#include <cairo.h>
#include <inttypes.h>
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/frame.h>
#include <libavutil/imgutils.h>
#include <libavutil/opt.h>
#include <libavutil/pixfmt.h>
#include <libavutil/timestamp.h>
#include <libswscale/swscale.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <visualize.h>

// Convert Cairo surface to video frame format
void _cairo_surface_to_frame(cairo_surface_t *surface, AVFrame *frame) {
  int in_linesizes[4];
  unsigned char *surface_buffer = cairo_image_surface_get_data(surface);

  AVFrame *in_frame = av_frame_alloc();
  in_frame->format = AV_PIX_FMT_BGRA;
  in_frame->width = cairo_image_surface_get_width(surface);
  in_frame->height = cairo_image_surface_get_height(surface);

  av_image_fill_arrays(in_frame->data, in_frame->linesize, surface_buffer, in_frame->format, in_frame->width, in_frame->height, 1);

  struct SwsContext *fooContext = sws_getContext(in_frame->width, in_frame->height, in_frame->format, frame->width, frame->height, frame->format, SWS_FAST_BILINEAR, NULL, NULL, NULL);
  //perform the conversion
  sws_scale(fooContext, (const uint8_t *const *)in_frame->data, in_frame->linesize, 0, frame->height, frame->data, frame->linesize);
  av_frame_free(&in_frame);
}

// Render particles to video frame via Cairo
void render_particles_to_encoding_context(struct Particle **z_curve_idx, struct Mesh *mesh, int c_size, EncodingContext *ctx) {
  int width = ctx->avcc->width;
  int height = ctx->avcc->height;

  // Set up properly sized surface
  cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
  cairo_t *cr = cairo_create(surface);

  // Clear background
  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_rectangle(cr, 0, 0, width, height);
  cairo_fill(cr);

  //Set blending mode and particle color
  cairo_set_source_rgba(cr, 0, 0, 0, 0.15);
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);

  int x_px,
      y_px;

  // Draw grid
  for (int i = 0; i < NX_PSI; i++) {
    for (int j = 0; j < NY_PSI - 1; j++) {
      x_px = (mesh->x_psi[j][i] / 1100000. * (double)width);
      y_px = ((1 - mesh->y_psi[j][i] / 800000.) * (double)height);
      cairo_move_to(cr, x_px, y_px);
      x_px = (mesh->x_psi[j + 1][i] / 1100000. * (double)width);
      y_px = ((1 - mesh->y_psi[j + 1][i] / 800000.) * (double)height);
      cairo_line_to(cr, x_px, y_px);
    }
  }

  for (int j = 0; j < NY_PSI; j++) {
    for (int i = 0; i < NX_PSI - 1; i++) {
      x_px = (mesh->x_psi[j][i] / 1100000. * (double)width);
      y_px = ((1 - mesh->y_psi[j][i] / 800000.) * (double)height);
      cairo_move_to(cr, x_px, y_px);
      x_px = (mesh->x_psi[j][i + 1] / 1100000. * (double)width);
      y_px = ((1 - mesh->y_psi[j][i + 1] / 800000.) * (double)height);
      cairo_line_to(cr, x_px, y_px);
    }
  }

  cairo_stroke(cr);
  cairo_set_source_rgba(cr, 0, 0, 0, 1);

  int p_size = z_curve_idx[c_size] - z_curve_idx[0];

  // Draw particles
  for (struct Particle *p = *z_curve_idx; p < *z_curve_idx + p_size; p++) {
    // Get grid index
    x_px = (p->x / 1100000. * (double)width);
    y_px = ((1 - p->y / 800000.) * (double)height);
    cairo_move_to(cr, x_px, y_px);
    cairo_line_to(cr, x_px, y_px);
  }

  cairo_stroke(cr);
  cairo_surface_flush(surface);

  // Set up frame
  AVFrame *frame = av_frame_alloc();
  frame->format = AV_PIX_FMT_YUV420P;
  frame->width = ctx->avcc->width;
  frame->height = ctx->avcc->height;
  av_frame_get_buffer(frame, 0);

  // Convert to frame and encode
  _cairo_surface_to_frame(surface, frame);
  _encode(ctx, frame);

  // Clean up
  av_frame_free(&frame);
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
}

// Create encoding context
struct EncodingContext *create_encoding_context(char *filename, int width, int height) {
  struct EncodingParams sp = {0};
  sp.video_codec = "libx265";
  sp.codec_priv_key = "x265-params";
  sp.codec_priv_value = "keyint=60:min-keyint=60:scenecut=0";

  struct EncodingContext *ctx = (struct EncodingContext *)calloc(1, sizeof(struct EncodingContext));
  ctx->filename = malloc(sizeof(char) * (strlen(filename) + 1));
  strcpy(ctx->filename, filename);

  avformat_alloc_output_context2(&ctx->avfc, NULL, NULL, ctx->filename);
  if (!ctx->avfc) {
    printf("could not allocate memory for output format");
    free(ctx);
    return NULL;
  }

  _prepare_video_encoder(ctx, width, height, sp);

  if (!(ctx->avfc->oformat->flags & AVFMT_NOFILE)) {
    if (avio_open(&ctx->avfc->pb, ctx->filename, AVIO_FLAG_WRITE) < 0) {
      printf("could not open the output file");
      free(ctx);
      return NULL;
    }
  }

  if (avformat_write_header(ctx->avfc, NULL) < 0) {
    printf("an error occurred when opening output file");
    free(ctx);
    return NULL;
  }

  return ctx;
}

// Set up encoder
int _prepare_video_encoder(struct EncodingContext *ctx, int width, int height, struct EncodingParams sp) {
  ctx->avs = avformat_new_stream(ctx->avfc, NULL);

  ctx->avc = avcodec_find_encoder_by_name(sp.video_codec);
  if (!ctx->avc) {
    printf("could not find the proper codec");
    return -1;
  }

  ctx->avcc = avcodec_alloc_context3(ctx->avc);
  if (!ctx->avcc) {
    printf("could not allocated memory for codec context");
    return -1;
  }

  av_opt_set(ctx->avcc->priv_data, "preset", "medium", 0);
  if (sp.codec_priv_key && sp.codec_priv_value)
    av_opt_set(ctx->avcc->priv_data, sp.codec_priv_key, sp.codec_priv_value, 0);

  ctx->avcc->height = height;
  ctx->avcc->width = width;
  ctx->avcc->pix_fmt = AV_PIX_FMT_YUV420P;

  ctx->avcc->bit_rate = 2 * 1000 * 1000 * 5;
  ctx->avcc->rc_buffer_size = 4 * 1000 * 1000 * 5;
  ctx->avcc->rc_max_rate = 2.5 * 1000 * 1000 * 5;
  ctx->avcc->rc_min_rate = 2 * 1000 * 1000 * 5;

  ctx->avcc->time_base = (AVRational){1, 30};
  ctx->avcc->framerate = (AVRational){30, 1};

  if (avcodec_open2(ctx->avcc, ctx->avc, NULL) < 0) {
    printf("could not open the codec");
    return -1;
  }
  avcodec_parameters_from_context(ctx->avs->codecpar, ctx->avcc);
  ctx->avs->codecpar->codec_tag = MKTAG('h', 'v', 'c', '1');
  return 0;
}

// Encode arbitrary frame to encoding context
int _encode(struct EncodingContext *ctx, AVFrame *input_frame) {
  if (input_frame) input_frame->pict_type = AV_PICTURE_TYPE_NONE;
  AVPacket *output_packet = av_packet_alloc();
  if (!output_packet) {
    printf("could not allocate memory for output packet");
    return -1;
  }
  if (input_frame != NULL)
    input_frame->pts = ctx->avcc->frame_number;

  int response = avcodec_send_frame(ctx->avcc, input_frame);

  while (response >= 0) {
    response = avcodec_receive_packet(ctx->avcc, output_packet);
    if (response == AVERROR(EAGAIN) || response == AVERROR_EOF) {
      break;
    } else if (response < 0) {
      printf("Error while receiving packet from encoder: %s", av_err2str(response));
      return -1;
    }

    if (output_packet->pts != AV_NOPTS_VALUE)
      output_packet->pts = av_rescale_q(output_packet->pts, ctx->avcc->time_base, ctx->avs->time_base);
    if (output_packet->dts != AV_NOPTS_VALUE)
      output_packet->dts = av_rescale_q(output_packet->dts, ctx->avcc->time_base, ctx->avs->time_base);

    response = av_interleaved_write_frame(ctx->avfc, output_packet);

    if (response != 0) {
      printf("Error %d while receiving packet from decoder: %s", response, av_err2str(response));
      return -1;
    }
  }
  av_packet_free(&output_packet);
  return 0;
}

// Clean up encoding context
int close_encoding_context(struct EncodingContext *ctx) {
  if (_encode(ctx, NULL)) return -1;

  av_write_trailer(ctx->avfc);

  avformat_free_context(ctx->avfc);
  ctx->avfc = NULL;

  free(ctx);
  ctx = NULL;

  return 0;
}