#ifndef VIS_H_INCLUDED
#define VIS_H_INCLUDED

#include <cairo.h>
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/frame.h>
#include <libavutil/imgutils.h>
#include <particle.h>
#include <stdint.h>

typedef struct EncodingParams {
  char *video_codec;
  char *codec_priv_key;
  char *codec_priv_value;
} EncodingParams;

typedef struct EncodingContext {
  AVFormatContext *avfc;
  AVCodec *avc;
  AVStream *avs;
  AVCodecContext *avcc;
  char *filename;
} EncodingContext;

int _prepare_video_encoder(struct EncodingContext *ctx, int width, int height, struct EncodingParams p);

struct EncodingContext *create_encoding_context(char *filename, int width, int height);
int close_encoding_context(struct EncodingContext *ctx);

int _encode(struct EncodingContext *ctx, AVFrame *input_frame);

void _cairo_surface_to_frame(cairo_surface_t *surface, AVFrame *frame);

void render_test(struct EncodingContext *ctx);
void render_particles_to_encoding_context(struct Particle **z_curve_idx, struct Mesh *mesh, int c_size, struct EncodingContext *ctx);

#endif