#include <interpolate.h>
#include <main.h>
#include <netcdf.h>
#include <particle.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <visualize.h>

#define GRD_PATH "/path/to/croco_grd.nc"
#define AVG_PATH "/path/to/croco_avg.nc"
#define VIS_PATH "/path/to/croco_vis.mp4"

double rand_in_range(double min, double max) {
  return ((double)rand() / (double)(RAND_MAX) * (max - min)) + min;
}

int main() {
  int grd_ncid, avg_ncid;
  uint64_t p_size = 5000, c_size;
  struct Mesh *mesh = malloc(sizeof(struct Mesh));
  struct Particle *particles = malloc(sizeof(struct Particle) * p_size);

  srand(time(NULL));

  // Open files and build mesh
  CHECK_NC_NOERR(nc_open(GRD_PATH, NC_WRITE, &grd_ncid));
  CHECK_NC_NOERR(nc_open(AVG_PATH, NC_WRITE, &avg_ncid));

  // Build mesh and load velocity field
  initialize_mesh(grd_ncid, mesh);
  update_mesh(avg_ncid, mesh, 0);

  // Sample particles from uniform distribution
  // This can be quite slow, especially in testing
  // One alternative is to generate particles once, and save them to a binary file
  for (struct Particle *p = particles; p < particles + p_size; p++) {
    while (1) {
      create_particle(p, mesh, p - particles, rand_in_range(0, 1100000), rand_in_range(0, 800000), rand_in_range(-1650, 1));
      if (_update_particle_index(p, mesh) != -1 && mesh->mask_rho[p->y_idx][p->x_idx] == 1)
        break;
    }
  }

  // Build z-curve index index
  struct Particle **z_curve_idx =
      create_index(particles, mesh, p_size, &c_size);

  // Create video encoding context for visualization
  EncodingContext *ctx = create_encoding_context(VIS_PATH, 2048, 1024);

  for (int i = 0; i < 720; i++) {
    // Update velocity field every 24 hours
    if (i % 24 == 0) {
      update_mesh(avg_ncid, mesh, i / 24);
    }
    // Step forward an hour in 4 minute increments
    for (int j = 0; j < 15; j++) {
      step_index(z_curve_idx, c_size, mesh, 240);
    }
    // Write frame to visualization
    render_particles_to_encoding_context(z_curve_idx, mesh, c_size, ctx);
  }
  // Clean up video encoding context
  close_encoding_context(ctx);
}