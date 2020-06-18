#include <interpolate.h>
#include <main.h>
#include <netcdf.h>
#include <particle.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define GRD_PATH "/Users/jonah/Documents/work/pptp/netcdf/data/croco_grd.nc"
#define AVG_PATH "/Users/jonah/Documents/work/pptp/netcdf/data/croco_avg.nc"

double rand_in_range(double min, double max) {
  return ((double)rand() / (double)(RAND_MAX) * (max - min)) + min;
}

int main() {
  int grd_ncid, avg_ncid;
  uint64_t p_size = 100000, c_size;
  struct Mesh mesh;
  struct Particle *particles = malloc(sizeof(struct Particle) * p_size);

  // Open files and build mesh
  CHECK_NC_NOERR(nc_open(GRD_PATH, NC_NOWRITE, &grd_ncid));
  CHECK_NC_NOERR(nc_open(AVG_PATH, NC_NOWRITE, &avg_ncid));
  build_mesh(grd_ncid, avg_ncid, &mesh);

  // Set random seed
  srand(42);

  for (struct Particle *p = particles; p < particles + p_size; p++) {
    // Generate random particles in region of gulf
    create_particle(p, &mesh, p - particles, rand_in_range(400000, 500000),
                    rand_in_range(200000, 300000), rand_in_range(0, 100));
  }

  // Create index
  struct Particle **z_curve_idx =
      create_index(particles, &mesh, p_size, &c_size);

  // One month with dt=1 hour
  for (int i = 0; i < 8760; i++) {
    step_index(z_curve_idx, c_size, &mesh, 3600.0);
  }
}