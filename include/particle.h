#include <main.h>
#include <mesh.h>
#include <stdint.h>
#include <stdio.h>

struct Particle {
  uint64_t id;
  uint64_t z_curve_idx;
  int x_idx, y_idx, z_idx;
  double x, y, z;
};

void create_particle(struct Particle *p, struct Mesh *mesh, uint64_t id, double x, double y, double z);

void _z_curve_to_xyz(uint64_t v, uint64_t *x, uint64_t *y, uint64_t *z);
uint64_t _xyz_to_z_curve(uint64_t x, uint64_t y, uint64_t z);

void _assign_vectors(const double *p_a, const double *p_b, const double *p_c, double *vec_a, double *vec_b);
int _cross_product_sign(const double *vec1, const double *vec2);

struct Particle **create_index(struct Particle *particles, struct Mesh *mesh, int p_size, uint64_t *c_size);
int step_index(struct Particle **z_curve_idx, int c_size, struct Mesh *mesh, double dt);
void _sort_index(struct Particle **z_curve_idx, int c_size);

void _step_particle(struct Particle *p, struct Mesh *mesh, double dt);
int _update_particle_index(struct Particle *particle, const struct Mesh *mesh);
