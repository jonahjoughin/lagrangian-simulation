#include <interpolate.h>
#include <main.h>
#include <math.h>
#include <particle.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void create_particle(struct Particle *p, struct Mesh *mesh, uint64_t id, double x, double y, double z) {
  // Assign basic properties
  p->id = id;
  p->x = x, p->y = y, p->z = z;

  p->x_idx = 1, p->y_idx = 1, p->z_idx = 0;
  // Assign z curve index
  _update_particle_index(p, mesh);
}

void _z_curve_to_xyz(uint64_t v, uint64_t *x, uint64_t *y, uint64_t *z) {
  *x = 0;
  *y = 0;
  *z = 0;

  for (int i = 0; i < 48; i += 3) {
    *x |= (v & ((uint64_t)1 << (i))) >> (0 + 2 * i / 3);
    *y |= ((v & ((uint64_t)1 << (i + 1))) >> (1 + 2 * i / 3));
    *z |= ((v & ((uint64_t)1 << (i + 2))) >> (2 + 2 * i / 3));
  }
}

uint64_t _xyz_to_z_curve(uint64_t x, uint64_t y, uint64_t z) {
  uint64_t masks[5] = {(uint64_t)(pow(2, 16)) - 1, 0xff0000ff, 0xf00f00f00f, 0xc30c30c30c3, 0x249249249249};

  x = x & masks[0];
  x = (x | x << 16) & masks[1];
  x = (x | x << 8) & masks[2];
  x = (x | x << 4) & masks[3];
  x = (x | x << 2) & masks[4];

  y = y & masks[0];
  y = (y | y << 16) & masks[1];
  y = (y | y << 8) & masks[2];
  y = (y | y << 4) & masks[3];
  y = (y | y << 2) & masks[4];

  z = z & masks[0];
  z = (z | z << 16) & masks[1];
  z = (z | z << 8) & masks[2];
  z = (z | z << 4) & masks[3];
  z = (z | z << 2) & masks[4];

  return x + (y << 1) + (z << 2);
}

void _a_min_c_b_min_a(const double *p_a, const double *p_b, const double *p_c, double *vec_a, double *vec_b) {
  for (int i = 0; i < 2; i++) {
    vec_a[i] = p_a[i] - p_c[i];
    vec_b[i] = p_b[i] - p_a[i];
  }
}

int _cross_product_sign(const double *vec1, const double *vec2) {
  double val = vec1[0] * vec2[1] - vec1[1] * vec2[0];
  if (val == 0) {
    return 0;
  } else {
    return val > 0 ? 1 : -1;
  }
}

double _square(double x) {
  return x * x;
}

struct Particle **create_index(struct Particle *particles, struct Mesh *mesh, int p_size, uint64_t *c_size) {
  // This will get very large if we increase NX, NY, NZ
  // Currently, for instance, the max z_curve value is ~24x larger than NX*NY*NZ
  // Will eventually need to switch to precomputing and renumbering z_indices
  // This will also cut down computation from repeated interleaving and deinterleaving of bits
  *c_size = NX_RHO * NY_RHO * NZ_RHO;
  // Allocate index memory
  struct Particle **z_curve_idx = malloc(sizeof(struct Particle *) * (*c_size + 1));
  // Set end pointers of index
  z_curve_idx[0] = particles;
  z_curve_idx[*c_size] = particles + p_size;
  // Sort index (will set all other pointers)
  _sort_index(z_curve_idx, *c_size);
  return z_curve_idx;
}

int step_index(struct Particle **z_curve_idx, int c_size, struct Mesh *mesh, double dt) {
  int p_size = z_curve_idx[c_size] - z_curve_idx[0];
  int err_count = 0;

  for (struct Particle *p = *z_curve_idx; p < *z_curve_idx + p_size; p++) {
    // Get grid index
    _step_particle(p, mesh, dt);
    err_count += _update_particle_index(p, mesh);
  }

  _sort_index(z_curve_idx, c_size);
  return err_count;
}

void _sort_index(struct Particle **z_curve_idx, int c_size) {
  int p_size = (z_curve_idx[c_size] - z_curve_idx[0]);
  int *count = (int *)calloc(c_size, sizeof(int));

  for (struct Particle *p = *z_curve_idx; p < *z_curve_idx + p_size; p++) {
    count[p->z_curve_idx]++;
  }

  int total = 0, tmp;
  for (int i = 0; i < c_size; i++) {
    tmp = count[i];
    count[i] = total;
    total += tmp;
  }

  for (int i = 0; i < c_size; i++) {
    z_curve_idx[i] = *z_curve_idx + count[i];
  }

  struct Particle *particles_copy = (struct Particle *)malloc(sizeof(struct Particle) * p_size);
  memcpy(particles_copy, *z_curve_idx, sizeof(struct Particle) * p_size);

  for (int i = 0; i < p_size; i++) {
    (*z_curve_idx)[count[particles_copy[i].z_curve_idx]] = particles_copy[i];
    count[particles_copy[i].z_curve_idx] += 1;
  }

  free(particles_copy);
  free(count);
}

void _step_particle(struct Particle *p, struct Mesh *mesh, double dt) {
  int x = p->x_idx, y = p->y_idx, z = p->z_idx;

  double dx_dxi, dx_deta, dy_dxi, dy_deta, det, dxi_dt, deta_dt, dx_dt, dy_dt;
  double U0, U1, V0, V1;
  double v1[2], v2[2], v3[2], v4[2], coord[2], normalized_coord[2];

  // // Calculate fluxes
  U0 = mesh->u[y][x - 1] * sqrt(_square(mesh->y_psi[y - 1][x - 1] - mesh->y_psi[y][x - 1]) + _square(mesh->x_psi[y - 1][x - 1] - mesh->x_psi[y][x - 1]));
  U1 = mesh->u[y][x] * sqrt(_square(mesh->y_psi[y - 1][x] - mesh->y_psi[y][x]) + _square(mesh->x_psi[y - 1][x] - mesh->x_psi[y][x]));
  V0 = mesh->v[y - 1][x] * sqrt(_square(mesh->y_psi[y - 1][x - 1] - mesh->y_psi[y - 1][x]) + _square(mesh->x_psi[y - 1][x - 1] - mesh->x_psi[y - 1][x]));
  V1 = mesh->v[y][x] * sqrt(_square(mesh->y_psi[y][x - 1] - mesh->y_psi[y][x]) + _square(mesh->x_psi[y][x - 1] - mesh->x_psi[y][x]));

  // Set up vertex vectors
  v1[0] = mesh->y_psi[y - 1][x - 1];
  v1[1] = mesh->x_psi[y - 1][x - 1];
  v2[0] = mesh->y_psi[y][x - 1];
  v2[1] = mesh->x_psi[y][x - 1];
  v3[0] = mesh->y_psi[y][x];
  v3[1] = mesh->x_psi[y][x];
  v4[0] = mesh->y_psi[y - 1][x];
  v4[1] = mesh->x_psi[y - 1][x];
  coord[0] = p->y;
  coord[1] = p->x;

  // Get normalized coordinates
  inverse_bilinear(coord, v1, v2, v3, v4, normalized_coord);
  // Calculate partial derivatives with respect to normalized coordinates
  dx_deta = -v1[1] + v2[1] + (v1[1] - v2[1] + v3[1] - v4[1]) * normalized_coord[1];
  dy_deta = -v1[0] + v2[0] + (v1[0] - v2[0] + v3[0] - v4[0]) * normalized_coord[1];
  dx_dxi = -v1[1] + v4[1] + (v1[1] - v2[1] + v3[1] - v4[1]) * normalized_coord[0];
  dy_dxi = -v1[0] + v4[0] + (v1[0] - v2[0] + v3[0] - v4[0]) * normalized_coord[0];
  // Calculate derivatives of normalized coordinates with respect to time
  det = (dx_dxi * dy_deta) - (dy_dxi * dx_deta);
  dxi_dt = ((1 - normalized_coord[1]) * U0 + normalized_coord[1] * U1) / det;
  deta_dt = ((1 - normalized_coord[0]) * V0 + normalized_coord[0] * V1) / det;
  // Get derivatives with respect to time
  dx_dt = dx_dxi * dxi_dt + dx_deta * deta_dt;
  dy_dt = dy_dxi * dxi_dt + dy_deta * deta_dt;
  // Simple update using timestep
  p->x += dx_dt * dt;
  p->y += dy_dt * dt;
}

int _update_particle_index(struct Particle *particle, const struct Mesh *mesh) {
  double v1[2], v2[2], v3[2], v4[2], p_yx[2];
  double vec1[2], vec2[2];
  uint64_t idx[3];
  // Set initial approximation to last known index
  idx[0] = particle->y_idx, idx[1] = particle->x_idx, idx[2] = particle->z_idx;
  // Locate new index
  // This should take very few iterations (probably one), unless particles are crossing multiple grid cells in a single step
  while (1) {
    // Out of bounds
    if (idx[0] < 1 || idx[1] < 1 || idx[0] >= NY_RHO - 1 || idx[1] >= NX_RHO - 1) {
      return -1;
    }
    // Set up vertex vectors
    v1[0] = mesh->y_psi[idx[0] - 1][idx[1] - 1];
    v1[1] = mesh->x_psi[idx[0] - 1][idx[1] - 1];
    v2[0] = mesh->y_psi[idx[0]][idx[1] - 1];
    v2[1] = mesh->x_psi[idx[0]][idx[1] - 1];
    v3[0] = mesh->y_psi[idx[0]][idx[1]];
    v3[1] = mesh->x_psi[idx[0]][idx[1]];
    v4[0] = mesh->y_psi[idx[0] - 1][idx[1]];
    v4[1] = mesh->x_psi[idx[0] - 1][idx[1]];
    p_yx[0] = particle->y;
    p_yx[1] = particle->x;

    // Check if particle outside of grid cell, and update cell accordingly

    _a_min_c_b_min_a(v1, v4, p_yx, vec1, vec2);
    if (_cross_product_sign(vec2, vec1) == -1) {
      idx[0] -= 1;
      continue;
    }

    _a_min_c_b_min_a(v4, v3, p_yx, vec1, vec2);
    if (_cross_product_sign(vec2, vec1) == -1) {
      idx[1] += 1;
      continue;
    }

    _a_min_c_b_min_a(v3, v2, p_yx, vec1, vec2);
    if (_cross_product_sign(vec2, vec1) == -1) {
      idx[0] += 1;
      continue;
    }

    _a_min_c_b_min_a(v2, v1, p_yx, vec1, vec2);
    if (_cross_product_sign(vec2, vec1) == -1) {
      idx[1] -= 1;
      continue;
    }

    // Copy correct index and return
    particle->x_idx = idx[1], particle->y_idx = idx[0], particle->z_idx = idx[2];
    particle->z_curve_idx = mesh->z_curve_idx_rho[idx[2]][idx[0]][idx[1]];
    return 0;
  }
}
