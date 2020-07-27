#include <interpolate.h>
#include <main.h>
#include <math.h>
#include <particle.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))

void create_particle(struct Particle *p, struct Mesh *mesh, uint64_t id, double x, double y, double z) {
  // Assign basic properties
  p->id = id;
  p->x = x, p->y = y, p->z = z;
  p->x_idx = -1, p->y_idx = -1, p->z_idx = -1;
  // Assign z curve index
  _update_particle_index(p, mesh);
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

// Helper function for stencil jumping algorithm
// v1 := p1 - p3, v2 := p2 - p1
void _a_min_c_b_min_a(const double *p_1, const double *p_2, const double *p_3, double *vec_1, double *vec_2) {
  for (int i = 0; i < 2; i++) {
    vec_1[i] = p_1[i] - p_3[i];
    vec_2[i] = p_2[i] - p_1[i];
  }
}

// Helper function for stencil jumping algorithm
// Computes the sign of the cross product (+/- 1)
int _cross_product_sign(const double *vec1, const double *vec2) {
  double val = vec1[0] * vec2[1] - vec1[1] * vec2[0];
  if (val == 0) {
    return 0;
  } else {
    return val > 0 ? 1 : -1;
  }
}

// Reinsert particle that has crossed boundary at random on the boundary
void reinsert_particle_at_open_boundary(struct Particle *p, struct Mesh *mesh) {

  struct Boundary *b = &(mesh->open_boundary);
  struct BoundarySegment *s;
  double d;

  // Continue until sampling process succeeds
  while (1) {
    // Rejection Sampling
    s = &(mesh->open_boundary.segments[rand() % b->segments_count]);
    d = ((double)rand()) / (RAND_MAX / b->max_flux);
    if (d < s->flux) {
      break;
    }
  }

  // Small perturbation to avoid placing particles precisely on boundary
  double eta = 0.001;
  double normalized_coord[2];
  double v1[2], v2[2], v3[3], v4[4];
  int x_idx, y_idx, z_idx = s->z_idx;

  // Sample particle uniformly along correct edge of normalized grid cell [0, 1] x [0,1]
  // Bilinear interpolation will be applied to transform this normalized coordinate into actual coordinate
  if (s->axis == 'u') {
    normalized_coord[0] = ((double)rand()) / (RAND_MAX);
    if (s->orientation == 1) {
      normalized_coord[1] = eta;
      x_idx = s->x_idx + 1;
      y_idx = s->y_idx;
    } else {
      normalized_coord[1] = 1 - eta;
      x_idx = s->x_idx;
      y_idx = s->y_idx;
    }
  } else if (s->axis == 'v') {
    normalized_coord[1] = ((double)rand()) / (RAND_MAX);
    if (s->orientation == 1) {
      normalized_coord[0] = eta;
      x_idx = s->x_idx;
      y_idx = s->y_idx + 1;
    } else {
      normalized_coord[0] = 1 - eta;
      x_idx = s->x_idx;
      y_idx = s->y_idx;
    }
  }

  // Assign vectors for bilinear interpolation
  v1[0] = mesh->y_psi[y_idx - 1][x_idx - 1];
  v1[1] = mesh->x_psi[y_idx - 1][x_idx - 1];
  v2[0] = mesh->y_psi[y_idx][x_idx - 1];
  v2[1] = mesh->x_psi[y_idx][x_idx - 1];
  v3[0] = mesh->y_psi[y_idx][x_idx];
  v3[1] = mesh->x_psi[y_idx][x_idx];
  v4[0] = mesh->y_psi[y_idx - 1][x_idx];
  v4[1] = mesh->x_psi[y_idx - 1][x_idx];

  // Sample z coordinate uniformly
  p->z = mesh->z_w[z_idx][y_idx][x_idx] + ((double)rand()) / (RAND_MAX / (mesh->z_w[z_idx + 1][y_idx][x_idx] - mesh->z_w[z_idx][y_idx][x_idx]));
  // Calculate x, y through bilinear interpolation
  p->y = bilinear(normalized_coord, v1[0], v2[0], v3[0], v4[0]);
  p->x = bilinear(normalized_coord, v1[1], v2[1], v3[1], v4[1]);

  // Assign new indices
  p->z_idx = z_idx;
  p->y_idx = y_idx;
  p->x_idx = x_idx;
  p->z_curve_idx = mesh->z_curve_idx_rho[z_idx][y_idx][x_idx];
}

// Create z-order index
struct Particle **create_index(struct Particle *particles, struct Mesh *mesh, int p_size, uint64_t *c_size) {
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

// Update all particles in a particle index
int step_index(struct Particle **z_curve_idx, int c_size, struct Mesh *mesh, double dt) {

  int p_size = z_curve_idx[c_size] - z_curve_idx[0];

  for (struct Particle *p = *z_curve_idx; p < *z_curve_idx + p_size; p++) {
    // Get grid index
    _step_particle(p, mesh, dt);
    _update_particle_index(p, mesh);

    // Check if particle has crossed open boundary and relocate if necessary
    if (mesh->mask_rho[p->y_idx][p->x_idx] == -1) {
      reinsert_particle_at_open_boundary(p, mesh);
    }
  }

  // Reorder index based on updated locations
  _sort_index(z_curve_idx, c_size);
  return 0;
}

void _sort_index(struct Particle **z_curve_idx, int c_size) {
  // O(n + m) sort
  // n: number of particles
  // m: number of values in index
  int p_size = (z_curve_idx[c_size] - z_curve_idx[0]);
  int *count = (int *)calloc(c_size, sizeof(int));

  // Construct histogram
  for (struct Particle *p = *z_curve_idx; p < *z_curve_idx + p_size; p++) {
    count[p->z_curve_idx]++;
  }

  // Construct cumulative histogram
  int total = 0, tmp;
  for (int i = 0; i < c_size; i++) {
    tmp = count[i];
    count[i] = total;
    total += tmp;
  }

  // Update index pointers
  for (int i = 0; i < c_size; i++) {
    z_curve_idx[i] = *z_curve_idx + count[i];
  }

  struct Particle *particles_copy = (struct Particle *)malloc(sizeof(struct Particle) * p_size);
  memcpy(particles_copy, *z_curve_idx, sizeof(struct Particle) * p_size);

  // Reorder particles
  for (int i = 0; i < p_size; i++) {
    (*z_curve_idx)[count[particles_copy[i].z_curve_idx]] = particles_copy[i];
    count[particles_copy[i].z_curve_idx] += 1;
  }

  free(particles_copy);
  free(count);
}

// Update individual particle
void _step_particle(struct Particle *p, struct Mesh *mesh, double dt) {

  int x = p->x_idx, y = p->y_idx, z = p->z_idx;
  double dx_dxi, dx_deta, dy_dxi, dy_deta, det, dxi_dt, deta_dt, dx_dt, dy_dt, dz_dt, zeta;
  double U0, U1, V0, V1;
  double v1[2], v2[2], v3[2], v4[2], coord[2], normalized_coord[2];

  // // Calculate fluxes
  U0 = mesh->u[z][y][x - 1] * sqrt(SQUARE(mesh->y_psi[y - 1][x - 1] - mesh->y_psi[y][x - 1]) + SQUARE(mesh->x_psi[y - 1][x - 1] - mesh->x_psi[y][x - 1]));
  U1 = mesh->u[z][y][x] * sqrt(SQUARE(mesh->y_psi[y - 1][x] - mesh->y_psi[y][x]) + SQUARE(mesh->x_psi[y - 1][x] - mesh->x_psi[y][x]));
  V0 = mesh->v[z][y - 1][x] * sqrt(SQUARE(mesh->y_psi[y - 1][x - 1] - mesh->y_psi[y - 1][x]) + SQUARE(mesh->x_psi[y - 1][x - 1] - mesh->x_psi[y - 1][x]));
  V1 = mesh->v[z][y][x] * sqrt(SQUARE(mesh->y_psi[y][x - 1] - mesh->y_psi[y][x]) + SQUARE(mesh->x_psi[y][x - 1] - mesh->x_psi[y][x]));

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
  // Calculate zeta coordinate
  zeta = (p->z - mesh->z_w[z][y][x]) / (mesh->z_w[z + 1][y][x] - mesh->z_w[z][y][x]);

  // Get derivatives with respect to time
  dx_dt = dx_dxi * dxi_dt + dx_deta * deta_dt;
  dy_dt = dy_dxi * dxi_dt + dy_deta * deta_dt;
  dz_dt = zeta * mesh->w[z][y][x] + (1 - zeta) * mesh->w[z + 1][y][x];

  // Simple update using timestep
  p->x += dx_dt * dt;
  p->y += dy_dt * dt;
  p->z += dz_dt * dt;
}

int _update_particle_index(struct Particle *particle, const struct Mesh *mesh) {
  int retval = 0;
  double v1[2], v2[2], v3[2], v4[2], p_yx[2];
  double vec1[2], vec2[2];
  int idx[3];
  // Set initial approximation to last known index
  idx[0] = MAX(particle->y_idx, 1), idx[1] = MAX(particle->x_idx, 1), idx[2] = MAX(particle->z_idx, 0);
  // Locate new index
  // This should take very few iterations (probably one), unless particles are crossing multiple grid cells in a single step
  while (1) {
    // Out of bounds
    if (idx[0] < 1 || idx[1] < 1 || idx[0] >= NY_RHO - 1 || idx[1] >= NX_RHO - 1) {
      particle->x_idx = -1;
      particle->y_idx = -1;
      particle->z_idx = -1;
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

    break;
  }

  // Find vertical index of particle
  // Initial test version, could be improved upon with binary search
  while (1) {
    // Restrict particle z-coordinate to vertical range
    if (idx[2] < 0) {
      particle->z = mesh->z_w[0][idx[0]][idx[1]] + .001;
      idx[2] = 0;
      retval = -1;
    } else if (idx[2] >= NZ_RHO) {
      particle->z = mesh->z_w[NZ_RHO][idx[0]][idx[1]] - .001;
      idx[2] = NZ_RHO - 1;
      retval = -1;
    }
    if (particle->z < mesh->z_w[idx[2]][idx[0]][idx[1]]) {
      idx[2] -= 1;
      continue;
    } else if (particle->z > mesh->z_w[idx[2] + 1][idx[0]][idx[1]]) {
      idx[2] += 1;
      continue;
    }
    break;
  }

  // Copy correct index and return
  particle->x_idx = idx[1];
  particle->y_idx = idx[0];
  particle->z_idx = idx[2];
  particle->z_curve_idx = mesh->z_curve_idx_rho[idx[2]][idx[0]][idx[1]];

  return retval;
}
