#include <main.h>
#include <math.h>
#include <netcdf.h>
#include <particle.h>
#include <stdio.h>
#include <stdlib.h>

// Used to sort grid cells by z-curve index
int z_curve_idx_comparator(const void *a_, const void *b_) {
  uint64_t *a = (uint64_t *)a_, *b = (uint64_t *)b_;
  if (*a < *b) {
    return -1;
  } else if (*a == *b) {
    return 0;
  } else {
    return 1;
  }
}

// Build boundary from mesh
void initialize_mesh_open_boundary(struct Mesh *mesh) {
  int s_count = 0, s_alloc_count = 1;
  struct BoundarySegment *segments = malloc(sizeof(struct BoundarySegment) * s_alloc_count);
  int i, j, k;

  // Loop through all u faces
  for (i = 0; i < NX_RHO - 1; i++) {
    for (j = 0; j < NY_RHO; j++) {
      for (k = 0; k < NZ_RHO; k++) {
        // Check if open boundary segment on a u face is found
        // A segment is considered part of an open boundary it is between a cell with mask value 1 and a cell with mask value -1
        if (mesh->mask_rho[j][i] * mesh->mask_rho[j][i + 1] == -1) {

          // Grow size of segments array as necessary
          if (s_count == s_alloc_count) {
            s_alloc_count *= 2;
            segments = realloc(segments, sizeof(struct BoundarySegment) * s_alloc_count);
          }

          // Initialize segment data
          segments[s_count].axis = 'u';
          segments[s_count].x_idx = i;
          segments[s_count].y_idx = j;
          segments[s_count].z_idx = k;
          segments[s_count].orientation = (mesh->mask_rho[j][i + 1] == 1) ? 1 : -1;

          s_count += 1;
        }
      }
    }
  }

  // Loop through all v faces
  for (i = 0; i < NX_RHO; i++) {
    for (j = 0; j < NY_RHO - 1; j++) {
      for (k = 0; k < NZ_RHO; k++) {
        // Check if open boundary segment on a v face is found
        // A segment is considered part of an open boundary it is between a cell with mask value 1 and a cell with mask value -1
        if (mesh->mask_rho[j][i] * mesh->mask_rho[j + 1][i] == -1) {

          // Grow size of segments array as necessary
          if (s_count == s_alloc_count) {
            s_alloc_count *= 2;
            segments = realloc(segments, sizeof(struct BoundarySegment) * s_alloc_count);
          }

          // Initialize segment data
          segments[s_count].axis = 'v';
          segments[s_count].x_idx = i;
          segments[s_count].y_idx = j;
          segments[s_count].z_idx = k;
          segments[s_count].orientation = (mesh->mask_rho[j + 1][i] == 1) ? 1 : -1;

          s_count += 1;
        }
      }
    }
  }
  // Resize array once exact length is known
  mesh->open_boundary.segments = realloc(segments, sizeof(struct Boundary) * s_count);
  mesh->open_boundary.segments_count = s_count;
}

int initialize_mesh(int grd_ncid, struct Mesh *mesh) {

  int y_psi_id, x_psi_id, mask_rho_id, y_u_id, x_u_id, y_v_id, x_v_id, z_w_id;

  // Get grid variable ids
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "y_psi", &y_psi_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "x_psi", &x_psi_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "mask_rho", &mask_rho_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "y_u", &y_u_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "x_u", &x_u_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "y_v", &y_v_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "x_v", &x_v_id));

  // Copy grid arrays
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, y_psi_id, &(mesh->y_psi[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, x_psi_id, &(mesh->x_psi[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, mask_rho_id, &(mesh->mask_rho[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, y_u_id, &(mesh->y_u[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, x_u_id, &(mesh->x_u[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, y_v_id, &(mesh->y_v[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, x_v_id, &(mesh->x_v[0][0])));

  // Set up boundary
  initialize_mesh_open_boundary(mesh);

  int i, j, k;
  uint64_t **idx = (uint64_t **)calloc(NX_RHO * NY_RHO * NZ_RHO, sizeof(uint64_t *));

  // Get interleaved z_curve codes
  for (i = 0; i < NX_RHO; i++) {
    for (j = 0; j < NY_RHO; j++) {
      for (k = 0; k < NZ_RHO; k++) {
        mesh->z_curve_idx_rho[k][j][i] = _xyz_to_z_curve(i, j, k);
        idx[k + j * NZ_RHO + i * NY_RHO * NZ_RHO] = &mesh->z_curve_idx_rho[k][j][i];
      }
    }
  }

  // Sort codes and reindex
  // This eliminates large number of unused codes, which speeds up sorting by an order of magnitude
  qsort((void *)idx, NX_RHO * NY_RHO * NZ_RHO, sizeof(uint64_t *), z_curve_idx_comparator);
  for (int i = 0; i < NX_RHO * NY_RHO * NZ_RHO; i++) {
    *(idx[i]) = i;
  }

  free(idx);
  return 0;
}

void update_mesh_open_boundary(struct Mesh *mesh) {

  struct Boundary *b = &(mesh->open_boundary);
  struct BoundarySegment *s;
  int x, y, z;
  double surface_area;

  // Reset max (inward) flux
  b->max_flux = 0;

  // Iterate throguh segments
  for (s = b->segments; s < b->segments + b->segments_count; s++) {

    x = s->x_idx;
    y = s->y_idx;
    z = s->z_idx;

    // Calculate updated flux values for each segment
    if (s->axis == 'u') {
      surface_area = sqrt(SQUARE(mesh->y_psi[y - 1][x] - mesh->y_psi[y][x]) + SQUARE(mesh->x_psi[y - 1][x] - mesh->x_psi[y][x])) * (mesh->z_w[z + 1][y][x] - mesh->z_w[z][y][x]);
      s->flux = surface_area * mesh->u[z][y][x] * s->orientation;
    } else if (s->axis == 'v') {
      surface_area = sqrt(SQUARE(mesh->y_psi[y][x - 1] - mesh->y_psi[y][x]) + SQUARE(mesh->x_psi[y][x - 1] - mesh->x_psi[y][x])) * (mesh->z_w[z + 1][y][x] - mesh->z_w[z][y][x]);
      s->flux = surface_area * mesh->v[z][y][x] * s->orientation;
    }
    // Update max flux if necessary
    if (s->flux > b->max_flux) b->max_flux = s->flux;
  }
}

int update_mesh(int avg_ncid, struct Mesh *mesh, int t_idx) {
  int u_id, v_id, w_id;
  // Get u, v, variable ids
  CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "u", &u_id));
  CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "v", &v_id));
  CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "w", &w_id));

  // Set up array slices based on t_idx (time)
  size_t u_start[4] = {t_idx, 0, 0, 0};
  size_t u_end[4] = {1, NZ_U, NY_U, NX_U};
  size_t v_start[4] = {t_idx, 0, 0, 0};
  size_t v_end[4] = {1, NZ_V, NY_V, NX_V};
  size_t w_start[4] = {t_idx, 0, 0, 0};
  size_t w_end[4] = {1, NZ_W - 1, NY_W, NX_W};

  // Copy slices of u, v, w arrays
  CHECK_NC_NOERR(nc_get_vara_double(avg_ncid, u_id, u_start, u_end, &(mesh->u[0][0][0])));
  CHECK_NC_NOERR(nc_get_vara_double(avg_ncid, v_id, v_start, v_end, &(mesh->v[0][0][0])));
  CHECK_NC_NOERR(nc_get_vara_double(avg_ncid, w_id, w_start, w_end, &(mesh->w[0][0][0])));

  int hc_id, s_w_id, h_id, Cs_w_id, zeta_id;
  double hc, *s_w, *h, *Cs_w, *zeta;

  s_w = malloc(sizeof(double) * NZ_W);
  Cs_w = malloc(sizeof(double) * NZ_W);
  h = malloc(sizeof(double) * NY_W * NX_W);
  zeta = malloc(sizeof(double) * NY_W * NX_W);

  CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "hc", &hc_id));
  CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "s_w", &s_w_id));
  CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "Cs_w", &Cs_w_id));
  CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "h", &h_id));
  CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "zeta", &zeta_id));

  CHECK_NC_NOERR(nc_get_var_double(avg_ncid, hc_id, &hc));
  CHECK_NC_NOERR(nc_get_var_double(avg_ncid, s_w_id, s_w));
  CHECK_NC_NOERR(nc_get_var_double(avg_ncid, Cs_w_id, Cs_w));
  CHECK_NC_NOERR(nc_get_var_double(avg_ncid, h_id, h));

  size_t zeta_start[4] = {t_idx, 0, 0};
  size_t zeta_end[4] = {1, NY_W, NX_W};

  CHECK_NC_NOERR(nc_get_vara_double(avg_ncid, zeta_id, zeta_start, zeta_end, zeta));

  int i, j, k;
  double z_;

  for (k = 0; k < NZ_W; k++) {
    for (j = 0; j < NY_W; j++) {
      for (i = 0; i < NX_W; i++) {
        z_ = hc * s_w[k] + (h[j * NX_W + i] - hc) * Cs_w[k];
        mesh->z_w[k][j][i] = z_ + zeta[j * NX_W + i] * (1 + z_ / h[j * NX_W + i]);
      }
    }
  }

  // Update flux on open boundary
  update_mesh_open_boundary(mesh);

  // Clean up
  free(s_w);
  free(Cs_w);
  free(h);
  free(zeta);

  return 0;
}