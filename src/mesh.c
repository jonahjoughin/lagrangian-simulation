#include <main.h>
#include <netcdf.h>
#include <particle.h>
#include <sort.h>
#include <stdio.h>
#include <stdlib.h>

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

int build_mesh(int grd_ncid, int avg_ncid, struct Mesh *mesh) {
  int y_psi_id, x_psi_id, y_rho_id, x_rho_id, mask_rho_id, y_u_id, x_u_id, y_v_id, x_v_id;
  int u_id, v_id;

  size_t u_start[4] = {0, 0, 0, 0};
  size_t u_end[4] = {1, 1, NY_U, NX_U};
  size_t v_start[4] = {0, 0, 0, 0};
  size_t v_end[4] = {1, 1, NY_V, NX_V};
  // Get grid variable ids
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "y_psi", &y_psi_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "x_psi", &x_psi_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "y_rho", &y_rho_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "x_rho", &x_rho_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "mask_rho", &mask_rho_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "y_u", &y_u_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "x_u", &x_u_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "y_v", &y_v_id));
  CHECK_NC_NOERR(nc_inq_varid(grd_ncid, "x_v", &x_v_id));
  // Get u, v, variable ids
  CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "u", &u_id));
  CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "v", &v_id));
  // Copy grid arrays
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, y_psi_id, &(mesh->y_psi[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, x_psi_id, &(mesh->x_psi[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, y_rho_id, &(mesh->y_rho[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, x_rho_id, &(mesh->x_rho[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, mask_rho_id, &(mesh->mask_rho[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, y_u_id, &(mesh->y_u[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, x_u_id, &(mesh->x_u[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, y_v_id, &(mesh->y_v[0][0])));
  CHECK_NC_NOERR(nc_get_var_double(grd_ncid, x_v_id, &(mesh->x_v[0][0])));
  // Copy u, v arrays
  CHECK_NC_NOERR(nc_get_vara_double(avg_ncid, u_id, u_start, u_end, &(mesh->u[0][0])));
  CHECK_NC_NOERR(nc_get_vara_double(avg_ncid, v_id, v_start, v_end, &(mesh->v[0][0])));

  uint64_t **idx = (uint64_t **)calloc(NX_RHO * NY_RHO * NZ_RHO, sizeof(uint64_t *));
  // Get interleaved z_curve codes
  for (int i = 0; i < NX_RHO; i++) {
    for (int j = 0; j < NY_RHO; j++) {
      for (int k = 0; k < NZ_RHO; k++) {
        mesh->z_curve_idx_rho[k][j][i] = _xyz_to_z_curve(i, j, k);
        idx[k + j * NZ_RHO + i * NY_RHO * NZ_RHO] = &mesh->z_curve_idx_rho[k][j][i];
      }
    }
  }
  // Sort codes and reindex
  // This eliminates large number of unused codes, which speeds up sorting by an order of magnitude
  merge_sort((void **)idx, z_curve_idx_comparator, 0, NX_RHO * NY_RHO * NZ_RHO - 1);
  for (int i = 0; i < NX_RHO * NY_RHO * NZ_RHO; i++) {
    *(idx[i]) = i;
  }

  return 0;
}