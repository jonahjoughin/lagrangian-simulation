#include <main.h>
#include <netcdf.h>
#include <stdio.h>

#include <mesh.h>

int build_mesh(int grd_ncid, int avg_ncid, struct Mesh *mesh) {

    int y_psi_id, x_psi_id, y_u_id, x_u_id, y_v_id, x_v_id, y_rho_id, x_rho_id;
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
    CHECK_NC_NOERR(nc_get_var_double(grd_ncid, y_u_id, &(mesh->y_u[0][0])));
    CHECK_NC_NOERR(nc_get_var_double(grd_ncid, x_u_id, &(mesh->x_u[0][0])));
    CHECK_NC_NOERR(nc_get_var_double(grd_ncid, y_v_id, &(mesh->y_v[0][0])));
    CHECK_NC_NOERR(nc_get_var_double(grd_ncid, x_v_id, &(mesh->x_v[0][0])));
    // Copy u, v arrays
    CHECK_NC_NOERR(nc_get_vara_double(avg_ncid, u_id, u_start, u_end, &(mesh->u[0][0])));
    CHECK_NC_NOERR(nc_get_vara_double(avg_ncid, v_id, v_start, v_end, &(mesh->v[0][0])));

    return 0;
}