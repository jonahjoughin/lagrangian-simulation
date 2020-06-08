#include <main.h>
#include <stdio.h>
#include <netcdf.h>

int build_coord_mesh(int ncid, int pm_id, int pn_id, int ny, int nx, float *coord_mesh) {
    
    float pm[ny][nx], pn[ny][nx];
    float (*mesh)[ny+1][nx+1] = (float (*)[ny+1][nx+1]) coord_mesh;

    CHECK_NC_NOERR(nc_get_var_float(ncid, pm_id, &pm[0][0]));
    CHECK_NC_NOERR(nc_get_var_float(ncid, pn_id, &pn[0][0]));
    
    // Fill in mesh coordinates from pn, pm
    for (int i = 0; i <= ny; i++) {
       for (int j = 0; j <= nx; j++) {
           if (i == 0) {
               mesh[0][i][j] = 0;
           } else if (j == 0) {
               mesh[0][i][j] = mesh[0][i-1][j] + 1/pn[i-1][j];
           } else if (j == nx) {
               mesh[0][i][j] = mesh[0][i-1][j] + 1/pn[i-1][j-1];
           } else {
               mesh[0][i][j] = mesh[0][i-1][j] + (1/pn[i-1][j] + 1/pn[i-1][j-1])/2;
           }
           if (j == 0) {
               mesh[1][i][j] = 0;
           } else if (i == 0) {
              mesh[1][i][j] = mesh[1][i][j-1] + 1/pm[i][j-1];
           } else if (i == ny) {
               mesh[1][i][j] = mesh[1][i][j-1] + 1/pm[i-1][j-1];
           } else {
               mesh[1][i][j] = mesh[1][i][j-1] + (1/pm[i][j-1] + 1/pm[i-1][j-1])/2;
           }
       }
    }

    return 0;
}

