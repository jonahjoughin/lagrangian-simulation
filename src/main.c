#include <main.h>
#include <locate.h>
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <mesh.h>
#include <time.h>

#define GRD_PATH "../data/croco_grd.nc"
#define AVG_PATH "../data/croco_avg.nc"
#define NX 164
#define NY 130
#define NP 1000000

int main() {

    int grd_ncid, avg_ncid, pm_id, pn_id;
    float coord_mesh[2][NY+1][NX+1];

    // Get file ids
    CHECK_NC_NOERR(nc_open(GRD_PATH, NC_NOWRITE, &grd_ncid));
    CHECK_NC_NOERR(nc_open(AVG_PATH, NC_NOWRITE, &avg_ncid));

    // GET var ids
    CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "pm", &pm_id));
    CHECK_NC_NOERR(nc_inq_varid(avg_ncid, "pn", &pn_id));

    // Build coordinate mesh from pm, pn variables
    build_coord_mesh(avg_ncid, pm_id, pn_id, NY, NX, (float *) coord_mesh);
    srand(time(NULL));

    clock_t t; 
    t = clock(); 

    float particle_coord[2];
    int particle_idx[2];

    for (int i = 0; i < NP; i++) {
        particle_coord[0] = (float) (rand() % 700000);
        particle_coord[1] = (float) (rand() % 700000);
        locate_particle((float *) coord_mesh, NY, NX, particle_coord, particle_idx);
    }

    t = clock() - t; 
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
  
    printf("located %d particles in %f seconds\n", NP, time_taken); 
    return 0;
}
