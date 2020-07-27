#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED

#include <stdint.h>

#define NX_RHO 164
#define NY_RHO 130
#define NZ_RHO 12

#define NX_U (NX_RHO - 1)
#define NY_U (NY_RHO)
#define NZ_U (NZ_RHO)

#define NX_V (NX_RHO)
#define NY_V (NY_RHO - 1)
#define NZ_V (NZ_RHO)

#define NX_W (NX_RHO)
#define NY_W (NY_RHO)
#define NZ_W (NZ_RHO + 1)

#define NX_PSI (NX_RHO - 1)
#define NY_PSI (NY_RHO - 1)

#define SQUARE(x) ({ typeof (x) _x = (x); _x * _x; })

typedef struct BoundarySegment {
  char axis;
  int x_idx, y_idx, z_idx;
  int orientation;
  double flux;
} BoundarySegment;

typedef struct Boundary {
  struct BoundarySegment *segments;
  int segments_count;
  double max_flux;
} Boundary;

typedef struct Mesh {
  // Coordinate arrays
  double x_psi[NY_PSI][NX_PSI];
  double y_psi[NY_PSI][NX_PSI];
  double y_u[NY_U][NX_U];
  double x_u[NY_U][NX_U];
  double y_v[NY_V][NX_V];
  double x_v[NY_V][NX_V];
  double z_w[NZ_W][NY_W][NX_W];

  double mask_rho[NY_RHO][NX_RHO];

  struct Boundary open_boundary;

  // u,v,w arrays
  double u[NZ_U][NY_U][NX_U];
  double v[NZ_V][NY_V][NX_V];
  double w[NZ_W - 1][NY_W][NX_W];

  uint64_t z_curve_idx_rho[NZ_RHO][NY_RHO][NX_RHO];
} Mesh;

int build_test_mesh(struct Mesh *mesh);
int build_mesh(int grd_ncid, int avg_ncid, struct Mesh *mesh);

int initialize_mesh(int grd_ncid, struct Mesh *mesh);
int update_mesh(int avg_ncid, struct Mesh *mesh, int t_idx);

#endif