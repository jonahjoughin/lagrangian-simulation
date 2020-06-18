#define NX_PSI 163
#define NY_PSI 129
#define NX_U NX_PSI
#define NY_U (NY_PSI + 1)
#define NX_V (NX_PSI + 1)
#define NY_V NY_PSI
#define NX_RHO (NX_PSI + 1)
#define NY_RHO (NY_PSI + 1)
#define NZ_RHO 1

struct Mesh {
  // Grid Coordinates
  double x_psi[NY_PSI][NX_PSI];
  double y_psi[NY_PSI][NX_PSI];

  double y_rho[NY_RHO][NX_RHO];
  double x_rho[NY_RHO][NX_RHO];
  double mask_rho[NY_RHO][NX_RHO];
  uint64_t z_curve_idx_rho[NZ_RHO][NY_RHO][NX_RHO];

  double y_u[NY_U][NX_U];
  double x_u[NY_U][NX_U];

  double y_v[NY_V][NX_V];
  double x_v[NY_V][NX_V];
  // UV values
  double u[NY_U][NX_U];
  double v[NY_V][NX_V];
};

int build_mesh(int grd_ncid, int avg_ncid, struct Mesh *mesh);