#include <math.h>

// Equals function with tolerance (to avoid floating point)
int _soft_eq(double a, double b, double epsilon) {
  return (a == b) || (a <= b + epsilon && a >= b - epsilon);
}

// Bilinear interpolation
double bilinear(const double* normalized_coord, double val1, const double val2, double val3, double val4) {
  // Bilinear weights
  double V[4];
  // Set up weights
  V[0] = val1;
  V[1] = -val1 + val2;
  V[2] = -val1 + val4;
  V[3] = val1 - val2 + val3 - val4;
  // Evaluate
  return V[0] + V[1] * normalized_coord[0] + V[2] * normalized_coord[1] + V[3] * normalized_coord[0] * normalized_coord[1];
}

// Inverse bilinear interpolation
// Given xy coordinates and the vertices of a quadrilateral, returns the normalized coordinates of xy in [0,1]x[0,1]
int inverse_bilinear(const double* coord, const double* vec1, const double* vec2, const double* vec3, const double* vec4, double* normalized_coord) {
  // Bilinear weights
  double V[4], W[4];
  // Polynomial variables
  double a, b, c, det;
  // Normalized coordinates
  double xi, eta;
  // Set up weights
  V[0] = vec1[0];
  V[1] = -vec1[0] + vec2[0];
  V[2] = -vec1[0] + vec4[0];
  V[3] = vec1[0] - vec2[0] + vec3[0] - vec4[0];
  W[0] = vec1[1];
  W[1] = -vec1[1] + vec2[1];
  W[2] = -vec1[1] + vec4[1];
  W[3] = vec1[1] - vec2[1] + vec3[1] - vec4[1];
  // Calculate polynomial coefficients
  a = V[3] * W[2] - V[2] * W[3];
  b = V[3] * W[0] - V[0] * W[3] + V[1] * W[2] - V[2] * W[1] + coord[0] * W[3] - coord[1] * V[3];
  c = V[1] * W[0] - V[0] * W[1] + coord[0] * W[1] - coord[1] * V[1];

  // Solve for xi, eta in [0,1]x[0,1]

  if (_soft_eq(a, 0, 1e-10)) {
    // Linear case
    xi = -c / b;
  } else {
    // Quadratic case
    det = sqrt(b * b - 4 * a * c);
    xi = (-b + det) / (2 * a);
  }

  eta = (coord[0] - V[0] - V[2] * xi) / (V[1] + V[3] * xi);

  // Assign normalized coordinates (eta - y, xi - x)
  normalized_coord[0] = eta;
  normalized_coord[1] = xi;

  return 0;
}

// Inverse bilinear interpolation
// Given xyz coordinates and the vertices of a cell, returns the normalized coordinates of xyz in [0,1]x[0,1]x[0,1]
// It is assumed that xy coordinates remain constant in the z direction
// Not yet completed
int inverse_trilinear(const double* coord, const double* v1, const double* v2, const double* v3, const double* v4, const double* v5, const double* v6, const double* v7, const double* v8, double* normalized_coord) {
  double z0, z1;

  inverse_bilinear(coord, v1, v2, v3, v4, normalized_coord);

  z0 = bilinear(normalized_coord, v1[2], v2[2], v3[2], v4[2]);
  z1 = bilinear(normalized_coord, v5[2], v6[2], v7[2], v8[2]);

  normalized_coord[2] = (coord[2] - z0) / (z1 - z0);

  return 0;
}
