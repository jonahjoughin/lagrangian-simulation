#ifndef INTERPOLATE_H_INCLUDED
#define INTERPOLATE_H_INCLUDED

double bilinear(const double* normalized_coord, double val1, const double val2, double val3, double val4);

int inverse_bilinear(const double* coord, const double* vec1, const double* vec2, const double* vec3, const double* vec4, double* normalized_coord);

#endif