#include <main.h>
#include <stdio.h>
#include <string.h>

void _assign_vectors(const float *p_a, const float *p_b, const float *p_c, float *vec_a, float *vec_b) {
    for (int i = 0; i < 2; i++) {
            vec_a[i] = p_a[i] - p_c[i];
            vec_b[i] = p_b[i] - p_a[i];
        }
}

int _cross_product_sign(const float* vec1, const float* vec2) {
    float val = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    if (val == 0) {
        return 0;
    } else {
        return val > 0 ? 1 : -1;
    }
}

int locate_particle(float *coord_mesh, int ny, int nx, float *particle_coord, int *idx) {

    float (*mesh)[ny+1][nx+1] = (float (*)[ny+1][nx+1]) coord_mesh;
    float p1[2], p2[2], p3[2], p4[2];
    float vec1[2], vec2[2];
    int _idx[2] = {0, 0};

    while (1) {
        if (_idx[0] < 0 || _idx[1] < 0 || _idx[0] >= ny || _idx[1] >= nx) {
            return -1;
        }

        // Assign grid vertices
        p1[0] = mesh[0][_idx[0]][_idx[1]];
        p1[1] = mesh[1][_idx[0]][_idx[1]];
        p2[0] = mesh[0][_idx[0]][_idx[1]+1];
        p2[1] = mesh[1][_idx[0]][_idx[1]+1];
        p3[0] = mesh[0][_idx[0]+1][_idx[1]+1];
        p3[1] = mesh[1][_idx[0]+1][_idx[1]+1];
        p4[0] = mesh[0][_idx[0]+1][_idx[1]];
        p4[1] = mesh[1][_idx[0]+1][_idx[1]];

        // Check cross product of vectors
        // If sign of one of the below is negative, particle is outside current box
        // See (https://en.wikipedia.org/wiki/Stencil_jumping) for details

        _assign_vectors(p1, p2, particle_coord, vec1, vec2);
        if (_cross_product_sign(vec2, vec1) == -1) {
            _idx[0] -= 1;
            continue;
        }

        _assign_vectors(p2, p3, particle_coord, vec1, vec2);
        if (_cross_product_sign(vec2, vec1) == -1) {
            _idx[1] += 1;
            continue;
        }

        _assign_vectors(p3, p4, particle_coord, vec1, vec2);
        if (_cross_product_sign(vec2, vec1) == -1) {
            _idx[0] += 1;
            continue;
        }

        _assign_vectors(p4, p1, particle_coord, vec1, vec2);
        if (_cross_product_sign(vec2, vec1) == -1) {
            _idx[1] -= 1;
            continue;
        }

        // Copy correct index and return
        memcpy (idx, _idx, sizeof(_idx));
        return 0;
    }

}

