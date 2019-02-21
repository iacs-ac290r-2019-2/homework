#include "gaussian_quadrature.h"

void gaussian_quadrature_1(float* xx, float* ww) {
    static float xx_[1] = {0.};
    static float ww_[1] = {2.};
    static size_t n = 1;
    memcpy(xx, xx_, n*sizeof(float));
    memcpy(ww, ww_, n*sizeof(float));
    return;
}

void gaussian_quadrature_2(float* xx, float* ww) {
    static float xx_[2] = {-SQRT3_INV, SQRT3_INV};
    static float ww_[2] = {1., 1.};
    static size_t n = 2;
    memcpy(xx, xx_, n*sizeof(float));
    memcpy(ww, ww_, n*sizeof(float));
    return;
}

void gaussian_quadrature(size_t n, float* xx, float* ww) {
    
    if (n == 1)
        gaussian_quadrature_1(xx, ww);
    else if (n == 2)
        gaussian_quadrature_2(xx, ww);
    return;
}

float local_N(float xx, size_t a) {
    if (a == 0) {
        return 0.5 * (1. + xx);
    } else if (a == 1) {
        return 0.5 * (1. - xx);
    } else {
        return 0.;
    }
}

float local_pN(float xx, size_t a) {
    if (a == 0) {
        return 0.5;
    } else if (a == 1) {
        return -0.5;
    } else {
        return 0.;
    }
}

float int_local_NN(size_t a, size_t b) {
    static size_t n = 2;
    float* xx = new float[n];
    float* ww = new float[n];
    gaussian_quadrature(n, xx, ww);
    float s = 0.;
    for (int i=0; i<n; i++) {
        s += local_N(xx[i], a) * local_N(xx[i], b) * ww[i];
    }
    delete [] xx;
    delete [] ww;
    return s;
}

float int_local_pNpN(size_t a, size_t b) {
    static size_t n = 1;
    float* xx = new float[n];
    float* ww = new float[n];
    gaussian_quadrature(n, xx, ww);
    float s = 0.;
    for (int i=0; i<0; i++) {
        s += local_pN(xx[i], a) * local_pN(xx[i], b) * ww[i];
    }
    delete [] xx;
    delete [] ww;
    return s;
}

