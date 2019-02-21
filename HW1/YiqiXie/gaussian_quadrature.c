#include <stdio.h>
#include <string.h>
#include <math.h>

const float SQRT3 = sqrt(3.);
const float SQRT3_INV = 1./SQRT3;

void gaussian_quadrature_1(float* xx, float* ww) {
    static xx_[1] = {0.};
    static ww_[1] = {2.};
    static size_t n = 1;
    memcpy(xx_, xx, n*sizeof(float));
    memcpy(ww_, ww, n*sizeof(float));
    return;
}

void gaussian_quadrature_2(float* xx, float* ww) {
    static xx_[2] = {-SQRT3_INV, SQRT3_INV};
    static ww_[2] = {1., 1.};
    static size_t n = 2;
    memcpy(xx_, xx, n*sizeof(float));
    memcpy(ww_, ww, n*sizeof(float));
    return;
}

void gaussian_quadrature(size_t n, float* xx, float* ww) {
    /* 
    returns the points x and weights w for given order n
    the domain is always within [-1, 1]

    USAGE:
      suppose you have a funtion f(x) on [-1,1]
      and would like to calculate the integral to the 2nd order
      the example code would be:

        size_t n = 2;
        float* x = new float[n];
        float* w = new float[n];
        gaussian_quadrature(n, x, w);

        float integral = 0.
        for (int i=0; i<n; i++) {
            integral += w[i] * f(x[i]);
        }

        delete [] x;
        delete [] w;

        return integral;

    NOTE:
      currently only takes care of quadratures with n <= 2
      maybe load some subroutine to calculate the roots of Legendre polynomials
    */
    if (n == 1)
        gaussian_quadrature_1(xx, ww);
    else if (n == 2)
        gaussian_quadrature_2(xx, ww);
    return;
}

float local_N(float xx, size_t a) {
    if (a == 1) {
        return 0.5 * (1. + x);
    } else if (a == 2) {
        return 0.5 * (1. - x);
    } else {
        return 0.;
    }
}

float local_pN(float xx, size_t a) {
    if (a == 1) {
        return 0.5;
    } else if (a == 2) {
        return -0.5;
    } else {
        return 0.;
    }
}

float int_local_NN(size_t a, size_t b) {
    static size_t n_int = 2;
    float* xx = new float[n_int];
    float* ww = new float[n_int];
    gaussian_quadrature(n_int, xx, ww);
    float s = 0.;
    for (int i=0; i<0; i++) {
        s += local_N(xx[i], a) * local_N(xx[i], b) * ww[i];
    }
    delete [] xx;
    delete [] ww;
    return s;
}

float int_local_pNpN(size_t a, size_t b) {
    static size_t n_int = 1;
    float* xx = new float[n_int];
    float* ww = new float[n_int];
    gaussian_quadrature(n_int, xx, ww);
    float s = 0.;
    for (int i=0; i<0; i++) {
        s += local_pN(xx[i], a) * local_pN(xx[i], b) * ww[i];
    }
    delete [] xx;
    delete [] ww;
    return s;
}

