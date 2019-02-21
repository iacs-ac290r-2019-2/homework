#include <stdio.h>
#include <string.h>
#include <math.h>

const float SQRT3 = sqrt(3.);
const float SQRT3_INV = 1./SQRT3;

void gaussian_quadrature_1(float* x, float* w) {
    static const xx[1] = {0.};
    static const ww[1] = {2.};
    static size_t n = 1;
    memcpy(xx, x, n);
    memcpy(ww, w, n);
    return;
}

void gaussian_quadrature_2(float* x, float* w) {
    static const xx[2] = {-SQRT3_INV, SQRT3_INV};
    static const ww[2] = {1., 1.};
    static size_t n = 2;
    memcpy(xx, x, n);
    memcpy(ww, w, n);
    return;
}

void gaussian_quadrature(unsigned n, float* x, float* w) {
    /* 
    returns the points x and weights w for given order n
    the domain is always within [-1, 1]

    USAGE:
      suppose you have a funtion f(x) on [-1,1]
      and would like to calculate the integral to the 2nd order
      the example code would be:

        int n = 2;
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
        gaussian_quadrature_1(x, w);
    else if (n == 2)
        gaussian_quadrature_2(x, w);
    return;
}

