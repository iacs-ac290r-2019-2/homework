#include <stdio.h>
#include <string.h>
#include <math.h>
#include "gaussian_quadrature.h"

#define SQRT3_INV 0.5773502691896257

// void gaussian_quadrature_1(float* xx, float* ww) {
//     static float xx_[1] = {0.};
//     static float ww_[1] = {2.};
//     static size_t n = 1;
//     memcpy(xx, xx_, n*sizeof(float));
//     memcpy(ww, ww_, n*sizeof(float));
//     return;
// }

// void gaussian_quadrature_2(float* xx, float* ww) {
//     static float xx_[2] = {-SQRT3_INV, SQRT3_INV};
//     static float ww_[2] = {1., 1.};
//     static size_t n = 2;
//     memcpy(xx, xx_, n*sizeof(float));
//     memcpy(ww, ww_, n*sizeof(float));
//     return;
// }

// void gaussian_quadrature(size_t n, float* xx, float* ww) {
//     if (n == 1)
//         gaussian_quadrature_1(xx, ww);
//     else if (n == 2)
//         gaussian_quadrature_2(xx, ww);
//     return;
// }

float N_local(size_t a, float xx) {
    if (a == 0) {
        return 0.5 * (1. - xx);
    } else if (a == 1) {
        return 0.5 * (1. + xx);
    } else {
        return 0.;
    }
}

float pN_local(size_t a, float xx) {
    if (a == 0) {
        return 0.5;
    } else if (a == 1) {
        return -0.5;
    } else {
        return 0.;
    }
}

// use gaussian quadrature to integrate the local Na*Nb
void NN_local_intmat(float* ddiag, float* ssub) {
    static size_t n = 2;
    static float xx[2] = {-SQRT3_INV, SQRT3_INV};
    static float ww[2] = {1., 1.};

    memset(ddiag, 0, 2*sizeof(float)); 
    memset(ssub,  0, 1*sizeof(float)); 
    for (size_t i=0; i<n; i++) {
        for (size_t a=0; a<2; a++) {
            ddiag[a] += ww[i] \
                        * N_local(a, xx[i]) \
                        * N_local(a, xx[i]);
        }
        for (size_t a=0; a<1; a++) {
            ssub[a]  += ww[i] \
                        * N_local(a, xx[i]) \
                        * N_local(a+1, xx[i]);
        }
    }
    
    return;
}

// use gaussian quadrature to integrate the local par(Na)*par(Nb)
void pNpN_local_intmat(float* ddiag, float* ssub) {
    static size_t n = 1;
    static float xx[1] = {0.};
    static float ww[1] = {2.};

    memset(ddiag, 0, 2*sizeof(float)); 
    memset(ssub,  0, 1*sizeof(float)); 
    for (size_t i=0; i<n; i++) {
        for (size_t a=0; a<2; a++) {
            ddiag[a] += ww[i] \
                        * pN_local(a, xx[i]) \
                        * pN_local(a, xx[i]);
        }
        for (size_t a=0; a<1; a++) {
            ssub[a]  += ww[i] \
                        * pN_local(a, xx[i]) \
                        * pN_local(a+1, xx[i]);
        }
    }
    return;
}

// assemble the NA*NB integral
void NN_global_intmat(size_t N, float* diag, float* sub) {
    static float ddiag[2];
    static float ssub[0];
    NN_local_intmat(ddiag, ssub);

    float coeff = 0.5 / N;
    diag[0] = coeff * ddiag[0];
    diag[N] = coeff * ddiag[1];
    for (size_t A=1; A<N; A++) {
        diag[A] = coeff * (ddiag[0] + ddiag[1]);
    }
    for (size_t A=0; A<N; A++) {
        sub[A]  = coeff * ssub[0];
    }
    return;
}

// assemble the par(NA)*par(NB) integral
void pNpN_global_intmat(size_t N, float* diag, float* sub) {
    static float ddiag[2];
    static float ssub[0];
    pNpN_local_intmat(ddiag, ssub);

    float coeff = 2. * N;
    diag[0] = coeff * ddiag[0];
    diag[N] = coeff * ddiag[1];
    for (size_t A=1; A<N; A++) {
        diag[A] = coeff * (ddiag[0] + ddiag[1]);
    }
    for (size_t A=0; A<N; A++) {
        sub[A]  = coeff * ssub[0];
    }
    return;
}

