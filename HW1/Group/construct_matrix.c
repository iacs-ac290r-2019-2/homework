#include <stdio.h>
#include <stdlib.h> 
#include <string.h> 
#include "gaussian_quadrature.h"

// float f(float x) { 
//     return 1.; 
// }

// void compute_K(int N, float *K_diag, float *K_subdiag)
// {
//     float h = 1.0/((float)N);

// 	for(int j=0; j<N-1; ++j){
// 		K_diag[j] = 2/h;
// 		K_subdiag[j] = -1/h;
// 	}
// 	K_diag[N-1] = 2/h;

//     return;

// }

// void compute_F(int N, float h_bound, float g_bound, float *F)
// {
//     // Some fake F!
//     // TODO: Use correctly computed F!

//     // for(int j=0; j<N; ++j){
//     //     F[j] = j * h_bound + g_bound;
// 	// }

//     float* diag_NN = malloc((N+1)*sizeof(float));
//     float* sub_NN  = malloc((N)*sizeof(float));
//     NN_global_intmat(N, diag_NN, sub_NN);

//     memset(F, 0, N*sizeof(float));
//     F[0] += h_bound;
//     for (size_t A=0; A<N; A++) {
//         F[A] += diag_NN[A] * f((float)A/N);
//     }
//     for (size_t A=0; A<N-1; A++) {
//         F[A]   += sub_NN[A] * f((float)(A+1)/N);
//         F[A+1] += sub_NN[A] * f((float)(A)/N);
//     }
//     // F[N-1] += -g_bound * (2*N * pNpN[0][1]);
//     float ssub;
//     pNpN_local_intmat(NULL, &ssub);
//     F[N-1] += -g_bound * (2.*N * ssub);

//     free(diag_NN);
//     free(sub_NN);

//     return;
// }

void compute_KF(int N, float h_bound, float g_bound, float* f, float* K_diag, float* K_sub, float* F) {
    float* NN_diag   = malloc((N+1)*sizeof(float));
    float* NN_sub    = malloc( N   *sizeof(float));
    float* pNpN_diag = malloc((N+1)*sizeof(float));
    float* pNpN_sub  = malloc( N   *sizeof(float));

    NN_global_intmat  (N, NN_diag,   NN_sub  );
    pNpN_global_intmat(N, pNpN_diag, pNpN_sub);

    memcpy(K_diag, pNpN_diag,  N   *sizeof(float));
    memcpy(K_sub,  pNpN_sub,  (N-1)*sizeof(float));

    memset(F, 0, N*sizeof(float));
    for (size_t A=1; A<N; A++) {
        F[A] += NN_sub[A] * f[A+1] + NN_diag[A] * f[A] + NN_sub[A-1] * f[A-1];
    }
    F[0]   +=  h_bound + NN_sub[0] * f[1] + NN_diag[0] * f[0];
    F[N-1] += -g_bound * pNpN_sub[N-1];

    free(NN_diag);
    free(NN_sub);
    free(pNpN_diag);
    free(pNpN_sub);
    return;
}