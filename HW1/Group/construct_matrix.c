#include <stdio.h>
#include <string.h> 
#include "gaussian_quadrature.h"

float f(float x) {
    return 0.;
}

void compute_K(int N, float *K_diag, float *K_subdiag)
{
    float h = 1.0/((float)N);

	for(int j=0; j<N-1; ++j){
		K_diag[j] = 2/h;
		K_subdiag[j] = -1/h;
	}
	K_diag[N-1] = 2/h;

}

void compute_F(int N, float h_bound, float g_bound, float *F)
{
 //    // Some fake F!
 //    // TODO: Use correctly computed F!

 //    for(int j=0; j<N; ++j){
 //        F[j] = j * h_bound + g_bound;
	// }
    float NN[2][2];
    float pNpN[2][2];
    for (int a=0; a<2; a++) {
        for (int b=0; b<2; b++) {
            NN[a][b] = int_local_NN(a,b);
            pNpN[a][b] = int_local_pNpN(a,b);
            printf("%f ", NN[a][b]);
            printf("%f ", pNpN[a][b]);
            printf("\n");
        }
    }

    memset(F, 0, N*sizeof(float));
    for (int A=0; A<N; A++) {
        for (int B=A-1; B<=A+1; B++) {
            float x = (B+1) * 1./N;
            if (B == A) {
                F[A] += f(x) * (2*N * NN[0][0] + 2*N * NN[1][1]);
            } else {
                F[A] += f(x) * (2*N * NN[0][1]);
            }
        }
    }
    F[0] += h_bound;
    F[N-1] += -g_bound * (2*N * pNpN[0][1]);

    return;
}