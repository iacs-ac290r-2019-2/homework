#include <stdio.h>

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
    // Some fake F!
    // TODO: Use correctly computed F!

    for(int j=0; j<N-1; ++j){
        F[j] = j * h_bound + g_bound;
	}
}