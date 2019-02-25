#include <stdio.h>

#include "tridiag_solver.h"
#include "construct_matrix.h"

#define N 100  // Number of nodes

int main(int argc, char* argv[])
{
    // -----------------------
    // Variable initialization
    // -----------------------

    // Problem specification. Can change by users.
    float h_bound = 1.0;
    float g_bound = 1.0;
    // f(x) is hard-coded elsewhere
    float f[N+1]; 
    for (int i=0; i<N+1; i++) {
        f[i] = 1.;
    }

    // values that will be computed later
    float K_diag[N]; // major diagonal
    float K_subdiag[N-1]; // sub diagonal
    float F[N]; // right hand side
    float u[N+1]; u[N]=g_bound; // solution on nodes

    // -----------------------
    // Construct global matrices K and F
    // -----------------------

    // Construct tri-diagonal matrix K
    // Input: N
    // Output: K_diag, K_subdiag
    // compute_K(N, K_diag, K_subdiag);

    // Construct right-hand side vector F
    // Input: N, h_bound, g_bound
    // Output: F
    // compute_F(N, h_bound, g_bound, F);

    compute_KF(N, h_bound, g_bound, f, K_diag, K_subdiag, F);

    // -----------------------
    // Construct matrices
    // -----------------------

    // Solve linear system
    // Input: N, K_diag, K_subdiag, F
    // Output: u
    tridiag_solver(N, K_diag, K_subdiag, F, u);
    
    // -----------------------
    // Show result
    // -----------------------

    printf("=================== \n");
    printf("Problem specification \n");
    printf("=================== \n");

    printf("Problem size N: %i \n", N);
    printf("Boundary condition: h_bound=%2.2f, g_bound= %2.2f \n", h_bound, g_bound);
    
    printf("=================== \n");
    printf("Value of matrices \n");
    printf("=================== \n");

    printf("K_diag: \n");
    for(int i = 0; i < N; i++) {
            printf("%f ", K_diag[i]);
        }
    printf("\n");

    printf("K_subdiag: \n");
    for(int i = 0; i < N-1; i++) {
            printf("%f ", K_subdiag[i]);
        }
    printf("\n");


    printf("RHS vector F: \n");
    for(int i = 0; i < N; i++) {
            printf("%f ", F[i]);
        }
    printf("\n");


    printf("Solution u: \n");
    for(int i = 0; i < N+1; i++) {
            printf("%f ", u[i]);
        }
    printf("\n");

    // TODO: write to file and plot in Python

}

