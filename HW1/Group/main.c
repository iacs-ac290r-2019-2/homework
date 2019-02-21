#include <stdio.h>

#define N 10  // Number of nodes

int main(int argc, char* argv[])
{
    // -----------------------
    // Variable initialization
    // -----------------------

    // Problem specification. Can change by users.
    float h_bound = 1.0;
    float g_bound = 0.0;
    // f(x) is hard-coded elsewhere

    // values that will be computed later
    float K_diag[N]; // major diagonal
    float K_subdiag[N-1]; // sub diagonal
    float F[N]; // right hand side
    float u[N]; // solution on nodes

    // -----------------------
    // Construct global matrices K and F
    // -----------------------

    // Construct tri-diagonal matrix K
    // Input: N
    // Output: K_diag, K_subdiag
    compute_K(N, K_diag, K_subdiag);

    // Construct right-hand side vector F
    // Input: N, h_bound, g_bound
    // Output: F
    compute_F(N, h_bound, g_bound, F);

    // -----------------------
    // Construct matrices
    // -----------------------

    // Solve linear system
    // Input: N, K_diag, K_subdiag, F
    // Output: u
    tridiag_solver(N, K_diag, K_subdiag, F, u);
    
    // -----------------------
    // Output result
    // -----------------------
    
    // Print result to check
    for(int i = 0; i < N; i++) {
            printf("%f ", u[i]);
        }
    printf("\n");

    // TODO: write to file and plot in Python

}

