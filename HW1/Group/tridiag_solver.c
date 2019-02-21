#include <lapacke.h>

// A tiny wrapper of LAPACK for single-column RHS
void tridiag_solver(int N, float *K_diag, float *K_subdiag, float *F, float *u)
{
    // LAPACK parameters
    lapack_int info,n,nrhs,ldb;
    n = N;
    nrhs = 1;  // only a single column on RHS
    ldb = N;

    // Temporary array
    float dl[N-1];
    float d[N];
    float du[N-1];

    // LAPACK routine will overwrite RHS vector by the solution
    // So we put RHS into solution vector u before calling LAPACK
    // It ensures that F is not modified
    for(int i = 0; i < N; i++) {
            u[i] = F[i];
        }

    // LAPACK will overwrite the LHS matrix as well
    // Use Temporary array to prevent K from being changed,
    // as we might want to print K later
    for(int i = 0; i < N; i++) {
        d[i] = K_diag[i];
    }

    for(int i = 0; i < N-1; i++) {
        du[i] = K_subdiag[i];
        dl[i] = K_subdiag[i];
    }

    // Standard LAPACK signature
    info = LAPACKE_sgtsv(LAPACK_COL_MAJOR,n,nrhs,dl,d,du,u,ldb);

}
