#include <lapacke.h>

// A tiny wrapper of LAPACK for single-column RHS
void tridiag_solver(int N, float* dl, float* d, float *du, float *f)
{
    // Solve linear system with LAPACK
    lapack_int info,n,nrhs,ldb;
    n = N;
    nrhs = 1;  // only a single column on RHS
    ldb = N;

    info = LAPACKE_sgtsv(LAPACK_COL_MAJOR,n,nrhs,dl,d,du,f,ldb);

}
