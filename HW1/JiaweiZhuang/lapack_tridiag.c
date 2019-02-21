#include <stdio.h>
#include <lapacke.h>

#define N 5

int main(int argc, char* argv[])
{
    // Initialize array
    // TODO: use input data
    double a[N-1] = {-1, -1, -1, -1}; // lower diagonal
    double b[N] = {2, 2, 2, 2, 2}; // major diagonal
    double c[N-1] = {-1, -1, -1, -1}; // upper diagonal

    double f[N] = {1, 2, 3, 4, 5}; // right hand side, will be overwritten by solution

    // Solve linear system with LAPACK
    lapack_int info,n,nrhs,ldb;
    n = N;
    nrhs = 1;  // only a single column on RHS
    ldb = N;

    info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,n,nrhs,a,b,c,f,ldb);

    // print result to check
    for(int i = 0; i < N; i++) { 
            printf("%f ", f[i]);
        }
    printf("\n");

}
