#include <stdio.h>
#include "tridiag_solver.h"

#define N 5

int main(int argc, char* argv[])
{
    // Initialize array
    // TODO: use input data
    float a[N-1] = {-1, -1, -1, -1}; // lower diagonal
    float b[N] = {2, 2, 2, 2, 2}; // major diagonal
    float c[N-1] = {-1, -1, -1, -1}; // upper diagonal

    float f[N] = {1, 2, 3, 4, 5}; // right hand side, will be overwritten by solution

    tridiag_solver(N, a, b, c, f);

    // print result to check
    for(int i = 0; i < N; i++) {
            printf("%f ", f[i]);
        }
    printf("\n");

    

}

