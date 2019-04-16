/* Daniel Willen, 2019
 *
 * Solve the transient heat conduction problem with homogeneous Dirichlet
 *  boundary conditions:
 *
 *    u(x={0,L}) = u(y={0,L}) = 0
 *
 *  and initial condition:
 *
 *    u(x,y,0) = sin(x) * sin(y)
 *
 *  on the domain 0 <= x,y <= L, with L = pi.
 *
 * This program solves the above problem on a single GPU with the Jacobi method.
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda.h>

#include <thrust/device_ptr.h>
#include <thrust/reduce.h>

#define PI 3.14159265358979323846
#define MAX_THREADS_DIM 18        // Note that this depends on the hardware

/* Note on the structure of this file:
 *  - Cuda device constant memory declarations are at the top
 *  - Functions definitions are in the middle. Functions include:
 *  - - parse_cmdline: Read command-line arguments for domain size
 *  - - jacobi_solver: Advance the soln to the next time step using Jacobi
 *  - - check_error:   Calculate the error b/t the numeric and analytic solns
 *  - The `main' function is at the bottom
 *
 *  Note that it is good practice to use header files and break functions out
 *   into separate files. This has not been done here for simplicity.
 */

/*** Auxiliary Functions ***/

/* Read the command line inputs */
// - argv[0] is the program name
// - argv[1] is the first input (number of points)
int parse_cmdline(int argc, char *argv[]) {
  int nx;
  if (argc >= 2) {
    nx = atoi(argv[1]); // Number of grid points
    printf("Grid is %d by %d\n\n", nx, nx);
  } else {
    printf("Input error. Run like: \n\n");
    printf("  $ ./parallel.c n\n\n");
    printf("  where n is the number of grid cells in one dimension\n");
    exit(EXIT_FAILURE);
  }
  return nx;
}

/*** GPU Constants ***/
__constant__ int _nx;
__constant__ int _ny;
__constant__ double _Lx;
__constant__ double _Ly;
__constant__ double _dx;
__constant__ double _dy;
__constant__ double _dt;
__constant__ double _D;
__constant__ double _pref;

/*******************************************************************************
 * Step IV: Launch the GPU kernel to advance to the next time step with the    *
 *          Jacobi method here.                                                *
 ******************************************************************************/

__global__ void jacobi_solver(double* u, double* u_new) {

  // int ti = blockIdx.x*blockDim.x + threadIdx.x;
  // int tj = blockIdx.y*blockDim.y + threadIdx.y;

  // if ((ti >= 1 && ti < (_nx-1)) && 
  //     (tj >= 1 && tj < (_ny-1))) {
  //   u_new[ti + tj*_nx] = 
  //       u[ti + tj*_nx] + _pref * (
  //         u[(ti+1) + tj*_nx] + 
  //         u[(ti-1) + tj*_nx] + 
  //         u[ti + (tj+1)*_nx] + 
  //         u[ti + (tj-1)*_nx] + 
  //         u[ti + tj*_nx] * (-4)
  //       );
  // }

  __shared__ double s_u[MAX_THREADS_DIM*MAX_THREADS_DIM];

  int si = threadIdx.x;
  int sj = threadIdx.y;
  int ti = blockIdx.x*(blockDim.x-2) + threadIdx.x;
  int tj = blockIdx.y*(blockDim.y-2) + threadIdx.y;
  
  if (ti < _nx && tj < _ny) {
    s_u[si + sj*blockDim.x] = u[ti + tj*_nx];
  }

  __syncthreads();

  if ((ti >= 1 && ti < (_nx-1)) && 
      (tj >= 1 && tj < (_ny-1)) &&
      (si > 0 && si < (blockDim.x-1)) &&
      (sj > 0 && sj < (blockDim.y-1))) {
    u_new[ti + tj*_nx] = 
        s_u[si + sj*blockDim.x] + _pref * (
          s_u[(si+1) + sj*blockDim.x] + 
          s_u[(si-1) + sj*blockDim.x] + 
          s_u[si + (sj+1)*blockDim.x] + 
          s_u[si + (sj-1)*blockDim.x] + 
          s_u[si + sj*blockDim.x] * (-4)
        );
  }

  return;
}

/******************************************************************************
 * Step V: Launch the GPU kernel to calculate the error at each grid point    *
 *         here.                                                              *
 *****************************************************************************/

__global__ void check_error(double* u, double* error, double time) {

  // int ti = blockIdx.x*blockDim.x + threadIdx.x;
  // int tj = blockIdx.y*blockDim.y + threadIdx.y;

  // if ((ti >= 1 && ti < (_nx-1)) && 
  //     (tj >= 1 && tj < (_ny-1))) {
  //   error[ti + tj*_nx] = u[ti + tj*_nx] / (sin(ti*_dx) * sin(tj*_dy) * exp(-2*_D*time)) - 1;
  // }

  int ti = blockIdx.x*(blockDim.x-2) + threadIdx.x;
  int tj = blockIdx.y*(blockDim.y-2) + threadIdx.y;

  if ((ti >= 1 && ti < (_nx-1)) && 
      (tj >= 1 && tj < (_ny-1)) &&
      (threadIdx.x > 0 && threadIdx.x < (blockDim.x-1)) &&
      (threadIdx.y > 0 && threadIdx.y < (blockDim.y-1))) {
    error[ti + tj*_nx] = u[ti + tj*_nx] / (sin(ti*_dx) * sin(tj*_dy) * exp(-2*_D*time)) - 1;
  }

  return;
}

/*** Main Function ***/
int main(int argc, char *argv[])
{
  /* Variable declaration */
  double Lx = PI;           // Domain length in x-direction
  double Ly = PI;           // Domain length in y-direction
  double D = 1.;            // Diffusion constant

  int nx, ny;               // Grid points (grid cells + 1)
  double dx, dy;            // Grid spacing
  double dt;                // Time step size
  double sim_time;          // Length of sim time, arbitrary for simplicity
  double pref;              // Pre-factor in the Jacobi method

  double error = 0.;        // Mean percent-difference at each grid point
  error = error;            // To prevent compiler warning

  /* Parse command-line for problem size */
  nx = parse_cmdline(argc, argv);
  ny = nx;                  // Assume a square grid

  /* Initialize variables */
  dx = Lx / (nx - 1);       // Cell width in x-direction
  dy = Ly / (ny - 1);       // Cell width in y-direction
  dt = 0.25*dx*dy/D;        // Limited by diffusive stability
  sim_time = Lx*Ly/D;       // Arbitrary simulation length
  pref = D*dt/(dx*dx);      // Jacobi pre-factor

  printf("Parameters\n");
  printf("---------------------------\n");
  printf("Lx = %.5lf\n", Lx); 
  printf("Lx = %.5lf\n", Ly); 
  printf("T  = %.5lf\n", sim_time); 
  printf("D  = %.5lf\n", D);
  printf("nx = %d\n", nx);
  printf("ny = %d\n", nx);
  printf("dx = %.5lf\n", dx);
  printf("dy = %.5lf\n", dy);
  printf("dt = %.5lf\n", dt);
  printf("\n");

  cudaMemcpyToSymbol(_nx, &nx, sizeof(int));
  cudaMemcpyToSymbol(_ny, &ny, sizeof(int));
  cudaMemcpyToSymbol(_Lx, &Lx, sizeof(double));
  cudaMemcpyToSymbol(_Ly, &Ly, sizeof(double));
  cudaMemcpyToSymbol(_dx, &dx, sizeof(double));
  cudaMemcpyToSymbol(_dy, &dy, sizeof(double));
  cudaMemcpyToSymbol(_dt, &dt, sizeof(double));
  cudaMemcpyToSymbol(_D,  &D,  sizeof(double));
  cudaMemcpyToSymbol(_pref, &pref, sizeof(double));

  /*****************************************************************************
   * Step I: Declare, allocate, and initialize memory for the field variable   *
   *         u on the CPU.                                                     *
   ****************************************************************************/

  double *u;
  u = (double*) malloc(nx*ny * sizeof(double));
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      u[i+j*nx] = sin(i*dx) * sin(j*dy);
    }
  }

  /*****************************************************************************
   * Step II: Declare and allocate GPU memory for _u, _u_new, and _error. Copy *
   *          the initial condition to the GPU.                                *
   ****************************************************************************/

  double *_u, *_u_new, *_error;
  cudaMalloc(&_u, nx*ny * sizeof(double));
  cudaMalloc(&_u_new, nx*ny * sizeof(double));
  cudaMalloc(&_error, nx*ny * sizeof(double));
  cudaMemcpy(_u, u, nx*ny * sizeof(double), cudaMemcpyHostToDevice);

  // Set the new soln and error to 0
  cudaMemset(_u_new, 0., nx*ny * sizeof(double));
  cudaMemset(_error, 0., nx*ny * sizeof(double));

  // Create thrust pointers to device memory for error calculation
  thrust::device_ptr<double> t_error(_error);

  /*****************************************************************************
   * Step III: Set up the kernel execution configuration for the domain based  *
   *           on the input domain size and the MAX_THREADS_DIM variable.      *
   ****************************************************************************/

  int threads_x = MAX_THREADS_DIM;
  int threads_y = MAX_THREADS_DIM;

  // int blocks_x = (int) ceil((double) nx / (double) threads_x);
  // int blocks_y = (int) ceil((double) ny / (double) threads_y);
  int blocks_x = (int) ceil((double) nx / (double) (threads_x - 2));
  int blocks_y = (int) ceil((double) ny / (double) (threads_y - 2));

  dim3 dim_blocks(threads_x, threads_y);
  dim3 num_blocks(blocks_x, blocks_y);

  printf("Parallelization\n");
  printf("---------------------------\n");
  printf("MAX_THREADS_DIM = %d\n", MAX_THREADS_DIM);
  printf("threads_x = %d\n", threads_x); 
  printf("threads_y = %d\n", threads_y); 
  printf("blocks_x  = %d\n", blocks_x); 
  printf("blocks_y  = %d\n", blocks_y); 
  printf("\n");

  /***************************/
  /* Main Time-Stepping Loop */
  /***************************/

  for (double time = 0.; time <= sim_time; time += dt) {
    /***************************************************************************
     * Step IV: Launch the GPU kernel to advance to the next time step with    *
     *          the Jacobi method here.                                        *
     **************************************************************************/

    jacobi_solver<<<num_blocks, dim_blocks>>>(_u, _u_new);

    /***************************************************************************
     * Step V: Launch the GPU kernel to calculate the error at each grid point *
     *         here.                                                           *
     **************************************************************************/

    check_error<<<num_blocks, dim_blocks>>>(_u, _error, time);

    // Use thrust to do a parallel reduction on the error
    error = thrust::reduce(t_error, t_error + nx*ny, 0., thrust::plus<double>());
    printf("Error at t* = %.5lf is %e\n", time*D/(Lx*Lx), error/(nx*ny));

    // Copy new soln to old. This also blocks to ensure computations are finished.
    cudaMemcpy(_u, _u_new, nx*ny * sizeof(double), cudaMemcpyDeviceToDevice);
  }

  /*****************************************************************************
   * Step VI: Copy the memory back to the CPU.                                 *
   ****************************************************************************/

  cudaMemcpy(u, _u, nx*ny * sizeof(double), cudaMemcpyDeviceToHost);

  /*****************************************************************************
   * Step I and Step II: Free the memory that you declared and allocated       *
   *                     earlier in the program.                               *
   ****************************************************************************/

  cudaFree(_u);
  cudaFree(_u_new);
  cudaFree(_error);
  free(u);

  return EXIT_SUCCESS;
}

