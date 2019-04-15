import numpy as np
from numba import jit
import xarray as xr

# == parse command line arguments
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--omega", type=float, default=1.0,
                    help="relaxation parameter, between [0, 2]")
parser.add_argument("--force", type=float, default=1e-5,
                    help="external forcing in x direction")
parser.add_argument("--wx", type=int, default=30,
                    help="width of the narrowing in region, between [0, 60]")
parser.add_argument("--nt", type=int, default=1000,
                    help="total number of time steps")
args = parser.parse_args()

print("Arguments:", args)

# == constants ==
nt = args.nt  # time steps
nx = 200  # points in x direction
ny = 60  # points in y direction
nq = 9  # D2Q9 LBM
omega = args.omega # relaxation parameters, between [0, 2]
force = args.force  # external forcing in x direction

cs2 = 1/3  # sound speed squared
cs = np.sqrt(cs2)  # sound speed
cs4 = cs2**2

# BC
wx = args.wx
lx = 50
wyh = int(ny - (ny - wx)/2)
wyl = int((ny - wx)/2)
b_left = int((nx - lx)/2)  # 74
b_right = int((nx - lx)/2 + lx)  #124

#  weights for D2Q9
w = np.array([4/9,
     1/9 ,1/9, 1/9, 1/9,
     1/36, 1/36, 1/36, 1/36])

# == prognostic variables ==
# Is (ny, nx, nq) a more efficient memory layout than (nq, ny, nx) ?
# +2 is for buffering regions
f = np.zeros((ny+2, nx+2, nq))  # particle population. 
feq = np.zeros((ny+2, nx+2, nq))  # The equilibrium distribution

ux = np.zeros((ny+2, nx+2))  #  macroscopic velocity in x
uy = np.zeros((ny+2, nx+2))  #  macroscopic velocity in y
rho = np.zeros((ny+2, nx+2))  # macroscopic density (sum of particle in all directions)


@jit
def apply_forcing(f, force):
    for j in range(1, ny+1):
        for i in range(1, nx+1):
            if j >= wyh or j <= wyl:
                if i <= b_right and i >= b_left:
                    continue
            f[j, i, 1] += force
            f[j, i, 5] += force
            f[j, i, 8] += force

            f[j, i, 3] -= force
            f[j, i, 6] -= force
            f[j, i, 7] -= force
         

@jit
def equilibrium(rho, ux, uy, feq):   
    for j in range(1, ny+1):
        for i in range(1, nx+1):
            if j >= wyh or j <= wyl:
                if i <= b_right and i >= b_left:
                    feq[j,i,1] = 0.0
                    feq[j,i,2] = 0.0
                    feq[j,i,3] = 0.0
                    feq[j,i,4] = 0.0
                    continue
            # == expressions that can be reused ==  
            # Comments are the corresponding variable names in Sauro Succi's Fortran code
            ux_d_cs2 = ux[j, i] / cs2  # ui
            uy_d_cs2 = uy[j, i] / cs2  # vi
            ux2_d_2cs4 = ux[j, i]**2 / (2*cs4)  # u2
            uy2_d_2cs4 = uy[j, i]**2 / (2*cs4)  # v2
            
            v2_d_2cs2 = (ux[j, i]**2 + uy[j, i]**2) / (2*cs2)  # sumsq (v = velocity magnitude)
            v2_d_sp = v2_d_2cs2 * (1.0 - cs2) / cs2  # sumsq2  (sp = special)
            uxuy_d_cs4 = ux_d_cs2 * uy_d_cs2  # uv
        
            # == computing equilibrium ==
            feq[j, i, 0] = rho[j, i]*w[0]*(1 - v2_d_2cs2)

            feq[j, i, 1] = rho[j, i]*w[1]*(1 - v2_d_2cs2 + ux_d_cs2 + ux2_d_2cs4)
            feq[j, i, 2] = rho[j, i]*w[2]*(1 - v2_d_2cs2 + uy_d_cs2 + uy2_d_2cs4)
            feq[j, i, 3] = rho[j, i]*w[3]*(1 - v2_d_2cs2 - ux_d_cs2 + ux2_d_2cs4)
            feq[j, i, 4] = rho[j, i]*w[4]*(1 - v2_d_2cs2 - uy_d_cs2 + uy2_d_2cs4)

            feq[j, i, 5] = rho[j, i]*w[5]*(1 + v2_d_sp + ux_d_cs2 + uy_d_cs2 + uxuy_d_cs4)
            feq[j, i, 6] = rho[j, i]*w[6]*(1 + v2_d_sp - ux_d_cs2 + uy_d_cs2 - uxuy_d_cs4)
            feq[j, i, 7] = rho[j, i]*w[7]*(1 + v2_d_sp - ux_d_cs2 - uy_d_cs2 + uxuy_d_cs4)
            feq[j, i, 8] = rho[j, i]*w[8]*(1 + v2_d_sp + ux_d_cs2 - uy_d_cs2 - uxuy_d_cs4)
            
            
@jit
def collision(f, feq):
    for j in range(1, ny+1):
        for i in range(1, nx+1):
            if j >= wyh or j <= wyl:
                if i <= b_right and i >= b_left:
                    continue
            for q in range(nq):
                f[j, i, q] = (1.0 - omega)*f[j, i, q] + omega*feq[j, i, q]


@jit
def streaming(f):
    f_temp = f.copy()  # copy can be avoided with write-back strategy    
    for j in range(1, ny+1):
        for i in range(1, nx+1):
            # first-order neighbours
            f[j, i, 1] = f_temp[j, i-1, 1]
            f[j, i, 2] = f_temp[j-1, i, 2]
            f[j, i, 3] = f_temp[j, i+1, 3]
            f[j, i, 4] = f_temp[j+1, i, 4]
            
            # second-order neighbours
            f[j, i, 5] = f_temp[j-1, i-1, 5]
            f[j, i, 6] = f_temp[j-1, i+1, 6]
            f[j, i, 7] = f_temp[j+1, i+1, 7]
            f[j, i, 8] = f_temp[j+1, i-1, 8]
            
            
@jit
def macroscopic(f, rho, ux, uy):
    
    for j in range(1, ny+1):
        for i in range(1, nx+1):
            if j >= wyh or j <= wyl:
                if i <= b_right and i >= b_left:
                    rho[j, i] = np.nan  #1.0
                    ux[j, i] = np.nan  #0.0
                    uy[j, i] = np.nan  #0.0
                    continue
            rho[j, i] = f[j, i, :].sum()
            
                # unroll loop explictly to avoid summing over zeros
            ux[j, i] = (f[j, i, 1] + f[j, i, 5] + f[j, i, 8] 
                           - f[j, i, 3] - f[j, i, 6] - f[j, i, 7]
                          )/rho[j, i]

            uy[j, i] = (f[j, i, 2] + f[j, i, 5] + f[j, i, 6] 
                           - f[j, i, 4] - f[j, i, 7] - f[j, i, 8]
                          )/rho[j, i]


@jit
def BC0(f):
    # north wall, bounce-back
    for i in range(1, nx+1):
        f[ny+1, i, 4] = f[ny, i, 2]
        f[ny+1, i, 8] = f[ny, i+1, 6]
        f[ny+1, i, 7] = f[ny, i-1, 5]
        
    # south wall, bounce-back
    for i in range(1, nx+1):
        f[0, i, 2] = f[1, i, 4]
        f[0, i, 6] = f[1, i-1, 8]
        f[0, i, 5] = f[1, i+1, 7]
    
    # west interface, periodic
    for j in range(1, ny+1):
        f[j, 0, 1] = f[j, nx, 1]
        f[j, 0, 5] = f[j, nx, 5]
        f[j, 0, 8] = f[j, nx, 8]
    
    # east interface, periodic
    for j in range(1, ny+1):
        f[j, nx+1, 3] = f[j, 1, 3]
        f[j, nx+1, 6] = f[j, 1, 6]
        f[j, nx+1, 7] = f[j, 1, 7] 
    
    
    for j in range(b_left, b_right+1):
        # top
        f[wyh, j, 4] = f[wyh+1, j, 2]
        f[wyh, j, 8] = f[wyh+1, j-1, 6]
        f[wyh, j, 7] = f[wyh+1, j+1, 5]
        # bottom 
        f[wyl, j, 2] = f[wyl-1, j, 4]
        f[wyl, j, 5] = f[wyl-1, j-1, 7]
        f[wyl, j, 6] = f[wyl-1, j+1, 8]
        
    for i in range(1, wyl+1):
        # left
        f[i, b_left, 3] = f[i, b_left+1, 1]
        f[i, b_left, 7] = f[i+1, b_left+1, 5]
        f[i, b_left, 6] = f[i-1, b_left+1, 8]
        # right
        f[i, b_right, 1] = f[i, b_right-1, 3]
        f[i, b_right, 5] = f[i-1, b_right-1, 7]
        f[i, b_right, 8] = f[i+1, b_right-1, 6]
    
    for i in range(wyh, ny+1):
        # left
        f[i, b_left, 3] = f[i, b_left+1, 1]
        f[i, b_left, 7] = f[i+1, b_left+1, 5]
        f[i, b_left, 6] = f[i-1, b_left+1, 8]
        #right
        f[i, b_right, 1] = f[i, b_right-1, 3]
        f[i, b_right, 5] = f[i-1, b_right-1, 7]
        f[i, b_right, 8] = f[i+1, b_right-1, 6]
        
    # square corners
    # top square
    f[ny, b_left, 7] = f[ny+1, b_left+1, 5] # north-west
    f[ny, b_right, 8] = f[ny+1, b_right-1, 6] # north-east
    f[wyh, b_right, 8] = f[wyh+1, b_right-1, 6] # south-west
    f[wyh, b_left, 7] = f[wyh+1, b_left+1, 5] # south-east
    # bottom square
    f[wyl, b_left, 6] = f[wyl-1, b_left+1, 8]
    f[wyl, b_right, 5] = f[wyh-1, b_right-1, 7]   
    f[0, b_left, 6] = f[1, b_left+1, 8]
    f[0, b_right, 5] = f[1, b_right-1, 7]

    # 4 corners, bounce-back
    f[ny+1, 0, 8] = f[ny, 1, 6]  # north-west
    f[0, 0, 5] = f[1, 1, 7]  # south-west
    f[ny+1, nx+1, 7] = f[ny, nx, 5]  # north-east
    f[0, nx+1, 6] = f[1, nx, 8]  # south-east
    
    
@jit
def initialize(f, feq, rho, ux, uy):
    '''
    Output
    ------
    f: 3D numpy array, population
    feq: 3D numpy array, equilibrium population
    rho: 2D numpy array, density
    ux, uy: 2D numpy array, velocity
    '''
    
    # uniform density
    rho[:] = 1.0
    for j in range(1, ny+1):
        for i in range(1, nx+1):
            if j >= wyh or j <= wyl:
                if i <= b_right and i >= b_left:
                    rho[j,i] = 0.0
                    continue
    ux[:] = 0.0
    uy[:] = 0.0
    
    # use equilibrium as initial population
    equilibrium(rho, ux, uy, feq)
    f[:] = feq


@jit
def lbm_solver(f, rho, ux, uy, nt):
    
    initialize(f, feq, rho, ux, uy)
    
    for it in range(0, nt):
        BC0(f)
        streaming(f)
        macroscopic(f, rho, ux, uy)
        equilibrium(rho, ux, uy, feq)
        collision(f, feq)
        apply_forcing(f, force)


if __name__ == "__main__":
    filename = 'lbm_w{0}_omega{1}_force{2}_nt{3}.nc'.format(wx, omega, force, nt)  # output file

    print('Important parameters:')
    print('Time steps:', nt)
    print('Narrowing width w:', wx)
    print('Forcing:', force)
    print('Lelaxation omega:', omega)
    print('viscosity: ', cs2*(1/omega-0.5))
    print('================ \n')

    print('running solver... \n')
    lbm_solver(f, rho, ux, uy, nt)

    # post processing
    coords={'x': np.arange(1, nx+1, dtype=float), 'y': np.arange(1, ny+1, dtype=float)}
    dr_ux = xr.DataArray(ux[1:-1, 1:-1], name='ux', dims=('y', 'x'), coords=coords)
    dr_uy = xr.DataArray(uy[1:-1, 1:-1], name='uy', dims=('y', 'x'), coords=coords)
    dr_rho = xr.DataArray(rho[1:-1, 1:-1], name='rho', dims=('y', 'x'), coords=coords)
    dr_p = dr_rho.rename('p') / 6  # pressure field 

    ds = xr.merge([dr_ux, dr_uy, dr_rho, dr_p])
    print('writing to file:', filename)
    ds.to_netcdf(filename)
