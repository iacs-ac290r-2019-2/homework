#include <stdio.h>
#include <math.h>

#define wx 30
#define lx 50
#define nx 200  // points in x direction
#define ny 60  // points in y direction
#define nq 9  // D2Q9 LBM

// global constants
int nt = 2000; // number of time steps

double omega = 1.0;  // relaxation parameters, between [0, 2]
double force = 1e-4;  // external forcing in x direction

double cs2 = 1.0/3;  // sound speed squared
double cs = 0.57735;  // sound speed 1/sqrt(3)
double cs4 = 1.0/9;

int wyh = ny - (ny - wx)/2;
int wyl = (ny - wx)/2;

int b_left = (nx - lx)/2; // 74
int b_right = (nx - lx)/2 + lx; //124


double w[nq] = {4.0/9, 1.0/9 ,1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36}; // weights for D2Q9

int BC1(double f[ny+2][nx+2][nq]){
    
    // north wall, bounce-back
    for (int i = 1; i < nx+1; i++) {

        if (i >= b_left && i <= b_right){ //  
            f[wyh+1][i][4] = f[wyh][i][2];
            f[wyh+1][i][8] = f[wyh][i+1][6];
            f[wyh+1][i][7] = f[wyh][i-1][5];
        } else {
            f[ny+1][i][4] = f[ny][i][2];
            f[ny+1][i][8] = f[ny][i+1][6];
            f[ny+1][i][7] = f[ny][i-1][5];
        }
    }
        
    // south wall, bounce-back
    for (int i = 1; i < nx+1; i++) {
        if (i >= b_left && i <= b_right){ //  
            f[wyl-1][i][2] = f[wyl][i][4];
            f[wyl-1][i][6] = f[wyl][i-1][8];
            f[wyl-1][i][5] = f[wyl][i+1][7];
        } else {
            f[0][i][2] = f[1][i][4];
            f[0][i][6] = f[1][i-1][8];
            f[0][i][5] = f[1][i+1][7];
        }
    }
    
    // west interface, periodic
    for (int j = 1; j < ny+1; j++) {
        // west wall, bounce-back
        if (j >= wyh || j <= wyl){
            f[j][b_right-1][6] = f[j][b_right][8];
            f[j][b_right-1][3] = f[j][b_right][1];
            f[j][b_right-1][7] = f[j][b_right][5];
        }
        // periodic
        f[j][0][1] = f[j][nx][1];
        f[j][0][5] = f[j][nx][5];
        f[j][0][8] = f[j][nx][8];
    }
    
    // east interface, periodic
    for (int j = 1; j < ny+1; j++) {
        // east wall, bounce-back
        if (j >= wyh || j <= wyl){
            f[j][b_left+1][8] = f[j][b_left][6];
            f[j][b_left+1][1] = f[j][b_left][3];
            f[j][b_left+1][5] = f[j][b_left][7];
        }
        f[j][nx+1][3] = f[j][1][3];
        f[j][nx+1][6] = f[j][1][6];
        f[j][nx+1][7] = f[j][1][7];
    }
    
    // 4 corners, bounce-back
    f[ny+1][0][8] = f[ny][1][6];  // north-west
    f[0][0][5] = f[1][1][7];  // south-west
    f[ny+1][nx+1][7] = f[ny][nx][5];  // north-east
    f[0][nx+1][6] = f[1][nx][8];  // south-east

    // middle square corners
    // top square
    f[ny+1][b_left-1][7] = f[ny][b_left][5];  //check this
    f[wyh+1][b_left-1][7] = f[wyh][b_left][5];

    f[ny+1][b_right][8] = f[ny][b_right-1][6];  
    f[wyh+1][b_right][8] = f[wyh][b_right-1][6];
    // bottom square
    f[1][b_left-1][6] = f[0][b_left][8];  
    f[wyl][b_left-1][6] = f[wyl][b_left][8];
    f[wyl][b_right-1][5] = f[wyl][b_right-1][7];  
    f[1][b_right-1][5] = f[0][b_right-1][7];

    return 0;

}

int fill_buffer(double f[ny+2][nx+2][nq]){
    
    // north wall, bounce-back
    for (int i = 1; i < nx+1; i++) {
        f[ny+1][i][4] = f[ny][i][2];
        f[ny+1][i][8] = f[ny][i+1][6];
        f[ny+1][i][7] = f[ny][i-1][5];
    }
        
    // south wall, bounce-back
    for (int i = 1; i < nx+1; i++) {
        f[0][i][2] = f[1][i][4];
        f[0][i][6] = f[1][i-1][8];
        f[0][i][5] = f[1][i+1][7];
    }
    
    // west interface, periodic
    for (int j = 1; j < ny+1; j++) {
        f[j][0][1] = f[j][nx][1];
        f[j][0][5] = f[j][nx][5];
        f[j][0][8] = f[j][nx][8];
    }
    
    // east interface, periodic
    for (int j = 1; j < ny+1; j++) {
        f[j][nx+1][3] = f[j][1][3];
        f[j][nx+1][6] = f[j][1][6];
        f[j][nx+1][7] = f[j][1][7];
    }
    
    // 4 corners, bounce-back
    f[ny+1][0][8] = f[ny][1][6];  // north-west
    f[0][0][5] = f[1][1][7];  // south-west
    f[ny+1][nx+1][7] = f[ny][nx][5];  // north-east
    f[0][nx+1][6] = f[1][nx][8];  // south-east

    return 0;

}



int streaming(double f[ny+2][nx+2][nq]){
    
    double f_temp[ny+2][nx+2][nq];  // copy can be avoided with write-back strategy

    // notice that the loop includes the buffering region (unlike other loops)!
    // it is important to also copy the buffer,
    // otherwise streaming at the boundary will be wrong!
    for (int j = 0; j < ny+2; j++) {
    for (int i = 0; i < nx+2; i++) {
        for (int q = 0; q < nq; q++) {
            f_temp[j][i][q] = f[j][i][q];
        }
    }
    }
   

    for (int j = 1; j < ny+1; j++) {
    for (int i = 1; i < nx+1; i++) {
            
        // first-order neighbours
        f[j][i][1] = f_temp[j][i-1][1];
        f[j][i][2] = f_temp[j-1][i][2];
        f[j][i][3] = f_temp[j][i+1][3];
        f[j][i][4] = f_temp[j+1][i][4];
        
        // second-order neighbours
        f[j][i][5] = f_temp[j-1][i-1][5];
        f[j][i][6] = f_temp[j-1][i+1][6];
        f[j][i][7] = f_temp[j+1][i+1][7];
        f[j][i][8] = f_temp[j+1][i-1][8];
    }
    }
    return 0;
}


int macroscopic(double f[ny+2][nx+2][nq], double rho[ny+2][nx+2], 
                double ux[ny+2][nx+2], double uy[ny+2][nx+2]){

    for (int j = 1; j < ny+1; j++) {
    for (int i = 1; i < nx+1; i++) {
            
            // numpy version: rho[j, i] = f[j, i, :].sum()
            rho[j][i] = 0.0;
            for (int q = 0; q < nq; q++) {
                rho[j][i] += f[j][i][q];
            }
            
            // unroll loop explictly to avoid summing over zeros
            ux[j][i] = (f[j][i][1] + f[j][i][5] + f[j][i][8] 
                       - f[j][i][3] - f[j][i][6] - f[j][i][7]
                      )/rho[j][i];
            
            uy[j][i] = (f[j][i][2] + f[j][i][5] + f[j][i][6] 
                       - f[j][i][4] - f[j][i][7] - f[j][i][8]
                      )/rho[j][i];

    }
    }
    return 0;
}


int equilibrium(double rho[ny+2][nx+2], double ux[ny+2][nx+2], 
                double uy[ny+2][nx+2], double feq[ny+2][nx+2][nq]){
    
    double ux_d_cs2, uy_d_cs2, ux2_d_2cs4, uy2_d_2cs4;
    double v2_d_2cs2, v2_d_sp, uxuy_d_cs4;

    for (int j = 1; j < ny+1; j++) {
    for (int i = 1; i < nx+1; i++) {
            
        // == expressions that can be reused ==  
        // Comments are the corresponding variable names in Sauro Succi's Fortran code
        ux_d_cs2 = ux[j][i] / cs2;  // ui
        uy_d_cs2 = uy[j][i] / cs2;  // vi
        ux2_d_2cs4 = ux[j][i]*ux[j][i] / (2*cs4);  // u2
        uy2_d_2cs4 = uy[j][i]*uy[j][i] / (2*cs4);  // v2
        
        v2_d_2cs2 = (ux[j][i]*ux[j][i] + uy[j][i]*uy[j][i]) / (2*cs2);  // sumsq (v = velocity magnitude)
        v2_d_sp = v2_d_2cs2 * (1.0 - cs2) / cs2;  // sumsq2  (sp = special)
        uxuy_d_cs4 = ux_d_cs2 * uy_d_cs2;  // uv
    
        // == computing equilibrium ==
        feq[j][i][0] = rho[j][i]*w[0]*(1 - v2_d_2cs2);

        feq[j][i][1] = rho[j][i]*w[1]*(1 - v2_d_2cs2 + ux_d_cs2 + ux2_d_2cs4);
        feq[j][i][2] = rho[j][i]*w[2]*(1 - v2_d_2cs2 + uy_d_cs2 + uy2_d_2cs4);
        feq[j][i][3] = rho[j][i]*w[3]*(1 - v2_d_2cs2 - ux_d_cs2 + ux2_d_2cs4);
        feq[j][i][4] = rho[j][i]*w[4]*(1 - v2_d_2cs2 - uy_d_cs2 + uy2_d_2cs4);

        feq[j][i][5] = rho[j][i]*w[5]*(1 + v2_d_sp + ux_d_cs2 + uy_d_cs2 + uxuy_d_cs4);
        feq[j][i][6] = rho[j][i]*w[6]*(1 + v2_d_sp - ux_d_cs2 + uy_d_cs2 - uxuy_d_cs4);
        feq[j][i][7] = rho[j][i]*w[7]*(1 + v2_d_sp - ux_d_cs2 - uy_d_cs2 + uxuy_d_cs4);
        feq[j][i][8] = rho[j][i]*w[8]*(1 + v2_d_sp + ux_d_cs2 - uy_d_cs2 - uxuy_d_cs4);

    }
    }
    return 0;        
}

int collision(double f[ny+2][nx+2][nq], double feq[ny+2][nx+2][nq]){

    for (int j = 1; j < ny+1; j++) {
    for (int i = 1; i < nx+1; i++) {
        for (int q = 0; q < nq; q++) {
                f[j][i][q] = (1.0 - omega)*f[j][i][q] + omega*feq[j][i][q];
        }
    }
    }
    return 0;
}


int apply_forcing(double f[ny+2][nx+2][nq], double force){
    
    for (int j = 1; j < ny+1; j++) {
    for (int i = 1; i < nx+1; i++) { 
            f[j][i][1] += force;
            f[j][i][5] += force;
            f[j][i][8] += force;

            f[j][i][3] -= force;
            f[j][i][6] -= force;
            f[j][i][7] -= force;
    }
    }
    return 0;
}

int initialize(double f[ny+2][nx+2][nq], double feq[ny+2][nx+2][nq], 
               double rho[ny+2][nx+2], double ux[ny+2][nx+2], double uy[ny+2][nx+2]){
    
    // uniform density
    for (int j = 1; j < ny+1; j++) {
    for (int i = 1; i < nx+1; i++) { 
        rho[j][i] = 1.0;
        ux[j][i] = 0.0;
        uy[j][i] = 0.0;
    }
    }
    
    // use equilibrium as initial population
    equilibrium(rho, ux, uy, feq);
    for (int j = 1; j < ny+1; j++) {
    for (int i = 1; i < nx+1; i++) {
        for (int q = 0; q < nq; q++) {
            f[j][i][q] = feq[j][i][q];
        }
    }
    }

    return 0;
}

int main(){

    // prognostic variables
    double f[ny+2][nx+2][nq] = {{{0}}}; // particle population. 
    double feq[ny+2][nx+2][nq] = {{{0}}};  // The equilibrium distribution

    double ux[ny+2][nx+2] = {{0}};  //  macroscopic velocity in x
    double uy[ny+2][nx+2] = {{0}};  //  macroscopic velocity in y
    double rho[ny+2][nx+2] = {{0}};   //macroscopic density (sum of particle in all directions)

    // LBM solver

    initialize(f, feq, rho, ux, uy);
    
    for (int it = 0; it < nt; it++) {
        BC1(f);
        streaming(f);
        macroscopic(f, rho, ux, uy);
        equilibrium(rho, ux, uy, feq);
        collision(f, feq);
        apply_forcing(f, force);
    }

    // print result to check
    for (int j = 0; j < ny+2; j++) {
        for (int i = 0; i < nx+2; i++) {
            printf("%f ", ux[j][i]);
        }
        printf("\n");
    }

    return 0;
}