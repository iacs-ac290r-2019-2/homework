#include <stdio.h>
#include <math.h>

#define nx 200  // points in x direction
#define ny 60  // points in y direction
#define nq 9  // D2Q9 LBM

#define wx 30
#define lx 50

// global constants
int wyh = ny - (ny - wx)/2;
int wyl = (ny - wx)/2;

int b_left = (nx - lx)/2 - 1;
int b_right = (nx - lx)/2 + lx - 1;


int nt = 10000; // number of time steps

double omega = 1.0;  // relaxation parameters, between [0, 2]
double force = 1e-10;  // external forcing in x direction

double cs2 = 1.0/3;  // sound speed squared
double cs = 0.57735;  // sound speed 1/sqrt(3)
double cs4 = 1.0/9;

double w[nq] = {4.0/9, 1.0/9 ,1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36}; // weights for D2Q9


int BC0(double f[ny+2][nx+2][nq]){
    
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
    
    
    // boxes in center
    for (int j = b_left; j < b_right+1; j++) {
        // top
        f[wyh][j][4] = f[wyh+1][j][2];
        f[wyh][j][8] = f[wyh+1][j-1][6];
        f[wyh][j][7] = f[wyh+1][j+1][5];
        // bottom
        f[wyl][j][2] = f[wyl-1][j][4];
        f[wyl][j][5] = f[wyl-1][j-1][7];
        f[wyl][j][6] = f[wyl-1][j+1][8];
    }


    for (int i = 2; i < wyl+1; i++) {
        // left
        f[i][b_left][3] = f[i][b_left+1][1];
        f[i][b_left][7] = f[i+1][b_left+1][5];
        f[i][b_left][6] = f[i-1][b_left+1][8];
        // right
        f[i][b_right][1] = f[i][b_right-1][3];
        f[i][b_right][5] = f[i-1][b_right-1][7];
        f[i][b_right][8] = f[i+1][b_right-1][6];
    }

    for (int i = wyh; i < ny; i++) {
        // left
        f[i][b_left][3] = f[i][b_left+1][1];
        f[i][b_left][7] = f[i+1][b_left+1][5];
        f[i][b_left][6] = f[i-1][b_left+1][8];
        // right
        f[i][b_right][1] = f[i][b_right-1][3];
        f[i][b_right][5] = f[i-1][b_right-1][7];
        f[i][b_right][8] = f[i+1][b_right-1][6];
    }
    
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
            
            if ((j >= wyh) || (j <= wyl)) {
                if ((i <= b_right) && (i >= b_left)){
                    rho[j][i] = 0.0;
                    ux[j][i] = 0.0;
                    uy[j][i] = 0.0;
                    continue;
                }
            }

            if(i==1){
                ux[j][i] = 1.0/36000.0;
                uy[j][i] = 0.0;
                rho[j][i] = 1.0;
                continue;
            }
            
            if (i == nx){
                ux[j][i] = 1.0/36000.0;
                uy[j][i] = 0.0;
                rho[j][i] = 1.0;
                continue;
            }
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
        if ((j >= wyh) || (j <= wyl)) {
                if ((i <= b_right) && (i >= b_left)){
                    feq[j][i][0] = 0.0;
                    feq[j][i][1] = 0.0;
                    feq[j][i][2] = 0.0;
                    feq[j][i][3] = 0.0;
                    feq[j][i][4] = 0.0;
                    feq[j][i][5] = 0.0;
                    feq[j][i][6] = 0.0;
                    feq[j][i][7] = 0.0;
                    feq[j][i][8] = 0.0;
                    continue;
                }
            }
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
                if ((j >= wyh) || (j <= wyl)) {
                    if ((i <= b_right) && (i >= b_left)){
                        continue;
                    }
                }
                f[j][i][q] = (1.0 - omega)*f[j][i][q] + omega*feq[j][i][q];
            }
        }
    }
    return 0;
}


int apply_forcing(double f[ny+2][nx+2][nq], double force){
    
    for (int j = 1; j < ny+1; j++) {
        for (int i = 1; i < nx+1; i++) { 
            if ((j >= wyh) || (j <= wyl)) {
                if ((i <= b_right) && (i >= b_left)){
                    continue;
                }
            }
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
        if ((j >= wyh) || (j <= wyl)) {
                if ((i <= b_right) && (i >= b_left)){
                    rho[j][i] = 0.0;
                    ux[j][i] = 0.0;
                    uy[j][i] = 0.0;
                    continue;
                }
            }
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

int main(int argc, char** argv){

    // prognostic variables
    double f[ny+2][nx+2][nq] = {{{0}}}; // particle population. 
    double feq[ny+2][nx+2][nq] = {{{0}}};  // The equilibrium distribution

    double ux[ny+2][nx+2] = {{0}};  //  macroscopic velocity in x
    double uy[ny+2][nx+2] = {{0}};  //  macroscopic velocity in y
    double rho[ny+2][nx+2] = {{0}};   //macroscopic density (sum of particle in all directions)

    // LBM solver

    initialize(f, feq, rho, ux, uy);
    
    for (int it = 0; it < nt; it++) {
        BC0(f);
        streaming(f);
        macroscopic(f, rho, ux, uy);
        equilibrium(rho, ux, uy, feq);
        collision(f, feq);
        // apply_forcing(f, force);
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