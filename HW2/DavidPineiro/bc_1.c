#include <stdio.h>

#define w 30
#define l 50
#define nx 200  // L; points in x direction
#define ny 60  // H; points in y direction


int wyh = ny - (ny - w)/2;
int wyl = (ny - w)/2;

int b_left = (nx - l)/2 - 1; // 74
int b_right = (nx - l)/2 + l - 1; //124

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
    		f[wyl][i][2] = f[wyl + 1][i][4];
        	f[wyl][i][6] = f[wyl + 1][i-1][8];
        	f[wyl][i][5] = f[wyl + 1][i+1][7];
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
    		f[j][b_right][6] = f[j][b_right + 1][8];
        	f[j][b_right][3] = f[j][b_right + 1][1];
        	f[j][b_right][7] = f[j][b_right + 1][5];
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
    		f[j][b_left][6] = f[j][b_left + 1][8];
        	f[j][b_left][3] = f[j][b_left + 1][1];
        	f[j][b_left][7] = f[j][b_left + 1][5];
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
    f[ny+1][bl+1][5] = f[ny][bl][7];  //check this
    f[wyh+1][bl+1][5] = f[wyh][bl][7];
    f[ny+1][br-1][6] = f[ny][br][8];  
    f[wyh+1][br-1][6] = f[wyh][bl][8];
    // bottom square
    f[0][bl+1][8] = f[ny][bl][6];  
    f[wyl-1][bl+1][8] = f[wyl][bl][6];
    f[wyl-1][br-1][7] = f[wyl][br][5];  
    f[0][br-1][7] = f[wyh][bl][5];

    return 0;

}