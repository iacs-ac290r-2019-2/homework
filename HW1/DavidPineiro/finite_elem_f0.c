/**

Finite Elements Method Code

**/
#include <stdio.h>

int N = 10;
int h_bound = 1;
int g_bound = 1;

// x will always be between 0 and 1
float N_A (float x, int A, int n){
	float X_A = A/n;
	float X_A_n1 = (A - 1)/n;
	float X_A_p1 = (A + 1)/n;

	if (x >= X_A_n1 && x <= X_A){
		return (x - X_A_n1)/(X_A - X_A_n1);
	} if (x >= X_A && x <= X_A_p1){
		return (X_A_p1 - x)/(X_A_p1 - X_A);
	} else{
		return 0;
	}
}

float partial_N_A(float x, int A, int n){
	float X_A = A/n;
	float X_A_n1 = (A - 1)/n;
	float X_A_p1 = (A + 1)/n;

	if (x >= X_A_n1 && x <= X_A){
		return 1;
	} if (x >= X_A && x <= X_A_p1){
		return -1;
	} else{
		return 0;
	}
}

float h_a(int A, int n){
	float X_A = (float)A / n;
	float X_A_n1 = (A - 1.0) / n;
	return X_A - X_A_n1;
}

void K_AB_matrix(int n, float *diag, float *t_diag){
	for(int j=0; j<n-1; ++j){
		float h = h_a(j, n);
		diag[j] = 2/h;
		t_diag[j] = -1/h;
	}
	diag[n-1] = 2/h_a(n-1, n);
}

float f(float x){
	return 0;
}


float F_a(float x){
	return 0;
}

float G_a(float x, int A){
	return g_bound*partial_N_A(x, A, N)*(1.0/h_a(A, N));
}


float H_a(float x, int A){
	return N_A(x, A, N)*h_bound; 
}


void F(float x){

}

int main(){
	// float a[5] = {1, 2, 3, 4, 5};
	// float diag[5];
	// float t_diag[4];
	// // float result = h_a(3, a);
	// K_AB_matrix(5, diag, t_diag);

	// for(int j=0; j<5; j++){
	// 	printf("%f ", diag[j]);
	// }
	// printf("neg one = %f \n", simpson(*f_1, 100, 0, 10));
	// printf("mix = %f \n", f_mixer(*f_1, *f_1)(3.0));
	// printf("hello world = %f\n", result);
}