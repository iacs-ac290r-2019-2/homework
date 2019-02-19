/**

Finite Elements Method Code

**/
#include <stdio.h>

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

float f_1(float x){
	return x*x;
}

/* Integration using Simpson's rule */ 
float simpson (float (*f)(float), int no, float min, float max){  
   int n;				 
   float interval, sum=0., x;
   interval = ((max -min) /(no-1));
   
   for (n=2; n<no; n+=2)                /* loop for odd points  */
   {
       x = interval * (n-1);
       sum += 4 * f(x);
   }
   for (n=3; n<no; n+=2)                /* loop for even points  */
   {
      x = interval * (n-1);
      sum += 2 * f(x);
   }   
   sum +=  f(min) + f(max);	 	/* add first and last value */          
   sum *= interval/3.;        		/* then multilpy by interval*/
   
   return (sum);
}  

// float Na_f(float x){
// 	return 
// }

// float N_A_F(float x){
// 	return f_1(x)*
// }

// float f_mixer(float (*f_1)(float), float (*f_2)(float)){
// 	return *float mixer(float x){return f_1(x)*f_2(x);};
// }


// void F_A(float min, float max, int n, float *A, float *M){

// 	for(int i; i<n, i++){

// 		M[i] = simpson(, 1000, min, max)
// 	}
// 	float fa = simpson()

// }


int main(){
	// float a[5] = {1, 2, 3, 4, 5};
	float diag[5];
	float t_diag[4];
	// float result = h_a(3, a);
	K_AB_matrix(5, diag, t_diag);

	for(int j=0; j<5; j++){
		printf("%f ", diag[j]);
	}
	// printf("neg one = %f \n", simpson(*f_1, 100, 0, 10));
	// printf("mix = %f \n", f_mixer(*f_1, *f_1)(3.0));
	// printf("hello world = %f\n", result);
}