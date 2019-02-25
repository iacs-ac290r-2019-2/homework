// void gaussian_quadrature_1(float* xx, float* ww);

// void gaussian_quadrature_2(float* xx, float* ww);

// void gaussian_quadrature(size_t n, float* xx, float* ww);

float N_local(size_t a, float xx);

float pN_local(size_t a, float xx);

void NN_local_intmat(float* ddiag, float* ssub);

void pNpN_local_intmat(float* ddiag, float* ssub);

void NN_global_intmat(size_t N, float* diag, float* sub);

void pNpN_global_intmat(size_t N, float* diag, float* sub);

