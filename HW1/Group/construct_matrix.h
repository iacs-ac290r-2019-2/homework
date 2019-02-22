// float f(float x);

// void compute_K(int N, float *K_diag, float *K_subdiag);

// void compute_F(int N, float h_bound, float g_bound, float *F);

void compute_KF(size_t N, float h_bound, float g_bound, float* f, float* K_diag, float* K_sub, float* F);