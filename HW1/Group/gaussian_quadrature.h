void gaussian_quadrature_1(float* xx, float* ww);

void gaussian_quadrature_2(float* xx, float* ww);

void gaussian_quadrature(size_t n, float* xx, float* ww);

float local_N(float xx, size_t a);

float local_pN(float xx, size_t a);

float int_local_NN(size_t a, size_t b);

float int_local_pNpN(size_t a, size_t b);

