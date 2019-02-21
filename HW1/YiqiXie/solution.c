float solution(float x, float* u, float g, int n) {
    float s = 0.;
    for (int i=0; i<n; i++) {
        s += u[i] * N_A(x, i+1, n);
    }
    s += g * N_A(x, n+1, n);
}