// build: gcc blas_utils.c -O2 -lopenblas -o gemm

#include <stdio.h>
#include <cblas.h>

int main(void) {
    const int M = 3, K = 2, N = 3;

    double A[3*2] = {
        1, 4, 7,  // col 0
        2, 5, 8   // col 1
    };

    double B[2*3] = {
        1, 2,  // col 0
        3, 4,  // col 1
        5, 6   // col 2
    };

    double C[3*3] = {0};

    const double alpha = 1.0, beta = 0.0;

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                M, N, K,
                alpha,
                A, M,
                B, K,
                beta,
                C, M);

    printf("C =\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            printf("%5.1f ", C[j*M + i]);
        }
        printf("\n");
    }
    return 0;
}
