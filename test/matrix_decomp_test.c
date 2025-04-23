#include <stdio.h>
#include "matrix_ops.h"
#include "matrix_decomp.h"


void test_lu_decomposition() {
    printf("=== LU Decomposition Test ===\n");
    Matrix A = matrix_create(3, 3);
    matrix_set(&A, 0, 0, 2); matrix_set(&A, 0, 1, -1); matrix_set(&A, 0, 2, -2);
    matrix_set(&A, 1, 0, -4); matrix_set(&A, 1, 1, 6); matrix_set(&A, 1, 2, 3);
    matrix_set(&A, 2, 0, -4); matrix_set(&A, 2, 1, -2); matrix_set(&A, 2, 2, 8);

    matrix_print(&A);

    Matrix L, U;
    matrix_LU_decomposition(&A, &L, &U);

    printf("\nL matrix:\n");
    matrix_print(&L);

    printf("\nU matrix:\n");
    matrix_print(&U);

    matrix_free(&A);
    matrix_free(&L);
    matrix_free(&U);
}

void test_qr_decomposition() {
    printf("\n=== QR Decomposition Test ===\n");
    Matrix A = matrix_create(3, 2);
    matrix_set(&A, 0, 0, 1); matrix_set(&A, 0, 1, 1);
    matrix_set(&A, 1, 0, 1); matrix_set(&A, 1, 1, -1);
    matrix_set(&A, 2, 0, 1); matrix_set(&A, 2, 1, 1);

    matrix_print(&A);

    Matrix Q, R;
    matrix_QR_decomposition(&A, &Q, &R);

    printf("\nQ matrix:\n");
    matrix_print(&Q);

    printf("\nR matrix:\n");
    matrix_print(&R);

    matrix_free(&A);
    matrix_free(&Q);
    matrix_free(&R);
}

int main() {
    test_lu_decomposition();
    test_qr_decomposition();
    return 0;
}
