#include <stdio.h>
#include <math.h>
#include "matrix_ops.h"

void test_matrix_ops01() {
    Matrix A = matrix_create(2, 2);
    matrix_set(&A, 0, 0, 1.0);
    matrix_set(&A, 0, 1, 2.0);
    matrix_set(&A, 1, 0, 3.0);
    matrix_set(&A, 1, 1, 4.0);

    Matrix B = matrix_identity(2);

    Matrix C = matrix_mul(&A, &B);
    matrix_print(&C);

    matrix_free(&A);
    matrix_free(&B);
    matrix_free(&C);
}

// 计算矩阵的行列式和逆矩阵
void test_matrix_ops02() {
    // 测试 2x2 单位矩阵
    Matrix mat1 = matrix_identity(2);
    double det1 = matrix_determinant(&mat1);
    printf("Determinant of identity matrix (2x2): %lf\n", det1);
    Matrix inv1 = matrix_inverse(&mat1);
    printf("Inverse of identity matrix (2x2):\n");
    matrix_print(&inv1);

    matrix_free(&inv1);
    matrix_free(&mat1);

    // 测试 2x2 正常矩阵
    Matrix mat2 = matrix_create(2, 2);
    matrix_set(&mat2, 0, 0, 1.0);
    matrix_set(&mat2, 0, 1, 2.0);
    matrix_set(&mat2, 1, 0, 3.0);
    matrix_set(&mat2, 1, 1, 4.0);
    
    double det2 = matrix_determinant(&mat2);
    printf("Determinant of matrix (2x2): %lf\n", det2);

    Matrix inv2 = matrix_inverse(&mat2);
    printf("Inverse of matrix (2x2):\n");
    matrix_print(&inv2);

    matrix_free(&inv2);
    matrix_free(&mat2);

    // 测试奇异矩阵（行列式为 0）
    Matrix mat3 = matrix_create(2, 2);
    matrix_set(&mat3, 0, 0, 1.0);
    matrix_set(&mat3, 0, 1, 2.0);
    matrix_set(&mat3, 1, 0, 2.0);
    matrix_set(&mat3, 1, 1, 4.0);
    
    double det3 = matrix_determinant(&mat3);
    printf("Determinant of singular matrix (2x2): %lf\n", det3);
    
    // 矩阵不可逆，应该返回NaN或0
    if (isnan(det3) || det3 == 0.0) {
        printf("Matrix is singular, cannot compute inverse.\n");
    } else {
        Matrix inv3 = matrix_inverse(&mat3);
        printf("Inverse of singular matrix (2x2):\n");
        matrix_print(&inv3);
        matrix_free(&inv3);
    }

    matrix_free(&mat3);
}

int main() {
    // test_matrix_ops01();
    test_matrix_ops02();
    return 0;
}
