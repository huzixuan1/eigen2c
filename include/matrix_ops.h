#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include <stddef.h>  // for size_t

typedef struct {
    size_t rows;
    size_t cols;
    double* data;
} Matrix;

// 创建矩阵
Matrix matrix_create(size_t rows, size_t cols);

// 释放矩阵
void matrix_free(Matrix* mat);

// 设置矩阵值
void matrix_set(Matrix* mat, size_t row, size_t col, double value);

// 获取矩阵值
double matrix_get(const Matrix* mat, size_t row, size_t col);

// 矩阵赋值
void matrix_copy(Matrix* dest, const Matrix* src);

// 矩阵加法
Matrix matrix_add(const Matrix* a, const Matrix* b);

// 矩阵减法
Matrix matrix_sub(const Matrix* a, const Matrix* b);

// 矩阵乘法
Matrix matrix_mul(const Matrix* a, const Matrix* b);

// 矩阵转置
Matrix matrix_transpose(const Matrix* mat);

// 求矩阵的逆
Matrix matrix_inverse(const Matrix* mat);

// 求矩阵的行列式
double matrix_determinant(const Matrix* mat);

// 单位矩阵
Matrix matrix_identity(size_t size);

// 打印矩阵（调试用）
void matrix_print(const Matrix* mat);


#endif  // MATRIX_OPS_H
