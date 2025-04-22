#include "matrix_ops.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define EPSILON 1e-10  // 判断是否为0的小数误差

Matrix matrix_create(size_t rows, size_t cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (double*)calloc(rows * cols, sizeof(double));
    return mat;
}

void matrix_free(Matrix* mat) {
    free(mat->data);
    mat->data = NULL;
    mat->rows = 0;
    mat->cols = 0;
}

void matrix_set(Matrix* mat, size_t row, size_t col, double value) {
    mat->data[row * mat->cols + col] = value;
}

double matrix_get(const Matrix* mat, size_t row, size_t col) {
    return mat->data[row * mat->cols + col];
}

void matrix_copy(Matrix* dest, const Matrix* src) {
    if (dest->rows != src->rows || dest->cols != src->cols) return;
    memcpy(dest->data, src->data, sizeof(double) * src->rows * src->cols);
}

Matrix matrix_add(const Matrix* a, const Matrix* b) {
    Matrix result = matrix_create(a->rows, a->cols);
    for (size_t i = 0; i < a->rows * a->cols; ++i) {
        result.data[i] = a->data[i] + b->data[i];
    }
    return result;
}

Matrix matrix_sub(const Matrix* a, const Matrix* b) {
    Matrix result = matrix_create(a->rows, a->cols);
    for (size_t i = 0; i < a->rows * a->cols; ++i) {
        result.data[i] = a->data[i] - b->data[i];
    }
    return result;
}

Matrix matrix_mul(const Matrix* a, const Matrix* b) {
    Matrix result = matrix_create(a->rows, b->cols);
    for (size_t i = 0; i < a->rows; ++i) {
        for (size_t j = 0; j < b->cols; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < a->cols; ++k) {
                sum += matrix_get(a, i, k) * matrix_get(b, k, j);
            }
            matrix_set(&result, i, j, sum);
        }
    }
    return result;
}

Matrix matrix_transpose(const Matrix* mat) {
    Matrix result = matrix_create(mat->cols, mat->rows);
    for (size_t i = 0; i < mat->rows; ++i) {
        for (size_t j = 0; j < mat->cols; ++j) {
            matrix_set(&result, j, i, matrix_get(mat, i, j));
        }
    }
    return result;
}

Matrix matrix_inverse(const Matrix* mat) {
    if (mat->rows != mat->cols) {
        // 非方阵无法求逆
        Matrix empty = {0, 0, NULL};
        return empty;
    }

    size_t n = mat->rows;
    Matrix augmented = matrix_create(n, n * 2);

    // 构建增广矩阵 [A | I]
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            matrix_set(&augmented, i, j, matrix_get(mat, i, j));
        }
        matrix_set(&augmented, i, n + i, 1.0);  // 单位矩阵部分
    }

    // 高斯-约旦消元
    for (size_t i = 0; i < n; ++i) {
        // 找到最大值进行主元交换
        double max_val = fabs(matrix_get(&augmented, i, i));
        size_t max_row = i;
        for (size_t k = i + 1; k < n; ++k) {
            double val = fabs(matrix_get(&augmented, k, i));
            if (val > max_val) {
                max_val = val;
                max_row = k;
            }
        }

        if (max_val < EPSILON) {
            matrix_free(&augmented);
            Matrix empty = {0, 0, NULL};
            return empty;  // 不可逆
        }

        // 行交换
        if (max_row != i) {
            for (size_t j = 0; j < 2 * n; ++j) {
                double tmp = matrix_get(&augmented, i, j);
                matrix_set(&augmented, i, j, matrix_get(&augmented, max_row, j));
                matrix_set(&augmented, max_row, j, tmp);
            }
        }

        // 将对角线归一化
        double diag = matrix_get(&augmented, i, i);
        for (size_t j = 0; j < 2 * n; ++j) {
            double val = matrix_get(&augmented, i, j) / diag;
            matrix_set(&augmented, i, j, val);
        }

        // 将其他行该列置为0
        for (size_t k = 0; k < n; ++k) {
            if (k == i) continue;
            double factor = matrix_get(&augmented, k, i);
            for (size_t j = 0; j < 2 * n; ++j) {
                double val = matrix_get(&augmented, k, j) - factor * matrix_get(&augmented, i, j);
                matrix_set(&augmented, k, j, val);
            }
        }
    }

    // 提取右边的逆矩阵
    Matrix inverse = matrix_create(n, n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            matrix_set(&inverse, i, j, matrix_get(&augmented, i, j + n));
        }
    }

    matrix_free(&augmented);
    return inverse;
}


double matrix_determinant(const Matrix* mat) {
    if (mat->rows != mat->cols) {
        // 非方阵无法计算行列式
        return NAN;  // 返回NaN表示错误
    }

    size_t n = mat->rows;
    Matrix augmented = matrix_create(n, n);

    // 创建矩阵副本，以避免改变原矩阵
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            matrix_set(&augmented, i, j, matrix_get(mat, i, j));
        }
    }

    double det = 1.0;

    // 高斯消元
    for (size_t i = 0; i < n; ++i) {
        // 找到最大值进行主元交换
        double max_val = fabs(matrix_get(&augmented, i, i));
        size_t max_row = i;
        for (size_t k = i + 1; k < n; ++k) {
            double val = fabs(matrix_get(&augmented, k, i));
            if (val > max_val) {
                max_val = val;
                max_row = k;
            }
        }

        if (max_val < EPSILON) {
            matrix_free(&augmented);
            return 0.0;  // 行列式为0，矩阵不可逆
        }

        // 行交换
        if (max_row != i) {
            for (size_t j = 0; j < n; ++j) {
                double tmp = matrix_get(&augmented, i, j);
                matrix_set(&augmented, i, j, matrix_get(&augmented, max_row, j));
                matrix_set(&augmented, max_row, j, tmp);
            }
            // 每次交换行时，行列式的符号要反转
            det = -det;
        }

        // 将对角线归一化
        double diag = matrix_get(&augmented, i, i);
        for (size_t j = 0; j < n; ++j) {
            double val = matrix_get(&augmented, i, j) / diag;
            matrix_set(&augmented, i, j, val);
        }

        // 将其他行该列置为0
        for (size_t k = i + 1; k < n; ++k) {
            double factor = matrix_get(&augmented, k, i);
            for (size_t j = 0; j < n; ++j) {
                double val = matrix_get(&augmented, k, j) - factor * matrix_get(&augmented, i, j);
                matrix_set(&augmented, k, j, val);
            }
        }
    }

    // 计算行列式为对角线元素的乘积
    for (size_t i = 0; i < n; ++i) {
        det *= matrix_get(&augmented, i, i);
    }

    matrix_free(&augmented);
    return det;
}

Matrix matrix_identity(size_t size) {
    Matrix mat = matrix_create(size, size);
    for (size_t i = 0; i < size; ++i) {
        matrix_set(&mat, i, i, 1.0);
    }
    return mat;
}

void matrix_print(const Matrix* mat) {
    for (size_t i = 0; i < mat->rows; ++i) {
        for (size_t j = 0; j < mat->cols; ++j) {
            printf("%8.4f ", matrix_get(mat, i, j));
        }
        printf("\n");
    }
}
