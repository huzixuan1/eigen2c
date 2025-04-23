#ifndef MATRIX_DECOMP_H
#define MATRIX_DECOMP_H

#include "matrix_ops.h"

double matrix_LU_decomposition(const Matrix* A, Matrix* L, Matrix* U);
double matrix_QR_decomposition(const Matrix* A, Matrix* Q, Matrix* R);
double matrix_cholesky_decomposition(const Matrix* A, Matrix* L);
double matrix_svd(const Matrix* A, Matrix* U, Matrix* S, Matrix* V);

#endif  // MATRIX_DECOMP_H

