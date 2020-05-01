#include <R.h>
#include <stddef.h>
#pragma once
void _mean(double* x, int *n, double* out);
void _var(double* x, int *n, double* out);
void _sd(double* x, int *n, double* out);
void _sum(double* x, int *n, double* out);
void _cumsum(double* x, int *n, double* out); 
void _cumprod(double *x, int* n, double* out);
void _prodcum(double *x, int* n, double* out);
void _rowSums(size_t nrow, size_t ncol, double A[nrow][ncol], double* out);
void _colSums(size_t nrow, size_t ncol, double A[nrow][ncol], double out[ncol]);
void _rowMeans(size_t nrow, size_t ncol, double A[nrow][ncol], 
	       double out[nrow]);
void _colMeans(size_t nrow, size_t ncol, double A[nrow][ncol], 
	       double out[ncol]);
void _rep(double* x, int* n, double* out);
void _REP(double* x, size_t nrow, size_t ncol, double out[nrow][ncol]);
void SUM(double* x, double* y, int* n, double* out);
void _rMultinom(size_t nrow, size_t ncol, double P[nrow][ncol], int* out);
void SampleReplace(int k, int n, int *y);
