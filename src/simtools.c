/*  ==========================================================================  *
 *  Miguel de Carvalho                                                          *
 *  Copyright (C) 2014                                                          *
 *  --------------------------------------------------------------------------  *
 *  This program is free software; you can redistribute it and/or modify        *
 *  it under the terms of the GNU General Public License as published by        *
 *  the Free Software Foundation; either version 2 of the License, or           *
 *  (at your option) any later version.                                         *
 *                                                                              *
 *  This program is distributed in the hope that it will be useful,             *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
 *  GNU General Public License for more details.                                *
 *                                                                              *
 *  You should have received a copy of the GNU General Public License           *
 *  along with this program; if not, a copy is available at                     *
 *  http://www.r-project.org/Licenses/                                          *
 *  ==========================================================================  *
 */
 
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "simtools.h" 

void _mean(double* x, int *n, double* out) {
  /* computes the mean of a vector */
  int i;
  *out = 0;
  for(i = 0; i < *n; i++) 
    *out += *x++;
  *out /= *n;
}
void _var(double* x, int* n, double* out) {
  /* computes the variance of a vector */
  int i;
  double mean;
  _mean(x, n, &mean);
  out[0] = 0;
  for(i = 0; i < *n; i++) 
    *out += pow(*x++, 2);
  *out = *out/(double)(*n - 1) - pow(mean, 2) *(*n / (double)(*n - 1));
}
void _sd(double* x, int *n, double* out) {
  /* computes the standard deviation of a vector */
  double var;
  _var(x, n, &var);
  *out = sqrt(var);
}
void _sum(double* x, int *n, double* out) {
  /* sum of vector elements */
  int i;
  *out = 0;
  for(i = 0; i < *n; i++) 
    *out += *x++;
}
void _cumsum(double* x, int *n, double* out) {
  /* cumulative sum of a vector's elements */
  int i;
  *out = x[0];
  for(i = 1; i < *n; i++) 
    out[i] = out[i - 1] + x[i];
}
void _cumprod(double* x, int* n, double* out) {
  /* cumulative product of the entries in a vector x
     of size n; output is saved in out*/
  int i;
  out[0] = x[0];
  for(i = 1; i < *n; i++)
    out[i] = out[i - 1] * x[i];
}
void _prodcum(double *x, int* n, double* out) {
  /* cumulative product of the one minus each entry of a vector x  
     of size n; output is saved in out*/
  int i;
  out[0] = 1 - x[0];
  for(i = 1; i < *n; i++)
    out[i] = out[i - 1] * (1 - x[i]);
}
void _rowSums(size_t nrow, size_t ncol, double A[nrow][ncol], 
	  double* out) {
  /* computes row sums of a matrix A; stores output in the 
     vector out */
  int i, j;
  for(i = 0; i < nrow ; i++) {
    out[i] = 0;
    for(j = 0; j < ncol; j++)
      out[i] += A[i][j];
  }
}
void _colSums(size_t nrow, size_t ncol, double A[nrow][ncol], 
	  double out[ncol]) {
  /* computes column sums of a matrix A; stores output in the 
     vector out */
  int i, j;
  for(i = 0; i < nrow ; i++) {
    out[i] = 0;
    for(j = 0; j < ncol; j++)
      out[i] += A[j][i];
  }
}
void _rowMeans(size_t nrow, size_t ncol, double A[nrow][ncol],
	  double out[nrow]) {
  /* computes column means of a matrix A; stores output in the
     vector out */
  int i, j;
  for(i = 0; i < nrow ; i++) {
    out[i] = 0;
    for(j = 0; j < ncol; j++)
      out[i] += A[i][j];
    out[i] /= (double)nrow;
  }
}
void _colMeans(size_t nrow, size_t ncol, double A[nrow][ncol],
	  double out[ncol]) {
  /* computes column means of a matrix A; stores output in the
     vector out */
  int i, j;
  for(i = 0; i < nrow ; i++) {
    out[i] = 0;
    for(j = 0; j < ncol; j++)
      out[i] += A[j][i];
    out[i] /= (double)nrow;
  }  
}
void _rep(double* x, int* n, double* out) {
  /* initializes a vector of size n with all entries equal to x*/
  int i;
  for (i = 0; i < *n; i++) 
    out[i] = *x;
}
void _REP(double* x, size_t nrow, size_t ncol, double out[nrow][ncol]) {
/* computes column means of a matrix A; stores output in the
     vector out */
  int i, j;
  for(i = 0; i < nrow ; i++) {
    for(j = 0; j < ncol; j++)
      out[i][j] = *x;
  }  
}
void SUM(double* x, double* y, int* n, double* out) {
  /* sums two vectors x and y of size n; output is saved in out*/
  int i;
  for (i = 0; i < *n; i++) 
    out[i] = x[i] + y[i];
}
void _rMultinom(size_t nrow, size_t ncol, double P[nrow][ncol],
		int* out) {
  /* this is a version of rMultinom from the R package Hmisc (m = 1) */
  int i, j;
  int count[nrow];
  int lev[ncol];
  double U[nrow][ncol];
  double un[nrow][ncol];
  
  for(i = 0; i < ncol; i++) 
    lev[i] = i + 1;
  for(i = 0; i < nrow; i++) {
    count[i] = 0;
    U[i][0] = P[i][0];
    un[i][0] = runif(0, 1);
    if(un[i][0] > U[i][0])
      count[i] += 1; 
    for(j = 1; j < ncol; j++) {
      U[i][j] = U[i][j - 1] + P[i][j];
      un[i][j] = un[i][0];
      if(un[i][j] > U[i][j])
	count[i] += 1; 
    }    
  }  
  for(i = 0; i < nrow; i++) {
    out[i] = lev[count[i]];  
  }
}
