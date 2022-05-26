/*  ==========================================================================  *
 *  Miguel de Carvalho                                                          *
 *  Copyright (C) 2022                                                          *
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
 *  --------------------------------------------------------------------------  *
 *  The componentwise adaptive Metropolis–Hastings implemented below follows    *
 *  the steps:                                                                  *
 *    1) Enumerate the indices of the Bernstein polynomial;                     *      
 *    2) Sample weights with componentwise adaptive MH (Haario et al);          *
 *    3) Compute trajectory of the target of interest (density/cdf).            *
 *  Step 2 proceeds as follows:                                                 * 
 *    1 <= gibbs <= burn: use static MH update;                                 *          
 *    burn < gibbs < 2 * burn: start building adaptive variance;                * 
 *    2 * burn < gibbs: use adaptive variance and harvest iterates;             *   
 *                                                                              *  
 *  For further details on the method see the references below.                 *
 *                                                                              *
 *  Main variables                                                              *
 *   - w: vector of pseudo-angles;                                              *
 *   - k: number of pseudo-angles;                                              *
 *   - c: hyperparameter of Dirichlet prior for weights ~ Dir(c, ..., c);       * 
 *   - p: number of parameters of p-dimensional Dirichlet kernel;               * 
 *   - m: number of basis functions;                                            * 
 *   - J: row sum of Bernstein polynomial index matrix (indices below);         *                        
 *   - T: number of MCMC iterations;                                            *
 *   - burn: number of burn in iterations;                                      * 
 *   - grid: grid on which the posterior density is to be evaluated;            *
 *   - lgrid: length of grid;                                                   * 
 *   - target: none (0), or angular density (1), or angular measure (2);        * 
 *                                                                              *    
 *  Outputs                                                                     *
 *  - weights: (T * m)-matrix with the posterior simulated weights;             *
 *  - indices: (m * p)-matrix with indices of Bernstein polynomials (vertex     *  
 *    coefficients appear on the first p rows)                                  *   
 *  - traj: (T * lgrid)-matrix with trajectories of densities evaluated on      *
 *    a grid of length lgrid;                                                   *
 *                                                                              *
 *  References                                                                  *
 *  - Haario, H., Eero, S. and Johanna, T. (2005) Componentwise adaptation for  *
 *    high dimensional MCMC. Computational Statistics, 20, 265-73.              *
 *                                                                              *
 *  - Hanson, T., de Carvalho, M. and Chen, Y. (2017) Bernstein polynomial      *
 *    angular densities of multivariate extreme value distributions,            *             
 *    Statistics and Probability Letters, 128, 60-66.                           *
 *                                                                              *
 *  ==========================================================================  *
 */

#include "header.h"

/* AUXILIARY FUNCTIONS */
double dir(int p, double alpha[p], double win[p]) {
  int i;
  double result;
  result = 0;
  for (i = 0; i < p; i++) 
    result += alpha[i];
  result = lgamma(result);
  for (i = 0; i < p; i++)
    result += -lgamma(alpha[i]) + (alpha[i] - (double)(1)) * log(win[i]);
  return(exp(result));
  return(result);
}

static void bern(double* w, int* k, double* c, int*p, int*m, int*J, int* T,
		 int* burn, double* grid, int* lgrid, double* weights,
		 double* indices, double* traj, int* target) {
  
  /* 0) ALLOCATE AND INITIALIZE VARIABLES */
  int aux, flag, gibbs, h, i, j, l, t; 
  int r[*p];
  int ind[*m][*p];
  double dm[*k][*m];
  double W[*k][*p];
  double Wr[*m], indr[*m], pip[*m], pipc[*m]; 
  double mean[*m - *p], v[*m - *p], var[*m - *p], tally[*m - *p];
  double adaptvar, ll, llold, llnew, po, s, sum, vold;
  double auxv[*k];
  
  llold = -DBL_MAX;  
  for (i = 0; i < *p; i++)
    r[i] = 1;  
  for (i = 0; i < *m; i++)
    for (j = 0; j < *p; j++)
      ind[i][j] = 0;  
  for (i = 0; i < *T * (*lgrid); i++)
    traj[i] = 0;    
  for (i = 0; i < *k; i++)
    for (j = 0; j < *p; j++)
      W[i][j] = w[i + j * *k];  
  for (i = 0; i < *m; i++) 
    pip[i] = 1;     
  for(i = 0; i < *m - *p; i++) {
    v[i] = 0; tally[i] = 0; mean[i] = 0; var[i] = 0;
  }
  
  /* 1) ENUMERATE INDICES OF BERNSTEIN POLYNOMIAL */
  r[0] = *J - *p + 1;
  t = *J - *p + 1;
  h = -1;
  for (j = 0; j < *p; j++)
    ind[0][j] = r[j];
  for (i = 1; i < *m; i++) {
    if (t != 2)
      h = -1;
    h++; t = r[h]; r[h] = 1; r[0] = t - 1; r[h + 1] = r[h + 1] + 1;
    for (j = 0; j < *p; j++) 
      ind[i][j] = r[j];
  }
  
  /* put vertex coefficients into first p slots */
  aux = 0; 			
  for (j = 1; j < *p; j++) {	
    aux += (int)(choose(*J - *p - 2 + j + 1, *J - *p - 1));
    for (l = 0; l < *p; l++) {
      r[l] = ind[j][l];
      ind[j][l] = ind[aux][l];
      ind[aux][l] = r[l];
    }
  }
  for (i = 0; i < *p; i++)
    for (j = 0; j < *m; j++)
      indices[*m * i  + j] = ind[j][i];
  
  /* 2) SAMPLE WEIGHTS WITH COMPONENTWISE ADAPTIVE MH (Haario et al) */
  /* Density Matrix */
  for (i = 0; i < *k; i++) {
    for (j = 0; j < *m; j++) {
      for (l = 0; l < *p; l++) {
	Wr[l] = W[i][l];
	indr[l] = (double)(ind[j][l]);
      }
      dm[i][j] = dir(*p, indr, Wr);
    }
  }   
  /* Sample */
  Rprintf("MCMC iterations:\n================\n");
  for (gibbs = 1; gibbs <= *T; gibbs++) {
    if (gibbs > *burn) {
      for (j = 0; j < *m - *p; j++) {
      	mean[j] += v[j];
	var[j] += pow(v[j], 2);
      }
    }
    for (j = 0; j < *m - *p; j++) {
      vold = v[j];
      if (gibbs > *burn * 2) {
	po = pow(mean[j] / (double)(gibbs - *burn), 2);	
	adaptvar = (double)(2.4) * 
	  (var[j] / (double)(gibbs - *burn) - po + (double)(0.001));
      } else {
	adaptvar = (double)(1);
      }
      v[j] += rnorm(0, sqrt(adaptvar));
      /* Evaluate log-likelihood */
      ll = 0; flag = 0;
      s = (double)(*p);
      for (i = 0; i < *m - *p; i++)
	s += exp(v[i]);
      for (i = *p; i < *m; i++)
	pip[i] = exp(v[i - *p]) / s;
      for (i = 0; i < *p; i++)
	pip[i] = (double)(1) / (double)(*p);
      
      for (i = *p; i < *m; i++)
	for (l = 0; l < *p; l++)
	  pip[l] += - pip[i] * ((double)(ind[i][l]) -
				(double)(1)) / (double)(*J - *p);
      
      for (i = 0; i < *p; i++)
	if (pip[i] < 0 || pip[i] > 1)
	  flag = 1;
      if (flag == 0) {
	ll = 0;
	for(i = 0; i < *m; i++)
	  ll += *c * log(pip[i]);
	for (i = 0; i < *k; i++) {
	  auxv[i] = 0;
	  for (l = 0; l < *m; l++)
	    auxv[i] += dm[i][l] * pip[l];	  
	}
	sum = 0;
	for (i = 0; i < *k; i++)
	  sum += log(auxv[i]);
	ll += sum;
      }      
      /* Metropolis step  */
      llnew = ll;
      if (llnew != 0) {
	if (log(runif(0, 1)) < (llnew - llold)) {
	  tally[j] = tally[j] + 1;
	  llold = llnew;
	  for (i = 0; i < *m; i++)
	    pipc[i] = pip[i];	
	} else {
	  v[j] = vold;
	}
      } else {
	v[j] = vold;
      }
    }
    if((gibbs + 1) % 1000 == 0)
      Rprintf("%d\n", (gibbs + 1));    
    for (i = 0; i < *m; i++) {
      weights[(gibbs - 1) + *T * i] = pipc[i];
    }
  }
  
  /* 3) COMPUTE DENSITY/ANGULAR MEASURE */
  /* Angular density */
  if (*target == 1) {
    Rprintf("Preparing outputs...\n");
    for (t = 0; t < *T; t++)
      for (j = 0; j < *lgrid; j++) 
	for (h = 0; h < *m; h++)
	  traj[t + *T * j] += weights[t + *T * h] *
	    dbeta(grid[j], (double)(ind[h][0]), (double)(ind[h][1]), 0);
  }
  /* Angular measure */
  if (*target == 2) {
    Rprintf("Preparing outputs...\n");
    for (t = 0; t < *T; t++)
      for (j = 0; j < *lgrid; j++) 
	for (h = 0; h < *m; h++)
	  traj[t + *T * j] += weights[t + *T * h] *
	    pbeta(grid[j], (double)(ind[h][0]), (double)(ind[h][1]), 1, 0);    
  }
  Rprintf("DONE\n");    
}

SEXP BERN(SEXP w, SEXP k, SEXP c, SEXP p, SEXP m, SEXP J, SEXP T,
	  SEXP burn, SEXP grid, SEXP lgrid, SEXP target) {
  int nprot = 0;
  int _k = asInteger(k);
  double _c = asReal(c);
  int _p = asInteger(p);
  int _m = asInteger(m); 
  int _J = asInteger(J);
  int _T = asInteger(T);
  int _burn = asInteger(burn);
  int _lgrid = asInteger(lgrid);
  int _target = asInteger(target);
  SEXP indices, weights, traj;
  
  w = PROTECT(coerceVector(w, REALSXP)); nprot++;
  grid = PROTECT(coerceVector(grid, REALSXP)); nprot++;
  PROTECT(weights = allocMatrix(REALSXP, _T, _m)); nprot++;
  PROTECT(indices = allocMatrix(REALSXP, _m, _p)); nprot++;  
  PROTECT(traj = allocMatrix(REALSXP, _T, _lgrid)); nprot++;
  
  double* weightsptr; weightsptr = REAL(weights);
  double* indptr; indptr = REAL(indices);
  double* trajptr; trajptr = REAL(traj);
  
  GetRNGstate();
  
  bern(REAL(w), &_k, &_c, &_p, &_m, &_J, &_T, &_burn, REAL(grid),
       &_lgrid, weightsptr, indptr, trajptr, &_target);
  
  PutRNGstate();
  
  SEXP ans = PROTECT(allocVector(VECSXP, 3)); nprot++;
  SET_VECTOR_ELT(ans, 0, weights);
  SET_VECTOR_ELT(ans, 1, indices);
  SET_VECTOR_ELT(ans, 2, traj);

  SEXP nm = allocVector(STRSXP, 3);
  setAttrib(ans, R_NamesSymbol, nm);
  SET_STRING_ELT(nm, 0, mkChar("weights"));
  SET_STRING_ELT(nm, 1, mkChar("indices"));
  SET_STRING_ELT(nm, 2, mkChar("traj"));
 
  UNPROTECT(nprot);
  return(ans);
}


