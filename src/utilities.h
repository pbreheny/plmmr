#include <RcppArmadillo.h>
#include "bigmemory/BigMatrix.h"
#include <time.h>
#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/bigmemoryDefines.h"

#include "plmmr_omp.h"

#ifndef UTILITIES_H
#define UTILITIES_H

using namespace Rcpp;
using namespace std;

double prod(XPtr<BigMatrix> X_, double *y, int i, int p);
double crossprod(XPtr<BigMatrix> X_, double *y, int j, int n);
NumericVector mean_sqsum(XPtr<BigMatrix> X_, int n, int p);
NumericVector colwise_l2mean(XPtr<BigMatrix> X_, int n, int p);
void rescale_cols(XPtr<BigMatrix> X_, int n, int p, NumericVector scales);
double g_loss(double *r, int n);
void standardize_and_get_residual(NumericVector &center, NumericVector &scale,
                                  int *p_keep_ptr, vector<int> &col_idx, //columns to keep, removing columns whose scale < 1e-6
                                  vector<double> &z, double *lambda_max_ptr,
                                  int *xmax_ptr, XPtr<BigMatrix> xMat, double *y,
                                  int *row_idx, double alpha, int n, int p);
#endif