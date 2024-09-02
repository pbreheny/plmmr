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

// double prod(XPtr<BigMatrix> X_, double *y, int i, int p);
double crossprod(XPtr<BigMatrix> X_, double *y, int j, int n);
NumericVector mean_sqsum(XPtr<BigMatrix> X_, int n, int p);
NumericVector col_means(XPtr<BigMatrix> X_, int n, int p);
NumericVector colwise_l2mean(XPtr<BigMatrix> X_, int n, int p);
void center_cols(XPtr<BigMatrix> X_, int n, int p, NumericVector centers);
void scale_cols(XPtr<BigMatrix> X_, int n, int p, NumericVector scales);
NumericVector sd(XPtr<BigMatrix> centered_X_, int n, int p);
double g_loss(double *r, int n);
// void standardize_and_get_residual(NumericVector &center, NumericVector &scale,
//                                   int *p_keep_ptr, vector<int> &col_idx, //columns to keep, removing columns whose scale < 1e-6
//                                   vector<double> &z, double *lambda_max_ptr,
//                                   int *xmax_ptr, XPtr<BigMatrix> xMat, double *y,
//                                   int *row_idx, double alpha, int n, int p);
// void fill_in(XPtr<BigMatrix> fill_into_,
//              int rows_into,
//              int cols_into,
//              NumericMatrix fill_from,
//              int rows_from,
//              int cols_from);
// void fill_in_filebacked(XPtr<BigMatrix> fill_into_,
//                         int rows_into,
//                         int cols_into,
//                         XPtr<BigMatrix> fill_from_,
//                         int rows_from,
//                         int cols_from,
//                         int skip);
#endif