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
#endif