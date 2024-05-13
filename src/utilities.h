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

double crossprod(double *X, double *y, int n, int j);

// Sum of squares of jth column of X
double sqsum(double *X, int n, int j);

double g_loss(double *r, int n);

#endif