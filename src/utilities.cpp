#include "utilities.h"

// Cross product of y with jth column of filebacked matrix X
double crossprod(XPtr<BigMatrix> X_, double *y_, int j, int n) {
  MatrixAccessor<double> X(*X_);
  double *xCol = X[j];
  double val=0.0;
  int i;

// #pragma omp parallel for private(i) schedule(static)

  for (i=0;i < n;i++) {
    val += xCol[i]*y_[i];
  }

  return(val);
}

// columnwise means of filebacked matrix X
NumericVector col_means(XPtr<BigMatrix> X_, int n, int p){
  MatrixAccessor<double> X(*X_);
  NumericVector center_vals(p);
  int j;

#pragma omp parallel for private(j) schedule(static)

  for (j=0;j<p;j++) {
    double *xCol = X[j];
    double sum_j = 0;
    for (int i=0;i<n;i++) sum_j += xCol[i];
    center_vals[j] = sum_j/n;
    xCol = nullptr;
  }

  return center_vals;
}

// columnwise sum of squares of filebacked matrix X
NumericVector colwise_l2mean(XPtr<BigMatrix> X_, int n, int p) {
  MatrixAccessor<double> X(*X_);
  NumericVector scale_vals(p);
  int j;

#pragma omp parallel for private(j) schedule(static)

  for (j=0;j<p;j++) {
    double *xCol = X[j];
    double sqsum = 0;
    for (int i=0;i<n;i++) sqsum += pow(xCol[i], 2);
    scale_vals[j] = sqrt(sqsum/n);
    xCol = nullptr;
  }

  return scale_vals;
}

// center columns of a filebacked matrix X
void center_cols(XPtr<BigMatrix> X_, int n, int p, NumericVector centers) {
  MatrixAccessor<double> X(*X_);
  int j;

  #pragma omp parallel for private(j) schedule(static)

  for (j=0;j<p;j++){
    double *xCol = X[j];
    for (int i=0;i<n;i++){
      xCol[i] = xCol[i] - centers[j];
    }
    xCol = nullptr;
  }
}

// scale columns of a filebacked matrix X
void scale_cols(XPtr<BigMatrix> X_, int n, int p, NumericVector scales) {
  MatrixAccessor<double> X(*X_);
  int j;

#pragma omp parallel for private(j) schedule(static)

  for (j=0;j<p;j++){
    double *xCol = X[j];
    for (int i=0;i<n;i++){
      if (scales[j] < 1e-3) {
        xCol[i] = 0;
      } else {
        xCol[i] = xCol[i]/scales[j];
      }
    }
    xCol = nullptr;
  }
}

// column-wise standard deviation of a (*centered*), filebacked matrix X
NumericVector sd(XPtr<BigMatrix> centered_X_, int n, int p){
  MatrixAccessor<double> X(*centered_X_);
  NumericVector sd_vals(p);
  int j;

//#pragma omp parallel for private(j) schedule(static)

  for (j=0;j<p;j++) {
    double *xCol = X[j];
    double sum_j = 0;
    for (int i=0;i<n;i++) sum_j += pow(xCol[i],2);
    sd_vals[j] = sqrt(sum_j/(n-1));
  }

  return(sd_vals);
}
