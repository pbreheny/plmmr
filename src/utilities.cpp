#include "utilities.h"


// Product filebacked matrix X with vector y
// ind_col and p_length let select elements of X be used
double prod(XPtr<BigMatrix> X_, double *y, int i, int p) {
  MatrixAccessor<double> X(*X_);
  double val=0;
  for (int j=0;j < p;j++) {
    val += X[i][j] * y[j];
  }
  return(val);
}

// Cross product of y with jth column of filebacked matrix X
// ind_col and p_length let select elements of X be used
double crossprod(XPtr<BigMatrix> X_, double *y_, int j, int n) {
  MatrixAccessor<double> X(*X_);
  double *col = X[j];
  double val=0.0;
  for (int i=0;i < n;i++) {
    val += col[i]*y_[i];
  }
 // Rprintf("crossprod's val is %f\n", val);

  return(val);
}

// column-wise sums of squares for filebacked matrix X
NumericVector mean_sqsum(XPtr<BigMatrix> X_, int n, int p){
  MatrixAccessor<double> X(*X_);
  NumericVector sqsum(p);

  for (int j=0;j<p;j++){
    for (int i=0;i<n;i++) sqsum[j] += pow(X[j][i], 2);
    sqsum[j] = sqsum[j]/n; // take mean
  }
  return sqsum;
}

// columnwise sum of squares of filebacked matrix X
NumericVector colwise_l2mean(XPtr<BigMatrix> X_, int n, int p) {
  MatrixAccessor<double> X(*X_);
  NumericVector sqsum(p);
  NumericVector scale_vals(p);
  for (int j=0;j<p;j++){
    for (int i=0;i<n;i++) sqsum[j] += pow(X[j][i], 2);
    scale_vals[j] = sqrt(sqsum[j]/n);
  }

  return scale_vals;
}

// scale columns of a filebacked matrix X
void rescale_cols(XPtr<BigMatrix> X_, int n, int p, NumericVector scales) {
  MatrixAccessor<double> X(*X_);
  for (int j=0;j<p;j++){
    for (int i=0;i<n;i++){
      X[j][i] = X[j][i]/scales[j];
    }
  }
}

// Gaussian loss
double g_loss(double *r, int n) {
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
}

// standardize a filebacked matrix
void standardize_and_get_residual(NumericVector &center, NumericVector &scale,
                                  int *p_keep_ptr, vector<int> &col_idx, //columns to keep, removing columns whose scale < 1e-6
                                  vector<double> &z, double *lambda_max_ptr,
                                  int *xmax_ptr, XPtr<BigMatrix> xMat, double *y,
                                  int *row_idx, double alpha, int n, int p) {
  MatrixAccessor<double> xAcc(*xMat);
  double *xCol;
  double sum_xy, sum_y;
  double zmax = 0.0, zj = 0.0;
  int i, j;

  for (j = 0; j < p; j++) {
    xCol = xAcc[j];
    sum_xy = 0.0;
    sum_y = 0.0;

    for (i = 0; i < n; i++) {
      center[j] += xCol[row_idx[i]];
      scale[j] += pow(xCol[row_idx[i]], 2);

      sum_xy = sum_xy + xCol[row_idx[i]] * y[i];
      sum_y = sum_y + y[i];
    }

    center[j] = center[j] / n; //center
    scale[j] = sqrt(scale[j] / n - pow(center[j], 2)); //scale

    if (scale[j] > 1e-6) {
      col_idx.push_back(j);
      zj = (sum_xy - center[j] * sum_y) / (scale[j] * n); //residual
      if (fabs(zj) > zmax) {
        zmax = fabs(zj);
        *xmax_ptr = j; // xmax_ptr is the index in the raw xMat, not index in col_idx!
      }
      z.push_back(zj);
    }
  }
  *p_keep_ptr = col_idx.size();
  *lambda_max_ptr = zmax / alpha;
}