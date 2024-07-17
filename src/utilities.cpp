#include "utilities.h"


// Product filebacked matrix X with vector y
// NOTE: this function has *not* been optimized for memory management -- needs work!
// double prod(XPtr<BigMatrix> X_, double *y, int i, int p) {
//   MatrixAccessor<double> X(*X_);
//   double val=0;
//   for (int j=0;j < p;j++) {
//     val += X[i][j] * y[j];
//   }
//   return(val);
// }

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
 // Rprintf("crossprod's val is %f\n", val);

  return(val);
}

// column-wise sums of squares for filebacked matrix X
NumericVector mean_sqsum(XPtr<BigMatrix> X_, int n, int p){
  MatrixAccessor<double> X(*X_);
  NumericVector sqsum(p);
  int j;

#pragma omp parallel for private(j) schedule(static)

  for (j=0;j<p;j++){
    double *xCol = X[j];
    for (int i=0;i<n;i++) sqsum[j] += pow(xCol[i], 2);
    sqsum[j] = sqsum[j]/n; // take mean
    xCol = nullptr;
  }
  return sqsum;
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
      xCol[i] = xCol[i]/scales[j];
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

// Gaussian loss
double g_loss(double *r, int n) {
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
}

// standardize a filebacked matrix
// NOTE: this function has *not* been optimized for memory management -- needs work!
// void standardize_and_get_residual(NumericVector &center, NumericVector &scale,
//                                   int *p_keep_ptr, vector<int> &col_idx, //columns to keep, removing columns whose scale < 1e-6
//                                   vector<double> &z, double *lambda_max_ptr,
//                                   int *xmax_ptr, XPtr<BigMatrix> xMat, double *y,
//                                   int *row_idx, double alpha, int n, int p) {
//   MatrixAccessor<double> xAcc(*xMat);
//   double *xCol;
//   double sum_xy, sum_y;
//   double zmax = 0.0, zj = 0.0;
//   int i, j;
//
//   for (j = 0; j < p; j++) {
//     xCol = xAcc[j];
//     sum_xy = 0.0;
//     sum_y = 0.0;
//
//     for (i = 0; i < n; i++) {
//       center[j] += xCol[row_idx[i]];
//       scale[j] += pow(xCol[row_idx[i]], 2);
//
//       sum_xy = sum_xy + xCol[row_idx[i]] * y[i];
//       sum_y = sum_y + y[i];
//     }
//
//     center[j] = center[j] / n; //center
//     scale[j] = sqrt(scale[j] / n - pow(center[j], 2)); //scale
//
//     if (scale[j] > 1e-6) {
//       col_idx.push_back(j);
//       zj = (sum_xy - center[j] * sum_y) / (scale[j] * n); //residual
//       if (fabs(zj) > zmax) {
//         zmax = fabs(zj);
//         *xmax_ptr = j; // xmax_ptr is the index in the raw xMat, not index in col_idx!
//       }
//       z.push_back(zj);
//     }
//   }
//   *p_keep_ptr = col_idx.size();
//   *lambda_max_ptr = zmax / alpha;
// }


// fill in the values of a filebacked matrix with values from a numeric matrix
// NOTE: this function has *not* been optimized for memory management -- needs work!
// void fill_in(XPtr<BigMatrix> fill_into_,
//              int rows_into,
//              int cols_into,
//              NumericMatrix fill_from,
//              int rows_from,
//              int cols_from){
//   MatrixAccessor<double> fill_into(*fill_into_);
//   Rprintf("\nMade it past fill_in() declarations");
//
// #pragma omp parallel for
//
//   for (int i = 0; i < rows_from; ++i) {
//     for (int j = 0; j < cols_from; ++j) {
//       double *intoCol = fill_into[j];
//       intoCol[i] = fill_from[j * rows_from + i];
//     }
//   }
//
// }

// fill in the values of one filebacked matrix with values from another filebacked matrix
// NOTE: this function has *not* been optimized for memory management -- needs work!
// void fill_in_filebacked(XPtr<BigMatrix> fill_into_,
//              int rows_into,
//              int cols_into,
//              XPtr<BigMatrix> fill_from_,
//              int rows_from,
//              int cols_from,
//              int skip){
//   MatrixAccessor<double> fill_into(*fill_into_);
//   MatrixAccessor<double> fill_from(*fill_from_);
// Rprintf("\nMade it past fill_in_filebacked() declarations");
//   for (int i=0;i<rows_from;i++){
//     for (int j=0;j<cols_from;j++){
//       fill_into[j + skip][i] = fill_from[j][i];
//     }
//   }
//
// }