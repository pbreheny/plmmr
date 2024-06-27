#include "utilities.h"

// rotating filebacked data
RcppExport SEXP rotate_filebacked(SEXP std_X_,
                                  SEXP WUt_,
                                  SEXP rot_X_,
                                  SEXP ncore_) {
  //Rprintf("\n declarations");
  // declarations: input
  XPtr<BigMatrix> std_X(std_X_); // this points to the filebacked matrix of standardized data
  NumericMatrix WUt(WUt_); // this is the in-memory n x r projection matrix
  XPtr<BigMatrix> rot_X(rot_X_); // this points to the *empty* filebacked matrix; will hold result
  int n = std_X->nrow(); // n = number of observations
  int p = std_X->ncol(); // p = number of features (including non-genomic predictors!)
  int r = WUt.nrow(); // r = number of eigenvalues used in rotation (r stands for 'rank')
  // Rprintf("\n nrow(std_X) %d\n", n);
  // Rprintf("\n ncol(std_X) %d\n", p);
  // Rprintf("\n nrow(WUt) %d\n", r);

 //Rprintf("\n initialize accessors");
  // initialize MatrixAccessors with std_X and rot_X
  MatrixAccessor<double> std_X_acc(*std_X); // initialize MatrixAccessor with std_X
  MatrixAccessor<double> rot_X_acc(*rot_X);

  // set up omp
  //Rprintf("\n set up OpenMP");
  int useCores = INTEGER(ncore_)[0];
#ifdef PLMMR_OMP_H_
  int haveCores = omp_get_num_procs();
  if(useCores < 1) {
    useCores = haveCores;
  }
  omp_set_dynamic(0);
  omp_set_num_threads(useCores);
#endif

  //Rprintf("\n multiplication");
  // matrix multiplication
  for (int k = 0; k < r; ++k) { // loops thru rows of WUt
    for (int j = 0; j < p; ++j) { // loops thru columns of std_X
      double sum = 0.0;
      // Rprintf("\n Multiplying the %d row of WUt ", k);
      // Rprintf("with the %d column of std_X", j);
      for (int i = 0; i < n; ++i) { // loops thru elements in the kth row of WUt/jth column of std_X
        sum += WUt[k * n + i] * std_X_acc[i][j];
      }
      rot_X_acc[k][j] = sum;
    }
  }

  // re-scale (analog to R function scale_varp())
  NumericVector scale_vals = colwise_l2mean(rot_X, r, p);
  scale_cols(rot_X, r, p, scale_vals);

  // save the means of the square values of each column (will pass to xtx in model fitting)
  NumericVector xtx = mean_sqsum(rot_X, r ,p);

  // cleanup and return
  Rprintf("\n return");
  Rcpp::List result;
  result["stdrot_X"] = rot_X; // the 'std' cues that this matrix has been standardized
  result["stdrot_X_scales"] = scale_vals;
  result["xtx"] = xtx;
  return result;
}
