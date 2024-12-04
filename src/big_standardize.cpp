#include "utilities.h"

RcppExport SEXP big_std(SEXP X_,
                        SEXP ncore_,
                        SEXP center_ = R_NilValue,
                        SEXP scale_ = R_NilValue){


  // *********************** NOTE ***********************
  // This function modifies the X that is currently being
  // used in the corresponding R session.
  //*****************************************************


  // declarations
  // Rprintf("\nDeclarations in big_std()");
  XPtr<BigMatrix> X(X_); // points to the filebacked matrix of rotated data
  // MatrixAccessor<double> X_acc(*X);
  int n = X->nrow();
  int p = X->ncol();


  // set up omp
  // Rprintf("\n set up OpenMP");
  int useCores = INTEGER(ncore_)[0];
#ifdef PLMMR_OMP_H_
  int haveCores = omp_get_num_procs();
  if(useCores < 1) {
    useCores = haveCores;
  }
  omp_set_dynamic(0);
  omp_set_num_threads(useCores);
#endif

  // re-center
  // Rprintf("\nCalling col_means()");
  NumericVector center_vals;
  if (Rf_isNull(center_)) {
    center_vals = col_means(X, n, p);
  } else {
    center_vals = as<NumericVector>(center_);
  }
  center_cols(X, n, p, center_vals);

  // re-scale
  // Rprintf("\nCalling scaling functions");
  NumericVector scale_vals;
  if (Rf_isNull(scale_)) {
    scale_vals = colwise_l2mean(X, n, p);
  } else {
    scale_vals = as<NumericVector>(scale_);
  }
  scale_cols(X, n, p, scale_vals);

  // save the means of the square values of each column (will pass to xtx in model fitting)
  //NumericVector xtx = mean_sqsum(X, r ,p);

  Rcpp::List result;
  result["std_X"] = X; // the 'std' cues that this matrix has been standardized
  result["std_X_center"] = center_vals;
  result["std_X_scale"] = scale_vals;
  //result["xtx"] = xtx;
  return result;
}
