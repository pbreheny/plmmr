#include "utilities.h"

// cross product of filebacked matrix X with vector y
RcppExport SEXP big_crossprod(SEXP X_,
                              SEXP y_,
                              SEXP ncore_){

  // declarations: input
  XPtr<BigMatrix> X(X_); // this points to the filebacked matrix of standardized data
  int n = X->nrow(); // n = number of observations
  int p = X->ncol(); // p = number of features
  double *y = REAL(y_);
  NumericVector res(p);

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

  // multiplication
  for (int j=0;j<p;j++){ // looping over rows of t(X) [columns of X]
    // Rprintf("multiplication with column %d\n", j);
    res[j] = crossprod(X, y, j, n);
    // Rprintf("crossprod's res is %f\n", res[j]);
  }
  Rcpp::List result;
  result["cp"] = res; // cp = 'cross product'
  return result;
}