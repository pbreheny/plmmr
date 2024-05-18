#include "utilities.h"

RcppExport SEXP big_rescale(SEXP rot_X_,
                            SEXP ncore_){

  // declarations
  XPtr<BigMatrix> rot_X(rot_X_); // points to the filebacked matrix of rotated data
  MatrixAccessor<double> rot_X_acc(*rot_X);
  int r = rot_X->nrow();
  int p = rot_X->ncol();


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

  // re-scale (analog to R function scale_varp())
  NumericVector scale_vals = colwise_l2mean(rot_X, r, p);
  rescale_cols(rot_X, r, p, scale_vals);

  // save the means of the square values of each column (will pass to xtx in model fitting)
  NumericVector xtx = mean_sqsum(rot_X, r ,p);

  Rcpp::List result;
  result["stdrot_X"] = rot_X; // the 'std' cues that this matrix has been standardized
  result["stdrot_X_scales"] = scale_vals;
  result["xtx"] = xtx;
  return result;
}