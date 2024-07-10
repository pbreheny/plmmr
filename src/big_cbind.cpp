#include "utilities.h"

// a version of cbind() for file-backed matrices
RcppExport SEXP big_cbind(SEXP A_, // in-memory data
                          SEXP B_, // file-backed data
                          SEXP C_, // file-backed placeholder for combined data
                          SEXP ncore_){



  // *********************** NOTE ***********************
  // This function modifies the C placeholder matrix
  // that is currently being used
  // in the corresponding R session.
  //*****************************************************

  Rprintf("\nDeclarations in big_cbind()");
  // declarations: input
  NumericMatrix A(A_);
  int rows_A = A.nrow();
  int cols_A = A.ncol();
  XPtr<BigMatrix> B(B_);
  XPtr<BigMatrix> C(C_);
  int rows_B = B->nrow();
  int cols_B = B->ncol();
  int rows_C = C->nrow();
  int cols_C = C->ncol();


  // setup OMP
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

  Rprintf("\nfill in first few columns of C with A");
  // fill in first few columns of C with A
  fill_in(C, rows_C, cols_C, A, rows_A, cols_A);

  Rprintf("\nfill in the rest of the columns of C with B");
  // fill in the rest of the columns of C with B
  fill_in_filebacked(C, rows_C, cols_C, B, rows_B, cols_B, cols_A);

  // return result
  Rcpp::List result;
  result["res"] = C;
  return result;
}
