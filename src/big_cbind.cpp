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


  // declarations: input
  NumericMatrix A(A_);
  int rows_A = A.nrow();
  int cols_A = A.ncol();
  XPtr<BigMatrix> Bptr(B_);
  MatrixAccessor<double> B(*Bptr);
  XPtr<BigMatrix> Cptr(C_);
  MatrixAccessor<double> C(*Cptr);
  int rows_B = Bptr->nrow();
  int cols_B = Bptr->ncol();
  int rows_C = Cptr->nrow();
  int cols_C = Cptr->ncol();


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


  // fill in first few columns of C with A
  for (int i=0;i<rows_A;i++){
    for (int j=0;j<cols_A;j++){
      C[j][i] = A[j*rows_A + i];
    }
  }

  // fill in the rest of the columns of C with B
  for (int i=0;i<rows_B;i++){
    for (int j=0;j<cols_B;j++){
      C[j + cols_A][i] = B[j][i];
    }
  }

  // return result
  Rcpp::List result;
  result["res"] = C;
  return result;
}
