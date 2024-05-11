#include <string.h>
//#include <math.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
//#include <Rmath.h>
#include <R_ext/Applic.h>
#include <RcppArmadillo.h>
#include "bigmemory/BigMatrix.h"
#include <time.h>
#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/bigmemoryDefines.h"
#include "plmmr_omp.h"

// rotating filebacked data 
SEXP rotate_filebacked(SEXP std_X_, SEXP WUt_, SEXP y_, SEXP ncore_, const char* backingfile) {
  
  // declarations: input
  XPtr<BigMatrix> xMat(std_X_); // this points to the filebacked matrix of standardized data
  double *WUt = REAL(WUt_); // this points to the in-memory n x r projection matrix
  double *y = REAL(y_); // this is the observed outcome (an n x 1 vector)
  int n = xMat->nrow(); // n = number of observations 
  int p = xMat->ncol(); // p = number of features (including non-genomic predictors!)
  int r = ncols(WUt_); // r = number of eigenvalues used in rotation (r stands for 'rank')
  MatrixAccessor<double> std_X_acc(*std_X_); // initialize MatrixAccessor with std_X
  
  // // declarations: output
  // SEXP rot_X;
  // PROTECT(rot_X = R_NilValue);
  
  // create the filebacked.big.matrix
  FilebackedBigMatrix* rot_X = new FilebackedBigMatrix(r, p, backingfile);
  
  // convert the C pointer to an R external pointer
  // rot_X = R_MakeExternalPtr(fbm, R_NilValue, R_NilValue);
  
  // set class attribute to indicate it's a filebacked.big.matrix
  // Rf_classgets(rot_X, Rf_mkString("filebacked.big.matrix"));
  
  // Initialize MatrixAccessor with rot_X
  MatrixAccessor<double> rot_X_acc(*rot_X);
  
  // set up omp
  int useCores = INTEGER(ncore_)[0];
#ifdef PLMMR_OMP_H_
  int haveCores = omp_get_num_procs();
  if(useCores < 1) {
    useCores = haveCores;
  }
  omp_set_dynamic(0);
  omp_set_num_threads(useCores);
#endif
  
  // matrix multiplication 
  for (int k = 0; k < r; ++k) { // loops thru rows of WUt
    for (int j = 0; j < p; ++j) { // loops thru columns of std_X
      double sum = 0.0;
      for (int i = 0; i < n; ++i) { // loops thru elements in the kth row of WUt/jth column of std_X
        sum += WUt[k * n + i] * std_X_acc[i][j];
      }
      rot_X_acc[k][j] = sum;
    }
  }
  
  // cleanup and return
  delete rot_X; // release memory for the filebacked.big.matrix
  return rot_X;
}