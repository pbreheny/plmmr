#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
//#include <math.h>
#include <string.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>

// rotate_filebacked
extern SEXP rotate_filebacked(SEXP std_X_, SEXP WUt_, SEXP y_, SEXP ncore_,
                              const char* backingfile);


// Cross product of y with jth column of X
double crossprod(double *X, double *y, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i];
  return(val);
}

// TODO: adjust what is below to handle filebacked case
// Sum of squares of jth column of X
double sqsum(double *X, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += pow(X[nn+i], 2);
  return(val);
}

// Gaussian loss
double g_loss(double *r, int n) {
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
}

static const R_CallMethodDef CallEntries[] = {
  {"rotate_filebacked", (DL_FUNC) &rotate_filebacked, 5},
  {NULL, NULL, 0}
};

void R_init_plmm(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
