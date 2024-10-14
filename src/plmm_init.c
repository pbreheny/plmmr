#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
// #include <R_ext/Applic.h>
#include <R_ext/Visibility.h>  // optional


// find column-wise standard deviation of a filebacked matrix
extern SEXP big_sd(SEXP X_, SEXP ncore_);

// standardize a matrix of filebacked data
extern SEXP big_std(SEXP X_, SEXP ncore_);

// standardize a matrix of in-memory data
extern SEXP in_mem_std(SEXP X_);

// big_crossprod
extern SEXP big_crossprod(SEXP X_,
                          SEXP y_,
                          SEXP ncore_);

static const R_CallMethodDef callMethods[] = {
  {"big_crossprod", (DL_FUNC) &big_crossprod, 4},
  {"big_std", (DL_FUNC) &big_std, 2},
  {"big_sd", (DL_FUNC) &big_sd, 2},
  {"in_mem_std", (DL_FUNC) &in_mem_std, 1},
  {NULL, NULL, 0}
};

void R_init_plmm(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
