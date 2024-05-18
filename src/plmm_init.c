#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
// #include <R_ext/Applic.h>
#include <R_ext/Visibility.h>  // optional

// re-scale a rotated matrix of data
extern SEXP big_rescale(SEXP rot_X_,
                       SEXP ncore_);

// rotate_filebacked
extern SEXP rotate_filebacked(SEXP std_X_,
                              SEXP WUt_,
                              SEXP rot_X_,
                              SEXP ncore_);

// big_crossprod
extern SEXP big_crossprod(SEXP X_,
                   SEXP y_,
                   SEXP ind_col_,
                   SEXP ncore_);

static const R_CallMethodDef callMethods[] = {
  {"rotate_filebacked", (DL_FUNC) &rotate_filebacked, 4},
  {"big_crossprod", (DL_FUNC) &big_crossprod, 4},
  {"big_rescale", (DL_FUNC) &big_rescale, 2},
  {NULL, NULL, 0}
};

void R_init_plmm(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}