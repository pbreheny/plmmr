#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
// #include <R_ext/Applic.h>
#include <R_ext/Visibility.h>  // optional

// rotate_filebacked
extern SEXP rotate_filebacked(SEXP std_X_,
                              SEXP WUt_, 
                              SEXP rot_X_, 
                              SEXP ncore_);




static const R_CallMethodDef callMethods[] = {
  {"rotate_filebacked", (DL_FUNC) &rotate_filebacked, 4},
  {NULL, NULL, 0}
};

void R_init_plmm(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
