#include "utilities.h"
// Note: this function is taken from the R package ncvreg;
// See https://github.com/pbreheny/ncvreg/blob/master/src/standardize.c

RcppExport SEXP in_mem_std(SEXP X_) {

  // Declarations
  NumericMatrix X = NumericMatrix(X_);
  int n = X.nrow();
  int p = X.ncol();
  NumericMatrix XX(n, p);
  NumericVector c(p);
  NumericVector s(p);

  for (int j = 0; j < p; j++) {
    // Center
    c[j] = 0;
    for (int i = 0; i < n; i++) {
      c[j] += X(i, j);
    }
    c[j] = c[j] / n;
    for (int i = 0; i < n; i++) {
      XX(i, j) = X(i, j) - c[j];
    }

    // Scale
    s[j] = 0;
    for (int i = 0; i < n; i++) {
      s[j] += pow(XX(i, j), 2);
    }
    s[j] = sqrt(s[j] / n);
    for (int i = 0; i < n; i++) {
      XX(i, j) = XX(i, j) / s[j];
    }
  }

  // Return list
  return List::create(Named("XX") = XX, Named("c") = c, Named("s") = s);
}
