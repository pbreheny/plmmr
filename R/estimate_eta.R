#' Estimate eta (to be used in rotating the data)
#'
#' This function is called internally by `plmm()`
#'
#' @param n The number of observations
#' @param s The non-zero eigenvalues of K, the realized relationship matrix
#' @param U The eigenvectors of K associated with s
#' @param y Continuous outcome vector
#' @param incpt_flag Logical: Does the model require fitting an intercept?
#'
#' @return a numeric value with the estimated value of eta, the variance parameter
#'
#' @keywords internal
#'
estimate_eta <- function(n, s, U, y, incpt_flag) {
  opt <- stats::optimize(
    f = log_lik,
    c(0.01, 0.99),
    n = n,
    s = s,
    U = U,
    y = y,
    incpt_flag = incpt_flag
  )

  opt$minimum
}
