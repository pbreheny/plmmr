#' Print method for summary.cv_plmm objects
#'
#' @param x An object of class \code{summary.cv_plmm}
#' @param digits The number of digits to use in formatting output
#' @param ... Not used
#'
#' @rdname print.summary.cv_plmm
#'
#' @returns Nothing is returned; instead, a message is printed to the console
#' summarizing the results of the cross-validated model fit.
#'
#'
#' @export
#'
#' @examples
#' admix_design <- create_design(X = admix$X, y = admix$y)
#' cv_fit <- cv_plmm(design = admix_design)
#' print(summary(cv_fit))
#'
print.summary.cv_plmm <- function(x, digits, ...){
  n <- nrow(x$fit$std_Xbeta)
  p <- nrow(x$fit$beta_vals)
  digits <- if (missing(digits)) digits <- c(2, 4, 2, 2, 3) else rep(digits, length.out=5)
  cat(x$fit$penalty, "-penalized model with n=", n, " and p=", p,"\n", sep="")
  cat("At minimum cross-validation error (lambda=", formatC(x$lambda_min, digits[2], format="f"), "):\n", sep="")
  cat("-------------------------------------------------\n")
  cat("  Nonzero coefficients: ", x$nvars, "\n", sep="")
  cat("  Cross-validation error (deviance): ", formatC(min(x$cve), digits[1], format="f"), "\n", sep="")
  cat("  Scale estimate (sigma): ", formatC(sqrt(x$cve[x$min]), digits[5], format="f"), "\n", sep="")
}
