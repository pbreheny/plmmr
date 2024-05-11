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
#' cv_fit <- cv_plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' print(summary(cv_fit))
#' 
print.summary.cv_plmm <- function(x, digits, ...){
  digits <- if (missing(digits)) digits <- c(2, 4, 2, 2, 3) else rep(digits, length.out=5)
  cat(x$fit$penalty, "-penalized model with n=", x$fit$std_X_n, " and p=", x$fit$p,"\n", sep="")
  cat("At minimum cross-validation error (lambda=", formatC(x$lambda.min, digits[2], format="f"), "):\n", sep="")
  cat("-------------------------------------------------\n")
  cat("  Nonzero coefficients: ", x$nvars, "\n", sep="")
  cat("  Cross-validation error (deviance): ", formatC(min(x$cve), digits[1], format="f"), "\n", sep="")
  cat("  Scale estimate (sigma): ", formatC(sqrt(x$cve[x$min]), digits[5], format="f"), "\n", sep="")
}
