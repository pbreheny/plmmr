#' A function to print the summary of a \code{plmm} model
#'
#' @param x A `summary.plmm` object
#' @param ... Not used
#' @rdname print.summary.plmm
#'
#' @returns Nothing is returned; instead, a message is printed to the console
#' summarizing the results of the model fit.
#'
#'
#' @export
#'
#' @examples
#' lam <- rev(seq(0.01, 1, length.out=20)) |> round(2) # for sake of example
#' fit <- plmm(X = admix$X, y = admix$y, lambda = lam)
#' fit2 <- plmm(X = admix$X, y = admix$y, penalty = "SCAD", lambda = lam)
#' print(summary(fit, idx = 12))
#' print(summary(fit2, lambda = 0.11))
print.summary.plmm <- function(x, ...){
  cat(x$penalty, "-penalized regression model with n=", x$std_X_n, ", p=", x$p, sep="")
  cat(" at lambda=", x$lambda_char, "\n", sep="")
  cat("-------------------------------------------------\n")

  # did the model converge?
  if(x$converged){
      cat("The model converged", "\n")
    } else {
      cat("The model did not converge - max. number of iterations reached", "\n")
    }

  cat("-------------------------------------------------\n")
  # nonzero coefficients
  cat("# of non-zero coefficients: ", x$nvars, "\n")
  cat("-------------------------------------------------\n")


}
