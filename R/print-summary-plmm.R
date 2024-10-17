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
#' admix_design <- create_design(X = admix$X, y = admix$y)
#' fit <- plmm(design = admix_design, lambda = lam)
#' fit2 <- plmm(design = admix_design, penalty = "SCAD", lambda = lam)
#' print(summary(fit, idx = 18))
#' print(summary(fit2, idx = 18))
print.summary.plmm <- function(x, ...){
  cat(x$penalty, "-penalized regression model with n=", x$n, ", p=", x$p, sep="")
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
