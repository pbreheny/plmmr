#' A summary method for the plmm objects
#'
#' @param object An object of class \code{plmm}
#' @param lambda The regularization parameter value at which inference should be reported.
#' @param idx Alternatively, `lambda` may be specified by an index; `idx=10` means:
#'      report inference for the 10th value of `lambda` along the regularization path. If both `lambda` and `idx` are specified, `lambda` takes precedence.
#' @param eps If lambda is given, eps is the tolerance for difference between the given lambda value and a lambda value from the object. Defaults to 0.0001 (1e-5)
#' @param ... Not used
#'
#' @return The return value is an object with S3 class `summary.plmm`. The class has its own print method and contains the following list elements:
#' * `penalty`: The penalty used by `plmm` (e.g. SCAD, MCP, lasso)
#' * `n`: Number of instances/observations
#' * `std_X_n`: the number of observations in the standardized data;
#'       the only time this would differ from 'n' is if data are from PLINK and the external data does not include all the same samples
#' * `p`: Number of regression coefficients (not including the intercept)
#' * `converged`: Logical indicator for whether the model converged
#' * `lambda`: The `lambda` value at which inference is being reported
#' * `lambda_char`: A formatted character string indicating the lambda value
#' * `nvars`: The number of nonzero coefficients (again, not including the intercept) at that value of `lambda`
#' * `nonzero`: The column names indicating the nonzero coefficients in the model at the specified value of `lambda`
#'
#' @rdname summary.plmm
#' @export
#'
#' @examples
#' admix_design <- create_design(X = admix$X, y = admix$y)
#' fit <- plmm(design = admix_design)
#' summary(fit, idx = 97)
summary.plmm <- function(object, lambda, idx, eps = 1e-5, ...) {

  # lambda/which
  if (missing(lambda) && missing(idx)) stop("One of the arguments 'lambda' or 'idx' must be provided.")
  if (missing(lambda)) lambda <- object$lambda[idx]
  if (missing(idx)) {
    idx <- which(abs(object$lambda - lambda) < eps)
    if (sum(idx) == 0) {
      warning("The user-specified lambda value is not within epsilon of any lambda values used to fit the model. Will proceed with the fitted model lambda value closest to the user-specified lambda.")
      idx <- which.min(abs(object$lambda - lambda))
    }

  }

  # nvars (tells number of non-zero coefficients)
  nvars <- predict(object, type = "nvars", lambda = lambda, idx = idx, ...)

  # error checking
  if (length(nvars) > 1) stop("You must specify a single model (i.e., a single value of lambda)", call. = FALSE)

  lambda_char <- colnames(object$beta_vals)[idx]


  # tells WHICH variables have non-zero coefficients
  nz <- which(object$beta_vals[, lambda_char] > .Machine$double.eps)
  nonzero <- rownames(object$beta_vals[nz, lambda_char, drop = FALSE])
  # don't drop, because we need the dimnames here ^


  out <- structure(list(
    penalty = object$penalty,
    n = nrow(object$K$U),
    std_X_n = object$std_X_n,
    p = nrow(object$beta_vals),
    converged = object$converged[idx],
    lambda = lambda,
    lambda_char = lambda_char,
    nvars = nvars,
    nonzero = nonzero
  ),
  class = "summary.plmm")

  return(out)

}
