#' Compute sequence of lambda values
#'
#' This function allows you compute a sequence of lambda values for plmm models.
#' @param X Rotated and standardized design matrix which *includes* the intercept column if present. May include clinical covariates and other non-SNP data. This can be either a 'matrix' or 'FBM' object.
#' @param y Continuous outcome vector.
#' @param intercept Logical: does X contain an intercept column? Defaults to TRUE.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda_min The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise. A value of lambda_min = 0 is not supported.
#' @param nlambda The desired number of lambda values in the sequence to be generated.
#' @param penalty_factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty_factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty_factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty_factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @keywords internal
#'
#' @returns a numeric vector of lambda values, equally spaced on the log scale
#'
setup_lambda <- function(X, y, alpha, lambda_min, nlambda, penalty_factor, intercept = TRUE) {

  # make sure alpha is neither missing nor zero
  if (is.na(alpha) | is.null(alpha) | alpha == 0) {
    stop("Must provide a non-zero value for alpha.")
  }
  # make sure user is not trying to use a lambda_min value of 0
  if (!missing(lambda_min)){
    if(lambda_min == 0){stop("User-specified value for lambda_min cannot be zero")}
  }

  # label dimensions of X
  n <- nrow(X)
  p <- ncol(X)

  # fit the model to unpenalized covariates (residuals will be used later)
  ind <- which(penalty_factor != 0)
  if (length(ind)!=p) {
    fit <- lm(y ~ X[, -ind])
  } else {
    fit <- lm(y~1)
  }

  # determine the maximum value for lambda
  if('matrix' %in% class(X)){
    decomp_backsolve <- abs(crossprod(X[,ind], fit$residuals))/penalty_factor[ind]
  } else {
    # NOTE: the following line is the reason that the penalized columns must be a
    # *contiguous* submatrix! We checked for this back in plmm_checks()
    pen_X <- bigmemory::sub.big.matrix(
      X,
      firstCol = ind[1],
      lastCol = ind[length(ind)])
    cprod_res <- .Call("big_crossprod",
                       pen_X@address,
                       fit$residuals,
                       as.integer(count_cores()),
                       PACKAGE = "plmmr")
    cprod <- cprod_res[[1]]
    decomp_backsolve <- abs(cprod)/penalty_factor[ind]
  }
  zmax <- max(stats::na.exclude(decomp_backsolve))/n
  lambda.max <- zmax/alpha

  # error check
  if(!is.finite(log(lambda.max))){stop("log(lambda.max) is not finite")}

  # Default is .001 if the number of observations is larger than the number of
  # covariates and .05 otherwise. A value of lambda_min = 0 is not supported.
  if(missing(lambda_min)){ # case 1: if the user does not specify a lambda_min value...
    if(n > p){
      lambda <- exp(seq(log(lambda.max), log(0.001*lambda.max), len = nlambda))
    } else {
      lambda <- exp(seq(log(lambda.max), log(0.05*lambda.max), len = nlambda))
    }
  } else { # case 2: the user specifies a (nonzero) lambda_min value
    lambda <- exp(seq(log(lambda.max), log(lambda_min*lambda.max), len = nlambda))
  }

}
