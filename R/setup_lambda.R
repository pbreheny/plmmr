#' Compute sequence of lambda values
#'
#' This function allows you compute a sequence of lambda values for plmm models.
#' @param X Design matrix. May include clinical covariates and other non-SNP data.
#' @param y Continuous outcome vector.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param standardize Logical flag for X variable standardization, prior to data transformation. The coefficients are always returned on the original scale. Default is TRUE. If variables are in the same units already, or manually standardized, you might not wish to standardize.
#' @param rotation Logical flag to indicate whether the weighted rotation of the data should be performed (TRUE), or not (FALSE). This is primarily for testing purposes and defaults to TRUE.
#' @export


setup_lambda <- function(X, y, alpha, lambda.min, nlambda, penalty.factor, standardize, rotation) {



  if (standardize & !rotation){

    lambda <- ncvreg::setupLambda(X, y, 'gaussian', alpha, lambda.min, nlambda, penalty.factor)

  } else { # rotated X not standardized so need to actually compute cov with residuals

    n <- nrow(X)
    p <- ncol(X)

    ind <- which(penalty.factor!=0)
    if (length(ind)!=p) {
      fit <- stats::glm(y~X[, -ind], family='gaussian')
    } else {
      fit <- stats::glm(y~1, family='gaussian')
    }
    # double check modifications for unpenalized things
    zmax <- max(stats::cov(X[, ind], fit$residuals) / penalty.factor[ind] * (n-1) / n)
    lambda.max <- zmax/alpha

    if (lambda.min==0) {
      lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)), 0)
    } else {
      lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len=nlambda))
    }

    if (length(ind)!=p) lambda[1] <- lambda[1] * 1.000001
    lambda
  }

}
