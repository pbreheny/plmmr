#' Fit a linear mixed model with non-convex regularization
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param X Design matrix. May include clinical covariates and other non-SNP data. If this is the case, X_for_K should be supplied witha  matrix containing only SNP data for computation of GRM.
#' @param y Continuous outcome vector.
#' @param X_for_K X matrix used to compute the similarity matrix, K. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param lambda A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.
#' @param eps Convergence threshhold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max.iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param convex Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param init Initial values for coefficients.
#' @param warn Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param returnX Return the standardized design matrix along with the fit? By default, this option is turned on if X is under 100 MB, but turned off for larger matrices to preserve memory.
#' @param standardize Logical flag for X variable standardization, prior to data transformation. The coefficients are always returned on the original scale. Default is TRUE. If variables are in the same units already, or manually standardized, you might not wish to standardize.
#' @param rotation Logical flag to indicate whether the weighted rotation of the data should be performed (TRUE), or not (FALSE). This is primarily for testing purposes and defaults to TRUE.
#' @param ... Not used.
#' @importFrom zeallot %<-%
#' @export

### need to create a non-testing version of this function

### should i check for things to be missing or null?

plmm <- function(X,
                 y,
                 X_for_K = X,
                 penalty = c("MCP", "SCAD", "lasso"),
                 gamma,
                 alpha = 1,
                 lambda.min = ifelse(n>p, 0.001, 0.05),
                 nlambda = 100,
                 lambda,
                 eps = 1e-04,
                 max.iter = 10000,
                 convex = TRUE,
                 dfmax = p + 1,
                 warn = TRUE,
                 penalty.factor = rep(1, ncol(X)),
                 init = rep(0, ncol(X)),
                 returnX = TRUE,
                 standardize = TRUE,
                 rotation = TRUE,
                 ...) {



  # Coersion
  S <- U <- eta <- NULL
  penalty <- match.arg(penalty)
  if (missing(gamma)) gamma <- switch(penalty, SCAD = 3.7, 3)
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (typeof(X)=="integer") storage.mode(X) <- "double"
  if (typeof(X)=="character") stop("X must be a numeric matrix", call.=FALSE)
  if (!is.double(y)) {
    op <- options(warn=2)
    on.exit(options(op))
    y <- tryCatch(
      error = function(cond) stop("y must be numeric or able to be coerced to numeric", call.=FALSE),
      as.double(y))
    options(op)
  }
  if (!is.double(penalty.factor)) penalty.factor <- as.double(penalty.factor)

  # Error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (length(penalty.factor)!=ncol(X)) stop("Dimensions of penalty.factor and X do not match", call.=FALSE)
  if (length(init)!=ncol(X)) stop("Dimensions of init and X do not match", call.=FALSE)
  if (length(y) != nrow(X)) stop("X and y do not have the same number of observations", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg", call.=FALSE)

  ## Deprecation support
  dots <- list(...)
  if ("n.lambda" %in% names(dots)) nlambda <- dots$n.lambda

  ## Set up XX, yy, lambda
  if (standardize){
    XX <- ncvreg::std(X)
    X_for_K <- ncvreg::std(X_for_K)
  } else {
    XX <- X
    attributes(XX)$center <- rep(0, ncol(XX))
    attributes(XX)$scale <- rep(1, ncol(XX))
    attributes(XX)$nonsingular <- 1:ncol(XX)
  }
  ns <- attr(XX, "nonsingular")
  penalty.factor <- penalty.factor[ns]
  yy <- y
  p <- ncol(XX)
  n <- length(yy)
  ##############################################################################
  ### get things working w/o rotation first

  ## Calculate eta
  if (rotation){
    c(S, U, eta) %<-% plmm_null(X_for_K, yy) # change name so it doesn't contain lasso
    W <- diag((eta * S + (1 - eta))^(-1/2))
  } else {
    U <- diag(n)
    W <- diag(n)
  }

  ## Rotate data
  SUX <- W %*% crossprod(U, cbind(1, XX)) # manual intercept
  SUy <- drop(W %*% crossprod(U, yy))

  ## Set up lambda
  if (missing(lambda)) {
    # Don't include the intercept in calculating lambda
    lambda <- ncvreg::setupLambda(SUX[,-1], SUy - mean(SUy), "gaussian", alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Placeholders for results
  init <- c(mean(yy), init) # add initial value for intercept
  b <- matrix(NA, ncol(SUX), nlambda)
  iter <- integer(nlambda)
  converged <- logical(nlambda)
  loss <- numeric(nlambda)
  # think about putting this loop in C
  for (ll in 1:nlambda){
    lam <- lambda[ll]
    res <- ncvreg::ncvfit(SUX, SUy, init, penalty, gamma, alpha, lam, eps, max.iter, c(0, penalty.factor), warn)
    b[, ll] <- init <- res$beta
    iter[ll] <- res$iter
    converged[ll] <- ifelse(res$iter < max.iter, TRUE, FALSE)
    loss[ll] <- res$loss
    # add break for saturation?
  }

  ### from ncvreg ##############################################################

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  a <- b[1, ind]
  b <- b[-1, ind, drop=FALSE]
  iter <- iter[ind]
  converged <- converged[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & sum(iter) == max.iter) warning("Maximum number of iterations reached")

  ## Local convexity?
  convex.min <- if (convex) convexMin(b, SUX[,-1], penalty, gamma, lambda*(1-alpha), family = 'gaussian', penalty.factor, a=a) else NULL

  ## Unstandardize
  beta <- matrix(0, nrow=(ncol(SUX)), ncol=length(lambda))
  bb <- b/attr(XX, "scale")[ns]
  beta[ns + 1,] <- bb
  beta[1,] <- a - crossprod(attr(XX, "center")[ns], bb)

  varnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)
  varnames <- c("(Intercept)", varnames)
  dimnames(beta) <- list(varnames, lamNames(lambda))

  ## Output
  val <- structure(list(beta = beta,
                        eta = eta,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        convex.min = convex.min,
                        loss = loss,
                        penalty.factor = penalty.factor,
                        n = n,
                        iter = iter,
                        converged = converged),
                        class = "plmm")
  if (missing(returnX)) {
    if (utils::object.size(SUX) > 1e8) {
      warning("Due to the large size of SUX (>100 Mb), returnX has been turned off.\nTo turn this message off, explicitly specify returnX=TRUE or returnX=FALSE).")
      returnX <- FALSE
    } else {
      returnX <- TRUE
    }
  }
  if (returnX) {
    val$X <- SUX
    val$y <- SUy
  }
  return(val)
}

