#' Fit a linear mixed model with non-convex regularization
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param X Design matrix. May include clinical covariates and other non-SNP data. If this is the case, X_for_K should be supplied witha  matrix containing only SNP data for computation of GRM.
#' @param y Continuous outcome vector.
#' @param X_for_K X matrix used to compute the similarity matrix, K. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @param init
#' @param penalty
#' @param penalty.factor
#' @param gamma
#' @param alpha
#' @param lambda.min
#' @param nlambda
#' @param lambda
#' @param eps
#' @param max.iter
#' @param convex
#' @param dfmax
#' @param warn
#' @param standardize
#' @param rotation Logical to indicate whether the weighted rotation of the data should be performed (TRUE), or not (FALSE).
#' This is primarily for testing purposes and defaults to TRUE.
#' If \code{FALSE} results should correspond to those from \code{ncvreg} or \code{glmnet}.
#' @param returnX
#' @param
#' @importFrom zeallot %<-%
#' @export

### need to create a non-testing version of this function

### should i check for things to be missing or null?

plmm <- function(X,
                 y,
                 X_for_K = X,
                 penalty = c("MCP", "SCAD", "lasso"),
                 penalty.factor = rep(1, ncol(X)),
                 gamma,
                 alpha = 1,
                 lambda.min = ifelse(n>p, 0.001, 0.05),
                 nlambda = 100,
                 lambda,
                 init = rep(0, ncol(X)),
                 eps = 1e-04,
                 max.iter = 1000,
                 convex = TRUE,
                 dfmax = p + 1,
                 warn = TRUE,
                 standardize = TRUE,
                 rotation = FALSE, # just for testing - get rid of this argument
                 returnX = TRUE,
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
    c(S, U, eta) %<-% lmm_lasso_null(X_for_K, yy) # change name so it doesn't contain lasso
    W <- diag((eta * S + (1 - eta))^(-1/2))
  } else {
    U <- diag(n)
    W <- diag(n)
  }

  ## Intercept
  # XXX <- cbind(1, XX)
  # attributes(XXX)$center <- c(0, attr(XX, "center"))
  # attributes(XXX)$scale <- c(1, attr(XX, "scale"))
  # attributes(XXX)$nonsingular <- c(1, attr(XX, "nonsingular") + 1)
  # ns <- attr(XXX, "nonsingular")
  # penalty.factor <- penalty.factor[ns]

  ## Rotate data
  SUX <- W %*% crossprod(U, cbind(1, XX))
  SUy <- drop(W %*% crossprod(U, yy))

  ## Set up lambda
  if (missing(lambda)) {
    # Don't include intercept in calculating lambda
    lambda <- ncvreg::setupLambda(SUX[,-1], SUy - mean(SUy), "gaussian", alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Placeholders for results
  init <- c(mean(SUy), init) # add initial value for intercept
  b <- matrix(NA, ncol(SUX), nlambda)
  iter <- integer(nlambda)
  loss <- numeric(nlambda)
  for (ll in 1:nlambda){
    lam <- lambda[ll]
    res <- ncvreg::ncvfit(SUX, SUy, init, penalty, gamma, alpha, lam, eps, max.iter, c(0, penalty.factor), warn)
    b[, ll] <- init <- res$beta
    iter[ll] <- res$iter
    loss[ll] <- res$loss
  }

  ### from ncvreg ##############################################################

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  a <- b[1, ind]
  b <- b[-1, ind, drop=FALSE]
  iter <- iter[ind]
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
                        iter = iter,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        convex.min = convex.min,
                        loss = loss,
                        penalty.factor = penalty.factor,
                        n = n),
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

### Need to add some checks for non-penalized covars!

