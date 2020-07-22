#' Fit a linear mixed model with non-convex regularization
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param X Design matrix of SNP data.
#' @param y Continuous outcome vector.
#' @param X0 Design matrix of non-SNP data.
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
                 X0,
                 X_for_K,
                 init,
                 penalty = c("MCP", "SCAD", "lasso"),
                 penalty.factor,
                 gamma,
                 alpha = 1,
                 lambda.min = ifelse(n>p, 0.001, 0.05),
                 nlambda = 100,
                 lambda,
                 eps = 1e-04,
                 max.iter = 1000,
                 convex = FALSE,
                 dfmax = p + 1,
                 warn = TRUE,
                 standardize = TRUE,
                 rotation = FALSE, # just for testing - get rid of this argument
                 returnX = TRUE,
                 ...) {



  # Coersion
  S <- U <- eta <- NULL
  penalty <- match.arg(penalty)
  if (missing(X0)){
    p0 <- 1 # intercept
  } else {
    p0 <- 1 + ncol(X0) # intercept + covars
  }
  if (missing(penalty.factor)){ # set up penalty factor so only SNPs are penalized
    penalty.factor <- rep(c(0, 1), times = c(p0, ncol(X)))
  } else {
    if (length(penalty.factor) != p0 + ncol(X)){
      penalty.factor <- rep(c(0, 1), times = c(p0, ncol(X)))
      warning(paste("penalty.factor must have length equal to supplied covariates + 1 (for intercept). Setting `penalty.factor <- rep(c(0, 1), times = c(", p0, ", ncol(X)))`."))
    }
  }
  if (missing(X_for_K)) X_for_K <- X # define SNP matrix used to construct RRM
  if (!missing(X0)) X <- cbind(X0, X) # set up full X matrix - SNPs + other covars
  if (missing(init)) init <- rep(0, 1 + ncol(X))
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
  # if (length(penalty.factor)!=ncol(X)) stop("Dimensions of penalty.factor and X do not match", call.=FALSE)
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
  # ns <- attr(XX, "nonsingular")
  # penalty.factor <- penalty.factor[ns]
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
  XXX <- cbind(1, XX)
  attributes(XXX)$center <- c(0, attr(XX, "center"))
  attributes(XXX)$scale <- c(1, attr(XX, "scale"))
  attributes(XXX)$nonsingular <- c(1, attr(XX, "nonsingular") + 1)
  ns <- attr(XXX, "nonsingular")
  penalty.factor <- penalty.factor[ns]

  ## Rotate data
  SUX <- W %*% crossprod(U, XXX)
  SUy <- drop(W %*% crossprod(U, yy))

  ## Set up lambda
  if (missing(lambda)) {
    # Don't include intercept in calculating lambda (?)
    lambda <- ncvreg::setupLambda(SUX[,-1], SUy, "gaussian", alpha, lambda.min, nlambda, penalty.factor[-1])
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Placeholders for results
  b <- matrix(NA, ncol(SUX), nlambda)
  iter <- integer(nlambda)
  loss <- numeric(nlambda)
  for (ll in 1:nlambda){
    lam <- lambda[ll]
    res <- ncvreg::ncvfit(SUX, SUy, init, penalty, gamma, alpha, lam, eps, max.iter, penalty.factor, warn)
    b[, ll] <- init <- res$beta
    iter[ll] <- res$iter
    loss[ll] <- res$loss
  }

  ### from ncvreg ##############################################################

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & sum(iter) == max.iter) warning("Maximum number of iterations reached")

  ## Local convexity? - later?
  # convex.min <- if (convex) convexMin(b, XX, penalty, gamma, lambda*(1-alpha), family, penalty.factor, a=a) else NULL

  ## Unstandardize
  beta <- matrix(0, nrow=(ncol(SUX)), ncol=length(lambda))
  bb <- b/attr(XXX, "scale")[ns]
  beta[ns[-1],] <- bb[-1,]
  a <- bb[1,]
  beta[1,] <- a - crossprod(attr(XXX, "center")[ns[-1]], bb[-1,])
  # if (intercept){ # manual intercept
  #   beta <- matrix(0, nrow=(ncol(XX)), ncol=length(lambda))
  #   bb <- b/attr(XX, "scale")[ns]
  #   beta[ns[-1],] <- bb[-1,]
  #   a <- bb[1,]
  #   beta[1,] <- a - crossprod(attr(XX, "center")[ns[-1]], bb[-1,])
  # } else { # no manual standardization, no manual intercept
  #   beta <- matrix(0, nrow=(ncol(XX) + 1), ncol=length(lambda))
  #   bb <- b/attr(XX, "scale")[ns]
  #   beta[ns+1,] <- bb
  #   beta[1,] <- a - crossprod(attr(XX, "center")[ns], bb)
  # }

  ## Names
  lamNames <- function(l) {
    if (length(l) > 1) {
      d <- ceiling(-log10(-max(diff(l))))
      d <- min(max(d,4), 10)
    } else {
      d <- 4
    }
    formatC(l, format="f", digits=d)
  }
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
                        # convex.min = convex.min,
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

