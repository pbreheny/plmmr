#' Fit a linear mixed model with non-convex regularization
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param X_for_K X matrix used to compute the similarity matrix, K. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @importFrom zeallot %<-%
#' @export



penalizedLMM <- function(X,
                       y,
                       X_for_K = NULL,
                       init = rep(0, ncol(X)),
                       penalty = c("MCP", "SCAD", "lasso"),
                       gamma = switch(penalty, SCAD = 3.7, 3),
                       alpha = 1,
                       lambda.min = ifelse(n>p,.001,.05),
                       nlambda = 100,
                       lambda,
                       eps = 1e-04,
                       max.iter = 1000,
                       convex = FALSE,
                       dfmax = p + 1,
                       penalty.factor = rep(1, ncol(X)),
                       warn = TRUE) {

  S <- U <- eta <- NULL

  ### from ncvreg ##############################################################
  # Coersion
  penalty <- match.arg(penalty)
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
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
  if (length(y) != nrow(X)) stop("X and y do not have the same number of observations", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg", call.=FALSE)

  ## Deprecation support
  dots <- list(...)
  if ("n.lambda" %in% names(dots)) nlambda <- dots$n.lambda

  ## Set up XX, yy, lambda
  XX <- std(X)
  ns <- attr(XX, "nonsingular")
  penalty.factor <- penalty.factor[ns]
  p <- ncol(XX)
  if (!is.null(X_for_K)) X_for_K <- ncvreg::std(X_for_K)
  yy <- y
  n <- length(yy)
  ##############################################################################
  ### get things working w/o rotation first

  ## Calculate eta
  if (is.null(X_for_K)){
    c(S, U, eta) %<-% lmm_lasso_null(XX, yy) ### change names so doesn't contain lasso
  } else {
    c(S, U, eta) %<-% lmm_lasso_null(X_for_K, yy)
  }
  W <- diag((eta * S + (1 - eta))^(-1/2))

  ## Add manual intercept - is this necessary or will centering SUy work?
  # XXX <- cbind(1, XX)
  # if (missing(penalty.factor) || length(penalty.factor) != ncol(XXX)){
  #   penalty.factor <- rep(c(0, 1), times = c(1, ncol(XX)))
  # }

  ## Rotate data
  SUX <- W %*% crossprod(U, XX)
  SUy_uncentered <- drop(W %*% crossprod(U, yy))
  SUy <- SUy_uncentered - mean(SUy_uncentered)

  ## Set up lambda
  if (missing(lambda)) {
    lambda <- ncvreg::setupLambda(SUX, SUy, "gaussian", alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  fit_ncv <- ncvreg(SUX, SUy)

  # fit <- ncvreg::ncvreg(std(SUX)[,-1], SUy_uncentered) # using this with std(SUX) and no intercept to check
  # fit_gl <- glmnet::glmnet(std(SUX), SUy, intercept = FALSE, standardize = FALSE, penalty.factor = penalty.factor, lambda = lambda)

  ## Set placeholders for results
  b <- matrix(NA, ncol(SUX), nlambda)
  a <- rep(mean(SUy_uncentered), nlambda)
  iter <- integer(nlambda)
  loss <- numeric(nlambda)
  for (ll in 1:nlambda){
    ### need to pass in previous beta values as move from one lambda to the next...duh
    lam <- lambda[ll]
    res <- ncvreg::ncvfit(SUX, SUy, init, penalty, gamma, alpha, lam, eps, max.iter, penalty.factor, warn)
    b[, ll] <- init <- res$beta # does this include intercept?
    iter[ll] <- res$iter
    loss[ll] <- res$loss
  }

  fit_gl <- glmnet(SUX, SUy_uncentered, intercept = TRUE, standardize = FALSE, lambda = lambda)
  fit_lam <- glmnet(SUX, SUy_uncentered, intercept = TRUE, standardize = FALSE) # do setupLambda and glmnet lambda get the same values?
  head(coefficients(fit_gl, s = lambda[16]))

  ### check against glmnet for a particular lambda...
  # checkFun <- function(l){
  #   these <- which(coefficients(fit_gl, s = lambda[l]) != 0)
  #   # print("glmnet:\n")
  #   r1 <- coefficients(fit_gl, s = lambda[l])[these]
  #   those <- which(beta[, l] != 0)
  #   # print("penalizedLMM:\n")
  #   r2 <- beta[those, l]
  #   list(glmnet_res = r1,
  #        plmm_res = r2)
  # }


  ### from ncvreg ##############################################################

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  a <- a[ind]
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & sum(iter) == max.iter) warning("Maximum number of iterations reached")

  ## Local convexity? - later?
  # convex.min <- if (convex) convexMin(b, XX, penalty, gamma, lambda*(1-alpha), family, penalty.factor, a=a) else NULL

  ## Unstandardize
  beta <- matrix(0, nrow=(ncol(XX) + 1), ncol=length(lambda))
  bb <- b/attr(XX, "scale")[ns]
  beta[ns+1,] <- bb
  beta[1,] <- a - crossprod(attr(XX, "center")[ns], bb)

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
                   class = "ncvreg") # penalizedLMM class?
  if (missing(returnX)) {
    if (utils::object.size(SUX) > 1e8) {
      warning("Due to the large size of SUX (>100 Mb), returnX has been turned off.\nTo turn this message off, explicitly specify returnX=TRUE or returnX=FALSE).")
      returnX <- FALSE
    } else {
      returnX <- TRUE
    }
  }
  if (returnX) {
    val$X <- SUX ### return standardized version?
    val$y <- SUy
  }
  return(val)
}


