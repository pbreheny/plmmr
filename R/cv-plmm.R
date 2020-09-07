#' Cross-validation for plmm
#'
#' Performs k-fold cross validation for lasso-, MCP-, or SCAD-penalized penalized linear mixed models over a grid of values for the regularization parameter lambda.
#' @param X Design matrix. May include clinical covariates and other non-SNP data. If this is the case, X_for_K should be supplied witha  matrix containing only SNP data for computation of GRM.
#' @param y Continuous outcome vector.
#' @param X_for_K X matrix used to compute the similarity matrix, K. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @param ... Additional arguments to plmm
#' @param cluster cv.plmm can be run in parallel across a cluster using the parallel package. The cluster must be set up in advance using the makeCluster function from that package. The cluster must then be passed to cv.plmm.
#' @param nfolds The number of cross-validation folds. Default is 10.
#' @param fold Which fold each observation belongs to. By default the observations are randomly assigned.
#' @param seed You may set the seed of the random number generator in order to obtain reproducible results.
#' @param returnY Should cv.plmm return the linear predictors from the cross-validation folds? Default is FALSE; if TRUE, this will return a matrix in which the element for row i, column j is the fitted value for observation i from the fold in which observation i was excluded from the fit, at the jth value of lambda.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each CV fold. Default is FALSE.
#' @export


cv.plmm <- function(X, y, X_for_K = X, ..., cluster, nfolds=10, seed, fold,
                    returnY=FALSE, trace=FALSE) {

  # Coercion
  if (!is.matrix(X)) {
    tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (!is.double(y)) {
    tmp <- try(y <- as.double(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be numeric or able to be coerced to numeric", call.=FALSE)
  }

  fit.args <- c(list(X = X, y = y, X_for_K = X_for_K), list(...))
  fit.args$returnX <- TRUE
  fit <- do.call('plmm', fit.args)
  # fit <- plmm(X=X, y=y, X_for_K=X_for_K, fit.args)

  # Extract rotated X and y and do CV on the rotated objects
  XX <- fit$X
  yy <- fit$y

  n <- length(y)
  E <- Y <- matrix(NA, nrow=n, ncol=length(fit$lambda))

  if (!missing(seed)) set.seed(seed)
  sde <- sqrt(.Machine$double.eps)
  if (missing(fold)) {
    fold <- sample(1:n %% nfolds)
    fold[fold==0] <- nfolds
  } else {
    nfolds <- max(fold)
  }

  cv.args <- list(...)
  cv.args$penalty.factor <- rep(c(0, 1), times = c(ncol(XX) - ncol(X_for_K), ncol(X_for_K)))
  cv.args$lambda <- fit$lambda
  cv.args$warn <- FALSE
  cv.args$convex <- FALSE
  cv.args$rotation <- FALSE
  cv.args$intercept <- FALSE
  cv.args$standardize <- FALSE # is this right? Do I need to do something special with standardization?
  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("fold","fit","X", "y", "cv.args"), envir=environment())
    parallel::clusterCall(cluster, function() library(penalizedLMM))
    fold.results <- parallel::parLapply(cl=cluster, X=1:nfolds, fun=cvf, X=XX, y=yy, fold=fold, cv.args=cv.args)
  }

  for (i in 1:nfolds) {
    if (!missing(cluster)) {
      res <- fold.results[[i]]
    } else {
      if (trace) cat("Starting CV fold #", i, sep="","\n")
      res <- cvf(i, XX, yy, fold, cv.args)
    }
    E[fold==i, 1:res$nl] <- res$loss
    Y[fold==i, 1:res$nl] <- res$yhat
  }

  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]

  ## Return
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, stats::sd) / sqrt(n)
  min <- which.min(cve)

  # Bias correction
  e <- sapply(1:nfolds, function(i) apply(E[fold==i, , drop=FALSE], 2, mean))
  Bias <- mean(e[min,] - apply(e, 2, min))

  val <- list(cve=cve, cvse=cvse, fold=fold, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min],
              null.dev=mean(loss.plmm(y, rep(mean(y), n))), Bias=Bias)
  if (returnY) val$Y <- Y
  structure(val, class="cv.plmm")
}