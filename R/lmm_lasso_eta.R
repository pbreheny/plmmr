#' Fit a linear mixed model with lasso regularization
#'
#' This function allows you to fit a linear mixed model via lasso-penalized maximum likelihood.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param standardize Should standardization be performed within \code{glmnet()}? Defaults to FALSE.
#' @param X_for_K X matrix used to compute the similarity matrix, K. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @importFrom zeallot %<-%
#' @export



lmm_lasso_eta <- function(X, y, p1, X_for_K = NULL, standardize = FALSE) {
  S <- U <- eta <- NULL
  if (is.null(X_for_K)){
    c(S, U, eta) %<-% lmm_lasso_eta_null(X, y)
  } else {
    c(S, U, eta) %<-% lmm_lasso_eta_null(X_for_K, y)
  }
  W <- diag((eta * S + (1 - eta))^(-1/2))
  SUX <- W %*% crossprod(U, X)
  SUy <- drop(W %*% crossprod(U, y))
  fit <- glmnet::glmnet(SUX, SUy, standardize = standardize)
  sel <- sapply(stats::predict(fit, type='nonzero'), length)
  coef <- coef(fit, min(fit$lambda[sel <= p1]))[-1]
  names(coef) <- colnames(X)
  return(list(fit = fit,
              nonzero = length(which(coef != 0)),
              coef = coef,
              delta = (1/eta) - 1,
              eta = eta))
}
