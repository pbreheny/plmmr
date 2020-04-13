#' Fit an LMM with lasso regularization
#'
#' This function allows you to fit an approximate LMM (using the rotation method) with lasso regularization.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param standardize Should standardization be performed within glmnet()? Defaults to FALSE.
#' @keywords
#' @export
#' @examples

lmm_lasso_eta <- function(X, y, p1, standardize = FALSE) {
  require(stats)
  require(zeallot)
  c(S, U, eta) %<-% lmm_lasso_eta_null(X, y)
  W <- diag((eta * S + (1 - eta))^(-1/2))
  SUX <- W %*% crossprod(U, X)
  SUy <- drop(W %*% crossprod(U, y))
  fit <- glmnet(SUX, SUy, standardize = standardize)
  sel <- sapply(predict(fit, type='nonzero'), length)
  coef <- coef(fit, min(fit$lambda[sel <= p1]))[-1]
  return(list(fit = fit,
              nonzero = length(which(coef != 0)),
              coef = coef,
              delta = (1/eta) - 1,
              eta = eta))
}
