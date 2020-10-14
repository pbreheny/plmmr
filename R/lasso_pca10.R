#' Fit a linear model with lasso regularization and PC adjutment
#'
#' This function allows you to fit a linear model via lasso-penalized maximum likelihood, adjusted for the first 10 principal components (unpenalized), with the same output values as those of lmm_lasso_eta. Primarily used for simulation and comparison with lmm_lasso_eta.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param standardize Should standardization be performed within \code{glmnet()}? Defaults to FALSE.
#' @importFrom zeallot %<-%
#' @export


lasso_pca10 <- function(X, y, p1, standardize = FALSE){
  S <- U <- NULL
  K <- tcrossprod(ncvreg::std(X))/ncol(X)
  c(S, U) %<-% methods::as(eigen(K), "list")
  pc <- U[, 1:10]
  XX <- cbind(pc, X)
  fit <- glmnet::glmnet(XX, y, penalty.factor = rep(c(0, 1), times=c(10, ncol(X))), standardize = standardize)
  sel <- sapply(stats::predict(fit, type='nonzero'), length) - 10 # remove unpenalized effects
  coef <- coef(fit, min(fit$lambda[sel <= p1]))[-c(1:11)] # remove intercept and unpenalized effects
  names(coef) <- colnames(X)
  coef_pred <- coef(fit, min(fit$lambda[sel <= p1]))
  names(coef_pred) <- c('(Intercept)', paste0('PC', 1:10), colnames(X))
  return(list(fit = fit,
              nonzero = length(which(coef != 0)),
              coef = coef,
              coef_pred = coef_pred,
              delta = NA,
              eta = NA))
}
