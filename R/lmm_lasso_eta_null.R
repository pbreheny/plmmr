#' Estimate eta using a null LmmLasso model
#'
#' This function allows you to estimate eta (the signal to noise ratio, or narrow-sense variability) under the assumption of a null model in order to estimate the variance parameter, eta.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @keywords
#' @export
#' @examples

lmm_lasso_eta_null<-function(X, y){
  require(stats)
  require(zeallot)
  K <- tcrossprod(std(X))/ncol(X)
  c(S, U, V) %<-% svd(K)
  Uy <- crossprod(U, y)
  opt <- optimize(f=lmm_lasso_eta_nll, c(0, 1), Uy=Uy, S=S)
  # opt <- optim(par = 0.5, fn = lmm_lasso_eta_nll, Uy = Uy, S = S, method = "L-BFGS-B", lower = -Inf, upper = 0)
  eta <- opt$minimum
  return(list(S=S, U=U, eta=eta))
}
