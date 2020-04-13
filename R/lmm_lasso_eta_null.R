#' Estimate eta using a null LmmLasso model
#'
#' This function allows you to estimate eta (the signal to noise ratio, or narrow-sense variability) under the assumption of a null model.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @keywords
#' @export
#' @examples

lmm_lasso_eta_null<-function(X, y){
  K <- tcrossprod(std(X))/ncol(X)
  c(S, U, V) %<-% svd(K)
  Uy <- crossprod(U, y)
  opt <- optimize(f=lmm_lasso_eta_nll, c(0, 1), Uy=Uy, S=S)
  eta <- opt$minimum
  return(list(S=S, U=U, eta=eta))
}
