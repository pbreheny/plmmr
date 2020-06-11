#' Estimate eta using a null LmmLasso model
#'
#' This function allows you to estimate eta (the signal to noise ratio, or narrow-sense variability) under the assumption of a null model.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @importFrom zeallot %<-%
#' @export

lmm_lasso_null<-function(X, y){
  S <- U <- V <- NULL
  K <- tcrossprod(ncvreg::std(X))/ncol(X)
  c(S, U) %<-% methods::as(eigen(K), "list")
  Uy <- crossprod(U, y)
  opt <- stats::optimize(f=lmm_lasso_nll, c(0.01, 0.99), Uy=Uy, S=S)
  eta <- opt$minimum
  return(list(S=S, U=U, eta=eta))
}
