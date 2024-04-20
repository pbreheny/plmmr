#' a function to create the estimated variance matrix from a PLMM fit 
#' @param fit An object returned by `plmm()`
#' @param K An optional matrix or list as returned by `choose_K()`
#' @param eta An optional numeric value between 0 and 1; if `fit` is not supplied, then this option must be specified.
#' @return Sigma_hat, a matrix representing the estimated variance 
#' 
#' @export
#' 
construct_variance <- function(fit, K = NULL, eta = NULL){
  
  # if 'fit' is given 
  if (!missing(fit)){
    # case 1: K is a list
    if(is.list(K)){
      K <- K$U %*% tcrossprod(diag(K$s), K$U) # TODO: make this computation more efficient
      Sigma_hat <- fit$eta * K + (1-fit$eta) * diag(fit$n)
    } else {
      # case 2: K is a matrix 
      Sigma_hat <- (fit$eta * fit$K) + ((1-fit$eta) * diag(nrow(fit$K)))
    }
  } else if (!is.null(K) & !is.null(eta)) {
    # case 1: K is a list
    if(is.list(K)){
      K <- K$U %*% tcrossprod(diag(K$s), K$U)
      Sigma_hat <- eta * K + (1-eta) * diag(n)
    } else {
      # case 2: K is a matrix 
      Sigma_hat <- (eta * K) + ((1-eta) * diag(nrow(K)))
    }
    
  } else {
    stop("\nOptions to construct_variance must be either (1) supply a plmm object to 
         fit or \n(2) supply both K and eta arguments.")
  }
  
  return(Sigma_hat)
}

