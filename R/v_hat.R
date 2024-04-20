#' a function to create the estimated variance matrix from a PLMM fit 
#' @param fit An object returned by `plmm()`
#' @param K An optional matrix or list as returned by `choose_K()`

#' @return Sigma_hat, a matrix representing the estimated variance 
#' 
#' @export
#' 
construct_variance <- function(fit, K = NULL, eta = NULL){
 
 
  # if K is supplied:
  if (!is.null(K) & !is.null(eta)) {
    # case 1: K is from K_svd list 
    if(is.list(K)){
      K <- K$U %*% tcrossprod(diag(K$s), K$U)
      Sigma_hat <- fit$eta * K + (1-fit$eta) * diag(fit$n)
    } else if (is.matrix(K)){
    # case 2: K is a matrix 
    Sigma_hat <- fit$eta * K + (1-fit$eta) * diag(fit$n)
    }
  }
  
  # if K is not supplied:
  if(is.null(K) & is.list(fit)){
    
    if ("K" %in% names(fit)){
      Sigma_hat <- fit$eta * fit$K + (1-fit$eta) * diag(fit$std_X_n)
    } else {
      stop("\nK is not supplied in either the fit or K arguments.")
    }
    
    
  }
  
  
    return(Sigma_hat)
}
 
