#' a function to create the estimated variance matrix from a PLMM fit 
#' @param fit An object returned by `plmm()`
#' @param K An optional matrix 
#' @param eta An optional numeric value between 0 and 1; if `fit` is not supplied, then this option must be specified.
#' @returns Sigma_hat, a matrix representing the estimated variance 
#' 
#' @export
#' 
construct_variance <- function(fit, K = NULL, eta = NULL){
  # TODO: check the inputs for this function...
  # if 'fit' is given 
  if (!missing(fit)){
      # note: K is a matrix 
      Sigma_hat <- (fit$eta * fit$K) + ((1-fit$eta) * diag(nrow(fit$K)))
  } else if (!is.null(K) & !is.null(eta)) {
    Sigma_hat <- (eta * K) + ((1-eta) * diag(nrow(K)))
  } else {
    stop("\nOptions to construct_variance must be either (1) supply a plmm object to 
         fit or \n(2) supply both K and eta arguments.")
  }
  
  return(Sigma_hat)
}

