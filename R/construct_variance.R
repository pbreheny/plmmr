#' a function to create the estimated variance matrix from a PLMM fit
#' @param fit An object returned by `plmm()`
#' @param K An optional matrix
#' @param eta An optional numeric value between 0 and 1; if `fit` is not supplied, then this option must be specified.
#'
#' @returns Sigma_hat, a matrix representing the estimated variance
#'
#' @keywords internal
#'
construct_variance <- function(fit, K = NULL, eta = NULL){

  # if 'fit' is given
  if (!missing(fit)){
    K <- fit$K
      # case 1: K is a matrix
    if (is.matrix(K)){
      Sigma_hat <- (fit$eta * K) + ((1-fit$eta) * diag(nrow(K)))
    } else {
      # case 2: K is a list with U,s
      SUt <- sweep(t(K$U), MARGIN = 1, STATS = K$s, FUN = "*")
      K <- K$U%*%SUt
      Sigma_hat <- (fit$eta * K) + ((1-fit$eta) * diag(nrow(K)))
    }
  } else if (!is.null(K) & !is.null(eta)) {
    if (is.matrix(K)){
      Sigma_hat <- (eta * K) + ((1-eta) * diag(nrow(K)))
    } else {
      SUt <- sweep(t(K$U), MARGIN = 1, STATS = K$s, FUN = "*")
      K_mat <- K$U%*%SUt
      Sigma_hat <- (eta * K_mat) + ((1-eta) * diag(nrow(K$U)))
    }

  } else {
    stop("\nOptions to construct_variance must be either: \n(1) supply a plmm object to
         fit or \n(2) supply both K and eta arguments.")
  }

  return(Sigma_hat)
}

