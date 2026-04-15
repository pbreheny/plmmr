#' a function to compute the BLUP
#' @param fit An object returned by `plmm()`
#' @param Xb Linear predictor
#' @param Sigma_21 Covariance matrix between the training and the testing data. Extracted from `estimated_Sigma` that is generated using all observations
#' @param idx Vector of indices of the penalty parameter \code{lambda} at which predictions are required. By default, all indices are returned.
#'
#' @returns Sigma_hat, a matrix representing the estimated variance
#'
#' @keywords internal
#'

compute_blup <- function(fit, Xb, Sigma_21, idx) {
  # if 'fit' is given
  if (!missing(fit)) {
    K <- fit$K
    eta <- fit$eta
    if (is.matrix(K)) {
      # case 1: K is a matrix
      decomp <- eigen(K, symmetric = TRUE)
      nz <- decomp$values > 1e-4
      U <- decomp$vectors[, nz, drop = FALSE]
      s <-  decomp$values[nz]
    } else {
      # case 2: K is a list with U,s
      U <- K$U
      s <- K$s
    }
  } else {
    stop("\nMust supply a plmm object.", call. = FALSE)
  }

  resid_old <- drop(fit$y) - fit$std_Xbeta[, idx]
  Ut_r   <- drop(t(U) %*% resid_old)
  proj_r <- U %*% Ut_r
  tmp <- U %*% (Ut_r / (eta * s + (1 - eta)))
  ranef <- Sigma_21 %*% tmp
  blup <- Xb + ranef

  return(blup)
}
