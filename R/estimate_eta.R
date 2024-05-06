#' Estimate eta (to be used in rotating the data)
#' This function is called internally by \code{plmm()}
#' @param n The number of observations 
#' @param s The singular values of K, the realized relationship matrix
#' @param U The left-singular vectors of the *standardized* design matrix
#' @param y Continuous outcome vector.
#' 
#' @keywords internal
estimate_eta <- function(n, s, U, y, eta_star){

  opt <- stats::optimize(f=log_lik,
                         c(0.01, 0.99),
                         n = n,
                         s = s,
                         U = U,
                         y = y)
  
  eta <- opt$minimum 
  
  
  return(eta)
}
