#' A function to take the eigendecomposition of K
#' Note: This is faster than taking SVD of X when p >> n 
#'
#' @param std_X The *standardized* design matrix 
#' @param p The number of columns in the *unstandardized* design matrix
#' @return A list with the eigenvectors and eigenvalues of K
#' @keywords internal
#'
eigen_K <- function(std_X, p){
  # idea: add @param k Optional integer argument indicating the number of singular values to use in a truncated SVD. Default (NULL) does full decomposition.
  K <- tcrossprod(std_X)/p
  decomp <- eigen(K)
  return(list(s = decomp$values,
              U = decomp$vectors,
              K = K))
}