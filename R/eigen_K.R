#' A function to take the eigendecomposition of K
#'
#' Note: This is faster than taking SVD of X when p >> n
#'
#' @param std_X The *standardized* design matrix.
#'
#' @return A list with three elements:
#' * `s`: The non-zero eigenvalues of K
#' * `U`: The eigenvectors of K associated with s
#' * `K`: The fully computed K matrix
#'
#' @keywords internal
#'
eigen_K <- function(std_X) {
  # Note: std_X has already been scaled, so no need to do that here
  # calculate K (which will be stored in memory regardless of how std_X is stored)
  K <- relatedness_mat(std_X, std = FALSE)

  # take eigendecomposition
  decomp <- eigen(K, symmetric = TRUE)
  nz <- decomp$values > 1e-4
  list(
    s = decomp$values[nz],
    U = decomp$vectors[, nz, drop = FALSE],
    K = K
  )
}
