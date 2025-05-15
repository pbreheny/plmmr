#' A function to take the eigendecomposition of K
#' Note: This is faster than taking SVD of X when p >> n
#'
#' @param std_X The *standardized* design matrix, stored as big.matrix object.
#' @param fbm_flag Logical: is std_X an FBM obejct? Passed from `plmm()`.
#'
#' @return A list with the eigenvectors and eigenvalues of K
#'
#' @keywords internal
#'
eigen_K <- function(std_X, fbm_flag) {
  # Note: std_X has already been scaled, so no need to do that here
  # calculate K (which will be stored in memory regardless of how std_X is stored)
  if (fbm_flag) {
    # TODO: explore ways to make this multiplication faster
    # TODO: explore ways to exclude specified columns from XX (e.g, in GWAS, maybe non-genomic predictors should be excluded)
    XX <- bigalgebra::dgemm(TRANSA = "N",
                            TRANSB = "T",
                            A = std_X,
                            B = std_X)
    K <- XX[,] / ncol(std_X)
  } else {
    XX <- tcrossprod(std_X)
    K <- XX / ncol(std_X)
  }

  # take eigendecomposition
  # TODO: explore ways to make this faster
  decomp <- eigen(K)
  return(list(s = decomp$values,
              U = decomp$vectors,
              K = K))
}
