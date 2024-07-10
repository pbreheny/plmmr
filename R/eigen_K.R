#' A function to take the eigendecomposition of K
#' Note: This is faster than taking SVD of X when p >> n
#'
#' @param std_X The *standardized* design matrix, stored as big.matrix object.
#' @param fbm_flag Logical: is std_X an FBM obejct? Passed from `plmm()`.
#' @param ... Optional additional arguments to `bigstatsr::big_tcrossprodSelf()`,
#' Of note, one additional argument is `ind.col`, which is used to *exclude* the
#' unpenalized additional, non-genomic covariates (like sex, age, etc.) from being
#' included in the calculation of kinship.
#' @return A list with the eigenvectors and eigenvalues of K
#' @keywords internal
#'
eigen_K <- function(std_X, fbm_flag, ...){
  # Note: std_X has already been scaled, so no need to do that here
  # calculate K (which will be stored in memory regardless of how std_X is stored)

  if (fbm_flag) {
    f <- paste0(bigmemory::dir.name(std_X), bigmemory::file.name(std_X))
    bk_file_path_sans_extension <- bigstatsr::sub_bk(f)
    fbm_X <- bigstatsr::FBM(nrow = nrow(std_X), ncol = ncol(std_X),
                            create_bk = FALSE,
                            backingfile = bk_file_path_sans_extension)
    K <- bigstatsr::big_tcrossprodSelf(fbm_X, ...)
    K <- K[,]/ncol(std_X)
  } else {
    # make std_X an FBM so that bigstatsr::big_tcrossprodSelf can be used...
    # need this for subsetting columns & passing additional arguments
    std_X <- bigstatsr::as_FBM(std_X)
    K <- bigstatsr::big_tcrossprodSelf(std_X, ...)
    K <- K[,]/ncol(std_X)
  }

  # take eigendecomposition
  # TODO: explore ways to make this faster
  decomp <- eigen(K)
  return(list(s = decomp$values,
              U = decomp$vectors,
              K = K))
}
