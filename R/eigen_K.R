#' A function to take the eigendecomposition of K
#' Note: This is faster than taking SVD of X when p >> n 
#'
#' @param std_X The *standardized* design matrix, stored as FBM.  
#' @param p The number of columns in the *unstandardized* design matrix. 
#' @param fbm_flag Logical: is std_X an FBM obejct? Passed from `plmm()`.
#' @param ... Optional additional arguments to `bigstatsr::big_tcrossprodSelf()`,
#' Of note, one additional argument is `ind.col`, which is used to *exclude* the 
#' unpenalized additional, non-genomic covariates (like sex, age, etc.) from being
#' included in the calculation of kinship. 
#' @return A list with the eigenvectors and eigenvalues of K
#' @keywords internal
#'
eigen_K <- function(std_X, p, fbm_flag, ...){
  # Note: std_X has already been scaled, so no need to do that here
  # calculate K (which will be stored in memory regardless of how std_X is stored)
  if(fbm_flag){
    K <- bigstatsr::big_tcrossprodSelf(std_X, ...) 
    # TODO: this function returns an FBM. Need to think about adding the option 
    # for K to be stored filebacked. For now, make K stay in memory.
     K <- K[,]/p
  } else {
    # make std_X an FBM so that bigstatsr::big_tcrossprodSelf can be used... 
    # need this for subsetting columns & passing additional arguments
    std_X <- bigstatsr::as_FBM(std_X) 
    K <- bigstatsr::big_tcrossprodSelf(std_X, ...) 
    K <- K[,]/p
  }
  
  # take eigendecomposition
  decomp <- eigen(K)
  return(list(s = decomp$values,
              U = decomp$vectors,
              K = K))
}