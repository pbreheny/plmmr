#' Calculate a relatedness matrix
#'
#' This function allows you to generate an n by n genetic relatedness matrix. If a numeric matrix is supplied, the RRM (Hayes, 2009) is used
#' and is computed XX'/p, where X is standardized. 
#' @param X A numeric matrix of genotypes (from *fully-imputed* data)
#' @param std Logical: should X be standardized? If you set this to FALSE, you should know exactly why that is appropriate... standardization is a best practice, and this will impact results.
#' @param fbm Logical: is X stored as an FBM? Defaults to FALSE
#' @param ns Optional vector of values indicating the indices of nonsingular features 
#' @param ... Other optional arguments to `bigstatsr::bigapply()` (like `ncores = ...`)
#' @export
#' 
#' @examples 
#' RRM <- relatedness_mat(X = admix$X)
relatedness_mat <- function(X, std = TRUE, fbm = FALSE, ns = NULL, ...){
  if (fbm){
    if (!std) stop("\nAt this time, standardization cannot be 'turned off' for filebacked data (because of the default behavior of process_plink()). 
                          \nWe are working to extend this to a more general filebacked case.")
    # cross product 
    cprod <- bigstatsr::big_tcrossprodSelf(X,
                                           # TODO: figure out why the line below makes the function return all NA values
                                           # fun.scaling = bigstatsr::as_scaling_fun(scale.col = rep(X$ncol, length(ns)),
                                           #                                         center.col = rep(0, length(ns))),
                                           ind.col = ns,
                                           ...)
    rrm <- bigstatsr::FBM(cprod$nrow, cprod$ncol)
    # scale by p, the number of columns in the design matrix (including constant features)
    bigstatsr::big_apply(X = cprod,
                         a.FUN = function(X, ind, p, res){
                           res[,ind] <- sweep(x = X[,ind],
                                              MARGIN = 2,
                                              STATS = p,
                                              FUN = "/")},
                         a.combine = cbind,
                         p = X$ncol,
                         res = rrm)
    
  } else {
    if (std){
      rrm <- tcrossprod(ncvreg::std(X))/ncol(X)
    } else {
      rrm <- tcrossprod(X)/ncol(X)
    }
    
  }
 
 return(rrm)
}

