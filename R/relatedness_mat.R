#' Calculate a relatedness matrix
#'
#' This function allows you to generate an n by n genetic relatedness matrix.
#' If a numeric matrix is supplied, the RRM [(Hayes, 2009)](https://doi.org/10.1017/S0016672308009981) is used
#' and is computed XX'/p, where X is standardized.
#'
#' @param X An n x p numeric matrix of genotypes (from *fully-imputed* data).
#' Note: This matrix should *not* include non-genetic features.
#' @param std Logical: should X be standardized? If you set this to FALSE (which can only be done if data are stored in memory), you should know exactly why that is appropriate... standardization is a best practice, and this will impact results.
#' @param fbm Logical: is X stored as an FBM? Defaults to FALSE
#' @param ns Optional vector of values indicating the indices of nonsingular features
#' @param ... Other optional arguments to `bigstatsr::bigapply()` (like `ncores = ...`)
#'
#' @returns An n x n numeric matrix capturing the genomic relatedness of the
#' samples represented in `X`. In our notation, we call this matrix K for 'kinship';
#' this matrix is also sometimes called the GRM for 'genomic relatedness matrix'.
#'
#' @export
#'
#' @examples
#' RRM <- relatedness_mat(X = admix$X)
#' RRM[1:5, 1:5]
relatedness_mat <- function(X, std = TRUE, fbm = FALSE, ns = NULL, ...){
  if (fbm){
    if (!std) stop("\nAt this time, standardization cannot be 'turned off' for
    filebacked data (because of the default behavior of process_plink()).")
    # cross product
    if (length(ns) %in% 1:(ncol(X)-1)){
      cprod <- bigstatsr::big_tcrossprodSelf(X,
                                             ind.col = ns,
                                             ...)
    } else {
      cprod <- bigstatsr::big_tcrossprodSelf(X,
                                             ...)
    }

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
                         res = rrm,
                         ncores = count_cores())

  } else {
    if (std){
      rrm <- tcrossprod(ncvreg::std(X))/ncol(X)
    } else {
      rrm <- tcrossprod(X)/ncol(X)
    }

  }

 return(rrm)
}

