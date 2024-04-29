#' big_std
#'A function to standardize a matrix stored as an FBM object 
#'
#' @param X An file-backed design matrix 
#' @param center A vector of centering values
#' @param scale A vector of scaling values 
#' @param ns A vector with the indices marking the Non-Singular columns of X. Defaults to NULL (meaning all columns are assumed to be nonsingular).
#' @param std_bk_extension The string specifying the backingfile for the standardized copy of X. Default uses a temporary directory.
#'
#' @return A standardized matrix of class FBM
#' @keywords internal
big_std <- function(X, center = NULL, scale, ns = NULL, std_bk_extension = NULL){
  
  if (is.null(ns)) {
    ns <- bigstatsr::cols_along(X)
  }
  
  if(is.null(center)){
    centered_X <- X
  } else {
    # allocate space 
    centered_X <- bigstatsr::FBM(X$nrow, X$ncol)
    bigstatsr::big_apply(X = X,
                         a.FUN = function(X, ind, center, res){
                           res[,ind] <- sweep(x = X[,ind],
                                              MARGIN = 2,
                                              STATS = center[ind],
                                              FUN = "-")},
                         a.combine = cbind,
                         ncores = bigstatsr::nb_cores(),
                         center = center,
                         res = centered_X)
  }
  
  
  if (!is.null(std_bk_extension)) { 
    # Note: this is the case used in `process_plink()` and its helper funs
    scaled_X <- bigstatsr::FBM(nrow = centered_X$nrow,
                               ncol = centered_X$ncol,
                               backingfile = std_bk_extension)
  } else {
    scaled_X <- bigstatsr::FBM(nrow = centered_X$nrow,
                               ncol = centered_X$ncol)
  }
  
  bigstatsr::big_apply(X = centered_X, 
                       a.FUN = function(X, ind, scale, res){
                         res[,ind] <- sweep(x = X[,ind],
                                            MARGIN = 2,
                                            STATS = scale[ind],
                                            FUN = "/")},
                       a.combine = cbind,
                       # NB: only scale the nonsingular columns
                       ind = ns, 
                       ncores = bigstatsr::nb_cores(),
                       scale = scale,
                       res = scaled_X)
  
  
  
  return(scaled_X)
  
}
