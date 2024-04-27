#' Helper function to index standardized data
#'
#' @param std_X The standardized matrix of data (which may be filebacked)
#' @param non_genomic Integer vector of columns in `std_X` representing non-genomic data. 
#'
#' @return list with indicies
#' @keywords internal
#'
index_std_X <- function(std_X, non_genomic){
  # designate dimensions of the standardized data 
  std_X_n <- nrow(std_X)
  std_X_p <- ncol(std_X)
  
  # keep track of non-genomic features for non-FBM case
  col_idx <- 1:std_X_p
  if (is.null(non_genomic)) {
    genomic <- col_idx
  } else {
    genomic <- col_idx[-non_genomic]
  }
  
  return(list(genomic = genomic, 
              std_X_n = std_X_n,
              std_X_p = std_X_p))
}