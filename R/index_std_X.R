#' Helper function to index standardized data
#'
#' @param std_X_p The number of features in the *standardized* matrix of data (which may be filebacked)
#' @param non_genomic Integer vector of columns in `std_X` representing non-genomic data.
#'
#' @return list with indices
#'
#' @keywords internal
#'
index_std_X <- function(std_X_p, non_genomic){
  # keep track of non-genomic features
  col_idx <- 1:std_X_p
  if (is.null(non_genomic)) {
    genomic <- col_idx
  } else {
    genomic <- col_idx[-non_genomic]
  }

  return(genomic)
}