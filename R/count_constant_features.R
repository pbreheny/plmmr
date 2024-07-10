#' A helper function to count constant features
#'
#' @param fbm An FBM object, usually the 'genotypes' component from a bigSNP list
#' @param ind.row Optional indicators for rows in fbm. Default includes all rows.
#' @param outfile String specifying name of log file
#' @param quiet Logical: should a message be printed to the console
#' @return ns A numeric vector with the indicesof the non-singular columns of the matrix associated with `counts`
#' @keywords internal
#'
count_constant_features <- function(fbm, ind.row = bigstatsr::rows_along(fbm), outfile, quiet){

  # NB: pruning out samples with incomplete phenotypes can make some features
  #   *become* constant!
  colstats <- bigstatsr::big_colstats(fbm, ind.row = ind.row)
  ns <- which(colstats$var > 1e-4)
  constants_idx <- sum(colstats$var < 1e-4)

  if(!quiet){
    cat("There are", sum(constants_idx), "constant features in the data\n")
  }

  cat("There are", sum(constants_idx), "constant features in the data\n",
      file = outfile, append = TRUE)

  return(ns)
}
