#' A helper function to count constant features 
#'
#' @param counts An object returned by `bigstatsr::big_counts()`
#' @param outfile String specifying name of log file 
#' @param quiet Logical: should a message be printed to the console ? 
#'
#' @return ns A numeric vector with the indicies of the non-singular columns of the matrix associated with `counts`
#' @keywords internal
#'
count_constant_features <- function(counts, outfile, quiet){
  # NB: pruning out samples with incomplete phenotypes can make some features
  #   *become* constant! 
  constants_idx <- apply(X = counts[1:3,],
                         MARGIN = 2,
                         # see which ~called~ features have all same value
                         FUN = function(c){sum(c == sum(c)) > 0})
  
  ns <- which(!constants_idx) # need this for analysis downstream
  
  
  if(!quiet){
    cat("\nThere are ", sum(constants_idx), " constant features in the data",
        file = outfile, append = TRUE)
    
    cat("\nThere are ", sum(constants_idx), " constant features in the data")
  }
  
  
  return(ns)
}