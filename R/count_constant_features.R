#' A helper function to count constant features
#'
#' @param fbm An FBM object, usually the 'genotypes' component from a bigSNP list
#' @param outfile String specifying name of log file
#' @param quiet Logical: should a message be printed to the console
#' @return ns A numeric vector with the indices of the non-singular columns of the matrix associated with `counts`
#' @keywords internal
#'
count_constant_features <- function(fbm, outfile, quiet){

  # NB: pruning out samples with incomplete phenotypes can make some features
  #   *become* constant!
  colstats <- .Call('big_sd',
                    fbm@address,
                    as.integer(bigstatsr::nb_cores()),
                    PACKAGE = 'plmmr')
  ns <- which(colstats$sd_vals > 1e-4)
  constants_idx <- sum(colstats$sd_vals < 1e-4)

  if(!quiet){
    cat("There are", sum(constants_idx), "constant features in the data\n")
  }

  cat("There are", sum(constants_idx), "constant features in the data\n",
      file = outfile, append = TRUE)

  return(ns)
}
