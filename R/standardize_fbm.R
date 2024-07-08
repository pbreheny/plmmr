#' A helper function to standardize an `FBM`
#'
#' @param subset_X An FBM
#' @param prefix   The data file name without the extension
#' @param rds_dir  The directory where the user wants to create the '.rds' and '.bk' files
#' @param ns       Integer vector with indices of non-singular columns of `subset_X`
#' @param non_gen  Integer vector of indices with non-genomic covariates
#' @param outfile  Optional: the name (character string) of the prefix of the logfile to be written.
#' @param quiet    Logical: should the messages printed to the console be silenced? Defaults to FALSE.
#'
#' @returns A list that includes:
#' * 'std_X' - the standardized FBM
#' * 'center' - a vector with  the values used to center the columns
#' * 'scale' - a vector with the values used to scale the columns
#'
#' @keywords internal
#'
standardize_fbm <- function(subset_X, prefix, rds_dir, ns, non_gen,
                           outfile, quiet){

  # standardization ------------------------------------------------
  if (!quiet) {cat("Column-standardizing the design matrix...\n")}
  # convert FBM pointer into a big.matrix pointer
  subset_X_bm <- subset_X |> fbm2bm()
  # centering & scaling
  std_res <- .Call("big_std",
                   subset_X_bm@address,
                   as.integer(bigstatsr::nb_cores()),
                   PACKAGE = "plmmr")

  std_X <- bigmemory::big.matrix(nrow = nrow(subset_X), ncol = ncol(subset_X))
  std_X@address <- std_res$std_X

  # label return object --------------------------------------------------------
  # naming these center and scale values so that I know they relate to the first
  # standardization; there will be another standardization after the rotation
  # in plmm_fit().
  ret <- list(
    std_X = bigmemory::describe(std_X),
    std_X_center = std_res$std_X_center,
    std_X_scale = std_res$std_X_scale,
    ns = ns,
    non_gen = non_gen) # save indices for non-genomic covariates

  if (!quiet){
    cat("Done with standardization. File formatting in progress.\n",
        file = outfile, append = TRUE)
  }

  return(ret)
}