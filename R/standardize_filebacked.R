#' A helper function to standardize a filebacked matrix
#'
#' @param X             A `big.matrix` object that has been subset &/or had any additional predictors appended as columns
#' @param outfile       Optional: the name (character string) of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet         Logical: should console messages be silenced? Defaults to FALSE
#' @param tocenter      Should the matrix be centered in addition to scaled? Defaults to TRUE.
#'
#' @return A list with a component called `std_X` - this is an FBM with column-standardized data.
#' List also includes several other indices/meta-data on the standardized matrix
#'
#' @keywords internal
#'
standardize_filebacked <- function(X, outfile, quiet, tocenter = TRUE) {

  # standardization ------------------------------------------------
  if (!quiet) {
    cat("Column-standardizing the design matrix...\n")
  }
  # centering & scaling
  # NOTE: this C++ call will change the .bk file so that its data are column-standardized
  std_res <- .Call("big_std",
                   X@address,
                   as.integer(count_cores()),
                   tocenter,
                   NULL, # no center values to pass here -- will calculate these
                   NULL, # no scaling values to pass -- will calculate these
                   PACKAGE = "plmmr")
  X@address <- std_res$std_X # saves standardized .bk

  if (!quiet) {
    cat("Standardization completed at", pretty_time())
  }

  cat("Standardization completed at", pretty_time(), file = outfile, append = TRUE)

  # label return object ------------------------------------------------
  # naming these center and scale values so that I know they relate to the first
  # standardization; there will be another standardization after the rotation
  # in plmm_fit().
  if (!quiet) {
    cat("Done with standardization. File formatting in progress\n")
  }

  list(
    std_X = bigmemory::describe(X),
    std_X_n = nrow(X),
    std_X_p = ncol(X),
    std_X_center = std_res$std_X_center,
    std_X_scale = std_res$std_X_scale)
}
