#' A helper function to standardize a `bigSNP`
#'
#' @param X           A list that includes:
#'                       (1) subset_X: a `big.matrix` object that has been subset &/or had any additional predictors appended as columns
#'                       (2) ns: a numeric vector indicating the indices of nonsingular columns in subset_X
#' @param new_file        The new_file (as a character string) of the bed/fam data files (e.g., `new_file = 'mydata'`)
#' @param rds_dir       The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`
#' @param outfile       Optional: the name (character string) of the new_file of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet         Logical: should messages be printed to the console? Defaults to FALSE (which leaves the print messages on...)
#' @param overwrite     Logical: if existing `.bk`/`.rds` files exist for the specified directory/new_file, should these be overwritten?
#'
#' @return A list with a new component of `obj` called 'std_X' - this is an FBM with column-standardized data.
#' List also includes several other indices/meta-data on the standardized matrix
#' @keywords internal
#'
standardize_bigsnp <- function(X, new_file, rds_dir, non_gen, complete_outcome, id_var,
                               outfile, quiet, overwrite){

  # standardization ------------------------------------------------
  if (!quiet) {cat("Column-standardizing the design matrix...\n")}

  # centering & scaling
  # NOTE: this C++ call will change the .bk file so that its data are column-standardized
  std_res <- .Call("big_std",
                   X@address,
                   as.integer(bigstatsr::nb_cores()),
                   PACKAGE = "plmmr")
X@address <- std_res$std_X # saves standardized .bk
# TODO: pick up here -- what to do about file names here... the same .bk is being modified..
  if (!quiet) {cat("Standardization completed at", pretty_time())}

  cat("Standardization completed at", pretty_time(), file = outfile, append = TRUE)

  # label return object ------------------------------------------------
  # naming these center and scale values so that I know they relate to the first
  # standardization; there will be another standardization after the rotation
  # in plmm_fit().
  ret <- list(
    std_X = bigmemory::describe(X),
    std_X_n = nrow(X),
    std_X_p = ncol(X),
    std_X_center = std_res$std_X_center,
    std_X_scale = std_res$std_X_scale
  )


  if (!quiet){
    cat("Done with standardization. File formatting in progress\n")
  }

  return(ret)
}
