#' A helper function to standardize a `bigSNP`
#'
#' @param obj           A list that includes:
#'                       (1) subset_X: a `big.matrix` object that has been subset &/or had any additional predictors appended as columns
#'                       (2) ns: a numeric vector indicating the indices of nonsingular columns in subset_X
#' @param prefix        The prefix (as a character string) of the bed/fam data files (e.g., `prefix = 'mydata'`)
#' @param rds_dir       The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`
#' @param non_gen       An integer vector that ranges from 1 to the number of added predictors. Example: if 2 predictors are added, non_gen = 1:2.
#' Note: this is typically passed from the result of `add_predictors()`
#' @param complete_phen Numeric vector with indicesmarking the rows of the original data which have a non-missing entry in the 6th column of the `.fam` file
#' @param id_var        String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID".
#' @param outfile       Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet         Logical: should messages be printed to the console? Defaults to TRUE
#' @param overwrite     Logical: if existing `.bk`/`.rds` files exist for the specified directory/prefix, should these be overwritten?
#'
#' @return A list with a new component of `obj` called 'std_X' - this is an FBM with column-standardized data.
#' List also includes several other indices/meta-data on the standardized matrix
#' @keywords internal
#'
standardize_bigsnp <- function(obj, prefix, rds_dir, non_gen, complete_phen, id_var,
                               outfile, quiet, overwrite){
  # standardization ------------------------------------------------
  if (!quiet) {cat("Column-standardizing the design matrix...\n")}
  # convert FBM pointer into a big.matrix pointer
  std_X <- obj$subset_X |> fbm2bm()

  # centering & scaling
  # NOTE: this C++ call will change the .bk file so that its data are column-standardized
  std_res <- .Call("big_std",
                   std_X@address,
                   as.integer(bigstatsr::nb_cores()),
                   PACKAGE = "plmmr")

  # add timestamp to log -- the standardization step could take a while
  if (!quiet) {cat("Standardization completed at", pretty_time())}

  cat("Standardization completed at", pretty_time(), file = outfile, append = TRUE)

  # label return object ------------------------------------------------
  # naming these center and scale values so that I know they relate to the first
  # standardization; there will be another standardization after the rotation
  # in plmm_fit().
  ret <- list(
    std_X = bigmemory::describe(std_X),
    std_X_n = nrow(std_X),
    std_X_p = ncol(std_X),
    std_X_center = std_res$std_X_center,
    std_X_scale = std_res$std_X_scale,
    ns = obj$ns
  )


  if (!quiet){
    cat("Done with standardization. File formatting in progress\n")
  }

  return(ret)
}
