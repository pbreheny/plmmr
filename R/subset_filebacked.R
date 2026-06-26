#' A helper function to subset `big.matrix` objects
#'
#' @param X                     A filebacked `big.matrix` with the to-be-standardized design matrix
#' @param new_file              Optional user-specified new file for the to-be-created .rds/.bk files.
#' @param complete_samples      Numeric vector with indices marking the rows of the original data which have a non-missing entry in the 6th column of the `.fam` file
#' @param ns                    Numeric vector with the indices of the non-singular columns
#' @param rds_dir               The path to the directory in which you want to create the new `.rds` and `.bk` files. Defaults to `data_dir`
#' @param outfile               Optional: the name (character string) of the logfile to be written.
#' @param quiet                 Logical: should console messages be silenced? Defaults to FALSE
#'
#' @return A list with two components. First, a `big.matrix` object, `subset_X`, representing a design matrix wherein:
#'  * rows are subset to those with complete phenotype information
#'  * columns are subset so that no constant features remain -- this is important for standardization downstream
#'
#' The list also includes the integer vector `ns` which marks which columns of the original matrix were 'non-singular' (i.e. *not* constant features).
#' The `ns` index plays an important role in `plmm_format()` and `untransform()`
#'
#' @keywords internal
#'
subset_filebacked <- function(X, new_file, complete_samples, ns, rds_dir, outfile, quiet) {
  # goal here is to subset the features so that constant features (monomorphic SNPs) are not
  # included in analysis
  # NB: this is also where we remove observations with missing phenotypes, if that was requested
  if (!quiet) {
    cat("Subsetting data to exclude constant features (e.g., monomorphic SNPs)\n")
  }
  cat(
    "Subsetting data to exclude constant features (e.g., monomorphic SNPs)\n",
    file = outfile,
    append = TRUE
  )

  subset_X <- bigmemory::deepcopy(
    x = X,
    row = complete_samples, # filters out samples not represented in both feature data and outcome/predictor data
    col = ns, # filters out singular (constant) features, including any constant predictors
    type = "double", # note object storage type -- this is key...
    backingfile = paste0(new_file, ".bk"),
    backingpath = rds_dir,
    descriptorfile = paste0(new_file, ".desc")
  )

  # save ns indices as part of our object
  list(subset_X = subset_X, ns = ns)
}
