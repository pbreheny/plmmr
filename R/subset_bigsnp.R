#' A helper function to subset `bigSNP` objects
#'
#' @param X                   A filebacked `big.matrix` with the to-be-standardized design matrix
#' @param new_file                Optional user-specified new_file for the to-be-created .rds/.bk files.
#' @param complete_outcome         Numeric vector with indicesmarking the rows of the original data which have a non-missing entry in the 6th column of the `.fam` file
#' @param non_gen               an integer vector that ranges from 1 to the number of added predictors. Example: if 2 predictors are added, non_gen = 1:2.
#'                              **Note**: this is typically passed from the result of `add_predictors()`
#' @param ns_genotypes          Numeric vector with the indices of the non-singular genotype (feature) columns
#'                              This vector is created in `handle_missingness()`
#' @param rds_dir               The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`
#' @param outfile               Optional: the name (character string) of the new_file of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet                 Logical: should messages be printed to the console? Defaults to FALSE (which leaves the print messages on...)
#'
#' @return A `bigSNP` object that includes a new component, 'subset_X', representing a design matrix wherein:
#' *  rows are subset according to user's specification in `handle_missing_phen`
#' *  columns are subset so that no constant features remain -- this is important for standardization downstream
#' The new 'obj' also includes the integer vector 'ns' which marks which columns of the original matrix were 'non-singular' (i.e. *not* constant features).
#' The 'ns' index plays an important role in `plmm_format()` and `untransform()` (both helper functions in model fitting)
#' @keywords internal
#'
subset_bigsnp <- function(X, complete_outcome, non_gen,
                          ns_genotypes, rds_dir, new_file, outfile, quiet){

  # goal here is to subset the features so that constant features (monomorphic SNPs) are not
  # included in analysis
  # NB: this is also where we remove observations with missing phenotypes, if that was requested
  if (!quiet){
    cat("Subsetting data to exclude constant features (e.g., monomorphic SNPs)\n")
  }
  cat("Subsetting data to exclude constant features (e.g., monomorphic SNPs)\n",
      file = outfile, append = TRUE)

    if (length(non_gen) > 0) {
      # TODO: adjust this option for more general case -- the issue here is not PLINK
      #   format, just the adjustment of the ns indices due to the additional covariates
      ns <- c(non_gen, ns_genotypes + length(non_gen))
      # Note: add_predictors() already scanned for constant features among the added predictors
    }

  subset_X <- bigmemory::deepcopy(x = X,
                                  row = complete_outcome, # filters out rows with missing phenotypes
                                  col = ns,
                                  type = "double", # this is key...
                                  backingfile = paste0(new_file, ".bk"),
                                  backingpath = rds_dir,
                                  descriptorfile = paste0(new_file, '.desc'))

  # save ns indices as part of our object
  return(list(subset_X = subset_X, ns = ns))

}
