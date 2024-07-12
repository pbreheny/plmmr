#' A helper function to handle missing genotype values from PLINK data
#' @param obj     a `bigSNP` object
#' @param na_phenotype_vals A vector of numeric values used to code NA values in the phenotype/outcome (this is the 'affection' column in a `bigSNP` object, or the last column of a `.fam` file). Defaults to -9 (matching PLINK conventions).
#' @param handle_missing_phen A string indicating how missing phenotypes should be handled:
#'  * "prune" (default): observations with missing phenotype are removed
#'  * "asis": leaves missing phenotypes as NA (this is fine if outcome will be supplied later from a separate file
#' @param outfile Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet   Logical: should messages be printed to the console? Defaults to TRUE
#'
#' @return a list of two components:
#' * na_idx: a logical vector indicating which columns have missing values
#' * prop_na: the proportion of missingness in each of the columns tagged by 'na_idx'
#' @keywords internal
#'
handle_missingness <- function(obj, na_phenotype_vals,
                               handle_missing_phen, outfile, quiet){
  # handle missing phenotypes ---------------------------------------
  # make missing phenotypes explicit (need both of the following because
  # bigstatsr::big_copy() does not handle negative indices)
  is_phen_missing <- obj$fam[,6] %in% na_phenotype_vals
  complete_phen <- which(!(is_phen_missing))
  na_phen <- which(is_phen_missing)
  names(na_phen) <- obj$fam$sample.ID[na_phen]

  if (handle_missing_phen == 'prune'){
    if(!quiet){
      cat("\nWill prune out", length(na_phen), "samples/observations with missing phenotype data.")
      # Note: the actual pruning happens in the 'subset' step
    }
    cat("\nWill prune out", length(na_phen), "samples/observations with missing phenotype data\n",
        file = outfile, append = TRUE)

  } else if (handle_missing_phen == 'asis'){
    if(!quiet){
      cat("\nWill mark", length(na_phen), "samples/observations as having missing phenotype data.")
    }
    cat("\nWill mark", length(na_phen), "samples/observations as having missing phenotype data\n",
        file = outfile, append = TRUE)

    obj$fam$affection[na_phen] <- NA_integer_

  }

  return(list(obj = obj,
              complete_phen = complete_phen,
              na_phen = na_phen))

}