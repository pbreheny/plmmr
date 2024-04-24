#' A helper function to handle missing genotype values from PLINK data 
#' @param obj     a `bigSNP` object 
#' @param counts  A dataframe as created by `bigstatsr::big_counts()` and returned by `name_and_count_bigsnp`
#' @param X       A matrix of genotype data as  returned by `name_and_count_bigsnp`
#' @param na_phenotype_vals A vector of numeric values used to code NA values in the phenotype/outcome (this is the 'affection' column in a `bigSNP` object, or the last column of a `.fam` file). Defaults to -9 (matching PLINK conventions).
#' @param handle_missing_phen A string indicating how missing phenotypes should be handled: 
#'  * "prune" (default): observations with missing phenotype are removed
#'  * "asis": leaves missing phenotypes as NA (this is fine if outcome will be supplied later from a separate file)
#'  * "median": impute missing phenotypes using the median (warning: this is overly simplistic in many cases).
#'  * "mean": impute missing phenotypes using the mean (warning: this is overly simplistic in many cases).
#' @param outfile Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet   Logical: should messages be printed to the console? Defaults to TRUE
#'
#' @return a list of two components: 
#' * na_idx: a logical vector indicating which columns have missing values
#' * prop_na: the proportion of missingness in each of the columns tagged by 'na_idx'
#' @keywords internal
#'
handle_missingness <- function(obj, X, counts, na_phenotype_vals, 
                               handle_missing_phen, outfile, quiet){
  # notify about missing (genotype) values ------------------------------
  na_idx <- counts[4,] > 0
  prop_na <- counts[4,]/nrow(X)
  
  cat("\nThere are a total of ", sum(na_idx), "SNPs with missing values",
      file = outfile, append = TRUE)
  cat("\nOf these, ", sum(prop_na > 0.5),
      " are missing in at least 50% of the samples",
      file = outfile, append = TRUE)
  if(!quiet){
    cat("\nThere are a total of ", sum(na_idx), "SNPs with missing values")
    cat("\nOf these, ", sum(prop_na > 0.5), " are missing in at least 50% of the samples")
  }
  
  # handle missing phenotypes ---------------------------------------
  # make missing phenotypes explicit (need both of the following because 
  # bigstatsr::big_cop() does not handle negative indices)
  complete_phen <- which(!(obj$fam$affection %in% na_phenotype_vals))
  na_phen <- which(obj$fam$affection %in% na_phenotype_vals)
  names(na_phen) <- obj$fam$sample.ID[na_phen]
  
  if (handle_missing_phen == 'prune'){
    if(!quiet){
      cat("\nWill prune out ", length(na_phen), " samples/observations with missing phenotype data.")
      # Note: the actual pruning happens in the 'subset' step 
    }
    
  } else if (handle_missing_phen == 'asis'){
    if(!quiet){
      cat("\nWill mark ", length(na_phen), " samples/observations as having missing phenotype data.")
    }
    obj$fam$affection[na_phen] <- NA_integer_
  } else {
    if(!quiet){
      cat("\nImputing phenotype data for ", length(na_phen), " samples/observations.")
    }
    obj$fam$affection[na_phen] <- switch(handle_missing_phen,
                                         median = median(obj$fam$affection[complete_phen]),
                                         mean = mean(obj$fam$affection[complete_phen]))
  }
  
  return(list(obj = obj, 
              na_idx = na_idx,
              prop_na = prop_na,
              complete_phen = complete_phen,
              na_phen = na_phen))
  
}