#' A helper function to handle missing genotype values from PLINK data
#' @param obj     a list with elements 'geno_plus_predictors' (a `big.matrix` object), fam (the fam file data), and ...
#' @param na_phenotype_vals A vector of numeric values used to code NA values in the phenotype/outcome (this is the 'affection' column in a `bigSNP` object, or the last column of a `.fam` file). Defaults to -9 (matching PLINK conventions).
#' @param handle_missing_phen A string indicating how missing phenotypes should be handled:
#'  * "prune" (default): observations with missing phenotype are removed
#'  * "asis": leaves missing phenotypes as NA (this is fine if outcome will be supplied later from a separate file
#' @param outfile Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet   Logical: should messages be printed to the console? Defaults to FALSE (which leaves the print messages on...)
#'
#' @return a list of two components:
#' * na_idx: a logical vector indicating which columns have missing values
#' * prop_na: the proportion of missingness in each of the columns tagged by 'na_idx'
#' @keywords internal
#'
handle_missingness <- function(obj, na_phenotype_vals,
                               handle_missing_phen, outfile, quiet){

  # check for constant features in genotypes also:
  ns_genotypes <- count_constant_features(fbm = obj$X,
                                          ind.row = complete_phen,
                                          outfile = outfile,
                                          quiet = quiet)

  return(list(complete_phen = complete_phen,
              na_phen = na_phen,
              ns_genotypes = ns_genotypes))

}