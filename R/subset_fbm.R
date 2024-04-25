#' A helper function to subset `FBM` objects 
#'
#' @param obj                   A `bigSNP` object 
#' @param counts                A dataframe as returned by `bigstatsr::big_counts()`
#' @param handle_missing_phen   A string indicating how missing phenotypes should be handled: 
#'  * "prune" (default): observations with missing phenotype are removed
#'  * "asis": leaves missing phenotypes as NA (this is fine if outcome will be supplied later from a separate file)
#'  * "median": impute missing phenotypes using the median (warning: this is overly simplistic in many cases).
#'  * "mean": impute missing phenotypes using the mean (warning: this is overly simplistic in many cases). 
#' @param complete_phen         Numeric vector with indicies marking the rows of the original data which have a non-missing entry in the 6th column of the `.fam` file 
#' @param non_gen               an integer vector that ranges from 1 to the number of added predictors. Example: if 2 predictors are added, non_gen = 1:2. 
#' Note: this is typically passed from the result of `add_predictors()`
#' @param data_dir              The path to the bed/bim/fam data files, *without* a trailing "/" (e.g., use `data_dir = '~/my_dir'`, **not** `data_dir = '~/my_dir/'`)
#' @param prefix                The prefix (as a character string) of the bed/fam data files (e.g., `prefix = 'mydata'`)
#' @param outfile               Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet                 Logical: should messages be printed to the console? Defaults to TRUE
#'
#' @return A `bigSNP` object that includes a new component, 'subset_X', representing a design matrix wherein:
#' *  rows are subset according to user's specification in `handle_missing_phen`
#' *  columns are subset so that no constant features remain -- this is important for standardization downstream
#' The new 'obj' also includes the integer vector 'ns' which marks which columns of the original matrix were 'non-singular' (i.e. *not* constant features).
#' The 'ns' index plays an important role in `plmm_format()` and `untransform()` (both helper functions in model fitting)
#' @keywords internal
#'
subset_fbm <- function(obj, counts, handle_missing_phen, complete_phen, non_gen,
                       data_dir, prefix, outfile, quiet){
  # goal here is to subset the features so that constant features (monomorphic SNPs) are not 
  # included in analysis
  # NB: this is also where we remove observations with missing phenotypes, if that was requested
  
  if (!quiet){
    cat("\nSubsetting data to exclude constant features (e.g., monomorphic SNPs)",
        file = outfile, append = TRUE)
  }
  
  # sub_bk_extension <- paste0("subset_", prefix) 

  if (handle_missing_phen == "prune"){
    if ("geno_plus_predictors" %in% names(obj)) {
      new_counts <- bigstatsr::big_counts(obj$genotypes,
                                          ind.row = complete_phen)
      ns_genotypes <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      ns <- c(non_gen, ns_genotypes + length(non_gen)) # Note: add_predictors() already scanned for constant features among the added predictors 
      obj$subset_X <- bigstatsr::big_copy(obj$geno_plus_predictors,
                                          # backingfile = paste0(data_dir,"/", sub_bk_extension),
                                          ind.row = complete_phen, # filters out rows with missing phenotypes
                                          ind.col = ns)
    } else {
      new_counts <- bigstatsr::big_counts(obj$genotypes,
                                          ind.row = complete_phen)
      ns <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      obj$subset_X <- bigstatsr::big_copy(obj$genotypes,
                                          # backingfile = paste0(data_dir,"/", sub_bk_extension),
                                          ind.row = complete_phen, # filters out rows with missing phenotypes
                                          ind.col = ns)
    }
    
  } else {
    if ("geno_plus_predictors" %in% names(obj)) {
      new_counts <- bigstatsr::big_counts(obj$genotypes)
      ns_genotypes <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      ns <- c(non_gen, ns_genotypes + length(non_gen)) # Note: add_predictors() already scanned for constant features among the added predictors 
      obj$subset_X <- bigstatsr::big_copy(obj$geno_plus_predictors,
                                          # backingfile = paste0(data_dir,"/", sub_bk_extension),
                                          ind.col = ns)
    } else {
      new_counts <- bigstatsr::big_counts(obj$genotypes)
      ns <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      obj$subset_X <- bigstatsr::big_copy(obj$genotypes,
                                          # backingfile = paste0(data_dir,"/", sub_bk_extension),
                                          ind.col = ns)
    }
  }
  
  # save ns indices as part of our object 
  obj$ns <- ns
  
  return(obj)
}