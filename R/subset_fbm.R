#' A helper function to subset `FBM` objects 
#'
#' @param obj 
#' @param counts 
#' @param handle_missing_phen 
#' @param complete_phen
#' @param non_gen 
#' @param data_dir              The path to the bed/bim/fam data files, *without* a trailing "/" (e.g., use `data_dir = '~/my_dir'`, **not** `data_dir = '~/my_dir/'`)
#' @param prefix 
#' @param outfile               Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet 
#'
#' @return
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
  
  sub_bk_extension <- paste0("subset_", prefix) 

  if (handle_missing_phen == "prune"){
    if ("geno_plus_predictors" %in% names(obj)) {
      new_counts <- bigstatsr::big_counts(obj$genotypes,
                                          ind.row = complete_phen)
      ns_genotypes <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      ns <- c(non_gen, ns_genotypes + length(non_gen)) # Note: add_predictors() already scanned for constant features among the added predictors 
      obj$subset_X <- bigstatsr::big_copy(obj$geno_plus_predictors,
                                          ind.row = complete_phen, # filters out rows with missing phenotypes
                                          ind.col = ns,
                                          backingfile = paste0(data_dir,"/", sub_bk_extension))
    } else {
      new_counts <- bigstatsr::big_counts(obj$genotypes,
                                          ind.row = complete_phen)
      ns <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      obj$subset_X <- bigstatsr::big_copy(obj$genotypes,
                                          ind.row = complete_phen, # filters out rows with missing phenotypes
                                          ind.col = ns,
                                          backingfile = paste0(data_dir,"/", sub_bk_extension))
    }
    
  } else {
    if ("geno_plus_predictors" %in% names(obj)) {
      new_counts <- bigstatsr::big_counts(obj$genotypes)
      ns_genotypes <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      ns <- c(non_gen, ns_genotypes + length(non_gen)) # Note: add_predictors() already scanned for constant features among the added predictors 
      obj$subset_X <- bigstatsr::big_copy(obj$geno_plus_predictors,
                                          ind.col = ns,
                                          backingfile = paste0(data_dir,"/", sub_bk_extension))
    } else {
      new_counts <- bigstatsr::big_counts(obj$genotypes)
      ns <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      obj$subset_X <- bigstatsr::big_copy(obj$genotypes,
                                          ind.col = ns,
                                          backingfile = paste0(data_dir,"/", sub_bk_extension))
    }
  }
  
  # save ns indices as part of our object 
  obj$ns <- ns
  
  return(obj)
}