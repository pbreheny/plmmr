#' Preprocess PLINK files using the `bigsnpr` package
#' 
#' @param data_dir The path to the bed/bim/fam data files, *without* a trailing "/" (e.g., use `data_dir = '~/my_dir'`, **not** `data_dir = '~/my_dir/'`)
#' @param prefix The prefix (as a character string) of the bed/fam data files (e.g., `prefix = 'mydata'`)
#' @param impute Logical: should data be imputed? Default to TRUE.
#' @param impute_method If 'impute' = TRUE, this argument will specify the kind of imputation desired. Options are: 
#'  * mode (default): Imputes the most frequent call. See `bigsnpr::snp_fastImputeSimple()` for details. 
#'  * random: Imputes sampling according to allele frequencies.
#'  * mean0: Imputes the rounded mean.
#'  * mean2: Imputes the mean rounded to 2 decimal places.
#'  * xgboost: Imputes using an algorithm based on local XGBoost models. See `bigsnpr::snp_fastImpute()` for details. Note: this can take several minutes, even for a relatively small data set. 
#' @param na_phenotype_vals A vector of numeric values used to code NA values in the phenotype/outcome (this is the 'affection' column in a `bigSNP` object, or the last column of a `.fam` file). Defaults to -9 (matching PLINK conventions).
#' @param id_var String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID". 
#' @param handle_missing_phen A string indicating how missing phenotypes should be handled: 
#'  * "prune" (default): observations with missing phenotype are removed
#'  * "asis": leaves missing phenotypes as NA (this is fine if outcome will be supplied later from a separate file)
#'  * "median": impute missing phenotypes using the median (warning: this is overly simplistic in many cases).
#'  * "mean": impute missing phenotypes using the mean (warning: this is overly simplistic in many cases).
#' @param quiet Logical: should messages be printed to the console? Defaults to TRUE
#' @param gz Logical: are the bed/bim/fam files g-zipped? Defaults to FALSE. NOTE: if TRUE, process_plink will unzip your zipped files.
#' @param outfile Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param overwrite Logical: if existing `.bk`/`.rds` files exist for the specified directory/prefix, should these be overwritten? Defaults to FALSE. Set to TRUE if you want to change the imputation method you're using, etc. 
#' @param ... Optional: additional arguments to `bigsnpr::snp_fastImpute()` (relevant only if impute_method = "xgboost")
#' @param add_predictor_fam Optional: if you want to include "sex" (the 5th column of `.fam` file) in the analysis, specify 'sex' here.
#' @param add_predictor_ext Optional: add additional covariates/predictors/features from an external file (i.e., not a PLINK file). 
#' This argument takes one of two kinds of arguments: 
#'  (a) a **named** numeric vector, where the names align with the sample IDs in the PLINK files. 
#' The names will be used to subset and align this external covariate with the supplied PLINK data.
#' 
#'  (b) a numeric matrix whose row names align with the sample IDs in the PLINK files. 
#'  The names will be used to subset and align this external covariate with the supplied PLINK data.
#' @returns Nothing is returned by this function; instead, files 'prefix.rds',
#'  'prefix.bk', and std_prefix.bk are created in the location specified by data_dir. Note that this 
#'  this function need only be run once; in subsequent data analysis/scripts, 
#'  `get_data()` will access the '.rds' file. 
#'    
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' process_plink(data_dir = plink_example(parent = T),
#'   prefix = "penncath_lite",
#'   gz = TRUE,
#'   outfile = "process_penncath",
#'   overwrite = TRUE,
#'   impute_method = "mode")
#' }
#' 
process_plink <- function(data_dir,
                          prefix,
                          impute = TRUE,
                          impute_method = 'mode',
                          na_phenotype_vals = c(-9),
                          id_var = "IID",
                          handle_missing_phen = "prune",
                          quiet = FALSE,
                          gz = FALSE,
                          outfile,
                          overwrite = FALSE,
                          add_predictor_fam = NULL,
                          add_predictor_ext = NULL,
                          ...){
  
  # start log ------------------------------------------
  if(missing(outfile)){
    outfile = paste0(data_dir, "/process_plink.log")
  } else {
    outfile = paste0(outfile, ".log")
  }
  log_con <- file(outfile)
  cat("### Processing PLINK files for PLMM ###", file = log_con)
  cat("\nLogging to ", outfile, file = outfile, append = TRUE)
  cat("\nPreprocessing", prefix, "data:", file = outfile, append = TRUE)
  
  if(!quiet){
    cat("\nLogging to", outfile)
    cat("\nPreprocessing", prefix, "data:")
  }

  # read in PLINK files --------------------------------
  step1_obj <- read_plink_files(data_dir, prefix, gz, outfile, overwrite, quiet)

  # name and count ------------------------------------
  step2 <- name_and_count_bigsnp(step1_obj, id_var, quiet)

  # chromosome check ---------------------------------
  # only consider SNPs on chromosomes 1-22
  if(step2$chr_range[1] < 1 | step2$chr_range[2] > 22){
    stop("\nplmmr only analyzes autosomes -- please remove variants on chromosomes outside 1-22.
         This can be done in PLINK 1.9; see the documentation in https://www.cog-genomics.org/plink/1.9/filter#chr")
  }
  
  # notify about missing genotypes & phenotypes ------------------
  step3 <- handle_missingness(obj = step2$obj, counts = step2$counts, X = step2$X,
                              na_phenotype_vals = na_phenotype_vals,
                              handle_missing_phen = handle_missing_phen,
                              outfile = outfile, quiet = quiet)

  # imputation -------------------------------------
  step4_obj <- impute_snp_data(step2$obj, step2$X, impute, impute_method,
                               outfile, quiet,...)

  # add predictors from external files ----------------
  step5 <- add_predictors(step4_obj, add_predictor_fam, add_predictor_ext,
                          id_var, step2$og_plink_ids, quiet)

  # subsetting -----------------------------------------
  step6_obj <- subset_fbm(step5$obj, step2$counts, handle_missing_phen,
                      step3$complete_phen, step5$non_gen, data_dir, prefix, outfile, quiet)

  # standardization ------------------------------------------------
  step7 <- standardize_fbm(step6_obj, prefix, step5$non_gen, step3$complete_phen,
                  id_var, outfile, quiet)
  
  if(!quiet){cat("\nDone with standardization. Processed files now saved as .rds object.")}
  close(log_con)
}