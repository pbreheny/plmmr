#' A function to impute SNP data
#'
#' @param obj a `bigSNP` object (as created by `read_plink_files()`)
#' @param X A matrix of genotype data as  returned by `name_and_count_bigsnp`
#' @param impute Logical: should data be imputed? Default to TRUE.
#' @param impute_method If 'impute' = TRUE, this argument will specify the kind of imputation desired. Options are:
#'  * mode (default): Imputes the most frequent call. See `bigsnpr::snp_fastImputeSimple()` for details.
#'  * random: Imputes sampling according to allele frequencies.
#'  * mean0: Imputes the rounded mean.
#'  * mean2: Imputes the mean rounded to 2 decimal places.
#'  * xgboost: Imputes using an algorithm based on local XGBoost models. See `bigsnpr::snp_fastImpute()` for details. Note: this can take several minutes, even for a relatively small data set.
#'  @param seed Numeric value to be passed as the seed for `impute_method = 'xgboost'`. Defaults to `as.numeric(Sys.Date())`
#' \@param outfile Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet Logical: should messages be printed to the console? Defaults to TRUE
#' @param ... Optional: additional arguments to `bigsnpr::snp_fastImpute()` (relevant only if impute_method = "xgboost")
#'
#' @return Nothing is returned, but the `obj$genotypes` is overwritten with the imputed version of the data
#' @keywords internal
#'
impute_snp_data <- function(obj, X, impute, impute_method,
                            outfile, quiet, seed = as.numeric(Sys.Date()),...){
  if(!quiet & impute){
    # catch for misspellings
    if(!(impute_method %in% c('mode', 'random', 'mean0', 'mean2', 'xgboost'))){
      stop("\nImpute method is misspecified or misspelled. Please use one of the
           \n5 options listed in the documentation.")
    }
    cat("\nImputing the missing (genotype) values using", impute_method, "method\n")
  }

  if(impute){
    cat("Imputing the missing values using", impute_method, "method\n",
        file = outfile, append = TRUE)

    if(impute_method %in% c('mode', 'random', 'mean0', 'mean2')){

      obj$genotypes <- bigsnpr::snp_fastImputeSimple(Gna = X,
                                                     ncores = bigstatsr::nb_cores(),
                                                     method = impute_method) # dots can pass other args

    } else if (impute_method == "xgboost"){

      obj$genotypes <- bigsnpr::snp_fastImpute(Gna = X,
                                               ncores = bigstatsr::nb_cores(),
                                               infos.chr = obj$chr,
                                               seed = seed,
                                               ...) # dots can pass other args

      cat("\n ***************** NOTE ********************************
          \nAugust 2023: With the xgboost imputation method, there
          \nhave been some issues (particularly on Mac OS) with warnings
          \nthat appear saying 'NA or NaN values in the resulting correlation matrix.'
          \nHowever, we (plmm authors) have
          \n not seen missing values appear in the results -- the imputed data
          \n does not show any NA or NaN values, and models fit on these data run without issue.
          \n We are actively investigating this warning message, and will
          \n make a note in a future release. If using xgboost, proceed with
          \n caution and file an issue if you notice any problems downstream.
          \n ********************************************************")

      # save imputed values (NB: will overwrite obj$genotypes)
      obj$genotypes$code256 <- bigsnpr::CODE_IMPUTE_PRED
    }

    # save the subset data
    obj <- bigsnpr::snp_save(obj)

    cat("Done with imputation.\n",
        file = outfile, append = TRUE)
  }

  return(obj)
}