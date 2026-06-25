#' A function to impute SNP data
#'
#' @param obj A `bigSNP` object (as created by `read_plink_files()`)
#' @param X A matrix of genotype data as returned by `name_and_count_bigsnp()`
#' @param chr A numeric vector of chromosomal locations of the SNPs.
#' @param impute Logical: should data be imputed? Defaults to TRUE.
#' @param impute_method If `impute = TRUE`, this argument will specify the kind of imputation desired. Options are:
#'  * `mode` (default): Imputes the most frequent call. See `bigsnpr::snp_fastImputeSimple()` for details.
#'  * `random`: Imputes sampling according to allele frequencies.
#'  * `mean0`: Imputes the rounded mean.
#'  * `mean2`: Imputes the mean rounded to 2 decimal places.
#'  * `xgboost`: Imputes using an algorithm based on local XGBoost models. See `bigsnpr::snp_fastImpute()` for details. Note: this can take several minutes, even for a relatively small data set.
#' @param parallel Logical: should the computations within this function be run in parallel? Defaults to TRUE. See `count_cores()` and `?bigparallelr::assert_cores` for more details.
#'                  In particular, the user should be aware that too much parallelization can make computations *slower*.
#' @param outfile Optional: the name (character string) of the prefix of the logfile to be written.
#' @param quiet Logical: should console messages be silenced? Defaults to FALSE
#' @param seed Numeric value to be passed as the seed for `impute_method = 'xgboost'`. Defaults to `as.numeric(Sys.Date())`
#' @param ... Optional: additional arguments to `bigsnpr::snp_fastImpute()` (relevant only if `impute_method = 'xgboost'`)
#'
#' @return Nothing is returned, but the `obj$genotypes` is overwritten with the imputed version of the data
#'
#' @keywords internal
#'
impute_snp_data <- function(
  obj,
  X,
  chr,
  impute,
  impute_method,
  parallel,
  outfile,
  quiet,
  seed = as.numeric(Sys.Date()),
  ...
) {
  if (!quiet && impute) {
    # catch for misspellings
    if (
      !(impute_method %in% c("mode", "random", "mean0", "mean2", "xgboost"))
    ) {
      stop(
        "\nImpute method is misspecified or misspelled. Please use one of the
           \n5 options listed in the documentation."
      )
    }
    cat(
      "Imputing the missing (genotype) values using",
      impute_method,
      "method...\n"
    )
  }

  if (impute) {
    cat(
      "Imputing the missing values using",
      impute_method,
      "method...\n",
      file = outfile,
      append = TRUE
    )

    if (impute_method %in% c("mode", "random", "mean0", "mean2")) {
      if (parallel) {
        obj$genotypes <- bigsnpr::snp_fastImputeSimple(
          Gna = X,
          ncores = count_cores(),
          method = impute_method
        )
      } else {
        obj$genotypes <- bigsnpr::snp_fastImputeSimple(
          Gna = X,
          method = impute_method
        )
      }
    } else if (impute_method == "xgboost") {
      withCallingHandlers(
        {
          if (parallel) {
            bigsnpr::snp_fastImpute(
              Gna = X,
              ncores = count_cores(),
              infos.chr = chr,
              seed = seed,
              ...
            ) # dots can pass other args
          } else {
            bigsnpr::snp_fastImpute(Gna = X, infos.chr = chr, seed = seed, ...) # dots can pass other args
          }
        },
        warning = function(w) {
          if (!.plmmr_env$xgboost_warning_shown) {
            cat(
              "\n*************************** NOTE ***************************\n",
              "With the xgboost imputation method, you may receive the warning",
              "\n'NA or NaN values in the resulting correlation matrix.' This",
              "\nis thrown by the bigsnpr package and likely results from several",
              "\nof your SNPs having extremely small minor allele counts. These",
              "\nvalues are excluded by the imputation algorithm, but if you would",
              "\nlike to remedy the warnings, it is easiest to filter your data by",
              "\nMAC (or MAF) using PLINK. For example:\n",
              "\n",
              "plink --bfile penncath_lite --mac 20 --make-bed --out penncath_filt",
              "\n***************************************************************"
            )
            .plmmr_env$xgboost_warning_shown <- TRUE
          }
        }
      )

      # save imputed values (NB: will overwrite obj$genotypes)
      X$code256 <- bigsnpr::CODE_IMPUTE_PRED
    }

    # save the imputed data
    obj <- bigsnpr::snp_save(obj)

    if (!quiet) {
      cat("Done with imputation.\n")
    }

    cat("Done with imputation.\n", file = outfile, append = TRUE)
  }

  obj
}
