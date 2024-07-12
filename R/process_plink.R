#' Preprocess PLINK files using the `bigsnpr` package
#'
#' @param data_dir              The path to the bed/bim/fam data files, *without* a trailing "/" (e.g., use `data_dir = '~/my_dir'`, **not** `data_dir = '~/my_dir/'`)
#' @param rds_dir               The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`
#' @param prefix                The prefix (as a character string) of the bed/fam data files (e.g., `prefix = 'mydata'`)
#' @param outfile               Optional: the name (character string) of the prefix of the logfile to be written.
#'                              Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile, created in the same directory as 'data_dir'.
#' @param impute                Logical: should data be imputed? Default to TRUE.
#' @param impute_method         If 'impute' = TRUE, this argument will specify the kind of imputation desired. Options are:
#'                                * mode (default): Imputes the most frequent call. See `bigsnpr::snp_fastImputeSimple()` for details.
#'                                * random: Imputes sampling according to allele frequencies.
#'                                * mean0: Imputes the rounded mean.
#'                                * mean2: Imputes the mean rounded to 2 decimal places.
#'                                * xgboost: Imputes using an algorithm based on local XGBoost models. See `bigsnpr::snp_fastImpute()` for details. Note: this can take several minutes, even for a relatively small data set.
#' @param na_phenotype_vals   A vector of numeric values used to code NA values in the phenotype/outcome (this is the 'affection' column in a `bigSNP` object, or the last column of a `.fam` file). Defaults to -9 (matching PLINK conventions).
#' @param id_var              String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID".
#' @param add_phen            Optional: A **data frame** with at least two columns: and ID column and a phenotype column
#' @param pheno_id            Optional: A string specifying the name of the ID column in `pheno`. MUST be specified if `add_phen` is specified.
#' @param pheno_name           Optional: A string specifying the name of the phenotype column in `pheno`.  MUST be specified if `add_phen` is specified. This column will be used as the default `y` argument to 'plmm()'.
#' @param handle_missing_phen A string indicating how missing phenotypes should be handled:
#'                                * "prune" (default): observations with missing phenotype are removed
#'                                * "median": impute missing phenotypes using the median (warning: this is overly simplistic in many cases).
#'                                * "mean": impute missing phenotypes using the mean (warning: this is overly simplistic in many cases).
#'                            Note: for data coming from PLINK, no missing values of the phenotype are allowed. You have to (1) supply phenotype from an external file,
#'                            (2) prune missing values, or (3) impute missing values.
#' @param quiet               Logical: should messages to be printed to the console be silenced? Defaults to FALSE
#' @param overwrite           Logical: if existing `.bk`/`.rds` files exist for the specified directory/prefix, should these be overwritten? Defaults to FALSE. Set to TRUE if you want to change the imputation method you're using, etc.
#'                            **Note**: If there are multiple `.rds` files with names that start with "std_prefix_...", **this will error out**.
#'                            To protect users from accidentally deleting files with saved results, only one `.rds` file can be removed with this option.
#' @param ...                 Optional: additional arguments to `bigsnpr::snp_fastImpute()` (relevant only if impute_method = "xgboost")
#'
#' @return Nothing is returned by this function, but (at least) two files are created in
#' the location specified by `rds_dir`:
#'
#' * 'std_prefix.rds': This is the `bigsnpr::bigSNP` object
#' that holds the PLINK data along with meta-data. See details for explanation of what
#' is included in this meta-data
#'
#' * 'std_prefix.bk': Created by the call to `standardize_fbm()`, this is the
#' backingfile that stores the numeric data of the standardized design matrix `std_X`
#'
#'  Note that `process_plink()` need only be run once for a given set of PLINK
#'  files; in subsequent data analysis/scripts, `get_data()` will access the '.rds' file.
#'
#' @export
#'
#' @details
#' The '.rds' object created by this function has the following elements:
#' * std_X: The file-backed design matrix, as an `FBM` object (see `bigstatsr` package for more info)
#' * fam: Data frame equivalent of PLINK '.bim' file
#' * map: Data frame equivalent of PLINK '.fam' file
#' * colnames: Character vector of column names for the original data (includes constant features, i.e. monomorphic SNPs)
#' * rownames: Character vector of row names for the original data.
#' * n: Number of rows (samples) in the original data
#' * p: Number of columns (features, SNPs, markers, ...) in the original data
#' * ns: Numeric vector of indices marking the non-singular columns of the original data
#' * std_X_center: Numeric vector of values used to center the non-singular columns of the data
#' * std_X_scale: Numeric vector of values used to scale the non-singular columns of the data
#' * std_X_colnames: Character vector of column names for the standardized data
#' * std_X_rownames: Character vector of row names for the standardized data.
#' * complete_phen: Numeric vector of indices marking the samples with a non-missing phenotype. Ony applicable if `handle_missing_phen = 'prune'`
#' * id_var: String specifying which ID column in the '.fam' file had the unique sample ID: 'FID' (1st column) or 'IID' (2nd column)
#'
#' @examples
#' \donttest{
#' temp_dir <- paste0(tempdir()) # using a temporary directory here
#' process_plink(data_dir = find_example_data(parent = TRUE), # reads data that ships with plmmr
#'               rds_dir = temp_dir,
#'               prefix = "penncath_lite",
#'               outfile = "process_penncath",
#'               overwrite = TRUE,
#'               impute_method = "mode")
#'
#'  # pen <- get_data(file.path(temp_dir, "std_penncath_lite"))
#'  # str(pen)
#'  }
#'
process_plink <- function(data_dir,
                          rds_dir = data_dir,
                          prefix,
                          outfile,
                          impute = TRUE,
                          impute_method = 'mode',
                          na_phenotype_vals = c(-9),
                          id_var = "IID",
                          handle_missing_phen = "prune",
                          add_phen = NULL,
                          pheno_id = NULL,
                          pheno_name = NULL,
                          quiet = FALSE,
                          overwrite = FALSE,
                          ...){

  # start log -----------------------------------------
  if(missing(outfile)){
    outfile = file.path(data_dir, "process_plink")
  }

  logfile <- create_log(outfile = outfile)

  if(!quiet){
    cat("\nLogging to", logfile)
    cat("\nPreprocessing", prefix, "data:")
  }
  cat("\nPreprocessing", prefix, "data\n", file = logfile, append = TRUE)


  # read in PLINK files --------------------------------
  plink_obj <- read_plink_files(data_dir = data_dir,
                                prefix = prefix,
                                rds_dir = rds_dir,
                                outfile = logfile,
                                overwrite = overwrite,
                                quiet = quiet)

  # add external phenotype, if needed -----------------
  if (id_var == "IID"){
    geno_id <- "sample.ID"
  } else if (id_var == "FID"){
    geno_id <- "family.ID"
  } else {
    stop("\nThe argument to id_var is misspecified. Must be one of 'IID' or 'FID', and
         the corresponding variable *must* be of type 'char'.")
  }

  if (!is.null(add_phen)){
    step1 <- add_external_phenotype(geno = plink_obj,
                                    geno_id = geno_id,
                                    pheno = add_phen,
                                    pheno_id = pheno_id,
                                    pheno_col = pheno_name,
                                    outfile = logfile)
    handle_missing_phen <- 'asis'

  } else {
    step1 <- plink_obj

  }

  # name and count ------------------------------------
  step2 <- name_and_count_bigsnp(obj = step1,
                                 id_var = id_var,
                                 quiet = quiet,
                                 outfile = logfile)

  # chromosome check ---------------------------------
  # only consider SNPs on chromosomes 1-22
  if(step2$chr_range[1] < 1 | step2$chr_range[2] > 22){
    stop("plmmr only analyzes autosomes -- please remove variants on
         chromosomes outside 1-22.
         This can be done in PLINK 1.9; see the documentation in
         https://www.cog-genomics.org/plink/1.9/filter#chr \n")
  }

  # notify about missing genotypes & phenotypes ---------------------------------
  step3 <- handle_missingness(obj = step2$obj,
                              counts = step2$counts,
                              X = step2$X,
                              na_phenotype_vals = na_phenotype_vals,
                              handle_missing_phen = handle_missing_phen,
                              outfile = logfile,
                              quiet = quiet)

  # imputation ------------------------------------------------------------------
  step4 <- impute_snp_data(step2$obj,
                           step2$X,
                           impute = impute,
                           impute_method = impute_method,
                           outfile = logfile,
                           quiet = quiet, ...)

  step4$complete_phen <- step3$complete_phen
  step4$handle_missing_phen <- handle_missing_phen
  bigsnpr::snp_save(step4)
  # cleanup --------------------------------------------------------------------
  # These steps remove intermediate rds/bk files created by the steps of the data management process
  list.files(rds_dir, pattern=paste0('^file.*.bk'), full.names=TRUE) |>
    file.remove()
  gc()

  if(!quiet){cat("\nprocess_plink() completed \nProcessed files now saved as .rds object.")}

  cat("\nprocess_plink() completed. \nProcessed files now saved as",
      file.path(rds_dir, paste0(prefix, ".rds")),
      "at", pretty_time(),
      file = logfile, append = TRUE)


  return(file.path(data_dir, paste0(prefix, ".rds")))

}



