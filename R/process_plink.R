#' Preprocess PLINK files using the `bigsnpr` package
#'
#' @param data_dir              The path to the bed/bim/fam data files, *without* a trailing "/" (e.g., use `data_dir = '~/my_dir'`, **not** `data_dir = '~/my_dir/'`)
#' @param data_prefix           The prefix (as a character string) of the bed/fam data files (e.g., `data_prefix = 'mydata'`)
#' @param rds_dir               The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`
#' @param rds_prefix            String specifying the user's preferred filename for the to-be-created .rds file (will be create insie `rds_dir` folder)
#'                              Note: 'rds_prefix' cannot be the same as 'data_prefix'
#' @param logfile               Optional: the name (character string) of the prefix of the logfile to be written in 'rds_dir'. Default to NULL (no log file written).
#'                              Note: if you supply a file path in this argument, it will error out with a "file not found" error. Only supply the string; e.g., if you want my_log.log, supply 'my_log', the my_log.log file will appear in rds_dir.
#' @param impute                Logical: should data be imputed? Default to TRUE.
#' @param impute_method         If 'impute' = TRUE, this argument will specify the kind of imputation desired. Options are:
#'                                * mode (default): Imputes the most frequent call. See `bigsnpr::snp_fastImputeSimple()` for details.
#'                                * random: Imputes sampling according to allele frequencies.
#'                                * mean0: Imputes the rounded mean.
#'                                * mean2: Imputes the mean rounded to 2 decimal places.
#'                                * xgboost: Imputes using an algorithm based on local XGBoost models. See `bigsnpr::snp_fastImpute()` for details. Note: this can take several minutes, even for a relatively small data set.

#' @param id_var              String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID"
#' @param quiet               Logical: should messages to be printed to the console be silenced? Defaults to FALSE
#' @param overwrite           Logical: if existing `.bk`/`.rds` files exist for the specified directory/prefix, should these be overwritten? Defaults to FALSE. Set to TRUE if you want to change the imputation method you're using, etc.
#'                            **Note**: If there are multiple `.rds` files with names that start with "std_prefix_...", **this will error out**.
#'                            To protect users from accidentally deleting files with saved results, only one `.rds` file can be removed with this option.
#' @param ...                 Optional: additional arguments to `bigsnpr::snp_fastImpute()` (relevant only if impute_method = "xgboost")
#'
#' @return The filepath to the '.rds' object created; see details for explanation.
#'
#' @details Three files are created in the location specified by `rds_dir`:
#'
#' * 'rds_prefix.rds': This is a list with three items:
#'    (1) `X`: the filebacked `bigmemory::big.matrix` object pointing to the imputed genotype data.
#'    This matrix has type 'double', which is important for downstream operations in `create_design()`
#'    (2) `map`: a data.frame with the PLINK 'bim' data (i.e., the variant information)
#'    (3) `fam`: a data.frame with the PLINK 'fam' data (i.e., the pedigree information)
#'
#' * 'prefix.bk': This is the
#' backingfile that stores the numeric data of the genotype matrix
#'
#' * 'rds_prefix.desc'" This is the description file, as needed by the
#'
#'  Note that `process_plink()` need only be run once for a given set of PLINK
#'  files; in subsequent data analysis/scripts, `get_data()` will access the '.rds' file.
#'
#' @export
#'
#' @details See vignette on processing PLINK files
#' @examples
#' \donttest{
#' temp_dir <- paste0(tempdir()) # using a temporary directory here
#' process_plink(data_dir = find_example_data(parent = TRUE), # reads data that ships with plmmr
#'               rds_dir = temp_dir,
#'               data_prefix = "penncath_lite",
#'               logfile = "process_penncath",
#'               overwrite = TRUE,
#'               impute_method = "mode")
#'
#'  # pen <- get_data(file.path(temp_dir, "std_penncath_lite"))
#'  # str(pen)
#'  }
#'
process_plink <- function(data_dir,
                          data_prefix,
                          rds_dir = data_dir,
                          rds_prefix,
                          logfile = NULL,
                          impute = TRUE,
                          impute_method = 'mode',
                          id_var = "IID",
                          quiet = FALSE,
                          overwrite = FALSE,
                          ...){

  if (identical(rds_prefix, data_prefix)) stop("rds_prefix cannot be the same as data_prefix. You need to change your choice of argument to rds_prefix.\n")

  # start log -----------------------------------------
  if(!is.null(logfile)){
    # TODO there seems to be a bug here -- needs investigation
    logfile <- create_log(file.path(rds_dir, logfile))
    cat("\nLogging to", logfile)
  } else {
    logfile <- tempfile()
  }


  if(!quiet){
    cat("\nPreprocessing", data_prefix, "data:")
  }
  cat("\nPreprocessing", data_prefix, "data\n", file = logfile, append = TRUE)


  # handle overwrite ----------------------------------
  if (overwrite) {
    to_remove <- list.files(file.path(rds_dir, rds_prefix))
    file.remove(to_remove)
  }

  # read in PLINK files --------------------------------
  step1 <- read_plink_files(data_dir = data_dir,
                            data_prefix = data_prefix,
                            rds_dir = rds_dir,
                            outfile = logfile,
                            overwrite = overwrite,
                            quiet = quiet)

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

  # notify about missing (genotype) values ------------------------------
  na_idx <- step2$na_counts > 0
  prop_na <- step2$na_counts/nrow(step2$X)

  cat("There are a total of", sum(na_idx), "SNPs with missing values\n",
      file = logfile, append = TRUE)
  cat("Of these,", sum(prop_na > 0.5),
      "are missing in at least 50% of the samples\n",
      file = logfile, append = TRUE)
  if(!quiet){
    cat("There are a total of", sum(na_idx), "SNPs with missing values\n")
    cat("Of these,", sum(prop_na > 0.5), "are missing in at least 50% of the samples\n")
  }

  # imputation ------------------------------------------------------------------
  step3 <- impute_snp_data(step2$obj,
                           step2$X,
                           impute = impute,
                           impute_method = impute_method,
                           outfile = logfile,
                           quiet = quiet, ...)

  bigsnpr::snp_save(step3) # save imputed data

  # format return object --------------------------------------------------------

  X  <- bigmemory::deepcopy(
    x = fbm2bm(step3$genotypes),
    row = 1:nrow(step3$genotypes),
    col = 1:ncol(step3$genotypes),
    type = 'double', # this is why we're making a copy -- need type 'double' downstream
    backingpath = rds_dir,
    backingfile = paste0(rds_prefix, '.bk'),
    descriptorfile = paste0(rds_prefix, '.desc')
  )

  ret <- list(X = describe(X),
              map = step3$map,
              fam = step3$fam,
              n = step3$n,
              p = step3$p)

  structure(ret, class = "processed_plink")

  rds_filename <- paste0(rds_prefix, ".rds")
  saveRDS(ret, file = file.path(rds_dir, rds_filename))

  # cleanup --------------------------------------------------------------------
  # These steps remove intermediate rds/bk files created by the steps of the data management process
  list.files(rds_dir, pattern=paste0('^file.*.bk'), full.names=TRUE) |>
    file.remove()
  list.files(rds_dir, pattern = paste0('^', data_prefix, ".rds"), full.names = TRUE) |>
    file.remove()
  list.files(rds_dir, pattern = paste0('^', data_prefix, ".bk"), full.names = TRUE) |>
    file.remove()

  gc()


  if(!quiet){cat("\nprocess_plink() completed \nProcessed files now saved as",
                 file.path(rds_dir, rds_filename))}

  cat("\nprocess_plink() completed. \nProcessed files now saved as",
      file.path(rds_dir, rds_filename),
      "at", pretty_time(),
      file = logfile, append = TRUE)

  # file.remove(tempfile()) # TODO: Not needed? Think about best way to handle logfile...
  return(file.path(rds_dir, rds_filename))

}



