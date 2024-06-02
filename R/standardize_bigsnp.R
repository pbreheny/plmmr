#' A helper function to standardize a `bigSNP`
#'
#' @param obj           A `bigSNP` object
#' @param prefix        The prefix (as a character string) of the bed/fam data files (e.g., `prefix = 'mydata'`)
#' @param rds_dir       The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`
#' @param non_gen       An integer vector that ranges from 1 to the number of added predictors. Example: if 2 predictors are added, non_gen = 1:2.
#' Note: this is typically passed from the result of `add_predictors()`
#' @param complete_phen Numeric vector with indicies marking the rows of the original data which have a non-missing entry in the 6th column of the `.fam` file
#' @param id_var        String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID".
#' @param outfile       Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet         Logical: should messages be printed to the console? Defaults to TRUE
#' @param overwrite     Logical: if existing `.bk`/`.rds` files exist for the specified directory/prefix, should these be overwritten?
#'
#' @return A list with a new component of `obj` called 'std_X' - this is an FBM with column-standardized data.
#' List also includes several other indices/meta-data on the standardized matrix
#' @keywords internal
#'
standardize_bigsnp <- function(obj, prefix, rds_dir, non_gen, complete_phen, id_var,
                               outfile, quiet, overwrite){
  # check for files to be overwritten
  if (overwrite){
    list.files(rds_dir, pattern=paste0('^std_.*.bk'), full.names=TRUE) |>
      file.remove()
    list.files(rds_dir, pattern=paste0('^std_.*.rds'), full.names=TRUE) |>
      file.remove()
  }


  # standardization ------------------------------------------------
  if (!quiet) {cat("\nColumn-standardizing the design matrix...")}
  # centering & scaling
  subset_X <- bigstatsr::big_copy(obj$subset_X,
                             type = "double", # this is key...
                             backingfile = file.path(rds_dir, "std_X"))
  subset_X_bm <- subset_X |> fbm2bm()

  std_res <- .Call("big_std",
                   subset_X_bm@address,
                   as.integer(bigstatsr::nb_cores()),
                   PACKAGE = "plmmr")

  std_X <- bigmemory::big.matrix(nrow = nrow(subset_X), ncol = ncol(subset_X))
  std_X@address <- std_res$std_X


  # label return object ------------------------------------------------
  # naming these center and scale values so that I know they relate to the first
  # standardization; there will be another standardization after the rotation
  # in plmm_fit().

  ret <- list(
    std_X = describe(std_X),
    std_X_center = std_res$std_X_center,
    std_X_scale = std_res$std_X_scale,
    std_X_colnames = obj$colnames[obj$ns],
    std_X_rownames = obj$rownames[complete_phen],
    fam = obj$fam,
    map = obj$map,
    non_gen = non_gen, # save indices for non-genomic covariates
    complete_phen = complete_phen, # save indices for which samples had complete phenotypes
    id_var = id_var # save ID variable - will need this downstream for analysis
  )


  if (!quiet){
    cat("\nDone with standardization. File formatting in progress.",
        file = outfile, append = TRUE)
  }

  return(ret)
}