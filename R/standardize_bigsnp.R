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


  # standardization ------------------------------------------------
  if (!quiet) {cat("Column-standardizing the design matrix...\n")}
  # convert FBM pointer into a big.matrix pointer
  subset_X_bm <- obj$subset_X |> fbm2bm()
  # centering & scaling
  std_res <- .Call("big_std",
                   subset_X_bm@address,
                   as.integer(bigstatsr::nb_cores()),
                   PACKAGE = "plmmr")

  std_X <- bigmemory::big.matrix(nrow = nrow(obj$subset_X), ncol = ncol(obj$subset_X))
  std_X@address <- std_res$std_X

  # checks ---------------------------------------------------------
  # std_X and obj$subset_X are two pointers, both pointing to the same backing file
  # obj$subset_X$backingfile; paste0(dir.name(std_X), file.name(std_X))

  # foo1 <- obj$subset_X[,]
  # foo2 <- std_X[,]
  # tinytest::expect_equivalent(foo1, foo2) #

  # label return object ------------------------------------------------
  # naming these center and scale values so that I know they relate to the first
  # standardization; there will be another standardization after the rotation
  # in plmm_fit().

  ret <- list(
    std_X = bigmemory::describe(std_X),
    std_X_center = std_res$std_X_center,
    std_X_scale = std_res$std_X_scale,
    ns = obj$ns,
    std_X_colnames = obj$colnames[obj$ns],
    std_X_rownames = obj$rownames[complete_phen],
    X_colnames = obj$colnames,
    X_rownames = obj$rownames,
    n = nrow(obj$fam),
    p = nrow(obj$map),
    fam = obj$fam,
    map = obj$map,
    non_gen = non_gen, # save indices for non-genomic covariates
    complete_phen = complete_phen, # save indices for which samples had complete phenotypes
    id_var = id_var # save ID variable - will need this downstream for analysis
  )


  if (!quiet){
    cat("Done with standardization. File formatting in progress\n")
  }

  return(ret)
}