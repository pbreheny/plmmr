#' A function to create a design matrix, outcome, and penalty factor to be passed to a model fitting function
#'
#' @param dat_file                Filepath to rds file of processed data (data from `process_plink()` or `process_delim()`)
#' @param rds_dir                 The path to the directory in which you want to create the new '.rds' and '.bk' files.
#' @param new_file                User-specified filename (*without .bk/.rds extension*) for the to-be-created .rds/.bk files. Must be different from any existing .rds/.bk files in the same folder.
#' @param feature_id                  This argument is the ID for the feature data. Two options here:
#'                                  - for PLINK data: a string specifying an ID column of the PLINK `.fam` file. Options are "IID" (default) and "FID"
#'                                  - for all other data: a character vector of unique identifiers (IDs) for each row of the feature data (i.e., the data processed with `process_delim()`)
#' @param add_outcome             A data frame or matrix with two columns: and ID column and a column with the outcome value (to be used as 'y' in the final design). IDs must be characters, outcome must be numeric.
#' @param outcome_id              A string specifying the name of the ID column in 'add_outcome'
#' @param outcome_col             A string specifying the name of the phenotype column in 'add_outcome'
#' @param na_outcome_vals         A vector of numeric values used to code NA values in the outcome. Defaults to `c(-9, NA_integer)` (the -9 matches PLINK conventions).
#' @param add_predictor       Optional: a matrix or data frame to be used for adding additional covariates/predictors/features from an external file (i.e., not a PLINK file).
#'                                This matrix must have one column that is an ID column; all other columns aside the ID will be used as covariates in the design matrix. Columns must be named.
#' @param predictor_id            A string specifying the name of the column in 'add_predictor' with sample IDs. Required if 'add_predictor' is supplied.
#'                                The names will be used to subset and align this external covariate with the supplied PLINK data.
#' @param overwrite               Logical: should existing .rds files be overwritten? Defaults to FALSE.
#' @param logfile                 Optional: name of the '.log' file to be written
#' @param quiet                   Logical: should messages to be printed to the console be silenced? Defaults to FALSE
#'
#' @export
#'
create_design <- function(dat_file,
                          rds_dir,
                          new_file,
                          feature_id,
                          add_outcome,
                          outcome_id,
                          outcome_col,
                          na_outcome_vals = c(-9, NA_integer_),
                          add_predictor = NULL,
                          predictor_id = NULL,
                          logfile = NULL,
                          overwrite = FALSE,
                          quiet = FALSE){

  # check for input errors ----------------------------------------

  if (any(add_outcome[,outcome_col] %in% na_outcome_vals)) {
    stop('It appears that you have some missing values in the outcome data.
         Please remove these samples; missing values are not permitted in the design.')
  }

  if (is.null(colnames(add_outcome))) {
    stop('The columns of "add_outcome" must be named.')
  }

  # additional checks for case where add_predictor is specified
  if (!is.null(add_predictor)) {

    if (is.null(predictor_id)) {
      stop("If add_predictor is specified, the user must also specify predictor_id")
    }

    if (is.null(colnames(add_predictor))){
      stop('The columns of "add_predictor" must be named.')
    }

    if (any(is.na(add_predictor[,]))) {
      stop('It appears that there are missing values in the predictor data.
         Please remove these samples; missing values are not permitted in the design.')
    }

    if (!identical(add_outcome[,outcome_id], add_predictor[,predictor_id])) {
      stop("Something is off in the supplied outcome and/or predictor data.
         Make sure the indicated ID columns are character type, represent the same samples, and have the same order.")
    }

  }

  # initial setup --------------------------------------------------
  existing_files <- list.files(rds_dir)

  if(!is.null(logfile)){
    logfile = file.path(rds_dir, logfile)
  }

  logfile <- create_log(outfile = logfile)

  # create list to be returned
  design <- list()

  # check for files to be overwritten---------------------------------
  if (overwrite){

    # remove files with name pattern
    to_remove <- paste0(file.path(rds_dir, new_file), c('.bk', '.rds', '.desc'))
    if (any(file.exists(to_remove))) {
      file.remove(to_remove)
    }


    # check for left over intermediate files
    if (file.exists(file.path(rds_dir,'unstd_design_matrix.bk'))) {
      file.remove(c(file.path(rds_dir,'unstd_design_matrix.bk'),
                    file.path(rds_dir,'unstd_design_matrix.desc')))
    }


  }

  # read in the processed data -------------------------------
  obj <- readRDS(dat_file)
  obj$X <- bigmemory::attach.big.matrix(obj$X)

  # save these original dim names
  design$X_colnames <- obj$colnames
  design$X_rownames <- obj$rownames

  # flag for data coming from plink
  is_plink <- any(grepl('fam', names(obj)))

  # determine which IDs to use ---------------------------------
  if (feature_id == "IID"){
    indiv_id <- "sample.ID"
    og_ids <- obj$fam[,indiv_id] |> as.character()

  } else if (feature_id == "FID"){
    indiv_id <- "family.ID"
    og_ids <- obj$fam[,indiv_id] |> as.character()

  } else if (is.character(feature_id) & length(feature_id) == nrow(obj$X)) {
    og_ids <- feature_id

  } else {
    stop("The feature_id argument is either misspecified or missing (see documentation for options).")
  }

  if (is.null(colnames(add_outcome))) {
    stop("The matrix supplied to add_outcome must have column names.")
  }

  # index samples for subsetting ------------
  # Note: this step uses the outcome (from external file) to determine which
  #   samples/observations should be pruned out; observations with no feature
  #   data will be removed from analysis
  sample_idx <- index_samples(obj = obj,
                              rds_dir = rds_dir,
                              indiv_id = og_ids,
                              add_outcome = add_outcome,
                              outcome_id = outcome_id,
                              outcome_col = outcome_col,
                              na_outcome_vals = na_outcome_vals,
                              outfile = logfile,
                              quiet = quiet)

  # save items to return
  design$outcome_idx <- sample_idx$outcome_idx # save indices of which rows in the feature data should be included in the design
  design$y <-  sample_idx$complete_samples[,..outcome_col]
  design$std_X_rownames <- sample_idx$complete_samples$ID

  if (!is.null(add_predictor)) {
    if (is.null(predictor_id)) stop('If add_predictor is supplied, the predictor_id argument must also be supplied')

    # align IDs between feature data and external data -------------------------
    aligned_add_predictor <- align_ids(id_var = predictor_id,
                                       quiet = quiet,
                                       add_predictor = add_predictor,
                                       og_ids = og_ids)

    # add predictors from external files --------------------------------------
    unstd_X <- add_predictors(obj = obj,
                              add_predictor = aligned_add_predictor,
                              id_var = feature_id,
                              rds_dir = rds_dir,
                              quiet = quiet)

    # save items to return
    design$non_gen <- unstd_X$non_gen # save indices for non-genomic covariates
    design$non_gen_colnames <- setdiff(colnames(add_predictor), predictor_id)

    if (is_plink){
      design$fam <- unstd_X$obj$fam
      design$map <- unstd_X$obj$map
    }

    # again, clean up to save space
    rm(obj)

  # index features for subsetting --------------------------------------------

  design$ns <- count_constant_features(fbm = unstd_X$design_matrix,
                                                 outfile = logfile,
                                                 quiet = quiet)

  # subsetting -----------------------------------------------------------------
  subset_res <- subset_bigsnp(X = unstd_X$design_matrix,
                              complete_samples = design$outcome_idx,
                              ns = design$ns,
                              rds_dir = rds_dir,
                              new_file = new_file,
                              outfile = logfile,
                              quiet = quiet)
  # clean up
  design$ns <- subset_res$ns
  design$std_X_colnames <- unstd_X$colnames[subset_res$ns]
  rm(unstd_X)

  # standardization ------------------------------------------------------------
  std_res <- standardize_bigsnp(X = subset_res$subset_X,
                                new_file = new_file,
                                rds_dir = rds_dir,
                                outfile = logfile,
                                quiet = quiet,
                                overwrite = overwrite)
  rm(subset_res)
  gc()

  # add meta data -------------------------------------------------------------
  design$std_X <- std_res$std_X
  design$std_X_n <- std_res$std_X_n
  design$std_X_p <- std_res$std_X_p
  design$std_X_center <- std_res$std_X_center
  design$std_X_scale <- std_res$std_X_scale
  design$penalty_factor <- c(rep(0, length(design$non_gen)),
                             rep(1, design$std_X_p - length(design$non_gen)))

  # cleanup -------------------------------------------------------------------
  list.files(rds_dir,
             pattern=paste0('^unstd_design.*'),
             full.names=TRUE) |> file.remove()

  # return -------------------------------------------------------------
  saveRDS(design, file.path(rds_dir, paste0(new_file, ".rds")))
  return(file.path(rds_dir, paste0(new_file, ".rds")))
  } else {
    stop('still working through the case where no predictors are added')
  }
}



