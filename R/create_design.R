#' A function to create a design matrix, outcome, and penalty factor to be passed to a model fitting function
#'
#' @param dat                     Filepath to rds file of processed data (data from `process_plink()` or `process_delim()`)
#' @param rds_dir                 The path to the directory in which you want to create the new '.rds' and '.bk' files.
#' @param new_file                User-specified filename (*without .bk/.rds extension*) for the to-be-created .rds/.bk files. Must be different from any existing .rds/.bk files in the same folder.
#' @param add_outcome             A data frame or matrix with two columns: and ID column and a column with the outcome value (to be used as 'y' in the final design). IDs must be characters, outcome must be numeric.
#' @param outcome_id              A string specifying the name of the ID column in 'add_outcome'
#' @param outcome_col             A string specifying the name of the phenotype column in 'add_outcome'
#' @param id_var                  This argument is the ID:
#'                                  - for PLINK data: a string specifying an ID column of the PLINK `.fam` file. Options are "IID" (default) and "FID"
#'                                  - for all other data: a character vector of unique identifiers (IDs) for each row of the data (i.e., the data processed with `process_delim()`)
#' @param na_outcome_vals         A vector of numeric values used to code NA values in the outcome. Defaults to `c(-9, NA_integer)` (the -9 matches PLINK conventions).
#' @param handle_missing_outcome  A string indicating how missing phenotypes should be handled:
#'                                  * "prune" (default): observations with missing phenotype are removed
#'                                  Note: for data coming from PLINK, no missing values of the phenotype are allowed. You have to (1) supply phenotype from an external file,
#'                                  or (2) prune missing values
#' @param add_predictor_ext       Optional: add additional covariates/predictors/features from an external file (i.e., not a PLINK file).
#'                                This argument takes one of two kinds of arguments:
#'                                 - a **named** numeric vector, where the names align with the sample IDs in the PLINK files. **No NA values** can be in this vector.
#'                                  The names will be used to subset and align this external covariate with the supplied PLINK data. **No NA values** can be in this matrix.
#'                                 - a numeric matrix whose row names align with the sample IDs in the PLINK files.
#'                                The names will be used to subset and align this external covariate with the supplied PLINK data.
#' @param overwrite               Logical: should existing .rds files be overwritten? Defaults to FALSE.
#' @param logfile                 Optional: name of the '.log' file to be written
#' @param quiet                   Logical: should messages to be printed to the console be silenced? Defaults to FALSE
#'
#' @export
#'
create_design <- function(dat,
                          rds_dir,
                          new_file,
                          add_outcome,
                          outcome_id,
                          outcome_col,
                          id_var,
                          na_outcome_vals = c(-9, NA_integer_),
                          add_predictor_ext = NULL,
                          logfile = NULL,
                          overwrite = FALSE,
                          quiet = FALSE){

  # initial setup --------------------------------------------------
  existing_files <- list.files(rds_dir)
  if (any(grepl(paste0(new_file, ".bk"), existing_files, fixed = TRUE))) {
    stop("The new_file you specified already exists in rds_dir. Please choose a different file name.\n")
  }

  if(!is.null(logfile)){
    logfile = file.path(rds_dir, logfile)
  }

  logfile <- create_log(outfile = logfile)

  # create list to be returned
  design <- list()

  # check for files to be overwritten---------------------------------
  if (overwrite){
    # double check how much this will erase
    rds_to_remove <-  list.files(rds_dir,
                                 pattern=paste0('^', new_file, '.*.rds'),
                                 full.names=TRUE)
    if (length(rds_to_remove) > 1) {
      stop("You set overwrite=TRUE, but this looks like it will overwrite multiiple .rds files
      that have the 'new_file' pattern. To save you from overwriting anything important, I will not erase anything yet.
      Please move any .rds files with this file name pattern to another directory.")
    }
    file.remove(rds_to_remove)
    gc()

    bk_to_remove <-  list.files(rds_dir,
                                pattern=paste0('^', new_file, '.*.bk'),
                                full.names=TRUE)
    file.remove(bk_to_remove)
    gc()
  }

  # read in the processed data -------------------------------
  obj <- readRDS(dat)
  obj$X <- bigmemory::attach.big.matrix(obj$X)

  # determine which IDs to use ---------------------------------
  if (id_var == "IID"){
    indiv_id <- "sample.ID"
    og_ids <- obj$fam[,indiv_id] |> as.character()

  } else if (id_var == "FID"){
    indiv_id <- "family.ID"
    og_ids <- obj$fam[,indiv_id] |> as.character()

  } else if (is.character(id_var) & length(id_var) == nrow(obj$X)) {
    og_ids <- id_var

  } else {
    stop("The id_var argument is either misspecified or missing (see documentation for options).")
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

  gc() # cleanup

  # save items to return
  design$std_X_rownames <- og_ids[sample_idx$outcome_present]
  design$outcome_present <- sample_idx$outcome_present # save indices for which samples had nonmissing outcome values
  design$outcome_idx <- sample_idx$outcome_idx # save indices of which rows in the feature data should be included in the design

  # index features for subsetting --------------------------------------------
  design$ns_genotypes <- count_constant_features(fbm = obj$X,
                                                 ind.row = sample_idx$outcome_idx,
                                                 outfile = logfile,
                                                 quiet = quiet)

  # add predictors from external files -----------------------------
  pred_X <- add_predictors(obj = obj,
                           add_predictor_ext = add_predictor_ext,
                           id_var = id_var,
                           og_ids = og_ids,
                           rds_dir = rds_dir,
                           quiet = quiet)

  # save items to return
  design$non_gen <- pred_X$non_gen # save indices for non-genomic covariates
  design$X_colnames <- pred_X$obj$colnames
  design$X_rownames <- pred_X$obj$rownames
  design$fam <- pred_X$obj$fam
  design$map <- pred_X$obj$map
  # design$y <- pred_X$obj$fam[,6] # TODO: fix this behavior

  # again, clean up to save space
  rm(obj); gc()



  # subsetting -----------------------------------------------------------------
  subset_X <- subset_bigsnp(obj = pred_X$obj,
                            handle_missing_outcome = handle_missing_outcome,
                            outcome_present = design$outcome_present,
                            non_gen = design$non_gen,
                            ns_genotypes = design$ns_genotypes,
                            rds_dir = rds_dir,
                            new_file = new_file,
                            outfile = logfile,
                            quiet = quiet)
  # clean up
  rm(pred_X); gc()
  design$ns <- subset_X$ns
  design$std_X_colnames <- design$X_colnames[subset_X$ns]

  # standardization ------------------------------------------------------------
  std_res <- standardize_bigsnp(obj = subset_X,
                                new_file = new_file,
                                rds_dir = rds_dir,
                                non_gen = design$non_gen,
                                outcome_present = design$outcome_present,
                                id_var = id_var,
                                outfile = logfile,
                                quiet = quiet,
                                overwrite = overwrite)
  rm(subset_X)
  gc()

  # add meta data -------------------------------------------------------------
  design$std_X <- std_res$std_X
  design$std_X_n <- std_res$std_X_n
  design$std_X_p <- std_res$std_X_p
  design$std_X_center <- std_res$std_X_center
  design$std_X_scale <- std_res$std_X_scale
  design$ns <- std_res$ns
  design$penalty_factor <- c(rep(0, length(design$non_gen)),
                             rep(1, design$std_X_p - length(design$non_gen)))

  # return -------------------------------------------------------------
  saveRDS(design, file.path(rds_dir, paste0(new_file, ".rds")))
  return(file.path(rds_dir, paste0(new_file, ".rds")))

}



