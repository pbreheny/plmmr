#' A function to create a design matrix to be passed to a model fitting function
#'
#' @param dat                 Filepath to rds file of processed data
#' @param rds_dir               The path to the directory in which you want to create the new '.rds' and '.bk' files.
#' @param prefix              Optional user-specified prefix for the to-be-created .rds/.bk files. Must be different from any existing .rds/.bk files in the same folder.
#' @param is_bigsnp
#' @param add_predictor_fam   Optional: if you want to include "sex" (the 5th column of `.fam` file) in the analysis, specify 'sex' here.
#' @param add_predictor_ext   Optional: add additional covariates/predictors/features from an external file (i.e., not a PLINK file).
#'                            This argument takes one of two kinds of arguments:
#'                              - a **named** numeric vector, where the names align with the sample IDs in the PLINK files. **No NA values** can be in this vector.
#'                            The names will be used to subset and align this external covariate with the supplied PLINK data. **No NA values** can be in this matrix.
#'                              - a numeric matrix whose row names align with the sample IDs in the PLINK files.
#'                            The names will be used to subset and align this external covariate with the supplied PLINK data.
#' @param id_var
#' @param overwrite
#' @param outfile
#' @param quiet               Logical: should messages to be printed to the console be silenced? Defaults to FALSE
#'
#' @return
#' @export
#'
#' @examples
create_design_matrix <- function(dat,
                                 rds_dir,
                                 prefix,
                                 is_bigsnp,
                                 add_predictor_fam = NULL,
                                 add_predictor_ext = NULL,
                                 id_var = "IID",
                                 overwrite = FALSE,
                                 outfile = NULL,
                                 quiet = FALSE){

  # initial setup --------------------------------------------------


  # if(missing(outfile)){
  #   outfile = file.path(data_dir, "process_plink")
  # }

  logfile <- create_log(outfile = outfile)

  # check for files to be overwritten---------------------------------
  if (overwrite){
    # double check how much this will erase
    rds_to_remove <-  list.files(rds_dir,
                                 pattern=paste0('^', prefix, '.*.rds'),
                                 full.names=TRUE)
    if (length(rds_to_remove) > 1) {
      stop("You set overwrite=TRUE, but this looks like it will overwrite multiiple .rds files
      that have the 'std_prefix' pattern.
           To save you from overwriting anything important, I will not erase anything yet.
           Please move any .rds files with this file name pattern to another directory.")
    }
    file.remove(rds_to_remove)
    gc()

    bk_to_remove <-  list.files(rds_dir,
                                pattern=paste0('^', prefix, '.*.bk'),
                                full.names=TRUE)
    file.remove(bk_to_remove)
    gc()
  } else {
    gc()
  }

  if (is_bigsnp) {
    obj <- bigsnpr::snp_attach(rdsfile = dat)
  } else {
    stop("This function currently only works for data coming from process_plink().
         We are working to add the option to pass other kinds of data. For now,
         keep using process_delim()")
  }

  if (id_var == "FID"){
    obj$rownames <- og_plink_ids <- as.character(obj$fam$family.ID)
  } else if (id_var == "IID") {
    obj$rownames <- og_plink_ids <- as.character(obj$fam$sample.ID)
  }
  # add predictors from external files -----------------------------
  pred_X <- add_predictors_to_bigsnp(obj = obj,
                                     add_predictor_fam = add_predictor_fam,
                                     add_predictor_ext = add_predictor_ext,
                                     id_var = id_var,
                                     og_plink_ids = og_plink_ids,
                                     rds_dir = rds_dir,
                                     quiet = quiet)
  gc()

  # subsetting -----------------------------------------------------------------
  subset_X <- subset_bigsnp(obj = pred_X$obj,
                            handle_missing_phen = obj$handle_missing_phen,
                            complete_phen = obj$complete_phen,
                            non_gen = pred_X$non_gen,
                            rds_dir = rds_dir,
                            prefix = prefix,
                            outfile = logfile,
                            quiet = quiet)
  gc()
  # standardization ------------------------------------------------------------
  std_X <- standardize_bigsnp(obj = subset_X,
                              prefix = prefix,
                              rds_dir = rds_dir,
                              non_gen = pred_X$non_gen,
                              complete_phen = obj$complete_phen,
                              id_var = id_var,
                              outfile = logfile,
                              quiet = quiet,
                              overwrite = overwrite)

  # cleanup and output ---------------------------------------------------------
  std_X$std_X_colnames <- pred_X$obj$colnames[std_X$ns]
  std_X$std_X_rownames <- pred_X$obj$rownames[pred_X$obj$complete_phen]
  std_X$X_colnames <- pred_X$obj$colnames
  std_X$X_rownames <- pred_X$obj$rownames
  std_X$n <- nrow(obj$fam)
  std_X$p <- nrow(obj$map)
  std_X$fam <- obj$fam
  std_X$map <- obj$map
  std_X$non_gen <- pred_X$non_gen # save indices for non-genomic covariates
  std_X$complete_phen <- pred_X$cobj$omplete_phen # save indices for which samples had complete phenotypes
  std_X$id_var <- id_var # save ID variable - will need this downstream for analysis

  saveRDS(std_X, file.path(rds_dir, paste0(prefix, ".rds")))
  return(file.path(rds_dir, paste0(prefix, ".rds")))
}



