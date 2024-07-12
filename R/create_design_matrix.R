#' A function to create a design matrix to be passed to a model fitting function
#'
#' @param X
#' @param add_predictor_fam   Optional: if you want to include "sex" (the 5th column of `.fam` file) in the analysis, specify 'sex' here.
#' @param add_predictor_ext   Optional: add additional covariates/predictors/features from an external file (i.e., not a PLINK file).
#'                            This argument takes one of two kinds of arguments:
#'                              - a **named** numeric vector, where the names align with the sample IDs in the PLINK files. **No NA values** can be in this vector.
#'                            The names will be used to subset and align this external covariate with the supplied PLINK data. **No NA values** can be in this matrix.
#'                              - a numeric matrix whose row names align with the sample IDs in the PLINK files.
#'                           The names will be used to subset and align this external covariate with the supplied PLINK data.
#' @param id_var
#' @param overwrite
#' @param outfile
#'
#' @return
#' @export
#'
#' @examples
create_design_matrix <- function(X,
                                 add_predictor_fam = NULL,
                                 add_predictor_ext = NULL,
                                 id_var = "IID",
                                 overwrite = FALSE,
                                 outfile = NULL){

  # check for files to be overwritten---------------------------------
  if (overwrite){
    gc()
    # double check how much this will erase
    rds_to_remove <-  list.files(rds_dir, pattern=paste0('^std_.*.rds'), full.names=TRUE)
    if (length(rds_to_remove) > 1) {
      stop("You set overwrite=TRUE, but this looks like it will overwrite multiiple .rds files
      that have the 'std_prefix' pattern.
           To save you from overwriting anything important, I will not erase anything yet.
           Please move any .rds files with this file name pattern to another directory.")
    }
    file.remove(rds_to_remove)
    gc()
    list.files(rds_dir, pattern=paste0('^std_.*.rds'), full.names=TRUE) |>
      file.remove()
    gc()
  } else {
    gc()
  }

  # add predictors from external files -----------------------------
  step5 <- add_predictors_to_bigsnp(obj = step4,
                                    add_predictor_fam = add_predictor_fam,
                                    add_predictor_ext = add_predictor_ext,
                                    id_var = id_var,
                                    og_plink_ids = step2$og_plink_ids,
                                    rds_dir = rds_dir,
                                    quiet = quiet)
  gc()

  # subsetting -----------------------------------------------------------------
  step6 <- subset_bigsnp(obj = step5$obj,
                         counts = step2$counts,
                         handle_missing_phen = handle_missing_phen,
                         complete_phen = step3$complete_phen,
                         non_gen = step5$non_gen,
                         data_dir = data_dir,
                         rds_dir = rds_dir,
                         prefix = prefix,
                         bk_filename = bk_filename,
                         outfile = logfile,
                         quiet = quiet)
  gc()
  # standardization ------------------------------------------------------------
  step7 <- standardize_bigsnp(obj = step6,
                              prefix = prefix,
                              rds_dir = rds_dir,
                              non_gen = step5$non_gen,
                              complete_phen = step3$complete_phen,
                              id_var = id_var,
                              outfile = logfile,
                              quiet = quiet,
                              overwrite = overwrite)

  # cleanup and output ---------------------------------------------------------
  saveRDS(step7, file.path(rds_dir, paste0("std_", prefix, ".rds")))

}