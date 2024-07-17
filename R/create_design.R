#' A function to create a design matrix to be passed to a model fitting function
#'
#' @param dat                 Filepath to rds file of processed data
#' @param rds_dir             The path to the directory in which you want to create the new '.rds' and '.bk' files.
#' @param new_file            User-specified filename (*without .bk/.rds extension*) for the to-be-created .rds/.bk files. Must be different from any existing .rds/.bk files in the same folder.
#' @param is_bigsnp           Logical: is the data coming from `process_plink()`? For now, this must be TRUE. Working to expand to accommodate other kinds of data.
#' @param na_phenotype_vals   A vector of numeric values used to code NA values in the phenotype/outcome (this is the 'affection' column in a `bigSNP` object, or the last column of a `.fam` file). Defaults to -9 (matching PLINK conventions).
#' @param add_phen            Optional: A **data frame** with at least two columns: and ID column and a phenotype column
#' @param pheno_id            Optional: A string specifying the name of the ID column in `pheno`. MUST be specified if `add_phen` is specified.
#' @param pheno_name           Optional: A string specifying the name of the phenotype column in `pheno`.  MUST be specified if `add_phen` is specified. This column will be used as the default `y` argument to 'plmm()'.
#' @param handle_missing_phen A string indicating how missing phenotypes should be handled:
#'                                * "prune" (default): observations with missing phenotype are removed
#'                            Note: for data coming from PLINK, no missing values of the phenotype are allowed. You have to (1) supply phenotype from an external file,
#'                            or (2) prune missing values
#' @param add_predictor_fam   Optional: if you want to include "sex" (the 5th column of `.fam` file) in the analysis, specify 'sex' here.
#' @param add_predictor_ext   Optional: add additional covariates/predictors/features from an external file (i.e., not a PLINK file).
#'                            This argument takes one of two kinds of arguments:
#'                              - a **named** numeric vector, where the names align with the sample IDs in the PLINK files. **No NA values** can be in this vector.
#'                            The names will be used to subset and align this external covariate with the supplied PLINK data. **No NA values** can be in this matrix.
#'                              - a numeric matrix whose row names align with the sample IDs in the PLINK files.
#'                            The names will be used to subset and align this external covariate with the supplied PLINK data.
#' @param id_var               String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID"
#' @param overwrite           Logical: should existing .rds files be overwritten? Defaults to FALSE.
#' @param logfile             Optional: name of the '.log' file to be written
#' @param quiet               Logical: should messages to be printed to the console be silenced? Defaults to FALSE
#'
#' @export
#'
create_design <- function(dat,
                          rds_dir,
                          new_file,
                          is_bigsnp,
                          na_phenotype_vals = c(-9),
                          handle_missing_phen = "prune",
                          add_phen = NULL,
                          pheno_id = NULL,
                          pheno_name = NULL,
                          add_predictor_fam = NULL,
                          add_predictor_ext = NULL,
                          id_var = "IID",
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
      that have the 'new_file' pattern.
           To save you from overwriting anything important, I will not erase anything yet.
           Please move any .rds files with this file name pattern to another directory.")
    }
    file.remove(rds_to_remove)
    gc()

    bk_to_remove <-  list.files(rds_dir,
                                pattern=paste0('^', new_file, '.*.bk'),
                                full.names=TRUE)
    file.remove(bk_to_remove)
    gc()
  } else {
    gc()
  }

  if (is_bigsnp) {
    obj <- readRDS(dat)
    obj$X <- bigmemory::attach.big.matrix(obj$X)
  } else {
    stop("This function currently only works for data coming from process_plink().
         We are working to add the option to pass other kinds of data. For now,
         keep using process_delim()")
  }

  # save items to return
  # design$n <- nrow(obj$fam)
  # design$p <- nrow(obj$map)

  # determine which IDs to use ---------------------------------
  if (id_var == "IID"){
    geno_id <- "sample.ID"
  } else if (id_var == "FID"){
    geno_id <- "family.ID"
  } else {
    stop("\nThe argument to id_var is misspecified. Must be one of 'IID' or 'FID', and
         the corresponding variable *must* be of type 'char'.")
  }

  if (!is.null(add_phen) & is.null(colnames(add_phen))) {
    stop("If you supply an external phenotype file to add_phen, this matrix must have
         column names.")
  }

  og_plink_ids <- obj$fam[,geno_id] |> as.character()

  # add phenotype from external files -----------------------------
  if (!is.null(add_phen)){
    step1 <- add_external_phenotype(geno = obj,
                                    rds_dir = rds_dir,
                                    geno_id = geno_id,
                                    pheno = add_phen,
                                    pheno_id = pheno_id,
                                    pheno_col = pheno_name,
                                    outfile = logfile,
                                    quiet = quiet)
  } else {
    step1 <- obj
  }
  gc() # cleanup

  # address missingness in phenotype values -----------------------
  step2 <- handle_missingness(obj = step1,
                              na_phenotype_vals = na_phenotype_vals,
                              handle_missing_phen = handle_missing_phen,
                              outfile = logfile,
                              quiet = quiet)



  # save items to return
  design$std_X_rownames <- og_plink_ids[step2$complete_phen]
  design$complete_phen <- step2$complete_phen # save indices for which samples had complete phenotypes
  design$ns_genotypes <- step2$ns_genotypes

  # add predictors from external files -----------------------------
  step1$obj$colnames <- obj$colnames # need these to carry to next step
 # if (!is.null(add_predictor_fam) | !is.null(add_predictor_ext)){
    pred_X <- add_predictors_to_bigsnp(obj = step1,
                                       add_predictor_fam = add_predictor_fam,
                                       add_predictor_ext = add_predictor_ext,
                                       id_var = id_var,
                                       og_plink_ids = og_plink_ids,
                                       rds_dir = rds_dir,
                                       quiet = quiet)

    # save items to return
    design$map <- pred_X$obj$map
    design$non_gen <- pred_X$non_gen # save indices for non-genomic covariates
    design$X_colnames <- pred_X$obj$colnames
    design$X_rownames <- pred_X$obj$rownames
    design$fam <- pred_X$obj$fam
    design$y <- pred_X$obj$fam[,6]

    # again, clean up to save space
    rm(step1); gc()

  #} else {
  #   pred_X <- step1
  #   design$non_gen <- NULL
  #
  #   # again, clean up to save space
  #   rm(step1); gc()
  # }


  # subsetting -----------------------------------------------------------------
  subset_X <- subset_bigsnp(obj = pred_X$obj,
                            handle_missing_phen = handle_missing_phen,
                            complete_phen = design$complete_phen,
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
browser()
  # standardization ------------------------------------------------------------
  std_res <- standardize_bigsnp(obj = subset_X,
                               new_file = new_file,
                               rds_dir = rds_dir,
                               non_gen = design$non_gen,
                               complete_phen = design$complete_phen,
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



  # cleanup -------------------------------------------------------------
  # These steps remove intermediate rds/bk files created by the steps of the data management process
  # list.files(rds_dir, pattern=paste0('^', new_file, '.*.rds'), full.names=TRUE) |>
  #   file.remove()
  # gc()
  # list.files(rds_dir, pattern=paste0('^', new_file, '.*.bk'), full.names=TRUE) |>
  #   file.remove()
  # gc()
  # list.files(rds_dir, pattern=paste0('^file.*.bk'), full.names=TRUE) |>
  #   file.remove()
  # gc()
  # Note: calls to gc() are important!


  # return -------------------------------------------------------------
  saveRDS(design, file.path(rds_dir, paste0(new_file, ".rds")))
  return(file.path(rds_dir, paste0(new_file, ".rds")))

}



