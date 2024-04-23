#' Preprocess PLINK files using the `bigsnpr` package
#' 
#' @param data_dir The path to the bed/bim/fam data files, *without* a trailing "/" (e.g., use `data_dir = '~/my_dir'`, **not** `data_dir = '~/my_dir/'`)
#' @param prefix The prefix (as a character string) of the bed/fam data files (e.g., `prefix = 'mydata'`)
#' @param impute Logical: should data be imputed? Default to TRUE.
#' @param impute_method If 'impute' = TRUE, this argument will specify the kind of imputation desired. Options are: 
#'  * mode (default): Imputes the most frequent call. See `bigsnpr::snp_fastImputeSimple()` for details. 
#'  * random: Imputes sampling according to allele frequencies.
#'  * mean0: Imputes the rounded mean.
#'  * mean2: Imputes the mean rounded to 2 decimal places.
#'  * xgboost: Imputes using an algorithm based on local XGBoost models. See `bigsnpr::snp_fastImpute()` for details. Note: this can take several minutes, even for a relatively small data set. 
#' @param na_phenotype_vals A vector of numeric values used to code NA values in the phenotype/outcome (this is the 'affection' column in a `bigSNP` object, or the last column of a `.fam` file). Defaults to -9 (matching PLINK conventions).
#' @param id_var String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID". 
#' @param handle_missing_phen A string indicating how missing phenotypes should be handled: 
#'  * "prune" (default): observations with missing phenotype are removed
#'  * "asis": leaves missing phenotypes as NA (this is fine if outcome will be supplied later from a separate file)
#'  * "median": impute missing phenotypes using the median (warning: this is overly simplistic in many cases).
#'  * "mean": impute missing phenotypes using the mean (warning: this is overly simplistic in many cases).
#' @param quiet Logical: should messages be printed to the console? Defaults to TRUE
#' @param gz Logical: are the bed/bim/fam files g-zipped? Defaults to FALSE. NOTE: if TRUE, process_plink will unzip your zipped files.
#' @param outfile Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param overwrite Logical: if existing `.bk`/`.rds` files exist for the specified directory/prefix, should these be overwritten? Defaults to FALSE. Set to TRUE if you want to change the imputation method you're using, etc. 
#' @param ... Optional: additional arguments to `bigsnpr::snp_fastImpute()` (relevant only if impute_method = "xgboost")
#' @param add_predictor_fam Optional: if you want to include "sex" (the 5th column of `.fam` file) in the analysis, specify 'sex' here.
#' @param add_predictor_ext Optional: add additional covariates/predictors/features from an external file (i.e., not a PLINK file). 
#' This argument takes one of two kinds of arguments: 
#'  (a) a **named** numeric vector, where the names align with the sample IDs in the PLINK files. 
#' The names will be used to subset and align this external covariate with the supplied PLINK data.
#' 
#'  (b) a numeric matrix whose row names align with the sample IDs in the PLINK files. 
#'  The names will be used to subset and align this external covariate with the supplied PLINK data.
#' @returns Nothing is returned by this function; instead, files 'prefix.rds',
#'  'prefix.bk', and std_prefix.bk are created in the location specified by data_dir. Note that this 
#'  this function need only be run once; in subsequent data analysis/scripts, 
#'  `get_data()` will access the '.rds' file. 
#'    
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' process_plink(data_dir = plink_example(parent = T),
#' prefix = "penncath_lite",
#' gz = TRUE,
#' outfile = "process_penncath",
#' overwrite = TRUE,
#' impute_method = "mode")
#' }
#' 
process_plink <- function(data_dir,
                          prefix,
                          impute = TRUE,
                          impute_method = 'mode',
                          na_phenotype_vals = c(-9),
                          id_var = "IID",
                          handle_missing_phen = "prune",
                          quiet = FALSE,
                          gz = FALSE,
                          outfile,
                          overwrite = FALSE,
                          add_predictor_fam = NULL,
                          add_predictor_ext = NULL,
                          ...){
  
  # start log ------------------------------------------
  if(missing(outfile)){
    outfile = paste0(data_dir, "/process_plink.log")
  } else {
    outfile = paste0(outfile, ".log")
  }
  log_con <- file(outfile)
  cat("### Processing PLINK files for PLMM ###", file = log_con)
  cat("\nLogging to ", outfile, file = outfile, append = TRUE)
  cat("\nPreprocessing", prefix, "data:", file = outfile, append = TRUE)
  
  if(!quiet){
    cat("\nLogging to", outfile)
    cat("\nPreprocessing", prefix, "data:")
  }
  
  # read in PLINK files --------------------------------
  path <- paste0(data_dir, "/", prefix, ".rds")
  bk_path <- paste0(data_dir, "/", prefix, ".bk")
  std_bk_path <- paste0(data_dir, "/std_", prefix, ".bk")
  sub_bk_path <- paste0(data_dir, "/subset_", prefix, ".bk")
  
  # check for overwrite: 
  if (file.exists(bk_path)){
    if (overwrite){
      # notify 
      cat("\nOverwriting existing files: ", prefix, ".bk/.rds\n",
          file = outfile, append = TRUE)
      
      if (!quiet){
        cat("\nOverwriting existing files: ", prefix, ".bk/.rds\n")
      }
      
      # overwrite existing files 
      system(paste0("rm ", bk_path))
      system(paste0("rm ", std_bk_path))
      system(paste0("rm ", sub_bk_path))
      system(paste0("rm ", path))
    } else {
      stop("\nThere are existing prefix.rds and prefix.bk files in the specified directory.  
           \nIf you want to overwrite these existing files, set 'overwrite = TRUE'. 
           \nOtherwise, choose a different prefix.")
    }
  }
  
  
  # create the RDS file first ------------------------
  cat("\nCreating ", prefix, ".rds\n", file = outfile, append = TRUE)
  if(!quiet){
    cat("\nCreating ", prefix, ".rds\n")
    
    # check for compressed files 
    if (gz){
      cat("\nUnzipping .gz files - this could take a second", file = outfile, append = TRUE)
      if (!quiet){cat("\nUnzipping .gz files - this could take a second")}
      system(paste0("gunzip -k ", file.path(data_dir, paste0(prefix, "*"))))
    }
    
    bigsnpr::snp_readBed(bedfile = paste0(data_dir, "/", prefix, ".bed"))
    obj <- bigsnpr::snp_attach(path)
  }
  
  # set object names --------------------------------
  obj$colnames <- obj$map$marker.ID
  
  if (id_var == "FID"){
    obj$rownames <- og_plink_ids <- as.character(obj$fam$family.ID)
  } else if (id_var == "IID") {
    obj$rownames <- og_plink_ids <- as.character(obj$fam$sample.ID)
  } else {
    stop("\nThe argument to id_var is misspecified. Must be one of 'IID' or 'FID'.")
  }

  chr <- obj$map$chromosome
  chr_range <- range(obj$map$chromosome)
  X <- obj$genotypes
  pos <- obj$map$physical.pos
  # save the dimensions of the *original* (pre-standardized) design matrix
  obj$n <- X$nrow
  obj$p <- X$ncol
  
  if(!quiet){
    cat("\nThere are ", obj$n, " observations and ",
        obj$p, " genomic features in the specified data files, representing chromosomes ",
        chr_range[1], " - ", chr_range[2])
  }
  
  # save these counts 
  counts <- bigstatsr::big_counts(X) # NB: this is a matrix 
  
  # chromosome check ---------------------------------
  # only consider SNPs on chromosomes 1-22

  if(chr_range[1] < 1 | chr_range[2] > 22){
    stop("\nplmmr only analyzes autosomes -- please remove variants on chromosomes outside 1-22.
         This can be done in PLINK 1.9; see the documentation in https://www.cog-genomics.org/plink/1.9/filter#chr")
  }
  
  #   # TODO: below is an idea for future development: 
  #   cat("\nplmmr only analyzes autosomes -- removing chromosomes outside 1-22")
  #   cat("\nplmmr only analyzes autosomes -- removing chromosomes outside 1-22",
  #       file = outfile, append = TRUE)
  #   
  #   original_dim <- dim(obj$genotypes)[2]
  #   chr_filtered <- bigsnpr::snp_subset(obj,
  #                                       ind.col = obj$map$chromosome %in% 1:22)
  #   obj <- bigsnpr::snp_attach(chr_filtered)
  #   new_dim <- dim(obj$genotypes)[2]
  #   
  #   cat("\nRemoved ", original_dim - new_dim, "SNPs that are outside of chromosomes 1-22.",
  #       file = outfile, append = TRUE)
  #   if(!quiet){
  #     cat("\nRemoved ", original_dim - new_dim, "SNPs that are outside of chromosomes 1-22.")
  #     
  #   }
  # }
  
  # TODO: figure out how to add a 'sexcheck' with bigsnpr functions
  # e.g., if sexcheck = TRUE, remove subjects with sex discrepancies
  
  # notify about missing (genotype) values ----------------------------
  na_idx <- counts[4,] > 0
  prop_na <- counts[4,]/nrow(X)
  
  cat("\nThere are a total of ", sum(na_idx), "SNPs with missing values",
      file = outfile, append = TRUE)
  cat("\nOf these, ", sum(prop_na > 0.5),
      " are missing in at least 50% of the samples",
      file = outfile, append = TRUE)
  if(!quiet){
    cat("\nThere are a total of ", sum(na_idx), "SNPs with missing values")
    cat("\nOf these, ", sum(prop_na > 0.5), " are missing in at least 50% of the samples")
  }
  
  # handle missing phenotypes ---------------------------------------
  # make missing phenotypes explicit (need both of the following because 
  # bigstatsr::big_cop() does not handle negative indices)
  complete_phen <- which(!(obj$fam$affection %in% na_phenotype_vals))
  na_phen <- which(obj$fam$affection %in% na_phenotype_vals)
  names(na_phen) <- obj$fam$sample.ID[na_phen]
  
  if (handle_missing_phen == 'prune'){
    if(!quiet){
      cat("\nWill prune out ", length(na_phen), " samples/observations with missing phenotype data.")
      # Note: the actual pruning happens in the 'subset' step 
    }
    
  } else if (handle_missing_phen == 'asis'){
    if(!quiet){
      cat("\nWill mark ", length(na_phen), " samples/observations as having missing phenotype data.")
    }
    obj$fam$affection[na_phen] <- NA_integer_
  } else {
    if(!quiet){
      cat("\nImputing phenotype data for ", length(na_phen), " samples/observations.")
    }
    obj$fam$affection[na_phen] <- switch(handle_missing_phen,
                                         median = median(obj$fam$affection[complete_phen]),
                                         mean = mean(obj$fam$affection[complete_phen]))
  }
  
  # imputation -------------------------------------------------
  if(!quiet & impute){
    # catch for misspellings
    if(!(impute_method %in% c('mode', 'random', 'mean0', 'mean2', 'xgboost'))){
      stop("\nImpute method is misspecified or misspelled. Please use one of the 
           \n5 options listed in the documentation.")
    }
    cat("\nImputing the missing (genotype) values using ", impute_method, " method\n")
  }
  
  if(impute){
    cat("\nImputing the missing values using ", impute_method, " method",
        file = outfile, append = TRUE)
    
    if(impute_method %in% c('mode', 'random', 'mean0', 'mean2')){ 
      
      obj$genotypes <- bigsnpr::snp_fastImputeSimple(Gna = X,
                                                 ncores = bigstatsr::nb_cores(),
                                                 method = impute_method) # dots can pass other args
      
    } else if (impute_method == "xgboost"){
      
      obj$genotypes <- bigsnpr::snp_fastImpute(Gna = X,
                                     ncores = bigstatsr::nb_cores(),
                                     infos.chr = chr,
                                     seed = as.numeric(Sys.Date()),
                                     ...) # dots can pass other args
      
      cat("\n ***************** NOTE ********************************
          \nAugust 2023: With the xgboost imputation method, there 
          \nhave been some issues (particularly on Mac OS) with warnings 
          \nthat appear saying 'NA or NaN values in the resulting correlation matrix.' 
          \nHowever, we (plmm authors) have 
          \n not seen missing values appear in the results -- the imputed data
          \n does not show any NA or NaN values, and models fit on these data run without issue. 
          \n We are actively investigating this warning message, and will
          \n make a note in a future release. If using xgboost, proceed with 
          \n caution and file an issue if you notice any problems downstream.
          \n ********************************************************")
      
      # save imputed values (NB: will overwrite obj$genotypes)
      obj$genotypes$code256 <- bigsnpr::CODE_IMPUTE_PRED
    }
    
    # save the imputed object 
    obj <- bigsnpr::snp_save(obj)
    
    cat("\nDone with imputation.",
        file = outfile, append = TRUE)
  }

  # add additional covariates -----------------------
  # first, set up some indices; even if no additional args are used, these NULL
  #   values are important for checks downstream
  non_gen <- NULL 
  ## covariates from .fam file ---------------------------
  if (!is.null(add_predictor_fam)) {
    if (!quiet) {
      cat("\nAdding predictors from .fam file.")
    }
    if (add_predictor_fam == "sex"){
      # add space for extra column
      obj$geno_plus_predictors <- bigstatsr::FBM(init = 0,
                                                 nrow = nrow(obj$fam), 
                                                 ncol = obj$genotypes$ncol + 1) 
      # fill in new matrix 
      bigstatsr::big_apply(obj$genotypes,
                           a.FUN = function(X, ind, res){
                             res[,1] <- obj$fam[add_predictor_fam]
                             res[,ind+1] <- X[,ind]
                           },
                           a.combine = cbind,
                           res = obj$geno_plus_predictors)
      
      # save non_gen: an index marking the first column as non-genomic predictor
      non_gen <- 1
      
    } # TODO: may add option to supply 6th column of .fam file here...
  }
  
  ## covariates from external file   ----------------------------------
  if (!is.null(add_predictor_ext)) {
    if (!quiet) {
      cat("\nAdding predictors from external data.")
    }
    if (is.vector(add_predictor_ext)) {
      ### vector case -------------------------------
      # make sure types match
      if (!is.numeric(add_predictor_ext)) {
        stop("\nThe vector supplied to the 'add_predictor_ext' argument must be numeric.")
      }
      names(add_predictor_ext) <- as.numeric(names(add_predictor_ext))
      
      if (var(add_predictor_ext) == 0) {
        stop("\nThe supplied argument to add_predictor_ext is constant (no variation).
             This would not be a meaningful predictor.")
      }
      
      # check for alignment 
      if (is.null(names(add_predictor_ext)) | 
          length(intersect(og_plink_ids, names(add_predictor_ext))) == 0) {
        stop("\nYou supplied an argument to 'add_predictor_ext', but the names of this
         vector either (a) do not exist or (b) do not align with either of the ID columns in the PLINK fam file.
         \nPlease create or align the names of this vector - alignment is essential for accurate analysis.")
      }
      
      add_predictor_ext <- align_famfile_ids(id_var = id_var,
                                             quiet = quiet,
                                             add_predictor = add_predictor_ext,
                                             og_plink_ids = og_plink_ids)
      
      # save non_gen: an index marking the first column as non-genomic predictor
      non_gen <- 1
      
      obj$geno_plus_predictors <- bigstatsr::FBM(init = 0,
                                                 nrow = nrow(obj$fam),
                                                 ncol = obj$genotypes$ncol + length(non_gen)) 
      # fill in new matrix 
      bigstatsr::big_apply(obj$genotypes,
                           a.FUN = function(X, ind, res){
                             res[,1:length(non_gen)] <- add_predictor_ext
                             res[,ind+length(non_gen)] <- X[,ind]
                           },
                           a.combine = cbind,
                           res = obj$geno_plus_predictors)
      
      # adjust colnames
      obj$colnames <- c(deparse(substitute(add_predictor_ext)), obj$colnames)
      
    } else if (is.matrix(add_predictor_ext) | is.data.frame(add_predictor_ext)) {
      ### matrix case --------------------------------
      if (is.data.frame(add_predictor_ext)) {
        add_predictor_ext <- as.matrix(add_predictor_ext)
      }
      # make sure types match
      if (!is.numeric(add_predictor_ext[,1])) {
        stop("\nThe matrix supplied to the 'add_predictor_ext' argument must have numeric values only.")
      }
      
      if (any(apply(add_predictor_ext, 2, var) == 0)) {
        stop("\nThe matrix supplied to the 'add_predictor_ext' argument has at least one
             constant column (a column that does not vary over the given samples).")
      }
      
      # check for alignment 
      if (is.null(rownames(add_predictor_ext)) | 
          length(intersect(og_plink_ids, rownames(add_predictor_ext))) == 0) {
        stop("\nYou supplied an argument to 'add_predictor_ext', but the row names of this
         matrix either (a) do not exist or (b) do not align with either of the ID columns in the PLINK fam file.
         \nPlease create or align the names of this matrix - alignment is essential for accurate analysis.")
      }
      
      add_predictor_ext <- align_famfile_ids(id_var = id_var,
                                             quiet = quiet,
                                             add_predictor = add_predictor_ext,
                                             og_plink_ids = og_plink_ids)

      # save non_gen: an index marking added columns as non-genomic predictors
      non_gen <- 1:ncol(add_predictor_ext)
      
      obj$geno_plus_predictors <- bigstatsr::FBM(init = 0,
                                                 nrow = nrow(obj$fam),
                                                 ncol = obj$genotypes$ncol + length(non_gen)) 
      # fill in new matrix 
      bigstatsr::big_apply(obj$genotypes,
                           a.FUN = function(X, ind, res){
                             res[,1:length(non_gen)] <- add_predictor_ext
                             res[,ind+length(non_gen)] <- X[,ind]
                           },
                           a.combine = cbind,
                           res = obj$geno_plus_predictors)
      
      # adjust colnames if applicable 
      if (!is.null(colnames(add_predictor_ext))){
        obj$colnames <- c(colnames(add_predictor_ext), obj$colnames)
      }
      
    }
  
  }
  
  # identify monomorphic SNPs --------------------------------
  
  
  # subsetting -----------------------------------------
  # goal here is to subset the features so that constant features (monomorphic SNPs) are not 
  # included in analysis
  # NB: this is also where we remove observations with missing phenotypes, if that was requested

  if (!quiet){
    cat("\nSubsetting data to exclude constant features (e.g., monomorphic SNPs)",
        file = outfile, append = TRUE)
  }
  
  sub_bk_extension <- paste0("subset_", prefix) 
  
  if (handle_missing_phen == "prune"){
    if ("geno_plus_predictors" %in% names(obj)) {
      new_counts <- bigstatsr::big_counts(obj$geno_plus_predictors,
                                          ind.row = complete_phen)
      ns <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      obj$subset_X <- bigstatsr::big_copy(obj$geno_plus_predictors,
                                          ind.row = complete_phen, # filters out rows with missing phenotypes
                                          ind.col = ns,
                                          backingfile = paste0(data_dir,"/", sub_bk_extension))
    } else {
      new_counts <- bigstatsr::big_counts(obj$genotypes,
                                          ind.row = complete_phen)
      ns <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      obj$subset_X <- bigstatsr::big_copy(obj$genotypes,
                                          ind.row = complete_phen, # filters out rows with missing phenotypes
                                          ind.col = ns,
                                          backingfile = paste0(data_dir,"/", sub_bk_extension))
    }
   
  } else {
    if ("geno_plus_predictors" %in% names(obj)) {
      new_counts <- bigstatsr::big_counts(obj$geno_plus_predictors)
      ns <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      obj$subset_X <- bigstatsr::big_copy(obj$geno_plus_predictors,
                                          ind.col = ns,
                                          backingfile = paste0(data_dir,"/", sub_bk_extension))
    } else {
      new_counts <- bigstatsr::big_counts(obj$genotypes)
      ns <- count_constant_features(new_counts, outfile = outfile, quiet = quiet)
      obj$subset_X <- bigstatsr::big_copy(obj$genotypes,
                                          ind.col = ns,
                                          backingfile = paste0(data_dir,"/", sub_bk_extension))
    }
  }
  
  # standardization ------------------------------------------------
  if (!quiet) {cat("\nColumn-standardizing the design matrix...")}

  # centering & scaling 
  scale_info <- bigstatsr::big_scale()(obj$subset_X)

  obj$std_X <- big_std(X = obj$subset_X,
                 center = scale_info$center,
                 scale = scale_info$scale) # leave ns = NULL; X is already subset
  
  std_bk_extension <- paste0("std_", prefix) 
  
  # now, save the indices of non-singular values in the *new* X; if additional 
  #   predictors have been added, these index values need to be updated! 
  if (is.null(non_gen)) {
    obj$ns <- ns
  } else {
    obj$ns <- c(non_gen, ns + length(non_gen))
  }
  
  
  # label return object ------------------------------------------------
  # naming these center and scale values so that I know they relate to the first
  # standardization; there will be another standardization after the rotation
  # in plmm_fit().
  obj$std_X_center <- scale_info$center
  obj$std_X_scale <- scale_info$scale
  obj$std_X_colnames <- obj$colnames[obj$ns]
  obj$std_X_rownames <- obj$rownames[complete_phen]
  obj$non_gen <- non_gen # save indices for non-genomic covariates
  obj$complete_phen <- complete_phen # save indices for which samples had complete phenotypes
  obj$id_var <- id_var # save ID variable - will need this downstream for analysis
  obj <- bigsnpr::snp_save(obj)
  
  
  if (!quiet){  
    cat("\nDone with standardization. File formatting in progress.",
        file = outfile, append = TRUE)
  }
  
  if(!quiet & impute){cat("\nDone with standardization. Processed files now saved as .rds object.")}
  close(log_con)
}



