#' A helper function to standardize an FBM
#'
#' @param obj 
#' @param prefix 
#' @param non_gen 
#' @param complete_phen 
#' @param id_var        String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID". 
#' @param outfile       Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet 
#'
#' @return
#' @keywords internal
#'
standardize_fbm <- function(obj, prefix, non_gen, complete_phen, id_var,
                            outfile, quiet){
  # standardization ------------------------------------------------
  if (!quiet) {cat("\nColumn-standardizing the design matrix...")}

  # centering & scaling 
  scale_info <- bigstatsr::big_scale()(obj$subset_X)
  
  obj$std_X <- big_std(X = obj$subset_X,
                       center = scale_info$center,
                       scale = scale_info$scale) # leave ns = NULL; X is already subset
  
  std_bk_extension <- paste0("std_", prefix) 
  
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
  
  return(obj)
}