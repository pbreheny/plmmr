#' A helper function to standardize an `FBM`

standardize_fbm <- function(subset_X, prefix, rds_dir, ns, non_gen, id_var,
                           outfile, quiet){
  # standardization ------------------------------------------------
  if (!quiet) {cat("\nColumn-standardizing the design matrix...")}

  # centering & scaling 
  scale_info <- bigstatsr::big_scale()(subset_X)
  
  std_X <- big_std(X = subset_X,
                       std_bk_extension = paste0(rds_dir, "/std_", prefix),
                       center = scale_info$center,
                       scale = scale_info$scale) # leave ns = NULL; X is already subset
  
  # label return object --------------------------------------------------------
  # naming these center and scale values so that I know they relate to the first
  # standardization; there will be another standardization after the rotation
  # in plmm_fit().
  ret <- list(
    std_X = std_X,
    std_X_center = scale_info$center,
    std_X_scale = scale_info$scale,
    non_gen = non_gen, # save indices for non-genomic covariates
    id_var = id_var # save ID variable - will need this downstream for analysis
  )
  
  if (!quiet){  
    cat("\nDone with standardization. File formatting in progress.",
        file = outfile, append = TRUE)
  }
  
  return(ret)
}