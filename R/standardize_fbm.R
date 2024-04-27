#' A helper function to standardize an `FBM`
#'
#' @param subset_X An FBM 
#' @param prefix   The data file name without the extension
#' @param rds_dir  The directory where the user wants to create the '.rds' and '.bk' files
#' @param ns       Integer vector with indices of non-singular columns of `subset_X`
#' @param non_gen  Integer vector of indices with non-genomic covariates 
#' @param outfile  Optional: the name (character string) of the prefix of the logfile to be written. 
#' @param quiet    Logical: should the messages printed to the console be silenced? Defaults to FALSE.
#' 
#' @returns A list that includes:
#' * 'std_X' - the standardized FBM
#' * 'center' - a vector with  the values used to center the columns 
#' * 'scale' - a vector with the values used to scale the columns 
#' 
#' @keywords internal
#' 
standardize_fbm <- function(subset_X, prefix, rds_dir, ns, non_gen,
                           outfile, quiet){
  # standardization ------------------------------------------------
  if (!quiet) {cat("\nColumn-standardizing the design matrix...")}

  # centering & scaling 
  scale_info <- bigstatsr::big_scale()(subset_X)
  
  std_X <- big_std(X = subset_X,
                       std_bk_extension = file.path(rds_dir, paste0("std_", prefix)),
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
    ns = ns,
    non_gen = non_gen) # save indices for non-genomic covariates
  
  if (!quiet){  
    cat("\nDone with standardization. File formatting in progress.",
        file = outfile, append = TRUE)
  }
  
  return(ret)
}