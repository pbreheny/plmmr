#' Read in processed data
#' This function is intended to be called after either `process_plink()` or `process_delim()` has been called once.
#'
#' @param path The file path to the RDS object containing the processed data. Do not add the '.rds' extension to the path.
#' @param returnX Logical: Should the design matrix be returned as a numeric matrix that will be stored in memory. By default, this will be FALSE.
#' @param trace Logical: Should trace messages be shown? Default is TRUE.
#'
#' @returns A list with these components:
#'  * std_X, the column-standardized design matrix as either (1) a numeric matrix or (2) a filebacked matrix (FBM). See `bigstatsr::FBM()` and `bigsnpr::bigSnp-class` documentation for details.
#'  * (if PLINK data) fam, a data frame containing the pedigree information (like a .fam file in PLINK)
#'  * (if PLINK data) map, a data frame containing the feature information (like a .bim file in PLINK)
#'  * ns: A vector indicating the which columns of X contain nonsingular features (i.e., features with variance != 0.
#'  * center: A vector of values for centering each column in X
#'  * scale: A vector of values for scaling each column in X
#'
#' @keywords internal
#'
#'
get_data <- function(path, returnX = FALSE, trace = TRUE){

 path <- check_for_file_extension(path)

  rds <- paste0(path, ".rds")
  bk <- paste0(path, ".bk") # .bk will be present if RDS was created with bigsnpr methods
  obj <- readRDS(rds)

  # attach std_X
  std_X_bm <- attach.big.matrix(obj$std_X)

  if(returnX){
    obj$std_X <- std_X_bm[,]
  if (trace){
    cat("Reminder: the X that is returned here is column-standardized, with constant features removed.\n")
  }
  } else {
    if (trace){
      cat("Note: The design matrix is being returned as a file-backed big.matrix object -- see bigmemory::big.matrix() documentation for details.\n")
      cat("Reminder: the X that is returned here is column-standardized\n")
    }
    obj$std_X <- std_X_bm
  }

  return(obj)

}
