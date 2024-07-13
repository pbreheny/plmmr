#' Read in processed data
#' This function is intended to be called after either `process_plink()` or `process_delim()` has been called once.
#'
#' @param path The file path to the RDS object containing the processed data. Do not add the '.rds' extension to the path.
#' @param returnX Logical: Should the design matrix be returned as a numeric matrix that will be stored in memory. By default, this will be FALSE if the object sizes exceeds 100 Mb.
#' @param trace Logical: Should trace messages be shown? Default is TRUE.
#' @returns A list with these components:
#'  * std_X, the column-standardized design matrix as either (1) a numeric matrix or (2) a filebacked matrix (FBM). See `bigstatsr::FBM()` and `bigsnpr::bigSnp-class` documentation for details.
#'  * fam, a data frame containing the pedigree information (like a .fam file in PLINK)
#'  * map, a data frame containing the feature information (like a .bim file in PLINK)
#'  * ns: A vector indicating the which columns of X contain nonsingular features (i.e., features with variance != 0.
#'  * center: A vector of values for centering each column in X
#'  * scale: A vector of values for scaling each column in X
#'
#' @export
#'
#' @examples
#' \donttest{
#' temp_dir <- paste0(tempdir()) # using a temporary directory here
#' process_plink(data_dir = find_example_data(parent = TRUE), # reads data that ships with plmmr
#'               rds_dir = temp_dir,
#'               prefix = "penncath_lite",
#'               gz = TRUE,
#'               outfile = "process_penncath",
#'               overwrite = TRUE,
#'               impute_method = "mode")
#'
#'  pen <- get_data(file.path(temp_dir, "std_penncath_lite"))
#'  str(pen)
#' }
#'
#' @details
#' The .rds object should have an 'std_X' element - this is what will be used as the design matrix for analysis. This design matrix should *not* include an intercept column (this will be added later in `plmm_fit`()).
#'
#' In the returned list, the `fam` data will be sorted by family and by individual, as in `dplyr::arrange(family.ID, sample.ID)`.
#' The rows of `X` will be sorted to align in the same order as in `fam`, where rownames of `X` will be sample ID.
#'
#'
get_data <- function(path, returnX, trace = TRUE){

 path <- check_for_file_extension(path)

  rds <- paste0(path, ".rds")
  bk <- paste0(path, ".bk") # .bk will be present if RDS was created with bigsnpr methods
  obj <- readRDS(rds)

  # return data in a tractable format
  if (missing(returnX)) {
    if (utils::object.size(obj$std_X) > 1e8) {
      warning("Due to the large size of X (>100 Mb), X has been returned as a file-backed matrix\n.
              To turn this message off, explicitly specify fbm=TRUE or fbm=FALSE)\n.")
      returnX <- FALSE
    } else {
      # if it fits, it ships
      returnX <- TRUE
    }
  }

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
