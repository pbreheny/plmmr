#' Read in processed data
#' This function is intended to be called after `process_plink()` has been called once. 
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
#' \dontrun{
#' pen <- get_data(path = "../tmp_files/penncath_lite", trace = TRUE)
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
  
  rds <- paste0(path, ".rds")
  bk <- paste0(path, ".bk") # .bk will be present if RDS was created with bigsnpr methods 
  
  if(file.exists(bk)){
    obj <- bigsnpr::snp_attach(rdsfile = rds)
  } else {
    obj <- readRDS(rds)
  }
  
  # return data in a tractable format 
  if (missing(returnX)) {
    if (utils::object.size(obj$std_X) > 1e8) {
      warning("\nDue to the large size of X (>100 Mb), X has been returned as a file-backed matrix (FBM).
              \nTo turn this message off, explicitly specify fbm=TRUE or fbm=FALSE).")
      returnX <- FALSE
    } else {
      # if it fits, it ships 
      returnX <- TRUE
    }
  }

  if(returnX){
    # get std_X as a matrix 
    std_X <- obj$std_X[,]
    
    # set row names 
    if(obj$id_var == "IID"){row_names <- as.character(obj$fam$sample.ID[obj$complete_phen])}
    if(obj$id_var == "FID"){row_names <- as.character(obj$fam$family.ID[obj$complete_phen])}
    
    dimnames(std_X) <- list(obj$std_X_rownames,
                            obj$std_X_colnames)
    
    cat("\nReminder: the X that is returned here is column-standardized, with constant features removed.
        \nA copy of the original data is available via the 'genotypes' matrix in the .rds object")
    return(list(n = obj$n,
                p = obj$p,
                X_colnames = obj$colnames,
                X_rownames = obj$rownames,
                std_X = std_X,
                std_X_center = obj$std_X_center,
                std_X_scale = obj$std_X_scale,
                fam = obj$fam,
                map = obj$map,
                ns = obj$ns,
                complete_phen = obj$complete_phen,
                id_var = obj$id_var,
                non_gen = obj$non_gen))
  } else {
    cat("Note: The design matrix is being returned as a file-backed matrix (FBM) -- see bigstatsr::FBM() for details.")
    
    cat("\nReminder: the X that is returned here is column-standardized.
        \nA copy of the original data is available via the 'genotypes' matrix in the .rds object")

    return(list(n = obj$n,
                p = obj$p,
                X_colnames = obj$colnames,
                X_rownames = obj$rownames,
                std_X = obj$std_X,
                std_X_center = obj$std_X_center,
                std_X_scale = obj$std_X_scale,
                std_X_rownames = obj$std_X_rownames,
                std_X_colnames = obj$std_X_colnames,
                fam = obj$fam,
                map = obj$map,
                ns = obj$ns,
                phen = obj$fam$affection[obj$complete_phen],
                complete_phen = obj$complete_phen,
                id_var = obj$id_var,
                non_gen = obj$non_gen))
  }
  
}