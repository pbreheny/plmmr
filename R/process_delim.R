#' A function to read in large data files as an FBM 
#'
#' @param file          The file to be read in, without the filepath. This should be a file of numeric values, no header row! 
#'                      Headers should be taken out of the file and supplied to the `col_names` argument in `plmm()`
#'                      Example: use `file = "myfile.txt"`, not `file = "~/mydirectory/myfile.txt"`
#' @param data_dir      The directory to the file
#' @param rds_dir       The directory where the user wants to create the '.rds' and '.bk' files
#'                      Defaults to `data_dir`
#' @param ind.col       Numeric vector of the columns to read in. Don't use negative indicies. 
#'                      If you're not sure how many columns are in your file, 
#'                      use something like `data.table::fread()` to examine the first row or two.
#' @param non_gen       Numeric vector with indices of columns that have non-genomic information. **Note**: if you plan on 
#'                      estimating the genomic relatedness among observations in your file, you 
#'                      **must** include this information in order for the estimates to be correct. 
#' @param outfile       Optional: the name (character string) of the prefix of the 
#'                      logfile to be written. Defaults to 'process_delim', i.e. you will get 'process_delim.log' as the outfile.
#' @param overwrite     Optional: the name (character string) of the prefix of the logfile to be written. 
#'                      Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param quiet         Logical: should the messages printed to the console be silenced? Defaults to FALSE.
#' 
#' @return Nothing is returned by this function, but (at least) two files are created in 
#' the location specified by `rds_dir`:
#' 
#' * 'std_prefix.rds': This is the `bigsnpr::bigSNP` object
#' that holds the PLINK data along with meta-data. See details for explanation of what 
#' is included in this meta-data
#' 
#' * 'std_prefix.bk': Created by the call to `standardize_fbm()`, this is the 
#' backingfile that stores the numeric data of the standardized design matrix `std_X`
#'  
#' @export
#'
#' @examples
#' temp_dir <- tempdir()
#' process_delim(file = "colon2.txt",
#'  data_dir = get_example_data(parent = TRUE),
#'  rds_dir = temp_dir,
#'   ind.col = 2:2002)
#'   
#' # colon2_rds <- readRDS(paste0(temp_dir, "std_colon2.rds"))
#' # str(colon2_rds)
#' 
process_delim <- function(file,
                      data_dir,
                      rds_dir = data_dir,
                      ind.col, 
                      non_gen = NULL,
                      outfile,
                      overwrite = FALSE,
                      quiet = FALSE){
  
  prefix <- unlist(strsplit(file, split = "\\."))[1]
  
  # start log ------------------------------------------
  if(missing(outfile)){
    outfile = paste0(rds_dir, "/process_data.log")
  } else {
    outfile = paste0(outfile, ".log")
  }
  
  log_con <- file(outfile)
  cat("### Processing files for PLMM ###", file = log_con)
  cat("\nLogging to ", outfile, file = outfile, append = TRUE)
  cat("\nPreprocessing", prefix, "data:", file = outfile, append = TRUE)
  
  if(!quiet){
    cat("\nLogging to", outfile)
    cat("\nPreprocessing", prefix, "data:")
  }
  
  # read in data files --------------------------------
 X <- read_data_files(file, data_dir, rds_dir, ind.col, outfile, overwrite, quiet)
  # note the original dimensions
  n <- nrow(X)
  p <- ncol(X)
  
  # notify about missing values ---------------------------------
  colstats <- bigstatsr::big_colstats(X)
  na_idx <- is.na(colstats$sum) # logical index
  
  cat("\nThere are a total of ", sum(na_idx), "features with missing values. 
      \nAt this time, plmmr::process_delim() does not impute missing values. We 
      are working to develop this feature. 
      \nFor now, plmmr::process_delim() will drop all columns of X with NA values.",
      file = outfile, append = TRUE)
  if(!quiet){
    cat("\nThere are a total of ", sum(na_idx), "features with missing values. 
      \nAt this time, plmmr::process_delim() does not impute missing values. We 
      are working to develop this feature. 
      \nFor now, plmmr::process_delim() will drop all columns of X with NA values.")
  }
 
  # subsetting -----------------------------------------------------------------
  # goal: subset columns to remove constant features 
  ns <- which(colstats$var > 1e-8) 
  # ns = index marking which columns have nonzero, NON-MISSING variance 
  # constant features and features with missing values are subset out! 
  if (length(ns) < ncol(X)){
    subset_X <- bigstatsr::big_copy(X,
                                    ind.col = ns)
  } else {
    subset_X <- X
  }
  
  # standardization ------------------------------------------------------------
  std_X_list <- standardize_fbm(subset_X, prefix, rds_dir, ns, non_gen, 
                           outfile, quiet)
  std_X_list$n <- n
  std_X_list$p <- p
  
  saveRDS(std_X_list, file.path(rds_dir, paste0("std_", prefix, ".rds")))
  
  # cleanup --------------------------------------------------------------------
    file.remove(paste0(rds_dir, "/", prefix, ".rds"))
    file.remove(paste0(rds_dir, "/", prefix, ".bk"))

  if(!quiet){cat("\nDone with standardization. 
                 Processed files now saved as .rds object.")}
  close(log_con)
  
  
}
