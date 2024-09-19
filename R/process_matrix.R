#' A function to process data from an in-memory matrix
#'
#' @param X           A numeric matrix of data. Note: rows should be observations (e.g., samples); columns should be features (e.g., predictors, SNPs, variants, etc.)
#' @param rds_dir     String specifying a filepath for a directory where output will be created
#' @param rds_prefix  String specifying the prefix of the '.rds' file to be created
#' @param quiet       Logical; if true, a messsage will be printed indicating the location of the processed data.
#'
#' @returns           The complete filepath (including '.rds' extension) of the newly created object
#' @export
#'
#' @examples
#' process_matrix(X = admix$X, rds_dir = tempdir(), rds_prefix = 'admix')
process_matrix <- function(X,
                           rds_dir,
                           rds_prefix,
                           quiet = FALSE){
  # check for missingness -------------------------------
  if (any(is.na(X))) {
      stop('The data supplied to "X" contains NA values
         Please take actions to address missing data.
         There can be no NA values in "X"\n')
  }

  # format output -------------------------------------

  # rows must be named
  if (is.null(colnames(X))) {
    colnames(X) <- paste0('feat', 1:ncol(X))
  }

  bm_X <- bigmemory::as.big.matrix(X, type = 'double')
  rds_filepath <- file.path(rds_dir, rds_prefix)
  rds <- list(X = bigmemory::describe(bm_X),
              n = nrow(X),
              p = ncol(X))
  saveRDS(rds, file = rds_filepath)
  if (!quiet) {
   cat("Output saved to", rds_filepath)
  }

  return(rds_filepath)
}