#' A helper function to add predictors to a filebacked matrix of data
#'
#' @param obj               A `bigSNP` object
#' @param add_predictor Optional: add additional covariates/predictors/features from an external file (i.e., not a PLINK file).
#' @param id_var            String specifying which column of the PLINK `.fam` file has the unique sample identifiers.
#' @param rds_dir           The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`(from `process_plink()` call)
#' @param quiet             Logical: should messages be printed to the console? Defaults to FALSE (which leaves the print messages on...)
#'
#' @return A list of 2 components:
#' * 'obj' - a `bigSNP` object with an added element representing the matrix that includes the additional predictors as the first few columns
#' * 'non_gen' - an integer vector that ranges from 1 to the number of added predictors. Example: if 2 predictors are added, non_gen = 1:2
#' @keywords internal
#'
add_predictors <- function(obj,
                           add_predictor,
                           id_var,
                           rds_dir,
                           quiet){


  # add additional covariates -----------------------
  # first, set up some indices; even if no additional args are used, these NULL
  #   values are important for checks downstream
  non_gen <- NULL

  if (!quiet) {
    cat("Adding predictors from external data.\n")
  }
  if (is.data.frame(add_predictor)) {
    add_predictor <- as.matrix(add_predictor)
  }
  # make sure types match
  if (!is.numeric(add_predictor[,1])) {
    stop("\nThe matrix supplied to the 'add_predictor' argument must have numeric values only.")
  }

  if (any(apply(add_predictor, 2, var) <1e-4)) {
    stop("\nThe matrix supplied to the 'add_predictor' argument has at least one
             constant column (a column that does not vary over the given samples).")
  }

  if (!quiet) cat('Aligning IDs between fam and predictor files\n')

  # save non_gen: an index marking added columns as non-genomic predictors
  non_gen <- 1:ncol(add_predictor)


  design_matrix <- big.matrix(nrow = nrow(obj$X),
                              ncol = ncol(obj$X) + length(non_gen),
                              type = 'double',
                              backingfile = "unstd_design_matrix.bk",
                              backingpath = rds_dir,
                              descriptorfile = "unstd_design_matrix.desc")

  design_matrix <- big_cbind(A = add_predictor,
                             B = obj$X,
                             C = design_matrix,
                             quiet = quiet)

  ret <- list(design_matrix = design_matrix, non_gen = non_gen)
  # adjust colnames if applicable
  if (!is.null(colnames(add_predictor))){
    ret$colnames <- c(colnames(add_predictor), obj$colnames)


    return(ret)

  }
}
