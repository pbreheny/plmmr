#' a function to create a design for PLMM modeling
#' @param data_file     For filebacked data (data from `process_plink()` or `process_delim()`), this is the
#'                      filepath to the processed data. Defaults to NULL (does not apply for in-memory data)
#' @param rds_dir       For filebacked data, this is the filepath to the directory/folder where you want the design to be saved.
#'                      **Note**: do not include/append the name you want for the to-be-created file -- the name is the argument `new_file`,
#'                      passed to `create_design_filebacked()`.
#' @param X             For in-memory data, this is the design matrix.
#' @param outcome_col   For in-memory data, this is the numeric vector representing the outcome.
#' @param ...           Additional arguments to pass to `create_design_filebacked()` or `create_design_in_memory()`.
#'                      See the documentation for those helper functions for details.
#'
#' @return An object of class `plmm_design`
#'
#' @export
#'
#' @details
#' This function is a wrapper for the other `create_design...()` inner functions; all arguments
#' included here are passed along to the `create_design...()` inner function that
#' matches the type of the data being supplied.
#'
#' @examples
#' ## Example 1: delimited data ##
#' # process delimited data
#' temp_dir <- tempdir()
#' colon_dat <- process_delim(data_file = "colon2.txt",
#'  data_dir = find_example_data(parent = TRUE), overwrite = TRUE,
#'  rds_dir = temp_dir, rds_prefix = "processed_colon2", sep = "\t", header = TRUE)
#'
#' # prepare outcome data
#' colon_outcome <- read.delim(find_example_data(path = "colon2_outcome.txt"))
#'
#' # create a design
#' colon_design <- create_design(data_file = colon_dat, rds_dir = temp_dir, new_file = "std_colon2",
#' add_outcome = colon_outcome, outcome_id = "ID", outcome_col = "y",
#' overwrite = TRUE, logfile = "test.log")
#'
#' # look at the results
#' colon_rds <- readRDS(colon_design)
#' str(colon_rds)
#'
#' ## Example 2: matrix data in-memory
#' admix_design <- create_design(X = admix$X, outcome_col = admix$y, unpen = "Snp1")
#'
create_design <- function(data_file = NULL,
                          rds_dir = NULL,
                          X = NULL,
                          outcome_col = NULL,
                          ...) {

  if (is.null(data_file)) { # case 1: in-memory matrix
    processed_matrix = create_design_in_memory(X,
                                               outcome_col,
                                               ...)
  } else { # case 2: filebacked data
    obj <- readRDS(data_file)
    switch (class(obj),
            processed_plink = create_design_filebacked(data_file = data_file,
                                                       rds_dir = rds_dir,
                                                       obj = obj,
                                                       outcome_col = outcome_col,
                                                       ...),
            processed_delim = create_design_filebacked(data_file = data_file,
                                                       rds_dir = rds_dir,
                                                       obj = obj,
                                                       outcome_col = outcome_col,
                                                       ...)
    )
  }
}