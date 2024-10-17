#' a function to create a design for PLMM modeling
#' @param data_file     For **filebacked data** (data from `process_plink()` or `process_delim()`), this is the
#'                      filepath to the processed data. Defaults to NULL (this argument does not apply for in-memory data).
#' @param rds_dir       For **filebacked data**, this is the filepath to the directory/folder where you want the design to be saved.
#'                      **Note**: do not include/append the name you want for the to-be-created file -- the name is the argument `new_file`,
#'                      passed to `create_design_filebacked()`. Defaults to NULL (this argument does not apply for in-memory data).
#' @param X             For **in-memory data (data in a matrix or data frame)**, this is the design matrix. Defaults to NULL (this argument does not apply for filebacked data).
#' @param y             For **in-memory data**, this is the numeric vector representing the outcome. Defaults to NULL (this argument does not apply for filebacked data).
#' @param ...           Additional arguments to pass to `create_design_filebacked()` or `create_design_in_memory()`.
#'                      See the documentation for those helper functions for details.
#'
#' @return A filepath to an object of class `plmm_design`, which is a named list with the design matrix,
#'  outcome, penalty factor vector, and other details needed for fitting a model. This list is stored as an .rds
#'  file for filebacked data, so in the filebacked case a string with the path to that file is returned. For in-memory data,
#'  the list itself is returned.
#'
#' @export
#'
#' @details
#' This function is a wrapper for the other `create_design...()` inner functions; all arguments
#' included here are passed along to the `create_design...()` inner function that
#' matches the type of the data being supplied. Note which arguments are optional
#' and which ones are not.
#'
#' Additional arguments for **all filebacked** data:
#'
#'    - **new_file**                User-specified filename (*without .bk/.rds extension*) for the to-be-created .rds/.bk files. Must be different from any existing .rds/.bk files in the same folder.
#'
#'    - **feature_id**              Optional: A string specifying the column in the data X (the feature data) with the row IDs (e.g., identifiers for each row/sample/participant/, etc.). No duplicates allowed.
#'                                  - for PLINK data: a string specifying an ID column of the PLINK `.fam` file. Options are "IID" (default) and "FID"
#'                                  - for all other filebacked data: a character vector of unique identifiers (IDs) for each row of the feature data (i.e., the data processed with `process_delim()`)
#'                                  - if left NULL (default), X is assumed to have the same row-order as add_outcome.
#'                                  **Note**: if this assumption is made in error, calculations downstream will be incorrect. Pay close attention here.
#'
#'    - **add_outcome**             A data frame or matrix with two columns: and ID column and a column with the outcome value (to be used as 'y' in the final design). IDs must be characters, outcome must be numeric.
#'
#'    - **outcome_id**              A string specifying the name of the ID column in 'add_outcome'
#'
#'    - **outcome_col**             A string specifying the name of the phenotype column in 'add_outcome'
#'
#'    - **na_outcome_vals**        Optional: a vector of numeric values used to code NA values in the outcome. Defaults to `c(-9, NA_integer)` (the -9 matches PLINK conventions).
#'
#'    - **overwrite**              Optional: logical - should existing .rds files be overwritten? Defaults to FALSE.
#'
#'    - **logfile**                Optional: name of the '.log' file to be written -- **Note:** do not append a `.log` to the filename; this is done automatically.
#'
#'    - **quiet**                  Optional: logical - should messages to be printed to the console be silenced? Defaults to FALSE
#'
#' Additional arguments specific to **PLINK** data:
#'    - **add_predictor**           Optional (for PLINK data only): a matrix or data frame to be used for adding additional **unpenalized** covariates/predictors/features from an external file (i.e., not a PLINK file).
#'                                This matrix must have one column that is an ID column; all other columns aside the ID will be used as covariates in the design matrix. Columns must be named.
#'
#'    - **predictor_id**            Optional (for PLINK data only): A string specifying the name of the column in 'add_predictor' with sample IDs. Required if 'add_predictor' is supplied.
#'                                The names will be used to subset and align this external covariate with the supplied PLINK data.
#'
#' Additional arguments specific to **delimited file** data:
#'    - **unpen**         Optional: an character vector with the names of columns to mark as unpenalized (i.e., these features would always be included in a model).
#'                      **Note**: if you choose to use this option, your delimited file **must** have column names.
#'
#' Additional arguments for **in-memory** data:
#'
#'    - **y**             A numeric vector representing the outcome for the model.
#'                      **Note**: it is the responsibility of the user to ensure that the y and X have the same row order!
#'    - **unpen**         Optional: an character vector with the names of columns to mark as unpenalized (i.e., these features would always be included in a model).
#'                      **Note**: if you choose to use this option, X must have column names.
#'
#' @examples
#'
#' ## Example 1: matrix data in-memory ##
#' admix_design <- create_design(X = admix$X, y = admix$y, unpen = "Snp1")
#'
#' ## Example 2: delimited data ##
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
#' add_outcome = colon_outcome, outcome_id = "ID", outcome_col = "y", unpen = "sex",
#' overwrite = TRUE, logfile = "test.log")
#'
#' # look at the results
#' colon_rds <- readRDS(colon_design)
#' str(colon_rds)
#'
#' ## Example 3: PLINK data ##
#' \donttest{
#' # process PLINK data
#' temp_dir <- tempdir()
#' unzip_example_data(outdir = temp_dir)
#'
#' plink_data <- process_plink(data_dir = temp_dir,
#'   data_prefix = "penncath_lite",
#'   rds_dir = temp_dir,
#'   rds_prefix = "imputed_penncath_lite",
#'   # imputing the mode to address missing values
#'   impute_method = "mode",
#'   # overwrite existing files in temp_dir
#'   # (you can turn this feature off if you need to)
#'   overwrite = TRUE,
#'   # turning off parallelization - leaving this on causes problems knitting this vignette
#'   parallel = FALSE)
#'
#' # get outcome data
#' penncath_pheno <- read.csv(find_example_data(path = 'penncath_clinical.csv'))
#'
#' outcome <- data.frame(FamID = as.character(penncath_pheno$FamID),
#'                   CAD = penncath_pheno$CAD)
#'
#' unpen_predictors <- data.frame(FamID = as.character(penncath_pheno$FamID),
#'                                sex = penncath_pheno$sex,
#'                                age = penncath_pheno$age)
#'
#'
#' # create design where sex and age are always included in the model
#' pen_design <- create_design(data_file = plink_data,
#'   feature_id = "FID",
#'   rds_dir = temp_dir,
#'   new_file = "std_penncath_lite",
#'   add_outcome = outcome,
#'   outcome_id = "FamID",
#'   outcome_col = "CAD",
#'   add_predictor = unpen_predictors,
#'   predictor_id = "FamID",
#'   logfile = "design",
#'   # again, overwrite if needed; use with caution
#'   overwrite = TRUE)
#'
#' # examine the design - notice the components of this object
#' pen_design_rds <- readRDS(pen_design)
#'
#'}
#'
#'
create_design <- function(data_file = NULL,
                          rds_dir = NULL,
                          X = NULL,
                          y = NULL,
                          ...) {

  if (is.null(data_file)) { # case 1: in-memory matrix
    processed_matrix = create_design_in_memory(X,
                                               y,
                                               ...)
  } else { # case 2: filebacked data
    obj <- readRDS(data_file)
    switch (class(obj),
            processed_plink = create_design_filebacked(data_file = data_file,
                                                       rds_dir = rds_dir,
                                                       obj = obj,
                                                       ...),

            processed_delim = create_design_filebacked(data_file = data_file,
                                                       rds_dir = rds_dir,
                                                       obj = obj,
                                                       ...)
    )
  }
}
