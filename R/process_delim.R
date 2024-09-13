#' A function to read in large data files as an FBM
#'
#' @param data_dir      The directory to the file.
#' @param data_file     The file to be read in, without the filepath. This should be a file of numeric values, no header row!
#'                      Headers should be taken out of the file and supplied to the `col_names` argument in `plmm()`
#'                      Example: use `data_file = "myfile.txt"`, not `data_file = "~/mydirectory/myfile.txt"`
#' @param feature_id    A string specifying the column in the data X (the feature data) with the row IDs (e.g., identifiers for each row/sample/participant/, etc.). No duplicates allowed.
#' @param rds_dir       The directory where the user wants to create the '.rds' and '.bk' files
#'                      Defaults to `data_dir`
#' @param rds_prefix    String specifying the user's preferred filename for the to-be-created .rds file (will be create insie `rds_dir` folder)
#'                      Note: 'rds_prefix' cannot be the same as 'data_prefix'
#' @param id_var        String specifying which column of the `data_file` has the unique sample identifiers.
#' @param outfile       Optional: the name (character string) of the prefix of the
#'                      logfile to be written. Defaults to 'process_delim', i.e. you will get 'process_delim.log' as the outfile.
#' @param overwrite     Optional: the name (character string) of the prefix of the logfile to be written.
#'                      Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#'                      **Note**: If there are multiple `.rds` files with names that start with "std_prefix_...", **this will error out**.
#'                      To protect users from accidentally deleting files with saved results, only one `.rds` file can be removed with this option.
#' @param quiet         Logical: should the messages printed to the console be silenced? Defaults to FALSE.
#' @param ...           Optional: other arguments to be passed to `bigmemory::read.big.matrix()`. Note: 'sep' is an option to pass here, as is 'header'.
#' @return The file path to the newly created '.rds' file
#'
#' @export
#'
#' @examples
#' temp_dir <- tempdir()
#' colon_dat <- process_delim(data_file = "colon2.txt",
#'  data_dir = find_example_data(parent = TRUE), overwrite = TRUE,
#'  rds_dir = temp_dir, rds_prefix = "processed_colon2", sep = "\t", header = TRUE)
#'
#' colon2 <- readRDS(colon_dat)
#' str(colon2)
#'
process_delim <- function(data_dir,
                          data_file,
                          feature_id,
                          rds_dir = data_dir,
                          rds_prefix,
                          logfile = NULL,
                          col_ind,
                          quiet = FALSE,
                          overwrite = FALSE,
                          ...){

  prefix <- unlist(strsplit(data_file, split = "\\."))[1]

  if (identical(rds_prefix, prefix)) stop("rds_prefix cannot be the same as data_prefix. You need to change your choice of argument to rds_prefix.\n")

  # start log ------------------------------------------
  if(!is.null(logfile)){
    logfile <- create_log(file.path(rds_dir, logfile))
    cat("\nLogging to", logfile)
    cat("Preprocessing", prefix, "data\n", file = logfile, append = TRUE)
  } else {
    # TODO: change this default so that there is an option to turn off log files
    logfile <- tempfile()
  }

  # read in data files --------------------------------
  X <- read_data_files(data_file = data_file,
                       data_dir = data_dir,
                       rds_dir = rds_dir,
                       rds_prefix = rds_prefix,
                       outfile = logfile,
                       overwrite = overwrite,
                       quiet = quiet,
                       ...)

  cat("There are", nrow(X), "observations and", ncol(X), "features in the specified data files.\n",
      file = logfile, append = TRUE)

  if (!quiet){
    cat("There are", nrow(X), "observations and", ncol(X),
        "features in the specified data files.\n")
  }

  # notify about missing values ---------------------------------
  if(!quiet){
    cat("At this time, plmmr::process_delim() does not not handle missing values in delimited data.
      Please make sure you have addressed missingness before you proceed.\n")
  }

  # create return object --------------------------------------------
  ret <- structure(list(X = bigmemory::describe(X),
              # save original dimensions
              n = nrow(X),
              p = ncol(X)), class = "processed_delim")

  rds_filename <- paste0(rds_prefix, ".rds")
  saveRDS(ret, file = file.path(rds_dir, rds_filename))

  # cleanup ----------------------------------------------------------
  # These steps remove intermediate rds/bk files created by the steps of the data management process
  list.files(rds_dir, pattern=paste0('^file.*.bk'), full.names=TRUE) |>
    file.remove()
  list.files(rds_dir, pattern = paste0('^', prefix, ".rds"), full.names = TRUE) |>
    file.remove()
  list.files(rds_dir, pattern = paste0('^', prefix, ".bk"), full.names = TRUE) |>
    file.remove()

  gc()

  if(!quiet){cat("\nprocess_plink() completed \nProcessed files now saved as",
                 file.path(rds_dir, rds_filename))}

  cat("\nprocess_plink() completed. \nProcessed files now saved as",
      file.path(rds_dir, rds_filename),
      "at", pretty_time(),
      file = logfile, append = TRUE)

  return(file.path(rds_dir, rds_filename))


}
