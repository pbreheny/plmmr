#' A function to read in large data files as an FBM
#'
#' @param file          The file to be read in, without the filepath. This should be a file of numeric values, no header row!
#'                      Headers should be taken out of the file and supplied to the `col_names` argument in `plmm()`
#'                      Example: use `file = "myfile.txt"`, not `file = "~/mydirectory/myfile.txt"`
#' @param data_dir      The directory to the file
#' @param rds_dir       The directory where the user wants to create the '.rds' and '.bk' files
#'                      Defaults to `data_dir`
#' @param bk_filename   Optional string to name the backingfile that will be created for the output data. Defaults to `paste0(std_`. `prefix`).
#'                      **Note**: Do NOT include a `.bk` extension in the filename.
#' @param ind.col       Numeric vector of the columns to read in. Don't use negative indicies.
#'                      If you're not sure how many columns are in your file,
#'                      use something like `data.table::fread()` to examine the first row or two.
#' @param unpen      Numeric vector with indices of columns that have non-genomic information. **Note**: if you plan on
#'                      estimating the genomic relatedness among observations in your file, you
#'                      **must** include this information in order for the estimates to be correct.
#' @param outfile       Optional: the name (character string) of the prefix of the
#'                      logfile to be written. Defaults to 'process_delim', i.e. you will get 'process_delim.log' as the outfile.
#' @param overwrite     Optional: the name (character string) of the prefix of the logfile to be written.
#'                      Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#'                      **Note**: If there are multiple `.rds` files with names that start with "std_prefix_...", **this will error out**.
#'                      To protect users from accidentally deleting files with saved results, only one `.rds` file can be removed with this option.
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
#'  data_dir = find_example_data(parent = TRUE),
#'  rds_dir = temp_dir,
#'   ind.col = 2:2002)
#'
#' # colon2_rds <- readRDS(paste0(temp_dir, "std_colon2.rds"))
#' # str(colon2_rds)
#'
process_delim <- function(file,
                          data_dir,
                          rds_dir = data_dir,
                          bk_filename,
                          ind.col,
                          unpen= NULL,
                          outfile,
                          overwrite = FALSE,
                          quiet = FALSE){

  prefix <- unlist(strsplit(file, split = "\\."))[1]

  if(missing(bk_filename)){
    bk_filename <- paste0("std_", prefix)
  }

  # start log ------------------------------------------
  if(missing(outfile)){
    outfile = file.path(data_dir, "process_data")
  }

  logfile <- create_log(outfile = outfile)

  if(!quiet){
    cat("Logging to", logfile, "\n")
    cat("Preprocessing", prefix, "data:\n")
  }

  cat("Preprocessing", prefix, "data\n", file = logfile, append = TRUE)

  # read in data files --------------------------------
  X <- read_data_files(file = file,
                       data_dir = data_dir,
                       rds_dir = rds_dir,
                       ind.col = ind.col,
                       outfile = logfile,
                       overwrite = overwrite,
                       quiet = quiet)

  # note the original dimensions
  n <- nrow(X)
  p <- ncol(X)

  cat("There are", n, "observations and", p, "features in the specified data files.\n",
      file = logfile, append = TRUE)

  if (!quiet){
    cat("There are", n, "observations and", p,
        "features in the specified data files.\n")

  }

  # notify about missing values ---------------------------------
  colstats <- bigstatsr::big_colstats(X)
  na_idx <- is.na(colstats$sum) # logical index

  cat("There are a total of ", sum(na_idx), "features with missing values.
      At this time, plmmr::process_delim() does not impute missing values. We
      are working to develop this feature.
      For now, plmmr::process_delim() will drop all columns of X with NA values.\n",
      file = logfile, append = TRUE)
  if(!quiet){
    cat("There are a total of ", sum(na_idx), "features with missing values.\n
      At this time, plmmr::process_delim() does not impute missing values. We
      are working to develop this feature.\n
      For now, plmmr::process_delim() will drop all columns of X with NA values.\n")
  }

  # check for files to be overwritten---------------------------------
  if (overwrite){
    gc()

    # double check how much this will erase
    rds_to_remove <-  list.files(rds_dir, pattern=paste0('^std_.*.rds'), full.names=TRUE)
    if (length(rds_to_remove) > 1) {
      stop("You set overwrite=TRUE, but this looks like it will overwrite multiiple .rds files
      that have the 'std_prefix' pattern.
           To save you from overwriting anything important, I will not erase anything yet.
           Please move any .rds files with this file name pattern to another directory.")
    }
    file.remove(rds_to_remove)
    gc()

    list.files(rds_dir, pattern=paste0('^std_.*.bk'), full.names=TRUE) |>
      file.remove()

    gc()  # this is important!
  }

  # subsetting -----------------------------------------------------------------
  # goal: subset columns to remove constant features
  ns <- which(colstats$var > 1e-8)
  # ns = index marking which columns have nonzero, NON-MISSING variance
  # constant features and features with missing values are subset out!
  if (length(ns) < ncol(X)){
    subset_X <- bigstatsr::big_copy(X, ind.col = ns,
                                    backingfile = file.path(rds_dir, bk_filename))
  } else {
    subset_X <- X
  }

  # standardization ------------------------------------------------------------
  std_X_list <- standardize_fbm(subset_X = subset_X,
                                prefix = prefix,
                                rds_dir = rds_dir,
                                ns = ns,
                                unpen= non_gen,
                                outfile = logfile,
                                quiet = quiet)
  std_X_list$n <- n
  std_X_list$p <- p


  # cleanup --------------------------------------------------------------------
  # These steps remove intermediate rds/bk files created by the steps of the data management process
  list.files(rds_dir, pattern=paste0('^', prefix, '.*.rds'), full.names=TRUE) |>
    file.remove()
  gc() # this is important!
  list.files(rds_dir, pattern=paste0('^', prefix, '.*.bk'), full.names=TRUE) |>
    file.remove()
  gc() # this is important!
  list.files(rds_dir, pattern=paste0('^file.*.bk'), full.names=TRUE) |>
    file.remove()
  gc()

  cat("Processed files now saved to:",
      file.path(rds_dir, paste0("std_", prefix, ".rds")),
      "at",
      pretty_time(),
      file = logfile,
      append = TRUE)

  if(!quiet){cat("Processed files now saved to:",
                 file.path(rds_dir, paste0("std_", prefix, ".rds")),
                 "at",
                 pretty_time())}

  saveRDS(std_X_list, file.path(rds_dir, paste0("std_", prefix, ".rds")))
  return(file.path(rds_dir, paste0("std_", prefix, ".rds")))
}
