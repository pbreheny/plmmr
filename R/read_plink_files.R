#' A function to read in PLINK files using `bigsnpr` methods
#'
#' @param data_dir        The path to the bed/bim/fam data files, *without* a trailing "/" (e.g., use `data_dir = '~/my_dir'`, **not** `data_dir = '~/my_dir/'`)
#' @param data_prefix     The prefix (as a character string) of the bed/fam data files (e.g., `prefix = 'mydata'`)
#' @param rds_dir         The path to the directory in which you want to create the new `.rds` and `.bk` files. Defaults to `data_dir`
#' @param rds_prefix      String specifying the user's preferred filename for the to-be-created `.rds` file (will be create inside `rds_dir` folder). If no rds_prefix is provided, the processed data files will be returned in memory.
#'                        Note: `rds_prefix` cannot be the same as `data_prefix`
#' @param outfile         Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param parallel        Logical: should the computations within this function be run in parallel? Defaults to TRUE. See `count_cores()` and `?bigparallelr::assert_cores` for more details.
#'                        In particular, the user should be aware that too much parallelization can make computations *slower*.
#' @param overwrite       Logical: if existing `.bk`/`.rds` files exist for the specified directory/prefix, should these be overwritten? Defaults to FALSE. Set to TRUE if you want to change the imputation method you're using, etc.
#' @param quiet           Logical: should messages be printed to the console? Defaults to TRUE
#'
#' @return `.rds` and `.bk` files are created in `data_dir`, and `obj` (a `bigSNP` object) is returned. See `bigsnpr` documentation for more info on the `bigSNP` class.
#'
#' @keywords internal
#'
read_plink_files <- function(data_dir, data_prefix, rds_dir, rds_prefix, outfile,
                             parallel, overwrite, quiet) {

  # check for compressed files
  if (!file.exists(file.path(data_dir, paste0(data_prefix, ".bed")))) {
    if (file.exists(file.path(data_dir, paste0(data_prefix, ".bed.gz")))) {
      cat("\nIt looks like your files are zipped -- please unzip them before calling process_plink().")
    } else {
      cat("\nThe PLINK files with the specified prefix do not appear in the provided data_dir folder.")
    }
  }

  to_remove <- lapply(c(file.path(rds_dir, data_prefix), file.path(rds_dir, rds_prefix)),
                      \(x) paste0(x, c(".rds", ".desc", ".bk"))) |> unlist()

  # check for overwrite:
  if (any(file.exists(to_remove))) {
    if (overwrite) {
      # notify
      cat("\nOverwriting existing files: ", data_prefix, ".bk/.rds and",
          rds_prefix, ".bk/.rds/.desc\n", sep = "", file = outfile, append = TRUE)

      if (!quiet) {
        cat("\nOverwriting existing files: ", data_prefix, ".bk/.rds and ",
            rds_prefix, ".bk/.rds/.desc\n", sep = "")
      }

      gc() # DO NOT REMOVE - unlink will fail on .bk files otherwise
      unlink(to_remove, force = TRUE)
    } else {
      stop("\nThere are existing .rds and .bk files in the specified directory with the given prefix.
           \nIf you want to overwrite these existing files, set 'overwrite = TRUE'.
           \nOtherwise, choose a different prefix.")
    }
  }

  # create the RDS file  ------------------------
  if (!quiet) cat("\nCreating ", data_prefix, ".rds\n", sep = "")

  if (parallel) {
    bigsnpr::snp_readBed2(bedfile = paste0(file.path(data_dir, data_prefix), ".bed"),
                          backingfile = file.path(rds_dir, data_prefix),
                          ncores = count_cores())
  } else {
    bigsnpr::snp_readBed2(bedfile = paste0(file.path(data_dir, data_prefix), ".bed"),
                          backingfile = file.path(rds_dir, data_prefix))
  }

  bigsnpr::snp_attach(paste0(file.path(rds_dir, data_prefix), ".rds"))
}
