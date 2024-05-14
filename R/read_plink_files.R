#' A function to read in PLINK files using `bigsnpr` methods
#'
#' @param data_dir  The path to the bed/bim/fam data files, *without* a trailing "/" (e.g., use `data_dir = '~/my_dir'`, **not** `data_dir = '~/my_dir/'`)
#' @param prefix    The prefix (as a character string) of the bed/fam data files (e.g., `prefix = 'mydata'`)
#' @param rds_dir   The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`
#' @param gz        Logical: are the bed/bim/fam files g-zipped? Defaults to FALSE. NOTE: if TRUE, process_plink will unzip your zipped files.
#' @param outfile   Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param overwrite Logical: if existing `.bk`/`.rds` files exist for the specified directory/prefix, should these be overwritten? Defaults to FALSE. Set to TRUE if you want to change the imputation method you're using, etc.
#' @param quiet Logical: should messages be printed to the console? Defaults to TRUE
#' @returns '.rds' and '.bk' files are created in `data_dir`, and `obj` (a `bigSNP` object) is returned. See `bigsnpr` documentation for more info on the `bigSNP` class.
#' @keywords internal
#'
read_plink_files <- function(data_dir, prefix, rds_dir, gz, outfile, overwrite, quiet){
  path <- file.path(rds_dir, paste0(prefix, ".rds"))
  bk_path <- file.path(rds_dir, paste0(prefix, ".bk"))
  std_bk_path <- file.path(rds_dir, paste0("std_", prefix, ".bk"))
  sub_bk_path <- file.path(rds_dir, paste0("subset_", prefix, ".bk"))

  # check for overwrite:
  if (file.exists(bk_path)){
    if (overwrite){
      # notify
      cat("\nOverwriting existing files: ", prefix, ".bk/.rds\n",
          file = outfile, append = TRUE)

      if (!quiet){
        cat("\nOverwriting existing files: ", prefix, ".bk/.rds\n")
      }

      # overwrite existing files
      system(paste0("rm ", bk_path))
      system(paste0("rm ", std_bk_path))
      system(paste0("rm ", sub_bk_path))
      system(paste0("rm ", path))
    } else {
      stop("\nThere are existing prefix.rds and prefix.bk files in the specified directory.
           \nIf you want to overwrite these existing files, set 'overwrite = TRUE'.
           \nOtherwise, choose a different prefix.")
    }
  }

  # create the RDS file first ------------------------
  cat("\nCreating ", prefix, ".rds\n", file = outfile, append = TRUE)
  if(!quiet){
    cat("\nCreating ", prefix, ".rds\n")
  }
  # check for compressed files
  if (file.exists(file.path(data_dir, paste0(prefix, ".bed.gz")))) {
    gz <- TRUE
  }
  if (gz){
    cat("\nUnzipping .gz files - this could take a second", file = outfile, append = TRUE)
    if (!quiet) cat("\nUnzipping .gz files - this could take a second")
    system(paste0("gunzip -k ", file.path(data_dir, paste0(prefix, "*"))))
  }

  bigsnpr::snp_readBed2(bedfile = paste0(data_dir, "/", prefix, ".bed"),
                        backingfile = paste0(rds_dir, "/", prefix),
                        ncores = bigstatsr::nb_cores())

  obj <- bigsnpr::snp_attach(paste0(rds_dir, "/", prefix, ".rds"))

  return(obj)
}
