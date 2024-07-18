#' A function to read in PLINK files using `bigsnpr` methods
#'
#' @param data_dir  The path to the bed/bim/fam data files, *without* a trailing "/" (e.g., use `data_dir = '~/my_dir'`, **not** `data_dir = '~/my_dir/'`)
#' @param data_prefix    The prefix (as a character string) of the bed/fam data files (e.g., `prefix = 'mydata'`)
#' @param rds_dir   The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`
#' @param outfile   Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param overwrite Logical: if existing `.bk`/`.rds` files exist for the specified directory/prefix, should these be overwritten? Defaults to FALSE. Set to TRUE if you want to change the imputation method you're using, etc.
#' @param quiet Logical: should messages be printed to the console? Defaults to TRUE
#' @returns '.rds' and '.bk' files are created in `data_dir`, and `obj` (a `bigSNP` object) is returned. See `bigsnpr` documentation for more info on the `bigSNP` class.
#' @keywords internal
#'
read_plink_files <- function(data_dir, data_prefix, rds_dir, outfile, overwrite, quiet){

  # check for compressed files
  if (!file.exists(file.path(data_dir, paste0(data_prefix, ".bed")))){
    if (file.exists(file.path(data_dir, paste0(data_prefix, ".bed.gz")))) {
      cat("\nIt looks like your files are zipped -- please unzip them before calling process_plink().")
    } else {
      cat("\nThe PLINK files with the specified prefix do not appear in the provided data_dir folder.")
    }
  }


  path <- file.path(rds_dir, paste0(data_prefix, ".rds"))
  bk_path <- file.path(rds_dir, paste0(data_prefix, ".bk"))

  # check for overwrite:
  if (file.exists(bk_path)){
    if (overwrite){
      # notify
      cat("\nOverwriting existing files: ", data_prefix, ".bk/.rds\n",
          file = outfile, append = TRUE)

      if (!quiet){
        cat("\nOverwriting existing files: ", data_prefix, ".bk/.rds\n")
      }

      # overwrite existing rds file
      gc()
      if (file.exists(path)) file.remove(path)
      gc()
      # remove old backingfile(s)
      list.files(rds_dir, pattern=paste0('^.*.bk'), full.names=TRUE) |>
        file.remove()
      gc()
    } else {
      stop("\nThere are existing data_prefix.rds and data_prefix.bk files in the specified directory.
           \nIf you want to overwrite these existing files, set 'overwrite = TRUE'.
           \nOtherwise, choose a different prefix.")
    }
  }

  # create the RDS file  ------------------------
  if (!quiet) cat("\nCreating ", data_prefix, ".rds\n", sep='')

  bigsnpr::snp_readBed2(bedfile = paste0(file.path(data_dir, data_prefix), ".bed"),
                        backingfile = file.path(rds_dir, data_prefix),
                        ncores = bigstatsr::nb_cores())

  obj <- bigsnpr::snp_attach(paste0(file.path(rds_dir, data_prefix), ".rds"))
  return(obj)
}
