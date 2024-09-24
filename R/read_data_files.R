#' A function to read in a large file as a numeric file-backed matrix (`FBM`)
#' Note: this function is a wrapper for `bigstatsr::big_read()`
#' @param data_file The name of the file to read, not including its directory. Directory should be specified in `data_dir`
#' @param data_dir  The path to the directory where 'file' is
#' @param rds_dir   The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`
#' @param rds_prefix    String specifying the user's preferred filename for the to-be-created .rds/.bk files (will be create insie `rds_dir` folder)
#'                      Note: 'rds_prefix' cannot be the same as 'data_file'
#' @param outfile   Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param overwrite Logical: if existing `.bk`/`.rds` files exist for the specified directory/prefix, should these be overwritten? Defaults to FALSE. Set to TRUE if you want to change the imputation method you're using, etc.
#' @param quiet     Logical: should messages be printed to the console? Defaults to TRUE
#' @param ...           Optional: other arguments to be passed to `bigmemory::read.big.matrix()`. Note: 'sep' is an option to pass here.
#'
#' @returns '.rds', '.bk', and '.desc' files are created in `data_dir`, and `obj` (a filebacked `bigmemory big.matrix` object) is returned. See `bigmemory` documentation for more info on the `big.matrix` class.
#'
#' @keywords internal
#'
read_data_files <- function(data_file, data_dir, rds_dir, rds_prefix,
                            outfile, overwrite, quiet, ...){

  rds_path <- file.path(rds_dir,  paste0(rds_prefix, ".rds"))
  bk_path <- file.path(rds_dir, paste0(rds_prefix, ".bk"))
  desc_path <- file.path(rds_dir, paste0(rds_prefix, ".desc"))

  # check for overwrite:
  if (file.exists(bk_path)){
    if (overwrite){
      # notify
      cat("\nOverwriting existing files:", rds_prefix, ".bk/.rds/.desc\n",
          file = outfile, append = TRUE, sep='')

      if (!quiet){
        cat("\nOverwriting existing files:", rds_prefix, ".bk/.rds/.desc\n", sep='')
      }
      file.remove(bk_path)
      file.remove(rds_path)
      file.remove(desc_path)
    } else {
      stop("\nThere are existing .rds and/or .bk files in the specified directory.
           \nIf you want to overwrite these existing files, set 'overwrite = TRUE'.
           \nOtherwise, choose a different prefix.")
    }
  }

  # create the bk file ------------------------
  obj <- bigmemory::read.big.matrix(filename = file.path(data_dir, data_file),
                             backingfile = paste0(rds_prefix, ".bk"),
                             backingpath = rds_dir,
                             descriptorfile = paste0(rds_prefix, ".desc"),
                             type = 'double',
                             ...)

  return(obj)

}
