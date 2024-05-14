#' A function to read in a large file as a numeric file-backed matrix (`FBM`)
#' Note: this function is a wrapper for `bigstatsr::big_read()`
#' @param file      The name of the file to read, not including its directory. Directory should be specified in `data_dir`
#' @param data_dir  The path to the directory where 'file' is 
#' @param rds_dir   The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`
#' @param ind.col   The indices (positive integers) that are sorted in ascending order. Passed to `bigstatsr::big_read()`
#' @param outfile   Optional: the name (character string) of the prefix of the logfile to be written. Defaults to 'process_plink', i.e. you will get 'process_plink.log' as the outfile.
#' @param overwrite Logical: if existing `.bk`/`.rds` files exist for the specified directory/prefix, should these be overwritten? Defaults to FALSE. Set to TRUE if you want to change the imputation method you're using, etc. 
#' @param quiet     Logical: should messages be printed to the console? Defaults to TRUE
#' 
#' @returns '.rds' and '.bk' files are created in `data_dir`, and `obj` (a `bigSNP` object) is returned. See `bigsnpr` documentation for more info on the `bigSNP` class.
#' @keywords internal
#'
read_data_files <- function(file, data_dir, rds_dir, ind.col, outfile, overwrite, quiet){
  
  prefix <- unlist(strsplit(file, split = "\\."))[1]
  path <- file.path(rds_dir,  paste0(prefix, ".rds"))
  bk_path <- file.path(rds_dir, paste0(prefix, ".bk"))
  std_bk_path <- file.path(rds_dir, paste0("std_", prefix, ".bk"))

  # check for overwrite: 
  if (file.exists(bk_path)){
    if (overwrite){
      # notify 
      cat("\nOverwriting existing files: ", prefix, ".bk/.rds\n",
          file = outfile, append = TRUE, sep='')
      
      if (!quiet){
        cat("\nOverwriting existing files: ", prefix, ".bk/.rds\n", sep='')
      }
      
      # overwrite existing files
      file.remove(bk_path)
      file.remove(std_bk_path)
      file.remove(path)
    } else {
      stop("\nThere are existing .rds and/or .bk files in the specified directory.  
           \nIf you want to overwrite these existing files, set 'overwrite = TRUE'. 
           \nOtherwise, choose a different prefix.")
    }
  }

  # create the RDS file ------------------------
  bigstatsr::big_read(file = file.path(data_dir, file),
                      select = ind.col,
                      backingfile = file.path(rds_dir, prefix))
  
  obj <- bigstatsr::big_attach(file.path(rds_dir, paste0(prefix, ".rds")))
  
  
  return(obj)
  
}
