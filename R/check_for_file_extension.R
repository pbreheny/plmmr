#' check_for_file_extension: a function to make our package 'smart' enough to
#' handle .rds file extensions
#' @param path A string specifying a file path that ends in a file name, e.g. "~/dir/my_file.rds"
#'
#' @return a string with a filepath *without* an extension, e.g. "~/dir/my_file"
#' @keywords internal

check_for_file_extension <- function(path){
  if (grepl('.rds', path)) {
    path <- unlist(strsplit(path, split = ".rds", fixed = TRUE))

  }

  if (grepl('.bk', string_path)) {
    path <- unlist(strsplit(path, split = ".bk", fixed = TRUE))
  }

  return(path)

}