#' helper function to get the file path of a file *without* the extension
#'
#' @param path The path of a file
#'
#' @return path_sans_ext The filepath without the extension
#'
#' @keywords internal
file_sans_ext <- function(path){
  f <- strsplit(x = path, ".", fixed = T) |> unlist()
  return(f[1])
}