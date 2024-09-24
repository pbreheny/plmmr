#' a function to format the time
#' @keywords internal
#' @returns A string with the formatted current date and time
pretty_time <- function(){
  format(Sys.time(), "%Y-%m-%d %H:%M:%S\n")
}