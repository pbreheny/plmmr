#' A function to format the time
#'
#' @return A string with the formatted current date and time
#'
#' @keywords internal
#'
pretty_time <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S\n")
}
