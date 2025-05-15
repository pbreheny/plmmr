#' a function to return the computer's host name
#'
#' @return String with hostname of current machine
#' @keywords internal
#'
get_hostname <- function() {
  return(as.character(Sys.info()["nodename"]))
}
