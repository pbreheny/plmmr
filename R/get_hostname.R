#' a function to return the computer's host name
#' @return String with hostname of current machine
#' @export
get_hostname <- function() {
  return(as.character(Sys.info()["nodename"]))
}


