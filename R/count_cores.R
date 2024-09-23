#' A helper function to count the number of cores available on the current machine
#'
#' @return A number of cores to use; if `parallel` is installed, this will be `parallel::detectCores()`. Otherwise, this returns a 1.
#' @keywords internal
#'
count_cores <- function(){
  where <- find.package(package = 'parallel', quiet = T)
  if (length(where) == 0) {
    ncores <- 1
  } else {
    ncores <- max(1, parallel::detectCores() - 1)
  }

  return(ncores)
}