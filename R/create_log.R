# Define the function to create the .log file
#' create_log_file
#' @param outfile   String specifying the name of the to-be-created file, *without* extension
#' @param ...       Not used
#' @keywords internal
#'
#' @returns Nothing is returned, intead a text file with the suffix .log is created.
create_log <- function(outfile, ...) {
  if (missing(outfile)) {
    stop("You must specify a name for the output file(s) via the outfile argument.")
  }

  logfile <- paste0(outfile, ".log")

  # open the log file for writing
  log_con <- file(logfile)

  # write header to the log file
  cat("### plmmr log file ###\n", file = logfile)
  cat("Logging to", logfile, "\n", file = logfile, append = TRUE)
  cat("Host:", get_hostname(), "\n", file = logfile, append = TRUE)
  cat("Current working directory:", getwd(), "\n", file = logfile, append = TRUE)
  cat("Start log at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S\n"),
      file = logfile, append = TRUE)

  # get the name of the calling function
  calling_function <- deparse(sys.call(-1))
  cat("Call:", trimws(calling_function), "\n", file = logfile, append = TRUE)

  # close the log file
  close(log_con)

  # return logfile name
  return(logfile)
}
