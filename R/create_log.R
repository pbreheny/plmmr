#' Create the `.log` file
#'
#' @param outfile String specifying the name of the to-be-created file, *without* extension
#'
#' @return Nothing is returned, instead a text file with the suffix `.log` is created.
#' If outfile is NULL, the path to the null device is returned.
#'
#' @keywords internal
#'
create_log <- function(outfile) {
  if (missing(outfile)) {
    stop(
      "You must specify a name for the output file(s) via the outfile argument."
    )
  }

  if (is.null(outfile)) {
    logfile <- nullfile()
  } else {
    logfile <- paste0(outfile, ".log")
  }

  # open the log file for writing
  log_con <- file(logfile)

  hostname <- as.character(Sys.info()["nodename"])

  # write header to the log file
  cat("### plmmr log file ###\n", file = logfile)
  cat("Logging to", logfile, "\n", file = logfile, append = TRUE)
  cat("Host:", hostname, "\n", file = logfile, append = TRUE)
  cat(
    "Current working directory:",
    getwd(),
    "\n",
    file = logfile,
    append = TRUE
  )
  cat(
    "Start log at:",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S\n"),
    file = logfile,
    append = TRUE
  )

  # get the name of the calling function
  syscall <- sys.calls()
  calling_function <- deparse(syscall[[1]])
  cat("Call:", trimws(calling_function), "\n", file = logfile, append = TRUE)

  # close the log file
  close(log_con)

  # return logfile name
  logfile
}
