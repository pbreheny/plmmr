#' A function to help with accessing example PLINK files.
#'
#' @param path Argument (string) specifying a path (filename) for an external data file in \code{extdata/}.
#' @param parent If the user wants the name of the parent directory where the example data is located, set \code{parent=TRUE}. Defaults to FALSE.
#'
#' @returns If \code{path=NULL}, a character vector of file names is returned. If path is given, then a character string
#' with the full file path.
#'
#' @export
#'
#' @examples
#' find_example_data(parent = TRUE)
#'
find_example_data <- function(path, parent = FALSE) {
  if (parent) {
    # if parent option selected, will return path to folder
    system.file("extdata", package = "plmmr", mustWork = TRUE)
  } else {
    # if parent option not selected, will return path to folder WITH file name
    system.file("extdata", path, package = "plmmr", mustWork = TRUE)
  }
}

#' Companion function to unzip the .gz files that ship with the `plmmr` package.
#'
#' @param outdir The file path to the directory to which the .gz files should be written.
#'
#' @details For an example of this function, look at `vignette('plink_files', package = "plmmr")`.
#'
#' @returns Nothing is returned; the PLINK files that ship with the `plmmr` package are stored in the directory specified by 'outdir'.
#'
#' @export
#'
unzip_example_data <- function(outdir) {
  # Set the input directory
  indir <- find_example_data(parent = TRUE)

  # Get a list of all .gz files in the input directory
  gz_files <- list.files(indir, pattern = "\\.gz$", full.names = TRUE)

  # Loop through and unzip each .gz file
  for(i in seq_along(gz_files)) {
    file_base <- tools::file_path_sans_ext(basename(gz_files[i]))
    R.utils::gunzip(gz_files[i],
                    destname = file.path(outdir, file_base),
                    overwrite = TRUE, remove = FALSE)
  }

  # Print a success message
  message("Unzipped files are saved in ", outdir)
}
