#' A function to help with accessing example PLINK files
#'
#' @param path Argument (string) specifying a path (filename) for an external data file in \code{extdata/}
#' @param parent If \code{path=TRUE} and the user wants the name of the parent directory where that file is located, set \code{parent=TRUE}. Defaults to FALSE.
#'
#' @return If \code{path=NULL}, a character vector of file names is returned. If path is given, then a character string
#' with the full file path
#'
#' @export
#'
find_example_data <- function(path, parent=FALSE) {
  if(parent){
    # if parent option selected, will return path to folder
    system.file("extdata", package = "plmmr", mustWork = TRUE)
  } else {
    # if parent option not selected, will return path to folder WITH file name
    system.file("extdata", path, package = "plmmr", mustWork = TRUE)
  }
}

#' a function to unzip the .gz files that ship with the package
#' *Note*: for Linux/Unix and MacOS only
#'
#' @param outdir The filepath to the directory to which the .gz files should be written
unzip_example_data <- function(outdir){

  # Check if the operating system is Windows
  if (.Platform$OS.type == "windows") {
    stop("\nThis function only works for Linux/Unix and MacOS users.
         \nWindows users may use find_example_data() to get the path to the
         PLINK files that ship with the plmmr package, and take it from there.")
  }


  # Set the input and output directories
  input_directory <- find_example_data(parent = TRUE)
  output_directory <- outdir

  # Get a list of all .gz files in the input directory
  gz_files <- list.files(input_directory, pattern = "\\.gz$", full.names = TRUE)

  # Loop through each .gz file
  for (gz_file in gz_files) {
    # Get the base name of the .gz file
    base_name <- basename(gz_file)

    # Remove the .gz extension
    output_file <- file.path(output_directory, sub("\\.gz$", "", base_name))

    # Unzip the .gz file and save it to the output directory
    system(paste("gzip -dc", shQuote(gz_file), ">", shQuote(output_file)))
  }

  # Print a success message
  cat("Unzipped files are saved in", output_directory, "\n")

}
