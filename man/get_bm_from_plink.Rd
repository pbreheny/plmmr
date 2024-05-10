% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_bm_from_plink.R
\name{get_bm_from_plink}
\alias{get_bm_from_plink}
\title{A function to convert the output of \code{process_plink()} into a \code{big.matrix} (such as could be passed to \code{biglasso::biglasso()})}
\usage{
get_bm_from_plink(path, standardize = TRUE)
}
\arguments{
\item{path}{The file path to the RDS object containing the processed data. Do not add the '.rds' extension to the path.}

\item{standardize}{Logical: should the data be returned in a column-standardized form? Defaults to TRUE. Please don't change this unless you are super confident about what you are doing...}
}
\value{
A list with elements:
\itemize{
\item 'X': a \code{big.matrix} object with the design matrix of the data. If standardize = TRUE, the data corresponding to \code{std_X} is returned.
Otherwise, the data corresponding to \code{subset_X} is returned
\item 'map': data frame with \code{.bim} file information
\item 'fam': data frame with \code{.fam} file information
\item 'ns': A vector indicating the which columns of X contain nonsingular features (i.e., features with variance != 0.
\item 'center': If standardize = TRUE, a vector of values for centering each column in X is returned
\item 'scale': If standardize = TRUE, a vector of values for scaling each column in X is returned
}
}
\description{
A function to convert the output of \code{process_plink()} into a \code{big.matrix} (such as could be passed to \code{biglasso::biglasso()})
}
\details{
This function is a wrapper combining \code{get_data()} and \code{fbm2bm()}
}
\examples{
\donttest{
process_plink(data_dir = get_example_data(parent = TRUE),
  prefix = "penncath_lite",
  gz = TRUE,
  outfile = "process_penncath",
  # overwrite = TRUE, # uncomment if needed 
  impute_method = "mode")
  
  my_path <- paste0(get_example_data(parent = TRUE), "/penncath_lite")
  bm_data <- get_bm_from_plink(my_path)
}
}