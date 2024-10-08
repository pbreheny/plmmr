% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize_matrix.R
\name{standardize_in_memory}
\alias{standardize_in_memory}
\title{A helper function to standardize matrices}
\usage{
standardize_in_memory(X)
}
\arguments{
\item{X}{a matrix}
}
\value{
a list with the standardized matrix, vectors with the centering/scaling values, and a vector with the indices of nonsingular columns
}
\description{
A helper function to standardize matrices
}
\details{
This function is adapted from https://github.com/pbreheny/ncvreg/blob/master/R/std.R
NOTE: this function returns a matrix \strong{in memory}. For standardizing filebacked
data, use \code{big_std()}  -- see src/big_standardize.cpp
}
\keyword{internal}
