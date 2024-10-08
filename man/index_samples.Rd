% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index_samples.R
\name{index_samples}
\alias{index_samples}
\title{A function to align genotype and phenotype data}
\usage{
index_samples(
  obj,
  rds_dir,
  indiv_id,
  add_outcome,
  outcome_id,
  outcome_col,
  na_outcome_vals,
  outfile,
  quiet
)
}
\arguments{
\item{obj}{An object created by \code{process_plink()}}

\item{rds_dir}{The path to the directory in which you want to create the new '.rds' and '.bk' files.}

\item{indiv_id}{A character string indicating the ID column name in the 'fam'
element of the genotype data list. Defaults to 'sample.ID', equivalent to 'IID' in PLINK. The other option is 'family.ID', equivalent to 'FID' in PLINK.}

\item{add_outcome}{A data frame with at least two columns: and ID column and a phenotype column}

\item{outcome_id}{A string specifying the name of the ID column in \code{pheno}}

\item{outcome_col}{A string specifying the name of the phenotype column in \code{pheno}. This column will be used as the default \code{y} argument to 'plmm()'.}

\item{na_outcome_vals}{A vector of numeric values used to code NA values in the outcome. Defaults to \code{c(-9, NA_integer)} (the -9 matches PLINK conventions).}

\item{outfile}{A string with the name of the filepath for the log file}

\item{quiet}{Logical: should messages be printed to the console? Defaults to FALSE (which leaves the print messages on...}
}
\value{
a list with two items:
\itemize{
\item a data.table with rows corresponding to the samples for which both genotype and phenotype are available.
\item a numeric vector with indices indicating which samples were 'complete' (i.e., which samples from add_outcome had corresponding data in the PLINK files)
}
}
\description{
A function to align genotype and phenotype data
}
\keyword{internal}
