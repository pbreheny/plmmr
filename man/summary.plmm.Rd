% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary-plmm.R
\name{summary.plmm}
\alias{summary.plmm}
\title{A summary method for the plmm objects}
\usage{
\method{summary}{plmm}(object, lambda, idx, eps = 1e-05, ...)
}
\arguments{
\item{object}{An object of class \code{plmm}}

\item{lambda}{The regularization parameter value at which inference should be reported.}

\item{idx}{Alternatively, \code{lambda} may be specified by an index; \code{idx=10} means:
report inference for the 10th value of \code{lambda} along the regularization path. If both \code{lambda} and \code{idx} are specified, \code{lambda} takes precedence.}

\item{eps}{If lambda is given, eps is the tolerance for difference between the given lambda value and a lambda value from the object. Defaults to 0.0001 (1e-5)}

\item{...}{Not used}
}
\value{
The return value is an object with S3 class \code{summary.plmm}. The class has its own print method and contains the following list elements:
\itemize{
\item \code{penalty}: The penalty used by \code{plmm} (e.g. SCAD, MCP, lasso)
\item \code{n}: Number of instances/observations
\item \code{std_X_n}: the number of observations in the standardized data;
the only time this would differ from 'n' is if data are from PLINK and the external data does not include all the same samples
\item \code{p}: Number of regression coefficients (not including the intercept)
\item \code{converged}: Logical indicator for whether the model converged
\item \code{lambda}: The \code{lambda} value at which inference is being reported
\item \code{lambda_char}: A formatted character string indicating the lambda value
\item \code{nvars}: The number of nonzero coefficients (again, not including the intercept) at that value of \code{lambda}
\item \code{nonzero}: The column names indicating the nonzero coefficients in the model at the specified value of \code{lambda}
}
}
\description{
A summary method for the plmm objects
}
\examples{
admix_design <- create_design(X = admix$X, y = admix$y)
fit <- plmm(design = admix_design)
summary(fit, idx = 97)
}
