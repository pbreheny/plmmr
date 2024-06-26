% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print-summary-cv-plmm.R
\name{print.summary.cv_plmm}
\alias{print.summary.cv_plmm}
\title{Print method for summary.cv_plmm objects}
\usage{
\method{print}{summary.cv_plmm}(x, digits, ...)
}
\arguments{
\item{x}{An object of class \code{summary.cv_plmm}}

\item{digits}{The number of digits to use in formatting output}

\item{...}{Not used}
}
\value{
Nothing is returned; instead, a message is printed to the console
summarizing the results of the cross-validated model fit.
}
\description{
Print method for summary.cv_plmm objects
}
\examples{
cv_fit <- cv_plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
print(summary(cv_fit))

}
