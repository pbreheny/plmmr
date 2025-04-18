% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-cv-plmm.R
\name{plot.cv_plmm}
\alias{plot.cv_plmm}
\title{Plot method for cv_plmm class}
\usage{
\method{plot}{cv_plmm}(
  x,
  log.l = TRUE,
  type = c("cve", "rsq", "scale", "snr", "all"),
  selected = TRUE,
  vertical.line = TRUE,
  col = "red",
  ...
)
}
\arguments{
\item{x}{An object of class cv_plmm}

\item{log.l}{Logical to indicate the plot should be returned on the natural log scale. Defaults to \code{log.l = FALSE}.}

\item{type}{Type of plot to return. Defaults to "cve."}

\item{selected}{Logical to indicate which variables should be plotted. Defaults to TRUE.}

\item{vertical.line}{Logical to indicate whether vertical line should be plotted at the minimum/maximum value. Defaults to TRUE.}

\item{col}{Color for vertical line, if plotted. Defaults to "red."}

\item{...}{Additional arguments.}
}
\value{
Nothing is returned; instead, a plot is drawn representing the relationship
between the tuning parameter 'lambda' value (x-axis) and the cross validation error (y-axis).
}
\description{
Plot method for cv_plmm class
}
\examples{
admix_design <- create_design(X = admix$X, y = admix$y)
cvfit <- cv_plmm(design = admix_design)
plot(cvfit)

}
