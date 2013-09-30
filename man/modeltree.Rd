\name{modeltree}
\alias{modeltree}
\title{Plots a flow chart with model comparisons}
\usage{
  modeltree(x, digits = 3, ...)
}
\arguments{
  \item{x}{A cRSA object (= output from the
  \code{\link{compare}} function)}

  \item{digits}{The number of digits to which numbers are
  rounded}

  \item{...}{Additional parameters (not used yet)}
}
\description{
  Plots a flow chart with model comparisons from a RSA
  object
}
\details{
  The plot can be either requested within the
  \code{compare} function: \code{compare(r1, plot=TRUE)} Or
  it can be plotted from a cRSA object (= output from the
  \code{\link{compare}} function): \code{c1 <- compare(r1)}
  \code{plot(c1)}
}
\seealso{
  \code{\link{RSA}}, \code{\link{compare}}
}

