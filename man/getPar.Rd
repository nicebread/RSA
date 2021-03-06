% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPar.R
\name{getPar}
\alias{getPar}
\title{Retrieves several variables from an RSA object}
\usage{
getPar(x, type = "coef", model = "full", digits = NA, ...)
}
\arguments{
\item{x}{RSA object}

\item{type}{One of: "syntax", "coef", "R2", "R2.adj", "free", "summary", "p.value"}

\item{model}{A string specifying the model; defaults to "full"}

\item{digits}{Number of digits the output is rounded to; if NA, digits are unconstrained}

\item{...}{Additional parameters passed to the extraction function}
}
\description{
Retrieves several variables from an RSA object
}
\details{
None so far.
}
\examples{

set.seed(0xBEEF)
n <- 300
err <- 2
x <- rnorm(n, 0, 5)
y <- rnorm(n, 0, 5)
df <- data.frame(x, y)
df <- within(df, {
	diff <- x-y
	absdiff <- abs(x-y)
	SD <- (x-y)^2
	z.sq <- SD + rnorm(n, 0, err)
})

r1 <- RSA(z.sq~x*y, df, models=c("full", "SSQD"))
getPar(r1, "syntax")
getPar(r1, "R2")
getPar(r1, "coef")
}
\seealso{
\code{\link{RSA}}
}
