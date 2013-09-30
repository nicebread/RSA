\name{bootRSA}
\alias{bootRSA}
\title{Compute bootstrap replications for model parameters}
\usage{
  bootRSA(x, model = "full", ...)
}
\arguments{
  \item{x}{RSA object}

  \item{model}{A string specifying the model; defaults to
  "full"}

  \item{...}{Additional parameters passed to the
  bootstrapLavaan function, e.g., \code{R} for the number
  of bootstrap replications (e.g., R=1000).}
}
\description{
  Compute bootstrap replications for model parameters
}
\details{
  None so far.
}
\examples{
\dontrun{
set.seed(0xBEEF)
n <- 300
err <- 2
x <- rnorm(n, 0, 5)
y <- rnorm(n, 0, 5)
df <- data.frame(x, y)
df <- within(df, {
	diff <- x-y
	absdiff <- abs(x-y)
	sqdiff <- (x-y)^2
	z.sq <- sqdiff + rnorm(n, 0, err)
})

r1 <- RSA(z.sq~x*y, df)
b1 <- bootRSA(r1, model="SSD", R=5000, parallel="multicore", ncpus=2)

# compute percentile confidence intervals and percentile p-value of the bootstrapped values
r1.boot.CI <- CI.boot(b1)
}
}
\seealso{
  \code{\link{RSA}}
}

