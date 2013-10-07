\name{confint.RSA}
\alias{confint}
\alias{confint.RSA}
\title{Computes confidence intervals for RSA parameters, standard or bootstrapped}
\usage{
  \method{confint}{RSA} (object, parm, level = 0.95, ...,
    model = "full", method = "standard", R = 5000)
}
\arguments{
  \item{object}{An RSA object}

  \item{parm}{Not used.}

  \item{level}{The confidence level required.}

  \item{model}{A string specifying the model; defaults to
  "full"}

  \item{method}{"standard" returns the CI for the lavaan
  object as it was computed. "boot" computes new percentile
  bootstrapped CIs.}

  \item{R}{If \code{method = "boot"}, R specifies the
  number of bootstrap samples}

  \item{...}{Additional parameters passed to the
  bootstrapLavaan function, e.g.,
  \code{parallel="multicore", ncpus=2}.}
}
\description{
  Computes confidence intervals for RSA parameters,
  standard or bootstrapped
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
	sqdiff <- (x-y)^2
	z.sq <- sqdiff + rnorm(n, 0, err)
})

r1 <- RSA(z.sq~x*y, df, models="SSD")
(c1 <- confint(r1, model="SSD"))

# Dummy example with 10 bootstrap replications - better use >= 5000!
(c2 <- confint(r1, model="SSD", method="boot", R=10))
\dontrun{
# multicore version
confint(r1, model="SSD", R=5000, parallel="multicore", ncpus=2)
}
}
\seealso{
  \code{\link{RSA}}
}

