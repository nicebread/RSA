\name{RSA}
\alias{RSA}
\title{Performs several RSA model tests on a data set with two predictors}
\usage{
  RSA(formula, data = NULL, center = FALSE, scale = FALSE,
    na.rm = FALSE, out.rm = TRUE, breakline = FALSE,
    models = "default", cubic = FALSE, verbose = TRUE,
    control.variables = c(), ...)
}
\arguments{
  \item{formula}{A formula in the form \code{z ~ x*y},
  specifying the variable names used from the data frame,
  where z is the name of the response variable, and x and y
  are the names of the predictor variables.}

  \item{data}{A data frame with the variables}

  \item{center}{Should predictor variables be centered on
  the sample mean before analyses?}

  \item{scale}{Should predictor variables be scaled to SD =
  1 before analyses?}

  \item{na.rm}{Remove missings before proceeding?}

  \item{out.rm}{Should outliers according to Bollen &
  Jackman (1980) criteria be excluded from analyses?}

  \item{breakline}{Should the breakline in the
  unconstrained absolute difference model be allowed (the
  breakline is possible from the model formulation, but
  empirically rather unrealistic ...)}

  \item{verbose}{Should additional information during the
  computation process be printed?}

  \item{models}{A vector with names of all models that
  should be computed. Should be any from c("absdiff",
  "absunc", "diff", "additive", "IA", "sqdiff", "RR",
  "SSD", "SRSD", "full", "null"). For \code{models="all"},
  all models are computed, for \code{models="default"} all
  models besides absolute difference models are computed.}

  \item{cubic}{Should a cubic model with the additional
  terms Y^3, XY^2, YX^2, and X^3 be included?}

  \item{control.variables}{A string vector with variable
  names from \code{data}. These variables are added as
  linear predictors to the model (in order "to control for
  them"). No interactions with the other variables are
  modeled.}

  \item{...}{Additional parameters passed to the lavaan sem
  function. For example: \code{se="boot"}}
}
\description{
  Performs several RSA model tests on a data set with two
  predictors
}
\details{
  You can also fit binary outcome variables with a probit
  link function. For that purpose, the response variable
  has to be defined as "ordered": \code{r1 <- RSA(Z.binary
  ~ X*Y, dat, ordered="Z.binary")} (for more details see
  the help file of the \code{sem} function in the
  \code{lavaan} package.). The results can also be plotted
  with probabilities on the z axis using the probit link
  function: \code{plot(r1, link="probit", zlim=c(0, 1),
  zlab="Probability")}. \code{lavaan} at the moment only
  supports a probit link function for binary outcomes, not
  a logit link.
}
\examples{
# Compute response surface from a fake data set
set.seed(0xBEEF)
n <- 300
err <- 15
x <- rnorm(n, 0, 5)
y <- rnorm(n, 0, 5)
df <- data.frame(x, y)
df <- within(df, {
	diff <- x-y
	absdiff <- abs(x-y)
	sqdiff <- (x-y)^2
	z.diff <- diff + rnorm(n, 0, err)
	z.abs <- absdiff + rnorm(n, 0, err)
	z.sq <- sqdiff + rnorm(n, 0, err)
	z.add <- diff + 0.4*x + rnorm(n, 0, err)
	z.complex <- 0.4*x + - 0.2*x*y + + 0.1*x^2 - 0.03*y^2 + rnorm(n, 0, err)
})
\dontrun{
r1 <- RSA(z.sq~x*y, df)
print(r1)
compare(r1)
plot(r1)
plot(r1, model="SRSD")
plot(r1, model="full", type="c")
getPar(r1, "coef")	# print model parameters including SE and CI
RSA.ST(r1)	# get surface parameters

# Motive congruency example
data(motcon)
r.m <- RSA(postVA~ePow*iPow, motcon)

# Get boostrapped CIs with 10 bootstrap samples (usually this should be set to 5000 or higher),
# only from the SSD model
c1 <- confint(r.m, model="SSD", method="boot", R=10)

# Plot the final model
plot(r.m, model="RR", xlab="Explicit power motive",
		ylab="Implicit power motive", zlab="Affective valence")
}
}
\seealso{
  \code{\link{demoRSA}}, \code{\link{plotRSA}},
  \code{\link{RSA.ST}}, \code{\link{confint.RSA}}
}

