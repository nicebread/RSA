\name{plotRSA}
\alias{plotRSA}
\title{Plots a response surface of a polynomial equation of second degree}
\usage{
  plotRSA(x = NULL, y = 0, x2 = 0, y2 = 0, xy = 0, w = 0,
    wx = 0, wy = 0, x3 = 0, xy2 = 0, x2y = 0, y3 = 0,
    b0 = 0, type = "3d", model = "full", xlim = NULL,
    ylim = NULL, zlim = NULL, xlab = NULL, ylab = NULL,
    zlab = NULL, surface = "predict", lambda = NULL,
    rotation = list(x = -63, y = 32, z = 15),
    label.rotation = list(x = 19, y = -40, z = 92),
    gridsize = 21, bw = FALSE, legend = TRUE, param = TRUE,
    axes = c("LOC", "LOIC", "PA1", "PA2"), project = FALSE,
    maxlines = FALSE, cex = 1.2,
    points = list(show = TRUE, value = "raw", jitter = 0, color = "black", cex = 0.5),
    demo = FALSE, fit = NULL, link = "identity",
    tck = c(1.5, 1.5, 1.5), distance = c(1.3, 1.3, 1.4),
    border = TRUE,
    contour = list(show = FALSE, color = "grey40", highlight = c()),
    hull = FALSE, SP.CI = FALSE, pal = NULL,
    pal.range = "box", pad = 0, ...)
}
\arguments{
  \item{x}{Either an RSA object (returned by the \code{RSA}
  function), or the coefficient for the X predictor}

  \item{y}{Y coefficient}

  \item{x2}{X^2 coefficient}

  \item{y2}{Y^2 coefficient}

  \item{xy}{XY interaction coefficient}

  \item{w}{W coefficient (for (un)constrained absolute
  difference model)}

  \item{wx}{WX coefficient (for (un)constrained absolute
  difference model)}

  \item{wy}{WY coefficient (for (un)constrained absolute
  difference model)}

  \item{y3}{Y^3 coefficient}

  \item{x3}{X^3 coefficient}

  \item{xy2}{XY^2 coefficient}

  \item{x2y}{X^2Y coefficient}

  \item{b0}{Intercept}

  \item{xlim}{Limits of the x axis}

  \item{ylim}{Limits of the y axis}

  \item{zlim}{Limits of the z axis}

  \item{xlab}{Label for x axis}

  \item{ylab}{Label for y axis}

  \item{zlab}{Label for z axis}

  \item{surface}{Method for the calculation of the surface
  z values. "predict" takes the predicted values from the
  model, "smooth" uses a thin plate smoother (function
  \code{Tps} from the \code{fields} package) of the raw
  data}

  \item{lambda}{lambda parameter for the smoother. Default
  (NULL) means that it is estimated by the smoother
  function. Small lambdas around 1 lead to rugged surfaces,
  big lambdas to very smooth surfaces.}

  \item{rotation}{Rotation of the 3d surface plot (when
  type == "3d")}

  \item{label.rotation}{Rotation of the axis labls (when
  type == "3d")}

  \item{gridsize}{Number of grid nodes in each dimension}

  \item{bw}{Print surface in black and white instead of
  colors?}

  \item{legend}{Print color legend for z values?}

  \item{cex}{Font size factor for axes labels}

  \item{type}{\code{3d} for 3d surface plot, \code{contour}
  for 2d contour plot, "interactive" for interactive
  rotatable plot. Shortcuts (i.e., first letter of string)
  are sufficient}

  \item{points}{A list of parameters which define the
  appearance of the raw scatter points: show = TRUE: Should
  the original data points be overplotted? value="raw":
  Plot the original z value, "predicted": plot the
  predicted z value. jitter=0: Amount of jitter for the raw
  data points. cex = .5: multiplication factor for point
  size.}

  \item{model}{If x is an RSA object: from which model
  should the response surface be computed?}

  \item{demo}{Do not change that parameter (internal use
  only)}

  \item{fit}{Do not change that parameter (internal use
  only)}

  \item{param}{Should the surface parameters a1 to a4 be
  shown on the plot? In case of a 3d plot a1 to a4 are
  printed on the upper left side; in case of a contour plot
  the principal axes are plotted.}

  \item{axes}{A vector of strings specifying the axes that
  should be plotted. Can be any combination of c("LOC",
  "LOIC", "PA1", "PA2"). LOC = line of congruence, LOIC =
  line of incongruence, PA1 = first principal axis, PA2 =
  second principal axis}

  \item{project}{Should the LOC, LOIC, etc. (as defined in
  parameter \code{axes}) be also plotted as a projection on
  the bottom of the cube?}

  \item{maxlines}{Should the maximum lines be plotted?
  (red: maximum X for a given Y, blue: maximum Y for a
  given X). Works only in type="3d"}

  \item{link}{Link function to transform the z axes.
  Implemented are "identity" (no transformation; default),
  "probit", and "logit"}

  \item{border}{Should a thicker border around the surface
  be plotted? Sometimes this border leaves the surrounding
  box, which does not look good. In this case the border
  can be suppressed by setting \code{border=FALSE}.}

  \item{contour}{A list defining the appearance of contour
  lines (aka. height lines). show=TRUE: Should the contour
  lines be plotted on the 3d wireframe plot? (Parameter
  only relevant for \code{type="3d"}). color = "grey40":
  Color of the contour lines. highlight = c(): A vector of
  heights which should be highlighted (i.e., printed in
  bold). Be careful: the highlighted line is not
  necessarily exactly at the specified height; instead the
  nearest height line is selected.}

  \item{hull}{Plot a bag plot on the surface (This is a
  bivariate extension of the boxplot. 50% of points are in
  the inner bag, 50% in the outer region). See Rousseeuw,
  Ruts, & Tukey (1999).}

  \item{SP.CI}{Plot the CI of the stationary point (only
  relevant for \code{type="contour"})}

  \item{distance}{A vector of three values defining the
  distance of labels to the axes}

  \item{tck}{A vector of three values defining the position
  of labels to the axes (see ?wireframe)}

  \item{pal}{A palette for shading}

  \item{pal.range}{Should the color range be scaled to the
  box (\code{pal.range = "box"}, default), or to the min
  and max of the surface (\code{pal.range = "surface"})? If
  set to "box", different surface plots can be compared
  along their color, as long as the zlim is the same for
  both.}

  \item{pad}{Pad controls the margin around the figure
  (positive numbers: larger margin, negative numbers:
  smaller margin)}

  \item{...}{Additional parameters passed to the plotting
  function (e.g., main="Title"). A useful title might be
  the R squared of the plotted model: main =
  as.expression(bquote(R^2==.(round(getPar(x, "r2",
  model="full"), 3))))}
}
\description{
  Plots an RSA object, or a response surface with specified
  parameters
}
\details{
  Each plot type has its distinctive advantages. The
  two-dimensional contour plot gives a clear view of the
  position of the principal axes and the stationary point.
  The 3d plot gives a three dimensional impression of the
  surface, allows overplotting of the original data points
  (in case an RSA object is provided), and allows the
  interactive adjustment of regression weights in the
  \code{\link{demoRSA}} function. The interactive plot
  allows rotating and exploring a three-dimensional surface
  with the mouse (nice for demonstration purposes). If you
  want to export publication-ready plots, it is recommended
  to export it with following commands: \code{p1 <-
  plot(r1, bw=TRUE) trellis.device(device="cairo_pdf",
  filename="RSA_plot.pdf") print(p1) dev.off()}
}
\examples{
# Plot response surfaces from known parameters
# example of Edwards (2002), Figure 3
plotRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628)
plotRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628, type="c")

# Plot response surface from an RSA object
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
	z.diff <- diff + rnorm(n, 0, err)
	z.abs <- absdiff + rnorm(n, 0, err)
	z.sq <- sqdiff + rnorm(n, 0, err)
	z.add <- diff + 0.4*x + rnorm(n, 0, err)
	z.complex <- 0.4*x + - 0.2*x*y + + 0.1*x^2 - 0.03*y^2 + rnorm(n, 0, err)
})

r1 <- RSA(z.sq~x*y, df, models=c("sqdiff", "full", "IA"))
plot(r1)	# default: model = "full"
plot(r1, model="sqdiff", points=list(show=TRUE, value="predicted"))
}
\references{
  Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The
  Bagplot: A Bivariate Boxplot. The American Statistician,
  53(4), 382-387. doi:10.1080/00031305.1999.10474494
}
\seealso{
  \code{\link{demoRSA}}, \code{\link{RSA}}
}

