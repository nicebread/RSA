#' @title Plots a response surface of a polynomial equation of second degree with interactive controls
#'
#' @description
#' Plots an RSA object, or a response surface with specified parameters, with interactive controls for coefficients.
#'
#' @details
#' No details so far. Just play around with the interface!
#'
#' @aliases demoSRR demoSRRR
#'
#' @export
#' @param x Either an RSA object (returned by the \code{RSA} function), or the coefficient for the X predictor
#' @param y Y coefficient
#' @param x2 X^2 coefficient
#' @param y2 Y^2 coefficient
#' @param xy XY interaction coefficient
#' @param y3 Y^3 coefficient
#' @param x3 X^3 coefficient
#' @param x2y X^2Y coefficient
#' @param xy2 XY^2 coefficient
#' @param w W coefficient (for (un)constrained absolute difference model)
#' @param wx WX coefficient (for (un)constrained absolute difference model)
#' @param wy WY coefficient (for (un)constrained absolute difference model)
#' @param b0 Intercept
#' @param xlim Limits of the x axis
#' @param ylim Limits of the y axis
#' @param zlim Limits of the z axis
#' @param xlab Label of the x axis
#' @param ylab Label of the y axis
#' @param zlab Label of the z axis

#' @param type \code{3d} for 3d surface plot, \code{contour} for 2d contour plot. Shortcuts (i.e., first letter of string) are sufficient; be careful: "contour" is very slow at the moment
#' @param points A list of parameters which define the appearance of the raw scatter points: show = TRUE: Should the original data points be overplotted? value="raw": Plot the original z value, "predicted": plot the predicted z value. jitter=0: Amount of jitter for the raw data points. cex = .5: multiplication factor for point size. See ?plotRSA for details.
#' @param project Which features should be projected on the floor? See ?plotRSA for details.
#' @param model If x is an RSA object: from which model should the response surface be computed?
#' @param extended Show additional controls (not implemented yet)
#' @param ... Other parameters passed through to plot.RSA (e.g., xlab, ylab, zlab, cex, legend)
#'
#'
#' @seealso \code{\link{plotRSA}}, \code{\link{RSA}}
#'
#' @examples
#' # Plot response surfaces from known parameters
#' # example of Edwards (2002), Figure 3
#' \dontrun{
#' demoRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628, type="3d")
#' demoRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628, legend=FALSE, type="c")
#' }
#'
#' # Plot response surface from an RSA object
#' \dontrun{
#' set.seed(0xBEEF)
#' n <- 300
#' err <- 2
#' x <- rnorm(n, 0, 5)
#' y <- rnorm(n, 0, 5)
#' df <- data.frame(x, y)
#' df <- within(df, {
#' 	diff <- x-y
#' 	absdiff <- abs(x-y)
#' 	SD <- (x-y)^2
#' 	z.diff <- diff + rnorm(n, 0, err)
#' 	z.abs <- absdiff + rnorm(n, 0, err)
#' 	z.sq <- SD + rnorm(n, 0, err)
#' 	z.add <- diff + 0.4*x + rnorm(n, 0, err)
#' 	z.complex <- 0.4*x + - 0.2*x*y + + 0.1*x^2 - 0.03*y^2 + rnorm(n, 0, err)
#' })
#' 
#' r1 <- RSA(z.sq~x*y, df)
#' demoRSA(r1)
#' demoRSA(r1, points=TRUE, model="SQD")
#' }

## TODO: Convert to Shiny app

#x=NULL; y=0; x2=0; y2=0; xy=0; w=0; wx=0; wy=0; x3=0; xy2=0; x2y=0; y3=0; b0=0; type="3d"; zlim=c(-2, 2); xlim=c(-2, 2); ylim=c(-2, 2); xlab=NULL; ylab=NULL; zlab=NULL; points = TRUE; model="full"; project=c("PA1", "PA2"); extended=FALSE

demoRSA <- function(x=NULL, y=0, x2=0, y2=0, xy=0, w=0, wx=0, wy=0, x3=0, xy2=0, x2y=0, y3=0, b0=0, type="3d", zlim=c(-2, 2), xlim=c(-2, 2), ylim=c(-2, 2), xlab=NULL, ylab=NULL, zlab=NULL, points = TRUE, model="full", project=c("PA1", "PA2"), extended=FALSE, ...) {

	type <- match.arg(type, c("interactive", "3d", "contour"))
	type2 <- type
	if (type2 == "interactive") stop("demoRSA only works with type == '3d' or 'contour'!")
		
	if (!requireNamespace("tkrplot", quietly = TRUE)) {
		stop('`tkrplot` package needed for demoRSA. Please install with install.packages("tkrplot")', call. = FALSE)
	}	
	

	# if model is provided: take its parameters as starting values
	if (is.null(x)) {
		x <- 0
		fit <- NULL
		points <- list(data=NULL, show=NA, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
		if (is.null(xlab)) {xlab <- "X"}
		if (is.null(ylab)) {ylab <- "Y"}
		if (is.null(zlab)) {zlab <- "Z"}
			
	} else if (!is.null(x) & !is.null(attr(x, "class"))) {
		if (attr(x, "class") == "RSA") {
			fit <- x
			C <- coef(fit$models[[model]])
			b0.0 <- b0 <- as.numeric(ifelse(is.na(C[paste0(fit$DV, "~1")]), b0, C[paste0(fit$DV, "~1")]))
			x.0 <- x <- as.numeric(ifelse(is.na(C["b1"]), 0, C["b1"]))
			y.0 <- y <- as.numeric(ifelse(is.na(C["b2"]), y, C["b2"]))
			x2.0 <- x2 <- as.numeric(ifelse(is.na(C["b3"]), x2, C["b3"]))
			y2.0 <- y2 <- as.numeric(ifelse(is.na(C["b5"]), y2, C["b5"]))
			xy.0 <- xy <- as.numeric(ifelse(is.na(C["b4"]), xy, C["b4"]))
			w.0 <- w <- as.numeric(ifelse(is.na(C["w1"]), w, C["w1"]))
			wx.0 <- wx <- as.numeric(ifelse(is.na(C["w2"]), wx, C["w2"]))
			wy.0 <- wy <- as.numeric(ifelse(is.na(C["w3"]), wy, C["w3"]))
			
			x3.0 <- x3 <- as.numeric(ifelse(is.na(C["b6"]), x3, C["b6"]))
			x2y.0 <- x2y <- as.numeric(ifelse(is.na(C["b7"]), x2y, C["b7"]))
			xy2.0 <- xy2 <- as.numeric(ifelse(is.na(C["b8"]), xy2, C["b8"]))
			y3.0 <- y3 <- as.numeric(ifelse(is.na(C["b9"]), y3, C["b9"]))
		
			xlim <- c(min(fit$data[, fit$IV1], na.rm=TRUE), max(fit$data[, fit$IV1], na.rm=TRUE))
			ylim <- c(min(fit$data[, fit$IV2], na.rm=TRUE), max(fit$data[, fit$IV2], na.rm=TRUE))
			
			if (is.null(xlab)) xlab <- fit$IV1
			if (is.null(ylab)) ylab <- fit$IV2
			if (is.null(zlab)) zlab <- fit$DV
				
			# expand range by 20% at each end
			xlim[1] <- xlim[1]*ifelse(xlim[1]<0, 1.1, 0.9)
			xlim[2] <- xlim[2]*ifelse(xlim[2]<0, 0.9, 1.1)
			ylim[1] <- ylim[1]*ifelse(ylim[1]<0, 1.1, 0.9)
			ylim[2] <- ylim[2]*ifelse(ylim[2]<0, 0.9, 1.1)
				
			# for the correct visual diagonal: same range for X and Y
			xlim[1] <- ylim[1] <- min(xlim[1], ylim[1])
			xlim[2] <- ylim[2] <- max(xlim[2], ylim[2])
		
			zlim <- c(min(fit$data[, fit$DV], na.rm=TRUE)*0.8, max(fit$data[, fit$DV], na.rm=TRUE)*1.2)
			
			# define the defaults
			if (is.null(points) || (typeof(points) == "logical" && points == TRUE)) {
				points <- list(show=TRUE, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
			}
			if (is.null(points) || (typeof(points) == "logical" && points == FALSE)) {
				points <- list(show=FALSE, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
			}
			if (is.null(points$out.mark)) points$out.mark <- FALSE

			if (points$out.mark == FALSE) {data.used <- fit$data[fit$data$out==FALSE, ]}
			if (points$out.mark == TRUE) {data.used <- fit$data}
			
			points$data <- data.used[, c(fit$IV1, fit$IV2, fit$DV, colnames(fit$data)[which(!colnames(fit$data) %in% c(fit$IV1, fit$IV2, fit$DV))])]
		}
	} else {
		fit <- NULL
		points <- list(data=NULL, show=NA, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
		if (is.null(xlab)) {xlab <- "X"}
		if (is.null(ylab)) {ylab <- "Y"}
		if (is.null(zlab)) {zlab <- "Z"}
	}
	

  TYPE <- tcltk::tclVar(); tcltk::tclvalue(TYPE) <- "full"
	
	B0 <- tcltk::tclVar(); tcltk::tclvalue(B0) <- b0
	X <- tcltk::tclVar(); tcltk::tclvalue(X) <- x
	Y <- tcltk::tclVar(); tcltk::tclvalue(Y) <- y
	X2 <- tcltk::tclVar(); tcltk::tclvalue(X2) <- x2
	Y2 <- tcltk::tclVar(); tcltk::tclvalue(Y2) <- y2
	XY <- tcltk::tclVar(); tcltk::tclvalue(XY) <- xy
	
	W <- tcltk::tclVar(); tcltk::tclvalue(W) <- w
	WX <- tcltk::tclVar(); tcltk::tclvalue(WX) <- wx
	WY <- tcltk::tclVar(); tcltk::tclvalue(WY) <- wy
	
	# rotation of the 3d-frame
	RX <- tcltk::tclVar(); tcltk::tclvalue(RX) <- -45
	RY <- tcltk::tclVar(); tcltk::tclvalue(RY) <- 45
	RZ <- tcltk::tclVar(); tcltk::tclvalue(RZ) <- 35
	
	# Dummy variables: Shift and rotation
	C <- tcltk::tclVar(); tcltk::tclvalue(C) <- 0
	S <- tcltk::tclVar(); tcltk::tclvalue(S) <- 1
	
	if (extended==TRUE) {
		X.Y2 <- tcltk::tclVar(); tcltk::tclvalue(X.Y2) <- 0
	}
	
	
	setAllBlack <- function() {
		sapply(list(X.lab, Y.lab, X2.lab, Y2.lab, XY.lab, W.lab, WX.lab, WY.lab), tcltk::tkconfigure, foreground="black")
	}

	update <- function(...) {
				
		# read parameters from sliders
        type <- as.character(tcltk::tclvalue(TYPE))
		b0 <- as.numeric(tcltk::tclvalue(B0))
		x <- as.numeric(tcltk::tclvalue(X))
		y <- as.numeric(tcltk::tclvalue(Y))
		x2 <- as.numeric(tcltk::tclvalue(X2))
		y2 <- as.numeric(tcltk::tclvalue(Y2))
		xy <- as.numeric(tcltk::tclvalue(XY))
		w <- as.numeric(tcltk::tclvalue(W))
		wx <- as.numeric(tcltk::tclvalue(WX))
		wy <- as.numeric(tcltk::tclvalue(WY))
		rx <- as.numeric(tcltk::tclvalue(RX))
		ry <- as.numeric(tcltk::tclvalue(RY))
		rz <- as.numeric(tcltk::tclvalue(RZ))
		c <- as.numeric(tcltk::tclvalue(C))
		s <- as.numeric(tcltk::tclvalue(S))
		if (extended==TRUE) {x.y2 <- as.numeric(tcltk::tclvalue(X.Y2))}
        
		setAllBlack()
		
		# set constraints
		if (type == "all") {
		}
		if (type == "poly") {
			tcltk::tclvalue(W) <- tcltk::tclvalue(WX) <- tcltk::tclvalue(WY) <- 0
			sapply(list(W.lab, WX.lab, WY.lab), tcltk::tkconfigure, foreground="grey40")
		}
		
		if (type == "diff") {
			tcltk::tclvalue(Y) <- -x
			tcltk::tclvalue(X2) <- tcltk::tclvalue(Y2) <- tcltk::tclvalue(XY) <- tcltk::tclvalue(W) <- tcltk::tclvalue(WX) <- tcltk::tclvalue(WY) <- 0
			sapply(list(X2.lab, Y2.lab, XY.lab, W.lab, WX.lab, WY.lab), tcltk::tkconfigure, foreground="grey40")
		}
		if (type == "SQD") {
			tcltk::tclvalue(Y) <- 0
			tcltk::tclvalue(X) <- 0
			tcltk::tclvalue(Y2) <- x2
			tcltk::tclvalue(XY) <- -2*x2
			tcltk::tclvalue(W) <- tcltk::tclvalue(WX) <- tcltk::tclvalue(WY) <- 0
			sapply(list(X.lab, Y.lab, Y2.lab, XY.lab, W.lab, WX.lab, WY.lab), tcltk::tkconfigure, foreground="grey40")
		}
		if (type == "sq.shift") {
			tcltk::tclvalue(Y) <- -x
			tcltk::tclvalue(Y2) <- x2
			tcltk::tclvalue(XY) <- -2*x2
			if (x2 != 0) {
				tcltk::tclvalue(B0) <- x^2 / (4*x2)				
				tcltk::tclvalue(C) <- x/(2*x2)
			}
			tcltk::tclvalue(W) <- tcltk::tclvalue(WX) <- tcltk::tclvalue(WY) <- 0
			sapply(list(Y.lab, Y2.lab, XY.lab, W.lab, WX.lab, WY.lab), tcltk::tkconfigure, foreground="grey40")
		}
		if (type == "sq.rot") {
			
			if (y2 != 0) {
				tcltk::tclvalue(X2) <- (xy^2) / (4*y2)
				x <- tcltk::tclvalue(X) <- (y*xy)/(2*y2)
			}
			
			if (y2 != 0 & y != 0) {
				tcltk::tclvalue(B0) <- y^2 / (4*y2)
				tcltk::tclvalue(C) <- -0.5*(y/y2)
				tcltk::tclvalue(S) <- -(x/y)
			}
			tcltk::tclvalue(W) <- tcltk::tclvalue(WX) <- tcltk::tclvalue(WY) <- 0
			sapply(list(X.lab, X2.lab, W.lab, WX.lab, WY.lab), tcltk::tkconfigure, foreground="grey40")
		}
		if (type == "IA") {
			tcltk::tclvalue(X2) <- 0
			tcltk::tclvalue(Y2) <- 0
			tcltk::tclvalue(W) <- tcltk::tclvalue(WX) <- tcltk::tclvalue(WY) <- 0
			sapply(list(X2.lab, Y2.lab, W.lab, WX.lab, WY.lab), tcltk::tkconfigure, foreground="grey40")
		}
		if (type == "absunc") {
			tcltk::tclvalue(X2) <- 0
			tcltk::tclvalue(Y2) <- 0
			tcltk::tclvalue(XY) <- 0
			sapply(list(X2.lab, Y2.lab, XY.lab), tcltk::tkconfigure, foreground="grey40")
		}
		if (type == "absdiff") {
			tcltk::tclvalue(X) <- 0
			tcltk::tclvalue(Y) <- 0
			tcltk::tclvalue(X2) <- 0
			tcltk::tclvalue(Y2) <- 0
			tcltk::tclvalue(XY) <- 0
			tcltk::tclvalue(W) <- 0
			tcltk::tclvalue(WY) <- -wx
			sapply(list(X.lab, Y.lab, X2.lab, Y2.lab, XY.lab, W.lab, WY.lab), tcltk::tkconfigure, foreground="grey40")
		}

    tkrplot::tkrreplot(img, hscale=1.5, vscale=1.5)
   }

    replot <- function() {
			# read parameters from sliders
	    type <- as.character(tcltk::tclvalue(TYPE))
			b0 <- as.numeric(tcltk::tclvalue(B0))
			x <- as.numeric(tcltk::tclvalue(X))
			y <- as.numeric(tcltk::tclvalue(Y))
			x2 <- as.numeric(tcltk::tclvalue(X2))
			y2 <- as.numeric(tcltk::tclvalue(Y2))
			xy <- as.numeric(tcltk::tclvalue(XY))
			w <- as.numeric(tcltk::tclvalue(W))
			wx <- as.numeric(tcltk::tclvalue(WX))
			wy <- as.numeric(tcltk::tclvalue(WY))
			rx <- as.numeric(tcltk::tclvalue(RX))
			ry <- as.numeric(tcltk::tclvalue(RY))
			rz <- as.numeric(tcltk::tclvalue(RZ))
		
			plot(plotRSA(x=x, y=y, x2=x2, y2=y2, xy=xy, w=w, wx=wx, wy=wy, b0=b0, rotation=list(x=rx, y=ry, z=rz), zlim=zlim, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, zlab=zlab, points=points, demo=TRUE, type=type2, fit=fit, project=project))
    }

	# define framework
    tt <- tcltk::tktoplevel()
    tcltk::tkwm.title(tt, "Response surface plot - polynomial model")

    img <- tkrplot::tkrplot(tt, replot, vscale=1.5, hscale=1.5)
    tcltk::tkpack(img, side='left')
	
	# define radiobuttons
	tcltk::tkpack(tfr <- tcltk::tkframe(tt, relief='groove', borderwidth=3), side='top')
	
	tcltk::tkpack(typebox <- tcltk::tkframe(tfr), side='top', fill='x')
    tcltk::tkpack(tcltk::tklabel(typebox,text='Constraints: '), side='left',anchor='s')
	tcltk::tkpack(tcltk::tkradiobutton(typebox, variable=TYPE, command=update, value="all", text="All parameters"))
	tcltk::tkpack(tcltk::tkradiobutton(typebox, variable=TYPE, command=update, value="poly", text="Full polynomial"))
	tcltk::tkpack(tcltk::tkradiobutton(typebox, variable=TYPE, command=update, value="IA", text="Interaction"))
	tcltk::tkpack(tcltk::tkradiobutton(typebox, variable=TYPE, command=update, value="SQD", text="Squared difference"))
	tcltk::tkpack(tcltk::tkradiobutton(typebox, variable=TYPE, command=update, value="sq.shift", text="Shifted squared difference"))
	tcltk::tkpack(tcltk::tkradiobutton(typebox, variable=TYPE, command=update, value="sq.rot", text="Shifted and rotated squared difference"))
	tcltk::tkpack(tcltk::tkradiobutton(typebox, variable=TYPE, command=update, value="diff", text="Difference score X-Y"))
	
	tcltk::tkpack(tcltk::tkradiobutton(typebox, variable=TYPE, command=update, value="absunc", text="Unconstrained absolute difference"))
	tcltk::tkpack(tcltk::tkradiobutton(typebox, variable=TYPE, command=update, value="absdiff", text="Absolute difference"))

	

	# define sliders: polynomial model
	tcltk::tkpack(tfr <- tcltk::tkframe(tt, relief='groove', borderwidth=3), side='left')
	tcltk::tkpack(fr0 <- tcltk::tkframe(tfr), side='top',fill='x')
	tcltk::tkpack(fr1 <- tcltk::tkframe(tfr), side='top',fill='x')
	tcltk::tkpack(fr2 <- tcltk::tkframe(tfr), side='top',fill='x')
	tcltk::tkpack(fr3 <- tcltk::tkframe(tfr), side='top',fill='x')
	tcltk::tkpack(fr4 <- tcltk::tkframe(tfr), side='top',fill='x')
	tcltk::tkpack(fr5 <- tcltk::tkframe(tfr), side='top',fill='x')
	tcltk::tkpack(fr6 <- tcltk::tkframe(tfr), side='top',fill='x')
	tcltk::tkpack(fr7 <- tcltk::tkframe(tfr), side='top',fill='x')
	tcltk::tkpack(fr8 <- tcltk::tkframe(tfr), side='top',fill='x')
	B0.lab <- tcltk::tklabel(fr0,text='Intercept: ')
	X.lab <- tcltk::tklabel(fr1,text='x: ')
	Y.lab <- tcltk::tklabel(fr2,text='y: ')
	XY.lab <- tcltk::tklabel(fr3,text='xy: ')
	X2.lab <- tcltk::tklabel(fr4,text='x2: ')
	Y2.lab <- tcltk::tklabel(fr5,text='y2: ')
	W.lab <- tcltk::tklabel(fr6,text='w: ')
	WX.lab <- tcltk::tklabel(fr7,text='wx: ')
	WY.lab <- tcltk::tklabel(fr8,text='wy: ')
	
    tcltk::tkpack(B0.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr0, variable=B0, orient='horizontal', command=update, from=-5, to=5, resolution=.1), side='left')
	
    tcltk::tkpack(X.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr1, variable=X, orient='horizontal', command=update, from=ifelse(is.null(fit), -5, -abs(x.0)*2), to=ifelse(is.null(fit), 5, abs(x.0)*2), resolution=0.01), side='left')

    tcltk::tkpack(Y.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr2, variable=Y, orient='horizontal', command=update, from=ifelse(is.null(fit), -5, -abs(y.0)*2), to=ifelse(is.null(fit), 5, abs(y.0)*2), resolution=0.01), side='left')
	
    tcltk::tkpack(XY.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr3, variable=XY, orient='horizontal', command=update, from=ifelse(is.null(fit), -3, -abs(xy.0)*2), to=ifelse(is.null(fit), 3, abs(xy.0)*2), resolution=0.01), side='left')
	
    tcltk::tkpack(X2.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr4, variable=X2, orient='horizontal', command=update, from=ifelse(is.null(fit), -3, -abs(x2.0)*2), to=ifelse(is.null(fit), 3, abs(x2.0)*2), resolution=0.01), side='left')

    tcltk::tkpack(Y2.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr5, variable=Y2, orient='horizontal', command=update, from=ifelse(is.null(fit), -3, -abs(y2.0)*2), to=ifelse(is.null(fit), 3, abs(y2.0)*2), resolution=0.01), side='left')
	
	# define sliders: absdiff model
	tcltk::tkpack(tfr <- tcltk::tkframe(tt, relief='groove', borderwidth=3), side='right')
	
    tcltk::tkpack(W.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr6, variable=W, orient='horizontal', command=update, from=ifelse(is.null(fit), -5, -abs(w.0)*2), to=ifelse(is.null(fit), 5, abs(w.0)*2), resolution=0.01), side='left')

    tcltk::tkpack(WX.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr7, variable=WX, orient='horizontal', command=update, from=ifelse(is.null(fit), -1, -abs(wx.0)*2), to=ifelse(is.null(fit), 1, abs(wx.0)*2), resolution=0.01), side='left')
	
    tcltk::tkpack(WY.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr8, variable=WY, orient='horizontal', command=update, from=ifelse(is.null(fit), -1, -abs(wy.0)*2), to=ifelse(is.null(fit), 1, abs(wy.0)*2), resolution=0.01), side='left')
	
	
	## Rotation of display
	tcltk::tkpack(tfr3d <- tcltk::tkframe(tt, relief='groove', borderwidth=3), side='right')
	tcltk::tkpack(fr3.1 <- tcltk::tkframe(tfr3d), side='top',fill='x')
	tcltk::tkpack(fr3.2 <- tcltk::tkframe(tfr3d), side='top',fill='x')
	tcltk::tkpack(fr3.3 <- tcltk::tkframe(tfr3d), side='top',fill='x')
	X3.lab <- tcltk::tklabel(fr3.1,text='x rotation: ')
	Y3.lab <- tcltk::tklabel(fr3.2,text='y rotation: ')
	Z3.lab <- tcltk::tklabel(fr3.3,text='z rotation: ')
    tcltk::tkpack(X3.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr3.1, variable=RX, orient='horizontal', command=update, from=-90, to=90, resolution=1), side='left')

    tcltk::tkpack(Y3.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr3.2, variable=RY, orient='horizontal', command=update, from=-90, to=90, resolution=1), side='left')
	
    tcltk::tkpack(Z3.lab, side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(fr3.3, variable=RZ, orient='horizontal', command=update, from=-90, to=90, resolution=1), side='left')
	
	
	
	## Extra (dummy) parameters
	tcltk::tkpack(frROT1 <- tcltk::tkframe(tfr3d), side='top',fill='x')
	tcltk::tkpack(frROT2 <- tcltk::tkframe(tfr3d), side='top',fill='x')
	
    tcltk::tkpack(tcltk::tklabel(frROT1,text='Shift (C): '), side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(frROT1, variable=C, orient='horizontal', command=update, from=-20, to=20, resolution=0.1), side='left')

    tcltk::tkpack(tcltk::tklabel(frROT2,text='Rotation (S): '), side='left',anchor='s')
	tcltk::tkpack(tcltk::tkscale(frROT2, variable=S, orient='horizontal', command=update, from=0, to=3, resolution=0.1), side='left')
	
    return(invisible(NULL))
}


#demoRSA()
#demoRSA(fit=r1, points=TRUE)

# Hack to please CRAN:
if(getRversion() >= "2.15.1")  {
	utils::globalVariables(c('tclVar', 'tclvalue', 'tkconfigure' , 'tkframe', 'tklabel', 'tkpack', 'tkradiobutton', 'tkscale', 'tktoplevel', 'tkwm.title', "tkrplot", "tkrreplot"))
}