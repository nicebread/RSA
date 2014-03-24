\name{RSA}
\title{News for package 'RSA'}
\encoding{UTF-8}


\section{Version 0.9.1}{


\subsection{New features}{
	\itemize{
		\item New models: onlyx and onlyx2, resp. onlyy and onlyy2, test whether one of both predictors can be dropped altogether
		\item New parameter \code{add} in \code{RSA} function: The user can pass additional lavaan syntax that is added at the end of the generated syntax. For example: \code{RSA(z ~ x*y, data=df, add="abs(a4) > abs(a2); a1==0")}
		\item Added the \code{compare2} function which compares two specific models, e.g. \code{compare2(r1, m1="SRR", m2="full")}. Attention: You must take care yourself that the compared models are nested! There is no automatic check, so you could, in principle, compare non-nested models. This is valid for AIC, BIC, CFI, TLI, and R2 indices, but *not* for the chi2-LR test!
		\item plotRSA gained a \code{main} parameter: The main title of the plot.
		\item New generic functions:
		\itemize{
			\item \code{fitted}: Return fitted values of a model, e.g.: \code{fitted(r1, model="RR")}
			\item \code{residuals} (or \code{resid}): Return residuals of a fitted model, e.g.: \code{resid(r1, model="RR")}. For example, use the \code{qqPlot} function from the \code{car} package to display a qq-plot for regression diagnostics: qqPlot(resid(r2))
		}
	}
}

\subsection{Improvements}{
	\itemize{
		\item RSA(..., cubic=TRUE) now works more reliable. getPar(x, model="cubic") now returns the cubic parameters x^3, y^3, xy^2, and x^2y.
		\item New function \code{summary}. Currently, it is simply an alias to \code{print}
		\item \code{cex} parameter now applies to all text output of plotRSA.
		\item If outliers are removed automatically, the \code{RSA} function prints the row numbers of the outliers.
		\item Plotting of outliers (with \code{points=list(out.mark=TRUE)}) is more reliable.
	}
}
}



\section{Version 0.9}{

\subsection{New Features}{
	\itemize{
		\item Two new models: "onlyx" (x + x^2) and "onlyy" (y + y^2). Also included in \code{compare()} and in \code{compare(... , plot=TRUE)}
		\item Arbitrary data points can be plotted in \code{plotRSA} - i.e., you can plot a fit with existing regression weights, and still plot raw data points. See \code{ plotRSA(..., points(data=df))}
		\item \code{plotRSA(..., points(value="predict"))} now also works in interactive plots \code{(type="i")}
		\item The hull in \code{plotRSA} now is much smoother and lies on the surface (before, long connection lines cut through the surface)
		\item \code{print.RSA} now prints the regression coefficients including standardized beta weights, and can print summaries for a submodel: \code{print(r1, model="RR")}. Nicer output with significance stars, formatted \emph{p}-values, etc.
	}
}

\subsection{Improvements}{
	\itemize{
		\item Refactored \code{plotRSA} and \code{plot.RSA} - much cleaner now
	}
}
}



\section{Version 0.8.2}{

\subsection{New Features}{
	\itemize{
		\item If outliers are removed, this is now printed as a warning (more salient)
		\item added parameter "out.mark" in the points parameter of plotRSA: If set to TRUE, outliers according to Bollen & Jackman (1980) are printed as a red X symbol. This option works regardless of whether the RSA function has set out.rm to TRUE or FALSE:
		\itemize{
			\item If \code{out.rm == TRUE} (in \code{RSA}) and \code{out.mark == FALSE} (in \code{plotRSA}), the outlier is removed from the model and *not plotted* in \code{plotRSA}. 
			\item If \code{out.rm == TRUE} (in RSA) and \code{out.mark == TRUE} (in \code{plotRSA}), the outlier is removed from the model but plotted and marked in \code{plotRSA}.
			\item If \code{out.rm == FALSE} (in RSA) and \code{out.mark == FALSE} (in \code{plotRSA}), the outlier is not removed from the model and plotted as a normal point in \code{plotRSA} (but not marked as outlier). 
			\item If \code{out.rm == FALSE} (in RSA) and \code{out.mark == TRUE} (in \code{plotRSA}), the outlier is not removed from the model, but plotted and marked in \code{plotRSA}.
			\item Example syntax: \code{plotRSA(r1, points=list(show=TRUE, out.mark=TRUE))}
		}
	}
}

}


\section{Version 0.8}{

\subsection{New Features}{
	\itemize{
		\item replaced \code{bootRSA} and \code{boot.CI} functions by a unified interface: \code{confint()}. This is now standard R behavior. Confidence intervals can be extracted directly from the lavaan model, or bootstrapped. See \code{?confint}
		\item New sanity check: false/implausible convergence of constrained models is detected and a warning is printed.
		\item Changed parameters for points and contours in \code{plotRSA}: now is a list with subelements, e.g.: points=list(show=TRUE, color="red", jitter=.2)
	}
}

\subsection{Improvements}{
	\itemize{
		\item General R 3.0.2 compatibility
		\item Improved non-linear constraints; now false convergences are more uncommon
	}
}
}



\section{Version 0.7.4}{

\subsection{New Features}{
	\itemize{
		\item new parameter in \code{plotRSA}: `points.value`. If set to "raw", the original response is splotted, if set to "predicted", the predicted value is plotted.
		\item new parameter in \code{plotRSA}: `jitter`. The amount of jitter of the raw data points (to avoid overplotting). For standardized predictors, values around .02 make sense.
		\item new parameter in \code{plotRSA}: `contour.highlight`. Provide a vector of values, which indicate contour lines that should be plotted in bold.
		\item new option in \code{getPar}: "r2.adj"
		\item added the null model (= intercept-only model) to \code{RSA} and \code{compare} (the null model can also be plotted and only shows the flat intercept surface)
		\item Added a graphical model comparison: \code{plot(compare(r1))}, or \code{compare(r1, plot=TRUE)} (requires package \code{qgraph})
		\item new parameter in \code{plotRSA}: `pal.range`. If set to "box", the color ramp in the 3d and contour plots is scaled to the zlim, if set to "surface" it is scaled to the range of the actual surface. If set to "box", different 3d surface plots can be compared along their color, as long as the zlim is the same for both.
		\item New parameter in \code{plotRSA}: `maxlines`. Plots the maximum X for a given Y and vice versa.
		\item New models:
			\itemize{
				\item added the "shifted rising ridge model" (SRR)
				\item added the "shifted and rotated rising ridge model" (SRRR)
			}
	}
}

\subsection{Improvements}{
	\itemize{
		\item Suppressed another (uncritical) lavaan warning: NaN in sqrt(b3 * b5).
		\item Improved \code{modeltree} plot: non-significant models (with R2.p > .05) are printed in gray
		\item Moved plyr, gridExtra, and RColorBrewer into the Imports section
	}
}
}


\section{Version 0.7}{

\subsection{New Features}{
	\itemize{
		\item Added the "rising ridge" model (label = "RR")
		\item Added a movie function \code{movieRSA}: Create rotating and surface changing videos of surface plots! (this function needs ffmpeg to be accessible from the command line. It works on Mac OS, not sure about Windows or Linux!)
	}
}

\subsection{Changes}{
	\itemize{
		\item The color ramp in the 3d and contour plots now is scaled to the zlim, not to the range of the actual surface (as before). That means: Different surface plots can be compared along their color, as long as the zlim is the same for both.
		\item Moved packages tkrplot and rgl into "Suggests"
		\item Renamed parameter "showSP" to "param" in the \code{plotRSA} function
		\item Renamed parameter "rot" to "rotation" and "label.rot" to "label.rotation" in \code{plotRSA}
		\item Changed default rotation of 2d plot (it's now more comparable with the contour plot)
	}
}
}


\section{Version 0.6}{

\subsection{New Features}{
	\itemize{
		\item Added parameter \code{cubic} to RSA function. If set to \code{TRUE}, also a cubic surface is fitted, adding four new parameters to the equation: X^3, Y^3, XY^2, and X^2Y
		\item Added parameter \code{project} to plotRSA. If \code{project=TRUE}, all lines which are defines by parameter \code{axes} are projected onto the floor of the cube.
		\item \code{plotRSA} now can plot a convex hull around the data points (set \code{hull=TRUE}); both in \code{type="3d"} and type=\code{"contour"}
		\item \code{plotRSA} now can plot the CI of the stationary point (set \code{SP.CI=TRUE}); only for \code{type="contour"}
		\item All surface parameters now are also calculated for IA and sqdiff models
		\item Added parameter \code{surface} to \code{plot.RSA}. \code{surface = "predict"} (default) plots the predicted surface. This is the known behavior of the function. \code{surface = "smooth"} plots a 2d-smoothed surface of the raw data (employs the \code{Tps} function of the \code{fields} package).
	}
}

\subsection{Changes}{
	\itemize{
		\item Changed parameter name \code{showBorder} to \code{border}
		\item Changed parameter name \code{showContour} to \code{contour}
	}
}
}