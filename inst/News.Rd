\name{RSA}
\title{News for package 'RSA'}
\encoding{UTF-8}

\section{Version 0.9.12}{


\subsection{New features}{
	\itemize{
		\item New parameter \code{axesStyles} in \code{plotRSA}: Define custom styles for LOC, LOIC, PA1, and PA2. Recognizes \code{lty}, \code{lwd}, and \code{col}. If you define a style for an axis, you have to provide all three parameters, otherwise a warning will be shown.
		\item The cubic models CA, RRCA, CL, and RRCL have been implemented (see Humberg, Nestler, Schönbrodt, & Back, in preparation, for details on how to use these models) and the plot() function has been adapted respectively. The new functions caRange() and clRange() have been implemented. They should be applied when testing cubic models to ensure that the model is interpreted for the whole range of realistic predictor combinations.
	}
}

\subsection{Bug fixes/Glitches}{
	\itemize{
		\item Changed the outlier detection, which was very liberal (i.e., flagging very few potential outliers). It now implements the recommendations from ... Please note: '[...] that the diagnostics cannot be employed mechanically. Instead, once a case is diagnosed as influential, remedial action requires a firm substantive grounding. Although case deletion may be warranted in some circumstances, it is an extreme remedy of last resort that should not be routinely followed. The more fruitful approach is to ask why a given case is influential' (Bollen & Jackman, 1985).
		\item The newly introduced a5 parameter was displayed as a second 'a1' in the 3d plot - this has been fixed.		
		\item When counting the rate of discrepant predictor combinatios, standardize the predictor difference at the predictors' grand standard deviation. This ensures that the percentage of (in)congruent predictors is independent of common linear transformations of the predictor variables. Accordingly, a predictor combination is now categorized as incongruent when X and Y differ by more than half a grand standard deviation from each other.
	}
}

}




\section{Version 0.9.11}{


\subsection{New features}{
	\itemize{
		\item New param in \code{plotRSA}: \code{suppress.grid}.
		\item If \code{link="probit"} in \code{plotRSA}, the 3d plot automatically scales the z-axis to c(0, 1)
		\item If \code{link="probit"} in \code{plotRSA}, the raw data points now can take the predicted probability (\code{plotRSA(..., points=list(value="predicted"))})
		\item \code{cex} params in plotRSA slightly changed: now we have three cex parameters to independently control the font size of the main title (\code{cex.main}), the axes labels (\code{cex.axesLabel}), and the tick labels (cex.tickLabel)
		\item New surface parameter for testing strict congruence patterns ("a5"). Thanks to Sarah Humberg for implementation.
	}
}

\subsection{Bug fixes}{
	\itemize{
		\item Fixed a wrong formula in the standard error of a4 (the formula was corrected in an erratum of Shanock, L. R., Baran, B. E., Gentry, W. A., & Pattison, S. C. (2014). Erratum to: Polynomial regression with response surface analysis: A powerful approach for examining moderation and overcoming limitations of difference scores. Journal of Business and Psychology, 29, 161. http://doi.org/10.1007/s10869-013-9317-6. Thanks to Sarah Humberg for pointing to the erratum).
		\item Fixed a bug in the curvature of the principal axes (closing issue 2 on Github; thanks to scith for reporting the bug).
	}
}

}


\section{Version 0.9.10}{



\subsection{Bug fixes}{
	\itemize{
		\item Fixed a wrong sign in the C parameter of the SRR model (thanks to Sarah Humberg for detecting the bug.)
	}
}

}


\section{Version 0.9.9}{


\subsection{New features and improvements}{
	\itemize{
		\item Better documentation for RSA.ST function.
	}
}

\subsection{Bug fixes}{
	\itemize{
		\item Made necessary changes due to lavaan update.
		\item Fixed issue due to ggplot2 update.
	}
}

}


\section{Version 0.9.8}{


\subsection{Bug fixes}{
	\itemize{
		\item Removed the beta version of the control.variables feature - this broke some functions.
		\item Removed some unnecessary warning that might confuse users.
	}
}

}


\section{Version 0.9.7}{


\subsection{New features and improvements}{
	\itemize{
		\item Improved aictab-plot; now also has a black/white option (\code{bw=TRUE})		
		\item Added evidence ratio to the \code{aictab} output: How much more likely is the best model (in first row) compared to each other model in the set?
		\item \code{compare}, \code{compare2}, and \code{aictab} gained a \code{digits} argument
		\item Reduced the set of models that is computed by default in the \code{RSA} function. Uncommon models, such as the \code{mean} model or the \code{diff} model have been removed.
		\item Removed eigenvalues from standard summary. Most users will never need them.
		\item When CIs and p values are bootstrapped, S is automatically tested against 1 (which is the null hypothesis).
		\item Improved contour plot (bag plot is now better visible; if xlim and ylim are set, the range is not extended anymore to include all points).
		\item New option in plotRSA: \code{showSP} (TRUE/FALSE) controls the display of the stationary point in contour plots, \code{showSP.CI} (TRUE/FALSE) controls the display of the stationary point's confidence interval.
		\item New demo data set: \code{motcon2}
	}
}

\subsection{Bug fixes}{
	\itemize{
		\item Fixed a bug that broke the \code{surface="smooth"} option in plotRSA
		\item Fixed a bug that \code{out.rm=FALSE} did not work any more.
		\item Fixed a bug where \code{getPar(..., type="r2.p")} returned a wrong p value (the output in \code{print} and model comparisons always was correct.)
		\item Fixed wrong label in \code{aictab} plot
	}
}

}



\section{Version 0.9.6}{

This version has been uploaded to CRAN.

\subsection{New features and improvements}{
	\itemize{
		\item Removed dependency \code{gridExtra}, added dependency \code{tcltk}
		\item Fixed a bug in \code{demoRSA}
		\item Added new demo: \code{demoSRRR}, which has surface parameteres as sliders.
	}
}
}


\section{Version 0.9.5}{

\subsection{New features and improvements}{
	\itemize{
		\item Robustness against non-normality: Default estimator now is "MLR", which provides robust standard errors, a robust scaled test statistic, and can handle missing values.
		\item Option \code{missing} in \code{RSA}. If missing values are present, by default \code{missing = "fiml"} is set, and a warning is issued that this is only valid if the data are missing completely at random (MCAR) or missing at random (MAR). Use \code{missing = "listwise"} to exclude missing values.
		\item In the 3d plot, addtionally data points and hull can be projected on the floor. Example: \code{plot(r1, project=c("PA1", "contour", "points", "hull"))}
		\item Use \code{pal="flip"} in plotRSA to reverse the default palette (green on bottom, red on top)
		\item Use \code{coefs=TRUE} in plotRSA to print the regression coefficients in the top right corner of a 3d plot.
		\item Added new demo: \code{demoSRR}, which has surface parameteres as sliders.
	}
}


\subsection{Bug fixes/ Improvements}{
	\itemize{
		\item Internal restructuring to please CRAN; less dependencies (more imports) etc.
		\item Added Travis CI to the Github repo
		\item Catched some border cases in plotRSA (e.g., a flat surface led to problems)
	}
}
}


\section{Version 0.9.4}{

\subsection{New features and improvements}{
	\itemize{
		\item The surface can be suppressed in the 3d plot: \code{suppress.surface=TRUE}. Useful for didactical purposes (first show empty cube, or cube with raw data points, then fade-in the regression surface.)
		\item \code{summary} now also reports the R^2 of the selected model and the intercept.
		\item Changed default for lavaan standard errors and test statistic to "robust"
	}
}
\subsection{Bug Fixes}{
	\itemize{
		\item \code{summary}: The standardized beta weights are now correct.
		\item \code{plotRSA}: The LOC and LOIC now are black in b/w mode
		\item \code{plotRSA}: When attempting to plot an RSA object with \code{plotRSA} (instead of simply \code{plot(r1)}), an informative warning is shown.
		\item New shortcut in \code{plotRSA}: Simply set \code{points=TRUE} (instead of defining the full list).
		\item Fixed a bug in \code{demoRSA} which prevented plotting of RSA-objects.
		\item Moved package \code{aplpack} into \code{Depends}, as plotting the hull now by default is turned on.
	}
}
}

\section{Version 0.9.3}{

\subsection{New features and improvements}{
	\itemize{
		\item Added function \code{compare2} to compare two specific models.
		\item \code{plotRSA}: parameter \code{project} can now plot axes ("LOC", "LOIC", "PA1", and "PA2") and the contour lines "contour" on the bottom of the cube. Example: \code{plot(r1, project=c("PA1", "contour"))}
		\item Significantly reduced memory footprint of RSA objects, some speed improvements.
		\item Plotting outliers is now more consistent.
		\item The hull is plotted by default (whenever raw data is provided, either in an RSA-fit object, or in the \code{points\$data} parameter)
	}
}
\subsection{Bug Fixes}{
	\itemize{
		\item Removed unnecessary warning when plotting an interactive plot.
	}
}
}

\section{Version 0.9.2}{

\subsection{User visible changes}{
	\itemize{
		\item Changed model abbreviations to be more consistent:
			\itemize{
				\item \code{sqdiff} to \code{SQD}
				\item \code{SSD} to \code{SSQD}
				\item \code{SRSD} to \code{SRSQD}
			}
	}
}

\subsection{New features}{
	\itemize{
		\item AICc (corrected AIC) is provided in addition in the \code{compare} function - this is better for small sample sizes
		\item An AIC model comparison table can be printed (based on code by Jarrett Byrnes).
	}
}

\subsection{Bug Fixes}{
	\itemize{
		\item Predicted value of stationary point now is reported in \code{summary}
	}
}
}




\section{Version 0.9.1}{


\subsection{New features}{
	\itemize{
		\item New models: onlyx and onlyx2, resp. onlyy and onlyy2, test whether one of both predictors can be dropped altogether. mean model: b1 == b2, i.e., x and y are averaged and taken as a single predictor.
		\item New parameter \code{add} in \code{RSA} function: The user can pass additional lavaan syntax that is added at the end of the generated syntax. For example: \code{RSA(z ~ x*y, data=df, add="abs(a4) > abs(a2); a1==0")}
		\item Added the \code{compare2} function which compares two specific models, e.g. \code{compare2(r1, m1="SRR", m2="full")}. Attention: You must take care yourself that the compared models are nested! There is no automatic check, so you could, in principle, compare non-nested models. This is valid for AIC, BIC, CFI, TLI, and R2 indices, but *not* for the chi2-LR test!
		\item plotRSA gained a \code{title} parameter: The main title of the plot.
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


\subsection{Bug Fixes}{
	\itemize{
		\item Removed the \code{error using packet 1 'x' and 'units' must have length > 0} bug in 3d plots
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
		\item All surface parameters now are also calculated for IA and SD models
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