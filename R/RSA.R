#' @title Performs several RSA model tests on a data set with two predictors
#'
#' @description
#' Performs several RSA model tests on a data set with two predictors
#'
#' @details
#' You can also fit binary outcome variables with a probit link function. For that purpose, the response variable has to be defined as "ordered": \code{r1 <- RSA(Z.binary ~ X*Y, dat, ordered="Z.binary")} (for more details see the help file of the \code{sem} function in the \code{lavaan} package.). The results can also be plotted with probabilities on the z axis using the probit link function: \code{plot(r1, link="probit", zlim=c(0, 1), zlab="Probability")}. \code{lavaan} at the moment only supports a probit link function for binary outcomes, not a logit link.
#'
#' @export
#' @import lavaan
#' @import ggplot2
#' @import gridExtra
#' @import lattice
#' @import RColorBrewer
#' @param formula A formula in the form \code{z ~ x*y}, specifying the variable names used from the data frame, where z is the name of the response variable, and x and y are the names of the predictor variables.
#' @param data A data frame with the variables
#' @param center Should predictor variables be centered on \emph{each variable's} sample mean before analyses? You should think carefully about this option, as different centering of the predictor variables can affect the commensurability of the predictor scales.
#' @param scale Should predictor variables be scales on the SD of \emph{each variable} before analyses? You should think carefully about this option, as different scaling of the predictor variables can affect the commensurability of the predictor scales.
#' @param na.rm Remove missings before proceeding?
#' @param add Additional syntax that is added to the lavaan model. Can contain, for example, additional constraints, like "p01 == 0; p11 == 0"
#' @param out.rm Should outliers according to Bollen & Jackman (1980) criteria be excluded from analyses?
#' @param breakline Should the breakline in the unconstrained absolute difference model be allowed (the breakline is possible from the model formulation, but empirically rather unrealistic ...). Defaults to \code{FALSE}
#' @param verbose Should additional information during the computation process be printed?
#' @param models A vector with names of all models that should be computed. Should be any from \code{c("absdiff", "absunc", "diff", "additive", "IA", "sqdiff", "RR", "SSD", "SRSD", "full", "null", "onlyx", "onlyy")}. For \code{models="all"}, all models are computed, for \code{models="default"} all models besides absolute difference models are computed.
#' @param cubic Should a cubic model with the additional terms Y^3, XY^2, YX^2, and X^3 be included?
#' @param control.variables A string vector with variable names from \code{data}. These variables are added as linear predictors to the model (in order "to control for them"). No interactions with the other variables are modeled.
#' @param ... Additional parameters passed to the \code{lavaan} \code{\link{sem}} function. For example, you can obtained bootstrapped standard errors by setting \code{se="boot"}.
#'
#'
#' @seealso \code{\link{demoRSA}}, \code{\link{plotRSA}}, \code{\link{RSA.ST}}, \code{\link{confint.RSA}}
#'
#' @examples
#' # Compute response surface from a fake data set
#' set.seed(0xBEEF)
#' n <- 300
#' err <- 15
#' x <- rnorm(n, 0, 5)
#' y <- rnorm(n, 0, 5)
#' df <- data.frame(x, y)
#' df <- within(df, {
#' 	diff <- x-y
#' 	absdiff <- abs(x-y)
#' 	sqdiff <- (x-y)^2
#' 	z.diff <- diff + rnorm(n, 0, err)
#' 	z.abs <- absdiff + rnorm(n, 0, err)
#' 	z.sq <- sqdiff + rnorm(n, 0, err)
#' 	z.add <- diff + 0.4*x + rnorm(n, 0, err)
#' 	z.complex <- 0.4*x + - 0.2*x*y + + 0.1*x^2 - 0.03*y^2 + rnorm(n, 0, err)
#' })
#' \dontrun{
#' r1 <- RSA(z.sq~x*y, df)
#' print(r1)
#' compare(r1)
#' plot(r1)
#' plot(r1, model="SRSD")
#' plot(r1, model="full", type="c")
#' getPar(r1, "coef")	# print model parameters including SE and CI
#' RSA.ST(r1)	# get surface parameters
#'
#' # Motive congruency example
#' data(motcon)
#' r.m <- RSA(postVA~ePow*iPow, motcon)
#'
#' # Get boostrapped CIs with 10 bootstrap samples (usually this should be set to 5000 or higher),
#' # only from the SSD model
#' c1 <- confint(r.m, model="SSD", method="boot", R=10)
#' 
#' # Plot the final model
#' plot(r.m, model="RR", xlab="Explicit power motive", 
#' 		ylab="Implicit power motive", zlab="Affective valence")
#' }

RSA <- function(formula, data=NULL, center=FALSE, scale=FALSE, na.rm=FALSE, 
	out.rm=TRUE, breakline=FALSE, models="default", cubic=FALSE, 
	verbose=TRUE, add = "", 
	control.variables=c(), ...) {

	validmodels <- c("absdiff", "absunc", "diff", "additive", "IA", "sqdiff", "SRRR", "SRR", "RR", "SSD", "SRSD", "full", "null", "onlyx", "onlyy")
	if (length(models)==1 & models[1]=="all") {models <- validmodels}
	if (length(models)==1 & models[1]=="default") {models <- c("diff", "additive", "IA", "sqdiff", "SRRR", "SRR", "RR", "SSD", "SRSD", "full", "null", "onlyx", "onlyy")}
	if (any(!models %in% validmodels))
		stop("Unknown model name provided in parameter 'models'.")
	
	# set all result objects to NULL as default
	s.NULL <- s.full <- s.IA <- s.diff <- s.absdiff <- s.additive <- s.sqdiff <- s.SSD <- s.SRSD <- s.absunc <- s.cubic <- s.RR <- s.SRR <- s.SRRR <- s.onlyx <- s.onlyy <- NULL
	SRSD.rot <- ""
	SRRR.rot <- ""
	
	add <- paste0("\n# User defined syntax:\n", add)
	
	DV <- all.vars(formula)[1]
	IV1 <- all.vars(formula)[2]
	IV2 <- all.vars(formula)[3]

	## Step 0a: Standardize values (if requested) and calculate higher order terms
	df <- data
	df[, IV1] <- scale(df[, IV1], center=center, scale=scale)
	df[, IV2] <- scale(df[, IV2], center=center, scale=scale)
		
	df <- add.variables(formula, data.frame(data.matrix(df)))
	
	# give warnings if the zero point is outside of data range
	if (0 < min(df[, IV1], na.rm=TRUE) | 0 > max(df[, IV1], na.rm=TRUE)) 
		warning(paste("The numerical zero point is outside of the range of variable", IV1, ". Please consider re-centering the variable."))
	if (0 < min(df[, IV2], na.rm=TRUE) | 0 > max(df[, IV2], na.rm=TRUE)) 
		warning(paste("The numerical zero point is outside of the range of variable", IV2, ". Please consider re-centering the variable."))
		
	# give warning if one variable has a much higher range than the other variable
	if ((max(df[, IV1], na.rm=TRUE) - min(df[, IV1], na.rm=TRUE)) / (max(df[, IV2], na.rm=TRUE) - min(df[, IV2], na.rm=TRUE)) > 2)
		warning("Predictor variables have a very different range (by factor 2)- please check scaling of variables.")
	
	IV12 <- paste0(IV1, "2")
	IV22 <- paste0(IV2, "2")
	IV13 <- paste0(IV1, "3")
	IV23 <- paste0(IV2, "3")
	IV_IA <- paste0(IV1, "_", IV2)
	IV_IA2 <- paste0(IV1, "_", IV2, "2")
	IV_IA3 <- paste0(IV1, "2", "_", IV2)
	W_IV1 <- paste0("W_", IV1)
	W_IV2 <- paste0("W_", IV2)

	# define control variable
	CV <- ifelse(length(control.variables > 0), paste0(" + ", paste(control.variables, collapse=" + ")), "")

	## Run polynomial regression as a OLS linear model
	addcubic <- ""
	if (cubic==TRUE) addcubic <- paste0(" + ", paste(IV13, IV23, IV_IA2, IV_IA3, sep=" + "))
	f <- paste0(paste0(DV, " ~ ", paste(IV1, IV2, IV12, IV_IA, IV22, sep=" + ")), addcubic, CV)
	rs <- lm(f, df)
	
	# Mark outliers and influential cases according to Bollen & Jackman, 1980
	df0 <- df	# df0 keeps the original data frame, before outliers are removed (we need this in plotRSA, if we want to show the outliers)
	outs <- c()	# keeps the row numbers of the outliers
	inf <- influence.measures(rs)
	outs <- which(apply(inf$is.inf[, c("dffit", "cook.d", "hat")], 1, sum) == 3)
	names(outs) <- NULL
	if (out.rm == TRUE) {
		if (verbose==TRUE & length(outs)>0) {
			warning(paste("Removed", length(outs), "multivariate outlier(s) according to Bollen & Jackman (1980) criteria. Outliers are in row(s):", paste(outs , collapse=", ")))
			df <- df[-outs, ]
		}
	}

	## Test all models
	
# suppress some types of lavaan warning, which cannot be ruled out analytically ...
withCallingHandlers({	
	
	# Standard full polynomial of second degree
	poly <- paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b3*", IV12, " + b4*", IV_IA, " + b5*", IV22, CV)
	
	if ("null" %in% models) {
		s.NULL <- sem(paste0(DV, "~ 1 + 0*", IV1, " + 0*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, CV), data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}
	
	if ("additive" %in% models) {
		if (verbose==TRUE) print("Computing additive model (additive) ...")
		m.additive <-  paste(poly,
			"b3==0",
			"b4==0",
			"b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
		add, sep="\n")
		s.additive <- sem(m.additive, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}
	
	if ("onlyx" %in% models) {
		if (verbose==TRUE) print("Computing x + x^2 model (onlyx) ...")
		m.onlyx <-  paste(poly,
			"b2==0",
			"b4==0",
			"b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
		add, sep="\n")
		s.onlyx <- sem(m.onlyx, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}
	
	if ("onlyy" %in% models) {
		if (verbose==TRUE) print("Computing y + y^2 model (onlyy) ...")
		m.onlyy <-  paste(poly,
			"b1==0",
			"b3==0",
			"b4==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
		add, sep="\n")
		s.onlyy <- sem(m.onlyy, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}

	if ("diff" %in% models) {
		if (verbose==TRUE) print("Computing difference model (diff) ...")
		m.diff <- paste(poly,
			"b3==0",
			"b4==0",
			"b5==0",
			"b1 == -b2",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			add, sep="\n")
			s.diff <- sem(m.diff, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
		#summary(s.diff, fit.measures=TRUE)
	}

	if ("IA" %in% models) {
		if (verbose==TRUE) print("Computing interaction model (IA)...")
		m.IA <- paste(poly,
			"b3==0",
			"b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"as1X := b1 + p11*b2 + b4*p10 + 2*b5*p10*p11",
			"as2X := b3 + b4*p11 + (p11^2)*b5",
			"as1Y := b1/p11 + b2 - (2*b3*p10)/p11^2 - (b4*p10)/p11",
			"as2Y := b3/p11^2 + b4/p11 + b5",
			"as3X := b1 + p21*b2 + b4*p20 + 2*b5*p20*p21",
			"as4X := b3 + b4*p21 + (p21^2)*b5",
			"as3Y := b1/p21 + b2 - (2*b3*p20)/p21^2 - (b4*p20)/p21",
			"as4Y := b3/p21^2 + b4/p21 + b5",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
		add, sep="\n")
		
			s.IA <- sem(m.IA, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}
	
	if ("sqdiff" %in% models) {
		if (verbose==TRUE) print("Computing squared difference model (sqdiff) ...")
		m.sqdiff <- paste(poly,
			"b1==0",
			"b2==0",
			"b3==b5",
			"b3+b4+b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")			
			s.sqdiff <- sem(m.sqdiff, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}
	
	if ("SSD" %in% models) {
		if (verbose==TRUE) print("Computing shifted squared difference model (SSD) ...")
		m.SSD <- paste(poly,
			"b1==-b2",
			"b3==b5",
			"b3+b4+b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"C := b1 / (2*b3)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")			
			s.SSD <- sem(m.SSD, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}
	
	if (any(models %in% c("RR"))) {
		if (verbose==TRUE) print("Computing rising ridge model (RR) ...")
		m.RR <- paste(poly,
			"b1==b2",
			"b3==b5",
			"b3+b4+b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"meaneffect := b1+b2",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			
			add, sep="\n")
			s.RR <- sem(m.RR, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}
	
	if (any(models %in% c("SRR"))) {
		if (verbose==TRUE) print("Computing shifted rising ridge model (SRR) ...")
		m.SRR <- paste(poly,
			"b3==b5",
			"b3+b4+b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"meaneffect := b1+b2",
			"C := (b2-b1) / (4*b3)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			
			add, sep="\n")
			s.SRR <- sem(m.SRR, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}
	
	
	if (any(models %in% c("SRRR"))) {
			if (verbose==TRUE) print("Computing rotated and shifted rising ridge model (SRRR), up ...")
			m.SRRR.up <- paste(paste(poly, " + start(0.01)*", IV12, " + start(0.01)*", IV22),
				"b3 > 0.000001",
				"b5 > 0.000001",
				"b4^2 == 4*b3*b5",
				"a1 := b1+b2",
				"a2 := b3+b4+b5",
				"a3 := b1-b2",
				"a4 := b3-b4+b5",
				"meaneffect := (b2*b4 - 2*b1*b5) / b4",
				"C := (-2*b1*b5 - b2*b4) / (4*b4*b5)",
				"S := (-b4) / (2*b5)",
				# eigenvalues
				"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
				"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
		
				add, sep="\n")
				s.SRRR.up <- sem(m.SRRR.up, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)	
				
			if (verbose==TRUE) print("Computing rotated and shifted rising ridge model (SRRR), down ...")
			m.SRRR.down <- paste(paste(poly, " + start(-0.01)*", IV12, " + start(-0.01)*", IV22),
			#m.SRRR <- paste(paste(poly, " + start(-0.001)*", IV22),
				"b3 < 0.000001",
				"b5 < 0.000001",
				"b4^2 == 4*b3*b5",
				"a1 := b1+b2",
				"a2 := b3+b4+b5",
				"a3 := b1-b2",
				"a4 := b3-b4+b5",
				"meaneffect := (b2*b4 - 2*b1*b5) / b4",
				"C := (-2*b1*b5 - b2*b4) / (4*b4*b5)",
				"S := (-b4) / (2*b5)",
				# eigenvalues
				"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
				"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
	
				add, sep="\n")
				s.SRRR.down <- sem(m.SRRR.down, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)	
				
			if (inspect(s.SRRR.up, "converged") == FALSE & inspect(s.SRRR.down, "converged") == TRUE) {
				SRRR.rot <- "down"
			} else 
			if (inspect(s.SRRR.up, "converged") == TRUE & inspect(s.SRRR.down, "converged") == FALSE) {
				SRRR.rot <- "up"
			} else 
			if (inspect(s.SRRR.up, "converged") == TRUE & inspect(s.SRRR.down, "converged") == TRUE) {
				SRRR.rot <- ifelse(fitMeasures(s.SRRR.up, "chisq") > fitMeasures(s.SRRR.down, "chisq"), "down", "up")
			} else {
				if (verbose==TRUE) print("Warning: SRRR model has not converged (neither up nor down curvature)")
			}
			if (SRRR.rot == "up") {
				s.SRRR <- s.SRRR.up
			} else 
			if (SRRR.rot == "down") {
				s.SRRR <- s.SRRR.down
			}
			if (verbose == TRUE) print(paste0("Direction of SRRR curvature: ", SRRR.rot))
			
	}
	
	
	if (any(models %in% c("SRSD"))) {
		if (verbose==TRUE) print("Computing rotated squared difference model (SRSD), up ...")
		m.SRSD.up <- paste(paste(poly, " + start(0.001)*", IV22),
			"b1 == (b2*b4)/(2*b5)",
			"b3 > 0.000001",
			"b5 > 0.000001",
			"b4^2 == 4*b3*b5",
			
			"C := -.5*(b2/b5)",
			"S := (-b4) / (2*b5)",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"as1X := b1 + p11*b2 + b4*p10 + 2*b5*p10*p11",
			"as2X := b3 + b4*p11 + (p11^2)*b5",
			"as1Y := b1/p11 + b2 - (2*b3*p10)/p11^2 - (b4*p10)/p11",
			"as2Y := b3/p11^2 + b4/p11 + b5",
			"as3X := b1 + p21*b2 + b4*p20 + 2*b5*p20*p21",
			"as4X := b3 + b4*p21 + (p21^2)*b5",
			"as3Y := b1/p21 + b2 - (2*b3*p20)/p21^2 - (b4*p20)/p21",
			"as4Y := b3/p21^2 + b4/p21 + b5",
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")
			s.SRSD.up <- sem(m.SRSD.up, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
			
			
		if (verbose==TRUE) print("Computing rotated squared difference model (SRSD), down ...")
		m.SRSD.down <- paste(paste(poly, " + start(-0.001)*", IV22),
			"b1 == (b2*b4)/(2*b5)",
			"b3 < 0.000001",
			"b5 < 0.000001",
			"b4^2 == 4*b3*b5",
		
			"C := -.5*(b2/b5)",
			"S := (-b4) / (2*b5)",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"as1X := b1 + p11*b2 + b4*p10 + 2*b5*p10*p11",
			"as2X := b3 + b4*p11 + (p11^2)*b5",
			"as1Y := b1/p11 + b2 - (2*b3*p10)/p11^2 - (b4*p10)/p11",
			"as2Y := b3/p11^2 + b4/p11 + b5",
			"as3X := b1 + p21*b2 + b4*p20 + 2*b5*p20*p21",
			"as4X := b3 + b4*p21 + (p21^2)*b5",
			"as3Y := b1/p21 + b2 - (2*b3*p20)/p21^2 - (b4*p20)/p21",
			"as4Y := b3/p21^2 + b4/p21 + b5",
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")
			s.SRSD.down <- sem(m.SRSD.down, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	
			if (inspect(s.SRSD.up, "converged") == FALSE & inspect(s.SRSD.down, "converged") == TRUE) {
				SRSD.rot <- "down"
			} else 
			if (inspect(s.SRSD.up, "converged") == TRUE & inspect(s.SRSD.down, "converged") == FALSE) {
				SRSD.rot <- "up"
			} else 
			if (inspect(s.SRSD.up, "converged") == TRUE & inspect(s.SRSD.down, "converged") == TRUE) {
				SRSD.rot <- ifelse(fitMeasures(s.SRSD.up, "chisq") > fitMeasures(s.SRSD.down, "chisq"), "down", "up")
			} else {
				if (verbose==TRUE) print("Warning: SRSD model has not converged (neither up nor down curvature)")
			}
			if (SRSD.rot == "up") {
				s.SRSD <- s.SRSD.up
			} else 
			if (SRSD.rot == "down") {
				s.SRSD <- s.SRSD.down
			}
			if (verbose == TRUE) print(paste0("Direction of SRSD curvature: ", SRSD.rot))

		
	}
	
	
	if ("full" %in% models) {
		if (verbose==TRUE) print("Computing polynomial model ...")
		m.full <-  paste(poly,
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"as1X := b1 + p11*b2 + b4*p10 + 2*b5*p10*p11",
			"as2X := b3 + b4*p11 + (p11^2)*b5",
			"as1Y := b1/p11 + b2 - (2*b3*p10)/p11^2 - (b4*p10)/p11",
			"as2Y := b3/p11^2 + b4/p11 + b5",
			"as3X := b1 + p21*b2 + b4*p20 + 2*b5*p20*p21",
			"as4X := b3 + b4*p21 + (p21^2)*b5",
			"as3Y := b1/p21 + b2 - (2*b3*p20)/p21^2 - (b4*p20)/p21",
			"as4Y := b3/p21^2 + b4/p21 + b5",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add,
			sep="\n"
		)
		s.full <- sem(m.full, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}
	
	
	if (cubic==TRUE) {
		if (verbose==TRUE) print("Computing full cubic model ...")
		m.cubic <-  paste(paste0(poly, " + b9*", IV13, " + b10*", IV_IA2, " + b11*", IV_IA3, " + b12*", IV23),
			"u1 := b1 + b2",				# linear part of LOC
			"u2 := b3 + b4 + b5",			# quadratic part of LOC
			"u3 := b9 + b10 + b11 + b12",	# cubic part of LOC
			"v1 := b1 - b2",				# linear part of LOIC
			"v2 := b3 - b4 + b5",			# quadratic part of LOIC
			"v3 := b9 + b10 - b11 - b12",	# cubic part of LOIC: If v3 != 0, then there is an enhancement effect (i.e., the slope is different on both sides of the optimum)
			add,
			sep="\n"
		)
		#print(m.cubic)
		s.cubic <- sem(m.cubic, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)		
	}
	
	#m.absdiff.JRE <-  paste(
	#	paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, " + 0*W.JRE + b7*W.JRE_", IV1, " + b8*W.JRE_", IV2),
	#	"b1 == -b2",
	#	"b7 == -b8",
	#	"b7 == -2*b1",
	#	add, sep="\n")
	#s.absdiff.JRE <-  sem(m.absdiff.JRE, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	#summary(s.absdiff.JRE, fit.measures=TRUE)
	
	# the unconstrained absolute difference model - Edwards (2002) formula
	#m.absunc.JRE <-  paste(
	#	paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, " + b6*W.JRE + b7*W.JRE_", IV1, " + b8*W.JRE_", IV2),
	#	add, sep="\n")
	#s.absunc.JRE <-  sem(m.absunc.JRE, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	#summary(s.absunc.JRE, fit.measures=TRUE)
	
	
	if ("absdiff" %in% models) {
		if (verbose==TRUE) print("Computing constrained absolute difference model ...")
		m.absdiff <-  paste(
			paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b6*W + b7*W_", IV1, " + b8*W_", IV2),
			"b1 == 0",
			"b2 == 0",
			"b6 == 0",
			"b7 == -b8",
			add, sep="\n")
			
			s.absdiff <- sem(m.absdiff, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}
	
	if ("absunc" %in% models) {
		# the unconstrained absolute difference model - new formula
		if (verbose==TRUE) print("Computing unconstrained absolute difference model ...")
		m.absunc <-  paste(
			paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b6*W + b7*W_", IV1, " + b8*W_", IV2),
			ifelse(breakline==FALSE, "b6==0", ""),
			add, sep="\n")
			
			s.absunc <- sem(m.absunc, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	}
	
},	  # end of "withCallingHandlers"

# suppress two types of warning
  warning=function(w) {
	   W <- as.character(w$call)
	   if (
		   (W[1] == "sqrt" & W[2] == "diag(def.cov)" & grepl("NaNs", w$message)) |
		   (W[1] == "sqrt" & W[2] == "b3 * b5") |
		   (W[1] == "nlminb" & W[2] == "x.par")
		  ) {invokeRestart("muffleWarning")}
} )


	# ---------------------------------------------------------------------
	# Sanity check: Check results for convergence problems
	# Sometimes the SRRR or the SRSD model find a bad solution and have a higher chi2 than their nested models (which is theoretically not possible)
	chisq1 <- plyr::ldply(list(full=s.full, SRRR=s.SRRR, SRR=s.SRR, RR=s.RR, sqdiff=s.sqdiff), function(x) {
		chi <- -1
		if (!is.null(x)) {
			if (inspect(x, "converged")==TRUE) chi <-  fitMeasures(x, "chisq")
		}
		return(chi)
	})
	chisq1 <- chisq1[chisq1[, 2]>=0, ]
	if (nrow(chisq1)>1) {
		chisq1$lag <- c(diff(chisq1[, 2], lag=1), NA)
		if (any(chisq1$lag < 0, na.rm=TRUE)) {
			warning(paste0("There are convergence problems with model ", chisq1[which(chisq1$lag < 0), ".id"], ". Its chi-square value is higher than that of a nested model, which is theoretically not possible. Please inspect the results with care, using the compare()-function"))
		}
	}
	
	chisq2 <- plyr::ldply(list(full=s.full, SRRR=s.SRRR, SRSD=s.SRSD, SSD=s.SSD, sqdiff=s.sqdiff), function(x) {
		chi <- -1
		if (!is.null(x)) {
			if (inspect(x, "converged")==TRUE) chi <-  fitMeasures(x, "chisq")
		}
		return(chi)
		
	})
	chisq2 <- chisq2[chisq2[, 2]>=0, ]
	
	if (nrow(chisq1)>1) {
		chisq2$lag <- c(diff(chisq2[, 2], lag=1), NA)
		if (any(chisq2$lag < 0, na.rm=TRUE)) {
			warning(paste0("There are convergence problems with model ", chisq2[which(chisq2$lag < 0), ".id"], ". Its chi-square value is higher than that of a nested model, which is theoretically not possible. Please inspect the results with care, using the compare()-function"))
		}
	}

	
	# ---------------------------------------------------------------------
	# Build results object
	res <- list(
		models = list(null=s.NULL, full=s.full, IA=s.IA, diff=s.diff, absdiff=s.absdiff, additive=s.additive, sqdiff=s.sqdiff, SRRR=s.SRRR, SRR=s.SRR, RR=s.RR, SSD=s.SSD, SRSD=s.SRSD, absunc=s.absunc, cubic=s.cubic, onlyx=s.onlyx, onlyy=s.onlyy), 
		SRSD.rot = SRSD.rot, SRRR.rot = SRRR.rot, LM=rs, formula=formula, 
		data=df, data.original=df0, out.rm = out.rm, outliers = outs, DV=DV, IV1=IV1, IV2=IV2, IV12=IV12, IV22=IV22, IV_IA=IV_IA, W_IV1=W_IV1, W_IV2=W_IV2, IV13=IV13, IV23=IV23, IV_IA2=IV_IA2, IV_IA3=IV_IA3, 
		r.squared = summary(rs)$r.squared)
	
	attr(res, "class") <- "RSA"
	return(res)
}



getfit <- function(x) {
	res <- plyr::laply(x$models, function(m) {fitMeasures(m)[c("aic", "bic", "cfi", "tli", "chisq", "df", "pvalue", "rmsea", "srmr")]})
	rownames(res) <- names(x$models)
	res <- data.frame(res)
	res$rel.aic <- (res$aic - min(res$aic))/(max(res$aic) - min(res$aic))
	res$rel.bic <- (res$bic - min(res$bic))/(max(res$bic) - min(res$bic))
	return(res)
}


bestmodel <- function(x) {
	
	F <- summary(x$LM)$fstatistic
	p.model <- 1-pf(F[1], F[2], F[3])
	if (p.model > .05) {
		return("Overall model is not significant.")
	}
	
	f <- getfit(x)
	
	# Choose family based on AIC and BIC
	
	# if only one model is selected by AIC: that's it!
	if (length(which(f$rel.aic < .01)) == 1) {
		m <- rownames(f)[which(f$rel.aic < .01)]
	} else 	
	# AIC is ambiguous: let BIC decide!
	if (length(which(f$rel.aic < .01)) > 1 & length(which(f$rel.bic < .01)) == 1) {
		m <- rownames(f)[which(f$rel.bic < .01)]
	} else 	
	# If AIC and BIC are ambiguous: let AIC decide
	if (length(which(f$rel.aic < .01)) > 1 & length(which(f$rel.bic < .01)) > 1) {
		m <- rownames(f)[which.min(f$rel.aic)]
	} else {
		warning("Could not determine best model! err1")
		return("err1")
	}
	
	
	## Let ALWAYS BIC decide
	#m <- rownames(f)[which.min(f$rel.bic)]
	
	if (m %in% c("full", "IA", "additive", "diff")) {
		# find the last model that is not significantly different from full model
		a1 <- with(x$models, {anova(diff, additive, IA, full)})
		M <- which(a1[, "Pr(>Chisq)"] < .05)
		if (length(M) > 0) {
			return(rownames(a1)[min(which(a1[, "Pr(>Chisq)"] < .05))-1])
			} else {
				return("diff")
			}
	} else 
	if (m %in% c("sqdiff")) {
		return("sqdiff")
	} else 
	if (m %in% c("absdiff", "absunc")) {
		a1 <- with(x$models, {anova(absdiff, absunc)})
		return(ifelse(a1[2, "Pr(>Chisq)"] < .05, "absunc", "absdiff"))
	} else {
		warning("Could not determine best model! err2")
		return("err2")
	}
}

