#' @title Performs several RSA model tests on a data set with two predictors
#'
#' @description
#' Performs several RSA model tests on a data set with two predictors
#'
#' @details
#' Even if the main variables of the model are normally distributed, their squared terms and interaction terms are necessarily non-normal. By default, the RSA function uses a scaled test statistic (\code{test="Satorra-Bentler"}) and robust standard errors (\code{se="robust"}), which are robust against violations of the normality assumption. 
#'
#' \emph{Why does my standard polynomial regression give different p-values and SEs than the RSA package? Shouldn't they be the same?} This is due to the robust standard errors employed in the RSA package. If you set \code{estimator="ML"} and \code{se="standard"}, you get p-values that are very close to the standard approach. (They might still not be identical because the standard regression approach usually uses an OLS estimator and RSA uses an ML estimator).
#'
#' Experimental feature (use with caution!): You can also fit \strong{binary outcome variables} with a probit link function. For that purpose, the response variable has to be defined as "ordered", and the \code{lavaan} estimator changed to "WLSMV": \code{r1.binary <- RSA(z.binary~x*y, df, ordered="z.binary", estimator="WLSMV", model="full")} (for more details see the help file of the \code{\link[lavaan:sem]{sem}} function in the \code{lavaan} package.). The results can also be plotted with probabilities on the z axis using the probit link function: \code{plot(r1.binary, link="probit", zlim=c(0, 1), zlab="Probability")}. For plotting, the binary outcome variable must be coded with 0 and 1 (not as a factor).
#' \code{lavaan} at the moment only supports a probit link function for binary outcomes, not a logit link. Please be aware that this experimental feature can fit the full model, but most other functions (such as model comparisons) might break and errors might show up.
#'
#' @export
#' @import lavaan
#' @import ggplot2
#' @import lattice
#' @import RColorBrewer
#' @param formula A formula in the form \code{z ~ x*y}, specifying the variable names used from the data frame, where z is the name of the response variable, and x and y are the names of the predictor variables.
#' @param data A data frame with the variables
#' @param center Method for centering the predictor variables before the analysis. Default option ("none") applies no centering. "pooled" centers the predictor variables on their \emph{pooled} sample mean, which preserves the commensurability of the predictor scales. "variablewise" centers the predictor variables on \emph{their respective} sample mean. You should think carefully before applying the "variablewise" option, as centering the predictor variables at different values (e.g., their respective means) can affect the commensurability of the predictor scales.
#' @param scale Method for scaling the predictor variables before the analysis. Default option ("none") applies no scaling. "pooled" scales the predictor variables on their \emph{pooled} sample SD, which preserves the commensurability of the predictor scales. "variablewise" scales the predictor variables on \emph{their respective} sample SD. You should think carefully before applying the "variablewise" option, as scaling the predictor variables at different values (e.g., their respective SDs) can affect the commensurability of the predictor scales.
#' @param na.rm Remove missings before proceeding?
#' @param add Additional syntax that is added to the lavaan model. Can contain, for example, additional constraints, like "p01 == 0; p11 == 0"
#' @param out.rm Should outliers according to Bollen & Jackman (1980) criteria be excluded from the analyses? In large data sets this analysis is the speed bottleneck. If you are sure that no outliers exist, set this option to FALSE for speed improvements.
#' @param breakline Should the breakline in the unconstrained absolute difference model be allowed (the breakline is possible from the model formulation, but empirically rather unrealistic ...). Defaults to \code{FALSE}
#' @param verbose Should additional information during the computation process be printed?
#' @param models A vector with names of all models that should be computed. Should be any from \code{c("absdiff", "absunc", "diff", "mean", "additive", "IA", "SQD", "RR", "SRR", "SRRR", "SSQD", "SRSQD", "full", "null", "onlyx", "onlyy", "onlyx2", "onlyy2", "cubic", "CA", "RRCA", "CL", "RRCL")}. For \code{models="all"}, all models are computed, for \code{models="default"} all models besides absolute difference models are computed.
#' @param cubic Should the cubic models with the additional terms Y^3, XY^2, YX^2, and X^3 be included?
#' @param control.variables A string vector with variable names from \code{data}. These variables are added as linear predictors to the model (in order "to control for them"). No interactions with the other variables are modeled.
#' @param center.control.variables Should the control variables be centered before analyses? This can improve interpretability of the intercept, which will then reflect the predicted outcome value at the point (X,Y)=(0,0) when all control variables take their respective \emph{average} values.
#' @param estimator Type of estimator that should be used by lavaan. Defaults to "MLR", which provides robust standard errors, a robust scaled test statistic, and can handle missing values. If you want to reproduce standard OLS estimates, use \code{estimator="ML"} and \code{se="standard"}
#' @param se Type of standard errors. This parameter gets passed through to the \code{sem} function of the \code{lavaan} package. See options there. By default, robust SEs are computed. If you use \code{se="boot"}, \code{lavaan} provides CIs and p-values based on the bootstrapped standard error. If you use \code{confint(..., method="boot")}, in contrast, you get CIs and p-values based on percentile bootstrap (see also \code{\link{confint.RSA}}).
#' @param missing Handling of missing values (this parameter is passed to the \code{lavaan} \code{sem} function). By default (\code{missing=NA}), Full Information Maximum Likelihood (FIML) is employed in case of missing values. If cases with missing values should be excluded, use \code{missing = "listwise"}.
#' @param ... Additional parameters passed to the \code{lavaan} \code{\link[lavaan:sem]{sem}} function.
#'
#'
#' @note 
#' For explanations of the meaning of the various different models that can be estimated, please see Schönbrodt (2016) for the second-order models (i.e., all models but "CA", "RRCA", "CL", "RRCL") and Humberg et al. (in press) for the third-order (cubic) models ("CA", "RRCA", "CL", "RRCL"). 
#' 
#' For most of the second-order models, several auxiliary parameters are computed from the estimated model coefficients (e.g., a1, ..., a5, p10, p11, p20, p21) and printed in the \code{summary} output. They can be used to guide interpretation by means of response surface methodology. Some references that explain how to use these parameters for interpretation are Edwards (2002; comprehensive overview of response surface methodology), Humberg et al. (2019; interpretation of a1, a2, a3, a4, p10, and p11, and how to use them to investigate congruence effects), Nestler et al. (2019; interpretation of a1, a2, a3, a4, and a5, and how to use them to investigate congruence effects, see in particular Appendix A for the introduction of a5), and Schönbrodt et al. (2018; interpretation of a1, ..., a5, see in particular Appendix A for a5).
#' 
#' The print function provides descriptive statistics about discrepancies in the predictors (with respect to numerical congruence). A cutpoint of |delta z| > 0.5 is used. The computation generally follows the idea of Shannock et al (2010) and Fleenor et al. (1996). However, in contrast to them, we standardize to the common mean and the common SD of both predictor variables. Otherwise we would break commensurability, and a person who has x=y in the unstandardized variable could become incongruent after variable-wise standardization. See also our discussion of commensurability and scale transformation in the cubic RSA paper (Humberg et al., in press; see pp. 35 - 37 in the preprint at \url{https://osf.io/preprints/psyarxiv/v6m35)}.
#' 
#' @references 
#' Edwards, J. R. (2002). Alternatives to difference scores: Polynomial regression analysis and response surface methodology. In F. Drasgow & N. W. Schmitt (Eds.), \emph{Advances in measurement and data analysis} (pp. 350–400). San Francisco, CA: Jossey-Bass.
#' 
#' Humberg, S., Nestler, S., & Back, M. D. (2019). Response Surface Analysis in Personality and Social Psychology: Checklist and Clarifications for the Case of Congruence Hypotheses. \emph{Social Psychological and Personality Science}, 10(3), 409–419. doi:10.1177/1948550618757600
#' 
#' Humberg, S., Schönbrodt, F. D., Back, M. D., & Nestler, S. (in press). Cubic response surface analysis: Investigating asymmetric and level-dependent congruence effects with third-order polynomial models. Psychological Methods. doi:10.1037/met0000352
#' 
#' Nestler, S., Humberg, S., & Schönbrodt, F. D. (2019). Response surface analysis with multilevel data: Illustration for the case of congruence hypotheses. \emph{Psychological Methods}, 24(3), 291–308. doi:10.1037/met0000199
#' 
#' Schönbrodt, F. D. (2016). \emph{Testing fit patterns with polynomial regression models.} Retrieved from osf.io/3889z
#' 
#' Schönbrodt, F. D., Humberg, S., & Nestler, S. (2018). Testing similarity effects with dyadic response surface analysis. \emph{European Journal of Personality}, 32(6), 627-641. doi:10.1002/per.2169
#' 
#'
#' @seealso \code{\link{demoRSA}}, \code{\link{plotRSA}}, \code{\link{RSA.ST}}, \code{\link{confint.RSA}}, \code{\link{compare}}
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
#' 	SD <- (x-y)^2
#' 	z.diff <- diff + rnorm(n, 0, err)
#' 	z.abs <- absdiff + rnorm(n, 0, err)
#' 	z.sq <- SD + rnorm(n, 0, err)
#' 	z.add <- diff + 0.4*x + rnorm(n, 0, err)
#' 	z.complex <- 0.4*x + - 0.2*x*y + + 0.1*x^2 - 0.03*y^2 + rnorm(n, 0, err)
#' })
#' df$z.binary <- as.numeric(df$z.sq < median(df$z.sq))
#'
#' \dontrun{
#' r1 <- RSA(z.sq~x*y, df)
#' summary(r1)
#' compare(r1)
#' plot(r1)
#' plot(r1, model="SRSQD")
#' plot(r1, model="full", type="c")
#' getPar(r1, "coef")	# print model parameters including SE and CI
#' RSA.ST(r1)	# get surface parameters
#'
#' # Example with binary outcome 
#' (probit regression, see Details; Experimental and a dirty workaround!).
#' # The standard summary output does not work; you have to access the full model directly:
#' r1.binary <- RSA(z.binary~x*y, df, ordered="z.binary", estimator="WLSMV", 
#'     model="full", se="standard")
#' # --> ignore the warning
#' summary(r1.binary$models[["full"]])
#' plot(r1.binary, link="probit", zlim=c(0, 1), zlab="Probability")
#' 
#' # Motive congruency example
#' data(motcon)
#' r.m <- RSA(postVA~ePow*iPow, motcon)
#'
#' # Get boostrapped CIs with 10 bootstrap samples (usually this should be set to 5000 or higher),
#' # only from the SSQD model
#' c1 <- confint(r.m, model="SSQD", method="boot", R=10)
#' 
#' # Plot the final model
#' plot(r.m, model="RR", xlab="Explicit power motive", 
#' 		ylab="Implicit power motive", zlab="Affective valence")
#' 		
#' # Inclusion of control variables: Fake data on self-other agreement
#' data(selfother)
#' r.c <- RSA(liking~IQ_self*IQ_friend, 
#'            center="pooled",
#'            control.variables=c("age", "int"),
#'            center.control.variables = TRUE,
#'            data=selfother)
#' summary(r.c)	
#' 		
#' }

RSA <- function(formula, data=NULL, center="none", scale="none", na.rm=FALSE, 
	out.rm=TRUE, breakline=FALSE, models="default", cubic=FALSE, 
	verbose=TRUE, add = "", estimator="MLR",
	se = "robust", missing=NA, control.variables=c(), center.control.variables=FALSE, ...) {

	validmodels <- c("absdiff", "absunc", "diff", "mean", "additive", "IA", "SQD", "SRRR", "SRR", "RR", "SSQD", "SRSQD", "full", "null", "onlyx", "onlyy", "onlyx2", "onlyy2", "weak", "strong", "cubic", "CA", "CL", "RRCA", "RRCL")
	
	if (length(models)==1 & models[1]=="all") {models <- validmodels}
	#if (length(models)==1 & models[1]=="default") {models <- c("diff", "mean", "additive", "IA", "SQD", "SRRR", "SRR", "RR", "SSQD", "SRSQD", "full", "null", "onlyx2", "onlyy2", "onlyx", "onlyy", "weak", "strong")}
	if (length(models)==1 & models[1]=="default" & cubic==FALSE) {models <- c("additive", "IA", "SQD", "SRRR", "SRR", "RR", "SSQD", "SRSQD", "full", "null", "onlyx2", "onlyy2", "onlyx", "onlyy")} 
	if (length(models)==1 & models[1]=="default" & cubic==TRUE) {models <- c("additive", "IA", "SQD", "SRRR", "SRR", "RR", "SSQD", "SRSQD", "full", "null", "onlyx2", "onlyy2", "onlyx", "onlyy", "cubic", "CA", "CL", "RRCA", "RRCL")}
	if (any(!models %in% validmodels))
		stop("Unknown model name provided in parameter 'models'.")
	
	# set cubic flag if any third-order model is contained in the models vector
	cubicmodels <- c("cubic", "CA", "CL", "RRCA", "RRCL")
	is.cubic <- any(models %in% cubicmodels)
	
	# set all result objects to NULL as default
	s.NULL <- s.full <- s.IA <- s.diff <- s.mean <- s.absdiff <- s.additive <- s.SQD <- s.SSQD <- s.SRSQD <- s.absunc <- s.cubic <- s.RR <- s.SRR <- s.SRRR <- s.onlyx <- s.onlyy <- s.onlyx2 <- s.onlyy2 <- s.weak <- s.strong <- s.CA <- s.CL <- s.RRCA <- s.RRCL <- NULL
	SRSQD.rot <- ""
	SRRR.rot <- ""
	
	add <- paste0("\n# User defined syntax:\n", add)
	
	DV <- all.vars(formula)[1]
	IV1 <- all.vars(formula)[2]
	IV2 <- all.vars(formula)[3]

	# reduce data frame to actually used variables
	df <- data[, c(DV, IV1, IV2, control.variables)]	
	
	
	## Step 0a: Standardize predictor variables (if requested)

	# when boolean options are used for center/scale (as was necessary in earlier package versions), restore behavior of earlier versions
	if(center==TRUE){center <- "variablewise"}
	if(scale==TRUE){scale <- "variablewise"}
	if (!center %in% c("none", "variablewise", "pooled")) {
		warning('Parameter `center` must be one of c("none", "variablewise", "pooled"), or TRUE (which corresponds to "variablewise") or FALSE (which corresponds to "none"). Setting to "none".')
		center <- "none"
	}
	
	# "variablewise": standardization at respective variable means/SDs
	
	if(center=="variablewise"){
	  warning("You specified 'variablewise' centering of the predictor variables. Make sure to check whether this is a good choice, as it can distort commensurability. Use center='pooled' if you want to center at the pooled mean instead.", call.=FALSE)
	  df[, IV1] <- scale(df[, IV1], center=TRUE, scale=FALSE)
	  df[, IV2] <- scale(df[, IV2], center=TRUE, scale=FALSE)
	}
	
	if(scale=="variablewise"){
	  warning("You specified 'variablewise' scaling of the predictor variables. Make sure to check whether this is a good choice, as it can distort commensurability. Use scale='pooled' if you want to scale at the pooled SD instead.", call.=FALSE)
	  df[, IV1] <- scale(df[, IV1], center=FALSE, scale=TRUE)
	  df[, IV2] <- scale(df[, IV2], center=FALSE, scale=TRUE)
	}
	
	# "pooled": standardization at pooled means/SDs
	
	if(center=="pooled"){
	  pooled.mean <-  mean( c(df[,IV1],df[,IV2]) , na.rm=T)
	  df[, IV1] <- df[, IV1] - pooled.mean
	  df[, IV2] <- df[, IV2] - pooled.mean
	}
	
	if(scale=="pooled"){
	  pooled.sd <-  sd( c(df[,IV1],df[,IV2]) , na.rm=T)
	  df[, IV1] <- df[, IV1] / pooled.sd
	  df[, IV2] <- df[, IV2] / pooled.sd
	}
	
	
	# calculate higher order terms
	df <- add.variables(formula, data.frame(data.matrix(df)))
	
	# give warnings if the zero point is outside of data range
	if (0 < min(df[, IV1], na.rm=TRUE) | 0 > max(df[, IV1], na.rm=TRUE)) 
		warning(paste("The numerical zero point is outside of the range of variable", IV1, ". Please consider re-centering the variable."))
	if (0 < min(df[, IV2], na.rm=TRUE) | 0 > max(df[, IV2], na.rm=TRUE)) 
		warning(paste("The numerical zero point is outside of the range of variable", IV2, ". Please consider re-centering the variable."))
		
	# give warning if one variable has a much higher range than the other variable
	if ((max(df[, IV1], na.rm=TRUE) - min(df[, IV1], na.rm=TRUE)) / (max(df[, IV2], na.rm=TRUE) - min(df[, IV2], na.rm=TRUE)) > 2)
		warning("Predictor variables have a very different range (by factor 2 or larger) - please check scaling of variables.")
	
	
	
	# set defaults for missing option
	if (is.na(missing)) {
		if (any(is.na(df))) {
			missing <- "fiml"
			warning("There are missing values in your data set. Model is computed with option `missing = 'fiml'`. This is only valid if the data are missing completely at random (MCAR) or missing at random (MAR)! If you want to exclude NAs, use `missing = 'listwise'`", call.=FALSE)
		} else {
			missing <- "listwise"
		}
	}
	
	
	# temporary workaround to avoid listwise deletion when the data contains missings
	#
	# lavaan versions < 0.6.3 will apply listwise deletion (because fixed.x = TRUE in the sem() models) and have no workaround implemented, so they should not be applied.
	if (any(is.na(df)) & missing=="fiml" & packageVersion("lavaan") < "0.6.3") {stop("Please install the latest version of 'lavaan'! lavaan versions < 0.6.3 will apply listwise deletion (because fixed.x = TRUE in the sem() models) and have no workaround implemented, so they should not be applied.")}
	#
	# lavaan versions >= 0.6.3 have the fiml.x workaround
	if (any(is.na(df)) & missing=="fiml" & packageVersion("lavaan") >= "0.6.3") {missing <- "fiml.x"}
	
	
	IV12 <- paste0(IV1, "2")
	IV22 <- paste0(IV2, "2")
	IV13 <- paste0(IV1, "3")
	IV23 <- paste0(IV2, "3")
	IV_IA <- paste0(IV1, "_", IV2)
	IV_IA2 <- paste0(IV1, "2", "_", IV2)
	IV_IA3 <- paste0(IV1, "_", IV2, "2")
	W_IV1 <- paste0("W_", IV1)
	W_IV2 <- paste0("W_", IV2)

	
	## control variables
	
	# control variables: set cv flag if control variables are involved
	is.cv <- length(control.variables) > 0
	
	# center control variables (if requested)
	if(is.cv & center.control.variables){
	  df[,control.variables] <- scale( df[,control.variables], center=TRUE, scale=FALSE )
	}
	
	# define control variable part of the regression model
	CV <- ifelse(is.cv, paste0(" + ", paste(control.variables, collapse=" + ")), "")

	
	## Run polynomial regression as a OLS linear model
	addcubic <- ""
	if (is.cubic) addcubic <- paste0(" + ", paste(IV13, IV_IA2, IV_IA3, IV23, sep=" + "))
	f <- paste0(paste0(DV, " ~ ", paste(IV1, IV2, IV12, IV_IA, IV22, sep=" + ")), addcubic, CV)
	lm.full <- lm(f, df, na.action=na.exclude)
	
	# ---------------------------------------------------------------------
	#  Mark outliers and influential cases 
	
	# define the defaults
	if (is.null(out.rm) || (typeof(out.rm) == "logical" && out.rm == TRUE)) {
		out.rm <- "bj1980"
	}
	if ((typeof(out.rm) == "logical" && out.rm == FALSE)) {
		out.rm <- "none"
	}
	out.rm <- match.arg(out.rm, c("bj1980", "robust", "none"))
	df$out <- FALSE

	#...according to Bollen & Jackman, 1980
	if (out.rm == "bj1980") {
		inf <- influence.measures(lm.full)
		df$out <- apply(inf$is.inf[, c("dffit", "cook.d", "hat")], 1, sum) == 3
		n.out <- sum(na.omit(df$out) == TRUE)
		if (verbose==TRUE & n.out>0) {
			warning(paste("Removed", n.out, "multivariate outlier(s) according to Bollen & Jackman (1980) criteria. Outliers are in row(s):", paste(which(df$out == TRUE) , collapse=", ")))
		}
	}
	#...outpro function from WRS
	if (out.rm == "robust") {
		stop("Robust outlier detection not implemented yet.")
		# out.pro <- outpro(df[, c(DV, IV1, IV2)])
# 		n.out <- out.pro$n.out
# 		df$out[out.pro$out.id] <- TRUE
# 		if (verbose==TRUE & n.out>0) {
# 			warning(paste("Removed", n.out, "multivariate outlier(s) using a robust procedure (Wilcox, 2012). Outliers are in row(s):", paste(which(df$out == TRUE) , collapse=", ")))
# 		}
	}
	
	# Rows with outliers have a NA in $out. Should they be kept or removed? I will keep them, otherwise it would not make sense to have fiml estimation, as all NAs are automatically processed as outliers.
	# FIXME: is this a good default choice?
	df$out[is.na(df$out)] <- FALSE


	## Test all models of second degree
	
# suppress some types of lavaan warning, which cannot be ruled out analytically ...
withCallingHandlers({	
	
	# Standard full polynomial of second degree
	poly <- paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b3*", IV12, " + b4*", IV_IA, " + b5*", IV22, CV)
	
	# Standard full polynomial of third degree
	polycubic <- paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b3*", IV12, " + b4*", IV_IA, " + b5*", IV22, " + b6*", IV13, " + b7*", IV_IA2, " + b8*", IV_IA3, " + b9*", IV23, CV) 
	
	
	# always estimate the null model (intercept + control variables only)
		m.null <- ifelse(is.cubic, 
		                 paste0(DV, "~ 1 + 0*", IV1, " + 0*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, " + 0*", IV13, " + 0*", IV_IA2, " + 0*", IV_IA3, " + 0*", IV23, CV),
		                 paste0(DV, "~ 1 + 0*", IV1, " + 0*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, CV))
	  s.NULL <- sem(m.null, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)

	  
	if ("additive" %in% models) {
	  if (verbose==TRUE) print("Computing additive model (additive) ...")
	  m.additive <-  paste(ifelse(is.cubic, polycubic, poly),
	                       "b3==0",
	                       "b4==0",
	                       "b5==0",
	                       ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
	                       "a1 := b1+b2",
	                       "a2 := b3+b4+b5",
	                       "a3 := b1-b2",
	                       "a4 := b3-b4+b5",
	                       "a5 := b3-b5",
	                       add, sep="\n")
	  s.additive <- sem(m.additive, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("onlyx2" %in% models) {
		if (verbose==TRUE) print("Computing x + x^2 model (onlyx2) ...")
		m.onlyx2 <-  paste(ifelse(is.cubic,polycubic,poly),
			"b2==0",
			"b4==0",
			"b5==0",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
		add, sep="\n")
		s.onlyx2 <- sem(m.onlyx2, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("onlyy2" %in% models) {
		if (verbose==TRUE) print("Computing y + y^2 model (onlyy2) ...")
		m.onlyy2 <-  paste(ifelse(is.cubic,polycubic,poly),
			"b1==0",
			"b3==0",
			"b4==0",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
		add, sep="\n")
		s.onlyy2 <- sem(m.onlyy2, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}

	if ("onlyx" %in% models) {
		if (verbose==TRUE) print("Computing x model (onlyx) ...")
		m.onlyx <-  paste(ifelse(is.cubic,polycubic,poly),
			"b2==0",
			"b3==0",
			"b4==0",
			"b5==0",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
		add, sep="\n")
		s.onlyx <- sem(m.onlyx, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("onlyy" %in% models) {
		if (verbose==TRUE) print("Computing y model (onlyy) ...")
		m.onlyy <-  paste(ifelse(is.cubic,polycubic,poly),
			"b1==0",
			"b3==0",
			"b4==0",
			"b5==0",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
		add, sep="\n")
		s.onlyy <- sem(m.onlyy, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	

	if ("diff" %in% models) {
		if (verbose==TRUE) print("Computing difference model (diff) ...")
		m.diff <- paste(ifelse(is.cubic,polycubic,poly),
			"b3==0",
			"b4==0",
			"b5==0",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"b1 == -b2",
			"a1 := b1+b2",
			"a2 := 0",
			"a3 := b1-b2",
			"a4 := 0",
			add, sep="\n")
			s.diff <- sem(m.diff, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
		#summary(s.diff, fit.measures=TRUE)
	}
	
	if ("mean" %in% models) {
		if (verbose==TRUE) print("Computing mean model (mean) ...")
		m.mean <- paste(ifelse(is.cubic,polycubic,poly),
			"b3==0",
			"b4==0",
			"b5==0",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"b1 == b2",
			"a1 := b1+b2",
			"a2 := 0",
			"a3 := b1-b2",
			"a4 := 0",
			add, sep="\n")
			s.mean <- sem(m.mean, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
		#summary(s.mean, fit.measures=TRUE)
	}

	if ("IA" %in% models) {
		if (verbose==TRUE) print("Computing interaction model (IA)...")
		m.IA <- paste(ifelse(is.cubic,polycubic,poly),
			"b3==0",
			"b5==0",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p10 := Y0 - p11*X0",			
			"p20 := Y0 - p21*X0",
			"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
		add, sep="\n")
		
			s.IA <- sem(m.IA, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("SQD" %in% models) {
		if (verbose==TRUE) print("Computing squared difference model (SQD) ...")
		m.SQD <- paste(ifelse(is.cubic,polycubic,poly),
			"b1==0",
			"b2==0",
			"b3==b5",
			"b3+b4+b5==0",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")			
			s.SQD <- sem(m.SQD, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("SSQD" %in% models) {
		if (verbose==TRUE) print("Computing shifted squared difference model (SSQD) ...")
		m.SSQD <- paste(ifelse(is.cubic,polycubic,poly),
			"b1==-b2",
			"b3==b5",
			"b3+b4+b5==0",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
			"C := b1 / (2*b3)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")			
			s.SSQD <- sem(m.SSQD, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if (any(models %in% c("RR"))) {
		if (verbose==TRUE) print("Computing rising ridge model (RR) ...")
		m.RR <- paste(ifelse(is.cubic,polycubic,poly),
			"b1==b2",
			"b3==b5",
			"b3+b4+b5==0",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
			"meaneffect := b1+b2",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			
			add, sep="\n")
			s.RR <- sem(m.RR, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if (any(models %in% c("SRR"))) {
		if (verbose==TRUE) print("Computing shifted rising ridge model (SRR) ...")
		m.SRR <- paste(ifelse(is.cubic,polycubic,poly),
			"b3==b5",
			"b3+b4+b5==0",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
			"meaneffect := a1",
			"C := (b1-b2) / (4*b3)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			
			add, sep="\n")
			s.SRR <- sem(m.SRR, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	
	if (any(models %in% c("SRRR"))) {
			if (verbose==TRUE) print("Computing rotated and shifted rising ridge model (SRRR), up ...")
			m.SRRR.up <- paste(paste(ifelse(is.cubic,polycubic,poly), " + start(0.01)*", IV12, " + start(0.01)*", IV22),
				"b3 > 0.000001",
				"b5 > 0.000001",
				"b4^2 == 4*b3*b5",
				ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
				"a1 := b1+b2",
				"a2 := b3+b4+b5",
				"a3 := b1-b2",
				"a4 := b3-b4+b5",
				"a5 := b3-b5",
				"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
				"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
				"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
				"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
				"meaneffect := (b2*b4 - 2*b1*b5) / b4",
				"C := (-2*b1*b5 - b2*b4) / (4*b4*b5)",
				"S := (-b4) / (2*b5)",
				"a4.rescaled := b3/S^2 - b4/S + b5",
				# eigenvalues
				"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
				"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
		
				add, sep="\n")
				s.SRRR.up <- sem(m.SRRR.up, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)	
				
			if (verbose==TRUE) print("Computing rotated and shifted rising ridge model (SRRR), down ...")
			m.SRRR.down <- paste(paste(ifelse(is.cubic,polycubic,poly), " + start(-0.01)*", IV12, " + start(-0.01)*", IV22),
			#m.SRRR <- paste(paste(ifelse(is.cubic,polycubic,poly), " + start(-0.001)*", IV22),
				"b3 < -0.000001",
				"b5 < -0.000001",
				"b4^2 == 4*b3*b5",
			  ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
				"a1 := b1+b2",
				"a2 := b3+b4+b5",
				"a3 := b1-b2",
				"a4 := b3-b4+b5",
				"a5 := b3-b5",
				"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
				"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
				"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
				"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
				"meaneffect := (b2*b4 - 2*b1*b5) / b4",
				"C := (-2*b1*b5 - b2*b4) / (4*b4*b5)",
				"S := (-b4) / (2*b5)",
				"a4.rescaled := b3/S^2 - b4/S + b5",
				# eigenvalues
				"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
				"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
	
				add, sep="\n")
				s.SRRR.down <- sem(m.SRRR.down, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)	
				
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
	
	
	if (any(models %in% c("SRSQD"))) {
		if (verbose==TRUE) print("Computing rotated squared difference model (SRSQD), up ...")
		m.SRSQD.up <- paste(paste(ifelse(is.cubic,polycubic,poly), " + start(0.001)*", IV22),
		  "2*b1*b5 == b2*b4",	# this is a different (but algebraically equivalent) formulation of the constraints
			"b3 > 0.000001",
			"b5 > 0.000001",
			"b4^2 == 4*b3*b5",	# this is a different (but algebraically equivalent) formulation of the constraints
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"C := -.5*(b2/b5)",
			"S := (-b4) / (2*b5)",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
			"a4.rescaled := b3/S^2 - b4/S + b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p10 := Y0 - p11*X0",			
			"p20 := Y0 - p21*X0",
			"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")
			s.SRSQD.up <- sem(m.SRSQD.up, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
			
			
		if (verbose==TRUE) print("Computing rotated squared difference model (SRSQD), down ...")
		m.SRSQD.down <- paste(paste(ifelse(is.cubic,polycubic,poly), " + start(-0.001)*", IV22),
			"2*b1*b5 == b2*b4",	# this is a different (but algebraically equivalent) formulation of the constraints
			"b3 < -0.000001",
			"b5 < -0.000001",
			"b4^2 == 4*b3*b5",
			ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"C := -.5*(b2/b5)",
			"S := (-b4) / (2*b5)",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
			"a4.rescaled := b3/S^2 - b4/S + b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")
			s.SRSQD.down <- sem(m.SRSQD.down, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	
			if (inspect(s.SRSQD.up, "converged") == FALSE & inspect(s.SRSQD.down, "converged") == TRUE) {
				SRSQD.rot <- "down"
			} else 
			if (inspect(s.SRSQD.up, "converged") == TRUE & inspect(s.SRSQD.down, "converged") == FALSE) {
				SRSQD.rot <- "up"
			} else 
			if (inspect(s.SRSQD.up, "converged") == TRUE & inspect(s.SRSQD.down, "converged") == TRUE) {
				SRSQD.rot <- ifelse(fitMeasures(s.SRSQD.up, "chisq") > fitMeasures(s.SRSQD.down, "chisq"), "down", "up")
			} else {
				if (verbose==TRUE) warning("Warning: SRSQD model has not converged (neither up nor down curvature)")
			}
			if (SRSQD.rot == "up") {
				s.SRSQD <- s.SRSQD.up
			} else 
			if (SRSQD.rot == "down") {
				s.SRSQD <- s.SRSQD.down
			}
			if (verbose == TRUE) print(paste0("Direction of SRSQD curvature: ", SRSQD.rot))

		
	}
	
	
	# always estimate the full second-order model (= global model for all second-order models)
		if (verbose==TRUE) print("Computing polynomial model (full) ...")
		m.full <-  paste(ifelse(is.cubic,polycubic,poly),
		  ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add,
			sep="\n"
		)
		s.full <- sem(m.full, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	
	
	if ("weak" %in% models) {
		if (verbose==TRUE) print("Computing weak fit pattern ...")
		m.weak <-  paste(ifelse(is.cubic,polycubic,poly),
		  ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			# constraints for weak pattern
			"b3*b5 > 0",	# must be > 0
			add,
			sep="\n"
		)
		s.weak <- sem(m.weak, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	
	if ("strong" %in% models) {
		if (verbose==TRUE) print("Computing strong fit pattern ...")
		m.strong <-  paste(ifelse(is.cubic,polycubic,poly),
		  ifelse(is.cubic, paste("b6==0","b7==0","b8==0","b9==0", sep="\n"), ""),
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a5 := b3-b5",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"PA1.curv := b3 + b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 + b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			# constraints for strong pattern
			"b3*b5 > 0.000001",			# must be > 0
			"(b2*b4) == 2*b1*b5",	# must be 0
			"4*b3*b5  == b4^2",	# must be 0
			add,
			sep="\n"
		)
		s.strong <- sem(m.strong, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	
	#m.absdiff.JRE <-  paste(
	#	paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, " + 0*W.JRE + w2*W.JRE_", IV1, " + w3*W.JRE_", IV2),
	#	"b1 == -b2",
	#	"w2 == -w3",
	#	"w2 == -2*b1",
	#	add, sep="\n")
	#s.absdiff.JRE <-  sem(m.absdiff.JRE, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	#summary(s.absdiff.JRE, fit.measures=TRUE)
	
	# the unconstrained absolute difference model - Edwards (2002) formula
	#m.absunc.JRE <-  paste(
	#	paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, " + w1*W.JRE + w2*W.JRE_", IV1, " + w3*W.JRE_", IV2),
	#	add, sep="\n")
	#s.absunc.JRE <-  sem(m.absunc.JRE, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	#summary(s.absunc.JRE, fit.measures=TRUE)
	
	
	if ("absdiff" %in% models) {
		if (verbose==TRUE) print("Computing constrained absolute difference model (absdiff) ...")
		m.absdiff <-  paste(
			paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + w1*W + w2*W_", IV1, " + w3*W_", IV2),
			"b1 == 0",
			"b2 == 0",
			"w1 == 0",
			"w2 == -w3",
			add, sep="\n")
			
			s.absdiff <- sem(m.absdiff, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("absunc" %in% models) {
		# the unconstrained absolute difference model - new formula
		if (verbose==TRUE) print("Computing unconstrained absolute difference model (absunc) ...")
		m.absunc <-  paste(
			paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + w1*W + w2*W_", IV1, " + w3*W_", IV2),
			ifelse(breakline==FALSE, "w1==0", ""),
			add, sep="\n")
			
			s.absunc <- sem(m.absunc, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}

	
	
	
	## Test all models of third degree

	# Full model of third degree
	if ("cubic" %in% models | is.cubic) {
		if (verbose==TRUE) print("Computing full cubic model (cubic) ...")
		m.cubic <-  paste(polycubic,
			# description of surface above LOC and LOIC
			"u1 := b1 + b2",				# linear term coefficient of LOC
			"u2 := b3 + b4 + b5",			# quadratic term coefficient of LOC
			"u3 := b6 + b7 + b8 + b9",	# cubic term coefficient of LOC
			"v1 := b1 - b2",				# linear term coefficient of LOIC
			"v2 := b3 - b4 + b5",			# quadratic term coefficient of LOIC
			"v3 := b6 - b7 + b8 - b9",	# cubic term coefficient of LOIC
			# user-defined syntax
			add,
			sep="\n"
		)
		s.cubic <- sem(m.cubic, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)		
	}
	
	# "Cubic Difference" models of third degree
	if ("CA" %in% models) {
		if (verbose==TRUE) print("Computing cubic asymmetry model (CA) ...")
		m.CA <-  paste(polycubic,
      # constraints
			"b1 == 0",
      "b2 == 0",
			"b4 == -2*b3",
      "b5 == b3",
			"b7 == -3*b6",
			"b8 == 3*b6",
			"b9 == -b6",
			# description of surface above LOC and LOIC
			"u1 := b1 + b2",
			"u2 := b3 + b4 + b5",
			"u3 := b6 + b7 + b8 + b9",
			"v1 := b1 - b2",
			"v2 := b3 - b4 + b5",
			"v3 := b6 - b7 + b8 - b9",
			# user-defined syntax
			add,
			sep="\n"
		)
		s.CA <- sem(m.CA, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)		
	}
	
	
	if ("CL" %in% models) {
	  if (verbose==TRUE) print("Computing cubic level model (CL) ...")
	  m.CL <-  paste(polycubic,
       # constraints
       "b1 == 0",
       "b2 == 0",
       "b4 == -2*b3",
       "b5 == b3",
       "b7 == -b6",
       "b8 == -b6",
       "b9 == b6",
       # description of surface above LOC and LOIC
       "u1 := b1 + b2",
       "u2 := b3 + b4 + b5",
       "u3 := b6 + b7 + b8 + b9",
       "v1 := b1 - b2",
       "v2 := b3 - b4 + b5",
       "v3 := b6 - b7 + b8 - b9",
       # user-defined syntax
       add,
       sep="\n"
	  )
	  s.CL <- sem(m.CL, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)		
	}

	
	if ("RRCA" %in% models) {
		if (verbose==TRUE) print("Computing rising ridge cubic asymmetry model (RRCA) ...")
		m.RRCA <-  paste(polycubic,
      # constraints
			"b1 == b2",
			"b4 == -2*b3",
      "b5 == b3",
			"b7 == -3*b6",
			"b8 == 3*b6",
			"b9 == -b6",
			# description of surface above LOC and LOIC
			"u1 := b1 + b2",
			"u2 := b3 + b4 + b5",
			"u3 := b6 + b7 + b8 + b9",
			"v1 := b1 - b2",
			"v2 := b3 - b4 + b5",
			"v3 := b6 - b7 + b8 - b9",
			# special parameters in this model
			"meaneffect := u1",
			# user-defined syntax
			add,
			sep="\n"
		)
		s.RRCA <- sem(m.RRCA, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)		
	}	
	
	
	if ("RRCL" %in% models) {
	  if (verbose==TRUE) print("Computing rising ridge cubic level model (RRCL) ...")
	  m.RRCL <-  paste(polycubic,
	                   # constraints
	                   "b1 == b2",
	                   "b4 == -2*b3",
	                   "b5 == b3",
	                   "b7 == -b6",
	                   "b8 == -b6",
	                   "b9 == b6",
	                   # description of surface above LOC and LOIC
	                   "u1 := b1 + b2",
	                   "u2 := b3 + b4 + b5",
	                   "u3 := b6 + b7 + b8 + b9",
	                   "v1 := b1 - b2",
	                   "v2 := b3 - b4 + b5",
	                   "v3 := b6 - b7 + b8 - b9",
	                   # special parameters in this model
	                   "meaneffect := u1",
	                   # user-defined syntax
	                   add,
	                   sep="\n"
	  )
	  s.RRCL <- sem(m.RRCL, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)		
	}	
	
	
	},	  # end of "withCallingHandlers"



# suppress several types of warning
  warning=function(w) {
	   W <- as.character(w$call)
	   if (
			 (W[1] == "sqrt") |
		   (W[1] == "nlminb" & W[2] == "x.par") |
		   (W[2] %in% c("m.SRRR.up", "m.SRRR.down", "m.SRSQD.up", "m.SRSQD.down") & grepl("model has NOT converged", w$message))
		  ) {invokeRestart("muffleWarning")}
} )


	# ---------------------------------------------------------------------
	# Sanity check: Check results for convergence problems
	# Sometimes the SRRR or the SRSQD model find a bad solution and have a higher chi2 than their nested models (which is theoretically not possible)
	chisq1 <- plyr::ldply(list(full=s.full, SRRR=s.SRRR, SRR=s.SRR, RR=s.RR, SQD=s.SQD), function(x) {
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
	
	chisq2 <- plyr::ldply(list(full=s.full, SRRR=s.SRRR, SRSQD=s.SRSQD, SSQD=s.SSQD, SQD=s.SQD), function(x) {
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
	modellist <- list(null=s.NULL, full=s.full, IA=s.IA, diff=s.diff, mean=s.mean, absdiff=s.absdiff, additive=s.additive, SQD=s.SQD, SRRR=s.SRRR, SRR=s.SRR, RR=s.RR, SSQD=s.SSQD, SRSQD=s.SRSQD, absunc=s.absunc, cubic=s.cubic, onlyx=s.onlyx, onlyy=s.onlyy, onlyx2=s.onlyx2, onlyy2=s.onlyy2, weak=s.weak, strong=s.strong, 
	CA=s.CA, CL=s.CL, RRCA=s.RRCA, RRCL=s.RRCL)
	
	res <- list(
		models = modellist, 
		SRSQD.rot = SRSQD.rot, SRRR.rot = SRRR.rot, LM=summary(lm.full), formula=formula, 
		data=df, out.rm = out.rm, outliers = which(df$out == TRUE), DV=DV, IV1=IV1, IV2=IV2, IV12=IV12, IV22=IV22, 
		IV_IA=IV_IA, W_IV1=W_IV1, W_IV2=W_IV2, IV13=IV13, IV_IA2=IV_IA2, IV_IA3=IV_IA3, IV23=IV23, 
		control.variables = control.variables, is.cv = is.cv, 
		r.squared = summary(lm.full)$r.squared, 
		is.cubic=is.cubic)
	
	attr(res, "class") <- "RSA"
	return(res)
}
