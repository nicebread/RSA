#' @title Show a table of AIC model comparisons
#' @aliases aictab
#'
#' @description
#' Show a table of AIC model comparisons
#'
#' @return
#' \describe{
#'  \item{Modnames}{Model names.}
#'  \item{K}{Number of estimated parameters (including the intercept, residual variance, and, if present in the model, control variables).}
#'  \item{LL}{Model log-likelihood.}
#'  \item{AICc}{Akaike Information Criterion (corrected).}
#'  \item{Delta_AICc}{Difference in AICc between this model and the best model.}
#'  \item{AICcWt}{The Akaike weights, also termed "model probabilities" by Burnham and Anderson (2002). Indicates the level of support (i.e., weight of evidence) of a model being the most parsimonious among the candidate model set.}
#'  \item{Cum.Wt}{Cumulative Akaike weight. One possible strategy is to restrict interpretation to the "confidence set" of models, that is, discard models with a Cum.Wt > .95 (see Burnham & Anderson, 2002, for details and alternatives).}
#'  \item{evidence.ratio}{Likelihood ratio of this model vs. the best model.}
#'  \item{cfi}{Comparative Fit Index (CFI).}
#'  \item{R2}{Coefficient of determination (R-squared).}
#'  \item{R2.adj}{Adjusted R-squared.}
#'  \item{R2.baseline}{Only provided if the model contains control variables. Difference in R-squared as compared to the baseline model with intercept and control variables (= the model "null"). This R^2 increment will typically be of interest because it refers to the amount of variance explained by the two predictors X and Y (plus their squared and interaction terms) in the RSA model.}
#'  \item{R2.baseline.p}{Only provided if the model contains control variables. p-value for the F-test of the model against the baseline model.}
#' }                       
#'
#' @references
#' Burnham, K. P., & Anderson, D. R. (2002). \emph{Model selection and multimodel inference: A practical information-theoretic approach.} Springer Science & Business Media.
#' 
#' @export
#' @param x An RSA object
#' @param models A vector with all model names of the candidate set. Defaults to all polynomial models in the RSA object.
#' @param plot Should a plot of the AICc table be plotted?
#' @param bw Should the plot be black & white?
#' @param digits The output is rounded to this number of digits. No rounding if NA (default).
#' 
#' @note This function is similar to the function \code{\link[AICcmodavg]{aictab}} in the \code{AICcmodavg} package.
#' 
#' @examples
#' \dontrun{
#' data(motcon)
#' r.m <- RSA(postVA~ePow*iPow, motcon, verbose=FALSE)
#' aictab(r.m, plot=TRUE)
#' }

aictab <- function(x, plot=FALSE, bw=FALSE, models=names(x$models)[!names(x$models) %in% c("absdiff", "absunc")], digits=NA) {
	cand.set.models <- x$models[models]
	
	# remove NULL models
	cand.set.models <- cand.set.models[!unlist(lapply(cand.set.models, is.null))]	
	
	# remove non-converged models
	cand.set.models <- cand.set.models[unlist(lapply(cand.set.models, inspect, "converged"))]
	
	a1 <- aictab.lavaan(cand.set.models, modnames=names(cand.set.models))
	class(a1) <- "data.frame"
	a1 <- a1[, c("Modnames", "K", "LL", "AICc", "Delta_AICc", "AICcWt", "Cum.Wt")]
	a1$evidence.ratio <- evidenceRatio(a1$Delta_AICc)
	a1$evidence.ratio[1] <- NA	

	c1 <- plyr::ldply(names(cand.set.models), function(m) {

	  Modnames <- m

	  # CFI
	  cfi <- fitmeasures(x$models[[m]], fit.measures = "cfi")

	  # R2 and F-test of this model versus intercept only model
	  R2.vs.interceptonly <- R2difftest(x, unrestricted=m, restricted="interceptonly")
	  R2 <- R2.vs.interceptonly$delta.R2
	  R2.p <- R2.vs.interceptonly$R2.p
	  
	  # number of parameters (without intercept and residual variance) for computation of adjusted R^2
	  object <- as.list(fitMeasures(x$models[[m]]))
	  df <- object$baseline.df - object$df

	  return( data.frame(Modnames=Modnames, cfi=as.vector(cfi), R2=as.vector(R2), R2.p=R2.p, df) )
	})
	
	# un-factor
	c1$Modnames <- as.character(levels(c1$Modnames))[c1$Modnames]
	
	# adjusted R2
	N <- lavaan::nobs(x$models[["full"]])
	c1$R2.adj <- 1 - ((1-c1$R2))*((N-1)/(N-c1$df-1))

	# merge aic and other fit indices
	a2 <- merge(a1, c1[, c("Modnames", "cfi", "R2", "R2.p", "R2.adj")], by="Modnames")

	# if control variables are present: include delta R2 and F-tests the models versus the baseline model (intercept + control variables)
	if(x$is.cv){

	  c1.cv <- plyr::ldply(names(cand.set.models), function(m) {

	    Modnames <- m

	    R2.vs.baseline <- R2difftest(x, unrestricted=m, restricted="null")
	    R2.baseline <- R2.vs.baseline$delta.R2
	    R2.baseline.p <- R2.vs.baseline$R2.p

	    return( data.frame(Modnames=Modnames, R2.baseline=as.vector(R2.baseline), R2.baseline.p=R2.baseline.p) )
	  })

	  # un-factor
	  c1.cv$Modnames <- as.character(levels(c1.cv$Modnames))[c1.cv$Modnames]

	  a2 <- merge(a2, c1.cv, by="Modnames")
	}

	# order by aic and round
	a2 <- a2[order(a2$Delta_AICc), ]
	if (!is.na(digits)) {
	  a2[, -c(1:2)] <- round(a2[, -c(1:2)], digits)
	}
	
	
	if (plot==TRUE) {
		a3 <- a2
		a3$Modnames <- factor(a3$Modnames, levels=a3$Modnames, ordered=TRUE)
		a3$ValidModels <- as.character(a3$Modnames)
		a3$ValidModels[a3$Delta_AICc > 2] <- ""
		a3$color <- ifelse(a3$Delta_AICc <= 2, "green", "yellow")
		a3$color[a3$Cum.Wt > .95] <- "red"
		a3$color[a3$Delta_AICc > 10] <- "red"
		
		if (bw==TRUE) {
			col_scale <- c("green"="grey10", "yellow"="grey10", "red"="grey10")
			shape_scale <- c("green"=19, "yellow"=1, "red"=4)
		} else {			
			col_scale <- c("green"="green2", "yellow"="darkgoldenrod1", "red"="firebrick3")
			shape_scale <- c("green"=19, "yellow"=19, "red"=19)
		}
		
		
		p1 <- ggplot(a3, aes_string(y="AICcWt", x="Delta_AICc", group=1)) + geom_line(color="grey30") + theme_bw() + geom_vline(xintercept=2, linetype="dotted", color="grey30") + geom_vline(xintercept=7, linetype="dotted", color="grey30") + geom_point(aes_string(color="color", shape="color"), size=3) + coord_cartesian(xlim=c(-0.5, max(a3$Delta_AICc)+1)) + xlab(bquote(Delta~"AICc")) + ylab("Model weight") + geom_text(aes_string(label="ValidModels"), size=3.7, hjust=-0.2) + scale_x_continuous(breaks=c(2, 7, 10)) + scale_color_manual(values=col_scale, guide="none") + scale_shape_manual(values=shape_scale, guide="none")		
		
		p1 <- p1 + 
			geom_hline(yintercept=a3$AICcWt[max(which(a3$Cum.Wt < .95))], linetype="dotted", color="grey30") + 
			annotate("text", label="Cumulative weight < .95", x=max(a3$Delta_AICc), y=a3$AICcWt[max(which(a3$Cum.Wt < .95))], size=3.3, hjust=1, vjust=-.4) + 
			annotate("text", label="Cumulative weight > .95", x=max(a3$Delta_AICc), y=a3$AICcWt[max(which(a3$Cum.Wt < .95))], size=3.3, hjust=1, vjust=1.4) 
		
		p1 <- p1 + 
			annotate("text", label="Practically equivalent models", x=0.5, y=0.03, size=3.3, hjust=0, vjust=0, angle=90) + 
			annotate("text", label="Implausible models", x=7.5, y=0.03, size=3.3, hjust=0, vjust=0, angle=90)
		
		print(p1)
	}

	return(a2)
}


evidenceRatio <- function(Delta.AICc) {1/exp(-0.5*Delta.AICc)}

# adapted from: http://byrneslab.net/classes/lavaan_materials/lavaan.modavg.R
AICc.lavaan<-function(object, second.ord=TRUE, c.hat = 1, return.K = FALSE){
	object <- as.list(fitMeasures(object))
	npar <- object$baseline.df - object$df + 2 # number of parameters, including coefficients of control variables (if there are any), intercept and residual variance (thus the +2)
	if(return.K == TRUE) return(npar)
	if(second.ord == FALSE && c.hat>1) return(-2*object$logl/c.hat+2*npar)
	if(second.ord == FALSE) return(-2*object$logl + 2*npar)
    if(c.hat>1) return( -2*object$logl/c.hat+2*npar + 2*( npar*(object$npar+1))/(object$ntotal-npar-1))
	  -2*object$logl + 2*npar + 2*( npar*(npar+1))/(object$ntotal-npar-1)
}
    
aictab.lavaan<-function(cand.set, modnames, sort = TRUE, c.hat = 1, second.ord = TRUE, nobs = NULL){
	if (is.null(modnames)) modnames<-1:length(cand.set)
	# check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
   # if (length(unique(check.resp)) > 1) 
   #     stop("You must use the same response variable for all models\n")
    Results <- NULL
    Results <- data.frame(Modnames = modnames)
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc.lavaan, return.K = TRUE, c.hat = c.hat, second.ord = second.ord))
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc.lavaan, return.K = FALSE, c.hat = c.hat, second.ord = second.ord))
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)
    Results$ModelLik <- exp(-0.5 * Results$Delta_AICc)
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)
    if (length(unique(Results$AICc)) != length(cand.set)) 
        warning("\nCheck model structure carefully as some models may be redundant\n")
    if (second.ord == TRUE && c.hat == 1) {
        Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
    }
    if (second.ord == TRUE && c.hat > 1) {
        colnames(Results) <- c("Modnames", "K", "QAICc", "Delta QAICc", 
            "ModelLik", "QAICcWt")
        LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
        Results$Quasi.LL <- LL/c.hat
        Results$c_hat <- c.hat
    }
    if (second.ord == FALSE && c.hat == 1) {
        colnames(Results) <- c("Modnames", "K", "AIC", "Delta AIC", 
            "ModelLik", "AICWt")
        Results$LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
    }
    if (second.ord == FALSE && c.hat > 1) {
        colnames(Results) <- c("Modnames", "K", "QAIC", "Delta QAIC", 
            "ModelLik", "QAICWt")
        LL <- unlist(lapply(X = cand.set, FUN = function(i) logLik(i)[1]))
        Results$Quasi.LL <- LL/c.hat
        Results$c_hat <- c.hat
    }
    if (sort) {
        Results <- Results[rev(order(Results[, 6])), ]
        Results$Cum.Wt <- cumsum(Results[, 6])
    }
    else {
        Results$Cum.Wt <- NULL
    }
    class(Results) <- c("aictab", "data.frame")
    return(Results)
	
}
