#' @title Show a table of AIC model comparisons
#' @aliases aictab
#'
#' @description
#' Show a table of AIC model comparisons
#'
#' @details
#' TODO
#'
#' @export
#' @param x An RSA object
#' @param cand.set A vector with all model names of the candidate set. Defaults to all polynomial models in the RSA object.
#' @param plot Should a plot of the AICc table be plotted? (Experimental)

aictab <- function(x, plot=FALSE, cand.set=names(x$models)[!names(x$models) %in% c("absdiff", "absunc")]) {
	cand.set.models <- x$models[cand.set]
	cand.set.models <- cand.set.models[!unlist(lapply(cand.set.models, is.null))]
	a1 <- aictab.lavaan(cand.set.models, modnames=names(cand.set.models))
	
	if (plot==TRUE) {
		a2 <- a1
		a2$Modnames <- factor(a2$Modnames, levels=a2$Modnames, ordered=TRUE)
		a2$ValidModels <- as.character(a2$Modnames)
		a2$ValidModels[a2$Delta_AICc > 2] <- ""
		a2$color <- ifelse(a2$Delta_AICc <= 2, "green", "orange")
		a2$color[a2$Cum.Wt > .95] <- "gray"
		a2$color <- factor(a2$color, levels=c("green", "orange", "gray"))
		
		p1 <- ggplot(a2, aes_string(y="AICcWt", x="Delta_AICc", group=1)) + geom_line() + theme_bw() + geom_vline(xintercept=2, linetype="dashed") + geom_vline(xintercept=7, linetype="dotted") + geom_point(aes_string(color="color"), size=3) + coord_cartesian(xlim=c(-0.5, max(a2$Delta_AICc)+1)) + xlab("Delta AIC") + ylab("Model weight") + geom_text(aes_string(label="ValidModels"), size=3, hjust=-0.3) + scale_x_continuous(breaks=c(2, 7, 10))
		
		p1 <- p1 + geom_hline(yintercept=a2$AICcWt[max(which(a2$Cum.Wt < .95))], linetype="solid", color="grey") + annotate("text", label="Cumulative weight > .95", x=max(a2$Delta_AICc), y=a2$AICcWt[max(which(a2$Cum.Wt < .95))], size=3, hjust=1, vjust=-.4) + annotate("text", label="Cumulative weight < .95", x=max(a2$Delta_AICc), y=a2$AICcWt[max(which(a2$Cum.Wt < .95))], size=3, hjust=1, vjust=1.4) 
		
		p1 <- p1 + annotate("text", label="Practically equivalent models", x=1, y=0, size=3, hjust=0, vjust=0, angle=90) + annotate("text", label="Implausible models", x=10, y=0, size=3, hjust=0, vjust=0, angle=90)
		
		# + scale_color_discrete("", levels=1:3, values=c("green", "orange", "gray"))
		
		print(p1)
	}
	return(a1)
}


# from: http://byrneslab.net/classes/lavaan_materials/lavaan.modavg.R
AICc.lavaan<-function(object, second.ord=TRUE, c.hat = 1, return.K = FALSE){
	object <- as.list(fitMeasures(object))
	npar<-object$baseline.df - object$df
	if(return.K == TRUE) return(object$npar)
	if(second.ord == FALSE && c.hat>1) return(-2*object$logl/c.hat+2*npar)
	if(second.ord == FALSE) return(object$aic)
    if(c.hat>1) return( -2*object$logl/c.hat+2*npar + 2*( npar*(object$npar+1))/(object$ntotal-npar-1))
    object$aic + 2*( npar*(npar+1))/(object$ntotal-npar-1)
}
    
aictab.lavaan<-function(cand.set, modnames, sort = TRUE, c.hat = 1, second.ord = TRUE, nobs = NULL){
	if (is.null(modnames)) modnames<-1:length(cand.set)
	# check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
   # if (length(unique(check.resp)) > 1) 
   #     stop("You must use the same response variable for all models\n")
    Results <- NULL
    Results <- data.frame(Modnames = modnames)
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc.lavaan, return.K = TRUE, c.hat = c.hat, second.ord = second.ord))
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc.lavaan, return.K = FALSE, c.hat = c.hat,second.ord = second.ord))
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
