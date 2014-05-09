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
#' @param plot Should a plot of the AICc table be plotted? (Experimental)

aictab <- function(x, plot=FALSE) {
	cand.set <- x$models[!names(x$models) %in% c("absdiff", "absunc")]
	cand.set <- cand.set[!unlist(lapply(cand.set, is.null))]
	a1 <- aictab.lavaan(cand.set, modnames=names(cand.set))
	
	if (plot==TRUE) {
		a2 <- a1
		a2 <- a2[order(a2$Delta_AICc, decreasing=TRUE), ]
		a2$Modnames <- factor(a2$Modnames, levels=a2$Modnames, ordered=TRUE)
		a2$color <- ifelse(a2$Delta_AICc <= 2, "green", "orange")
		a2$color[a2$Cum.Wt > .95] <- "gray"
		a2$color <- factor(a2$color, levels=c("green", "orange", "gray"))
		p1 <- ggplot(a2, aes_string(x="Modnames", y="Delta_AICc", group=1)) + geom_line() + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + geom_hline(yintercept=2, linetype="dashed") + geom_hline(yintercept=7, linetype="dotted") + geom_point(aes_string(color="color"), size=3) + geom_vline(xintercept=which(a2$Cum.Wt < .95)[1] - 0.5, linetype="solid", color="grey")  + annotate("text", label="Cumulative weight > .95", y=max(a2$Delta_AICc), x=which(a2$Cum.Wt < .95)[1]-1, size=3, hjust=1) + annotate("text", label="Practically equivalent models", y=1, x=0.5, size=3, hjust=0, angle=90) + annotate("text", label="Implausible models", y=10, x=0.5, size=3, hjust=0, angle=90) + annotate("text", label="Cumulative weight < .95", y=max(a2$Delta_AICc), x=which(a2$Cum.Wt < .95)[1], size=3, hjust=1) + coord_flip(ylim=c(-0.5, max(a2$Delta_AICc)+1)) + ylab("Delta AIC") + xlab("Model")
		# + scale_color_discrete("", levels=1:3, values=c("green", "orange", "gray"))
		
		print(p1)
	}
	return(a1)
}


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
	if(is.null(modnames)) modnames<-1:length(cand.set)
	# check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
   # if (length(unique(check.resp)) > 1) 
   #     stop("You must use the same response variable for all models\n")
    Results <- NULL
    Results <- data.frame(Modnames = modnames)
    Results$K <- unlist(lapply(X = cand.set, FUN = AICc.lavaan, 
        return.K = TRUE, c.hat = c.hat,second.ord = second.ord))
    Results$AICc <- unlist(lapply(X = cand.set, FUN = AICc.lavaan, 
        return.K = FALSE, c.hat = c.hat,second.ord = second.ord))
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
