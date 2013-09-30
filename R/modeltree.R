#' @title Plots a flow chart with model comparisons
#'
#' @description
#' Plots a flow chart with model comparisons from a RSA object
#'
#' @details
#' The plot can be either requested within the \code{compare} function:
#' \code{compare(r1, plot=TRUE)}
#' Or it can be plotted from a cRSA object (= output from the \code{\link{compare}} function):
#' \code{c1 <- compare(r1)}
#' \code{plot(c1)}
#'
#' @export
#' @param x A cRSA object (= output from the \code{\link{compare}} function)
#' @param digits The number of digits to which numbers are rounded
#' @param ... Additional parameters (not used yet)
#'
#' @seealso \code{\link{RSA}}, \code{\link{compare}}
#'
modeltree <- function(x, digits=3, ...) {
	library(qgraph)

	c1 <- x
	c1$fromto <- c("", paste0(c1$model[1:(nrow(c1)-1)], "_", c1$model[2:(nrow(c1))]))
	
	# define labels of boxes
	m <- c("full", "SRRR", "IA", "additive", "diff", "SRR", "RR", "sqdiff", "SRSD", "SSD", "null")

	# get values from compare-object
	R2.adj <- c()
	CFI <- c()
	R2.p <- c()
	for (i in 1:length(m)) {
		R2.adj <- c(R2.adj, c1[c1$model==m[i], "R2.adj"][1])
		CFI <- c(CFI, c1[c1$model==m[i], "cfi"][1])
		R2.p <- c(R2.p, c1[c1$model==m[i], "R2.p"][1])
	}

	pos <- matrix(ncol=4, byrow=TRUE, c(
	12, 0, 1, 0,
	13, 0, 2, 0,	
	14, 3, 6, 9,
	15, 4, 7,10,
	16, 5, 0, 8,
	17, 0,11, 0
	))

	# define edgelist, without weight
	eL <- matrix(ncol=2, byrow=TRUE, data=c(
	1,2,
	1,3,
	3,4,
	4,5,
	5,11,
	2,6,
	6,7,
	7,8,
	8,11,
	2,9,
	9,10,
	10,8,
	6,10,
	
	12,12,
	13,13,
	14,14,
	15,15,
	16,16,
	17,17
	))

	# define weights of edges: will be translated to color
	w <- c()
	for (i in 1:13) {
	    w <- c(w, c1[c1$fromto == paste0(m[eL[i, 1]], "_", m[eL[i, 2]]), "Pr(>Chisq)"][1])
	}
	w[is.na(w)] <- 0

	# compute box labels
	lab <- c("full", "SRRR", "Interaction", "Additive", "Difference", "SRR", "RR", "Squared difference", "SRSD", "SSD", "Intercept only", "k = 5", "k = 4", "k = 3", "k = 2", "k = 1", "k = 0")	
	lab2 <- list()
	lab.color <- c()
	
	for (i in 1:length(lab)) {
		if (i < 12) {
			lab2[[i]] <- paste(lab[i],
				paste("CFI = ", f2(CFI[i], digits)),
				paste("R^2[adj] = ", f2(R2.adj[i], digits)), sep="\n")
			lab.color <- c(lab.color, ifelse(R2.p < .05, "black", "grey"))	
		} else {
			lab2[[i]] <- lab[i]
			lab.color <- c(lab.color, "black")
		}
	}
	
	
	qgraph(eL, 
		edgeList	= TRUE,
		nNodes		= 17,
		layout 		= pos, 
	
		# define edges
		edge.labels = paste0("p = ", f2(w, digits)), 
		edge.color	= c(pRamp(w), rep("#FFFFFF", 6)), 
		edge.label.cex = 0.8,
	
		# define boxes
		labels		= lab2,
		label.color	= lab.color,
		shape		= "rectangle",
		border.color = c(rep("black", 11), rep("white", 6)),
		asize		= 3,		# size of arrowhead
		vsize 		= 14,		# horizontal size of boxes
		vsize2 		= 7,		# vertical size of boxes
		border.width = R2.adj
	)
}

#' @S3method plot cRSA
plot.cRSA <- function(x, ...) {
	plot.cRSA(x, ...)
}