#' @title Retrieves several variables from an RSA object
#'
#' @description
#' Retrieves several variables from an RSA object
#'
#' @details
#' None so far.
#'
#' @export
#' @param x RSA object
#' @param type One of: "syntax", "coef", "R2", "R2.adj", "free", "summary"
#' @param model A string specifying the model; defaults to "full"
#' @param digits Number of digits the output is rounded to; if NA, digits are unconstrained
#' @param ... Additional parameters passed to the extraction function
#'
#' @seealso \code{\link{RSA}}
#'
#' @examples
#'
#' set.seed(0xBEEF)
#' n <- 300
#' err <- 2
#' x <- rnorm(n, 0, 5)
#' y <- rnorm(n, 0, 5)
#' df <- data.frame(x, y)
#' df <- within(df, {
#' 	diff <- x-y
#' 	absdiff <- abs(x-y)
#' 	sqdiff <- (x-y)^2
#' 	z.sq <- sqdiff + rnorm(n, 0, err)
#' })
#' 
#' r1 <- RSA(z.sq~x*y, df, models=c("full", "SSD"))
#' getPar(r1, "syntax")
#' getPar(r1, "R2")
#' getPar(r1, "coef")


getPar <- function(x, type="coef", model="full", digits=NA, ...) {
	type <- tolower(type)
	type <- match.arg(type, c("syntax", "coef", "r2", "rsquared", "r.squared", "r2.p", "rsquared.p", "r.squared.p", "r2.adj", "rsquared.adj", "r.squared.adj", "npar", "free", "summary"))
	if (type=="syntax") {
		return(x$models[[model]]@Options$model)
	}
	if (type=="coef") {
		p1 <- parameterEstimates(x$models[[model]], ...)
		p1 <- data.frame(p1[p1$label != "", ])
		rownames(p1) <- p1$rhs
		p1 <- p1[, -c(1:3)]
		if (!is.na(digits))	p1[, -1] <- round(p1[, -1], digits)
		return(p1)
	}
	if (type %in% c("r2", "rsquared", "r.squared")) {
		return(inspect(x$models[[model]], "R2", ...))
	}
	if (type %in% c("r2.p", "rsquared.p", "r.squared.p")) {
		R <- inspect(x$models[[model]], "R2", ...)
		n <- nobs(x$models[[model]])
		k <- fitmeasures(x$models[[model]], "npar")
		return(pf(((n-k-1)*R)/(k*(1-R)), k, n-k-1, lower.tail=FALSE))
	}
	
	if (type %in% c("r2.adj", "rsquared.adj", "r.squared.adj")) {
		freeparam <- getFreeParameters(x$models[[model]]) - fitmeasures(x$models[[model]], "Df")
		r2.adj <- 1 - (1-inspect(x$models[[model]], "R2")) * ((nobs(x$models[[model]])-1)/(nobs(x$models[[model]]) - freeparam - 1))
		names(r2.adj) <- "r2.adj"
		return(r2.adj)
	}
	if (type %in% c("npar", "free")) {
		return(fitmeasures(x$models[[model]], "npar"))
	}
	if (type %in% c("summary")) {
		return(summary(x$models[[model]], ...))
	}
}
