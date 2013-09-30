#' @title Compute bootstrapped confidence interval and p value
#' @aliases CI.boot
#' @description
#' Compute bootstrapped confidence interval and p value
#'
#' @details
#' See \code{\link{bootRSA}} for more details and example.
#'
#' @export
#' @param r.boot Results object from the \code{\link{bootRSA}} function
#' @param CI Width of confidence interval
#'
#' @seealso \code{\link{RSA}}

CI.boot <- function(r.boot, CI=.95) {
	CIs <- apply(r.boot, 2, function(x) {
		p <- sum(x<0)/length(x)
		qu <- quantile(x, probs=c((1-CI)/2, 1-(1-CI)/2))
		res <- c(qu, p.value=min(p, 1-p)*2)	# *2 to make p-values two-sided
		return(res)
	})
	
	CIs <- t(CIs)
	return(CIs)
}