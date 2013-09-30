
# helpers.R

add.variables <- function(formula, df) {
	IV1 <- all.vars(formula)[2]
	IV2 <- all.vars(formula)[3]
	
	IV12 <- paste0(IV1, "2")
	IV22 <- paste0(IV2, "2")
	IV13 <- paste0(IV1, "3")
	IV23 <- paste0(IV2, "3")
	IV_IA <- paste0(IV1, "_", IV2)
	IV_IA2 <- paste0(IV1, "_", IV2, "2")
	IV_IA3 <- paste0(IV1, "2", "_", IV2)
	
	
	df[, IV12] <- df[, IV1]^2
	df[, IV22] <- df[, IV2]^2
	df[, IV_IA] <- df[, IV1]*df[, IV2]
	
	# three new variables for piecewise regression (test absolute difference score) - Edwards (2002) model
	df$W.JRE <- ifelse(df[, IV1] >= df[, IV2], 0, 1)
	df[, paste0("W.JRE_", IV1)] <- df$W.JRE*df[, IV1]
	df[, paste0("W.JRE_", IV2)] <- df$W.JRE*df[, IV2]
	
	# three new variables for piecewise regression (test absolute difference score) - new model Schoenbrodt 2012
	df$W <- ifelse(df[, IV1] >= df[, IV2], 1, -1)
	df$W[df[, IV1] == df[, IV2]] <- 0
	df[, paste0("W_", IV1)] <- df$W*df[, IV1]
	df[, paste0("W_", IV2)] <- df$W*df[, IV2]
	
	df$diff <- df[, IV2] - df[, IV1]
	df$sqdiff <- df$diff^2
	df$absdiff <- abs(df$diff)
	
	# cubic terms
	df[, IV13] <- df[, IV1]^3
	df[, IV_IA2] <- df[, IV1]*df[, IV2]^2
	df[, IV_IA3] <- df[, IV1]^2*df[, IV2]
	df[, IV23] <- df[, IV2]^3
	
	return(df)
}


sig2star <- function(val) {
	
	res <- val
	
	for (i in 1:length(val)) {
		res[i] <- ""
		if (is.na(val[i])) next();
		if (val[i] <= 0.1) res[i] <- "\U2020"
		if (val[i] <= 0.05) res[i] <- "*"
		if (val[i] <= 0.01) res[i] <- "**"
		if (val[i] <= 0.001) res[i] <- "***"
	}
	
	return(res)
}


# returns number of maximum free parameters of a regression model
getFreeParameters <- function(model) {
	VARS <- nrow(inspect(model, "free")$beta)	# number of variables
	df.max <- (VARS*(VARS+1))/2		# maximum df
	df.pred <- ((VARS-1)*(VARS))/2 + 1 # df bound in the predictors (i.e., (co)variances of the predictors & variance of DV)
	free.max <- df.max - df.pred	# maximum of free parameters
	return(free.max)
}




# helper function: takes a list of lavaan models (can include NULLs), and returns the usual anova object
anovaList <- function(modellist) {
	mods <- modellist[!sapply(modellist, function(x) is.null(x))]
	mods <- mods[!sapply(mods, function(x) !inspect(x, "converged"))]
	
	if (length(mods) == 0) {
		return(list(n.mods=0))
	}
	
    # put them in order (using df)
    DF <- sapply(mods, fitmeasures, "df")
    mods <- mods[order(DF, decreasing = FALSE)]
		
	pStr <- sapply(1:length(mods), function(x){ 
		if(x==1) {
			paste("mods[[",x,"]]",sep = "")
		} else {
			paste("force(mods[[",x,"]])",sep = "")
		}
	})
	pStr2 <- paste0("anova(", paste(pStr, collapse=", "), ")")
	
	a1 <- eval(parse(text = pStr2))
	
	if (length(mods) > 1) {
		rownames(a1) <- names(mods)
	}
	
	attr(a1, "n.mods") <- length(mods)
	return(list(ANOVA=a1, models=mods, n.mods=length(mods)))
}


# computes the coordinates of an arbitrary intersection of the surface,
# defined by a line on the X-Y plane (p0 = intercept, p1=slope)
getIntersect <- function(b0=0, x=0, y=0, x2=0, xy=0, y2=0, p0, p1, xlim=c(-2, 2), grid=21) {
	X <- seq(min(xlim), max(xlim), length.out=grid)
	Y <- p0 + p1*X
	n <- data.frame(X, Y)
	n2 <- add.variables(z~X+Y, n)
	n2$Z <- b0 + colSums(c(x, y, x2, y2, xy)*t(n2[, c(1:5)]))
	return(n2[, c("X", "Y", "Z")])
}


model <- function(x, model="full") x$models[[model]]

# transforms p-values to colors
pRamp <- function(p, sig=.05, borderline=.10, bias=.8) {
	# calculate bias that the color transition is at the borderline value
	bias2 <- .33/(borderline/(1 - sig))
	cR1 <- colorRamp(c("red", "red", "orange"), bias=bias, space="Lab")
	cR2 <- colorRamp(c("orange", "green", "green"), bias=bias2, space="Lab")
	
	p2 <- rep("#FFFFFF", length(p))
	if (length(p[p < sig])>0) {
		p2[p < sig] <- rgb(cR1(p[p < sig]/sig), maxColorValue=255)
	}
	if (length(p[p >= sig])>0) {
		p2[p >= sig] <- rgb(cR2((p[p >= sig] - sig) / (1 - sig)), maxColorValue=255)
	}
	return(p2)
}



# simple wrapper: formats a number in f.2 format
f2 <- function(x, digits=2, prepoint=0, skipZero=FALSE) {
	
	if (skipZero == TRUE) {zero <- "."} else {zero <- "0."}
	
	if (length(dim(x)) == 2) {
		apply(x, 2, function(x2) {gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x2) , fixed=TRUE)})
	} else {
		gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x) , fixed=TRUE)
	}
}



# helper function: find closest value in vector
f0 <- function (vec, target, unique = TRUE) {
    ret <- vec[sapply(target, function(x) which.min(abs(x - vec)))]
    if (unique) { ret <- unique(ret) }
    ret
}