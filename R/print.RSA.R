#' @S3method print RSA
print.RSA <- function(x, ..., model="full", digits=3) {
	summary.RSA(object=x, ..., model=model, digits=digits)
}

#' @S3method summary RSA
summary.RSA <- function(object, ..., model="full", digits=3) {
	x <- object
	with(x, {
		
		eff <- getPar(x, model=model, standardized=TRUE)
		
		if (!model %in% c("cubic", "absunc", "absdiff")) {
			ST <- RSA.ST(x, model=model)
		} else {
			ST <- NULL
		}
		
		
	## Step 1: Examine amount of discrepancy
	#--------------------------------------------------
	# Before conducting the polynomial regression analyses, it is important to inspect how many participants would be considered to have discrepancies between the two predictors so that you have an idea of the base rate of discrep- ancies in your sample.
	
	cat("Are there discrepancies in the predictors (with respect to numerical congruence)?\n----------------------------\n")
	D <- data[, IV2] - data[, IV1]
	Congruence <- cut(D, breaks=c(-Inf, -.5, .5, Inf), labels=c(paste0(IV2, " < ", IV1), "Congruence", paste0(IV2, " > ", IV1)))
	print(round(prop.table(table(Congruence)), 3)*100)
	
	
	cat("\nIs the full polynomial model significant?\n----------------------------\n")
	# --> is R2 significant?
	r2.model <- summary(LM)$r.squared
	
	F <- summary(LM)$fstatistic
	p.model <- 1-pf(F[1], F[2], F[3])
	cat(paste0("Test on model significance: R2 = ", round(r2.model, 3), ", ", p(p.model), "\n"))
		

	cat(paste0("\n\nRegression coefficients for model <", model, ">\n----------------------------\n"))
	if (model != "cubic") {
		coef.sel <- paste0("b", 1:5)
	} else {
		coef.sel <- paste0("b", c(1:5, 9:12))
	}
	RC <- eff[eff$label %in% coef.sel, c(1:3, 6:7)]
	RC[, 2:5] <- round(RC[, 2:5], digits)
	RC$beta <- round(eff[eff$label %in% coef.sel, "std.lv"], digits)
	RC$pvalue <- p(eff[eff$label %in% coef.sel, "pvalue"])
	RC$sig <- p2star(eff[eff$label %in% coef.sel, "pvalue"])
	print(RC)	
	
	if (!model %in% c("onlyx", "onlyy", "cubic")) {
		cat(paste0("\n\nEigenvalues and shape of surface for model <", model, ">\n----------------------------\n"))
		cat("Eigenvalues:\n")
		EV <- eff[eff$label %in% c("l1", "l2"), c("label", "est", "se", "ci.lower", "ci.upper")]
		EV[, 2:5] <- round(EV[, 2:5], digits)
		EV$pvalue <- p(eff[eff$label %in% c("l1", "l2"), "pvalue"])
		EV$sig <- p2star(eff[eff$label %in% c("l1", "l2"), "pvalue"])
		rownames(EV) <- NULL
		EV$label[1:2] <- c("lambda_1", "lambda_2")
		print(EV)
		shape <- "undefined"
		if (all(EV$est < 0)) shape <- "Dome shaped (stationary point = maximum response)"
		if (all(EV$est > 0)) shape <- "Bowl shaped (stationary point = minimum response)"
		if ((EV$est[1] > 0 & EV$est[2] < 0) | (EV$est[2] > 0 & EV$est[1] < 0)) shape <- "Saddle shaped"
		cat(paste0("--> ", shape))
	}
	
		
	if (!model %in% c("onlyx", "onlyy", "cubic")) {
		cat(paste0("\n\n\nSurface tests (a1 to a4) for model <", model, ">\n----------------------------\n"))
		as <- eff[eff$label %in% paste0("a", 1:4), c(1:3, 6:7)]
		as[, 2:5] <- round(as[, 2:5], digits)
		as$pvalue <- p(eff[eff$label %in% paste0("a", 1:4), "pvalue"])
		as$sig <- p2star(eff[eff$label %in% paste0("a", 1:4), "pvalue"])
		rownames(as) <- NULL
		print(as)
	
		# print interpretations:
		cat(paste0("\na1: Linear additive effect on line of congruence? ", ifelse(eff[eff$label %in% "a1", "pvalue"] <= .05, "YES", "NO"), "\n"))

		cat(paste0("a2: Is there curvature on the line of congruence? ", ifelse(eff[eff$label %in% "a2", "pvalue"] <= .05, 
		"YES", "NO"), "\n"))

		cat(paste0("a3: Is the ridge shifted away from the LOC? ", ifelse(eff[eff$label %in% "a3", "pvalue"] <= .05, "YES", "NO"), "\n"))	
	
		cat(paste0("a4: Is there a general effect of incongruence? ", ifelse(eff[eff$label %in% "a4", "pvalue"] <= .05, "YES", "NO"), "\n"))

	
		cat(paste0("\n\nLocation of stationary point: ", shape, " for model <", model, ">\n----------------------------\n"))
		cat(paste0(IV1, " = ", round(eff[eff$label %in% "X0", "est"], 3), "; ", IV2, " = ", round(eff[eff$label %in% "Y0", "est"], 3), "; predicted ", DV, " = ", round(ST$Z0, 3), "\n\n"))
	
		
		
		cat(paste0("\nPrincipal axes for model <", model, ">\n----------------------------\n"))
		
		PA <- eff[eff$label %in% c("p10", "p11", "p20", "p21"), c(1:3, 6:7)]
		PA[, 2:5] <- round(PA[, 2:5], digits)
		PA$pvalue <- p(eff[eff$label %in% c("p10", "p11", "p20", "p21"), "pvalue"])
		PA$sig <- p2star(eff[eff$label %in% c("p10", "p11", "p20", "p21"), "pvalue"])
		rownames(PA) <- NULL
		rownames(PA)[PA$label == "p10"] <- "Intercept of 1. PA"
		rownames(PA)[PA$label == "p11"] <- "Slope of 1. PA"
		rownames(PA)[PA$label == "p20"] <- "Intercept of 2. PA"
		rownames(PA)[PA$label == "p21"] <- "Slope of 2. PA"
		print(PA)
		
		C <- coef(models[["full"]], "all")
		
		cat("  --> Lateral shift of first PA from LOC at point (0; 0): C1 = ", round((-C["p10"])/(C["p11"] + 1), 3), "\n")
		cat("  --> Lateral shift of second PA from LOC at point (0; 0): C2 = ", round((-C["p20"])/(C["p21"] + 1), 3), "\n")
	}
	
	})
}

