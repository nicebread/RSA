#' @S3method print RSA
print.RSA <- function(x, ..., model="full", digits=3) {
	with(x, {
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
	cat(paste0("Test on model significance: R2 = ", round(r2.model, 3), " (p = ", round(p.model, 3), ")\n"))
		
	cat(paste0("\n\nRegression coefficients for model <", model, ">\n----------------------------\n"))
	RC <- getPar(x, model=model, standardized=TRUE)[paste0("b", 1:5), 1:8]
	rownames(RC) <- NULL
	colnames(RC)[8] <- "beta"
	RC[, 2:8] <- round(RC[, 2:8], digits)
	RC$sig <- sig2star(RC$pvalue)
	print(RC)
	
	
	cat(paste0("\n\nEigenvalues and shape of surface for model <", model, ">\n----------------------------\n"))
	
	ST <- RSA.ST(x, model=model)
	cat("Eigenvalues:\n")
	cat("lambda_1: ", round(ST$l[1], 3), "\n")
	cat("lambda_2: ", round(ST$l[2], 3), "\n")
	shape <- "undefined"
	if (all(ST$l < 0)) shape <- "Dome shaped (stationary point = maximum response)"
	if (all(ST$l > 0)) shape <- "Bowl shaped (stationary point = minimum response)"
	if ((ST$l[1] > 0 & ST$l[2] < 0) | (ST$l[2] > 0 & ST$l[1] < 0)) shape <- "Saddle shaped"
	cat(paste0("--> ", shape))
	
	
		
	cat(paste0("\n\n\nSurface tests (a1 to a4) for model <", model, ">\n----------------------------\n"))
	print(round(ST$SP, 3))
	# print interpretations:
	cat(paste0("\na1: Linear additive effect on line of congruence? ", ifelse(ST$SP$p.value[1] <= .05, "YES", "NO"), "\n"))

	cat(paste0("a2: Is there curvature on the line of congruence? ", ifelse(ST$SP$p.value[2] <= .05, 
	"YES", "NO"), "\n"))

	cat(paste0("a3: Is the ridge shifted away from the LOC? ", ifelse(ST$SP$p.value[3] <= .05, "YES", "NO"), "\n"))	
	
	cat(paste0("a4: Is there a general effect of incongruence? ", ifelse(ST$SP$p.value[4] <= .05, "YES", "NO"), "\n"))

	
	
	cat(paste0("\n\nLocation of stationary point (minimum, maximum, or saddle point response) for model <", model, ">\n----------------------------\n"))
	cat(paste0(IV1, " = ", round(ST$X0, 3), "; ", IV2, " = ", round(ST$Y0, 3), "; predicted ", DV, " = ", round(ST$Z0, 3), "\n\n"))
	
	cat(paste0("\nPrincipal axes for model <", model, ">\n----------------------------\n"))
	C <- coef(models[["full"]], "all")
	cat(paste("First principal axis: p11 =", round(C["p11"], 3), "; p10 =", round(C["p10"], 3), "\n"))
	#cat("  --> Deviation of p11 from 1 = Rotation of the first principal axis from X=Y (line of congruence)\n")
	cat("  --> Lateral shift from LOC at point (0; 0): C1 = ", round((-C["p10"])/(C["p11"] + 1), 3), "\n")
	cat(paste("Second principal axis: p21 =", round(C["p21"], 3),  "; p20 =", round(C["p20"], 3), "\n"))
	#cat("  --> Deviation of p21 from 1 = Rotation of the second principal axis from Y= -X (line of congruence)\n")
	cat("  --> Lateral shift from LOC at point (0; 0): C2 = ", round((-C["p20"])/(C["p21"] + 1), 3), "\n")
	
	})
}

