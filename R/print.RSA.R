#' @export
#' @method print RSA
print.RSA <- function(x, ..., model=NULL, digits=3) {
	summary.RSA(object=x, ..., model=model, digits=digits)
}

#' @export
#' @method summary RSA
summary.RSA <- function(object, ..., model=NULL, digits=3) {
	x <- object
	
	# No model provided? Default to the respective global model
	if (is.null(model)) {
	  
	  # get all model names that have been computed
	  available.models <- names(x$models)[!sapply(x$models, is.null)]
	  
	  if ("cubic" %in% available.models) {
	    model <- "cubic"
	  } else if ("full" %in% available.models) {
	    model <- "full"
	  } 
	  
	   message(paste0("No model has been specified - showing results for the full ", ifelse(model=="cubic", "third-order", "second-order")," polynomial model (<", model, ">)."))
	   
	}
	
	# is the model a cubic model?
	is.cubicmodel <- model %in% c("cubic","CA","RRCA","CL","RRCL")
	
	# Print model summary, also show package version
	if(!exists("meta") || is.null(meta)) meta <- packageDescription("RSA")
	cat(sprintf("RSA output (package version %s)", meta$Version))
	cat("\n===========================================\n\n")
	
	with(x, {
		
		eff <- getPar(x, model=model, standardized=TRUE)
		eff <- eff[order(eff$label), ]
		
		if (!model %in% c("absunc", "absdiff") & !is.cubicmodel) {
			ST <- RSA.ST(x, model=model)
		} else {
			ST <- NULL
		}
		
		
	## Step 1: Examine amount of discrepancy
	#--------------------------------------------------
	# When using polynomial regression analyses for the investigation of congruence effects, it is important to inspect how many participants would be considered to have discrepancies between the two predictors so that you have an idea of the base rate of discrep- ancies in your sample.
	
	cat("Are there discrepancies in the predictors (with respect to numerical congruence)?\n----------------------------\n")
	D <- (data[, IV2] - data[, IV1])/sd(c(data[, IV1], data[, IV2]), na.rm=TRUE)
	Congruence <- cut(D, breaks=c(-Inf, -.5, .5, Inf), labels=c(paste0(IV2, " < ", IV1), "congruent", paste0(IV2, " > ", IV1)))

	congruenceTab <- paste0(as.vector(table(Congruence)), "%")
	names(congruenceTab) <- names(table(Congruence))
	print(congruenceTab)
	
	
	cat(paste0("\nIs the corresponding global model (full ", ifelse(is.cubicmodel, "third-order", "second-order")," polynomial model) significant?\n----------------------------\n"))
	
	globalmodel <- ifelse(is.cubicmodel, "cubic", "full")
	
	R2.global <- R2difftest(x, unrestricted=globalmodel)
	cat(paste0("Test on ", ifelse(is.cv, "overall ", ""), "model significance: R^2 = ", round(R2.global$delta.R2, 3), ", ", p0(R2.global$R2.p), "\n"))
	
	if (is.cv) {
	  R2.global.cv <- R2difftest(x, unrestricted=globalmodel, restricted="null")
	  cat(paste0("R^2 increment above baseline model (with control variables): Delta_R^2 = ",  round(R2.global.cv$delta.R2, 3), ", ", p0(R2.global.cv$R2.p), "\n"))
	}
	

	if ( model != globalmodel ) {
		cat(paste0("\nIs the selected model <", model, "> significant?\n----------------------------\n"))
	  
	  R2.model <- R2difftest(x, unrestricted=model)
	  cat(paste0("Test on ", ifelse(is.cv, "overall ", ""), "model significance: R^2 = ", round(R2.model$delta.R2, 3), ", ", p0(R2.model$R2.p), "\n"))
	  
	  if (is.cv) {
	    R2.model.cv <- R2difftest(x, unrestricted=model, restricted="null")
	    cat(paste0("R^2 increment above baseline model (with control variables): Delta_R^2 = ",  round(R2.model.cv$delta.R2, 3), ", ", p0(R2.model.cv$R2.p), "\n"))
	  }
	
	}

	
	cat(paste0("\n\nNumber of observations: n = ", nobs(x$models[[model]]), "\n----------------------------\n"))

	cat(paste0("\n\nRegression coefficients for model <", model, ">\n----------------------------\n"))
	
	coef.sel <- paste0("b", 0:5)
	if (is.cubicmodel) { coef.sel <- c(coef.sel, paste0("b", 6:9)) }
	if (is.cv) { coef.sel <- c(coef.sel, paste0("cv", 1:length(control.variables))) }
	
	RC <- eff[eff$label %in% coef.sel, c(1:3, 6:7)]
	RC[, 2:5] <- round(RC[, 2:5], digits)
	RC$beta <- round(eff[eff$label %in% coef.sel, "std.all"], digits)
	RC$pvalue <- p(eff[eff$label %in% coef.sel, "pvalue"])
	RC$sig <- p2star(eff[eff$label %in% coef.sel, "pvalue"])
	print(RC)	
	

	if (!model %in% c("onlyx", "onlyy") & !is.cubicmodel) {
		cat(paste0("\n\n\nSurface tests (a1 to a5) for model <", model, ">\n----------------------------\n"))
		as <- eff[eff$label %in% paste0("a", 1:5), c(1:3, 6:7)]
		as[, 2:5] <- round(as[, 2:5], digits)
		as$pvalue <- p(eff[eff$label %in% paste0("a", 1:5), "pvalue"])
		as$sig <- p2star(eff[eff$label %in% paste0("a", 1:5), "pvalue"])
		rownames(as) <- NULL
		print(as)
	
		# print interpretations:
		cat(paste0("\na1: Linear additive effect on line of congruence? ", ifelse(eff[eff$label %in% "a1", "pvalue"] <= .05, "YES", "NO"), "\n"))

		cat(paste0("a2: Is there curvature on the line of congruence? ", ifelse(eff[eff$label %in% "a2", "pvalue"] <= .05, 
		"YES", "NO"), "\n"))

		cat(paste0("a3: Is the ridge shifted away from the LOC? ", ifelse(eff[eff$label %in% "a3", "pvalue"] <= .05, "YES", "NO"), "\n"))	
	
		cat(paste0("a4: Is there curvature on the line of incongruence? ", ifelse(eff[eff$label %in% "a4", "pvalue"] <= .05, "YES", "NO"), "\n"))

	
		PA <- eff[eff$label %in% c("p10", "p11", "p20", "p21"), c(1:3, 6:7)]
		
		if (nrow(eff[eff$label %in% c("X0", "Y0"), ]) > 0) {
			
			cat(paste0("\n\nLocation of stationary point for model <", model, ">\n----------------------------\n"))
			cat(paste0(IV1, " = ", round(eff[eff$label %in% "X0", "est"], 3), "; ", IV2, " = ", round(eff[eff$label %in% "Y0", "est"], 3), "; predicted ", DV, " = ", round(predictRSA(x, eff[eff$label %in% "X0", "est"], eff[eff$label %in% "Y0", "est"], model=model), 3), "\n\n"))
		}
		
		if (nrow(PA) > 0) {
		
			cat(paste0("\nPrincipal axes for model <", model, ">\n----------------------------\n"))
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
	}
	
	
	# ---------------------------------------------------------------------
	# Tests for fit patterns
	
#	if (all(c("full", "weak", "strong") %in% names(models))) {
# 		cat("\n\n\nTests for weak and strong fit patterns:\n----------------------------\n")
#
# 		eff.full <- getPar(x, model="full", standardized=TRUE)
# 		coef.sel <- paste0("b", c(3, 5))
# 		RC.full <- eff.full[eff.full$label %in% coef.sel, c(1:3, 6:7)]
# 		RC.full[, 2:5] <- round(RC.full[, 2:5], digits)
# 		RC.full$beta <- round(eff.full[eff.full$label %in% coef.sel, "std.all"], digits)
# 		RC.full$pvalue <- p(eff.full[eff.full$label %in% coef.sel, "pvalue"])
# 		RC.full$sig <- p2star(eff.full[eff.full$label %in% coef.sel, "pvalue"])
# 		print(RC.full)
#
# 		cat("\n\nWeak fit pattern:\n")
# 		aictab1 <- aictab(x, cand.set=c("full", "weak"))
# 		aicdiff.weak <- aictab1$AICc[aictab1$Modnames == "weak"] - aictab1$AICc[aictab1$Modnames == "full"]
#
# 		print(aictab1)
#
# 		is.weak <- aicdiff.weak <= 0.001
# 		cat(paste0("The data ", ifelse(is.weak==TRUE, "DO", "DO NOT"), " support a weak fit pattern."))
#
# 		if (is.weak == TRUE) {
# 			cat("\n\nStrong fit pattern:\n")
# 			aictab2 <- aictab(x, cand.set=c("weak", "strong"))
# 			ctab2 <- compare2(x, m1="weak", m2="strong", verbose=FALSE)
# 			print(round(ctab2[, 1:13], 3))
# 			aicdiff.weak <- aictab1$AICc[aictab1$Modnames == "weak"] - aictab1$AICc[aictab1$Modnames == "full"]
# 			cat(paste0("The data ", ifelse(ctab2["strong", "Pr(>Chisq)"] > .05, "DO", "DO NOT"), " support a strong fit pattern."))
# 		}
# 	}
	
	
	
	})
}

