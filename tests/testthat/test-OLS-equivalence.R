context('RSA - equivalence with OLS regression')

test_that("Standard models (additive, IA, full, cubic) are equivalent to standard OLS", {
	
	tol <- .00001
  
	# create a new random data set
	data(motcon)
	r1 <- RSA(postVA ~ ePow * iPow, data=motcon, estimator="ML", se="standard", models=c("additive", "IA", "full", "SQD"), verbose=FALSE)
	r.cubic <- RSA(postVA ~ ePow * iPow, data=motcon, estimator="ML", se="standard", models=c("full", "cubic"), cubic=TRUE, verbose=FALSE)
	
	# compute analogous OLS regressions
	lm.additive <- lm(postVA ~ ePow + iPow, data=motcon)
	lm.IA <- lm(postVA ~ ePow*iPow, data=motcon)
	lm.full <- lm(postVA ~ ePow + iPow + I(ePow^2) + ePow:iPow + I(iPow^2), data=motcon)
	lm.cubic <- lm(postVA ~ ePow + iPow + I(ePow^2) + I(ePow^3) + ePow:iPow + I(iPow^2) + I(iPow^3) + I(ePow^2):iPow + I(iPow^2):ePow, data=motcon)
	
	RSA.additive.coef <- getPar(r1, type="coef", model="additive")
	rownames(RSA.additive.coef) <- RSA.additive.coef$label
	
	RSA.IA.coef <- getPar(r1, type="coef", model="IA")
	rownames(RSA.IA.coef) <- RSA.IA.coef$label
	
	RSA.full.coef <- getPar(r1, type="coef", model="full")
	rownames(RSA.full.coef) <- RSA.full.coef$label
	
	RSA.cubic.coef <- getPar(r.cubic, type="coef", model="cubic")
	rownames(RSA.cubic.coef) <- RSA.cubic.coef$label
	
	# are parameter estimates the same between RSA and lm?
	expect_equal(as.vector(coef(lm.additive)), RSA.additive.coef[paste0("b", 0:2), "est"], tolerance = tol)
	expect_equal(as.vector(coef(lm.IA)), RSA.IA.coef[c("b0", "b1", "b2", "b4"), "est"], tolerance = tol)
	expect_equal(as.vector(coef(lm.full)), RSA.full.coef[c("b0", "b1", "b2", "b3", "b5", "b4"), "est"], tolerance = tol)
	expect_equal(as.vector(coef(lm.cubic)), RSA.cubic.coef[c("b0", "b1", "b2", "b3", "b6", "b5", "b9", "b4", "b7", "b8"), "est"], tolerance = tol)
   	
})

