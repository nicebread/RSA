
# helpers_cubic.R
# = functions needed for the interpretation of the cubic models


# ---------------------------------------------------------------------
# Range check function for strict and broad asymmetric congruence models
# ---------------------------------------------------------------------


# define function that, for a given point (x,y), computes the one-sided confidence interval of the model-prediction z at (x,y) according to the CA/RRCA model,

ci_pred <- function(obj, x, y, side, n, p, alpha, model){

  DV <- obj$DV
  
  # compute predicted outcome value at (x,y)
  z <- predictRSA(obj, x, y, model=model)

  # vector of predictor values
  r = c(1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3)
  
  # get covariances of the vector beta of estimated coefficients
  COV0 = vcov(obj$models[[model]])
  covBeta = COV0[c(paste0(DV,"~1"), "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9"),][,c(paste0(DV,"~1"), "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9")]
  
  # Compute standard error of the predicted value z
  var_z = t(r) %*% covBeta %*% r
  se_z = sqrt(var_z)
  
  # compute one-sided confidence interval of z. 
  if(side=="right"){
    lower_z = as.vector(-Inf)
    upper_z = as.vector(z + qt(1-alpha, df=n-p)*se_z)
  }
  
  if(side=="left"){
    lower_z = as.vector(z - qt(1-alpha, df=n-p)*se_z)
    upper_z = as.vector(Inf)
  }
  
  # output
  out <- list(z=z, lower_z=lower_z, upper_z=upper_z)
  return(out)
  }


# main function to identify data points behind E2 and test how many of them have outcome predictions that significantly differ from predictions for predictor combinations on E2 (that have the same level)
caRange <- function(object, alpha=0.05, verbose=TRUE, model="CA"){
      
  rsa = object
  
  # this function only makes sense for the models CA and RRCA
  if (!model %in% c("CA","RRCA")) stop("This function only makes sense for the models CA and RRCA.")

  # check whether model is contained in the RSA object and has converged
  if (!model %in% names(rsa$models[!unlist(lapply(rsa$models, is.null))])) stop("Please provide an RSA object that contains estimates for the asymmetric congruence model (CA or RRCA).")  
  
  if (inspect(rsa$models[[model]], "converged") == FALSE) {
    warning("The model has not converged!")
    return(NULL)
  }
  
  
  ## Extract information from RSA object
  
  # get variable names of outcome, x, and y
  DV <-  rsa$DV
  IV1 <- rsa$IV1
  IV2 <- rsa$IV2
  
  # get coefficient estimates of CA model
  c0 = as.numeric(coef(rsa$models[[model]])[paste0(DV,"~1")])
  c1 = as.numeric(coef(rsa$models[[model]])["b3"])
  c2 = as.numeric(coef(rsa$models[[model]])["b6"])
  
  # the present analysis makes no sense if the quadratic or cubic term coefficient is zero, because there is no second extremum line in these cases
  if (c1 == 0) stop("The quadratic term coefficient c1 is zero.")
  if (c2 == 0) stop("The cubic term coefficient c2 is zero.")
  
  # get raw data that was used to estimate the model
  data = rsa$data
  df = data[data$out==FALSE, ]
  
  # create variable with predicted outcome values
  df$z.pred <- NA
  for (k in 1:nrow(df)){df[k,"z.pred"] = predictRSA(rsa, df[k,IV1], df[k,IV2], model=model)}
  
                   
  
  # ---------------------------------------------------------------------
  # How many percent of the data points lie behind the second extremum line? 
  
  # define variable e = y - (x + 2c1/3c2) to be used in the next step
  df$e <- df[,IV2] - (df[,IV1] + 2*c1/(3*c2))
  
  if(c1*c2 < 0){
    
    # how many data points / what percentage of data points lie behind the second extremum line  <->   satisfy y < x + 2c1/3c2  <->  e = y - (x + 2c1/3c2) < 0
    df$behind <- df$e < 0
    df$behind[is.na(df$behind)] <- FALSE
    behind = sum(df$behind)
    p.behind = 100*(behind/nrow(df))
  }
  
  if(c1*c2 > 0){
    
    # how many data points / what percentage of data points lie behind the second extremum line  <->   satisfy y > x + 2c1/3c2  <->  e = y - (x + 2c1/3c2) > 0
    df$behind <- df$e > 0
    df$behind[is.na(df$behind)] <- FALSE
    behind = sum(df$behind)
    p.behind = 100*(behind/nrow(df))
  }  
  
  
  
  # ---------------------------------------------------------------------
  # How many percent of the data points behind the second extremum line have outcome predictions that differ significantly from the outcome prediction on the second extremum line?
  
  # In the following, determine the number of "bad" data points behind the second extremum line, separately for CA versus RRCA model.
  # Note that in fact, one could apply the procedure defined for the RRCA model (comparing the outcome prediction at each data point behind E2
  # to the outcome prediction at the point on E2 with the same level as the data point), as the CA model is a special case of the RR model. 
  # The results would be exactly the same as with the procedure defined for the CA model in the following.
  # Despite these redundancies, we strictly follow the descriptions in the manuscript here.
  
  if(model == "CA"){
    
    ## compute the confidence interval of the model prediction zr at the intersection point of E2 and the LOIC
    
    # Identify x value xr of intersection point of LOIC and second extremum line and compute the confidence interval around the predicted value zr at this point
    xr = -c1/(3*c2)
    yr = -xr

    # compute one-sided confidence interval of z. The direction of one-sidedness depends on whether the surface goes up versus down behind the second extremum line, which depends on the sign of c1.
    z_ci <- ci_pred(obj=rsa, x = xr, y = yr,
                    side=ifelse(c1<0, "right", "left"), 
                    n = nrow(df), p = 10, alpha=alpha,
                    model=model)
    
    zr <- z_ci$z

    ## Count points behind second extremum line whose outcome prediction differs significantly from zr
    df$lower_z <- z_ci$lower_z
    df$upper_z <- z_ci$upper_z
    df$badcase <- ifelse( (df$behind==T & df$z.pred < df$lower_z) | (df$behind==T & df$z.pred > df$upper_z) , TRUE, FALSE)
    
    bad = sum(df$badcase)
    p.bad = 100*(bad/nrow(df))
    
    response <- paste0(p.bad, "% of the data points lie behind the second extremum line AND have a significantly different outcome prediction than the points on the second extremum line.")
    
    # build results object for case of CA model
    res <- list(
      data.used = df,
      percentage.behind = p.behind,
      reversion_point_xy = paste0("(",round(xr,2),", ",round(yr,2), ")"),
      reversion_point_z = round(zr,2),
      reversion_point_z_lower = round(z_ci$lower_z,2),
      reversion_point_z_upper = round(z_ci$upper_z,2),
      reversion_point_z_ci = paste0("[", round(z_ci$lower_z,2),", ",round(z_ci$upper_z,2),"]"),
      which_bad_points = which(df$badcase==TRUE),
      how_many_bad_points = bad,
      percentage_bad_points = p.bad,
      response = response
      )
    
  } # end of computation for CA model
  
  

  if(model == "RRCA"){

    # for the RRCA model, the confidence interval to which the predicted value of each data point behind E2 must be compared is specific to the respective data point.

    # for each data point, compute the coordinates (xr,yr) of the point on E2 at the same level as the data point
    df$xr <- 0.5*(df[,IV1] + df[,IV2]) - (c1/(3*c2))
    df$yr <- 0.5*(df[,IV1] + df[,IV2]) + (c1/(3*c2))

    # create variables containing the zr value [ = predicted value at (xr,yr) ] and its confidence interval per data point
    df[,"zr"] <- apply(df, 1, function(data){ci_pred(obj=rsa, x=data["xr"], y=data["yr"],
                                                     side=ifelse(c1<0, "right", "left"),
                                                     n = nrow(df), p = 10, alpha=alpha,
                                                     model=model)$z})
    
    df[,"lower_zr"] <- apply(df, 1, function(data){ci_pred(obj=rsa, x=data["xr"], y=data["yr"],
                                                     side=ifelse(c1<0, "right", "left"),
                                                     n = nrow(df), p = 10, alpha=alpha,
                                                     model=model)$lower_z})
    
    df[,"upper_zr"] <- apply(df, 1, function(data){ci_pred(obj=rsa, x=data["xr"], y=data["yr"],
                                                     side=ifelse(c1<0, "right", "left"),
                                                     n = nrow(df), p = 10, alpha=alpha,
                                                     model=model)$upper_z})
    
    ## Count points behind second extremum line whose outcome prediction differs significantly from zr
    df$badcase <- ifelse( (df$behind==T & df$z.pred < df$lower_zr) | (df$behind==T & df$z.pred > df$upper_zr) , TRUE, FALSE)

    bad = sum(df$badcase)
    p.bad = 100*(bad/nrow(df))


    response <- paste0(p.bad, "% of the data points lie behind the second extremum line AND have a significantly different outcome prediction than the same-level point on the second extremum line.")
    
    # build results object for case of CA model
    res <- list(
      # data.used = df,
      percentage.behind = p.behind,
      which_bad_points = which(df$badcase==TRUE),
      how_many_bad_points = bad,
      percentage_bad_points = p.bad,
      response = response
    )


  } # end of computation for RRCA model



  # print results  
  if(verbose==T){
  print(response)
  }

  # return results
  return(res)
}



# ---------------------------------------------------------------------
# Range check function for strict and broad level-dependent congruence models
# ---------------------------------------------------------------------

# function to compute roots of the quadratic function ax^2 + bx + c
quad_roots <- function(a,b,c) {
  
  sqrt_term = b^2 - 4*a*c
  
  if (sqrt_term > 0) {
    pos_root <- ((-b) + sqrt(sqrt_term)) / (2*a)
    neg_root <- ((-b) - sqrt(sqrt_term)) / (2*a)
  } 
  
  if (sqrt_term == 0) {
    pos_root <- ((-b) + sqrt(sqrt_term)) / (2*a)
    neg_root <- NA
  } 
  
  if (sqrt_term < 0) {
    pos_root = neg_root = NA
  } 
  
  # output
  out <- c(pos_root,neg_root)
  return(out)
}


# function to compute the test statistic of the curvature of gk for a fixed value of k
pick_a_point <- function(c1,c3,cov,k){
  
  # compute curvature of gk
  eta = 4*c1 + 8*k*c3
  
  # variance of eta
  var = 16*cov["b3","b3"] + 64*cov["b3","b6"]*k + 64*cov["b6","b6"]*k^2
  
  # test statistic of eta
  teststat = eta/sqrt(var)
  
  # p value
  p = 2*pnorm(-abs(teststat))
  
  # output
  res = c(eta=eta, se=sqrt(var), teststat=teststat, p.value=p)
  return(res)
}



# main function for computing the regions of significance and testing their intersection with the data
clRange <- function(object, alpha=0.05, verbose=TRUE, model="CL"){
  
  rsa = object
  
  # this function only makes sense for the models CL and RRCL
  if (!model %in% c("CL","RRCL")) stop("This function only makes sense for the models CL and RRCL.")
  
  # check whether the model is contained in the RSA object and has converged
  if (!model %in% names(rsa$models[!unlist(lapply(rsa$models, is.null))])) stop("Please provide an RSA object that contains estimates for the level-dependent congruence model (CL or RRCL).")  
  
  if (inspect(rsa$models[[model]], "converged") == FALSE) {
    warning("The model has not converged!")
    return(NULL)
  }
  
  
  ## Extract information from RSA object
  
  # get variable names of outcome, x, and y
  DV <-  rsa$DV
  IV1 <- rsa$IV1
  IV2 <- rsa$IV2
  
  # get coefficient estimates of CL model
  c0 = as.numeric(coef(rsa$models[[model]])[paste0(DV,"~1")])
  c1 = as.numeric(coef(rsa$models[[model]])["b3"])
  c3 = as.numeric(coef(rsa$models[[model]])["b6"])
  
  # the present analysis makes no sense if the quadratic or cubic term coefficient is zero
  if (c1 == 0) stop("The quadratic term coefficient c1 is zero.")
  if (c3 == 0) stop("The coefficient c3 is zero.")
  
  # get covariances of the vector beta of estimated coefficients
  COV0 = vcov(rsa$models[[model]])

  # get raw data that was used to estimate the model
  data = rsa$data
  df = data[data$out==FALSE, ]
  

  # ---------------------------------------------------------------------
  ## Prepare result table
  
  tab = data.frame(region = NA, 
                   lower_bound = NA,
                   upper_bound = NA, 
                   data_points = NA,
                   percent_data = NA)
  tab <- tab[-1,]

  
  # ---------------------------------------------------------------------
  ## Determine the regions of negative/non/positive significance
  
  # compute k values where the curvature of gk changes from one significance state to another
  # = roots of the quadratic function ak^2 + bk + c, with a, b, c computed as follows
  z = qnorm(1-alpha/2)
  
  a = 64*c3^2 - 64*COV0["b6","b6"]*z^2
  b = 64*c1*c3 - 64*COV0["b3","b6"]*z^2
  c = 16*c1^2 - 16*COV0["b3","b3"]*z^2
  
  roots = quad_roots(a,b,c)
 
  
  # if there are two roots, the significance status changes at each of them. Determine the significance status of each interval.
  if (sum(is.na(roots)) == 0) {
    
    # name roots
    k1 <- lower_root <- min(roots)
    k2 <- higher_root <- max(roots)
    
    # check significance status in the interval ]-Inf, lower_root]
    row = nrow(tab)+1
    tab[row, c("lower_bound","upper_bound")] <- c(-Inf, lower_root)
    
    testvalue = lower_root - 10
    teststat = pick_a_point(c1, c3, cov=COV0, k=testvalue)["teststat"]
    
    if (teststat <= -z){ tab[row, "region"] <- "neg. sign."
    } else if (teststat >= z) { tab[row, "region"] <- "pos. sign."
    } else { tab[row, "region"] <- "nonsign."} 
    
    
    # check significance status in the interval [lower_root, higher_root]
    row = nrow(tab)+1
    tab[row, c("lower_bound","upper_bound")] <- c(lower_root, higher_root)
    
    testvalue = (lower_root+higher_root)/2
    teststat = pick_a_point(c1, c3, cov=COV0, k=testvalue)["teststat"]
    
    if (teststat <= -z){ tab[row, "region"] <- "neg. sign."
    } else if (teststat >= z) { tab[row, "region"] <- "pos. sign."
    } else { tab[row, "region"] <- "nonsign."} 
    
    
    # check significance status in the interval [higher_root, Inf[
    row = nrow(tab)+1
    tab[row, c("lower_bound","upper_bound")] <- c(higher_root, Inf)
    
    testvalue = higher_root + 10
    teststat = pick_a_point(c1, c3, cov=COV0, k=testvalue)["teststat"]
    
    if (teststat <= -z){ tab[row, "region"] <- "neg. sign."
    } else if (teststat >= z) { tab[row, "region"] <- "pos. sign."
    } else { tab[row, "region"] <- "nonsign."} 
    
  }
  

  # if there is one root, the significance status changes at it. Determine the significance status of each of the two intervals.
  if (sum(is.na(roots)) == 1) {
    
    # name roots
    k1 <- single_root <- roots[1]
    k2 <- NA

    # check significance status in the interval ]-Inf, single_root]
    row = nrow(tab)+1
    tab[row, c("lower_bound","upper_bound")] <- c(-Inf, single_root)
    
    testvalue = single_root - 10
    teststat = pick_a_point(c1, c3, cov=COV0, k=testvalue)["teststat"]
    
    if (teststat <= -z){ tab[row, "region"] <- "neg. sign."
    } else if (teststat >= z) { tab[row, "region"] <- "pos. sign."
    } else { tab[row, "region"] <- "nonsign."} 
    
    # check significance status in the interval [single_root, Inf[
    row = nrow(tab)+1
    tab[row, c("lower_bound","upper_bound")] <- c(single_root, Inf)
    
    testvalue = single_root + 10
    teststat = pick_a_point(c1, c3, cov=COV0, k=testvalue)["teststat"]
    
    if (teststat <= -z){ tab[row, "region"] <- "neg. sign."
    } else if (teststat >= z) { tab[row, "region"] <- "pos. sign."
    } else { tab[row, "region"] <- "nonsign."} 
    
  }
  
 
  # if there are no roots, the significance status is the same for all k. Determine which one it is.
  if (sum(is.na(roots)) == 2) {
    
    k1 <- k2 <- NA
    
    # check significance status for the whole surface
    row = nrow(tab)+1
    tab[row, c("lower_bound","upper_bound")] <- c(-Inf, Inf)
    
    teststat = pick_a_point(c1, c3, cov=COV0, k=0)["teststat"]
    
    if (teststat <= -z){ tab[row, "region"] <- "neg. sign."
    } else if (teststat >= z) { tab[row, "region"] <- "pos. sign."
    } else { tab[row, "region"] <- "nonsign."} 
    
  } 
  
 

  ## check how many data points lie in each region
  
  # for each data point (xi,yi), determine k so that (xi,yi) lies on gk
  df$k = (df[,IV1] + df[,IV2])/2

  # how many data points lie within each region? 
  for (i in 1:nrow(tab)){
    tab[i, "data_points"] <- length(which(df$k > tab[i,"lower_bound"] & df$k < tab[i,"upper_bound"]))
  }  
  
  tab$percent_data = 100*(tab$data_points/nrow(df))

  
  # ---------------------------------------------------------------------
  # build results object
  res <- list(
    data.used = df,
    k1 = k1,
    k2 = k2,
    regions = tab
  )
  
  # return results
  return(res)
}
