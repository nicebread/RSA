<!-- README.md is generated from README.Rmd. Please edit that file 
library(badger)
rmarkdown::render("README.rmd")
-->

[![Build
Status](https://travis-ci.org/nicebread/RSA.svg?branch=master)](https://travis-ci.org/nicebread/RSA)
[![Last-changedate](https://img.shields.io/badge/last%20change-2020--04--29-yellowgreen.svg)](/commits/master)
[![](https://www.r-pkg.org/badges/version/RSA?color=orange)](https://cran.r-project.org/package=RSA)
[![packageversion](https://img.shields.io/badge/Package%20version-0.10.0-orange.svg?style=flat-square)](commits/master)
![](http://cranlogs.r-pkg.org/badges/RSA)
[![](https://codecov.io/gh/nicebread/RSA/branch/master/graph/badge.svg)](https://codecov.io/gh/nicebread/RSA)

# RSA

An R package for Response Surface Analysis

## Installation

The stable CRAN version can be installed by:

    install.packages("RSA")

The current development version can be installed by:

    install.packages(c("devtools", "lavaan", "plyr", "ggplot2", "lattice", "tkrplot", "RColorBrewer", "rgl", "gridExtra", "aplpack", "fields", "qgraph", "AICcmodavg"))
    devtools::install_github("nicebread/RSA")

## Questions? Got to our Google mailing list

  - An email list for asking questions related to the RSA-package has
    been created at Google groups:
    <https://groups.google.com/forum/?fromgroups&hl=en#!forum/rsa-in-r>

## Demo script (with built-in data set)

    # if not already done: 
    # install the RSA package
    # (only has to be done once)
    # install.packages("RSA")
    
    # load the RSA package for
    # the active session
    library(RSA)
    
    # open the help page
    ?RSA
    
    ## Motive congruency example
    # load the built-in data set
    data(motcon)
    
    # Compute the RSA and save the result
    # into the new variable r1
    r1 <- RSA(postVA ~ ePow * iPow, data=motcon)
    
    # Show summary of the RSA
    print(r1)
    
    # Compare all models
    compare(r1, plot=TRUE)
    aictab(r1, plot=TRUE)
    
    # Show all RSA parameters with parametric
    # SEs, p values, and CIs
    getPar(r1, "coef", model="RR")
    
    # Standard CIs
    c1 <- confint(r1, model="RR")
    c1
    
    # Get bootstrapped confidence intervals
    # (5000 bootstrap replications), 
    # only from the RR model
    c2 <- confint(r1, model="RR", method="boot", R=5000)
    c2
    
    # Plot the final model
    plot(r1, model="RR", 
        xlab="Explicit power motive",
        ylab="Implicit power motive",
        zlab="Affective valence")
          
    
    ## Additional functions
    # contour plot
    plot(r1, model="RR", type="contour")
    
    # interactive, rotatable 3d plot
    plot(r1, model="RR", type="interactive")
    
    # open an interactive widget with control
    # sliders for regression weights
    demoRSA()

# Interactive plotting of polynomial surfaces

Not part of this package, but related: you can create polynomial surface
plots with this [interactive
app](http://shinyapps.org/showapp.php?app=https://tellmi.psy.lmu.de/felix/polySurface):

![](http://shinyapps.org/teaserpics/polynomialSurfaceExplorer.jpg)

# Maintainers

[Felix Schönbrodt](https://www.nicebread.de/) [Sarah
Humberg](https://www.uni-muenster.de/PsyIFP/AEBack/members/Sarah_Humberg.html)
