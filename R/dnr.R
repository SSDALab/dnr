#' @keywords internal
#' @import ergm
#' @import network
#' @importFrom graphics axis grid image
#' @importFrom stats as.formula binomial glm na.omit rbinom setNames simulate terms
#' @importFrom utils combn
##' @title dnr: A package for simulating dynamic networks using ERGM family models.

##' @description
##' This package provides functions for fitting lagged exponential family models on dynamic network data, simulation from the models and model diagnostics. \cr

##' @section note
##' This package was developed with help from ARO YIP award \#W911NF-14-1-0577. \cr

##' @section references
##' Abhirup Mallik and Zack W. Almquist (2017). "Stable Multiple Time Step Simulation/Prediction from Lagged Dynamic Network Regression Models." Working paper. University of Minnesota. \cr

##' Abhirup Mallik and Zack W. Almquist (2017). "An R Package for Dynamic Network Regression." Working paper. University of Minnesota. \cr

##' Zack W. Almquist and Carter T. Butts (forthcoming). "Dynamic Network Analysis with Missing Data: Theory and Methods." Statistica Sinica. doi:10.5705/ss.202016.0108. \cr

##' Zack W. Almquist and Carter T. Butts. (2013). "Dynamic Network Logistic Regression: A Logistic Choice Analysis of Inter- and Intra-group Blog Citation Dynamics in the 2004 US Presidential Election." Political Analysis, 21(4), 430-448. \cr

##' Zack W. Almquist and Carter T. Butts. (2014). "Bayesian Analysis of Dynamic Network Regression with Joint Edge/Vertex Dynamics." In Bayesian Inference in the Social Sciences. Ed. by I. Jeliazkov and X.-S. Yang. Hoboken, New Jersey: John Wiley & Sons. \cr

##' @docType package
##' @name dnr
NULL