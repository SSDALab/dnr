##' General purpose regression engine for the methods bayesglm, glm and glmnet
##' @param XYdata matrix with X and Y columns. First column is named as y, other columns are X.
##' @param method string among ("glm", "glmnet", "bayesglm").
##' @param regIntercept Logical. Should intercept be included in the model? 
##' @param lambda for method "glmnet".
##' @param alpha for "glmnet"
##' @return list with elements: coef, se, lambda, fit (Coefficients, SE, lambda, if used, fit object.)
##' @author Abhirup

regEngine <- function(XYdata,
                      method = "bayesglm",
                      regIntercept = FALSE,
                      lambda = NA,
                      alpha = 1){
    if(method == "glmnet"){
        if(is.na(lambda)){
            mSelect <- glmnet::cv.glmnet(
                                   data.matrix(XYdata[, -1]),
                                   as.vector(XYdata[, 1]),
                                   family = "binomial",
                                   alpha = alpha)
            lambda <- mSelect$lambda.min
        }
        mFit <- glmnet::glmnet(data.matrix(XYdata[, -1]),
                               as.vector(XYdata[, 1]),
                               family = "binomial",
                               alpha = alpha,
                               lambda = lambda,
                               intercept = regIntercept)
        mCoef <- setNames(as.vector(mFit$beta),
                          colnames(XYdata[, -1]))
        mSE <- NA
        out <- list(coef = mCoef,
                    se = mSE,
                    lambda = lambda,
                    fit = mFit)
    } else if (method == "glm") {
        if(regIntercept){
            mFit <- glm(y ~ ., data = XYdata,
                        family = binomial(link = "logit"))
            mSummery <- (summary(mFit))
            mCoef <- mFit$coefficients[-1]
            mSE <- numeric(length(mCoef))
            mSE[is.na(mCoef)] <- 0
            mSE[!is.na(mCoef)] <- mSummery$coefficients[-1 , 2]
            names(mSE) <- colnames(XYdata[, -1])
        } else{
            mFit <- glm(y ~ . -1, data = XYdata,
                        family = binomial(link = "logit"))
                        mSummery <- (summary(mFit))
            mCoef <- mFit$coefficients
            mSE <- numeric(length(mCoef))
            mSE[is.na(mCoef)] <- 0
            mSE[!is.na(mCoef)] <- mSummery$coefficients[ , 2]
            names(mSE) <- colnames(XYdata[, -1])
        }
        out <- list(coef = mCoef,
                    se = mSE,
                    lambda = lambda,
                    fit = mFit)
    } else if(method == "bayesglm"){
        if(regIntercept) {
            mFit <- arm::bayesglm(y ~., data = XYdata,
                                  family = binomial(link = "logit"))
            mSummery <- (summary(mFit))
            mCoef <- mFit$coefficients[-1]
            mSE <- numeric(length(mCoef))
            mSE[is.na(mCoef)] <- 0
            mSE[!is.na(mCoef)] <- mSummery$coefficients[-1 , 2]
            names(mSE) <- colnames(XYdata[, -1])
        } else {
            mFit <- arm::bayesglm(y ~ . -1, data = XYdata,
                                  family = binomial(link = "logit"))
            mSummery <- (summary(mFit))
            mCoef <- mFit$coefficients
            mSE <- numeric(length(mCoef))
            mSE[is.na(mCoef)] <- 0
            mSE[!is.na(mCoef)] <- mSummery$coefficients[ , 2]
            names(mSE) <- colnames(XYdata[, -1])
        }
        out <- list(coef = mCoef,
                    se = mSE,
                    lambda = lambda,
                    fit = mFit)
    }
    return(out)
}




ungvectorize <- function(x,nvertex,gmode){
  if(gmode == "digraph"){
    out <- diag(nvertex)
    tmp <- c(out)
    tmp[tmp==0] <- x
    out <- matrix(tmp,nvertex)
    diag(out) <- 0
  } else{
    out <- diag(nvertex)
    tmp <- c(out)
    tmpid <- c(lower.tri(out))
    tmp[tmpid] <- x
    out <- matrix(tmp, nvertex)
    out[upper.tri(out)] <- t(out)[upper.tri(out)]
    diag(out) <- 0
  }
  out
}


##'binaryPlot
##' @title binaryPlot
##' @param x matrix
##' @param axlabs Binary, should the axis labels be shown.
##' @param ... title, xlabs, ylabs.
##' @description Plot for binary matrices, especially adjacency matrices.
##' @export

binaryPlot <- function(x, axlabs = TRUE,...){
    xmin <- min(x)
    xmax <- max(x)
    yLabels <- rownames(x)
    xLabels <- colnames(x)
    title <-c()
    if( length(list(...)) ){
        Lst <- list(...)
        if( !is.null(Lst$zlim) ){
            min <- Lst$zlim[1]
            max <- Lst$zlim[2]
        }
        if( !is.null(Lst$yLabels) ){
            yLabels <- c(Lst$yLabels)
        }
        if( !is.null(Lst$xLabels) ){
            xLabels <- c(Lst$xLabels)
        }
        if( !is.null(Lst$title) ){
            title <- Lst$title
        }
    }
                                        # check for null values
    if( is.null(xLabels) ){
        xLabels <- c(1:ncol(x))
    }
    if( is.null(yLabels) ){
        yLabels <- c(1:nrow(x))
    }
                                        # Reverse x
    reverse <- nrow(x) : 1
    yLabels <- yLabels[reverse]
    x <- x[reverse,]
    image(1:length(xLabels), 1:length(yLabels), t(x),
          col=c("#cccccc", "#525252"), xlab="",
          ylab="", axes=FALSE, zlim=c(xmin,xmax))
    if( !is.null(title) ){
        title(main=title)
    }
    if(axlabs){
        axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
        axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
             cex.axis=0.7)
    }
    grid(NULL, NULL, col = "#969696", lty = 6)
}


ilogit <- function(x) 1/(1+exp(-x))
thresh <- function(x, hi = 0.95, lo = 0.05) {
    x <- ifelse(x>hi, hi, x)
    x <- ifelse(x<lo, lo, x)
}

## supporting functions:
## network size (that allows NA)
network.size.1 <- function(x){
  if(!network::is.network(x)){
    if(is.na(x)) return(0)
  } else return(network::network.size(x))
}

## Remove networks that are of size 0.
rmNAnetworks <- function(netlist){
  netlens <- unlist(lapply(netlist, network.size.1))
  toremove <- which(netlens < 1)
  if(length(toremove) > 0) {
    return(netlist[-toremove])
  } else {
    return(netlist)
  }
}

network.vertex.names.1 <- function(x) {
  if(!network::is.network(x)){
    if(is.na(x)) return (NA)
  } else {
    return (network::network.vertex.names(x))
  }
}
