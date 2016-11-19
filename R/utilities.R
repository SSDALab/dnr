ungvectorize <- function(x,nvertex,gmode){
  if(gmode == "digraph"){
    out <- diag(nvertex)
    tmp <- c(out)
    tmp[tmp==0] <- x
    out <- matrix(tmp,nvertex)
    diag(out) <- 0
  }
  out
}

#'binaryPlot
#' @title binaryPlot
#' @param x matrix
#' @param ... title, xlabs, ylabs.
#' @description Plot for binary matrices, especially adjacency matrices.

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

