##' Parameter estimation for Vertex model only for a list of dynamic networks.
##' @param InputNetwork Input network list.
##' @param VertexStatsvec Binary vector of size 9, indicating vertex model.
##' @param maxLag maximum lag.
##' @param VertexModelGroup Group term for vertex model.
##' @param VertexLag Binary vector of size maxLag, indicating Lag terms in the model.
##' @param VertexAttLag Vertex group term lag vector.
##' @param VertexLagMatrix Binary matrix indicating lagged vertex statistics in
##'     the model.
##' @param regMethod one of "glm", "glmnet", "bayesglm"
##' @return List of 3 elements:\cr
##' VertexFit: Output from regEngine. \cr
##' VertexStats: Subsetted vertex stats matrix. \cr
##' VertexStatsFull: Full matrix of vertex stats.
##' @author Abhirup
##' @export
##' @examples
##' nvertexstats <- 9
##' InputNetwork <- beach
##' maxLag <- 3
##' VertexStatsvec <- rep(1, nvertexstats)
##' VertexLag <- rep(1, maxLag)
##' regMethod <- "bayesglm"
##' VertexModelGroup <- "regular"
##' VertexLagMatrix <- matrix(0, maxLag, nvertexstats)
##' VertexLagMatrix[, c(4, 7)] <- 1
##' VertexLagMatrix[c(2,3),7] <- 0
##' Vout1 <- paramVertexOnlyGroup(InputNetwork = beach,
##'                           maxLag = maxLag,
##'                           VertexStatsvec = VertexStatsvec,
##'                           VertexModelGroup = VertexModelGroup,
##'                           VertexLag = VertexLag,
##'                           VertexLagMatrix = VertexLagMatrix)
##' summary(Vout1$VertexFit$fit)

paramVertexOnlyGroup <- function(InputNetwork,
                                 VertexStatsvec = rep(1, nvertexstats),
                                 maxLag,
                                 VertexModelGroup = NA,
                                 VertexLag = rep(1, maxLag),
                                 VertexAttLag = rep(1, maxLag),
                                 VertexLagMatrix = matrix(1, maxLag,
                                                          length(VertexStatsvec)),
                                 regMethod = "bayesglm"){


    InputNetwork <- rmNAnetworks(InputNetwork)
    netlength <- length(InputNetwork)
    Vunion <- unique(unlist(lapply(InputNetwork, network.vertex.names.1)))
    Vunion <- na.omit(Vunion)
    nv <- length(Vunion)

    ## is the network directed?
    ## we are using 1st network as representative network hoping these
    ## attributes wont change.
    isdirected <- network::get.network.attribute(InputNetwork[[1]], "directed")
    if(isdirected) {
        gmode <- "digraph"
    } else gmode <- "graph"

    net1.vstats <- vertexstats(InputNetwork[[1]], gmode = gmode)
    nvertexstats <- ncol(net1.vstats)
    ## Change this number when making modification in
    ## vertexStatsNew.R
    Verstats <- matrix(0, nrow = length(Vunion), ncol = nvertexstats)
    rownames(Verstats) <- Vunion
    Verstats[match(rownames(net1.vstats),
                   rownames(Verstats)),] <-
        net1.vstats
    y <- NULL
    x.Nets <- NULL
    x.vstats <- NULL
    x.vatts <- NULL
###############################
    for(i in seq_len(netlength - maxLag)){
        vstats.current <- NULL
        xlags.current <- NULL
        veratts.current <- NULL
###############################
        for(j in (maxLag:1)){
            x.current <- InputNetwork[[i + j - 1]]
            x.lag <- numeric(length(Vunion))
            current.vnames <- network.vertex.names.1(x.current)
            x.lag[match(current.vnames, Vunion)] <- 1
            xlags.current <- cbind(xlags.current, x.lag)
            current.vstats <- vertexstats(x.current, gmode = "digraph")
            verstats.lag <- matrix(0, nrow = length(Vunion), ncol = nvertexstats)
            rownames(verstats.lag) <- Vunion
            verstats.lag[match(rownames(current.vstats),
                               rownames(verstats.lag)), ] <- current.vstats
            vstats.current <- cbind(vstats.current, verstats.lag)

#####################
            if(!is.na(VertexModelGroup)){
                veratts.lag <- numeric(length(Vunion))
                veratts.lag[match(current.vnames, Vunion)] <-
                    network::get.vertex.attribute(x.current, VertexModelGroup)
                veratts.current <- cbind(veratts.current, veratts.lag)
            }
        }
        
###############################
        net.current <- InputNetwork[[i + maxLag]]
        y.current <- numeric(length(Vunion))
        current.vnames <- network.vertex.names.1(net.current)
        y.current[match(current.vnames, Vunion)] <- 1

        y <- c(y, y.current)
        x.Nets <- rbind(x.Nets, xlags.current)
        x.vstats <- rbind(x.vstats, vstats.current)
        if(!is.na(VertexModelGroup)){
            x.vatts <- rbind(x.vatts, veratts.current)
        }
    }


    for(i in seq_len(NCOL(x.Nets))){
        colnames(x.Nets)[i] <- paste0("lag", i, sep = "")
    }
    if(!is.na(VertexModelGroup)){
        for(i in seq_len(NCOL(x.vatts))){
            colnames(x.vatts)[i] <- paste0("attrib", i, sep = "")
        }
    }
    cnames <- numeric(ncol(x.vstats))
    for(i in seq_len(maxLag)){
        for(j in seq_len(nvertexstats)){
            cnames[(i - 1)*nvertexstats + j] <- 
                paste0("Vstat",j,"Lag",i, sep = ".")
        }
    }
    colnames(x.vstats) <- cnames
    if(!is.na(VertexModelGroup)){
        XYdata <- cbind(y, x.Nets, x.vatts, x.vstats)
    } else{
        XYdata <- cbind(y, x.Nets, x.vstats)
    }
    
    ## colnames(XYdata) <- NULL
    XYdata <- as.data.frame(XYdata)
    colnames(XYdata)[1] <- "y"
    ## subset
    VertexLagvec <- c(t(VertexLagMatrix))
    if(maxLag > 1) {
        if(!is.na(VertexModelGroup)){
            VertexLagvec <-
                c(VertexLag, VertexAttLag,VertexLagvec)
        } else {
            VertexLagvec <-
                c(VertexLag,VertexLagvec)
        }
    }
    VertexLagvec <- c(1, VertexLagvec)
    VertexLagvec <- VertexLagvec == 1
    
    VertexRegout <- regEngine(XYdata[, VertexLagvec], method = regMethod)
    
    ## end of vertex only regression ##
    out <- list(VertexFit = VertexRegout,
                VertexStats = XYdata[, VertexLagvec],
                VertexStatsFull = XYdata)

}
