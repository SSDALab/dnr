##' @title Simulation Engine for dynamic Vertex case.
##' @description Simulation engine for dynamic networks with variable number of vertices. 
##' Implements exponential family based hierachical model for vertice and the edges.
##' @param InputNetwork List of input networks
##' @param numSim number of time points to simulate
##' @param maxLag maximum Lag
##' @param VertexStatsvec Binary vector for vertex model.
##' @param VertexLag vector of lag for vertex
##' @param VertexLagMatrix matrix of lags for vertex stats.
##' @param VertexModelGroup Group term for vertex model.
##' @param VertexAttLag Lag vector for group term for vertex.
##' @param dayClassObserved Observed day class.
##' @param dayClassFuture Dayclass vector for future, must be of size numsim.
##' @param EdgeModelTerms Edge Model terms
##' @param EdgeModelFormula Edge model formula
##' @param EdgeGroup edge group term
##' @param EdgeIntercept edge intercept
##' @param EdgeNetparam edge network parameter name
##' @param EdgeExvar edge extraneous variable
##' @param EdgeLag edge Lag vector
##' @param EdgeLagMatrix edge lag matrix
##' @param regMethod regression method. "bayesglm" by default
##' @param paramout T/F on if regression needs to run.
##' @return List with following elements:
##' SimNetwork: Output Networks
##' EdgeParameterMat: Matrix of edge parameter
##' VertexParameterMat: Matrix of Vertex parameters.
##' @export
##'@examples
##'\dontrun{
##' nvertexstats <- 9
##' maxLag = 3
##' VertexLag = rep(1, maxLag)
##' VertexLagMatrix <- matrix(0, maxLag, nvertexstats)
##' VertexLagMatrix[, c(4, 7)] <- 1
##' VertexLagMatrix[c(2,3),7] <- 0
##' 
##' getWeekend <- function(z){
##'     weekends <- c("Saturday", "Sunday")
##'     if(!network::is.network(z)){
##'         if(is.na(z)) return(NA)
##'     } else {
##'          zDay <- get.network.attribute(z, attrname = "day")
##'          out <- ifelse(zDay %in% weekends, 1, 0)
##'          return(out)   
##'     }
##' }
##' 
##' dayClass <- numeric(length(beach))
##' for(i in seq_along(dayClass)) {
##'     dayClass[i] <- getWeekend(beach[[i]])
##' }
##' dayClass <- na.omit(dayClass)
##' simResult <- suppressWarnings(engineVertex(InputNetwork = beach,
##'                           numSim = 5,
##'                           maxLag = 3,
##'                           VertexStatsvec = rep(1, nvertexstats),
##'                           VertexModelGroup = "regular",
##'                           VertexAttLag = rep(1, maxLag),
##'                           VertexLag = rep(1, maxLag),
##'                           VertexLagMatrix = VertexLagMatrix,
##'                           dayClassObserved = dayClass,
##'                           dayClassFuture = c(1, 0, 0, 0, 0),
##'                           EdgeModelTerms = NA,
##'                           EdgeModelFormula = NA,
##'                           EdgeGroup = NA,
##'                           EdgeIntercept = c("edges"),
##'                           EdgeNetparam = c("logSize"),
##'                           EdgeExvar = NA,
##'                           EdgeLag = c(0, 1, 0),
##'                           paramout = TRUE
##'                           ))}
##' 
##' @author Abhirup




engineVertex <- function(InputNetwork,
                         numSim,
                         maxLag,
                         VertexStatsvec = rep(1, nvertexstats),
                         VertexLag = rep(1, maxLag),
                         VertexLagMatrix = matrix(1, maxLag,
                                                  length(VertexStatsvec)),
                         VertexModelGroup = NA,
                         VertexAttLag = rep(1, maxLag),
                         dayClassObserved = NA,
                         dayClassFuture = NA,
                         EdgeModelTerms,
                         EdgeModelFormula,
                         EdgeGroup = NA,
                         EdgeIntercept = c("edges"),
                         EdgeNetparam = NA,
                         EdgeExvar = NA,
                         EdgeLag = rep(1, maxLag),
                         EdgeLagMatrix = matrix(1, maxLag,
                                                length(EdgeModelTerms)),
                         regMethod = "bayesglm",
                         paramout = TRUE){
    InputNetwork <- rmNAnetworks(InputNetwork)
    dayClassObserved <- na.omit(dayClassObserved)
    Vunion <- unique(unlist(lapply(InputNetwork, network.vertex.names.1)))
    netlength <- length(InputNetwork)
    repfac <- netlength - maxLag
    nv <- length(Vunion)

    simulatedNetworks <- list()
    EdgeParameters <- list()
    VertexParameters <- list()
    
    isdirected <- network::get.network.attribute(InputNetwork[[1]], "directed")
    if(isdirected) {
        gmode <- "digraph"
    } else gmode <- "graph"
    
    
    for(simcount in seq_len(numSim)){
        print(simcount)
        Out.param <- paramVertex(InputNetwork = InputNetwork,
                                 VertexStatsvec = VertexStatsvec,
                                 maxLag = maxLag,
                                 VertexLag = VertexLag,
                                 VertexLagMatrix = VertexLagMatrix,
                                 VertexModelGroup = VertexModelGroup,
                                 VertexAttLag = VertexAttLag,
                                 dayClass = dayClassObserved,
                                 EdgeModelTerms = EdgeModelTerms,
                                 EdgeModelFormula = EdgeModelFormula,
                                 EdgeGroup = EdgeGroup,
                                 EdgeIntercept = EdgeIntercept,
                                 EdgeNetparam = EdgeNetparam,
                                 EdgeExvar = EdgeExvar,
                                 EdgeLag = EdgeLag,
                                 EdgeLagMatrix = EdgeLagMatrix,
                                 regMethod = regMethod,
                                 paramout = paramout)
        XYdata <- Out.param$Vstats[, -1]
        VertexCoeffs <- Out.param$VertexCoef
        ## check if nrows of XYdata is multiple of length(Vunion)
        if(nrow(XYdata) %% length(Vunion) != 0) stop("Wrong dimension of XYdata")
        Vstats.smooth <- matrix(0, nv, ncol(XYdata))
        for(i in seq_len(repfac)){
            Vstats.smooth <- Vstats.smooth + XYdata[(((i - 1)*nv + 1):(i*nv)),]
        }
        Vstats.smooth <- as.matrix(Vstats.smooth)/repfac
        Vertex.predictors <- Vstats.smooth %*% VertexCoeffs

        ## fix dayClass, if Day is present
        if(sum(!is.na(dayClassObserved)) > 0) {
            Vstats.smooth[, "Day"] <- dayClassFuture[simcount]
        }

        ## generate vertices:

        Vertex.new <- rbinom(nv, 1, ilogit(Vertex.predictors))
        names(Vertex.new) <- Vunion
        ## construct the empty networks with vertices give by
        
        ## Predict the edges
        InputMPLEmat <- as.matrix(Out.param$EdgePredictor0[, -1])
        nEdges <- ifelse(gmode == "digraph", nv*(nv -1), nv*(nv - 1)/2)
        repfac <- netlength - maxLag
        smmpleMat <- matrix(0, nEdges, ncol(InputMPLEmat))
        for(i in seq_len(repfac)){
            smmpleMat <- smmpleMat + InputMPLEmat[(((i-1)*nEdges+1):(i*nEdges)), ]
        }
        EdgeCoef <- Out.param$EdgeCoef
        ## name matching for updating the predictors.
        colnames(smmpleMat) <- names(EdgeCoef)
        ## list of predictors that needs updating: dayEffect
        if(sum(!is.na(dayClassFuture)) > 0){
            smmpleMat[, "dayEffect"] <- dayClassFuture[simcount]
        }
        smmpleMat <- as.matrix(smmpleMat/repfac)
        inputPred <- smmpleMat %*% EdgeCoef
        ## create an empty full sized network
        adjUnion <- matrix(0, nv, nv)
        rownames(adjUnion) <- Vunion
        colnames(adjUnion) <- Vunion
        dirTF <- network::get.network.attribute(InputNetwork[[1]], "directed")
        netUnion <- network(adjUnion, directed = dirTF)
        ## add vertex attributes:
        vatbs <- network::list.vertex.attributes(InputNetwork[[1]])
        vatbs <- vatbs[vatbs != "vertex.names"]
        for(vatb in vatbs){
            vatbUnion <- numeric(nv)
            names(vatbUnion) <- Vunion
            ##    vatbUnion <- ## iterate and match for vatb
            for(net in InputNetwork){
                netVnames <- network.vertex.names.1(net)
                vIdx <- match(netVnames, Vunion)
                vatbUnion[vIdx] <- network::get.vertex.attribute(net, vatb)
            }
            network::set.vertex.attribute(netUnion, vatb, vatbUnion)
        }

        edgeProbs <- ungvectorize(inputPred, nv, gmode)
        netUnion %n% "X" <- edgeProbs
        newNet.tmp <- simulate(ergm(netUnion ~ edgecov("X")))

        ## remove vertices that are absent in vertex.new
        ## get vertex names in newNet
        ## match with positions in Vunion
        vertexTodelete <- which(Vertex.new == 0)
        newNet <- network::delete.vertices(newNet.tmp, vertexTodelete)

        simulatedNetworks[[simcount]] <- newNet
        ## Let InputNetwork grow, no one dies.
        InputNetwork[[netlength + simcount]] <- newNet 
        dayClassObserved[netlength + simcount] <- dayClassFuture[simcount]
        ## parameters time series
        EdgeParameters[[simcount]] <- Out.param$EdgeCoef
        VertexParameters[[simcount]] <- Out.param$VertexCoef
    }
    EdgeParameterMat <- do.call(rbind, EdgeParameters)
    VertexParameterMat <- do.call(rbind, VertexParameters)

    return(list(SimNetwork = simulatedNetworks,
                EdgeParameterMat = EdgeParameterMat,
                VertexParameterMat = VertexParameterMat))
}



