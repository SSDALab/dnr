##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Simulation Engine for dynamic Vertex case.
##' @param InputNetwork List of input networks
##' @param numSim number of time points to simulate
##' @param maxLag maximum Lag
##' @param VertexLag vector of lag for vertex
##' @param VertexLagMatrix matrix of lags for vertex stats.
##' @param EdgeModelTerms Edge Model terms
##' @param EdgeModelFormula Edge model formula
##' @param EdgeGroup edge group term
##' @param EdgeIntercept edge intercept
##' @param EdgeExvar edge extraneous variable
##' @param EdgeLag edge Lag vector
##' @param EdgeLagMatrix edge lag matrix
##' @param regMethod regression method. "bayesglm" by default
##' @param paramout T/F on if regression needs to run.
##' @return List with following elements:
##' SimNetwork: Output Networks
##' EdgeParameterMat: Matrix of edge parameter
##' VertexParameterMat: Matrix of Vertex parameters.
##'@examples
##' InputNetwork <- beach
##' numSim <- 2
##' maxLag <- 3
##' EdgeModelTerms <- c("triadcensus.003", "triadcensus.012", "triadcensus.102",
##'                     "triadcensus.021D")
##' EdgeModelFormula <- net ~ triadcensus(0:3)
##' VertexLagMatrix <-  matrix(1, maxLag,
##'                            length(VertexStatsvec))
##' 
##' ## test function:
##' simResult <- engineVertex(InputNetwork = beach,
##'                           numSim = 5,
##'                           maxLag = 3,
##'                           EdgeModelTerms = c("triadcensus.003",
##'                                              "triadcensus.012",
##'                                              "triadcensus.102",
##'                                              "triadcensus.021D"),
##'                           EdgeModelFormula = net ~ triadcensus(0:3),
##'                           )
##' @author Abhirup
engineVertex <- function(InputNetwork,
                          numSim,
                          maxLag,
                          VertexLag = rep(1, maxLag),
                          VertexLagMatrix = matrix(1, maxLag,
                                                   length(VertexStatsvec)),
                          EdgeModelTerms,
                          EdgeModelFormula,
                          EdgeGroup = NA,
                          EdgeIntercept = c("edges"),
                          EdgeExvar = NA,
                          EdgeLag = rep(1, maxLag),
                          EdgeLagMatrix = matrix(1, maxLag,
                                                 length(EdgeModelTerms)),
                          regMethod = "bayesglm",
                          paramout = TRUE){
    InputNetwork <- rmNAnetworks(InputNetwork)
    Vunion <- unique(unlist(lapply(InputNetwork, network.vertex.names.1)))
    netlength <- length(InputNetwork)
    repfac <- netlength - maxLag
    nv <- length(Vunion)

    simulatedNetworks <- list()
    EdgeParameters <- list()
    VertexParameters <- list()
    
    for(simcount in seq_len(numSim)){
        print(simcount)
        Out.param <- paramVertex(InputNetwork = InputNetwork,
                                 VertexStatsvec = VertexStatsvec,
                                 maxLag = maxLag,
                                 VertexLag = VertexLag,
                                 VertexLagMatrix = VertexLagMatrix,
                                 EdgeModelTerms = EdgeModelTerms,
                                 EdgeModelFormula = EdgeModelFormula,
                                 EdgeGroup = EdgeGroup,
                                 EdgeIntercept = EdgeIntercept,
                                 EdgeExvar = EdgeExvar,
                                 EdgeLag = EdgeLag,
                                 EdgeLagMatrix = EdgeLagMatrix,
                                 regMethod = regMethod,
                                 paramout = paramout)
        XYdata <- Out.param$Vstats[, -1]
        InputMPLEmat <- Out.param$Edgemplemat
        VertexCoeffs <- Out.param$VertexCoef$coef.edge
        ## check if nrows of XYdata is multiple of length(Vunion)
        if(nrow(XYdata) %% length(Vunion) != 0) stop("Wrong dimension of XYdata")
        Vstats.smooth <- matrix(0, nv, ncol(XYdata))
        for(i in seq_len(repfac)){
            Vstats.smooth <- Vstats.smooth + XYdata[(((i - 1)*nv + 1):(i*nv)),]
        }
        Vstats.smooth <- as.matrix(Vstats.smooth)/repfac
        Vertex.predictors <- Vstats.smooth %*% VertexCoeffs

        ## generate vertices:

        Vertex.new <- rbinom(nv, 1, ilogit(Vertex.predictors))
        names(Vertex.new) <- Vunion
        ## construct the empty networks with vertices give by
        
        ## Predict the edges
        InputMPLEmat <- as.matrix(out$EdgePredictor0[, -1])
        nEdges <- ifelse(gmode == "digraph", nv*(nv -1), nv*(nv - 1)/2)
        repfac <- netlength - maxLag
        smmpleMat <- matrix(0, nEdges, ncol(InputMPLEmat))
        for(i in seq_len(repfac)){
            smmpleMat <- smmpleMat + InputMPLEmat[(((i-1)*nEdges+1):(i*nEdges)), ]
        }
        EdgeCoef <- out$EdgeCoef$coef.edge
        smmpleMat <- smmpleMat/repfac
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
        newNet <- delete.vertices(newNet.tmp, vertexTodelete)

        simulatedNetworks[[simcount]] <- newNet
        ## Let InputNetwork grow, no one dies.
        InputNetwork[[netlength + simcount]] <- newNet 

        ## parameters time series
        EdgeParameters[[simcount]] <- Out.param$EdgeCoef$coef.edge
        VertexParameters[[simcount]] <- Out.param$VertexCoef$coef.edge
    }
    EdgeParameterMat <- do.call(rbind, EdgeParameters)
    VertexParameterMat <- do.call(rbind, VertexParameters)

    return(list(SimNetwork = simulatedNetworks,
                EdgeParameterMat = EdgeParameterMat,
                VertexParameterMat = VertexParameterMat))
}

