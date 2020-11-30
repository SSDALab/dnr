##' @title Parameter estimation for Vertex dynamics
##' @description Parameter estimation fro dynamic vertex case. The interface remaining almost identical to the static vertex one.
##' @param InputNetwork list of networks.
##' @param VertexStatsvec binary vector of size 8.
##' @param maxLag maximum lag, numeric.
##' @param VertexLag binary vector of length maxLag.
##' @param VertexLagMatrix binary matrix of size maxLag x 8.
##' @param VertexModelGroup Grouping term for vertex model. Must be from vertex attribute list.
##' @param VertexAttLag Lag vector for vertex group terms. Of length maxLag.
##' @param dayClass Any network level present time attribute vector. Here used to indicate week/weekend as 0/1.
##' @param EdgeModelTerms Model terms in edge model.
##' @param EdgeModelFormula Model formula in edge model.
##' @param EdgeGroup Group terms in edge model.
##' @param EdgeIntercept Intercept for edge model.
##' @param EdgeNetparam Network level parameter for edge model (currently only supported parameter is current network size).
##' @param EdgeExvar Extraneous variable for edge model.
##' @param EdgeLag binary vector of length maxLag.
##' @param EdgeLagMatrix binary matrix of dim maxLag x length(EdgeModelTerms)
##' @param regMethod Regression method. default: "bayesglm"
##' @param paramout T/F Should the parameter estimates be returned?
##' @details The Vertex model parameter list is as follows (Freeman degree, In degree, Out degree, Eigen Centrality, Between centrality, Info centrality, Closeness centrality, log k cycles, log size). For more details about the definitions of the terms, please refer to the vertexstats.R file, which implements all of these. The definitions are in sna or igraph.
##' @return list with following elements: \cr
##' EdgeCoef: edge coefficients. \cr
##' Edgemplematfull: MPLE matrix from edges. \cr
##' Edgemplemat: Subsetted MPLE matrix. \cr
##' VertexCoef: Coefficients from vertex. \cr
##' Vstats: Vertex statistics matrix.\cr
##' EdgePredictor0: Edge predictors with imputations with 0.\cr
##' EdgePredictor1: Edge predictors with imputations with 1. \cr
##' EdgePredictorNA: Edge predictors with imputations with NA. \cr
##' EdgeFit: Edge model. \cr 
##' VertexStatsFull: Vertex statistics matrix, full. \cr
##' VertexFit: Vertex model. \cr
##' @author Abhirup
##' @export
##' @examples
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
##' 
##' 
##' out <- paramVertex(InputNetwork = beach,
##'                    maxLag = 3,
##'                    VertexStatsvec = rep(1, nvertexstats),
##'                    VertexModelGroup = "regular",
##'                    VertexLag = rep(1, maxLag),
##'                    VertexLagMatrix = VertexLagMatrix,
##'                    dayClass = dayClass,
##'                    EdgeModelTerms = NA,
##'                    EdgeModelFormula = NA,
##'                    EdgeGroup = NA,
##'                    EdgeIntercept = c("edges"),
##'                    EdgeNetparam = c("logSize"),
##'                    EdgeExvar = NA,
##'                    EdgeLag = c(1, 1, 0),
##'                    paramout = TRUE)

paramVertex <- function(InputNetwork,
                        VertexStatsvec = rep(1, nvertexstats),
                        maxLag,
                        VertexLag = rep(1, maxLag),
                        VertexLagMatrix = matrix(1, maxLag,
                                                 length(VertexStatsvec)),
                        VertexModelGroup = NA,
                        VertexAttLag = rep(1, maxLag),
                        dayClass = NA,
                        EdgeModelTerms,
                        EdgeModelFormula,
                        EdgeGroup,
                        EdgeIntercept = c("edges"),
                        EdgeNetparam = NA,
                        EdgeExvar = NA,
                        EdgeLag = rep(1, maxLag),
                        EdgeLagMatrix = matrix(1, maxLag,
                                               length(EdgeModelTerms)),
                        regMethod = "bayesglm",
                        paramout = FALSE){
    ## Section 1: Vertex model
    
    gengroup <- function(input_network,group,net1){
        Vmax = input_network;
        
        Vmax.label = network::get.vertex.attribute(Vmax[[1]], attrname = group);
        
        Vmax = network.vertex.names(Vmax[[1]])
        
        foo.index = which(Vmax%in%network.vertex.names(net1)==T)
        
        if (is.null(group)){
            grouping = c(rep(1, floor(length(Vmax)/2)), rep(2, length(Vmax)-floor(length(Vmax)/2)))
        } else {
            grouping = Vmax.label
        }
        grouping.sub = grouping[foo.index]
        dirTF <- network::get.network.attribute(input_network[[i]], "directed")
        if(dirTF){
            grouping.perm = expand.grid(unique(grouping.sub),unique(grouping.sub))
            grouping.perm$indicator = seq_along(grouping.perm[,1])
            
            
            grouping.perm.full = expand.grid(unique(grouping),unique(grouping))
            grouping.perm.full$indicator = seq_along(grouping.perm.full[,1])
            
            foo = matrix(NA,ncol=length(grouping.sub), nrow=length(grouping.sub))
            
            
            grouping.terms = sapply(1:nrow(grouping.perm.full),function(i) paste(group,grouping.perm.full[i,1],grouping.perm.full[i,2],sep=''))
            for (i in 1:length(grouping.sub)){
                for(j in 1:length(grouping.sub)){
                    foo[i,j]=grouping.perm[,3] [which(grouping.perm[,1]==grouping.sub[i]&grouping.perm[,2]==grouping.sub[j])]
                    #print(i,j)
                }
            }
            for(i in grouping.perm.full$indicator){
                pos <- 1
                envir = as.environment(pos)
                assign(grouping.terms[i], (foo==grouping.perm.full$indicator[i])*1, envir = envir)
            }
        } else{
            
            if(length(unique(grouping.sub))>1){
                grouping.perm = data.frame(rbind(t(sapply(unique(grouping.sub), function(x)c(x,x))), t(combn(unique(grouping.sub),2))))
            } else{
                grouping.perm = data.frame(t(c(unique(grouping.sub),unique(grouping.sub))))
            }
            
            
            grouping.perm$indicator = seq_along(grouping.perm[,1])
            
            if(length(unique(grouping))>1){
                grouping.perm.full = data.frame(rbind(t(sapply(unique(grouping), function(x)c(x,x))), t(combn(unique(grouping),2))))
            } else{
                grouping.perm.full = data.frame(t(c(unique(grouping),unique(grouping))))
            }
            
            
            grouping.perm.full$indicator = seq_along(grouping.perm.full[,1])
            
            foo = matrix(NA,ncol=length(grouping.sub), nrow=length(grouping.sub))
            
            
            grouping.terms = sapply(1:nrow(grouping.perm.full),function(i) paste(group,grouping.perm.full[i,1],grouping.perm.full[i,2],sep=''))
            
            
            for (i in 1:length(grouping.sub)){
                for(j in 1:length(grouping.sub)){
                    ind1 = apply(grouping.perm[,1:2],1,function(x) all(x == c(grouping.sub[i], grouping.sub[j]))) + apply(grouping.perm[,1:2],1,function(x) all(x == c(grouping.sub[j], grouping.sub[i])));
                    ind1[which(ind1>0)] = 1;
                    foo[i,j]=grouping.perm[,3][which(ind1!=0)];
                    
                    #foo[i,j]=grouping.perm[,3] [which((grouping.perm[,1]==grouping.sub[i]&grouping.perm[,2]==grouping.sub[j])|(grouping.perm[,1]==grouping.sub[j]&grouping.perm[,2]==grouping.sub[i]))]
                    #print(i,j)
                }
            }
            
            for(i in grouping.perm.full$indicator){
                pos <- 1
                envir = as.environment(pos)
                assign(grouping.terms[i], (foo==grouping.perm.full$indicator[i])*1, envir = envir)
            }
        }
        
        grouping.edgecov = sapply(1:nrow(grouping.perm.full), function(x) paste0('edgecov(',grouping.terms[x],')',sep=''))
        return(grouping.edgecov)
    }
    
    
    
    ## Function for upscaling networks in a window
    windowUnion <- function(net.window){
        nwindow <- length(net.window)
        Vunion.window <- unique(unlist(lapply(net.window, network.vertex.names.1)))
        adjunion <- matrix(0, nrow = length(Vunion.window),
                           ncol = length(Vunion.window))
        rownames(adjunion) <- Vunion.window
        colnames(adjunion) <- Vunion.window
        common.window <- list()
        for(i in 1:nwindow){
            adjLag <- as.matrix(net.window[[i]])
            adjLagIdx <- match(rownames(adjLag), Vunion.window) 
            adjunion[adjLagIdx, adjLagIdx] <- adjLag
            dirTF <- network::get.network.attribute(net.window[[i]], "directed")
            netUnion <- network(adjunion, directed = dirTF)
            
            ## add vertex attributes
            vatbs <- network::list.vertex.attributes(net.window[[i]])
            vatbs <- vatbs[vatbs != "vertex.names"]
            for (vatb in vatbs){
                vatbsCommon <- numeric(length(Vunion.window))
                vatbsCommon[adjLagIdx] <- network::get.vertex.attribute(net.window[[i]], vatb)
                network::set.vertex.attribute(netUnion, vatb, vatbsCommon)
            }
            
            ## add edge attributes
            eatbs <- network::list.edge.attributes(net.window[[i]])
            for(eatb in eatbs){
                netEatb <- network::get.edge.attribute(net.window[[i]], eatb)
                network::set.edge.attribute(netUnion, eatb, netEatb)
            }
            
            common.window[[i]] <- netUnion
        }
        names(common.window) <- names(net.window)
        return(common.window)
    }
    
    
    
    genintercept <- function(intercept,grouping.edgecov,netname){
        if(is.na(grouping.edgecov)||is.na(intercept)){
            if(is.na(grouping.edgecov)){
                formula <- as.formula(paste(paste0(netname,"~",sep=""),paste(c(intercept,"edgecov(m)"),collapse = "+")))
            } else {
                formula <- as.formula(paste(paste0(netname,"~",sep=""),paste(c(grouping.edgecov,"edgecov(m)"),collapse = "+")))
            }
        }else{
            formula <- as.formula(paste(paste0(netname,"~",sep=""),paste(c(intercept,grouping.edgecov,"edgecov(m)"),collapse = "+")))
        }
        return(formula)
    }
    
    genformula <- function(model.formula,netname="net1"){
        term_formula = terms(model.formula);
        term_formula = attr(term_formula, 'term.labels');
        formula1 = as.formula(paste(paste0(netname,"~",sep=""),
                                    paste(c(term_formula,"edgecov(m)"),
                                          collapse= "+")))
        return(formula1)
    }
    
    
    
    
    ## Remove networks from InputNetowrk that are NA.
    InputNetwork <- rmNAnetworks(InputNetwork)
    dayClass <- na.omit(dayClass)    
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
    ## construct the response and predictors
    for(i in seq_len(netlength - maxLag)) {
        vstats.current <- NULL
        xlags.current <- NULL
        veratts.current <- NULL
        for(j in (maxLag:1)) { ## count down
            x.current <- InputNetwork[[i + j - 1]]
            x.lag <- numeric(length(Vunion))
            current.vnames <- network.vertex.names.1(x.current)
            x.lag[match(current.vnames, Vunion)] <- 1
            x.lag <- matrix(x.lag)
            colnames(x.lag) <- paste0("Lag.",(1+maxLag-j),sep="")
            xlags.current <- cbind(xlags.current, x.lag)
            current.vstats <- vertexstats(x.current, gmode = "digraph")
            verstats.lag <- matrix(0, nrow = length(Vunion), ncol = nvertexstats)
            rownames(verstats.lag) <- Vunion
            verstats.lag[match(rownames(current.vstats),
                               rownames(verstats.lag)), ] <- current.vstats
            colnames(verstats.lag) <- paste0(vstatNames, ".Lag.", (1+maxLag-j),sep="")
            vstats.current <- cbind(vstats.current, verstats.lag)
            
            ###############################
            #             if(!is.na(VertexModelGroup)){
            #                 veratts.lag <- numeric(length(Vunion))
            #                 veratts.lag[match(current.vnames, Vunion)] <-
            #                     network::get.vertex.attribute(x.current, VertexModelGroup)
            #                 veratts.current <- cbind(veratts.current, veratts.lag)
            #             }
            if(!any(is.na(VertexModelGroup))){
                veratts.lag.combined = NULL
                for (vertAtt in VertexModelGroup){
                    veratts.lag <- numeric(length(Vunion))
                    veratts.lag[match(current.vnames, Vunion)] <-
                        network::get.vertex.attribute(x.current, vertAtt)
                    veratts.lag <- matrix(veratts.lag)
                    colnames(veratts.lag) <- paste0(vertAtt, ".Lag", j, sep="")
                    veratts.lag.combined <- cbind(veratts.lag.combined, veratts.lag)
                }
                veratts.current <- cbind(veratts.current, veratts.lag.combined)
            }
            
        }
        net.current <- InputNetwork[[i + maxLag]]
        if(sum(!is.na(dayClass)) > 0) {
            ## we add network level exogenous variable here.
            day.current <- dayClass[i + maxLag]
            day.current <- matrix(day.current, nrow(xlags.current), 1)
            colnames(day.current) <- "NetworkAttribute"
            xlags.current <- cbind(xlags.current, day.current)
        }
        
        y.current <- numeric(length(Vunion))
        current.vnames <- network.vertex.names.1(net.current)
        y.current[match(current.vnames, Vunion)] <- 1
        
        y <- c(y, y.current)
        x.Nets <- rbind(x.Nets, xlags.current)
        x.vstats <- rbind(x.vstats, vstats.current)
        ###############################
        if(!any(is.na(VertexModelGroup))){
            x.vatts <- rbind(x.vatts, veratts.current)
        }
        
    }
    
    #     for(i in seq_len(NCOL(x.Nets) - 1)){
    #         colnames(x.Nets)[i] <- paste0("lag", i, sep = "")
    #     }
    #     colnames(x.Nets)[ncol(x.Nets)] <- "Day"   
    ###############################
    #     if(!any(is.na(VertexModelGroup))){
    #         for(i in seq_len(NCOL(x.vatts))){
    #             colnames(x.vatts)[i] <- paste0("attribLag", i, sep = "")
    #         }
    #     }
    
    
    #     cnames <- numeric(ncol(x.vstats))
    #     for(i in seq_len(maxLag)){
    #         for(j in seq_len(nvertexstats)){
    #             cnames[(i - 1)*nvertexstats + j] <- 
    #                 paste0(vstatNames[j],".Lag",i, sep = "")
    #         }
    #     }
    #     colnames(x.vstats) <- cnames
    
    y <- matrix(y)
    colnames(y) <- "y"
    
    if(!any(is.na(VertexModelGroup))){
        XYdata <- cbind.data.frame(y, x.Nets, x.vatts, x.vstats)
    } else{
        XYdata <- cbind.data.frame(y, x.Nets, x.vstats)
    }
    
    
    ## subset
    # terms: from right: vertexStats, vertexAttributes, networkAttributes, Lags, Y
    # vertexStats:
    VertexLagvec <- c(t(VertexLagMatrix))
    #vertexAttributes:
    if(!any(is.na(VertexModelGroup))){
        for(i in VertexModelGroup){
            VertexLagvec <- c(1, VertexLagvec)
        }
    }
    #networkAttributes:
    if(sum(!is.na(dayClass)) > 0) VertexLagvec <- c(1, VertexLagvec)
    #lags:
    VertexLagvec <- c(VertexLag, VertexLagvec)
    # Y
    VertexLagvec <- c(1, VertexLagvec)
    
    VertexLagvec <- VertexLagvec == 1
    
    VertexRegout <- regEngine(XYdata[, VertexLagvec], method = regMethod)
    
    
    ## end of vertex only regression ###
    
    ## section 2: edge conditional on vertex:
    matout <- NULL
    fullPredictorStack0 <- NULL
    fullPredictorStack1 <- NULL
    fullPredictorStackNA <- NULL
    if(is.na(EdgeGroup)){
        grouping.edgecov <- NA
    } else{
        grouping.edgecov <- gengroup(InputNetwork, EdgeGroup, InputNetwork[[1]])
    }
    
    for(i in 1:(netlength - maxLag)){
        net.window <- InputNetwork[i : (i + maxLag)]
        common.window <- windowUnion(net.window)
        net.current <- common.window[[maxLag + 1]]
        
        netsize.window <- network.size.1(net.current)
        m <- matrix(1:(netsize.window*netsize.window), nrow = netsize.window, ncol = netsize.window)
        
        
        #Construction of formula##
        
        if(is.na(EdgeIntercept)&&is.na(EdgeGroup)){
            csmodel <- NULL
        } else{
            formula <- genintercept(EdgeIntercept,grouping.edgecov,
                                    netname="net.current")
            #TODO: Add support for exogenous variable later
            if(!is.na(EdgeExvar)){
                formula <- as.formula(paste(c(formula,paste0("edgecov(",EdgeExvar,")",sep="")),collapse = "+"))
            }
            mplemat  <-  ergm::ergmMPLE(formula, output="matrix")
            csintercept <- as.matrix(mplemat$predictor[,1:(ncol(mplemat$predictor)-1)])
            colnames(csintercept) <- colnames(mplemat$predictor)[1:(ncol(mplemat$predictor)-1)]
            csmodel <- csintercept
        }
        #fit model terms
        #Comment: We always count down!
        if(sum(is.na(EdgeModelFormula)) == 0){
            for(j in maxLag:1){
                formula <- genformula(EdgeModelFormula,netname = "common.window[[j]]")
                mplemat.tmp  <-  ergm::ergmMPLE(formula, output="matrix");
                edgeLag.tmp <- mplemat.tmp$predictor[,1:(ncol(mplemat.tmp$predictor)-1)]
                csmodel <- cbind(csmodel, edgeLag.tmp);
            }
        }
        ###############################
        ## fit EdgeNetparam terms (only supported term is logSize)
        ## In future, more terms can be added
        ## REMEMBER to update the subsetting section!!
        if(!is.na(EdgeNetparam)) {
            if("logSize" %in% EdgeNetparam) {
                if(gmode == "digraph"){
                    currNetSize <- rep(netsize.window,
                                       netsize.window*(netsize.window - 1))
                    logCurrNetSize <- log(currNetSize + 1e-10)
                } else {
                    currNetSize <- rep(netsize.window,
                                       netsize.window*(netsize.window - 1)/2)
                    logCurrNetSize <- log(currNetSize + 1e-10)
                }
                csmodel <- cbind(csmodel, logCurrNetSize)
            }
        }
        
        ## We add network level exogenous variables here
        ## Update the subsetting section also.
        if(sum(!is.na(dayClass)) > 0) {
            day.current <- dayClass[i + maxLag]
            if(gmode == "digraph") {
                dayEffect <- rep(day.current,
                                 netsize.window*(netsize.window - 1))
            } else {
                dayEffect <- rep(day.current,
                                 netsize.window*(netsize.window - 1)/2)
            }
            csmodel <- cbind(csmodel, dayEffect)
        }
        
        ###############################
        
        ## fit the lag terms
        ## Comment: Counting down.
        if(maxLag > 1){
            if(gmode=='digraph'){
                lagstats <- matrix(0, ncol=maxLag, nrow=netsize.window*(netsize.window-1))
            } else{
                lagstats <- matrix(0, ncol=maxLag, nrow=netsize.window*(netsize.window-1)/2)
            }
            lagnames <- rep("lag",(maxLag))
            for (j in (maxLag):1){
                lagstats[,j] <- sna::gvectorize(common.window[[j]][,],mode=gmode, censor.as.na=F)
                lagnames[j] <- paste("lag",j,sep = "")
            }
            colnames(lagstats) <- lagnames
            #response
            y <- sna::gvectorize(net.current[,],mode=gmode, censor.as.na=FALSE)
            mat <- data.frame(cbind(y,csmodel,lagstats))
        } else {
            y <- sna::gvectorize(net.current[,],mode=gmode, censor.as.na=FALSE)
            mat <- data.frame(cbind(y,csmodel))
        }
        
        ## imputed matrix1: with 0
        ## create empty matrix of 0 of full dim. (length(Vunion)).
        ## num of rows: num of edges. num of cols: cols in mat.
        ncolXY <- ncol(mat)
        nEdges <- ifelse(gmode == "digraph", nv*(nv-1), nv*(nv-1)/2)
        fullPredictor0 <- matrix(0, nEdges, ncolXY)
        fullPredictor1 <- matrix(1, nEdges, ncolXY)
        fullPredictorNA <- matrix(NA, nEdges,ncolXY)
        ## get indicators of vertices present out of Vunion
        vnames.window <- network.vertex.names(net.current)
        vind <- match(vnames.window, Vunion)
        ##vind <- na.omit(vind)
        AdjMatUnion <- matrix(0, nrow = nv, ncol = nv)
        AdjMatUnion[vind, vind] <- 1
        misPattern <- na.omit(sna::gvectorize(AdjMatUnion, mode = gmode))
        subId <- which(misPattern == 1)
        fullPredictor0[subId, ] <- as.matrix(mat)
        fullPredictor1[subId, ] <- as.matrix(mat)
        fullPredictorNA[subId, ] <- as.matrix(mat)
        
        matout <- rbind(matout,mat)
        fullPredictorStack0 <- rbind(fullPredictorStack0, fullPredictor0)
        fullPredictorStack1 <- rbind(fullPredictorStack1, fullPredictor1)
        fullPredictorStackNA <- rbind(fullPredictorStackNA, fullPredictorNA)
    }
    
    ## subsetting
    if(sum(!is.na(csintercept)) > 0){
        lagvec <- rep(1,NCOL(csintercept))
    }
    if(sum(!is.na(EdgeExvar)) > 0) {
        lagvec <- c(lagvec, 1)
    }
    if(sum(!is.na(EdgeNetparam)) > 0) {
        lagvec <- c(lagvec, 1)
    }
    if(sum(!is.na(dayClass)) > 0) {
        lagvec <- c(lagvec, 1)
    }
    if(sum(!is.na(EdgeModelTerms)) > 0) lagvec <- c(lagvec,t(EdgeLagMatrix))
    if(maxLag > 1) lagvec <- c(lagvec,EdgeLag)
    lagvec <- c(1,lagvec)
    lagvec <- lagvec==1
    #regression
    if(paramout){
        out <- regEngine(matout[,lagvec],regMethod)
    } else out <- NULL
    
    return(list(EdgeCoef=out$coef,
                EdgeFit = out,
                Edgemplematfull=matout,
                Edgemplemat=matout[,lagvec],
                EdgePredictor0 = fullPredictorStack0[, lagvec],
                EdgePredictor1 = fullPredictorStack1[, lagvec],
                EdgePredictorNA = fullPredictorStackNA[, lagvec],
                VertexCoef = VertexRegout$coef,
                Vstats = XYdata[, VertexLagvec],
                VertexStatsFull = XYdata,
                VertexFit = VertexRegout))
}

