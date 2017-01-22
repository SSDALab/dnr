##' .. content for \description{} (no empty lines) ..
##' 
##' .. content for \details{} ..
##' @title Parameter estimation for Vertex dynamics
##' @param InputNetwork list of networks.
##' @param VertexStatsvec binary vector of size 8.
##' @param maxLag maximum lag, numeric.
##' @param VertexLag binary vector of length maxLag.
##' @param VertexLagMatrix binary matrix of size maxLag x 8.
##' @param EdgeModelTerms Model terms in edge model.
##' @param EdgeModelFormula Model formula in edge model.
##' @param EdgeGroup Group terms in edge model.
##' @param EdgeIntercept Intercept for edge model.
##' @param EdgeExvar Extraneous variable for edge model.
##' @param EdgeLag binary vector of length maxLag.
##' @param EdgeLagMatrix binary matrix of dim maxLag x length(EdgeModelTerms)
##' @param regMethod Regression method. default: "bayesglm"
##' @param paramout T/F Should the parameter estimates be returned?
##' @return list with following elements:
##' EdgeCoef: edge coefficients.
##' Edgemplematfull: MPLE matrix from edges.
##' Edgemplemat: Subsetted MPLE matrix.
##' VertexCoef: Coefficients from vertex.
##' Vstats: Vertex statistics matrix.
##' @author Abhirup
##' @examples
##' InputNetwork <- beach
##' maxLag <- 3
##' regMethod <- "bayesglm"
##' 
##' lagVector <- rep(1, maxLag)
##' EdgeModelTerms <- c("triadcensus.003", "triadcensus.012", "triadcensus.102",
##'                     "triadcensus.021D")
##' EdgeModelFormula <- net ~ triadcensus(0:3)
##' EdgeGroup <- "group1"
##' EdgeIntercept <- c("edges")
##' EdgeExvar <- NA
##' EdgeLag <- rep(1, maxLag)
##' 
##' EdgeLagMatrix <- rbind(c(1, 1, 0, 1),
##'                        c(0, 0, 1, 1),
##'                        c(1, 0, 0, 0))
##' VertexLag <- c(1, 1, 0)
##' VertexLagMatrix <-  matrix(1, maxLag,
##'                            length(VertexStatsvec))
##' 
##' out <- paramVertex(InputNetwork = beach,
##'                    maxLag = 3,
##'                    VertexLag = VertexLag,
##'                    EdgeModelTerms = c("triadcensus.003", "triadcensus.012", "triadcensus.102","triadcensus.021D"),
##'                    EdgeModelFormula = net ~ triadcensus(0:3),
##'                    EdgeGroup = "group1",
##'                    EdgeIntercept = c("edges"),
##'                    paramout = TRUE)

paramVertex <- function(InputNetwork,
                        VertexStatsvec = rep(1, 8),
                        maxLag,
                        VertexLag = rep(1, maxLag),
                        VertexLagMatrix = matrix(1, maxLag,
                                                 length(VertexStatsvec)),
                        EdgeModelTerms,
                        EdgeModelFormula,
                        EdgeGroup,
                        EdgeIntercept = c("edges"),
                        EdgeExvar = NA,
                        EdgeLag = rep(1, maxLag),
                        EdgeLagMatrix = matrix(1, maxLag,
                                               length(EdgeModelTerms)),
                        regMethod = "bayesglm",
                        paramout = FALSE){
    ## Section 1: Vertex model
    ## supporting functions:
    ## network size (that allows NA)
    network.size.1 <- function(x){
        if(!is.network(x)){
            if(is.na(x)) return(0)
        } else return(network.size(x))
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
        if(!is.network(x)){
            if(is.na(x)) return (NA)
        } else {
            return (network.vertex.names(x))
        }
    }
    
    ## regression
    ## TODO: add options for models with intercept (regression without -1)
    
    edgecoeff <- function(output.edge,method='glmnet'){
        if(method == 'glmnet'){
            if(is.na(lambda)){
                blogfit.select.edge = cv.glmnet(data.matrix(output.edge[,-1]),as.vector(output.edge[,1]), family="binomial", alpha=alpha.glmnet)
                lambda = blogfit.select.edge$lambda.min
            }

            blogfit.sim.edge = glmnet(data.matrix(output.edge[,-1]),as.vector(output.edge[,1]),family="binomial", alpha=alpha.glmnet, lambda=lambda, intercept=F);

            coef.edge=setNames(as.vector(blogfit.sim.edge$beta), colnames(output.edge)[-1]);

            return(list(coef.edge=coef.edge, lambda=lambda, std.error=NA))
        } else if (method == 'glm'){
            blogfit.sim.edge = glm(y~.-1, data= output.edge,family=binomial(logit));

            coef.edge=setNames(as.vector(blogfit.sim.edge$coefficients), colnames(output.edge)[-1]);
            ## set NA to 0
            coef.edge[which(is.na(coef.edge))] = 0;
            return(list(coef.edge=coef.edge, lambda=NA, std.error=summary(blogfit.sim.edge)$coefficients[,2]))

        } else if (method == 'bayesglm'){
            blogfit.sim.edge = bayesglm(y~.-1, data= output.edge,family=binomial(logit));

            coef.edge=setNames(as.vector(blogfit.sim.edge$coefficients), colnames(output.edge)[-1]);
            ## set NA to 0
            coef.edge[which(is.na(coef.edge))] = 0;
            return(list(coef.edge=coef.edge,lambda=NA, std.error=summary(blogfit.sim.edge)$coefficients[,2]))
        }
    }

    ## supporting functions for edge case:

    ## Grouping term generator function
    #'Grouping term generator function
    #'@title gengroup
    #'@description: Generates the grouping terms from the formula group variable input by the user.
    #'@param input_network
    #'@param group: vector of strings, indicating grouping variables.
    #'@param net1: current network
    #'@return a vector of characters, that should be appended (+) at the begining of a formula to feed into ergm
    #'

    gengroup <- function(input_network,group,net1){
        Vmax = input_network;

        Vmax.label = get.vertex.attribute(Vmax[[1]], attrname = group);

        Vmax = network.vertex.names(Vmax[[1]])

        foo.index = which(Vmax%in%network.vertex.names(net1)==T)

        if (is.null(group)){
            grouping = c(rep(1, floor(length(Vmax)/2)), rep(2, length(Vmax)-floor(length(Vmax)/2)))
        } else {
            grouping = Vmax.label
        }
        grouping.sub = grouping[foo.index]

        if(directed){
            grouping.perm = expand.grid(unique(grouping.sub),unique(grouping.sub))
            grouping.perm$indicator = seq_along(grouping.perm[,1])


            grouping.perm.full = expand.grid(unique(grouping),unique(grouping))
            grouping.perm.full$indicator = seq_along(grouping.perm.full[,1])

            foo = matrix(,ncol=length(grouping.sub), nrow=length(grouping.sub))


            grouping.terms = sapply(1:nrow(grouping.perm.full),function(i) paste(group,grouping.perm.full[i,1],grouping.perm.full[i,2],sep=''))
            for (i in 1:length(grouping.sub)){
                for(j in 1:length(grouping.sub)){
                    foo[i,j]=grouping.perm[,3] [which(grouping.perm[,1]==grouping.sub[i]&grouping.perm[,2]==grouping.sub[j])]
                                        #print(i,j)
                }
            }
            for(i in grouping.perm.full$indicator){
                assign(grouping.terms[i], (foo==grouping.perm.full$indicator[i])*1, env = globalenv())
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

            foo = matrix(,ncol=length(grouping.sub), nrow=length(grouping.sub))


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
                assign(grouping.terms[i], (foo==grouping.perm.full$indicator[i])*1, env = globalenv())
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


    #'Generate Intercept
    #'@title genintercept
    #'@param intercept: intercept vector
    #'@param grouping.edgecov: grouping terms from gengroup
    #'@param netname: name of the network
    #'@return formula for just the intercept, to be fed to ergm
    #'
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

    #'Generate formula for feeding into ergm
    #'@title genformula
    #'@param model.formula: user given formula
    #'@param netname: left side network name to be used. For now, I use "net1"
    #'@return formula to be fed to ergm
    #'

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
    netlength <- length(InputNetwork)
    Vunion <- unique(unlist(lapply(InputNetwork, network.vertex.names.1)))
    Vunion <- na.omit(Vunion)
    nv <- length(Vunion)

    ## is the network directed?
    ## we are using 1st network as representative network hoping these
    ## attributes wont change.
    isdirected <- get.network.attribute(InputNetwork[[1]], "directed")
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
    ## construct the response and predictors
    for(i in seq_len(netlength - maxLag)) {
        vstats.current <- NULL
        xlags.current <- NULL
        for(j in (maxLag:1)) { ## count down
            x.current <- InputNetwork[[i + j - 1]]
            x.lag <- numeric(length(Vunion))
            current.vnames <- network.vertex.names.1(x.current)
            x.lag[match(current.vnames, Vunion)] <- 1
            xlags.current <- cbind(xlags.current, x.lag)
            current.vstats <- vertexstats(x.current, gmode = "digraph")
            verstats.lag <- matrix(0, nrow = length(Vunion), ncol = nvertstats)
            rownames(verstats.lag) <- Vunion
            verstats.lag[match(rownames(current.vstats),
                               rownames(verstats.lag)), ] <- current.vstats
            vstats.current <- cbind(vstats.current, verstats.lag)
        }
        net.current <- InputNetwork[[i + maxLag]]
        y.current <- numeric(length(Vunion))
        current.vnames <- network.vertex.names.1(net.current)
        y.current[match(current.vnames, Vunion)] <- 1

        y <- c(y, y.current)
        x.Nets <- rbind(x.Nets, xlags.current)
        x.vstats <- rbind(x.vstats, vstats.current)
    }

    for(i in seq_len(NCOL(x.Nets))){
        colnames(x.Nets)[i] <- paste0("lag", i, sep = "")
    }

    cnames <- numeric(ncol(x.vstats))
    for(i in seq_len(maxLag)){
        for(j in seq_len(nvertexstats)){
            cnames[(i - 1)*nvertexstats + j] <- 
                paste0("Vstat",j,"Lag",i, sep = ".")
        }
    }
    colnames(x.vstats) <- cnames
    XYdata <- cbind(y, x.Nets, x.vstats)

    
    ## colnames(XYdata) <- NULL
    XYdata <- as.data.frame(XYdata)
    colnames(XYdata)[1] <- "y"
    ## subset
    VertexLagvec <- c(t(VertexLagMatrix))
    if(maxLag > 1) VertexLagvec <- c(VertexLag, VertexLagvec)
    VertexLagvec <- c(1, VertexLagvec)
    VertexLagvec <- VertexLagvec == 1
    
    VertexRegout <- edgecoeff(XYdata[, VertexLagvec], method = regMethod)

    ## end of vertex only regression ###

    ## section 2: edge conditional on vertex:
    matout <- NULL
    fullPredictorStack0 <- NULL
    fullPredictorStack1 <- NULL
    fullPredictorStackNA <- NULL
    if(is.na(EdgeGroup)){
        grouping.edgecov <- NA
    } else{
        grouping.edgecov <- gengroup(InputNetowrk, EdgeGroup, InputNetwork[[1]])
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
                formula <- as.formula(paste(c(formula,paste0("edgecov(",exvar,")",sep="")),collapse = "+"))
            }
            mplemat  <-  ergmMPLE(formula, output="matrix")
            csintercept <- as.matrix(mplemat$predictor[,1:(ncol(mplemat$predictor)-1)])
            colnames(csintercept) <- colnames(mplemat$predictor)[1:(ncol(mplemat$predictor)-1)]
            csmodel <- csintercept
        }
                                        #fit model terms
                                        #Comment: We always count down!
        if(sum(is.na(model.formula)) == 0){
            for(j in maxLag:1){
                formula <- genformula(model.formula,netname = "common.window[[j]]")
                mplemat.tmp  <-  ergmMPLE(formula, output="matrix");
                edgeLag.tmp <- mplemat.tmp$predictor[,1:(ncol(mplemat.tmp$predictor)-1)]
                csmodel <- cbind(csmodel, edgeLag.tmp);
            }
        }
                                        #fit the lag terms
                                        #Comment: Counting down.
        if(maxLag > 1){
            if(gmode=='digraph'){
                lagstats <- matrix(0, ncol=maxLag, nrow=netsize.window*(netsize.window-1))
            } else{
                lagstats <- matrix(0, ncol=maxLag, nrow=netsize.window*(netsize.window-1)/2)
            }
            lagnames <- rep("lag",(maxLag))
            for (j in (maxLag):1){
                lagstats[,j] <- gvectorize(common.window[[j]][,],mode=gmode, censor.as.na=F)
                lagnames[j] <- paste("lag",j,sep = "")
            }
            colnames(lagstats) <- lagnames
                                        #response
            y <- gvectorize(net.current[,],mode=gmode, censor.as.na=FALSE)
            mat <- data.frame(cbind(y,csmodel,lagstats))
        } else {
            y <- gvectorize(net.current[,],mode=gmoded, censor.as.na=FALSE)
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
        misPattern <- na.omit(gvectorize(AdjMatUnion, mode = gmode))
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
    lagvec <- rep(1,NCOL(csintercept))
    if(sum(is.na(EdgeModelTerms)) != 0) lagvec <- c(lagvec,t(EdgeLagMatrix))
    if(maxLag > 1) lagvec <- c(lagvec,EdgeLag)
    lagvec <- c(1,lagvec)
    lagvec <- lagvec==1
                                        #regression
    if(paramout){
        out <- edgecoeff(matout[,lagvec],regMethod)
    } else out <- NULL

    return(list(EdgeCoef=out,
                Edgemplematfull=matout,
                Edgemplemat=matout[,lagvec],
                EdgePredictor0 = fullPredictorStack0,
                EdgePredictor1 = fullPredictorStack1,
                EdgePredictorNA = fullPredictorStackNA,
                VertexCoef = VertexRegout,
                Vstats = XYdata[, VertexLagvec]))
}
