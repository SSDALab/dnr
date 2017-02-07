##' Parameter estimation for Vertex model only for a list of dynamic networks.
##' @param InputNetwork Input network list.
##' @param VertexStatsvec Binary vector of size 8, indicating vertex model.
##' @param maxLag maximum lag.
##' @param VertexLag Binary vector of size maxLag, indicating Lag terms in the model.
##' @param VertexLagMatrix Binary matrix indicating lagged vertex statistics in
##'     the model.
##' @param regMethod one of "glm", "glmnet", "bayesglm"
##' @return List of 3 elements:\cr
##' VertexFit: Output from regEngine. \cr
##' VertexStats: Subsetted vertex stats matrix. \cr
##' VertexStatsFull: Full matrix of vertex stats.
##' @author Abhirup
##' @export
paramVertexOnly <- function(InputNetwork,
                            VertexStatsvec = rep(1, 8),
                            maxLag,
                            VertexLag = rep(1, maxLag),
                            VertexLagMatrix = matrix(1, maxLag,
                                                     length(VertexStatsvec)),
                            regMethod = "bayesglm"){
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
        dirTF <- network::get.network.attribute(input_network[[i]], "directed")
        if(dirTF){
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

    ## Remove networks from InputNetowrk that are NA.
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
    
    VertexRegout <- regEngine(XYdata[, VertexLagvec], method = regMethod)
    
    ## end of vertex only regression ##
    out <- list(VertexFit = VertexRegout,
                VertexStats = XYdata[, VertexLagvec],
                VertexStatsFull = XYdata)
}
