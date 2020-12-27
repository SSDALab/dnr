##' Parameter estimation for Vertex model only for a list of dynamic networks.
##' @param InputNetwork Input network list.
##' @param VertexStatsvec Binary vector of size 9, indicating vertex model.
##' @param maxLag maximum lag.
##' @param VertexLag Binary vector of size maxLag, indicating Lag terms in the model.
##' @param VertexLagMatrix Binary matrix indicating lagged vertex statistics in
##'     the model.
##' @param dayClass Any network level present time attribute vector. Here used to indicate week/weekend as 0/1.
##' @param regMethod one of "glm", "glmnet", "bayesglm"
##' @return List of 3 elements:\cr
##' VertexFit: Output from regEngine. \cr
##' VertexStats: Subsetted vertex stats matrix. \cr
##' VertexStatsFull: Full matrix of vertex stats.
##' @examples
##' nvertexstats <- 9
##' maxLag = 3
##' VertexLag = rep(1, maxLag)
##' VertexLagMatrix <- matrix(0, maxLag, nvertexstats)
##' VertexLagMatrix[, c(4, 7)] <- 1
##' VertexLagMatrix[c(2,3),7] <- 0
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
##' ## for(i in 1:31) print(getWeekend(beach[[i]]))
##' ## generate a vector of network level exogenous variable
##' dayClass <- numeric(length(beach))
##' for(i in seq_along(dayClass)) {
##'     dayClass[i] <- getWeekend(beach[[i]])
##' }
##' out <- paramVertexOnly(InputNetwork = beach,
##'                        maxLag = 3,
##'                        VertexStatsvec = rep(1, nvertexstats),
##'                        VertexLag = rep(1, maxLag),
##'                        VertexLagMatrix = VertexLagMatrix,
##'                        dayClass = dayClass)
##' @author Abhirup
##' @export

paramVertexOnly <- function(InputNetwork,
                            VertexStatsvec = rep(1, nvertexstats),
                            maxLag,
                            VertexLag = rep(1, maxLag),
                            VertexLagMatrix = matrix(1, maxLag,
                                                     length(VertexStatsvec)),
                            dayClass = NA,
                            regMethod = "bayesglm"){
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
            verstats.lag <- matrix(0, nrow = length(Vunion), ncol = nvertexstats)
            rownames(verstats.lag) <- Vunion
            verstats.lag[match(rownames(current.vstats),
                               rownames(verstats.lag)), ] <- current.vstats
            vstats.current <- cbind(vstats.current, verstats.lag)
        }
        net.current <- InputNetwork[[i + maxLag]]

        if(sum(!is.na(dayClass)) > 0) {
            ## we add network level exogenous variable here.
            day.current <- dayClass[i + maxLag]
            xlags.current <- cbind(xlags.current, day.current)
        }
        y.current <- numeric(length(Vunion))
        current.vnames <- network.vertex.names.1(net.current)
        y.current[match(current.vnames, Vunion)] <- 1

        y <- c(y, y.current)
        x.Nets <- rbind(x.Nets, xlags.current)
        x.vstats <- rbind(x.vstats, vstats.current)
    }

    for(i in seq_len(NCOL(x.Nets) - 1)){
        colnames(x.Nets)[i] <- paste0("lag", i, sep = "")
    }
    colnames(x.Nets)[ncol(x.Nets)] <- "Day"
    
    cnames <- numeric(ncol(x.vstats))
    for(i in seq_len(maxLag)){
        for(j in seq_len(nvertexstats)){
            cnames[(i - 1)*nvertexstats + j] <- 
                paste0(vstatNames[j],"Lag",i, sep = ".")
        }
    }
    colnames(x.vstats) <- cnames
    XYdata <- cbind(y, x.Nets, x.vstats)

    
    ## colnames(XYdata) <- NULL
    XYdata <- as.data.frame(XYdata)
    colnames(XYdata)[1] <- "y"
    ## subset
    VertexLagvec <- c(t(VertexLagMatrix))
    if(sum(!is.na(dayClass)) > 0) VertexLagvec <- c(1, VertexLagvec)
    if(maxLag > 1) VertexLagvec <- c(VertexLag, VertexLagvec)
    VertexLagvec <- c(1, VertexLagvec)
    VertexLagvec <- VertexLagvec == 1
    
    VertexRegout <- regEngine(XYdata[, VertexLagvec], method = regMethod)
    
    ## end of vertex only regression ##
    out <- list(VertexFit = VertexRegout,
                VertexStats = XYdata[, VertexLagvec],
                VertexStatsFull = XYdata)
}
