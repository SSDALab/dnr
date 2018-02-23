deg <- function(z, gmode){
    out <- cbind(sna::degree(z, gmode = gmode, cmode = "freeman"),
                 sna::degree(z, gmode = gmode, cmode = "indegree"),
                 sna::degree(z, gmode = gmode, cmode = "outdegree"))
    return(out)
}

eigenCentrality <- function(z, gmode){
    out <- igraph::eigen_centrality(igraph::graph_from_adjacency_matrix(z[,]))$vector
    out
}

BetCentrality <- function(z, gmode){
    out <- sna::betweenness(z, gmode = gmode)
    return(out)
}

InfoCentrality <- function(z, gmode){
    out <- tryCatch(sna::infocent(z, gmode = gmode, tol = 1e-40),
                    error = function(err){return(0)})
    return(out)
}

closenessCent <- function(z, gmode){
    out <- sna::closeness(z, gmode = gmode)
    return(out)
}

logkCycle <- function(z, gmode) {
    x <- sna::kcycle.census(z, mode = gmode, maxlen = 4)
    d <- log(x$cycle.count[ 3, 2:(network.size(z)+1)] + 1)
    return(d)
}

## log size
logSize <- function(z, gmode){
    if(!network::is.network(z)) {
        if(is.na(z)) return (NULL)
    } else if(sum(c(z[ , ])) == 0) {
        return (NULL)
    } else {
        nz <- network.size.1(z)
        d <- log(nz + 1e-10)
        d <- rep(d, nz)
    }
    return(d)
}


vertexstats <- function(z, gmode){
    if(!network::is.network(z)) {
        if(is.na(z)) return (NULL)
    } else if(sum(c(z[ , ])) == 0) {
        return (NULL)
    } else {
        out <- cbind(deg(z, gmode), eigenCentrality(z, gmode),
                     BetCentrality(z, gmode), InfoCentrality(z, gmode),
                     closenessCent(z, gmode), logkCycle(z, gmode),
                     logSize(z, gmode))
        colnames(out) <- c("degree", "InDeg", "OutDeg", "Eigencent",
                           "BetCent", "InfoCent", "CloseCent",
                           "LogCyc", "LogSize")
        rownames(out) <- network.vertex.names(z)
    }
    return(out)
}
nvertexstats <- 9
vstatNames <- c(
    "Degree",
    "InDegree",
    "OutDegree",
    "EigenCentrality",
    "BetweenCentrality",
    "InfoCentrality",
    "CloseCentrality",
    "LogCycle",
    "LogSize")
