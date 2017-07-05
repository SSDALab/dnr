###################################################
## Implementation of comparison metrics for performance evaluation of the
## simulation methods.
###################################################

##' Calculate number of triangles of a network
##'
##' This function calculates the number of triangles in a network given an
##'     adjacency matrix. We use igraph for this.
##' @title ntriangles
##' @param x square matrix (adjacency matrix)
##' @return scaler, number of triangles
##' @author Abhirup
##' @export
##' @examples
##' ntriangles(beach[[1]][, ])


ntriangles <- function(x){
    tmp <- matrix(igraph::triangles(igraph::graph_from_adjacency_matrix(x)),
                  nrow = 3)
    ncol(tmp)
}


##' Calculates the cluster coefficient from a network
##'
##' Given a network in the form of adjacency matrix, this calculates the cluster
##'     coefficient. For a definition of cluster coefficient, please refer to
##'     the igraph documentation.
##' @title clustCoef
##' @param x adjacency matrix
##' @return scaler
##' @author Abhirup
##' @export
##' @examples
##' clustCoef(beach[[1]][, ])

clustCoef <- function(x){
    igraph::transitivity(igraph::graph_from_adjacency_matrix(x))
}

##' Calculates the degree of each vertices.
##'
##' Given a network as adjacency matrix, calculate degree stats for each vertex.
##' @title vdegree
##' @param x Adjacency matrix.
##' @return vector of length number of vertices.
##' @author Abhirup
##' @export
##' @examples
##' vdegree(beach[[1]][, ])

vdegree <- function(x){
    igraph::degree(igraph::graph_from_adjacency_matrix(x))
}

##' Calculate the expectation of degree distribution of network
##'
##' Given a network in adjacency matrix form, this calculates the expected
##'     degree statsitic using igraph degree distribution function.
##' @title expdeg
##' @param x adjacency matrix
##' @return scaler
##' @author Abhirup
##' @export
##' @examples
##' expdeg(beach[[1]][, ])

expdeg <- function(x){
    tmp <- igraph::degree_distribution(igraph::graph_from_adjacency_matrix(x))
    sum(tmp*seq_along(tmp))
}
