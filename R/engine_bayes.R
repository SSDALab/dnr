##' Implementation of simulation engine for dynamic networks using smoothing estimates of change statistics.
##' @param start_network Initial list of networks 
##' @param inputcoeff coefficient vector
##' @param ns number of time points for simulation
##' @param model.terms model terms in formula
##' @param model.formula model formula (ergm)
##' @param graph_mode 'digraph' by default
##' @param group group terms
##' @param intercept intercept terms
##' @param exvar extraneous covariates
##' @param maxlag maximum lag
##' @param lagmat lag matrix
##' @param ylag lag vector for network lag terms
##' @param lambda NA
##' @param method 'bayesglm' by default
##' @param alpha.glmnet NA
##' @param paramout T/F parameter estimation is returned.
##' @param Theta = prior probability matrix.
##' @examples
##' \dontrun{
##' startNet <- rdNets[1:50]
##' model.terms=c("triadcensus.003", "triadcensus.012", "triadcensus.102", "triadcensus.021D", "gwesp")
##' model.formula = net~triadcensus(0:3)+gwesp(alpha=0, fixed=FALSE, cutoff=30)-1
##' graph_mode <- 'digraph'
##' group <- 'dnc'
##' alpha.glmnet <- 1
##' method <- 'bayesglm'
##' maxlag <- 3
##' lambda <- NA
##' intercept <- "edges"
##' cdim <- length(model.terms)
##' lagmat <- matrix(sample(c(0,1),(maxlag+1)*cdim,replace = TRUE),ncol = cdim)
##' ylag <- rep(1,maxlag)
##' lagmat[1,] <- rep(0,ncol(lagmat))
##' 
##' out.coef <- paramEdge(input_network = startNet,
##'                 model.terms = model.terms,
##'                 model.formula = model.formula,
##'                 graph_mode='digraph',
##'                 group=group,intercept = intercept,
##'                 exvar=NA,
##'                 maxlag = maxlag,
##'                 lagmat = lagmat,
##'                 ylag = ylag,
##'                 lambda = NA, method='bayesglm',
##'                 alpha.glmnet=1)
##' 
##' 
##' inputcoeff <- out.coef$coef$coef.edge
##' nvertex <- 47 ##find vertex here
##' ns <- 1
##' exvar <- NA
##' for(i in seq_along(startNet)) Theta <- Theta + startNet[[i]][,]
##' Theta <- Theta/length(startNet)
##' Theta <- thresh(Theta)
##' out.bayes <- engineEdgeBayes(start_network=startNet,
##' inputcoeff=inputcoeff,
##' ns=ns,
##' model.terms=model.terms,
##' model.formula=model.formula,
##' graph_mode=graph_mode,
##' group=group,intercept=intercept,
##' exvar=exvar,
##' maxlag=maxlag,
##' lagmat=lagmat,
##' ylag=ylag,
##' lambda = NA, method='bayesglm',
##' alpha.glmnet=alpha.glmnet,
##' Theta = Theta)
##'}
engineEdgeBayes <- function(start_network,inputcoeff,ns,
                         model.terms, model.formula,
                         graph_mode,group,intercept,
                         exvar,
                         maxlag,
                         lagmat,
                         ylag,
                         lambda = NA, method='bayesglm',
                         alpha.glmnet,
                         paramout = TRUE,
                         Theta = NA){
  nnetinput <- length(start_network)
  repfac <- nnetinput-maxlag
  nvertex <- network.size(start_network[[1]])
  nedges <- if(graph_mode=='digraph') nvertex*(nvertex-1)
  out_network <- list()
  coeflist <- list()
  lagmat[1,] <- rep(0,ncol(lagmat))
  for(ncount in 1:ns){
    print(ncount)
    out1 <- paramEdge(start_network,model.terms, model.formula,
                     graph_mode=graph_mode,group,intercept = intercept,
                     exvar=exvar,
                     maxlag = maxlag,
                     lagmat = lagmat,
                     ylag = ylag,
                     lambda = NA, method=method,
                     alpha.glmnet=alpha.glmnet,
                     paramout = paramout)
    inputmpleMat <- as.matrix(out1$mplemat[,-1])
    smmpleMat <- matrix(0,nedges,ncol(inputmpleMat))
    for(i in 1:repfac){
      smmpleMat <- smmpleMat + inputmpleMat[(((i-1)*nedges+1):(i*nedges)),]
    }
    smmpleMat <- smmpleMat/repfac
    inputpred <- smmpleMat%*%inputcoeff
    
    net.current <- start_network[[1]]
    X_t <- ungvectorize(inputpred,nvertex,graph_mode)
    Pobs <- ilogit(X_t)
    alpha.prior <- Theta/(1-Theta)
    P.post <- (Pobs + alpha.prior)/(alpha.prior/Theta + 1)
    net.current %n% "X" <- arm::logit(P.post)
    net.current <- simulate(ergm(net.current ~ edgecov("X")))
    out_network[[ncount]] <- net.current
    for(i in 1:(nnetinput-1)){
      start_network[[i]] <- start_network[[i+1]]
    }
    start_network[[nnetinput]] <- net.current
    coeflist[[ncount]] <- out1$coef$coef
  }
  coefmat <- do.call(rbind,coeflist)
  return(list(out_network=out_network,coefmat=coefmat))
}
