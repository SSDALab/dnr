##' @title Parameter estimation for static vertex case.
##' @description Parameter estimation for the static vertex case.
##' @param input_network Input network.
##' @param model.terms model terms, must be ERGM terms expanded.
##' @param model.formula ERGM formula for each time point.
##' @param graph_mode 'digraph' by default for bidirectional.
##' @param group grouping covariates for vertices.
##' @param intercept intercept terms.
##' @param exvar Extraneous variables
##' @param maxlag maximum lag.
##' @param lagmat Matrix of dimension (maxlag+1)x(length(model.terms))
##' @param ylag lag vectors of length=maxlag.
##' @param lambda NA
##' @param method Regression method, default is 'bayesglm'
##' @param alpha.glmnet if regularization is used. not needed for bayesglm.
##' @param paramout TRUE by default. if parameters are needed.
##' @return list with elements:
##'   coef: coefficients
##'   mplematfull: full martix of change statistics
##'   mplemat: subset of matrix of change statistics
##' @author Abhirup
##' @export
##' @examples 
##' input_network=rdNets[1:6]
##' model.terms=c("triadcensus.003", "triadcensus.012", "triadcensus.102", "triadcensus.021D", "gwesp");
##' model.formula = net~triadcensus(0:3)+gwesp(decay=0, fixed=FALSE, cutoff=30)-1;
##' graph_mode='digraph';
##' group='dnc';
##' alpha.glmnet=1
##' directed=TRUE;
##' method <- 'bayesglm'
##' maxlag <- 3
##' lambda=NA
##' intercept = c("edges")
##' cdim <- length(model.terms)
##' lagmat <- matrix(sample(c(0,1),(maxlag+1)*cdim,replace = TRUE),ncol = cdim)
##' ylag <- rep(1,maxlag)
##' exvar <- NA
##' out <- paramEdge(input_network,model.terms, model.formula,
##'                 graph_mode='digraph',group,intercept = c("edges"),exvar=NA,
##'                 maxlag = 3,
##'                 lagmat = matrix(sample(c(0,1),(maxlag+1)*cdim,
##'                                        replace = TRUE),ncol = cdim),
##'                 ylag = rep(1,maxlag),
##'                 lambda = NA, method='bayesglm',
##'                 alpha.glmnet=1)
##' 

paramEdge <- function(input_network,model.terms, model.formula,
                     graph_mode='digraph',group,intercept = c("edges"),exvar=NA,
                     maxlag = 3,
                     lagmat = matrix(
                       sample(c(0,1),
                              (maxlag+1)*length(model.terms),
                                            replace = T),ncol = length(model.terms)),
                     ylag = rep(1,maxlag),
                     lambda = NA, method='glmnet',
                     alpha.glmnet=1,paramout=TRUE){



  gengroup <- function(input_network,group,net1){
    Vmax <- input_network;

    Vmax.label <- network::get.vertex.attribute(Vmax[[1]], attrname = group);

    Vmax <- network.vertex.names(Vmax[[1]])

    foo.index <- which(Vmax%in%network.vertex.names(net1)==T)

    if (is.null(group)){
      grouping <- c(rep(1, floor(length(Vmax)/2)), rep(2, length(Vmax)-floor(length(Vmax)/2)))
    } else {
      grouping <- Vmax.label
    }
    grouping.sub <- grouping[foo.index]

    if(directed){
      grouping.perm <- expand.grid(unique(grouping.sub),unique(grouping.sub))
      grouping.perm$indicator <- seq_along(grouping.perm[,1])


      grouping.perm.full <- expand.grid(unique(grouping),unique(grouping))
      grouping.perm.full$indicator <- seq_along(grouping.perm.full[,1])

      foo <- matrix(NA,ncol=length(grouping.sub), nrow=length(grouping.sub))


      grouping.terms <- sapply(1:nrow(grouping.perm.full),function(i) paste(group,grouping.perm.full[i,1],grouping.perm.full[i,2],sep=''))
      for (i in 1:length(grouping.sub)){
        for(j in 1:length(grouping.sub)){
          foo[i,j]<-grouping.perm[,3] [which(grouping.perm[,1]==grouping.sub[i]&grouping.perm[,2]==grouping.sub[j])]
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


      grouping.perm$indicator <- seq_along(grouping.perm[,1])

      if(length(unique(grouping))>1){
        grouping.perm.full <- data.frame(rbind(t(sapply(unique(grouping), function(x)c(x,x))), t(combn(unique(grouping),2))))
      } else{
        grouping.perm.full <- data.frame(t(c(unique(grouping),unique(grouping))))
      }


      grouping.perm.full$indicator <- seq_along(grouping.perm.full[,1])

      foo <- matrix(NA,ncol=length(grouping.sub), nrow=length(grouping.sub))


      grouping.terms = sapply(1:nrow(grouping.perm.full),function(i) paste(group,grouping.perm.full[i,1],grouping.perm.full[i,2],sep=''))


      for (i in 1:length(grouping.sub)){
        for(j in 1:length(grouping.sub)){
          ind1 <- apply(grouping.perm[,1:2],1,function(x) all(x == c(grouping.sub[i], grouping.sub[j]))) + apply(grouping.perm[,1:2],1,function(x) all(x == c(grouping.sub[j], grouping.sub[i])));
          ind1[which(ind1>0)] <- 1;
          foo[i,j]<-grouping.perm[,3][which(ind1!=0)];

          #foo[i,j]<-grouping.perm[,3] [which((grouping.perm[,1]==grouping.sub[i]&grouping.perm[,2]==grouping.sub[j])|(grouping.perm[,1]==grouping.sub[j]&grouping.perm[,2]==grouping.sub[i]))]
          #print(i,j)
        }
      }

      for(i in grouping.perm.full$indicator){
        pos <- 1
        envir = as.environment(pos)
        assign(grouping.terms[i], (foo==grouping.perm.full$indicator[i])*1, envir = envir)
      }
    }

    grouping.edgecov <- sapply(1:nrow(grouping.perm.full), function(x) paste0('edgecov(',grouping.terms[x],')',sep=''))
    return(grouping.edgecov)
  }

  ########## End of gengroup ####



  genformula1 <- function(model.formula,grouping.edgecov,netname="net1"){
    term_formula = terms(model.formula);
    term_formula = attr(term_formula, 'term.labels');
    formula1 = as.formula(paste(paste0(netname,"~",sep=""), paste(c(grouping.edgecov,term_formula,"edgecov(m)"), collapse= "+")))
    return(formula1)
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


  excterm <- function(model.terms,exterm = 'logsize'){
    find.logsize = grep(exterm,model.terms);
    if (length(find.logsize) !=0 ){
      model.terms.new = model.terms[-find.logsize]
    } else{
      model.terms.new = model.terms
    }
    return(model.terms.new)
  }


  excform <- function(model.formula,exterm = 'logsize'){
    term_formula <- terms(model.formula)
    find.logsize = grep('logsize',term_formula);
    if (length(find.logsize) !=0 ){
      term_formula.new = term_formula[-find.logsize]
    } else{
      term_formula.new = term_formula
    }
    return(term_formula.new)
  }

  # term_formula = terms(model.formula);
  # term_formula = attr(term_formula, 'term.labels');
  #
  # formula = as.formula(paste("net1~ ", paste(c(grouping.edgecov,term_formula,"edgecov(m)"), collapse= "+")))
  #
  # mplemat  =  ergmMPLE(formula, output="matrix");
  # edgeLag = mplemat$predictor[,1:(dim(mplemat$predictor)[2]-1)];
  #

  ####### End of Formula manipulation ####

  #initialization
  k <- maxlag + 1
  matout <- NULL
  n.network <- length(input_network)
  #Legacy code from Yang
  if (graph_mode=='digraph'){
    gmode='digraph';
    mode='digraph'
    directed=TRUE;
    cmode="directed";
  } else{
    gmode='graph';
    mode='graph'
    directed=FALSE;
  }
  #sliding window
  for(i in 1:(n.network-k+1)){
    net.window <- input_network[i:(i+k-1)]
    net.current <- net.window[[k]]
    #fit intercept
    #construct the intercept
    size <- network.size(net.current);
    m <- matrix(1:(size*size),nrow=size,ncol=size);
    if(is.na(group)){
      grouping.edgecov <- NA
    }else{
      grouping.edgecov <- gengroup(input_network,group,net.current)
    }

    #Construction of formula##

    if(is.na(intercept)&&is.na(group)){
      csmodel <- NULL
    } else{
      formula <- genintercept(intercept,grouping.edgecov,
                              netname="net.current")
      #TODO: Add support for exogenous variable later
      if(!is.na(exvar)){
        formula <- as.formula(paste(c(formula,paste0("edgecov(",exvar,")",sep="")),collapse = "+"))
      }
      mplemat  <-  ergmMPLE(formula, output="matrix")
      csintercept <- as.matrix(mplemat$predictor[,1:(ncol(mplemat$predictor)-1)])
      colnames(csintercept) <- colnames(mplemat$predictor)[1:(ncol(mplemat$predictor)-1)]
      csmodel <- csintercept
    }
    #fit model terms
    #Comment: We always count down!
    if(sum(!is.na(model.formula)) > 0){
      for(j in k:1){
        formula <- genformula(model.formula,netname = "net.window[[j]]")
        mplemat.tmp  <-  ergmMPLE(formula, output="matrix");
        edgeLag.tmp <- mplemat.tmp$predictor[,1:(ncol(mplemat.tmp$predictor)-1)]
        csmodel <- cbind(csmodel, edgeLag.tmp);
      }
    }
    #fit the lag terms
    #Comment: Counting down.
    if(k > 1){
      if(graph_mode=='digraph'){
        lagstats <- matrix(0, ncol=k-1, nrow=size*(size-1))
      } else{
        lagstats <- matrix(0, ncol=k-1, nrow=size*(size-1)/2)
      }
      lagnames <- rep("lag",(k-1))
      for (j in (k-1):1){
        lagstats[,j] <- sna::gvectorize(net.window[[j]][,],mode=graph_mode, censor.as.na=F)
        lagnames[j] <- paste("lag",j,sep = "")
      }
      colnames(lagstats) <- lagnames
      #response
      y <- sna::gvectorize(net.current[,],mode=graph_mode, censor.as.na=FALSE)
      mat <- data.frame(cbind(y,csmodel,lagstats))
    } else {
      y <- sna::gvectorize(net.current[,],mode=graph_mode, censor.as.na=FALSE)
      mat <- data.frame(cbind(y,csmodel))
    }
    matout <- rbind(matout,mat)
  }
  #subsetting
  lagvec <- rep(1,NCOL(csintercept))
  if(sum(!is.na(model.formula)) > 0) lagvec <- c(lagvec,t(lagmat))
  if(k > 1) lagvec <- c(lagvec,ylag)
  lagvec <- c(1,lagvec)
  lagvec <- lagvec==1
  #regression
  if(paramout){
    out <- regEngine(matout[,lagvec],method)
  } else out <- NULL
  #results
  return(list(coef=out,mplematfull=matout,mplemat=matout[,lagvec]))
}
