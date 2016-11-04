context("paramest")

input_network=rdNets[1:6];
model.terms=c("triadcensus.003", "triadcensus.012", "triadcensus.102", "triadcensus.021D", "gwesp");
model.formula = net~triadcensus(0:3)+gwesp(alpha=0, fixed=FALSE, cutoff=30)-1;
graph_mode='digraph';
group='dnc';
alpha.glmnet=1
directed=TRUE;
method <- 'bayesglm'
maxlag <- 3
lambda=NA
intercept = c("edges")
cdim <- length(model.terms)
lagmat <- matrix(sample(c(0,1),(maxlag+1)*cdim,replace = T),ncol = cdim)
ylag <- rep(1,maxlag)
lagmat[1,] <- rep(0,ncol(lagmat))
out <- paramest(input_network,model.terms, model.formula,
                graph_mode='digraph',group,intercept = c("edges"),exvar=NA,
                maxlag = 3,
                lagmat = lagmat,
                ylag = rep(1,maxlag),
                lambda = NA, method='bayesglm',
                alpha.glmnet=1)

test_that("Paramest returns a list of length 3",{
  expect_equal(length(out), 3)
})