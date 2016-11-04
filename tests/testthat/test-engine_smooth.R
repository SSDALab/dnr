context("engine_smooth")
start_network <- input_network
inputcoeff <- out$coef$coef.edge
nvertex <- 47
ns <- 10
exvar <- NA
tmp <- engine_smooth(start_network=start_network,inputcoeff=inputcoeff,ns=ns,
                     model.terms=model.terms, model.formula=model.formula,
                     graph_mode=graph_mode,group=group,intercept=intercept,
                     exvar=exvar,
                     maxlag=maxlag,
                     lagmat=lagmat,
                     ylag=ylag,
                     lambda = NA, method='bayesglm',
                     alpha.glmnet=alpha.glmnet)

test_that("engine_smooth returns a list of length 2",{
  expect_equal(length(tmp), 2)
})

## TODO@abhirup: Add more tests