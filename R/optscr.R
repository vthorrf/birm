optscr <- function(x, levels=NULL, knot=NULL, method="VB",
                   Iters=500, Smpl=1000, Thin=1, A=500, seed=666){

  ### Start====
  require(LaplacesDemon)
  require(compiler)
  require(parallel)
  require(tidyr)
  CPUs = detectCores(all.tests = FALSE, logical = TRUE) - 1
  if(CPUs == 0) CPUs = 1

  ### Convert data to long format====
  if (is.null(levels)) levels <- max(x) + 1
  if (is.null(knot)) knot <- seq(0,1, len=if(ncol(x) > 10) 10 else ncol(x))
  base_score <- as.vector(scale(rowMeans(x), center=(levels - 1) * .5))
  var_bscore <- (pnorm(base_score) * (1 - pnorm(base_score)) ) * (levels - 1)
  basis <- length(knot)
  lonlong <- gather(data.frame(x), item, resp, colnames(x), factor_key=TRUE)
  data_long <- data.frame(ID=rep(1:nrow(x), times=ncol(x)),lonlong)

  ### Assemble data list====
  mon.names  <- "LP"
  parm.names <- as.parm.names(list( theta=rep(0,nrow(x)), b=rep(0,ncol(x) * basis) ))
  pos.theta  <- grep("theta", parm.names)
  pos.b      <- grep("b", parm.names)
  PGF <- function(Data) {
    theta <- rnorm(Data$n)
    b     <- rlaplace(Data$basis * Data$v)
    return(c(theta, b))
  }
  MyData <- list(parm.names=parm.names, mon.names=mon.names, levels=levels,
                 PGF=PGF, X=data_long, n=nrow(x), v=ncol(x), basis=basis,
                 base_score=base_score, var_bscore=var_bscore,
                 pos.theta=pos.theta, pos.b=pos.b, knot=knot)
  is.data(MyData)

  ### Model====
  Model <- function(parm, Data){

    ## Prior parameters
    theta <- parm[Data$pos.theta]
    b     <- parm[Data$pos.b]

    ### Log-Priors
    theta.prior <- sum(dnorm(theta, mean=Data$base_score, sd=sqrt(Data$var_bscore), log=T))
    b.prior     <- sum(dlaplace(b, location=0, scale=1, log=T))
    Lpp <- theta.prior + b.prior

    ### Log-Likelihood
    thetaLL <- pnorm( rep(theta, times=Data$v) )
    BsL     <- rep(b    , each=Data$n)
    bLL     <- matrix(BsL , nrow=nrow(Data$X), ncol=Data$basis)
    kLL     <- t(matrix(knot, nrow=Data$basis, ncol=nrow(Data$X)))
    IRF     <- plogis( rowSums((((thetaLL < kLL) * -2) + 1) * bLL) )
    LL      <- sum( dbinom(Data$X[,3], size=(Data$levels - 1), prob=IRF, log=T) )

    ### Log-Posterior
    LP <- LL + Lpp
    ### Estimates
    yhat <- rbinom(length(IRF), size=(Data$levels - 1), prob=IRF)
    ### Output
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
    return(Modelout)
  }
  Model <- compiler::cmpfun(Model)
  Initial.Values <- GIV(Model, MyData, PGF=T)
  is.model(Model, Initial.Values, MyData)
  is.bayesian(Model, Initial.Values, MyData)

  ### Run!====
  set.seed(seed)
  if (method=="VB") {
    ## Variational Bayes====
    #Iters=1000; Smpl=1000
    Iters=Iters; Smpl=Smpl
    Fit <- VariationalBayes(Model=Model, parm=Initial.Values, Data=MyData,
                            Covar=NULL, Interval=1e-6, Iterations=Iters,
                            Method="Salimans2", Samples=Smpl, sir=TRUE,
                            Stop.Tolerance=1e-5, CPUs=CPUs, Type="PSOCK")
  } else if (method=="LA") {
    ## Laplace Approximation====
    #Iters=100; Smpl=1000
    Iters=Iters; Smpl=Smpl
    Fit <- LaplaceApproximation(Model, parm=Initial.Values, Data=MyData,
                                Interval=1e-6, Iterations=Iters,
                                Method="SPG", Samples=Smpl, sir=TRUE,
                                CovEst="Identity", Stop.Tolerance=1e-5,
                                CPUs=CPUs, Type="PSOCK")
  } else if (method=="MCMC") {
    ## No-U-Turn Sampler====
    #Iters=10000; Status=100; Thin=10; Ad=500; delta=.6
    Iters=Iters; Status=Iters/10; Thin=Thin; Ad=A
    Fit <- LaplacesDemon(Model=Model, Data=MyData,
                         Initial.Values=Initial.Values,
                         Covar=NULL, Iterations=Iters,Status=Status,
                         Thinning=Thin, Algorithm="NUTS",
                         Specs=list(A=Ad,delta=0.6,epsilon=NULL,Lmax=Inf))
  } else if (method=="PMC") {
    ## Population Monte Carlo====
    #Iters=10; Thin=11; Smpl=1000
    Iters=Iters; Smpl=Smpl; Thin=Thin
    Fit <- PMC(Model=Model, Data=MyData, Initial.Values=Initial.Values,
               Covar=NULL, Iterations=Iters, Thinning=Thin, alpha=NULL,
               M=2, N=Smpl, nu=1e3, CPUs=CPUs, Type="PSOCK")
  } else if (method=="IQ") {
    ## Iterative Quadrature====
    #Iters=100; Smpl=1000
    Iters=Iters; Smpl=Smpl
    Fit <- IterativeQuadrature(Model=Model, parm=Initial.Values,
                               Data=MyData, Covar=NULL,
                               Iterations=Iters, Algorithm="CAGH",
                               Specs=list(N=3, Nmax=10, Packages=NULL,
                                          Dyn.libs=NULL),
                               Samples=Smpl, sir=T,
                               Stop.Tolerance=c(1e-5,1e-15),
                               Type="PSOCK", CPUs=CPUs)
  } else {stop('Unknown optimization method.')}

  ### Results====
  abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
  Base = matrix(Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1], nrow=ncol(x))
  rownames(Base) = colnames(x)
  colnames(Base) = paste("Beta",1:basis,sep="_")
  ICC <- lapply(1:nrow(Base), function(g) {
    Pr <- plogis(
      as.vector(crossprod(Base[g,],
                          t((((plogis(seq(-3,3,len=100)) <
                                 t(matrix(knot, nrow=basis, ncol=100))) * -2) +1)) )))
    theta <- seq(-3,3,len=100)
    return(data.frame(theta, Pr))
  })
  Dev  = Fit$Deviance
  DIC  = list(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)

  Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                  'abil'=abil,'icc'=ICC,'DIC'=DIC)
  return(Results)
}
