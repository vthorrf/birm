plm <- function(x, p=1, scaling=1.7, method="VB", Iters=500, Smpl=500,
                Status=100, Thin=1, Ad=500, delta=.6, seed=666){

  ### Start====
  require(LaplacesDemon)
  require(compiler)
  require(parallel)
  require(tidyr)
  CPUs = detectCores(all.tests = FALSE, logical = TRUE) - 1
  if(CPUs == 0) CPUs = 1

  ### Convert data to long format====
  lonlong <- gather(data.frame(x), item, resp, colnames(x), factor_key=TRUE)
  data_long <- data.frame(ID=rep(1:nrow(x), times=ncol(x)),lonlong)

  ### Choose Model====
  if (p == 1) {
    ## One-parameter logistic====
    # Assemble data list
    mon.names  <- "LP"
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)), b=rep(0,ncol(x)),
                                      Ds=rep(0,1) ))
    pos.theta  <- grep("theta", parm.names)
    pos.b      <- grep("b", parm.names)
    pos.Ds    <- grep("Ds", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n)
      b     <- rnorm(Data$v)
      Ds    <- rnorm(1)
      return(c(theta, b, Ds))
    }
    MyData <- list(parm.names=parm.names, mon.names=mon.names,
                   PGF=PGF, X=data_long, n=nrow(x), v=ncol(x), scaling=scaling,
                   pos.theta=pos.theta, pos.b=pos.b, pos.Ds=pos.Ds)
    is.data(MyData)

    # Model
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- parm[Data$pos.b]
      Ds    <- parm[Data$pos.Ds]

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=TRUE))
      b.prior     <- sum(dnorm(b, mean=0, sd=1, log=T))
      Ds.prior    <- sum(dnorm(Ds, mean=0, sd=1, log=T))
      Lpp <- theta.prior + b.prior + Ds.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      DLL     <- rep(Ds   , each=Data$v * Data$n)
      IRF     <- plogis( Data$scaling * DLL * (thetaLL - bLL) )
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- rbinom(length(IRF), size=1, prob=IRF)
      ### Output
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
      return(Modelout)
    }
    Model <- compiler::cmpfun(Model)
    Initial.Values <- GIV(Model, MyData, PGF=T)
    is.model(Model, Initial.Values, MyData)
    is.bayesian(Model, Initial.Values, MyData)

  } else if (p == 2) {
    ## Two-parameter logistic====
    # Assemble data list
    mon.names  <- "LP"
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)), b=rep(0,ncol(x)),
                                      Ds=rep(0,ncol(x)) ))
    pos.theta  <- grep("theta", parm.names)
    pos.b      <- grep("b", parm.names)
    pos.Ds    <- grep("Ds", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n)
      b     <- rnorm(Data$v)
      Ds    <- rnorm(Data$v)
      return(c(theta, b, Ds))
    }
    MyData <- list(parm.names=parm.names, mon.names=mon.names,
                   PGF=PGF, X=data_long, n=nrow(x), v=ncol(x),
                   pos.theta=pos.theta, pos.b=pos.b, pos.Ds=pos.Ds)
    is.data(MyData)

    # Model
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- parm[Data$pos.b]
      Ds    <- parm[Data$pos.Ds]

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=TRUE))
      b.prior     <- sum(dnorm(b, mean=0, sd=1, log=T))
      Ds.prior    <- sum(dnorm(Ds, mean=0, sd=1, log=T))
      Lpp <- theta.prior + b.prior + Ds.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      DLL     <- rep(Ds   , each=Data$n)
      IRF     <- plogis( DLL * (thetaLL - bLL) )
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- rbinom(length(IRF), size=1, prob=IRF)
      ### Output
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
      return(Modelout)
    }
    Model <- compiler::cmpfun(Model)
    Initial.Values <- GIV(Model, MyData, PGF=T)
    is.model(Model, Initial.Values, MyData)
    is.bayesian(Model, Initial.Values, MyData)

  } else if (p == 3) {
    ## Three-parameter logistic====
    # Assemble data list
    mon.names  <- "LP"
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)), b=rep(0,ncol(x)),
                                      Ds=rep(0,ncol(x)), c=rep(0,ncol(x)) ))
    pos.theta  <- grep("theta", parm.names)
    pos.b      <- grep("b", parm.names)
    pos.Ds     <- grep("Ds", parm.names)
    pos.c      <- grep("c", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n)
      b     <- rnorm(Data$v)
      Ds    <- rnorm(Data$v)
      c     <- rbeta(Data$v, 1, 1)
      return(c(theta, b, Ds, c))
    }
    MyData <- list(parm.names=parm.names, mon.names=mon.names,
                   PGF=PGF, X=data_long, n=nrow(x), v=ncol(x),
                   pos.theta=pos.theta, pos.b=pos.b, pos.c=pos.c, pos.Ds=pos.Ds)
    is.data(MyData)

    # Model
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- parm[Data$pos.b]
      Ds    <- parm[Data$pos.Ds]
      c     <- interval( parm[Data$pos.c], 1e-100, (1 - 1e-100) )
      parm[Data$pos.c] <- c

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=TRUE))
      b.prior     <- sum(dnorm(b, mean=0, sd=1, log=T))
      Ds.prior    <- sum(dnorm(Ds, mean=0, sd=1, log=T))
      c.prior     <- sum(dbeta(c, 1, 1, log=T))
      Lpp <- theta.prior + b.prior + Ds.prior + c.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      DLL     <- rep(Ds   , each=Data$n)
      cLL     <- rep(c    , each=Data$n)
      IRF     <- cLL + ( (1 - cLL) / (1 + exp(-DLL * ( thetaLL - bLL ))) )
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- rbinom(length(IRF), size=1, prob=IRF)
      ### Output
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
      return(Modelout)
    }
    Model <- compiler::cmpfun(Model)
    Initial.Values <- GIV(Model, MyData, PGF=T)
    is.model(Model, Initial.Values, MyData)
    is.bayesian(Model, Initial.Values, MyData)

  } else if (p == 4) {
    ## Four-parameter logistic====
    # Assemble data list
    mon.names  <- "LP"
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)), b=rep(0,ncol(x)),
                                      Ds=rep(0,ncol(x)), c=rep(0,ncol(x)),
                                      UA=rep(0,ncol(x)) ))
    pos.theta  <- grep("theta", parm.names)
    pos.b      <- grep("b", parm.names)
    pos.Ds     <- grep("Ds", parm.names)
    pos.c      <- grep("c", parm.names)
    pos.UA     <- grep("UA", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n)
      b     <- rnorm(Data$v)
      Ds    <- rnorm(Data$v)
      c     <- rbeta(Data$v, 1, 1)
      UA    <- rbeta(Data$v, 1, 1)
      return(c(theta, b, Ds, c, UA))
    }
    MyData <- list(parm.names=parm.names, mon.names=mon.names,
                   PGF=PGF, X=data_long, n=nrow(x), v=ncol(x),
                   pos.theta=pos.theta, pos.b=pos.b, pos.c=pos.c,
                   pos.UA=pos.UA, pos.Ds=pos.Ds)
    is.data(MyData)

    # Model
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- parm[Data$pos.b]
      Ds    <- parm[Data$pos.Ds]
      c     <- interval( parm[Data$pos.c], 1e-100, (1 - 1e-100) )
      parm[Data$pos.c] <- c
      UA    <- interval( parm[Data$pos.UA], 1e-100, (1 - 1e-100) )
      parm[Data$pos.UA] <- UA

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=TRUE))
      b.prior     <- sum(dnorm(b, mean=0, sd=1, log=T))
      Ds.prior    <- sum(dnorm(Ds, mean=0, sd=1, log=T))
      c.prior     <- sum(dbeta(c, 1, 1, log=T))
      UA.prior     <- sum(dbeta(UA, 1, 1, log=T))
      Lpp <- theta.prior + b.prior + Ds.prior + c.prior + UA.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      DLL     <- rep(Ds   , each=Data$n)
      cLL     <- rep(c    , each=Data$n)
      UALL     <- rep(UA   , each=Data$n)
      IRF     <- cLL + ( (UALL - cLL) / (1 + exp(-DLL * ( thetaLL - bLL ))) )
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- rbinom(length(IRF), size=1, prob=IRF)
      ### Output
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
      return(Modelout)
    }
    Model <- compiler::cmpfun(Model)
    Initial.Values <- GIV(Model, MyData, PGF=T)
    is.model(Model, Initial.Values, MyData)
    is.bayesian(Model, Initial.Values, MyData)

  } else warning("Unknow model :(")

  ### Run!====
  set.seed(seed)
  if (method=="VB") {
    ## Variational Bayes====
    #Iters=1000; Samples=1000
    Iters=Iters; Smpl=Smpl
    Fit <- VariationalBayes(Model=Model, parm=Initial.Values, Data=MyData,
                            Covar=NULL, Interval=1e-6, Iterations=Iters,
                            Method="Salimans2", Samples=Smpl, sir=TRUE,
                            Stop.Tolerance=1e-5, CPUs=CPUs, Type="PSOCK")
  } else if (method=="LA") {
    ## Laplace Approximation====
    #Iters=100; Samples=1000
    Iters=Iters; Smpl=Smpl
    Fit <- LaplaceApproximation(Model, parm=Initial.Values, Data=MyData,
                                Interval=1e-6, Iterations=Iters,
                                Method="SPG", Samples=Smpl, sir=TRUE,
                                CovEst="Identity", Stop.Tolerance=1e-5,
                                CPUs=CPUs, Type="PSOCK")
  } else if (method=="MCMC") {
    ## No-U-Turn Sampler====
    #Iters=10000; Status=100; Thin=10; Ad=500; delta=.6
    Iters=Iters; Status=Status; Thin=Thin; Ad=Ad; delta=delta
    Fit <- LaplacesDemon(Model=Model, Data=MyData,
                         Initial.Values=Initial.Values,
                         Covar=NULL, Iterations=Iters,Status=Status,
                         Thinning=Thin, Algorithm="NUTS",
                         Specs=list(A=Ad,delta=delta,epsilon=NULL,Lmax=Inf))
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
  } else {stop('Unknown estimation method.')}

  ### Results====
  if (p == 1) {
    ## One-parameter logistic====
    abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
    diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1]
    disc = Fit$Summary1[grep("Ds", rownames(Fit$Summary1), fixed=TRUE),1]
    Dev  = Fit$Deviance
    DIC  = list(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)

    Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                    'abil'=abil,'diff'=diff,"disc"=disc,'DIC'=DIC)
  } else if (p == 2) {
    ## Two-parameter logistic====
    abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
    diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1]
    disc = Fit$Summary1[grep("Ds", rownames(Fit$Summary1), fixed=TRUE),1]
    Dev  = Fit$Deviance
    DIC  = list(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)

    Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                    'abil'=abil,'diff'=diff,"disc"=disc,'DIC'=DIC)
  } else if (p == 3) {
    ## Three-parameter logistic====
    abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
    diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1]
    disc = Fit$Summary1[grep("Ds", rownames(Fit$Summary1), fixed=TRUE),1]
    gues = Fit$Summary1[grep("c", rownames(Fit$Summary1), fixed=TRUE),1]
    Dev  = Fit$Deviance
    DIC  = list(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)

    Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                    'abil'=abil,'diff'=diff,"disc"=disc,"gues"=gues,'DIC'=DIC)
  } else if (p == 4) {
    ## Four-parameter logistic====
    abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
    diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1]
    disc = Fit$Summary1[grep("Ds", rownames(Fit$Summary1), fixed=TRUE),1]
    gues = Fit$Summary1[grep("c", rownames(Fit$Summary1), fixed=TRUE),1]
    UpAs = Fit$Summary1[grep("UA", rownames(Fit$Summary1), fixed=TRUE),1]
    Dev  = Fit$Deviance
    DIC  = list(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)

    Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                    'abil'=abil,'diff'=diff,"disc"=disc,"gues"=gues,
                    "UpAs"=UpAs,'DIC'=DIC)

  } else warning("Can't return any result :P")

  return(Results)
}
