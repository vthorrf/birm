rasch <- function(x, interaction=F, method="LA", Iters=100, Smpl=1000,
                  Thin=1, a.s=0.234, temp=1e-2, tmax=1, algo="GA", seed=666){

  ### Start====
  #require(LaplacesDemon)
  #require(compiler)
  #require(parallel)
  #require(tidyr)
  CPUs = detectCores(all.tests = FALSE, logical = TRUE) - 1
  if(CPUs == 0) CPUs = 1

  ### Convert data to long format====
  lonlong <- gather(data.frame(x), "item", "resp", colnames(x), factor_key=TRUE)
  data_long <- data.frame(ID=rep(1:nrow(x), times=ncol(x)),lonlong)

  ### Choose a model====
  if (interaction == F) {
    ### Rasch model
    # Assemble data list
    if (method == "MAP") {
      mon.names  <- "LL"
    } else { mon.names  <- "LP" }
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)), b=rep(0,ncol(x)) ))
    pos.theta  <- grep("theta", parm.names)
    pos.b      <- grep("b", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n)
      b     <- rnorm(Data$v)
      return(c(theta, b))
    }
    MyData <- list(parm.names=parm.names, mon.names=mon.names,
                   PGF=PGF, X=data_long, n=nrow(x), v=ncol(x),
                   pos.theta=pos.theta, pos.b=pos.b)
    is.data(MyData)

    # Model
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- parm[Data$pos.b]

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=TRUE))
      b.prior     <- sum(dnorm(b    , mean=0, sd=1, log=T))
      Lpp <- theta.prior + b.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      IRF     <- plogis( thetaLL - bLL )
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

  } else if (interaction == T) {
    ### Rasch Interaction model
    # Assemble data list
    mon.names  <- "LP"
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)),
                                      b=rep(0,ncol(x)),
                                      delta=rep(0,1) ))
    pos.theta  <- grep("theta", parm.names)
    pos.b      <- grep("b", parm.names)
    pos.delta  <- grep("delta", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n)
      b     <- rnorm(Data$v)
      delta <- rnorm(1, 0, 1)
      return(c(theta, b, delta))
    }
    MyData <- list(parm.names=parm.names, mon.names=mon.names,
                   PGF=PGF, X=data_long, n=nrow(x), v=ncol(x),
                   pos.theta=pos.theta, pos.b=pos.b, pos.delta=pos.delta)
    is.data(MyData)

    # Model
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- parm[Data$pos.b]
      delta <- parm[Data$pos.delta]

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=T))
      b.prior     <- sum(dnorm(b, mean=0, sd=1, log=T))
      delta.prior <- sum(dnorm(delta, mean=0, sd=1, log=T))
      Lpp <- theta.prior + b.prior + delta.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      deltaLL <- rep(delta, each=Data$n * Data$v)
      IRF     <- plogis( thetaLL + bLL + (deltaLL * thetaLL * bLL) )
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

  } else stop("Unknow model :(")

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
    ## Hit-And-Run Metropolis====
    Iters=Iters; Status=Iters/10; Thin=Thin; A=a.s
    Fit <- LaplacesDemon(Model=Model, Data=MyData,
                         Initial.Values=Initial.Values,
                         Covar=NULL, Iterations=Iters,
                         Status=Status, Thinning=Thin,
                         Algorithm="HARM",
                         Specs=list(alpha.star=A, B=NULL))
  } else if (method=="PMC") {
    ## Population Monte Carlo====
    #Iters=10; Thin=1; Smpl=1000
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
  } else if (method=="MAP") {
    ## Maximum a Posteriori====
    #Iters=100; Smpl=1000
    Iters=Iters; Status=Iters/10
    Fit <- MAP(Model=Model, parm=Initial.Values, Data=MyData, algo=algo,
               maxit=Iters, temp=temp, tmax=tmax, REPORT=Status)
  } else {stop('Unknown optimization method.')}

  ### Results====
  if (method=="MAP") {
    abil = Fit$parm[pos.theta]
    diff = Fit$parm[pos.b]
    FI    = Fit$FI

    if (interaction == F) {
      Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                      'abil'=abil,'diff'=diff,'FitIndexes'=FI)
    } else {
      weight = Fit$parm[pos.delta]
      Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                      'abil'=abil,'diff'=diff,'weight'=weight,'FitIndexes'=FI)
    }
  } else {
    if (method=="PMC") {
      abil = Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
      diff = Fit$Summary[grep("b", rownames(Fit$Summary), fixed=TRUE),1]
    } else {
      abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
      diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1]
    }
    Dev  = Fit$Deviance
    DIC  = list(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)

    if (interaction == F) {
      Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                      'abil'=abil,'diff'=diff,'DIC'=DIC)
    } else {
      if (method=="PMC") {
        weight = Fit$Summary[grep("delta", rownames(Fit$Summary), fixed=TRUE),1]
      } else {
        weight = Fit$Summary1[grep("delta", rownames(Fit$Summary1), fixed=TRUE),1]
      }
      Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                      'abil'=abil,'diff'=diff,'weight'=weight,'DIC'=DIC)
    }
  }

  return(Results)
}
