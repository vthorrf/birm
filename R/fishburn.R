fishburn <- function(x, method="LA", model=1, Iters=100, Smpl=1000,
                     Thin=1, a.s=0.234, temp=1e-2, tmax=NULL,
                     algo="GA", seed=666, Interval=1e-8){

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

  ### Assemble data list====
  if (method == "MAP") {
    mon.names  <- "LL"
  } else { mon.names  <- "LP" }
  parm.names <- as.parm.names(list( theta=rep(0,nrow(x)), b=rep(0,ncol(x)) ))
  pos.theta  <- grep("theta", parm.names)
  pos.b      <- grep("b", parm.names)
  PGF <- function(Data) {
    theta <- rnorm(Data$n, mean=0, sd=1)
    b     <- rnorm(Data$v, mean=0, sd=1)
    return(c(theta, b))
  }
  MyData <- list(parm.names=parm.names, mon.names=mon.names,
                 PGF=PGF, X=data_long, n=nrow(x), v=ncol(x),
                 pos.theta=pos.theta, pos.b=pos.b)
  is.data(MyData)

  ### Model====
  if (model == 1) {
    ############################## 01 - Additive model
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- parm[Data$pos.b]

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=T))
      b.prior     <- sum(dnorm(b    , mean=0, sd=1, log=T))
      Lpp <- theta.prior + b.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      IRF     <- plogis( thetaLL + bLL )
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qbinom(rep(.5, length(IRF)), size=1, prob=IRF)
      ### Output
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
      return(Modelout)
    }
  } else if (model == 2) {
    ############################## 02 - Double-distributive
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- exp(parm[Data$pos.theta])
      b     <- exp(parm[Data$pos.b])

      ### Log-Priors
      theta.prior <- sum(dlnorm(theta, meanlog=0, sdlog=1, log=T))
      b.prior     <- sum(dlnorm(b    , meanlog=0, sdlog=1, log=T))
      Lpp <- theta.prior + b.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      #IRF     <- plogis( thetaLL * bLL )
      IRF     <- tanh( thetaLL * bLL )
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qbinom(rep(.5, length(IRF)), size=1, prob=IRF)
      ### Output
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
      return(Modelout)
    }
  } else if (model == 3) {
    ############################## 03 - Person-distributive
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- exp(parm[Data$pos.theta])
      b     <- parm[Data$pos.b]

      ### Log-Priors
      theta.prior <- sum(dlnorm(theta, meanlog=0, sdlog=1, log=T))
      b.prior     <- sum(dnorm(b    , mean=0, sd=1, log=T))
      Lpp <- theta.prior + b.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      IRF     <- plogis( log(thetaLL) + {thetaLL * bLL} )
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qbinom(rep(.5, length(IRF)), size=1, prob=IRF)
      ### Output
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
      return(Modelout)
    }
  } else if (model == 4) {
    ############################## 04 - Item-distributive
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- exp(parm[Data$pos.b])

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0   , sd=1, log=T))
      b.prior     <- sum(dlnorm(b   , meanlog=0, sdlog=1, log=T))
      Lpp <- theta.prior + b.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      IRF     <- plogis( log(bLL) + {thetaLL * bLL} )
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qbinom(rep(.5, length(IRF)), size=1, prob=IRF)
      ### Output
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
      return(Modelout)
    }
  } else if (model == 5) {
    ############################## 05 - Interaction
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- parm[Data$pos.b]

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=T))
      b.prior     <- sum(dnorm(b    , mean=0, sd=1, log=T))
      Lpp <- theta.prior + b.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      #IRF     <- plogis( {thetaLL + bLL} + {exp(thetaLL) * exp(bLL)} )
      IRF     <- plogis( {thetaLL + bLL} + {thetaLL * bLL} )
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qbinom(rep(.5, length(IRF)), size=1, prob=IRF)
      ### Output
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
      return(Modelout)
    }
  } else { stop("Unknown model! Check if expecification is correct.") }
  Model <- compiler::cmpfun(Model)
  Initial.Values <- GIV(Model, MyData, PGF=T)
  is.model(Model, Initial.Values, MyData)
  is.bayesian(Model, Initial.Values, MyData)

  ### Run!====
  set.seed(seed)
  if (method=="VB") {
    Iters=Iters; Smpl=Smpl
    Fit <- VariationalBayes(Model=Model, parm=Initial.Values, Data=MyData,
                            Covar=NULL, Interval=1e-6, Iterations=Iters,
                            Method="Salimans2", Samples=Smpl, sir=TRUE,
                            Stop.Tolerance=1e-5, CPUs=CPUs, Type="PSOCK")
  } else if (method=="LA") {
    Iters=Iters; Smpl=Smpl
    Fit <- LaplaceApproximation(Model, parm=Initial.Values, Data=MyData,
                                Interval=1e-6, Iterations=Iters,
                                Method="SPG", Samples=Smpl, sir=TRUE,
                                CovEst="Identity", Stop.Tolerance=1e-5,
                                CPUs=CPUs, Type="PSOCK")
  } else if (method=="MCMC") {
    ## Hit-And-Run Metropolis
    Iters=Iters; Status=Iters/10; Thin=Thin; A=a.s
    Fit <- LaplacesDemon(Model=Model, Data=MyData,
                         Initial.Values=Initial.Values,
                         Covar=NULL, Iterations=Iters,
                         Status=Status, Thinning=Thin,
                         Algorithm="HARM",
                         Specs=list(alpha.star=A, B=NULL))
  } else if (method=="PMC") {
    Iters=Iters; Smpl=Smpl; Thin=Thin
    Fit <- PMC(Model=Model, Data=MyData, Initial.Values=Initial.Values,
               Covar=NULL, Iterations=Iters, Thinning=Thin, alpha=NULL,
               M=2, N=Smpl, nu=1e3, CPUs=CPUs, Type="PSOCK")
  } else if (method=="IQ") {
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
    Fit <- MAP(Model=Model, parm=Initial.Values, Data=MyData, algo=algo, seed=seed,
               maxit=Iters, temp=temp, tmax=tmax, REPORT=Status, Interval=Interval)
  } else {stop('Unknown optimization method.')}

  ### Results====
  if (method=="MAP") {
    abil = Fit$parm[pos.theta]
    diff = Fit$parm[pos.b]
    FI    = Fit$FI

    Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                    'abil'=abil,'diff'=diff,'FitIndexes'=FI)

  } else {
    if (method=="PMC") {
      abil = Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
      diff = Fit$Summary[grep("b", rownames(Fit$Summary), fixed=TRUE),1]
    } else {
      abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
      diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1]
    }
    Dev    <- Fit$Deviance
    mDD    <- Dev - min(Dev)
    pDD    <- Dev[min(which(mDD < 100)):length(Dev)]
    pV     <- var(pDD)/2
    Dbar   <- mean(pDD)
    #Dbar = mean(Dev)
    #pV <- var(Dev)/2
    DIC  = list(DIC=Dbar + pV, Dbar=Dbar, pV=pV)

    Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                    'abil'=abil,'diff'=diff,'DIC'=DIC)

  }
  return(Results)
}
