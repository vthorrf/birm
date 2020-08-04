pcm <- function(x, levels=NULL, p=1, method="LA", Iters=100, Smpl=1000, Thin=1,
                a.s=0.234, temp=1e-2, tmax=1, algo="GA", seed=666){

  ### Start====
  #require(LaplacesDemon)
  #require(compiler)
  #require(parallel)
  #require(tidyr)
  #require(matrixStats)
  CPUs = detectCores(all.tests = FALSE, logical = TRUE) - 1
  if(CPUs == 0) CPUs = 1

  ### Convert data to long format====
  if (is.null(levels)) levels <- max(x)
  lonlong <- gather(data.frame(x), "item", "resp", colnames(x), factor_key=TRUE)
  data_long <- data.frame(ID=rep(1:nrow(x), times=ncol(x)),lonlong)

  ### Choose model====
  if (p == 1) {
    ## Partial Credit Model====
    # Assemble data list
    if (method == "MAP") {
      mon.names  <- "LL"
    } else { mon.names  <- "LP" }
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)), b=rep(0,ncol(x) * levels) ))
    pos.theta  <- grep("theta", parm.names)
    pos.b      <- grep("b", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n)
      b     <- rnorm(Data$v * Data$levels)
      return(c(theta, b))
    }
    MyData <- list(parm.names=parm.names, mon.names=mon.names,
                   PGF=PGF, X=data_long, n=nrow(x), v=ncol(x), levels=levels,
                   pos.theta=pos.theta, pos.b=pos.b)
    is.data(MyData)

    # Model
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- parm[Data$pos.b]

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=T))
      b.prior     <- sum(dnorm(b    , mean=0, sd=1, log=T))
      Lpp <- theta.prior + b.prior

      ### Log-Likelihood
      thetaLL  <- rep(theta, times=Data$v)
      BsL      <- rep(b    , each=Data$n)
      bLL      <- matrix(BsL , nrow=nrow(Data$X), ncol=Data$levels)
      eta      <- thetaLL - bLL
      exp.psum <- exp(matrixStats::rowCumsums(eta))
      IRF      <- exp.psum / rowSums(exp.psum)
      LL       <- sum( dcat(Data$X[,3], p=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qcat(rep(.5, nrow(IRF)), p=IRF)
      ### Output
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
      return(Modelout)
    }
    Model <- compiler::cmpfun(Model)
    Initial.Values <- GIV(Model, MyData, PGF=T)
    is.model(Model, Initial.Values, MyData)
    is.bayesian(Model, Initial.Values, MyData)

  } else if (p == 2) {
    ## Generalized Partial Credit Model====
    # Assemble data list
    mon.names  <- "LP"
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)), b=rep(0,ncol(x) * levels),
                                      Ds=rep(0,ncol(x)) ))
    pos.theta  <- grep("theta", parm.names)
    pos.b      <- grep("b", parm.names)
    pos.Ds    <- grep("Ds", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n)
      b     <- rnorm(Data$v * Data$levels)
      Ds    <- rnorm(Data$v)
      return(c(theta, b, Ds))
    }
    MyData <- list(parm.names=parm.names, mon.names=mon.names,
                   PGF=PGF, X=data_long, n=nrow(x), v=ncol(x), levels=levels,
                   pos.theta=pos.theta, pos.b=pos.b, pos.Ds=pos.Ds)
    is.data(MyData)

    # Model
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- parm[Data$pos.b]
      Ds    <- parm[Data$pos.Ds]

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=T))
      b.prior     <- sum(dnorm(b    , mean=0, sd=1, log=T))
      Ds.prior    <- sum(dnorm(Ds   , mean=0, sd=1, log=T))
      Lpp <- theta.prior + b.prior + Ds.prior

      ### Log-Likelihood
      thetaLL  <- rep(theta, times=Data$v)
      BsL      <- rep(b    , each=Data$n)
      bLL      <- matrix(BsL , nrow=nrow(Data$X), ncol=Data$levels)
      DLL      <- rep(Ds   , each=Data$n)
      eta      <- DLL * (thetaLL - bLL)
      exp.psum <- exp(matrixStats::rowCumsums(eta))
      IRF      <- exp.psum / rowSums(exp.psum)
      LL       <- sum( dcat(Data$X[,3], p=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qcat(rep(.5, nrow(IRF)), p=IRF)
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
  } else if (method=="MAP") {
    ## Maximum a Posteriori====
    #Iters=100; Smpl=1000
    Iters=Iters; Status=Iters/10
    Fit <- MAP(Model=Model, parm=Initial.Values, Data=MyData, algo=algo,
               maxit=Iters, temp=temp, tmax=tmax, REPORT=Status)
  } else {stop('Unknown optimization method.')}

  ### Results====
  if (p == 1) {

    if (method=="MAP") {
      abil = Fit$parm[pos.theta]
      diff = matrix(Fit$parm[pos.b], nrow=ncol(x))
      rownames(diff) = colnames(x)
      colnames(diff) = paste("Answer_Key",1:levels,sep="_")
      FI    = Fit$FI

      Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                      'abil'=abil,'diff'=diff,'FitIndexes'=FI)
    } else {
      if (method=="PMC") {
        abil = Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
        diff = matrix(Fit$Summary[grep("b", rownames(Fit$Summary), fixed=TRUE),1],
                      nrow=ncol(x))
      } else {
        abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
        diff = matrix(Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1],
                      nrow=ncol(x))
      }
      rownames(diff) = colnames(x)
      colnames(diff) = paste("Answer_Key",1:levels,sep="_")
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

  } else if (p == 2) {

    if (method=="MAP") {
      abil = Fit$parm[pos.theta]
      diff = matrix(Fit$parm[pos.b], nrow=ncol(x))
      rownames(diff) = colnames(x)
      colnames(diff) = paste("Answer_Key",1:levels,sep="_")
      disc = Fit$parm[pos.Ds]
      FI    = Fit$FI

      Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                      'abil'=abil,'diff'=diff,"disc"=disc,'FitIndexes'=FI)
    } else {
      if (method=="PMC") {
        abil = Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
        diff = matrix(Fit$Summary[grep("b", rownames(Fit$Summary), fixed=TRUE),1],
                      nrow=ncol(x))
        disc = Fit$Summary[grep("Ds", rownames(Fit$Summary), fixed=TRUE),1]
      } else {
        abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
        diff = matrix(Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1],
                      nrow=ncol(x))
        disc = Fit$Summary1[grep("Ds", rownames(Fit$Summary1), fixed=TRUE),1]
      }
      rownames(diff) = colnames(x)
      colnames(diff) = paste("Answer_Key",1:levels,sep="_")
      Dev    <- Fit$Deviance
      mDD    <- Dev - min(Dev)
      pDD    <- Dev[min(which(mDD < 100)):length(Dev)]
      pV     <- var(pDD)/2
      Dbar   <- mean(pDD)
      #Dbar = mean(Dev)
      #pV <- var(Dev)/2
      DIC  = list(DIC=Dbar + pV, Dbar=Dbar, pV=pV)

      Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                      'abil'=abil,'diff'=diff,"disc"=disc,'DIC'=DIC)
  }

  } else warning("Can't return any result :P")

  return(Results)
}
