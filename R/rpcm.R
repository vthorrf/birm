rpcm <- function(x, method="LA", ZI=FALSE, p=1, Iters=100, Smpl=1000,
                 Thin=1, a.s=0.234, temp=1e-2, tmax=NULL,
                 algo="GA", seed=666, Interval=1e-8){

  ### Start====
  set.seed(seed)
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
  if (p == 1) {
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)),
                                      sigma=rep(0,ncol(x)) ))
    pos.theta  <- grep("theta", parm.names)
    pos.sigma  <- grep("sigma", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n, mean=0, sd=1)
      sigma <- rnorm(Data$v, mean=0, sd=1)
      return(c(theta, sigma))
    }
    MyData <- list(parm.names=parm.names, mon.names=mon.names,
                   PGF=PGF, X=data_long, n=nrow(x), v=ncol(x),
                   pos.theta=pos.theta, pos.sigma=pos.sigma,
                   yB={data_long$resp != 0} * 1)
    is.data(MyData)

    ### Model====
    if(ZI == T) {
      MyData$X$resp[MyData$X$resp == 0] <- NA
      MyData$yB[MyData$yB == 1] <- NA
      Model <- function(parm, Data){
        # Parameters
        theta <- parm[grep("theta",Data$parm.names)]
        sigma <- parm[grep("sigma",Data$parm.names)]

        # Priors
        TP    <- sum(dnorm(theta, 0, 1, log=T))
        SP    <- sum(dnorm(sigma, 0, 1, log=T))
        Lpp   <- TP + SP

        # Likelihood
        thetaLL <- rep(theta, times=Data$v)
        sigmaLL <- rep(sigma, each= Data$n)
        mu      <- exp(thetaLL + sigmaLL)
        IRF     <- plogis(log(mu))
        IRF[which(IRF == 0)] <- 1e-7
        IRF[which(IRF == 1)] <- 1 - 1e-7
        LLB     <- sum( dbinom(Data$yB, 1, IRF, log=T), na.rm=T )
        LL      <- sum( dpois(Data$X$resp, lambda=mu, log=T), na.rm=T )

        ### Log-Posterior
        LP <- LL + LLB + Lpp
        ### Estimates
        bin  <- qbinom(rep(.5, length(mu)), size=1, prob=IRF)
        yhat <- qpois(rep(.5, length(mu)), mu) * bin
        ### Output
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
        return(Modelout)
      }
    } else if (ZI == F) {
      Model <- function(parm, Data){
        # Parameters
        theta <- parm[grep("theta",Data$parm.names)]
        sigma <- parm[grep("sigma",Data$parm.names)]

        # Priors
        TP    <- sum(dnorm(theta, 0, 1, log=T))
        SP    <- sum(dnorm(sigma, 0, 1, log=T))
        Lpp   <- TP + SP

        # Likelihood
        thetaLL <- rep(theta, times=Data$v)
        sigmaLL <- rep(sigma, each= Data$n)
        mu      <- exp(thetaLL + sigmaLL)
        LL      <- sum( dpois(Data$X$resp, lambda=mu, log=T), na.rm=T )

        ### Log-Posterior
        LP <- LL + Lpp
        ### Estimates
        yhat <- qpois(rep(.5, length(mu)), mu)
        ### Output
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
        return(Modelout)
      }
    } else {stop("Unknown model. ZI should be equal TRUE or FALSE (see documentation).")}

  } else if(p == 2) {
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)),
                                      sigma=rep(0,ncol(x)),
                                      alpha=rep(0,ncol(x)) ))
    pos.theta  <- grep("theta", parm.names)
    pos.sigma  <- grep("sigma", parm.names)
    pos.alpha  <- grep("alpha", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n, mean=0, sd=1)
      sigma <- rnorm(Data$v, mean=0, sd=1)
      alpha <- rep(0, Data$v)
      return(c(theta, sigma, alpha))
    }
    MyData <- list(parm.names=parm.names, mon.names=mon.names,
                   PGF=PGF, X=data_long, n=nrow(x), v=ncol(x),
                   pos.theta=pos.theta, pos.sigma=pos.sigma,
                   pos.alpha=pos.alpha, yB={data_long$resp != 0} * 1)
    is.data(MyData)

    ### Model====
    if (ZI == T) {
      MyData$X$resp[MyData$X$resp == 0] <- NA
      MyData$yB[MyData$yB == 1] <- NA
      Model <- function(parm, Data){
        # Parameters
        theta <- parm[grep("theta",Data$parm.names)]
        sigma <- parm[grep("sigma",Data$parm.names)]
        alpha <- exp(parm[grep("alpha",Data$parm.names)])

        # Priors
        TP    <- sum(dnorm(theta, 0, 1, log=T))
        SP    <- sum(dnorm(sigma, 0, 1, log=T))
        AP    <- sum(dlnorm(alpha, 0, 1, log=T))
        Lpp   <- TP + SP + AP

        # Likelihood
        thetaLL <- rep(theta, times=Data$v)
        sigmaLL <- rep(sigma, each= Data$n)
        alphaLL <- rep(alpha, each= Data$n)
        mu      <- exp(alphaLL * {thetaLL + sigmaLL})
        IRF     <- plogis(log(mu))
        IRF[which(IRF == 0)] <- 1e-7
        IRF[which(IRF == 1)] <- 1 - 1e-7
        LLB     <- sum( dbinom(Data$yB, 1, IRF, log=T), na.rm=T )
        LL      <- sum( dpois(Data$X$resp, lambda=mu, log=T), na.rm=T )

        ### Log-Posterior
        LP <- LL + LLB + Lpp
        ### Estimates
        bin  <- qbinom(rep(.5, length(mu)), size=1, prob=IRF)
        yhat <- qpois(rep(.5, length(mu)), mu) * bin
        ### Output
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
        return(Modelout)
      }
    } else if (ZI == F) {
      Model <- function(parm, Data){
        # Parameters
        theta <- parm[grep("theta",Data$parm.names)]
        sigma <- parm[grep("sigma",Data$parm.names)]
        alpha <- exp(parm[grep("alpha",Data$parm.names)])

        # Priors
        TP    <- sum(dnorm(theta, 0, 1, log=T))
        SP    <- sum(dnorm(sigma, 0, 1, log=T))
        AP    <- sum(dlnorm(alpha, 0, 1, log=T))
        Lpp   <- TP + SP + AP

        # Likelihood
        thetaLL <- rep(theta, times=Data$v)
        sigmaLL <- rep(sigma, each= Data$n)
        alphaLL <- rep(alpha, each= Data$n)
        mu      <- exp(alphaLL * {thetaLL + sigmaLL})
        LL      <- sum( dpois(Data$X$resp, lambda=mu, log=T), na.rm=T )

        ### Log-Posterior
        LP <- LL + Lpp
        ### Estimates
        yhat <- qpois(rep(.5, length(mu)), mu)
        ### Output
        Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
        return(Modelout)
      }
    } else {stop("Unknown model. ZI should be equal TRUE or FALSE (see documentation).")}
  } else {stop("Unknown model. p should be equal to 1 or 2 (see documentation).")}
  Model <- compiler::cmpfun(Model)
  Initial.Values <- GIV(Model, MyData, PGF=T)
  is.model(Model, Initial.Values, MyData)
  is.bayesian(Model, Initial.Values, MyData)

  ### Run!====
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
    abil = Fit[["Model"]]$parm[pos.theta]
    diff = Fit[["Model"]]$parm[pos.sigma]
    if(p == 2) disc = exp(Fit[["Model"]]$parm[pos.alpha])
    FI   = Fit$FI

    if(p == 1) {
      Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                      'abil'=abil,'diff'=diff,'FitIndexes'=FI)
    } else {
      Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                      'abil'=abil,'diff'=diff,'disc'=disc,
                      'FitIndexes'=FI)
    }

  } else {
    if (method=="PMC") {
      abil = Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
      diff = Fit$Summary[grep("sigma", rownames(Fit$Summary), fixed=TRUE),1]
      if(p == 2) disc = exp(Fit$Summary[grep("alpha", rownames(Fit$Summary), fixed=TRUE),1])
    } else {
      abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
      diff = Fit$Summary1[grep("sigma", rownames(Fit$Summary1), fixed=TRUE),1]
      if(p == 2) disc = exp(Fit$Summary1[grep("alpha", rownames(Fit$Summary1), fixed=TRUE),1])
    }
    Dev    <- Fit$Deviance
    mDD    <- Dev - min(Dev)
    pDD    <- Dev[min(which(mDD < 100)):length(Dev)]
    pV     <- var(pDD)/2
    Dbar   <- mean(pDD)
    #Dbar = mean(Dev)
    #pV <- var(Dev)/2
    DIC  = list(DIC=Dbar + pV, Dbar=Dbar, pV=pV)

    if (p == 1) {
      Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                      'abil'=abil,'diff'=diff,'DIC'=DIC)
    } else {
      Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                      'abil'=abil,'diff'=diff,'disc'=disc,
                      'DIC'=DIC)
    }
  }
  return(Results)
}
