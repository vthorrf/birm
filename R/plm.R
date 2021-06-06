plm <- function(x, p=1, scaling=1.7, method="LA",
                Iters=100, Smpl=1000, Thin=1,
                a.s=0.234, temp=1e-2, tmax=NULL, algo="GA",
                seed=666, Interval=1e-8){

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

  ### Choose Model====
  if (p == 1) {
    ## One-parameter logistic====
    # Assemble data list
    if (method == "MAP") {
      mon.names  <- "LL"
    } else { mon.names  <- "LP" }
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)), b=rep(0,ncol(x)),
                                      Ds=rep(0,1) ))
    pos.theta  <- grep("theta", parm.names)
    pos.b      <- grep("b", parm.names)
    pos.Ds    <- grep("Ds", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n)
      b     <- rnorm(Data$v)
      Ds    <- rnorm(1)
      #Ds    <- rlnorm(1)
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
      #Ds    <- interval( parm[Data$pos.Ds], 1e-100, Inf )
      #parm[Data$pos.Ds] <- Ds
      Ds    <- exp( parm[Data$pos.Ds] )

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=TRUE))
      b.prior     <- sum(dnorm(b, mean=0, sd=1, log=T))
      Ds.prior    <- sum(dlnorm(Ds, 0, 1, log=T))
      Lpp <- theta.prior + b.prior + Ds.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      DLL     <- rep(Ds   , each=Data$v * Data$n)
      IRF     <- plogis( Data$scaling * DLL * (thetaLL - bLL) )
      IRF[which(IRF == 1)] <- 1 - 1e-7
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qbinom(rep(.5, length(IRF)), size=1, prob=IRF)
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
      #Ds    <- rlnorm(Data$v)
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
      #Ds    <- interval( parm[Data$pos.Ds], 1e-100, Inf )
      #parm[Data$pos.Ds] <- Ds
      Ds    <- exp( parm[Data$pos.Ds] )

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=TRUE))
      b.prior     <- sum(dnorm(b, mean=0, sd=1, log=T))
      Ds.prior    <- sum(dlnorm(Ds, 0, 1, log=T))
      Lpp <- theta.prior + b.prior + Ds.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      DLL     <- rep(Ds   , each=Data$n)
      IRF     <- plogis( DLL * (thetaLL - bLL) )
      IRF[which(IRF == 1)] <- 1 - 1e-7
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qbinom(rep(.5, length(IRF)), size=1, prob=IRF)
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
      c     <- rnorm(Data$v)
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
      #Ds    <- interval( parm[Data$pos.Ds], 1e-100, Inf )
      Ds    <- exp( parm[Data$pos.Ds] )
      #parm[Data$pos.Ds] <- Ds
      #c     <- interval( parm[Data$pos.c], 1e-100, (1 - 1e-100) )
      c     <- plogis( parm[Data$pos.c] ) / 2
      #parm[Data$pos.c] <- c

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=TRUE))
      b.prior     <- sum(dnorm(b, mean=0, sd=1, log=T))
      Ds.prior    <- sum(dlnorm(Ds, 0, 1, log=T))
      c.prior     <- sum(dbeta(c, 1, 1, log=T))
      Lpp <- theta.prior + b.prior + Ds.prior + c.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      DLL     <- rep(Ds   , each=Data$n)
      cLL     <- rep(c    , each=Data$n)
      IRF     <- cLL + ( (1 - cLL) / (1 + exp(-DLL * ( thetaLL - bLL ))) )
      IRF[which(IRF == 1)] <- 1 - 1e-7
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qbinom(rep(.5, length(IRF)), size=1, prob=IRF)
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
      #Ds    <- rlnorm(Data$v)
      c     <- rnorm(Data$v)
      #c     <- rbeta(Data$v, 1, 1)
      UA    <- rnorm(Data$v)
      #UA    <- rbeta(Data$v, 1, 1)
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
      Ds    <- exp(parm[Data$pos.Ds])
      #Ds    <- interval( parm[Data$pos.Ds], 1e-100, Inf )
      #parm[Data$pos.Ds] <- Ds
      c     <- plogis(parm[Data$pos.c])/2
      #c     <- interval( parm[Data$pos.c], 1e-100, (1 - 1e-100) )
      #parm[Data$pos.c] <- c
      UA    <- {plogis(parm[Data$pos.UA]) + 1} / 2
      #UA    <- interval( parm[Data$pos.UA], 1e-100, (1 - 1e-100) )
      #parm[Data$pos.UA] <- UA

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=TRUE))
      b.prior     <- sum(dnorm(b, mean=0, sd=1, log=T))
      Ds.prior    <- sum(dlnorm(Ds, 0, 1, log=T))
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
      IRF[which(IRF == 1)] <- 1 - 1e-7
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qbinom(rep(.5, length(IRF)), size=1, prob=IRF)
      ### Output
      Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
      return(Modelout)
    }
    Model <- compiler::cmpfun(Model)
    Initial.Values <- GIV(Model, MyData, PGF=T)
    is.model(Model, Initial.Values, MyData)
    is.bayesian(Model, Initial.Values, MyData)

  } else  if (p == 5) {
    ## Five-parameter logistic====
    # Assemble data list
    mon.names  <- "LP"
    parm.names <- as.parm.names(list( theta=rep(0,nrow(x)), b=rep(0,ncol(x)),
                                      Ds=rep(0,ncol(x)), c=rep(0,ncol(x)),
                                      UA=rep(0,ncol(x)), AS=rep(0,ncol(x)) ))
    pos.theta  <- grep("theta", parm.names)
    pos.b      <- grep("b", parm.names)
    pos.Ds     <- grep("Ds", parm.names)
    pos.c      <- grep("c", parm.names)
    pos.UA     <- grep("UA", parm.names)
    pos.AS     <- grep("AS", parm.names)
    PGF <- function(Data) {
      theta <- rnorm(Data$n)
      b     <- rnorm(Data$v)
      Ds    <- rnorm(Data$v)
      #Ds    <- rlnorm(Data$v)
      c     <- rnorm(Data$v)
      #c     <- rbeta(Data$v, 1, 1)
      UA    <- rnorm(Data$v)
      #UA    <- rbeta(Data$v, 1, 1)
      AS    <- rnorm(Data$v)
      #AS    <- rlnorm(Data$v)
      return(c(theta, b, Ds, c, UA, AS))
    }
    MyData <- list(parm.names=parm.names, mon.names=mon.names,
                   PGF=PGF, X=data_long, n=nrow(x), v=ncol(x),
                   pos.theta=pos.theta, pos.b=pos.b, pos.c=pos.c,
                   pos.UA=pos.UA, pos.Ds=pos.Ds, pos.ASs=pos.AS)
    is.data(MyData)

    # Model
    Model <- function(parm, Data){

      ## Prior parameters
      theta <- parm[Data$pos.theta]
      b     <- parm[Data$pos.b]
      Ds    <- exp(parm[Data$pos.Ds])
      #Ds    <- interval( parm[Data$pos.Ds], 1e-100, Inf )
      #parm[Data$pos.Ds] <- Ds
      c     <- plogis(parm[Data$pos.c])/2
      #c     <- interval( parm[Data$pos.c], 1e-100, (1 - 1e-100) )
      #parm[Data$pos.c] <- c
      UA    <- {plogis(parm[Data$pos.UA]) + 1} / 2
      #UA    <- interval( parm[Data$pos.UA], 1e-100, (1 - 1e-100) )
      #parm[Data$pos.UA] <- UA
      AS    <- exp(parm[Data$pos.AS])
      #AS    <- interval( parm[Data$pos.AS], 1e-100, Inf )
      #parm[Data$pos.AS] <- AS

      ### Log-Priors
      theta.prior <- sum(dnorm(theta, mean=0, sd=1, log=TRUE))
      b.prior     <- sum(dnorm(b, mean=0, sd=1, log=T))
      Ds.prior    <- sum(dlnorm(Ds, 0, 1, log=T))
      c.prior     <- sum(dbeta(c, 1, 1, log=T))
      UA.prior    <- sum(dbeta(UA, 1, 1, log=T))
      AS.prior    <- sum(dlnorm(AS, 0, 1, log=T))
      Lpp <- theta.prior + b.prior + Ds.prior + c.prior + UA.prior + AS.prior

      ### Log-Likelihood
      thetaLL <- rep(theta, times=Data$v)
      bLL     <- rep(b    , each=Data$n)
      DLL     <- rep(Ds   , each=Data$n)
      cLL     <- rep(c    , each=Data$n)
      UALL    <- rep(UA   , each=Data$n)
      SLL     <- rep(AS   , each=Data$n)
      IRF     <- cLL + ( (UALL - cLL) / ((1 + exp(-DLL * ( thetaLL - bLL )))^SLL) )
      IRF[which(IRF == 1)] <- 1 - 1e-7
      LL      <- sum( dbinom(Data$X[,3], size=1, prob=IRF, log=T) )

      ### Log-Posterior
      LP <- LL + Lpp
      ### Estimates
      yhat <- qbinom(rep(.5, length(IRF)), size=1, prob=IRF)
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
    Fit <- MAP(Model=Model, parm=Initial.Values, Data=MyData, algo=algo, seed=seed,
               maxit=Iters, temp=temp, tmax=tmax, REPORT=Status, Interval=Interval)
  } else {stop('Unknown estimation method.')}

  ### Results====
  if (p == 1) {
    if (method=="MAP") {
      abil = Fit[["Model"]]$parm[pos.theta]
      diff = Fit[["Model"]]$parm[pos.b]
      disc = exp( Fit[["Model"]]$parm[pos.Ds] )
      FI    = Fit$FI

      Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                      'abil'=abil,'diff'=diff,"disc"=disc,'FitIndexes'=FI)
    } else {
    ## One-parameter logistic====
    if (method=="PMC") {
      abil = Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
      diff = Fit$Summary[grep("b", rownames(Fit$Summary), fixed=TRUE),1]
      disc = exp( Fit$Summary[grep("Ds", rownames(Fit$Summary), fixed=TRUE),1] )
    } else {
      abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
      diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1]
      disc = exp( Fit$Summary1[grep("Ds", rownames(Fit$Summary1), fixed=TRUE),1] )
    }
    Dev    <- Fit$Deviance
    mDD    <- Dev - min(Dev)
    pDD    <- Dev[min(which(mDD < 100)):length(Dev)]
    pV     <- var(pDD)/2
    Dbar   <- mean(pDD)
    #Dbar = mean(Dev)
    #pV <- var(Dev)/2
    DIC  = list(DIC=Dbar + pV, Dbar=Dbar, pV=pV)

    Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                    'abil'=abil,'diff'=diff,"disc"=disc,'DIC'=DIC)
    }
  } else if (p == 2) {
    if (method=="MAP") {
      abil = Fit[["Model"]]$parm[pos.theta]
      diff = Fit[["Model"]]$parm[pos.b]
      disc = exp( Fit[["Model"]]$parm[pos.Ds] )
      FI    = Fit$FI

      Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                      'abil'=abil,'diff'=diff,"disc"=disc,'FitIndexes'=FI)
    } else {
    ## Two-parameter logistic====
    if (method=="PMC") {
      abil = Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
      diff = Fit$Summary[grep("b", rownames(Fit$Summary), fixed=TRUE),1]
      disc = exp( Fit$Summary[grep("Ds", rownames(Fit$Summary), fixed=TRUE),1] )
    } else {
      abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
      diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1]
      disc = exp( Fit$Summary1[grep("Ds", rownames(Fit$Summary1), fixed=TRUE),1] )
    }
    Dev    <- Fit$Deviance
    mDD    <- Dev - min(Dev)
    pDD    <- Dev[min(which(mDD < 100)):length(Dev)]
    pV     <- var(pDD)/2
    Dbar   <- mean(pDD)
    #Dbar = mean(Dev)
    #pV <- var(Dev)/2
    DIC  = list(DIC=Dbar + pV, Dbar=Dbar, pV=pV)

    Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                    'abil'=abil,'diff'=diff,"disc"=disc,'DIC'=DIC)
    }
  } else if (p == 3) {
    if (method=="MAP") {
      abil = Fit[["Model"]]$parm[pos.theta]
      diff = Fit[["Model"]]$parm[pos.b]
      disc = exp(Fit[["Model"]]$parm[pos.Ds])
      gues = plogis(Fit[["Model"]]$parm[pos.c])/2
      FI    = Fit$FI

      Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                      'abil'=abil,'diff'=diff,"disc"=disc,"gues"=gues,'FitIndexes'=FI)
    } else {
    ## Three-parameter logistic====
    if (method=="PMC") {
      abil = Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
      diff = Fit$Summary[grep("b", rownames(Fit$Summary), fixed=TRUE),1]
      disc = exp(Fit$Summary[grep("Ds", rownames(Fit$Summary), fixed=TRUE),1])
      gues = plogis(Fit$Summary[grep("c", rownames(Fit$Summary), fixed=TRUE),1])/2
    } else {
      abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
      diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1]
      disc = exp(Fit$Summary1[grep("Ds", rownames(Fit$Summary1), fixed=TRUE),1])
      gues = plogis(Fit$Summary1[grep("c", rownames(Fit$Summary1), fixed=TRUE),1])/2
    }
    Dev    <- Fit$Deviance
    mDD    <- Dev - min(Dev)
    pDD    <- Dev[min(which(mDD < 100)):length(Dev)]
    pV     <- var(pDD)/2
    Dbar   <- mean(pDD)
    #Dbar = mean(Dev)
    #pV <- var(Dev)/2
    DIC  = list(DIC=Dbar + pV, Dbar=Dbar, pV=pV)

    Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                    'abil'=abil,'diff'=diff,"disc"=disc,"gues"=gues,'DIC'=DIC)
    }
  } else if (p == 4) {
    if (method=="MAP") {
      abil = Fit[["Model"]]$parm[pos.theta]
      diff = Fit[["Model"]]$parm[pos.b]
      disc = exp(Fit[["Model"]]$parm[pos.Ds])
      gues = plogis(Fit[["Model"]]$parm[pos.c])/2
      UpAs = {plogis(Fit[["Model"]]$parm[pos.UA])+1}/2
      FI    = Fit$FI

      Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                      'abil'=abil,'diff'=diff,"disc"=disc,"gues"=gues,
                      "UpAs"=UpAs,'FitIndexes'=FI)
    } else {
    ## Four-parameter logistic====
    if (method=="PMC") {
      abil = Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
      diff = Fit$Summary[grep("b", rownames(Fit$Summary), fixed=TRUE),1]
      disc = exp(Fit$Summary[grep("Ds", rownames(Fit$Summary), fixed=TRUE),1])
      gues = plogis(Fit$Summary[grep("c", rownames(Fit$Summary), fixed=TRUE),1])/2
      UpAs = {plogis(Fit$Summary[grep("UA", rownames(Fit$Summary), fixed=TRUE),1])+1}/2
    } else {
      abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
      diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1]
      disc = exp(Fit$Summary1[grep("Ds", rownames(Fit$Summary1), fixed=TRUE),1])
      gues = plogis(Fit$Summary1[grep("c", rownames(Fit$Summary1), fixed=TRUE),1])/2
      UpAs = {plogis(Fit$Summary1[grep("UA", rownames(Fit$Summary1), fixed=TRUE),1])+1}/2
    }
    Dev    <- Fit$Deviance
    mDD    <- Dev - min(Dev)
    pDD    <- Dev[min(which(mDD < 100)):length(Dev)]
    pV     <- var(pDD)/2
    Dbar   <- mean(pDD)
    #Dbar = mean(Dev)
    #pV <- var(Dev)/2
    DIC  = list(DIC=Dbar + pV, Dbar=Dbar, pV=pV)

    Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                    'abil'=abil,'diff'=diff,"disc"=disc,"gues"=gues,
                    "UpAs"=UpAs,'DIC'=DIC)
    }
  } else if (p == 5) {
    if (method=="MAP") {
      abil = Fit[["Model"]]$parm[pos.theta]
      diff = Fit[["Model"]]$parm[pos.b]
      disc = exp(Fit[["Model"]]$parm[pos.Ds])
      gues =  plogis(Fit[["Model"]]$parm[pos.c])/2
      UpAs = {plogis(Fit[["Model"]]$parm[pos.UA])+1}/2
      asym = exp(Fit[["Model"]]$parm[pos.AS])
      FI   = Fit$FI

      Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                      'abil'=abil,'diff'=diff,"disc"=disc,"gues"=gues,
                      "UpAs"=UpAs,"asym"=asym,'FitIndexes'=FI)
    } else {
      ## Five-parameter logistic====
      if (method=="PMC") {
        abil = Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
        diff = Fit$Summary[grep("b", rownames(Fit$Summary), fixed=TRUE),1]
        disc = exp(Fit$Summary[grep("Ds", rownames(Fit$Summary), fixed=TRUE),1])
        gues = plogis(Fit$Summary[grep("c", rownames(Fit$Summary), fixed=TRUE),1])/2
        UpAs = {plogis(Fit$Summary[grep("UA", rownames(Fit$Summary), fixed=TRUE),1])+1}/2
        asym = exp(Fit$Summary[grep("AS", rownames(Fit$Summary), fixed=TRUE),1])
      } else {
        abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
        diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1]
        disc = exp(Fit$Summary1[grep("Ds", rownames(Fit$Summary1), fixed=TRUE),1])
        gues = plogis(Fit$Summary1[grep("c", rownames(Fit$Summary1), fixed=TRUE),1])/2
        UpAs = {plogis(Fit$Summary1[grep("UA", rownames(Fit$Summary1), fixed=TRUE),1])+1}/2
        asym = exp(Fit$Summary1[grep("AS", rownames(Fit$Summary1), fixed=TRUE),1])
      }
      Dev    <- Fit$Deviance
      mDD    <- Dev - min(Dev)
      pDD    <- Dev[min(which(mDD < 100)):length(Dev)]
      pV     <- var(pDD)/2
      Dbar   <- mean(pDD)
      #Dbar = mean(Dev)
      #pV <- var(Dev)/2
      DIC  = list(DIC=Dbar + pV, Dbar=Dbar, pV=pV)

      Results <- list("Data"=MyData,"Model"=Model,"Fit"=Fit,
                      'abil'=abil,'diff'=diff,"disc"=disc,"gues"=gues,
                      "UpAs"=UpAs,"asym"=asym,'DIC'=DIC)
    }
  } else warning("Can't return any result :P")

  return(Results)
}
