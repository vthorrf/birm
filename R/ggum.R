ggum <- function(x, K=NULL, method="LA", Iters=100, Smpl=1000,
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
  if(is.null(K)) K <- length(table(unlist(c(as.data.frame(x)))))
  if (method == "MAP") {
    mon.names  <- "LL"
  } else { mon.names  <- "LP" }
  parm.names <- as.parm.names(list( theta=rep(0,nrow(x)),
                                    b=rep(0,ncol(x) * {K+1}),
                                    Ds=rep(0,ncol(x)) ))
  pos.theta  <- grep("theta", parm.names)
  pos.b      <- grep("b", parm.names)
  pos.Ds     <- grep("Ds", parm.names)
  PGF <- function(Data) {
    theta <- rnorm(Data$n)
    b     <- rnorm(Data$v * {Data$levels+1})
    Ds    <- rlnorm(Data$v)
    return(c(theta, b, Ds))
  }
  MyData <- list(parm.names=parm.names, mon.names=mon.names,
                 PGF=PGF, X=data_long, n=nrow(x), v=ncol(x), levels=K,
                 pos.theta=pos.theta, pos.b=pos.b, pos.Ds=pos.Ds)
  is.data(MyData)

  ### Model====
  Model <- function(parm, Data){

    ## Prior parameters
    theta <- parm[Data$pos.theta]
    b     <- parm[Data$pos.b]
    Ds    <- interval( parm[Data$pos.Ds], 1e-100, Inf )
    parm[Data$pos.Ds] <- Ds

    ### Log-Priors
    theta.prior <- sum(dnorm(theta,    mean=0,    sd=1, log=T))
    b.prior     <- sum(dnorm(b    ,    mean=0,    sd=1, log=T))
    Ds.prior    <- sum(dlnorm(Ds  , meanlog=0, sdlog=1, log=T))
    Lpp <- theta.prior + b.prior + Ds.prior

    ### Log-Likelihood
    thetaLL  <- rep(theta, times=Data$v)
    bLL      <- matrix(rep(b    , each=Data$n),
                       nrow=nrow(Data$X), ncol=Data$levels+1)
    DLL      <- matrix(rep(Ds    , each=Data$n),
                       nrow=nrow(Data$X), ncol=Data$levels)
    # Summation of tau
    delta    <- bLL[,1]
    tauk     <- matrixStats::rowCumsums(bLL[,-1])
    # Diffs
    diff <- sweep(matrix(thetaLL, nrow=nrow(Data$X), ncol=Data$levels) -
                  matrix(delta, nrow=nrow(Data$X), ncol=Data$levels),
                  2, c(1:Data$levels), "*")
    P1 <- exp(DLL * {diff - tauk})
    # Weight diffs
    K <- matrix(Data$levels, nrow=nrow(Data$X), ncol=Data$levels)
    k <- matrix(rep(c(1:Data$levels), each=nrow(Data$X)),
                nrow=nrow(Data$X), ncol=Data$levels)
    P2  <- exp(DLL * {{{{2*K} - k - 1} * {diff}} - tauk})
    IRF <- {P1 + P2} / matrix(rowSums(P1 + P2), nrow=nrow(Data$X), ncol=Data$levels)
    IRF[which(IRF == 1)] <- 1 - 1e-7
    LL  <- sum( dcat(Data$X[,3], p=IRF, log=T) )

    ### Log-Posterior
    LP <- LL + Lpp
    ### Estimates
    yhat <- tryCatch(qcat(rep(.5, nrow(IRF)), p=IRF),
                     error=function(e) {
                       qbinom(rep(.5, nrow(IRF)), Data$levels-1,
                              rowMeans(IRF)) + min(Data$X[,3])
                     })
    ### Output
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
    return(Modelout)
  }
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
    diff = Fit[["Model"]]$parm[pos.b][c(1:MyData$v)]
    tau = matrix(Fit[["Model"]]$parm[pos.b][-c(1:MyData$v)], nrow=ncol(x))
    rownames(tau) = colnames(x)
    colnames(tau) = paste("Answer_Key",1:K,sep="_")
    disc = Fit[["Model"]]$parm[pos.Ds]
    FI    = Fit$FI

    Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,'abil'=abil,
                    'diff'=diff,"tau"=tau,"disc"=disc,'FitIndexes'=FI)
  } else {
    if (method=="PMC") {
      abil = Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
      diff = Fit$Summary[grep("b", rownames(Fit$Summary), fixed=TRUE),1][c(1:MyData$v)]
      tau  = matrix(Fit$Summary[grep("b", rownames(Fit$Summary),
                                     fixed=TRUE),1][-c(1:MyData$v)],
                    nrow=ncol(x))
      disc = Fit$Summary[grep("Ds", rownames(Fit$Summary), fixed=TRUE),1]
    } else {
      abil = Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
      diff = Fit$Summary1[grep("b", rownames(Fit$Summary1), fixed=TRUE),1][c(1:MyData$v)]
      tau  = matrix(Fit$Summary1[grep("b", rownames(Fit$Summary1),
                                      fixed=TRUE),1][-c(1:MyData$v)],
                    nrow=ncol(x))
      disc = Fit$Summary1[grep("Ds", rownames(Fit$Summary1), fixed=TRUE),1]
    }
    rownames(tau) = colnames(x)
    colnames(tau) = paste("Answer_Key",1:K,sep="_")
    Dev    <- Fit$Deviance
    mDD    <- Dev - min(Dev)
    pDD    <- Dev[min(which(mDD < 100)):length(Dev)]
    pV     <- var(pDD)/2
    Dbar   <- mean(pDD)
    #Dbar = mean(Dev)
    #pV <- var(Dev)/2
    DIC  = list(DIC=Dbar + pV, Dbar=Dbar, pV=pV)

    Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,'abil'=abil,
                    'tau'=tau,'diff'=diff,"disc"=disc,'DIC'=DIC)
  }
  return(Results)
}
