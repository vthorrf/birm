optscr <- function(x, levels=NULL, M=5, basis="rademacher", err=NULL, knots=NULL,
                   degree=3, method="VB", Iters=500, Smpl=1000, Thin=1, A=500,
                   temp=1e-2, tmax=1, algo="GA", seed=666){

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
  if (is.null(levels)) {
    base_scor <- rowSums(x) / (max(x) * ncol(x))
  } else base_scor <- rowSums(x) / ((levels-1) * ncol(x))
  if (is.null(knots)) knots <- unique(sort(qtrunc(base_scor, "norm", a=-4, b=4)))
  if (is.null(err)) {
    err <- rep(1e-4,nrow(x))
  } else if(length(err) == nrow(x)) {
    err <- err
  } else stop("Error estimates for theta (err) are probably missing or are too many")

  ### Assemble data list====
  if (method == "MAP") {
    mon.names  <- "LL"
  } else { mon.names  <- "LP" }
  if(basis=="legendre") {
    parm.names <- as.parm.names(list( kappa=rep(0,(degree+2) * ncol(x)),
                                      theta=rep(0, nrow(x)) ))
  } else if(basis=="rademacher") {
    parm.names <- as.parm.names(list( kappa=rep(0,length(knots)*ncol(x)),
                                      theta=rep(0, nrow(x))  ))
  } else if(basis=="bs") {
    parm.names <- as.parm.names(list( kappa=rep(0,(length(knots)+
                                                   degree)*ncol(x)),
                                      theta=rep(0, nrow(x)) ))
  } else if(basis=="nn"){
    parm.names <- as.parm.names(list( kappa=rep(0, M * ncol(x) * 2),
                                      theta=rep(0, nrow(x)) ))

  }else stop("Unknow basis type :/")
  K          <- length(grep("kappa", parm.names))
  pos.kappa  <- grep("kappa", parm.names)
  pos.theta  <- grep("theta", parm.names)
  PGF <- function(Data) {
    kappa <- rlaplace(Data$K, 0, 1)
    theta <- rtrunc(Data$n, "norm", -4, 4, mean=Data$SS, sd=Data$err)
    return(c(kappa, theta))
  }
  LLfn   <- function(theta, kappa, knots, degree, n, v, SS, K, M) {
    if(basis=="legendre") {
      ### Legendre Orthogonal Polynomials basis expansion
      kLL   <- matrix(rep(kappa,each=n),ncol=degree+2)
      pol   <- matrix(rep(poly(theta, degree=degree), times=v),
                      ncol=degree, nrow=n * v)
      poll  <- cbind(1,theta,pol)
      W       <- rowSums( kLL * poll )
    } else if(basis=="rademacher") {
      ### Rademacher basis expansion
      kLL   <- matrix(rep(kappa,each=n),ncol=length(knots))
      thetaSS <- rep(theta, times=v)
      thetaLL <- matrix(rep(thetaSS, times=length(knots)),
                        ncol=length(knots))
      RadBas  <- sign(sweep(thetaLL,2,knots)) * kLL
      W       <- rowSums(RadBas)
    } else if(basis=="bs") {
      ### B-spline basis expansion
      kLL   <- matrix(rep(kappa,each=n),ncol=(length(knots)+
                                                     degree))
      thetaSS <- splines::bs(theta + SS, knots=knots,
                             degree=degree)
      BSP  <- matrix(rep(thetaSS, times=v),
                     ncol=degree + length(knots),
                     nrow=n * v, byrow=T)
      W       <- rowSums( kLL *BSP )
    } else if(basis=="nn"){
      ### Neural Network
      inLL  <- matrix(rep(kappa[1:(K/2)],each=n),ncol=M)
      pol   <- matrix(rep(theta, times=v), ncol=M,
                      nrow=n * v)
      otLL  <- matrix(rep(kappa[((K/2)+1):K],each=n),ncol=M)
      W       <- rowSums( plogis(inLL * pol) * otLL )
    }else stop("Unknow basis type :/")
    return(W)
  }
  MyData <- list(parm.names=parm.names, mon.names=mon.names,
                 PGF=PGF, X=data_long, n=nrow(x), v=ncol(x),
                 pos.kappa=pos.kappa, pos.theta=pos.theta, K=K,
                 SS=qtrunc(base_scor, "norm", a=-4, b=4), M=M,
                 knots=knots, degree=degree, err=err)
  is.data(MyData)

  ### Model====
  Model <- function(parm, Data){

    ## Prior parameters
    theta <- interval(parm[Data$pos.theta], -4, 4)
    parm[Data$pos.theta] <- theta
    kappa <- tanh(parm[Data$pos.kappa])

    ### Log-Priors
    theta.prior <- sum(dtrunc(theta, "norm",    -4, 4, mean=Data$SS,
                              sd=Data$err, log=T))
    kappa.prior <- sum(dlaplace(atanh(kappa), 0, 1, log=T))
    Lpp <- kappa.prior + theta.prior

    ### Log-Likelihood
    W       <- LLfn(theta, kappa, Data$knots, Data$degree,
                    Data$n, Data$v, Data$SS, Data$K, Data$M)
    IRF     <- 1 / ( 1 + exp(-W) )
    LL      <- sum( dbinom(Data$X[,3], size=max(Data$X[,3]), prob=IRF, log=T) )

    ### Log-Posterior
    LP <- LL + Lpp
    ### Estimates
    yhat <- qbinom(rep(.5, length(IRF)), size=max(Data$X[,3]), prob=IRF)
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
    Iters=Iters; Status=Iters/10; Thin=Thin; Ad=A
    Fit <- LaplacesDemon(Model=Model, Data=MyData,
                         Initial.Values=Initial.Values,
                         Covar=NULL, Iterations=Iters,Status=Status,
                         Thinning=Thin, Algorithm="NUTS",
                         Specs=list(A=Ad,delta=0.6,epsilon=NULL,Lmax=Inf))
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
    Fit <- MAP(Model=Model, parm=Initial.Values, Data=MyData, algo=algo,
               maxit=Iters, temp=temp, tmax=tmax, REPORT=Status)
  } else {stop('Unknown optimization method.')}

  ### Results====
  if (method=="MAP") {
    ppkp <- Fit$parm[pos.kappa]
    if(basis=="legendre") {
      kappa  = matrix(ppkp,ncol=(degree+2))
      rownames(kappa) <- colnames(x)
      colnames(kappa) <- paste("Kappa_",1:(degree+2),sep="")
    } else if(basis=="rademacher") {
      kappa  = matrix(ppkp,ncol=length(knots))
      rownames(kappa) <- colnames(x)
      colnames(kappa) <- paste("Kappa_",1:length(knots),sep="")
    } else if(basis=="bs") {
      kappa  = matrix(ppkp,ncol=(length(knots)+degree))
      rownames(kappa) <- colnames(x)
      colnames(kappa) <- paste("Kappa_",1:(length(knots)+degree),sep="")
    } else if(basis=="nn"){
      kappain  <- matrix(ppkp[1:(K/2)],ncol=M)
      kappaout <- matrix(ppkp[((K/2)+1):K],ncol=M)
      rownames(kappain) <- rownames(kappaout) <- colnames(x)
      colnames(kappain) <- paste("Kappa_IN_",1:M,sep="")
      colnames(kappaout) <- paste("Kappa_OUT_",1:M,sep="")
      kappa    <- list(kappain, kappaout)
    }else stop("Unknow basis type :/")
    abil <- Fit$parm[pos.theta]
    FI    = Fit$FI

    Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                    'abil'=abil, 'kappa'=kappa,'FitIndexes'=FI)

  } else {
    ppkp <- Fit$Summary1[grep("kappa", rownames(Fit$Summary1), fixed=TRUE),1]
    if(basis=="legendre") {
      kappa  = matrix(ppkp,ncol=(degree+2))
      rownames(kappa) <- colnames(x)
      colnames(kappa) <- paste("Kappa_",1:(degree+2),sep="")
    } else if(basis=="rademacher") {
      kappa  = matrix(ppkp,ncol=length(knots))
      rownames(kappa) <- colnames(x)
      colnames(kappa) <- paste("Kappa_",1:length(knots),sep="")
    } else if(basis=="bs") {
      kappa  = matrix(ppkp,ncol=(length(knots)+degree))
      rownames(kappa) <- colnames(x)
      colnames(kappa) <- paste("Kappa_",1:(length(knots)+degree),sep="")
    } else if(basis=="nn"){
      kappain  <- matrix(ppkp[1:(K/2)],ncol=M)
      kappaout <- matrix(ppkp[((K/2)+1):K],ncol=M)
      rownames(kappain) <- rownames(kappaout) <- colnames(x)
      colnames(kappain) <- paste("Kappa_IN_",1:M,sep="")
      colnames(kappaout) <- paste("Kappa_OUT_",1:M,sep="")
      kappa    <- list(kappain, kappaout)
    }else stop("Unknow basis type :/")
    abil <- Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
    Dev   = Fit$Deviance
    DIC   = list(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)

    Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                    'abil'=abil,'kappa'=kappa,'DIC'=DIC)

  }
  return(Results)
}
