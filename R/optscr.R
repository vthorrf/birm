optscr <- function(x, levels=NULL, M=5, basis="rademacher", err=NULL, knots=NULL,
                   degree=3, method="LA", Iters=100, Smpl=1000, Thin=1, a.s=0.234,
                   B=TRUE, temp=1e-2, tmax=NULL, algo="GA", seed=666, Interval=1e-8){

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
  if (is.null(levels)) {
    base_scor <- rowSums(x) / (max(x) * ncol(x))
  } else base_scor <- rowSums(x) / ((levels-1) * ncol(x))
  #if (is.null(knots)) knots <- unique(sort(qtrunc(base_scor, "norm", a=-2, b=2)))
  if (is.null(knots)) knots <- quantile(qtrunc(base_scor, "norm", a=-2, b=2), seq(0, .95, len=5))
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
    parm.names <- as.parm.names(list( kappa=rep(0,(degree) * ncol(x)),
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
    kappa <- rnorm(Data$K, 0, 1)
    #theta <- Data$SS
    #theta <- rtrunc(Data$n, "norm", -2, 2, mean=Data$SS, sd=Data$err)
    theta <- rnorm(Data$n, mean=Data$SS, sd=Data$err)
    return(c(kappa, theta))
  }
  relu   <- function(x) x * (x > 0)
  LLfn   <- function(theta, kappa, knots, degree, n, v, SS, K, M) {
    if(basis=="legendre") {
      ### Legendre Orthogonal Polynomials basis expansion
      kLL   <- matrix(rep(kappa,each=n),ncol=degree)
      #poll <- do.call("rbind", lapply(seq_len(v), function(g) {
      #  poly(theta, degree=degree)
      #}))
      poll  <- apply(poly(theta, degree=degree), 2, rep, times=v)
      W     <- rowSums( kLL * poll )
    } else if(basis=="rademacher") {
      ### Rademacher basis expansion
      kLL   <- matrix(rep(kappa,each=n),ncol=length(knots))
      thetaLL <- matrix(rep(rep(theta, times=v),
                            times=length(knots)),
                        ncol=length(knots))
      RadBas  <- sign(sweep(thetaLL,2,knots,"-")) * kLL
      W       <- rowSums(RadBas)
    } else if(basis=="bs") {
      ### B-spline basis expansion
      kLL   <- matrix(rep(kappa,each=n),
                      ncol=(length(knots)+degree))
      #BSP  <- do.call("rbind", lapply(seq_len(v), function(g) {
      #  splines::bs(theta, knots=knots, degree=degree,
      #              Boundary.knots=c(min(theta), max(theta)))
      #}))
      #BSP  <- apply(splines::bs(theta, knots=knots, degree=degree,
      #              Boundary.knots=c(min(theta), max(theta))), 2, rep, times=v)
      BSP  <- apply(splines::bs(theta, knots=knots, degree=degree,
                    Boundary.knots=c(min(knots), max(knots))), 2, rep, times=v)
      W       <- rowSums( kLL * BSP )
    } else if(basis=="nn"){
      ### Neural Network
      inLL  <- matrix(rep(kappa[1:(K/2)],each=n),ncol=M)
      pol   <- matrix(rep(theta, times=v), ncol=M,
                      nrow=n * v)
      otLL  <- matrix(rep(kappa[((K/2)+1):K],each=n),ncol=M)
      W     <- rowSums( relu(inLL * pol) * otLL )
    }else stop("Unknow basis type :/")
    return(W)
  }
  MyData <- list( parm.names=parm.names, mon.names=mon.names,
                  PGF=PGF, X=data_long, n=nrow(x), v=ncol(x),
                  pos.kappa=pos.kappa, pos.theta=pos.theta, K=K,
                  SS=((base_scor * 4) -2), M=M, knots=knots,
                  degree=degree, err=err, ssecdf=ecdf((base_scor * 4) -2) )
  is.data(MyData)

  ### Model====
  Model <- function(parm, Data){

    ## Prior parameters
    #theta <- interval(parm[Data$pos.theta], -2, 2)
    theta <- tanh(parm[Data$pos.theta]) * 2
    #parm[Data$pos.theta] <- theta
    kappa <- tanh(parm[Data$pos.kappa])

    ### Log-Priors
    tecdf       <- ecdf(theta)
    ss.tt.d     <- Data$ssecdf(Data$SS) - tecdf(theta)
    theta.prior <- sum(dlaplace(ss.tt.d, location=0, scale=1e-4, log=T))
    kappa.prior <- sum(dnorm(atanh(kappa), 0, 1, log=T))
    Lpp <- kappa.prior + theta.prior

    ### Log-Likelihood
    W       <- LLfn(theta, kappa, Data$knots, Data$degree,
                    Data$n, Data$v, Data$SS, Data$K, Data$M)
    IRF     <- 1 / ( 1 + exp(-W) )
    IRF[which(IRF == 1)] <- 1 - 1e-7
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
    if (B == T) {
      B=list(grep("theta",parm.names),
             seq_len(length(parm.names))[-grep("theta",parm.names)])
    } else if (B == F) {
      B = NULL
    } else { stop(paste("Don't know what you mean with B = ", B, sep="")) }
    Fit <- LaplacesDemon(Model=Model, Data=MyData,
                         Initial.Values=Initial.Values,
                         Covar=NULL, Iterations=Iters,
                         Status=Status, Thinning=Thin,
                         Algorithm="HARM",
                         Specs=list(alpha.star=A, B=B))
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
    ppkp <- Fit[["Model"]]$parm[pos.kappa]
    if(basis=="legendre") {
      kappa  = matrix(ppkp,ncol=(degree))
      rownames(kappa) <- colnames(x)
      colnames(kappa) <- paste("Kappa_",1:(degree),sep="")
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
    abil <- Fit[["Model"]]$parm[pos.theta]
    FI    = Fit$FI

    Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                    'abil'=abil, 'kappa'=kappa,'FitIndexes'=FI)

  } else {
    if (method=="PMC") {
      ppkp <- Fit$Summary[grep("kappa", rownames(Fit$Summary), fixed=TRUE),1]
      abil <- Fit$Summary[grep("theta", rownames(Fit$Summary), fixed=TRUE),1]
    } else {
      ppkp <- Fit$Summary1[grep("kappa", rownames(Fit$Summary1), fixed=TRUE),1]
      abil <- Fit$Summary1[grep("theta", rownames(Fit$Summary1), fixed=TRUE),1]
    }
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
    } else stop("Unknow basis type :/")
    Dev    <- Fit$Deviance
    mDD    <- Dev - min(Dev)
    pDD    <- Dev[min(which(mDD < 100)):length(Dev)]
    pV     <- var(pDD)/2
    Dbar   <- mean(pDD)
    #Dbar = mean(Dev)
    #pV <- var(Dev)/2
    DIC  = list(DIC=Dbar + pV, Dbar=Dbar, pV=pV)

    Results <- list("Data"=MyData,"Fit"=Fit,"Model"=Model,
                    'abil'=abil,'kappa'=kappa,'DIC'=DIC)
  }
  return(Results)
}
