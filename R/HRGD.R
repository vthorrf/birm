### Numerical gradient====
sHAR <- function(par, mod, ...) {
  theta <- rnorm(length(par))
  d     <- theta / sqrt(sum(theta*theta))
  u     <- optim(.5, steepA, Model=mod, df=d, ..., est=par,
                 method="Brent", upper=1, lower=1e-10)$par
  prop  <- par + {u * d}
  return(list("prop"=prop, "u"=u, "d"=d))
}

### Hit-and-Run Gradient Descent MAP====
HRGD <- function(Model, startvalue, Data, Interval=1e-8, maxit=100) {
  # Opening message
  cat("Hit-and-Run Gradient Descent Maximum a Posteriori estimation will run for ",
      maxit, " iterations.\n\n", sep="")
  startTime = proc.time()

  # Initial settings
  pb   <- txtProgressBar(min=0, max=maxit, style=3)
  par  <- df <- array(dim = c(maxit,length(startvalue)))
  f    <- vector("numeric", length=maxit)
  step <- vector("numeric", length=maxit)
  convergence = T

  # First step
  df[1,]   <- grad(Model, startvalue, Data, Interval=Interval)
  SS       <- suppressWarnings(unlist(optim(c(0,0), steep, Model=Model,
                                            Data=Data, df=df[1,],
                                            est=startvalue)$par))
  if (sum(abs(SS)) == 0) {
    SS       <- suppressWarnings(unlist(ucminf(c(0,0), steep, Model=Model,
                                               Data=Data, df=df[1,],
                                               est=startvalue)$par))
  }
  step[1]   <- mean(SS)
  par[1,]   <- startvalue +
               {SS[1]*df[1,]*{(abs(df[1,]) >  sd(df[1,])) * 1}} +
               {SS[2]*df[1,]*{(abs(df[1,]) <= sd(df[1,])) * 1}}
  f[1]      <- Model(par[1,], Data)[["LP"]]
  setTxtProgressBar(pb, 1)

  # Start estimation
  for(run in 2:maxit) {
    # Calculate gradient and estimate step parameter
    temp <- sHAR(par=par[run-1,], mod=Model, Data=Data)
    df[run,] <- temp$d
    step[run] <- temp$u
    par[run,] <- temp$prop
    f[run]    <- Model(par[run,], Data)[["LP"]]
    setTxtProgressBar(pb, run)
  }
  cat("\nConvergence achieved!")
  #if (run < maxit) {
  #  f    <- f[-c({run+1}:maxit)]
  #  par  <- par[-c({run+1}:maxit),]
  #  df   <- df[-c({run+1}:maxit),]
  #  step <- step[-c({run+1}:maxit),]
  #}
  close(pb)
  cat("\n")

  # Final messages
  stopTime = proc.time()
  elapsedTime = stopTime - startTime
  cat("It took ",round(elapsedTime[3],2)," secs for the run to finish. \n", sep="")

  # Return results
  Results <- list("LogPosterior"=f, "Estimates"=par, "Gradients"=df,
                  "MAP"=par[which.max(f),], "step"=step)
  return(Results)
}
