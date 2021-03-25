### Numerical gradient====
post <- function(Model, Data, par) {
  return(Model(par, Data)[["LP"]])
}
grad <- function(Model, par, Data, Interval=1e-8) {
  mat  <- matrix(par, nrow=length(par), ncol=length(par))
  diag(mat) <- par + Interval
  df <- {apply(mat, 2, post, Model=Model, Data=Data) - Model(par, Data)[["LP"]]} / Interval
  df[which(!is.finite(df))] <- 0
  return(df)
}
steep <- function(Model, par, Data, est, df) {
  -Model(est + {par[1]*df*{(abs(df) >  sd(df)) * 1}} +
               {par[2]*df*{(abs(df) <= sd(df)) * 1}}, Data)[["LP"]]
}

### Steepest N-group Gradient Descent MAP====
SNGD <- function(Model, startvalue, Data, Interval=1e-8, maxit=100) {
  # Opening message
  cat("Steepest N-group Gradient Descent Maximum a Posteriori estimation will run for ",
      maxit, " iterations.\n\n", sep="")
  startTime = proc.time()

  # Initial settings
  pb   <- txtProgressBar(min=0, max=maxit, style=3)
  par  <- df <- array(dim = c(maxit,length(startvalue)))
  f    <- vector("numeric", length=maxit)
  step <- array(dim = c(maxit,2))

  # First step
  df[1,]   <- grad(Model, startvalue, Data, Interval=Interval)
  step[1,] <- suppressWarnings(unlist(optim(c(0,0), steep, Model=Model,
                                            Data=Data, df=df[1,], est=startvalue,
                                            method="Nelder-Mead")$par))
  par[1,]   <- startvalue +
    {step[1,1]*df[1,]*{(abs(df[1,]) >  sd(df[1,])) * 1}} +
    {step[1,2]*df[1,]*{(abs(df[1,]) <= sd(df[1,])) * 1}}
  f[1]      <- Model(par[1,], Data)[["LP"]]
  setTxtProgressBar(pb, 1)

  # Start estimation
  for(run in 2:maxit) {
    # Calculate gradient and estimate step parameter
    df[run,] <- grad(Model, par[run-1,], Data, Interval=Interval)
    step[run,] <- suppressWarnings(unlist(optim(step[run-1,], steep, Model=Model,
                                                Data=Data, df=df[run,], est=par[run-1,],
                                                method="Nelder-Mead")$par))
    # Calculate the new parameters' and LogPosterior values
    par[run,] <- par[run-1,] +
      {step[run,1]*df[run,]*{(abs(df[run,]) >  sd(df[run,])) * 1}} +
      {step[run,2]*df[run,]*{(abs(df[run,]) <= sd(df[run,])) * 1}}
    f[run]    <- Model(par[run,], Data)[["LP"]]
    # Update progress bar
    setTxtProgressBar(pb, run)
  }
  close(pb)
  cat("\n")

  # Final messages
  stopTime = proc.time()
  elapsedTime = stopTime - startTime
  cat("It took ",round(elapsedTime[3],2)," secs for the run to finish. \n", sep="")
  plot.ts(f, ylab="LogPosterior", xlab="Iteration",
          xlim=c(1,maxit), ylim=c(min(f), max(f)))

  # Return results
  Results <- list("LogPosterior"=f, "Estimates"=par, "Gradients"=df,
                  "MAP"=par[which.max(f),], "step"=step)
  return(Results)
}
