### Best single step
steepA <- function(Model, par, Data, est, df) {
  -Model(est + {par*df}, Data)[["LP"]]
}
Gammas <- function(Model, par, Data, est, df, vi, si, iter, epsilon) {
  v    <- {{par[1] * vi} + {{1 - par[1]} * df     }}/{1 - {par[1] ^ iter}}
  s    <- {{par[2] * si} + {{1 - par[2]} * {df*df}}}/{1 - {par[2] ^ iter}}
  newV <- est + {{par[3]*v}/{epsilon+sqrt(s)}}
  return(-Model(newV, Data)[["LP"]])
}

### ADAM MAP====
SEGD <- function(Model, startvalue, Data, Interval=1e-6, maxit=100, tol=1e-3) {
  # Opening message
  cat("Steepest Adam Maximum a Posteriori estimation will run for ",
      maxit, " iterations at most.\n\n", sep="")
  startTime = proc.time()

  # Initial settings
  epsilon <- Interval
  pb      <- txtProgressBar(min=0, max=maxit, style=3)
  par     <- s <- v <- G <- array(dim = c(maxit,length(startvalue)))
  GG <- array(dim = c(maxit,3))
  f       <- vector("numeric", length=maxit)
  convergence = F

  # First step
  G[1,] <- grad(Model, startvalue, Data, Interval=Interval)
  alpha <- suppressWarnings(unlist(optim(0, steepA, Model=Model, df=G[1,],
                                         Data=Data, est=startvalue,
                                         method="Nelder-Mead")$par))
  GG[1,]   <- c(.5, .5, alpha)
  v[1,]    <- G[1,]*GG[1,1]
  s[1,]    <- G[1,]*GG[1,2]
  par[1,]  <- startvalue + {alpha*G[1,]}
  f[1]     <- Model(par[1,], Data)[["LP"]]
  setTxtProgressBar(pb, 1)

  # Run ADAM algorithm
  for(i in 2:maxit) {
    G[i,]  <- grad(Model, par[i-1,], Data, Interval=Interval)
    GG[i,] <- suppressWarnings(unlist(ucminf(par=c(0,0,0), fn=Gammas, Model=Model,
                                             df=G[i,], Data=Data, est=par[i-1,],
                                             vi=v[i-1,], si=s[i-1,], iter=i,
                                             epsilon=epsilon)$par))
    if (sum(G[i,]) == 0) {
      GG[i,] <- suppressWarnings(unlist(optim(par=c(0,0,0), fn=Gammas, Model=Model,
                                              df=G[i,], Data=Data, est=par[i-1,],
                                              vi=v[i-1,], si=s[i-1,], iter=i,
                                              epsilon=epsilon)$par))
    }
    gammav  <- GG[i,1]
    gammas  <- GG[i,2]
    alpha   <- GG[i,3]
    v[i,]   <- {{gammav * v[i-1,]} + {{1 - gammav} *  G[i,]       }}/{1 - {gammav^i}}
    s[i,]   <- {{gammas * s[i-1,]} + {{1 - gammas} * {G[i,]*G[i,]}}}/{1 - {gammas^i}}
    par[i,] <- par[i-1,] + {{alpha*v[i,]}/{epsilon+sqrt(s[i,])}}
    f[i]    <- Model(par[i,], Data)[["LP"]]
    if ({abs(c(G[i,] %*% G[i,])) < tol} |
        {abs(c(GG[i,] %*% GG[i,])) < {tol * tol}} |
        {mean(abs(c(G[i,] - G[i-1,]))) < tol}) {
      convergence = T
      setTxtProgressBar(pb, maxit)
      cat("\nConvergence achieved!")
      break
    } else { setTxtProgressBar(pb, i) }
  }
  if (convergence == F) {
    cat("\nConvergence may not have been achieved!")
  }
  if (i < maxit) {
    f   <- f[-c({i+1}:maxit)]
    par <- par[-c({i+1}:maxit),]
    G   <- G[-c({i+1}:maxit),]
    s   <- s[-c({i+1}:maxit),]
    v   <- v[-c({i+1}:maxit),]
    GG  <- GG[-c({i+1}:maxit),]
  }
  close(pb)
  cat("\n")

  # Final messages
  stopTime = proc.time()
  elapsedTime = stopTime - startTime
  cat("It took ",round(elapsedTime[3],2)," secs for the run to finish. \n", sep="")
  #plot.ts(f, ylab="LogPosterior", xlab="Iteration",
  #        xlim=c(1,length(f)), ylim=c(min(f), max(f)))

  # Return results
  Results <- list("LogPosterior"=f, "Estimates"=par, "Gradients"=G, "DSG"=s,
                  "Momentum"=v, "MAP"=par[which.max(f),], "steps"=GG)
  return(Results)
}
