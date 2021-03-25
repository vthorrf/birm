GA <- function(Model, parm, Data, maxit) {
  print("Initiating MAP estimation with GA algorithm")
  #want = c("GA")
  #have = want %in% rownames(installed.packages())
  #if ( any(!have) ) { install.packages( want[!have] ) }
  #suppressMessages(require(GA))

  LogPost <- function(para) {
    LP <- Model(para[1:length(Data$parm.names)], Data)
    LP$LP
  }
  estimates <- ga("real-valued", fitness = LogPost,
                  lower=rep(-5,length(Data$parm.names)),
                  upper=rep(5,length(Data$parm.names)),
                  names=Data$parm.names, monitor=T,
                  maxiter=maxit, suggestions=parm)
  return(Model(estimates@solution, Data))
}
