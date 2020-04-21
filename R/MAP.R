MAP <- function(Model, parm, Data, maxit=65000, temp=1e-2, tmax=1, REPORT=1000) {
  print("Initiating MAP estimation with SANN algorithm")
  Sys.sleep(.25)
  LogPost <- function(para, Model, Data) {
    LP <- Model(para[1:length(Data$parm.names)], Data)
    -LP$LP
  }
  estimates <- optim(parm, LogPost, Data=Data, Model=Model, method = "SANN",
                     control=list(maxit=maxit, temp=temp, tmax=tmax,
                                  trace=2, REPORT=REPORT))
  return(Model(estimates$par, Data))
}
