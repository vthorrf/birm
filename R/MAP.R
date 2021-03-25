MAP <- function(Model, parm, Data, algo=NULL,
                maxit=65000, temp=1e-2, tmax=1,
                REPORT=1000, Interval=1e-8) {

  if (is.null(algo)) algo = "GA"
  if (algo=="SANN") {
    Results <- SANN(Model, parm, Data, maxit=maxit, temp=temp, tmax=tmax, REPORT=REPORT)
  } else if (algo=="GA") {
    Results <- GA(Model, parm, Data, maxit=maxit)
  } else if (algo=="SD") {
    Results <- SD(Model, parm, Data, maxit=maxit, Interval=Interval)
  } else  stop("Unkown estimation algorithm")

  Results$FI <- list("AIC"=(2 * length(parm)) + Results$Dev,
                     "BIC"=(log(Data$n) * length(parm)) + Results$Dev,
                     "CAIC"=((log(Data$n) + 1) * length(parm)) + Results$Dev,
                     "SABIC"=(log((Data$n + 2)/24) * length(parm)) + Results$Dev)
  return(Results)
}
