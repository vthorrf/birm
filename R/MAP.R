MAP <- function(Model, parm, Data, algo=NULL, seed=666,
                maxit=65000, temp=1e-2, tmax=NULL,
                REPORT=1000, Interval=1e-8) {

  if (is.null(algo)) algo = "GA"
  if (algo=="SANN") {
    if (is.null(tmax)) tmax = 1
    Results <- SANN(Model, parm, Data, maxit=maxit, temp=temp, tmax=tmax, REPORT=REPORT)
  } else if (algo=="GA") {
    Results <- GA(Model, parm, Data, maxit=maxit)
  } else if (algo=="SD") {
    if (is.null(tmax)) tmax = 1e-3
    Results <- SD(Model, parm, Data, maxit=maxit, Interval=Interval, tol=tmax)
  } else if (algo=="ADAM") {
    if (is.null(tmax)) tmax = 1e-3
    Results <- SE(Model, parm, Data, maxit=maxit, Interval=Interval, tol=tmax)
  } else stop("Unkown estimation algorithm")

  Results$FI <- list("AIC"=(2 * length(parm)) + Results$Dev,
                     "BIC"=(log(Data$n) * length(parm)) + Results$Dev,
                     "CAIC"=((log(Data$n) + 1) * length(parm)) + Results$Dev,
                     "SABIC"=(log((Data$n + 2)/24) * length(parm)) + Results$Dev)
  return(Results)
}
