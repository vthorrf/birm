MAP <- function(Model, parm, Data, algo=NULL,
                maxit=65000, temp=1e-2, tmax=1,
                REPORT=1000) {

  if (is.null(algo)) algo = "SANN"
  if (algo=="SANN") {
    Results <- SANN(Model, parm, Data, maxit=maxit, temp=temp, tmax=tmax, REPORT=REPORT)
  } else if (algo=="GA") {
    Results <- GA(Model, Data)
  } else stop("Unkown estimation algorithm")

  return(Results)
}
