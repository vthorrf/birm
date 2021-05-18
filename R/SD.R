SD <- function(Model, startvalue, Data, Interval=1e-8, maxit=100, tol=1e-3) {

  estimates <- SNGD(Model=Model, startvalue=startvalue, Data=Data,
                    Interval=Interval, maxit=maxit, tol=tol)
  return(Model(estimates$MAP, Data))

}
