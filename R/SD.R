SD <- function(Model, startvalue, Data, Interval=1e-8, maxit=100) {

  estimates <- SNGD(Model=Model, startvalue=startvalue, Data=Data,
                    Interval=Interval, maxit=maxit)
  return(Model(estimates$MAP, Data))

}
