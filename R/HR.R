HR <- function(Model, startvalue, Data, Interval=1e-8, maxit=100) {

  estimates <- HRGD(Model=Model, startvalue=startvalue, Data=Data,
                    Interval=Interval, maxit=maxit)
  Results   <- list("Model"=Model(estimates$MAP, Data), "Fit"=estimates)
  return(Results)
}
