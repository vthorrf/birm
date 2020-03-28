simIRasch <- function(n,v,l,R=666,standard=F){

  set.seed(R)
  if (standard==F){
    sigma <- seq(length=n,.791,14.38); theta <- seq(length=v,.791,14.38)
  } else {
    sigma <- rgamma(n,shape=5,rate=1); theta <- rgamma(v,shape=5,rate=1)
  }

  E <- outer(sigma, theta, function(r,c) ((2*log(r))-(2*log(c))) )
  eta <- sapply(1:ncol(E), function(x) plogis(E[,x]))
  Rasch <- eta
  for (i in 1:nrow(Rasch)) {
    for (j in 1:ncol(Rasch)) {
      Rasch[i,j] <- rbinom(1,(l-1),eta[i,j])
    }
  }

  colnames(Rasch) <- sapply(1:ncol(Rasch), function(x) paste("Item_",x,sep=""))
  rownames(Rasch) <- sapply(1:nrow(Rasch), function(x) paste("Ind_",x,sep=""))

  Result <- list("data"=Rasch,"abil"=sigma,"diff"=theta)

  return(Result)
}
