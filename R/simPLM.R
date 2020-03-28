simTwo.PL <- function(n,v,l,R=666,d=1,standard=F){

  set.seed(R)
  if (standard==F){
    sigma <- seq(length=n,-3,3); theta <- seq(length=v,-3,3)
  } else {
    sigma <- rnorm(n,0,1); theta <- rnorm(v,0,1)
  }
  if (length(d) == 1) { discr <- rep(d,v) } else {discr <- d}

  E <- outer(sigma, theta, function(r,c) (r-c) )
  eta <- sapply(1:ncol(E), function(x) plogis(discr[x]*E[,x]))
  TPLM <- eta
  for (i in 1:nrow(TPLM)) {
    for (j in 1:ncol(TPLM)) {
      TPLM[i,j] <- rbinom(1,(l-1),eta[i,j])
    }
  }

  colnames(TPLM) <- sapply(1:ncol(TPLM), function(x) paste("Item_",x,sep=""))
  rownames(TPLM) <- sapply(1:nrow(TPLM), function(x) paste("Ind_",x,sep=""))

  Result <- list("data"=TPLM,"abil"=sigma,"diff"=theta)

  return(Result)
}
