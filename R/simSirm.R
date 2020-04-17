simSirm <- function(n, v, l=NULL, seed=666, sequence=F, dist="norm"){

  set.seed(seed)
  if (is.null(l)) l <- 2
  if (sequence==T){
    if (dist == "norm") {
      sigma <- exp(seq(length=n,-3,3)); theta <- exp(seq(length=v,-3,3))
    } else if (dist == "beta") {
      sigma <- seq(length=n, 1e-20, (1-(1e-20)))
      theta <- seq(length=v, 1e-20, (1-(1e-20)))
    } else stop("Unknow distribution for parameters :(")

  } else if (sequence == F) {
    if (dist == "norm") {
      sigma <- exp(rnorm(n, 0, 1)); theta <- exp(rnorm(v, 0, 1))
    } else if (dist == "beta") {
      sigma <- rbeta(n, 1, 1); theta <- rbeta(v, 1, 1)
    } else stop("Unknow distribution for parameters :(")
  } else stop("The argument 'sequence' has some problem, mate")

  eta <- outer(sigma, theta, function(r,c) r / (c + r) )
  Sirm <- eta
  for (i in 1:nrow(Sirm)) {
    for (j in 1:ncol(Sirm)) {
      Sirm[i,j] <- rbinom(1, (l-1), eta[i,j])
    }
  }

  colnames(Sirm) <- sapply(1:ncol(Sirm), function(x) paste("Item_",x,sep=""))
  rownames(Sirm) <- sapply(1:nrow(Sirm), function(x) paste("Ind_",x,sep=""))

  Result <- list("data"=Sirm,"abil"=sigma,"diff"=theta)

  return(Result)
}
