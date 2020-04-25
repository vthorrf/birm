simRasch <- function(n, v ,l=NULL, seed=666, sequence=F, dist="norm",
                     interaction=F, weight="normal"){

  set.seed(seed)
  if (is.null(l)) l <- 2
  if (sequence==T){
    if (dist == "norm") {
      sigma <- seq(length=n,-3,3); theta <- seq(length=v,-3,3)
    } else if (dist == "beta") {
      sigma <- seq(length=n, .0027,.9973)
      theta <- seq(length=v, .0027,.9973)
    } else stop("Unknow distribution for parameters :(")

  } else {
    if (dist == "norm") {
      sigma <- rnorm(n, 0, 1); theta <- rnorm(v, 0, 1)
    } else if (dist == "beta") {
      sigma <- rbeta(n, 1, 1); theta <- rbeta(v, 1, 1)
    } else stop("Unknow distribution for parameters :(")
  }

  if (interaction==F) {
    E <- outer(sigma, theta, function(r,c) (r - c) )
  } else {
    if (weight=="normal") {
      weight <- rnorm(1)
    } else if (weight=="pos") {
      weight <- abs(rnorm(1))
    } else if (weight=="neg") {
      weight <- -abs(rnorm(1))
    } else stop("Unknow procedure.")
    E <- outer(sigma, theta, function(r,c) (r + c + (weight*r*c)) )
  }
  eta <- sapply(1:ncol(E), function(x) plogis(E[,x]))
  Rasch <- eta
  for (i in 1:nrow(Rasch)) {
    for (j in 1:ncol(Rasch)) {
      Rasch[i,j] <- rbinom(1,(l-1),eta[i,j])
    }
  }

  colnames(Rasch) <- sapply(1:ncol(Rasch), function(x) paste("Item_",x,sep=""))
  rownames(Rasch) <- sapply(1:nrow(Rasch), function(x) paste("Ind_",x,sep=""))

  if (interaction == F) {
    Result <- list("data"=Rasch,"abil"=sigma,"diff"=theta)
  } else {
    Result <- list("data"=Rasch,"abil"=sigma,"diff"=theta,"weight"=weight)
  }

  return(Result)
}
