simPLM <- function(n, v, l=NULL, p=1, scaling=1.7, seed=666,
                   sequence=F, dist="norm"){

  set.seed(seed)
  if (is.null(l)) l <- 2
  if (sequence==T){
    if (dist == "norm") {
      sigma <- seq(length=n,-3,3); theta <- seq(length=v,-3,3)
    } else if (dist == "beta") {
      sigma <- seq(length=n, .0027,.9973)
      theta <- seq(length=v, .0027,.9973)
    } else stop("Unknow distribution for parameters :(")

  } else if (sequence == F) {
    if (dist == "norm") {
      sigma <- rnorm(n, 0, 1); theta <- rnorm(v, 0, 1)
    } else if (dist == "beta") {
      sigma <- rbeta(n, 1, 1); theta <- rbeta(v, 1, 1)
    } else stop("Unknow distribution for parameters :(")
  } else stop("The argument 'sequence' has some problem, mate")

  if (p == 1) {
    eta <- outer(sigma, theta, function(r,c) plogis(scaling * (r - c)) )
  } else if (p == 2) {
    disc <- rnorm(n, 0, 1)
    E <- outer(sigma, theta, function(r,c) r - c )
    eta <- sapply(1:ncol(E), function(x) plogis(disc[x] * E[,x]))
  } else if (p == 3) {
    disc <- rnorm(v, 0, 1)
    ch   <- rbeta(v, 1, 1)
    E <- outer(sigma, theta, function(r,c) r - c )
    eta <- sapply(1:ncol(E), function(x) plogis(
      ch[x] + ( (1 - ch[x]) / (1 + exp(-disc[x] * ( E[,x] ))) )
    ) )
  } else if (p == 4) {
    disc <- rnorm(v, 0, 1)
    ch   <- rbeta(v, 1, 1)
    UA   <- runif(v, .5, 1)
    E <- outer(sigma, theta, function(r,c) r - c )
    eta <- sapply(1:ncol(E), function(x) plogis(
      ch[x] + ( (UA[x] - ch[x]) / (1 + exp(-disc[x] * ( E[,x] ))) )
    ) )
  } else stop("Unknow model :(")

  PLM <- eta
  for (i in 1:nrow(PLM)) {
    for (j in 1:ncol(PLM)) {
      PLM[i,j] <- rbinom(1, (l-1), eta[i,j])
    }
  }

  colnames(PLM) <- sapply(1:ncol(PLM), function(x) paste("Item_",x,sep=""))
  rownames(PLM) <- sapply(1:nrow(PLM), function(x) paste("Ind_",x,sep=""))

  if (p == 1) {
    Result <- list("data"=PLM,"abil"=sigma,"diff"=theta, "scaling"=scaling)
  } else if (p == 2) {
    Result <- list("data"=PLM,"abil"=sigma,"diff"=theta, "disc"=disc)
  } else if (p == 3) {
    Result <- list("data"=PLM,"abil"=sigma,"diff"=theta, "disc"=disc,
                   "c"=ch)
  } else if (p == 4) {
    Result <- list("data"=PLM,"abil"=sigma,"diff"=theta, "disc"=disc,
                   "c"=ch, "UA"=UA)
  } else stop("Unknow value por 'p', bro")

  return(Result)
}
