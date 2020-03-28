simHIRM <- function(n,v,l,R=666,base=exp(1),standard=F){

  set.seed(R)
  base=base
  if (standard==F){
    psi <- seq(length=n,.0027,.9973); gamma <- seq(length=v,.0027,.9973)
  } else {
    psi <- rbeta(n,1,1); gamma <- rbeta(v,1,1)
  }

  fi <- function(base,x) {
    return( ((base^(x) - base^(-x)) / (base^(x) + base^(-x))) )
  }

  E <- outer(psi, gamma, function(r,c) (c-r) )
  eta <- sapply(1:ncol(E), function(x) ((1 - fi(base,E[,x]))/2) )
  hirm <- eta
  for (i in 1:nrow(hirm)) {
    for (j in 1:ncol(hirm)) {
      hirm[i,j] <- rbinom(1,(l-1),eta[i,j])
    }
  }

  colnames(hirm) <- sapply(1:ncol(hirm), function(x) paste("Item_",x,sep=""))
  rownames(hirm) <- sapply(1:nrow(hirm), function(x) paste("Ind_",x,sep=""))

  Result <- list("data"=hirm,"abil"=psi,"diff"=gamma)

  return(Result)
}
