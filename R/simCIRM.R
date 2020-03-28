simCIRM <- function(n,v,l,R=666,p=.5,standard=F){

  p=p
  set.seed(R)
  if (standard==F){
    psi <- seq(length=n,.0027,.9973); gamma <- seq(length=v,.0027,.9973)
  } else {
    psi <- rbeta(n,1,1); gamma <- rbeta(v,1,1)
  }

  eta <- outer(psi, gamma, function(r,c) (p*r)/((p*r)+((1-p)*c)) )
  cirm <- eta
  for (i in 1:nrow(cirm)) {
    for (j in 1:ncol(cirm)) {
      cirm[i,j] <- rbinom(1,(l-1),eta[i,j])
    }
  }

  colnames(cirm) <- sapply(1:ncol(cirm), function(x) paste("Item_",x,sep=""))
  rownames(cirm) <- sapply(1:nrow(cirm), function(x) paste("Ind_",x,sep=""))

  Result <- list("data"=cirm,"abil"=psi,"diff"=gamma)

  return(Result)
}
