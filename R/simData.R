simData <- function(n=NULL, v=NULL ,l=NULL, p=1, scaling=1.7,
                    seed=666, sequence=F, model="rasch",
                    dist="norm", interaction=F, weight="normal"){

  set.seed(seed)
  if (model == "rasch") {
    source("simRasch.R")
    Result <- simRasch(n=n, v=v, l=l, seed=seed, sequence=sequence,
                       dist=dist, interaction=interaction, weight=weight)
  } else if (model == "sirm") {
    source("simSirm.R")
    Result <- simSirm(n=n, v=v, l=l, seed=seed, sequence=sequence, dist=dist)
  } else if (model == "plm") {
    source("simPLM.R")
    Result <- simPLM(n=n, v=v, l=l, p=p, scaling=scaling,
                     seed=seed, sequence=sequence, dist=dist)

  } else stop("Unknow model :(")

  return(Result)
}
