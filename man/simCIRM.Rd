\name{scoring}
\alias{scoring}
\title{Scoring Methods to Change Observed Scores}
\description{
Creates a ordered matrix to be evaluated using Additive Conjoint Measurement Theory.
}
\usage{
scoring(rabil,rdiff,n,v,x,l,p=.5,type="rank")
}
\arguments{
   \item{rabil}{A vector of estimates of person abilities.}
   \item{rdiff}{A vector of estimates of item difficulty.}
   \item{n}{The number of individuals in the sample.}
   \item{v}{The number of items in the test.}
   \item{x}{A long format dataset. Only used for type "bayes".}
   \item{l}{The number of levels on the scale. Only used for type "bayes".}
   \item{p}{The probability of answering right when ability and difficulty match. Used for types "cond" and "condb".}
   \item{type}{The method of scoring. "rank" if the ranking method is prefered. "cond" if conditional probability method is prefered. "bayes" if Bayesian Rasch method is prefered. "condb" if Bayesian Conditional is prefered. "2PL" if Two-parameter logistic is prefered.}
}
\value{
A list containing the following component (for "rank" and "cond" methods):
  \item{scores}{A rescored data set.}
For "bayes", "condb" and "2PL" methods:
  \item{matrix}{A rescored data set.}
  \item{abil}{Aptitude estimates.}
  \item{easi}{Easiness estimates.}
  \item{dic}{DIC fit.}
  \item{full}{The full Bayesian run.}
For "2PL" method:
  \item{matrix}{A rescored data set.}
  \item{abil}{Aptitude estimates.}
  \item{easi}{Easiness estimates.}
  \item{disc}{Discrimination estimates.}
  \item{dic}{DIC fit.}
  \item{full}{The full Bayesian run.}
}
\examples{
### Random Data with 100 observations and 10 variables with 5 levels of response
n = 100; v = 10; l = 5
Data <- matrix(NA,nrow=n,ncol=v)
psi <- seq(length=n,-3,3); theta <- seq(length=v,-3,3); eta <- matrix(NA,ncol=v,nrow=n)
seed <- matrix(c(1:(n*v)),nrow=n,ncol=v)
for(i in 1:n)
{
  for(j in 1:v)
  {
    set.seed(seed[i,j])
    eta[i,j] <- exp(psi[i]-theta[j])/(1+exp(psi[i]-theta[j]))
    Data[i,j] <- rbinom(1,l-1,eta[i,j])
  }
}
rm(eta,seed,i,j)
Subj1 <- vector(mode="character",length=9); Subj2 <- vector(mode="character",length=n)
Item <- vector(mode="character",length=v)
for (i in 1:9) {Subj1[i] <- paste("S",0,i,sep="")}
for (i in 1:n) {Subj2[i] <- paste("S",i,sep="")}; Subj2 <- Subj2[10:n]
Subj <- c(Subj1,Subj2); rm(Subj1,Subj2)
for (i in 1:v) {Item[i] <- paste("I",i,sep="")}
DF <- data.frame(Subj,Data); colnames(DF) <- c("Subj",Item[1:v])
require(tidyr)
theData <- gather(DF, Subj, Correct, I1:I10, factor_key=TRUE)
colnames(theData) <- c("Subj","Item","Correct")
theData <- theData[with(theData, order(Subj)), ]
rownames(theData) <- c(1:(n*v)); head(theData)
theData$Subj = factor( theData$Subj )
theData$Item = factor( theData$Item )
theData$Correct = as.integer( theData$Correct )

### Create Scores
AbilF <- latent.vs(Data,l,what="person",method="freq")
Af <- AbilF$abil
EasiF <- latent.vs(Data,l,what="item",method="freq")
Ef <- EasiF$easi

### Rescore
scoring(n=n,v=v,type="rank")
scoring(rabil=Af,rdiff=1-Ef,type="cond")
scoring(n=n,v=v,x=theData,l=l,type="bayes")
scoring(n=n,v=v,x=theData,l=l,type="condb")
scoring(n=n,v=v,x=theData,l=l,type="2PL")
}
