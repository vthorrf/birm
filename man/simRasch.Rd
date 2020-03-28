\name{simRasch}
\alias{simRasch}
\title{Simulate Data From the Rasch model}
\description{
Simulate data from the Rasch model.
}
\usage{
simRasch(n, v, l, R=666, standard=F)
}
\arguments{
   \item{n}{Sample size.}
   \item{v}{The number of variables.}
   \item{l}{The number of response keys.}
   \item{R}{The random seed for guaranteeing replicable results. As most examples in LaplacesDemon, it defaults to 666.}
   \item{standard}{If the true scores and true difficulties should be a sequence (standard = F) or if they should be random values draw from 0 mean, 1 sd normal distributions (standard = T).}
}
\value{
A list containing the following components:
  \item{Data}{The dataset.}
  \item{Abil}{The true scores.}
  \item{Diff}{The true difficulty.}
}
\examples{
### Random Data with 100 observations and 10 variables with 5 levels of response
n = 100; v = 10; l = 2
RaschData <- simRasch(n, v, l)

Data <- RaschData$data
Abil <- RaschData$abil
Diff <- RaschData$diff

}
