\name{simData}
\alias{simData}
\title{Simulate Item Response Data}
\description{
Simulate data from many different IRMs.
}
\usage{
simData(n=NULL, v=NULL ,l=NULL, p=1, scaling=1.7, seed=666, sequence=F,
        model="rasch", dist="norm", interaction=F, weight="normal")
}
\arguments{
   \item{n}{Sample size.}
   \item{v}{The number of variables.}
   \item{l}{The number of response keys.}
   \item{p}{Indicates which parameter-logistic model will be fitted if the argument model is equal to "plm". When p = 1, the function returns the One-Parameter Logistic Model (1PLM). When p = 2, the function returns the 2PLM. p = 3 returns the 3PLM and p = 4 returns the 4PLM.}
   \item{scaling}{The scaling parameter for when simulating the one-parameter logistic model (p = 1 and model="plm").}
   \item{seed}{The random seed for guaranteeing replicable results. As most examples in LaplacesDemon, it defaults to 666.}
   \item{sequence}{If the true scores and true difficulties should be a sequence (sequence = T) or if they should be random values (sequence = F) according to the base distribution (the value of the "dist" argument).}
   \item{model}{Which model to use to simulate the data. It defaults to the rasch model, with value "rasch". It can also be the SIRM, "sirm", or the Parameter-Logistic Models, "plm".}
   \item{dist}{Which distribution to base the sequence or random draw of the parameters. It defaults to "norm", the normal distribution with mean 0 and sd 1. It can also be "beta", the beta distribution with parameters 1 and 1.}
   \item{interaction}{Set this argument to TRUE if a simulated dataset with an interaction term is desired. Defaults to FALSE and only used if model = "rasch".}
   \item{weight}{If interaction is set to TRUE, then a weight term will be added to the model. The default, "normal", draws the weight parameter from a 0 mean, 1 sd normal distribution. It can also be "pos" or "neg", returning a positive or negative value, respectively, draw from a truncated 0 mean, 1 sd normal distribution. Only used if model = "rasch".}
}
\value{
A list containing a dataset and the parameters used to create the dataset.
}
\examples{
### Random Data with 100 observations and 10 variables with 2 levels of response
n = 100; v = 10; l = 2
RaschData <- simData(n, v, l, model="rasch")

Data <- RaschData$data
Abil <- RaschData$abil
Diff <- RaschData$diff

}
