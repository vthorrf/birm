\name{optscr}
\alias{optscr}
\title{Bayesian Optimal Scoring}
\description{
An Bayesian implementation of the Optimal Scoring prodecure.
}
\usage{
optscr(x, levels=NULL, basis=NULL, knot=10, irf=NULL, method="VB",
       Iters=500, Smpl=1000, Thin=1, A=500, temp=1e-2, tmax=1, seed=666)
}
\arguments{
   \item{x}{A matrix in the wide format. You should certify that: the minimum possible score is 0; all variables are named; and that are only variables to be modeled.}
   \item{levels}{The number of responses keys in the scale. If no value is inserted, the function assumes that the maximum value in x is the number of levels.}
   \item{basis}{What basis expansion to use. It can be "b-spline" or "rademacher". Defaults to "b-spline".}
   \item{knot}{How many knots should the basis expansion have. Defaults to 10.}
   \item{irf}{The Item Response Function that converts the basis expansion to probabilities. Defaults to "logit", but can also be "sirm".}
   \item{method}{The Bayesian estimation method to be used. Defaults to "VB" Variational Bayes. It can also be Laplace Approximation ("LA"), MCMC No-U-Turn Sampler ("MCMC"), Population Monte Carlo ("PMC"), Iterative Quadrature ("IQ"), or Maximum a Posteriori with SANN ("MAP").}
   \item{Iters}{The maximum number of iterations. The default value, 500, is proposed for the default estimation method, "VB". For method "LA", we suggest 100 as the default. For method "MCMC", 20000. For method "PMC", 10. For method "IQ", 100. Finally, for method "MAP", 50000.}
   \item{Smpl}{This argument indicates the number of posterior samples to be taken with sampling importance resampling for the estimation methods "IQ","LA", and "VB". For the "PMC" method, this is the number of samples per mixture component. This argument is not used by both the "MCMC" and the "MAP" methods. The default value, 1000, is already fit for all the estimation methods that use this argument.}
   \item{Thin}{This is the number by which the posterior is thinned. Only "MCMC" and "PMC" estimation methods use this argument.}
   \item{A}{A is used only by the "MCMC" method. It is the number of initial, adaptive iterations to be discarded as burn-in. The value defaults to 500.}
   \item{temp}{Controls the "SANN" algorithm applied with the "MAP" method. It is the starting temperature for the cooling schedule. Defaults to 1e-2.}
   \item{tmax}{Is the number of function evaluations at each temperature for the "SANN" algorithm applied with the "MAP" method. Defaults to 1.}
   \item{seed}{The random seed for guaranteeing replicable results. As most examples in LaplacesDemon, it defaults to 666.}
}
\value{
A list containing the following components:
  \item{Data}{The information used for running the model.}
  \item{Model}{The LaplacesDemon model used.}
  \item{Fit}{All the information returned by LaplacesDemon estimation method regarding the fit of the model.}
  \item{abil}{The MAP estimates of the theta parameter.}
  \item{icc}{A list with the estimated ICCs of each item.}
  \item{DIC}{A list with the values of DIC, Dbar, and pV.}
}
\examples{
### Random Data with 100 observations and 10 variables with 5 levels of response
n = 100; v = 10; l = 5
Data <- simData(n, v, l)$data
fit  <- optscr(Data, levels=l, method="LA", Iters=100, Smpl=1000)
}
