\name{optscr}
\alias{optscr}
\title{Bayesian Optimal Scoring}
\description{
An Bayesian implementation of the Optimal Scoring prodecure.
}
\usage{
optscr(x, levels=NULL, M=5, basis="rademacher", err=NULL, knots=NULL,
       degree=3, method="LA", Iters=100, Smpl=1000, Thin=1, a.s=0.234,
       B=TRUE, temp=1e-2, tmax=1, algo="GA", seed=666, Interval=1e-8)
}
\arguments{
   \item{x}{A matrix in the wide format. You should certify that: the minimum possible score is 0; all variables are named; and that are only variables to be modeled.}
   \item{M}{The number of hidden nodes for the single layer when the argument "basis"" is equal to "nn". Defaults to 5.}
   \item{levels}{The number of responses keys in the scale. If no value is inserted, the function assumes that the maximum value in x is the number of levels.}
   \item{basis}{What basis expansion to use. It can be "bs", for b-spline, "rademacher", "legendre", or "nn", for neural network. Defaults to "rademacher".}
   \item{err}{A vector of error estimates for true scores. Defaults to 1e-4 for all scores.}
   \item{knots}{A vector of knots to be used with the basis expansion. If no knots are provided, a default vector, based on the unique values observed from the sum scores, is used.}
   \item{degree}{Degree of the piecewise ("bs") or orthogonal ("legendre") polynomial basis. Default is 3.}
   \item{method}{The Bayesian estimation method to be used. Defaults to "LA" Variational Bayes. It can also be Laplace Approximation ("LA"), MCMC Hit-And-Run Metropolis ("MCMC"), Population Monte Carlo ("PMC"), Iterative Quadrature ("IQ"), or Maximum a Posteriori with SANN ("MAP").}
   \item{Iters}{The maximum number of iterations. The default value, 100, is proposed for the default estimation method, "LA". For method "VB", we suggest 500 as the default. For method "MCMC", 20000. For method "PMC", 10. For method "IQ", 100. Finally, for method "MAP", if algo="SANN", than 50000, if algo="GA", than 500, and if algo="SD", than 100. However, ALWAYS check for convergence.}
   \item{Smpl}{This argument indicates the number of posterior samples to be taken with sampling importance resampling for the estimation methods "IQ","LA", and "VB". For the "PMC" method, this is the number of samples per mixture component. This argument is not used by both the "MCMC" and the "MAP" methods. The default value, 1000, is already fit for all the estimation methods that use this argument.}
   \item{Thin}{This is the number by which the posterior is thinned. Only "MCMC" and "PMC" estimation methods use this argument.}
   \item{a.s}{This is the alpha.star parameter, representing the target acceptance rate of the "MCMC" method. The value defaults to 0.234.}
   \item{B}{A logical parameter representing if blockwise sampling (i.e., true scores first, then item parameters) should be used instead of multivariate sampling when using the "MCMC" method. The value defaults to TRUE.}
   \item{temp}{Controls the "SANN" algorithm applied with the "MAP" method. It is the starting temperature for the cooling schedule. Defaults to 1e-2.}
   \item{tmax}{Is the number of function evaluations at each temperature for the "SANN" algorithm applied with the "MAP" method. Defaults to 1.}
   \item{algo}{The algorithm to be used if the estimation method is "MAP". Defaults to "GA", for a real-valued based genetic algorithm. The other possible values are "SANN", for Simulated Annealing, and "SD", for Steepest Gradient Descent.}
   \item{seed}{The random seed for guaranteeing replicable results. As most examples in LaplacesDemon, it defaults to 666.}
   \item{Interval}{This is the interval of numeric differencing when using the "SD" algorithm with the "MAP" method. Defaults to 1e-8.}
}
\value{
A list containing the following components:
  \item{Data}{The information used for running the model.}
  \item{Model}{The LaplacesDemon model used.}
  \item{Fit}{All the information returned by LaplacesDemon estimation method regarding the fit of the model.}
  \item{abil}{The MAP estimates of the theta parameter.}
  \item{icc}{A list with the estimated ICCs of each item.}
For all estimation methods but "MAP", it is also returned:
  \item{DIC}{A list with the values of DIC, Dbar, and pV.}
When "MAP" is used, instead of DIC, it is returned:
  \item{FitIndexes}{A list with the values of AIC, BIC, CAIC, and SABIC.}
}
\examples{
### Random Data with 100 observations and 10 variables with 5 levels of response
n = 100; v = 10; l = 5
Data <- simData(n, v, l)$data
fit  <- optscr(Data, levels=l, method="MAP", algo="GA", Iters=250)
}
