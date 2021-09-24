\name{fishburn}
\alias{fishburn}
\title{Bayesian Fishburn Models}
\description{
A Bayesian implementation of each of the five factorial utility models proposed by Franco et al. (under review). These models were inspired by Fishburn (1974).
}
\usage{
fishburn(x, method="LA", model=1, Iters=100, Smpl=1000, Thin=1, a.s=0.234,
         temp=1e-2, tmax=NULL, algo="GA", seed=666, Interval=1e-8)
}
\arguments{
   \item{x}{A matrix in the wide format. You should certify that: the score are binary and the minimum possible score is 0; all variables are named; and that are only variables to be modeled.}
   \item{method}{The Bayesian estimation method to be used. Defaults to "LA" Variational Bayes. It can also be Laplace Approximation ("LA"), MCMC Hit-And-Run Metropolis ("MCMC"), Population Monte Carlo ("PMC"), Iterative Quadrature ("IQ"), or Maximum a Posteriori ("MAP").}
   \item{model}{Indicates which model to be fitted. When p = 1, the Additive model is fitted. When p = 2, the Multiplicative model is fitted. p = 3 will fit the Person-Multiplicative model. p = 4 will fit the Item-Multiplicative model. Finally, p = 5 fits the Additive-Multiplicative model.}
   \item{Iters}{The maximum number of iterations. The default value, 100, is proposed for the default estimation method, "LA". For method "VB", we suggest 500 as the default. For method "MCMC", 20000. For method "PMC", 10. For method "IQ", 100. Finally, for method "MAP", if algo="SANN", than 50000, if algo="GA", than 500, and if algo="SD", than 100. However, ALWAYS check for convergence.}
   \item{Smpl}{This argument indicates the number of posterior samples to be taken with sampling importance resampling for the estimation methods "IQ","LA", and "VB". For the "PMC" method, this is the number of samples per mixture component. This argument is not used by both the "MCMC" and the "MAP" methods. The default value, 1000, is already fit for all the estimation methods that use this argument.}
   \item{Thin}{This is the number by which the posterior is thinned. Only "MCMC" and "PMC" estimation methods use this argument.}
   \item{a.s}{This is the alpha.star parameter, representing the target acceptance rate of the "MCMC" method. The value defaults to 0.234.}
   \item{temp}{Controls the "SANN" algorithm applied with the "MAP" method. It is the starting temperature for the cooling schedule. Defaults to 1e-2.}
   \item{tmax}{It is the number of function evaluations at each temperature for the "SANN" algorithm applied with the "MAP" method. Defaults to 1. It is also the relative convergence tolerance when algo is "SD" or "ADAM". Defaults to 1e-8.}
   \item{algo}{The algorithm to be used if the estimation method is "MAP". Defaults to "GA", for a real-valued based genetic algorithm. The other possible values are "SANN", for Simulated Annealing, "SD", for Steepest 2-Group Gradient Descent, and "ADAM", for Steepest Adaptive Moment Estimation.}
   \item{seed}{The random seed for guaranteeing replicable results. As most examples in LaplacesDemon, it defaults to 666.}
   \item{Interval}{This is the interval of numeric differencing when using the "SD" and "ADAM" algorithms with the "MAP" method. Defaults to 1e-8.}
}
\value{
A list containing the following components:
  \item{Data}{The information used for running the model.}
  \item{Model}{The LaplacesDemon model used.}
  \item{Fit}{All the information returned by LaplacesDemon estimation method regarding the fit of the model.}
  \item{abil}{The MAP estimates of the theta parameter.}
  \item{diff}{The MAP estimates of the difficulty parameter.}
For all estimation methods but "MAP", it is also returned:
  \item{DIC}{A list with the values of DIC, Dbar, and pV.}
When "MAP" is used, instead of DIC, it is returned:
  \item{FitIndexes}{A list with the values of AIC, BIC, CAIC, and SABIC.}
}
\examples{
### Random Data with 100 observations and 10 variables with 2 levels of response
n = 100; v = 10; l = 2
Data <- simData(n, v, l)$data
# 1 - AFM
fit1  <- fishburn(Data, model=1, method="MAP", algo="GA", Iters=250)
# 2 - MFM
fit2  <- fishburn(Data, model=2, method="MAP", algo="GA", Iters=250)
# 3 - P-MFM
fit3  <- fishburn(Data, model=3, method="MAP", algo="GA", Iters=250)
# 4 - I-MFM
fit4  <- fishburn(Data, model=4, method="MAP", algo="GA", Iters=250)
# 5 - AM-FM
fit5  <- fishburn(Data, model=5, method="MAP", algo="GA", Iters=250)
}
