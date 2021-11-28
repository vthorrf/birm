\name{plm}
\alias{plm}
\title{Bayesian of the N-Parameter Logistic Models}
\description{
A Bayesian implementation of the N-Parameter Logistic Models, with N from 1 to 5.
}
\usage{
plm(x, p=1, scaling=1.7, method="LA", Iters=100, Smpl=1000, Thin=1, a.s=0.234,
    temp=1e-2, tmax=NULL, algo="GA", seed=666, Interval=1e-8)
}
\arguments{
   \item{x}{A matrix in the wide format. You should certify that: the score are binary and the minimum possible score is 0; all variables are named; and that are only variables to be modeled.}
   \item{p}{Indicates the number of item parameters in the model, being the models cumulative (i.e., p = 4 has all the parameters of p = 3 plus a new one). When p = 1, the function returns only the difficulties estimate. When p = 2, also returns the estimates of discrimination per item. p = 3 will return the guessing estimates per item. p = 4 is the upper-asymptote of each item. Finally, p = 5 returns the curve-asymetry parameter for each item.}
   \item{scaling}{The scaling parameter for when fitting the one-parameter logistic model (i.e., when p = 1).}
   \item{method}{The Bayesian estimation method to be used. Defaults to Laplace Approximation, "LA". It can also be Variational Bayes ("VB"), MCMC Hit-And-Run Metropolis ("MCMC"), Population Monte Carlo ("PMC"), Iterative Quadrature ("IQ"), or Maximum a Posteriori ("MAP").}
   \item{Iters}{The maximum number of iterations. The default value, 100, is proposed for the default estimation method, "LA". For method "VB", we suggest 500 as the default. For method "MCMC", 20000. For method "PMC", 10. For method "IQ", 100. Finally, for method "MAP", if algo="SANN", than 50000, if algo="GA", than 500, and if algo="SD" or "ADAM", than 100. However, ALWAYS check for convergence.}
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
  \item{disc}{The MAP estimates of the discrimination parameter.}
  \item{gues}{The MAP estimates of the guessing parameter. Not returned when p < 3.}
  \item{UpAs}{The MAP estimates of the upper asymptote parameter. Not returned when p < 4.}
  \item{asym}{The MAP estimates of the curve asymetry parameter. Not returned when p < 5.}
For all estimation methods but "MAP", it is also returned:
  \item{DIC}{A list with the values of DIC, Dbar, and pV.}
When "MAP" is used, instead of DIC, it is returned:
  \item{FitIndexes}{A list with the values of AIC, BIC, CAIC, and SABIC.}
}
\examples{
### Random Data with 100 observations and 10 variables with 2 levels of response
n = 100; v = 10; l = 2
Data <- simData(n, v, l)$data
# 1PLM
fit1  <- plm(Data, p=1, scaling=1.7, method="MAP", algo="GA", Iters=50)
# 2PLM
fit2  <- plm(Data, p=2, method="MAP", algo="GA", Iters=50)
# 3PLM
fit3  <- plm(Data, p=3, method="MAP", algo="GA", Iters=50)
# 4PLM
fit4  <- plm(Data, p=4, method="MAP", algo="GA", Iters=50)
# 5PLM
fit5  <- plm(Data, p=5, method="MAP", algo="GA", Iters=50)
}
