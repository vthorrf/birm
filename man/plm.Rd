\name{plm}
\alias{plm}
\title{Bayesian of the N-Parameter Logistic Models}
\description{
An Bayesian implementation of the N-Parameter Logistic Models, with N from 1 to 4.
}
\usage{
plm(x, p=1, scaling=1.7, method="VB", Iters=500, Smpl=1000,
    Thin=1, A=500, temp=1e-2, tmax=1, algo="GA", seed=666)
}
\arguments{
   \item{x}{A matrix in the wide format. You should certify that: the score are binary and the minimum possible score is 0; all variables are named; and that are only variables to be modeled.}
   \item{p}{Indicates which model will be fitted. When p = 1, the function returns the One-Parameter Logistic Model (1PLM). When p = 2, the function returns the 2PLM. p = 3 returns the 3PLM and p = 4 returns the 4PLM.}
   \item{scaling}{The scaling parameter for when fitting the one-parameter logistic model (p = 1).}
   \item{method}{The Bayesian estimation method to be used. Defaults to "VB" Variational Bayes. It can also be Laplace Approximation ("LA"), MCMC No-U-Turn Sampler ("MCMC"), Population Monte Carlo ("PMC"), Iterative Quadrature ("IQ"), or Maximum a Posteriori with SANN ("MAP").}
   \item{Iters}{The maximum number of iterations. The default value, 500, is proposed for the default estimation method, "VB". For method "LA", we suggest 100 as the default. For method "MCMC", 20000. For method "PMC", 10. For method "IQ", 100. Finally, for method "MAP", if algo="SANN", than 50000, but if algo="GA", than 500.}
   \item{Smpl}{This argument indicates the number of posterior samples to be taken with sampling importance resampling for the estimation methods "IQ","LA", and "VB". For the "PMC" method, this is the number of samples per mixture component. This argument is not used by both the "MCMC" and the "MAP" methods. The default value, 1000, is already fit for all the estimation methods that use this argument.}
   \item{Thin}{This is the number by which the posterior is thinned. Only "MCMC" and "PMC" estimation methods use this argument.}
   \item{A}{A is used only by the "MCMC" method. It is the number of initial, adaptive iterations to be discarded as burn-in. The value defaults to 500.}
   \item{temp}{Controls the "SANN" algorithm applied with the "MAP" method. It is the starting temperature for the cooling schedule. Defaults to 1e-2.}
   \item{tmax}{Is the number of function evaluations at each temperature for the "SANN" algorithm applied with the "MAP" method. Defaults to 1.}
   \item{algo}{The algorithm to be used if the estimation method is "MAP". Defaults to "GA", for a real-valued based genetic algorithm. The other possible value is "SANN", Simulated Annealing.}
   \item{seed}{The random seed for guaranteeing replicable results. As most examples in LaplacesDemon, it defaults to 666.}
}
\value{
A list containing the following components:
  \item{Data}{The information used for running the model.}
  \item{Model}{The LaplacesDemon model used.}
  \item{Fit}{All the information returned by LaplacesDemon estimation method regarding the fit of the model.}
  \item{abil}{The MAP estimates of the theta parameter.}
  \item{diff}{The MAP estimates of the difficulty parameter.}
  \item{disc}{The MAP estimates of the discrimination parameter. Not returned when p = 1.}
  \item{gues}{The MAP estimates of the guessing parameter. Not returned when p < 3.}
  \item{UpAs}{The MAP estimates of the upper asymptote parameter. Not returned when p != 4.}
  \item{DIC}{A list with the values of DIC, Dbar, and pV.}
}
\examples{
### Random Data with 100 observations and 10 variables with 2 levels of response
n = 100; v = 10; l = 2
Data <- simData(n, v, l)$data
# 1PLM
fit1  <- plm(Data, p=1, scaling=1.7, method="MAP", algo="GA", Iters=250)
# 2PLM
fit2  <- plm(Data, p=2, method="MAP", algo="GA", Iters=250)
# 3PLM
fit3  <- plm(Data, p=3, method="MAP", algo="GA", Iters=250)
# 4PLM
fit4  <- plm(Data, p=4, method="MAP", algo="GA", Iters=250)
}
