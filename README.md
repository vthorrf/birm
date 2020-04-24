birm: Bayesian Item Response Models
=============

This package is intended to be an "easy to use tool" for Bayesian Item Response Modeling. Currently, specific IRMs are implemented, but future versions will aim on adding more flexibility to the users.

birm depends on LaplacesDemon.

This package should be considered experimental. Currently, the following models are implemented:

* Rasch model (`rasch` function)
* Partial and Generalized Partial Credit models (`pcm` function)
* Rating and Generalized Rating Scale models (`rsm` function)
* 1, 2, 3 and 4-Parameter Logistic models (`plm` function)
* Strict Item Response model (`sirm` function)
* Bayesian Optimal Scoring (`optscr` function)

All of the following estimation methods implemented in LaplacesDemon are available to use:

* Laplace Approximation (`method = "LA"` argument)
* Variational Bayes (`method = "VB"` argument)
* No-U-Turn Sampler (`method = "MCMC"` argument)
* Population Monte Carlo (`method = "PMC"` argument)
* Iterative Quadrature (`method = "IQ"` argument)
* Maximum a Posteriori estimation with Simulated Annealing or Genetic Algorithm (`method = "MAP"` and `algo=SANN` or `algo=GA`; not originally implemented in LaplacesDemon)

# Installation #

Using the 'remotes' package:

    install.packages("remotes")
    remotes::install_github("vthorrf/birm")
