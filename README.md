birm: Bayesian Item Response Models
=============

This package is intended to be an "play-and-plug" tool for Bayesian Item Response Modeling. Currently, specific unidimensional IRMs are implemented, but future versions will aim on adding more flexibility to the intermediary and advanced users.

birm depends on LaplacesDemon, an open-source statistical package that is intended to provide a complete environment for Bayesian inference.

This package should be considered experimental. The following models are implemented:

* Rasch model (`rasch` function)
* Partial and Generalized Partial Credit models (`pcm` function)
* Rating and Generalized Rating Scale models (`rsm` function)
* 1, 2, 3, 4, and 5-Parameter Logistic models (`plm` function)
* Strict Item Response model (`sirm` function)
* Bayesian Optimal Scoring (`optscr` function)

All of the following estimation methods implemented in LaplacesDemon are available to use:

* Laplace Approximation (`method = "LA"` argument)
* Variational Bayes (`method = "VB"` argument)
* Hit-And-Run Metropolis (`method = "MCMC"` argument)
* Population Monte Carlo (`method = "PMC"` argument)
* Iterative Quadrature (`method = "IQ"` argument)
* Maximum a Posteriori estimation with Simulated Annealing, Genetic Algorithm, or Steepest Gradient Descent (`method = "MAP"` and `algo="SANN"`, `algo="GA"`, or `algo="SD"`; not originally implemented in LaplacesDemon)

# Installation #

Using the 'remotes' package:

    install.packages("remotes")
    remotes::install_github("vthorrf/birm")
