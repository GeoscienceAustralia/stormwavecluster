Testing the sensitivity of the statistical modelling to data truncation
-----------------------------------------------------------------------

This folder contains code to help re-run the models in [../statistical_model_fit](../statistical_model_fit),
using 'perturbed' storm event summary statistical data. 

The reason for doing this is that the orginal data is subject to round-off
related truncations [e.g. storm duration is only recorded to the nearest hour,
because the underlying data are hourly]. In some situations, it is possible for
such rounding to lead to problems for statistical procedures for model
selection and parameter estimation, if those procedures implicitly assume continuous
input data [for which the probability of 'ties' should be zero]

Thus, it is a good idea to re-run our models with perturbed versions of the
input data, to check whether the fit is robust.


USAGE
-----

To run 100 models with different random perturbations to the data, call these
commands in succession [note: it is suggested to run on a multi-core linux machine,
since the code will only use parallel in that case]. 

    Rscript batch_run_storm_timings.R

    Rscript batch_run_univariate_distributions.R

    Rscript batch_run_vine_copula.R

Beware the run time will be around 100x as long as it took to run the code in
[../statistical_model_fit](../statistical_model_fit), unless you have a
multi-core machine running linux [in which case the work is distributed]. To save
time, the Bayesian model fits only run 1 chain when using perturbed data (since
our earlier checks indicated that the convergence of all chains was already quite good).

Code to check the outputs (and produce a markdown report) can be run from
within R with:
```r
library(knitr)
knit('check_storm_timings.Rmd')
```
Proceed similarly for the other check_XXXX.R routines.

The results are in [check_storm_timings.md](check_storm_timings.md); [check_univariate_dist.md](check_univariate_dist.md) ; [check_vine_copula.md](check_vine_copula.md)


