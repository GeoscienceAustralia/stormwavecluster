Vignette on the statistical modelling
-------------------------------------

This folder demonstrates the fit of the statistical model to the data

To run it, it is essential that all codes in [../preprocessing](../preprocessing) have
been run successfuly. Then, run the following codes in order:

* [statistical_model_storm_timings.Rmd](statistical_model_storm_timings.Rmd) Fitting of non-homogeneous poisson process to the event timngs
```r
    # Start R in a terminal, then do:
    library(knitr)
    knit('statistical_model_storm_timings.Rmd')
```
* [statistical_model_univariate_distributions.Rmd](statistical_model_univariate_distributions.Rmd) Fitting probability distributions to each storm summary statistic -- including consideration of seasonal and ENSO dependence
```r
    # Start R in a terminal, then do:
    library(knitr)
    knit('statistical_model_univariate_distributions.Rmd')
```
* [statistical_model_vine_copula.Rmd](statistical_model_vine_copula.Rmd) Fitting a vine copula to the dependencies in the storm summary statistics
```r
    # Start R in a terminal, then do:
    library(knitr)
    knit('statistical_model_vine_copula.Rmd')
```
