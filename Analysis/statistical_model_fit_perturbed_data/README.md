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
1. Copy the code [statistical_model_fit/statistical_model_storm_timings.Rmd](statistical_model_fit/statistical_model_storm_timings.Rmd) to the current directory
2. Within the copied file, find the line:
```r
    # Optionally remove ties in the event statistics by jittering
    break_ties_with_jitter = FALSE
```
and change it to read:
```r
    # Optionally remove ties in the event statistics by jittering
    break_ties_with_jitter = TRUE
```
* This change means that the code will perturb the `event_statistics` before
fitting the model. The idea is to run this code with many different perturbations
to the data, and check the extent to which the statisical model fit is affected.
3. Run the following from within R:
```r
    source('batch_run.R')
```
This assumes you are using a multi-core linux machine. It takes a few hours
on my 6 core ubuntu desktop.
