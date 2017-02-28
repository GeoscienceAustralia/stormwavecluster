Testing the sensitivity of the statistical modelling to data truncation
-----------------------------------------------------------------------

This folder contains code to help re-run the models in ../statistical_model,
using 'perturbed' storm event summary statistical data. 

The reason for doing this is that the orginal data is subject to round-off
related truncations [e.g.  duration to the nearest hour], and in some
situations, it is possible for such rounding to lead to problems for statistical
procedures for model selection and parameter estimation.

Thus, it is a good idea to re-run our models with perturbed versions of the
input data, to check whether the fit is robust.

