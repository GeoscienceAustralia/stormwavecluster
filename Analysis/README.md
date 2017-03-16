Storm wave clustering analysis for Old Bar
------------------------------------------

These folders contain code implementing the storm clustering model at Old Bar.

If you have not used them before, then note they should be worked through in order:

1. [preprocessing](preprocessing) contains initial data processing and
exploratory analyses. It must be run first
2. [statistical_model_fit](statistical_model_fit) contains the fit of the
statistical model
3. [statistical_model_fit_perturbed_data](statistical_model_fit_perturbed_data)
contains code to re-run the analysis in
[statistical_model_fit](statistical_model_fit) using 'perturbed' event
statistics data, which removes ties from the dataset. This is useful to check
whether data discretization artefacts are causing problems with the statistical
modelling. However, it is not necessary to run this code to understand the
analysis
