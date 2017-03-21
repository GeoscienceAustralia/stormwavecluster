Monte-carlo data preprocessing
------------------------------

Different perturbations are applied to the `hsig` and measured tide with each
run. In particular, this permits random changes to `hsig` to propagate through
to the event definition, in the event that perturbations are large enough to
change that.

Run with:

    Rscript batch_run_preprocessing.R
