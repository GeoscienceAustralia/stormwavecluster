**evmix_fit**
-------------

This code is a simple interface for fitting and doing MCMC on some models in
the `evmix` R package, as required in our storm wave analysis. 

The `evmix` package is an R package for fitting extreme value mixture models, which
was developed by Carl Scarrott and Yang Hu, University of Canterbury. See
https://CRAN.R-project.org/package=evmix

The code here also includes some (small) bug-fixes and speed-work-arounds to
`evmix` routines (thanks to evmix developer Carl Scarrot for assistence). We
have been informed the updates might be slow in coming to CRAN (the main
download site) -- hence why some routines are replicated here. 

This code only applies to fitting gamma-gpd or normal-gpd mixture models for
which the gpd models the upper tail. The latter models were sufficient for our
storm wave applications -- and we are not aiming to provide more generic tools.

If you need to fit extreme value mixture models, we suggest you consult the
`evmix` package directly, as well as the relevant literature, where the subject
is treated more generically than here. On the other hand, the codes here might
suit your needs if your application is very similar to ours.

Beware that often extreme value mixture models are hard to fit (e.g. see
discussion in Scarrot (2015) Univariate Extreme Value Mixture Modelling, in the
book 'Extreme Value Modeling and Risk Analysis: Methods and Applications'). In
our experience the initial log-likelihood optimization will often not converge
to a global maxima, and it may be necessary to tweak the optimization
parameters or constrain the parameter ranges to ensure optimal and sensible
results. Bayesian methods seem like a good option for exploring the likelihood
and bypassing issues of local optima, but again care is required to ensure
convergence. As an example of why parameter constraints may be needed, for some
datasets the optimal (ML) gpd threshold parameter occurs at a very high data
quantile, meaning very little data is used to fit the gpd, so the fit can be
erratic. The methods implemented herein have worked well with some examples in
our study, but care was required. 


**USAGE**

To install the code, you need to have the R packages `evmix` and `MCMCpack` installed.
If you don't have them already, this can be achieved by starting R and running
the following commands:

    install.packages(c('evmix', 'MCMCpack'))

The codes in 'evmix_fit.R' include some inline doxygen documentation. See the
codes in '../../Analysis' for examples of their usage in our storm wave
clustering code. See the codes in 'test_evmix_fit.R' for simpler examples.

**TESTS**

To test the code, run

    source('test_evmix_fit.R') 

from within R. It should print information on a number of tests with several PASS
statements, a few package startup messages, but no FAIL's or other errors. Some
figures are also produced.

Note that these tests take 10s of minutes on my multicore linux machine --
and will take longer on machines running windows, since they are not setup to
run in parallel with that OS. 
