**evmix_fit**
-------------

This code is a simple interface for fitting and doing MCMC on some models in
the exmix R package.

The code here also includes some (small) bug-fixes and speed-work-arounds to evmix
routines (thanks to evmix developer Carl Scarrot for assistence in fixing some
of the bugs).  We have been informed the updates might be slow in coming to
CRAN's evmix (the main download site) -- hence why some evmix routines are
replicated here. 

As well as evmix, this code relies on the packages MCMCpack and parallel.

If you need to fit extreme value mixture models we suggest you consult the
evmix package directly, as well as the relevant literature. Beware that often
these models are hard to fit (e.g. see discussion in Scarrot (2015) Univariate
Extreme Value Mixture Modelling, in the book 'Extreme Value Modeling and Risk
Analysis: Methods and Applications'). In our experience the initial
log-likelihood optimization will often not converge to a global maxima, and it
may be necessary to tweak the optimization parameters or constrain the
parameter ranges to ensure optimal and sensible results. Bayesian methods seem
like a good option for exploring the likelihood and bypassing issues of local
optima, but again care is required to ensure convergence. As an example of why
parameter constraints may be needed, for some datasets the optimal (ML) gpd
threshold parameter occurs at a very high data quantile, meaning very little
data is used to fit the gpd, so the fit can be erratic. The methods implemented
herein have worked well with some examples in our study, but care was required. 


**TESTS**

To test the code, run

    source('test_evmix_fit.R') 

from within R. It should print information on a number of tests with several PASS
statements, a few package startup messages, but no FAIL's or other errors. 
