**evmix_fit**
-------------

This code is a simple interface for fitting and doing MCMC (Bayesian
estimation) on some extreme value mixture models that are provided by the
`evmix` R package. 

The `evmix` package is an R package for fitting extreme value mixture models,
which was developed by Carl Scarrott and Yang Hu, University of Canterbury. For
more information on the evmix package, see
https://CRAN.R-project.org/package=evmix and
http://www.math.canterbury.ac.nz/~c.scarrott/evmix/ .

The code here includes some fixes and speed-work-arounds to `evmix` routines
(thanks to evmix developer Carl Scarrot for assistence), as at the time of
writing the updates were not available on CRAN.

The code here only applies to fitting Gamma-GPD or Normal-GPD mixture models,
with a GPD upper tail. Other models found in the `evmix` package are not used
herein - simply because the latter models were sufficient for our application.

Beware that often extreme value mixture models are hard to fit (e.g. see
discussion in *Scarrot (2015) Univariate Extreme Value Mixture Modelling, in
the book 'Extreme Value Modeling and Risk Analysis: Methods and
Applications'*). Our wrapper routines attempt to make fitting easier by
trialling many fits, and by using simple Bayesian techniques to sample the
posterior parameter distribution (with uniform priors). While these techniques
worked well in our applications, in general care is required, and you should always
use graphical methods to check the fit.

In some applications, for a reasonable model fit it may be necessary to use
the Bayesian priors to restrict the model parameters. For example, for some
datasets the optimal (ML) GPD `threshold` parameter occurs at a very high data
quantile, meaning very little data is used to fit the upper tail model. In this
case ML can lead to an erratic fit. As a work-around, one may use Bayesian
priors to 'force' the threshold parameter to be below some data quantile, thus
ensuring that the upper tail model is fit with sufficient data. 


**USAGE**

To install the code, you need to have the R packages `evmix` and `MCMCpack` installed.
If you don't have them already, this can be achieved by starting R and running
the following commands:
```r
    install.packages(c('evmix', 'MCMCpack'))
```

The codes in 'evmix_fit.R' include some inline doxygen documentation. See the
codes in '../../Analysis' for examples of their usage in our storm wave
clustering code. See the codes in [test_evmix_fit.R](test_evmix_fit.R) for
simpler examples.

**TESTS**

To test the code, run
```r
    source('test_evmix_fit.R') 
```

from within R. It should print information on a number of tests with several PASS
statements, a few package startup messages, but no FAIL's or other errors. Some
figures are also produced.

Note that these tests take 10s of minutes on my multicore linux machine --
and will take longer on machines running windows, since they are not setup to
run in parallel with that OS. 
