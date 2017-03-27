# **nhpoisp**
-----------

Code for fitting and simulating non-homogeneous Poisson processes using maximum
likelihood.

Series with gaps between events (e.g. to prevent storm overlap) are supported.
The event rate function can also depend on time and the time since the last event.

The main functions are `rnhpoisp` for simulating synthetic series, and
`fit_nhpoisp` for fitting models to data. 

## **Usage**
------------

An example illustrating some features is provided below. It shows how to:
* Define an event rate function, consisting of a constant rate, a sinusoidal seasonal component, and an exponential clustering term.
* Simulate a random time-series from the rate function, here of 50 years duration
* Estimate the model parameters from data

For more details/options see the in-line documentation in
[nhpoisp.R](nhpoisp.R) (this follows the doxygen style), and look at the tests.

**Below, we make a rate function `lambda(t, tlast=-Inf)` which depends on the time of 
year `t` and the time since the last event `tlast`**. It has a sinusoidal
variation through the year, with a greatly enhanced rate of events just after
an event occurs (often termed 'clustering').
```r
nhp = new.env()
source('nhpoisp.R', local=nhp)

# Define the rate function 'lambda(t, tlast)'. 
#
# To specify lambda, we need a 'rate equation' with parameters 'theta'
#
rate_equation = 'theta[1] + theta[2]*sin(2*pi*(t+theta[3])) + theta[4]*exp(-(t-tlast)/theta[5])'
theta_par = c(1, 2, 0.1, 15.0, 1/50)
#
# Note that lambda must always be >= 0, which we enforce with the
# 'minimum_rate' argument. (By default this is 0)
#
lambda = nhp$get_lambda_function(theta_par, rate_equation, minimum_rate=0)

# Example rate computation
lambda(t=1.3) 
#[1] 2.175571

# If 'tlast' is not provided, it is assumed to be -Inf
# In this example, the rate should increase if 'tlast' is close to 't'
lambda(t=1.3, tlast=1.29)
#[1] 11.27353
```

**Next we simulate a random synthetic timeseries using the above lambda function.**
The main function for this is `rnhpoisp` (for more information see
documentation in the function header).
```r  
set.seed(1) # Make the example reproducible

# 50 year series
series_duration = 50

synthetic_data = nhp$rnhpoisp(duration=series_duration, lambda = lambda)

# 'synthetic_data' is a vector of increasing times at which events occur. 
#
head(synthetic_data)
#[1] 0.1168604 0.1771514 0.8794946 2.2488395 2.2669781 4.1100736

# The clustering in the series is manifest in the distribution of the time
# between events. A considerable fraction of event spacings are much smaller
# than the median value, indicating sequences of 'clustered' events.
#
quantile(diff(synthetic_data), p=c(0.01, 0.1, 0.5, 0.9, 0.99))
#       1%         10%         50%         90%         99% 
# 0.001318402 0.015139788 0.381288058 1.364453604 2.532100262 
```

**Next we back-estimate the parameters of lambda from the synthetic_data series.**
Having good starting parameters is important for getting the fit to converge.
Note that if we were fitting real data, then `synthetic_data` would be read
from a file (e.g. using `scan` or `read.table`). However, here we fit to the
synthetic data simulated above, to show that we can back-estimate the known
paramters.
```r
model_fit = nhp$fit_nhpoisp(
    synthetic_data, 
    rate_equation=rate_equation,
    minimum_rate=0.,
    initial_theta=c(1., 0.2, 0., 10., 1/100),
    integration_dt = 1.0e-04,
    ##
    ## The arguments below control details of the optimization
    ## It can be difficult to fit these models with complex rate functions,
    ## so adjustment may be required.
    ## Consult the code documentation and see help on R's "optim"
    ## function for details.
    ##
    number_of_passes=1,
    enforce_nonnegative_theta=TRUE,
    optim_method=c('Nelder-Mead'),
    optim_control=list(maxit = 500, parscale=c(1, 1, 0.1, 1, 0.01)),
    verbose=FALSE)

# The fitted parameters should approximate the true ones,
# considering their standard errors
model_fit$par  
#[1] 1.17433808 2.51362821 0.11083784 8.44613070 0.01482651

nhp$get_fit_standard_errors(model_fit)
# [1] 0.19456223 0.45338702 0.01351391 5.62600650 0.01261502

# Note these standard errors are approximate only (based on inverting the
# hessian of the likelihood, so valid as the amount of data --> Inf). 

# For comparison, recall the true theta_par parameters were set above as
## theta_par = c(1, 2, 0.1, 15.0, 1/50)

```

**TESTS**

To run the tests, open R and do

```r
source('test_nhpoisp.R')
```

This will take awhile (say 30min on my 6 core linux box). It should
print information on each set of tests which are being run, along with
many 'PASS' statements.

**TO-DO**

This code could be converted to a stand-alone R package with relatively little work, 
and is generically useful. 

Note that the interface is flexible but a bit non-standard (particularly the way the rate
function is often passed as a character string). Consider revising.

It could be adapted to allow more complex rate functions (e.g. depending on the
time of many previous events -- which is relevant for earthquake clustering).
But this would require a bit of work.

**SEE ALSO**

There is an R package on CRAN called nhpoisson which is also for fitting
non-homogeneous Poisson processes. It was released after the code here was
developed, and seems to have different functionality to the current package
(e.g. when I last looked it did not seem to treat the clustering that we
required for the storm analysis, although that may have changed). It's
interface is very different to that of the current package. Anyway it is
certainly worth investigating if you are working on non-homogeneous Poisson
processes.
