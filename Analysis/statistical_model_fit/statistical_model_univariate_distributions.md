
# **Modelling the univariate distributions of storm event statistics**
--------------------------------------------------------------------------

*Gareth Davies, Geoscience Australia 2017*

# Introduction
------------------

This document follows on from
[statistical_model_storm_timings.md](statistical_model_storm_timings.md)
in describing our statistical analysis of storm waves at Old Bar. 

It illustrates the process of fitting probability distributions to the storm event summary statistics,
which are conditional on the time of year and ENSO.

It is essential that the code
[statistical_model_storm_timings.md](statistical_model_storm_timings.md) has
alread been run, and produced an Rdata file
*'Rimages/session_storm_timings_XXXX.Rdata'*, where XXXX contains information on 
whether data perturbation was applied.


Supposing the prerequisites have been run, and you have R installed,
along with all the packages required to run this code, and a copy of the
*stormwavecluster* git repository, then you should be able to re-run the
analysis here by simply copy-pasting the code. Alternatively, it can be run
with the `knit` command in the *knitr* package: 

```r
library(knitr)
knit('statistical_model_univariate_distributions.Rmd')
```

To run the code in tie-breaking mode, be sure to pass the a command-line
argument matching `break_ties` to R when starting, followed by an integer ID > 0,
e.g.

    R --args --break_ties 1234

or

    Rscript script_name_here.R --break_ties 1234

Running the above commands many times is facilitated by scripts in
[../statistical_model_fit_perturbed_data/README.md](../statistical_model_fit_perturbed_data/README.md)

The basic approach followed here is to:
* **Step 1: Load the previous session**
* **Step 2: Exploratory analysis of seasonal non-stationarity in event statistics**
* **Step 3: Model the distribution of each storm summary statistic, dependent on season (and mean annual SOI for wave direction)**

Later we will model the remaining joint dependence between these variables, and
simulate synthetic storm sequences. 

# **Step 1: Load the previous session and set some key parameters**
Here we re-load the session from the previous stage of the modelling. We also
set some parameters controlling the Monte-Carlo Markov-Chain (MCMC) computations 
further in the document. 
* The default parameter values should be appropriate for the analysis
herein. To save computational effort (for testing purposes) users might reduce
the `mcmc_chain_length`. To reduce memory usage, users can increase the
`mcmc_chain_thin` parameter. If using other datasets, it may be necessary to
increase the `mcmc_chain_length` to get convergence.
* The code is also setup to run using a previous session with data ties broken
at random. See
[../statistical_model_fit_perturbed_data/README.md](../statistical_model_fit_perturbed_data/README.md)
for information on how to do this.


```r
# Here we support multiple runs with random tie-breaking of the data
# If R was passed a commandline argument 'break_ties n' on startup (with n = integer),
# then read the n'th R session matching 'Rimages/session_storm_timings_TRUE_*.Rdata'.
# That session will correspond to one of the tie-breaking sessions
if( length(grep('break_ties', commandArgs(trailingOnly=TRUE))) > 0 ){
    # Found --break_ties command line argument

    break_ties_with_jitter = TRUE

    # Read one of the sessions with tie-breaking
    session_n = as.numeric(commandArgs(trailingOnly=TRUE)[2])

    if(session_n < 1) stop('Invalid run ID')

    # In this case, only run 1 mcmc chain on 1 core [since we will check many
    # tie-breaking sessions]
    mcmc_nchains = 1
    mcmc_ncores = 1

    mcmc_chain_length = 2e+06 # One long chain

}else{
    # No tie-breaking -- using 'raw' data

    break_ties_with_jitter = FALSE
    session_n = 0

    # In this case, run more chains in parallel.
    mcmc_nchains = 6
    mcmc_ncores = 6

    # However, the parallel framework used here does not work on windows,
    # so if running windows, only use 1 core
    if(.Platform$OS.type == 'windows') mcmc_ncores = 1

    # Length of each MCMC chain. Should be 'large' e.g 10^6, except for test runs 
    # We run multiple chains to enhance the likelihood of detecting non-convergence
    # since anyway this is cheap in parallel. These are pooled for final estimates,
    # but it is essential to manually check the convergence of the chains [e.g.
    # by comparing high return period confidence intervals].
    mcmc_chain_length = 1e+06 #1e+05 

}

# Make a 'title' which can appear in filenames to identify this run
run_title_id = paste0(break_ties_with_jitter, '_', session_n)

previous_R_session_file = paste0('Rimages/session_storm_timings_', run_title_id, 
    '.Rdata')
load(previous_R_session_file)

# To reduce the data size, we can throw away all but a fraction of the mcmc
# chains. This has computational (memory) benefits if the MCMC samples are
# strongly autocorrelated, but no other advantages.
mcmc_chain_thin = 20 
```

# **Step 2: Exploratory analysis of seasonal non-stationarity in event statistics**
----------------------------------------------------------------------

**Here we plot the distribution of each storm statistic by month.** This
highlights the seasonal non-stationarity. Below we will take some steps to
check the statistical significance of this, and later will use copula-based
techniques to make the modelled univariate distribution of each variable
conditional on the time of year.

```r
# Get month as 1, 2, ... 12
month_num = as.numeric(format(event_statistics$time, '%m'))
par(mfrow=c(3,2))
for(i in 1:5){
    boxplot(event_statistics[,i] ~ month_num, xlab='Month', 
        ylab=names(event_statistics)[i], names=month.abb,
        col='grey')
    title(main = names(event_statistics)[i], cex.main=2)
}

rm(month_num)
```

![plot of chunk monthplot](figure/monthplot-1.png)

To model the seasonal non-stationarity illustrated above, we define a seasonal
variable periodic in time, of the form `cos(2*pi*(t - offset))` where the time
`t` is in years. The `offset` is a phase variable which can be optimised for
each storm summary statistic separately, to give the 'best' cosine seasonal
pattern matching the data. One way to do this is to find the value of `offset`
which maximises the rank-correlation between each storm variable and the seasonal
variable.

**Below we compute the offset for each storm summary statistic, and also assess
it's statistical significance using a permutation test.** The figure shows the
rank correlation between each variable and a seasonal variable, for each value
of `offset` in [-0.5, 0.5] (which represents all possible values). Note the
`offset` value with the strongest rank correlation may be interpreted as the
optimal offset (*here we choose the `offset` with largest negative rank
correlation, so many `offset`'s are close to zero*). 


```r
# Store some useful statistics
stat_store = data.frame(var = rep(NA, 5), phi=rep(NA,5), cor = rep(NA, 5), 
    p = rep(NA, 5), cor_05=rep(NA, 5))
stat_store$var = names(event_statistics)[1:5]

# Test these values of the 'offset' parameter
phi_vals = seq(-0.5, 0.5, by=0.01)
par(mfrow=c(3,2))
for(i in 1:5){

    # Compute spearman correlation for all values of phi, for variable i
    corrs = phi_vals*0
    for(j in 1:length(phi_vals)){
        corrs[j] =  cor(event_statistics[,i], 
            cos(2*pi*(event_statistics$startyear - phi_vals[j])),
            method='s', use='pairwise.complete.obs')
    }

    plot(phi_vals, corrs, xlab='Offset', ylab='Spearman Correlation', 
        main=names(event_statistics)[i], cex.main=2,
        cex.lab=1.5)
    grid()
    abline(v=0, col='orange')

    # Save the 'best' result
    stat_store$phi[i] = phi_vals[which.min(corrs)]
    stat_store$cor[i] = min(corrs)

    # Function to compute the 'best' correlation of season with
    # permuted data, which by definition has no significant correlation with
    # the season. We can use this to assess the statistical significance of the
    # observed correlation between each variable and the season.
    cor_phi_function<-function(i0=i){
        # Resample the data
        d0 = sample(event_statistics[,i0], size=length(event_statistics[,i0]), 
            replace=TRUE)
        # Correlation function
        g<-function(phi){ 
            cor(d0, cos(2*pi*(event_statistics$startyear - phi)), 
                method='s', use='pairwise.complete.obs')
        }
        # Find best 'phi'
        best_phi = optimize(g, c(-0.5, 0.5), tol=1.0e-06)$minimum

        return(g(best_phi))
    }
   
    # Let's get statistical significance 
    cor_boot = replicate(5000, cor_phi_function())

    # Because our optimizer minimises, the 'strongest' correlations
    # it finds are negative. Of course if 0.5 is added to phi this is equivalent
    # to a positive correlation. 
    
    qcb = quantile(cor_boot, 0.05, type=6)
    stat_store$cor_05[i] = qcb
    stat_store$p[i] = mean(cor_boot < min(corrs))

    polygon(rbind( c(-1, -qcb), c(-1, qcb), c(1, qcb), c(1, -qcb)),
        col='brown', density=10)
}

dir.create('stat_store', showWarnings=FALSE)
write.table(stat_store, 
    file=paste0('stat_store/seasonal_correlation_statistics_', run_title_id, '.csv'), 
    sep="  &  ",
    quote=FALSE, row.names=FALSE)

rm(phi_vals, corrs)
```

![plot of chunk seasonphase1](figure/seasonphase1-1.png)

In the above figure, the shaded region represents a 95% interval for the best
correlation expected of 'random' data (i.e. a random sample of the original
data with an optimized offset).  Correlations outside the shaded interval are
unlikely to occur at random, and are intepreted as reflecting true seasonal
non-stationarity. 

Below we will make each storm summary statistic dependent on the seasonal
variable. For wave direction, the mean annual SOI value will also be treated.
Recall that relationships between mean annual SOI and storm wave direction
were established earlier (
[../preprocessing/extract_storm_events.md](../preprocessing/extract_storm_events.md),
[statistical_model_storm_timings.md](statistical_model_storm_timings.md) ). We
also found relationships between mean annual SOI and the rate of storms, and
MSL, which were treated in those sections (using the non-homogeneous poisson
process model, and the STL decomposition, respectively). Therefore, the latter
relationships are not treated in the section below, but they are included in
the overall model.


# **Step 3: Model the distribution of each storm summary statistic, dependent on season (and mean annual SOI for wave direction)**

In this section we model the distribution of each storm summary statistic, and
then make it conditional on the seasonal variable (and on mean annual SOI in
the case of wave direction only). 

The distributions of `hsig`, `duration` and `tideResid` are initially modelled
as extreme value mixture distributions. The distributions of `dir` and
`steepness` are initially modelled using non-parametric smoothing (based on the
log-spline method).

## Hsig

**Below we fit an extreme value mixture model to Hsig, using maximum
likelihood.** The model has a GPD upper tail, and a Gamma lower tail.

```r
# Get the exmix_fit routines in their own environment
evmix_fit = new.env()
source('../../R/evmix_fit/evmix_fit.R', local=evmix_fit, chdir=TRUE)

# Define the minimum possible value of hsig. The gamma distribution has a lower
# bound of 0, so we need to offset the data to match this. Note that even if hsig
# was perturbed, we account for that in defining hsig_threshold
hsig_offset = hsig_threshold

# Fit it
hsig_mixture_fit = evmix_fit$fit_gpd_mixture(
    data=event_statistics$hsig, 
    data_offset=hsig_offset, 
    bulk='gamma')
```

```
## [1] "  evmix fit NLLH: " "518.713140710555"  
## [1] "Warning: all iterations used"
## [1] "  fit_optim NLLH: " "518.713140596953"  
## [1] "  Bulk par estimate0: " "0.842162501418341"     
## [3] "1.01609622350921"       "1.27053769881576"      
## [5] "-0.214369193836053"    
## [1] "           estimate1: " "0.842159241658271"     
## [3] "1.01609240686169"       "1.27050586311419"      
## [5] "-0.214364006061649"    
## [1] "  Difference: "        "3.25976006920747e-06"  "3.81664751314403e-06" 
## [4] "3.18357015749449e-05"  "-5.18777440447482e-06"
## [1] "PASS: checked qfun and pfun are inverse functions"
```

```r
# Make a plot
DU$qqplot3(event_statistics$hsig, hsig_mixture_fit$qfun(runif(100000)), 
    main='Hsig QQ-plot')
abline(0, 1, col='red'); grid()
```

![plot of chunk hsig_fitA](figure/hsig_fitA-1.png)

The above code leads to print-outs of the maximum likelihood parameter fits
achieved by different methods, and the differences between them (which are only
a few parts per million in this case). Because fitting extreme value mixture
models can be challenging, internally the code tries many different fits.

During the fitting process, we also compute quantile and inverse quantile
functions for the fitted distribution. The code checks numerically that these
really are the inverse of each other, and will print information about whether
this was found to be true (*if not, there is a problem!*)

The quantile-quantile plot of the observed and fitted Hsig should fall close to
a straight line, if the fit worked. Poor fits are suggested by strong
deviations from the 1:1 line. While in this case the fit looks good, if the fit
is poor then further analysis is required. For example, it is possible that the
model fit did not converge, or that the statistical model is a poor choice for
the data.

Given that the above fit looks OK, **below we use Monte-Carlo-Markov-Chain
(MCMC) techniques to compute the Bayesian posterior distribution of the 4 model
parameters**. A few points about this process:
* The prior probability is uniform for each variable. Here we use
a very broad uniform distribution to represent an approximately
'non-informative' prior. The Gamma distribution parameters have uniform prior
over [0, 100 000 000]. The GPD threshold parameter prior is uniform
from zero to the 50th highest data point (to ensure that the tail
part of the model is fit using at least 50 data points). The GPD shape parameter
prior is uniform over [-1000 , 1000]. Note that for some other datasets, it
might be necessary to constrain the GPD shape parameter prior more strongly
than we do below, if it cannot be well estimated from the data (e.g. see the
literature). Overall we are aiming to make our priors reasonably
'non-informative', while still imposing pragmatic constraints required to
achieve a reasonable fit. 
* The routines update the object `hsig_mixture_fit`, so it contains
multiple chains, i.e. 'random walks' through the posterior parameter
distribution.
* Here we run 6 separate chains, with randomly chosen starting parameters, to
make it easier to detect non-convergence (i.e. to reduce the chance that a
single chain gets 'stuck' in part of the posterior distribution). The parameter
`mcmc_start_perturbation` defines the scale for that perturbation.
* It is possible that the randomly chosen start parameters are theoretically
impossible. In this case, the code will report that it had `Bad random start
parameters`, and will generate new ones.
* We use a burn-in of 1000 (i.e. the first 1000 entries in the chain are
discarded). This can assist with convergence.
* We make a simple diagnostic plot to check the MCMC convergence.
* The code runs in parallel, using 6 cores below. The parallel framework will
only work correctly on a shared memory linux machine.

```r
#' MCMC computations for later uncertainty characterisation

# Prevent the threshold parameter from exceeding the highest 50th data point
# Note that inside the fitting routine, Hsig was transformed to have lower
# bound of slightly above zero before fitting, since the Gamma distribution has
# a lower bound of zero. Hence we subtract hsig_offset here.
hsig_u_upper_limit = sort(event_statistics$hsig, decreasing=TRUE)[50] - hsig_offset
hsig_u_lower_limit = hsig_threshold - hsig_offset # Should be 0

# Compute the MCMC chains in parallel
hsig_mixture_fit = evmix_fit$mcmc_gpd_mixture(
    fit_env=hsig_mixture_fit, 
    par_lower_limits=c(0, 0, hsig_u_lower_limit, -1000.), 
    par_upper_limits=c(1e+08, 1.0e+08, hsig_u_upper_limit, 1000),
    mcmc_start_perturbation=c(0.4, 0.4, 2., 0.2), 
    mcmc_length=mcmc_chain_length,
    mcmc_thin=mcmc_chain_thin,
    mcmc_burnin=1000,
    mcmc_nchains=mcmc_nchains,
    mcmc_tune=c(1,1,1,1)*1,
    mc_cores=mcmc_ncores,
    annual_event_rate=mean(events_per_year_truncated))

# Graphical convergence check of one of the chains. 
plot(hsig_mixture_fit$mcmc_chains[[1]])
```

![plot of chunk hsigmixtureFitBayes](figure/hsigmixtureFitBayes-1.png)

**Below, we investigate the parameter estimates for each chain.** If all the
changes have converged, the quantiles of each parameter estimate should be
essentially the same (although if the underlying posterior distribution is
unbounded, then of course the min/max will not converge, although all other
quantiles eventually will). We also look at the 1/100 year event Hsig implied
by each chain, and make a return level plot.

```r
# Look at mcmc parameter estimates in each chain
lapply(hsig_mixture_fit$mcmc_chains, f<-function(x) summary(as.matrix(x)))
```

```
## [[1]]
##       var1             var2             var3              var4        
##  Min.   :0.6465   Min.   :0.7295   Min.   :0.02147   Min.   :-0.4565  
##  1st Qu.:0.8149   1st Qu.:0.9581   1st Qu.:1.07880   1st Qu.:-0.2499  
##  Median :0.8450   Median :1.0090   Median :1.34488   Median :-0.2012  
##  Mean   :0.8451   Mean   :1.0190   Mean   :1.33378   Mean   :-0.1979  
##  3rd Qu.:0.8753   3rd Qu.:1.0674   3rd Qu.:1.62905   3rd Qu.:-0.1486  
##  Max.   :1.0330   Max.   :1.9802   Max.   :2.14746   Max.   : 0.1811  
## 
## [[2]]
##       var1             var2             var3              var4        
##  Min.   :0.6218   Min.   :0.7518   Min.   :0.04965   Min.   :-0.4410  
##  1st Qu.:0.8151   1st Qu.:0.9578   1st Qu.:1.07903   1st Qu.:-0.2505  
##  Median :0.8452   Median :1.0090   Median :1.34879   Median :-0.2013  
##  Mean   :0.8453   Mean   :1.0189   Mean   :1.33636   Mean   :-0.1985  
##  3rd Qu.:0.8751   3rd Qu.:1.0666   3rd Qu.:1.63362   3rd Qu.:-0.1488  
##  Max.   :1.0456   Max.   :1.9460   Max.   :2.14749   Max.   : 0.1651  
## 
## [[3]]
##       var1             var2             var3              var4        
##  Min.   :0.6399   Min.   :0.7589   Min.   :0.01241   Min.   :-0.4417  
##  1st Qu.:0.8147   1st Qu.:0.9582   1st Qu.:1.08136   1st Qu.:-0.2500  
##  Median :0.8449   Median :1.0093   Median :1.34148   Median :-0.2008  
##  Mean   :0.8452   Mean   :1.0188   Mean   :1.33395   Mean   :-0.1980  
##  3rd Qu.:0.8755   3rd Qu.:1.0673   3rd Qu.:1.62674   3rd Qu.:-0.1488  
##  Max.   :1.0499   Max.   :1.9037   Max.   :2.14748   Max.   : 0.1339  
## 
## [[4]]
##       var1             var2             var3              var4        
##  Min.   :0.6102   Min.   :0.7430   Min.   :0.02186   Min.   :-0.4411  
##  1st Qu.:0.8144   1st Qu.:0.9581   1st Qu.:1.07462   1st Qu.:-0.2499  
##  Median :0.8449   Median :1.0092   Median :1.33967   Median :-0.2006  
##  Mean   :0.8448   Mean   :1.0209   Mean   :1.32984   Mean   :-0.1976  
##  3rd Qu.:0.8752   3rd Qu.:1.0674   3rd Qu.:1.62483   3rd Qu.:-0.1474  
##  Max.   :1.0476   Max.   :2.2568   Max.   :2.14748   Max.   : 0.1522  
## 
## [[5]]
##       var1             var2             var3             var4        
##  Min.   :0.6657   Min.   :0.7322   Min.   :0.0388   Min.   :-0.4493  
##  1st Qu.:0.8149   1st Qu.:0.9580   1st Qu.:1.0816   1st Qu.:-0.2504  
##  Median :0.8451   Median :1.0086   Median :1.3471   Median :-0.2018  
##  Mean   :0.8452   Mean   :1.0185   Mean   :1.3368   Mean   :-0.1985  
##  3rd Qu.:0.8754   3rd Qu.:1.0665   3rd Qu.:1.6288   3rd Qu.:-0.1494  
##  Max.   :1.0359   Max.   :1.7276   Max.   :2.1474   Max.   : 0.2883  
## 
## [[6]]
##       var1             var2             var3              var4        
##  Min.   :0.6375   Min.   :0.7173   Min.   :0.02051   Min.   :-0.4495  
##  1st Qu.:0.8152   1st Qu.:0.9570   1st Qu.:1.08814   1st Qu.:-0.2503  
##  Median :0.8453   Median :1.0084   Median :1.35076   Median :-0.2016  
##  Mean   :0.8453   Mean   :1.0189   Mean   :1.34052   Mean   :-0.1985  
##  3rd Qu.:0.8755   3rd Qu.:1.0660   3rd Qu.:1.63500   3rd Qu.:-0.1486  
##  Max.   :1.0706   Max.   :1.8848   Max.   :2.14749   Max.   : 0.2063
```

```r
# Look at ari 100 estimates
lapply(hsig_mixture_fit$ari_100_chains, 
    f<-function(x) quantile(x, p=c(0.025, 0.5, 0.975)))
```

```
## [[1]]
##     2.5%      50%    97.5% 
## 7.071025 7.566212 8.961925 
## 
## [[2]]
##     2.5%      50%    97.5% 
## 7.071145 7.562861 8.962399 
## 
## [[3]]
##     2.5%      50%    97.5% 
## 7.070005 7.565163 8.950241 
## 
## [[4]]
##     2.5%      50%    97.5% 
## 7.074583 7.571037 8.965667 
## 
## [[5]]
##     2.5%      50%    97.5% 
## 7.068235 7.565387 8.954716 
## 
## [[6]]
##     2.5%      50%    97.5% 
## 7.072355 7.565621 8.950204
```

```r
# Look at model prediction of the maximum observed value
# (supposing we observed the system for the same length of time as the data covers)
lapply(hsig_mixture_fit$ari_max_data_chains, 
    f<-function(x) quantile(x, p=c(0.025, 0.5, 0.975)))
```

```
## [[1]]
##     2.5%      50%    97.5% 
## 6.816331 7.204091 8.163641 
## 
## [[2]]
##     2.5%      50%    97.5% 
## 6.817394 7.200616 8.161220 
## 
## [[3]]
##     2.5%      50%    97.5% 
## 6.815820 7.201827 8.151903 
## 
## [[4]]
##     2.5%      50%    97.5% 
## 6.819930 7.205832 8.162963 
## 
## [[5]]
##     2.5%      50%    97.5% 
## 6.815356 7.202598 8.160739 
## 
## [[6]]
##     2.5%      50%    97.5% 
## 6.817280 7.203659 8.156988
```

```r
# If the chains are well behaved, we can combine all 
summary(hsig_mixture_fit$combined_chains)
```

```
##        V1               V2               V3                V4         
##  Min.   :0.6102   Min.   :0.7173   Min.   :0.01241   Min.   :-0.4565  
##  1st Qu.:0.8149   1st Qu.:0.9579   1st Qu.:1.08067   1st Qu.:-0.2502  
##  Median :0.8450   Median :1.0089   Median :1.34555   Median :-0.2012  
##  Mean   :0.8452   Mean   :1.0192   Mean   :1.33521   Mean   :-0.1982  
##  3rd Qu.:0.8753   3rd Qu.:1.0669   3rd Qu.:1.62969   3rd Qu.:-0.1486  
##  Max.   :1.0706   Max.   :2.2568   Max.   :2.14749   Max.   : 0.2883
```

```r
# If the chains are well behaved then we might want a merged 1/100 hsig
quantile(hsig_mixture_fit$combined_ari100, c(0.025, 0.5, 0.975))
```

```
##     2.5%      50%    97.5% 
## 7.071227 7.565993 8.957463
```

```r
# This is an alternative credible interval -- the 'highest posterior density' interval.
HPDinterval(as.mcmc(hsig_mixture_fit$combined_ari100))
```

```
##         lower    upper
## var1 6.974378 8.666495
## attr(,"Probability")
## [1] 0.95
```

```r
evmix_fit$mcmc_rl_plot(hsig_mixture_fit)
```

![plot of chunk hsigmixtureFitBayesB](figure/hsigmixtureFitBayesB-1.png)

**Here we use a different technique to compute the 1/100 AEP Hsig, as a
cross-check on the above analysis.** A simple Generalised Extreme Value model
fit to annual maxima is undertaken. While this technique is based on limited
data (i.e. only one observation per year), it is not dependent on our storm
event definition or choice of wave height threshold. In this sense it is quite
different to our peaks-over-threshold method above -- and thus serves as a
useful cross-check on the former results. 

```r
# Here we do an annual maximum analysis with a gev
# This avoids issues with event definition
annual_max_hsig = aggregate(event_statistics$hsig, 
    list(year=floor(event_statistics$startyear)), max)
# Remove the first and last years with incomplete data
keep_years = which(annual_max_hsig$year %in% 1986:2015)
library(ismev)
```

```
## Loading required package: mgcv
```

```
## Loading required package: nlme
```

```
## This is mgcv 1.8-12. For overview type 'help("mgcv-package")'.
```

```r
gev_fit_annual_max = gev.fit(annual_max_hsig[keep_years,2])
```

```
## $conv
## [1] 0
## 
## $nllh
## [1] 31.07656
## 
## $mle
## [1]  5.4962086  0.6511133 -0.2123074
## 
## $se
## [1] 0.1366988 0.1018133 0.1593190
```

```r
gev.prof(gev_fit_annual_max, m=100, xlow=6.5, xup=12, conf=0.95)
```

```
## If routine fails, try changing plotting interval
```

```r
title(main='Profile likehood confidence interval for 1/100 AEP Hsig \n using a GEV fit to annual maxima')
# Add vertical lines at the limits of the 95% interval
abline(v=c(6.97, 10.4), col='red', lty='dashed')
# Add vertical line at ML estimate
abline(v=7.4, col='orange')
```

![plot of chunk gevHsigFit](figure/gevHsigFit-1.png)

**Here we use copulas to determine a distribution for Hsig, conditional on the season**.
The computational details are wrapped up in a function that we source.
Essentially, the code:
* Finds the optimal seasonal `offset` for the chosen variable (`hsig`), and uses
this to create a function to compute the season statistic (which is `hsig`
specific) from the event time of year.
* Automatically chooses a copula family (based on AIC) to model dependence
between the chosen variable and the season variable, and fits the copula.
* Uses the copula to create new quantile and inverse quantile functions, for
which the user can pass conditional variables (i.e. to get the distribution,
given that the season variable attains a particular value).
* Test that the quantile and inverse quantile functions really are inverses of
each other (this can help catch user input errors)
* Make quantile-quantile plots of the data and model for a number of time
periods (here the first, middle and last thirds of the calendar year). The top
row shows the model with the distribution varying by season, and the bottom row
shows the model without seasonal effects. It is not particularly easy to
visually detect seasonal non-stationarities in these plots [compared, say, with
using monthly boxplots].  Their main purpose is compare the model and data
distribution at different times of year, and so detect poor model fits.
However, you might notice that the top row of plots 'hug' the 1:1 line slightly
better than corresponding plots in the bottom row in the central data values.
This reflects the modelled seasonal non-stationarities. *Note the tail behaviour
can be erratic, since the 'model' result is actually a random sample from the model.*

```r
# Get code to fit the conditional distribution
# Give a path that will also work if run from another directory inside Analysis.
source('../../Analysis/statistical_model_fit/make_conditional_distribution.R', chdir=TRUE)

# This returns an environment containing the conditional quantile and inverse
# quantile functions, among other information
hsig_fit_conditional = make_fit_conditional_on_season(
    event_statistics,
    var='hsig', 
    q_raw=hsig_mixture_fit$qfun, 
    p_raw=hsig_mixture_fit$pfun,
    startyear = 'startyear')
```

```
## [1] "Conditional p/q functions passed test: "
## [1] "  (Check plots to see if quantiles are ok)"
```

![plot of chunk fitCopulaHsig](figure/fitCopulaHsig-1.png)

```r
# What kind of copula was selected to model dependence between season and hsig?
print(hsig_fit_conditional$var_season_copula)
```

```
## Bivariate copula: Gaussian (par = -0.13, tau = -0.08)
```


## Duration

Here we model storm `duration`, using techniques virtually identical to those applied above.
As before:
* We first fit the univariate extreme value mixture distribution with maximum
likelihood; 
* Next we compute the posterior distribution of each parameter; 
* Finally we make the `duration` distribution conditional on the time of year,
using a seasonal variable that has been optimised to capture seasonality in the
storm `duration`.

**Here is the extreme value mixture model maximum likelihood fit**

```r
# Do the maximum likelihood fit. 
#
# Set the lower limit of the tail model to just below the lower limit of the
# data, in the event we perturb it by half an hour. [Initial lower limit = 1hr]
# However, if working with the 'raw' data, we should not set the lower limit to 1hr,
# since the latter is actually a data value, so some lower bound must be selected.
duration_offset = ifelse(break_ties_with_jitter, 
    1 - as.numeric(default_jitter_amounts['duration']), 
    0.5) 

duration_mixture_fit = evmix_fit$fit_gpd_mixture(
    data=event_statistics$duration, 
    data_offset=duration_offset, 
    bulk='gamma')
```

```
## Warning in FUN(X[[i]], ...): initial parameter values for threshold u = 0.5
## are invalid

## Warning in FUN(X[[i]], ...): initial parameter values for threshold u = 0.5
## are invalid
```

```
## [1] "  evmix fit NLLH: " "2752.84379433546"  
## [1] "  fit_optim NLLH: " "2752.84379427367"  
## [1] "  Bulk par estimate0: " "0.693823596763882"     
## [3] "35.9454038109696"       "48.0295299401044"      
## [5] "-0.163363982336055"    
## [1] "           estimate1: " "0.693829686546937"     
## [3] "35.9454043148681"       "48.0295307213169"      
## [5] "-0.163372906989212"    
## [1] "  Difference: "        "-6.08978305538521e-06" "-5.03898476722497e-07"
## [4] "-7.81212513345508e-07" "8.92465315702196e-06" 
## [1] "PASS: checked qfun and pfun are inverse functions"
```

```r
# Make a plot
DU$qqplot3(event_statistics$duration, duration_mixture_fit$qfun(runif(100000)), 
    main='Duration QQ-plot')
abline(0, 1, col='red'); grid()
```

![plot of chunk durationMixtureML](figure/durationMixtureML-1.png)

**Here is the extreme value mixture model posterior probability computation, using MCMC**
As before, note that we run a number of MCMC chains with random starting values, and in 
the event that the random starting parameters are invalid the code will simply try new ones.

```r
# MCMC computations for later uncertainty characterisation
#
# Prevent the threshold parameter from exceeding the highest 50th data point
duration_u_upper_limit = sort(event_statistics$duration, decreasing=TRUE)[50] - duration_offset
duration_u_lower_limit = 0

# Compute the MCMC chains in parallel.
duration_mixture_fit = evmix_fit$mcmc_gpd_mixture(
    fit_env=duration_mixture_fit, 
    par_lower_limits=c(0, 0, duration_u_lower_limit, -1000.), 
    par_upper_limits=c(1e+08, 1.0e+08, duration_u_upper_limit, 1000),
    mcmc_start_perturbation=c(0.4, 0.4, 2., 0.2), 
    mcmc_length=mcmc_chain_length,
    mcmc_thin=mcmc_chain_thin,
    mcmc_burnin=1000,
    mcmc_nchains=mcmc_nchains,
    mcmc_tune=c(1,1,1,1)*1,
    mc_cores=mcmc_ncores,
    annual_event_rate=mean(events_per_year_truncated))

# Graphical convergence check of one of the chains. 
plot(duration_mixture_fit$mcmc_chains[[1]])
```

![plot of chunk durationmixtureFitBayes](figure/durationmixtureFitBayes-1.png)

**Here we check the similarity of all the MCMC chains, and make a return-level
plot for storm `duration`**

```r
# Look at mcmc parameter estimates in each chain
lapply(duration_mixture_fit$mcmc_chains, f<-function(x) summary(as.matrix(x)))
```

```
## [[1]]
##       var1             var2            var3             var4         
##  Min.   :0.5260   Min.   :26.63   Min.   : 3.487   Min.   :-0.34612  
##  1st Qu.:0.6648   1st Qu.:34.53   1st Qu.:30.829   1st Qu.:-0.17470  
##  Median :0.6892   Median :36.45   Median :45.256   Median :-0.12025  
##  Mean   :0.6890   Mean   :37.07   Mean   :42.729   Mean   :-0.11564  
##  3rd Qu.:0.7137   3rd Qu.:38.81   3rd Qu.:55.835   3rd Qu.:-0.06102  
##  Max.   :0.8458   Max.   :72.17   Max.   :70.498   Max.   : 0.35519  
## 
## [[2]]
##       var1             var2            var3             var4         
##  Min.   :0.5179   Min.   :26.19   Min.   : 3.602   Min.   :-0.34868  
##  1st Qu.:0.6650   1st Qu.:34.50   1st Qu.:30.608   1st Qu.:-0.17449  
##  Median :0.6896   Median :36.45   Median :45.221   Median :-0.11960  
##  Mean   :0.6893   Mean   :37.02   Mean   :42.705   Mean   :-0.11511  
##  3rd Qu.:0.7139   3rd Qu.:38.77   3rd Qu.:55.969   3rd Qu.:-0.06028  
##  Max.   :0.8385   Max.   :69.04   Max.   :70.499   Max.   : 0.33582  
## 
## [[3]]
##       var1             var2            var3             var4         
##  Min.   :0.5169   Min.   :27.79   Min.   : 3.582   Min.   :-0.34286  
##  1st Qu.:0.6651   1st Qu.:34.51   1st Qu.:30.994   1st Qu.:-0.17456  
##  Median :0.6898   Median :36.43   Median :45.162   Median :-0.12010  
##  Mean   :0.6894   Mean   :37.01   Mean   :42.771   Mean   :-0.11555  
##  3rd Qu.:0.7139   3rd Qu.:38.74   3rd Qu.:55.798   3rd Qu.:-0.06074  
##  Max.   :0.8328   Max.   :65.35   Max.   :70.499   Max.   : 0.33824  
## 
## [[4]]
##       var1             var2            var3             var4         
##  Min.   :0.5181   Min.   :26.60   Min.   : 3.297   Min.   :-0.35807  
##  1st Qu.:0.6649   1st Qu.:34.52   1st Qu.:30.746   1st Qu.:-0.17484  
##  Median :0.6895   Median :36.46   Median :45.197   Median :-0.11981  
##  Mean   :0.6890   Mean   :37.09   Mean   :42.729   Mean   :-0.11576  
##  3rd Qu.:0.7137   3rd Qu.:38.78   3rd Qu.:55.916   3rd Qu.:-0.06072  
##  Max.   :0.8352   Max.   :74.75   Max.   :70.499   Max.   : 0.34363  
## 
## [[5]]
##       var1             var2            var3             var4         
##  Min.   :0.5321   Min.   :26.90   Min.   : 3.259   Min.   :-0.34948  
##  1st Qu.:0.6649   1st Qu.:34.51   1st Qu.:30.842   1st Qu.:-0.17461  
##  Median :0.6896   Median :36.44   Median :45.146   Median :-0.12014  
##  Mean   :0.6893   Mean   :37.03   Mean   :42.722   Mean   :-0.11574  
##  3rd Qu.:0.7142   3rd Qu.:38.79   3rd Qu.:55.889   3rd Qu.:-0.06139  
##  Max.   :0.8506   Max.   :66.95   Max.   :70.498   Max.   : 0.26321  
## 
## [[6]]
##       var1             var2            var3             var4         
##  Min.   :0.5219   Min.   :27.16   Min.   : 3.821   Min.   :-0.34901  
##  1st Qu.:0.6652   1st Qu.:34.52   1st Qu.:30.709   1st Qu.:-0.17411  
##  Median :0.6895   Median :36.43   Median :45.206   Median :-0.11960  
##  Mean   :0.6894   Mean   :37.01   Mean   :42.745   Mean   :-0.11543  
##  3rd Qu.:0.7139   3rd Qu.:38.76   3rd Qu.:55.910   3rd Qu.:-0.06089  
##  Max.   :0.8387   Max.   :64.35   Max.   :70.499   Max.   : 0.31240
```

```r
# Look at ari 100 estimates
lapply(duration_mixture_fit$ari_100_chains, 
    f<-function(x) quantile(x, p=c(0.025, 0.5, 0.975)))
```

```
## [[1]]
##     2.5%      50%    97.5% 
## 151.0460 175.8854 247.5436 
## 
## [[2]]
##     2.5%      50%    97.5% 
## 151.1405 175.9947 248.1162 
## 
## [[3]]
##     2.5%      50%    97.5% 
## 151.1846 175.8902 246.6152 
## 
## [[4]]
##     2.5%      50%    97.5% 
## 151.0658 175.8605 246.4530 
## 
## [[5]]
##     2.5%      50%    97.5% 
## 151.1150 175.8877 247.3284 
## 
## [[6]]
##     2.5%      50%    97.5% 
## 151.1149 176.0287 247.0162
```

```r
# Look at model prediction of the maximum observed value
# (supposing we observed the system for the same length of time as the data covers)
lapply(duration_mixture_fit$ari_max_data_chains, 
    f<-function(x) quantile(x, p=c(0.025, 0.5, 0.975)))
```

```
## [[1]]
##     2.5%      50%    97.5% 
## 138.6896 156.8650 202.8168 
## 
## [[2]]
##     2.5%      50%    97.5% 
## 138.7187 156.9115 203.0560 
## 
## [[3]]
##     2.5%      50%    97.5% 
## 138.7641 156.8510 202.2531 
## 
## [[4]]
##     2.5%      50%    97.5% 
## 138.7841 156.9096 202.1857 
## 
## [[5]]
##     2.5%      50%    97.5% 
## 138.6258 156.8185 202.4731 
## 
## [[6]]
##     2.5%      50%    97.5% 
## 138.7873 156.9765 202.4894
```

```r
# If the chains seem ok, we can combine all 
summary(duration_mixture_fit$combined_chains)
```

```
##        V1               V2              V3               V4          
##  Min.   :0.5169   Min.   :26.19   Min.   : 3.259   Min.   :-0.35807  
##  1st Qu.:0.6650   1st Qu.:34.52   1st Qu.:30.785   1st Qu.:-0.17456  
##  Median :0.6895   Median :36.44   Median :45.197   Median :-0.11992  
##  Mean   :0.6892   Mean   :37.04   Mean   :42.733   Mean   :-0.11554  
##  3rd Qu.:0.7139   3rd Qu.:38.78   3rd Qu.:55.886   3rd Qu.:-0.06085  
##  Max.   :0.8506   Max.   :74.75   Max.   :70.499   Max.   : 0.35519
```

```r
# If the chains are well behaved then we might want a merged 1/100 duration
quantile(duration_mixture_fit$combined_ari100, c(0.025, 0.5, 0.975))
```

```
##     2.5%      50%    97.5% 
## 151.1143 175.9267 247.1504
```

```r
HPDinterval(as.mcmc(duration_mixture_fit$combined_ari100))
```

```
##         lower    upper
## var1 146.5178 232.3902
## attr(,"Probability")
## [1] 0.95
```

```r
# Return level plot
evmix_fit$mcmc_rl_plot(duration_mixture_fit)
```

![plot of chunk durationMCMCcheck](figure/durationMCMCcheck-1.png)

**Finally we make the `duration` fit conditional on the time of year, using a
seasonal variable.** The seasonal QQ-plots below highlight the hourly
discretization in the `duration` data, which is prominent at low quantiles. It
is important to check whether this effects the analysis through randomization
of the data, which we undertake in the section
[../statistical_model_fit_pertubed_data](../statistical_model_fit_perturbed_data).

```r
# This returns an environment containing the conditional quantile and inverse
# quantile functions, among other information
duration_fit_conditional = make_fit_conditional_on_season(
    event_statistics,
    var='duration', 
    q_raw=duration_mixture_fit$qfun, 
    p_raw=duration_mixture_fit$pfun,
    startyear = 'startyear')
```

```
## [1] "Conditional p/q functions passed test: "
## [1] "  (Check plots to see if quantiles are ok)"
```

![plot of chunk fitCopulaDuration](figure/fitCopulaDuration-1.png)

```r
# What kind of copula was selected to model dependence between season and duration?
print(duration_fit_conditional$var_season_copula)
```

```
## Bivariate copula: Gaussian (par = -0.15, tau = -0.1)
```

## Tidal residual

Here we generally follow the steps implemented above for `hsig` and `duration`,
for modelling the tidal residual `tideResid`. An important change is that we
fit an extreme value mixture model with a normal lower tail (instead of a Gamma
lower tail). This is done because unlike storm `hsig` and `duration`, there is
no natural lower limit on the `tideResid` [e.g.  it can even be negative on
occasion].


```r
# Manually remove missing (NA) data before fitting
tideResid_mixture_fit = evmix_fit$fit_gpd_mixture(
    data=na.omit(event_statistics$tideResid),  
    bulk='normal',
    gpd_threshold_quantile_range=c(0.05, 0.85) # Numerical reasons
    )
```

```
## [1] "  evmix fit NLLH: " "-433.252459580517" 
## [1] "  fit_optim NLLH: " "-433.252460377451" 
## [1] "  Bulk par estimate0: " "0.115246548723832"     
## [3] "0.11463701571688"       "0.197500638710626"     
## [5] "-0.112013335868541"    
## [1] "           estimate1: " "0.115246500854982"     
## [3] "0.114636889799082"      "0.197499970552887"     
## [5] "-0.112079837783983"    
## [1] "  Difference: "       "4.7868849911703e-08"  "1.25917797341724e-07"
## [4] "6.68157738387132e-07" "6.65019154414553e-05"
## [1] "PASS: checked qfun and pfun are inverse functions"
```

```r
# Make a plot
DU$qqplot3(na.omit(event_statistics$tideResid), 
    tideResid_mixture_fit$qfun(runif(100000)))
abline(0, 1, col='red')
grid()
```

![plot of chunk tideResidExtremeValueMixture](figure/tideResidExtremeValueMixture-1.png)

Below is the MCMC computation of the posterior probability distribution.
As before, bad random starting parameters are rejected, with a warning.

```r
#' MCMC computations for later uncertainty characterisation
min_tr = min(event_statistics$tideResid, na.rm=TRUE)

tideResid_u_upper_limit = sort(event_statistics$tideResid, decreasing=TRUE)[50]
tideResid_u_lower_limit = min_tr

tideResid_mixture_fit = evmix_fit$mcmc_gpd_mixture(
    fit_env=tideResid_mixture_fit, 
    par_lower_limits=c(min_tr, 0, tideResid_u_lower_limit, -1000), 
    par_upper_limits=c(1e+08, 1e+08, tideResid_u_upper_limit, 1000),
    mcmc_start_perturbation=c(0.2, 0.2, 0.2, 0.3), 
    mcmc_length=mcmc_chain_length,
    mcmc_thin=mcmc_chain_thin,
    mcmc_burnin=1000,
    mcmc_nchains=mcmc_nchains,
    mcmc_tune=c(1,1,1,1)*1.,
    mc_cores=mcmc_ncores,
    annual_event_rate=mean(events_per_year_truncated))

# Graphical convergence check
plot(tideResid_mixture_fit$mcmc_chains[[1]])
```

![plot of chunk tideResidMCMC](figure/tideResidMCMC-1.png)

```r
# Clean up
rm(min_tr)
```

**Here we further investigate the behaviour of the MCMC chains for the tidal
residual fit, and make a return-level plot**

```r
# Look at mcmc parameter estimates in each chain
lapply(tideResid_mixture_fit$mcmc_chains, f<-function(x) summary(as.matrix(x)))
```

```
## [[1]]
##       var1              var2              var3             var4         
##  Min.   :0.09322   Min.   :0.09928   Min.   :0.1156   Min.   :-0.25155  
##  1st Qu.:0.11271   1st Qu.:0.11359   1st Qu.:0.2040   1st Qu.:-0.10523  
##  Median :0.11597   Median :0.11618   Median :0.2336   Median :-0.05436  
##  Mean   :0.11596   Mean   :0.11618   Mean   :0.2329   Mean   :-0.03773  
##  3rd Qu.:0.11921   3rd Qu.:0.11880   3rd Qu.:0.2654   3rd Qu.: 0.01083  
##  Max.   :0.13474   Max.   :0.13180   Max.   :0.2917   Max.   : 0.63688  
## 
## [[2]]
##       var1              var2             var3             var4         
##  Min.   :0.09453   Min.   :0.1004   Min.   :0.1120   Min.   :-0.24881  
##  1st Qu.:0.11271   1st Qu.:0.1136   1st Qu.:0.2041   1st Qu.:-0.10474  
##  Median :0.11593   Median :0.1162   Median :0.2340   Median :-0.05455  
##  Mean   :0.11595   Mean   :0.1162   Mean   :0.2332   Mean   :-0.03717  
##  3rd Qu.:0.11919   3rd Qu.:0.1188   3rd Qu.:0.2658   3rd Qu.: 0.01150  
##  Max.   :0.13501   Max.   :0.1332   Max.   :0.2917   Max.   : 0.67773  
## 
## [[3]]
##       var1              var2              var3             var4         
##  Min.   :0.09669   Min.   :0.09922   Min.   :0.1111   Min.   :-0.24503  
##  1st Qu.:0.11266   1st Qu.:0.11353   1st Qu.:0.2038   1st Qu.:-0.10456  
##  Median :0.11594   Median :0.11616   Median :0.2337   Median :-0.05432  
##  Mean   :0.11594   Mean   :0.11617   Mean   :0.2331   Mean   :-0.03765  
##  3rd Qu.:0.11918   3rd Qu.:0.11878   3rd Qu.:0.2655   3rd Qu.: 0.01130  
##  Max.   :0.13541   Max.   :0.13180   Max.   :0.2917   Max.   : 0.67869  
## 
## [[4]]
##       var1              var2             var3             var4         
##  Min.   :0.09617   Min.   :0.1003   Min.   :0.1092   Min.   :-0.24683  
##  1st Qu.:0.11266   1st Qu.:0.1136   1st Qu.:0.2038   1st Qu.:-0.10487  
##  Median :0.11593   Median :0.1162   Median :0.2341   Median :-0.05470  
##  Mean   :0.11593   Mean   :0.1162   Mean   :0.2332   Mean   :-0.03775  
##  3rd Qu.:0.11921   3rd Qu.:0.1188   3rd Qu.:0.2658   3rd Qu.: 0.01115  
##  Max.   :0.13689   Max.   :0.1320   Max.   :0.2917   Max.   : 0.87237  
## 
## [[5]]
##       var1              var2              var3             var4         
##  Min.   :0.09541   Min.   :0.09891   Min.   :0.1113   Min.   :-0.24535  
##  1st Qu.:0.11266   1st Qu.:0.11355   1st Qu.:0.2037   1st Qu.:-0.10565  
##  Median :0.11593   Median :0.11618   Median :0.2335   Median :-0.05432  
##  Mean   :0.11594   Mean   :0.11618   Mean   :0.2328   Mean   :-0.03814  
##  3rd Qu.:0.11923   3rd Qu.:0.11879   3rd Qu.:0.2654   3rd Qu.: 0.01136  
##  Max.   :0.13679   Max.   :0.13244   Max.   :0.2917   Max.   : 0.61219  
## 
## [[6]]
##       var1              var2             var3             var4         
##  Min.   :0.09632   Min.   :0.1010   Min.   :0.1169   Min.   :-0.23047  
##  1st Qu.:0.11272   1st Qu.:0.1136   1st Qu.:0.2040   1st Qu.:-0.10521  
##  Median :0.11597   Median :0.1162   Median :0.2339   Median :-0.05484  
##  Mean   :0.11596   Mean   :0.1162   Mean   :0.2332   Mean   :-0.03731  
##  3rd Qu.:0.11919   3rd Qu.:0.1188   3rd Qu.:0.2657   3rd Qu.: 0.01180  
##  Max.   :0.13783   Max.   :0.1321   Max.   :0.2917   Max.   : 0.66314
```

```r
# Look at ari 100 estimates
lapply(tideResid_mixture_fit$ari_100_chains, 
    f<-function(x) quantile(x, p=c(0.025, 0.5, 0.975)))
```

```
## [[1]]
##      2.5%       50%     97.5% 
## 0.5374044 0.6168566 0.8948476 
## 
## [[2]]
##      2.5%       50%     97.5% 
## 0.5371809 0.6164556 0.9013884 
## 
## [[3]]
##      2.5%       50%     97.5% 
## 0.5366296 0.6167381 0.8920260 
## 
## [[4]]
##      2.5%       50%     97.5% 
## 0.5375629 0.6163248 0.8898772 
## 
## [[5]]
##      2.5%       50%     97.5% 
## 0.5367417 0.6166927 0.8905517 
## 
## [[6]]
##      2.5%       50%     97.5% 
## 0.5370156 0.6164471 0.8981844
```

```r
# Look at model prediction of the maximum observed value
# (supposing we observed the system for the same length of time as the data covers)
lapply(tideResid_mixture_fit$ari_max_data_chains, 
    f<-function(x) quantile(x, p=c(0.025, 0.5, 0.975)))
```

```
## [[1]]
##      2.5%       50%     97.5% 
## 0.4876126 0.5438992 0.6909386 
## 
## [[2]]
##      2.5%       50%     97.5% 
## 0.4874693 0.5435498 0.6954991 
## 
## [[3]]
##      2.5%       50%     97.5% 
## 0.4871322 0.5436727 0.6926639 
## 
## [[4]]
##      2.5%       50%     97.5% 
## 0.4877267 0.5436781 0.6906689 
## 
## [[5]]
##      2.5%       50%     97.5% 
## 0.4873131 0.5439120 0.6912303 
## 
## [[6]]
##      2.5%       50%     97.5% 
## 0.4873816 0.5437372 0.6921664
```

```r
# If the chains seem ok, we can combine all 
summary(tideResid_mixture_fit$combined_chains)
```

```
##        V1                V2                V3               V4          
##  Min.   :0.09322   Min.   :0.09891   Min.   :0.1092   Min.   :-0.25155  
##  1st Qu.:0.11268   1st Qu.:0.11357   1st Qu.:0.2039   1st Qu.:-0.10504  
##  Median :0.11594   Median :0.11618   Median :0.2338   Median :-0.05452  
##  Mean   :0.11595   Mean   :0.11618   Mean   :0.2331   Mean   :-0.03763  
##  3rd Qu.:0.11920   3rd Qu.:0.11880   3rd Qu.:0.2656   3rd Qu.: 0.01132  
##  Max.   :0.13783   Max.   :0.13323   Max.   :0.2917   Max.   : 0.87237
```

```r
# If the chains are well behaved then we might want a merged 1/100 hsig
quantile(tideResid_mixture_fit$combined_ari100, c(0.025, 0.5, 0.975))
```

```
##      2.5%       50%     97.5% 
## 0.5370639 0.6165770 0.8946910
```

```r
HPDinterval(as.mcmc(tideResid_mixture_fit$combined_ari100))
```

```
##          lower     upper
## var1 0.5195418 0.8237584
## attr(,"Probability")
## [1] 0.95
```

```r
# Return level plot
evmix_fit$mcmc_rl_plot(tideResid_mixture_fit)
```

![plot of chunk tideResidMCMCcheck](figure/tideResidMCMCcheck-1.png)

Below we make the tidal residual distribution conditional on the time of year,
via a seasonal variable with phase optimized to model tidal residual seasonality.

```r
# This returns an environment containing the conditional quantile and inverse
# quantile functions, among other information
tideResid_fit_conditional = make_fit_conditional_on_season(
    event_statistics,
    var='tideResid', 
    q_raw=tideResid_mixture_fit$qfun, 
    p_raw=tideResid_mixture_fit$pfun,
    startyear = 'startyear')
```

```
## [1] "Conditional p/q functions passed test: "
## [1] "  (Check plots to see if quantiles are ok)"
```

![plot of chunk tideResidDependence](figure/tideResidDependence-1.png)

```r
# What kind of copula was selected to model dependence between season and
# tidal residual?
print(tideResid_fit_conditional$var_season_copula)
```

```
## Bivariate copula: Gaussian (par = -0.16, tau = -0.1)
```

## Steepness

**Here we model the distribution of wave `steepness` conditional on the time of
year.** The approach differs from that used above in that the density is
modelled with non-parametric techniques. This was done because for wave
`steepness` we are not so interested in extrapolating beyond the range of the
data, whereas such extrapolation is necessary for `hsig`, `duration`, and
`tideResid`. 

**The non-parametric fit of wave steepness is implemented here**

```r
library(logspline)

# Use the 'old' logspline density estimator
qsteepness0 = logspline::oldlogspline(
    na.omit(event_statistics$steepness), 
    lbound=min(event_statistics$steepness, na.rm=TRUE)-1.0e-04, 
    ubound=max(event_statistics$steepness, na.rm=TRUE)+1.0e-03)

# Make 'raw' quantile and inverse quantile functions
# With linear approximation we can be sure they are inverses of each other
ptmp = seq(0, 1, len=1000)
qtmp = logspline::qoldlogspline(ptmp, qsteepness0)
# Quantile function
qsteepness_raw = approxfun(ptmp, qtmp)
# Inverse quantile function
psteepness_raw = approxfun(qtmp, ptmp)  

rm(ptmp, qtmp, qsteepness0)

# Plot it
par(mfrow=c(1,2))
plot(ecdf(event_statistics$steepness), 
    main = 'Wave steepness empirical cumulative distribution',
    xlab='Steepness')
points(qsteepness_raw(seq(0,1,len=200)), seq(0,1,len=200), t='l', 
    col='red', lwd=2, lty='dashed')
legend('topleft', c('Data', 'Model'), col=c('black', 'red'), 
    lty=c('solid', 'dashed'), pch=c(19,NA), lwd=c(1,2))

x = seq(0,1,len=200)
hist(na.omit(event_statistics$steepness), freq=FALSE, breaks=30,
    main = 'Wave steepness density', xlab='Steepness')
points(qsteepness_raw(x),
    c(0, diff(x)/diff(qsteepness_raw(x))), t='l', col='red')
legend('topright', c('Model'), lty=1, col='red')
```

![plot of chunk steepnessdensity](figure/steepnessdensity-1.png)

**The conditional modelling is implemented here**. It follows the same
approach as used above.

```r
# This returns an environment containing the conditional quantile and inverse
# quantile functions, among other information
steepness_fit_conditional = make_fit_conditional_on_season(
    event_statistics,
    var='steepness', 
    q_raw=qsteepness_raw, 
    p_raw=psteepness_raw,
    startyear = 'startyear')
```

```
## [1] "Conditional p/q functions passed test: "
## [1] "  (Check plots to see if quantiles are ok)"
```

![plot of chunk steepnessConditional](figure/steepnessConditional-1.png)

```r
# What kind of copula was selected to model dependence between season and
# steepness?
print(steepness_fit_conditional$var_season_copula)
```

```
## Bivariate copula: Rotated Clayton 270 degrees (par = -0.19, tau = -0.09)
```

## Wave direction

Here we model the distribution of wave direction `dir`, dependent on both the
season, and the mean annual SOI. The latter was not treated above. 

**Below the density of `dir` is estimated with a logspline smoother**

```r
# Fit logspline to direction distribution
# Deliberately extend the range slightly beyond the data range
qdir0 = logspline::oldlogspline(
    na.omit(event_statistics$dir), 
    lbound=min(event_statistics$dir, na.rm=TRUE)-0.5, 
    ubound=max(event_statistics$dir, na.rm=TRUE)+0.5)

# Make 'raw' quantile and inverse quantile functions
# With linear approximation we can be sure they are inverses of each other
ptmp = seq(0, 1, len=1000)
qtmp = logspline::qoldlogspline(seq(0, 1, len=1000), qdir0)
# Quantile
qdir_raw = approxfun(ptmp, qtmp)
# Inverse quantile
pdir_raw = approxfun(qtmp, ptmp)  
# Cleanup
rm(ptmp, qtmp, qdir0)

# Plot it
par(mfrow=c(1,2))
plot(ecdf(event_statistics$dir), 
    main = 'Wave direction empirical cumulative distribution',
    xlab='Direction (degrees)')
points(qdir_raw(seq(0,1,len=200)), seq(0,1,len=200), t='l', 
    col='red', lwd=2, lty='dashed')
legend('topleft', c('Data', 'Model'), col=c('black', 'red'), 
    lty=c('solid', 'dashed'), pch=c(19,NA), lwd=c(1,2))

x = seq(0,1,len=200)
hist(na.omit(event_statistics$dir), freq=FALSE, breaks=30,
    main = 'Wave direction density', xlab='Direction (degrees)')
points(qdir_raw(x),
    c(0, diff(x)/diff(qdir_raw(x))), t='l', col='red')
legend('topleft', c('Model'), lty=1, col='red')
```

![plot of chunk directionDensity](figure/directionDensity-1.png)

**Below we derive the distribution of wave direction, conditional
on the mean annual SOI, and season**


```r
dir_fit_conditional = make_fit_conditional_on_soiA_and_season(
    event_statistics,
    var='dir', 
    q_raw=qdir_raw, 
    p_raw=pdir_raw,
    startyear = 'startyear',
    soiA = 'soiA')
```

```
## [1] "Conditional p/q functions passed test: Check quantiles are ok"
```

![plot of chunk directionConditional](figure/directionConditional-1.png)

```r
# What kind of copula was selected to model dependence between soiA and
# direction?
print(dir_fit_conditional$var_soiA_copula)
```

```
## Bivariate copula: Frank (par = -0.81, tau = -0.09)
```

```r
# What kind of copula was selected to model dependence between season and
# 'direction given soiA'?
print(dir_fit_conditional$vargivensoiA_season_copula)
```

```
## Bivariate copula: Frank (par = -0.74, tau = -0.08)
```

## Collating all the univariate distributions

For convenience later on, here we make a function to help generate random storm
properties conditional on the time of year and soiA. The function takes as
input a set of parameter percentiles in (0-1) (one for each storm statistic),
and a set of known values for the conditional variables (time of year and
soiA). The function uses this to generate values of the storm summary
statistics, conditional on the provided values of the conditional variables
(time of year and soiA). If the parameter percentiles are drawn from a uniform
distribution, then the random storm properties will correspond to the
distributions fit above.

```r
make_stormVarFun<-function(qduration=NULL, qhsig=NULL, qtideResid=NULL, qdir=NULL,
    qsteepness=NULL){

    qduration = qduration
    qhsig = qhsig
    qtideResid = qtideResid
    qdir = qdir
    qsteepness = qsteepness

    #' 
    #' Compute storm variables from inputs vector, which are all of the same length and
    #' in (0-1).
    #'
    #' The input vectors give the percentiles of values in the distribution of each
    #' variable.
    #' 
    #' @param duration The duration percentile, in (0-1)
    #' @param hsig The hsig percentile, in (0-1)
    #' @param steepness the steepness percentile, in (0,1)
    #' @param tideResid the tideResid percentile, in (0,1)
    #' @param conditional_variables list containing other variables required to evaluate
    #' the quantile functions (for example, complex models might use the time of
    #' year, or the soiA value)
    #'
    stormVarFun<-function(duration, hsig, dir, steepness, tideResid, 
        conditional_variables=NULL){

        # duration
        duration_vals = qduration(duration, conditional_variables)

        # hsig
        hsig_vals = qhsig(hsig, conditional_variables)

        # tidal residual
        tr_vals = qtideResid(tideResid, conditional_variables)

        # direction 
        dir_vals = qdir(dir, conditional_variables)

        # steepness
        steepness_vals = qsteepness(steepness, conditional_variables)
        stopifnot(all(steepness_vals > 0.))

        return(data.frame(duration=duration_vals, hsig=hsig_vals, 
            tideResid=tr_vals, dir=dir_vals, steepness = steepness_vals))
    }

    return(stormVarFun)
}

#' Get a function which transforms vectors in [0,1] to model quantiles
stormVarFun = make_stormVarFun(
    qduration = duration_fit_conditional$qfun,
    qhsig = hsig_fit_conditional$qfun,
    qtideResid = tideResid_fit_conditional$qfun,
    qdir = dir_fit_conditional$qfun,
    qsteepness = steepness_fit_conditional$qfun)
```


**Save output for use later.** 

We use the same `run_title_id` as was computed in the previous section
([statistical_model_storm_timings.md](statistical_model_storm_timings.md)).

```r
dir.create('Rimages', showWarnings=FALSE)
Rimage_title = paste0('Rimages/session_univariate_distributions_', run_title_id, '.Rdata')
save.image(Rimage_title)
```


## **Moving On**
The next part of this vignette begins at
[statistical_model_vine_copula.md](statistical_model_vine_copula.md).
