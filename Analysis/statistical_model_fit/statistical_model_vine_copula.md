
# **Statistical modelling of multivariate dependence in the storm event dataset**
------------------------------------------------------------------------------------

*Gareth Davies, Geoscience Australia 2017*

# Introduction
------------------

This document follows on from
[statistical_model_univariate_distributions.md](statistical_model_univariate_distributions.md)
in describing our statistical analysis of storm waves at Old Bar. 

It illustrates the process of modelling the joint distribution of all storm
variables, after their univariate distributions have already been fit. Recall
from the previous section that the univariate distributions are themselves
conditional on the event time of year, and mean annual SOI value. 

It is essential that the code in
[statistical_model_univariate_distributions.md](statistical_model_univariate_distributions.md)
has alread been run, and produced an Rdata file
*'Rimages/session_univariate_distributions_FALSE_0.Rdata'*. **To make sure, the
code below throws an error if the latter file does not exist.**

```r
# Check that the pre-requisites exist
if(!file.exists('Rimages/session_univariate_distributions_FALSE_0.Rdata')){
    stop('It appears you have not yet run the code in statistical_model_univariate_distributions.Rmd. That must be run before continuing')
}
```

Supposing the above did not generate any errors, and you have R installed,
along with all the packages required to run this code, and a copy of the
*stormwavecluster* git repository, then you should be able to re-run the
analysis here by simply copy-pasting the code. Alternatively, it can be run
with the `knit` command in the *knitr* package: 

```r
library(knitr)
knit('statistical_model_vine_copula.Rmd')
```

The basic approach followed here is to:
* **Step 1: Load the previous session**

Later we will **FIXME:POINT TO SUBSEQUENT STEPS**

# **Step 1: Load the previous session**
---------------------------------------

Here we load the session previously derived by 
[statistical_model_univariate_distributions.Rmd](statistical_model_univariate_distributions.Rmd).
As before, the code is adjusted to optionally use an analysis with random
perturbations to the data, to check the impact of ties and data discretization.


```r
# Need to re-load packages, as R does not automatically do this when re-loading
# a session
library(evmix)
```

```
## Loading required package: MASS
```

```
## Loading required package: splines
```

```
## Loading required package: gsl
```

```
## Loading required package: SparseM
```

```
## 
## Attaching package: 'SparseM'
```

```
## The following object is masked from 'package:base':
## 
##     backsolve
```

```r
library(logspline)
library(CDVine) # Used to set structure of C-Vine copula
```

```
## The CDVine package is no longer developed actively.
## Please consider using the more general VineCopula package
## (see https://CRAN.R-project.org/package=VineCopula),
## which extends and improves the functionality of CDVine.
```

```r
library(VineCopula) # Main copula fitting routine. 
```

```
## 
## Attaching package: 'VineCopula'
```

```
## The following objects are masked from 'package:CDVine':
## 
##     BiCopCDF, BiCopChiPlot, BiCopEst, BiCopHfunc, BiCopIndTest,
##     BiCopKPlot, BiCopLambda, BiCopMetaContour, BiCopName,
##     BiCopPar2TailDep, BiCopPar2Tau, BiCopPDF, BiCopSelect,
##     BiCopSim, BiCopTau2Par, BiCopVuongClarke
```

```r
# Here we support multiple runs with random tie-breaking of the data
# If R was passed a commandline argument 'break_ties n' on startup (with n = integer),
# then read the n'th R session matching 'Rimages/session_storm_timings_TRUE_*.Rdata'.
# That session will correspond to one of the tie-breaking sessions
if( length(grep('break_ties', commandArgs(trailingOnly=TRUE))) > 0 ){

    # Read one of the sessions with tie-breaking
    session_n = as.numeric(commandArgs(trailingOnly=TRUE)[2])
    previous_R_session_file = Sys.glob(
        'Rimages/session_univariate_distributions_TRUE_*.Rdata')[session_n]

}else{

    # Read the session that does not do any tie-breaking
    previous_R_session_file = 'Rimages/session_univariate_distributions_FALSE_0.Rdata'

}

load(previous_R_session_file)

# Definitions controlling the number of years in the synthetic series
nyears_synthetic_series = 1e+06 

print(previous_R_session_file)
```

```
## [1] "Rimages/session_storm_timings_FALSE_0.Rdata"
```

# **Step 2: Fit a copula to the remaining dependencies in the multivariate storm data**

Here we use a copula to model the remaining dependencies in the storm data. The
first step, implemented below, is to **convert our data to inverse quantile
values, conditional on the event time of year (`startyear`) and the `soiA`
value**. We only use rows which do not containing missing observations. Note
how the conditional inverse quantile values are defined using the previously
derived conditional distribution fits.


```r
# Extract observations which will be used in the copula. Ignore observations
# with any missing value, since the code cannot treat these.
skip = which(rowSums(is.na(
    event_statistics[c('duration', 'hsig', 'dir', 'tideResid', 'steepness', 'soiA')])) > 0)

# Transform the data into uniformly distributed margins, accounting
# for the value of any conditional variables
conditional_variables = list(startyear=event_statistics$startyear[-skip], 
    soiA=event_statistics$soiA[-skip])

events_conditional_copuladata = data.frame( 
    duration=duration_fit_conditional$pfun(
        event_statistics$duration[-skip], conditional_variables),
    hsig=hsig_fit_conditional$pfun(event_statistics$hsig[-skip], 
        conditional_variables),
    dir=dir_fit_conditional$pfun(event_statistics$dir[-skip], 
        conditional_variables),
    tideResid=tideResid_fit_conditional$pfun(event_statistics$tideResid[-skip], 
        conditional_variables),
    steepness=steepness_fit_conditional$pfun(event_statistics$steepness[-skip],
        conditional_variables)
    )

rm(skip, conditional_variables)


# Plot the transformed data
DU$nice_pairs(events_conditional_copuladata)
```

![plot of chunk empiricalCopulaData](figure/empiricalCopulaData-1.png)

We can also plot contours of the empirical copula. **The following plot depicts in a
non-parametric way the sorts of relationships that we will be modelling with
copulas.**

```r
# Convert our data to an object of type 'copuladata', to use the VineCopula
# functions
es_cop = as.copuladata(events_conditional_copuladata)
pairs(es_cop)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

**Here we repeat an earlier plot showing the dependencies in the full dataset**. This
is the most useful reference against which to compare the model, which is fit subsequently.

```r
# Sort the variables for easy comparison with the plot of the simulated data 
DU$nice_pairs(event_statistics[,c('hsig', 'duration', 'tideResid', 'steepness', 'dir')])
```

![plot of chunk pairsPlot22](figure/pairsPlot22-1.png)

**Below, we fit a C-Vine copula to the data, and plot the simulated results**.
The black points are the data, and the green points are the model. 

```r
# Try C Vine with most strongly correlated nodes fit first:
c_vine_node_order=c('hsig', 'duration', 'tideResid', 'steepness', 'dir')
c_vine_order = match(c_vine_node_order, names(es_cop))
es_cop_reorder = es_cop[,c_vine_order]

#
# Make a function which does the fit and provides a function to randomly
# sample from the copula, including transforming the modelled [0,1] values
# to the storm variable scales
#
make_Rvine_random_sampler<-function(es_cop_reorder, copula_fit=NULL, 
    plot=FALSE, fixed_variable_order=TRUE){
    # Fit an RVine copula to es_cop_reorder
    #
    # Restricting the familyset seems to help with the hsig/duration relation
    # (it can become too scattered otherwise)
    #
    # Optionally we can provide a fitted copula (must be of CVine type based 
    # on RVineStructureSelect). In this case the copula parameters will be
    # re-estimated with MLE to produce copula_fit_mle. This is useful for 
    # bootstrapping

    if(is.null(copula_fit)){
        # Choose the structure of the copula

        if(!fixed_variable_order){
            # Let VineCopula select the order of variables
            copula_fit = RVineStructureSelect(
                es_cop_reorder, 
                indeptest=TRUE, 
                type='RVine', 
                familyset=NA, 
                selectioncrit='AIC')
        }else{

            # Codes for all one-parameter copulas, and the t-copula
            one_par_copulas = c(1:6, 13:14, 16, 23:24, 26, 33:34, 36) 

            # Order the variables the same as our input data
            copula_fit = CDVineCopSelect(
                es_cop_reorder, 
                familyset=one_par_copulas, # 1:6 = all 1 parameter families, without rotations, and t-copula
                type='CVine', 
                indeptest=TRUE,
                selectioncrit="AIC")

            copula_fit = C2RVine(1:5, copula_fit$family, copula_fit$par, 
                copula_fit$par2)
        }

    }else{
        copula_fit = copula_fit
    }

    # Update parameters using maximum likelihood At the time of writing
    # (14/03/2017), the CRAN version of VineCopula has a bug in this routine,
    # but that seems fixed in the github version, which can be obtained using
    # the command:
    # devtools::install_github("tnagler/VineCopula")
    copula_fit_mle = RVineMLE(es_cop_reorder, copula_fit)

    if(plot){
        simdata2 = RVineSim(1e+04, RVM=copula_fit_mle$RVM)
        colnames(simdata2) = names(es_cop_reorder)
        simdata2 = as.data.frame(simdata2)
        sim_full2 = with(simdata2, 
            stormVarFun(duration, hsig, dir, steepness, tideResid))
       
        # Plot it 
        DU$nice_pairs(sim_full2[1:6000, names(es_cop_reorder)], 
            extra_data=event_statistics[names(es_cop_reorder)])
    }

    # Convenience function to sample randomly from the copula
    random_copula_samples<-function(n){
        out = RVineSim(n, RVM=copula_fit_mle$RVM)
        return(out)
    }

    return(environment())
}


copula_model = make_Rvine_random_sampler(es_cop_reorder, plot=TRUE)
```

```
## iter   10 value -475.478183
## final  value -475.488607 
## converged
```

![plot of chunk makeRvineRandomSampler](figure/makeRvineRandomSampler-1.png)

```r
random_copula_samples = copula_model$random_copula_samples 

# Print information
print(copula_model$copula_fit_mle)
```

```
## $value
## [1] 475.4886
## 
## $convergence
## [1] 0
## 
## $message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
## 
## $counts
## function gradient 
##       23       23 
## 
## $RVM
## C-vine copula with the following pair-copulas:
## Tree 1:
## 1,5  Independence 
## 1,4  Gaussian (par = 0.35, tau = 0.23) 
## 1,3  Gaussian (par = 0.54, tau = 0.37) 
## 1,2  Survival Gumbel (par = 2.45, tau = 0.59) 
## 
## Tree 2:
## 2,5;1  Frank (par = -1.3, tau = -0.14) 
## 2,4;1  Independence 
## 2,3;1  Gaussian (par = 0.21, tau = 0.14) 
## 
## Tree 3:
## 3,5;2,1  Frank (par = 1.55, tau = 0.17) 
## 3,4;2,1  Frank (par = -0.54, tau = -0.06) 
## 
## Tree 4:
## 4,5;3,2,1  Independence
```

**Below, we plot the empirical contours of the randomly simulated points, for
comparison with a similar plot of the data given earlier**

```r
# Order the plot to match names(events_conditional_copuladata
pairs(as.copuladata(
    copula_model$simdata2[,names(events_conditional_copuladata)]
    ))
```

![plot of chunk pairsOfModel](figure/pairsOfModel-1.png)

**Here we fit a more complex copula**. This one does not use a pre-specified
order, and allows for more copula families. *It requires use of the github
version of the `VineCopula` package, to work-around a bug in `VineCopula`*
This can be installed from github after the `devtools` package in installed.

```r
# You only need to do this once
devtools::install_github("tnagler/VineCopula")
library(VineCopula)
```


```r
copula_model2 = make_Rvine_random_sampler(es_cop_reorder, plot=TRUE, 
    fixed_variable_order=FALSE)
```

```
## iter   10 value -484.032091
## iter   20 value -484.052037
## iter   30 value -484.077278
## iter   40 value -484.084781
## iter   50 value -484.096117
## iter   60 value -484.101422
## iter   70 value -484.104479
## iter   80 value -484.105172
## final  value -484.105698 
## converged
```

![plot of chunk copula_alternative](figure/copula_alternative-1.png)

```r
# Print information on the fit
print(copula_model2$copula_fit_mle)
```

```
## $value
## [1] 484.1057
## 
## $convergence
## [1] 0
## 
## $message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
## 
## $counts
## function gradient 
##      100      100 
## 
## $RVM
## R-vine copula with the following pair-copulas:
## Tree 1:
## 1,4  Gaussian (par = 0.35, tau = 0.23) 
## 2,1  Survival BB8 (par = 5.93, par2 = 0.82, tau = 0.61) 
## 2,3  Gaussian (par = 0.53, tau = 0.36) 
## 5,2  Rotated Tawn type 2 90 degrees (par = -1.49, par2 = 0.15, tau = -0.09) 
## 
## Tree 2:
## 2,4;1  Independence 
## 3,1;2  Tawn  type 2 (par = 1.27, par2 = 0.38, tau = 0.12) 
## 5,3;2  Survival BB8 (par = 1.63, par2 = 0.88, tau = 0.17) 
## 
## Tree 3:
## 3,4;2,1  Frank (par = -0.6, tau = -0.07) 
## 5,1;3,2  Independence 
## 
## Tree 4:
## 5,4;3,2,1  Independence 
## 
## ---
## 1 <-> hsig,   2 <-> duration,
## 3 <-> tideResid,   4 <-> steepness,
## 5 <-> dir
```

As above, here we plot contours of psuedo-observations generated by the model

```r
# Order the plot to match names(events_conditional_copuladata
pairs(as.copuladata(
    copula_model2$simdata2[,names(events_conditional_copuladata)]
    ))
```

![plot of chunk copula_alternative2](figure/copula_alternative2-1.png)
