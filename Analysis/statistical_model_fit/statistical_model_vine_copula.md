
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
# If running via knitr, ensure knitr halts on error [do not use this command if
# copy-pasting the code]
opts_knit$set(stop_on_error=2L)

# Check that the pre-requisites exist
if(!file.exists('../statistical_model_fit/Rimages/session_univariate_distributions_FALSE_0.Rdata')){
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
The above command produces a .md file with the same name for viewing, which includes
updated figures and print-outs.

The basic approach followed here is to:
* **Step 1: Load the previous session**
* **Step 2: Fit a copula to the remaining dependencies in the multivariate storm data**
* **Step 3: Use the fitted model to simulate a long synthetic storm time-series**

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
## Loading required package: methods
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

#' Function to fit the copula
#'
#' Fit a copula to the psuedo-observations, and provides a function to randomly
#' sample from the copula (including transforming the modelled [0,1] values
#' to the storm variable scales)
#
#' @param es_cop_reorder copuladata.
#' @param copula_fit an existing copula fit, which will be updated using MLE on
#' es_cop_reorder. If NULL, then the copula structure is selected algorithmically
#' @param plot logical. Make a plot?
#' @param cvine_with_fixed_variable_order logical. If copula_fit = NULL, then if
#' this variable is TRUE, we select a C-Vine copula structure, with the variables
#' treated in the order they appear in es_cop_reorder, and a limited available set of
#' bivariate copula families. Otherwise, we use a (more generic) R-Vine copula.
#' @return the function environment
#'
make_Rvine_random_sampler<-function(es_cop_reorder, copula_fit=NULL, 
    plot=FALSE, cvine_with_fixed_variable_order=TRUE){
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
        # The copula structure was not provided. 
        # So choose the structure of the copula

        if(!cvine_with_fixed_variable_order){

            # Use a full 'RVine' copula to model the data
            #
            # This approach is more general than the CVine, but this also means
            # that more is demanded of the automated techniques to find the
            # best copula. 

            # Let VineCopula select the order of variables
            copula_fit = RVineStructureSelect(
                es_cop_reorder, 
                indeptest=TRUE, 
                type='RVine', 
                familyset=NA, 
                selectioncrit='AIC')

        }else{

            #
            # Use a restricted 'CVine' approach. 
            #
            # We specify the order of the variables [same order as columns in
            # es_cop_reorder]. We also only use a subset of the available
            # copula families. 
            #

            # Codes for all one-parameter copulas
            one_par_copulas = c(1, 3:6, 13:14, 16, 23:24, 26, 33:34, 36) 
            simple_copula_families =  c(1, 3:6) #= all 1 parameter families, without rotations

            # Copula selection using the package CDVine.
            # Order the variables the same as our input data
            copula_fit = CDVineCopSelect(
                es_cop_reorder, 
                familyset=one_par_copulas, #simple_copula_families
                type='CVine', 
                indeptest=TRUE,
                selectioncrit="AIC")

            # Convert to an object of type 'RVine', so we can use functions from
            # VineCopula
            copula_fit = C2RVine(1:5, copula_fit$family, copula_fit$par, 
                copula_fit$par2)

        }

    }else{
        # In this case the copula structure was already provided
        copula_fit = copula_fit
    }

    # Update parameters using maximum likelihood At the time of writing
    # (14/03/2017), the CRAN version of VineCopula has a bug in this routine,
    # but that is fixed in the current github version, which can be obtained using
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
## iter   10 value -478.207281
## final  value -478.211150 
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
## [1] 478.2112
## 
## $convergence
## [1] 0
## 
## $message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
## 
## $counts
## function gradient 
##       16       16 
## 
## $RVM
## C-vine copula with the following pair-copulas:
## Tree 1:
## 1,5  Independence 
## 1,4  Gaussian (par = 0.35, tau = 0.23) 
## 1,3  Gaussian (par = 0.54, tau = 0.37) 
## 1,2  Survival Gumbel (par = 2.46, tau = 0.59) 
## 
## Tree 2:
## 2,5;1  Frank (par = -1.29, tau = -0.14) 
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

**Here we plot the Vine structure of the copula**

```r
par(mfrow=c(2,3))
plot(copula_model$copula_fit_mle$RVM, type=2, edge.labels='family')
```

![plot of chunk treeplot](figure/treeplot-1.png)

**Here we statistically test for equality between a sample from the data
and a random sample from the model**. The computational demands of this
test grow rapidly with sample size (e.g. taking 1.3s for n=50, 20 s for n=100, 368 s for
n=200, ...), so here we only use a sample size of 100.

```r
# Find rows of event_statistics with no NA values, since we can't have NA's for this test
non_na_es = apply(!is.na(event_statistics[,c_vine_node_order]), 1, f<-function(x) all(x))
non_na_es = which(non_na_es)

library(TwoCop)
test_data_size = 1:100
if(!break_ties_with_jitter){
    # We must break ties in the data, if not already done, since otherwise the test
    # does not work [it's based on ranks]
    twocopula_test = TwoCop( 
        jitter(as.matrix(event_statistics[non_na_es[test_data_size], c_vine_node_order]), amount=1e-05),
        as.matrix(copula_model$sim_full2[test_data_size, c_vine_node_order]))
}else{
    # There are no ties
    twocopula_test = TwoCop( 
        as.matrix(event_statistics[non_na_es[test_data_size], c_vine_node_order]),
        as.matrix(copula_model$sim_full2[test_data_size, c_vine_node_order]))
}

# Print it out
print(twocopula_test)
```

```
## $pvalue
## [1] 0.69
## 
## $cvm
## [1] 0.02248761
## 
## $VaR
##        95% 
## 0.05884563 
## 
## $cvmsim
##   [1] 0.02739849 0.05875520 0.02310483 0.03240377 0.04133967 0.02196260
##   [7] 0.01804450 0.02454822 0.01978256 0.04955802 0.01960030 0.02399962
##  [13] 0.01851748 0.02229193 0.03094275 0.02064604 0.03211878 0.03065228
##  [19] 0.02053376 0.02119630 0.02331550 0.02163053 0.01919836 0.05133067
##  [25] 0.02454309 0.02486719 0.01992090 0.02732295 0.03339164 0.03766111
##  [31] 0.03870839 0.03235224 0.02187605 0.01867453 0.02787705 0.01859570
##  [37] 0.05878386 0.02336605 0.01542298 0.01797721 0.02279818 0.03229862
##  [43] 0.02613078 0.05179731 0.02764380 0.03025439 0.02941385 0.02454943
##  [49] 0.02648642 0.03891873 0.02407681 0.06459222 0.02292677 0.06232210
##  [55] 0.03667034 0.05156617 0.06001927 0.02657297 0.02900734 0.03030116
##  [61] 0.03769127 0.03083017 0.01999761 0.01822445 0.01797604 0.02226903
##  [67] 0.01856260 0.02116888 0.02313712 0.03031876 0.02297130 0.03354682
##  [73] 0.02938124 0.02189255 0.03928390 0.01773387 0.02639150 0.05423201
##  [79] 0.01939419 0.02725468 0.01843285 0.04575102 0.02792443 0.02211097
##  [85] 0.02199447 0.02554436 0.03282494 0.02748423 0.03990957 0.02041514
##  [91] 0.12044191 0.02529272 0.07479675 0.02075966 0.04729980 0.02432983
##  [97] 0.03416098 0.02268565 0.04879826 0.02345337
```

```r
if(twocopula_test$pvalue > 0.05){
    print('two-copula test DOES NOT REJECT null hypothesis at 5% level')
}else{
    print('two-copula test REJECTS null hypothesis at 5% level')
}
```

```
## [1] "two-copula test DOES NOT REJECT null hypothesis at 5% level"
```

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
    cvine_with_fixed_variable_order=FALSE)
```

```
## iter   10 value -484.261091
## iter   20 value -484.264838
## iter   30 value -484.288973
## iter   40 value -484.298080
## iter   50 value -484.328151
## iter   60 value -484.337849
## iter   70 value -484.339432
## iter   80 value -484.339928
## iter   90 value -484.340345
## final  value -484.340411 
## converged
```

![plot of chunk copula_alternative](figure/copula_alternative-1.png)

```r
# Print information on the fit
print(copula_model2$copula_fit_mle)
```

```
## $value
## [1] 484.3404
## 
## $convergence
## [1] 0
## 
## $message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
## 
## $counts
## function gradient 
##      104      104 
## 
## $RVM
## R-vine copula with the following pair-copulas:
## Tree 1:
## 1,4  Gaussian (par = 0.35, tau = 0.23) 
## 2,1  Survival BB8 (par = 5.82, par2 = 0.83, tau = 0.61) 
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

As above, here we plot contours of psuedo-observations generated by the more
complex model.

```r
# Order the plot to match names(events_conditional_copuladata
pairs(as.copuladata(
    copula_model2$simdata2[,names(events_conditional_copuladata)]
    ))
```

![plot of chunk copula_alternative2](figure/copula_alternative2-1.png)

**Here we plot the Vine structure of the more complex copula.**

```r
par(mfrow=c(2,3))
plot(copula_model2$copula_fit_mle$RVM, type=2, edge.labels='family')
```

![plot of chunk treeplot2](figure/treeplot2-1.png)

**Here we do the test above**

```r
if(!break_ties_with_jitter){
    # We must break ties in the data, if not already done, since otherwise the test
    # does not work (it's based on ranks)
    twocopula_testB = TwoCop( 
        jitter(as.matrix(event_statistics[non_na_es[test_data_size], c_vine_node_order]), amount=1e-05),
        as.matrix(copula_model2$sim_full2[test_data_size, c_vine_node_order]))
}else{
    # There are no ties
    twocopula_testB = TwoCop( 
        as.matrix(event_statistics[non_na_es[test_data_size], c_vine_node_order]),
        as.matrix(copula_model2$sim_full2[test_data_size, c_vine_node_order]))
}

# Print it out
print(twocopula_testB)
```

```
## $pvalue
## [1] 0.5
## 
## $cvm
## [1] 0.02560812
## 
## $VaR
##        95% 
## 0.04977946 
## 
## $cvmsim
##   [1] 0.03602990 0.03873566 0.02060185 0.04471292 0.01666976 0.03798539
##   [7] 0.01882996 0.02261145 0.02397355 0.04028588 0.02065314 0.02949590
##  [13] 0.04958978 0.05535840 0.02126548 0.02114225 0.02428030 0.03207685
##  [19] 0.02814026 0.02381844 0.02870713 0.04235719 0.02083069 0.02253649
##  [25] 0.01650364 0.01999840 0.02691909 0.03530230 0.03399995 0.02985775
##  [31] 0.02129359 0.03541424 0.04317802 0.04883129 0.04132254 0.01827753
##  [37] 0.06719158 0.02544377 0.02933784 0.02039111 0.02304718 0.02034132
##  [43] 0.01701568 0.02389673 0.01543827 0.02571380 0.01746617 0.02302817
##  [49] 0.02598294 0.02287361 0.01273191 0.06111369 0.03517542 0.02032321
##  [55] 0.02125867 0.03663157 0.02029448 0.01626011 0.03535930 0.02704746
##  [61] 0.03312930 0.02731390 0.03368083 0.04406091 0.04227164 0.03141644
##  [67] 0.02141410 0.06375702 0.01873181 0.01999440 0.03081604 0.01963916
##  [73] 0.01445480 0.02772478 0.03280734 0.01253570 0.02199733 0.02414705
##  [79] 0.03309409 0.01819545 0.04240291 0.01810118 0.01598667 0.02425463
##  [85] 0.02370813 0.03964503 0.05338341 0.02642563 0.04372382 0.02297220
##  [91] 0.03516655 0.02169600 0.03237378 0.02370677 0.03462840 0.03005128
##  [97] 0.02462447 0.01527409 0.02045260 0.04557866
```

```r
if(twocopula_testB$pvalue > 0.05){
    print('two-copula test DOES NOT REJECT null hypothesis at 5% level')
}else{
    print('two-copula test REJECTS null hypothesis at 5% level')
}
```

```
## [1] "two-copula test DOES NOT REJECT null hypothesis at 5% level"
```
## Save the Rimage for later use

We use the same `run_title_id` as was computed in the previous 2 sections
([statistical_model_storm_timings.md](statistical_model_storm_timings.md), and
[statistical_model_univariate_distributions.md](statistical_model_univariate_distributions.md)).

```r
dir.create('Rimages', showWarnings=FALSE)
Rimage_title = paste0('Rimages/session_vine_copula_', run_title_id, '.Rdata')
save.image(Rimage_title)
```



