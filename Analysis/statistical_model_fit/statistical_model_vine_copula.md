
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
## iter   10 value -478.761494
## final  value -478.761554 
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
## [1] 478.7616
## 
## $convergence
## [1] 0
## 
## $message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
## 
## $counts
## function gradient 
##       14       14 
## 
## $RVM
## C-vine copula with the following pair-copulas:
## Tree 1:
## 1,5  Independence 
## 1,4  Gaussian (par = 0.35, tau = 0.23) 
## 1,3  Gaussian (par = 0.54, tau = 0.36) 
## 1,2  Survival Gumbel (par = 2.46, tau = 0.59) 
## 
## Tree 2:
## 2,5;1  Frank (par = -1.25, tau = -0.14) 
## 2,4;1  Independence 
## 2,3;1  Gaussian (par = 0.21, tau = 0.13) 
## 
## Tree 3:
## 3,5;2,1  Frank (par = 1.53, tau = 0.17) 
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
## [1] 0.45
## 
## $cvm
## [1] 0.02568118
## 
## $VaR
##        95% 
## 0.04531514 
## 
## $cvmsim
##   [1] 0.03345984 0.02451727 0.03221835 0.01746202 0.02768248 0.02496632
##   [7] 0.03178810 0.02351762 0.01947479 0.03699932 0.02615755 0.03205686
##  [13] 0.02418097 0.02557383 0.04703020 0.03232856 0.01995323 0.02162552
##  [19] 0.01959892 0.01583444 0.03165481 0.02100903 0.03047798 0.02637523
##  [25] 0.02126720 0.02688547 0.03211258 0.01374130 0.02528536 0.01365627
##  [31] 0.03799811 0.02121901 0.01888583 0.01200331 0.01752408 0.03424945
##  [37] 0.01491876 0.05300675 0.04016222 0.11848757 0.01474694 0.03801706
##  [43] 0.03194228 0.02454352 0.01696684 0.02653506 0.01868602 0.03158114
##  [49] 0.04261334 0.02113892 0.01871162 0.02526573 0.05248468 0.02882710
##  [55] 0.02248453 0.01612417 0.03693360 0.01685173 0.02172421 0.01886014
##  [61] 0.02150356 0.02722538 0.01748865 0.03441657 0.01629330 0.02889078
##  [67] 0.03269380 0.01868022 0.02127569 0.02237868 0.02668643 0.02500072
##  [73] 0.02153187 0.02141166 0.01899468 0.01852368 0.02926237 0.01758371
##  [79] 0.03173380 0.03523403 0.01662757 0.03189625 0.01650450 0.01833167
##  [85] 0.02501182 0.04197092 0.04530583 0.02832155 0.02137582 0.02259909
##  [91] 0.01811316 0.02437129 0.02739377 0.03632064 0.02232444 0.02881066
##  [97] 0.01792063 0.03482033 0.03688827 0.04549211
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
## iter   10 value -487.772637
## iter   20 value -487.777201
## iter   30 value -487.798429
## iter   40 value -487.849119
## iter   50 value -487.850526
## final  value -487.850561 
## converged
```

![plot of chunk copula_alternative](figure/copula_alternative-1.png)

```r
# Print information on the fit
print(copula_model2$copula_fit_mle)
```

```
## $value
## [1] 487.8506
## 
## $convergence
## [1] 0
## 
## $message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
## 
## $counts
## function gradient 
##       62       62 
## 
## $RVM
## R-vine copula with the following pair-copulas:
## Tree 1:
## 1,4  Gaussian (par = 0.35, tau = 0.23) 
## 2,1  Survival BB8 (par = 6, par2 = 0.82, tau = 0.61) 
## 2,3  Gaussian (par = 0.53, tau = 0.36) 
## 5,2  Rotated Tawn type 2 90 degrees (par = -1.5, par2 = 0.16, tau = -0.1) 
## 
## Tree 2:
## 2,4;1  Independence 
## 3,1;2  Tawn  type 2 (par = 1.26, par2 = 0.38, tau = 0.12) 
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
        jitter(as.matrix(event_statistics[non_na_es[test_data_size], c_vine_node_order]), 
            amount=1e-05),
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
## [1] 0.57
## 
## $cvm
## [1] 0.02558276
## 
## $VaR
##        95% 
## 0.04170438 
## 
## $cvmsim
##   [1] 0.02043238 0.02707432 0.02188912 0.04002636 0.02925130 0.02524008
##   [7] 0.02936260 0.06360737 0.01547817 0.02966861 0.03612714 0.02675758
##  [13] 0.03569780 0.01935043 0.01982385 0.02684007 0.02985134 0.03420338
##  [19] 0.02537852 0.03121681 0.02657861 0.03437038 0.06863088 0.03351921
##  [25] 0.02615569 0.03150437 0.03402662 0.02375997 0.03414153 0.05734665
##  [31] 0.04092620 0.02142769 0.02518727 0.03156873 0.03158061 0.03064303
##  [37] 0.02913235 0.02124264 0.03012918 0.02312832 0.01526617 0.02456140
##  [43] 0.01983034 0.02737334 0.03397487 0.02816920 0.06612399 0.02118778
##  [49] 0.03752198 0.01825508 0.03246849 0.03244686 0.02723314 0.03870146
##  [55] 0.04073818 0.02211299 0.03364747 0.02279555 0.02502559 0.02081524
##  [61] 0.03274997 0.02347108 0.03274747 0.02013345 0.01852908 0.02991977
##  [67] 0.05648982 0.02163212 0.03282419 0.02224332 0.01467618 0.01778747
##  [73] 0.03910438 0.03822065 0.02510423 0.02727314 0.01829525 0.02171495
##  [79] 0.01438497 0.02300971 0.02411677 0.01553040 0.03243996 0.02859000
##  [85] 0.02639537 0.03010746 0.02454911 0.02074967 0.01870788 0.02025841
##  [91] 0.03035207 0.03084688 0.02310682 0.03536892 0.01810688 0.03833294
##  [97] 0.02595916 0.01913847 0.02588524 0.02202059
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



