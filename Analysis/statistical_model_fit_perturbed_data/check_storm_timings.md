# How is the storm timing model fit affected by data perturbation?
------------------------------------------------------------------------

The code below can be run after
[batch_run_storm_timings.R](batch_run_storm_timings.R) has successfully been
run, and produced 100 Rdata files containing fits with perturbed input data.

For our data, the code below should show that the optimal model type (`best_nhp_model$rate_equation`)
for the storm event times model is not affected by data perturbation.

Further, the parameter values (`best_nhp_model$par`) only change by (at most) a
few parts in 1000, and often by much less.

Thus the discretization of our data does not seem to be significantly affecting
the storm timing model fit.


```r
#
# Check the results
#
original_storm_timing_session = new.env()
load('../statistical_model_fit/Rimages/session_storm_timings_FALSE_0.Rdata', 
    envir=original_storm_timing_session)

# Read each R session into its own environment
all_storm_timing_sessions = Sys.glob('Rimages/session_storm_timings_TRUE_*.Rdata')
list_envs = list()
for(i in 1:length(all_storm_timing_sessions)){
    session_file = all_storm_timing_sessions[i]
    list_envs[[session_file]] = new.env()
    load(session_file, envir = list_envs[[session_file]])
}

# Check the 'best fit' rate equations in each R session
all_rate_eqns = unlist(lapply(list_envs, f<-function(x) x$best_nhp_model$rate_equation))

# If the chosen model is unaffected by data perturbations, then only one
# equation should appear here.
print('The table of best fit lambda models (only one is expected)')
```

```
## [1] "The table of best fit lambda models (only one is expected)"
```

```r
print(table(all_rate_eqns))
```

```
## all_rate_eqns
##                                                                theta[1] + 0*t+theta[1 + 1]*sin(2*pi*(t - theta[1 + 2]))+0 
##                                                                                                                         5 
## theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+0 
##                                                                                                                        88 
##           theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*sin(2*pi*(t - theta[2 + 2]))+0 
##                                                                                                                         7
```

```r
if(length(unique(all_rate_eqns)) > 1){
    msg = paste0('More than one rate model identified with perturbed data.\n',
        ' The code below must be changed to deal with this case')
    stop(msg)
}else{
    stopifnot(all(all_rate_eqns == original_storm_timing_session$best_nhp_model$rate_equation))
    print('All perturbed fits have the same best fit lambda model as obtained from the original data')
}
```

```
## Error in eval(expr, envir, enclos): More than one rate model identified with perturbed data.
##  The code below must be changed to deal with this case
```

```r
# Check variations in the model parameters.
# The following only works if all_rate_eqns are identical
all_rate_par= matrix(
    unlist(lapply(list_envs, f<-function(x) x$best_nhp_model$par)), 
    ncol=4, byrow=TRUE)
```

```
## Warning in matrix(unlist(lapply(list_envs, f <- function(x) x
## $best_nhp_model$par)), : data length [395] is not a sub-multiple or
## multiple of the number of rows [99]
```

```r
# Coefficient of variation of estimates. Seems to be very small (e.g. 1/1000)
all_rate_CoV = apply(all_rate_par, 2, sd)/apply(all_rate_par, 2, mean)

print('Coefficient of variation of all perturbed model parameters: ')
```

```
## [1] "Coefficient of variation of all perturbed model parameters: "
```

```r
print(all_rate_CoV)
```

```
## [1] 1.0454750 0.9355830 1.0549814 0.9127257
```

```r
all_rate_err = all_rate_par
for(i in 1:ncol(all_rate_err)){ 
    orig_par = original_storm_timing_session$best_nhp_model$par[i]
    all_rate_err[,i] = (all_rate_err[,i] - orig_par)/orig_par
}

print('Summary of [perturbed - original]/original for best fit lambda model parameters')
```

```
## [1] "Summary of [perturbed - original]/original for best fit lambda model parameters"
```

```r
summary(all_rate_err)
```

```
##        V1                   V2                  V3           
##  Min.   :-0.9859492   Min.   : -0.01957   Min.   :-0.987008  
##  1st Qu.:-0.9856313   1st Qu.:  0.00078   1st Qu.:-0.979240  
##  Median :-0.9694841   Median : 28.98113   Median :-0.971860  
##  Mean   :-0.4991508   Mean   : 37.57359   Mean   :-0.538251  
##  3rd Qu.: 0.0002978   3rd Qu.: 74.42760   3rd Qu.:-0.001102  
##  Max.   : 0.5418070   Max.   :106.33387   Max.   : 0.421852  
##        V4          
##  Min.   :-0.53838  
##  1st Qu.:-0.00009  
##  Median :31.76921  
##  Mean   :18.25931  
##  3rd Qu.:31.78792  
##  Max.   :49.52843
```

```r
print('Summary of ABS[perturbed - original]/original for best fit lambda model parameters')
```

```
## [1] "Summary of ABS[perturbed - original]/original for best fit lambda model parameters"
```

```r
summary(abs(all_rate_err))
```

```
##        V1                  V2                  V3           
##  Min.   :0.0000118   Min.   :  0.00009   Min.   :0.0000125  
##  1st Qu.:0.0004984   1st Qu.:  0.00216   1st Qu.:0.0012915  
##  Median :0.9694841   Median : 28.98113   Median :0.9718604  
##  Mean   :0.5534357   Mean   : 37.57452   Mean   :0.5549797  
##  3rd Qu.:0.9856313   3rd Qu.: 74.42760   3rd Qu.:0.9792399  
##  Max.   :0.9859492   Max.   :106.33387   Max.   :0.9870077  
##        V4          
##  Min.   : 0.00000  
##  1st Qu.: 0.00011  
##  Median :31.76921  
##  Mean   :18.43055  
##  3rd Qu.:31.78792  
##  Max.   :49.52843
```
