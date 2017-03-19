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
## theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+0 
##                                                                                                                       100
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
## [1] "All perturbed fits have the same best fit lambda model as obtained from the original data"
```

```r
# Check variations in the model parameters.
# The following only works if all_rate_eqns are identical
all_rate_par= matrix(
    unlist(lapply(list_envs, f<-function(x) x$best_nhp_model$par)), 
    ncol=4, byrow=TRUE)
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
## [1] 2.108547e-04 1.277276e-03 4.030207e-04 7.146903e-05
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
##        V1                   V2                   V3            
##  Min.   :-0.0001612   Min.   :-0.0043260   Min.   :-1.678e-03  
##  1st Qu.: 0.0001557   1st Qu.:-0.0012330   1st Qu.:-1.098e-03  
##  Median : 0.0003181   Median :-0.0004562   Median :-7.716e-04  
##  Mean   : 0.0003152   Mean   :-0.0004547   Mean   :-7.786e-04  
##  3rd Qu.: 0.0004624   3rd Qu.: 0.0002504   3rd Qu.:-4.327e-04  
##  Max.   : 0.0007224   Max.   : 0.0035070   Max.   : 1.145e-07  
##        V4            
##  Min.   :-2.997e-04  
##  1st Qu.:-7.485e-05  
##  Median :-2.776e-05  
##  Mean   :-3.539e-05  
##  3rd Qu.: 6.544e-06  
##  Max.   : 1.617e-04
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
##  Min.   :3.143e-07   Min.   :3.324e-06   Min.   :1.145e-07  
##  1st Qu.:1.640e-04   1st Qu.:3.975e-04   1st Qu.:4.327e-04  
##  Median :3.181e-04   Median :8.922e-04   Median :7.716e-04  
##  Mean   :3.220e-04   Mean   :1.046e-03   Mean   :7.786e-04  
##  3rd Qu.:4.624e-04   3rd Qu.:1.508e-03   3rd Qu.:1.098e-03  
##  Max.   :7.224e-04   Max.   :4.326e-03   Max.   :1.678e-03  
##        V4           
##  Min.   :4.606e-07  
##  1st Qu.:1.996e-05  
##  Median :4.330e-05  
##  Mean   :5.888e-05  
##  3rd Qu.:8.676e-05  
##  Max.   :2.997e-04
```
