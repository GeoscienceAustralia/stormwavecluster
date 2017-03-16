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

## If the chosen model is unaffected by data perturbations, then only one
## equation should appear here.
# print('The table of best fit lambda models (only one is expected)')
# print(table(all_rate_eqns))

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
## [1] 2.484223e-04 1.016987e-03 4.619618e-04 7.763279e-05
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
##  Min.   :-0.0003300   Min.   :-2.590e-03   Min.   :-0.0015596  
##  1st Qu.: 0.0000609   1st Qu.:-1.275e-03   1st Qu.:-0.0011230  
##  Median : 0.0001983   Median :-7.787e-04   Median :-0.0007685  
##  Mean   : 0.0002620   Mean   :-6.003e-04   Mean   :-0.0007324  
##  3rd Qu.: 0.0004976   3rd Qu.: 3.918e-05   3rd Qu.:-0.0004373  
##  Max.   : 0.0011129   Max.   : 2.061e-03   Max.   : 0.0004262  
##        V4            
##  Min.   :-1.836e-04  
##  1st Qu.:-7.892e-05  
##  Median :-3.949e-05  
##  Mean   :-3.133e-05  
##  3rd Qu.: 1.744e-05  
##  Max.   : 2.846e-04
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
##  Min.   :4.812e-06   Min.   :3.827e-05   Min.   :1.594e-06  
##  1st Qu.:7.403e-05   1st Qu.:4.997e-04   1st Qu.:4.394e-04  
##  Median :2.056e-04   Median :1.007e-03   Median :7.685e-04  
##  Mean   :2.796e-04   Mean   :1.002e-03   Mean   :7.555e-04  
##  3rd Qu.:4.976e-04   3rd Qu.:1.519e-03   3rd Qu.:1.123e-03  
##  Max.   :1.113e-03   Max.   :2.590e-03   Max.   :1.560e-03  
##        V4           
##  Min.   :2.320e-06  
##  1st Qu.:2.673e-05  
##  Median :5.622e-05  
##  Mean   :6.704e-05  
##  3rd Qu.:9.398e-05  
##  Max.   :2.846e-04
```
