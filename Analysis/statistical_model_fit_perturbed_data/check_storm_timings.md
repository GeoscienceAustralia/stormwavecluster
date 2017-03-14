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
## [1] 2.088343e-04 1.074498e-03 4.237163e-04 6.997662e-05
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
##  Min.   :-0.0001997   Min.   :-3.515e-03   Min.   :-1.678e-03  
##  1st Qu.: 0.0001276   1st Qu.:-1.342e-03   1st Qu.:-1.084e-03  
##  Median : 0.0002506   Median :-8.517e-04   Median :-7.571e-04  
##  Mean   : 0.0002853   Mean   :-7.137e-04   Mean   :-7.669e-04  
##  3rd Qu.: 0.0004616   3rd Qu.: 6.316e-05   3rd Qu.:-4.483e-04  
##  Max.   : 0.0007078   Max.   : 1.676e-03   Max.   : 2.137e-05  
##        V4            
##  Min.   :-2.482e-04  
##  1st Qu.:-9.198e-05  
##  Median :-4.909e-05  
##  Mean   :-4.916e-05  
##  3rd Qu.:-7.085e-06  
##  Max.   : 1.701e-04
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
##  Min.   :5.423e-07   Min.   :3.367e-06   Min.   :3.712e-07  
##  1st Qu.:1.358e-04   1st Qu.:5.193e-04   1st Qu.:4.483e-04  
##  Median :2.506e-04   Median :9.655e-04   Median :7.571e-04  
##  Mean   :2.988e-04   Mean   :1.065e-03   Mean   :7.673e-04  
##  3rd Qu.:4.616e-04   3rd Qu.:1.378e-03   3rd Qu.:1.084e-03  
##  Max.   :7.078e-04   Max.   :3.515e-03   Max.   :1.678e-03  
##        V4           
##  Min.   :2.044e-06  
##  1st Qu.:3.606e-05  
##  Median :5.669e-05  
##  Mean   :6.942e-05  
##  3rd Qu.:9.782e-05  
##  Max.   :2.482e-04
```
