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
## [1] 2.292026e-04 1.045050e-03 4.424752e-04 7.425761e-05
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
##  Min.   :-3.300e-04   Min.   :-3.515e-03   Min.   :-0.0016782  
##  1st Qu.: 8.576e-05   1st Qu.:-1.320e-03   1st Qu.:-0.0010920  
##  Median : 2.303e-04   Median :-8.193e-04   Median :-0.0007634  
##  Mean   : 2.737e-04   Mean   :-6.570e-04   Mean   :-0.0007496  
##  3rd Qu.: 4.799e-04   3rd Qu.: 4.534e-05   3rd Qu.:-0.0004432  
##  Max.   : 1.113e-03   Max.   : 2.061e-03   Max.   : 0.0004262  
##        V4            
##  Min.   :-2.482e-04  
##  1st Qu.:-8.466e-05  
##  Median :-4.487e-05  
##  Mean   :-4.024e-05  
##  3rd Qu.:-8.857e-07  
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
##  Min.   :5.423e-07   Min.   :3.367e-06   Min.   :3.712e-07  
##  1st Qu.:1.022e-04   1st Qu.:5.098e-04   1st Qu.:4.432e-04  
##  Median :2.327e-04   Median :1.003e-03   Median :7.634e-04  
##  Mean   :2.892e-04   Mean   :1.034e-03   Mean   :7.614e-04  
##  3rd Qu.:4.799e-04   3rd Qu.:1.447e-03   3rd Qu.:1.092e-03  
##  Max.   :1.113e-03   Max.   :3.515e-03   Max.   :1.678e-03  
##        V4           
##  Min.   :2.044e-06  
##  1st Qu.:3.395e-05  
##  Median :5.622e-05  
##  Mean   :6.823e-05  
##  3rd Qu.:9.642e-05  
##  Max.   :2.846e-04
```
