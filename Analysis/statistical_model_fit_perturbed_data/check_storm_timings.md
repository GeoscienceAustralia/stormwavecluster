# How is the storm timing model fit affected by data perturbation?
------------------------------------------------------------------------

The code below can be run after
[batch_run_storm_timings.R](batch_run_storm_timings.R) has successfully been
run, and produced 100 Rdata files containing fits with perturbed input data.

It shows the range of 'best-fit' storm timing models, and summarises some of their
parameters.


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

# Check that all fits passed
number_of_fitted_models = sapply(list_envs, f<-function(x) length(x$exhaustive_AICs))
print(summary(number_of_fitted_models))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##      16      16      16      16      16      16
```

```r
if(!( all(number_of_fitted_models == length(original_storm_timing_session$exhaustive_AICs)) )){
    print('WARNING: Some sessions did not fit all models')
}else{
    print('SUCCESS: All sessions fit all models')
}
```

```
## [1] "SUCCESS: All sessions fit all models"
```

```r
# Summary of number of events
nevents = sapply(list_envs, f<-function(x) length(x$event_statistics[,1]))
print(summary(nevents))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   692.0   708.0   712.0   712.1   717.0   730.0
```

```r
# Number of events in original data
print(length(original_storm_timing_session$event_statistics[,1]))
```

```
## [1] 678
```

```r
# Check the 'best fit' rate equations in each R session
all_rate_eqns = unlist(lapply(list_envs, f<-function(x) x$best_nhp_model$rate_equation))

# Look at the range of Lambda models
print('The table of best fit lambda models')
```

```
## [1] "The table of best fit lambda models"
```

```r
print(table(all_rate_eqns))
```

```
## all_rate_eqns
##                                                                                                            theta[1] + 0*t+theta[1 + 1]*abs(2/pi*asin(cos(pi*(t-theta[1+2]))))+0 
##                                                                                                                                                                              23 
##                                                      theta[1] + 0*t+theta[1 + 1]*abs(2/pi*asin(cos(pi*(t-theta[1+2]))))+theta[1 + 2 + 1]*exp((tlast - t)*abs(theta[1 + 2 + 2])) 
##                                                                                                                                                                               2 
##                                                       theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+0 
##                                                                                                                                                                              66 
## theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+theta[2 + 2 + 1]*exp((tlast - t)*abs(theta[2 + 2 + 2])) 
##                                                                                                                                                                               9
```

```r
unique_rate_eqns = unique(all_rate_eqns)
match_rate_eqns = match(all_rate_eqns, unique_rate_eqns)

#
# Check variations in the model parameters.
# Loop over each rate equation separately
#
for(i in 1:length(unique_rate_eqns)){

    # Find runs which selected this rate equation
    keep = which(match_rate_eqns == i)

    print('')
    print('################')
    print(paste0('Equation : ', unique_rate_eqns[i]))
    print(paste0(' .... was selected by: ', length(keep), ' runs'))

    # Get the fitted parameters
    all_rate_par= matrix(
        unlist(lapply(list_envs[keep], f<-function(x) x$best_nhp_model$par)),
        nrow=length(keep), byrow=TRUE)

    # Coefficient of variation of estimates. Seems to be very small (e.g. 1/1000)
    all_rate_CoV = apply(all_rate_par, 2, sd)/apply(all_rate_par, 2, mean)

    print('')
    print('.... Summary of fitted parameters')
    print(summary(all_rate_par))

    print('')
    print('.... Coefficient of variation of perturbed model parameters: ')
    print(all_rate_CoV)

    print('')
    print('.... Parameters / approximate_standard_errors')
    all_rate_se = matrix(
        unlist(lapply(
            list_envs[keep], f<-function(x){
                rate_se =  try(x$nhp$get_fit_standard_errors(x$best_nhp_model))
                if(class(rate_se) == 'try_error'){
                    rate_se = x$best_nhp_model$par * NA
                }
                return(rate_se)
                }
            )),
        nrow=length(keep), byrow=TRUE)
    print(summary(all_rate_par/all_rate_se))
}
```

```
## [1] ""
## [1] "################"
## [1] "Equation : theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+theta[2 + 2 + 1]*exp((tlast - t)*abs(theta[2 + 2 + 2]))"
## [1] " .... was selected by: 9 runs"
## [1] ""
## [1] ".... Summary of fitted parameters"
##        V1              V2               V3              V4        
##  Min.   :17.37   Min.   :0.1917   Min.   :16.67   Min.   :0.4950  
##  1st Qu.:18.03   1st Qu.:0.2128   1st Qu.:18.14   1st Qu.:0.5055  
##  Median :18.24   Median :0.2313   Median :18.45   Median :0.5061  
##  Mean   :18.29   Mean   :0.2271   Mean   :18.47   Mean   :0.5080  
##  3rd Qu.:18.40   3rd Qu.:0.2426   3rd Qu.:19.30   3rd Qu.:0.5098  
##  Max.   :19.31   Max.   :0.2585   Max.   :19.91   Max.   :0.5239  
##        V5                V6       
##  Min.   :-637.53   Min.   : 1891  
##  1st Qu.:-206.51   1st Qu.: 2598  
##  Median : -38.15   Median : 5453  
##  Mean   :-167.68   Mean   :11962  
##  3rd Qu.: -24.27   3rd Qu.:19943  
##  Max.   : -20.54   Max.   :28900  
## [1] ""
## [1] ".... Coefficient of variation of perturbed model parameters: "
## [1]  0.03385630  0.09671043  0.05709182  0.01555927 -1.25028785  0.90523635
## [1] ""
## [1] ".... Parameters / approximate_standard_errors"
##        V1               V2              V3              V4       
##  Min.   : 9.725   Min.   :1.433   Min.   :4.777   Min.   :15.53  
##  1st Qu.: 9.840   1st Qu.:1.575   1st Qu.:5.118   1st Qu.:18.42  
##  Median : 9.992   Median :1.727   Median :5.320   Median :18.94  
##  Mean   :10.047   Mean   :1.701   Mean   :5.283   Mean   :21.06  
##  3rd Qu.:10.115   3rd Qu.:1.816   3rd Qu.:5.463   3rd Qu.:24.46  
##  Max.   :10.514   Max.   :1.956   Max.   :5.732   Max.   :28.28  
##        V5                V6        
##  Min.   :-2.8254   Min.   :0.7818  
##  1st Qu.:-2.5707   1st Qu.:1.2361  
##  Median :-1.9096   Median :1.5428  
##  Mean   :-1.8046   Mean   :1.9094  
##  3rd Qu.:-1.0407   3rd Qu.:2.5578  
##  Max.   :-0.8672   Max.   :4.2497  
## [1] ""
## [1] "################"
## [1] "Equation : theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+0"
## [1] " .... was selected by: 66 runs"
## [1] ""
## [1] ".... Summary of fitted parameters"
##        V1              V2               V3              V4        
##  Min.   :16.88   Min.   :0.1916   Min.   :15.19   Min.   :0.4948  
##  1st Qu.:17.89   1st Qu.:0.2201   1st Qu.:17.61   1st Qu.:0.5011  
##  Median :18.27   Median :0.2333   Median :18.39   Median :0.5089  
##  Mean   :18.26   Mean   :0.2371   Mean   :18.31   Mean   :0.5078  
##  3rd Qu.:18.56   3rd Qu.:0.2583   3rd Qu.:19.02   3rd Qu.:0.5121  
##  Max.   :19.48   Max.   :0.3398   Max.   :20.88   Max.   :0.5269  
## [1] ""
## [1] ".... Coefficient of variation of perturbed model parameters: "
## [1] 0.02878948 0.11916824 0.05871899 0.01786373
## [1] ""
## [1] ".... Parameters / approximate_standard_errors"
```

```
## Warning in sqrt(diag(solve(fit$hessian))): NaNs produced
```

```
## [1] "Warning: standard errors could not be computed from raw Hessian."
## [1] ".... Trying nearPD to get nearest positive definite matrix "
##        V1               V2                 V3              V4           
##  Min.   : 3.594   Min.   :0.001976   Min.   :4.416   Min.   : 0.000406  
##  1st Qu.: 9.871   1st Qu.:1.630985   1st Qu.:5.058   1st Qu.:18.996892  
##  Median :10.023   Median :1.747321   Median :5.259   Median :21.979852  
##  Mean   : 9.938   Mean   :1.750274   Mean   :5.236   Mean   :21.486140  
##  3rd Qu.:10.184   3rd Qu.:1.944316   3rd Qu.:5.433   3rd Qu.:24.495685  
##  Max.   :10.478   Max.   :2.525993   Max.   :6.007   Max.   :30.463352  
## [1] ""
## [1] "################"
## [1] "Equation : theta[1] + 0*t+theta[1 + 1]*abs(2/pi*asin(cos(pi*(t-theta[1+2]))))+theta[1 + 2 + 1]*exp((tlast - t)*abs(theta[1 + 2 + 2]))"
## [1] " .... was selected by: 2 runs"
## [1] ""
## [1] ".... Summary of fitted parameters"
##        V1              V2              V3               V4       
##  Min.   :17.83   Min.   :18.53   Min.   :0.5019   Min.   :-3598  
##  1st Qu.:17.96   1st Qu.:18.59   1st Qu.:0.5040   1st Qu.:-2954  
##  Median :18.09   Median :18.64   Median :0.5061   Median :-2311  
##  Mean   :18.09   Mean   :18.64   Mean   :0.5061   Mean   :-2311  
##  3rd Qu.:18.22   3rd Qu.:18.70   3rd Qu.:0.5083   3rd Qu.:-1667  
##  Max.   :18.35   Max.   :18.76   Max.   :0.5104   Max.   :-1023  
##        V5       
##  Min.   :36136  
##  1st Qu.:39347  
##  Median :42559  
##  Mean   :42559  
##  3rd Qu.:45770  
##  Max.   :48982  
## [1] ""
## [1] ".... Coefficient of variation of perturbed model parameters: "
## [1]  0.020313074  0.008570504  0.011898054 -0.787814505  0.213424097
## [1] ""
## [1] ".... Parameters / approximate_standard_errors"
##        V1               V2              V3              V4         
##  Min.   : 9.911   Min.   :5.258   Min.   :16.69   Min.   :-0.8773  
##  1st Qu.: 9.951   1st Qu.:5.287   1st Qu.:18.11   1st Qu.:-0.7692  
##  Median : 9.991   Median :5.317   Median :19.54   Median :-0.6610  
##  Mean   : 9.991   Mean   :5.317   Mean   :19.54   Mean   :-0.6610  
##  3rd Qu.:10.032   3rd Qu.:5.346   3rd Qu.:20.96   3rd Qu.:-0.5528  
##  Max.   :10.072   Max.   :5.375   Max.   :22.39   Max.   :-0.4446  
##        V5       
##  Min.   :1.723  
##  1st Qu.:2.372  
##  Median :3.022  
##  Mean   :3.022  
##  3rd Qu.:3.672  
##  Max.   :4.322  
## [1] ""
## [1] "################"
## [1] "Equation : theta[1] + 0*t+theta[1 + 1]*abs(2/pi*asin(cos(pi*(t-theta[1+2]))))+0"
## [1] " .... was selected by: 23 runs"
## [1] ""
## [1] ".... Summary of fitted parameters"
##        V1              V2              V3        
##  Min.   :16.13   Min.   :16.94   Min.   :0.4946  
##  1st Qu.:17.59   1st Qu.:18.06   1st Qu.:0.4962  
##  Median :17.91   Median :18.54   Median :0.5019  
##  Mean   :17.86   Mean   :18.55   Mean   :0.5039  
##  3rd Qu.:18.27   3rd Qu.:18.82   3rd Qu.:0.5110  
##  Max.   :18.77   Max.   :21.15   Max.   :0.5141  
## [1] ""
## [1] ".... Coefficient of variation of perturbed model parameters: "
## [1] 0.03365180 0.05247720 0.01449679
## [1] ""
## [1] ".... Parameters / approximate_standard_errors"
##        V1               V2              V3       
##  Min.   : 9.296   Min.   :4.795   Min.   :15.92  
##  1st Qu.: 9.841   1st Qu.:5.181   1st Qu.:21.03  
##  Median : 9.972   Median :5.302   Median :24.04  
##  Mean   : 9.930   Mean   :5.321   Mean   :24.52  
##  3rd Qu.:10.061   3rd Qu.:5.437   3rd Qu.:29.06  
##  Max.   :10.201   Max.   :6.177   Max.   :31.75
```

```r
print('Original model equation and parameters: ')
```

```
## [1] "Original model equation and parameters: "
```

```r
print(original_storm_timing_session$best_nhp_model$rate_equation)
```

```
## [1] "theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+0"
```

```r
print(original_storm_timing_session$best_nhp_model$par)
```

```
## [1] 16.7987242  0.2412990 18.2172477  0.5126168
```
