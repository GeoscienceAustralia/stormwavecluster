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
##     666     666     666     666     666     666
```

```r
# Number of events in original data
print(length(original_storm_timing_session$event_statistics[,1]))
```

```
## [1] 666
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
## theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+0 
##                                                                                                                       100
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
## [1] "Equation : theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+0"
## [1] " .... was selected by: 100 runs"
## [1] ""
## [1] ".... Summary of fitted parameters"
##        V1              V2               V3              V4        
##  Min.   :16.53   Min.   :0.2186   Min.   :17.57   Min.   :0.5221  
##  1st Qu.:16.54   1st Qu.:0.2192   1st Qu.:17.61   1st Qu.:0.5222  
##  Median :16.55   Median :0.2194   Median :17.62   Median :0.5238  
##  Mean   :16.55   Mean   :0.2195   Mean   :17.61   Mean   :0.5234  
##  3rd Qu.:16.55   3rd Qu.:0.2197   3rd Qu.:17.62   3rd Qu.:0.5239  
##  Max.   :16.56   Max.   :0.2227   Max.   :17.65   Max.   :0.5239  
## [1] ""
## [1] ".... Coefficient of variation of perturbed model parameters: "
## [1] 0.0003391353 0.0025008745 0.0007453631 0.0014414920
## [1] ""
## [1] ".... Parameters / approximate_standard_errors"
##        V1              V2              V3              V4       
##  Min.   :9.530   Min.   :1.695   Min.   :5.254   Min.   :19.05  
##  1st Qu.:9.532   1st Qu.:1.700   1st Qu.:5.263   1st Qu.:19.70  
##  Median :9.534   Median :1.702   Median :5.266   Median :22.03  
##  Mean   :9.534   Mean   :1.702   Mean   :5.266   Mean   :21.38  
##  3rd Qu.:9.536   3rd Qu.:1.704   3rd Qu.:5.268   3rd Qu.:22.13  
##  Max.   :9.540   Max.   :1.728   Max.   :5.276   Max.   :22.29
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
## [1] 16.5435353  0.2195488 17.6198233  0.5238601
```
