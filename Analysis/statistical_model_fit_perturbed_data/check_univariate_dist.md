Sensitivity of the fitted univariate distributions to random tie-breaking perturbations of the data
---------------------------------------------------------------------------------------------------

The code below investigates how fitted model parameters and AIC-based choice of
copula can change due to data perturbation.

In general, the changes to the model parameters are very small (often < 1%), as
compared with fitting the model to the un-perturbed data. They are not regarded
as significant.

In some cases we find that the automatically chosen copula family differs
randomly due to the perturbation. This reflects that multiple copula may fit
the data quite well, and so have similar AIC values (close enough that the
perturbation can change the optimal family). With the exception of the duration
variable, there is usually one dominant family, with a few cases having other
families. However, for duration, the automatically chosen copula varies
substantially between the Frank (large minority) and Gaussian (majority)
families. This simply reflects that both provide nearly equal fits to the data,
so random perturbations can lead to one or the other having best AIC.



```r
all_ud = Sys.glob('Rimages/session_univariate_distributions_TRUE_*.Rdata')

#
# The session images are too large load all at once.
# So get the important variables here. Do it in parallel (although possibly
# speed is mainly limited by disk reads?) 
#
get_Rimage_data<-function(all_ud_i){
    
    random_env = new.env()
    load(all_ud_i, envir=random_env)

    # Store key variables from random_env in a list, which we will later put in
    # store_var_list
    output = list()

    # Get break_ties_with_jitter
    output$break_ties_with_jitter = random_env$break_ties_with_jitter

    # Get the perturbed event statistics
    output$event_statistics = random_env$event_statistics

    # Hsig info -- parameters, and (Bayesian) quantiles of 1/100 AEP, and seasonal copula type, and seasonal phase
    output$hsig_mixture_fit_par = random_env$hsig_mixture_fit$fit_optim$par
    output$hsig_aep100_quantiles = quantile(random_env$hsig_mixture_fit$ari_100_chains[[1]], 
        p=c(0.025, 0.5, 0.975))
    output$hsig_season_phi = random_env$hsig_fit_conditional$season_phi_value$minimum
    output$hsig_season_copula = random_env$hsig_fit_conditional$var_season_copula

    # Duration info -- parameters, and (Bayesian) quantiles of 1/100 AEP, and seasonal copula type, and seasonal phase
    output$duration_mixture_fit_par = random_env$duration_mixture_fit$fit_optim$par
    output$duration_aep100_quantiles = quantile(random_env$duration_mixture_fit$ari_100_chains[[1]], 
        p=c(0.025, 0.5, 0.975))
    output$duration_season_phi = random_env$duration_fit_conditional$season_phi_value$minimum
    output$duration_season_copula = random_env$duration_fit_conditional$var_season_copula

    # TideResid info -- parameters, and (Bayesian) quantiles of 1/100 AEP, and seasonal copula type, and seasonal phase
    output$tideResid_mixture_fit_par = random_env$tideResid_mixture_fit$fit_optim$par
    output$tideResid_aep100_quantiles = quantile(random_env$tideResid_mixture_fit$ari_100_chains[[1]], 
        p=c(0.025, 0.5, 0.975))
    output$tideResid_season_phi = random_env$tideResid_fit_conditional$season_phi_value$minimum
    output$tideResid_season_copula = random_env$tideResid_fit_conditional$var_season_copula

    # steepness
    output$steepness_season_phi = random_env$steepness_fit_conditional$season_phi_value$minimum
    output$steepness_season_copula = random_env$steepness_fit_conditional$var_season_copula 

    # direction
    output$dir_season_phi = random_env$dir_fit_conditional$season_phi_value$minimum
    output$dir_season_copula = random_env$dir_fit_conditional$vargivensoiA_season_copula 
    output$dir_soiA_copula = random_env$dir_fit_conditional$var_soiA_copula

    # Store the result and move on
    #store_var_list[[all_ud[i]]] = output
    return(output)
}

# Read all images
store_var_list = lapply(as.list(all_ud), get_Rimage_data)

# Read the original fit (based on un-perturbed data)
original_var_list = get_Rimage_data(
    '../statistical_model_fit/Rimages/session_univariate_distributions_FALSE_0.Rdata')

# Check that all the perturbed data sessions do jittering
stopifnot(all(sapply(store_var_list, f<-function(x) x$break_ties_with_jitter)))
# Check the original fit does not do jittering
stopifnot(original_var_list$break_ties_with_jitter == FALSE)

# Check that the event_statistics is unique in every session [i.e. the perturbed
# sessions really do randomly perturb event_statistics. We perturbed hsig, duration,
# tp1, and dir [tideResid was already unique].
#
# To do the check, compute the column sums of all event statistics. They should
# all be unique
#
max_es_vals = sapply(store_var_list, 
    f<-function(x) colSums(x$event_statistics[,1:4], na.rm=TRUE))
stopifnot(length(unique(max_es_vals)) == length(max_es_vals))


#
# Useful function
#
relative_error_summary<-function(variable_name){
    variable_differences = sapply(store_var_list,
        f<-function(x){
            num = x[[variable_name]] - original_var_list[[variable_name]]
            denom = original_var_list[[variable_name]]
            return(num/denom)
        }
    )

    variable_differences = t(variable_differences)
    print(summary(variable_differences))
    return(invisible())
}

#
# HSIG MODEL CHECKS
#

# Check how the hsig_mixture_fit parameters vary due to jittering
#
# Errors typically O(1/1000), with extrema of about 1%
relative_error_summary('hsig_mixture_fit_par')
```

```
##        V1                  V2                  V3            
##  Min.   :-0.011693   Min.   :-0.016745   Min.   :-0.0024026  
##  1st Qu.:-0.001048   1st Qu.:-0.004978   1st Qu.:-0.0002168  
##  Median : 0.002034   Median :-0.002351   Median : 0.0006817  
##  Mean   : 0.001495   Mean   :-0.001768   Mean   : 0.0024109  
##  3rd Qu.: 0.004352   3rd Qu.: 0.001229   3rd Qu.: 0.0015281  
##  Max.   : 0.014942   Max.   : 0.013552   Max.   : 0.0158706  
##        V4           
##  Min.   :-0.014714  
##  1st Qu.:-0.004210  
##  Median :-0.001039  
##  Mean   :-0.000516  
##  3rd Qu.: 0.004086  
##  Max.   : 0.019716
```

```r
#
# Errors typically O(1/10000)
relative_error_summary('hsig_aep100_quantiles')
```

```
##       2.5%                 50%                 97.5%           
##  Min.   :-5.889e-04   Min.   :-0.0006708   Min.   :-0.0001268  
##  1st Qu.:-1.556e-04   1st Qu.:-0.0002324   1st Qu.: 0.0020165  
##  Median : 8.623e-05   Median : 0.0001026   Median : 0.0033010  
##  Mean   : 9.514e-05   Mean   : 0.0001258   Mean   : 0.0033217  
##  3rd Qu.: 3.219e-04   3rd Qu.: 0.0004196   3rd Qu.: 0.0043866  
##  Max.   : 1.074e-03   Max.   : 0.0020585   Max.   : 0.0112398
```

```r
# Small errors in optimal season phi [units of years -- typical value O(1 day)]
hsig_season_phi_err = sapply(store_var_list, 
    f<-function(x) x$hsig_season_phi%%1 - original_var_list$hsig_season_phi%%1)
print(summary(hsig_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0224500  0.0002547  0.0002783  0.0003373  0.0022280  0.0129100
```

```r
# Look at the automatically chosen copula family. Mostly Frank, with
# with occasional alternatives
print(original_var_list$hsig_season_copula$familyname)
```

```
## [1] "Frank"
```

```r
hsig_copula_type = sapply(store_var_list, f<-function(x) x$hsig_season_copula$familyname)
print(table(hsig_copula_type))
```

```
## hsig_copula_type
##                          Frank                       Gaussian 
##                             94                              3 
## Rotated Tawn type 2 90 degrees 
##                              3
```

```r
#
# DURATION MODEL CHECKS
#

# Check how the duration_mixture_fit parameters vary due to jittering
#
# Errors typically a few percent
relative_error_summary('duration_mixture_fit_par')
```

```
##        V1                  V2                  V3           
##  Min.   :-0.025404   Min.   :-0.006086   Min.   :-0.028341  
##  1st Qu.:-0.013573   1st Qu.: 0.007539   1st Qu.:-0.018860  
##  Median :-0.010814   Median : 0.011989   Median :-0.016317  
##  Mean   :-0.010329   Mean   : 0.011674   Mean   :-0.012794  
##  3rd Qu.:-0.006789   3rd Qu.: 0.015427   3rd Qu.:-0.003627  
##  Max.   : 0.006985   Max.   : 0.027493   Max.   : 0.010442  
##        V4           
##  Min.   :-0.034481  
##  1st Qu.: 0.005988  
##  Median : 0.014010  
##  Mean   : 0.013806  
##  3rd Qu.: 0.022406  
##  Max.   : 0.055536
```

```r
#
# Errors typically a few parts per thousand 
relative_error_summary('duration_aep100_quantiles')
```

```
##       2.5%                 50%                 97.5%           
##  Min.   :-0.0029876   Min.   :-0.0046000   Min.   :-0.0179435  
##  1st Qu.:-0.0018137   1st Qu.:-0.0024460   1st Qu.:-0.0063712  
##  Median :-0.0009194   Median :-0.0014128   Median :-0.0035568  
##  Mean   :-0.0009550   Mean   :-0.0014550   Mean   :-0.0033549  
##  3rd Qu.:-0.0001476   3rd Qu.:-0.0003593   3rd Qu.:-0.0008886  
##  Max.   : 0.0011102   Max.   : 0.0018253   Max.   : 0.0083458
```

```r
# Small errors in optimal season phi [units of years -- typical value O(1 day)]
duration_season_phi_err = sapply(store_var_list, 
    f<-function(x) x$duration_season_phi%%1 - original_var_list$duration_season_phi%%1)
print(summary(duration_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -7.513e-03 -4.577e-03 -1.105e-03 -1.139e-03 -8.049e-05  2.207e-02
```

```r
# Look at the automatically chosen copula family. Not very stable to
# perturbations. Both Frank and Gaussian come up, but Gaussian is most common. 
print(original_var_list$duration_season_copula$familyname)
```

```
## [1] "Gaussian"
```

```r
duration_copula_type = sapply(store_var_list, f<-function(x) x$duration_season_copula$familyname)
print(table(duration_copula_type))
```

```
## duration_copula_type
##                  Frank               Gaussian Rotated BB1 90 degrees 
##                     42                     57                      1
```

```r
#
# TIDERESID MODEL CHECKS
#
# These should not show error, because we didn't jitter tideResid.
# However, there can be small errors due to MCMC (since it uses random
# numbers)
#

# Check how the tideResid_mixture_fit parameters vary due to jittering
#
# No errors
relative_error_summary('tideResid_mixture_fit_par')
```

```
##        V1          V2          V3          V4   
##  Min.   :0   Min.   :0   Min.   :0   Min.   :0  
##  1st Qu.:0   1st Qu.:0   1st Qu.:0   1st Qu.:0  
##  Median :0   Median :0   Median :0   Median :0  
##  Mean   :0   Mean   :0   Mean   :0   Mean   :0  
##  3rd Qu.:0   3rd Qu.:0   3rd Qu.:0   3rd Qu.:0  
##  Max.   :0   Max.   :0   Max.   :0   Max.   :0
```

```r
#
# Errors typically < 1/1000 -- this is purely due to MCMC
relative_error_summary('tideResid_aep100_quantiles')
```

```
##       2.5%                 50%                 97.5%           
##  Min.   :-8.367e-05   Min.   :-0.0005413   Min.   :-0.0004299  
##  1st Qu.:-8.367e-05   1st Qu.:-0.0005413   1st Qu.: 0.0007583  
##  Median :-8.367e-05   Median :-0.0005413   Median : 0.0007583  
##  Mean   :-7.654e-05   Mean   :-0.0005320   Mean   : 0.0009042  
##  3rd Qu.:-8.367e-05   3rd Qu.:-0.0005413   3rd Qu.: 0.0007583  
##  Max.   : 2.285e-04   Max.   : 0.0001433   Max.   : 0.0086477
```

```r
# No errors
tideResid_season_phi_err = sapply(store_var_list, 
    f<-function(x) x$tideResid_season_phi%%1 - original_var_list$tideResid_season_phi%%1)
print(summary(tideResid_season_phi_err))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       0       0       0       0       0       0
```

```r
# Look at the automatically chosen copula family. 
# Always Gaussian
print(original_var_list$tideResid_season_copula$familyname)
```

```
## [1] "Gaussian"
```

```r
tideResid_copula_type = sapply(store_var_list, f<-function(x) x$tideResid_season_copula$familyname)
print(table(tideResid_copula_type))
```

```
## tideResid_copula_type
## Gaussian 
##      100
```

```r
#
# Steepness model checks
#
# We use a non-parametric distribution, so here just check the copula and
# seasonal phase
#

# Errors O(1 day)
steepness_season_phi_err = sapply(store_var_list, 
    f<-function(x) (x$steepness_season_phi%%1 - original_var_list$steepness_season_phi%%1))
print(summary(steepness_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0039540  0.0009908  0.0010330  0.0011030  0.0014690  0.0037230
```

```r
# Mostly Clayton copula, like in the original data
print(original_var_list$steepness_season_copula$familyname)
```

```
## [1] "Rotated Clayton 270 degrees"
```

```r
steepness_copula_type = sapply(store_var_list, f<-function(x) x$steepness_season_copula$familyname)
print(table(steepness_copula_type))
```

```
## steepness_copula_type
##     Rotated BB1 270 degrees Rotated Clayton 270 degrees 
##                           7                          93
```

```r
#
# Direction model checks
#
# We use a non-parametric distribution, so here just check the copula (soiA +
# seasonal) and seasonal phase
#

# Error O(1 day) or less
dir_season_phi_err = sapply(store_var_list, 
    f<-function(x) (x$dir_season_phi%%1 - original_var_list$dir_season_phi%%1))
print(summary(dir_season_phi_err))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## -0.005543  0.005377  0.006421  0.005470  0.006793  0.019980
```

```r
# Mostly Frank, like with the original data 
print(original_var_list$dir_season_copula$familyname)
```

```
## [1] "Frank"
```

```r
dir_copula_type = sapply(store_var_list, f<-function(x) x$dir_season_copula$familyname)
print(table(dir_copula_type))
```

```
## dir_copula_type
##                      Frank Rotated Clayton 90 degrees 
##                         93                          5 
## Rotated Gumbel 270 degrees 
##                          2
```

```r
# Always Frank, like with the original data
print(original_var_list$dir_soiA_copula$familyname)
```

```
## [1] "Frank"
```

```r
dir_copula_type = sapply(store_var_list, f<-function(x) x$dir_soiA_copula$familyname)
print(table(dir_copula_type))
```

```
## dir_copula_type
## Frank 
##   100
```
