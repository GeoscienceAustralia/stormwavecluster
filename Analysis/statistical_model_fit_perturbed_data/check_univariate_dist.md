# Sensitivity of the fitted univariate distributions to random perturbations of the data
---------------------------------------------------------------------------------------

The code below investigates how fitted model parameters and AIC-based choice of
copula can change due to data perturbation.

**Summary**: For the most part the changes to the model parameters are
small (often < 1%), as compared with fitting the model to the un-perturbed
data. However, we see that the maximum likelihood fit of the storm duration
data is sensitive to data perturbation.

Further, in some cases we find that the automatically chosen copula family differs
randomly due to the perturbation. This reflects that multiple copula may fit
the data quite well, and so have similar AIC values (close enough that the
perturbation can change the optimal family). With the exception of the duration
variable, there is usually one dominant family, with a few cases having other
families. 

However, for duration, the automatically chosen copula varies more
substantially, with both the Frank and Gaussian occurring often. This further
suggests that our fit to the duration data is affected by the data
discretization, to a greater extent than the other variables.

In all cases, the most common family derived from the perturbed data is the
same as the family selected for the original data.

**Read the model results, and do some basic checks**

```r
source('get_Rimage_data_univariate_distributions.R', local=TRUE)

# Read the summary statistics -- saved earlier for speed
store_var_list = readRDS('univariate_runs_summary_statistics.RDS')

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
            # Handle errors gracefully
            if( (length(num) != length(denom)) || any(is.na(num)) || 
                (class(num) != class(denom))){
                num = denom*NA
            }
            return(num/abs(denom))
        }
    )

    variable_differences = t(variable_differences)
    print(summary(variable_differences))
    return(invisible())
}
```



**Hsig model checks**

The hsig model is not very effected by data perturbation, although visually a
slight thickening of the upper credible interval is noticable for high return
periods.

```r
# Make a plot
plot_ylim = range(c(original_var_list$hsig_mixture_fit_bayes_lower_q,
    original_var_list$hsig_mixture_fit_bayes_upper_q))

xs = original_var_list$hsig_mixture_fit_bayes_rates
plot(xs, original_var_list$hsig_mixture_fit_ml, 
    t='l', ylim=plot_ylim, main='Hsig fit (all overplotted)', xlim=c(max(xs), min(xs)),
    log='xy', xlab='Event rate (number per year)', ylab='m')
points(xs, original_var_list$hsig_mixture_fit_bayes_median_q, col='orange', t='l')
points(xs, original_var_list$hsig_mixture_fit_bayes_upper_q, col='red', t='l')
points(xs, original_var_list$hsig_mixture_fit_bayes_lower_q, col='red', t='l')
for(i in 1:length(store_var_list)){
    xs = store_var_list[[i]]$hsig_mixture_fit_bayes_rates
    points(xs, store_var_list[[i]]$hsig_mixture_fit_ml, t='l')
    points(xs, store_var_list[[i]]$hsig_mixture_fit_bayes_median_q, col='orange', t='l')
    points(xs, store_var_list[[i]]$hsig_mixture_fit_bayes_upper_q, col='red', t='l')
    points(xs, store_var_list[[i]]$hsig_mixture_fit_bayes_lower_q, col='red', t='l')
}
points(xs, original_var_list$hsig_mixture_fit_ml, col='blue', t='p', pch=1, lty='dashed')
legend('topleft', 
    c('Maximum Likelihood', 'Bayesian Median', 'Bayesian 97.5%', 'Bayesian 2.5%', 'Unperturbed data ML'),
    col=c('black', 'orange', 'red', 'red', 'blue'), lwd=c(1,1,1,1,NA), pch=c(NA, NA, NA, NA, 1))
```

![plot of chunk hsig](figure/hsig-1.png)

```r
# Check how the hsig_mixture_fit parameters vary due to jittering
#
#
print(original_var_list$hsig_mixture_fit_par)
```

```
## [1]  0.8500381  1.0100971  1.2916491 -0.2194104
```

```r
#
# Errors typically O(1/1000), with extrema of about 1%
relative_error_summary('hsig_mixture_fit_par')
```

```
##        V1                   V2                   V3            
##  Min.   :-0.0240411   Min.   :-0.0538487   Min.   :-0.0188498  
##  1st Qu.:-0.0002164   1st Qu.:-0.0234349   1st Qu.:-0.0123787  
##  Median : 0.0036912   Median :-0.0038379   Median : 0.0007689  
##  Mean   : 0.0114793   Mean   :-0.0112137   Mean   : 0.0117865  
##  3rd Qu.: 0.0244210   3rd Qu.: 0.0006029   3rd Qu.: 0.0285843  
##  Max.   : 0.0529665   Max.   : 0.0253095   Max.   : 0.1109334  
##        V4           
##  Min.   :-0.022255  
##  1st Qu.: 0.001194  
##  Median : 0.007209  
##  Mean   : 0.006515  
##  3rd Qu.: 0.012237  
##  Max.   : 0.042043
```

```r
#
# Errors typically O(1/10000)
relative_error_summary('hsig_aep100_quantiles')
```

```
##       2.5%                 50%                 97.5%           
##  Min.   :-0.0010682   Min.   :-0.0008690   Min.   :-0.0041211  
##  1st Qu.:-0.0006229   1st Qu.:-0.0001589   1st Qu.:-0.0008875  
##  Median :-0.0003631   Median : 0.0003062   Median : 0.0001698  
##  Mean   :-0.0003792   Mean   : 0.0003829   Mean   : 0.0011156  
##  3rd Qu.:-0.0001472   3rd Qu.: 0.0007939   3rd Qu.: 0.0030575  
##  Max.   : 0.0003287   Max.   : 0.0023992   Max.   : 0.0111256
```

```r
# Small errors in optimal season phi [units of years -- typical value O(1 day)]
hsig_season_phi_err = unlist(sapply(store_var_list, 
    f<-function(x) x$hsig_season_phi%%1 - original_var_list$hsig_season_phi%%1))
print(summary(hsig_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0195400  0.0002563  0.0002773  0.0008247  0.0022280  0.0129100
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
##                             95                              2 
## Rotated Tawn type 2 90 degrees 
##                              3
```

**Duration model checks**

Duration is most significantly affected by its hourly discretization, among all
the cases reviewed here. The maximum likelihood fit seems to vary between two
main branches (one of which includes the model for the unperturbed data). The
Bayesian median is more stable. The variation in the credible intervals is also
visually evident (they appear thick when overplotted).

This can also be seen in the variations in the maximum likelihood model
parameters. Further, the copula type is less stable for duration, than for
other variables.

```r
# Make a plot
plot_ylim = range(c(original_var_list$duration_mixture_fit_bayes_lower_q,
    original_var_list$duration_mixture_fit_bayes_upper_q))

xs = original_var_list$duration_mixture_fit_bayes_rates
plot(xs, original_var_list$duration_mixture_fit_ml, 
    t='l', ylim=plot_ylim, main='Duration fit (all overplotted)', xlim=c(max(xs), min(xs)),
    log='xy', xlab='Event rate (number per year)', ylab='m')
points(xs, original_var_list$duration_mixture_fit_bayes_median_q, col='orange', t='l')
points(xs, original_var_list$duration_mixture_fit_bayes_upper_q, col='red', t='l')
points(xs, original_var_list$duration_mixture_fit_bayes_lower_q, col='red', t='l')
for(i in 1:length(store_var_list)){
    xs = store_var_list[[i]]$duration_mixture_fit_bayes_rates
    points(xs, store_var_list[[i]]$duration_mixture_fit_ml, t='l')
    points(xs, store_var_list[[i]]$duration_mixture_fit_bayes_median_q, col='orange', t='l')
    points(xs, store_var_list[[i]]$duration_mixture_fit_bayes_upper_q, col='red', t='l')
    points(xs, store_var_list[[i]]$duration_mixture_fit_bayes_lower_q, col='red', t='l')
}
points(xs, original_var_list$duration_mixture_fit_ml, col='blue', pch=1)
legend('topleft', 
    c('Maximum Likelihood', 'Bayesian Median', 'Bayesian 97.5%', 'Bayesian 2.5%', 'Unperturbed data ML'),
    col=c('black', 'orange', 'red', 'red', 'blue'), lwd=c(1,1,1,1, NA), pch=c(NA, NA, NA, NA, 1))
```

![plot of chunk durationplot](figure/durationplot-1.png)

```r
# Check how the duration_mixture_fit parameters vary due to jittering
#
print(original_var_list$duration_mixture_fit_par)
```

```
## [1]  0.7876761 31.8814876 51.3829002 -0.1393468
```

```r
#
# Errors 
relative_error_summary('duration_mixture_fit_par')
```

```
##        V1                V2               V3                 V4         
##  Min.   :-0.3117   Min.   :0.1656   Min.   :-0.87701   Min.   :-0.3009  
##  1st Qu.:-0.2591   1st Qu.:0.2007   1st Qu.:-0.83460   1st Qu.:-0.2576  
##  Median :-0.2324   Median :0.4990   Median :-0.77328   Median : 0.6176  
##  Mean   :-0.2232   Mean   :0.4811   Mean   :-0.55465   Mean   : 0.3446  
##  3rd Qu.:-0.1692   3rd Qu.:0.6840   3rd Qu.:-0.10041   3rd Qu.: 0.7227  
##  Max.   :-0.1468   Max.   :1.0393   Max.   :-0.08175   Max.   : 0.8161  
##  NA's   :2         NA's   :2        NA's   :2          NA's   :2
```

```r
#
# Errors in 1/100 AEP are small
relative_error_summary('duration_aep100_quantiles')
```

```
##       2.5%               50%                97.5%          
##  Min.   :0.004077   Min.   :-0.002154   Min.   :-0.048331  
##  1st Qu.:0.007401   1st Qu.: 0.006369   1st Qu.:-0.035300  
##  Median :0.008840   Median : 0.014390   Median :-0.028039  
##  Mean   :0.009561   Mean   : 0.016687   Mean   :-0.026744  
##  3rd Qu.:0.010842   3rd Qu.: 0.022366   3rd Qu.:-0.021420  
##  Max.   :0.023378   Max.   : 0.061010   Max.   : 0.008759  
##  NA's   :2          NA's   :2           NA's   :2
```

```r
# Small errors in optimal season phi [units of years -- typical value O(1 day)]
duration_season_phi_err = unlist(sapply(store_var_list, 
    f<-function(x) x$duration_season_phi%%1 - original_var_list$duration_season_phi%%1))
print(summary(duration_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -7.513e-03 -4.581e-03 -1.089e-03 -1.131e-03 -5.759e-05  7.281e-03
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
duration_copula_type = unlist(sapply(store_var_list, f<-function(x) x$duration_season_copula$familyname))
print(table(duration_copula_type))
```

```
## duration_copula_type
##                  Frank               Gaussian Rotated BB1 90 degrees 
##                     50                     47                      1
```

**Tidal residual model checks**

No perturbation was applied to the tidal residual data, because there were no
ties. Therefore, any variation in the model fits *only* reflects the use of a
finite number of MCMC samples to characterise the Bayesian model (note only one
chain was used). The figure suggests that the MCMC related variations are
small, consistent with our earlier analysis. 

```r
# Make a plot
plot_ylim = range(c(original_var_list$tideResid_mixture_fit_bayes_lower_q,
    original_var_list$tideResid_mixture_fit_bayes_upper_q))

xs = original_var_list$tideResid_mixture_fit_bayes_rates
plot(xs, original_var_list$tideResid_mixture_fit_ml, 
    t='l', ylim=plot_ylim, main='Tidal residual fit (all overplotted)', xlim=c(max(xs), min(xs)),
    log='xy', xlab='Event rate (number per year)', ylab='m')
points(xs, original_var_list$tideResid_mixture_fit_bayes_median_q, col='orange', t='l')
points(xs, original_var_list$tideResid_mixture_fit_bayes_upper_q, col='red', t='l')
points(xs, original_var_list$tideResid_mixture_fit_bayes_lower_q, col='red', t='l')
for(i in 1:length(store_var_list)){
    xs = store_var_list[[i]]$tideResid_mixture_fit_bayes_rates
    points(xs, store_var_list[[i]]$tideResid_mixture_fit_ml, t='l')
    points(xs, store_var_list[[i]]$tideResid_mixture_fit_bayes_median_q, col='orange', t='l')
    points(xs, store_var_list[[i]]$tideResid_mixture_fit_bayes_upper_q, col='red', t='l')
    points(xs, store_var_list[[i]]$tideResid_mixture_fit_bayes_lower_q, col='red', t='l')
}
points(xs, original_var_list$tideResid_mixture_fit_ml, col='blue', t='p', pch=1)
legend('topleft', 
    c('Maximum Likelihood', 'Bayesian Median', 'Bayesian 97.5%', 'Bayesian 2.5%', 'Unperturbed data ML'),
    col=c('black', 'orange', 'red', 'red', 'blue'), lwd=c(1,1,1,1,NA), pch=c(NA, NA, NA, NA, 1))
```

![plot of chunk tidal_residualA](figure/tidal_residualA-1.png)

```r
#
# These should not show error, because we didn't jitter tideResid.
# However, there can be small errors due to MCMC (since it uses random
# numbers)
#

# Check how the tideResid_mixture_fit parameters vary due to jittering
#
print(original_var_list$tideResid_mixture_fit_par)
```

```
## [1]  0.1146811  0.1135725  0.1855709 -0.1256922
```

```r
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
##  Min.   :-1.715e-03   Min.   :-1.502e-03   Min.   :-0.0120964  
##  1st Qu.:-3.260e-04   1st Qu.:-4.727e-04   1st Qu.:-0.0059022  
##  Median :-2.252e-05   Median :-2.061e-04   Median :-0.0031346  
##  Mean   : 1.172e-05   Mean   :-2.180e-04   Mean   :-0.0028974  
##  3rd Qu.: 3.554e-04   3rd Qu.: 5.136e-05   3rd Qu.: 0.0006905  
##  Max.   : 1.328e-03   Max.   : 8.716e-04   Max.   : 0.0079544
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


**Steepness model checks**

Steepness is modelled non-parametrically. The non-parametric model shows only
slight variation due to data perturbation. 

```r
#
# We use a non-parametric distribution, so here just check the copula and
# seasonal phase
#
plot(seq(0,1,len=200), original_var_list$steepness_quantiles_seq, t='l',
    main='Steepness', xlab='Non-exceedance probability for one event', ylab='Quantile')
for(i in 1:length(store_var_list)){
    points(seq(0,1,len=200), store_var_list[[i]]$steepness_quantiles_seq, t='l')
}
points(seq(0,1,len=200), original_var_list$steepness_quantiles_seq, pch=1, col='red')

legend('topleft', c('Non-parametric model', 'Unperturbed data model'), col=c('black', 'red'), 
    lwd=c(1,NA), pch=c(NA,1))
```

![plot of chunk steepnessModel](figure/steepnessModel-1.png)

```r
# Errors O(1 day)
steepness_season_phi_err = sapply(store_var_list, 
    f<-function(x) (x$steepness_season_phi%%1 - original_var_list$steepness_season_phi%%1))
print(summary(steepness_season_phi_err))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## -0.003998  0.001002  0.001039  0.001157  0.001393  0.003723
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
##                           2                          98
```

**Direction model checks**

Direction is modelled non-parametrically. The non-parametric model shows only
slight variation due to data perturbation. 

```r
#
# We use a non-parametric distribution, so here just check the copula (soiA +
# seasonal) and seasonal phase
#
plot(seq(0,1,len=200), original_var_list$dir_quantiles_seq, t='l',
    main='Direction', xlab='Non-exceedance probability for one event', ylab='Quantile')
for(i in 1:length(store_var_list)){
    points(seq(0,1,len=200), store_var_list[[i]]$dir_quantiles_seq, t='l')
}
points(seq(0,1,len=200), original_var_list$dir_quantiles_seq, pch=1, col='red')
legend('topleft', c('Non-parametric model', 'Unperturbed data model'), col=c('black', 'red'), 
    lwd=c(1,NA), pch=c(NA,1))
```

![plot of chunk directionModel](figure/directionModel-1.png)

```r
# Error O(1 day) or less
dir_season_phi_err = sapply(store_var_list, 
    f<-function(x) (x$dir_season_phi%%1 - original_var_list$dir_season_phi%%1))
print(summary(dir_season_phi_err))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## -0.005153  0.005390  0.006489  0.005374  0.006794  0.016810
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
##                         97                          3
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
