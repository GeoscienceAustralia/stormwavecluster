# Sensitivity of the fitted univariate distributions to random perturbations of the data
---------------------------------------------------------------------------------------

The code below investigates how fitted model parameters and AIC-based choice of
copula can change due to data perturbation.

**Read the model results, and do some basic checks**

```r
source('get_Rimage_data_univariate_distributions.R', local=TRUE)

# Read the summary statistics -- saved earlier for speed
store_var_list = readRDS('univariate_runs_summary_statistics.RDS')

# Read the original fit (based on un-perturbed data)
original_var_list = get_Rimage_data(
    '../statistical_model_fit/Rimages/session_univariate_distributions_FALSE_0.Rdata')
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
## 
## Attaching package: 'SparseM'
```

```
## The following object is masked from 'package:base':
## 
##     backsolve
```

```r
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
perturbed_summary<-function(variable_name){

    variable_vals = sapply(store_var_list,
        f<-function(x){
            num = x[[variable_name]] 
            denom = original_var_list[[variable_name]]
            # Handle errors gracefully
            if( (length(num) != length(denom)) || any(is.na(num)) || 
                (class(num) != class(denom))){
                num = denom*NA
            }
            return(num)
        }
    )

    variable_vals = t(variable_vals)
    print(summary(variable_vals))
    return(invisible())

}

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


```r
# Make a plot
plot_ylim = range(c(original_var_list$hsig_mixture_fit_bayes_lower_q,
    original_var_list$hsig_mixture_fit_bayes_upper_q))

xs = original_var_list$hsig_mixture_fit_bayes_rates
plot(xs, original_var_list$hsig_mixture_fit_ml, 
    t='l', ylim=plot_ylim, main='Hsig fit (all overplotted)', 
    xlim=c(max(xs), min(xs)),
    log='xy', xlab='Event rate (number per year)', ylab='m')
points(xs, original_var_list$hsig_mixture_fit_bayes_median_q, col='orange', t='l', lty='dotted')
points(xs, original_var_list$hsig_mixture_fit_bayes_upper_q, col='red', t='l', lty='dotted')
points(xs, original_var_list$hsig_mixture_fit_bayes_lower_q, col='red', t='l', lty='dotted')
for(i in 1:length(store_var_list)){
    xs = store_var_list[[i]]$hsig_mixture_fit_bayes_rates
    points(xs, store_var_list[[i]]$hsig_mixture_fit_ml, t='l')
    points(xs, store_var_list[[i]]$hsig_mixture_fit_bayes_median_q, col='orange', t='l')
    points(xs, store_var_list[[i]]$hsig_mixture_fit_bayes_upper_q, col='red', t='l')
    points(xs, store_var_list[[i]]$hsig_mixture_fit_bayes_lower_q, col='red', t='l')
}
points(xs, original_var_list$hsig_mixture_fit_ml, col='blue', t='p', pch=1, lty='dashed')
grid()
legend('topleft', 
    c('Maximum Likelihood', 'Bayesian Median', 'Bayesian 97.5%', 'Bayesian 2.5%', 
        'Unperturbed data ML'),
    col=c('black', 'orange', 'red', 'red', 'blue'), lwd=c(1,1,1,1,NA), 
    pch=c(NA, NA, NA, NA, 1))
```

![plot of chunk hsig](figure/hsig-1.png)

```r
# Check how the hsig_mixture_fit parameters vary due to jittering
#
#
print(original_var_list$hsig_mixture_fit_par)
```

```
## [1]  0.8424030  1.0204286  1.2726887 -0.2198767
```

```r
#
perturbed_summary('hsig_mixture_fit_par')
```

```
##        V1              V2                V3                V4          
##  Min.   :1.007   Min.   :0.01189   Min.   :0.01693   Min.   :-0.20142  
##  1st Qu.:1.775   1st Qu.:0.10557   1st Qu.:0.03229   1st Qu.:-0.04391  
##  Median :2.054   Median :0.15239   Median :0.04098   Median :-0.03941  
##  Mean   :2.312   Mean   :0.16789   Mean   :0.05445   Mean   :-0.04103  
##  3rd Qu.:2.407   3rd Qu.:0.23120   3rd Qu.:0.04696   3rd Qu.:-0.03530  
##  Max.   :7.013   Max.   :0.86687   Max.   :1.49498   Max.   :-0.02568  
##  NA's   :1       NA's   :1         NA's   :1         NA's   :1
```

```r
relative_error_summary('hsig_mixture_fit_par')
```

```
##        V1               V2                V3                V4         
##  Min.   :0.1959   Min.   :-0.9883   Min.   :-0.9867   Min.   :0.08393  
##  1st Qu.:1.1076   1st Qu.:-0.8965   1st Qu.:-0.9746   1st Qu.:0.80028  
##  Median :1.4388   Median :-0.8507   Median :-0.9678   Median :0.82074  
##  Mean   :1.7449   Mean   :-0.8355   Mean   :-0.9572   Mean   :0.81339  
##  3rd Qu.:1.8578   3rd Qu.:-0.7734   3rd Qu.:-0.9631   3rd Qu.:0.83947  
##  Max.   :7.3247   Max.   :-0.1505   Max.   : 0.1747   Max.   :0.88323  
##  NA's   :1        NA's   :1         NA's   :1         NA's   :1
```

```r
#
print(original_var_list$hsig_aep100_quantiles)
```

```
##     2.5%      50%    97.5% 
## 7.068316 7.542978 8.877614
```

```r
perturbed_summary('hsig_aep100_quantiles')
```

```
##       2.5%            50%            97.5%       
##  Min.   :7.025   Min.   :7.593   Min.   : 9.268  
##  1st Qu.:7.056   1st Qu.:7.674   1st Qu.: 9.670  
##  Median :7.079   Median :7.741   Median : 9.936  
##  Mean   :7.083   Mean   :7.845   Mean   :10.093  
##  3rd Qu.:7.097   3rd Qu.:7.953   3rd Qu.:10.408  
##  Max.   :7.205   Max.   :9.188   Max.   :12.584  
##  NA's   :8       NA's   :8       NA's   :8
```

```r
# Errors typically O(1/10000)
relative_error_summary('hsig_aep100_quantiles')
```

```
##       2.5%                50%               97.5%        
##  Min.   :-0.006089   Min.   :0.006617   Min.   :0.04394  
##  1st Qu.:-0.001764   1st Qu.:0.017343   1st Qu.:0.08930  
##  Median : 0.001509   Median :0.026253   Median :0.11922  
##  Mean   : 0.002078   Mean   :0.040034   Mean   :0.13690  
##  3rd Qu.: 0.004059   3rd Qu.:0.054373   3rd Qu.:0.17244  
##  Max.   : 0.019346   Max.   :0.218112   Max.   :0.41752  
##  NA's   :8           NA's   :8          NA's   :8
```

```r
# Small errors in optimal season phi [units of years -- typical value O(1 day)]
hsig_season_phi_err = unlist(sapply(store_var_list, 
    f<-function(x) x$hsig_season_phi%%1 - original_var_list$hsig_season_phi%%1))
print(summary(hsig_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0244000 -0.0002318  0.0002733 -0.0019780  0.0009050  0.0129300
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
hsig_copula_type = sapply(store_var_list, f<-function(x) as.character(x$hsig_season_copula$familyname))
print(table(hsig_copula_type))
```

```
## Error in table(hsig_copula_type): all arguments must have the same length
```

**Duration model checks**


```r
# Make a plot
plot_ylim = range(c(original_var_list$duration_mixture_fit_bayes_lower_q,
    original_var_list$duration_mixture_fit_bayes_upper_q))

xs = original_var_list$duration_mixture_fit_bayes_rates
plot(xs, original_var_list$duration_mixture_fit_ml, 
    t='l', ylim=plot_ylim, main='Duration fit (all overplotted)', xlim=c(max(xs), min(xs)),
    log='xy', xlab='Event rate (number per year)', ylab='m')
points(xs, original_var_list$duration_mixture_fit_bayes_median_q, col='orange', t='l', lty='dotted')
points(xs, original_var_list$duration_mixture_fit_bayes_upper_q, col='red', t='l', lty='dotted')
points(xs, original_var_list$duration_mixture_fit_bayes_lower_q, col='red', t='l', lty='dotted')
for(i in 1:length(store_var_list)){
    xs = store_var_list[[i]]$duration_mixture_fit_bayes_rates
    points(xs, store_var_list[[i]]$duration_mixture_fit_ml, t='l')
    points(xs, store_var_list[[i]]$duration_mixture_fit_bayes_median_q, col='orange', t='l')
    points(xs, store_var_list[[i]]$duration_mixture_fit_bayes_upper_q, col='red', t='l')
    points(xs, store_var_list[[i]]$duration_mixture_fit_bayes_lower_q, col='red', t='l')
}
points(xs, original_var_list$duration_mixture_fit_ml, col='blue', pch=1)
```

```
## Error in xy.coords(x, y): 'x' and 'y' lengths differ
```

```r
grid()
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
## [1]  0.6812525 36.5813149 47.3420093 -0.1675303
```

```r
#
# Errors 
perturbed_summary('duration_mixture_fit_par')
```

```
##        V1            V2            V3            V4     
##  Min.   : NA   Min.   : NA   Min.   : NA   Min.   : NA  
##  1st Qu.: NA   1st Qu.: NA   1st Qu.: NA   1st Qu.: NA  
##  Median : NA   Median : NA   Median : NA   Median : NA  
##  Mean   :NaN   Mean   :NaN   Mean   :NaN   Mean   :NaN  
##  3rd Qu.: NA   3rd Qu.: NA   3rd Qu.: NA   3rd Qu.: NA  
##  Max.   : NA   Max.   : NA   Max.   : NA   Max.   : NA  
##  NA's   :100   NA's   :100   NA's   :100   NA's   :100
```

```r
relative_error_summary('duration_mixture_fit_par')
```

```
##        V1            V2            V3            V4     
##  Min.   : NA   Min.   : NA   Min.   : NA   Min.   : NA  
##  1st Qu.: NA   1st Qu.: NA   1st Qu.: NA   1st Qu.: NA  
##  Median : NA   Median : NA   Median : NA   Median : NA  
##  Mean   :NaN   Mean   :NaN   Mean   :NaN   Mean   :NaN  
##  3rd Qu.: NA   3rd Qu.: NA   3rd Qu.: NA   3rd Qu.: NA  
##  Max.   : NA   Max.   : NA   Max.   : NA   Max.   : NA  
##  NA's   :100   NA's   :100   NA's   :100   NA's   :100
```

```r
#
# Errors in 1/100 AEP are small
perturbed_summary('duration_aep100_quantiles')
```

```
##       2.5%          50%          97.5%    
##  Min.   : NA   Min.   : NA   Min.   : NA  
##  1st Qu.: NA   1st Qu.: NA   1st Qu.: NA  
##  Median : NA   Median : NA   Median : NA  
##  Mean   :NaN   Mean   :NaN   Mean   :NaN  
##  3rd Qu.: NA   3rd Qu.: NA   3rd Qu.: NA  
##  Max.   : NA   Max.   : NA   Max.   : NA  
##  NA's   :100   NA's   :100   NA's   :100
```

```r
relative_error_summary('duration_aep100_quantiles')
```

```
##       2.5%          50%          97.5%    
##  Min.   : NA   Min.   : NA   Min.   : NA  
##  1st Qu.: NA   1st Qu.: NA   1st Qu.: NA  
##  Median : NA   Median : NA   Median : NA  
##  Mean   :NaN   Mean   :NaN   Mean   :NaN  
##  3rd Qu.: NA   3rd Qu.: NA   3rd Qu.: NA  
##  Max.   : NA   Max.   : NA   Max.   : NA  
##  NA's   :100   NA's   :100   NA's   :100
```

```r
# Small errors in optimal season phi [units of years -- typical value O(1 day)]
duration_season_phi_err = unlist(sapply(store_var_list, 
    f<-function(x) x$duration_season_phi%%1 - original_var_list$duration_season_phi%%1))
print(summary(duration_season_phi_err))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 
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
## < table of extent 0 >
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
points(xs, original_var_list$tideResid_mixture_fit_bayes_median_q, col='orange', t='l', lty='dotted')
points(xs, original_var_list$tideResid_mixture_fit_bayes_upper_q, col='red', t='l', lty='dotted')
points(xs, original_var_list$tideResid_mixture_fit_bayes_lower_q, col='red', t='l', lty='dotted')
for(i in 1:length(store_var_list)){
    xs = store_var_list[[i]]$tideResid_mixture_fit_bayes_rates
    points(xs, store_var_list[[i]]$tideResid_mixture_fit_ml, t='l')
    points(xs, store_var_list[[i]]$tideResid_mixture_fit_bayes_median_q, col='orange', t='l')
    points(xs, store_var_list[[i]]$tideResid_mixture_fit_bayes_upper_q, col='red', t='l')
    points(xs, store_var_list[[i]]$tideResid_mixture_fit_bayes_lower_q, col='red', t='l')
}
points(xs, original_var_list$tideResid_mixture_fit_ml, col='blue', t='p', pch=1)
legend('topleft', 
    c('Maximum Likelihood', 'Bayesian Median', 'Bayesian 97.5%', 
        'Bayesian 2.5%', 'Unperturbed data ML'),
    col=c('black', 'orange', 'red', 'red', 'blue'), lwd=c(1,1,1,1,NA), 
    pch=c(NA, NA, NA, NA, 1))
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
##        V1                   V2                   V3           
##  Min.   :-3.397e-03   Min.   :-5.528e-03   Min.   :-0.022940  
##  1st Qu.: 2.745e-05   1st Qu.: 3.909e-05   1st Qu.: 0.002292  
##  Median : 9.438e-04   Median : 1.744e-03   Median : 0.013804  
##  Mean   : 1.340e-03   Mean   : 1.929e-03   Mean   : 0.017867  
##  3rd Qu.: 2.710e-03   3rd Qu.: 3.720e-03   3rd Qu.: 0.028752  
##  Max.   : 1.096e-02   Max.   : 1.581e-02   Max.   : 0.178794  
##        V4            
##  Min.   :-0.0618716  
##  1st Qu.:-0.0001238  
##  Median : 0.0277509  
##  Mean   : 0.0304730  
##  3rd Qu.: 0.0488003  
##  Max.   : 0.3404370
```

```r
#
# Errors typically < 1/1000 -- this is purely due to MCMC
relative_error_summary('tideResid_aep100_quantiles')
```

```
##       2.5%                 50%                 97.5%          
##  Min.   :-0.0039623   Min.   :-0.0085172   Min.   :-0.015499  
##  1st Qu.:-0.0016827   1st Qu.:-0.0029981   1st Qu.: 0.001251  
##  Median :-0.0001182   Median :-0.0009441   Median : 0.007571  
##  Mean   :-0.0000747   Mean   :-0.0005440   Mean   : 0.007729  
##  3rd Qu.: 0.0013958   3rd Qu.: 0.0017748   3rd Qu.: 0.013771  
##  Max.   : 0.0038761   Max.   : 0.0065260   Max.   : 0.032681
```

```r
# No errors
tideResid_season_phi_err = sapply(store_var_list, 
    f<-function(x) x$tideResid_season_phi%%1 - original_var_list$tideResid_season_phi%%1)
print(summary(tideResid_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -5.275e-03 -4.445e-03 -2.820e-03 -2.917e-03 -9.427e-04  7.757e-05
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
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0042970  0.0009921  0.0010250  0.0009423  0.0010670  0.0036640
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
##                    Gaussian Rotated Clayton 270 degrees 
##                           4                          96
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
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0051810  0.0003676  0.0057140  0.0045230  0.0067150  0.0168400
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
##                         96                          2 
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
