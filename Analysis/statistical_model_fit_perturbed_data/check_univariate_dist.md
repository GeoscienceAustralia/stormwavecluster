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
# Original parameters
print(original_var_list$hsig_mixture_fit_par)
```

```
## [1]  0.8421592  1.0160924  1.2705059 -0.2143640
```

```r
# Perturbed parameter summary
perturbed_summary('hsig_mixture_fit_par')
```

```
##        V1               V2              V3              V4         
##  Min.   :0.8383   Min.   :1.014   Min.   :1.269   Min.   :-0.2158  
##  1st Qu.:0.8413   1st Qu.:1.016   1st Qu.:1.270   1st Qu.:-0.2147  
##  Median :0.8418   Median :1.017   Median :1.270   Median :-0.2145  
##  Mean   :0.8418   Mean   :1.017   Mean   :1.270   Mean   :-0.2145  
##  3rd Qu.:0.8423   3rd Qu.:1.017   3rd Qu.:1.271   3rd Qu.:-0.2143  
##  Max.   :0.8440   Max.   :1.021   Max.   :1.271   Max.   :-0.2136
```

```r
# Relative error in parameters summary
relative_error_summary('hsig_mixture_fit_par')
```

```
##        V1                   V2                   V3            
##  Min.   :-0.0046218   Min.   :-0.0024134   Min.   :-1.099e-03  
##  1st Qu.:-0.0009976   1st Qu.:-0.0001996   1st Qu.:-2.303e-04  
##  Median :-0.0003884   Median : 0.0004163   Median :-9.223e-05  
##  Mean   :-0.0004498   Mean   : 0.0004818   Mean   :-1.063e-04  
##  3rd Qu.: 0.0001810   3rd Qu.: 0.0011059   3rd Qu.: 3.896e-05  
##  Max.   : 0.0021966   Max.   : 0.0051084   Max.   : 4.971e-04  
##        V4            
##  Min.   :-0.0066776  
##  1st Qu.:-0.0014171  
##  Median :-0.0005245  
##  Mean   :-0.0006189  
##  3rd Qu.: 0.0002528  
##  Max.   : 0.0034474
```

```r
#

# Original 1/100 quantile
print(original_var_list$hsig_aep100_quantiles)
```

```
##     2.5%      50%    97.5% 
## 7.071127 7.563662 8.957316
```

```r
# Summary of perturbed 1/100 quantile
perturbed_summary('hsig_aep100_quantiles')
```

```
##       2.5%            50%            97.5%      
##  Min.   :7.069   Min.   :7.558   Min.   :8.927  
##  1st Qu.:7.071   1st Qu.:7.561   1st Qu.:8.944  
##  Median :7.072   Median :7.562   Median :8.949  
##  Mean   :7.072   Mean   :7.562   Mean   :8.949  
##  3rd Qu.:7.072   3rd Qu.:7.563   3rd Qu.:8.955  
##  Max.   :7.074   Max.   :7.565   Max.   :8.972
```

```r
# Summary of relative error
relative_error_summary('hsig_aep100_quantiles')
```

```
##       2.5%                 50%                 97.5%           
##  Min.   :-2.778e-04   Min.   :-7.229e-04   Min.   :-0.0033312  
##  1st Qu.:-3.671e-05   1st Qu.:-3.205e-04   1st Qu.:-0.0014879  
##  Median : 6.039e-05   Median :-1.764e-04   Median :-0.0009047  
##  Mean   : 6.036e-05   Mean   :-2.132e-04   Mean   :-0.0008767  
##  3rd Qu.: 1.821e-04   3rd Qu.:-9.373e-05   3rd Qu.:-0.0002986  
##  Max.   : 3.475e-04   Max.   : 1.204e-04   Max.   : 0.0016582
```

```r
# Variability in optimal season phi [units of years -- typical value O(1 day)]
hsig_season_phi_err = unlist(sapply(store_var_list, 
    f<-function(x) x$hsig_season_phi%%1 - original_var_list$hsig_season_phi%%1))
print(summary(hsig_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -3.312e-04  2.948e-05  1.580e-04  2.434e-03  2.546e-03  1.352e-02
```

```r
# Look at the automatically chosen copula family -- possibly more than one is identified
print(original_var_list$hsig_season_copula$familyname)
```

```
## [1] "Gaussian"
```

```r
hsig_copula_type = sapply(store_var_list, f<-function(x) as.character(x$hsig_season_copula$familyname))
print(table(hsig_copula_type))
```

```
## hsig_copula_type
##                     Frank                  Gaussian 
##                        15                        84 
## Rotated Gumbel 90 degrees 
##                         1
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
## [1]  0.6938297 35.9454043 48.0295307 -0.1633729
```

```r
#
# Perturbed parameters
perturbed_summary('duration_mixture_fit_par')
```

```
##        V1               V2              V3               V4          
##  Min.   :0.5644   Min.   :36.46   Min.   : 7.232   Min.   :-0.18032  
##  1st Qu.:0.6025   1st Qu.:37.30   1st Qu.: 8.082   1st Qu.:-0.17131  
##  Median :0.6202   Median :46.10   Median :11.956   Median :-0.05013  
##  Mean   :0.6314   Mean   :44.88   Mean   :24.981   Mean   :-0.09260  
##  3rd Qu.:0.6715   3rd Qu.:51.85   3rd Qu.:47.445   3rd Qu.:-0.03467  
##  Max.   :0.6853   Max.   :58.73   Max.   :48.399   Max.   :-0.01620
```

```r
relative_error_summary('duration_mixture_fit_par')
```

```
##        V1                 V2                V3           
##  Min.   :-0.18650   Min.   :0.01425   Min.   :-0.849430  
##  1st Qu.:-0.13160   1st Qu.:0.03757   1st Qu.:-0.831727  
##  Median :-0.10607   Median :0.28244   Median :-0.751063  
##  Mean   :-0.08993   Mean   :0.24862   Mean   :-0.479880  
##  3rd Qu.:-0.03215   3rd Qu.:0.44257   3rd Qu.:-0.012163  
##  Max.   :-0.01233   Max.   :0.63373   Max.   : 0.007683  
##        V4          
##  Min.   :-0.10373  
##  1st Qu.:-0.04856  
##  Median : 0.69314  
##  Mean   : 0.43318  
##  3rd Qu.: 0.78780  
##  Max.   : 0.90086
```

```r
#
# 1/100 AEP
print(original_var_list$duration_mixture_fit_par)
```

```
## [1]  0.6938297 35.9454043 48.0295307 -0.1633729
```

```r
perturbed_summary('duration_aep100_quantiles')
```

```
##       2.5%            50%            97.5%      
##  Min.   :151.2   Min.   :176.1   Min.   :246.5  
##  1st Qu.:151.4   1st Qu.:177.2   1st Qu.:247.8  
##  Median :151.7   Median :178.4   Median :249.6  
##  Mean   :151.8   Mean   :178.7   Mean   :249.7  
##  3rd Qu.:152.0   3rd Qu.:179.7   3rd Qu.:251.0  
##  Max.   :153.4   Max.   :185.5   Max.   :257.5
```

```r
relative_error_summary('duration_aep100_quantiles')
```

```
##       2.5%                50%                97.5%          
##  Min.   :0.0006586   Min.   :0.0008929   Min.   :-0.010343  
##  1st Qu.:0.0025281   1st Qu.:0.0069984   1st Qu.:-0.005222  
##  Median :0.0042821   Median :0.0139828   Median : 0.001922  
##  Mean   :0.0047055   Mean   :0.0155280   Mean   : 0.002542  
##  3rd Qu.:0.0063036   3rd Qu.:0.0210464   3rd Qu.: 0.007850  
##  Max.   :0.0155917   Max.   :0.0541647   Max.   : 0.033860
```

```r
# Errors in optimal season phi [units of years]
duration_season_phi_err = unlist(sapply(store_var_list, 
    f<-function(x) x$duration_season_phi%%1 - original_var_list$duration_season_phi%%1))
print(summary(duration_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0072590 -0.0007415  0.0001093  0.0010370  0.0021640  0.0283200
```

```r
# Look at the automatically chosen copula family. Possibly a range of families are selected?
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
##    Frank Gaussian 
##       28       72
```

**Tidal residual model checks**



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
# Check how the tideResid_mixture_fit parameters vary due to jittering
#
print(original_var_list$tideResid_mixture_fit_par)
```

```
## [1]  0.1152465  0.1146369  0.1975000 -0.1120798
```

```r
perturbed_summary('tideResid_mixture_fit_par')
```

```
##        V1               V2               V3               V4          
##  Min.   :0.1146   Min.   :0.1138   Min.   :0.1856   Min.   :-0.12766  
##  1st Qu.:0.1151   1st Qu.:0.1143   1st Qu.:0.1918   1st Qu.:-0.11976  
##  Median :0.1154   Median :0.1146   Median :0.1961   Median :-0.11271  
##  Mean   :0.1156   Mean   :0.1149   Mean   :0.2050   Mean   :-0.09792  
##  3rd Qu.:0.1161   3rd Qu.:0.1156   3rd Qu.:0.2158   3rd Qu.:-0.09009  
##  Max.   :0.1180   Max.   :0.1183   Max.   :0.3660   Max.   : 0.35644
```

```r
relative_error_summary('tideResid_mixture_fit_par')
```

```
##        V1                   V2                   V3           
##  Min.   :-0.0058208   Min.   :-0.0072074   Min.   :-0.060453  
##  1st Qu.:-0.0009931   1st Qu.:-0.0027166   1st Qu.:-0.028872  
##  Median : 0.0013453   Median :-0.0004212   Median :-0.007269  
##  Mean   : 0.0029443   Mean   : 0.0023400   Mean   : 0.038025  
##  3rd Qu.: 0.0073228   3rd Qu.: 0.0082975   3rd Qu.: 0.092482  
##  Max.   : 0.0235728   Max.   : 0.0321720   Max.   : 0.853162  
##        V4           
##  Min.   :-0.138969  
##  1st Qu.:-0.068505  
##  Median :-0.005654  
##  Mean   : 0.126327  
##  3rd Qu.: 0.196197  
##  Max.   : 4.180276
```

```r
#
# 1/100 AEP 
print(original_var_list$tideResid_aep100_quantiles)
```

```
##      2.5%       50%     97.5% 
## 0.5365498 0.6166083 0.8987738
```

```r
perturbed_summary('tideResid_aep100_quantiles')
```

```
##       2.5%             50%             97.5%       
##  Min.   :0.5343   Min.   :0.6114   Min.   :0.8731  
##  1st Qu.:0.5362   1st Qu.:0.6148   1st Qu.:0.8857  
##  Median :0.5369   Median :0.6163   Median :0.8920  
##  Mean   :0.5369   Mean   :0.6162   Mean   :0.8929  
##  3rd Qu.:0.5378   3rd Qu.:0.6178   3rd Qu.:0.8996  
##  Max.   :0.5392   Max.   :0.6206   Max.   :0.9182
```

```r
relative_error_summary('tideResid_aep100_quantiles')
```

```
##       2.5%                 50%                 97.5%          
##  Min.   :-0.0041817   Min.   :-0.0084982   Min.   :-0.028559  
##  1st Qu.:-0.0006293   1st Qu.:-0.0029380   1st Qu.:-0.014563  
##  Median : 0.0007440   Median :-0.0005009   Median :-0.007573  
##  Mean   : 0.0006931   Mean   :-0.0005972   Mean   :-0.006590  
##  3rd Qu.: 0.0023959   3rd Qu.: 0.0018980   3rd Qu.: 0.000910  
##  Max.   : 0.0049482   Max.   : 0.0064121   Max.   : 0.021636
```

```r
# Changes in seasonal phase
tideResid_season_phi_err = sapply(store_var_list, 
    f<-function(x) x$tideResid_season_phi%%1 - original_var_list$tideResid_season_phi%%1)
print(summary(tideResid_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -5.625e-06  7.181e-03  1.182e-02  9.808e-03  1.216e-02  1.640e-02
```

```r
# Look at the automatically chosen copula family. 
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

Steepness is modelled non-parametrically.  

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
# Variability in seasonal phase
steepness_season_phi_err = sapply(store_var_list, 
    f<-function(x) (x$steepness_season_phi%%1 - original_var_list$steepness_season_phi%%1))
print(summary(steepness_season_phi_err))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## -0.007146 -0.002175 -0.002152 -0.001973 -0.001690  0.002909
```

```r
# Variability in copula type
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
## Rotated Clayton 270 degrees 
##                         100
```

**Direction model checks**

Direction is modelled non-parametrically.  

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
# Variability in seasonal phase
dir_season_phi_err = sapply(store_var_list, 
    f<-function(x) (x$dir_season_phi%%1 - original_var_list$dir_season_phi%%1))
print(summary(dir_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0039490  0.0001312  0.0105900  0.0087030  0.0120300  0.0232900
```

```r
# Variability in direction/season copula family
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
##                         90                          9 
## Rotated Gumbel 270 degrees 
##                          1
```

```r
# Variability in direction_soiA copula family
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
