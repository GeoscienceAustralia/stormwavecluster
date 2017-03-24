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
# Function to summarise parameters from perturbed runs
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

#
# Function to summarise "relative errors" in perturbed runs
# i.e. (perturbed - original)/abs(original)
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
##  Min.   :0.8399   Min.   :1.014   Min.   :1.270   Min.   :-0.2152  
##  1st Qu.:0.8414   1st Qu.:1.016   1st Qu.:1.270   1st Qu.:-0.2147  
##  Median :0.8419   Median :1.016   Median :1.270   Median :-0.2145  
##  Mean   :0.8418   Mean   :1.017   Mean   :1.270   Mean   :-0.2145  
##  3rd Qu.:0.8424   3rd Qu.:1.017   3rd Qu.:1.271   3rd Qu.:-0.2143  
##  Max.   :0.8437   Max.   :1.019   Max.   :1.271   Max.   :-0.2137
```

```r
# Relative error in parameters summary
relative_error_summary('hsig_mixture_fit_par')
```

```
##        V1                   V2                   V3            
##  Min.   :-0.0027110   Min.   :-0.0020800   Min.   :-6.771e-04  
##  1st Qu.:-0.0008945   1st Qu.:-0.0002465   1st Qu.:-2.054e-04  
##  Median :-0.0002992   Median : 0.0003185   Median :-6.104e-05  
##  Mean   :-0.0003835   Mean   : 0.0004308   Mean   :-8.769e-05  
##  3rd Qu.: 0.0002303   3rd Qu.: 0.0010237   3rd Qu.: 6.098e-05  
##  Max.   : 0.0018790   Max.   : 0.0029767   Max.   : 4.667e-04  
##        V4            
##  Min.   :-0.0037486  
##  1st Qu.:-0.0014571  
##  Median :-0.0005010  
##  Mean   :-0.0005783  
##  3rd Qu.: 0.0003216  
##  Max.   : 0.0029474
```

```r
#

# Original 1/100 quantile
print(original_var_list$hsig_aep100_quantiles)
```

```
##     2.5%      50%    97.5% 
## 7.073745 7.566690 8.964898
```

```r
# Summary of perturbed 1/100 quantile
perturbed_summary('hsig_aep100_quantiles')
```

```
##       2.5%            50%            97.5%      
##  Min.   :7.068   Min.   :7.558   Min.   :8.918  
##  1st Qu.:7.071   1st Qu.:7.561   1st Qu.:8.941  
##  Median :7.071   Median :7.562   Median :8.947  
##  Mean   :7.071   Mean   :7.562   Mean   :8.947  
##  3rd Qu.:7.072   3rd Qu.:7.563   3rd Qu.:8.955  
##  Max.   :7.074   Max.   :7.566   Max.   :8.964
```

```r
# Summary of relative error
relative_error_summary('hsig_aep100_quantiles')
```

```
##       2.5%                 50%                 97.5%           
##  Min.   :-8.245e-04   Min.   :-1.093e-03   Min.   :-0.0052769  
##  1st Qu.:-4.531e-04   1st Qu.:-8.093e-04   1st Qu.:-0.0026401  
##  Median :-3.403e-04   Median :-6.576e-04   Median :-0.0019418  
##  Mean   :-3.537e-04   Mean   :-6.574e-04   Mean   :-0.0019494  
##  3rd Qu.:-2.599e-04   3rd Qu.:-5.174e-04   3rd Qu.:-0.0010779  
##  Max.   :-1.265e-05   Max.   :-2.534e-05   Max.   :-0.0001129
```

```r
# Variability in optimal season phi [units of years -- typical value O(1 day)]
hsig_season_phi_err = unlist(sapply(store_var_list, 
    f<-function(x) x$hsig_season_phi%%1 - original_var_list$hsig_season_phi%%1))
print(summary(hsig_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -3.120e-04  2.618e-05  6.018e-05  2.094e-03  1.986e-03  1.352e-02
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
##                        14                        84 
## Rotated Gumbel 90 degrees 
##                         2
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
##  Min.   :0.5759   Min.   :36.00   Min.   : 7.127   Min.   :-0.18134  
##  1st Qu.:0.6128   1st Qu.:37.35   1st Qu.:11.302   1st Qu.:-0.17210  
##  Median :0.6225   Median :45.68   Median :11.998   Median :-0.05064  
##  Mean   :0.6359   Mean   :43.88   Mean   :27.062   Mean   :-0.10065  
##  3rd Qu.:0.6699   3rd Qu.:47.06   3rd Qu.:47.566   3rd Qu.:-0.04705  
##  Max.   :0.6932   Max.   :56.93   Max.   :48.539   Max.   :-0.02255
```

```r
relative_error_summary('duration_mixture_fit_par')
```

```
##        V1                   V2                 V3           
##  Min.   :-0.1700023   Min.   :0.001448   Min.   :-0.851620  
##  1st Qu.:-0.1168128   1st Qu.:0.038987   1st Qu.:-0.764679  
##  Median :-0.1027397   Median :0.270842   Median :-0.750194  
##  Mean   :-0.0834969   Mean   :0.220694   Mean   :-0.436558  
##  3rd Qu.:-0.0345601   3rd Qu.:0.309347   3rd Qu.:-0.009646  
##  Max.   :-0.0009158   Max.   :0.583670   Max.   : 0.010608  
##        V4          
##  Min.   :-0.10999  
##  1st Qu.:-0.05344  
##  Median : 0.69005  
##  Mean   : 0.38390  
##  3rd Qu.: 0.71203  
##  Max.   : 0.86195
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
##  Min.   :151.1   Min.   :175.6   Min.   :246.0  
##  1st Qu.:151.5   1st Qu.:177.3   1st Qu.:247.8  
##  Median :151.7   Median :178.3   Median :249.6  
##  Mean   :151.7   Mean   :178.6   Mean   :249.6  
##  3rd Qu.:152.0   3rd Qu.:179.4   3rd Qu.:250.8  
##  Max.   :152.9   Max.   :183.9   Max.   :255.7
```

```r
relative_error_summary('duration_aep100_quantiles')
```

```
##       2.5%                 50%                97.5%           
##  Min.   :-0.0006114   Min.   :-0.001799   Min.   :-0.0069794  
##  1st Qu.: 0.0018944   1st Qu.: 0.007463   1st Qu.: 0.0005652  
##  Median : 0.0031132   Median : 0.013374   Median : 0.0077083  
##  Mean   : 0.0035608   Mean   : 0.014723   Mean   : 0.0077081  
##  3rd Qu.: 0.0050792   3rd Qu.: 0.019730   3rd Qu.: 0.0126087  
##  Max.   : 0.0114902   Max.   : 0.045189   Max.   : 0.0321220
```

```r
# Errors in optimal season phi [units of years]
duration_season_phi_err = unlist(sapply(store_var_list, 
    f<-function(x) x$duration_season_phi%%1 - original_var_list$duration_season_phi%%1))
print(summary(duration_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0072540 -0.0006843  0.0001347  0.0015950  0.0056200  0.0283600
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
##                     Frank                  Gaussian 
##                        36                        62 
## Rotated Gumbel 90 degrees 
##                         2
```

**Tidal residual model checks**



```r
# Make a plot
plot_ylim = range(c(original_var_list$tideResid_mixture_fit_bayes_lower_q,
    original_var_list$tideResid_mixture_fit_bayes_upper_q))
```

```
## Warning in min(x, na.rm = na.rm): no non-missing arguments to min;
## returning Inf
```

```
## Warning in max(x, na.rm = na.rm): no non-missing arguments to max;
## returning -Inf
```

```r
xs = original_var_list$tideResid_mixture_fit_bayes_rates
plot(xs, original_var_list$tideResid_mixture_fit_ml, 
    t='l', ylim=plot_ylim, main='Tidal residual fit (all overplotted)', xlim=c(max(xs), min(xs)),
    log='xy', xlab='Event rate (number per year)', ylab='m')
```

```
## Error in plot.window(...): need finite 'ylim' values
```

![plot of chunk tidal_residualA](figure/tidal_residualA-1.png)

```r
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

#
# Check how the tideResid_mixture_fit parameters vary due to jittering
#
print(original_var_list$tideResid_mixture_fit_par)
```

```
## [1]  0.1150170  0.1148036  0.1974718 -0.1132414
```

```r
perturbed_summary('tideResid_mixture_fit_par')
```

```
##        V1               V2               V3               V4          
##  Min.   :0.1147   Min.   :0.1137   Min.   :0.1836   Min.   :-0.13056  
##  1st Qu.:0.1151   1st Qu.:0.1143   1st Qu.:0.1915   1st Qu.:-0.12059  
##  Median :0.1154   Median :0.1146   Median :0.1961   Median :-0.11459  
##  Mean   :0.1155   Mean   :0.1148   Mean   :0.2039   Mean   :-0.09936  
##  3rd Qu.:0.1160   3rd Qu.:0.1154   3rd Qu.:0.2138   3rd Qu.:-0.09114  
##  Max.   :0.1178   Max.   :0.1182   Max.   :0.3664   Max.   : 0.37464
```

```r
relative_error_summary('tideResid_mixture_fit_par')
```

```
##        V1                  V2                   V3           
##  Min.   :-0.002514   Min.   :-0.0092246   Min.   :-0.070425  
##  1st Qu.: 0.001036   1st Qu.:-0.0044007   1st Qu.:-0.030089  
##  Median : 0.003292   Median :-0.0021547   Median :-0.006833  
##  Mean   : 0.004573   Mean   : 0.0002522   Mean   : 0.032326  
##  3rd Qu.: 0.008956   3rd Qu.: 0.0054620   3rd Qu.: 0.082929  
##  Max.   : 0.024465   Max.   : 0.0297268   Max.   : 0.855522  
##        V4          
##  Min.   :-0.15291  
##  1st Qu.:-0.06485  
##  Median :-0.01195  
##  Mean   : 0.12261  
##  3rd Qu.: 0.19517  
##  Max.   : 4.30836
```

```r
#
# 1/100 AEP 
print(original_var_list$tideResid_aep100_quantiles)
```

```
##      2.5%       50%     97.5% 
## 0.5366994 0.6168376 0.8978648
```

```r
perturbed_summary('tideResid_aep100_quantiles')
```

```
##       2.5%             50%             97.5%       
##  Min.   :0.5349   Min.   :0.6115   Min.   :0.8700  
##  1st Qu.:0.5363   1st Qu.:0.6150   1st Qu.:0.8865  
##  Median :0.5371   Median :0.6162   Median :0.8920  
##  Mean   :0.5371   Mean   :0.6162   Mean   :0.8924  
##  3rd Qu.:0.5379   3rd Qu.:0.6176   3rd Qu.:0.8972  
##  Max.   :0.5395   Max.   :0.6217   Max.   :0.9247
```

```r
relative_error_summary('tideResid_aep100_quantiles')
```

```
##       2.5%                 50%                 97.5%           
##  Min.   :-0.0032930   Min.   :-0.0086723   Min.   :-0.0310866  
##  1st Qu.:-0.0008351   1st Qu.:-0.0029125   1st Qu.:-0.0126576  
##  Median : 0.0008260   Median :-0.0009671   Median :-0.0065477  
##  Mean   : 0.0006699   Mean   :-0.0009823   Mean   :-0.0060552  
##  3rd Qu.: 0.0021690   3rd Qu.: 0.0011894   3rd Qu.:-0.0007056  
##  Max.   : 0.0051453   Max.   : 0.0078647   Max.   : 0.0298562
```

```r
# Changes in seasonal phase
tideResid_season_phi_err = sapply(store_var_list, 
    f<-function(x) x$tideResid_season_phi%%1 - original_var_list$tideResid_season_phi%%1)
print(summary(tideResid_season_phi_err))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 1.365e-05 9.797e-03 1.196e-02 1.034e-02 1.269e-02 1.455e-02
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
## -0.006675 -0.002191 -0.002152 -0.001969 -0.001734  0.001092
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
## -0.0048750  0.0004312  0.0115000  0.0101300  0.0210700  0.0233000
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
##                         87                         12 
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
