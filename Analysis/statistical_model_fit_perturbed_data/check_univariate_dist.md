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

# Function to interpolate rate/value curve at rate = 1/50,
# with log-transform of values
rate_1_in_50<-function(rate, value){
    exp(approx(log(rate), log(value), xout=log(1/50))$y)
}

ml_1_in_50<-function(variable_name){
    


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
##  1st Qu.:0.8415   1st Qu.:1.016   1st Qu.:1.270   1st Qu.:-0.2146  
##  Median :0.8420   Median :1.016   Median :1.270   Median :-0.2144  
##  Mean   :0.8419   Mean   :1.016   Mean   :1.270   Mean   :-0.2145  
##  3rd Qu.:0.8423   3rd Qu.:1.017   3rd Qu.:1.271   3rd Qu.:-0.2143  
##  Max.   :0.8436   Max.   :1.019   Max.   :1.271   Max.   :-0.2138
```

```r
# Relative error in parameters summary
relative_error_summary('hsig_mixture_fit_par')
```

```
##        V1                   V2                   V3            
##  Min.   :-0.0026812   Min.   :-0.0019291   Min.   :-6.728e-04  
##  1st Qu.:-0.0008067   1st Qu.:-0.0002382   1st Qu.:-1.812e-04  
##  Median :-0.0001407   Median : 0.0001628   Median :-4.343e-05  
##  Mean   :-0.0002812   Mean   : 0.0003094   Mean   :-6.306e-05  
##  3rd Qu.: 0.0002182   3rd Qu.: 0.0008843   3rd Qu.: 5.466e-05  
##  Max.   : 0.0017302   Max.   : 0.0029132   Max.   : 3.961e-04  
##        V4            
##  Min.   :-0.0038367  
##  1st Qu.:-0.0012270  
##  Median :-0.0001839  
##  Mean   :-0.0004272  
##  3rd Qu.: 0.0002547  
##  Max.   : 0.0025782
```

```r
#

# Original 1/100 quantile
print(original_var_list$hsig_aep100_quantiles)
```

```
##     2.5%      50%    97.5% 
## 7.071843 7.568026 8.974433
```

```r
# Summary of perturbed 1/100 quantile
perturbed_summary('hsig_aep100_quantiles')
```

```
##       2.5%            50%            97.5%      
##  Min.   :7.069   Min.   :7.559   Min.   :8.925  
##  1st Qu.:7.071   1st Qu.:7.561   1st Qu.:8.941  
##  Median :7.071   Median :7.562   Median :8.948  
##  Mean   :7.071   Mean   :7.562   Mean   :8.948  
##  3rd Qu.:7.072   3rd Qu.:7.563   3rd Qu.:8.954  
##  Max.   :7.074   Max.   :7.565   Max.   :8.977
```

```r
# Summary of relative error
relative_error_summary('hsig_aep100_quantiles')
```

```
##       2.5%                 50%                 97.5%           
##  Min.   :-3.949e-04   Min.   :-0.0012372   Min.   :-0.0054722  
##  1st Qu.:-1.754e-04   1st Qu.:-0.0009442   1st Qu.:-0.0037407  
##  Median :-5.987e-05   Median :-0.0008399   Median :-0.0029923  
##  Mean   :-6.447e-05   Mean   :-0.0008301   Mean   :-0.0029992  
##  3rd Qu.: 4.534e-05   3rd Qu.:-0.0007070   3rd Qu.:-0.0022433  
##  Max.   : 3.508e-04   Max.   :-0.0003546   Max.   : 0.0003249
```

```r
# ML rates of 1/50 event
orig_1_in_50 = rate_1_in_50(original_var_list$hsig_mixture_fit_bayes_rates, 
    original_var_list$hsig_mixture_fit_ml)
print(orig_1_in_50)
```

```
## [1] 7.234595
```

```r
perturb_1_in_50 = lapply(store_var_list, f<-function(x){
        rate_1_in_50(x$hsig_mixture_fit_bayes_rates, x$hsig_mixture_fit_ml)})
print(summary(unlist(perturb_1_in_50) - orig_1_in_50))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0010350 -0.0003212 -0.0001364 -0.0001339  0.0001271  0.0007237
```

```r
# Variability in optimal season phi [units of years -- typical value O(1 day)]
hsig_season_phi_err = unlist(sapply(store_var_list, 
    f<-function(x) x$hsig_season_phi%%1 - original_var_list$hsig_season_phi%%1))
print(summary(hsig_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -7.414e-04  2.307e-05  7.288e-05  2.417e-03  2.512e-03  1.360e-02
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
##                        16                        80 
## Rotated Gumbel 90 degrees 
##                         4
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
##  Min.   :0.5661   Min.   :36.63   Min.   : 7.084   Min.   :-0.17720  
##  1st Qu.:0.5979   1st Qu.:37.45   1st Qu.: 8.261   1st Qu.:-0.17248  
##  Median :0.6202   Median :46.20   Median :11.891   Median :-0.04938  
##  Mean   :0.6306   Mean   :44.96   Mean   :25.220   Mean   :-0.09381  
##  3rd Qu.:0.6682   3rd Qu.:51.93   3rd Qu.:47.408   3rd Qu.:-0.03296  
##  Max.   :0.6817   Max.   :59.00   Max.   :48.378   Max.   :-0.02275
```

```r
relative_error_summary('duration_mixture_fit_par')
```

```
##        V1                 V2                V3           
##  Min.   :-0.18403   Min.   :0.01905   Min.   :-0.852500  
##  1st Qu.:-0.13823   1st Qu.:0.04186   1st Qu.:-0.827999  
##  Median :-0.10613   Median :0.28527   Median :-0.752413  
##  Mean   :-0.09108   Mean   :0.25072   Mean   :-0.474904  
##  3rd Qu.:-0.03692   3rd Qu.:0.44462   3rd Qu.:-0.012934  
##  Max.   :-0.01742   Max.   :0.64138   Max.   : 0.007257  
##        V4          
##  Min.   :-0.08463  
##  1st Qu.:-0.05576  
##  Median : 0.69775  
##  Mean   : 0.42580  
##  3rd Qu.: 0.79827  
##  Max.   : 0.86076
```

```r
#
# 1/100 AEP
print(original_var_list$duration_aep100_quantiles)
```

```
##     2.5%      50%    97.5% 
## 151.3099 175.9107 248.3765
```

```r
perturbed_summary('duration_aep100_quantiles')
```

```
##       2.5%            50%            97.5%      
##  Min.   :151.1   Min.   :176.3   Min.   :246.3  
##  1st Qu.:151.5   1st Qu.:177.4   1st Qu.:248.4  
##  Median :151.7   Median :178.5   Median :249.5  
##  Mean   :151.8   Mean   :178.8   Mean   :249.9  
##  3rd Qu.:152.1   3rd Qu.:179.7   3rd Qu.:251.0  
##  Max.   :153.4   Max.   :185.2   Max.   :257.5
```

```r
relative_error_summary('duration_aep100_quantiles')
```

```
##       2.5%                50%               97.5%           
##  Min.   :-0.001185   Min.   :0.002184   Min.   :-0.0082461  
##  1st Qu.: 0.001036   1st Qu.:0.008309   1st Qu.:-0.0000616  
##  Median : 0.002814   Median :0.014822   Median : 0.0043501  
##  Mean   : 0.003265   Mean   :0.016616   Mean   : 0.0062148  
##  3rd Qu.: 0.004940   3rd Qu.:0.021335   3rd Qu.: 0.0106906  
##  Max.   : 0.013571   Max.   :0.052894   Max.   : 0.0367782
```

```r
# ML rates of 1/50 event
orig_1_in_50 = rate_1_in_50(original_var_list$duration_mixture_fit_bayes_rates, 
    original_var_list$duration_mixture_fit_ml)
print(orig_1_in_50)
```

```
## [1] 157.1726
```

```r
perturb_1_in_50 = lapply(store_var_list, f<-function(x){
        rate_1_in_50(x$duration_mixture_fit_bayes_rates, x$duration_mixture_fit_ml)})
print(summary(unlist(perturb_1_in_50) - orig_1_in_50))
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -0.93970 -0.09341 14.68000 10.09000 18.32000 22.30000
```

```r
# Bayesian rates, lower
orig_1_in_50 = rate_1_in_50(original_var_list$duration_mixture_fit_bayes_rates, 
    original_var_list$duration_mixture_fit_bayes_lower_q)
print(orig_1_in_50)
```

```
## [1] 144.4226
```

```r
perturb_1_in_50 = lapply(store_var_list, f<-function(x){
        rate_1_in_50(x$duration_mixture_fit_bayes_rates, x$duration_mixture_fit_bayes_lower_q)})
print(summary(unlist(perturb_1_in_50) - orig_1_in_50))
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -0.02123  0.25340  0.43560  0.47950  0.67210  1.51700
```

```r
# Bayesian rates, upper
orig_1_in_50 = rate_1_in_50(original_var_list$duration_mixture_fit_bayes_rates, 
    original_var_list$duration_mixture_fit_bayes_upper_q)
print(orig_1_in_50)
```

```
## [1] 221.7531
```

```r
perturb_1_in_50 = lapply(store_var_list, f<-function(x){
        rate_1_in_50(x$duration_mixture_fit_bayes_rates, x$duration_mixture_fit_bayes_upper_q)})
print(summary(unlist(perturb_1_in_50) - orig_1_in_50))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -1.3280  0.1169  1.0540  1.3520  2.2140  7.1100
```

```r
# Errors in optimal season phi [units of years]
duration_season_phi_err = unlist(sapply(store_var_list, 
    f<-function(x) x$duration_season_phi%%1 - original_var_list$duration_season_phi%%1))
print(summary(duration_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0075940 -0.0006901  0.0001174  0.0009504  0.0055980  0.0222800
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
##                        29                        70 
## Rotated Gumbel 90 degrees 
##                         1
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
## [1]  0.1154248  0.1144187  0.1973124 -0.1100628
```

```r
perturbed_summary('tideResid_mixture_fit_par')
```

```
##        V1               V2               V3               V4          
##  Min.   :0.1148   Min.   :0.1138   Min.   :0.1850   Min.   :-0.12754  
##  1st Qu.:0.1151   1st Qu.:0.1143   1st Qu.:0.1900   1st Qu.:-0.12083  
##  Median :0.1153   Median :0.1146   Median :0.1950   Median :-0.11496  
##  Mean   :0.1155   Mean   :0.1148   Mean   :0.2005   Mean   :-0.10858  
##  3rd Qu.:0.1161   3rd Qu.:0.1155   3rd Qu.:0.2150   3rd Qu.:-0.09010  
##  Max.   :0.1163   Max.   :0.1160   Max.   :0.2222   Max.   :-0.08048
```

```r
relative_error_summary('tideResid_mixture_fit_par')
```

```
##        V1                   V2                  V3          
##  Min.   :-0.0058395   Min.   :-0.005787   Min.   :-0.06235  
##  1st Qu.:-0.0027071   1st Qu.:-0.001384   1st Qu.:-0.03684  
##  Median :-0.0008415   Median : 0.001184   Median :-0.01147  
##  Mean   : 0.0006113   Mean   : 0.003133   Mean   : 0.01628  
##  3rd Qu.: 0.0055840   3rd Qu.: 0.009555   3rd Qu.: 0.08973  
##  Max.   : 0.0080013   Max.   : 0.013474   Max.   : 0.12615  
##        V4          
##  Min.   :-0.15883  
##  1st Qu.:-0.09783  
##  Median :-0.04448  
##  Mean   : 0.01343  
##  3rd Qu.: 0.18138  
##  Max.   : 0.26880
```

```r
#
# 1/100 AEP 
print(original_var_list$tideResid_aep100_quantiles)
```

```
##      2.5%       50%     97.5% 
## 0.5371929 0.6163213 0.8986268
```

```r
perturbed_summary('tideResid_aep100_quantiles')
```

```
##       2.5%             50%             97.5%       
##  Min.   :0.5347   Min.   :0.6114   Min.   :0.8679  
##  1st Qu.:0.5360   1st Qu.:0.6147   1st Qu.:0.8854  
##  Median :0.5370   Median :0.6163   Median :0.8925  
##  Mean   :0.5370   Mean   :0.6162   Mean   :0.8919  
##  3rd Qu.:0.5380   3rd Qu.:0.6181   3rd Qu.:0.8981  
##  Max.   :0.5393   Max.   :0.6203   Max.   :0.9157
```

```r
relative_error_summary('tideResid_aep100_quantiles')
```

```
##       2.5%                 50%                 97.5%           
##  Min.   :-0.0047293   Min.   :-0.0079290   Min.   :-0.0341791  
##  1st Qu.:-0.0022660   1st Qu.:-0.0026516   1st Qu.:-0.0147044  
##  Median :-0.0002779   Median :-0.0001122   Median :-0.0067825  
##  Mean   :-0.0003945   Mean   :-0.0002329   Mean   :-0.0075196  
##  3rd Qu.: 0.0014716   3rd Qu.: 0.0028649   3rd Qu.:-0.0005776  
##  Max.   : 0.0039224   Max.   : 0.0064233   Max.   : 0.0190254
```

```r
# Changes in seasonal phase
tideResid_season_phi_err = sapply(store_var_list, 
    f<-function(x) x$tideResid_season_phi%%1 - original_var_list$tideResid_season_phi%%1)
print(summary(tideResid_season_phi_err))
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -2.627e-06  7.243e-03  1.080e-02  9.624e-03  1.208e-02  1.454e-02
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
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -6.674e-03 -2.182e-03 -2.142e-03 -1.977e-03 -1.700e-03 -6.176e-05
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
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## -0.004037  0.005411  0.011640  0.010550  0.012050  0.023470
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
##                         91                          9
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
