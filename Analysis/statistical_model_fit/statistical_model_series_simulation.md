
# **Statistical modelling of multivariate dependence in the storm event dataset**
------------------------------------------------------------------------------------

*Gareth Davies, Geoscience Australia 2017*

# Introduction
------------------

This document follows on from
[statistical_model_vine_copula.md](statistical_model_vine_copula.md)
in describing our statistical analysis of storm waves at Old Bar. 

It illustrates the process of simulating synthetic series from the fitted model. 

It is essential that the code in
[statistical_model_vine_copula.md](statistical_model_vine_copula.md)
has alread been run, and produced an Rdata file
*'Rimages/session_univariate_distributions_XXXX.Rdata'*. Here XXXX is a flag related to
whether or not perturbations were applied to the data.


Supposing the above did not generate any errors, and you have R installed,
along with all the packages required to run this code, and a copy of the
*stormwavecluster* git repository, then you should be able to re-run the
analysis here by simply copy-pasting the code. Alternatively, it can be run
with the `knit` command in the *knitr* package: 

```r
library(knitr)
knit('statistical_model_series_simulation.Rmd')
```
The above command produces a .md file with the same name for viewing, which includes
updated figures and print-outs.

To run the code in tie-breaking mode, be sure to pass the a command-line
argument matching `break_ties` to R when starting, followed by an integer ID > 0,
e.g.

    R --args --break_ties 1234

or

    Rscript script_name_here.R --break_ties 1234

Running the above commands many times is facilitated by scripts in
[../statistical_model_fit_perturbed_data/README.md](../statistical_model_fit_perturbed_data/README.md)

The basic approach followed here is to:
* **Step 1: Load the previous session**
* **Step 2: Simulate random storm properties with all dependencies**
* **Step 3: Appending `msl` and `tp1` to the synthetic series**
* **Step 4: Graphical checks of the simulated results**
* **Step 5: Generate synthetic series for bootstrapping**

# **Step 1: Load the previous session**
---------------------------------------

Here we load the session previously derived by 
[statistical_model_univariate_distributions.Rmd](statistical_model_univariate_distributions.Rmd).
As before, the code is adjusted to optionally use an analysis with random
perturbations to the data, to check the impact of ties and data discretization.


```r
# Need to re-load packages, as R does not automatically do this when re-loading
# a session
library(evmix)
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
## Loading required package: methods
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
library(logspline)
library(CDVine) # Used to set structure of C-Vine copula
```

```
## The CDVine package is no longer developed actively.
## Please consider using the more general VineCopula package
## (see https://CRAN.R-project.org/package=VineCopula),
## which extends and improves the functionality of CDVine.
```

```r
library(VineCopula) # Main copula fitting routine. 
```

```
## 
## Attaching package: 'VineCopula'
```

```
## The following objects are masked from 'package:CDVine':
## 
##     BiCopCDF, BiCopChiPlot, BiCopEst, BiCopHfunc, BiCopIndTest,
##     BiCopKPlot, BiCopLambda, BiCopMetaContour, BiCopName,
##     BiCopPar2TailDep, BiCopPar2Tau, BiCopPDF, BiCopSelect,
##     BiCopSim, BiCopTau2Par, BiCopVuongClarke
```

```r
# Here we support multiple runs with random tie-breaking of the data
# If R was passed a commandline argument 'break_ties n' on startup (with n = integer),
# then read the n'th R session matching 'Rimages/session_storm_timings_TRUE_*.Rdata'.
# That session will correspond to one of the tie-breaking sessions
if( length(grep('break_ties', commandArgs(trailingOnly=TRUE))) > 0 ){

    break_ties_with_jitter = TRUE

    # Read one of the sessions with tie-breaking
    session_n = as.numeric(commandArgs(trailingOnly=TRUE)[2])
    if(session_n < 1) stop('Invalid input ID')

    # Definitions controlling the number of years in the synthetic series
    nyears_synthetic_series = 1e+03 

}else{

    break_ties_with_jitter = FALSE
    session_n = 0

    # Definitions controlling the number of years in the synthetic series
    nyears_synthetic_series = 1e+05

}

# Make a 'title' which can appear in filenames to identify this run
run_title_id = paste0(break_ties_with_jitter, '_', session_n)

previous_R_session_file = paste0('Rimages/session_vine_copula_', run_title_id, '.Rdata')

load(previous_R_session_file)

print(previous_R_session_file)
```

```
## [1] "Rimages/session_storm_timings_FALSE_0.Rdata"
```


# **Step 2: Simulate random storm properties with all dependencies**
--------------------------------------------------------------------

Here we simulate a synthetic timeseries of events, using the event magnitude
and timings suggested before. To do this we need a function to generate the
event properties, given the time `t`. Because storm wave events cannot overlap
in time, the randomly generated event `duration` at time `t` will affect the
subsequent event timings (and therefore cannot be decoupled from the event
temporal modelling). As the event properties are affected by `soiA` we include
that in the model. `soiA` has a weak autocorrelation structure, which is treated
below.

**Below we make a function which creates an environment containing data needed
to simulate the storm properties, as well as a function which can do the
simulation.** Later will pass the latter function to the non-homogeneous
Poisson process simulation code.

```r
#' Function to initialise the 'event_creator' environment.
#'
#' This contains a function that can simulate the event properties
#' including soiA
#'
#' @param random_copula_samples function to generate random samples (in [0,1]) from the multivariate copula
#' @param stormVarFun function to compute marginals given random copula samples and conditional variables
#' @param observation_start_time time that the simulation can begin (in years)
#' @param lambda_rate_equation equation for rate model lambda (in text)
#' @param lamda_theta_par value of vector theta referred to in lambda_rate_equation
#' @param plot_soiA logical -- make plots of soiA acf/pacf?
#' @param nyears_synthetic_series -- number of years in the series we will simulate. 
#' @return an environment, which includes a function to simulate event properties
#'
build_event_creator<-function(random_copula_samples, stormVarFun, 
    observation_start_time, lambda_rate_equation, lambda_theta_par, 
    plot_soiA=FALSE, nyears_synthetic_series = nyears_synthetic_series){

    event_creator = new.env()

    event_creator$random_copula_samples = random_copula_samples
    event_creator$stormVarFun = stormVarFun
    event_creator$observation_start_time = observation_start_time
    event_creator$lambda_rate_equation = lambda_rate_equation
    event_creator$lambda_theta_par = lambda_theta_par
    event_creator$plot_soiA = plot_soiA
    event_creator$nyears_synthetic_series = nyears_synthetic_series

    event_creator$year2hours = 365.25 * 24

    with(event_creator, 
        {

            nhp = new.env()
            source('../../R/nhpoisp/nhpoisp.R', local=nhp)

            DU = new.env()
            source('../preprocessing/data_utilities.R', local=DU)

            # Include autoregressive model of soiA here
            CI = DU$read_climate_indices(
                '../../Data/Climate_Indices/ENSO/SOI_BOM.txt',
                '../../Data/Climate_Indices/AAO/AAO.txt')

            yearly_soiA = aggregate(CI$soi$index, 
                list(year=as.numeric(format(CI$soi$time, '%Y'))), 
                mean)

            # (partial) Autocorrelation plots for soi
            # Based on this we fit an auto-regressive model of order 2
            if(plot_soiA){
                par(mfrow=c(2,1))
                acf(yearly_soiA[,2], main='Autocorrelation of soiA')
                pacf(yearly_soiA[,2], main='Partial Autocorrelation of soiA')
            }

            # Fit the model and simulate a series longer than what we need
            suppressPackageStartupMessages(library(forecast))
            yearly_soiA_model = forecast::Arima(
                yearly_soiA[,2], 
                order=c(2,0,0),
                include.mean=FALSE)
            yearly_soiA_sim = forecast::simulate.Arima(
                yearly_soiA_model, 
                n=nyears_synthetic_series+10)

            # Trick to index into yearly_soiA_sim
            start_year = floor(observation_start_time)

            # lambda function
            # This will use the yearly_soiA_sim defined above
            lambda = nhp$get_lambda_function(
                lambda_theta_par,
                rate_equation = lambda_rate_equation,
                minimum_rate = 0)

            #' Pre-compute random samples from copula (to speed up code)
            #'
            #' @param N size of the random copula table
            #' @return A function which generates a single random sample
            #' from the copula table. It updates the table as required
            precompute_random_copula_samples<-function(N){

                copula_counter = 0

                random_copula_samples_precomputed = random_copula_samples(N)

                # Return a function which gives a single sample
                # and updates the random sample table as required
                get_random_copula_sample<-function(){

                    copula_counter <<- copula_counter+1

                    # Check if we need to update the table
                    if(copula_counter > dim(random_copula_samples_precomputed)[1]){
                        random_copula_samples_precomputed <<- random_copula_samples(N)
                        copula_counter <<- 1
                    }

                    return(random_copula_samples_precomputed[copula_counter,])
                }

                return(get_random_copula_sample)
            }

            # Function to get a copula sample
            get_random_copula_sample = precompute_random_copula_samples(1e+03)
        
            
            #' This function is passed to rnhpoisp
            #'
            #' @param t Event start time 
            #' @return Event properties as a named list (single row of data.frame)
            #'
            event_properties_function<-function(t){

                soiA = yearly_soiA_sim[floor(t) - start_year+1]

                # Use copula
                rand5 = get_random_copula_sample() 
            
                output = stormVarFun(
                    duration=rand5['duration'], hsig=rand5['hsig'], 
                    dir=rand5['dir'], steepness=rand5['steepness'], 
                    tideResid=rand5['tideResid'], 
                    conditional_variables=list(startyear=t, soiA=soiA)
                    )

                #
                # NOTE!
                # We must return duration in units 'years' for the timings to work!
                #
                output$duration = output$duration/(year2hours)
                output$soiA = soiA

                return(output)
            }

        }
    ) # end event_creator

    return(event_creator)

} # end function
```

**Below we use the function above to make the `event_creator`**, which will be
passed to the random non-homogeneous Poisson process simulation function. When
doing this, we have to be careful to edit the `best_nhp_model` rate equation so
that it uses the simulated `soiA` values, instead of the data values.

```r
# Make our event creator
#
# Recall that when we fit the lambda model initially, we made the soi equation
# loop repeatedly over the historical measurements of soiA.
#
# In contrast, we want the lambda model used for simulation to be like the best
# fit model.  BUT if the soi annual component is being used, we should replace
# it with:
# "theta[1] + theta[2] * yearly_soiA_sim[floor(t) - start_year + 1]"
# since this is what the event_creator environment uses to get soi as a
# function of time [see the code above where the variables yearly_soiA_sim and
# start_year are defined]
#
simulation_rate_equation = best_nhp_model$rate_equation
simulation_rate_equation = gsub(
    annual_rate_equations$soi$eqn, 
    'theta[1] + theta[2]*yearly_soiA_sim[floor(t) - start_year + 1]', 
    simulation_rate_equation,
    fixed=TRUE)
print(simulation_rate_equation)
```

```
## [1] "theta[1] + theta[2]*yearly_soiA_sim[floor(t) - start_year + 1]+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+0"
```

```r
event_creator = build_event_creator(
    random_copula_samples, 
    stormVarFun, 
    observation_start_time=obs_start_time, 
    lambda_rate_equation = simulation_rate_equation, 
    lambda_theta_par = best_nhp_model$par, 
    plot_soiA=TRUE,
    nyears_synthetic_series=nyears_synthetic_series)
```

![plot of chunk makeEventCreator](figure/makeEventCreator-1.png)

```r
# Simulate the series
# Note we add an extra gap of 'duration_gap_hours/(24*365.24)' between events
# -- since with the current definition of events, no other event can occur
# within this length of time. 
synthetic_series = nhp$rnhpoisp(
    duration = nyears_synthetic_series, 
    lambda = event_creator$lambda, 
    event_properties_function = event_creator$event_properties_function,
    observation_start_time=obs_start_time,
    extra_duration_gap=duration_gap_hours/(year2hours),
    print_progress = nyears_synthetic_series)
```

```
## [1] 1e+05
## [1] "2017-04-10 12:06:05 AEST"
## [1] 2e+05
## [1] "2017-04-10 12:12:48 AEST"
## [1] 3e+05
## [1] "2017-04-10 12:19:27 AEST"
## [1] 4e+05
## [1] "2017-04-10 12:26:06 AEST"
## [1] 5e+05
## [1] "2017-04-10 12:32:50 AEST"
## [1] 6e+05
## [1] "2017-04-10 12:39:30 AEST"
## [1] 7e+05
## [1] "2017-04-10 12:46:11 AEST"
## [1] 8e+05
## [1] "2017-04-10 12:52:48 AEST"
## [1] 9e+05
## [1] "2017-04-10 12:59:24 AEST"
## [1] 1e+06
## [1] "2017-04-10 13:06:01 AEST"
## [1] 1100000
## [1] "2017-04-10 13:12:45 AEST"
## [1] 1200000
## [1] "2017-04-10 13:19:38 AEST"
## [1] 1300000
## [1] "2017-04-10 13:26:17 AEST"
## [1] 1400000
## [1] "2017-04-10 13:32:41 AEST"
## [1] 1500000
## [1] "2017-04-10 13:39:09 AEST"
## [1] 1600000
## [1] "2017-04-10 13:45:41 AEST"
## [1] 1700000
## [1] "2017-04-10 13:52:10 AEST"
## [1] 1800000
## [1] "2017-04-10 13:58:36 AEST"
## [1] 1900000
## [1] "2017-04-10 14:05:07 AEST"
## [1] 2e+06
## [1] "2017-04-10 14:11:39 AEST"
## [1] 2100000
## [1] "2017-04-10 14:18:15 AEST"
## [1] 2200000
## [1] "2017-04-10 14:24:44 AEST"
```

```r
# Extract the information in a more convenient format
synthetic_attr = as.data.frame(attr(synthetic_series, 'event_properties'))
synthetic_attr$startyear = as.numeric(synthetic_series)

# Print the top few rows (note duration is still in years)
head(synthetic_attr)
```

```
##       duration     hsig   tideResid      dir  steepness      soiA
## 1 0.0023715065 3.579068  0.07683052 162.9246 0.02220782 0.3987766
## 2 0.0029527107 4.192566  0.16773847 165.5280 0.02445249 0.3987766
## 3 0.0044079952 3.841130  0.23527678 183.2093 0.03304199 0.3987766
## 4 0.0005577913 3.054618  0.15904881 177.5396 0.01539019 0.3987766
## 5 0.0002172042 3.042827  0.18534659 185.6018 0.03015004 0.3987766
## 6 0.0027642636 3.677255 -0.01905451 140.8115 0.01708789 3.2826634
##   startyear
## 1  1985.781
## 2  1985.839
## 3  1985.856
## 4  1985.886
## 5  1985.939
## 6  1986.050
```

**Here we graphically compare a few years of the model with a few years of data**

```r
# Plot a few years
plot_ylim = c(0, 9)
par(mfrow=c(2,1))
par(mar=c(2,4,2,1))
plot(synthetic_attr$startyear, synthetic_attr$hsig, t='h', 
    xlim=c(obs_start_time, obs_start_time+20),
    main='20 years of the synthetic series', 
    xlab='Year', ylab=bquote(H[sig]~'(m)'),
    ylim = plot_ylim)
grid()
plot(event_statistics$startyear, event_statistics$hsig, t='h',
    xlim=c(obs_start_time, obs_start_time+20),
    main='20 years of data',
    xlab='Year', ylab=bquote(H[sig]~'(m)'),
    ylim = plot_ylim)
grid()
```

![plot of chunk quickplot1](figure/quickplot1-1.png)

# **Step 3: Appending `msl` and `tp1` to the synthetic series**
---------------------------------------------------------------

We need to add in a smoothly varying MSL value to the series, as earlier we
established this was related to `soiA`. This is done using a simple approach. We
simulate an MSL value for each year based on the models fitted earlier
(`soi_SL_lm`), and then interpolate linearly between these values (assumed to
occur mid-year) to get the slowly-varying component of MSL for every storm. We
also add the monthly sea level variation to this, based on the previously
computed `smooth_tideResid_fun_stl_monthly`.


```r
#' Compute MSL given time/soi of synthetic storms
#'
#' For each year, we randomly generate an annual MSL value from the MSL vs soiA
#' regression. This is assumed to represent MSL halfway through the year.
#' Then we linearly interpolate from that series to get annual MSL for all events,
#' and add to this a seasonal component derived from the STL decomposition earlier.
#' 
#' @param output_times storm times (year)
#' @param output_soiA annual soiA during the storm
#' @param include_MSL_rise. If FALSE, only use soiA when predicting MSL. In that
#' case the year is set to 1999 (which gives MSL~0 for current data Old bar). Otherwise use the
#' output_time (rounded mid-year) to get the year
#' @return A numeric vector giving MSL during every storm
#'
compute_soi_MSL_perturbation<-function(output_times, output_soiA, 
    include_MSL_rise=FALSE){

    # Round all times to mid-year
    rounded_times = floor(output_times) + 0.5

    # Keep a single soiA, time for each year
    lrt = length(rounded_times)
    keepers = which(rounded_times[1:(lrt-1)] != rounded_times[2:lrt]) + 1
    if(rounded_times[1] != rounded_times[keepers[1]]) keepers = c(1, keepers)

    yearly_soiA_series = cbind(rounded_times[keepers], output_soiA[keepers])

    # Predict a single MSL for each year. This includes a random
    # component and an soiA related deterministic component
    #
    # To do this we use predict.lm to generative predictive intervals
    # for every point in yearly_soiA_series, **with a uniform random confidence level**.
    #
    # The underlying regression model was made in the `preprocessing` stage of the analysis
    #
    # We then randomly sample the upper or lower limit of this interval to
    # get our final value.
    lk = length(yearly_soiA_series[,1])

    if(include_MSL_rise){
        regression_year = yearly_soiA_series[,1]
    }else{
        default_year = 1999
        regression_year = rep(default_year, length(yearly_soiA_series[,1]))
    }

    random_predictions = predict(
        soi_SL_lm, 
        newdata=data.frame(soiA = yearly_soiA_series[,2], year=regression_year), 
        interval='prediction', 
        level = runif(lk))

    # Randomly sample upper or lower limit
    random_index = sample(c(0,1), size=lk, replace=TRUE) 
    random_SL = random_index * random_predictions[,'lwr'] + 
        (1-random_index) * random_predictions[,'upr']

    # Compute values of SL for each synthetic storm
    output_sl = approx(yearly_soiA_series[,1], random_SL, 
        xout=output_times, rule=2)$y

    # Add on monthly component
    # Since smooth_tideResid_fun_stl_monthly is defined from 1985-2014 and is
    # seasonally periodic, we just convert all output times to [0-1] and say
    # they are in a single year (year 2000, although any choice of year inside
    # our data coverage would be fine)
    output_sl = output_sl + smooth_tideResid_fun_stl_monthly(output_times%%1 + 2000)

    # Basic logical check
    stopifnot(length(output_sl) == lrt)

    return(output_sl)
}

# Compute an MSL series
output_sl = compute_soi_MSL_perturbation(
        output_times = synthetic_attr$startyear, 
        output_soiA = synthetic_attr$soiA,
        include_MSL_rise=FALSE)

# Append 
synthetic_attr = cbind(synthetic_attr, data.frame(msl=output_sl))

# Print the top few rows
head(synthetic_attr)
```

```
##       duration     hsig   tideResid      dir  steepness      soiA
## 1 0.0023715065 3.579068  0.07683052 162.9246 0.02220782 0.3987766
## 2 0.0029527107 4.192566  0.16773847 165.5280 0.02445249 0.3987766
## 3 0.0044079952 3.841130  0.23527678 183.2093 0.03304199 0.3987766
## 4 0.0005577913 3.054618  0.15904881 177.5396 0.01539019 0.3987766
## 5 0.0002172042 3.042827  0.18534659 185.6018 0.03015004 0.3987766
## 6 0.0027642636 3.677255 -0.01905451 140.8115 0.01708789 3.2826634
##   startyear         msl
## 1  1985.781 -0.03221010
## 2  1985.839 -0.03802943
## 3  1985.856 -0.04014978
## 4  1985.886 -0.03948911
## 5  1985.939 -0.02992279
## 6  1986.050 -0.04028984
```

**Here we back-calculate the wave period from the Airy wave dispersion relation**

```r
# Get the code 
wavedisp = new.env()
source('../../R/wave_dispersion/wave_dispersion_relation.R', local=wavedisp)

# Wavelength = hsig / (steepness) [since the latter = (hsig/(hsig/wavelength)) ]
synthetic_attr$tp1 = wavedisp$airy_period(
    lambda=synthetic_attr$hsig/synthetic_attr$steepness, 
    h=buoy_depth)
```


# **Step 4: Graphical checks of the simulated results**
-------------------------------------------------------


**Here we show a pairwise scatterplot comparisons between the model and data.**
Note that the synthetic `msl` includes seasonal variation, and so is weakly
correlated with most other storm variables, which also show seasonal variation.

```r
# Make pairwise scatterplots (ignoring time variables)
# Don't use all points (takes too long to plot!)
plot_inds = 1:5000

# Get observations with duration in years for plot
tmp_data = event_statistics[names(synthetic_attr)]
tmp_data$duration = tmp_data$duration/year2hours

DU$nice_pairs(synthetic_attr[plot_inds,], extra_data=tmp_data)
```

![plot of chunk plot1](figure/plot1-1.png)

```r
rm(tmp_data)
```

**Check whether the boundary in the relation between `hsig` and `tp1` is represented**

```r
plot_inds = 1:5000
plot(event_statistics$hsig, event_statistics$tp1, 
    xlim=c(2.5, 10), ylim=c(5, 20), 
    pch = 19, cex=0.5,
    xlab='Hsig', ylab='TP1', 
    main='Hsig vs TP1 in data (black) and model (red)')

points(synthetic_attr$hsig[plot_inds], synthetic_attr$tp1[plot_inds], 
    col=rgb(1,0,0,alpha=0.5), pch=19, cex=0.1)
```

![plot of chunk hsigtp1](figure/hsigtp1-1.png)


**Quantile-Quantile plots of data and model**

Here we run QQ-plots of the key data and model variables. If the fit is good,
the model and data should be close to the 1:1 line. We also include commented
out code to do KS tests comparing the model and data, although for efficiency
reasons that is not run here. Note that we provide our own QQ-plot, since the
default one in R always matches the extreme quantiles of both samples -- which
introduces strong bias if the sample sizes are very unequal, which is the case
here.

```r
par(mfrow=c(3,2))

# Hsig
DU$qqplot3(synthetic_attr$hsig, event_statistics$hsig, main='Hsig (m)', 
    xlab='Model', ylab='Data')
abline(0,1, col='red')

#suppressPackageStartupMessages(library(Matching))
#ks.boot(synthetic_attr$hsig[plot_inds], event_statistics$hsig)

# Duration in years
DU$qqplot3(synthetic_attr$duration*(year2hours), event_statistics$duration, 
    main='Duration (hour)', xlab='Model', ylab='Data')
abline(0,1, col='red')

#ks.boot(synthetic_attr$duration[plot_inds]*(year2hours), event_statistics$duration)

# TP1
DU$qqplot3(synthetic_attr$tp1, event_statistics$tp1, main='Tp1 (s)', 
    xlab='Model', ylab='Data')
abline(0,1, col='red')

#ks.boot(synthetic_attr$tp1[plot_inds], event_statistics$tp1)

# DIRECTION
DU$qqplot3(synthetic_attr$dir, event_statistics$dir, main='Direction (deg)', 
    xlab='Model', ylab='Data')
abline(0,1, col='red')

#ks.boot(synthetic_attr$dir[plot_inds], event_statistics$dir)

# Tidal residual (adjusted for annual MSL changes)
DU$qqplot3(synthetic_attr$tideResid, event_statistics$tideResid, 
    main='Tidal Residual (m)', xlab='Model', ylab='Data')
abline(0,1, col='red')

#ks.boot(synthetic_attr$tideResid[plot_inds], event_statistics$tideResid)
```

![plot of chunk check22](figure/check22-1.png)

**Comparison of time-between-events in the model and data**

```r
# Histogram bin breaks in hours
binBreaks = c(0, 3, 6, 9, 12, 18, 24, 36, 48, 
              seq(3*24, 50*24, len=48), 0.5*year2hours) /(year2hours)

hist(diff(event_time), breaks=binBreaks, freq=FALSE, 
    main='Time between events: Histogram (data) \n and smoothed densities (model + data)')

# Add model density
theoretical_dens = density(diff(synthetic_series), adjust=0.5, from=0)

#' R's density function produces densities which don't integrate
#' to zero when we apply the 'from' and 'to' arguments.
#' Fix that here
correct_density_integral<-function(theoretical_dens){
    density_integral = sum(theoretical_dens$y)*diff(theoretical_dens$x[1:2])

    theoretical_dens$y = theoretical_dens$y/density_integral

    return(theoretical_dens)
}

theoretical_dens = correct_density_integral(theoretical_dens)
points(theoretical_dens, t='l', col='red')

# Add data density
empirical_dens = correct_density_integral(density(diff(event_time), from=0))
points(empirical_dens, t='l', col='blue')

legend('topright', c('Model', 'Data'), lty=c(1,1), col=c('red', 'blue'), 
    title='Smoothed densities')
```

![plot of chunk time_between_events](figure/time_between_events-1.png)

**Comparison of the event time-of-year distribution in the model and data**

```r
data_tme = event_time - floor(event_time)
model_tme = synthetic_series - floor(synthetic_series)

hist(data_tme, freq=FALSE, density=10, 
    main='Event time of year in data (black) and model (red)', 
    xlab='Time of year ([0-1])', ylim=c(0, 1.4))

hist(model_tme, add=T, freq=FALSE, col='red', density=0)
```

![plot of chunk eventsTimeOfYear](figure/eventsTimeOfYear-1.png)

**Comparison of the distributions of number-of-events-each-year**

```r
hist(events_per_year_truncated, freq=FALSE, density=10, n=10,
    main='Number of events each year in data (black) and model (red)', 
    xlab='Number of events each year',
    xlim=c(0, 50))

model_counts = aggregate(synthetic_series, list(floor(synthetic_series)), 
    length)

hist(model_counts[,2], add=TRUE, freq=FALSE, n=20, col='red', density=0)
```

![plot of chunk numEventsPerYear](figure/numEventsPerYear-1.png)

# **Step 5: Simulate a number of synthetic series with the same duration as the data**

Here we simulate 10 synthetic series with the same duration as the data. Each of them
can be used to produce a single bootstrapped parameter estimate. With about 1000 such
series reasonable bootstrap type confidence intervals may be produced. These tend to work
quite well for non-extreme quantities, but may not perform as well for extremes.

**A first step is to make a function which can simulate a series.** In doing
this, it is important to update the random `soiA` and `msl` perturbations. This
is implemented below.

```r
#
#' Function to simulate data from the fitted model, with new random soiA and MSL values
#
simulate_data_from_fitted_model<-function(synthetic_series_duration){

    #
    # Re-set event creator (with a new random soiA series)
    # Add a few extra years so we have enough random soiA values
    #
    event_creator = build_event_creator(
        random_copula_samples, 
        stormVarFun, 
        observation_start_time=obs_start_time, 
        lambda_rate_equation = simulation_rate_equation, 
        lambda_theta_par = best_nhp_model$par, 
        plot_soiA=FALSE,
        nyears_synthetic_series=synthetic_series_duration + 10)

    #
    # Generate series
    #
    synthetic_series = nhp$rnhpoisp(
        duration = synthetic_series_duration,
        lambda = event_creator$lambda,
        event_properties_function = event_creator$event_properties_function,
        observation_start_time=obs_start_time,
        extra_duration_gap=duration_gap_hours/(year2hours))

    # Extract the information in a more convenient format
    synthetic_attr = as.data.frame(attr(synthetic_series, 'event_properties'))
    synthetic_attr$startyear = as.numeric(synthetic_series)

    # Convert duration to hours (not years) so it is consistent with 
    # event_statistics based on the data
    synthetic_attr$duration = synthetic_attr$duration*year2hours

    # Get wave period
    synthetic_attr$tp1 = wavedisp$airy_period(
        lambda=synthetic_attr$hsig/synthetic_attr$steepness, 
        h=buoy_depth)

    # Make synthetic MSL
    output_sl = compute_soi_MSL_perturbation(
            output_times = synthetic_attr$startyear, 
            output_soiA = synthetic_attr$soiA,
            include_MSL_rise = FALSE)
    synthetic_attr = cbind(synthetic_attr, data.frame(msl=output_sl))

    # Write out
    #output_file = paste0(synthetic_series_fitted_model_dir, '/', 'series_', 
    #    nseries*10 + i, '.csv')
    #write.table(synthetic_attr, file = output_file, sep=",", row.names=FALSE)
    return(synthetic_attr)
}
```

**Now we compute the synthetic series.** Only ten series are generated here,
all having the same length as the data, although for proper uncertainty
quantification many more series would be needed (e.g. 1000).

```r
n_series = 10
simulated_series_list = vector(mode='list', length=n_series) 
for(i in 1:n_series){
    simulated_series_list[[i]] = simulate_data_from_fitted_model(
        synthetic_series_duration=data_duration_years)
}
```


**Save outputs for later use**

```r
dir.create('Rimages', showWarnings=FALSE)
Rimage_title = paste0('Rimages/session_series_simulation_', run_title_id, '.Rdata')
save.image(Rimage_title)

# Print R version and package info
print(sessionInfo())
```

```
## R version 3.3.1 (2016-06-21)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.5 LTS
## 
## locale:
##  [1] LC_CTYPE=en_AU.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_AU.UTF-8        LC_COLLATE=en_AU.UTF-8    
##  [5] LC_MONETARY=en_AU.UTF-8    LC_MESSAGES=en_AU.UTF-8   
##  [7] LC_PAPER=en_AU.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_AU.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] methods   splines   stats     graphics  grDevices utils     datasets 
## [8] base     
## 
## other attached packages:
## [1] forecast_8.0          VineCopula_2.1.2.9000 CDVine_1.4           
## [4] logspline_2.1.9       evmix_2.6             SparseM_1.7          
## [7] gsl_1.9-10.3          MASS_7.3-45           knitr_1.15.1         
## 
## loaded via a namespace (and not attached):
##  [1] pcaPP_1.9-61      Rcpp_0.12.9       highr_0.6        
##  [4] plyr_1.8.4        iterators_1.0.8   tseries_0.10-38  
##  [7] tools_3.3.1       evaluate_0.10     gtable_0.2.0     
## [10] lattice_0.20-33   pspline_1.0-17    Matrix_1.2-6     
## [13] foreach_1.4.3     igraph_1.0.1      parallel_3.3.1   
## [16] mvtnorm_1.0-6     copula_0.999-16   stringr_1.1.0    
## [19] lmtest_0.9-35     stats4_3.3.1      grid_3.3.1       
## [22] nnet_7.3-12       ADGofTest_0.3     ggplot2_2.1.0    
## [25] magrittr_1.5      scales_0.4.0      codetools_0.2-14 
## [28] stabledist_0.7-1  timeDate_3012.100 colorspace_1.3-2 
## [31] fracdiff_1.4-2    numDeriv_2016.8-1 quadprog_1.5-5   
## [34] stringi_1.1.2     network_1.13.0    doParallel_1.0.10
## [37] munsell_0.4.3     zoo_1.7-13
```

