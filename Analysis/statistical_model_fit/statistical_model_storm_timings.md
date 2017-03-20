
# **Statistical modelling of the storm event dataset: Storm event timings**
-------------------------------------------------------------

*Gareth Davies, Geoscience Australia 2017*

# Introduction
------------------

This document follows on from
[../preprocessing/extract_storm_events.md](../preprocessing/extract_storm_events.md)
in describing our statistical analysis of storm waves at Old Bar. 

It illustrates the process of fitting the storm event timing statistical model to the data.

It is essential that the scripts in [../preprocessing](../preprocessing) (or
equivalent in
[../preprocessing_perturbed_data](../preprocessing_perturbed_data) ) have
alread been run, and produced an RDS file
*'../preprocessing/Derived_data/event_statistics_XXXX.RDS'* (or equivalent)
where XXXX is a flag for recording the particulars of any data perturbation
that was applied. 


Supposing the prerequisites have been run, and you have R installed,
along with all the packages required to run this code, and a copy of the
*stormwavecluster* git repository, then you should be able to re-run the
analysis here by simply copy-pasting the code. Alternatively, it can be run
with the `knit` command in the *knitr* package: 

```r
library(knitr)
knit('statistical_model_storm_timings.Rmd')
```

To run the code in tie-breaking mode, be sure to pass the a command-line
argument matching `break_ties` to R when starting, followed by an integer ID > 0,
e.g.

    R --args --break_ties 1234

or

    Rscript script_name_here.R --break_ties 1234

Running the above commands many times is facilitated by scripts in
[../statistical_model_fit_perturbed_data/README.md](../statistical_model_fit_perturbed_data/README.md)


The basic approach followed here is to:
* **Step 1: Get the preprocessed data and optionally break ties in the data**
* **Step 2: Compute the wave steepness**
* **Step 3: Exploration of the annual number of storms**
* **Step 4: Modelling the storm event timings as a non-homogeneous Poisson process**

Later we will use the statistical model to simulate synthetic storm event
time-series.

# **Step 1: Get the preprocessed data and optionally break ties in the data**
----------------------------------------------------------------------


**Below we read the data from earlier steps of the analysis**

```r
# Get the data_utilities functions (in their own environment to keep the
# namespace clean)
DU = new.env()
source('../preprocessing/data_utilities.R', local=DU, chdir=TRUE) 

# Useful number to convert from years to hours (ignoring details of leap-years)
year2hours = 365.25*24

# Optionally remove ties in the event statistics by jittering.
# To do this, a commandline argument matching 'break_ties' must
# have been passed to R before running. Otherwise no jittering is applied. 
# We design this way to facilitate running many random jobs with scripts, using
# the same code.
if( length(grep('break_ties', commandArgs(trailingOnly=TRUE))) > 0){
    # Apply jittering
    break_ties_with_jitter = TRUE
    
    # Get the ID for this run
    session_n = as.numeric(commandArgs(trailingOnly=TRUE)[2])

    if(session_n < 1) stop('Invalid tie-breaking ID value')
    
}else{

    # No jittering -- use the raw data
    break_ties_with_jitter = FALSE

    session_n = 0

}

# Make a 'title' which can appear in filenames to identify this run
run_title_id = paste0(break_ties_with_jitter, '_', session_n)

#
# Read the last session
#
if(break_ties_with_jitter){
    # Read data produced by pre-processing scripts
    event_statistics_list = readRDS(paste0(
        '../preprocessing_perturbed_data/Derived_data/event_statistics_', 
        run_title_id, '.RDS'))
}else{
    # Read data produced by pre-processing scripts
    event_statistics_list = readRDS(paste0(
        '../preprocessing/Derived_data/event_statistics_', 
        run_title_id, '.RDS'))
}

# Extract variables that we need from previous analysis
for(varname in names(event_statistics_list)){
    assign(varname, event_statistics_list[[varname]])
}

# No need to keep event statistics_list, as we have extracted all required
# variables
rm(event_statistics_list)


# Look at the variables we have
ls()
```

```
##  [1] "break_ties_with_jitter"           "CI_annual_fun"                   
##  [3] "data_duration_years"              "DU"                              
##  [5] "duration_gap_hours"               "duration_offset_hours"           
##  [7] "duration_threshold_hours"         "event_statistics"                
##  [9] "hsig_threshold"                   "obs_start_time_strptime"         
## [11] "run_title_id"                     "session_n"                       
## [13] "smooth_tideResid_fun_stl_monthly" "soi_SL_lm"                       
## [15] "varname"                          "year2hours"
```

```r
# Record whether or not we are breaking ties
print(break_ties_with_jitter)
```

```
## [1] FALSE
```

If our event statistics are subject to rounding (introducing 'ties' or repeated
values into the data), then some statistical methods designed for continuous
data may perform badly. For instance, our storm duration data is always in
multiples of one hour (because we use hourly data), and so there are many
storms with durations of 1, 2, 3... hours. These 'ties' lead to ambiguity in
the definition of data ranks, can sometimes result in poor performance for
statistical methods which assume continuous data, because for continuous data,
ties have probability zero, so data ranks are always well defined. While this
is not always a problem, it needs to be checked.

Therefore, **below we optionally perturb the `event_statistics` to remove
ties**. Recall that earlier in the analysis, the optional perturbation was
applied to the raw `hsig` values and the raw tidal observations, so there is no
need to perturb `hsig` or `tideResid` again below.  The remaining variables
have ties due to limited measurement resolution. We thus apply a perturbation
of 1/2 hour for `duration` and `startyear`, and 1/2 degree for `dir`, as these
represent half of the bin-width of the measured data we have (1 hour / 1 degree
increments).  For `tp1` (which has the most ties, and only 40 unique values),
the bins are irregularly spaced without an obvious pattern. The median distance
between unique `tp1` values after sorting is 0.25, with a maximum of 1.06, and
a minimum of 0.01.  Therefore, a uniform perturbation of plus/minus 0.1 second
is applied to `tp1`. 

```r
#
# Jitter of variables described in the text above -- recalling that hsig and tidal measurement 
# were already jittered
#
default_jitter_vars = c( 'hsig', 'duration', 'tideResid', 'dir', 'tp1')
default_jitter_amounts = c( 0.0,        0.5,         0.0,   0.5,   0.1) 
names(default_jitter_amounts) = default_jitter_vars


#' Make a function which will return a jittered version of the original
#' event_statistics
#'
#' See comments above regarding default jitter variables and amounts
#'
#' @param event_statistics_orig data.frame with the 'raw' event statistics
#' @param default_jitter_vars character vector with names of variables to
#' jitter. 'duration' will lead to both 'duration' and 'startyear' being
#' changed, with care to ensure consistency
#' @param default_jitter_amounts numeric vector with jitter amounts for each variable.
#' Note that for 'hsig', the jitter is interpreted as a fraction of the original value,
#' and all other jitters are absolute.
#' @return function which generates jittered event statistics
#'
make_jitter_event_statistics_function<-function(
    event_statistics_orig,
    default_jitter_vars,
    default_jitter_amounts){

    # Bring arguments into this environment
    event_statistics_orig = event_statistics_orig
    default_jitter_vars = default_jitter_vars
    default_jitter_amounts = as.numeric(default_jitter_amounts) # Strip names if provided

    # Function that will jitter the original event_statistics
    jitter_event_statistics_function<-function(
        jitter_vars = default_jitter_vars,
        jitter_amounts = default_jitter_amounts
        ){

        event_statistics = event_statistics_orig

        # Jitter
        for(i in 1:length(jitter_vars)){ 

            # Duration/startyear require special treatment
            if(jitter_vars[i] == 'duration'){

                # We must jitter both startyear and duration, while ensuring
                # that 'event start time' + (duration + duration_gap_hours -
                # duration_offset_hours) does not overtake the next event start
                # time -- since that was a key feature of our event definition,
                # which is used in fitting the storm timing model [i.e. we know
                # that events were merged if the 'gap' was less than
                # duration_gap_hours, and must preserve that]
                #
                # For our data, such invalid jitter values are rare, but it can
                # happen with an unlucky jitter. In that case we just generate
                # a new random jitter until there are no invalid cases.
                #
                # Since it is such a rare event with our data, we do not expect
                # bias with this approach [as only a few values will be affected]
                #
                is_not_consistent = TRUE
                n = length(event_statistics[,1]) 
                jittered_duration = event_statistics$duration * NA
                jittered_startyear = event_statistics$startyear * NA
                ci = 1:n # Indices to change
                while(is_not_consistent){

                    # Trial new values of duration, startyear, endyear
                    jittered_duration[ci] = jitter(event_statistics$duration[ci], 
                        amount=jitter_amounts[i])
                    jittered_startyear[ci] = jitter(event_statistics$startyear[ci], 
                        amount=jitter_amounts[i]/year2hours)

                    new_endyear = jittered_startyear + jittered_duration/year2hours
                  
                    # Find values where the jitter 
                    ci = which(new_endyear[1:(n-1)] + 
                        (duration_gap_hours - duration_offset_hours)/year2hours > 
                        jittered_startyear[2:n])

                    if(length(ci) == 0) is_not_consistent = FALSE
                }

                event_statistics$duration = jittered_duration
                event_statistics$startyear = jittered_startyear
                event_statistics$endyear = new_endyear
    
            }else{

                # Beware -- for compatibility with S, the 'jitter' function
                # interprets the argument 'amount = 0' as 'use the default
                # jitter', rather than 'give no jitter' which would seem more
                # intuitive. Hence, we do not call the function if
                # jitter_amounts[i] == 0
                if(jitter_amounts[i] > 0){
                    event_statistics[[jitter_vars[i]]] = 
                        jitter(event_statistics[[jitter_vars[i]]], 
                            amount = jitter_amounts[i])
                }
            }
        }

        return(event_statistics)
    }
    return(jitter_event_statistics_function)
}

# Function that will return a jitter of the original event_statistics
event_statistics_orig = event_statistics

# Make function which can return jitter event statistics
jitter_event_statistics_function = make_jitter_event_statistics_function(
    event_statistics_orig,
    default_jitter_vars,
    default_jitter_amounts
    )

if(break_ties_with_jitter){
    # Jitter the event statistics
    event_statistics = jitter_event_statistics_function()
    summary(event_statistics)
}
```



# **Step 2: Compute the wave steepness**
----------------------------------------

In our analysis, we choose to model the wave steepness (a function of the wave
height and period) instead of working directly with wave period tp1 (since in
initial investigations, this led to better performance of later stages of the
statistical modelling). Below the wave steepness is computed using the Airy
wave dispersion relation, assuming the water depth is 80m (which is appropriate
for the MHL wave buoys).

```r
wavedisp = new.env()
source('../../R/wave_dispersion/wave_dispersion_relation.R', local=wavedisp)
buoy_depth = 80 # Depth in m. 
wavelengths = wavedisp$airy_wavelength(period=event_statistics$tp1, h=buoy_depth)
event_statistics$steepness = event_statistics$hsig/wavelengths
rm(wavelengths)
summary(event_statistics$steepness)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.008328 0.015950 0.019660 0.020650 0.024470 0.051790
```

```r
hist(event_statistics$steepness, main='Histogram of computed wave steepness values')
```

![plot of chunk computeWaveSteepness](figure/computeWaveSteepness-1.png)

# **Step 3: Exploration of the annual number of storms, the time between storms, and seasonal patterns**
--------------------------------------------------------------------------------------------------------

**Here we plot the number of events each year in the data.** Notice
how there are few events in the first and last year. This reflects
that our data is incomplete during those years.

```r
events_per_year = table(format(event_statistics$time, '%Y'))
# Hack into numeric type
year = as.numeric(names(events_per_year))
events_per_year = as.numeric(events_per_year)

# Plot
plot(year, as.numeric(events_per_year), t='h', lwd=3, lend=1, 
    main='Number of storm events per year',
    xlab='Year',
    ylab='Events per year', 
    ylim=c(0, max(events_per_year)))
```

![plot of chunk temporal_spacing1](figure/temporal_spacing1-1.png)

```r
# Clean up
rm(year)
```

Before considering detailed modelling of the time-series, we informally check
whether the annual counts behave like a Poisson distribution. This would
imply the annual count data has a mean that is equal to its variance (though
because of finite sampling, this is not expected to hold exactly). We check
this by simulating a large number of samples from a Poisson distribution with
the same mean as our data, and computing their variance. The first and last years
of the data are removed to avoid the artefacts mentioned above.

If the data were truely Poisson, then we would expect the data variance to fall
well within the samples from the simulation -- which it does. **The code below implements
this check.**

```r
l = length(events_per_year)
events_per_year_truncated = events_per_year[-c(1,l)]

# For Poisson, mean = variance (within sampling variability)
sample_mean = mean(events_per_year_truncated)
sample_mean
```

```
## [1] 21.96667
```

```r
sample_var = var(events_per_year_truncated)
sample_var
```

```
## [1] 15.41264
```

```r
# Simulate
n = length(events_per_year_truncated)
nsim = 10000
simulated_variance = replicate(nsim, var( rpois( n, lambda=sample_mean)))
empirical_distribution = ecdf(simulated_variance)

# What fraction of the empirical samples have a variance less than our data sample?
empirical_distribution(sample_var)
```

```
## [1] 0.1191
```

```r
hist(simulated_variance, breaks=60, freq=FALSE,
    main=paste0('Distribution of the sample variance for count data from a \n',
                'Poisson distribution (with the same mean and sample size as our data) \n',
                ' (Red line is our data variance)'))
abline(v=sample_var, col='red')
```

![plot of chunk checkpoisson1](figure/checkpoisson1-1.png)

```r
# Clean up
rm(l, n, nsim, simulated_variance, sample_var, empirical_distribution)
```

The above graphical and statistical checks do not suggests any particularly
strong deviation from the Poisson model *for the annual count data*. **Below, we
examine the data on shorter time-scales by plotting the distribution of times
between events, and the number of events each season.**

```r
num_events = length(event_statistics$startyear)
time_between_events = event_statistics$startyear[2:num_events] - 
                      event_statistics$endyear[1:(num_events-1)]

par(mfrow=c(2,1))
hist(time_between_events, freq=FALSE, xlab='Time between events (years)', 
    breaks=40, main='Histogram of time between events', cex.main=1.5)

# Add an exponential decay curve with rate based on the mean (corresponding to
# ideal Poisson Process)
xs = seq(0,1,len=100)
points(xs, exp(-xs/mean(time_between_events))/mean(time_between_events), 
    t='l', col=2)
grid()

# Compute the fraction of events which occur in each month
events_per_month = aggregate(
    rep(1/num_events, num_events), 
    list(month=format(event_statistics$time, '%m')), 
    FUN=sum)

# Get month labels as Jan, Feb, Mar, ...
month_label = month.abb[as.numeric(events_per_month$month)]
barplot(events_per_month$x, names.arg = month_label, 
        main='Fraction of events occuring in each calendar month',
        cex.main=1.5)
grid()
```

![plot of chunk timebetween_seasonal](figure/timebetween_seasonal-1.png)

```r
# Clean up
rm(time_between_events, xs, events_per_month, month_label, num_events)
```


# **Step 4: Modelling the storm event timings as a non-homogeneous Poisson process**
------------------------------------------------------------------------------------

Below we use the nhpoisp.R script to fit various non-homogeneous Poisson
process models to the storm event timings. The code below is somewhat complex,
since it automates the fit of a range of different models. **If you are trying
to learn to fit these models, it is strongly suggested you consult the tests
and introductory illustrations contained in the
folder[../../R/nhpoisp/](../../R/nhpoisp/), rather than starting with the code
below.**

```r
nhp = new.env()
source('../../R/nhpoisp/nhpoisp.R', local=nhp)

#
# Prepare the data for modelling
#

# Get event start time in years
event_time = event_statistics$startyear
# Get event duration in years 
#
# NOTE: Because we demand a 'gap' between events, it makes sense to add 
#       (duration_gap_hours - duration_offset_hours) to the
#       event duration, since **by our definitions** no other event is possible
#       in this time.
#       The term 'duration_offset_hours' occurs since we use inclusive counting
#       when deriving initial storm duratins [e.g. this is related to the fact
#       that if there are 2 hsig observations above the threshold, then the
#       storm duration is taken as 2 hours, whereas the time difference between
#       these observations is only 1 hour.]
#

event_duration_years = (
    event_statistics$duration + 
    duration_gap_hours - duration_offset_hours
    )/year2hours 
obs_start_time = DU$time_to_year(obs_start_time_strptime)


# We cannot use all the data for these fits which include SOI terms,
# because we don't have annual SOI for 2016+
bulk_fit_indices = which(event_time < 2016) 

#
# Define the 'annual' component of all equations we test. 
#
annual_rate_equations = list(
    # The most basic model. -- need to keep reference to 't' so it vectorises
    # Since there are about 22 events/year, a reasonable order-of-magnitude
    # starting guess for theta[1] is '30'. Using a value somewhat higher than
    # '22' also reduces the chance of lambda being negative (and clipped to zero)
    # at the starding parameters, which makes fitting numerically easier
    constant=list(eqn = 'theta[1] + 0*t', npar = 1, start_par = 30, par_scale=1),

    # A model where soi matters. Because we only have annual soi until end of
    # 2015, we make the function loop over values from 1985 - 2015 inclusive
    # (31 years inclusive).
    # No values are 'looped' for the actual model fit so the parameters are
    # unaffected - however, for plotting the loop helps, since we have to
    # simulate many values.
    #
    soi = list(
        eqn='theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))',
        npar=2, start_par = c(30, 0.1), par_scale=c(1,1/10))
    #
    # NOTE: Later, if simulating a series with synthetic soiA data, we will
    # have to change this equation so it no-longer 'loops' over only 31 years.
    # But for fitting, this is most convenient.
)

#
# Define the 'seasonal' component of all equations we test. 
#
seasonal_rate_equations = list(
    # The most basic model,
    constant=list(eqn='0', npar=0, start_par=c(), par_scale=c()),

    # Simple sinusoid,
    single_freq = list(
        eqn='theta[n_annual_par + 1]*sin(2*pi*(t - theta[n_annual_par + 2]))', 
        npar=2, start_par=c(5, 0.1), par_scale = c(1, 1/10)),
        
    # Double sinusoid
    double_freq = list(eqn=paste0(
        'theta[n_annual_par + 1]*sin(2*pi*(t - theta[n_annual_par + 2])) + ',
        'theta[n_annual_par + 3]*sin(4*pi*(t - theta[n_annual_par + 4]))'),
        npar=4, start_par=c(5, 0.1, 1, 0.1), par_scale=c(1, 1/10, 1, 1/10)),

    # Sawtooth
    sawtooth = list(
        eqn='theta[n_annual_par + 1]*abs(2/pi*asin(cos(pi*(t-theta[n_annual_par+2]))))',
        npar=2, start_par=c(5, 0.5), par_scale=c(1, 1/10))
)

#
# Define the 'cluster' component of all equations we test
#
cluster_rate_equations = list(
    # The most basic model
    constant=list(eqn="0", npar=0, start_par=c(), par_scale=c()),

    # A clustering model
    # Note the cluster time-scale must be positive to keep it sensible
    cluster = list(
        eqn=paste0('theta[n_annual_par + n_seasonal_par + 1]*',
            'exp((tlast - t)*abs(theta[n_annual_par + n_seasonal_par + 2]))'),
        npar=2, start_par=c(1, 100), par_scale=c(1, 100))
)

# We will loop over all combinations of the above equations
annual_rate_names = names(annual_rate_equations)
seasonal_rate_names = names(seasonal_rate_equations)
cluster_rate_names = names(cluster_rate_equations)

# Prepare for loop -- use these variables to store results
counter = 0
exhaustive_lambda_model_fits = list()

# Options to control the loop over all models
print_info = TRUE 
fit_model = TRUE # If FALSE, just test we can construct the eqn's correctly
use_previous_fit_for_cluster_pars = FALSE # Trick to improve starting parameter guess for models with clustering

# Loop over all rate equations (representing all combinations of
# annual/seasona/cluster models)
for(ar_name in annual_rate_names){
    for(sr_name in seasonal_rate_names){
        for(cr_name in cluster_rate_names){

            counter = counter + 1
            # Make starting parameters
            if( (cr_name == 'constant') | 
                (!use_previous_fit_for_cluster_pars) | 
                !fit_model ){
                start_par = c(annual_rate_equations[[ar_name]]$start_par, 
                              seasonal_rate_equations[[sr_name]]$start_par,
                              cluster_rate_equations[[cr_name]]$start_par)
            }else{
                # Better to use the parameters from before
                start_par = c(local_fit$par, cluster_rate_equations[[cr_name]]$start_par)
            }

            # Make scaling parameters for model parameters -- having these a similar magnitude
            # to the true value of the model parameters can help with convergence.
            par_scale = c(annual_rate_equations[[ar_name]]$par_scale, 
                          seasonal_rate_equations[[sr_name]]$par_scale,
                          cluster_rate_equations[[cr_name]]$par_scale)

            # Make preliminary equation
            rate_equation = paste0(
                annual_rate_equations[[ar_name]]$eqn, '+', 
                seasonal_rate_equations[[sr_name]]$eqn, '+', 
                cluster_rate_equations[[cr_name]]$eqn)
            
            # Sub in the required parameters 
            rate_equation = gsub('n_annual_par', annual_rate_equations[[ar_name]]$npar, 
                rate_equation)
            rate_equation = gsub('n_seasonal_par', seasonal_rate_equations[[sr_name]]$npar, 
                rate_equation)
           
            # Print out information 
            if(print_info){
                print('')
                print('')
                print(c('Annual model: ', ar_name))
                print(c('Seasonal model: ', sr_name))
                print(c('Clustering model: ', cr_name))
                print(c('Rate equation: ', rate_equation))
                print(c('Starting par: ', start_par))
                print(c('parscale for optim: ', par_scale))
            }

            if(fit_model){

                # For all single parameter models, just use BFGS. For others,
                # do two Nelder-Mead fits, followed by BFGS
                if(length(start_par) == 1){
                    optim_method_sequence = 'BFGS'
                }else{
                    optim_method_sequence = c('Nelder-Mead', 'Nelder-Mead', 'BFGS')
                }
            
                # This is the main fitting routine    
                local_fit =  nhp$fit_nhpoisp(event_time[bulk_fit_indices],
                    rate_equation=rate_equation,
                    minimum_rate=0.0,
                    initial_theta=start_par,
                    x0 = obs_start_time,
                    event_durations = event_duration_years[bulk_fit_indices],
                    number_of_passes = length(optim_method_sequence),
                    optim_method=optim_method_sequence,
                    enforce_nonnegative_theta=FALSE,
                    optim_control=list(parscale = par_scale),
                    use_optim2=FALSE,
                    use_numDeriv_hessian=TRUE)

                # Store the result in a list
                exhaustive_lambda_model_fits[[counter]] = local_fit  
        
                
                if(print_info) {
                    print('...Fit...')
                    print(c('...parameters...:', local_fit$par))
                    print(c('...standard errors (approximate)...:', nhp$get_fit_standard_errors(local_fit)))
                    print(c('...convergence flag...:', local_fit$convergence))
                    print(c('...negative log likelihood...:', local_fit$value))
                }
            }
        }
    }
}
```

```
## [1] ""
## [1] ""
## [1] "Annual model: " "constant"      
## [1] "Seasonal model: " "constant"        
## [1] "Clustering model: " "constant"          
## [1] "Rate equation: "    "theta[1] + 0*t+0+0"
## [1] "Starting par: " "30"            
## [1] "parscale for optim: " "1"                   
## [1] "...Fit..."
## [1] "...parameters...:" "24.9048999590305" 
## [1] "...standard errors (approximate)...:"
## [2] "0.967226167176675"                   
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1468.59276230447"             
## [1] ""
## [1] ""
## [1] "Annual model: " "constant"      
## [1] "Seasonal model: " "constant"        
## [1] "Clustering model: " "cluster"           
## [1] "Rate equation: "                                                         
## [2] "theta[1] + 0*t+0+theta[1 + 0 + 1]*exp((tlast - t)*abs(theta[1 + 0 + 2]))"
## [1] "Starting par: " "30"             "1"              "100"           
## [1] "parscale for optim: " "1"                    "1"                   
## [4] "100"                 
## [1] "...Fit..."
## [1] "...parameters...:" "24.6613469336676"  "4.41412134306482" 
## [4] "423.176142248569" 
## [1] "...standard errors (approximate)...:"
## [2] "1.04283094769225"                    
## [3] "7.42070768339896"                    
## [4] "774.483896813472"                    
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1468.84516691669"             
## [1] ""
## [1] ""
## [1] "Annual model: " "constant"      
## [1] "Seasonal model: " "single_freq"     
## [1] "Clustering model: " "constant"          
## [1] "Rate equation: "                                           
## [2] "theta[1] + 0*t+theta[1 + 1]*sin(2*pi*(t - theta[1 + 2]))+0"
## [1] "Starting par: " "30"             "5"              "0.1"           
## [1] "parscale for optim: " "1"                    "1"                   
## [4] "0.1"                 
## [1] "...Fit..."
## [1] "...parameters...:" "25.0999431159312"  "7.01586446990227" 
## [4] "0.245997453398378"
## [1] "...standard errors (approximate)...:"
## [2] "0.975538253586096"                   
## [3] "1.36346577413517"                    
## [4] "0.0308062726939097"                  
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1481.84016810313"             
## [1] ""
## [1] ""
## [1] "Annual model: " "constant"      
## [1] "Seasonal model: " "single_freq"     
## [1] "Clustering model: " "cluster"           
## [1] "Rate equation: "                                                                                                 
## [2] "theta[1] + 0*t+theta[1 + 1]*sin(2*pi*(t - theta[1 + 2]))+theta[1 + 2 + 1]*exp((tlast - t)*abs(theta[1 + 2 + 2]))"
## [1] "Starting par: " "30"             "5"              "0.1"           
## [5] "1"              "100"           
## [1] "parscale for optim: " "1"                    "1"                   
## [4] "0.1"                  "1"                    "100"                 
## [1] "...Fit..."
## [1] "...parameters...:" "24.9897360849056"  "6.98730114009642" 
## [4] "0.24702119748545"  "2.55975915634248"  "559.524701716886" 
## [1] "...standard errors (approximate)...:"
## [2] "1.02765530339647"                    
## [3] "1.36679664489064"                    
## [4] "0.0310246243838514"                  
## [5] "7.73170221121681"                    
## [6] "1433.08930713862"                    
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1481.90629337346"             
## [1] ""
## [1] ""
## [1] "Annual model: " "constant"      
## [1] "Seasonal model: " "double_freq"     
## [1] "Clustering model: " "constant"          
## [1] "Rate equation: "                                                                                       
## [2] "theta[1] + 0*t+theta[1 + 1]*sin(2*pi*(t - theta[1 + 2])) + theta[1 + 3]*sin(4*pi*(t - theta[1 + 4]))+0"
## [1] "Starting par: " "30"             "5"              "0.1"           
## [5] "1"              "0.1"           
## [1] "parscale for optim: " "1"                    "1"                   
## [4] "0.1"                  "1"                    "0.1"                 
## [1] "...Fit..."
## [1] "...parameters...:" "25.0989162319808"  "7.10217832415091" 
## [4] "0.241673685824822" "0.9490501440895"   "0.269307965925712"
## [1] "...standard errors (approximate)...:"
## [2] "0.975554773049773"                   
## [3] "1.3871072089424"                     
## [4] "0.0308159104604589"                  
## [5] "1.35446448299282"                    
## [6] "0.114938563398429"                   
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1482.08488371701"             
## [1] ""
## [1] ""
## [1] "Annual model: " "constant"      
## [1] "Seasonal model: " "double_freq"     
## [1] "Clustering model: " "cluster"           
## [1] "Rate equation: "                                                                                                                                             
## [2] "theta[1] + 0*t+theta[1 + 1]*sin(2*pi*(t - theta[1 + 2])) + theta[1 + 3]*sin(4*pi*(t - theta[1 + 4]))+theta[1 + 4 + 1]*exp((tlast - t)*abs(theta[1 + 4 + 2]))"
## [1] "Starting par: " "30"             "5"              "0.1"           
## [5] "1"              "0.1"            "1"              "100"           
## [1] "parscale for optim: " "1"                    "1"                   
## [4] "0.1"                  "1"                    "0.1"                 
## [7] "1"                    "100"                 
## [1] "...Fit..."
## [1] "...parameters...:" "24.980080414316"   "7.07039402049601" 
## [4] "0.242807238234928" "0.951835354721623" "0.267203707534672"
## [7] "2.59708291922587"  "546.088744253807" 
## [1] "...standard errors (approximate)...:"
## [2] "1.02815259350692"                    
## [3] "1.38962246498701"                    
## [4] "0.0310560492122841"                  
## [5] "1.35428272906393"                    
## [6] "0.114269177753343"                   
## [7] "7.64231231426186"                    
## [8] "1337.47400762545"                    
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1482.1571144146"              
## [1] ""
## [1] ""
## [1] "Annual model: " "constant"      
## [1] "Seasonal model: " "sawtooth"        
## [1] "Clustering model: " "constant"          
## [1] "Rate equation: "                                                     
## [2] "theta[1] + 0*t+theta[1 + 1]*abs(2/pi*asin(cos(pi*(t-theta[1+2]))))+0"
## [1] "Starting par: " "30"             "5"              "0.5"           
## [1] "parscale for optim: " "1"                    "1"                   
## [4] "0.1"                 
## [1] "...Fit..."
## [1] "...parameters...:" "16.3051496023983"  "17.5858795883073" 
## [4] "0.51952099068441" 
## [1] "...standard errors (approximate)...:"
## [2] "1.7197909963374"                     
## [3] "3.33951286834653"                    
## [4] "0.0313336301636397"                  
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1482.28086747763"             
## [1] ""
## [1] ""
## [1] "Annual model: " "constant"      
## [1] "Seasonal model: " "sawtooth"        
## [1] "Clustering model: " "cluster"           
## [1] "Rate equation: "                                                                                                           
## [2] "theta[1] + 0*t+theta[1 + 1]*abs(2/pi*asin(cos(pi*(t-theta[1+2]))))+theta[1 + 2 + 1]*exp((tlast - t)*abs(theta[1 + 2 + 2]))"
## [1] "Starting par: " "30"             "5"              "0.5"           
## [5] "1"              "100"           
## [1] "parscale for optim: " "1"                    "1"                   
## [4] "0.1"                  "1"                    "100"                 
## [1] "...Fit..."
## [1] "...parameters...:" "16.2217674454589"  "17.4852169737216" 
## [4] "0.519815087723505" "2.9338814913764"   "-540.3761287584"  
## [1] "...standard errors (approximate)...:"
## [2] "1.73936704610301"                    
## [3] "3.34560165474652"                    
## [4] "0.0316537126717003"                  
## [5] "7.6426887942624"                     
## [6] "1217.11270369906"                    
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1482.37362320853"             
## [1] ""
## [1] ""
## [1] "Annual model: " "soi"           
## [1] "Seasonal model: " "constant"        
## [1] "Clustering model: " "constant"          
## [1] "Rate equation: "                                                        
## [2] "theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+0+0"
## [1] "Starting par: " "30"             "0.1"           
## [1] "parscale for optim: " "1"                    "0.1"                 
## [1] "...Fit..."
## [1] "...parameters...:" "25.1669932097885"  "0.22106465153444" 
## [1] "...standard errors (approximate)...:"
## [2] "0.989635621802474"                   
## [3] "0.13141875041298"                    
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1470.01268662021"             
## [1] ""
## [1] ""
## [1] "Annual model: " "soi"           
## [1] "Seasonal model: " "constant"        
## [1] "Clustering model: " "cluster"           
## [1] "Rate equation: "                                                                                                              
## [2] "theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+0+theta[2 + 0 + 1]*exp((tlast - t)*abs(theta[2 + 0 + 2]))"
## [1] "Starting par: " "30"             "0.1"            "1"             
## [5] "100"           
## [1] "parscale for optim: " "1"                    "0.1"                 
## [4] "1"                    "100"                 
## [1] "...Fit..."
## [1] "...parameters...:" "24.9769043872307"  "0.215120622332752"
## [4] "3.81169538612404"  "491.32178867225"  
## [1] "...standard errors (approximate)...:"
## [2] "1.05617668512625"                    
## [3] "0.131926023092234"                   
## [4] "7.70759065377083"                    
## [5] "990.863281049692"                    
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1470.17797609399"             
## [1] ""
## [1] ""
## [1] "Annual model: " "soi"           
## [1] "Seasonal model: " "single_freq"     
## [1] "Clustering model: " "constant"          
## [1] "Rate equation: "                                                                                                
## [2] "theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*sin(2*pi*(t - theta[2 + 2]))+0"
## [1] "Starting par: " "30"             "0.1"            "5"             
## [5] "0.1"           
## [1] "parscale for optim: " "1"                    "0.1"                 
## [4] "1"                    "0.1"                 
## [1] "...Fit..."
## [1] "...parameters...:" "25.3489207315751"  "0.212002716299742"
## [4] "7.01293258463332"  "0.2482775625759"  
## [1] "...standard errors (approximate)...:"
## [2] "0.996720719965785"                   
## [3] "0.129143393364926"                   
## [4] "1.36706039239213"                    
## [5] "0.03062895070311"                    
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1483.19164488692"             
## [1] ""
## [1] ""
## [1] "Annual model: " "soi"           
## [1] "Seasonal model: " "single_freq"     
## [1] "Clustering model: " "cluster"           
## [1] "Rate equation: "                                                                                                                                                      
## [2] "theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*sin(2*pi*(t - theta[2 + 2]))+theta[2 + 2 + 1]*exp((tlast - t)*abs(theta[2 + 2 + 2]))"
## [1] "Starting par: " "30"             "0.1"            "5"             
## [5] "0.1"            "1"              "100"           
## [1] "parscale for optim: " "1"                    "0.1"                 
## [4] "1"                    "0.1"                  "1"                   
## [7] "100"                 
## [1] "...Fit..."
## [1] "...parameters...:" "25.2644702465053"  "0.209946195312091"
## [4] "6.99657251831648"  "0.248968755217749" "2.09383406005446" 
## [7] "632.277024205568" 
## [1] "...standard errors (approximate)...:"
## [2] "1.0440488103676"                     
## [3] "0.129304314435757"                   
## [4] "1.36919792315664"                    
## [5] "0.0308279357887612"                  
## [6] "8.10398074296041"                    
## [7] "2042.44995605656"                    
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1483.23248507946"             
## [1] ""
## [1] ""
## [1] "Annual model: " "soi"           
## [1] "Seasonal model: " "double_freq"     
## [1] "Clustering model: " "constant"          
## [1] "Rate equation: "                                                                                                                                            
## [2] "theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*sin(2*pi*(t - theta[2 + 2])) + theta[2 + 3]*sin(4*pi*(t - theta[2 + 4]))+0"
## [1] "Starting par: " "30"             "0.1"            "5"             
## [5] "0.1"            "1"              "0.1"           
## [1] "parscale for optim: " "1"                    "0.1"                 
## [4] "1"                    "0.1"                  "1"                   
## [7] "0.1"                 
## [1] "...Fit..."
## [1] "...parameters...:" "25.3439463185761"  "0.210528456409697"
## [4] "7.07606048708399"  "0.244114579107015" "0.895572386217896"
## [7] "0.252003942604594"
## [1] "...standard errors (approximate)...:"
## [2] "0.996575445329289"                   
## [3] "0.129711185757779"                   
## [4] "1.38424749403921"                    
## [5] "0.0309362771125661"                  
## [6] "1.36159605302444"                    
## [7] "0.120852927894629"                   
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1483.4095571846"              
## [1] ""
## [1] ""
## [1] "Annual model: " "soi"           
## [1] "Seasonal model: " "double_freq"     
## [1] "Clustering model: " "cluster"           
## [1] "Rate equation: "                                                                                                                                                                                                  
## [2] "theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*sin(2*pi*(t - theta[2 + 2])) + theta[2 + 3]*sin(4*pi*(t - theta[2 + 4]))+theta[2 + 4 + 1]*exp((tlast - t)*abs(theta[2 + 4 + 2]))"
## [1] "Starting par: " "30"             "0.1"            "5"             
## [5] "0.1"            "1"              "0.1"            "1"             
## [9] "100"           
## [1] "parscale for optim: " "1"                    "0.1"                 
## [4] "1"                    "0.1"                  "1"                   
## [7] "0.1"                  "1"                    "100"                 
## [1] "...Fit..."
## [1] "...parameters...:" "25.2497291289696"  "0.209879588470234"
## [4] "7.03634095980815"  "0.245360632915104" "0.91392475041329" 
## [7] "0.250496414143282" "2.24059816548656"  "634.211751243086" 
## [1] "...standard errors (approximate)...:"
## [2] "1.04610587854082"                    
## [3] "0.129839598148953"                   
## [4] "1.38581180220226"                    
## [5] "0.031187507834233"                   
## [6] "1.36146420680171"                    
## [7] "0.119187414092851"                   
## [8] "8.02263485840419"                    
## [9] "1977.71527031315"                    
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1483.45654995779"             
## [1] ""
## [1] ""
## [1] "Annual model: " "soi"           
## [1] "Seasonal model: " "sawtooth"        
## [1] "Clustering model: " "constant"          
## [1] "Rate equation: "                                                                                                          
## [2] "theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+0"
## [1] "Starting par: " "30"             "0.1"            "5"             
## [5] "0.5"           
## [1] "parscale for optim: " "1"                    "0.1"                 
## [4] "1"                    "0.1"                 
## [1] "...Fit..."
## [1] "...parameters...:" "16.5435352722584"  "0.219548780906959"
## [4] "17.6198233434521"  "0.523860117374763"
## [1] "...standard errors (approximate)...:"
## [2] "1.73546515299518"                    
## [3] "0.12895041841605"                    
## [4] "3.34550076308239"                    
## [5] "0.0237126874889317"                  
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1483.7273022572"              
## [1] ""
## [1] ""
## [1] "Annual model: " "soi"           
## [1] "Seasonal model: " "sawtooth"        
## [1] "Clustering model: " "cluster"           
## [1] "Rate equation: "                                                                                                                                                                
## [2] "theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+theta[2 + 2 + 1]*exp((tlast - t)*abs(theta[2 + 2 + 2]))"
## [1] "Starting par: " "30"             "0.1"            "5"             
## [5] "0.5"            "1"              "100"           
## [1] "parscale for optim: " "1"                    "0.1"                 
## [4] "1"                    "0.1"                  "1"                   
## [7] "100"                 
## [1] "...Fit..."
## [1] "...parameters...:" "16.5046739024132"  "0.219884464274771"
## [4] "17.613983456309"   "0.523982434879071" "6.16523842934303" 
## [7] "-3633.94742358632"
## [1] "...standard errors (approximate)...:"
## [2] "1.73734670754271"                    
## [3] "0.128869826949164"                   
## [4] "3.34330749479561"                    
## [5] "0.0240212865457841"                  
## [6] "19.1687109733373"                    
## [7] "8073.72603430121"                    
## [1] "...convergence flag...:" "0"                      
## [1] "...negative log likelihood...:" "-1483.79632661922"
```
* For every model, the above print-outs show the best fit parameters, approximate
standard errors, convergence flags (0 indicates convergence), and other information.

* For our (unperturbed) data, all models should converge.

* For models with clustering, we should see that the clustering terms always
have approximate standard errors which are larger than the parameter estimates,
indicating that those parameters are not well constrained. More advanced methods
(e.g. profile likelihood) would be required to precicely quantify the uncertainties,
but for our purposes this is not required.

* With different data, it is possible for models not to converge. If that
happens, then it is necessary to try different starting values and/or different
optimization options. Even if the optimization algorithm reports a `convergence
flag` of `0` (which indicates successful convergence), it is possible that a
local (but not global) optimum was found. Graphical checks of the model fit
(such as diagnostics shown below) can help identify situations when the fit is
poor, in which case both numerical convergence and the suitability of the
lambda model itself should be checked.

* With other data, it is also possible for the model to appear to converge, but
for the numerically computed Hessian of the likelihood at the optimum to be
non-positive definite. The latter situation either means that the model has not
really converged, or that minor numerical errors in the computation of the
Hessian have made it non-positive-definite. In this instance the code will
print a warning, and compute approximate standard errors using a Hessian matrix
which is forced to be positive definite. However, it is highly advisable to
double check the fit if the model is to be used subsequently.

**Next we compute the AIC for each model fit above, in order to select the most
parsimonious model.** This provides a means of ranking the models, based on
their parsimony. We use the 'corrected AIC' which slightly adjusts the standard
AIC for sample size. This correction has no impact on the model rankings here
however, because we have a relatively large sample compared with the number of
model parameters, which means the correction is always small.

```r
# Quick check that we have results for every model [e.g. convergence failures could break this]
expected_num_models = length(annual_rate_names) * length(seasonal_rate_names) * length(cluster_rate_names)
if(counter != expected_num_models){
    print('WARNING: Some storm timing model fits failed')
}

# corrected AIC for every model trialled above 
exhaustive_AICs =  unlist(lapply(exhaustive_lambda_model_fits, 
    f<-function(x) nhp$compute_fit_AIC_BIC(x, correct_AIC=TRUE)$AIC))

# Choose most parsimonious (i.e. model having min AIC)
best_nhp_model = exhaustive_lambda_model_fits[[which.min(exhaustive_AICs)]]
# What is the equation of the best model?
print(best_nhp_model$rate_equation)
```

```
## [1] "theta[1] + theta[2]*CI_annual_fun$soi(floor((t - 1985)%%31 + 1985))+theta[2 + 1]*abs(2/pi*asin(cos(pi*(t-theta[2+2]))))+0"
```

```r
# Get lambda function from the best fit
lambda = nhp$get_lambda_function(
    best_nhp_model$par, 
    rate_equation=best_nhp_model$rate_equation, 
    minimum_rate=0.)

# Make a plot comparing the model and data -- integrating over all observed SOI
# The plot scripts also compare the empirical distribution of the data and 
# simulations from the fitted model, using a number of KS-test based statistics
nhp$plot_nhpoisson_diagnostics(event_time[bulk_fit_indices], 
    event_duration_years[bulk_fit_indices], lambda, nbins=25)
```

```
## Loading required package: Matching
```

```
## Loading required package: MASS
```

```
## ## 
## ##  Matching (Version 4.9-2, Build Date: 2015-12-25)
## ##  See http://sekhon.berkeley.edu/matching for additional documentation.
## ##  Please cite software as:
## ##   Jasjeet S. Sekhon. 2011. ``Multivariate and Propensity Score Matching
## ##   Software with Automated Balance Optimization: The Matching package for R.''
## ##   Journal of Statistical Software, 42(7): 1-52. 
## ##
```

```
## [1] "KS TEST OF THE EVENTS TIME-OF-YEAR"
## $ks.boot.pvalue
## [1] 0.908
## 
## $ks
## 
## 	Two-sample Kolmogorov-Smirnov test
## 
## data:  Tr and Co
## D = 0.021498, p-value = 0.922
## alternative hypothesis: two-sided
## 
## 
## $nboots
## [1] 1000
## 
## attr(,"class")
## [1] "ks.boot"
```

```
## [1] "KS TEST OF THE TIME BETWEEN EVENTS"
## $ks.boot.pvalue
## [1] 0.957
## 
## $ks
## 
## 	Two-sample Kolmogorov-Smirnov test
## 
## data:  Tr and Co
## D = 0.019428, p-value = 0.9656
## alternative hypothesis: two-sided
## 
## 
## $nboots
## [1] 1000
## 
## attr(,"class")
## [1] "ks.boot"
## [1] "KS TEST OF THE NUMBER OF EVENTS EACH YEAR"
## $ks.boot.pvalue
## [1] 0.782
## 
## $ks
## 
## 	Two-sample Kolmogorov-Smirnov test
## 
## data:  Tr and Co
## D = 0.090558, p-value = 0.9629
## alternative hypothesis: two-sided
## 
## 
## $nboots
## [1] 1000
## 
## attr(,"class")
## [1] "ks.boot"
```

![plot of chunk computeCorrectedAIC](figure/computeCorrectedAIC-1.png)


**Save output for use later**

```r
dir.create('Rimages', showWarnings=FALSE)
Rimage_title = paste0('Rimages/session_storm_timings_', run_title_id, '.Rdata')
save.image(Rimage_title)
```


## **Moving On**
The next steps of the vignette begin at
[statistical_model_univariate_distributions.md](statistical_model_univariate_distributions.md).
