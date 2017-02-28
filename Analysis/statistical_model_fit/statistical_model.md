
# **Statistical modelling of the storm event dataset**
-------------------------------------------------------------

*Gareth Davies, Geoscience Australia 2017*

# Introduction
------------------

This document follows on from
[../preprocessing/extract_storm_events.md](../preprocessing/extract_storm_events.md)
in describing our statistical analysis of storm waves at Old Bar. 

It illustrates the process of fitting the statistical model to the data 

It is essential that the scripts in [../preprocessing](../preprocessing) have
alread been run, and produced an
RDS file *'../preprocessing/Derived_data/event_statistics.RDS'*. **To make sure, the
code below throws an error if the latter file does not exist.**

```r
# Check that the pre-requisites exist
if(!file.exists('../preprocessing/Derived_data/event_statistics.RDS')){
    stop('It appears you have not yet run all codes in ../preprocessing. They must be run before continuing')
}
```

Supposing the above did not generate any errors, and you have R installed,
along with all the packages required to run this code, and a copy of the
*stormwavecluster* git repository, then you should be able to re-run the
analysis here by simply copy-pasting the code. Alternatively, it can be run
with the `knit` command in the *knitr* package: 

```r
library(knitr)
knit('statistical_model.Rmd')
```

The basic approach followed here is to:
* **Step 1: Get the preprocessed data and address rounding artefacts**

Later we will use the statistical model to simulate synthetic storm event
time-series.

# **Step 1: Get the preprocessed data and address rounding artefacts**
-----------------------------------------------------------------


**Below we read

```r
# Get the data_utilities (in their own environment to keep the namespace clean)
DU = new.env()
source('../preprocessing/data_utilities.R', local=DU) 

# Read data saved by processing scripts
event_statistics_list = readRDS('../preprocessing/Derived_data/event_statistics.RDS')

# Extract variables that we need from previous analysis
for(varname in names(event_statistics_list)){
    assign(varname, event_statistics_list[[varname]])
}

# Clean up
rm(event_statistics_list)

# Definitions controlling the synthetic series creation
nyears_synthetic_series = 1e+06 #1e+03

# Length of each MCMC chain. Should be 'large' e.g 10^6, except for test runs 
# We run multiple chains to enhance the likelihood of detecting non-convergence
# since anyway this is cheap in parallel. These are pooled for final estimates,
# but it is essential to manually check the convergence of the chains [e.g.
# by comparing high return period confidence intervals].
mcmc_chain_length = 1e+06 #1e+05 
# To reduce the data size, we can throw away all but a fraction of the mcmc
# chains. This has computational (memory) benefits if the MCMC samples are
# strongly autocorrelated, but no other advantages.
mcmc_chain_thin = 20 

# Useful number to convert from years to hours (ignoring details of leap-years)
year2hours = 365.25*24

# Optionally remove ties in the event statistics by jittering
break_ties_with_jitter = FALSE

# Look at the variables we have
ls()
```

```
##  [1] "break_ties_with_jitter"           "CI_annual_fun"                   
##  [3] "data_duration_years"              "DU"                              
##  [5] "duration_gap_hours"               "duration_offset_hours"           
##  [7] "duration_threshold_hours"         "event_statistics"                
##  [9] "event_statistics_orig"            "hsig_threshold"                  
## [11] "jitter_event_statistics"          "make_jitter_event_statistics"    
## [13] "mcmc_chain_length"                "mcmc_chain_thin"                 
## [15] "nyears_synthetic_series"          "obs_start_time_strptime"         
## [17] "smooth_tideResid_fun_stl_monthly" "soi_SL_lm"                       
## [19] "varname"                          "year2hours"
```

If our event statistics are subject to rounding (introducing 'ties' or repeated
values into the data), then it is possible for some statistical methods used
here to perform badly (since they assume continuous data, which has probability
zero of ties).

**Below we optionally perturb the `event_statistics` to remove ties**. The perturbation
size is related to the resolution of the data, which is 1 cm for Hsig, 1 hour for
duration, and 1 degree for direction. For tp1 (which has the most ties) and only
40 unique values, the bins are irregularly spaced without an obvious pattern.
The median distance between unique tp1 values after sorting is 0.25, with a maximum of
1.06, and a minimum of 0.01. Below a uniform perturbation of 0.1 is applied.


```r
# Make a function which will return a jittered version of the original
# event_statistics
make_jitter_event_statistics<-function(event_statistics_orig){
    # Save original data
    event_statistics_orig = event_statistics_orig

    # Function that will jitter the original event_statistics
    jitter_event_statistics<-function(
        jitter_vars = c('hsig', 'duration', 'dir', 'tp1'),
        # Jitter amounts = half of the bin size [see comments above regarding TP1]
        jitter_amounts = c(0.005, 0.5, 0.5, 0.1)
        ){

        event_statistics = event_statistics_orig

        # Jitter
        for(i in 1:length(jitter_vars)){ 
            event_statistics[[jitter_vars[i]]] = 
                jitter(event_statistics[[jitter_vars[i]]], 
                    amount = jitter_amounts[i])
        }

        # But hsig must be above the threshold
        kk = which(event_statistics$hsig <= hsig_threshold)
        if(length(kk) > 0){
            event_statistics$hsig[kk] = event_statistics_orig$hsig[kk]
        }

        return(event_statistics)
    }
    return(jitter_event_statistics)
}

# Function that will return a jitter of the original event_statistics
event_statistics_orig = event_statistics
jitter_event_statistics = make_jitter_event_statistics(event_statistics_orig)

if(break_ties_with_jitter){
    # Jitter the event statistics
    event_statistics = jitter_event_statistics()
    summary(event_statistics)
}
```