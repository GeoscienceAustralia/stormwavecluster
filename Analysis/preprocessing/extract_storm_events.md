
# **Extract storm events and conduct an exploratory data analysis**
-------------------------------------------------------------

*Gareth Davies, Geoscience Australia 2017*

# Introduction
------------------

This document follows on from [preprocess_data.md](preprocess_data.md) in describing
our analysis of storm waves at Old Bar. It illustrates the extraction
of storm wave events from the "Old Bar" time-series created in the earlier script, and 
some preliminary analyses of the data.

It is essential that the [preprocess_data.md](preprocess_data.md) code has
already been run, and produced an RDS file
*'Rimages/Session_data_processing_clean_XXXX.Rdata'*, where XXXX is the
`run_title_id`.


Supposing the above did not generate any errors, and you have R installed,
along with all the packages required to run this code, and a copy of the
*stormwavecluster* git repository, then you should be able to re-run the
analysis here by simply copy-pasting the code. Alternatively, it can be run
with the `knit` command in the *knitr* package: 

```r
library(knitr)
knit('extract_storm_events.Rmd')
```

To run the code in tie-breaking mode, be sure to pass the a commandline
argument matching `break_ties` to R when starting, followed by an integer ID > 0,
e.g.

    R --args --break_ties 1234

or

    Rscript knit_preprocess_data.R --break_ties 1234

You must have already run the dependencies mentioned above with the same
commandline arguments.

The basic approach followed here is to:
* **Step 1**: Extract storm events from the "Old Bar" time-series created earlier
* **Step 2**: Compute summary statistics for each storm event
* **Step 3**: Study changes in monthly mean sea level, and remove seasonal and inter-annual trends from our tidal residual, to better estimate the storm surge component.
* **Step 4**: Perform bias correction of the wave directions obtained from stations other than Crowdy Head.
* **Step 5**: Do some exploratory analysis of the data

Later we will develop the statistical analysis of the storm events.


# **Step 0: Process commandline arguments**
-------------------------------------------

This is where we read the commandline arguments relevant for running with
perturbed data. See the corresponding section in
[preprocess_data.md](preprocess_data.md) for more information. 


```r
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

    # A dummy id for the run with no data perturbation
    session_n = 0

}

# Make a 'title' which can appear in filenames to identify this run
run_title_id = paste0(break_ties_with_jitter, '_', session_n)

print(c('run_title_id: ', run_title_id))
```

```
## [1] "run_title_id: " "FALSE_0"
```


# **Step 1: Extract storm events from the "Old Bar" time-series**
-----------------------------------------------------------------

Here we follow on from the analysis in [preprocess_data.md](preprocess_data.md)
by loading the associated R session.

```r
# This will only work if the 'preprocess_data' code has already been run
input_file = paste0('Rimages/Session_data_processing_clean_', run_title_id, '.Rdata')
if(!file.exists(input_file)){
    print('Looking for file: ')
    print(input_file)
    stop('Input file does not exist')
}else{
    load(input_file)
}
```

```
## [1] "Looking for file: "
## [1] "Rimages/Session_data_processing_clean_FALSE_0.Rdata"
```

```
## Error in eval(expr, envir, enclos): Input file does not exist
```

**Now, we extract storm events from the `full_data` time-series.** We initially
define storm events as periods in which the significant wave height exceeds its
95th percentile. We also merge storm events which are separated by less than 24
hours.

```r
# Event extraction
hsig_threshold = quantile(full_data$hsig, p=0.95, na.rm=T)
```

```
## Error in quantile(full_data$hsig, p = 0.95, na.rm = T): object 'full_data' not found
```

```r
hsig_threshold
```

```
## Error in eval(expr, envir, enclos): object 'hsig_threshold' not found
```

```r
duration_threshold_hours = 0 # Suggest 0 -- remove short events later if needed 
duration_gap_hours = 24 # Must always be >= data time spacing (1 in this case)

# See the routine extract_events in data_utilities.R for details
event_set = DU$extract_events(full_data, 
    hsig_threshold = hsig_threshold,
    duration_threshold_hours = duration_threshold_hours, 
    duration_gap_hours = duration_gap_hours,
    events_to_combine = NULL)
```

```
## Error in eval(expr, envir, enclos): object 'DU' not found
```

**Here we check how many storms were identified, and show an example plot of one**

```r
# event_set is a list containing the start/end index of each event, and another
# list called 'data' which contains the event data (a separate entry for each
# event)
names(event_set)
```

```
## Error in eval(expr, envir, enclos): object 'event_set' not found
```

```r
# How many events?
length(event_set$start_index)
```

```
## Error in eval(expr, envir, enclos): object 'event_set' not found
```

```r
# Start should always be <= end!
stopifnot(all(event_set$start_index <= event_set$end_index))
```

```
## Error in stopifnot(all(event_set$start_index <= event_set$end_index)): object 'event_set' not found
```

```r
# How many datasets do we have (should be 1 for each event)
num_events = length(event_set$data)
```

```
## Error in eval(expr, envir, enclos): object 'event_set' not found
```

```r
num_events
```

```
## Error in eval(expr, envir, enclos): object 'num_events' not found
```

```r
# Plot one event [ change the index to get a good one, number 442 was good for me ]
event_example = event_set$data[[442]]
```

```
## Error in eval(expr, envir, enclos): object 'event_set' not found
```

```r
par(mar=c(3,4.5,1,1)) # Change plot margins, make it look better
DU$plot_single_storm_event(event_example)
```

```
## Error in eval(expr, envir, enclos): object 'DU' not found
```

```r
# Clean up
rm(event_example)
```

```
## Warning in rm(event_example): object 'event_example' not found
```

**Here we plot the entire full_data time-series, with events overlain.** This
goes to a separate PDF file under FIG, because it is very large.


```r
# Plot will go in this directory
dir.create('FIG', showWarnings=FALSE)

# Function to plot all wave data, with storm events overlain
multi_year_pdf_plot<-function(site, event_data=NULL){

    # Open a pdf file to plot to
    pdf(paste0('FIG/yearly_timeseries_', site, '.pdf'), width=45, height=8)

    # Compute the axis limits for the plot
    max_hsig = max(wd[[site]]$hsig, na.rm=TRUE)
    max_tp1 = max(wd[[site]]$tp1, na.rm=TRUE)

    # Make a one page plot for each year
    for(year in 1985:2016){
        DU$wave_data_single_year_plot(year, site, wd, max_hsig, max_tp1, event_data,
            add_days=TRUE, add_event_start_lines=TRUE)
    }

    dev.off()
}

multi_year_pdf_plot('full_data', event_set$data)
```

```
## Error in multi_year_pdf_plot("full_data", event_set$data): object 'wd' not found
```

# **Step 2: Compute summary statistics for each storm event**
---------------------------------------------------------------

**The statistical analysis of event magnitude / frequency is performed on storm
event summary statistics -- and the latter are defined and extracted in the
following code**. Look at the computational routine in *data_utilities.R* for
more information on the extraction. At the time of writing we extract the
maximum significant wave height (m); TP1 (s) and direction (degrees from North)
at the time of peak significant wave height; the maximum tidal residual (m),
the event duration in hours, and the start time and end time as decimal years
(this is useful for some later analysis). The event summary statistics should
obviously reflect our event definition, which imposes a lower limit on the
significant wave height (`> hsig_threshold`), and the time between storms (`> 24 hours`).

```r
duration_offset_hours = 1.0 # Duration (hrs) for single-point event. Must be <= duration gap hours
stopifnot(duration_offset_hours <= duration_gap_hours)

event_statistics = DU$extract_event_statistics(event_set$data, median_tp1_dir=FALSE, 
    duration_offset_hours=duration_offset_hours)
```

```
## Error in eval(expr, envir, enclos): object 'DU' not found
```

```r
# Get a basic summary
summary(event_statistics, digits=6)
```

```
## Error in summary(event_statistics, digits = 6): object 'event_statistics' not found
```

```r
# Make pairwise scatterplots (ignoring time variables)
DU$nice_pairs(event_statistics[c('duration', 'hsig', 'tp1', 'dir', 'tideResid')])
```

```
## Error in eval(expr, envir, enclos): object 'DU' not found
```

# **Step 3: Study changes in monthly mean sea level, and remove seasonal and inter-annual trends from the tidal residual, to better estimate the storm surge**
-------------------------------------------------------------------------------------------------------------------------------------------------

**At this point, various exploratory analyses are undertaken to understand
non-stationarities and inhomogeneities in the data.** In particular,
we remove MSL related non-stationarities in the tidal residual, and also
perform an adjustment to the wave direction.



**Here we show that the tidal-residual is non-stationary**, with an obvious
increasing trend.

```r
# Tidal residual
scatter.smooth(event_statistics$startyear, event_statistics$tideResid, 
    xlab='Year', ylab='Tidal Residual (m)', col='blue', 
    cex = 0.5, main='Tidal residual over time')
```

```
## Error in xy.coords(x, y, xlabel, ylabel): object 'event_statistics' not found
```

```r
# Fit a linear regression
tidal_resid_vs_startyear = lm(event_statistics$tideResid ~ event_statistics$startyear)
```

```
## Error in eval(expr, envir, enclos): object 'event_statistics' not found
```

```r
abline(coef(tidal_resid_vs_startyear), col='red')
```

```
## Error in coef(tidal_resid_vs_startyear): object 'tidal_resid_vs_startyear' not found
```

```r
grid(col='brown')
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
```

```r
# Look at regression coefficients
summary(tidal_resid_vs_startyear)
```

```
## Error in summary(tidal_resid_vs_startyear): object 'tidal_resid_vs_startyear' not found
```

```r
# Clean up
rm(tidal_resid_vs_startyear)
```

```
## Warning in rm(tidal_resid_vs_startyear): object 'tidal_resid_vs_startyear'
## not found
```

```r
# A different check -- basic spearman correlation
cortest_startyear_tideResid_A = cor.test(event_statistics$startyear, 
    event_statistics$tideResid, method='s')
```

```
## Error in cor.test(event_statistics$startyear, event_statistics$tideResid, : object 'event_statistics' not found
```

```r
print(cortest_startyear_tideResid_A)
```

```
## Error in print(cortest_startyear_tideResid_A): object 'cortest_startyear_tideResid_A' not found
```

```r
# Can remove the warnings about p-values with ties by 'jittering' the data
# Qualitative result should be unchanged
cortest_startyear_tideResid_B = cor.test( 
    jitter(event_statistics$startyear), jitter(event_statistics$tideResid), 
    method='s')
```

```
## Error in jitter(event_statistics$startyear): object 'event_statistics' not found
```

```r
print(cortest_startyear_tideResid_B)
```

```
## Error in print(cortest_startyear_tideResid_B): object 'cortest_startyear_tideResid_B' not found
```

**The increasing trend in the surge might be reflective of changes in MSL**
(e.g.  due to climate change, among other forcings). To further investigate
this, we study monthly mean sea levels at Tomaree -- with the idea of removing
the long-term and seasonal trends in mean-sea-level from the surge residual we
have just defined.

```r
# Compute annual sea levels
yearly_SL = aggregate(wd$full_data$tide, list(floor(wd$full_data$year)), 
    f<-function(x) mean(x, na.rm=T))
```

```
## Error in aggregate(wd$full_data$tide, list(floor(wd$full_data$year)), : object 'wd' not found
```

```r
# Also compute shorter time-scale sea levels
monthly_SL = aggregate(wd$full_data$tide, 
    list(month=as.numeric(format(wd$full_data$time, '%m')), 
         year=as.numeric(format(wd$full_data$time, '%Y'))),
    f<-function(x) mean(x, na.rm=TRUE))
```

```
## Error in aggregate(wd$full_data$tide, list(month = as.numeric(format(wd$full_data$time, : object 'wd' not found
```

```r
monthly_tidal_residual = aggregate(wd$full_data$tideResid, 
    list(month=as.numeric(format(wd$full_data$time, '%m')), 
         year=as.numeric(format(wd$full_data$time, '%Y'))),
    f<-function(x) mean(x, na.rm=TRUE))
```

```
## Error in aggregate(wd$full_data$tideResid, list(month = as.numeric(format(wd$full_data$time, : object 'wd' not found
```

```r
# Convert the time of the monthly sea level /residual to decimal year, and
# gap-fill NA values
# Associate the monthly value to the 15th of the month
monthly_SL_year = DU$time_to_year(
    strptime(paste0(monthly_SL$year, '-', monthly_SL$month, '-15 00:00:00'), 
    format='%Y-%m-%d %H:%M:%S', tz='GMT'))
```

```
## Error in eval(expr, envir, enclos): object 'DU' not found
```

```r
monthly_SL_filled = monthly_SL$x
```

```
## Error in eval(expr, envir, enclos): object 'monthly_SL' not found
```

```r
monthly_tidal_residual_filled = monthly_tidal_residual$x
```

```
## Error in eval(expr, envir, enclos): object 'monthly_tidal_residual' not found
```

```r
nan_fix = which(is.nan(monthly_SL_filled)) 
```

```
## Error in which(is.nan(monthly_SL_filled)): object 'monthly_SL_filled' not found
```

```r
# The missing values of the residual and the monthly SL should be the same
# Confirm they are!
stopifnot(all(nan_fix == which(is.nan(monthly_tidal_residual$x))))
```

```
## Error in stopifnot(all(nan_fix == which(is.nan(monthly_tidal_residual$x)))): object 'nan_fix' not found
```

```r
# This will gap fill, assuming no consecutive nan values within the data, which
# is true for our current data (but check before applying to other data!)
monthly_SL_filled[nan_fix] = 0.5 * 
    (monthly_SL_filled[nan_fix - 1] + monthly_SL_filled[nan_fix + 1])
```

```
## Error in eval(expr, envir, enclos): object 'monthly_SL_filled' not found
```

```r
monthly_tidal_residual_filled[nan_fix] = 0.5 * 
    (monthly_tidal_residual_filled[nan_fix - 1] + 
     monthly_tidal_residual_filled[nan_fix + 1])
```

```
## Error in eval(expr, envir, enclos): object 'monthly_tidal_residual_filled' not found
```

```r
# Remove trailing nan
nan_fix = which(is.nan(monthly_SL_filled)) 
```

```
## Error in which(is.nan(monthly_SL_filled)): object 'monthly_SL_filled' not found
```

```r
monthly_SL_filled = monthly_SL_filled[-nan_fix]
```

```
## Error in eval(expr, envir, enclos): object 'monthly_SL_filled' not found
```

```r
monthly_SL_year = monthly_SL_year[-nan_fix]
```

```
## Error in eval(expr, envir, enclos): object 'monthly_SL_year' not found
```

```r
monthly_tidal_residual_filled = monthly_tidal_residual_filled[-nan_fix]
```

```
## Error in eval(expr, envir, enclos): object 'monthly_tidal_residual_filled' not found
```

```r
#
plot(monthly_SL_year, monthly_SL_filled, t='l')
```

```
## Error in plot(monthly_SL_year, monthly_SL_filled, t = "l"): object 'monthly_SL_year' not found
```

```r
points(monthly_SL_year, monthly_tidal_residual_filled, t='l', col='red')
```

```
## Error in points(monthly_SL_year, monthly_tidal_residual_filled, t = "l", : object 'monthly_SL_year' not found
```

```r
legend('topleft', c('Monthly sea level', 'Monthly tidal residual'), 
    col=c('black', 'red'), lty=c(1,1), pch=c(NA, NA))
```

```
## Error in strwidth(legend, units = "user", cex = cex, font = text.font): plot.new has not been called yet
```

The above figure suggests significant non-stationarities in the monthly mean
sea level.

**Here we apply the STL method to model non-stationarities in the monthly tidal
residual series** (Cleveland et al 1990). We decompose the sea level into an
seasonal periodic component and an annual trend. This is later used to correct
the derived non-astronimical tidal surge, so that the latter becomes more
obviously related to the storm related surge component.

```r
# Make a seasonal timeseries. 
monthly_tidal_residual_filled_ts = ts(monthly_tidal_residual_filled, 
    start=c(1985, 10), frequency=12)
```

```
## Error in is.data.frame(data): object 'monthly_tidal_residual_filled' not found
```

```r
# Stl smoothing
monthly_tidal_residual_filled_stl = stl(monthly_tidal_residual_filled_ts, 
    s.window='periodic')
```

```
## Error in as.ts(x): object 'monthly_tidal_residual_filled_ts' not found
```

```r
plot(monthly_tidal_residual_filled_stl)
```

```
## Error in plot(monthly_tidal_residual_filled_stl): object 'monthly_tidal_residual_filled_stl' not found
```

```r
## Convert to function
smooth_tideResid_fun_stl = approxfun(monthly_SL_year, 
    rowSums(monthly_tidal_residual_filled_stl$time.series[,1:2]),
    rule=2)
```

```
## Error in is.data.frame(x): object 'monthly_tidal_residual_filled_stl' not found
```

```r
# For convenience later, we also keep separate annual and monthly functions
smooth_tideResid_fun_stl_annual = approxfun(monthly_SL_year, 
    monthly_tidal_residual_filled_stl$time.series[,2],
    rule=2)
```

```
## Error in xy.coords(x, y): object 'monthly_tidal_residual_filled_stl' not found
```

```r
smooth_tideResid_fun_stl_monthly = approxfun(monthly_SL_year, 
    monthly_tidal_residual_filled_stl$time.series[,1],
    rule=2)
```

```
## Error in xy.coords(x, y): object 'monthly_tidal_residual_filled_stl' not found
```

Below we investigate relationships between MSL and annual mean SOI, as well as
long-term MSL changes. To enable further analyses we read climate variables,
computing various averages, and appending them to the event statistics to
support later analysis.

```r
# Get climate index info, along with a smoothed soi with df ~= number of years
CI = DU$read_climate_indices(
    soi_file = '../../Data/Climate_Indices/ENSO/SOI_BOM.txt',
    aao_file = '../../Data/Climate_Indices/AAO/AAO.txt')
```

```
## Error in eval(expr, envir, enclos): object 'DU' not found
```

```r
CI_annual = lapply(CI, 
    f<-function(x) aggregate(x[,2], list(year=as.numeric(format(x[,1], '%Y'))), mean ))
```

```
## Error in lapply(CI, f <- function(x) aggregate(x[, 2], list(year = as.numeric(format(x[, : object 'CI' not found
```

```r
CI_annual_fun = lapply(CI_annual, 
    f<-function(x) approxfun(x[,1], x[,2], method='constant') )
```

```
## Error in lapply(CI_annual, f <- function(x) approxfun(x[, 1], x[, 2], : object 'CI_annual' not found
```

```r
# In the table, soiA = annual mean soi, aaoA = annual mean aao, etc
for(nm in names(CI_annual_fun)){
    es_name = paste0(nm, 'A')
    event_statistics[[es_name]] = CI_annual_fun[[nm]](floor(event_statistics$startyear))
}
```

```
## Error in eval(expr, envir, enclos): object 'CI_annual_fun' not found
```

```r
soi_time = DU$time_to_year(CI$soi$time)
```

```
## Error in eval(expr, envir, enclos): object 'DU' not found
```

```r
soi_fun = approxfun(soi_time, CI$soi$index)
```

```
## Error in xy.coords(x, y): object 'CI' not found
```

```r
soi_NA = which(is.na(CI$soi$index))
```

```
## Error in which(is.na(CI$soi$index)): object 'CI' not found
```

```r
if(length(soi_NA) > 0){
    smooth_soi = smooth.spline(soi_time[-soi_NA], CI$soi$index[-soi_NA], 
        df = diff(range(soi_time[-soi_NA])))
}else{
    smooth_soi = smooth.spline(soi_time[], CI$soi$index[], 
        df = diff(range(soi_time[])))
}
```

```
## Error in eval(expr, envir, enclos): object 'soi_NA' not found
```

**Question: Is annual averaged SOI related to annual averaged sea level?**

```r
# Is annual average SOI related to annual average sea level?
# This is suggested in White et al (2014)
yearly_soi = CI_annual$soi
```

```
## Error in eval(expr, envir, enclos): object 'CI_annual' not found
```

```r
mm = match(yearly_SL[,1], yearly_soi[,1])
```

```
## Error in match(yearly_SL[, 1], yearly_soi[, 1]): object 'yearly_SL' not found
```

```r
cortest_SL_soi_A = cor.test(yearly_SL[,2], yearly_soi[mm,2]) 
```

```
## Error in cor.test(yearly_SL[, 2], yearly_soi[mm, 2]): object 'yearly_SL' not found
```

```r
print(cortest_SL_soiA) # Yes, positive relationship
```

```
## Error in print(cortest_SL_soiA): object 'cortest_SL_soiA' not found
```

**Question: How does the correlation between MSL and mean annual SOI change if
we assume recent sea level rise ~ 1.8mm/year (White et al., 2014)?**

```r
cortest_SL_soiB = cor.test(yearly_SL[,2] - 0.0018*(yearly_SL[,1] - 1985), 
    yearly_soi[mm,2])
```

```
## Error in cor.test(yearly_SL[, 2] - 0.0018 * (yearly_SL[, 1] - 1985), yearly_soi[mm, : object 'yearly_SL' not found
```

```r
print(cortest_SL_soiB) # Not much change compared to the result above
```

```
## Error in print(cortest_SL_soiB): object 'cortest_SL_soiB' not found
```

```r
# Here we make a simple linear regression of MSL, related to SOI and time
soi_SL_yearly = data.frame(soiA=yearly_soi[mm,2], sl=yearly_SL[,2], year = yearly_soi[mm,1])
```

```
## Error in data.frame(soiA = yearly_soi[mm, 2], sl = yearly_SL[, 2], year = yearly_soi[mm, : object 'yearly_soi' not found
```

```r
soi_SL_lm = lm(sl ~ soiA + year, data=soi_SL_yearly)
```

```
## Error in is.data.frame(data): object 'soi_SL_yearly' not found
```

```r
summary(soi_SL_lm) # Should suggest relations between MSL, SOI, and time
```

```
## Error in summary(soi_SL_lm): object 'soi_SL_lm' not found
```

Below we make some plots of the sea level information, and the tidal residual
after monthly SL residuals are removed. **The definition of the tidal residual
is changed in the following code** to remove the component related to
inter-annual and seasonal mean sea level. The idea is that A) the adjusted
tidal residual should more strongly reflect storm type processes. Further, B)
we can still model the influence of seasons and SOI/sea-level-rise on MSL, and
integrate that into our analysis, so nothing is lost by this modelling approach. 

```r
# Plot MSL over time, with SOI
par(mfrow=c(3,1))
plot(yearly_SL[,1]+0.5, yearly_SL[,2], xlab = 'Year', ylab='Annual mean sea level',
    main='Annual MSL over time with a smoother (red) and smooth scaled SOI (green)',
    cex.main=2)
```

```
## Error in plot(yearly_SL[, 1] + 0.5, yearly_SL[, 2], xlab = "Year", ylab = "Annual mean sea level", : object 'yearly_SL' not found
```

```r
points(monthly_SL_year, smooth_tideResid_fun_stl_annual(monthly_SL_year), t='l', col='red')
```

```
## Error in points(monthly_SL_year, smooth_tideResid_fun_stl_annual(monthly_SL_year), : object 'monthly_SL_year' not found
```

```r
points(smooth_soi$x, smooth_soi$y/200, t='l', col='green')
```

```
## Error in points(smooth_soi$x, smooth_soi$y/200, t = "l", col = "green"): object 'smooth_soi' not found
```

```r
points(yearly_soi[,1]+0.5, yearly_soi[,2]/200, col='green')
```

```
## Error in points(yearly_soi[, 1] + 0.5, yearly_soi[, 2]/200, col = "green"): object 'yearly_soi' not found
```

```r
abline(v=1985:2015, col='grey')
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
```

```r
grid(col='brown')
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
```

```r
plot(monthly_SL_year, monthly_SL_filled, t='l', 
    xlab='Year', ylab='Monthly mean sea level',
    main='Monthly sea level with smoother (red), and scaled SOI (green)',
    cex.main=2)
```

```
## Error in plot(monthly_SL_year, monthly_SL_filled, t = "l", xlab = "Year", : object 'monthly_SL_year' not found
```

```r
points(monthly_SL_year, smooth_tideResid_fun_stl(monthly_SL_year), 
    col='red', t='l')
```

```
## Error in points(monthly_SL_year, smooth_tideResid_fun_stl(monthly_SL_year), : object 'monthly_SL_year' not found
```

```r
points(DU$time_to_year(CI$soi$time), CI$soi$index/200, t='l', 
    col='green')
```

```
## Error in points(DU$time_to_year(CI$soi$time), CI$soi$index/200, t = "l", : object 'DU' not found
```

```r
grid(col='brown')
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
```

```r
abline(v=1985:2016, col='grey')
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
```

```r
# Use the STL decomposition to adjust the tidal residuals
SL_adjustment_events = smooth_tideResid_fun_stl(event_statistics$startyear)
```

```
## Error in eval(expr, envir, enclos): could not find function "smooth_tideResid_fun_stl"
```

```r
#
# REMOVE THE MSL COMPONENT FROM THE EVENT STATISTICS
#
event_statistics$tideResid = event_statistics$tideResid - SL_adjustment_events
```

```
## Error in eval(expr, envir, enclos): object 'event_statistics' not found
```

```r
# Append msl to the event_statistics table
event_statistics$msl = SL_adjustment_events
```

```
## Error in eval(expr, envir, enclos): object 'SL_adjustment_events' not found
```

```r
# Add in a plot of the NEW tidal residual
plot(event_statistics$startyear, event_statistics$tideResid, col='blue',
    xlab='Time', ylab='Adjusted Tidal Residual', 
    main='Tidal Residual, adjusted for annual and monthly drifts in the sea level',
    cex.main=2)
```

```
## Error in plot(event_statistics$startyear, event_statistics$tideResid, : object 'event_statistics' not found
```

```r
grid(col='brown'); abline(h=0, col='red')
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
```

**Here we perform a few statistical checks on the adjusted tidal residual, to see
whether long-term and ENSO related trends have been removed**

```r
# Is there still a small increasing trend in the tidal residual? It's possible,
# but the result is statistically borderline, and is much weaker than the result
# prior to the removal of the monthly sea level trends (reported above)
cortest_startyear_tr_new = cor.test(event_statistics$startyear, event_statistics$tideResid, method='s')
```

```
## Error in cor.test(event_statistics$startyear, event_statistics$tideResid, : object 'event_statistics' not found
```

```r
print(cortest_startyear_tr_new)
```

```
## Error in print(cortest_startyear_tr_new): object 'cortest_startyear_tr_new' not found
```

```r
# Here check the same thing as above, with linear regression. 
lm_startyear_tr_new = lm(event_statistics$tideResid ~ event_statistics$startyear)
```

```
## Error in eval(expr, envir, enclos): object 'event_statistics' not found
```

```r
# Suggests a small increasing trend, with borderline statistical significance
summary(lm_startyear_tr_new)
```

```
## Error in summary(lm_startyear_tr_new): object 'lm_startyear_tr_new' not found
```

```r
# What about a relation between the surge and soiA ?
cortest_soiA_tr_new = cor.test(event_statistics$soiA, event_statistics$tideResid, method='s')
```

```
## Error in cor.test(event_statistics$soiA, event_statistics$tideResid, method = "s"): object 'event_statistics' not found
```

```r
print(cortest_soiA_tr_new) # Seems to have been convincingly removed.
```

```
## Error in print(cortest_soiA_tr_new): object 'cortest_soiA_tr_new' not found
```

```r
# Any warnings about 'cannot compute exact p-values with ties' can be cross-checked
# by removing ties from the data with a jitter. e.g.
cortest_soiA_tr_new2 = cor.test(jitter(event_statistics$soiA), 
    jitter(event_statistics$tideResid), method='s') 
```

```
## Error in jitter(event_statistics$soiA): object 'event_statistics' not found
```

```r
print(cortest_soiA_tr_new2) # Should be very similar to the last one
```

```
## Error in print(cortest_soiA_tr_new2): object 'cortest_soiA_tr_new2' not found
```

```r
# As an alternative to correlation, let's check by fitting a linear regression. 
lm_tr_soiA_new = summary(lm(event_statistics$tideResid ~ event_statistics$soiA))
```

```
## Error in eval(expr, envir, enclos): object 'event_statistics' not found
```

```r
# Should not suggest a significant trend [since we removed SOI from the tidal
# residual already]
print(summary(lm_tr_soiA_new)) 
```

```
## Error in summary(lm_tr_soiA_new): object 'lm_tr_soiA_new' not found
```























