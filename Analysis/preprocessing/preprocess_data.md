
# **Creation of storm wave event summary statistics**
------------------------------------------------------

*Gareth Davies, Geoscience Australia 2017*

# Introduction
------------------

This document illustrates the extraction of storm wave event summary statistics
from observational wave and tide data, using R code. It was generated from the
file [preprocess_data.Rmd](preprocess_data.Rmd), using the literate programming R package *knitr*. 

If you have R installed, along with all the packages required to run this code,
and a copy of the *stormwavecluster* git repository, then you should be able to
re-run the analysis here by simply copy-pasting the code (or using *knitr* to
extract it to a script).

The basic approach followed here is to:
* **Step 1:** Parse relevant wave time-series data at a number of sites (all near-ish to Old Bar), and convert them to a single time-series representing waves at Old Bar. 
* **Step 2:** Parse tidal observations, and astronomical tidal predictions, for a site near Old Bar.
* **Step 3:** Extract storm-events from these time-series, based on some precise definition of storms.
* **Step 4:** Extract *summary-statistics* describing each storm event.

The above summary statistics will later be subject to statistical modelling.


# **Step 1: Parse the wave time-series data**
--------------------------------------------

# *Parsing the measured wave time-series at Crowdy Head, Coffs Harbour and Sydney*
------------------------------------------------------------------------------------------

Here we parse measured wave time-series from Crowdy Head (since that is near
Old Bar, our primary site of interest). We also parse similar data measured at
Coffs Harbour and Sydney, for the purposes of filling gaps in the Crowdy Head
observations. 

First we read a script 'data_utilities.R' which contains functions for things
like 'reading the MHL wave data' into an environment. By using an environment,
we avoid the possibility of accidentally over-writing variables inside
'data_utilities.R'. You will have to read the 'data_utilities.R' code and comments
to understand each function. 


```r
DU = new.env()
source('data_utilities.R', local=DU) 
```

Wave data for Crowdy Head, Coffs Harbour and Sydney was kindly provided by
Manly Hydraulics Laboratory, and is included in this repository. First we get
data measured from 1985-2014.


```r
mhl_wave_dir = '../../Data/NSW_Waves/'

# There are zipped csv files for each site, for the data from 1985-2014
mhl_wave_files = Sys.glob(paste0(mhl_wave_dir, '*1985*.csv.zip'))
mhl_wave_files
```

```
## [1] "../../Data/NSW_Waves/COFHOW 10-10-1985 to 1-11-2014.csv.zip"
## [2] "../../Data/NSW_Waves/CRHDOW 10-10-1985 to 1-11-2014.csv.zip"
## [3] "../../Data/NSW_Waves/SYDDOW 10-10-1985 to 1-11-2014.csv.zip"
```

```r
# Make a short name for each site
mhl_wave_sites = substr(basename(mhl_wave_files), 1, 4)
mhl_wave_sites
```

```
## [1] "COFH" "CRHD" "SYDD"
```

```r
# Read the data
# wd = "wave data"
wd = DU$parse_MHL_wave_buoy_data(mhl_wave_sites, mhl_wave_files)
```

Next, we append data from 2015 to the above 1985-2014 data. We also read data
from another Sydney station which includes hindcast wave directions (based on
meteorological charts).

**Append the more recent wave data to `wd`**

```r
# Update: In early 2016 we received updated wave-buoy data
#
# Append the updated MHL wave buoy data
mhl_wave_files = Sys.glob(paste0(mhl_wave_dir, '*2016*.csv.zip'))
mhl_wave_sites = substr(basename(mhl_wave_files), 1, 4)
wd_update = DU$parse_MHL_wave_buoy_data(mhl_wave_sites, mhl_wave_files)

for(nm in names(wd)){
    matchInds = match(wd_update[[nm]]$time, wd[[nm]]$time)

    ## Lots of checks here

    # From graphical checks the overlapping data seems identical
    # So matchInds should be a sequence of consecutive integers, followed
    # by a sequence of NA's
    stopifnot( (max(which(!is.na(matchInds))) + 1) == min(which(is.na(matchInds))) )
    stopifnot( max(diff(matchInds), na.rm=TRUE)  == 1)
    stopifnot( min(diff(matchInds), na.rm=TRUE)  == 1)

    # Ensure names and overlapping data are identical
    stopifnot(all(names(wd_update[[nm]]) == names(wd[[nm]])))
    for(varname in names(wd_update[[nm]])){
        stopifnot(
            all(range(wd_update[[nm]][[varname]] - wd[[nm]][[varname]][matchInds], 
                    na.rm=TRUE) == 
                c(0, 0))
            )
    }

    # If all those tests have passed, we can update the data with an rbind 
    ll = matchInds[1] - 1
    wd[[nm]] = rbind(wd[[nm]][1:ll,], wd_update[[nm]])
}


# Get the other sydney data. It's in a different format, so we parse it here,
# and append it to the 'wd' list

syd1 = read.table(
    unz(description = paste0(mhl_wave_dir, 'SYDNOW 17-7-1987 to 4-10-2000.csv.zip'), 
        filename = 'SYDNOW 17-7-1987 to 4-10-2000.csv'),
    skip=7, header=TRUE, sep=",")

# Add a time variable
syd1$time = strptime(syd1$Date.Time, format='%d-%B-%Y %H:%M', tz='Etc/GMT-10')
# Add names to data.frame
names(syd1) = c('datetime', 'hsig', 'hmax', 'tz', 'tsig', 'tp1', 'dir', 'time')
# Add a year variable
syd1$year = DU$time_to_year(syd1$time)
# Append to 'wd' list under the name 'SYDL'
wd$SYDL = syd1
```

As a result of the above operations, we have a list `wd` containing data.frames
for each station. The station names and data dimensions are:

```r
# Station names
names(wd)
```

```
## [1] "COFH" "CRHD" "SYDD" "SYDL"
```

```r
# Number of rows/columns for each station
lapply(wd, dim)
```

```
## $COFH
## [1] 229880      7
## 
## $CRHD
## [1] 229983      7
## 
## $SYDD
## [1] 174953      7
## 
## $SYDL
## [1] 107378      9
```
The data format looks like this (using the example of Crowdy Head):

```r
wd$CRHD[1:10,]
```

```
##                   time  hsig  hmax  tz tp1 dir     year
## 1  1985-10-10 08:00:00 0.953 1.723 5.8 9.5  NA 1985.774
## 2  1985-10-10 09:00:00 0.827 1.422 5.8 9.5  NA 1985.774
## 3  1985-10-10 10:00:00 0.913 1.657 6.5 9.5  NA 1985.774
## 4  1985-10-10 11:00:00 0.895 1.513 6.2 9.5  NA 1985.774
## 5  1985-10-10 12:00:00 0.911 1.561 6.3 8.2  NA 1985.774
## 6  1985-10-10 13:00:00 1.051 1.764 6.7 8.2  NA 1985.774
## 7  1985-10-10 14:00:00 0.863 1.493 6.0 8.8  NA 1985.774
## 8  1985-10-10 15:00:00 0.951 1.852 6.5 8.8  NA 1985.774
## 9  1985-10-10 16:00:00 0.926 1.459 5.8 9.5  NA 1985.774
## 10 1985-10-10 17:00:00 0.968 1.816 5.7 8.8  NA 1985.775
```
The data for the SYDL station actually has a few more columns. That does not
affect the analysis so long as it also contains columns with the same names and
data as observed at the other stations.

```r
wd$SYDL[1:10,]
```

```
##             datetime  hsig  hmax  tz tsig  tp1 dir                time
## 1  17-JUL-1987 13:00 1.001 2.185 5.1  9.0 11.1 123 1987-07-17 13:00:00
## 2  17-JUL-1987 14:00 1.063 1.982 5.9 10.0 11.1 123 1987-07-17 14:00:00
## 3  17-JUL-1987 15:00 0.994 1.762 6.3 10.4 11.1 122 1987-07-17 15:00:00
## 4  17-JUL-1987 16:00 0.912 1.736 5.8  9.7 11.1 122 1987-07-17 16:00:00
## 5  17-JUL-1987 17:00 0.915 1.721 5.6  9.8 11.1 121 1987-07-17 17:00:00
## 6  17-JUL-1987 18:00 0.910 1.767 5.8  9.9 10.2 121 1987-07-17 18:00:00
## 7  17-JUL-1987 19:00 0.844 1.701 5.7 10.0 10.2 120 1987-07-17 19:00:00
## 8  17-JUL-1987 20:00 0.869 1.710 5.9 10.0 11.1 120 1987-07-17 20:00:00
## 9  17-JUL-1987 21:00 0.945 1.703 5.7 10.0 11.1 119 1987-07-17 21:00:00
## 10 17-JUL-1987 22:00 1.031 1.970 5.8  9.9 11.1 119 1987-07-17 22:00:00
##        year
## 1  1987.541
## 2  1987.541
## 3  1987.541
## 4  1987.542
## 5  1987.542
## 6  1987.542
## 7  1987.542
## 8  1987.542
## 9  1987.542
## 10 1987.542
```

The data_utilities.R script includes a function to plot the station data.
Here's an example. Note there are gaps in the data, which will be filled at
a later stage of the analysis.

```r
# Plot 2013 for Crowdy head
DU$wave_data_single_year_plot(year=2013, site='CRHD', wd=wd, max_hsig=8, max_tp1=15)
```

![plot of chunk plot2013](figure/plot2013-1.png)

# *Converting the wave observations to a single 'Old Bar' time-series*
----------------------------------------------------------------------

Here we create a 'gap-filled' wave time-series which we will treat as
representative of wave conditions at Old Bar. 

Because Crowdy Head is much closer to Old Bar than the other wave measuring
sites, it is natural to take measurements at Crowdy Head as a 'first
preference' representation of waves at Old Bar. When Crowdy Head data was
missing, we decided to gap-fill preferentially with data from Coffs Harbour. If
the latter were missing, we used Sydney observations, firstly from the `SYDD` site
(the directional wave-rider buoy, measured since 1992), and secondly from the
`SYDL` site (the long-reef wave-rider measurements complemented with hind-cast
wave-directional data). 

These preferences were justified by comparison of observed wave height and
direction at the gap-filling stations with those from Crowdy Head during
high-wave events. A few plots comparing wave heights and directions are shown
below. 

**Wave directions during storms @ Crowdy Head, in 2013**

```r
# Compare wave directions in a year, when waves at Crowdy Head exceed hsig_thresh
hsig_thresh = 2.9
year2compare = 2013

par(mfrow=c(1,2))
DU$check_station_correlations(year2compare, 'CRHD', 'SYDD', wd, 'dir', 
    site1_restriction = (wd$CRHD$hsig > hsig_thresh))
```

```
## [1] 0.8254451
```

```r
title(main='Wave direction at Crowdy Head vs Sydney', line=0.5)
DU$check_station_correlations(year2compare, 'CRHD', 'COFH', wd, 'dir', 
    site1_restriction = (wd$CRHD$hsig > hsig_thresh))
```

```
## [1] 0.899565
```

```r
title(main='Wave direction at Crowdy Head vs Coffs Harbour', line=0.5)
```

![plot of chunk direction_compare](figure/direction_compare-1.png)

```r
# Note -- we do not have SYDL direction data from this year
```

**Significant wave height during storms @ Crowdy Head, in 2013**

```r
# Compare 2013
par(mfrow=c(2,2))
DU$check_station_correlations(year2compare, 'CRHD', 'SYDD', wd, 'hsig',
    site1_restriction = (wd$CRHD$hsig > hsig_thresh))
```

```
## [1] 0.3872359
```

```r
title(main=bquote(H[sig] ~ 'at Crowdy Head vs Sydney'), line=0.5)
DU$check_station_correlations(year2compare, 'CRHD', 'COFH', wd, 'hsig',
    site1_restriction = (wd$CRHD$hsig > hsig_thresh))
```

```
## [1] 0.6368558
```

```r
title(main=bquote(H[sig] ~ 'at Crowdy Head vs Coffs Harbour'), line=0.5)
DU$check_station_correlations(1990, 'CRHD', 'SYDL', wd, 'hsig',
    site1_restriction = (wd$CRHD$hsig > hsig_thresh))
```

```
## [1] 0.3298255
```

```r
title(main=bquote(H[sig] ~ 'at Crowdy Head vs Long Reef'), line=0.5)
```

![plot of chunk hsig_compare](figure/hsig_compare-1.png)


**Make the gap-filled data, stored in a variable `full_data`**

```r
# Get times to interpolate at
len_crhd = length(wd$CRHD$time)
desired_times = seq(wd$CRHD$time[1], wd$CRHD$time[len_crhd], by='hour')

# Get the interpolated 'full' data
site_preference_order = c('CRHD', 'COFH', 'SYDD', 'SYDL')
full_data = DU$gap_fill_wave_data(desired_times, site_preference_order, wd)
head(full_data)
```

```
##                  time  hsig  hmax  tz tp1 dir     year waves_site dir_site
## 1 1985-10-10 08:00:00 0.953 1.723 5.8 9.5  NA 1985.774       CRHD     SYDL
## 2 1985-10-10 09:00:00 0.827 1.422 5.8 9.5  NA 1985.774       CRHD     SYDL
## 3 1985-10-10 10:00:00 0.913 1.657 6.5 9.5  NA 1985.774       CRHD     SYDL
## 4 1985-10-10 11:00:00 0.895 1.513 6.2 9.5  NA 1985.774       CRHD     SYDL
## 5 1985-10-10 12:00:00 0.911 1.561 6.3 8.2  NA 1985.774       CRHD     SYDL
## 6 1985-10-10 13:00:00 1.051 1.764 6.7 8.2  NA 1985.774       CRHD     SYDL
```

```r
tail(full_data)
```

```
##                       time  hsig hmax   tz  tp1 dir     year waves_site
## 265692 2016-01-31 19:00:00 0.970 2.07 5.09 9.77  98 2016.084       CRHD
## 265693 2016-01-31 20:00:00 0.899 1.72 5.02 8.17  84 2016.084       CRHD
## 265694 2016-01-31 21:00:00 0.936 1.73 5.24 8.52  89 2016.084       CRHD
## 265695 2016-01-31 22:00:00 1.008 1.75 4.99 9.32  95 2016.084       CRHD
## 265696 2016-01-31 23:00:00 0.981 1.76 5.14 9.77  94 2016.085       CRHD
## 265697 2016-02-01 00:00:00 0.946 1.77 5.21 8.17  71 2016.085       CRHD
##        dir_site
## 265692     CRHD
## 265693     CRHD
## 265694     CRHD
## 265695     CRHD
## 265696     CRHD
## 265697     CRHD
```

```r
# Append the 'full_data' to wd, and plot it
wd$full_data = full_data
```

The following plots show that most wave directions in `full_data` originate
from Crowdy Head, whereas most wave directions originate from the Sydney
waverider buoy. This is inevitable, because wave direction was only measured at
Crowdy Head and Coffs Harbour after ~ 2011, while measurements have been taken at
Sydney since 1992.


```r
par(mfrow=c(1,2))
pie(table(full_data$waves_site), main='Source of wave data in Old Bar wave series')
pie(table(full_data$dir_site), main='Source of direction data in Old Bar wave series')
```

![plot of chunk gap_filling_check](figure/gap_filling_check-1.png)

