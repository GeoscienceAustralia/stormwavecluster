
# **Creation of storm wave event summary statistics**
------------------------------------------------------

*Gareth Davies, Geoscience Australia 2017*

# Introduction
------------------

This document illustrates the extraction of storm wave event summary statistics
from observational wave and tide data, using R code. It was generated from the
file 'preprocess_data.Rmd', using the literate programming R package *knitr*. 

The basic approach followed here is to:
* **Step 1:** Parse relevant time-series data (from models or measurements) describing storm waves, sea levels, and climatic indices of interest. 
* **Step 2:** Convert the observational data to a single time-series representing observations at Crowdy Head. 
* **Step 3:** Extract storm-events from these time-series, based on some precise definition of storms.
* **Step 4:** Extract *summary-statistics* describing each storm event.

The above summary statistics will later be subject to statistical modelling.


# **Step 1: Parse the time-series data**
----------------------------------------

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


