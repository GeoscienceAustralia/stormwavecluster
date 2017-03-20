# Various utility functions for the storm clustering analysis
#
# AUTHORS
# Gareth Davies, Geoscience Australia 2014-15, gareth.davies.ga.code@gmail.com
# (Contributors add name here)

###############################################################################
#'
#' Return the current value of .Random.seed 
#'
#' Utility function to get the random seed -- and make it first if required
#'
#' @export
#' @examples
#'  x = get_random_seed()
#'
get_random_seed<-function(){

    if(!exists('.Random.seed', where=.GlobalEnv)){
        # Create the random seed by calling a random number generator
        x = rnorm(1)
    }

    return(.Random.seed)
}


#' Determine if a year is a leap year
is_leap_year<-function(year){

    (year%%4 == 0)*( (year%%100 != 0) + (year%%100 == 0)*(year%%400 == 0))
}

.test_is_leap_year<-function(){

    years =      c(1800, 1804, 1900, 1952, 2000, 2004, 2005, 2200)

    leap_years = c(F,    T   , F   , T   , T   , T   , F   , F   )

    if(any(is_leap_year(years)!=leap_years)){
        print('FAIL')
    }else{
        print('PASS')
    }

    days_in_many_yr = strptime('9999-01-01', format='%Y-%m-%d', tz='GMT') - 
                      strptime('00-01-01', format='%Y-%m-%d', tz='GMT')

    days_in_many_yr_B = sum(is_leap_year(0:9998))*366 + 
                        sum(!is_leap_year(0:9998))*365

    if(days_in_many_yr == days_in_many_yr_B){
        print('PASS')
    }else{
        print('FAIL')
    }

}

#' Convert a decimal year to a strptime object. 
#' @param year time as a decimal year
#' @param tz timezone as recognized by strptime
year_to_time<-function(year, tz='Etc/GMT-10'){

    ## Max 'year' allowed by strptime is 9999
    #if(year >= 1e+04){
    #    extra_ten_thousands = floor(year/1e+04)
    #    year = year - 1e+04*extra_ten_thousands

    #    extra_days = as.difftime(366, units='days')*2500 + 
    #                 as.difftime(365, units='days')*7500
    #}else{
    #    extra_days = as.difftime(0, units='days')
    #}

    extra_ten_thousands = floor(year/1e+04)

    year = year - 1e+04*extra_ten_thousands

    tenK_days = sum(is_leap_year(0:9999))*366 + sum(!is_leap_year(0:9999))*365

    extra_days = as.difftime(extra_ten_thousands*tenK_days, units='days')

    days_in_year = 365 + 1*is_leap_year(floor(year))

    day_decimal = (year - floor(year))*days_in_year

    year_strp = strptime(paste0(floor(year), '-01-01'), format='%Y-%m-%d', 
        tz=tz)

    day_difftime = as.difftime(day_decimal, units='days')

    output_date = year_strp + day_difftime + extra_days

    return(output_date)
}

.test_year_to_time<-function(){
    
    # Test 1
    year = 2400
    year_strp = year_to_time(year)

    if(as.character(year_strp)=='2400-01-01'){
        print('PASS')
    }else{
        print('FAIL')
    }    

    # Test 2
    year = 2400.5
    year_strp = year_to_time(year)

    if(as.character(year_strp)=='2400-07-02'){
        print('PASS')
    }else{
        print('FAIL')
    }    

    # Test 3
   
    year = 2000.184734 
    year_strp = year_to_time(year)
    year_back = time_to_year(year_strp)

    if(abs(year - year_back) < 1/(365*24*60*60)){
        # Error < 1s
        print('PASS')
    }else{
        print('FAIL')
    }

}

###############################################################################
#'
#' Convert times to a decimal year IN THEIR OWN TIME ZONE (going down to seconds)
#'
time_to_year<-function(strptime_vec){

    strptime_year = format(strptime_vec, '%Y')

    # Figure out if it is a leap year
    days_in_year = 365 + 1*is_leap_year(as.numeric(strptime_year))

    day_of_year = format(strptime_vec, '%j')
    hour_of_day = format(strptime_vec, '%H')
    minute_of_hour = format(strptime_vec, '%M')
    second_of_minute = format(strptime_vec, '%S')

    output_time = as.numeric(strptime_year) +
                  (as.numeric(day_of_year)-1)/days_in_year +
                  as.numeric(hour_of_day)/(24*days_in_year) +
                  as.numeric(minute_of_hour)/(60*24*days_in_year) +
                  as.numeric(second_of_minute)/(60*60*24*days_in_year)
                    

    return(output_time)
}

#'
.test_time_to_year<-function(){

    Y = strptime('1985-01-01 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='GMT')
    Y_frac = time_to_year(Y)
    stopifnot(Y_frac==1985)
    print('PASS')

    # Check non trivial timezones (but we report in our own timezone)
    Y = strptime('1985-01-01 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT-10')
    Y_frac = time_to_year(Y)
    stopifnot(Y_frac==(1985))
    print('PASS')

    Y = strptime('1985-01-01 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT+10')
    Y_frac = time_to_year(Y)
    stopifnot(Y_frac==(1985))
    print('PASS')


    Y = strptime('1985-02-01 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='GMT')
    Y_frac = time_to_year(Y)
    stopifnot(all.equal(Y_frac,1985 + (31/365)))
    print('PASS')

    # More complex
    Y = strptime('1985-03-02 02:30:05', format='%Y-%m-%d %H:%M:%S', tz='GMT')
    Y_frac = time_to_year(Y)
    expected_time = 1985 + (31+28+1)/365 + 2/(24*365) + 30/(60*24*365) + 5/(60*60*24*365) 
    stopifnot(all.equal(Y_frac,expected_time))
    print('PASS')
 
    # Test with leap year   
    Y = strptime('1988-03-01 02:30:05', format='%Y-%m-%d %H:%M:%S', tz='GMT')
    Y_frac = time_to_year(Y)
    stopifnot(all.equal(Y_frac,1988 + ((31+29)/366+ (2*60*60+30*60+5)/(366*24*60*60))))
    print('PASS')

    # Test with leap year
    Y = strptime('1988-03-01 23:30:05', format='%Y-%m-%d %H:%M:%S', tz='GMT')
    Y_frac = time_to_year(Y)
    stopifnot(all.equal(Y_frac,1988 + ((31+29)/366+ (23*60*60+30*60+5)/(366*24*60*60))))
    print('PASS')

    # Weird date with time of 24:00:00
    Y = strptime('1988-12-31 24:00:00', format='%Y-%m-%d %H:%M:%S', tz='GMT')
    Y_frac = time_to_year(Y)
    stopifnot(all.equal(Y_frac,1989))
    print('PASS')

}

#'
#' Function to read MHL wave buoy data
#' @param MHL_wave_sites A character vector giving a name for each site
#' @param MHL_wave_files A character vector giving the filenames corresponding to each site name
#' @return A named list containing one data.frame for each site. The data.frame
#'   contains time, hsig, hmax, tz, tp1, dir. See the MHL website for explanation.
parse_MHL_wave_buoy_data<-function(
    MHL_wave_sites=mhl_wave_sites, 
    MHL_wave_files=mhl_wave_files){

    # Store all sites in lists
    wd = list() # wave data
    for(i in 1:length(MHL_wave_sites)){
        site = MHL_wave_sites[i]
        # Read csv
        if(grepl('\\.zip', MHL_wave_files[i])){
            wave_data = read.csv(
                unz(description=MHL_wave_files[i], 
                    filename=basename(gsub('\\.zip', '', MHL_wave_files[i]))), 
                skip=9, header=TRUE)
        }else{
            wave_data = read.csv(MHL_wave_files[i], skip=9, header=TRUE)
        }
        # Correct erronious final column, caused by ',' in header
        wave_data = wave_data[,1:6]
        # Adjust column names
        names(wave_data) = c('time', 'hsig', 'hmax', 'tz', 'tp1', 'dir')

        # Convert time to date-time object -- according to MHL website, waverider
        # time-zone is UTC+10 (beware R denotes this with a minus, as 'Etc/GMT-10')
        # e.g. http://www.mhl.nsw.gov.au/htbin/wave_data_plot.com?Location=Sydney
        time_data = strptime(wave_data[,1], 
                             format = '%d-%b-%Y %H:%M', 
                             tz = 'Etc/GMT-10')

        # Also store time as a decimal year -- useful
        year_data = DU$time_to_year(time_data)

        # Append the time as a decimal year to wd
        wd[[site]] = cbind(wave_data, year=year_data)
        # Replace the character time with the date-time object
        wd[[site]]$time = time_data
    }

    # Put these variables in the global environment
    return(wd)
}


#'
#' Parse the MHL tidal gauge data
#' 
#' @param read_MHL_csv_tide_gauge
#' @return A data.frame containing the time, julian time (days since 1970), tidal level, and status
read_MHL_csv_tide_gauge<-function(gauge_file){

    if(grep('\\.zip', gauge_file)){
        tomaree = read.csv(
            unz(gauge_file, filename=basename(gsub('\\.zip', '', gauge_file))), 
            skip=17,header=T,
            na.strings=c('NA', 'n/a'), stringsAsFactors=FALSE)
    }else{
        tomaree = read.csv(gauge_file, skip=17,header=T,
            na.strings=c('NA', 'n/a'), stringsAsFactors=FALSE)
    }

    JULIAN_TIME_ORIGIN = strptime("1970-01-01 00:00:00", format='%Y-%m-%d %H:%M:%S', tz = "Etc/GMT-10")
    tomaree_time = strptime(paste(tomaree[,1], tomaree[,2]), 
        format='%d/%m/%Y %H:%M:%S', tz='Etc/GMT-10')
    tomaree_julian_time = julian(tomaree_time, 
        origin = JULIAN_TIME_ORIGIN)

    output = data.frame(time=tomaree_time, 
                        julian_time = as.numeric(tomaree_julian_time),
                        tide=tomaree[,3], status=tomaree[,4])
    return(output)
}


#'
#' Plot a single year of the  MHL wave buoy data, for a particular site
#'
#' @param year Integer year to plot
#' @param site site name which appears in wd
#' @param wd output of parse_MHL_wave_buoy_data
#' @param max_hsig Upper limit on hsig for the plot
#' @param max_tp1 Upper limit on tp1 for the plot
#' @param event_data if not null, a list of data.frames each containing events with the same variables as wd[[site]]
#' @param add_days Add day tickmarks to the plot
#' @param add_event_start_lines Add lines showing event starts
#' @return Doesn't return anything, but makes a plot
wave_data_single_year_plot<-function(year, site, wd, max_hsig, max_tp1, 
    event_data=NULL, add_days=FALSE, add_event_start_lines=FALSE){

    # Get data subset
    data_year = format(wd[[site]]$time, '%Y')
    keep = which(data_year == year)
    if(length(keep) == 0){
        return()
    }
     
    # Get the data-subset (convenience)
    ds = wd[[site]][keep,]
    
    if(add_days){
        days = unique(round(ds$time, 'day'))
        days_tick = as.numeric(days)

        days_lab = days[which(format(days, '%d') %in% c('07', '14', '21', '28' ))]

    }


    if('tide'%in%names(ds)){
        # Plot tides as well
        par(mfrow=c(4,1))
        par(mar=c(3,4,1,1)) # Adjust margins 
    }else{
        par(mfrow=c(3,1)) # 3 plots per page, arranged vertically
        par(mar=c(3,3,2,1)) # Adjust margins 
    }

    # Hsig plot
    plot(ds$time, ds$hsig,t='h', 
         col='orange', lend=1,
         #ylab='Hsig', main=paste0('Hsig: ', year,' ', site),
         ylab='Hsig (m)', main=paste0('Significant Wave Height: ', year),
         ylim=c(0, max_hsig))
    points(ds$time, ds$hsig,t='l')
    grid()

    if(!is.null(event_data)){
        for(i in 1:length(event_data)){
            points(event_data[[i]]$time, event_data[[i]]$hsig, t='l', col='red')
            if(add_event_start_lines){
                abline(v = as.numeric(event_data[[i]]$time[1]), col='pink', lty='solid')
                #text(as.numeric(event_data[[i]]$time[1]), max_hsig, label=as.character(i), cex=0.5)
            }
        }
    }
    if(add_days){
        axis(side=1, at=days_tick, labels=FALSE, tick=TRUE)
        axis(side=1, at=as.numeric(days_lab), labels=format(days_lab,'%d'))
    }

    # TP1 plot
    plot(ds$time, ds$tp1,t='h', 
         col='blue', lend=1, 
         #ylab='TP1', main=paste0('T P1: ', year,' ', site), 
         ylab='Tp1 (s)', main=paste0('Wave period: ', year),
         ylim=c(0, max_tp1))
    points(ds$time, ds$tp1,t='l')
    grid()
    if(!is.null(event_data)){
        for(i in 1:length(event_data)){
            points(event_data[[i]]$time, event_data[[i]]$tp1, t='l', col='red')
            if(add_event_start_lines) abline(v = as.numeric(event_data[[i]]$time[1]), col='pink', lty='solid')
        }
    }
    if(add_days){
        axis(side=1, at=days_tick, labels=FALSE, tick=TRUE)
        axis(side=1, at=as.numeric(days_lab), labels=format(days_lab,'%d'))
    }

    # Wave direction plot
    plot(ds$time, ds$dir,t='h', 
         col='green', lend=1, 
         ylab='Direction (degrees)', main=paste0('Wave Direction: ', year), 
         ylim=c(0, 360))
    points(ds$time, ds$dir,t='l')
    grid()
    if(!is.null(event_data)){
        for(i in 1:length(event_data)){
            points(event_data[[i]]$time, event_data[[i]]$dir, t='l', col='red')
            if(add_event_start_lines) abline(v = as.numeric(event_data[[i]]$time[1]), col='pink', lty='solid')
        }

    }
    if(add_days){
        axis(side=1, at=days_tick, labels=FALSE, tick=TRUE)
        axis(side=1, at=as.numeric(days_lab), labels=format(days_lab,'%d'))
    }

    if('tide'%in%names(ds) && sum(!is.na(ds$tide)) > 0){
        # Tidal plot
        plot(ds$time, ds$tide,t='l', 
             col='blue', lend=1, 
             ylab='Tide (m)', main=paste0('Tide (with astronomical tidal residual): ', year))
        if(!is.null(event_data)){
            for(i in 1:length(event_data)){
                points(event_data[[i]]$time, event_data[[i]]$tide, t='l', col='red')
                if(add_event_start_lines) abline(v = as.numeric(event_data[[i]]$time[1]), col='pink', lty='solid')
            }
        }

        points(ds$time, ds$tideResid, t='l', col='pink')
        grid(col='brown')


    }    
    if(add_days){
        axis(side=1, at=days_tick, labels=FALSE, tick=TRUE)
        axis(side=1, at=as.numeric(days_lab), labels=format(days_lab,'%d'))
    }

    return(invisible())
}



###############################################################################
#'
#' Check correlations of Hsig between different sites, by year
#' @param year The year to compare
#' @param site1 A site in the wd list
#' @param site2 Another site in the wd list
#' @param wd output of parse_MHL_wave_buoy_data
#' @param variable name (which variable to correlate)
#' @param site1_restriction. Logical vector with the same length as
#'        wd[[site1]]$year. Comparison only occurs where it has the value TRUE
check_station_correlations<-function(year, site1, site2, wd, variable_name='hsig', 
    site1_restriction = TRUE){
   
    # Find indices corresponding to 'year' at site1
    year_indices = which( (floor(wd[[site1]]$year) == year)&(site1_restriction))

    if(length(year_indices)==0){
        plot(c(0,1),c(0,1), main='Insufficient data')
        return(NA)
    }

    # Find matching indices at site2
    matching_indices = match(wd[[site1]]$year[year_indices], wd[[site2]]$year)
    
    if(all(is.na(year_indices))){
        plot(c(0,1),c(0,1), main='Insufficient data')
        return(NA)
    }

    
    hsig_site1 = wd[[site1]][[variable_name]]
    hsig_site2 = wd[[site2]][[variable_name]]
    
    plot(hsig_site1[year_indices], hsig_site2[matching_indices],
        xlab=paste0(site1, ': ', variable_name),
        ylab=paste0(site2, ': ', variable_name),
        main = year)

    local_cor = cor(hsig_site1[year_indices], hsig_site2[matching_indices], 
        use='pairwise.complete')
    title(sub=paste0('Correlation coefficient: ', round(local_cor,2)), col.sub='red')

    abline(0,1,col='red')
    grid()

    return(local_cor)
}

#' Given a series which possibly contains NA values, use
#' linear interpolation to replace NA values which are in 
#' a group of less than n values.
#'
#' e.g. Consider:
#'  series = c(1,4,5,NA,7,NA,NA,5,4,3,NA,NA,NA, 8)
#' 
#' Calling the function with n = 3 will replace the first 3 NAs, but not the remaining
#' ones which are in a group with size >= 3. Calling the function with n=1 will
#' only replace the first NA.
#'
#' Note that NA values outside the bounds of the nonNA series will stay as NA, because
#' we cannot interpolate them!
#'
interpolate_NA_runs_with_length_less_than_n<-function(series, n){

    # Cannot interpolate with only 2 non-missing values
    if(sum(is.na(series)) > (length(series) - 2)) stop('Too many NA values to gap fill')
    
    # Find runs of missing values 
    missing_values = is.na(series)
    # Runs of missing values have 'values = TRUE'
    missing_values_runs = rle(missing_values)
   
    # Set runs of length > n to value = FALSE.
    # The remaining TRUE values need interpolation 
    interpolate_over = missing_values_runs
    interpolate_over$values[interpolate_over$lengths >= n] = FALSE
    interpolate_over_logic = inverse.rle(interpolate_over)

    if(any(interpolate_over_logic)){
        sn = 1:length(series)
        interp_inds = sn[which(interpolate_over_logic)]
        
        rn = sn[which(!interpolate_over_logic)]

        interp_out = approx(rn, series[rn], xout = interp_inds)
        
        series[interp_inds] = interp_out$y
    }

    return(series)
}

#' Quick test code for the above function
.test_interpolate_NA_runs_with_length_less_than_n<-function(){
    
    series = c(1,4,5,NA,7,NA,NA,4,4,3,NA,NA,NA, 7)
    
    s1 = interpolate_NA_runs_with_length_less_than_n(series, n=1)
    stopifnot(identical(s1, series))

    s1 = interpolate_NA_runs_with_length_less_than_n(series, n=2)

    stopifnot(all.equal(s1[4], 6))
    stopifnot(all(which(is.na(s1)) == c(6,7, 11, 12, 13)))

    s1 = interpolate_NA_runs_with_length_less_than_n(series, n=3)
    stopifnot( all.equal(s1[4], 6))
    stopifnot( all.equal(s1[6:7], c(6,5)))
    stopifnot(all(which(is.na(s1)) == c(11,12,13)))
    
    s1 = interpolate_NA_runs_with_length_less_than_n(series, n=4)
    stopifnot( all.equal(s1[4], 6))
    stopifnot( all.equal(s1[6:7], c(6,5)))
    stopifnot(all.equal(s1[11:13], c(4,5,6)))
    
}


#################################################################################
#'
#' Replace missing data at one station with another station
#'
#' @param desired_times Vector of strptime type with desired times for the output data.frame
#' @param site_preference_order Character vector giving the preference order for the output data
#' @param wd Output from parse_MHL_wave_buoy_data
#' @param use_interpolation_for_gaps_less_than integer Once we have filled missing values with data
#' from one station, we check for remaining gaps of <= 'use_interpolation_for_gaps_less_than' sequential
#' values. We then fill those gaps with interpolation of the current data. The idea is that for 'short' gaps, 
#' it is better to interpolate, rather than use data from other stations. Note that wave direction is rounded
#' to the nearest integer for consistency with other values
#' @return A data.frame containing data from each site in order of preference.
gap_fill_wave_data<-function(desired_times, site_preference_order, wd, 
    use_interpolation_for_gaps_less_than=0){

    site1 = site_preference_order[1]
    ncol = length(wd[[site1]][1,])

    # Make output data.frame
    output = matrix(NA, ncol=ncol, nrow = length(desired_times))
    output = as.data.frame(output)
    names(output) = names(wd[[site1]][1,])

    # Store the site too
    output_waves_site = rep(NA,length(desired_times))
    output_dir_site = rep(NA,length(desired_times))

    for(site in site_preference_order){

        #
        # Find rows where hsig is missing, and fill in those we can
        #

        # Find rows that need filling
        missing_waves = which(is.na(output$hsig))

        # Find the indices of the same times at site
        matches = match(as.character(desired_times[missing_waves]),
            as.character(wd[[site]]$time))

        # Fill cols 2:ncol (everything except time)
        output$hsig[missing_waves] = wd[[site]]$hsig[matches] 
        output$hmax[missing_waves] = wd[[site]]$hmax[matches] 
        output$tz[missing_waves] = wd[[site]]$tz[matches] 
        output$tp1[missing_waves] = wd[[site]]$tp1[matches] 
        output$year[missing_waves] = wd[[site]]$year[matches]

        # Store the site name too
        output_waves_site[missing_waves] = site

        #
        # Interpolate over runs of NA's with length 1, 2, or 3
        #
        n = use_interpolation_for_gaps_less_than
        output$hsig = interpolate_NA_runs_with_length_less_than_n(output$hsig, n=n)
        output$hmax = interpolate_NA_runs_with_length_less_than_n(output$hmax, n=n)
        output$tz = interpolate_NA_runs_with_length_less_than_n(output$tz, n=n)
        output$tp1 = interpolate_NA_runs_with_length_less_than_n(output$tp1, n=n)
        output$year = interpolate_NA_runs_with_length_less_than_n(output$year, n=n)

        #
        # Find rows where wave dir is missing, and fill in those we can
        #

        missing_dir = which(is.na(output$dir))

        # Find the indices of the same times at site
        matches = match(as.character(desired_times[missing_dir]),
            as.character(wd[[site]]$time))

        output$dir[missing_dir] = wd[[site]]$dir[matches]

        newdir = interpolate_NA_runs_with_length_less_than_n(output$dir, n=n)
        # For direction, round to integer
        newdir = round(newdir)
        output$dir = newdir

        # Store the site name too
        output_dir_site[missing_dir] = site

        #
        # Interpolate over gaps of 3 or less hours
        #

    }

    output$time = desired_times

    # Append the site names as a final column to the output
    output$waves_site = output_waves_site
    output$dir_site = output_dir_site

    return(output)

}


#' Quick code to plot a single event set
#' 
#' @param event_data a subset of the full data.frame which contains an event
plot_single_storm_event<-function(event_data){

    event_example = event_data

    par(mfrow=c(4,1))
    par(oma=c(0,0,1,0))

    plot(event_example$time, event_example$hsig, t='h', xlab='Time', 
        ylab=bquote(H[sig]~'(m)'), col='yellow')
    points(event_example$time, event_example$hsig, t='l')
    grid()

    plot(event_example$time, event_example$tp1, t='h', xlab='Time', 
        ylab='TP1 (s)', col='blue')
    points(event_example$time, event_example$tp1, t='l')
    grid()

    plot(event_example$time, event_example$dir, t='h', xlab='Time',
        ylab='Direction', col='green')
    points(event_example$time, event_example$dir, t='l')
    grid()

    plot(event_example$time, event_example$tide, t='l', xlab='Time',
        ylab='Tide (m above MSL)', col='blue')
    points(event_example$time, event_example$tideResid, t='l', col='red')
    grid()

    title(main=paste0('Wave event starting ', event_example$time[1]), outer=TRUE, 
        cex.main=1.4)
}



#' Get event start and end times from wd timeseries.  
#'
#' Events are defined as having hsig > hsig_threshold, having a duration >
#' duration_threshold, and a spacing of at least duration_gap_hours since the last event.
#'
#' @param site_data data.frame with the site data
#' @param hsig_threshold 
#' @param duration_threshold_hours event duration threshold in hours
#' @param duration_gap_hours minimum gap between distinct events in hours. Must
#'        be >= smallest time increment in the data.
#' @param events_to_combine If not NULL, then a list giving indices of initially defined events to merge.
#'        The idea is that you initially run the code with this set to NULL, then manually identify
#'        the events to merge (preliminary event IDs are just 1,2,3,4.... in order).
#'        Then you re-run this function with the same settings, but passing events_to_combine, and
#'        the code will combine the events
#' @return A list with the start_indices , end_indices for each event (indices
#'         correspond to wd$site), and a list named data, containing the sub-set
#'         of the wd$site data for each event
extract_events<-function(
    site_data, 
    hsig_threshold, 
    duration_threshold_hours, 
    duration_gap_hours,
    events_to_combine=NULL
    ){

    # Get all 'runs' of < threshold and > threshold
    high_waves = which(site_data$hsig > hsig_threshold)

    lhw = length(high_waves)

    # Time between high-wave data records
    high_wave_gap = as.numeric(
        difftime(site_data$time[high_waves[2:lhw]], 
                 site_data$time[high_waves[1:(lhw-1)]], units='hours')
        )

    # 'Events' start following a sufficiently long time since the previous
    # exceedence, and end similarly
    event_start_index = high_waves[ c(1, which(high_wave_gap >
        duration_gap_hours) + 1) ]
    event_end_index = high_waves[ c(which(high_wave_gap > duration_gap_hours),
                               length(high_waves))]

    # At this stage, events are defined by event_start_index, event_end_index,
    # and we know that between events there is a 'time gap' of > duration_gap_hours

    # Since start & end occur when Hsig>threshold, can have events with 'zero'
    # duration
    event_duration_hours = as.numeric(difftime(site_data$time[event_end_index],
        site_data$time[event_start_index], units='hours'))

    les = length(event_start_index)
    event_gap_hours = as.numeric(difftime(
        site_data$time[event_start_index[2:les]],
        site_data$time[event_end_index[1:(les-1)]], units='hours'))

    # Remove events which are not long enough
    remove_events = which(event_duration_hours < duration_threshold_hours)

    if(length(remove_events)>0){
        event_start_index = event_start_index[-remove_events]
        event_end_index = event_end_index[-remove_events]
    }

    # Manual merge here
    if(!is.null(events_to_combine)){
        print('Combining events based on manual definitions')
        start_index_removal = c()
        end_index_removal = c()
        for(i in 1:length(events_to_combine)){
            ll = length(events_to_combine[[i]])
            stopifnot(all(diff(events_to_combine[[i]]) == 1))
            start_index_removal = c(start_index_removal, events_to_combine[[i]][2:ll])
            end_index_removal = c(end_index_removal, events_to_combine[[i]][1:(ll-1)])
        }

        event_start_index = event_start_index[-start_index_removal]
        event_end_index = event_end_index[-end_index_removal]
    }

    output = list(start_index = event_start_index,
                  end_index = event_end_index) 

    # Append data to the output
    output$data = list()
    for(i in 1:length(event_start_index)){
        output$data[[i]] = site_data[event_start_index[i]:event_end_index[i],]
    }

    return(output)
}


#' Convenience function to extract the event stats
#'
#' @param event_data a list containing event data frames, like in the $data
#'        entry of extract_events
#' @param median_tp1_dir Logical. If TRUE, tp1 and dir are equal to the median
#'        storm value. Otherwise the value at peak hsig is used.
#' @param dir_min Minimum value for direction. If the event wave direction is smaller
#'        than this, we add 360 to it. This can help avoid an arbitrary split in 
#'        direction data for sites which receive waves form the north.
#' @param duration_offset_hours. The event duration = (end_time - start_time + duration_offset_hours). 
#'        This means events with only 1 data point have duration=duration_offset_hours
#' @return a data.frame containing summary statistics for each event
extract_event_statistics<-function(event_data, median_tp1_dir=TRUE, dir_min=0, 
    duration_offset_hours=0){

    event_statistics = list()

    for(i in 1:length(event_data)){
        data = event_data[[i]]

        if(dim(data)[1] > 1){
            dt = difftime(max(data$time), min(data$time), units='hours') + duration_offset_hours
        }else{
            dt = duration_offset_hours
        }
        hsig_max = max(data$hsig, na.rm=T)

        if(median_tp1_dir){
            tp1_median = median(data$tp1, na.rm=T)
            dir_median = median(data$dir, na.rm=T)
        }else{
            # Use values at peak hsig (use first if the maxima is non unique)
            hsig_ml = which(data$hsig == hsig_max)[1]
            tp1_median = data$tp1[hsig_ml]
            dir_median = data$dir[hsig_ml]
        }

        if(!is.na(dir_median)){

            if(dir_median < dir_min) dir_median = dir_median + 360
        }
        

        # If we don't have any tidal data, just set the value to NA
        if(('tideResid' %in% names(data)) && sum(!is.na(data$tideResid)) > 0){
            tidalresid_max = max(data$tideResid, na.rm=T)
        }else{
            tidalresid_max = NA
        }
        start_time = DU$time_to_year(data$time[1])
        end_time = DU$time_to_year(max(data$time))

        event_statistics[[i]] = c(as.numeric(dt), hsig_max, tp1_median, 
            dir_median, tidalresid_max, start_time, end_time)
    } 

    # Convert to data.frame for output
    output = matrix(unlist(event_statistics), ncol=7, byrow=T)
    output = data.frame(duration=output[,1], hsig = output[,2],
                        tp1 = output[,3], dir = output[,4], 
                        tideResid = output[,5], startyear = output[,6], 
                        endyear=output[,7])

    # Append a strptime time object for the start time (can be useful)
    timevec = rep(event_data[[1]]$time[1], length(output[,1]))
    for(i in 1:length(event_data)){
        timevec[i] = event_data[[i]]$time[1]
    }
    output$time = timevec
    return(output)
}


#'
#' Function to extract clustered event sets from a simulated timeseries
#'
#' See instead the function event_cluster_counter (in cluster_counting_functions),
#' which is better tested
#'
#' @param synthetic_attr data.frame or list. The event set with attributes startyear, hsig,
#'        duration, ... etc
#' @param N integer. The number of events in the cluster
#' @param time_envelope numeric. The maximal time over which the clustered event set
#'        can occur (units of years)
#' @param thresholds list. A list with names corresponding to column names of
#'        synthetic_attr, and values corresponding to thresholds that the
#'        attribute must exceed to be counted in the cluster.
#' @param time_envelope_type character. 'start-to-end' or 'end-to-start'. In the
#' former case 'time-envelope' is measured from the 'start' of the first event until the 'end'
#'  of the last event.
#' @return a 2 column matrix: First column gives the starting index of each
#'         cluster, 2nd column gives the ending index
#'
event_cluster_indices_nonindependent<-function(
    synthetic_attr, 
    N = 2, 
    time_envelope = 60/365.25, 
    thresholds = list(hsig=0, duration=0/365.25),
    time_envelope_type='start-to-end'){

    # Check input
    stopifnot(time_envelope_type %in% c('start-to-end', 'end-to-start'))
    stopifnot(N >= 1)
    stopifnot(time_envelope >= 0)

    # Intialise events_to_keep --- then update later
    events_to_keep = rep(1, length(synthetic_attr[,1]))
    
    # Ensure only events exceeding the thresholds stay in events_to_keep   
    threshold_names = names(thresholds)
    if(length(threshold_names) > 0){
        for(nm in threshold_names){
            stopifnot(nm %in% names(synthetic_attr))
            events_to_keep = events_to_keep*(synthetic_attr[[nm]] >= thresholds[[nm]])
        }
    }

    # Ensure we have > 0 events
    le = sum(events_to_keep)
    if(le < N) stop('No events matching these criteria')

    # Compute the time for the next N events matching the criteria
    events_to_keep_inds = which(events_to_keep == 1)
    clusters0  = synthetic_attr$startyear[events_to_keep_inds]
    durations0 = synthetic_attr$duration[events_to_keep_inds]

    # Time interval to have N events
    if(time_envelope_type == 'start-to-end'){
        # Duration goes from start of first event to end of last event
        duration_to_n_events = ( (clusters0[N:le] + durations0[N:le]) - 
            clusters0[1:(le-(N-1))])

    }else if(time_envelope_type == 'end-to-start'){
        # Duration goes from end of first event to start of last event

        duration_to_n_events = ( clusters0[N:le] - 
            (clusters0[1:(le-(N-1))] + durations0[1:(le-(N-1))]))

    }else{

        stop(paste0('Unrecognized value of "time_envelope_type = "', time_envelope_type))
    }

    # Apply constraint
    clusters1 = which( duration_to_n_events <= time_envelope)

    # Return clusters (which might overlap)
    cluster_start = events_to_keep_inds[clusters1]
    cluster_end = events_to_keep_inds[clusters1 + (N-1)]

    # Record whether the event overlaps with the previous or next event
    lc = length(cluster_start)
    overlap_with_next = c((cluster_end[1:(lc-1)] >= cluster_start[2:lc]), 0)
    #non_independence_flag = 
    #non_independence_flag = (c(0, non_independence_flag) + c(non_independence_flag, 0))>0

    return(cbind(cluster_start, cluster_end, overlap_with_next))
}


#' Convenience function to add a histogram to a 'pairs' plot on the diagonal
#' See the examples in ?pairs
#' @param x A numeric vector
panel.hist <- function(x, ...){
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(usr[1:2], 0, 1.5) )
         h <- hist(x, plot = FALSE)
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         rect(breaks[-nB], 0, breaks[-1], y, col = "brown", ...)
}


#' Convenience function to plot a median smoother on the data, assuming it's periodig
periodic_annual_smooth<-function(time_of_year, value, k=11, smoothcol='black', ...){

    stopifnot(max(time_of_year) <= 1)
    stopifnot(min(time_of_year) >= 0)
    time_of_year_range = range(time_of_year)

    # Sort inputs
    time_of_year_order = sort(time_of_year, index.return=TRUE)$ix
    time_of_year = time_of_year[time_of_year_order]
    value = value[time_of_year_order]

    periodic_time = c(time_of_year, time_of_year+1, time_of_year+2)
    periodic_var = c(value, value, value)

    var_na = which(is.na(periodic_var))
    if(length(var_na) > 0){
        periodic_time = periodic_time[-var_na]
        periodic_var = periodic_var[-var_na]
    }

    smooth00 = runmed(periodic_var, k=k)
    plot(time_of_year, value, ...)
    points(periodic_time-1, smooth00, t='o', col=smoothcol)

}


#' Make an empirical quantile function, and inverse quantile function
#'
#' The function uses R's "density" to evaluate the density x,y values (being careful to
#' normalise so it integrates to 1). The user can provide from/to arguments to
#' determine the domain of the density, arguments to control the smoother, etc.
#' See ?density for information on those options. 
#'
#' Once we have computed the density, we integrate it with the trapozoidal rule, 
#' and fit linear approxfuns [quantile function takes integral(y) to x, 
#' inverse takes x to integral(y)]. Note that by using linear approxfuns, we can ensure
#' that the inverse relation holds.
#' 
#'
#' @param data a vector of data from which the empirical quantile function will be made
#' @param ... further arguments to density
#' @return a list with the quantile and inverse quantile functions
#'
make_empirical_quantile_and_inverse_quantile_functions<-function(data, ...){
    
    data_dens = density(data, ...)

    # Trapezoidal integral of density
    l = length(data_dens$y)
    dx = data_dens$x[2] - data_dens$x[1]
    dens_integral = cumsum(( c(0, data_dens$y[2:l]) + c(0, data_dens$y[1:(l-1)])))*dx*0.5

    quant_q = data_dens$x
    quant_p = dens_integral/max(dens_integral)
    stopifnot(quant_p[1]==0)
    stopifnot(all.equal(quant_p[l],1))

    # Quantile function
    qempirical = approxfun(quant_p, quant_q) 
    # Inverse quantile function
    pempirical = approxfun(quant_q, quant_p)

    return(list(qempirical = qempirical, 
                pempirical = pempirical))
}


.test_make_empirical_quantile_and_inverse_quantile_functions<-function(){

    # Some data
    data = rnorm(300)

    # Make the quantile / inverse quantile functions
    tmp = make_empirical_quantile_and_inverse_quantile_functions(data)
    qemp = tmp$qempirical
    pemp = tmp$pempirical

    # Logical checks on probabilities
    probs = pemp(data)
    stopifnot(min(probs)>=0)
    stopifnot(max(probs)<=1)

    # Check that quantile really is the inverse of probability
    q_probs = qemp(probs)
    stopifnot(all.equal(data, q_probs))

    print('PASS')

}

test_all<-function(){
    .test_make_empirical_quantile_and_inverse_quantile_functions()
    .test_year_to_time()
    .test_time_to_year()
    .test_is_leap_year()
}

#'
#' Read the climate indices files
#' Keep this in a function to avoid clutter elsewhere
#'
read_climate_indices<-function(soi_file, aao_file){
    
    CI = list() # For outputs

    #
    # Read SOI into a 2-column (time, value) format
    #
    soi_temp = read.table(soi_file, skip=6, na.strings='*', header=TRUE) 
    soi_values = soi_temp[,2:13]
    soi_years = rep(soi_temp[,1], each=12, len=12*length(soi_temp[,1]))
    soi_months = rep(month.name, length(soi_temp[,1]))
    # Assign the time to the 15th of the month
    soi_time = strptime( paste(soi_years, soi_months, '15', sep="-"), 
        format="%Y-%B-%d", tz='Etc/GMT')
    CI$soi = data.frame(time=soi_time, index = c(t(as.matrix(soi_values))) )
    rm(soi_temp, soi_values, soi_years, soi_months, soi_time)

    # The definition of SOI involves multiplication by 10. Let's undo that so
    # the scale is more similar to the other series.
    #CI$soi[,2] = CI$soi[,2]/10 


    ##
    ## Read Nino3.4 into a 2-column (time-value) format
    #nino34_temp = read.table(nino34_file, skip=1, header=TRUE)
    #nino34_values = nino34_temp$NINO3.4
    ## Normalise
    #nino34_values = (nino34_values - mean(nino34_values))/sd(nino34_values)
    #nino34_time = strptime(
    #    paste0(nino34_temp$YR, '-', nino34_temp$MON, '-', 15, sep=""),
    #    format='%Y-%m-%d', tz='Etc/GMT')
    #CI$nino34 = data.frame(time = nino34_time, index = nino34_values)
    

    ##
    ## Read PDO into a 2-column (time, value) format
    ##
    #pdo_temp = read.table(pdo_file, skip=30, na.strings='*', header=FALSE, 
    #    nrows = 116, fill=TRUE) 
    #pdo_values = pdo_temp[,2:13]
    #pdo_years = rep(gsub('*', '', as.character(pdo_temp[,1]), fixed=TRUE),
    #    each = 12, len=12*length(pdo_values[,2]))
    #pdo_months = rep(month.name, len=12*length(pdo_values[,2]))
    ## Assign pdo to the 15th of the month
    #pdo_time = strptime( paste(pdo_years, pdo_months, '15', sep="-"), 
    #    format='%Y-%B-%d', tz='Etc/GMT')
    #CI$pdo = data.frame(time = pdo_time, index = c(t(as.matrix(pdo_values))))
    #rm(pdo_temp, pdo_values, pdo_years, pdo_months, pdo_time)

    #
    # Read Antarctic Oscillation into a 2-column (time, value) format
    #
    aao_temp = read.table(aao_file, header=FALSE, skip=1, fill=TRUE) 
    aao_values = aao_temp[,2:13]
    aao_years = rep(aao_temp[,1], each=12, len=12*length(aao_values[,1]))
    aao_months = rep(month.name, len=12*length(aao_values[,1]))
    # Assign aao to the 15th of the month
    aao_time = strptime( paste(aao_years, aao_months, '15', sep="-"), format='%Y-%B-%d', tz='Etc/GMT')
    CI$aao = data.frame(time=aao_time, index=c(t(as.matrix(aao_values))))
    rm(aao_temp, aao_values, aao_years, aao_months, aao_time)

    return(CI)
}



read_Adelaide_wave_data<-function(
    Adelaide_wave_datasets = Adelaide_wave_datasets,
    Adelaide_wave_files = Adelaide_wave_files){

    library(ncdf4)

    # Store all datasets in lists
    wd = list() # wave data

    for(i in 1:length(Adelaide_wave_datasets)){

        site = Adelaide_wave_datasets[i]

        # Read wave dataset .nc
        wave_file = nc_open(Adelaide_wave_files)

        time_Julian = ncvar_get(wave_file, varid = 'time')
        date_base = as.POSIXct('1990-01-01 00:00:00', tz = 'UTC')
        date_base_offset = 9.5*60*60 # Use this to convert from UTC to ACST
        date_time = format(date_base + time_Julian*24*60*60 + date_base_offset, 
            "%Y-%m-%d %H:%M:%S")
        #date_time = format(date_base + time_Julian*24*60*60, "%Y-%m-%d %H:%M")

        hsig = as.numeric(ncvar_get(wave_file, varid = 'hs'))
        tm0 = as.numeric(ncvar_get(wave_file, varid = 'tm01'))
        tp1 = as.numeric(1/ncvar_get(wave_file, varid = 'fp'))
        dir = as.numeric(ncvar_get(wave_file, varid = 'dir'))

        wave_data = data.frame(#time_Julian = time_Julian,
                               time = date_time,
                               hsig = hsig,
                               #tm0 = tm0,
                               tp1 = tp1,
                               dir = dir,
                               stringsAsFactors=FALSE)

         
        wave_data$tp1[is.infinite(wave_data$tp1)] <- NA
        # wave_data$dir[wave_data$dir<150] <- 150

        # Convert time to date-time object -- Australian Central Standard Time without daylight.
        time_data = strptime(wave_data$time,
                            format = '%Y-%m-%d %H:%M:%S',
                             tz = 'Australia/Darwin')

        # We are sometimes 1 second off. Fix that here
        time_data = round(time_data, 'mins')

        # Also store time as a decimal year -- useful
         year_data = DU$time_to_year(time_data)
        # Append the time as a decimal year to wave_data
         wd[[site]] = cbind(wave_data, year=year_data)

        # Replace with the date-time object
         wd[[site]]$time = time_data
    }
    # Put these variables in the global environment
    return(wd)
}


#' Parse the Adelaide tidal gauge data
#' 
#' @param read_Adelaide_csv_tide_gauge
#' @return A data.frame containing the time, julian time (days since 1970), tidal level, and status
read_Adelaide_csv_tide_gauge<-function(gauge_file){

    JULIAN_TIME_ORIGIN = strptime("1970-01-01 00:00:00", format='%Y-%m-%d %H:%M:%S', tz = "Australia/Darwin")
    Adelaide = read.csv(gauge_file, skip=0,header=T,
        na.strings=c('NA', 'n/a', '-9.999'), stringsAsFactors=FALSE)
    Adelaide_time = strptime(paste(Adelaide[,1], Adelaide[,2]),
        format='%d/%m/%Y %H:%M:%S', tz='Australia/Darwin')
    Adelaide_julian_time = julian(Adelaide_time,
        origin = JULIAN_TIME_ORIGIN)

    output = data.frame(time=Adelaide_time,
                        julian_time = as.numeric(Adelaide_julian_time),
                        tide=Adelaide[,4], status=Adelaide[,5])
    return(output)
}

read_Adelaide_csv_tide_pred<-function(prediction_file){
    JULIAN_TIME_ORIGIN = strptime("1970-01-01 00:00:00", format='%Y-%m-%d %H:%M:%S', tz = "Australia/Darwin")
    Adelaide = read.csv(prediction_file, skip=0, header=T,
        na.strings=c('NA', 'n/a'), stringsAsFactors=FALSE)
    Adelaide_time = strptime(Adelaide[,1],
        format='%d/%m/%Y %H:%M', tz='Australia/Darwin')
    Adelaide_julian_time = julian(Adelaide_time,
        origin = JULIAN_TIME_ORIGIN)

    output = data.frame(time=Adelaide_time,
                        julian_time = as.numeric(Adelaide_julian_time),
                        tide=Adelaide[,2])
    return(output)
}

#'
#' Nice looking alternative to 'pairs'
#'
nice_pairs<-function( mydata, extra_data=NULL ){
    # Prettier pairs plot
    pairs(mydata,
        pch='.', cex=3, 
        upper.panel=function(x,y,...){ 
            points(x,y, col='green', ...); 
            grid(col='orange')
            if(!is.null(extra_data)){
                # Hack to identify the data columns
                icol = get('i', pos=parent.frame(n=2))
                jcol = get('j', pos=parent.frame(n=2))
                points(extra_data[,jcol], extra_data[,icol], col='black', pch='.', cex=2)
            }

            spearman_cortest = try(suppressWarnings(cor.test(x, y, method='s')))
            if(class(spearman_cortest) != 'try-error'){
                spearman_cor = spearman_cortest$estimate
                title_word = as.character(round(spearman_cor, 3))
                if(spearman_cortest$p.value < 0.05){
                    title_word = paste0(title_word, '*')
                }
                if(spearman_cortest$p.value < 0.005){
                    title_word = paste0(title_word, '*')
                }

                title(main=title_word, line=-1, col.main='blue', cex.main=1.5) 
            }
        }, 
        diag.panel=DU$panel.hist,
        lower.panel=function(x,y,...){ 
            points(x,y, col=0); 

            spearman_cortest = try(suppressWarnings(cor.test(x, y, method='s')))
            if(class(spearman_cortest) != 'try-error'){
                spearman_cor = spearman_cortest$estimate
                title_word = as.character(round(spearman_cor, 3))
                font_size = 1.0
                if(spearman_cortest$p.value < 0.05){
                    title_word = paste0(title_word, '*')
                    font_size = font_size + 0.5
                }
                if(spearman_cortest$p.value < 0.005){
                    title_word = paste0(title_word, '*')
                    font_size = font_size + 0.5
                }

                title(main=title_word, line=-2, col.main='red', cex.main=font_size) 
            }
        })
}

#' Improved alternative to qqplot
qqplot3<-function(x, y, plot.it = TRUE, xlab = deparse(substitute(x)), 
    ylab = deparse(substitute(y)), ...){
    sx <- sort(x)
    sy <- sort(y)
    lenx <- length(sx)
    leny <- length(sy)
    if (leny < lenx) 
        sx <- quantile(sx, p = (1:leny)/(leny+1), type=6) #approx(1L:lenx, sx, n = leny)$y
    if (leny > lenx) 
        sy <- quantile(sy, p = (1:lenx)/(lenx+1), type=6) #approx(1L:leny, sy, n = lenx)$y
    if (plot.it) 
        plot(sx, sy, xlab = xlab, ylab = ylab, ...)
    invisible(list(x = sx, y = sy))
}

#' Copy an environment (rather than just pass by reference)
#'
#' Currently the recursive option will not work (only relevant if the
#' environment contains an environment), though that should be an easy fix
#'
#' @param fit_env An environment
#' @param recursive If the environment contains other environments, should they
#' be deep copied as well?
#' @return A copy of fit_env
#'
clone_environment<-function(env, recursive=TRUE){

    if(any(eapply(env, class) == 'environment')){
        if(recursive){
            stop('Recursive environment clone not yet implemented')
        }
    }

    return(list2env(as.list(env, all.names=TRUE)))
}

