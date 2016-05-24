######################################################################################################
#' Read data from a csv file containing the gauge info
read_csv_gauge<-function(gauge_file){

    JULIAN_TIME_ORIGIN = strptime("1970-01-01 00:00:00", format='%Y-%m-%d %H:%M:%S', tz = "Etc/GMT-10")
    tomaree = read.csv(gauge_file, skip=17,header=T,
        na.strings=c('NA', 'n/a'), stringsAsFactors=FALSE)
    tomaree_time = strptime(paste(tomaree[,1], tomaree[,2]), 
        format='%d/%m/%Y %H:%M:%S', tz='Etc/GMT-10')
    tomaree_julian_time = julian(tomaree_time, 
        origin = JULIAN_TIME_ORIGIN)

    output = data.frame(time=tomaree_time, 
                        julian_time = as.numeric(tomaree_julian_time),
                        tide=tomaree[,3], status=tomaree[,4])
    return(output)
}

coffsh = read_csv_gauge(Sys.glob('Coffs*.csv'))

########################################################################################################
# Get some predictions
source('../predict_tide.R', chdir=TRUE)

site_name = 'Coffs_harbour' # (No spaces in name)
# Site coordinates as decimal degrees (longitude, latitude)
site_coordinates = c(153 + 48/60 + 45.82/(60*60), -(30 + 18/60 + 10.29/(60*60)))

# Timezone MUST be GMT / UTC. Note that -10 means 'ahead 10 hours', 
start_time = coffsh$time[1]
end_time   = coffsh$time[length(coffsh$time)]
time_interval = '15 min'

coffs_pred = get_tidal_prediction(site_name, site_coordinates, 
    start_time, end_time, time_interval)

## Compare

# Times should be identical
stopifnot(all(coffsh$time == coffs_pred[,1]))

# Tidal residual
tidal_residual = coffsh$tide - coffs_pred[,2]
mean_offset = mean(coffsh$tide, na.rm=T)

compute_30day_moving_average = FALSE
if(compute_30day_moving_average){
    # Have seen some analyses which take the difference between
    # measurements and 30 day moving average, and compare that to
    # tidal predictions, so that the longer term sea level anomalies
    # are removed. We can't apply this in the current work though (since we
    # have no way to predict the anomaly separately)

    coffsh_30 = coffsh$tide*NA
    filter_length = 4*24*30/2
    series_length = length(coffsh_30)
    for(i in 1:series_length){
        # ('filter' is faster than this loop, but doesn't treat NA)
        # filter(coffsh$tide, rep(1, 4*24*30)/(4*24*30)) # 30 day running mean
        coffsh_30[i] = mean(coffsh$tide[max(i-filter_length, 1):min(i+filter_length, series_length)], na.rm=T)
    }
    tidal_anomoly = coffsh$tide - coffsh_30 - coffs_pred[,2]
}


# Looks ok?
for( year in 1996:2014){
    #year = 1997
    year_inds = which(format(coffsh$time, '%Y')==as.character(year))

    if(length(year_inds)==0) next
    if(all(is.na(coffsh$tide[year_inds]))) next

    png(paste0('Tidal_plot_', year, '.png'), width=12, height=7, units='in', res=100)
        plot(coffsh$time[year_inds], coffsh$tide[year_inds] - mean_offset, t='l',
            xlab='Time', ylab='Tide (m MSL)', xaxs='i', main=year)
        points(coffs_pred[,1], coffs_pred[,2], t='l',col='red')
        #plot(coffsh$time, tidal_residual-mean_offset,t='l')
        #abline(h=0, col='red')
        points(coffs_pred[,1], tidal_residual - mean_offset,t='l',col=3)
        grid(col='brown')
        legend('bottomright', legend=c('Data', 'Prediction (TPXO72)', 'Residual'), 
            lty=c(1,1,1), col = c('black', 'red', 'green'), bg='white')
    dev.off()
}
