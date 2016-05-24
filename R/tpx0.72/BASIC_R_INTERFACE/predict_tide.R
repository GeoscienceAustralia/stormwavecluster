## Crude R interface to use of predict_tide
## Uses system calls, but takes care of creating input files, calling, etc.

## The actual case specific code is not here -- but these routines can be used to do the work

# This should be the directory containing the 'predict_tide' program
#.OTPS_directory = '/home/gareth/Code_Experiments/TIDAL_PREDICTION/TPX072/OTPS'
.OTPS_directory = normalizePath('../OTPS')

###################################################################################
#
# BELOW HERE THINGS SHOULDN'T BE CHANGED
#
###################################################################################

#' Make the lat_lon_time file for predict_tide
make_lat_lon_time_file<-function(site_coordinates, prediction_times, 
    output_file = NA, return_output = FALSE, verbose=TRUE){

    if(verbose) print('Making lat_lon_time_file...')
   
    if(return_output == TRUE){ 
        
        stop('This option is broken, FIXME')

        lats = rep(site_coordinates[2], len=length(years))
        lons = rep(site_coordinates[1], len=length(years))
        years = format(prediction_times, '%Y') 
        months = format(prediction_times, '%m')
        days = format(prediction_times, '%d')
        hours = format(prediction_times, '%H')
        minutes = format(prediction_times, '%M')
        seconds = format(prediction_times, '%S')

        output_data = data.frame(lats, lons, years, months, days, hours, minutes, 
            seconds, stringsAsFactors=FALSE)

        if(!is.na(output_file)){
            write.table(output_data, file=output_file, sep=" ", quote=FALSE, 
                row.names=FALSE, col.names=FALSE)
        }

    }else{
        # Low memory version
        lp = length(prediction_times)

        time_info = gsub("-", " ", prediction_times)
        time_info = gsub(":", " ", time_info)
        #browser()
        output_data = paste(
            rep(site_coordinates[2], len=lp),
            rep(site_coordinates[1], len=lp),
            time_info,
            #format(prediction_times, '%Y'),
            #format(prediction_times, '%m'),
            #format(prediction_times, '%d'),
            #format(prediction_times, '%H'),
            #format(prediction_times, '%M'),
            #format(prediction_times, '%S'), 
            sep=" ")

        if(!is.na(output_file)){
            cat(output_data, file=output_file, sep="\n")
        }
    }

    if(return_output) return(output_data)

    return()
}


#' Make the setup.inp file for predict_tide
make_setup.inp<-function(lat_lon_time_file, output_tide_file, setup_file, 
    verbose=TRUE){
    if(verbose) print('Making setup.inp ...')
    # First line = model data source
    file_lines = 'DATA/Model_tpxo7'
    # Second = file with lat/lon/time
    file_lines = c(file_lines, basename(lat_lon_time_file))
    # Third = What to get
    file_lines = c(file_lines, 'z')
    # 4th/5th = blank
    file_lines = c(file_lines, '', '')
    # 5th = height datum
    file_lines = c(file_lines, 'oce')
    # 6th = option which we must set to zero
    file_lines = c(file_lines, '0')
    # 7th = output file
    file_lines = c(file_lines, basename(output_tide_file))

    # Write the output
    cat(file_lines, file=setup_file, sep="\n")
}

#' Make the command line call
run_predict_tide<-function(setup_file, OTPS_directory, verbose=TRUE){

    if(verbose) print('Calling predict_tide ... (its output printed below) ...')
    # Go back to the current directory when this function finishes
    mydir = getwd()
    on.exit(setwd(mydir))

    setwd(OTPS_directory)

    # Run the predict_tide program
    system_command = paste0('./predict_tide <', basename(setup_file))
    system(system_command)

}

#' Get rid of temp files used to run predict_tide
cleanup_intermediate_files<-function(lat_lon_time_file, output_tide_file, 
    setup_file){
    file.remove(lat_lon_time_file)
    file.remove(setup_file)
    file.remove(output_tide_file)
}


#' Determine if a year is a leap year
is_leap_year<-function(year){

    (year%%4 == 0)*( (year%%100 != 0) + (year%%100 == 0)*(year%%400 == 0))
}

TEN_THOUSAND_DAYS = sum(is_leap_year(0:9999))*366 + sum(!is_leap_year(0:9999))*365

#' Timezone conversion + deal with years > 9999
#'
#' @param start_time time to convert (strptime object)
#' @param tz Timezone (e.g. 'Etc/GMT-10')
#' @return Strptime object with start_time converted to the new timezone
#'
convert_timezone<-function(start_time, tz){
    
        time_format_string = '%Y-%m-%d %H:%M:%S'

        year = as.numeric(format(start_time, '%Y'))

        if(year >= 1e+04){

            extra_ten_thousands = floor(year/1e+04)

            extra_days = as.difftime(extra_ten_thousands*TEN_THOUSAND_DAYS, units='days')

            start_time = start_time - extra_days

        }else{

            extra_days = as.difftime(0, units='days')
        }

        start_time_char = format(as.POSIXct(start_time), tz=tz, 
            format=time_format_string)

        new_start_time = strptime(start_time_char, format=time_format_string,
            tz=tz)

        new_start_time = new_start_time + extra_days

        return(new_start_time)
}

#' Given a list containing vectors all of the same time, unpack into a
#' large vector, without problematic coercion
unpack_list_to_vector<-function(big_list){

    lpt = sum(unlist(lapply(big_list, length))) # Length of all
    all_prediction_times = rep(big_list[[1]][1], len=lpt)
    
    counter = 0
    for(i in 1:length(big_list)){
        inds = counter + (1:length(big_list[[i]]))
        all_prediction_times[inds] = big_list[[i]]
        big_list[[i]] = NULL
        counter = inds[length(inds)]
    }

    return(all_prediction_times)
}

###############################################################################
#
# Main program below here
#
###############################################################################

#' Get tides at a particular site using TPX072
#' @param site_name Name for the site (just used for files)
#' @param site_coordinates vector of length 2 giving lon,lat in decimal degrees
#' @param start_time start time (strptime object). Must have correct timezone
#' @param end_time end time (strptime object). Must have correct timezone
#' @param time_interval. Time interval accepted by 'seq' (e.g. '1 hour' or '15 min' or '5 days' or '30 sec')
#' @param OPTS_directory location of the TPX072 OPTS directory
get_tidal_prediction<-function(site_name, site_coordinates, 
    start_time, end_time, time_interval, OTPS_directory=.OTPS_directory){

    #on.exit(browser())

    if(!file.exists(OTPS_directory)){
        msg1 = paste0('Cannot find OTPS directory ', OTPS_directory)
        msg2 = 'You might need to pass the variable OTPS_directory to get_tidal_prediction'
        msg3 = "Alternatively, adjust the variable name '.OTPS_directory' in predict_tide.R"
        print(c(msg1, msg2, msg3))
        stop()
    }

    stopifnot(length(start_time) == length(end_time))

    # Key filenames
    lat_lon_time_file = paste0(OTPS_directory, '/' , 'lat_lon_time_R_', 
        site_name)
    output_tide_file = paste0(OTPS_directory, '/' , 'predictions_', site_name,
        '.out')
    setup_file = paste0(OTPS_directory, '/', 'setup_', site_name, '.inp')


    prediction_times = vector(mode='list', len = length(start_time))
    gmt_prediction_times = vector(mode='list', len = length(start_time))
    tz = attr(start_time[1], 'tzone') #format(start_time[1], '%Z') 
    format_string = '%Y-%m-%d %H:%M:%S'
    GMT_tz = 'Etc/GMT'
    apply_timezone_conversion = (tz != GMT_tz)
    for(i in 1:length(start_time)){

        if(i%%1e+04 == 0) print(i)

        if(apply_timezone_conversion){
            # Compute times to get output
            pred_seq = seq(start_time[i], end_time[i], by=time_interval)
            prediction_times[[i]] = format(pred_seq, '%s') #strftime(pred_seq, format=format_string, tz=tz) #as.character(pred_seq)

            lps = length(pred_seq)

            # Change timezone
            gmt_start_time = convert_timezone(pred_seq[1], tz=GMT_tz)
            gmt_end_time = convert_timezone(pred_seq[lps], tz=GMT_tz)

            gmt_prediction_times[[i]] = strftime(seq(gmt_start_time, gmt_end_time, len=lps), format=format_string, tz=GMT_tz)

            stopifnot(length(prediction_times[[i]]) == length(gmt_prediction_times[[i]]))
        }else{
            # Compute times to get output
            pred_seq = seq(start_time[i], end_time[i], by=time_interval)
            prediction_times[[i]] = format(pred_seq, '%s') #strftime(pred_seq, format=format_string, tz=tz)
        }
    }

    print('Unpack 1 ...')
    # Pack into a single vector
    prediction_times = strptime((unlist(prediction_times)), format='%s', tz=tz) #unpack_list_to_vector(prediction_times)
    #browser()
    print('Unpack 2 ...')
    if(tz == GMT_tz){
        gmt_prediction_times = strftime(prediction_times, format=format_string, tz=GMT_tz)
    }else{
        gmt_prediction_times = unlist(gmt_prediction_times) #unpack_list_to_vector(prediction_times)
    }

    # Make the lat_lon_time file
    make_lat_lon_time_file(site_coordinates, prediction_times=gmt_prediction_times, 
        output_file=lat_lon_time_file)

    # Make setup.inp
    make_setup.inp(lat_lon_time_file, output_tide_file, setup_file)

    # Call predict_tide
    run_predict_tide(setup_file, OTPS_directory)

    # Read the output of predict_tide
    output_text = read.table(output_tide_file, sep="", colClasses='character', 
        stringsAsFactors=FALSE, skip=6)

    output_data = data.frame(time=prediction_times, tide=as.numeric(output_text[,5]))

    cleanup_intermediate_files(lat_lon_time_file, setup_file, output_tide_file)

    return(output_data)
}


