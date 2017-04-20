#
# Read all the RDS files produced when bootstrapping, and compute various
# summary datasets, and export the data.
#


library(methods)
library(evmix)
library(copula)
library(VineCopula)
library(CDVine)

nhp = new.env()
source('../../R/nhpoisp/nhpoisp.R', local=nhp)

DU = new.env()
source('../preprocessing/data_utilities.R', local=DU)

evmix_fit = new.env()
source('../../R/evmix_fit/evmix_fit.R', local=evmix_fit, chdir=TRUE)

wavedisp = new.env()
source('../../R/wave_dispersion/wave_dispersion_relation.R', local=wavedisp, chdir=TRUE)

# 
input_RDS = Sys.glob('Synthetic_series_bootstrap/bootstrap_fit_series_*.RDS')

series_duration_years = 1000

allRDS = vector(mode='list', length = length(input_RDS))
names(allRDS) = input_RDS

for(i in 1:length(input_RDS)) allRDS[[input_RDS[i]]] = readRDS(input_RDS[i])

# Get original fit for comparison [we expect the boostrap uncertainties to
# surround this in some sense]. FORCE IT TO BE THE UNPERTURBED FIT
initial_fit_env = new.env()
load('../statistical_model_fit/Rimages/session_series_simulation_FALSE_0.Rdata', envir=initial_fit_env)

###################################################################################
#
# Export timeseries for coastal modelling
#

export_series_to_csv = FALSE
if(export_series_to_csv){

    csv_dir = 'Synthetic_series_bootstrap_CSV'
    dir.create(csv_dir, showWarnings=FALSE)

    mydir = getwd()

    format_output_data<-function(output_data){
        # Reduce the size for output
        var_keep = c('duration', 'hsig', 'tideResid', 'dir', 'tp1', 'startyear', 'soiA', 'msl')
        output_data = output_data[,var_keep]

        # Reduce significant figures
        output_data$duration = formatC(output_data$duration, digits=4, format='e')
        output_data$hsig = formatC(output_data$hsig, digits=4, format='e')
        output_data$tideResid = formatC(output_data$tideResid, digits=4, format='e')
        output_data$msl = formatC(output_data$msl, digits=4, format='e')
        output_data$dir = formatC(output_data$dir, digits=6, format='e')
        output_data$tp1 = formatC(output_data$tp1, digits=4, format='e')
        output_data$startyear = formatC(output_data$startyear, digits=11, format='e')
        output_data$soiA = formatC(output_data$soiA, digits=4, format='e')

        # Append units to names
        names(output_data) = paste(var_keep, 
            c('(hrs)', '(m)', '(m)', '(deg)', '(s)', '(yr)', '(soi)', '(m)'), 
            sep="")

        return(output_data)
    }
    

    output_compressed_bootstrap<-function(input_RDS_i){
        # Get the synthetic time-series of events
        output_data = allRDS[[ input_RDS_i ]]$synthetic_attr
    
        output_filename = paste0(csv_dir, '/', basename(gsub('RDS', 'csv', input_RDS_i)))

        output_data = format_output_data(output_data)

        write.csv(output_data, output_filename, row.names=FALSE, quote=FALSE)
    
        # Change dir to csv_dir so that 7z does not contain folder structure
        setwd(csv_dir)

        #compressed_file = gsub('csv', '7z', basename(output_filename))
        #compressed_command = paste0('7z a ', compressed_file, ' ', basename(output_filename))

        compressed_file = paste0(basename(output_filename), '.bz2')
        compressed_command = paste0('bzip2 -f ', basename(output_filename))
        system(compressed_command)
        setwd(mydir)
        
        # Return compressed file so we can operate on it later
        compressed_filename = paste0(dirname(output_filename), '/', compressed_file)
        return(compressed_filename)

    }



    # Run the creation command in parallel
    library(parallel)
    compressed_filenames = parallel::mclapply(as.list(input_RDS), output_compressed_bootstrap, mc.cores=12)
    compressed_filenames = unlist(compressed_filenames)

    # Move files to 'clean' directory
    final_output_dir = 'Storm_series_bootstrap_simulations'
    dir.create(final_output_dir, showWarnings=FALSE)
    target_files = compressed_filenames #gsub('RDS', '7z', paste0(csv_dir, '/', basename(input_RDS)) )
    target_files_new = paste0(final_output_dir, '/', target_files)
    dir.create(dirname(target_files_new[1]), showWarnings=FALSE)
    file.rename(target_files, target_files_new)

    # Save the synthetic series for the fitted model to the same folder
    # structure
    # Do it in an environment to keep everything clean
    #tmp = new.env()
    #load('Session_end.Rdata', envir=tmp)
    output_data = initial_fit_env$synthetic_attr

    # Convert years to hours
    output_data$duration = output_data$duration * initial_fit_env$year2hours
    #rm(tmp)
    output_data = format_output_data(output_data) 
    output_data_file = paste0(final_output_dir, '/Fitted_model_series.csv')
    write.csv(output_data, output_data_file,
        row.names=FALSE, quote=FALSE)
    setwd(final_output_dir)

    #seven_zip_command = paste0('7z a ', gsub('csv', '7z', basename(output_data_file)), ' ', basename(output_data_file))
    ## Bzip2 is only slightly worse, bit it's much faster 
    compressed_command = paste0('bzip2 -f ', ' ', basename(output_data_file))
    system(compressed_command)
    setwd(mydir)
    #file.remove(output_data_file)
    
}

#################################################################################
#
# Check parameter uncertainties
#
#################################################################################

pdf('methodology_parametric_bootstrap_checks_diagnostic_plots.pdf', width=12, height=10)

# Get event timing model fits

lnhp = length(initial_fit_env$best_nhp_model$par)
nhp_par = matrix( unlist(lapply(allRDS, f<-function(x) x$best_nhp_model$par)), 
    ncol=lnhp, byrow=TRUE)

par(mfrow=c(2,2))
for(i in 1:lnhp){
    # Visually check that initial fit is consistent with bootstrap
    hist(nhp_par[,i], n=20, main=paste0('NHP model parameter ', i))
    abline(v=initial_fit_env$best_nhp_model$par[i], col='red')
}

# Get convergence flag (should be zero)
nhp_conv = unlist(lapply(allRDS, f<-function(x) x$best_nhp_model$convergence))

range(nhp_conv)


# Get Hsig random parameters
hsig_random_par = matrix(unlist(lapply(allRDS, f<-function(x) x$hsig_random_par)), 
    ncol=length(initial_fit_env$hsig_mixture_fit$fit_optim$par), byrow=TRUE)

par(mfrow = c(2,2))
for(i in 1:4){
    # Visually check that distributions seem 'the-same'
    initial_fit_env$DU$qqplot3(
        initial_fit_env$hsig_mixture_fit$combined_chains[,i],
        hsig_random_par[,i])
    title(main=paste0('QQ plot hsig: parameter ', i))
    abline(0, 1, col='red')
    grid()
}


# Get Hsig bootstrap parameters
hsig_boot_par = matrix(unlist(lapply(allRDS, f<-function(x) x$bootstrap_hsig_mixture_fit)), 
    ncol=length(initial_fit_env$hsig_mixture_fit$fit_optim$par), byrow=TRUE)

par(mfrow=c(2,2))
for(i in 1:4){
    hist(hsig_boot_par[,i], main=paste0('Hsig bootstrap parameter ', i))
    abline(v=initial_fit_env$hsig_mixture_fit$fit_optim$par[i], col='red')
}

# Get duration random parameters
duration_random_par = matrix(unlist(lapply(allRDS, f<-function(x) x$duration_random_par)), 
    ncol=length(initial_fit_env$duration_mixture_fit$fit_optim$par), byrow=TRUE)

par(mfrow = c(2,2))
for(i in 1:4){
    # Visually check that distributions seem 'the-same'
    initial_fit_env$DU$qqplot3(
        initial_fit_env$duration_mixture_fit$combined_chains[,i],
        duration_random_par[,i])
    title(main=paste0('QQ plot duration: parameter ', i))
    abline(0, 1, col='red')
    grid()
}


# Get duration bootstrap parameters
duration_boot_par = matrix(unlist(lapply(allRDS, f<-function(x) x$bootstrap_duration_mixture_fit)), 
    ncol=length(initial_fit_env$duration_mixture_fit$fit_optim$par), byrow=TRUE)

par(mfrow=c(2,2))
for(i in 1:4){
    hist(duration_boot_par[,i], main=paste0('Duration bootstrap parameter ', i))
    abline(v=initial_fit_env$duration_mixture_fit$fit_optim$par[i], col='red')
}


# Get tideResid random parameters
tideResid_random_par = matrix(unlist(lapply(allRDS, f<-function(x) x$tideResid_random_par)), 
    ncol=length(initial_fit_env$tideResid_mixture_fit$fit_optim$par), byrow=TRUE)

par(mfrow = c(2,2))
for(i in 1:4){
    # Visually check that distributions seem 'the-same'
    initial_fit_env$DU$qqplot3(
        initial_fit_env$tideResid_mixture_fit$combined_chains[,i],
        tideResid_random_par[,i])
    title(main=paste0('QQ plot tideResid: parameter ', i))
    abline(0, 1, col='red')
    grid()
}


# Get tideResid bootstrap parameters
tideResid_boot_par = matrix(unlist(lapply(allRDS, f<-function(x) x$bootstrap_tideResid_mixture_fit)), 
    ncol=length(initial_fit_env$tideResid_mixture_fit$fit_optim$par), byrow=TRUE)

par(mfrow=c(2,2))
for(i in 1:4){
    hist(tideResid_boot_par[,i], main=paste0('tideResid bootstrap par ', i))
    abline(v=initial_fit_env$tideResid_mixture_fit$fit_optim$par[i], col='red')
}


# copula parameters
# Need to account for the fact that the family is varying in the
# bootstrap/data perturbation case
cop_par = array(
    unlist(lapply(allRDS, f<-function(x) x$copula_model$copula_fit_mle$RVM$par)), 
    dim = c(dim(allRDS[[1]]$copula_model$copula_fit_mle$RVM$par), length(allRDS)))
cop_fam = array(
    unlist(lapply(allRDS, f<-function(x) x$copula_model$copula_fit_mle$RVM$family)), 
    dim = c(dim(allRDS[[1]]$copula_model$copula_fit_mle$RVM$family), length(allRDS)))

# Only the lower triangular entries of cop_par and cop_fam are meaningful
# Here that means sum(1:4) plots
par(mfrow=c(4,3))
for(i in 1:nrow(cop_par[,,1])){
    for(j in 1:ncol(cop_par[,,1])){
        # Ignore main diagonal + upper-triangle of matrix
        if(i <= j) next

        # Note the independence copulas will plot as a single bin.
        fams = cop_fam[i,j,]
        unique_fams = unique(fams)
        vals = cop_par[i,j,]
        boxplot(vals ~ fams, col=1 + (1:length(unique_fams)))
        abline(h = initial_fit_env$copula_model$copula_fit_mle$RVM$par[i,j], 
            col=1 + match(initial_fit_env$copula_model$copula_fit_mle$RVM$family[i,j], unique_fams))
        title(paste0(i, ',', j))
    }
}


# Length of synthetic series
len_synth = unlist(lapply(allRDS, f<-function(x) length(x$synthetic_attr[,1])))
mean_rate = len_synth/series_duration_years

par(mfrow=c(1,1))
hist(mean_rate, main='Annual storm rate')
original_mean_rate = length(initial_fit_env$synthetic_attr[,1])/initial_fit_env$nyears_synthetic_series
abline(v=original_mean_rate, col='red')


# Direction
par(mfrow=c(1,1))
tmp = hist(initial_fit_env$dir_fit_conditional$q_raw(runif(1e+04)), n=50, 
    freq=FALSE, col='red', main='DIR initial fit and synthetic series')
for(i in 1:length(allRDS)){
    tmp2 = hist(allRDS[[i]]$dir_fit_conditional$q_raw(runif(1e+04)), n=50, 
        plot=FALSE)
    points(tmp2$mids, tmp2$density, t='l', col='grey')
}
points(tmp$mids, tmp$density, t='o', col='black', lwd=4)

# steepness
par(mfrow=c(1,1))
tmp = hist(initial_fit_env$steepness_fit_conditional$q_raw(runif(1e+04)), n=10, 
    freq=FALSE, col='red', main='Steepness initial fit and synthetic series')
for(i in 1:length(allRDS)){
    tmp2 = hist(allRDS[[i]]$steepness_fit_conditional$q_raw(runif(1e+04)), n=10, 
        plot=FALSE)
    points(tmp2$mids, tmp2$density, t='l', col='grey')
}
points(tmp$mids, tmp$density, t='o', col='black', lwd=4)


# Direction in el-nino

par(mfrow=c(1,1))
kk = which(initial_fit_env$synthetic_attr$soiA < -5)

elnino_dir = hist(initial_fit_env$synthetic_attr$dir[kk], n=50, 
    freq=FALSE, col='red', main='DIR El-Nino initial fit and synthetic series')

kk = which(initial_fit_env$synthetic_attr$soiA > 5)
lanina_dir = hist(initial_fit_env$synthetic_attr$dir[kk], n=50, 
    plot=FALSE)

for(i in 1:length(allRDS)){
    kk = which(allRDS[[i]]$synthetic_attr$soiA < -5)
    tmp2 = hist(allRDS[[i]]$synthetic_attr$dir[kk], n=50, 
        plot=FALSE)
    points(tmp2$mids, tmp2$density, t='l', col='grey')
}
points(elnino_dir$mids, elnino_dir$density, t='o', col='black', lwd=4)
points(lanina_dir$mids, lanina_dir$density, t='l', col='blue')

# Direction in la-nina

par(mfrow=c(1,1))
kk = which(initial_fit_env$synthetic_attr$soiA > 5)
lanina_dir = hist(initial_fit_env$synthetic_attr$dir[kk], n=50, 
    freq=FALSE, col='red', main='DIR La-Nina initial fit and synthetic series')
for(i in 1:length(allRDS)){
    kk = which(allRDS[[i]]$synthetic_attr$soiA > 5)
    tmp2 = hist(allRDS[[i]]$synthetic_attr$dir[kk], n=50, 
        plot=FALSE)
    points(tmp2$mids, tmp2$density, t='l', col='grey')
}
points(lanina_dir$mids, lanina_dir$density, t='o', col='black', lwd=4)
points(elnino_dir$mids, elnino_dir$density, t='l', col='blue')

# More Direction/ENSO comparison
elnino_120_rate = unlist(lapply(allRDS, f<-function(x){
    keep = which(x$synthetic_attr$soiA < -5)
    events = sum(x$synthetic_attr$dir[keep] < 120)
    years = length(unique(floor(x$synthetic_attr$startyear[keep])))
    return(events/years)
    }
))

quantile(elnino_120_rate, c(0.025, 0.5, 0.975))

elnino_170_rate = unlist(lapply(allRDS, f<-function(x){
    keep = which(x$synthetic_attr$soiA < -5)
    events = sum(x$synthetic_attr$dir[keep] > 170)
    years = length(unique(floor(x$synthetic_attr$startyear[keep])))
    return(events/years)
    }
))
quantile(elnino_170_rate, c(0.025, 0.5, 0.975))

lanina_120_rate = unlist(lapply(allRDS, f<-function(x){
    keep = which(x$synthetic_attr$soiA > 5)
    events = sum(x$synthetic_attr$dir[keep] < 120)
    years = length(unique(floor(x$synthetic_attr$startyear[keep])))
    return(events/years)
    }
))
quantile(lanina_120_rate, c(0.025, 0.5, 0.975))

lanina_170_rate = unlist(lapply(allRDS, f<-function(x){
    keep = which(x$synthetic_attr$soiA > 5)
    events = sum(x$synthetic_attr$dir[keep] > 170)
    years = length(unique(floor(x$synthetic_attr$startyear[keep])))
    return(events/years)
    }
))
quantile(lanina_170_rate, c(0.025, 0.5, 0.975))

# Seasonality in all parameters

season_par = c('hsig', 'duration', 'tp1', 'tideResid', 'dir')
for(sp in season_par){
    year_to_10 = ceiling(12 * (initial_fit_env$synthetic_attr$startyear - 
        floor(initial_fit_env$synthetic_attr$startyear)))
    fit1 = aggregate(initial_fit_env$synthetic_attr[[sp]], by=list(month=year_to_10), 
        median)

    # Deal with units for duration
    if(sp == 'duration'){
        scaler = initial_fit_env$year2hours
    }else{
        scaler = 1
    }

    plot(fit1[,1], fit1[,2]*scaler, t='o')

    for(i in 1:length(allRDS)){
        year_to_10 = ceiling(12 * (allRDS[[i]]$synthetic_attr$startyear - 
            floor(allRDS[[i]]$synthetic_attr$startyear)))
        tmp = aggregate(allRDS[[i]]$synthetic_attr[[sp]], by=list(month=year_to_10), 
            median)
        points(tmp[,1], tmp[,2], t='l', col='grey')
    }
    title(main=paste0('Monthly variation in median ', sp))
    points(fit1[,1], fit1[,2]*scaler, t='o', lwd=4)
}


###############################################################################
#
# Check that distribution of 'event-gaps',  is appropriately uniform once
# transformed with the theoretical model
#
###############################################################################
check_time_between_events<-function(fit_env, duration_in_years=FALSE){

    # Find the maximum gap and the time that it starts
    if(duration_in_years){
        scaler = 1
    }else{
        scaler = initial_fit_env$year2hours
    }
 
    start_of_gap = fit_env$synthetic_attr$startyear + 
        fit_env$synthetic_attr$duration/scaler + 
        initial_fit_env$duration_gap_hours/initial_fit_env$year2hours

    lx = length(start_of_gap)

    stopifnot(all(start_of_gap[1:(lx-1)] < fit_env$synthetic_attr$startyear[2:lx]))

    # Compute the probability of the gap

    soiAfun_years = rle(floor(fit_env$synthetic_attr$startyear))$values
    soiAfun_soiA = rle(fit_env$synthetic_attr$soiA)$values
    stopifnot(length(soiAfun_years) == length(soiAfun_soiA))
    soiA_fun_for_lambda = approxfun(soiAfun_years, soiAfun_soiA)

    source('../../R/nhpoisp/nhpoisp.R', local=TRUE)

    lambda_fun = get_lambda_function(
        theta = fit_env$best_nhp_model$par, 
        rate_equation = fit_env$best_nhp_model$rate_equation,
        minimum_rate = 0)

    gap_prob = rep(NA, (lx-1))
    for(i in 1:(lx-1)){

        #gap_prob[i] = integrate(lambda_fun, start_of_gap[i], fit_env$synthetic_attr$startyear[i+1], 
        #    rel.tol=1.0e-12, abs.tol=1.0e-12)$value
        gap_prob[i] = integrate(lambda_fun, start_of_gap[i], fit_env$synthetic_attr$startyear[i+1], 
            rel.tol=1.0e-6, abs.tol=1.0e-6)$value
    }
   
    # This should be uniformly distributed 
    gap_prob = (1 - exp(-gap_prob))
    
    ks_test = ks.test(gap_prob, 'punif')

    return(list(gap_prob=gap_prob, ks_test = ks_test, max_gap_prob = max(gap_prob), 
        min_gap_prob = min(gap_prob), theoretical_max = lx/(lx+1)))
}

# Check the distribution of the gap F value,  vs the theoretical uniform distribution
event_gap_summary = parallel::mclapply(allRDS, f<-function(x) check_time_between_events(x), mc.cores=12)

# Distribution of the 'most extreme' gap
empirical_dist = unlist(lapply(event_gap_summary, f<-function(x) x$max_gap_prob))
theoretical_dist = replicate(10000, {max(runif(original_mean_rate * 1000))})
initial_fit_env$DU$qqplot3(theoretical_dist, empirical_dist, 
    main='QQ-plot of maximum event gap nonexceedence probability vs theoretical value')
abline(0, 1, col='red')

# Distribution of ks.test results
empirical_dist = unlist(lapply(event_gap_summary, f<-function(x) x$ks_test$p.value))
hist(empirical_dist, main='p-values of ks.test results (should be uniform)')

# Distribution of the 'least extreme' gap
empirical_dist = unlist(lapply(event_gap_summary, f<-function(x) x$min_gap_prob))
theoretical_dist = replicate(10000, {min(runif(original_mean_rate * 1000))})
initial_fit_env$DU$qqplot3(theoretical_dist, empirical_dist, 
    main='QQ-plot of minimum event gap nonexceedence probability vs theoretical value')
abline(0, 1, col='red')

###############################################################################
#
# Hundred-year hsig, duration, tideResid
#
###############################################################################

# To get percentiles corresponding to the desired ARI, need to adjust for
# the mean-rate != unity
ari_100_percentile = 1 - 1/(mean_rate*100)

# Try another approach based on pure table lookup
hsig_100b = unlist(mapply(
    f<-function(x, ari_100_percentile){
        quantile(x$synthetic_attr$hsig, p=ari_100_percentile)
    },
    x=allRDS, ari_100_percentile = ari_100_percentile))

hist(initial_fit_env$hsig_mixture_fit$combined_ari100, freq=FALSE, n=20,
    main='ARI 100 hsig in initial bayesian fit (black), and empirically from synthetic series')
hist(hsig_100b, freq=FALSE, add=T, col='red', density=20)


# Do comparison with duration
duration_100b = unlist(mapply(
    f<-function(x, ari_100_percentile){
        quantile(x$synthetic_attr$duration, p=ari_100_percentile)
    },
    x=allRDS, ari_100_percentile = ari_100_percentile))

hist(initial_fit_env$duration_mixture_fit$combined_ari100, freq=FALSE, n=20,
    main='ARI 100 duration in initial bayesian fit (black), and empirically from synthetic series')
hist(duration_100b, freq=FALSE, add=T, col='red', density=20)

# Do comparison with tideResid
tideResid_100b = unlist(mapply(
    f<-function(x, ari_100_percentile){
        quantile(x$synthetic_attr$tideResid, p=ari_100_percentile)
    },
    x=allRDS, ari_100_percentile = ari_100_percentile))

hist(initial_fit_env$tideResid_mixture_fit$combined_ari100, freq=FALSE, n=20,
    main='ARI 100 tideResid in initial bayesian fit (black), and empirically from synthetic series')
hist(tideResid_100b, freq=FALSE, add=T, col='red', density=20)


#################################################################################
#
# Maxima in synthetic series and relation to shape parameter
#
#################################################################################
par(mfrow=c(2,2))
hsig_max_synthetic_series = unlist(lapply(allRDS, f<-function(x) max(x$synthetic_attr$hsig)))
plot(hsig_random_par[,4], hsig_max_synthetic_series, 
    main='hsig: Synthetic series peak vs GPD shape parameter',
    xlab='GPD Shape parameter', 
    ylab='Max hsig in 1000 year synthetic series', log='y')

duration_max_synthetic_series = unlist(lapply(allRDS, f<-function(x) max(x$synthetic_attr$duration)))
plot(duration_random_par[,4], duration_max_synthetic_series, 
    main='duration: Synthetic series peak vs GPD shape parameter',
    xlab='GPD Shape parameter', ylab='Max duration in 1000 year synthetic series', log='y')

tideResid_max_synthetic_series = unlist(lapply(allRDS, f<-function(x) max(x$synthetic_attr$tideResid)))
plot(tideResid_random_par[,4], tideResid_max_synthetic_series, 
    main='tideResid: Synthetic series peak vs GPD shape parameter',
    xlab='GPD Shape parameter', ylab='Max tideResid in 1000 year synthetic series', log='y')

# How do bootstrap values compare?
# Get APPROXIMATE maximum hsig that would occur if we used a bootstrap instead.
par(mfrow=c(2,3))
hist(hsig_max_synthetic_series, n=50)
abline(v=max(initial_fit_env$event_statistics$hsig), col='red')
hist(duration_max_synthetic_series, n=50)
abline(v=max(initial_fit_env$event_statistics$duration), col='red')
hist(tideResid_max_synthetic_series, n=50)
abline(v=max(initial_fit_env$event_statistics$tideResid, na.rm=TRUE), col='red')

prob_val = 1 - 1/(1000*original_mean_rate)
hsig_max_using_bootstrap = initial_fit_env$hsig_mixture_fit$qf(prob_val, 
    hsig_boot_par) + initial_fit_env$hsig_threshold
hist(hsig_max_using_bootstrap, n=50)
abline(v=max(initial_fit_env$event_statistics$hsig), col='red')

duration_max_using_bootstrap = initial_fit_env$hsig_mixture_fit$qf(prob_val,
    duration_boot_par)
hist(duration_max_using_bootstrap, n=50)
abline(v=max(initial_fit_env$event_statistics$duration), col='red')

tideResid_max_using_bootstrap = initial_fit_env$hsig_mixture_fit$qf(prob_val, 
    tideResid_boot_par)
hist(tideResid_max_using_bootstrap, n=50)
abline(v=max(initial_fit_env$event_statistics$tideResid, na.rm=TRUE), col='red')


## Plot of bootstrap vs Bayes uncertainties
## Beware -- here 'aep' is not really AEP, it is a rate
## AEP = 1 - exp(-aep)

aeps = c(20, 15, 10, 7, 5, 3, 2, 1, 1/2, 1/5, 1/7, 1/10, 1/15, 1/20, 1/40, 1/50, 1/75, 1/100, 1/200, 1/500, 1/1000)
hsig_data = initial_fit_env$event_statistics$hsig
empirical_aep = (1 - (rank(hsig_data)/(length(hsig_data)+1)))*original_mean_rate

par(mfrow=c(1,1))
plot(empirical_aep, hsig_data, log='xy', xlim=c(20, 0.01), ylim=c(2.8, 10), pch=19, cex=1)

bayes_aep = matrix(NA, ncol=3, nrow = length(aeps))
boot_aep = matrix(NA, ncol=3, nrow = length(aeps))
for(i in 1:length(aeps)){
    print(i)
    aep = aeps[i]
    prob_val = (1 - aep/original_mean_rate)
    hsig_max_using_bootstrap = initial_fit_env$hsig_mixture_fit$qf(prob_val, 
        hsig_boot_par) + initial_fit_env$hsig_threshold
    hsig_max_using_bayes = initial_fit_env$hsig_mixture_fit$qf(prob_val, 
        hsig_random_par) + initial_fit_env$hsig_threshold

    bayes_aep[i,] = quantile(hsig_max_using_bayes, c(0.025, 0.5, 0.975))
    boot_aep[i,] = quantile(hsig_max_using_bootstrap, c(0.025, 0.5, 0.975))
}

for(j in c(1, 3)) points(aeps, bayes_aep[,j], t='l', col='red')
points(aeps, bayes_aep[,2], t='l', col='red', lwd=4)
for(j in 1:3) points(aeps, boot_aep[,j], t='l', col='blue')
points(aeps, boot_aep[,2], t='l', col='blue', lwd=4)
grid()
ml_hsig = initial_fit_env$hsig_mixture_fit$qfun(1-aeps/original_mean_rate)
points(aeps, ml_hsig, t='l', 
    lwd=4, col='black', lty='dashed')
legend('topleft', 
    c('Bayesian (0.025, 0.5, 0.975)', 'Parametric Bootstrap (0.025, 0.5, 0.975)', 'Maximum Likelihood Estimator'),
    col=c('red', 'blue', 'black'), lty=c('solid', 'solid', 'dashed'))
title('Comparison Hsig RL Bayesian and Bootstrap')

output_data = data.frame(aep=aeps, bayes_lower = bayes_aep[,1], bayes_med = bayes_aep[,2], bayes_upper = bayes_aep[,3],
    boot_lower = boot_aep[,1], boot_med = boot_aep[,2], boot_upper = boot_aep[,3], ml_hsig = ml_hsig)

write.table(output_data, 'hsig_return_periods_detailed.csv', sep=",", row.names=FALSE)
# Finish diagnostic plots
dev.off()
stop('Deliberate stop here')

###################################################################################
#
# Better plot like above, for hsig, duration, and tidal residual
#
#

# Make a nice plot
png('hsig_return_period.png', width=8, height=6, units='in', res=200)

true_aep = 1 - exp(-output_data$aep)
plot(true_aep, output_data$ml_hsig, t='l', lwd=2, xlim=c(1, 1/1000), 
    ylim=c(2.9, 10.5), log='x', frame.plot=TRUE, axes=FALSE, 
    xlab='Annual Exceedance Probability', ylab="")
axis(side=2)
mtext(side=2, bquote(paste(H[sig], ' (m)')), line=2, cex=1.5)
axis(side=1, at = c(1, 1/10, 1/100, 1/1000), labels=c('1', '1/10', '1/100', '1/1000'))
grid(col='brown')

points(true_aep,  output_data$bayes_lower, t='l', col='red', lwd=2)
points(true_aep,  output_data$bayes_upper, t='l', col='red', lwd=2)
#points(true_aep,  output_data$bayes_med, t='l', col='red')

points(true_aep, output_data$boot_lower, t='l', col='green', lty='dashed', lwd=2)
points(true_aep, output_data$boot_upper, t='l', col='green', lty='dashed', lwd=2)

points(1 - exp(-empirical_aep), hsig_data, pch=19, cex=0.6)

arrows(x0 = 1/100, x1 = 1/100, y0=7.4, y1=10.4, angle=90, col='brown', length=0.1)
arrows(x0 = 1/100, x1 = 1/100, y0=7.4, y1=6.97, angle=90, col='brown', length=0.1)
points(1/100, 7.4, pch=8, col='brown', cex=2)


legend('bottomright', c('Data', 'Maximum likelihood estimator', 'Bayesian 95% interval', 'Bootstrap 95% interval', 'GEV 1/100 95% Interval'),
    lty=c(NA, 'solid', 'solid', 'dashed', 'solid'), col=c('black', 'black', 'red', 'green', 'brown'), lwd=2, bg='white', pch=c(19, NA, NA, NA, 8))
dev.off()


################################################################################
#
# Recurrence of 1974 storm
#
###############################################################################
source('../statistical_model_fit/cluster_counting_functions.R')
library(parallel)

# Example definitions 'like the 1974 storm'
# The event_cluster_counter is currently slow, but could be made very fast in Rcpp
rates_5_3 = mcmapply( 
    f<-function(x){
        event_cluster_counter(x$synthetic_attr, hsig_threshold=5, nevents=3, 
            time_window=6*7/365.25)
    },
    x=allRDS,
    mc.cores=6, SIMPLIFY=TRUE, USE.NAMES=FALSE
    )
rates_5_3_max6 = mcmapply( 
    f<-function(x){
        event_cluster_counter(x$synthetic_attr, hsig_threshold=5, nevents=3, 
            time_window=6*7/365.25, max_hsig_threshold=6)
    },
    x=allRDS,
    mc.cores=6, SIMPLIFY=TRUE, USE.NAMES=FALSE
    )

## Get confidence interval for 1000 year series
quantile(rates_5_3/1000, p = c(0.025, 0.5, 0.975))
quantile(rates_5_3_max6/1000, p = c(0.025, 0.5, 0.975))


# Compute some statistics for all bootstrap events
rates_55_2 = mcmapply(
    f<-function(x){event_cluster_counter(x$synthetic_attr, 
        hsig_threshold=5.5, nevents=2, time_window=6*7/365.25)},
    x = allRDS, mc.cores=6, SIMPLIFY=TRUE, USE.NAMES=FALSE)

rates_5_2 = mcmapply(
    f<-function(x){event_cluster_counter(x$synthetic_attr, 
        hsig_threshold=5., nevents=2, time_window=6*7/365.25)},
    x = allRDS, mc.cores=6, SIMPLIFY=TRUE, USE.NAMES=FALSE)

rates_45_2 = mcmapply(
    f<-function(x){event_cluster_counter(x$synthetic_attr, 
        hsig_threshold=4.5, nevents=2, time_window=6*7/365.25)},
    x = allRDS, mc.cores=6, SIMPLIFY=TRUE, USE.NAMES=FALSE)

rates_4_2 = mcmapply(
    f<-function(x){event_cluster_counter(x$synthetic_attr, 
        hsig_threshold=4., nevents=2, time_window=6*7/365.25)},
    x = allRDS, mc.cores=6, SIMPLIFY=TRUE, USE.NAMES=FALSE)




# Compare events above with empirical estimates from the data, assuming the
# cluster events in the data can be approximated as poisson in time
obs = initial_fit_env$event_statistics #read.csv('event_statistics_out.csv')

obs_5_3 = event_cluster_counter(obs, hsig_threshold=5., nevents=3, time_window=6*7/365.25)
poisson.test(obs_5_3, T=(2015-1985))
quantile(rates_5_3/1000, p=c(0.025, 0.5, 0.975))

obs_55_2 = event_cluster_counter(obs, hsig_threshold=5.5, nevents=2, time_window=6*7/365.25)
poisson.test(obs_55_2, T=(2015-1985))
quantile(rates_55_2/1000, p=c(0.025, 0.5, 0.975))

obs_5_2 = event_cluster_counter(obs, hsig_threshold=5., nevents=2, time_window=6*7/365.25)
poisson.test(obs_5_2, T=(2015-1985))
quantile(rates_5_2/1000, p=c(0.025, 0.5, 0.975))

obs_45_2 = event_cluster_counter(obs, hsig_threshold=4.5, nevents=2, time_window=6*7/365.25)
poisson.test(obs_45_2, T=(2015-1985))
quantile(rates_45_2/1000, p=c(0.025, 0.5, 0.975))

obs_4_2 = event_cluster_counter(obs, hsig_threshold=4., nevents=2, time_window=6*7/365.25)
poisson.test(obs_4_2, T=(2015-1985))
quantile(rates_4_2/1000, p=c(0.025, 0.5, 0.975))

# Compare the model/data predictions for event sequences
#' Function to do the work
model_data_cluster_compare = function(hsig_threshold, nevents, time_window){

    model_rates = mcmapply(
        f<-function(x){
            event_cluster_counter(x$synthetic_attr, 
                hsig_threshold=hsig_threshold, 
                nevents=nevents, time_window=time_window)},
            x = allRDS, mc.cores=6, SIMPLIFY=TRUE, USE.NAMES=FALSE)

    data_rates = event_cluster_counter(
        initial_fit_env$event_statistics, 
        hsig_threshold=hsig_threshold,
        nevents = nevents,
        time_window = time_window)

    return(list(model_rates = model_rates, data_rates = data_rates))

}

## Rate of 2 events in 4 weeks with hsig > threshold
hsig_threshold_seq = seq(4.0, 7.0, by = 0.25)
model_store = list()
for(i in 1:length(hsig_threshold_seq)){
    model_store[[i]] = model_data_cluster_compare(hsig_threshold_seq[i], 2, 4*7/365.25)
    names(model_store)[i] = as.character(hsig_threshold_seq[i])
}

tmp = rep(NA, length(hsig_threshold_seq))
output = data.frame(obs_rate = tmp, model_rate = tmp, model_lower=tmp, model_upper=tmp, hsig_threshold = tmp)
for(i in 1:length(model_store)){
    output$hsig_threshold[i] = hsig_threshold_seq[i]
    tmp = quantile(model_store[[i]]$model_rates, c(0.025, 0.5, 0.975))/1000
    output$model_lower[i] = tmp[1]
    output$model_rate[i] = tmp[2]
    output$model_upper[i] = tmp[3]
    output$obs_rate[i] = model_store[[i]]$data_rates/initial_fit_env$data_duration_years
}

pdf('clustered_event_rate_test.pdf', width=8, height=6)
plot(output$hsig_threshold, output$model_rate, ylim=c(0.001, 2), t='l',
    xlab='Clustered event definition threshold wave height', ylab='Clustered event rate (#/year)',
    log='y')
points(output$hsig_threshold, output$model_lower + 1.0e-012, t='l', col='red')
points(output$hsig_threshold, output$model_upper, t='l', col='red')
points(output$hsig_threshold, output$obs_rate + 1.0e-12, pch=19, col='blue', cex=2)
grid(col='brown')
legend('topright', c('Model', 'Model 95% Confidence Interval', 'Empirical Rate'),
    lty=c(1, 1, NA), col=c('black', 'red', 'blue'), pch=c(NA, NA, 19), bg='white',
    cex=1.)
dev.off()



#
# 
# Plot of seasonal enso impacts, using perturbed data fits
#
# This should be run in the statistical_model_fit_perturbed_data folder
#
stopifnot(basename(getwd()) == 'statistical_model_fit_perturbed_data')

# Make synthetic_attr by combining series from all the perturbed data fits [but
# not bootstrap fits]
session_series_simulation = Sys.glob('Rimages/session_series_simulation_TRUE_*.Rdata')
synthetic_series_list = vector(mode='list', length=length(session_series_simulation))
for(i in 1:length(session_series_simulation)){
    print(i)
    tmp_env = new.env()
    load(session_series_simulation[i], envir=tmp_env)
    synthetic_series_list[[i]] = tmp_env$synthetic_attr
}

synthetic_attr = synthetic_series_list[[1]] * NA

synthetic_attr = do.call(rbind, synthetic_series_list)


month_days = c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
month_breaks = (c(0, cumsum(month_days)))/sum(month_days)

plotpar = c('hsig', 'duration', 'tideResid', 'dir', 'tp1')
plotpar_title = c('Hsig', 'D', 'R', expression(theta), 'T')
plot_ylabs = c('m', 'hours', 'm', 'degrees', 's')

dec<-function(x) x - floor(x)
month_cat = findInterval(dec(synthetic_attr$startyear), month_breaks)

png('modelled_seasonal_enso_impacts.png', width=7, height=7, units='in', res=200)
par(mfrow=c(3,2))
par(mar=c(4,4,2,1))

dir1 = hist(synthetic_attr$dir, plot=FALSE, breaks=seq(55, 205, by=5))
plot(dir1$mids, dir1$density, t='l', ylim=c(0.001, 0.035), lwd=2, col='red',
    xlab='Direction', ylab='Density (note log scale)', 
    main = 'Storm wave direction density ', log='y', las=1)    
kk = which(synthetic_attr$soiA < -5)
dir1 = hist(synthetic_attr$dir[kk], plot=FALSE, breaks=seq(55, 205, by=5))
mean(synthetic_attr$dir[kk] < 120)
points(dir1$mids, dir1$density, t='l', col='orange', lty='dashed', lwd=2)
kk = which(synthetic_attr$soiA > 5)
dir1 = hist(synthetic_attr$dir[kk], plot=FALSE, breaks=seq(55, 205, by=5))
mean(synthetic_attr$dir[kk] < 120)
points(dir1$mids, dir1$density, t='l', col='blue', lty='dashed', lwd=2)
grid()
legend('topleft', c('Model', 'La-Nina', 'El-Nino'), 
    lty=c('solid', 'dashed', 'dashed'), lwd=2, col=c('red', 'blue', 'orange'))


for(i in 1:5){
    if(plotpar[i] == 'duration'){
        boxplot((synthetic_attr[[plotpar[i]]][1:1e+05]*365.25*24) ~ month_cat[1:1e+05],
            main=plotpar_title[i], names=month.abb, col='grey', las=1,
            ylab=plot_ylabs[i])
    }else if(plotpar[i] == 'tideResid'){
        boxplot((synthetic_attr$tideResid + synthetic_attr$msl)[1:1e+05] ~ month_cat[1:1e+05],
            main='Tidal residual (including MSL)', names=month.abb, col='grey', las=1,
            ylab=plot_ylabs[i])
    }else{
        boxplot(synthetic_attr[[plotpar[i]]][1:1e+05] ~ month_cat[1:1e+05],
            main=plotpar_title[i], names=month.abb, col='grey', las=1,
            ylab=plot_ylabs[i])
    }
    grid()
}
dev.off()

#
#
# Compare model [over all data perturbations] with data
#
#

nice_pairs_paper<-function( mydata, labels, extra_data=NULL, mydata_col='blue', add_signif_stars=TRUE, copula_type=NULL){
    # Prettier pairs plot
    pairs(mydata, labels,
        pch='.', cex=3, 
        upper.panel=function(x,y,...){ 
            points(x,y, col=mydata_col, ...); 
            grid(col='orange')
            if(!is.null(extra_data)){
                # Hack to identify the data columns
                icol = get('i', pos=parent.frame(n=2))
                jcol = get('j', pos=parent.frame(n=2))
                points(extra_data[,jcol], extra_data[,icol], col='black', pch='.', cex=2)
            }

        }, 
        diag.panel=DU$panel.hist,
        lower.panel=function(x,y,...){ 
            points(x,y, col=0); 

            spearman_cortest = try(cor.test(x, y, method='s'))
            if(class(spearman_cortest) != 'try-error'){
                spearman_cor = spearman_cortest$estimate
                title_word = as.character(round(spearman_cor, 3))
                font_size = 2.0

                if(add_signif_stars){
                    if(spearman_cortest$p.value < 0.05){
                        title_word = paste0(title_word, '*')
                        font_size = font_size + 0.5
                    }
                    if(spearman_cortest$p.value < 0.005){
                        title_word = paste0(title_word, '*')
                        font_size = font_size + 0.5
                    }
                }else{
                    font_size = font_size + 0.5 
                }

                if(!is.null(extra_data)){
                    # Hack to identify the data columns
                    icol = get('i', pos=parent.frame(n=2))
                    jcol = get('j', pos=parent.frame(n=2))
                    #points(extra_data[,jcol], extra_data[,icol], col='black', pch='.', cex=2)
                    extra_cortest = try(cor.test(extra_data[,jcol], extra_data[,icol], method='s'))
                    
                    title_word = paste0(title_word, '\n(', as.character(round(extra_cortest$estimate,3)), 
                        ') \n')#, copula_type[icol, jcol])
                }

                title(main=title_word, line=-8, col.main='black', cex.main=font_size) 

                if(!is.null(extra_data)){
                    title(copula_type[icol, jcol], line=-12.3, col.main='black', cex.main=font_size/1.5)
                }
            }
        })
}



event_statistics = initial_fit_env$event_statistics
plotvar = c('hsig', 'duration', 'tideResid', 'tp1', 'dir') #c('duration', 'hsig', 'tp1', 'dir', 'tideResid')
plotvar_names = c('Hsig (m)', 'D (hours)', 'R (m)', 'T (s)', expression(theta ~ (degrees) )) #c('D', 'Hsig', 'T', 'theta', 'h')
event_stats_tmp = event_statistics[,plotvar]
names(event_stats_tmp) = plotvar_names

png('data_pairs.png', width=10, height=10, units='in', res=200)
nice_pairs_paper(event_stats_tmp, labels=plotvar_names)
dev.off()


## Make matrix recording copula types for plot
vc_summary_stats = readRDS('vine_copula_runs_summary_statistics.RDS')
stopifnot(basename(getwd()) == 'statistical_model_fit_perturbed_data')
#copula_type = matrix(
#    VineCopula::BiCopName(initial_fit_env$copula_model$copula_fit_mle$RVM$family, short=FALSE), 
#    ncol=5, nrow=5)
##copula_type = matrix("", ncol=5, nrow=5) # We are pooling data, so should not do this
#copula_type = t(copula_type[5:1,5:1])
#copula_type[4,] = paste(copula_type[4,], '*', sep="")
#copula_type[5,4] = paste(copula_type[5,4], '*', sep="")
#copula_type[upper.tri(copula_type, diag=TRUE)] = ""
#copula_type = gsub('Independence', 'Indep', copula_type)
#copula_type = gsub('Survival Gumbel', '180_Gumbel', copula_type)
#
copula_type = matrix("", ncol=5, nrow=5)
for(i in 1:4){
    for(j in (i+1):5){
    
        all_cops = unlist(
            lapply(vc_summary_stats, f<-function(x) x$copula_model$copula_fit_mle$RVM$family[6-i,6-j])
            )
        uac = unique(all_cops)
        for(k in 1:length(uac)){
            nm = BiCopName(uac[k], short=FALSE)
            if(i == 4 | j == 4) nm = paste0(nm, '*')
            nm_frac = paste0(sum(uac[k] == all_cops), '%')
            copula_type[j,i] = paste0(copula_type[j,i], nm, ' ', nm_frac, ' \n ')
        }
    }
}

copula_type = gsub('Independence', 'Indep', copula_type)
copula_type = gsub('Survival Gumbel', '180_Gumbel', copula_type) 


png('model_data_pairs.png', width=10, height=10, res=200, units='in')
ri = sample(1:length(synthetic_attr[,1]), size=10000, replace=FALSE)
synthetic_attr_tmp = synthetic_attr[ri, plotvar]
synthetic_attr_tmp$duration = synthetic_attr_tmp$duration * 365.25 * 24
names(synthetic_attr_tmp) = plotvar_names
nice_pairs_paper(synthetic_attr_tmp, labels=plotvar_names, #names(synthetic_attr_tmp), 
    extra_data = event_stats_tmp, mydata_col='green', add_signif_stars=FALSE,
    copula_type = copula_type)
dev.off()


#
# Copula homogeneity test
#

set.seed(1)
s2 = (which(!is.na(event_statistics$dir)& !is.na(event_statistics$tideResid)))
s1 = sample(1:length(synthetic_attr[,1]), size=length(s2))
c_vine_node_order = initial_fit_env$c_vine_node_order
m1 = as.matrix(synthetic_attr[s1, c_vine_node_order])
m1[,'duration'] = m1[,'duration'] * 365.25 * 24 # duration in hours
m2 = as.matrix(event_statistics[s2, c_vine_node_order])
m2 = jitter(m2, amount=1e-05) # Break ties in data

TwoCop(m1, m2)

## > TwoCop(m1, m2)
## $pvalue
## [1] 0.22
## 
## $cvm
## [1] 0.02980479
## 
## $VaR
##        95% 
## 0.04793045 
## 
## $cvmsim
##   [1] 0.02711472 0.03046887 0.02111260 0.01237741 0.05125965 0.01681944
##   [7] 0.02147249 0.02220483 0.03060476 0.02668288 0.02480595 0.01887465
##  [13] 0.02249961 0.01606233 0.03939016 0.06398929 0.02046512 0.02209881
##  [19] 0.02881696 0.02394589 0.02194467 0.02764813 0.03255751 0.02113861
##  [25] 0.02974173 0.02500367 0.02087241 0.03553133 0.03317558 0.01735588
##  [31] 0.02981508 0.02606398 0.01552303 0.02224072 0.03622861 0.01526160
##  [37] 0.03425128 0.04808110 0.02191122 0.02675620 0.01867642 0.02051027
##  [43] 0.01809407 0.02771165 0.04030497 0.01368774 0.03317211 0.02098451
##  [49] 0.05000661 0.06027332 0.01171297 0.02214880 0.02255617 0.02059420
##  [55] 0.02159111 0.02407127 0.01630921 0.02260625 0.02787017 0.02246802
##  [61] 0.01865251 0.02962439 0.03906682 0.02379714 0.02454988 0.02068883
##  [67] 0.02015369 0.02196137 0.04777760 0.01989220 0.02476762 0.02750582
##  [73] 0.04570753 0.02451394 0.02705929 0.02805862 0.01923541 0.01598564
##  [79] 0.02784126 0.02275817 0.02050807 0.04509234 0.02660505 0.02770208
##  [85] 0.04792252 0.01942857 0.01943109 0.01468725 0.01582234 0.01380013
##  [91] 0.02487943 0.02647437 0.01602166 0.04784928 0.02665370 0.02939497
##  [97] 0.01398869 0.02324431 0.01194160 0.02069712

png('Empirical_copula_plot.png', width=12, height=10, units='in', res=200)
nn = length(m1[,1]) + 1
contour_levels = seq(0.05, 1.8, by=0.05)
contour_col = rev(rainbow(length(contour_levels)+3)[1:length(contour_levels)])
par(oma=c(0,0,0,8))
par(mfrow=c(5, 5))
par(mar=c(2,2,2,1))
for(i in 1:5){
    for(j in 1:5){
        if(i > j){
    
            #BiCopKDE(rank(m1[,i])/nn, rank(m1[,j])/nn, margins='unif', kde.pars=list(method='T', renorm.iter=10),
            #    levels=contour_levels, col=contour_col, lwd=2)
    
            v1 = bkde2D(cbind(rank(m1[,i])/nn, rank(m1[,j])/nn), bandwidth=c(0.1,0.1))
            contour(seq(0,1,len=51), seq(0,1,len=51), v1$fhat, levels=contour_levels, col=contour_col)
            points(rank(m1[,i])/nn, rank(m1[,j])/nn, col='black', pch=19, cex=0.5)
            title('Model', line=0.3, cex.main=1.8)
            #BiCopKDE(rank(m1[,i])/nn, rank(m1[,j])/nn, col='red', margins='unif')
            #BiCopKDE(rank(m2[,i])/nn, rank(m2[,j])/nn, add=TRUE, margins='unif')
        }
        if(j > i){
            v1 = bkde2D(cbind(rank(m2[,i])/nn, rank(m2[,j])/nn), bandwidth=c(0.1,0.1))
            contour(seq(0,1,len=51), seq(0,1,len=51), v1$fhat, levels=contour_levels, col=contour_col)
            points(rank(m2[,i])/nn, rank(m2[,j])/nn, col='black', pch=19, cex=0.5)
            title('Data', line=0.3, cex.main=1.8)
            #plot(rank(m1[,i])/nn, rank(m1[,j])/nn, col='red', pch=19, cex=0.5)
            #points(rank(m2[,i])/nn, rank(m2[,j])/nn, col='black', pch=19, cex=0.5)
        }

        if(i == j){
            if(i != 4){
                plot(c(0,1), c(0,1), col='white', axes=FALSE, frame.plot=TRUE, 
                    xlab="", ylab="", main = plotvar_names[i], cex.main=2.5, line=-7)
            }else{

                plot(c(0,1), c(0,1), col='white', axes=FALSE, frame.plot=TRUE, 
                    xlab="", ylab="", main = expression('S (m/m)'), cex.main=2.5, line=-7)
            }
        }
    }
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

rec_left = 0.95
rec_right = 1.0 
rect(rep(rec_left , len=length(contour_levels)), (contour_levels-0.025) - 1, 
     rep(rec_right, len=length(contour_levels)), (contour_levels+0.025) - 1, 
    col=contour_col)
text(0.5*(rec_left+rec_right), max(contour_levels)+0.025 - 1, 
    'Empirical \n Copula Density', pos=3, cex=1.5)
lab_levels = c(0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8) 
text(1.00, lab_levels-1, labels=lab_levels,
    pos=4, cex=1.5)
dev.off()


#smoothscatter(rank(m1[,i])/nn, rank(m1[,j])/nn),
#    levels=contour_levels, col=contour_col, lwd=2)
#
#v1 = bkde2D(cbind(rank(m1[,i])/nn, rank(m1[,j])/nn), bandwidth=c(0.1, 0.1))


