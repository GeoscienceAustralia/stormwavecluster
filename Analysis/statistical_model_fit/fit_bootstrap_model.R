# Code to fit the statistical model to series simulated from the fitted model.
# This is used for uncertainty quantification.


#############################
#
# LOAD THE PREVIOUS SESSION
#
#############################

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

    MC_CORES = 1 # No parallel for tie-breaking cases at this level -- instead run many cases in parallel

}else{

    break_ties_with_jitter = FALSE
    session_n = 0

    # Definitions controlling the number of years in the synthetic series
    nyears_synthetic_series = 1e+03 

    MC_CORES = parallel::detectCores()

}

# Useful constant
year2hours = 365.25 * 24

# Make a 'title' which can appear in filenames to identify this run
run_title_id = paste0(break_ties_with_jitter, '_', session_n)

# Put the previous R session in its own environment, rather than loading it directly
previous_R_session_file = paste0('Rimages/session_series_simulation_', run_title_id, '.Rdata')
initial_fit_env = new.env()
load(previous_R_session_file, envir=initial_fit_env)

#
# Read some key variables from the derived data
#

if(break_ties_with_jitter){
    derived_data_file = paste0('../preprocessing_perturbed_data/Derived_data/event_statistics_', run_title_id, '.RDS')
}else{
    derived_data_file = paste0('../preprocessing/Derived_data/event_statistics_', run_title_id, '.RDS')
}
event_statistics_list = readRDS(derived_data_file)

# Extract variables that we need from previous analysis
for(varname in names(event_statistics_list)){
    assign(varname, event_statistics_list[[varname]])
}

# Replace event_statistics with a synthetic version, generated from the fitted
# model
rm(event_statistics_list)
rm(event_statistics)


# Make an output directory for the bootstrap info, if it doesn't already exist
bootstrap_output_dir = 'Synthetic_series_bootstrap'
dir.create(bootstrap_output_dir, showWarnings=FALSE)

#
# Make a function which does the work so its easy to run in parallel
#
single_bootstrap<-function(event_statistics, iterate_number){

    # Need to re-load packages, as R does not automatically do this when re-loading
    # a session
    library(evmix)
    library(logspline)
    library(CDVine) # Used to set structure of C-Vine copula
    library(VineCopula) # Main copula fitting routine. 
    library(Matching)

    nhp = new.env()
    source('../../R/nhpoisp/nhpoisp.R', local=nhp)

    DU = new.env()
    source('../preprocessing/data_utilities.R', local=DU)

    evmix_fit = new.env()
    source('../../R/evmix_fit/evmix_fit.R', local=evmix_fit, chdir=TRUE)

    wavedisp = new.env()
    source('../../R/wave_dispersion/wave_dispersion_relation.R', local=wavedisp, chdir=TRUE)

    #################################################################################
    #
    # Fit the event timing model (FIXME: copied from other code)
    #
    #################################################################################

    # Get event start time in years
    event_time = event_statistics$startyear

    # Get event duration in years 
    event_duration_years = (event_statistics$duration + duration_gap_hours - 
        duration_offset_hours)/(year2hours)
    # NOTE: Because we demand a 'gap' between events, it makes sense to add this
    #       time to event duration, since no other event is possible in this time.

    # Get the observation start time in years
    obs_start_time = DU$time_to_year(obs_start_time_strptime)

    # Define the model
    # We need to provide soiA from the simulated series and fix up the soi function
    # 
    soiA_fun_unique_years = rle(floor(event_statistics$startyear))$values
    soiA_fun_unique_soiA = rle(event_statistics$soiA)$values
    if(length(soiA_fun_unique_years) != length(soiA_fun_unique_soiA)){
        stop('annual soiA series not correctly created')
    }
    # It is possible (though rare) there are no events in 1985. This will make
    # the lambda function undefined in 1985, which will cause a problem in the
    # fit since the observation start time is 1985. So we deal with that here
    if(min(soiA_fun_unique_years) != floor(obs_start_time)){
        print('Appending artificial first year to soiA series (should occur rarely)')
        soiA_fun_unique_years = c(floor(obs_start_time), soiA_fun_unique_years)
        soiA_fun_unique_soiA = c(0, soiA_fun_unique_soiA)
    }

    soiA_fun_data = cbind(soiA_fun_unique_years, soiA_fun_unique_soiA)
    soiA_fun_for_lambda = approxfun(soiA_fun_unique_years, soiA_fun_unique_soiA, 
        method='constant')

    rate_equation = initial_fit_env$best_nhp_model$rate_equation
    simulation_rate_equation = gsub(
        initial_fit_env$annual_rate_equations$soi$eqn, 
        'theta[1] + theta[2]*soiA_fun_for_lambda(floor(t))', 
        rate_equation,
        fixed=TRUE)
    rate_equation = simulation_rate_equation
    rm(simulation_rate_equation)

    # NOTE: Use initial_theta values based on the main fit.
    # Ensure that initial theta values will not give any event sequence 0
    # probability (since if the dataset is viewed as impossible given the
    # initial parameters, the optimization will fail).
    # The use of parscale helps the optimization converge in this specific case
    # but would need to be changed for other situations. See ?optim for
    # information.
    best_nhp_model = nhp$fit_nhpoisp(
        event_time,
        rate_equation=rate_equation,
        minimum_rate=0.0,
        initial_theta=initial_fit_env$best_nhp_model$par, 
        x0 = obs_start_time,
        event_durations = event_duration_years,
        number_of_passes = 3,
        optim_method=c('Nelder-Mead', 'Nelder-Mead', 'BFGS'),
        enforce_nonnegative_theta=FALSE,
        optim_control = initial_fit_env$best_nhp_model$optim_control)

    # If the above fit fails, exit gracefully
    if(is.na(best_nhp_model$par[1])){
        return(list(paste0('Fit failed on ', input_file)))
    }

    ###############################################################################
    #
    # Bootstrap marginal models
    #
    # Note that currently we use bootstrap models for uncertainty in copula's and
    # empirical distributions, but we use the original bayesian MCMC for
    # uncertainty in the mixture models.
    #
    ###############################################################################

    hsig_mixture_fit = DU$clone_environment(initial_fit_env$hsig_mixture_fit, 
        recursive=FALSE)
    duration_mixture_fit = DU$clone_environment(initial_fit_env$duration_mixture_fit, 
        recursive=FALSE)
    tideResid_mixture_fit = DU$clone_environment(initial_fit_env$tideResid_mixture_fit, 
        recursive=FALSE)

    #' 
    #' Get boostrap fits for hsig, duration, and tideResid
    #' Use the boostrapped data in event_statistics
    #'
    #' For starting values, use the best-fit parameters in initial_fit_env
    #' because these were used to generate the boostrap data so should never
    #' evaluate to Inf.... unless their threshold limit is exceeded by the 
    #' initial fit threshold, which is possible. So we change the latter if required
    #
    # Hsig
    #
    hsig_offset = hsig_threshold
    
    hsig_u_upper_limit = sort(event_statistics$hsig, decreasing=TRUE)[50] - hsig_offset
    hsig_u_lower_limit = hsig_threshold - hsig_offset # Should be 0

    start_par = initial_fit_env$hsig_mixture_fit$fit_optim$par
    if(start_par[3] > hsig_u_upper_limit) start_par[3] = hsig_u_upper_limit -0.01
    if(start_par[3] < hsig_u_lower_limit) start_par[3] = hsig_u_lower_limit +0.01
    evmix_fit$constrained_fit_gpd_mixture(
        hsig_mixture_fit,
        lower_bounds = c(-Inf, -Inf, hsig_u_lower_limit, -1000),
        upper_bounds = c( Inf,  Inf, hsig_u_upper_limit,  1000),
        start_par = start_par,
        data = event_statistics$hsig)

    #
    # Duration
    #
    duration_offset = ifelse(break_ties_with_jitter, 
        1 - as.numeric(initial_fit_env$default_jitter_amounts['duration']), 
        0.5) 

    duration_u_upper_limit = sort(event_statistics$duration, decreasing=TRUE)[50] - duration_offset
    duration_u_lower_limit = 0

    # Lower /upper bounds of parameter ranges
    lb = c(0, 0, duration_u_lower_limit, -Inf)
    ub = c(Inf, Inf, duration_u_upper_limit, Inf)

    start_par = initial_fit_env$duration_mixture_fit$fit_optim$par
    if( (start_par[3] > duration_u_upper_limit) | (start_par[3] < duration_u_lower_limit)){
        start_par[3] = 0.5*(duration_u_lower_limit + duration_u_upper_limit)
        start_par[4] = 0
    }

    evmix_fit$constrained_fit_gpd_mixture(
        duration_mixture_fit,
        lower_bounds = lb,
        upper_bounds = ub,
        start_par = start_par,
        data = event_statistics$duration)

    #   
    # tideResid 
    #

    # Limit the threshold
    tideResid_u_upper_limit = sort(event_statistics$tideResid, decreasing=TRUE)[50]
    min_tr = min(event_statistics$tideResid, na.rm=TRUE)
    tideResid_u_lower_limit = min_tr

    # Lower /upper bounds of parameter ranges
    lb = c(0, 0, tideResid_u_lower_limit, -Inf)
    ub = c(Inf, Inf, tideResid_u_upper_limit, Inf)
    # Starting parameters
    start_par = tideResid_mixture_fit$fit_optim$par
    if( (start_par[3] > tideResid_u_upper_limit) | (start_par[3] < tideResid_u_lower_limit)){
        start_par[3] = 0.5*(tideResid_u_lower_limit + tideResid_u_upper_limit)
        start_par[4] = 0
    }

    evmix_fit$constrained_fit_gpd_mixture(
        tideResid_mixture_fit,
        lower_bounds = lb,
        upper_bounds = ub,
        start_par = start_par,
        data = na.omit(event_statistics$tideResid))

    #
    # Steepness
    #

    # Use the 'old' logspline density estimator
    qsteepness0 = logspline::oldlogspline(
        na.omit(event_statistics$steepness), 
        lbound=min(event_statistics$steepness, na.rm=TRUE)-1.0e-04, 
        ubound=max(event_statistics$steepness, na.rm=TRUE)+1.0e-03)

    # Make 'raw' quantile and inverse quantile functions
    # With linear approximation we can be sure they are inverses of each other
    ptmp = seq(0, 1, len=1000)
    qtmp = logspline::qoldlogspline(ptmp, qsteepness0)
    # Quantile function
    qsteepness_raw = approxfun(ptmp, qtmp)
    # Inverse quantile function
    psteepness_raw = approxfun(qtmp, ptmp)  

    rm(ptmp, qtmp, qsteepness0)

    #
    # Direction
    # 

    # Fit logspline to direction distribution
    # Deliberately extend the range slightly beyond the data range
    qdir0 = logspline::oldlogspline(
        na.omit(event_statistics$dir), 
        lbound=min(event_statistics$dir, na.rm=TRUE)-0.5, 
        ubound=max(event_statistics$dir, na.rm=TRUE)+0.5)

    # Make 'raw' quantile and inverse quantile functions
    # With linear approximation we can be sure they are inverses of each other
    ptmp = seq(0, 1, len=1000)
    qtmp = logspline::qoldlogspline(seq(0, 1, len=1000), qdir0)
    # Quantile
    qdir_raw = approxfun(ptmp, qtmp)
    # Inverse quantile
    pdir_raw = approxfun(qtmp, ptmp)  
    # Cleanup
    rm(ptmp, qtmp, qdir0)
    
    # Get conditional fits
    #
    # BEWARE: In XXX_fit_conditional as modified below, q_raw and p_raw are no longer
    # inverse functions! 
    #
    # This is because the conditional q_raw uses is based on the randomly
    # sampled bayesian parameters for hsig, duration, and tideResid. OTOH p_raw
    # is based on a bootstrap fit [since herein it is only used to transform
    # data to (0-1) for copula estimation]
    #
    # This is appropriate for what I do below, since p_raw is ONLY used for
    # copula fits, and q_raw is ONLY used for the marginFun simulation (non-copula part).
    #
    # But BE CAREFUL WITH ANY OTHER USES, which would probably expect q_raw and
    # p_raw to be inverse functions!
    #

    # Hsig
    evmix_fit$make_random_mcmc_qfun(hsig_mixture_fit)

    hsig_fit_conditional = initial_fit_env$make_fit_conditional_on_season(
        event_statistics, 
        var = 'hsig',
        q_raw=hsig_mixture_fit$qfun_random, 
        p_raw=hsig_mixture_fit$pfun_constrained, 
        test_pqfun=FALSE, plot_quantiles=FALSE)

    # Duration
    evmix_fit$make_random_mcmc_qfun(duration_mixture_fit)

    duration_fit_conditional = initial_fit_env$make_fit_conditional_on_season(
        event_statistics, 
        var='duration',
        q_raw=duration_mixture_fit$qfun_random, 
        p_raw=duration_mixture_fit$pfun_constrained, 
        test_pqfun=FALSE, plot_quantiles=FALSE)

    # tideResid
    evmix_fit$make_random_mcmc_qfun(tideResid_mixture_fit)

    tideResid_fit_conditional = initial_fit_env$make_fit_conditional_on_season(
        event_statistics, 
        q_raw=tideResid_mixture_fit$qfun_random, 
        p_raw=tideResid_mixture_fit$pfun_constrained, 
        test_pqfun=FALSE, plot_quantiles=FALSE)

    # dir -- since the marginal is nonparametric we entirely base this
    # on the bootstrap data
    dir_fit_conditional = initial_fit_env$make_fit_conditional_on_soiA_and_season(
        event_statistics, 
        var = 'dir',
        q_raw = qdir_raw,
        p_raw = pdir_raw,
        test_pqfun=FALSE, plot_quantiles=FALSE)

    # steepness -- non parametric also
    steepness_fit_conditional = initial_fit_env$make_fit_conditional_on_season(
        event_statistics, 
        var = 'steepness',
        q_raw = qsteepness_raw,
        p_raw = psteepness_raw,
        test_pqfun=FALSE, plot_quantiles=FALSE)


    # Function to make the marginal distributions
    stormVarFun = initial_fit_env$make_stormVarFun(
        qduration = duration_fit_conditional$qfun,
        qhsig = hsig_fit_conditional$qfun,
        qtideResid = tideResid_fit_conditional$qfun,
        qdir = dir_fit_conditional$qfun,
        qsteepness = steepness_fit_conditional$qfun)


    ###############################################################################
    #
    # Fit Copula
    #
    ###############################################################################

    # Make ranked data with bootstrap
    conditional_variables = list(
        startyear=event_statistics$startyear, 
        soiA=event_statistics$soiA)

    # Convert to conditional probability values.
    # NOTE: This does not use the random bayesian p function -- 
    # It uses a the original fits p-function
    es_01 = data.frame( 
        duration=duration_fit_conditional$pfun(
            event_statistics$duration, conditional_variables),
        hsig=hsig_fit_conditional$pfun(event_statistics$hsig, 
            conditional_variables),
        dir=dir_fit_conditional$pfun(event_statistics$dir, 
            conditional_variables),
        tideResid=tideResid_fit_conditional$pfun(event_statistics$tideResid, 
            conditional_variables),
        steepness=steepness_fit_conditional$pfun(event_statistics$steepness,
            conditional_variables)
        )
    rm(conditional_variables)

    # Adapt format for CVine copula
    es_cop = as.copuladata(es_01)
    c_vine_node_order=c('hsig', 'duration', 'tideResid', 'steepness', 'dir')
    c_vine_order = match(c_vine_node_order, names(es_cop))
    es_cop_reorder = es_cop[,c_vine_order]

    # Update the copula parameters while keeping the originally selected structure.
    copula_model = initial_fit_env$make_Rvine_random_sampler(es_cop_reorder, 
        copula_fit=initial_fit_env$copula_model$copula_fit, plot=FALSE)
    random_copula_samples = copula_model$random_copula_samples


    ###############################################################################
    #
    # Fit the model
    #
    ###############################################################################

    event_creator = initial_fit_env$build_event_creator(
        random_copula_samples, 
        stormVarFun, 
        observation_start_time=obs_start_time, 
        lambda_rate_equation = initial_fit_env$simulation_rate_equation, 
        lambda_theta_par = best_nhp_model$par, 
        plot_soiA=FALSE,
        nyears_synthetic_series = nyears_synthetic_series)

    # Simulate the series

    # Note we add an extra gap of 'duration_gap_hours/(24*365.24)' between events
    # -- since with the current definition of events, no other event can occur
    # within this length of time. 
    synthetic_series = nhp$rnhpoisp(
        duration = nyears_synthetic_series, 
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

    synthetic_attr$tp1 = wavedisp$airy_period(
        lambda=synthetic_attr$hsig/synthetic_attr$steepness,
        h=initial_fit_env$buoy_depth)

    # Append sea level
    output_sl = initial_fit_env$compute_soi_MSL_perturbation(
        output_times = synthetic_attr$startyear, 
        output_soiA = synthetic_attr$soiA)

    synthetic_attr = cbind(synthetic_attr, data.frame(msl=output_sl))

    # Write out
    output_file = paste0(bootstrap_output_dir, '/', 'bootstrap_fit_series_', 
        iterate_number, '_', run_title_id, '.RDS')
    

    # Keep this just to save later
    #lambda_final = nhp$get_lambda_function(
    #        best_nhp_model$par, 
    #        rate_equation=best_nhp_model$rate_equation, 
    #        minimum_rate=0.)

    saveRDS(
        list(
            hsig_random_par=hsig_mixture_fit$random_par, 
            #hsig_fit_conditional=hsig_fit_conditional,
            bootstrap_hsig_mixture_fit = hsig_mixture_fit$fit_constrained$par,
            duration_random_par = duration_mixture_fit$random_par, 
            #duration_fit_conditional = duration_fit_conditional, 
            bootstrap_duration_mixture_fit = duration_mixture_fit$fit_constrained$par,
            tideResid_random_par = tideResid_mixture_fit$random_par, 
            #tideResid_fit_conditional = tideResid_fit_conditional, 
            bootstrap_tideResid_mixture_fit = tideResid_mixture_fit$fit_constrained$par,
            dir_fit_conditional = dir_fit_conditional, 
            steepness_fit_conditional = steepness_fit_conditional, 
            copula_model = copula_model, 
            #random_copula_samples = random_copula_samples,
            best_nhp_model = best_nhp_model, 
            soiA_fun_data = soiA_fun_data,
            synthetic_attr = synthetic_attr, 
            output_file = output_file, 
            event_statistics = event_statistics
        ),
        file = output_file
    )

    return(list(NULL))

}


#
#
# Run the bootstrap fits
#
#
#stop()

library(parallel)

# Set up random numbers to work in parallel
RNGkind("L'Ecuyer-CMRG")
#set.seed(12345)

# Wrap the 'single_series_bootstrap' function in a 'try'
# so we can catch errors without global failure
try_single_bootstrap<-function(random_series, iterate_number){

    output = try(single_bootstrap(random_series, iterate_number))    

    if(class(output) == 'try-error'){
        return(list(output, list(head(random_series)), list(iterate_number)))
    }else{
        return(output)
    }

}

X = mcmapply(try_single_bootstrap, 
    initial_fit_env$simulated_series_list, 
    as.list(1:length(initial_fit_env$simulated_series_list)),
    mc.cores=MC_CORES, mc.preschedule=TRUE)

# This should give us messages about errors
print(unlist(X))

#
#kk = which(is.na(initial_fit_env$event_statistics$soiA))
#single_bootstrap(initial_fit_env$event_statistics[-kk,], 99999)
