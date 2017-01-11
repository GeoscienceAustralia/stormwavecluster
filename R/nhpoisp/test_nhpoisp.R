# Source the nhpoisp.R code into its own environment, to keep the namespace
# clean. Functions in nhpoisp.R can now be called with the prefix 'nhp$'

nhp = new.env(parent = parent.env(.GlobalEnv))
source('nhpoisp.R', local=nhp)

#' Workhorse test function that rpoisp is correct by checking the statistics of
#' simulated events
.test_rpoisp<-function(
    duration = 10000,
    lambda = 1,
    reps = 1000
    ){

    #print(' Generating synthetic samples ...')
    length_rpoisp = replicate(reps, length(nhp$rpoisp(duration, lambda)) )

    # Expected mean number of events = duration*lambda
    err = (mean(length_rpoisp) - duration*lambda)
    abs_err = abs(err)

    #cat(c(' mean difference: ', err, '\n'))
    mean_difference_se = 1/sqrt(reps)*sqrt(lambda)*sqrt(duration)
    #cat(c(' expected mean difference standard error: ', 
    #      mean_difference_se, '\n'))

    # If everything works, this check should fail only 2/100000 times
    # (statistically)
    if(abs_err > abs(qnorm(1.0e-05))*mean_difference_se){
        success='FAIL -- should only occur ~ 2/1e+05 times by chance'
        #print('Issues, going into browser')
        #browser()
    }else{
        success='PASS'
    }

    # Variance of poisson = mean, so
    # sd(number_of_events) = sqrt(duration*expected_mean_number_of_events)
    expected_sd = sqrt(lambda*duration)
    #cat(c(' observed standard deviation of length_rpoisp ', 
    #      sd(length_rpoisp), '\n'))
    #cat(c(' expected standard deviation of length_rpoisp ', 
    #      expected_sd, '\n'))
   
    if(FALSE){ 
        # Make a plot of time-spacings, to check reasonableness
        dlr = diff(nhp$rpoisp(duration, lambda))
        if(length(dlr)>1){
            theoretical_x = seq(0, max(dlr)*1., len=100)
            theoretical_y = lambda*exp(-lambda*theoretical_x)
            dlr_hist = hist(dlr, plot=FALSE)

            plot(theoretical_x, theoretical_y ,
                t='l', col='red', 
                ylim=c(0,1.2)*max(c(theoretical_y, dlr_hist$dens)),
                main='Theoretical (red) and empirical from large random sample (black)')
            points(dlr_hist$mids, dlr_hist$dens,lwd=3,lend=1,type='h')
            title(sub=paste0('Duration: ', duration, ' Rate: ', lambda))
            abline(h=0,col='grey',lty='dotted')
        }else{
            plot(c(0,1),c(0,1),col=0,
                main=' Less than 3 points in sample -- normal for low (rate * duration)',
                xlab="", ylab="")
            title(sub=paste0('Duration: ', duration, ' Rate: ', lambda, 
                '\n rate*duration: ', duration*lambda))
        }
    }

    print(success)

}

#' To run some particular tests or rpoisp, call this
test_rpoisp<-function(){
    #par(mfrow=c(2,2))
    .test_rpoisp()
    .test_rpoisp(duration=5, lambda=1000, reps=500)
    .test_rpoisp(duration=50, lambda=40, reps=500)
    .test_rpoisp(duration=50000, lambda=0.3, reps=500)
}


#' Workhorse test function for rnhpoisp
#' 
#' Checks that the statistics of the simulated data really do behave as
#' expected, by comparing synthetic and analytical rates
#'
.test_rnhpoisp<-function(
    lambda_function){

    # Test parameters
    nSim = 1000
    series_duration=100
    mean_store = rep(NA,nSim)
    seasonal_breaks = seq(0,1,by=1/30) # Bin boundaries for time-of-year
    seasonal_mean_store = matrix(NA,nrow = nSim, 
        ncol = length(seasonal_breaks)-1)

    event_properties_function<-function(t) return(list(duration=0*t))

    #for(i in 1:nSim){
    #    if(i%%100==0) cat('.')
    single_run<-function(i){
        x = nhp$rnhpoisp(
                duration=series_duration, 
                lambda=lambda_function,
                event_properties_function = event_properties_function
            )

        # Basic tests
        stopifnot(max(x)<series_duration)
        stopifnot(min(x)>0)

        mean_store[i] = length(x)
        
        # Bin the rates
        x_frac = x - floor(x)
        time_of_year = cut(x_frac, breaks=seasonal_breaks, labels=FALSE)
        seasonal_mean_store_local = tabulate(time_of_year, 
            nbins=length(seasonal_breaks)-1)
        return(list(seasonal_mean_store=seasonal_mean_store_local))
    }
    
    library(parallel)
    all_runs = mclapply(as.list(1:nSim), single_run, mc.cores=12, 
        mc.preschedule=TRUE)
    for(i in 1:length(all_runs)){
        seasonal_mean_store[i,] = all_runs[[i]]$seasonal_mean_store
    }

    # Compute expected value of rates
    expected_rates = rep(NA, length(seasonal_breaks)-1)
    for(i in 1:length(expected_rates)){
        expected_rates[i] = integrate(lambda_function, seasonal_breaks[i], 
           seasonal_breaks[i+1])$value
    }

    # Visually compare simulated / analytical rates
    if(FALSE){
        boxplot(seasonal_mean_store)
        points(1:length(expected_rates), expected_rates*series_duration,
            col='red',pch=19)
        points(1:length(expected_rates), colMeans(seasonal_mean_store),
            col='green', pch=4)
        legend('topright', legend=c('Theoretical counts', 'Simulated mean counts'), 
            col=c('red', 'green'), pch=c(19,4))
    }

    # Test significance of deviations from expected rates
    # Since we do multiple tests, we want to adjust the p value
    desired_overall_p = 0.001
    f<-function(p) (p.adjust(p, 'hommel', 30) - desired_overall_p)**2
    required_p_single_test = optimize(f,lower=0,upper=0.01, tol=1.0e-08)$minimum

    for(i in 1:length(expected_rates)){
        # Test to see whether the synthetic results are consistent with the
        #test_stat = t.test(seasonal_mean_store[,i], 
        #    mu = expected_rates[i]*series_duration)
        test_stat = poisson.test(sum(seasonal_mean_store[,i]),
            r = expected_rates[i]*series_duration*nSim)

        # Check if the p.value is ok
        if(test_stat$p.value<required_p_single_test){
            print(i)
            print(test_stat)
            stop('mean number of events seems wrong (but there is a very small chance of a statistical fluke)')
        }else{
            print('PASS')
        }
    }
    
    #return(environment())
}

#'
#' Say we have a rate function lambda0(t) defining the background rate, but each
#' event has a duration=gap_size during which no other events can occur. We want
#' to compute the true CDF of the time between events accounting for these gaps. 
#' Note this is no-longer a poisson process.
#'
#' METHOD:
#'     The CDF of time between events for the 'raw' rate function is:
#' F(t_i | t_(i-1)) = 1 - exp{-integral(lambda0 from t_(i-1) to t_i)}
#'
#'     When there are gaps, this becomes:
#' F(t_i | t_(i-1)) 
#' = 1 - exp{-integral(lambda0 from t_(i-1)+GAP to t_i} %% When t_i>GAP+t_(i-1)
#' = 0 %% When t_i <= t_(i-1)+GAP
#'
#'
.test_rnhpoisp_constantgap<-function(){

    lambda<-function(t, tlast) return(t*0+1)
    gap = 100/365
    x = nhp$rnhpoisp(10000, lambda=lambda, 
        event_properties_function=function(t) list(duration=t*0+gap))

    delays = diff(x) - gap
    if(min(delays)<0){
        stop('time between events was < event_duration')
    }    

    compare1 = sort(rexp(length(delays)))
    sorted_delays = sort(delays)
    if(FALSE){
        qqplot(sorted_delays, compare1, log='xy')
        abline(0,1,col='red')
        title('Q-Q plot of post-gap delays vs an exponential random variable \n grey points are simulations')
    }

    sim_min = rep(Inf, length(delays))
    sim_max = rep(0, length(delays))
    nsim = 2000
    for(i in 1:nsim){
        compare_local = sort(rexp(length(delays)))
        sim_min = pmin(sim_min, compare_local) 
        sim_max = pmax(sim_max, compare_local)
        if(i<nsim/10){
            if(FALSE) points(compare_local, compare1, col='grey',pch=19,cex=0.2) 
        }
    }
    if(FALSE) points(sort(delays), compare1,col='red',pch=19,cex=0.2)

    if(all((sorted_delays - sim_min)*(sim_max - sorted_delays) > 0)){
        print('PASS')
    }else{
        print('FAIL')
        print(paste0('delays outside the range of ', nsim, 
            ' synthetic qqplots'))
    }

}

#' Run various tests
test_rnhpoisp<-function(){

    print('## Test 1 ## ')
    lambda = nhp$get_lambda_function(c(1,1, 0.1), 
        rate_equation='theta[1]+theta[2]*sin(2*pi*(t+theta[3]))',
        minimum_rate=0.0)
    test_result = .test_rnhpoisp(lambda_function = lambda)
    
    print('## Test 2 ## ')
    lambda = nhp$get_lambda_function(c(1,2,0.3), 
        rate_equation='theta[1]+theta[2]*sin(2*pi*(t+theta[3]))',
        minimum_rate=0.0)
    test_result = .test_rnhpoisp(lambda_function = lambda)
    
    print('## Test 3 ## ')
    lambda = nhp$get_lambda_function(c(1,2,0.3, 1, 0.5), 
        rate_equation='theta[1]+theta[2]*sin(2*pi*(t+theta[3]))+theta[4]*sin(4*pi*(t+theta[5]))',
        minimum_rate=0.0)
    test_result = .test_rnhpoisp(lambda_function = lambda)

    print('## Test 4 ##')
    .test_rnhpoisp_constantgap()

}

#' Convenience function to get the current .Random.seed (or create it if it
#' doesn't exist)
#'
get_random_seed<-function(){
    if(exists('.Random.seed')){
        return(.Random.seed)
    }else{
        # Need to generate a random number to create the seed
        x = runif(1)
        return(.Random.seed)
    }
}

#'
#' Test the method for fitting non-homogeneous poisson. We simulate the process, 
#' then back-estimate the parameters. Tweaks were required (selection of
#' optim_method, number_of_passes, etc) to get it to work, as is not suprising
#' for nonlinear optimization
#'
#' Changing the series duration / comparison tol will change the pass / failure
#' In practical applications, you need to tweak the initial_theta, optim
#' method, etc
#'
#' @param series_duration length(years) of sythetic series we test
#' @param seed Reproducible parameter for set.seed(seed)
#' @return The function environment
#'
test_fit_nhpoisp<-function(series_duration=100, seed = 1, cases=1:7){

    # Reproducible random numbers
    original_seed = get_random_seed()
    set.seed(seed)

    if(1 %in% cases){
        ########################################################################### 
        # Make some synthetic data 
        theta_par = c(1, 2, 0.1)
        rate_equation = 'theta[1]+theta[2]*sin(2*pi*(t+theta[3]))'

        lambda = nhp$get_lambda_function(theta_par, rate_equation, minimum_rate=0)
        observed_data = nhp$rnhpoisp(duration=series_duration, lambda)

        # Fit the model
        # Initial theta must be consistent with the data (i.e. finite neg log
        # likelihood)
        fit1 = nhp$fit_nhpoisp(observed_data, 
            rate_equation=rate_equation,
            minimum_rate=0.,
            initial_theta=c(1.2, 2, 0.1), 
            x0=0,
            number_of_passes=1,
            enforce_nonnegative_theta=TRUE,
            verbose=FALSE)
           
        # Check that values are within 10% 
        fitted_par = fit1$par
        # The phase variables can be altered modulo 1 without changing things
        # so for the test we adjust 
        fitted_par[3] = fitted_par[3]%%1

        err = fitted_par - theta_par

        rel_err = abs(err)/abs(nhp$get_fit_standard_errors(fit1))

        if(any(rel_err > 3)){
            # Compare the negative log likelihood of the theta_par and fitted_par
            cat('FAIL case 1 \n') 
        }else{
            cat('PASS\n')
        }

    }

    if(2 %in% cases){
        ###########################################################################
        #
        # Another example below with more harmonics
        #
        ###########################################################################

        # Make some synthetic data 
        theta_par = c(1, 2, 0.1, 0.4, 0.2)
        rate_equation = 'theta[1]+theta[2]*sin(2*pi*(t+theta[3])) + theta[4]*sin(4*pi*(t+theta[5]))'

        lambda = nhp$get_lambda_function(theta_par, rate_equation, minimum_rate=0)
        observed_data = nhp$rnhpoisp(duration=series_duration, lambda = lambda)

        # Fit the model
        # Initial theta must be consistent with the data (i.e. finite neg log
        # likelihood)
        fit2 = nhp$fit_nhpoisp(observed_data, 
            rate_equation=rate_equation,
            minimum_rate=0.,
            initial_theta=c(1.0, 2.0, 0.1, 0.4, 0.2),
            x0=0,
            number_of_passes=2,
            enforce_nonnegative_theta=TRUE,
            optim_method=c('Nelder-Mead','BFGS'),
            verbose=FALSE)
        
        # Check that values are within comparison_tol
        fitted_par = fit2$par
        # The phase variables can vary modulo 1 or 0.5, so adjust them for the test
        fitted_par[3] = fitted_par[3]%%1
        fitted_par[5] = fitted_par[5]%%0.5

        err = fitted_par - theta_par

        rel_err = abs(err)/nhp$get_fit_standard_errors(fit2)

        if(any(rel_err>3)){
            # Compare the negative log likelihood of the theta_par and fitted_par
            cat('FAIL on case 2 \n') 
        }else{
            cat('PASS \n')
        }

    }


    if(3 %in% cases){
        ###########################################################################
        #
        # Another example below -- with clustering
        #
        ###########################################################################

        # Make some synthetic data 
        theta_par = c(1, 2, 15.0, 1/36)
        rate_equation = 'theta[1]+theta[2]*sin(2*pi*(t)) + theta[3]*exp(-(t-tlast)/theta[4])'

        lambda = nhp$get_lambda_function(theta_par, rate_equation, minimum_rate=0.)
        observed_data = nhp$rnhpoisp(duration=series_duration, lambda = lambda)

        # Fit the model
        # Initial theta must be consistent with the data (i.e. finite neg log
        # likelihood)
        fit3 = nhp$fit_nhpoisp(observed_data, 
            rate_equation=rate_equation,
            minimum_rate=0.,
            initial_theta=c(1., 2.0, 15.0, 1/36),
            x0=0,
            number_of_passes=2,
            enforce_nonnegative_theta=TRUE,
            optim_method=c('Nelder-Mead', 'BFGS'),
            optim_control=list(maxit = 500),
            verbose=FALSE)
        
        # Check that values are within comparison_tol
        fitted_par = fit3$par

        err = fitted_par - theta_par

        rel_err = abs(err)/nhp$get_fit_standard_errors(fit3)

        if(any(rel_err > 3)){
                cat('FAIL on case 3\n') 
        }else{
            cat('PASS\n')
        }
    }
   
    if(4 %in% cases){ 

        ###########################################################################
        #
        # Another more complex clustering example
        #
        ###########################################################################

        # Make some synthetic data 
        theta_par = c(1, 2, 0.1, 0.7, 0.2, 15.0, 1/50)
        rate_equation = 'theta[1]+theta[2]*sin(2*pi*(t+theta[3])) + theta[4]*sin(4*pi*(t+theta[5])) + theta[6]*exp(-(t-tlast)/theta[7])'

        lambda = nhp$get_lambda_function(theta_par, rate_equation, minimum_rate=0)
        observed_data = nhp$rnhpoisp(duration=series_duration, lambda = lambda)

        # Fit the model
        # Initial theta must be consistent with the data (i.e. finite neg log
        # likelihood)
        fit4 = nhp$fit_nhpoisp(observed_data, 
            rate_equation=rate_equation,
            minimum_rate=0.,
            initial_theta=c(1., 2.0, 0.1, 0.7, 0.2, 15., 1/50),
            x0=0,
            number_of_passes=2,
            enforce_nonnegative_theta=TRUE,
            optim_method=c('Nelder-Mead', 'BFGS'),
            optim_control=list(maxit = 500),
            verbose=FALSE)
        
        # Check that values are within comparison_tol
        fitted_par = fit4$par
        # The phase variables can vary modulo 1 or 0.5, so adjust them for the test
        fitted_par[3] = fitted_par[3]%%1
        fitted_par[5] = fitted_par[5]%%0.5

        err = fitted_par - theta_par
        rel_err = abs(err)/nhp$get_fit_standard_errors(fit4)

        if(any(rel_err > 3)){
            cat('FAIL on case 4\n')
        }else{
            cat('PASS\n')
        }
    }

    if(5 %in% cases){
        ###########################################################################
        #
        # Example with clustering and gaps
        #
        ###########################################################################

        # Make some synthetic data 
        theta_par = c(1, 2, 0.1, 15.0, 1/50)
        rate_equation = 'theta[1]+theta[2]*sin(2*pi*(t+theta[3])) + theta[4]*exp(-(t-tlast)/theta[5])'

        lambda = nhp$get_lambda_function(theta_par, rate_equation, minimum_rate=0)

        # Make the event duration vary uniformly
        event_properties_function<-function(t){
            return(list(duration = runif(n=length(t))*0.013))
        }
        
        observed_data = nhp$rnhpoisp(duration=series_duration, lambda = lambda, 
            event_properties_function = event_properties_function)

        event_durations = attr(observed_data, 'event_properties')[,1]

        # Fit the model
        # Initial theta must be consistent with the data (i.e. finite neg log
        # likelihood)
        fit5 = nhp$fit_nhpoisp(observed_data, 
            rate_equation=rate_equation,
            minimum_rate=0.,
            initial_theta=c(1., 2.0, 0.1, 15.0, 1/50),
            x0=0,
            event_durations=event_durations,
            number_of_passes=2,
            enforce_nonnegative_theta=TRUE,
            optim_method=c('Nelder-Mead', 'BFGS'),
            optim_control=list(maxit = 500),
            verbose=FALSE)
        
        # Check that values are within comparison_tol
        fitted_par = fit5$par
        # The phase variables can vary modulo 1 or 0.5, so adjust them for the test
        fitted_par[3] = fitted_par[3]%%1

        err = fitted_par - theta_par
        rel_err = abs(err)/nhp$get_fit_standard_errors(fit5)

        if(any(rel_err > 3)){
            cat('FAIL on case 5\n')
        }else{
            cat('PASS\n')
        }
    }

    if(6 %in% cases){
        ###########################################################################
        #
        # Example like the previous one, with bigger gaps, more likely to show any
        # problems with gap implementation
        #
        ###########################################################################

        # Make some synthetic data 
        theta_par = c(1, 2, 0.1, 15.0, 1/50)
        rate_equation = 'theta[1]+theta[2]*sin(2*pi*(t+theta[3])) + theta[4]*exp(-(t-tlast)/theta[5])'

        lambda = nhp$get_lambda_function(theta_par, rate_equation, minimum_rate=0)

        # Make the event duration vary uniformly
        event_properties_function<-function(t){
            return(list(duration = runif(n=length(t))*0.5))
        }
        
        observed_data = nhp$rnhpoisp(duration=series_duration, lambda = lambda, 
            event_properties_function = event_properties_function)

        event_durations = attr(observed_data, 'event_properties')[,1]

        # Fit the model
        # Initial theta must be consistent with the data (i.e. finite neg log
        # likelihood)
        fit6 = nhp$fit_nhpoisp(observed_data, 
            rate_equation=rate_equation,
            minimum_rate=0.,
            initial_theta=c(1., 2.0, 0.1, 15., 1/50),
            x0=0,
            event_durations=event_durations,
            number_of_passes=2,
            enforce_nonnegative_theta=TRUE,
            optim_method=c('Nelder-Mead', 'BFGS'),
            optim_control=list(maxit = 500),
            verbose=FALSE)
        
        # Check that values are within comparison_tol
        fitted_par = fit6$par
        # The phase variables can vary modulo 1 or 0.5, so adjust them for the test
        fitted_par[3] = fitted_par[3]%%1

        err = fitted_par - theta_par
        rel_err = abs(err)/nhp$get_fit_standard_errors(fit6)

        if(any(rel_err > 3)){
            cat('FAIL on case 6\n')
        }else{
            cat('PASS\n')
        }

    }

    if(7 %in% cases){
        ###########################################################################
        #
        # Example like a dataset we fit
        #
        ###########################################################################

        # Make some synthetic data 
        theta_par = c(26, 10.5, 0.5, 29.0, 165.0)

        triangle_eqn = 'abs(2/pi * asin( cos(pi*(t + 0.5 - theta[3])) ))'
        cluster_eqn = '- theta[4]*exp((tlast - t)*theta[5])'
        rate_equation = paste0('theta[1] - theta[2]*', triangle_eqn, cluster_eqn)

        lambda = nhp$get_lambda_function(theta_par, rate_equation, minimum_rate=0)

        # Make the event duration vary uniformly
        event_properties_function<-function(t){
            return(list(duration = runif(n=length(t))*6/365.25))
        }
        
        observed_data = nhp$rnhpoisp(duration=series_duration, lambda = lambda, 
            event_properties_function = event_properties_function)

        event_durations = attr(observed_data, 'event_properties')[,1]

        # Fit the model
        # Initial theta must be consistent with the data (i.e. finite neg log
        # likelihood)
        fit7 = nhp$fit_nhpoisp(observed_data, 
            rate_equation=rate_equation,
            minimum_rate=0.,
            initial_theta=c(25., 10, 0.5, 0., 140.0),
            x0=0,
            event_durations=event_durations,
            number_of_passes=2,
            enforce_nonnegative_theta=TRUE,
            optim_method=c('Nelder-Mead', 'BFGS'),
            optim_control=list(maxit = 500, parscale=c(1, 1, 0.1, 1, 10)),
            verbose=FALSE)
        
        # Check that values are within comparison_tol
        fitted_par = fit7$par
        # The phase variables can vary modulo 1 or 0.5, so adjust them for the test
        fitted_par[3] = fitted_par[3]%%1

        err = fitted_par - theta_par
        rel_err = abs(err)/nhp$get_fit_standard_errors(fit7)

        if(any(rel_err > 3)){
            cat('FAIL on case 7\n')
        }else{
            cat('PASS\n')
        }

    }

    # Reset the random seed
    .Random.seed = original_seed

    #return(environment())
}

#' Apply the test_fit_nhpoisp function in parallel
test_cases_parallel<-function(series_duration = 1000, seed = 1, cases = 1:7, 
    MC_CORES = 6){

    library(parallel)

    test_case<-function(i){
       test_fit_nhpoisp(series_duration = series_duration, 
                        seed = seed, 
                        cases=i) 
    }

    results = mclapply(as.list(cases), test_case, mc.cores = MC_CORES, 
        mc.preschedule=TRUE) 

    #return(results)
}

#' This test simulates many series, then back-fits their parameters, and
#' checks that the accuracy is consistent with expectations from likelihood
#' theory. Because it is a simulation test there is a small chance of failure
#' even if everything is working correctly (around 1% for each check). If it
#' fails, run it again repeatedly to confirm it fails more than expected before
#' assuming there are problems. 
#'
#'
#' Note we need to use the joint confidence region (2 degrees of
#' freedom) to cover the 'real' synthetic lambda parameters with desired
#' frequency (as expected)
#'
#' However, for confidence intervals on any single parameter, a single degree of
#' freedom should be ok
test_fit_simulate_and_confidence_interval<-function(){

    # The model we fit
    synthetic_lambda_parameters = c(2.6, 2.4)
    rate_equation = "theta[1]+theta[2]*sin(2*pi*t)"

    df = length(synthetic_lambda_parameters) # 1  
    desired_coverage = 0.95
    
    # Parameter estimates seem to be more bias for low series
    # duration
    # Makes sense as the justification of ML is asymptotic
    series_duration = 100 

    nsim = 300
    CI_covered = rep(NA, nsim)
    store_fitpar = matrix(NA, ncol=df, nrow=nsim)

    # Function to run the test (to ease parallel treatment)
    single_test<-function(i){

        # Make a synthetic series
        lambda_function = nhp$get_lambda_function(
            theta=synthetic_lambda_parameters,
            rate_equation=rate_equation)
        x_ann = nhp$rnhpoisp(
            duration=series_duration, lambda = lambda_function)

        # Estimate its parameters with a bad initial guess
        fit = nhp$fit_nhpoisp(
            observed_data = x_ann, 
            rate_equation = rate_equation,
            initial_theta = c(1, 1))
            
        store_fitpar = fit$par

        # Check whether the 'real' parameter values are within the confidence
        # interval (1) or not (0)
        loglik_region_max = fit$value + qchisq(desired_coverage, df)/2.
        if(loglik_region_max > 
            nhp$negloglik_from_theta(synthetic_lambda_parameters, 
                observed_data=x_ann)
            ){
            CI_covered = 1
        }else{
            CI_covered = 0
        }

        return(list(store_fitpar=store_fitpar, CI_covered=CI_covered))

    }
    
    library(parallel)
    all_fits = mclapply(as.list(1:nsim), single_test, mc.cores=12, 
        mc.preschedule=TRUE)

    store_fitpar = matrix(
        unlist(lapply(all_fits, f<-function(x) x$store_fitpar)),
        ncol=2, byrow=TRUE)

    CI_covered = unlist(lapply(all_fits, f<-function(x) x$CI_covered))

    # How often do we cover? Should be close to desired_coverage
    coverage_test = 
        binom.test(sum(CI_covered), length(CI_covered), p=desired_coverage,
        conf.level=0.99)

    if(coverage_test$p.value > 0.01){
        print('PASS')
    }else{
        print("FAIL (on CI coverage).")
        print(coverage_test)
    }

    # Check it

    par1_test = t.test(store_fitpar[,1], mu = synthetic_lambda_parameters[1],
        conf.level=0.99)
    if(par1_test$p.value > 0.01){
        print('PASS')
    }else{
        print('FAIL (on fit of par1)')
        print(par1_test)
    }

    par2_test = t.test(store_fitpar[,2], mu = synthetic_lambda_parameters[2],
        conf.level=0.99)
    if(par2_test$p.value > 0.01){
        print('PASS')
    }else{
        print('FAIL (on fit of par2)')
        print(par2_test)
    }
    
}

#' Run all tests
test_all<-function(){

    #
    RNGkind("L'Ecuyer-CMRG")
    set.seed(12345)

    cat('\n')
    cat('Testing rpoisp....\n')
    cat('\n')
    test_rpoisp()

    cat('\n')
    cat('Testing rnhpoisp....\n')
    cat('(this takes a few minutes on a 12-core linux box)\n')
    cat('(It assumes you are running on a shared memory linux machine)\n')
    cat('\n')
    test_rnhpoisp()

    cat('\n')
    cat('Testing we can back-fit simulated series with expected accuracy ...\n')
    cat('(this takes a few minutes on a 12-core linux box)\n')
    cat('(It assumes you are running on a shared memory linux machine)\n')
    cat('')

    test_fit_simulate_and_confidence_interval()
    
    cat('\n')
    cat('Testing back-fitting of simulated series with complex rate functions ...\n')
    cat('(this takes a few minutes on a 12-core linux box)\n')
    cat('(It assumes you are running on a shared memory linux machine)\n')
    cat('')

    test_cases_parallel(series_duration=100, MC_CORES=7)

}

test_all()
