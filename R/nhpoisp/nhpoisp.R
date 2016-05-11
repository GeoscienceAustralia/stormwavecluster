# 
# Various functions for simulating / fitting non-homogeneous poisson processes
#
# AUTHORS
# Gareth Davies, Geoscience Australia 2014-15, gareth.davies.ga.code@gmail.com
# (Contributors add name here)
# 


#' Simulate a homogeneous poisson process with constant rate lambda (mean
#' occurrances per "year") for duration "years". 
#'
#' Note "year" can be any other time unit, and the result with be relative to
#' this
#'
#' @param duration number of "years" of simulation
#' @param lambda mean number of events per "year"
#' @return vector with times of synthetic events
#' @export
#'
rpoisp <-function(
    duration = 100, 
    lambda=1
    ){

    if(lambda<=0 | duration <= 0){
        stop('duration and lambda must be positive')
    }

    # We expect duration*lambda events over time duration, on average
    # Simulate more events than this to make it likely we make enough for
    # 'duration'
    event_count_multiplier = 1.5 
    n_sim_events = max(ceiling(duration*lambda*event_count_multiplier), 10)

    # Simulate the events
    spacing = rexp(n_sim_events, rate = lambda)
    events = cumsum(spacing)

    while(max(events) < duration){
        # Need to simulate more events to cover duration. We want to aviod this
        # for efficiency, hence why 'event_count_multiplier' is used
        spacing = c(spacing, rexp(n_sim_events, rate = lambda))
        events = cumsum(spacing)
    }

    if( sum(events < duration) == 0 ){
        # We can have no events within 'duration'. Cleanly deal with this by
        # returning a vector of length 0.
        return(numeric(0))
    }else{
        return(events[events < duration])
    }
}

###############################################################################
#'
#' Make the rate function lambda(t, tlast=-Inf) giving the instantaneous mean
#' number of events per 'year'
#'
#' The resulting lambda function takes an optional argument tlast, giving the
#' time of the last event
#' 
#' The function can be specified with text + a vector of parameters. 
#' This interface is useful for maximum likelihood estimation of parameters for
#' lambda(t). A minimum_rate can be used to prevent negative rates
#'
#' @param theta vector of parameters for rate function
#' @param rate_equation text of code to make the rate function lambda(t)
#'        (except clipping is applied to ensure the minimum rate is minimum_rate)  
#'        It can include references to t, to theta[1], theta[2], ... (up to the 
#'        length of theta) and also tlast, which is the time of the last event
#' @param minimum_rate minimum allowed rate. Any rate_equation values < this 
#'        are clipped to minimum_rate. Beware if this is too small
#'        and you have data occurring during a minimum_rate time, then
#'        it will be interpreted as an EXTREMELY unlikely event
#' @return the function lambda(t, tlast=-Inf)
#' @export
#'
get_lambda_function<-function(
    theta, 
    rate_equation='theta[1]+theta[2]*sin(2*pi*t)', 
    minimum_rate=0.0e-100
    ){

    rate_equation_expression = parse(text=rate_equation, keep.source=FALSE)
    lambda<-function(t, tlast=-Inf){
        raw_rate = eval(rate_equation_expression)
        #return( pmax(raw_rate, minimum_rate) )
        noclip = (raw_rate > minimum_rate)
        value = raw_rate * noclip + minimum_rate * (1 - noclip) 
        return( value )
    }
    
    # Store the actual equation with the function
    attr(lambda, 'rate_equation') = rate_equation
    attr(lambda, 'minimum_rate') = minimum_rate
    attr(lambda, 'theta') = theta

    return(lambda)
}

###############################################################################
#' Simulate a non-homogeneous poisson process with rate function lambda(t)
#'
#' The first event starts at t(1) when
#' F(t(1) ) = 1-exp(-integral(lambda(t) from observation_start_time to t(1))) 
#'              = uniform_random_number
#' Subsequent events start at t(i) when 
#' F(t(i) | t(i-1)+event_duration(i-1))
#' = 1-exp(-integral(lambda(t, tlast = t(i-1)+event_duration(i-1)) from (t(i-1)+event_duration(i-1)) to t(i))) 
#' = uniform_random_number
#' The integral is computed with the trapezoidal rule to solve for t(i)
#' 
#' @param duration the number of 'years' (or other time units) to simulate
#' @param lambda a function lambda(t) giving the rate at time t
#' @param event_properties_function a function which can compute a 
#'        (possibly random) event duration based on the time t. This allows
#'        events to have a finite duration, as has been suggested for coastal
#'        storm analysis. The values of this function will be attached to the 
#'        output as an attribute. Its default value gives a zero event duration.
#'        The function MUST return a list containing all numeric values, with a
#'        'duration' attribute. It can also have other attributes like Hsig,
#'        etc, which are probably correlated with duration.
#' @param integration_dt The trapezoidal integration time-step (in the same
#'        units as duration)
#' @param observation_start_time the time we begin observing. Until the first
#'        event occurs we treat tlast as -Inf when evaluating lambda (no clustering),
#'        so if this is important to you, simulate a 'dummy' bit of series initially
#'        and then cut it off.
#' @param extra_duration_gap This number is added to the event_duration computed from
#'        event_properties_function before all computations, so no event is 
#'        permitted before this additional time has elapsed. 
#'        However, the actual duration attribute of the returned value is
#'        not changed (i.e. it is still the result of event_properties_function).
#'        This is useful for datasets where events are distinguished to include
#'        a minimum 'gap' between events. 
#' @param print_progress Integer. Print the time every print_progress event
#' @return The synthetic event timeseries -- a sorted vector of random times 
#'         when events occurred (started). It has an attribute
#'         named 'event_properties', which is a numeric matrix with column names
#'         corresponding to the event properties
#' @export
#'
rnhpoisp<-function(
    duration = 100,
    lambda=get_lambda_function(
        c(1,1), 
        rate_equation='theta[1]+theta[2]*sin(2*pi*t)',
        minimum_rate=0.0),
    event_properties_function=function(t) { return(list(duration=0.0*t)) },
    integration_dt=1.0e-04,
    observation_start_time=0.,
    extra_duration_gap=0.,
    print_progress = 1e+100
    ){

    if(!is.function(lambda)){
        stop('lambda must be a function of time lambda(t, tlast = -Inf)')
    }

    # Make a sequence tseq to estimate the integral of lambda (for root finding)
    # The maxima of tseq is chosen so it is likely that we get an event
    # before the end of tseq (otherwise we iterate in a while loop)
    yearly_mean_lambda = integrate(lambda, observation_start_time, 
        observation_start_time + 1)$value
    dt = integration_dt
    tseq = seq(0, 1.0/yearly_mean_lambda, by=dt) 
    len_tseq = length(tseq)

    # Vector to store event times -- we don't know how many there will be so
    # make it 'likely big enough' to avoid later memory reallocation
    preallocation_storage_size = ceiling(yearly_mean_lambda*duration*2)
    event_storage_filler = rep(NA, len = preallocation_storage_size)
    event_storage = event_storage_filler

    # List to store results of event_properties_function -- pre-allocate some data
    # of the same size as event_storage
    event_properties = event_properties_function(observation_start_time)
    tmp = lapply(event_properties, f<-function(x) x*NA)
    event_properties_storage_filler = matrix(as.numeric(tmp), 
        nrow = preallocation_storage_size, ncol = length(tmp), byrow=TRUE) 
    colnames(event_properties_storage_filler) = names(tmp)
    event_properties_storage = event_properties_storage_filler

    # Assume there was no previous event at the start of the series
    number_of_events_kept = 0
    last_event_start_time = -Inf 
    last_event_end_time = -Inf

    while(last_event_end_time < observation_start_time + duration){

        # The next event occurs when 
        #
        # 1 - exp( - {integral of lambda(t, tlast = t[i-1] + event_duration[i-1]) from (t[i-1] + event_duration[i-1]) to t[i]}) 
        #     = uniform_random_number
        #
        # where t[i-1] is the time of the last event (or the start of
        # observatiions), event_duration[i-1] is its duration, and t[i] is the
        # time of the next event. 
        #
        # We need to find t(i) such that this is satisfied
        #
        # NOTE: t[i-1] + event_duration[i-1] = last_event_end_time in this code
        #
        # Equivalently, we need to find t[i] such that 
        # lambda_integeral = 0
        # where:
        # lambda_integral = 
        #    {integral of lambda(t) from (t[i-1]+event_duration[i-1]) to t[i]} + 
        #    log(1-uniform_random_number) 
        #
        # METHOD:
        # Compute lambda_integral along a section of the x-axis
        # of range(tseq), as a function of the x position
        #
        # If it becomes > 0., then find the root = t[i], which is the time of
        # the next event. 
        # Otherwise, we need to integrate further along the x-axis until we
        # find the root. Increase tseq by the max(tseq), add the integral along
        # this new section to the lambda_integral.  Proceed until we can find
        # the root

        uniform_random_number = runif(1)
        log_unif = log(1-uniform_random_number)
   
        # These are updated in the while loop
        prior_integral_time = 0
        lambda_integral = log_unif # < 0

        while(lambda_integral < 0.){

            # Evaluate lambda over a timeslice with duration/spacing from tseq
            # Starting time = [last event time + time-inverval over which we
            # have previously integrated in this while loop]
            if(number_of_events_kept == 0){
                # There was no previous event
                lambda_vals = lambda(observation_start_time + prior_integral_time + tseq)
            }else{
                lambda_vals = 
                    lambda(last_event_end_time + prior_integral_time + tseq, 
                        tlast = last_event_end_time)
            }

            if(min(lambda_vals)<0.){
                stop('Negative lambda values --> invalid lambda function')
            }

            # Compute trapezoidal approximation to integral for all t in tseq
            trapezoidal_darea = 
                (0.5*dt)*(lambda_vals[1:(len_tseq-1)] + lambda_vals[2:len_tseq])
            trapezoidal_integral = c(0, cumsum(trapezoidal_darea)) + lambda_integral

            # Update lambda_integral to include the newly integrated domain
            lambda_integral = trapezoidal_integral[len_tseq]

            if(lambda_integral < 0){
                # Update the prior integral time, and go back to the start of
                # the while loop 
                prior_integral_time = tseq[len_tseq] + prior_integral_time

            }else{
                # An event has occurred
                # Find the root. We will break out of the inner while loop

                if(lambda_integral == 0){
                    # Unlikely case that the integral is exactly 0 at the end
                    # of tseq

                    event_ind = max(which(lambda_vals>0))
                    if(event_ind == -Inf) stop('event_ind < -Inf')

                    event_spacing = tseq[event_ind] + prior_integral_time

                }else{

                    # Interpolate linearly over the change in sign to find root
                    # Note: m can never be 1 since 
                    # trapezoidal_integral[1] = 
                    #     0 + the previous value of lambda_integral

                    m = sum(trapezoidal_integral <= 0)

                    event_spacing = tseq[m] + prior_integral_time + 
                        dt*(-trapezoidal_integral[m])/
                        (trapezoidal_integral[m+1]-trapezoidal_integral[m])

                }

                # Record the start time of the event
                if(number_of_events_kept == 0){
                    last_event_start_time = event_spacing + observation_start_time
                }else{
                    last_event_start_time = last_event_end_time + event_spacing
                }
                number_of_events_kept = number_of_events_kept + 1

                if (number_of_events_kept %% print_progress == 0){
                    print(number_of_events_kept)
                    print(Sys.time())
                }

                # Get the end time of the event
                event_properties = event_properties_function(last_event_start_time)
                event_duration = event_properties$duration
                if(event_duration<0) stop('event duration < 0')

                # Here we denote the 'end' as the time after any extra_duration_gap
                last_event_end_time = last_event_start_time + event_duration + extra_duration_gap
    
                # Allocate memory for event storage if needed
                if(length(event_storage) < number_of_events_kept){
                    event_storage = c(event_storage, event_storage_filler)
                    event_properties_storage = rbind(event_properties_storage, event_properties_storage_filler)
                }

                # Store the start time + event properties
                event_storage[number_of_events_kept] = last_event_start_time
                event_properties_storage[number_of_events_kept, ] = 
                    as.numeric(event_properties)

            }
        }
    }

    # Event storage was pre-allocated, and might be too long (with extra NA
    # values). Remove NA's
    keepers = which(!is.na(event_storage))
    event_storage = event_storage[keepers]
    event_properties_storage = event_properties_storage[keepers, ,drop=FALSE]

    if(sum(event_storage < duration + observation_start_time) == 0){

        return(numeric(0))

    }else{

        # Return all events which occurred before 'duration', with an attribute
        # giving the event_properties
        keepers = which(event_storage < duration + observation_start_time)
        event_storage = event_storage[keepers]
        event_properties_storage = event_properties_storage[keepers, , drop=FALSE]
       
        # Convert event_properties to a matrix 
        #nc = length(event_properties_storage[[1]])
        #column_names = names(event_properties_storage[[1]])
        #event_properties_storage = matrix(unlist(event_properties_storage),
        #    ncol = nc, byrow=T)
        #colnames(event_properties_storage) = column_names
       
        # Return as an attribute 
        attr(event_storage, 'event_properties') = event_properties_storage

        return(event_storage)
    }
}

###############################################################################
#' 
#' Probability density function for a non-homogeneous possion process 
#'
#' If lambda is ever evaluated as < 0 at points x or in-between then
#' a message is printed, and density values of 0 (or -Inf for log=TRUE) are
#' returned [reflecting the 'invalid' lambda function]
#'
#' @param x the times at which events were observed, x[1] < x[2] < ... < x[n]
#' @param lambda function for the rate. lambda(t, tlast) is a function of time,
#'        and the (end) time of the last event 
#' @param x0 The starting time of observations, x0 <= x[1]
#' @param event_durations numeric vector of length(x) giving the durations of 
#'        the events. x[i] + event_durations[i] gives the end time of the ith 
#'        event. For a typical poisson process, event_durations[i] = 0. 
#' @param log TRUE/FALSE return the log pdf. Mainly convenient for computing 
#'        the log-likelihood
#' @param integration_dt The trapezoidal integration time-step in the same time-units as x
#' @return a vector of length = (length(x)). The i'th entry gives the 
#'         probability density of x[i] given x[i-1]. When i=1, we integrate
#'         from x0 but set tlast = -Inf for tlast in lambda 
#'         [i.e. lambda(t, tlast=-Inf)] 
#' 
#' @details The probability density of x[i] given x[i-1] and lambda is:
#'          lambda(x[i+1])*exp(-[integral_{x[i]+event_durations[i]}^{x[i+1]} lambda(t) dt])
#'          See e.g. Luceno et al., (2006) The effect of temporal dependence on
#'          the estimation of the frequency of extreme ocean climate events,
#'          Proc. R. Soc. A (2006) 462, 1683–1697
#' @export
#'
dnhpoisp<-function(
    x, 
    lambda=function(t, tlast=-Inf){NA},
    x0=0,
    event_durations = rep(0,length(x)),
    log=FALSE,
    integration_dt=1.0e-03
    ){

    lx = length(x)

    if(is.na(lambda(x0))){
        stop('Must provide a lambda function (and which does not evaluate to NA)')
    }

    # Various checks on data
    if(all(is.na(diff(x)))){
        stop('x must have at least 2 consecutive non NA values')
    }

    if(x0 > x[1]){
        stop('starting time x0 must be <= x[1]')
    }

    if(min(diff(x), na.rm=TRUE)<=0){
        stop('x must be increasing')
    }

    if(is.na(x[1]) | is.na(x[lx])){
        stop('x cannot have first/last values being NA')
    }

    if(length(event_durations)!=lx){
        stop('event_durations must have the same length as x')
    }

    if(max(event_durations)>0){
        if(any(diff(x) < event_durations[1:(lx-1)])){
            stop('Event duration cannot be < time to next event')
        }
    }

    if(min(event_durations)<0) stop('Cannot have negative event duration')

    dens = NA*x 
    # Append starting time to times (since the math works nicely that way)
    xnew = c(x0, x)
    event_durations_new = c(0, event_durations)
    for(i in 1:lx){ 

        # Skip 'NA' time values
        if(is.na(xnew[i]) | is.na(xnew[i+1]) ){
            next
        }

        # Compute 'last event time' for models with clustering
        if(i==1){
            # No previous value
            tlast = -Inf
        }else{
            tlast = xnew[i] + event_durations_new[i]
        }

        # Compute the integral of lambda from x[i] to x[i+1]
        trapezoidal_integral=TRUE
        if(trapezoidal_integral){
            # Trapezoidal integration on a fine grid 'xseq'
            min_xseq = xnew[i] + event_durations_new[i]
            max_xseq = xnew[i+1]
            len_xseq = ceiling((max_xseq-min_xseq)/integration_dt+1)
            xseq = seq(from=min_xseq, to=max_xseq, len=len_xseq)

            lambda_xseq = lambda(xseq, tlast = tlast)

            # Give impossible lambda functions a density of 0
            if(min(lambda_xseq)<0.){
                if(log==FALSE){ 
                    print('Negative lambda, returning 0')
                    return(rep(0, length(dens)))
                }else{
                    print('Negative lambda, returning -Inf')
                    return(rep(-Inf, length(dens)))
                }
            }

            dxseq = xseq[2]-xseq[1]
            rate_integral = sum(lambda_xseq)*dxseq - 
                0.5*dxseq*(lambda_xseq[1] + lambda_xseq[len_xseq])
        }else{
            rate_integral = 
                integrate(lambda, lower=xnew[i], upper=xnew[i+1], 
                    tlast=tlast)$value
        }

        # Compute log density for numerical stability
        dens[i] = log(lambda(xnew[i+1], tlast = tlast)) - rate_integral

    }

    if(log==FALSE) dens = exp(dens)

    return(dens)
}

###############################################################################
#'
#' Compute the negative log likelihood function for a non-homogeneous poisson
#' process with rate lambda(t)
#' 
#' @param x vector of times at which events occurred. It may contain NA values.
#'        If x[i] is NA, then the time-intervals x[i-1]:x[i] and x[i]:x[i+1] 
#'        are assigned an NA density (which is ignored in the summed negative
#'        log likelihood))
#' @param lambda function lambda(t) giving the rate at time t
#' @param x0 The observation start time, x0 < x[1]
#' @param integration_dt The increment used for numerical integration 
#'        in dnhpoisp
#' @return negative log likelihood 
#' @export
#'
negloglik_nhpoisp<-function(
    x,
    lambda=get_lambda_function(c(1,1)),
    x0=0,
    event_durations = rep(0,length(x)),
    integration_dt=1.0e-03
    ){

    negloglik = - sum(dnhpoisp(x, lambda=lambda, x0=x0, 
        event_durations = event_durations, log=TRUE, 
        integration_dt=integration_dt), na.rm=TRUE)
    # NOTE: the na.rm=TRUE above suggests the code can work
    # with some NA x values. This also occurs in dnhpoisp.
    # However, this is not well tested, and it is not
    # clear in what situation we would think it reasonable
    # to pass an NA time value. Maybe to denote a block of missing data? 
    # Consider revising.

    return(negloglik)
}


###############################################################################
#
#' Negative log likelihood as a function of the parameters theta
#' 
#' This form is useful for minimization & maximum likelihood estimation,
#' since the first argument is theta, as required by various nonlinear
#' optimizers
#'
#' @param theta vector of model parameters
#' @param observed_data sorted times of events
#' @param x0 the starting time of the observations
#' @param event_durations FIXME
#' @param rate_equation rate_equation used in get_lambda_function
#' @param minimum_rate minimum_rate used in get_lambda_function
#' @param enforce_nonnegative_theta if TRUE, any negative theta parameters
#'        cause -Inf to be returned
#' @param integration_dt the traezoidal integration dt used in dnhpoisson
#' @return negative log likelihood (theta | observed_data)
#' @export
negloglik_from_theta<-function(
    theta, 
    observed_data, 
    x0=0,
    event_durations = rep(0,length(observed_data)),
    rate_equation='theta[1] + theta[2]*sin(2*pi*t)',
    minimum_rate=0.0,
    enforce_nonnegative_theta=FALSE,
    integration_dt=1.0e-03){

    if(enforce_nonnegative_theta){
        # Prevent negative theta's
        if(any(theta<0)) return(-Inf)
    }

    lambda = get_lambda_function(theta, rate_equation=rate_equation, 
        minimum_rate = minimum_rate)

    return(negloglik_nhpoisp(observed_data, lambda=lambda, x0=x0, 
        event_durations=event_durations,
        integration_dt=integration_dt))
    
}

###############################################################################
#'
#' Fit a non-homogeneous possion model to the data
#'
#' @param observed_data vector of event times (sorted to be increasing)
#' @param rate_equation String of code defining the non-homogeneous rate
#'        lambda, in terms of 'theta' (a vector of parameters) and tlast (the
#'        time of the previous event)
#' @param minimum_rate If rate_equation allows values < minimum_rate, they are 
#'        clipped to minimum_rate. Note all rates must be >= 0.
#' @param initial_theta initial vector theta to begin minimisation. 
#'        theta contains parameters of rate_equation
#' @param x0 starting time of the observations
#' @param number_of_passes If greater than 1, then the optimization method is
#'        run this many times, with the theta starting values after the first
#'        run coming from the fitted parameters of the previous optimization.
#' @param first_pass_data_length if number_of_passes > 1, then on the first pass
#'        only fit the model to the first 'first_pass_data_length' number of data points
#'        This might help to get reasonable starting values for the next optimization?
#' @param enforce_nonnegative_theta If TRUE then any negative theta parameters
#'        cause -Inf to be returned in the negative log likelihood, so it is
#'        unlikely the optimizer will settle on negative values.
#' @param optim_method The optimization method to use, see ?optim
#' @param optimization control parameters, see ?optim
#' @return The result of a call to optim
#' @export
#'
fit_nhpoisp<-function(
    observed_data,
    rate_equation='theta[1] + theta[2]*sin(2*pi*t)',
    minimum_rate=0.0,
    initial_theta=c(1,1),
    x0=0,
    event_durations = rep(0, length(observed_data)),
    number_of_passes = 1,
    first_pass_data_length = Inf,
    enforce_nonnegative_theta=FALSE,
    optim_method='Nelder-Mead',
    optim_control=list(),
    verbose=FALSE,
    integration_dt=1.0e-04,
    use_optim2=FALSE){

    if(length(optim_method)>1){
        if(length(optim_method)!=number_of_passes){
            stop(paste0('optim_method can either be a single optimization \n',
                        ' method or a vector of length = number of passes'))
        }
    }else{
        optim_method = rep(optim_method, number_of_passes)
    }

    # Option to use the optim2 interface (experimental)
    if(use_optim2){
        optim_fun = optim2
    }else{
        optim_fun = optim
    }
    

    for(i in 1:number_of_passes){
        
        if(i>1) initial_theta = fit$par

        if(verbose){
            cat('\n')
            cat(c('    Iteration: ', i, '\n'))
            cat(c('    Optim method: ', optim_method[i], '\n'))
            cat(c('    Initial theta:', initial_theta, '\n'))
            cat('    Fitting....\n')
        }

        # If we use > 1 pass, then we probably have a hard-to-fit model.
        # Try using only a data subset on the first pass. 
        # Maybe this will help get a better starting guess for the 2nd iteration?
        if(number_of_passes>1 & i==1){
            ll = min(first_pass_data_length, length(observed_data))
            local_observed_data = observed_data[1:ll]
            local_event_durations = event_durations[1:ll]
        }else{
            local_observed_data = observed_data
            local_event_durations = event_durations
        }


        fit = try(optim_fun(initial_theta,
            fn = negloglik_from_theta, 
            # Other parameters for negloglik_from_theta below
            observed_data=local_observed_data,
            x0=x0,
            event_durations=local_event_durations,
            rate_equation=rate_equation,
            minimum_rate=minimum_rate,
            enforce_nonnegative_theta=enforce_nonnegative_theta,
            integration_dt=integration_dt,
            # Parameters to control optimization
            method=optim_method[i],
            control=optim_control,
            hessian=FALSE))

        if(class(fit)=='try-error'){
            fit = list(par=NA)
        }

        if(verbose){
            cat(c('    Fit ', i, '\n'))
            cat(c('    ', fit$par, '\n'))
        }
        
    }

    # Only add the hessian at the end
    fit_hessian = try( 
        optimHess(
            fit$par, 
            fn = negloglik_from_theta,
            gr = NULL,
            observed_data=observed_data,
            x0=x0,
            event_durations=event_durations,
            rate_equation=rate_equation,
            minimum_rate=minimum_rate,
            enforce_nonnegative_theta=enforce_nonnegative_theta,
            integration_dt=integration_dt
        )
    )

    if(class(fit_hessian)=='try-error'){
        fit$hessian = NA
    }else{
        fit$hessian = fit_hessian
    }

    fit$datalength = length(observed_data)
    fit$rate_equation = rate_equation
    fit$optim_control = optim_control

    return(fit)

}
   
############################################################################### 
#
#' Extract APPROXIMATE standard errors from a maximum likelihood fit, by
#' inverting the hessian. A practical approach for problems with many
#' parameters
#'
#' @param fit The output from fit_nhpoisp (or from optim)
#' @return a vector of standard errors (one for each fit$par)
#' @export
get_fit_standard_errors<-function(fit){
    ses = try(sqrt(diag(solve(fit$hessian))))

    if( (class(ses) == 'try-error') | any(is.na(ses))){
        print('Invalid standard errors produced: Use a more advanced method or improve the fit')
        return(NA)
    }

    return(ses)
} 



## #' Plot of likelihood
## plot_negloglik<-function(
##     x_ann,
##     theta_min,
##     theta_max,
##     nr = 19,
##     nc = 21,
##     lambda_function_minimum_rate = 0.,
##     lambda_function_rate_equation='theta[1]+theta[2]*sin(2*pi*t)'
##     ){
## 
##     # Presently only supports 2 parameter model
##     stopifnot(length(theta_min)==2)
## 
##     store_negloglik = matrix(NA,ncol=nc,nrow=nr)
## 
##     # Evaluate the negative log likelihood on this grid
##     consts = seq(theta_min[1], theta_max[1], len=nc)
##     amps  = seq(theta_min[2], theta_max[2], len=nr)
## 
##     for(fi in 1:nr){
##         print(fi)
##         for(ci in 1:nc){
##             lambda = get_lambda_function(c(consts[ci], amps[fi]), 
##                 rate_equation = lambda_function_rate_equation,
##                 minimum_rate = lambda_function_minimum_rate)
##             store_negloglik[fi,ci] = negloglik_nhpoisp(x_ann, lambda = lambda)
##         }
##     }
## 
##     # plot it
##     image(amps, consts, store_negloglik, col=rainbow(50))
##     contour(amps, consts, store_negloglik,add=T, nlevels=30)
##     # add a rough confidence ellipse
##     negloglikmin = min(store_negloglik)
##     loglik_CI_range = qchisq(0.95, 2)/2 # Joint confidence interval
##     #loglik_CI_range = qchisq(0.95, 1)/2 # confidence interval
##     contour(amps, consts, store_negloglik, add=T, 
##             levels = negloglikmin+ loglik_CI_range, lwd=3)
## }


##########################################################################

#' plot diagnostics for nhpoisson fit
#'
#' @param event_time Times of observations, in units of years
#' @param event_duration Durations of observations, in units of years. 
#' @param fitted_lambda The 'best fit' lambda function to compare with the data
#' @param num_simulated_duration Simulate this many series to compare with the
#'        data. Event durations for the simulated series will be created by resampling
#'        from the provided event_duration
#' @param nbins number of histogram bins
#' @return Nothing, but make a nice diagnostic plot
plot_nhpoisson_diagnostics<-function(
    event_time, 
    event_durations, 
    fitted_lambda, 
    num_simulated_series = 100, 
    nbins = 20){

    require(Matching)

    l = length(event_time)
    stopifnot(all(diff(event_time) >= event_durations[1:(l-1)]))

    par(mfrow = c(2,2))

    ########################################################################
    # Distribution of events within the year
    #
    yearly_breaks = seq(0, 1, len=nbins+1)
    # Compute hist so we can extend the ylim range
    data_hist = hist(event_time - floor(event_time), breaks=yearly_breaks, plot=FALSE)
    hist(event_time-floor(event_time), freq=FALSE, breaks=data_hist$breaks,
        main='Within year event frequencies \n (estimate of theoretical model in red, \n points from synthetic datasets)', 
        xlab='Time of year (as a decimal, 0 = start, 1 = end)',
        col='grey', ylim=c(0, max(data_hist$density)*1.5))

    starttime = min(event_time)
    sampleduration = diff(range(event_time))
    # Simulate to get the theoretical frequencies. To simulate, we need to know the
    # event durations, so we randomly sample these from the data
    synthetic_series = rnhpoisp(
        duration=sampleduration*num_simulated_series, 
        lambda=fitted_lambda, 
        event_properties_function=function(t) { 
            # This needs to return a list with a duration attribute
            durations = sample(event_durations, size=length(t), replace=TRUE)
            return(list(duration=durations))
        },
        observation_start_time = starttime 
        )

    # ks test that the distributions are similar
    print('KS TEST OF THE EVENTS TIME-OF-YEAR')
    #print(ks.test(event_time - floor(event_time), synthetic_series - floor(synthetic_series)))
    print(ks.boot(event_time - floor(event_time), synthetic_series - floor(synthetic_series)))

    # Split up the simulated series, and overplot the histogram
    for(i in 1:num_simulated_series){
        inds = which((synthetic_series >= starttime +(i-1)*sampleduration)&
                     (synthetic_series <= starttime + i*sampleduration))
        ss = synthetic_series[inds]
        synthetic_hist = hist(ss - floor(ss), breaks=data_hist$breaks, plot=FALSE)
        points(synthetic_hist$mids, synthetic_hist$density, pch='.', col='red')
    }

    # Do full curve
    #synthetic_hist = density(synthetic_series - floor(synthetic_series), 
    #    from=0, to=1) 
    #points(synthetic_hist$x, synthetic_hist$y, t='l', col='red')
    synthetic_hist = hist(synthetic_series - floor(synthetic_series), 
        breaks = data_hist$breaks, plot=FALSE) 
    #dx = diff(synthetic_hist$mids)[1]
    points(synthetic_hist$breaks, c(synthetic_hist$density,0), t='s', col='red')

   
    #########################################################
    # Time between events
 
    # Compute hist so we can extend the ylim range
    yearly_breaks_max = max(max(diff(event_time)), max(diff(synthetic_series)))
    yearly_breaks = seq(0, yearly_breaks_max, len=nbins+1)
    data_hist = hist(diff(event_time), breaks=yearly_breaks, plot=FALSE)

    hist(diff(event_time), freq=FALSE, breaks=data_hist$breaks, 
        main='Time between events \n (estimate of theoretical model in red, \n points from synthetic datasets)', 
        xlab='Time between event starting times (units of years)',
        col='grey',
        ylim=c(0, max(data_hist$density*1.5)))

    for(i in 1:num_simulated_series){
        inds = which((synthetic_series >= starttime +(i-1)*sampleduration)&
                     (synthetic_series <= starttime + i*sampleduration))
        ss = synthetic_series[inds]
        synthetic_hist = hist(diff(ss), breaks=data_hist$breaks, plot=FALSE)
        points(synthetic_hist$mids, synthetic_hist$density, pch = '.', col='red')
    }
    # Do full curve
    #synthetic_hist = density(diff(synthetic_series), from=0) 
    #points(synthetic_hist$x, synthetic_hist$y, t='l', col='red')
    synthetic_hist = hist(diff(synthetic_series), 
        breaks = data_hist$breaks, plot=FALSE) 
    points(synthetic_hist$breaks, c(synthetic_hist$density, 0.), t='s', col='red')

    # ks test that the distributions are similar
    print('KS TEST OF THE TIME BETWEEN EVENTS')
    #print(ks.test(diff(event_time), diff(synthetic_series)))
    print(ks.boot(diff(event_time), diff(synthetic_series)))

    ###########################################################
    # Count number of events in each year

    # For synthetic series
    obs_years = diff(range(floor(synthetic_series)))+1
    first_year = floor(min(synthetic_series))
    synthetic_events_per_year = tabulate(floor(synthetic_series) - first_year, 
        nbins=obs_years)

    ms = 0:max(synthetic_events_per_year)
    # FIXME: Must be a 1-line version of this -- count how many times each
    # integer occurs in synthetic_events_per_year
    msC = ms*0
    for(i in 1:length(ms)){
        msC[i] = sum(synthetic_events_per_year==ms[i])/length(synthetic_events_per_year)
    } 
    #synthetic_dens = density(synthetic_events_per_year, bw='SJ', from=0, to=m)

    # For data
    obs_years = diff(range(floor(event_time)))+1
    first_year = floor(min(event_time))
    events_per_year = tabulate(floor(event_time) - first_year + 1, 
        nbins=obs_years)

    print('KS TEST OF THE NUMBER OF EVENTS EACH YEAR')
    #print(ks.test(as.numeric(events_per_year), as.numeric(synthetic_events_per_year)))
    print(ks.boot(as.numeric(events_per_year), as.numeric(synthetic_events_per_year)))

    md = 0:max(events_per_year)
    mdC = md*0
    for(i in 1:length(mdC)){
        mdC[i] = sum(events_per_year==md[i])/length(events_per_year)
    }
    
    ylim_plot = c(0, max(mdC, msC))
    plot(ms-0.5, msC , t='s', xlab='Number of events per year',
        ylab='Density', main='Number of events each year (density)',
        ylim=ylim_plot, col='red')
    points(md-0.5, mdC, t='s',col='black')
    grid(col='brown')
    legend('topright', c('Data', 'Synthetic data (large sample)'), col=c('black', 'red'), 
        lty=c(1,1), bg='white')

    # Raw lambda plot
    t = seq(first_year, first_year + 3, len=1000)
    plot( t, lambda(t), t='l', xlab='Year')
    grid()
    title('Fitted lambda (tlast = -Inf)')
    
}

#' Compute the AIC and BIC for a fitted model. 
#'
#' @param fit Output of fit_nhpoisp
#' @return List with AIC and BIC
compute_fit_AIC_BIC<-function(fit){
   
    npar = length(fit$par) #+ 1
    ndata = fit$datalength
    negloglik = fit$value
   
    AIC = 2*negloglik + 2*npar 

    BIC = 2*negloglik + log(ndata)*npar

    return(list(AIC=AIC, BIC=BIC))
}



#'
#' Evaluate a lamdba function along a given timeseries
#' @param lambda The lambda function
#' @param times A sequence of times to evaluate lambda
#' @param observed_times A sequence of observed times
#' @param observed_durations The duration of events at observed_times
#' @return for all times[i], we return lambda(times[i], tlast=max((observed_times + observed_durations)[observed_times<times[i]])) 
#' 
#' @examples 
#' 
#' par(mfrow=c(3,1))
#' for(year in 1989:1991){
#'     mytimes = seq(year, year+1, by=1/5000)
#'     mylambda = nhp$evaluate_lambda_function(lambda, 
#'         times=mytimes, 
#'         observed_times = event_time, 
#'         observed_durations = event_duration)
#'     plot(mytimes, mylambda, t='l', lwd=3,
#'         xlab='Year', ylab='lambda')
#'     abline(v=event_time, col='red')
#'     abline(v=event_time + event_duration, col='blue')
#'     grid(col='grey', lty='dashed')
#' }
#' 
#' 
evaluate_lambda_function<-function(lambda, times, observed_times, observed_durations){

    lambda_vals = times*0.0

    if(min(diff(observed_times))<=0) stop('observed_times not increasing')
    if(min(diff(times))<= 0) stop('times not increasing')
    if(min(observed_durations)<0) stop('observed_durations are not >= 0')
    if(length(observed_times) != length(observed_durations)) stop('length(observed_times) is not = length(observed_durations)')

    # Special case when start < first observation
    if(times[1] < observed_times[1]){
        inds = which(times <= observed_times[1])
        if(length(inds) > 0) lambda_vals[inds] = lambda(times[inds])
       
        # Treat event gaps 
        inds = which((times > observed_times[1])&(times < observed_times[1] + observed_durations[1]))
        if(length(inds) > 0) lambda_vals[inds] = 0.
    }


    if(length(observed_times)==1) return(lambda_vals)

    # Typical case
    for(i in 2:length(observed_times)){
        inds = which((times >= observed_times[i-1] + observed_durations[i-1])&(times <= observed_times[i]))
        if(length(inds) > 0){
            lambda_vals[inds] = lambda(times[inds], tlast=observed_times[i-1] + observed_durations[i-1])
        }

        inds = which((times > observed_times[i]) & (times < observed_times[i] + observed_durations[i]))
        if(length(inds) > 0){
            lambda_vals[inds] = 0.
        }
    }
    
    return(lambda_vals)
}


#' Replacement for optim so we can also use optimx methods
#' Arguments are the same as for 'optim', except method
#' can be anything supported by optimx, and if method=NA
#' then all optimx methods are run
#' 
optim2<-function(par, fn, gr = NULL, ...,
           method = "Nelder-Mead",
           lower = -Inf, upper = Inf,
           control = list(), hessian = FALSE){

    library(optimx)

    if(is.na(method)){
        method<-c('Nelder-Mead',
             'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf',
             'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin')
    }

    optimx_fit = optimx::optimx(par=par, fn=fn, gr=gr, lower=lower, 
        upper=upper, method=method, hessian=hessian, control=control, ...)

    if(length(method) > 1){
        keep_method = which.min(optimx_fit$value)
    }else{
        keep_method = 1
    }
    
    par = as.numeric(coef(optimx_fit[keep_method,]))
    value = optimx_fit[keep_method,'value']
    convergence = optimx_fit[keep_method, 'convcode']
    if(hessian){
        hessian = optimHess(par, fn, gr, ..., control=control)
    }

    return(as.list(environment()))
}

