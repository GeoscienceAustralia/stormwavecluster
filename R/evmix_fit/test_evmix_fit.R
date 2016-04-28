source('evmix_fit.R', local=TRUE)

#' Simple test code for fit_gpd_mixture
#' This is problematic without a VERY large sample size,
#' probably because the 'u' parameter is hard to constrain 
#' from the data alone (since the tail gpd is a similar shape to the
#' bulk gamma). 
#test_fit_gpd_mixture<-function(){
#
#    set.seed(1)
#
#    gshape= 1
#    gscale= 1
#    xi = 0.05
#    u = qgamma(0.8, gshape, 1/gscale)
#
#    ## Example like one of my hsig fits
#    #gshape= 0.84
#    #gscale= 1.01
#    #xi = -0.22
#    #u = 1.27
#    
#    # The sample size needs to be large for the gpd part to be estimated ok.
#    #sample_size = 5e+05 # Parameter estimates ok
#    sample_size = 5e+04 
#    #sample_size = 500 # Threshold very biased downward, gpd shape parameter not great
#    
#    test_data = qgammagpdcon(runif(sample_size), gshape = gshape, gscale = gscale, 
#        u = u, xi = xi, phiu = TRUE)
#
#    test_fit = fit_gpd_mixture(data = test_data, 
#        gpd_threshold_quantile_range=c(0.5, 0.95), verbose=FALSE)
#
#    true_par = c(gshape, gscale, u, xi)
#
#    if(all(abs(test_fit$fit_optim$par - true_par) < 
#               c(0.01, 0.01, 0.2, 0.3)*true_par)){
#        print('PASS')
#    }else{
#        print('FAIL: Fit is too far from true parameters')
#    }
#
#}


#' Function to fit some random data and make a plot
test_mcmc_gpd_mixture<-function(myseed=1){

    set.seed(myseed) 

    gshape=1
    gscale=1
    xi = 0.05
    u = qgamma(0.8, gshape, 1/gscale)

    true_par = c(gshape, gscale, u, xi)

    # The sample size needs to be large for the gpd part to be estimated ok.
    # If it's not large, we expect large uncertainties -- that's ok
    # But for a test it's cleaner / clearer if we use lots of data
    test_data = qgammagpdcon(runif(50000), gshape = gshape, gscale = gscale, 
        u = u, xi = xi)

    # Fit the data
    test_fit = fit_gpd_mixture(data = test_data, ntrial_u = 20, 
        gpd_threshold_quantile_range = c(0.6, 0.9), ntrial_iterations=10,
        verbose=FALSE)

    mcmc_gpd_mixture(test_fit, 
        par_lower_limits = c(0, 0, quantile(test_data, 0.5), -0.5), 
        par_upper_limits = c(1e+08, 1e+08, quantile(test_data, 0.9), 0.5), 
        mcmc_start_perturbation=0.05, mcmc_length=1e+05, mcmc_thin=100,
        mcmc_burnin=1000, mcmc_nchains=1, mcmc_tune=1,
        mc_cores=1, annual_event_rate=1)

    # Make a plot
    png('test_mcmc_gpd_mixture.png', width=10, height=8, res=200, units='in')

    mcmc_rl_plot(test_fit)

    # Add the 'true' parameters
    true_quantiles = test_fit$qf(1-test_fit$desired_rates, 
        c(gshape, gscale, u, xi))
    points(test_fit$desired_rates, true_quantiles, t='l', lwd=3, col='purple')

    dev.off()

    # Test that the true values are within the confidence limits
    enveloped = ((true_quantiles > test_fit$desired_lower_q) &
        (true_quantiles < test_fit$desired_upper_q))

    # With a seed of 1 this should pass, although of course if we tried
    # enough other seeds it should fail at some stage
    if(all(enveloped)){
        print('PASS')
    }else{
        print('FAIL')
    }

}

#' A few test cases for the gpd_mixture fit routine
#' We repeatedly generate synthetic data with known parameters,
#' then estimate the parameters with the fit routine. Return
#' plots comparing the estimated and true parameters, and a basic pass/fail test
#'
test_fit_gpd_mixture_B<-function(test_case=3){

    if(test_case == 1){ 
        ## Hypothetical example. Requires very large sample size to accurately
        ## estimate 'u' and 'xi' [e.g. 5e+05, but not 500 or 5000] 
        ## Probably this is because the gamma distribution is well approximated
        ## in its tail with the GPD distribution having a small shape parameter,
        ## so the fit is still 'good' if the fitted GPD threshold creeps into
        ## the 'true' gamma distribution
        gshape = 1
        gscale = 1
        xi = 0.05
        u = qgamma(0.9, gshape, 1/gscale)
        sample_size = 500
    }
    
    if(test_case == 2){
        # Modification of test_case 1
        # Question: Can we more easily estimate u and xi if the true GPD shape
        # parameter is not similar to the 'tail' gpd-shape parameter that the
        # tail of a gamma distribution should converge to?
        gshape = 1
        gscale = 1
        xi = 0.25
        u = qgamma(0.9, gshape, 1/gscale)
        sample_size = 500
    }

    if(test_case == 3){
        # Example like one of my hsig fits
        # Can be fit reasonably well with 500 data points
        gshape= 0.84
        gscale= 1.01
        u = 1.27
        xi = -0.22
        sample_size = 500
    }

    if(test_case == 4){
        # Example like a duration fit
        gshape= 0.6
        gscale=  48
        u = 11
        xi = -0.05
        sample_size = 500
    }

    # Simulate random data and fit it
    parfun<-function(i){
        
        test_data = qgammagpdcon(runif(sample_size), gshape=gshape, 
            gscale=gscale, u=u, xi=xi)
        
        test_fit = fit_gpd_mixture(data=test_data, ntrial_u=20)

        return(test_fit)
    }

    RNGkind("L'Ecuyer-CMRG")

    # Some of the optimizations will fail, reflecting that these models
    # are hard to fit. The code expects this -- suppress those warnings.
    outputs = suppressWarnings(
        parallel::mclapply(as.list(1:120), parfun, mc.preschedule=FALSE, 
            mc.cores=12, mc.silent=TRUE)
        )

    # Extract all fitted parameters
    fitted_par = matrix(as.numeric(unlist(
        lapply(outputs, 
            f<-function(x){ 
                if(class(x) == 'try-error'){ 
                    return(rep(NA,4))
                }else{
                    return(x$fit_optim$par)
                }
            }
        ))),
        ncol=4, byrow=TRUE)

    png(paste0('test_gpd_mixture_fit_', test_case, '.png'), width=14, height=10,
        res=200, units='in')

    # Add plot of predicted/true parameters
    par(mfrow=c(3,2))
    for(i in 1:4){
        hist(fitted_par[,i], n=20, main=paste0('Par ', i))
        true_par = c(gshape, gscale, u, xi)[i] 
        abline(v=true_par, col='red')
        if( (quantile(fitted_par[,i], 0.95, na.rm=TRUE) > true_par) &
            (quantile(fitted_par[,i], 0.05, na.rm=TRUE) < true_par)){
            print('PASS')
        }else{
            print('FAIL')
        }
    }

    # It is useful later to have a quantile function that
    # we can pass the 'true' parameters to. This will be 
    # in 'qf' inside a successful fit
    for(i in 1:length(outputs)){
        if(class(outputs[[i]]) != 'try-error'){
            full_quantile_function = outputs[[i]]$qf
            break
        }
    }

    # Add plots of predicted / true quantiles
    quantiles_to_plot = c(0.99, 0.999)
    for(i in 1:length(quantiles_to_plot)){
        quantile_to_plot = quantiles_to_plot[i]
        value_q = unlist(as.numeric(lapply(outputs,
            f<-function(x){
                if(class(x) == 'try-error'){
                    return(NA)
                }else{
                    return(x$qfun(quantile_to_plot))
                }
            }
        )))

        hist(value_q, main=paste0(quantile_to_plot, ' quantile'), n=20)
        true_quantile = full_quantile_function(quantile_to_plot, 
            c(gshape, gscale, u, xi))
        abline(v=true_quantile, col='red')

        # Make enveloping the quantile into a test
        if( (quantile(value_q, 0.95, na.rm=TRUE) > true_quantile ) &
            (quantile(value_q, 0.05, na.rm=TRUE) < true_quantile )){
            print('PASS')
        }else{
            print('FAIL')
        }
    }

    dev.off()

}


test_all<-function(){

    set.seed(1)
    
    print('')
    print('Test gpd_mixture fit on a case with realistic parameters ...')
    test_fit_gpd_mixture_B(test_case=3)
    print('')

    print('')
    print('Test gpd_mixture fit on another case with realistic parameters ...')
    test_fit_gpd_mixture_B(test_case=4)
    print('')
    
    
    print('')
    print('Test of mcmc uncertainty for gpd mixture fit (also makes a png plot) ...')
    test_mcmc_gpd_mixture()
    print()

}

test_all()
