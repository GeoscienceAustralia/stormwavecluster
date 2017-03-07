library(VineCopula)
#library(copula)
#library(vines)

source('../preprocessing/data_utilities.R', local=TRUE, chdir=TRUE)

#' Make a univariate distribution conditional on a seasonal variable 
#'
#'
#' @param event_statistics data.frame containing the event summary statistics
#' @param var character name of variable in event_statistics which we will make
#' conditional on season
#' @param q_raw function giving the (unconditional) quantiles of 'var', e.g.
#' q_raw(0.7) should give the quantile value that is not exceeded by 70% of the
#' population of 'var'. q_raw must accept vectorised input arguments
#' @param p_raw function giving the (unconditional) inverse quantiles of 'var', e.g.
#' p_raw(2.5) should give the fraction of the population of 'var' below  2.5.
#' p_raw must accept vectorised input arguments.
#' @param test_pqfun Logical. If TRUE, then run a test to check that the derived 
#' conditional quantile and conditional inverse quantile functions really are the
#' inverse of each other. 
#' @param plot_quantiles Logical. If TRUE, then make a plot of conditional/unconditional
#' quantiles at different times of the year.
#' @param startyear character Name of the variable in event_statistics which contains
#' the event time of year (as a decimal year)
#' @return the function environment. The most important
#' variables contained in this environment are 'qfun' and 'pfun' (the quantile
#' and inverse quantile functions, which allow conditional variables to be provided).
#'
make_fit_conditional_on_season<-function(
    event_statistics, 
    var = 'hsig', 
    q_raw=NULL, 
    p_raw=NULL, 
    test_pqfun=TRUE, 
    plot_quantiles=TRUE,
    startyear = 'startyear'){

    # Bring key function arguments into the function environment (since we
    # return the function environment)
    var = var
    startyear = startyear
    # 'Raw' quantile and inverse quantile functions
    q_raw = q_raw
    p_raw = p_raw

    # Remove event_statisics rows which are missing data for var
    is_na_var = which(is.na(event_statistics[[var]]))
    if(length(is_na_var) > 0){
        event_statistics = event_statistics[-is_na_var,]
    }

    # Best value of the phase in 'season'
    season_phi_value = optimize(
        f<-function(x){ 
            cor(event_statistics[[var]], 
                cos(2*pi*(event_statistics[[startyear]] - x)), 
                method='s')
        }, 
        interval=c(-1,1))

    # Season function optimized for var
    seasonfun<-function(t) cos(2*pi*(t - season_phi_value$minimum))

    #
    # Copula relating 'season' and var
    #
    lx = length(event_statistics[[var]])
    copula_data = cbind(rank(event_statistics[[var]])/(lx+1),
        rank(seasonfun(event_statistics[[startyear]]))/(lx+1) )

    # The type of copula is chosen using BiCopSelect
    var_season_copula = BiCopSelect(copula_data[,1], copula_data[,2])

    #
    # Distribution of 'season' using season_phi_value appropriate for var
    # This is an approximation from the data -- but should be ok for our
    # purposes
    #
    pseason = approxfun( 
        c(-1, sort(seasonfun(event_statistics[[startyear]])), 1),
        (0:(lx+1))/(lx+1)
    ) 

    # Quantile function
    qfun<-function(p, conditional_variables=NULL){
        # If needed, adjust p based on the time of year
        if(!is.null(conditional_variables)){
            if(length(p) != length(conditional_variables[[startyear]])){
                stop('qfun inputs must be of same length') 
            }
            season_p = pseason(seasonfun(conditional_variables[[startyear]]))
            #new_p = vines::hinverse(var_season_copula, p, season_p)
            #new_p = VineCopula::BiCopHinv(p, season_p, var_season_copula)[[2]] 
            new_p = BiCopHinv2(p, season_p, var_season_copula)
        }else{
            new_p = p
        }

        # Get quantile value from fitted distribution
        qvar_vals = q_raw(new_p)

        return(qvar_vals)  
    }

    # Inverse quantile function
    pfun<-function(q, conditional_variables=NULL){

        # Get inverse cdf value from fitted distribution
        p_vals = p_raw(q)

        # If needed, adjust p based on the time of year
        if(!is.null(conditional_variables)){
            if(length(q) != length(conditional_variables[[startyear]])){
                stop('pfun inputs must be of same length') 
            }
            season_p = pseason(seasonfun(conditional_variables[[startyear]]))
            #new_p = vines::h(var_season_copula, p_vals, season_p)
            #new_p = VineCopula::BiCopHfunc(pvar_vals, season_p, var_season_copula)[[2]]
            new_p = BiCopHfunc2(p_vals, season_p, var_season_copula)
        }else{
            new_p = p_vals
        }

        return(new_p)
    }

    # Random number generator
    rfun<-function(n, conditional_variables=NULL){
        random_p = runif(n)
        qfun(random_p, conditional_variables=conditional_variables)
    }


    # Code to test that pfun and qfun work
    test_qfun_pfun<-function(){
        # Check that qfun is the inverse of pfun, for a range of p and
        # startyear values, including extremes of each

        p = c(0.0, 1.0e-06, 0.01, 0.3, 0.5, 0.7, 0.9, 0.99, 1.-1e-06, 1)
        startyear_vals = seq(0, 1, len=20)

        pall = expand.grid(p, startyear_vals)

        conditional_var = list()
        conditional_var[[startyear]] = pall[,2]
        q1 = qfun(pall[,1], conditional_variables=conditional_var)
        p1 = pfun(q1, conditional_variables=conditional_var)

        if(!isTRUE(all.equal(pall[,1], p1, tol=5*sqrt(.Machine$double.eps)))){
            print(cbind(p1, pall[,1], q1))
            stop()
        }
        print('Conditional p/q functions passed test: ')
        print('  (Check plots to see if quantiles are ok)')
    }

    # Plot to check that fitted/data quantiles are similar 
    diagnostic_plot<-function(){
    
        # Number large enough so empirical quantiles are close to exact
        NN = 100000 
        model_conditional_var = list(
            startyear=rep(event_statistics[[startyear]], 
                length.out=NN))
        conditional_quantiles = qfun(runif(NN), 
            conditional_variables=model_conditional_var)
        unconditional_quantiles = qfun(runif(NN))
        #conditional_quantiles = qfun(pfun(event_statistics[[var]]), 
        #    conditional_variables=model_conditional_var)
        #unconditional_quantiles = qfun(pfun(event_statistics[[var]])

        # Split the year into categories and do qqplot for each
        year_boundaries = (1:3)/ 3
        year_boundary_labels = paste0(c('First', 'Middle', 'Last'), ' third of year')

        data_cat = rep(1, length(event_statistics[[startyear]]))
        model_cat = rep(1, NN)

        dec<-function(x) x - trunc(x)

        for(yb in year_boundaries){
            data_cat = data_cat + (dec(event_statistics[[startyear]]) > yb)
            model_cat = model_cat + (dec(model_conditional_var[[startyear]]) > yb)
        }

        # Setup the plot    
        layout(matrix(c(1:6), ncol=3, byrow=TRUE))

        for(i in 1:length(year_boundaries)){

            # Only make log-log plot when there are non-positive data values
            if(any(event_statistics[[var]][data_cat==i] <= 0) | 
                any(conditional_quantiles[data_cat==i] <= 0)){
                logplot_par = ''
            }else{
                logplot_par = 'xy'
            }

            # This qq plot has better treatment of extreme quantiles when
            # datasets are of very different size
            qqplot3(event_statistics[[var]][data_cat==i], 
                conditional_quantiles[data_cat==i], 
                xlab='Data', ylab='Model', 
                main=paste0('Non-stationary model (', var, ')\n', 
                    year_boundary_labels[i], sep=""), 
                cex.main=1.5,
                log=logplot_par)
            abline(0, 1, col='red'); grid()
        }

        for(i in 1:length(year_boundaries)){

            # Only make log-log plot when there are non-positive data values
            if(any(event_statistics[[var]][data_cat==i] <= 0) | 
                any(conditional_quantiles[data_cat==i] <= 0)){
                logplot_par = ''
            }else{
                logplot_par = 'xy'
            }
            
            # This qq plot has better treatment of extreme quantiles when
            # datasets are of very different size
            qqplot3(event_statistics[[var]][data_cat==i], 
                unconditional_quantiles[data_cat==i], 
                xlab='Data', ylab='Model', 
                main=paste0('Stationary model (', var, ')\n', 
                    year_boundary_labels[i], sep=""),
                cex.main=1.5,
                log=logplot_par)
            abline(0, 1, col='red'); grid()
        }

    }


    # Run tests when we make it
    if(test_pqfun) test_qfun_pfun()
    if(plot_quantiles) diagnostic_plot()

    # Return the function environment
    return(environment())
}

#hsig_fit_conditional = make_fit_conditional_on_season(
#    event_statistics,
#    var='hsig', 
#    q_raw=hsig_mixture_fit$qfun, 
#    p_raw=hsig_mixture_fit$pfun,
#    startyear = 'startyear')
