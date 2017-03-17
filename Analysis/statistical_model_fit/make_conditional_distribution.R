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
    # Restrict to 1 parameter families, which are indexed by these integers
    one_par_copulas = c(1, 3:6, 13:14, 16, 23:24, 26, 33:34, 36) 
    var_season_copula = BiCopSelect(copula_data[,1], copula_data[,2], 
        familyset=one_par_copulas)

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

            # Only make log-log plot when all data are positive
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

#' Make a univariate distribution conditional on the mean annual SOI, and a
#' seasonal variable.
#'
#' Assumes the mean annual SOI is independent of the seasonal variable, and
#' has distribution which is well modelled by a normal distribution. Both hold
#' for our dataset
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
#' @param soiA character Name of the variable in event_statistics which contains the
#' mean annual SOI at the time of the event.
#' @return the function environment. The most important
#' variables contained in this environment are 'qfun' and 'pfun' (the quantile
#' and inverse quantile functions, which allow conditional variables to be provided).
#'
make_fit_conditional_on_soiA_and_season<-function(
    event_statistics, 
    var = 'dir', 
    q_raw=NULL, 
    p_raw=NULL, 
    test_pqfun=TRUE, 
    plot_quantiles=TRUE,
    startyear = 'startyear',
    soiA = 'soiA'){

    # Bring key function arguments into the function environment (since we
    # return the function environment)
    var = var
    startyear = startyear
    soiA = soiA

    # 'Raw' quantile and inverse quantile functions
    q_raw = q_raw
    p_raw = p_raw

    # Remove event_statisics rows which are missing data for var, or where soiA
    # is missing [i.e. incomplete years]
    is_na_var = which(is.na(event_statistics[[var]]) | is.na(event_statistics[[soiA]]))
    if(length(is_na_var) > 0){
        event_statistics = event_statistics[-is_na_var,]
    }


    # Need a soiA percentile function
    # Assume soiA normal distribution (graphically this is a good fit)
    unique_soiA = unique(event_statistics[[soiA]])
    unique_soiA_mean = mean(unique_soiA)
    unique_soiA_sd = sd(unique_soiA)
    psoiA<-function(q) pnorm(q, mean=unique_soiA_mean, sd=unique_soiA_sd)
    qsoiA<-function(p) qnorm(p, mean=unique_soiA_mean, sd=unique_soiA_sd)


    #
    # Copula relating 'soiA' and dir. Breaking ties in soiA at random does not
    # significantly change the copula fit.
    #
    copula_data_soiA = cbind(
        rank(event_statistics[[var]])/
            (length(event_statistics[[var]])+1),
        psoiA(event_statistics[[soiA]])
        #psoiA(jitter(event_statistics[[soiA]]))
        )

    # Restrict to 1 parameter families, which are indexed by these integers
    one_par_copulas = c(1, 3:6, 13:14, 16, 23:24, 26, 33:34, 36) 
    var_soiA_copula = BiCopSelect(copula_data_soiA[,1], copula_data_soiA[,2], 
        familyset=one_par_copulas)

    #
    # Make the 'h' and 'hinverse' functions, giving percentiles of dir adjusted
    # for soiA, and the inverse
    #
    # This gives the conditional non-exceedance probability corresponding to
    # an unconditional non-exceedance probability p (i.e. u_1 value), with the
    # given soiA value (i.e. u_2 value)
    h_soiA<-function(p, soiA){
        #new_p = vines::h(var_soiA_copula, p, psoiA(soiA), eps=1.0e-10)
        new_p = BiCopHfunc2(p, psoiA(soiA), var_soiA_copula)
        return(new_p)
    }
    
    # This gives the unconditional non-exceedance probability (i.e. u_1 value) corresponding to
    # the conditional non-exceedance probability p associated with the given soiA value (i.e. u_2 value)
    hinverse_soiA<-function(p, soiA){
        #new_p = vines::hinverse(var_soiA_copula, p, psoiA(soiA), eps=1.0e-10)
        new_p = BiCopHinv2(p, psoiA(soiA), var_soiA_copula)
        return(new_p)
    }

    ###########################################################################
    #
    # Below we make (var | soiA) conditional on season as well
    #
    ###########################################################################

    # Compute empirical percentiles for dir given soiA
    p_empirical_given_soiA = h_soiA(
        rank(event_statistics[[var]])/(length(event_statistics[[var]])+1),
        event_statistics[[soiA]])

    # Best value of the phase in 'season'
    season_phi_value = optimize(
        f<-function(x){ 
            cor(p_empirical_given_soiA, 
                cos(2*pi*(event_statistics[[startyear]] - x)), 
                method='s')
        }, 
        interval=c(-1,1))

    # Season function optimized for var
    seasonfun<-function(t) cos(2*pi*(t - season_phi_value$minimum))
    #
    # Distribution of 'season' using season_phi_value appropriate for 'var' 
    # Enforce lower bound of -1 (p = 0) and upper bound of 1 (p = 1) implied
    # by a single sinusoidal component in seasonfun
    #
    lx = length(event_statistics[[startyear]])
    pseason = approxfun(
        c(-1, sort(seasonfun(event_statistics[[startyear]])), 1),
        (0:(lx+1))/(lx+1)) 

    # Make copula data with 'var' conditional on soi, to fit with 'season'
    copula_data_season_given_soiA = cbind(
        p_empirical_given_soiA, 
        pseason(seasonfun(event_statistics[[startyear]])))

    # Fit the copula
    vargivensoiA_season_copula = BiCopSelect(copula_data_season_given_soiA[,1], 
        copula_data_season_given_soiA[,2], familyset=one_par_copulas) 

    # h function, which transforms 'soiA corrected p' into 'soiA and season
    # corrected p' 
    h_season_given_soiA<-function(p_adjusted_for_soiA, season){
        #new_p = vines::h(invdir_season_given_soiA_copula, 
        #    p_adjusted_for_soiA, pseason(season), eps=1.0e-10)
        new_p = BiCopHfunc2(p_adjusted_for_soiA, pseason(season), 
            vargivensoiA_season_copula)
        return(new_p)
    }
    
    hinverse_season_given_soiA<-function(p_adjusted_for_soiA, season){
        #new_p = vines::hinverse(invdir_season_given_soiA_copula, 
        #    p_adjusted_for_soiA, pseason(season), eps=1.0e-10)
        new_p = BiCopHinv2(p_adjusted_for_soiA, pseason(season), 
            vargivensoiA_season_copula)
        return(new_p)
    }

    # Quantile function
    qfun<-function(p, conditional_variables=NULL){
        # If needed, adjust p based on the time of year
        if(!is.null(conditional_variables)){

            if(length(p) != length(conditional_variables[[soiA]])){
                stop('qfun inputs must be of same length') 
            }
            if(length(p) != length(conditional_variables[[startyear]])){
                stop('qfun inputs must be of same length') 
            }

            # Adjust for soiA 
            p_given_soiA = hinverse_soiA(p, conditional_variables[[soiA]])
   
            # Adjust for season 
            p_given_soiA_season = hinverse_season_given_soiA(p_given_soiA, 
                seasonfun(conditional_variables[[startyear]]))
            new_p = p_given_soiA_season
        }else{
            new_p = p
        }

        # Get quantile value from fitted distribution
        qfun_vals = q_raw(new_p)

        return(qfun_vals)  
    }


    # Inverse quantile function
    pfun<-function(q, conditional_variables=NULL){

        p_given_soiA_season = p_raw(q)

        # If needed adjust based on time of year and soiA
        if(!is.null(conditional_variables)){
            # Adjust for season
            p_given_soiA = h_season_given_soiA(p_given_soiA_season, 
                seasonfun(conditional_variables[[startyear]]))

            # Adjust for soiA
            new_p = h_soiA(p_given_soiA, conditional_variables[[soiA]])
        }else{
            new_p = p_given_soiA_season
        }

        return(new_p)
    }


    #' Test that pfun, qfun are inverses of each other
    test_qfun_pfun<-function(){
        
        conditional_vars = expand.grid(
            p = c(0.0, 1.0e-06, 0.01, 0.3, 0.5, 0.7, 0.9, 0.99, 1.-1e-06, 1),
            startyear=seq(0,1,len=10), 
            soiA = qsoiA(seq(0.001, 0.999, len=10))
            )
           
        p = conditional_vars$p 
        
        q1 = qfun(p, conditional_variables=conditional_vars)

        p1 = pfun(q1, conditional_variables=conditional_vars)

        if(isTRUE(all.equal(p1, p))){
            print('Conditional p/q functions passed test: Check quantiles are ok')
        }else{
            print(cbind(p, p1, p-p1, q1))
        }
    }

    # For direction, better to focus on ENSO in diagnostic?
    diagnostic_plot<-function(){
    
        # Number large enough so empirical quantiles are close to exact
        NN = 1e+05
        model_conditional_var = list(startyear=rep(event_statistics[[startyear]], 
            length.out=NN), soiA = rep(event_statistics[[soiA]], length.out=NN))
        conditional_quantiles = qfun(runif(NN), 
            conditional_variables=model_conditional_var)
        unconditional_quantiles = qfun(runif(NN))

        #browser()

        # Split soiA up into 3 categories
        soiA_boundaries = c(-5, 5, Inf)
        soiA_labels = c('soiA <= -5', '-5 < soiA <= 5', 'soiA > 5')

        data_cat = rep(1, length(event_statistics[[var]]))
        model_cat = rep(1, NN)
    
        # Value of 1 --> soiA < -5
        # Value of 2 -->  -5 < soiA < 5
        # Value of 3 --> soiA > 5
        for(yb in soiA_boundaries){
            data_cat = data_cat + (event_statistics[[soiA]] > yb)
            model_cat = model_cat + (model_conditional_var[[soiA]] > yb)
        }

        par(mfrow=c(2, length(soiA_boundaries)))

        for(i in 1:length(soiA_boundaries)){
            qqplot3(event_statistics[[var]][data_cat==i], 
                conditional_quantiles[data_cat==i], 
                xlab='Data', ylab='Model', 
                main=paste0('Conditional model (', var, ') \n', soiA_labels[i] ),
                cex.main=1.5)
            abline(0, 1, col='red'); grid()
        }

        for(i in 1:length(soiA_boundaries)){
            
            qqplot3(event_statistics[[var]][data_cat==i], 
                unconditional_quantiles[data_cat==i], 
                xlab='Data', ylab='Model', 
                main=paste0('Unconditional model (', var, ') \n', soiA_labels[i]),
                cex.main=1.5)
            abline(0, 1, col='red'); grid()
        }
    }
    
    if(test_pqfun) test_qfun_pfun()
    if(plot_quantiles) diagnostic_plot()

    return(environment())

}
