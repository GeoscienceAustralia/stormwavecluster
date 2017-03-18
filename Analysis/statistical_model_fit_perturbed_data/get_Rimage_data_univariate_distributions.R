
#
# The session images are too large load all at once.
# So get the important variables here. 
#
get_Rimage_data<-function(all_ud_i){
    library(evmix)
    library(logspline)
    library(VineCopula)
    
    random_env = new.env()
    load(all_ud_i, envir=random_env)

    # Store key variables from random_env in a list, which we will later put in
    # store_var_list
    output = list()

    # Get break_ties_with_jitter
    output$break_ties_with_jitter = random_env$break_ties_with_jitter

    # Get the perturbed event statistics
    output$event_statistics = random_env$event_statistics

    # Hsig info -- parameters, and (Bayesian) quantiles of 1/100 AEP, and
    # seasonal copula type, and seasonal phase
    output$hsig_mixture_fit_par = random_env$hsig_mixture_fit$fit_optim$par
    output$hsig_aep100_quantiles = quantile(random_env$hsig_mixture_fit$ari_100_chains[[1]], 
        p=c(0.025, 0.5, 0.975))
    output$hsig_season_phi = random_env$hsig_fit_conditional$season_phi_value$minimum
    output$hsig_season_copula = random_env$hsig_fit_conditional$var_season_copula

    output$hsig_mixture_fit_bayes_rates = random_env$hsig_mixture_fit$desired_rates * 
        random_env$hsig_mixture_fit$annual_event_rate
    output$hsig_mixture_fit_bayes_median_q = random_env$hsig_mixture_fit$desired_median_q
    output$hsig_mixture_fit_bayes_lower_q = random_env$hsig_mixture_fit$desired_lower_q
    output$hsig_mixture_fit_bayes_upper_q = random_env$hsig_mixture_fit$desired_upper_q
    output$hsig_mixture_fit_ml = random_env$hsig_mixture_fit$desired_best_q

    # Duration info -- parameters, and (Bayesian) quantiles of 1/100 AEP, and
    # seasonal copula type, and seasonal phase
    output$duration_mixture_fit_par = random_env$duration_mixture_fit$fit_optim$par
    output$duration_aep100_quantiles = quantile(random_env$duration_mixture_fit$ari_100_chains[[1]], 
        p=c(0.025, 0.5, 0.975))
    output$duration_season_phi = random_env$duration_fit_conditional$season_phi_value$minimum
    output$duration_season_copula = random_env$duration_fit_conditional$var_season_copula

    output$duration_mixture_fit_bayes_rates = random_env$duration_mixture_fit$desired_rates * 
        random_env$duration_mixture_fit$annual_event_rate
    output$duration_mixture_fit_bayes_median_q = random_env$duration_mixture_fit$desired_median_q
    output$duration_mixture_fit_bayes_lower_q = random_env$duration_mixture_fit$desired_lower_q
    output$duration_mixture_fit_bayes_upper_q = random_env$duration_mixture_fit$desired_upper_q
    output$duration_mixture_fit_ml = random_env$duration_mixture_fit$desired_best_q

    # TideResid info -- parameters, and (Bayesian) quantiles of 1/100 AEP, and
    # seasonal copula type, and seasonal phase
    output$tideResid_mixture_fit_par = random_env$tideResid_mixture_fit$fit_optim$par
    output$tideResid_aep100_quantiles = quantile(random_env$tideResid_mixture_fit$ari_100_chains[[1]], 
        p=c(0.025, 0.5, 0.975))
    output$tideResid_season_phi = random_env$tideResid_fit_conditional$season_phi_value$minimum
    output$tideResid_season_copula = random_env$tideResid_fit_conditional$var_season_copula

    output$tideResid_mixture_fit_bayes_rates = random_env$tideResid_mixture_fit$desired_rates * 
        random_env$tideResid_mixture_fit$annual_event_rate
    output$tideResid_mixture_fit_bayes_median_q = random_env$tideResid_mixture_fit$desired_median_q
    output$tideResid_mixture_fit_bayes_lower_q = random_env$tideResid_mixture_fit$desired_lower_q
    output$tideResid_mixture_fit_bayes_upper_q = random_env$tideResid_mixture_fit$desired_upper_q
    output$tideResid_mixture_fit_ml = random_env$tideResid_mixture_fit$desired_best_q

    # steepness
    output$steepness_season_phi = random_env$steepness_fit_conditional$season_phi_value$minimum
    output$steepness_season_copula = random_env$steepness_fit_conditional$var_season_copula 
    output$steepness_quantiles_seq = random_env$steepness_fit_conditional$q_raw(seq(0,1,len=200))

    # direction
    output$dir_season_phi = random_env$dir_fit_conditional$season_phi_value$minimum
    output$dir_season_copula = random_env$dir_fit_conditional$vargivensoiA_season_copula 
    output$dir_soiA_copula = random_env$dir_fit_conditional$var_soiA_copula
    output$dir_quantiles_seq = random_env$dir_fit_conditional$q_raw(seq(0,1,len=200))

    # Store the result and move on
    #store_var_list[[all_ud[i]]] = output
    return(output)
}

