# Look at coverage of:
#   1) Nonparametric bootstrap, and
#   2) Bayesian intervals with uniform priors
# See related code 'gpd_bootstrap_coverage.R' for tests of parametric bootstrap and L-moments

fit_gpd_eva<-function(data){

    fit = eva::gpdFit(data, threshold=0, npp=1)

    # Return as 'mle' for compatibility with ismev::gpd.fit
    fit$mle = as.numeric(fit$par.ests)

    return(fit)
}


#' Compute nonparametric bootstrap based CI's for a return level
gpd_Rl_CI_npbootstrap<-function(data, level, conf, Nboot){
    
    data_fit = fit_gpd_eva(data)

    quantile_store = rep(NA, Nboot)
    for(i in 1:Nboot){
        random_data = sample(data, size=length(data), replace=TRUE)
        random_data_fit = fit_gpd_eva(random_data)
        quantile_store[i] = eva::qgpd(level, loc=0, scale=random_data_fit$mle[1], 
            shape=random_data_fit$mle[2])
    }

    CI = quantile(quantile_store, c((1-conf)/2, 1 - (1-conf)/2), type=6)

    return(CI)
}


gpd_local = new.env()
source('gpd.R', local=gpd_local)

#' Compute bayesian credible intervals with uniform priors for a return level
gpd_Rl_CI_bayes<-function(data, level, conf){

    ll_fun<-function(par){
        # Finite scale
        if(par[1] <= 0){
            return(-Inf)
        }else{
            sum(gpd_local$dgpd(data, location=0, scale=par[1], shape=par[2], log=TRUE))
        }
    }

    gpd_mcmc = MCMCpack::MCMCmetrop1R(ll_fun, theta.init=c(1, 0), burnin=1000, mcmc=10000)

    CI = quantile(gpd_local$qgpd(rep(level, length(gpd_mcmc[,1])), 
        location=0, scale=gpd_mcmc[,1], shape=gpd_mcmc[,2]), 
        probs = c((1-conf)/2, 1 - (1-conf)/2), type=6)

    return(CI)
}

#' 
gpd_CIs<-function(data, level, conf, Nboot=1000){
    CI_bootstrap = gpd_Rl_CI_npbootstrap(data, level, conf, Nboot)
    CI_bayes = gpd_Rl_CI_bayes(data, level, conf)
    return(list(boot=CI_bootstrap, bayes=CI_bayes))
}

#' Compute confidence intervals for random data.
#' Here 'i' is not used, but lets us run the function with mclapply
gpd_CIs_random_data<-function(i, level=0.99, conf=0.95, Nboot=1000, sam_size=100, 
    shape_par=0.1){

    random_data = eva::rgpd(sam_size, shape=shape_par, scale=1, loc=0)

    random_CIs = gpd_CIs(random_data, level, conf, Nboot)

    return(random_CIs)
}


#' Run a coverage test (takes quite some time even in parallel)
coverage_test<-function(sam_size, shape_par, N, N_boot, myquant){

    gpd_local = new.env()
    source('gpd.R', local=gpd_local)
    quantile_true = gpd_local$qgpd(myquant, shape=shape_par)

    parallel_fun<-function(i){
        output = try(gpd_CIs_random_data(i, level=myquant, conf=0.95, Nboot=N_boot, 
            sam_size=sam_size, shape_par=shape_par))

        return(output)
    }

    library(parallel)
    RNGkind("L'Ecuyer-CMRG")
    all_gpd_tests = mclapply(as.list(1:N), parallel_fun, mc.cores=12)

    coverage_boot = (unlist(lapply(all_gpd_tests, 
        f<-function(x){ 
            if('boot' %in% names(x)){
                out = (x$boot[1]<quantile_true)*(x$boot[2]>quantile_true)
            }else{
                out = NA
            }
            return(out)
            }
        )))
    mean_coverage_boot = mean(coverage_boot, na.rm=TRUE)

    coverage_bayes = (unlist(lapply(all_gpd_tests, 
        f<-function(x){ 
            if('bayes' %in% names(x)){
                out = (x$bayes[1]<quantile_true)*(x$bayes[2]>quantile_true)
            }else{
                out = NA
            }
            return(out)
            }
        )))
    mean_coverage_bayes = mean(coverage_bayes, na.rm=TRUE)

    # What about the upper tail performance?

    coverage_upper_boot = (unlist(lapply(all_gpd_tests, 
        f<-function(x){ 
            if('boot' %in% names(x)){
                out = (x$boot[2]>quantile_true)
            }else{
                out = NA
            }
            return(out)
            }
        )))
    mean_coverage_upper_boot = mean(coverage_upper_boot, na.rm=TRUE)

    coverage_upper_bayes = (unlist(lapply(all_gpd_tests, 
        f<-function(x){ 
            if('bayes' %in% names(x)){
                out = (x$bayes[2]>quantile_true)
            }else{
                out = NA
            }
            return(out)
            }
        )))
    mean_coverage_upper_bayes = mean(coverage_upper_bayes, na.rm=TRUE)

    return(environment())
}



coverage_list = list()
for(shape_par in c(-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4)){

    sam_size = 150 # 
    annual_event_rate = sam_size/30 # 30 years of data
    shape_par = shape_par
    N = 1000
    N_boot = 1000
    myquant = 1 - 1/(100 * annual_event_rate)

    coverage_list[[as.character(shape_par)]] = 
        coverage_test(sam_size, shape_par, N, N_boot, myquant)

}

# Store the results
save.image('gpd_CI_coverage2.Rdata')


# Write out results
tmp = rep(NA, length(coverage_list))
# FIXME: After running this I realised I didn't fix the names in the next line of code (only)
# The effect is that the output table has 4 columns of NA's at the start, followed by 4 columns
# of data that I want. The names are all good, so no real problem, it's just unclean.
output_data = data.frame(shape=tmp, bootstrap_coverage=tmp, 
    bootstrap_upper_coverage=tmp, likelihood_coverage=tmp,
    likelihood_upper_coverage=tmp)
for(i in 1:length(coverage_list)){
    output_data$shape[i] = names(coverage_list)[i]
    output_data$npbootstrap_coverage[i] = coverage_list[[i]]$mean_coverage_boot
    output_data$npbootstrap_upper_coverage[i] = coverage_list[[i]]$mean_coverage_upper_boot
    output_data$bayes_coverage[i] = coverage_list[[i]]$mean_coverage_bayes
    output_data$bayes_upper_coverage[i] = coverage_list[[i]]$mean_coverage_upper_bayes
    }

# Round
output_data$npbootstrap_coverage = round(output_data$npbootstrap_coverage, 2)
output_data$npbootstrap_upper_coverage = round(output_data$npbootstrap_upper_coverage, 2)
output_data$bayes_coverage = round(output_data$bayes_coverage, 2)
output_data$bayes_upper_coverage = round(output_data$bayes_upper_coverage, 2)

write.table(output_data, 'gpd_coverage_experiments2.txt', sep=" & ", 
    row.names=FALSE, quote=FALSE)
