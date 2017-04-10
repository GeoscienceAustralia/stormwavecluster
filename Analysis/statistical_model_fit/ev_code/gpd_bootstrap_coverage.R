# Look at parametric bootstrap CI coverage for
# 1) ML, and
# 2) Lmoments.
#


## #' My function to fit a GPD (cross-check with other methods) 
## fit_gpd_GD<-function(data){
## 
##     # scale = exp(par[1]) to prevent negative scale optimization issues
##     gpd_nll_local<-function(par){ 
##         -sum(gpd_local$dgpd(data, shape=par[2], scale=exp(par[1]), threshold=0,
##             log=TRUE))
##     }
## 
##     fit = optim(c(1, 0), gpd_nll_local, method='BFGS', 
##         control=list(parscale=c(1, 0.1)))
## 
##     # Correct scale. Return as 'mle' for compatibility with ismev::gpd.fit
##     fit$mle = c(exp(fit$par[1]), fit$par[2])
## 
##     return(fit)
## }

#' Fit a GPD using the eva package. This has a nice method for profile
#' likelihood CI's of return levels
fit_gpd_eva<-function(data){

    fit = eva::gpdFit(data, threshold=0, npp=1)

    # Return as 'mle' for compatibility with ismev::gpd.fit
    fit$mle = as.numeric(fit$par.ests)

    return(fit)
}

#' Compute parametric bootstrap based CI's for a return level
gpd_Rl_CI_bootstrap<-function(data, level, conf, Nboot){
    
    data_fit = fit_gpd_eva(data)

    quantile_store = rep(NA, Nboot)
    for(i in 1:Nboot){
        random_data = eva::rgpd(n=length(data), loc=0, scale=data_fit$mle[1], 
            shape=data_fit$mle[2])
        random_data_fit = fit_gpd_eva(random_data)
        quantile_store[i] = eva::qgpd(level, loc=0, scale=random_data_fit$mle[1], 
            shape=random_data_fit$mle[2])
    }

    CI = quantile(quantile_store, c((1-conf)/2, 1 - (1-conf)/2), type=6)

    return(CI)
}

#' Compute profile likelihood based CIs for a return level
gpd_Rl_CI_proflik<-function(data, level, conf){
    data_fit = fit_gpd_eva(data)
    
    level_CI = eva::gpdRl(data_fit, period=1/(1-level), conf=conf, 
        method='profile')

    return(level_CI$CI)
}

#' 
gpd_CIs<-function(data, level, conf, Nboot=1000){
    CI_bootstrap = gpd_Rl_CI_bootstrap(data, level, conf, Nboot)
    CI_proflik = gpd_Rl_CI_proflik(data, level, conf)
    return(list(boot=CI_bootstrap, prof=CI_proflik))
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

    coverage_prof = (unlist(lapply(all_gpd_tests, 
        f<-function(x){ 
            if('prof' %in% names(x)){
                out = (x$prof[1]<quantile_true)*(x$prof[2]>quantile_true)
            }else{
                out = NA
            }
            return(out)
            }
        )))
    mean_coverage_prof = mean(coverage_prof, na.rm=TRUE)

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

    coverage_upper_prof = (unlist(lapply(all_gpd_tests, 
        f<-function(x){ 
            if('prof' %in% names(x)){
                out = (x$prof[2]>quantile_true)
            }else{
                out = NA
            }
            return(out)
            }
        )))
    mean_coverage_upper_prof = mean(coverage_upper_prof, na.rm=TRUE)

    return(environment())
}


##################################################################
#
# GPD CI coverage with maximum likelihood
#
##################################################################

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
save.image('gpd_CI_coverage.Rdata')

stop('Deliberate stop here')

# Write out results
tmp = rep(NA, length(coverage_list))
output_data = data.frame(shape=tmp, bootstrap_coverage=tmp, 
    bootstrap_upper_coverage=tmp, likelihood_coverage=tmp,
    likelihood_upper_coverage=tmp)
for(i in 1:length(coverage_list)){
    output_data$shape[i] = names(coverage_list)[i]
    output_data$bootstrap_coverage[i] = coverage_list[[i]]$mean_coverage_boot
    output_data$bootstrap_upper_coverage[i] = coverage_list[[i]]$mean_coverage_upper_boot
    output_data$likelihood_coverage[i] = coverage_list[[i]]$mean_coverage_prof
    output_data$likelihood_upper_coverage[i] = coverage_list[[i]]$mean_coverage_upper_prof
    }

# Round
output_data$bootstrap_coverage = round(output_data$bootstrap_coverage, 2)
output_data$bootstrap_upper_coverage = round(output_data$bootstrap_upper_coverage, 2)
output_data$likelihood_coverage = round(output_data$likelihood_coverage, 2)
output_data$likelihood_upper_coverage = round(output_data$likelihood_upper_coverage, 2)

write.table(output_data, 'gpd_coverage_experiments.txt', sep=" & ", 
    row.names=FALSE, quote=FALSE)

######################################################################
#
# Here we run a variant of the above test with lmoments, to see if
# we can get the same results as Kysely (2010). 
#
# Summary: For the tests I did we get similar results to Kysely, and my checks
# also suggest the parametric bootstrap has good coverage using L-moment based
# fitting
# 
######################################################################

if(FALSE){
    # Kysley, sample size of 100, shape par of 0.1, 100 year return level
    # from 100 data points
    sam_size = 100 # 
    shape_par = 0.1
    N = 1000
    N_boot = 1000
    myquant = 1 - 1/100 #1 - 1/(150 * 3)
}else{
    # This case has a bounded upper tail. The coverage is still good.
    sam_size = 100 # 
    shape_par = -0.1
    N = 1000
    N_boot = 1000
    myquant = 1 - 1/100 #1 - 1/(150 * 3)
}


gpd_local = new.env()
source('gpd.R', local=gpd_local)
quantile_true = gpd_local$qgpd(myquant, shape=shape_par)

library(lmomco)
gpd_bootstrap_test_lmom<-function(i){
    # Lmomco uses the reverse shape parameter
    theoretical_fit = list(
        type='gpa', 
        para=c(0, 1, -shape_par), 
        zeta=1, 
        source='pargpd')
   
    # Generate random data and fit to it 
    tmp_data = rlmomco(sam_size, theoretical_fit)
    tmp_fit = pargpa(lmom.ub(tmp_data), xi=0)

    quantile_store = qlmomco(myquant, tmp_fit)

    # Do N_boot bootstraps
    inner_quantile_store = rep(NA, N_boot)
    inner_fit_par = matrix(NA, ncol=2, nrow=N_boot)
    for(j in 1:N_boot){
        tmp_data_inner =  rlmomco(sam_size, tmp_fit)
        tmp_fit_inner = pargpa(lmom.ub(tmp_data_inner), xi=0)
        
        inner_fit_par[j,] = tmp_fit_inner$par[2:3]
        inner_quantile_store[j] = qlmomco(myquant, tmp_fit_inner)
    }
    return(environment())

}

library(parallel)
RNGkind("L'Ecuyer-CMRG")
all_gpd_tests_lmom = mclapply(as.list(1:N), gpd_bootstrap_test_lmom, mc.cores=12)

# How often would a confidence interval contain the true value?
all_gpd_quantiles_lmom_upper = unlist(lapply(all_gpd_tests_lmom, 
    f<-function(x) quantile(x$inner_quantile_store, 0.95, type=6)))
all_gpd_quantiles_lmom_lower = unlist(lapply(all_gpd_tests_lmom, 
    f<-function(x) quantile(x$inner_quantile_store, 0.05, type=6)))
mean((all_gpd_quantiles_lmom_upper > quantile_true)*(all_gpd_quantiles_lmom_lower < quantile_true))

#
bootstrap_maxima_below_data_maxima_fraction_lmom = unlist(lapply(
    all_gpd_tests_lmom, 
    f<-function(x){
        # Maximum of GPD is -scale/shape if shape<0 and threshold=0
        # Otherwise maximum is Inf
        # But lmom has reversed shape!
        boot_maxima = x$inner_fit_par[,1]/x$inner_fit_par[,2]
        boot_maxima[x$inner_fit_par[,2] < 0] = Inf
        return(mean(boot_maxima < max(x$tmp_data)))
    }
))

