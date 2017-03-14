#
# Run multiple analyses using perturbed data
#

# Copy the code from statistical_model_fit here
system('cp ../statistical_model_fit/statistical_model_storm_timings.Rmd .')

# Function to run a single model with randomly perturbed data, in a new shell.
run_random_model<-function(i){
    system('Rscript knit_storm_timings.R --break_ties ')
}

# Do the above 100 times
# This takes a few hours on my 12 core desktop.
nsim = 100
library(parallel)
mclapply(as.list(1:nsim), run_random_model, mc.cores=detectCores())

#
# Check the results
#

# Read each R session into its own environment
all_storm_timing_sessions = Sys.glob('Rimages/*.Rdata')
list_envs = list()
for(i in 1:length(all_storm_timing_sessions)){
    session_file = all_storm_timing_sessions[i]
    list_envs[[session_file]] = new.env()
    load(session_file, envir = list_envs[[session_file]])
}

# Check the 'best fit' rate equations in each R session
all_rate_eqns = unlist(lapply(list_envs, f<-function(x) x$best_nhp_model$rate_equation))

# If the chosen model is unaffected by data perturbations, then only one
# equation should appear here.
print(table(all_rate_eqns))

if(length(unique(all_rate_eqns)) > 1){
    msg = paste0('More than one rate model identified with perturbed data.\n',
        ' The code below must be changed to deal with this case')
    stop(msg)
}

# Check variations in the model parameters.
# The following only works if all_rate_eqns are identical
all_rate_par= matrix(
    unlist(lapply(list_envs, f<-function(x) x$best_nhp_model$par)), 
    ncol=4, byrow=TRUE)
# Coefficient of variation of estimates. Seems to be very small (e.g. 1/1000)
all_rate_cov = apply(all_rate_par, 2, sd)/apply(all_rate_par, 2, mean)
