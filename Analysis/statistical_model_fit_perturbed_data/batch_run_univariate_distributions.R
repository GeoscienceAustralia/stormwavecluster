#
# Run multiple analyses using perturbed data
#

# Copy the code from statistical_model_fit here
system('cp ../statistical_model_fit/statistical_model_univariate_distributions.Rmd .')

# Function to run a single model with randomly perturbed data, in a new shell.
# The actual model run corrresponds to the 'ith' model from knit_univariate_distributions.R
run_random_model<-function(i){
    command_to_run = paste0('Rscript knit_univariate_distributions.R --break_ties ', i)
    system(command_to_run)
}

# Do the above 100 times
nsim = 100
library(parallel)
mclapply(as.list(1:nsim), run_random_model, mc.cores=detectCores())


# Save an object containing summary statistics [makes later analysis convenient]

all_ud = Sys.glob('Rimages/session_univariate_distributions_TRUE_*.Rdata')

source('get_Rimage_data_univariate_distributions.R', local=TRUE)

# Read all images
store_var_list = lapply(as.list(all_ud), get_Rimage_data)

saveRDS(store_var_list, 'univariate_runs_summary_statistics.RDS')
