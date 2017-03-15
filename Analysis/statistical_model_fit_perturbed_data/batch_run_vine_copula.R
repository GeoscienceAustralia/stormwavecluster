#
# Run multiple analyses using perturbed data
#

# Copy the code from statistical_model_fit here
system('cp ../statistical_model_fit/statistical_model_vine_copula.Rmd .')

# Function to run a single model with randomly perturbed data, in a new shell.
# The actual model run corrresponds to the 'ith' model from knit_univariate_distributions.R
run_random_model<-function(i){
    command_to_run = paste0('Rscript knit_vine_copula.R --break_ties ', i)
    system(command_to_run)
}

# Do the above 100 times
nsim = 10
library(parallel)
mclapply(as.list(1:nsim), run_random_model, mc.cores=detectCores())

