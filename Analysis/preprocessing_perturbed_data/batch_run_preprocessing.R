#
# Run multiple analyses using perturbed data
#

# Copy the code from statistical_model_fit here
system('cp ../preprocessing/*.Rmd ../preprocessing/*.R .')

# Function to run a single model with randomly perturbed data, in a new shell.
run_random_model<-function(i){
    system(paste0('Rscript knit_preprocessing.R --break_ties ', i))
    system(paste0('Rscript knit_extract_storm_events.R --break_ties ', i))
}

# Do the above 100 times
nsim = 100
library(parallel)
mclapply(as.list(1:nsim), run_random_model, mc.cores=detectCores())


