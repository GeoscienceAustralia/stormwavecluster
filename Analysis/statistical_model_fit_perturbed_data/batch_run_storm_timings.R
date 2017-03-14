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


