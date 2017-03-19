# Simulate the sampling error in the hsig calculation, assuming wave heights
# have a rayleigh distribution with a given number of waves each hour.

library(VGAM) # Package to get Rayleigh distribution

# Rayleigh parameter -- this is equal to the [ distribution_mean / (sqrt(pi/2)) ]
sigmas = seq(1,7,len=100)

# Compute the (approximate) number of waves measured during each hour. The
# sampling variability in hsig will vary inversely with the square root of
# this, so will not be strongly sensitive to moderate changes in the assumed
# wave period.
nominal_wave_period = 11
nwaves = 34 * 60/ nominal_wave_period # The buoys record a 34 min burst #

# Use monte-carlo sampling to compute the sampling variability expected
# of hsig measurement under these conditions.
nreps = 10000

# Storage matrices
store_hsig = matrix(NA, nrow=nreps, ncol=length(sigmas))
store_hmax = store_hsig

# Loop over all 'sigma' values. Each of these will correspond to a single 'true' hsig,
# but measurements of 'hsig' will vary because they are only based on 'nwaves/3' waves.
for(i in 1:length(sigmas)){

    # Loop over each monte-carlo iterate
    for(j in 1:nreps){

        # Generate random heights for 'nwaves' random waves
        rwaves = rrayleigh(nwaves, scale=sigmas[i])

        # Sort, and find the mean of the top third [= hsig], and the max [ = hmax]
        rwaves_sorted = sort(rwaves, decreasing=TRUE)
        top_third = 1:(nwaves/3)
        store_hsig[j,i] = mean(rwaves_sorted[top_third])
        store_hmax[j,i] = rwaves_sorted[1]
    }
}

# For each sigma, compute the associated mean hsig [which will be effectively
# constant]
mean_hsig = colMeans(store_hsig)
# For each sigma, compute the standard deviation of the computed hsig -- this
# gives the expected 'sampling error'. It is proportional to mean_hsig
sd_mean_hsig = apply(store_hsig, 2, sd)

# Note that for a given sigma, the hsig will be very close to normally distributed
# thanks to the central limit theorem. 
summary(sd_mean_hsig / mean_hsig)
