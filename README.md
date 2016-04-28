stormwavestat
-------------

This is a collection of R scripts used for our statistical modelling of coastal
storm waves.

It is currently in development and we do not encourage its usage unless you are
collaborating with the developer. 


**Structure**

The subfolders are:

./Analysis - 'driver' scripts used to perform the high level analysis

./Data - location for data. Not all of the data underlying our analysis is open source, so not all dependencies can currently be provided.

./R - self-contained R scripts used for the analysis. These include:

    * nhpoisp -- Code for fitting non-homogeneous poisson processes

    * nearest_index_sorted -- Fast Rcpp code for doing a nearest-neighbour search on a sorted numeric vector

    * wave_dispersion -- Code to solve the airy wave dispersion relation

    * evmix_fit -- Convenience wrappers for fitting extreme value mixture models

Information on using and testing the above codes is provided within their directories.


**Advice on adapting the analysis to another site**

*Note: Currently the Analysis code has not been pushed to the repository, so these comments may not apply*

The Analysis code is not expected to be applied to another site without
modification. This is because it includes choices which may not generalise well to
other locations (e.g. which probability distributions best fit a site), and
hard-coded site specific details (e.g. time-zones, assumptions about structure
of source data, etc). 

It is provided here to be useful as:
    * a source of examples using routines in ./R, which are more generic
    * an example analysis which can be adapted to another site by an experienced user
    * a demonstration of using R for modelling of storm wave statistics


**Bugs and contributions**

Bugs relating to code inside ./R can be raised on the 'issues' page of the
github site. Consider making a github pull request with any fixes. However, if 
you would like to make a major contribution, you should discuss this with the package
maintainer first (gareth.davies.ga.code@gmail.com) to ensure it will be accepted. 

