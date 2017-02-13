stormwavecluster
----------------

This is a collection of R scripts used for our statistical modelling of coastal
storm waves.

It is currently in development and we do not encourage its use unless you are
collaborating with the developer.


**Structure**

The sub-folders are:

[./Analysis](Analysis) - 'driver' scripts used to perform the high level analysis

[./Data](Data) - location for data. Not all of the data underlying our analysis is open source, so not all dependencies can currently be provided.

[./R](R) - self-contained R scripts used for the analysis. These include:

* [nhpoisp](R/nhpoisp) -- Code for fitting non-homogeneous Poisson processes

* [wave_dispersion](R/wave_dispersion) -- Code to solve the airy wave dispersion relation

* [nearest_index_sorted](R/nearest_index_sorted) -- Fast Rcpp code for doing a
nearest-neighbour search on a sorted numeric vector. We were unable to find a
similarly fast and memory efficient way of doing this in native R.

* [evmix_fit](R/evmix_fit) -- Convenience wrappers for fitting extreme value mixture models

* [tpx07.2](R/tpx07.2) -- Convenience interface for the tpx07.2 tidal prediction software. The latter must be installed separately.

Information on using and testing the above codes is provided within their directories.

**Testing**

To run all the tests for the 'R' folder, start R and execute the
[test_all.R](test_all.R) script in the current directory:

    source('test_all.R')

This may take awhile to run (currently the extreme value mixture model tests
take tens of minutes on my machine).

If you are only interested in a single piece of code in the R folder, then look
their for a corresponding test script.

**Advice on adapting the analysis to another site**

*Note: Currently the Analysis code has not been pushed to the repository, so these comments may not apply*

The Analysis code is not expected to be applied to another site without
modification. This is because it includes choices which may not generalise well to
other locations (e.g. which probability distributions best fit a site), and
hard-coded site specific details (e.g. time-zones, assumptions about structure
of source data, etc.). 

It is provided here to be useful as:

* a source of examples using routines in ./R, which are more generic
* an example analysis which can be adapted to another site by an experienced user
* a demonstration of using R for modelling of storm wave statistics
* to make our work more reproducible and transparent


**Bugs, maintainence and contributions**

Bugs relating to code inside ./R can be raised on the 'issues' page of the
github site. Consider making a github pull request with any fixes. However, if 
you would like to make a major contribution, you should discuss this with the package
maintainer first (gareth.davies.ga.code@gmail.com) to ensure it will be accepted. 

