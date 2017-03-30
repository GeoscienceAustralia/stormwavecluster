stormwavecluster
----------------

This is a collection of R scripts and data used for our statistical modelling of coastal
storm waves.

**Structure**

The sub-folders are:

[./Analysis](Analysis) - 'driver' scripts used to perform the high level analysis, and an associated vignette.

[./Data](Data) - location for data used in the analysis.

[./R](R) - self-contained R scripts used for the analysis. These include:

* [nhpoisp](R/nhpoisp) -- Code for fitting non-homogeneous Poisson processes

* [wave_dispersion](R/wave_dispersion) -- Code to solve the airy wave dispersion relation

* [nearest_index_sorted](R/nearest_index_sorted) -- Fast Rcpp code for doing a
nearest-neighbour search on a sorted numeric vector. This is a simple task, but
we were unable to find a similarly fast and memory efficient way of doing this
in native R.

* [evmix_fit](R/evmix_fit) -- Convenience wrappers for fitting extreme value mixture models

* [tpxo7.2](R/tpxo7.2) -- Convenience interface for the tpxo7.2 tidal prediction software. The latter must be installed separately.

Information on using and testing the above codes is provided within their directories.

**Dependencies**

The tpxo7.2 package should be installed to use the tpxo7.2 interface.  However,
the remainder of the code and analysis will run fine even if this is not done,
so if you just want to look at the statistical analysis, this can be skipped.

The code relies on many different R packages, which can be obtained within R
using the command:
```r

install.packages(c( 
    'knitr', 'devtools', 'Rcpp',  'forecast',
    'CDVine', 'evmix', 'MCMCpack', 'numDeriv', 
    'optimx', 'logspline', 'ismev', 'VineCopula',
    'Matching', 'TwoCop')) 

```
At the time of writing there is a bug (of significance for our analysis) in the
version of `VineCopula` available on CRAN. This version was downloaded with the
above command. The bug is fixed in a version of VineCopula on github, which can
be obtained within R using the command:
```r
devtools::install_github("tnagler/VineCopula")
```
The above requires that you are in an environment where you can build R
packages. This is usually automatically supported on linux and mac, but may require
pre-installation of Rtools on windows, see:
https://cran.r-project.org/bin/windows/Rtools/ .

**Testing**

To run all the tests for the 'R' folder, start R and execute the
[test_all.R](test_all.R) script in the current directory:

```r
    source('test_all.R')
```

This may take awhile to run (currently the extreme value mixture model tests
take tens of minutes on my machine, as do the tests of non-homogeneous-Poisson process fitting).

If you are only interested in a single piece of code in the R folder, then look
there for a corresponding test script.

**Advice on adapting the analysis to another site**

The Analysis code is not expected to be applied to another site without
modification. This is because it includes choices which may not generalise well to
other locations (e.g. which probability distributions best fit a site), and
hard-coded site specific details (e.g. time-zones, assumptions about structure
of source data, etc.). Also, it relies heavily on non-linear optimization techniques,
which may need tweaking to converge with other datasets.

Therefore, we strongly advise against its use as a black-box code. Nonetheless,
it is provided here to be useful as:

* a source of examples using routines in ./R, which are more generic
* an example analysis which can be adapted to another site by an experienced user
* to make our work more reproducible and transparent


**Bugs, maintainence and contributions**

Bugs relating can be raised on the 'issues' page of the github site. Consider
making a github pull request with any fixes. However, if you would like to make
a major contribution, you should discuss this with the package maintainer first
(gareth.davies.ga.code@gmail.com) to ensure it will be accepted. 


**Acknowledgements**

Development of this code was supported through the [Bushfire and Natural
Hazards Cooperative Research Centre](https://www.bnhcrc.com.au/), as part of
the project *Resilience to clustered disaster events on the coast: Storm Surge*.

This project would be impractical without the efforts of the community who
develop and maintain the R software, and its large ecosystem of packages.
