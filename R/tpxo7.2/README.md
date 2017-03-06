R interface to tpxo7.2
----------------------

This code is a basic interface to the tpxo7.2 program. 

To run it, you need to have tpxo7.2 installed.

Installation
-------------
The tpxo7.2 program must be installed and compiled to run this package.

Currently (24/05/2016), tpxo7.2 can be obtained here: http://volkov.oce.orst.edu/tides/global.html

Once tpxo7.2 is installed on your system, the variable `.OPTS_directory` in the
file [OTPS_directory_name.R](BASIC_R_INTERFACE/OTPS_directory_name.R) must be changed to
correspond to the OPTS directory of the tpxo7.2 installation. Otherwise, the code
will not be able to find tpxo7.2, and will throw an error.

Usage
-----

See the script [tidal_computations.R](BASIC_R_INTERFACE/tidal_computations.R) for an example of usage.

Essentially, the user has to provide a site name, site coordinates (lon,lat),
the times at which predictions are desired (with appropriate timezone
information). The script takes care of the rest. 

Testing
--------

A test program is in
[BASIC_R_INTERFACE/test/test_tides.R](BASIC_R_INTERFACE/test/test_tides.R). It
compares the tidal predictions to measured data at a site, and checks that the
residual standard deviation is sufficiently small. It also makes a plot. To run
the test program, open R in the same directory as the [test_tides.R](BASIC_R_INTERFACE/test/test_tides.R)
script, and then run:

```r
    source('test_tides.R')
```

The observational data used for this test was kindly provided by Manly
Hydraulics Laboratory.
