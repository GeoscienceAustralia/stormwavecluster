R interface to tpx0.72
----------------------

This code is a basic interface to the tpx0.72 program. 

To run it, you need to have tpx0.72 installed.

Installation
-------------
The tpx0.72 program must be installed and compiled to run this package.

Currently (24/05/2016), tpx0.72 can be obtained here: http://volkov.oce.orst.edu/tides/global.html

Once tpx0.72 is installed on your system, the variable `.OPTS_directory` in the
file [OTPS_directory_name.R](BASIC_R_INTERFACE/OTPS_directory_name.R) must be changed to
correspond to the OPTS directory of the tpx072 installation. Otherwise, the code
will not be able to find tpx072, and will throw an error.

Usage
-----

See the script [tidal_computations.R](BASIC_R_INTERFACE/tidal_computations.R) for an example of usage.

Essentially, the user has to provide a site name, site coordinates (lon,lat),
the times at which predictions are desired (with appropriate timezone
information). The script takes care of the rest. 

Testing
--------

We have tested the code by comparison with a measured gauge, but cannot provide
the test data currently as it is not open source. However, the scripts to run
the test (which make many plots comparing observations and measurements) are in
BASIC_R_INTERFACE/test.
