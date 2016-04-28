**nhpoisp**
-----------

Code for fitting and simulating non-homogeneous poisson processes. 

Series with gaps between events (e.g. to prevent storm overlap) are supported.
The event rate function can also depend on time and the time since the last event.
Additional covariates can be introduced with some care.

For usage information see the in-line documentation in nhpoisp.R
(this follows the doxygen style), and look at the tests.

**USAGE**

FIXME


**TESTS**

To run the tests, open R and do

    source('test_nhpoisp.R')

This will take awhile (say 30min on my 6 core linux box). It should
print information on each set of tests which are being run, along with
many 'PASS' statements.

**TO-DO**

This code could be converted to a stand-alone R package with relatively little work, 
and is generically useful. 

Note that the interface is flexible but a bit non-standard (particularly the way the rate
function is often passed as a character string). Consider revising.

**SEE ALSO**

There is an R package on CRAN called nhpoisson which is also for fitting
non-homogeneous poisson processes. It was released after the code here was
developed, and seems to have different functionality to the current package
(e.g. when I last looked it did not seem to treat the clustering that we
required for the storm analysis, although that may have changed). It's
interface is very different to that of the current package. Anyway it is
certainly worth investigating if you are working on non-homogeneous poisson
processes.
