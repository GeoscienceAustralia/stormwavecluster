**nearest_index_sorted_cpp**
----------------------------

Code to efficiently find the 'nearest neighbour' index in a sorted (monotonic
non-decreasing) numeric vector, using a simple binary search. While this is a simple
operation, a very fast and memory efficient approach seems to be lacking in R.

It takes 3 arguments as input:

    1. A sorted numeric vector X

    2. A real number (or a vector of them) Y. We want to find values in X that are close to Y.

    3. A flag indicating whether we should check that the first argument really is sorted.

The code returns an integer vector I (with the same length as Y) such that X[I]
are the nearest-neighbours of Y in X.

This code relies on the Rcpp package being installed.


**USAGE**

    source('nearest_index_sorted_cpp.R')
   
    # Example sorted data 
    x = c(-10, 1, 2,3,4,7)

    # Find the index in x nearest to 1.3 (should be 2, since x[2] = 1, which is closest to 1.3)

    nearest_index_sorted_cpp(x, 1.3, check_is_sorted=1)
    # [1] 2

    # We can do multiple numbers at once
    nearest_index_sorted_cpp(x, c(1.3, 3.2), check_is_sorted=1)
    # [1] 2 4
   
    # Since check_is_sorted=1 the above command first checks that x is sorted
    # (monotonic non-decreasing).
    # If x were not sorted (i.e. increasing) then the above command would fail

    y = c(x, -20) # y is not increasing
    nearest_index_sorted_cpp(y, 1.3, check_is_sorted=1)
    # Error: First argument must be a monotonic non-decreasing vector

    # If you are doing repeated searches through a large x, it is potentially
    # very inefficient to check that x is sorted every time. Thus you can suppress
    # the check by setting check_is_sorted=0. 

    nearest_index_sorted_cpp(x, 1.3, check_is_sorted=0) # Faster for large x, but dangerous!
    # [1] 2

**BEWARE:** 

If x is not monotonic non-decreasing and check_is_sorted=0, then the code may
produce unreliable results without warning. Thus you should always use
check_is_sorted=1 on the first lookup (to catch unexpected errors), and only
use check_is_sorted=0 on subsequent lookups if speed is an issue.

**TESTS:**

To test the code, open R and run:

    source('test_nearest_index_sorted_cpp.R')

If it repeatedly prints 'PASS' and gives no errors then all is well

