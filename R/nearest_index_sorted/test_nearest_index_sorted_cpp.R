#######################################################################################
#
# Run a simple test
#
#######################################################################################
source('nearest_index_sorted_cpp.R', local=TRUE)

test_nearest_index_sorted_cpp<-function(){
    x = c(-10, 1, 2,3,4,7)

    # Typical case
    stopifnot(nearest_index_sorted_cpp(x, 1.3, check_is_sorted = 1) == 2)
    print('PASS')

    # Lower bound
    stopifnot(nearest_index_sorted_cpp(x, -20, check_is_sorted = 1) == 1)
    print('PASS')
    
    # Upper bound
    stopifnot(nearest_index_sorted_cpp(x, 20, check_is_sorted = 1) == 6)
    print('PASS')

    # Vectorized case
    stopifnot(all( nearest_index_sorted_cpp(x, c(-20, 1.3, 2.6, 20), check_is_sorted = 1) == c(1, 2, 4, 6)))
    print('PASS')

    # Check we fail when we should
    y = c(x, -20)
    z = try(nearest_index_sorted_cpp(y, 0, check_is_sorted=1), silent=TRUE)
    stopifnot(class(z) == 'try-error')
    print('PASS')

    
    # Case with a small vector to search through
    x2 = c(-10, 10)
    z = nearest_index_sorted_cpp(x2, c(-100, -10, 0.001, 10, 100), check_is_sorted=1)
    stopifnot(all(z == c(1, 1, 2, 2, 2)))

}

test_nearest_index_sorted_cpp()

