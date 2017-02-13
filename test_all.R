# Run all the R/*/test_* test codes

test_codes = c(Sys.glob('R/*/test_*.R'), Sys.glob('R/*/*/*/test_*.R'))

for(test_code in test_codes){

    cat('\n')
    cat('##################################################################\n')
    print(paste0('Running ', test_code))

    # Run the test in its own environment to avoid interactions of variables   
    test_env = new.env() 
    source(test_code, chdir=TRUE, local=test_env)
}
