## Code for generalised weibul distribution

#dweibull(x, shape, scale = 1, log = FALSE)
#pweibull(q, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
#qweibull(p, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
#rweibull(n, shape, scale = 1)

dwbl<-function(x, shape, scale=1, location = 0, log=FALSE){
    dweibull(x - location, shape, scale, log)
}

qwbl<-function(p, shape, scale = 1, location = 0, lower.tail = TRUE, log.p = FALSE){
    qweibull(p, shape, scale, lower.tail, log.p) + location
}

pwbl<-function(q, shape, scale = 1, location=0, lower.tail = TRUE, log.p = FALSE){
    pweibull(q - location, shape, scale, lower.tail, log.p)
}

rwbl<-function(n, shape, scale = 1, location=0){
    rweibull(n, shape, scale) + location
}


test_wbl<-function(){

    x=seq(0,1000,len=100)
    p=seq(0,1,len=101)
    #library(fExtremes) # Should be identical to this
    
    xd = dwbl(x, shape=0.2, location=200,scale=100)
    xq = qwbl(p, shape=0.2, location=200,scale=100)
    xp = pwbl(x, shape=0.2, location=200,scale=100)

    xd_compare = dweibull(x-200, shape=0.2, scale=100)

    stopifnot(max(abs(xd-xd_compare)) < 1.0e-12)
    print('PASS')
       
    xp_compare = pweibull(x - 200, shape=0.2, scale=100) 
    stopifnot(max(abs(xp-xp_compare)) < 1.0e-12)
    print('PASS')
   

    xq_compare = qweibull(p, shape=0.2, scale=100) + 200 
    stopifnot(isTRUE(all.equal(xq, xq_compare)))
    print('PASS')


    # Inverse quantile of a random sample should be uniformly distributed
    random_sam = rwbl(1000, shape=0.2, scale=100, location=200)
    random_p = pwbl(random_sam, shape=0.2, scale=100, location=200)
    # Avoid warnings about ties, which can happen around very very small random values.
    stopifnot(suppressWarnings(ks.test(random_p, punif)$p.value) > 1.0e-04)
    print('PASS')

    # Check that qwbl is inverse of pwbl
    p2 = pwbl(xq, shape=0.2, location=200, scale=100)
    stopifnot(isTRUE(all(abs(p2 - p) <= (1.0e-05 * p))))
    print('PASS')

}
