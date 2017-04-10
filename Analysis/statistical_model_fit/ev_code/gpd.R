## All this code is slightly adapted from fExtremes
## so I can distribute gpd code without the
## many dependencies of fExtremes

dgpd<-function (x, location = 0, scale = 1, shape = 0, log = FALSE) 
{
    stopifnot(min(scale) > 0)
    d <- (x - location)/scale

    nn <- max(c(length(d), length(location), length(scale), length(shape)))
    scale <- rep(scale, length.out = nn)
    shape <- rep(shape, length.out = nn)
    index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)

    if(max(index) == 1){

        d[index] = ifelse(shape[index] == 0,
            log(1/scale[index]) - d[index],
            log(1/scale[index]) - (1/shape[index] + 1) * log(1 + shape[index] * d[index]) 
            )
    }

    if(min(index) == 0){
        d[!index] <- -Inf
    }

    if (!log) 
        d <- exp(d)
    #attr(d, "control") = data.frame(location = location[1], scale = scale[1], 
    #    shape = shape[1], log = log, row.names = "")
    return(d)
}


qgpd<-function (p, location = 0, scale = 1, shape = 0, lower.tail = TRUE){
    stopifnot(min(scale) > 0)
    #stopifnot(length(shape) == 1)
    stopifnot(min(p, na.rm = TRUE) >= 0)
    stopifnot(max(p, na.rm = TRUE) <= 1)
    
    nn <- max(c(length(p), length(scale), length(shape), length(location)))
    scale <- rep(scale, length.out = nn)
    shape <- rep(shape, length.out = nn)
    location <- rep(location, length.out = nn)

    if (lower.tail) 
        p <- 1 - p

    index = which(shape == 0)

    q = p * 0

    if (length(index) > 0) {
        q[index] = location[index] - scale[index] * log(p[index])
    }

    index = which(shape != 0)
    if(length(index) > 0){
        q[index] = location[index] + scale[index] * (p[index]^(-shape[index]) - 1)/shape[index]
    }
    #attr(q, "control") = data.frame(location = location[1], scale = scale[1], 
    #    shape = shape[1], lower.tail = lower.tail, row.names = "")
    return(q)
}

pgpd<-function (q, location = 0, scale = 1, shape = 0, lower.tail = TRUE){
    stopifnot(min(scale) > 0)
    #stopifnot(length(shape) == 1)

    nn <- max(c(length(q), length(location), length(scale), length(shape)))
    scale <- rep(scale, length.out = nn)
    shape <- rep(shape, length.out = nn)
    location <- rep(location, length.out = nn)

    q <- pmax(q - location, 0)/scale

    p = q * 0

    index = which(shape == 0)
    if(length(index) > 0){
        p[index] <- 1 - exp(-q[index])
    }

    index = which(shape != 0)
    if(length(index) > 0){
        p[index] <- pmax(1 + shape[index] * q[index], 0)
        p[index] <- 1 - p[index]^(-1/shape[index])
    }

    if (!lower.tail) p <- 1 - p

    return(p)
}


rgpd<-function (n, location = 0, scale = 1, shape = 0) 
{
    stopifnot(min(scale) > 0)

    if(max(c(length(location), length(scale), length(shape))) > n){
        stop('n is less than the length of a parameter vector')
    }

    scale = rep(scale, length.out=n)
    shape = rep(shape, length.out=n)
    location = rep(location, length.out=n)

    zero_shp = which(shape == 0)
    non_zero_shp = which(shape != 0)

    r = rep(NA, n)

    if(length(zero_shp) > 0){
        r[zero_shp]  = location[zero_shp] + scale[zero_shp] * rexp(n)
    }

    if (length(non_zero_shp) > 0) {
        r[non_zero_shp] = location[non_zero_shp] + scale[non_zero_shp] * (runif(n)^(-shape[non_zero_shp]) - 1)/shape[non_zero_shp]
    }
    #attr(r, "control") = data.frame(location = location[1], scale = scale[1], 
    #    shape = shape[1], row.names = "")
    return(r)
}

#############################
test_gpd_code<-function(){

    x=seq(0,1000,len=100)
    p=seq(0,1,len=101)
    #library(fExtremes) # Should be identical to this
    
    xd = dgpd(x, shape=0.2, location=200,scale=100)
    xq = qgpd(p, shape=0.2, location=200,scale=100)
    xp = pgpd(x, shape=0.2, location=200,scale=100)
    set.seed(1) 
    xr=rgpd(100, shape=0.2, location=200,scale=100)


    library(fExtremes)

    yd=fExtremes::dgpd(x,xi=0.2,mu=200,beta=100)
    
    stopifnot(all(xd==yd))
    print('PASS - same as dgpd fExtremes')

    yq=fExtremes::qgpd(p,xi=0.2,mu=200,beta=100)
    stopifnot(all(xq==yq))
    print('PASS - same as qgpd fExtremes')

    yp=fExtremes::pgpd(x,xi=0.2,mu=200,beta=100)
    stopifnot(all(xp==yp))
    print('PASS - same as pgpd fExtremes')
   
    set.seed(1) 
    yr=fExtremes::rgpd(100,xi=0.2,mu=200,beta=100)
    stopifnot(all(xr==yr))
    print('PASS - same as rgpd fExtremes')

}
