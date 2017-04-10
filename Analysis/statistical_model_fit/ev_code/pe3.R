##############################################################################

#' Pearson type 3 distribution density function
#'
#' Supports both positive and negative scale parameters
#'
#' @param x data
#' @param shape shape parameter
#' @param scale scale parameter, can be negative
#' @param location threshold or location parameter
#' @param log TRUE/FALSE return log density?
#' @return density value, or log density if log=TRUE
#'
dpe3<-function(x, shape, scale, location=0, log=FALSE){
    fx = dgamma( (x-location)/scale, shape=shape, scale=1, log=TRUE) - log(abs(scale)) 
    if(log == FALSE){ fx = exp(fx) }
    return(fx)
}

#' Pearson type 3 cumulative distribution function
#'
#' Supports both positive and negative scale parameters
#'
#' @param q quantile values
#' @param shape shape parameter
#' @param scale scale parameter, can be negative
#' @param location threshold or location parameter
#' @param lower.tail TRUE/FALSE return lower tail or upper tail area
#' @param log.p TRUE/FALSE return log of p?
#' @return CFD value (=area under pdf with domain < q, if lower.tail=TRUE)
#'
ppe3<-function(q, shape, scale, location=0,lower.tail=TRUE, log.p=FALSE){
    # CDF for pearson type 3 distribution
    #
    # Idea: sign(scale)*[x-location] has a gamma distribution with parameters shape,abs(scale)
    #
    # Idea: pgamma( q, shape, scale) = pgamma(scale*q, shape, 1)

    # Ensure shape, scale,location have length=length(x)
    scale = q*0+scale
    shape = q*0+shape
    location = q*0+location

    # Find where scale >0
    mm = (scale>=0)
    log_p = numeric(length(scale))

    # Cases with positive scale parameter
    if(sum(mm)>0){
        indz=which(mm)
        # Compute log(p)
        log_p[indz]=pgamma( (q[indz] - location[indz]), shape=shape[indz], scale=scale[indz], 
            lower.tail=lower.tail, log.p=TRUE)
    
        if(log.p==FALSE) log_p[indz]=exp(log_p)
    }

    # Cases with negative scale parameter
    if(sum(mm) < length(scale)){
        indz = which(!mm)
        # Need to convert p to 1-p, due to the reversal of sign caused by the
        # sign(scale) operation, which reverses the limits in the cdf
        # integral
        # Since we work directly with log(p), use expm1 to get 1-p
        #
        # NOTE: expm1 = exp(x) - 1, and is accurate for x close to 0

        log_p[indz]= 1- pgamma(sign(scale[indz])*(q[indz]-location[indz]), 
                               shape=shape[indz], scale=abs(scale[indz]),
                               lower.tail=lower.tail,log.p=FALSE)
        if(log.p==TRUE) log_p[indz] = log(log_p)
    }

    return(log_p)
}

#' Pearson type 3 inverse cumulative distribution function
#'
#' Supports both positive and negative scale parameters
#'
#' @param p probabilities 
#' @param shape shape parameter
#' @param scale scale parameter, can be negative
#' @param location threshold or location parameter
#' @param lower.tail TRUE/FALSE return lower tail or upper tail area
#' @param log.p TRUE/FALSE assume the input p is actually log(p)?
#' @return quantile value 
#'
qpe3<-function(p, shape, scale, location=0, lower.tail=TRUE, log.p=FALSE){
    # Inverse CDF for log-pearson type 3 distribution
    #
    # Idea: sign(scale)*[x-location] has a gamma distribution with parameters shape,abs(scale)
    #
    # Because of the sign(scale) operation, we must replace p with 1-p for
    # negative scale
    #
    
    # Ensure shape, scale,location have length=length(p)
    scale = p*0 + scale
    shape = p*0 + shape
    location = p*0 + location 

    # Find where scale >0
    mm = (scale >= 0)
  
    qua_nolog = numeric(length(scale)) 

    if(sum(mm) > 0){ 
        # Cases with positive scale
        indz = which(mm)
        qua_nolog[indz] = qgamma(p[indz], shape=shape[indz], scale=abs(scale[indz]), lower.tail=lower.tail, log.p=log.p)
    }

    if(sum(mm) < length(scale)){
        # Cases with negative scale
        # WARNING: Can be subject to round-off error -- could potentially be improved.
        indz = which(!mm)
        if(log.p == FALSE){
            qua_nolog[indz] = qgamma(p[indz], shape=shape[indz], scale=abs(scale[indz]), lower.tail=!lower.tail, log.p=log.p)
        }else{
            # Take exp of probabilities 
            qua_nolog[indz] = qgamma(exp(p[indz]), shape=shape[indz], scale=abs(scale[indz]), lower.tail=!lower.tail, log.p=FALSE)
        }

    }
    
    qua_nolog = (qua_nolog/sign(scale)+location)

    return(qua_nolog)
}

#' Pearson type 3 distribution random numbers
#'
#' Supports both positive and negative scale parameters
#'
#' @param n number of values
#' @param shape shape parameter
#' @param scale scale parameter, can be negative
#' @param location threshold or location parameter
#' @return n random numbers
#'
rpe3<-function(n, shape, scale, location=0){
    # Random numbers for log-pearson type 3 distribution
    #
    # Idea: sign(scale)*[x-location] has a gamma distribution with parameters shape,abs(scale)
    #
    myrand = (rgamma(n, shape=shape, scale = abs(scale))/sign(scale) + location)
    return(myrand)
}


test_pe3<-function(){
    # Test against lmomco if available
    require(lmomco)
  
    # Preliminary parameters
    shape = 0.5
    scale = 9.5
    location = 13
    
    data0 = rgamma(1e+06, shape=shape, scale=abs(scale))*sign(scale) + location

    # Lmoment parameters for testing
    data0_lmom_fit = parpe3(lmom.ub(data0))

    # scale, shape, location for testing (similar but not exactly = preliminary parameters)
    shape = as.numeric( 4/(data0_lmom_fit$para['gamma'])**2 )
    scale = as.numeric( 0.5 * data0_lmom_fit$para['sigma'] * abs(data0_lmom_fit$para['gamma']) )
    location = as.numeric( data0_lmom_fit$para['mu'] - 2*data0_lmom_fit$para['sigma']/data0_lmom_fit$para['gamma'] )

    # quantile
    p = c(0.01, 0.1, 0.5, 0.9, 0.99)
    checkQ = all.equal( quape3(p, data0_lmom_fit) , qpe3(p, shape, scale, location) )
    if(!isTRUE(checkQ)){
        stop('quantile function failed')
    }else{
        print('pass')
    }

    # Probability
    q = c(1, 5, 10, 100)
    checkp = all.equal( as.numeric(cdfpe3(q, data0_lmom_fit)) , ppe3(q, shape, scale, location) )
    if(!isTRUE(checkp)){
        stop('quantile function failed')
    }else{
        print('pass')
    }

    # Density
    x = seq(1,500,len=500)
    lmom_vals = pdfpe3(x, data0_lmom_fit) 
    our_vals = dpe3(x, shape, scale, location)
    # The lmoments package can give density values of NaN when they should be zero
    # Fix this here
    lmom_vals[is.na(lmom_vals)] = 0
    checkd = all.equal(lmom_vals, our_vals)
    if(!isTRUE(checkd)){
        stop('density function failed')
    }else{
        print('pass')
    }

    # Random numbers
    #print('Making random numbers ...')
    N = 10000000
    m1 = rpe3(N, shape, scale, location)
    m1_dens = density(m1, n=N/4, from = min(m1), to=max(m1))
    plotshift = min(m1_dens$x)+1

    # These should agree with each other except perhaps for edge effects
    #plot(m1_dens$x+plotshift, m1_dens$y, t='l', log='xy', 
    #    main='Empirical density from rpe3 vs analytical density, positive scale')
    #points(m1_dens$x+plotshift, dpe3(m1_dens$x, shape, scale, location), t='l', col=2)
    #legend('bottomleft', c('Empirical Density from random deviates','Density'), 
    #    lwd=c(1,1), col=c(1,2))

    p = c(0.01, 0.1, 0.4, 0.7, 0.9, 0.99)
    random_quants = quantile(m1, p=p)
    pe3_quants = qpe3(p, shape, scale, location)

    err = (random_quants - pe3_quants)/pe3_quants
    if(max(abs(err)) < 0.01){
        print('pass')
    }else{
        stop('random number generation failed')
    }

    ################################################################################
    #
    # Check for accuracy / consistency with negative scale
    shape = 187
    scale = -9.6
    location = 13

    # Check ppe3 is the integral of dpe3
    probs = c(0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
    xx = qpe3(probs, shape, scale, location) 
    pxx = ppe3(xx,shape,scale,location)
    passTest = rep(0,length(xx))
    for(i in 1:length(pxx)){
        myInt = integrate(dpe3, -Inf, xx[i], shape=shape, scale=scale, location=location, rel.tol=1.0e-12)
        passTest[i] = isTRUE(all.equal(myInt$value, pxx[i]))
        if(!passTest[i]){
            stop('CDF fails with negative scale parameter')
        }
    }
    print('pass')
    

    # Check that exp( ppe3(..., log.p=TRUE) ) = ppe3(...)
    pxxLog=ppe3(xx,shape,scale,location,log.p=TRUE)
    if(!isTRUE(all.equal(exp(pxxLog),pxx))){
        stop('CDF with log.p=TRUE not consistent with log.p=FALSE')
    }else{
        print('pass')
    }
    
    # Check that 1-ppe3(..., lower.tail=FALSE) ) = ppe3(...)
    pxxUpper=ppe3(xx,shape,scale,location,lower.tail=FALSE)
    if(!isTRUE(all.equal(1-pxxUpper,pxx))){
        stop('CDF with lower.tail=FALSE not consistent with default')
    }else{
        print('pass')
    }
    
    # Check that 1-exp(ppe3(..., lower.tail=FALSE,log.p=TRUE) )) = ppe3(...)
    pxxUpperLog = ppe3(xx,shape,scale,location,lower.tail=FALSE,log.p=TRUE)
    if(!isTRUE(all.equal(1-exp(pxxUpperLog),pxx))){
        stop('CDF with lower.tail=FALSE and log.p=TRUE not consistent with default')
    }else{
        print('pass')
    }

    # Use qpe3 to get quantile values at pxx, assuming it has passed
    quant_pxx=qpe3(pxxLog, shape,scale,location,log.p=TRUE)
    quant_pxx2=qpe3(pxx, shape,scale,location,log.p=FALSE)
    quant_pxx3=qpe3(pxxUpper, shape,scale,location,lower.tail=FALSE,log.p=FALSE)
    quant_pxx4=qpe3(pxxUpperLog, shape,scale,location,lower.tail=FALSE,log.p=TRUE)
    if((!isTRUE(all.equal(quant_pxx,quant_pxx2))) |
       (!isTRUE(all.equal(quant_pxx,quant_pxx3))) |
       (!isTRUE(all.equal(quant_pxx,quant_pxx4))) 
        ){
        stop('Inv CDF not consistent using some combinations of log.p=TRUE/FALSE and lower.tail=TRUE/FALSE')
    }else{
        print('pass')
    }

    if(!isTRUE(all.equal(quant_pxx,xx, tol=1.0e-03))){
        stop('Quantile function fails with negative scale parameter')
    }else{
        print('pass')
    }

    #print('Making random numbers ...')
    N=1000000
    m1=rpe3(N, shape, scale, location)
    
    p = c(0.01, 0.1, 0.4, 0.7, 0.9, 0.99)
    random_quants = quantile(m1, p=p)
    pe3_quants = qpe3(p, shape, scale, location)
    
    err = (random_quants - pe3_quants)/pe3_quants
    if(max(abs(err)) < 0.01){
        print('pass')
    }else{
        stop('random number generation failed')
    }

}
