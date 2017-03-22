suppressPackageStartupMessages(library(evmix))
# Override some evmix functions with local versions
source('fgammagpd.r', local=TRUE)

library(parallel)
# Here we use mclapply for parallel computation. That doesn't work
# on windows with > 1 core, hence....
if(.Platform$OS.type == 'windows'){
    .DEFAULT_MC_CORES = 1
}else{
    .DEFAULT_MC_CORES = detectCores()
}


#'
#' Copy of evmix::lnormgpdcon, with some parameters checks disabled to enhance
#' speed.
#'
lnormgpdcon<-function (x, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd), 
    xi = 0, phiu = TRUE, log = TRUE){

    
    if(FALSE){
        check.quant(x, allowna = TRUE, allowinf = TRUE)
        check.param(nmean)
        check.param(nsd)
        check.param(u)
        check.param(xi)
        check.phiu(phiu, allowfalse = TRUE)
        check.logic(log)
        if (any(!is.finite(x))) {
            warning("non-finite cases have been removed")
            x = x[is.finite(x)]
        }
        check.quant(x)
    }

    n = length(x)
    check.inputn(c(length(nmean), length(nsd), length(u), length(xi), 
        length(phiu)), allowscalar = TRUE)
    xu = x[which(x > u)]
    nu = length(xu)
    xb = x[which(x <= u)]
    nb = length(xb)
    if (n != nb + nu) {
        stop("total non-finite sample size is not equal to those above threshold and those below or equal to it")
    }
    if ((nsd <= 0) | (u <= min(x)) | (u >= max(x))) {
        l = -Inf
    }
    else {
        if (is.logical(phiu)) {
            pu = pnorm(u, nmean, nsd)
            if (phiu) {
                phiu = 1 - pu
            }
            else {
                phiu = nu/n
            }
        }
        phib = (1 - phiu)/pu
        du = dnorm(u, nmean, nsd)
        sigmau = phiu/(phib * du)
        syu = 1 + xi * (xu - u)/sigmau
        yb = (xb - nmean)/nsd
        if ((min(syu) <= 0) | (sigmau <= 0) | (du < .Machine$double.eps) | 
            (phiu <= 0) | (phiu >= 1) | (pu <= 0) | (pu >= 1)) {
            l = -Inf
        }
        else {
            l = lgpd(xu, u, sigmau, xi, phiu)
            l = l - nb * log(2 * pi * nsd^2)/2 - sum(yb^2)/2 + 
                nb * log(phib)
        }
    }
    if (!log) 
        l = exp(l)
    l
}

#' Copy of evmix::nlnormgpdcon with parameter checks disabled for speed
nlnormgpdcon<-function (pvector, x, phiu = TRUE, finitelik = FALSE){
    np = 4
    if(FALSE){
        check.nparam(pvector, nparam = np)
        check.quant(x, allowna = TRUE, allowinf = TRUE)
        check.phiu(phiu, allowfalse = TRUE)
        check.logic(finitelik)
    }
    nmean = pvector[1]
    nsd = pvector[2]
    u = pvector[3]
    xi = pvector[4]
    nllh = -lnormgpdcon(x, nmean, nsd, u, xi, phiu)
    if (finitelik & is.infinite(nllh)) {
        nllh = sign(nllh) * 1e+06
    }
    nllh
}

#' Copy of evmix::lgpd with parameter checks disabled for speed
lgpd<-function (x, u = 0, sigmau = 1, xi = 0, phiu = 1, log = TRUE){
    if(FALSE){
        check.quant(x, allowna = TRUE, allowinf = TRUE)
        check.param(u)
        check.param(sigmau)
        check.param(xi)
        check.prob(phiu)
        check.logic(log)
        check.inputn(c(length(u), length(sigmau), length(xi), length(phiu)), 
            allowscalar = TRUE)
    }

    if (any(!is.finite(x))) {
        warning("non-finite cases have been removed")
        x[!is.finite(x)] = NA
    }
    xu = x[which(x > u)]
    nu = length(xu)
    yu = (xu - u)/sigmau
    syu = 1 + xi * yu
    if ((min(syu) <= 0) | (sigmau <= 0) | (phiu <= 0) | (phiu > 
        1)) {
        l = -Inf
    }
    else {
        if (abs(xi) < 1e-06) {
            l = -nu * log(sigmau) - sum(yu) + nu * log(phiu)
        }
        else {
            l = -nu * log(sigmau) - (1/xi + 1) * sum(log(syu)) + 
                nu * log(phiu)
        }
    }
    if (!log) 
        l = exp(l)
    l
}

#' Copy of evmix::lgammagpdcon with parameter checks turned off for speed
lgammagpdcon<-function (x, gshape = 1, gscale = 1, u = qgamma(0.9, gshape, 
    1/gscale), xi = 0, phiu = TRUE, log = TRUE){

    if(FALSE){
        check.quant(x, allowna = TRUE, allowinf = TRUE)
        check.param(gshape)
        check.param(gscale)
        check.param(u)
        check.param(xi)
        check.phiu(phiu, allowfalse = TRUE)
        check.logic(log)
    }

    if (any(!is.finite(x))) {
        warning("non-finite cases have been removed")
        x = x[is.finite(x)]
    }
    if (any(x <= 0)) {
        warning("non-positive values have been removed")
        x = x[x > 0]
    }
    check.quant(x)
    n = length(x)

    if(FALSE){
        check.inputn(c(length(gshape), length(gscale), length(u), 
            length(xi), length(phiu)), allowscalar = TRUE)
    }

    xu = x[which(x > u)]
    nu = length(xu)
    xb = x[which(x <= u)]
    nb = length(xb)
    if (n != nb + nu) {
        stop("total non-finite sample size is not equal to those above threshold and those below or equal to it")
    }
    if ((gscale <= 0) | (gshape <= 0) | (u <= 0) | (u <= min(x)) | 
        (u >= max(x))) {
        l = -Inf
    }
    else {
        if (is.logical(phiu)) {
            pu = pgamma(u, gshape, scale = gscale)
            if (phiu) {
                phiu = 1 - pu
            }
            else {
                phiu = nu/n
            }
        }
        phib = (1 - phiu)/pu
        du = dgamma(u, gshape, scale = gscale)
        sigmau = phiu/(phib * du)
        syu = 1 + xi * (xu - u)/sigmau
        if ((min(syu) <= 0) | (sigmau <= 0) | (du < .Machine$double.eps) | 
            (phiu <= 0) | (phiu >= 1) | (pu <= 0) | (pu >= 1)) {
            l = -Inf
        }
        else {
            l = lgpd(xu, u, sigmau, xi, phiu)
            l = l + (gshape - 1) * sum(log(xb)) - sum(xb)/gscale - 
                nb * gshape * log(gscale) - nb * lgamma(gshape) + 
                nb * log(phib)
        }
    }
    if (!log) 
        l = exp(l)
    l
}


#' Copy of nlgammagpdcon without parameter checks for speed
nlgammagpdcon<-function (pvector, x, phiu = TRUE, finitelik = FALSE){
    np = 4
    if(FALSE){
        check.nparam(pvector, nparam = np)
        check.quant(x, allowna = TRUE, allowinf = TRUE)
        check.phiu(phiu, allowfalse = TRUE)
        check.logic(finitelik)
    }
    gshape = pvector[1]
    gscale = pvector[2]
    u = pvector[3]
    xi = pvector[4]
    nllh = -lgammagpdcon(x, gshape, gscale, u, xi, phiu)
    if (finitelik & is.infinite(nllh)) {
        nllh = sign(nllh) * 1e+06
    }
    nllh
}

#' Version of qnormgpd with checks removed (for more speed)
qnormgpd<-function (p, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd), 
    sigmau = nsd, xi = 0, phiu = TRUE, lower.tail = TRUE) 
{

    if(FALSE){
        check.prob(p, allowna = TRUE)
        check.param(nmean, allowvec = TRUE)
        check.posparam(nsd, allowvec = TRUE)
        check.param(u, allowvec = TRUE)
        check.posparam(sigmau, allowvec = TRUE)
        check.param(xi, allowvec = TRUE)
        check.phiu(phiu, allowvec = TRUE)
        check.logic(lower.tail)
    }

    n = check.inputn(c(length(p), length(nmean), length(nsd), 
        length(u), length(sigmau), length(xi), length(phiu)), 
        allowscalar = TRUE)
    if (!lower.tail) 
        p = 1 - p
    p = rep(p, length.out = n)
    nmean = rep(nmean, length.out = n)
    nsd = rep(nsd, length.out = n)
    u = rep(u, length.out = n)
    sigmau = rep(sigmau, length.out = n)
    xi = rep(xi, length.out = n)
    pu = pnorm(u, nmean, nsd)
    if (is.logical(phiu)) {
        phiu = 1 - pu
    }
    else {
        phiu = rep(phiu, length.out = n)
    }
    phib = (1 - phiu)/pu
    q = p
    whichb = which(p <= (1 - phiu))
    nb = length(whichb)
    whichu = which(p > (1 - phiu))
    nu = length(whichu)
    if (nb > 0) 
        q[whichb] = qnorm(p[whichb]/phib[whichb], nmean[whichb], 
            nsd[whichb])
    if (nu > 0) 
        q[whichu] = evmix::qgpd(p[whichu], u[whichu], sigmau[whichu], 
            xi[whichu], phiu[whichu])
    q
}

#' Version of qnormgpdcon with checks removed for speed
qnormgpdcon<-function (p, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd), 
    xi = 0, phiu = TRUE, lower.tail = TRUE) 
{

    if(FALSE){
        check.prob(p, allowna = TRUE)
        check.param(nmean, allowvec = TRUE)
        check.posparam(nsd, allowvec = TRUE)
        check.param(u, allowvec = TRUE)
        check.param(xi, allowvec = TRUE)
        check.phiu(phiu, allowvec = TRUE)
        check.logic(lower.tail)
    }
    n = check.inputn(c(length(p), length(nmean), length(nsd), 
        length(u), length(xi), length(phiu)), allowscalar = TRUE)
    p = rep(p, length.out = n)
    nmean = rep(nmean, length.out = n)
    nsd = rep(nsd, length.out = n)
    u = rep(u, length.out = n)
    xi = rep(xi, length.out = n)
    pu = pnorm(u, nmean, nsd)
    if (is.logical(phiu)) {
        phiu = 1 - pu
    }
    else {
        phiu = rep(phiu, length.out = n)
    }
    phib = (1 - phiu)/pu
    sigmau = phiu/(phib * dnorm(u, nmean, nsd))
    if(FALSE) check.posparam(sigmau, allowvec = TRUE)
    qnormgpd(p, nmean, nsd, u, sigmau, xi, phiu, lower.tail)
}



#' Version of qgammagpd with checks removed (for more speed)
qgammagpd<-function (p, gshape = 1, gscale = 1, u = qgamma(0.9, gshape, 
    1/gscale), sigmau = sqrt(gshape) * gscale, xi = 0, phiu = TRUE, 
    lower.tail = TRUE) 
{

    if(FALSE){
        check.prob(p, allowna = TRUE)
        check.posparam(gshape, allowvec = TRUE)
        check.posparam(gscale, allowvec = TRUE)
        check.posparam(u, allowvec = TRUE)
        check.posparam(sigmau, allowvec = TRUE)
        check.param(xi, allowvec = TRUE)
        check.phiu(phiu, allowvec = TRUE)
        check.logic(lower.tail)
    }
    n = check.inputn(c(length(p), length(gshape), length(gscale), 
        length(u), length(sigmau), length(xi), length(phiu)), 
        allowscalar = TRUE)
    if (!lower.tail) 
        p = 1 - p
    p = rep(p, length.out = n)
    gshape = rep(gshape, length.out = n)
    gscale = rep(gscale, length.out = n)
    u = rep(u, length.out = n)
    sigmau = rep(sigmau, length.out = n)
    xi = rep(xi, length.out = n)
    pu = pgamma(u, gshape, scale = gscale)
    if (is.logical(phiu)) {
        phiu = 1 - pu
    }
    else {
        phiu = rep(phiu, length.out = n)
    }
    phib = (1 - phiu)/pu
    q = p
    whichb = which(p <= (1 - phiu))
    nb = length(whichb)
    whichu = which(p > (1 - phiu))
    nu = length(whichu)
    if (nb > 0) 
        q[whichb] = qgamma(p[whichb]/phib[whichb], gshape[whichb], 
            scale = gscale[whichb])
    if (nu > 0) 
        q[whichu] = evmix::qgpd(p[whichu], u[whichu], sigmau[whichu], 
            xi[whichu], phiu[whichu])
    q
}

#' Version of evmix:pgammagpdcon. For speed some parameter checks are removed
pgammagpdcon<-function (q, gshape = 1, gscale = 1, u = qgamma(0.9, gshape, 
    1/gscale), xi = 0, phiu = TRUE, lower.tail = TRUE){

    # Removing these checks speeds it up significantly
    if(FALSE){
        check.quant(q, allowna = TRUE, allowinf = TRUE)
        check.posparam(gshape, allowvec = TRUE)
        check.posparam(gscale, allowvec = TRUE)
        check.posparam(u, allowvec = TRUE)
        check.param(xi, allowvec = TRUE)
        check.phiu(phiu, allowvec = TRUE)
        check.logic(lower.tail)
    }
    n = check.inputn(c(length(q), length(gshape), length(gscale), 
        length(u), length(xi), length(phiu)), allowscalar = TRUE)
    if (any(is.infinite(q))) 
        warning("infinite quantiles set to NA")
    q[is.infinite(q)] = NA
    q = rep(q, length.out = n)
    gshape = rep(gshape, length.out = n)
    gscale = rep(gscale, length.out = n)
    u = rep(u, length.out = n)
    xi = rep(xi, length.out = n)
    pu = pgamma(u, gshape, scale = gscale)
    if (is.logical(phiu)) {
        phiu = 1 - pu
    }
    else {
        phiu = rep(phiu, length.out = n)
    }
    phib = (1 - phiu)/pu
    sigmau = phiu/(phib * dgamma(u, gshape, scale = gscale))
    if(FALSE) check.posparam(sigmau, allowvec = TRUE)
    pgammagpd(q, gshape, gscale, u, sigmau, xi, phiu, lower.tail)
}


#' This fixes a bug in evmix (the argument 'scale=' is required in the dgamma call)
#' Carl Scarrot identified the fix, it will go into the evmix package sometime
#' in the coming months. Also, for speed some parameter checks are removed
qgammagpdcon<-function(p, gshape = 1, gscale = 1, u = qgamma(0.9, gshape, 
    1/gscale), xi = 0, phiu = TRUE, lower.tail = TRUE){

    # Removing these checks speeds it up significantly
    if(FALSE){
        check.prob(p, allowna = TRUE)
        check.posparam(gshape, allowvec = TRUE)
        check.posparam(gscale, allowvec = TRUE)
        check.posparam(u, allowvec = TRUE)
        check.param(xi, allowvec = TRUE)
        check.phiu(phiu, allowvec = TRUE)
        check.logic(lower.tail)
    }
    n = check.inputn(c(length(p), length(gshape), length(gscale), 
        length(u), length(xi), length(phiu)), allowscalar = TRUE)
    p = rep(p, length.out = n)
    gshape = rep(gshape, length.out = n)
    gscale = rep(gscale, length.out = n)
    u = rep(u, length.out = n)
    xi = rep(xi, length.out = n)
    pu = pgamma(u, gshape, scale = gscale)
    if (is.logical(phiu)) {
        phiu = 1 - pu
    }
    else {
        phiu = rep(phiu, length.out = n)
    }
    phib = (1 - phiu)/pu
    sigmau = phiu/(phib * dgamma(u, gshape, scale=gscale))
    if(FALSE) check.posparam(sigmau, allowvec = TRUE)
    qgammagpd(p, gshape, gscale, u, sigmau, xi, phiu, lower.tail)

}


#'
#' Code to fit gamma-gpd or normal-gpd mixture models where
#' the gpd models the upper tail
#'
#' @param fit_env environment. Most code is ran in this environment and it is
#' returned by the function.
#' @param data Numeric vector. The data to fit
#' @param data_offset Numeric constant. Before fitting this is subtracted from
#' data. Typically required to ensure the lower bound of the gamma-gpd is zero
#' @param starting_par Vector of starting parameters for the model. If not
#' provided a decent guess is made
#' @param gpd_threshold_quantile_range Range of the INITIAL search for the gpd
#' threshold parameter, in terms of percentiles of (data - data_offset). Note this
#' does not constrain the final fit (see constrained_fit_gpd_mixture for that)
#' @param bulk 'gamma' or 'normal'
#' @param phiu TRUE/FALSE. TRUE corresponds to the 'bulk' tail fraction approach
#' of Scarrot (2015), whereas FALSE corresponds to the 'parameterised' tail fraction
#' @param ntrial_u Number of thresholds to trial in initial fit
#' @param ntrial_iterations. Repeatedly apply optim 'ntrial_iterations' times
#' using different methods to make convergence more likely
#' @return a modification of fit_env
#'
fit_gpd_mixture<-function(
    fit_env = new.env(), 
    data = NULL, 
    data_offset = 0, 
    starting_par = NULL,
    gpd_threshold_quantile_range=c(0.05, 0.95),
    bulk='gamma',
    phiu=TRUE,
    continuous=TRUE,
    ntrial_u = 20,
    ntrial_iterations=10,
    verbose=TRUE){

    if(is.null(data)) stop('Must provide data')

    if(!is.logical(phiu)) stop('phiu must be logical.')

    with(fit_env, {

        # 
        data = data
        data_offset = data_offset
        starting_par = starting_par
        gpd_threshold_quantile_range = gpd_threshold_quantile_range
        bulk = bulk
        phiu = phiu
        ntrial_u = ntrial_u
        verbose = verbose
       
        # Gamma must have lower bound of zero
        data_trans = data - data_offset

        # Get a consistent interface for all model families
        # fitter_evmix, nll_fun, qf, pf
        if(bulk == 'gamma' & continuous){
            fitter_evmix = evmix::fgammagpdcon

            # Use the local nllfun (with parameter checks disabled) for speed
            nll_fun<-function(par, x) nlgammagpdcon(pvector=par[1:4], x=x, phiu=phiu)
        
            qf<-function(p, par, ...){
                if(is.null(dim(par))){
                    qgammagpdcon(p, par[1], par[2], par[3], par[4], ...)
                }else{
                    qgammagpdcon(p, par[,1], par[,2], par[,3], par[,4], ...)
                }
            }

            pf<-function(q, par, ...){
                if(is.null(dim(par))){
                    pgammagpdcon(q, par[1], par[2], par[3], par[4], ...)
                }else{
                    pgammagpdcon(q, par[,1], par[,2], par[,3], par[,4], ...)
                }
            }
        
        }else if(bulk == 'normal' & continuous){
            fitter_evmix = evmix::fnormgpdcon

            # Use the local nllfun (with parameter checks disabled) for speed
            nll_fun<-function(par, x) nlnormgpdcon(pvector=par[1:4], x=x, phiu=phiu)

            qf<-function(p, par, ...){
                if(is.null(dim(par))){
                    qnormgpdcon(p, par[1], par[2], par[3], par[4], ...)
                }else{
                    qnormgpdcon(p, par[,1], par[,2], par[,3], par[,4], ...)
                }
            }
            pf<-function(q, par, ...){
                if(is.null(dim(par))){
                    pnormgpdcon(q, par[1], par[2], par[3], par[4], ...)
                }else{
                    pnormgpdcon(q, par[,1], par[,2], par[,3], par[,4], ...)
                }
            }
        
        }else if(bulk == 'gamma' & !continuous){
            fitter_evmix = fgammagpd

            nll_fun<-function(par, x) nlgammagpd(pvector=par[1:5], x=x, phiu=phiu)
        
            qf<-function(p, par, ...){
                if(is.null(dim(par))){
                    qgammagpd(p, par[1], par[2], par[3], par[4], par[5], ...)
                }else{
                    qgammagpd(p, par[,1], par[,2], par[,3], par[,4], par[,5], ...)
                }
            }
            pf<-function(q, par, ...){
                if(is.null(dim(par))){
                    pgammagpd(q, par[1], par[2], par[3], par[4], par[5], ...)
                }else{
                    pgammagpd(q, par[,1], par[,2], par[,3], par[,4], par[,5], ...)
                }
            }
        
        }else if(bulk == 'normal' & !continuous){
            fitter_evmix = evmix::fnormgpd

            nll_fun<-function(par, x) evmix::nlnormgpd(pvector=par[1:5], x=x, phiu=phiu)

            qf<-function(p, par, ...){
                if(is.null(dim(par))){
                    qnormgpd(p, par[1], par[2], par[3], par[4], par[5], ...)
                }else{
                    qnormgpd(p, par[,1], par[,2], par[,3], par[,4], par[,5], ...)
                }
            }
            pf<-function(q, par, ...){
                if(is.null(dim(par))){
                    pnormgpd(q, par[1], par[2], par[3], par[4], par[5], ...)
                }else{
                    pnormgpd(q, par[,1], par[,2], par[,3], par[,4], par[,5], ...)
                }
            }
        }

        # Fit a bulk model
        useq_values = quantile( 
            data_trans, 
            p = seq(gpd_threshold_quantile_range[1], gpd_threshold_quantile_range[2], 
                len=ntrial_u))

        # Initially, use a grid search over a range of threshold values
        fit_evmix = fitter_evmix(x=data_trans, phiu=phiu, useq=useq_values)

        # Repeat optimization to dodge local minima
        for(i in 1:ntrial_iterations){
            init_fit = fit_evmix

            fit_evmix = suppressWarnings(fitter_evmix(x=data_trans, phiu=phiu, pvector=init_fit$mle, 
                method='Nelder-Mead'))
            # Make sure it really did improve things
            if(fit_evmix$nllh > init_fit$nllh) fit_evmix = init_fit
            init_fit = fit_evmix

            fit_evmix = suppressWarnings(fitter_evmix(x=data_trans, phiu=phiu, pvector=init_fit$mle, 
                method='BFGS'))
            # Make sure it really did improve things
            if(fit_evmix$nllh > init_fit$nllh) fit_evmix = init_fit

            if(max(abs(fit_evmix$mle - init_fit$mle))< 1.0e-06) break
        }
        #rm(init_fit)

        if(verbose) print(c('  evmix fit NLLH: ', fit_evmix$nllh))

        if(is.null(starting_par)) starting_par = fit_evmix$mle

        # Check the continuous bulk fit
        fit_optim = optim(starting_par, nll_fun, x=data_trans, method='Nelder-Mead')
        
        # Repeat optimization to dodge local minima
        for(i in 1:ntrial_iterations){
            init_fit = fit_optim

            fit_optim = optim(init_fit$par, nll_fun, x=data_trans, method='Nelder-Mead')
            # Check it did improve things
            if(fit_optim$value > init_fit$value) fit_optim = init_fit

            init_fit = fit_optim

            fit_optim = optim(init_fit$par, nll_fun, x=data_trans, method='BFGS')
            # Check it did improve things
            if(fit_optim$value > init_fit$value) fit_optim = init_fit

            if(max(abs(fit_optim$par - init_fit$par))< 1.0e-06) break
        }
        if(i == ntrial_iterations & verbose) print('Warning: all iterations used')

        if(verbose){
            print(c('  fit_optim NLLH: ', fit_optim$value)) 
            print(c('  Bulk par estimate0: ', fit_evmix$mle))
            print(c('           estimate1: ', fit_optim$par))
            print(c('  Difference: ', fit_evmix$mle - fit_optim$par))
        }

        # fit_optim seems to be 'always' better or equivalent to using evmix's
        # routines directly, if we use the same starting parameters. Ensure it is
        if(fit_optim$value > fit_evmix$nllh + 1.0e-06){
            stop('fit_optim is not converged')
        }

        # For quantile and inverse quantile functions, phiu cannot be false
        if(isTRUE(phiu)){
            phiu2 = phiu
        }else{
            # Scarrot shows this is equal to the mean number of exceedances
            phiu2 = mean(data_trans > fit_optim$par[3]) 
        }

        # Key functions -- use the local (debugged) quantile/inverse quantile
        qfun <-function(p){
            qf(p, fit_optim$par, phiu=phiu2) + data_offset
        }

        pfun <-function(q){
            pf(q - data_offset, fit_optim$par, phiu=phiu2)
        }

        # Check that it seems to be working
        test_qfun_pfun<-function(){
            pseq = c(1e-12, 1e-06, 0.01, 0.1, 0.5, 0.9, 0.99, 1-1e-06, 1-1e-12)
            
            q0 = qfun(pseq)
            p0 = pfun(q0)

            #rel_err = abs(pseq - p0)/pseq

            l = length(pseq)
            
            # Use all.equal to check for round-off, except when we get right near zero
            # Direct use of all.equal is more conservative for the pseq[1] case, and was
            # sometimes making undesirable failures for our application.
            success = (isTRUE(all.equal(p0[2:l], pseq[2:l])) & (abs(pseq[1] - p0[1]) < 2.0e-12)) 

            if(!success){
            #if(any(rel_err > 1.0e-06)){
                print(cbind(pseq, p0, q0, pseq-p0, (pseq-p0)/pseq), digits=12)
                stop('ERROR in quantile functions: (pfun(qfun(p)) != p) to within the desired tolerance')
            }else{
                if(verbose) print('PASS: checked qfun and pfun are inverse functions')
            }
        }

        # Run tests when we make this
        test_qfun_pfun()

    })

    return(fit_env)
}

#' 
#' Constrained maximum likelihood fit for a gpd mixture model
#' 
#' @param fit_env environment. output of fit_gpd_mixture which fits unconstrained gpd
#' mixture models wite ML
#' @param lower_bounds vector of lower bounds on parameters. length = length(fit_env$fit_optim$par)
#' @param upper_bounds vector of upper bounds on parameters. length = length(fit_env$fit_optim$par)
#' @param start_par starting parameters for constrained optimization. If NULL, parameters from
#' previous (unconstrained) optimization in fit_env will be used. If these do not satisfy
#' the constraints there will be problems.
#' @param data data for the constrained fit. If NULL, then existing data in fit_env is used
#' @param data_offset data_offset for the constrained fit. If NULL, then
#' existing data_offset in fit_env is used.
#' @return Nothing but fit_env is modified to include the constrained
#'
constrained_fit_gpd_mixture<-function(fit_env, lower_bounds, upper_bounds, 
    start_par=NULL, data=NULL, data_offset=NULL, ntrial_iterations = 10){

    if(!exists('nll_fun', envir=fit_env)){
        stop('Environment fit_env does not have an nll_fun')
    }

    # Put the function arguments into fit_env

    fit_env$constrained_ML_lower_bounds = lower_bounds
    fit_env$constrained_ML_upper_bounds = upper_bounds
    fit_env$constrained_ML_start_par = start_par

    if(!is.null(data)){
        fit_env$constrained_ML_data = data
    }else{
        fit_env$constrained_ML_data = fit_env$data
    }

    if(!is.null(data_offset)){
        fit_env$constrained_ML_data_offset = data_offset
    }else{
        fit_env$constrained_ML_data_offset = fit_env$data_offset
    }

    fit_env$constrained_ML_data_trans = fit_env$constrained_ML_data - 
        fit_env$constrained_ML_data_offset

    fit_env$constrained_ML_ntrial_iterations = ntrial_iterations
        
    with(fit_env, {

        #browser()

        # Negative log likelihood with constraints
        nll_fun_constrained<-function(par, x){
            if(any(par < constrained_ML_lower_bounds) | 
               any(par > constrained_ML_upper_bounds)){
                return(Inf)
            }else{
                return(nll_fun(par, x))
            }
        }

        # If start par are not provided, use earlier fit optim
        if(is.null(constrained_ML_start_par)){
            constrained_ML_start_par = fit_optim$par
        }

        if(any(constrained_ML_start_par < constrained_ML_lower_bounds) |
           any(constrained_ML_start_par > constrained_ML_upper_bounds)){
            stop('error in constrained_fit_gpd_mixture: Starting parameters do not satisfy constraints')
        }


        # Initial fit
        fit_constrained = optim(constrained_ML_start_par, nll_fun_constrained, 
            x=constrained_ML_data_trans)

        # Repeatedly optimize to make it more likely we converge
        for(i in 1:(2*constrained_ML_ntrial_iterations)){
            fit_init = fit_constrained

            fit_constrained = optim(fit_init$par, nll_fun_constrained, 
                x=constrained_ML_data_trans)

            # Check it really did improve things
            if(fit_constrained$value > fit_init$value) fit_constrained = fit_init

            # Cannot use BFGS
            if(max(abs(fit_constrained$par - fit_init$par)) < 1.0e-06) break

        }

        if(i == 2*constrained_ML_ntrial_iterations){ 
            print('Warning: all iterations used')
        }

        # Make constrained p/q fun
        qfun_constrained<-function(p){
            qf(p, fit_constrained$par, phiu=phiu2) + constrained_ML_data_offset
        }

        pfun_constrained<-function(q){
            pf(q - constrained_ML_data_offset, fit_constrained$par, phiu=phiu2)
        }

    })

}


#'
#' Compute multiple MCMC chains to characterise uncertainty in a fitted
#' gpd mixture model
#'
#' @param fit_env environment output from fit_gpd_mixture
#' @param par_lower_limits numeric vector. Lower bound on uniform prior parameters
#' @param par_upper_limits numeric vector. Upper bound on uniform prior parameters
#' @param mcmc_start_perturbation numeric vector. MCMC chains are initialised
#' from the optimal value of fit_env$fit_optim plus a uniformly distributed
#' perturbation ranging over +- this value.
#' @param mcmc_length The length of each mcmc chain before thinning
#' @param mcmc_thin Only output every mcmc_thin'th result of the chain. While
#' the optimal value is 1, for highly correlated chains higher values can reduce
#' the output size without significantly compromising accuracy.
#' @param mcmc_burnin Length of the burnin for each mcmc chain
#' @param mcmc_nchains number of mcmc chains to run
#' @param mcmc_tune value of tune passed to mcmcpack::MCMCmetrop1R
#' @param mc_cores Number of cores to use. Chains are spread over cores
#' @param annual_event_rate Number of events per year on average (used to compute ari statistics)
#' @return the environment fit_env updated with the mcmc information. Note that if
#' the mcmc search finds a better optimum than found by fit_gpd_mixture,
#' then it modifies fit_optim to reflect this.
#'
mcmc_gpd_mixture<-function(
    fit_env, 
    par_lower_limits = -Inf, 
    par_upper_limits = Inf, 
    mcmc_start_perturbation=0.05, mcmc_length=1.0e+04, mcmc_thin=1,
    mcmc_burnin=1000, mcmc_nchains=12, mcmc_tune=1,
    mc_cores=.DEFAULT_MC_CORES, annual_event_rate=1, verbose=TRUE){

    # Copy input parameters to fit_env
    fit_env$par_lower_limits = par_lower_limits
    fit_env$par_upper_limits = par_upper_limits
    fit_env$mcmc_start_perturbation = mcmc_start_perturbation
    fit_env$mcmc_length = mcmc_length
    fit_env$mcmc_thin = mcmc_thin
    fit_env$mcmc_burnin = mcmc_burnin
    fit_env$mcmc_nchains = mcmc_nchains
    fit_env$mcmc_tune = mcmc_tune
    fit_env$mc_cores = mc_cores
    fit_env$annual_event_rate = annual_event_rate
    fit_env$mcmc_verbose = verbose

    # Run the main code inside fit_env, and return it as an environment
    with(fit_env, {
        
        loglik_fun <-function(par) {
            if(any(par < par_lower_limits) | any(par > par_upper_limits)){
                return(-Inf)
            }else{
                return(-nll_fun(par, x=data_trans))
            }
        }

        # Simulate multiple MCMC chains (as a check on convergence)
        suppressPackageStartupMessages(library(parallel))
        suppressPackageStartupMessages(library(MCMCpack))
        RNGkind("L'Ecuyer-CMRG")
        mcmc_chain_maker<-function(x){

            # If starting values are 'bad' then mcmc will fail
            # Get around this with a 'while' loop
            finished = FALSE
            counter = 0
            while(finished == FALSE){
                start_vals = fit_optim$par * 1.0 + 
                    (runif(length(fit_optim$par))-0.5)*mcmc_start_perturbation

                mcmcsam = try(MCMCpack::MCMCmetrop1R(loglik_fun, start_vals, 
                    mcmc=mcmc_length, thin=mcmc_thin, burnin=mcmc_burnin,
                    tune=mcmc_tune), silent=TRUE)

                if(class(mcmcsam) != 'try-error' | counter==100){
                    if(counter == 100) stop('Too many bad random start parameters')
                    finished = TRUE
                }else{
                    counter = counter+1
                    if(mcmc_verbose) print('-- Bad random start parameters, trying again..')
                }
            }

            return(mcmcsam)
        }

        if(.Platform$OS.type == 'windows'){
            # mclapply only works with 1 core on windows
            mc_cores = 1
        }
        mcmc_chains = parallel::mclapply(1:mcmc_nchains, mcmc_chain_maker, 
            mc.cores=mc_cores)


        ## Compute some quantiles
        ## Make sure qgammagpdcon has been updated with the bugfix
        upper_bound <-function(par){
            qf(1, par, phiu=phiu2) + data_offset
        }

        ari_max_data<-function(par){
            qf(1-1/(length(data_trans)+1), par, phiu=phiu2) + data_offset
        }

        ari_100<-function(par){
            qf(1 - 1/(100*annual_event_rate), par, phiu=phiu2) + data_offset 
        }

        # Get important statistics for all chains -- so we can check that they
        # have converged
        upper_bound_chains = parallel::mclapply(1:mcmc_nchains, 
            f<-function(x) upper_bound(mcmc_chains[[x]]),
            mc.cores=mc_cores)
        ari_max_data_chains = parallel::mclapply(1:mcmc_nchains, 
            f<-function(x) ari_max_data(mcmc_chains[[x]]), 
            mc.cores=mc_cores)
        ari_100_chains = parallel::mclapply(1:mcmc_nchains, 
            f<-function(x) ari_100(mcmc_chains[[x]]), 
            mc.cores=mc_cores)
        loglik_chains = parallel::mclapply(1:mcmc_nchains, 
            f<-function(x) apply(mcmc_chains[[x]], 1, loglik_fun), 
            mc.cores=mc_cores)

        loglik_chains_maximum = max(unlist(lapply(loglik_chains, max)))

        if(loglik_chains_maximum > -fit_optim$value){
            if(mcmc_verbose) print('Warning: Original fit was not optimal')
            chain_max = which.max(unlist(lapply(loglik_chains, max)))
            chain_ind = which.max(loglik_chains[[chain_max]])
            if(mcmc_verbose) print('Better parameters are:')
            if(mcmc_verbose) print(mcmc_chains[[chain_max]][chain_ind,])
            if(mcmc_verbose) print('Adjusting fit_optim to reflect this')
            fit_optim$par = mcmc_chains[[chain_max]][chain_ind,]
            fit_optim$value = nll_fun(fit_optim$par, x=data_trans)
        }

        # Combine chains for later analysis
        # Should not use this unless we think the chains have converged
        combined_chains = do.call(rbind, mcmc_chains)
        combined_loglik = do.call(c, loglik_chains)
        combined_ari100 = do.call(c, ari_100_chains)

    })
    #return(environment())
    return(fit_env)
}


#' Function to make a new 'qfun' in an existing mixture_fit environment,
#' using a random sample from the bayesian distribution generated
#' with mcmc_gpd_mixture
#'
#' This can be preferable to fitting the new bootstrap sample
#' (for the same reasons that Bayesian uncertainty is often 
#' preferrable to a bootstrap, if it is doable)
#'
make_random_mcmc_qfun<-function(fit_env){

    with(fit_env,{
        # Get random parameter values from the MCMC distribution
        random_par_index = sample(1:length(combined_chains[,1]), size=1)
        random_par = combined_chains[random_par_index,]

        # New qfun using the random parameters
        # See exmix_fit.R for definitions of the other parameters in the
        # environment
        qfun_random<-function(p) qf(p, random_par, phiu=phiu2) + data_offset
    })

}


#' Return level plot 
#'
#' @param fit_env environment which has mcmc_gpd_mixture applied
#' @param ci level of Bayesian quantile confidence limits which are plotted
#' @param xlim x limits for the plot. If null, these are provided. Note you
#' probably want the left limit to be larger than the right limit (since
#' we plot annual exceedance rate on the x axis)
#' @param log log parameter passed to the plot
#' @return Nothing, but fit_env is modified to include new variables required for the analysis,
#' and a nice plot is made
mcmc_rl_plot<-function(fit_env, ci=0.95, xlim=NULL, log='x'){

    if(is.null(xlim)){
        xlim = c(1 - 0.5*(length(fit_env$data)/(length(fit_env$data)+1)), 
                 0.5*(1 - length(fit_env$data)/(length(fit_env$data)+1)) )
        xlim = xlim*fit_env$annual_event_rate
    }

    fit_env$rl_plot_par = list(ci=ci, xlim=xlim, log=log)

    with(fit_env, {
        #browser()

        # Get data ready for plotting
        data_vals = sort(data)

        data_exceed_rate = 1 - (1:length(data_vals))/(length(data_vals)+1)

        mean_exceedences_per_year = data_exceed_rate * annual_event_rate

        
        ## Pr(one or more events in time T) = 1 - Pr(no events in time T)
        #data_annual_exceed_rate = 1 - exp(-data_exceed_rate * annual_event_rate)

        # Get model ready for plotting. 

        # Compute rates extending outside the range of the data
        min_rate = min(rl_plot_par$xlim)/annual_event_rate
        max_rate = max(rl_plot_par$xlim)/annual_event_rate

        desired_rates = seq(log(min_rate), log(max_rate), 
            len=50)
        desired_rates = exp(rev(desired_rates))
       
        # Save info and cleanup 
        rm(min_rate, max_rate) 

        desired_quantile_mcmc = matrix(NA, ncol = length(combined_chains), 
            nrow=length(desired_rates))
        for(i in 1:length(desired_rates)){
            desired_quantile_mcmc[i,] = qf(1-desired_rates[i], combined_chains, 
                phiu=phiu2) + data_offset
        }

        desired_lower_q = apply(desired_quantile_mcmc, 1, 
            f<-function(x) quantile(x, (1-rl_plot_par$ci)/2))
        desired_upper_q = apply(desired_quantile_mcmc, 1, 
            f<-function(x) quantile(x, 1 -0.5*(1-rl_plot_par$ci)))
        desired_median_q = apply(desired_quantile_mcmc, 1, 
            f<-function(x) quantile(x, 0.5))


        desired_best_q = qfun(1-desired_rates)

        # Plot data
        plot(mean_exceedences_per_year, data_vals, 
            xlim=rl_plot_par$xlim, 
            ylim=c(min(desired_lower_q), max(desired_upper_q)),
            t='p', cex=0.5, pch=19, log=rl_plot_par$log, 
            xlab = 'Annual exceedance rate', ylab="")

        # Fitted curve
        points(desired_rates*annual_event_rate, desired_best_q, t='l', 
            col='blue', lwd=2)
        points(desired_rates*annual_event_rate, desired_lower_q, t='l', 
            col='red', lty='dashed')
        points(desired_rates*annual_event_rate, desired_upper_q, t='l', 
            col='red', lty='dashed')
        points(desired_rates*annual_event_rate, desired_median_q, t='l',
            col='orange', lty='dashed', lwd=2)

        legend('topleft', 
            c('Data', 'ML fit', 'Bayesian median', 
                paste0('Bayesian credible interval (', rl_plot_par$ci, ')', sep="")),
            col=c('black', 'blue', 'orange', 'red'), 
            lty=c(NA, 'solid', 'dashed', 'dashed', 'dashed'),
            lwd = c(NA, 2, 2, 1, 1), 
            pch=c(19, NA, NA, NA, NA))

    })
    
}

