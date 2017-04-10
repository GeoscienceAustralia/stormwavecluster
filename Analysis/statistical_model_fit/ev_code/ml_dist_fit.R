# Code to assist with fitting scale-shape-location models with covariates

# Get gpd and pe3
source('pe3.R', local=TRUE)
source('gpd.R', local=TRUE)
source('wbl.R', local=TRUE)


#
# Negative log likelihood functions with a particular interface
#
gpd_nll<-function(scale, shape, location, data){

    if(any(is.nan(scale) | is.nan(shape) | is.nan(location))) return(Inf)
    if(any(scale <= 0)) return(Inf)

    -sum(dgpd(data, location=location, scale=scale, shape=shape, log=TRUE))
}

pe3_nll<-function(scale, shape, location, data){
    if(any(is.nan(scale) | is.nan(shape) | is.nan(location))) return(Inf)

    -sum(dpe3(data, location=location, scale=scale, shape=shape, log=TRUE))
}

wbl_nll<-function(scale, shape, location, data){
    if(any(is.nan(scale) | is.nan(shape) | is.nan(location))) return(Inf)

    -sum(dwbl(data, scale=scale, shape=shape, location=location, log=TRUE))
}


## General function to pass parameters to a scale/shape/location model,
## with linear predictors
#' @param par Vector of parameters. = c(scale_parameters, shape_parameters)
#' @param location fixed value for location parameter
#' @param data data to which we fit model
#' @param preds a matrix with columns containing predictors
#' @param scale_pred vector of column indices to use as linear predictors for scale
#' @param shape_pred vector of column indices to use a linear predictors for shape
#' @param scale_link. Function. Scale is modelled as scale_link(predictors%*%scale_pred)
#' @param shape_link. Function. Shape is modelled as shape_link(predictors%*%shape_pred)
#' @param nll_fun a function returning the negative log likelihood with the same
#' interface as gpd_nll
#' @param return_nll if TRUE return the negative log likelihood. Otherwise return
#' the scale and shape parameters
#' @param reparameterise_scale. Numeric constant. Transform scale to 'scale + reparamterise_scale * shape'
#' Theory suggests then the parameterization should be more invarient to changes
#' in the threshold, if reparameterise_scale = new_threshold - old_threshold
#' @return
model_nll<-function(par, location, data=NULL, preds=NULL, scale_pred=NULL, shape_pred=NULL, 
    scale_link = identity, shape_link = identity, nll_fun=gpd_nll, return_nll=TRUE, 
    reparameterise_scale=0){

    stopifnot(length(par) == (2+length(scale_pred) + length(shape_pred)))

    len_scale_pred = length(scale_pred)
    len_shape_pred = length(shape_pred)

    # Scale predictor
    if(is.null(scale_pred)){
        scale = scale_link(par[1])
    }else{
        scale = scale_link(par[1] + par[2:(2+len_scale_pred-1)]%*%t(preds[,scale_pred]))
    }

    # Shape predictor
    if(is.null(shape_pred)){
        shape = shape_link(par[2+len_scale_pred])
    }else{
        shape = shape_link(par[2+len_scale_pred] + par[(2+len_scale_pred)+(1:len_shape_pred)]%*%t(preds[,shape_pred]))
    }

    scale = scale + reparameterise_scale*shape 

    if(return_nll){
        output = nll_fun(scale, shape, location, data)
        return(output)
    }else{
        return(data.frame(scale=c(scale), shape=c(shape), location=c(shape)*0 + location))
    }
}

#' Convenience wrapper to optimize model_nll
#'
#' Returns optim object with a few other useful bits of info
#' See model_nll for documentation
#'
fit_ssl_distr<-function(par, location, data, preds, scale_pred=NULL, shape_pred=NULL, 
    scale_link = identity, shape_link = identity, distr='gpd', skip_hessian=FALSE, 
    reparameterise_scale=0., ...){
    
    if(any(data < location)){
        print('Warning: Removing data < location')
        keep = which(data > location)
        if(length(keep) == 0) stop('location > max(data)')
        data = data[keep]
        if(length(dim(preds))>1){
            preds = preds[keep,, drop=FALSE]
        }else{
            preds = preds[keep]
        }
    }else{
        keep = 1:length(data)
    }
  
    # Get negative log likelihood 
    nll_fun_name = paste0(distr, '_nll') 
    nll_fun = get(nll_fun_name)
   
    negloglik_fun <-function(par){
        model_nll(par, location=location, data=data, preds=preds, 
           scale_pred = scale_pred, shape_pred = shape_pred, scale_link=scale_link, 
           shape_link = shape_link, nll_fun=nll_fun, 
           reparameterise_scale = reparameterise_scale)  
        }

    loglik_fun<-function(par) -negloglik_fun(par)

    #model_fit = optim(par, fn=model_nll, location=location, data=data, 
    #    preds=preds, scale_pred=scale_pred, shape_pred=shape_pred, scale_link=scale_link, 
    #    shape_link = shape_link, nll_fun=nll_fun, hessian=!skip_hessian, 
    #    reparameterise_scale=reparameterise_scale, ...)
    model_fit = optim(par, fn=negloglik_fun, hessian=!skip_hessian, ...)

    model_fit$data = data
    model_fit$keep = keep

    model_fit$location = location
    model_fit$nll_fun = nll_fun
    model_fit$distr = distr
    
    model_fit$negloglik_fun = negloglik_fun
    model_fit$loglik_fun = loglik_fun

    # Make a function which will turn a matrix of predictors into scale/shape parameters
    model_fit$data_pars_function<-function(preds, par=model_fit$par){
        model_nll(par, location, data, preds, 
            scale_pred = scale_pred, shape_pred = shape_pred, scale_link=scale_link, 
            shape_link = shape_link, nll_fun=nll_fun, return_nll=FALSE,
            reparameterise_scale = reparameterise_scale)  
    }
    
    model_fit$data_pars = model_fit$data_pars_function(preds)

    #
    model_fit$aic = 2*length(model_fit$par) + 2*model_fit$value

    if(!skip_hessian){
        model_fit$se = sqrt(diag(solve(model_fit$hessian)))
    }else{
        model_fit$se = model_fit$par * NA
    }

    if(model_fit$convergence != 0) warning('Not converged')
    

    return(model_fit)
}

#'
#' Automated test for ml_dist_fit. Compare against the results if 'ismev'
#' which fits a smaller range of models.
#'
test_ml_dist_fit<-function(){

    library(ismev)

    mydat = rgpd(1000, scale=10, location=10, shape=0.2)

    rand_dat = runif(1000)

    for(threshold in c(10, 11)){

        ismev_fit = gpd.fit(mydat, threshold=threshold, ydat=data.frame(rand_dat), sigl=1,shl=1, show=FALSE)

        newfit = fit_ssl_distr(c(8, 0, 0, 0), location=threshold, data=mydat, preds=matrix(rand_dat, ncol=1), 
            shape_pred=1, scale_pred=1, method='Nelder-Mead')

        tester = isTRUE(all.equal(newfit$par, ismev_fit$mle, tol=1.0e-02))

        if(tester){
            print('PASS')
            print(newfit$value)
            print(ismev_fit$nllh)
        }else{
            if(newfit$value < ismev_fit$nllh){
                print('Better than gpd.fit')
            }else{
                print('Worse than gpd.fit')
            }
            print(ismev_fit$mle)
            print(newfit$par)
            print(ismev_fit$nllh)
            print(newfit$value)
        }
    }

}

#' QQplot
#'
#' Alternative to R's qqplot, which avoids an artefact for datasets with very
#' unequal size. The key change is commented in the code.
#' 
qqplot2<-function(x, y, plot.it = TRUE, xlab = deparse(substitute(x)), ylab = deparse(substitute(y)), ...){
    sx <- sort(x)
    sy <- sort(y)
    lenx <- length(sx)
    leny <- length(sy)
    if (leny < lenx) 
        sx <- quantile(sx, p = (1:leny)/(leny+1), type=6) #approx(1L:lenx, sx, n = leny)$y
    if (leny > lenx) 
        sy <- quantile(sy, p = (1:lenx)/(lenx+1), type=6) #approx(1L:leny, sy, n = lenx)$y
    if (plot.it) 
        plot(sx, sy, xlab = xlab, ylab = ylab, ...)
    invisible(list(x = sx, y = sy))
}


#' Useful diagnostic plots for a model
#'
#' The distribution must have a random number generator rdist (e.g. for gpd, it is rgpd)
#' with the standard scale, shape, location interface [e.g. rgpd(n, scale=scale,
#' shape=shape, location=location)]
#'
#' @param best_hsig_model A model fit using 'fit_ssl_model'
#' @param input_data_time A vector of times corresponding to the data that fit_ssl_model was fit to.
#' Times for all data provided to that function should be included, even if they correspond to data 
#' values < location which would have been removed
#' @param nrand Use random sample this many times larger than the data to
#' compute key model predictions without complex integration
#' @return Nothing, but make a nice diagnostic plot
#'
plot_ssl_fit<-function(best_hsig_model, input_data_time,
    nrand=1000){

    # Get random number generator
    rdistr = get(paste0('r', best_hsig_model$distr))
    
    # Diagnostic plots to compare the model and data
    par(mfrow=c(3,4))

    # QQ-plot of the 'best' fitting model
    # Avoid the interpolation artefact of the standard qqplot
    # Make a long random dataset to compare against
    random_data = rdistr(length(best_hsig_model$data)*nrand, 
        location=best_hsig_model$data_pars[,3], 
        scale=best_hsig_model$data_pars[,1], 
        shape=best_hsig_model$data_pars[,2])
    qqplot2(best_hsig_model$data, random_data)
    abline(0, 1, col='red')
    grid(col='brown')
    title(main='QQ-plot (large model sample) ')

    # Look at what our fitted model suggests the maximum data value could be.
    max_vals = apply(matrix(random_data, ncol=nrand, byrow=FALSE), 2, max)
    print(paste0(
        'Distribution of max data value according to the model (nrand = ', 
        nrand, '):'))
    print(summary(max_vals))
    hist(max_vals, main='Maximum observed value \n as predicted by model', n=50)
    abline(v=max(best_hsig_model$data), col='red')


    # Shape/scale parameters over time
    years = strptime(
        paste(format(input_data_time, '%Y'), '-01-01 00:00:00', sep=""), 
        format='%Y-%m-%d %H:%M:%S')
    years = unique(years)

    t0 = input_data_time[best_hsig_model$keep]

    plot(t0, rep(best_hsig_model$data_pars[,1], length.out=length(t0)), t='l', 
        main='Fitted scale parameter')
    for(i in 1:length(years)){
        points(c(years[i], years[i]), c(-100, 100), t='l', col='red')
    }
    plot(t0, rep(best_hsig_model$data_pars[,2], length.out=length(t0)), t='l', 
        main='Fitted shape parameter')
    for(i in 1:length(years)){
        points(c(years[i], years[i]), c(-100, 100), t='l', col='red')
    }

    for(i in 1:4){
        random_data = rdistr(
                length(best_hsig_model$data)*1, 
                location=best_hsig_model$data_pars[,3], 
                scale=best_hsig_model$data_pars[,1], 
                shape=best_hsig_model$data_pars[,2])
        
        qqplot(best_hsig_model$data, random_data,
            xlab='Data', ylab='Random sample',
            main=paste0('QQplot: Data vs random sample ', i))
        abline(0, 1, col='red')
        grid(col='brown')

        # Monthly boxplots. For this we need the times of the original data
        data_times = input_data_time[best_hsig_model$keep]
        random_data_times = rep(input_data_time, 
            length.out = length(random_data))
        boxplot(
            best_hsig_model$data ~ as.numeric(format(data_times, '%m')), 
            col=rgb(1, 0, 0, alpha=0.5), border='red', 
            ylim=c(min(best_hsig_model$data_pars[,3]), 
                   max(c(best_hsig_model$data, random_data))))
        boxplot(random_data ~ as.numeric(format(random_data_times, '%m')),
            add=TRUE, col=rgb(0, 0, 1, alpha=0.5), border='blue')
        title(main='Hsig ~ month (data: red, model sample: blue)')
    }
}
