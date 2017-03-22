#
# Code for cross-validation of copula selection
#


#' Use xvCopula_2D for copula selection, with an interface similar to BiCopSelect
#'
#' @param u1 vector of data which will be transformed to psuedo-observations
#' @param u2 vector of data which will be transformed to psuedo-observations
#' @param familyset vector of integers corresponding to families to test, see ?BiCopName
#' @param k perform k-fold cross-validation. If NULL, do leave-one-out cv (slow)
#' @param xv_only Logical. If TRUE, only return the families with their xvCopula_2D score. Otherwise,
#' return the best fitted copula
#' @return if xv_only, return the best fit copula. Otherwise, return a 2 column
#' matrix giving the families, and their score.
#'
BiCopSelect_xv<-function(u1, u2, familyset=NULL, k=NULL, verbose = FALSE, xv_only=FALSE){

    all_families = c(0:10, 13, 14, 16:20, 23:24, 26:30, 33:34, 36:40, 104, 114, 124, 134, 
        204, 214, 224, 234)
   
    if(is.null(familyset)){
        familyset = all_families
    }

    if( is.character(familyset)){
        if(familyset == 'one_parameter_models'){
            familyset = c(0, 1, 3:6, 13:14, 16, 23:24, 26, 33:34, 36) 
        }else{
            stop('unrecognized character familyset')
        }
    }

    mm = match(familyset, all_families)
    if( any( is.na(mm) )){
        print(c('Cannot match copula family: ', familyset[which(is.na(mm))]))
    }

    family_xv_value = familyset * NA 

    for(i in 1:length(familyset)){
        family_xv_value[i] = xvCopula_2D(familyset[i], cbind(u1, u2), k=k, verbose=verbose)
    }

    if(xv_only){
        output = cbind(familyset, family_xv_value)
        colnames(output) = c('family', 'xvCopulaValue')
    }else{
        output = familyset[which.max(family_xv_value)]
        fitted_copula = BiCopEst(u1, u2, family=output)
        output = fitted_copula
    }

    return(output)

}


#' Edit of copula::xvCopula to work with bivariate copulas from VineCopula
#' 
#' The main changes are A) Making sure the data dimension = 2, and B) using
#' BiCopEst and BiCopPDF inplace of 'copula' package equivalent functions.
#' Note from the test function .test_BiCopSelect_xv, that I don't get exactly
#' the same answer as using 'copula' routines. This seems to be because the
#' fitting and/or PDF computations are not identical in 'copula' and 'VineCopula'.
#' However, I have not seen significant differences.
#' 
#' @param copula_family integer corresponding to a bivariate copula family, see '?BiCopName'
#' @param x matrix with 2 columns in (0,1), giving 2D data which will be converted to psuedo observations.
#' @param k if not NULL, use 'leave-one-out' cross validation, otherwise do k-fold cross-validation. In the
#' latter case, the result will depend on the value of the random seed.
#' @param verbose logical. Print progress info
xvCopula_2D<-function (copula_family, x, k = NULL, verbose = interactive(), ...){

    if (!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        x <- as.matrix(x)
        stopifnot(is.matrix(x))
    }

    d <- 2
    n <- nrow(x)
    stopifnot(is.numeric(x), (d == 2), (n > 0), ncol(x) == 2) #, dim(copula) == d)

    if(is.null(k)){
        k = n
    }else{
        k = as.integer(k)
    }
    
    stopifnot(k >= 2L, n%/%k >= 1)

    if (k < n) x <- x[sample.int(n), ]

    if (verbose) {
        pb <- txtProgressBar(max = k, style = if (isatty(stdout())) 
            3
        else 1)
        on.exit(close(pb))
    }

    xv <- 0

    p <- n%/%k

    m <- rep(p, k)

    r <- n - k * p

    if (r > 0) 
        m[seq_len(r)] <- p + 1

    b <- c(0, cumsum(m))

    v <- matrix(NA, p + 1, d)

    for (i in seq_len(k)) {

        sel <- (b[i] + 1):b[i + 1]

        x.not.s <- x[-sel, , drop = FALSE]
        u <- pobs(x.not.s)
        #u <- x.not.s

        #copula <- fitCopula(copula, u, estimate.variance = FALSE, 
        #    ...)@copula

        ## Testing
        ##copula <- fitCopula(normalCopula(), u, estimate.variance = FALSE, 
        ##    ...)@copula

        ## Try to suppress printing
        tmpchar<-capture.output({
            copula <- try(BiCopEst(u[,1], u[,2], family = copula_family), silent=TRUE)
        })
        # Some copulas will fail, e.g. if trying to fit a copula which only has
        # negative rank correlation, to data with positive rank correlation. 
        # Setting the output to -Inf should ensure we don't select this one!
        if( (class(copula) == 'try-error') || 
            ( (length(tmpchar) > 0) & 
              (length(grep('cannot be used', tmpchar)) > 0) )){
            return(-Inf)
        }

        imi <- seq_len(m[i])

        v.i <- v[imi, , drop = FALSE]

        x.sel <- x[sel, , drop = FALSE]

        nmi <- n - m[i]

        for (j in seq_len(d)) {
            vals <- unique.default(xj <- sort(x.not.s[, j]))
            v.i[, j] <- approxfun(vals, cumsum(tabulate(match(xj, 
                vals))), method = "constant", ties = "ordered", 
                yleft = 0, yright = nmi)(x.sel[, j])
        }

        v.i <- (v.i + 1/2)/(nmi + 1)

        #xv <- xv + mean(dCopula(v.i, copula, log = TRUE))
        xv <- xv + mean(log(BiCopPDF(v.i[,1], v.i[,2], copula)))

        if (verbose) 
            setTxtProgressBar(pb, i)
    }
    output = xv/k * n
    return(output)
}

#' Check that it works by comparison with similar functions in 'copula'
#'
#' We do not get exact agreement for most datasets. This seems to be because
#' of differences in parameter estimation or log-likelihood computation
#' between VineCopula and 'copula'. [Based on tests where I hacked the code to
#' use copula's fitting routines]
#' 
.test_BiCopSelect_xv<-function(){

    # Check a gaussian copula, by comparison with 'copula' package
    m1 = rnorm(100)
    m2 = rnorm(100) + 0.2*m1

    fit1 = BiCopSelect_xv(m1, m2, familyset = 1, k=NULL, verbose=FALSE, xv_only=TRUE)

    fit2 = copula::xvCopula(copula::normalCopula(), cbind(m1, m2), k=NULL)


    if(abs(fit1[,2] - fit2) < abs(fit2) * 0.001){
        print('PASS')
    }else{
        print('FAIL')
    }    

   
    # Check a frank copula with the same data 
    fit1 = BiCopSelect_xv(m1, m2, familyset = 5, k=NULL, verbose=FALSE, xv_only=TRUE)

    fit2 = copula::xvCopula(copula::frankCopula(), cbind(m1, m2), k=NULL)

    if(abs(fit1[,2] - fit2) < abs(fit2) * 0.001){
        print('PASS')
    }else{
        print('FAIL')
    }    
}
