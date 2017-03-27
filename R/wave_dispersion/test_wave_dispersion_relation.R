source('wave_dispersion_relation.R', local=TRUE)

#' Test code
test_airy_wavelength_period<-function(){

    ## Shallow water limit
    wl = 20000
    h = 20
    period0 = airy_period(wl, h)
    wavelength0 = airy_wavelength(period0, h)

    stopifnot(abs(wl - wavelength0) < 1.0e-05*wavelength0)

    stopifnot(abs(period0 - wl/sqrt(9.81 * h)) < 1.0e-05 * period0)
    print('PASS')

    ## Deep water limit
    wl = 2
    h = 200
    
    period0 = airy_period(wl, h)
    wavelength0 = airy_wavelength(period0, h)

    stopifnot(abs(wl - wavelength0) < 1.0e-05)

    stopifnot(abs(wl - period0**2 * 9.81/(2*pi)) < 1.0e-05 * wl)
    print('PASS')

    ## Intermediate
    
    wl = 200
    h = 20
    
    period0 = airy_period(wl, h)
    wavelength0 = airy_wavelength(period0, h)

    stopifnot(abs(wl - wavelength0) < 1.0e-05)
    print('PASS')

    # Check that the dispersion relation is satisfied to within an expected tolerance
    resid = wavelength0 - period0**2 * 9.81 /(2*pi) * tanh( h * 2 * pi / wavelength0)
    stopifnot(abs(resid) < 1.0e-06)
    print('PASS')

}

test_airy_wavelength_period()
