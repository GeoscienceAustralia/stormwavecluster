**wave_dispersion_relation**
----------------------------

Code to solve the wave dispersion relation for airy wave theory (i.e. linear
wave theory).

In practice this means computing the wave period given the wavelength and depth,
and computing the wavelength given the wave period and depth.

The key routines are: 

    airy_wavelength(period, h, g=9.81, full_output=FALSE)
    airy_period(wavelength, h, g=9.81, full_output=FALSE)

They accept vectorized period/h/wavelength arguments.


**USAGE**

Here we give an example of computing the wave period given the wavelength and
depth, for a wave in-between the shallow and deep water limits. We then back-calculate
the wavelength from the computed period, and check that it is identical to the input value
(to within errors tolerated by the numerical optimization).


    source('wave_dispersion_relation.R')


    wl = 200 # Wave length
    h = 20 # Water depth
   
    # Compute the wave period 
    period0 = airy_period(wl, h)

    # Back-calculate the wavelength (should be 200)
    wavelength0 = airy_wavelength(period0, h)

    # Check
    stopifnot(abs(wl - wavelength0) < 1.0e-05)


**TESTS**

From within R, run:

    source('test_wave_dispersion_relation.R')

It should print PASS a few times (with no FAIL's)
