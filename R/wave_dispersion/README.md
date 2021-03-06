**wave_dispersion_relation**
----------------------------

Code to solve the wave dispersion relation for airy wave theory (i.e. linear
wave theory).

The solution satisfies the implicit equation:

    lambda = g/(2 pi) * T^2 * tanh ( 2 * pi * h / lambda)

where `lambda` is the wavelength (m), `T` is the period (s), `h` is the water depth (m), and `g` is gravity (m/s/s).

In practice we want to compute the wave period given the wavelength and depth,
and/or compute the wavelength given the wave period and depth.

The key routines are: 
```r
airy_wavelength(period, h, g=9.81, full_output=FALSE)
airy_period(wavelength, h, g=9.81, full_output=FALSE)
```
They accept vectorized period/h/wavelength arguments. Note that `full_output` can be used
to give detailed information on the numerical minimization (probably only required for
debugging)


**USAGE**

Here we give an example of computing the wave period given the wavelength and
depth, for a wave in-between the shallow and deep water limits. We then back-calculate
the wavelength from the computed period, and check that it is identical to the input value
(to within errors tolerated by the numerical optimization). To run this code you should
be in the same directory as the file
[wave_dispersion_relation.R](wave_dispersion_relation.R), or else change the filename in the
`source` command below to include the full directory path to the latter script.
```r
source('wave_dispersion_relation.R')


wl = 200 # Wave length
h = 20 # Water depth

# Compute the wave period 
period0 = airy_period(wl, h)

# Back-calculate the wavelength (should be 200)
wavelength0 = airy_wavelength(period0, h)

# Check
stopifnot(abs(wl - wavelength0) < 1.0e-05)
```

**TESTS**

Open R in the current directory, and run:
```r
source('test_wave_dispersion_relation.R')
```

It should print PASS a few times (with no FAIL's)
