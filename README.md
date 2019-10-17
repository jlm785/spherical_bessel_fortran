# spherical_bessel_fortran

Spherical Bessel funstions in fortran90

## function sbessj(n,x)

Calculates the spherical bessel function of first kind j_n(x)
Formulas 10.1.2, 10.1.19 and 10.1.59 of Abramowitz and Stegun
SphericalBesselJ[n,x] of Mathematica

I have been using this function for decades.  It is a compromise in code-size, precision and reliability.
In the course of upgrading the code from f77 to f90, cleaned it and added a test for the parameter range where the results are reliable, in the sense that it provides at least 8 correct significant figures or terminates with an error message.  For all c
