# spherical_bessel_fortran

Spherical Bessel funstions in fortran90

## function sbessj(n,x)

Calculates the spherical bessel function of first kind j_n(x)
Formulas 10.1.2, 10.1.19 and 10.1.59 of Abramowitz and Stegun
SphericalBesselJ[n,x] of Mathematica

I have been using this function for decades.  It is a compromise in code-size, precision and reliability.
In the course of upgrading the code from f77 to f90, I cleaned it, tested it against the results of Mathematica and other software. I also added a test for the parameter range where the results are reliable, in the sense that it provides at least 8 correct significant figures or terminates with an error message.  For all commom cases it as 15 digits accuracy.

To test the subroutine just compile
F90 test_bessel.f90 sbessj.f90 -o test.exe
and run test.exe

It will compare the present subroutine with a subroutine from Shanjie Zhang, Jianming Jin
