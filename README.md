# spherical_bessel_fortran (real and complex)

Spherical Bessel funstions in fortran90

## function sbessj(n,x)

Calculates the spherical Bessel function of first kind j_n(x) for real x
Formulas 10.1.2, 10.1.19 and 10.1.59 of Abramowitz and Stegun
SphericalBesselJ[n,x] of Mathematica

I have been using this function for decades.  It is a compromise in code-size, precision and reliability.
In the course of upgrading the code from f77 to f90, I cleaned it, tested it against the results of Mathematica and other software. I also added a test for the parameter range where the results are reliable, in the sense that it provides at least 8 correct significant figures or terminates with an error message.  For all commom cases it as 15 digits accuracy.

The exception are the zeroes of the spherical Bessel function where you have absolute precision but loose
the relative precision

## function zsbessj(n,z)

Calculates the spherical Bessel function of first kind j_n(z) for complex z

It is an adaptation of the previous real version discarding the compromises for efficiency.
It is new, so caveat emptor, but it was tested for -100 < Re(z) < 100, -100 < Im(z) < 100, n=0,1,2,3,...,100,
always returning values with 14 decimal places of absolute precision or relative precision.

Relative precision is relevant  for |j_n(z)| > 1 (|Im(z)|>~6),
absolute precision is relevant for |j_n(z)| < 1 (near zeroes of j_n (z) on the real axis).

## how to run

To test the sbessj subroutine just compile

F90 test_bessel.f90 sbessj.f90 -o test.exe

where F90 is your compiler and run test.exe.
It will compare the present subroutine with a subroutine from Shanjie Zhang, Jianming Jin

To test the zsbessj subroutine just compile

F90 test_zsbessel_one.f90 zsbessj.f90 zsbessj_pow_quad.f90 -o test.exe

and run test.exe.  It will use a quad precision subroutine to compare with the double precision subroutine.
