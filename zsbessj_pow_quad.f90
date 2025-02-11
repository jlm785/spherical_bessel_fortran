!------------------------------------------------------------------------------!
!  This file is free software: you can redistribute it and/or modify           !
!  it under the terms of the GNU Lesser General Public License as published by !
!  the Free Software Foundation, either version 3 of the License, or           !
!  (at your option) any later version.                                         !
!                                                                              !
!  This file is distributed in the hope that it will be useful,                !
!  but WITHOUT ANY WARRANTY; without even the implied warranty of              !
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               !
!  GNU Lesser General Public License for more details.                         !
!                                                                              !
!  You should have received a copy of the GNU Lesser General Public License    !
!  along with this program.  If not, see <https://www.gnu.org/licenses/>.      !
!------------------------------------------------------------------------------!

!>  Calculates the spherical bessel function of first kind j_n(z)
!>  for a complex argument with a power series around z = 0.
!>  Formula 10.1.2 of Abramowitz and Stegun.
!>
!>  It provides an (under?) estimate of the number of accurate digits.
!>
!>  \author       Jose Luis Martins
!>  \version      0.01
!>  \date         19 April 2018, 10 February 2025.
!>  \copyright    GNU LGPL v3

subroutine zsbessj_pow_quad(n, z, zsb, acc)

!  Written 19 April 2018 from old atom code, Siesta code, and NumSBF code.

!  Relative errors are relevant for Abs(Re(z)) ~ 0.75 n, n > 30 and
!  Abs(Re(z)) >> Abs(Im(z))

  implicit none

  integer, parameter          :: REAL128 = selected_real_kind(28)

! input

  integer, intent(in)                  ::  n                           !<  n >= 0 order of function
  complex(REAL128), intent(in)          ::  z                           !<  argument

! output

  complex(REAL128), intent(out)         ::  zsb                         !<  result
  real(REAL128) , intent(out)           ::  acc                         !<  (under)estimation of number of accurate digits)

! local variables

  complex(REAL128)               ::  pref                               !  prefactor of series expansion
  complex(REAL128)               ::  z2                                 !  0.5*z^2
  complex(REAL128)               ::  sumz, fac                          !  series expansion
  logical                        ::  fail
  real(REAL128)                  ::  facmax
  integer                        ::  jmax

! counter

  integer                    ::  j

! parameters

  real(REAL128), parameter           ::  UM = 1.0_REAL128
  real(REAL128), parameter           ::  EPS = epsilon(UM)

! series expansion

  pref = UM
  if(n > 0) then
    do j = 1,n
      pref = pref*z/(2*j+1)
    enddo
  endif

! maximum j can be estimated from Stirling formula

  z2 = z*z/2
  sumz = UM
  fac = UM

  fail = .TRUE.
  facmax = fac

  jmax = 40 + nint(1.5*abs(z))

  do j = 1,jmax
    fac = -fac*z2 / (j*(2*n+2*j+1))
    sumz = sumz + fac
    if(abs(fac) > facmax)  facmax = abs(fac)

    if(abs(fac) < EPS .and. abs(z2) < real(4*j*j)) then
       fail = .FALSE.

       exit

     endif
  enddo

  zsb = sumz*pref
  acc = facmax / (abs(sumz)+EPS)
  acc = -log10(acc*EPS)

  if(fail) acc = 0

  return

end subroutine zsbessj_pow_quad
