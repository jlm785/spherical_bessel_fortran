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
!>  for a complex argument with a forward series from n = 0,1.
!>  Formula 10.1.19 of Abramowitz and Stegun.
!>
!>  It provides an estimate of the number of accurate digits.
!>  This estimate may be optimist for large values of z.
!>
!>  \author       Jose Luis Martins
!>  \version      0.02
!>  \date         19 April 2018, 10 February 2025.
!>  \copyright    GNU LGPL v3

subroutine zsbessj_fwd(n, z, zsb, acc)

!  Written 19 April 2018 from old atom code, Siesta code, and NumSBF code.

!  Relative errors are relevant for Abs(Re(z)) ~ 0.75 n, n > 30 and
!  Abs(Re(z)) >> Abs(Im(z))

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                  ::  n                             !<  n >= 0 order of function
  complex(REAL64), intent(in)          ::  z                             !<  argument

! output

  complex(REAL64), intent(out)         ::  zsb                           !<  result
  real(REAL64) , intent(out)           ::  acc                           !<  (over)estimation of number of accurate digits)

! local variables

  complex(REAL64)               ::  by, bym, byp, uz                     !  recurrence variables
  real(REAL64)                  ::  bymax

! counter

  integer                    ::  j

! parameters

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64
  real(REAL64), parameter           ::  UM = 1.0_REAL64
  real(REAL64), parameter           ::  EPS = epsilon(UM)

  uz = UM / z

! recursion formula

  if(n == 0) then
    zsb = sin(z) * uz
    acc = UM
  else
    by = sin(z) * uz
    bym = cos(z) * uz
    byp = uz*by - bym
    bymax = abs(uz*by)
    if(abs(bym) > bymax) bymax = abs(bym)
    bym = by
    by = byp
    if( n > 1) then
      do j = 1,n-1
        bymax = (2*j+1)*bymax*abs(uz)
        byp = (2*j+1)*uz*by - bym
        if(abs((2*j+1)*uz*by) > bymax) bymax = abs((2*j+1)*uz*by)
        if(abs(bym) > bymax) bymax = abs(bym)
        bym = by
        by = byp
        if(abs(by) > bymax)  bymax = abs(by)
      enddo
    endif

    zsb = by

    acc = bymax / (abs(by)+EPS)

  endif

  acc = -log10(acc*EPS)

  return

end subroutine zsbessj_fwd
