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
!>  for a complex argument with a backward series1.
!>  Formula 10.1.19 of Abramowitz and Stegun.
!>
!>  Accuracy is estimated from the value of j_0.
!>
!>  \author       Jose Luis Martins
!>  \version      0.02
!>  \date         11 February 2025.
!>  \copyright    GNU LGPL v3

subroutine zsbessj_bwd(n, z, zsb, acc)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                  ::  n                             !<  n >= 0 order of function
  complex(REAL64), intent(in)          ::  z                             !<  argument

! output

  complex(REAL64), intent(out)         ::  zsb                           !<  result
  real(REAL64) , intent(out)           ::  acc                           !<  estimation of number of accurate digits)

! local variables

  integer                       ::  nup                                  !  starting point nup

  complex(REAL64)               ::  by, bym, byp, uz                     !  recurrence variables
  complex(REAL64)               ::  zsb0

! counter

  integer                    ::  j

! parameters

  real(REAL64), parameter           ::  UM = 1.0_REAL64

  uz = UM / z

! finds starting n

  nup = nint(abs(z) / 0.78)
  if (nup < n+5) nup = n+5

! recursion formula

  call zsbessj_pow(nup+1, z, byp, acc)
  call zsbessj_pow(nup, z, by, acc)

  do j = nup,1,-1
    bym = (2*j+1)*uz*by - byp
    byp = by
    by = bym
    if(j-1 == n) zsb = by
  enddo

  zsb0 = sin(z) * uz
  zsb = zsb * zsb0/by

  acc = -log10(abs(abs(zsb0/by)-UM))


  return

end subroutine zsbessj_bwd
