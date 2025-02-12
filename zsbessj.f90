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
!>  for a complex argument.
!>  Formulas 10.1.2 and 10.1.19 of Abramowitz and Stegun.
!>
!>  For reasonable values of n and z it has almost machine precision.
!>
!>  \author       Jose Luis Martins
!>  \version      0.02
!>  \date         12 February 2025.
!>  \copyright    GNU LGPL v3

function zsbessj(n, z)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                  ::  n                             !<  n >= 0 order of function
  complex(REAL64), intent(in)          ::  z                             !<  argument

! output

  complex(REAL64)                      ::  zsbessj                       !<  result

! local variables

  real(REAL64)                  ::  acc
  integer                       ::  xpow
  real(REAL64)                  ::  ax


  xpow = min(10.0,n+0.01)
  ax = abs(real(z,REAL64))

  if(ax < xpow) then

    call zsbessj_pow(n, z, zsbessj, acc)

  elseif(n < 10) then

    call zsbessj_fwd(n, z, zsbessj, acc)

  else

    call zsbessj_bwd(n, z, zsbessj)

  endif



end function zsbessj
