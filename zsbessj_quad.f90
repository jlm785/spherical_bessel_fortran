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

subroutine zsbessj_quad(n, z, zsb, acc)


  implicit none

  integer, parameter          :: REAL128 = selected_real_kind(28)

! input

  integer, intent(in)                  ::  n                             !<  n >= 0 order of function
  complex(REAL128), intent(in)         ::  z                             !<  argument

! output

  complex(REAL128), intent(out)        ::  zsb                           !<  result
  real(REAL128), intent(out)           ::  acc                           !<  estimate of accuracy

! local variables

  real(REAL128)                  ::  xpow
  real(REAL128)                  ::  ax

  real(REAL128)                  ::  acc_fwd
  complex(REAL128)               ::  zsb_fwd

  real(REAL128)                  ::  acc_pow
  complex(REAL128)               ::  zsb_pow



  xpow = min(20.0,n+0.01)
  ax = abs(real(z,REAL128))

  if(ax < xpow) then

    call zsbessj_pow_quad(n, z, zsb, acc)

  elseif(n < 30) then

    call zsbessj_fwd_quad(n, z, zsb, acc)

  else

    call zsbessj_bwd_quad(n, z, zsb, acc)

    if(acc < 20)then

      if(ax > real(n)) then

        call zsbessj_fwd_quad(n, z, zsb_fwd, acc_fwd)

        if(acc_fwd > acc) zsb = zsb_fwd

      else

        call zsbessj_pow_quad(n, z, zsb_pow, acc_pow)
        call zsbessj_fwd_quad(n, z, zsb_fwd, acc_fwd)


        if(acc_fwd > acc .and. acc_fwd >= acc_pow) zsb = zsb_fwd
        if(acc_pow > acc .and. acc_pow >  acc_fwd) zsb = zsb_pow

      endif
    endif

  endif

end subroutine zsbessj_quad




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

  integer, intent(in)                   ::  n                           !<  n >= 0 order of function
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




!>  Calculates the spherical bessel function of first kind j_n(z)
!>  for a complex argument with a forward series from n = 0,1.
!>  Formula 10.1.19 of Abramowitz and Stegun.
!>
!>  Quad precision version
!>
!>  It provides an estimate of the number of accurate digits.
!>  This estimate may be optimist for large values of z.
!>
!>  \author       Jose Luis Martins
!>  \version      0.02
!>  \date         19 April 2018, 10 February 2025.
!>  \copyright    GNU LGPL v3

subroutine zsbessj_fwd_quad(n, z, zsb, acc)

!  Written 19 April 2018 from old atom code, Siesta code, and NumSBF code.

!  Relative errors are relevant for Abs(Re(z)) ~ 0.75 n, n > 30 and
!  Abs(Re(z)) >> Abs(Im(z))

  implicit none

  integer, parameter          :: REAL128 = selected_real_kind(28)

! input

  integer, intent(in)                  ::  n                             !<  n >= 0 order of function
  complex(REAL128), intent(in)         ::  z                             !<  argument

! output

  complex(REAL128), intent(out)        ::  zsb                           !<  result
  real(REAL128) , intent(out)          ::  acc                           !<  (over)estimation of number of accurate digits)

! local variables

  complex(REAL128)              ::  by, bym, byp, uz                     !  recurrence variables
  real(REAL128)                 ::  bymax

! counter

  integer                    ::  j

! parameters

  real(REAL128), parameter          ::  ZERO = 0.0_REAL128
  real(REAL128), parameter          ::  UM = 1.0_REAL128
  real(REAL128), parameter          ::  EPS = epsilon(UM)

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

end subroutine zsbessj_fwd_quad



!>  Calculates the spherical bessel function of first kind j_n(z)
!>  for a complex argument with a backward series1.
!>  Formula 10.1.19 of Abramowitz and Stegun.
!>
!>  Quad precision version
!>
!>  Accuracy is estimated from the correction yo the value of j_0.
!>  The accuracy of j_n may be quite higher.
!>
!>  \author       Jose Luis Martins
!>  \version      0.02
!>  \date         11 February 2025.
!>  \copyright    GNU LGPL v3

subroutine zsbessj_bwd_quad(n, z, zsb, acc)


  implicit none

  integer, parameter          :: REAL128 = selected_real_kind(28)

! input

  integer, intent(in)                  ::  n                             !<  n >= 0 order of function
  complex(REAL128), intent(in)         ::  z                             !<  argument

! output

  complex(REAL128), intent(out)        ::  zsb                           !<  result
  real(REAL128) , intent(out)          ::  acc                           !<  number of accurate digits of j_0.

! local variables

  integer                       ::  nup                                  !  starting point nup

  complex(REAL128)              ::  by, bym, byp, uz                     !  recurrence variables
  complex(REAL128)              ::  zsb0

! counter

  integer                    ::  j

! parameters

  real(REAL128), parameter          ::  UM = 1.0_REAL128
  real(REAL128), parameter          ::  EPS = epsilon(UM)

  uz = UM / z

! finds starting n

  nup = nint(abs(z) / 0.78)
  if (nup < n+50) nup = n+50

! recursion formula

  call zsbessj_pow_quad(nup+1, z, byp, acc)
  call zsbessj_pow_quad(nup, z, by, acc)

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

end subroutine zsbessj_bwd_quad


