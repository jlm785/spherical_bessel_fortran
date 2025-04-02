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
!>
!>  Calls power, forward or backwards algorithms depending on the values of n and z
!>
!>  For reasonable values of n and z it has almost machine precision.
!>
!>  \author       Jose Luis Martins
!>  \version      0.1
!>  \date         12 February 2025.
!>  \copyright    GNU LGPL v3

function zsbessj(n, z)

!  use, intrinsic :: ieee_arithmetic

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
! input

  integer, intent(in)                  ::  n                             !<  n >= 0 order of function
  complex(REAL64), intent(in)          ::  z                             !<  argument

! output

  complex(REAL64)                      ::  zsbessj                       !<  result

! local variables

  real(REAL64)                  ::  acc
  real(REAL64)                  ::  xpow
  real(REAL64)                  ::  zr, zi, zm
  integer                       ::  nfwd
  real(REAL64)                  ::  zimax

  complex(REAL64)               ::  zsb                                  !  partial result

! constants

  real(REAL64), parameter           ::  UM = 1.0_REAL64
  real(REAL64), parameter           ::  VERYBIG = huge(UM)


  zr = abs(real(z,REAL64))
  zi = abs(aimag(z))
  zm = abs(z)

  zimax = log(VERYBIG) - 5*UM

! Exceptional cases

  if(n < 0) then
    write(6,*)
    write(6,*) '   Stopped in zsbessj, n = ',n,' < 0'
    write(6,*)

    stop

! if you do not want to stop use the following lines instead (valid without debugging options)

!     zr = 0.01 / VERYBIG
!     zi = 0.02 / VERYBIG
!     zsbessj = zr / zi
!
!     return

!   or the following lines (using ieee_arithmetic module)

!    zsbessj = ieee_value( UM, ieee_quiet_nan )
!
!    return

  endif

  if(zi > zimax) then

    write(6,*)
    write(6,*) '   Stopped in zsbessj, result is computer infinity'
    write(6,*)

    stop

! if you do not want to stop use the following lines instead (valid without debugging options)

!     zsbessj = 10.0*VERYBIG
!
!     return

!   or the following lines (using ieee_arithmetic module)

!    zsbessj = ieee_value( UM, ieee_positive_inf )
!
!    return

  endif

! heuristic range where power is safer

  if(n >= 75) then
!   instability for large numbers
    xpow = 0.6*(n-25)
  elseif(n >= 15) then
!   bridge the two extremes
    xpow = 0.31*n + 6.754
  else
!   2/3 times the asymptotic expansion of the first maxima of j_n(z) (
!   Formula 9.5 of Abramowitz and Ategun
    xpow = 0.65*(n+0.5 + 0.8086*(n+0.5)**0.3333333 + 0.07249/(n+0.5)**0.3333333)
  endif

! forward is safe for small n values
! can use forward for other ranges...

  nfwd = 1

  if(zm < xpow) then

    call zsbessj_pow(n, z, zsb, acc)

  elseif(n <= nfwd) then

    call zsbessj_fwd(n, z, zsb, acc)

  else

    call zsbessj_bwd(n, z, zsb, acc)

  endif

  zsbessj = zsb

end function zsbessj





!>  Calculates the spherical bessel function of first kind j_n(z)
!>  for a complex argument with a backward series1.
!>  Formula 10.1.19 of Abramowitz and Stegun.
!>
!>  Accuracy is estimated from the correction to the value of j_0 or j_1.
!>  The accuracy of j_n may be quite higher.
!>
!>  \author       Jose Luis Martins
!>  \version      0.1
!>  \date         11 February 2025. 27 March 2025.
!>  \copyright    GNU LGPL v3

subroutine zsbessj_bwd(n, z, zsb, acc)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                  ::  n                             !<  n >= 0 order of function
  complex(REAL64), intent(in)          ::  z                             !<  argument

! output

  complex(REAL64), intent(out)         ::  zsb                           !<  result
  real(REAL64) , intent(out)           ::  acc                           !<  number of accurate digits of j_0.

! local variables

  integer                       ::  nup                                  !  starting point nup

  complex(REAL64)               ::  by, bym, byp, uz                     !  recurrence variables
  complex(REAL64)               ::  zsb0, zsb1
  complex(REAL64)               ::  renorm0, renorm1                     !  renormalizations (avoid overflow, do not "simplify" the code!)
  real(REAL64)                  ::  zm, zr, zi

! counter

  integer                    ::  j

! parameters

  real(REAL64), parameter           ::  UM = 1.0_REAL64
  real(REAL64), parameter           ::  EPS = epsilon(UM)

  uz = UM / z
  zm = abs(z)
  zr = abs(real(z,REAL64))
  zi = abs(aimag(z))

! finds starting nup, by and byp

  if(zm < 5*n .or. (n > 15 .and. zr < 5*n + zi)) then

    if(abs(z) < 10*UM) then
      nup = nint(40*sqrt(abs(z)) / sqrt(10.0)) + 10
    elseif(abs(z) < 20*UM / (2 - exp(UM)/2)) then
      nup = nint(2*abs(z)+30)
    else
      nup = nint(exp(UM)*abs(z)/2 + 50)
   endif

    if(nup < n+10) nup = n+10

    call zsbessj_pow(nup+1, z, byp, acc)
    call zsbessj_pow(nup, z, by, acc)

  else

    nup = n + 10

    call zsbessj_fwd(nup+1, z, byp, acc)
    call zsbessj_fwd(nup, z, by, acc)

  endif

! recursion formula

  do j = nup,1,-1
    bym = (2*j+1)*uz*by - byp
    byp = by
    by = bym
    if(j-1 == n) zsb = by
  enddo

! renormalizes from the final value of j_0 or j_1

  zsb0 = sin(z) * uz
  zsb1 = (zsb0-cos(z)) * uz
  renorm0 = zsb0/by
  renorm1 = zsb1/byp
  if(abs(zsb0) > abs(zsb1)) then
    zsb = zsb * renorm0
    acc = -log10(abs(abs(renorm0)-UM)+EPS)
  else
    zsb = zsb * renorm1
    acc = -log10(abs(abs(renorm1)-UM)+EPS)
  endif

  return

end subroutine zsbessj_bwd





!>  Calculates the spherical bessel function of first kind j_n(z)
!>  for a complex argument with a power series around z = 0.
!>  Formula 10.1.2 of Abramowitz and Stegun.
!>
!>  It provides an (under?) estimate of the number of accurate digits.
!>
!>  \author       Jose Luis Martins
!>  \version      0.1
!>  \date         19 April 2018, 27 March 2025.
!>  \copyright    GNU LGPL v3

subroutine zsbessj_pow(n, z, zsb, acc)

!  Code basis was adapted on 19 April 2018 from old atom code, Siesta code, and NumSBF code.

!  Relative errors are relevant for Abs(Re(z)) ~ 0.75 n, n > 30 and
!  Abs(Re(z)) >> Abs(Im(z))

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                  ::  n                             !<  n >= 0 order of function
  complex(REAL64), intent(in)          ::  z                             !<  argument

! output

  complex(REAL64), intent(out)         ::  zsb                           !<  result
  real(REAL64) , intent(out)           ::  acc                           !<  (under)estimation of number of accurate digits)

! local variables

  complex(REAL64)               ::  pref                                 !  prefactor of series expansion
  complex(REAL64)               ::  z2                                   !  0.5*z^2
  complex(REAL64)               ::  sumz, fac                            !  series expansion
  logical                       ::  fail
  real(REAL64)                  ::  facmax
  integer                       ::  jmax, kmax

! counter

  integer                    ::  j, k

! parameters

  real(REAL64), parameter           ::  UM = 1.0_REAL64
  real(REAL64), parameter           ::  EPS = epsilon(UM)


! Series expansion.  Reordering of products delays onset of overflow.

  pref = UM
  if(n > 0) then
!     do j = 1,n
!       pref = pref*z/(2*j+1)
!     enddo
    do j = 1,n/2
      pref = pref*z/(2*j+1)
      pref = pref*z/(2*(n-j)+3)
    enddo
    if(mod(n,2) == 1)  pref = pref*z/(n+2)
  endif

! maximum j can be estimated from Stirling formula, but it is simpler
! to use a overestimate.

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
       kmax = j

       exit

     endif
  enddo

  zsb = sumz*pref

! The improvement of the next few lines is very very minor.
! They may be commented for efficiency.

  kmax = kmax + 1
  do k = kmax,1,-1
    fac = UM - z2 * fac / ((k * (2*n + 2*k +1))*UM)
  enddo

  zsb = fac*pref

  acc = facmax / (abs(sumz)+EPS)
  acc = -log10(acc*EPS)

  if(fail) acc = 0

  return

end subroutine zsbessj_pow




!>  Calculates the spherical bessel function of first kind j_n(z)
!>  for a complex argument with a forward series from n = 0,1.
!>  Formula 10.1.19 of Abramowitz and Stegun.
!>
!>  It provides an estimate of the number of accurate digits.
!>  This estimate may be optimist for large values of z.
!>
!>  \author       Jose Luis Martins
!>  \version      0.1
!>  \date         19 April 2018, 27 March 2025.
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




!>  Calculates the spherical bessel function of first kind j_n(z)
!>  for a complex argument with an asymptotic series.
!>  Formula 10.1.8 of Abramowitz and Stegun.
!>
!>  Can be used for Re(z) >> n
!>
!>  A guess of accuracy is provided.
!>
!>  In the asymptotic region the forward formula is usualy more
!>  accurate.  It is not used in zsbessj, it is included in the
!>  eventuality it may be useful...
!>
!>  \author       Jose Luis Martins
!>  \version      0.1
!>  \date         11 February 2025. 27 March 2025.
!>  \copyright    GNU LGPL v3


subroutine zsbessj_asy(n, z, zsb, acc)



  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                  ::  n                             !<  n >= 0 order of function
  complex(REAL64), intent(in)          ::  z                             !<  argument

! output

  complex(REAL64), intent(out)         ::  zsb                           !<  result
  real(REAL64) , intent(out)           ::  acc                           !<  heuristic estimation of number of accurate digits

! local variables

  complex(REAL64)               ::  pz, qz                               !  prefactor of series expansion
  complex(REAL64)               ::  uz                                   !  1/z
  complex(REAL64)               ::  fac
  integer                       ::  kmod, nmod

! counter

  integer                    ::  k

! parameters

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64
  real(REAL64), parameter           ::  UM = 1.0_REAL64
  complex(REAL64), parameter        ::  C_UM = cmplx(UM,ZERO,REAL64)
  complex(REAL64), parameter        ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  real(REAL64), parameter           ::  EPS = epsilon(UM)


  uz = UM / z

! recursion formula

  if(n == 0) then
    zsb = sin(z) * uz
    acc = -log10(EPS)
  elseif(n == 1) then
    zsb = (sin(z)*uz - cos(z)) * uz
    acc = -log10(abs(sin(z)*uz - cos(z))*EPS)
  else

    pz = C_UM
    qz = C_ZERO
    fac = C_UM

    do k = 1,n

      fac = uz*(fac*(n+k)*(n-k+1)) / (2*k)

      kmod = mod(k,4)

      if(kmod == 0) then
        pz = pz + fac
      elseif(kmod == 1) then
        qz = qz + fac
      elseif(kmod == 2) then
        pz = pz - fac
      elseif(kmod == 3) then
        qz = qz - fac
      endif

    enddo

    nmod = mod(n,4)
    if(nmod == 0) then
      zsb = uz * ( pz*sin(z) + qz*cos(z) )
    elseif(nmod == 1) then
      zsb =-uz * ( pz*cos(z) - qz*sin(z) )
    elseif(nmod == 2) then
      zsb =-uz * ( pz*sin(z) + qz*cos(z) )
    elseif(nmod == 3) then
      zsb = uz * ( pz*cos(z) - qz*sin(z) )
    endif
    acc = (n+1)*n*abs(uz)/2

!   heuristic guess

    acc = 15.5 - acc/10 - acc*acc/100

  endif

  return

end subroutine zsbessj_asy

