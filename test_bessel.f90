       program test_bessel

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

       integer                     ::  n
       real(REAL64)                ::  x
       real(REAL64),allocatable    ::  sj(:),dj(:)
       integer                     ::  nm
       real(REAL64)                ::  y,z
       real(REAL64), external      ::  sbessj

       integer                     ::  i

       write(6,*)
       write(6,*)  '  This program compares 2 two codes for the'
       write(6,*)  '  spherical Bessel function j_n(x)'
       write(6,*)

       do i=1,100        
         write(6,*) 'enter n,x'
         read(5,*) n,x
         allocate(sj(0:n),dj(0:n))
         call sphj(n,x,nm,sj,dj)
         y = sj(n)
         z = sbessj(n,x)
         write(6,'(8x,"sbessj",22x,"sphj",23x,"difference",10x,          &
     &             "relative difference")')
         write(6,'(g24.15,2x,g24.15,2x,f24.15,2x,f24.15)') z,y,z-y,      &
     &        abs(z-y)/max(abs(z),abs(y))
         deallocate(sj,dj)
       enddo

       stop
       end

subroutine sphj ( n, x, nm, sj, dj )

!*****************************************************************************80
!
!! SPHJ computes spherical Bessel functions jn(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    12 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin.
!    Modifications suggested by Vincent Lagage, 12 January 2016.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, real ( kind = 8 ) SJ(0:N), the values of jn(x).
!
!    Output, real ( kind = 8 ) DJ(0:N), the values of jn'(x).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cs
  real ( kind = 8 ) dj(0:n)
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) msta1
  integer ( kind = 4 ) msta2
  integer ( kind = 4 ) nm
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) sj(0:n)
  real ( kind = 8 ) x

  nm = n
!
!  Original code.
!
  if ( .true. ) then

    if ( abs ( x ) <= 1.0D-100 ) then
      do k = 0, n
        sj(k) = 0.0D+00
        dj(k) = 0.0D+00
      end do
      sj(0) = 1.0D+00
      dj(1) = 0.3333333333333333D+00
      return
    end if
!
!  Updated code.
!
  else

    if ( abs ( x ) <= 1.0D-16 ) then
      do k = 0, n
        sj(k) = 0.0D+00
        dj(k) = 0.0D+00
      end do
      sj(0) = 1.0D+00
      if ( 0 < n ) then
        do k = 1, n
          sj(k) = sj(k-1) * x / real ( 2 * k + 1, kind = 8 )
        end do
        dj(1) = 1.0D+00 / 3.0D+00
      end if
      return
    end if

  end if

  sj(0) = sin ( x ) / x
  sj(1) = ( sj(0) - cos ( x ) ) / x

  if ( 2 <= n ) then

    sa = sj(0)
    sb = sj(1)
    m = msta1 ( x, 200 )
    if ( m < n ) then
      nm = m
    else
      m = msta2 ( x, n, 15 )
    end if

    f0 = 0.0D+00
    f1 = 1.0D+00-100
    do k = m, 0, -1
      f = ( 2.0D+00 * k + 3.0D+00 ) * f1 / x - f0
      if ( k <= nm ) then
        sj(k) = f
      end if
      f0 = f1
      f1 = f
    end do

    if ( abs ( sa ) <= abs ( sb ) ) then
      cs = sb / f0
    else
      cs = sa / f
    end if

    do k = 0, nm
      sj(k) = cs * sj(k)
    end do

  end if      

  dj(0) = ( cos(x) - sin(x) / x ) / x
  do k = 1, nm
    dj(k) = sj(k-1) - ( k + 1.0D+00 ) * sj(k) / x
  end do

  return
end
function msta1 ( x, mp )

!*****************************************************************************80
!
!! MSTA1 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for backward  
!    recurrence such that the magnitude of    
!    Jn(x) at that point is about 10^(-MP).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) MP, the negative logarithm of the 
!    desired magnitude.
!
!    Output, integer ( kind = 4 ) MSTA1, the starting point.
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) envj
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  integer ( kind = 4 ) it
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) msta1
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nn
  real ( kind = 8 ) x

  a0 = abs ( x )
  n0 = int ( 1.1D+00 * a0 ) + 1
  f0 = envj ( n0, a0 ) - mp
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - mp
  do it = 1, 20       
    nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )                  
    f = envj ( nn, a0 ) - mp
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta1 = nn

  return
end
function msta2 ( x, n, mp )

!*****************************************************************************80
!
!! MSTA2 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for a backward
!    recurrence such that all Jn(x) has MP significant digits.
!
!    Jianming Jin supplied a modification to this code on 12 January 2016.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of Jn(x).
!
!    Input, integer ( kind = 4 ) N, the order of Jn(x).
!
!    Input, integer ( kind = 4 ) MP, the number of significant digits.
!
!    Output, integer ( kind = 4 ) MSTA2, the starting point.
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) ejn
  real ( kind = 8 ) envj
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) hmp
  integer ( kind = 4 ) it
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) msta2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nn
  real ( kind = 8 ) obj
  real ( kind = 8 ) x

  a0 = abs ( x )
  hmp = 0.5D+00 * mp
  ejn = envj ( n, a0 )

  if ( ejn <= hmp ) then
    obj = mp
!
!  Original code:
!
!   n0 = int ( 1.1D+00 * a0 )
!
!  Updated code:
!
    n0 = int ( 1.1D+00 * a0 ) + 1
  else
    obj = hmp + ejn
    n0 = n
  end if

  f0 = envj ( n0, a0 ) - obj
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - obj

  do it = 1, 20
    nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )
    f = envj ( nn, a0 ) - obj
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta2 = nn + 10

  return
end
function envj ( n, x )

!*****************************************************************************80
!
!! ENVJ is a utility function used by MSTA1 and MSTA2.
!
!  Discussion:
!
!    ENVJ estimates -log(Jn(x)) from the estimate
!    Jn(x) approx 1/sqrt(2*pi*n) * ( e*x/(2*n))^n
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!    Modifications suggested by Vincent Lafage, 11 January 2016.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Bessel function.
!
!    Input, real ( kind = 8 ) X, the absolute value of the argument.
!
!    Output, real ( kind = 8 ) ENVJ, the value.
!
  implicit none

  real ( kind = 8 ) envj
  real ( kind = 8 ) logten
  integer ( kind = 4 ) n
  real ( kind = 8 ) n_r8
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) x
!
!  Original code
!
  if ( .true. ) then

    envj = 0.5D+00 * log10 ( 6.28D+00 * n ) &
      - n * log10 ( 1.36D+00 * x / n )
!
!  Modification suggested by Vincent Lafage.
!
  else

    n_r8 = real ( n, kind = 8 )
    logten = log ( 10.0D+00 )
    envj = r8_gamma_log ( n_r8 + 1.0D+00 ) / logten - n_r8 * log10 ( x )

  end if

  return
end
 
