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

!>  Tests the complex spherical bessel functions
!>
!>  \author       Jose Luis Martins
!>  \version      0.01
!>  \date         February 2025.
!>  \copyright    GNU LGPL v3



program test_bessel

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: REAL128 = selected_real_kind(28)

  integer                     ::  n

  complex(REAL64)             ::  z

  complex(REAL64)             ::  zsb_p
  real(REAL64)                ::  acc_p

  complex(REAL128)            ::  z_q
  complex(REAL128)            ::  zsb_q
  real(REAL128)               ::  acc_q

  complex(REAL64)             ::  zsb_f
  real(REAL64)                ::  acc_f

  complex(REAL64)             ::  zsb_b
  real(REAL64)                ::  acc_b

  complex(REAL64)             ::  zsb
  complex(REAL64), external   ::  zsbessj

  real(REAL64)                ::  acc

!   real(REAL64)                ::  x,y
!   real(REAL64)                ::  relerror, abserror
!
!
!   real(REAL64), parameter     ::  EPS = epsilon(1.0_REAL64)
!
!   real(REAL64), parameter           ::  ZERO = 0.0_REAL64
!   real(REAL64), parameter           ::  UM = 1.0_REAL64
!   complex(REAL64), parameter        ::  C_UM = cmplx(UM,ZERO,REAL64)
!   complex(REAL64), parameter        ::  C_I = cmplx(ZERO,UM,REAL64)

  integer                     ::  i
  real(REAL64), parameter           ::  UM = 1.0_REAL64
  real(REAL64), parameter           ::  EPS = epsilon(UM)

  write(6,*)
  write(6,*)  '  This compares spherical Bessel function j_n(x)'
  write(6,*)

  do i = 1,100

    write(6,*) 'enter n,x'
    read(5,*) n,z
    write(6,*)

    z_q = z

    if(n < 0) stop

    call zsbessj_pow(n, z, zsb_p, acc_p)

    call zsbessj_fwd(n, z, zsb_f, acc_f)

    call zsbessj_bwd(n, z, zsb_b, acc_b)

    call zsbessj_quad(n, z_q, zsb_q, acc_q)

    zsb = zsbessj(n, z)

    if(acc_q < 17) then
      write(6,'(2g24.16,10x,f10.3,"    quad    ")') zsb_q, acc_q
      write(6,'(2g24.16                         )') zsb
      write(6,*)
      write(6,'(2g24.16,10x,f10.3,"    power   ")') zsb_p, acc_p
      write(6,'(2g24.16,10x,f10.3,"    forward ")') zsb_f, acc_f
      write(6,'(2g24.16,10x,f10.3,"    backward")') zsb_b, acc_b
      write(6,*)
    else
      write(6,'(2g24.16,10x,f10.3,"    quad    ")') zsb_q, acc_q
      acc = -log10(abs(zsb  -zsb_q)/(abs(zsb_q)+EPS)+EPS)
      write(6,'(2g24.16                         ,f10.3)') zsb, acc
      write(6,*)
      acc = -log10(abs(zsb_p-zsb_q)/(abs(zsb_q)+EPS)+EPS)
      write(6,'(2g24.16,10x,f10.3,"    power   ",f10.3)') zsb_p, acc_p, acc
      acc = -log10(abs(zsb_f-zsb_q)/(abs(zsb_q)+EPS)+EPS)
      write(6,'(2g24.16,10x,f10.3,"    forward ",f10.3)') zsb_f, acc_f, acc
      acc = -log10(abs(zsb_b-zsb_q)/(abs(zsb_q)+EPS)+EPS)
      write(6,'(2g24.16,10x,f10.3,"    backward",f10.3)') zsb_b, acc_b, acc
      write(6,*)
    endif

  enddo

  stop

end




