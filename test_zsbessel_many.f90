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

  complex(REAL128)            ::  z_q
  complex(REAL128)            ::  zsb_q
  real(REAL128)               ::  acc_q

  complex(REAL64)             ::  zsb
  complex(REAL64), external   ::  zsbessj

  real(REAL64)                ::  accabs
  real(REAL64)                ::  accrel

  real(REAL64)                ::  accabsmin
  complex(REAL64)             ::  zabsmin

  real(REAL64)                ::  accrelmin
  complex(REAL64)             ::  zrelmin

  real(REAL64)                ::  acut

  real(REAL64)                ::  alogmod

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

  integer                     ::  i, j

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64
  real(REAL64), parameter           ::  UM = 1.0_REAL64
  real(REAL64), parameter           ::  EPS = epsilon(UM)
  complex(REAL64), parameter        ::  C_UM = cmplx(UM,ZERO,REAL64)
  complex(REAL64), parameter        ::  C_I = cmplx(ZERO,UM,REAL64)

  write(6,*)
  write(6,*)  '  Testing spherical Bessel function j_n(z)'
  write(6,*)

! only prints for acc < acut

  acut = 14.0


  do n = 0,100

    accabsmin = 100.0
    accrelmin = 100.0

!    write(6,'("  testing n = ",i5)') n

    do i = 0,100
    do j = 0,100

      z = i*C_UM + j*C_I
      z = z / 1
      z_q = z

      call zsbessj_quad(n, z_q, zsb_q, acc_q)

      zsb = zsbessj(n, z)

      accrel = -log10(abs(zsb  -zsb_q)/(abs(zsb_q)+EPS)+EPS)
      accabs = -log10(abs(zsb  -zsb_q)+EPS)

      alogmod = log10(abs(zsb)+1000*tiny(UM))

      if(acc_q > 17) then
        if(accrelmin > accrel .and. alogmod > -4 ) then
          accrelmin = accrel
          zrelmin = z
        endif
        if(accabsmin > accabs .and. alogmod < 4) then
          accabsmin = accabs
          zabsmin = z
        endif
        if(accrel < acut .and. accabs < acut) then
          write(75,*)
          write(75,'("  n = ",i5,"     z = ",2f12.4,"    acc = ",2f10.2)')  n, z, accrel, accabs
          write(75,'("  j_n(z) = ",2g24.16)') zsb
        endif
        if(accrel < acut .and. alogmod > -4) then
          write(76,*)
          write(76,'("  n = ",i5,"     z = ",2f12.4,"    acc = ",f10.2)')  n, z, accrel
          write(76,'("  j_n(z) = ",2g24.16)') zsb
        endif
        if(accabs < acut .and. alogmod < -4) then
            write(77,*)
            write(77,'("  n = ",i5,"     z = ",2f12.4,"    acc = ",f10.2)')  n, z, accabs
            write(77,'("  j_n(z) = ",2g24.16)') zsb
        endif
      else
        write(78,'("  n = ",i5,"     z = ",2f12.4,"    acc = ",f10.2)')  n, z, acc_q
      endif

    enddo
    enddo

    write(6,*)
    write(6,'("  worst absolute accuracy for n = ",i5,"  z = ",2f10.4,"   acc = ",f10.4)')   &
                n, zabsmin, accabsmin
    write(6,'("  worst relative accuracy for n = ",i5,"  z = ",2f10.4,"   acc = ",f10.4)')   &
                n, zrelmin, accrelmin

  enddo

  stop

end




