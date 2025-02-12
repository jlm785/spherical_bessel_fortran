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

    call zsbessj_bwd(n, z, zsb_b)

    call zsbessj_pow_quad(n, z_q, zsb_q, acc_q)

    zsb = zsbessj(n, z)

    if(acc_q < 20) then
      write(6,'(2g24.16,10x,f10.3,"    quad    ")') zsb_q, acc_q
      write(6,'(2g24.16                         )') zsb
      write(6,*)
      write(6,'(2g24.16,10x,f10.3,"    power   ")') zsb_p, acc_p
      write(6,'(2g24.16,10x,f10.3,"    forward ")') zsb_f, acc_f
      write(6,'(2g24.16,20x      ,"    backward")') zsb_b
      write(6,*)
    else
      write(6,'(2g24.16,10x,f10.3,"    quad    ")') zsb_q, acc_q
      acc = -log10(abs(zsb  -zsb_q)/(abs(zsb_q)+EPS*abs(zsb  -zsb_q)))
      write(6,'(2g24.16                         ,f10.3)') zsb, acc
      write(6,*)
      acc = -log10(abs(zsb_p-zsb_q)/(abs(zsb_q)+EPS*abs(zsb_p-zsb_q)))
      write(6,'(2g24.16,10x,f10.3,"    power   ",f10.3)') zsb_p, acc_p, acc
      acc = -log10(abs(zsb_f-zsb_q)/(abs(zsb_q)+EPS*abs(zsb_f-zsb_q)))
      write(6,'(2g24.16,10x,f10.3,"    forward ",f10.3)') zsb_f, acc_f, acc
      acc = -log10(abs(zsb_b-zsb_q)/(abs(zsb_q)+EPS*abs(zsb_b-zsb_q)))
      write(6,'(2g24.16,20x      ,"    backward",f10.3)') zsb_b, acc
      write(6,*)
    endif

  enddo

  stop

end


