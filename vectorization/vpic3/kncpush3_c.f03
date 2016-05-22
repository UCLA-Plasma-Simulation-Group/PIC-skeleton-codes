!-----------------------------------------------------------------------
! Fortran2003 interface file for kncpush3.c
      module kncpush3_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine ckncxiscan2(isdata,nths) bind(C,name='ckncxiscan2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nths
         type (c_ptr), value :: isdata
         end subroutine
      end interface
!
      interface
         subroutine ckncgpush3lt(s_part,s_fxyz,qbm,dt,s_ek,idimp,nop,npe&
     &,nx,ny,nz,nxv,nyv,nzv,ipbc) bind(C,name='ckncgpush3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nop, npe, nx, ny, nz
         integer(c_int), value :: nxv, nyv, nzv, ipbc
         real(c_float), value :: qbm, dt
         type (c_ptr), value :: s_part, s_fxyz, s_ek
         end subroutine
      end interface
!
      interface
         subroutine ckncgpost3lt(s_part,s_q,qm,nop,npe,idimp,nxv,nyv,nzv&
     &) bind(C,name='ckncgpost3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nop, npe, idimp, nxv, nyv, nzv
         real(c_float), value :: qm
         type (c_ptr), value :: s_part, s_q
         end subroutine
      end interface
!
      interface
         subroutine cknc2gpost3lt(s_part,s_q,qm,nop,npe,idimp,nxv,nyv,  &
     &nzv) bind(C,name='cknc2gpost3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nop, npe, idimp, nxv, nyv, nzv
         real(c_float), value :: qm
         type (c_ptr), value :: s_part, s_q
         end subroutine
      end interface
!
      interface
         subroutine ckncdsortp3yzlt(s_parta,s_partb,s_npic,idimp,nop,npe&
     &,ny1,nyz1) bind(C,name='ckncdsortp3yzlt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nop, npe, ny1, nyz1
         type (c_ptr), value :: s_parta, s_partb, s_npic
         end subroutine
      end interface
!
      interface
         subroutine cknccguard3l(s_fxyz,nx,ny,nz,nxe,nye,nze)           &
     &bind(C,name='cknccguard3l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nz, nxe, nye, nze
         type (c_ptr), value :: s_fxyz
         end subroutine
      end interface
!
      interface
         subroutine ckncaguard3l(s_q,nx,ny,nz,nxe,nye,nze)              &
     &bind(C,name='ckncaguard3l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nz, nxe, nye, nze
         type (c_ptr), value :: s_q
         end subroutine
      end interface
!
      interface
         subroutine ckncpois33(s_q,s_fxyz,isign,s_ffc,ax,ay,az,affp,s_we&
     &,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd) bind(C,name='ckncpois33')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer(c_int), value :: nxhd, nyhd, nzhd
         real(c_float), value :: ax, ay, az, affp
         type (c_ptr), value :: s_q, s_fxyz, s_ffc, s_we
         end subroutine
      end interface
!
      interface
         subroutine ckncwfft3rvx(s_f,isign,mixup,sct,indx,indy,indz,nxhd&
     &,nyd,nzd,nxhyzd,nxyzhd) bind(C,name='ckncwfft3rvx')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, indz, nxhd, nyd
         integer(c_int), value :: nzd, nxhyzd, nxyzhd
         type (c_ptr), value :: s_f, mixup, sct
         end subroutine
      end interface
!
      interface
         subroutine ckncwfft3rv3(s_f,isign,mixup,sct,indx,indy,indz,nxhd&
     &,nyd,nzd,nxhyzd,nxyzhd) bind(C,name='ckncwfft3rv3')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, indz, nxhd, nyd
         integer(c_int), value :: nzd, nxhyzd, nxyzhd
         type (c_ptr), value :: s_f, mixup, sct
         end subroutine
      end interface
!
      end module
