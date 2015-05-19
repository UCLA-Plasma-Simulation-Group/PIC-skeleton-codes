!-----------------------------------------------------------------------
! Fortran2003 interface file for ssebpush2.c
      module ssebpush2_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine csse2xiscan2(isdata,nths)                           &
     &bind(C,name='csse2xiscan2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nths
         type (c_ptr), value :: isdata
         end subroutine
      end interface
!
      interface
         subroutine csse2gbpush23lt(s_part,s_fxy,s_bxy,qbm,dt,dtc,s_ek, &
     &idimp,nop,npe,nx,ny,nxv,nyv,ipbc) bind(C,name='csse2gbpush23lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nop, npe, nx, ny, nxv, nyv
         integer(c_int), value :: ipbc
         real(c_float), value :: qbm, dt, dtc
         type (c_ptr), value :: s_part, s_fxy, s_bxy, s_ek
         end subroutine
      end interface
!
      interface
         subroutine csse2grbpush23lt(s_part,s_fxy,s_bxy,qbm,dt,dtc,ci,  &
     &s_ek,idimp,nop,npe,nx,ny,nxv,nyv,ipbc)                            &
     &bind(C,name='csse2grbpush23lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nop, npe, nx, ny, nxv, nyv
         integer(c_int), value :: ipbc
         real(c_float), value :: qbm, dt, dtc, ci
         type (c_ptr), value :: s_part, s_fxy, s_bxy, s_ek
         end subroutine
      end interface
!
      interface
         subroutine csse2gpost2lt(s_part,s_q,qm,nop,npe,idimp,nxv,nyv)  &
     &bind(C,name='csse2gpost2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nop, npe, idimp, nxv, nyv
         real(c_float), value :: qm
         type (c_ptr), value :: s_part, s_q
         end subroutine
      end interface
!
      interface
         subroutine csse2gjpost2lt(s_part,s_cu,qm,dt,nop,npe,idimp,nx,ny&
     &,nxv,nyv,ipbc) bind(C,name='csse2gjpost2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nop, npe, idimp, nx, ny, nxv, nyv
         integer(c_int), value :: ipbc
         real(c_float), value :: qm, dt
         type (c_ptr), value :: s_part, s_cu
         end subroutine
      end interface
!
      interface
         subroutine csse2grjpost2lt(s_part,s_cu,qm,dt,ci,nop,npe,idimp, &
     &nx,ny,nxv,nyv,ipbc) bind(C,name='csse2grjpost2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nop, npe, idimp, nx, ny, nxv, nyv
         integer(c_int), value :: ipbc
         real(c_float), value :: qm, dt, ci
         type (c_ptr), value :: s_part, s_cu
         end subroutine
      end interface
!
      interface
         subroutine csse2dsortp2ylt(s_parta,s_partb,s_npic,idimp,nop,npe&
     &,ny1) bind(C,name='csse2dsortp2ylt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nop, npe, ny1
         type (c_ptr), value :: s_parta, s_partb, s_npic
         end subroutine
      end interface
!
      interface
         subroutine csse2bguard2l(s_bxy,nx,ny,nxe,nye)                  &
     &bind(C,name='csse2bguard2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxe, nye
         type (c_ptr), value :: s_bxy
         end subroutine
      end interface
!
      interface
         subroutine csse2acguard2l(s_cu,nx,ny,nxe,nye)                  &
     &bind(C,name='csse2acguard2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxe, nye
         type (c_ptr), value :: s_cu
         end subroutine
      end interface
!
      interface
         subroutine csse2aguard2l(s_q,nx,ny,nxe,nye)                    &
     &bind(C,name='csse2aguard2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxe, nye
         type (c_ptr), value :: s_q
         end subroutine
      end interface
!
      interface
         subroutine csse2pois23(s_q,s_fxy,isign,s_ffc,ax,ay,affp,s_we,nx&
     &,ny,nxvh,nyv,nxhd,nyhd) bind(C,name='csse2pois23')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real(c_float), value :: ax, ay, affp
         type (c_ptr), value :: s_q, s_fxy, s_ffc, s_we
         end subroutine
      end interface
!
      interface
         subroutine csse2cuperp2(s_cu,nx,ny,nxvh,nyv)                   &
     &bind(C,name='csse2cuperp2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxvh, nyv
         type (c_ptr), value :: s_cu
         end subroutine
      end interface
!
      interface
         subroutine csse2ibpois23(s_cu,s_bxy,s_ffc,ci,s_wm,nx,ny,nxvh,  &
     &nyv,nxhd,nyhd) bind(C,name='csse2ibpois23')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxvh, nyv, nxhd, nyhd
         real(c_float), value :: ci
         type (c_ptr), value :: s_cu, s_bxy, s_ffc, s_wm
         end subroutine
      end interface
!
      interface
         subroutine csse2maxwel2(s_exy,s_bxy,s_cu,s_ffc,ci,dt,s_wf,s_wm,&
     &nx,ny,nxvh,nyv,nxhd,nyhd) bind(C,name='csse2maxwel2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxvh, nyv, nxhd, nyhd
         real(c_float), value :: ci, dt
         type (c_ptr), value :: s_exy, s_bxy, s_cu, s_ffc, s_wf, s_wm
         end subroutine
      end interface
!
      interface
         subroutine csse2emfield2(s_fxy,s_exy,s_ffc,isign,nx,ny,nxvh,nyv&
     &,nxhd,nyhd) bind(C,name='csse2emfield2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         type (c_ptr), value :: s_fxy, s_exy, s_ffc
         end subroutine
      end interface
!
      interface
         subroutine csse2wfft2rx(s_f,isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd) bind(C,name='csse2wfft2rx')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, nxhd, nyd, nxhyd
          integer(c_int), value :: nxyhd
         type (c_ptr), value :: s_f, mixup, sct
         end subroutine
      end interface
!
      interface
         subroutine csse2wfft2r3(s_f,isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd) bind(C,name='csse2wfft2r3')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, nxhd, nyd, nxhyd
         integer(c_int), value :: nxyhd
         type (c_ptr), value :: s_f, mixup, sct
         end subroutine
      end interface
!
      end module
