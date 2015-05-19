!-----------------------------------------------------------------------
! Fortran2003 interface file for ssepush2.c
      module ssepush2_c
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
         subroutine csse2gpush2lt(s_part,s_fxy,qbm,dt,s_ek,idimp,nop,npe&
     &,nx,ny,nxv,nyv,ipbc) bind(C,name='csse2gpush2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nop, npe, nx, ny, nxv, nyv
         integer(c_int), value :: ipbc
         real(c_float), value :: qbm, dt
         type (c_ptr), value :: s_part, s_fxy, s_ek
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
         subroutine csse2cguard2l(s_fxy,nx,ny,nxe,nye)                  &
     &bind(C,name='csse2cguard2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxe, nye
         type (c_ptr), value :: s_fxy
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
         subroutine csse2pois22(s_q,s_fxy,isign,s_ffc,ax,ay,affp,s_we,nx&
     &,ny,nxvh,nyv,nxhd,nyhd) bind(C,name='csse2pois22')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real(c_float), value :: ax, ay, affp
         type (c_ptr), value :: s_q, s_fxy, s_ffc, s_we
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
         subroutine csse2wfft2r2(s_f,isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd) bind(C,name='csse2wfft2r2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, nxhd, nyd, nxhyd
         integer(c_int), value :: nxyhd
         type (c_ptr), value :: s_f, mixup, sct
         end subroutine
      end interface
!
      end module
