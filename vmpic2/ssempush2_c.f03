!-----------------------------------------------------------------------
! Fortran2003 interface file for ssempush2.c
      module ssempush2_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine csse2gppush2lt(s_ppart,s_fxy,s_kpic,qbm,dt,s_ek,    &
     &idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)                    &
     &bind(C,name='csse2gppush2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my
         integer(c_int), value :: nxv, nyv, mx1, mxy1, ipbc
         real(c_float), value :: qbm, dt
         type (c_ptr), value :: s_ppart, s_fxy, s_kpic, s_ek
         end subroutine
      end interface
!
      interface
         subroutine csse2gppushf2lt(s_ppart,s_fxy,s_kpic,s_ncl,s_ihole, &
     &qbm,dt,s_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)   &
     &bind(C,name='csse2gppushf2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ntmax
         integer(c_int) :: irc
         real(c_float), value :: qbm, dt
         type (c_ptr), value :: s_ppart, s_fxy, s_kpic, s_ncl, s_ihole
         type (c_ptr), value :: s_ek
         end subroutine
      end interface
!
      interface
         subroutine csse2gppost2lt(s_ppart,s_q,s_kpic,qm,nppmx,idimp,mx,&
     &my,nxv,nyv,mx1,mxy1) bind(C,name='csse2gppost2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1
         real(c_float), value :: qm
         type (c_ptr), value :: s_ppart, s_q, s_kpic
         end subroutine
      end interface
!
      interface
         subroutine csse2pporder2lt(s_ppart,s_ppbuff,s_kpic,s_ncl,      &
     &s_ihole,idimp,nppmx,nx,ny,mx,my,mx1,my1,npbmx,ntmax,irc)        &
     &bind(C,name='csse2pporder2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my
         integer(c_int), value :: mx1, my1, npbmx, ntmax
         integer(c_int) :: irc
         type (c_ptr), value :: s_ppart, s_ppbuff, s_kpic
         type (c_ptr), value :: s_ncl, s_ihole
         end subroutine
      end interface
!
      interface
         subroutine csse2pporderf2lt(s_ppart,s_ppbuff,s_kpic,s_ncl,     &
     &s_ihole,idimp,nppmx,mx1,my1,npbmx,ntmax,irc)                    &
     &bind(C,name='csse2pporderf2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, mx1, my1, npbmx, ntmax
         integer(c_int) :: irc
         type (c_ptr), value :: s_ppart, s_ppbuff, s_kpic
         type (c_ptr), value :: s_ncl, s_ihole
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
         subroutine csse2mpois22(s_q,s_fxy,isign,s_ffc,ax,ay,affp,s_we, &
     &nx,ny,nxvh,nyv,nxhd,nyhd) bind(C,name='csse2mpois22')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real(c_float), value :: ax, ay, affp
         type (c_ptr), value :: s_q, s_fxy, s_ffc, s_we
         end subroutine
      end interface
!
      interface
         subroutine csse2wfft2rmx(s_f,isign,mixup,sct,indx,indy,nxhd,nyd&
     &,nxhyd,nxyhd) bind(C,name='csse2wfft2rmx')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, nxhd, nyd, nxhyd
         integer(c_int), value :: nxyhd
         type (c_ptr), value :: s_f, mixup, sct
         end subroutine
      end interface
!
      interface
         subroutine csse2wfft2rm2(s_f,isign,mixup,sct,indx,indy,nxhd,nyd&
     &,nxhyd,nxyhd) bind(C,name='csse2wfft2rm2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, nxhd, nyd, nxhyd
         integer(c_int), value :: nxyhd
         type (c_ptr), value :: s_f, mixup, sct
         end subroutine
      end interface
!
      end module
