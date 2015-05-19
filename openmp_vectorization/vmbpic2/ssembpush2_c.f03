!-----------------------------------------------------------------------
! Fortran2003 interface file for ssembpush2.c
      module ssembpush2_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine csse2gbppush23lt(s_ppart,s_fxy,s_bxy,s_kpic,qbm,dt, &
     &dtc,s_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)           &
     &bind(C,name='csse2gbppush23lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ipbc
         real(c_float), value :: qbm, dt, dtc
         type (c_ptr), value :: s_ppart, s_fxy, s_bxy, s_kpic, s_ek
         end subroutine
      end interface
!
      interface
         subroutine csse2gbppushf23lt(s_ppart,s_fxy,s_bxy,s_kpic,s_ncl, &
     &s_ihole,qbm,dt,dtc,s_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1, &
     &ntmax,irc) bind(C,name='csse2gbppushf23lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ntmax
         integer(c_int) :: irc
         real(c_float), value :: qbm, dt, dtc
         type (c_ptr), value :: s_ppart, s_fxy, s_bxy, s_kpic, s_ncl
         type (c_ptr), value :: s_ihole, s_ek
         end subroutine
      end interface
!
      interface
         subroutine csse2grbppush23lt(s_ppart,s_fxy,s_bxy,s_kpic,qbm,dt,&
     &dtc,ci,s_ek,idimp, nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)       &
     &bind(C,name='csse2grbppush23lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ipbc
         real(c_float), value :: qbm, dt, dtc, ci
         type (c_ptr), value :: s_ppart, s_fxy, s_bxy, s_kpic, s_ek
         end subroutine
      end interface
!
      interface
         subroutine csse2grbppushf23lt(s_ppart,s_fxy,s_bxy,s_kpic,s_ncl,&
     &s_ihole,qbm,dt,dtc,ci,s_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,   &
     &mxy1,ntmax,irc) bind(C,name='csse2grbppushf23lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ntmax
         integer(c_int) :: irc
         real(c_float), value :: qbm, dt, dtc, ci
         type (c_ptr), value :: s_ppart, s_fxy, s_bxy, s_kpic, s_ncl
         type (c_ptr), value :: s_ihole, s_ek
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
         subroutine csse2gjppost2lt(s_ppart,s_cu,s_kpic,qm,dt,nppmx,    &
     &idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)                          &
     &bind(C,name='csse2gjppost2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ipbc
         real(c_float), value :: qm, dt
         type (c_ptr), value :: s_ppart, s_cu, s_kpic
         end subroutine
      end interface
!
      interface
         subroutine csse2gjppostf2lt(s_ppart,s_cu,s_kpic,s_ncl,s_ihole, &
     &qm,dt,nppmx,idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)         &
     &bind(C,name='csse2gjppostf2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ntmax
         integer(c_int) :: irc
         real(c_float), value :: qm, dt
         type (c_ptr), value :: s_ppart, s_cu, s_kpic, s_ncl, s_ihole
         end subroutine
      end interface
!
      interface
         subroutine csse2grjppost2lt(s_ppart,s_cu,s_kpic,qm,dt,ci,nppmx,&
     &idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)                          &
     &bind(C,name='csse2grjppost2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ipbc
         real(c_float), value :: qm, dt, ci
         type (c_ptr), value :: s_ppart, s_cu, s_kpic
         end subroutine
      end interface
!
      interface
         subroutine csse2grjppostf2lt(s_ppart,s_cu,s_kpic,s_ncl,s_ihole,&
     &qm,dt,ci,nppmx,idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)      &
     &bind(C,name='csse2grjppostf2lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer(c_int), value ::  mx1, mxy1, ntmax
         integer(c_int) :: irc
         real(c_float), value :: qm, dt, ci
         type (c_ptr), value :: s_ppart, s_cu, s_kpic, s_ncl, s_ihole
         end subroutine
      end interface
!
      interface
         subroutine csse2pporder2lt(s_ppart,s_ppbuff,s_kpic,s_ncl,      &
     &s_ihole,idimp,nppmx,nx,ny,mx,my,mx1,my1,npbmx,ntmax,irc)          &
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
     &s_ihole,idimp,nppmx,mx1,my1,npbmx,ntmax,irc)                      &
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
         subroutine csse2mpois23(s_q,s_fxy,isign,s_ffc,ax,ay,affp,s_we, &
     &nx,ny,nxvh,nyv,nxhd,nyhd) bind(C,name='csse2mpois23')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real(c_float), value :: ax, ay, affp
         type (c_ptr), value :: s_q, s_fxy, s_ffc, s_we
         end subroutine
      end interface
!
      interface
         subroutine csse2mcuperp2(s_cu,nx,ny,nxvh,nyv)                  &
     &bind(C,name='csse2mcuperp2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxvh, nyv
         type (c_ptr), value :: s_cu
         end subroutine
      end interface
!
      interface
         subroutine csse2mibpois23(s_cu,s_bxy,s_ffc,ci,s_wm,nx,ny,nxvh, &
     &nyv,nxhd,nyhd) bind(C,name='csse2mibpois23')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxvh, nyv, nxhd, nyhd
         real(c_float), value :: ci
         type (c_ptr), value :: s_cu, s_bxy, s_ffc, s_wm
         end subroutine
      end interface
!
      interface
         subroutine csse2mmaxwel2(s_exy,s_bxy,s_cu,s_ffc,ci,dt,s_wf,s_wm&
     &,nx,ny,nxvh,nyv,nxhd,nyhd) bind(C,name='csse2mmaxwel2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxvh, nyv, nxhd, nyhd
         real(c_float), value :: ci, dt
         type (c_ptr), value :: s_exy, s_bxy, s_cu, s_ffc, s_wf, s_wm
         end subroutine
      end interface
!
      interface
         subroutine csse2memfield2(s_fxy,s_exy,s_ffc,isign,nx,ny,nxvh,  &
     &nyv,nxhd,nyhd) bind(C,name='csse2memfield2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         type (c_ptr), value :: s_fxy, s_exy, s_ffc
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
      interface
         subroutine csse2wfft2rm3(s_f,isign,mixup,sct,indx,indy,nxhd,nyd&
     &,nxhyd,nxyhd) bind(C,name='csse2wfft2rm3')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, nxhd, nyd, nxhyd
         integer(c_int), value :: nxyhd
         type (c_ptr), value :: s_f, mixup, sct
         end subroutine
      end interface
!
      end module
