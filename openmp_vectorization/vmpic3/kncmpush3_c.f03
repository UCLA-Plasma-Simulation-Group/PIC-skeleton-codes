!-----------------------------------------------------------------------
! Fortran2003 interface file for kncmpush3.c
      module kncmpush3_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine ckncgppush3lt(s_ppart,s_fxyz,s_kpic,qbm,dt,ek,idimp,&
     &nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)           &
     &bind(C,name='ckncgppush3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real(c_float), value :: qbm, dt
         real(c_float) :: ek
         type (c_ptr), value :: s_ppart, s_fxyz, s_kpic
         end subroutine
      end interface
!
      interface
         subroutine ckncgppushf3lt(s_ppart,s_fxyz,s_kpic,s_ncl,s_ihole, &
     &qbm,dt,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,&
     &ntmax,irc) bind(C,name='ckncgppushf3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer(c_int) :: irc
         real(c_float), value :: qbm, dt
         real(c_float) :: ek
         type (c_ptr), value :: s_ppart, s_fxyz, s_kpic, s_ncl, s_ihole
         end subroutine
      end interface
!
      interface
         subroutine ckncgppost3lt(s_ppart,s_q,s_kpic,qm,nppmx,idimp,mx, &
     &my,mz,nxv,nyv,nzv,mx1,my1,mxyz1) bind(C,name='ckncgppost3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1
         real(c_float), value :: qm
         type (c_ptr), value :: s_ppart, s_q, s_kpic
         end subroutine
      end interface
!
      interface
         subroutine cknc2gppost3lt(s_ppart,s_q,s_kpic,qm,nppmx,idimp,mx,&
     &my,mz,nxv,nyv,nzv,mx1,my1,mxyz1) bind(C,name='cknc2gppost3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1
         real(c_float), value :: qm
         type (c_ptr), value :: s_ppart, s_q, s_kpic
         end subroutine
      end interface
!
      interface
         subroutine ckncpporder3lt(s_ppart,s_ppbuff,s_kpic,s_ncl,s_ihole&
     &,idimp,nppmx,nx,ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)       &
     &bind(C,name='ckncpporder3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer(c_int), value :: mx1, my1, mz1, npbmx, ntmax
         integer(c_int) :: irc
         type (c_ptr), value :: s_ppart, s_ppbuff
         type (c_ptr), value :: s_kpic, s_ncl, s_ihole
         end subroutine
      end interface
!
      interface
         subroutine ckncpporderf3lt(s_ppart,s_ppbuff,s_kpic,s_ncl,      &
     &s_ihole,idimp,nppmx,mx1,my1,mz1,npbmx,ntmax,irc)                  &
     &bind(C,name='ckncpporderf3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, mx1, my1, mz1
         integer(c_int), value :: npbmx, ntmax
         integer(c_int) :: irc
         type (c_ptr), value :: s_ppart, s_ppbuff
         type (c_ptr), value :: s_kpic, s_ncl, s_ihole
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
         subroutine ckncmpois33(s_q,s_fxyz,isign,s_ffc,ax,ay,az,affp,we,&
     &nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd) bind(C,name='ckncmpois33')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer(c_int), value :: nxhd, nyhd, nzhd
         real(c_float), value :: ax, ay, az, affp
         real(c_float) :: we
         type (c_ptr), value :: s_q, s_fxyz, s_ffc
         end subroutine
      end interface
!
      interface
         subroutine ckncwfft3rmx(s_f,isign,s_mixup,s_sct,indx,indy,indz,&
     &nxhd,nyd,nzd,nxhyzd,nxyzhd) bind(C,name='ckncwfft3rmx')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, indz, nxhd, nyd
         integer(c_int), value :: nzd, nxhyzd, nxyzhd
         type (c_ptr), value :: s_f, s_mixup, s_sct
         end subroutine
      end interface
!
      interface
         subroutine ckncwfft3rm3(s_f,isign,s_mixup,s_sct,indx,indy,indz,&
     &nxhd,nyd,nzd,nxhyzd,nxyzhd) bind(C,name='ckncwfft3rm3')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, indz, nxhd, nyd
         integer(c_int), value :: nzd, nxhyzd, nxyzhd
         type (c_ptr), value :: s_f, s_mixup, s_sct
         end subroutine
      end interface
!
      end module
