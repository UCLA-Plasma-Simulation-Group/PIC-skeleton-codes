!-----------------------------------------------------------------------
! Fortran2003 interface file for kncmbpush3.c
      module kncmbpush3_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine ckncgbppush3lt(s_ppart,s_fxyz,s_bxyz,s_kpic,qbm,dt, &
     &dtc,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,   &
     &ipbc) bind(C,name='ckncgbppush3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real(c_float), value :: qbm, dt, dtc
         real(c_float) :: ek
         type (c_ptr), value :: s_ppart, s_fxyz, s_bxyz, s_kpic
         end subroutine
      end interface
!
      interface
         subroutine ckncgbppushf3lt(s_ppart,s_fxyz,s_bxyz,s_kpic,s_ncl, &
     &s_ihole,qbm,dt,dtc,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,  &
     &mx1,my1,mxyz1,ntmax,irc) bind(C,name='ckncgbppushf3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer(c_int) :: irc
         real(c_float), value :: qbm, dt, dtc
         real(c_float) :: ek
         type (c_ptr), value :: s_ppart, s_fxyz, s_bxyz, s_kpic
         type (c_ptr), value :: s_ncl, s_ihole
         end subroutine
      end interface
!
      interface
         subroutine ckncgrbppush3lt(s_ppart,s_fxyz,s_bxyz,s_kpic,qbm,dt,&
     &dtc,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,&
     ipbc) bind(C,name='ckncgrbppush3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real(c_float), value :: qbm, dt, dtc, ci
         real(c_float) :: ek
         type (c_ptr), value :: s_ppart, s_fxyz, s_bxyz, s_kpic
         end subroutine
      end interface
!
      interface
         subroutine ckncgrbppushf3lt(s_ppart,s_fxyz,s_bxyz,s_kpic,s_ncl,&
     &s_ihole,qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv&
     &,mx1,my1,mxyz1,ntmax,irc) bind(C,name='ckncgrbppushf3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer(c_int) :: irc
         real(c_float), value :: qbm, dt, dtc, ci
         real(c_float) :: ek
         type (c_ptr), value :: s_ppart, s_fxyz, s_bxyz, s_kpic
         type (c_ptr), value :: s_ncl, s_ihole
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
         subroutine ckncgjppost3lt(s_ppart,s_cu,s_kpic,qm,dt,nppmx,idimp&
     &,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)                &
     &bind(C,name='ckncgjppost3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real(c_float), value :: qm, dt
         type (c_ptr), value :: s_ppart, s_cu, s_kpic
         end subroutine
      end interface
!
      interface
         subroutine ckncgjppostf3lt(s_ppart,s_cu,s_kpic,s_ncl,s_ihole,qm&
     &,dt,nppmx,idimp,nx,ny,nz,mx, my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax&
     &,irc) bind(C,name='ckncgjppostf3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer(c_int) :: irc
         real(c_float), value :: qm, dt
         type (c_ptr), value :: s_ppart, s_cu, s_kpic, s_ncl, s_ihole
         end subroutine
      end interface
!
      interface
         subroutine ckncgrjppost3lt(s_ppart,s_cu,s_kpic,qm,dt,ci,nppmx, &
     &idimp,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)           &
     &bind(C,name='ckncgrjppost3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real(c_float), value :: qm, dt, ci
         type (c_ptr), value :: s_ppart, s_cu, s_kpic
         end subroutine
      end interface
!
      interface
         subroutine ckncgrjppostf3lt(s_ppart,s_cu,s_kpic,s_ncl,s_ihole, &
     &qm,dt,ci,nppmx,idimp,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1, &
     &ntmax,irc) bind(C,name='ckncgrjppostf3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer(c_int) :: irc
         real(c_float), value :: qm, dt, ci
         type (c_ptr), value :: s_ppart, s_cu, s_kpic, s_ncl, s_ihole
         end subroutine
      end interface
!
      interface
         subroutine cknc2gjppost3lt(s_ppart,s_cu,s_kpic,qm,dt,nppmx,    &
     &idimp,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)           &
     &bind(C,name='cknc2gjppost3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real(c_float), value :: qm, dt
         type (c_ptr), value :: s_ppart, s_cu, s_kpic
         end subroutine
      end interface
!
      interface
         subroutine cknc2grjppost3lt(s_ppart,s_cu,s_kpic,qm,dt,ci,nppmx,&
     &idimp,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)           &
     &bind(C,name='cknc2grjppost3lt')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer(c_int), value :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real(c_float), value :: qm, dt, ci
         type (c_ptr), value :: s_ppart, s_cu, s_kpic
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
         subroutine ckncacguard3l(s_cu,nx,ny,nz,nxe,nye,nze)            &
     &bind(C,name='ckncacguard3l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nz, nxe, nye, nze
         type (c_ptr), value :: s_cu
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
         subroutine ckncmcuperp3(s_cu,nx,ny,nz,nxvh,nyv,nzv)            &
        &bind(C,name='ckncmcuperp3')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nz, nxvh, nyv, nzv
         type (c_ptr), value :: s_cu
         end subroutine
      end interface
!
      interface
         subroutine ckncmibpois33(s_cu,s_bxyz,s_ffc,ci,wm,nx,ny,nz,nxvh,&
     &nyv,nzv,nxhd,nyhd,nzhd) bind(C,name='ckncmibpois33')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nz, nxvh, nyv, nzv
         integer(c_int), value :: nxhd, nyhd, nzhd
         real(c_float), value :: ci
         real(c_float) :: wm
         type (c_ptr), value :: s_cu, s_bxyz, s_ffc
         end subroutine
      end interface
!
      interface
         subroutine ckncmmaxwel3(s_exyz,s_bxyz,s_cu,s_ffc,ci,dt,wf,wm,nx&
     &,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd) bind(C,name='ckncmmaxwel3')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nz, nxvh, nyv, nzv
         integer(c_int), value :: nxhd, nyhd, nzhd
         real(c_float), value :: ci, dt
         real(c_float) :: wf, wm
         type (c_ptr), value :: s_exyz, s_bxyz, s_cu, s_ffc
         end subroutine
      end interface
!
      interface
         subroutine ckncmemfield3(s_fxyz,s_exyz,s_ffc,isign,nx,ny,nz,   &
     &nxvh,nyv,nzv,nxhd,nyhd,nzhd) bind(C,name='ckncmemfield3')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer(c_int), value :: nxhd, nyhd, nzhd
         type (c_ptr), value :: s_fxyz, s_exyz, s_ffc
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
