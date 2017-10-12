!-----------------------------------------------------------------------
! Interface file for kncmbpush3.c
      module kncmbpush3_h
      implicit none
!
      interface
         subroutine ckncgbppush3lt(ppart,fxyz,bxyz,kpic,qbm,dt,dtc,ek,  &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(4,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine ckncgbppushf3lt(ppart,fxyz,bxyz,kpic,ncl,ihole,qbm, &
     &dt,dtc,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,&
     &ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(4,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine ckncgrbppush3lt(ppart,fxyz,bxyz,kpic,qbm,dt,dtc,ci, &
     &ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(4,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine ckncgrbppushf3lt(ppart,fxyz,bxyz,kpic,ncl,ihole,qbm,&
     &dt,dtc,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,   &
     &mxyz1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(4,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine ckncgppost3lt(ppart,q,kpic,qm,nppmx,idimp,mx,my,mz, &
     &nxv,nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, my, mz, nxv, nyv, nzv
         integer, intent(in) :: mx1, my1, mxyz1
         real, intent(in) :: qm
         real, dimension(nppmx,idimp,mxyz1), intent(in) :: ppart
         real, dimension(nxv,nyv,nzv), intent(inout) :: q
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cknc2gppost3lt(ppart,q,kpic,qm,nppmx,idimp,mx,my,mz,&
     &nxv,nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, my, mz, nxv, nyv, nzv
         integer, intent(in) :: mx1, my1, mxyz1
         real, intent(in) :: qm
         real, dimension(nppmx,idimp,mxyz1), intent(in) :: ppart
         real, dimension(nxv,nyv,nzv), intent(inout) :: q
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine ckncgjppost3lt(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny&
     &,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qm, dt
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(4,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine ckncgrjppost3lt(ppart,cu,kpic,qm,dt,ci,nppmx,idimp, &
     &nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine ckncpporder3lt(ppart,ppbuff,kpic,ncl,ihole,idimp,   &
     &nppmx,nx,ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: mx1, my1, mz1, npbmx, ntmax
         integer, intent(inout) :: irc
         real, dimension(nppmx,idimp,mx1*my1*mz1), intent(inout) ::     &
     &ppart
         real, dimension(npbmx,idimp,mx1*my1*mz1), intent(inout) ::     &
    &ppbuff
         integer, dimension(mx1*my1*mz1), intent(inout) :: kpic
         integer, dimension(26,mx1*my1*mz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine ckncpporderf3lt(ppart,ppbuff,kpic,ncl,ihole,idimp,  &
     &nppmx,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, my1, mz1, npbmx
         integer, intent(in) :: ntmax
         integer, intent(inout) :: irc
         real, dimension(nppmx,idimp,mx1*my1*mz1), intent(inout) ::     &
     &ppart
         real, dimension(npbmx,idimp,mx1*my1*mz1), intent(inout) ::     &
    &ppbuff
         integer, dimension(mx1*my1*mz1), intent(inout) :: kpic
         integer, dimension(26,mx1*my1*mz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1), intent(in) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine cknccguard3l(fxyz,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, dimension(4,nxe,nye,nze), intent(inout) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine ckncacguard3l(cu,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, dimension(4,nxe,nye,nye), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine ckncaguard3l(q,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, dimension(nxe,nye,nze), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine ckncmpois33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,&
     &nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ax, ay, az, affp
         real, intent(inout) :: we
         real, dimension(2*nxvh,nyv,nzv), intent(in) :: q
         real, dimension(4,2*nxvh,nyv,nzv), intent(inout) :: fxyz
         complex, dimension(nxhd,nyhd,nzhd), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cknccuperp3(cu,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         real, dimension(4,2*nxvh,nyv,nzv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine ckncibpois33(cu,bxyz,ffc,ci,wm,nx,ny,nz,nxvh,nyv,nzv&
     &,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(4,2*nxvh,nyv,nzv), intent(in) :: cu
         complex, dimension(4,nxvh,nyv,nzv), intent(inout) :: bxyz
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ckncmaxwel3(exyz,bxyz,cu,ffc,ci,dt,wf,wm,nx,ny,nz,  &
     &nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(4,nxvh,nyv,nzv), intent(inout) :: exyz, bxyz
         real, dimension(4,2*nxvh,nyv,nzv), intent(in) :: cu
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ckncemfield3(fxyz,exyz,ffc,isign,nx,ny,nz,nxvh,nyv, &
     &nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, dimension(4,2*nxvh,nyv,nzv), intent(inout) :: fxyz
         complex, dimension(4,nxvh,nyv,nzv), intent(in) :: exyz
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ckncwfft3rmx(f,isign,mixup,sct,indx,indy,indz,nxhd, &
     &nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer, intent(in) :: nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine ckncwfft3rm3(f,isign,mixup,sct,indx,indy,indz,nxhd, &
     &nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer, intent(in) :: nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
          subroutine ckncfft3rmxy(f,isign,mixup,sct,indx,indy,indz,nzi, &
     &nzp,nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine ckncfft3rmz(f,isign,mixup,sct,indx,indy,indz,nyi,nyp&
     &,nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine ckncfft3rm3xy(f,isign,mixup,sct,indx,indy,indz,nzi, &
     &nzp,nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine ckncfft3rm3z(f,isign,mixup,sct,indx,indy,indz,nyi,  &
     &nyp,nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer , intent(in):: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      end module
