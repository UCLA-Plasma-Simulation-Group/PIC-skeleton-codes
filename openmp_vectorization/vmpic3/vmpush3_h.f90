!-----------------------------------------------------------------------
! Interface file for vmpush3.f
      module vmpush3_h
      implicit none
!
      interface
         subroutine DISTR3(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,    &
     &idimp,nop,nx,ny,nz,ipbc)
         implicit none
         integer :: npx, npy, npz, idimp, nop, nx, ny, nz, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
!
      interface
         subroutine DBLKP3L(part,kpic,nppmx,idimp,nop,mx,my,mz,mx1,my1, &
     &mxyz1,irc)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mz, mx1, my1, mxyz1, irc
         real, dimension(idimp,nop) :: part
         integer, dimension(mxyz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN3LT(part,ppart,kpic,nppmx,idimp,nop,mx,my,mz,&
     &mx1,my1,mxyz1,irc)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mz, mx1, my1, mxyz1, irc
         real, dimension(idimp,nop) :: part
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         integer, dimension(mxyz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN3LTP(part,ppart,kpic,kp,nppmx,idimp,nop,mx,my&
     &,mz,mx1,my1,mxyz1,irc)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mz, mx1, my1, mxyz1, irc
         real, dimension(idimp,nop) :: part
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         integer, dimension(mxyz1) :: kpic
         integer, dimension(nppmx,mxyz1) :: kp
         end subroutine
      end interface
!
      interface
         subroutine PPCHECK3LT(ppart,kpic,idimp,nppmx,nx,ny,nz,mx,my,mz,&
     &mx1,my1,mz1,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, my1, mz1
         integer :: irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         integer, dimension(mx1*my1*mz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GPPUSH3LT(ppart,fxyz,kpic,qbm,dt,ek,idimp,nppmx,nx, &
     &ny, nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nyv, nzv
         integer :: mx1, my1, mxyz1, ipbc
         real :: qbm, dt, ek
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         real, dimension(4,nxv,nyv,nzv) :: fxyz
         integer, dimension(mxyz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GPPUSHF3LT(ppart,fxyz,kpic,ncl,ihole,qbm,dt,ek,idimp&
     &,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nyv, nzv
         integer :: mx1, my1, mxyz1, ntmax, irc
         real :: qbm, dt, ek
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         real, dimension(4,nxv,nyv,nzv) :: fxyz
         integer, dimension(mxyz1) :: kpic
         integer, dimension(26,mxyz1) :: ncl
         integer, dimension(2,ntmax+1,mxyz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VGPPUSH3LT(ppart,fxyz,kpic,qbm,dt,ek,idimp,nppmx,nx,&
     &ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nyv, nzv
         integer :: mx1, my1, mxyz1, ipbc
         real :: qbm, dt, ek
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         real, dimension(4,nxv,nyv,nzv) :: fxyz
         integer, dimension(mxyz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGPPUSHF3LT(ppart,fxyz,kpic,ncl,ihole,qbm,dt,ek,    &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nyv, nzv
         integer :: mx1, my1, mxyz1, ntmax, irc
         real :: qbm, dt, ek
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         real, dimension(4,nxv,nyv,nzv) :: fxyz
         integer, dimension(mxyz1) :: kpic
         integer, dimension(26,mxyz1) :: ncl
         integer, dimension(2,ntmax+1,mxyz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GPPOST3LT(ppart,q,kpic,qm,nppmx,idimp,mx,my,mz,nxv, &
     &nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer :: nppmx, idimp, mx, my, mz, nxv, nyv, nzv, mx1, my1
         integer :: mxyz1
         real :: qm
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         real, dimension(nxv,nyv,nzv) :: q
         integer, dimension(mxyz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGPPOST3LT(ppart,q,kpic,qm,nppmx,idimp,mx,my,mz,nxv,&
     &nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer :: nppmx, idimp, mx, my, mz, nxv, nyv, nzv, mx1, my1
         integer :: mxyz1
         real :: qm
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         real, dimension(nxv,nyv,nzv) :: q
         integer, dimension(mxyz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPORDER3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &nx,ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, my1, mz1
         integer :: npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1*mz1) :: ppbuff
         integer, dimension(mx1*my1*mz1) :: kpic
         integer, dimension(26,mx1*my1*mz1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPORDERF3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,&
     &mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, my1, mz1, npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1*mz1) :: ppbuff
         integer, dimension(mx1*my1*mz1) :: kpic
         integer, dimension(26,mx1*my1*mz1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VPPORDER3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,&
     &nx,ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, my1, mz1
         integer :: npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1*mz1) :: ppbuff
         integer, dimension(mx1*my1*mz1) :: kpic
         integer, dimension(26,mx1*my1*mz1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VPPORDERF3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx&
     &,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, my1, mz1, npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1*mz1) :: ppbuff
         integer, dimension(mx1*my1*mz1) :: kpic
         integer, dimension(26,mx1*my1*mz1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine V2PPORDERF3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,    &
     &nppmx,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, my1, mz1, npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         real, dimension(idimp,npbmx,mx1*my1*mz1) :: ppbuff
         integer, dimension(mx1*my1*mz1) :: kpic
         integer, dimension(26,mx1*my1*mz1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine CGUARD3L(fxyz,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer :: nx, ny, nz, nxe, nye, nze
         real, dimension(4,nxe,nye,nze) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine AGUARD3L(q,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer :: nx, ny, nz, nxe, nye, nze
         real, dimension(nxe,nye,nze) :: q
         end subroutine
      end interface
!
      interface
         subroutine VMPOIS33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,&
     &nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer :: isign, nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
         real :: ax, ay, az, affp, we
         real, dimension(2*nxvh,nyv,nzv) :: q
         real, dimension(4,2*nxvh,nyv,nzv) :: fxyz
         complex, dimension(nxhd,nyhd,nzhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
         implicit none
         integer :: indx, indy, indz, nxhyzd, nxyzhd
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RVMX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,&
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RVM3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,&
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RVMXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RVMXZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RVM3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,&
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RVM3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine SET_SZERO3(q,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer :: mx, my, mz, nxv, nyv, nzv, mx1, my1, mxyz1
         real, dimension(nxv,nyv,nzv) :: q
         end subroutine
      end interface
!
      interface
         subroutine SET_VZERO3(cu,mx,my,mz,ndim,nxv,nyv,nzv,mx1,my1,    &
     &mxyz1)
         implicit none
         integer :: mx, my, mz, ndim, nxv, nyv, nzv, mx1, my1, mxyz1
         real, dimension(ndim,nxv,nyv,nzv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine SET_CVZERO3(exyz,nx,ny,nz,ndim,nxvh,nyv,nzv)
         implicit none
         integer :: nx, ny, nz, ndim, nxvh, nyv, nzv
         complex, dimension(ndim,nxvh,nyv,nzv) :: exyz
         end subroutine
      end interface
!
      interface
         function ranorm()
         implicit none
         double precision :: ranorm
         end function
      end interface
!
      end module
