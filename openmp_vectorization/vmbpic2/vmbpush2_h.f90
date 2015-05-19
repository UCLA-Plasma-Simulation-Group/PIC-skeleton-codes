!-----------------------------------------------------------------------
! Interface file for vmbpush2.f
      module vmbpush2_h
      implicit none
!
      interface
         subroutine DISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp, &
     &nop,nx,ny,ipbc)
         implicit none
         integer :: npx, npy, idimp, nop, nx, ny, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
!
      interface
         subroutine DBLKP2L(part,kpic,nppmx,idimp,nop,mx,my,mx1,mxy1,irc&
     &)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mx1, mxy1, irc
         real, dimension(idimp,nop) :: part
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN2LT(part,ppart,kpic,nppmx,idimp,nop,mx,my,mx1&
     &,mxy1,irc)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mx1, mxy1, irc
         real, dimension(idimp,nop) :: part
         real, dimension(nppmx,idimp,mxy1) :: ppart
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN2LTP(part,ppart,kpic,kp,nppmx,idimp,nop,mx,my&
     &,mx1,mxy1,irc)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mx1, mxy1, irc
         real, dimension(idimp,nop) :: part
         real, dimension(nppmx,idimp,mxy1) :: ppart
         integer, dimension(mxy1) :: kpic
         integer, dimension(nppmx,mxy1) :: kp
         end subroutine
      end interface
!
      interface
         subroutine PPCHECK2LT(ppart,kpic,idimp,nppmx,nx,ny,mx,my,mx1,  &
     &my1,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, mx1, my1, irc
         real, dimension(nppmx,idimp,mx1*my1) :: ppart
         integer, dimension(mx1*my1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSH23LT(ppart,fxy,bxy,kpic,qbm,dt,dtc,ek,idimp, &
     &nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: fxy, bxy
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSHF23LT(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,dtc&
     &,ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qbm, dt, dtc, ek
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: fxy, bxy
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GRBPPUSH23LT(ppart,fxy,bxy,kpic,qbm,dt,dtc,ci,ek,   &
     &idimp, nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: fxy, bxy
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRBPPUSHF23LT(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,  &
     &dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: fxy, bxy
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VGBPPUSH23LT(ppart,fxy,bxy,kpic,qbm,dt,dtc,ek,idimp,&
     &nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: fxy, bxy
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGBPPUSHF23LT(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,  &
     &dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qbm, dt, dtc, ek
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: fxy, bxy
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VGRBPPUSH23LT(ppart,fxy,bxy,kpic,qbm,dt,dtc,ci,ek,  &
     &idimp, nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: fxy, bxy
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGRBPPUSHF23LT(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt, &
     &dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: fxy, bxy
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GPPOST2LT(ppart,q,kpic,qm,nppmx,idimp,mx,my,nxv,nyv,&
     &mx1,mxy1)
         implicit none
         integer :: nppmx, idimp, mx, my, nxv, nyv, mx1, mxy1
         real :: qm
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(nxv,nyv) :: q
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGPPOST2LT(ppart,q,kpic,qm,nppmx,idimp,mx,my,nxv,nyv&
     &,mx1,mxy1)
         implicit none
         integer :: nppmx, idimp, mx, my, nxv, nyv, mx1, mxy1
         real :: qm
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(nxv,nyv) :: q
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GJPPOST2LT(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny,mx,&
     &my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qm, dt
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: cu
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GJPPOSTF2LT(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,    &
     &idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qm, dt
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: cu
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GRJPPOST2LT(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx,ny&
     &,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qm, dt, ci
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: cu
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRJPPOSTF2LT(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx,&
     &idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qm, dt, ci
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: cu
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VGJPPOST2LT(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny,mx&
     &,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qm, dt
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: cu
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGJPPOSTF2LT(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,   &
     &idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qm, dt
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: cu
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VGRJPPOST2LT(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx, &
     &ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qm, dt, ci
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: cu
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGRJPPOSTF2LT(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx&
     &,idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qm, dt, ci
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(4,nxv*nyv) :: cu
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPORDER2LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &nx,ny,mx,my,mx1,my1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, mx1, my1, npbmx, ntmax
         integer :: irc
         real, dimension(nppmx,idimp,mx1*my1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1) :: ppbuff
         integer, dimension(mx1*my1) :: kpic
         integer, dimension(8,mx1*my1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPORDERF2LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,&
     &mx1,my1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, my1, npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1) :: ppbuff
         integer, dimension(mx1*my1) :: kpic
         integer, dimension(8,mx1*my1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VPPORDER2LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,&
     &nx,ny,mx,my,mx1,my1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, mx1, my1, npbmx, ntmax
         integer :: irc
         real, dimension(nppmx,idimp,mx1*my1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1) :: ppbuff
         integer, dimension(mx1*my1) :: kpic
         integer, dimension(8,mx1*my1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VPPORDERF2LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx&
     &,mx1,my1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, my1, npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1) :: ppbuff
         integer, dimension(mx1*my1) :: kpic
         integer, dimension(8,mx1*my1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine BGUARD2L(bxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(4,nxe,nye) :: bxy
         end subroutine
      end interface
!
      interface
         subroutine ACGUARD2L(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(4,nxe,nye) :: cu
         end subroutine
      end interface
!
      interface
         subroutine AGUARD2L(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
!
      interface
         subroutine VMPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,  &
     &nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ax, ay, affp, we
         real, dimension(2*nxvh,nyv) :: q
         real, dimension(4,2*nxvh,nyv) :: fxy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MCUPERP2(cu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
         real, dimension(4,2*nxvh,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine VMIBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd&
     &)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ci, wm
         real, dimension(4,2*nxvh,nyv) :: cu
         complex, dimension(4,nxvh,nyv) :: bxy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine VMMAXWEL2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,&
     &nxhd,nyhd)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ci, dt, wf, wm
         complex, dimension(4,nxvh,nyv) :: exy, bxy
         real, dimension(4,2*nxvh,nyv) :: cu
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine VMEMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,   &
     &nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, dimension(4,2*nxvh,nyv) :: fxy
         complex, dimension(4,nxvh,nyv) :: exy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer :: indx, indy, nxhyd, nxyhd
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RVMX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd&
     &,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         real, dimension(2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RVM3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd&
     &,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         real, dimension(4,2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RVMXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd, &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
         real, dimension(2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
         real, dimension(2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RVM3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd, &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
         real, dimension(4,2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RVM3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd, &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
         real, dimension(4,2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
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
