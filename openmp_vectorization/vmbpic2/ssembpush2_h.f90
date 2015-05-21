!-----------------------------------------------------------------------
! Interface file for ssembpush2.c
      module ssembpush2_h
      implicit none
!
      interface
         subroutine csse2gbppush23lt(ppart,fxy,bxy,kpic,qbm,dt,dtc,ek,  &
     &idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
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
         subroutine csse2gbppushf23lt(ppart,fxy,bxy,kpic,ncl,ihole,qbm, &
     &dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
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
         subroutine csse2grbppush23lt(ppart,fxy,bxy,kpic,qbm,dt,dtc,ci, &
     &ek,idimp, nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
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
         subroutine csse2grbppushf23lt(ppart,fxy,bxy,kpic,ncl,ihole,qbm,&
     &dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
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
         subroutine csse2gppost2lt(ppart,q,kpic,qm,nppmx,idimp,mx,my,nxv&
     &,nyv,mx1,mxy1)
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
         subroutine csse2gjppost2lt(ppart,cu,kpic,qm,dt,nppmx,idimp,nx, &
     &ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
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
         subroutine csse2gjppostf2lt(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx&
     &,idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
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
         subroutine csse2grjppost2lt(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,&
     &nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
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
         subroutine csse2grjppostf2lt(ppart,cu,kpic,ncl,ihole,qm,dt,ci, &
     &nppmx,idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
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
         subroutine csse2pporder2lt(ppart,ppbuff,kpic,ncl,ihole,idimp,  &
     &nppmx,nx,ny,mx,my,mx1,my1,npbmx,ntmax,irc)
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
         subroutine csse2pporderf2lt(ppart,ppbuff,kpic,ncl,ihole,idimp, &
     &nppmx,mx1,my1,npbmx,ntmax,irc)
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
         subroutine csse2bguard2l(bxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(4,nxe,nye) :: bxy
         end subroutine
      end interface
!
      interface
         subroutine csse2acguard2l(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(4,nxe,nye) :: cu
         end subroutine
      end interface
!
      interface
         subroutine csse2aguard2l(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
!
      interface
         subroutine csse2mpois23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,   &
     &nxvh,nyv,nxhd,nyhd)
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
         subroutine csse2mcuperp2(cu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
         real, dimension(4,2*nxvh,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine csse2mibpois23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,&
     &nyhd)
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
         subroutine csse2mmaxwel2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,&
     &nyv,nxhd,nyhd)
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
         subroutine csse2memfield2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd&
     &,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, dimension(4,2*nxvh,nyv) :: fxy
         complex, dimension(4,nxvh,nyv) :: exy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine csse2wfft2rmx(f,isign,mixup,sct,indx,indy,nxhd,nyd, &
     &nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         real, dimension(2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine csse2wfft2rm3(f,isign,mixup,sct,indx,indy,nxhd,nyd, &
     &nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         real, dimension(4,2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      end module
