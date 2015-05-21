!-----------------------------------------------------------------------
! Interface file for ssempush2.c
      module ssempush2_h
      implicit none
!
      interface
         subroutine csse2gppush2lt(ppart,fxy,kpic,qbm,dt,ek,idimp,nppmx,&
     &nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qbm, dt, ek
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(2,nxv*nyv) :: fxy
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine csse2gppushf2lt(ppart,fxy,kpic,ncl,ihole,qbm,dt,ek, &
     &idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qbm, dt, ek
         real, dimension(nppmx,idimp,mxy1) :: ppart
         real, dimension(2,nxv*nyv) :: fxy
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
         subroutine csse2cguard2l(fxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: fxy
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
         subroutine csse2mpois22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,   &
     &nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ax, ay, affp, we
         real, dimension(2*nxvh,nyv) :: q
         real, dimension(2,2*nxvh,nyv) :: fxy
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
         subroutine csse2wfft2rm2(f,isign,mixup,sct,indx,indy,nxhd,nyd, &
     &nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         real, dimension(2,2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      end module
