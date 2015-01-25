!-----------------------------------------------------------------------
! Interface file for mpush2.f
      module mpush2_h
      implicit none
!
      interface
         subroutine DISTR2(part,vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx,ny,&
     &ipbc)
         implicit none
         integer :: npx, npy, idimp, nop, nx, ny, ipbc
         real :: vtx, vty, vdx, vdy
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
         subroutine PPMOVIN2L(part,ppart,kpic,nppmx,idimp,nop,mx,my,mx1,&
     &mxy1,irc)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mx1, mxy1, irc
         real, dimension(idimp,nop) :: part
         real, dimension(idimp,nppmx,mxy1) :: ppart
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPCHECK2L(ppart,kpic,idimp,nppmx,nx,ny,mx,my,mx1,my1&
     &,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, mx1, my1, irc
         real, dimension(idimp,nppmx,mx1*my1) :: ppart
         integer, dimension(mx1*my1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GPPUSH2L(ppart,fxy,kpic,qbm,dt,ek,idimp,nppmx,nx,ny,&
     &mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nppmx,mxy1) :: ppart
         real, dimension(2,nxv,nyv) :: fxy
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GPPUSHF2L(ppart,fxy,kpic,ncl,ihole,qbm,dt,ek,idimp, &
     &nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qbm, dt, ek
         real, dimension(idimp,nppmx,mxy1) :: ppart
         real, dimension(2,nxv,nyv) :: fxy
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GPPOST2L(ppart,q,kpic,qm,nppmx,idimp,mx,my,nxv,nyv, &
     &mx1,mxy1)
         implicit none
         integer :: nppmx, idimp, mx, my, nxv, nyv, mx1, mxy1
         real :: qm
         real, dimension(idimp,nppmx,mxy1) :: ppart
         real, dimension(nxv,nyv) :: q
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPORDER2L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx&
     &,ny,mx,my,mx1,my1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, mx1, my1, npbmx, ntmax
         integer :: irc
         real, dimension(idimp,nppmx,mx1*my1) :: ppart
         real, dimension(idimp,npbmx,mx1*my1) :: ppbuff
         integer, dimension(mx1*my1) :: kpic
         integer, dimension(8,mx1*my1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPORDERF2L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &mx1,my1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, my1, npbmx, ntmax, irc
         real, dimension(idimp,nppmx,mx1*my1) :: ppart
         real, dimension(idimp,npbmx,mx1*my1) :: ppbuff
         integer, dimension(mx1*my1) :: kpic
         integer, dimension(8,mx1*my1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine CGUARD2L(fxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: fxy
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
         subroutine MPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv&
     &,nxhd,nyhd)
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
         subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer :: indx, indy, nxhyd, nxyhd
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RMX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         real, dimension(2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RM2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         real, dimension(2,2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  &
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
         subroutine FFT2RM2X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
         real, dimension(2,2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RM2Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
         real, dimension(2,2*nxhd,nyd) :: f
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