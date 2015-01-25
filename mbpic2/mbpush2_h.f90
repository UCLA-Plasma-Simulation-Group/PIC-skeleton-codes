!-----------------------------------------------------------------------
! Interface file for mbpush2.f
      module mbpush2_h
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
         subroutine GBPPUSH23L(ppart,fxy,bxy,kpic,qbm,dt,dtc,ek,idimp,  &
     &nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nppmx,mxy1) :: ppart
         real, dimension(3,nxv,nyv) :: fxy, bxy
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,dtc,&
     &ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nppmx,mxy1) :: ppart
         real, dimension(3,nxv,nyv) :: fxy, bxy
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GRBPPUSH23L(ppart,fxy,bxy,kpic,qbm,dt,dtc,ci,ek,    &
     &idimp, nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nppmx,mxy1) :: ppart
         real, dimension(3,nxv,nyv) :: fxy, bxy
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,dtc&
     &,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nppmx,mxy1) :: ppart
         real, dimension(3,nxv,nyv) :: fxy, bxy
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
         subroutine GJPPOST2L(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny,mx, &
     &my, nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qm, dt
         real, dimension(idimp,nppmx,mxy1) :: ppart
         real, dimension(3,nxv,nyv) :: cu
         integer, dimension(mxy1) :: kpic

         end subroutine
      end interface
!
      interface
         subroutine GJPPOSTF2L(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,idimp&
     &,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qm, dt
         real, dimension(idimp,nppmx,mxy1) :: ppart
         real, dimension(3,nxv,nyv) :: cu
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GRJPPOST2L(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx,ny,&
     &mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nppmx,mxy1) :: ppart
         real, dimension(3,nxv,nyv) :: cu
         integer, dimension(mxy1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRJPPOSTF2L(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx, &
     &idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax, irc
         real :: qm, dt, ci
         real, dimension(idimp,nppmx,mxy1) :: ppart
         real, dimension(3,nxv,nyv) :: cu
         integer, dimension(mxy1) :: kpic
         integer, dimension(8,mxy1) :: ncl
         integer, dimension(2,ntmax+1,mxy1) :: ihole
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
         subroutine BGUARD2L(bxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: bxy
         end subroutine
      end interface
!
      interface
         subroutine ACGUARD2L(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cu
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
         subroutine MPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv&
     &,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ax, ay, affp, we
         real, dimension(2*nxvh,nyv) :: q
         real, dimension(3,2*nxvh,nyv) :: fxy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MCUPERP2(cu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
         real, dimension(3,2*nxvh,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine MIBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ci, wm
         real, dimension(3,2*nxvh,nyv) :: cu
         complex, dimension(3,nxvh,nyv) :: bxy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MMAXWEL2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv, &
     &nxhd,nyhd)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ci, dt, wf, wm
         complex, dimension(3,nxvh,nyv) :: exy, bxy
         real, dimension(3,2*nxvh,nyv) :: cu
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MEMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd&
     &)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, dimension(3,2*nxvh,nyv) :: fxy
         complex, dimension(3,nxvh,nyv) :: exy
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
         subroutine WFFT2RM3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd) :: f
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
         subroutine FFT2RM3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RM3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd) :: f
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
