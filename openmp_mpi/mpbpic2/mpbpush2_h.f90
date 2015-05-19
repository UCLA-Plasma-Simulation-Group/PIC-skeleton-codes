!-----------------------------------------------------------------------
! Interface file for mpbpush2.f
      module mpbpush2_h
      implicit none
!
      interface
         subroutine PDICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,  &
     &idps)
         implicit none
         integer :: nyp, noff, nypmx, nypmn, ny, kstrt, nvp, idps
         real, dimension(idps) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PDISTR2H(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,&
     &npx,npy,nx,ny,idimp,npmax,idps,ipbc,ierr)
         implicit none
         integer :: npp, nps, npx, npy, nx, ny, idimp, npmax, idps, ipbc
         integer :: ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax) :: part
         real, dimension(idps) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my&
     &,mx1,mxyp1,irc)
         implicit none
         integer :: nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
         integer :: npp, noff
         real, dimension(idimp,npmax) :: part
         integer, dimension(mxyp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,    &
     &npmax,mx,my,mx1,mxyp1,irc)
         implicit none
         integer :: nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
         integer :: npp, noff
         real, dimension(idimp,npmax) :: part
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         integer, dimension(mxyp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,&
     &mx1,myp1,irc)
         implicit none
         integer :: idimp, nppmx, nx, mx, my, mx1, myp1, irc
         integer :: noff, nyp
         real, dimension(idimp,nppmx,mx1*myp1) :: ppart
         integer, dimension(mx1*myp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,&
     &ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer :: noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
         integer :: mx1, mxyp1, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(3,nxv,nypmx) :: fxy, bxy
         integer, dimension(mxyp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,&
     &qbm,dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,  &
     &irc)
         implicit none
         integer :: noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
         integer :: mx1, mxyp1, ntmax, irc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(3,nxv,nypmx) :: fxy, bxy
         integer, dimension(mxyp1) :: kpic
         integer, dimension(8,mxyp1) :: ncl
         integer, dimension(2,ntmax+1,mxyp1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc&
     &,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer :: noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
         integer :: mx1, mxyp1, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(3,nxv,nypmx) :: fxy, bxy
         integer, dimension(mxyp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp&
     &,qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,    &
     &ntmax,irc)
         implicit none
         integer :: noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
         integer :: mx1, mxyp1, ntmax, irc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(3,nxv,nypmx) :: fxy, bxy
         integer, dimension(mxyp1) :: kpic
         integer, dimension(8,mxyp1) :: ncl
         integer, dimension(2,ntmax+1,mxyp1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,  &
     &nxv, nypmx,mx1,mxyp1)
         implicit none
         integer :: noff, idimp, nppmx, mx, my, nxv, nypmx, mx1, mxyp1
         real :: qm
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(nxv,nypmx) :: q
         integer, dimension(mxyp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGJPPOST2L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,&
     &ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer :: noff, nppmx, idimp, nx, ny, mx, my, nxv, nypmx
         integer :: mx1, mxyp1, ipbc
         real :: qm, dt
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(3,nxv,nypmx) :: cu
         integer, dimension(mxyp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,&
     &nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
         implicit none
         integer :: noff, nyp, nppmx, idimp, nx, ny, mx, my, nxv, nypmx
         integer :: mx1, mxyp1, ntmax, irc
         real :: qm, dt
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(3,nxv,nypmx) :: cu
         integer, dimension(mxyp1) :: kpic
         integer, dimension(8,mxyp1) :: ncl
         integer, dimension(2,ntmax+1,mxyp1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGRJPPOST2L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp&
     &,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer :: noff, nppmx, idimp, nx, ny, mx, my, nxv, nypmx
         integer :: mx1, mxyp1, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(3,nxv,nypmx) :: cu
         integer, dimension(mxyp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGRJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt&
     &,ci,nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
         implicit none
         integer :: noff, nyp, nppmx, idimp, nx, ny, mx, my, nxv, nypmx
         integer :: mx1, mxyp1, ntmax, irc
         real :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(3,nxv,nypmx) :: cu
         integer, dimension(mxyp1) :: kpic
         integer, dimension(8,mxyp1) :: ncl
         integer, dimension(2,ntmax+1,mxyp1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,&  
     &ncll,nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,  &
     &nbmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, mx1, myp1
         integer :: npbmx, ntmax, nbmax, irc
         integer :: noff, nyp
         real, dimension(idimp,nppmx,mx1*myp1) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1) :: ppbuff
         real, dimension(idimp,nbmax) :: sbufl, sbufr
         integer, dimension(mx1*myp1) :: kpic
         integer, dimension(8,mx1*myp1) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1) :: ihole
         integer, dimension(3,mx1) :: ncll, nclr
         end subroutine
      end interface
!
      interface
         subroutine PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll&
     &,nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, myp1, npbmx, ntmax, nbmax, irc
         real, dimension(idimp,nppmx,mx1*myp1) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1) :: ppbuff
         real, dimension(idimp,nbmax) :: sbufl, sbufr
         integer, dimension(8,mx1*myp1) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1) :: ihole
         integer, dimension(3,mx1) :: ncll, nclr
         end subroutine
      end interface
!
      interface
         subroutine PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,&
     &mcll,mclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, myp1, npbmx, ntmax, nbmax, irc
         real, dimension(idimp,nppmx,mx1*myp1) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1) :: ppbuff
         real, dimension(idimp,nbmax) :: rbufl, rbufr
         integer, dimension(mx1*myp1) :: kpic
         integer, dimension(8,mx1*myp1) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1) :: ihole
         integer, dimension(3,mx1) :: mcll, mclr
         end subroutine
      end interface
!
      interface
         subroutine PPCGUARD2XL(fxy,nyp,nx,ndim,nxe,nypmx)
         implicit none
         integer :: nyp, nx, ndim, nxe, nypmx
         real, dimension(ndim,nxe,nypmx) :: fxy
         end subroutine
      end interface
!
      interface
         subroutine PPAGUARD2XL(q,nyp,nx,nxe,nypmx)
         implicit none
         integer :: nyp, nx, nxe, nypmx
         real, dimension(nxe,nypmx) :: q
         end subroutine
      end interface
!
      interface
         subroutine PPACGUARD2XL(cu,nyp,nx,ndim,nxe,nypmx)
         implicit none
         integer :: nyp, nx, ndim, nxe, nypmx
         real, dimension(ndim,nxe,nypmx) :: cu
         end subroutine
      end interface
!
      interface
         subroutine MPPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt, &
     &nyv,kxp,nyhd)
         implicit none
         integer :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real :: ax, ay, affp, we
         complex, dimension(nyv,kxp) :: q
         complex, dimension(3,nyv,kxp) :: fxy
         complex, dimension(nyhd,kxp) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPCUPERP2(cu,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp) :: cu
         end subroutine
      end interface
!
      interface
         subroutine MIPPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,  &
     &nyhd)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, nyhd
         real :: ci, wm
         complex, dimension(3,nyv,kxp) :: cu
         complex, dimension(3,nyv,kxp) :: bxy
         complex, dimension(nyhd,kxp) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,   &
     &kstrt,nyv,kxp,nyhd)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, nyhd
         real :: affp, ci, dt, wf, wm
         complex, dimension(3,nyv,kxp) :: exy, bxy, cu
         complex, dimension(nyhd,kxp) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,  &
     &nyhd)
         implicit none
         integer :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         complex, dimension(3,nyv,kxp) :: fxy, exy
         complex, dimension(nyhd,kxp) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine WPFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer :: indx, indy, nxhyd, nxyhd
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WPPFFT2RM(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
         integer :: kxp, kyp, kypd, nxhyd, nxyhd
         real :: ttp
         real, dimension(2*nxvh,kypd) :: f
         complex, dimension(nyv,kxp) :: g
         complex, dimension(kxp,kyp) :: bs, br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WPPFFT2RM2(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx&
     &,indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
         integer :: kxp, kyp, kypd, nxhyd, nxyhd
         real :: ttp
         real, dimension(2,2*nxvh,kypd) :: f
         complex, dimension(2,nyv,kxp) :: g
         complex, dimension(2,kxp,kyp) :: bs, br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RMXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,  &
     &kypp,nxvh,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, nxvh, kypi, kypp, kypd
         integer :: nxhyd, nxyhd
         real, dimension(2*nxvh,kypd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RMXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,  &
     &kxpp,nyv,kxp,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxp
         integer :: nxhyd, nxyhd
         complex, dimension(nyv,kxp) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RM3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, nxvh, kypi, kypp, kypd
         integer :: nxhyd, nxyhd
         real, dimension(3,2*nxvh,kypd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RM3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxp,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxp
         integer :: nxhyd, nxyhd
         complex, dimension(3,nyv,kxp) :: g
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