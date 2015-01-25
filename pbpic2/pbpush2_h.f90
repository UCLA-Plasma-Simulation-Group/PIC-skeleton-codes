!-----------------------------------------------------------------------
! Interface file for pbpush2.f
      module pbpush2_h
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
         subroutine PPGBPUSH23L(part,fxy,bxy,edges,npp,noff,ihole,qbm,dt&
     &,dtc,ek,nx,ny,idimp,npmax,nxv,nypmx,idps,ntmax,ipbc)
         implicit none
         integer :: npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv
         integer :: nypmx, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax) :: part
         real, dimension(3,nxv,nypmx) :: fxy, bxy
         real, dimension(idps) :: edges
         integer, dimension(ntmax+1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGRBPUSH23L(part,fxy,bxy,edges,npp,noff,ihole,qbm, &
     &dt,dtc,ci,ek,nx,ny,idimp,npmax,nxv,nypmx,idps,ntmax,ipbc)
         implicit none
         integer :: npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv
         integer :: nypmx, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax) :: part
         real, dimension(3,nxv,nypmx) :: fxy, bxy
         real, dimension(idps) :: edges
         integer, dimension(ntmax+1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGPOST2L(part,q,npp,noff,qm,idimp,npmax,nxv,nypmx)
         implicit none
         integer :: npp, noff, idimp, npmax, nxv, nypmx
         real :: qm
         real, dimension(idimp,npmax) :: part
         real, dimension(nxv,nypmx) :: q
         end subroutine
      end interface
!
      interface
         subroutine PPGJPOST2L(part,cu,edges,npp,noff,ihole,qm,dt,nx,ny,&
     &idimp,npmax,nxv,nypmx,idps,ntmax,ipbc)
         implicit none
         integer :: npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv
         integer :: nypmx, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax) :: part
         real, dimension(3,nxv,nypmx) :: cu
         real, dimension(idps) :: edges
         integer, dimension(ntmax+1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGRJPOST2L(part,cu,edges,npp,noff,ihole,qm,dt,ci,nx&
     &,ny,idimp,npmax,nxv,nypmx,idps,ntmax,ipbc)
         implicit none
         integer :: npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv
         integer :: nypmx, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax) :: part
         real, dimension(3,nxv,nypmx) :: cu
         real, dimension(idps) :: edges
         integer, dimension(ntmax+1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPDSORTP2YL(parta,partb,npic,npp,noff,nyp,idimp,    &
     &npmax,nypm1)
         implicit none
         integer :: npp, noff, nyp, idimp, npmax, nypm1
         real, dimension(idimp,npmax) :: parta, partb
         integer, dimension(nypm1) :: npic
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
         subroutine PPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,  &
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
         subroutine PPCUPERP2(cu,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp) :: cu
         end subroutine
      end interface
!
      interface
         subroutine IPPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,   &
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
         subroutine PPMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,    &
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
         subroutine PPEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,   &
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
         subroutine WPPFFT2R(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx, &
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
         subroutine WPPFFT2R3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
         integer :: kxp, kyp, kypd, nxhyd, nxyhd
         real :: ttp
         real, dimension(3,2*nxvh,kypd) :: f
         complex, dimension(3,nyv,kxp) :: g
         complex, dimension(3,kxp,kyp) :: bs, br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,   &
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
         subroutine PPFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,   &
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
         subroutine PPFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,  &
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
         subroutine PPFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,  &
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