!-----------------------------------------------------------------------
! Interface file for ppush2.f
      module ppush2_h
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
         subroutine PDISTR2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy, &
     &nx,ny,idimp,npmax,idps,ipbc,ierr)
         implicit none
         integer :: npp, nps, npx, npy, nx, ny, idimp, npmax, idps, ipbc
         integer :: ierr
         real :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax) :: part
         real, dimension(idps) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PPGPUSH2L(part,fxy,edges,npp,noff,ihole,qbm,dt,ek,nx&
     &,ny,idimp,npmax,nxv,nypmx,idps,ntmax,ipbc)
         implicit none
         integer :: npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv
         integer :: nypmx, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax) :: part
         real, dimension(2,nxv,nypmx) :: fxy
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
         subroutine PPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,  &
     &nyv,kxp,nyhd)
         implicit none
         integer :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real :: ax, ay, affp, we
         complex, dimension(nyv,kxp) :: q
         complex, dimension(2,nyv,kxp) :: fxy
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
         subroutine WPPFFT2R2(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
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
         subroutine PPFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,  &
     &kypp,nxvh,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, nxvh, kypi, kypp, kypd
         integer :: nxhyd, nxyhd
         real, dimension(2,2*nxvh,kypd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,  &
     &kxpp,nyv,kxp,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxp
         integer :: nxhyd, nxyhd
         complex, dimension(2,nyv,kxp) :: g
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