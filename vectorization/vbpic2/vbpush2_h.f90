!-----------------------------------------------------------------------
! Interface file for vbpush2.f
      module vbpush2_h
      implicit none
!
      interface
         subroutine DISTR2HT(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp,&
     &npe,nx,ny,ipbc)
         implicit none
         integer :: npx, npy, idimp, npe, nx, ny, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(npe,idimp) :: part
         end subroutine
      end interface
!
      interface
         subroutine GBPUSH23LT(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,npe,&
     &nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, npe, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine GRBPUSH23LT(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &npe,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, npe, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine VGBPUSH23LT(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,npe&
     &,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, npe, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine VGRBPUSH23LT(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop&
     &,npe,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, npe, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine GPOST2LT(part,q,qm,nop,npe,idimp,nxv,nyv)
         implicit none
         integer :: nop, npe, idimp, nxv, nyv
         real :: qm
         real, dimension(npe,idimp) :: part
         real, dimension(nxv,nyv) :: q
         end subroutine
      end interface
!
      interface
         subroutine VGPOST2LT(part,q,qm,nop,npe,idimp,nxv,nyv)
         implicit none
         integer :: nop, npe, idimp, nxv, nyv
         real :: qm
         real, dimension(npe,idimp) :: part
         real, dimension(nxv,nyv) :: q
         end subroutine
      end interface
!
      interface
         subroutine GJPOST2LT(part,cu,qm,dt,nop,npe,idimp,nx,ny,nxv,nyv,&
     &ipbc)
         implicit none
         integer :: nop, npe, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine GRJPOST2LT(part,cu,qm,dt,ci,nop,npe,idimp,nx,ny,nxv,&
     &nyv,ipbc)
         implicit none
         integer :: nop, npe, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt, ci
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine VGJPOST2LT(part,cu,qm,dt,nop,npe,idimp,nx,ny,nxv,nyv&
     &,ipbc)
         implicit none
         integer :: nop, npe, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine VGRJPOST2LT(part,cu,qm,dt,ci,nop,npe,idimp,nx,ny,nxv&
     &,nyv,ipbc)
         implicit none
         integer :: nop, npe, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt, ci
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine DSORTP2YLT(parta,partb,npic,idimp,nop,npe,ny1)
         implicit none
         integer :: idimp, nop, npe, ny1
         real, dimension(npe,idimp) :: parta, partb
         integer, dimension(ny1) :: npic
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
         subroutine VPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv&
     &,nxhd,nyhd)
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
         subroutine CUPERP2(cu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
         real, dimension(4,2*nxvh,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine VIBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
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
         subroutine VMAXWEL2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv, &
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
         subroutine VEMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd&
     &)
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
         subroutine WFFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd, &
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
         subroutine WFFT2RVX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
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
         subroutine WFFT2RV3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         real, dimension(4,2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RVXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  &
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
         subroutine FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd&
     &,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
         real, dimension(2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RV3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  &
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
         subroutine FFT2RV3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,  &
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
      interface
         function randum()
         implicit none
         double precision :: randum
         end function
      end interface
!
      end module
