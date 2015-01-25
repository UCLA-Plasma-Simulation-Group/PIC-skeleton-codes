!-----------------------------------------------------------------------
! Interface file for bpush2.f
      module bpush2_h
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
         subroutine GBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine GRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop, &
     &nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine GPOST2L(part,q,qm,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,nyv) :: q
         end subroutine
      end interface
!
      interface
         subroutine GJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine GRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv, &
     &ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine DSORTP2YL(parta,partb,npic,idimp,nop,ny1)
         implicit none
         integer :: idimp, nop, ny1
         real, dimension(idimp,nop) :: parta, partb
         integer, dimension(ny1) :: npic
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
         subroutine POIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,&
     &nxhd,nyhd)
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
         subroutine CUPERP2(cu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
         real, dimension(3,2*nxvh,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine IBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
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
         subroutine MAXWEL2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,  &
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
         subroutine EMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
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
         subroutine WFFT2R3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd, &
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
         subroutine FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd&
     &,nxhyd,nxyhd)
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
         subroutine FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd&
     &,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd&
     &,nxhyd,nxyhd)
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
      interface
         function randum()
         implicit none
         double precision :: randum
         end function
      end interface
!
      end module
