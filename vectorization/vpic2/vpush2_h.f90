!-----------------------------------------------------------------------
! Interface file for vpush2.f
      module vpush2_h
      implicit none
!
      interface
         subroutine DISTR2T(part,vtx,vty,vdx,vdy,npx,npy,idimp,npe,nx,ny&
     &,ipbc)
         implicit none
         integer :: npx, npy, idimp, npe, nx, ny, ipbc
         real :: vtx, vty, vdx, vdy
         real, dimension(npe,idimp) :: part
         end subroutine
      end interface
!
      interface
         subroutine GPUSH2LT(part,fxy,qbm,dt,ek,idimp,nop,npe,nx,ny,nxv,&
     &nyv,ipbc)
         implicit none
         integer :: idimp, nop, npe, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(npe,idimp) :: part
         real, dimension(2,nxv,nyv) :: fxy
         end subroutine
      end interface
!
      interface
         subroutine VGPUSH2LT(part,fxy,qbm,dt,ek,idimp,nop,npe,nx,ny,nxv&
     &,nyv,ipbc)
         implicit none
         integer :: idimp, nop, npe, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(npe,idimp) :: part
         real, dimension(2,nxv,nyv) :: fxy
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
         subroutine DSORTP2YLT(parta,partb,npic,idimp,nop,npe,ny1)
         implicit none
         integer :: idimp, nop, npe, ny1
         real, dimension(npe,idimp) :: parta, partb
         integer, dimension(ny1) :: npic
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
         subroutine VPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv&
     &,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ax, ay, affp ,we
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
         subroutine WFFT2RV2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
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
         subroutine FFT2RV2X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  &
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
         subroutine FFT2R2Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd&
     &,nxhyd,nxyhd)
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
      interface
         function randum()
         implicit none
         double precision :: randum
         end function
      end interface
!
      end module
