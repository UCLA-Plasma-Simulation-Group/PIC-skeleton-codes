!-----------------------------------------------------------------------
! Interface file for push2.f
      module push2_h
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
         subroutine GPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv, &
     &ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
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
         subroutine DSORTP2YL(parta,partb,npic,idimp,nop,ny1)
         implicit none
         integer :: idimp, nop, ny1
         real, dimension(idimp,nop) :: parta, partb
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
         subroutine POIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,&
     &nxhd,nyhd)
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
         subroutine WFFT2R2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd, &
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
         subroutine FFT2R2X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd&
     &,nxhyd,nxyhd)
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
