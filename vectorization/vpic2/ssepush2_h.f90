!-----------------------------------------------------------------------
! Interface file for ssepush2.c
      module ssepush2_h
      implicit none
!
      interface
         subroutine csse2xiscan2(isdata,nths)
         implicit none
         integer :: nths
         integer, dimension(nths) :: isdata
         end subroutine
      end interface
!
      interface
         subroutine csse2gpush2lt(part,fxy,qbm,dt,ek,idimp,nop,npe,nx,ny&
     &,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, npe, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(npe,idimp) :: part
         real, dimension(2,nxv,nyv) :: fxy
         end subroutine
      end interface
!
      interface
         subroutine csse2gpost2lt(part,q,qm,nop,npe,idimp,nxv,nyv)
         implicit none
         integer :: nop, npe, idimp, nxv, nyv
         real :: qm
         real, dimension(npe,idimp) :: part
         real, dimension(nxv,nyv) :: q
         end subroutine
      end interface
!
      interface
         subroutine csse2dsortp2ylt(parta,partb,npic,idimp,nop,npe,ny1)
         implicit none
         integer :: idimp, nop, npe, ny1
         real, dimension(npe,idimp) :: parta, partb
         integer, dimension(ny1) :: npic
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
         subroutine csse2pois22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh&
     &,nyv,nxhd,nyhd)
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
         subroutine csse2wfft2rx(f,isign,mixup,sct,indx,indy,nxhd,nyd,  &
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
         subroutine csse2wfft2r2(f,isign,mixup,sct,indx,indy,nxhd,nyd,  &
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
