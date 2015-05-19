!-----------------------------------------------------------------------
! Interface file for ssebpush2.c
      module ssebpush2_h
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
         subroutine csse2gbpush23lt(part,fxy,bxy,qbm,dt,dtc,ek,idimp,   &
     &nop,npe,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, npe, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine csse2grbpush23lt(part,fxy,bxy,qbm,dt,dtc,ci,ek,     &
     &idimp,nop,npe,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, npe, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: fxy, bxy
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
         subroutine csse2gjpost2lt(part,cu,qm,dt,nop,npe,idimp,nx,ny,nxv&
    &,nyv,ipbc)
         implicit none
         integer :: nop, npe, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine csse2grjpost2lt(part,cu,qm,dt,ci,nop,npe,idimp,nx,ny&
     &,nxv,nyv,ipbc)
         implicit none
         integer :: nop, npe, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt, ci
         real, dimension(npe,idimp) :: part
         real, dimension(4,nxv,nyv) :: cu
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
         subroutine csse2bguard2l(bxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(4,nxe,nye) :: bxy
         end subroutine
      end interface
!
      interface
         subroutine csse2acguard2l(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(4,nxe,nye) :: cu
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
         subroutine csse2pois23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh&
     &,nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ax, ay, affp ,we
         real, dimension(2*nxvh,nyv) :: q
         real, dimension(4,2*nxvh,nyv) :: fxy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine csse2cuperp2(cu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
         real, dimension(4,2*nxvh,nyv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine csse2ibpois23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd, &
     &nyhd)
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
         subroutine csse2maxwel2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh, &
     &nyv,nxhd,nyhd)
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
         subroutine csse2emfield2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,&
     &nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, dimension(4,2*nxvh,nyv) :: fxy
         complex, dimension(4,nxvh,nyv) :: exy
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
         subroutine csse2wfft2r3(f,isign,mixup,sct,indx,indy,nxhd,nyd,  &
     &nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         real, dimension(4,2*nxhd,nyd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      end module
