!-----------------------------------------------------------------------
! Interface file for bfield3.f
      module bfield3_h
      implicit none
!
      interface
         subroutine POTP3(q,pot,ffc,we,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd, &
     &nzhd)
         implicit none
         integer :: nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
         real :: we
         real, dimension(2*nxvh,nyv,nzv) :: q
         complex, dimension(nxvh,nyv,nzv) :: pot
         complex, dimension(nxhd,nyhd,nzhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine DIVF3(f,df,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer :: nx, ny, nz, nxvh, nyv, nzv
         complex, dimension(3,nxvh,nyv,nzv) :: f
         complex, dimension(nxvh,nyv,nzv) :: df
         end subroutine
      end interface
!
      interface
         subroutine GRADF3(df,f,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer :: nx, ny, nz, nxvh, nyv, nzv
         complex, dimension(nxvh,nyv,nzv) :: df
         complex, dimension(3,nxvh,nyv,nzv) :: f
         end subroutine
      end interface
!
      interface
         subroutine CURLF3(f,g,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer :: nx, ny, nz, nxvh, nyv, nzv
         complex, dimension(3,nxvh,nyv,nzv) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine AVPOT33(bxyz,axyz,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer :: nx, ny, nz, nxvh, nyv, nzv
         complex, dimension(3,nxvh,nyv,nzv) :: bxyz, axyz
         end subroutine
      end interface
!
      interface
         subroutine AVRPOT33(axyz,bxyz,ffc,ci,nx,ny,nz,nxvh,nyv,nzv,nxhd&
     &,nyhd,nzhd)
         implicit none
         integer :: nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
         real :: ci
         complex, dimension(3,nxvh,nyv,nzv) :: axyz, bxyz
         complex, dimension(nxhd,nyhd,nzhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine SMOOTH3(q,qs,ffc,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,   &
     &nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, dimension(2*nxvh,nyv,nzv), intent(in) :: q
         complex, dimension(nxvh,nyv,nzv), intent(inout) :: qs
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine SMOOTH33(cu,cus,ffc,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,&
     &nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, dimension(3,2*nxvh,nyv,nzv), intent(in) :: cu
         complex, dimension(3,nxvh,nyv,nzv), intent(inout) :: cus
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine RDMODES3(pot,pott,nx,ny,nz,modesx,modesy,modesz,nxvh&
     &,nyv,nzv,modesxd,modesyd,modeszd)
         implicit none
         integer :: nx, ny, nz, modesx, modesy, modesz, nxvh, nyv, nzv
         integer :: modesxd, modesyd, modeszd
         complex, dimension(nxvh,nyv,nzv) :: pot
         complex, dimension(modesxd,modesyd,modeszd) :: pott
         end subroutine
      end interface
!
      interface
         subroutine WRMODES3(pot,pott,nx,ny,nz,modesx,modesy,modesz,nxvh&
     &,nyv,nzv,modesxd,modesyd,modeszd)
         implicit none
         integer :: nx, ny, nz, modesx, modesy, modesz, nxvh, nyv, nzv
         integer :: modesxd, modesyd, modeszd
         complex, dimension(nxvh,nyv,nzv) :: pot
         complex, dimension(modesxd,modesyd,modeszd) :: pott
         end subroutine
      end interface
!
      interface
         subroutine RDVMODES3(vpot,vpott,nx,ny,nz,modesx,modesy,modesz, &
     &ndim,nxvh,nyv,nzv,modesxd,modesyd,modeszd)
         implicit none
         integer :: nx, ny, nz, modesx, modesy, modesz, ndim
         integer :: nxvh, nyv, nzv, modesxd, modesyd, modeszd
         complex, dimension(ndim,nxvh,nyv,nzv) :: vpot
         complex, dimension(ndim,modesxd,modesyd,modeszd) :: vpott
         end subroutine
      end interface
!
      interface
         subroutine WRVMODES3(vpot,vpott,nx,ny,nz,modesx,modesy,modesz, &
     &ndim,nxvh,nyv,nzv,modesxd,modesyd,modeszd)
         implicit none
         integer :: nx, ny, nz, modesx, modesy, modesz, ndim
         integer :: nxvh, nyv, nzv, modesxd, modesyd, modeszd
         complex, dimension(ndim,nxvh,nyv,nzv) :: vpot
         complex, dimension(ndim,modesxd,modesyd,modeszd) :: vpott
         end subroutine
      end interface
!
      end module
