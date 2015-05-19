!-----------------------------------------------------------------------
! Interface file for field3.f
      module field3_h
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
      end module
