!-----------------------------------------------------------------------
! Interface file for field1.f
      module field1_h
      implicit none
!
      interface
         subroutine POTP1(q,pot,ffc,we,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(inout) :: we
         real, dimension(2*nxvh), intent(in) :: q
         complex, dimension(nxvh), intent(inout) :: pot
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine DIVF1(f,df,nx,ndim,nxvh)
         implicit none
         integer, intent(in) :: nx, ndim, nxvh
         complex, dimension(ndim,nxvh), intent(in) :: f
         complex, dimension(nxvh), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine GRADF1(df,f,nx,ndim,nxvh)
         implicit none
         integer, intent(in) :: nx, ndim, nxvh
         complex, dimension(nxvh), intent(in) :: df
         complex, dimension(ndim,nxvh), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine SMOOTH1(q,qs,ffc,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, dimension(2*nxvh), intent(in) :: q
         complex, dimension(nxvh), intent(inout) :: qs
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine RDMODES1(pot,pott,nx,modesx,nxvh,modesxd)
         implicit none
         integer, intent(in) :: nx, modesx, nxvh, modesxd
         complex, dimension(nxvh), intent(in) :: pot
         complex, dimension(modesxd), intent(inout) :: pott
         end subroutine
      end interface
!
      interface
         subroutine WRMODES1(pot,pott,nx,modesx,nxvh,modesxd)
         implicit none
         integer, intent(in) :: nx, modesx, nxvh, modesxd
         complex, dimension(nxvh), intent(inout) :: pot
         complex, dimension(modesxd), intent(in) :: pott
         end subroutine
      end interface
!
      end module
