!-----------------------------------------------------------------------
! Interface file for field2.f
      module field2_h
      implicit none
!
      interface
         subroutine POTP2(q,pot,ffc,we,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(inout) :: we
         real, dimension(2*nxvh,nyv), intent(in) :: q
         complex, dimension(nxvh,nyv), intent(inout) :: pot
         complex, dimension(nxhd,nyhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine DIVF2(f,df,nx,ny,ndim,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, ndim, nxvh, nyv
         complex, dimension(3,nxvh,nyv), intent(in) :: f
         complex, dimension(nxvh,nyv), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine GRADF2(df,f,nx,ny,ndim,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, ndim, nxvh, nyv
         complex, dimension(nxvh,nyv), intent(in) :: df
         complex, dimension(3,nxvh,nyv), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine SMOOTH2(q,qs,ffc,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv, nxhd, nyhd
         real, dimension(2*nxvh,nyv), intent(in) :: q
         complex, dimension(nxvh,nyv), intent(inout) :: qs
         complex, dimension(nxhd,nyhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine RDMODES2(pot,pott,nx,ny,modesx,modesy,nxvh,nyv,     &
     &modesxd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, nxvh, nyv
         integer, intent(in) :: modesxd, modesyd
         complex, dimension(nxvh,nyv), intent(in) :: pot
         complex, dimension(modesxd,modesyd), intent(inout) :: pott
         end subroutine
      end interface
!
      interface
         subroutine WRMODES2(pot,pott,nx,ny,modesx,modesy,nxvh,nyv,     &
     &modesxd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, nxvh, nyv
         integer, intent(in) :: modesxd, modesyd
         complex, dimension(nxvh,nyv), intent(inout) :: pot
         complex, dimension(modesxd,modesyd), intent(in) :: pott
         end subroutine
      end interface
!
      end module
