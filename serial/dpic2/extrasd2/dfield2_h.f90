!-----------------------------------------------------------------------
! Interface file for dfield2.f
      module dfield2_h
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
         subroutine CURLF2(f,g,nx,ny,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv
         complex, dimension(3,nxvh,nyv), intent(in) :: f
         complex, dimension(3,nxvh,nyv), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine APOTP23(cu,axy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(3,2*nxvh,nyv), intent(in) :: cu
         complex, dimension(3,nxvh,nyv), intent(inout) :: axy
         complex, dimension(nxhd,nyhd), intent(in) :: ffc
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
         subroutine SMOOTH23(cu,cus,ffc,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv, nxhd, nyhd
         real, dimension(3,2*nxvh,nyv), intent(in) :: cu
         complex, dimension(3,nxvh,nyv), intent(inout) :: cus
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
      interface
         subroutine RDVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,nxvh, &
     &nyv,modesxd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, ndim, nxvh, nyv
         integer, intent(in) :: modesxd, modesyd
         complex, dimension(ndim,nxvh,nyv), intent(in) :: vpot
         complex, dimension(ndim,modesxd,modesyd), intent(inout) ::     &
     &   vpott
         end subroutine
      end interface
!
      interface
         subroutine WRVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,nxvh, &
     &nyv,modesxd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, ndim, nxvh, nyv
         integer, intent(in) :: modesxd, modesyd
         complex, dimension(ndim,nxvh,nyv), intent(inout) :: vpot
         complex, dimension(ndim,modesxd,modesyd), intent(in) :: vpott
         end subroutine
      end interface
!
      end module
