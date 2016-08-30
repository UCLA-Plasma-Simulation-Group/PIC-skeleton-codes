!-----------------------------------------------------------------------
! Interface file for dfield1.f
      module dfield1_h
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
         subroutine CURLF1(f,g,nx,nxvh)
         implicit none
         integer, intent(in) :: nx, nxvh
         complex, dimension(2,nxvh), intent(in) :: f
         complex, dimension(2,nxvh), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine APOTP13(cu,ayz,ffc,ci,wm,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(2,2*nxvh), intent(in) :: cu
         complex, dimension(2,nxvh), intent(inout) :: ayz
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ETFIELD13(dcu,eyz,ffe,ci,wf,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci
         real, intent(inout) :: wf
         real, dimension(3,2*nxvh), intent(in) :: dcu
         complex, dimension(3,nxvh), intent(inout) :: eyz
         complex, dimension(nxhd), intent(in) :: ffe
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
         subroutine SMOOTH13(cu,cus,ffc,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, dimension(2,2*nxvh), intent(in) :: cu
         complex, dimension(2,nxvh), intent(inout) :: cus
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
      interface
         subroutine RDVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
         implicit none
         integer, intent(in) :: nx, modesx, ndim, nxvh, modesxd
         complex, dimension(ndim,nxvh), intent(in) :: vpot
         complex, dimension(ndim,modesxd), intent(inout) ::     &
     &   vpott
         end subroutine
      end interface
!
      interface
         subroutine WRVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
         implicit none
         integer, intent(in) :: nx, modesx, ndim, nxvh, modesxd
         complex, dimension(ndim,nxvh), intent(inout) :: vpot
         complex, dimension(ndim,modesxd), intent(in) :: vpott
         end subroutine
      end interface
!
      end module
