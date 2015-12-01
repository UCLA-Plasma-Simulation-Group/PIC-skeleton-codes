!-----------------------------------------------------------------------
! Interface file for dfield3.f
      module dfield3_h
      implicit none
!
      interface
         subroutine POTP3(q,pot,ffc,we,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd, &
     &nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(inout) :: we
         real, dimension(2*nxvh,nyv,nzv), intent(in) :: q
         complex, dimension(nxvh,nyv,nzv), intent(inout) :: pot
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine DIVF3(f,df,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         complex, dimension(3,nxvh,nyv,nzv), intent(in) :: f
         complex, dimension(nxvh,nyv,nzv), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine GRADF3(df,f,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         complex, dimension(nxvh,nyv,nzv), intent(in) :: df
         complex, dimension(3,nxvh,nyv,nzv), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine CURLF3(f,g,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         complex, dimension(3,nxvh,nyv,nzv), intent(in) :: f
         complex, dimension(3,nxvh,nyv,nzv), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine APOTP33(cu,axyz,ffc,ci,wm,nx,ny,nz,nxvh,nyv,nzv,nxhd&
     &,nyhd,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(3,2*nxvh,nyv,nzv), intent(in) :: cu
         complex, dimension(3,nxvh,nyv,nzv), intent(inout) :: axyz
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ETFIELD33(dcu,exyz,ffe,ci,wf,nx,ny,nz,nxvh,nyv,nzv, &
     &nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ci
         real, intent(inout) :: wf
         real, dimension(3,2*nxvh,nyv,nzv), intent(in) :: dcu
         complex, dimension(3,nxvh,nyv,nzv), intent(inout) :: exyz
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffe
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
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
         integer, intent(in) :: nxvh, nyv, nzv
         integer, intent(in) :: modesxd, modesyd, modeszd
         complex, dimension(nxvh,nyv,nzv), intent(in) :: pot
         complex, dimension(modesxd,modesyd,modeszd), intent(inout) ::  &
     &   pott
         end subroutine
      end interface
!
      interface
         subroutine WRMODES3(pot,pott,nx,ny,nz,modesx,modesy,modesz,nxvh&
     &,nyv,nzv,modesxd,modesyd,modeszd)
         implicit none
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
         integer, intent(in) :: nxvh, nyv, nzv
         integer, intent(in) :: modesxd, modesyd, modeszd
         complex, dimension(nxvh,nyv,nzv), intent(inout) :: pot
         complex, dimension(modesxd,modesyd,modeszd), intent(in) :: pott
         end subroutine
      end interface
!
      interface
         subroutine RDVMODES3(vpot,vpott,nx,ny,nz,modesx,modesy,modesz, &
     &ndim,nxvh,nyv,nzv,modesxd,modesyd,modeszd)
         implicit none
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz, ndim
         integer, intent(in) :: nxvh, nyv, nzv
         integer, intent(in) :: modesxd, modesyd, modeszd
         complex, dimension(ndim,nxvh,nyv,nzv), intent(in) :: vpot
         complex, dimension(ndim,modesxd,modesyd,modeszd), intent(inout)&
     & :: vpott
         end subroutine
      end interface
!
      interface
         subroutine WRVMODES3(vpot,vpott,nx,ny,nz,modesx,modesy,modesz, &
     &ndim,nxvh,nyv,nzv,modesxd,modesyd,modeszd)
         implicit none
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz, ndim
         integer, intent(in) :: nxvh, nyv, nzv, modesxd, modesyd, modeszd
         complex, dimension(ndim,nxvh,nyv,nzv), intent(inout) :: vpot
         complex, dimension(ndim,modesxd,modesyd,modeszd), intent(in) ::&
     &   vpott
         end subroutine
      end interface
!
      end module
