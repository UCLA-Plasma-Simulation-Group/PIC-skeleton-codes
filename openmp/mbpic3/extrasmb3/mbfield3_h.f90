!-----------------------------------------------------------------------
! Interface file for mbfield3.f
      module mbfield3_h
      implicit none
!
      interface
         subroutine MPOTP3(q,pot,ffc,we,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,&
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
         subroutine MDIVF3(f,df,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         complex, dimension(3,nxvh,nyv,nzv), intent(in) :: f
         complex, dimension(nxvh,nyv,nzv), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine MGRADF3(df,f,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         complex, dimension(nxvh,nyv,nzv), intent(in) :: df
         complex, dimension(3,nxvh,nyv,nzv), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine MCURLF3(f,g,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         complex, dimension(3,nxvh,nyv,nzv), intent(in) :: f
         complex, dimension(3,nxvh,nyv,nzv), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine MAVPOT33(bxyz,axyz,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         complex, dimension(3,nxvh,nyv,nzv), intent(in) :: bxyz
         complex, dimension(3,nxvh,nyv,nzv), intent(inout) :: axyz
         end subroutine
      end interface
!
      interface
         subroutine MAVRPOT33(axyz,bxyz,ffc,ci,nx,ny,nz,nxvh,nyv,nzv,   &
     &nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ci
         complex, dimension(3,nxvh,nyv,nzv), intent(inout) :: axyz
         complex, dimension(3,nxvh,nyv,nzv), intent(in) :: bxyz
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MSMOOTH3(q,qs,ffc,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,  &
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
         subroutine MSMOOTH33(cu,cus,ffc,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd&
     &,nzhd)
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
