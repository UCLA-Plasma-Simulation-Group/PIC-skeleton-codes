!-----------------------------------------------------------------------
! Interface file for pfield2.f
      module pfield2_h
      implicit none
!
      interface
         subroutine PPOTP2(q,pot,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(inout) :: we
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(nyv,kxp), intent(inout) :: pot
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine PPDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp
         complex, dimension(ndim,nyv,kxp), intent(in) :: f
         complex, dimension(nyv,kxp), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine PPGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp
         complex, dimension(nyv,kxp), intent(in) :: df
         complex, dimension(ndim,nyv,kxp), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine PPSMOOTH2(q,qs,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(nyv,kxp), intent(inout) :: qs
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine PPRDMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,  &
     &kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, kstrt, nyv, kxp
         integer, intent(in) :: modesxpd, modesyd
         complex, dimension(nyv,kxp), intent(in) :: pot
         complex, dimension(modesyd,modesxpd), intent(inout) :: pott
         end subroutine
      end interface
!
      interface
         subroutine PPWRMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,  &
     &kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, kstrt, nyv, kxp
         integer, intent(in) :: modesxpd, modesyd
         complex, dimension(nyv,kxp), intent(inout) :: pot
         complex, dimension(modesyd,modesxpd), intent(in) :: pott
         end subroutine
      end interface
!
      end module
