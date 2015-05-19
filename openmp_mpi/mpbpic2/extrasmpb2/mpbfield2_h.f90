!-----------------------------------------------------------------------
! Interface file for mpbfield2.f
      module mpbfield2_h
      implicit none
!
      interface
         subroutine MPPOTP2(q,pot,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
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
         subroutine MPPDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp
         complex, dimension(ndim,nyv,kxp), intent(in) :: f
         complex, dimension(nyv,kxp), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine MPPGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp
         complex, dimension(nyv,kxp), intent(in) :: df
         complex, dimension(ndim,nyv,kxp), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine MPPCURLF2(f,g,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(in) :: f
         complex, dimension(3,nyv,kxp), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine MPPAVPOT23(bxy,axy,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(in) :: bxy
         complex, dimension(3,nyv,kxp), intent(inout) :: axy
         end subroutine
      end interface
!
      interface
         subroutine MPPAVRPOT23(axy,bxy,ffc,affp,ci,nx,ny,kstrt,nyv,kxp,&
     &nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: affp, ci
         complex, dimension(3,nyv,kxp), intent(inout) :: axy
         complex, dimension(3,nyv,kxp), intent(in) :: bxy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTH2(q,qs,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(nyv,kxp), intent(inout) :: qs
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTH23(cu,cus,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         complex, dimension(3,nyv,kxp), intent(in) :: cu
         complex, dimension(3,nyv,kxp), intent(inout) :: cus
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
      interface
         subroutine PPRDVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,    &
     &kstrt,nyv,kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, ndim, kstrt
         integer, intent(in) :: nyv, kxp, modesxpd, modesyd
         complex, dimension(ndim,nyv,kxp), intent(in) :: vpot
         complex, dimension(ndim,modesyd,modesxpd), intent(inout) ::     &
     &   vpott
         end subroutine
      end interface
!
      interface
         subroutine PPWRVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,    &
     &kstrt,nyv,kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, ndim, kstrt
         integer, intent(in) :: nyv, kxp, modesxpd, modesyd
         complex, dimension(ndim,nyv,kxp), intent(inout) :: vpot
         complex, dimension(ndim,modesyd,modesxpd), intent(in) :: vpott
         end subroutine
      end interface
!
      end module
