!-----------------------------------------------------------------------
! Interface file for mpbfield3.f
      module mpbfield3_h
      implicit none
!
      interface
         subroutine MPPOTP32(q,pot,ffc,we,nx,ny,nz,kstrt,nvpy,nvpz,nzv, &
     &kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp, nzhd
         real, intent(inout) :: we
         complex, dimension(nzv,kxyp,kyzp), intent(in) :: q
         complex, dimension(nzv,kxyp,kyzp), intent(inout) :: pot
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPDIVF32(f,df,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,   &
     &kyzp)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: f
         complex, dimension(nzv,kxyp,kyzp), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine MPPGRADF32(df,f,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,  &
     &kyzp)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp
         complex, dimension(nzv,kxyp,kyzp), intent(in) :: df
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine MPPCURLF32(f,g,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,   &
     &kyzp)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp,kyzp
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: f
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine MPPAVPOT332(bxyz,axyz,nx,ny,nz,kstrt,nvpy,nvpz,nzv, &
     &kxyp,kyzp)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: bxyz
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: axyz
         end subroutine
      end interface
!
      interface
         subroutine MPPAVRPOT332(axyz,bxyz,ffc,affp,ci,nx,ny,nz,kstrt,  &
     &nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz
         integer, intent(in) :: nzv, kxyp, kyzp, nzhd
         real, intent(in) :: affp, ci
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: axyz
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: bxyz
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTH32(q,qs,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,  &
     &kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp, nzhd
         complex, dimension(nzv,kxyp,kyzp), intent(in) :: q
         complex, dimension(nzv,kxyp,kyzp), intent(inout) :: qs
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTH332(cu,cus,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv&
     &,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp, nzhd
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: cu
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: cus
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine PPRDMODES32(pot,pott,nx,ny,nz,modesx,modesy,modesz, &
     &kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
         implicit none
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
         integer, intent(in) :: kstrt, nvpy, nvpz, nzv, kxyp, kyzp
         integer, intent(in) :: modesxpd, modesypd, modeszd
         complex, dimension(nzv,kxyp,kyzp), intent(in) :: pot
         complex, dimension(modeszd,modesxpd,modesypd), intent(inout) ::&
     & pott
         end subroutine
      end interface
!
      interface
         subroutine PPWRMODES32(pot,pott,nx,ny,nz,modesx,modesy,modesz, &
     &kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
         implicit none
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
         integer, intent(in) :: kstrt, nvpy, nvpz, nzv, kxyp, kyzp
         integer, intent(in) :: modesxpd, modesypd, modeszd
         complex, dimension(nzv,kxyp,kyzp), intent(inout) :: pot
         complex, dimension(modeszd,modesxpd,modesypd), intent(in) ::   &
     & pott
         end subroutine
      end interface
!
      interface
         subroutine PPRDVMODES32(vpot,vpott,nx,ny,nz,modesx,modesy,     &
     &modesz,ndim,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,      &
     &modeszd)
         implicit none
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
         integer, intent(in) :: ndim, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
         integer, intent(in) :: modesxpd, modesypd, modeszd
         complex, dimension(ndim,nzv,kxyp,kyzp), intent(in) :: vpot
         complex, dimension(ndim,modeszd,modesxpd,modesypd),            &
     &intent(inout) :: vpott
         end subroutine
      end interface
!
      interface
         subroutine PPWRVMODES32(vpot,vpott,nx,ny,nz,modesx,modesy,     &
     &modesz,ndim,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd       &
     &,modeszd)
         implicit none
         integer, intent(in) :: nx, ny, nz, modesx, modesy, modesz
         integer, intent(in) :: ndim, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
         integer, intent(in) :: modesxpd, modesypd, modeszd
         complex, dimension(ndim,nzv,kxyp,kyzp), intent(inout) :: vpot
         complex, dimension(ndim,modeszd,modesxpd,modesypd), intent(in) &
     &:: vpott
         end subroutine
      end interface
!
      end module
