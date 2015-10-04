!-----------------------------------------------------------------------
! Interface file for pdfield3.f
      module pdfield3_h
      implicit none
!
      interface
         subroutine PPOTP32(q,pot,ffc,we,nx,ny,nz,kstrt,nvpy,nvpz,nzv,  &
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
         subroutine PPDIVF32(f,df,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp&
     &)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: f
         complex, dimension(nzv,kxyp,kyzp), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine PPGRADF32(df,f,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,   &
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
         subroutine PPCURLF32(f,g,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp&
     &)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp,kyzp
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: f
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine PPAPOTP32(cu,axyz,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,nvpz&
     &,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp, nzhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: cu
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: axyz
         complex, dimension(nzhd,kxyp,kyzp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine PPETFIELD332(dcu,exyz,ffe,affp,ci,wf,nx,ny,nz,kstrt,& 
     &nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kstrt, nvpy, nvpz, nzv
         integer, intent(in) :: kxyp, kyzp, nzhd
         real, intent(in) :: affp, ci
         real, intent(inout) :: wf
         complex, dimension(3,nzv,kxyp,kyzp), intent(in) :: dcu
         complex, dimension(3,nzv,kxyp,kyzp), intent(inout) :: exyz
         complex, dimension(nzhd,kxyp,kyzp), intent(inout) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine PPSMOOTH32(q,qs,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,   &
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
         subroutine PPSMOOTH332(cu,cus,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,&
     &kxyp,kyzp,nzhd)
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
