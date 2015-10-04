!-----------------------------------------------------------------------
! Interface file for pplib3.f
      module pplib3_h
      implicit none
!
      interface
         subroutine PPINIT2(idproc,nvp)
         implicit none
         integer, intent(inout) :: idproc, nvp
         end subroutine
      end interface
!
      interface
         subroutine PPEXIT
         implicit none
         end subroutine
      end interface
!
      interface
         subroutine PPABORT
         implicit none
         end subroutine
      end interface
!
      interface
         subroutine PWTIMERA(icntrl,time,dtime)
         implicit none
         integer, intent(in) :: icntrl
         real, intent(inout) :: time
         double precision, intent(inout) :: dtime
         end subroutine
      end interface
!
      interface
         subroutine PPSUM(f,g,nxp)
         implicit none
         integer, intent(in) :: nxp
         real, dimension(nxp), intent(inout) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine PPDSUM(f,g,nxp)
         implicit none
         integer, intent(in) :: nxp
         double precision, dimension(nxp), intent(inout) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine PPIMAX(if,ig,nxp)
         implicit none
         integer, intent(in) :: nxp
         integer, dimension(nxp), intent(inout) :: if, ig
         end subroutine
      end interface
!
      interface
         subroutine PPDMAX(f,g,nxp)
         implicit none
         integer, intent(in) :: nxp
         double precision, dimension(nxp), intent(inout) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine PPNCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,  &
     &nzpmx,idds)
         implicit none
         integer, intent(in) :: kstrt, nvpy, nvpz, nxv, nypmx, nzpmx
         integer, intent(in) :: idds
         real, dimension(nxv,nypmx,nzpmx), intent(inout) :: f
         real, dimension(nxv,nzpmx,2), intent(inout) :: scs
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPNAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv, &
     &nypmx,nzpmx,idds)
         implicit none
         integer, intent(in) :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
         integer, intent(in) :: idds
         real, dimension(nxv,nypmx,nzpmx), intent(inout) :: f
         real, dimension(nxv,nzpmx,2), intent(inout) :: scs
         real, dimension(nxv,nypmx), intent(inout) :: scr
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPNACGUARD32L(f,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx&
     &,nxv,nypmx,nzpmx,idds)
         implicit none
         integer, intent(in) :: ndim, kstrt ,nvpy, nvpz, nx, nxv, nypmx
         integer, intent(in) :: nzpmx, idds
         real, dimension(ndim,nxv,nypmx,nzpmx), intent(inout) :: f
         real, dimension(ndim,nxv,nzpmx,2), intent(inout) :: scs
         real, dimension(ndim,nxv,nypmx), intent(inout) :: scr
         integer, dimension(idds), intent(in) :: nyzp
         end subroutine
      end interface
!
      interface
         subroutine PPTPOS3A(f,g,s,t,nx,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,  &
     &nxv,nyv,kxypd,kypd,kzpd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kxyp, kyp, kzp, kstrt, nvpy
         integer, intent(in) :: nxv, nyv, kxypd, kypd, kzpd
         real, dimension(2*nxv,kypd,kzpd), intent(in) :: f
         complex, dimension(nyv,kxypd,kzpd), intent(inout) :: g
         complex, dimension(kxyp*kyp*kzp), intent(inout) :: s, t
         end subroutine
      end interface
!
      interface
         subroutine PPTPOS3B(g,h,s,t,nx,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy, &
     &nvpz,nyv,nzv,kxypd,kyzpd,kzpd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kxyp, kyzp, kzp, kstrt, nvpy
         integer, intent(in) :: nvpz, nyv, nzv, kxypd, kyzpd, kzpd
         complex, dimension(nyv,kxypd,kzpd), intent(inout) :: g
         complex, dimension(nzv,kxypd,kyzpd), intent(inout) :: h
         complex, dimension(kyzp*kxyp*kzp), intent(inout) :: s, t
         end subroutine
      end interface
!
      interface
         subroutine PPNTPOS3A(f,g,s,t,nx,ny,nz,kxyp,kyp,kzp,kstrt,nvpy, &
     &ndim,nxv,nyv,kxypd,kypd,kzpd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kxyp, kyp, kzp, kstrt, nvpy
         integer, intent(in) :: ndim, nxv, nyv, kxypd, kypd, kzpd
         real, dimension(ndim,2*nxv,kypd,kzpd), intent(in) :: f
         complex, dimension(ndim,nyv,kxypd,kzpd), intent(inout) :: g
         complex, dimension(ndim,kxyp*kyp*kzp), intent(inout) :: s, t
         end subroutine
      end interface
!
      interface
         subroutine PPNTPOS3B(g,h,s,t,nx,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,&
     &nvpz,ndim,nyv,nzv,kxypd,kyzpd,kzpd)
         implicit none
         integer, intent(in) :: nx, ny, nz, kxyp, kyzp, kzp, kstrt, nvpy
         integer, intent(in) :: nvpz, ndim, nyv, nzv, kxypd, kyzpd, kzpd
         complex, dimension(ndim,nyv,kxypd,kzpd), intent(inout) :: g
         complex, dimension(ndim,nzv,kxypd,kyzpd), intent(inout) :: h
         complex, dimension(ndim,kyzp*kxyp*kzp), intent(inout) :: s, t
         end subroutine
      end interface
!
      interface
         subroutine PPMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,    &
     &ihole,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,idps,nbmax,ntmax,info)
         implicit none
         integer, intent(in) :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer, intent(in) :: idps, nbmax, ntmax
         integer, intent(inout) :: npp
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idps), intent(in) :: edges
         real, dimension(idimp,nbmax), intent(inout) :: sbufr, sbufl
         real, dimension(idimp,nbmax), intent(inout) :: rbufr, rbufl
         integer, dimension(ntmax+1,2) , intent(inout):: ihole
         integer, dimension(7), intent(inout) :: info
         end subroutine
      end interface
!
      interface
         subroutine PPMOVEG32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,   &
     &ihole,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,idps,nbmax,ntmax,info)
         implicit none
         integer, intent(in) :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer, intent(in) :: idps, nbmax, ntmax
         integer, intent(inout) :: npp
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idps), intent(in) :: edges
         real, dimension(idimp,nbmax), intent(inout) :: sbufr, sbufl
         real, dimension(idimp,nbmax), intent(inout) :: rbufr, rbufl
         integer, dimension(ntmax+1,2) , intent(inout):: ihole
         integer, dimension(7), intent(inout) :: info
         end subroutine
      end interface
!
      end module

