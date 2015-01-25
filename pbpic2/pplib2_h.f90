!-----------------------------------------------------------------------
! Interface file for pplib2.f
      module pplib2_h
      implicit none
!
      interface
         subroutine PPINIT2(idproc,nvp)
         implicit none
         integer :: idproc, nvp
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
         integer :: icntrl
         real :: time
         double precision :: dtime
         end subroutine
      end interface
!
      interface
         subroutine PPSUM(f,g,nxp)
         implicit none
         integer :: nxp
         real, dimension(nxp) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine PPDSUM(f,g,nxp)
         implicit none
         integer :: nxp
         double precision, dimension(nxp) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine PPIMAX(if,ig,nxp)
         implicit none
         integer :: nxp
         integer, dimension(nxp) :: if, ig
         end subroutine
      end interface
!
      interface
         subroutine PPNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
         implicit none
         integer :: nyp, kstrt, nvp, nxv, nypmx
         real, dimension(nxv,nypmx) :: f
         end subroutine
      end interface
!
      interface
         subroutine PPNAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
         implicit none
         integer :: nyp, kstrt, nvp, nx, nxv, nypmx
         real, dimension(nxv,nypmx) :: f
         real, dimension(nxv) :: scr
         end subroutine
      end interface
!
      interface
         subroutine PPTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd&
     &,kypd)
         implicit none
         integer :: nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv, kxpd, kypd
         real, dimension(2*nxv,kypd) :: f
         complex, dimension(nyv,kxpd) :: g
         complex, dimension(kxp*kyp) :: s, t
         end subroutine
      end interface
!
      interface
         subroutine PPNTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,  &
     &nyv,kxpd,kypd)
         implicit none
         integer :: nx, ny, kxp, kyp, kstrt, nvp, ndim, nxv, nyv
         integer :: kxpd, kypd
         real, dimension(ndim,2*nxv,kypd) :: f
         complex, dimension(ndim,nyv,kxpd) :: g
         complex, dimension(ndim,kxp*kyp) :: s, t
         end subroutine
      end interface
!
      interface
         subroutine PPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole&
     &,ny,kstrt,nvp,idimp,npmax,idps,nbmax,ntmax,info)
         implicit none
         integer :: ny, kstrt, nvp, idimp, npmax, idps, nbmax, ntmax
         integer :: npp
         real, dimension(idimp,npmax) :: part
         real, dimension(idps) :: edges
         real, dimension(idimp,nbmax) :: sbufr, sbufl, rbufr, rbufl
         integer, dimension(ntmax+1) :: ihole
         integer, dimension(5) :: info
         end subroutine
      end interface
!
      end module

