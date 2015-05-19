!-----------------------------------------------------------------------
! Interface file for mpplib2.f
!     module mpplib2_h
!     implicit none
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
!        real, dimension(nxp), intent(inout) :: f, g
         real, dimension(*), intent(inout) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine PPDSUM(f,g,nxp)
         implicit none
         integer, intent(in) :: nxp
!        double precision, dimension(nxp), intent(inout) :: f, g
         double precision, dimension(*), intent(inout) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine PPIMAX(if,ig,nxp)
         implicit none
         integer, intent(in) :: nxp
!        integer, dimension(nxp), intent(inout) :: if, ig
         integer, dimension(*), intent(inout) :: if, ig
         end subroutine
      end interface
!
      interface
         subroutine PPDMAX(f,g,nxp)
         implicit none
         integer, intent(in) :: nxp
!        double precision, dimension(nxp), intent(inout) :: f, g
         double precision, dimension(*), intent(inout) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine PPNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
         implicit none
         integer, intent(in) :: nyp, kstrt, nvp, nxv, nypmx
!        real, dimension(nxv,nypmx), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine PPNAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
         implicit none
         integer, intent(in) :: nyp, kstrt, nvp, nx, nxv, nypmx
!        real, dimension(nxv,nypmx), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        real, dimension(nxv), intent(inout) :: scr
         real, dimension(*), intent(inout) :: scr
         end subroutine
      end interface
!
      interface
         subroutine PPNACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
         implicit none
         integer, intent(in) :: ndim, nyp, kstrt, nvp, nx, nxv, nypmx
!        real, dimension(ndim,nxv,nypmx), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        real, dimension(ndim,nxv), intent(inout) :: scr
         real, dimension(*), intent(inout) :: scr
         end subroutine
      end interface
!
      interface
         subroutine PPTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd&
     &,kypd)
         implicit none
         integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv
         integer, intent(in) :: kxpd, kypd
!        real, dimension(2*nxv,kypd), intent(in) :: f
         real, dimension(*), intent(in) :: f
!        complex, dimension(nyv,kxpd), intent(inout) :: g
         complex, dimension(*), intent(inout) :: g
!        complex, dimension(kxp*kyp), intent(inout) :: s, t
         complex, dimension(*), intent(inout) :: s, t
         end subroutine
      end interface
!
      interface
         subroutine PPNTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,  &
     &nyv,kxpd,kypd)
         implicit none
         integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
         integer, intent(in) :: nxv, nyv, kxpd, kypd
!        real, dimension(ndim,2*nxv,kypd), intent(in) :: f
         real, dimension(*), intent(in) :: f
!        complex, dimension(ndim,nyv,kxpd), intent(inout) :: g
         complex, dimension(*), intent(inout) :: g
!        complex, dimension(ndim,kxp*kyp), intent(inout) :: s, t
         complex, dimension(*), intent(inout) :: s, t
         end subroutine
      end interface
!
      interface
         subroutine PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr&
     &,kstrt,nvp,idimp,nbmax,mx1)
         implicit none
         integer, intent(in) :: kstrt, nvp, idimp, nbmax, mx1
!        real, dimension(idimp,nbmax), intent(in) :: sbufr, sbufl
         real, dimension(*), intent(in) :: sbufr, sbufl
!        real, dimension(idimp,nbmax), intent(inout) :: rbufr, rbufl
         real, dimension(*), intent(inout) :: rbufr, rbufl
!        real, dimension(3,mx1), intent(in) :: ncll, nclr
         real, dimension(*), intent(in) :: ncll, nclr
!        real, dimension(3,mx1), intent(inout) :: mcll, mclr
         real, dimension(*), intent(inout) :: mcll, mclr
         end subroutine
      end interface
!
!     end module
      end

