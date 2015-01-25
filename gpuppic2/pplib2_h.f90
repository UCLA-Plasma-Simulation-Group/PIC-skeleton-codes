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
         subroutine PPFNDGRP(locl,kstrt,nvp,idev,ndev)
         implicit none
         integer :: kstrt, nvp, idev, ndev
         integer, dimension(nvp) :: locl
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
         subroutine PPPNCGUARD2L(scs,scr,kstrt,nvp,nxv)
         implicit none
         integer :: kstrt, nvp, nxv
         real, dimension(nxv) :: scs, scr
         end subroutine
      end interface
!
      interface
         subroutine PPPNAGUARD2L(scs,scr,kstrt,nvp,nxv)
         implicit none
         integer :: kstrt, nvp, nxv
         real, dimension(nxv) :: scs, scr
         end subroutine
      end interface
!
      interface
         subroutine PPPTPOSE(sm,tm,nx,ny,kxp,kyp,kstrt,nvp)
         implicit none
         integer :: nx, ny, kxp, kyp, kstrt, nvp
         complex, dimension(kxp*kyp,nvp) :: sm, tm
         end subroutine
      end interface
!
      interface
         subroutine PPPTPOSEN(sm,tm,nx,ny,kxp,kyp,kstrt,nvp,ndim)
         implicit none
         integer :: nx, ny, kxp, kyp, kstrt, nvp, ndim
         complex, dimension(kxp*ndim*kyp,nvp) :: sm, tm
         end subroutine
      end interface
!
      interface
         subroutine ACSNDREC(stm,idproc,nsize,ntag,mode)
         implicit none
         integer :: idproc, nsize, ntag, mode
         complex, dimension(nsize) :: stm
         end subroutine
      end interface
!
      interface
         subroutine PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr&
     &,kstrt,nvp,idimp,nbmax,mx1)
         implicit none
         integer :: kstrt, nvp, idimp, nbmax, mx1
         real, dimension(idimp,nbmax) :: sbufr, sbufl, rbufr, rbufl
         integer, dimension(3,mx1) :: ncll, nclr, mcll, mclr
         end subroutine
      end interface
!
      end module

