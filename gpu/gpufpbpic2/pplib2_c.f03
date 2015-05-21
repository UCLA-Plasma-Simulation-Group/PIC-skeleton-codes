!-----------------------------------------------------------------------
! Fortran2003 Basic parallel PIC library for MPI communications
! interface to C
! written by Viktor K. Decyk, UCLA
      module pplib2_c
      use iso_c_binding
      use pplib2
      implicit none
      private
!
      contains
!
!-----------------------------------------------------------------------
      subroutine fppinit2(idproc,nvp,argc,argv) bind(C,name='cppinit2')
      implicit none
      integer(c_int) :: idproc, nvp
      integer(c_int), value :: argc
      type (c_ptr) :: argv
      call PPINIT2(idproc,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fppfndgrp(locl,kstrt,nvp,idev,ndev)                    &
     &bind(C,name='cppfndgrp')
      implicit none
      integer(c_int), value :: kstrt, nvp
      integer(c_int) :: idev, ndev
      type (c_ptr), value :: locl
! local data
      integer, dimension(:), pointer :: kocl
      call c_f_pointer(locl,kocl,(/nvp/))
      call PPFNDGRP(kocl,kstrt,nvp,idev,ndev)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fppexit() bind(C,name='cppexit')
      implicit none
      call PPEXIT()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fppabort() bind(C,name='cppabort')
      implicit none
      call PPABORT()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fpwtimera(icntrl,time,dtime) bind(C,name='cpwtimera')
      implicit none
      integer(c_int), value :: icntrl
      real(c_float) :: time
      real(c_double) :: dtime
      call PWTIMERA(icntrl,time,dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fppsum(f,g,nxp) bind(C,name='cppsum')
      implicit none
      integer(c_int), value :: nxp
      type (c_ptr), value :: f, g
! local data
      real, dimension(:), pointer :: a, b
      call c_f_pointer(f,a,(/nxp/))
      call c_f_pointer(g,b,(/nxp/))
      call PPSUM(a,b,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fppdsum(f,g,nxp) bind(C,name='cppdsum')
      implicit none
      integer(c_int), value :: nxp
      type (c_ptr), value :: f, g
! local data
      double precision, dimension(:), pointer :: a, b
      call c_f_pointer(f,a,(/nxp/))
      call c_f_pointer(g,b,(/nxp/))
      call PPDSUM(a,b,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fppimax(if,ig,nxp) bind(C,name='cppimax')
      implicit none
      integer(c_int), value :: nxp
      type (c_ptr), value :: if, ig
! local data
      integer, dimension(:), pointer :: ia, ib
      call c_f_pointer(if,ia,(/nxp/))
      call c_f_pointer(ig,ib,(/nxp/))
      call PPIMAX(ia,ib,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fpppcncguard2l(scs,scr,kstrt,nvp,nxvh)                 &
     &bind(C,name='cpppcncguard2l')
      implicit none
      integer(c_int), value :: kstrt, nvp, nxvh
      type (c_ptr), value :: scs, scr
! local data
      complex, dimension(:), pointer :: rcs, rcr
      call c_f_pointer(scs,rcs,(/nxvh/))
      call c_f_pointer(scr,rcr,(/nxvh/))
      call PPPCNCGUARD2L(rcs,rcr,kstrt,nvp,nxvh)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fpppcnaguard2l(scs,scr,kstrt,nvp,nxvh)                 &
     &bind(C,name='cpppcnaguard2l')
      implicit none
      integer(c_int), value :: kstrt, nvp, nxvh
      type (c_ptr), value :: scs, scr
! local data
      complex, dimension(:), pointer :: rcs, rcr
      call c_f_pointer(scs,rcs,(/nxvh/))
      call c_f_pointer(scr,rcr,(/nxvh/))
      call PPPCNAGUARD2L(rcs,rcr,kstrt,nvp,nxvh)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fppptpose(sm,tm,nx,ny,kxp,kyp,kstrt,nvp)               &
     &bind(C,name='cppptpose')
      implicit none
      integer(c_int), value :: nx, ny, kxp, kyp, kstrt, nvp
      type (c_ptr), value :: sm, tm
! local data
      complex, dimension(:,:), pointer :: um, vm
      call c_f_pointer(sm,um,(/kxp*kyp,nvp/))
      call c_f_pointer(tm,vm,(/kxp*kyp,nvp/))
      call PPPTPOSE(um,vm,nx,ny,kxp,kyp,kstrt,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fppptposen(sm,tm,nx,ny,kxp,kyp,kstrt,nvp,ndim)         &
     &bind(C,name='cppptposen')
      implicit none
      integer(c_int), value :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      type (c_ptr), value :: sm, tm
! local data
      complex, dimension(:,:), pointer :: um, vm
      call c_f_pointer(sm,um,(/kxp*ndim*kyp,nvp/))
      call c_f_pointer(tm,vm,(/kxp*ndim*kyp,nvp/))
      call PPPTPOSEN(um,vm,nx,ny,kxp,kyp,kstrt,nvp,ndim)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine facsndrec(stm,idproc,nsize,ntag,mode)                  &
     &bind(C,name='cacsndrec')
      implicit none
      integer(c_int), value :: idproc, nsize, ntag, mode
      type (c_ptr), value :: stm
! local data
      complex, dimension(:), pointer :: rtm
      call c_f_pointer(stm,rtm,(/nsize/))
      call ACSNDREC(rtm,idproc,nsize,ntag,mode)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine fpppmove2(sbufr,sbufl,rbufr, rbufl,ncll,nclr,mcll,mclr,&
     &kstrt,nvp,idimp,nbmax,mx1) bind(C,name='cpppmove2')
      implicit none
      integer(c_int), value :: kstrt, nvp, idimp, nbmax, mx1
      type (c_ptr), value :: sbufr, sbufl, rbufr, rbufl
      type (c_ptr), value :: ncll, nclr, mcll, mclr
! local data
      real, dimension(:), pointer :: tbufr, tbufl, ubufr, ubufl
      integer, dimension(:,:), pointer :: jcll, jclr, kcll, kclr
      call c_f_pointer(sbufr,tbufr,(/idimp*nbmax/))
      call c_f_pointer(sbufl,tbufl,(/idimp*nbmax/))
      call c_f_pointer(rbufr,ubufr,(/idimp*nbmax/))
      call c_f_pointer(rbufl,ubufl,(/idimp*nbmax/))
      call c_f_pointer(ncll,jcll,(/3,mx1/))
      call c_f_pointer(nclr,jclr,(/3,mx1/))
      call c_f_pointer(mcll,kcll,(/3,mx1/))
      call c_f_pointer(mclr,kclr,(/3,mx1/))
      call PPPMOVE2(tbufr,tbufl,ubufr,ubufl,jcll,jclr,kcll,kclr,kstrt,  &
     &nvp,idimp,nbmax,mx1)
         end subroutine
!
      end module

