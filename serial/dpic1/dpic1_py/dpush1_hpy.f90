!-----------------------------------------------------------------------
! Interface file for dpush1.f
!     module dpush1_h
!     implicit none
!
      interface
         subroutine DISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,idimp,nop, &
     &nx,ipbc)
         implicit none
         integer, intent(in) :: npx, idimp, nop, nx, ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine GBPUSH13L(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,nop,&
     &nx,nxv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, nxv, ipbc
         real, intent(in) :: omx, qbm, dt, dtc
         real, intent(inout) :: ek
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(*), intent(in) :: fxyz
!        real, dimension(2,nxv), intent(in) :: byz
         real, dimension(*), intent(in) :: byz
         end subroutine
      end interface
!
      interface
         subroutine GPOST1L(part,q,qm,nop,idimp,nxv)
         implicit none
         integer, intent(in) :: nop, idimp, nxv
         real, intent(in) :: qm
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        real, dimension(nxv), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine GJPOST1L(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
         implicit none
         integer, intent(in) :: nop, idimp, nx, nxv,ipbc
         real, intent(in) :: qm, dt
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(2,nxv), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine GMJPOST1L(part,amu,qm,nop,idimp,nxv)
         implicit none
         integer, intent(in) :: nop, idimp, nxv
         real, intent(in) :: qm
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        real, dimension(2,nxv), intent(inout) :: amu
         real, dimension(*), intent(inout) :: amu
         end subroutine
      end interface
!
      interface
         subroutine GDJPOST1L(part,fxyz,byz,dcu,amu,omx,qm,qbm,dt,idimp,&
     &nop,nxv)
         implicit none
         integer, intent(in) :: idimp, nop, nxv
         real, intent(in) :: omx, qm, qbm, dt
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(*), intent(in) :: fxyz
!        real, dimension(2,nxv), intent(in) :: byz
         real, dimension(*), intent(in) :: byz
!        real, dimension(2,nxv), intent(inout) :: dcu, amu
         real, dimension(*), intent(inout) :: dcu, amu
         end subroutine
      end interface
!
      interface
         subroutine GDCJPOST1L(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt,  &
     &idimp,nop,nxv)
         implicit none
         integer, intent(in) :: idimp, nop, nxv
         real, intent(in) :: omx, qm, qbm, dt
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(*), intent(in) :: fxyz
!        real, dimension(2,nxv), intent(in) :: byz
         real, dimension(*), intent(in) :: byz
!        real, dimension(2,nxv), intent(inout) :: cu, dcu, amu
         real, dimension(*), intent(inout) :: cu, dcu, amu
         end subroutine
      end interface
!
      interface
         subroutine DSORTP1XL(parta,partb,npic,idimp,nop,nx1)
         implicit none
         integer, intent(in) :: idimp, nop, nx1
!        real, dimension(idimp,nop), intent(in) :: parta
         real, dimension(*), intent(in) :: parta
!        real, dimension(idimp,nop), intent(inout) :: partb
         real, dimension(*), intent(inout) :: partb
!        integer, dimension(nx1), intent(inout) :: npic
         integer, dimension(*), intent(inout) :: npic
         end subroutine
      end interface
!
      interface
         subroutine DGUARD1L(fx,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
!        real, dimension(nxe), intent(inout) :: fx
         real, dimension(*), intent(inout) :: fx
         end subroutine
      end interface
!
      interface
         subroutine CGUARD1L(byz,nx,nxe)
         implicit none
         integer , intent(in):: nx, nxe
!        real, dimension(2,nxe), intent(inout) :: byz
         real, dimension(*), intent(inout) :: byz
         end subroutine
      end interface
!
      interface
         subroutine ACGUARD1L(cu,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
!        real, dimension(2,nxe), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine AGUARD1L(q,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
!        real, dimension(nxe), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine ASCFGUARD1L(dcu,cus,q2m0,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: q2m0
!        real, dimension(2,nxe), intent(inout) :: dcu
         real, dimension(*), intent(inout) :: dcu
!        real, dimension(2,nxe), intent(in) :: cus
         real, dimension(*), intent(in) :: cus
         end subroutine
      end interface
!
      interface
         subroutine FWPMINMX1(qe,qbme,wpmax,wpmin,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: qbme
         real, intent(inout) :: wpmax, wpmin
!        real, dimension(nxe), intent(in) :: qe
         real, dimension(*), intent(in) :: qe
         end subroutine
      end interface
!
      interface
         subroutine POIS1(q,fx,isign,ffc,ax,affp,we,nx)
         implicit none
         integer, intent(in) :: isign, nx
         real, intent(in) :: ax, affp
         real, intent(inout) :: we
!        real, dimension(nx), intent(in) :: q
         real, dimension(*), intent(in) :: q
!        real, dimension(nx), intent(inout) :: fx
         real, dimension(*), intent(inout) :: fx
!        complex, dimension(nx/2), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine BBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
!        real, dimension(2,2*nxvh), intent(in) :: cu
         real, dimension(*), intent(in) :: cu
!        real, dimension(2,2*nxvh), intent(inout) :: byz
         real, dimension(*), intent(inout) :: byz
!        complex, dimension(nxhd), intent(in) :: ffc
         complex, dimension(*), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine BADDEXT1(byz,omy,omz,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: omy, omz
!        real, dimension(2,nxe), intent(inout) :: byz
         real, dimension(*), intent(inout) :: byz
         end subroutine
      end interface
!
      interface
         subroutine DCUPERP13(dcu,amu,nx,nxvh)
         implicit none
         integer, intent(in) :: nx, nxvh
!        real, dimension(2,2*nxvh), intent(inout) :: dcu
         real, dimension(*), intent(inout) :: dcu
!        real, dimension(2,2*nxvh), intent(in) :: amu
         real, dimension(*), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine ADCUPERP13(dcu,amu,nx,nxvh)
         implicit none
         integer, intent(in) :: nx, nxvh
!        real, dimension(2,2*nxvh), intent(inout) :: dcu
         real, dimension(*), intent(inout) :: dcu
!        real, dimension(2,2*nxvh), intent(in) :: amu
         real, dimension(*), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine EPOIS13(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,nxvh,&
     &nxhd)
         implicit none
         integer, intent(in) :: isign, nx, nxvh, nxhd
         real, intent(in) :: ax, affp, wp0, ci
         real, intent(inout) :: wf
!        real, dimension(2,2*nxvh), intent(in) :: dcu
         real, dimension(*), intent(in) :: dcu
!        real, dimension(2,2*nxvh), intent(inout) :: eyz
         real, dimension(*), intent(inout) :: eyz
!        complex, dimension(nxhd), intent(inout) :: ffe
         complex, dimension(*), intent(inout) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine ADDVRFIELD13(fxyze,eyze,fxe,nxe)
         implicit none
         integer, intent(in) :: nxe
!        real, dimension(3,nxe), intent(inout) :: fxyze
         real, dimension(*), intent(inout) :: fxyze
!        real, dimension(2,nxe), intent(in) :: eyze
         real, dimension(*), intent(in) :: eyze
!        real, dimension(nxe), intent(in) :: fxe
         real, dimension(*), intent(in) :: fxe
         end subroutine
      end interface
!
      interface
         subroutine WFFT1RINIT(mixup,sct,indx,nxhd)
         implicit none
         integer, intent(in) :: indx, nxhd
!        integer, dimension(nxhd), intent(inout) :: mixup
         integer, dimension(*), intent(inout) :: mixup
!        complex, dimension(nxhd), intent(inout) :: sct
         complex, dimension(*), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer, intent(in) :: isign, indx, nxd, nxhd
!        real, dimension(nxd), intent(inout) :: f, t
         real, dimension(*), intent(inout) :: f, t
!        integer, dimension(nxhd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT1R2X(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer, intent(in) :: isign, indx, nxd, nxhd
!        real, dimension(3,nxd), intent(inout) :: f, t
         real, dimension(*), intent(inout) :: f, t
!        integer, dimension(nxhd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         function ranorm()
         implicit none
         double precision :: ranorm
         end function
      end interface
!
      interface
         function randum()
         implicit none
         double precision :: randum
         end function
      end interface
!
!     end module
      end
