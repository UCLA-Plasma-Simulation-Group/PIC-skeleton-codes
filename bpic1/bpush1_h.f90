!-----------------------------------------------------------------------
! Interface file for bpush1.f
      module bpush1_h
      implicit none
!
      interface
         subroutine DISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,idimp,nop, &
     &nx,ipbc)
         implicit none
         integer, intent(in) :: npx, idimp, nop, nx, ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop), intent(inout) :: part
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
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         end subroutine
      end interface
!
      interface
         subroutine GRBPUSH13L(part,fxyz,byz,omx,qbm,dt,dtc,ci,ek,idimp,&
     &nop,nx,nxv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, nxv, ipbc
         real, intent(in) :: omx, qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         end subroutine
      end interface
!
      interface
         subroutine GPOST1L(part,q,qm,nop,idimp,nxv)
         implicit none
         integer, intent(in) :: nop, idimp, nxv
         real, intent(in) :: qm
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(nxv), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine GJPOST1L(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
         implicit none
         integer, intent(in) :: nop, idimp, nx, nxv,ipbc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(2,nxv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine GRJPOST1L(part,cu,qm,dt,ci,nop,idimp,nx,nxv,ipbc)
         implicit none
         integer, intent(in) :: nop, idimp, nx, nxv, ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(2,nxv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine DSORTP1XL(parta,partb,npic,idimp,nop,nx1)
         implicit none
         integer, intent(in) :: idimp, nop, nx1
         real, dimension(idimp,nop), intent(in) :: parta
         real, dimension(idimp,nop), intent(inout) :: partb
         integer, dimension(nx1), intent(inout) :: npic
         end subroutine
      end interface
!
      interface
         subroutine CGUARD1L(byz,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(2,nxe), intent(inout) :: byz
         end subroutine
      end interface
!
      interface
         subroutine BGUARD1L(fxyz,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(3,nxe), intent(inout) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine ACGUARD1L(cu,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(2,nxe), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine AGUARD1L(q,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(nxe), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine POIS1(q,fx,isign,ffc,ax,affp,we,nx)
         implicit none
         integer, intent(in) :: isign, nx
         real, intent(in) :: ax, affp
         real, intent(inout) :: we
         real, dimension(nx), intent(in) :: q
         real, dimension(nx), intent(inout) :: fx
         complex, dimension(nx/2), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine IBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(2,2*nxvh), intent(in) :: cu
         complex, dimension(2,nxvh), intent(inout) :: byz
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MAXWEL1(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(2,nxvh), intent(inout) :: eyz, byz
         real, dimension(2,2*nxvh), intent(in) :: cu
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine EMFIELD1(fxyz,fx,eyz,ffc,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, dimension(3,2*nxvh), intent(inout) :: fxyz
         real, dimension(2*nxvh), intent(in) :: fx
         complex, dimension(3,nxvh), intent(in) :: eyz
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine BMFIELD1(fyz,eyz,ffc,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, dimension(2,2*nxvh), intent(inout) :: fyz
         complex, dimension(2,nxvh), intent(in) :: eyz
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine WFFT1RINIT(mixup,sct,indx,nxhd)
         implicit none
         integer, intent(in) :: indx, nxhd
         integer, dimension(nxhd), intent(inout) :: mixup
         complex, dimension(nxhd), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer, intent(in) :: isign, indx, nxd, nxhd
         real, dimension(nxd), intent(inout) :: f, t
         integer, dimension(nxhd), intent(in) :: mixup
         complex, dimension(nxhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT1R2X(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer, intent(in) :: isign, indx, nxd, nxhd
         real, dimension(2,nxd), intent(inout) :: f, t
         integer, dimension(nxhd), intent(in) :: mixup
         complex, dimension(nxhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT1R3X(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer , intent(in):: isign, indx, nxd, nxhd
         real, dimension(3,nxd), intent(inout) :: f, t
         integer, dimension(nxhd), intent(in) :: mixup
         complex, dimension(nxhd), intent(in) :: sct
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
      end module
