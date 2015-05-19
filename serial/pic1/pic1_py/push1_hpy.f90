!-----------------------------------------------------------------------
! Interface file for push1.f
!     module push1_h
!     implicit none
!
      interface
         subroutine DISTR1(part,vtx,vdx,npx,idimp,nop,nx,ipbc)
         implicit none
         integer, intent(in) :: npx, idimp, nop, nx, ipbc
         real, intent(in) :: vtx, vdx
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine GPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, nxv, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(nxv), intent(in) :: fx
         real, dimension(*), intent(in) :: fx
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
         subroutine CGUARD1L(fx,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
!        real, dimension(nxe), intent(inout) :: fx
         real, dimension(*), intent(inout) :: fx
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

