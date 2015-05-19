!-----------------------------------------------------------------------
! Interface file for push2.c
!     module cpush2_h
!     implicit none
!
      interface
         subroutine cdistr2(part,vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx,ny&
     &,ipbc)
         implicit none
         integer, intent(in) :: npx, npy, idimp, nop, nx, ny, ipbc
         real, intent(in) :: vtx, vty, vdx, vdy
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine cgpush2l(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,&
     &ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(2,nxv,nyv), intent(in) :: fxy
         real, dimension(*), intent(in) :: fxy
         end subroutine
      end interface
!
      interface
         subroutine cgpost2l(part,q,qm,nop,idimp,nxv,nyv)
         implicit none
         integer, intent(in) :: nop, idimp, nxv, nyv
         real, intent(in) :: qm
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        real, dimension(nxv,nyv), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine cdsortp2yl(parta,partb,npic,idimp,nop,ny1)
         implicit none
         integer, intent(in) :: idimp, nop, ny1
!        real, dimension(idimp,nop), intent(in) :: parta
         real, dimension(*), intent(in) :: parta
!        real, dimension(idimp,nop), intent(inout) :: partb
         real, dimension(*), intent(inout) :: partb
!        integer, dimension(ny1), intent(inout) :: npic
         integer, dimension(*), intent(inout) :: npic
         end subroutine
      end interface
!
      interface
         subroutine ccguard2l(fxy,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
!        real, dimension(2,nxe,nye), intent(inout) :: fxy
         real, dimension(*), intent(inout) :: fxy
         end subroutine
      end interface
!
      interface
         subroutine caguard2l(q,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
!        real, dimension(nxe,nye), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine cpois22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv&
     &,nxhd,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
!        real, dimension(2*nxvh,nyv), intent(in) :: q
         real, dimension(*), intent(in) :: q
!        real, dimension(2,2*nxvh,nyv), intent(inout) :: fxy
         real, dimension(*), intent(inout) :: fxy
!        complex, dimension(nxhd,nyhd), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cwfft2rinit(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: indx, indy, nxhyd, nxyhd
!        integer, dimension(nxhyd), intent(inout) :: mixup
         integer, dimension(*), intent(inout) :: mixup
!        complex, dimension(nxyhd), intent(inout) :: sct
         complex, dimension(*), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cwfft2rx(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
!        real, dimension(2*nxhd,nyd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(inout) :: mixup
         integer, dimension(*), intent(inout) :: mixup
!        complex, dimension(nxyhd), intent(inout) :: sct
         complex, dimension(*), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cwfft2r2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
!        real, dimension(2,2*nxhd,nyd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(inout) :: mixup
         integer, dimension(*), intent(inout) :: mixup
!        complex, dimension(nxyhd), intent(inout) :: sct
         complex, dimension(*), intent(inout) :: sct
         end subroutine
      end interface
!
!     end module
      end

