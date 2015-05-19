!-----------------------------------------------------------------------
! Interface file for bpush2.c
!     module cbpush2_h
!     implicit none
!
      interface
         subroutine cdistr2h(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp,&
     &nop,nx,ny,ipbc)
         implicit none
         integer, intent(in) :: npx, npy, idimp, nop, nx, ny, ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine cgbpush23l(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx, &
     &ny,nxv,nyv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         real, dimension(*), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine cgrbpush23l(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ipbc)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         real, dimension(*), intent(in) :: fxy, bxy
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
         subroutine cgjpost2l(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc&
     &)
         implicit none
         integer, intent(in) :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qm, dt
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(3,nxv,nyv), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine cgrjpost2l(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,&
     &ipbc)
         implicit none
         integer, intent(in) :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qm, dt, ci
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(3,nxv,nyv), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
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
         subroutine cbguard2l(bxy,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
!        real, dimension(3,nxe,nye), intent(inout)  :: bxy
         real, dimension(*), intent(inout)  :: bxy
         end subroutine
      end interface
!
      interface
         subroutine cacguard2l(cu,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
!        real, dimension(3,nxe,nye), intent(inout)  :: cu
         real, dimension(*), intent(inout)  :: cu
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
         subroutine cpois23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv&
     &,nxhd,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
!        real, dimension(2*nxvh,nyv), intent(in) :: q
         real, dimension(*), intent(in) :: q
!        real, dimension(3,2*nxvh,nyv), intent(inout) :: fxy
         real, dimension(*), intent(inout) :: fxy
!        complex, dimension(nxhd,nyhd), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ccuperp2(cu,nx,ny,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv
!        real, dimension(3,2*nxvh,nyv), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine cibpois23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
!        real, dimension(3,2*nxvh,nyv), intent(in) :: cu
         real, dimension(*), intent(in) :: cu
!        complex, dimension(3,nxvh,nyv), intent(inout) :: bxy
         complex, dimension(*), intent(inout) :: bxy
!        complex, dimension(nxhd,nyhd), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cmaxwel2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv, &
     &nxhd,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
!        complex, dimension(3,nxvh,nyv), intent(inout) :: exy, bxy
         complex, dimension(*), intent(inout) :: exy, bxy
!        real, dimension(3,2*nxvh,nyv), intent(in) :: cu
         real, dimension(*), intent(in) :: cu
!        complex, dimension(nxhd,nyhd), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cemfield2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd&
     &)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
!        real, dimension(3,2*nxvh,nyv), intent(inout) :: fxy
         real, dimension(*), intent(inout) :: fxy
!        complex, dimension(3,nxvh,nyv), intent(in) :: exy
         complex, dimension(*), intent(in) :: exy
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
         integer, intent(in) :: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
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
         subroutine cwfft2r3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
!        real, dimension(3,2*nxhd,nyd), intent(inout) :: f
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
