!-----------------------------------------------------------------------
! Interface file for bpush2.f
      module bpush2_h
      implicit none
!
      interface
         subroutine DISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp, &
     &nop,nx,ny,ipbc)
         implicit none
         integer, intent(in) :: npx, npy, idimp, nop, nx, ny, ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine GBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine GRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop, &
     &nx,ny,nxv,nyv,ipbc)
         implicit none
         integer , intent(in):: idimp, nop, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         end subroutine
      end interface
!
      interface
         subroutine GPOST2L(part,q,qm,nop,idimp,nxv,nyv)
         implicit none
         integer, intent(in) :: nop, idimp, nxv, nyv
         real, intent(in) :: qm
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(nxv,nyv), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine GJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer, intent(in) :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv,nyv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine GRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv, &
     &ipbc)
         implicit none
         integer, intent(in) :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv,nyv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine DSORTP2YL(parta,partb,npic,idimp,nop,ny1)
         implicit none
         integer, intent(in) :: idimp, nop, ny1
         real, dimension(idimp,nop), intent(in) :: parta
         real, dimension(idimp,nop), intent(inout) :: partb
         integer, dimension(ny1), intent(inout) :: npic
         end subroutine
      end interface
!
      interface
         subroutine BGUARD2L(bxy,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye), intent(inout)  :: bxy
         end subroutine
      end interface
!
      interface
         subroutine ACGUARD2L(cu,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye), intent(inout)  :: cu
         end subroutine
      end interface
!
      interface
         subroutine AGUARD2L(q,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
         real, dimension(nxe,nye), intent(inout)  :: q
         end subroutine
      end interface
!
      interface
         subroutine POIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,&
     &nxhd,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         real, dimension(2*nxvh,nyv), intent(in) :: q
         real, dimension(3,2*nxvh,nyv), intent(inout) :: fxy
         complex, dimension(nxhd,nyhd), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine CUPERP2(cu,nx,ny,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv
         real, dimension(3,2*nxvh,nyv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine IBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(3,2*nxvh,nyv), intent(in) :: cu
         complex, dimension(3,nxvh,nyv), intent(inout) :: bxy
         complex, dimension(nxhd,nyhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MAXWEL2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,  &
     &nxhd,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(3,nxvh,nyv), intent(inout) :: exy, bxy
         real, dimension(3,2*nxvh,nyv), intent(in) :: cu
         complex, dimension(nxhd,nyhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine EMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, dimension(3,2*nxvh,nyv), intent(inout) :: fxy
         complex, dimension(3,nxvh,nyv), intent(in) :: exy
         complex, dimension(nxhd,nyhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: indx, indy, nxhyd, nxyhd
         integer, dimension(nxhyd), intent(inout) :: mixup
         complex, dimension(nxyhd), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd, &
     &nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2R3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd, &
     &nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd&
     &,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd&
     &,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd&
     &,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd&
     &,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
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
