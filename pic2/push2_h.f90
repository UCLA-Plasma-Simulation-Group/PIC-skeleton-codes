!-----------------------------------------------------------------------
! Interface file for push2.f
      module push2_h
      implicit none
!
      interface
         subroutine DISTR2(part,vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx,ny,&
     &ipbc)
         implicit none
         integer, intent(in) :: npx, npy, idimp, nop, nx, ny, ipbc
         real, intent(in) :: vtx, vty, vdx, vdy
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine GPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv, &
     &ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(2,nxv,nyv), intent(in) :: fxy
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
         subroutine CGUARD2L(fxy,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye), intent(inout) :: fxy
         end subroutine
      end interface
!
      interface
         subroutine AGUARD2L(q,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
         real, dimension(nxe,nye), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine POIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,&
     &nxhd,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         real, dimension(2*nxvh,nyv), intent(in) :: q
         real, dimension(2,2*nxvh,nyv), intent(inout) :: fxy
         complex, dimension(nxhd,nyhd), intent(inout) :: ffc
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
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         real, dimension(2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2R2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd, &
     &nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2,2*nxhd,nyd), intent(inout) :: f
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
         subroutine FFT2R2X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd&
     &,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2R2Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd&
     &,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2,2*nxhd,nyd), intent(inout) :: f
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
