!-----------------------------------------------------------------------
! Interface file for push3.f
!     module push3_h
!     implicit none
!
      interface
         subroutine DISTR3(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,    &
     &idimp,nop,nx,ny,nz,ipbc)
         implicit none
         integer, intent(in) :: npx, npy, npz, idimp, nop, nx, ny, nz
         integer, intent(in) :: ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine GPUSH3L(part,fxyz,qbm,dt,ek,idimp,nop,nx,ny,nz,nxv, &
     &nyv,nzv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nz, nxv, nyv, nzv
         integer, intent(in) :: ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz
         real, dimension(*), intent(in) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine GPOST3L(part,q,qm,nop,idimp,nxv,nyv,nzv)
         implicit none
         integer, intent(in) :: nop, idimp, nxv, nyv, nzv
         real, intent(in) :: qm
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        real, dimension(nxv,nyv,nzv), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine DSORTP3YZL(parta,partb,npic,idimp,nop,ny1,nyz1)
         implicit none
         integer, intent(in) :: idimp, nop, ny1, nyz1
!        real, dimension(idimp,nop), intent(in) :: parta
         real, dimension(*), intent(in) :: parta
!        real, dimension(idimp,nop), intent(inout) :: partb
         real, dimension(*), intent(inout) :: partb
!        integer, dimension(nyz1), intent(inout) :: npic
         integer, dimension(*), intent(inout) :: npic
         end subroutine
      end interface
!
      interface
         subroutine CGUARD3L(fxyz,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
!        real, dimension(3,nxe,nye,nze), intent(inout) :: fxyz
         real, dimension(*), intent(inout) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine AGUARD3L(q,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
!        real, dimension(nxe,nye,nze), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine POIS33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,  &
     &nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ax, ay, az, affp
         real, intent(inout) :: we
!        real, dimension(2*nxvh,nyv,nzv), intent(in) :: q
         real, dimension(*), intent(in) :: q
!        real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: fxyz
         real, dimension(*), intent(inout) :: fxyz
!        complex, dimension(nxhd,nyhd,nzhd), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: indx, indy, indz, nxhyzd, nxyzhd
!        integer, dimension(nxhyzd), intent(inout) :: mixup
         integer, dimension(*), intent(inout) :: mixup
!        complex, dimension(nxyzhd), intent(inout) :: sct
         complex, dimension(*), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,  &
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer, intent(in) :: nxhyzd, nxyzhd
!        real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3R3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,  &
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer, intent(in) :: nxhyzd, nxyzhd
!        real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
          subroutine FFT3RXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,  &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
!        real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RXZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,   &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
!        real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3R3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,  &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
!        real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3R3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,   &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer , intent(in):: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
!        real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
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
