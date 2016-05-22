!-----------------------------------------------------------------------
! Interface file for vpush3.f
      module vpush3_h
      implicit none
!
      interface
         subroutine DISTR3T(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,   &
     &idimp,npe,nx,ny,nz,ipbc)
         implicit none
         integer, intent(in) :: npx, npy, npz, idimp, npe, nx, ny, nz
         integer, intent(in) :: ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(npe,idimp), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine GPUSH3LT(part,fxyz,qbm,dt,ek,idimp,nop,npe,nx,ny,nz,&
     &nxv,nyv,nzv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, npe, nx, ny, nz
         integer, intent(in) :: nxv, nyv, nzv, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
         real, dimension(npe,idimp), intent(inout) :: part
         real, dimension(4,nxv*nyv*nzv), intent(in) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine VGPUSH3LT(part,fxyz,qbm,dt,ek,idimp,nop,npe,nx,ny,nz&
     &,nxv,nyv,nzv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, npe, nx, ny, nz
         integer, intent(in) :: nxv, nyv, nzv, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
         real, dimension(npe,idimp), intent(inout) :: part
         real, dimension(4,nxv*nyv*nzv), intent(in) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine GPOST3LT(part,q,qm,nop,npe,idimp,nxv,nyv,nzv)
         implicit none
         integer, intent(in) :: nop, npe, idimp, nxv, nyv, nzv
         real, intent(in) :: qm
         real, dimension(npe,idimp), intent(in) :: part
         real, dimension(nxv*nyv*nzv), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine VGPOST3LT(part,q,qm,nop,npe,idimp,nxv,nyv,nzv)
         implicit none
         integer, intent(in) :: nop, npe, idimp, nxv, nyv, nzv
         real, intent(in) :: qm
         real, dimension(npe,idimp), intent(in) :: part
         real, dimension(nxv*nyv*nzv), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine DSORTP3YZLT(parta,partb,npic,idimp,nop,npe,ny1,nyz1)
         implicit none
         integer, intent(in) :: idimp, nop, npe, ny1, nyz1
         real, dimension(npe,idimp), intent(in) :: parta
         real, dimension(npe,idimp), intent(inout) :: partb
         integer, dimension(nyz1), intent(inout) :: npic
         end subroutine
      end interface
!
      interface
         subroutine CGUARD3L(fxyz,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, dimension(4,nxe,nye,nze), intent(inout) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine AGUARD3L(q,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, dimension(nxe,nye,nze), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine VPOIS33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz, &
     &nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ax, ay, az, affp
         real, intent(inout) :: we
         real, dimension(2*nxvh,nyv,nzv), intent(in) :: q
         real, dimension(4,2*nxvh,nyv,nzv), intent(inout) :: fxyz
         complex, dimension(nxhd,nyhd,nzhd), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: indx, indy, indz, nxhyzd, nxyzhd
         integer, dimension(nxhyzd), intent(inout) :: mixup
         complex, dimension(nxyzhd), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RVX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd, &
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer, intent(in) :: nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RV3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd, &
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer, intent(in) :: nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
          subroutine FFT3RVXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RXZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,   &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RV3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RV3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,  &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer , intent(in):: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
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
