!-----------------------------------------------------------------------
! Interface file for bpush3.f
      module bpush3_h
      implicit none
!
      interface
         subroutine DISTR3(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,    &
     &idimp,nop,nx,ny,nz,ipbc)
         implicit none
         integer, intent(in) :: npx, npy, npz, idimp, nop, nx, ny, nz
         integer, intent(in) :: ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine GBPUSH3L(part,fxyz,bxyz,qbm,dt,dtc,ek,idimp,nop,nx, &
     &ny,nz,nxv,nyv,nzv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nz, nxv, nyv, nzv
         integer, intent(in) :: ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         end subroutine
      end interface
!
      interface
         subroutine GRBPUSH3L(part,fxyz,bxyz,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nz,nxv,nyv,nzv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, nx, ny, nz, nxv, nyv, nzv
         integer, intent(in) :: ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         end subroutine
      end interface
!
      interface
         subroutine GPOST3L(part,q,qm,nop,idimp,nxv,nyv,nzv)
         implicit none
         integer, intent(in) :: nop, idimp, nxv, nyv, nzv
         real, intent(in) :: qm
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(nxv,nyv,nzv), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine GJPOST3L(part,cu,qm,dt,nop,idimp,nx,ny,nz,nxv,nyv,  &
     &nzv,ipbc)
         implicit none
         integer, intent(in) :: nop, idimp, nx, ny, nz, nxv, nyv, nzv
         integer, intent(in) :: ipbc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine GRJPOST3L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nz,nxv,  &
     &nyv,nzv,ipbc)
         implicit none
         integer, intent(in) :: nop, idimp, nx, ny, nz, nxv, nyv, nzv
         integer, intent(in) :: ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine DSORTP3YZL(parta,partb,npic,idimp,nop,ny1,nyz1)
         implicit none
         integer, intent(in) :: idimp, nop, ny1, nyz1
         real, dimension(idimp,nop), intent(in) :: parta
         real, dimension(idimp,nop), intent(inout) :: partb
         integer, dimension(nyz1), intent(inout) :: npic
         end subroutine
      end interface
!
      interface
         subroutine CGUARD3L(fxyz,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, dimension(3,nxe,nye,nze), intent(inout) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine ACGUARD3L(cu,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, dimension(3,nxe,nye,nye), intent(inout) :: cu
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
         subroutine POIS33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,  &
     &nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ax, ay, az, affp
         real, intent(inout) :: we
         real, dimension(2*nxvh,nyv,nzv), intent(in) :: q
         real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: fxyz
         complex, dimension(nxhd,nyhd,nzhd), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine CUPERP3(cu,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine IBPOIS33(cu,bxyz,ffc,ci,wm,nx,ny,nz,nxvh,nyv,nzv,   &
     &nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(3,2*nxvh,nyv,nzv), intent(in) :: cu
         complex, dimension(3,nxvh,nyv,nzv), intent(inout) :: bxyz
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MAXWEL3(exyz,bxyz,cu,ffc,ci,dt,wf,wm,nx,ny,nz,nxvh, &
     &nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(3,nxvh,nyv,nzv), intent(inout) :: exyz, bxyz
         real, dimension(3,2*nxvh,nyv,nzv), intent(in) :: cu
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine EMFIELD3(fxyz,exyz,ffc,isign,nx,ny,nz,nxvh,nyv,nzv, &
     &nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: fxyz
         complex, dimension(3,nxvh,nyv,nzv), intent(in) :: exyz
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
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
         subroutine WFFT3RX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,  &
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
         subroutine WFFT3R3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,  &
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer, intent(in) :: nxhyzd, nxyzhd
         real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
          subroutine FFT3RXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,  &
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
         subroutine FFT3R3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,  &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3R3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,   &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
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
