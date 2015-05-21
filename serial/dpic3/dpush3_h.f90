!-----------------------------------------------------------------------
! Interface file for dpush3.f
      module dpush3_h
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
         subroutine GMJPOST3L(part,amu,qm,nop,idimp,nxv,nyv,nzv)
         implicit none
         integer, intent(in) :: nop, idimp, nxv, nyv, nzv
         real, intent(in) :: qm
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(6,nxv,nyv,nzv), intent(inout) :: amu
         end subroutine
      end interface
!
      interface
         subroutine GDJPOST3L(part,fxyz,bxyz,dcu,amu,qm,qbm,dt,idimp,nop&
     &,nxv,nyv,nzv)
         implicit none
         integer, intent(in) :: idimp, nop, nxv, nyv, nzv
         real, intent(in) :: qm, qbm, dt
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: dcu
         real, dimension(6,nxv,nyv,nzv), intent(inout) :: amu
         end subroutine
      end interface
!
      interface
         subroutine GDCJPOST3L(part,fxyz,bxyz,cu,dcu,amu,qm,qbm,dt,idimp&
     &,nop,nxv,nyv,nzv)
         implicit none
         integer, intent(in) :: idimp, nop, nxv, nyv, nzv
         real, intent(in) :: qm, qbm, dt
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu, dcu
         real, dimension(6,nxv,nyv,nzv), intent(inout) :: amu
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
         subroutine AMCGUARD3L(amu,nx,ny,nz,nxe,nye,nze,ndim)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze, ndim
         real, dimension(ndim,nxe,nye,nze), intent(inout) :: amu
         end subroutine
      end interface
!
      interface
         subroutine ASCFGUARD3L(dcu,cus,q2m0,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, intent(in) :: q2m0
         real, dimension(3,nxe,nye,nze), intent(inout) :: dcu
         real, dimension(3,nxe,nye,nze), intent(in) :: cus
         end subroutine
      end interface
!
      interface
         subroutine FWPMINMX3(qe,qbme,wpmax,wpmin,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, intent(in) :: qbme
         real, intent(inout) :: wpmax, wpmin
         real, dimension(nxe,nye,nze), intent(in) :: qe
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
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
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
         subroutine BBPOIS33(cu,bxyz,ffc,ci,wm,nx,ny,nz,nxvh,nyv,nzv,   &
     &nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(3,2*nxvh,nyv,nzv), intent(in) :: cu
         real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: bxyz
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine BADDEXT3(bxyz,omx,omy,omz,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, intent(in) :: omx, omy, omz
         real, dimension(3,nxe,nye,nze), intent(inout) :: bxyz
         end subroutine
      end interface
!
      interface
         subroutine DCUPERP3(dcu,amu,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: dcu
         real, dimension(6,2*nxvh,nyv,nzv), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine ADCUPERP3(dcu,amu,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: dcu
         real, dimension(6,2*nxvh,nyv,nzv), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine EPOIS33(dcu,exyz,isign,ffe,ax,ay,az,affp,wp0,ci,wf, &
     &nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ax, ay, az, affp, wp0, ci
         real, intent(inout) ::  wf
         real, dimension(3,2*nxvh,nyv,nzv), intent(in) :: dcu
         real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: exyz
         complex, dimension(nxhd,nyhd), intent(inout) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine ADDVRFIELD3(a,b,c,ndim,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nxe, nye, nze, ndim
         real, dimension(ndim,nxe,nye,nze), intent(inout) :: a
         real, dimension(ndim,nxe,nye,nze), intent(in) :: b, c
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
         subroutine WFFT3RN(f,ss,isign,mixup,sct,indx,indy,indz,nxhd,nyd&
     &,nzd,ndim,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer, intent(in) :: ndim, nxhyzd, nxyzhd
         real, dimension(ndim,2*nxhd,nyd,nzd), intent(inout) :: f
         complex, dimension(ndim,nxhd), intent(inout) :: ss
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
         subroutine FFT3RNXY(f,ss,isign,mixup,sct,indx,indy,indz,nzi,nzp&
     &,nxhd,nyd,nzd,ndim,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, ndim, nxhyzd, nxyzhd
         real, dimension(ndim,2*nxhd,nyd,nzd), intent(inout) :: f
         complex, dimension(ndim,nxhd), intent(inout) :: ss
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RNZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,   &
     &nxhd,nyd,nzd,ndim,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, ndim, nxhyzd, nxyzhd
         real, dimension(ndim,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine SWAP3CN(f,s,isign,nxh,ny,nzi,nzt,nxhd,nyd,nzd,ndim)
         implicit none
         integer, intent(in) :: isign, nxh, ny, nzi, nzt, nxhd, nyd, nzd, ndim
         real, dimension(ndim,2*nxhd,nyd,nzd), intent(inout) :: f
         complex, dimension(ndim*nxhd), intent(inout) :: s
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
