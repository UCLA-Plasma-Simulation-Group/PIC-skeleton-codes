!-----------------------------------------------------------------------
! Interface file for kncbpush3.c
      module kncbpush3_h
      implicit none
!
      interface
         subroutine ckncgbpush3lt(part,fxyz,bxyz,qbm,dt,dtc,ek,idimp,nop&
     &,npe,nx,ny,nz,nxv,nyv,nzv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, npe, nx, ny, nz
         integer, intent(in) :: nxv, nyv, nzv, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(npe,idimp), intent(inout) :: part
         real, dimension(4,nxv*nyv*nzv), intent(in) :: fxyz, bxyz
         end subroutine
      end interface
!
      interface
         subroutine ckncgrbpush3lt(part,fxyz,bxyz,qbm,dt,dtc,ci,ek,idimp&
     &,nop,npe,nx,ny,nz,nxv,nyv,nzv,ipbc)
         implicit none
         integer, intent(in) :: idimp, nop, npe, nx, ny, nz
         integer, intent(in) :: nxv, nyv, nzv, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(npe,idimp), intent(inout) :: part
         real, dimension(4,nxv*nyv*nzv), intent(in) :: fxyz, bxyz
         end subroutine
      end interface
!
      interface
         subroutine ckncgpost3lt(part,q,qm,nop,npe,idimp,nxv,nyv,nzv)
         implicit none
         integer, intent(in) :: nop, npe, idimp, nxv, nyv, nzv
         real, intent(in) :: qm
         real, dimension(npe,idimp), intent(in) :: part
         real, dimension(nxv*nyv*nzv), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine cknc2gpost3lt(part,q,qm,nop,npe,idimp,nxv,nyv,nzv)
         implicit none
         integer, intent(in) :: nop, npe, idimp, nxv, nyv, nzv
         real, intent(in) :: qm
         real, dimension(npe,idimp), intent(in) :: part
         real, dimension(nxv*nyv*nzv), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine ckncgjpost3lt(part,cu,qm,dt,nop,npe,idimp,nx,ny,nz, &
     &nxv,nyv,nzv,ipbc)
         implicit none
         integer, intent(in) :: nop, npe, idimp, nx, ny, nz
         integer, intent(in) :: nxv, nyv, nzv, ipbc
         real, intent(in) :: qm, dt
         real, dimension(npe,idimp), intent(inout) :: part
         real, dimension(4,nxv*nyv*nzv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine ckncgrjpost3lt(part,cu,qm,dt,ci,nop,npe,idimp,nx,ny,&
     &nz,nxv,nyv,nzv,ipbc)
         implicit none
         integer, intent(in) :: nop, npe, idimp, nx, ny, nz
         integer, intent(in) :: nxv, nyv, nzv, ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(npe,idimp), intent(inout) :: part
         real, dimension(4,nxv*nyv*nzv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine ckncdsortp3yzlt(parta,partb,npic,idimp,nop,npe,ny1, &
     &nyz1)
         implicit none
         integer, intent(in) :: idimp, nop, npe, ny1, nyz1
         real, dimension(npe,idimp), intent(in) :: parta
         real, dimension(npe,idimp), intent(inout) :: partb
         integer, dimension(nyz1), intent(inout) :: npic
         end subroutine
      end interface
!
      interface
         subroutine cknccguard3l(fxyz,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, dimension(4,nxe,nye,nze), intent(inout) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine ckncacguard3l(cu,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, dimension(4,nxe,nye,nye), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine ckncaguard3l(q,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
         real, dimension(nxe,nye,nze), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine ckncpois33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny, &
     &nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd)
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
         subroutine cknccuperp3(cu,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         real, dimension(4,2*nxvh,nyv,nzv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine ckncibpois33(cu,bxyz,ffc,ci,wm,nx,ny,nz,nxvh,nyv,nzv&
     &,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(4,2*nxvh,nyv,nzv), intent(in) :: cu
         complex, dimension(4,nxvh,nyv,nzv), intent(inout) :: bxyz
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ckncmaxwel3(exyz,bxyz,cu,ffc,ci,dt,wf,wm,nx,ny,nz,  &
     &nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
         complex, dimension(4,nxvh,nyv,nzv), intent(inout) :: exyz, bxyz
         real, dimension(4,2*nxvh,nyv,nzv), intent(in) :: cu
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ckncemfield3(fxyz,exyz,ffc,isign,nx,ny,nz,nxvh,nyv, &
     &nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, dimension(4,2*nxvh,nyv,nzv), intent(inout) :: fxyz
         complex, dimension(4,nxvh,nyv,nzv), intent(in) :: exyz
         complex, dimension(nxhd,nyhd,nzhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine ckncwfft3rvx(f,isign,mixup,sct,indx,indy,indz,nxhd, &
     &nyd,nzd,nxhyzd,nxyzhd)
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
         subroutine ckncwfft3rv3(f,isign,mixup,sct,indx,indy,indz,nxhd, &
     &nyd,nzd,nxhyzd,nxyzhd)
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
          subroutine ckncfft3rvxy(f,isign,mixup,sct,indx,indy,indz,nzi, &
     &nzp,nxhd,nyd,nzd,nxhyzd,nxyzhd)
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
         subroutine ckncfft3rxz(f,isign,mixup,sct,indx,indy,indz,nyi,nyp&
     &,nxhd,nyd,nzd,nxhyzd,nxyzhd)
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
         subroutine ckncfft3rv3xy(f,isign,mixup,sct,indx,indy,indz,nzi, &
     &nzp,nxhd,nyd,nzd,nxhyzd,nxyzhd)
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
         subroutine ckncfft3rv3z(f,isign,mixup,sct,indx,indy,indz,nyi,  &
     &nyp,nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer , intent(in):: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd), intent(inout) :: f
         integer, dimension(nxhyzd), intent(in) :: mixup
         complex, dimension(nxyzhd), intent(in) :: sct
         end subroutine
      end interface
!
      end module
