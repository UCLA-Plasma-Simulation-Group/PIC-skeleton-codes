!-----------------------------------------------------------------------
! Interface file for vmbpush3.f
      module vmbpush3_h
      implicit none
!
      interface
         subroutine DISTR3(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,    &
     &idimp,nop,nx,ny,nz,ipbc)
         implicit none
         integer :: npx, npy, npz, idimp, nop, nx, ny, nz, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
!
      interface
         subroutine DBLKP3L(part,kpic,nppmx,idimp,nop,mx,my,mz,mx1,my1, &
     &mxyz1,irc)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mz, mx1, my1, mxyz1, irc
         real, dimension(idimp,nop) :: part
         integer, dimension(mxyz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN3LT(part,ppart,kpic,nppmx,idimp,nop,mx,my,mz,&
     &mx1,my1,mxyz1,irc)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mz, mx1, my1, mxyz1, irc
         real, dimension(idimp,nop) :: part
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         integer, dimension(mxyz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN3LTP(part,ppart,kpic,kp,nppmx,idimp,nop,mx,my&
     &,mz,mx1,my1,mxyz1,irc)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mz, mx1, my1, mxyz1, irc
         real, dimension(idimp,nop) :: part
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         integer, dimension(mxyz1) :: kpic
         integer, dimension(nppmx,mxyz1) :: kp
         end subroutine
      end interface
!
      interface
         subroutine PPCHECK3LT(ppart,kpic,idimp,nppmx,nx,ny,nz,mx,my,mz,&
     &mx1,my1,mz1,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, my1, mz1
         integer :: irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         integer, dimension(mx1*my1*mz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSH3LT(ppart,fxyz,bxyz,kpic,qbm,dt,dtc,ek,idimp,&
     &nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSHF3LT(ppart,fxyz,bxyz,kpic,ncl,ihole,qbm,dt,  &
     &dtc,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,   &
     &ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GRBPPUSH3LT(ppart,fxyz,bxyz,kpic,qbm,dt,dtc,ci,ek,  &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRBPPUSHF3LT(ppart,fxyz,bxyz,kpic,ncl,ihole,qbm,dt, &
     &dtc,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,&
     &ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VGBPPUSH3LT(ppart,fxyz,bxyz,kpic,qbm,dt,dtc,ek,idimp&
     &,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGBPPUSHF3LT(ppart,fxyz,bxyz,kpic,ncl,ihole,qbm,dt, &
     &dtc,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,   &
     &ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VGRBPPUSH3LT(ppart,fxyz,bxyz,kpic,qbm,dt,dtc,ci,ek, &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGRBPPUSHF3LT(ppart,fxyz,bxyz,kpic,ncl,ihole,qbm,dt,&
     &dtc,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,&
     &ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine V2GBPPUSH3LT(ppart,fxyz,bxyz,kpic,qbm,dt,dtc,ek,    &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine V2GBPPUSHF3LT(ppart,fxyz,bxyz,kpic,ncl,ihole,qbm,dt,&
     &dtc,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,   &
     &ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine V2GRBPPUSH3LT(ppart,fxyz,bxyz,kpic,qbm,dt,dtc,ci,ek,&
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine V2GRBPPUSHF3LT(ppart,fxyz,bxyz,kpic,ncl,ihole,qbm,dt&
     &,dtc,ci,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1&
     &,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GPPOST3LT(ppart,q,kpic,qm,nppmx,idimp,mx,my,mz,nxv, &
     &nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer :: nppmx, idimp, mx, my, mz, nxv, nyv, nzv, mx1, my1
         integer :: mxyz1
         real :: qm
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         real, dimension(nxv,nyv,nzv) :: q
         integer, dimension(mxyz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGPPOST3LT(ppart,q,kpic,qm,nppmx,idimp,mx,my,mz,nxv,&
     &nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer :: nppmx, idimp, mx, my, mz, nxv, nyv, nzv, mx1, my1
         integer :: mxyz1
         real :: qm
         real, dimension(nppmx,idimp,mxyz1) :: ppart
         real, dimension(nxv,nyv,nzv) :: q
         integer, dimension(mxyz1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GJPPOST3LT(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny,nz,&
     &mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qm, dt
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GJPPOSTF3LT(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,    &
     &idimp,nx,ny,nz,mx, my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GRJPPOST3LT(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx,ny&
     &,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GRJPPOSTF3LT(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx,&
     &idimp,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VGJPPOST3LT(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny,nz&
     &,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qm, dt
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGJPPOSTF3LT(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,   &
     &idimp,nx,ny,nz,mx, my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt
         real, dimension(nppmx,idimp,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VGRJPPOST3LT(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx, &
     &ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine VGRJPPOSTF3LT(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx&
     &,idimp,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine V2GRJPPOST3LT(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx,&
     &ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qm, dt, ci
         real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPORDER3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &nx,ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, my1, mz1
         integer :: npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1*mz1) :: ppbuff
         integer, dimension(mx1*my1*mz1) :: kpic
         integer, dimension(26,mx1*my1*mz1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPORDERF3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,&
     &mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, my1, mz1, npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1*mz1) :: ppbuff
         integer, dimension(mx1*my1*mz1) :: kpic
         integer, dimension(26,mx1*my1*mz1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VPPORDER3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,&
     &nx,ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, my1, mz1
         integer :: npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1*mz1) :: ppbuff
         integer, dimension(mx1*my1*mz1) :: kpic
         integer, dimension(26,mx1*my1*mz1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine VPPORDERF3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx&
     &,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, my1, mz1, npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         real, dimension(npbmx,idimp,mx1*my1*mz1) :: ppbuff
         integer, dimension(mx1*my1*mz1) :: kpic
         integer, dimension(26,mx1*my1*mz1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine V2PPORDERF3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,    &
     &nppmx,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, my1, mz1, npbmx, ntmax, irc
         real, dimension(nppmx,idimp,mx1*my1*mz1) :: ppart
         real, dimension(idimp,npbmx,mx1*my1*mz1) :: ppbuff
         integer, dimension(mx1*my1*mz1) :: kpic
         integer, dimension(26,mx1*my1*mz1) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine CGUARD3L(fxyz,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer :: nx, ny, nz, nxe, nye, nze
         real, dimension(4,nxe,nye,nze) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine ACGUARD3L(cu,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer , intent(in):: nx, ny, nz, nxe, nye, nze
         real, dimension(3,nxe,nye,nye), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine AGUARD3L(q,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer :: nx, ny, nz, nxe, nye, nze
         real, dimension(nxe,nye,nze) :: q
         end subroutine
      end interface
!
      interface
         subroutine VMPOIS33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,&
     &nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer :: isign, nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
         real :: ax, ay, az, affp, we
         real, dimension(2*nxvh,nyv,nzv) :: q
         real, dimension(4,2*nxvh,nyv,nzv) :: fxyz
         complex, dimension(nxhd,nyhd,nzhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MCUPERP3(cu,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine VMIBPOIS33(cu,bxyz,ffc,ci,wm,nx,ny,nz,nxvh,nyv,nzv, &
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
         subroutine VMMAXWEL3(exyz,bxyz,cu,ffc,ci,dt,wf,wm,nx,ny,nz,nxvh&
     &,nyv,nzv,nxhd,nyhd,nzhd)
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
         subroutine VMEMFIELD3(fxyz,exyz,ffc,isign,nx,ny,nz,nxvh,nyv,nzv&
     &,nxhd,nyhd,nzhd)
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
         integer :: indx, indy, indz, nxhyzd, nxyzhd
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RVMX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,&
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RVM3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,&
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RVMXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RVMXZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RVM3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,&
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RVM3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer :: isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
         integer :: nxhyzd, nxyzhd
         real, dimension(4,2*nxhd,nyd,nzd) :: f
         integer, dimension(nxhyzd) :: mixup
         complex, dimension(nxyzhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine SET_SZERO3(q,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer :: mx, my, mz, nxv, nyv, nzv, mx1, my1, mxyz1
         real, dimension(nxv,nyv,nzv) :: q
         end subroutine
      end interface
!
      interface
         subroutine SET_VZERO3(cu,mx,my,mz,ndim,nxv,nyv,nzv,mx1,my1,    &
     &mxyz1)
         implicit none
         integer :: mx, my, mz, ndim, nxv, nyv, nzv, mx1, my1, mxyz1
         real, dimension(ndim,nxv,nyv,nzv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine SET_CVZERO3(exyz,nx,ny,nz,ndim,nxvh,nyv,nzv)
         implicit none
         integer :: nx, ny, nz, ndim, nxvh, nyv, nzv
         complex, dimension(ndim,nxvh,nyv,nzv) :: exyz
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
      end module
