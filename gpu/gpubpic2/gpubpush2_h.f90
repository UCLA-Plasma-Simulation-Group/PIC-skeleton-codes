!-----------------------------------------------------------------------
! Fortran interface to CUDA Library for Skeleton 2D Electrostatic GPU
! PIC Code */
! written by Viktor K. Decyk, UCLA
      module gpubpush2_h
      implicit none
!
      interface
         subroutine cgpubppush23l(gp_ppart,gp_fxy,gp_bxy,gp_kpic,qbm,dt,&
     &dtc,gp_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qbm, dt, dtc
         integer, dimension(2) :: gp_ppart, gp_fxy, gp_bxy, gp_kpic
         integer, dimension(2) :: gp_ek
         end subroutine
      end interface
!
      interface
         subroutine cgpubppushf23l(gp_ppart,gp_fxy,gp_bxy,gp_kpic,      &
     &gp_ncl,gp_ihole,qbm,dt,dtc,gp_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv, &
     &mx1,mxy1,ntmax,gp_irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax
         real :: qbm, dt, dtc
         integer, dimension(2) :: gp_ppart, gp_fxy, gp_bxy, gp_kpic
         integer, dimension(2) :: gp_ncl, gp_ihole, gp_ek, gp_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpurbppush23l(gp_ppart,gp_fxy,gp_bxy,gp_kpic,qbm,dt&
     &,dtc,ci,gp_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qbm, dt, dtc, ci
         integer, dimension(2) :: gp_ppart, gp_fxy, gp_bxy, gp_kpic
         integer, dimension(2) :: gp_ek
         end subroutine
      end interface
!
      interface
         subroutine cgpurbppushf23l(gp_ppart,gp_fxy,gp_bxy,gp_kpic,     &
     &gp_ncl,gp_ihole,qbm,dt,dtc,ci,gp_ek,idimp,nppmx,nx,ny,mx,my,nxv,  &
     &nyv,mx1,mxy1,ntmax,gp_irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax
         real :: qbm, dt, dtc, ci
         integer, dimension(2) :: gp_ppart, gp_fxy, gp_bxy, gp_kpic
         integer, dimension(2) :: gp_ncl, gp_ihole, gp_ek, gp_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpu2ppost2l(gp_ppart,gp_q,gp_kpic,qm,nppmx,idimp,  &
     &mx,my,nxv,nyv,mx1,mxy1) 
         implicit none
         integer :: nppmx, idimp, mx, my, nxv, nyv, mx1, mxy1
         real :: qm
         integer, dimension(2) :: gp_ppart, gp_q, gp_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpu2jppost2l(gp_ppart,gp_cu,gp_kpic,qm,dt,nppmx,   &
     &idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qm, dt
         integer, dimension(2) :: gp_ppart, gp_cu, gp_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpu2jppostf2l(gp_ppart,gp_cu,gp_kpic,gp_ncl,       &
     &gp_ihole,qm,dt,nppmx,idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,    &
     &gp_irc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax
         real :: qm, dt
         integer, dimension(2) :: gp_ppart, gp_cu, gp_kpic, gp_ncl
         integer, dimension(2) :: gp_ihole, gp_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpu2rjppost2l(gp_ppart,gp_cu,gp_kpic,qm,dt,ci,nppmx&
     &,idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qm, dt, ci
         integer, dimension(2) :: gp_ppart, gp_cu, gp_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpu2rjppostf2l(gp_ppart,gp_cu,gp_kpic,gp_ncl,      &
     &gp_ihole,qm,dt,ci,nppmx,idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax, &
     &gp_irc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax
         real :: qm, dt, ci
         integer, dimension(2) :: gp_ppart, gp_cu, gp_kpic, gp_ncl
         integer, dimension(2) :: gp_ihole, gp_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpucaguard2l(gp_qc,gp_q,nx,ny,nxe,nye,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxe, nye, nxvh, nyv
         integer, dimension(2) :: gp_qc, gp_q
         end subroutine
      end interface
!
      interface
         subroutine cgpucacguard2l(gp_cuc,gp_cu,nx,ny,nxe,nye,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxe, nye, nxvh, nyv
         integer, dimension(2) :: gp_cuc, gp_cu
         end subroutine
      end interface
!
      interface
         subroutine cgpucbguard2l(gp_bxyc,gp_bxy,nx,ny,nxe,nye,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxe, nye, nxvh, nyv
         integer, dimension(2) :: gp_bxyc, gp_bxy
         end subroutine
      end interface
!
      interface
         subroutine cgpuppord2l(gp_ppart,gp_ppbuff,gp_kpic,gp_ncl,      &
     &gp_ihole,idimp,nppmx,nx,ny,mx,my,mx1,my1,npbmx,ntmax,gp_irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, mx1, my1, npbmx, ntmax
         integer, dimension(2) :: gp_ppart, gp_ppbuff, gp_kpic, gp_ncl
         integer, dimension(2) :: gp_ihole, gp_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpuppordf2l(gp_ppart,gp_ppbuff,gp_kpic,gp_ncl,     &
     &gp_ihole,idimp,nppmx,mx1,my1,npbmx,ntmax,gp_irc)
         implicit none
         integer :: idimp, nppmx, mx1, my1, npbmx, ntmax
         integer, dimension(2) :: gp_ppart, gp_ppbuff, gp_kpic, gp_ncl
         integer, dimension(2) :: gp_ihole, gp_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpupois23t(gp_qt,gp_fxyt,gp_ffct,gp_we,nx,ny,nxvh, &
     &nyv,nxhd,nyhd)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         integer, dimension(2) :: gp_qt, gp_fxyt, gp_ffct, gp_we
         end subroutine
      end interface
!
      interface
         subroutine cgpucuperp2t(gp_cut,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
         integer, dimension(2) :: gp_cut
         end subroutine
      end interface
!
      interface
         subroutine cgpuibpois23t(gp_cut,gp_bxyt,gp_ffct,ci,gp_wm,nx,ny,&
     &nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ci
         integer, dimension(2) :: gp_cut, gp_bxyt, gp_ffct, gp_wm
         end subroutine
      end interface
!
      interface
         subroutine cgpumaxwel2t(gp_exyt,gp_bxyt,gp_cut,gp_ffct,ci,dt,  &
     &gp_wf,gp_wm,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ci, dt
         integer, dimension(2) :: gp_exyt, gp_bxyt, gp_cut, gp_ffct
         integer, dimension(2) :: gp_wf, gp_wm
         end subroutine
      end interface
!
      interface
         subroutine cgpuemfield2t(gp_fxyt,gp_exyt,gp_ffct,isign,nx,ny,  &
     &nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         integer, dimension(2) :: gp_fxyt, gp_exyt, gp_ffct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwfft2rcs(gp_f,gp_g,isign,gp_mixup,gp_sct,indx,  &
     &indy,nxhd,nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         integer, dimension(2) :: gp_f, gp_g, gp_mixup, gp_sct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwfft2rcsn(gp_fn,gp_gn,isign,gp_mixup,gp_sct,indx&
     &,indy,ndim,nxhd,nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, ndim, nxhd, nyd, nxhyd, nxyhd
         integer, dimension(2) :: gp_fn, gp_gn, gp_mixup, gp_sct
         end subroutine
      end interface
!
      interface
         subroutine cgpusum2(gp_a,gp_sa,nx)
         implicit none
         integer :: nx
         integer, dimension(2) :: gp_a, gp_sa
         end subroutine
      end interface
!
      end module


