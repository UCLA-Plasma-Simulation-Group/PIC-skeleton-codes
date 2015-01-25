!-----------------------------------------------------------------------
! Fortran interface to CUDA Library for Skeleton 2-1/2D Electromagnetic
! GPU-MPI PIC Code */
! written by Viktor K. Decyk, UCLA
      module gpupbpush2_h
      implicit none
!
      interface
         subroutine cgpuppgbppush23l(gp_ppart,gp_fxy,gp_bxy,gp_kpic,noff&
     &,nyp,qbm,dt,dtc,gp_ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,&
     &ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nypmx, mx1, mxyp1
         integer :: ipbc
         integer :: noff, nyp
         real :: qbm, dt, dtc
         integer, dimension(2) :: gp_ppart, gp_fxy, gp_bxy, gp_kpic
         integer, dimension(2) :: gp_ek
         end subroutine
      end interface
!
      interface
         subroutine cgpuppgrbppush23l(gp_ppart,gp_fxy,gp_bxy,gp_kpic,   &
     &noff,nyp,qbm,dt,dtc,ci,gp_ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1&
     &,mxyp1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nypmx, mx1, mxyp1
         integer :: ipbc
         integer :: noff, nyp
         real :: qbm, dt, dtc, ci
         integer, dimension(2) :: gp_ppart, gp_fxy, gp_bxy, gp_kpic
         integer, dimension(2) :: gp_ek
         end subroutine
      end interface
!
      interface
         subroutine cgpu2ppgppost2l(gp_ppart,gp_q,gp_kpic,noff,qm,idimp,&
     &nppmx,mx,my,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer :: idimp, nppmx, mx, my, nxv, nypmx, mx1, mxyp1
         integer :: noff
         real :: qm
         integer, dimension(2) :: gp_ppart, gp_q, gp_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpu2ppjppost2l(gp_ppart,gp_cu,gp_kpic,noff,qm,dt,  &
     &nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1, mxyp1
         integer :: ipbc
         integer :: noff
         real :: qm, dt
         integer, dimension(2) :: gp_ppart, gp_cu, gp_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpu2pprjppost2l(gp_ppart,gp_cu,gp_kpic,noff,qm,dt, &
     &ci,nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer :: nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1, mxyp1
         integer :: ipbc
         integer :: noff
         real :: qm, dt, ci
         integer, dimension(2) :: gp_ppart, gp_cu, gp_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcaguard2xl(gp_qc,gp_scs,gp_q,nyp,nx,nxe,nypmx,&
     &nxvh,kypd)
         implicit none
         integer :: nx, nxe, nypmx, nxvh, kypd
         integer :: nyp
         integer, dimension(2) :: gp_qc, gp_scs, gp_q
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcaguard2yl(gp_qc,gp_scr,nx,nxvh,kypd)
         implicit none
         integer :: nx, nxvh, kypd
         integer, dimension(2) :: gp_qc, gp_scr
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcacguard2xl(gp_cuc,gp_scs,gp_cu,nyp,nx,nxe,   &
     &nypmx,nxvh,kypd)
         implicit none
         integer :: nx, nxe, nypmx, nxvh, kypd
         integer :: nyp
         integer, dimension(2) :: gp_cuc, gp_scs, gp_cu
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcacguard2yl(gp_fvc,gp_scr,nx,nxvh,kypd)
         implicit none
         integer :: nx, nxvh, kypd
         integer, dimension(2) :: gp_fvc, gp_scr
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcbguard2xl(gp_fxyc,gp_scs,gp_fxy,nyp,nx,nxe,  &
     &nypmx,nxvh,kypd)
         implicit none
         integer :: nx, nxe, nypmx, nxvh, kypd
         integer :: nyp
         integer, dimension(2) :: gp_fxyc, gp_scs, gp_fxy
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcbguard2yl(gp_fxy,gp_scr,nyp,nx,nxe,nxvh,nypmx&
     &)
         implicit none
         integer :: nx, nxe, nxvh, nypmx
         integer :: nyp
         integer, dimension(2) :: gp_fxy, gp_scr
         end subroutine
      end interface
!
      interface
         subroutine cgpupppord2la(gp_ppart,gp_ppbuff,gp_sbufl,gp_sbufr, &
     &gp_kpic,gp_ncl,gp_ihole,gp_ncll,gp_nclr,noff,nyp,idimp,nppmx,nx,ny&
     &,mx,my,mx1,myp1,npbmx,ntmax,nbmax,gp_irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, mx1, myp1, npbmx
         integer :: ntmax, nbmax
         integer :: noff, nyp
         integer, dimension(2) :: gp_ppart, gp_ppbuff
         integer, dimension(2) :: gp_sbufl, gp_sbufr
         integer, dimension(2) :: gp_kpic, gp_ncl, gp_ihole
         integer, dimension(2) :: gp_ncll, gp_nclr, gp_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpupppord2lb(gp_ppart,gp_ppbuff,gp_rbufl,gp_rbufr, &
     &gp_kpic,gp_ncl,gp_ihole,gp_mcll,gp_mclr,idimp,nppmx,mx1,myp1,npbmx&
     &,ntmax,nbmax,gp_irc)
         implicit none
         integer :: idimp, nppmx, mx1, myp1, npbmx, ntmax, nbmax
         integer, dimension(2) :: gp_ppart, gp_ppbuff
         integer, dimension(2) :: gp_rbufl, gp_rbufr
         integer, dimension(2) :: gp_kpic, gp_ncl, gp_ihole
         integer, dimension(2) :: gp_mcll, gp_mclr, gp_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpuppois23t(gp_qt,gp_fxyt,gp_ffct,gp_we,nx,ny,kstrt&
     &,nyv,kxp1,nyhd)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp1, nyhd
         integer, dimension(2) :: gp_qt, gp_fxyt, gp_ffct, gp_we
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcuperp2t(gp_cut,nx,ny,kstrt,nyv,kxp1)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp1, nyhd
         integer, dimension(2) :: gp_cut
         end subroutine
      end interface
!
      interface
         subroutine cgpuippbpoisp23t(gp_cut,gp_bxyt,gp_ffct,ci,gp_wm,nx,&
     &ny,kstrt,nyv,kxp1,nyhd)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp1, nyhd
         real :: ci
         integer, dimension(2) :: gp_cut, gp_bxyt, gp_ffct, gp_wm
         end subroutine
      end interface
!
      interface
         subroutine cgpuppmaxwel2t(gp_exyt,gp_bxyt,gp_cut,gp_ffct,affp, &
     &ci,dt,gp_wf,gp_wm,nx,ny,kstrt,nyv,kxp1,nyhd)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp1, nyhd
         real :: affp, ci, dt
         integer, dimension(2) :: gp_exyt, gp_bxyt, gp_cut, gp_ffct
         integer, dimension(2) :: gp_wf, gp_wm
         end subroutine
      end interface
!
      interface
         subroutine cgpuppemfield2t(gp_fxyt,gp_exyt,gp_ffct,isign,nx,ny,&
     &kstrt,nyv,kxp1,nyhd)
         implicit none
         integer :: isign, nx, ny, kstrt, nyv, kxp1, nyhd
         integer, dimension(2) :: gp_fxyt, gp_exyt, gp_ffct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwppfft2rcsx(gp_f,gp_bsm,isign,gp_mixup,gp_sct,  &
     &indx,indy,kstrt,nvp,kxp1,kyp,nxhd,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, nvp, kxp1, kyp
         integer :: nxhd, kypd, nxhyd, nxyhd
         integer, dimension(2) :: gp_f, gp_bsm, gp_mixup, gp_sct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwppfft2rcsy(gp_g,gp_brm,isign,gp_mixup,gp_sct,  &
     &indx,indy,kstrt,nvp,kxp1,kyp,nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, nvp, kxp1, kyp
         integer :: nyd, nxhyd, nxyhd
         integer, dimension(2) :: gp_g, gp_brm, gp_mixup, gp_sct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwppfft2rcsxn(gp_fn,gp_bsm,isign,gp_mixup,gp_sct,&
     &indx,indy,ndim,kstrt,nvp,kxp1,kyp,nxhd,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, ndim, kstrt, nvp, kxp1, kyp
         integer :: nxhd, kypd, nxhyd, nxyhd
         integer, dimension(2) :: gp_fn, gp_bsm, gp_mixup, gp_sct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwppfft2rcsyn(gp_gn,gp_brm,isign,gp_mixup,gp_sct,&
     &indx,indy,ndim,kstrt,nvp,kxp1,kyp,nyd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, ndim, kstrt, nvp, kxp1, kyp
         integer :: nyd, nxhyd, nxyhd
         integer, dimension(2) :: gp_gn, gp_brm, gp_mixup, gp_sct
         end subroutine
      end interface
!
      interface
         subroutine cgpuppltpose(gp_f,gp_g,nx,ny,kxp,kyp,kstrt,nxv,nyv)
         implicit none
         integer :: nx, ny, kxp, kyp, kstrt, nxv, nyv
         integer, dimension(2) :: gp_f, gp_g
         end subroutine
      end interface
!
      interface
         subroutine cgpuppltposen(gp_fn,gp_gn,nx,ny,kxp,kyp,kstrt,ndim, &
     &nxv,nyv)
         implicit none
         integer :: nx, ny, kxp, kyp, kstrt, ndim, nxv, nyv
         integer, dimension(2) :: gp_fn, gp_gn
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


