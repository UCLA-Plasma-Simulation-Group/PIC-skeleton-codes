!-----------------------------------------------------------------------
! Fortran2003 interface to CUDA Library for Skeleton 2-1/2D
! Electromagnetic GPU-MPI PIC Code */
! written by Viktor K. Decyk, UCLA
      module gpupbpush2_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine cgpuppgbppush23l(g_ppart,g_fxy,g_bxy,g_kpic,noff,nyp&
     &,qbm,dt,dtc,g_ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)&
     & bind(C,name='cgpuppgbppush23l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv
         integer(c_int), value :: nypmx, mx1, mxyp1, ipbc
         integer(c_int), value :: noff, nyp
         real(c_float), value :: qbm, dt, dtc
         type (c_ptr), value :: g_ppart, g_fxy, g_bxy, g_kpic, g_ek
         end subroutine
      end interface
!
      interface
         subroutine cgpuppgrbppush23l(g_ppart,g_fxy, g_bxy,g_kpic,noff, &
     &nyp,qbm,dt,dtc,ci,g_ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1&
     &,ipbc) bind(C,name='cgpuppgrbppush23l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv
         integer(c_int), value :: nypmx, mx1, mxyp1, ipbc
         integer(c_int), value :: noff, nyp
         real(c_float), value :: qbm, dt, dtc, ci
         type (c_ptr), value :: g_ppart, g_fxy, g_bxy, g_kpic, g_ek
         end subroutine
      end interface
!
      interface
         subroutine cgpu2ppgppost2l(g_ppart,g_q,g_kpic,noff,qm,idimp,   &
     &nppmx,mx,my,nxv,nypmx,mx1,mxyp1) bind(C,name='cgpu2ppgppost2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, mx, my, nxv, nypmx, mx1
         integer(c_int), value :: mxyp1
         integer(c_int), value :: noff
         real(c_float), value :: qm
         type (c_ptr), value :: g_ppart, g_q, g_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpu2ppjppost2l(g_ppart,g_cu,g_kpic,noff,qm,dt,nppmx&
     &,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)                      &
     &bind(C,name='cgpu2ppjppost2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, mx, my, nxv
         integer(c_int), value :: nypmx, mx1, mxyp1, ipbc
         integer(c_int), value :: noff
         real(c_float), value :: qm, dt
         type (c_ptr), value :: g_ppart, g_cu, g_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpu2pprjppost2l(g_ppart,g_cu,g_kpic,noff,qm,dt,ci, &
     &nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)                 &
     &bind(C,name='cgpu2pprjppost2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, mx, my, nxv
         integer(c_int), value :: nypmx, mx1, mxyp1, ipbc
         integer(c_int), value :: noff
         real(c_float), value :: qm, dt, ci
         type (c_ptr), value :: g_ppart, g_cu, g_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcaguard2xl(g_qc,g_scs,g_q,nyp,nx,nxe,nypmx,   &
     &nxvh,kypd) bind(C,name='cgpuppcaguard2xl')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, nxe, nypmx, nxvh, kypd
         integer(c_int), value :: nyp
         type (c_ptr), value :: g_qc, g_scs, g_q
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcaguard2yl(g_qc,g_scr,nx,nxvh,kypd)           &
     & bind(C,name='cgpuppcaguard2yl')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, nxvh, kypd
         type (c_ptr), value :: g_qc, g_scr
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcacguard2xl(g_cuc,g_scs,g_cu,nyp,nx,nxe,nypmx,&
     &nxvh,kypd) bind(C,name='cgpuppcacguard2xl')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, nxe, nypmx, nxvh, kypd
         integer(c_int), value :: nyp
         type (c_ptr), value :: g_cuc, g_scs, g_cu
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcacguard2yl(g_fvc,g_scr,nx,nxvh,kypd)         &
     & bind(C,name='cgpuppcacguard2yl')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, nxvh, kypd
         type (c_ptr), value :: g_fvc, g_scr
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcbguard2xl(g_fxyc,g_scs,g_fxy,nyp,nx,nxe,nypmx&
     &,nxvh,kypd) bind(C,name='cgpuppcbguard2xl')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, nxe, nypmx, nxvh, kypd
         integer(c_int), value :: nyp
         type (c_ptr), value :: g_fxyc, g_scs, g_fxy
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcbguard2yl(g_fxy,g_scr,nyp,nx,nxe,nxvh,nypmx) &
     &bind(C,name='cgpuppcbguard2yl')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, nxe, nxvh, nypmx
         integer(c_int), value :: nyp
         type (c_ptr), value :: g_fxy, g_scr
         end subroutine
      end interface
!
      interface
         subroutine cgpupppord2la(g_ppart,g_ppbuff,g_sbufl,g_sbufr,     &
     &g_kpic,g_ncl,g_ihole,g_ncll,g_nclr,noff,nyp,idimp,nppmx,nx,ny,mx, &
     &my,mx1,myp1,npbmx,ntmax,nbmax,g_irc) bind(C,name='cgpupppord2la')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, mx1
         integer(c_int), value ::myp1, npbmx, ntmax, nbmax
         integer(c_int), value:: noff, nyp
         type (c_ptr), value :: g_ppart, g_ppbuff, g_sbufl, g_sbufr
         type (c_ptr), value :: g_kpic, g_ncl, g_ihole, g_ncll, g_nclr
         type (c_ptr), value :: g_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpupppord2lb(g_ppart,g_ppbuff,g_rbufl,g_rbufr,     &
     &g_kpic,g_ncl,g_ihole,g_mcll,g_mclr,idimp,nppmx,mx1,myp1,npbmx,    &
     &ntmax,nbmax,g_irc) bind(C,name='cgpupppord2lb')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, mx1, myp1, npbmx, ntmax
         integer(c_int), value :: nbmax
         type (c_ptr), value :: g_ppart, g_ppbuff, g_rbufl, g_rbufr
         type (c_ptr), value :: g_kpic, g_ncl, g_ihole, g_mcll, g_mclr
         type (c_ptr), value :: g_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpuppois23t(g_qt,g_fxyt,g_ffct,g_we,nx,ny,kstrt,nyv&
     &,kxp1,nyhd) bind(C,name='cgpuppois23t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, kstrt, nyv, kxp1, nyhd
         type (c_ptr), value :: g_qt, g_fxyt, g_ffct, g_we
         end subroutine
      end interface
!
      interface
         subroutine cgpuppcuperp2t(g_cut,nx,ny,kstrt,nyv,kxp1)          &
     &bind(C,name='cgpuppcuperp2t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, kstrt, nyv, kxp1
         type (c_ptr), value :: g_cut
         end subroutine
      end interface
!
      interface
         subroutine cgpuippbpoisp23t(g_cut,g_bxyt,g_ffct,ci,g_wm,nx,ny, &
     &kstrt,nyv,kxp1,nyhd) bind(C,name='cgpuippbpoisp23t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, kstrt, nyv, kxp1, nyhd
         real(c_float), value :: ci
         type (c_ptr), value :: g_cut, g_bxyt, g_ffct, g_wm
         end subroutine
      end interface
!
      interface
         subroutine cgpuppmaxwel2t(g_exyt,g_bxyt,g_cut,g_ffct,affp,ci,dt&
     &,g_wf,g_wm,nx,ny,kstrt,nyv,kxp1,nyhd)                             &
     &bind(C,name='cgpuppmaxwel2t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, kstrt, nyv, kxp1, nyhd
         real(c_float), value :: affp, ci, dt
         type (c_ptr), value :: g_exyt, g_bxyt, g_cut, g_ffct
         type (c_ptr), value :: g_wf, g_wm
         end subroutine
      end interface
!
      interface
         subroutine cgpuppemfield2t(g_fxyt,g_exyt,g_ffct,isign,nx,ny,   &
     &kstrt,nyv,kxp1,nyhd) bind(C,name='cgpuppemfield2t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, nx, ny, kstrt, nyv, kxp1, nyhd
         type (c_ptr), value :: g_fxyt, g_exyt, g_ffct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwppfft2rcsx(g_f,g_bsm,isign,g_mixup,g_sct,indx, &
     &indy,kstrt,nvp,kxp1,kyp,nxhd,kypd,nxhyd,nxyhd)                    &
     &bind(C,name='cgpuwppfft2rcsx')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, kstrt, nvp, kxp1
         integer(c_int), value :: kyp, nxhd, kypd, nxhyd, nxyhd
         type (c_ptr), value :: g_f, g_bsm, g_mixup, g_sct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwppfft2rcsy(g_g,g_brm,isign,g_mixup,g_sct,indx, &
     &indy,kstrt,nvp,kxp1,kyp,nyd,nxhyd,nxyhd)                          &
     &bind(C,name='cgpuwppfft2rcsy')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, kstrt, nvp, kxp1
         integer(c_int), value :: kyp, nyd, nxhyd, nxyhd
         type (c_ptr), value :: g_g, g_brm, g_mixup, g_sct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwppfft2rcsxn(g_fn,g_bsm,isign,g_mixup,g_sct,indx&
     &,indy,ndim,kstrt,nvp,kxp1,kyp,nxhd,kypd,nxhyd,nxyhd)              &
     &bind(C,name='cgpuwppfft2rcsxn')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, ndim, kstrt, nvp
         integer(c_int), value :: kxp1, kyp, nxhd, kypd, nxhyd, nxyhd
         type (c_ptr), value :: g_fn, g_bsm, g_mixup, g_sct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwppfft2rcsyn(g_gn,g_brm,isign,g_mixup,g_sct,indx&
     &,indy,ndim,kstrt,nvp,kxp1,kyp,nyd,nxhyd,nxyhd)                    &
     &bind(C,name='cgpuwppfft2rcsyn')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, ndim, kstrt, nvp
         integer(c_int), value :: kxp1, kyp, nyd, nxhyd, nxyhd
         type (c_ptr), value :: g_gn, g_brm, g_mixup, g_sct
         end subroutine
      end interface
!
      interface
         subroutine cgpuppltpose(g_f,g_g,nx,ny,kxp,kyp,kstrt,nxv,nyv)   &
     &bind(C,name='cgpuppltpose')
         use iso_c_binding
         implicit none
         integer(c_int), value:: nx, ny, kxp, kyp, kstrt, nxv, nyv
         type (c_ptr), value :: g_f, g_g
         end subroutine
      end interface
!
      interface
         subroutine cgpuppltposen(g_fn,g_gn,nx,ny,kxp,kyp,kstrt,ndim,nxv&
     &,nyv) bind(C,name='cgpuppltposen')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, kxp, kyp, kstrt, ndim, nxv
         integer(c_int), value :: nyv
         type (c_ptr), value :: g_fn, g_gn
         end subroutine
      end interface
!
      interface
         subroutine cgpusum2(g_a,g_sa,nx) bind(C,name='cgpusum2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx
         type (c_ptr), value :: g_a, g_sa
         end subroutine
      end interface
!
      end module


