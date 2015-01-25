!-----------------------------------------------------------------------
! Fortran2003 interface to CUDA Library for Skeleton 2D Electrostatic
! GPU-MPI PIC Code */
! written by Viktor K. Decyk, UCLA
      module gpuppush2_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine cgpuppgppush2l(g_ppart,g_fxy,g_kpic,noff,nyp,qbm,dt,&
     &g_ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)            &
     &bind(C,name='cgpuppgppush2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, mx, my, idimp, nppmx, nxv
         integer(c_int), value :: nypmx, mx1, mxyp1, ipbc
         integer(c_int), value :: noff, nyp
         real(c_float), value :: qbm, dt
         type (c_ptr), value :: g_ppart, g_fxy, g_kpic, g_ek
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
         subroutine cgpuppccguard2xl(g_fxyc,g_scs,g_fxy,nyp,nx,nxe,nypmx&
     &,nxvh,kypd) bind(C,name='cgpuppccguard2xl')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, nxe, nypmx, nxvh, kypd
         integer(c_int), value :: nyp
         type (c_ptr), value :: g_fxyc, g_scs, g_fxy
         end subroutine
      end interface
!
      interface
         subroutine cgpuppccguard2yl(g_fxy,g_scr,nyp,nx,nxe,nxvh,nypmx) &
     &bind(C,name='cgpuppccguard2yl')
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
         subroutine cgpuppois22t(g_qt, g_fxyt,g_ffct,g_we,nx,ny,kstrt,  &
     &nyv,kxp1,nyhd) bind(C,name='cgpuppois22t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, kstrt, nyv, kxp1, nyhd
         type (c_ptr), value :: g_qt, g_fxyt, g_ffct, g_we
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


