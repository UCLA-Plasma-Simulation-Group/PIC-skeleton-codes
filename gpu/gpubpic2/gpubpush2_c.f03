!-----------------------------------------------------------------------
! Fortran2003 interface to CUDA Library for Skeleton 2-1/2D
! Electromagnetic GPU PIC Code
! written by Viktor K. Decyk, UCLA
      module gpubpush2_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine cgpubppush23l(g_ppart,g_fxy,g_bxy,g_kpic,qbm,dt,dtc,&
     &g_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)               &
     &bind(C,name='cgpubppush23l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ipbc
         real(c_float), value :: qbm, dt, dtc
         type (c_ptr), value :: g_ppart, g_fxy, g_bxy, g_kpic, g_ek
         end subroutine
      end interface
!
      interface
         subroutine cgpubppushf23l(g_ppart,g_fxy,g_bxy,g_kpic,g_ncl,    &
     &g_ihole,qbm,dt,dtc,g_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1, &
     &ntmax,g_irc) bind(C,name='cgpubppushf23l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ntmax
         real(c_float), value :: qbm, dt, dtc
         type (c_ptr), value :: g_ppart, g_fxy, g_bxy, g_kpic, g_ncl
         type (c_ptr), value :: g_ihole, g_ek, g_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpurbppush23l(g_ppart,g_fxy,g_bxy,g_kpic,qbm,dt,dtc&
     &,ci,g_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)           &
     &bind(C,name='cgpurbppush23l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ipbc
         real(c_float), value :: qbm, dt, dtc, ci
         type (c_ptr), value :: g_ppart, g_fxy, g_bxy, g_kpic, g_ek
         end subroutine
      end interface
!
      interface
         subroutine cgpurbppushf23l(g_ppart,g_fxy,g_bxy,g_kpic,g_ncl,   &
     &g_ihole,qbm,dt,dtc,ci,g_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,   &
     &mxy1,ntmax,g_irc) bind(C,name='cgpurbppushf23l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ntmax
         real(c_float), value :: qbm, dt, dtc, ci
         type (c_ptr), value :: g_ppart, g_fxy, g_bxy, g_kpic, g_ncl
         type (c_ptr), value :: g_ihole, g_ek, g_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpu2ppost2l(g_ppart,g_q,g_kpic,qm,nppmx,idimp,mx,my&
     &,nxv,nyv,mx1,mxy1) bind(C,name='cgpu2ppost2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1
         real(c_float), value :: qm
         type (c_ptr), value :: g_ppart, g_q, g_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpu2jppost2l(g_ppart,g_cu,g_kpic,qm,dt,nppmx,idimp,&
     &nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc) bind(C,name='cgpu2jppost2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ipbc
         real(c_float), value :: qm, dt
         type (c_ptr), value :: g_ppart, g_cu, g_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpu2jppostf2l(g_ppart,g_cu,g_kpic,g_ncl,g_ihole,qm,&
     &dt,nppmx,idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,g_irc)          &
     &bind(C,name='cgpu2jppostf2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ntmax
         real(c_float), value :: qm, dt
         type (c_ptr), value :: g_ppart, g_cu, g_kpic, g_ncl, g_ihole
         type (c_ptr), value :: g_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpu2rjppost2l(g_ppart,g_cu,g_kpic,qm,dt,ci,nppmx,  &
     &idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)                          &
     &bind(C,name='cgpu2rjppost2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ipbc
         real(c_float), value :: qm, dt, ci
         type (c_ptr), value :: g_ppart, g_cu, g_kpic
         end subroutine
      end interface
!
      interface
         subroutine cgpu2rjppostf2l(g_ppart,g_cu,g_kpic,g_ncl,g_ihole,qm&
     &,dt,ci,nppmx,idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,g_irc)      &
     &bind(C,name='cgpu2rjppostf2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ntmax
         real(c_float), value :: qm, dt, ci
         type (c_ptr), value :: g_ppart, g_cu, g_kpic, g_ncl, g_ihole
         type (c_ptr), value :: g_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpucaguard2l(g_qc,g_q,nx,ny,nxe,nye,nxvh,nyv)      &
     &bind(C,name='cgpucaguard2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxe, nye, nxvh, nyv
         type (c_ptr), value :: g_qc, g_q
         end subroutine
      end interface
!
      interface
         subroutine cgpucacguard2l(g_cuc,g_cu,nx,ny,nxe,nye,nxvh,nyv)   &
     &bind(C,name='cgpucacguard2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxe, nye, nxvh, nyv
         type (c_ptr), value :: g_cuc, g_cu
         end subroutine
      end interface
!
      interface
         subroutine cgpucbguard2l(g_bxyc,g_bxy,nx,ny,nxe,nye,nxvh,nyv)  &
     &bind(C,name='cgpucbguard2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxe, nye, nxvh, nyv
         type (c_ptr), value :: g_bxyc, g_bxy
         end subroutine
      end interface
!
      interface
         subroutine cgpuppord2l(g_ppart,g_ppbuff,g_kpic,g_ncl,g_ihole,  &
     &idimp,nppmx,nx,ny,mx,my,mx1,my1,npbmx,ntmax,g_irc)                &
     &bind(C,name='cgpuppord2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, mx1, my1
         integer(c_int), value :: npbmx, ntmax
         type (c_ptr), value :: g_ppart, g_ppbuff, g_kpic, g_ncl
         type (c_ptr), value :: g_ihole, g_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpuppordf2l(g_ppart,g_ppbuff,g_kpic,g_ncl,g_ihole, &
     &idimp,nppmx,mx1,my1,npbmx,ntmax,g_irc) bind(C,name='gpuppordf2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, mx1, my1, npbmx, ntmax
         type (c_ptr), value :: g_ppart, g_ppbuff, g_kpic, g_ncl
         type (c_ptr), value :: g_ihole, g_irc
         end subroutine
      end interface
!
      interface
         subroutine cgpupois23t(g_qt,g_fxyt,g_ffct,g_we,nx,ny,nxvh,nyv, &
     &nxhd,nyhd) bind(C,name='cgpupois23t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxvh, nyv, nxhd, nyhd
         type (c_ptr), value :: g_qt, g_fxyt, g_ffct, g_we
         end subroutine
      end interface
!
      interface
         subroutine cgpucuperp2t(g_cut,nx,ny,nxvh,nyv)                  &
     &bind(C,name='cgpucuperp2t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxvh, nyv
         type (c_ptr), value :: g_cut
         end subroutine
      end interface
!
      interface
         subroutine cgpuibpois23t(g_cut,g_bxyt,g_ffct,ci,g_wm,nx,ny,nxvh&
     &,nyv,nxhd,nyhd) bind(C,name='cgpuibpois23t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxvh, nyv, nxhd, nyhd
         real(c_float), value :: ci
         type (c_ptr), value :: g_cut, g_bxyt, g_ffct, g_wm
         end subroutine
      end interface
!
      interface
         subroutine cgpumaxwel2t(g_exyt,g_bxyt,g_cut,g_ffct,ci,dt,g_wf, &
     &g_wm,nx,ny,nxvh,nyv,nxhd,nyhd) bind(C,name='cgpumaxwel2t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxvh, nyv, nxhd, nyhd
         real(c_float), value :: ci, dt
         type (c_ptr), value :: g_exyt, g_bxyt, g_cut, g_ffct
         type (c_ptr), value :: g_wf, g_wm
         end subroutine
      end interface
!
      interface
         subroutine cgpuemfield2t(g_fxyt,g_exyt,g_ffct,isign,nx,ny,nxvh,&
     &nyv,nxhd,nyhd) bind(C,name='cgpuemfield2t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         type (c_ptr), value :: g_fxyt, g_exyt, g_ffct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwfft2rcs(g_f,g_g,isign,g_mixup,g_sct,indx,indy, &
     &nxhd,nyd,nxhyd,nxyhd) bind(C,name='cgpuwfft2rcs')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, nxhd, nyd
         integer(c_int), value :: nxhyd, nxyhd
         type (c_ptr), value :: g_f, g_g, g_mixup, g_sct
         end subroutine
      end interface
!
      interface
         subroutine cgpuwfft2rcsn(g_fn,g_gn,isign,g_mixup,g_sct,indx,   &
     &indy,ndim,nxhd,nyd,nxhyd,nxyhd) bind(C,name='cgpuwfft2rcsn')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, ndim, nxhd, nyd
         integer(c_int), value :: nxhyd, nxyhd
         type (c_ptr), value :: g_fn, g_gn, g_mixup, g_sct
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


