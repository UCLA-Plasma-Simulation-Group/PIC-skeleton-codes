!-----------------------------------------------------------------------
! Fortran2003 interface to CUDA Library for Skeleton 2D Electrostatic
! GPU PIC Code */
! written by Viktor K. Decyk, UCLA
      module gpupush2_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine cgpuppush2l(g_ppart,g_fxy,g_kpic,qbm,dt,g_ek,idimp, &
     &nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc) bind(C,name='cgpuppush2l'&
     &)
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ipbc
         real(c_float), value :: qbm, dt
         type (c_ptr), value :: g_ppart, g_fxy, g_kpic, g_ek
         end subroutine
      end interface
!
      interface
         subroutine cgpuppushf2l(g_ppart,g_fxy,g_kpic,g_ncl,g_ihole,qbm,&
     &dt,g_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,g_irc)     &
     &bind(C,name='cgpuppushf2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer(c_int), value :: mx1, mxy1, ntmax
         real(c_float), value :: qbm, dt
         type (c_ptr), value :: g_ppart, g_fxy, g_kpic, g_ncl, g_ihole
         type (c_ptr), value :: g_ek, g_irc
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
         subroutine cgpuccguard2l(g_fxyc,g_fxy,nx,ny,nxe,nye,nxvh,nyv)  &
     &bind(C,name='cgpuccguard2l')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxe, nye, nxvh, nyv
         type (c_ptr), value :: g_fxyc, g_fxy
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
         subroutine cgpupois22t(g_qt,g_fxyt,g_ffct,g_we,nx,ny,nxvh,nyv, &
     &nxhd,nyhd) bind(C,name='cgpupois22t')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, nxvh, nyv, nxhd, nyhd
         type (c_ptr), value :: g_qt, g_fxyt, g_ffct, g_we
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


