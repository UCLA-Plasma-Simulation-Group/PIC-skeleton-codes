!-----------------------------------------------------------------------
! Fortran interface to CUDA Library for Skeleton 2D Electrostatic GPU
! PIC Code */
! written by Viktor K. Decyk, UCLA
      module gpupush2_h
      implicit none
!
      interface
         subroutine cgpuppush2l(gp_ppart,gp_fxy,gp_kpic,qbm,dt,gp_ek,   &
     &idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ipbc
         real :: qbm, dt
         integer, dimension(2) :: gp_ppart, gp_fxy, gp_kpic, gp_ek
         end subroutine
      end interface
!
      interface
         subroutine cgpuppushf2l(gp_ppart,gp_fxy,gp_kpic,gp_ncl,gp_ihole&
     &,qbm,dt,gp_ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,     &
     &gp_irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
         integer :: ntmax
         real :: qbm, dt
         integer, dimension(2) :: gp_ppart, gp_fxy, gp_kpic, gp_ncl
         integer, dimension(2) :: gp_ihole, gp_ek, gp_irc
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
         subroutine cgpucaguard2l(gp_qc,gp_q,nx,ny,nxe,nye,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxe, nye, nxvh, nyv
         integer, dimension(2) :: gp_qc, gp_q
         end subroutine
      end interface
!
      interface
         subroutine cgpuccguard2l(gp_fxyc,gp_fxy,nx,ny,nxe,nye,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxe, nye, nxvh, nyv
         integer, dimension(2) :: gp_fxyc, gp_fxy
         end subroutine
      end interface
!
      interface
         subroutine cgpuppord2l(gp_ppart,gp_ppbuff,gp_kpic,gp_ncl,      &
     &gp_ihole,idimp,nppmx,nx,ny,mx,my,mx1,my1,npbmx,ntmax,gp_irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, mx1, my1, npbmx
         integer :: ntmax
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
         subroutine cgpupois22t(gp_qt,gp_fxyt,gp_ffct,gp_we,nx,ny,nxvh, &
     &nyv,nxhd,nyhd) 
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         integer, dimension(2) :: gp_qt, gp_fxyt, gp_ffct, gp_we
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


