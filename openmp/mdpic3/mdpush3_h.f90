!-----------------------------------------------------------------------
! Interface file for mdpush3.f
      module mdpush3_h
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
         subroutine DBLKP3L(part,kpic,nppmx,idimp,nop,mx,my,mz,mx1,my1, &
     &mxyz1,irc)
         implicit none
         integer, intent(in) :: idimp, nop, mx, my, mz, mx1, my1, mxyz1
         integer, intent(inout) :: nppmx, irc
         real, dimension(idimp,nop), intent(in) :: part
         integer, dimension(mxyz1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN3L(part,ppart,kpic,nppmx,idimp,nop,mx,my,mz, &
     &mx1,my1,mxyz1,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nop, mx, my, mz
         integer, intent(in) :: mx1, my1, mxyz1
         integer, intent(inout) :: irc
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         integer, dimension(mxyz1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPCHECK3L(ppart,kpic,idimp,nppmx,nx,ny,nz,mx,my,mz, &
     &mx1,my1,mz1,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: mx1, my1, mz1
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*my1*mz1), intent(in) :: ppart
         integer, dimension(mx1*my1*mz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSH3L(ppart,fxyz,bxyz,kpic,qbm,dt,dtc,ek,idimp, &
     &nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSHF3L(ppart,fxyz,bxyz,kpic,ncl,ihole,qbm,dt,dtc&
     &,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax,&
     &irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GPPOST3L(ppart,q,kpic,qm,nppmx,idimp,mx,my,mz,nxv,  &
     &nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, my, mz, nxv, nyv, nzv
         integer, intent(in) :: mx1, my1, mxyz1
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mxyz1), intent(in) :: ppart
         real, dimension(nxv,nyv,nzv), intent(inout) :: q
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GJPPOST3L(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny,nz, &
     &mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GMJPPOST3L(ppart,amu,kpic,qm,nppmx,idimp,mx,my,mz,  &
     &nxv,nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, my, mz, nxv, nyv, nzv
         integer, intent(in) :: mx1, my1, mxyz1
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mxyz1), intent(in) :: ppart
         real, dimension(6,nxv,nyv,nzv), intent(inout) :: amu
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GDJPPOST3L(ppart,fxyz,bxyz,kpic,dcu,amu,qm,qbm,dt,  &
     &idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1
         real, intent(in) :: qm, qbm, dt
         real, dimension(idimp,nppmx,mxyz1), intent(in) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: dcu
         real, dimension(6,nxv,nyv,nzv), intent(inout) :: amu
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GDCJPPOST3L(ppart,fxyz,bxyz,kpic,cu,dcu,amu,qm,qbm, &
     &dt,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1
         real, intent(in) :: qm, qbm, dt
         real, dimension(idimp,nppmx,mxyz1), intent(in) :: ppart
         real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz, bxyz
         real, dimension(3,nxv,nyv,nzv), intent(inout) :: cu, dcu
         real, dimension(6,nxv,nyv,nzv), intent(inout) :: amu
         integer, dimension(mxyz1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPORDER3L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx&
     &,ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: mx1, my1, mz1, npbmx, ntmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*my1*mz1), intent(inout) ::     &
     &ppart
         real, dimension(idimp,npbmx,mx1*my1*mz1), intent(inout) ::     &
     &ppbuff
         integer, dimension(mx1*my1*mz1), intent(inout) :: kpic
         integer, dimension(26,mx1*my1*mz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1), intent(inout) ::    &
     &ihole
         end subroutine
      end interface
!
      interface
         subroutine PPORDERF3L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in):: idimp, nppmx, mx1, my1, mz1, npbmx, ntmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*my1*mz1), intent(inout) ::     &
     &ppart
         real, dimension(idimp,npbmx,mx1*my1*mz1), intent(inout) ::     &
     &ppbuff
         integer, dimension(mx1*my1*mz1), intent(inout) :: kpic
         integer, dimension(26,mx1*my1*mz1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1*mz1), intent(in) :: ihole
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
         integer, intent(in):: nx, ny, nz, nxe, nye, nze
         real, intent(in) :: qbme
         real, intent(inout) :: wpmax, wpmin
         real, dimension(nxe,nye,nze), intent(in) :: qe
         end subroutine
      end interface
!
      interface
         subroutine MPOIS33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz, &
     &nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ax, ay, az, affp
         real, intent(inout) :: we
         real, dimension(2*nxvh,nyv,nzv), intent(in) :: q
         real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: fxyz
         complex, dimension(nxhd,nyhd,nzhd), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MCUPERP3(cu,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         real, dimension(3,2*nxvh,nyv,nzv) :: cu
         end subroutine
      end interface
!
      interface
         subroutine MBBPOIS33(cu,bxyz,ffc,ci,wm,nx,ny,nz,nxvh,nyv,nzv,  &
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
         subroutine MDCUPERP3(dcu,amu,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: dcu
         real, dimension(6,2*nxvh,nyv,nzv), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MADCUPERP3(dcu,amu,nx,ny,nz,nxvh,nyv,nzv)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxvh, nyv, nzv
         real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: dcu
         real, dimension(6,2*nxvh,nyv,nzv), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MEPOIS33(dcu,exyz,isign,ffe,ax,ay,az,affp,wp0,ci,wf,&
     &nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ax, ay, az, affp, wp0, ci
         real, intent(inout) :: wf
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
         subroutine WFFT3RMX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd, &
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
         subroutine WFFT3RM3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd, &
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
         subroutine WFFT3RMN(f,ss,isign,mixup,sct,indx,indy,indz,nxhd,  &
     &nyd,nzd,ndim,nxhyzd,nxyzhd)
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
          subroutine FFT3RMXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp, &
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
         subroutine FFT3RMXZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,  &
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
         subroutine FFT3RM3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp, &
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
         subroutine FFT3RM3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,  &
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
         subroutine FFT3RMNXY(f,ss,isign,mixup,sct,indx,indy,indz,nzi,  &
     &nzp,nxhd,nyd,nzd,ndim,nxhyzd,nxyzhd)
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
         subroutine FFT3RMNZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,  &
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
         subroutine MSWAP3CN(f,s,isign,nxh,ny,nzi,nzt,nxhd,nyd,nzd,ndim)
         implicit none
         integer, intent(in) :: isign, nxh, ny, nzi, nzt, nxhd, nyd, nzd
         integer, intent(in) :: ndim
         real, dimension(ndim,2*nxhd,nyd,nzd), intent(inout) :: f
         complex, dimension(ndim*nxhd,nzd), intent(inout) :: s
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

