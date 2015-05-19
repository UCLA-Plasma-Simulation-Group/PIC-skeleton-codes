!-----------------------------------------------------------------------
! Interface file for mdpush2.f
      module mdpush2_h
      implicit none
!
      interface
         subroutine DISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp, &
     &nop,nx,ny,ipbc)
         implicit none
         integer, intent(in) :: npx, npy, idimp, nop, nx, ny, ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine DBLKP2L(part,kpic,nppmx,idimp,nop,mx,my,mx1,mxy1,irc&
     &)
         implicit none
         integer, intent(in) :: idimp, nop, mx, my, mx1, mxy1
         integer, intent(inout) :: nppmx, irc
         real, dimension(idimp,nop), intent(in) :: part
         integer, dimension(mxy1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN2L(part,ppart,kpic,nppmx,idimp,nop,mx,my,mx1,&
     &mxy1,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nop, mx, my, mx1, mxy1
         integer, intent(inout) :: irc
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         integer, dimension(mxy1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPCHECK2L(ppart,kpic,idimp,nppmx,nx,ny,mx,my,mx1,my1&
     &,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, mx1, my1
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*my1), intent(in) :: ppart
         integer, dimension(mx1*my1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSH23L(ppart,fxy,bxy,kpic,qbm,dt,dtc,ek,idimp,  &
     &nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         integer, dimension(mxy1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,dtc,&
     &ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         integer, dimension(mxy1), intent(in) :: kpic
         integer, dimension(8,mxy1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxy1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GPPOST2L(ppart,q,kpic,qm,nppmx,idimp,mx,my,nxv,nyv, &
     &mx1,mxy1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mxy1), intent(in) :: ppart
         real, dimension(nxv,nyv), intent(inout) :: q
         integer, dimension(mxy1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GJPPOST2L(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny,mx, &
     &my, nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ipbc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(3,nxv,nyv), intent(inout) :: cu
         integer, dimension(mxy1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GMJPPOST2L(ppart,amu,kpic,qm,nppmx,idimp,mx,my,nxv, &
     &nyv,mx1,mxy1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mxy1), intent(in) :: ppart
         real, dimension(4,nxv,nyv), intent(inout) :: amu
         integer, dimension(mxy1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,qm,qbm,dt,    &
     &idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1
         real, intent(in) :: qm, qbm, dt
         real, dimension(idimp,nppmx,mxy1), intent(in) :: ppart
         real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         real, dimension(3,nxv,nyv), intent(inout) :: dcu
         real, dimension(4,nxv,nyv), intent(inout) :: amu
         integer, dimension(mxy1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,qm,qbm,dt,&
     &idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1
         real, intent(in) :: qm, qbm, dt
         real, dimension(idimp,nppmx,mxy1), intent(in) :: ppart
         real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         real, dimension(3,nxv,nyv), intent(inout) :: cu, dcu
         real, dimension(4,nxv,nyv), intent(inout) :: amu
         integer, dimension(mxy1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPORDER2L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx&
     &,ny,mx,my,mx1,my1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, mx1, my1
         integer, intent(in) :: npbmx, ntmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*my1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1*my1), intent(inout) :: ppbuff
         integer, dimension(mx1*my1), intent(inout) :: kpic
         integer, dimension(8,mx1*my1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPORDERF2L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &mx1,my1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, my1, npbmx, ntmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*my1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1*my1), intent(inout) :: ppbuff
         integer, dimension(mx1*my1), intent(inout) :: kpic
         integer, dimension(8,mx1*my1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*my1), intent(in) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine BGUARD2L(bxy,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye), intent(inout) :: bxy
         end subroutine
      end interface
!
      interface
         subroutine ACGUARD2L(cu,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine AGUARD2L(q,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
         real, dimension(nxe,nye), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine AMCGUARD2L(amu,nx,ny,nxe,nye,ndim)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye, ndim
         real, dimension(ndim,nxe,nye), intent(inout) :: amu
         end subroutine
      end interface
!
      interface
         subroutine ASCFGUARD2L(dcu,cus,q2m0,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
         real, intent(in) :: q2m0
         real, dimension(3,nxe,nye), intent(inout) :: dcu
         real, dimension(3,nxe,nye), intent(in) :: cus
         end subroutine
      end interface
!
      interface
         subroutine FWPMINMX2(qe,qbme,wpmax,wpmin,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
         real, intent(in) :: qbme
         real, intent(inout) :: wpmax, wpmin
         real, dimension(nxe,nye), intent(in) :: qe
         end subroutine
      end interface
!
      interface
         subroutine MPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv&
     &,nxhd,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         real, dimension(2*nxvh,nyv), intent(in) :: q
         real, dimension(3,2*nxvh,nyv), intent(inout) :: fxy
         complex, dimension(nxhd,nyhd), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MCUPERP2(cu,nx,ny,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv
         real, dimension(3,2*nxvh,nyv), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine MBBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(3,2*nxvh,nyv), intent(in) :: cu
         real, dimension(3,2*nxvh,nyv), intent(inout) :: bxy
         complex, dimension(nxhd,nyhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine BADDEXT2(bxy,omx,omy,omz,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
         real, intent(in) :: omx, omy, omz
         real, dimension(3,nxe,nye), intent(inout) :: bxy
         end subroutine
      end interface
!
      interface
         subroutine MDCUPERP23(dcu,amu,nx,ny,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv
         real, dimension(3,2*nxvh,nyv), intent(inout) :: dcu
         real, dimension(4,2*nxvh,nyv), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MADCUPERP23(dcu,amu,nx,ny,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv
         real, dimension(3,2*nxvh,nyv), intent(inout) :: dcu
         real, dimension(4,2*nxvh,nyv), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MEPOIS23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx, &
     &ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ax, ay, affp, wp0, ci
         real, intent(inout) :: wf
         real, dimension(3,2*nxvh,nyv), intent(in) :: dcu
         real, dimension(3,2*nxvh,nyv), intent(inout) :: exy
         complex, dimension(nxhd,nyhd), intent(inout) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine ADDVRFIELD2(a,b,c,ndim,nxe,nye)
         implicit none
         integer, intent(in) :: nxe, nye, ndim
         real, dimension(ndim,nxe,nye), intent(inout) :: a
         real, dimension(ndim,nxe,nye), intent(in) :: b, c
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: indx, indy, nxhyd, nxyhd
         integer, dimension(nxhyd), intent(inout) :: mixup
         complex, dimension(nxyhd), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RMX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RM3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RMN(f,ss,isign,mixup,sct,indx,indy,nxhd,nyd,   &
     &ndim,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxhd, nyd, ndim
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(ndim,2*nxhd,nyd), intent(inout) :: f
         complex, dimension(ndim*nxhd,nyd), intent(inout) :: ss
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RM3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RM3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,  &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
         real, dimension(3,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RMNX(f,ss,isign,mixup,sct,indx,indy,nyi,nyp,nxhd&
     &,nyd,ndim,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: ndim, nxhyd, nxyhd
         real, dimension(ndim,2*nxhd,nyd), intent(inout) :: f
         complex, dimension(ndim*nxhd,nyd), intent(inout) :: ss
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT2RMNY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,  &
     &nyd,ndim,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: ndim, nxhyd, nxyhd
         real, dimension(ndim,2*nxhd,nyd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine MSWAPC2N(f,s,isign,nxh,nyi,nyt,nxhd,nyd,ndim)
         implicit none
         integer, intent(in) :: isign, nxh, nyi, nyt, nxhd, nyd, ndim
         real, dimension(ndim,2*nxhd,nyd), intent(inout) :: f
         complex, dimension(ndim*nxhd), intent(inout) :: s
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
