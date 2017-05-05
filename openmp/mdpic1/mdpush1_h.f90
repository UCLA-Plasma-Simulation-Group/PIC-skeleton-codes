!-----------------------------------------------------------------------
! Interface file for mdpush1.f
      module mdpush1_h
      implicit none
!
      interface
         subroutine DISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,idimp,nop, &
     &nx,ipbc)
         implicit none
         integer, intent(in) :: npx, idimp, nop, nx, ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,nop), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine DBLKP1L(part,kpic,nppmx,idimp,nop,mx,mx1,irc)
         implicit none
         integer, intent(in) :: idimp, nop, mx, mx1
         integer, intent(inout) :: nppmx, irc
         real, dimension(idimp,nop), intent(in) :: part
         integer, dimension(mx1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN1L(part,ppart,kpic,nppmx,idimp,nop,mx,mx1,irc&
     &)
         implicit none
         integer, intent(in) :: nppmx, idimp, nop, mx, mx1
         integer, intent(inout) :: irc
         real, dimension(idimp,nop), intent(in) :: part
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         integer, dimension(mx1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPCHECK1L(ppart,kpic,idimp,nppmx,nx,mx,mx1,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, mx1
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,   &
     &idimp,nppmx,nx,mx,nxv,mx1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: omx, qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt&
     &,dtc,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: omx, qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GPPOST1L(ppart,q,kpic,qm,nppmx,idimp,mx,nxv,mx1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, nxv, mx1
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(nxv), intent(inout) :: q
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GJPPOST1L(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,mx,nxv,&
     &mx1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(2,nxv), intent(inout) :: cu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GMJPPOST1L(ppart,amu,kpic,qm,nppmx,idimp,mx,nxv,mx1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, nxv, mx1
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(2,nxv), intent(inout) :: amu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GDJPPOST1L(ppart,fxyz,byz,dcu,amu,kpic,omx,qm,qbm,dt&
     &,idimp,nppmx,nx,mx,nxv,mx1)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1
         real, intent(in) :: omx, qm, qbm, dt
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         real, dimension(2,nxv), intent(inout) :: dcu, amu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GDCJPPOST1L(ppart,fxyz,byz,cu,dcu,amu,kpic,omx,qm,  &
     &qbm,dt,idimp,nppmx,nx,mx,nxv,mx1)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1
         real, intent(in) :: omx, qm, qbm, dt
         real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(3,nxv), intent(in) :: fxyz
         real, dimension(2,nxv), intent(in) :: byz
         real, dimension(2,nxv), intent(inout) :: cu, dcu, amu
         integer, dimension(mx1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPORDER1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx&
     &,mx,mx1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, mx1, npbmx, ntmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1), intent(inout) :: ppbuff
         integer, dimension(mx1), intent(inout) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPORDERF1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &mx1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, npbmx, ntmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1), intent(inout) :: ppbuff
         integer, dimension(mx1), intent(inout) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(in) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine DGUARD1L(fx,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(nxe), intent(inout) :: fx
         end subroutine
      end interface
!
      interface
         subroutine CGUARD1L(byz,nx,nxe)
         implicit none
         integer , intent(in):: nx, nxe
         real, dimension(2,nxe), intent(inout) :: byz
         end subroutine
      end interface
!
      interface
         subroutine ACGUARD1L(cu,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(2,nxe), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine AGUARD1L(q,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, dimension(nxe), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine ASCFGUARD1L(dcu,cus,q2m0,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: q2m0
         real, dimension(2,nxe), intent(inout) :: dcu
         real, dimension(2,nxe), intent(in) :: cus
         end subroutine
      end interface
!
      interface
         subroutine FWPMINMX1(qe,qbme,wpmax,wpmin,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: qbme
         real, intent(inout) :: wpmax, wpmin
         real, dimension(nxe), intent(in) :: qe
         end subroutine
      end interface
!
      interface
         subroutine POIS1(q,fx,isign,ffc,ax,affp,we,nx)
         implicit none
         integer, intent(in) :: isign, nx
         real, intent(in) :: ax, affp
         real, intent(inout) :: we
         real, dimension(nx), intent(in) :: q
         real, dimension(nx), intent(inout) :: fx
         complex, dimension(nx/2), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine BBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
         implicit none
         integer, intent(in) :: nx, nxvh, nxhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(2,2*nxvh), intent(in) :: cu
         real, dimension(2,2*nxvh), intent(inout) :: byz
         complex, dimension(nxhd), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine BADDEXT1(byz,omy,omz,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
         real, intent(in) :: omy, omz
         real, dimension(2,nxe), intent(inout) :: byz
         end subroutine
      end interface
!
      interface
         subroutine DCUPERP13(dcu,amu,nx,nxvh)
         implicit none
         integer, intent(in) :: nx, nxvh
         real, dimension(2,2*nxvh), intent(inout) :: dcu
         real, dimension(2,2*nxvh), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine ADCUPERP13(dcu,amu,nx,nxvh)
         implicit none
         integer, intent(in) :: nx, nxvh
         real, dimension(2,2*nxvh), intent(inout) :: dcu
         real, dimension(2,2*nxvh), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine EPOIS13(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,nxvh,&
     &nxhd)
         implicit none
         integer, intent(in) :: isign, nx, nxvh, nxhd
         real, intent(in) :: ax, affp, wp0, ci
         real, intent(inout) :: wf
         real, dimension(2,2*nxvh), intent(in) :: dcu
         real, dimension(2,2*nxvh), intent(inout) :: eyz
         complex, dimension(nxhd), intent(inout) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine ADDVRFIELD13(fxyze,eyze,fxe,nxe)
         implicit none
         integer, intent(in) :: nxe
         real, dimension(3,nxe), intent(inout) :: fxyze
         real, dimension(2,nxe), intent(in) :: eyze
         real, dimension(nxe), intent(in) :: fxe
         end subroutine
      end interface
!
      interface
         subroutine WFFT1RINIT(mixup,sct,indx,nxhd)
         implicit none
         integer, intent(in) :: indx, nxhd
         integer, dimension(nxhd), intent(inout) :: mixup
         complex, dimension(nxhd), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer, intent(in) :: isign, indx, nxd, nxhd
         real, dimension(nxd), intent(inout) :: f, t
         integer, dimension(nxhd), intent(in) :: mixup
         complex, dimension(nxhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT1R2X(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer, intent(in) :: isign, indx, nxd, nxhd
         real, dimension(2,nxd), intent(inout) :: f, t
         integer, dimension(nxhd), intent(in) :: mixup
         complex, dimension(nxhd), intent(in) :: sct
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
