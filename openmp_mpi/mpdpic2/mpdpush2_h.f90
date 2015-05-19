!-----------------------------------------------------------------------
! Interface file for mpdpush2.f
      module mpdpush2_h
      implicit none
!
      interface
         subroutine PDICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,  &
     &idps)
         implicit none
         integer, intent(in) :: ny, kstrt, nvp, idps
         integer, intent(inout) :: nyp, noff, nypmx, nypmn
         real, dimension(idps), intent(inout) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PDISTR2H(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,&
     &npx,npy,nx,ny,idimp,npmax,idps,ipbc,ierr)
         implicit none
         integer, intent(in) :: nps, npx, npy, nx, ny, idimp, npmax
         integer, intent(in) :: idps, ipbc
         integer, intent(inout) :: npp, ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idps), intent(in) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my&
     &,mx1,mxyp1,irc)
         implicit none
         integer, intent(in) :: idimp, npmax, mx, my, mx1, mxyp1, npp
         integer, intent(in) :: noff
         integer, intent(inout) :: nppmx, irc
         real, dimension(idimp,npmax), intent(in) :: part
         integer, dimension(mxyp1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,    &
     &npmax,mx,my,mx1,mxyp1,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, npmax, mx, my, mx1, mxyp1
         integer, intent(in) :: npp, noff
         integer, intent(inout) :: irc
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         integer, dimension(mxyp1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,&
     &mx1,myp1,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, my, mx1, myp1
         integer, intent(in) :: noff, nyp
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(in) :: ppart
         integer, dimension(mx1*myp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,&
     &ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,&
     &qbm,dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,  &
     &irc)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout)  :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(8,mxyp1), intent(inout)  :: ncl
         integer, dimension(2,ntmax+1,mxyp1), intent(inout)  :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,  &
     &nxv, nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, idimp, nppmx, mx, my, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nxv,nypmx), intent(inout) :: q
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGJPPOST2L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,&
     &ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: qm, dt
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(inout) :: cu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGMJPPOST2L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx, &
     &my,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, mx, my, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1
         real, intent(in) ::  qm
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(4,nxv,nypmx), intent(inout) :: amu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,noff,nyp,qm,&
     &qbm,dt,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1
         real, intent(in) ::  qm, qbm, dt
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(3,nxv,nypmx), intent(inout) :: dcu
         real, dimension(4,nxv,nypmx), intent(inout) :: amu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,noff,nyp&
     &,qm,qbm,dt,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1
         real, intent(in) ::  qm, qbm, dt
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(3,nxv,nypmx), intent(inout) :: cu, dcu
         real, dimension(4,nxv,nypmx), intent(inout) :: amu
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,&  
     &ncll,nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,  &
     &nbmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, mx1, myp1
         integer, intent(in) :: npbmx, ntmax, nbmax, noff, nyp
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1), intent(inout) :: ppbuff
         real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr
         integer, dimension(mx1*myp1), intent(in) :: kpic
         integer, dimension(8,mx1*myp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1), intent(inout) :: ihole
         integer, dimension(3,mx1), intent(inout) :: ncll, nclr
         end subroutine
      end interface
!
      interface
         subroutine PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll&
     &,nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, myp1, npbmx, ntmax
         integer, intent(in) :: nbmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1), intent(inout) :: ppbuff
         real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr
         integer, dimension(8,mx1*myp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1), intent(in) :: ihole
         integer, dimension(3,mx1), intent(inout) :: ncll, nclr
         end subroutine
      end interface
!
      interface
         subroutine PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,&
     &mcll,mclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, myp1, npbmx, ntmax
         integer, intent(in) :: nbmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1), intent(in) :: ppbuff
         real, dimension(idimp,nbmax), intent(in) :: rbufl, rbufr
         integer, dimension(mx1*myp1), intent(inout) :: kpic
         integer, dimension(8,mx1*myp1), intent(in) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1), intent(in) :: ihole
         integer, dimension(3,mx1), intent(in) :: mcll, mclr
         end subroutine
      end interface
!
      interface
         subroutine PPCGUARD2XL(fxy,nyp,nx,ndim,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, ndim, nxe, nypmx
         real, dimension(ndim,nxe,nypmx), intent(inout) :: fxy
         end subroutine
      end interface
!
      interface
         subroutine PPAGUARD2XL(q,nyp,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, dimension(nxe,nypmx), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine PPACGUARD2XL(cu,nyp,nx,ndim,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, ndim, nxe, nypmx
         real, dimension(ndim,nxe,nypmx), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine PPASCFGUARD2L(dcu,cus,nyp,q2m0,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, intent(in) :: q2m0
         real, dimension(3,nxe,nypmx), intent(inout) :: dcu
         real, dimension(3,nxe,nypmx), intent(in) :: cus
         end subroutine
      end interface
!
      interface
         subroutine PPFWPMINMX2(qe,nyp,qbme,wpmax,wpmin,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, intent(in) :: qbme
         real, intent(inout) :: wpmax, wpmin
         real, dimension(nxe,nypmx), intent(in) :: qe
         end subroutine
      end interface
!
      interface
         subroutine MPPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt, &
     &nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(3,nyv,kxp), intent(inout) :: fxy
         complex, dimension(nyhd,kxp), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPCUPERP2(cu,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine MPPBBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,  &
     &nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nyv,kxp), intent(in) :: cu
         complex, dimension(3,nyv,kxp), intent(inout) :: bxy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine PPBADDEXT2(bxy,nyp,omx,omy,omz,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, intent(in) :: omx, omy, omz
         real, dimension(3,nxe,nypmx), intent(inout) :: bxy
         end subroutine
      end interface
!
      interface
         subroutine MPPDCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(inout) :: dcu
         complex, dimension(4,nyv,kxp), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MPPADCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(inout) :: dcu
         complex, dimension(4,nyv,kxp), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MPPEPOISP23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf, &
     &nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ax, ay, affp, wp0, ci
         real, intent(inout) :: wf
         complex, dimension(3,nyv,kxp), intent(in) :: dcu
         complex, dimension(3,nyv,kxp), intent(inout) :: exy
         complex, dimension(nyhd,kxp), intent(inout) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine PPADDVRFIELD2(a,b,c,ndim,nxe,nypmx)
         implicit none
         integer, intent(in) :: ndim, nxe, nypmx
         real, dimension(ndim,nxe,nypmx), intent(inout) :: a
         real, dimension(ndim,nxe,nypmx), intent(in) :: b, c
         end subroutine
      end interface
!
      interface
         subroutine WPFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: indx, indy, nxhyd, nxyhd
         integer, dimension(nxhyd), intent(inout) :: mixup
         complex, dimension(nxyhd), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WPPFFT2RM(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp, kyp, kypd, nxhyd, nxyhd
         real, intent(inout) :: ttp
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         complex, dimension(nyv,kxp), intent(inout) :: g
         complex, dimension(kxp,kyp), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WPPFFT2RM3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx&
     &,indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp, kyp, kypd, nxhyd, nxyhd
         real, intent(inout) :: ttp
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         complex, dimension(3,nyv,kxp), intent(inout) :: g
         complex, dimension(3,kxp,kyp), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WPPFFT2RMN(f,g,bs,br,ss,isign,ntpose,mixup,sct,ttp, &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,ndim,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp, kyp, kypd, ndim
         integer, intent(in) :: nxhyd, nxyhd
         real, intent(inout) :: ttp
         real, dimension(ndim,2*nxvh,kypd), intent(inout) :: f
         complex, dimension(ndim,nyv,kxp), intent(inout) :: g
         complex, dimension(ndim,kxp,kyp), intent(inout) :: bs, br
         complex, dimension(ndim*nxvh,kypd), intent(inout) :: ss
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RMXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,  &
     &kypp,nxvh,kypd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kypi, kypp
         integer, intent(in) :: nxvh, kypd, nxhyd, nxyhd
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RMXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,  &
     &kxpp,nyv,kxp,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kxpi, kxpp
         integer, intent(in) :: nyv, kxp, nxhyd, nxyhd
         complex, dimension(nyv,kxp), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RM3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kypi, kypp
         integer, intent(in) :: nxvh, kypd, nxhyd, nxyhd
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RM3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxp,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kxpi, kxpp
         integer, intent(in) :: nyv, kxp, nxhyd, nxyhd
         complex, dimension(3,nyv,kxp), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RMNXX(f,ss,isign,mixup,sct,indx,indy,kstrt,   &
     &kypi,kypp,nxvh,kypd,ndim,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kypi, kypp
         integer, intent(in) :: nxvh, kypd, ndim, nxhyd, nxyhd
         real, dimension(ndim,2*nxvh,kypd), intent(inout) :: f
         complex, dimension(ndim*nxvh,kypd), intent(inout) :: ss
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RMNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxp,ndim,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kxpi, kxpp
         integer, intent(in) :: nyv, kxp, ndim, nxhyd, nxyhd
         complex, dimension(ndim,nyv,kxp), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyhd), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine MPPSWAPC2N(f,s,isign,nxh,kypi,kypt,nxvh,kypd,ndim)
         implicit none
         integer, intent(in) :: isign, nxh, kypi, kypt, nxvh, kypd, ndim
         real, dimension(ndim,2*nxvh,kypd), intent(inout) :: f
         complex, dimension(ndim*nxvh,kypd), intent(inout) :: s
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