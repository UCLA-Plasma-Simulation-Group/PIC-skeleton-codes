!-----------------------------------------------------------------------
! Interface file for mppush2.f
      module mppush2_h
      implicit none
!
      interface
         subroutine PDICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,  &
     &idps)
         implicit none
         integer :: nyp, noff, nypmx, nypmn, ny, kstrt, nvp, idps
         real, dimension(idps) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PDISTR2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy, &
     &nx,ny,idimp,npmax,idps,ipbc,ierr)
         implicit none
         integer :: npp, nps, npx, npy, nx, ny, idimp, npmax, idps, ipbc
         integer :: ierr
         real :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax) :: part
         real, dimension(idps) :: edges
         end subroutine
      end interface
!
      interface
         subroutine PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my&
     &,mx1,mxyp1,irc)
         implicit none
         integer :: nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
         integer :: npp, noff
         real, dimension(idimp,npmax) :: part
         integer, dimension(mxyp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,    &
     &npmax,mx,my,mx1,mxyp1,irc)
         implicit none
         integer :: nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
         integer :: npp, noff
         real, dimension(idimp,npmax) :: part
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         integer, dimension(mxyp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,&
     &mx1,myp1,irc)
         implicit none
         integer :: idimp, nppmx, nx, mx, my, mx1, myp1, irc
         integer :: noff, nyp
         real, dimension(idimp,nppmx,mx1*myp1) :: ppart
         integer, dimension(mx1*myp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGPPUSH2L(ppart,fxy,kpic,noff,nyp,qbm,dt,ek,nx,ny, &
     &mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer :: nx, ny, mx, my, idimp, nppmx, nxv, nypmx, mx1, mxyp1
         integer :: ipbc
         integer :: noff, nyp
         real :: qbm, dt, ek
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(2,nxv,nypmx) :: fxy
         integer, dimension(mxyp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGPPUSHF2L(ppart,fxy,kpic,ncl,ihole,noff,nyp,qbm,dt&
     &,ek,nx,ny,mx,my,idimp,nppmx,nxv,nypmx,mx1,mxyp1,ntmax,irc)
         implicit none
         integer :: nx, ny, mx, my, idimp, nppmx, nxv, nypmx, mx1, mxyp1
         integer :: ntmax, irc
         integer :: noff, nyp
         real :: qbm, dt, ek
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(2,nxv,nypmx) :: fxy
         integer, dimension(mxyp1) :: kpic
         integer, dimension(8,mxyp1) :: ncl
         integer, dimension(2,ntmax+1,mxyp1) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,  &
     &nxv,nypmx,mx1,mxyp1)
         implicit none
         integer :: idimp, nppmx, mx, my, nxv, nypmx, mx1, mxyp1
         integer :: noff
         real :: qm
         real, dimension(idimp,nppmx,mxyp1) :: ppart
         real, dimension(nxv,nypmx) :: q
         integer, dimension(mxyp1) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,&  
     &ncll,nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,  &
     &nbmax,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, mx1, myp1
         integer :: npbmx, ntmax, nbmax, irc
         integer :: noff, nyp
         real, dimension(idimp,nppmx,mx1*myp1) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1) :: ppbuff
         real, dimension(idimp,nbmax) :: sbufl, sbufr
         integer, dimension(mx1*myp1) :: kpic
         integer, dimension(8,mx1*myp1) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1) :: ihole
         integer, dimension(3,mx1) :: ncll, nclr
         end subroutine
      end interface
!
      interface
         subroutine PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll&
     &,nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, myp1, npbmx, ntmax, nbmax, irc
         real, dimension(idimp,nppmx,mx1*myp1) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1) :: ppbuff
         real, dimension(idimp,nbmax) :: sbufl, sbufr
         integer, dimension(8,mx1*myp1) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1) :: ihole
         integer, dimension(3,mx1) :: ncll, nclr
         end subroutine
      end interface
!
      interface
         subroutine PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,&
     &mcll,mclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer :: idimp, nppmx, mx1, myp1, npbmx, ntmax, nbmax, irc
         real, dimension(idimp,nppmx,mx1*myp1) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1) :: ppbuff
         real, dimension(idimp,nbmax) :: rbufl, rbufr
         integer, dimension(mx1*myp1) :: kpic
         integer, dimension(8,mx1*myp1) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1) :: ihole
         integer, dimension(3,mx1) :: mcll, mclr
         end subroutine
      end interface
!
      interface
         subroutine PPCGUARD2XL(fxy,nyp,nx,ndim,nxe,nypmx)
         implicit none
         integer :: nyp, nx, ndim, nxe, nypmx
         real, dimension(ndim,nxe,nypmx) :: fxy
         end subroutine
      end interface
!
      interface
         subroutine PPAGUARD2XL(q,nyp,nx,nxe,nypmx)
         implicit none
         integer :: nyp, nx, nxe, nypmx
         real, dimension(nxe,nypmx) :: q
         end subroutine
      end interface
!
      interface
         subroutine MPPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt, &
     &nyv,kxp,nyhd)
         implicit none
         integer :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real :: ax, ay, affp, we
         complex, dimension(nyv,kxp) :: q
         complex, dimension(2,nyv,kxp) :: fxy
         complex, dimension(nyhd,kxp) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine WPFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer :: indx, indy, nxhyd, nxyhd
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WPPFFT2RM(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
         integer :: kxp, kyp, kypd, nxhyd, nxyhd
         real :: ttp
         real, dimension(2*nxvh,kypd) :: f
         complex, dimension(nyv,kxp) :: g
         complex, dimension(kxp,kyp) :: bs, br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WPPFFT2RM2(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx&
     &,indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
         integer :: kxp, kyp, kypd, nxhyd, nxyhd
         real :: ttp
         real, dimension(2,2*nxvh,kypd) :: f
         complex, dimension(2,nyv,kxp) :: g
         complex, dimension(2,kxp,kyp) :: bs, br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RMXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,  &
     &kypp,nxvh,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, nxvh, kypi, kypp, kypd
         integer :: nxhyd, nxyhd
         real, dimension(2*nxvh,kypd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RMXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,  &
     &kxpp,nyv,kxp,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxp
         integer :: nxhyd, nxyhd
         complex, dimension(nyv,kxp) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RM2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, nxvh, kypi, kypp, kypd
         integer :: nxhyd, nxyhd
         real, dimension(3,2*nxvh,kypd) :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
!
      interface
         subroutine PPFFT2RM2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxp,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxp
         integer :: nxhyd, nxyhd
         complex, dimension(3,nyv,kxp) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
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