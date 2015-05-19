!-----------------------------------------------------------------------
! Interface file for mpbpush2.c
!
      interface
         subroutine cpdicomp2l(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp, &
     &idps)
         implicit none
         integer, intent(in) :: ny, kstrt, nvp, idps
         integer, intent(inout) :: nyp, noff, nypmx, nypmn
!        real, dimension(idps), intent(inout) :: edges
         real, dimension(*), intent(inout) :: edges
         end subroutine
      end interface
!
      interface
         subroutine cpdistr2h(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz&
     &,npx,npy,nx,ny,idimp,npmax,idps,ipbc,ierr)
         implicit none
         integer, intent(in) :: nps, npx, npy, nx, ny, idimp, npmax
         integer, intent(in) :: idps, ipbc
         integer, intent(inout) :: npp, ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
!        real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(idps), intent(in) :: edges
         real, dimension(*), intent(in) :: edges
         end subroutine
      end interface
!
      interface
         subroutine cppdblkp2l(part,kpic,npp,noff,nppmx,idimp,npmax,mx, &
     &my,mx1,mxyp1,irc)
         implicit none
         integer, intent(in) :: idimp, npmax, mx, my, mx1, mxyp1, npp
         integer, intent(in) :: noff
         integer, intent(inout) :: nppmx, irc
!        real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        integer, dimension(mxyp1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cpppmovin2l(part,ppart,kpic,npp,noff,nppmx,idimp,   &
     &npmax,mx,my,mx1,mxyp1,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, npmax, mx, my, mx1, mxyp1
         integer, intent(in) :: npp, noff
         integer, intent(inout) :: irc
!        real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        integer, dimension(mxyp1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cpppcheck2l(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my&
     &,mx1,myp1,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, my, mx1, myp1
         integer, intent(in) :: noff, nyp
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1*myp1), intent(in) :: ppart
         real, dimension(*), intent(in) :: ppart
!        integer, dimension(mx1*myp1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cppgbppush23l(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc&
     &,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(*), intent(in) :: fxy, bxy
!        integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cppgbppushf23l(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp&
     &,qbm,dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax, &
     &irc)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxyp1), intent(inout)  :: ppart
         real, dimension(*), intent(inout)  :: ppart
!        real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(*), intent(in) :: fxy, bxy
!        integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
!        integer, dimension(8,mxyp1), intent(inout)  :: ncl
         integer, dimension(*), intent(inout)  :: ncl
!        integer, dimension(2,ntmax+1,mxyp1), intent(inout)  :: ihole
         integer, dimension(*), intent(inout)  :: ihole
         end subroutine
      end interface
!
      interface
         subroutine cppgrbppush23l(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,  &
     &dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(*), intent(in) :: fxy, bxy
!        integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cppgrbppushf23l(ppart,fxy,bxy,kpic,ncl,ihole,noff,  &
     &nyp,qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1, &
     &ntmax,irc)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ntmax, irc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         real, dimension(*), intent(in) :: fxy, bxy
!        integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
!        integer, dimension(8,mxyp1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mxyp1), intent(inout) :: ihole
         integer, dimension(*), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine cppgppost2l(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my, &
     &nxv, nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: noff, idimp, nppmx, mx, my, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1
         real, intent(in) :: qm
!        real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(*), intent(in) :: ppart
!        real, dimension(nxv,nypmx), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
!        integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cppgjppost2l(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,  &
     &nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: qm, dt
!        real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nypmx), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
!        integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cppgjppostf2l(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt&
     &,nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
         implicit none
         integer, intent(in) :: noff, nyp, nppmx, idimp, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ntmax, irc
         real, intent(in) :: qm, dt
!        real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nypmx), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
!        integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
!        integer, dimension(8,mxyp1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mxyp1), intent(inout) :: ihole
         integer, dimension(*), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine cppgrjppost2l(ppart,cu,kpic,noff,qm,dt,ci,nppmx,    &
     &idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
         implicit none
         integer, intent(in) :: noff, nppmx, idimp, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: qm, dt, ci
!        real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nypmx), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
!        integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cppgrjppostf2l(ppart,cu,kpic,ncl,ihole,noff,nyp,qm, &
     &dt,ci,nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
         implicit none
         integer, intent(in) :: noff, nyp, nppmx, idimp, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ntmax, irc
         real, intent(in) :: qm, dt, ci
!        real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nypmx), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
!        integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
!        integer, dimension(8,mxyp1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mxyp1), intent(inout) :: ihole
         integer, dimension(*), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine cppporder2la(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole&
     &,ncll,nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax, &
     &nbmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, mx1, myp1
         integer, intent(in) :: npbmx, ntmax, nbmax, noff, nyp
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(idimp,npbmx,mx1*myp1), intent(inout) :: ppbuff
         real, dimension(*), intent(inout) :: ppbuff
!        real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr
         real, dimension(*), intent(inout) :: sbufl, sbufr
!        integer, dimension(mx1*myp1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
!        integer, dimension(8,mx1*myp1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mx1*myp1), intent(inout) :: ihole
         integer, dimension(*), intent(inout) :: ihole
!        integer, dimension(3,mx1), intent(inout) :: ncll, nclr
         integer, dimension(*), intent(inout) :: ncll, nclr
         end subroutine
      end interface
!
      interface
         subroutine cppporderf2la(ppart,ppbuff,sbufl,sbufr,ncl,ihole,   &
     &ncll,nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, myp1, npbmx, ntmax
         integer, intent(in) :: nbmax
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(idimp,npbmx,mx1*myp1), intent(inout) :: ppbuff
         real, dimension(*), intent(inout) :: ppbuff
!        real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr
         real, dimension(*), intent(inout) :: sbufl, sbufr
!        integer, dimension(8,mx1*myp1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mx1*myp1), intent(in) :: ihole
         integer, dimension(*), intent(in) :: ihole
!        integer, dimension(3,mx1), intent(inout) :: ncll, nclr
         integer, dimension(*), intent(inout) :: ncll, nclr
         end subroutine
      end interface
!
      interface
         subroutine cppporder2lb(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole&
     &,mcll,mclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, myp1, npbmx, ntmax
         integer, intent(in) :: nbmax
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(idimp,npbmx,mx1*myp1), intent(in) :: ppbuff
         real, dimension(*), intent(in) :: ppbuff
!        real, dimension(idimp,nbmax), intent(in) :: rbufl, rbufr
         real, dimension(*), intent(in) :: rbufl, rbufr
!        integer, dimension(mx1*myp1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
!        integer, dimension(8,mx1*myp1), intent(in) :: ncl
         integer, dimension(*), intent(in) :: ncl
!        integer, dimension(2,ntmax+1,mx1*myp1), intent(in) :: ihole
         integer, dimension(*), intent(in) :: ihole
!        integer, dimension(3,mx1), intent(in) :: mcll, mclr
         integer, dimension(*), intent(in) :: mcll, mclr
         end subroutine
      end interface
!
      interface
         subroutine cppcguard2xl(fxy,nyp,nx,ndim,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, ndim, nxe, nypmx
!        real, dimension(ndim,nxe,nypmx), intent(inout) :: fxy
         real, dimension(*), intent(inout) :: fxy
         end subroutine
      end interface
!
      interface
         subroutine cppaguard2xl(q,nyp,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
!        real, dimension(nxe,nypmx), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine cppacguard2xl(cu,nyp,nx,ndim,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, ndim, nxe, nypmx
!        real, dimension(ndim,nxe,nypmx), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine cmppois23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,&
     &nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
!        complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(*), intent(in) :: q
!        complex, dimension(3,nyv,kxp), intent(inout) :: fxy
         complex, dimension(*), intent(inout) :: fxy
!        complex, dimension(nyhd,kxp), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cmppcuperp2(cu,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
!        complex, dimension(3,nyv,kxp), intent(inout) :: cu
         complex, dimension(*), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine cmippbpoisp23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp, &
     &nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
!        complex, dimension(3,nyv,kxp), intent(in) :: cu
         complex, dimension(*), intent(in) :: cu
!        complex, dimension(3,nyv,kxp), intent(inout) :: bxy
         complex, dimension(*), intent(inout) :: bxy
!        complex, dimension(nyhd,kxp), intent(in) :: ffc
         complex, dimension(*), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cmppmaxwel2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,  &
     &kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: affp, ci, dt
         real, intent(inout) :: wf, wm
!        complex, dimension(3,nyv,kxp), intent(inout) :: exy, bxy
         complex, dimension(*), intent(inout) :: exy, bxy
!        complex, dimension(3,nyv,kxp), intent(in)  :: cu
         complex, dimension(*), intent(in)  :: cu
!        complex, dimension(nyhd,kxp), intent(in) :: ffc
         complex, dimension(*), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cmppemfield2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp, &
     &nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
!        complex, dimension(3,nyv,kxp), intent(inout) :: fxy
         complex, dimension(*), intent(inout) :: fxy
!        complex, dimension(3,nyv,kxp), intent(in) :: exy
         complex, dimension(*), intent(in) :: exy
!        complex, dimension(nyhd,kxp), intent(in) :: ffc
         complex, dimension(*), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cwpfft2rinit(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: indx, indy, nxhyd, nxyhd
!        integer, dimension(nxhyd), intent(inout) :: mixup
         integer, dimension(*), intent(inout) :: mixup
!        complex, dimension(nxyhd), intent(inout) :: sct
         complex, dimension(*), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cwppfft2rm(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx&
     &,indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp, kyp, kypd, nxhyd, nxyhd
         real, intent(inout) :: ttp
!        real, dimension(2*nxvh,kypd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        complex, dimension(nyv,kxp), intent(inout) :: g
         complex, dimension(*), intent(inout) :: g
!        complex, dimension(kxp,kyp), intent(inout) :: bs, br
         complex, dimension(*), intent(inout) :: bs, br
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cwppfft2rm3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,   &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp, kyp, kypd, nxhyd, nxyhd
         real, intent(inout) :: ttp
!        real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        complex, dimension(3,nyv,kxp), intent(inout) :: g
         complex, dimension(*), intent(inout) :: g
!        complex, dimension(3,kxp,kyp), intent(inout) :: bs, br
         complex, dimension(*), intent(inout) :: bs, br
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cppfft2rmxx(f,isign,mixup,sct,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kypi, kypp
         integer, intent(in) :: nxvh, kypd, nxhyd, nxyhd
!        real, dimension(2*nxvh,kypd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cppfft2rmxy(g,isign,mixup,sct,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxp,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kxpi, kxpp
         integer, intent(in) :: nyv, kxp, nxhyd, nxyhd
!        complex, dimension(nyv,kxp), intent(inout) :: g
         complex, dimension(*), intent(inout) :: g
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cppfft2rm3xx(f,isign,mixup,sct,indx,indy,kstrt,kypi,&
     &kypp,nxvh,kypd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kypi, kypp
         integer, intent(in) :: nxvh, kypd, nxhyd, nxyhd
!        real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cppfft2rm3xy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,&
     &kxpp,nyv,kxp,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kxpi, kxpp
         integer, intent(in) :: nyv, kxp, nxhyd, nxyhd
!        complex, dimension(3,nyv,kxp), intent(inout) :: g
         complex, dimension(*), intent(inout) :: g
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
!-----------------------------------------------------------------------
! Interface file for mpplib2.c
!
      interface
         subroutine cppinit2(idproc,nvp,argc,argv)
         implicit none
         integer, intent(inout) :: idproc, nvp
         integer, intent(in) :: argc
         character, dimension(*), intent(in) :: argv
         end subroutine
      end interface
!
      interface
         subroutine cppexit
         implicit none
         end subroutine
      end interface
!
      interface
         subroutine cppabort
         implicit none
         end subroutine
      end interface
!
      interface
         subroutine cpwtimera(icntrl,time,dtime)
         implicit none
         integer, intent(in) :: icntrl
         real, intent(inout) :: time
         double precision, intent(inout) :: dtime
         end subroutine
      end interface
!
      interface
         subroutine cppsum(f,g,nxp)
         implicit none
         integer, intent(in) :: nxp
!        real, dimension(nxp), intent(inout) :: f, g
         real, dimension(*), intent(inout) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine cppdsum(f,g,nxp)
         implicit none
         integer, intent(in) :: nxp
!        double precision, dimension(nxp), intent(inout) :: f, g
         double precision, dimension(*), intent(inout) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine cppimax(if,ig,nxp)
         implicit none
         integer, intent(in) :: nxp
!        integer, dimension(nxp), intent(inout) :: if, ig
         integer, dimension(*), intent(inout) :: if, ig
         end subroutine
      end interface
!
      interface
         subroutine cppdmax(f,g,nxp)
         implicit none
         integer, intent(in) :: nxp
!        double precision, dimension(nxp), intent(inout) :: f, g
         double precision, dimension(*), intent(inout) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine cppncguard2l(f,nyp,kstrt,nvp,nxv,nypmx)
         implicit none
         integer, intent(in) :: nyp, kstrt, nvp, nxv, nypmx
!        real, dimension(nxv,nypmx), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine cppnaguard2l(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
         implicit none
         integer, intent(in) :: nyp, kstrt, nvp, nx, nxv, nypmx
!        real, dimension(nxv,nypmx), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        real, dimension(nxv), intent(inout) :: scr
         real, dimension(*), intent(inout) :: scr
         end subroutine
      end interface
!
      interface
         subroutine cppnacguard2l(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
         implicit none
         integer, intent(in) :: ndim, nyp, kstrt, nvp, nx, nxv, nypmx
!        real, dimension(ndim,nxv,nypmx), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        real, dimension(ndim,nxv), intent(inout) :: scr
         real, dimension(*), intent(inout) :: scr
         end subroutine
      end interface
!
      interface
         subroutine cpptpose(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,   &
     &kxpd,kypd)
         implicit none
         integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv
         integer, intent(in) :: kxpd, kypd
!        real, dimension(2*nxv,kypd), intent(in) :: f
         real, dimension(*), intent(in) :: f
!        complex, dimension(nyv,kxpd), intent(inout) :: g
         complex, dimension(*), intent(inout) :: g
!        complex, dimension(kxp*kyp), intent(inout) :: s, t
         complex, dimension(*), intent(inout) :: s, t
         end subroutine
      end interface
!
      interface
         subroutine cppntpose(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv, &
     &nyv,kxpd,kypd)
         implicit none
         integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
         integer, intent(in) :: nxv, nyv, kxpd, kypd
!        real, dimension(ndim,2*nxv,kypd), intent(in) :: f
         real, dimension(*), intent(in) :: f
!        complex, dimension(ndim,nyv,kxpd), intent(inout) :: g
         complex, dimension(*), intent(inout) :: g
!        complex, dimension(ndim,kxp*kyp), intent(inout) :: s, t
         complex, dimension(*), intent(inout) :: s, t
         end subroutine
      end interface
!
      interface
         subroutine cpppmove2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,   &
     &mclr,kstrt,nvp,idimp,nbmax,mx1)
         implicit none
         integer, intent(in) :: kstrt, nvp, idimp, nbmax, mx1
!        real, dimension(idimp,nbmax), intent(in) :: sbufr, sbufl
         real, dimension(*), intent(in) :: sbufr, sbufl
!        real, dimension(idimp,nbmax), intent(inout) :: rbufr, rbufl
         real, dimension(*), intent(inout) :: rbufr, rbufl
!        integer, dimension(3,mx1), intent(in) :: ncll, nclr
         integer, dimension(*), intent(in) :: ncll, nclr
!        integer, dimension(3,mx1), intent(inout) :: mcll, mclr
         integer, dimension(*), intent(inout) :: mcll, mclr
         end subroutine
      end interface
!
      end
