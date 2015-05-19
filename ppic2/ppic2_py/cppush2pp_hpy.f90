!-----------------------------------------------------------------------
! Interface file for ppush2.c
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
         subroutine cpdistr2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,&
     &nx,ny,idimp,npmax,idps,ipbc,ierr)
         implicit none
         integer, intent(in) :: nps, npx, npy, nx, ny, idimp, npmax
         integer, intent(in) :: idps, ipbc
         integer, intent(inout) :: npp, ierr
         real, intent(in) :: vtx, vty, vdx, vdy
!        real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(idps), intent(in) :: edges
         real, dimension(*), intent(in) :: edges
         end subroutine
      end interface
!
      interface
         subroutine cppgpush2l(part,fxy,edges,npp,noff,ihole,qbm,dt,ek, &
     &nx,ny,idimp,npmax,nxv,nypmx,idps,ntmax,ipbc)
         implicit none
         integer, intent(in) :: npp, noff, nx, ny, idimp, npmax, idps
         integer, intent(in) :: ntmax, nxv, nypmx, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
!        real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(2,nxv,nypmx), intent(in) :: fxy
         real, dimension(*), intent(in) :: fxy
!        real, dimension(idps), intent(in) :: edges
         real, dimension(*), intent(in) :: edges
!        integer, dimension(ntmax+1), intent(inout) :: ihole
         integer, dimension(*), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine cppgpost2l(part,q,npp,noff,qm,idimp,npmax,nxv,nypmx)
         implicit none
         integer, intent(in) :: npp, noff, idimp, npmax, nxv, nypmx
         real, intent(in) :: qm
!        real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        real, dimension(nxv,nypmx), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine cppdsortp2yl(parta,partb,npic,npp,noff,nyp,idimp,   &
     &npmax,nypm1)
         implicit none
         integer, intent(in) :: npp, noff, nyp, idimp, npmax, nypm1
!        real, dimension(idimp,npmax), intent(in) :: parta
         real, dimension(*), intent(in) :: parta
!        real, dimension(idimp,npmax), intent(inout) :: partb
         real, dimension(*), intent(inout) :: partb
!        integer, dimension(nypm1), intent(inout) :: npic
         integer, dimension(*), intent(inout) :: npic
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
         subroutine cppois22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt, &
     &nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
!        complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(*), intent(in) :: q
!        complex, dimension(2,nyv,kxp), intent(inout) :: fxy
         complex, dimension(*), intent(inout) :: fxy
!        complex, dimension(nyhd,kxp), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
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
         subroutine cwppfft2r(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
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
         subroutine cwppfft2r2(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx&
     &,indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp, kyp, kypd, nxhyd, nxyhd
         real, intent(inout) :: ttp
!        real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        complex, dimension(2,nyv,kxp), intent(inout) :: g
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
         subroutine cppfft2rxx(f,isign,mixup,sct,indx,indy,kstrt,kypi,  &
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
         subroutine cppfft2rxy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,  &
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
         subroutine cppfft2r2xx(f,isign,mixup,sct,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, nxvh
         integer, intent(in) :: kypi, kypp, kypd, nxhyd, nxyhd
!        real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cppfft2r2xy(g,isign,mixup,sct,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxp,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt, kxpi, kxpp
         integer, intent(in) :: nyv, kxp, nxhyd, nxyhd
!        complex, dimension(2,nyv,kxp), intent(inout) :: g
         complex, dimension(*), intent(inout) :: g
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
!-----------------------------------------------------------------------
! Interface file for pplib2.c
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
         integer, intent(in) :: nyp, ndim, kstrt, nvp, nx, nxv, nypmx
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
         subroutine cppmove2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,    &
     &ihole,ny,kstrt,nvp,idimp,npmax,idps,nbmax,ntmax,info)
         implicit none
         integer, intent(in) :: ny, kstrt, nvp, idimp, npmax, idps
         integer, intent(in) :: nbmax, ntmax
         integer, intent(inout) :: npp
!        real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
!        real, dimension(idps), intent(in) :: edges
         real, dimension(*), intent(in) :: edges
!        real, dimension(idimp,nbmax), intent(inout) :: sbufr, sbufl
         real, dimension(*), intent(inout) :: sbufr, sbufl
!        real, dimension(idimp,nbmax), intent(inout) :: rbufr, rbufl
         real, dimension(*), intent(inout) :: rbufr, rbufl
!        integer, dimension(ntmax+1) , intent(inout):: ihole
         integer, dimension(*) , intent(inout):: ihole
!        integer, dimension(5), intent(inout) :: info
         integer, dimension(*), intent(inout) :: info
         end subroutine
      end interface
!
      end
