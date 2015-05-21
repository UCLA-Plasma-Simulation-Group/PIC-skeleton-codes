!-----------------------------------------------------------------------
! Interface file for ppush2.f
      module ppush2_h
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
     &nx, ny,idimp,npmax,idps,ipbc,ierr)
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
         integer :: nppmx, idimp, npmax, mx, my, mx1, mxyp1
         integer :: npp, noff
         real, dimension(idimp,npmax) :: part
         integer, dimension(mxyp1) :: kpic
         integer, dimension(1) :: irc
         end subroutine
      end interface
!
      interface
         subroutine PPPMOVIN2LT(part,ppart,kpic,npp,noff,nppmx,idimp,   &
     &npmax,mx,my,mx1,mxyp1,irc)
         implicit none
         integer :: nppmx, idimp, npmax, mx, my, mx1, mxyp1
         integer :: npp, noff
         real, dimension(idimp,npmax) :: part
         real, dimension(nppmx,idimp,mxyp1) :: ppart
         integer, dimension(mxyp1) :: kpic
         integer, dimension(1) :: irc
         end subroutine
      end interface
!
      interface
         subroutine PPPCHECK2LT(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my&
     &,mx1,myp1,irc)
         implicit none
         integer :: noff, nyp, idimp, nppmx, nx, mx, my, mx1, myp1
         real, dimension(nppmx,idimp,mx1*myp1) :: ppart
         integer, dimension(mx1*myp1) :: kpic
         integer, dimension(1) :: irc
         end subroutine
      end interface
!
      interface
         subroutine PPOIS22T(qt,fxyt,isign,ffct,ax,ay,affp,we,nx,ny,    &
     &kstrt, nyv,kxp1,nyhd)
         implicit none
         integer :: isign, nx, ny, kstrt, nyv, kxp1, nyhd
         real :: ax, ay, affp, we
         complex, dimension(nyv,kxp1) :: qt
         complex, dimension(nyv,2,kxp1) :: fxyt
         complex, dimension(nyhd,kxp1) :: ffct
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
         function ranorm()
         implicit none
         double precision :: ranorm
         end function
      end interface
!
      end module