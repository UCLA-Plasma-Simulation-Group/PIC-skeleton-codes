!-----------------------------------------------------------------------
! Interface file for push2.f
      module push2_h
      implicit none
!
      interface
         subroutine DISTR2(part,vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx,ny,&
     &ipbc)
         implicit none
         integer :: npx, npy, idimp, nop, nx, ny, ipbc
         real :: vtx, vty, vdx, vdy
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
!
      interface
         subroutine DBLKP2L(part,kpic,nppmx,idimp,nop,mx,my,mx1,mxy1,irc&
     &)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mx1, mxy1
         real, dimension(idimp,nop) :: part
         integer, dimension(mxy1) :: kpic
         integer, dimension(1) :: irc
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN2LT(part,ppart,kpic,nppmx,idimp,nop,mx,my,mx1&
     &,mxy1,irc)
         implicit none
         integer :: nppmx, idimp, nop, mx, my, mx1, mxy1
         real, dimension(idimp,nop) :: part
         real, dimension(nppmx,idimp,mxy1) :: ppart
         integer, dimension(mxy1) :: kpic
         integer, dimension(1) :: irc
         end subroutine
      end interface
!
      interface
         subroutine PPCHECK2LT(ppart,kpic,idimp,nppmx,nx,ny,mx,my,mx1,  &
     &my1,irc)
         implicit none
         integer :: idimp, nppmx, nx, ny, mx, my, mx1, my1
         real, dimension(nppmx,idimp,mx1*my1) :: ppart
         integer, dimension(mx1*my1) :: kpic
         integer, dimension(1) :: irc
         end subroutine
      end interface
!
      interface
         subroutine POIS22T(qt,fxyt,isign,ffct,ax,ay,affp,we,nx,ny,nxvh,&
     &nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ax, ay, affp, we
         complex, dimension(nyv,nxvh) :: qt
         complex, dimension(nyv,2,nxvh) :: fxyt
         complex, dimension(nyhd,nxhd) :: ffct
         end subroutine
      end interface
!
      interface
         subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
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