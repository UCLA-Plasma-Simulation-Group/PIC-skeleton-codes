!-----------------------------------------------------------------------
! Interface file for mpush2.c
!
      interface
         subroutine cdistr2(part,vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx,ny&
     &,ipbc)
         implicit none
         integer, intent(in) :: npx, npy, idimp, nop, nx, ny, ipbc
         real, intent(in) :: vtx, vty, vdx, vdy
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine cdblkp2l(part,kpic,nppmx,idimp,nop,mx,my,mx1,mxy1,  &
     &irc)
         implicit none
         integer, intent(in) :: idimp, nop, mx, my, mx1, mxy1
         integer, intent(inout) :: nppmx, irc
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        integer, dimension(mxy1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cppmovin2l(part,ppart,kpic,nppmx,idimp,nop,mx,my,mx1&
     &,mxy1,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nop, mx, my, mx1, mxy1
         integer, intent(inout) :: irc
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        integer, dimension(mxy1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cppcheck2l(ppart,kpic,idimp,nppmx,nx,ny,mx,my,mx1,  &
     &my1,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, mx1, my1
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1*my1), intent(in) :: ppart
         real, dimension(*), intent(in) :: ppart
!        integer, dimension(mx1*my1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cgppush2l(ppart,fxy,kpic,qbm,dt,ek,idimp,nppmx,nx,ny&
     &,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(2,nxv,nyv), intent(in) :: fxy
         real, dimension(*), intent(in) :: fxy
!        integer, dimension(mxy1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cgppushf2l(ppart,fxy,kpic,ncl,ihole,qbm,dt,ek,idimp,&
     &nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(2,nxv,nyv), intent(in) :: fxy
         real, dimension(*), intent(in) :: fxy
!        integer, dimension(mxy1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
!        integer, dimension(8,mxy1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mxy1), intent(inout) :: ihole
         integer, dimension(*), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine cgppost2l(ppart,q,kpic,qm,nppmx,idimp,mx,my,nxv,nyv,&
     &mx1,mxy1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1
         real, intent(in) :: qm
!        real, dimension(idimp,nppmx,mxy1), intent(in) :: ppart
         real, dimension(*), intent(in) :: ppart
!        real, dimension(nxv,nyv), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
!        integer, dimension(mxy1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cpporder2l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &nx,ny,mx,my,mx1,my1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, mx1, my1
         integer, intent(in) :: npbmx, ntmax
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1*my1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(idimp,npbmx,mx1*my1), intent(inout) :: ppbuff
         real, dimension(*), intent(inout) :: ppbuff
!        integer, dimension(mx1*my1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
!        integer, dimension(8,mx1*my1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mx1*my1), intent(inout) :: ihole
         integer, dimension(*), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine cpporderf2l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,&
     &mx1,my1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, my1, npbmx, ntmax
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1*my1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(idimp,npbmx,mx1*my1), intent(inout) :: ppbuff
         real, dimension(*), intent(inout) :: ppbuff
!        integer, dimension(mx1*my1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
!        integer, dimension(8,mx1*my1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mx1*my1), intent(in) :: ihole
         integer, dimension(*), intent(in) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine ccguard2l(fxy,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
!        real, dimension(2,nxe,nye), intent(inout) :: fxy
         real, dimension(*), intent(inout) :: fxy
         end subroutine
      end interface
!
      interface
         subroutine caguard2l(q,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
!        real, dimension(nxe,nye), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine cmpois22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,  &
     &nyv,nxhd,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
!        real, dimension(2*nxvh,nyv), intent(in) :: q
         real, dimension(*), intent(in) :: q
!        real, dimension(2,2*nxvh,nyv), intent(inout) :: fxy
         real, dimension(*), intent(inout) :: fxy
!        complex, dimension(nxhd,nyhd), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cwfft2rinit(mixup,sct,indx,indy,nxhyd,nxyhd)
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
         subroutine cwfft2rmx(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd&
     &,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
!        real, dimension(2*nxhd,nyd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cwfft2rm2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd&
     &,nxyhd)
         implicit none
         integer , intent(in):: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
!        real, dimension(2,2*nxhd,nyd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cfft2rmxx(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd, &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer , intent(in):: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
!        real, dimension(2*nxhd,nyd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cfft2rmxy(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd, &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer , intent(in):: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
!        real, dimension(2*nxhd,nyd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cfft2rm2x(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd, &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
!        real, dimension(2,2*nxhd,nyd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cfft2rm2y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd, &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
!        real, dimension(2,2*nxhd,nyd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      end
