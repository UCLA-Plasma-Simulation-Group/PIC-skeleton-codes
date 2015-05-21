!-----------------------------------------------------------------------
! Interface file for mpush1.f
!     module mpush1_h
!     implicit none
!
      interface
         subroutine DISTR1(part,vtx,vdx,npx,idimp,nop,nx,ipbc)
         implicit none
         integer, intent(in) :: npx, idimp, nop, nx, ipbc
         real, intent(in) :: vtx, vdx
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine DBLKP1L(part,kpic,nppmx,idimp,nop,mx,mx1,irc)
         implicit none
         integer, intent(in) :: idimp, nop, mx, mx1
         integer, intent(inout) :: nppmx, irc
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        integer, dimension(mx1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN1L(part,ppart,kpic,nppmx,idimp,nop,mx,mx1,irc&
     &)
         implicit none
         integer, intent(in) :: nppmx, idimp, nop, mx, mx1
         integer, intent(inout) :: irc
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        integer, dimension(mx1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPCHECK1L(ppart,kpic,idimp,nppmx,nx,mx,mx1,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, mx1
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(*), intent(in) :: ppart
!        integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GPPUSH1L(ppart,fx,kpic,qbm,dt,ek,idimp,nppmx,nx,mx, &
     &nxv,mx1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(nxv), intent(in) :: fx
         real, dimension(*), intent(in) :: fx
!        integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GPPUSHF1L(ppart,fx,kpic,ncl,ihole,qbm,dt,ek,idimp,  &
     &nppmx,nx,mx,nxv,mx1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, nxv, mx1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(nxv), intent(in) :: fx
         real, dimension(*), intent(in) :: fx
!        integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
!        integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         integer, dimension(*), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GPPOST1L(ppart,q,kpic,qm,nppmx,idimp,mx,nxv,mx1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, nxv, mx1
         real, intent(in) :: qm
!        real, dimension(idimp,nppmx,mx1), intent(in) :: ppart
         real, dimension(*), intent(in) :: ppart
!        real, dimension(nxv), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
!        integer, dimension(mx1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPORDER1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx&
     &,mx,mx1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, mx1, npbmx, ntmax
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(idimp,npbmx,mx1), intent(inout) :: ppbuff
         real, dimension(*), intent(inout) :: ppbuff
!        integer, dimension(mx1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
!        integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         integer, dimension(*), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPORDERF1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &mx1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, npbmx, ntmax
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(idimp,npbmx,mx1), intent(inout) :: ppbuff
         real, dimension(*), intent(inout) :: ppbuff
!        integer, dimension(mx1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
!        integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mx1), intent(in) :: ihole
         integer, dimension(*), intent(in) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine CGUARD1L(fx,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
!        real, dimension(nxe), intent(inout) :: fx
         real, dimension(*), intent(inout) :: fx
         end subroutine
      end interface
!
      interface
         subroutine AGUARD1L(q,nx,nxe)
         implicit none
         integer, intent(in) :: nx, nxe
!        real, dimension(nxe), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine POIS1(q,fx,isign,ffc,ax,affp,we,nx)
         implicit none
         integer, intent(in) :: isign, nx
         real, intent(in) :: ax, affp
         real, intent(inout) :: we
!        real, dimension(nx), intent(in) :: q
         real, dimension(*), intent(in) :: q
!        real, dimension(nx), intent(inout) :: fx
         real, dimension(*), intent(inout) :: fx
!        complex, dimension(nx/2), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine WFFT1RINIT(mixup,sct,indx,nxhd)
         implicit none
         integer, intent(in) :: indx, nxhd
!        integer, dimension(nxhd), intent(inout) :: mixup
         integer, dimension(*), intent(inout) :: mixup
!        complex, dimension(nxhd), intent(inout) :: sct
         complex, dimension(*), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer, intent(in) :: isign, indx, nxd, nxhd
!        real, dimension(nxd), intent(inout) :: f, t
         real, dimension(*), intent(inout) :: f, t
!        integer, dimension(nxhd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
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
!     end module
      end
