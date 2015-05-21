!-----------------------------------------------------------------------
! Interface file for mbpush2.c
!
      interface
         subroutine cdistr2h(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp,&
     &nop,nx,ny,ipbc)
         implicit none
         integer, intent(in) :: npx, npy, idimp, nop, nx, ny, ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
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
         subroutine cgbppush23l(ppart,fxy,bxy,kpic,qbm,dt,dtc,ek,idimp, &
     &nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ipbc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         real, dimension(*), intent(in) :: fxy, bxy
!        integer, dimension(mxy1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cgbppushf23l(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,   &
     &dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         real, dimension(*), intent(in) :: fxy, bxy
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
         subroutine cgrbppush23l(ppart,fxy,bxy,kpic,qbm,dt,dtc,ci,ek,   &
     &idimp, nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxy1), intent(inout)  :: ppart
         real, dimension(*), intent(inout)  :: ppart
!        real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         real, dimension(*), intent(in) :: fxy, bxy
!        integer, dimension(mxy1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cgrbppushf23l(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,  &
     &dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nyv), intent(in) :: fxy, bxy
         real, dimension(*), intent(in) :: fxy, bxy
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
         subroutine cgjppost2l(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny,mx,&
     &my, nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ipbc
         real, intent(in) :: qm, dt
!        real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nyv), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
!        integer, dimension(mxy1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cgjppostf2l(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,    &
     &idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt
!        real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nyv), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
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
         subroutine cgrjppost2l(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx,ny&
     &,mx,my,nxv,nyv,mx1,mxy1,ipbc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ipbc
         real, intent(in) :: qm, dt, ci
!        real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nyv), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
!        integer, dimension(mxy1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine cgrjppostf2l(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx,&
     &idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nx, ny, mx, my, nxv, nyv
         integer, intent(in) :: mx1, mxy1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qm, dt, ci
!        real, dimension(idimp,nppmx,mxy1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nyv), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
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
         subroutine cbguard2l(bxy,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
!        real, dimension(3,nxe,nye), intent(inout) :: bxy
         real, dimension(*), intent(inout) :: bxy
         end subroutine
      end interface
!
      interface
         subroutine cacguard2l(cu,nx,ny,nxe,nye)
         implicit none
         integer, intent(in) :: nx, ny, nxe, nye
!        real, dimension(3,nxe,nye), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
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
         subroutine cmpois23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,  &
     &nyv,nxhd,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
!        real, dimension(2*nxvh,nyv), intent(in) :: q
         real, dimension(*), intent(in) :: q
!        real, dimension(3,2*nxvh,nyv), intent(inout) :: fxy
         real, dimension(*), intent(inout) :: fxy
!        complex, dimension(nxhd,nyhd), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cmcuperp2(cu,nx,ny,nxvh,nyv)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv
!        real, dimension(3,2*nxvh,nyv), intent(inout) :: cu
         real, dimension(*), intent(inout) :: cu
         end subroutine
      end interface
!
      interface
         subroutine cmibpois23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,    &
     &nyhd)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
!        real, dimension(3,2*nxvh,nyv), intent(in) :: cu
         real, dimension(*), intent(in) :: cu
!        complex, dimension(3,nxvh,nyv), intent(inout) :: bxy
         complex, dimension(*), intent(inout) :: bxy
!        complex, dimension(nxhd,nyhd), intent(in) :: ffc
         complex, dimension(*), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cmmaxwel2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,&
     &nxhd,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, nxvh, nyv, nxhd, nyhd
         real, intent(in) :: ci, dt
         real, intent(inout) :: wf, wm
!        complex, dimension(3,nxvh,nyv), intent(inout) :: exy, bxy
         complex, dimension(*), intent(inout) :: exy, bxy
!        real, dimension(3,2*nxvh,nyv), intent(in) :: cu
         real, dimension(*), intent(in) :: cu
!        complex, dimension(nxhd,nyhd), intent(in) :: ffc
         complex, dimension(*), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine cmemfield2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,   &
     &nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
!        real, dimension(3,2*nxvh,nyv), intent(inout) :: fxy
         real, dimension(*), intent(inout) :: fxy
!        complex, dimension(3,nxvh,nyv), intent(in) :: exy
         complex, dimension(*), intent(in) :: exy
!        complex, dimension(nxhd,nyhd), intent(in) :: ffc
         complex, dimension(*), intent(in) :: ffc
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
         subroutine cwfft2rm3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd&
     &,nxyhd)
         implicit none
         integer , intent(in):: isign, indx, indy, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
!        real, dimension(3,2*nxhd,nyd), intent(inout) :: f
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
         integer, intent(in) :: isign, indx, indy, nyi, nyp, nxhd, nyd
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
         integer, intent(in) :: isign, indx, indy, nxi, nxp, nxhd, nyd
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
         subroutine cfft2rm3x(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd, &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nyi, nyp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
!        real, dimension(3,2*nxhd,nyd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine cfft2rm3y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd, &
     &nyd,nxhyd,nxyhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, nxi, nxp, nxhd, nyd
         integer, intent(in) :: nxhyd, nxyhd
!        real, dimension(3,2*nxhd,nyd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      end
