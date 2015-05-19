!-----------------------------------------------------------------------
! Interface file for mpush3.f
!     module mpush3_h
!     implicit none
!
      interface
         subroutine DISTR3(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,    &
     &idimp,nop,nx,ny,nz,ipbc)
         implicit none
         integer, intent(in) :: npx, npy, npz, idimp, nop, nx, ny, nz
         integer, intent(in) :: ipbc
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
!        real, dimension(idimp,nop), intent(inout) :: part
         real, dimension(*), intent(inout) :: part
         end subroutine
      end interface
!
      interface
         subroutine DBLKP3L(part,kpic,nppmx,idimp,nop,mx,my,mz,mx1,my1, &
     &mxyz1,irc)
         implicit none
         integer, intent(in) :: idimp, nop, mx, my, mz, mx1, my1, mxyz1
         integer, intent(inout) :: nppmx, irc
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        integer, dimension(mxyz1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPMOVIN3L(part,ppart,kpic,nppmx,idimp,nop,mx,my,mz, &
     &mx1,my1,mxyz1,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, nop, mx, my, mz
         integer, intent(in) :: mx1, my1, mxyz1
         integer, intent(inout) :: irc
!        real, dimension(idimp,nop), intent(in) :: part
         real, dimension(*), intent(in) :: part
!        real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        integer, dimension(mxyz1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPCHECK3L(ppart,kpic,idimp,nppmx,nx,ny,nz,mx,my,mz, &
     &mx1,my1,mz1,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: mx1, my1, mz1
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1*my1*mz1), intent(in) :: ppart
         real, dimension(*), intent(in) :: ppart
!        integer, dimension(mx1*my1*mz1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GPPUSH3L(ppart,fxyz,kpic,qbm,dt,ek,idimp,nppmx,nx,ny&
     &,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ipbc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz
         real, dimension(*), intent(in) :: fxyz
!        integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine GPPUSHF3L(ppart,fxyz,kpic,ncl,ihole,qbm,dt,ek,idimp,&
     &nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: nxv, nyv, nzv, mx1, my1, mxyz1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt
         real, intent(inout) :: ek
!        real, dimension(idimp,nppmx,mxyz1), intent(inout) :: ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(3,nxv,nyv,nzv), intent(in) :: fxyz
         real, dimension(*), intent(in) :: fxyz
!        integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
!        integer, dimension(26,mxyz1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mxyz1), intent(inout) :: ihole
         integer, dimension(*), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine GPPOST3L(ppart,q,kpic,qm,nppmx,idimp,mx,my,mz,nxv,  &
     &nyv,nzv,mx1,my1,mxyz1)
         implicit none
         integer, intent(in) :: nppmx, idimp, mx, my, mz, nxv, nyv, nzv
         integer, intent(in) :: mx1, my1, mxyz1
         real, intent(in) :: qm
!        real, dimension(idimp,nppmx,mxyz1), intent(in) :: ppart
         real, dimension(*), intent(in) :: ppart
!        real, dimension(nxv,nyv,nzv), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
!        integer, dimension(mxyz1), intent(in) :: kpic
         integer, dimension(*), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPORDER3L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx&
     &,ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, nz, mx, my, mz
         integer, intent(in) :: mx1, my1, mz1, npbmx, ntmax
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1*my1*mz1), intent(inout) ::     &
!    &ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(idimp,npbmx,mx1*my1*mz1), intent(inout) ::     &
!    &ppbuff
         real, dimension(*), intent(inout) :: ppbuff
!        integer, dimension(mx1*my1*mz1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
!        integer, dimension(26,mx1*my1*mz1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mx1*my1*mz1), intent(inout) ::    &
!    &ihole
         integer, dimension(*), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPORDERF3L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &mx1,my1,mz1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in):: idimp, nppmx, mx1, my1, mz1, npbmx, ntmax
         integer, intent(inout) :: irc
!        real, dimension(idimp,nppmx,mx1*my1*mz1), intent(inout) ::     &
!    &ppart
         real, dimension(*), intent(inout) :: ppart
!        real, dimension(idimp,npbmx,mx1*my1*mz1), intent(inout) ::     &
!    &ppbuff
         real, dimension(*), intent(inout) :: ppbuff
!        integer, dimension(mx1*my1*mz1), intent(inout) :: kpic
         integer, dimension(*), intent(inout) :: kpic
!        integer, dimension(26,mx1*my1*mz1), intent(inout) :: ncl
         integer, dimension(*), intent(inout) :: ncl
!        integer, dimension(2,ntmax+1,mx1*my1*mz1), intent(in) :: ihole
         integer, dimension(*), intent(in) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine CGUARD3L(fxyz,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
!        real, dimension(3,nxe,nye,nze), intent(inout) :: fxyz
         real, dimension(*), intent(inout) :: fxyz
         end subroutine
      end interface
!
      interface
         subroutine AGUARD3L(q,nx,ny,nz,nxe,nye,nze)
         implicit none
         integer, intent(in) :: nx, ny, nz, nxe, nye, nze
!        real, dimension(nxe,nye,nze), intent(inout) :: q
         real, dimension(*), intent(inout) :: q
         end subroutine
      end interface
!
      interface
         subroutine MPOIS33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz, &
     &nxvh,nyv,nzv,nxhd,nyhd,nzhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, nz, nxvh, nyv, nzv
         integer, intent(in) :: nxhd, nyhd, nzhd
         real, intent(in) :: ax, ay, az, affp
         real, intent(inout) :: we
!        real, dimension(2*nxvh,nyv,nzv), intent(in) :: q
         real, dimension(*), intent(in) :: q
!        real, dimension(3,2*nxvh,nyv,nzv), intent(inout) :: fxyz
         real, dimension(*), intent(inout) :: fxyz
!        complex, dimension(nxhd,nyhd,nzhd), intent(inout) :: ffc
         complex, dimension(*), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: indx, indy, indz, nxhyzd, nxyzhd
!        integer, dimension(nxhyzd), intent(inout) :: mixup
         integer, dimension(*), intent(inout) :: mixup
!        complex, dimension(nxyzhd), intent(inout) :: sct
         complex, dimension(*), intent(inout) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RMX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd, &
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer, intent(in) :: nxhyzd, nxyzhd
!        real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine WFFT3RM3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd, &
     &nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nxhd, nyd, nzd
         integer, intent(in) :: nxhyzd, nxyzhd
!        real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
          subroutine FFT3RMXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
!        real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RMXZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,  &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
!        real, dimension(2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RM3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp, &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nzi, nzp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
!        real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
         complex, dimension(*), intent(in) :: sct
         end subroutine
      end interface
!
      interface
         subroutine FFT3RM3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,  &
     &nxhd,nyd,nzd,nxhyzd,nxyzhd)
         implicit none
         integer, intent(in) :: isign, indx, indy, indz, nyi, nyp
         integer, intent(in) :: nxhd, nyd, nzd, nxhyzd, nxyzhd
!        real, dimension(3,2*nxhd,nyd,nzd), intent(inout) :: f
         real, dimension(*), intent(inout) :: f
!        integer, dimension(nxhyzd), intent(in) :: mixup
         integer, dimension(*), intent(in) :: mixup
!        complex, dimension(nxyzhd), intent(in) :: sct
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
