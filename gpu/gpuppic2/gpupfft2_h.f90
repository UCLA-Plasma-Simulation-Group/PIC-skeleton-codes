!-----------------------------------------------------------------------
! Fortran interface to CUDA FFT Library
! written by Viktor K. Decyk, UCLA
      module gpupfft2_h
      implicit none
!
      interface
         subroutine gpufft2rrcuinit(nx,ny,ndim)
         implicit none
         integer :: nx, ny, ndim
         end subroutine
      end interface
!
      interface
         subroutine gpufft2cuinit(nx,ny,ndim)
         implicit none
         integer :: nx, ny, ndim
         end subroutine
      end interface
!
      interface
         subroutine gpufft2rrcudel()
         implicit none
         end subroutine
      end interface
!
      interface
         subroutine gpufft2cudel()
         implicit none
         end subroutine
      end interface
!
      interface
         subroutine gpupfft2rrcux(gp_f,gp_bsm,isign,indx,indy,kstrt,nvp,&
     &kxp1,kyp,nxh1d,kypd)
         implicit none
         integer :: isign, indx, indy, kstrt, nvp, kxp1, kyp
         integer :: nxh1d, kypd
         integer, dimension(2) :: gp_f, gp_bsm
         end subroutine
      end interface
!
      interface
         subroutine gpupfft2rrcuy(gp_g,gp_brm,isign,indx,indy,kstrt,nvp,&
     &kxp1,kyp,nyd)
         implicit none
         integer :: isign, indx, indy, kstrt, nvp, kxp1, kyp, nyd
         integer, dimension(2) :: gp_g, gp_brm
         end subroutine
      end interface
!
      interface
         subroutine gpupfft2rrcuxn(gp_fn,gp_bsm,isign,indx,indy,ndim,   &
     &kstrt,nvp,kxp1,kyp,nxh1d,kypd)
         implicit none
         integer :: isign, indx, indy, ndim, kstrt, nvp, kxp1, kyp
         integer :: nxh1d, kypd
         integer, dimension(2) :: gp_fn, gp_bsm
         end subroutine
      end interface
!
      interface
         subroutine gpupfft2rrcuyn(gp_gn,gp_brm,isign,indx,indy,ndim,   &
     &kstrt,nvp,kxp1,kyp,nyd)
         implicit none
         integer :: isign, indx, indy, ndim, kstrt, nvp, kxp1, kyp, nyd
         integer, dimension(2) :: gp_gn, gp_brm
         end subroutine
      end interface
!
      interface
         subroutine cgpuppsltpose(gp_f,gp_g,ani,nx,ny,kxp,kyp,kstrt,nxv,&
     &nyv)
         implicit none
         integer :: nx, ny, kxp, kyp, kstrt, nxv, nyv
         real :: ani
         integer, dimension(2) :: gp_f, gp_g
         end subroutine
      end interface
!
      interface
         subroutine cgpuppsltposen(gp_fn,gp_gn,ani,nx,ny,kxp,kyp,kstrt, &
     &ndim,nxv,nyv)
         implicit none
         integer :: nx, ny, kxp, kyp, kstrt, ndim, nxv, nyv
         real :: ani
         integer, dimension(2) :: gp_fn, gp_gn
         end subroutine
      end interface
!
      end module


