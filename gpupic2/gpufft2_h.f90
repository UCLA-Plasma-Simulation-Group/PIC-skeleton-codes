!-----------------------------------------------------------------------
! Fortran interface to CUDA FFT Library
! written by Viktor K. Decyk, UCLA
      module gpufft2_h
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
         subroutine gpufft2rrcu(gp_f,gp_g,isign,indx,indy,nxh1d,nyd)
         implicit none
         integer :: isign, indx, indy, nxh1d, nyd
         integer, dimension(2) :: gp_f, gp_g
         end subroutine
      end interface
!
      interface
         subroutine gpufft2rrcun(gp_fn,gp_gn,isign,indx,indy,ndim,nxh1d,&
     &nyd)
         implicit none
         integer :: isign, indx, indy, ndim, nxh1d, nyd
         integer, dimension(2) :: gp_fn, gp_gn
         end subroutine
      end interface
!
      end module


