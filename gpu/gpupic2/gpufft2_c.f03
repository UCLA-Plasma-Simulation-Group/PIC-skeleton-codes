!-----------------------------------------------------------------------
! Fortran2003 interface to CUDA FFT Library
! written by Viktor K. Decyk, UCLA
      module gpufft2_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine gpufft2rrcuinit(nx,ny,ndim)                         &
     &bind(C,name='gpufft2rrcuinit')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, ndim
         end subroutine
      end interface
!
      interface
         subroutine gpufft2cuinit(nx,ny,ndim)                           &
     &bind(C,name='gpufft2cuinit')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, ndim
         end subroutine
      end interface
!
      interface
         subroutine gpufft2rrcudel() bind(C,name='gpufft2rrcudel')
         implicit none
         end subroutine
      end interface
!
      interface
         subroutine gpufft2cudel() bind(C,name='gpufft2cudel')
         implicit none
         end subroutine
      end interface
!
      interface
         subroutine gpufft2rrcu(g_f,g_g,isign,indx,indy,nxh1d,nyd)      &
     &bind(C,name='gpufft2rrcu')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, nxh1d, nyd
         type (c_ptr), value :: g_f, g_g
         end subroutine
      end interface
!
      interface
         subroutine gpufft2rrcun(g_fn,g_gn,isign,indx,indy,ndim,nxh1d,  &
     &nyd) bind(C,name='gpufft2rrcun')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, ndim, nxh1d, nyd
         type (c_ptr), value :: g_fn, g_gn
         end subroutine
      end interface
!
      end module


