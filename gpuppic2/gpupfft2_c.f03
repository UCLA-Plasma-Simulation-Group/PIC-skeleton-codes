!-----------------------------------------------------------------------
! Fortran2003 interface to CUDA FFT Library
! written by Viktor K. Decyk, UCLA
      module gpupfft2_c
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
         subroutine gpupfft2rrcux(g_f,g_bsm,isign,indx,indy,kstrt,nvp,  &
     &kxp1,kyp,nxh1d,kypd) bind(C,name='gpupfft2rrcux')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, kstrt, nvp, kxp1
         integer(c_int), value :: kyp, nxh1d, kypd
         type (c_ptr), value :: g_f, g_bsm
         end subroutine
      end interface
!
      interface
         subroutine gpupfft2rrcuy(g_g,g_brm,isign,indx,indy,kstrt,nvp,  &
     &kxp1,kyp,nyd) bind(C,name='gpupfft2rrcuy')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, kstrt, nvp, kxp1
         integer(c_int), value :: kyp, nyd
         type (c_ptr), value :: g_g, g_brm
         end subroutine
      end interface
!
      interface
         subroutine gpupfft2rrcuxn(g_fn,g_bsm,isign,indx,indy,ndim,kstrt&
     &,nvp,kxp1,kyp,nxh1d,kypd) bind(C,name='gpupfft2rrcuxn')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, ndim, kstrt, nvp
         integer(c_int), value :: kxp1, kyp, nxh1d, kypd
         type (c_ptr), value :: g_fn, g_bsm
         end subroutine
      end interface
!
      interface
         subroutine gpupfft2rrcuyn(g_gn,g_brm,isign,indx,indy,ndim,kstrt&
     &,nvp,kxp1,kyp,nyd) bind(C,name='gpupfft2rrcuyn')
         use iso_c_binding
         implicit none
         integer(c_int), value :: isign, indx, indy, ndim, kstrt, nvp
         integer(c_int), value :: kxp1, kyp, nyd
         type (c_ptr), value :: g_gn, g_brm
         end subroutine
      end interface
!
      interface
         subroutine cgpuppsltpose(g_f,g_g,ani,nx,ny,kxp,kyp,kstrt,nxv,  &
     &nyv) bind(C,name='cgpuppsltpose')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, kxp, kyp, kstrt, nxv, nyv
         real(c_float), value :: ani
         type (c_ptr), value :: g_f, g_g
         end subroutine
      end interface
!
      interface
         subroutine cgpuppsltposen(g_fn,g_gn,ani,nx,ny,kxp,kyp,kstrt,   &
     &ndim,nxv,nyv) bind(C,name='cgpuppsltposen')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nx, ny, kxp, kyp, kstrt, ndim, nxv
         integer(c_int), value :: nyv
         real(c_float), value :: ani
         type (c_ptr), value :: g_fn, g_gn
         end subroutine
      end interface
!
      end module


