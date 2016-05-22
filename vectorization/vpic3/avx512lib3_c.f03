!-----------------------------------------------------------------------
! Fortran2003 interface file for avx512lib3.c
      module avx512lib3_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine avx512_fallocate(s_f,nsize,irc)                     &
     &bind(C,name='avx512_fallocate')
         use iso_c_binding
         implicit none
         type (c_ptr) :: s_f
         integer(c_int), value :: nsize
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine avx512_callocate(s_c,nsize,irc)                     &
     &bind(C,name='avx512_callocate')
         use iso_c_binding
         implicit none
         type (c_ptr) :: s_c
         integer(c_int), value :: nsize
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine avx512_iallocate(s_i,nsize,irc)                     &
     &bind(C,name='avx512_iallocate')
         use iso_c_binding
         implicit none
         type (c_ptr) :: s_i
         integer(c_int), value :: nsize
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine avx512_deallocate(s_d,irc)                          &
     &bind(C,name='avx512_deallocate')
         use iso_c_binding
         implicit none
         type (c_ptr), value :: s_d
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine ckncxiscan2(isdata,nths) bind(C,name='ckncxiscan2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nths
         type (c_ptr), value :: isdata
         end subroutine
      end interface
!
      contains
!
!-----------------------------------------------------------------------
      subroutine avx512_f2allocate(s_f,nx,ny,irc)
      implicit none 
      integer, intent(in) :: nx, ny
      integer :: irc
      real, dimension(:,:), pointer :: s_f
! local data
      type (c_ptr) :: fptr
      if (digits(s_f(1,1)) > 24) then
         call avx512_fallocate(fptr,2*nx*ny,irc)
      else
         call avx512_fallocate(fptr,nx*ny,irc)
      endif
      call c_f_pointer(fptr,s_f,(/nx,ny/))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine avx512_f3allocate(s_f,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: nx, ny, nz
      integer :: irc
      real, dimension(:,:,:), pointer :: s_f
! local data
      type (c_ptr) :: fptr
      if (digits(s_f(1,1,1)) > 24) then
         call avx512_fallocate(fptr,2*nx*ny*nz,irc)
      else
         call avx512_fallocate(fptr,nx*ny*nz,irc)
      endif
      call c_f_pointer(fptr,s_f,(/nx,ny,nz/))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine avx512_f4allocate(s_f,ndim,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: ndim, nx, ny, nz
      integer :: irc
      real, dimension(:,:,:,:), pointer :: s_f
! local data
      type (c_ptr) :: fptr
      if (digits(s_f(1,1,1,1)) > 24) then
         call avx512_fallocate(fptr,2*ndim*nx*ny*nz,irc)
      else
         call avx512_fallocate(fptr,ndim*nx*ny*nz,irc)
      endif
      call c_f_pointer(fptr,s_f,(/ndim,nx,ny,nz/))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine avx512_c2allocate(s_f,nx,ny,irc)
      implicit none 
      integer, intent(in) :: nx, ny
      integer :: irc
      complex, dimension(:,:), pointer :: s_f
! local data
      type (c_ptr) :: cptr
      if (digits(real(s_f(1,1))) > 24) then
         call avx512_callocate(cptr,2*nx*ny,irc)
      else
         call avx512_callocate(cptr,nx*ny,irc)
      endif
      call c_f_pointer(cptr,s_f,(/nx,ny/))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine avx512_c3allocate(s_f,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: nx, ny, nz
      integer :: irc
      complex, dimension(:,:,:), pointer :: s_f
! local data
      type (c_ptr) :: cptr
      if (digits(real(s_f(1,1,1))) > 24) then
         call avx512_callocate(cptr,2*nx*ny*nz,irc)
      else
         call avx512_callocate(cptr,nx*ny*nz,irc)
      endif
      call c_f_pointer(cptr,s_f,(/nx,ny,nz/))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine avx512_i1allocate(s_f,nx,irc)
      implicit none 
      integer, intent(in) :: nx
      integer :: irc
      integer, dimension(:), pointer :: s_f
! local data
      type (c_ptr) :: iptr
      if (digits(s_f(1)) > 31) then
         call avx512_iallocate(iptr,2*nx,irc)
      else
         call avx512_iallocate(iptr,nx,irc)
      endif
      call c_f_pointer(iptr,s_f,(/nx/))
      end subroutine
!
      end module
