!-----------------------------------------------------------------------
! Fortran2003 interface file for sselib2.c
      module sselib2_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine sse_fallocate(s_f,nsize,irc)                        &
     &bind(C,name='sse_fallocate')
         use iso_c_binding
         implicit none
         type (c_ptr) :: s_f
         integer(c_int), value :: nsize
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine sse_callocate(s_c,nsize,irc)                        &
     &bind(C,name='sse_callocate')
         use iso_c_binding
         implicit none
         type (c_ptr) :: s_c
         integer(c_int), value :: nsize
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine sse_iallocate(s_i,nsize,irc)                        &
     &bind(C,name='sse_iallocate')
         use iso_c_binding
         implicit none
         type (c_ptr) :: s_i
         integer(c_int), value :: nsize
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine sse_deallocate(s_d,irc) bind(C,name='sse_deallocate'&
     &)
         use iso_c_binding
         implicit none
         type (c_ptr), value :: s_d
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine csse2iscan2(isdata,nths) bind(C,name='csse2iscan2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nths
         type (c_ptr), value :: isdata
         end subroutine
      end interface
!
      interface
         integer (c_int) function check_sse2() bind(C,name='check_sse2')
         use iso_c_binding
         implicit none
         end function
      end interface
!
      interface
         integer (c_int) function check_avx() bind(C,name='check_avx')
         use iso_c_binding
         implicit none
         end function
      end interface
!
      contains
!
!-----------------------------------------------------------------------
      subroutine sse_f2allocate(s_f,nx,ny,irc)
      implicit none 
      integer, intent(in) :: nx, ny
      integer :: irc
      real, dimension(:,:), pointer :: s_f
! local data
      type (c_ptr) :: fptr
      if (digits(s_f(1,1)) > 24) then
         call sse_fallocate(fptr,2*nx*ny,irc)
      else
         call sse_fallocate(fptr,nx*ny,irc)
      endif
      call c_f_pointer(fptr,s_f,(/nx,ny/))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine sse_f3allocate(s_f,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: nx, ny, nz
      integer :: irc
      real, dimension(:,:,:), pointer :: s_f
! local data
      type (c_ptr) :: fptr
      if (digits(s_f(1,1,1)) > 24) then
         call sse_fallocate(fptr,2*nx*ny*nz,irc)
      else
         call sse_fallocate(fptr,nx*ny*nz,irc)
      endif
      call c_f_pointer(fptr,s_f,(/nx,ny,nz/))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine sse_c2allocate(s_f,nx,ny,irc)
      implicit none 
      integer, intent(in) :: nx, ny
      integer :: irc
      complex, dimension(:,:), pointer :: s_f
! local data
      type (c_ptr) :: cptr
      if (digits(real(s_f(1,1))) > 24) then
         call sse_callocate(cptr,2*nx*ny,irc)
      else
         call sse_callocate(cptr,nx*ny,irc)
      endif
      call c_f_pointer(cptr,s_f,(/nx,ny/))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine sse_c3allocate(s_f,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: nx, ny, nz
      integer :: irc
      complex, dimension(:,:,:), pointer :: s_f
! local data
      type (c_ptr) :: cptr
      if (digits(real(s_f(1,1,1))) > 24) then
         call sse_callocate(cptr,2*nx*ny*nz,irc)
      else
         call sse_callocate(cptr,nx*ny*nz,irc)
      endif
      call c_f_pointer(cptr,s_f,(/nx,ny,nz/))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine sse_i1allocate(s_f,nx,irc)
      implicit none 
      integer, intent(in) :: nx
      integer :: irc
      integer, dimension(:), pointer :: s_f
! local data
      type (c_ptr) :: iptr
      if (digits(s_f(1)) > 31) then
         call sse_iallocate(iptr,2*nx,irc)
      else
         call sse_iallocate(iptr,nx,irc)
      endif
      call c_f_pointer(iptr,s_f,(/nx/))
      end subroutine
!
      end module
