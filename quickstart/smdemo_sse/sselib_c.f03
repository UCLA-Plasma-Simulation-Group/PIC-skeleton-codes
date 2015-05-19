!-----------------------------------------------------------------------
! Fortran2003 interface file for sselib.c
      module sselib_c
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
         subroutine ssadd(s_a,s_b,s_c,nx) bind(C,name='ssadd')
         use iso_c_binding
         implicit none
         type (c_ptr), value :: s_a, s_b, s_c
         integer(c_int), value :: nx
         end subroutine
      end interface
!
      contains
!
      subroutine sse_f1allocate(s_f,nx,irc)
      implicit none 
      integer, intent(in) :: nx
      integer, intent(inout) :: irc
      real, dimension(:), pointer :: s_f
! local data
      type (c_ptr) :: fptr
      call sse_fallocate(fptr,nx,irc)
      call c_f_pointer(fptr,s_f,(/nx/))
      end subroutine
!
      subroutine fssadd(s_a,s_b,s_c,nx)
! Fortran interface to SSE function ssadd
      implicit none 
      integer, intent(in) :: nx
      real, dimension(:), pointer :: s_a, s_b, s_c
      call ssadd(c_loc(s_a(1)),c_loc(s_b(1)),c_loc(s_c(1)),nx)
      end subroutine
!
      end module
