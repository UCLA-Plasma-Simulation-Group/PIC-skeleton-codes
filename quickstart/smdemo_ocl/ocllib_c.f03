!-----------------------------------------------------------------------
! Fortran2003 interface file for gpulib_cl.c
      module ocllib_c
      use iso_c_binding
      implicit none
!
      interface
         subroutine gpadd(g_a,g_b,g_c,nx) bind(C,name='gpadd')
         use iso_c_binding
         implicit none
         type (c_ptr), value :: g_a, g_b, g_c
         integer(c_int), value :: nx
         end subroutine
      end interface
!
      interface
         subroutine init_cl(platf,dev,irc) bind(C,name='init_cl')
         use iso_c_binding
         implicit none
         integer(c_int), value :: platf, dev
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine end_cl(irc) bind(C,name='end_cl')
         use iso_c_binding
         implicit none
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine gpu_setgbsize(nblock) bind(C,name='gpu_setgbsize')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nblock
         end subroutine
      end interface
!
      interface
         subroutine gpu_fallocate(g_f,nsize,irc)                        &
     &bind(C,name='gpu_fallocate')
         use iso_c_binding
         implicit none
         type (c_ptr) :: g_f
         integer(c_int), value :: nsize
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine gpu_deallocate(g_d,irc) bind(C,name='gpu_deallocate'&
     &)
         use iso_c_binding
         implicit none
         type (c_ptr) :: g_d
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine gpu_fcopyin(f,g_f,nsize) bind(C,name='gpu_fcopyin')
         use iso_c_binding
         implicit none
         type (c_ptr), value :: f, g_f
         integer(c_int), value :: nsize
         end subroutine
      end interface
!
      interface
         subroutine gpu_fcopyout(f,g_f,nsize) bind(C,name='gpu_fcopyout'&
     &)
         use iso_c_binding
         implicit none
         type (c_ptr), value :: f, g_f
         integer(c_int), value :: nsize
         end subroutine
      end interface
!
      end module
