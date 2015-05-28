!-----------------------------------------------------------------------
! Fortran2003 interface to CUDA Library for GPU Tutorial
! written by Viktor K. Decyk, UCLA
      module gpulib2_c
      use iso_c_binding
      implicit none
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
         function getmmcc() bind(C,name='getmmcc')
         use iso_c_binding
         implicit none
         integer(c_int) :: getmmcc
         end function
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
         subroutine gpu_iallocate(g_i,nsize,irc)                        &
     &bind(C,name='gpu_iallocate')
         use iso_c_binding
         implicit none
         type (c_ptr) :: g_i
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
         type (c_ptr), value :: g_d
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
      interface
         subroutine emptykernel() bind(C,name='emptykernel')
         implicit none
         end subroutine
      end interface
!
      interface
         subroutine init_cu(dev,irc) bind(C,name='init_cu')
         use iso_c_binding
         implicit none
         integer(c_int), value :: dev
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine end_cu() bind(C,name='end_cu')
         implicit none
         end subroutine
      end interface
!
      interface
         subroutine gpu_copy1(gp_a,gp_b,mx,nx) bind(C,name='gpu_copy1')
         use iso_c_binding
         implicit none
         integer(c_int), value :: mx, nx
         type (c_ptr), value :: gp_a, gp_b
         end subroutine
      end interface
!
      interface
         subroutine gpu_copy2a(gp_a,gp_b,mx,nx,ny)                      &
     &bind(C,name='gpu_copy2a')
         use iso_c_binding
         implicit none
         integer(c_int), value :: mx, nx, ny
         type (c_ptr), value :: gp_a, gp_b
         end subroutine
      end interface
!
      interface
         subroutine gpu_copy2b(gp_a,gp_b,mx,nx,ny)                      &
     &bind(C,name='gpu_copy2b')
         use iso_c_binding
         implicit none
         integer(c_int), value :: mx, nx, ny
         type (c_ptr), value :: gp_a, gp_b
         end subroutine
      end interface
!
      interface
         subroutine gpu_copy3(gp_a,gp_b,mx,my,nx,ny)                    &
     &bind(C,name='gpu_copy3')
         use iso_c_binding
         implicit none
         integer(c_int), value :: mx, my, nx, ny
         type (c_ptr), value :: gp_a, gp_b
         end subroutine
      end interface
!
      interface
         subroutine gpu_transpose2(gp_a,gp_b,mx,nx,ny)                  &
     &bind(C,name='gpu_transpose2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: mx, nx, ny
         type (c_ptr), value :: gp_a, gp_b
         end subroutine
      end interface
!
      interface
         subroutine gpu_sum1(gp_a,gp_sa,mx,nx) bind(C,name='gpu_sum1')
         use iso_c_binding
         implicit none
         integer(c_int), value :: mx, nx
         type (c_ptr), value :: gp_a, gp_sa
         end subroutine
      end interface
!
      interface
         subroutine gpu_sum2(gp_a,gp_d,mx,nx) bind(C,name='gpu_sum2')
         use iso_c_binding
         implicit none
         integer(c_int), value :: mx, nx
         type (c_ptr), value :: gp_a, gp_d
         end subroutine
      end interface
!
      interface
         subroutine gpu_sum3(gp_a,gp_d,gp_sa,mx,nx)                     &
     &bind(C,name='gpu_sum3')
         use iso_c_binding
         implicit none
         integer(c_int), value :: mx, nx
         type (c_ptr), value :: gp_a, gp_d, gp_sa
         end subroutine
      end interface
!
      end module


