!-----------------------------------------------------------------------
! Fortran2003 interface to CUDA utility Library
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
         subroutine gpu_callocate(g_c,nsize,irc)                        &
     &bind(C,name='gpu_callocate')
         use iso_c_binding
         implicit none
         type (c_ptr) :: g_c
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
         subroutine gpu_icopyin(f,g_f,nsize) bind(C,name='gpu_icopyin')
         use iso_c_binding
         implicit none
         type (c_ptr), value :: f, g_f
         integer(c_int), value :: nsize
         end subroutine
      end interface
!
      interface
         subroutine gpu_icopyout(f,g_f,nsize) bind(C,name='gpu_icopyout'&
     &)
         use iso_c_binding
         implicit none
         type (c_ptr), value :: f, g_f
         integer(c_int), value :: nsize
         end subroutine
      end interface
!
      interface
         subroutine gpu_ccopyin(f,g_f,nsize) bind(C,name='gpu_ccopyin')
         use iso_c_binding
         implicit none
         type (c_ptr), value :: f, g_f
         integer(c_int), value :: nsize
         end subroutine
      end interface
!
      interface
         subroutine gpu_ccopyout(f,g_f,nsize) bind(C,name='gpu_ccopyout'&
     &)
         use iso_c_binding
         implicit none
         type (c_ptr), value :: f, g_f
         integer(c_int), value :: nsize
         end subroutine
      end interface
!
      interface
         subroutine gpu_zfmem(g_f,nsize) bind(C,name='gpu_zfmem')
         use iso_c_binding
         implicit none
         type (c_ptr), value :: g_f
         integer(c_int), value :: nsize
         end subroutine
      end interface
!
      interface
         subroutine gpu_zcmem(g_f,nsize) bind(C,name='gpu_zcmem')
         use iso_c_binding
         implicit none
         type (c_ptr), value :: g_f
         integer(c_int), value :: nsize
         end subroutine
      end interface
!
      interface
         subroutine gpu_set_cache_size(nscache)                         &
     &bind(C,name='gpu_set_cache_size')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nscache
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
      end module


