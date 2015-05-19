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
         subroutine hpl_fallocate(h_f,nsize,irc)                        &
     &bind(C,name='hpl_fallocate')
         use iso_c_binding
         implicit none 
         type (c_ptr) :: h_f
         integer(c_int), value :: nsize
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine hpl_callocate(h_c,nsize,irc)                        &
     &bind(C,name='hpl_callocate')
         use iso_c_binding
         implicit none 
         type (c_ptr) :: h_c
         integer(c_int), value :: nsize
         integer(c_int) :: irc
         end subroutine
      end interface
!
      interface
         subroutine hpl_deallocate(h_d,irc) bind(C,name='hpl_deallocate'&
     &)
         use iso_c_binding
         implicit none
         type (c_ptr), value :: h_d
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
         subroutine gpu_initstream(nstream)                             &
     &bind(C,name='gpu_initstream')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nstream
         end subroutine
      end interface
!
      interface
         subroutine gpu_delstream(nstream) bind(C,name='gpu_delstream')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nstream
         end subroutine
      end interface
!
      interface
         subroutine waitstream(nstream) bind(C,name='waitstream')
         use iso_c_binding
         implicit none
         integer(c_int), value :: nstream
         end subroutine
      end interface
!
      interface
         subroutine gpu_cascopyin(f,g_f,noff,nsize,nstream)             &
     &bind(C,name='gpu_cascopyin')
         use iso_c_binding
         implicit none
         type (c_ptr), value :: f, g_f
         integer(c_int), value :: noff, nsize, nstream
         end subroutine
      end interface
!
      interface
         subroutine gpu_cascopyout(f,g_f,noff,nsize,nstream)            &
     &bind(C,name='gpu_cascopyout')
         use iso_c_binding
         implicit none
         type (c_ptr), value :: f, g_f
         integer(c_int), value :: noff, nsize, nstream
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
      contains
!
!-----------------------------------------------------------------------
      subroutine hpl_f1allocate(h_f,nx,irc)
      implicit none 
      integer, intent(in) :: nx
      integer :: irc
      real, dimension(:), pointer :: h_f
! local data
      type (c_ptr) :: fptr
      call hpl_fallocate(fptr,nx,irc)
      call c_f_pointer(fptr,h_f,(/nx/))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine hpl_f2allocate(h_f,nx,ny,irc)
      implicit none 
      integer, intent(in) :: nx, ny
      integer :: irc
      real, dimension(:,:), pointer :: h_f
! local data
      type (c_ptr) :: fptr
      call hpl_fallocate(fptr,nx*ny,irc)
      call c_f_pointer(fptr,h_f,(/nx,ny/))
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine hpl_c2allocate(h_f,nx,ny,irc)
      implicit none 
      integer, intent(in) :: nx, ny
      integer :: irc
      complex, dimension(:,:), pointer :: h_f
! local data
      type (c_ptr) :: cptr
      call hpl_callocate(cptr,nx*ny,irc)
      call c_f_pointer(cptr,h_f,(/nx,ny/))
      end subroutine
!
      end module


