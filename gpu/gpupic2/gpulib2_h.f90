!-----------------------------------------------------------------------
! Fortran interface to CUDA utility Library
! written by Viktor K. Decyk, UCLA
      module gpulib2_h
      implicit none
!
      interface
         subroutine gpu_setgbsize(nblock)
         implicit none
         integer :: nblock
         end subroutine
      end interface
!
      interface
         integer function getmmcc()
         implicit none
         end function
      end interface
!
      interface
         subroutine gpu_fallocate(gp_f,nsize,irc)
         implicit none
         integer :: nsize, irc
         integer, dimension(2) :: gp_f
         end subroutine
      end interface
!
      interface
         subroutine gpu_iallocate(gp_i,nsize,irc)
         implicit none
         integer :: nsize, irc
         integer, dimension(2) :: gp_i
         end subroutine
      end interface
!
      interface
         subroutine gpu_callocate(gp_c,nsize,irc)
         implicit none
         integer :: nsize, irc
         integer, dimension(2) :: gp_c
         end subroutine
      end interface
!
      interface
         subroutine gpu_deallocate(gp_d,irc)
         implicit none
         integer :: irc
         integer, dimension(2) :: gp_d
         end subroutine
      end interface
!
      interface
         subroutine gpu_fcopyin(f,gp_f,nsize)
         implicit none
         integer :: nsize
         real, dimension(*) :: f
         integer, dimension(2) :: gp_f
         end subroutine
      end interface
!
      interface
         subroutine gpu_fcopyout(f,gp_f,nsize)
         implicit none
         integer :: nsize
         real, dimension(*) :: f
         integer, dimension(2) :: gp_f
         end subroutine
      end interface
!
      interface
         subroutine gpu_icopyin(f,gp_f,nsize)
         implicit none
         integer :: nsize
         integer, dimension(*) :: f
         integer, dimension(2) :: gp_f
         end subroutine
      end interface
!
      interface
         subroutine gpu_icopyout(f,gp_f,nsize)
         implicit none
         integer :: nsize
         integer, dimension(*) :: f
         integer, dimension(2) :: gp_f
         end subroutine
      end interface
!
      interface
         subroutine gpu_ccopyin(f,gp_f,nsize)
         implicit none
         integer :: nsize
         complex, dimension(*) :: f
         integer, dimension(2) :: gp_f
         end subroutine
      end interface
!
      interface
         subroutine gpu_ccopyout(f,gp_f,nsize)
         implicit none
         integer :: nsize
         complex, dimension(*) :: f
         integer, dimension(2) :: gp_f
         end subroutine
      end interface
!
      interface
         subroutine gpu_zfmem(gp_f,nsize)
         implicit none
         integer :: nsize
         integer, dimension(2) :: gp_f
         end subroutine
      end interface
!
      interface
         subroutine gpu_set_cache_size(nscache)
         implicit none
         integer :: nscache
         end subroutine
      end interface
!
      interface
         subroutine emptykernel()
         implicit none
         end subroutine
      end interface
!
      interface
         subroutine init_cu(dev,irc)
         implicit none
         integer :: dev, irc
         end subroutine
      end interface
!
      interface
         subroutine end_cu()
         implicit none
         end subroutine
      end interface
!
      end module


