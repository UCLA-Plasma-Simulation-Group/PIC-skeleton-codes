!-----------------------------------------------------------------------
! Fortran interface to CUDA Library for GPU Tutorial
! written by Viktor K. Decyk, UCLA
      module gpulib2
      implicit none

      interface
         subroutine setgbsize(nblock)
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
      interface
         subroutine gpu_copy1(gp_a,gp_b,mx,nx)
         implicit none
         integer :: mx, nx
         integer, dimension(2) :: gp_a, gp_b
         end subroutine
      end interface
!
      interface
         subroutine gpu_copy2a(gp_a,gp_b,mx,nx,ny)
         implicit none
         integer :: mx, nx, ny
         integer, dimension(2) :: gp_a, gp_b
         end subroutine
      end interface
!
      interface
         subroutine gpu_copy2b(gp_a,gp_b,mx,nx,ny)
         implicit none
         integer :: mx, nx, ny
         integer, dimension(2) :: gp_a, gp_b
         end subroutine
      end interface
!
      interface
         subroutine gpu_copy3(gp_a,gp_b,mx,my,nx,ny)
         implicit none
         integer :: mx, my, nx, ny
         integer, dimension(2) :: gp_a, gp_b
         end subroutine
      end interface
!
      interface
         subroutine gpu_transpose2(gp_a,gp_b,mx,nx,ny)
         implicit none
         integer :: mx, my, nx, ny
         integer, dimension(2) :: gp_a, gp_b
         end subroutine
      end interface
!
      interface
         subroutine gpu_sum1(gp_a,gp_sa,mx,nx)
         implicit none
         integer :: mx, nx
         integer, dimension(2) :: gp_a, gp_sa
         end subroutine
      end interface
!
      interface
         subroutine gpu_sum2(gp_a,gp_d,mx,nx)
         implicit none
         integer :: mx, nx
         integer, dimension(2) :: gp_a, gp_d
         end subroutine
      end interface
!
      interface
         subroutine gpu_sum3(gp_a,gp_d,gp_sa,mx,nx)
         implicit none
         integer :: mx, nx
         integer, dimension(2) :: gp_a, gp_d, gp_sa
         end subroutine
      end interface
!
      end module


