!-----------------------------------------------------------------------
! Interface file for avx512flib3.c
      module avx512flib3_h
      implicit none
!
      interface
         subroutine avx512_f3allocatex(sp_f,nx,ny,nz,nd,irc)
         implicit none
         integer :: nx, ny, nz, nd, irc
         real, dimension(:,:,:), pointer :: sp_f
         end subroutine
      end interface
!
      interface
         subroutine avx512_f4allocatex(sp_f,ndim,nx,ny,nz,nd,irc)
         implicit none
         integer :: ndim, nx, ny, nz, nd, irc
         real, dimension(:,:,:,:), pointer :: sp_f
         end subroutine
      end interface
!
      interface
         subroutine avx512_c3allocatex(sp_c,nx,ny,nz,nd,irc)
         implicit none
         integer :: nx, ny, nz, nd, irc
         complex, dimension(:,:,:), pointer :: sp_c
         end subroutine
      end interface
!
      interface
         subroutine avx512_c4allocatex(sp_c,ndim,nx,ny,nz,nd,irc)
         implicit none
         integer :: ndim, nx, ny, nz, nd, irc
         complex, dimension(:,:,:,:), pointer :: sp_c
         end subroutine
      end interface
!
      contains
!
!-----------------------------------------------------------------------
      subroutine avx512_f3allocate(s_f,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: nx, ny, nz
      integer :: irc
      real, dimension(:,:,:), pointer :: s_f
      if (digits(s_f(1,1,1)) > 24) then
         call avx512_f3allocatex(s_f,nx,ny,nz,2,irc)
      else
         call avx512_f3allocatex(s_f,nx,ny,nz,1,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine avx512_f4allocate(s_f,ndim,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: ndim, nx, ny, nz
      integer :: irc
      real, dimension(:,:,:,:), pointer :: s_f
      if (digits(s_f(1,1,1,1)) > 24) then
         call avx512_f4allocatex(s_f,ndim,nx,ny,nz,2,irc)
      else
         call avx512_f4allocatex(s_f,ndim,nx,ny,nz,1,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine avx512_c3allocate(s_f,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: nx, ny, nz
      integer :: irc
      complex, dimension(:,:,:), pointer :: s_f
      if (digits(real(s_f(1,1,1))) > 24) then
         call avx512_c3allocatex(s_f,nx,ny,nz,2,irc)
      else
         call avx512_c3allocatex(s_f,nx,ny,nz,1,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine avx512_c4allocate(s_f,ndim,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: ndim, nx, ny, nz
      integer :: irc
      complex, dimension(:,:,:,:), pointer :: s_f
      if (digits(real(s_f(1,1,1,1))) > 24) then
         call avx512_c4allocatex(s_f,ndim,nx,ny,nz,2,irc)
      else
         call avx512_c4allocatex(s_f,ndim,nx,ny,nz,1,irc)
      endif
      end subroutine
! 
   end module
