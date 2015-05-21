!-----------------------------------------------------------------------
! Interface file for sseflib3.c
      module sseflib3_h
      use sselib3_h
      implicit none
!
      interface
         subroutine sse_f3allocatex(sp_f,nx,ny,nz,nd,irc)
         implicit none
         integer :: nx, ny, nz, nd, irc
         real, dimension(:,:,:), pointer :: sp_f
         end subroutine
      end interface
!
      interface
         subroutine sse_f4allocatex(sp_f,ndim,nx,ny,nz,nd,irc)
         implicit none
         integer :: ndim, nx, ny, nz, nd, irc
         real, dimension(:,:,:,:), pointer :: sp_f
         end subroutine
      end interface
!
      interface
         subroutine sse_c3allocatex(sp_c,nx,ny,nz,nd,irc)
         implicit none
         integer :: nx, ny, nz, nd, irc
         complex, dimension(:,:,:), pointer :: sp_c
         end subroutine
      end interface
!
      contains
!
!-----------------------------------------------------------------------
      subroutine sse_f3allocate(s_f,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: nx, ny, nz
      integer :: irc
      real, dimension(:,:,:), pointer :: s_f
      if (digits(s_f(1,1,1)) > 24) then
         call sse_f3allocatex(s_f,nx,ny,nz,2,irc)
      else
         call sse_f3allocatex(s_f,nx,ny,nz,1,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine sse_f4allocate(s_f,ndim,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: ndim, nx, ny, nz
      integer :: irc
      real, dimension(:,:,:,:), pointer :: s_f
      if (digits(s_f(1,1,1,1)) > 24) then
         call sse_f4allocatex(s_f,ndim,nx,ny,nz,2,irc)
      else
         call sse_f4allocatex(s_f,ndim,nx,ny,nz,1,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine sse_c3allocate(s_f,nx,ny,nz,irc)
      implicit none 
      integer, intent(in) :: nx, ny, nz
      integer :: irc
      complex, dimension(:,:,:), pointer :: s_f
      if (digits(real(s_f(1,1,1))) > 24) then
         call sse_c3allocatex(s_f,nx,ny,nz,2,irc)
      else
         call sse_c3allocatex(s_f,nx,ny,nz,1,irc)
      endif
      end subroutine
! 
   end module
