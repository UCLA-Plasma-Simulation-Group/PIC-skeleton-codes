!-----------------------------------------------------------------------
! Interface file for sseflib2.c
      module sseflib2_h
      use sselib2_h
      implicit none
!
      interface
         subroutine sse_f2allocatex(sp_f,nx,ny,nd,irc)
         implicit none
         integer :: nx, ny, nd, irc
         real, dimension(:,:), pointer :: sp_f
         end subroutine
      end interface
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
         subroutine sse_c2allocatex(sp_c,nx,ny,nd,irc)
         implicit none
         integer :: nx, ny, nd, irc
         complex, dimension(:,:), pointer :: sp_c
         end subroutine
      end interface
!
      interface
         subroutine sse_i1allocatex(sp_i,nx,nd,irc)
         implicit none
         integer :: nx, nd, irc
         integer, dimension(:), pointer :: sp_i
         end subroutine
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
      if (digits(s_f(1,1)) > 24) then
         call sse_f2allocatex(s_f,nx,ny,2,irc)
      else
         call sse_f2allocatex(s_f,nx,ny,1,irc)
      endif
      end subroutine
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
      subroutine sse_c2allocate(s_f,nx,ny,irc)
      implicit none 
      integer, intent(in) :: nx, ny
      integer :: irc
      complex, dimension(:,:), pointer :: s_f
      if (digits(real(s_f(1,1))) > 24) then
         call sse_c2allocatex(s_f,nx,ny,2,irc)
      else
         call sse_c2allocatex(s_f,nx,ny,1,irc)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine sse_i1allocate(s_f,nx,irc)
      implicit none 
      integer, intent(in) :: nx
      integer :: irc
      integer, dimension(:), pointer :: s_f
      if (digits(s_f(1)) > 31) then
         call sse_i1allocatex(s_f,nx,2,irc)
      else
         call sse_i1allocatex(s_f,nx,1,irc)
      endif
      end subroutine
! 
   end module
