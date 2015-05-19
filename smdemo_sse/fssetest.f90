!-----------------------------------------------------------------------
! SSE vector add test program
! written by Viktor K. Decyk, UCLA
! nx = size of array
      program ssetest
      use iso_c_binding
      use sselib_c
      implicit none
      integer :: nx = 1048576
      integer :: j, irc = 0
      real :: eps, epsmax
! timing data
      double precision :: dtime
      integer, dimension(4) :: itime
! data for Fortran
      real, dimension(:), allocatable :: a, b, c
! data for SSE
      real, dimension(:), pointer :: s_a => null()
      real, dimension(:), pointer :: s_b => null()
      real, dimension(:), pointer :: s_c => null()
!
! initialize Fortran data
      allocate(a(nx),b(nx),c(nx))
! initialize vectors
      a = 0.0
      do j = 1, nx
         b(j) = j
         c(j) = 2*j
      enddo
! allocate aligned 1d array for SSE
      call sse_f1allocate(s_a,nx,irc)
      call sse_f1allocate(s_b,nx,irc)
      call sse_f1allocate(s_c,nx,irc)
      if (irc /= 0) then
         write (*,*) 'SSE allocate error!'
         stop
      endif
! Copy initial data for SSE
      s_b = b; s_c = c
!
! First execute in Fortran: a = b + c
      call dtimer(dtime,itime,-1)
      call fadd(a,b,c,nx)
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran add time=', real(dtime)
!
! Execute on SSE: s_a = s_b + s_c
      call dtimer(dtime,itime,-1)
      call fssadd(s_a,s_b,s_c,nx)
      call dtimer(dtime,itime,1)
      write (*,*) 'SSE add time=', real(dtime)
!
! Check for correctness: compare a and s_a
      epsmax = 0.0
      do j = 1, nx
         eps = abs(a(j)-s_a(j))
         if (eps > epsmax) epsmax = eps
      enddo
      write (*,*) 'maximum difference = ', epsmax
!
! deallocate memory for SSE
      call sse_deallocate(c_loc(s_a(1)),irc); nullify(s_a)
      call sse_deallocate(c_loc(s_b(1)),irc); nullify(s_b)
      call sse_deallocate(c_loc(s_c(1)),irc); nullify(s_c)
! deallocate Fortran memory
      deallocate(a,b,c)
!
      end program
!
      subroutine fadd(a,b,c,nx)
! Fortran add procedure
      integer :: nx
      real, dimension(nx) :: a, b, c
      integer :: j
      do j = 1, nx
         a(j) = b(j) + c(j)
      enddo
      end subroutine

