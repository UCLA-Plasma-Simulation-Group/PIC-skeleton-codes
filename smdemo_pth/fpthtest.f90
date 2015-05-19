!-----------------------------------------------------------------------
      program pthtest
! Pthreads vector add test program
! written by Viktor K. Decyk, UCLA
! nx = size of array, nthreads = number of threads
      implicit none
      integer :: nx = 1048576, nthreads = 1
      integer :: j, irc
      real :: time = 0.0
      real :: eps, epsmax
      double precision :: dtime
      real, dimension(:), pointer :: a, b, c, d
      integer, dimension(4) :: itime
      real, dimension(:), pointer :: p_a, p_b, p_c
!
! initialize Host data
      allocate(a(nx),b(nx),c(nx),d(nx))
      do j = 1, nx
         b(j) = j
         c(j) = 2*j
      enddo
      a = 0.0; d = -1.0
! set up pthreads
      irc = 0
      call init_pt(0,irc)
      if (irc /= 0) then
         write (*,*) 'Pthreads initialization error!'
         stop
      endif
!     call setnthsize(nthreads)
! allocate data for pthreads
      allocate(p_a(nx),p_b(nx),p_c(nx))
!
! First execute on Host in Fortran: a = b + c
      call dtimer(dtime,itime,-1)
      call fadd(a,b,c,nx)
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran add time=', real(dtime)
! Copy data for pthreads
      call dtimer(dtime,itime,-1)
      p_b = b
      p_c = c
      call dtimer(dtime,itime,1)
      time = time + real(dtime)
      write (*,*) 'Copyin time=', real(dtime)
! Execute with pthreads: g_a = g_b + g_c
      call dtimer(dtime,itime,-1)
      call ptadd(p_a,p_b,p_c,nx,irc)
      call dtimer(dtime,itime,1)
      time = time + real(dtime)
      if (irc /= 0) write (*,*) 'ptadd error: irc=',irc
      write (*,*) 'pthreads add time=', real(dtime)
! Copy data from pthreads: d = g_a
      call dtimer(dtime,itime,-1)
      d = p_a
      call dtimer(dtime,itime,1)
      time = time + real(dtime)
      write (*,*) 'Copyout time=', real(dtime)
      write (*,*) 'Total pthreads time=',time
!
! Check for correctness: compare a and d
      epsmax = 0.0
      do j = 1, nx
         eps = abs(a(j)-d(j))
         if (eps > epsmax) epsmax = eps
      enddo
      write (*,*) 'maximum difference = ', epsmax
!
! deallocate memory for pthreads
      deallocate(p_a,p_b,p_c)
! close down pthreads
      call end_pt()
! deallocate Host memory
      deallocate(a,b,c,d)
!
      end program
!
      subroutine fadd(a,b,c,nx)
      integer :: nx
      real, dimension(nx) :: a, b, c
      integer :: j
      do j = 1, nx
         a(j) = b(j) + c(j)
      enddo
      end subroutine
      
