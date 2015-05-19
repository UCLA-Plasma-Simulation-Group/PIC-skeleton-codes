!-----------------------------------------------------------------------
      program omptest
! OpenMP vector add test program
! written by Viktor K. Decyk, UCLA
! nx = size of array, nthreads = number of threads
      implicit none
      integer :: nx = 1048576, nthreads = 1
      integer :: j, irc
      real :: time = 0.0
      real :: eps, epsmax
! timing data
      double precision :: dtime
      integer, dimension(4) :: itime
! data for Fortran
      real, dimension(:), pointer :: a, b, c
! data for OpenMP
      real, dimension(:), pointer :: p_a, p_b, p_c
!
! initialize Host data
      allocate(a(nx),b(nx),c(nx))
! initialize vectors
      do j = 1, nx
         b(j) = j
         c(j) = 2*j
      enddo
      a = 0.0
! set up OpenMP
      call init_omp(0)
!     call setnthsize(nthreads)
! allocate data for OpenMP
      allocate(p_a(nx),p_b(nx),p_c(nx))
! Copy initial data for OpenMP
      p_b = b; p_c = c
!
! First execute on Host in Fortran: a = b + c
      call dtimer(dtime,itime,-1)
      call fadd(a,b,c,nx)
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran add time=', real(dtime)
!
! Execute with OpenMP: p_a = p_b + p_c
      call dtimer(dtime,itime,-1)
      call mpadd(p_a,p_b,p_c,nx)
      call dtimer(dtime,itime,1)
      time = time + real(dtime)
      write (*,*) 'OpenMP add time=', real(dtime)
!
! Check for correctness: compare a and p_a
      epsmax = 0.0
      do j = 1, nx
         eps = abs(a(j)-p_a(j))
         if (eps > epsmax) epsmax = eps
      enddo
      write (*,*) 'maximum difference = ', epsmax
!
! deallocate memory for OpenMP
      deallocate(p_a,p_b,p_c)
! deallocate Host memory
      deallocate(a,b,c)
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
      
