!-----------------------------------------------------------------------
      program gputest
! GPU vector add test program
! written by Viktor K. Decyk, UCLA
! nx = size of array, nblock = block size on GPU
      use gpulib_cuf
      implicit none
! nx = size of array, nblock = block size on GPU
      integer :: nx = 1048576, nblock = 64
      integer :: j, irc
      real :: eps, epsmax
! timing data
      double precision :: dtime
      integer, dimension(4) :: itime
! data for Fortran Host
      real, dimension(:), pointer :: a => null(), b => null()
      real, dimension(:), pointer :: c => null()
! data for GPU
      real, device, dimension(:), allocatable :: g_a, g_b, g_c
!
! initialize Fortran data on Host
      allocate(a(nx),b(nx),c(nx))
! initialize vectors
      do j = 1, nx
         b(j) = j
         c(j) = 2*j
      enddo
      a = 0.0
! set up GPU
      irc = 0
      call setgbsize(nblock)
      call init_cuf(0,irc)
      if (irc /= 0) then
         write (*,*) 'CUDA initialization error!'
         stop
      endif
! allocate data on GPU
      allocate(g_a(nx),g_b(nx),g_c(nx))
! Copy initial data to GPU
      g_b = b; g_c = c
!
! First execute on Host in Fortran: a = b + c
      call dtimer(dtime,itime,-1)
      call fadd(a,b,c,nx)
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran add time=', real(dtime)
!
! Execute on GPU: g_a = g_b + g_c
      call dtimer(dtime,itime,-1)
      call gpadd(g_a,g_b,g_c,nx)
      call dtimer(dtime,itime,1)
      write (*,*) 'GPU add time=', real(dtime)
!
! Check for correctness: compare a and g_a
! Copy g_a from GPU to c on Host, then compare a with c
      c = g_a
      epsmax = 0.0
      do j = 1, nx
         eps = abs(a(j)-c(j))
         if (eps > epsmax) epsmax = eps
      enddo
      write (*,*) 'maximum difference = ', epsmax
!
! deallocate memory on GPU
      deallocate(g_a,g_b,g_c)
! close down GPU
      call end_cuf()
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
      
