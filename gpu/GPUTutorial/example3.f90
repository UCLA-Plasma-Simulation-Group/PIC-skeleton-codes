!-----------------------------------------------------------------------
! CUDA Fortran GPU Tutorial: Reduction
! written by Viktor K. Decyk, UCLA
      program example3
      use redux
      use gpuflib2
      implicit none
      integer, parameter :: nx = 3000, mx = 128, nbx = (nx - 1)/mx + 1
      integer, parameter :: nbxs = (nbx - 1)/mx + 1
      integer :: nblock = 64
      integer :: j, k, irc
      real :: s, t, eps, epsmax
      double precision :: dtime
      integer, dimension(4) :: itime
      real, dimension(nx) :: a1
      real, dimension(nbx) :: d1
      real, device, dimension(:), allocatable :: g_a1, g_d1, g_s
!
! set up GPU
      irc = 0
      call setgbsize(nblock)
      call init_cu(0,irc)
      if (irc /= 0) then
         write (*,*) 'CUDA initialization error!'
         stop
      endif
!
! allocate 1d data on GPU
      allocate(g_a1(nx),g_d1(nbx),g_s(2*nbxs))
!
! initialize 1d data on host
      do j = 1, nx
         a1(j) = real(j)
      enddo
      d1 = 0.0
      g_d1 = d1
      s = 0.0; t = 0.0
!
! measure overhead time by running empty kernel
      call dtimer(dtime,itime,-1)
      call empty_kernel()
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran empty kernel time=', real(dtime)
!
! segmented 1d sum on host with block size mx
      call dtimer(dtime,itime,-1)
!     call sum0(a1,s)
!     call sum1(a1,s,mx)
      call sum2(a1,d1,mx)
      call sum1(d1,s,mx)
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran 1d sum time=', real(dtime)
!
! 1d sum on GPU with block size mx
      g_a1 = a1
      call dtimer(dtime,itime,-1)
!     call gpu_sum1(g_a1,g_s,mx)
!     call gpu_sum2(g_a1,g_d1,mx)
!     call gpu_sum1(g_d1,g_s,mx)
      call gpu_sum3(g_a1,g_d1,g_s,mx)
      call dtimer(dtime,itime,1)
      write (*,*) 'GPU 1d sum time=', real(dtime)
      a1(1:nbx) = g_d1; t = g_s(1)
!
! Check for correctness: compare d1 and g_d1
      epsmax = 0.0
      do j = 1, nbx
         eps = abs(d1(j)-a1(j))
         if (eps > epsmax) epsmax = eps
      enddo
      write (*,*) '1d sum maximum difference = ', epsmax
      write (*,*) 's,t = ', s, t
!
! deallocate memory on GPU
      deallocate(g_a1,g_d1,g_s)
! close down GPU
      call end_cu()
!
      end program