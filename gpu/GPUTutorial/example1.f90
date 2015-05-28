!-----------------------------------------------------------------------
! CUDA Fortran GPU Tutorial: Copy
! written by Viktor K. Decyk, UCLA
      program example1
      use copy
      use gpuflib2
      implicit none
! nx, ny = size of array
      integer, parameter :: nx = 3000, ny = 600
! nblock = block size on GPU
      integer :: nblock = 64
! mx, my = data block size
      integer :: mx, my
      real :: s = 0.5
      integer :: j, k, irc
      real :: eps, epsmax
! timing data
      double precision :: dtime
      integer, dimension(4) :: itime
! data for Fortran Host
      real, dimension(nx) :: a1, b1
      real, dimension(nx,ny) :: a2, b2
! data for GPU
      real, device, dimension(:), allocatable :: g_a1, g_b1
      real, device, dimension(:,:), allocatable :: g_a2, g_b2
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
      allocate(g_a1(nx),g_b1(nx))
!
! allocate 2d data on GPU
      allocate(g_a2(nx,ny),g_b2(nx,ny))
!
! initialize 1d data on Host
      do j = 1, nx
         b1(j) = real(j)
      enddo
      a1 = 0.0
! initialize 2d data on Host
      do k = 1, ny
      do j = 1, nx
         b2(j,k) = real(j + nx*(k-1))
      enddo
      enddo
      a2 = 0.0
! copy data to GPU
      g_a1 = a1
      g_b1 = b1
      g_a2 = a2
      g_b2 = b2
!
! measure overhead time by running empty kernel
      call dtimer(dtime,itime,-1)
      call empty_kernel()
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran empty kernel time=', real(dtime)
!
! 1D copy
      mx = 128
!
! segmented 1d copy on Host with block size mx
      call dtimer(dtime,itime,-1)
!     call copy0(a1,b1)
      call copy1(a1,b1,mx)
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran 1d copy time=', real(dtime)
!
! 1d copy on GPU with block size mx
      call dtimer(dtime,itime,-1)
      call gpu_copy1(g_a1,g_b1,mx)
      call dtimer(dtime,itime,1)
      write (*,*) 'GPU 1d copy time=', real(dtime)
!
! copy data from GPU
      b1 = g_a1
!
! Check for correctness: compare a1 and g_a1
      epsmax = 0.0
      do j = 1, nx
         eps = abs(a1(j)-b1(j))
         if (eps > epsmax) epsmax = eps
      enddo
      write (*,*) '1d copy maximum difference = ', epsmax
!
! 2D copy
      mx = 16; my = 16
!
! segmented 2d copy on host with block size mx, my
      call dtimer(dtime,itime,-1)
!     call copy2(a2,b2,mx)
!     call saxpy2(a2,b2,s,mx)
      call copy3(a2,b2,mx,my)
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran 2d copy time=', real(dtime)
!
! 2d copy on GPU with block size mx, my
      call dtimer(dtime,itime,-1)
!     call gpu_copy2a(g_a2,g_b2,mx)
!     call gpu_copy2b(g_a2,g_b2,mx)
!     call gpu_saxpy2(g_a2,g_b2,s,mx)
      call gpu_copy3(g_a2,g_b2,mx,my)
      call dtimer(dtime,itime,1)
      write (*,*) 'GPU 2d copy time=', real(dtime)
!
! copy data from GPU
      b2 = g_a2
!
! Check for correctness: compare a2 and g_a2
      epsmax = 0.0
      do k = 1, ny
         do j = 1, nx
            eps = abs(a2(j,k)-b2(j,k))
            if (eps > epsmax) epsmax = eps
         enddo
      enddo
      write (*,*) '2d copy maximum difference = ', epsmax
!
! deallocate memory on GPU
      deallocate(g_a1,g_b1)
      deallocate(g_a2,g_b2)
! close down GPU
      call end_cu()
!
      end program