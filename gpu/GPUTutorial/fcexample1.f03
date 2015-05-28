!-----------------------------------------------------------------------
! Fortran2003 with CUDA C GPU Tutorial: Copy
! written by Viktor K. Decyk, UCLA
      program fcexample1
      use iso_c_binding
      use copy
      use gpulib2_c
      implicit none
! nx, ny = size of array
      integer, parameter :: nx = 3000, ny = 600
! nblock = block size on GPU
      integer :: nblock = 64
! mx, my = data block size
      integer :: mx, my
      real :: s = 0.5
      integer :: j, k, irc
! timing data
      real :: eps, epsmax
      double precision :: dtime
      integer, dimension(4) :: itime
! data for Fortran Host
      real, dimension(nx), target :: a1, b1
      real, dimension(nx,ny), target :: a2, b2
! data for GPU
      type (c_ptr) :: g_a1 = c_null_ptr, g_b1 = c_null_ptr
      type (c_ptr) :: g_a2 = c_null_ptr, g_b2 = c_null_ptr
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
      call gpu_fallocate(g_a1,nx,irc)
      call gpu_fallocate(g_b1,nx,irc)
!
! allocate 2d data on GPU
      call gpu_fallocate(g_a2,nx*ny,irc)
      call gpu_fallocate(g_b2,nx*ny,irc)
!
      if (irc /= 0) then
         write (*,*) 'GPU allocate error!'
         stop
      endif
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
      call gpu_fcopyin(c_loc(a1),g_a1,nx)
      call gpu_fcopyin(c_loc(b1),g_b1,nx)
      call gpu_fcopyin(c_loc(a2),g_a2,nx*ny)
      call gpu_fcopyin(c_loc(b2),g_b2,nx*ny)
!
! measure overhead time by running empty kernel
      call dtimer(dtime,itime,-1)
      call emptykernel()
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
      call gpu_copy1(g_a1,g_b1,mx,nx)
      call dtimer(dtime,itime,1)
      write (*,*) 'GPU 1d copy time=', real(dtime)
!
! copy data from GPU
      call gpu_fcopyout(c_loc(b1),g_a1,nx);
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
      mx = 32; my = 16
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
!     call gpu_copy2a(g_a2,g_b2,mx,nx,ny)
!     call gpu_copy2b(g_a2,g_b2,mx,nx,ny)
!     call gpu_saxpy2(g_a2,g_b2,s,mx,nx,ny)
      call gpu_copy3(g_a2,g_b2,mx,my,nx,ny)
      call dtimer(dtime,itime,1)
      write (*,*) 'GPU 2d copy time=', real(dtime)
!
! copy data from GPU
      call gpu_fcopyout(c_loc(b2),g_a2,nx*ny)
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
      call gpu_deallocate(g_a1,irc)
      call gpu_deallocate(g_b1,irc)
      call gpu_deallocate(g_a2,irc)
      call gpu_deallocate(g_b2,irc)
      if (irc /= 0) then
         write (*,*) 'GPU deallocate error!'
         stop
      endif
! close down GPU
      call end_cu()
!
      end program