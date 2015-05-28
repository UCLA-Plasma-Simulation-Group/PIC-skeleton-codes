!-----------------------------------------------------------------------
! Fortran2003 with CUDA C GPU Tutorial: Transpose
! written by Viktor K. Decyk, UCLA
      program fcexample2
      use iso_c_binding
      use transpose
      use gpulib2_c
      implicit none
! nx, ny = size of array
! mx, my = data block size
      integer, parameter :: nx = 512, ny = 512, mx = 16, my = 16
! nblock = block size on GPU
      integer :: nblock = 64
      integer :: j, k, irc
      real :: eps, epsmax
! timing data
      double precision :: dtime
      integer, dimension(4) :: itime
! data for Fortran Host
      real, dimension(nx,ny), target :: b2
      real, dimension(ny,nx), target :: a2, c2
! data for GPU
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
! allocate 2d data on GPU
      call gpu_fallocate(g_a2,ny*nx,irc)
      call gpu_fallocate(g_b2,nx*ny,irc)

      if (irc /= 0) then
         write (*,*) 'GPU allocate error!'
         stop
      endif
!
! initialize 2d data on Host
      do k = 1, ny
      do j = 1, nx
         b2(j,k) = real(j + nx*(k-1))
      enddo
      enddo
      a2 = 0.0
! copy data to GPU
      call gpu_fcopyin(c_loc(a2),g_a2,ny*nx)
      call gpu_fcopyin(c_loc(b2),g_b2,nx*ny)
!
! measure overhead time by running empty kernel
      call dtimer(dtime,itime,-1)
      call emptykernel()
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran empty kernel time=', real(dtime)
!
! segmented 2d transpose on Host with block size mx, my
      call dtimer(dtime,itime,-1)
!     call transpose0(a2,b2)
      call transpose2(a2,b2,mx,my)
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran 2d transpose time=', real(dtime)
!
! 2d transpose on GPU with block size mx, mx
      call dtimer(dtime,itime,-1)
      call gpu_transpose2(g_a2,g_b2,mx,nx,ny)
      call dtimer(dtime,itime,1)
      write (*,*) 'GPU 2d transpose time=', real(dtime)
!
! copy data from GPU
      call gpu_fcopyout(c_loc(c2),g_a2,ny*nx)
!
! Check for correctness: compare a2 and g_a2
      epsmax = 0.0
      do j = 1, nx
         do k = 1, ny
            eps = abs(a2(k,j)-c2(k,j))
            if (eps > epsmax) epsmax = eps
         enddo
      enddo
      write (*,*) '2d transpose maximum difference = ', epsmax
!
! deallocate memory on GPU
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