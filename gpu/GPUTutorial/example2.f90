!-----------------------------------------------------------------------
! Fortran GPU Tutorial: Transpose
! written by Viktor K. Decyk, UCLA
      program example2
      use  transpose
      use gpuflib2
      implicit none
      integer, parameter :: nx = 512, ny = 512, mx = 16, my = 16
      integer :: nblock = 64
      integer :: j, k, irc
      real :: eps, epsmax
      double precision :: dtime
      integer, dimension(4) :: itime
      real, dimension(nx,ny) :: b2
      real, dimension(ny,nx) :: a2, c2
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
! allocate 2d data on GPU
      allocate(g_a2(ny,nx),g_b2(nx,ny))
!
! initialize 2d data
      do k = 1, ny
      do j = 1, nx
         b2(j,k) = real(j + nx*(k-1))
      enddo
      enddo
      a2 = 0.0
      g_a2 = a2
!
! measure overhead time by running empty kernel
      call dtimer(dtime,itime,-1)
      call empty_kernel()
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran empty kernel time=', real(dtime)
!
! segmented 2d transpose on host with block size mx, my
      call dtimer(dtime,itime,-1)
!     call transpose0(a2,b2)
      call transpose2(a2,b2,mx,my)
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran 2d transpose time=', real(dtime)
!
! 2d transpose on GPU with block size mx, mx
      g_b2 = b2
      call dtimer(dtime,itime,-1)
      call gpu_transpose2(g_a2,g_b2,mx)
      call dtimer(dtime,itime,1)
      write (*,*) 'GPU 2d transpose time=', real(dtime)
      c2 = g_a2
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
!
! deallocate memory on GPU
      deallocate(g_a2,g_b2)
! close down GPU
      call end_cu()
!
      end program