!-----------------------------------------------------------------------
! Fortran with CUDA C GPU Tutorial: Transpose
! written by Viktor K. Decyk, UCLA
      program example2
      use  transpose
      use gpulib2
      implicit none
      integer, parameter :: nx = 512, ny = 512, mx = 16, my = 16
      integer :: nblock = 64
      integer :: j, k, irc
      real :: eps, epsmax
      double precision :: dtime
      integer, dimension(4) :: itime
      real, dimension(nx,ny) :: b2
      real, dimension(ny,nx) :: a2, c2
      integer, dimension(2) :: g_a2 = 0.0, g_b2 = 0.0
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
! initialize 2d data
      do k = 1, ny
      do j = 1, nx
         b2(j,k) = real(j + nx*(k-1))
      enddo
      enddo
      a2 = 0.0
      call gpu_fcopyin(a2,g_a2,ny*nx)
!
! measure overhead time by running empty kernel
      call dtimer(dtime,itime,-1)
      call emptykernel()
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
      call gpu_fcopyin(b2,g_b2,nx*ny)
      call dtimer(dtime,itime,-1)
      call gpu_transpose2(g_a2,g_b2,mx,nx,ny)
      call dtimer(dtime,itime,1)
      write (*,*) 'GPU 2d transpose time=', real(dtime)
      call gpu_fcopyout(c2,g_a2,ny*nx)
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