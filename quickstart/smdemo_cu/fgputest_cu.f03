!-----------------------------------------------------------------------
! GPU vector add test program
! written by Viktor K. Decyk, UCLA
      program gputest
      use iso_c_binding
      use gpulib_c
      implicit none
! nx = size of array, nblock = block size on GPU
      integer :: nx = 1048576, nblock = 64
      integer :: j, irc
      real :: eps, epsmax
! timing data
      double precision :: dtime
      integer, dimension(4) :: itime
! data for Fortran
      real, dimension(:), pointer :: a => null(), b => null()
      real, dimension(:), pointer :: c => null()
! data for GPU
      type (c_ptr) :: g_a = c_null_ptr, g_b = c_null_ptr
      type (c_ptr) :: g_c = c_null_ptr
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
      call init_cu(0,irc)
      if (irc /= 0) then
         write (*,*) 'CUDA initialization error!'
         stop
      endif
! allocate data on GPU, return address to Fortran
      call gpu_fallocate(g_a,nx,irc)
      call gpu_fallocate(g_b,nx,irc)
      call gpu_fallocate(g_c,nx,irc)
      if (irc /= 0) then
         write (*,*) 'GPU allocate error!'
         stop
      endif
! Copy initial data to GPU
      call gpu_fcopyin(c_loc(b(1)),g_b,nx)
      call gpu_fcopyin(c_loc(c(1)),g_c,nx)
!
! First execute on Host in Fortran: a = b + c
      call dtimer(dtime,itime,-1)
      call fadd(a,b,c,nx)
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran add time=', real(dtime)

! Execute on GPU: g_a = g_b + g_c
      call dtimer(dtime,itime,-1)
      call gpadd(g_a,g_b,g_c,nx)
      call dtimer(dtime,itime,1)
      write (*,*) 'GPU add time=', real(dtime)
!
! Check for correctness: compare a and g_a
! Copy g_a from GPU to c on Host, then compare a with c
      call gpu_fcopyout(c_loc(c(1)),g_a,nx)
      epsmax = 0.0
      do j = 1, nx
         eps = abs(a(j)-c(j))
         if (eps > epsmax) epsmax = eps
      enddo
      write (*,*) 'maximum difference = ', epsmax
!
! deallocate memory on GPU
      call gpu_deallocate(g_a,irc)
      call gpu_deallocate(g_b,irc)
      call gpu_deallocate(g_c,irc)
! close down GPU
      call end_cu()
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
      
