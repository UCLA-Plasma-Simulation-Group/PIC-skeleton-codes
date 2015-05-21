!-----------------------------------------------------------------------
      program gputest
! GPU vector add test program
! written by Viktor K. Decyk, UCLA
! nx = size of array, nblock = block size on GPU
      use gpumain_cuf
      implicit none
      integer :: nx = 1048576, nblock = 64
      integer :: j, irc
      real :: time = 0.0
      real :: eps, epsmax
      double precision :: dtime
      real, dimension(:), pointer :: a, b, c, d
      integer, dimension(4) :: itime
      real, device, dimension(:), allocatable :: g_a, g_b, g_c
!
! initialize Host data
      allocate(a(nx),b(nx),c(nx),d(nx))
      do j = 1, nx
         b(j) = j
         c(j) = 2*j
      enddo
      a = 0.0; d = -1.0
! set up GPU
      irc = 0
      call setgbsize(nblock)
      call init_cu(0,irc)
      if (irc /= 0) then
         write (*,*) 'CUDA initialization error!'
         stop
      endif
! allocate data on GPU, using CUDA functions
!     call g_fallocate(g_a,nx,irc)
!     call g_fallocate(g_b,nx,irc)
!     call g_fallocate(g_c,nx,irc)
!     if (irc /= 0) then
!        write (*,*) 'GPU allocate error!'
!        stop
!     endif
! allocate data on GPU, using Fortran90 array syntax
      allocate(g_a(nx),g_b(nx),g_c(nx))
!
! First execute on Host in Fortran: a = b + c
      call dtimer(dtime,itime,-1)
      call fadd(a,b,c,nx)
      call dtimer(dtime,itime,1)
      write (*,*) 'Fortran add time=', real(dtime)
! Copy data to GPU
      call dtimer(dtime,itime,-1)
! using CUDA functions
!     call copyin_gmemptr(b,g_b,nx)
!     call copyin_gmemptr(c,g_c,nx)
! using Fortran90 array syntax
      g_b = b; g_c = c
      call dtimer(dtime,itime,1)
      time = time + real(dtime)
      write (*,*) 'Copyin time=', real(dtime)
! Execute on GPU: g_a = g_b + g_c
      call dtimer(dtime,itime,-1)
      call gpadd(g_a,g_b,g_c,nx)
      call dtimer(dtime,itime,1)
      time = time + real(dtime)
      write (*,*) 'GPU add time=', real(dtime)
! Copy data from GPU: d = g_a
      call dtimer(dtime,itime,-1)
! using CUDA functions
!     call copyout_gmemptr(d,g_a,nx)
! using Fortran90 array syntax
      d = g_a
      call dtimer(dtime,itime,1)
      time = time + real(dtime)
      write (*,*) 'Copyout time=', real(dtime)
      write (*,*) 'Total GPU time=',time
!
! Check for correctness: compare a and d
      epsmax = 0.0
      do j = 1, nx
         eps = abs(a(j)-d(j))
         if (eps > epsmax) epsmax = eps
      enddo
      write (*,*) 'maximum difference = ', epsmax
!
! deallocate memory on GPU, using CUDA functions
!     call g_deallocate(g_a,irc)
!     call g_deallocate(g_b,irc)
!     call g_deallocate(g_c,irc)
! deallocate memory on GPU, using Fortran90 array syntax
      deallocate(g_a,g_b,g_c)
! close down GPU
      call end_cu()
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
      
