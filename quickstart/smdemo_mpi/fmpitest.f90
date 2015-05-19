!-----------------------------------------------------------------------
      program mpitest
! mpi vector add test program
! written by Viktor K. Decyk, UCLA
! nx = size of array
      implicit none
      integer :: nx = 1048576
! idproc = processor id
! nvp = number of processors obtained
      integer :: idproc, nvp
      integer :: j, nxp, nxps, irc
      real :: time = 0.0
      real :: eps, epsmax
      double precision :: dtime
      real, dimension(:), pointer :: a, b, c, d
      integer, dimension(4) :: itime
      real, dimension(:), pointer :: p_a, p_b, p_c
!
! set up mpi
      irc = 0
      call init_mpi(idproc,nvp,irc)
      if (irc /= 0) then
         write (*,*) idproc, 'MPI initialization error!'
         stop
      endif
      if (idproc==0) write (*,*) 'mpi nodes available = ', nvp
! initialize Host data
      if (idproc==0) then
         allocate(a(nx),b(nx),c(nx),d(nx))
         do j = 1, nx
            b(j) = j
            c(j) = 2*j
         enddo
         a = 0.0; d = -1.0
      endif
! nxp = maximum size of array on each processor
      nxp = (nx - 1)/nvp + 1
! allocate data for mpi
      allocate(p_a(nxp),p_b(nxp),p_c(nxp))
! nxps = size of actual data, if nx is not an exact multiple of nvp
      nxps = min(nxp,max(0,nx-nxp*idproc))
!
! First execute on Host in Fortran: a = b + c
      if (idproc==0) then
         call dtimer(dtime,itime,-1)
         call fadd(a,b,c,nx)
         call dtimer(dtime,itime,1)
         write (*,*) 'Fortran add time=', real(dtime)
      endif
! Copy data for mpi
      call dtimer(dtime,itime,-1)
      call vscatter(b,p_b,nx,nxp)
      call vscatter(c,p_c,nx,nxp)
      call dtimer(dtime,itime,1)
      time = time + real(dtime)
      if (idproc==0) write (*,*) 'Copyin time=', real(dtime)
! Execute with mpi: g_a = g_b + g_c
      call dtimer(dtime,itime,-1)
      call mpadd(p_a,p_b,p_c,nxps)
      call dtimer(dtime,itime,1)
      time = time + real(dtime)
      if (idproc==0) write (*,*) 'mpi add time=', real(dtime)
! Copy data from mpi: d = g_a
      call dtimer(dtime,itime,-1)
      call vgather(d,p_a,nx,nxp)
      call dtimer(dtime,itime,1)
      time = time + real(dtime)
      if (idproc==0) then
         write (*,*) 'Copyout time=', real(dtime)
         write (*,*) 'Total mpi time=',time
      endif
!
! Check for correctness: compare a and d
      if (idproc==0) then
         epsmax = 0.0
         do j = 1, nx
            eps = abs(a(j)-d(j))
            if (eps > epsmax) epsmax = eps
         enddo
         write (*,*) 'maximum difference = ', epsmax
      endif
!
! deallocate memory for mpi
      deallocate(p_a,p_b,p_c)
! close down mpi
      call end_mpi()
! deallocate Host memory
      if (idproc==0) deallocate(a,b,c,d)
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
      
