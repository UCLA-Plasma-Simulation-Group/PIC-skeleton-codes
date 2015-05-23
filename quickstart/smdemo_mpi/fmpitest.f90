!-----------------------------------------------------------------------
      program mpitest
! MPI vector add test program
! written by Viktor K. Decyk, UCLA
! nx = size of array
      implicit none
      integer :: nx = 1048576
! idproc = processor id
! nvp = number of processors obtained
      integer :: idproc, nvp
      integer :: j, joff, nxp, nxps, irc
      real :: time = 0.0
      real :: eps, epsmax
      real, dimension(1) :: epsm, work
! timing data
      double precision :: dtime
      integer, dimension(4) :: itime
! global data for Fortran
      real, dimension(:), pointer :: a, b, c
! local data for MPI
      real, dimension(:), pointer :: p_a, p_b, p_c
!
! initialize global Host data
      allocate(a(nx),b(nx),c(nx))
! initialize the same global vectors on each MPI node
      do j = 1, nx
         b(j) = j
         c(j) = 2*j
      enddo
      a = 0.0;
! set up MPI
! idproc contains local processor id, different for each MPI node
      irc = 0
      call init_mpi(idproc,nvp,irc)
      if (irc /= 0) then
         write (*,*) idproc, 'MPI initialization error!'
         stop
      endif
      if (idproc==0) write (*,*) 'MPI nodes available = ', nvp
!
! vector data will be partitioned into local arrays among MPI nodes
! nxp = maximum size of local array on each MPI node
      nxp = (nx - 1)/nvp + 1
! allocate local data for MPI
      allocate(p_a(nxp),p_b(nxp),p_c(nxp))
! joff = starting location in global array for local array
! different value on each MPI node, since it depends on idproc
      joff = nxp*idproc
! nxps = actual size of local array, if nx not an exact multiple of nvp
! possibly different value on each MPI node, since it depends on joff
      nxps = min(nxp,max(0,nx-joff))
! Copy part of initial vector data to local arrays for MPI
      do j = 1, nxps
         p_b(j) = b(j+joff)
         p_c(j) = c(j+joff)
      enddo
! alternative copy data for MPI:
! copy from MPI node 0 to other nodes if only node 0 has initial data
!     call vscatter(b,p_b,nx,nxp)
!     call vscatter(c,p_c,nx,nxp)
!
! First execute on Host in Fortran: a = b + c
! Each processor executes the same add of global data
      call dtimer(dtime,itime,-1)
      call fadd(a,b,c,nx)
      call dtimer(dtime,itime,1)
      if (idproc==0) write (*,*) 'Fortran add time=', real(dtime)

! Execute with MPI: p_a = p_b + p_c
! Each processor adds only its local data
      call dtimer(dtime,itime,-1)
      call mpadd(p_a,p_b,p_c,nxps)
      call dtimer(dtime,itime,1)
      time = time + real(dtime)
      if (idproc==0) write (*,*) 'MPI add time=', real(dtime)
!
! Check for correctness: compare a and p_a
! Each MPI node compares global result with local array result
      epsmax = 0.0
      do j = 1, nxps
         eps = abs(a(j+joff)-p_a(j))
         if (eps > epsmax) epsmax = eps
      enddo
! find global maximum error for each of the local maximum errors
      epsm(1) = epsmax
      call PPMAX(epsm,work,1)
      epsmax = epsm(1)
      if (idproc==0) write (*,*) 'maximum difference = ', epsmax
!
! Alternate check for correctness: compare a and p_a
! copy to c from other MPI nodes to node 0 if only node 0 has answer
!     call vgather(c,p_a,nx,nxp)
!     if (idproc==0) then
!        epsmax = 0.0
!        do j = 1, nx
!           eps = abs(a(j)-c(j))
!           if (eps > epsmax) epsmax = eps
!        enddo
!        write (*,*) 'maximum difference = ', epsmax
!     endif
!
! deallocate memory for MPI
      deallocate(p_a,p_b,p_c)
! close down MPI
      call end_mpi()
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
      
