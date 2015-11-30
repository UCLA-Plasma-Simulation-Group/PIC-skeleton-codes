!-----------------------------------------------------------------------
! Basic parallel PIC library for MPI communications with OpenMP
! pplib3.f90 contains basic communications procedures for 1d partitions:
! PPINIT2 initializes parallel processing for Fortran90, returns
!         number of processors and processor id.
! PPEXIT terminates parallel processing.
! PPABORT aborts parallel processing.
! PWTIMERA performs parallel local wall clock timing.
! PPSUM performs parallel sum of a real vector.
! PPDSUM performs parallel sum of a double precision vector.
! PPIMAX performs parallel maximum of an integer vector.
! PPDMAX performs parallel maximum of a double precision vector.
! PPNCGUARD32L copies data to guard cells in y and z for scalar data,
!              linear interpolation, and distributed data with 2D
!              non-uniform partition.
! PPNAGUARD32L adds guard cells in y and z for scalar array, linear
!              interpolation, and distributed data with 2D non-uniform
!              partition.
! PPNACGUARD32L adds guard cells in y and z for vector array, linear
!               interpolation, and distributed data with 2D non-uniform
!               partition.
! PPTPOS3A performs a transpose of a complex scalar array, distributed
!          in y and z, to a complex scalar array, distributed in x and z
! PPTPOS3B performs a transpose of a complex scalar array, distributed
!          in x and z, to a complex scalar array, distributed in x and y
! PPNTPOS3A performs a transpose of an n component complex vector array,
!           distributed in y and z, to an n component complex vector
!           array, distributed in x and z.
! PPNTPOS3B performs a transpose of an n component complex vector array,
!           distributed in x and z, to an n component complex vector
!           array, distributed in x and y.
! PPPMOVE32 moves particles in y/z into appropriate spatial regions for
!           tiled distributed data with 2D spatial decomposition
! written by viktor k. decyk, ucla
! copyright 1995, regents of the university of california
! update: october 27, 2015
      module mpplib3
      use mpi
      implicit none
!
! common data for parallel processing
! lstat = length of status array
      integer, parameter :: lstat = MPI_STATUS_SIZE
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! mreal = default datatype for reals
! mint = default datatype for integers
! mcplx = default datatype for complex type
! mdouble = default double precision type
! lworld = MPI_COMM_WORLD communicator
      integer :: nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! msum = MPI_SUM
! mmax = MPI_MAX
      integer :: msum, mmax
      save
!
      private
      public :: PPINIT2, PPEXIT, PPABORT, PWTIMERA
      public :: PPSUM, PPDSUM, PPIMAX, PPDMAX
      public :: PPNCGUARD32L, PPNAGUARD32L, PPNACGUARD32L
      public :: PPTPOS3A, PPTPOS3B, PPNTPOS3A, PPNTPOS3B, PPPMOVE32
!
      contains
!
!-----------------------------------------------------------------------
      subroutine PPINIT2(idproc,nvp)
! this subroutine initializes parallel processing
! lgrp communicator = MPI_COMM_WORLD
! output: idproc, nvp
! idproc = processor id in lgrp communicator
! nvp = number of real or virtual processors obtained
      implicit none
      integer, intent(inout) :: idproc, nvp
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! mreal = default datatype for reals
! mint = default datatype for integers
! mcplx = default datatype for complex type
! mdouble = default double precision type
! lworld = MPI_COMM_WORLD communicator
! msum = MPI_SUM
! mmax = MPI_MAX
! local data
      integer :: ierror, ndprec, idprec
      integer :: iprec
      logical :: flag
      real :: prec
! ndprec = (0,1) = (no,yes) use (normal,autodouble) precision
      if (digits(prec) > 24) then
         ndprec = 1
      else
         ndprec = 0
      endif
! idprec = (0,1) = (no,yes) use (normal,autodouble) integer precision
      if (digits(iprec) > 31) then
         idprec = 1
      else
         idprec = 0
      endif
! this segment is used for mpi computers
! indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (.not.flag) then
! initialize the MPI execution environment
         call MPI_INIT(ierror)
         if (ierror /= 0) stop
      endif
      lworld = MPI_COMM_WORLD
      lgrp = lworld
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierror)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nproc,ierror)
! set default datatypes
      mint = MPI_INTEGER
      mdouble = MPI_DOUBLE_PRECISION
! single precision real
      if (ndprec==0) then
         mreal = MPI_REAL
         mcplx = MPI_COMPLEX
! double precision real
      else
         mreal = MPI_DOUBLE_PRECISION
         mcplx = MPI_DOUBLE_COMPLEX
      endif
! single precision integer
!     if (idprec==0) then
!        mint = MPI_INTEGER
! double precision integer
!     else
!        mint = MPI_INTEGER8
!     endif
! operators
      msum = MPI_SUM
      mmax = MPI_MAX
      nvp = nproc
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPEXIT()
! this subroutine terminates parallel processing
      implicit none
! lworld = MPI_COMM_WORLD communicator
! local data
      integer :: ierror
      logical :: flag
! indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (flag) then
! synchronize processes
         call MPI_BARRIER(lworld,ierror)
! terminate MPI execution environment
         call MPI_FINALIZE(ierror)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPABORT()
! this subroutine aborts parallel processing
      implicit none
! lworld = MPI_COMM_WORLD communicator
! local data
      integer :: errorcode, ierror
      logical :: flag
! indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (flag) then
         errorcode = 1
! terminate MPI execution environment
         call MPI_ABORT(lworld,errorcode,ierror)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
! this subroutine performs local wall clock timing
! input: icntrl, dtime
! icntrl = (-1,0,1) = (initialize,ignore,read) clock
! clock should be initialized before it is read!
! time = elapsed time in seconds
! dtime = current time
! written for mpi
      implicit none
      integer, intent(in) :: icntrl
      real, intent(inout) :: time
      double precision, intent(inout) :: dtime
! local data
      double precision :: jclock
! initialize clock
      if (icntrl==(-1)) then
         dtime = MPI_WTIME()
! read clock and write time difference from last clock initialization
      else if (icntrl==1) then
         jclock = dtime
         dtime = MPI_WTIME()
         time = real(dtime - jclock)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPSUM(f,g,nxp)
! this subroutine performs a parallel sum of a vector, that is:
! f(j,k) = sum over k of f(j,k)
! at the end, all processors contain the same summation.
! f = input and output real data
! g = scratch real array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mreal = default datatype for reals
! msum = MPI_SUM
! local data
      integer :: j, ierr
! perform sum
      call MPI_ALLREDUCE(f,g,nxp,mreal,msum,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSUM(f,g,nxp)
! this subroutine performs a parallel sum of a vector, that is:
! f(j,k) = sum over k of f(j,k)
! at the end, all processors contain the same summation.
! f = input and output double precision data
! g = scratch double precision array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mdouble = default double precision type
! msum = MPI_SUM
! local data
      integer :: j, ierr
! perform sum
      call MPI_ALLREDUCE(f,g,nxp,mdouble,msum,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPIMAX(if,ig,nxp)
! this subroutine finds parallel maximum for each element of a vector
! that is, if(j,k) = maximum as a function of k of if(j,k)
! at the end, all processors contain the same maximum.
! if = input and output integer data
! ig = scratch integer array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if, ig
! lgrp = current communicator
! mint = default datatype for integers
! mmax = MPI_MAX
! local data
      integer :: j, ierr
! find maximum
      call MPI_ALLREDUCE(if,ig,nxp,mint,mmax,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         if(j) = ig(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDMAX(f,g,nxp)
! this subroutine finds parallel maximum for each element of a vector
! that is, f(j,k) = maximum as a function of k of f(j,k)
! at the end, all processors contain the same maximum.
! f = input and output double precision data
! g = scratch double precision array
! nxp = number of data values in vector
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
! lgrp = current communicator
! mdouble = default double precision type
! mmax = MPI_MAX
! local data
      integer j, ierr
! find maximum
      call MPI_ALLREDUCE(f,g,nxp,mdouble,mmax,lgrp,ierr)
! copy output from scratch array
      do j = 1, nxp
         f(j) = g(j)
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx&
     &,idds)
! this subroutine copies data to guard cells in non-uniform partitions
! f(j,k,l) = real data for grid j,k,l in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! scs(j,k) = scratch array for particle partition
! nyzp(1:2) = number of primary gridpoints in y/z in particle partition
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition
! linear interpolation, for distributed data,
! with 2D spatial decomposition
      implicit none
      integer, intent(in) :: kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, idds
      real, dimension(nxv,nypmx,nzpmx), intent(inout) :: f
      real, dimension(nxv,nzpmx,2), intent(inout) :: scs
      integer, dimension(idds), intent(in) :: nyzp
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: j, k, js, ks, noff, kr, kl
      integer :: nxvz, nxvzs, nxvy, nxvys, nyp1, nzp1
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      nyp1 = nyzp(1) + 1
      nzp1 = nyzp(2) + 1
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      noff = nypmx*nzpmx
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
! special case for one processor in y
      if (nvpy==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nyzp(2)
            do j = 1, nxv
               f(j,nyp1,k) = f(j,1,k)
            enddo
         enddo
!$OMP END PARALLEL DO
      else
! buffer data in y
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nyzp(2)
            do j = 1, nxv
               scs(j,k,1) = f(j,1,k)
            enddo
         enddo
!$OMP END PARALLEL DO
! copy to guard cells in y
         nxvzs = nxv*nyzp(2)
         kr = js + 1
         if (kr >= nvpy) kr = kr - nvpy
         kl = js - 1
         if (kl < 0) kl = kl + nvpy
         kr = kr + nvpy*ks
         kl = kl + nvpy*ks
! this segment is used for mpi computers
         call MPI_IRECV(scs(1,1,2),nxvz,mreal,kr,noff+3,lgrp,msid,ierr)
         call MPI_SEND(scs,nxvzs,mreal,kl,noff+3,lgrp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
! copy guard cells
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nyzp(2)
            do j = 1, nxv
               f(j,nyp1,k) = scs(j,k,2)
           enddo
         enddo
!$OMP END PARALLEL DO
      endif
! special case for one processor in z
      if (nvpz==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nyp1
            do j = 1, nxv
               f(j,k,nzp1) = f(j,k,1)
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! copy to guard cells in z
      nxvys = nxv*nyp1
      kr = ks + 1
      if (kr >= nvpz) kr = kr - nvpz
      kl = ks - 1
      if (kl < 0) kl = kl + nvpz
      kr = js + nvpy*kr
      kl = js + nvpy*kl
! this segment is used for mpi computers
      call MPI_IRECV(f(1,1,nzp1),nxvy,mreal,kr,noff+4,lgrp,msid,ierr)
      call MPI_SEND(f,nxvys,mreal,kl,noff+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,    &
     &nypmx,nzpmx,idds)
! this subroutine adds data from guard cells in non-uniform partitions
! f(j,k,l) = real data for grid j,k,l in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! scs/scr = scratch arrays for particle partition
! nyzp(1:2) = number of primary gridpoints in y/z in particle partition
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nx = system length in x direction
! nxv = first dimension of f, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition
! linear interpolation, for distributed data
! with 2D spatial decomposition
      implicit none
      integer, intent(in) :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
      integer, intent(in) :: idds
      real, dimension(nxv,nypmx,nzpmx), intent(inout) :: f
      real, dimension(nxv,nzpmx,2), intent(inout) :: scs
      real, dimension(nxv,nypmx), intent(inout) :: scr
      integer, dimension(idds), intent(in) :: nyzp
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: j, k, js, ks, noff, kr, kl
      integer :: nx1, nxvz, nxvzs, nxvy, nxvys, nyp1, nzp1
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      nx1 = nx + 1
      nyp1 = nyzp(1) + 1
      nzp1 = nyzp(2) + 1
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      noff = nypmx*nzpmx
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
! special case for one processor in y
      if (nvpy==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nzp1
            do j = 1, nx1
               f(j,1,k) = f(j,1,k) + f(j,nyp1,k)
               f(j,nyp1,k) = 0.0
            enddo
         enddo
!$OMP END PARALLEL DO
      else
! buffer data in y
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nzp1 
            do j = 1, nxv
               scs(j,k,1) = f(j,nyp1,k)
            enddo
         enddo
!$OMP END PARALLEL DO
! add guard cells in y
         nxvzs = nxv*nzp1
         kr = js + 1
         if (kr >= nvpy) kr = kr - nvpy
         kl = js - 1
         if (kl < 0) kl = kl + nvpy
         kr = kr + nvpy*ks
         kl = kl + nvpy*ks
! this segment is used for mpi computers
         call MPI_IRECV(scs(1,1,2),nxvz,mreal,kl,noff+1,lgrp,msid,ierr)
         call MPI_SEND(scs,nxvzs,mreal,kr,noff+1,lgrp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
! add up the guard cells
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nzp1
            do j = 1, nx1
               f(j,1,k) = f(j,1,k) + scs(j,k,2)
               f(j,nyp1,k) = 0.0
            enddo
         enddo
!$OMP END PARALLEL DO
      endif
! special case for one processor in z
      if (nvpz==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, nyp1
            do j = 1, nx1
               f(j,k,1) = f(j,k,1) + f(j,k,nzp1)
               f(j,k,nzp1) = 0.0
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! add guard cells in z
      nxvys = nxv*nyp1
      kr = ks + 1
      if (kr >= nvpz) kr = kr - nvpz
      kl = ks - 1
      if (kl < 0) kl = kl + nvpz
      kr = js + nvpy*kr
      kl = js + nvpy*kl
! this segment is used for mpi computers
      call MPI_IRECV(scr,nxvy,mreal,kl,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,nzp1),nxvys,mreal,kr,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
! add up the guard cells
!$OMP PARALLEL DO PRIVATE(j,k)
      do k = 1, nyp1
         do j = 1, nx1
            f(j,k,1) = f(j,k,1) + scr(j,k)
            f(j,k,nzp1) = 0.0
         enddo
      enddo
!$OMP END PARALLEL DO
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNACGUARD32L(f,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,  &
     &nxv,nypmx,nzpmx,idds)
! this subroutine adds data from guard cells in non-uniform partitions
! f(ndim,j,k,l) = real data for grid j,k,l in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! scs/scr = scratch arrays for particle partition
! nyzp(1:2) = number of primary gridpoints in y/z in particle partition
! ndim = leading dimension of array f
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nx = system length in x direction
! nxv = second dimension of f, must be >= nx+1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! idds = dimensionality of domain decomposition
! linear interpolation, for distributed data
! with 2D spatial decomposition
      implicit none
      integer, intent(in) :: ndim, kstrt, nvpy, nvpz, nx, nxv
      integer, intent(in) :: nypmx, nzpmx, idds
      real, dimension(ndim,nxv,nypmx,nzpmx), intent(inout) :: f
      real, dimension(ndim,nxv,nzpmx,2), intent(inout) :: scs
      real, dimension(ndim,nxv,nypmx), intent(inout) :: scr
      integer, dimension(idds), intent(in) :: nyzp
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: i, j, k, js, ks, noff, kr, kl
      integer :: nx1, nxvz, nxvzs, nxvy, nxvys, nyp1, nzp1
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      nx1 = nx + 1
      nyp1 = nyzp(1) + 1
      nzp1 = nyzp(2) + 1
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      noff = ndim*nypmx*nzpmx
      nxvz = ndim*nxv*nzpmx
      nxvy = ndim*nxv*nypmx
! special case for one processor in y
      if (nvpy==1) then
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, nzp1
            do j = 1, nx1
               do i = 1, ndim
                  f(i,j,1,k) = f(i,j,1,k) + f(i,j,nyp1,k)
                  f(i,j,nyp1,k) = 0.0
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
      else
! buffer data in y
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, nzp1 
            do j = 1, nxv
               do i = 1, ndim
                  scs(i,j,k,1) = f(i,j,nyp1,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
! add guard cells in y
         nxvzs = ndim*nxv*nzp1
         kr = js + 1
         if (kr >= nvpy) kr = kr - nvpy
         kl = js - 1
         if (kl < 0) kl = kl + nvpy
         kr = kr + nvpy*ks
         kl = kl + nvpy*ks
! this segment is used for mpi computers
         call MPI_IRECV(scs(1,1,1,2),nxvz,mreal,kl,noff+1,lgrp,msid,ierr&
     &)
         call MPI_SEND(scs,nxvzs,mreal,kr,noff+1,lgrp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
! add up the guard cells
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, nzp1
           do j = 1, nx1
              do i = 1, ndim
                  f(i,j,1,k) = f(i,j,1,k) + scs(i,j,k,2)
                  f(i,j,nyp1,k) = 0.0
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
      endif
! special case for one processor in z
      if (nvpz==1) then
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, nyp1
            do j = 1, nx1
               do i = 1, ndim
                  f(i,j,k,1) = f(i,j,k,1) + f(i,j,k,nzp1)
                  f(i,j,k,nzp1) = 0.0
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! add guard cells in z
      nxvys = ndim*nxv*nyp1
      kr = ks + 1
      if (kr >= nvpz) kr = kr - nvpz
      kl = ks - 1
      if (kl < 0) kl = kl + nvpz
      kr = js + nvpy*kr
      kl = js + nvpy*kl
! this segment is used for mpi computers
      call MPI_IRECV(scr,nxvy,mreal,kl,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,nzp1),nxvys,mreal,kr,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
! add up the guard cells
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = 1, nyp1
         do j = 1, nx1
            do i = 1, ndim
               f(i,j,k,1) = f(i,j,k,1) + scr(i,j,k)
               f(i,j,k,nzp1) = 0.0
            enddo
         enddo
      enddo
!$OMP END PARALLEL DO
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPTPOS3A(f,g,s,t,nx,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,nxv, &
     &nyv,kxypd,kypd,kzpd)
! this subroutine performs a transpose of a matrix f, distributed in y
! and z to a matrix g, distributed in x and z, that is,
! g(k+kyp*(m-1),j,l) = f(j+kxyp*(n-1),k,l), where
! 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
! 1 <= n <= nx/kxyp, 1 <= m <= ny/kyp
! and where indices n and m can be distributed across processors.
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! f = complex input array
! g = complex output array
! s, t = complex scratch arrays
! nx/ny/nz = number of points in x/y/z
! kxyp/kyp/kzp = number of data values per block in x/y/z
! kstrt = starting data block number
! nvpy = number of real or virtual processors in y
! nxv/nyv = first dimension of f/g
! kypd/kxypd = second dimension of f/g
! kzpd = third dimension of f and g
      implicit none
      integer, intent(in) :: nx, ny, nz, kxyp, kyp, kzp, kstrt, nvpy
      integer, intent(in) :: nxv, nyv, kxypd, kypd, kzpd
      complex, dimension(nxv,kypd,kzpd), intent(in) :: f
      complex, dimension(nyv,kxypd,kzpd), intent(inout) :: g
      complex, dimension(kxyp*kyp*kzp), intent(inout) :: s, t
! lgrp = current communicator
! mcplx = default datatype for complex
! local data
      integer :: n, j, k, l, js, ks, kxyps, kyps, kzps, id, joff, koff
      integer :: ld, jd, kxyzp
      integer :: ierr, msid, ll
      integer, dimension(lstat) :: istatus
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxyps = min(kxyp,max(0,nx-kxyp*js))
      kyps = min(kyp,max(0,ny-kyp*js))
      kzps = min(kzp,max(0,nz-kzp*ks))
      kxyzp = kxyp*kyp*kzp
! special case for one processor
      if (nvpy==1) then
!$OMP PARALLEL DO PRIVATE(j,k,l,ll)
         do ll = 1, kyp*kzp
            l = (ll - 1)/kyp
            k = ll - kyp*l
            l = l + 1
            do j = 1, kxyp
               g(k,j,l) = f(j,k,l)
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! this segment is used for shared memory computers
!     do l = 1, nz
!        do m = 1, min(ny,nvpy)
!           koff = kyp*(m - 1)
!           do i = 1, min(nx,nvpy)
!           joff = kxyp*(i - 1)
!              do k = 1, min(kyp,max(0,ny-koff))
!                 do j = 1, min(kxyp,max(0,nx-joff))
!                    g(k+koff,j+joff,l) = f(j+joff,k+koff,l)
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
! this segment is used for mpi computers
      do n = 1, nvpy
         id = n - js - 1
         if (id < 0) id = id + nvpy
! extract data to send
         joff = kxyp*id
         ld = min(kxyp,max(0,nx-joff))
!$OMP PARALLEL DO PRIVATE(j,k,l,ll,koff)
         do ll = 1, kyps*kzps
            l = (ll - 1)/kyps
            k = ll - kyps*l
            l = l + 1
            koff = kyps*(l - 1) - 1
            do j = 1, ld
               s(j+ld*(k+koff)) = f(j+joff,k,l)
            enddo
         enddo
!$OMP END PARALLEL DO
         jd = id + nvpy*ks
         ld = ld*kyps*kzps
! post receive
         call MPI_IRECV(t,kxyzp,mcplx,jd,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mcplx,jd,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         koff = kyp*id
         ld = min(kyp,max(0,ny-koff))
!$OMP PARALLEL DO PRIVATE(j,k,l,ll,joff)
         do ll = 1, ld*kzps
            l = (ll - 1)/ld
            k = ll - ld*l
            l = l + 1
            joff = ld*(l - 1) - 1
            do j = 1, kxyps
               g(k+koff,j,l) = t(j+kxyps*(k+joff))
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPTPOS3B(g,h,s,t,nx,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,nvpz&
     &,nyv,nzv,kxypd,kyzpd,kzpd)
! this subroutine performs a transpose of a matrix g, distributed in x
! and z to a matrix h, distributed in x and y, that is,
! h(l+kzp*(n-1),j,k) = g(k+kyzp*(m-1),j,l), where
! 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
! 1 <= m <= ny/kyzp, 1 <= n <= nz/kzp
! and where indices n and m can be distributed across processors.
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! g = complex input array
! h = complex output array
! s, t = complex scratch arrays
! nx/ny/nz = number of points in x/y/z
! kxyp/kyzp/kzp = number of data values per block in x/y/z
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! nyv/nzv = first dimension of g/h
! kxypd = second dimension of g and h
! kzpd/kyzpd = third dimension of g/h
      implicit none
      integer, intent(in) :: nx, ny, nz, kxyp, kyzp, kzp, kstrt
      integer, intent(in) :: nvpy, nvpz, nyv, nzv, kxypd, kyzpd, kzpd
      complex, dimension(nyv,kxypd,kzpd), intent(inout) :: g
      complex, dimension(nzv,kxypd,kyzpd), intent(inout) :: h
      complex, dimension(kyzp*kxyp*kzp), intent(inout) :: s, t
! lgrp = current communicator
! mcplx = default datatype for complex
! local data
      integer :: n, j, k, l, js, ks, kxyps, kyzps, kzps, id, koff, loff
      integer :: ld, jd, kxyzp, ll
      integer :: ierr, msid
      integer, dimension(lstat) :: istatus
! js/ks = processor co-ordinates in x/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxyps = min(kxyp,max(0,nx-kxyp*js))
      kyzps = min(kyzp,max(0,ny-kyzp*ks))
      kzps = min(kzp,max(0,nz-kzp*ks))
      kxyzp = kxyp*kyzp*kzp
! special case for one processor
      if (nvpz==1) then
!$OMP PARALLEL DO PRIVATE(j,k,l,ll)
         do ll = 1, kxyp*kzp
            l = (ll - 1)/kxyp
            j = ll - kxyp*l
            l = l + 1
            do k = 1, kyzp
               h(l,j,k) = g(k,j,l)
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! this segment is used for shared memory computers
!     do i = 1, min(nz,nvpz)
!        loff = kzp*(i - 1)
!        do m = 1, min(ny,nvpz)
!           koff = kyzp*(m - 1)
!           do l = 1, min(kzp,max(0,nz-loff))
!              do j = 1, nx
!                 do k = 1, min(kyzp,max(0,ny-koff))
!                    h(l+loff,j,k+koff) = g(k+koff,j,l+loff)
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
! this segment is used for mpi computers
      do n = 1, nvpz
         id = n - ks - 1
         if (id < 0) id = id + nvpz
! extract data to send
         koff = kyzp*id
         ld = min(kyzp,max(0,ny-koff))
!$OMP PARALLEL DO PRIVATE(j,k,l,ll,loff)
         do ll = 1, kxyps*kzps
            l = (ll - 1)/kxyps
            j = ll - kxyps*l
            l = l + 1
            loff = kxyps*(l - 1) - 1
            do k = 1, ld
               s(k+ld*(j+loff)) = g(k+koff,j,l)
            enddo
         enddo
!$OMP END PARALLEL DO
         jd = js + nvpy*id
         ld = ld*kxyps*kzps
! post receive
         call MPI_IRECV(t,kxyzp,mcplx,jd,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mcplx,jd,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         loff = kzp*id
         ld = min(kzp,max(0,nz-loff))
!$OMP PARALLEL DO PRIVATE(j,k,l,ll,koff)
         do ll = 1, kxyps*ld
            l = (ll - 1)/kxyps
            j = ll - kxyps*l
            l = l + 1
            koff = kxyps*(l - 1) - 1
            do k = 1, kyzps
               h(l+loff,j,k) = t(k+kyzps*(j+koff))
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNTPOS3A(f,g,s,t,nx,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,ndim&
     &,nxv,nyv,kxypd,kypd,kzpd)
! this subroutine performs a transpose of a matrix f, distributed in y
! and z to a matrix g, distributed in x and z, that is,
! g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxyp*(n-1),k,l), where
! 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
! 1 <= n <= nx/kxyp, 1 <= m <= ny/kyp
! and where indices n and m can be distributed across processors.
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! f = complex input array
! g = complex output array
! s, t = complex scratch arrays
! nx/ny/nz = number of points in x/y/z
! kxyp/kyp/kzp = number of data values per block in x/y/z
! kstrt = starting data block number
! nvpy = number of real or virtual processors in y
! ndim = leading dimension of arrays f and g
! nxv/nyv = first dimension of f/g
! kypd/kxypd = second dimension of f/g
! kzpd = third dimension of f and g
      implicit none
      integer, intent(in) :: nx, ny, nz, kxyp, kyp, kzp, kstrt, nvpy
      integer, intent(in) :: ndim, nxv, nyv, kxypd, kypd, kzpd
      complex, dimension(ndim,nxv,kypd,kzpd), intent(in) :: f
      complex, dimension(ndim,nyv,kxypd,kzpd), intent(inout) :: g
      complex, dimension(ndim,kxyp*kyp*kzp), intent(inout) :: s, t
! lgrp = current communicator
! mcplx = default datatype for complex
! local data
      integer :: i, n, j, k, l, js, ks, kxyps, kyps, kzps, id
      integer :: joff, koff, ld, jd, kxyzp, ll
      integer :: ierr, msid
      integer, dimension(lstat) :: istatus
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxyps = min(kxyp,max(0,nx-kxyp*js))
      kyps = min(kyp,max(0,ny-kyp*js))
      kzps = min(kzp,max(0,nz-kzp*ks))
      kxyzp = ndim*kxyp*kyp*kzp
! special case for one processor
      if (nvpy==1) then
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll)
         do ll = 1, kyp*kzp
            l = (ll - 1)/kyp
            k = ll - kyp*l
            l = l + 1
            do j = 1, kxyp
               do i = 1, ndim
                  g(i,k,j,l) = f(i,j,k,l)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! this segment is used for shared memory computers
!     do l = 1, nz
!        do m = 1, min(ny,nvpy)
!           koff = kyp*(m - 1)
!           do i = 1, min(nx,nvpy)
!              joff = kxyp*(i - 1)
!              do k = 1, min(kyp,max(0,ny-koff))
!                 do j = 1, min(kxyp,max(0,nx-joff))
!                    do n = 1, ndim
!                       g(n,k+koff,j+joff,l) = f(n,j+joff,k+koff,l)
!                    enddo
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
! this segment is used for mpi computers
      do n = 1, nvpy
         id = n - js - 1
         if (id < 0) id = id + nvpy
! extract data to send
         joff = kxyp*id
         ld = min(kxyp,max(0,nx-joff))
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,koff)
         do ll = 1, kyps*kzps
            l = (ll - 1)/kyps
            k = ll - kyps*l
            l = l + 1
            koff = kyps*(l - 1) - 1
            do j = 1, ld
               do i = 1, ndim
                  s(i,j+ld*(k+koff)) = f(i,j+joff,k,l)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         jd = id + nvpy*ks
         ld = ndim*ld*kyps*kzps
! post receive
         call MPI_IRECV(t,kxyzp,mcplx,jd,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mcplx,jd,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         koff = kyp*id
         ld = min(kyp,max(0,ny-koff))
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,joff)
         do ll = 1, ld*kzps
            l = (ll - 1)/ld
            k = ll - ld*l
            l = l + 1
            joff = ld*(l - 1) - 1
            do j = 1, kxyps
               do i = 1, ndim
                  g(i,k+koff,j,l) = t(i,j+kxyps*(k+joff))
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNTPOS3B(g,h,s,t,nx,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,   &
     &nvpz,ndim,nyv,nzv,kxypd,kyzpd,kzpd)
! this subroutine performs a transpose of a matrix g, distributed in x
! and z to a matrix h, distributed in x and y, that is,
! h(1:ndim,l+kzp*(n-1),j,k) = g(1:ndim,k+kyzp*(m-1),j,l), where
! 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
! 1 <= m <= ny/kyzp, 1 <= n <= nz/kzp
! and where indices n and m can be distributed across processors.
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! g = complex input array
! h = complex output array
! s, t = complex scratch arrays
! nx/ny/nz = number of points in x/y/z
! kxyp/kyzp/kzp = number of data values per block in x/y/z
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! ndim = leading dimension of arrays g and h
! nyv/nzv = first dimension of g/h
! kxypd = second dimension of g and h
! kzpd/kyzpd = third dimension of g/h
      implicit none
      integer, intent(in) :: nx, ny, nz, kxyp, kyzp, kzp, kstrt, nvpy
      integer, intent(in) :: nvpz, ndim, nyv, nzv, kxypd, kyzpd, kzpd
      complex, dimension(ndim,nyv,kxypd,kzpd), intent(inout) :: g
      complex, dimension(ndim,nzv,kxypd,kyzpd), intent(inout) :: h
      complex, dimension(ndim,kyzp*kxyp*kzp), intent(inout) :: s, t
! lgrp = current communicator
! mcplx = default datatype for complex
! local data
      integer :: i, n, j, k, l, js, ks, kxyps, kyzps, kzps, id
      integer :: koff, loff, ld, jd, kxyzp, ll
      integer :: ierr, msid
      integer, dimension(lstat) :: istatus
! js/ks = processor co-ordinates in x/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxyps = min(kxyp,max(0,nx-kxyp*js))
      kyzps = min(kyzp,max(0,ny-kyzp*ks))
      kzps = min(kzp,max(0,nz-kzp*ks))
      kxyzp = ndim*kxyp*kyzp*kzp
! special case for one processor
      if (nvpz==1) then
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll)
         do ll = 1, kxyp*kzp
            l = (ll - 1)/kxyp
            j = ll - kxyp*l
            l = l + 1
            do k = 1, kyzp
               do i = 1, ndim
                  h(i,l,j,k) = g(i,k,j,l)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! this segment is used for shared memory computers
!     do i = 1, min(nz,nvpz)
!        loff = kzp*(i - 1)
!        do m = 1, min(ny,nvpz)
!           koff = kyzp*(m - 1)
!           do l = 1, min(kzp,max(0,nz-loff))
!              do j = 1, nx
!                 do k = 1, min(kyzp,max(0,ny-koff))
!                    do n = 1, ndim
!                       h(n,l+loff,j,k+koff) = g(n,k+koff,j,l+loff)
!                    enddo
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
! this segment is used for mpi computers
      do n = 1, nvpz
         id = n - ks - 1
         if (id < 0) id = id + nvpz
! extract data to send
         koff = kyzp*id
         ld = min(kyzp,max(0,ny-koff))
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,loff)
         do ll = 1, kxyps*kzps
            l = (ll - 1)/kxyps
            j = ll - kxyps*l
            l = l + 1
            loff = kxyps*(l - 1) - 1
            do k = 1, ld
               do i = 1, ndim
                  s(i,k+ld*(j+loff)) = g(i,k+koff,j,l)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         jd = js + nvpy*id
         ld = ndim*ld*kxyps*kzps
! post receive
         call MPI_IRECV(t,kxyzp,mcplx,jd,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mcplx,jd,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         loff = kzp*id
         ld = min(kzp,max(0,nz-loff))
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,koff)
         do ll = 1, kxyps*ld
            l = (ll - 1)/kxyps
            j = ll - kxyps*l
            l = l + 1
            koff = kxyps*(l - 1) - 1
            do k = 1, kyzps
               do i = 1, ndim
                  h(i,l+loff,j,k) = t(i,k+kyzps*(j+koff))
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPMOVE32(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr, &
     &mcls,kstrt,nvpy,nvpz,idimp,nbmax,mx1,myp1,mzp1,mxzyp1,irc)
! this subroutine moves particles into appropriate spatial regions in y
! for distributed data, with 2d domain decomposition in y/z.
! tiles are assumed to be arranged in 3D linear memory
! output: rbufr, rbufl, mcll, mclr
! sbufl = buffer for particles being sent to lower/back processor
! sbufr = buffer for particles being sent to upper/forward processor
! rbufl = buffer for particles being received from lower/back processor
! rbufr = buffer for particles being received from upper/forward
! processor
! ncll = particle number offsets sent to lower/back processor
! nclr = particle number offsets sent to upper/forward processor
! mcll = particle number offsets received from lower/back processor
! mclr = particle number offsets received from upper/forward processor
! mcls = particle number ofsets received from corner processors
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! idimp = size of phase space = 4 or 5
! nbmax =  size of buffers for passing particles between processors
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! mzp1 = (partition length in z direction - 1)/mz + 1
! mxzyp1 = mx1*max(myp1,mzp1)
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer, intent(in) :: kstrt, nvpy, nvpz, idimp, nbmax
      integer, intent(in) :: mx1, myp1, mzp1, mxzyp1
      integer, intent(inout) :: irc
      real, dimension(idimp,nbmax,2), intent(in) :: sbufr, sbufl
      real, dimension(idimp,nbmax,2), intent(inout) :: rbufr, rbufl
      integer, dimension(3,mxzyp1,3,2), intent(inout) :: ncll, nclr
      integer, dimension(3,mxzyp1,3,2), intent(inout) :: mcll, mclr
      integer, dimension(3,mx1+1,4), intent(inout) :: mcls
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: ierr, js, ks, kl, kr, i, j, k, n, jsl, jsr
      integer :: m, ll, lr, krr, krl, klr, kll, jsrr, jsrl, jslr, jsll
      integer :: nr, nl, mr, ml, nbr, nbl
      integer :: mxyp1, mxzp1, nbsize, ncsize, nsize
      integer, dimension(8) :: msid
      integer, dimension(12) :: itg
      integer, dimension(lstat) :: istatus
      integer, dimension(1) :: nb, iwork
      data itg /3,4,5,6,7,8,9,10,11,12,13,14/
! js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      mxyp1 = mx1*myp1
      mxzp1 = mx1*mzp1
      nbsize = idimp*nbmax
      ncsize = 9*mxzyp1
! copy particle buffers in y:
! update rbufl(:,1), rbufr(:,1), mcll(:,1), mclr(:,1)
! special case for one processor
      if ((nvpy*nvpz)==1) then
!$OMP PARALLEL
!$OMP DO PRIVATE(i,j,k,ll)
         do ll = 1, 3*mxzp1
            k = (ll - 1)/mxzp1
            j = ll - mxzp1*k
            k = k + 1
            do i = 1, 3
               mcll(i,j,k,1) = nclr(i,j,k,1)
               mclr(i,j,k,1) = ncll(i,j,k,1)
            enddo
         enddo
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i,j)
         do j = 1, nclr(3,mxzp1,3,1)
            do i = 1, idimp
               rbufl(i,j,1) = sbufr(i,j,1)
            enddo
         enddo
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i,j)
         do j = 1, ncll(3,mxzp1,3,1)
            do i = 1, idimp
               rbufr(i,j,1) = sbufl(i,j,1)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
! get particles from corners
         n = mx1*(mzp1 - 1)
! zero out base addresses in prefix scans
         if (n.gt.0) then
            nr = nclr(3,n,3,1)
            nl = ncll(3,n,3,1)
         else
            nr = nclr(3,mxzp1,2,1)
            nl = ncll(3,mxzp1,2,1)
         endif
         do j = 1, mx1
            do i = 1, 3
               nclr(i,j,2,1) = nclr(i,j,2,1) - nclr(3,mxzp1,1,1)
               nclr(i,n+j,3,1) = nclr(i,n+j,3,1) - nr
               ncll(i,j,2,1) = ncll(i,j,2,1) - ncll(3,mxzp1,1,1)
               ncll(i,n+j,3,1) = ncll(i,n+j,3,1) - nl
            enddo
         enddo
! add new base addresses in prefix scans
         ml = mcll(3,mxzp1,3,1)
         mr = mclr(3,mxzp1,3,1)
         do j = 1, mx1
            do i = 1, 3
               mcls(i,j,1) = nclr(i,j,2,1) + ml
               mcls(i,j,3) = ncll(i,j,2,1) + mr
            enddo
         enddo
         mcls(1,mx1+1,1) = ml
         mcls(1,mx1+1,3) = mr
! append corner particles to end of buffers
         k = nclr(3,mx1,2,1)
         m = nclr(3,mxzp1,1,1)
         do j = 1, k
            do i = 1, idimp
               rbufl(i,j+ml,1) = sbufr(i,j+m,1)
            enddo
         enddo
         ml = ml + k
         k = ncll(3,mx1,2,1)
         m = ncll(3,mxzp1,1,1)
         do j = 1, k
            do i = 1, idimp
               rbufr(i,j+mr,1) = sbufl(i,j+m,1)
            enddo
         enddo
         mr = mr + k
! add new base addresses in prefix scans
         do j = 1, mx1
            do i = 1, 3
               mcls(i,j,2) = nclr(i,n+j,3,1) + ml
               mcls(i,j,4) = ncll(i,n+j,3,1) + mr
            enddo
         enddo
         mcls(1,mx1+1,2) = ml
         mcls(1,mx1+1,4) = mr
! append more corner particles to end of buffers
         do j = 1, nclr(3,n+mx1,3,1)
            do i = 1, idimp
               rbufl(i,j+ml,1) = sbufr(i,j+nr,1)
            enddo
         enddo
         do j = 1, ncll(3,n+mx1,3,1)
            do i = 1, idimp
               rbufr(i,j+mr,1) = sbufl(i,j+nl,1)
            enddo
         enddo
! this segment is used for mpi computers
      else
! get particles from below and above
         kr = js + 1
         if (kr.ge.nvpy) kr = kr - nvpy
         kl = js - 1
         if (kl.lt.0) kl = kl + nvpy
         kr = kr + nvpy*ks
         kl = kl + nvpy*ks
! post receives
         call MPI_IRECV(mcll(1,1,1,1),ncsize,mint,kl,itg(1),lgrp,msid(1)&
     &,ierr)
         call MPI_IRECV(mclr(1,1,1,1),ncsize,mint,kr,itg(2),lgrp,msid(2)&
     &,ierr)
         call MPI_IRECV(rbufl(1,1,1),nbsize,mreal,kl,itg(3),lgrp,msid(3)&
     &,ierr)
         call MPI_IRECV(rbufr(1,1,1),nbsize,mreal,kr,itg(4),lgrp,msid(4)&
     &,ierr)
! send particle number offsets
         call MPI_ISEND(nclr(1,1,1,1),ncsize,mint,kr,itg(1),lgrp,msid(5)&
     &,ierr)
         call MPI_ISEND(ncll(1,1,1,1),ncsize,mint,kl,itg(2),lgrp,msid(6)&
     &,ierr)
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_WAIT(msid(2),istatus,ierr)
! send particles
         jsr = idimp*nclr(3,mxzp1,3,1)
         call MPI_ISEND(sbufr(1,1,1),jsr,mreal,kr,itg(3),lgrp,msid(7),  &
     &ierr)
         jsl = idimp*ncll(3,mxzp1,3,1)
         call MPI_ISEND(sbufl(1,1,1),jsl,mreal,kl,itg(4),lgrp,msid(8),  &
     &ierr)
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
! make sure sbufr, sbufl, ncll, and nclr have been sent
         do i = 1, 4
            call MPI_WAIT(msid(i+4),istatus,ierr)
         enddo
! get particles from corners
         kr = js + 1
         if (kr.ge.nvpy) kr = kr - nvpy
         kl = js - 1
         if (kl.lt.0) kl = kl + nvpy
         lr = ks + 1
         if (lr.ge.nvpz) lr = lr - nvpz
         ll = ks - 1
         if (ll.lt.0) ll = ll + nvpz
         krl = kr + nvpy*ll
         krr = kr + nvpy*lr
         kll = kl + nvpy*ll
         klr = kl + nvpy*lr
         nsize = 3*mx1
         n = mx1*(mzp1 - 1)
! zero out base addresses in prefix scans
         if (n.gt.0) then
            nr = nclr(3,n,3,1)
            nl = ncll(3,n,3,1)
         else
            nr = nclr(3,mxzp1,2,1)
            nl = ncll(3,mxzp1,2,1)
         endif
         do j = 1, mx1
            do i = 1, 3
               nclr(i,j,2,1) = nclr(i,j,2,1) - nclr(3,mxzp1,1,1)
               nclr(i,n+j,3,1) = nclr(i,n+j,3,1) - nr
               ncll(i,j,2,1) = ncll(i,j,2,1) - ncll(3,mxzp1,1,1)
               ncll(i,n+j,3,1) = ncll(i,n+j,3,1) - nl
            enddo
         enddo
         n = n + 1
! post receives
         call MPI_IRECV(mcls(1,1,1),nsize,mint,klr,itg(5),lgrp,msid(1), &
     &ierr)
         call MPI_IRECV(mcls(1,1,2),nsize,mint,kll,itg(6),lgrp,msid(2), &
     &ierr)
         call MPI_IRECV(mcls(1,1,3),nsize,mint,krr,itg(7),lgrp,msid(3), &
     &ierr)
         call MPI_IRECV(mcls(1,1,4),nsize,mint,krl,itg(8),lgrp,msid(4), &
     &ierr)
! send particle number offsets
         call MPI_ISEND(nclr(1,1,2,1),nsize,mint,krl,itg(5),lgrp,msid(5)&
     &,ierr)
         call MPI_ISEND(nclr(1,n,3,1),nsize,mint,krr,itg(6),lgrp,msid(6)&
     &,ierr)
         call MPI_ISEND(ncll(1,1,2,1),nsize,mint,kll,itg(7),lgrp,msid(7)&
     &,ierr)
         call MPI_ISEND(ncll(1,n,3,1),nsize,mint,klr,itg(8),lgrp,msid(8)&
     &,ierr)
! make sure particle offsets have been sent to and received from corners
         do i = 1, 8
            call MPI_WAIT(msid(i),istatus,ierr)
         enddo
! check for overflow errors
         ml = mcll(3,mxzp1,3,1)
         mr = mclr(3,mxzp1,3,1)
         nbl = nbmax - (ml + (mcls(3,mx1,1) + mcls(3,mx1,2)))
         nbr = nbmax - (mr + (mcls(2,mx1,3) + mcls(3,mx1,4)))
         nb(1) = min(-nbl,-nbr)
         call PPIMAX(nb,iwork,1)
         if (nb(1) > 0) then
            write (2,*) 'corner buffer overflow error = ', nb(1)
            irc = nb(1)
            return
         endif
         nbl = idimp*nbl
         nbr = idimp*nbr
! add new base addresses in prefix scans
         do j = 1, mx1
            do i = 1, 3
               mcls(i,j,1) = mcls(i,j,1) + ml
               mcls(i,j,3) = mcls(i,j,3) + mr
            enddo
         enddo
         mcls(1,mx1+1,1) = ml
         mcls(1,mx1+1,3) = mr
! post first part of particle receives, append to end
         ml = ml + 1
         call MPI_IRECV(rbufl(1,ml,1),nbl,mreal,klr,itg(9),lgrp,msid(1),&
     &ierr)
         mr = mr + 1
         call MPI_IRECV(rbufr(1,mr,1),nbr,mreal,krr,itg(11),lgrp,msid(3)&
     &,ierr)
! send first part of particles
         m = nclr(3,mxzp1,1,1) + 1
         jsrl = idimp*nclr(3,mx1,2,1)
         call MPI_ISEND(sbufr(1,m,1),jsrl,mreal,krl,itg(9),lgrp,msid(5),&
     &ierr)
         m = ncll(3,mxzp1,1,1) + 1
         jsll = idimp*ncll(3,mx1,2,1)
         call MPI_ISEND(sbufl(1,m,1),jsll,mreal,kll,itg(11),lgrp,msid(7)&
     &,ierr)
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,m,ierr)
         nbl = nbl - m
         m = m/idimp
         ml = ml + m - 1
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,m,ierr)
         nbr = nbr - idimp*m
         m = m/idimp
         mr = mr + m - 1
! add new base addresses in prefix scans
         do j = 1, mx1
            do i = 1, 3
               mcls(i,j,2) = mcls(i,j,2) + ml
               mcls(i,j,4) = mcls(i,j,4) + mr
            enddo
         enddo
         mcls(1,mx1+1,2) = ml
         mcls(1,mx1+1,4) = mr
! post second part of particle receives, append to end
         ml = ml + 1
         call MPI_IRECV(rbufl(1,ml,1),nbl,mreal,kll,itg(10),lgrp,msid(2)&
     &,ierr)
         mr = mr + 1
         call MPI_IRECV(rbufr(1,mr,1),nbr,mreal,krl,itg(12),lgrp,msid(4)&
     &,ierr)
! send second part of particles
         jsrr = idimp*nclr(3,n+mx1-1,3,1)
         m = nr + 1
         call MPI_ISEND(sbufr(1,m,1),jsrr,mreal,krr,itg(10),lgrp,msid(6)&
     &,ierr)
         jslr = idimp*ncll(3,n+mx1-1,3,1)
         m = nl + 1
         call MPI_ISEND(sbufl(1,m,1),jslr,mreal,klr,itg(12),lgrp,msid(8)&
     &,ierr)
         call MPI_WAIT(msid(2),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
! make sure sbufr and sbufl have been sent to corners
         do i = 1, 4
            call MPI_WAIT(msid(i+4),istatus,ierr)
         enddo
      endif
! copy particle buffers in z:
! update rbufl(:,2), rbufr(:,2), mcll(:,2), mclr(:,2)
! special case for one processor
      if ((nvpy*nvpz)==1) then
!$OMP PARALLEL
!$OMP DO PRIVATE(i,j,k,ll)
         do ll = 1, 3*mxyp1
            k = (ll - 1)/mxyp1
            j = ll - mxyp1*k
            k = k + 1
            do i = 1, 3
               mcll(i,j,k,2) = nclr(i,j,k,2)
               mclr(i,j,k,2) = ncll(i,j,k,2)
            enddo
         enddo
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i,j)
         do j = 1, nclr(3,mxyp1,3,2)
            do i = 1, idimp
               rbufl(i,j,2) = sbufr(i,j,2)
            enddo
         enddo
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i,j)
         do j = 1, ncll(3,mxyp1,3,2)
            do i = 1, idimp
               rbufr(i,j,2) = sbufl(i,j,2)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
! this segment is used for mpi computers
      else
! get particles from back and front
         kr = ks + 1
         if (kr.ge.nvpz) kr = kr - nvpz
         kl = ks - 1
         if (kl.lt.0) kl = kl + nvpz
         kr = js + nvpy*kr
         kl = js + nvpy*kl
! post receives
         call MPI_IRECV(mcll(1,1,1,2),ncsize,mint,kl,itg(1),lgrp,msid(1)&
     &,ierr)
         call MPI_IRECV(mclr(1,1,1,2),ncsize,mint,kr,itg(2),lgrp,msid(2)&
     &,ierr)
         call MPI_IRECV(rbufl(1,1,2),nbsize,mreal,kl,itg(3),lgrp,msid(3)&
     &,ierr)
         call MPI_IRECV(rbufr(1,1,2),nbsize,mreal,kr,itg(4),lgrp,msid(4)&
     &,ierr)
! send particle number offsets
         call MPI_ISEND(nclr(1,1,1,2),ncsize,mint,kr,itg(1),lgrp,msid(5)&
     &,ierr)
         call MPI_ISEND(ncll(1,1,1,2),ncsize,mint,kl,itg(2),lgrp,msid(6)&
     &,ierr)
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_WAIT(msid(2),istatus,ierr)
! send particles
         jsr = idimp*nclr(3,mxyp1,3,2)
         call MPI_ISEND(sbufr(1,1,2),jsr,mreal,kr,itg(3),lgrp,msid(7),  &
     &ierr)
         jsl = idimp*ncll(3,mxyp1,3,2)
         call MPI_ISEND(sbufl(1,1,2),jsl,mreal,kl,itg(4),lgrp,msid(8),  &
     &ierr)
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
! make sure sbufr, sbufl, ncll, and nclr have been sent
         do i = 1, 4
            call MPI_WAIT(msid(i+4),istatus,ierr)
         enddo
      endif
      end subroutine
!
      end module
!
! Make functions callable by Fortran77
!
!-----------------------------------------------------------------------
      subroutine PPINIT2(idproc,nvp)
      use mpplib3, only: SUB => PPINIT2
      implicit none
      integer, intent(inout) :: idproc, nvp
      call SUB(idproc,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPEXIT
      use mpplib3, only: SUB => PPEXIT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPABORT
      use mpplib3, only: SUB => PPABORT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
      use mpplib3, only: SUB => PWTIMERA
      implicit none
      integer, intent(in) :: icntrl
      real, intent(inout) :: time
      double precision, intent(inout) :: dtime
      call SUB(icntrl,time,dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPSUM(f,g,nxp)
      use mpplib3, only: SUB => PPSUM
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSUM(f,g,nxp)
      use mpplib3, only: SUB => PPDSUM
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPIMAX(if,ig,nxp)
      use mpplib3, only: SUB => PPIMAX
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if, ig
      call SUB(if,ig,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDMAX(f,g,nxp)
      use mpplib3, only: SUB => PPDMAX
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx&
     &,idds)
      use mpplib3, only: SUB => PPNCGUARD32L
      implicit none
      integer, intent(in) :: kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, idds
      real, dimension(nxv,nypmx,nzpmx), intent(inout) :: f
      real, dimension(nxv,nzpmx,2), intent(inout) :: scs
      integer, dimension(idds), intent(in) :: nyzp
      call SUB(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,idds)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,    &
     &nypmx,nzpmx,idds)
      use mpplib3, only: SUB => PPNAGUARD32L
      implicit none
      integer, intent(in) :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx
      integer, intent(in) :: idds
      real, dimension(nxv,nypmx,nzpmx), intent(inout) :: f
      real, dimension(nxv,nzpmx,2), intent(inout) :: scs
      real, dimension(nxv,nypmx), intent(inout) :: scr
      integer, dimension(idds), intent(in) :: nyzp
      call SUB(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,idds)
      end subroutine
!
!-----------------------------------------------------------------------
       subroutine PPNACGUARD32L(f,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,  &
     &nxv,nypmx,nzpmx,idds)
      use mpplib3, only: SUB => PPNACGUARD32L
      implicit none
      integer, intent(in) :: ndim, kstrt, nvpy, nvpz, nx, nxv
      integer, intent(in) :: nypmx, nzpmx, idds
      real, dimension(ndim,nxv,nypmx,nzpmx), intent(inout) :: f
      real, dimension(ndim,nxv,nzpmx,2), intent(inout) :: scs
      real, dimension(ndim,nxv,nypmx), intent(inout) :: scr
      integer, dimension(idds), intent(in) :: nyzp
      call SUB(f,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,  &
     &idds)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPTPOS3A(f,g,s,t,nx,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,nxv, &
     &nyv,kxypd,kypd,kzpd)
      use mpplib3, only: SUB => PPTPOS3A
      implicit none
      integer, intent(in) :: nx, ny, nz, kxyp, kyp, kzp, kstrt, nvpy
      integer, intent(in) :: nxv, nyv, kxypd, kypd, kzpd
      complex, dimension(nxv,kypd,kzpd), intent(in) :: f
      complex, dimension(nyv,kxypd,kzpd), intent(inout) :: g
      complex, dimension(kxyp*kyp*kzp), intent(inout) :: s, t
      call SUB(f,g,s,t,nx,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,nxv,nyv,kxypd,  &
     &kypd,kzpd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPTPOS3B(g,h,s,t,nx,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,nvpz&
     &,nyv,nzv,kxypd,kyzpd,kzpd)
      use mpplib3, only: SUB => PPTPOS3B
      implicit none
      integer, intent(in) :: nx, ny, nz, kxyp, kyzp, kzp, kstrt
      integer, intent(in) :: nvpy, nvpz, nyv, nzv, kxypd, kyzpd, kzpd
      complex, dimension(nyv,kxypd,kzpd), intent(inout) :: g
      complex, dimension(nzv,kxypd,kyzpd), intent(inout) :: h
      complex, dimension(kyzp*kxyp*kzp), intent(inout) :: s, t
      call SUB(g,h,s,t,nx,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,nvpz,nyv,nzv,  &
     &kxypd,kyzpd,kzpd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNTPOS3A(f,g,s,t,nx,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,ndim&
     &,nxv,nyv,kxypd,kypd,kzpd)
      use mpplib3, only: SUB => PPNTPOS3A
      implicit none
      integer, intent(in) :: nx, ny, nz, kxyp, kyp, kzp, kstrt, nvpy
      integer, intent(in) :: ndim, nxv, nyv, kxypd, kypd, kzpd
      complex, dimension(ndim,nxv,kypd,kzpd), intent(in) :: f
      complex, dimension(ndim,nyv,kxypd,kzpd), intent(inout) :: g
      complex, dimension(ndim,kxyp*kyp*kzp), intent(inout) :: s, t
      call SUB(f,g,s,t,nx,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,ndim,nxv,nyv,   &
     &kxypd,kypd,kzpd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNTPOS3B(g,h,s,t,nx,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,   &
     &nvpz,ndim,nyv,nzv,kxypd,kyzpd,kzpd)
      use mpplib3, only: SUB => PPNTPOS3B
      implicit none
      integer, intent(in) :: nx, ny, nz, kxyp, kyzp, kzp, kstrt, nvpy
      integer, intent(in) :: nvpz, ndim, nyv, nzv, kxypd, kyzpd, kzpd
      complex, dimension(ndim,nyv,kxypd,kzpd), intent(inout) :: g
      complex, dimension(ndim,nzv,kxypd,kyzpd), intent(inout) :: h
      complex, dimension(ndim,kyzp*kxyp*kzp), intent(inout) :: s, t
      call SUB(g,h,s,t,nx,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,nvpz,ndim,nyv, &
     &nzv,kxypd,kyzpd,kzpd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPMOVE32(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr, &
     &mcls,kstrt,nvpy,nvpz,idimp,nbmax,mx1,myp1,mzp1,mxzyp1,irc)
      use mpplib3, only: SUB => PPPMOVE32
      implicit none
      integer, intent(in) :: kstrt, nvpy, nvpz, idimp, nbmax
      integer, intent(in) :: mx1, myp1, mzp1, mxzyp1
      integer, intent(inout) :: irc
      real, dimension(idimp,nbmax,2), intent(in) :: sbufr, sbufl
      real, dimension(idimp,nbmax,2), intent(inout) :: rbufr, rbufl
      integer, dimension(3,mxzyp1,3,2), intent(inout) :: ncll, nclr
      integer, dimension(3,mxzyp1,3,2), intent(inout) :: mcll, mclr
      integer, dimension(3,mx1+1,4), intent(inout) :: mcls
      call SUB(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,mcls,kstrt,  &
     &nvpy,nvpz,idimp,nbmax,mx1,myp1,mzp1,mxzyp1,irc)
      end subroutine
