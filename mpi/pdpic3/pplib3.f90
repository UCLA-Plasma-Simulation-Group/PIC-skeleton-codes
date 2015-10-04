!-----------------------------------------------------------------------
! Basic parallel PIC library for MPI communications
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
! PPMOVE32 moves particles into appropriate spatial regions with
!          periodic boundary conditions and 2D spatial decomposition.
!          Assumes ihole list has been found.
! written by viktor k. decyk, ucla
! copyright 1995, regents of the university of california
! update: august 29, 2015
      module pplib3
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
      public :: PPTPOS3A, PPTPOS3B, PPNTPOS3A, PPNTPOS3B, PPMOVE32
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
         do k = 1, nyzp(2)
            do j = 1, nxv
               f(j,nyp1,k) = f(j,1,k)
            enddo
         enddo
      else
! buffer data in y
         do k = 1, nyzp(2)
            do j = 1, nxv
               scs(j,k,1) = f(j,1,k)
            enddo
         enddo
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
         do k = 1, nyzp(2)
            do j = 1, nxv
               f(j,nyp1,k) = scs(j,k,2)
           enddo
         enddo
      endif
! special case for one processor in z
      if (nvpz==1) then
         do k = 1, nyp1
            do j = 1, nxv
               f(j,k,nzp1) = f(j,k,1)
            enddo
         enddo
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
         do k = 1, nzp1
            do j = 1, nx1
               f(j,1,k) = f(j,1,k) + f(j,nyp1,k)
               f(j,nyp1,k) = 0.0
            enddo
         enddo
      else
! buffer data in y
         do k = 1, nzp1 
            do j = 1, nxv
               scs(j,k,1) = f(j,nyp1,k)
            enddo
         enddo
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
         do k = 1, nzp1
            do j = 1, nx1
               f(j,1,k) = f(j,1,k) + scs(j,k,2)
               f(j,nyp1,k) = 0.0
            enddo
         enddo
      endif
! special case for one processor in z
      if (nvpz==1) then
         do k = 1, nyp1
            do j = 1, nx1
               f(j,k,1) = f(j,k,1) + f(j,k,nzp1)
               f(j,k,nzp1) = 0.0
            enddo
         enddo
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
      do k = 1, nyp1
         do j = 1, nx1
            f(j,k,1) = f(j,k,1) + scr(j,k)
            f(j,k,nzp1) = 0.0
         enddo
      enddo
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
         do k = 1, nzp1
            do j = 1, nx1
               do i = 1, ndim
                  f(i,j,1,k) = f(i,j,1,k) + f(i,j,nyp1,k)
                  f(i,j,nyp1,k) = 0.0
               enddo
            enddo
         enddo
      else
! buffer data in y
         do k = 1, nzp1 
            do j = 1, nxv
               do i = 1, ndim
                  scs(i,j,k,1) = f(i,j,nyp1,k)
               enddo
            enddo
         enddo
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
         do k = 1, nzp1
           do j = 1, nx1
              do i = 1, ndim
                  f(i,j,1,k) = f(i,j,1,k) + scs(i,j,k,2)
                  f(i,j,nyp1,k) = 0.0
               enddo
            enddo
         enddo
      endif
! special case for one processor in z
      if (nvpz==1) then
         do k = 1, nyp1
            do j = 1, nx1
               do i = 1, ndim
                  f(i,j,k,1) = f(i,j,k,1) + f(i,j,k,nzp1)
                  f(i,j,k,nzp1) = 0.0
               enddo
            enddo
         enddo
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
      do k = 1, nyp1
         do j = 1, nx1
            do i = 1, ndim
               f(i,j,k,1) = f(i,j,k,1) + scr(i,j,k)
               f(i,j,k,nzp1) = 0.0
            enddo
         enddo
      enddo
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
      integer :: ierr, msid
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
         do l = 1, kzp
            do k = 1, kyp
               do j = 1, kxyp
                  g(k,j,l) = f(j,k,l)
               enddo
            enddo
         enddo
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
         do l = 1, kzps
            koff = kyps*(l - 1) - 1
            do k = 1, kyps
               do j = 1, ld
                  s(j+ld*(k+koff)) = f(j+joff,k,l)
               enddo
            enddo
         enddo
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
         do l = 1, kzps
            joff = ld*(l - 1) - 1
            do k = 1, ld
               do j = 1, kxyps
                  g(k+koff,j,l) = t(j+kxyps*(k+joff))
               enddo
            enddo
         enddo
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
      integer :: ld, jd, kxyzp
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
         do l = 1, kzp
            do j = 1, kxyp
               do k = 1, kyzp
                  h(l,j,k) = g(k,j,l)
               enddo
            enddo
         enddo
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
         do l = 1, kzps
            loff = kxyps*(l - 1) - 1
            do j = 1, kxyps
               do k = 1, ld
                  s(k+ld*(j+loff)) = g(k+koff,j,l)
               enddo
            enddo
         enddo
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
         do l = 1, ld
            koff = kxyps*(l - 1) - 1
            do j = 1, kxyps
               do k = 1, kyzps
                  h(l+loff,j,k) = t(k+kyzps*(j+koff))
               enddo
            enddo
         enddo
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
      integer :: joff, koff, ld, jd, kxyzp
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
         do l = 1, kzp
            do k = 1, kyp
               do j = 1, kxyp
                  do i = 1, ndim
                     g(i,k,j,l) = f(i,j,k,l)
                  enddo
               enddo
            enddo
         enddo
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
         do l = 1, kzps
            koff = kyps*(l - 1) - 1
            do k = 1, kyps
               do j = 1, ld
                  do i = 1, ndim
                     s(i,j+ld*(k+koff)) = f(i,j+joff,k,l)
                  enddo
               enddo
            enddo
         enddo
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
         do l = 1, kzps
            joff = ld*(l - 1) - 1
            do k = 1, ld
               do j = 1, kxyps
                  do i = 1, ndim
                     g(i,k+koff,j,l) = t(i,j+kxyps*(k+joff))
                  enddo
               enddo
            enddo
         enddo
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
      integer :: koff, loff, ld, jd, kxyzp
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
         do l = 1, kzp
            do j = 1, kxyp
               do k = 1, kyzp
                  do i = 1, ndim
                     h(i,l,j,k) = g(i,k,j,l)
                  enddo
               enddo
            enddo
         enddo
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
         do l = 1, kzps
            loff = kxyps*(l - 1) - 1
            do j = 1, kxyps
               do k = 1, ld
                  do i = 1, ndim
                     s(i,k+ld*(j+loff)) = g(i,k+koff,j,l)
                  enddo
               enddo
            enddo
         enddo
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
         do l = 1, ld
            koff = kxyps*(l - 1) - 1
            do j = 1, kxyps
               do k = 1, kyzps
                  do i = 1, ndim
                     h(i,l+loff,j,k) = t(i,k+kyzps*(j+koff))
                  enddo
               enddo
            enddo
         enddo
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole, &
     &ny,nz,kstrt,nvpy,nvpz,idimp,npmax,idps,nbmax,ntmax,info)
! this subroutine moves particles into appropriate spatial regions
! ihole array is calculated in particle push procedure
! with periodic boundary conditions and 2D spatial decomposition
! output: part, ihole, npp, sbufr, sbufl, rbufr, rbufl, info
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = position z of particle n in partition
! part(4,n) = velocity vx of particle n in partition
! part(5,n) = velocity vy of particle n in partition
! part(6,n) = velocity vz of particle n in partition m
! edges(1:2) = lower/upper boundary in y of particle partition
! edges(3:4) = back/front boundary in z of particle partition
! npp = number of particles in partition
! sbufl = buffer for particles being sent to back processor
! sbufr = buffer for particles being sent to front processor
! rbufl = buffer for particles being received from back processor
! rbufr = buffer for particles being received from front processor
! ihole = location of holes left in particle arrays
! ny/nz = system length in y/z direction
! kstrt = starting data block number
! nvpy/nvpz = number of real or virtual processors in y/z
! idimp = size of phase space = 6
! npmax = maximum number of particles in each partition.
! idps = number of particle partition boundaries = 4
! nbmax =  size of buffers for passing particles between processors
! ntmax =  size of hole array for particles leaving processors
! info = status information
! info(1) = ierr = (0,N) = (no,yes) error condition exists
! info(2) = maximum number of particles per processor
! info(3) = minimum number of particles per processor
! info(4:5) = maximum number of buffer overflows in y/z
! info(6:7) = maximum number of particle passes required in y/z
      implicit none
      integer, intent(in) :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
      integer, intent(in) :: idps, nbmax, ntmax
      integer, intent(inout) :: npp
      real, dimension(idimp,npmax), intent(inout) :: part
      real, dimension(idps), intent(in) :: edges
      real, dimension(idimp,nbmax), intent(inout) :: sbufr, sbufl
      real, dimension(idimp,nbmax), intent(inout) :: rbufr, rbufl
      integer, dimension(ntmax+1,2) , intent(inout):: ihole
      integer, dimension(7), intent(inout) :: info
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
! iy/iz = partitioned co-ordinates
      integer, parameter :: iy = 2, iz = 3
      integer :: i, j, n, js, ks, ic, nvp, iter, nps, kl, kr, j1, j2
      integer :: ih, jh, nh, j3, joff, jin, nbsize, nter, mter, itermax
      integer :: itg, ierr
      real :: an, xt
      integer, dimension(4) :: msid
      integer, dimension(lstat) :: istatus
      integer, dimension(2) :: kb, jsl, jsr, jss
      integer, dimension(5) :: ibflg, iwork
! js/ks = processor co-ordinates in y/z=> idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      nbsize = idimp*nbmax
      info(1) = 0
      info(6) = 0
      info(7) = 0
      itermax = 2000
! buffer outgoing particles, first in y then in z direction
         do n = 1, 2
         if (n==1) then
            ic = iy
            nvp = nvpy
            an = real(ny)
            jh = ihole(1,n+1)
            ibflg(5) = 0
      else if (n==2) then
            ic = iz
            nvp = nvpz
            an = real(nz)
         endif
         iter = 2
         nter = 0
         ih = ihole(1,n)
         joff = 1
         jin = 1
! ih = number of particles extracted from holes
! joff = next hole location for extraction
! jss(1) = number of holes available to be filled
! jin = next hole location to be filled
! start loop
   10    mter = 0
         nps = 0
         kb(1) = js
         kb(2) = ks
! buffer outgoing particles
         jsl(1) = 0
         jsr(1) = 0
! load particle buffers
         do j = 1, ih
            j1 = ihole(j+joff,n)
            xt = part(ic,j1)
! particles going down or backward
            if (xt < edges(2*n-1)) then
               if (kb(n)==0) xt = xt + an
               if (jsl(1) < nbmax) then
                  jsl(1) = jsl(1) + 1
                  do i = 1, idimp
                     sbufl(i,jsl(1)) = part(i,j1)
                  enddo
                  sbufl(ic,jsl(1)) = xt
               else
                  nps = 1
                  exit
               endif
! particles going up or forward
            else
               if (kb(n)==(nvp-1)) xt = xt - an
               if (jsr(1) < nbmax) then
                  jsr(1) = jsr(1) + 1
                  do i = 1, idimp
                     sbufr(i,jsr(1)) = part(i,j1)
                  enddo
                  sbufr(ic,jsr(1)) = xt
               else
                  nps = 1
                  exit
               endif
            endif
         enddo
         jss(1) = jsl(1) + jsr(1)
         joff = joff + jss(1)
         ih = ih - jss(1)
! check for full buffer condition
         ibflg(3) = nps
! copy particle buffers
   60    iter = iter + 2
         mter = mter + 1
! special case for one processor
         if (nvp==1) then
            jsl(2) = jsr(1)
            do j = 1, jsl(2)
               do i = 1, idimp
                  rbufl(i,j) = sbufr(i,j)
               enddo
            enddo
            jsr(2) = jsl(1)
            do j = 1, jsr(2)
               do i = 1, idimp
                  rbufr(i,j) = sbufl(i,j)
               enddo
            enddo
! this segment is used for mpi computers
         else
! get particles from below and above or back and front
            kb(1) = js
            kb(2) = ks
            kl = kb(n)
            kb(n) = kl + 1
            if (kb(n) >= nvp) kb(n) = kb(n) - nvp
            kr = kb(1) + nvpy*kb(2)
            kb(n) = kl - 1
            if (kb(n) < 0) kb(n) = kb(n) + nvp
            kl = kb(1) + nvpy*kb(2)
! post receive
            itg = iter - 1
            call MPI_IRECV(rbufl,nbsize,mreal,kl,itg,lgrp,msid(1),ierr)
            call MPI_IRECV(rbufr,nbsize,mreal,kr,iter,lgrp,msid(2),ierr)
! send particles
            jsr(1) = idimp*jsr(1)
            call MPI_ISEND(sbufr,jsr(1),mreal,kr,itg,lgrp,msid(3),ierr)
            jsl(1) = idimp*jsl(1)
            call MPI_ISEND(sbufl,jsl(1),mreal,kl,iter,lgrp,msid(4),ierr)
! wait for particles to arrive
            call MPI_WAIT(msid(1),istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,nps,ierr)
            jsl(2) = nps/idimp
            call MPI_WAIT(msid(2),istatus,ierr)
            call MPI_GET_COUNT(istatus,mreal,nps,ierr)
            jsr(2) = nps/idimp
         endif
! check if particles must be passed further
! check if any particles coming from above or front belong here
         jsl(1) = 0
         jsr(1) = 0
         jss(2) = 0
         do j = 1, jsr(2)
            if (rbufr(ic,j) < edges(2*n-1)) jsl(1) = jsl(1) + 1
            if (rbufr(ic,j) >= edges(2*n)) jsr(1) = jsr(1) + 1
         enddo
!        if (jsr(1) /= 0) then
!           if (n==1) then
!              write (2,*) kb+1,'Info: particles returning above'
!           else if (n==2) then
!              write (2,*) kb+1,'Info: particles returning front'
!           endif
!       endif
! check if any particles coming from below or back belong here
         do j = 1, jsl(2)
            if (rbufl(ic,j) >= edges(2*n)) jsr(1) = jsr(1) + 1
            if (rbufl(ic,j) < edges(2*n-1)) jss(2) = jss(2) + 1
         enddo
!        if (jss(2) /= 0) then
!           if (n==1) then
!              write (2,*) kb+1,'Info: particles returning below'
!           else if (n==2) then
!              write (2,*) kb+1,'Info: particles returning back'
!           endif
!        endif
         nps = jsl(1) + jsr(1) + jss(2)
         ibflg(2) = nps
! make sure sbufr and sbufl have been sent
         if (nvp /= 1) then
            call MPI_WAIT(msid(3),istatus,ierr)
            call MPI_WAIT(msid(4),istatus,ierr)
         endif
         if (nps==0) go to 180
! remove particles which do not belong here
         kb(1) = js
         kb(2) = ks
! first check particles coming from above or front
         jsl(1) = 0
         jsr(1) = 0
         jss(2) = 0
         do j = 1, jsr(2)
            xt = rbufr(ic,j)
! particles going down or back
            if (xt < edges(2*n-1)) then
               jsl(1) = jsl(1) + 1
               if (kb(n)==0) xt = xt + an
               rbufr(ic,j) = xt
               do i = 1, idimp
                  sbufl(i,jsl(1)) = rbufr(i,j)
               enddo
! particles going up or front, should not happen
            else if (xt >= edges(2*n)) then
               jsr(1) = jsr(1) + 1
               if (kb(n)==(nvp-1)) xt = xt - an
               rbufr(ic,j) = xt
               do i = 1, idimp
                  sbufr(i,jsr(1)) = rbufr(i,j)
               enddo
! particles staying here
            else
               jss(2) = jss(2) + 1
               do i = 1, idimp
                  rbufr(i,jss(2)) = rbufr(i,j)
               enddo
            endif
         enddo
         jsr(2) = jss(2)
! next check particles coming from below or back
         jss(2) = 0
         do j = 1, jsl(2)
            xt = rbufl(ic,j)
! particles going up or front
            if (xt >= edges(2*n)) then
               if (jsr(1) < nbmax) then
                  jsr(1) = jsr(1) + 1
                  if (kb(n)==(nvp-1)) xt = xt - an
                  rbufl(ic,j) = xt
                  do i = 1, idimp
                     sbufr(i,jsr(1)) = rbufl(i,j)
                  enddo
               else
                  jss(2) = 2*npmax
                  exit
               endif
! particles going down or back, should not happen
            else if (xt < edges(2*n-1)) then
               if (jsl(1) < nbmax) then
                  jsl(1) = jsl(1) + 1
                  if (kb(n)==0) xt = xt + an
                  rbufl(ic,j) = xt
                  do i = 1, idimp
                     sbufl(i,jsl(1)) = rbufl(i,j)
                  enddo
               else
                  jss(2) = 2*npmax
                  exit
               endif
! particles staying here
            else
               jss(2) = jss(2) + 1
               do i = 1, idimp
                  rbufl(i,jss(2)) = rbufl(i,j)
               enddo
            endif
         enddo
         jsl(2) = jss(2)
! check if move would overflow particle array
  180    nps = npp + jsl(2) + jsr(2) - jss(1)
         ibflg(1) = nps
         ibflg(4) = -min0(npmax,nps)
         call PPIMAX(ibflg,iwork,5)
         info(2) = ibflg(1)
         info(3) = -ibflg(4)
         ierr = ibflg(1) - npmax
         if (ierr > 0) then
            write (2,*) 'particle overflow error, ierr = ', ierr
            info(1) = ierr
            return
         endif
! check for ihole overflow condition
         ierr = ibflg(5) - ntmax
         if (ierr > 0) then
            write (2,*) 'ihole overflow error, ierr = ', ierr
            info(1) = -ierr
            return
         endif
! distribute incoming particles from buffers
         nh = 0
! distribute particles coming from below or back into holes
         jss(2) = min0(jss(1),jsl(2))
         do j = 1, jss(2)
            j1 = ihole(j+jin,n)
            do i = 1, idimp
               part(i,j1) = rbufl(i,j)
            enddo
! check if incoming particle is also out of bounds in z
            if (n==1) then
               xt = part(iz,j1)
! if so, add it to list of particles in z
               if ((xt < edges(2*n+1)).or.(xt >= edges(2*n+2))) then
                  jh = jh + 1
                  if (jh <= ntmax) then
                     ihole(jh+1,n+1) = j1
                  else
                     nh = 1
                 endif
               endif
            endif
         enddo
         jin = jin + jss(2)
         if (jss(1) > jsl(2)) then
            jss(2) = min0(jss(1)-jsl(2),jsr(2))
         else
            jss(2) = jsl(2) - jss(1)
         endif
         do j = 1, jss(2)
! no more particles coming from below or back
! distribute particles coming from above or front into holes
            if (jss(1) > jsl(2)) then
               j1 = ihole(j+jin,n)
               do i = 1, idimp
                  part(i,j1) = rbufr(i,j)
               enddo
! check if incoming particle is also out of bounds in z
               if (n==1) then
                  xt = part(iz,j1)
! if so, add it to list of particles in z
                  if ((xt < edges(2*n+1)).or.(xt >= edges(2*n+2))) then
                     jh = jh + 1
                     if (jh <= ntmax) then
                        ihole(jh+1,n+1) = j1
                     else
                        nh = 1
                     endif
                  endif
             endif
! no more holes
! distribute remaining particles from below or back into bottom
            else
               do i = 1, idimp
                  part(i,j+npp) = rbufl(i,j+jss(1))
               enddo
! check if incoming particle is also out of bounds in z
               if (n==1) then
                  xt = part(iz,j+npp)
! if so, add it to list of particles in z
                  if ((xt < edges(2*n+1)).or.(xt >= edges(2*n+2))) then
                     jh = jh + 1
                     if (jh <= ntmax) then
                        ihole(jh+1,n+1) = j + npp
                     else
                        nh = 1
                     endif
                  endif
               endif
            endif
         enddo
         if (jss(1) > jsl(2)) jin = jin + jss(2)
         nps = jsl(2) + jsr(2)
         if (jss(1) <= jsl(2)) then
            npp = npp + (jsl(2) - jss(1))
            jss(1) = jsl(2)
         endif
! no more holes
! distribute remaining particles from above into bottom
         jsr(2) = max0(0,nps-jss(1))
         jss(1) = jss(1) - jsl(2)
         do j = 1, jsr(2)
            do i = 1, idimp
               part(i,j+npp) = rbufr(i,j+jss(1))
            enddo
! check if incoming particle is also out of bounds in z
            if (n==1) then
               xt = part(iz,j+npp)
! if so, add it to list of particles in z
               if ((xt < edges(2*n+1)).or.(xt >= edges(2*n+2))) then
                  jh = jh + 1
                  if (jh <= ntmax) then
                     ihole(jh+1,n+1) = j + npp
                  else
                     nh = 1
                  endif
               endif
            endif

         enddo
         npp = npp + jsr(2)
! check for ihole overflow condition
         if ((n==1).and.(nh > 0)) ibflg(5) = jh
! holes left over
! fill up remaining holes in particle array with particles from bottom
         if (ih==0) then
            jsr(2) = max0(0,ihole(1,n)-jin+1)
            nh = 0
! holes are stored in increasing value
            if (n==1) then
               do j = 1, jsr(2)
                  j1 = npp - j + 1
                  j2 = ihole(jsr(2)-j+jin+1,n)
                  if (j1 > j2) then
! move particle only if it is below current hole
                     do i = 1, idimp
                        part(i,j2) = part(i,j1)
                     enddo
! check if this move makes the ihole list for z invalid
                     xt = part(iz,j1)
                     if ((xt<edges(2*n+1)).or.(xt>=edges(2*n+2))) then
                        i = jh + 1
! if so, adjust the list of holes
                        j3 = ihole(i,n+1)
                        do while (j3 /= j1)
                           i = i - 1
                           if (i==1) then
                              write (2,*) kstrt,                        &
     &                        'cannot find particle:n,j1=', n, j1
                              nh = 1
                              exit
                           endif
                           j3 = ihole(i,n+1)
                        enddo
! update ihole list to use new location
                        ihole(i,n+1) = j2
                     endif
                  endif
               enddo
! holes may not be stored in increasing value
            else
               do j = 1, jsr(2)
                  j1 = npp - j + 1
                  j2 = ihole(jsr(2)-j+jin+1,n)
                  xt = part(iz,j1)
! determine if particle at location j1 represents an unfilled hole
                  if ((xt < edges(2*n-1)).or.(xt >= edges(2*n))) then
                     i = jh + 2 - j
! if so, adjust the list of holes
                     j3 = ihole(i,n)
                     do while (j3 /= j1)
                        i = i - 1
                        if (i==1) then
                           write (2,*) kstrt,                           &
     &                     'cannot find particle:n,j1=', n, j1
                           nh = 1
                           exit
                        endif
                        j3 = ihole(i,n)
                     enddo
! update ihole list to use new location
                     ihole(i,n) = j2
                  else if (j1 > j2) then
! move particle only if it is below current hole
                     do i = 1, idimp
                        part(i,j2) = part(i,j1)
                     enddo
                  endif
               enddo
! check for lost particle error
               if (nh > 0) call PPABORT
            endif
            jin = jin + jsr(2)
            npp = npp - jsr(2)
         endif
         jss(1) = 0
! check if any particles have to be passed further
         if (ibflg(3) > 0) ibflg(3) = 1
         info(5+n) = max0(info(5+n),mter)
         if (ibflg(2) > 0) then
!           write (2,*) 'Info: particles being passed further = ',      &
!    &                   ibflg(2)
            if (iter < itermax) go to 60
            ierr = -((iter-2)/2)
            if (kstrt==1) write (2,*) 'Iteration overflow, iter = ',    &
     &ierr
            info(1) = ierr
            return
         endif
! check if buffer overflowed and more particles remain to be checked
         if (ibflg(3) > 0) then
            nter = nter + 1
            info(3+n) = nter
            go to 10
         endif
         if (nter > 0) then
            if (kstrt==1) then
               write (2,*) 'Info: ',nter,' buffer overflows, nbmax=',   &
     &                     nbmax
            endif
         endif
! update ihole number in z
         if (n==1) ihole(1,2) = jh
      enddo
      end subroutine
!
      end module
!
! Make functions callable by Fortran77
!
!-----------------------------------------------------------------------
      subroutine PPINIT2(idproc,nvp)
      use pplib3, only: SUB => PPINIT2
      implicit none
      integer, intent(inout) :: idproc, nvp
      call SUB(idproc,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPEXIT
      use pplib3, only: SUB => PPEXIT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPABORT
      use pplib3, only: SUB => PPABORT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
      use pplib3, only: SUB => PWTIMERA
      implicit none
      integer, intent(in) :: icntrl
      real, intent(inout) :: time
      double precision, intent(inout) :: dtime
      call SUB(icntrl,time,dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPSUM(f,g,nxp)
      use pplib3, only: SUB => PPSUM
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSUM(f,g,nxp)
      use pplib3, only: SUB => PPDSUM
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPIMAX(if,ig,nxp)
      use pplib3, only: SUB => PPIMAX
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if, ig
      call SUB(if,ig,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDMAX(f,g,nxp)
      use pplib3, only: SUB => PPDMAX
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx&
     &,idds)
      use pplib3, only: SUB => PPNCGUARD32L
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
      use pplib3, only: SUB => PPNAGUARD32L
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
      use pplib3, only: SUB => PPNACGUARD32L
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
      use pplib3, only: SUB => PPTPOS3A
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
      use pplib3, only: SUB => PPTPOS3B
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
      use pplib3, only: SUB => PPNTPOS3A
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
      use pplib3, only: SUB => PPNTPOS3B
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
      subroutine PPMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole, &
     &ny,nz,kstrt,nvpy,nvpz,idimp,npmax,idps,nbmax,ntmax,info)
      use pplib3, only: SUB => PPMOVE32
      implicit none
      integer, intent(in) :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
      integer, intent(in) :: idps, nbmax, ntmax
      integer, intent(inout) :: npp
      real, dimension(idimp,npmax), intent(inout) :: part
      real, dimension(idps), intent(in) :: edges
      real, dimension(idimp,nbmax), intent(inout) :: sbufr, sbufl
      real, dimension(idimp,nbmax), intent(inout) :: rbufr, rbufl
      integer, dimension(ntmax+1,2) , intent(inout):: ihole
      integer, dimension(7), intent(inout) :: info
      call SUB(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,ny,nz,kstrt,&
     &nvpy,nvpz,idimp,npmax,idps,nbmax,ntmax,info)
      end subroutine

