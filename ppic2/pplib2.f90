!-----------------------------------------------------------------------
! Basic parallel PIC library for MPI communications
! pplib2.f90 contains basic communications procedures for 1d partitions:
! PPINIT2 initializes parallel processing for Fortran90, returns
!         number of processors and processor id.
! PPEXIT terminates parallel processing.
! PPABORT aborts parallel processing.
! PWTIMERA performs parallel local wall clock timing.
! PPSUM performs parallel sum of a real vector.
! PPDSUM performs parallel sum of a double precision vector.
! PPIMAX performs parallel maximum of an integer vector.
! PPDMAX performs parallel maximum of a double precision vector.
! PPNCGUARD2L copies data to guard cells in y for scalar data, linear
!             interpolation, and distributed data with non-uniform
!             partition.
! PPNAGUARD2L adds guard cells in y for scalar array, linear
!             interpolation, and distributed data with non-uniform
!             partition.
! PPNACGUARD2L adds guard cells in y for vector array, linear
!              interpolation, and distributed data with non-uniform
!              partition.
! PPTPOSE performs a transpose of a complex scalar array, distributed
!         in y, to a complex scalar array, distributed in x.
! PPNTPOSE performs a transpose of an n component complex vector array,
!          distributed in y, to an n component complex vector array,
!          distributed in x.
! PPMOVE2 moves particles into appropriate spatial regions with periodic
!         boundary conditions.  Assumes ihole list has been found.
! written by viktor k. decyk, ucla
! copyright 1995, regents of the university of california
! update: april 23, 2015
      module pplib2
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
      public :: PPNCGUARD2L, PPNAGUARD2L, PPNACGUARD2L
      public :: PPTPOSE, PPNTPOSE, PPMOVE2
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
      subroutine PPNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
! this subroutine copies data to guard cells in non-uniform partitions
! f(j,k) = real data for grid j,k in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! output: f
! nyp = number of primary gridpoints in field partition
! it is assumed the nyp > 0.
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of field partition, including guard cell.
! linear interpolation, for distributed data
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: j, ks, moff, kl, kr
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
! special case for one processor
      if (nvp==1) then
         do j = 1, nxv
            f(j,nyp+1) = f(j,1)
         enddo
         return
      endif
      ks = kstrt - 1
      moff = nypmx*nvp + 2
! copy to guard cells
      kr = ks + 1
      if (kr >= nvp) kr = kr - nvp
      kl = ks - 1
      if (kl < 0)  kl = kl + nvp
      ks = nyp + 1
! this segment is used for mpi computers
      call MPI_IRECV(f(1,ks),nxv,mreal,kr,moff,lgrp,msid,ierr)
      call MPI_SEND(f,nxv,mreal,kl,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
! this subroutine adds data from guard cells in non-uniform partitions
! f(j,k) = real data for grid j,k in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! output: f, scr
! scr(j) = scratch array for particle partition
! nyp = number of primary gridpoints in particle partition
! it is assumed the nyp > 0.
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nx = system length in x direction
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of field partition, including guard cells.
! linear interpolation, for distributed data
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      real, dimension(nxv), intent(inout) :: scr
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: j, nx1, ks, moff, kl, kr
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      nx1 = nx + 1
! special case for one processor
      if (nvp==1) then
         do j = 1, nx1
            f(j,1) = f(j,1) + f(j,nyp+1)
            f(j,nyp+1) = 0.
         enddo
         return
      endif
      ks = kstrt - 1
      moff = nypmx*nvp + 1
! add guard cells
      kr = ks + 1
      if (kr >= nvp) kr = kr - nvp
      kl = ks - 1
      if (kl < 0) kl = kl + nvp
      ks = nyp + 1
! this segment is used for mpi computers
      call MPI_IRECV(scr,nxv,mreal,kl,moff,lgrp,msid,ierr)
      call MPI_SEND(f(1,ks),nxv,mreal,kr,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
! add up the guard cells
      do j = 1, nx1
         f(j,1) = f(j,1) + scr(j)
         f(j,nyp+1) = 0.0
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
! this subroutine adds data from guard cells in non-uniform partitions
! f(ndim,j,k) = real data for grid j,k in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! output: f, scr
! scr(ndim,j) = scratch array for particle partition
! nyp = number of primary gridpoints in particle partition
! it is assumed the nyp > 0.
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nx = system length in x direction
! ndim = leading dimension of array f
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of field partition, including guard cells.
! linear interpolation, for distributed data
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx, ndim, nxv, nypmx
      real, dimension(ndim,nxv,nypmx), intent(inout) :: f
      real, dimension(ndim,nxv), intent(inout) :: scr
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: j, n, nx1, ks, moff, kl, kr
      integer :: nnxv
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      nx1 = nx + 1
! special case for one processor
      if (nvp==1) then
         do j = 1, nx1
            do n = 1, ndim
               f(n,j,1) = f(n,j,1) + f(n,j,nyp+1)
               f(n,j,nyp+1) = 0.0
            enddo
         enddo
         return
      endif
      ks = kstrt - 1
      moff = nypmx*nvp + 1
      nnxv = ndim*nxv
! add guard cells
      kr = ks + 1
      if (kr >= nvp) kr = kr - nvp
      kl = ks - 1
      if (kl < 0) kl = kl + nvp
      ks = nyp + 1
! this segment is used for mpi computers
      call MPI_IRECV(scr,nnxv,mreal,kl,moff,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,ks),nnxv,mreal,kr,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
! add up the guard cells
      do j = 1, nx1
         do n = 1, ndim
            f(n,j,1) = f(n,j,1) + scr(n,j)
            f(n,j,nyp+1) = 0.0
         enddo
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd,  &
     &kypd)
! this subroutine performs a transpose of a matrix f, distributed in y,
! to a matrix g, distributed in x, that is,
! g(k+kyp*(m-1),j,l) = f(j+kxp*(l-1),k,m), where
! 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
! and where indices l and m can be distributed across processors.
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! f = complex input array
! g = complex output array
! s, t = complex scratch arrays
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nxv/nyv = first dimension of f/g
! kypd/kxpd = second dimension of f/g
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv
      integer, intent(in) :: kxpd, kypd
      complex, dimension(nxv,kypd), intent(in) :: f
      complex, dimension(nyv,kxpd), intent(inout) :: g
      complex, dimension(kxp*kyp), intent(inout) :: s, t
! lgrp = current communicator
! mcplx = default datatype for complex
! local data
      integer :: n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld
      integer :: ierr, msid
      integer, dimension(lstat) :: istatus
      ks = kstrt - 1
      kxps = min(kxp,max(0,nx-kxp*ks))
      kyps = min(kyp,max(0,ny-kyp*ks))
      kxyp = kxp*kyp
! special case for one processor
      if (nvp==1) then
         do k = 1, kyp
            do j = 1, kxp
               g(k,j) = f(j,k)
            enddo
         enddo
         return
      endif
! this segment is used for shared memory computers
!     do m = 1, min(ny,nvp)
!        koff = kyp*(m - 1)
!        do k = 1, min(kyp,max(0,ny-koff))
!           do l = 1, min(nx,nvp)
!              joff = kxp*(l - 1)
!              do j = 1, min(kxp,max(0,nx-joff))
!                 g(k+koff,j+joff) = f(j+joff,k+koff)
!              enddo
!           enddo
!        enddo
!     enddo
! this segment is used for mpi computers
      do n = 1, nvp
         id = n - ks - 1
         if (id < 0) id = id + nvp
! extract data to send
         joff = kxp*id
         ld = min(kxp,max(0,nx-joff))
         do k = 1, kyps
            do j = 1, ld
               s(j+ld*(k-1)) = f(j+joff,k)
            enddo
         enddo
         ld = ld*kyps
! post receive
         call MPI_IRECV(t,kxyp,mcplx,id,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mcplx,id,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         koff = kyp*id
         ld = min(kyp,max(0,ny-koff))
         do k = 1, ld
            do j = 1, kxps
               g(k+koff,j) = t(j+kxps*(k-1))
            enddo
         enddo
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv, &
     &kxpd,kypd)
! this subroutine performs a transpose of a matrix f, distributed in y,
! to a matrix g, distributed in x, that is,
! g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxp*(l-1),k,m), where
! 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
! and where indices l and m can be distributed across processors.
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! f = complex input array
! g = complex output array
! s, t = complex scratch arrays
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
! ndim = leading dimension of arrays f and g
! nxv/nyv = first dimension of f/g
! kypd/kxpd = second dimension of f/g
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      integer, intent(in) :: nxv, nyv, kxpd, kypd
      complex, dimension(ndim,nxv,kypd), intent(in) :: f
      complex, dimension(ndim,nyv,kxpd), intent(inout) :: g
      complex, dimension(ndim,kxp*kyp), intent(inout) :: s, t
! lgrp = current communicator
! mcplx = default datatype for complex
! local data
      integer :: i, n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld
      integer :: ierr, msid
      integer, dimension(lstat) :: istatus
      ks = kstrt - 1
      kxps = min(kxp,max(0,nx-kxp*ks))
      kyps = min(kyp,max(0,ny-kyp*ks))
      kxyp = ndim*kxp*kyp
! special case for one processor
      if (nvp==1) then
         do k = 1, kyp
            do j = 1, kxp
               do i = 1, ndim
                  g(i,k,j) = f(i,j,k)
               enddo
            enddo
         enddo
         return
      endif
! this segment is used for shared memory computers
!     do m = 1, min(ny,nvp)
!        koff = kyp*(m - 1)
!        do k = 1, min(kyp,max(0,ny-koff))
!           do l = 1, min(nx,nvp)
!              joff = kxp*(l - 1)
!              do j = 1, min(kxp,max(0,nx-joff))
!                 do i = 1, ndim
!                    g(i,k+koff,j+joff) = f(i,j+joff,k+koff)
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
! this segment is used for mpi computers
      do n = 1, nvp
         id = n - ks - 1
         if (id.lt.0) id = id + nvp
! extract data to send
         joff = kxp*id
         ld = min(kxp,max(0,nx-joff))
         do k = 1, kyps
            do j = 1, ld
               do i = 1, ndim
                  s(i,j+ld*(k-1)) = f(i,j+joff,k)
               enddo
            enddo
         enddo
         ld = ndim*ld*kyps
! post receive
         call MPI_IRECV(t,kxyp,mcplx,id,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mcplx,id,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         koff = kyp*id
         ld = min(kyp,max(0,ny-koff))
         do k = 1, ld
            do j = 1, kxps
               do i = 1, ndim
                  g(i,k+koff,j) = t(i,j+kxps*(k-1))
               enddo
            enddo
         enddo
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,ny&
     &,kstrt,nvp,idimp,npmax,idps,nbmax,ntmax,info)
! this subroutine moves particles into appropriate spatial regions
! periodic boundary conditions
! output: part, npp, sbufr, sbufl, rbufr, rbufl, info
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! part(3,n) = velocity vx of particle n in partition
! part(4,n) = velocity vy of particle n in partition
! edges(1:2) = lower:lower boundary of particle partition
! npp = number of particles in partition
! sbufl = buffer for particles being sent to lower processor
! sbufr = buffer for particles being sent to upper processor
! rbufl = buffer for particles being received from lower processor
! rbufr = buffer for particles being received from upper processor
! ihole = location of holes left in particle arrays
! ny = system length in y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 4
! npmax = maximum number of particles in each partition.
! idps = number of partition boundaries
! nbmax =  size of buffers for passing particles between processors
! ntmax =  size of hole array for particles leaving processors
! info = status information
! info(1) = ierr = (0,N) = (no,yes) error condition exists
! info(2) = maximum number of particles per processor
! info(3) = minimum number of particles per processor
! info(4) = maximum number of buffer overflows
! info(5) = maximum number of particle passes required
      implicit none
      integer, intent(in) :: ny, kstrt, nvp, idimp, npmax, idps
      integer, intent(in) :: nbmax, ntmax
      integer, intent(inout) :: npp
      real, dimension(idimp,npmax), intent(inout) :: part
      real, dimension(idps), intent(in) :: edges
      real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr, rbufl, rbufr
      integer, dimension(ntmax+1), intent(inout) :: ihole
      integer, dimension(5), intent(inout) :: info
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
! iy = partitioned co-ordinate
      integer, parameter :: iy = 2
      integer :: ierr, ks, ih, iter, nps, itg, kl, kr, j, j1, j2, i
      integer :: joff, jin, nbsize, nter, mter, itermax
      real :: any, yt
      integer, dimension(4) :: msid
      integer, dimension(lstat) :: istatus
      integer, dimension(2) :: jsl, jsr, jss
      integer, dimension(4) :: ibflg, iwork
      any = real(ny)
      ks = kstrt - 1
      nbsize = idimp*nbmax
      iter = 2
      nter = 0
      info(1) = 0
      info(5) = 0
      ih = ihole(1)
      joff = 1
      jin = 1
      itermax = 2000
! ih = number of particles extracted from holes
! joff = next hole location for extraction
! jss(1) = number of holes available to be filled
! jin = next hole location to be filled
! start loop
   10 mter = 0
      nps = 0
! buffer outgoing particles
      jsl(1) = 0
      jsr(1) = 0
! load particle buffers
      do j = 1, ih
         j1 = ihole(j+joff)
         yt = part(iy,j1)
! particles going down
         if (yt < edges(1)) then
            if (ks==0) yt = yt + any
            if (jsl(1) < nbmax) then
               jsl(1) = jsl(1) + 1
               do i = 1, idimp
                  sbufl(i,jsl(1)) = part(i,j1)
               enddo
               sbufl(iy,jsl(1)) = yt
            else
               nps = 1
               exit
            endif
! particles going up
         else
            if (ks==(nvp-1)) yt = yt - any
            if (jsr(1) < nbmax) then
               jsr(1) = jsr(1) + 1
               do i = 1, idimp
                  sbufr(i,jsr(1)) = part(i,j1)
               enddo
               sbufr(iy,jsr(1)) = yt
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
   60 iter = iter + 2
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
! get particles from below and above
         kr = ks + 1
         if (kr >= nvp) kr = kr - nvp
         kl = ks - 1
         if (kl < 0) kl = kl + nvp
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
! check if any particles coming from above belong here
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do j = 1, jsr(2)
         if (rbufr(iy,j) < edges(1)) jsl(1) = jsl(1) + 1
         if (rbufr(iy,j) >= edges(2)) jsr(1) = jsr(1) + 1
      enddo
!     if (jsr(1) /= 0) write (2,*) ks+1, 'Info: particles returning up'
! check if any particles coming from below belong here
      do j = 1, jsl(2)
         if (rbufl(iy,j) >= edges(2)) jsr(1) = jsr(1) + 1
         if (rbufl(iy,j) < edges(1)) jss(2) = jss(2) + 1
      enddo
!     if (jss(2) /= 0) write (2,*) ks+1,'Info: particles returning down'
      nps = jsl(1) + jsr(1) + jss(2)
      ibflg(2) = nps
! make sure sbufr and sbufl have been sent
      if (nvp /= 1) then
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
      endif
      if (nps==0) go to 180
! remove particles which do not belong here
! first check particles coming from above
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do j = 1, jsr(2)
         yt = rbufr(iy,j)
! particles going down
         if (yt < edges(1)) then
            jsl(1) = jsl(1) + 1
            if (ks==0) yt = yt + any
            rbufr(iy,j) = yt
            do i = 1, idimp
               sbufl(i,jsl(1)) = rbufr(i,j)
            enddo
! particles going up, should not happen
         else if (yt >= edges(2)) then
            jsr(1) = jsr(1) + 1
            if (ks==(nvp-1)) yt = yt - any
            rbufr(iy,j) = yt
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
! next check particles coming from below
      jss(2) = 0
      do j = 1, jsl(2)
         yt = rbufl(iy,j)
! particles going up
         if (yt >= edges(2)) then
            if (jsr(1) < nbmax) then
               jsr(1) = jsr(1) + 1
               if (ks==(nvp-1)) yt = yt - any
               rbufl(iy,j) = yt
               do i = 1, idimp
                  sbufr(i,jsr(1)) = rbufl(i,j)
               enddo
            else
               jss(2) = 2*npmax
               exit
            endif
! particles going down, should not happen
         else if (yt < edges(1)) then
            if (jsl(1) < nbmax) then
               jsl(1) = jsl(1) + 1
               if (ks==0) yt = yt + any
               rbufl(iy,j) = yt
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
  180 nps = npp + jsl(2) + jsr(2) - jss(1)
      ibflg(1) = nps
      ibflg(4) = -min0(npmax,nps)
      call PPIMAX(ibflg,iwork,4)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr > 0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
! distribute incoming particles from buffers
! distribute particles coming from below into holes
      jss(2) = min0(jss(1),jsl(2))
      do j = 1, jss(2)
         j1 = ihole(j+jin)
         do i = 1, idimp
            part(i,j1) = rbufl(i,j)
         enddo
      enddo
      jin = jin + jss(2)
      if (jss(1) > jsl(2)) then
         jss(2) = min0(jss(1)-jsl(2),jsr(2))
      else
         jss(2) = jsl(2) - jss(1)
      endif
      do j = 1, jss(2)
! no more particles coming from below
! distribute particles coming from above into holes
         if (jss(1) > jsl(2)) then
            j1 = ihole(j+jin)
            do i = 1, idimp
               part(i,j1) = rbufr(i,j)
            enddo
! no more holes
! distribute remaining particles from below into bottom
         else
            do i = 1, idimp
               part(i,j+npp) = rbufl(i,j+jss(1))
            enddo
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
      enddo
      npp = npp + jsr(2)
! holes left over
! fill up remaining holes in particle array with particles from bottom
      if (ih.eq.0) then
         jsr(2) = max0(0,ihole(1)-jin+1)
         do j = 1, jsr(2)
            j1 = npp - j + 1
            j2 = ihole(jsr(2)-j+jin+1)
            if (j1 > j2) then
! move particle only if it is below current hole
               do i = 1, idimp
                  part(i,j2) = part(i,j1)
               enddo
            endif
         enddo
         jin = jin + jsr(2)
         npp = npp - jsr(2)
      endif
      jss(1) = 0
! check if any particles have to be passed further
      if (ibflg(3) > 0) ibflg(3) = 1
      info(5) = max0(info(5),mter)
      if (ibflg(2) > 0) then
!        write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (iter < itermax) go to 60
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         return
      endif
! check if buffer overflowed and more particles remain to be checked
      if (ibflg(3) > 0) then
         nter = nter + 1
         info(4) = nter
         go to 10
      endif
      if (nter > 0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      end subroutine
!
      end module
!
! Make functions callable by Fortran77
!
!-----------------------------------------------------------------------
      subroutine PPINIT2(idproc,nvp)
      use pplib2, only: SUB => PPINIT2
      implicit none
      integer, intent(inout) :: idproc, nvp
      call SUB(idproc,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPEXIT
      use pplib2, only: SUB => PPEXIT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPABORT
      use pplib2, only: SUB => PPABORT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
      use pplib2, only: SUB => PWTIMERA
      implicit none
      integer, intent(in) :: icntrl
      real, intent(inout) :: time
      double precision, intent(inout) :: dtime
      call SUB(icntrl,time,dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPSUM(f,g,nxp)
      use pplib2, only: SUB => PPSUM
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSUM(f,g,nxp)
      use pplib2, only: SUB => PPDSUM
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPIMAX(if,ig,nxp)
      use pplib2, only: SUB => PPIMAX
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if, ig
      call SUB(if,ig,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDMAX(f,g,nxp)
      use pplib2, only: SUB => PPDMAX
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
      use pplib2, only: SUB => PPNCGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      call SUB(f,nyp,kstrt,nvp,nxv,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
      use pplib2, only: SUB => PPNAGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      real, dimension(nxv), intent(inout) :: scr
      call SUB(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
      use pplib2, only: SUB => PPNACGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx, ndim, nxv, nypmx
      real, dimension(ndim,nxv,nypmx), intent(inout) :: f
      real, dimension(ndim,nxv), intent(inout) :: scr
      call SUB(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd,  &
     &kypd)
      use pplib2, only: SUB => PPTPOSE
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv
      integer, intent(in) :: kxpd, kypd
      complex, dimension(nxv,kypd), intent(in) :: f
      complex, dimension(nyv,kxpd), intent(inout) :: g
      complex, dimension(kxp*kyp), intent(inout) :: s, t
      call SUB(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd,kypd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv, &
     &kxpd,kypd)
      use pplib2, only: SUB => PPNTPOSE
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      integer, intent(in) :: nxv, nyv, kxpd, kypd
      complex, dimension(ndim,nxv,kypd), intent(in) :: f
      complex, dimension(ndim,nyv,kxpd), intent(inout) :: g
      complex, dimension(ndim,kxp*kyp), intent(inout) :: s, t
      call SUB(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv,kxpd,kypd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,ny&
     &,kstrt,nvp,idimp,npmax,idps,nbmax,ntmax,info)
      use pplib2, only: SUB => PPMOVE2
      implicit none
      integer, intent(in) :: ny, kstrt, nvp, idimp, npmax, idps
      integer, intent(in) :: nbmax, ntmax
      integer, intent(inout) :: npp
      real, dimension(idimp,npmax), intent(inout) :: part
      real, dimension(idps), intent(in) :: edges
      real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr
      real, dimension(idimp,nbmax), intent(inout) :: rbufl, rbufr
      integer, dimension(ntmax+1), intent(inout) :: ihole
      integer, dimension(5), intent(inout) :: info
      call SUB(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,ny,kstrt,nvp&
     &,idimp,npmax,idps,nbmax,ntmax,info)
      end subroutine

