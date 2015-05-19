!-----------------------------------------------------------------------
! Basic parallel PIC library for MPI communications with OpenMP
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
! PPPMOVE2 moves particles into appropriate spatial regions for tiled
!          distributed data.
! written by viktor k. decyk, ucla
! copyright 1995, regents of the university of california
! update: may 9, 2015
      module mpplib2
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
      public :: PPTPOSE, PPNTPOSE, PPPMOVE2
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
      integer, intent(in) ::  nyp, kstrt, nvp, nx, ndim, nxv, nypmx
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
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, kyp
            do j = 1, kxp
               g(k,j) = f(j,k)
            enddo
         enddo
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, kyps
            do j = 1, ld
               s(j+ld*(k-1)) = f(j+joff,k)
            enddo
         enddo
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, ld
            do j = 1, kxps
               g(k+koff,j) = t(j+kxps*(k-1))
            enddo
         enddo
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, kyp
            do j = 1, kxp
               do i = 1, ndim
                  g(i,k,j) = f(i,j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, kyps
            do j = 1, ld
               do i = 1, ndim
                  s(i,j+ld*(k-1)) = f(i,j+joff,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, ld
            do j = 1, kxps
               do i = 1, ndim
                  g(i,k+koff,j) = t(i,j+kxps*(k-1))
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,  &
     &kstrt,nvp,idimp,nbmax,mx1)
! this subroutine moves particles into appropriate spatial regions
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! output: rbufr, rbufl, mcll, mclr
! sbufl = buffer for particles being sent to lower processor
! sbufr = buffer for particles being sent to upper processor
! rbufl = buffer for particles being received from lower processor
! rbufr = buffer for particles being received from upper processor
! ncll = particle number being sent to lower processor
! nclr = particle number being sent to upper processor
! mcll = particle number being received from lower processor
! mclr = particle number being received from upper processor
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 4 or 5
! nbmax =  size of buffers for passing particles between processors
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer, intent(in) :: kstrt, nvp, idimp, nbmax, mx1
      real, dimension(idimp,nbmax), intent(in) :: sbufl, sbufr
      real, dimension(idimp,nbmax), intent(inout) :: rbufl, rbufr
      integer, dimension(3,mx1), intent(in) :: ncll, nclr
      integer, dimension(3,mx1), intent(inout) :: mcll, mclr
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: ierr, ks, kl, kr, i, j, jsl, jsr
      integer :: nbsize, ncsize
      integer, dimension(8) :: msid
      integer, dimension(4) :: itg
      integer, dimension(lstat) :: istatus
      data itg /3,4,5,6/
      ks = kstrt - 1
      nbsize = idimp*nbmax
      ncsize = 3*mx1
! copy particle buffers: update rbufl, rbufr, mcll, mclr
! special case for one processor
      if (nvp==1) then
         do j = 1, mx1
            do i = 1, 3
               mcll(i,j) = nclr(i,j)
            enddo
         continue
         enddo
         do j = 1, mx1
            do i = 1, 3
               mclr(i,j) = ncll(i,j)
            enddo
         enddo
         do j = 1, nclr(3,mx1)
            do i = 1, idimp
               rbufl(i,j) = sbufr(i,j)
            enddo
         enddo
         do j = 1, ncll(3,mx1)
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
! post receives
         call MPI_IRECV(mcll,ncsize,mint,kl,itg(1),lgrp,msid(1),ierr)
         call MPI_IRECV(mclr,ncsize,mint,kr,itg(2),lgrp,msid(2),ierr)
         call MPI_IRECV(rbufl,nbsize,mreal,kl,itg(3),lgrp,msid(3),ierr)
         call MPI_IRECV(rbufr,nbsize,mreal,kr,itg(4),lgrp,msid(4),ierr)
! send particle number offsets
         call MPI_ISEND(nclr,ncsize,mint,kr,itg(1),lgrp,msid(5),ierr)
         call MPI_ISEND(ncll,ncsize,mint,kl,itg(2),lgrp,msid(6),ierr)
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_WAIT(msid(2),istatus,ierr)
! send particles
         jsr = idimp*nclr(3,mx1)
         call MPI_ISEND(sbufr,jsr,mreal,kr,itg(3),lgrp,msid(7),ierr)
         jsl = idimp*ncll(3,mx1)
         call MPI_ISEND(sbufl,jsl,mreal,kl,itg(4),lgrp,msid(8),ierr)
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
      endif
! make sure sbufr, sbufl, ncll, and nclr have been sent
      if (nvp /= 1) then
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
      use mpplib2, only: SUB => PPINIT2
      implicit none
      integer, intent(inout) :: idproc, nvp
      call SUB(idproc,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPEXIT
      use mpplib2, only: SUB => PPEXIT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPABORT
      use mpplib2, only: SUB => PPABORT
      implicit none
      call SUB()
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
      use mpplib2, only: SUB => PWTIMERA
      implicit none
      integer, intent(in) :: icntrl
      real, intent(inout) :: time
      double precision, intent(inout) :: dtime
      call SUB(icntrl,time,dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPSUM(f,g,nxp)
      use mpplib2, only: SUB => PPSUM
      implicit none
      integer, intent(in) :: nxp
      real, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSUM(f,g,nxp)
      use mpplib2, only: SUB => PPDSUM
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPIMAX(if,ig,nxp)
      use mpplib2, only: SUB => PPIMAX
      implicit none
      integer, intent(in) :: nxp
      integer, dimension(nxp), intent(inout) :: if, ig
      call SUB(if,ig,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDMAX(f,g,nxp)
      use mpplib2, only: SUB => PPDMAX
      implicit none
      integer, intent(in) :: nxp
      double precision, dimension(nxp), intent(inout) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
      use mpplib2, only: SUB => PPNCGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      call SUB(f,nyp,kstrt,nvp,nxv,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
      use mpplib2, only: SUB => PPNAGUARD2L
      implicit none
      integer, intent(in) :: nyp, kstrt, nvp, nx, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
      real, dimension(nxv), intent(inout) :: scr
      call SUB(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPNACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
      use mpplib2, only: SUB => PPNACGUARD2L
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
      use mpplib2, only: SUB => PPTPOSE
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
      use mpplib2, only: SUB => PPNTPOSE
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
      subroutine PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,  &
     &kstrt,nvp,idimp,nbmax,mx1)
      use mpplib2, only: SUB => PPPMOVE2
      implicit none
      integer, intent(in) :: kstrt, nvp, idimp, nbmax, mx1
      real, dimension(idimp,nbmax), intent(in) :: sbufl, sbufr
      real, dimension(idimp,nbmax), intent(inout) :: rbufl, rbufr
      integer, dimension(3,mx1), intent(in) :: ncll, nclr
      integer, dimension(3,mx1), intent(inout) :: mcll, mclr
      call SUB(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,nvp,   &
     &idimp,nbmax,mx1)
      end subroutine


