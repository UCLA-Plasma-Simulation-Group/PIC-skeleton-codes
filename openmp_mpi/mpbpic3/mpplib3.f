c-----------------------------------------------------------------------
c Basic parallel PIC library for MPI communications with OpenMP
c pplib3.f contains basic communications procedures for 2d partitions:
c PPINIT2 initializes parallel processing for Fortran90, returns
c         number of processors and processor id.
c PPEXIT terminates parallel processing.
c PPABORT aborts parallel processing.
c PWTIMERA performs parallel local wall clock timing.
c PPSUM performs parallel sum of a real vector.
c PPDSUM performs parallel sum of a double precision vector.
c PPIMAX performs parallel maximum of an integer vector.
c PPDMAX performs parallel maximum of a double precision vector.
c PPNCGUARD32L copies data to guard cells in y and z for scalar data,
c              linear interpolation, and distributed data with 2D
c              non-uniform partition.
c PPNAGUARD32L adds guard cells in y and z for scalar array, linear
c              interpolation, and distributed data with 2D non-uniform
c              partition.
c PPNACGUARD32L adds guard cells in y and z for vector array, linear
c               interpolation, and distributed data with 2D non-uniform
c               partition.
c PPTPOS3A performs a transpose of a complex scalar array, distributed
c          in y and z, to a complex scalar array, distributed in x and z
c PPTPOS3B performs a transpose of a complex scalar array, distributed
c          in x and z, to a complex scalar array, distributed in x and y
c PPNTPOS3A performs a transpose of an n component complex vector array,
c           distributed in y and z, to an n component complex vector
c           array, distributed in x and z.
c PPNTPOS3B performs a transpose of an n component complex vector array,
c           distributed in x and z, to an n component complex vector
c           array, distributed in x and y.
c PPPMOVE32 moves particles in y/z into appropriate spatial regions for
c           tiled distributed data with 2D spatial decomposition
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: october 27, 2015
c-----------------------------------------------------------------------
      function vresult(prec)
      implicit none
      real prec, vresult
      vresult = prec
      return
      end
c-----------------------------------------------------------------------
      function iresult(iprec)
      implicit none
      integer iprec, iresult
      iresult = iprec
      return
      end
c-----------------------------------------------------------------------
      subroutine PPINIT2(idproc,nvp)
c this subroutine initializes parallel processing
c lgrp communicator = MPI_COMM_WORLD
c output: idproc, nvp
c idproc = processor id in lgrp communicator
c nvp = number of real or virtual processors obtained
      implicit none
      integer idproc, nvp
c get definition of MPI constants
      include 'mpif.h'
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
      integer msum, mmax
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
c mint = default datatype for integers
c mcplx = default datatype for complex type
c mdouble = default double precision type
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c msum = MPI_SUM
c mmax = MPI_MAX
      common /PPARMSX/ msum, mmax
c local data
      integer ierror, ndprec, idprec
      integer ibig, iprec, iresult
      logical flag
      real small, prec, vresult
      save /PPARMS/, /PPARMSX/
      data small /1.0e-12/
      data ibig /2147483647/
      prec = 1.0 + small
      iprec = ibig + 1
c ndprec = (0,1) = (no,yes) use (normal,autodouble) precision
      if (vresult(prec).gt.1.0) then
         ndprec = 1
      else
         ndprec = 0
      endif
c idprec = (0,1) = (no,yes) use (normal,autodouble) integer precision
      if (iresult(iprec).gt.0) then
         idprec = 1
      else
         idprec = 0
      endif
c this segment is used for mpi computers
      if (MPI_STATUS_SIZE.gt.lstat) then
         write (2,*) ' status size too small, actual/required = ', lstat
     1, MPI_STATUS_SIZE
         stop
      endif
c indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (.not.flag) then
c initialize the MPI execution environment
         call MPI_INIT(ierror)
         if (ierror.ne.0) stop
      endif
      lworld = MPI_COMM_WORLD
      lgrp = lworld
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierror)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nproc,ierror)
c set default datatypes
      mint = MPI_INTEGER
      mdouble = MPI_DOUBLE_PRECISION
c single precision real
      if (ndprec.eq.0) then
         mreal = MPI_REAL
         mcplx = MPI_COMPLEX
c double precision real
      else
         mreal = MPI_DOUBLE_PRECISION
         mcplx = MPI_DOUBLE_COMPLEX
      endif
c single precision integer
c     if (idprec.eq.0) then
c        mint = MPI_INTEGER
c double precision integer
c     else
c        mint = MPI_INTEGER8
c     endif
c operators
      msum = MPI_SUM
      mmax = MPI_MAX
      nvp = nproc
      return
      end
c-----------------------------------------------------------------------
      subroutine PPEXIT
c this subroutine terminates parallel processing
      implicit none
c common block for parallel processing
      integer nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierror
      logical flag
c indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (flag) then
c synchronize processes
         call MPI_BARRIER(lworld,ierror)
c terminate MPI execution environment
         call MPI_FINALIZE(ierror)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPABORT
c this subroutine aborts parallel processing
      implicit none
c common block for parallel processing
      integer nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer errorcode, ierror
      logical flag
c indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (flag) then
         errorcode = 1
c terminate MPI execution environment
         call MPI_ABORT(lworld,errorcode,ierror)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
c this subroutine performs local wall clock timing
c input: icntrl, dtime
c icntrl = (-1,0,1) = (initialize,ignore,read) clock
c clock should be initialized before it is read!
c time = elapsed time in seconds
c dtime = current time
c written for mpi
      implicit none
      integer icntrl
      real time
      double precision dtime
c local data
      double precision jclock
      double precision MPI_WTIME
      external MPI_WTIME
c initialize clock
      if (icntrl.eq.(-1)) then
         dtime = MPI_WTIME()
c read clock and write time difference from last clock initialization
      else if (icntrl.eq.1) then
         jclock = dtime
         dtime = MPI_WTIME()
         time = real(dtime - jclock)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPSUM(f,g,nxp)
c this subroutine performs a parallel sum of a vector, that is:
c f(j,k) = sum over k of f(j,k)
c at the end, all processors contain the same summation.
c f = input and output real data
c g = scratch real array
c nxp = number of data values in vector
      implicit none
      real f, g
      integer nxp
      dimension f(nxp), g(nxp)
c common blocks for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
      integer msum, mmax
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c msum = MPI_SUM
      common /PPARMSX/ msum, mmax
c local data
      integer j, ierr
c perform sum
      call MPI_ALLREDUCE(f,g,nxp,mreal,msum,lgrp,ierr)
c copy output from scratch array
      do 10 j = 1, nxp
      f(j) = g(j)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPDSUM(f,g,nxp)
c this subroutine performs a parallel sum of a vector, that is:
c f(j,k) = sum over k of f(j,k)
c at the end, all processors contain the same summation.
c f = input and output double precision data
c g = scratch double precision array
c nxp = number of data values in vector
      implicit none
      double precision f, g
      integer nxp
      dimension f(nxp), g(nxp)
c common blocks for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
      integer msum, mmax
      parameter(lstat=10)
c lgrp = current communicator
c mdouble = default double precision type
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c msum = MPI_SUM
      common /PPARMSX/ msum, mmax
c local data
      integer j, ierr
c perform sum
      call MPI_ALLREDUCE(f,g,nxp,mdouble,msum,lgrp,ierr)
c copy output from scratch array
      do 10 j = 1, nxp
      f(j) = g(j)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPIMAX(if,ig,nxp)
c this subroutine finds parallel maximum for each element of a vector
c that is, if(j,k) = maximum as a function of k of if(j,k)
c at the end, all processors contain the same maximum.
c if = input and output integer data
c ig = scratch integer array
c nxp = number of data values in vector
      implicit none
      integer if, ig
      integer nxp
      dimension if(nxp), ig(nxp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
      integer msum, mmax
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c mmax = MPI_MAX
      common /PPARMSX/ msum, mmax
c local data
      integer j, ierr
c find maximum
      call MPI_ALLREDUCE(if,ig,nxp,mint,mmax,lgrp,ierr)
c copy output from scratch array
      do 10 j = 1, nxp
      if(j) = ig(j)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPDMAX(f,g,nxp)
c this subroutine finds parallel maximum for each element of a vector
c that is, f(j,k) = maximum as a function of k of f(j,k)
c at the end, all processors contain the same maximum.
c f = input and output double precision data
c g = scratch double precision array
c nxp = number of data values in vector
      implicit none
      double precision f, g
      integer nxp
      dimension f(nxp), g(nxp)
c common blocks for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
      integer msum, mmax
      parameter(lstat=10)
c lgrp = current communicator
c mdouble = default double precision type
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c mmax = MPI_MAX
      common /PPARMSX/ msum, mmax
c local data
      integer j, ierr
c find maximum
      call MPI_ALLREDUCE(f,g,nxp,mdouble,mmax,lgrp,ierr)
c copy output from scratch array
      do 10 j = 1, nxp
      f(j) = g(j)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPNCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx
     1,idds)
c this subroutine copies data to guard cells in non-uniform partitions
c f(j,k,l) = real data for grid j,k,l in particle partition.
c the grid is non-uniform and includes one extra guard cell.
c scs(j,k) = scratch array for particle partition
c nyzp(1:2) = number of primary gridpoints in y/z in particle partition
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data,
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, idds
      integer nyzp
      real f, scs
      dimension nyzp(idds)
      dimension f(nxv,nypmx,nzpmx), scs(nxv,nzpmx,2)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer j, k, js, ks, noff, kr, kl
      integer nxvz, nxvzs, nxvy, nxvys, nyp1, nzp1
      dimension istatus(lstat)
      nyp1 = nyzp(1) + 1
      nzp1 = nyzp(2) + 1
c js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      noff = nypmx*nzpmx
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do 20 k = 1, nyzp(2)
         do 10 j = 1, nxv
         f(j,nyp1,k) = f(j,1,k)
   10    continue
   20    continue
!$OMP END PARALLEL DO
         go to 70
      endif
c buffer data in y
!$OMP PARALLEL DO PRIVATE(j,k)
      do 40 k = 1, nyzp(2)
      do 30 j = 1, nxv
      scs(j,k,1) = f(j,1,k)
   30 continue
   40 continue
!$OMP END PARALLEL DO
c copy to guard cells in y
      nxvzs = nxv*nyzp(2)
      kr = js + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      kl = js - 1
      if (kl.lt.0) kl = kl + nvpy
      kr = kr + nvpy*ks
      kl = kl + nvpy*ks
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,2),nxvz,mreal,kr,noff+3,lgrp,msid,ierr)
      call MPI_SEND(scs,nxvzs,mreal,kl,noff+3,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c copy guard cells
!$OMP PARALLEL DO PRIVATE(j,k)
      do 60 k = 1, nyzp(2)
      do 50 j = 1, nxv
      f(j,nyp1,k) = scs(j,k,2)
   50 continue
   60 continue
!$OMP END PARALLEL DO
c special case for one processor in z
   70 if (nvpz.eq.1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do 90 k = 1, nyp1
         do 80 j = 1, nxv
         f(j,k,nzp1) = f(j,k,1)
   80    continue
   90    continue
!$OMP END PARALLEL DO
         return
      endif
c copy to guard cells in z
      nxvys = nxv*nyp1
      kr = ks + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = ks - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = js + nvpy*kr
      kl = js + nvpy*kl
c this segment is used for mpi computers
      call MPI_IRECV(f(1,1,nzp1),nxvy,mreal,kr,noff+4,lgrp,msid,ierr)
      call MPI_SEND(f,nxvys,mreal,kl,noff+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      return
      end
c-----------------------------------------------------------------------
      subroutine PPNAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,    
     1nypmx,nzpmx,idds)
c this subroutine adds data from guard cells in non-uniform partitions
c f(j,k,l) = real data for grid j,k,l in particle partition.
c the grid is non-uniform and includes one extra guard cell.
c scs/scr = scratch arrays for particle partition
c nyzp(1:2) = number of primary gridpoints in y/z in particle partition
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, idds
      integer nyzp
      real f, scs, scr
      dimension nyzp(idds)
      dimension f(nxv,nypmx,nzpmx), scs(nxv,nzpmx,2), scr(nxv,nypmx)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer j, k, js, ks, noff, kr, kl
      integer nx1, nxvz, nxvzs, nxvy, nxvys, nyp1, nzp1
      dimension istatus(lstat)
      nx1 = nx + 1
      nyp1 = nyzp(1) + 1
      nzp1 = nyzp(2) + 1
c js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      noff = nypmx*nzpmx
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do 20 k = 1, nzp1
         do 10 j = 1, nx1
         f(j,1,k) = f(j,1,k) + f(j,nyp1,k)
         f(j,nyp1,k) = 0.0
   10    continue
   20    continue
!$OMP END PARALLEL DO
         go to 70
      endif
c buffer data in y
!$OMP PARALLEL DO PRIVATE(j,k)
      do 40 k = 1, nzp1 
      do 30 j = 1, nxv
      scs(j,k,1) = f(j,nyp1,k)
   30 continue
   40 continue
!$OMP END PARALLEL DO
c add guard cells in y
      nxvzs = nxv*nzp1
      kr = js + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      kl = js - 1
      if (kl.lt.0) kl = kl + nvpy
      kr = kr + nvpy*ks
      kl = kl + nvpy*ks
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,2),nxvz,mreal,kl,noff+1,lgrp,msid,ierr)
      call MPI_SEND(scs,nxvzs,mreal,kr,noff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
!$OMP PARALLEL DO PRIVATE(j,k)
      do 60 k = 1, nzp1
      do 50 j = 1, nx1
      f(j,1,k) = f(j,1,k) + scs(j,k,2)
      f(j,nyp1,k) = 0.0
   50 continue
   60 continue
!$OMP END PARALLEL DO
c special case for one processor in z
   70 if (nvpz.eq.1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do 90 k = 1, nyp1
         do 80 j = 1, nx1
         f(j,k,1) = f(j,k,1) + f(j,k,nzp1)
         f(j,k,nzp1) = 0.0
   80    continue
   90    continue
!$OMP END PARALLEL DO
         return
      endif
c add guard cells in z
      nxvys = nxv*nyp1
      kr = ks + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = ks - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = js + nvpy*kr
      kl = js + nvpy*kl
c this segment is used for mpi computers
      call MPI_IRECV(scr,nxvy,mreal,kl,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,nzp1),nxvys,mreal,kr,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
!$OMP PARALLEL DO PRIVATE(j,k)
      do 110 k = 1, nyp1
      do 100 j = 1, nx1
      f(j,k,1) = f(j,k,1) + scr(j,k)
      f(j,k,nzp1) = 0.0
  100 continue
  110 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPNACGUARD32L(f,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,  
     1nxv,nypmx,nzpmx,idds)
c this subroutine adds data from guard cells in non-uniform partitions
c f(ndim,j,k,l) = real data for grid j,k,l in particle partition.
c the grid is non-uniform and includes one extra guard cell.
c scs/scr = scratch arrays for particle partition
c nyzp(1:2) = number of primary gridpoints in y/z in particle partition
c ndim = leading dimension of array f
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = second dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer ndim, kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, idds
      integer nyzp
      real f, scs, scr
      dimension nyzp(idds)
      dimension f(ndim,nxv,nypmx,nzpmx)
      dimension scs(ndim,nxv,nzpmx,2), scr(ndim,nxv,nypmx)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer i, j, k, js, ks, noff, kr, kl
      integer nx1, nxvz, nxvzs, nxvy, nxvys, nyp1, nzp1
      dimension istatus(lstat)
      nx1 = nx + 1
      nyp1 = nyzp(1) + 1
      nzp1 = nyzp(2) + 1
c js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      noff = ndim*nypmx*nzpmx
      nxvz = ndim*nxv*nzpmx
      nxvy = ndim*nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do 30 k = 1, nzp1
         do 20 j = 1, nx1
         do 10 i = 1, ndim
         f(i,j,1,k) = f(i,j,1,k) + f(i,j,nyp1,k)
         f(i,j,nyp1,k) = 0.0
   10    continue
   20    continue
   30    continue
!$OMP END PARALLEL DO
         go to 100
      endif
c buffer data in y
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do 60 k = 1, nzp1 
      do 50 j = 1, nxv
      do 40 i = 1, ndim
      scs(i,j,k,1) = f(i,j,nyp1,k)
   40 continue
   50 continue
   60 continue
!$OMP END PARALLEL DO
c add guard cells in y
      nxvzs = ndim*nxv*nzp1
      kr = js + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      kl = js - 1
      if (kl.lt.0) kl = kl + nvpy
      kr = kr + nvpy*ks
      kl = kl + nvpy*ks
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,1,2),nxvz,mreal,kl,noff+1,lgrp,msid,ierr)
      call MPI_SEND(scs,nxvzs,mreal,kr,noff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do 90 k = 1, nzp1
      do 80 j = 1, nx1
      do 70 i = 1, ndim
      f(i,j,1,k) = f(i,j,1,k) + scs(i,j,k,2)
      f(i,j,nyp1,k) = 0.0
   70 continue
   80 continue
   90 continue
!$OMP END PARALLEL DO
c special case for one processor in z
  100 if (nvpz.eq.1) then
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do 130 k = 1, nyp1
         do 120 j = 1, nx1
         do 110 i = 1, ndim
         f(i,j,k,1) = f(i,j,k,1) + f(i,j,k,nzp1)
         f(i,j,k,nzp1) = 0.0
  110    continue
  120    continue
  130    continue
!$OMP END PARALLEL DO
         return
      endif
c add guard cells in z
      nxvys = ndim*nxv*nyp1
      kr = ks + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = ks - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = js + nvpy*kr
      kl = js + nvpy*kl
c this segment is used for mpi computers
      call MPI_IRECV(scr,nxvy,mreal,kl,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,nzp1),nxvys,mreal,kr,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do 160 k = 1, nyp1
      do 150 j = 1, nx1
      do 140 i = 1, ndim
      f(i,j,k,1) = f(i,j,k,1) + scr(i,j,k)
      f(i,j,k,nzp1) = 0.0
  140 continue
  150 continue
  160 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPTPOS3A(f,g,s,t,nx,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,nxv, 
     1nyv,kxypd,kypd,kzpd)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(k+kyp*(m-1),j,l) = f(j+kxyp*(n-1),k,l), where
c 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
c 1 <= n <= nx/kxyp, 1 <= m <= ny/kyp
c and where indices n and m can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kxyp/kyp/kzp = number of data values per block in x/y/z
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nxv/nyv = first dimension of f/g
c kypd/kxypd = second dimension of f/g
c kzpd = third dimension of f and g
      implicit none
      integer nx, ny, nz, kxyp, kyp, kzp, kstrt, nvpy, nxv, nyv
      integer kxypd, kypd, kzpd
      complex f, g, s, t
      dimension f(nxv,kypd,kzpd), g(nyv,kxypd,kzpd)
      dimension s(kxyp*kyp*kzp), t(kxyp*kyp*kzp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer n, j, k, l, js, ks, kxyps, kyps, kzps, id, joff, koff, ld
      integer jd, kxyzp, ll
      integer ierr, msid, istatus
      dimension istatus(lstat)
c js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxyps = min(kxyp,max(0,nx-kxyp*js))
      kyps = min(kyp,max(0,ny-kyp*js))
      kzps = min(kzp,max(0,nz-kzp*ks))
      kxyzp = kxyp*kyp*kzp
c special case for one processor
      if (nvpy.eq.1) then
!$OMP PARALLEL DO PRIVATE(j,k,l,ll)
         do 20 ll = 1, kyp*kzp
         l = (ll - 1)/kyp
         k = ll - kyp*l
         l = l + 1
         do 10 j = 1, kxyp
         g(k,j,l) = f(j,k,l)
   10    continue
   20    continue
!$OMP END PARALLEL DO
         return
      endif
c this segment is used for shared memory computers
c     do 100 l = 1, nz
c     do 90 m = 1, min(ny,nvpy)
c     koff = kyp*(m - 1)
c     do 80 i = 1, min(nx,nvpy)
c     joff = kxyp*(i - 1)
c     do 70 k = 1, min(kyp,max(0,ny-koff))
c     do 60 j = 1, min(kxyp,max(0,nx-joff))
c     g(k+koff,j+joff,l) = f(j+joff,k+koff,l)
c  60 continue
c  70 continue
c  80 continue
c  90 continue
c 100 continue
c this segment is used for mpi computers
      do 70 n = 1, nvpy
      id = n - js - 1
      if (id.lt.0) id = id + nvpy
c extract data to send
      joff = kxyp*id
      ld = min(kxyp,max(0,nx-joff))
!$OMP PARALLEL DO PRIVATE(j,k,l,ll,koff)
      do 40 ll = 1, kyps*kzps
      l = (ll - 1)/kyps
      k = ll - kyps*l
      l = l + 1
      koff = kyps*(l - 1) - 1
      do 30 j = 1, ld
      s(j+ld*(k+koff)) = f(j+joff,k,l)
   30 continue
   40 continue
!$OMP END PARALLEL DO
      jd = id + nvpy*ks
      ld = ld*kyps*kzps
c post receive
      call MPI_IRECV(t,kxyzp,mcplx,jd,n,lgrp,msid,ierr)
c send data
      call MPI_SEND(s,ld,mcplx,jd,n,lgrp,ierr)
c receive data
      call MPI_WAIT(msid,istatus,ierr)
c insert data received
      koff = kyp*id
      ld = min(kyp,max(0,ny-koff))
!$OMP PARALLEL DO PRIVATE(j,k,l,ll,joff)
      do 60 ll = 1, ld*kzps
      l = (ll - 1)/ld
      k = ll - ld*l
      l = l + 1
      joff = ld*(l - 1) - 1
      do 50 j = 1, kxyps
      g(k+koff,j,l) = t(j+kxyps*(k+joff))
   50 continue
   60 continue
!$OMP END PARALLEL DO
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPTPOS3B(g,h,s,t,nx,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,nvpz
     1,nyv,nzv,kxypd,kyzpd,kzpd)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(l+kzp*(n-1),j,k) = g(k+kyzp*(m-1),j,l), where
c 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
c 1 <= m <= ny/kyzp, 1 <= n <= nz/kzp
c and where indices n and m can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c g = complex input array
c h = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kxyp/kyzp/kzp = number of data values per block in x/y/z
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nyv/nzv = first dimension of g/h
c kxypd = second dimension of g and h
c kzpd/kyzpd = third dimension of g/h
      implicit none
      integer nx, ny, nz, kxyp, kyzp, kzp, kstrt, nvpy, nvpz, nyv, nzv
      integer kxypd, kyzpd, kzpd
      complex g, h, s, t
      dimension g(nyv,kxypd,kzpd), h(nzv,kxypd,kyzpd)
      dimension s(kyzp*kxyp*kzp), t(kyzp*kxyp*kzp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer n, j, k, l, js, ks, kxyps, kyzps, kzps, id, koff, loff, ld
      integer jd, kxyzp, ll
      integer ierr, msid, istatus
      dimension istatus(lstat)
c js/ks = processor co-ordinates in x/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxyps = min(kxyp,max(0,nx-kxyp*js))
      kyzps = min(kyzp,max(0,ny-kyzp*ks))
      kzps = min(kzp,max(0,nz-kzp*ks))
      kxyzp = kxyp*kyzp*kzp
c special case for one processor
      if (nvpz.eq.1) then
!$OMP PARALLEL DO PRIVATE(j,k,l,ll)
         do 20 ll = 1, kxyp*kzp
         l = (ll - 1)/kxyp
         j = ll - kxyp*l
         l = l + 1
         do 10 k = 1, kyzp
         h(l,j,k) = g(k,j,l)
   10    continue
   20    continue
!$OMP END PARALLEL DO
         return
      endif
c this segment is used for shared memory computers
c     do 100 i = 1, min(nz,nvpz)
c     loff = kzp*(i - 1)
c     do 90 m = 1, min(ny,nvpz)
c     koff = kyzp*(m - 1)
c     do 80 l = 1, min(kzp,max(0,nz-loff))
c     do 70 j = 1, nx
c     do 60 k = 1, min(kyzp,max(0,ny-koff))
c     h(l+loff,j,k+koff) = g(k+koff,j,l+loff)
c  60 continue
c  70 continue
c  80 continue
c  90 continue
c 100 continue
c this segment is used for mpi computers
      do 70 n = 1, nvpz
      id = n - ks - 1
      if (id.lt.0) id = id + nvpz
c extract data to send
      koff = kyzp*id
      ld = min(kyzp,max(0,ny-koff))
!$OMP PARALLEL DO PRIVATE(j,k,l,ll,loff)
      do 40 ll = 1, kxyps*kzps
      l = (ll - 1)/kxyps
      j = ll - kxyps*l
      l = l + 1
      loff = kxyps*(l - 1) - 1
      do 30 k = 1, ld
      s(k+ld*(j+loff)) = g(k+koff,j,l)
   30 continue
   40 continue
!$OMP END PARALLEL DO
      jd = js + nvpy*id
      ld = ld*kxyps*kzps
c post receive
      call MPI_IRECV(t,kxyzp,mcplx,jd,n,lgrp,msid,ierr)
c send data
      call MPI_SEND(s,ld,mcplx,jd,n,lgrp,ierr)
c receive data
      call MPI_WAIT(msid,istatus,ierr)
c insert data received
      loff = kzp*id
      ld = min(kzp,max(0,nz-loff))
!$OMP PARALLEL DO PRIVATE(j,k,l,ll,koff)
      do 60 ll = 1, kxyps*ld
      l = (ll - 1)/kxyps
      j = ll - kxyps*l
      l = l + 1
      koff = kxyps*(l - 1) - 1
      do 50 k = 1, kyzps
      h(l+loff,j,k) = t(k+kyzps*(j+koff))
   50 continue
   60 continue
!$OMP END PARALLEL DO
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPNTPOS3A(f,g,s,t,nx,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,ndim
     1,nxv,nyv,kxypd,kypd,kzpd)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxyp*(n-1),k,l), where
c 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
c 1 <= n <= nx/kxyp, 1 <= m <= ny/kyp
c and where indices n and m can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kxyp/kyp/kzp = number of data values per block in x/y/z
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c ndim = leading dimension of arrays f and g
c nxv/nyv = first dimension of f/g
c kypd/kxypd = second dimension of f/g
c kzpd = third dimension of f and g
      implicit none
      integer nx, ny, nz, kxyp, kyp, kzp, kstrt, nvpy, ndim, nxv, nyv
      integer kxypd, kypd, kzpd
      complex f, g, s, t
      dimension f(ndim,nxv,kypd,kzpd), g(ndim,nyv,kxypd,kzpd)
      dimension s(ndim,kxyp*kyp*kzp), t(ndim,kxyp*kyp*kzp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer i, n, j, k, l, js, ks, kxyps, kyps, kzps, id, joff, koff
      integer ld, jd, kxyzp, ll
      integer ierr, msid, istatus
      dimension istatus(lstat)
c js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxyps = min(kxyp,max(0,nx-kxyp*js))
      kyps = min(kyp,max(0,ny-kyp*js))
      kzps = min(kzp,max(0,nz-kzp*ks))
      kxyzp = ndim*kxyp*kyp*kzp
c special case for one processor
      if (nvpy.eq.1) then
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll)
         do 30 ll = 1, kyp*kzp
         l = (ll - 1)/kyp
         k = ll - kyp*l
         l = l + 1
         do 20 j = 1, kxyp
         do 10 i = 1, ndim
         g(i,k,j,l) = f(i,j,k,l)
   10    continue
   20    continue
   30    continue
!$OMP END PARALLEL DO
         return
      endif
c this segment is used for shared memory computers
c     do 110 l = 1, nz
c     do 100 m = 1, min(ny,nvpy)
c     koff = kyp*(m - 1)
c     do 90 i = 1, min(nx,nvpy)
c     joff = kxyp*(i - 1)
c     do 80 k = 1, min(kyp,max(0,ny-koff))
c     do 70 j = 1, min(kxyp,max(0,nx-joff))
c     do 60 n = 1, ndim
c     g(n,k+koff,j+joff,l) = f(n,j+joff,k+koff,l)
c  60 continue
c  70 continue
c  80 continue
c  90 continue
c 100 continue
c 110 continue
c this segment is used for mpi computers
      do 100 n = 1, nvpy
      id = n - js - 1
      if (id.lt.0) id = id + nvpy
c extract data to send
      joff = kxyp*id
      ld = min(kxyp,max(0,nx-joff))
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,koff)
      do 60 ll = 1, kyps*kzps
      l = (ll - 1)/kyps
      k = ll - kyps*l
      l = l + 1
      koff = kyps*(l - 1) - 1
      do 50 j = 1, ld
      do 40 i = 1, ndim
      s(i,j+ld*(k+koff)) = f(i,j+joff,k,l)
   40 continue
   50 continue
   60 continue
!$OMP END PARALLEL DO
      jd = id + nvpy*ks
      ld = ndim*ld*kyps*kzps
c post receive
      call MPI_IRECV(t,kxyzp,mcplx,jd,n,lgrp,msid,ierr)
c send data
      call MPI_SEND(s,ld,mcplx,jd,n,lgrp,ierr)
c receive data
      call MPI_WAIT(msid,istatus,ierr)
c insert data received
      koff = kyp*id
      ld = min(kyp,max(0,ny-koff))
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,joff)
      do 90 ll = 1, ld*kzps
      l = (ll - 1)/ld
      k = ll - ld*l
      l = l + 1
      joff = ld*(l - 1) - 1
      do 80 j = 1, kxyps
      do 70 i = 1, ndim
      g(i,k+koff,j,l) = t(i,j+kxyps*(k+joff))
   70 continue
   80 continue
   90 continue
!$OMP END PARALLEL DO
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPNTPOS3B(g,h,s,t,nx,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,   
     1nvpz,ndim,nyv,nzv,kxypd,kyzpd,kzpd)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(1:ndim,l+kzp*(n-1),j,k) = g(1:ndim,k+kyzp*(m-1),j,l), where
c 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
c 1 <= m <= ny/kyzp, 1 <= n <= nz/kzp
c and where indices n and m can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c g = complex input array
c h = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kxyp/kyzp/kzp = number of data values per block in x/y/z
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c ndim = leading dimension of arrays g and h
c nyv/nzv = first dimension of g/h
c kxypd = second dimension of g and h
c kzpd/kyzpd = third dimension of g/h
      implicit none
      integer nx, ny, nz, kxyp, kyzp, kzp, kstrt, nvpy, nvpz, ndim
      integer nyv, nzv, kxypd, kyzpd, kzpd
      complex g, h, s, t
      dimension g(ndim,nyv,kxypd,kzpd), h(ndim,nzv,kxypd,kyzpd)
      dimension s(ndim,kyzp*kxyp*kzp), t(ndim,kyzp*kxyp*kzp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer i, n, j, k, l, js, ks, kxyps, kyzps, kzps, id, koff, loff
      integer ld, jd, kxyzp, ll
      integer ierr, msid, istatus
      dimension istatus(lstat)
c js/ks = processor co-ordinates in x/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxyps = min(kxyp,max(0,nx-kxyp*js))
      kyzps = min(kyzp,max(0,ny-kyzp*ks))
      kzps = min(kzp,max(0,nz-kzp*ks))
      kxyzp = ndim*kxyp*kyzp*kzp
c special case for one processor
      if (nvpz.eq.1) then
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll)
         do 30 ll = 1, kxyp*kzp
         l = (ll - 1)/kxyp
         j = ll - kxyp*l
         l = l + 1
         do 20 k = 1, kyzp
         do 10 i = 1, ndim
         h(i,l,j,k) = g(i,k,j,l)
   10    continue
   20    continue
   30    continue
!$OMP END PARALLEL DO
         return
      endif
c this segment is used for shared memory computers
c     do 100 i = 1, min(nz,nvpz)
c     loff = kzp*(i - 1)
c     do 100 m = 1, min(ny,nvpz)
c     koff = kyzp*(m - 1)
c     do 90 l = 1, min(kzp,max(0,nz-loff))
c     do 80 j = 1, nx
c     do 70 k = 1, min(kyzp,max(0,ny-koff))
c     do 60 n = 1, ndim
c     h(n,l+loff,j,k+koff) = g(n,k+koff,j,l+loff)
c  60 continue
c  70 continue
c  80 continue
c  90 continue
c 100 continue
c 110 continue
c this segment is used for mpi computers
      do 100 n = 1, nvpz
      id = n - ks - 1
      if (id.lt.0) id = id + nvpz
c extract data to send
      koff = kyzp*id
      ld = min(kyzp,max(0,ny-koff))
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,loff)
      do 60 ll = 1, kxyps*kzps
      l = (ll - 1)/kxyps
      j = ll - kxyps*l
      l = l + 1
      loff = kxyps*(l - 1) - 1
      do 50 k = 1, ld
      do 40 i = 1, ndim
      s(i,k+ld*(j+loff)) = g(i,k+koff,j,l)
   40 continue
   50 continue
   60 continue
!$OMP END PARALLEL DO
      jd = js + nvpy*id
      ld = ndim*ld*kxyps*kzps
c post receive
      call MPI_IRECV(t,kxyzp,mcplx,jd,n,lgrp,msid,ierr)
c send data
      call MPI_SEND(s,ld,mcplx,jd,n,lgrp,ierr)
c receive data
      call MPI_WAIT(msid,istatus,ierr)
c insert data received
      loff = kzp*id
      ld = min(kzp,max(0,nz-loff))
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,koff)
      do 90 ll = 1, kxyps*ld
      l = (ll - 1)/kxyps
      j = ll - kxyps*l
      l = l + 1
      koff = kxyps*(l - 1) - 1
      do 80 k = 1, kyzps
      do 70 i = 1, ndim
      h(i,l+loff,j,k) = t(i,k+kyzps*(j+koff))
   70 continue
   80 continue
   90 continue
!$OMP END PARALLEL DO
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPMOVE32(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr, 
     1mcls,kstrt,nvpy,nvpz,idimp,nbmax,mx1,myp1,mzp1,mxzyp1,irc)
c this subroutine moves particles into appropriate spatial regions in y
c for distributed data, with 2d domain decomposition in y/z.
c tiles are assumed to be arranged in 3D linear memory
c output: rbufr, rbufl, mcll, mclr
c sbufl = buffer for particles being sent to lower/back processor
c sbufr = buffer for particles being sent to upper/forward processor
c rbufl = buffer for particles being received from lower/back processor
c rbufr = buffer for particles being received from upper/forward
c processor
c ncll = particle number offsets sent to lower/back processor
c nclr = particle number offsets sent to upper/forward processor
c mcll = particle number offsets received from lower/back processor
c mclr = particle number offsets received from upper/forward processor
c mcls = particle number ofsets received from corner processors
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c idimp = size of phase space = 4 or 5
c nbmax =  size of buffers for passing particles between processors
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c mzp1 = (partition length in z direction - 1)/mz + 1
c mxzyp1 = mx1*max(myp1,mzp1)
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kstrt, nvpy, nvpz, idimp, nbmax, mx1, myp1, mzp1, mxzyp1
      integer irc
      real sbufr, sbufl, rbufr, rbufl
      integer ncll, nclr, mcll, mclr, mcls
      dimension sbufl(idimp,nbmax,2), sbufr(idimp,nbmax,2)
      dimension rbufl(idimp,nbmax,2), rbufr(idimp,nbmax,2)
      dimension ncll(3,mxzyp1,3,2), nclr(3,mxzyp1,3,2)
      dimension mcll(3,mxzyp1,3,2), mclr(3,mxzyp1,3,2), mcls(3,mx1+1,4)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, js, ks, kl, kr, i, j, k, n, jsl, jsr
      integer m, ll, lr, krr, krl, klr, kll, jsrr, jsrl, jslr, jsll
      integer nr, nl, mr, ml, nbr, nbl
      integer mxyp1, mxzp1, nbsize, ncsize, nsize
      integer nb, iwork
      integer msid, itg, istatus
      dimension msid(8), itg(12), istatus(lstat)
      dimension nb(1), iwork(1)
      data itg /3,4,5,6,7,8,9,10,11,12,13,14/
c js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      mxyp1 = mx1*myp1
      mxzp1 = mx1*mzp1
      nbsize = idimp*nbmax
      ncsize = 9*mxzyp1
c copy particle buffers in y:
c update rbufl(:,1), rbufr(:,1), mcll(:,1), mclr(:,1)
c special case for one processor
      if ((nvpy*nvpz).eq.1) then
!$OMP PARALLEL
!$OMP DO PRIVATE(i,j,k,ll)
         do 20 ll = 1, 3*mxzp1
         k = (ll - 1)/mxzp1
         j = ll - mxzp1*k
         k = k + 1
         do 10 i = 1, 3
         mcll(i,j,k,1) = nclr(i,j,k,1)
         mclr(i,j,k,1) = ncll(i,j,k,1)
   10    continue
   20    continue
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i,j)
         do 40 j = 1, nclr(3,mxzp1,3,1)
         do 30 i = 1, idimp
         rbufl(i,j,1) = sbufr(i,j,1)
   30    continue
   40    continue
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i,j)
         do 60 j = 1, ncll(3,mxzp1,3,1)
         do 50 i = 1, idimp
         rbufr(i,j,1) = sbufl(i,j,1)
   50    continue
   60    continue
!$OMP END DO
!$OMP END PARALLEL
c get particles from corners
         n = mx1*(mzp1 - 1)
c zero out base addresses in prefix scans
         if (n.gt.0) then
            nr = nclr(3,n,3,1)
            nl = ncll(3,n,3,1)
         else
            nr = nclr(3,mxzp1,2,1)
            nl = ncll(3,mxzp1,2,1)
         endif
         do 80 j = 1, mx1
         do 70 i = 1, 3
         nclr(i,j,2,1) = nclr(i,j,2,1) - nclr(3,mxzp1,1,1)
         nclr(i,n+j,3,1) = nclr(i,n+j,3,1) - nr
         ncll(i,j,2,1) = ncll(i,j,2,1) - ncll(3,mxzp1,1,1)
         ncll(i,n+j,3,1) = ncll(i,n+j,3,1) - nl
   70    continue
   80    continue
c add new base addresses in prefix scans
         ml = mcll(3,mxzp1,3,1)
         mr = mclr(3,mxzp1,3,1)
         do 100 j = 1, mx1
         do 90 i = 1, 3
         mcls(i,j,1) = nclr(i,j,2,1) + ml
         mcls(i,j,3) = ncll(i,j,2,1) + mr
   90    continue
  100    continue
         mcls(1,mx1+1,1) = ml
         mcls(1,mx1+1,3) = mr
c append corner particles to end of buffers
         k = nclr(3,mx1,2,1)
         m = nclr(3,mxzp1,1,1)
         do 120 j = 1, k
         do 110 i = 1, idimp
         rbufl(i,j+ml,1) = sbufr(i,j+m,1)
  110    continue
  120    continue
         ml = ml + k
         k = ncll(3,mx1,2,1)
         m = ncll(3,mxzp1,1,1)
         do 140 j = 1, k
         do 130 i = 1, idimp
         rbufr(i,j+mr,1) = sbufl(i,j+m,1)
  130    continue
  140    continue
         mr = mr + k
c add new base addresses in prefix scans
         do 160 j = 1, mx1
         do 150 i = 1, 3
         mcls(i,j,2) = nclr(i,n+j,3,1) + ml
         mcls(i,j,4) = ncll(i,n+j,3,1) + mr
  150    continue
  160    continue
         mcls(1,mx1+1,2) = ml
         mcls(1,mx1+1,4) = mr
c append more corner particles to end of buffers
         do 180 j = 1, nclr(3,n+mx1,3,1)
         do 170 i = 1, idimp
         rbufl(i,j+ml,1) = sbufr(i,j+nr,1)
  170    continue
  180    continue
         do 200 j = 1, ncll(3,n+mx1,3,1)
         do 190 i = 1, idimp
         rbufr(i,j+mr,1) = sbufl(i,j+nl,1)
  190    continue
  200    continue
c this segment is used for mpi computers
      else
c get particles from below and above
         kr = js + 1
         if (kr.ge.nvpy) kr = kr - nvpy
         kl = js - 1
         if (kl.lt.0) kl = kl + nvpy
         kr = kr + nvpy*ks
         kl = kl + nvpy*ks
c post receives
         call MPI_IRECV(mcll(1,1,1,1),ncsize,mint,kl,itg(1),lgrp,msid(1)
     1,ierr)
         call MPI_IRECV(mclr(1,1,1,1),ncsize,mint,kr,itg(2),lgrp,msid(2)
     1,ierr)
         call MPI_IRECV(rbufl(1,1,1),nbsize,mreal,kl,itg(3),lgrp,msid(3)
     1,ierr)
         call MPI_IRECV(rbufr(1,1,1),nbsize,mreal,kr,itg(4),lgrp,msid(4)
     1,ierr)
c send particle number offsets
         call MPI_ISEND(nclr(1,1,1,1),ncsize,mint,kr,itg(1),lgrp,msid(5)
     1,ierr)
         call MPI_ISEND(ncll(1,1,1,1),ncsize,mint,kl,itg(2),lgrp,msid(6)
     1,ierr)
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_WAIT(msid(2),istatus,ierr)
c send particles
         jsr = idimp*nclr(3,mxzp1,3,1)
         call MPI_ISEND(sbufr(1,1,1),jsr,mreal,kr,itg(3),lgrp,msid(7),  
     1ierr)
         jsl = idimp*ncll(3,mxzp1,3,1)
         call MPI_ISEND(sbufl(1,1,1),jsl,mreal,kl,itg(4),lgrp,msid(8),  
     1ierr)
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
c make sure sbufr, sbufl, ncll, and nclr have been sent
         do 210 i = 1, 4
         call MPI_WAIT(msid(i+4),istatus,ierr)
  210    continue
c get particles from corners
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
c zero out base addresses in prefix scans
         if (n.gt.0) then
            nr = nclr(3,n,3,1)
            nl = ncll(3,n,3,1)
         else
            nr = nclr(3,mxzp1,2,1)
            nl = ncll(3,mxzp1,2,1)
         endif
         do 230 j = 1, mx1
         do 220 i = 1, 3
         nclr(i,j,2,1) = nclr(i,j,2,1) - nclr(3,mxzp1,1,1)
         nclr(i,n+j,3,1) = nclr(i,n+j,3,1) - nr
         ncll(i,j,2,1) = ncll(i,j,2,1) - ncll(3,mxzp1,1,1)
         ncll(i,n+j,3,1) = ncll(i,n+j,3,1) - nl
  220    continue
  230    continue
         n = n + 1
c post receives
         call MPI_IRECV(mcls(1,1,1),nsize,mint,klr,itg(5),lgrp,msid(1), 
     1ierr)
         call MPI_IRECV(mcls(1,1,2),nsize,mint,kll,itg(6),lgrp,msid(2), 
     1ierr)
         call MPI_IRECV(mcls(1,1,3),nsize,mint,krr,itg(7),lgrp,msid(3), 
     1ierr)
         call MPI_IRECV(mcls(1,1,4),nsize,mint,krl,itg(8),lgrp,msid(4), 
     1ierr)
c send particle number offsets
         call MPI_ISEND(nclr(1,1,2,1),nsize,mint,krl,itg(5),lgrp,msid(5)
     1,ierr)
         call MPI_ISEND(nclr(1,n,3,1),nsize,mint,krr,itg(6),lgrp,msid(6)
     1,ierr)
         call MPI_ISEND(ncll(1,1,2,1),nsize,mint,kll,itg(7),lgrp,msid(7)
     1,ierr)
         call MPI_ISEND(ncll(1,n,3,1),nsize,mint,klr,itg(8),lgrp,msid(8)
     1,ierr)
c make sure particle offsets have been sent to and received from corners
         do 240 i = 1, 8
         call MPI_WAIT(msid(i),istatus,ierr)
  240    continue
c check for overflow errors
         ml = mcll(3,mxzp1,3,1)
         mr = mclr(3,mxzp1,3,1)
         nbl = nbmax - (ml + (mcls(3,mx1,1) + mcls(3,mx1,2)))
         nbr = nbmax - (mr + (mcls(2,mx1,3) + mcls(3,mx1,4)))
         nb(1) = min(-nbl,-nbr)
         call PPIMAX(nb,iwork,1)
         if (nb(1).gt.0) then
            write (2,*) 'corner buffer overflow error = ', nb(1)
            irc = nb(1)
            return
         endif
         nbl = idimp*nbl
         nbr = idimp*nbr
c add new base addresses in prefix scans
         do 260 j = 1, mx1
         do 250 i = 1, 3
         mcls(i,j,1) = mcls(i,j,1) + ml
         mcls(i,j,3) = mcls(i,j,3) + mr
  250    continue
  260    continue
         mcls(1,mx1+1,1) = ml
         mcls(1,mx1+1,3) = mr
c post first part of particle receives, append to end
         ml = ml + 1
         call MPI_IRECV(rbufl(1,ml,1),nbl,mreal,klr,itg(9),lgrp,msid(1),
     1ierr)
         mr = mr + 1
         call MPI_IRECV(rbufr(1,mr,1),nbr,mreal,krr,itg(11),lgrp,msid(3)
     1,ierr)
c send first part of particles
         m = nclr(3,mxzp1,1,1) + 1
         jsrl = idimp*nclr(3,mx1,2,1)
         call MPI_ISEND(sbufr(1,m,1),jsrl,mreal,krl,itg(9),lgrp,msid(5),
     1ierr)
         m = ncll(3,mxzp1,1,1) + 1
         jsll = idimp*ncll(3,mx1,2,1)
         call MPI_ISEND(sbufl(1,m,1),jsll,mreal,kll,itg(11),lgrp,msid(7)
     1,ierr)
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
c add new base addresses in prefix scans
         do 280 j = 1, mx1
         do 270 i = 1, 3
         mcls(i,j,2) = mcls(i,j,2) + ml
         mcls(i,j,4) = mcls(i,j,4) + mr
  270    continue
  280    continue
         mcls(1,mx1+1,2) = ml
         mcls(1,mx1+1,4) = mr
c post second part of particle receives, append to end
         ml = ml + 1
         call MPI_IRECV(rbufl(1,ml,1),nbl,mreal,kll,itg(10),lgrp,msid(2)
     1,ierr)
         mr = mr + 1
         call MPI_IRECV(rbufr(1,mr,1),nbr,mreal,krl,itg(12),lgrp,msid(4)
     1,ierr)
c send second part of particles
         jsrr = idimp*nclr(3,n+mx1-1,3,1)
         m = nr + 1
         call MPI_ISEND(sbufr(1,m,1),jsrr,mreal,krr,itg(10),lgrp,msid(6)
     1,ierr)
         jslr = idimp*ncll(3,n+mx1-1,3,1)
         m = nl + 1
         call MPI_ISEND(sbufl(1,m,1),jslr,mreal,klr,itg(12),lgrp,msid(8)
     1,ierr)
         call MPI_WAIT(msid(2),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
c make sure sbufr and sbufl have been sent to corners
         do 290 i = 1, 4
         call MPI_WAIT(msid(i+4),istatus,ierr)
  290    continue
      endif
c copy particle buffers in z:
c update rbufl(:,2), rbufr(:,2), mcll(:,2), mclr(:,2)
c special case for one processor
      if ((nvpy*nvpz).eq.1) then
!$OMP PARALLEL
!$OMP DO PRIVATE(i,j,k,ll)
         do 310 ll = 1, 3*mxyp1
         k = (ll - 1)/mxyp1
         j = ll - mxyp1*k
         k = k + 1
         do 300 i = 1, 3
         mcll(i,j,k,2) = nclr(i,j,k,2)
         mclr(i,j,k,2) = ncll(i,j,k,2)
  300    continue
  310    continue
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i,j)
         do 330 j = 1, nclr(3,mxyp1,3,2)
         do 320 i = 1, idimp
         rbufl(i,j,2) = sbufr(i,j,2)
  320    continue
  330    continue
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i,j)
         do 350 j = 1, ncll(3,mxyp1,3,2)
         do 340 i = 1, idimp
         rbufr(i,j,2) = sbufl(i,j,2)
  340    continue
  350    continue
!$OMP END DO
!$OMP END PARALLEL
c this segment is used for mpi computers
      else
c get particles from back and front
         kr = ks + 1
         if (kr.ge.nvpz) kr = kr - nvpz
         kl = ks - 1
         if (kl.lt.0) kl = kl + nvpz
         kr = js + nvpy*kr
         kl = js + nvpy*kl
c post receives
         call MPI_IRECV(mcll(1,1,1,2),ncsize,mint,kl,itg(1),lgrp,msid(1)
     1,ierr)
         call MPI_IRECV(mclr(1,1,1,2),ncsize,mint,kr,itg(2),lgrp,msid(2)
     1,ierr)
         call MPI_IRECV(rbufl(1,1,2),nbsize,mreal,kl,itg(3),lgrp,msid(3)
     1,ierr)
         call MPI_IRECV(rbufr(1,1,2),nbsize,mreal,kr,itg(4),lgrp,msid(4)
     1,ierr)
c send particle number offsets
         call MPI_ISEND(nclr(1,1,1,2),ncsize,mint,kr,itg(1),lgrp,msid(5)
     1,ierr)
         call MPI_ISEND(ncll(1,1,1,2),ncsize,mint,kl,itg(2),lgrp,msid(6)
     1,ierr)
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_WAIT(msid(2),istatus,ierr)
c send particles
         jsr = idimp*nclr(3,mxyp1,3,2)
         call MPI_ISEND(sbufr(1,1,2),jsr,mreal,kr,itg(3),lgrp,msid(7),  
     1ierr)
         jsl = idimp*ncll(3,mxyp1,3,2)
         call MPI_ISEND(sbufl(1,1,2),jsl,mreal,kl,itg(4),lgrp,msid(8),  
     1ierr)
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
c make sure sbufr, sbufl, ncll, and nclr have been sent
         do 360 i = 1, 4
         call MPI_WAIT(msid(i+4),istatus,ierr)
  360    continue
      endif
      return
      end
