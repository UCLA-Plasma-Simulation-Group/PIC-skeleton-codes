c-----------------------------------------------------------------------
c Basic parallel PIC library for MPI communications
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
c PPMOVE32 moves particles into appropriate spatial regions with
c          periodic boundary conditions and 2D spatial decomposition.
c          Assumes ihole list has been found.
c PPMOVEG32 moves particles into appropriate spatial regions with
c           periodic boundary conditions and 2D spatial decomposition.
c           ihole list is calculated from particles co-ordinates
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: august 29, 2015
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
         do 20 k = 1, nyzp(2)
         do 10 j = 1, nxv
         f(j,nyp1,k) = f(j,1,k)
   10    continue
   20    continue
         go to 70
      endif
c buffer data in y
      do 40 k = 1, nyzp(2)
      do 30 j = 1, nxv
      scs(j,k,1) = f(j,1,k)
   30 continue
   40 continue
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
      do 60 k = 1, nyzp(2)
      do 50 j = 1, nxv
      f(j,nyp1,k) = scs(j,k,2)
   50 continue
   60 continue
c special case for one processor in z
   70 if (nvpz.eq.1) then
         do 90 k = 1, nyp1
         do 80 j = 1, nxv
         f(j,k,nzp1) = f(j,k,1)
   80    continue
   90    continue
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
         do 20 k = 1, nzp1
         do 10 j = 1, nx1
         f(j,1,k) = f(j,1,k) + f(j,nyp1,k)
         f(j,nyp1,k) = 0.0
   10    continue
   20    continue
         go to 70
      endif
c buffer data in y
      do 40 k = 1, nzp1 
      do 30 j = 1, nxv
      scs(j,k,1) = f(j,nyp1,k)
   30 continue
   40 continue
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
      do 60 k = 1, nzp1
      do 50 j = 1, nx1
      f(j,1,k) = f(j,1,k) + scs(j,k,2)
      f(j,nyp1,k) = 0.0
   50 continue
   60 continue
c special case for one processor in z
   70 if (nvpz.eq.1) then
         do 90 k = 1, nyp1
         do 80 j = 1, nx1
         f(j,k,1) = f(j,k,1) + f(j,k,nzp1)
         f(j,k,nzp1) = 0.0
   80    continue
   90    continue
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
      do 110 k = 1, nyp1
      do 100 j = 1, nx1
      f(j,k,1) = f(j,k,1) + scr(j,k)
      f(j,k,nzp1) = 0.0
  100 continue
  110 continue
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
         do 30 k = 1, nzp1
         do 20 j = 1, nx1
         do 10 i = 1, ndim
         f(i,j,1,k) = f(i,j,1,k) + f(i,j,nyp1,k)
         f(i,j,nyp1,k) = 0.0
   10    continue
   20    continue
   30    continue
         go to 100
      endif
c buffer data in y
      do 60 k = 1, nzp1 
      do 50 j = 1, nxv
      do 40 i = 1, ndim
      scs(i,j,k,1) = f(i,j,nyp1,k)
   40 continue
   50 continue
   60 continue
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
      do 90 k = 1, nzp1
      do 80 j = 1, nx1
      do 70 i = 1, ndim
      f(i,j,1,k) = f(i,j,1,k) + scs(i,j,k,2)
      f(i,j,nyp1,k) = 0.0
   70 continue
   80 continue
   90 continue
c special case for one processor in z
  100 if (nvpz.eq.1) then
         do 130 k = 1, nyp1
         do 120 j = 1, nx1
         do 110 i = 1, ndim
         f(i,j,k,1) = f(i,j,k,1) + f(i,j,k,nzp1)
         f(i,j,k,nzp1) = 0.0
  110    continue
  120    continue
  130    continue
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
      do 160 k = 1, nyp1
      do 150 j = 1, nx1
      do 140 i = 1, ndim
      f(i,j,k,1) = f(i,j,k,1) + scr(i,j,k)
      f(i,j,k,nzp1) = 0.0
  140 continue
  150 continue
  160 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNACGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypm
     1x,nzpmx,mblok,nblok,ngds,idds)
c this subroutine adds data from guard cells in non-uniform partitions
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes one extra guard cell.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds
      integer nyzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(3,nxv,nypmx,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, kr, kl, mnblok
      integer nx1, nxvz, nxvzs, nyzp1, nxvy, nxvys, m, my, mz, j, k, n
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
         do 40 m = 1, mnblok
         do 30 k = 1, nyzp(2,m)+1
         do 20 j = 1, nx1
         do 10 n = 1, 3
         f(n,j,1,k,m) = f(n,j,1,k,m) + f(n,j,nyzp(1,m)+1,k,m)
         f(n,j,nyzp(1,m)+1,k,m) = 0.
   10    continue
   20    continue
   30    continue
   40    continue
         go to 170
      endif
c buffer data in y
      do 80 m = 1, mnblok
      do 70 k = 1, nyzp(2,m)+1
      do 60 j = 1, nxv
      do 50 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,nyzp(1,m)+1,k,m)
   50 continue
   60 continue
   70 continue
   80 continue
c add guard cells in y
      do 160 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 150 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      nxvzs = nxv*nyzp1
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      kl = ky - 1
      if (kl.lt.0) kl = kl + nvpy
      kr = kr + kz
      kl = kl + kz
c this segment is used for shared memory computers
c     do 110 k = 1, nyzp1
c     do 100 j = 1, nxv
c     do 90 n = 1, 3
c     scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c  90 continue
c 100 continue
c 110 continue
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl-1,noff+1,lgrp,msid,i
     1err)
      call MPI_SEND(scs(1,1,1,1,m),3*nxvzs,mreal,kr-1,noff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 140 k = 1, nyzp1
      do 130 j = 1, nx1
      do 120 n = 1, 3
      f(n,j,1,k,m) = f(n,j,1,k,m) + scs(n,j,k,2,m)
      f(n,j,nyzp(1,m)+1,k,m) = 0.
  120 continue
  130 continue
  140 continue
  150 continue
  160 continue
c special case for one processor in z
  170 if (nvpz.eq.1) then
         do 210 m = 1, mnblok
         do 200 k = 1, nyzp(1,m)+1
         do 190 j = 1, nx1
         do 180 n = 1, 3
         f(n,j,k,1,m) = f(n,j,k,1,m) + f(n,j,k,nyzp(2,m)+1,m)
         f(n,j,k,nyzp(2,m)+1,m) = 0.
  180    continue
  190    continue
  200    continue
  210    continue
         return
      endif
c add guard cells in z
      do 290 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 280 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(1,m) + 1
      nxvys = nxv*nyzp1
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 240 k = 1, nyzp1
c     do 230 j = 1, nxv
c     do 220 n = 1, 3
c     scr(n,j,k,m) = f(n,j,k,nyzp(2,m)+1,kl)
c 220 continue
c 230 continue
c 240 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,3*nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,nyzp(2,m)+1,m),3*nxvys,mreal,kr-1,noff+2,lgr
     1p,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 270 k = 1, nyzp1
      do 260 j = 1, nx1
      do 250 n = 1, 3
      f(n,j,k,1,m) = f(n,j,k,1,m) + scr(n,j,k,m)
      f(n,j,k,nyzp(2,m)+1,m) = 0.
  250 continue
  260 continue
  270 continue
  280 continue
  290 continue
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
      integer jd, kxyzp
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
         do 30 l = 1, kzp
         do 20 k = 1, kyp
         do 10 j = 1, kxyp
         g(k,j,l) = f(j,k,l)
   10    continue
   20    continue
   30    continue
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
      do 100 n = 1, nvpy
      id = n - js - 1
      if (id.lt.0) id = id + nvpy
c extract data to send
      joff = kxyp*id
      ld = min(kxyp,max(0,nx-joff))
      do 60 l = 1, kzps
      koff = kyps*(l - 1) - 1
      do 50 k = 1, kyps
      do 40 j = 1, ld
      s(j+ld*(k+koff)) = f(j+joff,k,l)
   40 continue
   50 continue
   60 continue
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
      do 90 l = 1, kzps
      joff = ld*(l - 1) - 1
      do 80 k = 1, ld
      do 70 j = 1, kxyps
      g(k+koff,j,l) = t(j+kxyps*(k+joff))
   70 continue
   80 continue
   90 continue
  100 continue
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
      integer jd, kxyzp
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
         do 30 l = 1, kzp
         do 20 j = 1, kxyp
         do 10 k = 1, kyzp
         h(l,j,k) = g(k,j,l)
   10    continue
   20    continue
   30    continue
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
      do 100 n = 1, nvpz
      id = n - ks - 1
      if (id.lt.0) id = id + nvpz
c extract data to send
      koff = kyzp*id
      ld = min(kyzp,max(0,ny-koff))
      do 60 l = 1, kzps
      loff = kxyps*(l - 1) - 1
      do 50 j = 1, kxyps
      do 40 k = 1, ld
      s(k+ld*(j+loff)) = g(k+koff,j,l)
   40 continue
   50 continue
   60 continue
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
      do 90 l = 1, ld
      koff = kxyps*(l - 1) - 1
      do 80 j = 1, kxyps
      do 70 k = 1, kyzps
      h(l+loff,j,k) = t(k+kyzps*(j+koff))
   70 continue
   80 continue
   90 continue
  100 continue
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
      integer ld, jd, kxyzp
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
         do 40 l = 1, kzp
         do 30 k = 1, kyp
         do 20 j = 1, kxyp
         do 10 i = 1, ndim
         g(i,k,j,l) = f(i,j,k,l)
   10    continue
   20    continue
   30    continue
   40    continue
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
      do 130 n = 1, nvpy
      id = n - js - 1
      if (id.lt.0) id = id + nvpy
c extract data to send
      joff = kxyp*id
      ld = min(kxyp,max(0,nx-joff))
      do 80 l = 1, kzps
      koff = kyps*(l - 1) - 1
      do 70 k = 1, kyps
      do 60 j = 1, ld
      do 50 i = 1, ndim
      s(i,j+ld*(k+koff)) = f(i,j+joff,k,l)
   50 continue
   60 continue
   70 continue
   80 continue
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
      do 120 l = 1, kzps
      joff = ld*(l - 1) - 1
      do 110 k = 1, ld
      do 100 j = 1, kxyps
      do 90 i = 1, ndim
      g(i,k+koff,j,l) = t(i,j+kxyps*(k+joff))
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
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
      integer ld, jd, kxyzp
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
         do 40 l = 1, kzp
         do 30 j = 1, kxyp
         do 20 k = 1, kyzp
         do 10 i = 1, ndim
         h(i,l,j,k) = g(i,k,j,l)
   10    continue
   20    continue
   30    continue
   40    continue
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
      do 130 n = 1, nvpz
      id = n - ks - 1
      if (id.lt.0) id = id + nvpz
c extract data to send
      koff = kyzp*id
      ld = min(kyzp,max(0,ny-koff))
      do 80 l = 1, kzps
      loff = kxyps*(l - 1) - 1
      do 70 j = 1, kxyps
      do 60 k = 1, ld
      do 50 i = 1, ndim
      s(i,k+ld*(j+loff)) = g(i,k+koff,j,l)
   50 continue
   60 continue
   70 continue
   80 continue
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
      do 120 l = 1, ld
      koff = kxyps*(l - 1) - 1
      do 110 j = 1, kxyps
      do 100 k = 1, kyzps
      do 90 i = 1, ndim
      h(i,l+loff,j,k) = t(i,k+kyzps*(j+koff))
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole, 
     1ny,nz,kstrt,nvpy,nvpz,idimp,npmax,idps,nbmax,ntmax,info)
c this subroutine moves particles into appropriate spatial regions
c ihole array is calculated in particle push procedure
c with periodic boundary conditions and 2D spatial decomposition
c output: part, ihole, npp, sbufr, sbufl, rbufr, rbufl, info
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = position z of particle n in partition
c part(4,n) = velocity vx of particle n in partition
c part(5,n) = velocity vy of particle n in partition
c part(6,n) = velocity vz of particle n in partition m
c edges(1:2) = lower/upper boundary in y of particle partition
c edges(3:4) = back/front boundary in z of particle partition
c npp = number of particles in partition
c sbufl = buffer for particles being sent to back processor
c sbufr = buffer for particles being sent to front processor
c rbufl = buffer for particles being received from back processor
c rbufr = buffer for particles being received from front processor
c ihole = location of holes left in particle arrays
c ny/nz = system length in y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition.
c idps = number of particle partition boundaries = 4
c nbmax =  size of buffers for passing particles between processors
c ntmax =  size of hole array for particles leaving processors
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4:5) = maximum number of buffer overflows in y/z
c info(6:7) = maximum number of particle passes required in y/z
      implicit none
      integer npp, ny, nz, kstrt, nvpy, nvpz, idimp, npmax, idps, nbmax
      integer ntmax
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer ihole, info
      dimension part(idimp,npmax)
      dimension edges(idps)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension rbufl(idimp,nbmax), rbufr(idimp,nbmax)
      dimension ihole(ntmax+1,2)
      dimension info(7)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
c iy/iz = partitioned co-ordinates
      integer iy, iz
      parameter(iy=2,iz=3)
      integer i, j, n, js, ks, ic, nvp, iter, nps, itg, kl, kr, j1, j2
      integer ih, jh, nh, j3, joff, jin, nbsize, nter, mter, itermax
      integer ierr
      integer msid, istatus
      integer kb, jsl, jsr, jss, ibflg, iwork
      real an, xt
      dimension msid(4), istatus(lstat)
      dimension kb(2), jsl(2), jsr(2), jss(2), ibflg(5), iwork(5)
c js/ks = processor co-ordinates in y/z=> idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      nbsize = idimp*nbmax
      info(1) = 0
      info(6) = 0
      info(7) = 0
      itermax = 2000
c buffer outgoing particles, first in y then in z direction
      do 320 n = 1, 2
      if (n.eq.1) then
         ic = iy
         nvp = nvpy
         an = real(ny)
         jh = ihole(1,n+1)
         ibflg(5) = 0
      else if (n.eq.2) then
         ic = iz
         nvp = nvpz
         an = real(nz)
      endif
      iter = 2
      nter = 0
      ih = ihole(1,n)
      joff = 1
      jin = 1
c ih = number of particles extracted from holes
c joff = next hole location for extraction
c jss(1) = number of holes available to be filled
c jin = next hole location to be filled
c start loop
   10 mter = 0
      nps = 0
      kb(1) = js
      kb(2) = ks
c buffer outgoing particles
      jsl(1) = 0
      jsr(1) = 0
c load particle buffers
      do 40 j = 1, ih
      j1 = ihole(j+joff,n)
      xt = part(ic,j1)
c particles going down or backward
      if (xt.lt.edges(2*n-1)) then
         if (kb(n).eq.0) xt = xt + an
         if (jsl(1).lt.nbmax) then
            jsl(1) = jsl(1) + 1
            do 20 i = 1, idimp
            sbufl(i,jsl(1)) = part(i,j1)
   20       continue
            sbufl(ic,jsl(1)) = xt
         else
            nps = 1
            go to 50
         endif
c particles going up or forward
      else
         if (kb(n).eq.(nvp-1)) xt = xt - an
         if (jsr(1).lt.nbmax) then
            jsr(1) = jsr(1) + 1
            do 30 i = 1, idimp
            sbufr(i,jsr(1)) = part(i,j1)
   30       continue
            sbufr(ic,jsr(1)) = xt
         else
            nps = 1
            go to 50
         endif
      endif
   40 continue
   50 jss(1) = jsl(1) + jsr(1)
      joff = joff + jss(1)
      ih = ih - jss(1)
c check for full buffer condition
      ibflg(3) = nps
c copy particle buffers
   60 iter = iter + 2
      mter = mter + 1
c special case for one processor
      if (nvp.eq.1) then
         jsl(2) = jsr(1)
         do 64 j = 1, jsl(2)
         do 62 i = 1, idimp
         rbufl(i,j) = sbufr(i,j)
   62    continue
   64    continue
         jsr(2) = jsl(1)
         do 68 j = 1, jsr(2)
         do 66 i = 1, idimp
         rbufr(i,j) = sbufl(i,j)
   66    continue
   68    continue
c this segment is used for mpi computers
      else
c get particles from below and above or back and front
         kb(1) = js
         kb(2) = ks
         kl = kb(n)
         kb(n) = kl + 1
         if (kb(n).ge.nvp) kb(n) = kb(n) - nvp
         kr = kb(1) + nvpy*kb(2)
         kb(n) = kl - 1
         if (kb(n).lt.0) kb(n) = kb(n) + nvp
         kl = kb(1) + nvpy*kb(2)
c post receive
         itg = iter - 1
         call MPI_IRECV(rbufl,nbsize,mreal,kl,itg,lgrp,msid(1),ierr)
         call MPI_IRECV(rbufr,nbsize,mreal,kr,iter,lgrp,msid(2),ierr)
c send particles
         jsr(1) = idimp*jsr(1)
         call MPI_ISEND(sbufr,jsr(1),mreal,kr,itg,lgrp,msid(3),ierr)
         jsl(1) = idimp*jsl(1)
         call MPI_ISEND(sbufl,jsl(1),mreal,kl,iter,lgrp,msid(4),ierr)
c wait for particles to arrive
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsl(2) = nps/idimp
         call MPI_WAIT(msid(2),istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsr(2) = nps/idimp
      endif
c check if particles must be passed further
c check if any particles coming from above or front belong here
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do 70 j = 1, jsr(2)
      if (rbufr(ic,j).lt.edges(2*n-1)) jsl(1) = jsl(1) + 1
      if (rbufr(ic,j).ge.edges(2*n)) jsr(1) = jsr(1) + 1
   70 continue
c     if (jsr(1).ne.0) then
c        if (n.eq.1) then
c           write (2,*) kb+1,'Info:',jsr(1),' particles returning above'
c        else if (n.eq.2) then
c           write (2,*) kb+1,'Info:',jsr(1),' particles returning front'
c        endif
c     endif
c check if any particles coming from below or back belong here
      do 80 j = 1, jsl(2)
      if (rbufl(ic,j).ge.edges(2*n)) jsr(1) = jsr(1) + 1
      if (rbufl(ic,j).lt.edges(2*n-1)) jss(2) = jss(2) + 1
   80 continue
c     if (jss(2).ne.0) then
c        if (n.eq.1) then
c           write (2,*) kb+1,'Info:',jss(2),' particles returning below'
c        else if (n.eq.2) then
c           write (2,*) kb+1,'Info:',jss(2),' particles returning back'
c        endif
c     endif
      nps = jsl(1) + jsr(1) + jss(2)
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      if (nvp.ne.1) then
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
      endif
      if (nps.eq.0) go to 180
c remove particles which do not belong here
      kb(1) = js
      kb(2) = ks
c first check particles coming from above or front
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do 120 j = 1, jsr(2)
      xt = rbufr(ic,j)
c particles going down or back
      if (xt.lt.edges(2*n-1)) then
         jsl(1) = jsl(1) + 1
         if (kb(n).eq.0) xt = xt + an
         rbufr(ic,j) = xt
         do 90 i = 1, idimp
         sbufl(i,jsl(1)) = rbufr(i,j)
   90    continue
c particles going up or front, should not happen
      else if (xt.ge.edges(2*n)) then
         jsr(1) = jsr(1) + 1
         if (kb(n).eq.(nvp-1)) xt = xt - an
         rbufr(ic,j) = xt
         do 100 i = 1, idimp
         sbufr(i,jsr(1)) = rbufr(i,j)
  100    continue
c particles staying here
      else
         jss(2) = jss(2) + 1
         do 110 i = 1, idimp
         rbufr(i,jss(2)) = rbufr(i,j)
  110    continue
      endif
  120 continue
      jsr(2) = jss(2)
c next check particles coming from below or back
      jss(2) = 0
      do 160 j = 1, jsl(2)
      xt = rbufl(ic,j)
c particles going up or front
      if (xt.ge.edges(2*n)) then
         if (jsr(1).lt.nbmax) then
            jsr(1) = jsr(1) + 1
            if (kb(n).eq.(nvp-1)) xt = xt - an
            rbufl(ic,j) = xt
            do 130 i = 1, idimp
            sbufr(i,jsr(1)) = rbufl(i,j)
  130       continue
         else
            jss(2) = 2*npmax
            go to 170
         endif
c particles going down or back, should not happen
      else if (xt.lt.edges(2*n-1)) then
         if (jsl(1).lt.nbmax) then
            jsl(1) = jsl(1) + 1
            if (kb(n).eq.0) xt = xt + an
            rbufl(ic,j) = xt
            do 140 i = 1, idimp
            sbufl(i,jsl(1)) = rbufl(i,j)
  140       continue
         else
            jss(2) = 2*npmax
            go to 170
         endif
c particles staying here
      else
         jss(2) = jss(2) + 1
         do 150 i = 1, idimp
         rbufl(i,jss(2)) = rbufl(i,j)
  150    continue
      endif
  160 continue
  170 jsl(2) = jss(2)
c check if move would overflow particle array
  180 nps = npp + jsl(2) + jsr(2) - jss(1)
      ibflg(1) = nps
      ibflg(4) = -min0(npmax,nps)
      call PPIMAX(ibflg,iwork,5)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
c check for ihole overflow condition
      ierr = ibflg(5) - ntmax
      if (ierr.gt.0) then
         write (2,*) 'ihole overflow error, ierr = ', ierr
         info(1) = -ierr
         return
      endif
c distribute incoming particles from buffers
      nh = 0
c distribute particles coming from below or back into holes
      jss(2) = min0(jss(1),jsl(2))
      do 200 j = 1, jss(2)
      j1 = ihole(j+jin,n)
      do 190 i = 1, idimp
      part(i,j1) = rbufl(i,j)
  190 continue
c check if incoming particle is also out of bounds in z
      if (n.eq.1) then
         xt = part(iz,j1)
c if so, add it to list of particles in z
         if ((xt.lt.edges(2*n+1)).or.(xt.ge.edges(2*n+2))) then
            jh = jh + 1
            if (jh.le.ntmax) then
               ihole(jh+1,n+1) = j1
            else
               nh = 1
            endif
         endif
      endif
  200 continue
      jin = jin + jss(2)
      if (jss(1).gt.jsl(2)) then
         jss(2) = min0(jss(1)-jsl(2),jsr(2))
      else
         jss(2) = jsl(2) - jss(1)
      endif
      do 230 j = 1, jss(2)
c no more particles coming from below or back
c distribute particles coming from above or front into holes
      if (jss(1).gt.jsl(2)) then
         j1 = ihole(j+jin,n)
         do 210 i = 1, idimp
         part(i,j1) = rbufr(i,j)
  210    continue
c check if incoming particle is also out of bounds in z
         if (n.eq.1) then
            xt = part(iz,j1)
c if so, add it to list of particles in z
            if ((xt.lt.edges(2*n+1)).or.(xt.ge.edges(2*n+2))) then
               jh = jh + 1
               if (jh.le.ntmax) then
                  ihole(jh+1,n+1) = j1
               else
                  nh = 1
               endif
            endif
        endif
c no more holes
c distribute remaining particles from below or back into bottom
      else
         do 220 i = 1, idimp
         part(i,j+npp) = rbufl(i,j+jss(1))
  220    continue
c check if incoming particle is also out of bounds in z
         if (n.eq.1) then
            xt = part(iz,j+npp)
c if so, add it to list of particles in z
            if ((xt.lt.edges(2*n+1)).or.(xt.ge.edges(2*n+2))) then
               jh = jh + 1
               if (jh.le.ntmax) then
                  ihole(jh+1,n+1) = j + npp
               else
                  nh = 1
               endif
            endif
         endif
      endif
  230 continue
      if (jss(1).gt.jsl(2)) jin = jin + jss(2)
      nps = jsl(2) + jsr(2)
      if (jss(1).le.jsl(2)) then
         npp = npp + (jsl(2) - jss(1))
         jss(1) = jsl(2)
      endif
c no more holes
c distribute remaining particles from above into bottom
      jsr(2) = max0(0,nps-jss(1))
      jss(1) = jss(1) - jsl(2)
      do 250 j = 1, jsr(2)
      do 240 i = 1, idimp
      part(i,j+npp) = rbufr(i,j+jss(1))
  240 continue
c check if incoming particle is also out of bounds in z
      if (n.eq.1) then
         xt = part(iz,j+npp)
c if so, add it to list of particles in z
         if ((xt.lt.edges(2*n+1)).or.(xt.ge.edges(2*n+2))) then
            jh = jh + 1
            if (jh.le.ntmax) then
               ihole(jh+1,n+1) = j + npp
            else
               nh = 1
            endif
         endif
      endif
  250 continue
      npp = npp + jsr(2)
c check for ihole overflow condition
      if ((n.eq.1).and.(nh.gt.0)) ibflg(5) = jh
c holes left over
c fill up remaining holes in particle array with particles from bottom
      if (ih.eq.0) then
         jsr(2) = max0(0,ihole(1,n)-jin+1)
         nh = 0
c holes are stored in increasing value
         if (n.eq.1) then
            do 280 j = 1, jsr(2)
            j1 = npp - j + 1
            j2 = ihole(jsr(2)-j+jin+1,n)
            if (j1.gt.j2) then
c move particle only if it is below current hole
               do 260 i = 1, idimp
               part(i,j2) = part(i,j1)
  260          continue
c check if this move makes the ihole list for z invalid
               xt = part(iz,j1)
               if ((xt.lt.edges(2*n+1)).or.(xt.ge.edges(2*n+2))) then
                  i = jh + 1
c if so, adjust the list of holes
  270             j3 = ihole(i,n+1)
                  if (j3.ne.j1) then
                     i = i - 1
                     if (i.gt.1) go to 270
                     write (2,*) kstrt,'cannot find particle:n,j1=',n,j1
                     nh = 1
                  endif
c update ihole list to use new location
                  ihole(i,n+1) = j2
               endif
            endif
  280       continue
c holes may not be stored in increasing value
         else
            do 310 j = 1, jsr(2)
            j1 = npp - j + 1
            j2 = ihole(jsr(2)-j+jin+1,n)
            xt = part(iz,j1)
c determine if particle at location j1 represents an unfilled hole
            if ((xt.lt.edges(2*n-1)).or.(xt.ge.edges(2*n))) then
               i = jh + 2 - j
c if so, adjust the list of holes
  290          j3 = ihole(i,n)
               if (j3.ne.j1) then
                  i = i - 1
                  if (i.gt.1) go to 290
                  write (2,*) kstrt,'cannot find particle:n,j1=',n,j1
                  nh = 1
               endif
c update ihole list to use new location
               ihole(i,n) = j2
            else if (j1.gt.j2) then
c move particle only if it is below current hole
               do 300 i = 1, idimp
               part(i,j2) = part(i,j1)
  300          continue
            endif
  310       continue
c check for lost particle error
            if (nh.gt.0) call PPABORT
         endif
         jin = jin + jsr(2)
         npp = npp - jsr(2)
      endif
      jss(1) = 0
c check if any particles have to be passed further
      if (ibflg(3).gt.0) ibflg(3) = 1
      info(5+n) = max0(info(5+n),mter)
      if (ibflg(2).gt.0) then
c        write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (iter.lt.itermax) go to 60
         ierr = -((iter-2)/2)
         if (kstrt.eq.1) write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         return
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(3+n) = nter
         go to 10
      endif
      if (nter.gt.0) then
         if (kstrt.eq.1) then
            write (2,*) 'Info: ',nter,' buffer overflows, nbmax=', nbmax
         endif
      endif
c update ihole number in z
      if (n.eq.1) ihole(1,2) = jh
  320 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPMOVEG32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,
     1ny,nz,kstrt,nvpy,nvpz,idimp,npmax,idps,nbmax,ntmax,info)
c this subroutine moves particles into appropriate spatial regions
c ihole array is calculated from particles co-ordinates
c with periodic boundary conditions and 2D spatial decomposition
c output: part, ihole, npp, sbufr, sbufl, rbufr, rbufl, info
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = position z of particle n in partition
c part(4,n) = velocity vx of particle n in partition
c part(5,n) = velocity vy of particle n in partition
c part(6,n) = velocity vz of particle n in partition m
c edges(1:2) = lower/upper boundary in y of particle partition
c edges(3:4) = back/front boundary in z of particle partition
c npp = number of particles in partition
c sbufl = buffer for particles being sent to back processor
c sbufr = buffer for particles being sent to front processor
c rbufl = buffer for particles being received from back processor
c rbufr = buffer for particles being received from front processor
c ihole = location of holes left in particle arrays
c ny/nz = system length in y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition.
c idps = number of particle partition boundaries = 4
c nbmax =  size of buffers for passing particles between processors
c ntmax =  size of hole array for particles leaving processors
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4:5) = maximum number of buffer overflows in y/z
c info(6:7) = maximum number of particle passes required in y/z
      implicit none
      integer npp, ny, nz, kstrt, nvpy, nvpz, idimp, npmax, idps, nbmax
      integer ntmax
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer ihole, info
      dimension part(idimp,npmax)
      dimension edges(idps)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension rbufl(idimp,nbmax), rbufr(idimp,nbmax)
      dimension ihole(ntmax+1)
      dimension info(7)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
c iy/iz = partitioned co-ordinates
      integer iy, iz
      parameter(iy=2,iz=3)
      integer i, j, n, js, ks, ic, nvp, iter, nps, itg, kl, kr, j1, j2
      integer ih, joff, jin, nbsize, nter, mter, itermax, ierr
      integer msid, istatus
      integer kb, jsl, jsr, jss, ibflg, iwork
      real an, xt
      dimension msid(4), istatus(lstat)
      dimension kb(2), jsl(2), jsr(2), jss(2), ibflg(4), iwork(4)
c js/ks = processor co-ordinates in y/z=> idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      nbsize = idimp*nbmax
      info(1) = 0
      info(6) = 0
      info(7) = 0
      itermax = 2000
c buffer outgoing particles, first in y then in z direction
      do 280 n = 1, 2
      if (n.eq.1) then
         ic = iy
         nvp = nvpy
         an = real(ny)
      else if (n.eq.2) then
         ic = iz
         nvp = nvpz
         an = real(nz)
      endif
      iter = 2
      nter = 0
      joff = 1
c ih = number of particles extracted from holes
c joff = next hole location for extraction
c jss(1) = number of holes available to be filled
c jin = next hole location to be filled
c start loop
   10 mter = 0
      nps = 0
      jin = 1
      kb(1) = js
      kb(2) = ks
c buffer outgoing particles
      jsl(1) = 0
      jsr(1) = 0
c load particle buffers
      do 40 j = 1, npp
      xt = part(ic,j)
c particles going down or backward
      if (xt.lt.edges(2*n-1)) then
         if (kb(n).eq.0) xt = xt + an
         if (jsl(1).lt.nbmax) then
            jsl(1) = jsl(1) + 1
            do 20 i = 1, idimp
            sbufl(i,jsl(1)) = part(i,j)
   20       continue
            sbufl(ic,jsl(1)) = xt
            ihole(jsl(1)+jsr(1)+1) = j
         else
            nps = 1
            go to 50
         endif
c particles going up or forward
      else if (xt.ge.edges(2*n)) then
         if (kb(n).eq.(nvp-1)) xt = xt - an
         if (jsr(1).lt.nbmax) then
            jsr(1) = jsr(1) + 1
            do 30 i = 1, idimp
            sbufr(i,jsr(1)) = part(i,j)
   30       continue
            sbufr(ic,jsr(1)) = xt
            ihole(jsl(1)+jsr(1)+1) = j
         else
            nps = 1
            go to 50
         endif
      endif
   40 continue
   50 jss(1) = jsl(1) + jsr(1)
      joff = joff + jss(1)
      ihole(1) = jss(1)
      ih = 0
c check for full buffer condition
      ibflg(3) = nps
c copy particle buffers
   60 iter = iter + 2
      mter = mter + 1
c special case for one processor
      if (nvp.eq.1) then
         jsl(2) = jsr(1)
         do 64 j = 1, jsl(2)
         do 62 i = 1, idimp
         rbufl(i,j) = sbufr(i,j)
   62    continue
   64    continue
         jsr(2) = jsl(1)
         do 68 j = 1, jsr(2)
         do 66 i = 1, idimp
         rbufr(i,j) = sbufl(i,j)
   66    continue
   68    continue
c this segment is used for mpi computers
      else
c get particles from below and above or back and front
         kb(1) = js
         kb(2) = ks
         kl = kb(n)
         kb(n) = kl + 1
         if (kb(n).ge.nvp) kb(n) = kb(n) - nvp
         kr = kb(1) + nvpy*kb(2)
         kb(n) = kl - 1
         if (kb(n).lt.0) kb(n) = kb(n) + nvp
         kl = kb(1) + nvpy*kb(2)
c post receive
         itg = iter - 1
         call MPI_IRECV(rbufl,nbsize,mreal,kl,itg,lgrp,msid(1),ierr)
         call MPI_IRECV(rbufr,nbsize,mreal,kr,iter,lgrp,msid(2),ierr)
c send particles
         jsr(1) = idimp*jsr(1)
         call MPI_ISEND(sbufr,jsr(1),mreal,kr,itg,lgrp,msid(3),ierr)
         jsl(1) = idimp*jsl(1)
         call MPI_ISEND(sbufl,jsl(1),mreal,kl,iter,lgrp,msid(4),ierr)
c wait for particles to arrive
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsl(2) = nps/idimp
         call MPI_WAIT(msid(2),istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsr(2) = nps/idimp
      endif
c check if particles must be passed further
c check if any particles coming from above or front belong here
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do 70 j = 1, jsr(2)
      if (rbufr(ic,j).lt.edges(2*n-1)) jsl(1) = jsl(1) + 1
      if (rbufr(ic,j).ge.edges(2*n)) jsr(1) = jsr(1) + 1
   70 continue
c     if (jsr(1).ne.0) then
c        if (n.eq.1) then
c           write (2,*) kb+1,'Info:',jsr(1),' particles returning above'
c        else if (n.eq.2) then
c           write (2,*) kb+1,'Info:',jsr(1),' particles returning front'
c        endif
c     endif
c check if any particles coming from below or back belong here
      do 80 j = 1, jsl(2)
      if (rbufl(ic,j).ge.edges(2*n)) jsr(1) = jsr(1) + 1
      if (rbufl(ic,j).lt.edges(2*n-1)) jss(2) = jss(2) + 1
   80 continue
c     if (jss(2).ne.0) then
c        if (n.eq.1) then
c           write (2,*) kb+1,'Info:',jss(2),' particles returning below'
c        else if (n.eq.2) then
c           write (2,*) kb+1,'Info:',jss(2),' particles returning back'
c        endif
c     endif
      nps = jsl(1) + jsr(1) + jss(2)
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      if (nvp.ne.1) then
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
      endif
      if (nps.eq.0) go to 180
c remove particles which do not belong here
      kb(1) = js
      kb(2) = ks
c first check particles coming from above or front
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do 120 j = 1, jsr(2)
      xt = rbufr(ic,j)
c particles going down or back
      if (xt.lt.edges(2*n-1)) then
         jsl(1) = jsl(1) + 1
         if (kb(n).eq.0) xt = xt + an
         rbufr(ic,j) = xt
         do 90 i = 1, idimp
         sbufl(i,jsl(1)) = rbufr(i,j)
   90    continue
c particles going up or front, should not happen
      else if (xt.ge.edges(2*n)) then
         jsr(1) = jsr(1) + 1
         if (kb(n).eq.(nvp-1)) xt = xt - an
         rbufr(ic,j) = xt
         do 100 i = 1, idimp
         sbufr(i,jsr(1)) = rbufr(i,j)
  100    continue
c particles staying here
      else
         jss(2) = jss(2) + 1
         do 110 i = 1, idimp
         rbufr(i,jss(2)) = rbufr(i,j)
  110    continue
      endif
  120 continue
      jsr(2) = jss(2)
c next check particles coming from below or back
      jss(2) = 0
      do 160 j = 1, jsl(2)
      xt = rbufl(ic,j)
c particles going up or front
      if (xt.ge.edges(2*n)) then
         if (jsr(1).lt.nbmax) then
            jsr(1) = jsr(1) + 1
            if (kb(n).eq.(nvp-1)) xt = xt - an
            rbufl(ic,j) = xt
            do 130 i = 1, idimp
            sbufr(i,jsr(1)) = rbufl(i,j)
  130       continue
         else
            jss(2) = 2*npmax
            go to 170
         endif
c particles going down or back, should not happen
      else if (xt.lt.edges(2*n-1)) then
         if (jsl(1).lt.nbmax) then
            jsl(1) = jsl(1) + 1
            if (kb(n).eq.0) xt = xt + an
            rbufl(ic,j) = xt
            do 140 i = 1, idimp
            sbufl(i,jsl(1)) = rbufl(i,j)
  140       continue
         else
            jss(2) = 2*npmax
            go to 170
         endif
c particles staying here
      else
         jss(2) = jss(2) + 1
         do 150 i = 1, idimp
         rbufl(i,jss(2)) = rbufl(i,j)
  150    continue
      endif
  160 continue
  170 jsl(2) = jss(2)
c check if move would overflow particle array
  180 nps = npp + jsl(2) + jsr(2) - jss(1)
      ibflg(1) = nps
      ibflg(4) = -min0(npmax,nps)
      call PPIMAX(ibflg,iwork,4)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
c distribute incoming particles from buffers
c distribute particles coming from below or back into holes
      jss(2) = min0(jss(1),jsl(2))
      do 200 j = 1, jss(2)
      j1 = ihole(j+jin)
      do 190 i = 1, idimp
      part(i,j1) = rbufl(i,j)
  190 continue
  200 continue
      jin = jin + jss(2)
      if (jss(1).gt.jsl(2)) then
         jss(2) = min0(jss(1)-jsl(2),jsr(2))
      else
         jss(2) = jsl(2) - jss(1)
      endif
      do 230 j = 1, jss(2)
c no more particles coming from below or back
c distribute particles coming from above or front into holes
      if (jss(1).gt.jsl(2)) then
         j1 = ihole(j+jin)
         do 210 i = 1, idimp
         part(i,j1) = rbufr(i,j)
  210    continue
c no more holes
c distribute remaining particles from below or back into bottom
      else
         do 220 i = 1, idimp
         part(i,j+npp) = rbufl(i,j+jss(1))
  220    continue
      endif
  230 continue
      if (jss(1).gt.jsl(2)) jin = jin + jss(2)
      nps = jsl(2) + jsr(2)
      if (jss(1).le.jsl(2)) then
         npp = npp + (jsl(2) - jss(1))
         jss(1) = jsl(2)
      endif
c no more holes
c distribute remaining particles from above into bottom
      jsr(2) = max0(0,nps-jss(1))
      jss(1) = jss(1) - jsl(2)
      do 250 j = 1, jsr(2)
      do 240 i = 1, idimp
      part(i,j+npp) = rbufr(i,j+jss(1))
  240 continue
  250 continue
      npp = npp + jsr(2)
c holes left over
c fill up remaining holes in particle array with particles from bottom
      if (ih.eq.0) then
         jsr(2) = max0(0,ihole(1)-jin+1)
         do 270 j = 1, jsr(2)
         j1 = npp - j + 1
         j2 = ihole(jsr(2)-j+jin+1)
         if (j1.gt.j2) then
c move particle only if it is below current hole
            do 260 i = 1, idimp
            part(i,j2) = part(i,j1)
  260       continue
         endif
  270    continue
         jin = jin + jsr(2)
         npp = npp - jsr(2)
      endif
      jss(1) = 0
c check if any particles have to be passed further
      if (ibflg(3).gt.0) ibflg(3) = 1
      info(5+n) = max0(info(5+n),mter)
      if (ibflg(2).gt.0) then
c        write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (iter.lt.itermax) go to 60
         ierr = -((iter-2)/2)
         if (kstrt.eq.1) write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         return
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(3+n) = nter
         go to 10
      endif
      if (nter.gt.0) then
         if (kstrt.eq.1) then
            write (2,*) 'Info: ',nter,' buffer overflows, nbmax=', nbmax
         endif
      endif
  280 continue
      return
      end
