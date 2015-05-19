c-----------------------------------------------------------------------
c Basic parallel PIC library for MPI communications
c pplib2.f contains basic communications procedures for 1d partitions:
c PPINIT2 initializes parallel processing for Fortran90, returns
c         number of processors and processor id.
c PPEXIT terminates parallel processing.
c PPABORT aborts parallel processing.
c PWTIMERA performs parallel local wall clock timing.
c PPSUM performs parallel sum of a real vector.
c PPDSUM performs parallel sum of a double precision vector.
c PPIMAX performs parallel maximum of an integer vector.
c PPDMAX performs parallel maximum of a double precision vector.
c PPNCGUARD2L copies data to guard cells in y for scalar data, linear
c             interpolation, and distributed data with non-uniform
c             partition.
c PPNAGUARD2L adds guard cells in y for scalar array, linear
c             interpolation, and distributed data with non-uniform
c             partition.
c PPNACGUARD2L adds guard cells in y for vector array, linear
c              interpolation, and distributed data with non-uniform
c              partition.
c PPTPOSE performs a transpose of a complex scalar array, distributed
c         in y, to a complex scalar array, distributed in x.
c PPNTPOSE performs a transpose of an n component complex vector array,
c          distributed in y, to an n component complex vector array,
c          distributed in x.
c PPMOVE2 moves particles into appropriate spatial regions with periodic
c         boundary conditions.  Assumes ihole list has been found.
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: april 19, 2015
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
      subroutine PPNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx)
c this subroutine copies data to guard cells in non-uniform partitions
c f(j,k) = real data for grid j,k in particle partition.
c the grid is non-uniform and includes one extra guard cell.
c output: f
c nyp = number of primary gridpoints in field partition
c it is assumed the nyp > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cell.
c linear interpolation, for distributed data
      implicit none
      integer nyp, kstrt, nvp, nxv, nypmx
      real f
      dimension f(nxv,nypmx)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer j, ks, moff, kl, kr
      integer istatus, msid, ierr
      dimension istatus(lstat)
c special case for one processor
      if (nvp.eq.1) then
         do 10 j = 1, nxv
         f(j,nyp+1) = f(j,1)
   10    continue
         return
      endif
      ks = kstrt - 1
      moff = nypmx*nvp + 2
c copy to guard cells
      kr = ks + 1
      if (kr.ge.nvp) kr = kr - nvp
      kl = ks - 1
      if (kl.lt.0)  kl = kl + nvp
      ks = nyp + 1
c this segment is used for mpi computers
      call MPI_IRECV(f(1,ks),nxv,mreal,kr,moff,lgrp,msid,ierr)
      call MPI_SEND(f,nxv,mreal,kl,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      return
      end
c-----------------------------------------------------------------------
      subroutine PPNAGUARD2L(f,scr,nyp,nx,kstrt,nvp,nxv,nypmx)
c this subroutine adds data from guard cells in non-uniform partitions
c f(j,k) = real data for grid j,k in particle partition.
c the grid is non-uniform and includes one extra guard cell.
c output: f, scr
c scr(j) = scratch array for particle partition
c nyp = number of primary gridpoints in particle partition
c it is assumed the nyp > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c linear interpolation, for distributed data
      implicit none
      integer nyp, kstrt, nvp, nx, nxv, nypmx
      real f, scr
      dimension f(nxv,nypmx), scr(nxv)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer j, nx1, ks, moff, kl, kr
      integer istatus, msid, ierr
      dimension istatus(lstat)
      nx1 = nx + 1
c special case for one processor
      if (nvp.eq.1) then
         do 10 j = 1, nx1
         f(j,1) = f(j,1) + f(j,nyp+1)
         f(j,nyp+1) = 0.
   10    continue
         return
      endif
      ks = kstrt - 1
      moff = nypmx*nvp + 1
c add guard cells
      kr = ks + 1
      if (kr.ge.nvp) kr = kr - nvp
      kl = ks - 1
      if (kl.lt.0) kl = kl + nvp
      ks = nyp + 1
c this segment is used for mpi computers
      call MPI_IRECV(scr,nxv,mreal,kl,moff,lgrp,msid,ierr)
      call MPI_SEND(f(1,ks),nxv,mreal,kr,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 20 j = 1, nx1
      f(j,1) = f(j,1) + scr(j)
      f(j,nyp+1) = 0.0
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPNACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx)
c this subroutine adds data from guard cells in non-uniform partitions
c f(ndim,j,k) = real data for grid j,k in particle partition.
c the grid is non-uniform and includes one extra guard cell.
c output: f, scr
c scr(ndim,j) = scratch array for particle partition
c nyp = number of primary gridpoints in particle partition
c it is assumed the nyp > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c ndim = leading dimension of array f
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c linear interpolation, for distributed data
      implicit none
      integer nyp, kstrt, nvp, nx, ndim, nxv, nypmx
      real f, scr
      dimension f(ndim,nxv,nypmx), scr(ndim,nxv)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer j, n, nx1, ks, moff, kl, kr
      integer nnxv
      integer istatus, msid, ierr
      dimension istatus(lstat)
      nx1 = nx + 1
c special case for one processor
      if (nvp.eq.1) then
         do 20 j = 1, nx1
         do 10 n = 1, ndim
         f(n,j,1) = f(n,j,1) + f(n,j,nyp+1)
         f(n,j,nyp+1) = 0.0
   10    continue
   20    continue
         return
      endif
      ks = kstrt - 1
      moff = nypmx*nvp + 1
      nnxv = ndim*nxv
c add guard cells
      kr = ks + 1
      if (kr.ge.nvp) kr = kr - nvp
      kl = ks - 1
      if (kl.lt.0) kl = kl + nvp
      ks = nyp + 1
c this segment is used for mpi computers
      call MPI_IRECV(scr,nnxv,mreal,kl,moff,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,ks),nnxv,mreal,kr,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 40 j = 1, nx1
      do 30 n = 1, ndim
      f(n,j,1) = f(n,j,1) + scr(n,j)
      f(n,j,nyp+1) = 0.0
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd,  
     1kypd)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(k+kyp*(m-1),j,l) = f(j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny = number of points in x/y
c kxp/kyp = number of data values per block in x/y
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nxv/nyv = first dimension of f/g
c kypd/kxpd = second dimension of f/g
      implicit none
      integer nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv, kxpd, kypd
      complex f, g, s, t
      dimension f(nxv,kypd), g(nyv,kxpd)
      dimension s(kxp*kyp), t(kxp*kyp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld
      integer ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 1
      kxps = min(kxp,max(0,nx-kxp*ks))
      kyps = min(kyp,max(0,ny-kyp*ks))
      kxyp = kxp*kyp
c special case for one processor
      if (nvp.eq.1) then
         do 20 k = 1, kyp
         do 10 j = 1, kxp
         g(k,j) = f(j,k)
   10    continue
   20    continue
         return
      endif
c this segment is used for shared memory computers
c     do 60 m = 1, min(ny,nvp)
c     koff = kyp*(m - 1)
c     do 50 k = 1, min(kyp,max(0,ny-koff))
c     do 40 l = 1, min(nx,nvp)
c     joff = kxp*(l - 1)
c     do 30 j = 1, min(kxp,max(0,nx-joff))
c     g(k+koff,j+joff) = f(j+joff,k+koff)
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      do 70 n = 1, nvp
      id = n - ks - 1
      if (id.lt.0) id = id + nvp
c extract data to send
      joff = kxp*id
      ld = min(kxp,max(0,nx-joff))
      do 40 k = 1, kyps
      do 30 j = 1, ld
      s(j+ld*(k-1)) = f(j+joff,k)
   30 continue
   40 continue
      ld = ld*kyps
c post receive
      call MPI_IRECV(t,kxyp,mcplx,id,n,lgrp,msid,ierr)
c send data
      call MPI_SEND(s,ld,mcplx,id,n,lgrp,ierr)
c receive data
      call MPI_WAIT(msid,istatus,ierr)
c insert data received
      koff = kyp*id
      ld = min(kyp,max(0,ny-koff))
      do 60 k = 1, ld
      do 50 j = 1, kxps
      g(k+koff,j) = t(j+kxps*(k-1))
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPNTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv, 
     1kxpd,kypd)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny = number of points in x/y
c kxp/kyp = number of data values per block in x/y
c kstrt = starting data block number
c nvp = number of real or virtual processors
c ndim = leading dimension of arrays f and g
c nxv/nyv = first dimension of f/g
c kypd/kxpd = second dimension of f/g
      implicit none
      integer nx, ny, kxp, kyp, kstrt, nvp, ndim, nxv, nyv, kxpd, kypd
      complex f, g, s, t
      dimension f(ndim,nxv,kypd), g(ndim,nyv,kxpd)
      dimension s(ndim,kxp*kyp), t(ndim,kxp*kyp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer i, n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld
      integer ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 1
      kxps = min(kxp,max(0,nx-kxp*ks))
      kyps = min(kyp,max(0,ny-kyp*ks))
      kxyp = ndim*kxp*kyp
c special case for one processor
      if (nvp.eq.1) then
         do 30 k = 1, kyp
         do 20 j = 1, kxp
         do 10 i = 1, ndim
         g(i,k,j) = f(i,j,k)
   10    continue
   20    continue
   30    continue
         return
      endif
c this segment is used for shared memory computers
c     do 80 m = 1, min(ny,nvp)
c     koff = kyp*(m - 1)
c     do 70 k = 1, min(kyp,max(0,ny-koff))
c     do 60 l = 1, min(nx,nvp)
c     joff = kxp*(l - 1)
c     do 50 j = 1, min(kxp,max(0,nx-joff))
c     do 40 i = 1, ndim
c     g(i,k+koff,j+joff) = f(i,j+joff,k+koff)
c  40 continue
c  50 continue
c  60 continue
c  70 continue
c  80 continue
c this segment is used for mpi computers
      do 100 n = 1, nvp
      id = n - ks - 1
      if (id.lt.0) id = id + nvp
c extract data to send
      joff = kxp*id
      ld = min(kxp,max(0,nx-joff))
      do 60 k = 1, kyps
      do 50 j = 1, ld
      do 40 i = 1, ndim
      s(i,j+ld*(k-1)) = f(i,j+joff,k)
   40 continue
   50 continue
   60 continue
      ld = ndim*ld*kyps
c post receive
      call MPI_IRECV(t,kxyp,mcplx,id,n,lgrp,msid,ierr)
c send data
      call MPI_SEND(s,ld,mcplx,id,n,lgrp,ierr)
c receive data
      call MPI_WAIT(msid,istatus,ierr)
c insert data received
      koff = kyp*id
      ld = min(kyp,max(0,ny-koff))
      do 90 k = 1, ld
      do 80 j = 1, kxps
      do 70 i = 1, ndim
      g(i,k+koff,j) = t(i,j+kxps*(k-1))
   70 continue
   80 continue
   90 continue
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,ny
     1,kstrt,nvp,idimp,npmax,idps,nbmax,ntmax,info)
c this subroutine moves particles into appropriate spatial regions
c periodic boundary conditions
c output: part, npp, sbufr, sbufl, rbufr, rbufl, info
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = velocity vx of particle n in partition
c part(4,n) = velocity vy of particle n in partition
c edges(1:2) = lower:lower boundary of particle partition
c npp = number of particles in partition
c sbufl = buffer for particles being sent to lower processor
c sbufr = buffer for particles being sent to upper processor
c rbufl = buffer for particles being received from lower processor
c rbufr = buffer for particles being received from upper processor
c ihole = location of holes left in particle arrays
c ny = system length in y direction
c kstrt = starting data block number
c nvp = number of real or virtual processors
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition.
c idps = number of partition boundaries
c nbmax =  size of buffers for passing particles between processors
c ntmax =  size of hole array for particles leaving processors
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4) = maximum number of buffer overflows
c info(5) = maximum number of particle passes required
      implicit none
      integer ny, kstrt, nvp, idimp, npmax, idps, nbmax, ntmax
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, info
      dimension part(idimp,npmax)
      dimension edges(idps)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension rbufl(idimp,nbmax), rbufr(idimp,nbmax)
      dimension ihole(ntmax+1)
      dimension info(5)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
c iy = partitioned co-ordinate
      integer iy
      parameter(iy=2)
      integer ierr, ks, ih, iter, nps, itg, kl, kr, j, j1, j2, i
      integer joff, jin, nbsize, nter, mter, itermax
      integer msid, istatus
      integer jsl, jsr, jss, ibflg, iwork
      real any, yt
      dimension msid(4), istatus(lstat)
      dimension jsl(2), jsr(2), jss(2), ibflg(4), iwork(4)
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
c ih = number of particles extracted from holes
c joff = next hole location for extraction
c jss(1) = number of holes available to be filled
c jin = next hole location to be filled
c start loop
   10 mter = 0
      nps = 0
c buffer outgoing particles
      jsl(1) = 0
      jsr(1) = 0
c load particle buffers
      do 40 j = 1, ih
      j1 = ihole(j+joff)
      yt = part(iy,j1)
c particles going down
      if (yt.lt.edges(1)) then
         if (ks.eq.0) yt = yt + any
         if (jsl(1).lt.nbmax) then
            jsl(1) = jsl(1) + 1
            do 20 i = 1, idimp
            sbufl(i,jsl(1)) = part(i,j1)
   20       continue
            sbufl(iy,jsl(1)) = yt
         else
            nps = 1
            go to 50
         endif
c particles going up
      else
         if (ks.eq.(nvp-1)) yt = yt - any
         if (jsr(1).lt.nbmax) then
            jsr(1) = jsr(1) + 1
            do 30 i = 1, idimp
            sbufr(i,jsr(1)) = part(i,j1)
   30       continue
            sbufr(iy,jsr(1)) = yt
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
c get particles from below and above
         kr = ks + 1
         if (kr.ge.nvp) kr = kr - nvp
         kl = ks - 1
         if (kl.lt.0) kl = kl + nvp
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
c check if any particles coming from above belong here
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do 70 j = 1, jsr(2)
      if (rbufr(iy,j).lt.edges(1)) jsl(1) = jsl(1) + 1
      if (rbufr(iy,j).ge.edges(2)) jsr(1) = jsr(1) + 1
   70 continue
c     if (jsr(1).ne.0) write (2,*) ks+1, 'Info: particles returning up'
c check if any particles coming from below belong here
      do 80 j = 1, jsl(2)
      if (rbufl(iy,j).ge.edges(2)) jsr(1) = jsr(1) + 1
      if (rbufl(iy,j).lt.edges(1)) jss(2) = jss(2) + 1
   80 continue
c     if (jss(2).ne.0) write (2,*) ks+1,'Info: particles returning down'
      nps = jsl(1) + jsr(1) + jss(2)
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      if (nvp.ne.1) then
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
      endif
      if (nps.eq.0) go to 180
c remove particles which do not belong here
c first check particles coming from above
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do 120 j = 1, jsr(2)
      yt = rbufr(iy,j)
c particles going down
      if (yt.lt.edges(1)) then
         jsl(1) = jsl(1) + 1
         if (ks.eq.0) yt = yt + any
         rbufr(iy,j) = yt
         do 90 i = 1, idimp
         sbufl(i,jsl(1)) = rbufr(i,j)
   90    continue
c particles going up, should not happen
      else if (yt.ge.edges(2)) then
         jsr(1) = jsr(1) + 1
         if (ks.eq.(nvp-1)) yt = yt - any
         rbufr(iy,j) = yt
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
c next check particles coming from below
      jss(2) = 0
      do 160 j = 1, jsl(2)
      yt = rbufl(iy,j)
c particles going up
      if (yt.ge.edges(2)) then
         if (jsr(1).lt.nbmax) then
            jsr(1) = jsr(1) + 1
            if (ks.eq.(nvp-1)) yt = yt - any
            rbufl(iy,j) = yt
            do 130 i = 1, idimp
            sbufr(i,jsr(1)) = rbufl(i,j)
  130       continue
         else
            jss(2) = 2*npmax
            go to 170
         endif
c particles going down, should not happen
      else if (yt.lt.edges(1)) then
         if (jsl(1).lt.nbmax) then
            jsl(1) = jsl(1) + 1
            if (ks.eq.0) yt = yt + any
            rbufl(iy,j) = yt
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
c distribute particles coming from below into holes
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
c no more particles coming from below
c distribute particles coming from above into holes
      if (jss(1).gt.jsl(2)) then
         j1 = ihole(j+jin)
         do 210 i = 1, idimp
         part(i,j1) = rbufr(i,j)
  210    continue
c no more holes
c distribute remaining particles from below into bottom
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
      info(5) = max0(info(5),mter)
      if (ibflg(2).gt.0) then
c        write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (iter.lt.itermax) go to 60
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         return
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(4) = nter
         go to 10
      endif
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end

