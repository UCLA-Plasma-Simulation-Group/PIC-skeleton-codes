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
c PPNCGUARD2L copies data to guard cells in y for scalar data, linear
c             interpolation, and distributed data with non-uniform
c             partition.
c PPNAGUARD2L adds guard cells in y for scalar array, linear
c             interpolation, and distributed data with non-uniform
c             partition.
c PPTPOSE performs a transpose of a complex scalar array, distributed
c         in y, to a complex scalar array, distributed in x.
c PPNTPOSE performs a transpose of an n component complex vector array,
c          distributed in y, to an n component complex vector array,
c          distributed in x.
c PPPMOVE2 moves particles into appropriate spatial regions for tiled
c          distributed data.
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: April 10, 2014
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
!$OMP PARALLEL DO PRIVATE(j,k)
         do 20 k = 1, kyp
         do 10 j = 1, kxp
         g(k,j) = f(j,k)
   10    continue
   20    continue
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(j,k)
      do 40 k = 1, kyps
      do 30 j = 1, ld
      s(j+ld*(k-1)) = f(j+joff,k)
   30 continue
   40 continue
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(j,k)
      do 60 k = 1, ld
      do 50 j = 1, kxps
      g(k+koff,j) = t(j+kxps*(k-1))
   50 continue
   60 continue
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do 30 k = 1, kyp
         do 20 j = 1, kxp
         do 10 i = 1, ndim
         g(i,k,j) = f(i,j,k)
   10    continue
   20    continue
   30    continue
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do 60 k = 1, kyps
      do 50 j = 1, ld
      do 40 i = 1, ndim
      s(i,j+ld*(k-1)) = f(i,j+joff,k)
   40 continue
   50 continue
   60 continue
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do 90 k = 1, ld
      do 80 j = 1, kxps
      do 70 i = 1, ndim
      g(i,k+koff,j) = t(i,j+kxps*(k-1))
   70 continue
   80 continue
   90 continue
!$OMP END PARALLEL DO
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,  
     1kstrt,nvp,idimp,nbmax,mx1)
c this subroutine moves particles into appropriate spatial regions
c for distributed data, with 1d domain decomposition in y.
c tiles are assumed to be arranged in 2D linear memory
c output: rbufr, rbufl, mcll, mclr
c sbufl = buffer for particles being sent to lower processor
c sbufr = buffer for particles being sent to upper processor
c rbufl = buffer for particles being received from lower processor
c rbufr = buffer for particles being received from upper processor
c ncll = particle number being sent to lower processor
c nclr = particle number being sent to upper processor
c mcll = particle number being received from lower processor
c mclr = particle number being received from upper processor
c kstrt = starting data block number
c nvp = number of real or virtual processors
c idimp = size of phase space = 4 or 5
c nbmax =  size of buffers for passing particles between processors
c mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer kstrt, nvp, idimp, nbmax, mx1
      real sbufr, sbufl, rbufr, rbufl
      integer ncll, nclr, mcll, mclr
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension rbufl(idimp,nbmax), rbufr(idimp,nbmax)
      dimension ncll(3,mx1), nclr(3,mx1), mcll(3,mx1), mclr(3,mx1)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, ks, kl, kr, i, j, jsl, jsr
      integer nbsize, ncsize
      integer msid, itg, istatus
      dimension msid(8), itg(4), istatus(lstat)
      data itg /3,4,5,6/
      ks = kstrt - 1
      nbsize = idimp*nbmax
      ncsize = 3*mx1
c copy particle buffers: update rbufl, rbufr, mcll, mclr
c special case for one processor
      if (nvp.eq.1) then
         do 20 j = 1, mx1
         do 10 i = 1, 3
         mcll(i,j) = nclr(i,j)
   10    continue
   20    continue
         do 40 j = 1, mx1
         do 30 i = 1, 3
         mclr(i,j) = ncll(i,j)
   30    continue
   40    continue
         do 60 j = 1, nclr(3,mx1)
         do 50 i = 1, idimp
         rbufl(i,j) = sbufr(i,j)
   50    continue
   60    continue
         do 80 j = 1, ncll(3,mx1)
         do 70 i = 1, idimp
         rbufr(i,j) = sbufl(i,j)
   70    continue
   80    continue
c this segment is used for mpi computers
      else
c get particles from below and above
         kr = ks + 1
         if (kr.ge.nvp) kr = kr - nvp
         kl = ks - 1
         if (kl.lt.0) kl = kl + nvp
c post receives
         call MPI_IRECV(mcll,ncsize,mint,kl,itg(1),lgrp,msid(1),ierr)
         call MPI_IRECV(mclr,ncsize,mint,kr,itg(2),lgrp,msid(2),ierr)
         call MPI_IRECV(rbufl,nbsize,mreal,kl,itg(3),lgrp,msid(3),ierr)
         call MPI_IRECV(rbufr,nbsize,mreal,kr,itg(4),lgrp,msid(4),ierr)
c send particle number offsets
         call MPI_ISEND(nclr,ncsize,mint,kr,itg(1),lgrp,msid(5),ierr)
         call MPI_ISEND(ncll,ncsize,mint,kl,itg(2),lgrp,msid(6),ierr)
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_WAIT(msid(2),istatus,ierr)
c send particles
         jsr = idimp*nclr(3,mx1)
         call MPI_ISEND(sbufr,jsr,mreal,kr,itg(3),lgrp,msid(7),ierr)
         jsl = idimp*ncll(3,mx1)
         call MPI_ISEND(sbufl,jsl,mreal,kl,itg(4),lgrp,msid(8),ierr)
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
      endif
c make sure sbufr, sbufl, ncll, and nclr have been sent
      if (nvp.ne.1) then
         do 90 i = 1, 4
         call MPI_WAIT(msid(i+4),istatus,ierr)
   90    continue
      endif
      return
      end
