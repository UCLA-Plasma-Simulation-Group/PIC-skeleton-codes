c vector add test program for MPI
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine MPADD(a,b,c,nx)
      integer nx
      real a, b, c
      dimension a(nx), b(nx), c(nx)
      integer j
      do 10 j = 1, nx
      a(j) = b(j) + c(j)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine INIT_MPI(idproc,nvp,irc)
c this subroutine initializes parallel processing
c lgrp communicator = MPI_COMM_WORLD
c output: idproc, nvp
c idproc = processor id in lgrp communicator
c nvp = number of real or virtual processors obtained
c error code is modified only if there is an error
      implicit none
      integer idproc, nvp, irc
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
      integer ndprec
c     integer idprec
      integer ierror
      logical flag
      save /PPARMS/, /PPARMSX/
c ndprec = (0,1) = (no,yes) use (normal,autodouble) precision
      ndprec = 0
c idprec = (0,1) = (no,yes) use (normal,autodouble) integer precision
c     idprec = 0
c this segment is used for mpi computers
      if (MPI_STATUS_SIZE.gt.lstat) then
         write (2,*) ' status size too small, actual/required = ', lstat
     1, MPI_STATUS_SIZE
         irc = 1
         return
      endif
c indicate whether MPI_INIT has been called
      call MPI_INITIALIZED(flag,ierror)
      if (.not.flag) then
c initialize the MPI execution environment
         call MPI_INIT(ierror)
         if (ierror.ne.0) then
            irc = 1
            return
         endif
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
      subroutine PPMAX(f,g,nxp)
c this subroutine finds parallel maximum for each element of a vector
c that is, f(j,k) = maximum as a function of k of f(j,k)
c at the end, all processors contain the same maximum.
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
c mreal = default real type
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c mmax = MPI_MAX
      common /PPARMSX/ msum, mmax
c local data
      integer j, ierr
c find maximum
      call MPI_ALLREDUCE(f,g,nxp,mreal,mmax,lgrp,ierr)
c copy output from scratch array
      do 10 j = 1, nxp
      f(j) = g(j)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine VSCATTER(f,g,nx,nxp)
c this subroutine performs a scatter of a real vector, that is,
c successive nxp blocks of source f on processor 0, are sent to
c destination g on successive processors:
c g(1:nxp) = f(1+nxp*idproc:min(nxp*(idproc+1),nx))
c f = input real data
c g = output real data
c nx = size of source data
c nxp = size of destination data
      implicit none
      integer nx, nxp
      real f, g
      dimension f(nx), g(nxp)
c common blocks for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer i, j, idproc, nvp, is, nxps, koff, ierr
      integer istatus
      dimension istatus(lstat)
c find processor id
      call MPI_COMM_RANK(lgrp,idproc,ierr)
c find number of processors
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
c copy directly on processor 0
      if (idproc.eq.0) then
         do 10 j = 1, nxp
         g(j) = f(j)
   10    continue
c send to remaining processors
         do 20 i = 2, nvp
         is = i - 1
         nxps = nxp*is
         koff = min(nxps,nx-1)
         nxps = min(nxp,max(0,nx-nxps))
         call MPI_SEND(f(koff+1),nxps,mreal,is,1,lgrp,ierr)
   20    continue
c receive from processor 0
      else
         nxps = min(nxp,max(0,nx-nxp*idproc))
         call MPI_RECV(g,nxps,mreal,0,1,lgrp,istatus,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine VGATHER(f,g,nx,nxp)
c this subroutine performs a gather of a real vector, that is,
c successive nxp blocks of destination f on processor 0, are received
c from source g on successive processors:
c f(1+nxp*idproc:min(nxp*(idproc+1),nx)) = g(1:nxp)
c f = input real data
c g = output real data
c nx = size of destination data
c nxp = size of source data
      implicit none
      integer nx, nxp
      real f, g
      dimension f(nx), g(nxp)
c common blocks for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer i, j, idproc, nvp, is, nxps, koff, ierr
      integer istatus
      dimension istatus(lstat)
c find processor id
      call MPI_COMM_RANK(lgrp,idproc,ierr)
c find number of processors
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
c copy directly on processor 0
      if (idproc.eq.0) then
         do 10 j = 1, nxp
         f(j) = g(j)
   10    continue
c receive from remaining processors
         do 20 i = 2, nvp
         is = i - 1
         nxps = nxp*is
         koff = min(nxps,nx-1)
         nxps = min(nxp,max(0,nx-nxps))
         call MPI_RECV(f(koff+1),nxps,mreal,is,2,lgrp,istatus,ierr)
   20    continue
c send to processor 0
      else
         nxps = min(nxp,max(0,nx-nxp*idproc))
         call MPI_SEND(g,nxps,mreal,0,2,lgrp,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine END_MPI
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
