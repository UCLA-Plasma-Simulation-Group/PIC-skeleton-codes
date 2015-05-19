c-----------------------------------------------------------------------
c Basic parallel PIC library for MPI communications
c pplib2.f contains basic communications procedures for 1d partitions:
c PPINIT2 initializes parallel processing for Fortran90, returns
c         number of processors and processor id.
c PPFNDGRP finds which MPI nodes have the same host name, creates
c          ordered list of such nodes, and returns where current node is
c          located in that list 
c PPEXIT terminates parallel processing.
c PPABORT aborts parallel processing.
c PWTIMERA performs parallel local wall clock timing.
c PPSUM performs parallel sum of a real vector.
c PPDSUM performs parallel sum of a double precision vector.
c PPIMAX performs parallel maximum of an integer vector.
c PPPCNCGUARD2L sends/receives guard cells in y for complex scalar array
c               linear interpolation, and distributed data with
c               non-uniform partition.  for copying guard cells.
c PPPCNAGUARD2L sends/receives guard cells in y for complex scalar array
c               linear interpolation, and distributed data with
c               non-uniform partition.  for adding guard cells.
c PPPTPOSE performs a transpose of a complex scalar array, distributed
c          in y, to a complex scalar array, distributed in x.
c          optimized for GPU
c PPPTPOSEN performs a transpose of an n component complex vector array,
c           distributed in y, to an n component complex vector array,
c           distributed in x.  optimized for GPU with vector data
c ACSNDREC helps perform asynchronous transpose between GPUS
c PPPMOVE2 moves particles into appropriate spatial regions for tiled
c         distributed data.
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: may 7, 2014
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
      subroutine PPFNDGRP(locl,kstrt,nvp,idev,ndev)
c this subroutine finds which MPI nodes have the same host name
c creates ordered list of such nodes and returns where current node is
c located in that list 
c used to ensure that different GPUs on the same host have different ids
c input: all except locl, idev, ndev, output: locl, idev, ndev
c locl = ordered list of ranks on same host
c kstrt = starting data block number
c nvp = number of real or virtual processors
c idev = location in rank list with value kstrt - 1 (-1 if not found)
c ndev = number of nodes found on host
      implicit none
      integer kstrt, nvp, idev, ndev
      integer locl
      dimension locl(nvp)
c get definition of MPI constants
      include 'mpif.h'
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
c lname = length of hostname character
      parameter(lstat=10)
c lgrp = current communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer j, l, m, n, nn, ks, id, lname, mchar
      integer ierr, msid, istatus
      dimension istatus(lstat)
      character*(MPI_MAX_PROCESSOR_NAME) :: hname, name
      lname = MPI_MAX_PROCESSOR_NAME
      mchar = MPI_CHARACTER
      ks = kstrt - 1
      call MPI_GET_PROCESSOR_NAME(hname,l,ierr)
      nn = 0
c this segment is used for mpi computers
c find and save ranks with same host name
      do 10 n = 1, nvp
      id = n - ks - 1
      if (id.lt.0) id = id + nvp
      if (id.ne.ks) then
c post receive
         call MPI_IRECV(name,lname,mchar,id,n,lgrp,msid,ierr)
c send data
         call MPI_SEND(hname,lname,mchar,id,n,lgrp,ierr)
c receive data
         call MPI_WAIT(msid,istatus,ierr)
      else
         name = hname
      endif
c save rank if remote name equals local name
      if (name.eq.hname) then
         nn = nn + 1
         locl(nn) = id
      endif
   10 continue
c order rank list
      do 30 j = 1, nn-1
      l = j
      m = locl(j)
c find minimum value and location
      do 20 n = j+1, nn
      id = locl(n)
      if (id.lt.m) then
         m = id
         l = n
      endif
   20 continue
c swap minimum to beginning of array
      if (l.gt.j) then
         id = locl(j)
         locl(j) = m
         locl(l) = id
      endif
   30 continue
c find location in rank list with value kstrt - 1
      idev = 0
      do 40 j = 1, nn
      if (ks.eq.locl(j)) idev = j
   40 continue
      idev = idev - 1
c return number of nodes found on host
      ndev = nn
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
      subroutine PPPCNCGUARD2L(scs,scr,kstrt,nvp,nxvh)
c this subroutine sends/receives guard cells in y for scalar array
c for copying guard cells, sends to left processor, receives from right
c output: scr
c scs(j) = input data to send
c scr(j) = received data
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nxvh = size of complex array to send
c linear interpolation, for distributed data
      implicit none
      integer kstrt, nvp, nxvh
      complex scs, scr
      dimension scs(nxvh), scr(nxvh)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, moff, kl, kr
      integer istatus, msid, ierr
      dimension istatus(lstat)
      ks = kstrt - 1
      moff = 2*nxvh*nvp + 2
c copy to guard cells
      kr = ks + 1
      if (kr.ge.nvp) kr = kr - nvp
      kl = ks - 1
      if (kl.lt.0)  kl = kl + nvp
c this segment is used for mpi computers
      call MPI_IRECV(scr,nxvh,mcplx,kr,moff,lgrp,msid,ierr)
      call MPI_SEND(scs,nxvh,mcplx,kl,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPCNAGUARD2L(scs,scr,kstrt,nvp,nxvh)
c this subroutine sends/receives guard cells in y for scalar array
c for adding guard cells, sends to right processor, receives from left
c output: scr
c scs(j) = input data to send
c scr(j) = received data
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nxvh = size of complex array to send
c linear interpolation, for distributed data
      implicit none
      integer kstrt, nvp, nxvh
      complex scs, scr
      dimension scs(nxvh), scr(nxvh)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, moff, kl, kr
      integer istatus, msid, ierr
      dimension istatus(lstat)
c special case for one processor
      ks = kstrt - 1
      moff = 2*nxvh*nvp + 1
c add guard cells
      kr = ks + 1
      if (kr.ge.nvp) kr = kr - nvp
      kl = ks - 1
      if (kl.lt.0) kl = kl + nvp
c this segment is used for mpi computers
      call MPI_IRECV(scr,nxvh,mcplx,kl,moff,lgrp,msid,ierr)
      call MPI_SEND(scs,nxvh,mcplx,kr,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPTPOSE(sm,tm,nx,ny,kxp,kyp,kstrt,nvp)
c this subroutine sends and receives data between MPI nodes to perform
c a transpose of a matrix distributed in y, to another matrix
c distributed in x.  one message is sent and received at a time.
c optimized for GPU
c ss/tm = complex buffers on host to be sent/received
c nx/ny = number of points in x/y
c kxp/kyp = number of data values per block in x/y
c kstrt = starting data block number
c nvp = number of real or virtual processors
      implicit none
      integer nx, ny, kxp, kyp, kstrt, nvp
      complex sm, tm
      dimension sm(kxp*kyp,nvp), tm(kxp*kyp,nvp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer j, n, nn, ks, kyps, kxyp, id, joff, ld
      integer ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 1
      kyps = min(kyp,max(0,ny-kyp*ks))
      kxyp = kxp*kyp
c special case for one processor
      if (nvp.eq.1) then
         do 10 j = 1, kxyp
         tm(j,1) = sm(j,1)
   10    continue
         return
      endif
      nn = 0
c this segment is used for mpi computers
      do 20 n = 1, nvp
      id = n - ks - 1
      if (id.lt.0) id = id + nvp
      if (id.ne.ks) then
c adjust counter to omit data sent to oneself
         nn = nn + 1
c calculate length of data to send
         joff = kxp*id
         ld = kyps*min(kxp,max(0,nx-joff))
c post receive
         call MPI_IRECV(tm(1,nn),kxyp,mcplx,id,n,lgrp,msid,ierr)
c send data
         call MPI_SEND(sm(1,nn),ld,mcplx,id,n,lgrp,ierr)
c receive data
         call MPI_WAIT(msid,istatus,ierr)
      endif
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPTPOSEN(sm,tm,nx,ny,kxp,kyp,kstrt,nvp,ndim)
c this subroutine sends and receives data between MPI nodes to perform
c a transpose of an n component matrix distributed in y, to another 
c n component matrix distributed in x.
c one message is sent and received at a time.
c optimized for GPU with vector data
c ss/tm = complex buffers on host to be sent/received
c nx/ny = number of points in x/y
c kxp/kyp = number of data values per block in x/y
c kstrt = starting data block number
c nvp = number of real or virtual processors
      implicit none
      integer nx, ny, kxp, kyp, kstrt, nvp, ndim
      complex sm, tm
      dimension sm(kxp*ndim*kyp,nvp), tm(kxp*ndim*kyp,nvp)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer j, n, nn, ks, kyps, kxyp, id, joff, ld
      integer ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 1
      kyps = ndim*min(kyp,max(0,ny-kyp*ks))
      kxyp = kxp*ndim*kyp
c special case for one processor
      if (nvp.eq.1) then
         do 10 j = 1, kxyp
         tm(j,1) = sm(j,1)
   10    continue
         return
      endif
      nn = 0
c this segment is used for mpi computers
      do 20 n = 1, nvp
      id = n - ks - 1
      if (id.lt.0) id = id + nvp
      if (id.ne.ks) then
c adjust counter to omit data sent to oneself
         nn = nn + 1
c calculate length of data to send
         joff = kxp*id
         ld = kyps*min(kxp,max(0,nx-joff))
c post receive
         call MPI_IRECV(tm(1,nn),kxyp,mcplx,id,n,lgrp,msid,ierr)
c send data
         call MPI_SEND(sm(1,nn),ld,mcplx,id,n,lgrp,ierr)
c receive data
         call MPI_WAIT(msid,istatus,ierr)
      endif
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine ACSNDREC(stm,idproc,nsize,ntag,mode)
c this subroutine is part of a family of procedures for performing
c a transpose by sending and receiving asynchronous data between GPUs
c stm = receive/send buffer
c idproc = processor id for sending/receiving
c nsize = size of data packet in words
c ntag = MPI tag
c mode = (1,2,3) = (post receive, post send, wait for send/receive)
c modes 1 and 2 should be called before mode 3 is called.
      implicit none
      integer idproc, nsize, ntag, mode
      complex stm
      dimension stm(nsize)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, mrid, istatus
      dimension istatus(lstat)
      save msid, mrid
      if (mode.eq.1) then
         call MPI_IRECV(stm,nsize,mcplx,idproc,ntag,lgrp,mrid,ierr)
      else if (mode.eq.2) then
         call MPI_ISEND(stm,nsize,mcplx,idproc,ntag,lgrp,msid,ierr)
      else if (mode.eq.3) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(mrid,istatus,ierr)
      endif
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
      dimension sbufl(idimp*nbmax), sbufr(idimp*nbmax)
      dimension rbufl(idimp*nbmax), rbufr(idimp*nbmax)
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
         do 50 j = 1, idimp*nclr(3,mx1)
         rbufl(j) = sbufr(j)
   50    continue
         do 60 j = 1, idimp*ncll(3,mx1)
         rbufr(j) = sbufl(j)
   60    continue
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

