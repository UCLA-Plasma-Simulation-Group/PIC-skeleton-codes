!-----------------------------------------------------------------------
! Basic parallel PIC library for MPI communications
! pplib2.f90 contains basic communications procedures for 1d partitions:
! PPINIT2 initializes parallel processing for Fortran90, returns
!         number of processors and processor id.
! PPFNDGRP finds which MPI nodes have the same host name, creates
!          ordered list of such nodes, and returns where current node is
!          located in that list 
! PPEXIT terminates parallel processing.
! PPABORT aborts parallel processing.
! PWTIMERA performs parallel local wall clock timing.
! PPSUM performs parallel sum of a real vector.
! PPDSUM performs parallel sum of a double precision vector.
! PPIMAX performs parallel maximum of an integer vector.
! PPPNCGUARD2L sends/receives guard cells in y for scalar array, linear
!              interpolation, and distributed data with non-uniform
!              partition.  for copying guard cells.
! PPPNAGUARD2L sends/receives guard cells in y for scalar array, linear
!              interpolation, and distributed data with non-uniform
!              partition.  for adding guard cells.
! PPPTPOSE performs a transpose of a complex scalar array, distributed
!          in y, to a complex scalar array, distributed in x.
!          optimized for GPU
! PPPTPOSEN performs a transpose of an n component complex vector array,
!           distributed in y, to an n component complex vector array,
!           distributed in x.  optimized for GPU with vector data
! ACSNDREC helps perform asynchronous transpose between GPUS
! PPPMOVE2 moves particles into appropriate spatial regions for tiled
!         distributed data.
! written by viktor k. decyk, ucla
! copyright 1995, regents of the university of california
! update: may 7, 2014
      module pplib2
      use mpi
      implicit none
!     include 'mpif.h'
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
      public :: PPINIT2, PPFNDGRP, PPEXIT, PPABORT, PWTIMERA
      public :: PPSUM, PPDSUM, PPIMAX
      public :: PPPNCGUARD2L, PPPNAGUARD2L
      public :: PPPTPOSE, PPPTPOSEN, PPPMOVE2, ACSNDREC
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
      integer :: idproc, nvp
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
      subroutine PPFNDGRP(locl,kstrt,nvp,idev,ndev)
! this subroutine finds which MPI nodes have the same host name
! creates ordered list of such nodes andreturns where current node is
! located in that list 
! used to ensure that different GPUs on the same host have different ids
! input: all except locl, idev, ndev, output: locl, idev, ndev
! locl = ordered list of ranks on same host
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idev = location in rank list with value kstrt - 1 (-1 if not found)
! ndev = number of nodes found on host
      implicit none
      integer, intent(in) :: kstrt, nvp
      integer :: idev, ndev
      integer, dimension(nvp) :: locl
! lname = length of hostname character
! lgrp = current communicator
! local data
      integer :: j, l, m, n, nn, ks, id, lname, mchar
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      character(len=MPI_MAX_PROCESSOR_NAME) :: hname, name
      lname = MPI_MAX_PROCESSOR_NAME
      mchar = MPI_CHARACTER
      ks = kstrt - 1
      call MPI_GET_PROCESSOR_NAME(hname,l,ierr)
      nn = 0
! this segment is used for mpi computers
! find and save ranks with same host name
      do n = 1, nvp
         id = n - ks - 1
         if (id < 0) id = id + nvp
         if (id /= ks) then
! post receive
            call MPI_IRECV(name,lname,mchar,id,n,lgrp,msid,ierr)
! send data
            call MPI_SEND(hname,lname,mchar,id,n,lgrp,ierr)
! receive data
            call MPI_WAIT(msid,istatus,ierr)
         else
            name = hname
         endif
! save rank if remote name equals local name
         if (name==hname) then
            nn = nn + 1
            locl(nn) = id
         endif
      enddo
! order rank list
      do j = 1, nn-1
         l = j
         m = locl(j)
! find minimum value and location
         do n = j+1, nn
            id = locl(n)
            if (id < m) then
               m = id
               l = n
            endif
         enddo
! swap minimum to beginning of array
         if (l > j) then
            id = locl(j)
            locl(j) = m
            locl(l) = id
         endif
      enddo
! find location in rank list with value kstrt - 1
      idev = 0
      do j = 1, nn
         if (ks.eq.locl(j)) idev = j
      enddo
      idev = idev - 1
! return number of nodes found on host
      ndev = nn
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
      real :: time
      double precision :: dtime
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
      real, dimension(nxp) :: f, g
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
      double precision, dimension(nxp) :: f, g
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
      integer, dimension(nxp) :: if, ig
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
      subroutine PPPNCGUARD2L(scs,scr,kstrt,nvp,nxv)
! this subroutine sends/receives guard cells in y for scalar array
! for copying guard cells, sends to left processor, receives from right
! output: scr
! scs(j) = input data to send
! scr(j) = received data
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nxv = size of array to send
! linear interpolation, for distributed data
      implicit none
      integer, intent(in) :: kstrt, nvp, nxv
      real, dimension(nxv) :: scs, scr
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: ks, moff, kl, kr
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      ks = kstrt - 1
      moff = nxv*nvp + 2
! copy to guard cells
      kr = ks + 1
      if (kr >= nvp) kr = kr - nvp
      kl = ks - 1
      if (kl < 0)  kl = kl + nvp
! this segment is used for mpi computers
      call MPI_IRECV(scr,nxv,mreal,kr,moff,lgrp,msid,ierr)
      call MPI_SEND(scs,nxv,mreal,kl,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPNAGUARD2L(scs,scr,kstrt,nvp,nxv)
! this subroutine sends/receives guard cells in y for scalar array
! for adding guard cells, sends to right processor, receives from left
! output: scr
! scs(j) = input data to send
! scr(j) = received data
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nxv = size of array to send
! linear interpolation, for distributed data
      implicit none
      integer, intent(in) :: kstrt, nvp, nxv
      real, dimension(nxv) :: scs, scr
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: ks, moff, kl, kr
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
! special case for one processor
      ks = kstrt - 1
      moff = nxv*nvp + 1
! add guard cells
      kr = ks + 1
      if (kr >= nvp) kr = kr - nvp
      kl = ks - 1
      if (kl < 0) kl = kl + nvp
! this segment is used for mpi computers
      call MPI_IRECV(scr,nxv,mreal,kl,moff,lgrp,msid,ierr)
      call MPI_SEND(scs,nxv,mreal,kr,moff,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPTPOSE(sm,tm,nx,ny,kxp,kyp,kstrt,nvp)
! this subroutine sends and receives data between MPI nodes to perform
! a transpose of a matrix distributed in y, to another matrix
! distributed in x.  one message is sent and received at a time.
! optimized for GPU
! ss/tm = complex buffers on host to be sent/received
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
      implicit none
      integer, intent(in) ::  nx, ny, kxp, kyp, kstrt, nvp
      complex, dimension(kxp*kyp,nvp) :: sm, tm
! lgrp = current communicator
! mcplx = default datatype for complex
! local data
      integer :: j, n, nn, ks, kyps, kxyp, id, joff, ld
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      ks = kstrt - 1
      kyps = min(kyp,max(0,ny-kyp*ks))
      kxyp = kxp*kyp
! special case for one processor
      if (nvp==1) then
         do j = 1, kxyp
            tm(j,1) = sm(j,1)
         enddo
         return
      endif
      nn = 0
! this segment is used for mpi computers
      do n = 1, nvp
         id = n - ks - 1
         if (id < 0) id = id + nvp
         if (id /= ks) then
! adjust counter to omit data sent to oneself
            nn = nn + 1
! calculate length of data to send
            joff = kxp*id
            ld = kyps*min(kxp,max(0,nx-joff))
! post receive
            call MPI_IRECV(tm(1,nn),kxyp,mcplx,id,n,lgrp,msid,ierr)
! send data
            call MPI_SEND(sm(1,nn),ld,mcplx,id,n,lgrp,ierr)
! receive data
            call MPI_WAIT(msid,istatus,ierr)
         endif
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPTPOSEN(sm,tm,nx,ny,kxp,kyp,kstrt,nvp,ndim)
! this subroutine sends and receives data between MPI nodes to perform
! a transpose of an n component matrix distributed in y, to another 
! n component matrix distributed in x.
! one message is sent and received at a time.
! optimized for GPU with vector data
! ss/tm = complex buffers on host to be sent/received
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      complex, dimension(kxp*ndim*kyp,nvp) :: sm, tm
! lgrp = current communicator
! mcplx = default datatype for complex
! local data
      integer :: j, n, nn, ks, kyps, kxyp, id, joff, ld
      integer :: msid, ierr
      integer, dimension(lstat) :: istatus
      ks = kstrt - 1
      kyps = ndim*min(kyp,max(0,ny-kyp*ks))
      kxyp = kxp*ndim*kyp
! special case for one processor
      if (nvp==1) then
         do j = 1, kxyp
            tm(j,1) = sm(j,1)
         enddo
         return
      endif
      nn = 0
! this segment is used for mpi computers
      do n = 1, nvp
         id = n - ks - 1
         if (id < 0) id = id + nvp
         if (id /= ks) then
! adjust counter to omit data sent to oneself
            nn = nn + 1
! calculate length of data to send
            joff = kxp*id
            ld = kyps*min(kxp,max(0,nx-joff))
! post receive
            call MPI_IRECV(tm(1,nn),kxyp,mcplx,id,n,lgrp,msid,ierr)
! send data
            call MPI_SEND(sm(1,nn),ld,mcplx,id,n,lgrp,ierr)
! receive data
            call MPI_WAIT(msid,istatus,ierr)
         endif
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine ACSNDREC(stm,idproc,nsize,ntag,mode)
! this subroutine is part of a family of procedures for performing
! a transpose by sending and receiving asynchronous data between GPUs
! stm = receive/send buffer
! idproc = processor id for sending/receiving
! nsize = size of data packet in words
! ntag = MPI tag
! mode = (1,2,3) = (post receive, post send, wait for send/receive)
! modes 1 and 2 should be called before mode 3 is called.
      implicit none
      integer, intent(in) :: idproc, nsize, ntag, mode
      complex, dimension(nsize) :: stm
! lgrp = current communicator
! mcplx = default datatype for complex
! local data
      integer, save :: msid, mrid
      integer :: ierr
      integer, dimension(lstat) :: istatus
      if (mode==1) then
         call MPI_IRECV(stm,nsize,mcplx,idproc,ntag,lgrp,mrid,ierr)
      else if (mode==2) then
         call MPI_ISEND(stm,nsize,mcplx,idproc,ntag,lgrp,msid,ierr)
      else if (mode==3) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(mrid,istatus,ierr)
      endif
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
      real, dimension(idimp*nbmax) :: sbufl, sbufr, rbufl, rbufr
      integer, dimension(3,mx1) :: ncll, nclr, mcll, mclr
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
         do j = 1, idimp*nclr(3,mx1)
            rbufl(j) = sbufr(j)
         enddo
         do j = 1, idimp*ncll(3,mx1)
            rbufr(j) = sbufl(j)
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
      use pplib2, only: SUB => PPINIT2
      implicit none
      integer :: idproc, nvp
      call SUB(idproc,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPFNDGRP(locl,kstrt,nvp,idev,ndev)
      use pplib2, only: SUB => PPFNDGRP
      implicit none
      integer :: kstrt, nvp, idev, ndev
      integer, dimension(nvp) :: locl
      call SUB(locl,kstrt,nvp,idev,ndev)
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
      integer :: icntrl
      real :: time
      double precision :: dtime
      call SUB(icntrl,time,dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPSUM(f,g,nxp)
      use pplib2, only: SUB => PPSUM
      implicit none
      integer :: nxp
      real, dimension(nxp) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPDSUM(f,g,nxp)
      use pplib2, only: SUB => PPDSUM
      implicit none
      integer :: nxp
      double precision, dimension(nxp) :: f, g
      call SUB(f,g,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPIMAX(if,ig,nxp)
      use pplib2, only: SUB => PPIMAX
      implicit none
      integer :: nxp
      integer, dimension(nxp) :: if, ig
      call SUB(if,ig,nxp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPNCGUARD2L(scs,scr,kstrt,nvp,nxv)
      use pplib2, only: SUB => PPPNCGUARD2L
      implicit none
      integer :: kstrt, nvp, nxv
      real, dimension(nxv) :: scs, scr
      call SUB(scs,scr,kstrt,nvp,nxv)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPNAGUARD2L(scs,scr,kstrt,nvp,nxv)
      use pplib2, only: SUB => PPPNAGUARD2L
      implicit none
      integer:: kstrt, nvp, nxv
      real, dimension(nxv) :: scs, scr
      call SUB(scs,scr,kstrt,nvp,nxv)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPTPOSE(sm,tm,nx,ny,kxp,kyp,kstrt,nvp)
      use pplib2, only: SUB => PPPTPOSE
      implicit none
      integer ::  nx, ny, kxp, kyp, kstrt, nvp
      complex, dimension(kxp*kyp,nvp) :: sm, tm
      call SUB(sm,tm,nx,ny,kxp,kyp,kstrt,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPTPOSEN(sm,tm,nx,ny,kxp,kyp,kstrt,nvp,ndim)
      use pplib2, only: SUB => PPPTPOSEN
      implicit none
      integer :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      complex, dimension(kxp*ndim*kyp,nvp) :: sm, tm
      call SUB(sm,tm,nx,ny,kxp,kyp,kstrt,nvp,ndim)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine ACSNDREC(stm,idproc,nsize,ntag,mode)
      use pplib2, only: SUB => ACSNDREC
      integer :: idproc, nsize, ntag, mode
      complex, dimension(nsize) :: stm
      call SUB(stm,idproc,nsize,ntag,mode)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,  &
     &kstrt,nvp,idimp,nbmax,mx1)
      use pplib2, only: SUB => PPPMOVE2
      implicit none
      integer, intent(in) :: kstrt, nvp, idimp, nbmax, mx1
      real, dimension(idimp,nbmax) :: sbufl, sbufr, rbufl, rbufr
      integer, dimension(3,mx1) :: ncll, nclr, mcll, mclr
      call SUB(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,nvp,   &
     &idimp,nbmax,mx1)
      end subroutine

