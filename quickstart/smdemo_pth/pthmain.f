c vector add test program for Pthreads
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine PADD(a,b,c,nx)
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
      subroutine INIT_PT(nth,irc)
c initialize multi-tasking library
c use nth threads if nth > 0; otherwise, use the number found
c error code is modified only if there is an error
      implicit none
      integer nth, irc
c common block for parallel processing
      integer maxthreads, nthreads, idtask
c maxthreads = maximum number of threads supported
      parameter(maxthreads=16)
      dimension idtask(maxthreads)
      common /ptdata/ nthreads, idtask
      save /ptdata/
c local data
      integer ncpus
c determine how many processors are available
      call MP_INIT(ncpus,irc)
      nthreads = ncpus
      if (nthreads.eq.0) nthreads = 1
      write (*,*) 'number of cpus found = ', nthreads
      if (nth.gt.0) nthreads = nth
      write (*,*) 'using ',  nthreads, ' cpu(s)'
      return
      end
c-----------------------------------------------------------------------
      subroutine SETNTHSIZE(nth)
c set number of threads
      implicit none
      integer nth
c common block for parallel processing
      integer maxthreads, nthreads, idtask
      parameter(maxthreads=16)
      dimension idtask(maxthreads)
      common /ptdata/ nthreads, idtask
      if (nth.gt.0) nthreads = nth
      return
      end
c-----------------------------------------------------------------------
      subroutine PTADD(a,b,c,nx,irc)
c multitasking vector add
c irc = ierror indicator (0 = no error)
      implicit none
      integer nx, irc
      real a, b, c
      dimension a(nx), b(nx), c(nx)
c common block for parallel processing
      integer maxthreads, nthreads, idtask
      parameter(maxthreads=16)
      dimension idtask(maxthreads)
      common /ptdata/ nthreads, idtask
c local data
      integer i, nmt, nxp, nxo
      external PADD
      nmt = nthreads - 1
      nxp = (nx - 1)/nthreads + 1
c start add tasks
      do 10 i = 1, nmt
      nxo = nxp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),PADD,4,a(nxo),b(nxo),c(nxo),nxp)
c check for errors
      if (idtask(i).eq.0) then
        irc = -1
        return
      endif
   10 continue
c add remaining data
      nxo = nxp*nmt + 1
      nxp = nx - nxp*nmt
      call PADD(a(nxo),b(nxo),c(nxo),nxp)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) irc = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine END_PT()
      implicit none
c terminate pthreads library
      call MP_END
      return
      end
