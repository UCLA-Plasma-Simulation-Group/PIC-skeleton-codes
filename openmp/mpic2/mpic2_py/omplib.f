c OpenMP utility library
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine INIT_OMP(nth)
c initialize openmp library
c use nth threads if nth > 0; otherwise, use the number found
      implicit none
      integer nth
c get definition of OpenMP functions
      include 'omp_lib.h'
c common block for parallel processing
      integer nthreads
      common /omdata/ nthreads
      save /omdata/
c local data
      integer ncpus
c determine how many processors are available
      ncpus = omp_get_num_procs()
      write (*,*) 'number of cpus found = ', ncpus
      nthreads = omp_get_max_threads()
      write (*,*) 'maximum number of threads = ', nthreads
      if (nth.gt.0) nthreads = nth
      call omp_set_num_threads(nthreads)
      write (*,*) 'using ',  nthreads, ' thread(s)'
      return
      end
c-----------------------------------------------------------------------
      subroutine SETNTHSIZE(nth)
c set number of threads
      implicit none
      integer nth
c common block for parallel processing
      integer nthreads
      common /omdata/ nthreads
      if (nth.gt.0) nthreads = nth
      call omp_set_num_threads(nthreads)
      return
      end
c-----------------------------------------------------------------------
      function GETNTHSIZE()
c get number of threads
      implicit none
      integer GETNTHSIZE
c common block for parallel processing
      integer nthreads
      common /omdata/ nthreads
      GETNTHSIZE = nthreads
      return
      end
