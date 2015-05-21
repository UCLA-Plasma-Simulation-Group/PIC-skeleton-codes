!-----------------------------------------------------------------------
! Skeleton 1D Electrostatic OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpic1
      use mpush1_h
      use omplib_h
      implicit none
! indx = exponent which determines grid points in x direction:
! nx = 2**indx.
      integer, parameter :: indx =   9
! npx = number of electrons distributed in x direction.
      integer, parameter :: npx = 18432
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
! vtx = thermal velocity of electrons in x direction
! vx0 = drift velocity of electrons in x direction.
      real, parameter :: vtx = 1.0, vx0 = 0.0
! ax = smoothed particle size in x direction
      real :: ax = .912871
! idimp = number of particle coordinates = 2
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 2, ipbc = 1
! wke/we/wt = particle kinetic/electric field/total energy
      real :: wke = 0.0, we = 0.0, wt = 0.0
! mx = number of grids in x in sorting tiles
      integer :: mx = 32
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! declare scalars for standard code
      integer :: np, nx, nxh, nxe
      integer :: mx1, ntime, nloop, isign
      real :: qbme, affp
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, ntmax, npbmx, irc
      integer :: nvp
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), pointer :: part
! qe = electron charge density with guard cells
      real, dimension(:), pointer :: qe
! fxe = smoothed electric field with guard cells
      real, dimension(:), pointer :: fxe
! ffc = form factor array for poisson solver
      complex, dimension(:), pointer :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct
!
! declare arrays for OpenMP (tiled) code:
! ppart = tiled particle array
! ppbuff = buffer array for reordering tiled particle array
      real, dimension(:,:,:), pointer :: ppart, ppbuff
! kpic = number of particles in each tile
      integer, dimension(:), pointer :: kpic
! ncl = number of particles departing tile in each direction
      integer, dimension(:,:), pointer :: ncl
! ihole = location/destination of each particle departing tile
      integer, dimension(:,:,:), pointer :: ihole
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
      real :: tpush = 0.0, tsort = 0.0
      double precision :: dtime
!
      irc = 0
! nvp = number of shared memory nodes (0=default)
      nvp = 0
!     write (*,*) 'enter number of nodes:'
!     read (5,*) nvp
! initialize for shared memory parallel processing
      call INIT_OMP(nvp)
!
! initialize scalars for standard code
! np = total number of particles in simulation
! nx = number of grid points in x direction
      np = npx; nx = 2**indx; nxh = nx/2
      nxe = nx + 2
! mx1 = number of tiles in x direction
      mx1 = (nx - 1)/mx + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx)/real(np)
!
! allocate data for standard code
      allocate(part(idimp,np))
      allocate(qe(nxe),fxe(nxe))
      allocate(ffc(nxh),mixup(nxh),sct(nxh))
      allocate(kpic(mx1))
!
! prepare fft tables
      call WFFT1RINIT(mixup,sct,indx,nxh)
! calculate form factors
      isign = 0
      call POIS1(qe,fxe,isign,ffc,ax,affp,we,nx)
! initialize electrons
      call DISTR1(part,vtx,vx0,npx,idimp,np,nx,ipbc)
!
! find number of particles in each of mx, tiles: updates kpic, nppmx
      call DBLKP1L(part,kpic,nppmx,idimp,np,mx,mx1,irc)
      if (irc /= 0) then
         write (*,*) 'DBLKP1L error, irc=', irc
         stop
      endif
! allocate vector particle data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmax = xtras*nppmx
      npbmx = xtras*nppmx
      allocate(ppart(idimp,nppmx0,mx1))
      allocate(ppbuff(idimp,npbmx,mx1))
      allocate(ncl(2,mx1))
      allocate(ihole(2,ntmax+1,mx1))
! copy ordered particle data for OpenMP: updates ppart and kpic
      call PPMOVIN1L(part,ppart,kpic,nppmx0,idimp,np,mx,mx1,irc)
      if (irc /= 0) then
         write (*,*) 'PPMOVIN1L overflow error, irc=', irc
         stop
      endif
! sanity check
      call PPCHECK1L(ppart,kpic,idimp,nppmx0,nx,mx,mx1,irc)
      if (irc /= 0) then
         write (*,*) 'PPCHECK1L error: irc=', irc
         stop
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     write (*,*) 'ntime = ', ntime
!
! deposit charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call GPPOST1L(ppart,qe,kpic,qme,nppmx0,idimp,mx,nxe,mx1)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      call AGUARD1L(qe,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with standard procedure:
! updates qe, fxe
      call dtimer(dtime,itime,-1)
      isign = -1
      call FFT1RXX(qe,fxe,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! calculate force/charge in fourier space with standard procedure:
! updates fxe, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call POIS1(qe,fxe,isign,ffc,ax,affp,we,nx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform force to real space with standard procedure: updates fxe, qe
      call dtimer(dtime,itime,-1)
      isign = 1
      call FFT1RXX(fxe,qe,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with standard procedure: updates fxe
      call dtimer(dtime,itime,-1)
      call CGUARD1L(fxe,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! push particles with OpenMP:
      wke = 0.0
      call dtimer(dtime,itime,-1)
! updates part, wke
!     call GPPUSH1L(ppart,fxe,kpic,qbme,dt,wke,idimp,nppmx0,nx,mx,nxe,  &
!    &mx1,ipbc)
! updates ppart, ncl, ihole, wke, irc
      call GPPUSHF1L(ppart,fxe,kpic,ncl,ihole,qbme,dt,wke,idimp,nppmx0, &
     &nx,mx,nxe,mx1,ntmax,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
      if (irc /= 0) then
         write (*,*) 'GPPUSHF1L error: irc=', irc
         stop
      endif
!
! reorder particles by tile with OpenMP:
      call dtimer(dtime,itime,-1)
! updates ppart, ppbuff, kpic, ncl, ihole, and irc
!     call PPORDER1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,mx,mx1,&
!    &npbmx,ntmax,irc)
! updates ppart, ppbuff, kpic, ncl, and irc
      call PPORDERF1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,npbmx&
     &,ntmax,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tsort = tsort + time
      if (irc /= 0) then
         write (*,*) 'PPORDERF1L error: ntmax, irc=', ntmax, irc
         stop
      endif
!
      if (ntime==0) then
         write (*,*) 'Initial Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') we, wke, wke + we
      endif
      ntime = ntime + 1
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
      write (*,*) 'ntime = ', ntime
      write (*,*) 'Final Field, Kinetic and Total Energies:'
      write (*,'(3e14.7)') we, wke, wke + we
!
      write (*,*)
      write (*,*) 'deposit time = ', tdpost
      write (*,*) 'guard time = ', tguard
      write (*,*) 'solver time = ', tfield
      write (*,*) 'fft time = ', tfft
      write (*,*) 'push time = ', tpush
      write (*,*) 'sort time = ', tsort
      tfield = tfield + tguard + tfft
      write (*,*) 'total solver time = ', tfield
      time = tdpost + tpush + tsort
      write (*,*) 'total particle time = ', time
      wt = time + tfield
      write (*,*) 'total time = ', wt
      write (*,*)
!
      wt = 1.0e+09/(real(nloop)*real(np))
      write (*,*) 'Push Time (nsec) = ', tpush*wt
      write (*,*) 'Deposit Time (nsec) = ', tdpost*wt
      write (*,*) 'Sort Time (nsec) = ', tsort*wt
      write (*,*) 'Total Particle Time (nsec) = ', time*wt
      write (*,*)
!
      stop
      end program
