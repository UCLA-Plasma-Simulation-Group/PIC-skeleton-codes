!-----------------------------------------------------------------------
! Skeleton 2D Electrostatic OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpic2
      use mpush2_h
      use omplib_h
      implicit none
! indx/indy = exponent which determines grid points in x/y direction:
! nx = 2**indx, ny = 2**indy.
      integer, parameter :: indx =   9, indy =   9
! npx/npy = number of electrons distributed in x/y direction.
      integer, parameter :: npx =  3072, npy =   3072
! ndim = number of velocity coordinates = 2
      integer, parameter :: ndim = 2
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
! vtx/vty = thermal velocity of electrons in x/y direction
! vx0/vy0 = drift velocity of electrons in x/y direction.
      real, parameter :: vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0
! ax/ay = smoothed particle size in x/y direction
      real :: ax = .912871, ay = .912871
! idimp = number of particle coordinates = 4
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 4, ipbc = 1
! wke/we/wt = particle kinetic/electric field/total energy
      real :: wke = 0.0, we = 0.0, wt = 0.0
! mx/my = number of grids in x/y in sorting tiles
      integer :: mx = 16, my = 16
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! declare scalars for standard code
      integer :: np, nx, ny, nxh, nyh, nxe, nye, nxeh, nxyh, nxhy
      integer :: mx1, my1, mxy1, ntime, nloop, isign
      real :: qbme, affp
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, ntmax, npbmx, irc
      integer :: nvp
!
! declare arrays for standard code:
! part = original particle array
      real, dimension(:,:), pointer :: part
! qe = electron charge density with guard cells
      real, dimension(:,:), pointer :: qe
! fxye = smoothed electric field with guard cells
      real, dimension(:,:,:), pointer :: fxye
! ffc = form factor array for poisson solver
      complex, dimension(:,:), pointer :: ffc
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
! nx/ny = number of grid points in x/y direction
      np = npx*npy; nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = ny/2
      nxe = nx + 2; nye = ny + 1; nxeh = nxe/2
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny)
! mx1/my1 = number of tiles in x/y direction
      mx1 = (nx - 1)/mx + 1; my1 = (ny - 1)/my + 1; mxy1 = mx1*my1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx*ny)/real(np)
!
! allocate data for standard code
      allocate(part(idimp,np)) 
      allocate(qe(nxe,nye),fxye(ndim,nxe,nye))
      allocate(ffc(nxh,nyh),mixup(nxhy),sct(nxyh))
      allocate(kpic(mxy1))
!
! prepare fft tables
      call WFFT2RINIT(mixup,sct,indx,indy,nxhy,nxyh)
! calculate form factors
      isign = 0
      call MPOIS22(qe,fxye,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,  &
     &nyh)
! initialize electrons
      call DISTR2(part,vtx,vty,vx0,vy0,npx,npy,idimp,np,nx,ny,ipbc)
!
! find number of particles in each of mx, my tiles: updates kpic, nppmx
      call DBLKP2L(part,kpic,nppmx,idimp,np,mx,my,mx1,mxy1,irc)
      if (irc /= 0) then
         write (*,*) 'DBLKP2L error, irc=', irc
         stop
      endif
! allocate vector particle data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmax = xtras*nppmx
      npbmx = xtras*nppmx
      allocate(ppart(idimp,nppmx0,mxy1))
      allocate(ppbuff(idimp,npbmx,mxy1))
      allocate(ncl(8,mxy1))
      allocate(ihole(2,ntmax+1,mxy1))
! copy ordered particle data for OpenMP: updates ppart and kpic
      call PPMOVIN2L(part,ppart,kpic,nppmx0,idimp,np,mx,my,mx1,mxy1,irc)
      if (irc /= 0) then
         write (*,*) 'PPMOVIN2L overflow error, irc=', irc
         stop
      endif
! sanity check
      call PPCHECK2L(ppart,kpic,idimp,nppmx0,nx,ny,mx,my,mx1,my1,irc)
      if (irc /= 0) then
         write (*,*) 'PPCHECK2L error: irc=', irc
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
      call GPPOST2L(ppart,qe,kpic,qme,nppmx0,idimp,mx,my,nxe,nye,mx1,   &
     &mxy1)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      call AGUARD2L(qe,nx,ny,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      isign = -1
      call WFFT2RMX(qe,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! calculate force/charge in fourier space with OpenMP: updates fxye, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call MPOIS22(qe,fxye,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,  &
     &nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform force to real space with OpenMP: updates fxye
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT2RM2(fxye,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with OpenMP: updates fxye
      call dtimer(dtime,itime,-1)
      call CGUARD2L(fxye,nx,ny,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! push particles with OpenMP:
      wke = 0.0
      call dtimer(dtime,itime,-1)
! updates ppart, wke
!     call GPPUSH2L(ppart,fxye,kpic,qbme,dt,wke,idimp,nppmx0,nx,ny,mx,my&
!    &,nxe,nye,mx1,mxy1,ipbc)
! updates ppart, ncl, ihole, wke, irc
      call GPPUSHF2L(ppart,fxye,kpic,ncl,ihole,qbme,dt,wke,idimp,nppmx0,&
     &nx,ny,mx,my,nxe,nye,mx1,mxy1,ntmax,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
      if (irc /= 0) then
         write (*,*) 'GPPUSHF2L error: irc=', irc
         stop
      endif
!
! reorder particles by tile with OpenMP:
      call dtimer(dtime,itime,-1)
! updates ppart, ppbuff, kpic, ncl, ihole, and irc
!     call PPORDER2L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny,mx, &
!    &my,mx1,my1,npbmx,ntmax,irc)
! updates ppart, ppbuff, kpic, ncl, and irc
      call PPORDERF2L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,my1, &
     &npbmx,ntmax,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tsort = tsort + time
      if (irc /= 0) then
         write (*,*) 'PPORDERF2L error: ntmax, irc=', ntmax, irc
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
