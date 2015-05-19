!-----------------------------------------------------------------------
! Skeleton 1-2/2D Darwin OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mdpic1
      use mdpush1_h
      use omplib_h
      implicit none
! indx = exponent which determines grid points in x direction:
! nx = 2**indx.
      integer, parameter :: indx =   9
! npx = number of electrons distributed in x direction.
      integer, parameter :: npx =  18432
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
! vtx/vty = thermal velocity of electrons in x/y direction
! vx0/vy0 = drift velocity of electrons in x/y direction.
      real, parameter :: vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0
! vtx/vz0 = thermal/drift velocity of electrons in z direction
      real, parameter :: vtz = 1.0, vz0 = 0.0
! ax = smoothed particle size in x direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ci = 0.1
! idimp = number of particle coordinates = 4
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 4, ipbc = 1
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.4, omy = 0.0, omz = 0.0
! ndc = number of corrections in darwin iteration
      integer :: ndc = 1
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
      real :: zero = 0.0
! mx = number of grids in x in sorting tiles
      integer :: mx = 32
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! declare scalars for standard code
      integer :: k
      integer :: np, nx, nxh, nxe, nxeh
      integer :: mx1, ntime, nloop, isign
      real :: qbme, affp, q2m0, wpm, wpmax, wpmin
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, ntmax, npbmx, irc
      integer :: nvp
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), pointer :: part
! qe = electron charge density with guard cells
! fxe = smoothed longitudinal electric field with guard cells
      real, dimension(:), pointer :: qe, fxe
! cue = electron current density with guard cells
! dcu = acceleration density with guard cells
! cus = transverse electric field
! amu = momentum flux with guard cells
      real, dimension(:,:), pointer :: cue, dcu, cus, amu
! exyze/byze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:), pointer :: exyze, byze
! ffc, ffe = form factor arrays for poisson solvers
      complex, dimension(:), pointer :: ffc, ffe
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct
! gxe, gyze = scratch arrays for fft
      real, dimension(:), pointer :: gxe
      real, dimension(:,:), pointer :: gyze
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
      real :: tdjpost = 0.0, tdcjpost = 0.0, tpush = 0.0, tsort = 0.0
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
      nxe = nx + 2; nxeh = nxe/2
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
      allocate(cue(2,nxe),dcu(2,nxe),cus(2,nxe),amu(2,nxe))
      allocate(exyze(3,nxe),byze(2,nxe))
      allocate(ffc(nxh),ffe(nxh),mixup(nxh),sct(nxh))
      allocate(kpic(mx1))
      allocate(gxe(nxe),gyze(2,nxe))
!
! prepare fft tables
      call WFFT1RINIT(mixup,sct,indx,nxh)
! calculate form factor: ffc
      isign = 0
      call POIS1(qe,fxe,isign,ffc,ax,affp,we,nx)
! initialize electrons
      call DISTR1H(part,vtx,vty,vtz,vx0,vy0,vz0,npx,idimp,np,nx,ipbc)
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
! find maximum and minimum initial electron density
      qe = 0.0
      call GPPOST1L(ppart,qe,kpic,qme,nppmx0,idimp,mx,nxe,mx1)
      call AGUARD1L(qe,nx,nxe)
      call FWPMINMX1(qe,qbme,wpmax,wpmin,nx,nxe)
      wpm = 0.5*(wpmax + wpmin)*affp
! accelerate convergence: update wpm
      if (wpm <= 10.0) wpm = 0.75*wpm
      write (*,*) 'wpm=',wpm
      q2m0 = wpm/affp
! calculate form factor: ffe
      isign = 0
      call EPOIS13(dcu,cus,isign,ffe,ax,affp,wpm,ci,wf,nx,nxeh,nxh)
!
! initialize transverse electric field
      cus = 0.0
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     write (*,*) 'ntime = ', ntime
!
! deposit current with OpenMP: updates cue
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call GJPPOST1L(ppart,cue,kpic,qme,zero,nppmx0,idimp,nx,mx,nxe,mx1,&
     &ipbc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdjpost = tdjpost + time
!
! deposit charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call GPPOST1L(ppart,qe,kpic,qme,nppmx0,idimp,mx,nxe,mx1)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with standard procedure: updates qe, cue
      call dtimer(dtime,itime,-1)
      call AGUARD1L(qe,nx,nxe)
      call ACGUARD1L(cue,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with standard procedure:
! updates qe, gxe
      call dtimer(dtime,itime,-1)
      isign = -1
      call FFT1RXX(qe,gxe,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! calculate longitudinal force/charge in fourier space with standard
! procedure: updates fxe, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call POIS1(qe,fxe,isign,ffc,ax,affp,we,nx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform longitudinal electric force to real space with standard
! procedure: updates fxe, gxe
      call dtimer(dtime,itime,-1)
      isign = 1
      call FFT1RXX(fxe,gxe,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! transform current to fourier space with standard procedure:
! updates cue, gyze
      call dtimer(dtime,itime,-1)
      isign = -1
      call FFT1R2X(cue,gyze,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! calculate magnetic field in fourier space with standard procedure:
! updates byze, wm
      call dtimer(dtime,itime,-1)
      call BBPOIS13(cue,byze,ffc,ci,wm,nx,nxeh,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform magnetic force to real space with standard procedure:
! updates byze, gyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call FFT1R2X(byze,gyze,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! add constant to magnetic field with standard procedure: updates byze
      call dtimer(dtime,itime,-1)
      call BADDEXT1(byze,omy,omz,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! copy guard cells with standard procedure: updates fxe, byze
      call dtimer(dtime,itime,-1)
      call DGUARD1L(fxe,nx,nxe)
      call CGUARD1L(byze,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and old transverse electric fields with standard
! procedure: updates exyze
      call dtimer(dtime,itime,-1)
      call ADDVRFIELD13(exyze,cus,fxe,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! deposit electron acceleration density and momentum flux with OpenMP:
! updates dcu, amu
      call dtimer(dtime,itime,-1)
      dcu = 0.0; amu = 0.0
      call GDJPPOST1L(ppart,exyze,byze,dcu,amu,kpic,omx,qme,qbme,dt,    &
     &idimp,nppmx0,nx,mx,nxe,mx1)
! add old scaled electric field with standard procedure: updates dcu
      call ASCFGUARD1L(dcu,cus,q2m0,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdcjpost = tdcjpost + time
!
! add guard cells with standard procedure: updates dcu, amu
      call dtimer(dtime,itime,-1)
      call ACGUARD1L(dcu,nx,nxe)
      call ACGUARD1L(amu,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform acceleration density and momentum flux to fourier space
! with standard procedure: updates dcu, amu, gyze
      call dtimer(dtime,itime,-1)
      isign = -1
      call FFT1R2X(dcu,gyze,isign,mixup,sct,indx,nxe,nxh)
      call FFT1R2X(amu,gyze,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of time derivative of current with standard
! procedure: updates dcu
      call dtimer(dtime,itime,-1)
      call ADCUPERP13(dcu,amu,nx,nxeh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate transverse electric field with standard procedure:
! updates cus, wf
      call dtimer(dtime,itime,-1)
      isign = -1
      call EPOIS13(dcu,cus,isign,ffe,ax,affp,wpm,ci,wf,nx,nxeh,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform transverse electric field to real space with standard
! procedure: updates cus, gyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call FFT1R2X(cus,gyze,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with standard procedure: updates cus
      call dtimer(dtime,itime,-1)
      call CGUARD1L(cus,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and transverse electric fields with standard
! procedure: exyze = cus + fxe, updates exyze
! cus needs to be retained for next time step
      call dtimer(dtime,itime,-1)
      call ADDVRFIELD13(exyze,cus,fxe,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! inner iteration loop
      do k = 1, ndc
!
! deposit electron current and acceleration density and momentum flux
! with OpenMP: updates cue, dcu, amu
      call dtimer(dtime,itime,-1)
      cue = 0.0; dcu = 0.0; amu = 0.0
      call GDCJPPOST1L(ppart,exyze,byze,cue,dcu,amu,kpic,omx,qme,qbme,dt&
     &,idimp,nppmx0,nx,mx,nxe,mx1)
! add scaled electric field with standard procedure: updates dcu
      call ASCFGUARD1L(dcu,cus,q2m0,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdcjpost = tdcjpost + time
!
! add guard cells for current, acceleration density, and momentum flux
! with standard procedure: updates cue, dcu, amu
      call dtimer(dtime,itime,-1)
      call ACGUARD1L(cue,nx,nxe)
      call ACGUARD1L(dcu,nx,nxe)
      call ACGUARD1L(amu,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform current to fourier space with standard procedure:
! update cue, gyze
      call dtimer(dtime,itime,-1)
      isign = -1
      call FFT1R2X(cue,gyze,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! calculate magnetic field in fourier space with standard procedure:
! updates byze, wm
      call dtimer(dtime,itime,-1)
      call BBPOIS13(cue,byze,ffc,ci,wm,nx,nxeh,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform magnetic force to real space with standard procedure:
! updates byze, gyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call FFT1R2X(byze,gyze,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! add constant to magnetic field with standard procedure: updates bzye
      call dtimer(dtime,itime,-1)
      call BADDEXT1(byze,omy,omz,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform acceleration density and momentum flux to fourier space
! with standard procedure: updates dcu, amu, gyze
      call dtimer(dtime,itime,-1)
      isign = -1
      call FFT1R2X(dcu,gyze,isign,mixup,sct,indx,nxe,nxh)
      call FFT1R2X(amu,gyze,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of time derivative of current with standard
! procedure: updates dcu
      call dtimer(dtime,itime,-1)
      call ADCUPERP13(dcu,amu,nx,nxeh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate transverse electric field with standard procedure:
! updates cus, wf
      call dtimer(dtime,itime,-1)
      isign = -1
      call EPOIS13(dcu,cus,isign,ffe,ax,affp,wpm,ci,wf,nx,nxeh,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform transverse electric field to real space with standard
! procedure: updates cus, gyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call FFT1R2X(cus,gyze,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with standard procedure: updates byze, cus
      call dtimer(dtime,itime,-1)
      call CGUARD1L(byze,nx,nxe)
      call CGUARD1L(cus,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and transverse electric fields with standard
! procedure: exyze = cus + fxe, updates exyze
! cus needs to be retained for next time step
      call dtimer(dtime,itime,-1)
      call ADDVRFIELD13(exyze,cus,fxe,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
      enddo
!
! push particles with OpenMP:
      wke = 0.0
      call dtimer(dtime,itime,-1)
! updates ppart, wke
!     call GBPPUSH13L(ppart,exyze,byze,kpic,omx,qbme,dt,dt,wke,idimp,   &
!    &nppmx0,nx,mx,nxe,mx1,ipbc)
! updates ppart, ncl, ihole, wke, irc
      call GBPPUSHF13L(ppart,exyze,byze,kpic,ncl,ihole,omx,qbme,dt,dt,  &
     &wke,idimp,nppmx0,nx,mx,nxe,mx1,ntmax,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
      if (irc /= 0) then
         write (*,*) 'GBPPUSHF13L error: irc=', irc
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
         wt = we + wm
         write (*,*) 'Initial Total Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') wt, wke, wke + wt
         write (*,*) 'Initial Electrostatic, Transverse Electric and Mag&
     &netic Field Energies:'
         write (*,'(3e14.7)') we, wf, wm
      endif
      ntime = ntime + 1
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
      write (*,*) 'ntime, ndc = ', ntime, ndc
      wt = we + wm
      write (*,*) 'Final Total Field, Kinetic and Total Energies:'
      write (*,'(3e14.7)') wt, wke, wke + wt
      write (*,*) 'Final Electrostatic, Transverse Electric and Magnetic&
     & Field Energies:'
      write (*,'(3e14.7)') we, wf, wm
!
      write (*,*)
      write (*,*) 'deposit time = ', tdpost
      write (*,*) 'current deposit time = ', tdjpost
      write (*,*) 'current derivative deposit time = ', tdcjpost
      tdpost = tdpost + tdjpost + tdcjpost
      write (*,*) 'total deposit time = ', tdpost
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
!
      stop
      end program
