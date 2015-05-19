!-----------------------------------------------------------------------
! Skeleton 1-2/2D Darwin PIC code
! written by Viktor K. Decyk, UCLA
      program dpic1
      use dpush1_h
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
! sortime = number of time steps between standard electron sorting
      integer :: idimp = 4, ipbc = 1, sortime = 50
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.4, omy = 0.0, omz = 0.0
! ndc = number of corrections in darwin iteration
      integer :: ndc = 1
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
      real :: zero = 0.0
! declare scalars for standard code
      integer :: k
      integer :: np, nx, nxh, nxe, nxeh
      integer :: nx1, ntime, nloop, isign
      real :: qbme, affp, q2m0, wpm, wpmax, wpmin
!
! declare arrays for standard code:
! part, part2 = particle arrays
      real, dimension(:,:), pointer :: part, part2, tpart
! qe = electron charge density with guard cells
! fxe = smoothed longitudinal electric field with guard cells
      real, dimension(:), pointer :: qe, fxe
! cue = electron current density with guard cells
! dcu = acceleration density with guard cells
! cus = transverse electric field with guard cells
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
! npic = scratch array for reordering particles
      integer, dimension(:), pointer :: npic
! gxe, gyze = scratch arrays for fft
      real, dimension(:), pointer :: gxe
      real, dimension(:,:), pointer :: gyze
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
      real :: tdjpost = 0.0, tdcjpost = 0.0, tpush = 0.0, tsort = 0.0
      double precision :: dtime
!
! initialize scalars for standard code
! np = total number of particles in simulation
! nx = number of grid points in x direction
      np = npx; nx = 2**indx; nxh = nx/2
      nxe = nx + 2; nxeh = nxe/2; nx1 = nx + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx)/real(np)
!
! allocate data for standard code
      allocate(part(idimp,np))
      if (sortime > 0) allocate(part2(idimp,np))
      allocate(qe(nxe),fxe(nxe))
      allocate(cue(2,nxe),dcu(2,nxe),cus(2,nxe),amu(2,nxe))
      allocate(exyze(3,nxe),byze(2,nxe))
      allocate(ffc(nxh),ffe(nxh),mixup(nxh),sct(nxh))
      allocate(npic(nx1))
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
! find maximum and minimum initial electron density
      qe = 0.0
      call GPOST1L(part,qe,qme,np,idimp,nxe)
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
! deposit current with standard procedure: updates cue
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call GJPOST1L(part,cue,qme,zero,np,idimp,nx,nxe,ipbc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdjpost = tdjpost + time
!
! deposit charge with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call GPOST1L(part,qe,qme,np,idimp,nxe)
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
! deposit electron acceleration density and momentum flux with standard
! procedure: updates dcu, amu
      call dtimer(dtime,itime,-1)
      dcu = 0.0; amu = 0.0
      call GDJPOST1L(part,exyze,byze,dcu,amu,omx,qme,qbme,dt,idimp,np,  &
     &nxe)
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
! with standard procedure: updates cue, dcu, amu
      call dtimer(dtime,itime,-1)
      cue = 0.0; dcu = 0.0; amu = 0.0
      call GDCJPOST1L(part,exyze,byze,cue,dcu,amu,omx,qme,qbme,dt,idimp,&
     &np,nxe)
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
! procedure: exyze = cus + fxyze, updates exyze
! cus needs to be retained for next time step
      call dtimer(dtime,itime,-1)
      call ADDVRFIELD13(exyze,cus,fxe,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
      enddo
!
! push particles with standard procedure: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
      call GBPUSH13L(part,exyze,byze,omx,qbme,dt,dt,wke,idimp,np,nx,nxe,&
     &ipbc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
!
! sort particles by cell for standard procedure
      if (sortime > 0) then
         if (mod(ntime,sortime)==0) then
            call dtimer(dtime,itime,-1)
            call DSORTP1XL(part,part2,npic,idimp,np,nx1)
! exchange pointers
            tpart => part
            part => part2
            part2 => tpart
            call dtimer(dtime,itime,1)
            time = real(dtime)
            tsort = tsort + time
         endif
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
