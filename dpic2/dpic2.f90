!-----------------------------------------------------------------------
! Skeleton 2-1/2D Darwin PIC code
! written by Viktor K. Decyk, UCLA
      program dpic2
      use dpush2_h
      implicit none
! indx/indy = exponent which determines grid points in x/y direction:
! nx = 2**indx, ny = 2**indy.
      integer, parameter :: indx =   9, indy =   9
! npx/npy = number of electrons distributed in x/y direction.
      integer, parameter :: npx =  3072, npy =   3072
! ndim = number of velocity coordinates = 3
      integer, parameter :: ndim = 3
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
! vtx/vty = thermal velocity of electrons in x/y direction
! vx0/vy0 = drift velocity of electrons in x/y direction.
      real, parameter :: vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0
! vtx/vz0 = thermal/drift velocity of electrons in z direction
      real, parameter :: vtz = 1.0, vz0 = 0.0
! ax/ay = smoothed particle size in x/y direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ay = .912871, ci = 0.1
! idimp = number of particle coordinates = 5
! ipbc = particle boundary condition: 1 = periodic
! sortime = number of time steps between standard electron sorting
      integer :: idimp = 5, ipbc = 1, sortime = 50
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
      integer :: np, nx, ny, nxh, nyh, nxe, nye, nxeh, nxyh, nxhy
      integer :: mdim, ny1, ntime, nloop, isign
      real :: qbme, affp, q2m0, wpm, wpmax, wpmin
!
! declare arrays for standard code:
! part, part2 = particle arrays
      real, dimension(:,:), pointer :: part, part2, tpart
! qe = electron charge density with guard cells
      real, dimension(:,:), pointer :: qe
! cue = electron current density with guard cells
! dcu = acceleration density with guard cells
! cus = smoothed transverse electric field with guard cells
! amu = momentum flux with guard cells
      real, dimension(:,:,:), pointer :: cue, dcu, cus, amu
! exyze = smoothed total electric field with guard cells
! fxyze = smoothed longitudinal electric field with guard cells
! bxyze = smoothed magnetic field with guard cells
      real, dimension(:,:,:), pointer :: fxyze, exyze, bxyze
! ffc, ffe = form factor arrays for poisson solvers
      complex, dimension(:,:), pointer :: ffc, ffe
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct
! npicy = scratch array for reordering particles
      integer, dimension(:), pointer :: npicy
! ss = scratch array for WFFT2RN
      complex, dimension(:,:), pointer :: ss
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
! nx/ny = number of grid points in x/y direction
      np = npx*npy; nx = 2**indx; ny = 2**indy
      nxh = nx/2; nyh = max(1,ny/2)
      nxe = nx + 2; nye = ny + 1; nxeh = nxe/2
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny); ny1 = ny + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
! mdim = dimension of amu array
      mdim = 2*ndim - 2
      qbme = qme
      affp = real(nx*ny)/real(np)
!
! allocate data for standard code
      allocate(part(idimp,np))
      if (sortime > 0) allocate(part2(idimp,np))
      allocate(qe(nxe,nye),fxyze(ndim,nxe,nye),exyze(ndim,nxe,nye))
      allocate(cue(ndim,nxe,nye),dcu(ndim,nxe,nye),cus(ndim,nxe,nye))
      allocate(amu(mdim,nxe,nye),bxyze(ndim,nxe,nye))
      allocate(ffc(nxh,nyh),ffe(nxh,nyh),mixup(nxhy),sct(nxyh))
      allocate(npicy(ny1),ss(mdim,nxeh))
!
! prepare fft tables
      call WFFT2RINIT(mixup,sct,indx,indy,nxhy,nxyh)
! calculate form factor: ffc
      isign = 0
      call POIS23(qe,fxyze,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,  &
     &nyh)
! initialize electrons
      call DISTR2H(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,idimp,np,nx,ny, &
     &ipbc)
!
! find maximum and minimum initial electron density
      qe = 0.0
      call GPOST2L(part,qe,qme,np,idimp,nxe,nye)
      call AGUARD2L(qe,nx,ny,nxe,nye)
      call FWPMINMX2(qe,qbme,wpmax,wpmin,nx,ny,nxe,nye)
      wpm = 0.5*(wpmax + wpmin)*affp
! accelerate convergence: update wpm
      if (wpm <= 10.0) wpm = 0.75*wpm
      write (*,*) 'wpm=',wpm
      q2m0 = wpm/affp
! calculate form factor: ffe
      isign = 0
      call EPOIS23(dcu,cus,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,nxeh,nye&
     &,nxh,nyh)
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
      call GJPOST2L(part,cue,qme,zero,np,idimp,nx,ny,nxe,nye,ipbc)
!     call GSJPOST2L(part,cue,qme,zero,np,idimp,nx,ny,nxe,nxe*nye,ipbc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdjpost = tdjpost + time
!
! deposit charge with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call GPOST2L(part,qe,qme,np,idimp,nxe,nye)
!     call GSPOST2L(part,qe,qme,np,idimp,nxe,nxe*nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with standard procedure: updates qe, cue
      call dtimer(dtime,itime,-1)
      call AGUARD2L(qe,nx,ny,nxe,nye)
      call ACGUARD2L(cue,nx,ny,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      isign = -1
      call WFFT2RX(qe,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! calculate longitudinal force/charge in fourier space with standard
! procedure: updates fxyze, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call POIS23(qe,fxyze,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,  &
     &nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform longitudinal electric force to real space with standard
! procedure: updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT2R3(fxyze,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! transform current to fourier space with standard procedure: update cue
      call dtimer(dtime,itime,-1)
      isign = -1
      call WFFT2R3(cue,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of current with standard procedure: updates cue
      call dtimer(dtime,itime,-1)
      call CUPERP2(cue,nx,ny,nxeh,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate magnetic field in fourier space with standard procedure:
! updates bxyze, wm
      call dtimer(dtime,itime,-1)
      call BBPOIS23(cue,bxyze,ffc,ci,wm,nx,ny,nxeh,nye,nxh,nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform magnetic force to real space with standard procedure:
! updates bxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT2R3(bxyze,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! add constant to magnetic field with standard procedure: updates bxyze
      call dtimer(dtime,itime,-1)
      call BADDEXT2(bxyze,omx,omy,omz,nx,ny,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! copy guard cells with standard procedure: updates fxyze, bxyze
      call dtimer(dtime,itime,-1)
      call BGUARD2L(fxyze,nx,ny,nxe,nye)
      call BGUARD2L(bxyze,nx,ny,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and old transverse electric fields with standard
! procedure: updates exyze
      call dtimer(dtime,itime,-1)
      call ADDVRFIELD2(exyze,cus,fxyze,ndim,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! deposit electron acceleration density and momentum flux with standard
! procedure: updates dcu, amu
      call dtimer(dtime,itime,-1)
      dcu = 0.0; amu = 0.0
      call GDJPOST2L(part,exyze,bxyze,dcu,amu,qme,qbme,dt,idimp,np,nxe, &
     &nye)
!     call GSDJPOST2L(part,exyze,bxyze,dcu,amu,qme,qbme,dt,idimp,np,nxe,&
!    &nxe*nye)
! add old scaled electric field with standard procedure: updates dcu
      call ASCFGUARD2L(dcu,cus,q2m0,nx,ny,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdcjpost = tdcjpost + time
!
! add guard cells with standard procedure: updates dcu, amu
      call dtimer(dtime,itime,-1)
      call ACGUARD2L(dcu,nx,ny,nxe,nye)
      call AMCGUARD2L(amu,nx,ny,nxe,nye,mdim)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform acceleration density and momentum flux to fourier space
! with standard procedure: updates dcu, amu
      call dtimer(dtime,itime,-1)
      isign = -1
      call WFFT2R3(dcu,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call WFFT2RN(amu,ss,isign,mixup,sct,indx,indy,nxeh,nye,mdim,nxhy, &
     &nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of time derivative of current with standard
! procedure: updates dcu
      call dtimer(dtime,itime,-1)
      call ADCUPERP23(dcu,amu,nx,ny,nxeh,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate transverse electric field with standard procedure:
! updates cus, wf
      call dtimer(dtime,itime,-1)
      isign = -1
      call EPOIS23(dcu,cus,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,nxeh,nye&
     &,nxh,nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform transverse electric field to real space with standard
! procedure: updates cus
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT2R3(cus,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with standard procedure: updates cus
      call dtimer(dtime,itime,-1)
      call BGUARD2L(cus,nx,ny,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and transverse electric fields with standard
! procedure: exyze = cus + fxyze, updates exyze
! cus needs to be retained for next time step
      call dtimer(dtime,itime,-1)
      call ADDVRFIELD2(exyze,cus,fxyze,ndim,nxe,nye)
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
      call GDCJPOST2L(part,exyze,bxyze,cue,dcu,amu,qme,qbme,dt,idimp,np,&
     &nxe,nye)
!     call GSDCJPOST2L(part,exyze,bxyze,cue,dcu,amu,qme,qbme,dt,idimp,np&
!    &,nxe,nxe*nye)
! add scaled electric field with standard procedure: updates dcu
      call ASCFGUARD2L(dcu,cus,q2m0,nx,ny,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdcjpost = tdcjpost + time
!
! add guard cells for current, acceleration density, and momentum flux
! with standard procedure: updates cue, dcu, amu
      call dtimer(dtime,itime,-1)
      call ACGUARD2L(cue,nx,ny,nxe,nye)
      call ACGUARD2L(dcu,nx,ny,nxe,nye)
      call AMCGUARD2L(amu,nx,ny,nxe,nye,mdim)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform current to fourier space with standard procedure: update cue
      call dtimer(dtime,itime,-1)
      isign = -1
      call WFFT2R3(cue,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of current with standard procedure: updates cue
      call dtimer(dtime,itime,-1)
      call CUPERP2(cue,nx,ny,nxeh,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate magnetic field in fourier space with standard procedure:
! updates bxyze, wm
      call dtimer(dtime,itime,-1)
      call BBPOIS23(cue,bxyze,ffc,ci,wm,nx,ny,nxeh,nye,nxh,nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform magnetic force to real space with standard procedure:
! updates bxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT2R3(bxyze,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! add constant to magnetic field with standard procedure: updates bxzye
      call dtimer(dtime,itime,-1)
      call BADDEXT2(bxyze,omx,omy,omz,nx,ny,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform acceleration density and momentum flux to fourier space
! with standard procedure: updates dcu and amu
      call dtimer(dtime,itime,-1)
      isign = -1
      call WFFT2R3(dcu,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call WFFT2RN(amu,ss,isign,mixup,sct,indx,indy,nxeh,nye,mdim,nxhy, &
     &nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of time derivative of current with standard
! procedure: updates dcu
      call dtimer(dtime,itime,-1)
      call ADCUPERP23(dcu,amu,nx,ny,nxeh,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate transverse electric field with standard procedure:
! updates cus, wf
      call dtimer(dtime,itime,-1)
      isign = -1
      call EPOIS23(dcu,cus,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,nxeh,nye&
     &,nxh,nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform transverse electric field to real space with standard
! procedure: updates cus
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT2R3(cus,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with standard procedure: updates bxyze, cus
      call dtimer(dtime,itime,-1)
      call BGUARD2L(bxyze,nx,ny,nxe,nye)
      call BGUARD2L(cus,nx,ny,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and transverse electric fields with standard
! procedure: exyze = cus + fxyze, updates exyze
! cus needs to be retained for next time step
      call dtimer(dtime,itime,-1)
      call ADDVRFIELD2(exyze,cus,fxyze,ndim,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
      enddo
!
! push particles with standard procedure: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
      call GBPUSH23L(part,exyze,bxyze,qbme,dt,dt,wke,idimp,np,nx,ny,nxe,&
     &nye,ipbc)
!     call GSBPUSH23L(part,exyze,bxyze,qbme,dt,dt,wke,idimp,np,nx,ny,nxe&
!    &,nxe*nye,ipbc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
!
! sort particles by cell for standard procedure
      if (sortime > 0) then
         if (mod(ntime,sortime)==0) then
            call dtimer(dtime,itime,-1)
            call DSORTP2YL(part,part2,npicy,idimp,np,ny1)
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
!
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
