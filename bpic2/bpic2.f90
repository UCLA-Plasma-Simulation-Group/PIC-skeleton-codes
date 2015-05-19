!-----------------------------------------------------------------------
! Skeleton 2-1/2D Electromagnetic PIC code
! written by Viktor K. Decyk, UCLA
      program bpic2
      use bpush2_h
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
      real, parameter :: tend = 10.0, dt = 0.04, qme = -1.0
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
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 5, ipbc = 1, sortime = 50, relativity = 1
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
! declare scalars for standard code
      integer :: np, nx, ny, nxh, nyh, nxe, nye, nxeh, nxyh, nxhy
      integer :: ny1, ntime, nloop, isign
      real :: qbme, affp, dth
!
! declare arrays for standard code:
! part, part2 = particle arrays
      real, dimension(:,:), pointer :: part, part2, tpart
! qe = electron charge density with guard cells
      real, dimension(:,:), pointer :: qe
! cue = electron current density with guard cells
! fxyze/bxyze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:,:), pointer :: cue, fxyze, bxyze
! exyz/bxyz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:,:), pointer :: exyz, bxyz
! ffc = form factor array for poisson solver
      complex, dimension(:,:), pointer :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct
! npicy = scratch array for reordering particles
      integer, dimension(:), pointer :: npicy
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
      real :: tdjpost = 0.0, tpush = 0.0, tsort = 0.0
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
      qbme = qme
      affp = real(nx*ny)/real(np)
      dth = 0.0
!
! allocate data for standard code
      allocate(part(idimp,np))
      if (sortime > 0) allocate(part2(idimp,np))
      allocate(qe(nxe,nye),fxyze(ndim,nxe,nye))
      allocate(cue(ndim,nxe,nye),bxyze(ndim,nxe,nye))
      allocate(exyz(ndim,nxeh,nye),bxyz(ndim,nxeh,nye))
      allocate(ffc(nxh,nyh),mixup(nxhy),sct(nxyh))
      allocate(npicy(ny1))
!
! prepare fft tables
      call WFFT2RINIT(mixup,sct,indx,indy,nxhy,nxyh)
! calculate form factors
      isign = 0
      call POIS23(qe,fxyze,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,  &
     &nyh)
! initialize electrons
      call DISTR2H(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,idimp,np,nx,ny, &
     &ipbc)
!
! initialize transverse electromagnetic fields
      exyz = cmplx(0.0,0.0)
      bxyz = cmplx(0.0,0.0)
!
      if (dt > 0.45*ci) then
         write (*,*) 'Warning: Courant condition may be exceeded!'
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     write (*,*) 'ntime = ', ntime
!
! deposit current with standard procedure: updates part, cue
      call dtimer(dtime,itime,-1)
      cue = 0.0
      if (relativity==1) then
         call GRJPOST2L(part,cue,qme,dth,ci,np,idimp,nx,ny,nxe,nye,ipbc)
!        call GSRJPOST2L(part,cue,qme,dth,ci,np,idimp,nx,ny,nxe,nxe*nye,&
!    &ipbc)
      else
         call GJPOST2L(part,cue,qme,dth,np,idimp,nx,ny,nxe,nye,ipbc)
!        call GSJPOST2L(part,cue,qme,dth,np,idimp,nx,ny,nxe,nxe*nye,ipbc&
!    &)
      endif
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
! add guard cells with standard procedure: updates cue, qe
      call dtimer(dtime,itime,-1)
      call ACGUARD2L(cue,nx,ny,nxe,nye)
      call AGUARD2L(qe,nx,ny,nxe,nye)
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
! calculate electromagnetic fields in fourier space with standard
! procedure: updates exyz, bxyz, wf, wm
      call dtimer(dtime,itime,-1)
      if (ntime==0) then
         call IBPOIS23(cue,bxyz,ffc,ci,wm,nx,ny,nxeh,nye,nxh,nyh)
         wf = 0.0
         dth = 0.5*dt
      else
         call MAXWEL2(exyz,bxyz,cue,ffc,ci,dt,wf,wm,nx,ny,nxeh,nye,nxh, &
     &nyh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate force/charge in fourier space with standard procedure:
! updates fxyze, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call POIS23(qe,fxyze,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,  &
     &nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! add longitudinal and transverse electric fields with standard
! procedure: updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call EMFIELD2(fxyze,exyz,ffc,isign,nx,ny,nxeh,nye,nxh,nyh)
! copy magnetic field with standard procedure: updates bxyze
      isign = -1
      call EMFIELD2(bxyze,bxyz,ffc,isign,nx,ny,nxeh,nye,nxh,nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform electric force to real space with standard procedure:
! updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT2R3(fxyze,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
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
! copy guard cells with standard procedure: updates fxyze, bxyze
      call dtimer(dtime,itime,-1)
      call BGUARD2L(fxyze,nx,ny,nxe,nye)
      call BGUARD2L(bxyze,nx,ny,nxe,nye)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! push particles with standard procedure: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
      if (relativity==1) then
         call GRBPUSH23L(part,fxyze,bxyze,qbme,dt,dth,ci,wke,idimp,np,nx&
     &,ny,nxe,nye,ipbc)
!        call GSRBPUSH23L(part,fxyze,bxyze,qbme,dt,dth,ci,wke,idimp,np, &
!    &nx,ny,nxe,nxe*nye,ipbc)
      else
         call GBPUSH23L(part,fxyze,bxyze,qbme,dt,dth,wke,idimp,np,nx,ny,&
     &nxe,nye,ipbc)
!        call GSBPUSH23L(part,fxyze,bxyze,qbme,dt,dth,wke,idimp,np,nx,ny&
!    &,nxe,nxe*nye,ipbc)
      endif
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
         wt = we + wf + wm
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
      write (*,*) 'ntime, relativity = ', ntime, relativity
      wt = we + wf + wm
      write (*,*) 'Final Total Field, Kinetic and Total Energies:'
      write (*,'(3e14.7)') wt, wke, wke + wt
      write (*,*) 'Final Electrostatic, Transverse Electric and Magnetic&
     & Field Energies:'
      write (*,'(3e14.7)') we, wf, wm
!
      write (*,*)
      write (*,*) 'deposit time = ', tdpost
      write (*,*) 'current deposit time = ', tdjpost
      tdpost = tdpost + tdjpost
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
