!-----------------------------------------------------------------------
! Skeleton 1-2/2D Electromagnetic PIC code
! written by Viktor K. Decyk, UCLA
      program bpic1
      use bpush1_h
      implicit none
! indx = exponent which determines grid points in x direction:
! nx = 2**indx.
      integer, parameter :: indx =   9
! npx = number of electrons distributed in x direction.
      integer, parameter :: npx =  18432
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.05, qme = -1.0
! vtx/vty = thermal velocity of electrons in x/y direction
! vx0/vy0 = drift velocity of electrons in x/y direction.
      real, parameter :: vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0
! vtx/vz0 = thermal/drift velocity of electrons in z direction
      real, parameter :: vtz = 1.0, vz0 = 0.0
! omx = magnetic field electron cyclotron frequency in x
      real :: omx = 0.0
! ax = smoothed particle size in x direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ci = 0.1
! idimp = number of particle coordinates = 4
! ipbc = particle boundary condition: 1 = periodic
! sortime = number of time steps between standard electron sorting
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 4, ipbc = 1, sortime = 50, relativity = 1
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
! declare scalars for standard code
      integer :: np, nx, nxh, nxe, nxeh
      integer :: nx1, ntime, nloop, isign
      real :: qbme, affp, dth
!
! declare arrays for standard code:
! part, part2 = particle arrays
      real, dimension(:,:), pointer :: part, part2, tpart
! qe = electron charge density with guard cells
! fxe = smoothed longitudinal electric field with guard cells
      real, dimension(:), pointer :: qe, fxe
! cue = electron current density with guard cells
! fxyze/byze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:), pointer :: cue, fxyze, byze
! eyz/byz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:), pointer :: eyz, byz
! ffc = form factor array for poisson solver
      complex, dimension(:), pointer :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct
! npic = scratch array for reordering particles
      integer, dimension(:), pointer :: npic
! gxyze = scratch array for fft
      real, dimension(:,:), pointer :: gxyze
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
! nx = number of grid points in x direction
      np = npx; nx = 2**indx; nxh = nx/2
      nxe = nx + 2; nxeh = nxe/2; nx1 = nx + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx)/real(np)
      dth = 0.0
!
! allocate data for standard code
      allocate(part(idimp,np))
      if (sortime > 0) allocate(part2(idimp,np))
      allocate(qe(nxe),fxe(nxe))
      allocate(fxyze(3,nxe),cue(2,nxe),byze(2,nxe))
      allocate(eyz(2,nxeh),byz(2,nxeh))
      allocate(ffc(nxh),mixup(nxh),sct(nxh))
      allocate(npic(nx1))
      allocate(gxyze(3,nxe))
!
! prepare fft tables
      call WFFT1RINIT(mixup,sct,indx,nxh)
! calculate form factors
      isign = 0
      call POIS1(qe,fxe,isign,ffc,ax,affp,we,nx)
! initialize electrons
      call DISTR1H(part,vtx,vty,vtz,vx0,vy0,vz0,npx,idimp,np,nx,ipbc)
!
! initialize transverse electromagnetic fields
      eyz = cmplx(0.0,0.0)
      byz = cmplx(0.0,0.0)
!
      if (dt > 0.64*ci) then
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
         call GRJPOST1L(part,cue,qme,dth,ci,np,idimp,nx,nxe,ipbc)
      else
         call GJPOST1L(part,cue,qme,dth,np,idimp,nx,nxe,ipbc)
      endif
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
! add guard cells with standard procedure: updates cue, qe
      call dtimer(dtime,itime,-1)
      call ACGUARD1L(cue,nx,nxe)
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
! transform current to fourier space with standard procedure:
! updates cue, byze
      call dtimer(dtime,itime,-1)
      isign = -1
      call FFT1R2X(cue,byze,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! calculate electromagnetic fields in fourier space with standard
! procedure: updates eyz, byz
      call dtimer(dtime,itime,-1)
      if (ntime==0) then
         call IBPOIS13(cue,byz,ffc,ci,wm,nx,nxeh,nxh)
         wf = 0.0
         dth = 0.5*dt
      else
         call MAXWEL1(eyz,byz,cue,ffc,ci,dt,wf,wm,nx,nxeh,nxh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
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
! add longitudinal and transverse electric fields with standard
! procedure: updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call EMFIELD1(fxyze,fxe,eyz,ffc,nx,nxeh,nxh)
! copy magnetic field with standard procedure: updates byze
      isign = -1
      call BMFIELD1(byze,byz,ffc,nx,nxeh,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform electric force to real space with standard procedure:
! updates fxyze, gxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call FFT1R3X(fxyze,gxyze,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! transform magnetic force to real space with standard procedure:
! updates byze, cue
      call dtimer(dtime,itime,-1)
      isign = 1
      call FFT1R2X(byze,cue,isign,mixup,sct,indx,nxe,nxh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with standard procedure: updates fxyze, byze
      call dtimer(dtime,itime,-1)
      call BGUARD1L(fxyze,nx,nxe)
      call CGUARD1L(byze,nx,nxe)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! push particles with standard procedure: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
      if (relativity==1) then
         call GRBPUSH13L(part,fxyze,byze,omx,qbme,dt,dth,ci,wke,idimp,np&
     &,nx,nxe,ipbc)
      else
         call GBPUSH13L(part,fxyze,byze,omx,qbme,dt,dth,wke,idimp,np,nx,&
     &nxe,ipbc)
      endif
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
