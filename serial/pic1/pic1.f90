!-----------------------------------------------------------------------
! Skeleton 1D Electrostatic PIC code
! written by Viktor K. Decyk, UCLA
      program pic1
      use push1_h
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
! sortime = number of time steps between standard electron sorting
      integer :: idimp = 2, ipbc = 1, sortime = 50
! wke/we/wt = particle kinetic/electric field/total energy
      real :: wke = 0.0, we = 0.0, wt = 0.0
! declare scalars for standard code
      integer :: np, nx, nxh, nxe
      integer :: nx1, ntime, nloop, isign
      real :: qbme, affp
!
! declare arrays for standard code:
! part, part2 = particle arrays
      real, dimension(:,:), pointer :: part, part2, tpart
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
! npic = scratch array for reordering particles
      integer, dimension(:), pointer :: npic
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
      real :: tpush = 0.0, tsort = 0.0
      double precision :: dtime
!
! initialize scalars for standard code
! np = total number of particles in simulation
! nx = number of grid points in x direction
      np = npx; nx = 2**indx; nxh = nx/2
      nxe = nx + 2; nx1 = nx + 1
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
      allocate(ffc(nxh),mixup(nxh),sct(nxh))
      allocate(npic(nx1))
!
! prepare fft tables
      call WFFT1RINIT(mixup,sct,indx,nxh)
! calculate form factors
      isign = 0
      call POIS1(qe,fxe,isign,ffc,ax,affp,we,nx)
! initialize electrons
      call DISTR1(part,vtx,vx0,npx,idimp,np,nx,ipbc)
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     write (*,*) 'ntime = ', ntime
!
! deposit charge with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call GPOST1L(part,qe,qme,np,idimp,nxe)
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
! push particles with standard procedure: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
      call GPUSH1L(part,fxe,qbme,dt,wke,idimp,np,nx,nxe,ipbc)
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
!
      stop
      end program
