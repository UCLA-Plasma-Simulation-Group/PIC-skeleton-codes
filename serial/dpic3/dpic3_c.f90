!-----------------------------------------------------------------------
! Skeleton 3D Darwin PIC code
! written by Viktor K. Decyk, UCLA
      program dpic3
! #include "dpush3.h"
      implicit none
! indx/indy/indz = exponent which determines grid points in x/y/z
! direction: nx = 2**indx, ny = 2**indy, nz = 2**indz.
      integer, parameter :: indx =   7, indy =   7, indz =   7
! npx/npy/npz = number of electrons distributed in x/y/z direction.
      integer, parameter :: npx =  384, npy =   384, npz =   384
! ndim = number of velocity coordinates = 3
      integer, parameter :: ndim = 3
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real, parameter :: vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real, parameter :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
! ax/ay/az = smoothed particle size in x/y/z direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ay = .912871, az = .912871, ci = 0.1
! idimp = number of particle coordinates = 6
! ipbc = particle boundary condition: 1 = periodic
! sortime = number of time steps between standard electron sorting
      integer :: idimp = 6, ipbc = 1, sortime = 20
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
      integer :: np, nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh
      integer :: mdim, nxyzh, nxhyz, ny1, nyz1, ntime, nloop, isign
      real :: qbme, affp, q2m0, wpm, wpmax, wpmin
!
! declare arrays for standard code:
! part, part2 = particle arrays
      real, dimension(:,:), pointer :: part, part2, tpart
! qe = electron charge density with guard cells
      real, dimension(:,:,:), pointer :: qe
! cue = electron current density with guard cells
! dcu = acceleration density with guard cells
! cus = smoothed transverse electric field with guard cells
! amu = momentum flux with guard cells
      real, dimension(:,:,:,:), pointer :: cue, dcu, cus, amu
! exyze = smoothed total electric field with guard cells
! fxyze = smoothed longitudinal electric field with guard cells
! bxyze = smoothed magnetic field with guard cells
      real, dimension(:,:,:,:), pointer :: fxyze, exyze, bxyze
! ffc, ffe = form factor arrays for poisson solvers
      complex, dimension(:,:,:), pointer :: ffc, ffe
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct
! npic = scratch array for reordering particles
      integer, dimension(:), pointer :: npic
! ss = scratch array for WFFT3RN
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
! nx/ny/nz = number of grid points in x/y/z direction
      np = npx*npy*npz; nx = 2**indx; ny = 2**indy; nz = 2**indz
      nxh = nx/2; nyh = max(1,ny/2); nzh = max(1,nz/2)
      nxe = nx + 2; nye = ny + 1; nze = nz + 1; nxeh = nxe/2
      nxyzh = max(nx,ny,nz)/2; nxhyz = max(nxh,ny,nz)
      ny1 = ny + 1; nyz1 = ny1*(nz + 1)
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
! mdim = dimension of amu array
      mdim = 2*ndim
      qbme = qme
      affp = real(nx)*real(ny)*real(nz)/real(np)
!
! allocate data for standard code
      allocate(part(idimp,np))
      if (sortime > 0) allocate(part2(idimp,np))
      allocate(qe(nxe,nye,nze),fxyze(ndim,nxe,nye,nze))
      allocate(cue(ndim,nxe,nye,nze),dcu(ndim,nxe,nye,nze))
      allocate(cus(ndim,nxe,nye,nze),amu(mdim,nxe,nye,nze))
      allocate(exyze(ndim,nxe,nye,nze),bxyze(ndim,nxe,nye,nze))
      allocate(ffc(nxh,nyh,nzh),ffe(nxh,nyh,nzh))
      allocate(mixup(nxhyz),sct(nxyzh))
      allocate(npic(nyz1),ss(mdim,nxeh))
!
! prepare fft tables
      call CWFFT3RINIT(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
! calculate form factor: ffc
      isign = 0
      call CPOIS33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,nye&
     &,nze,nxh,nyh,nzh)
! initialize electrons
      call CDISTR3(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,idimp,np,nx,&
     &ny,nz,ipbc)
!
! find maximum and minimum initial electron density
      qe = 0.0
      call CGPOST3L(part,qe,qme,np,idimp,nxe,nye,nze)
      call CAGUARD3L(qe,nx,ny,nz,nxe,nye,nze)
      call CFWPMINMX3(qe,qbme,wpmax,wpmin,nx,ny,nz,nxe,nye,nze)
      wpm = 0.5*(wpmax + wpmin)*affp
! accelerate convergence: update wpm
      if (wpm <= 10.0) wpm = 0.75*wpm
      write (*,*) 'wpm=',wpm
      q2m0 = wpm/affp
! calculate form factor: ffe
      isign = 0
      call CEPOIS33(dcu,cus,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,nz, &
     &nxeh,nye,nze,nxh,nyh,nzh)
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
      call CGJPOST3L(part,cue,qme,zero,np,idimp,nx,ny,nz,nxe,nye,nze,   &
     &ipbc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdjpost = tdjpost + time
!
! deposit charge with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call CGPOST3L(part,qe,qme,np,idimp,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with standard procedure: updates qe, cue
      call dtimer(dtime,itime,-1)
      call CAGUARD3L(qe,nx,ny,nz,nxe,nye,nze)
      call CACGUARD3L(cue,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      isign = -1
      call CWFFT3RX(qe,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz&
     &,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! calculate longitudinal force/charge in fourier space with standard
! procedure: updates fxyze, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call CPOIS33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,nye&
     &,nze,nxh,nyh,nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform longitudinal electric force to real space with standard
! procedure: updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call CWFFT3R3(fxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,  &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! transform current to fourier space with standard procedure: update cue
      call dtimer(dtime,itime,-1)
      isign = -1
      call CWFFT3R3(cue,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of current with standard procedure: updates cue
      call dtimer(dtime,itime,-1)
      call CCUPERP3(cue,nx,ny,nz,nxeh,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate magnetic field in fourier space with standard procedure:
! updates bxyze, wm
      call dtimer(dtime,itime,-1)
      call CBBPOIS33(cue,bxyze,ffc,ci,wm,nx,ny,nz,nxeh,nye,nze,nxh,nyh, &
     &nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform magnetic force to real space with standard procedure:
! updates bxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call CWFFT3R3(bxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,  &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! add constant to magnetic field with standard procedure: updates bxyze
      call dtimer(dtime,itime,-1)
      call CBADDEXT3(bxyze,omx,omy,omz,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! copy guard cells with standard procedure: updates fxyze, bxyze
      call dtimer(dtime,itime,-1)
      call CCGUARD3L(fxyze,nx,ny,nz,nxe,nye,nze)
      call CCGUARD3L(bxyze,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and old transverse electric fields with standard
! procedure: updates exyze
      call dtimer(dtime,itime,-1)
      call CADDVRFIELD3(exyze,cus,fxyze,ndim,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! deposit electron acceleration density and momentum flux with standard
! procedure: updates dcu, amu
      call dtimer(dtime,itime,-1)
      dcu = 0.0; amu = 0.0
      call CGDJPOST3L(part,exyze,bxyze,dcu,amu,qme,qbme,dt,idimp,np,nxe,&
     &nye,nze)
! add old scaled electric field with standard procedure: updates dcu
      call CASCFGUARD3L(dcu,cus,q2m0,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdcjpost = tdcjpost + time
!
! add guard cells with standard procedure: updates dcu, amu
      call dtimer(dtime,itime,-1)
      call CACGUARD3L(dcu,nx,ny,nz,nxe,nye,nze)
      call CAMCGUARD3L(amu,nx,ny,nz,nxe,nye,nze,mdim)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform acceleration density and momentum flux to fourier space
! with standard procedure: updates dcu, amu
      call dtimer(dtime,itime,-1)
      isign = -1
      call CWFFT3R3(dcu,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call CWFFT3RN(amu,ss,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze, &
     &mdim,nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of time derivative of current with standard
! procedure: updates dcu
      call dtimer(dtime,itime,-1)
      call CADCUPERP3(dcu,amu,nx,ny,nz,nxeh,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate transverse electric field with standard procedure:
! updates cus, wf
      call dtimer(dtime,itime,-1)
      isign = -1
      call CEPOIS33(dcu,cus,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,nz, &
     &nxeh,nye,nze,nxh,nyh,nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform transverse electric field to real space with standard
! procedure: updates cus
      call dtimer(dtime,itime,-1)
      isign = 1
      call CWFFT3R3(cus,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with standard procedure: updates cus
      call dtimer(dtime,itime,-1)
      call CCGUARD3L(cus,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and transverse electric fields with standard
! procedure: exyze = cus + fxyze, updates exyze
! cus needs to be retained for next time step
      call dtimer(dtime,itime,-1)
      call CADDVRFIELD3(exyze,cus,fxyze,ndim,nxe,nye,nze)
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
      call CGDCJPOST3L(part,exyze,bxyze,cue,dcu,amu,qme,qbme,dt,idimp,np&
     &,nxe,nye,nze)
! add caled electric field with standard procedure: updates dcu
      call CASCFGUARD3L(dcu,cus,q2m0,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdcjpost = tdcjpost + time
!
! add guard cells for current, acceleration density, and momentum flux
! with standard procedure: updates cue, dcu, amu
      call dtimer(dtime,itime,-1)
      call CACGUARD3L(cue,nx,ny,nz,nxe,nye,nze)
      call CACGUARD3L(dcu,nx,ny,nz,nxe,nye,nze)
      call CAMCGUARD3L(amu,nx,ny,nz,nxe,nye,nze,mdim)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform current to fourier space with standard procedure: update cue
      call dtimer(dtime,itime,-1)
      isign = -1
      call CWFFT3R3(cue,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of current with standard procedure: updates cue
      call dtimer(dtime,itime,-1)
      call CCUPERP3(cue,nx,ny,nz,nxeh,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate magnetic field in fourier space with standard procedure:
! updates bxyze, wm
      call dtimer(dtime,itime,-1)
      call CBBPOIS33(cue,bxyze,ffc,ci,wm,nx,ny,nz,nxeh,nye,nze,nxh,nyh, &
     &nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform magnetic force to real space with standard procedure:
! updates bxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call CWFFT3R3(bxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,  &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! add constant to magnetic field with standard procedure: updates bxyze
      call dtimer(dtime,itime,-1)
      call CBADDEXT3(bxyze,omx,omy,omz,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform acceleration density and momentum flux to fourier space
! with standard procedure: updates dcu, amu
      call dtimer(dtime,itime,-1)
      isign = -1
      call CWFFT3R3(dcu,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call CWFFT3RN(amu,ss,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze, &
     &mdim,nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of time derivative of current with standard
! procedure: updates dcu
      call dtimer(dtime,itime,-1)
      call CADCUPERP3(dcu,amu,nx,ny,nz,nxeh,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate transverse electric field with standard procedure:
! updates cus, wf
      call dtimer(dtime,itime,-1)
      isign = -1
      call CEPOIS33(dcu,cus,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,nz, &
     &nxeh,nye,nze,nxh,nyh,nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform transverse electric field to real space with standard
! procedure: updates cus
      call dtimer(dtime,itime,-1)
      isign = 1
      call CWFFT3R3(cus,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with standard procedure: updates bxyze, cus
      call dtimer(dtime,itime,-1)
      call CCGUARD3L(bxyze,nx,ny,nz,nxe,nye,nze)
      call CCGUARD3L(cus,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and transverse electric fields with standard
! procedure: exyze = cus + fxyze, updates exyze
! cus needs to be retained for next time step
      call dtimer(dtime,itime,-1)
      call CADDVRFIELD3(exyze,cus,fxyze,ndim,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
      enddo
!
! push particles with standard procedure: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
      call CGBPUSH3L(part,exyze,bxyze,qbme,dt,dt,wke,idimp,np,nx,ny,nz, &
     &nxe,nye,nze,ipbc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
!
! sort particles by cell for standard procedure
      if (sortime > 0) then
         if (mod(ntime,sortime)==0) then
            call dtimer(dtime,itime,-1)
            call CDSORTP3YZL(part,part2,npic,idimp,np,ny1,nyz1)
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
