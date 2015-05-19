!-----------------------------------------------------------------------
! Skeleton 3D Darwin OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mdpic3
      use mdpush3_h
      use omplib_h
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
      integer :: idimp = 6, ipbc = 1
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.4, omy = 0.0, omz = 0.0
! ndc = number of corrections in darwin iteration
      integer :: ndc = 1
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
      real :: zero = 0.0
! mx/my/mz = number of grids in x/y/z in sorting tiles
      integer :: mx = 8, my = 8, mz = 8
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! declare scalars for standard code
      integer :: k
      integer :: np, nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh
      integer :: mdim, nxyzh, nxhyz, mx1, my1, mz1, mxyz1
      integer :: ntime, nloop, isign
      real :: qbme, affp, q2m0, wpm, wpmax, wpmin
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, ntmax, npbmx, irc
      integer :: nvp
!
! declare arrays for standard code:
! part = original particle array
      real, dimension(:,:), pointer :: part
! qe = electron charge density with guard cells
      real, dimension(:,:,:), pointer :: qe
! cue = electron current density with guard cells
! dcu = acceleration density with guard cells
! cus = transverse electric field with guard cells
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
! ss = scratch array for WFFT3RMN
      complex, dimension(:,:), pointer :: ss
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
! nx/ny/nz = number of grid points in x/y/z direction
      np = npx*npy*npz; nx = 2**indx; ny = 2**indy; nz = 2**indz
      nxh = nx/2; nyh = max(1,ny/2); nzh = max(1,nz/2)
      nxe = nx + 2; nye = ny + 1; nze = nz + 1; nxeh = nxe/2
      nxyzh = max(nx,ny,nz)/2; nxhyz = max(nxh,ny,nz)
! mx1/my1/mz1 = number of tiles in x/y/z direction
      mx1 = (nx - 1)/mx + 1; my1 = (ny - 1)/my + 1
      mz1 = (nz - 1)/mz + 1; mxyz1 = mx1*my1*mz1
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
      allocate(qe(nxe,nye,nze),fxyze(ndim,nxe,nye,nze))
      allocate(cue(ndim,nxe,nye,nze),dcu(ndim,nxe,nye,nze))
      allocate(cus(ndim,nxe,nye,nze),amu(mdim,nxe,nye,nze))
      allocate(exyze(ndim,nxe,nye,nze),bxyze(ndim,nxe,nye,nze))
      allocate(ffc(nxh,nyh,nzh),ffe(nxh,nyh,nzh))
      allocate(mixup(nxhyz),sct(nxyzh))
      allocate(ss(mdim*nxeh,nze))
      allocate(kpic(mxyz1))
!
! prepare fft tables
      call WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
! calculate form factor: ffc
      isign = 0
      call MPOIS33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,nye&
     &,nze,nxh,nyh,nzh)
! initialize electrons
      call DISTR3(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,idimp,np,nx, &
     &ny,nz,ipbc)
!
! find number of particles in each of mx, my, mz, tiles:
! updates kpic, nppmx
      call DBLKP3L(part,kpic,nppmx,idimp,np,mx,my,mz,mx1,my1,mxyz1,irc)
      if (irc /= 0) then
         write (*,*) 'DBLKP3L error, irc=', irc
         stop
      endif
! allocate vector particle data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmax = xtras*nppmx
      npbmx = xtras*nppmx
      allocate(ppart(idimp,nppmx0,mxyz1))
      allocate(ppbuff(idimp,npbmx,mxyz1))
      allocate(ncl(26,mxyz1))
      allocate(ihole(2,ntmax+1,mxyz1))
! copy ordered particle data for OpenMP: updates ppart and kpic
      call PPMOVIN3L(part,ppart,kpic,nppmx0,idimp,np,mx,my,mz,mx1,my1,  &
     &mxyz1,irc)
      if (irc /= 0) then
         write (*,*) 'PPMOVIN3L overflow error, irc=', irc
         stop
      endif
! sanity check
      call PPCHECK3L(ppart,kpic,idimp,nppmx0,nx,ny,nz,mx,my,mz,mx1,my1, &
     &mz1,irc)
      if (irc /= 0) then
         write (*,*) 'PPCHECK3L error: irc=', irc
         stop
      endif
!
! find maximum and minimum initial electron density
      qe = 0.0
      call GPPOST3L(ppart,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,nye,nze,&
     &mx1,my1,mxyz1)
      call AGUARD3L(qe,nx,ny,nz,nxe,nye,nze)
      call FWPMINMX3(qe,qbme,wpmax,wpmin,nx,ny,nz,nxe,nye,nze)
      wpm = 0.5*(wpmax + wpmin)*affp
! accelerate convergence: update wpm
      if (wpm <= 10.0) wpm = 0.75*wpm
      write (*,*) 'wpm=',wpm
      q2m0 = wpm/affp
! calculate form factor: ffe
      isign = 0
      call MEPOIS33(dcu,cus,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,nz, &
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
! deposit current with OpenMP: updates cue
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call GJPPOST3L(ppart,cue,kpic,qme,zero,nppmx0,idimp,nx,ny,nz,mx,my&
     &,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdjpost = tdjpost + time
!
! deposit charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call GPPOST3L(ppart,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,nye,nze,&
     &mx1,my1,mxyz1)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with OpenMP: updates qe, cue
      call dtimer(dtime,itime,-1)
      call AGUARD3L(qe,nx,ny,nz,nxe,nye,nze)
      call ACGUARD3L(cue,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      isign = -1
      call WFFT3RMX(qe,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz&
     &,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! calculate longitudinal force/charge in fourier space with OpenMP:
! updates fxyze, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call MPOIS33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,nye&
     &,nze,nxh,nyh,nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform longitudinal electric force to real space with OpenMP:
! updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT3RM3(fxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,  &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! transform current to fourier space with OpenMP: update cue
      call dtimer(dtime,itime,-1)
      isign = -1
      call WFFT3RM3(cue,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of current with OpenMP: updates cue
      call dtimer(dtime,itime,-1)
      call MCUPERP3(cue,nx,ny,nz,nxeh,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate magnetic field in fourier space with OpenMP:
! updates bxyze, wm
      call dtimer(dtime,itime,-1)
      call MBBPOIS33(cue,bxyze,ffc,ci,wf,nx,ny,nz,nxeh,nye,nze,nxh,nyh, &
     &nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform magnetic force to real space with OpenMP: updates bxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT3RM3(bxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,  &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! add constant to magnetic field with OpenMP: updates bxyze
      call dtimer(dtime,itime,-1)
      call BADDEXT3(bxyze,omx,omy,omz,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! copy guard cells with OpenMP: updates fxyze, bxyze
      call dtimer(dtime,itime,-1)
      call CGUARD3L(fxyze,nx,ny,nz,nxe,nye,nze)
      call CGUARD3L(bxyze,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and old transverse electric fields with OpenMP:
! updates exyze
      call dtimer(dtime,itime,-1)
      call ADDVRFIELD3(exyze,cus,fxyze,ndim,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! deposit electron acceleration density and momentum flux with OpenMP:
! updates dcu, amu
      call dtimer(dtime,itime,-1)
      dcu = 0.0; amu = 0.0
      call GDJPPOST3L(ppart,exyze,bxyze,kpic,dcu,amu,qme,qbme,dt,idimp, &
     &nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1)
! add old scaled electric field with OpenMP: updates dcu
      call ASCFGUARD3L(dcu,cus,q2m0,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdcjpost = tdcjpost + time
!
! add guard cells with OpenMP: updates dcu, amu
      call dtimer(dtime,itime,-1)
      call ACGUARD3L(dcu,nx,ny,nz,nxe,nye,nze)
      call AMCGUARD3L(amu,nx,ny,nz,nxe,nye,nze,mdim)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform acceleration density and momentum flux to fourier space
! with OpenMP: updates dcu, amu
      call dtimer(dtime,itime,-1)
      isign = -1
      call WFFT3RM3(dcu,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call WFFT3RMN(amu,ss,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze, &
     &mdim,nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of time derivative of current with OpenMP:
! updates dcu
      call dtimer(dtime,itime,-1)
      call MADCUPERP3(dcu,amu,nx,ny,nz,nxeh,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate transverse electric field with OpenMP: updates cus, wf
      call dtimer(dtime,itime,-1)
      isign = -1
      call MEPOIS33(dcu,cus,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,nz, &
     &nxeh,nye,nze,nxh,nyh,nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform transverse electric field to real space with OpenMP:
! updates cus
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT3RM3(cus,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with OpenMP: updates cus
      call dtimer(dtime,itime,-1)
      call CGUARD3L(cus,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and transverse electric fields with OpenMP:
! exyze = cus + fxyze, updates exyze
! cus needs to be retained for next time step
      call dtimer(dtime,itime,-1)
      call ADDVRFIELD3(exyze,cus,fxyze,ndim,nxe,nye,nze)
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
      call GDCJPPOST3L(ppart,exyze,bxyze,kpic,cue,dcu,amu,qme,qbme,dt,  &
     &idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1)
! add scaled electric field with OpenMP: updates dcu
      call ASCFGUARD3L(dcu,cus,q2m0,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdcjpost = tdcjpost + time
!
! add guard cells for current, acceleration density, and momentum flux
! with OpenMP: updates cue, dcu, amu
      call dtimer(dtime,itime,-1)
      call ACGUARD3L(cue,nx,ny,nz,nxe,nye,nze)
      call ACGUARD3L(dcu,nx,ny,nz,nxe,nye,nze)
      call AMCGUARD3L(amu,nx,ny,nz,nxe,nye,nze,mdim)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform current to fourier space with OpenMP: update cue
      call dtimer(dtime,itime,-1)
      isign = -1
      call WFFT3RM3(cue,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of current with OpenMP: updates cue
      call dtimer(dtime,itime,-1)
      call MCUPERP3(cue,nx,ny,nz,nxeh,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate magnetic field in fourier space with OpenMP:
! updates bxyze, wm
      call dtimer(dtime,itime,-1)
      call MBBPOIS33(cue,bxyze,ffc,ci,wm,nx,ny,nz,nxeh,nye,nze,nxh,nyh, &
     &nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform magnetic force to real space with OpenMP: updates bxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT3RM3(bxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,  &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! add constant to magnetic field with OpenMP: updates bxyze
      call dtimer(dtime,itime,-1)
      call BADDEXT3(bxyze,omx,omy,omz,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform acceleration density and momentum flux to fourier space
! with OpenMP: updates dcu, amu
      call dtimer(dtime,itime,-1)
      isign = -1
      call WFFT3RM3(dcu,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call WFFT3RMN(amu,ss,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze, &
     &mdim,nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of time derivative of current with OpenMP:
! updates dcu
      call dtimer(dtime,itime,-1)
      call MADCUPERP3(dcu,amu,nx,ny,nz,nxeh,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate transverse electric field with OpenMP: updates cus, wf
      call dtimer(dtime,itime,-1)
      isign = -1
      call MEPOIS33(dcu,cus,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,nz, &
     &nxeh,nye,nze,nxh,nyh,nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform transverse electric field to real space with OpenMP:
! updates cus
      call dtimer(dtime,itime,-1)
      isign = 1
      call WFFT3RM3(cus,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,    &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with OpenMP: updates bxyze, cus
      call dtimer(dtime,itime,-1)
      call CGUARD3L(bxyze,nx,ny,nz,nxe,nye,nze)
      call CGUARD3L(cus,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and transverse electric fields with OpenMP:
! exyze = cus + fxyze, updates exyze
! cus needs to be retained for next time step
      call dtimer(dtime,itime,-1)
      call ADDVRFIELD3(exyze,cus,fxyze,ndim,nxe,nye,nze)
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
!     call GBPPUSH3L(ppart,exyze,bxyze,kpic,qbme,dt,dt,wke,idimp,nppmx0,&
!    &nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
! updates ppart, ncl, ihole, wke, irc
      call GBPPUSHF3L(ppart,exyze,bxyze,kpic,ncl,ihole,qbme,dt,dt,wke,  &
     &idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
      if (irc /= 0) then
         write (*,*) 'GBPPUSHF3L error: irc=', irc
         stop
      endif
!
! reorder particles by tile with OpenMP:
      call dtimer(dtime,itime,-1)
! updates ppart, ppbuff, kpic, ncl, ihole, and irc
!     call PPORDER3L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny,nz, &
!    &mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
! updates ppart, ppbuff, kpic, ncl, and irc
      call PPORDERF3L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,my1, &
     &mz1,npbmx,ntmax,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tsort = tsort + time
      if (irc /= 0) then
         write (*,*) 'PPORDERF3L error: ntmax, irc=', ntmax, irc
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
      write (*,*)
!
      stop
      end program
