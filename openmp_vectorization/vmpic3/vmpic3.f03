!-----------------------------------------------------------------------
! Skeleton 3D Electrostatic OpenMP/Vector PIC code
! written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE
      program vmpic3
      use iso_c_binding
      use avx512lib3_c
      use kncmpush3_c
      use vmpush3_h
      use omplib_h
      implicit none
! indx/indy/indz = exponent which determines grid points in x/y/z
! direction: nx = 2**indx, ny = 2**indy, nz = 2**indz.
      integer, parameter :: indx =   7, indy =   7, indz =   7
! npx/npy/npz = number of electrons distributed in x/y/z direction.
      integer, parameter :: npx =  384, npy =   384, npz =   384
! ndim = number of velocity coordinates = 3
      integer, parameter :: ndim = 4
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real, parameter :: vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real, parameter :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
! ax/ay/az = smoothed particle size in x/y/z direction
      real :: ax = .912871, ay = .912871, az = .912871
! idimp = number of particle coordinates = 6
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 6, ipbc = 1
! wke/we/wt = particle kinetic/electric field/total energy
      real, target :: wke = 0.0, we = 0.0
      real :: wt = 0.0
! mx/my/mz = number of grids in x/y/z in sorting tiles
      integer :: mx = 8, my = 8, mz = 8
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! kvec = (1,2) = run (autovector,KNC) version
      integer :: kvec = 1
!
! declare scalars for standard code
      integer :: np, nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh
      integer :: nxyzh, nxhyz, mx1, my1, mz1, mxyz1
      integer :: ntime, nloop, isign, lvect
      integer :: irc = 0
      real :: qbme, affp
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, ntmax, npbmx
      integer :: nvp
!
! declare arrays for standard code:
! part = particle arrays
      real, dimension(:,:), pointer :: part => null()
! qe = electron charge density with guard cells
      real, dimension(:,:,:), pointer :: qe => null()
! fxyze = smoothed electric field with guard cells
      real, dimension(:,:,:,:), pointer :: fxyze => null()
! ffc = form factor array for poisson solver
      complex, dimension(:,:,:), pointer :: ffc => null()
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup => null()
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct => null()
!
! declare arrays for OpenMP (tiled) code:
! ppartt = tiled particle array
      real, dimension(:,:,:), pointer :: ppartt => null()
! ppbuff = buffer array for reordering tiled particle array
      real, dimension(:,:,:), pointer :: ppbuff => null()
! kpic = number of particles in each tile
      integer, dimension(:), pointer :: kpic => null()
! ncl = number of particles departing tile in each direction
      integer, dimension(:,:), pointer :: ncl => null()
! ihole = location/destination of each particle departing tile
      integer, dimension(:,:,:), pointer :: ihole => null()
! kp = original location of reordered particle
      integer, dimension(:,:), pointer :: kp => null()
!
! declare and initialize timing data
      real :: time
      type, bind(C) :: timeval
         integer(c_long) :: tv_sec
         integer(c_long) :: tv_usec
      end type
      type (timeval) :: itime
      double precision :: dtime
      real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
      real :: tpush = 0.0, tsort = 0.0
!
      interface
         subroutine dtimer(time,itime,icntrl) bind(C,name='dtimer')
         use iso_c_binding
         import :: timeval
         implicit none
         real(c_double) :: time
         type (timeval) :: itime
         integer(c_int), value :: icntrl
         end subroutine
      end interface
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
      qbme = qme
      affp = real(nx)*real(ny)*real(nz)/real(np)
!
! allocate data for standard code
      allocate(part(idimp,np))
      allocate(mixup(nxhyz),sct(nxyzh))
      allocate(kpic(mxyz1))
!
      lvect = 16
! allocate vector field data
      nxe = lvect*((nxe - 1)/lvect + 1)
      nxeh = nxe/2
      call avx512_f3allocate(qe,nxe,nye,nze,irc)
      call avx512_f4allocate(fxyze,ndim,nxe,nye,nze,irc)
      call avx512_c3allocate(ffc,nxh,nyh,nzh,irc)
      if (irc /= 0) then
         write (*,*) 'aligned field allocation error: irc = ', irc
      endif
!
! prepare fft tables
      call WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
! calculate form factors
      isign = 0
      call VMPOIS33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,  &
     &nye,nze,nxh,nyh,nzh)
! initialize electrons
      call DISTR3(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,idimp,np,nx, &
     &ny,nz,ipbc)
!
! find number of particles in each of mx, my mz, tiles:
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
! align data for Vector Processor
      nppmx0 = lvect*((nppmx0 - 1)/lvect + 1)
      ntmax = lvect*(ntmax/lvect + 1)
      npbmx = lvect*((npbmx - 1)/lvect + 1)
      call avx512_f3allocate(ppartt,nppmx0,idimp,mxyz1,irc)
      call avx512_f3allocate(ppbuff,npbmx,idimp,mxyz1,irc)
      allocate(ncl(26,mxyz1))
      allocate(ihole(2,ntmax+1,mxyz1))
      allocate(kp(nppmx0,mxyz1))
!
! copy ordered particle data for OpenMP: updates ppartt, kpic, and kp
      call PPMOVIN3LTP(part,ppartt,kpic,kp,nppmx0,idimp,np,mx,my,mz,mx1,&
     &my1,mxyz1,irc)
      if (irc /= 0) then
         write (*,*) 'PPMOVIN3LTP overflow error, irc=', irc
         stop
      endif
! sanity check
      call PPCHECK3LT(ppartt,kpic,idimp,nppmx0,nx,ny,nz,mx,my,mz,mx1,my1&
     &,mz1,irc)
      if (irc /= 0) then
         write (*,*) 'PPCHECK3LT error: irc=', irc
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
      call SET_SZERO3(qe,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1)
      if (kvec==1) then
!        call GPPOST3LT(ppartt,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,   &
!    &nye,nze,mx1,my1,mxyz1)
         call VGPPOST3LT(ppartt,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,  &
     &nye,nze,mx1,my1,mxyz1)
! KNC function
      else if (kvec==2) then
         call cknc2gppost3lt(c_loc(ppartt(1,1,1)),c_loc(qe(1,1,1)),     &
     &c_loc(kpic(1)),qme,nppmx0,idimp,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1&
     &)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      if (kvec==1) then
         call AGUARD3L(qe,nx,ny,nz,nxe,nye,nze)
! KNC function
      else if (kvec==2) then
         call ckncaguard3l(c_loc(qe(1,1,1)),nx,ny,nz,nxe,nye,nze)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      isign = -1
      if (kvec==1) then
         call WFFT3RVMX(qe,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze, &
     &nxhyz,nxyzh)
! KNC function
      else if (kvec==2) then
         call ckncwfft3rmx(c_loc(qe(1,1,1)),isign,c_loc(mixup(1)),      &
     &c_loc(sct(1)),indx,indy,indz,nxeh,nye,nze,nxhyz,nxyzh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! calculate force/charge in fourier space with OpenMP: updates fxyze, we
      call dtimer(dtime,itime,-1)
      isign = -1
      if (kvec==1) then
         call VMPOIS33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh&
     &,nye,nze,nxh,nyh,nzh)
! KNC function
      else if (kvec==2) then
         call ckncmpois33(c_loc(qe(1,1,1)),c_loc(fxyze(1,1,1,1)),isign, &
     &c_loc(ffc(1,1,1)),ax,ay,az,affp,we,nx,ny,nz,nxeh,nye,nze,nxh,nyh, &
     &nzh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform force to real space with OpenMP: updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      if (kvec==1) then
         call WFFT3RVM3(fxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,  &
     &nze,nxhyz,nxyzh)
! KNC function
      else if (kvec==2) then
         call ckncwfft3rm3(c_loc(fxyze(1,1,1,1)),isign,c_loc(mixup(1)), &
     &c_loc(sct(1)),indx,indy,indz,nxeh,nye,nze,nxhyz,nxyzh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with OpenMP: updates fxyze
      call dtimer(dtime,itime,-1)
      if (kvec==1) then
         call CGUARD3L(fxyze,nx,ny,nz,nxe,nye,nze)
! KNC function
      else if (kvec==2) then
         call cknccguard3l(c_loc(fxyze(1,1,1,1)),nx,ny,nz,nxe,nye,nze)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! push particles with OpenMP:
      wke = 0.0
      call dtimer(dtime,itime,-1)
      if (kvec==1) then
! updates ppart, wke
!        call GPPUSH3LT(ppartt,fxyze,kpic,qbme,dt,wke,idimp,nppmx0,nx,ny&
!    &,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
         call VGPPUSH3LT(ppartt,fxyze,kpic,qbme,dt,wke,idimp,nppmx0,nx, &
     &ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
!        call V2GPPUSH3LT(ppartt,fxyze,kpic,qbme,dt,wke,idimp,nppmx0,nx,&
!    &ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
! updates ppart, ncl, ihole, wke, irc
!        call GPPUSHF3LT(ppartt,fxyze,kpic,ncl,ihole,qbme,dt,wke,idimp, &
!    &nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,irc)
!        call VGPPUSHF3LT(ppartt,fxyze,kpic,ncl,ihole,qbme,dt,wke,idimp,&
!    &nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,irc)
!        call V2GPPUSHF3LT(ppartt,fxyze,kpic,ncl,ihole,qbme,dt,wke,idimp&
!    &,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,irc)
! KNC function
      else if (kvec==2) then
! updates ppart, wke
         call ckncgppush3lt(c_loc(ppartt(1,1,1)),c_loc(fxyze(1,1,1,1)), &
     &c_loc(kpic(1)),qbme,dt,wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,&
     &nze,mx1,my1,mxyz1,ipbc)
! updates ppart, ncl, ihole, wke, irc
!        call ckncgppushf3lt(c_loc(ppartt(1,1,1)),c_loc(fxyze(1,1,1,1)),&
!    &c_loc(kpic(1)),c_loc(ncl(1,1)),c_loc(ihole(1,1,1)),qbme,dt,wke,   &
!    &idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,irc&
!    &)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
      if (irc /= 0) then
         write (*,*) 'VGPPUSHF3LT error: irc=', irc
         stop
      endif
!
! reorder particles by tile with OpenMP:
      call dtimer(dtime,itime,-1)
     if (kvec==1) then
! updates ppartt, ppbuff, kpic, ncl, ihole, and irc
!        call PPORDER3LT(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny&
!    &,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
         call VPPORDER3LT(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx, &
     &ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
! updates ppartt, ppbuff, kpic, ncl, and irc
!        call PPORDERF3LT(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,&
!    &my1,mz1,npbmx,ntmax,irc)
!        call VPPORDERF3LT(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,   &
!    &mx1,my1,mz1,npbmx,ntmax,irc)
!        call V2PPORDERF3LT(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,  &
!    &mx1,my1,mz1,npbmx,ntmax,irc)
! KNC function
      else if (kvec==2) then
! updates ppartt, ppbuff, kpic, ncl, ihole, and irc
         call ckncpporder3lt(c_loc(ppartt(1,1,1)),c_loc(ppbuff(1,1,1)), &
     &c_loc(kpic(1)),c_loc(ncl(1,1)),c_loc(ihole(1,1,1)),idimp,nppmx0,  &
     &nx,ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
! updates ppartt, ppbuff, kpic, ncl, and irc
!        call ckncpporderf3lt(c_loc(ppartt(1,1,1)),c_loc(ppbuff(1,1,1)),&
!    &c_loc(kpic(1)),c_loc(ncl(1,1)),c_loc(ihole(1,1,1)),idimp,nppmx0,  &
!    &mx1,my1,mz1,npbmx,ntmax,irc)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tsort = tsort + time
      if (irc /= 0) then
         write (*,*) 'VPPORDERF3LT error: ntmax, irc=', ntmax, irc
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
      write (*,*) 'ntime = ', ntime, 'kvec = ', kvec
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
      call avx512_deallocate(c_loc(ppartt(1,1,1)),irc); nullify(ppartt)
      call avx512_deallocate(c_loc(ppbuff(1,1,1)),irc); nullify(ppbuff)
      call avx512_deallocate(c_loc(ffc(1,1,1)),irc); nullify(ffc)
      call avx512_deallocate(c_loc(fxyze(1,1,1,1)),irc); nullify(fxyze)
      call avx512_deallocate(c_loc(qe(1,1,1)),irc); nullify(qe)
!
      stop
      end program
