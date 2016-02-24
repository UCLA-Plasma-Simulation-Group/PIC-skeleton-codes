!-----------------------------------------------------------------------
! Skeleton 3D Electrostatic MPI/OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mppic3
      use mppush3_h
      use mpplib3       ! use with mpplib3.f90
!     use mpplib3_h     ! use with mpplib3.f
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
      real :: ax = .912871, ay = .912871, az = .912871
! idimp = number of particle coordinates = 6
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 6, ipbc = 1
! idps = number of partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
      integer :: idps = 4, idds =    2
! wke/we/wt = particle kinetic/electric field/total energy
      real :: wke = 0.0, we = 0.0, wt = 0.0
! mx/my/mz = number of grids in x/y/z in sorting tiles
! sorting tiles, should be less than or equal to 16
      integer :: mx = 8, my = 8, mz = 8
! fraction of extra particles needed for particle management
      real :: xtras = 0.2
! declare scalars for standard code
      integer :: nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh, nnxe
      integer :: nxyzh, nxhyz, mx1, ntime, nloop, isign, ierr
      real :: qbme, affp
      double precision :: np
!
! declare scalars for MPI code
      integer :: ntpose = 1
      integer :: nvpy, nvpz, nvp, idproc, kstrt, npmax, kyp, kzp
      integer :: kxyp, kyzp, kzyp, nypmx, nzpmx, nypmn, nzpmn
      integer :: npp, nps, myp1, mzp1, mxyzp1, mxzyp1
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, nbmaxp, ntmaxp, npbmx, irc
      integer :: nvpp
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), pointer :: part
! qe = electron charge density with guard cells
      real, dimension(:,:,:), pointer :: qe
! fxyze = smoothed electric field with guard cells
      real, dimension(:,:,:,:), pointer :: fxyze
! qt, qs = scalar charge density field arrays in fourier space
      complex, dimension(:,:,:), pointer :: qt, qs
! fxyzt, fxyzs = vector electric field arrays in fourier space
      complex, dimension(:,:,:,:), pointer :: fxyzt, fxyzs
! ffc = form factor array for poisson solver
      complex, dimension(:,:,:), pointer :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct
      double precision, dimension(4) :: wtot, work
!
! declare arrays for MPI code:
! bs/br = complex send/receive buffers for data transpose
      complex, dimension(:,:,:), pointer :: bs, br
! sbufl/sbufr = particle buffers sent to nearby processors
! rbufl/rbufr = particle buffers received from nearby processors
      real, dimension(:,:,:), pointer :: sbufl, sbufr, rbufl, rbufr
! edges(1:2) = lower:upper y boundaries of particle partition
! edges(3:4) = back:front z boundaries of particle partition
      real, dimension(:), pointer  :: edges
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! noff(1:2) = lowermost global gridpoint in y/z
      integer, dimension(:), pointer :: nyzp, noff
! scr/scs = guard cell buffers received/sent from nearby processors
      real, dimension(:,:), pointer  :: scr, scs
!
! declare arrays for OpenMP code
! ppart = tiled particle array
! ppbuff = buffer array for reordering tiled particle array
      real, dimension(:,:,:), pointer :: ppart, ppbuff
! kpic = number of particles in each tile
      integer, dimension(:), pointer :: kpic
! ncl = number of particles departing tile in each direction
      integer, dimension(:,:), pointer :: ncl
! iholep = location/destination of each particle departing tile
      integer, dimension(:,:,:), pointer :: iholep
! ncll/nclr/mcll/mclr = number offsets send/received from processors
      integer, dimension(:,:,:,:), pointer :: ncll, nclr, mcll, mclr
! mcls = number offsets received from corner processors
      integer, dimension(:,:,:), pointer :: mcls
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0
      real :: tpush = 0.0, tsort = 0.0, tmov = 0.0
      real, dimension(2) :: tfft = 0.0
      double precision :: dtime
!
      irc = 0
! nvpp = number of shared memory nodes (0=default)
      nvpp = 0
!     write (*,*) 'enter number of nodes:'
!     read (5,*) nvpp
! initialize for shared memory parallel processing
      call INIT_OMP(nvpp)
!
! initialize scalars for standard code
! np = total number of particles in simulation
      np =  dble(npx)*dble(npy)*dble(npz)
! nx/ny/nz = number of grid points in x/y/z direction
      nx = 2**indx; ny = 2**indy; nz = 2**indz
      nxh = nx/2; nyh = max(1,ny/2); nzh = max(1,nz/2)
      nxe = nx + 2; nye = ny + 2; nze = nz + 2
      nxeh = nxe/2; nnxe = ndim*nxe
      nxyzh = max(nx,ny,nz)/2; nxhyz = max(nxh,ny,nz)
! mx1 = number of tiles in x direction
      mx1 = (nx - 1)/mx + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = dble(nx)*dble(ny)*dble(nz)/np
!      
! nvp = number of MPI ranks
! initialize for distributed memory parallel processing
      call PPINIT2(idproc,nvp)
      kstrt = idproc + 1
! obtain 2D partition (nvpy,nvpz) from nvp:
! nvpy/nvpz = number of processors in y/z
      call FCOMP32(nvp,nx,ny,nz,nvpy,nvpz,ierr)
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'FCOMP32 error: nvp,nvpy,nvpz=', nvp, nvpy, nvpz
         endif
         go to 3000
      endif
!
! initialize data for MPI code
      allocate(edges(idps),nyzp(idds),noff(idds))
! calculate partition variables:
! edges, nyzp, noff, nypmx, nzpmx, nypmn, nzpmn
! edges(1:2) = lower:upper boundary of particle partition in y
! edges(3:4) = back:front boundary of particle partition in z
! nyzp(1:2) = number of primary (complete) gridpoints in y/z
! noff(1:2) = lowermost global gridpoint in y/z in particle partition
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
! nypmn = minimum value of nyzp(1)
! nzpmn = minimum value of nyzp(2)
      call PDICOMP32L(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,ny,nz,    &
     &kstrt,nvpy,nvpz,idps,idds)
      if (kstrt==1) then
         if (nypmn < 1) then
            write (*,*) 'combination not supported nvpy, ny =',nvpy,ny
         endif
         if (nzpmn < 1) then
            write (*,*) 'combination not supported nvpz, nz =',nvpz,nz
         endif
      endif
      if ((nypmn < 1).or.(nzpmn < 1)) go to 3000
!
! initialize additional scalars for MPI code
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvpy + 1
! kzp = number of complex grids in each field partition in z direction
      kzp = (nz - 1)/nvpz + 1
! kxyp = number of complex grids in each field partition in x direction
! in transposed data
      kxyp = (nxh - 1)/nvpy + 1
! kyzp = number of complex grids in each field partition in y direction,
! in transposed data
      kyzp = (ny - 1)/nvpz + 1; kzyp = max(kyzp,kyp)
! npmax = maximum number of electrons in each partition
      npmax = (np/nvp)*1.25
! myp1/mzp1 = number of tiles in y/z direction
      myp1 = (nyzp(1) - 1)/my + 1; mzp1 = (nyzp(2) - 1)/mz + 1
! mxzyp1 = mx1*max(max(mzp1),max(myp1))
      mxzyp1 = mx1*max((nzpmx-2)/mz+1,(nypmx-2)/my+1)
      mxyzp1 = mx1*myp1*mzp1
!
! allocate data for standard code
      allocate(part(idimp,npmax))
      allocate(qe(nxe,nypmx,nzpmx),fxyze(ndim,nxe,nypmx,nzpmx))
      allocate(qt(nze,kxyp,kyzp),qs(nye,kxyp,nzpmx))
      allocate(fxyzt(ndim,nze,kxyp,kyzp),fxyzs(ndim,nye,kxyp,nzpmx))
      allocate(ffc(nzh,kxyp,kyzp),mixup(nxhyz),sct(nxyzh))
      allocate(kpic(mxyzp1))
!
! allocate data for MPI code
      allocate(bs(ndim,kxyp*kzyp,kzp),br(ndim,kxyp*kzyp,kzp))
      allocate(scr(nxe,nypmx),scs(nnxe,2*nzpmx))
!
! prepare fft tables
      call WPFFT32RINIT(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
! calculate form factors
      isign = 0
      call MPPOIS332(qt,fxyzt,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,kstrt,&
     &nvpy,nvpz,nze,kxyp,kyzp,nzh)
!
! initialize electrons
      nps = 1
      npp = 0
      call PDISTR32(part,edges,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy, &
     &npz,nx,ny,nz,idimp,npmax,idps,ipbc,ierr)
! check for particle initialization error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'particle initialization error: ierr=', ierr
         endif
         go to 3000
      endif
!
! find number of particles in each of mx, my, mz tiles:
! updates kpic, nppmx
      call PPDBLKP3L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,mz,mx1, &
     &myp1,mxyzp1,idds,irc)
      if (irc /= 0) then
         write (*,*) 'PPDBLKP3L error, irc=', irc
         call PPABORT()
         stop
      endif
! allocate vector particle data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmaxp = xtras*nppmx
      npbmx = xtras*nppmx
      nbmaxp = 0.125*mxzyp1*npbmx
      allocate(sbufl(idimp,nbmaxp,2),sbufr(idimp,nbmaxp,2))
      allocate(rbufl(idimp,nbmaxp,2),rbufr(idimp,nbmaxp,2))
      allocate(ppart(idimp,nppmx0,mxyzp1))
      allocate(ppbuff(idimp,npbmx,mxyzp1))
      allocate(ncl(26,mxyzp1))
      allocate(iholep(2,ntmaxp+1,mxyzp1))
      allocate(ncll(3,mxzyp1,3,2),nclr(3,mxzyp1,3,2))
      allocate(mcll(3,mxzyp1,3,2),mclr(3,mxzyp1,3,2))
      allocate(mcls(3,mx1+1,4))
!
! copy ordered particle data for OpenMP: updates ppart and kpic
      call PPPMOVIN3L(part,ppart,kpic,npp,noff,nppmx0,idimp,npmax,mx,my,&
     &mz,mx1,myp1,mxyzp1,idds,irc)
      if (irc /= 0) then
         write (*,*) kstrt, 'PPPMOVIN3L overflow error, irc=', irc
         call PPABORT()
         stop
      endif
! sanity check
      call PPPCHECK3L(ppart,kpic,noff,nyzp,idimp,nppmx0,nx,mx,my,mz,mx1,&
     &myp1,mzp1,idds,irc)
      if (irc /= 0) then
         write (*,*) kstrt, 'PPPCHECK3L error: irc=', irc
         call PPABORT()
         stop
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     if (kstrt==1) write (*,*) 'ntime = ', ntime
!
! deposit charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call PPGPPOST32L(ppart,qe,kpic,noff,qme,nppmx0,idimp,mx,my,mz,nxe,&
     &nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      call PPAGUARD32XL(qe,nyzp,nx,nxe,nypmx,nzpmx,idds)
      call PPNAGUARD32L(qe,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxe,nypmx,   &
     &nzpmx,idds)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with OpenMP: updates qt, modifies qe
      call dtimer(dtime,itime,-1)
      isign = -1
      call WPPFFT32RM(qe,qs,qt,bs,br,isign,ntpose,mixup,sct,ttp,indx,   &
     &indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,kxyp,    &
     &nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! calculate force/charge in fourier space with OpenMP: updates fxyzt, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call MPPOIS332(qt,fxyzt,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,kstrt,&
     &nvpy,nvpz,nze,kxyp,kyzp,nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform force to real space with OpenMP: updates fxyze
! modifies fxyzt
      call dtimer(dtime,itime,-1)
      isign = 1
      call WPPFFT32RM3(fxyze,fxyzs,fxyzt,bs,br,isign,ntpose,mixup,sct,  &
     &ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,&
     &kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! copy guard cells with OpenMP: updates fxyze
      call dtimer(dtime,itime,-1)
      call PPNCGUARD32L(fxyze,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,&
     &idds)
      call PPCGUARD32XL(fxyze,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! push particles with OpenMP:
      call dtimer(dtime,itime,-1)
      wke = 0.0
! updates ppart and wke
!     call PPGPPUSH32L(ppart,fxyze,kpic,noff,nyzp,qbme,dt,wke,idimp,    &
!    &nppmx0,nx,ny,nz,mx,my,mz,nxe,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc&
!    &)
! updates ppart, wke, ncl, iholep, irc
      call PPGPPUSHF32L(ppart,fxyze,kpic,ncl,iholep,noff,nyzp,qbme,dt,  &
     &wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nypmx,nzpmx,mx1,myp1,mxyzp1&
     &,ntmaxp,idds,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
      if (irc /= 0) then
         write (*,*) kstrt, 'PPGPPUSHF32L error: irc=', irc
         call PPABORT()
         stop
      endif
!
! reorder particles by tile with OpenMP
! first part of particle reorder on x, y and z cell
! with mx, my, mz tiles:
      call dtimer(dtime,itime,-1)
! updates ppart, ppbuff, sbufl, sbufr, ncl, iholep, ncll, nclr, irc
!     call PPPORDER32LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,iholep,ncll,  &
!    &nclr,noff,nyzp,idimp,nppmx0,nx,ny,nz,mx,my,mz,mx1,myp1,mzp1,mxzyp1&
!    &,npbmx,ntmaxp,nbmaxp,idds,irc)
! updates: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
      call PPPORDERF32LA(ppart,ppbuff,sbufl,sbufr,ncl,iholep,ncll,nclr, &
     &idimp,nppmx0,mx1,myp1,mzp1,mxzyp1,npbmx,ntmaxp,nbmaxp,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tsort = tsort + time
      if (irc /= 0) then
         write (*,*) kstrt,'PPPORDERF32LA error: ntmaxp,irc=',ntmaxp,irc
         call PPABORT()
         stop
      endif
! move particles into appropriate spatial regions:
! updates rbufr, rbufl, mcll, mclr, mcls
      call dtimer(dtime,itime,-1)
      call PPPMOVE32(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,mcls,  &
     &kstrt,nvpy,nvpz,idimp,nbmaxp,mx1,myp1,mzp1,mxzyp1,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tmov = tmov + time
      if (irc /= 0) then
         if (kstrt==1) then
            write (*,*) kstrt,'PPPMOVE32 error: nbmaxp, irc=',nbmaxp,irc
         endif
         go to 3000
      endif
! second part of particle reorder on x and y cell with mx, my, mz tiles:
! updates ppart, kpic
      call dtimer(dtime,itime,-1)
      call PPPORDER32LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,iholep,mcll,  &
     &mclr,mcls,idimp,nppmx0,mx1,myp1,mzp1,mxzyp1,npbmx,ntmaxp,nbmaxp,  &
     &irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tsort = tsort + time
      if (irc /= 0) then
         write (*,*) kstrt,'PPPORDER32LB error: nppmx0, irc=',nppmx0,irc
         call PPABORT()
         stop
      endif
!
! energy diagnostic
      wtot(1) = we
      wtot(2) = wke
      wtot(3) = 0.0
      wtot(4) = we + wke
      call PPDSUM(wtot,work,4)
      we = wtot(1)
      wke = wtot(2)
      if (ntime==0) then
         if (kstrt==1) then
            write (*,*) 'Initial Field, Kinetic and Total Energies:'
            write (*,'(3e14.7)') we, wke, wke + we
         endif
      endif
      ntime = ntime + 1
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
      if (kstrt==1) then
         write (*,*) 'ntime = ', ntime
         write (*,*) 'MPI nodes nvpy, nvpz = ', nvpy, nvpz
         write (*,*) 'Final Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') we, wke, wke + we
!
         write (*,*)
         write (*,*) 'deposit time = ', tdpost
         write (*,*) 'guard time = ', tguard
         write (*,*) 'solver time = ', tfield
         write (*,*) 'fft and transpose time = ', tfft(1), tfft(2)
         write (*,*) 'push time = ', tpush
         write (*,*) 'particle move time = ', tmov
         write (*,*) 'sort time = ', tsort
         tfield = tfield + tguard + tfft(1)
         write (*,*) 'total solver time = ', tfield
         tsort = tsort + tmov
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
      endif
!
 3000 continue
      call PPEXIT()
      end program
