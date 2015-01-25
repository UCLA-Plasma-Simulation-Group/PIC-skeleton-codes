!-----------------------------------------------------------------------
! Skeleton 2D Electrostatic GPU PIC code
! written by Viktor K. Decyk, UCLA
      program gpupic2
      use gpulib2_h
      use gpupush2_h
      use gpufft2_h
      use push2_h
      implicit none
! indx/indy = exponent which determines grid points in x/y direction:
! nx = 2**indx, ny = 2**indy.
      integer, parameter :: indx =   9, indy =   9
! npx/npy = number of electrons distributed in x/y direction.
      integer, parameter :: npx =  3072, npy =   3072
! ndim = number of velocity coordinates = 2
      integer, parameter :: ndim = 2
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
! vtx/vty = thermal velocity of electrons in x/y direction
! vx0/vy0 = drift velocity of electrons in x/y direction.
      real, parameter :: vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0
! ax/ay = smoothed particle size in x/y direction
      real :: ax = .912871, ay = .912871
! idimp = number of particle coordinates = 4
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 4, ipbc = 1
! wke/we/wt = particle kinetic/electric field/total energy
      real :: wke = 0.0, we = 0.0, wt = 0.0
! mx/my = number of grids in x/y in sorting tiles
      integer :: mx = 16, my = 16
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! declare scalars for standard code
      integer :: np, nx, ny, nxh, nyh, nxh1, nxe, nye, nxeh, nxyh, nxhy
      integer :: mx1, my1, mxy1, ntime, nloop, isign
      real :: qbme, affp
      real, dimension(1) :: sum
!
! declare scalars for GPU code
      integer :: nblock = 128
! nscache = (0,1,2) = (no,small,big) cache size
      integer :: nscache = 1
      integer :: mmcc, nppmx, nppmx0, ntmax, npbmx
      integer :: nxhd
      integer, dimension(1) :: irc
!
! declare arrays for standard code:
! part = original particle array
      real, dimension(:,:), pointer :: part
! ffct = form factor array for poisson solver
      complex, dimension(:,:), pointer :: ffct
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct
!
! declare arrays for GPU code:
! g_qe = electron charge density with guard cells
! g_fxye = smoothed electric field with guard cells
      integer, dimension(2) :: g_qe = 0.0, g_fxye = 0.0
! g_ffct = form factor array for poisson solver
      integer, dimension(2) :: g_ffct = 0.0
! g_mixup = bit reverse table for FFT
! g_sct = sine/cosine table for FFT
      integer, dimension(2) :: g_mixup = 0.0, g_sct = 0.0
! g_q = scalar charge density field array in real space
! g_qt = scalar charge density field array in fourier space
      integer, dimension(2) :: g_q = 0.0, g_qt = 0.0
! g_fxy = vector electric field array in real space
! g_fxyt = vector electric field array in fourier space
      integer, dimension(2) :: g_fxy = 0.0, g_fxyt = 0.0
! g_wke = particle kinetic energy
! g_we = electric field energy
      integer, dimension(2) :: g_wke = 0.0, g_we = 0.0
! g_ppart = tiled particle array
! g_ppbuff = buffer array for reordering tiled particle array
      integer, dimension(2) :: g_ppart = 0.0, g_ppbuff = 0.0
! g_kpic = number of particles in each tile
      integer, dimension(2) :: g_kpic = 0.0
! g_ncl = number of particles departing tile in each direction
! g_ihole = location/destination of each particle departing tile
      integer, dimension(2) :: g_ncl = 0.0, g_ihole = 0.0
! g_sum = scratch array for energy sum reductions
      integer, dimension(2) :: g_sum = 0.0
! g_irc = error code (returned only if error occurs)
      integer, dimension(2) :: g_irc = 0.0
! qt = charge density array in fourier space on host
      complex, dimension(:,:), pointer :: qt
! fxyt = electric field array in fourier space on host
      complex, dimension(:,:,:), pointer :: fxyt
! ppart = tiled particle array on host
      real, dimension(:,:,:), pointer :: ppart
! kpic = number of particles in each tile on host
      integer, dimension(:), pointer :: kpic
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      double precision :: dtime
      real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
      real :: tpush = 0.0, tsort = 0.0
!
! initialize scalars for standard code
! np = total number of particles in simulation
! nx/ny = number of grid points in x/y direction
      np = npx*npy; nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = ny/2
      nxh1 = nxh + 1; nxe = nx + 2; nye = ny + 1; nxeh = nxe/2
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny)
! mx1/my1 = number of tiles in x/y direction
      mx1 = (nx - 1)/mx + 1; my1 = (ny - 1)/my + 1; mxy1 = mx1*my1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx*ny)/real(np)
! set size for FFT arrays
      nxhd = nxh1
!
! allocate data for standard code
      allocate(part(idimp,np))
      allocate(ffct(nyh,nxh))
      allocate(mixup(nxhy),sct(nxyh))
      allocate(kpic(mxy1))
      allocate(qt(ny,nxh1),fxyt(ny,ndim,nxh1))
!
! set up GPU
      irc = 0
      call gpu_setgbsize(nblock)
      call init_cu(0,irc(1))
      if (irc(1) /= 0) then
         write (*,*) 'CUDA initialization error!'
         stop
      endif
! obtain compute capability
      mmcc = getmmcc()
      if (mmcc < 20) then
         write (*,*) 'compute capability 2.x or higher required'
         stop
      endif
! set cache size
      call gpu_set_cache_size(nscache)
!
! allocate data for GPU code
      call gpu_fallocate(g_qe,nxe*nye,irc(1))
      call gpu_fallocate(g_fxye,ndim*nxe*nye,irc(1))
      call gpu_callocate(g_ffct,nyh*nxh,irc(1))
      call gpu_iallocate(g_mixup,nxhy,irc(1))
      call gpu_callocate(g_sct,nxyh,irc(1))
      call gpu_callocate(g_q,nxhd*ny,irc(1))
      call gpu_callocate(g_qt,ny*nxh1,irc(1))
      call gpu_callocate(g_fxy,nxhd*ndim*ny,irc(1))
      call gpu_callocate(g_fxyt,ny*ndim*nxh1,irc(1))
      call gpu_fallocate(g_wke,mxy1,irc(1))
      call gpu_fallocate(g_we,nxh1,irc(1))
      call gpu_fallocate(g_sum,1,irc(1))
      if (irc(1) /= 0) then
         write (*,*) 'GPU allocate error!'
         stop
      endif
!
! prepare fft tables
      call WFFT2RINIT(mixup,sct,indx,indy,nxhy,nxyh)
! prepare NVIDIA ffts
      call gpufft2rrcuinit(nx,ny,ndim)
      call gpufft2cuinit(nx,ny,ndim)
! calculate form factors
      isign = 0
      call POIS22T(qt,fxyt,isign,ffct,ax,ay,affp,we,nx,ny,nxh1,ny,nxh,  &
     &nyh)
! copy in solver arrays to GPU
      call gpu_icopyin(mixup,g_mixup,nxhy)
      call gpu_ccopyin(sct,g_sct,nxyh)
      call gpu_ccopyin(ffct,g_ffct,nyh*nxh)
! initialize electrons
      call DISTR2(part,vtx,vty,vx0,vy0,npx,npy,idimp,np,nx,ny,ipbc)
!
! find number of particles in each of mx, my tiles: updates kpic, nppmx
      call DBLKP2L(part,kpic,nppmx,idimp,np,mx,my,mx1,mxy1,irc)
      if (irc(1) /= 0) then
         write (*,*) 'DBLKP2L error, irc=', irc
         stop
      endif
! allocate vector particle data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmax = xtras*nppmx
      npbmx = xtras*nppmx
! align data to warp size
      nppmx0 = 32*((nppmx0 - 1)/32 + 1)
      ntmax = 32*(ntmax/32 + 1)
      npbmx = 32*((npbmx - 1)/32 + 1)
      call gpu_fallocate(g_ppart,nppmx0*idimp*mxy1,irc(1))
      call gpu_fallocate(g_ppbuff,npbmx*idimp*mxy1,irc(1))
      call gpu_iallocate(g_kpic,mxy1,irc(1))
      call gpu_iallocate(g_ncl,8*mxy1,irc(1))
      call gpu_iallocate(g_ihole,2*(ntmax+1)*mxy1,irc(1))
      call gpu_iallocate(g_irc,1,irc(1))
      if (irc(1) /= 0) then
         write (*,*) 'GPU allocate error!'
         stop
      endif
      allocate(ppart(nppmx0,idimp,mxy1))
!
! copy ordered particle data for GPU code: updates ppart and kpic
      call PPMOVIN2LT(part,ppart,kpic,nppmx0,idimp,np,mx,my,mx1,mxy1,irc&
     &)
      if (irc(1) /= 0) then
         write (*,*) 'PPMOVIN2LT overflow error, irc=', irc
         stop
      endif
! sanity check
      call PPCHECK2LT(ppart,kpic,idimp,nppmx0,nx,ny,mx,my,mx1,my1,irc)
      if (irc(1) /= 0) then
         write (*,*) 'PPCHECK2LT error: irc=', irc
         stop
      endif
! copy to GPU
      call gpu_icopyin(irc,g_irc,1)
      call gpu_fcopyin(ppart,g_ppart,nppmx0*idimp*mxy1)
      call gpu_icopyin(kpic,g_kpic,mxy1)
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     write (*,*) 'ntime = ', ntime
!
! deposit charge with GPU code: updates g_qe
      call dtimer(dtime,itime,-1)
      call gpu_zfmem(g_qe,nxe*nye)
      call cgpu2ppost2l(g_ppart,g_qe,g_kpic,qme,nppmx0,idimp,mx,my,nxe, &
     &nye,mx1,mxy1)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add and copy guard cells with GPU code: updates g_q
      call dtimer(dtime,itime,-1)
      call cgpucaguard2l(g_q,g_qe,nx,ny,nxe,nye,nxhd,ny)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with GPU code: updates g_q, g_qt
      call dtimer(dtime,itime,-1)
      isign = -1
      call cgpuwfft2rcs(g_q,g_qt,isign,g_mixup,g_sct,indx,indy,nxhd,ny, &
     &nxhy,nxyh)
! NVIDIA fft
!     call gpufft2rrcu(g_q,g_qt,isign,indx,indy,nxhd,ny)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! calculate force/charge in fourier space with GPU code:
! updates g_fxyt, g_we
      call dtimer(dtime,itime,-1)
      call cgpupois22t(g_qt,g_fxyt,g_ffct,g_we,nx,ny,nxh1,ny,nxh,nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform force to real space with GPU code: updates g_fxy, g_fxyt
      call dtimer(dtime,itime,-1)
      isign = 1
      call cgpuwfft2rcsn(g_fxy,g_fxyt,isign,g_mixup,g_sct,indx,indy,ndim&
     &,nxhd,ny,nxhy,nxyh)
! NVIDIA fft
!     call gpufft2rrcun(g_fxy,g_fxyt,isign,indx,indy,ndim,nxhd,ny)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with GPU code: updates g_fxye
      call dtimer(dtime,itime,-1)
      call cgpuccguard2l(g_fxy,g_fxye,nx,ny,nxe,nye,nxhd,ny)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! push particles with GPU code:
      call dtimer(dtime,itime,-1)
! updates g_ppart, g_wke
      call cgpuppush2l(g_ppart,g_fxye,g_kpic,qbme,dt,g_wke,idimp,nppmx0,&
     &nx,ny,mx,my,nxe,nye,mx1,mxy1,ipbc)
! updates g_ppart, g_ncl, g_ihole, g_wke, g_irc
!     call cgpuppushf2l(g_ppart,g_fxye,g_kpic,g_ncl,g_ihole,qbme,dt,    &
!    &g_wke,idimp,nppmx0,nx,ny,mx,my,nxe,nye,mx1,mxy1,ntmax,g_irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
!
! reorder particles by tile with GPU code:
      call dtimer(dtime,itime,-1)
! updates g_ppart, g_ppbuff, g_kpic, g_ncl, g_ihole, and g_irc
      call cgpuppord2l(g_ppart,g_ppbuff,g_kpic,g_ncl,g_ihole,idimp,     &
     &nppmx0,nx,ny,mx,my,mx1,my1,npbmx,ntmax,g_irc)
! updates g_ppart, g_ppbuff, g_kpic, g_ncl, and g_irc
!     call cgpuppordf2l(g_ppart,g_ppbuff,g_kpic,g_ncl,g_ihole,idimp,    &
!    &nppmx0,mx1,my1,npbmx,ntmax,g_irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tsort = tsort + time
!
! sanity check
      call gpu_icopyout(irc,g_irc,1)
      if (irc(1) /= 0) then
         write (*,*) 'push or reorder error: ntmax, irc=', ntmax, irc
         stop
      endif
!
! energy diagnostic
      if (ntime==0) then
         call gpu_zfmem(g_sum,1)
         call cgpusum2(g_we,g_sum,nxh1)
         call gpu_fcopyout(sum,g_sum,1); we = sum(1)
         call gpu_zfmem(g_sum,1)
         call cgpusum2(g_wke,g_sum,mxy1)
         call gpu_fcopyout(sum,g_sum,1); wke = sum(1)
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
! energy diagnostic
      call gpu_zfmem(g_sum,1)
      call cgpusum2(g_we,g_sum,nxh1)
      call gpu_fcopyout(sum,g_sum,1); we = sum(1)
      call gpu_zfmem(g_sum,1)
      call cgpusum2(g_wke,g_sum,mxy1)
      call gpu_fcopyout(sum,g_sum,1); wke = sum(1)
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
! close down NVIDIA fft
      call gpufft2cudel()
      call gpufft2rrcudel()
! deallocate memory on GPU
      call gpu_deallocate(g_irc,irc(1))
      call gpu_deallocate(g_ihole,irc(1))
      call gpu_deallocate(g_ncl,irc(1))
      call gpu_deallocate(g_kpic,irc(1))
      call gpu_deallocate(g_sum,irc(1))
      call gpu_deallocate(g_we,irc(1))
      call gpu_deallocate(g_wke,irc(1))
      call gpu_deallocate(g_fxyt,irc(1))
      call gpu_deallocate(g_fxy,irc(1))
      call gpu_deallocate(g_qt,irc(1))
      call gpu_deallocate(g_q,irc(1))
      call gpu_deallocate(g_sct,irc(1))
      call gpu_deallocate(g_mixup,irc(1))
      call gpu_deallocate(g_ffct,irc(1))
      call gpu_deallocate(g_ppbuff,irc(1))
      call gpu_deallocate(g_ppart,irc(1))
      call gpu_deallocate(g_fxye,irc(1))
      call gpu_deallocate(g_qe,irc(1))
! close down GPU
      call end_cu()
!
      stop
      end program
