!-----------------------------------------------------------------------
! Skeleton 2-1/2D Electromagnetic GPU-MPI PIC code
! written by Viktor K. Decyk, UCLA
      program gpupbpic2
      use iso_c_binding
      use gpulib2_c
      use gpupbpush2_c
      use gpupfft2_c
      use pbpush2_h
      use pplib2       ! use with pplib2.f90
!     use pplib2_h     ! use with pplib2.f
      use gpplib2_f03
      use dtimer_c
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
! ci = reciprical of velocity of light.
      real :: ax = .912871, ay = .912871, ci = 0.1
! idimp = number of particle coordinates = 5
! ipbc = particle boundary condition: 1 = periodic
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 5, ipbc = 1, relativity = 1
! idps = number of partition boundaries
      integer :: idps = 2
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
! mx/my = number of grids in x/y in sorting tiles
      integer :: mx = 16, my = 16
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! declare scalars for standard code
      integer :: nx, ny, nxh, nyh, nxh1, nxe, nye, nxeh
      integer :: nxyh, nxhy, mx1, ntime, nloop, isign, ierr
      real :: qbme, affp, dth
      real, dimension(1), target :: ssum
      double precision :: np
!
! declare scalars for MPI code
      integer :: nvp, idproc, kstrt, npmax, kxp, kyp, nypmx, nypmn
      integer :: nyp, noff, npp, nps
      integer :: myp1, mxyp1, kxpp, kypp
!
! declare scalars for GPU code
      integer :: nblock = 128
! nscache = (0,1,2) = (no,small,big) cache size
      integer :: nscache = 1
      integer :: mmcc, nppmx, nppmx0, nbmaxp, ntmaxp, npbmx
      integer :: nxhd, kxpd, idev, ndev
      integer, dimension(1), target :: irc
!
! declare arrays for standard code:
! part = original particle array
      real, dimension(:,:), pointer :: part => null()
! ffct = form factor array for poisson solver
      complex, dimension(:,:), pointer :: ffct => null()
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup => null()
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct => null()
      real, dimension(7) :: wtot, work
!
! declare arrays for MPI code:
! sbufl/sbufr = particle buffers sent to nearby processors
      real, dimension(:,:), pointer :: sbufl => null(), sbufr => null()
! rbufl/rbufr = particle buffers received from nearby processors
      real, dimension(:,:), pointer :: rbufl => null(), rbufr => null()
! edges(1:2) = lower:upper y boundaries of particle partition
      real, dimension(:), pointer  :: edges => null()
! scs/scr = guard cell buffers sent to/received from nearby processors
      real, dimension(:), pointer  :: scs => null(), scr => null()
! locl = ordered list of MPI ranks on the same host
      integer, dimension(:), pointer :: locl => null()
!
! declare arrays for GPU code:
! g_qe = electron charge density with guard cells
      type (c_ptr) :: g_qe = c_null_ptr
! g_cue = electron current density with guard cells
      type (c_ptr) :: g_cue = c_null_ptr
! g_fxyze/g_bxyze = smoothed electric/magnetic field with guard cells
      type (c_ptr) :: g_fxyze = c_null_ptr, g_bxyze = c_null_ptr
! g_ffct = form factor array for poisson solver
      type (c_ptr) :: g_ffct = c_null_ptr
! g_mixup = bit reverse table for FFT
! g_sct = sine/cosine table for FFT
      type (c_ptr) :: g_mixup = c_null_ptr, g_sct = c_null_ptr
! g_q = scalar charge density field array in real space
! g_cu = vector current density field array in real space
      type (c_ptr) :: g_q = c_null_ptr, g_cu = c_null_ptr
! g_qt = scalar charge density field array in fourier space
! g_cut = vector current density field array in fourier space
      type (c_ptr) :: g_qt = c_null_ptr, g_cut = c_null_ptr
! g_fxyz/g_hxyz = vector electric/magnetic field in real space
      type (c_ptr) :: g_fxyz = c_null_ptr, g_hxyz = c_null_ptr
! g_fxyzt/g_hxyzt = vector electric/magnetic field in fourier space
      type (c_ptr) :: g_fxyzt = c_null_ptr, g_hxyzt = c_null_ptr
! g_exyzt/g_bxyzt = transverse electric/magnetic field in fourier space
      type (c_ptr) :: g_exyzt = c_null_ptr, g_bxyzt = c_null_ptr
! g_wke/g_we = particle kinetic/electrostatic field energy
      type (c_ptr) :: g_wke = c_null_ptr, g_we = c_null_ptr
! g_wf/g_wm = magnetic field/transverse electric field energy
      type (c_ptr) :: g_wf = c_null_ptr, g_wm = c_null_ptr
! g_ppart = tiled particle array
! g_ppbuff = buffer array for reordering tiled particle array
      type (c_ptr) :: g_ppart = c_null_ptr, g_ppbuff = c_null_ptr
! g_kpic = number of particles in each tile
      type (c_ptr) :: g_kpic = c_null_ptr
! g_sbufl/g_sbufr = buffers for sending/receiving particle data
      type (c_ptr) :: g_sbufl = c_null_ptr, g_sbufr = c_null_ptr
! g_ncl = number of particles departing tile in each direction
! g_ihole = location/destination of each particle departing tile
      type (c_ptr) :: g_ncl = c_null_ptr, g_ihole = c_null_ptr
! g_ncll/g_nclr = buffers for sending/receiving particle tile info
      type (c_ptr) :: g_ncll = c_null_ptr, g_nclr = c_null_ptr
! g_bsm/g_brm = complex send/receive buffers for data transpose
      type (c_ptr) :: g_bsm = c_null_ptr, g_brm = c_null_ptr
! g_scs = buffer for sending/receiving guard cell data
      type (c_ptr) :: g_scs = c_null_ptr
! g_sum = scratch array for energy sum reductions
      type (c_ptr) :: g_sum = c_null_ptr
! g_irc = error code (returned only if error occurs)
      type (c_ptr) :: g_irc = c_null_ptr
! qt = charge density array in fourier space on host
      complex, dimension(:,:), pointer :: qt => null()
! fxyzt = electric field array in fourier space on host
      complex, dimension(:,:,:), pointer :: fxyzt => null()
! ppart = tiled particle array on host
      real, dimension(:,:,:), pointer :: ppart => null()
! kpic = number of particles in each tile on host
      integer, dimension(:), pointer :: kpic => null()
! ncll/nclr = particle tile info being sent to nearby processors
      integer, dimension(:,:), pointer :: ncll => null(), nclr => null()
! mcll/mclr = particle tile info being received from nearby processors
      integer, dimension(:,:), pointer :: mcll => null(), mclr => null()
! bsm/brm = complex send/receive buffers for data transpose
      complex, dimension(:,:), pointer :: bsm => null(), brm => null()
!
! declare and initialize timing data
      real :: time
      type (timeval) :: itime
      double precision :: dtime
      real :: tdpost = 0.0, tguard = 0.0, tsort = 0.0
      real :: tfield = 0.0, tdjpost = 0.0, tpush = 0.0
      real :: tmov = 0.0
      real, dimension(2) :: tmsort = 0.0
      real, dimension(2) :: tfft = 0.0
!
! initialize scalars for standard code
! np = total number of particles in simulation
      np =  dble(npx)*dble(npy)
! nx/ny = number of grid points in x/y direction
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = ny/2
      nxh1 = nxh + 1
      nxe = nx + 2; nye = ny + 2; nxeh = nxe/2
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny)
! mx1 = number of tiles in x direction
      mx1 = (nx - 1)/mx + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = dble(nx)*dble(ny)/np
      dth = 0.0
! set size for FFT arrays
      nxhd = nxh1
!      
! nvp = number of MPI ranks
! initialize for distributed memory parallel processing
      call PPINIT2(idproc,nvp)
      kstrt = idproc + 1
! check if too many processors
      if (nvp > ny) then
         if (kstrt==1) then
         write (*,*) 'Too many processors requested: ny, nvp=', ny, nvp
         endif
         go to 3000
      endif
!
! initialize data for MPI code
      allocate(edges(idps))
! calculate partition variables: edges, nyp, noff, nypmx
! edges(1:2) = lower:upper boundary of particle partition
! nyp = number of primary (complete) gridpoints in particle partition
! noff = lowermost global gridpoint in particle partition
! nypmx = maximum size of particle partition, including guard cells
! nypmn = minimum value of nyp
      call PDICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,idps)
      if (nypmn < 1) then
         if (kstrt==1) then
            write (*,*) 'combination not supported nvp, ny =',nvp,ny
         endif
         go to 3000
      endif
! initialize additional scalars for MPI code
! kxp = number of complex grids in each field partition in x direction
      kxp = (nxh - 1)/nvp + 1
! set size for FFT arrays
      kxpd = (nxh1 - 1)/nvp + 1
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvp + 1
! npmax = maximum number of electrons in each partition
      npmax = (np/nvp)*1.25
! myp1 = number of tiles in y direction
      myp1 = (nyp - 1)/my + 1; mxyp1 = mx1*myp1
! kxpp/kypp = actual size of GPU field partition
      kxpp = min(kxpd,max(0,nxhd-kxpd*idproc))
      kypp = min(kyp,max(0,ny-kyp*idproc))
!
! allocate data for standard code
      allocate(part(idimp,npmax))
      allocate(ffct(nyh,kxpd))
      allocate(mixup(nxhy),sct(nxyh))
      allocate(kpic(mxyp1))
      allocate(qt(ny,kxpd),fxyzt(ny,ndim,kxpd))
!
! allocate data for MPI code
      allocate(locl(nvp))
!
! set up GPU
      irc = 0
! get unique GPU device ids
      call PPFNDGRP(locl,kstrt,nvp,idev,ndev)
      if (idev < 0) then
         write (*,*) kstrt,'GPU device id error!'
         call PPABORT()
         stop
      endif
      call gpu_setgbsize(nblock)
      call init_cu(idev,irc(1))
      if (irc(1) /= 0) then
         write (*,*) kstrt,'CUDA initialization error!'
         call PPABORT()
         stop
      endif
! obtain compute capability
      mmcc = getmmcc()
      if (mmcc < 20) then
         write (*,*) kstrt, 'compute capability 2.x or higher required'
         call PPABORT()
         stop
      endif
! set cache size
      call gpu_set_cache_size(nscache)
! create asynchronous streams
      call gpu_initstream(1)
      call gpu_initstream(2)
      call gpu_initstream(3)
! allocate data for GPU code
      call gpu_fallocate(g_qe,nxe*nypmx,irc(1))
      call gpu_fallocate(g_cue,ndim*nxe*nypmx,irc(1))
      call gpu_fallocate(g_fxyze,ndim*nxe*nypmx,irc(1))
      call gpu_fallocate(g_bxyze,ndim*nxe*nypmx,irc(1))
      call gpu_callocate(g_ffct,nyh*kxpd,irc(1))
      call gpu_iallocate(g_mixup,nxhy,irc(1))
      call gpu_callocate(g_sct,nxyh,irc(1))
      call gpu_callocate(g_q,nxhd*kyp,irc(1))
      call gpu_callocate(g_cu,nxhd*ndim*kyp,irc(1))
      call gpu_callocate(g_qt,ny*kxpd,irc(1))
      call gpu_callocate(g_cut,ny*ndim*kxpd,irc(1))
      call gpu_callocate(g_fxyz,nxhd*ndim*kyp,irc(1))
      call gpu_callocate(g_hxyz,nxhd*ndim*kyp,irc(1))
      call gpu_callocate(g_fxyzt,ny*ndim*kxpd,irc(1))
      call gpu_callocate(g_hxyzt,ny*ndim*kxpd,irc(1))
      call gpu_callocate(g_exyzt,ny*ndim*kxpd,irc(1))
      call gpu_callocate(g_bxyzt,ny*ndim*kxpd,irc(1))
      call gpu_fallocate(g_wke,mxyp1,irc(1))
      call gpu_fallocate(g_we,kxpd,irc(1))
      call gpu_fallocate(g_wf,kxpd,irc(1))
      call gpu_fallocate(g_wm,kxpd,irc(1))
      call gpu_fallocate(g_sum,1,irc(1))
      if (irc(1) /= 0) then
         write (*,*) kstrt, 'GPU allocate error!'
         call PPABORT()
         stop
      endif
!
! prepare fft tables
      call WPFFT2RINIT(mixup,sct,indx,indy,nxhy,nxyh)
! prepare NVIDIA ffts
      call gpupfft2rrcuinit(nx,kypp,ndim)
      call gpupfft2cuinit(kxpp,ny,ndim)
! calculate form factors
      isign = 0
      call PPOIS23T(qt,fxyzt,isign,ffct,ax,ay,affp,we,nx,ny,kstrt,ny,   &
     &kxpd,nyh)
! copy in solver arrays to GPU
      call gpu_icopyin(c_loc(mixup(1)),g_mixup,nxhy)
      call gpu_ccopyin(c_loc(sct(1)),g_sct,nxyh)
      call gpu_ccopyin(c_loc(ffct(1,1)),g_ffct,nyh*kxpd)
! initialize electrons
      nps = 1
      npp = 0
      call PDISTR2H(part,edges,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy, &
     &nx,ny,idimp,npmax,idps,ipbc,ierr)
! check for particle initialization error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'particle initialization error: ierr=', ierr
         endif
         go to 3000
      endif
!
! initialize transverse electromagnetic fields
      call gpu_zcmem(g_exyzt,ny*ndim*kxpd)
      call gpu_zcmem(g_bxyzt,ny*ndim*kxpd)
!
! find number of particles in each of mx, my tiles: updates kpic, nppmx
      call PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,mx1,    &
     &mxyp1,irc(1))
      if (irc(1) /= 0) then
         write (*,*) kstrt, 'PPDBLKP2L error, irc=', irc(1)
         call PPABORT()
         stop
      endif
! allocate vector particle data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmaxp = 0.5*xtras*nppmx
      npbmx = 0.5*xtras*nppmx
      nbmaxp = 0.25*mx1*npbmx
! align data to warp size
      nppmx0 = 32*((nppmx0 - 1)/32 + 1)
      ntmaxp = 32*(ntmaxp/32 + 1)
      npbmx = 32*((npbmx - 1)/32 + 1)
      nbmaxp = 32*((nbmaxp - 1)/32 + 1)
      call gpu_fallocate(g_ppart,nppmx0*idimp*mxyp1,irc(1))
      call gpu_fallocate(g_ppbuff,npbmx*idimp*mxyp1,irc(1))
      call gpu_iallocate(g_kpic,mxyp1+1,irc(1))
      call gpu_fallocate(g_sbufl,nbmaxp*idimp,irc(1))
      call gpu_fallocate(g_sbufr,nbmaxp*idimp,irc(1))
      call gpu_iallocate(g_ncl,8*mxyp1,irc(1))
      call gpu_iallocate(g_ihole,2*(ntmaxp+1)*mxyp1,irc(1))
      call gpu_iallocate(g_ncll,3*mx1,irc(1))
      call gpu_iallocate(g_nclr,3*mx1,irc(1))
      call gpu_callocate(g_bsm,kxpd*ndim*kyp*nvp,irc(1))
      call gpu_callocate(g_brm,kxpd*ndim*kyp*nvp,irc(1))
      call gpu_fallocate(g_scs,nxe*ndim,irc(1))
      call gpu_iallocate(g_irc,1,irc(1))
      if (irc(1) /= 0) then
         write (*,*) kstrt, 'GPU allocate error!'
         call PPABORT()
         stop
      endif
      allocate(ppart(nppmx0,idimp,mxyp1))
      allocate(ncll(3,mx1),nclr(3,mx1),mcll(3,mx1),mclr(3,mx1))
!
! allocate data for GPU-MPI buffers
!     allocate(scs(nxe*ndim),scr(nxe*ndim))
!     allocate(sbufl(idimp,nbmaxp),sbufr(idimp,nbmaxp))
!     allocate(rbufl(idimp,nbmaxp),rbufr(idimp,nbmaxp))
!     allocate(bsm(kxpd*ndim*kyp,nvp),brm(kxpd*ndim*kyp,nvp))
! allocate host page-locked memory for GPU-MPI buffers
      call hpl_f1allocate(scs,nxe*ndim,irc(1))
      call hpl_f1allocate(scr,nxe*ndim,irc(1))
      call hpl_f2allocate(sbufl,idimp,nbmaxp,irc(1))
      call hpl_f2allocate(sbufr,idimp,nbmaxp,irc(1))
      call hpl_f2allocate(rbufl,idimp,nbmaxp,irc(1))
      call hpl_f2allocate(rbufr,idimp,nbmaxp,irc(1))
      call hpl_c2allocate(bsm,kxpd*ndim*kyp,nvp,irc(1))
      call hpl_c2allocate(brm,kxpd*ndim*kyp,nvp,irc(1))
      if (irc(1) /= 0) then
         write (*,*) kstrt, 'hpl allocate error, irc=', irc(1)
         call PPABORT()
         stop
      endif
!
! copy ordered particle data for GPU code: updates ppart, kpic
      call PPPMOVIN2LT(part,ppart,kpic,npp,noff,nppmx0,idimp,npmax,mx,my&
     &,mx1,mxyp1,irc(1))
      if (irc(1) /= 0) then
         write (*,*) kstrt, 'PPPMOVIN2LT overflow error, irc=', irc(1)
         call PPABORT()
         stop
      endif
! sanity check
      call PPPCHECK2LT(ppart,kpic,noff,nyp,idimp,nppmx0,nx,mx,my,mx1,   &
     &myp1,irc(1))
      if (irc(1) /= 0) then
         write (*,*) kstrt, 'PPPCHECK2LT error: irc=', irc(1)
         call PPABORT()
         stop
      endif
! copy to GPU
      call gpu_icopyin(c_loc(irc),g_irc,1)
      call gpu_fcopyin(c_loc(ppart(1,1,1)),g_ppart,nppmx0*idimp*mxyp1)
      call gpu_icopyin(c_loc(kpic(1)),g_kpic,mxyp1)
      call gpu_zfmem(g_we,kxpd)
      call gpu_zfmem(g_wf,kxpd)
      call gpu_zfmem(g_wm,kxpd)
!
      if (dt > 0.45*ci) then
         if (kstrt==1) then
            write (*,*) 'Warning: Courant condition may be exceeded!'
         endif
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     if (kstrt==1) write (*,*) 'ntime = ', ntime
!
! deposit current with GPU code: updates g_ppart, g_cue
      call dtimer(dtime,itime,-1)
      call gpu_zfmem(g_cue,ndim*nxe*nypmx)
      if (relativity==1) then
         call cgpu2pprjppost2l(g_ppart,g_cue,g_kpic,noff,qme,dth,ci,    &
     &nppmx0,idimp,nx,ny,mx,my,nxe,nypmx,mx1,mxyp1,ipbc)
      else
         call cgpu2ppjppost2l(g_ppart,g_cue,g_kpic,noff,qme,dth,nppmx0, &
     &idimp,nx,ny,mx,my,nxe,nypmx,mx1,mxyp1,ipbc)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdjpost = tdjpost + time
!
! reorder particles by tile with GPU code:
! updates g_ppart, g_ppbuff, g_kpic, g_ncl, g_ihole, and g_irc,
! as well as various buffers
      call GPPORDER2L(g_ppart,g_ppbuff,g_sbufl,g_sbufr,g_kpic,g_ncl,    &
     &g_ihole,g_ncll,g_nclr,sbufl,sbufr,rbufl,rbufr,ncll,nclr,mcll,mclr,&
     &tmsort,noff,nyp,kstrt,nvp,idimp,nppmx0,nx,ny,mx,my,mx1,myp1,npbmx,&
     &ntmaxp,nbmaxp,g_irc)
      tsort = tmsort(1)
      tmov = tmsort(2)
!
! deposit charge with GPU code: updates g_qe
      call dtimer(dtime,itime,-1)
      call gpu_zfmem(g_qe,nxe*nypmx)
      call cgpu2ppgppost2l(g_ppart,g_qe,g_kpic,noff,qme,idimp,nppmx0,mx,&
     &my,nxe,nypmx,mx1,mxyp1)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add and copy guard cells with GPU code: updates g_cu, g_q
      call dtimer(dtime,itime,-1)
      call GPPCACGUARD2L(g_cu,g_cue,g_scs,scs,scr,nx,nyp,kstrt,nvp,ndim,&
     &nxe,nypmx,nxhd,kyp)
      call GPPCAGUARD2L(g_q,g_qe,g_scs,scs,scr,nx,nyp,kstrt,nvp,nxe,    &
     &nypmx,nxhd,kyp)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with GPU code: updates g_q, g_qt,
! as well as various buffers
      isign = -1
      call WAPPFFT2RCS(g_q,g_qt,g_bsm,g_brm,bsm,brm,isign,g_mixup,g_sct,&
     &tfft,indx,indy,kstrt,nvp,kxpd,kyp,nxhd,ny,kyp,nxhy,nxyh)
! NVIDIA fft
!     call GPUPPFFT2RRCU(g_q,g_qt,g_bsm,g_brm,bsm,brm,isign,tfft,indx,  &
!    &indy,kstrt,nvp,kxpd,kyp,nxhd,ny,kyp)
!
! transform current to fourier space with GPU code: updates g_cu, g_cut,
! as well as various buffers
      isign = -1
      call WAPPFFT2RCSN(g_cu,g_cut,g_bsm,g_brm,bsm,brm,isign,g_mixup,   &
     &g_sct,tfft,indx,indy,kstrt,nvp,ndim,kxpd,kyp,nxhd,ny,kyp,nxhy,nxyh&
     &)
! NVIDIA fft
!     call GPUPPFFT2RRCUN(g_cu,g_cut,g_bsm,g_brm,bsm,brm,isign,tfft,indx&
!    &,indy,kstrt,nvp,ndim,kxpd,kyp,nxhd,ny,kyp)
!
! take transverse part of current with GPU code: updates g_cut
      call dtimer(dtime,itime,-1)
      call cgpuppcuperp2t(g_cut,nx,ny,kstrt,ny,kxpd)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate electromagnetic fields in fourier space with GPU code:
! updates g_exyzt, g_bxyzt, g_wf, g_wm
      call dtimer(dtime,itime,-1)
      if (ntime==0) then
         call cgpuippbpoisp23t(g_cut,g_bxyzt,g_ffct,ci,g_wm,nx,ny,kstrt,&
     &ny,kxpd,nyh)
         wf = 0.0
         dth = 0.5*dt
      else
         call cgpuppmaxwel2t(g_exyzt,g_bxyzt,g_cut,g_ffct,affp,ci,dt,   &
     &g_wf,g_wm,nx,ny,kstrt,ny,kxpd,nyh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate force/charge in fourier space with GPU code:
! updates g_fxyzt, g_we
      call dtimer(dtime,itime,-1)
      call cgpuppois23t(g_qt,g_fxyzt,g_ffct,g_we,nx,ny,kstrt,ny,kxpd,nyh&
     &)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! add longitudinal and transverse electric fields with with GPU code:
! updates g_fxyzt
      call dtimer(dtime,itime,-1)
      isign = 1
      call cgpuppemfield2t(g_fxyzt,g_exyzt,g_ffct,isign,nx,ny,kstrt,ny, &
     &kxpd,nyh)
! copy magnetic field with GPU code: updates g_hxyzt
      isign = -1
      call cgpuppemfield2t(g_hxyzt,g_bxyzt,g_ffct,isign,nx,ny,kstrt,ny, &
     &kxpd,nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform force to real space with GPU code: updates g_fxyz, g_fxyzt
! as well as various buffers
      isign = 1
      call WAPPFFT2RCSN(g_fxyz,g_fxyzt,g_bsm,g_brm,bsm,brm,isign,g_mixup&
     &,g_sct,tfft,indx,indy,kstrt,nvp,ndim,kxpd,kyp,nxhd,ny,kyp,nxhy,   &
     &nxyh)
! NVIDIA fft
!     call GPUPPFFT2RRCUN(g_fxyz,g_fxyzt,g_bsm,g_brm,bsm,brm,isign,tfft,&
!    &indx,indy,kstrt,nvp,ndim,kxpd,kyp,nxhd,ny,kyp)
!
! transform magnetic field to fourier space with GPU code:
! updates g_hxyz, g_hxyzt, as well as various buffers
      isign = 1
      call WAPPFFT2RCSN(g_hxyz,g_hxyzt,g_bsm,g_brm,bsm,brm,isign,g_mixup&
     &,g_sct,tfft,indx,indy,kstrt,nvp,ndim,kxpd,kyp,nxhd,ny,kyp,nxhy,   &
     &nxyh)
! NVIDIA fft
!     call GPUPPFFT2RRCUN(g_hxyz,g_hxyzt,g_bsm,g_brm,bsm,brm,isign,tfft,&
!    &indx,indy,kstrt,nvp,ndim,kxpd,kyp,nxhd,ny,kyp)
!
! copy guard cells with GPU code: updates g_fxyze, g_bxyze
      call dtimer(dtime,itime,-1)
      call GPPCBGUARD2L(g_fxyz,g_fxyze,g_scs,scs,scr,nx,nyp,kstrt,nvp,  &
     &ndim,nxe,nypmx,nxhd,kyp)
      call GPPCBGUARD2L(g_hxyz,g_bxyze,g_scs,scs,scr,nx,nyp,kstrt,nvp,  &
     &ndim,nxe,nypmx,nxhd,kyp)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! push particles with GPU code: updates g_ppart, g_wke
      call dtimer(dtime,itime,-1)
      if (relativity==1) then
         call cgpuppgrbppush23l(g_ppart,g_fxyze,g_bxyze,g_kpic,noff,nyp,&
     &qbme,dt,dth,ci,g_wke,idimp,nppmx0,nx,ny,mx,my,nxe,nypmx,mx1,mxyp1,&
     &ipbc)
      else
         call cgpuppgbppush23l(g_ppart,g_fxyze,g_bxyze,g_kpic,noff,nyp, &
     &qbme,dt,dth,g_wke,idimp,nppmx0,nx,ny,mx,my,nxe,nypmx,mx1,mxyp1,   &
     &ipbc)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
!
! reorder particles by tile with GPU code:
! updates g_ppart, g_ppbuff, g_kpic, g_ncl, g_ihole, and g_irc,
! as well as various buffers
      call GPPORDER2L(g_ppart,g_ppbuff,g_sbufl,g_sbufr,g_kpic,g_ncl,    &
     &g_ihole,g_ncll,g_nclr,sbufl,sbufr,rbufl,rbufr,ncll,nclr,mcll,mclr,&
     &tmsort,noff,nyp,kstrt,nvp,idimp,nppmx0,nx,ny,mx,my,mx1,myp1,npbmx,&
     &ntmaxp,nbmaxp,g_irc)
      tsort = tmsort(1)
      tmov = tmsort(2)
!
! sanity check
      call gpu_icopyout(c_loc(irc),g_irc,1)
      if (irc(1) /= 0) then
         write (*,*) kstrt, 'GPPORDER2L error: irc=', irc(1)
         call PPABORT()
         stop
      endif
!
! energy diagnostic
      if (ntime==0) then
         call gpu_zfmem(g_sum,1)
         call cgpusum2(g_we,g_sum,kxpd)
         call gpu_fcopyout(c_loc(ssum),g_sum,1); we = ssum(1)
         call gpu_zfmem(g_sum,1)
         call cgpusum2(g_wf,g_sum,kxpd)
         call gpu_fcopyout(c_loc(ssum),g_sum,1); wf = ssum(1)
         call gpu_zfmem(g_sum,1)
         call cgpusum2(g_wm,g_sum,kxpd)
         call gpu_fcopyout(c_loc(ssum),g_sum,1); wm = ssum(1)
         call gpu_zfmem(g_sum,1)
         call cgpusum2(g_wke,g_sum,mxyp1)
         call gpu_fcopyout(c_loc(ssum),g_sum,1); wke = ssum(1)
         wt = we + wf + wm
         wtot(1) = wt
         wtot(2) = wke
         wtot(3) = 0.0
         wtot(4) = wke + wt
         wtot(5) = we
         wtot(6) = wf
         wtot(7) = wm
         call PPSUM(wtot,work,7)
         wke = wtot(2)
         we = wtot(5)
         wf = wtot(6)
         wm = wtot(7)
         if (kstrt==1) then
            wt = we + wf + wm
            write (*,*) 'Initial Total Field, Kinetic and Total Energies&
     &:'
            write (*,'(3e14.7)') wt, wke, wke + wt
            write (*,*) 'Initial Electrostatic, Transverse Electric and &
     &Magnetic Field Energies:'
            write (*,'(3e14.7)') we, wf, wm
         endif
      endif
      ntime = ntime + 1
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
! energy diagnostic
      call gpu_zfmem(g_sum,1)
      call cgpusum2(g_we,g_sum,kxpd)
      call gpu_fcopyout(c_loc(ssum),g_sum,1); we = ssum(1)
      call gpu_zfmem(g_sum,1)
      call cgpusum2(g_wf,g_sum,kxpd)
      call gpu_fcopyout(c_loc(ssum),g_sum,1); wf = ssum(1)
      call gpu_zfmem(g_sum,1)
      call cgpusum2(g_wm,g_sum,kxpd)
      call gpu_fcopyout(c_loc(ssum),g_sum,1); wm = ssum(1)
      call gpu_zfmem(g_sum,1)
      call cgpusum2(g_wke,g_sum,mxyp1)
      call gpu_fcopyout(c_loc(ssum),g_sum,1); wke = ssum(1)
      wt = we + wf + wm
      wtot(1) = wt
      wtot(2) = wke
      wtot(3) = 0.0
      wtot(4) = wke + wt
      wtot(5) = we
      wtot(6) = wf
      wtot(7) = wm
      call PPSUM(wtot,work,7)
      wke = wtot(2)
      we = wtot(5)
      wf = wtot(6)
      wm = wtot(7)
      if (kstrt==1) then
         write (*,*) 'ntime, relativity = ', ntime, relativity
         write (*,*) 'MPI nodes nvp = ', nvp, ', GPUs per host = ', ndev
         wt = we + wf + wm
         write (*,*) 'Final Total Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') wt, wke, wke + wt
         write (*,*) 'Final Electrostatic, Transverse Electric and Magne&
     &tic Field Energies:'
         write (*,'(3e14.7)') we, wf, wm
         write (*,*)
!
         write (*,*) 'deposit time = ', tdpost
         write (*,*) 'current deposit time = ', tdjpost
         tdpost = tdpost + tdjpost
         write (*,*) 'total deposit time = ', tdpost
         write (*,*) 'guard time = ', tguard
         write (*,*) 'solver time = ', tfield
         write (*,*) 'fft times = ', sum(tfft), tfft
         write (*,*) 'push time = ', tpush
         write (*,*) 'move time = ', tmov
         write (*,*) 'sort time = ', tsort
         tfield = tfield + tguard + sum(tfft)
         write (*,*) 'total solver time = ', tfield
         time = tdpost + tpush + tsort + tmov
         write (*,*) 'total particle time = ', time
         wt = time + tfield
         write (*,*) 'total time = ', wt
         write (*,*)
!
         wt = 1.0e+09/(real(nloop)*real(np))
         write (*,*) 'Push Time (nsec) = ', tpush*wt
         write (*,*) 'Deposit Time (nsec) = ', tdpost*wt
         write (*,*) 'Sort Time (nsec) = ', tsort*wt
         write (*,*) 'Move Time (nsec) = ', tmov*wt
         write (*,*) 'Total Particle Time (nsec) = ', time*wt
      endif
!
! close down NVIDIA fft
      call gpupfft2cudel()
      call gpupfft2rrcudel()
! deallocate memory on GPU
      call gpu_deallocate(g_irc,irc(1))
      call gpu_deallocate(g_scs,irc(1))
      call gpu_deallocate(g_brm,irc(1))
      call gpu_deallocate(g_bsm,irc(1))
      call gpu_deallocate(g_nclr,irc(1))
      call gpu_deallocate(g_ncll,irc(1))
      call gpu_deallocate(g_ihole,irc(1))
      call gpu_deallocate(g_ncl,irc(1))
      call gpu_deallocate(g_sbufr,irc(1))
      call gpu_deallocate(g_sbufl,irc(1))
      call gpu_deallocate(g_kpic,irc(1))
      call gpu_deallocate(g_ppbuff,irc(1))
      call gpu_deallocate(g_ppart,irc(1))
      call gpu_deallocate(g_sum,irc(1))
      call gpu_deallocate(g_wm,irc(1))
      call gpu_deallocate(g_wf,irc(1))
      call gpu_deallocate(g_we,irc(1))
      call gpu_deallocate(g_wke,irc(1))
      call gpu_deallocate(g_bxyzt,irc(1))
      call gpu_deallocate(g_exyzt,irc(1))
      call gpu_deallocate(g_hxyzt,irc(1))
      call gpu_deallocate(g_fxyzt,irc(1))
      call gpu_deallocate(g_hxyz,irc(1))
      call gpu_deallocate(g_fxyz,irc(1))
      call gpu_deallocate(g_cut,irc(1))
      call gpu_deallocate(g_qt,irc(1))
      call gpu_deallocate(g_cu,irc(1))
      call gpu_deallocate(g_q,irc(1))
      call gpu_deallocate(g_sct,irc(1))
      call gpu_deallocate(g_mixup,irc(1))
      call gpu_deallocate(g_ffct,irc(1))
      call gpu_deallocate(g_bxyze,irc(1))
      call gpu_deallocate(g_fxyze,irc(1))
      call gpu_deallocate(g_cue,irc(1))
      call gpu_deallocate(g_qe,irc(1))
! deallocate host page-locked memory
      call hpl_deallocate(c_loc(scs(1)),irc(1)); nullify(scs)
      call hpl_deallocate(c_loc(scr(1)),irc(1)); nullify(scr)
      call hpl_deallocate(c_loc(sbufl(1,1)),irc(1)); nullify(sbufl)
      call hpl_deallocate(c_loc(sbufr(1,1)),irc(1)); nullify(sbufr)
      call hpl_deallocate(c_loc(rbufl(1,1)),irc(1)); nullify(rbufl)
      call hpl_deallocate(c_loc(rbufr(1,1)),irc(1)); nullify(rbufr)
      call hpl_deallocate(c_loc(bsm(1,1)),irc(1)); nullify(bsm)
      call hpl_deallocate(c_loc(brm(1,1)),irc(1)); nullify(brm)
 3000 continue
!
! delete asynchronous streams
      call gpu_delstream(3)
      call gpu_delstream(2)
      call gpu_delstream(1)
! close down GPU
      call end_cu()
! close down MPI
      call PPEXIT()
!
      stop
      end program
!
! Only used by Fortran90, not needed in Fortran2003
      subroutine getfcptr(cref,carray,nx) 
      end subroutine
!
      subroutine getf2cptr(cref,carray,nx,ny) 
      end subroutine
!
      subroutine getc2cptr(cref,carray,nx,ny) 
      end subroutine
