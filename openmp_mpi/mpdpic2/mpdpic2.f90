!-----------------------------------------------------------------------
! Skeleton 2-1/2D Darwin MPI/OpenMP PIC code
! written by Viktor K. Decyk, UCLA
      program mpdpic2
      use mpdpush2_h
      use mpplib2       ! use with pplib2.f90
!     use mpplib2_h     ! use with pplib2.f
      use omplib_h
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
! idimp = dimension of phase space = 5
! ipbc = particle boundary condition: 1 = periodic
      integer :: idimp = 5, ipbc = 1
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.4, omy = 0.0, omz = 0.0
! ndc = number of corrections in darwin iteration
      integer :: ndc = 1
! idps = number of partition boundaries
      integer :: idps = 2
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
      real :: zero = 0.0
! sorting tiles, should be less than or equal to 32
      integer :: mx = 16, my = 16
! fraction of extra particles needed for particle management
      real :: xtras = 0.2
! declare scalars for standard code
      integer :: k
      integer :: nx, ny, nxh, nyh, nxe, nye, nxeh, nnxe, nxyh, nxhy
      integer :: mdim, mx1, ntime, nloop, isign, ierr
      real :: qbme, affp, q2m0, wpm, wpmax, wpmin
      double precision :: np
!
! declare scalars for MPI code
      integer :: ntpose = 1
      integer :: nvp, idproc, kstrt, npmax, kxp, kyp, nypmx, nypmn
      integer :: nyp, noff, npp, nps, myp1, mxyp1
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, nbmaxp, ntmaxp, npbmx, irc
      integer :: nvpp
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), pointer :: part
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
! ss = scratch array for WPPFFT2RMN
      complex, dimension(:,:), pointer :: ss
! qt = scalar charge density field array in fourier space
      complex, dimension(:,:), pointer :: qt
! cut = vector current density in fourier space
! dcut = vector acceleration density in fourier space
! cur = vector transverse electric field in fourier space
! amut = tensor momentum flux in fourier space
      complex, dimension(:,:,:), pointer :: cut, dcut, cur, amut
! exyt = vector total electric field in fourier space
! fxyt = vector longitudinal electric field in fourier space
! bxyt = vector magnetic field in fourier space
      complex, dimension(:,:,:), pointer :: exyt, fxyt, bxyt
! ffc, ffe = form factor arrays for poisson solvers
      complex, dimension(:,:), pointer :: ffc, ffe
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct
      double precision, dimension(7) :: wtot, work
!
! declare arrays for MPI code:
! bs/br = complex send/receive buffers for data transpose
      complex, dimension(:,:,:), pointer :: bs, br
! sbufl/sbufr = particle buffers sent to nearby processors
! rbufl/rbufr = particle buffers received from nearby processors
      real, dimension(:,:), pointer :: sbufl, sbufr, rbufl, rbufr
! edges(1:2) = lower:upper y boundaries of particle partition
      real, dimension(:), pointer  :: edges
! scr = guard cell buffer received from nearby processors
      real, dimension(:), pointer  :: scr
!
! declare arrays for OpenMP code:
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
      integer, dimension(:,:), pointer :: ncll, nclr, mcll, mclr
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0
      real :: tdjpost = 0.0, tdcjpost = 0.0, tpush = 0.0, tsort = 0.0
      real :: tmov = 0.0
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
      np =  dble(npx)*dble(npy)
! nx/ny = number of grid points in x/y direction
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = max(1,ny/2)
      nxe = nx + 2; nye = ny + 2; nxeh = nxe/2; nnxe = ndim*nxe
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny)
! mx1 = number of tiles in x direction
      mx1 = (nx - 1)/mx + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
! mdim = dimension of amu array
      mdim = 2*ndim - 2
      qbme = qme
      affp = dble(nx)*dble(ny)/np
!      
! nvp = number of distributed memory nodes
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
!
! initialize additional scalars for MPI code
! kxp = number of complex grids in each field partition in x direction
      kxp = (nxh - 1)/nvp + 1
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvp + 1
! npmax = maximum number of electrons in each partition
      npmax = (np/nvp)*1.25
! myp1 = number of tiles in y direction
      myp1 = (nyp - 1)/my + 1; mxyp1 = mx1*myp1
!
! allocate and initialize data for standard code
      allocate(part(idimp,npmax))
      allocate(qe(nxe,nypmx),fxyze(ndim,nxe,nypmx))
      allocate(exyze(ndim,nxe,nypmx),cue(ndim,nxe,nypmx))
      allocate(dcu(ndim,nxe,nypmx),cus(ndim,nxe,nypmx))
      allocate(amu(mdim,nxe,nypmx),bxyze(ndim,nxe,nypmx))
      allocate(qt(nye,kxp),fxyt(ndim,nye,kxp),exyt(ndim,nye,kxp))
      allocate(cut(ndim,nye,kxp),dcut(ndim,nye,kxp),cur(ndim,nye,kxp))
      allocate(amut(mdim,nye,kxp),bxyt(ndim,nye,kxp))
      allocate(ffc(nyh,kxp),ffe(nyh,kxp),mixup(nxhy),sct(nxyh))
      allocate(ss(mdim*nxeh,nypmx))
      allocate(kpic(mxyp1))
!
! allocate and initialize data for MPI code
      allocate(bs(mdim,kxp,kyp),br(mdim,kxp,kyp))
      allocate(scr(mdim*nxe))
!
! prepare fft tables
      call WPFFT2RINIT(mixup,sct,indx,indy,nxhy,nxyh)
! calculate form factor: ffc
      isign = 0
      call MPPOIS23(qt,fxyt,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nye,kxp,&
     &nyh)
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
! find number of particles in each of mx, my tiles: updates kpic, nppmx
      call PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,mx1,    &
     &mxyp1,irc)
      if (irc /= 0) then
         write (*,*) 'PPDBLKP2L error, irc=', irc
         call PPABORT()
         stop
      endif
! allocate vector particle data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmaxp = xtras*nppmx
      npbmx = xtras*nppmx
      nbmaxp = 0.25*mx1*npbmx
      allocate(sbufl(idimp,nbmaxp),sbufr(idimp,nbmaxp))
      allocate(rbufl(idimp,nbmaxp),rbufr(idimp,nbmaxp))
      allocate(ppart(idimp,nppmx0,mxyp1))
      allocate(ppbuff(idimp,npbmx,mxyp1))
      allocate(ncl(8,mxyp1))
      allocate(iholep(2,ntmaxp+1,mxyp1))
      allocate(ncll(3,mx1),nclr(3,mx1),mcll(3,mx1),mclr(3,mx1))
!
! copy ordered particle data for OpenMP
      call PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx0,idimp,npmax,mx,my,&
     &mx1,mxyp1,irc)
      if (irc /= 0) then
         write (*,*) kstrt, 'PPPMOVIN2L overflow error, irc=', irc
         call PPABORT()
         stop
      endif
! sanity check
      call PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx0,nx,mx,my,mx1,myp1&
     &,irc)
      if (irc /= 0) then
         write (*,*) kstrt, 'PPPCHECK2L error: irc=', irc
         call PPABORT()
         stop
      endif
!
! find maximum and minimum initial electron density
      qe = 0.0
      call PPGPPOST2L(ppart,qe,kpic,noff,qme,idimp,nppmx0,mx,my,nxe,    &
     &nypmx,mx1,mxyp1)
      call PPAGUARD2XL(qe,nyp,nx,nxe,nypmx)
      call PPNAGUARD2L(qe,scr,nyp,nx,kstrt,nvp,nxe,nypmx)
      call PPFWPMINMX2(qe,nyp,qbme,wpmax,wpmin,nx,nxe,nypmx)
      wtot(1) = wpmax
      wtot(2) = -wpmin
      call PPDMAX(wtot,work,2)
      wpmax = wtot(1)
      wpmin = -wtot(2)
      wpm = 0.5*(wpmax + wpmin)*affp
! accelerate convergence: update wpm
      if (wpm <= 10.0) wpm = 0.75*wpm
      if (kstrt==1) write (*,*) 'wpm=',wpm
      q2m0 = wpm/affp
! calculate form factor: ffe
      isign = 0
      call MPPEPOISP23(dcut,cur,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,   &
     &kstrt,nye,kxp,nyh)
!
! initialize transverse electric field
      cus = 0.0
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     if (kstrt==1) write (*,*) 'ntime = ', ntime
!
! deposit current with OpenMP: updates cue
      call dtimer(dtime,itime,-1)
      cue = 0.0
      call PPGJPPOST2L(ppart,cue,kpic,noff,qme,zero,nppmx0,idimp,nx,ny, &
     &mx,my,nxe,nypmx,mx1,mxyp1,ipbc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdjpost = tdjpost + time
!
! deposit charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call PPGPPOST2L(ppart,qe,kpic,noff,qme,idimp,nppmx0,mx,my,nxe,    &
     &nypmx,mx1,mxyp1)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with OpenMP: updates cue, qe
      call dtimer(dtime,itime,-1)
      call PPAGUARD2XL(qe,nyp,nx,nxe,nypmx)
      call PPNAGUARD2L(qe,scr,nyp,nx,kstrt,nvp,nxe,nypmx)
      call PPACGUARD2XL(cue,nyp,nx,ndim,nxe,nypmx)
      call PPNACGUARD2L(cue,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with OpenMP: updates qt, modifies qe
      call dtimer(dtime,itime,-1)
      isign = -1
      call WPPFFT2RM(qe,qt,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,  &
     &kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! calculate longitudinal force/charge in fourier space with OpenMP:
! updates fxyt, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call MPPOIS23(qt,fxyt,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nye,kxp,&
     &nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform longitudinal electric force to real space with OpenMP:
! updates fxyze, modifies fxyt
      call dtimer(dtime,itime,-1)
      isign = 1
      call WPPFFT2RM3(fxyze,fxyt,bs,br,isign,ntpose,mixup,sct,ttp,indx, &
     &indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! transform current to fourier space with OpenMP: updates cut
! modifies cue
      call dtimer(dtime,itime,-1)
      isign = -1
      call WPPFFT2RM3(cue,cut,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy&
     &,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! take transverse part of current with OpenMP: updates cut
      call dtimer(dtime,itime,-1)
      call MPPCUPERP2(cut,nx,ny,kstrt,nye,kxp)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate magnetic field in fourier space with standard procedure:
! updates bxyt, wm
      call dtimer(dtime,itime,-1)
      call MPPBBPOISP23(cut,bxyt,ffc,ci,wm,nx,ny,kstrt,nye,kxp,nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform magnetic field to real space with OpenMP: updates bxyze
! modifies bxyt
      call dtimer(dtime,itime,-1)
      isign = 1
      call WPPFFT2RM3(bxyze,bxyt,bs,br,isign,ntpose,mixup,sct,ttp,indx, &
     &indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! add constant to magnetic field with standard procedure: updates bxyze
      call dtimer(dtime,itime,-1)
      call PPBADDEXT2(bxyze,nyp,omx,omy,omz,nx,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! copy guard cells with OpenMP: updates fxyze, bxyze
      call dtimer(dtime,itime,-1)
      call PPNCGUARD2L(fxyze,nyp,kstrt,nvp,nnxe,nypmx)
      call PPCGUARD2XL(fxyze,nyp,nx,ndim,nxe,nypmx)
      call PPNCGUARD2L(bxyze,nyp,kstrt,nvp,nnxe,nypmx)
      call PPCGUARD2XL(bxyze,nyp,nx,ndim,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and old transverse electric fields with standard
! procedure: updates exyze
      call dtimer(dtime,itime,-1)
      call PPADDVRFIELD2(exyze,cus,fxyze,ndim,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! deposit electron acceleration density and momentum flux with standard
! procedure: updates dcu, amu
      call dtimer(dtime,itime,-1)
      dcu = 0.0; amu = 0.0
      call PPGDJPPOST2L(ppart,exyze,bxyze,dcu,amu,kpic,noff,nyp,qme,qbme&
     &,dt,idimp,nppmx0,nx,mx,my,nxe,nypmx,mx1,mxyp1)
! add old scaled electric field with standard procedure: updates dcu
      call PPASCFGUARD2L(dcu,cus,nyp,q2m0,nx,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdcjpost = tdcjpost + time
!
! add guard cells with standard procedure: updates dcu, amu
      call dtimer(dtime,itime,-1)
      call PPACGUARD2XL(dcu,nyp,nx,ndim,nxe,nypmx)
      call PPNACGUARD2L(dcu,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx)
      call PPACGUARD2XL(amu,nyp,nx,mdim,nxe,nypmx)
      call PPNACGUARD2L(amu,scr,nyp,nx,mdim,kstrt,nvp,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform acceleration density and momentum flux to fourier space
! with OpenMP: updates dcut, amut, modifies dcu, amu
      call dtimer(dtime,itime,-1)
      isign = -1
      call WPPFFT2RM3(dcu,dcut,bs,br,isign,ntpose,mixup,sct,ttp,indx,   &
     &indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      tfft(2) = tfft(2) + ttp
      call WPPFFT2RMN(amu,amut,bs,br,ss,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,mdim,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! take transverse part of time derivative of current with standard
! procedure: updates dcut
      call dtimer(dtime,itime,-1)
      call MPPADCUPERP23(dcut,amut,nx,ny,kstrt,nye,kxp)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate transverse electric field with standard procedure:
! updates cur, wf
      call dtimer(dtime,itime,-1)
      isign = -1
      call MPPEPOISP23(dcut,cur,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,   &
     &kstrt,nye,kxp,nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform transverse electric field to real space with OpenMP:
! updates cus, modifies cur
      call dtimer(dtime,itime,-1)
      isign = 1
      call WPPFFT2RM3(cus,cur,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy&
     &,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! copy guard cells with OpenMP: updates cus
      call dtimer(dtime,itime,-1)
      call PPNCGUARD2L(cus,nyp,kstrt,nvp,nnxe,nypmx)
      call PPCGUARD2XL(cus,nyp,nx,ndim,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and transverse electric fields with standard
! procedure: exyze = cus + fxyze, updates exyze
! cus needs to be retained for next time step
      call dtimer(dtime,itime,-1)
      call PPADDVRFIELD2(exyze,cus,fxyze,ndim,nxe,nypmx)
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
      call PPGDCJPPOST2L(ppart,exyze,bxyze,cue,dcu,amu,kpic,noff,nyp,qme&
     &,qbme,dt,idimp,nppmx0,nx,mx,my,nxe,nypmx,mx1,mxyp1)
! add scaled electric field with standard procedure: updates dcu
      call PPASCFGUARD2L(dcu,cus,nyp,q2m0,nx,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdcjpost = tdcjpost + time
!
! add guard cells for current, acceleration density, and momentum flux
! with OpenMP: updates cue, dcu, amu
      call dtimer(dtime,itime,-1)
      call PPACGUARD2XL(cue,nyp,nx,ndim,nxe,nypmx)
      call PPNACGUARD2L(cue,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx)
      call PPACGUARD2XL(dcu,nyp,nx,ndim,nxe,nypmx)
      call PPNACGUARD2L(dcu,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx)
      call PPACGUARD2XL(amu,nyp,nx,mdim,nxe,nypmx)
      call PPNACGUARD2L(amu,scr,nyp,nx,mdim,kstrt,nvp,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform current to fourier space with OpenMP: updates cut
! modifies cue
      call dtimer(dtime,itime,-1)
      isign = -1
      call WPPFFT2RM3(cue,cut,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy&
     &,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! take transverse part of current with OpenMP: updates cut
      call dtimer(dtime,itime,-1)
      call MPPCUPERP2(cut,nx,ny,kstrt,nye,kxp)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate magnetic field in fourier space with standard procedure:
! updates bxyt, wm
      call dtimer(dtime,itime,-1)
      call MPPBBPOISP23(cut,bxyt,ffc,ci,wm,nx,ny,kstrt,nye,kxp,nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform magnetic field to real space with OpenMP: updates bxyze
! modifies bxyt
      call dtimer(dtime,itime,-1)
      isign = 1
      call WPPFFT2RM3(bxyze,bxyt,bs,br,isign,ntpose,mixup,sct,ttp,indx, &
     &indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! add constant to magnetic field with standard procedure: updates bxyze
      call dtimer(dtime,itime,-1)
      call PPBADDEXT2(bxyze,nyp,omx,omy,omz,nx,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform acceleration density and momentum flux to fourier space
! with OpenMP: updates dcut, amut, modifies dcu, amu
      call dtimer(dtime,itime,-1)
      isign = -1
      call WPPFFT2RM3(dcu,dcut,bs,br,isign,ntpose,mixup,sct,ttp,indx,   &
     &indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      tfft(2) = tfft(2) + ttp
      call WPPFFT2RMN(amu,amut,bs,br,ss,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,mdim,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! take transverse part of time derivative of current with standard
! procedure: updates dcut
      call dtimer(dtime,itime,-1)
      call MPPADCUPERP23(dcut,amut,nx,ny,kstrt,nye,kxp)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate transverse electric field with standard procedure:
! updates cur, wf
      call dtimer(dtime,itime,-1)
      isign = -1
      call MPPEPOISP23(dcut,cur,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,   &
     &kstrt,nye,kxp,nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform transverse electric field to real space with OpenMP:
! updates cus, modifies cur
      call dtimer(dtime,itime,-1)
      isign = 1
      call WPPFFT2RM3(cus,cur,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy&
     &,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! copy guard cells with OpenMP: updates bxyze, cus
      call dtimer(dtime,itime,-1)
      call PPNCGUARD2L(bxyze,nyp,kstrt,nvp,nnxe,nypmx)
      call PPCGUARD2XL(bxyze,nyp,nx,ndim,nxe,nypmx)
      call PPNCGUARD2L(cus,nyp,kstrt,nvp,nnxe,nypmx)
      call PPCGUARD2XL(cus,nyp,nx,ndim,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! add longitudinal and transverse electric fields with standard
! procedure: exyze = cus + fxyze, updates exyze
! cus needs to be retained for next time step
      call dtimer(dtime,itime,-1)
      call PPADDVRFIELD2(exyze,cus,fxyze,ndim,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
      enddo
!
! push particles with OpenMP:
      call dtimer(dtime,itime,-1)
      wke = 0.0
! updates ppart and wke
!     call PPGBPPUSH23L(ppart,exyze,bxyze,kpic,noff,nyp,qbme,dt,dt,wke, &
!    &idimp,nppmx0,nx,ny,mx,my,nxe,nypmx,mx1,mxyp1,ipbc)
! updates ppart, ncl, iholep, ek, irc
      call PPGBPPUSHF23L(ppart,exyze,bxyze,kpic,ncl,iholep,noff,nyp,qbme&
     &,dt,dt,wke,idimp,nppmx0,nx,ny,mx,my,nxe,nypmx,mx1,mxyp1,ntmaxp,irc&
     &)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
      if (irc /= 0) then
         write (*,*) kstrt, 'PPGBPPUSHF23L error: irc=', irc
         call PPABORT()
         stop
      endif
!
! reorder particles by tile with OpenMP
! first part of particle reorder on x and y cell with mx, my tiles:
      call dtimer(dtime,itime,-1)
! updates ppart, ppbuff, sbufl, sbufr, ncl, iholep, ncll, nclr, irc
!     call PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,iholep,ncll,   &
!    &nclr,noff,nyp,idimp,nppmx0,nx,ny,mx,my,mx1,myp1,npbmx,ntmaxp,     &
!    &nbmaxp,irc)
! updates ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
      call PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,iholep,ncll,nclr,  &
     &idimp,nppmx0,mx1,myp1,npbmx,ntmaxp,nbmaxp,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tsort = tsort + time
      if (irc /= 0) then
         write (*,*) kstrt,'PPPORDERF2LA error: ntmaxp, irc=',ntmaxp,irc
         call PPABORT()
         stop
      endif
! move particles into appropriate spatial regions:
! updates rbufr, rbufl, mcll, mclr
      call dtimer(dtime,itime,-1)
      call PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,  &
     &nvp,idimp,nbmaxp,mx1)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tmov = tmov + time
! second part of particle reorder on x and y cell with mx, my tiles:
! updates ppart, kpic
      call dtimer(dtime,itime,-1)
      call PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,iholep,mcll,   &
     &mclr,idimp,nppmx0,mx1,myp1,npbmx,ntmaxp,nbmaxp,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tsort = tsort + time
      if (irc /= 0) then
         write (*,*) kstrt,'PPPORDER2LB error: nppmx0, irc=',nppmx0,irc
         call PPABORT()
         stop
      endif
!
! energy diagnostic
      wt = we + wm
      wtot(1) = wt
      wtot(2) = wke
      wtot(3) = 0.0
      wtot(4) = wke + wt
      wtot(5) = we
      wtot(6) = wf
      wtot(7) = wm
      call PPDSUM(wtot,work,7)
      wke = wtot(2)
      we = wtot(5)
      wf = wtot(6)
      wm = wtot(7)
      if (ntime==0) then
         if (kstrt.eq.1) then
            wt = we + wm
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
      if (kstrt.eq.1) then
         write (*,*) 'ntime, ndc = ', ntime, ndc
         write (*,*) 'MPI nodes nvp = ', nvp
         wt = we + wm
         write (*,*) 'Final Total Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') wt, wke, wke + wt
         write (*,*) 'Final Electrostatic, Transverse Electric and Magne&
     &tic Field Energies:'
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
      stop
      end program
