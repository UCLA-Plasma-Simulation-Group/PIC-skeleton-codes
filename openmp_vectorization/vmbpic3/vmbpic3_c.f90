!-----------------------------------------------------------------------
! Skeleton 3D Electromagnetic OpenMP/Vector PIC code
! written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE
      program vmbpic3
      use avx512flib3_h
! #include "vmbpush3.h"
! #include "omplib.h"
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
      real, parameter :: tend = 10.0, dt = 0.035, qme = -1.0
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real, parameter :: vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real, parameter :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
! ax/ay/az = smoothed particle size in x/y/z direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ay = .912871, az = .912871, ci = 0.1
! idimp = number of particle coordinates = 6
! ipbc = particle boundary condition: 1 = periodic
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 6, ipbc = 1, relativity = 1
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
! mx/my/mz = number of grids in x/y/z in sorting tiles
      integer :: mx = 8, my = 8, mz = 8
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! declare scalars for standard code
      integer :: np, nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh
      integer :: nxyzh, nxhyz, mx1, my1, mz1, mxyz1
      integer :: ntime, nloop, isign, lvect
      integer :: irc = 0
      real :: qbme, affp, dth
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, ntmax, npbmx
      integer :: nvp
!
! declare arrays for standard code:
! part = original particle array
      real, dimension(:,:), pointer :: part
! qe = electron charge density with guard cells
      real, dimension(:,:,:), pointer :: qe
! cue = electron current density with guard cells
! fxyze/bxyze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:,:,:), pointer :: cue, fxyze, bxyze
! exyz/bxyz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:,:,:), pointer :: exyz, bxyz
! ffc = form factor array for poisson solver
      complex, dimension(:,:,:), pointer :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct
!
! declare arrays for OpenMP (tiled) code:
! ppartt = tiled particle array
! ppbuff = buffer array for reordering tiled particle array
      real, dimension(:,:,:), pointer :: ppartt, ppbuff
! kpic = number of particles in each tile
      integer, dimension(:), pointer :: kpic
! ncl = number of particles departing tile in each direction
      integer, dimension(:,:), pointer :: ncl
! ihole = location/destination of each particle departing tile
      integer, dimension(:,:,:), pointer :: ihole
! kp = original location of reordered particle
      integer, dimension(:,:), pointer :: kp
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
      real :: tdjpost = 0.0, tpush = 0.0, tsort = 0.0
      double precision :: dtime
!
      irc = 0
! nvp = number of shared memory nodes (0=default)
      nvp = 0
!     write (*,*) 'enter number of nodes:'
!     read (5,*) nvp
! initialize for shared memory parallel processing
      call cinit_omp(nvp)
!
! initialize scalars for standard code
! np = total number of particles in simulation
! nx/ny/nz = number of grid points in x/y direction
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
      dth = 0.0
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
      call avx512_f4allocate(cue,ndim,nxe,nye,nze,irc)
      call avx512_f4allocate(fxyze,ndim,nxe,nye,nze,irc)
      call avx512_f4allocate(bxyze,ndim,nxe,nye,nze,irc)
      call avx512_c4allocate(exyz,ndim,nxeh,nye,nze,irc)
      call avx512_c4allocate(bxyz,ndim,nxeh,nye,nze,irc)
      call avx512_c3allocate(ffc,nxh,nyh,nzh,irc)
      if (irc /= 0) then
         write (*,*) 'aligned field allocation error: irc = ', irc
      endif
!
! prepare fft tables
      call cwfft3rinit(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
! calculate form factors
      isign = 0
      call cvmpois33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh, &
     &nye,nze,nxh,nyh,nzh)
! initialize electrons
      call cdistr3(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,idimp,np,nx,&
     &ny,nz,ipbc)
!
! find number of particles in each of mx, my, mz, tiles:
! updates kpic, nppmx
      call cdblkp3l(part,kpic,nppmx,idimp,np,mx,my,mz,mx1,my1,mxyz1,irc)
      if (irc /= 0) then
         write (*,*) 'cdblkp3l error, irc=', irc
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
      call cppmovin3ltp(part,ppartt,kpic,kp,nppmx0,idimp,np,mx,my,mz,mx1&
     &,my1,mxyz1,irc)
      if (irc /= 0) then
         write (*,*) 'cppmovin3ltp overflow error, irc=', irc
         stop
      endif
! sanity check
      call cppcheck3lt(ppartt,kpic,idimp,nppmx0,nx,ny,nz,mx,my,mz,mx1,  &
     &my1,mz1,irc)
      if (irc /= 0) then
         write (*,*) 'cppcheck3lt error: irc=', irc
         stop
      endif
!
! initialize transverse electromagnetic fields
      call cset_cvzero3(exyz,nx,ny,nz,ndim,nxeh,nye,nze)
      call cset_cvzero3(bxyz,nx,ny,nz,ndim,nxeh,nye,nze)
!
      if (dt > 0.37*ci) then
         write (*,*) 'Warning: Courant condition may be exceeded!'
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     write (*,*) 'ntime = ', ntime
!
! deposit current with OpenMP: 
      call dtimer(dtime,itime,-1)
      call cset_vzero3(cue,mx,my,mz,ndim,nxe,nye,nze,mx1,my1,mxyz1)
      if (relativity==1) then
! updates ppart, cue
!        call cgrjppost3lt(ppartt,cue,kpic,qme,dth,ci,nppmx0,idimp,nx,ny&
!    &,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
!        call cvgrjppost3lt(ppartt,cue,kpic,qme,dth,ci,nppmx0,idimp,nx, &
!    &ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
! updates ppart, cue, ncl, ihole, irc
!        call cgrjppostf3lt(ppartt,cue,kpic,ncl,ihole,qme,dth,ci,nppmx0,&
!    &idimp,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,irc)
         call cvgrjppostf3lt(ppartt,cue,kpic,ncl,ihole,qme,dth,ci,nppmx0&
     &,idimp,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,irc)
      else
! updates ppart, cue
!        call cgjppost3lt(ppartt,cue,kpic,qme,dth,nppmx0,idimp,nx,ny,nz,&
!    &mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
!        call cvgjppost3lt(ppartt,cue,kpic,qme,dth,nppmx0,idimp,nx,ny,nz&
!    &,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
! updates ppart, cue, ncl, ihole, irc
!        call cgjppostf3lt(ppartt,cue,kpic,ncl,ihole,qme,dth,nppmx0,    &
!    &idimp,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,irc)
         call cvgjppostf3lt(ppartt,cue,kpic,ncl,ihole,qme,dth,nppmx0,   &
     &idimp,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,irc)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdjpost = tdjpost + time
      if (irc /= 0) then
         if (relativity==1) then
            write (*,*) 'cvgjppostf3lt error: irc=', irc
         else
            write (*,*) 'cvgrjppostf3lt error: irc=', irc
         endif
         stop
      endif
!
! reorder particles by tile with OpenMP:
      call dtimer(dtime,itime,-1)
! updates ppartt, ppbuff, kpic, ncl, ihole, and irc
!     call cpporder3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny, &
!    &nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
!     call cvpporder3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny,&
!    &nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
! updates ppartt, ppbuff, kpic, ncl, and irc
!     call cpporderf3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,  &
!    &my1,mz1,npbmx,ntmax,irc)
      call cvpporderf3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1, &
     &my1,mz1,npbmx,ntmax,irc)
!     call cv2pporderf3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,&
!    &my1,mz1,npbmx,ntmax,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tsort = tsort + time
      if (irc /= 0) then
         write (*,*) 'cvpporderf3lt error: ntmax, irc=', ntmax, irc
         stop
      endif
!
! deposit charge with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      call cset_szero3(qe,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1)
!     call cgppost3lt(ppartt,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,nye, &
!    &nze,mx1,my1,mxyz1)
      call cvgppost3lt(ppartt,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,nye,&
     &nze,mx1,my1,mxyz1)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with OpenMP: updates cue, qe
      call dtimer(dtime,itime,-1)
      call cacguard3l(cue,nx,ny,nz,nxe,nye,nze)
      call caguard3l(qe,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with OpenMP: updates qe
      call dtimer(dtime,itime,-1)
      isign = -1
      call cwfft3rvmx(qe,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,   &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! transform current to fourier space with OpenMP: update cue
      call dtimer(dtime,itime,-1)
      isign = -1
      call cwfft3rvm3(cue,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,  &
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of current with OpenMP: updates cue
      call dtimer(dtime,itime,-1)
      call cmcuperp3(cue,nx,ny,nz,nxeh,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate electromagnetic fields in fourier space with OpenMP:
! updates exyz, bxyz, wf, wm
      call dtimer(dtime,itime,-1)
      if (ntime==0) then
         call cvmibpois33(cue,bxyz,ffc,ci,wm,nx,ny,nz,nxeh,nye,nze,nxh, &
     &nyh,nzh)
         wf = 0.0
         dth = 0.5*dt
      else
         call cvmmaxwel3(exyz,bxyz,cue,ffc,ci,dt,wf,wm,nx,ny,nz,nxeh,nye&
     &,nze,nxh,nyh,nzh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate force/charge in fourier space with OpenMP:
! updates fxyze, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call cvmpois33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh, &
     &nye,nze,nxh,nyh,nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! add longitudinal and transverse electric fields with OpenMP:
! updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call cvmemfield3(fxyze,exyz,ffc,isign,nx,ny,nz,nxeh,nye,nze,nxh,nyh&
     &,nzh)
! copy magnetic field with OpenMP: updates bxyze
      isign = -1
      call cvmemfield3(bxyze,bxyz,ffc,isign,nx,ny,nz,nxeh,nye,nze,nxh,nyh&
     &,nzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform electric force to real space with OpenMP: updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call cwfft3rvm3(fxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,&
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! transform magnetic force to real space with OpenMP: updates bxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call cwfft3rvm3(bxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,&
     &nxhyz,nxyzh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with OpenMP: updates fxyze, bxyze
      call dtimer(dtime,itime,-1)
      call ccguard3l(fxyze,nx,ny,nz,nxe,nye,nze)
      call ccguard3l(bxyze,nx,ny,nz,nxe,nye,nze)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! push particles with OpenMP: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
      if (relativity==1) then
! updates ppart, wke
!        call cgrbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,ci,wke,  &
!    &idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
!        call cvgrbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,ci,wke, &
!    &idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
!        call cv2grbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,ci,wke,&
!    &idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
! updates ppart, ncl, ihole, wke, irc
!        call cgrbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt,  &
!    &dth,ci,wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,    &
!    &mxyz1,ntmax,irc)
         call cvgrbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt, &
     &dth,ci,wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,    &
     &mxyz1,ntmax,irc)
!        call cv2grbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt,&
!    &dth,ci,wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,    &
!    &mxyz1,ntmax,irc)
      else
! updates ppart, wke
!        call cgbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,wke,idimp,&
!    &nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
!        call cvgbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,wke,idimp&
!    &,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
!        call cv2gbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,wke,    &
!    &idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
! updates ppart, ncl, ihole, wke, irc
!        call cgbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt,dth&
!    &,wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,    &
!    &ntmax,irc)
         call cvgbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt,  &
     &dth,wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1, &
     &ntmax,irc)
!        call cv2gbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt, &
!    &dth,wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1, &
!    &ntmax,irc)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
      if (irc /= 0) then
         if (relativity==1) then
            write (*,*) 'cvgrbppushf3lt error: irc=', irc
         else
            write (*,*) 'cvgbppushf3lt error: irc=', irc
         endif
         stop
      endif
!
! reorder particles by tile with OpenMP:
      call dtimer(dtime,itime,-1)
! updates ppartt, ppbuff, kpic, ncl, ihole, and irc
!     call cpporder3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny, &
!    &nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
!     call cvpporder3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny,&
!    &nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
! updates ppartt, ppbuff, kpic, ncl, and irc
!     call cpporderf3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,  &
!    &my1,mz1,npbmx,ntmax,irc)
!     call cvpporderf3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1, &
!    &my1,mz1,npbmx,ntmax,irc)
      call cv2pporderf3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,&
     &my1,mz1,npbmx,ntmax,irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tsort = tsort + time
      if (irc /= 0) then
         write (*,*) 'VPPORDERF3LT error: ntmax, irc=', ntmax, irc
         stop
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
      write (*,*)
!
      call avx512_deallocate(ppartt); nullify(ppartt)
      call avx512_deallocate(ppbuff); nullify(ppbuff)
      call avx512_deallocate(ffc); nullify(ffc)
      call avx512_deallocate(bxyz); nullify(bxyz)
      call avx512_deallocate(exyz); nullify(exyz)
      call avx512_deallocate(bxyze); nullify(bxyze)
      call avx512_deallocate(fxyze); nullify(fxyze)
      call avx512_deallocate(cue); nullify(cue)
      call avx512_deallocate(qe); nullify(qe)
!
      stop
      end program
!
! Procedures to create Fortran90 pointers for data allocated in C.
! For details see V. K. Decyk, ACM Fortran Forum, vol. 27, no. 2 (2008).
      subroutine getf3cptr(cref,carray,nx,ny,nz) 
! set reference to C data in 3d real Fortran pointer object 
      implicit none 
      integer :: nx, ny, nz
      real, dimension(nx,ny,nz), target :: carray 
      real, dimension(:,:,:), pointer :: cref 
      cref => carray 
      end subroutine
!
      subroutine getf4cptr(cref,carray,ndim,nx,ny,nz) 
! set reference to C data in 4d real Fortran pointer object 
      implicit none 
      integer :: ndim, nx, ny, nz
      real, dimension(ndim,nx,ny,nz), target :: carray 
      real, dimension(:,:,:,:), pointer :: cref 
      cref => carray 
      end subroutine
!
      subroutine getc3cptr(cref,carray,nx,ny,nz) 
! set reference to C data in 3d complex Fortran pointer object 
      implicit none 
      integer :: nx, ny, nz
      complex, dimension(nx,ny,nz), target :: carray 
      complex, dimension(:,:,:), pointer :: cref 
      cref => carray 
      end subroutine
!
      subroutine getc4cptr(cref,carray,ndim,nx,ny,nz) 
! set reference to C data in vector 3d complex Fortran pointer object 
      implicit none 
      integer :: ndim, nx, ny, nz
      complex, dimension(ndim,nx,ny,nz), target :: carray 
      complex, dimension(:,:,:,:), pointer :: cref 
      cref => carray 
      end subroutine
