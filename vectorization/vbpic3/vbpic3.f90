!-----------------------------------------------------------------------
! Skeleton 3D Electromagnetic Vector PIC code
! written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE
      program vbpic3
      use avx512flib3_h
      use kncbpush3_h
      use vbpush3_h
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
! sortime = number of time steps between standard electron sorting
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 6, ipbc = 1, sortime = 20, relativity = 1
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
! kvec = (1,2) = run (autovector,KNC) version
      integer :: kvec = 1
!
! declare scalars for standard code
      integer :: np, nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh
      integer :: nxyzh, nxhyz, npe, ny1, nyz1, ntime, nloop, isign
      integer :: irc = 0
      real :: qbme, affp, dth
!
! declare arrays for standard code:
! partt, partt2 = transposed particle arrays
      real, dimension(:,:), pointer :: partt, partt2, tpartt
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
! npic = scratch array for reordering particles
      integer, dimension(:), pointer :: npic
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
! nx/ny/nz = number of grid points in x/y direction
      np = npx*npy*npz; nx = 2**indx; ny = 2**indy; nz = 2**indz
      nxh = nx/2; nyh = max(1,ny/2); nzh = max(1,nz/2)
      nxe = nx + 2; nye = ny + 1; nze = nz + 1; nxeh = nxe/2
      nxyzh = max(nx,ny,nz)/2; nxhyz = max(nxh,ny,nz)
      ny1 = ny + 1; nyz1 = ny1*(nz + 1)
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx)*real(ny)*real(nz)/real(np)
      dth = 0.0
!
! allocate data for standard code
      allocate(mixup(nxhyz),sct(nxyzh))
!
! align memory for KNC
      npe = 16*((np - 1)/16 + 1)
      nxe = 16*((nxe - 1)/16 + 1)
      nxeh = nxe/2
      call avx512_f2allocate(partt,npe,idimp,irc)
      if (sortime > 0) then
         call avx512_f2allocate(partt2,npe,idimp,irc)
      endif
      call avx512_f3allocate(qe,nxe,nye,nze,irc)
      call avx512_f4allocate(cue,ndim,nxe,nye,nze,irc)
      call avx512_f4allocate(fxyze,ndim,nxe,nye,nze,irc)
      call avx512_f4allocate(bxyze,ndim,nxe,nye,nze,irc)
      call avx512_c4allocate(exyz,ndim,nxeh,nye,nze,irc)
      call avx512_c4allocate(bxyz,ndim,nxeh,nye,nze,irc)
      call avx512_c3allocate(ffc,nxh,nyh,nzh,irc)
      call avx512_i1allocate(npic,nyz1,irc)
      if (irc /= 0) then
         write (*,*) 'aligned allocation error: irc = ', irc
      endif
!
! prepare fft tables
      call WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
! calculate form factors
      isign = 0
      call VPOIS33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,nye&
     &,nze,nxh,nyh,nzh)
! initialize electrons
      call DISTR3T(partt,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,idimp,npe, &
     &nx,ny,nz,ipbc)
!
! initialize transverse electromagnetic fields
      exyz = cmplx(0.0,0.0)
      bxyz = cmplx(0.0,0.0)
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
! deposit current with standard procedure: updates part, cue
      call dtimer(dtime,itime,-1)
      cue = 0.0
      if (relativity==1) then
         if (kvec==1) then
!           call GRJPOST3LT(partt,cue,qme,dth,ci,np,npe,idimp,nx,ny,nz, &
!    &nxe,nye,nze,ipbc)
            call VGRJPOST3LT(partt,cue,qme,dth,ci,np,npe,idimp,nx,ny,nz,&
     &nxe,nye,nze,ipbc)
! KNC function
         else if (kvec==2) then
            call ckncgrjpost3lt(partt,cue,qme,dth,ci,np,npe,idimp,nx,ny,&
     &nz,nxe,nye,nze,ipbc)
         endif
      else
         if (kvec==1) then
!           call GJPOST3LT(partt,cue,qme,dth,np,npe,idimp,nx,ny,nz,nxe, &
!    &nye,nze,ipbc)
            call VGJPOST3LT(partt,cue,qme,dth,np,npe,idimp,nx,ny,nz,nxe,&
     &nye,nze,ipbc)
! KNC function
         else if (kvec==2) then
            call ckncgjpost3lt(partt,cue,qme,dth,np,npe,idimp,nx,ny,nz, &
     &nxe,nye,nze,ipbc)
         endif
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdjpost = tdjpost + time
!
! deposit charge with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      if (kvec==1) then
         call GPOST3LT(partt,qe,qme,np,npe,idimp,nxe,nye,nze)
!        call VGPOST3LT(partt,qe,qme,np,npe,idimp,nxe,nye,nze)
! KNC function
      else if (kvec==2) then
         call cknc2gpost3lt(partt,qe,qme,np,npe,idimp,nxe,nye,nze)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with standard procedure: updates cue, qe
      call dtimer(dtime,itime,-1)
      if (kvec==1) then
         call ACGUARD3L(cue,nx,ny,nz,nxe,nye,nze)
         call AGUARD3L(qe,nx,ny,nz,nxe,nye,nze)
! KNC function
      else if (kvec==2) then
         call ckncacguard3l(cue,nx,ny,nz,nxe,nye,nze)
         call ckncaguard3l(qe,nx,ny,nz,nxe,nye,nze)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      isign = -1
      if (kvec==1) then
         call WFFT3RVX(qe,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,  &
     &nxhyz,nxyzh)
! KNC function
      else if (kvec==2) then
         call ckncwfft3rvx(qe,isign,mixup,sct,indx,indy,indz,nxeh,nye,  &
     &nze,nxhyz,nxyzh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! transform current to fourier space with standard procedure: update cue
      call dtimer(dtime,itime,-1)
      isign = -1
      if (kvec==1) then
         call WFFT3RV3(cue,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze, &
     &nxhyz,nxyzh)
! KNC function
      else if (kvec==2) then
         call ckncwfft3rv3(cue,isign,mixup,sct,indx,indy,indz,nxeh,nye, &
      &nze,nxhyz,nxyzh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! take transverse part of current with standard procedure: updates cue
      call dtimer(dtime,itime,-1)
      if (kvec==1) then
         call CUPERP3(cue,nx,ny,nz,nxeh,nye,nze)
! KNC function
      else if (kvec==2) then
         call cknccuperp3(cue,nx,ny,nz,nxeh,nye,nze)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time 
!
! calculate electromagnetic fields in fourier space with standard
! procedure: updates exyz, bxyz, wf, wm
      call dtimer(dtime,itime,-1)
      if (ntime==0) then
         if (kvec==1) then
            call VIBPOIS33(cue,bxyz,ffc,ci,wm,nx,ny,nz,nxeh,nye,nze,nxh,&
     &nyh,nzh)
! KNC function
         else if (kvec==2) then
            call ckncibpois33(cue,bxyz,ffc,ci,wm,nx,ny,nz,nxeh,nye,nze, &
     &nxh,nyh,nzh)
         endif
         wf = 0.0
         dth = 0.5*dt
      else
         if (kvec==1) then
            call VMAXWEL3(exyz,bxyz,cue,ffc,ci,dt,wf,wm,nx,ny,nz,nxeh,  &
     &nye,nze,nxh,nyh,nzh)
! KNC function
         else if (kvec==2) then
            call ckncmaxwel3(exyz,bxyz,cue,ffc,ci,dt,wf,wm,nx,ny,nz,nxeh&
     &,nye,nze,nxh,nyh,nzh)
         endif
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! calculate force/charge in fourier space with standard procedure:
! updates fxyze, we
      call dtimer(dtime,itime,-1)
      isign = -1
      if (kvec==1) then
         call VPOIS33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,&
     &nye,nze,nxh,nyh,nzh)
! KNC function
      else if (kvec==2) then
         call ckncpois33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,  &
     &nxeh,nye,nze,nxh,nyh,nzh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! add longitudinal and transverse electric fields with standard
! procedure: updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      if (kvec==1) then
         call VEMFIELD3(fxyze,exyz,ffc,isign,nx,ny,nz,nxeh,nye,nze,nxh, &
     &nyh,nzh)
! KNC function
      else if (kvec==2) then
         call ckncemfield3(fxyze,exyz,ffc,isign,nx,ny,nz,nxeh,nye,nze,  &
     &nxh,nyh,nzh)
      endif
! copy magnetic field with standard procedure: updates bxyze
      isign = -1
      if (kvec==1) then
         call VEMFIELD3(bxyze,bxyz,ffc,isign,nx,ny,nz,nxeh,nye,nze,nxh, &
     &nyh,nzh)
! KNC function
      else if (kvec==2) then
         call ckncemfield3(bxyze,bxyz,ffc,isign,nx,ny,nz,nxeh,nye,nze,  &
     &nxh,nyh,nzh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform electric force to real space with standard procedure:
! updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      if (kvec==1) then
         call WFFT3RV3(fxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze&
     &,nxhyz,nxyzh)
! KNC function
      else if (kvec==2) then
         call ckncwfft3rv3(fxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye&
     &,nze,nxhyz,nxyzh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! transform magnetic force to real space with standard procedure:
! updates bxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      if (kvec==1) then
         call WFFT3RV3(bxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze&
     &,nxhyz,nxyzh)
! KNC function
      else if (kvec==2) then
         call ckncwfft3rv3(bxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye&
     &,nze,nxhyz,nxyzh)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft = tfft + time
!
! copy guard cells with standard procedure: updates fxyze, bxyze
      call dtimer(dtime,itime,-1)
      if (kvec==1) then
         call CGUARD3L(fxyze,nx,ny,nz,nxe,nye,nze)
         call CGUARD3L(bxyze,nx,ny,nz,nxe,nye,nze)
! KNC function
      else if (kvec==2) then
         call cknccguard3l(fxyze,nx,ny,nz,nxe,nye,nze)
         call cknccguard3l(bxyze,nx,ny,nz,nxe,nye,nze)
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! push particles with standard procedure: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
      if (relativity==1) then
         if (kvec==1) then
            call GRBPUSH3LT(partt,fxyze,bxyze,qbme,dt,dth,ci,wke,idimp, &
     &np,npe,nx,ny,nz,nxe,nye,nze,ipbc)
!           call VGRBPUSH3LT(partt,fxyze,bxyze,qbme,dt,dth,ci,wke,idimp,&
!    &np,npe,nx,ny,nz,nxe,nye,nze,ipbc)
! KNC function
         else if (kvec==2) then
            call ckncgrbpush3lt(partt,fxyze,bxyze,qbme,dt,dth,ci,wke,   &
     &idimp,np,npe,nx,ny,nz,nxe,nye,nze,ipbc)
         endif
      else
         if (kvec==1) then
            call GBPUSH3LT(partt,fxyze,bxyze,qbme,dt,dth,wke,idimp,np,  &
     &npe,nx,ny,nz,nxe,nye,nze,ipbc)
!           call VGBPUSH3LT(partt,fxyze,bxyze,qbme,dt,dth,wke,idimp,np, &
!    &npe,nx,ny,nz,nxe,nye,nze,ipbc)
! KNC function
         else if (kvec==2) then
            call ckncgbpush3lt(partt,fxyze,bxyze,qbme,dt,dth,wke,idimp, &
     &np,npe,nx,ny,nz,nxe,nye,nze,ipbc)
         endif
      endif
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
!
! sort particles by cell for standard procedure
      if (sortime > 0) then
         if (mod(ntime,sortime)==0) then
            call dtimer(dtime,itime,-1)
            if (kvec==1) then
               call DSORTP3YZLT(partt,partt2,npic,idimp,np,npe,ny1,nyz1)
! KNC function
            else if (kvec==2) then
               call ckncdsortp3yzlt(partt,partt2,npic,idimp,np,npe,ny1, &
     &nyz1)
            endif
! exchange pointers
            tpartt => partt
            partt => partt2
            partt2 => tpartt
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
      write (*,*) 'kvec = ', kvec
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
      call avx512_deallocate(npic); nullify(npic)
      call avx512_deallocate(ffc); nullify(ffc)
      call avx512_deallocate(bxyz); nullify(bxyz)
      call avx512_deallocate(exyz); nullify(exyz)
      call avx512_deallocate(bxyze); nullify(bxyze)
      call avx512_deallocate(fxyze); nullify(fxyze)
      call avx512_deallocate(cue); nullify(cue)
      call avx512_deallocate(qe); nullify(qe)
      if (sortime > 0) then
         call avx512_deallocate(partt2); nullify(partt2)
      endif
      call avx512_deallocate(partt); nullify(partt)
      stop
      end program
!
! Procedures to create Fortran90 pointers for data allocated in C.
! For details see V. K. Decyk, ACM Fortran Forum, vol. 27, no. 2 (2008).
      subroutine getf2cptr(cref,carray,nx,ny) 
! set reference to C data in 2d real Fortran pointer object 
      implicit none 
      integer :: nx, ny
      real, dimension(nx,ny), target :: carray 
      real, dimension(:,:), pointer :: cref 
      cref => carray 
      end subroutine
!
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
! set reference to C data in vector 3d real Fortran pointer object 
      implicit none 
      integer :: ndim, nx, ny, nz
      real, dimension(ndim,nx,ny,nz), target :: carray 
      real, dimension(:,:,:,:), pointer :: cref 
      cref => carray 
      end subroutine
!
      subroutine getc2cptr(cref,carray,nx,ny) 
! set reference to C data in 2d complex Fortran pointer object 
      implicit none 
      integer :: nx, ny
      complex, dimension(nx,ny), target :: carray 
      complex, dimension(:,:), pointer :: cref 
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
!
      subroutine geti1cptr(cref,carray,nx) 
! set reference to C data in 1d integer Fortran pointer object 
      implicit none 
      integer :: nx
      integer, dimension(nx), target :: carray 
      integer, dimension(:), pointer :: cref 
      cref => carray 
      end subroutine
!
