!-----------------------------------------------------------------------
! Skeleton 2D Electrostatic MPI PIC code which uses C libraries
! written by Viktor K. Decyk, UCLA
      program ppic2
! #include "ppush2.h"
! #include "pplib2.h"
      implicit none
      integer, parameter :: indx =   9, indy =   9
      integer, parameter :: npx =  3072, npy =   3072
      integer, parameter :: ndim = 2
      real, parameter :: tend = 10.0, dt = 0.1, qme = -1.0
      real, parameter :: vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0
      real :: ax = .912871, ay = .912871
! idimp = dimension of phase space = 4
! sortime = number of time steps between standard electron sorting
      integer :: idimp = 4, ipbc = 1, sortime = 50
! idps = number of partition boundaries
      integer :: idps = 2
      real :: wke = 0.0, we = 0.0, wt = 0.0
! declare scalars for standard code
      integer :: nx, ny, nxh, nyh, nxe, nye, nxeh, nnxe, nxyh, nxhy
      integer :: ny1, ntime, nloop, isign, ierr
      real :: qbme, affp
      double precision :: np
!
! declare scalars for MPI code
      integer :: ntpose = 1
      integer :: nvp, idproc, kstrt, npmax, kxp, kyp, nypmx, nypmn
      integer :: nyp, noff, npp, nps, nbmax, ntmax
!
! declare arrays for standard code
      real, dimension(:,:), pointer :: part, part2, tpart
      real, dimension(:,:), pointer :: qe
      real, dimension(:,:,:), pointer :: fxye
      complex, dimension(:,:), pointer :: qt
      complex, dimension(:,:,:), pointer :: fxyt
      complex, dimension(:,:), pointer :: ffc
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      integer, dimension(:), pointer :: ihole
      integer, dimension(:), pointer :: npic
      double precision, dimension(4) :: wtot, work
      integer, dimension(7) :: info
!
! declare arrays for MPI code
      complex, dimension(:,:,:), pointer :: bs, br
      real, dimension(:,:), pointer :: sbufl, sbufr, rbufl, rbufr
      real, dimension(:), pointer  :: edges
      real, dimension(:), pointer  :: scr
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0
      real :: tpush = 0.0, tsort = 0.0, tmov = 0.0
      real, dimension(2) :: tfft = 0.0
      double precision :: dtime
!
! initialize scalars for standard code
      np =  dble(npx)*dble(npy)
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = max(1,ny/2)
      nxe = nx + 2; nye = ny + 2; nxeh = nxe/2; nnxe = ndim*nxe
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny); ny1 = ny + 1
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = dble(nx)*dble(ny)/np
!      
! nvp = number of distributed memory nodes
! initialize for distributed memory parallel processing
      call CPPINIT2(idproc,nvp,0,0)
      kstrt = idproc + 1
! check if too many processors
      if (nvp > ny) then
         if (kstrt==1) then
         write (*,*) 'Too many processors requested: ny, nvp=', ny, nvp
         call CPPEXIT()
         stop
         endif
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
      call CPDICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,idps)
      if (nypmn < 1) then
         if (kstrt==1) then
            write (*,*) 'combination not supported nvp, ny =',nvp,ny
         endif
         go to 3000
      endif
! initialize additional scalars for MPI code
! kxp = number of complex grids in each field partition in x direction
      kxp = (nxh - 1)/nvp + 1
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvp + 1
! npmax = maximum number of electrons in each partition
      npmax = (np/nvp)*1.25
! nbmax = size of buffer for passing particles between processors
      nbmax = 0.1*npmax
! ntmax = size of ihole buffer for particles leaving processor
      ntmax = 2*nbmax
!
! allocate and initialize data for standard code
      allocate(part(idimp,npmax),part2(idimp,npmax))
      allocate(qe(nxe,nypmx),fxye(ndim,nxe,nypmx))
      allocate(qt(nye,kxp),fxyt(ndim,nye,kxp))
      allocate(ffc(nyh,kxp),mixup(nxhy),sct(nxyh))
      allocate(ihole(ntmax+1),npic(nypmx))
!
! allocate and initialize data for MPI code
      allocate(bs(ndim,kxp,kyp),br(ndim,kxp,kyp))
      allocate(sbufl(idimp,nbmax),sbufr(idimp,nbmax))
      allocate(rbufl(idimp,nbmax),rbufr(idimp,nbmax))
      allocate(scr(nxe))
!
! prepare fft tables
      call CWPFFT2RINIT(mixup,sct,indx,indy,nxhy,nxyh)
! calculate form factors
      isign = 0
      call CPPOIS22(qt,fxyt,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nye,kxp,&
     &nyh)
! initialize electrons
      nps = 1
      npp = 0
      call CPDISTR2(part,edges,npp,nps,vtx,vty,vx0,vy0,npx,npy,nx,ny,   &
     &idimp,npmax,idps,ipbc,ierr)
! check for particle initialization error
      if (ierr /= 0) then
         if (kstrt==1) then
            write (*,*) 'particle initialization error: ierr=', ierr
         endif
         go to 3000
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
!     if (kstrt==1) write (*,*) 'ntime = ', ntime
!
! deposit charge with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call CPPGPOST2L(part,qe,npp,noff,qme,idimp,npmax,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tdpost = tdpost + time
!
! add guard cells with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      call CPPAGUARD2XL(qe,nyp,nx,nxe,nypmx)
      call CPPNAGUARD2L(qe,scr,nyp,nx,kstrt,nvp,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! transform charge to fourier space with standard procedure: updates qt
! modifies qe
      call dtimer(dtime,itime,-1)
      isign = -1
      call CWPPFFT2R(qe,qt,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,  &
     &kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! calculate force/charge in fourier space with standard procedure:
! updates fxyt
      call dtimer(dtime,itime,-1)
      isign = -1
      call CPPOIS22(qt,fxyt,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nye,kxp,&
     &nyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfield = tfield + time
!
! transform force to real space with standard procedure: updates fxye
! modifies fxyt
      call dtimer(dtime,itime,-1)
      isign = 1
      call CWPPFFT2R2(fxye,fxyt,bs,br,isign,ntpose,mixup,sct,ttp,indx,  &
     &indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tfft(1) = tfft(1) + time
      tfft(2) = tfft(2) + ttp
!
! copy guard cells with standard procedure: updates fxye
      call dtimer(dtime,itime,-1)
      call CPPNCGUARD2L(fxye,nyp,kstrt,nvp,nnxe,nypmx)
      call CPPCGUARD2XL(fxye,nyp,nx,ndim,nxe,nypmx)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tguard = tguard + time
!
! push particles: updates part, wke, and ihole
      call dtimer(dtime,itime,-1)
      wke = 0.0
      call CPPGPUSH2L(part,fxye,edges,npp,noff,ihole,qbme,dt,wke,nx,ny, &
     &idimp,npmax,nxe,nypmx,idps,ntmax,ipbc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tpush = tpush + time
! check for ihole overflow error
      if (ihole(1) < 0) then
         ierr = -ihole(1)
         write (*,*) kstrt,'ihole overflow error: ntmax,ih=', ntmax,ierr
         call CPPABORT()
         go to 3000
      endif
!
! move electrons into appropriate spatial regions: updates part, npp
      call dtimer(dtime,itime,-1)
      call CPPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,ny,    &
     &kstrt,nvp,idimp,npmax,idps,nbmax,ntmax,info)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      tmov = tmov + time
! check for particle manager error
      if (info(1) /= 0) then
         ierr = info(1)
         if (kstrt==1) then
            write (*,*) 'particle manager error: ierr=', ierr
         endif
         go to 3000
      endif
!
! sort particles for standard code: updates part
      if (sortime > 0) then
         if (mod(ntime,sortime)==0) then
            call dtimer(dtime,itime,-1)
            call CPPDSORTP2YL(part,part2,npic,npp,noff,nyp,idimp,npmax, &
     &nypmx)
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
! energy diagnostic
      wtot(1) = we
      wtot(2) = wke
      wtot(3) = 0.0
      wtot(4) = we + wke
      call CPPDSUM(wtot,work,4)
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
      call CPPEXIT()
      stop
      end
