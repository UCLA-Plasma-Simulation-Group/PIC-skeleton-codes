c Fortran Library for Skeleton 3D Electrostatic MPI/OpenMP PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine PDICOMP32L(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,ny, 
     1nz,kstrt,nvpy,nvpz,idps,idds)
c this subroutine determines spatial boundaries for uniform particle
c decomposition, calculates number of grid points in each spatial
c region, and the offset of these grid points from the global address
c nvpy must be < ny and nvpz must be < nz.
c some combinations of ny and nvpy and nz and nvpz result in a zero
c value of nyzp.  this is not supported.
c input: ny, nz, kstrt, nvpy, nvpz, idps, idds
c output: edges, nyzp, noff, nypmx, nzpmx, nypmn, nzpmn
c for 2D spatial decomposition
c edges(1:2) = lower/upper boundary in y of particle partition
c edges(3:4) = back/front boundary in z of particle partition
c nyzp(1:2) = number of primary (complete) gridpoints in y/z
c noff(1) = lowermost global gridpoint in y in particle partition
c noff(2) = backmost global gridpoint in z in particle partition
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c nypmn = minimum value of nyzp(1)
c nzpmn = minimum value of nyzp(2)
c ny/nz = system length in y/z direction
c kstrt = starting data block number (processor id + 1)
c nvpy/nvpz = number of real or virtual processors in y/z
c idps = number of particle partition boundaries = 4
c idds = dimensionality of domain decomposition = 2
      implicit none
      integer nypmx, nzpmx, nypmn, nzpmn, ny, nz, kstrt, nvpy, nvpz
      integer idps, idds
      integer nyzp, noff
      real edges
      dimension nyzp(idds), noff(idds)
      dimension edges(idps)
c local data
      integer jb, kb, kyp, kzp
      real at1, at2, any, anz
      integer myzpm, iwork4
      dimension myzpm(4), iwork4(4)
      any = real(ny)
      anz = real(nz)
c determine decomposition
c find processor id in y/z
      kb = (kstrt - 1)/nvpy
      jb = kstrt - nvpy*kb - 1
c boundaries in y
      kyp = (ny - 1)/nvpy + 1
      at1 = real(kyp)
      edges(1) = at1*real(jb)
      if (edges(1).gt.any) edges(1) = any
      noff(1) = edges(1)
      edges(2) = at1*real(jb + 1)
      if (edges(2).gt.any) edges(2) = any
      jb = edges(2)
      nyzp(1) = jb - noff(1)
c boundaries in z
      kzp = (nz - 1)/nvpz + 1
      at2 = real(kzp)
      edges(3) = at2*real(kb)
      if (edges(3).gt.anz) edges(3) = anz
      noff(2) = edges(3)
      edges(4) = at2*real(kb + 1)
      if (edges(4).gt.anz) edges(4) = anz
      kb = edges(4)
      nyzp(2) = kb - noff(2)
c find maximum/minimum partition size in y and z
      myzpm(1) = nyzp(1)
      myzpm(2) = -nyzp(1)
      myzpm(3) = nyzp(2)
      myzpm(4) = -nyzp(2)
      call PPIMAX(myzpm,iwork4,4)
      nypmx = myzpm(1) + 1
      nypmn = -myzpm(2)
      nzpmx = myzpm(3) + 1
      nzpmn = -myzpm(4)
      return
      end
c-----------------------------------------------------------------------
      subroutine FCOMP32(nvp,nx,ny,nz,nvpy,nvpz,ierr)
c determines optimal partition for nvp processors
c input: nvp, number of processors, nx, ny, nz = number of grids
c output: nvpy, nvpz, processors in y, z direction, ierr = error code
c nvp = number of real or virtual processors obtained
c nx/ny/nz = system length in x/y/z direction
c nvpy/nvpz = number of real or virtual processors in y/z
c ierr = (0,1) = (no,yes) error condition exists
      implicit none
      integer nvp, nx, ny, nz, nvpy, nvpz, ierr
c local data
      integer nxh, lvp
      double precision dt1
      nxh = nx/2
      ierr = 0
c algorithm 1: prefer equal number of grids in y and z partitions
      dt1 = sqrt(dble(nvp)*dble(ny)/dble(nz))
c algorithm 2: prefer equal number of grids in x*y and y*z partitions
c     dt1 = sqrt(nvp*sqrt(dble(nxh)/dble(nz)))
c return total number of processors in y and z
      nvpy = real(dt1)
      if (nvpy.lt.1) nvpy = 1
      nvpz = nvp/nvpy
      lvp = nvpy*nvpz
      if (lvp.gt.nvp) then
         write (*,*) 'invalid partition:nvpy,nvpz,nvp=', nvpy, nvpz, nvp
         ierr = 1
         return
      endif
   10 if (lvp.ne.nvp) then
         nvpy = nvpy - 1
         nvpz = nvp/nvpy
         lvp = nvpy*nvpz
         go to 10
      endif
      nvp = lvp
      return
      end
c-----------------------------------------------------------------------
      subroutine PDISTR32(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx
     1,npy,npz,nx,ny,nz,idimp,npmax,idps,ipbc,ierr)
c for 3d code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c for distributed data with 2D spatial decomposition
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = position z of particle n in partition
c part(4,n) = velocity vx of particle n in partition
c part(5,n) = velocity vy of particle n in partition
c part(6,n) = velocity vz of particle n in partition
c edges(1) = lower boundary in y of particle partition
c edges(2) = upper boundary in y of particle partition
c edges(3) = back boundary in z of particle partition
c edges(4) = front boundary in z of particle partition
c npp = number of particles in partition
c nps = starting address of particles in partition
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c idps = number of partition boundaries
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
c ierr = (0,1) = (no,yes) error condition exists
c ranorm = gaussian random number with zero mean and unit variance
c with 2D spatial decomposition
      implicit none
      integer npp, nps, npx, npy, npz, nx, ny, nz, idimp, npmax, idps
      integer ipbc, ierr
      real vtx, vty, vtz, vdx, vdy, vdz
      real part, edges
      dimension part(idimp,npmax), edges(idps)
c local data
      integer j, k, l, npt, npxyzp
      real edgelx, edgely, edgelz, at1, at2, at3
      real xt, yt, zt, vxt, vyt, vzt
      double precision dnpxy, dnpxyz, dt1
      integer ierr1, iwork1
      double precision sum4, work4
      dimension ierr1(1), iwork1(1), sum4(4), work4(4)
      double precision ranorm
      ierr = 0
c particle distribution constant
      dnpxy = dble(npx)*dble(npy)
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      at3 = real(nz)/real(npz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
         at3 = real(nz-2)/real(npz)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 0.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
         at3 = real(nz)/real(npz)
      endif
      do 30 l = 1, npz
      zt = edgelz + at3*(real(l) - 0.5)
      do 20 k = 1, npy
      yt = edgely + at2*(real(k) - 0.5)
      do 10 j = 1, npx
c uniform density profile
      xt = edgelx + at1*(real(j) - 0.5)
c maxwellian velocity distribution
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      vzt = vtz*ranorm()
      if ((yt.ge.edges(1)).and.(yt.lt.edges(2))) then
         if ((zt.ge.edges(3)).and.(zt.lt.edges(4))) then
            npt = npp + 1
            if (npt.le.npmax) then
               part(1,npt) = xt
               part(2,npt) = yt
               part(3,npt) = zt
               part(4,npt) = vxt
               part(5,npt) = vyt
               part(6,npt) = vzt
               npp = npt
            else
               ierr = ierr + 1
            endif
         endif
      endif
   10 continue
   20 continue
   30 continue
      npxyzp = 0
c add correct drift
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 40 j = nps, npp
      npxyzp = npxyzp + 1
      sum4(1) = sum4(1) + part(4,j)
      sum4(2) = sum4(2) + part(5,j)
      sum4(3) = sum4(3) + part(6,j)
   40 continue
      sum4(4) = npxyzp
      call PPDSUM(sum4,work4,4)
      dnpxyz = sum4(4)
      ierr1(1) = ierr
      call PPIMAX(ierr,iwork1,1)
      ierr = ierr1(1)
      dt1 = 1.0d0/dnpxyz
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 50 j = nps, npp
      part(4,j) = part(4,j) - sum4(1)
      part(5,j) = part(5,j) - sum4(2)
      part(6,j) = part(6,j) - sum4(3)
   50 continue
c process errors
      dnpxyz = dnpxyz - dnpxy*dble(npz)
      if (dnpxyz.ne.0.0d0) ierr = dnpxyz
      return
      end
c-----------------------------------------------------------------------
      subroutine PPDBLKP3L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,mz
     1,mx1,myp1,mxyzp1,idds,irc)
c this subroutine finds the maximum number of particles in each tile of
c mx, my, mz to calculate size of segmented particle array ppart
c linear interpolation, spatial decomposition in y/z direction
c input: all except kpic, nppmx, output: kpic, nppmx
c part = input particle array
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = position z of particle n in partition
c kpic = output number of particles per tile
c npp = number of particles in partition
c noff(1) = lowermost global gridpoint in y in particle partition
c noff(2) = backmost global gridpoint in z in particle partition
c nppmx = return maximum number of particles in tile
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mx/my/mz = number of grids in sorting cell in x, y and z
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c mxyzp1 = mx1*myp1*mzp1,
c where mzp1 = (partition length in z direction - 1)/mz + 1
c idds = dimensionality of domain decomposition
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer npp, nppmx, idimp, npmax, mx, my, mz, mx1, myp1, mxyzp1
      integer idds, irc
      integer kpic, noff
      real part
      dimension part(idimp,npmax), kpic(mxyzp1)
      dimension noff(idds)
c local datal, 
      integer j, k, n, m, l, mnoff, lnoff, mxyp1, isum, ist, npx, ierr
      mnoff = noff(1)
      lnoff = noff(2)
      ierr = 0
      mxyp1 = mx1*myp1
c clear counter array
      do 10 k = 1, mxyzp1
      kpic(k) = 0
   10 continue
c find how many particles in each tile
      do 20 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      l = part(3,j)
      l = (l - lnoff)/mz
      m = n + mx1*m + mxyp1*l
      if (m.le.mxyzp1) then
         kpic(m) = kpic(m) + 1
      else
         ierr = max(ierr,m-mxyzp1)
      endif
   20 continue
c find maximum
      isum = 0
      npx = 0
      do 30 k = 1, mxyzp1
      ist = kpic(k)
      npx = max(npx,ist)
      isum = isum + ist
   30 continue
      nppmx = npx
c check for errors
      if (ierr.gt.0) then
         irc = ierr
      else if (isum.ne.npp) then
         irc = -1
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPMOVIN3L(part,ppart,kpic,npp,noff,nppmx,idimp,npmax, 
     1mx,my,mz,mx1,myp1,mxyzp1,idds,irc)
c this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
c and copies to segmented array ppart
c linear interpolation, spatial decomposition in y/z direction
c input: all except ppart, kpic, output: ppart, kpic
c part/ppart = input/output particle arrays
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = position z of particle n in partition
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c ppart(3,n,m) = position z of particle n in tile m
c ppart(4,n,m) = velocity vx of particle n in tile m
c ppart(5,n,m) = velocity vy of particle n in tile m
c ppart(6,n,m) = velocity vz of particle n in tile m
c kpic = output number of particles per tile
c npp = number of particles in partition
c noff(1) = lowermost global gridpoint in y in particle partition
c noff(2) = backmost global gridpoint in z in particle partition
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mx/my/mz = number of grids in sorting cell in x, y and z
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c mxyzp1 = mx1*myp1*mzp1,
c where mzp1 = (partition length in z direction - 1)/mz + 1
c idds = dimensionality of domain decomposition
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer npp, nppmx, idimp, npmax, mx, my, mz, mx1, myp1, mxyzp1
      integer idds, irc
      integer kpic, noff
      real part, ppart
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyzp1)
      dimension kpic(mxyzp1), noff(idds)
c local data
      integer i, j, k, n, m, l, mnoff, lnoff, mxyp1, ip, ierr
      mnoff = noff(1)
      lnoff = noff(2)
      ierr = 0
      mxyp1 = mx1*myp1
c clear counter array
      do 10 k = 1, mxyzp1
      kpic(k) = 0
   10 continue
c find addresses of particles at each tile and reorder particles
      do 30 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      l = part(3,j)
      l = (l - lnoff)/mz
      m = n + mx1*m + mxyp1*l
      ip = kpic(m) + 1
      if (ip.le.nppmx) then
         do 20 i = 1, idimp
         ppart(i,ip,m) = part(i,j)
   20    continue
      else
         ierr = max(ierr,ip-nppmx)
      endif
      kpic(m) = ip
   30 continue
      if (ierr.gt.0) irc = ierr
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPCHECK3L(ppart,kpic,noff,nyzp,idimp,nppmx,nx,mx,my,mz
     1,mx1,myp1,mzp1,idds,irc)
c this subroutine performs a sanity check to make sure particles sorted
c by x,y,z grid in tiles of mx, my, mz, are all within bounds.
c tiles are assumed to be arranged in 3D linear memory
c input: all except irc
c output: irc
c ppart(1,n,l) = position x of particle n in tile l
c ppart(2,n,l) = position y of particle n in tile l
c ppart(3,n,l) = position a of particle n in tile l
c kpic(l) = number of reordered output particles in tile l
c noff(1) = lowermost global gridpoint in y in particle partition
c noff(2) = backmost global gridpoint in z in particle partition
c nyzp(1:2) = number of primary (complete) gridpoints in y/z
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx = system length in x
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c mzp1 = (partition length in z direction - 1)/mz + 1
c idds = dimensionality of domain decomposition
c irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, mx, my, mz, mx1, myp1, mzp1, idds, irc
      real ppart
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mx1*myp1*mzp1)
      dimension kpic(mx1*myp1*mzp1), noff(idds), nyzp(idds)
c local data
      integer mxyp1, mxyzp1, noffp, moffp, loffp, nppp
      integer j, k, l, nn, mm, ll, ist
      real edgelx, edgely, edgelz, edgerx, edgery, edgerz, dx, dy, dz
      mxyp1 = mx1*myp1
      mxyzp1 = mxyp1*mzp1
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,l,noffp,moffp,loffp,nppp,nn,mm,ll,ist,edgelx,edgely, 
!$OMP& edgelz,edgerx,edgery,edgerz,dx,dy,dz)
      do 20 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      nn = min(mx,nx-noffp)
      mm = min(my,nyzp(1)-moffp)
      ll = min(mz,nyzp(2)-loffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff(1) + moffp
      edgery = noff(1) + moffp + mm
      edgelz = noff(2) + loffp
      edgerz = noff(2) + loffp + ll
c loop over particles in tile
      do 10 j = 1, nppp
      dx = ppart(1,j,l)
      dy = ppart(2,j,l)
      dz = ppart(3,j,l)
c find particles going out of bounds
      ist = 0
      if (dx.lt.edgelx) ist = 1
      if (dx.ge.edgerx) ist = 2
      if (dy.lt.edgely) ist = ist + 3
      if (dy.ge.edgery) ist = ist + 6
      if (dz.lt.edgelz) ist = ist + 9
      if (dz.ge.edgerz) ist = ist + 18
      if (ist.gt.0) irc = l
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGPPUSH32L(ppart,fxyz,kpic,noff,nyzp,qbm,dt,ek,idimp, 
     1nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions.
c for distributed data, with 2D spatial decomposition
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 90 flops/particle, 30 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c ppart(3,n,m) = position z of particle n in partition in tile m
c ppart(4,n,m) = velocity vx of particle n in partition in tile m
c ppart(5,n,m) = velocity vy of particle n in partition in tile m
c ppart(6,n,m) = velocity vz of particle n in partition in tile m
c fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,ll)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
c that is, convolution of electric field over particle shape
c where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
c kpic = number of particles per tile
c noff(1) = lowermost global gridpoint in y in particle partition
c noff(2) = backmost global gridpoint in z in particle partition
c nyzp(1:2) = number of primary (complete) gridpoints in y/z
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = system length in x/y/z direction
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = second dimension of field array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c mxyzp1 = mx1*myp1*mzp1,
c where mzp1 = (partition length in z direction - 1)/mz + 1
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds, ipbc
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds), nyzp(idds)
c local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll
      real qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real sfxyz
      dimension sfxyz(3,MXV,MYV,MZV)
c     dimension sfxyz(3,mx+1,my+1,mz+1)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
      qtm = qbm*dt
      sum2 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      edgerz = real(nz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
         edgerz = real(nz-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      endif
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,nn,mm,ll,x,y,z
!$OMP& ,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,sum1,sfxyz)
!$OMP& REDUCTION(+:sum2)
      do 50 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
c load local fields from global array
      do 30 k = 1, min(mz,nyzp(2)-loffp)+1
      do 20 j = 1, min(my,nyzp(1)-moffp)+1
      do 10 i = 1, min(mx,nx-noffp)+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 40 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find acceleration
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)  
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)  
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)  
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            
     1             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            
     1             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            
     1             + dx1*sfxyz(3,nn+1,mm+1,ll))
      vx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      vy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      vz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,nn,mm+1,ll+1)                     
     1                  + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(vy + dyp*sfxyz(2,nn,mm+1,ll+1)                     
     1                  + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(vz + dyp*sfxyz(3,nn,mm+1,ll+1)                     
     1                  + dx1*sfxyz(3,nn+1,mm+1,ll+1))
c new velocity
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      dx = vx + qtm*dx
      dy = vy + qtm*dy
      dz = vz + qtm*dz
c average kinetic energy
      vx = vx + dx
      vy = vy + dy
      vz = vz + dz
      sum1 = sum1 + (vx*vx + vy*vy+ vz*vz)
      ppart(4,j,l) = dx
      ppart(5,j,l) = dy
      ppart(6,j,l) = dz
c new position
      dx = x + dx*dt
      dy = y + dy*dt
      dz = z + dz*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -ppart(4,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -ppart(5,j,l)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            ppart(6,j,l) = -ppart(6,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(4,j,l) = -ppart(4,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(5,j,l) = -ppart(5,j,l)
         endif
      endif
c set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
   40 continue
      sum2 = sum2 + sum1
   50 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGPPUSHF32L(ppart,fxyz,kpic,ncl,ihole,noff,nyzp,qbm,dt
     1,ek,idimp,nppmx,nx,ny,nz,mx,my,mz,nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,
     2ntmax,idds,irc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c for distributed data, with 2D spatial decomposition
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 90 flops/particle, 30 loads, 6 stores
c input: all except ncl, ihole, irc, output: part, ncl, ihole, ek, irc
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c ppart(3,n,m) = position z of particle n in partition in tile m
c ppart(4,n,m) = velocity vx of particle n in partition in tile m
c ppart(5,n,m) = velocity vy of particle n in partition in tile m
c ppart(6,n,m) = velocity vz of particle n in partition in tile m
c fxyz(1,j,k,l) = x component of force/charge at grid (j,kk,ll)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,kk,ll)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,kk,ll)
c that is, convolution of electric field over particle shape
c where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
c kpic(l) = number of particles in tile l
c ncl(i,l) = number of particles going to destination i, tile l
c ihole(1,:,l) = location of hole in array left by departing particle
c ihole(2,:,l) = direction destination of particle leaving hole
c all for tile l
c ihole(1,1,l) = ih, number of holes left (error, if negative)
c noff(1) = lowermost global gridpoint in y in particle partition
c noff(2) = backmost global gridpoint in z in particle partition
c nyzp(1:2) = number of primary (complete) gridpoints in y/z
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = system length in x/y/z direction
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = second dimension of field array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c mxyzp1 = mx1*myp1*mzp1,
c where mzp1 = (partition length in z direction - 1)/mz + 1
c ntmax = size of hole array for particles leaving tiles
c idds = dimensionality of domain decomposition
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, ntmax, idds, irc
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic, ncl, ihole, noff, nyzp
      dimension ppart(idimp,nppmx,mxyzp1), fxyz(3,nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), ncl(26,mxyzp1), ihole(2,ntmax+1,mxyzp1)
      dimension noff(idds), nyzp(idds)
c local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, ih, nh, nn, mm, ll
      real qtm, x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real sfxyz
      dimension sfxyz(3,MXV,MYV,MZV)
c     dimension sfxyz(3,mx+1,my+1,mz+1)
      double precision sum1, sum2
      mxyp1 = mx1*myp1
      qtm = qbm*dt
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
      sum2 = 0.0d0
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,ih,nh,nn,mm,ll
!$OMP& ,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,edgelx,     
!$OMP& edgely,edgelz,edgerx,edgery,edgerz,sum1,sfxyz)
!$OMP& REDUCTION(+:sum2)
      do 60 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      nn = min(mx,nx-noffp)
      mm = min(my,nyzp(1)-moffp)
      ll = min(mz,nyzp(2)-loffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff(1) + moffp
      edgery = noff(1) + moffp + mm
      edgelz = noff(2) + loffp
      edgerz = noff(2) + loffp + ll
      ih = 0
      nh = 0
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
c load local fields from global array
      do 30 k = 1, ll+1
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxyz(1,i,j,k) = fxyz(1,i+noffp,j+moffp,k+loffp)
      sfxyz(2,i,j,k) = fxyz(2,i+noffp,j+moffp,k+loffp)
      sfxyz(3,i,j,k) = fxyz(3,i+noffp,j+moffp,k+loffp)
   10 continue
   20 continue
   30 continue
c clear counters
      do 40 j = 1, 26
      ncl(j,l) = 0
   40 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 50 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find acceleration
      dx = amx*sfxyz(1,nn,mm,ll) + amy*sfxyz(1,nn+1,mm,ll)  
      dy = amx*sfxyz(2,nn,mm,ll) + amy*sfxyz(2,nn+1,mm,ll)  
      dz = amx*sfxyz(3,nn,mm,ll) + amy*sfxyz(3,nn+1,mm,ll)  
      dx = amz*(dx + dyp*sfxyz(1,nn,mm+1,ll)                            
     1             + dx1*sfxyz(1,nn+1,mm+1,ll))
      dy = amz*(dy + dyp*sfxyz(2,nn,mm+1,ll)                            
     1             + dx1*sfxyz(2,nn+1,mm+1,ll))
      dz = amz*(dz + dyp*sfxyz(3,nn,mm+1,ll)                            
     1             + dx1*sfxyz(3,nn+1,mm+1,ll))
      vx = amx*sfxyz(1,nn,mm,ll+1) + amy*sfxyz(1,nn+1,mm,ll+1)
      vy = amx*sfxyz(2,nn,mm,ll+1) + amy*sfxyz(2,nn+1,mm,ll+1)
      vz = amx*sfxyz(3,nn,mm,ll+1) + amy*sfxyz(3,nn+1,mm,ll+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,nn,mm+1,ll+1)                     
     1                  + dx1*sfxyz(1,nn+1,mm+1,ll+1))
      dy = dy + dzp*(vy + dyp*sfxyz(2,nn,mm+1,ll+1)                     
     1                  + dx1*sfxyz(2,nn+1,mm+1,ll+1))
      dz = dz + dzp*(vz + dyp*sfxyz(3,nn,mm+1,ll+1)                     
     1                  + dx1*sfxyz(3,nn+1,mm+1,ll+1))
c new velocity
      vx = ppart(4,j,l)
      vy = ppart(5,j,l)
      vz = ppart(6,j,l)
      dx = vx + qtm*dx
      dy = vy + qtm*dy
      dz = vz + qtm*dz
c average kinetic energy
      vx = vx + dx
      vy = vy + dy
      vz = vz + dz
      sum1 = sum1 + (vx*vx + vy*vy+ vz*vz)
      ppart(4,j,l) = dx
      ppart(5,j,l) = dy
      ppart(6,j,l) = dz
c new position
      dx = x + dx*dt
      dy = y + dy*dt
      dz = z + dz*dt
c find particles going out of bounds
      mm = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) dz = dz - anz
         mm = mm + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               mm = mm + 9
            else
               dz = 0.0
            endif
         else
            mm = mm + 9
         endif
      endif
c set new position
      ppart(1,j,l) = dx
      ppart(2,j,l) = dy
      ppart(3,j,l) = dz
c increment counters
      if (mm.gt.0) then
         ncl(mm,l) = ncl(mm,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j
            ihole(2,ih+1,l) = mm
         else
            nh = 1
         endif
      endif
   50 continue
      sum2 = sum2 + sum1
c set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
   60 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGPPOST32L(ppart,q,kpic,noff,qm,nppmx,idimp,mx,my,mz, 
     1nxv,nypmx,nzpmx,mx1,myp1,mxyzp1,idds)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c for distributed data, with 2D spatial decomposition
c OpenMP version using guard cells
c data deposited in tiles
c particles stored segmented array
c 33 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c ppart(3,n,m) = position z of particle n in partition in tile m
c q(j,k,l) = charge density at grid point (j,kk,ll),
c where kk = k + noff(1) - 1, and ll = l + noff(2) - 1
c kpic = number of particles per tile
c noff(1) = lowermost global gridpoint in y in particle partition
c noff(2) = backmost global gridpoint in z in particle partition
c qm = charge on particle, in units of e
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 6
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = first dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c mxyzp1 = mx1*myp1*mzp1,
c where mzp1 = (partition length in z direction - 1)/mz + 1
c idds = dimensionality of domain decomposition
      implicit none
      integer nppmx, idimp, mx, my, mz, nxv, nypmx, nzpmx
      integer mx1, myp1, mxyzp1, idds
      real qm
      real ppart, q
      integer kpic, noff
      dimension ppart(idimp,nppmx,mxyzp1), q(nxv,nypmx,nzpmx)
      dimension kpic(mxyzp1), noff(idds)
c local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxyp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, mnoff, lnoff, nn, mm, ll, nm, lm
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1
      real sq
c     dimension sq(MXV,MYV,MZV)
      dimension sq(mx+1,my+1,mz+1)
      mxyp1 = mx1*myp1
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noffp,moffp,loffp,nppp,mnoff,lnoff,nn,mm,ll,nm,lm
!$OMP& ,x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,sq)
      do 150 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      mnoff = moffp + noff(1) - 1
      lnoff = loffp + noff(2) - 1
c zero out local accumulator
      do 30 k = 1, mz+1
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      sq(i,j,k) = 0.0
   10 continue
   20 continue
   30 continue
c loop over particles in tile
      do 40 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,l)
      y = ppart(2,j,l)
      z = ppart(3,j,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noffp + 1
      mm = mm - mnoff
      ll = ll - lnoff
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c deposit charge within tile to local accumulator
      x = sq(nn,mm,ll) + amx*amz
      y = sq(nn+1,mm,ll) + amy*amz
      sq(nn,mm,ll) = x
      sq(nn+1,mm,ll) = y
      x = sq(nn,mm+1,ll) + dyp*amz
      y = sq(nn+1,mm+1,ll) + dx1*amz
      sq(nn,mm+1,ll) = x
      sq(nn+1,mm+1,ll) = y
      x = sq(nn,mm,ll+1) + amx*dzp
      y = sq(nn+1,mm,ll+1) + amy*dzp
      sq(nn,mm,ll+1) = x
      sq(nn+1,mm,ll+1) = y
      x = sq(nn,mm+1,ll+1) + dyp*dzp
      y = sq(nn+1,mm+1,ll+1) + dx1*dzp
      sq(nn,mm+1,ll+1) = x
      sq(nn+1,mm+1,ll+1) = y
   40 continue
c deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      ll = min(mz,nzpmx-loffp)
      do 70 k = 2, ll
      do 60 j = 2, mm
      do 50 i = 2, nn
      q(i+noffp,j+moffp,k+loffp) = q(i+noffp,j+moffp,k+loffp)
     1 + sq(i,j,k)
   50 continue
   60 continue
   70 continue
c deposit charge to edge points in global array
      lm = min(mz+1,nzpmx-loffp)
      do 90 j = 2, mm
      do 80 i = 2, nn
!$OMP ATOMIC
      q(i+noffp,j+moffp,1+loffp) = q(i+noffp,j+moffp,1+loffp)
     1 + sq(i,j,1)
      if (lm > mz) then
!$OMP ATOMIC
         q(i+noffp,j+moffp,lm+loffp) = q(i+noffp,j+moffp,lm+loffp)
     1   + sq(i,j,lm)
      endif
   80 continue
   90 continue
      nm = min(mx+1,nxv-noffp)
      mm = min(my+1,nypmx-moffp)
      do 120 k = 1, ll
      do 100 i = 2, nn
!$OMP ATOMIC
      q(i+noffp,1+moffp,k+loffp) = q(i+noffp,1+moffp,k+loffp)
     1 + sq(i,1,k)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noffp,mm+moffp,k+loffp) = q(i+noffp,mm+moffp,k+loffp)
     1   + sq(i,mm,k)
      endif
  100 continue
      do 110 j = 1, mm
!$OMP ATOMIC
      q(1+noffp,j+moffp,k+loffp) = q(1+noffp,j+moffp,k+loffp)
     1 + sq(1,j,k)
      if (nm > mx) then
!$OMP ATOMIC
         q(nm+noffp,j+moffp,k+loffp) = q(nm+noffp,j+moffp,k+loffp) 
     1   + sq(nm,j,k)
      endif
  110 continue
  120 continue
      if (lm > mz) then
         do 130 i = 2, nn
!$OMP ATOMIC
         q(i+noffp,1+moffp,lm+loffp) = q(i+noffp,1+moffp,lm+loffp)
     1   + sq(i,1,lm)
         if (mm > my) then
!$OMP ATOMIC
            q(i+noffp,mm+moffp,lm+loffp) = q(i+noffp,mm+moffp,lm+loffp)
     1      + sq(i,mm,lm)
         endif
  130    continue
         do 140 j = 1, mm
!$OMP ATOMIC
         q(1+noffp,j+moffp,lm+loffp) = q(1+noffp,j+moffp,lm+loffp)
     1   + sq(1,j,lm)
         if (nm > mx) then
!$OMP ATOMIC
            q(nm+noffp,j+moffp,lm+loffp) = q(nm+noffp,j+moffp,lm+loffp) 
     1      + sq(nm,j,lm)
         endif
  140    continue
      endif
  150 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPORDER32LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,  
     1ncll,nclr,noff,nyzp,idimp,nppmx,nx,ny,nz,mx,my,mz,mx1,myp1,mzp1,
     2mxzyp1,npbmx,ntmax,nbmax,idds,irc)
c this subroutine performs first part of a particle sort by x,y,z grid
c in tiles of mx, my, mz
c linear interpolation, with periodic boundary conditions
c for distributed data, with 2d domain decomposition in y/z.
c tiles are assumed to be arranged in 3D linear memory
c this part of the algorithm has 3 steps.  first, one finds particles
c leaving tile and stores their number in each directon, location, and
c destination in ncl and ihole.  then, a prefix scan of ncl is performed
c and departing particles are buffered in ppbuff in direction order.
c finally, we buffer particles leaving the processor in y/z direction in
c sbufl and sbufr, and store particle number offsets in ncll and nclr.
c input: all except ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
c output: ppart, ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c ppart(3,n,m) = position z of particle n in tile m
c ppbuff(i,n,l) = i co-ordinate of particle n in tile l
c sbufl = buffer for particles being sent to lower/back processor
c sbufr = buffer for particles being sent to upper/forward processor
c kpic(l) = number of particles in tile l
c ncl(i,l) = number of particles going to destination i, tile l
c ihole(1,:,l) = location of hole in array left by departing particle
c ihole(2,:,l) = direction destination of particle leaving hole
c all for tile l
c ihole(1,1,l) = ih, number of holes left (error, if negative)
c ncll = number offset being sent to lower/back processor
c nclr = number offset being sent to upper/forward processor
c noff(1) = lowermost global gridpoint in y in particle partition
c noff(2) = backmost global gridpoint in z in particle partition
c nyzp(1:2) = number of primary (complete) gridpoints in y/z
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = system length in x/y/z direction
c mx/my/mz = number of grids in sorting cell in x/y/z
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c mzp1 = (partition length in z direction - 1)/mz + 1
c mxzyp1 = mx1*max(myp1,mzp1)
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c nbmax =  size of buffers for passing particles between processors
c idds = dimensionality of domain decomposition
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, myp1, mzp1
      integer mxzyp1, npbmx, ntmax, nbmax, idds, irc
      real ppart, ppbuff, sbufl, sbufr
      integer kpic, noff, nyzp, ncl, ihole, ncll, nclr
      dimension ppart(idimp,nppmx,mx1*myp1*mzp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1*mzp1)
      dimension sbufl(idimp,nbmax,2), sbufr(idimp,nbmax,2)
      dimension kpic(mx1*myp1*mzp1), noff(idds), nyzp(idds)
      dimension ncl(26,mx1*myp1*mzp1), ihole(2,ntmax+1,mx1*myp1*mzp1)
      dimension ncll(3,mxzyp1,3,2), nclr(3,mxzyp1,3,2)
c local data
      integer mxyp1, mxzp1, mxyzp1, noffp, moffp, loffp, nppp
      integer i, j, k, l, n, ii, ih, nh, ist, nn, mm, ll, isum
      integer ip, j1, j2, k1, k2, kk, nr, nl
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dx, dy, dz
      mxyp1 = mx1*myp1
      mxyzp1 = mxyp1*mzp1
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
c find and count particles leaving tiles and determine destination
c update ppart, ihole, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,l,noffp,moffp,loffp,nppp,nn,mm,ll,ih,nh,ist,dx,dy,dz,
!$OMP& edgelx,edgely,edgelz,edgerx,edgery,edgerz)
      do 30 l = 1, mxyzp1
      loffp = (l - 1)/mxyp1
      k = l - mxyp1*loffp
      loffp = mz*loffp
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(l)
      nn = min(mx,nx-noffp)
      mm = min(my,nyzp(1)-moffp)
      ll = min(mz,nyzp(2)-loffp)
      ih = 0
      nh = 0
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff(1) + moffp
      edgery = noff(1) + moffp + mm
      edgelz = noff(2) + loffp
      edgerz = noff(2) + loffp + ll
c clear counters
      do 10 j = 1, 26
      ncl(j,l) = 0
   10 continue
c loop over particles in tile
      do 20 j = 1, nppp
      dx = ppart(1,j,l)
      dy = ppart(2,j,l)
      dz = ppart(3,j,l)
c find particles going out of bounds
      ist = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(1,j,l) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(1,j,l) = dx
         else
            ist = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(2,j,l) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(2,j,l) = dy
         else
            ist = ist + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) ppart(3,j,l) = dz - anz
         ist = ist + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               ist = ist + 9
            else
               dz = 0.0
            endif
            ppart(3,j,l) = dz
         else
            ist = ist + 9
         endif
      endif
      if (ist.gt.0) then
         ncl(ist,l) = ncl(ist,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j
            ihole(2,ih+1,l) = ist
         else
            nh = 1
         endif
      endif
   20 continue
c set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
   30 continue
!$OMP END PARALLEL DO
c ihole overflow
      if (irc.gt.0) return
c
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,isum,ist,nh,ip,j1,ii)
      do 70 l = 1, mxyzp1
c find address offset for ordered ppbuff array
      isum = 0
      do 40 j = 1, 26
      ist = ncl(j,l)
      ncl(j,l) = isum
      isum = isum + ist
   40 continue
      nh = ihole(1,1,l)
      ip = 0
c loop over particles leaving tile
      do 60 j = 1, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,l)
      ist = ihole(2,j+1,l)
      ii = ncl(ist,l) + 1
      if (ii.le.npbmx) then
         do 50 i = 1, idimp
         ppbuff(i,ii,l) = ppart(i,j1,l)
   50    continue
      else
         ip = 1
      endif
      ncl(ist,l) = ii
   60 continue
c set error
      if (ip.gt.0) irc = ncl(26,l)
   70 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c buffer particles and their number leaving the node up or down:
c update sbufl(:,1), sbufr(:,1), ncll(:,1), nclr(:,1)
      mxzp1 = mx1*mzp1
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k,l,ll)
      do 80 ll = 1, mxzp1
      l = (ll - 1)/mx1
      k = ll + (mxyp1 - mx1)*l
c going straight up and down
      ncll(1,ll,1,1) = ncl(5,k) - ncl(2,k)
      nclr(1,ll,1,1) = ncl(8,k+kk) - ncl(5,k+kk)
c particles going up and back
      ncll(1,ll,2,1) = ncl(14,k) - ncl(11,k)
c particles going down and back
      nclr(1,ll,2,1) = ncl(17,k+kk) - ncl(14,k+kk)
c particles going up and forward
      ncll(1,ll,3,1) = ncl(23,k) - ncl(20,k)
c particles going down and forward
      nclr(1,ll,3,1) = ncl(26,k+kk) - ncl(23,k+kk)
   80 continue
!$OMP END PARALLEL DO
c perform prefix scan
      kk = 1
   90 if (kk.ge.mxzp1) go to 110
!$OMP PARALLEL DO PRIVATE(ll,ii,nn,mm)
      do 100 ll = 1, mxzp1
      ii = (ll - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + ll + kk
      if (nn.le.mxzp1) then
         ncll(1,nn,1,1) = ncll(1,nn,1,1) + ncll(1,mm+1,1,1)
         nclr(1,nn,1,1) = nclr(1,nn,1,1) + nclr(1,mm+1,1,1)
         ncll(1,nn,2,1) = ncll(1,nn,2,1) + ncll(1,mm+1,2,1)
         nclr(1,nn,2,1) = nclr(1,nn,2,1) + nclr(1,mm+1,2,1)
         ncll(1,nn,3,1) = ncll(1,nn,3,1) + ncll(1,mm+1,3,1)
         nclr(1,nn,3,1) = nclr(1,nn,3,1) + nclr(1,mm+1,3,1)
      endif
  100 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 90
  110 kk = mx1*(myp1 - 1)
      j1 = ncll(1,mxzp1,1,1)
      k1 = j1 + ncll(1,mxzp1,2,1)
      j2 = nclr(1,mxzp1,1,1)
      k2 = j2 + nclr(1,mxzp1,2,1)
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,ii,nn,mm)
      do 300 ll = 1, mxzp1
      l = (ll - 1)/mx1
      k = ll + (mxyp1 - mx1)*l
c particles going straight up
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,ll,1,1) - ii
      do 130 j = 1, min(ii,nbmax-nn)
      do 120 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(2,k),k)
  120 continue
  130 continue
      do 140 i = 1, 3
      ncll(i,ll,1,1) = ncl(i+2,k) - ncl(2,k) + nn
  140 continue
c particles going straight down
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,ll,1,1) - ii
      do 160 j = 1, min(ii,nbmax-mm)
      do 150 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  150 continue
  160 continue
      do 170 i = 1, 3
      nclr(i,ll,1,1) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  170 continue
c particles going up and back
      ii = ncl(14,k) - ncl(11,k)
      nn = j1 + ncll(1,ll,2,1) - ii
      do 190 j = 1, min(ii,nbmax-nn)
      do 180 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(11,k),k)
  180 continue
  190 continue
      do 200 i = 1, 3
      ncll(i,ll,2,1) = ncl(i+11,k) - ncl(11,k) + nn
  200 continue
c particles going down and back
      ii = ncl(17,k+kk) - ncl(14,k+kk)
      mm = j2 + nclr(1,ll,2,1) - ii
      do 220 j = 1, min(ii,nbmax-mm)
      do 210 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(14,k+kk),k+kk)
  210 continue
  220 continue
      do 230 i = 1, 3
      nclr(i,ll,2,1) = ncl(i+14,k+kk) - ncl(14,k+kk) + mm
  230 continue
c particles going up and forward
      ii = ncl(23,k) - ncl(20,k)
      nn = k1 + ncll(1,ll,3,1) - ii
      do 250 j = 1, min(ii,nbmax-nn)
      do 240 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(20,k),k)
  240 continue
  250 continue
      do 260 i = 1, 3
      ncll(i,ll,3,1) = ncl(i+20,k) - ncl(20,k) + nn
  260 continue
c particles going down and forward, to different node
      ii = ncl(26,k+kk) - ncl(23,k+kk)
      mm = k2 + nclr(1,ll,3,1) - ii
      do 280 j = 1, min(ii,nbmax-mm)
      do 270 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(23,k+kk),k+kk)
  270 continue
  280 continue
      do 290 i = 1, 3
      nclr(i,ll,3,1) = ncl(i+23,k+kk) - ncl(23,k+kk) + mm
  290 continue
  300 continue
!$OMP END PARALLEL DO
c
c buffer particles and their number leaving the node back or forward:
c update sbufl(:,2), sbufr(:,2), ncll(:,2), nclr(:,2)
      kk = mxyp1*(mzp1 - 1)
!$OMP PARALLEL DO PRIVATE(k,ll)
      do 310 ll = 1, mxyp1
      k = ll
c going straight back or forward
      ncll(1,ll,1,2) = ncl(11,k) - ncl(8,k)
      nclr(1,ll,1,2) = ncl(20,k+kk) - ncl(17,k+kk)
c particles going back and up
      ncll(1,ll,2,2) = ncl(14,k) - ncl(11,k)
c particles going forward and up
      nclr(1,ll,2,2) = ncl(23,k+kk) - ncl(20,k+kk)
c particles going back and down
      ncll(1,ll,3,2) = ncl(17,k) - ncl(14,k)
c particles going forward and down
      nclr(1,ll,3,2) = ncl(26,k+kk) - ncl(23,k+kk)
  310 continue
!$OMP END PARALLEL DO
c perform prefix scan
      kk = 1
  320 if (kk.ge.mxyp1) go to 340
!$OMP PARALLEL DO PRIVATE(ll,ii,nn,mm)
      do 330 ll = 1, mxyp1
      ii = (ll - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + ll + kk
      if (nn.le.mxyp1) then
         ncll(1,nn,1,2) = ncll(1,nn,1,2) + ncll(1,mm+1,1,2)
         nclr(1,nn,1,2) = nclr(1,nn,1,2) + nclr(1,mm+1,1,2)
         ncll(1,nn,2,2) = ncll(1,nn,2,2) + ncll(1,mm+1,2,2)
         nclr(1,nn,2,2) = nclr(1,nn,2,2) + nclr(1,mm+1,2,2)
         ncll(1,nn,3,2) = ncll(1,nn,3,2) + ncll(1,mm+1,3,2)
         nclr(1,nn,3,2) = nclr(1,nn,3,2) + nclr(1,mm+1,3,2)
      endif
  330 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 320
  340 kk = mxyp1*(mzp1 - 1)
      j1 = ncll(1,mxyp1,1,2)
      k1 = j1 + ncll(1,mxyp1,2,2)
      j2 = nclr(1,mxyp1,1,2)
      k2 = j2 + nclr(1,mxyp1,2,2)
!$OMP PARALLEL DO PRIVATE(i,j,k,ll,ii,nn,mm)
      do 530 ll = 1, mxyp1
      k = ll
c particles going straight up
      ii = ncl(11,k) - ncl(8,k)
      nn = ncll(1,ll,1,2) - ii
      do 360 j = 1, min(ii,nbmax-nn)
      do 350 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(8,k),k)
  350 continue
  360 continue
      do 370 i = 1, 3
      ncll(i,ll,1,2) = ncl(i+8,k) - ncl(8,k) + nn
  370 continue
c particles going straight down
      ii = ncl(20,k+kk) - ncl(17,k+kk)
      mm = nclr(1,ll,1,2) - ii
      do 390 j = 1, min(ii,nbmax-mm)
      do 380 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(17,k+kk),k+kk)
  380 continue
  390 continue
      do 400 i = 1, 3
      nclr(i,ll,1,2) = ncl(i+17,k+kk) - ncl(17,k+kk) + mm
  400 continue
c particles going up and back
      ii = ncl(14,k) - ncl(11,k)
      nn = j1 + ncll(1,ll,2,2) - ii
      do 420 j = 1, min(ii,nbmax-nn)
      do 410 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(11,k),k)
  410 continue
  420 continue
      do 430 i = 1, 3
      ncll(i,ll,2,2) = ncl(i+11,k) - ncl(11,k) + nn
  430 continue
c particles going down and back
      ii = ncl(23,k+kk) - ncl(20,k+kk)
      mm = j2 + nclr(1,ll,2,2) - ii
      do 450 j = 1, min(ii,nbmax-mm)
      do 440 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(20,k+kk),k+kk)
  440 continue
  450 continue
      do 460 i = 1, 3
      nclr(i,ll,2,2) = ncl(i+20,k+kk) - ncl(20,k+kk) + mm
  460 continue
c particles going up and forward
      ii = ncl(17,k) - ncl(14,k)
      nn = k1 + ncll(1,ll,3,2) - ii
      do 480 j = 1, min(ii,nbmax-nn)
      do 470 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(14,k),k)
  470 continue
  480 continue
      do 490 i = 1, 3
      ncll(i,ll,3,2) = ncl(i+14,k) - ncl(14,k) + nn
  490 continue
c particles going down and forward, to different node
      ii = ncl(26,k+kk) - ncl(23,k+kk)
      mm = k2 + nclr(1,ll,3,2) - ii
      do 510 j = 1, min(ii,nbmax-mm)
      do 500 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(23,k+kk),k+kk)
  500 continue
  510 continue
      do 520 i = 1, 3
      nclr(i,ll,3,2) = ncl(i+23,k+kk) - ncl(23,k+kk) + mm
  520 continue
  530 continue
!$OMP END PARALLEL DO
c sbufl or sbufr overflow
      kk = max(ncll(3,mxzp1,3,1),nclr(3,mxzp1,3,1))
      ll = max(ncll(3,mxyp1,3,2),nclr(3,mxyp1,3,2))
      ii = max(kk,ll)
c corners overflow
      nn = nclr(3,mx1,2,1) - nclr(3,mxzp1,1,1)
      mm = ncll(3,mx1,2,1) - ncll(3,mxzp1,1,1)
      n = mx1*(mzp1 - 1)
      if (n.gt.0) then
         nr = nclr(3,n,3,1)
         nl = ncll(3,n,3,1)
      else
         nr = nclr(3,mxzp1,2,1)
         nl = ncll(3,mxzp1,2,1)
      endif
      kk = nclr(3,n+mx1,3,1) - nr
      ll = ncll(3,n+mx1,3,1) - nl
c total overflow: result valid only for one processor case
      ii = ii + max(nn+kk,mm+ll)
      if (ii.gt.nbmax) then
         irc = ii
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPORDERF32LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll, 
     1nclr,idimp,nppmx,mx1,myp1,mzp1,mxzyp1,npbmx,ntmax,nbmax,irc)
c this subroutine performs first part of a particle sort by x,y,z grid
c in tiles of mx, my, mz
c linear interpolation, with periodic boundary conditions
c for distributed data, with 2d domain decomposition in y/z.
c tiles are assumed to be arranged in 3D linear memory
c this part of the algorithm has 2 steps.  first, a prefix scan of ncl
c is performed and departing particles are buffered in ppbuff in
c direction order. then, we buffer particles leaving the processor in
c sbufl and sbufr, and store particle number offsets in ncll and nclr.
c it assumes that the number, location, and destination of particles 
c leaving a tile have been previously stored in ncl and ihole by the
c PPGPPUSHF32L subroutine.
c input: all except ppbuff, sbufl, sbufr, ncll, nclr, irc
c output: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c ppart(3,n,m) = position z of particle n in tile m
c ppbuff(i,n,l) = i co-ordinate of particle n in tile l
c sbufl = buffer for particles being sent to lower/back processor
c sbufr = buffer for particles being sent to upper/forward processor
c ncl(i,l) = number of particles going to destination i, tile l
c ihole(1,:,l) = location of hole in array left by departing particle
c ihole(2,:,l) = direction destination of particle leaving hole
c all for tile l
c ihole(1,1,l) = ih, number of holes left (error, if negative)
c ncll = number offset being sent to lower/back processor
c nclr = number offset being sent to upper/forward processor
c noff(1) = lowermost global gridpoint in y in particle partition
c noff(2) = backmost global gridpoint in z in particle partition
c nyzp(1:2) = number of primary (complete) gridpoints in y/z
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c mzp1 = (partition length in z direction - 1)/mz + 1
c mxzyp1 = mx1*max(myp1,mzp1)
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c nbmax =  size of buffers for passing particles between processors
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, myp1, mzp1, mxzyp1, npbmx, ntmax, nbmax
      integer irc
      real ppart, ppbuff, sbufl, sbufr
      integer ncl, ihole, ncll, nclr
      dimension ppart(idimp,nppmx,mx1*myp1*mzp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1*mzp1)
      dimension sbufl(idimp,nbmax,2), sbufr(idimp,nbmax,2)
      dimension ncl(26,mx1*myp1*mzp1), ihole(2,ntmax+1,mx1*myp1*mzp1)
      dimension ncll(3,mxzyp1,3,2), nclr(3,mxzyp1,3,2)
c local data
      integer mxyp1, mxzp1, mxyzp1
      integer i, j, k, l, n, ii, nh, ist, nn, mm, ll, isum
      integer ip, j1, j2, k1, k2, kk, nr, nl
      mxyp1 = mx1*myp1
      mxyzp1 = mxyp1*mzp1
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,isum,ist,nh,ip,j1,ii)
      do 40 l = 1, mxyzp1
c find address offset for ordered ppbuff array
      isum = 0
      do 10 j = 1, 26
      ist = ncl(j,l)
      ncl(j,l) = isum
      isum = isum + ist
   10 continue
      nh = ihole(1,1,l)
      ip = 0
c loop over particles leaving tile
      do 30 j = 1, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,l)
      ist = ihole(2,j+1,l)
      ii = ncl(ist,l) + 1
      if (ii.le.npbmx) then
         do 20 i = 1, idimp
         ppbuff(i,ii,l) = ppart(i,j1,l)
   20    continue
      else
         ip = 1
      endif
      ncl(ist,l) = ii
   30 continue
c set error
      if (ip.gt.0) irc = ncl(26,l)
   40 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c buffer particles and their number leaving the node up or down:
c update sbufl(:,1), sbufr(:,1), ncll(:,1), nclr(:,1)
      mxzp1 = mx1*mzp1
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k,l,ll)
      do 50 ll = 1, mxzp1
      l = (ll - 1)/mx1
      k = ll + (mxyp1 - mx1)*l
c going straight up and down
      ncll(1,ll,1,1) = ncl(5,k) - ncl(2,k)
      nclr(1,ll,1,1) = ncl(8,k+kk) - ncl(5,k+kk)
c particles going up and back
      ncll(1,ll,2,1) = ncl(14,k) - ncl(11,k)
c particles going down and back
      nclr(1,ll,2,1) = ncl(17,k+kk) - ncl(14,k+kk)
c particles going up and forward
      ncll(1,ll,3,1) = ncl(23,k) - ncl(20,k)
c particles going down and forward
      nclr(1,ll,3,1) = ncl(26,k+kk) - ncl(23,k+kk)
   50 continue
!$OMP END PARALLEL DO
c perform prefix scan
      kk = 1
   60 if (kk.ge.mxzp1) go to 80
!$OMP PARALLEL DO PRIVATE(ll,ii,nn,mm)
      do 70 ll = 1, mxzp1
      ii = (ll - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + ll + kk
      if (nn.le.mxzp1) then
         ncll(1,nn,1,1) = ncll(1,nn,1,1) + ncll(1,mm+1,1,1)
         nclr(1,nn,1,1) = nclr(1,nn,1,1) + nclr(1,mm+1,1,1)
         ncll(1,nn,2,1) = ncll(1,nn,2,1) + ncll(1,mm+1,2,1)
         nclr(1,nn,2,1) = nclr(1,nn,2,1) + nclr(1,mm+1,2,1)
         ncll(1,nn,3,1) = ncll(1,nn,3,1) + ncll(1,mm+1,3,1)
         nclr(1,nn,3,1) = nclr(1,nn,3,1) + nclr(1,mm+1,3,1)
      endif
   70 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 60
   80 kk = mx1*(myp1 - 1)
      j1 = ncll(1,mxzp1,1,1)
      k1 = j1 + ncll(1,mxzp1,2,1)
      j2 = nclr(1,mxzp1,1,1)
      k2 = j2 + nclr(1,mxzp1,2,1)
!$OMP PARALLEL DO PRIVATE(i,j,k,l,ll,ii,nn,mm)
      do 270 ll = 1, mxzp1
      l = (ll - 1)/mx1
      k = ll + (mxyp1 - mx1)*l
c particles going straight up
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,ll,1,1) - ii
      do 100 j = 1, min(ii,nbmax-nn)
      do 90 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(2,k),k)
   90 continue
  100 continue
      do 110 i = 1, 3
      ncll(i,ll,1,1) = ncl(i+2,k) - ncl(2,k) + nn
  110 continue
c particles going straight down
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,ll,1,1) - ii
      do 130 j = 1, min(ii,nbmax-mm)
      do 120 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  120 continue
  130 continue
      do 140 i = 1, 3
      nclr(i,ll,1,1) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  140 continue
c particles going up and back
      ii = ncl(14,k) - ncl(11,k)
      nn = j1 + ncll(1,ll,2,1) - ii
      do 160 j = 1, min(ii,nbmax-nn)
      do 150 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(11,k),k)
  150 continue
  160 continue
      do 170 i = 1, 3
      ncll(i,ll,2,1) = ncl(i+11,k) - ncl(11,k) + nn
  170 continue
c particles going down and back
      ii = ncl(17,k+kk) - ncl(14,k+kk)
      mm = j2 + nclr(1,ll,2,1) - ii
      do 190 j = 1, min(ii,nbmax-mm)
      do 180 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(14,k+kk),k+kk)
  180 continue
  190 continue
      do 200 i = 1, 3
      nclr(i,ll,2,1) = ncl(i+14,k+kk) - ncl(14,k+kk) + mm
  200 continue
c particles going up and forward
      ii = ncl(23,k) - ncl(20,k)
      nn = k1 + ncll(1,ll,3,1) - ii
      do 220 j = 1, min(ii,nbmax-nn)
      do 210 i = 1, idimp
      sbufl(i,j+nn,1) = ppbuff(i,j+ncl(20,k),k)
  210 continue
  220 continue
      do 230 i = 1, 3
      ncll(i,ll,3,1) = ncl(i+20,k) - ncl(20,k) + nn
  230 continue
c particles going down and forward, to different node
      ii = ncl(26,k+kk) - ncl(23,k+kk)
      mm = k2 + nclr(1,ll,3,1) - ii
      do 250 j = 1, min(ii,nbmax-mm)
      do 240 i = 1, idimp
      sbufr(i,j+mm,1) = ppbuff(i,j+ncl(23,k+kk),k+kk)
  240 continue
  250 continue
      do 260 i = 1, 3
      nclr(i,ll,3,1) = ncl(i+23,k+kk) - ncl(23,k+kk) + mm
  260 continue
  270 continue
!$OMP END PARALLEL DO
c
c buffer particles and their number leaving the node back or forward:
c update sbufl(:,2), sbufr(:,2), ncll(:,2), nclr(:,2)
      kk = mxyp1*(mzp1 - 1)
!$OMP PARALLEL DO PRIVATE(k,ll)
      do 280 ll = 1, mxyp1
      k = ll
c going straight back or forward
      ncll(1,ll,1,2) = ncl(11,k) - ncl(8,k)
      nclr(1,ll,1,2) = ncl(20,k+kk) - ncl(17,k+kk)
c particles going back and up
      ncll(1,ll,2,2) = ncl(14,k) - ncl(11,k)
c particles going forward and up
      nclr(1,ll,2,2) = ncl(23,k+kk) - ncl(20,k+kk)
c particles going back and down
      ncll(1,ll,3,2) = ncl(17,k) - ncl(14,k)
c particles going forward and down
      nclr(1,ll,3,2) = ncl(26,k+kk) - ncl(23,k+kk)
  280 continue
!$OMP END PARALLEL DO
c perform prefix scan
      kk = 1
  290 if (kk.ge.mxyp1) go to 310
!$OMP PARALLEL DO PRIVATE(ll,ii,nn,mm)
      do 300 ll = 1, mxyp1
      ii = (ll - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + ll + kk
      if (nn.le.mxyp1) then
         ncll(1,nn,1,2) = ncll(1,nn,1,2) + ncll(1,mm+1,1,2)
         nclr(1,nn,1,2) = nclr(1,nn,1,2) + nclr(1,mm+1,1,2)
         ncll(1,nn,2,2) = ncll(1,nn,2,2) + ncll(1,mm+1,2,2)
         nclr(1,nn,2,2) = nclr(1,nn,2,2) + nclr(1,mm+1,2,2)
         ncll(1,nn,3,2) = ncll(1,nn,3,2) + ncll(1,mm+1,3,2)
         nclr(1,nn,3,2) = nclr(1,nn,3,2) + nclr(1,mm+1,3,2)
      endif
  300 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 290
  310 kk = mxyp1*(mzp1 - 1)
      j1 = ncll(1,mxyp1,1,2)
      k1 = j1 + ncll(1,mxyp1,2,2)
      j2 = nclr(1,mxyp1,1,2)
      k2 = j2 + nclr(1,mxyp1,2,2)
!$OMP PARALLEL DO PRIVATE(i,j,k,ll,ii,nn,mm)
      do 500 ll = 1, mxyp1
      k = ll
c particles going straight up
      ii = ncl(11,k) - ncl(8,k)
      nn = ncll(1,ll,1,2) - ii
      do 330 j = 1, min(ii,nbmax-nn)
      do 320 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(8,k),k)
  320 continue
  330 continue
      do 340 i = 1, 3
      ncll(i,ll,1,2) = ncl(i+8,k) - ncl(8,k) + nn
  340 continue
c particles going straight down
      ii = ncl(20,k+kk) - ncl(17,k+kk)
      mm = nclr(1,ll,1,2) - ii
      do 360 j = 1, min(ii,nbmax-mm)
      do 350 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(17,k+kk),k+kk)
  350 continue
  360 continue
      do 370 i = 1, 3
      nclr(i,ll,1,2) = ncl(i+17,k+kk) - ncl(17,k+kk) + mm
  370 continue
c particles going up and back
      ii = ncl(14,k) - ncl(11,k)
      nn = j1 + ncll(1,ll,2,2) - ii
      do 390 j = 1, min(ii,nbmax-nn)
      do 380 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(11,k),k)
  380 continue
  390 continue
      do 400 i = 1, 3
      ncll(i,ll,2,2) = ncl(i+11,k) - ncl(11,k) + nn
  400 continue
c particles going down and back
      ii = ncl(23,k+kk) - ncl(20,k+kk)
      mm = j2 + nclr(1,ll,2,2) - ii
      do 420 j = 1, min(ii,nbmax-mm)
      do 410 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(20,k+kk),k+kk)
  410 continue
  420 continue
      do 430 i = 1, 3
      nclr(i,ll,2,2) = ncl(i+20,k+kk) - ncl(20,k+kk) + mm
  430 continue
c particles going up and forward
      ii = ncl(17,k) - ncl(14,k)
      nn = k1 + ncll(1,ll,3,2) - ii
      do 450 j = 1, min(ii,nbmax-nn)
      do 440 i = 1, idimp
      sbufl(i,j+nn,2) = ppbuff(i,j+ncl(14,k),k)
  440 continue
  450 continue
      do 460 i = 1, 3
      ncll(i,ll,3,2) = ncl(i+14,k) - ncl(14,k) + nn
  460 continue
c particles going down and forward, to different node
      ii = ncl(26,k+kk) - ncl(23,k+kk)
      mm = k2 + nclr(1,ll,3,2) - ii
      do 480 j = 1, min(ii,nbmax-mm)
      do 470 i = 1, idimp
      sbufr(i,j+mm,2) = ppbuff(i,j+ncl(23,k+kk),k+kk)
  470 continue
  480 continue
      do 490 i = 1, 3
      nclr(i,ll,3,2) = ncl(i+23,k+kk) - ncl(23,k+kk) + mm
  490 continue
  500 continue
!$OMP END PARALLEL DO
c sbufl or sbufr overflow
      kk = max(ncll(3,mxzp1,3,1),nclr(3,mxzp1,3,1))
      ll = max(ncll(3,mxyp1,3,2),nclr(3,mxyp1,3,2))
      ii = max(kk,ll)
c corners overflow
      nn = nclr(3,mx1,2,1) - nclr(3,mxzp1,1,1)
      mm = ncll(3,mx1,2,1) - ncll(3,mxzp1,1,1)
      n = mx1*(mzp1 - 1)
      if (n.gt.0) then
         nr = nclr(3,n,3,1)
         nl = ncll(3,n,3,1)
      else
         nr = nclr(3,mxzp1,2,1)
         nl = ncll(3,mxzp1,2,1)
      endif
      kk = nclr(3,n+mx1,3,1) - nr
      ll = ncll(3,n+mx1,3,1) - nl
c total overflow: result valid only for one processor case
      ii = ii + max(nn+kk,mm+ll)
      if (ii.gt.nbmax) then
         irc = ii
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPORDER32LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,  
     1mcll,mclr,mcls,idimp,nppmx,mx1,myp1,mzp1,mxzyp1,npbmx,ntmax,nbmax,
     2irc)
c this subroutine performs second part of a particle sort by x,y,z grid
c in tiles of mx, my, mz
c linear interpolation, with periodic boundary conditions
c for distributed data, with 2d domain decomposition in y/z.
c tiles are assumed to be arranged in 3D linear memory
c incoming particles from other tiles are copied from ppbuff, rbufl, and
c rbufr into ppart
c input: all except ppart, kpic, irc
c output: ppart, kpic, irc
c ppart(1,n,k) = position x of particle n in tile k
c ppart(2,n,k) = position y of particle n in tile k
c ppart(3,n,k) = position z of particle n in tile k
c ppbuff(i,n,k) = i co-ordinate of particle n in tile k
c rbufl = buffer for particles being received from lower/back processor
c rbufr = buffer for particles being received from upper/forward
c processor
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = direction destination of particle leaving hole
c all for tile k
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c mcll = number offset being received from lower/back processor
c mclr = number offset being received from upper/forward processor
c mcls = number ofsets received from corner processors
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c mzp1 = (partition length in z direction - 1)/mz + 1
c mxzyp1 = mx1*max(myp1,mzp1)
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c nbmax =  size of buffers for passing particles between processors
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, myp1, mzp1, mxzyp1, npbmx, ntmax, nbmax
      integer irc
      real ppart, ppbuff, rbufl, rbufr
      integer kpic, ncl, ihole, mcll, mclr, mcls
      dimension ppart(idimp,nppmx,mx1*myp1*mzp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1*mzp1)
      dimension rbufl(idimp,nbmax,2), rbufr(idimp,nbmax,2)
      dimension kpic(mx1*myp1*mzp1), ncl(26,mx1*myp1*mzp1)
      dimension ihole(2,ntmax+1,mx1*myp1*mzp1)
      dimension mcll(3,mxzyp1,3,2), mclr(3,mxzyp1,3,2), mcls(3,mx1+1,4)
c local data
      integer mxyp1, mxyzp1, nppp, ncoff, joff, koff
      integer i, j, k, l, ii, kx, ky, kz, ih, nh, ist, lorr
      integer ip, j1, j2, kxl, kxr, kk, kl, kr, ll, lk, lr, mm, kzs
      logical inside
      integer ks
      dimension ks(26)
      mxyp1 = mx1*myp1
      mxyzp1 = mxyp1*mzp1
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ii,kk,nppp,kx,ky,kz,kl,kr,kxl,kxr,lk,ll,lr,mm,kzs
!$OMP& ,ih,nh,ncoff,joff,koff,ist,j1,j2,ip,lorr,inside,ks)
      do 180 l = 1, mxyzp1
      nppp = kpic(l)
      kz = (l - 1)/mxyp1
      k = l - mxyp1*kz
      kzs = kz*mx1
      kz = kz + 1
c loop over tiles in z
      lk = (kz - 1)*mxyp1
c find tile behind
      ll = (kz - 2)*mxyp1
c find tile in front
      lr = kz*mxyp1
      ky = (k - 1)/mx1 + 1
c loop over tiles in y
      kk = (ky - 1)*mx1
c find tile above
      kl = (ky - 2)*mx1
c find tile below
      kr = ky*mx1
c loop over tiles in x, assume periodic boundary conditions
      kx = k - (ky - 1)*mx1
      kxl = kx - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = kx + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
c find tile number for different directions
      ks(1) = kxr + kk + lk
      ks(2) = kxl + kk + lk
      if (ky.eq.myp1) then
         ks(3) = -kx
         ks(4) = -kxr
         ks(5) = -kxl
      else
         ks(3) = kx + kr + lk
         ks(4) = kxr + kr + lk
         ks(5) = kxl + kr + lk
      endif
      if (ky.eq.1) then
         ks(6) = -kx
         ks(7) = -kxr 
         ks(8) = -kxl
      else
         ks(6) = kx + kl + lk
         ks(7) = kxr + kl + lk
         ks(8) = kxl + kl + lk
      endif
      if (kz.eq.mzp1) then
         ks(9) = -kx
         ks(10) = -kxr
         ks(11) = -kxl
      else
         ks(9) = kx + kk + lr
         ks(10) = kxr + kk + lr
         ks(11) = kxl + kk + lr
      endif
      if ((ky.eq.myp1).or.(kz.eq.mzp1)) then
         ks(12) = -kx
         ks(13) = -kxr
         ks(14) = -kxl
      else
         ks(12) = kx + kr + lr
         ks(13) = kxr + kr + lr
         ks(14) = kxl + kr + lr
      endif
      if ((ky.eq.1).or.(kz.eq.mzp1)) then
         ks(15) = -kx
         ks(16) = -kxr
         ks(17) = -kxl
      else
         ks(15) = kx + kl + lr
         ks(16) = kxr + kl + lr
         ks(17) = kxl + kl + lr
      endif
      if (kz.eq.1) then
         ks(18) = -kx
         ks(19) = -kxr 
         ks(20) = -kxl
      else
         ks(18) = kx + kk + ll
         ks(19) = kxr + kk + ll
         ks(20) = kxl + kk + ll
      endif
      if ((ky.eq.myp1).or.(kz.eq.1)) then
         ks(21) = -kx
         ks(22) = -kxr
         ks(23) = -kxl
      else
         ks(21) = kx + kr + ll
         ks(22) = kxr + kr + ll
         ks(23) = kxl + kr + ll
      endif
      if ((ky.eq.1).or.(kz.eq.1)) then
         ks(24) = -kx
         ks(25) = -kxr
         ks(26) = -kxl
      else
         ks(24) = kx + kl + ll
         ks(25) = kxr + kl + ll
         ks(26) = kxl + kl + ll
      endif
c identify interior
      if ((ky.gt.1).and.(ky.lt.myp1).and.(kz.gt.1).and.(kz.lt.mzp1))    
     1 then
         inside = .true.
      else
         inside = .false.
      endif
c loop over directions
      nh = ihole(1,1,l)
      joff = 0
      koff = 0
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 150 ii = 1, 26
      lorr = 0
      ip = -1
c ip = number of particles coming from direction ii
c interior layers
      if (inside) then
         if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
         ip = ncl(ii,ks(ii)) - ncoff
c edges
      else
c top layer
         if (ky.eq.1) then
            if ((ii.ge.6).and.(ii.le.8)) then
               lorr = -1
               joff = 0
               mm = kzs - ks(ii)
               if (ii.eq.6) then
                  if (mm.gt.1) joff = mcll(3,mm-1,1,1)
               else
                  joff = mcll(ii-6,mm,1,1)
               endif
               ip = mcll(ii-5,mm,1,1) - joff
            else if ((ii.ge.15).and.(ii.le.17)) then
               lorr = -2
               joff = mcll(3,mx1*mzp1,1,1)
               if (kz.lt.mzp1) then
                  mm = kzs + mx1 - ks(ii)
                  if (ii.eq.15) then
                     if (mm.gt.1) joff = mcll(3,mm-1,2,1)
                  else
                     joff = mcll(ii-15,mm,2,1)
                  endif
                  ip = mcll(ii-14,mm,2,1) - joff
c corner data, (ky=1,kz=mzp1)
               else
                  mm = -ks(ii)
                  if (ii.eq.15) then
                     if (mm.eq.1) joff = mcls(1,mx1+1,1)
                     if (mm.gt.1) joff = mcls(3,mm-1,1)
                  else
                     joff = mcls(ii-15,mm,1)
                  endif
                  ip = mcls(ii-14,mm,1) - joff
               endif
            else if ((ii.ge.24).and.(ii.le.26)) then
               lorr = -3
               joff = mcll(3,mx1*mzp1,2,1)
               if (kz.gt.1) then
                  mm = kzs - mx1 - ks(ii)
                  if (ii.eq.24) then
                     if (mm.gt.1) joff = mcll(3,mm-1,3,1)
                  else
                     joff = mcll(ii-24,mm,3,1)
                  endif
                  ip = mcll(ii-23,mm,3,1) - joff
c corner data, (ky=1,kz=1)
               else
                  mm = -ks(ii)
                  if (ii.eq.24) then
                     if (mm.eq.1) joff = mcls(1,mx1+1,2)
                     if (mm.gt.1) joff = mcls(3,mm-1,2)
                  else
                     joff = mcls(ii-24,mm,2)
                  endif
                  ip = mcls(ii-23,mm,2) - joff
               endif
c internal data
            else
               if (ks(ii).gt.0) then
                  if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
                  ip = ncl(ii,ks(ii)) - ncoff
               endif
            endif
         endif
c bottom layer
         if (ky.eq.myp1) then
            if ((ii.ge.3).and.(ii.le.5)) then
               lorr = 1
               joff = 0
               mm = kzs - ks(ii)
               if (ii.eq.3) then
                  if (mm.gt.1) joff = mclr(3,mm-1,1,1)
               else
                  joff = mclr(ii-3,mm,1,1)
               endif
               ip = mclr(ii-2,mm,1,1) - joff
            else if ((ii.ge.12).and.(ii.le.14)) then
               lorr = 2
               joff = mclr(3,mx1*mzp1,1,1)
               if (kz.lt.mzp1) then
                  mm = kzs + mx1 - ks(ii)
                  if (ii.eq.12) then
                     if (mm.gt.1) joff = mclr(3,mm-1,2,1)
                  else
                     joff = mclr(ii-12,mm,2,1)
                  endif
                  ip = mclr(ii-11,mm,2,1) - joff
c corner data, (ky=myp1,kz=mzp1)
               else
                  mm = -ks(ii)
                  if (ii.eq.12) then
                     if (mm.eq.1) joff = mcls(1,mx1+1,3)
                     if (mm.gt.1) joff = mcls(3,mm-1,3)
                  else
                     joff = mcls(ii-12,mm,3)
                  endif
                  ip = mcls(ii-11,mm,3) - joff
               endif
            else if ((ii.ge.21).and.(ii.le.23)) then
               lorr = 3
               joff = mclr(3,mx1*mzp1,2,1)
               if (kz.gt.1) then
                  mm = kzs - mx1 - ks(ii)
                  if (ii.eq.21) then
                    if (mm.gt.1) joff = mclr(3,mm-1,3,1)
                  else
                     joff = mclr(ii-21,mm,3,1)
                  endif
                  ip = mclr(ii-20,mm,3,1) - joff
c corner data, (ky=myp1,kz=1)
               else
                  mm = -ks(ii)
                  if (ii.eq.21) then
                     if (mm.eq.1) joff = mcls(1,mx1+1,4)
                     if (mm.gt.1) joff = mcls(3,mm-1,4)
                  else
                     joff = mcls(ii-21,mm,4)
                  endif
                  ip = mcls(ii-20,mm,4) - joff
               endif
c internal data
            else
               if (ks(ii).gt.0) then
                  if (ky.gt.1) then
                     if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
                     ip = ncl(ii,ks(ii)) - ncoff
                  endif
               endif
            endif
         endif
c front layer
         if (kz.eq.1) then
            if ((ii.ge.18).and.(ii.le.20)) then
               koff = 0
               mm = kk - ks(ii)
               if (ii.eq.18) then
                  if (mm.gt.1) koff = mcll(3,mm-1,1,2)
               else
                  koff = mcll(ii-18,mm,1,2)
               endif
               ip = mcll(ii-17,mm,1,2) - koff
            else if ((ii.ge.21).and.(ii.le.23)) then
               koff = mcll(3,mx1*myp1,1,2)
               if (ky.lt.myp1) then
                  mm = kr - ks(ii)
                  if (ii.eq.21) then
                     if (mm.gt.1) koff = mcll(3,mm-1,2,2)
                  else
                     koff = mcll(ii-21,mm,2,2)
                  endif
                  ip = mcll(ii-20,mm,2,2) - koff
c corner data, already done
c              else
               endif
            else if ((ii.ge.24).and.(ii.le.26)) then
               koff = mcll(3,mx1*myp1,2,2)
               if (ky.gt.1) then
                  mm = kl - ks(ii)
                  if (ii.eq.24) then
                     if (mm.gt.1) koff = mcll(3,mm-1,3,2)
                  else
                     koff = mcll(ii-24,mm,3,2)
                  endif
                  ip = mcll(ii-23,mm,3,2) - koff
c corner data, already done
c              else
               endif
c internal data
            else
               if (ks(ii).gt.0) then
                  if ((ky.gt.1).and.(ky.lt.myp1)) then
                     if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
                     ip = ncl(ii,ks(ii)) - ncoff
                  endif
               endif
            endif
         endif
c back layer
         if (kz.eq.mzp1) then
            if ((ii.ge.9).and.(ii.le.11)) then
               koff = 0
               mm = kk - ks(ii)
               if (ii.eq.9) then
                  if (mm.gt.1) koff = mclr(3,mm-1,1,2)
               else
                  koff = mclr(ii-9,mm,1,2)
               endif
               ip = mclr(ii-8,mm,1,2) - koff
            else if ((ii.ge.12).and.(ii.le.14)) then
               koff = mclr(3,mx1*myp1,1,2)
               if (ky.lt.myp1) then
                  mm = kr - ks(ii)
                  if (ii.eq.12) then
                     if (mm.gt.1) koff = mclr(3,mm-1,2,2)
                  else
                     koff = mclr(ii-12,mm,2,2)
                  endif
                  ip = mclr(ii-11,mm,2,2) - koff
c corner data, already done
c              else
               endif
            else if ((ii.ge.15).and.(ii.le.17)) then
               koff = mclr(3,mx1*myp1,2,2)
               if (ky.gt.1) then
                  mm = kl - ks(ii)
                  if (ii.eq.15) then
                     if (mm.gt.1) koff = mclr(3,mm-1,3,2)
                  else
                     koff = mclr(ii-15,mm,3,2)
                  endif
                  ip = mclr(ii-14,mm,3,2) - koff
c corner data, already done
c              else
               endif
c internal data
            else
               if (ks(ii).gt.0) then
                  if ((ky.gt.1).and.(ky.lt.myp1)) then
                     if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
                     ip = ncl(ii,ks(ii)) - ncoff
                  endif
               endif
            endif
         endif
      endif
c
      if (ip.lt.0) write (*,*) 'help, ip undefined:l,ii=',l,ii
c copy incoming particles
      do 140 j = 1, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,l)
c place overflow at end of array
      else
         j1 = nppp + 1
         nppp = j1
      endif
      if (j1.le.nppmx) then
c interior layers
         if (inside) then
            do 10 i = 1, idimp
            ppart(i,j1,l) = ppbuff(i,j+ncoff,ks(ii))
   10       continue
c edges
         else
c top layer
            if (ky.eq.1) then
c external data
               if (lorr.lt.0) then
                  do 20 i = 1, idimp
                  ppart(i,j1,l) = rbufl(i,j+joff,1)
   20             continue
c internal data
               else if (ks(ii).gt.0) then
                  do 30 i = 1, idimp
                  ppart(i,j1,l) = ppbuff(i,j+ncoff,ks(ii))
   30             continue
               endif
            endif
c bottom layer
            if (ky.eq.myp1) then
c external data
               if (lorr.gt.0) then
                  do 40 i = 1, idimp
                  ppart(i,j1,l) = rbufr(i,j+joff,1)
   40             continue
c internal data
               else if (ks(ii).gt.0) then
                  if (ky.gt.1) then
                     do 50 i = 1, idimp
                     ppart(i,j1,l) = ppbuff(i,j+ncoff,ks(ii))
   50                continue
                  endif
               endif
            endif
c front layer
            if (kz.eq.1) then
              if ((ii.ge.18).and.(ii.le.20)) then
                  do 60 i = 1, idimp
                  ppart(i,j1,l) = rbufl(i,j+koff,2)
   60             continue
               else if ((ii.ge.21).and.(ii.le.23)) then
                  if (ky.lt.myp1) then
                     do 70 i = 1, idimp
                     ppart(i,j1,l) = rbufl(i,j+koff,2)
   70                continue
                  endif
               else if ((ii.ge.24).and.(ii.le.26)) then
                  if (ky.gt.1) then
                     do 80 i = 1, idimp
                     ppart(i,j1,l) = rbufl(i,j+koff,2)
   80                continue
                  endif
c internal data
               else if (ks(ii).gt.0) then
                  if ((ky.gt.1).and.(ky.lt.myp1)) then
                     do 90 i = 1, idimp
                     ppart(i,j1,l) = ppbuff(i,j+ncoff,ks(ii))
   90                continue
                  endif
               endif
            endif
c back layer
            if (kz.eq.mzp1) then
               if ((ii.ge.9).and.(ii.le.11)) then
                  do 100 i = 1, idimp
                  ppart(i,j1,l) = rbufr(i,j+koff,2)
  100             continue
               else if ((ii.ge.12).and.(ii.le.14)) then
                  if (ky.lt.myp1) then
                     do 110 i = 1, idimp
                     ppart(i,j1,l) = rbufr(i,j+koff,2)
  110                continue
                  endif
               else if ((ii.ge.15).and.(ii.le.17)) then
                  if (ky.gt.1) then
                     do 120 i = 1, idimp
                     ppart(i,j1,l) = rbufr(i,j+koff,2)
  120                continue
                  endif
c internal data
               else if (ks(ii).gt.0) then
                  if ((ky.gt.1).and.(ky.lt.myp1)) then
                     do 130 i = 1, idimp
                     ppart(i,j1,l) = ppbuff(i,j+ncoff,ks(ii))
  130                continue
                  endif
               endif
            endif
         endif
      else
         ist = 1
      endif
  140 continue
  150 continue
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
      if (ih.lt.nh) then
         ip = nh - ih
         do 170 j = 1, ip
         j1 = nppp - j + 1
         j2 = ihole(1,nh-j+2,l)
         if (j1.gt.j2) then
c move particle only if it is below current hole
            do 160 i = 1, idimp
            ppart(i,j2,l) = ppart(i,j1,l)
  160       continue
         endif
  170    continue
         nppp = nppp - ip
      endif
      kpic(l) = nppp
  180 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPCGUARD32XL(fxyz,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
c replicate extended periodic vector field in x direction
c linear interpolation, for distributed data with 2D decomposition
c nyzp(1:2) = number of primary (complete) gridpoints in y/z
c nx = system length in xz direction
c ndim = leading dimension of field array fxyz
c nxe = first dimension of field array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition = 2
      implicit none
      integer nx, ndim, nxe, nypmx, nzpmx, idds
      integer nyzp
      real fxyz
      dimension fxyz(ndim,nxe,nypmx,nzpmx), nyzp(idds)
c local data
      integer i, k, l, myp1, mzp1
c replicate edges of extended field
      myp1 = nyzp(1) + 1
      mzp1 = nyzp(2) + 1
!$OMP PARALLEL DO PRIVATE(i,k,l)
      do 30 l = 1, mzp1
      do 20 k = 1, myp1
      do 10 i = 1, ndim
      fxyz(i,nx+1,k,l) = fxyz(i,1,k,l)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPAGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds)
c accumulate extended periodic scalar field in x direction
c linear interpolation, for distributed data with 2D decomposition
c nyzp(1:2) = number of primary (complete) gridpoints in y/z
c nx = system length in xz direction
c nxe = first dimension of field array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition = 2
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds
      integer nyzp
      real q
      dimension q(nxe,nypmx,nzpmx), nyzp(idds)
      integer k, l, myp1, mzp1
c accumulate edges of extended field
      myp1 = nyzp(1) + 1
      mzp1 = nyzp(2) + 1
!$OMP PARALLEL DO PRIVATE(k,l)
      do 20 l = 1, mzp1
      do 10 k = 1, myp1
      q(1,k,l) = q(1,k,l) + q(nx+1,k,l)
      q(nx+1,k,l) = 0.0
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPOIS332(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,  
     1kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions for distributed data,
c with 2D spatial decomposition and OpenMP
c for isign = 0, input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,nvpy,nvpz,
c                       kxyp,kyzp,nzhd
c output: ffc
c for isign =/ 0, input: q,ffc,isign,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,
c                        kyzp,nzhd
c output: fxyz,we
c approximate flop count is:
c 62*nxc*nyc*nzc + 33*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the equation used is:
c fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
c fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
c fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
c fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
c q(l,j,k) = complex charge density for fourier mode jj-1,kk-1,l-1
c fxyz(1,l,j,k) = x component of force/charge
c fxyz(2,l,j,k) = y component of force/charge
c fxyz(3,l,j,k) = z component of force/charge
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
c kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
c if isign = 0, form factor array is prepared
c aimag(ffc(l,j,k)) = finite-size particle shape factor s
c real(ffc(l,j,k)) = potential green's function g
c for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
c kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer isign, nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      integer nzhd
      real ax, ay, az, affp, we
      complex q, fxyz, ffc
      dimension q(nzv,kxyp,kyzp), fxyz(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
c local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz, at1, at2, at3, at4, at5, at6
      complex zero, zt1, zt2
      double precision wp, sum1, sum2
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
c find processor id and offsets in y/z
c js/ks = processor co-ordinates in x/y => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      if (isign.ne.0) go to 40
      if (kstrt.gt.(nvpy*nvpz)) return
c prepare form factor array
      do 30 k = 1, kyzps
      k1 = k + koff
      if (k1.gt.nyh) k1 = k1 - ny
      dky = dny*real(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyps
      dkx = dnx*real(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nzh
      dkz = dnz*real(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.0) then
         ffc(l,j,k) = cmplx(affp,1.0)
      else
         ffc(l,j,k) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 sum1 = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 160
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,at1,at2,at3,at4,zt1,zt2,wp)
!$OMP& REDUCTION(+:sum1)
      do 60 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 50 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,j,k))*aimag(ffc(l,j,k))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*real(l - 1)*at1
            zt1 = cmplx(aimag(q(l,j,k)),-real(q(l,j,k)))
            zt2 = cmplx(aimag(q(l1,j,k)),-real(q(l1,j,k)))
            fxyz(1,l,j,k) = at2*zt1
            fxyz(2,l,j,k) = at3*zt1
            fxyz(3,l,j,k) = at4*zt1
            fxyz(1,l1,j,k) = at2*zt2
            fxyz(2,l1,j,k) = at3*zt2
            fxyz(3,l1,j,k) = -at4*zt2
            wp = wp + at1*(q(l,j,k)*conjg(q(l,j,k))
     1              + q(l1,j,k)*conjg(q(l1,j,k)))
   50       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = real(ffc(1,j,k))*aimag(ffc(1,j,k))
            at2 = dkx*at1
            at3 = dky*at1
            zt1 = cmplx(aimag(q(1,j,k)),-real(q(1,j,k)))
            fxyz(1,1,j,k) = at2*zt1
            fxyz(2,1,j,k) = at3*zt1
            fxyz(3,1,j,k) = zero
            fxyz(1,l1,j,k) = zero
            fxyz(2,l1,j,k) = zero
            fxyz(3,l1,j,k) = zero
            wp = wp + at1*(q(1,j,k)*conjg(q(1,j,k)))
         endif
      endif
      sum1 = sum1 + wp
   60 continue
!$OMP END DO
!$OMP END PARALLEL
c mode numbers kx = 0, nx/2
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,at1,at3,at4,zt1,zt2,wp)
!$OMP& REDUCTION(+:sum2)
      do 90 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         wp = 0.0d0
         if (js.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 70 l = 2, nzh
               l1 = nz2 - l
               at1 = real(ffc(l,1,k))*aimag(ffc(l,1,k))
               at3 = dky*at1
               at4 = dnz*real(l - 1)*at1
               zt1 = cmplx(aimag(q(l,1,k)),-real(q(l,1,k)))
               zt2 = cmplx(aimag(q(l1,1,k)),-real(q(l1,1,k)))
               fxyz(1,l,1,k) = zero
               fxyz(2,l,1,k) = at3*zt1
               fxyz(3,l,1,k) = at4*zt1
               fxyz(1,l1,1,k) = zero
               fxyz(2,l1,1,k) = at3*zt2
               fxyz(3,l1,1,k) = -at4*zt2
               wp = wp + at1*(q(l,1,k)*conjg(q(l,1,k))
     1                 + q(l1,1,k)*conjg(q(l1,1,k)))
   70          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = real(ffc(1,1,k))*aimag(ffc(1,1,k))
               at3 = dky*at1
               zt1 = cmplx(aimag(q(1,1,k)),-real(q(1,1,k)))
               fxyz(1,1,1,k) = zero
               fxyz(2,1,1,k) = at3*zt1
               fxyz(3,1,1,k) = zero
               fxyz(1,l1,1,k) = zero
               fxyz(2,l1,1,k) = zero
               fxyz(3,l1,1,k) = zero
               wp = wp + at1*(q(1,1,k)*conjg(q(1,1,k)))
c throw away kx = nx/2
            else
               do 80 l = 1, nz
               fxyz(1,l,1,k) = zero
               fxyz(2,l,1,k) = zero
               fxyz(3,l,1,k) = zero
   80          continue
            endif
         endif
         sum2 = sum2 + wp
      endif
   90 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum2
c mode numbers ky = 0, ny/2
      sum2 = 0.0d0
c keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,at1,at2,at4,zt1,zt2,wp)
!$OMP& REDUCTION(+:sum2)
         do 110 j = 1, kxyps
         dkx = dnx*real(j + joff)
         wp = 0.0d0
         if ((j+joff).gt.0) then
            do 100 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,j,1))*aimag(ffc(l,j,1))
            at2 = dkx*at1
            at4 = dnz*real(l - 1)*at1
            zt1 = cmplx(aimag(q(l,j,1)),-real(q(l,j,1)))
            zt2 = cmplx(aimag(q(l1,j,1)),-real(q(l1,j,1)))
            fxyz(1,l,j,1) = at2*zt1
            fxyz(2,l,j,1) = zero
            fxyz(3,l,j,1) = at4*zt1
            fxyz(1,l1,j,1) = at2*zt2
            fxyz(2,l1,j,1) = zero
            fxyz(3,l1,j,1) = -at4*zt2
            wp = wp + at1*(q(l,j,1)*conjg(q(l,j,1))
     1              + q(l1,j,1)*conjg(q(l1,j,1)))
  100       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = real(ffc(1,j,1))*aimag(ffc(1,j,1))
            at2 = dkx*at1
            zt1 = cmplx(aimag(q(1,j,1)),-real(q(1,j,1)))
            fxyz(1,1,j,1) = at2*zt1
            fxyz(2,1,j,1) = zero
            fxyz(3,1,j,1) = zero
            fxyz(1,l1,j,1) = zero
            fxyz(2,l1,j,1) = zero
            fxyz(3,l1,j,1) = zero
            wp = wp + at1*(q(1,j,1)*conjg(q(1,j,1)))
         endif
         sum2 = sum2 + wp
  110    continue
!$OMP END PARALLEL DO
         wp = 0.0d0
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 120 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,1,1))*aimag(ffc(l,1,1))
            at4 = dnz*real(l - 1)*at1
            zt1 = cmplx(aimag(q(l,1,1)),-real(q(l,1,1)))
            fxyz(1,l,1,1) = zero
            fxyz(2,l,1,1) = zero
            fxyz(3,l,1,1) = at4*zt1
            fxyz(1,l1,1,1) = zero
            fxyz(2,l1,1,1) = zero
            fxyz(3,l1,1,1) = zero
            wp = wp + at1*(q(l,1,1)*conjg(q(l,1,1)))
  120       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            fxyz(1,1,1,1) = zero
            fxyz(2,1,1,1) = zero
            fxyz(3,1,1,1) = zero
            fxyz(1,l1,1,1) = zero
            fxyz(2,l1,1,1) = zero
            fxyz(3,l1,1,1) = zero
         endif
         sum2 = sum2 + wp
      endif
c throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 140 j = 1, kxyps
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 130 l = 1, nz
            fxyz(1,l,j,k1) = zero
            fxyz(2,l,j,k1) = zero
            fxyz(3,l,j,k1) = zero
  130       continue
         endif
  140    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 150 l = 1, nz
            fxyz(1,l,1,k1) = zero
            fxyz(2,l,1,k1) = zero
            fxyz(3,l,1,k1) = zero
  150       continue
         endif
      endif
      sum1 = sum1 + sum2
  160 continue
      we = real(nx)*real(ny)*real(nz)*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT32RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
c this subroutine calculates tables needed by a three dimensional
c real to complex fast fourier transform and its inverse.
c input: indx, indy, indz, nxhyzd, nxyzhd
c output: mixup, sct
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, indz, nxhyzd, nxyzhd
      integer mixup
      complex sct
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, ny, nz, nxyz, nxhyz, nxyzh
      integer j, k, lb, ll, jb, it
      real dnxyz, arg
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhyz
      lb = j - 1
      ll = 0
      do 10 k = 1, ndx1yz
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles 2*n*pi/nxyz
      nxyzh = nxyz/2
      dnxyz = 6.28318530717959/real(nxyz)
      do 30 j = 1, nxyzh
      arg = dnxyz*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPPFFT32RM(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx,
     1indy,indz,kstrt,nvpy,nvpz,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,   
     2kypd,kyzpd,kzpd,kzyp,nxhyzd,nxyzhd)
c wrapper function for 3d real to complex fft, with packed data
c parallelized with MPI/OpenMP
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nvpy, nvpz
      integer nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
      integer kxypd, kypd, kyzpd, kzpd, kzyp, nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f, g, h, bs, br, sct
      dimension f(nxvh,kypd,kzpd), g(nyv,kxypd,kzpd), h(nzv,kxypd,kyzpd)
      dimension bs(kxyp*kzyp,kzp), br(kxyp*kzyp,kzp)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer nxh, ny, nz, kypi, kxypi, js, ks, kxypp, kypp, kzpp, nvp
      real tp, tf
      double precision dtime
      data kypi, kxypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      nz = 2**indz
c js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxypp = min(kxyp,max(0,nxh-kxyp*js))
      kypp = min(kyp,max(0,ny-kyp*js))
      kzpp = min(kzp,max(0,nz-kzp*ks))
      nvp = nvpy*nvpz
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PPFFT32RMXX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp,   
     1kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOS3A(f,g,bs,br,nxh,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,nxvh,
     1nyv,kxypd,kypd,kzpd)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PPFFT32RMXY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy,  
     1nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PPTPOS3B(g,h,bs,br,nxh,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,nvpz
     1,nyv,nzv,kxypd,kyzpd,kzpd)
         call PWTIMERA(1,tp,dtime)
c perform z fft
         call PPFFT32RMXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy,  
     1nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
c transpose h array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOS3B(h,g,br,bs,nxh,nz,ny,kxyp,kzp,kyzp,kstrt,nvpy, 
     1nvpz,nzv,nyv,kxypd,kzpd,kyzpd)
            call PPTPOS3A(g,f,br,bs,ny,nxh,nz,kyp,kxyp,kzp,kstrt,nvpy,  
     1nyv,nxvh,kypd,kxypd,kzpd)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOS3A(f,g,bs,br,nxh,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,  
     1nxvh,nyv,kxypd,kypd,kzpd)
            call PPTPOS3B(g,h,bs,br,nxh,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy, 
     1nvpz,nyv,nzv,kxypd,kyzpd,kzpd)
            call PWTIMERA(1,tf,dtime)
         endif
c perform z fft
         call PPFFT32RMXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy,  
     1nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PPTPOS3B(h,g,br,bs,nxh,nz,ny,kxyp,kzp,kyzp,kstrt,nvpy,nvpz
     1,nzv,nyv,kxypd,kzpd,kyzpd)
         call PWTIMERA(1,tp,dtime)
c perform y fft
         call PPFFT32RMXY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy,  
     1nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOS3A(g,f,br,bs,ny,nxh,nz,kyp,kxyp,kzp,kstrt,nvpy,nyv, 
     1nxvh,kypd,kxypd,kzpd)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PPFFT32RMXX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp,   
     1kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPPFFT32RM3(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx
     1,indy,indz,kstrt,nvpy,nvpz,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,  
     2kypd,kyzpd,kzpd,kzyp,nxhyzd,nxyzhd)
c wrapper function for 3 3d real to complex ffts, with packed data
c c parallelized with MPI/OpenMP
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nvpy, nvpz
      integer nxvh, nyv, nzv, kxyp, kyp, kyzp, kzp
      integer kxypd, kypd, kyzpd, kzpd, kzyp, nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f, g, h, bs, br, sct
      dimension f(3,nxvh,kypd,kzpd), g(3,nyv,kxypd,kzpd)
      dimension h(3,nzv,kxypd,kyzpd)
      dimension bs(3,kxyp*kzyp,kzp), br(3,kxyp*kzyp,kzp)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer nxh, ny, nz, kypi, kxypi, js, ks, kxypp, kypp, kzpp, nvp
      real tp, tf
      double precision dtime
      data kypi, kxypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      nz = 2**indz
c js/ks = processor co-ordinates in y/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kxypp = min(kxyp,max(0,nxh-kxyp*js))
      kypp = min(kyp,max(0,ny-kyp*js))
      kzpp = min(kzp,max(0,nz-kzp*ks))
      nvp = nvpy*nvpz
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PPFFT32RM3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp,  
     1kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOS3A(f,g,bs,br,nxh,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,3,  
     1nxvh,nyv,kxypd,kypd,kzpd)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PPFFT32RM3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, 
     1nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PPNTPOS3B(g,h,bs,br,nxh,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,   
     1nvpz,3,nyv,nzv,kxypd,kyzpd,kzpd)
         call PWTIMERA(1,tp,dtime)
c perform z fft
         call PPFFT32RM3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, 
     1nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
c transpose h array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOS3B(h,g,br,bs,nxh,nz,ny,kxyp,kzp,kyzp,kstrt,nvpy,
     1nvpz,3,nzv,nyv,kxypd,kzpd,kyzpd)
            call PPNTPOS3A(g,f,br,bs,ny,nxh,nz,kyp,kxyp,kzp,kstrt,nvpy,3
     1,nyv,nxvh,kypd,kxypd,kzpd)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOS3A(f,g,bs,br,nxh,ny,nz,kxyp,kyp,kzp,kstrt,nvpy,3
     1,nxvh,nyv,kxypd,kypd,kzpd)
            call PPNTPOS3B(g,h,bs,br,nxh,ny,nz,kxyp,kyzp,kzp,kstrt,nvpy,
     1nvpz,3,nyv,nzv,kxypd,kyzpd,kzpd)
            call PWTIMERA(1,tf,dtime)
         endif
c perform z fft
         call PPFFT32RM3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, 
     1nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PPNTPOS3B(h,g,br,bs,nxh,nz,ny,kxyp,kzp,kyzp,kstrt,nvpy,   
     1nvpz,3,nzv,nyv,kxypd,kzpd,kyzpd)
         call PWTIMERA(1,tp,dtime)
c perform y fft
         call PPFFT32RM3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy, 
     1nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOS3A(g,f,br,bs,ny,nxh,nz,kyp,kxyp,kzp,kstrt,nvpy,3,  
     1nyv,nxvh,kypd,kxypd,kzpd)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PPFFT32RM3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp,  
     1kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT32RMXX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp,
     1kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
c this subroutine performs the x part of a three dimensional real to
c complex fast fourier transform and its inverse for a subset of y and z
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition and OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c f(n,k,i) = (1/nx*ny*nz)*sum(f(j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, a forward fourier transform is performed
c f(n,k,i) = sum(f(j,k,i)*exp(sqrt(-1)*2pi*n*j/nx))
c kstrt = starting data block number
c nvp = number of real or virtual processors
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kzpp = number of z indices used
c kypd = second dimension of f
c kzpd = third dimension of f
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c final fourier coefficients are stored as follows:
c h(l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js, kk = k + kyzp*ks
c and MPI rank idproc = js + nvpy*ks
c 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
c h(l,1,k) = mode nx/2,kk-1,l-1, where ny/2+2 <= kk <= ny, 1 <= l <= nz,
c the following are located on node js = 0 and ks = 0:
c h(l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
c imag(h(1,1,1)) = real part of mode nx/2,0,0
c imag(h(nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
c the following are located on node js = 0 and ks = nyh/kyzp:
c h(l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
c imag(h(1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
c imag(h(nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvp, kypi, kypp, nxvh
      integer kzpp, kypd, kzpd, nxhyzd, nxyzhd
      integer mixup
      complex f, sct
      dimension f(nxvh,kypd,kzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nz, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, nrxb
      integer nrx, kypt, nn
      real ani
      complex s, t, t1
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kypt = kypi + kypp - 1
      if (kstrt.gt.nvp) return
      if (isign.gt.0) go to 70
c inverse fourier transform
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t,t1)
      do 60 nn = 1, kypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kypi
c bit-reverse array elements in x
      do 10 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t = f(j1,i,n)
         f(j1,i,n) = f(j,i,n)
         f(j,i,n) = t
      endif
   10 continue
c then transform in x
      do 40 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*f(j2,i,n)
      f(j2,i,n) = f(j1,i,n) - t
      f(j1,i,n) = f(j1,i,n) + t
   20 continue
   30 continue
   40 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      do 50 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,i,n))
      s = f(j,i,n) + t
      t = (f(j,i,n) - t)*t1
      f(j,i,n) = ani*(s + t)
      f(nxh2-j,i,n) = ani*conjg(s - t)
   50 continue
      f(1,i,n) = 2.0*ani*cmplx(real(f(1,i,n)) + aimag(f(1,i,n)),
     1                         real(f(1,i,n)) - aimag(f(1,i,n)))
      if (nxhh.gt.0) f(nxhh+1,i,n) = 2.0*ani*conjg(f(nxhh+1,i,n))
   60 continue
!$OMP END PARALLEL DO
      return
c forward fourier transform
   70 nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t,t1)
      do 130 nn = 1, kypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kypi
c scramble coefficients
      kmr = nxyz/nx
      do 80 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,i,n))
      s = f(j,i,n) + t
      t = (f(j,i,n) - t)*t1
      f(j,i,n) = s + t
      f(nxh2-j,i,n) = conjg(s - t)
   80 continue
      f(1,i,n) = cmplx(real(f(1,i,n)) + aimag(f(1,i,n)),
     1                 real(f(1,i,n)) - aimag(f(1,i,n)))
      if (nxhh.gt.0) f(nxhh+1,i,n) = 2.0*conjg(f(nxhh+1,i,n))
c bit-reverse array elements in x
      do 90 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t = f(j1,i,n)
         f(j1,i,n) = f(j,i,n)
         f(j,i,n) = t
      endif
   90 continue
c then transform in x
      do 120 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*f(j2,i,n)
      f(j2,i,n) = f(j1,i,n) - t
      f(j1,i,n) = f(j1,i,n) + t
  100 continue
  110 continue
  120 continue
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT32RMXY(g,isign,mixup,sct,indx,indy,indz,kstrt,nvpy
     1,nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
c this subroutine performs the y part of a three dimensional real to
c complex fast fourier transform and its inverse for a subset of x and z
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition and OpenMP
c for isign = (-1,1), input: all, output: g
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c g(m,j,i) = sum(g(k,j,i)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(m,j,i) = sum(g(k,j,i)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c kxypi = initial x index used
c kxypp = number of x indices used
c nyv = first dimension of g
c kzpp = number of z indices used
c kxypd = second dimension of g
c kzpd = third dimension of g
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c final fourier coefficients are stored as follows:
c h(l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js, kk = k + kyzp*ks
c and MPI rank idproc = js + nvpy*ks
c 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
c h(l,1,k) = mode nx/2,kk-1,l-1, where ny/2+2 <= kk <= ny, 1 <= l <= nz,
c the following are located on node js = 0 and ks = 0:
c h(l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
c imag(h(1,1,1)) = real part of mode nx/2,0,0
c imag(h(nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
c the following are located on node js = 0 and ks = nyh/kyzp:
c h(l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
c imag(h(1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
c imag(h(nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvpy, nvpz, kxypi, kxypp
      integer nyv, kzpp, kxypd, kzpd, nxhyzd, nxyzhd
      integer mixup
      complex g, sct
      dimension g(nyv,kxypd,kzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, ny2, nz, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, nryb
      integer js, ks, nry, kxypt, nn
      complex s, t
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = max(1,ny/2)
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxypt = kxypi + kxypp - 1
c js/ks = processor co-ordinates in x/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      if (kstrt.gt.(nvpy*nvpz)) return
      if (isign.gt.0) go to 80
c inverse fourier transform
      nryb = nxhyz/ny
      nry = nxyz/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t)
      do 50 nn = 1, kxypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kxypi
c bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t = g(k1,i,n)
         g(k1,i,n) = g(k,i,n)
         g(k,i,n) = t
      endif
   10 continue
c then transform in y
      do 40 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*g(j2,i,n)
      g(j2,i,n) = g(j1,i,n) - t
      g(j1,i,n) = g(j1,i,n) + t
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
!$OMP PARALLEL DO PRIVATE(k,n,s)
         do 70 n = 1, kzpp
         do 60 k = 2, nyh
         s = g(ny2-k,1,n)
         g(ny2-k,1,n) = 0.5*cmplx(aimag(g(k,1,n) + s),
     1                            real(g(k,1,n) - s))
         g(k,1,n) = 0.5*cmplx(real(g(k,1,n) + s),aimag(g(k,1,n) - s))
   60    continue
   70    continue
!$OMP END PARALLEL DO
      endif
      return
c forward fourier transform
   80 nryb = nxhyz/ny
      nry = nxyz/ny
c scramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
!$OMP PARALLEL DO PRIVATE(k,n,s)
         do 100 n = 1, kzpp
         do 90 k = 2, nyh
         s = cmplx(aimag(g(ny2-k,1,n)),real(g(ny2-k,1,n)))
         g(ny2-k,1,n) = conjg(g(k,1,n) - s)
         g(k,1,n) = g(k,1,n) + s
   90    continue
  100    continue
!$OMP END PARALLEL DO
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t)
      do 150 nn = 1, kxypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kxypi
c bit-reverse array elements in y
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t = g(k1,i,n)
         g(k1,i,n) = g(k,i,n)
         g(k,i,n) = t
      endif
  110 continue
c then transform in y
      do 140 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*g(j2,i,n)
      g(j2,i,n) = g(j1,i,n) - t
      g(j1,i,n) = g(j1,i,n) + t
  120 continue
  130 continue
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT32RMXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,nvpy
     1,nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
c this subroutine performs the z part of a three dimensional real to
c complex fast fourier transform and its inverse for a subset of x and y
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition and OpenMP
c for isign = (-1,1), input: all, output: h
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c h(l,n,m) = sum(h(i,j,k)*exp(-sqrt(-1)*2pi*l*i/nz))
c if isign = 1, a forward fourier transform is performed
c h(l,n,m) = sum(h(i,j,k)*exp(sqrt(-1)*2pi*ll*ii/nz))
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c kxypi = initial x index used
c kxypp = number of x indices used
c nzv = first dimension of h
c kyzpp = number of y indices used
c kxypd = second dimension of h
c kyzpd = third dimension of h
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c final fourier coefficients are stored as follows:
c h(l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js, kk = k + kyzp*ks
c and MPI rank idproc = js + nvpy*ks
c 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
c h(l,1,k) = mode nx/2,kk-1,l-1, where ny/2+2 <= kk <= ny, 1 <= l <= nz,
c the following are located on node js = 0 and ks = 0:
c h(l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
c imag(h(1,1,1)) = real part of mode nx/2,0,0
c imag(h(nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
c the following are located on node js = 0 and ks = nyh/kyzp:
c h(l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
c imag(h(1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
c imag(h(nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvpy, nvpz, kxypi, kxypp
      integer nzv, kyzp, kxypd, kyzpd, nxhyzd, nxyzhd
      integer mixup
      complex h, sct
      dimension h(nzv,kxypd,kyzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, nz, nzh, nz2, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, nrzb
      integer l1, js, ks, nrz, kxypt, kyzpp, kyzb, nn
      complex s, t
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = max(1,ny/2)
      nz = 2**indz
      nzh = max(1,nz/2)
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxypt = kxypi + kxypp - 1
c js/ks = processor co-ordinates in x/y => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kyzpp = min(kyzp,max(0,ny-kyzp*ks))
      if (kstrt.gt.(nvpy*nvpz)) return
      if (isign.gt.0) go to 80
c inverse fourier transform
      nrzb = nxhyz/nz
      nrz = nxyz/nz
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,l1,s,t)
      do 50 nn = 1, kxypp*kyzpp
      i = (nn - 1)/kyzpp
      n = nn - kyzpp*i
      i = i +  kxypi
c bit-reverse array elements in z
      do 10 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         t = h(l1,i,n)
         h(l1,i,n) = h(l,i,n)
         h(l,i,n) = t
      endif
   10 continue
c finally transform in z
      do 40 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*h(j2,i,n)
      h(j2,i,n) = h(j1,i,n) - t
      h(j1,i,n) = h(j1,i,n) + t
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
         if (ks.eq.0) then
            do 60 n = 2, nzh
            s = h(nz2-n,1,1)
            h(nz2-n,1,1) = 0.5*cmplx(aimag(h(n,1,1) + s),
     1                               real(h(n,1,1) - s))
            h(n,1,1) = 0.5*cmplx(real(h(n,1,1) + s),aimag(h(n,1,1) - s))
   60       continue
         endif
         kyzb = nyh/kyzp
         if (ks.eq.kyzb) then
            k1 = nyh - kyzb*kyzp + 1
            do 70 n = 2, nzh
            s = h(nz2-n,1,k1)
            h(nz2-n,1,k1) = 0.5*cmplx(aimag(h(n,1,k1) + s),
     1                                real(h(n,1,k1) - s))
            h(n,1,k1) = 0.5*cmplx(real(h(n,1,k1) + s),
     1                            aimag(h(n,1,k1) - s))
   70       continue
        endif
      endif
      return
c forward fourier transform
   80 nrzb = nxhyz/nz
      nrz = nxyz/nz
c scramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
         if (ks.eq.0) then
            do 90 n = 2, nzh
            s = cmplx(aimag(h(nz2-n,1,1)),real(h(nz2-n,1,1)))
            h(nz2-n,1,1) = conjg(h(n,1,1) - s)
            h(n,1,1) = h(n,1,1) + s
   90       continue
         endif
         kyzb = nyh/kyzp
         if (ks.eq.kyzb) then
            k1 = nyh - kyzb*kyzp + 1
            do 100 n = 2, nzh
            s = cmplx(aimag(h(nz2-n,1,k1)),real(h(nz2-n,1,k1)))
            h(nz2-n,1,k1) = conjg(h(n,1,k1) - s)
            h(n,1,k1) = h(n,1,k1) + s
  100       continue
         endif
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,l1,s,t)
      do 150 nn = 1, kxypp*kyzpp
      i = (nn - 1)/kyzpp
      n = nn - kyzpp*i
      i = i +  kxypi
c bit-reverse array elements in z
      do 110 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         t = h(l1,i,n)
         h(l1,i,n) = h(l,i,n)
         h(l,i,n) = t
      endif
  110 continue
c first transform in z
      do 140 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*h(j2,i,n)
      h(j2,i,n) = h(j1,i,n) - t
      h(j1,i,n) = h(j1,i,n) + t
  120 continue
  130 continue
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT32RM3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,nvp
     1,kypi,kypp,nxvh,kzpp,kypd,kzpd,nxhyzd,nxyzhd)
c this subroutine performs the x part of 3 three dimensional real to
c complex fast fourier transforms and their inverses for a subset of
c y and z, using complex arithmetic,
c for data which is distributed in blocks, with 2D spatial decomposition
c with OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c f(1:3,n,k,i) = (1/nx*ny*nz)*sum(f(1:3,j,k,i)*
c                                 exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, a forward fourier transform is performed
c f(1:3,n,k,i) = sum(f(1:3,j,k,i)*exp(sqrt(-1)*2pi*n*j/nx))
c kstrt = starting data block number
c nvp = number of real or virtual processors
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kzpp = number of z indices used
c kypd = second dimension of f
c kzpd = third dimension of f
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c final fourier coefficients are stored as follows:
c h(1:3,l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js,
c kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
c 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
c h(1:3,l,1,k) = mode nx/2,kk-1,l-1,
c where ny/2+2 <= kk <= ny, 1 <= l <= nz,
c the following are located on node js = 0 and ks = 0:
c h(1:3,l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
c imag(h(1:3,1,1,1)) = real part of mode nx/2,0,0
c imag(h(1:3,nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
c the following are located on node js = 0 and ks = nyh/kyzp:
c h(1:3,l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
c imag(h(1:3,1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
c imag(h(1:3,nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvp, kypi, kypp, nxvh
      integer kzpp, kypd, kzpd, nxhyzd, nxyzhd
      integer mixup
      complex f, sct
      dimension f(3,nxvh,kypd,kzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nz, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj, nrxb
      integer nrx, kypt, nn
      real ani, at1, at2
      complex s, t, t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kypt = kypi + kypp - 1
      if (kstrt.gt.nvp) return
      if (isign.gt.0) go to 100
c inverse fourier transform
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,at2,s,t,t1,
!$OMP& t2,t3)
      do 90 nn = 1, kypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kypi
c swap complex components
      do 10 j = 1, nxh
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(real(f(2,j,i,n)),aimag(f(3,j,i,n)))
      at2 = aimag(f(2,j,i,n))
      f(2,j,i,n) = cmplx(aimag(f(1,j,i,n)),at1)
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
   10 continue
c bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i,n)
         t2 = f(2,j1,i,n)
         t3 = f(3,j1,i,n)
         f(1,j1,i,n) = f(1,j,i,n)
         f(2,j1,i,n) = f(2,j,i,n)
         f(3,j1,i,n) = f(3,j,i,n)
         f(1,j,i,n) = t1
         f(2,j,i,n) = t2
         f(3,j,i,n) = t3
      endif
   20 continue
c then transform in x
      do 50 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t1 = s*f(1,j2,i,n)
      t2 = s*f(2,j2,i,n)
      t3 = s*f(3,j2,i,n)
      f(1,j2,i,n) = f(1,j1,i,n) - t1
      f(2,j2,i,n) = f(2,j1,i,n) - t2
      f(3,j2,i,n) = f(3,j1,i,n) - t3
      f(1,j1,i,n) = f(1,j1,i,n) + t1
      f(2,j1,i,n) = f(2,j1,i,n) + t2
      f(3,j1,i,n) = f(3,j1,i,n) + t3
   30 continue
   40 continue
   50 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      do 70 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 60 jj = 1, 3
      t = conjg(f(jj,nxh2-j,i,n))
      s = f(jj,j,i,n) + t
      t = (f(jj,j,i,n) - t)*t1
      f(jj,j,i,n) = ani*(s + t)
      f(jj,nxh2-j,i,n) = ani*conjg(s - t)
   60 continue
   70 continue
      do 80 jj = 1, 3
      f(jj,1,i,n) = 
     1             2.0*ani*cmplx(real(f(jj,1,i,n)) + aimag(f(jj,1,i,n)),
     2                           real(f(jj,1,i,n)) - aimag(f(jj,1,i,n)))
      if (nxhh.gt.0) f(jj,nxhh+1,i,n) = 2.0*ani*conjg(f(jj,nxhh+1,i,n))
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
c forward fourier transform
  100 nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,at2,s,t,t1,
!$OMP& t2,t3)
      do 190 nn = 1, kypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kypi
c scramble coefficients
      kmr = nxyz/nx
      do 120 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 jj = 1, 3
      t = conjg(f(jj,nxh2-j,i,n))
      s = f(jj,j,i,n) + t
      t = (f(jj,j,i,n) - t)*t1
      f(jj,j,i,n) = s + t
      f(jj,nxh2-j,i,n) = conjg(s - t)
  110 continue
  120 continue
      do 130 jj = 1, 3
      f(jj,1,i,n) = cmplx(real(f(jj,1,i,n)) + aimag(f(jj,1,i,n)),
     1                    real(f(jj,1,i,n)) - aimag(f(jj,1,i,n)))
      if (nxhh.gt.0) f(jj,nxhh+1,i,n) = 2.0*conjg(f(jj,nxhh+1,i,n))
  130 continue
c bit-reverse array elements in x
      do 140 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i,n)
         t2 = f(2,j1,i,n)
         t3 = f(3,j1,i,n)
         f(1,j1,i,n) = f(1,j,i,n)
         f(2,j1,i,n) = f(2,j,i,n)
         f(3,j1,i,n) = f(3,j,i,n)
         f(1,j,i,n) = t1
         f(2,j,i,n) = t2
         f(3,j,i,n) = t3
      endif
  140 continue
c finally transform in x
      do 170 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 160 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 150 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*f(1,j2,i,n)
      t2 = s*f(2,j2,i,n)
      t3 = s*f(3,j2,i,n)
      f(1,j2,i,n) = f(1,j1,i,n) - t1
      f(2,j2,i,n) = f(2,j1,i,n) - t2
      f(3,j2,i,n) = f(3,j1,i,n) - t3
      f(1,j1,i,n) = f(1,j1,i,n) + t1
      f(2,j1,i,n) = f(2,j1,i,n) + t2
      f(3,j1,i,n) = f(3,j1,i,n) + t3
  150 continue
  160 continue
  170 continue
c swap complex components
      do 180 j = 1, nxh
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(aimag(f(2,j,i,n)),aimag(f(3,j,i,n)))
      at2 = real(f(2,j,i,n))
      f(2,j,i,n) = cmplx(at1,aimag(f(1,j,i,n)))
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT32RM3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,   
     1nvpy,nvpz,kxypi,kxypp,nyv,kzpp,kxypd,kzpd,nxhyzd,nxyzhd)
c this subroutine performs the y part of 3 three dimensional real to
c complex fast fourier transforms and their inverses for a subset of
c x and z, using complex arithmetic,
c for data which is distributed in blocks, with 2D spatial decomposition
c with OpenMP
c for isign = (-1,1), input: all, output: g
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c g(1:3,m,j,i) = sum(g(1:3,k,j,i)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(1:3,m,j,i) = sum(g(1:3,k,j,i)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c kxypi = initial x index used
c kxypp = number of x indices used
c nyv = first dimension of g
c kzpp = number of z indices used
c kxypd = second dimension of g
c kzpd = third dimension of g
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c final fourier coefficients are stored as follows:
c h(1:3,l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js,
c kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
c 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
c h(1:3,l,1,k) = mode nx/2,kk-1,l-1,
c where ny/2+2 <= kk <= ny, 1 <= l <= nz,
c the following are located on node js = 0 and ks = 0:
c h(1:3,l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
c imag(h(1:3,1,1,1)) = real part of mode nx/2,0,0
c imag(h(1:3,nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
c the following are located on node js = 0 and ks = nyh/kyzp:
c h(1:3,l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
c imag(h(1:3,1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
c imag(h(1:3,nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvpy, nvpz, kxypi, kxypp
      integer nyv, kzpp, kxypd, kzpd, nxhyzd, nxyzhd
      integer mixup
      complex g, sct
      dimension g(3,nyv,kxypd,kzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, ny2, nz, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj, nryb
      integer js, ks, nry, kxypt, nn
      complex s, t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = max(1,ny/2)
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxypt = kxypi + kxypp - 1
c js/ks = processor co-ordinates in x/z => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      if (kstrt.gt.(nvpy*nvpz)) return
      if (isign.gt.0) go to 90
c inverse fourier transform
      nryb = nxhyz/ny
      nry = nxyz/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1,t2,t3)
      do 50 nn = 1, kxypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kxypi
c bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = g(1,k1,i,n)
         t2 = g(2,k1,i,n)
         t3 = g(3,k1,i,n)
         g(1,k1,i,n) = g(1,k,i,n)
         g(2,k1,i,n) = g(2,k,i,n)
         g(3,k1,i,n) = g(3,k,i,n)
         g(1,k,i,n) = t1
         g(2,k,i,n) = t2
         g(3,k,i,n) = t3
      endif
   10 continue
c then transform in y
      do 40 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t1 = s*g(1,j2,i,n)
      t2 = s*g(2,j2,i,n)
      t3 = s*g(3,j2,i,n)
      g(1,j2,i,n) = g(1,j1,i,n) - t1
      g(2,j2,i,n) = g(2,j1,i,n) - t2
      g(3,j2,i,n) = g(3,j1,i,n) - t3
      g(1,j1,i,n) = g(1,j1,i,n) + t1
      g(2,j1,i,n) = g(2,j1,i,n) + t2
      g(3,j1,i,n) = g(3,j1,i,n) + t3
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
!$OMP PARALLEL DO PRIVATE(k,n,jj,s)
         do 80 n = 1, kzpp
         do 70 k = 2, nyh
         do 60 jj = 1, 3
         s = g(jj,ny2-k,1,n)
         g(jj,ny2-k,1,n) = 0.5*cmplx(aimag(g(jj,k,1,n) + s),
     1                               real(g(jj,k,1,n) - s))
         g(jj,k,1,n) = 0.5*cmplx(real(g(jj,k,1,n) + s),
     1                           aimag(g(jj,k,1,n) - s))
   60    continue
   70    continue
   80    continue
!$OMP END PARALLEL DO
      endif
      return
c forward fourier transform
  90  nryb = nxhyz/ny
      nry = nxyz/ny
c scramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
!$OMP PARALLEL DO PRIVATE(k,n,jj,s)
         do 120 n = 1, kzpp
         do 110 k = 2, nyh
         do 100 jj = 1, 3
         s = cmplx(aimag(g(jj,ny2-k,1,n)),real(g(jj,ny2-k,1,n)))
         g(jj,ny2-k,1,n) = conjg(g(jj,k,1,n) - s)
         g(jj,k,1,n) = g(jj,k,1,n) + s
  100    continue
  110    continue
  120    continue
!$OMP END PARALLEL DO
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1,t2,t3)
      do 170 nn = 1, kxypp*kzpp
      i = (nn - 1)/kzpp
      n = nn - kzpp*i
      i = i + kxypi
c bit-reverse array elements in y
      do 130 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = g(1,k1,i,n)
         t2 = g(2,k1,i,n)
         t3 = g(3,k1,i,n)
         g(1,k1,i,n) = g(1,k,i,n)
         g(2,k1,i,n) = g(2,k,i,n)
         g(3,k1,i,n) = g(3,k,i,n)
         g(1,k,i,n) = t1
         g(2,k,i,n) = t2
         g(3,k,i,n) = t3
      endif
  130 continue
c then transform in y
      do 160 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 150 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 140 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*g(1,j2,i,n)
      t2 = s*g(2,j2,i,n)
      t3 = s*g(3,j2,i,n)
      g(1,j2,i,n) = g(1,j1,i,n) - t1
      g(2,j2,i,n) = g(2,j1,i,n) - t2
      g(3,j2,i,n) = g(3,j1,i,n) - t3
      g(1,j1,i,n) = g(1,j1,i,n) + t1
      g(2,j1,i,n) = g(2,j1,i,n) + t2
      g(3,j1,i,n) = g(3,j1,i,n) + t3
  140 continue
  150 continue
  160 continue
  170 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT32RM3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,   
     1nvpy,nvpz,kxypi,kxypp,nzv,kyzp,kxypd,kyzpd,nxhyzd,nxyzhd)
c this subroutine performs the z part of 3 three dimensional real to
c complex fast fourier transforms and their inverses for a subset of
c x and y, using complex arithmetic,
c for data which is distributed in blocks, with 2D spatial decomposition
c with OpenMP
c for isign = (-1,1), input: all, output: h
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c h(1:3,l,n,m) = sum(h(1:3,i,j,k)*exp(-sqrt(-1)*2pi*l*i/nz))
c if isign = 1, a forward fourier transform is performed
c h(1:3,l,n,m) = sum(h(1:3,i,j,k)*exp(sqrt(-1)*2pi*ll*ii/nz))
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c kxypi = initial x index used
c kxypp = number of x indices used
c nzv = first dimension of h
c kyzpp = number of y indices used
c kxypd = second dimension of h
c kyzpd = third dimension of h
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c h(1:3,l,j,k) = mode jj-1,kk-1,l, where jj = j + kxyp*js,
c kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
c 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
c h(1:3,l,1,k) = mode nx/2,kk-1,l-1,
c where ny/2+2 <= kk <= ny, 1 <= l <= nz,
c the following are located on node js = 0 and ks = 0:
c h(1:3,l,1,1) = mode nx/2,0,l-1, where 2 <= l <= nz/2
c imag(h(1:3,1,1,1)) = real part of mode nx/2,0,0
c imag(h(1:3,nz/2+1,1,1)) = real part of mode nx/2,0,nz/2
c the following are located on node js = 0 and ks = nyh/kyzp:
c h(1:3,l,1,ny/2+1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
c imag(h(1:3,1,1,ny/2+1)) = real part of mode nx/2,ny/2,0
c imag(h(1:3,nz/2+1,1,ny/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, nvpy, nvpz, kxypi, kxypp
      integer nzv, kyzp, kxypd, kyzpd, nxhyzd, nxyzhd
      integer mixup
      complex h, sct
      dimension h(3,nzv,kxypd,kyzpd)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, nz, nzh, nz2, nxyz, nxhyz
      integer j, k, l, i, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj, nrzb
      integer l1, js, ks, nrz, kxypt, kyzpp, kyzb, nn
      complex s, t1, t2, t3
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = max(1,ny/2)
      nz = 2**indz
      nzh = max(1,nz/2)
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxypt = kxypi + kxypp - 1
c js/ks = processor co-ordinates in x/y => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      kyzpp = min(kyzp,max(0,ny-kyzp*ks))
      if (kstrt.gt.(nvpy*nvpz)) return
      if (isign.gt.0) go to 100
c inverse fourier transform
      nrzb = nxhyz/nz
      nrz = nxyz/nz
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,l1,s,t1,t2,t3)
      do 50 nn = 1, kxypp*kyzpp
      i = (nn - 1)/kyzpp
      n = nn - kyzpp*i
      i = i +  kxypi
c bit-reverse array elements in z
      do 10 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         t1 = h(1,l1,i,n)
         t2 = h(2,l1,i,n)
         t3 = h(3,l1,i,n)
         h(1,l1,i,n) = h(1,l,i,n)
         h(2,l1,i,n) = h(2,l,i,n)
         h(3,l1,i,n) = h(3,l,i,n)
         h(1,l,i,n) = t1
         h(2,l,i,n) = t2
         h(3,l,i,n) = t3
      endif
   10 continue
c finally transform in z
      do 40 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t1 = s*h(1,j2,i,n)
      t2 = s*h(2,j2,i,n)
      t3 = s*h(3,j2,i,n)
      h(1,j2,i,n) = h(1,j1,i,n) - t1
      h(2,j2,i,n) = h(2,j1,i,n) - t2
      h(3,j2,i,n) = h(3,j1,i,n) - t3
      h(1,j1,i,n) = h(1,j1,i,n) + t1
      h(2,j1,i,n) = h(2,j1,i,n) + t2
      h(3,j1,i,n) = h(3,j1,i,n) + t3
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
         if (ks.eq.0) then
            do 70 n = 2, nzh
            do 60 jj = 1, 3
            s = h(jj,nz2-n,1,1)
            h(jj,nz2-n,1,1) = 0.5*cmplx(aimag(h(jj,n,1,1) + s),
     1                                  real(h(jj,n,1,1) - s))
            h(jj,n,1,1) = 0.5*cmplx(real(h(jj,n,1,1) + s),
     1                              aimag(h(jj,n,1,1) - s))
   60       continue
   70       continue
         endif
         kyzb = nyh/kyzp
         if (ks.eq.kyzb) then
            k1 = nyh - kyzb*kyzp + 1
            do 90 n = 2, nzh
            do 80 jj = 1, 3
            s = h(jj,nz2-n,1,k1)
            h(jj,nz2-n,1,k1) = 0.5*cmplx(aimag(h(jj,n,1,k1) + s),
     1                                   real(h(jj,n,1,k1) - s))

            h(jj,n,1,k1) = 0.5*cmplx(real(h(jj,n,1,k1) + s),
     1                               aimag(h(jj,n,1,k1) - s))
   80       continue
   90       continue
        endif
      endif
      return
c forward fourier transform
  100 nrzb = nxhyz/nz
      nrz = nxyz/nz
c scramble modes kx = 0, nx/2
      if ((js.eq.0).and.(kxypi.eq.1)) then
         if (ks.eq.0) then
            do 120 n = 2, nzh
            do 110 jj = 1, 3
            s = cmplx(aimag(h(jj,nz2-n,1,1)),real(h(jj,nz2-n,1,1)))
            h(jj,nz2-n,1,1) = conjg(h(jj,n,1,1) - s)
            h(jj,n,1,1) = h(jj,n,1,1) + s
  110       continue
  120       continue
         endif
         kyzb = nyh/kyzp
         if (ks.eq.kyzb) then
            k1 = nyh - kyzb*kyzp + 1
            do 140 n = 2, nzh
            do 130 jj = 1, 3
            s = cmplx(aimag(h(jj,nz2-n,1,k1)),real(h(jj,nz2-n,1,k1)))
            h(jj,nz2-n,1,k1) = conjg(h(jj,n,1,k1) - s)
            h(jj,n,1,k1) = h(jj,n,1,k1) + s
  130       continue
  140       continue
         endif
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,nn,ns,ns2,km,kmr,k1,k2,j1,j2,l1,s,t1,t2,t3)
      do 190 nn = 1, kxypp*kyzpp
      i = (nn - 1)/kyzpp
      n = nn - kyzpp*i
      i = i +  kxypi
c bit-reverse array elements in z
      do 150 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         t1 = h(1,l1,i,n)
         t2 = h(2,l1,i,n)
         t3 = h(3,l1,i,n)
         h(1,l1,i,n) = h(1,l,i,n)
         h(2,l1,i,n) = h(2,l,i,n)
         h(3,l1,i,n) = h(3,l,i,n)
         h(1,l,i,n) = t1
         h(2,l,i,n) = t2
         h(3,l,i,n) = t3
      endif
  150 continue
c first transform in z
      do 180 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*h(1,j2,i,n)
      t2 = s*h(2,j2,i,n)
      t3 = s*h(3,j2,i,n)
      h(1,j2,i,n) = h(1,j1,i,n) - t1
      h(2,j2,i,n) = h(2,j1,i,n) - t2
      h(3,j2,i,n) = h(3,j1,i,n) - t3
      h(1,j1,i,n) = h(1,j1,i,n) + t1
      h(2,j1,i,n) = h(2,j1,i,n) + t2
      h(3,j1,i,n) = h(3,j1,i,n) + t3
  160 continue
  170 continue
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      function ranorm()
c this program calculates a random number y from a gaussian distribution
c with zero mean and unit variance, according to the method of
c mueller and box:
c    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
c    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
c where x is a random number uniformly distributed on (0,1).
c written for the ibm by viktor k. decyk, ucla
      implicit none
      integer iflg,isc,i1,r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end
c-----------------------------------------------------------------------
      function randum()
c this is a version of the random number generator dprandom due to
c c. bingham and the yale computer center, producing numbers
c in the interval (0,1).  written for the sun by viktor k. decyk, ucla
      implicit none
      integer isc,i1,r1,r2
      double precision randum,h1l,h1u,r0,r3,asc,bsc
      save r1,r2,h1l,h1u
      data r1,r2 /1271199957,1013501921/
      data h1l,h1u /65533.0d0,32767.0d0/
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      randum = (dble(r1) + dble(r2)*asc)*asc
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPCOPYOUT(part,ppart,kpic,npp,npmax,nppmx,idimp,mxyzp1
     1,irc)
c for 3d code, this subroutine copies segmented particle data ppart to
c the array part with original tiled layout
c spatial decomposition in y/z direction
c input: all except part, npp, output: part, npp
c part(i,j) = i-th coordinate for particle j in partition
c ppart(i,j,k) = i-th coordinate for particle j in partition in tile k
c kpic = number of particles per tile
c npp = number of particles in partition
c npmax = maximum number of particles in each partition
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 6
c mxyzp1 = total number of tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer npp, npmax, nppmx, idimp, mxyzp1, irc
      real part, ppart
      integer kpic
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyzp1)
      dimension kpic(mxyzp1)
c local data
      integer i, j, k, npoff, nppp, ne, ierr
      npoff = 0
      ierr = 0
c loop over tiles
      do 30 k = 1, mxyzp1
      nppp = kpic(k)
      ne = nppp + npoff
      if (ne.gt.npmax) ierr = max(ierr,ne-npmax)
      if (ierr.gt.0) nppp = 0
c loop over particles in tile
      do 20 j = 1, nppp
      do 10 i = 1, idimp
      part(i,j+npoff) = ppart(i,j,k)
   10 continue
   20 continue
      npoff = npoff + nppp
   30 continue
      npp = npoff
      if (ierr.gt.0) irc = ierr
      return
      end
