c Fortran Library for Skeleton 2-1/2D Electromagnetic MPI/OpenMP
c PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine PDICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,idps)
c this subroutine determines spatial boundaries for uniform particle
c decomposition, calculates number of grid points in each spatial
c region, and the offset of these grid points from the global address
c nvp must be < ny.  some combinations of ny and nvp result in a zero
c value of nyp.  this is not supported.
c integer boundaries are set.
c input: ny, kstrt, nvp, idps, output: edges, nyp, noff, nypmx, nypmn
c edges(1) = lower boundary of particle partition
c edges(2) = upper boundary of particle partition
c nyp = number of primary (complete) gridpoints in particle partition
c noff = lowermost global gridpoint in particle partition
c nypmx = maximum size of particle partition, including guard cells
c nypmn = minimum value of nyp
c ny = system length in y direction
c kstrt = starting data block number (processor id + 1)
c nvp = number of real or virtual processors
c idps = number of partition boundaries
      implicit none
      integer nyp, noff, nypmx, nypmn, ny, kstrt, nvp, idps
      real edges
      dimension edges(idps)
c local data
      integer kb, kyp
      real at1, any
      integer mypm, iwork2
      dimension mypm(2), iwork2(2)
      any = real(ny)
c determine decomposition
      kb = kstrt - 1
      kyp = (ny - 1)/nvp + 1
      at1 = real(kyp)
      edges(1) = at1*real(kb)
      if (edges(1).gt.any) edges(1) = any
      noff = edges(1)
      edges(2) = at1*real(kb + 1)
      if (edges(2).gt.any) edges(2) = any
      kb = edges(2)
      nyp = kb - noff
c find maximum/minimum partition size
      mypm(1) = nyp
      mypm(2) = -nyp
      call PPIMAX(mypm,iwork2,2)
      nypmx = mypm(1) + 1
      nypmn = -mypm(2)
      return
      end
c-----------------------------------------------------------------------
      subroutine PDISTR2H(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx
     1,npy,nx,ny,idimp,npmax,idps,ipbc,ierr)
c for 2-1/2d code, this subroutine calculates initial particle
c co-ordinates and velocities with uniform density and maxwellian
c velocity with drift for distributed data.
c input: all except part, ierr, output: part, npp, ierr
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = velocity vx of particle n in partition
c part(4,n) = velocity vy of particle n in partition
c part(5,n) = velocity vz of particle n in partition
c edges(1) = lower boundary of particle partition
c edges(2) = upper boundary of particle partition
c npp = number of particles in partition
c nps = starting address of particles in partition
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy = initial number of particles distributed in x/y direction
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c idps = number of partition boundaries
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ierr = (0,1) = (no,yes) error condition exists
c ranorm = gaussian random number with zero mean and unit variance
c with spatial decomposition
      implicit none
      integer npp, nps, npx, npy, nx, ny, idimp, npmax, idps, ipbc, ierr
      real vtx, vty, vtz, vdx, vdy, vdz
      real part, edges
      dimension part(idimp,npmax), edges(idps)
c local data
      integer j, k, npt, npxyp
      real edgelx, edgely, at1, at2, xt, yt, vxt, vyt, vzt
      double precision dnpx, dnpxy, dt1
      integer ierr1, iwork1
      double precision sum4, work4
      dimension ierr1(1), iwork1(1), sum4(4), work4(4)
      double precision ranorm
      ierr = 0
c particle distribution constant
      dnpx = dble(npx)
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
c uniform density profile
      do 20 k = 1, npy
      yt = edgely + at2*(real(k) - 0.5)
      do 10 j = 1, npx
      xt = edgelx + at1*(real(j) - 0.5)
c maxwellian velocity distribution
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      vzt = vtz*ranorm()
      if ((yt.ge.edges(1)).and.(yt.lt.edges(2))) then
         npt = npp + 1
         if (npt.le.npmax) then
            part(1,npt) = xt
            part(2,npt) = yt
            part(3,npt) = vxt
            part(4,npt) = vyt
            part(5,npt) = vzt
            npp = npt
         else
            ierr = ierr + 1
         endif
      endif
   10 continue
   20 continue
      npxyp = 0
c add correct drift
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 30 j = nps, npp
      npxyp = npxyp + 1
      sum4(1) = sum4(1) + part(3,j)
      sum4(2) = sum4(2) + part(4,j)
      sum4(3) = sum4(3) + part(5,j)
   30 continue
      sum4(4) = npxyp
      call PPDSUM(sum4,work4,4)
      dnpxy = sum4(4)
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      dt1 = 1.0d0/dnpxy
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 40 j = nps, npp
      part(3,j) = part(3,j) - sum4(1)
      part(4,j) = part(4,j) - sum4(2)
      part(5,j) = part(5,j) - sum4(3)
   40 continue
c process errors
      dnpxy = dnpxy - dnpx*dble(npy)
      if (dnpxy.ne.0.0d0) ierr = dnpxy
      return
      end
c-----------------------------------------------------------------------
      subroutine PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,  
     1mx1,mxyp1,irc)
c this subroutine finds the maximum number of particles in each tile of
c mx, my to calculate size of segmented particle array ppart
c linear interpolation, spatial decomposition in y direction
c input: all except kpic, nppmx, output: kpic, nppmx
c part = input particle array
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c kpic = output number of particles per tile
c nppmx = return maximum number of particles in tile
c npp = number of particles in partition
c noff = backmost global gridpoint in particle partition
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c mx/my = number of grids in sorting cell in x and y
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
      integer kpic, npp, noff
      real part
      dimension part(idimp,npmax)
      dimension kpic(mxyp1)
c local data
      integer j, k, n, m, mnoff, isum, ist, npx, ierr
      mnoff = noff
      ierr = 0
c clear counter array
      do 10 k = 1, mxyp1
      kpic(k) = 0
   10 continue
c find how many particles in each tile
      do 20 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      m = n + mx1*m
      if (m.le.mxyp1) then
         kpic(m) = kpic(m) + 1
      else
         ierr = max(ierr,m-mxyp1)
      endif
   20 continue
c find maximum
      isum = 0
      npx = 0
      do 30 k = 1, mxyp1
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
      subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,npmax, 
     1mx,my,mx1,mxyp1,irc)
c this subroutine sorts particles by x,y grid in tiles of
c mx, my and copies to segmented array ppart
c linear interpolation, spatial decomposition in y direction
c input: all except ppart, kpic, output: ppart, kpic
c part/ppart = input/output particle arrays
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c kpic = output number of particles per tile
c nppmx = maximum number of particles in tile
c npp = number of particles in partition
c noff = backmost global gridpoint in particle partition
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c mx/my = number of grids in sorting cell in x and y
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
      integer kpic, npp, noff
      real part, ppart
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1)
c local data
      integer i, j, k, n, m, mnoff, ip, ierr
      mnoff = noff
      ierr = 0
c clear counter array
      do 10 k = 1, mxyp1
      kpic(k) = 0
   10 continue
c find addresses of particles at each tile and reorder particles
      do 30 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      m = n + mx1*m
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
      subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,mx1
     1,myp1,irc)
c this subroutine performs a sanity check to make sure particles sorted
c by x,y grid in tiles of mx, my, are all within bounds.
c tiles are assumed to be arranged in 2D linear memory
c input: all except irc
c output: irc
c ppart(1,n,k) = position x of particle n in tile k
c ppart(2,n,k) = position y of particle n in tile k
c kpic(k) = number of reordered output particles in tile k
c noff = lowermost global gridpoint in particle partition.
c nyp = number of primary (complete) gridpoints in particle partition
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx = system length in x direction
c mx/my = number of grids in sorting cell in x/y
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, mx1, myp1, irc
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension kpic(mx1*myp1)
c local data
      integer mxyp1, noffp, moffp, nppp, j, k, ist, nn, mm
      real edgelx, edgely, edgerx, edgery, dx, dy
      mxyp1 = mx1*myp1
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noffp,moffp,nppp,nn,mm,ist,edgelx,edgely,edgerx,
!$OMP& edgery,dx,dy)
      do 20 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
c loop over particles in tile
      do 10 j = 1, nppp
      dx = ppart(1,j,k)
      dy = ppart(2,j,k)
c find particles going out of bounds
      ist = 0
      if (dx.lt.edgelx) ist = 1
      if (dx.ge.edgerx) ist = 2
      if (dy.lt.edgely) ist = ist + 3
      if (dy.ge.edgery) ist = ist + 6
      if (ist.gt.0) irc = k
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ek,
     1idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c OpenMP version using guard cells, for distributed data
c data read in tiles
c particles stored segmented array
c 119 flops/particle, 1 divide, 29 loads, 5 stores
c input: all, output: ppart, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c ppart(3,n,m) = velocity vx of particle n in partition in tile m
c ppart(4,n,m) = velocity vy of particle n in partition in tile m
c ppart(5,n,m) = velocity vz of particle n in partition in tile m
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c fxy(3,j,k) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape,
c where kk = k + noff - 1
c bxy(1,j,k) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape,
c where kk = k + noff - 1
c kpic = number of particles per tile
c noff = lowermost global gridpoint in particle partition.
c nyp = number of primary (complete) gridpoints in particle partition
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ipbc
      real qbm, dt, dtc, ek
      real ppart, fxy, bxy
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension kpic(mxyp1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y
      real sfxy, sbxy
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
c     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      sum2 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,dxp,dyp,amx,amy,  
!$OMP& dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,
!$OMP& rot3,rot4,rot5,rot6,rot7,rot8,rot9,sum1,sfxy,sbxy)
!$OMP& REDUCTION(+:sum2)
      do 60 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
c load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 50 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
      ppart(5,j,k) = dz
c new position
      dx = x + dx*dtc
      dy = y + dy*dtc
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
      endif
c set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   50 continue
      sum2 = sum2 + sum1
   60 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm
     1,dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells, for distributed data
c data read in tiles
c particles stored segmented array
c 119 flops/particle, 1 divide, 29 loads, 5 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c ppart(3,n,m) = velocity vx of particle n in partition in tile m
c ppart(4,n,m) = velocity vy of particle n in partition in tile m
c ppart(5,n,m) = velocity vz of particle n in partition in tile m
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c fxy(3,j,k) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape,
c where kk = k + noff - 1
c bxy(1,j,k) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape,
c where kk = k + noff - 1
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c noff = lowermost global gridpoint in particle partition.
c nyp = number of primary (complete) gridpoints in particle partition
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ntmax, irc
      real qbm, dt, dtc, ek
      real ppart, fxy, bxy
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, ih, nh, nn, mm
      real qtmh, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
      real x, y
      real sfxy, sbxy
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
c     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx,
!$OMP& amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1, 
!$OMP& rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,
!$OMP& edgery,sum1,sfxy,sbxy)
!$OMP& REDUCTION(+:sum2)
      do 70 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
      ih = 0
      nh = 0
      mnoff = moffp + noff - 1
c load local fields from global arrays
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm+1
      do 30 i = 1, nn+1
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
c clear counters
      do 50 j = 1, 8
      ncl(j,k) = 0
   50 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 60 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
      ppart(5,j,k) = dz
c new position
      dx = x + dx*dtc
      dy = y + dy*dtc
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
c set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
c increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
   60 continue
      sum2 = sum2 + sum1
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   70 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci
     1,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c OpenMP version using guard cells, for distributed data
c data read in tiles
c particles stored segmented array
c 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: ppart, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c ppart(3,n,m) = momentum vx of particle n in partition in tile m
c ppart(4,n,m) = momentum vy of particle n in partition in tile m
c ppart(5,n,m) = momentum vz of particle n in partition in tile m
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c fxy(3,j,k) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape,
c where kk = k + noff - 1
c bxy(1,j,k) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape,
c where kk = k + noff - 1
c kpic = number of particles per tile
c noff = lowermost global gridpoint in particle partition.
c nyp = number of primary (complete) gridpoints in particle partition
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprical of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 5
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ipbc
      real qbm, dt, dtc, ci, ek
      real ppart, fxy, bxy
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension kpic(mxyp1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y
      real sfxy, sbxy
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
c     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum2 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,dxp,dyp,amx,amy,
!$OMP& dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,
!$OMP& rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,dtg,sum1,sfxy,
!$OMP& sbxy)
!$OMP& REDUCTION(+:sum2)
      do 60 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
c load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 50 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new momentum
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
      ppart(5,j,k) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = x + dx*dtg
      dy = y + dy*dtg
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
      endif
c set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   50 continue
      sum2 = sum2 + sum1
   60 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,  
     1qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax
     2,irc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells, for distributed data
c data read in tiles
c particles stored segmented array
c 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c ppart(3,n,m) = momentum vx of particle n in partition in tile m
c ppart(4,n,m) = momentum vy of particle n in partition in tile m
c ppart(5,n,m) = momentum vz of particle n in partition in tile m
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c fxy(3,j,k) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape,
c where kk = k + noff - 1
c bxy(1,j,k) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape,
c where kk = k + noff - 1
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c noff = lowermost global gridpoint in particle partition.
c nyp = number of primary (complete) gridpoints in particle partition
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprical of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 5
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ntmax, irc
      real qbm, dt, dtc, ci, ek
      real ppart, fxy, bxy
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, ih, nh, nn, mm
      real qtmh, ci2, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
      real x, y
      real sfxy, sbxy
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
c     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx,
!$OMP& amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,
!$OMP& rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,
!$OMP& edgery,p2,gami,qtmg,dtg,sum1,sfxy,sbxy)
!$OMP& REDUCTION(+:sum2)
      do 70 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
      ih = 0
      nh = 0
      mnoff = moffp + noff - 1
c load local fields from global arrays
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm+1
      do 30 i = 1, nn+1
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
c clear counters
      do 50 j = 1, 8
      ncl(j,k) = 0
   50 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 60 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new momentum
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
      ppart(5,j,k) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = x + dx*dtg
      dy = y + dy*dtg
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
c set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
c increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
   60 continue
      sum2 = sum2 + sum1
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   70 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv, 
     1nypmx,mx1,mxyp1)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c OpenMP version using guard cells, for distributed data
c data deposited in tiles
c particles stored segmented array
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c q(j,k) = charge density at grid point (j,kk),
c where kk = k + noff - 1
c kpic = number of particles per tile
c noff = lowermost global gridpoint in particle partition.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, idimp, nppmx, mx, my, nxv, nypmx, mx1, mxyp1
      real qm
      real ppart, q
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), q(nxv,nypmx), kpic(mxyp1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real x, y, dxp, dyp, amx, amy
      real sq
c     dimension sq(MXV,MYV)
      dimension sq(mx+1,my+1)
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,mnoff,nn,mm,x,y,dxp,dyp,amx,amy,
!$OMP& sq)
      do 80 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
c zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      sq(i,j) = 0.0
   10 continue
   20 continue
c loop over particles in tile
      do 30 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit charge within tile to local accumulator
      x = sq(nn,mm) + amx*amy
      y = sq(nn+1,mm) + dxp*amy
      sq(nn,mm) = x
      sq(nn+1,mm) = y
      x = sq(nn,mm+1) + amx*dyp
      y = sq(nn+1,mm+1) + dxp*dyp
      sq(nn,mm+1) = x
      sq(nn+1,mm+1) = y
   30 continue
c deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 50 j = 2, mm
      do 40 i = 2, nn
      q(i+noffp,j+moffp) = q(i+noffp,j+moffp) + sq(i,j)
   40 continue
   50 continue
c deposit charge to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 60 i = 2, nn
!$OMP ATOMIC
      q(i+noffp,1+moffp) = q(i+noffp,1+moffp) + sq(i,1)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noffp,mm+moffp) = q(i+noffp,mm+moffp) + sq(i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noffp)
      do 70 j = 1, mm
!$OMP ATOMIC
      q(1+noffp,j+moffp) = q(1+noffp,j+moffp) + sq(1,j)
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noffp,j+moffp) = q(nn+noffp,j+moffp) + sq(nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGJPPOST2L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,ny,
     1mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c OpenMP version using guard cells, for distributed data
c data deposited in tiles
c particles stored segmented array
c 41 flops/particle, 17 loads, 14 stores
c input: all, output: ppart, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c ppart(3,n,m) = x velocity of particle n in partition in tile m
c ppart(4,n,m) = y velocity of particle n in partition in tile m
c ppart(5,n,m) = z velocity of particle n in partition in tile m
c cu(i,j,k) = ith component of current density at grid point (j,kk),
c where kk = k + noff - 1
c kpic = number of particles per tile
c noff = lowermost global gridpoint in particle partition.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1, mxyp1
      integer ipbc
      real qm, dt
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv,nypmx)
      dimension kpic(mxyp1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real scu
      dimension scu(3,MXV,MYV)
c     dimension scu(3,mx+1,my+1)
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,dxp,dyp,amx,amy,dx
!$OMP& ,dy,vx,vy,vz,scu)
      do 80 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
c zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   10 continue
   20 continue
c loop over particles in tile
      do 30 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      scu(1,nn,mm) = scu(1,nn,mm) + vx*dx
      scu(2,nn,mm) = scu(2,nn,mm) + vy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + vz*dx
      dx = amx*dyp
      scu(1,nn+1,mm) = scu(1,nn+1,mm) + vx*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + vy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + vz*dy
      dy = dxp*dyp
      scu(1,nn,mm+1) = scu(1,nn,mm+1) + vx*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + vy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + vz*dx
      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + vx*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + vy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
      endif
c set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   30 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 50 j = 2, mm
      do 40 i = 2, nn
      cu(1,i+noffp,j+moffp) = cu(1,i+noffp,j+moffp) + scu(1,i,j)
      cu(2,i+noffp,j+moffp) = cu(2,i+noffp,j+moffp) + scu(2,i,j)
      cu(3,i+noffp,j+moffp) = cu(3,i+noffp,j+moffp) + scu(3,i,j)
   40 continue
   50 continue
c deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 60 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noffp)
      do 70 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,   
     1nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells, for distributed data
c data deposited in tiles
c particles stored segmented array
c 41 flops/particle, 17 loads, 14 stores
c input: all except ncl, ihole, irc,
c output: ppart, cu, ncl, ihole, irc
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c ppart(3,n,m) = x velocity of particle n in partition in tile m
c ppart(4,n,m) = y velocity of particle n in partition in tile m
c ppart(5,n,m) = z velocity of particle n in partition in tile m
c cu(i,j,k) = ith component of current density at grid point (j,kk),
c where kk = k + noff - 1
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c noff = lowermost global gridpoint in particle partition.
c nyp = number of primary (complete) gridpoints in particle partition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer noff, nyp, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1
      integer mxyp1, ntmax, irc
      real qm, dt
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv,nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, ih, nh, nn, mm
      real dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real anx, any, edgelx, edgely, edgerx, edgery
      real scu
      dimension scu(3,MXV,MYV)
c     dimension scu(3,mx+1,my+1)
      anx = real(nx)
      any = real(ny)
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx,
!$OMP& amy,dx,dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,scu)
      do 90 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
      ih = 0
      nh = 0
      mnoff = moffp + noff - 1
c zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   10 continue
   20 continue
c clear counters
      do 30 j = 1, 8
      ncl(j,k) = 0
   30 continue
c loop over particles in tile
      do 40 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      scu(1,nn,mm) = scu(1,nn,mm) + vx*dx
      scu(2,nn,mm) = scu(2,nn,mm) + vy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + vz*dx
      dx = amx*dyp
      scu(1,nn+1,mm) = scu(1,nn+1,mm) + vx*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + vy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + vz*dy
      dy = dxp*dyp
      scu(1,nn,mm+1) = scu(1,nn,mm+1) + vx*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + vy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + vz*dx
      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + vx*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + vy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
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
c set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
c increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
   40 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 60 j = 2, mm
      do 50 i = 2, nn
      cu(1,i+noffp,j+moffp) = cu(1,i+noffp,j+moffp) + scu(1,i,j)
      cu(2,i+noffp,j+moffp) = cu(2,i+noffp,j+moffp) + scu(2,i,j)
      cu(3,i+noffp,j+moffp) = cu(3,i+noffp,j+moffp) + scu(3,i,j)
   50 continue
   60 continue
c deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 70 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
   70 continue
      nn = min(mx+1,nxv-noffp)
      do 80 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
   80 continue
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   90 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGRJPPOST2L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp,nx
     1,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c OpenMP version using guard cells, for distributed data
c data deposited in tiles
c particles stored segmented array
c 47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
c input: all, output: ppart, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c ppart(3,n,m) = x momentum of particle n in partition in tile m
c ppart(4,n,m) = y momentum of particle n in partition in tile m
c ppart(5,n,m) = z momentum of particle n in partition in tile m
c cu(i,j,k) = ith component of current density at grid point (j,kk),
c where kk = k + noff - 1
c kpic = number of particles per tile
c noff = lowermost global gridpoint in particle partition.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprical of velocity of light
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1, mxyp1
      integer ipbc
      real qm, dt, ci
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv,nypmx)
      dimension kpic(mxyp1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, p2, gami
      real scu
      dimension scu(3,MXV,MYV)
c     dimension scu(3,mx+1,my+1)
      ci2 = ci*ci
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,dxp,dyp,amx,amy,dx
!$OMP& ,dy,vx,vy,vz,p2,gami,scu)
      do 80 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
c zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   10 continue
   20 continue
c loop over particles in tile
      do 30 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
c find inverse gamma
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      scu(1,nn,mm) = scu(1,nn,mm) + vx*dx
      scu(2,nn,mm) = scu(2,nn,mm) + vy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + vz*dx
      dx = amx*dyp
      scu(1,nn+1,mm) = scu(1,nn+1,mm) + vx*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + vy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + vz*dy
      dy = dxp*dyp
      scu(1,nn,mm+1) = scu(1,nn,mm+1) + vx*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + vy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + vz*dx
      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + vx*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + vy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
      endif
c set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   30 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 50 j = 2, mm
      do 40 i = 2, nn
      cu(1,i+noffp,j+moffp) = cu(1,i+noffp,j+moffp) + scu(1,i,j)
      cu(2,i+noffp,j+moffp) = cu(2,i+noffp,j+moffp) + scu(2,i,j)
      cu(3,i+noffp,j+moffp) = cu(3,i+noffp,j+moffp) + scu(3,i,j)
   40 continue
   50 continue
c deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 60 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noffp)
      do 70 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGRJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci
     1,nppmx,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells, for distributed data
c data deposited in tiles
c particles stored segmented array
c 47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
c input: all except ncl, ihole, irc,
c output: ppart, cu, ncl, ihole, irc
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c ppart(1,n,m) = position x of particle n in partition in tile m
c ppart(2,n,m) = position y of particle n in partition in tile m
c ppart(3,n,m) = x momentum of particle n in partition in tile m
c ppart(4,n,m) = y momentum of particle n in partition in tile m
c ppart(5,n,m) = z momentum of particle n in partition in tile m
c cu(i,j,k) = ith component of current density at grid point (j,kk),
c where kk = k + noff - 1
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c noff = lowermost global gridpoint in particle partition.
c nyp = number of primary (complete) gridpoints in particle partition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprical of velocity of light
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer noff, nyp, nppmx, idimp, nx, ny, mx, my, nxv, nypmx, mx1
      integer mxyp1, ntmax, irc
      real qm, dt, ci
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1), cu(3,nxv,nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm, ih, nh
      real ci2, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, p2, gami
      real anx, any, edgelx, edgely, edgerx, edgery
      real scu
      dimension scu(3,MXV,MYV)
c     dimension scu(3,mx+1,my+1)
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx,
!$OMP& amy,dx,dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,p2,gami,scu)
      do 90 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
      ih = 0
      nh = 0
      mnoff = moffp + noff - 1
c zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   10 continue
   20 continue
c clear counters
      do 30 j = 1, 8
      ncl(j,k) = 0
   30 continue
c loop over particles in tile
      do 40 j = 1, nppp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
c find inverse gamma
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      scu(1,nn,mm) = scu(1,nn,mm) + vx*dx
      scu(2,nn,mm) = scu(2,nn,mm) + vy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + vz*dx
      dx = amx*dyp
      scu(1,nn+1,mm) = scu(1,nn+1,mm) + vx*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + vy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + vz*dy
      dy = dxp*dyp
      scu(1,nn,mm+1) = scu(1,nn,mm+1) + vx*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + vy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + vz*dx
      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + vx*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + vy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
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
c set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
c increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
   40 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 60 j = 2, mm
      do 50 i = 2, nn
      cu(1,i+noffp,j+moffp) = cu(1,i+noffp,j+moffp) + scu(1,i,j)
      cu(2,i+noffp,j+moffp) = cu(2,i+noffp,j+moffp) + scu(2,i,j)
      cu(3,i+noffp,j+moffp) = cu(3,i+noffp,j+moffp) + scu(3,i,j)
   50 continue
   60 continue
c deposit current to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 70 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
   70 continue
      nn = min(mx+1,nxv-noffp)
      do 80 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
   80 continue
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   90 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,   
     1ncll,nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,  
     2nbmax,irc)
c this subroutine performs first part of a particle sort by x,y grid
c in tiles of mx, my
c linear interpolation, with periodic boundary conditions
c for distributed data, with 1d domain decomposition in y.
c tiles are assumed to be arranged in 2D linear memory
c this part of the algorithm has 3 steps.  first, one finds particles
c leaving tile and stores their number in each directon, location, and
c destination in ncl and ihole.  then, a prefix scan of ncl is performed
c and departing particles are buffered in ppbuff in direction order.
c finally, we buffer particles leaving the processor in sbufl and sbufr,
c and store particle number offsets in ncll and nclr.
c input: all except ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
c output: ppart, ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
c ppart(1,n,k) = position x of particle n in tile k
c ppart(2,n,k) = position y of particle n in tile k
c ppbuff(i,n,k) = i co-ordinate of particle n in tile k
c sbufl = buffer for particles being sent to lower processor
c sbufr = buffer for particles being sent to upper processor
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = direction destination of particle leaving hole
c all for tile k
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c ncll = number offset being sent to lower processor
c nclr = number offset being sent to upper processor
c noff = lowermost global gridpoint in particle partition.
c nyp = number of primary (complete) gridpoints in particle partition
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c nbmax =  size of buffers for passing particles between processors
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, mx1, myp1, npbmx
      integer ntmax, nbmax, irc
      real ppart, ppbuff, sbufl, sbufr
      integer kpic, ncl, ihole, ncll, nclr
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension kpic(mx1*myp1), ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension ncll(3,mx1), nclr(3,mx1)
c local data
      integer mxyp1, noffp, moffp, nppp
      integer i, j, k, ii, ih, nh, ist, nn, mm, isum, ip, j1, kk
      real anx, any, edgelx, edgely, edgerx, edgery, dx, dy
      mxyp1 = mx1*myp1
      anx = real(nx)
      any = real(ny)
c find and count particles leaving tiles and determine destination
c update ppart, ihole, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noffp,moffp,nppp,nn,mm,ih,nh,ist,dx,dy,edgelx,edgely,
!$OMP& edgerx,edgery)
      do 30 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      ih = 0
      nh = 0
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
c clear counters
      do 10 j = 1, 8
      ncl(j,k) = 0
   10 continue
c loop over particles in tile
      do 20 j = 1, nppp
      dx = ppart(1,j,k)
      dy = ppart(2,j,k)
c find particles going out of bounds
      ist = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(1,j,k) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(1,j,k) = dx
         else
            ist = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(2,j,k) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(2,j,k) = dy
         else
            ist = ist + 3
         endif
      endif
      if (ist.gt.0) then
         ncl(ist,k) = ncl(ist,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = ist
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
      ihole(1,1,k) = ih
   30 continue
!$OMP END PARALLEL DO
c ihole overflow
      if (irc.gt.0) return
c
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 70 k = 1, mxyp1
c find address offset for ordered ppbuff array
      isum = 0
      do 40 j = 1, 8
      ist = ncl(j,k)
      ncl(j,k) = isum
      isum = isum + ist
   40 continue
      nh = ihole(1,1,k)
      ip = 0
c loop over particles leaving tile
      do 60 j = 1, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = ncl(ist,k) + 1
      if (ii.le.npbmx) then
         do 50 i = 1, idimp
         ppbuff(i,ii,k) = ppart(i,j1,k)
   50    continue
      else
         ip = 1
      endif
      ncl(ist,k) = ii
   60 continue
c set error
      if (ip.gt.0) irc = ncl(8,k)
   70 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c buffer particles and their number leaving the node:
c update sbufl, sbufr, ncll, nclr
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k)
      do 80 k = 1, mx1
      ncll(1,k) = ncl(5,k) - ncl(2,k)
      nclr(1,k) = ncl(8,k+kk) - ncl(5,k+kk)
   80 continue
!$OMP END PARALLEL DO
c perform prefix scan
      kk = 1
   90 if (kk.ge.mx1) go to 110
!$OMP PARALLEL DO PRIVATE(k,ii,nn,mm)
      do 100 k = 1, mx1
      ii = (k - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + k + kk
      if (nn.le.mx1) then
         ncll(1,nn) = ncll(1,nn) + ncll(1,mm+1)
         nclr(1,nn) = nclr(1,nn) + nclr(1,mm+1)
      endif
  100 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 90
  110 kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(i,j,k,ii,nn,mm)
      do 180 k = 1, mx1
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,k) - ii
      do 130 j = 1, min(ii,nbmax-nn)
      do 120 i = 1, idimp
      sbufl(i,j+nn) = ppbuff(i,j+ncl(2,k),k)
  120 continue
  130 continue
      do 140 i = 1, 3
      ncll(i,k) = ncl(i+2,k) - ncl(2,k) + nn
  140 continue
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,k) - ii
      do 160 j = 1, min(ii,nbmax-mm)
      do 150 i = 1, idimp
      sbufr(i,j+mm) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  150 continue
  160 continue
      do 170 i = 1, 3
      nclr(i,k) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  170 continue
  180 continue
!$OMP END PARALLEL DO
c sbufl or sbufr overflow
      ii = max(ncll(3,mx1),nclr(3,mx1))
      if (ii.gt.nbmax) then
         irc = ii
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,  
     1nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
c this subroutine performs first part of a particle sort by x,y grid
c in tiles of mx, my
c linear interpolation, with periodic boundary conditions
c for distributed data, with 1d domain decomposition in y.
c tiles are assumed to be arranged in 2D linear memory
c this part of the algorithm has 2 steps.  first, a prefix scan of ncl
c is performed and departing particles are buffered in ppbuff in
c direction order. then, we buffer particles leaving the processor in
c sbufl and sbufr, and store particle number offsets in ncll and nclr.
c it assumes that the number, location, and destination of particles 
c leaving a tile have been previously stored in ncl and ihole by the
c PPGPPUSHF2L subroutine.
c input: all except ppbuff, sbufl, sbufr, ncll, nclr, irc
c output: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
c ppart(1,n,k) = position x of particle n in tile k
c ppart(2,n,k) = position y of particle n in tile k
c ppbuff(i,n,k) = i co-ordinate of particle n in tile k
c sbufl = buffer for particles being sent to lower processor
c sbufr = buffer for particles being sent to upper processor
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = direction destination of particle leaving hole
c all for tile k
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c ncll = number offset being sent to lower processor
c nclr = number offset being sent to upper processor
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c nbmax =  size of buffers for passing particles between processors
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, myp1, npbmx, ntmax, nbmax, irc
      real ppart, ppbuff, sbufl, sbufr
      integer ncl, ihole, ncll, nclr
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension ncll(3,mx1), nclr(3,mx1)
c local data
      integer mxyp1
      integer i, j, k, ii, nh, ist, nn, mm, isum, ip, j1, kk
      mxyp1 = mx1*myp1
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 40 k = 1, mxyp1
c find address offset for ordered ppbuff array
      isum = 0
      do 10 j = 1, 8
      ist = ncl(j,k)
      ncl(j,k) = isum
      isum = isum + ist
   10 continue
      nh = ihole(1,1,k)
      ip = 0
c loop over particles leaving tile
      do 30 j = 1, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = ncl(ist,k) + 1
      if (ii.le.npbmx) then
         do 20 i = 1, idimp
         ppbuff(i,ii,k) = ppart(i,j1,k)
   20    continue
      else
         ip = 1
      endif
      ncl(ist,k) = ii
   30 continue
c set error
      if (ip.gt.0) irc = ncl(8,k)
   40 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c buffer particles and their number leaving the node:
c update sbufl, sbufr, ncll, nclr
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k)
      do 50 k = 1, mx1
      ncll(1,k) = ncl(5,k) - ncl(2,k)
      nclr(1,k) = ncl(8,k+kk) - ncl(5,k+kk)
   50 continue
!$OMP END PARALLEL DO
c perform prefix scan
      kk = 1
   60 if (kk.ge.mx1) go to 80
!$OMP PARALLEL DO PRIVATE(k,ii,nn,mm)
      do 70 k = 1, mx1
      ii = (k - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + k + kk
      if (nn.le.mx1) then
         ncll(1,nn) = ncll(1,nn) + ncll(1,mm+1)
         nclr(1,nn) = nclr(1,nn) + nclr(1,mm+1)
      endif
   70 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 60
   80 kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(i,j,k,ii,nn,mm)
      do 150 k = 1, mx1
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,k) - ii
      do 100 j = 1, min(ii,nbmax-nn)
      do 90 i = 1, idimp
      sbufl(i,j+nn) = ppbuff(i,j+ncl(2,k),k)
   90 continue
  100 continue
      do 110 i = 1, 3
      ncll(i,k) = ncl(i+2,k) - ncl(2,k) + nn
  110 continue
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,k) - ii
      do 130 j = 1, min(ii,nbmax-mm)
      do 120 i = 1, idimp
      sbufr(i,j+mm) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  120 continue
  130 continue
      do 140 i = 1, 3
      nclr(i,k) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  140 continue
  150 continue
!$OMP END PARALLEL DO
c sbufl or sbufr overflow
      ii = max(ncll(3,mx1),nclr(3,mx1))
      if (ii.gt.nbmax) then
         irc = ii
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,   
     1mcll,mclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
c this subroutine performs second part of a particle sort by x,y grid
c in tiles of mx, my
c linear interpolation, with periodic boundary conditions
c for distributed data, with 1d domain decomposition in y.
c tiles are assumed to be arranged in 2D linear memory
c incoming particles from other tiles are copied from ppbuff, rbufl, and
c rbufr into ppart
c input: all except ppart, kpic, irc
c output: ppart, kpic, irc
c ppart(1,n,k) = position x of particle n in tile k
c ppart(2,n,k) = position y of particle n in tile k
c ppbuff(i,n,k) = i co-ordinate of particle n in tile k
c rbufl = buffer for particles being received from lower processor
c rbufr = buffer for particles being received from upper processor
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = direction destination of particle leaving hole
c all for tile k
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c mcll = number offset being received from lower processor
c mclr = number offset being received from upper processor
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c nbmax =  size of buffers for passing particles between processors
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, myp1, npbmx
      integer ntmax, nbmax, irc
      real ppart, ppbuff, rbufl, rbufr
      integer kpic, ncl, ihole, mcll, mclr
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension rbufl(idimp,nbmax), rbufr(idimp,nbmax)
      dimension kpic(mx1*myp1), ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension mcll(3,mx1), mclr(3,mx1)
c local data
      integer mxyp1, nppp, ncoff, noff, moff
      integer i, j, k, ii, kx, ky, ih, nh, ist
      integer ip, j1, j2, kxl, kxr, kk, kl, kr
      integer ks
      dimension ks(8)
      mxyp1 = mx1*myp1
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,kk,nppp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,noff,
!$OMP& moff,ist,j1,j2,ip,ks)
      do 200 k = 1, mxyp1
      nppp = kpic(k)
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
      ks(1) = kxr + kk
      ks(2) = kxl + kk
      ks(3) = kx + kr
      ks(4) = kxr + kr
      ks(5) = kxl + kr
      ks(6) = kx + kl
      ks(7) = kxr + kl
      ks(8) = kxl + kl
c loop over directions
      nh = ihole(1,1,k)
      noff = 0
      moff = 0
      if (ky.eq.1) then
         if (kx.gt.1) noff = mcll(3,kx-1)
      endif
      if (ky.eq.myp1) then
         if (kx.gt.1) moff = mclr(3,kx-1)
      endif
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 170 ii = 1, 8
c ip = number of particles coming from direction ii
      if (ks(ii).le.0) then
         if (ii.gt.6) noff = mcll(ii-6,ks(ii)+mx1)
         ip = mcll(ii-5,ks(ii)+mx1) - noff
      else if (ks(ii).gt.mxyp1) then
         if (ii.gt.3) moff = mclr(ii-3,ks(ii)-mxyp1)
         ip = mclr(ii-2,ks(ii)-mxyp1) - moff
      else
         if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
         ip = ncl(ii,ks(ii)) - ncoff
      endif
      do 160 j = 1, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,k)
c place overflow at end of array
      else
         j1 = nppp + 1
         nppp = j1
      endif
      if (j1.le.nppmx) then
         if (ks(ii).le.0) then
            do 130 i = 1, idimp
            ppart(i,j1,k) = rbufl(i,j+noff)
  130       continue
         else if (ks(ii).gt.mxyp1) then
            do 140 i = 1, idimp
            ppart(i,j1,k) = rbufr(i,j+moff)
  140       continue
         else
            do 150 i = 1, idimp
            ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
  150       continue
         endif
      else
         ist = 1
      endif
  160 continue
  170 continue
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
      if (ih.lt.nh) then
         ip = nh - ih
         do 190 j = 1, ip
         j1 = nppp - j + 1
         j2 = ihole(1,nh-j+2,k)
         if (j1.gt.j2) then
c move particle only if it is below current hole
            do 180 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
  180       continue
         endif
  190    continue
         nppp = nppp - ip
      endif
      kpic(k) = nppp
  200 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPCGUARD2XL(fxy,nyp,nx,ndim,nxe,nypmx)
c replicate extended periodic vector field in x direction
c linear interpolation, for distributed data
c nyp = number of primary (complete) gridpoints in particle partition
c nx = system length in x direction
c ndim = leading dimension of array fxy
c nxe = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells
      implicit none
      integer nyp, nx, ndim, nxe, nypmx
      real fxy
      dimension fxy(ndim,nxe,nypmx)
c local data
      integer i, k, myp1
c replicate edges of extended field
      myp1 = nyp + 1
      do 20 k = 1, myp1
      do 10 i = 1, ndim
      fxy(i,nx+1,k) = fxy(i,1,k)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPAGUARD2XL(q,nyp,nx,nxe,nypmx)
c accumulate extended periodic scalar field in x direction
c linear interpolation, for distributed data
c nyp = number of primary (complete) gridpoints in particle partition
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells
      implicit none
      integer nyp, nx, nxe, nypmx
      real q
      dimension q(nxe,nypmx)
c local data
      integer k, myp1
c accumulate edges of extended field
      myp1 = nyp + 1
      do 10 k = 1, myp1
      q(1,k) = q(1,k) + q(nx+1,k)
      q(nx+1,k) = 0.0
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPACGUARD2XL(cu,nyp,nx,ndim,nxe,nypmx)
c accumulate extended periodic vector field in x direction
c linear interpolation, for distributed data
c nyp = number of primary (complete) gridpoints in particle partition
c nx = system length in x direction
c ndim = leading dimension of array fxy
c nxe = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells
      implicit none
      real cu
      integer nyp, nx, ndim, nxe, nypmx
      dimension cu(ndim,nxe,nypmx)
c local data
      integer i, k, myp1
c accumulate edges of extended field
      myp1 = nyp + 1
      do 20 k = 1, myp1
      do 10 i = 1, ndim
      cu(i,1,k) = cu(i,1,k) + cu(i,nx+1,k)
      cu(i,nx+1,k) = 0.0
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,
     1kxp,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.  Zeros out z component.
c for distributed data.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp,nyhd,
c output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,kstrt,nyv,kxp,nyhd,
c output: fxy,we
c approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c q(k,j) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j) = x component of complex force/charge,
c fxy(2,k,j) = y component of complex force/charge,
c fxy(3,k,j) = zero,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c kxp = number of data values per block
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated.
c aimag(ffc(k,j)) = finite-size particle shape factor s
c real(ffc(k,j)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp, nyhd
      real ax, ay, affp, we
      complex q, fxy, ffc
      dimension q(nyv,kxp), fxy(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
      if (kstrt.gt.nxh) return
c prepare form factor array
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.0) then
         ffc(k,j) = cmplx(affp,1.0)
      else
         ffc(k,j) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 70
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,at1,at2,at3,zt1,zt2,wp)
!$OMP& REDUCTION(+:sum1)
      do 50 j = 1, kxps
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 40 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j))*aimag(ffc(k,j))
         at2 = dkx*at1
         at3 = dny*real(k - 1)*at1
         zt1 = cmplx(aimag(q(k,j)),-real(q(k,j)))
         zt2 = cmplx(aimag(q(k1,j)),-real(q(k1,j)))
         fxy(1,k,j) = at2*zt1
         fxy(2,k,j) = at3*zt1
         fxy(3,k,j) = zero
         fxy(1,k1,j) = at2*zt2
         fxy(2,k1,j) = -at3*zt2
         fxy(3,k1,j) = zero
         wp = wp + at1*(q(k,j)*conjg(q(k,j)) + q(k1,j)*conjg(q(k1,j)))
   40    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffc(1,j))*aimag(ffc(1,j))
         at3 = dkx*at1
         zt1 = cmplx(aimag(q(1,j)),-real(q(1,j)))
         fxy(1,1,j) = at3*zt1
         fxy(2,1,j) = zero
         fxy(3,1,j) = zero
         fxy(1,k1,j) = zero
         fxy(2,k1,j) = zero
         fxy(3,k1,j) = zero
         wp = wp + at1*(q(1,j)*conjg(q(1,j)))
      endif
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 60 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,1))*aimag(ffc(k,1))
         at2 = dny*real(k - 1)*at1
         zt1 = cmplx(aimag(q(k,1)),-real(q(k,1)))
         fxy(1,k,1) = zero
         fxy(2,k,1) = at2*zt1
         fxy(3,k,1) = zero
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
         fxy(3,k1,1) = zero
         wp = wp + at1*(q(k,1)*conjg(q(k,1)))
   60    continue
         k1 = nyh + 1
         fxy(1,1,1) = zero
         fxy(2,1,1) = zero
         fxy(3,1,1) = zero
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
         fxy(3,k1,1) = zero
      endif
      sum1 = sum1 + wp
   70 continue
      we = real(nx)*real(ny)*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPCUPERP2(cu,nx,ny,kstrt,nyv,kxp)
c this subroutine calculates the transverse current in fourier space
c input: all, output: cu
c approximate flop count is: 36*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
c and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
c cu(i,k,j) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp = number of data values per block
      implicit none
      integer nx, ny, kstrt, nyv, kxp
      complex cu
      dimension cu(3,nyv,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky, dkx2, at1
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c calculate transverse part of current
      if (kstrt.gt.nxh) return
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dkx2,dky,at1,zt1)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = 1.0/(dky*dky + dkx2)
         zt1 = at1*(dkx*cu(1,k,j) + dky*cu(2,k,j))
         cu(1,k,j) = cu(1,k,j) - dkx*zt1
         cu(2,k,j) = cu(2,k,j) - dky*zt1
         zt1 = at1*(dkx*cu(1,k1,j) - dky*cu(2,k1,j))
         cu(1,k1,j) = cu(1,k1,j) - dkx*zt1
         cu(2,k1,j) = cu(2,k1,j) + dky*zt1
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         cu(1,1,j) = zero
         cu(1,k1,j) = zero
         cu(2,k1,j) = zero
      endif
   20 continue
!$OMP END PARALLEL DO
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         cu(2,k,1) = zero
         cu(1,k1,1) = zero
         cu(2,k1,1) = zero
   30    continue
         k1 = nyh + 1
         cu(1,1,1) = zero
         cu(2,1,1) = zero
         cu(1,k1,1) = zero
         cu(2,k1,1) = zero
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MIPPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field with periodic boundary conditions for distributed data.
c input: cu,ffc,ci,nx,ny,kstrt,nyv,kxp,jblok,nyhd, output: bxy,wm
c approximate flop count is: 85*nxc*nyc + 36*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0,
c bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
c cu(i,k,j) = i-th component of complex current density and
c bxy(i,k,j) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffc(k,j)) = finite-size particle shape factor s
c real(ffc(k,j)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real ci, wm
      complex cu, bxy, ffc
      dimension cu(3,nyv,kxp), bxy(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real ci2, dnx, dny, dkx, dky, at1, at2, at3
      complex zero, zt1, zt2, zt3
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 40
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,at1,at2,at3,zt1,zt2,zt3,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = ci2*real(ffc(k,j))
         at2 = dky*at1
         at3 = dkx*at1
         at1 = at1*aimag(ffc(k,j))
         zt1 = cmplx(-aimag(cu(3,k,j)),real(cu(3,k,j)))
         zt2 = cmplx(-aimag(cu(2,k,j)),real(cu(2,k,j)))
         zt3 = cmplx(-aimag(cu(1,k,j)),real(cu(1,k,j)))
         bxy(1,k,j) = at2*zt1
         bxy(2,k,j) = -at3*zt1
         bxy(3,k,j) = at3*zt2 - at2*zt3
         zt1 = cmplx(-aimag(cu(3,k1,j)),real(cu(3,k1,j)))
         zt2 = cmplx(-aimag(cu(2,k1,j)),real(cu(2,k1,j)))
         zt3 = cmplx(-aimag(cu(1,k1,j)),real(cu(1,k1,j)))
         bxy(1,k1,j) = -at2*zt1
         bxy(2,k1,j) = -at3*zt1
         bxy(3,k1,j) = at3*zt2 + at2*zt3
         wp = wp + at1*(cu(1,k,j)*conjg(cu(1,k,j))                      
     1   + cu(2,k,j)*conjg(cu(2,k,j)) + cu(3,k,j)*conjg(cu(3,k,j))   
     2   + cu(1,k1,j)*conjg(cu(1,k1,j)) + cu(2,k1,j)*conjg(cu(2,k1,j))  
     3   + cu(3,k1,j)*conjg(cu(3,k1,j)))
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = ci2*real(ffc(1,j))
         at2 = dkx*at1
         at1 = at1*aimag(ffc(1,j))
         zt1 = cmplx(-aimag(cu(3,1,j)),real(cu(3,1,j)))
         zt2 = cmplx(-aimag(cu(2,1,j)),real(cu(2,1,j)))
         bxy(1,1,j) = zero
         bxy(2,1,j) = -at2*zt1
         bxy(3,1,j) = at2*zt2
         bxy(1,k1,j) = zero
         bxy(2,k1,j) = zero
         bxy(3,k1,j) = zero
         wp = wp + at1*(cu(1,1,j)*conjg(cu(1,1,j))                      
     1   + cu(2,1,j)*conjg(cu(2,1,j)) + cu(3,1,j)*conjg(cu(3,1,j)))
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = ci2*real(ffc(k,1))
         at2 = dky*at1
         at1 = at1*aimag(ffc(k,1))
         zt1 = cmplx(-aimag(cu(3,k,1)),real(cu(3,k,1)))
         zt2 = cmplx(-aimag(cu(1,k,1)),real(cu(1,k,1)))
         bxy(1,k,1) = at2*zt1
         bxy(2,k,1) = zero
         bxy(3,k,1) = -at2*zt2
         bxy(1,k1,1) = zero
         bxy(2,k1,1) = zero
         bxy(3,k1,1) = zero
         wp = wp + at1*(cu(1,k,1)*conjg(cu(1,k,1))                      
     1   + cu(2,k,1)*conjg(cu(2,k,1)) + cu(3,k,1)*conjg(cu(3,k,1)))
   30    continue
         k1 = nyh + 1
         bxy(1,1,1) = zero
         bxy(2,1,1) = zero
         bxy(3,1,1) = zero
         bxy(1,k1,1) = zero
         bxy(2,k1,1) = zero
         bxy(3,k1,1) = zero
      endif
      sum1 = sum1 + wp
   40 continue
      wm = real(nx)*real(ny)*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,kstrt,
     1nyv,kxp,nyhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions.
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ez(ky=pi) = 0,
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c and similarly for bx, by, bz.
c cu(i,k,j) = i-th component of complex current density and
c exy(i,k,j) = i-th component of complex electric field,
c bxy(i,k,j) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c aimag(ffc(k,j)) = finite-size particle shape factor s
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
c nx/ny = system length in x/y direction
c kxp = number of data values per block
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real affp, ci, dt, wf, wm
      complex exy, bxy, cu, ffc
      dimension exy(3,nyv,kxp), bxy(3,nyv,kxp), cu(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dth, c2, cdt, adt, anorm, dkx, dky, afdt
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      double precision wp, ws, sum1, sum2
      if (ci.le.0.0) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dth = 0.5*dt
      c2 = 1.0/(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      sum1 = 0.0d0
      sum2 = 0.0d0
      if (kstrt.gt.nxh) go to 40
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,afdt,zt1,zt2,zt3,zt4,zt5,zt6,
!$OMP& zt7,zt8,zt9,ws,wp)
!$OMP& REDUCTION(+:sum1,sum2)
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      ws = 0.0d0
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         afdt = adt*aimag(ffc(k,j))
c update magnetic field half time step, ky > 0
         zt1 = cmplx(-aimag(exy(3,k,j)),real(exy(3,k,j)))
         zt2 = cmplx(-aimag(exy(2,k,j)),real(exy(2,k,j)))
         zt3 = cmplx(-aimag(exy(1,k,j)),real(exy(1,k,j)))
         zt4 = bxy(1,k,j) - dth*(dky*zt1)
         zt5 = bxy(2,k,j) + dth*(dkx*zt1)
         zt6 = bxy(3,k,j) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt2 = cmplx(-aimag(zt5),real(zt5))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,k,j) + cdt*(dky*zt1) - afdt*cu(1,k,j)
         zt8 = exy(2,k,j) - cdt*(dkx*zt1) - afdt*cu(2,k,j)
         zt9 = exy(3,k,j) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,k,j)
c update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt2 = cmplx(-aimag(zt8),real(zt8))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         exy(1,k,j) = zt7
         exy(2,k,j) = zt8
         exy(3,k,j) = zt9
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               
     1                  + zt9*conjg(zt9))
         zt4 = zt4 - dth*(dky*zt1)
         zt5 = zt5 + dth*(dkx*zt1)
         zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
         bxy(1,k,j) = zt4
         bxy(2,k,j) = zt5
         bxy(3,k,j) = zt6
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               
     1                  + zt6*conjg(zt6))
c update magnetic field half time step, ky < 0
         zt1 = cmplx(-aimag(exy(3,k1,j)),real(exy(3,k1,j)))
         zt2 = cmplx(-aimag(exy(2,k1,j)),real(exy(2,k1,j)))
         zt3 = cmplx(-aimag(exy(1,k1,j)),real(exy(1,k1,j)))
         zt4 = bxy(1,k1,j) + dth*(dky*zt1)
         zt5 = bxy(2,k1,j) + dth*(dkx*zt1)
         zt6 = bxy(3,k1,j) - dth*(dkx*zt2 + dky*zt3)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt2 = cmplx(-aimag(zt5),real(zt5))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,k1,j) - cdt*(dky*zt1) - afdt*cu(1,k1,j)
         zt8 = exy(2,k1,j) - cdt*(dkx*zt1) - afdt*cu(2,k1,j)
         zt9 = exy(3,k1,j) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,k1,j)
c update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt2 = cmplx(-aimag(zt8),real(zt8))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         exy(1,k1,j) = zt7
         exy(2,k1,j) = zt8
         exy(3,k1,j) = zt9
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               
     1                  + zt9*conjg(zt9))
         zt4 = zt4 + dth*(dky*zt1)
         zt5 = zt5 + dth*(dkx*zt1)
         zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
         bxy(1,k1,j) = zt4
         bxy(2,k1,j) = zt5
         bxy(3,k1,j) = zt6
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               
     1                  + zt6*conjg(zt6))
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         afdt = adt*aimag(ffc(1,j))
c update magnetic field half time step
         zt1 = cmplx(-aimag(exy(3,1,j)),real(exy(3,1,j)))
         zt2 = cmplx(-aimag(exy(2,1,j)),real(exy(2,1,j)))
         zt5 = bxy(2,1,j) + dth*(dkx*zt1)
         zt6 = bxy(3,1,j) - dth*(dkx*zt2)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt2 = cmplx(-aimag(zt5),real(zt5))
         zt8 = exy(2,1,j) - cdt*(dkx*zt1) - afdt*cu(2,1,j)
         zt9 = exy(3,1,j) + cdt*(dkx*zt2) - afdt*cu(3,1,j)
c update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt2 = cmplx(-aimag(zt8),real(zt8))
         exy(1,1,j) = zero
         exy(2,1,j) = zt8
         exy(3,1,j) = zt9
         ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
         zt5 = zt5 + dth*(dkx*zt1)
         zt6 = zt6 - dth*(dkx*zt2)
         bxy(1,1,j) = zero
         bxy(2,1,j) = zt5
         bxy(3,1,j) = zt6
         wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
         bxy(1,k1,j) = zero
         bxy(2,k1,j) = zero
         bxy(3,k1,j) = zero
         exy(1,k1,j) = zero
         exy(2,k1,j) = zero
         exy(3,k1,j) = zero
      endif
      sum1 = sum1 + ws
      sum2 = sum2 + wp
   20 continue
!$OMP END PARALLEL DO
      ws = 0.0d0
      wp = 0.0d0
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         afdt = adt*aimag(ffc(k,1))
c update magnetic field half time step
         zt1 = cmplx(-aimag(exy(3,k,1)),real(exy(3,k,1)))
         zt3 = cmplx(-aimag(exy(1,k,1)),real(exy(1,k,1)))
         zt4 = bxy(1,k,1) - dth*(dky*zt1)
         zt6 = bxy(3,k,1) + dth*(dky*zt3)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,k,1) + cdt*(dky*zt1) - afdt*cu(1,k,1)
         zt9 = exy(3,k,1) - cdt*(dky*zt3) - afdt*cu(3,k,1)
c update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         exy(1,k,1) = zt7
         exy(2,k,1) = zero
         exy(3,k,1) = zt9
         ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
         zt4 = zt4 - dth*(dky*zt1)
         zt6 = zt6 + dth*(dky*zt3)
         bxy(1,k,1) = zt4
         bxy(2,k,1) = zero
         bxy(3,k,1) = zt6
         wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
         bxy(1,k1,1) = zero
         bxy(2,k1,1) = zero
         bxy(3,k1,1) = zero
         exy(1,k1,1) = zero
         exy(2,k1,1) = zero
         exy(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         bxy(1,1,1) = zero
         bxy(2,1,1) = zero
         bxy(3,1,1) = zero
         exy(1,1,1) = zero
         exy(2,1,1) = zero
         exy(3,1,1) = zero
         bxy(1,k1,1) = zero
         bxy(2,k1,1) = zero
         bxy(3,k1,1) = zero
         exy(1,k1,1) = zero
         exy(2,k1,1) = zero
         exy(3,k1,1) = zero
      endif
      sum1 = sum1 + ws
      sum2 = sum2 + wp
   40 continue
      wf = real(nx)*real(ny)*sum1
      wm = real(nx)*real(ny)*c2*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,nyhd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign < 0
c includes additional smoothing
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp, nyhd
      complex fxy, exy, ffc
      dimension fxy(3,nyv,kxp), exy(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
c local data
      integer i, nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real at1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      if (kstrt.gt.nxh) return
c add the fields
      if (isign.gt.0) then
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(i,j,k,k1,at1)
         do 40 j = 1, kxps
         do 20 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j))
         do 10 i = 1, 3
         fxy(i,k,j) = fxy(i,k,j) + exy(i,k,j)*at1
         fxy(i,k1,j) = fxy(i,k1,j) + exy(i,k1,j)*at1
   10    continue
   20    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j))
         do 30 i = 1, 3
         fxy(i,1,j) = fxy(i,1,j) + exy(i,1,j)*at1
         fxy(i,k1,j) = fxy(i,k1,j) + exy(i,k1,j)*at1
   30    continue
   40    continue
!$OMP END PARALLEL DO
c copy the fields
      else if (isign.lt.0) then
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(i,j,k,k1,at1)
         do 80 j = 1, kxps
         do 60 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j))
         do 50 i = 1, 3
         fxy(i,k,j) = exy(i,k,j)*at1
         fxy(i,k1,j) = exy(i,k1,j)*at1
   50    continue
   60    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j))
         do 70 i = 1, 3
         fxy(i,1,j) = exy(i,1,j)*at1
         fxy(i,k1,j) = exy(i,k1,j)*at1
   70    continue
   80    continue
!$OMP END PARALLEL DO
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
c this subroutine calculates tables needed by a two dimensional
c real to complex fast fourier transform and its inverse.
c input: indx, indy, nxhyd, nxyhd
c output: mixup, sct
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, nxhyd, nxyhd
      integer mixup
      complex sct
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nxy, nxhy, nxyh
      integer j, k, lb, ll, jb, it
      real dnxy, arg
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/real(nxy)
      do 30 j = 1, nxyh
      arg = dnxy*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPPFFT2RM(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,   
     1indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
c parallelized with OpenMP
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv, kxp, kyp
      integer kypd, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(nxvh,kypd), g(nyv,kxp)
      dimension bs(kxp,kyp), br(kxp,kyp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi, ks, kxpp, kypp
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      ks = kstrt - 1
      kxpp = min(kxp,max(0,nxh-kxp*ks))
      kypp = min(kyp,max(0,ny-kyp*ks))
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PPFFT2RMXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,   
     1nxvh,kypd,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp,  
     1kypd)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PPFFT2RMXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv
     1,kxp,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,   
     1kypd,kxp)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp
     1,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PPFFT2RMXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv
     1,kxp,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,kypd, 
     1kxp)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PPFFT2RMXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,   
     1nxvh,kypd,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPPFFT2RM3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,  
     1indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
c parallelized with OpenMP
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv, kxp, kyp
      integer kypd, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(3,nxvh,kypd), g(3,nyv,kxp)
      dimension bs(3,kxp,kyp), br(3,kxp,kyp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi, ks, kxpp, kypp
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      ks = kstrt - 1
      kxpp = min(kxp,max(0,nxh-kxp*ks))
      kypp = min(kyp,max(0,ny-kyp*ks))
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PPFFT2RM3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,  
     1nxvh,kypd,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,3,nxvh,nyv,kxp
     1,kypd)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PPFFT2RM3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,  
     1nyv,kxp,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,3,nyv,nxvh,
     1kypd,kxp)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,3,nxvh,nyv,
     1kxp,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PPFFT2RM3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,  
     1nyv,kxp,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,3,nyv,nxvh,   
     1kypd,kxp)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PPFFT2RM3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,  
     1nxvh,kypd,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2RMXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,
     1nxvh,kypd,nxhyd,nxyhd)
c this subroutine performs the x part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic, with OpenMP,
c for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kypd = second dimension of f
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k) = mode j-1,kk-1, where kk = k + kyp*(kstrt - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1) = real part of mode nx/2,0 on mode kstrt=0
c imaginary part of f(1,1) = real part of mode nx/2,ny/2
c on mode kstrt=(ny/2)/kyp
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, kstrt, nxvh, kypi, kypp, kypd
      integer nxhyd, nxyhd, mixup
      complex f, sct
      dimension f(nxvh,kypd)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, kypt, j, k, nrx
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nrxb
      real ani
      complex s, t, t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 70
c inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t,t1)
      do 60 i = kypi, kypt
c bit-reverse array elements in x
      do 10 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t = f(j1,i)
         f(j1,i) = f(j,i)
         f(j,i) = t
      endif
   10 continue
c then transform in x
      do 40 m = 1, indx1
      ns = 2**(m - 1)
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
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
   20 continue
   30 continue
   40 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 50 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,i))
      s = f(j,i) + t
      t = (f(j,i) - t)*t1
      f(j,i) = ani*(s + t)
      f(nxh2-j,i) = ani*conjg(s - t)
   50 continue
      f(1,i) = 2.0*ani*cmplx(real(f(1,i)) + aimag(f(1,i)),              
     1                       real(f(1,i)) - aimag(f(1,i)))
      if (nxhh.gt.0) f(nxhh+1,i) = 2.0*ani*conjg(f(nxhh+1,i))
   60 continue
!$OMP END PARALLEL DO
      return
c forward fourier transform
   70 nrxb = nxhy/nxh
      nrx = nxhy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t,t1)
      do 130 i = kypi, kypt
c scramble coefficients
      kmr = nxy/nx
      do 80 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,i))
      s = f(j,i) + t
      t = (f(j,i) - t)*t1
      f(j,i) = s + t
      f(nxh2-j,i) = conjg(s - t)
   80 continue
      f(1,i) = cmplx(real(f(1,i)) + aimag(f(1,i)),                      
     1               real(f(1,i)) - aimag(f(1,i)))
      if (nxhh.gt.0) f(nxhh+1,i) = 2.0*conjg(f(nxhh+1,i))
c bit-reverse array elements in x
      do 90 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t = f(j1,i)
         f(j1,i) = f(j,i)
         f(j,i) = t
      endif
   90 continue
c then transform in x
      do 120 m = 1, indx1
      ns = 2**(m - 1)
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
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
  100 continue
  110 continue
  120 continue
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2RMXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,
     1nyv,kxp,nxhyd,nxyhd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic, with OpenMP,
c for data which is distributed in blocks
c for isign = (-1,1), input: all, output: g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c g(m,n) = sum(g(k,j)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(k,j) = sum(g(m,n)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c kxp = number of x indices per block
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g
c kxp = number of data values per block in x
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(k,j) = mode jj-1,k-1, where jj = j + kxp*(kstrt - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1) = real part of mode nx/2,ny/2
c on node kstrt=0
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxp
      integer nxhyd, nxyhd, mixup
      complex g, sct
      dimension g(nyv,kxp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nryb
      complex s, t
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 70
c inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t)
      do 50 i = kxpi, kxpt
c bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t = g(k1,i)
         g(k1,i) = g(k,i)
         g(k,i) = t
      endif
   10 continue
c then transform in y
      do 40 m = 1, indy
      ns = 2**(m - 1)
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
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 60 k = 2, nyh
         s = g(ny2-k,1)
         g(ny2-k,1) = 0.5*cmplx(aimag(g(k,1) + s),real(g(k,1) - s))
         g(k,1) = 0.5*cmplx(real(g(k,1) + s),aimag(g(k,1) - s))
   60    continue
      endif
      return
c forward fourier transform
   70 nryb = nxhy/ny
      nry = nxy/ny
c scramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 80 k = 2, nyh
         s = cmplx(aimag(g(ny2-k,1)),real(g(ny2-k,1)))
         g(ny2-k,1) = conjg(g(k,1) - s)
         g(k,1) = g(k,1) + s
   80    continue
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t)
      do 130 i = kxpi, kxpt
c bit-reverse array elements in y
      do 90 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t = g(k1,i)
         g(k1,i) = g(k,i)
         g(k,i) = t
      endif
   90 continue
c then transform in y
      do 120 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
  100 continue
  110 continue
  120 continue
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2RM3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp
     1,nxvh,kypd,nxhyd,nxyhd)
c this subroutine performs the x part of 3 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic, with OpenMP,
c for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kypd = second dimension of f
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(1:3,j,k) = mode j-1,kk-1, where kk = k + kyp*(kstrt - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1:3,1,k) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1:3,1,1) = real part of mode nx/2,0
c on mode kstrt=0
c imaginary part of f(1:3,1,1) = real part of mode nx/2,ny/2
c on mode kstrt=(ny/2)/kyp
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, nxvh, kypi, kypp
      integer kypd, nxhyd, nxyhd
      complex f, sct
      dimension f(3,nxvh,kypd)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, kypt, j, k, nrx
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nrxb
      real ani, at1, at2
      complex s, t, t1, t2, t3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 100
c inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,s,t,t1,t2,t3)
      do 90 i = kypi, kypt
c swap complex components
      do 10 j = 1, nxh
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(real(f(2,j,i)),aimag(f(3,j,i)))
      at2 = aimag(f(2,j,i))
      f(2,j,i) = cmplx(aimag(f(1,j,i)),at1)
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
   10 continue
c bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i)
         t2 = f(2,j1,i)
         t3 = f(3,j1,i)
         f(1,j1,i) = f(1,j,i)
         f(2,j1,i) = f(2,j,i)
         f(3,j1,i) = f(3,j,i)
         f(1,j,i) = t1
         f(2,j,i) = t2
         f(3,j,i) = t3
      endif
   20 continue
c then transform in x
      do 50 m = 1, indx1
      ns = 2**(m - 1)
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
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      t3 = s*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(3,j2,i) = f(3,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
      f(3,j1,i) = f(3,j1,i) + t3
   30 continue
   40 continue
   50 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 70 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 60 k = 1, 3
      t = conjg(f(k,nxh2-j,i))
      s = f(k,j,i) + t
      t = (f(k,j,i) - t)*t1
      f(k,j,i) = ani*(s + t)
      f(k,nxh2-j,i) = ani*conjg(s - t)
   60 continue
   70 continue
      do 80 k = 1, 3
      f(k,1,i) = 2.0*ani*cmplx(real(f(k,1,i)) + aimag(f(k,1,i)),        
     1                         real(f(k,1,i)) - aimag(f(k,1,i)))
      if (nxhh.gt.0) f(k,nxhh+1,i) = 2.0*ani*conjg(f(k,nxhh+1,i))
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
c forward fourier transform
  100 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,s,t,t1,t2,t3)
      do 190 i = kypi, kypt
c scramble coefficients
      kmr = nxy/nx
      do 120 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 k = 1, 3
      t = conjg(f(k,nxh2-j,i))
      s = f(k,j,i) + t
      t = (f(k,j,i) - t)*t1
      f(k,j,i) = s + t
      f(k,nxh2-j,i) = conjg(s - t)
  110 continue
  120 continue
      do 130 k = 1, 3
      f(k,1,i) = cmplx(real(f(k,1,i)) + aimag(f(k,1,i)),                
     1                 real(f(k,1,i)) -aimag(f(k,1,i)))
      if (nxhh.gt.0) f(k,nxhh+1,i) = 2.0*conjg(f(k,nxhh+1,i))
  130 continue
c bit-reverse array elements in x
      do 140 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i)
         t2 = f(2,j1,i)
         t3 = f(3,j1,i)
         f(1,j1,i) = f(1,j,i)
         f(2,j1,i) = f(2,j,i)
         f(3,j1,i) = f(3,j,i)
         f(1,j,i) = t1
         f(2,j,i) = t2
         f(3,j,i) = t3
      endif
  140 continue
c then transform in x
      do 170 m = 1, indx1
      ns = 2**(m - 1)
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
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      t3 = s*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(3,j2,i) = f(3,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
      f(3,j1,i) = f(3,j1,i) + t3
  150 continue
  160 continue
  170 continue
c swap complex components
      do 180 j = 1, nxh
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(aimag(f(2,j,i)),aimag(f(3,j,i)))
      at2 = real(f(2,j,i))
      f(2,j,i) = cmplx(at1,aimag(f(1,j,i)))
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2RM3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp
     1,nyv,kxp,nxhyd,nxyhd)
c this subroutine performs the y part of 3 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c using complex arithmetic, with OpenMP,
c for data which is distributed in blocks
c for isign = (-1,1), input: all, output: g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c g(1:3,m,n) = sum(g(1:3,k,j)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(1:3,k,j) = sum(g(1:3,m,n)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g
c kxp = number of data values per block in x
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(1:3,k,j) = mode jj-1,k-1, where jj = j + kxp*(kstrt - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(1:3,k,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1:3,1,1) = real part of mode nx/2,0 and
c imaginary part of g(1:3,ny/2+1,1) = real part of mode nx/2,ny/2
c on node kstrt=0
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxpi, kxpp, nyv
      integer kxp, nxhyd, nxyhd
      complex g, sct
      dimension g(3,nyv,kxp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nryb
      complex s, t1, t2, t3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 80
c inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1,t2,t3)
      do 50 i = kxpi, kxpt
c bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = g(1,k1,i)
         t2 = g(2,k1,i)
         t3 = g(3,k1,i)
         g(1,k1,i) = g(1,k,i)
         g(2,k1,i) = g(2,k,i)
         g(3,k1,i) = g(3,k,i)
         g(1,k,i) = t1
         g(2,k,i) = t2
         g(3,k,i) = t3
      endif
   10 continue
c then transform in y
      do 40 m = 1, indy
      ns = 2**(m - 1)
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
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      t3 = s*g(3,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(3,j2,i) = g(3,j1,i) - t3
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
      g(3,j1,i) = g(3,j1,i) + t3
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 70 k = 2, nyh
         do 60 j = 1, 3
         s = g(j,ny2-k,1)
         g(j,ny2-k,1) = 0.5*cmplx(aimag(g(j,k,1) + s),                  
     1                            real(g(j,k,1) - s))
         g(j,k,1) = 0.5*cmplx(real(g(j,k,1) + s),aimag(g(j,k,1) - s))
   60    continue
   70    continue
      endif
      return
c forward fourier transform
   80 nryb = nxhy/ny
      nry = nxy/ny
c scramble modes kx = 0, nx/2
      if ((ks.eq.0).and.(kxpi.eq.1)) then
         do 100 k = 2, nyh
         do 90 j = 1, 3
         s = cmplx(aimag(g(j,ny2-k,1)),real(g(j,ny2-k,1)))
         g(j,ny2-k,1) = conjg(g(j,k,1) - s)
         g(j,k,1) = g(j,k,1) + s
   90    continue
  100    continue
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,s,t1,t2,t3)
      do 150 i = kxpi, kxpt
c bit-reverse array elements in y
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = g(1,k1,i)
         t2 = g(2,k1,i)
         t3 = g(3,k1,i)
         g(1,k1,i) = g(1,k,i)
         g(2,k1,i) = g(2,k,i)
         g(3,k1,i) = g(3,k,i)
         g(1,k,i) = t1
         g(2,k,i) = t2
         g(3,k,i) = t3
      endif
  110 continue
c then transform in y
      do 140 m = 1, indy
      ns = 2**(m - 1)
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
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      t3 = s*g(3,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(3,j2,i) = g(3,j1,i) - t3
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
      g(3,j1,i) = g(3,j1,i) + t3
  120 continue
  130 continue
  140 continue
  150 continue
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
      subroutine PPPCOPYOUT(part,ppart,kpic,npp,npmax,nppmx,idimp,mxyp1,
     1irc)
c for 2d code, this subroutine copies segmented particle data ppart to
c the array part with original tiled layout
c spatial decomposition in y direction
c input: all except part, npp, irc, output: part, npp, irc
c part(i,j) = i-th coordinate for particle j in partition
c ppart(i,j,k) = i-th coordinate for particle j in partition in tile k
c kpic = number of particles per tile
c npp = number of particles in partition
c npmax = maximum number of particles in each partition
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c mxyp1 = total number of tiles in partition
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer npp, npmax, nppmx, idimp, mxyp1, irc
      real part, ppart
      integer kpic
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1)
c local data
      integer i, j, k, npoff, nppp, ne, ierr
      npoff = 0
      ierr = 0
c loop over tiles
      do 30 k = 1, mxyp1
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
