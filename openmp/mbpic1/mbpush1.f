c Fortran Library for Skeleton 1-2/2D Electromagnetic PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,idimp,nop,nx, 
     1ipbc)
c for 1-2/2d code, this subroutine calculates initial particle
c co-ordinate and velocity, with uniform density and maxwellian
c velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx = number of particles distributed in x direction
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, idimp, nop, nx, ipbc
      real part, vtx, vty, vtz, vdx, vdy,v dz
      dimension part(idimp,nop)
c local data
      integer j
      real edgelx, at1, sum1, sum2, sum3
      double precision dsum1, dsum2, dsum3
      double precision ranorm
c set boundary values
      edgelx = 0.0
      at1 = real(nx)/real(npx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
c uniform density profile
      do 10 j = 1, npx
      part(1,j) = edgelx + at1*(real(j) - .5)
   10 continue
c maxwellian velocity distribution
      do 20 j = 1, npx
      part(2,j) = vtx*ranorm()
      part(3,j) = vty*ranorm()
      part(4,j) = vtz*ranorm()
   20 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 30 j = 1, npx
      dsum1 = dsum1 + part(2,j)
      dsum2 = dsum2 + part(3,j)
      dsum3 = dsum3 + part(4,j)
   30 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1.0/real(npx)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 40 j = 1, npx
      part(2,j) = part(2,j) - sum1
      part(3,j) = part(3,j) - sum2
      part(4,j) = part(4,j) - sum3
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DBLKP1L(part,kpic,nppmx,idimp,nop,mx,mx1,irc)
c this subroutine finds the maximum number of particles in each tile of
c mx to calculate size of segmented particle array ppart
c linear interpolation
c input: all except kpic, nppmx, output: kpic, nppmx
c part = input particle array
c part(1,n) = position x of particle n
c kpic = output number of particles per tile
c nppmx = return maximum number of particles in tile
c idimp = size of phase space = 2
c nop = number of particles
c mx = number of grids in sorting cell in x
c mx1 = (system length in x direction - 1)/mx + 1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kpic, nppmx, idimp, nop, mx, mx1, irc
      real part
      dimension part(idimp,nop), kpic(mx1)
c local data
      integer j, k, n, isum, ist, npx, ierr
      ierr = 0
c clear counter array
      do 10 k = 1, mx1
      kpic(k) = 0
   10 continue
c find how many particles in each tile
      do 20 j = 1, nop
      n = part(1,j)
      n = n/mx + 1
      if (n.le.mx1) then
         kpic(n) = kpic(n) + 1
      else
         ierr = max(ierr,n-mx1)
      endif
   20 continue
c find maximum
      isum = 0
      npx = 0
      do 30 k = 1, mx1
      ist = kpic(k)
      npx = max(npx,ist)
      isum = isum + ist
   30 continue
      nppmx = npx
c check for errors
      if (ierr.gt.0) then
         irc = ierr
      else if (isum.ne.nop) then
         irc = -1
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPMOVIN1L(part,ppart,kpic,nppmx,idimp,nop,mx,mx1,irc)
c this subroutine sorts particles by x grid in tiles of mx and copies
c to segmented array ppart
c linear interpolation
c input: all except ppart, kpic, output: ppart, kpic, irc
c part/ppart = input/output particle arrays
c part(1,n) = position x of particle n
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = velocity vx of particle n in tile m
c kpic = output number of particles per tile
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 2
c nop = number of particles
c mx = number of grids in sorting cell in x
c mx1 = (system length in x direction - 1)/mx + 1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kpic, nppmx, idimp, nop, mx, mx1, irc
      real part, ppart
      dimension part(idimp,nop), ppart(idimp,nppmx,mx1)
      dimension kpic(mx1)
c local data
      integer i, j, k, n, ip, ierr
      ierr = 0
c clear counter array
      do 10 k = 1, mx1
      kpic(k) = 0
   10 continue
c find addresses of particles at each tile and reorder particles
      do 30 j = 1, nop
      n = part(1,j)
      n = n/mx + 1
      ip = kpic(n) + 1
      if (ip.le.nppmx) then
         do 20 i = 1, idimp
         ppart(i,ip,n) = part(i,j)
   20    continue
      else
         ierr = max(ierr,ip-nppmx)
      endif
      kpic(n) = ip
   30 continue
      if (ierr.gt.0) irc = ierr
      return
      end
c-----------------------------------------------------------------------
      subroutine PPCHECK1L(ppart,kpic,idimp,nppmx,nx,mx,mx1,irc)
c this subroutine performs a sanity check to make sure particles sorted
c by x grid in tiles of mx, are all within bounds.
c tiles are assumed to be arranged in 1D linear memory
c input: all except irc
c output: irc
c ppart(1,n,k) = position x of particle n in tile k
c kpic(k) = number of reordered output particles in tile k
c idimp = size of phase space = 2
c nppmx = maximum number of particles in tile
c nx = system length in x direction
c mx = number of grids in sorting cell in x
c mx1 = (system length in x direction - 1)/mx + 1
c irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, mx, mx1, irc
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension kpic(mx1)
c local data
      integer noff, npp, j, k, ist, nn
      real edgelx, edgerx, dx
c loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,noff,npp,nn,ist,edgelx,edgerx,dx)
      do 20 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
c loop over particles in tile
      do 10 j = 1, npp
      dx = ppart(1,j,k)
c find particles going out of bounds
      ist = 0
      if (dx.lt.edgelx) ist = 1
      if (dx.ge.edgerx) ist = 2
      if (ist.gt.0) irc = k
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ek,idimp,
     1nppmx,nx,mx,nxv,mx1,ipbc)
c for 1-2/2d code, this subroutine updates particle co-ordinate and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 78 flops/particle, 1 divide, 14 loads, 4 stores
c input: all, output: ppart, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fx(x(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fy(x(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fz(x(t))*dt)
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
c omy = (q/m)*by(x(t)), and omz = (q/m)*bz(x(t)).
c position equations used are:
c x(t+dt) = x(t) + vx(t+dt/2)*dt
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
c where n = nearest grid point and dx = x-n
c similarly for fy(x), fz(x), by(x), bz(x)
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = velocity vx of particle n in tile m
c ppart(3,n,m) = velocity vy of particle n in tile m
c ppart(4,n,m) = velocity vz of particle n in tile m
c fxyz(1,j) = x component of force/charge at grid (j)
c fxyz(2,j) = y component of force/charge at grid (j)
c fxyz(3,j) = z component of force/charge at grid (j)
c that is, convolution of electric field over particle shape
c byz(1,j) = y component of magnetic field at grid (j)
c byz(2,j) = z component of magnetic field at grid (j)
c that is, the convolution of magnetic field over particle shape
c kpic = number of particles per tile
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx = system length in x direction
c mx = number of grids in sorting cell in x
c nxv = second dimension of field arrays, must be >= nx+1
c mx1 = (system length in x direction - 1)/mx + 1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ipbc
      real omx, qbm, dt, dtc, ek
      real ppart, fxyz, byz
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension kpic(mx1)
c local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real qtmh, edgelx, edgerx, dxp, amx, x, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real sfxyz, sbyz
      dimension sfxyz(3,MXV), sbyz(2,MXV)
c     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      sum2 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c error if local array is too small
c     if (mx.ge.MXV) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy,acz, 
!$OMP& omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,
!$OMP& rot9,sum1,sfxyz,sbyz)
!$OMP& REDUCTION(+:sum2)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
c load local fields from global array
      nn = min(mx,nx-noff) + 1
      do 10 j = 1, nn
      sfxyz(1,j) = fxyz(1,j+noff)
      sfxyz(2,j) = fxyz(2,j+noff)
      sfxyz(3,j) = fxyz(3,j+noff)
   10 continue
      do 20 j = 1, nn
      sbyz(1,j) = byz(1,j+noff)
      sbyz(2,j) = byz(2,j+noff)
   20 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 30 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
c find electric field
      dx = amx*sfxyz(1,nn) + dxp*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + dxp*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + dxp*sfxyz(3,nn+1)
c find magnetic field
      ox = omx
      oy = amx*sbyz(1,nn) + dxp*sbyz(1,nn+1)
      oz = amx*sbyz(2,nn) + dxp*sbyz(2,nn+1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(2,j,k) + dx
      acy = ppart(3,j,k) + dy
      acz = ppart(4,j,k) + dz
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
      ppart(2,j,k) = dx
      ppart(3,j,k) = dy
      ppart(4,j,k) = dz
c new position
      dx = x + dx*dtc
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(2,j,k) = -ppart(2,j,k)
         endif
      endif
c set new position
      ppart(1,j,k) = dx
   30 continue
      sum2 = sum2 + sum1
   40 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + .5*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine GBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt,  
     1dtc,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
c for 1-2/2d code, this subroutine updates particle co-ordinate and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 78 flops/particle, 1 divide, 14 loads, 4 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fx(x(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fy(x(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fz(x(t))*dt)
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
c omy = (q/m)*by(x(t)), and omz = (q/m)*bz(x(t)).
c position equations used are:
c x(t+dt) = x(t) + vx(t+dt/2)*dt
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
c where n = nearest grid point and dx = x-n
c similarly for fy(x), fz(x), by(x), bz(x)
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = velocity vx of particle n in tile m
c ppart(3,n,m) = velocity vy of particle n in tile m
c ppart(4,n,m) = velocity vz of particle n in tile m
c fxyz(1,j) = x component of force/charge at grid (j)
c fxyz(2,j) = y component of force/charge at grid (j)
c fxyz(3,j) = z component of force/charge at grid (j)
c that is, convolution of electric field over particle shape
c byz(1,j) = y component of magnetic field at grid (j)
c byz(2,j) = z component of magnetic field at grid (j)
c that is, the convolution of magnetic field over particle shape
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx = system length in x direction
c mx = number of grids in sorting cell in x
c nxv = second dimension of field arrays, must be >= nx+1
c mx1 = (system length in x direction - 1)/mx + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ntmax, irc
      real omx, qbm, dt, dtc, ek
      real ppart, fxyz, byz
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
c local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real qtmh, anx, edgelx, edgerx, dxp, amx, x, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real sfxyz, sbyz
      dimension sfxyz(3,MXV), sbyz(2,MXV)
c     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      anx = real(nx)
      sum2 = 0.0d0
c error if local array is too small
c     if (mx.ge.MXV) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy
!$OMP& ,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,
!$OMP& rot8,rot9,edgelx,edgerx,sum1,sfxyz,sbyz)
!$OMP& REDUCTION(+:sum2)
      do 50 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
c load local fields from global array
      nn = nn + 1
      do 10 j = 1, nn
      sfxyz(1,j) = fxyz(1,j+noff)
      sfxyz(2,j) = fxyz(2,j+noff)
      sfxyz(3,j) = fxyz(3,j+noff)
   10 continue
      do 20 j = 1, nn
      sbyz(1,j) = byz(1,j+noff)
      sbyz(2,j) = byz(2,j+noff)
   20 continue
c clear counters
      do 30 j = 1, 2
      ncl(j,k) = 0
   30 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 40 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
c find electric field
      dx = amx*sfxyz(1,nn) + dxp*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + dxp*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + dxp*sfxyz(3,nn+1)
c find magnetic field
      ox = omx
      oy = amx*sbyz(1,nn) + dxp*sbyz(1,nn+1)
      oz = amx*sbyz(2,nn) + dxp*sbyz(2,nn+1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(2,j,k) + dx
      acy = ppart(3,j,k) + dy
      acz = ppart(4,j,k) + dz
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
      ppart(2,j,k) = dx
      ppart(3,j,k) = dy
      ppart(4,j,k) = dz
c new position
      dx = x + dx*dtc
c find particles going out of bounds
      nn = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c nn = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         nn = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               nn = 1
            else
               dx = 0.0
            endif
         else
            nn = 1
         endif
      endif
c set new position
      ppart(1,j,k) = dx
c increment counters
      if (nn.gt.0) then
         ncl(nn,k) = ncl(nn,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = nn
         else
            nh = 1
         endif
      endif
   40 continue
      sum2 = sum2 + sum1
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + .5*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine GRBPPUSH13L(ppart,fxyz,byz,kpic,omx,qbm,dt,dtc,ci,ek,  
     1idimp,nppmx,nx,mx,nxv,mx1,ipbc)
c for 1-2/2d code, this subroutine updates particle co-ordinate and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 90 flops/particle, 4 divides, 2 sqrts, 14 loads, 4 stores
c input: all, output: ppart, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
c    .5*(q/m)*fx(x(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
c    .5*(q/m)*fy(x(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
c    .5*(q/m)*fz(x(t))*dt)
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
c omy = (q/m)*by(x(t))*gami and omz = (q/m)*bz(x(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equation used is:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
c where n = nearest grid point and dx = x-n
c similarly for fy(x), fz(x), by(x), bz(x)
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = momentum px of particle n in tile m
c ppart(3,n,m) = momentum py of particle n in tile m
c ppart(4,n,m) = momentum pz of particle n in tile m
c fxyz(1,j) = x component of force/charge at grid (j)
c fxyz(2,j) = y component of force/charge at grid (j)
c fxyz(3,j) = z component of force/charge at grid (j)
c that is, convolution of electric field over particle shape
c byz(1,j) = y component of magnetic field at grid (j)
c byz(2,j) = z component of magnetic field at grid (j)
c that is, the convolution of magnetic field over particle shape
c kpic = number of particles per tile
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx = system length in x direction
c mx = number of grids in sorting cell in x
c nxv = second dimension of field arrays, must be >= nx+1
c mx1 = (system length in x direction - 1)/mx + 1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ipbc
      real omx, qbm, dt, dtc, ci, ek
      real ppart, fxyz, byz
      integer kpic
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension kpic(mx1)
c local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real qtmh, ci2, edgelx, edgerx, dxp, amx, x, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real p2, gami, qtmg, dtg
      real sfxyz, sbyz
      dimension sfxyz(3,MXV), sbyz(2,MXV)
c     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum2 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c error if local array is too small
c     if (mx.ge.MXV) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy,acz, 
!$OMP& omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,
!$OMP& rot9,p2,gami,qtmg,dtg,sum1,sfxyz,sbyz)
!$OMP& REDUCTION(+:sum2)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
c load local fields from global array
      nn = min(mx,nx-noff) + 1
      do 10 j = 1, nn
      sfxyz(1,j) = fxyz(1,j+noff)
      sfxyz(2,j) = fxyz(2,j+noff)
      sfxyz(3,j) = fxyz(3,j+noff)
   10 continue
      do 20 j = 1, nn
      sbyz(1,j) = byz(1,j+noff)
      sbyz(2,j) = byz(2,j+noff)
   20 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 30 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
c find electric field
      dx = amx*sfxyz(1,nn) + dxp*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + dxp*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + dxp*sfxyz(3,nn+1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(2,j,k) + dx
      acy = ppart(3,j,k) + dy
      acz = ppart(4,j,k) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = omx
      oy = amx*sbyz(1,nn) + dxp*sbyz(1,nn+1)
      oz = amx*sbyz(2,nn) + dxp*sbyz(2,nn+1)
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
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(2,j,k) = dx
      ppart(3,j,k) = dy
      ppart(4,j,k) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = x + dx*dtg
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(2,j,k) = -ppart(2,j,k)
         endif
      endif
c set new position
      ppart(1,j,k) = dx
   30 continue
      sum2 = sum2 + sum1
   40 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine GRBPPUSHF13L(ppart,fxyz,byz,kpic,ncl,ihole,omx,qbm,dt, 
     1dtc,ci,ek,idimp,nppmx,nx,mx,nxv,mx1,ntmax,irc)
c for 1-2/2d code, this subroutine updates particle co-ordinate and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover, with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 90 flops/particle, 4 divides, 2 sqrts, 14 loads, 4 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
c    .5*(q/m)*fx(x(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
c    .5*(q/m)*fy(x(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) +
c    .5*(q/m)*fz(x(t))*dt)
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
c omy = (q/m)*by(x(t))*gami and omz = (q/m)*bz(x(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equation used is:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
c where n = nearest grid point and dx = x-n
c similarly for fy(x), fz(x), by(x), bz(x)
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = momentum px of particle n in tile m
c ppart(3,n,m) = momentum py of particle n in tile m
c ppart(4,n,m) = momentum pz of particle n in tile m
c fxyz(1,j) = x component of force/charge at grid (j)
c fxyz(2,j) = y component of force/charge at grid (j)
c fxyz(3,j) = z component of force/charge at grid (j)
c that is, convolution of electric field over particle shape
c byz(1,j) = y component of magnetic field at grid (j)
c byz(2,j) = z component of magnetic field at grid (j)
c that is, the convolution of magnetic field over particle shape
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx = system length in x direction
c mx = number of grids in sorting cell in x
c nxv = second dimension of field arrays, must be >= nx+1
c mx1 = (system length in x direction - 1)/mx + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer idimp, nppmx, nx, mx, nxv, mx1, ntmax, irc
      real omx, qbm, dt, dtc, ci, ek
      real ppart, fxyz, byz
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
c local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real qtmh, ci2, anx, edgelx, edgerx, dxp, amx, x, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real p2, gami, qtmg, dtg
      real sfxyz, sbyz
      dimension sfxyz(3,MXV), sbyz(2,MXV)
c     dimension sfxyz(3,mx+1), sbyz(2,mx+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      anx = real(nx)
      sum2 = 0.0d0
c error if local array is too small
c     if (mx.ge.MXV) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,dy,dz,ox,oy,oz,acx,acy
!$OMP& ,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,
!$OMP& rot8,rot9,p2,gami,qtmg,dtg,edgelx,edgerx,sum1,sfxyz,sbyz)
!$OMP& REDUCTION(+:sum2)
      do 50 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
c load local fields from global array
      nn = nn + 1
      do 10 j = 1, nn
      sfxyz(1,j) = fxyz(1,j+noff)
      sfxyz(2,j) = fxyz(2,j+noff)
      sfxyz(3,j) = fxyz(3,j+noff)
   10 continue
      do 20 j = 1, nn
      sbyz(1,j) = byz(1,j+noff)
      sbyz(2,j) = byz(2,j+noff)
   20 continue
c clear counters
      do 30 j = 1, 2
      ncl(j,k) = 0
   30 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 40 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = x - real(nn)
      nn = nn - noff + 1
      amx = 1.0 - dxp
c find electric field
      dx = amx*sfxyz(1,nn) + dxp*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + dxp*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + dxp*sfxyz(3,nn+1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(2,j,k) + dx
      acy = ppart(3,j,k) + dy
      acz = ppart(4,j,k) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = omx
      oy = amx*sbyz(1,nn) + dxp*sbyz(1,nn+1)
      oz = amx*sbyz(2,nn) + dxp*sbyz(2,nn+1)
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
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(2,j,k) = dx
      ppart(3,j,k) = dy
      ppart(4,j,k) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = x + dx*dtg
c find particles going out of bounds
      nn = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c nn = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         nn = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               nn = 1
            else
               dx = 0.0
            endif
         else
            nn = 1
         endif
      endif
c set new position
      ppart(1,j,k) = dx
c increment counters
      if (nn.gt.0) then
         ncl(nn,k) = ncl(nn,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = nn
         else
            nh = 1
         endif
      endif
   40 continue
      sum2 = sum2 + sum1
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine GJPPOST1L(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,mx,nxv,mx1
     1,ipbc)
c for 1-2/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c OpenMP version using guard cells
c data deposited in tiles
c particles stored segmented array
c 14 flops/particle, 8 loads, 5 stores
c input: all, output: ppart, cu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vi, where i = y,z
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = velocity vx of particle n in tile m
c ppart(3,n,m) = velocity vy of particle n in tile m
c ppart(4,n,m) = velocity vz of particle n in tile m
c cu(i,j) = ith component of current density at grid point j
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 4
c nx = system length in x direction
c mx = number of grids in sorting cell in x
c nxv = second dimension of current array, must be >= nx+1
c mx1 = (system length in x direction - 1)/mx + 1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer nppmx, idimp, nx, mx, nxv, mx1, ipbc
      real qm, dt
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mx1), cu(2,nxv)
      dimension kpic(mx1)
c local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real edgelx, edgerx, dxp, amx, vy, vz, x, dx
      real scu
      dimension scu(2,MXV)
c     dimension scu(2,mx+1)
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c error if local array is too small
c     if (mx.ge.MXV) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,vy,vz,scu)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
c zero out local accumulator
      do 10 j = 1, mx+1
      scu(1,j) = 0.0
      scu(2,j) = 0.0
   10 continue
c loop over particles in tile
      do 20 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = qm*(x - real(nn))
      nn = nn - noff + 1
      amx = qm - dxp
c deposit current within tile to local accumulator
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      scu(1,nn) = scu(1,nn) + vy*amx
      scu(2,nn) = scu(2,nn) + vz*amx
      scu(1,nn+1) = scu(1,nn+1) + vy*dxp
      scu(2,nn+1) = scu(2,nn+1) + vz*dxp
c advance position half a time-step
      dx = x + ppart(2,j,k)*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(2,j,k) = -ppart(2,j,k)
         endif
      endif
c set new position
      ppart(1,j,k) = dx
   20 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      do 30 j = 2, nn
      cu(1,j+noff) = cu(1,j+noff) + scu(1,j)
      cu(2,j+noff) = cu(2,j+noff) + scu(2,j)
   30 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      cu(1,1+noff) = cu(1,1+noff) + scu(1,1)
!$OMP ATOMIC
      cu(2,1+noff) = cu(2,1+noff) + scu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff) = cu(1,nn+noff) + scu(1,nn)
!$OMP ATOMIC
         cu(2,nn+noff) = cu(2,nn+noff) + scu(2,nn)
      endif
   40 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GJPPOSTF1L(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,idimp,nx
     1,mx,nxv,mx1,ntmax,irc)
c for 1-2/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c with periodic boundary conditions
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells
c data deposited in tiles
c particles stored segmented array
c 14 flops/particle, 8 loads, 5 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, cu, irc
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vi, where i = y,z
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = velocity vx of particle n in tile m
c ppart(3,n,m) = velocity vy of particle n in tile m
c ppart(4,n,m) = velocity vz of particle n in tile m
c cu(i,j) = ith component of current density at grid point j
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 4
c nx = system length in x direction
c mx = number of grids in sorting cell in x
c nxv = second dimension of current array, must be >= nx+1
c mx1 = (system length in x direction - 1)/mx + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer nppmx, idimp, nx, mx, nxv, mx1, ntmax, irc
      real qm, dt
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1), cu(2,nxv)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
c local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real anx, edgelx, edgerx, dxp, amx, vy, vz, x, dx
      real scu
      dimension scu(2,MXV)
c     dimension scu(2,mx+1)
      anx = real(nx)
c error if local array is too small
c     if (mx.ge.MXV) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,vy,vz,edgelx,edgerx,
!$OMP& scu)
      do 50 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
c zero out local accumulator
      do 10 j = 1, mx+1
      scu(1,j) = 0.0
      scu(2,j) = 0.0
   10 continue
c clear counters
      do 20 j = 1, 2
      ncl(j,k) = 0
   20 continue
c loop over particles in tile
      do 30 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = qm*(x - real(nn))
      nn = nn - noff + 1
      amx = qm - dxp
c deposit current within tile to local accumulator
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      scu(1,nn) = scu(1,nn) + vy*amx
      scu(2,nn) = scu(2,nn) + vz*amx
      scu(1,nn+1) = scu(1,nn+1) + vy*dxp
      scu(2,nn+1) = scu(2,nn+1) + vz*dxp
c advance position half a time-step
      dx = x + ppart(2,j,k)*dt
c find particles going out of bounds
      nn = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c nn = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         nn = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               nn = 1
            else
               dx = 0.0
            endif
         else
            nn = 1
         endif
      endif
c set new position
      ppart(1,j,k) = dx
c increment counters
      if (nn.gt.0) then
         ncl(nn,k) = ncl(nn,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = nn
         else
            nh = 1
         endif
      endif
   30 continue
c deposit current  to interior points in global array
      nn = min(mx,nxv-noff)
      do 40 j = 2, nn
      cu(1,j+noff) = cu(1,j+noff) + scu(1,j)
      cu(2,j+noff) = cu(2,j+noff) + scu(2,j)
   40 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      cu(1,1+noff) = cu(1,1+noff) + scu(1,1)
!$OMP ATOMIC
      cu(2,1+noff) = cu(2,1+noff) + scu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff) = cu(1,nn+noff) + scu(1,nn)
!$OMP ATOMIC
         cu(2,nn+noff) = cu(2,nn+noff) + scu(2,nn)
      endif
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GRJPPOST1L(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx,mx,nxv
     1,mx1,ipbc)
c for 1-2/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c OpenMP version using guard cells
c data deposited in tiles
c particles stored segmented array
c 24 flops/particle, 1 divide, 1 sqrt, 9 loads, 5 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*pi*gami, where i = y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = x momentum of particle n in tile m
c ppart(3,n,m) = y momentum of particle n in tile m
c ppart(4,n,m) = z momentum of particle n in tile m
c cu(i,j) = ith component of current density at grid point j
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 4
c nx = system length in x direction
c mx = number of grids in sorting cell in x
c nxv = second dimension of current array, must be >= nx+1
c mx1 = (system length in x direction - 1)/mx + 1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer nppmx, idimp, nx, mx, nxv, mx1, ipbc
      real qm, dt, ci
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mx1), cu(2,nxv)
      dimension kpic(mx1)
c local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real edgelx, ci2, edgerx, dxp, amx, vx, vy, vz, x, dx, p2, gami
      real scu
      dimension scu(2,MXV)
c     dimension scu(2,mx+1)
      ci2 = ci*ci
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,x,dxp,amx,dx,vx,vy,vz,p2,gami,scu)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
c zero out local accumulator
      do 10 j = 1, mx+1
      scu(1,j) = 0.0
      scu(2,j) = 0.0
   10 continue
c loop over particles in tile
      do 20 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = qm*(x - real(nn))
      nn = nn - noff + 1
      amx = qm - dxp
c find inverse gamma
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c deposit current within tile to local accumulator
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      scu(1,nn) = scu(1,nn) + vy*amx
      scu(2,nn) = scu(2,nn) + vz*amx
      scu(1,nn+1) = scu(1,nn+1) + vy*dxp
      scu(2,nn+1) = scu(2,nn+1) + vz*dxp
c advance position half a time-step
      dx = x + vx*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(2,j,k) = -ppart(2,j,k)
         endif
      endif
c set new position
      ppart(1,j,k) = dx
   20 continue
c deposit current  to interior points in global array
      nn = min(mx,nxv-noff)
      do 30 j = 2, nn
      cu(1,j+noff) = cu(1,j+noff) + scu(1,j)
      cu(2,j+noff) = cu(2,j+noff) + scu(2,j)
   30 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      cu(1,1+noff) = cu(1,1+noff) + scu(1,1)
!$OMP ATOMIC
      cu(2,1+noff) = cu(2,1+noff) + scu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff) = cu(1,nn+noff) + scu(1,nn)
!$OMP ATOMIC
         cu(2,nn+noff) = cu(2,nn+noff) + scu(2,nn)
      endif
   40 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GRJPPOSTF1L(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx,    
     1idimp,nx,mx,nxv,mx1,ntmax,irc)
c for 1-2/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells
c data deposited in tiles
c particles stored segmented array
c 24 flops/particle, 1 divide, 1 sqrt, 9 loads, 5 stores
c input: all except ncl, ihole, irc, output: part, ncl, ihole, cu, irc
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*pi*gami, where i = y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = x momentum of particle n in tile m
c ppart(3,n,m) = y momentum of particle n in tile m
c ppart(4,n,m) = z momentum of particle n in tile m
c cu(i,j) = ith component of current density at grid point j
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 4
c nx = system length in x direction
c mx = number of grids in sorting cell in x
c nxv = second dimension of current array, must be >= nx+1
c mx1 = (system length in x direction - 1)/mx + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer nppmx, idimp, nx, mx, nxv, mx1, ntmax, irc
      real qm, dt, ci
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1), cu(2,nxv)
      dimension kpic(mx1), ncl(2,mx1), ihole(2,ntmax+1,mx1)
c local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, ih, nh, nn
      real edgelx, ci2, edgerx, dxp, amx, vx, vy, vz, x, dx, p2, gami
      real anx
      real scu
      dimension scu(2,MXV)
c     dimension scu(2,mx+1)
      ci2 = ci*ci
      anx = real(nx)
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,x,dxp,amx,dx,vx,vy,vz,p2,gami,     
!$OMP& edgelx,edgerx,scu)
      do 50 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      edgelx = noff
      edgerx = noff + nn
      ih = 0
      nh = 0
c zero out local accumulator
      do 10 j = 1, mx+1
      scu(1,j) = 0.0
      scu(2,j) = 0.0
   10 continue
c clear counters
      do 20 j = 1, 2
      ncl(j,k) = 0
   20 continue
c loop over particles in tile
      do 30 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dxp = qm*(x - real(nn))
      nn = nn - noff + 1
      amx = qm - dxp
c find inverse gamma
      vx = ppart(2,j,k)
      vy = ppart(3,j,k)
      vz = ppart(4,j,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c deposit current within tile to local accumulator
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      scu(1,nn) = scu(1,nn) + vy*amx
      scu(2,nn) = scu(2,nn) + vz*amx
      scu(1,nn+1) = scu(1,nn+1) + vy*dxp
      scu(2,nn+1) = scu(2,nn+1) + vz*dxp
c advance position half a time-step
      dx = x + vx*dt
c find particles going out of bounds
      nn = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c nn = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         nn = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               nn = 1
            else
               dx = 0.0
            endif
         else
            nn = 1
         endif
      endif
c set new position
      ppart(1,j,k) = dx
c increment counters
      if (nn.gt.0) then
         ncl(nn,k) = ncl(nn,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = nn
         else
            nh = 1
         endif
      endif
   30 continue
c deposit current  to interior points in global array
      nn = min(mx,nxv-noff)
      do 40 j = 2, nn
      cu(1,j+noff) = cu(1,j+noff) + scu(1,j)
      cu(2,j+noff) = cu(2,j+noff) + scu(2,j)
   40 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      cu(1,1+noff) = cu(1,1+noff) + scu(1,1)
!$OMP ATOMIC
      cu(2,1+noff) = cu(2,1+noff) + scu(2,1)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff) = cu(1,nn+noff) + scu(1,nn)
!$OMP ATOMIC
         cu(2,nn+noff) = cu(2,nn+noff) + scu(2,nn)
      endif
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GPPOST1L(ppart,q,kpic,qm,nppmx,idimp,mx,nxv,mx1)
c for 1d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c OpenMP version using guard cells
c data deposited in tiles
c particles stored segmented array
c 7 flops/particle, 3 loads, 3 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n)=qm*(1.-dx) and q(n+1)=qm*dx
c where n = nearest grid point and dx = x-n
c ppart(1,n,m) = position x of particle n in tile m
c q(j) = charge density at grid point j
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 2
c mx = number of grids in sorting cell in x
c nxv = first dimension of charge array, must be >= nx+1
c mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer nppmx, idimp, mx, nxv, mx1
      real qm
      real ppart, q
      integer kpic
      dimension ppart(idimp,nppmx,mx1), q(nxv)
      dimension kpic(mx1)
c local data
      integer MXV
      parameter(MXV=129)
      integer noff, npp
      integer j, k, nn
      real x, dx
      real sq
c     dimension sq(MXV)
      dimension sq(mx+1)
c error if local array is too small
c     if (mx.ge.MXV) return
c loop over tiles
!$OMP PARALLEL DO PRIVATE(j,k,noff,npp,nn,x,dx,sq)
      do 40 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
c zero out local accumulator
      do 10 j = 1, mx+1
      sq(j) = 0.0
   10 continue
c loop over particles in tile
      do 20 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      nn = x
      dx = qm*(x - real(nn))
      nn = nn - noff + 1
c deposit charge within tile to local accumulator
      sq(nn) = sq(nn) + (qm - dx)
      sq(nn+1) = sq(nn+1) + dx
   20 continue
c deposit charge to interior points in global array
      nn = min(mx,nxv-noff)
      do 30 j = 2, nn
      q(j+noff) = q(j+noff) + sq(j)
   30 continue
      nn = min(mx+1,nxv-noff)
!$OMP ATOMIC
      q(1+noff) = q(1+noff) + sq(1)
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noff) = q(nn+noff) + sq(nn)
      endif
   40 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPORDER1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx,mx
     1,mx1,npbmx,ntmax,irc)
c this subroutine sorts particles by x grid in tiles of mx, 
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 1D linear memory
c algorithm has 3 steps.  first, one finds particles leaving tile and
c stores their number in each directon, location, and destination in ncl
c and ihole.  second, a prefix scan of ncl is performed and departing
c particles are buffered in ppbuff in direction order.  finally, we copy
c the incoming particles from other tiles into ppart.
c input: all except ppbuff, ncl, ihole, irc
c output: ppart, ppbuff, kpic, ncl, ihole, irc
c ppart(1,n,k) = position x of particle n in tile k
c ppbuff(i,n,k) = i co-ordinate of particle n in tile k
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = direction destination of particle leaving hole
c all for tile k
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c idimp = size of phase space = 2
c nppmx = maximum number of particles in tile
c nx = system length in x direction
c mx = number of grids in sorting cell in x
c mx1 = (system length in x direction - 1)/mx + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, mx, mx1, npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1), ppbuff(idimp,npbmx,mx1)
      dimension kpic(mx1), ncl(2,mx1)
      dimension ihole(2,ntmax+1,mx1)
c local data
      integer noff, npp, ncoff
      integer i, j, k, ii, ih, nh, ist, nn, isum
      integer ip, j1, j2, kxl, kxr
      real anx, edgelx, edgerx, dx
      integer ks
      dimension ks(2)
      anx = real(nx)
c find and count particles leaving tiles and determine destination
c update ppart, ihole, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,npp,nn,ih,nh,ist,dx,edgelx,edgerx)
      do 30 k = 1, mx1
      noff = mx*(k - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      ih = 0
      nh = 0
      edgelx = noff
      edgerx = noff + nn
c clear counters
      do 10 j = 1, 2
      ncl(j,k) = 0
   10 continue
c loop over particles in tile
      do 20 j = 1, npp
      dx = ppart(1,j,k)
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
      do 70 k = 1, mx1
c find address offset for ordered ppbuff array
      isum = 0
      do 40 j = 1, 2
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
      if (ip.gt.0) irc = ncl(2,k)
   70 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,npp,kxl,kxr,ih,nh,ncoff,ist,j1,j2,ip,ks)
      do 130 k = 1, mx1
      npp = kpic(k)
c loop over tiles in x, assume periodic boundary conditions
      kxl = k - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = k + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
c find tile number for different directions
      ks(1) = kxr
      ks(2) = kxl
c loop over directions
      nh = ihole(1,1,k)
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 100 ii = 1, 2
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
c ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
      do 90 j = 1, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,k)
c place overflow at end of array
      else
         j1 = npp + 1
         npp = j1
      endif
      if (j1.le.nppmx) then
         do 80 i = 1, idimp
         ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
   80    continue
      else
         ist = 1
      endif
   90 continue
  100 continue
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
      if (ih.lt.nh) then
         ip = nh - ih
         do 120 j = 1, ip
         j1 = npp - j + 1
         j2 = ihole(1,nh-j+2,k)
         if (j1.gt.j2) then
c move particle only if it is below current hole
            do 110 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
  110       continue
         endif
  120    continue
         npp = npp - ip
      endif
      kpic(k) = npp
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPORDERF1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,mx1,
     1npbmx,ntmax,irc)
c this subroutine sorts particles by x grid in tiles of mx, 
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 1D linear memory.
c the algorithm has 2 steps.  first, a prefix scan of ncl is performed
c and departing particles are buffered in ppbuff in direction order.
c then we copy the incoming particles from other tiles into ppart.
c it assumes that the number, location, and destination of particles 
c leaving a tile have been previously stored in ncl and ihole by the
c GPPUSHF1L subroutine.
c input: all except ppbuff, irc
c output: ppart, ppbuff, kpic, ncl, irc
c ppart(1,n,k) = position x of particle n in tile k
c ppbuff(i,n,k) = i co-ordinate of particle n in tile k
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = direction destination of particle leaving hole
c all for tile k
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c idimp = size of phase space = 2
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1), ppbuff(idimp,npbmx,mx1)
      dimension kpic(mx1), ncl(2,mx1)
      dimension ihole(2,ntmax+1,mx1)
c local data
      integer npp, ncoff
      integer i, j, k, ii, ih, nh, ist, isum
      integer ip, j1, j2, kxl, kxr
      integer ks
      dimension ks(2)
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 40 k = 1, mx1
c find address offset for ordered ppbuff array
      isum = 0
      do 10 j = 1, 2
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
      if (ip.gt.0) irc = ncl(2,k)
   40 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,npp,kxl,kxr,ih,nh,ncoff,ist,j1,j2,ip,ks)
      do 100 k = 1, mx1
      npp = kpic(k)
c loop over tiles in x, assume periodic boundary conditions
      kxl = k - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = k + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
c find tile number for different directions
      ks(1) = kxr
      ks(2) = kxl
c loop over directions
      nh = ihole(1,1,k)
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 70 ii = 1, 2
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
c ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
      do 60 j = 1, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,k)
c place overflow at end of array
      else
         j1 = npp + 1
         npp = j1
      endif
      if (j1.le.nppmx) then
         do 50 i = 1, idimp
         ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
   50    continue
      else
         ist = 1
      endif
   60 continue
   70 continue
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
      if (ih.lt.nh) then
         ip = nh - ih
         do 90 j = 1, ip
         j1 = npp - j + 1
         j2 = ihole(1,nh-j+2,k)
         if (j1.gt.j2) then
c move particle only if it is below current hole
            do 80 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
   80       continue
         endif
   90    continue
         npp = npp - ip
      endif
      kpic(k) = npp
  100 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine CGUARD1L(byz,nx,nxe)
c replicate extended periodic vector field byz
c linear interpolation
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      real byz
      integer nx, nxe
      dimension byz(2,nxe)
c local data
      integer i
c copy edges of extended field
      do 10 i = 1, 2
      byz(i,nx+1) = byz(i,1)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine BGUARD1L(fxyz,nx,nxe)
c replicate extended periodic vector field fxyz
c linear interpolation
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      real fxyz
      integer nx, nxe
      dimension fxyz(3,nxe)
c local data
      integer i
c copy edges of extended field
      do 10 i = 1, 3
      fxyz(i,nx+1) = fxyz(i,1)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine ACGUARD1L(cu,nx,nxe)
c accumulate extended periodic vector field cu
c linear interpolation
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      real cu
      integer nx, nxe
      dimension cu(2,nxe)
c local data
      integer i
c accumulate edges of extended field
      do 10 i = 1, 2
      cu(i,1) = cu(i,1) + cu(i,nx+1)
      cu(i,nx+1) = 0.0
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine AGUARD1L(q,nx,nxe)
c accumulate extended periodic scalar field q
c linear interpolation
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      real q
      integer nx, nxe
      dimension q(nxe)
c accumulate edges of extended field
      q(1) = q(1) + q(nx+1)
      q(nx+1) = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine POIS1(q,fx,isign,ffc,ax,affp,we,nx)
c this subroutine solves 1d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,affp,nx, output: ffc
c for isign  /= 0, input: q,ffc,isign,nx, output: fx,we
c approximate flop count is: 6*nx
c equation used is:
c fx(k) = -sqrt(-1)*k*g(k)*s(k)*q(k), where k = 2pi*j/nx, j=fourier mode,
c g(k) = (affp/k**2)*s(k), and s(k) = exp(-(k*ax)**2/2), except for
c fx(k=0) = fx(k=pi) = 0.
c cmplx(q(2*j-1),q(2*j)) = complex charge density for fourier mode j-1
c cmplx(fx(2*j-1),fx(2*j)) = complex force/charge for fourier mode j-1
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c ffc(2*j) = finite-size particle shape factor s for fourier mode j-1
c ffc(2*j-1) = potential green's function g for fourier mode j-1
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np = number of particles
c electric field energy is also calculated, using
c we = nx*sum((affp/k**2)*|q(k)*s(k)|**2)
c nx = system length in x direction
      implicit none
      integer isign, nx
      real ax, affp, we
      real q, fx, ffc
      dimension q(nx), fx(nx), ffc(nx)
c local data
      integer j, nxh
      real dnx, dkx, at1, at2
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      if (isign.ne.0) go to 20
c prepare form factor array
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      ffc(2*j) = exp(-.5*(dkx*ax)**2)
      ffc(2*j-1) = affp*ffc(2*j)/(dkx*dkx)
   10 continue
      ffc(1) = affp
      ffc(2) = 1.0
      return
c calculate force/charge and sum field energy
   20 wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 30 j = 2, nxh
      at1 = ffc(2*j-1)*ffc(2*j)
      at2 = dnx*real(j - 1)*at1
      fx(2*j-1) = at2*q(2*j)
      fx(2*j) = -at2*q(2*j-1)
      wp = wp + at1*(q(2*j-1)**2 + q(2*j)**2)
   30 continue
c mode number kx = 0
      fx(1) = 0.
      fx(2) = 0.
      we = real(nx)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine IBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
c this subroutine solves 1-2/2d poisson's equation in fourier space for
c magnetic field, with periodic boundary conditions.
c input: cu,ffc,ci,nx,nxv, output: byz,wm
c approximate flop count is: 29*nxc
c where nxc = nx/2 - 1
c the magnetic field is calculated using the equations:
c by(kx) = -ci*ci*sqrt(-1)*g(kx)*kx*cuz(kx),
c bz(kx) = ci*ci*sqrt(-1)*g(kx)*kx*cuy(kx),
c where kx = 2pi*j/nx, and j = fourier mode number,
c g(kx) = (affp/kx**2)*s(kx),
c s(kx) = exp(-((kx*ax)**2)/2), except for
c by(kx=pi) = bz(kx=pi) = 0, and by(kx=0) = bz(kx=0) = 0.
c cu(i,j) = complex current density for fourier mode (j-1)
c byz(i,j) = i component of complex magnetic field
c all for fourier mode (j-1)
c aimag(ffc(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffc(j)) = potential green's function g
c for fourier mode (j-1)
c ci = reciprocal of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2), where
c affp = normalization constant = nx/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
c nxhd = second dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, wm
      complex cu, byz, ffc
      dimension cu(2,nxvh), byz(2,nxvh), ffc(nxhd)
c local data
      integer j, nxh
      real dnx, ci2, at1, at2
      complex zero, zt1, zt2
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at1 = ci2*real(ffc(j))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffc(j))
      zt1 = cmplx(-aimag(cu(2,j)),real(cu(2,j)))
      zt2 = cmplx(-aimag(cu(1,j)),real(cu(1,j)))
      byz(1,j) = -at2*zt1
      byz(2,j) = at2*zt2
      wp = wp + at1*(cu(1,j)*conjg(cu(1,j)) + cu(2,j)*conjg(cu(2,j)))
   10 continue
      byz(1,1) = zero
      byz(2,1) = zero
      wm = real(nx)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine MAXWEL1(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions.
c input: all, output: wf, wm, eyz, byz
c approximate flop count is: 87*nxc
c where nxc = nx/2 - 1
c the magnetic field is first updated half a step using the equations:
c by(kx) = by(kx) + .5*dt*sqrt(-1)*kx*ez(kx)
c bz(kx) = bz(kx) - .5*dt*sqrt(-1)*kx*ey(kx)
c the electric field is then updated a whole step using the equations:
c ey(kx) = ey(kx) - c2*dt*sqrt(-1)*kx*bz(kx) - affp*dt*cuy(kx)*s(kx)
c ez(kx) = ez(kx) + c2*dt*sqrt(-1)*kx*by(kx) - affp*dt*cuz(kx)*s(kx)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, c2 = 1./(ci*ci)
c and s(kx) = exp(-((kx*ax)**2)
c j = fourier mode numbers, except for
c ey(kx=pi) = ez(kx=pi) = 0, and ey(kx=0) = ez(kx=0) = 0.
c and similarly for by, bz.
c cu(i,j) = complex current density
c eyz(i,j) = complex transverse electric field
c byz(i,j) = complex magnetic field
c for component i, all for fourier mode (j-1)
c real(ffc(1)) = affp = normalization constant = nx/np,
c where np=number of particles
c aimag(ffc(j)) = finite-size particle shape factor s,
c s(kx) = exp(-((kx*ax)**2)/2)
c for fourier mode (j-1)
c ci = reciprocal of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*sum((1/affp)*|eyz(kx)|**2)
c magnetic field energy is also calculated, using
c wm = nx*sum((c2/affp)*|byz(kx)|**2)
c nx = system length in x direction
c nxvh = first dimension of field arrays, must be >= nxh
c nxhd = first dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, dt, wf, wm
      complex eyz, byz, cu, ffc
      dimension eyz(2,nxvh), byz(2,nxvh), cu(2,nxvh), ffc(nxhd)
c local data
      integer j, nxh
      real dnx, dth, c2, cdt, affp, adt, anorm, dkx, afdt
      complex zero, zt1, zt2, zt5, zt6, zt8, zt9
      double precision wp, ws
      if (ci.le.0.0) return
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      dth = 0.5*dt
      c2 = 1.0/(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1))
      adt = affp*dt
      zero = cmplx(0.0,0.0)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j))
c update magnetic field half time step
      zt1 = cmplx(-aimag(eyz(2,j)),real(eyz(2,j)))
      zt2 = cmplx(-aimag(eyz(1,j)),real(eyz(1,j)))
      zt5 = byz(1,j) + dth*(dkx*zt1)
      zt6 = byz(2,j) - dth*(dkx*zt2)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt8 = eyz(1,j) - cdt*(dkx*zt1) - afdt*cu(1,j)
      zt9 = eyz(2,j) + cdt*(dkx*zt2) - afdt*cu(2,j)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      eyz(1,j) = zt8
      eyz(2,j) = zt9
      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      byz(1,j) = zt5
      byz(2,j) = zt6
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
   10 continue
c mode numbers kx = 0, nx/2
c     eyz(1,1) = zero
c     eyz(2,1) = zero
c     byz(1,1) = zero
c     byz(2,1) = zero
      wf = real(nx)*ws
      wm = real(nx)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine EMFIELD1(fxyz,fx,eyz,ffc,nx,nxvh,nxhd)
c this subroutine merges complex vector fields
c includes additional smoothing
      implicit none
      integer nx, nxvh, nxhd
      complex fxyz, fx, eyz, ffc
      dimension fxyz(3,nxvh), fx(nxvh), eyz(2,nxvh)
      dimension ffc(nxhd)
c local data
      integer j, nxh
      real at1
      nxh = nx/2
c add the fields
      do 10 j = 1, nxh
      at1 = aimag(ffc(j))
      fxyz(1,j) = fx(j)
      fxyz(2,j) = eyz(1,j)*at1
      fxyz(3,j) = eyz(2,j)*at1
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine BMFIELD1(fyz,eyz,ffc,nx,nxvh,nxhd)
c this subroutine copies complex vector fields
c includes additional smoothing
      implicit none
      integer nx, nxvh, nxhd
      complex fyz, eyz, ffc
      dimension fyz(2,nxvh), eyz(2,nxvh)
      dimension ffc(nxhd)
c local data
      integer j, nxh
      real at1
      nxh = nx/2
      do 10 j = 1, nxh
      at1 = aimag(ffc(j))
      fyz(1,j) = eyz(1,j)*at1
      fyz(2,j) = eyz(2,j)*at1
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT1RINIT(mixup,sct,indx,nxhd)
c this subroutine calculates tables needed by a one dimensional
c real to complex fast fourier transform and its inverse.
c input: indx, nxhd
c output: mixup, sct
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c indx = exponent which determines length in x direction,
c where nx=2**indx
c nxhd = nx/2
c written by viktor k. decyk, ucla
      implicit none
      integer indx, nxhd
      integer mixup
      complex sct
      dimension mixup(nxhd), sct(nxhd)
c local data
      integer indx1, nx, nxh
      integer j, k, lb, ll, jb, it
      real dnx, arg
      indx1 = indx - 1
      nx = 2**indx
      nxh = nx/2
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxh
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles 2*n*pi/nx
      dnx = 6.28318530717959/real(nx)
      do 30 j = 1, nxh
      arg = dnx*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
c this subroutine performs a one dimensional real to complex fast
c fourier transform and its inverse, using complex arithmetic
c for isign = (-1,1), input: all except t, output: f, t
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = nx/2
c f = input and output data
c t = complex scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = -1, an inverse fourier transform is performed
c f(n) = (1/nx)*sum(f(j)*exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, a forward fourier transform is performed
c f(j) = sum(f(n)*exp(sqrt(-1)*2pi*n*j/nx))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c fourier coefficients are stored as follows:
c f(j) = real,imaginary part of mode j-1, 1 <= j <= nx/2
c except for
c real(f(1)) = real part of mode 0,
c aimag(f(1)) = real part of mode nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sct
      dimension f(nxd), t(nxhd), mixup(nxhd), sct(nxhd)
c local data
      integer indx1, nx, nxh, nxhh, nxh2, j, k, l, j1
      integer nxs, nxs2, km, km2, k1, k2
      real ani
      complex t1, t2
      if (isign.eq.0) return
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      indx1 = indx - 1
      if (isign.gt.0) go to 70
c inverse fourier transform
c bit-reverse array elements to complex temporary
      do 10 j = 1, nxh
      j1 = mixup(j)
      t(j) = cmplx(f(2*j1-1),f(2*j1))
   10 continue
c transform
      do 40 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 30 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 20 j = 1, nxs
      t1 = sct(1+km2*(j-1))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
   20 continue
   30 continue
   40 continue
c unscramble coefficients and normalize
      ani = 1.0/real(2*nx)
      do 50 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sct(j)),-real(sct(j)))
      t(j) = ani*(t1 + t2)
      t(nxh2-j) = ani*conjg(t1 - t2)
   50 continue
      ani = 2.*ani
      t(nxhh+1) = ani*conjg(t(nxhh+1))
      t(1) = ani*cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(1))
     1)
c move to real destination
      do 60 j = 1, nxh
      f(2*j-1) = real(t(j))
      f(2*j) = aimag(t(j))
   60 continue
      return
c forward fourier transform
c move to complex temporary
   70 do 80 j = 1, nxh
      t(j) = cmplx(f(2*j-1),f(2*j))
   80 continue
c scramble coefficients
      do 90 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sct(j)),real(sct(j)))
      t(j) = t1 + t2
      t(nxh2-j) = conjg(t1 - t2)
   90 continue
      t(nxhh+1) = 2.*conjg(t(nxhh+1))
      t(1) = cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(1)))
c bit-reverse array elements to real destination
      do 100 j = 1, nxh
      j1 = mixup(j)
      f(2*j-1) = real(t(j1))
      f(2*j) = aimag(t(j1))
  100 continue
c move back to complex temporary
      do 110 j = 1, nxh
      t(j) = cmplx(f(2*j-1),f(2*j))
  110 continue
c transform
      do 140 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 130 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 120 j = 1, nxs
      t1 = conjg(sct(1+km2*(j-1)))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
  120 continue
  130 continue
  140 continue
c move to real destination
      do 150 j = 1, nxh
      f(2*j-1) = real(t(j))
      f(2*j) = aimag(t(j))
  150 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT1R2X(f,t,isign,mixup,sct,indx,nxd,nxhd)
c this subroutine performs two one dimensional real to complex fast
c fourier transforms and their inverses, using complex arithmetic
c for isign = (-1,1), input: all except t, output: f, t
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = nx/2
c f = input and output data
c t = complex scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c f(1:2,n) = (1/nx)*sum(f(1:2,j)*exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, a forward fourier transform is performed
c f(1:2,j) = sum(f(1:2,n)*exp(sqrt(-1)*2pi*n*j/nx))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c fourier coefficients are stored as follows:
c f(1:2,j) = real,imaginary part of mode j-1, 1 <= j <= nx/2
c except for
c real(f(1:2,1)) = real part of mode 0,
c aimag(f(1:2,1)) = real part of mode nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sct
      dimension f(2,nxd), t(2,nxhd), mixup(nxhd), sct(nxhd)
c local data
      integer indx1, nx, nxh, nxhh, nxh2, j, k, l, j1
      integer nxs, nxs2, km, km2, k1, k2, jj
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      indx1 = indx - 1
      if (isign.gt.0) go to 90
c inverse fourier transform
c bit-reverse array elements to complex temporary
      do 10 j = 1, nxh
      j1 = mixup(j)
      t(1,j) = cmplx(f(1,2*j1-1),f(1,2*j1))
      t(2,j) = cmplx(f(2,2*j1-1),f(2,2*j1))
   10 continue
c transform
      do 40 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 30 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 20 j = 1, nxs
      t1 = sct(1+km2*(j-1))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
   20 continue
   30 continue
   40 continue
c unscramble coefficients and normalize
      ani = 1.0/real(2*nx)
      do 60 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),-real(sct(j)))
      do 50 jj = 1, 2
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = ani*(t1 + t2)
      t(jj,nxh2-j) = ani*conjg(t1 - t2)
   50 continue
   60 continue
      ani = 2.*ani
      do 70 jj = 1, 2
      t(jj,nxhh+1) = ani*conjg(t(jj,nxhh+1))
      t(jj,1) = ani*cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)) -
     1 aimag(t(jj,1)))
   70 continue
c move to complex destination
      do  80 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(2,2*j-1) = aimag(t(1,j))
      f(1,2*j) = real(t(2,j))
      f(2,2*j) = aimag(t(2,j))
   80 continue
      return
c forward fourier transform
c move complex source to complex temporary
   90 do 100 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(2,2*j-1))
      t(2,j) = cmplx(f(1,2*j),f(2,2*j))
  100 continue
c scramble coefficients
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),real(sct(j)))
      do 110 jj = 1, 2
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = t1 + t2
      t(jj,nxh2-j) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 jj = 1, 2
      t(jj,nxhh+1) = 2.*conjg(t(jj,nxhh+1))
      t(jj,1) = cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)) - aim
     1ag(t(jj,1)))
  130 continue
c bit-reverse array elements to real destination
      do 140 j = 1, nxh
      j1 = mixup(j)
      f(1,2*j-1) = real(t(1,j1))
      f(1,2*j) = aimag(t(1,j1))
      f(2,2*j-1) = real(t(2,j1))
      f(2,2*j) = aimag(t(2,j1))
  140 continue
c move back to complex temporary
      do 150 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(1,2*j))
      t(2,j) = cmplx(f(2,2*j-1),f(2,2*j))
  150 continue
c transform
      do 180 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 170 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 160 j = 1, nxs
      t1 = conjg(sct(1+km2*(j-1)))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
  160 continue
  170 continue
  180 continue
c move to real destination
      do 190 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(1,2*j) = aimag(t(1,j))
      f(2,2*j-1) = real(t(2,j))
      f(2,2*j) = aimag(t(2,j))
  190 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT1R3X(f,t,isign,mixup,sct,indx,nxd,nxhd)
c this subroutine performs three one dimensional real to complex fast
c fourier transforms and their inverses, using complex arithmetic
c for isign = (-1,1), input: all except t, output: f, t
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = nx/2
c f = input and output data
c t = complex scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c f(1:3,n) = (1/nx)*sum(f(1:3,j)*exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, a forward fourier transform is performed
c f(1:3,j) = sum(f(1:3,n)*exp(sqrt(-1)*2pi*n*j/nx))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c fourier coefficients are stored as follows:
c f(1) = real part of mode 0, f(2) = real part of mode nx/2
c f(2*j-1),f(2*j) = real,imaginary part of mode j-1, 0 < j < nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sct
      dimension f(3,nxd), t(3,nxhd), mixup(nxhd), sct(nxhd)
c local data
      integer indx1, nx, nxh, nxhh, nxh2, j, k, l, j1
      integer nxs, nxs2, km, km2, k1, k2, jj
      real ani
      complex t1, t2, t3, t4
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      indx1 = indx - 1
      if (isign.gt.0) go to 90
c inverse fourier transform
c bit-reverse array elements to complex temporary
      do 10 j = 1, nxh
      j1 = mixup(j)
      t(1,j) = cmplx(f(1,2*j1-1),f(1,2*j1))
      t(2,j) = cmplx(f(2,2*j1-1),f(2,2*j1))
      t(3,j) = cmplx(f(3,2*j1-1),f(3,2*j1))
   10 continue
c transform
      do 40 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 30 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 20 j = 1, nxs
      t1 = sct(1+km2*(j-1))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t4 = t1*t(3,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(3,j+k2) = t(3,j+k1) - t4
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
      t(3,j+k1) = t(3,j+k1) + t4
   20 continue
   30 continue
   40 continue
c unscramble coefficients and normalize result
      ani = 1./real(2*nx)
      do 60 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),-real(sct(j)))
      do 50 jj = 1, 3
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = ani*(t1 + t2)
      t(jj,nxh2-j) = ani*conjg(t1 - t2)
   50 continue
   60 continue
      ani = 2.0*ani
      do 70 jj = 1, 3
      t(jj,nxhh+1) = ani*conjg(t(jj,nxhh+1))
      t(jj,1) = ani*cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)) -
     1 aimag(t(jj,1)))
   70 continue
c move to complex destination
      do 80 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(2,2*j-1) = aimag(t(1,j))
      f(3,2*j-1) = real(t(2,j))
      f(1,2*j) = aimag(t(2,j))
      f(2,2*j) = real(t(3,j))
      f(3,2*j) = aimag(t(3,j))
   80 continue
      return
c forward fourier transform
c move complex source to complex temporary
   90 do 100 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(2,2*j-1))
      t(2,j) = cmplx(f(3,2*j-1),f(1,2*j))
      t(3,j) = cmplx(f(2,2*j),f(3,2*j))
  100 continue
c scramble coefficients
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),real(sct(j)))
      do 110 jj = 1, 3
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = t1 + t2
      t(jj,nxh2-j) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 jj = 1, 3
      t(jj,nxhh+1) = 2.*conjg(t(jj,nxhh+1))
      t(jj,1) = cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)) - aim
     1ag(t(jj,1)))
  130 continue
c bit-reverse array elements to real destination
      do 140 j = 1, nxh
      j1 = mixup(j)
      f(1,2*j-1) = real(t(1,j1))
      f(1,2*j) = aimag(t(1,j1))
      f(2,2*j-1) = real(t(2,j1))
      f(2,2*j) = aimag(t(2,j1))
      f(3,2*j-1) = real(t(3,j1))
      f(3,2*j) = aimag(t(3,j1))
  140 continue
c move back to complex temporary
      do 150 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(1,2*j))
      t(2,j) = cmplx(f(2,2*j-1),f(2,2*j))
      t(3,j) = cmplx(f(3,2*j-1),f(3,2*j))
  150 continue
c transform
      do 180 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 170 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 160 j = 1, nxs
      t1 = conjg(sct(1+km2*(j-1)))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t4 = t1*t(3,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(3,j+k2) = t(3,j+k1) - t4
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
      t(3,j+k1) = t(3,j+k1) + t4
  160 continue
  170 continue
  180 continue
c move to real destination
      do 190 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(1,2*j) = aimag(t(1,j))
      f(2,2*j-1) = real(t(2,j))
      f(2,2*j) = aimag(t(2,j))
      f(3,2*j-1) = real(t(3,j))
      f(3,2*j) = aimag(t(3,j))
  190 continue
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
      subroutine PPCOPYOUT(part,ppart,kpic,nop,nppmx,idimp,mx1,irc)
c for 1d code, this subroutine copies segmented particle data ppart to
c the array part with original tiled layout
c input: all except part, output: part
c part(i,j) = i-th coordinate for particle j
c ppart(i,j,k) = i-th coordinate for particle j in tile k
c kpic = number of particles per tile
c nop = number of particles
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 4
c mx1 = total number of tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nop, nppmx, idimp, mx1, irc
      real part, ppart
      integer kpic
      dimension part(idimp,nop), ppart(idimp,nppmx,mx1)
      dimension kpic(mx1)
c local data
      integer i, j, k, npoff, npp, ne, ierr
      npoff = 0
      ierr = 0
c loop over tiles
      do 30 k = 1, mx1
      npp = kpic(k)
      ne = npp + npoff
      if (ne.gt.nop) ierr = max(ierr,ne-nop)
      if (ierr.gt.0) npp = 0
c loop over particles in tile
      do 20 j = 1, npp
      do 10 i = 1, idimp
      part(i,j+npoff) = ppart(i,j,k)
   10 continue
   20 continue
      npoff = npoff + npp
   30 continue
      if (ierr.gt.0) irc = ierr
      return
      end
