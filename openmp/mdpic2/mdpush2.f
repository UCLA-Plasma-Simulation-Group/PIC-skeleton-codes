c Fortran Library for Skeleton 2-1/2D Darwin OpenMP PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp,nop,
     1nx,ny,ipbc)
c for 2-1/2d code, this subroutine calculates initial particle
c co-ordinates and velocities with uniform density and maxwellian
c velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, idimp, nop, nx, ny, ipbc
      real vtx, vty, vtz, vdx, vdy, vdz
      real part
      dimension part(idimp,nop)
c local data
      integer j, k, k1, npxy
      real edgelx, edgely, at1, at2, at3, sum1, sum2, sum3
      double precision dsum1, dsum2, dsum3
      double precision ranorm
      npxy = npx*npy
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
      k1 = npx*(k - 1)
      at3 = edgely + at2*(real(k) - 0.5)
      do 10 j = 1, npx
      part(1,j+k1) = edgelx + at1*(real(j) - 0.5)
      part(2,j+k1) = at3
   10 continue
   20 continue
c maxwellian velocity distribution
      do 30 j = 1, npxy
      part(3,j) = vtx*ranorm()
      part(4,j) = vty*ranorm()
      part(5,j) = vtz*ranorm()
   30 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 40 j = 1, npxy
      dsum1 = dsum1 + part(3,j)
      dsum2 = dsum2 + part(4,j)
      dsum3 = dsum3 + part(5,j)
   40 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1./real(npxy)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 50 j = 1, npxy
      part(3,j) = part(3,j) - sum1
      part(4,j) = part(4,j) - sum2
      part(5,j) = part(5,j) - sum3
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DBLKP2L(part,kpic,nppmx,idimp,nop,mx,my,mx1,mxy1,irc)
c this subroutine finds the maximum number of particles in each tile of
c mx, my to calculate size of segmented particle array ppart
c linear interpolation
c input: all except kpic, nppmx, output: kpic, nppmx
c part = input particle array
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c kpic = output number of particles per tile
c nppmx = return maximum number of particles in tile
c idimp = size of phase space = 4
c nop = number of particles
c mx/my = number of grids in sorting cell in x and y
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kpic, nppmx, idimp, nop, mx, my, mx1, mxy1, irc
      real part
      dimension part(idimp,nop), kpic(mxy1)
c local data
      integer j, k, n, m, isum, ist, npx, ierr
      ierr = 0
c clear counter array
      do 10 k = 1, mxy1
      kpic(k) = 0
   10 continue
c find how many particles in each tile
      do 20 j = 1, nop
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = m/my
      m = n + mx1*m
      if (m.le.mxy1) then
         kpic(m) = kpic(m) + 1
      else
         ierr = max(ierr,m-mxy1)
      endif
   20 continue
c find maximum
      isum = 0
      npx = 0
      do 30 k = 1, mxy1
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
      subroutine PPMOVIN2L(part,ppart,kpic,nppmx,idimp,nop,mx,my,mx1,   
     1mxy1,irc)
c this subroutine sorts particles by x,y grid in tiles of mx, my
c and copies to segmented array ppart
c linear interpolation
c input: all except ppart, kpic, output: ppart, kpic
c part/ppart = input/output particle arrays
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c ppart(3,n,m) = velocity vx of particle n in tile m
c ppart(4,n,m) = velocity vy of particle n in tile m
c kpic = output number of particles per tile
c nppmx = rmaximum number of particles in tile
c idimp = size of phase space = 4
c nop = number of particles
c mx/my = number of grids in sorting cell in x and y
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kpic, nppmx, idimp, nop, mx, my, mx1, mxy1, irc
      real part, ppart
      dimension part(idimp,nop), ppart(idimp,nppmx,mxy1)
      dimension kpic(mxy1)
c local data
      integer i, j, k, n, m, ip, ierr
      ierr = 0
c clear counter array
      do 10 k = 1, mxy1
      kpic(k) = 0
   10 continue
c find addresses of particles at each tile and reorder particles
      do 30 j = 1, nop
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = m/my
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
      subroutine PPCHECK2L(ppart,kpic,idimp,nppmx,nx,ny,mx,my,mx1,my1,  
     1irc)
c this subroutine performs a sanity check to make sure particles sorted
c by x,y grid in tiles of mx, my, are all within bounds.
c tiles are assumed to be arranged in 2D linear memory
c input: all except irc
c output: irc
c ppart(1,n,k) = position x of particle n in tile k
c ppart(2,n,k) = position y of particle n in tile k
c kpic(k) = number of reordered output particles in tile k
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, ny, mx, my, mx1, my1, irc
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mx1*my1)
      dimension kpic(mx1*my1)
c local data
      integer mxy1, noff, moff, npp, j, k, ist, nn, mm
      real edgelx, edgely, edgerx, edgery, dx, dy
      mxy1 = mx1*my1
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,moff,npp,nn,mm,ist,edgelx,edgely,edgerx,edgery,
!$OMP& dx,dy)
      do 20 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
c loop over particles in tile
      do 10 j = 1, npp
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
      subroutine GBPPUSH23L(ppart,fxy,bxy,kpic,qbm,dt,dtc,ek,idimp,nppmx
     1,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c OpenMP version using guard cells
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
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c ppart(3,n,m) = velocity vx of particle n in tile m
c ppart(4,n,m) = velocity vy of particle n in tile m
c ppart(5,n,m) = velocity vz of particle n in tile m
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c kpic = number of particles per tile
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
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ipbc
      real qbm, dt, dtc, ek
      real ppart, fxy, bxy
      integer kpic
      dimension ppart(idimp,nppmx,mxy1)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm
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
      edgely = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,x,y,dxp,dyp,amx,amy,dx,dy,dz,ox
!$OMP& ,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,
!$OMP& rot5,rot6,rot7,rot8,rot9,sum1,sfxy,sbxy)
!$OMP& REDUCTION(+:sum2)
      do 60 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c load local fields from global arrays
      nn = min(mx,nx-noff) + 1
      mm = min(my,ny-moff) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noff,j+moff)
      sfxy(2,i,j) = fxy(2,i+noff,j+moff)
      sfxy(3,i,j) = fxy(3,i+noff,j+moff)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noff,j+moff)
      sbxy(2,i,j) = bxy(2,i+noff,j+moff)
      sbxy(3,i,j) = bxy(3,i+noff,j+moff)
   30 continue
   40 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 50 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1
      mm = mm - moff + 1
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
      subroutine GBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,dtc,ek,
     1idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 119 flops/particle, 1 divide, 29 loads, 5 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, irc, ek
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
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c ppart(3,n,m) = velocity vx of particle n in tile m
c ppart(4,n,m) = velocity vy of particle n in tile m
c ppart(5,n,m) = velocity vz of particle n in tile m
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
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
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ntmax
      integer irc
      real qbm, dt, dtc, ek
      real ppart, fxy, bxy
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxy1)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension kpic(mxy1), ncl(8,mxy1)
      dimension ihole(2,ntmax+1,mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, ih, nh, nn, mm
      real qtmh, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
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
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,ih,nh,x,y,dxp,dyp,amx,amy,dx,dy
!$OMP& ,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,
!$OMP& rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,edgery,sum1,
!$OMP& sfxy,sbxy)
!$OMP& REDUCTION(+:sum2)
      do 70 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
      ih = 0
      nh = 0
c load local fields from global arrays
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxy(1,i,j) = fxy(1,i+noff,j+moff)
      sfxy(2,i,j) = fxy(2,i+noff,j+moff)
      sfxy(3,i,j) = fxy(3,i+noff,j+moff)
   10 continue
   20 continue
      do 40 j = 1, mm+1
      do 30 i = 1, nn+1
      sbxy(1,i,j) = bxy(1,i+noff,j+moff)
      sbxy(2,i,j) = bxy(2,i+noff,j+moff)
      sbxy(3,i,j) = bxy(3,i+noff,j+moff)
   30 continue
   40 continue
c clear counters
      do 50 j = 1, 8
      ncl(j,k) = 0
   50 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 60 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1
      mm = mm - moff + 1
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
      subroutine GPPOST2L(ppart,q,kpic,qm,nppmx,idimp,mx,my,nxv,nyv,mx1,
     1mxy1)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c OpenMP version using guard cells
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
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c q(j,k) = charge density at grid point j,k
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 4
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
      implicit none
      integer nppmx, idimp, mx, my, nxv, nyv, mx1, mxy1
      real qm
      real ppart, q
      integer kpic
      dimension ppart(idimp,nppmx,mxy1), q(nxv,nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm
      real x, y, dxp, dyp, amx, amy
      real sq
c     dimension sq(MXV,MYV)
      dimension sq(mx+1,my+1)
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,x,y,dxp,dyp,amx,amy,sq)
      do 80 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      sq(i,j) = 0.0
   10 continue
   20 continue
c loop over particles in tile
      do 30 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noff + 1
      mm = mm - moff + 1
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
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 50 j = 2, mm
      do 40 i = 2, nn
      q(i+noff,j+moff) = q(i+noff,j+moff) + sq(i,j)
   40 continue
   50 continue
c deposit charge to edge points in global array
      mm = min(my+1,nyv-moff)
      do 60 i = 2, nn
!$OMP ATOMIC
      q(i+noff,1+moff) = q(i+noff,1+moff) + sq(i,1)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noff,mm+moff) = q(i+noff,mm+moff) + sq(i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noff)
      do 70 j = 1, mm
!$OMP ATOMIC
      q(1+noff,j+moff) = q(1+noff,j+moff) + sq(1,j)
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noff,j+moff) = q(nn+noff,j+moff) + sq(nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GJPPOST2L(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny,mx,my, 
     1nxv,nyv,mx1,mxy1,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c OpenMP version using guard cells
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
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c ppart(3,n,m) = x velocity of particle n in tile m
c ppart(4,n,m) = y velocity of particle n in tile m
c ppart(5,n,m) = z velocity of particle n in tile m
c cu(i,j,k) = ith component of current density at grid point j,k
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of current array, must be >= nx+1
c nyv = second dimension of current array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ipbc
      real qm, dt
      real ppart, cu
      integer kpic
      dimension ppart(idimp,nppmx,mxy1), cu(3,nxv,nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real scu
      dimension scu(3,MXV,MYV)
c     dimension scu(3,mx+1,my+1)
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,x,y,dxp,dyp,amx,amy,dx,dy,vx,vy
!$OMP& ,vz,scu)
      do 80 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   10 continue
   20 continue
c loop over particles in tile
      do 30 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noff + 1
      mm = mm - moff + 1
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
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 50 j = 2, mm
      do 40 i = 2, nn
      cu(1,i+noff,j+moff) = cu(1,i+noff,j+moff) + scu(1,i,j)
      cu(2,i+noff,j+moff) = cu(2,i+noff,j+moff) + scu(2,i,j)
      cu(3,i+noff,j+moff) = cu(3,i+noff,j+moff) + scu(3,i,j)
   40 continue
   50 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 60 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noff,1+moff) = cu(1,i+noff,1+moff) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noff,1+moff) = cu(2,i+noff,1+moff) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noff,1+moff) = cu(3,i+noff,1+moff) + scu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noff,mm+moff) = cu(1,i+noff,mm+moff) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noff,mm+moff) = cu(2,i+noff,mm+moff) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noff,mm+moff) = cu(3,i+noff,mm+moff) + scu(3,i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noff)
      do 70 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noff,j+moff) = cu(1,1+noff,j+moff) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noff,j+moff) = cu(2,1+noff,j+moff) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noff,j+moff) = cu(3,1+noff,j+moff) + scu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff,j+moff) = cu(1,nn+noff,j+moff) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noff,j+moff) = cu(2,nn+noff,j+moff) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noff,j+moff) = cu(3,nn+noff,j+moff) + scu(3,nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GMJPPOST2L(ppart,amu,kpic,qm,nppmx,idimp,mx,my,nxv,nyv,
     1mx1,mxy1)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using first-order spline interpolation
c OpenMP version using guard cells
c data deposited in tiles
c particles stored segmented array
c 51 flops/particle, 21 loads, 16 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c ppart(1,n,m) = position x of particle n in tile m at t
c ppart(2,n,m) = position y of particle n in tile m at t
c ppart(3,n,m) = x velocity of particle n in tile m at t - dt/2
c ppart(4,n,m) = y velocity of particle n in tile m at t - dt/2
c ppart(5,n,m) = z velocity of particle n in tile m at t - dt/2
c amu(i,j,k) = ith component of momentum flux at grid point j,k
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c mx/my = number of grids in sorting cell in x/y
c nxv = second dimension of flux array, must be >= nx+1
c nyv = third dimension of flux array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
      implicit none
      integer nppmx, idimp, mx, my, nxv, nyv, mx1, mxy1
      real qm
      real ppart, amu
      integer kpic
      dimension ppart(idimp,nppmx,mxy1), amu(4,nxv,nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm
      real dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, v1, v2, v3, v4
      real samu
      dimension samu(4,MXV,MYV)
c     dimension samu(4,mx+1,my+1)
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,x,y,dxp,dyp,amx,amy,dx,dy,vx,vy
!$OMP& ,vz,v1,v2,v3,v4,samu)
      do 80 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      samu(1,i,j) = 0.0
      samu(2,i,j) = 0.0
      samu(3,i,j) = 0.0
      samu(4,i,j) = 0.0
   10 continue
   20 continue
c loop over particles in tile
      do 30 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noff + 1
      mm = mm - moff + 1
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit momentum flux
      dx = amx*amy
      dy = dxp*amy
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      samu(1,nn,mm) = samu(1,nn,mm) + v1*dx
      samu(2,nn,mm) = samu(2,nn,mm) + v2*dx
      samu(3,nn,mm) = samu(3,nn,mm) + v3*dx
      samu(4,nn,mm) = samu(4,nn,mm) + v4*dx
      dx = amx*dyp
      samu(1,nn+1,mm) = samu(1,nn+1,mm) + v1*dy
      samu(2,nn+1,mm) = samu(2,nn+1,mm) + v2*dy
      samu(3,nn+1,mm) = samu(3,nn+1,mm) + v3*dy
      samu(4,nn+1,mm) = samu(4,nn+1,mm) + v4*dy
      dy = dxp*dyp
      samu(1,nn,mm+1) = samu(1,nn,mm+1) + v1*dx
      samu(2,nn,mm+1) = samu(2,nn,mm+1) + v2*dx
      samu(3,nn,mm+1) = samu(3,nn,mm+1) + v3*dx
      samu(4,nn,mm+1) = samu(4,nn,mm+1) + v4*dx
      samu(1,nn+1,mm+1) = samu(1,nn+1,mm+1) + v1*dy
      samu(2,nn+1,mm+1) = samu(2,nn+1,mm+1) + v2*dy
      samu(3,nn+1,mm+1) = samu(3,nn+1,mm+1) + v3*dy
      samu(4,nn+1,mm+1) = samu(4,nn+1,mm+1) + v4*dy
   30 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 50 j = 2, mm
      do 40 i = 2, nn
      amu(1,i+noff,j+moff) = amu(1,i+noff,j+moff) + samu(1,i,j)
      amu(2,i+noff,j+moff) = amu(2,i+noff,j+moff) + samu(2,i,j)
      amu(3,i+noff,j+moff) = amu(3,i+noff,j+moff) + samu(3,i,j)
      amu(4,i+noff,j+moff) = amu(4,i+noff,j+moff) + samu(4,i,j)
   40 continue
   50 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 60 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noff,1+moff) = amu(1,i+noff,1+moff) + samu(1,i,1)
!$OMP ATOMIC
      amu(2,i+noff,1+moff) = amu(2,i+noff,1+moff) + samu(2,i,1)
!$OMP ATOMIC
      amu(3,i+noff,1+moff) = amu(3,i+noff,1+moff) + samu(3,i,1)
!$OMP ATOMIC
      amu(4,i+noff,1+moff) = amu(4,i+noff,1+moff) + samu(4,i,1)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noff,mm+moff) = amu(1,i+noff,mm+moff) + samu(1,i,mm)
!$OMP ATOMIC
         amu(2,i+noff,mm+moff) = amu(2,i+noff,mm+moff) + samu(2,i,mm)
!$OMP ATOMIC
         amu(3,i+noff,mm+moff) = amu(3,i+noff,mm+moff) + samu(3,i,mm)
!$OMP ATOMIC
         amu(4,i+noff,mm+moff) = amu(4,i+noff,mm+moff) + samu(4,i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noff)
      do 70 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noff,j+moff) = amu(1,1+noff,j+moff) + samu(1,1,j)
!$OMP ATOMIC
      amu(2,1+noff,j+moff) = amu(2,1+noff,j+moff) + samu(2,1,j)
!$OMP ATOMIC
      amu(3,1+noff,j+moff) = amu(3,1+noff,j+moff) + samu(3,1,j)
!$OMP ATOMIC
      amu(4,1+noff,j+moff) = amu(4,1+noff,j+moff) + samu(4,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noff,j+moff) = amu(1,nn+noff,j+moff) + samu(1,nn,j)
!$OMP ATOMIC
         amu(2,nn+noff,j+moff) = amu(2,nn+noff,j+moff) + samu(2,nn,j)
!$OMP ATOMIC
         amu(3,nn+noff,j+moff) = amu(3,nn+noff,j+moff) + samu(3,nn,j)
!$OMP ATOMIC
         amu(4,nn+noff,j+moff) = amu(4,nn+noff,j+moff) + samu(4,nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GDJPPOST2L(ppart,fxy,bxy,dcu,amu,kpic,qm,qbm,dt,idimp, 
     1nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c and acceleration density using first-order spline interpolation.
c OpenMP version using guard cells
c data read/written in tiles
c particles stored segmented array
c 194 flops/particle, 1 divide, 57 loads, 28 stores
c input: all, output: dcu, amu
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
c dcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
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
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c ppart(1,n,m) = position x of particle n in tile m at t
c ppart(2,n,m) = position y of particle n in tile m at t
c ppart(3,n,m) = velocity vx of particle n in tile m at t - dt/2
c ppart(4,n,m) = velocity vy of particle n in tile m at t - dt/2
c ppart(5,n,m) = velocity vz of particle n in tile m at t - dt/2
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c dcu(i,j,k) = ith component of acceleration density
c at grid point j,k for i = 1, 3
c amu(i,j,k) = ith component of momentum flux
c at grid point j,k for i = 1, 4
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 5
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
      implicit none
      integer idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
      real qm, qbm, dt
      real ppart, fxy, bxy, dcu, amu
      integer kpic
      dimension ppart(idimp,nppmx,mxy1)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension dcu(3,nxv,nyv), amu(4,nxv,nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm
      real qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, v1, v2, v3, v4
      real sfxy, sbxy, sdcu, samu
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
      dimension sdcu(3,MXV,MYV), samu(4,MXV,MYV)
c     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
c     dimension sdcu(3,mx+1,my+1), samu(4,mx+1,my+1)
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,x,y,vx,vy,vz,v1,v2,v3,v4,dxp,  
!$OMP& dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,    
!$OMP& anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sfxy,sbxy,sdcu
!$OMP& ,samu)
      do 120 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c load local fields from global arrays
      nn = min(mx,nx-noff) + 1
      mm = min(my,ny-moff) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noff,j+moff)
      sfxy(2,i,j) = fxy(2,i+noff,j+moff)
      sfxy(3,i,j) = fxy(3,i+noff,j+moff)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noff,j+moff)
      sbxy(2,i,j) = bxy(2,i+noff,j+moff)
      sbxy(3,i,j) = bxy(3,i+noff,j+moff)
   30 continue
   40 continue
c zero out local accumulators
      do 60 j = 1, my+1
      do 50 i = 1, mx+1
      samu(1,i,j) = 0.0
      samu(2,i,j) = 0.0
      samu(3,i,j) = 0.0
      samu(4,i,j) = 0.0
      sdcu(1,i,j) = 0.0
      sdcu(2,i,j) = 0.0
      sdcu(3,i,j) = 0.0
   50 continue
   60 continue
c loop over particles in tile
      do 70 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1
      mm = mm - moff + 1
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
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
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
c deposit momentum flux and acceleration density
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = amx*amy
      dy = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      samu(1,nn,mm) = samu(1,nn,mm) + v1*dx
      samu(2,nn,mm) = samu(2,nn,mm) + v2*dx
      samu(3,nn,mm) = samu(3,nn,mm) + v3*dx
      samu(4,nn,mm) = samu(4,nn,mm) + v4*dx
      sdcu(1,nn,mm) = sdcu(1,nn,mm) + vx*dx
      sdcu(2,nn,mm) = sdcu(2,nn,mm) + vy*dx
      sdcu(3,nn,mm) = sdcu(3,nn,mm) + vz*dx
      dx = amx*dyp
      samu(1,nn+1,mm) = samu(1,nn+1,mm) + v1*dy
      samu(2,nn+1,mm) = samu(2,nn+1,mm) + v2*dy
      samu(3,nn+1,mm) = samu(3,nn+1,mm) + v3*dy
      samu(4,nn+1,mm) = samu(4,nn+1,mm) + v4*dy
      sdcu(1,nn+1,mm) = sdcu(1,nn+1,mm) + vx*dy
      sdcu(2,nn+1,mm) = sdcu(2,nn+1,mm) + vy*dy
      sdcu(3,nn+1,mm) = sdcu(3,nn+1,mm) + vz*dy
      dy = dxp*dyp
      samu(1,nn,mm+1) = samu(1,nn,mm+1) + v1*dx
      samu(2,nn,mm+1) = samu(2,nn,mm+1) + v2*dx
      samu(3,nn,mm+1) = samu(3,nn,mm+1) + v3*dx
      samu(4,nn,mm+1) = samu(4,nn,mm+1) + v4*dx
      sdcu(1,nn,mm+1) = sdcu(1,nn,mm+1) + vx*dx
      sdcu(2,nn,mm+1) = sdcu(2,nn,mm+1) + vy*dx
      sdcu(3,nn,mm+1) = sdcu(3,nn,mm+1) + vz*dx
      samu(1,nn+1,mm+1) = samu(1,nn+1,mm+1) + v1*dy
      samu(2,nn+1,mm+1) = samu(2,nn+1,mm+1) + v2*dy
      samu(3,nn+1,mm+1) = samu(3,nn+1,mm+1) + v3*dy
      samu(4,nn+1,mm+1) = samu(4,nn+1,mm+1) + v4*dy
      sdcu(1,nn+1,mm+1) = sdcu(1,nn+1,mm+1) + vx*dy
      sdcu(2,nn+1,mm+1) = sdcu(2,nn+1,mm+1) + vy*dy
      sdcu(3,nn+1,mm+1) = sdcu(3,nn+1,mm+1) + vz*dy
   70 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 90 j = 2, mm
      do 80 i = 2, nn
      amu(1,i+noff,j+moff) = amu(1,i+noff,j+moff) + samu(1,i,j)
      amu(2,i+noff,j+moff) = amu(2,i+noff,j+moff) + samu(2,i,j)
      amu(3,i+noff,j+moff) = amu(3,i+noff,j+moff) + samu(3,i,j)
      amu(4,i+noff,j+moff) = amu(4,i+noff,j+moff) + samu(4,i,j)
      dcu(1,i+noff,j+moff) = dcu(1,i+noff,j+moff) + sdcu(1,i,j)
      dcu(2,i+noff,j+moff) = dcu(2,i+noff,j+moff) + sdcu(2,i,j)
      dcu(3,i+noff,j+moff) = dcu(3,i+noff,j+moff) + sdcu(3,i,j)
   80 continue
   90 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 100 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noff,1+moff) = amu(1,i+noff,1+moff) + samu(1,i,1)
!$OMP ATOMIC
      amu(2,i+noff,1+moff) = amu(2,i+noff,1+moff) + samu(2,i,1)
!$OMP ATOMIC
      amu(3,i+noff,1+moff) = amu(3,i+noff,1+moff) + samu(3,i,1)
!$OMP ATOMIC
      amu(4,i+noff,1+moff) = amu(4,i+noff,1+moff) + samu(4,i,1)
!$OMP ATOMIC
      dcu(1,i+noff,1+moff) = dcu(1,i+noff,1+moff) + sdcu(1,i,1)
!$OMP ATOMIC
      dcu(2,i+noff,1+moff) = dcu(2,i+noff,1+moff) + sdcu(2,i,1)
!$OMP ATOMIC
      dcu(3,i+noff,1+moff) = dcu(3,i+noff,1+moff) + sdcu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noff,mm+moff) = amu(1,i+noff,mm+moff) + samu(1,i,mm)
!$OMP ATOMIC
         amu(2,i+noff,mm+moff) = amu(2,i+noff,mm+moff) + samu(2,i,mm)
!$OMP ATOMIC
         amu(3,i+noff,mm+moff) = amu(3,i+noff,mm+moff) + samu(3,i,mm)
!$OMP ATOMIC
         amu(4,i+noff,mm+moff) = amu(4,i+noff,mm+moff) + samu(4,i,mm)
!$OMP ATOMIC
         dcu(1,i+noff,mm+moff) = dcu(1,i+noff,mm+moff) + sdcu(1,i,mm)
!$OMP ATOMIC
         dcu(2,i+noff,mm+moff) = dcu(2,i+noff,mm+moff) + sdcu(2,i,mm)
!$OMP ATOMIC
         dcu(3,i+noff,mm+moff) = dcu(3,i+noff,mm+moff) + sdcu(3,i,mm)
      endif
  100 continue
      nn = min(mx+1,nxv-noff)
      do 110 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noff,j+moff) = amu(1,1+noff,j+moff) + samu(1,1,j)
!$OMP ATOMIC
      amu(2,1+noff,j+moff) = amu(2,1+noff,j+moff) + samu(2,1,j)
!$OMP ATOMIC
      amu(3,1+noff,j+moff) = amu(3,1+noff,j+moff) + samu(3,1,j)
!$OMP ATOMIC
      amu(4,1+noff,j+moff) = amu(4,1+noff,j+moff) + samu(4,1,j)
!$OMP ATOMIC
      dcu(1,1+noff,j+moff) = dcu(1,1+noff,j+moff) + sdcu(1,1,j)
!$OMP ATOMIC
      dcu(2,1+noff,j+moff) = dcu(2,1+noff,j+moff) + sdcu(2,1,j)
!$OMP ATOMIC
      dcu(3,1+noff,j+moff) = dcu(3,1+noff,j+moff) + sdcu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noff,j+moff) = amu(1,nn+noff,j+moff) + samu(1,nn,j)
!$OMP ATOMIC
         amu(2,nn+noff,j+moff) = amu(2,nn+noff,j+moff) + samu(2,nn,j)
!$OMP ATOMIC
         amu(3,nn+noff,j+moff) = amu(3,nn+noff,j+moff) + samu(3,nn,j)
!$OMP ATOMIC
         amu(4,nn+noff,j+moff) = amu(4,nn+noff,j+moff) + samu(4,nn,j)
!$OMP ATOMIC
         dcu(1,nn+noff,j+moff) = dcu(1,nn+noff,j+moff) + sdcu(1,nn,j)
!$OMP ATOMIC
         dcu(2,nn+noff,j+moff) = dcu(2,nn+noff,j+moff) + sdcu(2,nn,j)
!$OMP ATOMIC
         dcu(3,nn+noff,j+moff) = dcu(3,nn+noff,j+moff) + sdcu(3,nn,j)
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GDCJPPOST2L(ppart,fxy,bxy,cu,dcu,amu,kpic,qm,qbm,dt,   
     1idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density and current density using first-order spline
c interpolation.
c OpenMP version using guard cells
c data read/written in tiles
c particles stored segmented array
c 218 flops/particle, 1 divide, 69 loads, 40 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
c dcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
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
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c ppart(1,n,m) = position x of particle n in tile m at t
c ppart(2,n,m) = position y of particle n in tile m at t
c ppart(3,n,m) = velocity vx of particle n in tile m at t - dt/2
c ppart(4,n,m) = velocity vy of particle n in tile m at t - dt/2
c ppart(5,n,m) = velocity vz of particle n in tile m at t - dt/2
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,j,k) = ith component of current density
c at grid point j,k for i = 1, 3
c dcu(i,j,k) = ith component of acceleration density
c at grid point j,k for i = 1, 3
c amu(i,j,k) = ith component of momentum flux
c at grid point j,k for i = 1, 4
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 5
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
      implicit none
      integer idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1
      real qm, qbm, dt
      real ppart, fxy, bxy, cu, dcu, amu
      integer kpic
      dimension ppart(idimp,nppmx,mxy1)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension cu(3,nxv,nyv), dcu(3,nxv,nyv), amu(4,nxv,nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm
      real qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, v1, v2, v3, v4
      real sfxy, sbxy, scu, sdcu, samu
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
      dimension scu(3,MXV,MYV), sdcu(3,MXV,MYV), samu(4,MXV,MYV)
c     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
c     dimension scu(3,mx+1,my+1,), sdcu(3,mx+1,my+1), samu(4,mx+1,my+1)
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,x,y,vx,vy,vz,v1,v2,v3,v4,dxp,  
!$OMP& dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,    
!$OMP& anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sfxy,sbxy,scu,
!$OMP& sdcu,samu)
      do 120 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c load local fields from global arrays
      nn = min(mx,nx-noff) + 1
      mm = min(my,ny-moff) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noff,j+moff)
      sfxy(2,i,j) = fxy(2,i+noff,j+moff)
      sfxy(3,i,j) = fxy(3,i+noff,j+moff)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noff,j+moff)
      sbxy(2,i,j) = bxy(2,i+noff,j+moff)
      sbxy(3,i,j) = bxy(3,i+noff,j+moff)
   30 continue
   40 continue
c zero out local accumulators
      do 60 j = 1, my+1
      do 50 i = 1, mx+1
      samu(1,i,j) = 0.0
      samu(2,i,j) = 0.0
      samu(3,i,j) = 0.0
      samu(4,i,j) = 0.0
      sdcu(1,i,j) = 0.0
      sdcu(2,i,j) = 0.0
      sdcu(3,i,j) = 0.0
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   50 continue
   60 continue
c loop over particles in tile
      do 70 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1
      mm = mm - moff + 1
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
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      vz = ppart(5,j,k)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
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
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = amx*amy
      dy = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      samu(1,nn,mm) = samu(1,nn,mm) + v1*dx
      samu(2,nn,mm) = samu(2,nn,mm) + v2*dx
      samu(3,nn,mm) = samu(3,nn,mm) + v3*dx
      samu(4,nn,mm) = samu(4,nn,mm) + v4*dx
      sdcu(1,nn,mm) = sdcu(1,nn,mm) + vx*dx
      sdcu(2,nn,mm) = sdcu(2,nn,mm) + vy*dx
      sdcu(3,nn,mm) = sdcu(3,nn,mm) + vz*dx
      scu(1,nn,mm) = scu(1,nn,mm) + ox*dx
      scu(2,nn,mm) = scu(2,nn,mm) + oy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + oz*dx
      dx = amx*dyp
      samu(1,nn+1,mm) = samu(1,nn+1,mm) + v1*dy
      samu(2,nn+1,mm) = samu(2,nn+1,mm) + v2*dy
      samu(3,nn+1,mm) = samu(3,nn+1,mm) + v3*dy
      samu(4,nn+1,mm) = samu(4,nn+1,mm) + v4*dy
      sdcu(1,nn+1,mm) = sdcu(1,nn+1,mm) + vx*dy
      sdcu(2,nn+1,mm) = sdcu(2,nn+1,mm) + vy*dy
      sdcu(3,nn+1,mm) = sdcu(3,nn+1,mm) + vz*dy
      scu(1,nn+1,mm) = scu(1,nn+1,mm) + ox*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + oy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + oz*dy
      dy = dxp*dyp
      samu(1,nn,mm+1) = samu(1,nn,mm+1) + v1*dx
      samu(2,nn,mm+1) = samu(2,nn,mm+1) + v2*dx
      samu(3,nn,mm+1) = samu(3,nn,mm+1) + v3*dx
      samu(4,nn,mm+1) = samu(4,nn,mm+1) + v4*dx
      sdcu(1,nn,mm+1) = sdcu(1,nn,mm+1) + vx*dx
      sdcu(2,nn,mm+1) = sdcu(2,nn,mm+1) + vy*dx
      sdcu(3,nn,mm+1) = sdcu(3,nn,mm+1) + vz*dx
      scu(1,nn,mm+1) = scu(1,nn,mm+1) + ox*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + oy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + oz*dx
      samu(1,nn+1,mm+1) = samu(1,nn+1,mm+1) + v1*dy
      samu(2,nn+1,mm+1) = samu(2,nn+1,mm+1) + v2*dy
      samu(3,nn+1,mm+1) = samu(3,nn+1,mm+1) + v3*dy
      samu(4,nn+1,mm+1) = samu(4,nn+1,mm+1) + v4*dy
      sdcu(1,nn+1,mm+1) = sdcu(1,nn+1,mm+1) + vx*dy
      sdcu(2,nn+1,mm+1) = sdcu(2,nn+1,mm+1) + vy*dy
      sdcu(3,nn+1,mm+1) = sdcu(3,nn+1,mm+1) + vz*dy
      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + ox*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + oy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + oz*dy
   70 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 90 j = 2, mm
      do 80 i = 2, nn
      amu(1,i+noff,j+moff) = amu(1,i+noff,j+moff) + samu(1,i,j)
      amu(2,i+noff,j+moff) = amu(2,i+noff,j+moff) + samu(2,i,j)
      amu(3,i+noff,j+moff) = amu(3,i+noff,j+moff) + samu(3,i,j)
      amu(4,i+noff,j+moff) = amu(4,i+noff,j+moff) + samu(4,i,j)
      dcu(1,i+noff,j+moff) = dcu(1,i+noff,j+moff) + sdcu(1,i,j)
      dcu(2,i+noff,j+moff) = dcu(2,i+noff,j+moff) + sdcu(2,i,j)
      dcu(3,i+noff,j+moff) = dcu(3,i+noff,j+moff) + sdcu(3,i,j)
      cu(1,i+noff,j+moff) = cu(1,i+noff,j+moff) + scu(1,i,j)
      cu(2,i+noff,j+moff) = cu(2,i+noff,j+moff) + scu(2,i,j)
      cu(3,i+noff,j+moff) = cu(3,i+noff,j+moff) + scu(3,i,j)
   80 continue
   90 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 100 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noff,1+moff) = amu(1,i+noff,1+moff) + samu(1,i,1)
!$OMP ATOMIC
      amu(2,i+noff,1+moff) = amu(2,i+noff,1+moff) + samu(2,i,1)
!$OMP ATOMIC
      amu(3,i+noff,1+moff) = amu(3,i+noff,1+moff) + samu(3,i,1)
!$OMP ATOMIC
      amu(4,i+noff,1+moff) = amu(4,i+noff,1+moff) + samu(4,i,1)
!$OMP ATOMIC
      dcu(1,i+noff,1+moff) = dcu(1,i+noff,1+moff) + sdcu(1,i,1)
!$OMP ATOMIC
      dcu(2,i+noff,1+moff) = dcu(2,i+noff,1+moff) + sdcu(2,i,1)
!$OMP ATOMIC
      dcu(3,i+noff,1+moff) = dcu(3,i+noff,1+moff) + sdcu(3,i,1)
!$OMP ATOMIC
      cu(1,i+noff,1+moff) = cu(1,i+noff,1+moff) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noff,1+moff) = cu(2,i+noff,1+moff) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noff,1+moff) = cu(3,i+noff,1+moff) + scu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noff,mm+moff) = amu(1,i+noff,mm+moff) + samu(1,i,mm)
!$OMP ATOMIC
         amu(2,i+noff,mm+moff) = amu(2,i+noff,mm+moff) + samu(2,i,mm)
!$OMP ATOMIC
         amu(3,i+noff,mm+moff) = amu(3,i+noff,mm+moff) + samu(3,i,mm)
!$OMP ATOMIC
         amu(4,i+noff,mm+moff) = amu(4,i+noff,mm+moff) + samu(4,i,mm)
!$OMP ATOMIC
         dcu(1,i+noff,mm+moff) = dcu(1,i+noff,mm+moff) + sdcu(1,i,mm)
!$OMP ATOMIC
         dcu(2,i+noff,mm+moff) = dcu(2,i+noff,mm+moff) + sdcu(2,i,mm)
!$OMP ATOMIC
         dcu(3,i+noff,mm+moff) = dcu(3,i+noff,mm+moff) + sdcu(3,i,mm)
!$OMP ATOMIC
         cu(1,i+noff,mm+moff) = cu(1,i+noff,mm+moff) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noff,mm+moff) = cu(2,i+noff,mm+moff) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noff,mm+moff) = cu(3,i+noff,mm+moff) + scu(3,i,mm)
      endif
  100 continue
      nn = min(mx+1,nxv-noff)
      do 110 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noff,j+moff) = amu(1,1+noff,j+moff) + samu(1,1,j)
!$OMP ATOMIC
      amu(2,1+noff,j+moff) = amu(2,1+noff,j+moff) + samu(2,1,j)
!$OMP ATOMIC
      amu(3,1+noff,j+moff) = amu(3,1+noff,j+moff) + samu(3,1,j)
!$OMP ATOMIC
      amu(4,1+noff,j+moff) = amu(4,1+noff,j+moff) + samu(4,1,j)
!$OMP ATOMIC
      dcu(1,1+noff,j+moff) = dcu(1,1+noff,j+moff) + sdcu(1,1,j)
!$OMP ATOMIC
      dcu(2,1+noff,j+moff) = dcu(2,1+noff,j+moff) + sdcu(2,1,j)
!$OMP ATOMIC
      dcu(3,1+noff,j+moff) = dcu(3,1+noff,j+moff) + sdcu(3,1,j)
!$OMP ATOMIC
      cu(1,1+noff,j+moff) = cu(1,1+noff,j+moff) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noff,j+moff) = cu(2,1+noff,j+moff) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noff,j+moff) = cu(3,1+noff,j+moff) + scu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noff,j+moff) = amu(1,nn+noff,j+moff) + samu(1,nn,j)
!$OMP ATOMIC
         amu(2,nn+noff,j+moff) = amu(2,nn+noff,j+moff) + samu(2,nn,j)
!$OMP ATOMIC
         amu(3,nn+noff,j+moff) = amu(3,nn+noff,j+moff) + samu(3,nn,j)
!$OMP ATOMIC
         amu(4,nn+noff,j+moff) = amu(4,nn+noff,j+moff) + samu(4,nn,j)
!$OMP ATOMIC
         dcu(1,nn+noff,j+moff) = dcu(1,nn+noff,j+moff) + sdcu(1,nn,j)
!$OMP ATOMIC
         dcu(2,nn+noff,j+moff) = dcu(2,nn+noff,j+moff) + sdcu(2,nn,j)
!$OMP ATOMIC
         dcu(3,nn+noff,j+moff) = dcu(3,nn+noff,j+moff) + sdcu(3,nn,j)
!$OMP ATOMIC
         cu(1,nn+noff,j+moff) = cu(1,nn+noff,j+moff) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noff,j+moff) = cu(2,nn+noff,j+moff) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noff,j+moff) = cu(3,nn+noff,j+moff) + scu(3,nn,j)
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPORDER2L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx,ny
     1,mx,my,mx1,my1,npbmx,ntmax,irc)
c this subroutine sorts particles by x,y grid in tiles of mx, my
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 2D linear memory
c algorithm has 3 steps.  first, one finds particles leaving tile and
c stores their number in each directon, location, and destination in ncl
c and ihole.  second, a prefix scan of ncl is performed and departing
c particles are buffered in ppbuff in direction order.  finally, we copy
c the incoming particles from other tiles into ppart.
c input: all except ppbuff, ncl, ihole, irc
c output: ppart, ppbuff, kpic, ncl, ihole, irc
c ppart(1,n,k) = position x of particle n in tile k
c ppart(2,n,k) = position y of particle n in tile k
c ppbuff(i,n,k) = i co-ordinate of particle n in tile k
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = direction destination of particle leaving hole
c all for tile k
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, ny, mx, my, mx1, my1, npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1*my1), ppbuff(idimp,npbmx,mx1*my1)
      dimension kpic(mx1*my1), ncl(8,mx1*my1)
      dimension ihole(2,ntmax+1,mx1*my1)
c local data
      integer mxy1, noff, moff, npp, ncoff
      integer i, j, k, ii, kx, ky, ih, nh, ist, nn, mm, isum
      integer ip, j1, j2, kxl, kxr, kk, kl, kr
      real anx, any, edgelx, edgely, edgerx, edgery, dx, dy
      integer ks
      dimension ks(8)
      mxy1 = mx1*my1
      anx = real(nx)
      any = real(ny)
c find and count particles leaving tiles and determine destination
c update ppart, ihole, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,moff,npp,nn,mm,ih,nh,ist,dx,dy,edgelx,edgely,
!$OMP& edgerx,edgery)
      do 30 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      ih = 0
      nh = 0
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
c clear counters
      do 10 j = 1, 8
      ncl(j,k) = 0
   10 continue
c loop over particles in tile
      do 20 j = 1, npp
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
      do 70 k = 1, mxy1
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
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,kk,npp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,ist,j1,j2
!$OMP& ,ip,ks)
      do 130 k = 1, mxy1
      npp = kpic(k)
      ky = (k - 1)/mx1 + 1
c loop over tiles in y, assume periodic boundary conditions
      kk = (ky - 1)*mx1
c find tile above
      kl = ky - 1 
      if (kl.lt.1) kl = kl + my1
      kl = (kl - 1)*mx1
c find tile below
      kr = ky + 1
      if (kr.gt.my1) kr = kr - my1
      kr = (kr - 1)*mx1
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
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 100 ii = 1, 8
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
      subroutine PPORDERF2L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,mx1,
     1my1,npbmx,ntmax,irc)
c this subroutine sorts particles by x,y grid in tiles of mx, my
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 2D linear memory.
c the algorithm has 2 steps.  first, a prefix scan of ncl is performed
c and departing particles are buffered in ppbuff in direction order.
c then we copy the incoming particles from other tiles into ppart.
c it assumes that the number, location, and destination of particles 
c leaving a tile have been previously stored in ncl and ihole by the
c GPPUSHF2L subroutine.
c input: all except ppbuff, irc
c output: ppart, ppbuff, kpic, ncl, irc
c ppart(1,n,k) = position x of particle n in tile k
c ppart(2,n,k) = position y of particle n in tile k
c ppbuff(i,n,k) = i co-ordinate of particle n in tile k
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = direction destination of particle leaving hole
c all for tile k
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, my1, npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1*my1), ppbuff(idimp,npbmx,mx1*my1)
      dimension kpic(mx1*my1), ncl(8,mx1*my1)
      dimension ihole(2,ntmax+1,mx1*my1)
c local data
      integer mxy1, npp, ncoff
      integer i, j, k, ii, kx, ky, ih, nh, ist, isum
      integer ip, j1, j2, kxl, kxr, kk, kl, kr
      integer ks
      dimension ks(8)
      mxy1 = mx1*my1
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 40 k = 1, mxy1
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
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,kk,npp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,ist,j1,j2
!$OMP& ,ip,ks)
      do 100 k = 1, mxy1
      npp = kpic(k)
      ky = (k - 1)/mx1 + 1
c loop over tiles in y, assume periodic boundary conditions
      kk = (ky - 1)*mx1
c find tile above
      kl = ky - 1 
      if (kl.lt.1) kl = kl + my1
      kl = (kl - 1)*mx1
c find tile below
      kr = ky + 1
      if (kr.gt.my1) kr = kr - my1
      kr = (kr - 1)*mx1
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
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 70 ii = 1, 8
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
      subroutine BGUARD2L(bxy,nx,ny,nxe,nye)
c replicate extended periodic vector field bxy
c linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nxe = second dimension of field arrays, must be >= ny+1
      implicit none
      real bxy
      integer nx, ny, nxe, nye
      dimension bxy(3,nxe,nye)
c local data
      integer j, k
c copy edges of extended field
      do 10 k = 1, ny
      bxy(1,nx+1,k) = bxy(1,1,k)
      bxy(2,nx+1,k) = bxy(2,1,k)
      bxy(3,nx+1,k) = bxy(3,1,k)
   10 continue
      do 20 j = 1, nx
      bxy(1,j,ny+1) = bxy(1,j,1)
      bxy(2,j,ny+1) = bxy(2,j,1)
      bxy(3,j,ny+1) = bxy(3,j,1)
   20 continue
      bxy(1,nx+1,ny+1) = bxy(1,1,1)
      bxy(2,nx+1,ny+1) = bxy(2,1,1)
      bxy(3,nx+1,ny+1) = bxy(3,1,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine ACGUARD2L(cu,nx,ny,nxe,nye)
c accumulate extended periodic vector field cu
c linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nxe = second dimension of field arrays, must be >= ny+1
      implicit none
      real cu
      integer nx, ny, nxe, nye
      dimension cu(3,nxe,nye)
c local data
      integer i, j, k
c accumulate edges of extended field
      do 20 k = 1, ny
      do 10 i = 1, 3
      cu(i,1,k) = cu(i,1,k) + cu(i,nx+1,k)
      cu(i,nx+1,k) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, 3
      cu(i,j,1) = cu(i,j,1) + cu(i,j,ny+1)
      cu(i,j,ny+1) = 0.0
   30 continue
   40 continue
      do 50 i = 1, 3
      cu(i,1,1) = cu(i,1,1) + cu(i,nx+1,ny+1)
      cu(i,nx+1,ny+1) = 0.0
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine AGUARD2L(q,nx,ny,nxe,nye)
c accumulate extended periodic scalar field q
c linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nxe = second dimension of field arrays, must be >= ny+1
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k
c accumulate edges of extended field
      do 10 k = 1, ny
      q(1,k) = q(1,k) + q(nx+1,k)
      q(nx+1,k) = 0.0
   10 continue
      do 20 j = 1, nx
      q(j,1) = q(j,1) + q(j,ny+1)
      q(j,ny+1) = 0.0
   20 continue
      q(1,1) = q(1,1) + q(nx+1,ny+1)
      q(nx+1,ny+1) = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine AMCGUARD2L(amu,nx,ny,nxe,nye,ndim)
c accumulate extended periodic tensor field
c linear interpolation
      implicit none
      real amu
      integer nx, ny, nxe, nye, ndim
      dimension amu(ndim,nxe,nye)
c local data
      integer i, j, k
c accumulate edges of extended field
      do 20 k = 1, ny
      do 10 i = 1, ndim
      amu(i,1,k) = amu(i,1,k) + amu(i,nx+1,k)
      amu(i,nx+1,k) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, ndim
      amu(i,j,1) = amu(i,j,1) + amu(i,j,ny+1)
      amu(i,j,ny+1) = 0.0
   30 continue
   40 continue
      do 50 i = 1, ndim
      amu(i,1,1) = amu(i,1,1) + amu(i,nx+1,ny+1)
      amu(i,nx+1,ny+1) = 0.0
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine ASCFGUARD2L(dcu,cus,q2m0,nx,ny,nxe,nye)
c add scaled field to extended periodic field
      implicit none
      real dcu, cus, q2m0
      integer nx, ny, nxe, nye
      dimension dcu(3,nxe,nye), cus(3,nxe,nye)
c local data
      integer i, j, k
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do 30 k = 1, ny
      do 20 j = 1, nx
      do 10 i = 1, 3
      dcu(i,j,k) = dcu(i,j,k) - q2m0*cus(i,j,k)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine FWPMINMX2(qe,qbme,wpmax,wpmin,nx,ny,nxe,nye)
c calculates maximum and minimum plasma frequency.  assumes guard cells
c have already been added
c qe = charge density for electrons
c qbme = charge/mass ratio for electrons
c wpmax/wpmin = maximum/minimum plasma frequency
c nx/ny = system length in x/y direction
c nxe = first dimension of charge arrays, nxe must be >= nx
c nye = second dimension of charge arrays, nye must be >= ny
      implicit none
      real qe, qbme, wpmax, wpmin
      integer nx, ny, nxe, nye
      dimension qe(nxe,nye)
c local data
      integer j, k
      real at1
      wpmax = qbme*qe(1,1)
      wpmin = wpmax
!$OMP PARALLEL DO PRIVATE(j,k)
!$OMP& REDUCTION(max:wpmax), REDUCTION(min:wpmin)
      do 20 k = 1, ny
      do 10 j = 1, nx
      at1 = qbme*qe(j,k)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine MPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,  
     1nxhd,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.  Zeros out z component.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,nxvh,nyhd, output: fxy,we
c approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c fxy(3,j,k) = zero,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, we
      complex q, fxy, ffc
      dimension q(nxvh,nyv), fxy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*real(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*real(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.0) then
         ffc(j,k) = cmplx(affp,1.0)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 sum1 = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,at1,at2,at3,zt1,zt2,wp)
!$OMP& REDUCTION(+:sum1)
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      wp = 0.0d0
      do 40 j = 2, nxh
      at1 = real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(3,j,k) = zero
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
   40 continue
c mode numbers kx = 0, nx/2
      at1 = real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dny*real(k - 1)*at1
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at3*zt1
      fxy(3,1,k) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nxh
      at1 = real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*real(j - 1)*at1
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      fxy(1,j,1) = at2*zt1
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   60 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      sum1 = sum1 + wp
      we = real(nx*ny)*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine MCUPERP2(cu,nx,ny,nxvh,nyv)
c this subroutine calculates the transverse current in fourier space
c input: all, output: cu
c approximate flop count is: 36*nxc*nyc
c and nxc*nyc divides
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
c and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = first dimension of current array, must be >= nxh
c nyv = second dimension of current array, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex cu
      dimension cu(3,nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dkx, dky, dky2, at1
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c calculate transverse part of current
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,dky2,dkx,at1,zt1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      zt1 = at1*(dkx*cu(1,j,k) + dky*cu(2,j,k))
      cu(1,j,k) = cu(1,j,k) - dkx*zt1
      cu(2,j,k) = cu(2,j,k) - dky*zt1
      zt1 = at1*(dkx*cu(1,j,k1) - dky*cu(2,j,k1))
      cu(1,j,k1) = cu(1,j,k1) - dkx*zt1
      cu(2,j,k1) = cu(2,j,k1) + dky*zt1
   10 continue
c mode numbers kx = 0, nx/2
      cu(2,1,k) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      cu(1,j,1) = zero
      cu(1,j,k1) = zero
      cu(2,j,k1) = zero
   30 continue
      cu(1,1,1) = zero
      cu(2,1,1) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MBBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c with periodic boundary conditions.
c input: cu,ffc,ci,nx,ny,nxvh,nyhd, output: bxy,wm
c approximate flop count is: 90*nxc*nyc + 40*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c            s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = bx(ky=pi) = by(ky=pi) = bz(ky=pi) 
c = 0, and bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c bxy(1,j,k) = x component of complex magnetic field
c bxy(2,j,k) = y component of complex magnetic field
c bxy(3,j,k) = z component of complex magnetic field
c all for fourier mode (j-1,k-1)
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ci = reciprocal of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
c    |cu(kx,ky)*s(kx,ky)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real ci, wm
      complex cu, bxy, ffc
      dimension cu(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dky, ci2, at1, at2, at3
      complex zero, zt1, zt2, zt3
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
c calculate smoothed magnetic field and sum field energy
      sum1 = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,at1,at2,at3,zt1,zt2,zt3,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      wp = 0.0d0
      do 10 j = 2, nxh
      at1 = ci2*real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      bxy(1,j,k) = at3*zt1
      bxy(2,j,k) = -at2*zt1
      bxy(3,j,k) = at2*zt2 - at3*zt3
      zt1 = cmplx(-aimag(cu(3,j,k1)),real(cu(3,j,k1)))
      zt2 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
      zt3 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
      bxy(1,j,k1) = -at3*zt1
      bxy(2,j,k1) = -at2*zt1
      bxy(3,j,k1) = at2*zt2 + at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k))                         
     1   + cu(2,j,k)*conjg(cu(2,j,k)) + cu(3,j,k)*conjg(cu(3,j,k))      
     2   + cu(1,j,k1)*conjg(cu(1,j,k1)) + cu(2,j,k1)*conjg(cu(2,j,k1))  
     3   + cu(3,j,k1)*conjg(cu(3,j,k1)))
   10 continue
c mode numbers kx = 0, nx/2
      at1 = ci2*real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dny*real(k - 1)*at1
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bxy(1,1,k) = at3*zt1
      bxy(2,1,k) = zero
      bxy(3,1,k) = -at3*zt3
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k))                         
     1   + cu(2,1,k)*conjg(cu(2,1,k)) + cu(3,1,k)*conjg(cu(3,1,k)))
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      at1 = ci2*real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*real(j - 1)*at1
      zt1 = cmplx(-aimag(cu(3,j,1)),real(cu(3,j,1)))
      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
      bxy(1,j,1) = zero
      bxy(2,j,1) = -at2*zt1
      bxy(3,j,1) = at2*zt2
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1))                         
     1   + cu(2,j,1)*conjg(cu(2,j,1)) + cu(3,j,1)*conjg(cu(3,j,1)))
   30 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      sum1 = sum1 + wp
      wm = real(nx*ny)*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine BADDEXT2(bxy,omx,omy,omz,nx,ny,nxe,nye)
c adds constant to magnetic field for 2-1/2d code
c bxy = magnetic field
c omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
c nx/ny = system length in x/y direction
c nxe = second dimension of magnetic field array, nxe must be >= nx
c nye = third dimension of magnetic field array, nye must be >= ny
      implicit none
      real bxy, omx, omy, omz
      integer nx, ny, nxe, nye
      dimension bxy(3,nxe,nye)
c local data
      integer j, k
!$OMP PARALLEL DO PRIVATE(j,k)
      do 20 k = 1, ny
      do 10 j = 1, nx
      bxy(1,j,k) = bxy(1,j,k) + omx
      bxy(2,j,k) = bxy(2,j,k) + omy
      bxy(3,j,k) = bxy(3,j,k) + omz
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine MDCUPERP23(dcu,amu,nx,ny,nxvh,nyv)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux
c in 2-1/2d with periodic boundary conditions.
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = -sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = -sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c dcu(3,kx,ky) = -sqrt(-1)*(kx*vx*vz+ky*vy*vz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) =  dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on output:
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c amu(1,j,k) = xx-yy component of complex momentum flux
c amu(2,j,k) = xy component of complex momentum flux
c amu(3,j,k) = zx component of complex momentum flux
c amu(4,j,k) = zy component of complex momentum flux
c all for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex dcu, amu
      dimension dcu(3,nxvh,nyv), amu(4,nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,k1,dky,dky2,dkx,dkx2,dkxy,dkxy2,at1,zt1,zt2,zt3)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      dkx2 = dkx*dkx
      dkxy = dkx*dky
      dkxy2 = dky2 - dkx2
      at1 = 1.0/(dkx2 + dky2)
      zt1 = cmplx(aimag(amu(1,j,k)),-real(amu(1,j,k)))
      zt2 = cmplx(aimag(amu(2,j,k)),-real(amu(2,j,k)))
      zt3 = at1*(dkxy*zt1 + dkxy2*zt2)
      dcu(1,j,k) = dky*zt3
      dcu(2,j,k) = -dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k)),-real(amu(3,j,k)))
      zt2 = cmplx(aimag(amu(4,j,k)),-real(amu(4,j,k)))
      dcu(3,j,k) = dkx*zt1 + dky*zt2
      zt1 = cmplx(aimag(amu(1,j,k1)),-real(amu(1,j,k1)))
      zt2 = cmplx(aimag(amu(2,j,k1)),-real(amu(2,j,k1)))
      zt3 = at1*(dkxy*zt1 - dkxy2*zt2)
      dcu(1,j,k1) = dky*zt3
      dcu(2,j,k1) = dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k1)),-real(amu(3,j,k1)))
      zt2 = cmplx(aimag(amu(4,j,k1)),-real(amu(4,j,k1)))
      dcu(3,j,k1) = dkx*zt1 - dky*zt2
   10 continue
c mode numbers kx = 0, nx/2
      zt2 = cmplx(aimag(amu(2,1,k)),-real(amu(2,1,k)))
      dcu(1,1,k) = dky*zt2
      dcu(2,1,k) = zero
      zt2 = cmplx(aimag(amu(4,1,k)),-real(amu(4,1,k)))
      dcu(3,1,k) = dky*zt2
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1)),-real(amu(2,j,1)))
      dcu(1,j,1) = zero
      dcu(2,j,1) = dkx*zt2
      zt1 = cmplx(aimag(amu(3,j,1)),-real(amu(3,j,1)))
      dcu(3,j,1) = dkx*zt1
      dcu(1,j,k1) = zero
      dcu(2,j,k1) = zero
      dcu(3,j,k1) = zero
   30 continue
      dcu(1,1,1) = zero
      dcu(2,1,1) = zero
      dcu(3,1,1) = zero
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MADCUPERP23(dcu,amu,nx,ny,nxvh,nyv)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux and acceleration density
c in 2-1/2d with periodic boundary conditions.
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = dcu(1,kx,ky)-sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = dcu(2,kx,ky)-sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c dcu(3,kx,ky) = dcu(3,kx,ky)-sqrt(-1)*(kx*vx*vz+ky*vy*vz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) =  dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on input:
c dcu(i,j,k) = complex acceleration density for fourier mode (j-1,k-1)
c on output:
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c amu(1,j,k) = xx-yy component of complex momentum flux
c amu(2,j,k) = xy component of complex momentum flux
c amu(3,j,k) = zx component of complex momentum flux
c amu(4,j,k) = zy component of complex momentum flux
c all for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex dcu, amu
      dimension dcu(3,nxvh,nyv), amu(4,nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,k1,dky,dky2,dkx,dkx2,dkxy,dkxy2,at1,zt1,zt2,zt3)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      dkx2 = dkx*dkx
      dkxy = dkx*dky
      dkxy2 = dky2 - dkx2
      at1 = 1.0/(dkx2 + dky2)
      zt1 = cmplx(aimag(amu(1,j,k)),-real(amu(1,j,k)))
      zt2 = cmplx(aimag(amu(2,j,k)),-real(amu(2,j,k)))
      zt3 = at1*(dky*dcu(1,j,k) - dkx*dcu(2,j,k) + dkxy*zt1 + dkxy2*zt2)
      dcu(1,j,k) = dky*zt3
      dcu(2,j,k) = -dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k)),-real(amu(3,j,k)))
      zt2 = cmplx(aimag(amu(4,j,k)),-real(amu(4,j,k)))
      dcu(3,j,k) = dcu(3,j,k) + dkx*zt1 + dky*zt2
      zt1 = cmplx(aimag(amu(1,j,k1)),-real(amu(1,j,k1)))
      zt2 = cmplx(aimag(amu(2,j,k1)),-real(amu(2,j,k1)))
      zt3 = at1*(dky*dcu(1,j,k1) + dkx*dcu(2,j,k1) + dkxy*zt1 - dkxy2*zt
     12)
      dcu(1,j,k1) = dky*zt3
      dcu(2,j,k1) = dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k1)),-real(amu(3,j,k1)))
      zt2 = cmplx(aimag(amu(4,j,k1)),-real(amu(4,j,k1)))
      dcu(3,j,k1) = dcu(3,j,k1) + dkx*zt1 - dky*zt2
   10 continue
c mode numbers kx = 0, nx/2
      zt2 = cmplx(aimag(amu(2,1,k)),-real(amu(2,1,k)))
      dcu(1,1,k) = dcu(1,1,k) + dky*zt2
      dcu(2,1,k) = zero
      zt2 = cmplx(aimag(amu(4,1,k)),-real(amu(4,1,k)))
      dcu(3,1,k) = dcu(3,1,k) + dky*zt2
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1)),-real(amu(2,j,1)))
      dcu(1,j,1) = zero
      dcu(2,j,1) = dcu(2,j,1) + dkx*zt2
      zt1 = cmplx(aimag(amu(3,j,1)),-real(amu(3,j,1)))
      dcu(3,j,1) = dcu(3,j,1) + dkx*zt1
      dcu(1,j,k1) = zero
      dcu(2,j,k1) = zero
      dcu(3,j,k1) = zero
   30 continue
      dcu(1,1,1) = zero
      dcu(2,1,1) = zero
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MEPOIS23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny, 
     1nxvh,nyv,nxhd,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c transverse electric field (or convolution of transverse electric field
c over particle shape), with periodic boundary conditions.
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c for isign = 0, input: isign,ax,ay,affp,wp0,nx,ny,nxvh,nyhd, output:ffe
c for isign /= 0, input: dcu,ffe,isign,ci,nx,ny,nxvh,nyv,nxhd,nyhd,
c output: exy,wf
c approximate flop count is: 68*nxc*nyc + 33*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign = 0, form factor array is prepared
c if isign = -1, smoothed transverse electric field is calculated
c using the equation:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)*s(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2+wp0*ci2*s(kx,ky)**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = ex(ky=pi) = ey(ky=pi) = ez(ky=pi) 
c = 0, and ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c if isign = 1, unsmoothed transverse electric field is calculated
c using the equation:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)
c ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c exy(1,j,k) = x component of complex transverse electric field
c exy(2,j,k) = y component of complex transverse electric field
c exy(3,j,k) = z component of complex transverse electric field
c all for fourier mode (j-1,k-1)
c aimag(ffe(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffe(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c wp0 = normalized total plasma frequency squared
c ci = reciprocal of velocity of light
c transverse electric field energy is also calculated, using
c wf = nx*ny*sum((affp/((kx**2+ky**2)*ci*ci)**2)
c    |dcu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the derivative of current is
c divergence-free
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, wp0, ci, wf
      complex dcu, exy, ffe
      dimension dcu(3,nxvh,nyv), exy(3,nxvh,nyv)
      dimension ffe(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, ci2, wpc, dkx, dky, at1, at2, at3, at4
      complex zero
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
      wpc = wp0*ci2
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*real(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*real(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.0) then
         ffe(j,k) = cmplx(affp,1.0)
      else
         ffe(j,k) = cmplx(affp*at4/(at3 + wpc*at4*at4),at4)
      endif
   10 continue
   20 continue
      return
c calculate smoothed transverse electric field and sum field energy
   30 if (isign.gt.0) go to 70
      sum1 = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1,at2,wp)
!$OMP& REDUCTION(+:sum1)
      do 50 k = 2, nyh
      k1 = ny2 - k
      wp = 0.0d0
      do 40 j = 2, nxh
      at2 = -ci2*real(ffe(j,k))
      at1 = at2*aimag(ffe(j,k))
      at2 = at2*at2
      exy(1,j,k) = at1*dcu(1,j,k)
      exy(2,j,k) = at1*dcu(2,j,k)
      exy(3,j,k) = at1*dcu(3,j,k)
      exy(1,j,k1) = at1*dcu(1,j,k1)
      exy(2,j,k1) = at1*dcu(2,j,k1)
      exy(3,j,k1) = at1*dcu(3,j,k1)
      wp = wp + at2*(dcu(1,j,k)*conjg(dcu(1,j,k))                       
     1   + dcu(2,j,k)*conjg(dcu(2,j,k)) + dcu(3,j,k)*conjg(dcu(3,j,k))  
     2   + dcu(1,j,k1)*conjg(dcu(1,j,k1))                               
     3   + dcu(2,j,k1)*conjg(dcu(2,j,k1))                               
     4   + dcu(3,j,k1)*conjg(dcu(3,j,k1)))
   40 continue
c mode numbers kx = 0, nx/2
      at2 = -ci2*real(ffe(1,k))
      at1 = at2*aimag(ffe(1,k))
      at2 = at2*at2
      exy(1,1,k) = at1*dcu(1,1,k)
      exy(2,1,k) = at1*dcu(2,1,k)
      exy(3,1,k) = at1*dcu(3,1,k)
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wp = wp + at2*(dcu(1,1,k)*conjg(dcu(1,1,k))                       
     1   + dcu(2,1,k)*conjg(dcu(2,1,k)) + dcu(3,1,k)*conjg(dcu(3,1,k)))
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nxh
      at2 = -ci2*real(ffe(j,1))
      at1 = at2*aimag(ffe(j,1))
      at2 = at2*at2
      exy(1,j,1) = at1*dcu(1,j,1)
      exy(2,j,1) = at1*dcu(2,j,1)
      exy(3,j,1) = at1*dcu(3,j,1)
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
      wp = wp + at2*(dcu(1,j,1)*conjg(dcu(1,j,1))                       
     1   + dcu(2,j,1)*conjg(dcu(2,j,1)) + dcu(3,j,1)*conjg(dcu(3,j,1)))
   60 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      sum1 = sum1 + wp
      wf = real(nx*ny)*sum1/real(ffe(1,1))
      return
c calculate unsmoothed transverse electric field and sum field energy
   70 sum1 = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1,at2,wp)
!$OMP& REDUCTION(+:sum1)
      do 90 k = 2, nyh
      k1 = ny2 - k
      wp = 0.0d0
      do 80 j = 2, nxh
      at2 = -ci2*real(ffe(j,k))
      at1 = at2*at2
      exy(1,j,k) = at2*dcu(1,j,k)
      exy(2,j,k) = at2*dcu(2,j,k)
      exy(3,j,k) = at2*dcu(3,j,k)
      exy(1,j,k1) = at2*dcu(1,j,k1)
      exy(2,j,k1) = at2*dcu(2,j,k1)
      exy(3,j,k1) = at2*dcu(3,j,k1)
      wp = wp + at1*(dcu(1,j,k)*conjg(dcu(1,j,k))                       
     1   + dcu(2,j,k)*conjg(dcu(2,j,k)) + dcu(3,j,k)*conjg(dcu(3,j,k))  
     2   + dcu(1,j,k1)*conjg(dcu(1,j,k1))                               
     3   + dcu(2,j,k1)*conjg(dcu(2,j,k1))                               
     4   + dcu(3,j,k1)*conjg(dcu(3,j,k1)))
   80 continue
c mode numbers kx = 0, nx/2
      at2 = -ci2*real(ffe(1,k))
      at1 = at2*at2
      exy(1,1,k) = at2*dcu(1,1,k)
      exy(2,1,k) = at2*dcu(2,1,k)
      exy(3,1,k) = at2*dcu(3,1,k)
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wp = wp + at1*(dcu(1,1,k)*conjg(dcu(1,1,k))                       
     1   + dcu(2,1,k)*conjg(dcu(2,1,k)) + dcu(3,1,k)*conjg(dcu(3,1,k)))
      sum1 = sum1 + wp
   90 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 100 j = 2, nxh
      at2 = -ci2*real(ffe(j,1))
      at1 = at2*at2
      exy(1,j,1) = at2*dcu(1,j,1)
      exy(2,j,1) = at2*dcu(2,j,1)
      exy(3,j,1) = at2*dcu(3,j,1)
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
      wp = wp + at1*(dcu(1,j,1)*conjg(dcu(1,j,1))                       
     1   + dcu(2,j,1)*conjg(dcu(2,j,1)) + dcu(3,j,1)*conjg(dcu(3,j,1)))
  100 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      sum1 = sum1 + wp
      wf = real(nx*ny)*sum1/real(ffe(1,1))
      return
      end
c-----------------------------------------------------------------------
      subroutine ADDVRFIELD2(a,b,c,ndim,nxe,nye)
c this subroutine calculates a = b + c for real vector fields
      implicit none
      integer ndim, nxe, nye
      real a, b, c
      dimension a(ndim,nxe,nye), b(ndim,nxe,nye), c(ndim,nxe,nye)
c local data
      integer i, j, k
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do 30 k = 1, nye
      do 20 j = 1, nxe
      do 10 i = 1, ndim
      a(i,j,k) = b(i,j,k) + c(i,j,k)
   10 continue
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
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
      subroutine WFFT2RMX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,   
     1nxyhd)
c wrapper function for real to complex fft, with packed data
c parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
c perform y fft
         call FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c perform x fft
         call FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT2RM3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,   
     1nxyhd)
c wrapper function for 3 2d real to complex ffts
c parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2RM3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
c perform y fft
         call FFT2RM3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RM3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c perform x fft
         call FFT2RM3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT2RMN(f,ss,isign,mixup,sct,indx,indy,nxhd,nyd,ndim, 
     1nxhyd,nxyhd)
c wrapper function for multiple 2d real to complex ffts
      implicit none
      complex f, ss, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, ndim, nxhyd, nxyhd
      dimension f(ndim,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      dimension ss(ndim*nxhd,nyd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2RMNX(f,ss,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,  
     1ndim,nxhyd,nxyhd)
c perform y fft
         call FFT2RMNY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,ndim
     1,nxhyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RMNY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,ndim
     1,nxhyd,nxyhd)
c perform x fft
         call FFT2RMNX(f,ss,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,  
     1ndim,nxhyd,nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd, 
     1nxhyd,nxyhd)
c this subroutine performs the x part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic, with OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform in x is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, a forward fourier transform in x is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2
c nyd = second dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c imag(f(1,1)) = real part of mode nx/2,0 and
c imag(f(1,ny/2+1) ) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, nrxb
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 70
c inverse fourier transform
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,ani,t1,t2,t3)
      do 60 i = nyi, nyt
c bit-reverse array elements in x
      do 10 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(j1,i)
         f(j1,i) = f(j,i)
         f(j,i) = t1
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
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
   20 continue
   30 continue
   40 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 0.5/(real(nx)*real(ny))
      do 50 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      t2 = conjg(f(nxh2-j,i))
      t1 = f(j,i) + t2
      t2 = (f(j,i) - t2)*t3
      f(j,i) = ani*(t1 + t2)
      f(nxh2-j,i) = ani*conjg(t1 - t2)
   50 continue
      ani = 2.0*ani
      f(nxhh+1,i) = ani*conjg(f(nxhh+1,i))
      f(1,i) = ani*cmplx(real(f(1,i)) + aimag(f(1,i)),                  
     1                   real(f(1,i)) - aimag(f(1,i)))
   60 continue
!$OMP END PARALLEL DO
      return
c forward fourier transform
   70 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,t1,t2,t3)
      do 130 i = nyi, nyt
c scramble coefficients
      kmr = nxy/nx
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      t2 = conjg(f(nxh2-j,i))
      t1 = f(j,i) + t2
      t2 = (f(j,i) - t2)*t3
      f(j,i) = t1 + t2
      f(nxh2-j,i) = conjg(t1 - t2)
   80 continue
      f(nxhh+1,i) = 2.0*conjg(f(nxhh+1,i))
      f(1,i) = cmplx(real(f(1,i)) + aimag(f(1,i)),                      
     1               real(f(1,i)) - aimag(f(1,i)))
c bit-reverse array elements in x
      do 90 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(j1,i)
         f(j1,i) = f(j,i)
         f(j,i) = t1
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
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
  100 continue
  110 continue
  120 continue
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd, 
     1nxhyd,nxyhd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic, with OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform in y is performed
c f(n,m) = sum(f(j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform in y is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = first dimension of f >= nx/2
c nyd = second dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c imag(f(1,1)) = real part of mode nx/2,0 and
c imag(f(1,ny/2+1) ) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, nryb
      complex t1, t2
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 70
c inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,t1,t2)
      do 50 i = nxi, nxt
c bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(i,k1)
         f(i,k1) = f(i,k)
         f(i,k) = t1
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
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 60 k = 2, nyh
         t1 = f(1,ny2-k)
         f(1,ny2-k) = 0.5*cmplx(aimag(f(1,k) + t1),real(f(1,k) - t1))
         f(1,k) = 0.5*cmplx(real(f(1,k) + t1),aimag(f(1,k) - t1))
   60    continue
      endif
      return
c forward fourier transform
   70 nryb = nxhy/ny
      nry = nxy/ny
c scramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 80 k = 2, nyh
         t1 = cmplx(aimag(f(1,ny2-k)),real(f(1,ny2-k)))
         f(1,ny2-k) = conjg(f(1,k) - t1)
         f(1,k) = f(1,k) + t1
   80    continue
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,t1,t2)
      do 130 i = nxi, nxt
c bit-reverse array elements in y
      do 90 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(i,k1)
         f(i,k1) = f(i,k)
         f(i,k) = t1
      endif
   90 continue
c then transform in y
      do 120 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
  100 continue
  110 continue
  120 continue
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RM3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd, 
     1nxhyd,nxyhd)
c this subroutine performs the x part of 3 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c y, using complex arithmetic, with OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, three inverse fourier transforms are performed
c f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, three forward fourier transforms are performed
c f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = second dimension of f >= nx/2
c nyd = third dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(1:3,j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1:3,1,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c imag(f(1:3,1,1)) = real part of mode nx/2,0 and
c imag(f(1:3,1,ny/2+1) ) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      integer nrxb
      real at1, at2, ani
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 100
c inverse fourier transform
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,at2,ani,t1,t2,t3
!$OMP& ,t4)
      do 90 i = nyi, nyt
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
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
   30 continue
   40 continue
   50 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 0.5/(real(nx)*real(ny))
      do 70 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 60 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,i))
      t1 = f(jj,j,i) + t2
      t2 = (f(jj,j,i) - t2)*t3
      f(jj,j,i) = ani*(t1 + t2)
      f(jj,nxh2-j,i) = ani*conjg(t1 - t2)
   60 continue
   70 continue
      ani = 2.0*ani
      do 80 jj = 1, 3
      f(jj,nxhh+1,i) = ani*conjg(f(jj,nxhh+1,i))
      f(jj,1,i) = ani*cmplx(real(f(jj,1,i)) + aimag(f(jj,1,i)),         
     1                      real(f(jj,1,i)) - aimag(f(jj,1,i)))
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
c forward fourier transform
  100 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,at2,t1,t2,t3,t4)
      do 190 i = nyi, nyt
c scramble coefficients
      kmr = nxy/nx
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,i))
      t1 = f(jj,j,i) + t2
      t2 = (f(jj,j,i) - t2)*t3
      f(jj,j,i) = t1 + t2
      f(jj,nxh2-j,i) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 jj = 1, 3
      f(jj,nxhh+1,i) = 2.0*conjg(f(jj,nxhh+1,i))
      f(jj,1,i) = cmplx(real(f(jj,1,i)) + aimag(f(jj,1,i)),             
     1                  real(f(jj,1,i)) - aimag(f(jj,1,i)))
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
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
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
      subroutine FFT2RM3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd, 
     1nxhyd,nxyhd)
c this subroutine performs the y part of 3 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c x, using complex arithmetic, with OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, three inverse fourier transforms are performed
c f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, three forward fourier transforms are performed
c f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = second dimension of f >= nx/2
c nyd = third dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(1:3,j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1:3,1,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c imag(f(1:3,1,1)) = real part of mode nx/2,0 and
c imag(f(1:3,1,ny/2+1) ) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      integer nryb
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 80
c inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,t1,t2,t3,t4)
      do 50 i = nxi, nxt
c bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(1,i,k1)
         t2 = f(2,i,k1)
         t3 = f(3,i,k1)
         f(1,i,k1) = f(1,i,k)
         f(2,i,k1) = f(2,i,k)
         f(3,i,k1) = f(3,i,k)
         f(1,i,k) = t1
         f(2,i,k) = t2
         f(3,i,k) = t3
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
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 70 k = 2, nyh
         do 60 jj = 1, 3
         t1 = f(jj,1,ny2-k)
         f(jj,1,ny2-k) = 0.5*cmplx(aimag(f(jj,1,k) + t1),               
     1                             real(f(jj,1,k) - t1))
         f(jj,1,k) = 0.5*cmplx(real(f(jj,1,k) + t1),                    
     1                         aimag(f(jj,1,k) - t1))
   60    continue
   70    continue
      endif
      return
c forward fourier transform
   80 nryb = nxhy/ny
      nry = nxy/ny
c scramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 100 k = 2, nyh
         do 90 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
         f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
         f(jj,1,k) = f(jj,1,k) + t1
   90    continue
  100    continue
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,t1,t2,t3,t4)
      do 150 i = nxi, nxt
c bit-reverse array elements in y
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(1,i,k1)
         t2 = f(2,i,k1)
         t3 = f(3,i,k1)
         f(1,i,k1) = f(1,i,k)
         f(2,i,k1) = f(2,i,k)
         f(3,i,k1) = f(3,i,k)
         f(1,i,k) = t1
         f(2,i,k) = t2
         f(3,i,k) = t3
      endif
  110 continue
c first transform in y
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
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
  120 continue
  130 continue
  140 continue
  150 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RMNX(f,ss,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,  
     1nyd,ndim,nxhyd,nxyhd)
c this subroutine performs the x part of N two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c y, using complex arithmetic, where N = ndim, with OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)
c for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)
c where M = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, N inverse fourier transforms are performed
c f(1:N,n,m) = (1/nx*ny)*sum(f(1:N,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, N forward fourier transforms are performed
c f(1:N,j,k) = sum(f(1:N,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c ss = scratch array
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = second dimension of f
c nyd = third dimension of f
c ndim = leading dimension of array f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(1:N,j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1:N,1,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c imag(f(1:N,1,1)) = real part of mode nx/2,0 and
c imag(f(1:N,1,ny/2+1)) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, ndim, nxhyd, nxyhd
      complex f, ss, sct
      integer mixup
      dimension f(ndim,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      dimension ss(ndim*nxhd,nyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, ns, ns2, km, kmr, i, j, k, l, k1, k2, j1, j2, jj
      integer nrxb
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 110
c inverse fourier transform
      nrxb = nxhy/nxh
      nrx = nxy/nxh
c swap complex components
      call MSWAPC2N(f,ss,isign,nxh,nyi,nyt,nxhd,nyd,ndim)
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,ani,t1,t2,t3)
      do 100 i = nyi, nyt
c bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 10 jj = 1, ndim
         t1 = f(jj,j1,i)
         f(jj,j1,i) = f(jj,j,i)
         f(jj,j,i) = t1
   10    continue
      endif
   20 continue
c first transform in x
      do 60 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 jj = 1, ndim
      t2 = t1*f(jj,j2,i)
      f(jj,j2,i) = f(jj,j1,i) - t2
      f(jj,j1,i) = f(jj,j1,i) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1.0/real(2*nx*ny)
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 jj = 1, ndim
      t2 = conjg(f(jj,nxh2-j,i))
      t1 = f(jj,j,i) + t2
      t2 = (f(jj,j,i) - t2)*t3
      f(jj,j,i) = ani*(t1 + t2)
      f(jj,nxh2-j,i) = ani*conjg(t1 - t2)
   70 continue
   80 continue
      ani = 2.0*ani
      do 90 jj = 1, ndim
      f(jj,nxhh+1,i) = ani*conjg(f(jj,nxhh+1,i))
      f(jj,1,i) = ani*cmplx(real(f(jj,1,i)) + aimag(f(jj,1,i)),
     1                      real(f(jj,1,i)) - aimag(f(jj,1,i)))
   90 continue
  100 continue
!$OMP END PARALLEL DO
      return
c forward fourier transform
  110 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,t1,t2,t3)
      do 210 i = nyi, nyt
c scramble coefficients
      kmr = nxy/nx
      do 130 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 120 jj = 1, ndim
      t2 = conjg(f(jj,nxh2-j,i))
      t1 = f(jj,j,i) + t2
      t2 = (f(jj,j,i) - t2)*t3
      f(jj,j,i) = t1 + t2
      f(jj,nxh2-j,i) = conjg(t1 - t2)
  120 continue
  130 continue
      do 140 jj = 1, ndim
      f(jj,nxhh+1,i) = 2.0*conjg(f(jj,nxhh+1,i))
      f(jj,1,i) = cmplx(real(f(jj,1,i)) + aimag(f(jj,1,i)),
     1                  real(f(jj,1,i)) - aimag(f(jj,1,i)))
  140 continue
c bit-reverse array elements in x
      do 160 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 150 jj = 1, ndim
         t1 = f(jj,j1,i)
         f(jj,j1,i) = f(jj,j,i)
         f(jj,j,i) = t1
  150    continue
      endif
  160 continue
c then transform in x
      do 200 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 190 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 180 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 170 jj = 1, ndim
      t2 = t1*f(jj,j2,i)
      f(jj,j2,i) = f(jj,j1,i) - t2
      f(jj,j1,i) = f(jj,j1,i) + t2
  170 continue
  180 continue
  190 continue
  200 continue
  210 continue
!$OMP END PARALLEL DO
c swap complex components
      call MSWAPC2N(f,ss,isign,nxh,nyi,nyt,nxhd,nyd,ndim)
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RMNY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd, 
     1ndim,nxhyd,nxyhd)
c this subroutine performs the y part of N two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c x, using complex arithmetic, where N = ndim, with OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)
c for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)
c where M = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, N inverse fourier transforms are performed
c f(1:N,n,m) = (1/nx*ny)*sum(f(1:N,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, N forward fourier transforms are performed
c f(1:N,j,k) = sum(f(1:N,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = second dimension of f
c nyd = third dimension of f
c ndim = leading dimension of array f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(1:N,j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1:N,1,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c imag(f(1:N,1,1)) = real part of mode nx/2,0 and
c imag(f(1:N,1,ny/2+1)) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, ndim, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(ndim,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, ns, ns2, km, kmr, i, j, k, l, k1, k2, j1, j2, jj
      integer nryb
      complex t1, t2
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 100
c inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,t1,t2)
      do 70 i = nxi, nxt
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 10 jj = 1, ndim
         t1 = f(jj,i,k1)
         f(jj,i,k1) = f(jj,i,k)
         f(jj,i,k) = t1
   10    continue
      endif
   20 continue
c then transform in y
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 jj = 1, ndim
      t2 = t1*f(jj,i,j2)
      f(jj,i,j2) = f(jj,i,j1) - t2
      f(jj,i,j1) = f(jj,i,j1) + t2
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      do 90 k = 2, nyh
      if (nxi.eq.1) then
         do 80 jj = 1, ndim
         t1 = f(jj,1,ny2-k)
         f(jj,1,ny2-k) = 0.5*cmplx(aimag(f(jj,1,k) + t1),
     1                             real(f(jj,1,k) - t1))
         f(jj,1,k) = 0.5*cmplx(real(f(jj,1,k) + t1),
     1                         aimag(f(jj,1,k) - t1))
   80    continue
      endif
   90 continue
      return
c forward fourier transform
  100 nryb = nxhy/ny
      nry = nxy/ny
c scramble modes kx = 0, nx/2
      do 120 k = 2, nyh
      if (nxi.eq.1) then
         do 110 jj = 1, ndim
         t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
         f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
         f(jj,1,k) = f(jj,1,k) + t1
  110    continue
      endif
  120 continue
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,t1,t2)
      do 190 i = nxi, nxt
c bit-reverse array elements in y
      do 140 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 130 jj = 1, ndim
         t1 = f(jj,i,k1)
         f(jj,i,k1) = f(jj,i,k)
         f(jj,i,k) = t1
  130    continue
      endif
  140 continue
c first transform in y
      do 180 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 150 jj = 1, ndim
      t2 = t1*f(jj,i,j2)
      f(jj,i,j2) = f(jj,i,j1) - t2
      f(jj,i,j1) = f(jj,i,j1) + t2
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine MSWAPC2N(f,s,isign,nxh,nyi,nyt,nxhd,nyd,ndim)
c this subroutine swaps components for multiple ffts
c f = input array
c s = scratch array
c isign = (-1,1) = swap (real-to-complex,complex-to-real)
c nxh = complex dimension in x direction
c nyi/nyt = initial/final y index used
c nxhd = half of the second dimension of f
c nyd = third dimension of f
c ndim = leading dimension of array f
      implicit none
      integer isign, nxh, nyi, nyt, nxhd, nyd, ndim
      real f, s
      dimension f(ndim,2*nxhd,nyd), s(2*ndim*nxhd,nyd)
c local data
      integer i, j, k, ioff
c swap complex components
c real to complex
      if (isign.lt.0) then
!$OMP PARALLEL DO PRIVATE(i,j,k,ioff)
         do 60 k = nyi, nyt
         do 20 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 10 i = 1, ndim
         s(2*i+ioff-1,k) = f(i,2*j-1,k)
         s(2*i+ioff,k) = f(i,2*j,k)
   10    continue
   20    continue
         do 50 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 30 i = 1, ndim
         f(i,2*j-1,k) = s(i+ioff,k)
   30    continue
         ioff = ioff + ndim
         do 40 i = 1, ndim
         f(i,2*j,k) = s(i+ioff,k)
   40    continue
   50    continue
   60    continue
!$OMP END PARALLEL DO
c complex to real
      else if (isign.gt.0) then
!$OMP PARALLEL DO PRIVATE(i,j,k,ioff)
         do 120 k = nyi, nyt
         do 90 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 70 i = 1, ndim
         s(i+ioff,k) = f(i,2*j-1,k)
   70    continue
         ioff = ioff + ndim
         do 80 i = 1, ndim
         s(i+ioff,k) = f(i,2*j,k)
   80    continue
   90    continue
         do 110 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 100 i = 1, ndim
         f(i,2*j-1,k) = s(2*i+ioff-1,k)
         f(i,2*j,k) = s(2*i+ioff,k)
  100    continue
  110    continue
  120    continue
!$OMP END PARALLEL DO
      endif
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
      subroutine PPCOPYOUT(part,ppart,kpic,nop,nppmx,idimp,mxy1,irc)
c for 2d code, this subroutine copies segmented particle data ppart to
c the array part with original tiled layout
c input: all except part, output: part
c part(i,j) = i-th coordinate for particle j
c ppart(i,j,k) = i-th coordinate for particle j in tile k
c kpic = number of particles per tile
c nop = number of particles
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c mxy1 = total number of tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nop, nppmx, idimp, mxy1, irc
      real part, ppart
      integer kpic
      dimension part(idimp,nop), ppart(idimp,nppmx,mxy1)
      dimension kpic(mxy1)
c local data
      integer i, j, k, npoff, npp, ne, ierr
      npoff = 0
      ierr = 0
c loop over tiles
      do 30 k = 1, mxy1
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
