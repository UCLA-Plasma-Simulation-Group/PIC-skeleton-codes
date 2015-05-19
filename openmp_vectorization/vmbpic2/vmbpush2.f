c Fortran Library for Skeleton 2-1/2D Electromagnetic OpenMP/Vector
c PIC Code
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
      subroutine PPMOVIN2LT(part,ppart,kpic,nppmx,idimp,nop,mx,my,mx1,  
     1mxy1,irc)
c this subroutine sorts particles by x,y grid in tiles of mx, my
c and copies to segmented array ppart
c linear interpolation
c input: all except ppart, kpic, output: ppart, kpic
c part/ppart = input/output particle arrays
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
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
      dimension part(idimp,nop), ppart(nppmx,idimp,mxy1)
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
         ppart(ip,i,m) = part(i,j)
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
      subroutine PPMOVIN2LTP(part,ppart,kpic,kp,nppmx,idimp,nop,mx,my,  
     1mx1,mxy1,irc)
c this subroutine sorts particles by x,y grid in tiles of mx, my
c and copies to segmented array ppart
c designed for NUMA architectures, where memory is associated with the
c processor which first writes a memory location.
c linear interpolation
c input: all except ppart, kpic, output: ppart, kpic
c part/ppart = input/output particle arrays
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c kpic = output number of particles per tile
c kp = original location of reordered particle
c nppmx = rmaximum number of particles in tile
c idimp = size of phase space = 4
c nop = number of particles
c mx/my = number of grids in sorting cell in x and y
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kpic, kp, nppmx, idimp, nop, mx, my, mx1, mxy1, irc
      real part, ppart
      dimension part(idimp,nop), ppart(nppmx,idimp,mxy1)
      dimension kpic(mxy1), kp(nppmx,mxy1)
c local data
      integer i, j, k, n, m, ip, npp, ierr
      ierr = 0
c clear counter array
      do 10 k = 1, mxy1
      kpic(k) = 0
   10 continue
c find addresses of particles at each tile to reorder particles
      do 20 j = 1, nop
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = m/my
      m = n + mx1*m
      ip = kpic(m) + 1
      if (ip.le.nppmx) then
          kp(ip,m) = j
      else
         ierr = max(ierr,ip-nppmx)
      endif
      kpic(m) = ip
   20 continue
c check for overflow
      if (ierr.gt.0) then
         irc = ierr
         return
      endif
c copy reordered particles
!$OMP PARALLEL DO PRIVATE(i,j,k,m,npp)
      do 50 k = 1, mxy1
      npp = kpic(k)
      do 40 j = 1, npp
      m = kp(j,k)
      do 30 i = 1, idimp
      ppart(j,i,k) = part(i,m)
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPCHECK2LT(ppart,kpic,idimp,nppmx,nx,ny,mx,my,mx1,my1, 
     1irc)
c this subroutine performs a sanity check to make sure particles sorted
c by x,y grid in tiles of mx, my, are all within bounds.
c tiles are assumed to be arranged in 2D linear memory, and transposed
c input: all except irc
c output: irc
c ppart(n,1,k) = position x of particle n in tile k
c ppart(n,2,k) = position y of particle n in tile k
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
      dimension ppart(nppmx,idimp,mx1*my1)
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
      dx = ppart(j,1,k)
      dy = ppart(j,2,k)
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
      subroutine GBPPUSH23LT(ppart,fxy,bxy,kpic,qbm,dt,dtc,ek,idimp,    
     1nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = velocity vx of particle n in tile m
c ppart(n,4,m) = velocity vy of particle n in tile m
c ppart(n,5,m) = velocity vz of particle n in tile m
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
      dimension ppart(nppmx,idimp,mxy1)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm, lxv
      real qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz
      real sfxy, sbxy
      dimension sfxy(4,MXV*MYV), sbxy(4,MXV*MYV)
c     dimension sfxy(4,(mx+1)*(my+1)), sbxy(4,(mx+1)*(my+1))
      double precision sum1, sum2
      lxv = mx + 1
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
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,x,y,vx,vy,vz,dxp,dyp,amx,amy,dx
!$OMP& ,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,  
!$OMP& rot3,rot4,rot5,rot6,rot7,rot8,rot9,sum1,sfxy,sbxy)
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
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noff+nxv*(j+moff-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noff+nxv*(j+moff-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noff+nxv*(j+moff-1))
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noff+nxv*(j+moff-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noff+nxv*(j+moff-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noff+nxv*(j+moff-1))
   30 continue
   40 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 50 j = 1, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn)
      dy = amx*sfxy(2,nn)
      dz = amx*sfxy(3,nn)
      dx = amy*(dxp*sfxy(1,nn+1) + dx)
      dy = amy*(dxp*sfxy(2,nn+1) + dy)
      dz = amy*(dxp*sfxy(3,nn+1) + dz)
      acx = amx*sfxy(1,nn+lxv)
      acy = amx*sfxy(2,nn+lxv)
      acz = amx*sfxy(3,nn+lxv)
      dx = dx + dyp*(dxp*sfxy(1,nn+1+lxv) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1+lxv) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1+lxv) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn)
      oy = amx*sbxy(2,nn)
      oz = amx*sbxy(3,nn)
      ox = amy*(dxp*sbxy(1,nn+1) + ox)
      oy = amy*(dxp*sbxy(2,nn+1) + oy)
      oz = amy*(dxp*sbxy(3,nn+1) + oz)
      acx = amx*sbxy(1,nn+lxv)
      acy = amx*sbxy(2,nn+lxv)
      acz = amx*sbxy(3,nn+lxv)
      ox = ox + dyp*(dxp*sbxy(1,nn+1+lxv) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1+lxv) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1+lxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(j,3,k) + dx
      acy = ppart(j,4,k) + dy
      acz = ppart(j,5,k) + dz
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c new position
      dx = x + vx*dtc
      dy = y + vy*dtc
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
      endif
c set new position
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
c set new velocity
      ppart(j,3,k) = vx
      ppart(j,4,k) = vy
      ppart(j,5,k) = vz
   50 continue
      sum2 = sum2 + sum1
   60 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine GBPPUSHF23LT(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,dtc,ek
     1,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = velocity vx of particle n in tile m
c ppart(n,4,m) = velocity vy of particle n in tile m
c ppart(n,5,m) = velocity vz of particle n in tile m
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
      dimension ppart(nppmx,idimp,mxy1)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
      dimension kpic(mxy1), ncl(8,mxy1)
      dimension ihole(2,ntmax+1,mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, ih, nh, nn, mm, lxv
      real qtmh, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
      real x, y, vx, vy, vz
      real sfxy, sbxy
      dimension sfxy(4,MXV*MYV), sbxy(4,MXV*MYV)
c     dimension sfxy(4,(mx+1)*(my+1)), sbxy(4,(mx+1)*(my+1))
      double precision sum1, sum2
      lxv = mx + 1
      qtmh = 0.5*qbm*dt
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,ih,nh,x,y,vx,vy,vz,dxp,dyp,amx,
!$OMP& amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1, 
!$OMP& rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,    
!$OMP& edgery,sum1,sfxy,sbxy)
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
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noff+nxv*(j+moff-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noff+nxv*(j+moff-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noff+nxv*(j+moff-1))
   10 continue
   20 continue
      do 40 j = 1, mm+1
      do 30 i = 1, nn+1
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noff+nxv*(j+moff-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noff+nxv*(j+moff-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noff+nxv*(j+moff-1))
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
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn)
      dy = amx*sfxy(2,nn)
      dz = amx*sfxy(3,nn)
      dx = amy*(dxp*sfxy(1,nn+1) + dx)
      dy = amy*(dxp*sfxy(2,nn+1) + dy)
      dz = amy*(dxp*sfxy(3,nn+1) + dz)
      acx = amx*sfxy(1,nn+lxv)
      acy = amx*sfxy(2,nn+lxv)
      acz = amx*sfxy(3,nn+lxv)
      dx = dx + dyp*(dxp*sfxy(1,nn+1+lxv) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1+lxv) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1+lxv) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn)
      oy = amx*sbxy(2,nn)
      oz = amx*sbxy(3,nn)
      ox = amy*(dxp*sbxy(1,nn+1) + ox)
      oy = amy*(dxp*sbxy(2,nn+1) + oy)
      oz = amy*(dxp*sbxy(3,nn+1) + oz)
      acx = amx*sbxy(1,nn+lxv)
      acy = amx*sbxy(2,nn+lxv)
      acz = amx*sbxy(3,nn+lxv)
      ox = ox + dyp*(dxp*sbxy(1,nn+1+lxv) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1+lxv) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1+lxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(j,3,k) + dx
      acy = ppart(j,4,k) + dy
      acz = ppart(j,5,k) + dz
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c new position
      dx = x + vx*dtc
      dy = y + vy*dtc
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
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
c set new velocity
      ppart(j,3,k) = vx
      ppart(j,4,k) = vy
      ppart(j,5,k) = vz
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
      subroutine GRBPPUSH23LT(ppart,fxy,bxy,kpic,qbm,dt,dtc,ci,ek,idimp,
     1nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = momentum vx of particle n in tile m
c ppart(n,4,m) = momentum vy of particle n in tile m
c ppart(n,5,m) = momentum vz of particle n in tile m
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      real qbm, dt, dtc, ci, ek
      real ppart, fxy, bxy
      integer kpic
      dimension ppart(nppmx,idimp,mxy1)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm, lxv
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz
      real sfxy, sbxy
      dimension sfxy(4,MXV*MYV), sbxy(4,MXV*MYV)
c     dimension sfxy(4,(mx+1)*(my+1)), sbxy(4,(mx+1)*(my+1))
      double precision sum1, sum2
      lxv = mx + 1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
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
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,x,y,vx,vy,vz,dxp,dyp,amx,amy,dx
!$OMP& ,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,  
!$OMP& rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,dtg,sum1,sfxy,   
!$OMP& sbxy)
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
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noff+nxv*(j+moff-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noff+nxv*(j+moff-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noff+nxv*(j+moff-1))
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noff+nxv*(j+moff-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noff+nxv*(j+moff-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noff+nxv*(j+moff-1))
   30 continue
   40 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 50 j = 1, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn)
      dy = amx*sfxy(2,nn)
      dz = amx*sfxy(3,nn)
      dx = amy*(dxp*sfxy(1,nn+1) + dx)
      dy = amy*(dxp*sfxy(2,nn+1) + dy)
      dz = amy*(dxp*sfxy(3,nn+1) + dz)
      acx = amx*sfxy(1,nn+lxv)
      acy = amx*sfxy(2,nn+lxv)
      acz = amx*sfxy(3,nn+lxv)
      dx = dx + dyp*(dxp*sfxy(1,nn+1+lxv) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1+lxv) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1+lxv) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn)
      oy = amx*sbxy(2,nn)
      oz = amx*sbxy(3,nn)
      ox = amy*(dxp*sbxy(1,nn+1) + ox)
      oy = amy*(dxp*sbxy(2,nn+1) + oy)
      oz = amy*(dxp*sbxy(3,nn+1) + oz)
      acx = amx*sbxy(1,nn+lxv)
      acy = amx*sbxy(2,nn+lxv)
      acz = amx*sbxy(3,nn+lxv)
      ox = ox + dyp*(dxp*sbxy(1,nn+1+lxv) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1+lxv) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1+lxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(j,3,k) + dx
      acy = ppart(j,4,k) + dy
      acz = ppart(j,5,k) + dz
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
c new velocity
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = x + vx*dtg
      dy = y + vy*dtg
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
      endif
c set new position
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
c set new momentum
      ppart(j,3,k) = vx
      ppart(j,4,k) = vy
      ppart(j,5,k) = vz
   50 continue
      sum2 = sum2 + sum1
   60 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine GRBPPUSHF23LT(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,dtc, 
     1ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c with periodic boundary conditions.
c Using the Boris Mover.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, irc, ek
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = momentum vx of particle n in tile m
c ppart(n,4,m) = momentum vy of particle n in tile m
c ppart(n,5,m) = momentum vz of particle n in tile m
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      real qbm, dt, dtc, ci, ek
      real ppart, fxy, bxy
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mxy1)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
      dimension kpic(mxy1), ncl(8,mxy1)
      dimension ihole(2,ntmax+1,mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, ih, nh, nn, mm, lxv
      real qtmh, ci2, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, p2, gami, qtmg, dtg, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
      real x, y, vx, vy, vz
      real sfxy, sbxy
      dimension sfxy(4,MXV*MYV), sbxy(4,MXV*MYV)
c     dimension sfxy(4,(mx+1)*(my+1)), sbxy(4,(mx+1)*(my+1))
      double precision sum1, sum2
      lxv = mx + 1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,ih,nh,x,y,vx,vy,vz,dxp,dyp,amx,
!$OMP& amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1, 
!$OMP& rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,    
!$OMP& edgery,p2,gami,qtmg,dtg,sum1,sfxy,sbxy)
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
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noff+nxv*(j+moff-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noff+nxv*(j+moff-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noff+nxv*(j+moff-1))
   10 continue
   20 continue
      do 40 j = 1, mm+1
      do 30 i = 1, nn+1
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noff+nxv*(j+moff-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noff+nxv*(j+moff-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noff+nxv*(j+moff-1))
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
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn)
      dy = amx*sfxy(2,nn)
      dz = amx*sfxy(3,nn)
      dx = amy*(dxp*sfxy(1,nn+1) + dx)
      dy = amy*(dxp*sfxy(2,nn+1) + dy)
      dz = amy*(dxp*sfxy(3,nn+1) + dz)
      acx = amx*sfxy(1,nn+lxv)
      acy = amx*sfxy(2,nn+lxv)
      acz = amx*sfxy(3,nn+lxv)
      dx = dx + dyp*(dxp*sfxy(1,nn+1+lxv) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1+lxv) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1+lxv) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn)
      oy = amx*sbxy(2,nn)
      oz = amx*sbxy(3,nn)
      ox = amy*(dxp*sbxy(1,nn+1) + ox)
      oy = amy*(dxp*sbxy(2,nn+1) + oy)
      oz = amy*(dxp*sbxy(3,nn+1) + oz)
      acx = amx*sbxy(1,nn+lxv)
      acy = amx*sbxy(2,nn+lxv)
      acz = amx*sbxy(3,nn+lxv)
      ox = ox + dyp*(dxp*sbxy(1,nn+1+lxv) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1+lxv) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1+lxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(j,3,k) + dx
      acy = ppart(j,4,k) + dy
      acz = ppart(j,5,k) + dz
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = x + vx*dtg
      dy = y + vy*dtg
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
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
c set new momentum
      ppart(j,3,k) = vx
      ppart(j,4,k) = vy
      ppart(j,5,k) = vz
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
      subroutine VGBPPUSH23LT(ppart,fxy,bxy,kpic,qbm,dt,dtc,ek,idimp,   
     1nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c vectorizable/OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = velocity vx of particle n in tile m
c ppart(n,4,m) = velocity vy of particle n in tile m
c ppart(n,5,m) = velocity vz of particle n in tile m
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
      dimension ppart(nppmx,idimp,mxy1)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noff, moff, npp, ipp, joff, nps
      integer i, j, k, m, nn, mm, lxv
      real qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz
      real sfxy, sbxy
      dimension sfxy(4,MXV*MYV), sbxy(4,MXV*MYV)
c     dimension sfxy(4,(mx+1)*(my+1)), sbxy(4,(mx+1)*(my+1))
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,2)
      double precision sum1, sum2
      lxv = mx + 1
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
!$OMP& PRIVATE(i,j,k,m,noff,moff,npp,ipp,joff,nps,nn,mm,x,y,vx,vy,vz,
!$OMP& dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,
!$OMP& anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sum1,sfxy,
!$OMP& sbxy,n,s1,s2,t)
!$OMP& REDUCTION(+:sum2)
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
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noff+nxv*(j+moff-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noff+nxv*(j+moff-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noff+nxv*(j+moff-1))
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noff+nxv*(j+moff-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noff+nxv*(j+moff-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noff+nxv*(j+moff-1))
   30 continue
   40 continue
      sum1 = 0.0d0
      ipp = npp/npblk
c outer loop over number of full blocks
      do 100 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 50 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noff + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
   50 continue
c find acceleration
      do 70 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
      do 60 i = 1, lvect
      if (i.gt.2) nn = mm
      dx = dx + sfxy(1,i+nn)*s1(j,i)
      dy = dy + sfxy(2,i+nn)*s1(j,i)
      dz = dz + sfxy(3,i+nn)*s1(j,i)
      ox = ox + sbxy(1,i+nn)*s1(j,i)
      oy = oy + sbxy(2,i+nn)*s1(j,i)
      oz = oz + sbxy(3,i+nn)*s1(j,i)
   60 continue
      s1(j,1) = dx
      s1(j,2) = dy
      s1(j,3) = dz
      s2(j,1) = ox
      s2(j,2) = oy
      s2(j,3) = oz
   70 continue
c new velocity
      do 80 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
c calculate half impulse
      dx = qtmh*s1(j,1)
      dy = qtmh*s1(j,2)
      dz = qtmh*s1(j,3)
c half acceleration
      acx = ppart(j+joff,3,k) + dx
      acy = ppart(j+joff,4,k) + dy
      acz = ppart(j+joff,5,k) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*s2(j,1)
      omyt = qtmh*s2(j,2)
      omzt = qtmh*s2(j,3)
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c new position
      s1(j,1) = x + vx*dtc
      s1(j,2) = y + vy*dtc
      s2(j,1) = vx
      s2(j,2) = vy
      s2(j,3) = vz
   80 continue
! check boundary conditions
!dir$ novector
      do 90 j = 1, npblk
      dx = s1(j,1)
      dy = s1(j,2)
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
      endif
c set new position
      ppart(j+joff,1,k) = dx
      ppart(j+joff,2,k) = dy
c set new velocity
      ppart(j+joff,3,k) = vx
      ppart(j+joff,4,k) = vy
      ppart(j+joff,5,k) = vz
   90 continue
  100 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 110 j = nps, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn)
      dy = amx*sfxy(2,nn)
      dz = amx*sfxy(3,nn)
      dx = amy*(dxp*sfxy(1,nn+1) + dx)
      dy = amy*(dxp*sfxy(2,nn+1) + dy)
      dz = amy*(dxp*sfxy(3,nn+1) + dz)
      acx = amx*sfxy(1,nn+lxv)
      acy = amx*sfxy(2,nn+lxv)
      acz = amx*sfxy(3,nn+lxv)
      dx = dx + dyp*(dxp*sfxy(1,nn+1+lxv) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1+lxv) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1+lxv) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn)
      oy = amx*sbxy(2,nn)
      oz = amx*sbxy(3,nn)
      ox = amy*(dxp*sbxy(1,nn+1) + ox)
      oy = amy*(dxp*sbxy(2,nn+1) + oy)
      oz = amy*(dxp*sbxy(3,nn+1) + oz)
      acx = amx*sbxy(1,nn+lxv)
      acy = amx*sbxy(2,nn+lxv)
      acz = amx*sbxy(3,nn+lxv)
      ox = ox + dyp*(dxp*sbxy(1,nn+1+lxv) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1+lxv) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1+lxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(j,3,k) + dx
      acy = ppart(j,4,k) + dy
      acz = ppart(j,5,k) + dz
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c new position
      dx = x + vx*dtc
      dy = y + vy*dtc
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
      endif
c set new position
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
c set new velocity
      ppart(j,3,k) = vx
      ppart(j,4,k) = vy
      ppart(j,5,k) = vz
  110 continue
      sum2 = sum2 + sum1
  120 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine VGBPPUSHF23LT(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,dtc, 
     1ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c vectorizable/OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = velocity vx of particle n in tile m
c ppart(n,4,m) = velocity vy of particle n in tile m
c ppart(n,5,m) = velocity vz of particle n in tile m
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
      dimension ppart(nppmx,idimp,mxy1)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
      dimension kpic(mxy1), ncl(8,mxy1)
      dimension ihole(2,ntmax+1,mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noff, moff, npp, ipp, joff, nps
      integer i, j, k, m, ih, nh, nn, mm, lxv
      real qtmh, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
      real x, y, vx, vy, vz
      real sfxy, sbxy
      dimension sfxy(4,MXV*MYV), sbxy(4,MXV*MYV)
c     dimension sfxy(4,(mx+1)*(my+1)), sbxy(4,(mx+1)*(my+1))
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,2)
      double precision sum1, sum2
      lxv = mx + 1
      qtmh = 0.5*qbm*dt
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,noff,moff,npp,ipp,joff,nps,nn,mm,ih,nh,x,y,vx,vy,
!$OMP& vz,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt, 
!$OMP& omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,   
!$OMP& edgely,edgerx,edgery,sum1,sfxy,sbxy,n,s1,s2,t)
!$OMP& REDUCTION(+:sum2)
      do 130 k = 1, mxy1
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
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noff+nxv*(j+moff-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noff+nxv*(j+moff-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noff+nxv*(j+moff-1))
   10 continue
   20 continue
      do 40 j = 1, mm+1
      do 30 i = 1, nn+1
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noff+nxv*(j+moff-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noff+nxv*(j+moff-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noff+nxv*(j+moff-1))
   30 continue
   40 continue
c clear counters
      do 50 j = 1, 8
      ncl(j,k) = 0
   50 continue
      sum1 = 0.0d0
      ipp = npp/npblk
c outer loop over number of full blocks
      do 110 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 60 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noff + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
   60 continue
c find acceleration
      do 80 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
      do 70 i = 1, lvect
      if (i.gt.2) nn = mm
      dx = dx + sfxy(1,i+nn)*s1(j,i)
      dy = dy + sfxy(2,i+nn)*s1(j,i)
      dz = dz + sfxy(3,i+nn)*s1(j,i)
      ox = ox + sbxy(1,i+nn)*s1(j,i)
      oy = oy + sbxy(2,i+nn)*s1(j,i)
      oz = oz + sbxy(3,i+nn)*s1(j,i)
   70 continue
      s1(j,1) = dx
      s1(j,2) = dy
      s1(j,3) = dz
      s2(j,1) = ox
      s2(j,2) = oy
      s2(j,3) = oz
   80 continue
c new velocity
      do 90 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
c calculate half impulse
      dx = qtmh*s1(j,1)
      dy = qtmh*s1(j,2)
      dz = qtmh*s1(j,3)
c half acceleration
      acx = ppart(j+joff,3,k) + dx
      acy = ppart(j+joff,4,k) + dy
      acz = ppart(j+joff,5,k) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*s2(j,1)
      omyt = qtmh*s2(j,2)
      omzt = qtmh*s2(j,3)
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c new position
      s1(j,1) = x + vx*dtc
      s1(j,2) = y + vy*dtc
      s2(j,1) = vx
      s2(j,2) = vy
      s2(j,3) = vz
   90 continue
! check boundary conditions
!dir$ novector
      do 100 j = 1, npblk
      dx = s1(j,1)
      dy = s1(j,2)
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
      ppart(j+joff,1,k) = dx
      ppart(j+joff,2,k) = dy
c set new velocity
      ppart(j+joff,3,k) = s2(j,1)
      ppart(j+joff,4,k) = s2(j,2)
      ppart(j+joff,5,k) = s2(j,3)
c increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j + joff
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
  100 continue
  110 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 120 j = nps, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn)
      dy = amx*sfxy(2,nn)
      dz = amx*sfxy(3,nn)
      dx = amy*(dxp*sfxy(1,nn+1) + dx)
      dy = amy*(dxp*sfxy(2,nn+1) + dy)
      dz = amy*(dxp*sfxy(3,nn+1) + dz)
      acx = amx*sfxy(1,nn+lxv)
      acy = amx*sfxy(2,nn+lxv)
      acz = amx*sfxy(3,nn+lxv)
      dx = dx + dyp*(dxp*sfxy(1,nn+1+lxv) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1+lxv) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1+lxv) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn)
      oy = amx*sbxy(2,nn)
      oz = amx*sbxy(3,nn)
      ox = amy*(dxp*sbxy(1,nn+1) + ox)
      oy = amy*(dxp*sbxy(2,nn+1) + oy)
      oz = amy*(dxp*sbxy(3,nn+1) + oz)
      acx = amx*sbxy(1,nn+lxv)
      acy = amx*sbxy(2,nn+lxv)
      acz = amx*sbxy(3,nn+lxv)
      ox = ox + dyp*(dxp*sbxy(1,nn+1+lxv) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1+lxv) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1+lxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(j,3,k) + dx
      acy = ppart(j,4,k) + dy
      acz = ppart(j,5,k) + dz
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c new position
      dx = x + vx*dtc
      dy = y + vy*dtc
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
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
c set new velocity
      ppart(j,3,k) = vx
      ppart(j,4,k) = vy
      ppart(j,5,k) = vz
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
  120 continue
      sum2 = sum2 + sum1
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
  130 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine VGRBPPUSH23LT(ppart,fxy,bxy,kpic,qbm,dt,dtc,ci,ek,idimp
     1,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c vectorizable/OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = momentum vx of particle n in tile m
c ppart(n,4,m) = momentum vy of particle n in tile m
c ppart(n,5,m) = momentum vz of particle n in tile m
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      real qbm, dt, dtc, ci, ek
      real ppart, fxy, bxy
      integer kpic
      dimension ppart(nppmx,idimp,mxy1)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noff, moff, npp, ipp, joff, nps
      integer i, j, k, m, nn, mm, lxv
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz
      real sfxy, sbxy
      dimension sfxy(4,MXV*MYV), sbxy(4,MXV*MYV)
c     dimension sfxy(4,(mx+1)*(my+1)), sbxy(4,(mx+1)*(my+1))
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,2)
      double precision sum1, sum2
      lxv = mx + 1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
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
!$OMP& PRIVATE(i,j,k,m,noff,moff,npp,nn,mm,ipp,joff,nps,x,y,vx,vy,vz,
!$OMP& dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,
!$OMP& anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,
!$OMP& dtg,sum1,sfxy,sbxy,n,s1,s2,t)
!$OMP& REDUCTION(+:sum2)
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
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noff+nxv*(j+moff-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noff+nxv*(j+moff-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noff+nxv*(j+moff-1))
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noff+nxv*(j+moff-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noff+nxv*(j+moff-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noff+nxv*(j+moff-1))
   30 continue
   40 continue
      sum1 = 0.0d0
      ipp = npp/npblk
c outer loop over number of full blocks
      do 100 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 50 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noff + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
   50 continue
c find acceleration
      do 70 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
      do 60 i = 1, lvect
      if (i.gt.2) nn = mm
      dx = dx + sfxy(1,i+nn)*s1(j,i)
      dy = dy + sfxy(2,i+nn)*s1(j,i)
      dz = dz + sfxy(3,i+nn)*s1(j,i)
      ox = ox + sbxy(1,i+nn)*s1(j,i)
      oy = oy + sbxy(2,i+nn)*s1(j,i)
      oz = oz + sbxy(3,i+nn)*s1(j,i)
   60 continue
      s1(j,1) = dx
      s1(j,2) = dy
      s1(j,3) = dz
      s2(j,1) = ox
      s2(j,2) = oy
      s2(j,3) = oz
   70 continue
c new momentum
      do 80 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
c calculate half impulse
      dx = qtmh*s1(j,1)
      dy = qtmh*s1(j,2)
      dz = qtmh*s1(j,3)
c half acceleration
      acx = ppart(j+joff,3,k) + dx
      acy = ppart(j+joff,4,k) + dy
      acz = ppart(j+joff,5,k) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*s2(j,1)
      omyt = qtmg*s2(j,2)
      omzt = qtmg*s2(j,3)
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      s1(j,1) = x + vx*dtg
      s1(j,2) = y + vy*dtg
      s2(j,1) = vx
      s2(j,2) = vy
      s2(j,3) = vz
   80 continue
! check boundary conditions
!dir$ novector
      do 90 j = 1, npblk
      dx = s1(j,1)
      dy = s1(j,2)
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
      endif
c set new position
      ppart(j+joff,1,k) = dx
      ppart(j+joff,2,k) = dy
c set new momentum
      ppart(j+joff,3,k) = vx
      ppart(j+joff,4,k) = vy
      ppart(j+joff,5,k) = vz
   90 continue
  100 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 110 j = nps, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn)
      dy = amx*sfxy(2,nn)
      dz = amx*sfxy(3,nn)
      dx = amy*(dxp*sfxy(1,nn+1) + dx)
      dy = amy*(dxp*sfxy(2,nn+1) + dy)
      dz = amy*(dxp*sfxy(3,nn+1) + dz)
      acx = amx*sfxy(1,nn+lxv)
      acy = amx*sfxy(2,nn+lxv)
      acz = amx*sfxy(3,nn+lxv)
      dx = dx + dyp*(dxp*sfxy(1,nn+1+lxv) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1+lxv) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1+lxv) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn)
      oy = amx*sbxy(2,nn)
      oz = amx*sbxy(3,nn)
      ox = amy*(dxp*sbxy(1,nn+1) + ox)
      oy = amy*(dxp*sbxy(2,nn+1) + oy)
      oz = amy*(dxp*sbxy(3,nn+1) + oz)
      acx = amx*sbxy(1,nn+lxv)
      acy = amx*sbxy(2,nn+lxv)
      acz = amx*sbxy(3,nn+lxv)
      ox = ox + dyp*(dxp*sbxy(1,nn+1+lxv) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1+lxv) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1+lxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(j,3,k) + dx
      acy = ppart(j,4,k) + dy
      acz = ppart(j,5,k) + dz
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
c new velocity
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = x + vx*dtg
      dy = y + vy*dtg
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
      endif
c set new position
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
c set new momentum
      ppart(j,3,k) = vx
      ppart(j,4,k) = vy
      ppart(j,5,k) = vz
  110 continue

      sum2 = sum2 + sum1
  120 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine VGRBPPUSHF23LT(ppart,fxy,bxy,kpic,ncl,ihole,qbm,dt,dtc,
     1ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c with periodic boundary conditions.
c Using the Boris Mover.
c also determines list of particles which are leaving this tile
c vectorizable/OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, irc, ek
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = momentum vx of particle n in tile m
c ppart(n,4,m) = momentum vy of particle n in tile m
c ppart(n,5,m) = momentum vz of particle n in tile m
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      real qbm, dt, dtc, ci, ek
      real ppart, fxy, bxy
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mxy1)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
      dimension kpic(mxy1), ncl(8,mxy1)
      dimension ihole(2,ntmax+1,mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noff, moff, npp, ipp, joff, nps
      integer i, j, k, m, ih, nh, nn, mm, lxv
      real qtmh, ci2, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, p2, gami, qtmg, dtg, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
      real x, y, vx, vy, vz
      real sfxy, sbxy
      dimension sfxy(4,MXV*MYV), sbxy(4,MXV*MYV)
c     dimension sfxy(4,(mx+1)*(my+1)), sbxy(4,(mx+1)*(my+1))
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,2)
      double precision sum1, sum2
      lxv = mx + 1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,noff,moff,npp,nn,mm,ipp,joff,nps,ih,nh,x,y,vx,vy,
!$OMP& vz,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt, 
!$OMP& omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,   
!$OMP& edgely,edgerx,edgery,p2,gami,qtmg,dtg,sum1,sfxy,sbxy,n,s1,s2,t)
!$OMP& REDUCTION(+:sum2)
      do 130 k = 1, mxy1
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
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noff+nxv*(j+moff-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noff+nxv*(j+moff-1))
      sfxy(3,i+lxv*(j-1)) = fxy(3,i+noff+nxv*(j+moff-1))
   10 continue
   20 continue
      do 40 j = 1, mm+1
      do 30 i = 1, nn+1
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noff+nxv*(j+moff-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noff+nxv*(j+moff-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noff+nxv*(j+moff-1))
   30 continue
   40 continue
c clear counters
      do 50 j = 1, 8
      ncl(j,k) = 0
   50 continue
      sum1 = 0.0d0
      ipp = npp/npblk
c outer loop over number of full blocks
      do 110 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 60 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noff + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
   60 continue
c find acceleration
      do 80 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
      do 70 i = 1, lvect
      if (i.gt.2) nn = mm
      dx = dx + sfxy(1,i+nn)*s1(j,i)
      dy = dy + sfxy(2,i+nn)*s1(j,i)
      dz = dz + sfxy(3,i+nn)*s1(j,i)
      ox = ox + sbxy(1,i+nn)*s1(j,i)
      oy = oy + sbxy(2,i+nn)*s1(j,i)
      oz = oz + sbxy(3,i+nn)*s1(j,i)
   70 continue
      s1(j,1) = dx
      s1(j,2) = dy
      s1(j,3) = dz
      s2(j,1) = ox
      s2(j,2) = oy
      s2(j,3) = oz
   80 continue
c new momentum
      do 90 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
c calculate half impulse
      dx = qtmh*s1(j,1)
      dy = qtmh*s1(j,2)
      dz = qtmh*s1(j,3)
c half acceleration
      acx = ppart(j+joff,3,k) + dx
      acy = ppart(j+joff,4,k) + dy
      acz = ppart(j+joff,5,k) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*s2(j,1)
      omyt = qtmg*s2(j,2)
      omzt = qtmg*s2(j,3)
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      s1(j,1) = x + vx*dtg
      s1(j,2) = y + vy*dtg
      s2(j,1) = vx
      s2(j,2) = vy
      s2(j,3) = vz
   90 continue
! check boundary conditions
!dir$ novector
      do 100 j = 1, npblk
      dx = s1(j,1)
      dy = s1(j,2)
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
      ppart(j+joff,1,k) = dx
      ppart(j+joff,2,k) = dy
c set new momentum
      ppart(j+joff,3,k) = s2(j,1)
      ppart(j+joff,4,k) = s2(j,2)
      ppart(j+joff,5,k) = s2(j,3)
c increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j + joff
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
  100 continue
  110 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 120 j = nps, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn)
      dy = amx*sfxy(2,nn)
      dz = amx*sfxy(3,nn)
      dx = amy*(dxp*sfxy(1,nn+1) + dx)
      dy = amy*(dxp*sfxy(2,nn+1) + dy)
      dz = amy*(dxp*sfxy(3,nn+1) + dz)
      acx = amx*sfxy(1,nn+lxv)
      acy = amx*sfxy(2,nn+lxv)
      acz = amx*sfxy(3,nn+lxv)
      dx = dx + dyp*(dxp*sfxy(1,nn+1+lxv) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1+lxv) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1+lxv) + acz)
c find magnetic field
      ox = amx*sbxy(1,nn)
      oy = amx*sbxy(2,nn)
      oz = amx*sbxy(3,nn)
      ox = amy*(dxp*sbxy(1,nn+1) + ox)
      oy = amy*(dxp*sbxy(2,nn+1) + oy)
      oz = amy*(dxp*sbxy(3,nn+1) + oz)
      acx = amx*sbxy(1,nn+lxv)
      acy = amx*sbxy(2,nn+lxv)
      acz = amx*sbxy(3,nn+lxv)
      ox = ox + dyp*(dxp*sbxy(1,nn+1+lxv) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1+lxv) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1+lxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = ppart(j,3,k) + dx
      acy = ppart(j,4,k) + dy
      acz = ppart(j,5,k) + dz
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
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c update inverse gamma
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = x + vx*dtg
      dy = y + vy*dtg
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
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
c set new momentum
      ppart(j,3,k) = vx
      ppart(j,4,k) = vy
      ppart(j,5,k) = vz
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
  120 continue
      sum2 = sum2 + sum1
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
  130 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine GPPOST2LT(ppart,q,kpic,qm,nppmx,idimp,mx,my,nxv,nyv,mx1
     1,mxy1)
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
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
      dimension ppart(nppmx,idimp,mxy1), q(nxv,nyv)
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
      x = ppart(j,1,k)
      y = ppart(j,2,k)
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
      subroutine VGPPOST2LT(ppart,q,kpic,qm,nppmx,idimp,mx,my,nxv,nyv,  
     1mx1,mxy1)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c vectorizable/OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
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
      dimension ppart(nppmx,idimp,mxy1), q(nxv*nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noff, moff, npp, ipp, joff, nps
      integer i, j, k, m, nn, mm, lxv
      real x, y, dxp, dyp, amx, amy
      real sq
c     dimension sq(MXV*MYV)
      dimension sq((mx+1)*(my+1))
c scratch arrays
      integer n
      real s
      dimension n(npblk), s(npblk,lvect)
      lxv = mx + 1
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,noff,moff,npp,ipp,joff,nps,nn,mm,x,y,dxp,dyp,amx,
!$OMP& amy,sq,n,s)
      do 110 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c zero out local accumulator
      do 10 j = 1, (mx+1)*(my+1)
      sq(j) = 0.0
   10 continue
c loop over particles in tile
      ipp = npp/npblk
c outer loop over number of full blocks
      do 50 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 20 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      n(j) = nn - noff + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
   20 continue
c deposit charge within tile to local accumulator
      do 40 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
!dir$ ivdep
      do 30 i = 1, lvect
      if (i.gt.2) nn = mm
      sq(i+nn) = sq(i+nn) + s(j,i)
   30 continue
   40 continue
   50 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 60 j = nps, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit charge within tile to local accumulator
      x = sq(nn) + amx*amy
      y = sq(nn+1) + dxp*amy
      sq(nn) = x
      sq(nn+1) = y
      x = sq(nn+lxv) + amx*dyp
      y = sq(nn+1+lxv) + dxp*dyp
      sq(nn+lxv) = x
      sq(nn+1+lxv) = y
   60 continue
c deposit charge to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 80 j = 2, mm
      do 70 i = 2, nn
      q(i+noff+nxv*(j+moff-1)) = q(i+noff+nxv*(j+moff-1)) +             
     1sq(i+lxv*(j-1))
   70 continue
   80 continue
c deposit charge to edge points in global array
      mm = min(my+1,nyv-moff)
      do 90 i = 2, nn
!$OMP ATOMIC
      q(i+noff+nxv*moff) = q(i+noff+nxv*moff) + sq(i)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noff+nxv*(mm+moff-1)) = q(i+noff+nxv*(mm+moff-1)) +        
     1sq(i+lxv*(mm-1))
      endif
   90 continue
      nn = min(mx+1,nxv-noff)
      do 100 j = 1, mm
!$OMP ATOMIC
      q(1+noff+nxv*(j+moff-1)) = q(1+noff+nxv*(j+moff-1)) +             
     1sq(1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noff+nxv*(j+moff-1)) = q(nn+noff+nxv*(j+moff-1)) +        
     1sq(nn+lxv*(j-1))
      endif
  100 continue
  110 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GJPPOST2LT(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny,mx,my,
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = x velocity of particle n in tile m
c ppart(n,4,m) = y velocity of particle n in tile m
c ppart(n,5,m) = z velocity of particle n in tile m
c cu(i,j,k) = ith component of current density at grid point j,k
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ipbc
      real qm, dt
      real ppart, cu
      integer kpic
      dimension ppart(nppmx,idimp,mxy1), cu(4,nxv*nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm, lxv
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real scu
      dimension scu(4,MXV*MYV)
c     dimension scu(4,(mx+1)*(my+1))
      lxv = mx + 1
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
      do 70 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c zero out local accumulator
      nn = lxv*(my + 1)
      do 10 i = 1, nn
      scu(1,i) = 0.0
      scu(2,i) = 0.0
      scu(3,i) = 0.0
   10 continue
c loop over particles in tile
      do 20 j = 1, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ppart(j,3,k)
      vy = ppart(j,4,k)
      vz = ppart(j,5,k)
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j,3,k) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(j,4,k) = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j,3,k) = -vx
         endif
      endif
c set new position
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
   20 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 40 j = 2, mm
      do 30 i = 2, nn
      cu(1,i+noff+nxv*(j+moff-1)) = cu(1,i+noff+nxv*(j+moff-1))         
     1 + scu(1,i+lxv*(j-1))
      cu(2,i+noff+nxv*(j+moff-1)) = cu(2,i+noff+nxv*(j+moff-1))         
     1 + scu(2,i+lxv*(j-1))
      cu(3,i+noff+nxv*(j+moff-1)) = cu(3,i+noff+nxv*(j+moff-1))         
     1 + scu(3,i+lxv*(j-1))
   30 continue
   40 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 50 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noff+nxv*moff) = cu(1,i+noff+nxv*moff) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noff+nxv*moff) = cu(2,i+noff+nxv*moff) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noff+nxv*moff) = cu(3,i+noff+nxv*moff) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noff+nxv*(mm+moff-1)) = cu(1,i+noff+nxv*(mm+moff-1))    
     1    + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noff+nxv*(mm+moff-1)) = cu(2,i+noff+nxv*(mm+moff-1))    
     1    + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noff+nxv*(mm+moff-1)) = cu(3,i+noff+nxv*(mm+moff-1))    
     1    + scu(3,i+lxv*(mm-1))
      endif
   50 continue
      nn = min(mx+1,nxv-noff)
      do 60 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noff+nxv*(j+moff-1)) = cu(1,1+noff+nxv*(j+moff-1))         
     1 + scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noff+nxv*(j+moff-1)) = cu(2,1+noff+nxv*(j+moff-1))         
     1 + scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noff+nxv*(j+moff-1)) = cu(3,1+noff+nxv*(j+moff-1))         
     1 + scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff+nxv*(j+moff-1)) = cu(1,nn+noff+nxv*(j+moff-1))    
     1    + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noff+nxv*(j+moff-1)) = cu(2,nn+noff+nxv*(j+moff-1))    
     1    + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noff+nxv*(j+moff-1)) = cu(3,nn+noff+nxv*(j+moff-1))    
     1    + scu(3,nn+lxv*(j-1))
      endif
   60 continue
   70 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GJPPOSTF2LT(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,idimp, 
     1nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = x velocity of particle n in tile m
c ppart(n,4,m) = y velocity of particle n in tile m
c ppart(n,5,m) = z velocity of particle n in tile m
c cu(i,j,k) = ith component of current density at grid point j,k
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ntmax
      integer irc
      real qm, dt
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mxy1), cu(4,nxv*nyv)
      dimension kpic(mxy1), ncl(8,mxy1)
      dimension ihole(2,ntmax+1,mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, ih, nh, nn, mm, lxv
      real dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real anx, any, edgelx, edgely, edgerx, edgery
      real scu
      dimension scu(4,MXV*MYV)
c     dimension scu(4,(mx+1)*(my+1))
      lxv = mx + 1
      anx = real(nx)
      any = real(ny)
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,ih,nh,x,y,dxp,dyp,amx,amy,dx,dy
!$OMP& ,vx,vy,vz,edgelx,edgely,edgerx,edgery,scu)
      do 80 k = 1, mxy1
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
c zero out local accumulator
      nn = lxv*(my + 1)
      do 10 i = 1, nn
      scu(1,i) = 0.0
      scu(2,i) = 0.0
      scu(3,i) = 0.0
   10 continue
c clear counters
      do 20 j = 1, 8
      ncl(j,k) = 0
   20 continue
c loop over particles in tile
      do 30 j = 1, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ppart(j,3,k)
      vy = ppart(j,4,k)
      vz = ppart(j,5,k)
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
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
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
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
   30 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 50 j = 2, mm
      do 40 i = 2, nn
      cu(1,i+noff+nxv*(j+moff-1)) = cu(1,i+noff+nxv*(j+moff-1))         
     1 + scu(1,i+lxv*(j-1))
      cu(2,i+noff+nxv*(j+moff-1)) = cu(2,i+noff+nxv*(j+moff-1))         
     1 + scu(2,i+lxv*(j-1))
      cu(3,i+noff+nxv*(j+moff-1)) = cu(3,i+noff+nxv*(j+moff-1))         
     1 + scu(3,i+lxv*(j-1))
   40 continue
   50 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 60 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noff+nxv*moff) = cu(1,i+noff+nxv*moff) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noff+nxv*moff) = cu(2,i+noff+nxv*moff) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noff+nxv*moff) = cu(3,i+noff+nxv*moff) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noff+nxv*(mm+moff-1)) = cu(1,i+noff+nxv*(mm+moff-1))    
     1    + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noff+nxv*(mm+moff-1)) = cu(2,i+noff+nxv*(mm+moff-1))    
     1    + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noff+nxv*(mm+moff-1)) = cu(3,i+noff+nxv*(mm+moff-1))    
     1    + scu(3,i+lxv*(mm-1))
      endif
   60 continue
      nn = min(mx+1,nxv-noff)
      do 70 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noff+nxv*(j+moff-1)) = cu(1,1+noff+nxv*(j+moff-1))         
     1 + scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noff+nxv*(j+moff-1)) = cu(2,1+noff+nxv*(j+moff-1))         
     1 + scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noff+nxv*(j+moff-1)) = cu(3,1+noff+nxv*(j+moff-1))         
     1 + scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff+nxv*(j+moff-1)) = cu(1,nn+noff+nxv*(j+moff-1))    
     1    + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noff+nxv*(j+moff-1)) = cu(2,nn+noff+nxv*(j+moff-1))    
     1    + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noff+nxv*(j+moff-1)) = cu(3,nn+noff+nxv*(j+moff-1))    
     1    + scu(3,nn+lxv*(j-1))
      endif
   70 continue
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   80 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GRJPPOST2LT(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx,ny,mx
     1,my,nxv,nyv,mx1,mxy1,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = x momentum of particle n in tile m
c ppart(n,4,m) = y momentum of particle n in tile m
c ppart(n,5,m) = z momentum of particle n in tile m
c cu(i,j,k) = ith component of current density at grid point j,k
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ipbc
      real qm, dt, ci
      real ppart, cu
      integer kpic
      dimension ppart(nppmx,idimp,mxy1), cu(4,nxv*nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm, lxv
      real ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami
      real scu
      dimension scu(4,MXV*MYV)
c     dimension scu(4,(mx+1)*(my+1))
      lxv = mx + 1
      ci2 = ci*ci
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
!$OMP& ,vz,ux,uy,uz,p2,gami,scu)
      do 70 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c zero out local accumulator
      nn = lxv*(my + 1)
      do 10 i = 1, nn
      scu(1,i) = 0.0
      scu(2,i) = 0.0
      scu(3,i) = 0.0
   10 continue
c loop over particles in tile
      do 20 j = 1, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
c find inverse gamma
      ux = ppart(j,3,k)
      uy = ppart(j,4,k)
      uz = ppart(j,5,k)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ux*gami
      vy = uy*gami
      vz = uz*gami
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j,3,k) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(j,4,k) = -uy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j,3,k) = -ux
         endif
      endif
c set new position
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
   20 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 40 j = 2, mm
      do 30 i = 2, nn
      cu(1,i+noff+nxv*(j+moff-1)) = cu(1,i+noff+nxv*(j+moff-1))         
     1 + scu(1,i+lxv*(j-1))
      cu(2,i+noff+nxv*(j+moff-1)) = cu(2,i+noff+nxv*(j+moff-1))         
     1 + scu(2,i+lxv*(j-1))
      cu(3,i+noff+nxv*(j+moff-1)) = cu(3,i+noff+nxv*(j+moff-1))         
     1 + scu(3,i+lxv*(j-1))
   30 continue
   40 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 50 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noff+nxv*moff) = cu(1,i+noff+nxv*moff) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noff+nxv*moff) = cu(2,i+noff+nxv*moff) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noff+nxv*moff) = cu(3,i+noff+nxv*moff) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noff+nxv*(mm+moff-1)) = cu(1,i+noff+nxv*(mm+moff-1))    
     1    + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noff+nxv*(mm+moff-1)) = cu(2,i+noff+nxv*(mm+moff-1))    
     1    + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noff+nxv*(mm+moff-1)) = cu(3,i+noff+nxv*(mm+moff-1))    
     1    + scu(3,i+lxv*(mm-1))
      endif
   50 continue
      nn = min(mx+1,nxv-noff)
      do 60 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noff+nxv*(j+moff-1)) = cu(1,1+noff+nxv*(j+moff-1))         
     1 + scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noff+nxv*(j+moff-1)) = cu(2,1+noff+nxv*(j+moff-1))         
     1 + scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noff+nxv*(j+moff-1)) = cu(3,1+noff+nxv*(j+moff-1))         
     1 + scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff+nxv*(j+moff-1)) = cu(1,nn+noff+nxv*(j+moff-1))    
     1    + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noff+nxv*(j+moff-1)) = cu(2,nn+noff+nxv*(j+moff-1))    
     1    + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noff+nxv*(j+moff-1)) = cu(3,nn+noff+nxv*(j+moff-1))    
     1    + scu(3,nn+lxv*(j-1))
      endif
   60 continue
   70 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GRJPPOSTF2LT(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx,   
     1idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = x momentum of particle n in tile m
c ppart(n,4,m) = y momentum of particle n in tile m
c ppart(n,5,m) = z momentum of particle n in tile m
c cu(i,j,k) = ith component of current density at grid point j,k
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ntmax
      integer irc
      real qm, dt, ci
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mxy1), cu(4,nxv*nyv)
      dimension kpic(mxy1), ncl(8,mxy1)
      dimension ihole(2,ntmax+1,mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, ih, nh, nn, mm, lxv
      real ci2, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami
      real anx, any, edgelx, edgely, edgerx, edgery
      real scu
      dimension scu(4,MXV*MYV)
c     dimension scu(4,(mx+1)*(my+1))
      lxv = mx + 1
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,ih,nh,x,y,dxp,dyp,amx,amy,dx,dy
!$OMP& ,vx,vy,vz,ux,uy,uz,edgelx,edgely,edgerx,edgery,p2,gami,scu)
      do 80 k = 1, mxy1
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
c zero out local accumulator
      nn = lxv*(my + 1)
      do 10 i = 1, nn
      scu(1,i) = 0.0
      scu(2,i) = 0.0
      scu(3,i) = 0.0
   10 continue
c clear counters
      do 20 j = 1, 8
      ncl(j,k) = 0
   20 continue
c loop over particles in tile
      do 30 j = 1, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
c find inverse gamma
      ux = ppart(j,3,k)
      uy = ppart(j,4,k)
      uz = ppart(j,5,k)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ux*gami
      vy = uy*gami
      vz = uz*gami
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
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
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
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
   30 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 50 j = 2, mm
      do 40 i = 2, nn
      cu(1,i+noff+nxv*(j+moff-1)) = cu(1,i+noff+nxv*(j+moff-1))         
     1 + scu(1,i+lxv*(j-1))
      cu(2,i+noff+nxv*(j+moff-1)) = cu(2,i+noff+nxv*(j+moff-1))         
     1 + scu(2,i+lxv*(j-1))
      cu(3,i+noff+nxv*(j+moff-1)) = cu(3,i+noff+nxv*(j+moff-1))         
     1 + scu(3,i+lxv*(j-1))
   40 continue
   50 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 60 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noff+nxv*moff) = cu(1,i+noff+nxv*moff) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noff+nxv*moff) = cu(2,i+noff+nxv*moff) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noff+nxv*moff) = cu(3,i+noff+nxv*moff) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noff+nxv*(mm+moff-1)) = cu(1,i+noff+nxv*(mm+moff-1))    
     1    + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noff+nxv*(mm+moff-1)) = cu(2,i+noff+nxv*(mm+moff-1))    
     1    + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noff+nxv*(mm+moff-1)) = cu(3,i+noff+nxv*(mm+moff-1))    
     1    + scu(3,i+lxv*(mm-1))
      endif
   60 continue
      nn = min(mx+1,nxv-noff)
      do 70 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noff+nxv*(j+moff-1)) = cu(1,1+noff+nxv*(j+moff-1))         
     1 + scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noff+nxv*(j+moff-1)) = cu(2,1+noff+nxv*(j+moff-1))         
     1 + scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noff+nxv*(j+moff-1)) = cu(3,1+noff+nxv*(j+moff-1))         
     1 + scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff+nxv*(j+moff-1)) = cu(1,nn+noff+nxv*(j+moff-1))    
     1    + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noff+nxv*(j+moff-1)) = cu(2,nn+noff+nxv*(j+moff-1))    
     1    + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noff+nxv*(j+moff-1)) = cu(3,nn+noff+nxv*(j+moff-1))    
     1    + scu(3,nn+lxv*(j-1))
      endif
   70 continue
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   80 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine VGJPPOST2LT(ppart,cu,kpic,qm,dt,nppmx,idimp,nx,ny,mx,my
     1,nxv,nyv,mx1,mxy1,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c vectorizable/OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = x velocity of particle n in tile m
c ppart(n,4,m) = y velocity of particle n in tile m
c ppart(n,5,m) = z velocity of particle n in tile m
c cu(i,j,k) = ith component of current density at grid point j,k
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ipbc
      real qm, dt
      real ppart, cu
      integer kpic
      dimension ppart(nppmx,idimp,mxy1), cu(4,nxv*nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noff, moff, npp, ipp, joff, nps
      integer i, j, k, m, nn, mm, lxv
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real scu
      dimension scu(4,MXV*MYV)
c     dimension scu(4,(mx+1)*(my+1))
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,2)
      lxv = mx + 1
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
!$OMP& PRIVATE(i,j,k,m,noff,moff,npp,ipp,joff,nps,nn,mm,x,y,dxp,dyp,amx,
!$OMP& amy,dx,dy,vx,vy,vz,scu,n,s1,s2,t)
      do 120 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c zero out local accumulator
      nn = lxv*(my + 1)
      do 10 i = 1, nn
      scu(1,i) = 0.0
      scu(2,i) = 0.0
      scu(3,i) = 0.0
   10 continue
      ipp = npp/npblk
c outer loop over number of full blocks
      do 60 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 20 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      n(j) = nn - noff + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
      s2(j,1) = ppart(j+joff,3,k)
      s2(j,2) = ppart(j+joff,4,k)
      s2(j,3) = ppart(j+joff,5,k)
   20 continue
c deposit current
      do 40 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
!dir$ ivdep
      do 30 i = 1, lvect
      if (i.gt.2) nn = mm
      scu(1,i+nn) = scu(1,i+nn) + vx*s1(j,i)
      scu(2,i+nn) = scu(2,i+nn) + vy*s1(j,i)
      scu(3,i+nn) = scu(3,i+nn) + vz*s1(j,i)
   30 continue
   40 continue
c advance position half a time-step
!dir$ novector
      do 50 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      vx = s2(j,1)
      vy = s2(j,2)
      dx = x + vx*dt
      dy = y + vy*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j+joff,3,k) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(j+joff,4,k) = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j+joff,3,k) = -vx
         endif
      endif
c set new position
      ppart(j+joff,1,k) = dx
      ppart(j+joff,2,k) = dy
   50 continue
   60 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 70 j = nps, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ppart(j,3,k)
      vy = ppart(j,4,k)
      vz = ppart(j,5,k)
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j,3,k) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(j,4,k) = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j,3,k) = -vx
         endif
      endif
c set new position
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
   70 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 90 j = 2, mm
      do 80 i = 2, nn
      cu(1,i+noff+nxv*(j+moff-1)) = cu(1,i+noff+nxv*(j+moff-1))         
     1 + scu(1,i+lxv*(j-1))
      cu(2,i+noff+nxv*(j+moff-1)) = cu(2,i+noff+nxv*(j+moff-1))         
     1 + scu(2,i+lxv*(j-1))
      cu(3,i+noff+nxv*(j+moff-1)) = cu(3,i+noff+nxv*(j+moff-1))         
     1 + scu(3,i+lxv*(j-1))
   80 continue
   90 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 100 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noff+nxv*moff) = cu(1,i+noff+nxv*moff) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noff+nxv*moff) = cu(2,i+noff+nxv*moff) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noff+nxv*moff) = cu(3,i+noff+nxv*moff) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noff+nxv*(mm+moff-1)) = cu(1,i+noff+nxv*(mm+moff-1))    
     1    + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noff+nxv*(mm+moff-1)) = cu(2,i+noff+nxv*(mm+moff-1))    
     1    + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noff+nxv*(mm+moff-1)) = cu(3,i+noff+nxv*(mm+moff-1))    
     1    + scu(3,i+lxv*(mm-1))
      endif
  100 continue
      nn = min(mx+1,nxv-noff)
      do 110 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noff+nxv*(j+moff-1)) = cu(1,1+noff+nxv*(j+moff-1))         
     1 + scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noff+nxv*(j+moff-1)) = cu(2,1+noff+nxv*(j+moff-1))         
     1 + scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noff+nxv*(j+moff-1)) = cu(3,1+noff+nxv*(j+moff-1))         
     1 + scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff+nxv*(j+moff-1)) = cu(1,nn+noff+nxv*(j+moff-1))    
     1    + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noff+nxv*(j+moff-1)) = cu(2,nn+noff+nxv*(j+moff-1))    
     1    + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noff+nxv*(j+moff-1)) = cu(3,nn+noff+nxv*(j+moff-1))    
     1    + scu(3,nn+lxv*(j-1))
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine VGJPPOSTF2LT(ppart,cu,kpic,ncl,ihole,qm,dt,nppmx,idimp,
     1nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c vectorizable/OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = x velocity of particle n in tile m
c ppart(n,4,m) = y velocity of particle n in tile m
c ppart(n,5,m) = z velocity of particle n in tile m
c cu(i,j,k) = ith component of current density at grid point j,k
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ntmax
      integer irc
      real qm, dt
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mxy1), cu(4,nxv*nyv)
      dimension kpic(mxy1), ncl(8,mxy1)
      dimension ihole(2,ntmax+1,mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noff, moff, npp, ipp, joff, nps
      integer i, j, k, m, ih, nh, nn, mm, lxv
      real dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
      real anx, any, edgelx, edgely, edgerx, edgery
      real scu
      dimension scu(4,MXV*MYV)
c     dimension scu(4,(mx+1)*(my+1))
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,2)
      lxv = mx + 1
      anx = real(nx)
      any = real(ny)
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,noff,moff,npp,ipp,joff,nps,nn,mm,ih,nh,x,y,dxp,  
!$OMP& dyp,amx,amy,dx,dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,scu,n,s1, 
!$OMP& s2,t)
      do 130 k = 1, mxy1
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
c zero out local accumulator
      nn = lxv*(my + 1)
      do 10 i = 1, nn
      scu(1,i) = 0.0
      scu(2,i) = 0.0
      scu(3,i) = 0.0
   10 continue
c clear counters
      do 20 j = 1, 8
      ncl(j,k) = 0
   20 continue
      ipp = npp/npblk
c outer loop over number of full blocks
      do 70 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 30 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      n(j) = nn - noff + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
      s2(j,1) = ppart(j+joff,3,k)
      s2(j,2) = ppart(j+joff,4,k)
      s2(j,3) = ppart(j+joff,5,k)
   30 continue
c deposit current
      do 50 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
!dir$ ivdep
      do 40 i = 1, lvect
      if (i.gt.2) nn = mm
      scu(1,i+nn) = scu(1,i+nn) + vx*s1(j,i)
      scu(2,i+nn) = scu(2,i+nn) + vy*s1(j,i)
      scu(3,i+nn) = scu(3,i+nn) + vz*s1(j,i)
   40 continue
   50 continue
c advance position half a time-step
!dir$ novector
      do 60 j = 1, npblk
      dx = t(j,1) + s2(j,1)*dt
      dy = t(j,2) + s2(j,2)*dt
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
      ppart(j+joff,1,k) = dx
      ppart(j+joff,2,k) = dy
c increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j + joff
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
   60 continue
   70 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 80 j = nps, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ppart(j,3,k)
      vy = ppart(j,4,k)
      vz = ppart(j,5,k)
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
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
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
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
   80 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 100 j = 2, mm
      do 90 i = 2, nn
      cu(1,i+noff+nxv*(j+moff-1)) = cu(1,i+noff+nxv*(j+moff-1))         
     1 + scu(1,i+lxv*(j-1))
      cu(2,i+noff+nxv*(j+moff-1)) = cu(2,i+noff+nxv*(j+moff-1))         
     1 + scu(2,i+lxv*(j-1))
      cu(3,i+noff+nxv*(j+moff-1)) = cu(3,i+noff+nxv*(j+moff-1))         
     1 + scu(3,i+lxv*(j-1))
   90 continue
  100 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 110 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noff+nxv*moff) = cu(1,i+noff+nxv*moff) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noff+nxv*moff) = cu(2,i+noff+nxv*moff) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noff+nxv*moff) = cu(3,i+noff+nxv*moff) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noff+nxv*(mm+moff-1)) = cu(1,i+noff+nxv*(mm+moff-1))    
     1    + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noff+nxv*(mm+moff-1)) = cu(2,i+noff+nxv*(mm+moff-1))    
     1    + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noff+nxv*(mm+moff-1)) = cu(3,i+noff+nxv*(mm+moff-1))    
     1    + scu(3,i+lxv*(mm-1))
      endif
  110 continue
      nn = min(mx+1,nxv-noff)
      do 120 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noff+nxv*(j+moff-1)) = cu(1,1+noff+nxv*(j+moff-1))         
     1 + scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noff+nxv*(j+moff-1)) = cu(2,1+noff+nxv*(j+moff-1))         
     1 + scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noff+nxv*(j+moff-1)) = cu(3,1+noff+nxv*(j+moff-1))         
     1 + scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff+nxv*(j+moff-1)) = cu(1,nn+noff+nxv*(j+moff-1))    
     1    + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noff+nxv*(j+moff-1)) = cu(2,nn+noff+nxv*(j+moff-1))    
     1    + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noff+nxv*(j+moff-1)) = cu(3,nn+noff+nxv*(j+moff-1))    
     1    + scu(3,nn+lxv*(j-1))
      endif
  120 continue
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine VGRJPPOST2LT(ppart,cu,kpic,qm,dt,ci,nppmx,idimp,nx,ny, 
     1mx,my,nxv,nyv,mx1,mxy1,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c vectorizable/OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = x momentum of particle n in tile m
c ppart(n,4,m) = y momentum of particle n in tile m
c ppart(n,5,m) = z momentum of particle n in tile m
c cu(i,j,k) = ith component of current density at grid point j,k
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ipbc
      real qm, dt, ci
      real ppart, cu
      integer kpic
      dimension ppart(nppmx,idimp,mxy1), cu(4,nxv*nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noff, moff, npp, ipp, joff, nps
      integer i, j, k, m, nn, mm, lxv
      real ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami
      real scu
      dimension scu(4,MXV*MYV)
c     dimension scu(4,(mx+1)*(my+1))
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,4)
      lxv = mx + 1
      ci2 = ci*ci
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
!$OMP& PRIVATE(i,j,k,m,noff,moff,npp,ipp,joff,nps,nn,mm,x,y,dxp,dyp,amx,
!$OMP& amy,dx,dy,vx,vy,vz,ux,uy,uz,p2,gami,scu,n,s1,s2,t)
      do 120 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c zero out local accumulator
      nn = lxv*(my + 1)
      do 10 i = 1, nn
      scu(1,i) = 0.0
      scu(2,i) = 0.0
      scu(3,i) = 0.0
   10 continue
      ipp = npp/npblk
c outer loop over number of full blocks
      do 60 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 20 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      n(j) = nn - noff + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
c find inverse gamma
      ux = ppart(j+joff,3,k)
      uy = ppart(j+joff,4,k)
      uz = ppart(j+joff,5,k)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      s2(j,1) = ux*gami
      s2(j,2) = uy*gami
      s2(j,3) = uz*gami
      t(j,3) = ux
      t(j,4) = uy
   20 continue
c deposit current
      do 40 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
!dir$ ivdep
      do 30 i = 1, lvect
      if (i.gt.2) nn = mm
      scu(1,i+nn) = scu(1,i+nn) + vx*s1(j,i)
      scu(2,i+nn) = scu(2,i+nn) + vy*s1(j,i)
      scu(3,i+nn) = scu(3,i+nn) + vz*s1(j,i)
   30 continue
   40 continue
c advance position half a time-step
!dir$ novector
      do 50 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      vx = s2(j,1)
      vy = s2(j,2)
      ux = t(j,3)
      uy = t(j,4)
      dx = x + vx*dt
      dy = y + vy*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j+joff,3,k) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(j+joff,4,k) = -uy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j+joff,3,k) = -ux
         endif
      endif
c set new position
      ppart(j+joff,1,k) = dx
      ppart(j+joff,2,k) = dy
   50 continue
   60 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 70 j = nps, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
c find inverse gamma
      ux = ppart(j,3,k)
      uy = ppart(j,4,k)
      uz = ppart(j,5,k)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ux*gami
      vy = uy*gami
      vz = uz*gami
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j,3,k) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            ppart(j,4,k) = -uy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            ppart(j,3,k) = -ux
         endif
      endif
c set new position
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
   70 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 90 j = 2, mm
      do 80 i = 2, nn
      cu(1,i+noff+nxv*(j+moff-1)) = cu(1,i+noff+nxv*(j+moff-1))         
     1 + scu(1,i+lxv*(j-1))
      cu(2,i+noff+nxv*(j+moff-1)) = cu(2,i+noff+nxv*(j+moff-1))         
     1 + scu(2,i+lxv*(j-1))
      cu(3,i+noff+nxv*(j+moff-1)) = cu(3,i+noff+nxv*(j+moff-1))         
     1 + scu(3,i+lxv*(j-1))
   80 continue
   90 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 100 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noff+nxv*moff) = cu(1,i+noff+nxv*moff) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noff+nxv*moff) = cu(2,i+noff+nxv*moff) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noff+nxv*moff) = cu(3,i+noff+nxv*moff) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noff+nxv*(mm+moff-1)) = cu(1,i+noff+nxv*(mm+moff-1))    
     1    + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noff+nxv*(mm+moff-1)) = cu(2,i+noff+nxv*(mm+moff-1))    
     1    + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noff+nxv*(mm+moff-1)) = cu(3,i+noff+nxv*(mm+moff-1))    
     1    + scu(3,i+lxv*(mm-1))
      endif
  100 continue
      nn = min(mx+1,nxv-noff)
      do 110 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noff+nxv*(j+moff-1)) = cu(1,1+noff+nxv*(j+moff-1))         
     1 + scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noff+nxv*(j+moff-1)) = cu(2,1+noff+nxv*(j+moff-1))         
     1 + scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noff+nxv*(j+moff-1)) = cu(3,1+noff+nxv*(j+moff-1))         
     1 + scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff+nxv*(j+moff-1)) = cu(1,nn+noff+nxv*(j+moff-1))    
     1    + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noff+nxv*(j+moff-1)) = cu(2,nn+noff+nxv*(j+moff-1))    
     1    + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noff+nxv*(j+moff-1)) = cu(3,nn+noff+nxv*(j+moff-1))    
     1    + scu(3,nn+lxv*(j-1))
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine VGRJPPOSTF2LT(ppart,cu,kpic,ncl,ihole,qm,dt,ci,nppmx,  
     1idimp,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c vectorizable/OpenMP version using guard cells
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = x momentum of particle n in tile m
c ppart(n,4,m) = y momentum of particle n in tile m
c ppart(n,5,m) = z momentum of particle n in tile m
c cu(i,j,k) = ith component of current density at grid point j,k
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer nppmx, idimp, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ntmax
      integer irc
      real qm, dt, ci
      real ppart, cu
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mxy1), cu(4,nxv*nyv)
      dimension kpic(mxy1), ncl(8,mxy1)
      dimension ihole(2,ntmax+1,mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer noff, moff, npp, ipp, joff, nps
      integer i, j, k, m, ih, nh, nn, mm, lxv
      real ci2, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, p2, gami
      real anx, any, edgelx, edgely, edgerx, edgery
      real scu
      dimension scu(4,MXV*MYV)
c     dimension scu(4,(mx+1)*(my+1))
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,2)
      lxv = mx + 1
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,noff,moff,npp,ipp,joff,nps,nn,mm,ih,nh,x,y,dxp,  
!$OMP& dyp,amx,amy,dx,dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,p2,gami,  
!$OMP& scu,n,s1,s2,t)
      do 130 k = 1, mxy1
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
c zero out local accumulator
      nn = lxv*(my + 1)
      do 10 i = 1, nn
      scu(1,i) = 0.0
      scu(2,i) = 0.0
      scu(3,i) = 0.0
   10 continue
c clear counters
      do 20 j = 1, 8
      ncl(j,k) = 0
   20 continue
      ipp = npp/npblk
c outer loop over number of full blocks
      do 70 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 30 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      n(j) = nn - noff + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
c find inverse gamma
      vx = ppart(j+joff,3,k)
      vy = ppart(j+joff,4,k)
      vz = ppart(j+joff,5,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      s2(j,1) = vx*gami
      s2(j,2) = vy*gami
      s2(j,3) = vz*gami
   30 continue
c deposit current
      do 50 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
!dir$ ivdep
      do 40 i = 1, lvect
      if (i.gt.2) nn = mm
      scu(1,i+nn) = scu(1,i+nn) + vx*s1(j,i)
      scu(2,i+nn) = scu(2,i+nn) + vy*s1(j,i)
      scu(3,i+nn) = scu(3,i+nn) + vz*s1(j,i)
   40 continue
   50 continue
c advance position half a time-step
!dir$ novector
      do 60 j = 1, npblk
      dx = t(j,1) + s2(j,1)*dt
      dy = t(j,2) + s2(j,2)*dt
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
      ppart(j+joff,1,k) = dx
      ppart(j+joff,2,k) = dy
c increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j + joff
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
   60 continue
   70 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 80 j = nps, npp
c find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
c find inverse gamma
      vx = ppart(j,3,k)
      vy = ppart(j,4,k)
      vz = ppart(j,5,k)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn - noff + 1 + lxv*(mm - moff)
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      scu(1,nn) = scu(1,nn) + vx*dx
      scu(2,nn) = scu(2,nn) + vy*dx
      scu(3,nn) = scu(3,nn) + vz*dx
      dx = amx*dyp
      scu(1,nn+1) = scu(1,nn+1) + vx*dy
      scu(2,nn+1) = scu(2,nn+1) + vy*dy
      scu(3,nn+1) = scu(3,nn+1) + vz*dy
      dy = dxp*dyp
      scu(1,nn+lxv) = scu(1,nn+lxv) + vx*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + vy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + vz*dx
      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + vx*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + vy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + vz*dy
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
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
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
   80 continue
c deposit current to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 100 j = 2, mm
      do 90 i = 2, nn
      cu(1,i+noff+nxv*(j+moff-1)) = cu(1,i+noff+nxv*(j+moff-1))         
     1 + scu(1,i+lxv*(j-1))
      cu(2,i+noff+nxv*(j+moff-1)) = cu(2,i+noff+nxv*(j+moff-1))         
     1 + scu(2,i+lxv*(j-1))
      cu(3,i+noff+nxv*(j+moff-1)) = cu(3,i+noff+nxv*(j+moff-1))         
     1 + scu(3,i+lxv*(j-1))
   90 continue
  100 continue
c deposit current to edge points in global array
      mm = min(my+1,nyv-moff)
      do 110 i = 2, nn
!$OMP ATOMIC
      cu(1,i+noff+nxv*moff) = cu(1,i+noff+nxv*moff) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noff+nxv*moff) = cu(2,i+noff+nxv*moff) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noff+nxv*moff) = cu(3,i+noff+nxv*moff) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
         cu(1,i+noff+nxv*(mm+moff-1)) = cu(1,i+noff+nxv*(mm+moff-1))    
     1    + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(2,i+noff+nxv*(mm+moff-1)) = cu(2,i+noff+nxv*(mm+moff-1))    
     1    + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
         cu(3,i+noff+nxv*(mm+moff-1)) = cu(3,i+noff+nxv*(mm+moff-1))    
     1    + scu(3,i+lxv*(mm-1))
      endif
  110 continue
      nn = min(mx+1,nxv-noff)
      do 120 j = 1, mm
!$OMP ATOMIC
      cu(1,1+noff+nxv*(j+moff-1)) = cu(1,1+noff+nxv*(j+moff-1))         
     1 + scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noff+nxv*(j+moff-1)) = cu(2,1+noff+nxv*(j+moff-1))         
     1 + scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noff+nxv*(j+moff-1)) = cu(3,1+noff+nxv*(j+moff-1))         
     1 + scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         cu(1,nn+noff+nxv*(j+moff-1)) = cu(1,nn+noff+nxv*(j+moff-1))    
     1    + scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(2,nn+noff+nxv*(j+moff-1)) = cu(2,nn+noff+nxv*(j+moff-1))    
     1    + scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
         cu(3,nn+noff+nxv*(j+moff-1)) = cu(3,nn+noff+nxv*(j+moff-1))    
     1    + scu(3,nn+lxv*(j-1))
      endif
  120 continue
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine viscan2(isdata,mb,nths)
c performs vectorizable prefix reduction of integer data
c using binary tree method.
      implicit none
      integer nths
      integer isdata, mb
      dimension isdata(nths), mb(nths/2)
c local data
      integer j, kxs, lb, ns
      ns = nths/2
      do 10 j = 1, ns
      mb(j) = j - 1
   10 continue
      kxs = 1
   20 if (kxs.lt.nths) then
!dir$ ivdep
         do 30 j = 1, ns
         lb = kxs*mb(j)
         if ((j+lb+kxs).le.nths) then
            isdata(j+lb+kxs) = isdata(j+lb+kxs) + isdata(2*lb+kxs)
         endif
         mb(j) = mb(j)/2
   30    continue
         kxs = kxs + kxs
         go to 20
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPORDER2LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx, 
     1ny,mx,my,mx1,my1,npbmx,ntmax,irc)
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppbuff(n,i,k) = i co-ordinate of particle n in tile k
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
      dimension ppart(nppmx,idimp,mx1*my1), ppbuff(npbmx,idimp,mx1*my1)
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
      dx = ppart(j,1,k)
      dy = ppart(j,2,k)
c find particles going out of bounds
      ist = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(j,1,k) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(j,1,k) = dx
         else
            ist = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(j,2,k) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(j,2,k) = dy
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
         ppbuff(ii,i,k) = ppart(j1,i,k)
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
!$OMP& PRIVATE(i,j,k,ii,kk,npp,kx,ky,kl,kr,kxl,kxr,ih,nh,nn,ncoff,ist,j1
!$OMP& ,j2,ip,ks)
      do 140 k = 1, mxy1
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
         ppart(j1,i,k) = ppbuff(j+ncoff,i,ks(ii))
   80    continue
      else
         ist = 1
      endif
   90 continue
  100 continue
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
c holes with locations great than npp-ip do not need to be filled
      if (ih.lt.nh) then
         ip = nh - ih
         ii = nh + 1
         nn = ihole(1,ii,k)
         ih = ih + 2
         j2 = ihole(1,ih,k)
c move particles from end into remaining holes
c holes are processed in increasing order
         do 130 j = 1, ip
         j1 = npp - j + 1
         if (j1.eq.nn) then
            ii = ii - 1
            nn = ihole(1,ii,k)
         else
            do 120 i = 1, idimp
            ppart(j2,i,k) = ppart(j1,i,k)
  120       continue
            ih = ih + 1
            j2 = ihole(1,ih,k)
         endif
  130    continue
         npp = npp - ip
      endif
      kpic(k) = npp
  140 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPORDERF2LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,mx1
     1,my1,npbmx,ntmax,irc)
c this subroutine sorts particles by x,y grid in tiles of mx, my
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 2D linear memory.
c the algorithm has 2 steps.  first, a prefix scan of ncl is performed
c and departing particles are buffered in ppbuff in direction order.
c then we copy the incoming particles from other tiles into ppart.
c it assumes that the number, location, and destination of particles 
c leaving a tile have been previously stored in ncl and ihole by the
c GPPUSHF2LT subroutine.
c input: all except ppbuff, irc
c output: ppart, ppbuff, kpic, ncl, irc
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppbuff(n,i,k) = i co-ordinate of particle n in tile k
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
      dimension ppart(nppmx,idimp,mx1*my1), ppbuff(npbmx,idimp,mx1*my1)
      dimension kpic(mx1*my1), ncl(8,mx1*my1)
      dimension ihole(2,ntmax+1,mx1*my1)
c local data
      integer mxy1, npp, ncoff
      integer i, j, k, ii, kx, ky, ih, nh, ist, nn, isum
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
         ppbuff(ii,i,k) = ppart(j1,i,k)
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
!$OMP& PRIVATE(i,j,k,ii,kk,npp,kx,ky,kl,kr,kxl,kxr,ih,nh,nn,ncoff,ist,j1
!$OMP& ,j2,ip,ks)
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
         ppart(j1,i,k) = ppbuff(j+ncoff,i,ks(ii))
   50    continue
      else
         ist = 1
      endif
   60 continue
   70 continue
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
c holes with locations great than npp-ip do not need to be filled
      if (ih.lt.nh) then
         ip = nh - ih
         ii = nh + 1
         nn = ihole(1,ii,k)
         ih = ih + 2
         j2 = ihole(1,ih,k)
c move particles from end into remaining holes
c holes are processed in increasing order
         do 90 j = 1, ip
         j1 = npp - j + 1
         if (j1.eq.nn) then
            ii = ii - 1
            nn = ihole(1,ii,k)
         else
            do 80 i = 1, idimp
            ppart(j2,i,k) = ppart(j1,i,k)
   80       continue
            ih = ih + 1
            j2 = ihole(1,ih,k)
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
      subroutine VPPORDER2LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx,
     1ny,mx,my,mx1,my1,npbmx,ntmax,irc)
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
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppbuff(n,i,k) = i co-ordinate of particle n in tile k
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
      dimension ppart(nppmx,idimp,mx1*my1), ppbuff(npbmx,idimp,mx1*my1)
      dimension kpic(mx1*my1), ncl(8,mx1*my1)
      dimension ihole(2,ntmax+1,mx1*my1)
c local data
      integer npblk
      parameter(npblk=16)
      integer mxy1, noff, moff, npp, ipp, joff, nps, ncoff
      integer i, j, k, m, ii, kx, ky, ih, nh, ist, nn, mm, in
      integer ip, j1, j2, kxl, kxr, kk, kl, kr, lb, kxs
      real anx, any, edgelx, edgely, edgerx, edgery, dx, dy
      integer sncl, ks
      dimension sncl(8), ks(8)
c scratch arrays
      integer n
      dimension n(npblk,3)
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
      dx = ppart(j,1,k)
      dy = ppart(j,2,k)
c find particles going out of bounds
      ist = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(j,1,k) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(j,1,k) = dx
         else
            ist = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(j,2,k) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(j,2,k) = dy
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
!$OMP& PRIVATE(i,j,k,m,kxs,lb,ist,nh,ip,ipp,nps,joff,j1,ii,sncl,ks,n)
      do 160 k = 1, mxy1
c find address offset for ordered ppbuff array
      do 40 j = 1, 8
      sncl(j) = ncl(j,k)
      ks(j) = j - 1
   40 continue
      kxs = 1
   50 if (kxs.lt.8) then
!dir$ ivdep
         do 60 j = 1, 4
         lb = kxs*ks(j)
         sncl(j+lb+kxs) = sncl(j+lb+kxs) + sncl(2*lb+kxs)
         ks(j) = ks(j)/2
   60    continue
         kxs = kxs + kxs
         go to 50
      endif
      do 70 j = 1, 8
      sncl(j) = sncl(j) - ncl(j,k)
   70 continue
      nh = ihole(1,1,k)
      ip = 0
c buffer particles that are leaving tile, in direction order
c loop over particles leaving tile
      ipp = nh/npblk
c outer loop over number of full blocks
      do 120 m = 1, ipp
      joff = npblk*(m - 1) + 1
c inner loop over particles in block
      do 80 j = 1, npblk
      n(j,1) = ihole(1,j+joff,k)
      n(j,2) = ihole(2,j+joff,k)
   80 continue
c calculate offsets
      do 90 j = 1, npblk
      ist = n(j,2)
      ii = sncl(ist) + 1
      n(j,2) = ii
      sncl(ist) = ii
   90 continue
c buffer particles that are leaving tile, in direction order
      do 110 i = 1, idimp
      do 100 j = 1, npblk
      j1 = n(j,1)
      ii = n(j,2)
      if (ii.le.npbmx) then
         ppbuff(ii,i,k) = ppart(j1,i,k)
      else
         ip = 1
      endif
  100 continue
  110 continue
  120 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 140 j = nps, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = sncl(ist) + 1
      if (ii.le.npbmx) then
         do 130 i = 1, idimp
         ppbuff(ii,i,k) = ppart(j1,i,k)
  130    continue
      else
         ip = 1
      endif
      sncl(ist) = ii
  140 continue
      do 150 j = 1, 8
      ncl(j,k) = sncl(j)
  150 continue
c set error
      if (ip.gt.0) irc = ncl(8,k)
  160 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ii,kk,in,npp,ipp,joff,nps,kx,ky,kl,kr,kxl,kxr,ih,
!$OMP& nh,nn,mm,ncoff,ist,j1,j2,ip,ks,n)
      do 310 k = 1, mxy1
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
      do 230 ii = 1, 8
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
c ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
c loop over particles coming from direction ii
      ipp = ip/npblk
c outer loop over number of full blocks
      do 200 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 170 j = 1, npblk
c insert incoming particles into holes
      if ((j+ih).le.nh) then
         j1 = ihole(1,j+ih+1,k)
c place overflow at end of array
      else
         j1 = npp + j + ih - nh
      endif
      n(j,1) = j1
  170 continue
      do 190 i = 1, idimp
      do 180 j = 1, npblk
      j1 = n(j,1)
      if (j1.le.nppmx) then
         ppart(j1,i,k) = ppbuff(j+joff+ncoff,i,ks(ii))
      else
         ist = 1
      endif
  180 continue
  190 continue
      ih = ih + npblk
  200 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 220 j = nps, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,k)
c place overflow at end of array
      else
         j1 = npp + ih - nh
      endif
      if (j1.le.nppmx) then
         do 210 i = 1, idimp
         ppart(j1,i,k) = ppbuff(j+ncoff,i,ks(ii))
  210    continue
      else
         ist = 1
      endif
  220 continue
  230 continue
      if (ih > nh) npp = npp + ih - nh
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
c holes with locations great than npp-ip do not need to be filled
      if (ih.lt.nh) then
         ip = nh - ih
c move particles from end into remaining holes
c holes are processed in increasing order
         ii = nh + 1
         ipp = ip/npblk
c outer loop over number of full blocks
         do 280 m = 1, ipp
         joff = npblk*(m - 1)
c inner loop over particles in block
         do 240 j = 1, npblk
         n(j,2) = ihole(1,ih+j+1,k)
         n(j,3) = ihole(1,ii-j+1,k)
  240    continue
         in = 1
         mm = 1
         nn = n(in,3)
         do 250 j = 1, npblk
         j1 = npp - j - joff + 1
         n(j,1) = n(mm,2)
         if (j1.eq.nn) then
            in = in + 1
            nn = n(in,3)
            n(j,1) = -1
         else
            mm = mm + 1
         endif
  250    continue
         do 270 i = 1, idimp
!dir$ ivdep
         do 260 j = 1, npblk
         j1 = npp - j - joff + 1
         j2 = n(j,1)
         if (j2.gt.0) ppart(j2,i,k) = ppart(j1,i,k)
  260    continue
  270    continue
         ii = ii - in + 1
         ih = ih + mm - 1
  280    continue
         nps = npblk*ipp + 1
         nn = ihole(1,ii,k)
         ih = ih + 2
         j2 = ihole(1,ih,k)
c loop over remaining particles
         do 300 j = nps, ip
         j1 = npp - j + 1
         if (j1.eq.nn) then
            ii = ii - 1
            nn = ihole(1,ii,k)
         else
            do 290 i = 1, idimp
            ppart(j2,i,k) = ppart(j1,i,k)
  290       continue
            ih = ih + 1
            j2 = ihole(1,ih,k)
         endif
  300    continue
         npp = npp - ip
      endif
      kpic(k) = npp
  310 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine VPPORDERF2LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,  
     1mx1,my1,npbmx,ntmax,irc)
c this subroutine sorts particles by x,y grid in tiles of mx, my
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 2D linear memory.
c the algorithm has 2 steps.  first, a prefix scan of ncl is performed
c and departing particles are buffered in ppbuff in direction order.
c then we copy the incoming particles from other tiles into ppart.
c it assumes that the number, location, and destination of particles 
c leaving a tile have been previously stored in ncl and ihole by the
c GPPUSHF2LT subroutine.
c input: all except ppbuff, irc
c output: ppart, ppbuff, kpic, ncl, irc
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppbuff(n,i,k) = i co-ordinate of particle n in tile k
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
      dimension ppart(nppmx,idimp,mx1*my1), ppbuff(npbmx,idimp,mx1*my1)
      dimension kpic(mx1*my1), ncl(8,mx1*my1)
      dimension ihole(2,ntmax+1,mx1*my1)
c local data
      integer npblk
      parameter(npblk=16)
      integer mxy1, npp, ncoff
      integer i, j, k, ii, kx, ky, ih, nh, ist, nn, mm, in
      integer ip, j1, j2, kxl, kxr, kk, kl, kr
      integer lb, kxs, m, ipp, nps, joff
      integer sncl, ks
      dimension sncl(8), ks(8)
c scratch arrays
      integer n
      dimension n(npblk,3)
      mxy1 = mx1*my1
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,kxs,lb,ist,nh,ip,ipp,nps,joff,j1,ii,sncl,ks,n)
      do 130 k = 1, mxy1
c find address offset for ordered ppbuff array
      do 10 j = 1, 8
      sncl(j) = ncl(j,k)
      ks(j) = j - 1
   10 continue
      kxs = 1
   20 if (kxs.lt.8) then
!dir$ ivdep
         do 30 j = 1, 4
         lb = kxs*ks(j)
         sncl(j+lb+kxs) = sncl(j+lb+kxs) + sncl(2*lb+kxs)
         ks(j) = ks(j)/2
   30    continue
         kxs = kxs + kxs
         go to 20
      endif
      do 40 j = 1, 8
      sncl(j) = sncl(j) - ncl(j,k)
   40 continue
      nh = ihole(1,1,k)
      ip = 0
c buffer particles that are leaving tile, in direction order
c loop over particles leaving tile
      ipp = nh/npblk
c outer loop over number of full blocks
      do 90 m = 1, ipp
      joff = npblk*(m - 1) + 1
c inner loop over particles in block
      do 50 j = 1, npblk
      n(j,1) = ihole(1,j+joff,k)
      n(j,2) = ihole(2,j+joff,k)
   50 continue
c calculate offsets
      do 60 j = 1, npblk
      ist = n(j,2)
      ii = sncl(ist) + 1
      n(j,2) = ii
      sncl(ist) = ii
   60 continue
c buffer particles that are leaving tile, in direction order
      do 80 i = 1, idimp
      do 70 j = 1, npblk
      j1 = n(j,1)
      ii = n(j,2)
      if (ii.le.npbmx) then
         ppbuff(ii,i,k) = ppart(j1,i,k)
      else
         ip = 1
      endif
   70 continue
   80 continue
   90 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 110 j = nps, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = sncl(ist) + 1
      if (ii.le.npbmx) then
         do 100 i = 1, idimp
         ppbuff(ii,i,k) = ppart(j1,i,k)
  100    continue
      else
         ip = 1
      endif
      sncl(ist) = ii
  110 continue
      do 120 j = 1, 8
      ncl(j,k) = sncl(j)
  120 continue
c set error
      if (ip.gt.0) irc = ncl(8,k)
  130 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ii,kk,in,npp,ipp,joff,nps,kx,ky,kl,kr,kxl,kxr,ih,
!$OMP& nh,nn,mm,ncoff,ist,j1,j2,ip,ks,n)
      do 280 k = 1, mxy1
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
      do 200 ii = 1, 8
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
c ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
c loop over particles coming from direction ii
      ipp = ip/npblk
c outer loop over number of full blocks
      do 170 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 140 j = 1, npblk
c insert incoming particles into holes
      if ((j+ih).le.nh) then
         j1 = ihole(1,j+ih+1,k)
c place overflow at end of array
      else
         j1 = npp + j + ih - nh
      endif
      n(j,1) = j1
  140 continue
      do 160 i = 1, idimp
      do 150 j = 1, npblk
      j1 = n(j,1)
      if (j1.le.nppmx) then
         ppart(j1,i,k) = ppbuff(j+joff+ncoff,i,ks(ii))
      else
         ist = 1
      endif
  150 continue
  160 continue
      ih = ih + npblk
  170 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 190 j = nps, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,k)
c place overflow at end of array
      else
         j1 = npp + ih - nh
      endif
      if (j1.le.nppmx) then
         do 180 i = 1, idimp
         ppart(j1,i,k) = ppbuff(j+ncoff,i,ks(ii))
  180    continue
      else
         ist = 1
      endif
  190 continue
  200 continue
      if (ih > nh) npp = npp + ih - nh
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
c holes with locations great than npp-ip do not need to be filled
      if (ih.lt.nh) then
         ip = nh - ih
c move particles from end into remaining holes
c holes are processed in increasing order
         ii = nh + 1
c loop over particles coming from direction ii
         ipp = ip/npblk
c outer loop over number of full blocks
         do 250 m = 1, ipp
         joff = npblk*(m - 1)
c inner loop over particles in block
         do 210 j = 1, npblk
         n(j,2) = ihole(1,ih+j+1,k)
         n(j,3) = ihole(1,ii-j+1,k)
  210    continue
         in = 1
         mm = 1
         nn = n(in,3)
         do 220 j = 1, npblk
         j1 = npp - j - joff + 1
         n(j,1) = n(mm,2)
         if (j1.eq.nn) then
            in = in + 1
            nn = n(in,3)
            n(j,1) = -1
         else
            mm = mm + 1
         endif
  220    continue
         do 240 i = 1, idimp
!dir$ ivdep
         do 230 j = 1, npblk
         j1 = npp - j - joff + 1
         j2 = n(j,1)
         if (j2.gt.0) ppart(j2,i,k) = ppart(j1,i,k)
  230    continue
  240    continue
         ii = ii - in + 1
         ih = ih + mm - 1
  250    continue
         nps = npblk*ipp + 1
         nn = ihole(1,ii,k)
         ih = ih + 2
         j2 = ihole(1,ih,k)
c loop over remaining particles
         do 270 j = nps, ip
         j1 = npp - j + 1
         if (j1.eq.nn) then
            ii = ii - 1
            nn = ihole(1,ii,k)
         else
            do 260 i = 1, idimp
            ppart(j2,i,k) = ppart(j1,i,k)
  260       continue
            ih = ih + 1
            j2 = ihole(1,ih,k)
         endif
  270    continue
         npp = npp - ip
      endif
      kpic(k) = npp
  280 continue
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
      dimension bxy(4,nxe,nye)
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
      dimension cu(4,nxe,nye)
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
      subroutine VMPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv, 
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
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, we
      complex q, fxy, ffc
      dimension q(nxvh,nyv), fxy(4,nxvh,nyv)
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
!dir$ ivdep
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
      at1 = at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
      wp = wp + dble(at1)
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
      at1 = at1*(q(1,k)*conjg(q(1,k)))
      wp = wp + dble(at1)
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
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
      at1 = at1*(q(j,1)*conjg(q(j,1)))
      wp = wp + dble(at1)
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
c nxvh = second dimension of current array, must be >= nxh
c nyv = third dimension of current array, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex cu
      dimension cu(4,nxvh,nyv)
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
!dir$ ivdep
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
!dir$ ivdep
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
      subroutine VMIBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with periodic boundary conditions.
c input: cu,ffc,ci,nx,ny,nxv,nyhd, output: bxy,wm
c approximate flop count is: 90*nxc*nyc + 40*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = bx(ky=pi) = by(ky=pi) = bz(ky=pi) 
c = 0, and bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c bxy(i,j,k) = i component of complex magnetic field
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
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real ci, wm
      complex cu, bxy, ffc
      dimension cu(4,nxvh,nyv), bxy(4,nxvh,nyv)
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
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      sum1 = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,at1,at2,at3,zt1,zt2,zt3,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      wp = 0.0d0
!dir$ ivdep
      do 10 j = 2, nxh
      at1 = ci2*real(ffc(j,k))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      at1 = at1*aimag(ffc(j,k))
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
      at1 = at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j,k))
     1   + cu(3,j,k)*conjg(cu(3,j,k)) + cu(1,j,k1)*conjg(cu(1,j,k1))    
     2   + cu(2,j,k1)*conjg(cu(2,j,k1)) + cu(3,j,k1)*conjg(cu(3,j,k1)))
      wp = wp + dble(at1)
   10 continue
c mode numbers kx = 0, nx/2
      at1 = ci2*real(ffc(1,k))
      at3 = dny*real(k - 1)*at1
      at1 = at1*aimag(ffc(1,k))
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bxy(1,1,k) = at3*zt1
      bxy(2,1,k) = zero
      bxy(3,1,k) = -at3*zt3
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      at1 = at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1,k))
     1   + cu(3,1,k)*conjg(cu(3,1,k)))
      wp = wp + dble(at1)
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
cdir$ ivdep
      do 30 j = 2, nxh
      at1 = ci2*real(ffc(j,1))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffc(j,1))
      zt1 = cmplx(-aimag(cu(3,j,1)),real(cu(3,j,1)))
      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
      bxy(1,j,1) = zero
      bxy(2,j,1) = -at2*zt1
      bxy(3,j,1) = at2*zt2
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      at1 = at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j,1))
     1   + cu(3,j,1)*conjg(cu(3,j,1)))
      wp = wp + dble(at1)
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
      subroutine VMMAXWEL2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,   
     1nxhd,nyhd)
c this subroutine solves 2-1/2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions.
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
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
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c and similarly for bx, by, bz.
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffc(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprocal of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real ci, dt, wf, wm
      complex exy, bxy, cu, ffc
      dimension exy(4,nxvh,nyv), bxy(4,nxvh,nyv), cu(4,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dth, c2, cdt, affp, anorm, dkx, dky, afdt, adt
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      real at1
      double precision wp, ws, sum1, sum2
      if (ci.le.0.0) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dth = 0.5*dt
      c2 = 1.0/(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1,1))
      adt = affp*dt
      zero = cmplx(0.0,0.0)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      sum1 = 0.0d0
      sum2 = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,dkx,afdt,at1,zt1,zt2,zt3,zt4,zt5,  
!$OMP& zt6,zt7,zt8,zt9,ws,wp)
!$OMP& REDUCTION(+:sum1,sum2)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      ws = 0.0d0
      wp = 0.0d0
!dir$ ivdep
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j,k))
c update magnetic field half time step, ky > 0
      zt1 = cmplx(-aimag(exy(3,j,k)),real(exy(3,j,k)))
      zt2 = cmplx(-aimag(exy(2,j,k)),real(exy(2,j,k)))
      zt3 = cmplx(-aimag(exy(1,j,k)),real(exy(1,j,k)))
      zt4 = bxy(1,j,k) - dth*(dky*zt1)
      zt5 = bxy(2,j,k) + dth*(dkx*zt1)
      zt6 = bxy(3,j,k) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k) + cdt*(dky*zt1) - afdt*cu(1,j,k)
      zt8 = exy(2,j,k) - cdt*(dkx*zt1) - afdt*cu(2,j,k)
      zt9 = exy(3,j,k) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,j,k)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,j,k) = zt7
      exy(2,j,k) = zt8
      exy(3,j,k) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 - dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
      bxy(1,j,k) = zt4
      bxy(2,j,k) = zt5
      bxy(3,j,k) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
c update magnetic field half time step, ky < 0
      zt1 = cmplx(-aimag(exy(3,j,k1)),real(exy(3,j,k1)))
      zt2 = cmplx(-aimag(exy(2,j,k1)),real(exy(2,j,k1)))
      zt3 = cmplx(-aimag(exy(1,j,k1)),real(exy(1,j,k1)))
      zt4 = bxy(1,j,k1) + dth*(dky*zt1)
      zt5 = bxy(2,j,k1) + dth*(dkx*zt1)
      zt6 = bxy(3,j,k1) - dth*(dkx*zt2 + dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k1) - cdt*(dky*zt1) - afdt*cu(1,j,k1)
      zt8 = exy(2,j,k1) - cdt*(dkx*zt1) - afdt*cu(2,j,k1)
      zt9 = exy(3,j,k1) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,j,k1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,j,k1) = zt7
      exy(2,j,k1) = zt8
      exy(3,j,k1) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 + dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
      bxy(1,j,k1) = zt4
      bxy(2,j,k1) = zt5
      bxy(3,j,k1) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
   10 continue
c mode numbers kx = 0, nx/2
      afdt = adt*aimag(ffc(1,k))
c update magnetic field half time step
      zt1 = cmplx(-aimag(exy(3,1,k)),real(exy(3,1,k)))
      zt3 = cmplx(-aimag(exy(1,1,k)),real(exy(1,1,k)))
      zt4 = bxy(1,1,k) - dth*(dky*zt1)
      zt6 = bxy(3,1,k) + dth*(dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,1,k) + cdt*(dky*zt1) - afdt*cu(1,1,k)
      zt9 = exy(3,1,k) - cdt*(dky*zt3) - afdt*cu(3,1,k)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,1,k) = zt7
      exy(2,1,k) = zero
      exy(3,1,k) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 - dth*(dky*zt1)
      zt6 = zt6 + dth*(dky*zt3)
      bxy(1,1,k) = zt4
      bxy(2,1,k) = zero
      bxy(3,1,k) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
      wp = wp + dble(at1)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      sum1 = sum1 + ws
      sum2 = sum2 + wp
   20 continue
!$OMP END PARALLEL DO
      ws = 0.0d0
      wp = 0.0d0
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 30 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j,1))
c update magnetic field half time step
      zt1 = cmplx(-aimag(exy(3,j,1)),real(exy(3,j,1)))
      zt2 = cmplx(-aimag(exy(2,j,1)),real(exy(2,j,1)))
      zt5 = bxy(2,j,1) + dth*(dkx*zt1)
      zt6 = bxy(3,j,1) - dth*(dkx*zt2)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt8 = exy(2,j,1) - cdt*(dkx*zt1) - afdt*cu(2,j,1)
      zt9 = exy(3,j,1) + cdt*(dkx*zt2) - afdt*cu(3,j,1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      exy(1,j,1) = zero
      exy(2,j,1) = zt8
      exy(3,j,1) = zt9
      at1 = anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      bxy(1,j,1) = zero
      bxy(2,j,1) = zt5
      bxy(3,j,1) = zt6
      at1 = anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
   30 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      sum1 = sum1 + ws
      sum2 = sum2 + wp
      wf = real(nx*ny)*sum1
      wm = real(nx*ny)*c2*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine VMEMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign < 0
c includes additional smoothing
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      complex fxy, exy, ffc
      dimension fxy(4,nxvh,nyv), exy(4,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer j, k, nxh, nyh, ny2, k1
      real at1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
c add the fields
      if (isign.gt.0) then
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1)
         do 20 k = 2, nyh
         k1 = ny2 - k
!dir$ ivdep
         do 10 j = 1, nxh
         at1 = aimag(ffc(j,k))
         fxy(1,j,k) = fxy(1,j,k) + exy(1,j,k)*at1
         fxy(2,j,k) = fxy(2,j,k) + exy(2,j,k)*at1
         fxy(3,j,k) = fxy(3,j,k) + exy(3,j,k)*at1
         fxy(1,j,k1) = fxy(1,j,k1) + exy(1,j,k1)*at1
         fxy(2,j,k1) = fxy(2,j,k1) + exy(2,j,k1)*at1
         fxy(3,j,k1) = fxy(3,j,k1) + exy(3,j,k1)*at1
   10    continue
   20    continue
!$OMP END PARALLEL DO
         k1 = nyh + 1
!dir$ ivdep
         do 30 j = 1, nxh
         at1 = aimag(ffc(j,1))
         fxy(1,j,1) = fxy(1,j,1) + exy(1,j,1)*at1
         fxy(2,j,1) = fxy(2,j,1) + exy(2,j,1)*at1
         fxy(3,j,1) = fxy(3,j,1) + exy(3,j,1)*at1
         fxy(1,j,k1) = fxy(1,j,k1) + exy(1,j,k1)*at1
         fxy(2,j,k1) = fxy(2,j,k1) + exy(2,j,k1)*at1
         fxy(3,j,k1) = fxy(3,j,k1) + exy(3,j,k1)*at1
   30    continue
c copy the fields
      else if (isign.lt.0) then
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1)
         do 50 k = 2, nyh
         k1 = ny2 - k
!dir$ ivdep
         do 40 j = 1, nxh
         at1 = aimag(ffc(j,k))
         fxy(1,j,k) = exy(1,j,k)*at1
         fxy(2,j,k) = exy(2,j,k)*at1
         fxy(3,j,k) = exy(3,j,k)*at1
         fxy(1,j,k1) = exy(1,j,k1)*at1
         fxy(2,j,k1) = exy(2,j,k1)*at1
         fxy(3,j,k1) = exy(3,j,k1)*at1
   40    continue
   50    continue
!$OMP END PARALLEL DO
         k1 = nyh + 1
!dir$ ivdep
         do 60 j = 1, nxh
         at1 = aimag(ffc(j,1))
         fxy(1,j,1) = exy(1,j,1)*at1
         fxy(2,j,1) = exy(2,j,1)*at1
         fxy(3,j,1) = exy(3,j,1)*at1
         fxy(1,j,k1) = exy(1,j,k1)*at1
         fxy(2,j,k1) = exy(2,j,k1)*at1
         fxy(3,j,k1) = exy(3,j,k1)*at1
   60    continue
      endif
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
      subroutine WFFT2RVMX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,  
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
         call FFT2RVMXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,    
     1nxhyd,nxyhd)
c perform y fft
         call FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c perform x fft
         call FFT2RVMXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,    
     1nxhyd,nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT2RVM3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,  
     1nxyhd)
c wrapper function for 3 2d real to complex ffts
c parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(4,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2RVM3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,    
     1nxhyd,nxyhd)
c perform y fft
         call FFT2RVM3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,   
     1nxhyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RVM3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,   
     1nxhyd,nxyhd)
c perform x fft
         call FFT2RVM3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,    
     1nxhyd,nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RVMXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd,
     1nxhyd,nxyhd)
c this subroutine performs the x part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic, with Vector/OpenMP
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
      integer nrx, i, j, k, l, j1, k1, k2, ns, ns2, km, kmr, nrxb
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
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,ani,t1,t2,t3)
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
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(j+k2,i)
      f(j+k2,i) = f(j+k1,i) - t2
      f(j+k1,i) = f(j+k1,i) + t2
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
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,t1,t2,t3)
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
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(j+k2,i)
      f(j+k2,i) = f(j+k1,i) - t2
      f(j+k1,i) = f(j+k1,i) + t2
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
      subroutine FFT2RVM3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd,
     1nxhyd,nxyhd)
c this subroutine performs the x part of 3 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c y, using complex arithmetic, with Vector/OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, three inverse fourier transforms are performed
c f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
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
      dimension f(4,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, jj, j1, k1, k2, ns, ns2, km, kmr
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
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,at1,at2,ani,t1,t2,t3,t4
!$OMP& )
      do 90 i = nyi, nyt
c swap complex components
      do 10 j = 1, nxh
      at1 = aimag(f(3,j,i))
      at2 = real(f(3,j,i))
      f(3,j,i) = cmplx(real(f(2,j,i)),real(f(4,j,i)))
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
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(1,j+k2,i)
      t3 = t1*f(2,j+k2,i)
      t4 = t1*f(3,j+k2,i)
      f(1,j+k2,i) = f(1,j+k1,i) - t2
      f(2,j+k2,i) = f(2,j+k1,i) - t3
      f(3,j+k2,i) = f(3,j+k1,i) - t4
      f(1,j+k1,i) = f(1,j+k1,i) + t2
      f(2,j+k1,i) = f(2,j+k1,i) + t3
      f(3,j+k1,i) = f(3,j+k1,i) + t4
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
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,at1,at2,t1,t2,t3,t4)
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
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(1,j+k2,i)
      t3 = t1*f(2,j+k2,i)
      t4 = t1*f(3,j+k2,i)
      f(1,j+k2,i) = f(1,j+k1,i) - t2
      f(2,j+k2,i) = f(2,j+k1,i) - t3
      f(3,j+k2,i) = f(3,j+k1,i) - t4
      f(1,j+k1,i) = f(1,j+k1,i) + t2
      f(2,j+k1,i) = f(2,j+k1,i) + t3
      f(3,j+k1,i) = f(3,j+k1,i) + t4
  150 continue
  160 continue
  170 continue
c swap complex components
      do 180 j = 1, nxh
      f(4,j,i) = cmplx(aimag(f(3,j,i)),aimag(f(4,j,i)))
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(aimag(f(1,j,i)),aimag(f(2,j,i)))
      at2 = real(f(2,j,i))
      f(2,j,i) = cmplx(at1,0.0)
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RVM3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd,
     1nxhyd,nxyhd)
c this subroutine performs the y part of 3 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c x, using complex arithmetic, with Vector/OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, three inverse fourier transforms are performed
c f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
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
      dimension f(4,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
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
