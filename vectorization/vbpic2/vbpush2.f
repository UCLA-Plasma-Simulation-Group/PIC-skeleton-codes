c Fortran Library for Skeleton 2-1/2D Electromagnetic Vector PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR2HT(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp,npe
     1,nx,ny,ipbc)
c for 2-1/2d code, this subroutine calculates initial particle
c co-ordinates and velocities with uniform density and maxwellian
c velocity with drift
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = velocity vx of particle n
c part(n,4) = velocity vy of particle n
c part(n,5) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 5
c npe = first dimension of particle array
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, idimp, npe, nx, ny, ipbc
      real vtx, vty, vtz, vdx, vdy, vdz
      real part
      dimension part(npe,idimp)
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
      part(j+k1,1) = edgelx + at1*(real(j) - 0.5)
      part(j+k1,2) = at3
   10 continue
   20 continue
c maxwellian velocity distribution
      do 30 j = 1, npxy
      part(j,3) = vtx*ranorm()
      part(j,4) = vty*ranorm()
      part(j,5) = vtz*ranorm()
   30 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 40 j = 1, npxy
      dsum1 = dsum1 + part(j,3)
      dsum2 = dsum2 + part(j,4)
      dsum3 = dsum3 + part(j,5)
   40 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1./real(npxy)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 50 j = 1, npxy
      part(j,3) = part(j,3) - sum1
      part(j,4) = part(j,4) - sum2
      part(j,5) = part(j,5) - sum3
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GBPUSH23LT(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,npe,nx,
     1ny,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells
c 119 flops/particle, 1 divide, 29 loads, 5 stores
c input: all, output: part, ek
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
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = velocity vx of particle n
c part(n,4) = velocity vy of particle n
c part(n,5) = velocity vz of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nop = number of particles
c npe = first dimension of particle array
c nx/ny = system length in x/y direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, npe, nx, ny, nxv, nyv, ipbc
      real qbm, dt, dtc, ek
      real part, fxy, bxy
      dimension part(npe,idimp)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
c local data
      integer j, nn, mm
      real qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz
      double precision sum1
      qtmh = 0.5*qbm*dt
      sum1 = 0.0d0
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
      do 10 j = 1, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn + nxv*mm + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*fxy(1,nn)
      dy = amx*fxy(2,nn)
      dz = amx*fxy(3,nn)
      dx = amy*(dxp*fxy(1,nn+1) + dx)
      dy = amy*(dxp*fxy(2,nn+1) + dy)
      dz = amy*(dxp*fxy(3,nn+1) + dz)
      acx = amx*fxy(1,nn+nxv)
      acy = amx*fxy(2,nn+nxv)
      acz = amx*fxy(3,nn+nxv)
      dx = dx + dyp*(dxp*fxy(1,nn+1+nxv) + acx) 
      dy = dy + dyp*(dxp*fxy(2,nn+1+nxv) + acy)
      dz = dz + dyp*(dxp*fxy(3,nn+1+nxv) + acz)
c find magnetic field
      ox = amx*bxy(1,nn)
      oy = amx*bxy(2,nn)
      oz = amx*bxy(3,nn)
      ox = amy*(dxp*bxy(1,nn+1) + ox)
      oy = amy*(dxp*bxy(2,nn+1) + oy)
      oz = amy*(dxp*bxy(3,nn+1) + oz)
      acx = amx*bxy(1,nn+nxv)
      acy = amx*bxy(2,nn+nxv)
      acz = amx*bxy(3,nn+nxv)
      ox = ox + dyp*(dxp*bxy(1,nn+1+nxv) + acx) 
      oy = oy + dyp*(dxp*bxy(2,nn+1+nxv) + acy)
      oz = oz + dyp*(dxp*bxy(3,nn+1+nxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(j,3) + dx
      acy = part(j,4) + dy
      acz = part(j,5) + dz
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
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
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
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
c set new velocity
      part(j,3) = vx
      part(j,4) = vy
      part(j,5) = vz
   10 continue
c normalize kinetic energy
      ek = ek + 0.5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine GRBPUSH23LT(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,npe
     1,nx,ny,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells
c 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: part, ek
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
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = momentum px of particle n
c part(n,4) = momentum py of particle n
c part(n,5) = momentum pz of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 5
c nop = number of particles
c npe = first dimension of particle array
c nx/ny = system length in x/y direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, npe, nx, ny, nxv, nyv, ipbc
      real qbm, dt, dtc, ci, ek
      real part, fxy, bxy
      dimension part(npe,idimp)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
c local data
      integer j, nn, mm
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz
      double precision sum1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum1 = 0.0d0
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
      do 10 j = 1, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn + nxv*mm + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*fxy(1,nn)
      dy = amx*fxy(2,nn)
      dz = amx*fxy(3,nn)
      dx = amy*(dxp*fxy(1,nn+1) + dx)
      dy = amy*(dxp*fxy(2,nn+1) + dy)
      dz = amy*(dxp*fxy(3,nn+1) + dz)
      acx = amx*fxy(1,nn+nxv)
      acy = amx*fxy(2,nn+nxv)
      acz = amx*fxy(3,nn+nxv)
      dx = dx + dyp*(dxp*fxy(1,nn+1+nxv) + acx) 
      dy = dy + dyp*(dxp*fxy(2,nn+1+nxv) + acy)
      dz = dz + dyp*(dxp*fxy(3,nn+1+nxv) + acz)
c find magnetic field
      ox = amx*bxy(1,nn)
      oy = amx*bxy(2,nn)
      oz = amx*bxy(3,nn)
      ox = amy*(dxp*bxy(1,nn+1) + ox)
      oy = amy*(dxp*bxy(2,nn+1) + oy)
      oz = amy*(dxp*bxy(3,nn+1) + oz)
      acx = amx*bxy(1,nn+nxv)
      acy = amx*bxy(2,nn+nxv)
      acz = amx*bxy(3,nn+nxv)
      ox = ox + dyp*(dxp*bxy(1,nn+1+nxv) + acx) 
      oy = oy + dyp*(dxp*bxy(2,nn+1+nxv) + acy)
      oz = oz + dyp*(dxp*bxy(3,nn+1+nxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(j,3) + dx
      acy = part(j,4) + dy
      acz = part(j,5) + dz
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
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
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
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
c set new momentum
      part(j,3) = vx
      part(j,4) = vy
      part(j,5) = vz
   10 continue
c normalize kinetic energy
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine VGBPUSH23LT(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,npe,nx
     1,ny,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c vectorizable version using guard cells
c 119 flops/particle, 1 divide, 29 loads, 5 stores
c input: all, output: part, ek
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
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = velocity vx of particle n
c part(n,4) = velocity vy of particle n
c part(n,5) = velocity vz of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nop = number of particles
c npe = first dimension of particle array
c nx/ny = system length in x/y direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, npe, nx, ny, nxv, nyv, ipbc
      real qbm, dt, dtc, ek
      real part, fxy, bxy
      dimension part(npe,idimp)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer i, j, k, ipp, joff, nps, nn, mm
      real qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,2)
      double precision sum1
      qtmh = 0.5*qbm*dt
      sum1 = 0.0d0
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
      ipp = nop/npblk
c outer loop over number of full blocks
      do 60 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn + nxv*mm
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
   10 continue
c find acceleration
      do 30 j = 1, npblk
      nn = n(j)
      mm = nn + nxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
      do 20 i = 1, lvect
      if (i.gt.2) nn = mm
      dx = dx + fxy(1,i+nn)*s1(j,i)
      dy = dy + fxy(2,i+nn)*s1(j,i)
      dz = dz + fxy(3,i+nn)*s1(j,i)
      ox = ox + bxy(1,i+nn)*s1(j,i)
      oy = oy + bxy(2,i+nn)*s1(j,i)
      oz = oz + bxy(3,i+nn)*s1(j,i)
   20 continue
      s1(j,1) = dx
      s1(j,2) = dy
      s1(j,3) = dz
      s2(j,1) = ox
      s2(j,2) = oy
      s2(j,3) = oz
   30 continue
c new velocity
      do 40 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
c calculate half impulse
      dx = qtmh*s1(j,1)
      dy = qtmh*s1(j,2)
      dz = qtmh*s1(j,3)
c half acceleration
      acx = part(j+joff,3) + dx
      acy = part(j+joff,4) + dy
      acz = part(j+joff,5) + dz
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
   40 continue
! check boundary conditions
!dir$ novector
      do 50 j = 1, npblk
      dx = s1(j,1)
      dy = s1(j,2)
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
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
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j+joff,1) = dx
      part(j+joff,2) = dy
c set new velocity
      part(j+joff,3) = vx
      part(j+joff,4) = vy
      part(j+joff,5) = vz
   50 continue
   60 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 70 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn + nxv*mm + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*fxy(1,nn)
      dy = amx*fxy(2,nn)
      dz = amx*fxy(3,nn)
      dx = amy*(dxp*fxy(1,nn+1) + dx)
      dy = amy*(dxp*fxy(2,nn+1) + dy)
      dz = amy*(dxp*fxy(3,nn+1) + dz)
      acx = amx*fxy(1,nn+nxv)
      acy = amx*fxy(2,nn+nxv)
      acz = amx*fxy(3,nn+nxv)
      dx = dx + dyp*(dxp*fxy(1,nn+1+nxv) + acx) 
      dy = dy + dyp*(dxp*fxy(2,nn+1+nxv) + acy)
      dz = dz + dyp*(dxp*fxy(3,nn+1+nxv) + acz)
c find magnetic field
      ox = amx*bxy(1,nn)
      oy = amx*bxy(2,nn)
      oz = amx*bxy(3,nn)
      ox = amy*(dxp*bxy(1,nn+1) + ox)
      oy = amy*(dxp*bxy(2,nn+1) + oy)
      oz = amy*(dxp*bxy(3,nn+1) + oz)
      acx = amx*bxy(1,nn+nxv)
      acy = amx*bxy(2,nn+nxv)
      acz = amx*bxy(3,nn+nxv)
      ox = ox + dyp*(dxp*bxy(1,nn+1+nxv) + acx) 
      oy = oy + dyp*(dxp*bxy(2,nn+1+nxv) + acy)
      oz = oz + dyp*(dxp*bxy(3,nn+1+nxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(j,3) + dx
      acy = part(j,4) + dy
      acz = part(j,5) + dz
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
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
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
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
c set new velocity
      part(j,3) = vx
      part(j,4) = vy
      part(j,5) = vz
   70 continue
c normalize kinetic energy
      ek = ek + 0.5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine VGRBPUSH23LT(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,  
     1npe,nx,ny,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c vectorizable version using guard cells
c 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: part, ek
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
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = momentum px of particle n
c part(n,4) = momentum py of particle n
c part(n,5) = momentum pz of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 5
c nop = number of particles
c npe = first dimension of particle array
c nx/ny = system length in x/y direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, npe, nx, ny, nxv, nyv, ipbc
      real qbm, dt, dtc, ci, ek
      real part, fxy, bxy
      dimension part(npe,idimp)
      dimension fxy(4,nxv*nyv), bxy(4,nxv*nyv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer i, j, k, ipp, joff, nps, nn, mm
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,2)
      double precision sum1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum1 = 0.0d0
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
      ipp = nop/npblk
c outer loop over number of full blocks
      do 60 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn + nxv*mm
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
   10 continue
c find acceleration
      do 30 j = 1, npblk
      nn = n(j)
      mm = nn + nxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
      do 20 i = 1, lvect
      if (i.gt.2) nn = mm
      dx = dx + fxy(1,i+nn)*s1(j,i)
      dy = dy + fxy(2,i+nn)*s1(j,i)
      dz = dz + fxy(3,i+nn)*s1(j,i)
      ox = ox + bxy(1,i+nn)*s1(j,i)
      oy = oy + bxy(2,i+nn)*s1(j,i)
      oz = oz + bxy(3,i+nn)*s1(j,i)
   20 continue
      s1(j,1) = dx
      s1(j,2) = dy
      s1(j,3) = dz
      s2(j,1) = ox
      s2(j,2) = oy
      s2(j,3) = oz
   30 continue
c new momentum
      do 40 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
c calculate half impulse
      dx = qtmh*s1(j,1)
      dy = qtmh*s1(j,2)
      dz = qtmh*s1(j,3)
c half acceleration
      acx = part(j+joff,3) + dx
      acy = part(j+joff,4) + dy
      acz = part(j+joff,5) + dz
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
   40 continue
! check boundary conditions
!dir$ novector
      do 50 j = 1, npblk
      dx = s1(j,1)
      dy = s1(j,2)
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
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
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j+joff,1) = dx
      part(j+joff,2) = dy
c set new momentum
      part(j+joff,3) = vx
      part(j+joff,4) = vy
      part(j+joff,5) = vz
   50 continue
   60 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 70 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn + nxv*mm + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*fxy(1,nn)
      dy = amx*fxy(2,nn)
      dz = amx*fxy(3,nn)
      dx = amy*(dxp*fxy(1,nn+1) + dx)
      dy = amy*(dxp*fxy(2,nn+1) + dy)
      dz = amy*(dxp*fxy(3,nn+1) + dz)
      acx = amx*fxy(1,nn+nxv)
      acy = amx*fxy(2,nn+nxv)
      acz = amx*fxy(3,nn+nxv)
      dx = dx + dyp*(dxp*fxy(1,nn+1+nxv) + acx) 
      dy = dy + dyp*(dxp*fxy(2,nn+1+nxv) + acy)
      dz = dz + dyp*(dxp*fxy(3,nn+1+nxv) + acz)
c find magnetic field
      ox = amx*bxy(1,nn)
      oy = amx*bxy(2,nn)
      oz = amx*bxy(3,nn)
      ox = amy*(dxp*bxy(1,nn+1) + ox)
      oy = amy*(dxp*bxy(2,nn+1) + oy)
      oz = amy*(dxp*bxy(3,nn+1) + oz)
      acx = amx*bxy(1,nn+nxv)
      acy = amx*bxy(2,nn+nxv)
      acz = amx*bxy(3,nn+nxv)
      ox = ox + dyp*(dxp*bxy(1,nn+1+nxv) + acx) 
      oy = oy + dyp*(dxp*bxy(2,nn+1+nxv) + acy)
      oz = oz + dyp*(dxp*bxy(3,nn+1+nxv) + acz)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(j,3) + dx
      acy = part(j,4) + dy
      acz = part(j,5) + dz
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
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
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
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
c set new momentum
      part(j,3) = vx
      part(j,4) = vy
      part(j,5) = vz
   70 continue
c normalize kinetic energy
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine GPOST2LT(part,q,qm,nop,npe,idimp,nxv,nyv)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c q(j,k) = charge density at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 4
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
      implicit none
      integer nop, npe, idimp, nxv, nyv
      real qm
      real part, q
      dimension part(npe,idimp), q(nxv,nyv)
c local data
      integer j, nn, mm
      real x, y, dxp, dyp, amx, amy
      do 10 j = 1, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit charge
      x = q(nn,mm) + amx*amy
      y = q(nn+1,mm) + dxp*amy
      q(nn,mm) = x
      q(nn+1,mm) = y
      x = q(nn,mm+1) + amx*dyp
      y = q(nn+1,mm+1) + dxp*dyp
      q(nn,mm+1) = x
      q(nn+1,mm+1) = y
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine VGPOST2LT(part,q,qm,nop,npe,idimp,nxv,nyv)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c vectorizable version using guard cells
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c q(j,k) = charge density at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 4
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
      implicit none
      integer nop, npe, idimp, nxv, nyv
      real qm
      real part, q
      dimension part(npe,idimp), q(nxv*nyv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer i, j, k, ipp, joff, nps, nn, mm
      real x, y, dxp, dyp, amx, amy
c scratch arrays
      integer n
      real s
      dimension n(npblk), s(npblk,lvect)
      ipp = nop/npblk
c outer loop over number of full blocks
      do 40 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      n(j) = nn + nxv*mm
      amx = qm - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
   10 continue
c deposit charge
      do 30 j = 1, npblk
      nn = n(j)
      mm = nn + nxv - 2
!dir$ ivdep
      do 20 i = 1, lvect
      if (i.gt.2) nn = mm
      q(i+nn) = q(i+nn) + s(j,i)
   20 continue
   30 continue
   40 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 50 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn + nxv*mm + 1
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit charge
      x = q(nn) + amx*amy
      y = q(nn+1) + dxp*amy
      q(nn) = x
      q(nn+1) = y
      x = q(nn+nxv) + amx*dyp
      y = q(nn+nxv+1) + dxp*dyp
      q(nn+nxv) = x
      q(nn+nxv+1) = y
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GJPOST2LT(part,cu,qm,dt,nop,npe,idimp,nx,ny,nxv,nyv,   
     1ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 41 flops/particle, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = x velocity of particle n
c part(n,4) = y velocity of particle n
c part(n,5) = z velocity of particle n
c cu(i,j,k) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, npe, idimp, nx, ny, nxv, nyv, ipbc
      real qm, dt
      real part, cu
      dimension part(npe,idimp), cu(4,nxv*nyv)
c local data
      integer j, nn, mm
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
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
      do 10 j = 1, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn + nxv*mm + 1
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = part(j,3)
      vy = part(j,4)
      vz = part(j,5)
      cu(1,nn) = cu(1,nn) + vx*dx
      cu(2,nn) = cu(2,nn) + vy*dx
      cu(3,nn) = cu(3,nn) + vz*dx
      dx = amx*dyp
      cu(1,nn+1) = cu(1,nn+1) + vx*dy
      cu(2,nn+1) = cu(2,nn+1) + vy*dy
      cu(3,nn+1) = cu(3,nn+1) + vz*dy
      dy = dxp*dyp
      cu(1,nn+nxv) = cu(1,nn+nxv) + vx*dx
      cu(2,nn+nxv) = cu(2,nn+nxv) + vy*dx
      cu(3,nn+nxv) = cu(3,nn+nxv) + vz*dx
      cu(1,nn+1+nxv) = cu(1,nn+1+nxv) + vx*dy
      cu(2,nn+1+nxv) = cu(2,nn+1+nxv) + vy*dy
      cu(3,nn+1+nxv) = cu(3,nn+1+nxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,3) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,4) = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,3) = -vx
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GRJPOST2LT(part,cu,qm,dt,ci,nop,npe,idimp,nx,ny,nxv,nyv
     1,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = x momentum of particle n
c part(n,4) = y momentum of particle n
c part(n,5) = z momentum of particle n
c cu(i,j,k) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, npe, idimp, nx, ny, nxv, nyv, ipbc
      real qm, dt, ci
      real part, cu
      dimension part(npe,idimp), cu(4,nxv*nyv)
c local data
      integer j, nn, mm
      real ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami
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
      do 10 j = 1, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
c find inverse gamma
      ux = part(j,3)
      uy = part(j,4)
      uz = part(j,5)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn + nxv*mm + 1
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ux*gami
      vy = uy*gami
      vz = uz*gami
      cu(1,nn) = cu(1,nn) + vx*dx
      cu(2,nn) = cu(2,nn) + vy*dx
      cu(3,nn) = cu(3,nn) + vz*dx
      dx = amx*dyp
      cu(1,nn+1) = cu(1,nn+1) + vx*dy
      cu(2,nn+1) = cu(2,nn+1) + vy*dy
      cu(3,nn+1) = cu(3,nn+1) + vz*dy
      dy = dxp*dyp
      cu(1,nn+nxv) = cu(1,nn+nxv) + vx*dx
      cu(2,nn+nxv) = cu(2,nn+nxv) + vy*dx
      cu(3,nn+nxv) = cu(3,nn+nxv) + vz*dx
      cu(1,nn+1+nxv) = cu(1,nn+1+nxv) + vx*dy
      cu(2,nn+1+nxv) = cu(2,nn+1+nxv) + vy*dy
      cu(3,nn+1+nxv) = cu(3,nn+1+nxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,3) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,4) = -uy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,3) = -ux
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine VGJPOST2LT(part,cu,qm,dt,nop,npe,idimp,nx,ny,nxv,nyv,  
     1ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c vectorizable version using guard cells
c 41 flops/particle, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = x velocity of particle n
c part(n,4) = y velocity of particle n
c part(n,5) = z velocity of particle n
c cu(i,j,k) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, npe, idimp, nx, ny, nxv, nyv, ipbc
      real qm, dt
      real part, cu
      dimension part(npe,idimp), cu(4,nxv*nyv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer i, j, k, ipp, joff, nps, nn, mm
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,2)
      ipp = nop/npblk
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
c outer loop over number of full blocks
      do 50 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      n(j) = nn + nxv*mm
      amx = qm - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
      s2(j,1) = part(j+joff,3)
      s2(j,2) = part(j+joff,4)
      s2(j,3) = part(j+joff,5)
   10 continue
c deposit current
      do 30 j = 1, npblk
      nn = n(j)
      mm = nn + nxv - 2
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
!dir$ ivdep
      do 20 i = 1, lvect
      if (i.gt.2) nn = mm
      cu(1,i+nn) = cu(1,i+nn) + vx*s1(j,i)
      cu(2,i+nn) = cu(2,i+nn) + vy*s1(j,i)
      cu(3,i+nn) = cu(3,i+nn) + vz*s1(j,i)
   20 continue
   30 continue
c advance position half a time-step
!dir$ novector
      do 40 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      vx = s2(j,1)
      vy = s2(j,2)
      dx = x + vx*dt
      dy = y + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j+joff,3) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j+joff,4) = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j+joff,3) = -vx
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j+joff,1) = dx
      part(j+joff,2) = dy
   40 continue
   50 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 60 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn + nxv*mm + 1
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = part(j,3)
      vy = part(j,4)
      vz = part(j,5)
      cu(1,nn) = cu(1,nn) + vx*dx
      cu(2,nn) = cu(2,nn) + vy*dx
      cu(3,nn) = cu(3,nn) + vz*dx
      dx = amx*dyp
      cu(1,nn+1) = cu(1,nn+1) + vx*dy
      cu(2,nn+1) = cu(2,nn+1) + vy*dy
      cu(3,nn+1) = cu(3,nn+1) + vz*dy
      dy = dxp*dyp
      cu(1,nn+nxv) = cu(1,nn+nxv) + vx*dx
      cu(2,nn+nxv) = cu(2,nn+nxv) + vy*dx
      cu(3,nn+nxv) = cu(3,nn+nxv) + vz*dx
      cu(1,nn+1+nxv) = cu(1,nn+1+nxv) + vx*dy
      cu(2,nn+1+nxv) = cu(2,nn+1+nxv) + vy*dy
      cu(3,nn+1+nxv) = cu(3,nn+1+nxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,3) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,4) = -vy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,3) = -vx
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine VGRJPOST2LT(part,cu,qm,dt,ci,nop,npe,idimp,nx,ny,nxv,  
     1nyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c vectorizable version using guard cells
c 47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = x momentum of particle n
c part(n,4) = y momentum of particle n
c part(n,5) = z momentum of particle n
c cu(i,j,k) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, npe, idimp, nx, ny, nxv, nyv, ipbc
      real qm, dt, ci
      real part, cu
      dimension part(npe,idimp), cu(4,nxv*nyv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer i, j, k, ipp, joff, nps, nn, mm
      real ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,4)
      ci2 = ci*ci
      ipp = nop/npblk
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
c outer loop over number of full blocks
      do 50 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      n(j) = nn + nxv*mm
      amx = qm - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
c find inverse gamma
      ux = part(j+joff,3)
      uy = part(j+joff,4)
      uz = part(j+joff,5)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      s2(j,1) = ux*gami
      s2(j,2) = uy*gami
      s2(j,3) = uz*gami
      t(j,3) = ux
      t(j,4) = uy
   10 continue
c deposit current
      do 30 j = 1, npblk
      nn = n(j)
      mm = nn + nxv - 2
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
!dir$ ivdep
      do 20 i = 1, lvect
      if (i.gt.2) nn = mm
      cu(1,i+nn) = cu(1,i+nn) + vx*s1(j,i)
      cu(2,i+nn) = cu(2,i+nn) + vy*s1(j,i)
      cu(3,i+nn) = cu(3,i+nn) + vz*s1(j,i)
   20 continue
   30 continue
c advance position half a time-step
!dir$ novector
      do 40 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      vx = s2(j,1)
      vy = s2(j,2)
      ux = t(j,3)
      uy = t(j,4)
      dx = x + vx*dt
      dy = y + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j+joff,3) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j+joff,4) = -uy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j+joff,3) = -ux
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j+joff,1) = dx
      part(j+joff,2) = dy
   40 continue
   50 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 60 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
c find inverse gamma
      ux = part(j,3)
      uy = part(j,4)
      uz = part(j,5)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn + nxv*mm + 1
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = ux*gami
      vy = uy*gami
      vz = uz*gami
      cu(1,nn) = cu(1,nn) + vx*dx
      cu(2,nn) = cu(2,nn) + vy*dx
      cu(3,nn) = cu(3,nn) + vz*dx
      dx = amx*dyp
      cu(1,nn+1) = cu(1,nn+1) + vx*dy
      cu(2,nn+1) = cu(2,nn+1) + vy*dy
      cu(3,nn+1) = cu(3,nn+1) + vz*dy
      dy = dxp*dyp
      cu(1,nn+nxv) = cu(1,nn+nxv) + vx*dx
      cu(2,nn+nxv) = cu(2,nn+nxv) + vy*dx
      cu(3,nn+nxv) = cu(3,nn+nxv) + vz*dx
      cu(1,nn+1+nxv) = cu(1,nn+1+nxv) + vx*dy
      cu(2,nn+1+nxv) = cu(2,nn+1+nxv) + vy*dy
      cu(3,nn+1+nxv) = cu(3,nn+1+nxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,3) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,4) = -uy
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,3) = -ux
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DSORTP2YLT(parta,partb,npic,idimp,nop,npe,ny1)
c this subroutine sorts particles by y grid
c linear interpolation
c parta/partb = input/output particle arrays
c parta(n,2) = position y of particle n
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c npe = first dimension of particle array
c ny1 = system length in y direction + 1
      implicit none
      integer npic, idimp, nop, npe, ny1
      real parta, partb
      dimension parta(npe,idimp), partb(npe,idimp), npic(ny1)
c local data
      integer i, j, k, m, isum, ist, ip
c clear counter array
      do 10 k = 1, ny1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      m = parta(j,2)
      m = m + 1
      npic(m) = npic(m) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, ny1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, nop
      m = parta(j,2)
      m = m + 1
      ip = npic(m) + 1
      do 40 i = 1, idimp
      partb(ip,i) = parta(j,i)
   40 continue
      npic(m) = ip
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine BGUARD2L(bxy,nx,ny,nxe,nye)
c replicate extended periodic vector field bxy
c linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
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
c nye = second dimension of field arrays, must be >= ny+1
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
c nye = second dimension of field arrays, must be >= ny+1
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
      subroutine VPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,  
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
c vectorizable version
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
      double precision wp
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
   30 wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
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
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
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
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 70 j = 2, nxh
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
   70 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      we = real(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine CUPERP2(cu,nx,ny,nxvh,nyv)
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
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      cu(2,1,k) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      cu(1,j,1) = zero
      cu(1,j,k1) = zero
      cu(2,j,k1) = zero
   40 continue
      cu(1,1,1) = zero
      cu(2,1,1) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine VIBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
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
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
c vectorizable version
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
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
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
     1    + cu(3,j,k)*conjg(cu(3,j,k)) + cu(1,j,k1)*conjg(cu(1,j,k1))   
     2    + cu(2,j,k1)*conjg(cu(2,j,k1)) + cu(3,j,k1)*conjg(cu(3,j,k1)))
      wp = wp + dble(at1)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
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
     1    + cu(3,1,k)*conjg(cu(3,1,k)))
      wp = wp + dble(at1)
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 40 j = 2, nxh
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
     1    + cu(3,j,1)*conjg(cu(3,j,1)))
      wp = wp + dble(at1)
   40 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wm = real(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine VMAXWEL2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,nxhd
     1,nyhd)
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
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
c vectorizable version
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
      double precision wp, ws
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
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
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
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
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
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 40 j = 2, nxh
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
   40 continue
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
      wf = real(nx*ny)*ws
      wm = real(nx*ny)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine VEMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
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
      subroutine WFFT2RVX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,   
     1nxyhd)
c wrapper function for real to complex fft, with packed data
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
         call FFT2RVXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
c perform y fft
         call FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c perform x fft
         call FFT2RVXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT2RV3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,   
     1nxyhd)
c wrapper function for 3 2d real to complex ffts
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
         call FFT2RV3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
c perform y fft
         call FFT2RV3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RV3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c perform x fft
         call FFT2RV3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RVXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd, 
     1nxhyd,nxyhd)
c this subroutine performs the x part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic.
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2
c nyd = second dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, j1, k1, k2, ns, ns2, km, kmr
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
      if (isign.gt.0) go to 100
c inverse fourier transform
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 k = nyi, nyt
      t1 = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t1
   10 continue
   20 continue
c first transform in x
      nrx = nxy/nxh
      do 60 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 i = nyi, nyt
      do 30 j = 1, ns
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(j+k2,i)
      f(j+k2,i) = f(j+k1,i) - t2
      f(j+k1,i) = f(j+k1,i) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1.0/real(2*nx*ny)
      do 80 k = nyi, nyt
      do 70 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      t2 = conjg(f(nxh2-j,k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = ani*(t1 + t2)
      f(nxh2-j,k) = ani*conjg(t1 - t2)
   70 continue
   80 continue
      ani = 2.0*ani
      do 90 k = nyi, nyt
      f(nxhh+1,k) = ani*conjg(f(nxhh+1,k))
      f(1,k) = ani*cmplx(real(f(1,k)) + aimag(f(1,k)),                  
     1                   real(f(1,k)) - aimag(f(1,k)))
   90 continue
      return
c forward fourier transform
c scramble coefficients
  100 kmr = nxy/nx
      do 120 k = nyi, nyt
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      t2 = conjg(f(nxh2-j,k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = t1 + t2
      f(nxh2-j,k) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 k = nyi, nyt
      f(nxhh+1,k) = 2.0*conjg(f(nxhh+1,k))
      f(1,k) = cmplx(real(f(1,k)) + aimag(f(1,k)),                      
     1               real(f(1,k)) - aimag(f(1,k)))
  130 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 150 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 150
      do 140 k = nyi, nyt
      t1 = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t1
  140 continue
  150 continue
c then transform in x
      nrx = nxy/nxh
      do 190 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 i = nyi, nyt
      do 160 j = 1, ns
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(j+k2,i)
      f(j+k2,i) = f(j+k1,i) - t2
      f(j+k1,i) = f(j+k1,i) + t2
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd,  
     1nxhyd,nxyhd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = first dimension of f >= nx/2
c nyd = second dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr
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
      if (isign.gt.0) go to 80
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t1
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
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
      do 30 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      do 70 k = 2, nyh
      if (nxi.eq.1) then
         t1 = f(1,ny2-k)
         f(1,ny2-k) = 0.5*cmplx(aimag(f(1,k) + t1),real(f(1,k) - t1))
         f(1,k) = 0.5*cmplx(real(f(1,k) + t1),aimag(f(1,k) - t1))
      endif
   70 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
   80 do 90 k = 2, nyh
      if (nxi.eq.1) then
         t1 = cmplx(aimag(f(1,ny2-k)),real(f(1,ny2-k)))
         f(1,ny2-k) = conjg(f(1,k) - t1)
         f(1,k) = f(1,k) + t1
      endif
   90 continue
c bit-reverse array elements in y
      nry = nxhy/ny
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 110
      do 100 j = nxi, nxt
      t1 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t1
  100 continue
  110 continue
c first transform in y
      nry = nxy/ny
      do 150 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 120 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
  120 continue
  130 continue
  140 continue
  150 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RV3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd, 
     1nxhyd,nxyhd)
c this subroutine performs the x part of 3 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c y, using complex arithmetic
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
      if (isign.gt.0) go to 140
c inverse fourier transform
c swap complex components
      do 20 i = nyi, nyt
      do 10 j = 1, nxh
      at1 = aimag(f(3,j,i))
      at2 = real(f(3,j,i))
      f(3,j,i) = cmplx(real(f(2,j,i)),real(f(4,j,i)))
      f(2,j,i) = cmplx(aimag(f(1,j,i)),at1)
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
   10 continue
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 k = nyi, nyt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      t3 = f(3,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(3,j1,k) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
   30 continue
   40 continue
c first transform in x
      nrx = nxy/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 i = nyi, nyt
      do 50 j = 1, ns
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
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1.0/real(2*nx*ny)
      do 110 k = nyi, nyt
      do 100 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 90 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,nxh2-j,k) = ani*conjg(t1 - t2)
   90 continue
  100 continue
  110 continue
      ani = 2.*ani
      do 130 k = nyi, nyt
      do 120 jj = 1, 3
      f(jj,nxhh+1,k) = ani*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = ani*cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),         
     1                      real(f(jj,1,k)) - aimag(f(jj,1,k)))
  120 continue
  130 continue
      return
c forward fourier transform
c scramble coefficients
  140 kmr = nxy/nx
      do 170 k = nyi, nyt
      do 160 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 150 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,nxh2-j,k) = conjg(t1 - t2)
  150 continue
  160 continue
  170 continue
      do 190 k = nyi, nyt
      do 180 jj = 1, 3
      f(jj,nxhh+1,k) = 2.0*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),             
     1                  real(f(jj,1,k)) - aimag(f(jj,1,k)))
  180 continue
  190 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 210 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 210
      do 200 k = nyi, nyt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      t3 = f(3,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(3,j1,k) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  200 continue
  210 continue
c then transform in x
      nrx = nxy/nxh
      do 250 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 i = nyi, nyt
      do 220 j = 1, ns
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
  220 continue
  230 continue
  240 continue
  250 continue
c swap complex components
      do 270 i = nyi, nyt
      do 260 j = 1, nxh
      f(4,j,i) = cmplx(aimag(f(3,j,i)),aimag(f(4,j,i)))
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(aimag(f(1,j,i)),aimag(f(2,j,i)))
      at2 = real(f(2,j,i))
      f(2,j,i) = cmplx(at1,0.0)
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RV3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd, 
     1nxhyd,nxyhd)
c this subroutine performs the y part of 3 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c x, using complex arithmetic
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
      if (isign.gt.0) go to 90
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
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
      do 30 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      do 80 k = 2, nyh
      if (nxi.eq.1) then
         do 70 jj = 1, 3
         t1 = f(jj,1,ny2-k)
         f(jj,1,ny2-k) = 0.5*cmplx(aimag(f(jj,1,k) + t1),               
     1                             real(f(jj,1,k) - t1))
         f(jj,1,k) = 0.5*cmplx(real(f(jj,1,k) + t1),                    
     1                         aimag(f(jj,1,k) - t1))
   70    continue
      endif
   80 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
   90 do 110 k = 2, nyh
      if (nxi.eq.1) then
         do 100 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
         f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
         f(jj,1,k) = f(jj,1,k) + t1
  100    continue
      endif
  110 continue
c bit-reverse array elements in y
      nry = nxhy/ny
      do 130 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 130
      do 120 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  120 continue
  130 continue
c first transform in y
      nry = nxy/ny
      do 170 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 160 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 150 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 140 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
  140 continue
  150 continue
  160 continue
  170 continue
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
