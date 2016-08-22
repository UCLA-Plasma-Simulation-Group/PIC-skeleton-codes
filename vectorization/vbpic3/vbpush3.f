c Fortran Library for Skeleton 3D Electromagnetic Vector PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR3T(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,idimp,
     1npe,nx,ny,nz,ipbc)
c for 3d code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c part(n,4) = velocity vx of particle n
c part(n,5) = velocity vy of particle n
c part(n,6) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c idimp = size of phase space = 6
c npe = first dimension of particle array
c nx/ny/nz = system length in x/y/z direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, npz, idimp, npe, nx, ny, nz, ipbc
      real vtx, vty, vtz, vdx, vdy, vdz
      real part
      dimension part(npe,idimp)
c local data
      integer j, k, l, k1, l1, npxy, npxyz
      real edgelx, edgely, edgelz, at1, at2, at3, at4, at5
      real sum1, sum2, sum3
      double precision dsum1, dsum2, dsum3
      double precision ranorm
      npxy = npx*npy
      npxyz = npxy*npz
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
      endif
c uniform density profile
      do 30 l = 1, npz
      l1 = npxy*(l - 1)
      at5 = edgelz + at3*(real(l) - .5)
      do 20 k = 1, npy
      k1 = npx*(k - 1) + l1
      at4 = edgely + at2*(real(k) - .5)
      do 10 j = 1, npx
      part(j+k1,1) = edgelx + at1*(real(j) - .5)
      part(j+k1,2) = at4
      part(j+k1,3) = at5
   10 continue
   20 continue
   30 continue
c maxwellian velocity distribution
      do 40 j = 1, npxyz
      part(j,4) = vtx*ranorm()
      part(j,5) = vty*ranorm()
      part(j,6) = vtz*ranorm()
   40 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 50 j = 1, npxyz
      dsum1 = dsum1 + part(j,4)
      dsum2 = dsum2 + part(j,5)
      dsum3 = dsum3 + part(j,6)
   50 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1.0/real(npxyz)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 60 j = 1, npxyz
      part(j,4) = part(j,4) - sum1
      part(j,5) = part(j,5) - sum2
      part(j,6) = part(j,6) - sum3
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GBPUSH3LT(part,fxyz,bxyz,qbm,dt,dtc,ek,idimp,nop,npe,nx
     1,ny,nz,nxv,nyv,nzv,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field.  Using the Boris Mover.
c scalar version using guard cells
c 190 flops/particle, 1 divide, 54 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c z(t+dt)=z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c part(n,4) = velocity vx of particle n
c part(n,5) = velocity vy of particle n
c part(n,6) = velocity vz of particle n
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c bxyz(1,j,k,l) = x component of magnetic field at grid (j,k,l)
c bxyz(2,j,k,l) = y component of magnetic field at grid (j,k,l)
c bxyz(3,j,k,l) = z component of magnetic field at grid (j,k,l)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive force calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nop = number of particles
c npe = first dimension of particle array
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nop, npe, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qbm, dt, dtc, ek
      real part, fxyz, bxyz
      dimension part(npe,idimp)
      dimension fxyz(4,nxv*nyv*nzv), bxyz(4,nxv*nyv*nzv)
c local data
      integer j, nn, mm, ll, nxyv
      real qtmh, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz
      double precision sum1
      nxyv = nxv*nyv
      qtmh = 0.5*qbm*dt
      sum1 = 0.0d0
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
      do 10 j = 1, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find electric field
      dx = amx*fxyz(1,nn) + amy*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + amy*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + amy*fxyz(3,nn+1)
      dx = amz*(dx + dyp*fxyz(1,nn+nxv) + dx1*fxyz(1,nn+1+nxv))
      dy = amz*(dy + dyp*fxyz(2,nn+nxv) + dx1*fxyz(2,nn+1+nxv))
      dz = amz*(dz + dyp*fxyz(3,nn+nxv) + dx1*fxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*fxyz(1,mm) + amy*fxyz(1,mm+1)
      acy = amx*fxyz(2,mm) + amy*fxyz(2,mm+1)
      acz = amx*fxyz(3,mm) + amy*fxyz(3,mm+1)
      dx = dx + dzp*(acx + dyp*fxyz(1,mm+nxv) + dx1*fxyz(1,mm+1+nxv))
      dy = dy + dzp*(acy + dyp*fxyz(2,mm+nxv) + dx1*fxyz(2,mm+1+nxv))
      dz = dz + dzp*(acz + dyp*fxyz(3,mm+nxv) + dx1*fxyz(3,mm+1+nxv))
c find magnetic field
      ox = amx*bxyz(1,nn) + amy*bxyz(1,nn+1)
      oy = amx*bxyz(2,nn) + amy*bxyz(2,nn+1)
      oz = amx*bxyz(3,nn) + amy*bxyz(3,nn+1)
      ox = amz*(ox + dyp*bxyz(1,nn+nxv) + dx1*bxyz(1,nn+1+nxv))
      oy = amz*(oy + dyp*bxyz(2,nn+nxv) + dx1*bxyz(2,nn+1+nxv))
      oz = amz*(oz + dyp*bxyz(3,nn+nxv) + dx1*bxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*bxyz(1,mm) + amy*bxyz(1,mm+1)
      acy = amx*bxyz(2,mm) + amy*bxyz(2,mm+1)
      acz = amx*bxyz(3,mm) + amy*bxyz(3,mm+1)
      ox = ox + dzp*(acx + dyp*bxyz(1,mm+nxv) + dx1*bxyz(1,mm+1+nxv))
      oy = oy + dzp*(acy + dyp*bxyz(2,mm+nxv) + dx1*bxyz(2,mm+1+nxv))
      oz = oz + dzp*(acz + dyp*bxyz(3,mm+nxv) + dx1*bxyz(3,mm+1+nxv))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(j,4) + dx
      acy = part(j,5) + dy
      acz = part(j,6) + dz
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
      dz = z + vz*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
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
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
      part(j,3) = dz
c set new velocity
      part(j,4) = vx
      part(j,5) = vy
      part(j,6) = vz
   10 continue
c normalize kinetic energy
      ek = ek + 0.5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine GRBPUSH3LT(part,fxyz,bxyz,qbm,dt,dtc,ci,ek,idimp,nop,  
     1npe,nx,ny,nz,nxv,nyv,nzv,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells
c 202 flops/particle, 4 divides, 2 sqrts, 54 loads, 6 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
c omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
c omy = (q/m)*by(x(t),y(t),z(t))*gami,
c omz = (q/m)*bz(x(t),y(t),z(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c part(n,4) = momentum px of particle n
c part(n,5) = momentum py of particle n
c part(n,6) = momentum pz of particle n
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c bxyz(1,j,k,l) = x component of magnetic field at grid (j,k,l)
c bxyz(2,j,k,l) = y component of magnetic field at grid (j,k,l)
c bxyz(3,j,k,l) = z component of magnetic field at grid (j,k,l)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive force calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 6
c nop = number of particles
c npe = first dimension of particle array
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nop, npe, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qbm, dt, dtc, ci, ek
      real part, fxyz, bxyz
      dimension part(npe,idimp)
      dimension fxyz(4,nxv*nyv*nzv), bxyz(4,nxv*nyv*nzv)
c local data
      integer j, nn, mm, ll, nxyv
      real qtmh, ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, p2, gami, qtmg, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9, dtg
      real x, y, z, vx, vy, vz
      double precision sum1
      nxyv = nxv*nyv
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum1 = 0.0d0
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
      do 10 j = 1, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find electric field
      dx = amx*fxyz(1,nn) + amy*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + amy*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + amy*fxyz(3,nn+1)
      dx = amz*(dx + dyp*fxyz(1,nn+nxv) + dx1*fxyz(1,nn+1+nxv))
      dy = amz*(dy + dyp*fxyz(2,nn+nxv) + dx1*fxyz(2,nn+1+nxv))
      dz = amz*(dz + dyp*fxyz(3,nn+nxv) + dx1*fxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*fxyz(1,mm) + amy*fxyz(1,mm+1)
      acy = amx*fxyz(2,mm) + amy*fxyz(2,mm+1)
      acz = amx*fxyz(3,mm) + amy*fxyz(3,mm+1)
      dx = dx + dzp*(acx + dyp*fxyz(1,mm+nxv) + dx1*fxyz(1,mm+1+nxv))
      dy = dy + dzp*(acy + dyp*fxyz(2,mm+nxv) + dx1*fxyz(2,mm+1+nxv))
      dz = dz + dzp*(acz + dyp*fxyz(3,mm+nxv) + dx1*fxyz(3,mm+1+nxv))
c find magnetic field
      ox = amx*bxyz(1,nn) + amy*bxyz(1,nn+1)
      oy = amx*bxyz(2,nn) + amy*bxyz(2,nn+1)
      oz = amx*bxyz(3,nn) + amy*bxyz(3,nn+1)
      ox = amz*(ox + dyp*bxyz(1,nn+nxv) + dx1*bxyz(1,nn+1+nxv))
      oy = amz*(oy + dyp*bxyz(2,nn+nxv) + dx1*bxyz(2,nn+1+nxv))
      oz = amz*(oz + dyp*bxyz(3,nn+nxv) + dx1*bxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*bxyz(1,mm) + amy*bxyz(1,mm+1)
      acy = amx*bxyz(2,mm) + amy*bxyz(2,mm+1)
      acz = amx*bxyz(3,mm) + amy*bxyz(3,mm+1)
      ox = ox + dzp*(acx + dyp*bxyz(1,mm+nxv) + dx1*bxyz(1,mm+1+nxv))
      oy = oy + dzp*(acy + dyp*bxyz(2,mm+nxv) + dx1*bxyz(2,mm+1+nxv))
      oz = oz + dzp*(acz + dyp*bxyz(3,mm+nxv) + dx1*bxyz(3,mm+1+nxv))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(j,4) + dx
      acy = part(j,5) + dy
      acz = part(j,6) + dz
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
      dz = z + vz*dtg
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
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
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
      part(j,3) = dz
c set new momentum
      part(j,4) = vx
      part(j,5) = vy
      part(j,6) = vz
   10 continue
c normalize kinetic energy
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine VGBPUSH3LT(part,fxyz,bxyz,qbm,dt,dtc,ek,idimp,nop,npe, 
     1nx,ny,nz,nxv,nyv,nzv,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field.  Using the Boris Mover.
c vectorizable version using guard cells
c 190 flops/particle, 1 divide, 54 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c z(t+dt)=z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c part(n,4) = velocity vx of particle n
c part(n,5) = velocity vy of particle n
c part(n,6) = velocity vz of particle n
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c bxyz(1,j,k,l) = x component of magnetic field at grid (j,k,l)
c bxyz(2,j,k,l) = y component of magnetic field at grid (j,k,l)
c bxyz(3,j,k,l) = z component of magnetic field at grid (j,k,l)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive force calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nop = number of particles
c npe = first dimension of particle array
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nop, npe, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qbm, dt, dtc, ek
      real part, fxyz, bxyz
      dimension part(npe,idimp)
      dimension fxyz(4,nxv*nyv*nzv), bxyz(4,nxv*nyv*nzv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer i, j, k, l, ipp, joff, nps, nn, mm, ll, nxyv
      real qtmh, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,3), t(npblk,3)
      double precision sum1
      nxyv = nxv*nyv
      qtmh = 0.5*qbm*dt
      sum1 = 0.0d0
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
c loop over particles
      ipp = nop/npblk
c outer loop over number of full blocks
      do 60 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      z = part(j+joff,3)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn + nxv*mm + nxyv*ll
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s1(j,1) = amx*amz
      s1(j,2) = amy*amz
      s1(j,3) = dyp*amz
      s1(j,4) = dx1*amz
      s1(j,5) = amx*dzp
      s1(j,6) = amy*dzp
      s1(j,7) = dyp*dzp
      s1(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   10 continue
c find acceleration
      do 30 j = 1, npblk
      nn = n(j)
      mm = nn + nxv - 2
      ll = nn + nxyv - 4
      l = ll + nxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
      do 20 i = 1, lvect
      if (i.gt.6) then
         nn = l
      else if (i.gt.4) then
         nn = ll
      else if (i.gt.2) then
         nn = mm
      endif
      dx = dx + fxyz(1,i+nn)*s1(j,i)
      dy = dy + fxyz(2,i+nn)*s1(j,i)
      dz = dz + fxyz(3,i+nn)*s1(j,i)
      ox = ox + bxyz(1,i+nn)*s1(j,i)
      oy = oy + bxyz(2,i+nn)*s1(j,i)
      oz = oz + bxyz(3,i+nn)*s1(j,i)
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
      z = t(j,3)
c calculate half impulse
      dx = qtmh*s1(j,1)
      dy = qtmh*s1(j,2)
      dz = qtmh*s1(j,3)
c half acceleration
      acx = part(j+joff,4) + dx
      acy = part(j+joff,5) + dy
      acz = part(j+joff,6) + dz
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
      s1(j,3) = z + vz*dtc
      s2(j,1) = vx
      s2(j,2) = vy
      s2(j,3) = vz
   40 continue
c check boundary conditions
      do 50 j = 1, npblk
      dx = s1(j,1)
      dy = s1(j,2)
      dz = s1(j,3)
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
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
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = t(j,3)
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j+joff,1) = dx
      part(j+joff,2) = dy
      part(j+joff,3) = dz
c set new velocity
      part(j+joff,4) = vx
      part(j+joff,5) = vy
      part(j+joff,6) = vz
   50 continue
   60 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 70 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find electric field
      dx = amx*fxyz(1,nn) + amy*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + amy*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + amy*fxyz(3,nn+1)
      dx = amz*(dx + dyp*fxyz(1,nn+nxv) + dx1*fxyz(1,nn+1+nxv))
      dy = amz*(dy + dyp*fxyz(2,nn+nxv) + dx1*fxyz(2,nn+1+nxv))
      dz = amz*(dz + dyp*fxyz(3,nn+nxv) + dx1*fxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*fxyz(1,mm) + amy*fxyz(1,mm+1)
      acy = amx*fxyz(2,mm) + amy*fxyz(2,mm+1)
      acz = amx*fxyz(3,mm) + amy*fxyz(3,mm+1)
      dx = dx + dzp*(acx + dyp*fxyz(1,mm+nxv) + dx1*fxyz(1,mm+1+nxv))
      dy = dy + dzp*(acy + dyp*fxyz(2,mm+nxv) + dx1*fxyz(2,mm+1+nxv))
      dz = dz + dzp*(acz + dyp*fxyz(3,mm+nxv) + dx1*fxyz(3,mm+1+nxv))
c find magnetic field
      ox = amx*bxyz(1,nn) + amy*bxyz(1,nn+1)
      oy = amx*bxyz(2,nn) + amy*bxyz(2,nn+1)
      oz = amx*bxyz(3,nn) + amy*bxyz(3,nn+1)
      ox = amz*(ox + dyp*bxyz(1,nn+nxv) + dx1*bxyz(1,nn+1+nxv))
      oy = amz*(oy + dyp*bxyz(2,nn+nxv) + dx1*bxyz(2,nn+1+nxv))
      oz = amz*(oz + dyp*bxyz(3,nn+nxv) + dx1*bxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*bxyz(1,mm) + amy*bxyz(1,mm+1)
      acy = amx*bxyz(2,mm) + amy*bxyz(2,mm+1)
      acz = amx*bxyz(3,mm) + amy*bxyz(3,mm+1)
      ox = ox + dzp*(acx + dyp*bxyz(1,mm+nxv) + dx1*bxyz(1,mm+1+nxv))
      oy = oy + dzp*(acy + dyp*bxyz(2,mm+nxv) + dx1*bxyz(2,mm+1+nxv))
      oz = oz + dzp*(acz + dyp*bxyz(3,mm+nxv) + dx1*bxyz(3,mm+1+nxv))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(j,4) + dx
      acy = part(j,5) + dy
      acz = part(j,6) + dz
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
      dz = z + vz*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
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
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
      part(j,3) = dz
c set new velocity
      part(j,4) = vx
      part(j,5) = vy
      part(j,6) = vz
   70 continue
c normalize kinetic energy
      ek = ek + 0.5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine VGRBPUSH3LT(part,fxyz,bxyz,qbm,dt,dtc,ci,ek,idimp,nop, 
     1npe,nx,ny,nz,nxv,nyv,nzv,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c vectorizable version using guard cells
c 202 flops/particle, 4 divides, 2 sqrts, 54 loads, 6 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
c omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
c omy = (q/m)*by(x(t),y(t),z(t))*gami,
c omz = (q/m)*bz(x(t),y(t),z(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c part(n,4) = momentum px of particle n
c part(n,5) = momentum py of particle n
c part(n,6) = momentum pz of particle n
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c bxyz(1,j,k,l) = x component of magnetic field at grid (j,k,l)
c bxyz(2,j,k,l) = y component of magnetic field at grid (j,k,l)
c bxyz(3,j,k,l) = z component of magnetic field at grid (j,k,l)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive force calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 6
c nop = number of particles
c npe = first dimension of particle array
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nop, npe, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qbm, dt, dtc, ci, ek
      real part, fxyz, bxyz
      dimension part(npe,idimp)
      dimension fxyz(4,nxv*nyv*nzv), bxyz(4,nxv*nyv*nzv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer i, j, k, l, ipp, joff, nps, nn, mm, ll, nxyv
      real qtmh, ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, p2, gami, qtmg, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9, dtg
      real x, y, z, vx, vy, vz
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,3), t(npblk,3)
      double precision sum1
      nxyv = nxv*nyv
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum1 = 0.0d0
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
c loop over particles
      ipp = nop/npblk
c outer loop over number of full blocks
      do 60 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      z = part(j+joff,3)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn + nxv*mm + nxyv*ll
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s1(j,1) = amx*amz
      s1(j,2) = amy*amz
      s1(j,3) = dyp*amz
      s1(j,4) = dx1*amz
      s1(j,5) = amx*dzp
      s1(j,6) = amy*dzp
      s1(j,7) = dyp*dzp
      s1(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   10 continue
c find acceleration
      do 30 j = 1, npblk
      nn = n(j)
      mm = nn + nxv - 2
      ll = nn + nxyv - 4
      l = ll + nxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
      do 20 i = 1, lvect
      if (i.gt.6) then
         nn = l
      else if (i.gt.4) then
         nn = ll
      else if (i.gt.2) then
         nn = mm
      endif
      dx = dx + fxyz(1,i+nn)*s1(j,i)
      dy = dy + fxyz(2,i+nn)*s1(j,i)
      dz = dz + fxyz(3,i+nn)*s1(j,i)
      ox = ox + bxyz(1,i+nn)*s1(j,i)
      oy = oy + bxyz(2,i+nn)*s1(j,i)
      oz = oz + bxyz(3,i+nn)*s1(j,i)
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
      z = t(j,3)
c calculate half impulse
      dx = qtmh*s1(j,1)
      dy = qtmh*s1(j,2)
      dz = qtmh*s1(j,3)
c half acceleration
      acx = part(j+joff,4) + dx
      acy = part(j+joff,5) + dy
      acz = part(j+joff,6) + dz
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
      s1(j,3) = z + vz*dtg
      s2(j,1) = vx
      s2(j,2) = vy
      s2(j,3) = vz
   40 continue
c check boundary conditions
      do 50 j = 1, npblk
      dx = s1(j,1)
      dy = s1(j,2)
      dz = s1(j,3)
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
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
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = t(j,3)
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j+joff,1) = dx
      part(j+joff,2) = dy
      part(j+joff,3) = dz
c set new momentum
      part(j+joff,4) = vx
      part(j+joff,5) = vy
      part(j+joff,6) = vz
   50 continue
   60 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 70 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find electric field
      dx = amx*fxyz(1,nn) + amy*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + amy*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + amy*fxyz(3,nn+1)
      dx = amz*(dx + dyp*fxyz(1,nn+nxv) + dx1*fxyz(1,nn+1+nxv))
      dy = amz*(dy + dyp*fxyz(2,nn+nxv) + dx1*fxyz(2,nn+1+nxv))
      dz = amz*(dz + dyp*fxyz(3,nn+nxv) + dx1*fxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*fxyz(1,mm) + amy*fxyz(1,mm+1)
      acy = amx*fxyz(2,mm) + amy*fxyz(2,mm+1)
      acz = amx*fxyz(3,mm) + amy*fxyz(3,mm+1)
      dx = dx + dzp*(acx + dyp*fxyz(1,mm+nxv) + dx1*fxyz(1,mm+1+nxv))
      dy = dy + dzp*(acy + dyp*fxyz(2,mm+nxv) + dx1*fxyz(2,mm+1+nxv))
      dz = dz + dzp*(acz + dyp*fxyz(3,mm+nxv) + dx1*fxyz(3,mm+1+nxv))
c find magnetic field
      ox = amx*bxyz(1,nn) + amy*bxyz(1,nn+1)
      oy = amx*bxyz(2,nn) + amy*bxyz(2,nn+1)
      oz = amx*bxyz(3,nn) + amy*bxyz(3,nn+1)
      ox = amz*(ox + dyp*bxyz(1,nn+nxv) + dx1*bxyz(1,nn+1+nxv))
      oy = amz*(oy + dyp*bxyz(2,nn+nxv) + dx1*bxyz(2,nn+1+nxv))
      oz = amz*(oz + dyp*bxyz(3,nn+nxv) + dx1*bxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*bxyz(1,mm) + amy*bxyz(1,mm+1)
      acy = amx*bxyz(2,mm) + amy*bxyz(2,mm+1)
      acz = amx*bxyz(3,mm) + amy*bxyz(3,mm+1)
      ox = ox + dzp*(acx + dyp*bxyz(1,mm+nxv) + dx1*bxyz(1,mm+1+nxv))
      oy = oy + dzp*(acy + dyp*bxyz(2,mm+nxv) + dx1*bxyz(2,mm+1+nxv))
      oz = oz + dzp*(acz + dyp*bxyz(3,mm+nxv) + dx1*bxyz(3,mm+1+nxv))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(j,4) + dx
      acy = part(j,5) + dy
      acz = part(j,6) + dz
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
      dz = z + vz*dtg
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
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
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
      part(j,3) = dz
c set new momentum
      part(j,4) = vx
      part(j,5) = vy
      part(j,6) = vz
   70 continue
c normalize kinetic energy
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine V2GBPUSH3LT(part,fxyz,bxyz,qbm,dt,dtc,ek,idimp,nop,npe,
     1nx,ny,nz,nxv,nyv,nzv,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field.  Using the Boris Mover.
c vectorizable version using guard cells
c 190 flops/particle, 1 divide, 54 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c z(t+dt)=z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c part(n,4) = velocity vx of particle n
c part(n,5) = velocity vy of particle n
c part(n,6) = velocity vz of particle n
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c bxyz(1,j,k,l) = x component of magnetic field at grid (j,k,l)
c bxyz(2,j,k,l) = y component of magnetic field at grid (j,k,l)
c bxyz(3,j,k,l) = z component of magnetic field at grid (j,k,l)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive force calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nop = number of particles
c npe = first dimension of particle array
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nop, npe, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qbm, dt, dtc, ek
      real part, fxyz, bxyz
      dimension part(npe,idimp)
      dimension fxyz(4,nxv*nyv*nzv), bxyz(4,nxv*nyv*nzv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer i, j, k, ipp, joff, nps, nn, mm, ll, nxyv
      real qtmh, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, z, vx, vy, vz
c scratch arrays
      integer n, m
      real s1, s2, t
      dimension n(npblk), m(lvect)
      dimension s1(npblk,lvect), s2(npblk,3), t(npblk,3)
      double precision sum1
      nxyv = nxv*nyv
      m(1) = 0
      m(2) = 1
      m(3) = nxv
      m(4) = nxv + 1
      m(5) = nxyv
      m(6) = nxyv + 1
      m(7) = nxyv + nxv
      m(8) = nxyv + nxv + 1
      qtmh = 0.5*qbm*dt
      sum1 = 0.0d0
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
c loop over particles
      ipp = nop/npblk
c outer loop over number of full blocks
      do 60 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      z = part(j+joff,3)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn + nxv*mm + nxyv*ll + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s1(j,1) = amx*amz
      s1(j,2) = amy*amz
      s1(j,3) = dyp*amz
      s1(j,4) = dx1*amz
      s1(j,5) = amx*dzp
      s1(j,6) = amy*dzp
      s1(j,7) = dyp*dzp
      s1(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   10 continue
c find acceleration
      do 30 j = 1, npblk
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
!dir$ ivdep
      do 20 i = 1, lvect
      dx = dx + fxyz(1,n(j)+m(i))*s1(j,i)
      dy = dy + fxyz(2,n(j)+m(i))*s1(j,i)
      dz = dz + fxyz(3,n(j)+m(i))*s1(j,i)
      ox = ox + bxyz(1,n(j)+m(i))*s1(j,i)
      oy = oy + bxyz(2,n(j)+m(i))*s1(j,i)
      oz = oz + bxyz(3,n(j)+m(i))*s1(j,i)
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
      z = t(j,3)
c calculate half impulse
      dx = qtmh*s1(j,1)
      dy = qtmh*s1(j,2)
      dz = qtmh*s1(j,3)
c half acceleration
      acx = part(j+joff,4) + dx
      acy = part(j+joff,5) + dy
      acz = part(j+joff,6) + dz
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
      s1(j,3) = z + vz*dtc
      s2(j,1) = vx
      s2(j,2) = vy
      s2(j,3) = vz
   40 continue
c check boundary conditions
      do 50 j = 1, npblk
      dx = s1(j,1)
      dy = s1(j,2)
      dz = s1(j,3)
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
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
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = t(j,3)
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j+joff,1) = dx
      part(j+joff,2) = dy
      part(j+joff,3) = dz
c set new velocity
      part(j+joff,4) = vx
      part(j+joff,5) = vy
      part(j+joff,6) = vz
   50 continue
   60 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 70 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find electric field
      dx = amx*fxyz(1,nn) + amy*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + amy*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + amy*fxyz(3,nn+1)
      dx = amz*(dx + dyp*fxyz(1,nn+nxv) + dx1*fxyz(1,nn+1+nxv))
      dy = amz*(dy + dyp*fxyz(2,nn+nxv) + dx1*fxyz(2,nn+1+nxv))
      dz = amz*(dz + dyp*fxyz(3,nn+nxv) + dx1*fxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*fxyz(1,mm) + amy*fxyz(1,mm+1)
      acy = amx*fxyz(2,mm) + amy*fxyz(2,mm+1)
      acz = amx*fxyz(3,mm) + amy*fxyz(3,mm+1)
      dx = dx + dzp*(acx + dyp*fxyz(1,mm+nxv) + dx1*fxyz(1,mm+1+nxv))
      dy = dy + dzp*(acy + dyp*fxyz(2,mm+nxv) + dx1*fxyz(2,mm+1+nxv))
      dz = dz + dzp*(acz + dyp*fxyz(3,mm+nxv) + dx1*fxyz(3,mm+1+nxv))
c find magnetic field
      ox = amx*bxyz(1,nn) + amy*bxyz(1,nn+1)
      oy = amx*bxyz(2,nn) + amy*bxyz(2,nn+1)
      oz = amx*bxyz(3,nn) + amy*bxyz(3,nn+1)
      ox = amz*(ox + dyp*bxyz(1,nn+nxv) + dx1*bxyz(1,nn+1+nxv))
      oy = amz*(oy + dyp*bxyz(2,nn+nxv) + dx1*bxyz(2,nn+1+nxv))
      oz = amz*(oz + dyp*bxyz(3,nn+nxv) + dx1*bxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*bxyz(1,mm) + amy*bxyz(1,mm+1)
      acy = amx*bxyz(2,mm) + amy*bxyz(2,mm+1)
      acz = amx*bxyz(3,mm) + amy*bxyz(3,mm+1)
      ox = ox + dzp*(acx + dyp*bxyz(1,mm+nxv) + dx1*bxyz(1,mm+1+nxv))
      oy = oy + dzp*(acy + dyp*bxyz(2,mm+nxv) + dx1*bxyz(2,mm+1+nxv))
      oz = oz + dzp*(acz + dyp*bxyz(3,mm+nxv) + dx1*bxyz(3,mm+1+nxv))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(j,4) + dx
      acy = part(j,5) + dy
      acz = part(j,6) + dz
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
      dz = z + vz*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
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
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
      part(j,3) = dz
c set new velocity
      part(j,4) = vx
      part(j,5) = vy
      part(j,6) = vz
   70 continue
c normalize kinetic energy
      ek = ek + 0.5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine V2GRBPUSH3LT(part,fxyz,bxyz,qbm,dt,dtc,ci,ek,idimp,nop,
     1npe,nx,ny,nz,nxv,nyv,nzv,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c vectorizable version using guard cells
c 202 flops/particle, 4 divides, 2 sqrts, 54 loads, 6 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
c omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
c omy = (q/m)*by(x(t),y(t),z(t))*gami,
c omz = (q/m)*bz(x(t),y(t),z(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c part(n,4) = momentum px of particle n
c part(n,5) = momentum py of particle n
c part(n,6) = momentum pz of particle n
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c bxyz(1,j,k,l) = x component of magnetic field at grid (j,k,l)
c bxyz(2,j,k,l) = y component of magnetic field at grid (j,k,l)
c bxyz(3,j,k,l) = z component of magnetic field at grid (j,k,l)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive force calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 6
c nop = number of particles
c npe = first dimension of particle array
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nop, npe, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qbm, dt, dtc, ci, ek
      real part, fxyz, bxyz
      dimension part(npe,idimp)
      dimension fxyz(4,nxv*nyv*nzv), bxyz(4,nxv*nyv*nzv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer i, j, k, ipp, joff, nps, nn, mm, ll, nxyv
      real qtmh, ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, p2, gami, qtmg, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9, dtg
      real x, y, z, vx, vy, vz
c scratch arrays
      integer n, m
      real s1, s2, t
      dimension n(npblk), m(lvect)
      dimension s1(npblk,lvect), s2(npblk,3), t(npblk,3)
      double precision sum1
      nxyv = nxv*nyv
      m(1) = 0
      m(2) = 1
      m(3) = nxv
      m(4) = nxv + 1
      m(5) = nxyv
      m(6) = nxyv + 1
      m(7) = nxyv + nxv
      m(8) = nxyv + nxv + 1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum1 = 0.0d0
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
c loop over particles
      ipp = nop/npblk
c outer loop over number of full blocks
      do 60 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      z = part(j+joff,3)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn + nxv*mm + nxyv*ll + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s1(j,1) = amx*amz
      s1(j,2) = amy*amz
      s1(j,3) = dyp*amz
      s1(j,4) = dx1*amz
      s1(j,5) = amx*dzp
      s1(j,6) = amy*dzp
      s1(j,7) = dyp*dzp
      s1(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   10 continue
c find acceleration
      do 30 j = 1, npblk
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
!dir$ ivdep
      do 20 i = 1, lvect
      dx = dx + fxyz(1,n(j)+m(i))*s1(j,i)
      dy = dy + fxyz(2,n(j)+m(i))*s1(j,i)
      dz = dz + fxyz(3,n(j)+m(i))*s1(j,i)
      ox = ox + bxyz(1,n(j)+m(i))*s1(j,i)
      oy = oy + bxyz(2,n(j)+m(i))*s1(j,i)
      oz = oz + bxyz(3,n(j)+m(i))*s1(j,i)
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
      z = t(j,3)
c calculate half impulse
      dx = qtmh*s1(j,1)
      dy = qtmh*s1(j,2)
      dz = qtmh*s1(j,3)
c half acceleration
      acx = part(j+joff,4) + dx
      acy = part(j+joff,5) + dy
      acz = part(j+joff,6) + dz
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
      s1(j,3) = z + vz*dtg
      s2(j,1) = vx
      s2(j,2) = vy
      s2(j,3) = vz
   40 continue
c check boundary conditions
      do 50 j = 1, npblk
      dx = s1(j,1)
      dy = s1(j,2)
      dz = s1(j,3)
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
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
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = t(j,3)
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j+joff,1) = dx
      part(j+joff,2) = dy
      part(j+joff,3) = dz
c set new momentum
      part(j+joff,4) = vx
      part(j+joff,5) = vy
      part(j+joff,6) = vz
   50 continue
   60 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 70 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find electric field
      dx = amx*fxyz(1,nn) + amy*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + amy*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + amy*fxyz(3,nn+1)
      dx = amz*(dx + dyp*fxyz(1,nn+nxv) + dx1*fxyz(1,nn+1+nxv))
      dy = amz*(dy + dyp*fxyz(2,nn+nxv) + dx1*fxyz(2,nn+1+nxv))
      dz = amz*(dz + dyp*fxyz(3,nn+nxv) + dx1*fxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*fxyz(1,mm) + amy*fxyz(1,mm+1)
      acy = amx*fxyz(2,mm) + amy*fxyz(2,mm+1)
      acz = amx*fxyz(3,mm) + amy*fxyz(3,mm+1)
      dx = dx + dzp*(acx + dyp*fxyz(1,mm+nxv) + dx1*fxyz(1,mm+1+nxv))
      dy = dy + dzp*(acy + dyp*fxyz(2,mm+nxv) + dx1*fxyz(2,mm+1+nxv))
      dz = dz + dzp*(acz + dyp*fxyz(3,mm+nxv) + dx1*fxyz(3,mm+1+nxv))
c find magnetic field
      ox = amx*bxyz(1,nn) + amy*bxyz(1,nn+1)
      oy = amx*bxyz(2,nn) + amy*bxyz(2,nn+1)
      oz = amx*bxyz(3,nn) + amy*bxyz(3,nn+1)
      ox = amz*(ox + dyp*bxyz(1,nn+nxv) + dx1*bxyz(1,nn+1+nxv))
      oy = amz*(oy + dyp*bxyz(2,nn+nxv) + dx1*bxyz(2,nn+1+nxv))
      oz = amz*(oz + dyp*bxyz(3,nn+nxv) + dx1*bxyz(3,nn+1+nxv))
      mm = nn + nxyv
      acx = amx*bxyz(1,mm) + amy*bxyz(1,mm+1)
      acy = amx*bxyz(2,mm) + amy*bxyz(2,mm+1)
      acz = amx*bxyz(3,mm) + amy*bxyz(3,mm+1)
      ox = ox + dzp*(acx + dyp*bxyz(1,mm+nxv) + dx1*bxyz(1,mm+1+nxv))
      oy = oy + dzp*(acy + dyp*bxyz(2,mm+nxv) + dx1*bxyz(2,mm+1+nxv))
      oz = oz + dzp*(acz + dyp*bxyz(3,mm+nxv) + dx1*bxyz(3,mm+1+nxv))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(j,4) + dx
      acy = part(j,5) + dy
      acz = part(j,6) + dz
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
      dz = z + vz*dtg
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
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
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
      part(j,3) = dz
c set new momentum
      part(j,4) = vx
      part(j,5) = vy
      part(j,6) = vz
   70 continue
c normalize kinetic energy
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine GPOST3LT(part,q,qm,nop,npe,idimp,nxv,nyv,nzv)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells
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
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c q(j,k,l) = charge density at grid point j,k,l
c qm = charge on particle, in units of e
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 6
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
c nzv = third dimension of charge array, must be >= nz+1
      implicit none
      integer nop, npe, idimp, nxv, nyv, nzv
      real qm
      real part, q
      dimension part(npe,idimp), q(nxv*nyv*nzv)
c local data
      integer j, nn, mm, ll, nxyv
      real x, y, z, w, dx1, dxp, dyp, dzp, amx, amy, amz
      nxyv = nxv*nyv
      do 10 j = 1, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c deposit charge
      x = q(nn) + amx*amz
      y = q(nn+1) + amy*amz
      z = q(nn+nxv) + dyp*amz
      w = q(nn+1+nxv) + dx1*amz
      q(nn) = x
      q(nn+1) = y
      q(nn+nxv) = z
      q(nn+1+nxv) = w
      mm = nn + nxyv
      x = q(mm) + amx*dzp
      y = q(mm+1) + amy*dzp
      z = q(mm+nxv) + dyp*dzp
      w = q(mm+1+nxv) + dx1*dzp
      q(mm) = x
      q(mm+1) = y
      q(mm+nxv) = z
      q(mm+1+nxv) = w
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine VGPOST3LT(part,q,qm,nop,npe,idimp,nxv,nyv,nzv)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c vectorizable version using guard cells
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
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c q(j,k,l) = charge density at grid point j,k,l
c qm = charge on particle, in units of e
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 6
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
c nzv = third dimension of charge array, must be >= nz+1
      implicit none
      integer nop, npe, idimp, nxv, nyv, nzv
      real qm
      real part, q
      dimension part(npe,idimp), q(nxv*nyv*nzv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer i, j, k, ipp, joff, nps, nn, mm, ll, nxyv
      real x, y, z, w, dx1, dxp, dyp, dzp, amx, amy, amz
c scratch arrays
      integer n, m
      real s
      dimension n(npblk), m(lvect), s(npblk,lvect)
      nxyv = nxv*nyv
      m(1) = 0
      m(2) = 1
      m(3) = nxv
      m(4) = nxv + 1
      m(5) = nxyv
      m(6) = nxyv + 1
      m(7) = nxyv + nxv
      m(8) = nxyv + nxv + 1
c loop over particles
      ipp = nop/npblk
c outer loop over number of full blocks
      do 40 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      z = part(j+joff,3)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn + nxv*mm + nxyv*ll + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
   10 continue
c deposit charge
      do 30 j = 1, npblk
!dir$ ivdep
      do 20 i = 1, lvect
      q(n(j)+m(i)) = q(n(j)+m(i)) + s(j,i)
   20 continue
   30 continue
   40 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 50 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c deposit charge
      x = q(nn) + amx*amz
      y = q(nn+1) + amy*amz
      z = q(nn+nxv) + dyp*amz
      w = q(nn+1+nxv) + dx1*amz
      q(nn) = x
      q(nn+1) = y
      q(nn+nxv) = z
      q(nn+1+nxv) = w
      mm = nn + nxyv
      x = q(mm) + amx*dzp
      y = q(mm+1) + amy*dzp
      z = q(mm+nxv) + dyp*dzp
      w = q(mm+1+nxv) + dx1*dzp
      q(mm) = x
      q(mm+1) = y
      q(mm+nxv) = z
      q(mm+1+nxv) = w
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GJPOST3LT(part,cu,qm,dt,nop,npe,idimp,nx,ny,nz,nxv,nyv,
     1nzv,ipbc)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 69 flops/particle, 30 loads, 27 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c part(n,4) = x velocity of particle n
c part(n,5) = y velocity of particle n
c part(n,6) = z velocity of particle n
c cu(i,j,k,l) = ith component of current density at grid point j,k,l
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 6
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c nzv = fourth dimension of current array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer nop, npe, idimp, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qm, dt
      real part, cu
      dimension part(npe,idimp), cu(4,nxv*nyv*nzv)
c local data
      integer j, nn, mm, ll, nxyv
      real edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz
      real x, y, z
      nxyv = nxv*nyv
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
c find interpolation weights
      do 10 j = 1, nop
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c deposit current
      dx = amx*amz
      dy = amy*amz
      vx = part(j,4)
      vy = part(j,5)
      vz = part(j,6)
      cu(1,nn) = cu(1,nn) + vx*dx
      cu(2,nn) = cu(2,nn) + vy*dx
      cu(3,nn) = cu(3,nn) + vz*dx
      dx = dyp*amz
      cu(1,nn+1) = cu(1,nn+1) + vx*dy
      cu(2,nn+1) = cu(2,nn+1) + vy*dy
      cu(3,nn+1) = cu(3,nn+1) + vz*dy
      dy = dx1*amz
      cu(1,nn+nxv) = cu(1,nn+nxv) + vx*dx
      cu(2,nn+nxv) = cu(2,nn+nxv) + vy*dx
      cu(3,nn+nxv) = cu(3,nn+nxv) + vz*dx
      dx = amx*dzp
      cu(1,nn+1+nxv) = cu(1,nn+1+nxv) + vx*dy
      cu(2,nn+1+nxv) = cu(2,nn+1+nxv) + vy*dy
      cu(3,nn+1+nxv) = cu(3,nn+1+nxv) + vz*dy
      mm = nn + nxyv
      dy = amy*dzp
      cu(1,mm) = cu(1,mm) + vx*dx
      cu(2,mm) = cu(2,mm) + vy*dx
      cu(3,mm) = cu(3,mm) + vz*dx
      dx = dyp*dzp
      cu(1,mm+1) = cu(1,mm+1) + vx*dy
      cu(2,mm+1) = cu(2,mm+1) + vy*dy
      cu(3,mm+1) = cu(3,mm+1) + vz*dy
      dy = dx1*dzp
      cu(1,mm+nxv) = cu(1,mm+nxv) + vx*dx
      cu(2,mm+nxv) = cu(2,mm+nxv) + vy*dx
      cu(3,mm+nxv) = cu(3,mm+nxv) + vz*dx
      cu(1,mm+1+nxv) = cu(1,mm+1+nxv) + vx*dy
      cu(2,mm+1+nxv) = cu(2,mm+1+nxv) + vy*dy
      cu(3,mm+1+nxv) = cu(3,mm+1+nxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,4) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,5) = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            part(j,6) = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,4) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,5) = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
      part(j,3) = dz
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GRJPOST3LT(part,cu,qm,dt,ci,nop,npe,idimp,nx,ny,nz,nxv,
     1nyv,nzv,ipbc)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 79 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c part(n,4) = x momentum of particle n
c part(n,5) = y momentum of particle n
c part(n,6) = z momentum of particle n
c cu(i,j,k,l) = ith component of current density at grid point j,k,l
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 6
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c nzv = fourth dimension of current array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer nop, npe, idimp, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qm, dt, ci
      real part, cu
      dimension part(npe,idimp), cu(4,nxv*nyv*nzv)
c local data
      integer j, nn, mm, ll, nxyv
      real ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz
      real x, y, z, ux, uy, uz
      real p2, gami
      nxyv = nxv*nyv
      ci2 = ci*ci
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
c find interpolation weights
      do 10 j = 1, nop
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
c find inverse gamma
      ux = part(j,4)
      uy = part(j,5)
      uz = part(j,6)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c deposit current
      dx = amx*amz
      dy = amy*amz
      vx = ux*gami
      vy = uy*gami
      vz = uz*gami
      cu(1,nn) = cu(1,nn) + vx*dx
      cu(2,nn) = cu(2,nn) + vy*dx
      cu(3,nn) = cu(3,nn) + vz*dx
      dx = dyp*amz
      cu(1,nn+1) = cu(1,nn+1) + vx*dy
      cu(2,nn+1) = cu(2,nn+1) + vy*dy
      cu(3,nn+1) = cu(3,nn+1) + vz*dy
      dy = dx1*amz
      cu(1,nn+nxv) = cu(1,nn+nxv) + vx*dx
      cu(2,nn+nxv) = cu(2,nn+nxv) + vy*dx
      cu(3,nn+nxv) = cu(3,nn+nxv) + vz*dx
      dx = amx*dzp
      cu(1,nn+1+nxv) = cu(1,nn+1+nxv) + vx*dy
      cu(2,nn+1+nxv) = cu(2,nn+1+nxv) + vy*dy
      cu(3,nn+1+nxv) = cu(3,nn+1+nxv) + vz*dy
      mm = nn + nxyv
      dy = amy*dzp
      cu(1,mm) = cu(1,mm) + vx*dx
      cu(2,mm) = cu(2,mm) + vy*dx
      cu(3,mm) = cu(3,mm) + vz*dx
      dx = dyp*dzp
      cu(1,mm+1) = cu(1,mm+1) + vx*dy
      cu(2,mm+1) = cu(2,mm+1) + vy*dy
      cu(3,mm+1) = cu(3,mm+1) + vz*dy
      dy = dx1*dzp
      cu(1,mm+nxv) = cu(1,mm+nxv) + vx*dx
      cu(2,mm+nxv) = cu(2,mm+nxv) + vy*dx
      cu(3,mm+nxv) = cu(3,mm+nxv) + vz*dx
      cu(1,mm+1+nxv) = cu(1,mm+1+nxv) + vx*dy
      cu(2,mm+1+nxv) = cu(2,mm+1+nxv) + vy*dy
      cu(3,mm+1+nxv) = cu(3,mm+1+nxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,4) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,5) = -uy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            part(j,6) = -uz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,4) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,5) = -uy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
      part(j,3) = dz
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine VGJPOST3LT(part,cu,qm,dt,nop,npe,idimp,nx,ny,nz,nxv,nyv
     1,nzv,ipbc)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c vectorizable version using guard cells
c 69 flops/particle, 30 loads, 27 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c part(n,4) = x velocity of particle n
c part(n,5) = y velocity of particle n
c part(n,6) = z velocity of particle n
c cu(i,j,k,l) = ith component of current density at grid point j,k,l
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 6
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c nzv = fourth dimension of current array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer nop, npe, idimp, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qm, dt
      real part, cu
      dimension part(npe,idimp), cu(4,nxv*nyv*nzv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer i, j, k, ipp, joff, nps, nn, mm, ll, nxyv
      real edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz
      real x, y, z
c scratch arrays
      integer n, m
      real s1, s2, t
      dimension n(npblk), m(lvect)
      dimension s1(npblk,lvect), s2(npblk,3), t(npblk,3)
      nxyv = nxv*nyv
      m(1) = 0
      m(2) = 1
      m(3) = nxv
      m(4) = nxv + 1
      m(5) = nxyv
      m(6) = nxyv + 1
      m(7) = nxyv + nxv
      m(8) = nxyv + nxv + 1
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
c loop over particles
      ipp = nop/npblk
c outer loop over number of full blocks
      do 50 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      z = part(j+joff,3)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn + nxv*mm + nxyv*ll + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s1(j,1) = amx*amz
      s1(j,2) = amy*amz
      s1(j,3) = dyp*amz
      s1(j,4) = dx1*amz
      s1(j,5) = amx*dzp
      s1(j,6) = amy*dzp
      s1(j,7) = dyp*dzp
      s1(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
      s2(j,1) = part(j+joff,4)
      s2(j,2) = part(j+joff,5)
      s2(j,3) = part(j+joff,6)
   10 continue
c deposit charge
      do 30 j = 1, npblk
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
!dir$ ivdep
      do 20 i = 1, lvect
      cu(1,n(j)+m(i)) = cu(1,n(j)+m(i)) + vx*s1(j,i)
      cu(2,n(j)+m(i)) = cu(2,n(j)+m(i)) + vy*s1(j,i)
      cu(3,n(j)+m(i)) = cu(3,n(j)+m(i)) + vz*s1(j,i)
   20 continue
   30 continue
c advance position half a time-step
      do 40 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j+joff,4) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j+joff,5) = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            part(j+joff,6) = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j+joff,4) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j+joff,5) = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j+joff,1) = dx
      part(j+joff,2) = dy
      part(j+joff,3) = dz
   40 continue
   50 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 60 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c deposit current
      dx = amx*amz
      dy = amy*amz
      vx = part(j,4)
      vy = part(j,5)
      vz = part(j,6)
      cu(1,nn) = cu(1,nn) + vx*dx
      cu(2,nn) = cu(2,nn) + vy*dx
      cu(3,nn) = cu(3,nn) + vz*dx
      dx = dyp*amz
      cu(1,nn+1) = cu(1,nn+1) + vx*dy
      cu(2,nn+1) = cu(2,nn+1) + vy*dy
      cu(3,nn+1) = cu(3,nn+1) + vz*dy
      dy = dx1*amz
      cu(1,nn+nxv) = cu(1,nn+nxv) + vx*dx
      cu(2,nn+nxv) = cu(2,nn+nxv) + vy*dx
      cu(3,nn+nxv) = cu(3,nn+nxv) + vz*dx
      dx = amx*dzp
      cu(1,nn+1+nxv) = cu(1,nn+1+nxv) + vx*dy
      cu(2,nn+1+nxv) = cu(2,nn+1+nxv) + vy*dy
      cu(3,nn+1+nxv) = cu(3,nn+1+nxv) + vz*dy
      mm = nn + nxyv
      dy = amy*dzp
      cu(1,mm) = cu(1,mm) + vx*dx
      cu(2,mm) = cu(2,mm) + vy*dx
      cu(3,mm) = cu(3,mm) + vz*dx
      dx = dyp*dzp
      cu(1,mm+1) = cu(1,mm+1) + vx*dy
      cu(2,mm+1) = cu(2,mm+1) + vy*dy
      cu(3,mm+1) = cu(3,mm+1) + vz*dy
      dy = dx1*dzp
      cu(1,mm+nxv) = cu(1,mm+nxv) + vx*dx
      cu(2,mm+nxv) = cu(2,mm+nxv) + vy*dx
      cu(3,mm+nxv) = cu(3,mm+nxv) + vz*dx
      cu(1,mm+1+nxv) = cu(1,mm+1+nxv) + vx*dy
      cu(2,mm+1+nxv) = cu(2,mm+1+nxv) + vy*dy
      cu(3,mm+1+nxv) = cu(3,mm+1+nxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,4) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,5) = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            part(j,6) = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,4) = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,5) = -vy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
      part(j,3) = dz
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine VGRJPOST3LT(part,cu,qm,dt,ci,nop,npe,idimp,nx,ny,nz,nxv
     1,nyv,nzv,ipbc)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c vectorizable version using guard cells
c 79 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c part(n,4) = x momentum of particle n
c part(n,5) = y momentum of particle n
c part(n,6) = z momentum of particle n
c cu(i,j,k,l) = ith component of current density at grid point j,k,l
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 6
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c nzv = fourth dimension of current array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer nop, npe, idimp, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qm, dt, ci
      real part, cu
      dimension part(npe,idimp), cu(4,nxv*nyv*nzv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer i, j, k, ipp, joff, nps, nn, mm, ll, nxyv
      real ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz
      real x, y, z, ux, uy, uz
      real p2, gami
c scratch arrays
      integer n, m
      real s1, s2, t
      dimension n(npblk), m(lvect)
      dimension s1(npblk,lvect), s2(npblk,3), t(npblk,6)
      nxyv = nxv*nyv
      m(1) = 0
      m(2) = 1
      m(3) = nxv
      m(4) = nxv + 1
      m(5) = nxyv
      m(6) = nxyv + 1
      m(7) = nxyv + nxv
      m(8) = nxyv + nxv + 1
      ci2 = ci*ci
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
c loop over particles
      ipp = nop/npblk
c outer loop over number of full blocks
      do 50 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      z = part(j+joff,3)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn + nxv*mm + nxyv*ll + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s1(j,1) = amx*amz
      s1(j,2) = amy*amz
      s1(j,3) = dyp*amz
      s1(j,4) = dx1*amz
      s1(j,5) = amx*dzp
      s1(j,6) = amy*dzp
      s1(j,7) = dyp*dzp
      s1(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
c find inverse gamma
      ux = part(j+joff,4)
      uy = part(j+joff,5)
      uz = part(j+joff,6)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      s2(j,1) = ux*gami
      s2(j,2) = uy*gami
      s2(j,3) = uz*gami
      t(j,4) = ux
      t(j,5) = uy
      t(j,6) = uz
   10 continue
c deposit current
      do 30 j = 1, npblk
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
!dir$ ivdep
      do 20 i = 1, lvect
      cu(1,n(j)+m(i)) = cu(1,n(j)+m(i)) + vx*s1(j,i)
      cu(2,n(j)+m(i)) = cu(2,n(j)+m(i)) + vy*s1(j,i)
      cu(3,n(j)+m(i)) = cu(3,n(j)+m(i)) + vz*s1(j,i)
   20 continue
   30 continue
c advance position half a time-step
      do 40 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
      vx = s2(j,1)
      vy = s2(j,2)
      vz = s2(j,3)
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
      ux = t(j,4)
      uy = t(j,5)
      uz = t(j,6)
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j+joff,4) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j+joff,5) = -uy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            part(j+joff,6) = -uz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j+joff,4) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j+joff,5) = -uy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j+joff,1) = dx
      part(j+joff,2) = dy
      part(j+joff,3) = dz
   40 continue
   50 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 60 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
c find inverse gamma
      ux = part(j,4)
      uy = part(j,5)
      uz = part(j,6)
      p2 = ux*ux + uy*uy + uz*uz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c deposit current
      dx = amx*amz
      dy = amy*amz
      vx = ux*gami
      vy = uy*gami
      vz = uz*gami
      cu(1,nn) = cu(1,nn) + vx*dx
      cu(2,nn) = cu(2,nn) + vy*dx
      cu(3,nn) = cu(3,nn) + vz*dx
      dx = dyp*amz
      cu(1,nn+1) = cu(1,nn+1) + vx*dy
      cu(2,nn+1) = cu(2,nn+1) + vy*dy
      cu(3,nn+1) = cu(3,nn+1) + vz*dy
      dy = dx1*amz
      cu(1,nn+nxv) = cu(1,nn+nxv) + vx*dx
      cu(2,nn+nxv) = cu(2,nn+nxv) + vy*dx
      cu(3,nn+nxv) = cu(3,nn+nxv) + vz*dx
      dx = amx*dzp
      cu(1,nn+1+nxv) = cu(1,nn+1+nxv) + vx*dy
      cu(2,nn+1+nxv) = cu(2,nn+1+nxv) + vy*dy
      cu(3,nn+1+nxv) = cu(3,nn+1+nxv) + vz*dy
      mm = nn + nxyv
      dy = amy*dzp
      cu(1,mm) = cu(1,mm) + vx*dx
      cu(2,mm) = cu(2,mm) + vy*dx
      cu(3,mm) = cu(3,mm) + vz*dx
      dx = dyp*dzp
      cu(1,mm+1) = cu(1,mm+1) + vx*dy
      cu(2,mm+1) = cu(2,mm+1) + vy*dy
      cu(3,mm+1) = cu(3,mm+1) + vz*dy
      dy = dx1*dzp
      cu(1,mm+nxv) = cu(1,mm+nxv) + vx*dx
      cu(2,mm+nxv) = cu(2,mm+nxv) + vy*dx
      cu(3,mm+nxv) = cu(3,mm+nxv) + vz*dx
      cu(1,mm+1+nxv) = cu(1,mm+1+nxv) + vx*dy
      cu(2,mm+1+nxv) = cu(2,mm+1+nxv) + vy*dy
      cu(3,mm+1+nxv) = cu(3,mm+1+nxv) + vz*dy
c advance position half a time-step
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,4) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,5) = -uy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            part(j,6) = -uz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(j,4) = -ux
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(j,5) = -uy
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(j,1) = dx
      part(j,2) = dy
      part(j,3) = dz
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DSORTP3YZLT(parta,partb,npic,idimp,nop,npe,ny1,nyz1)
c this subroutine sorts particles by y,z grid
c linear interpolation
c part = particle array
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c npic = address offset for reordering particles
c idimp = size of phase space = 6
c nop = number of particles
c ny1 = system length in y direction + 1
c nyz1 = ny1*nz1, where nz1 = system length in z direction + 1
      implicit none
      integer npic, idimp, nop, npe, ny1, nyz1
      real parta, partb
      dimension parta(npe,idimp), partb(npe,idimp), npic(nyz1)
c local data
      integer i, j, k, m, l, isum, ist, ip
c clear counter array
      do 10 k = 1, nyz1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      m = parta(j,2)
      l = parta(j,3)
      l = m + ny1*l + 1
      npic(l) = npic(l) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyz1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, nop
      m = parta(j,2)
      l = parta(j,3)
      l = m + ny1*l + 1
      ip = npic(l) + 1
      do 40 i = 1, idimp
      partb(ip,i) = parta(j,i)
   40 continue
      npic(l) = ip
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine CGUARD3L(fxyz,nx,ny,nz,nxe,nye,nze)
c replicate extended periodic vector field fxyz
c linear interpolation
c nx/ny/nz = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
c nze = third dimension of field arrays, must be >= nz+1
      implicit none
      real fxyz
      integer nx, ny, nz, nxe, nye, nze
      dimension fxyz(4,nxe,nye,nze)
c local data
      integer j, k, l
      do 30 l = 1, nz
      do 10 k = 1, ny
      fxyz(1,nx+1,k,l) = fxyz(1,1,k,l)
      fxyz(2,nx+1,k,l) = fxyz(2,1,k,l)
      fxyz(3,nx+1,k,l) = fxyz(3,1,k,l)
   10 continue
      do 20 j = 1, nx
      fxyz(1,j,ny+1,l) = fxyz(1,j,1,l)
      fxyz(2,j,ny+1,l) = fxyz(2,j,1,l)
      fxyz(3,j,ny+1,l) = fxyz(3,j,1,l)
   20 continue
      fxyz(1,nx+1,ny+1,l) = fxyz(1,1,1,l)
      fxyz(2,nx+1,ny+1,l) = fxyz(2,1,1,l)
      fxyz(3,nx+1,ny+1,l) = fxyz(3,1,1,l)
   30 continue
      do 50 k = 1, ny
      do 40 j = 1, nx
      fxyz(1,j,k,nz+1) = fxyz(1,j,k,1)
      fxyz(2,j,k,nz+1) = fxyz(2,j,k,1)
      fxyz(3,j,k,nz+1) = fxyz(3,j,k,1)
   40 continue
      fxyz(1,nx+1,k,nz+1) = fxyz(1,1,k,1)
      fxyz(2,nx+1,k,nz+1) = fxyz(2,1,k,1)
      fxyz(3,nx+1,k,nz+1) = fxyz(3,1,k,1)
   50 continue
      do 60 j = 1, nx
      fxyz(1,j,ny+1,nz+1) = fxyz(1,j,1,1)
      fxyz(2,j,ny+1,nz+1) = fxyz(2,j,1,1)
      fxyz(3,j,ny+1,nz+1) = fxyz(3,j,1,1)
   60 continue
      fxyz(1,nx+1,ny+1,nz+1) = fxyz(1,1,1,1)
      fxyz(2,nx+1,ny+1,nz+1) = fxyz(2,1,1,1)
      fxyz(3,nx+1,ny+1,nz+1) = fxyz(3,1,1,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine ACGUARD3L(cu,nx,ny,nz,nxe,nye,nze)
c accumulate extended periodic field cu
c linear interpolation
c nx/ny/nz = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
c nze = third dimension of field arrays, must be >= nz+1
      implicit none
      integer nx, ny, nz, nxe, nye, nze
      real cu
      dimension cu(4,nxe,nye,nze)
c local data
      integer i, j, k, l
c accumulate edges of extended field
      do 60 l = 1, nz
      do 20 k = 1, ny
      do 10 i = 1, 3
      cu(i,1,k,l) = cu(i,1,k,l) + cu(i,nx+1,k,l)
      cu(i,nx+1,k,l) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, 3
      cu(i,j,1,l) = cu(i,j,1,l) + cu(i,j,ny+1,l)
      cu(i,j,ny+1,l) = 0.0
   30 continue
   40 continue
      do 50 i = 1, 3
      cu(i,1,1,l) = cu(i,1,1,l) + cu(i,nx+1,ny+1,l)
      cu(i,nx+1,ny+1,l) = 0.0
   50 continue
   60 continue
      do 100 k = 1, ny
      do 80 j = 1, nx
      do 70 i = 1, 3
      cu(i,j,k,1) = cu(i,j,k,1) + cu(i,j,k,nz+1)
      cu(i,j,k,nz+1) = 0.0
   70 continue
   80 continue
      do 90 i = 1, 3
      cu(i,1,k,1) = cu(i,1,k,1) + cu(i,nx+1,k,nz+1)
      cu(i,nx+1,k,nz+1) = 0.0
   90 continue
  100 continue
      do 120 j = 1, nx
      do 110 i = 1, 3
      cu(i,j,1,1) = cu(i,j,1,1) + cu(i,j,ny+1,nz+1)
      cu(i,j,ny+1,nz+1) = 0.0
  110 continue
  120 continue
      do 130 i = 1, 3
      cu(i,1,1,1) = cu(i,1,1,1) + cu(i,nx+1,ny+1,nz+1)
      cu(i,nx+1,ny+1,nz+1) = 0.0
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine AGUARD3L(q,nx,ny,nz,nxe,nye,nze)
c accumulate extended periodic scalar field q
c linear interpolation
c nx/ny/nz = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
c nze = third dimension of field arrays, must be >= nz+1
      implicit none
      real q
      integer nx, ny, nz, nxe, nye, nze
      dimension q(nxe,nye,nze)
      integer j, k, l
c accumulate edges of extended field
      do 30 l = 1, nz
      do 10 k = 1, ny
      q(1,k,l) = q(1,k,l) + q(nx+1,k,l)
      q(nx+1,k,l) = 0.0
   10 continue
      do 20 j = 1, nx
      q(j,1,l) = q(j,1,l) + q(j,ny+1,l)
      q(j,ny+1,l) = 0.0
   20 continue
      q(1,1,l) = q(1,1,l) + q(nx+1,ny+1,l)
      q(nx+1,ny+1,l) = 0.0
   30 continue
      do 50 k = 1, ny
      do 40 j = 1, nx
      q(j,k,1) = q(j,k,1) + q(j,k,nz+1)
      q(j,k,nz+1) = 0.0
   40 continue
      q(1,k,1) = q(1,k,1) + q(nx+1,k,nz+1)
      q(nx+1,k,nz+1) = 0.0
   50 continue
      do 60 j = 1, nx
      q(j,1,1) = q(j,1,1) + q(j,ny+1,nz+1)
      q(j,ny+1,nz+1) = 0.0
   60 continue
      q(1,1,1) = q(1,1,1) + q(nx+1,ny+1,nz+1)
      q(nx+1,ny+1,nz+1) = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine VPOIS33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxvh
     1,nyv,nzv,nxhd,nyhd,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.
c for isign = 0, output: ffc
c input: isign,ax,ay,az,affp,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
c for isign = -1, output: fxyz, we
c input: q,ffc,isign,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
c approximate flop count is:
c 59*nxc*nyc*nzc + 26*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c equation used is:
c fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*s(kx,ky,kz),
c fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*s(kx,ky,kz),
c fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
c fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
c fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
c fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
c q(j,k,l) = complex charge density for fourier mode (j-1,k-1,l-1)
c fxyz(1,j,k,l) = x component of complex force/charge
c fxyz(2,j,k,l) = y component of complex force/charge
c fxyz(3,j,k,l) = z component of complex force/charge
c all for fourier mode (j-1,k-1,l-1)
c aimag(ffc(j,k,l)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1,l-1)
c real(ffc(j,k,l)) = potential green's function g
c for fourier mode (j-1,k-1,l-1)
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nzv = third dimension of field arrays, must be >= nz
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
c nzhd = third dimension of form factor array, must be >= nzh
c vectorizable version
      implicit none
      integer isign, nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
      real ax, ay, az, affp, we
      complex q, fxyz, ffc
      dimension q(nxvh,nyv,nzv), fxyz(4,nxvh,nyv,nzv)
      dimension ffc(nxhd,nyhd,nzhd)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real dnx, dny, dnz, dkx, dky, dkz, at1, at2, at3, at4, at5, at6
      complex zero, zt1, zt2
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 40
c prepare form factor array
      do 30 l = 1, nzh
      dkz = dnz*real(l - 1)
      at1 = dkz*dkz
      at2 = (dkz*az)**2
      do 20 k = 1, nyh
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = (dky*ay)**2 + at2
      do 10 j = 1, nxh
      dkx = dnx*real(j - 1)
      at5 = dkx*dkx + at3
      at6 = exp(-.5*((dkx*ax)**2 + at4))
      if (at5.eq.0.) then
         ffc(j,k,l) = cmplx(affp,1.0)
      else
         ffc(j,k,l) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 90 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
      do 60 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
!dir$ ivdep
      do 50 j = 2, nxh
      at1 = real(ffc(j,k,l))*aimag(ffc(j,k,l))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      at4 = dkz*at1
      zt1 = cmplx(aimag(q(j,k,l)),-real(q(j,k,l)))
      zt2 = cmplx(aimag(q(j,k1,l)),-real(q(j,k1,l)))
      fxyz(1,j,k,l) = at2*zt1
      fxyz(2,j,k,l) = at3*zt1
      fxyz(3,j,k,l) = at4*zt1
      fxyz(1,j,k1,l) = at2*zt2
      fxyz(2,j,k1,l) = -at3*zt2
      fxyz(3,j,k1,l) = at4*zt2
      zt1 = cmplx(aimag(q(j,k,l1)),-real(q(j,k,l1)))
      zt2 = cmplx(aimag(q(j,k1,l1)),-real(q(j,k1,l1)))
      fxyz(1,j,k,l1) = at2*zt1
      fxyz(2,j,k,l1) = at3*zt1
      fxyz(3,j,k,l1) = -at4*zt1
      fxyz(1,j,k1,l1) = at2*zt2
      fxyz(2,j,k1,l1) = -at3*zt2
      fxyz(3,j,k1,l1) = -at4*zt2
      at1 = at1*(q(j,k,l)*conjg(q(j,k,l)) + q(j,k1,l)*conjg(q(j,k1,l))  
     1    + q(j,k,l1)*conjg(q(j,k,l1)) + q(j,k1,l1)*conjg(q(j,k1,l1)))
      wp = wp + dble(at1)
   50 continue
   60 continue
c mode numbers kx = 0, nx/2
!dir$ ivdep
      do 70 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffc(1,k,l))*aimag(ffc(1,k,l))
      at3 = dny*real(k - 1)*at1
      at4 = dkz*at1
      zt1 = cmplx(aimag(q(1,k,l)),-real(q(1,k,l)))
      zt2 = cmplx(aimag(q(1,k,l1)),-real(q(1,k,l1)))
      fxyz(1,1,k,l) = zero
      fxyz(2,1,k,l) = at3*zt1
      fxyz(3,1,k,l) = at4*zt1
      fxyz(1,1,k1,l) = zero
      fxyz(2,1,k1,l) = zero
      fxyz(3,1,k1,l) = zero
      fxyz(1,1,k,l1) = zero
      fxyz(2,1,k,l1) = at3*zt2
      fxyz(3,1,k,l1) = -at4*zt2
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      at1 = at1*(q(1,k,l)*conjg(q(1,k,l)) + q(1,k,l1)*conjg(q(1,k,l1)))
      wp = wp + dble(at1)
   70 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 80 j = 2, nxh
      at1 = real(ffc(j,1,l))*aimag(ffc(j,1,l))
      at2 = dnx*real(j - 1)*at1
      at4 = dkz*at1
      zt1 = cmplx(aimag(q(j,1,l)),-real(q(j,1,l)))
      zt2 = cmplx(aimag(q(j,1,l1)),-real(q(j,1,l1)))
      fxyz(1,j,1,l) = at2*zt1
      fxyz(2,j,1,l) = zero
      fxyz(3,j,1,l) = at4*zt1
      fxyz(1,j,k1,l) = zero
      fxyz(2,j,k1,l) = zero
      fxyz(3,j,k1,l) = zero
      fxyz(1,j,1,l1) = at2*zt2
      fxyz(2,j,1,l1) = zero
      fxyz(3,j,1,l1) = -at4*zt2
      fxyz(1,j,k1,l1) = zero
      fxyz(2,j,k1,l1) = zero
      fxyz(3,j,k1,l1) = zero
      at1 = at1*(q(j,1,l)*conjg(q(j,1,l)) + q(j,1,l1)*conjg(q(j,1,l1)))
      wp = wp + dble(at1)
   80 continue
c mode numbers kx = 0, nx/2
      at1 = real(ffc(1,1,l))*aimag(ffc(1,1,l))
      at4 = dkz*at1
      fxyz(1,1,1,l) = zero
      fxyz(2,1,1,l) = zero
      fxyz(3,1,1,l) = at4*cmplx(aimag(q(1,1,l)),-real(q(1,1,l)))
      fxyz(1,1,k1,l) = zero
      fxyz(2,1,k1,l) = zero
      fxyz(3,1,k1,l) = zero
      fxyz(1,1,1,l1) = zero
      fxyz(2,1,1,l1) = zero
      fxyz(3,1,1,l1) = zero
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      at1 = at1*(q(1,1,l)*conjg(q(1,1,l)))
      wp = wp + dble(at1)
   90 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 110 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
!dir$ ivdep
      do 100 j = 2, nxh
      at1 = real(ffc(j,k,1))*aimag(ffc(j,k,1))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k,1)),-real(q(j,k,1)))
      zt2 = cmplx(aimag(q(j,k1,1)),-real(q(j,k1,1)))
      fxyz(1,j,k,1) = at2*zt1
      fxyz(2,j,k,1) = at3*zt1
      fxyz(3,j,k,1) = zero
      fxyz(1,j,k1,1) = at2*zt2
      fxyz(2,j,k1,1) = -at3*zt2
      fxyz(3,j,k1,1) = zero
      fxyz(1,j,k,l1) = zero
      fxyz(2,j,k,l1) = zero
      fxyz(3,j,k,l1) = zero
      fxyz(1,j,k1,l1) = zero
      fxyz(2,j,k1,l1) = zero
      fxyz(3,j,k1,l1) = zero
      at1 = at1*(q(j,k,1)*conjg(q(j,k,1)) + q(j,k1,1)*conjg(q(j,k1,1)))
      wp = wp + dble(at1)
  100 continue
  110 continue
c mode numbers kx = 0, nx/2
!dir$ ivdep
      do 120 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffc(1,k,1))*aimag(ffc(1,k,1))
      at3 = dny*real(k - 1)*at1
      fxyz(1,1,k,1) = zero
      fxyz(2,1,k,1) = at3*cmplx(aimag(q(1,k,1)),-real(q(1,k,1)))
      fxyz(3,1,k,1) = zero
      fxyz(1,1,k1,1) = zero
      fxyz(2,1,k1,1) = zero
      fxyz(3,1,k1,1) = zero
      fxyz(1,1,k,l1) = zero
      fxyz(2,1,k,l1) = zero
      fxyz(3,1,k,l1) = zero
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      at1 = at1*(q(1,k,1)*conjg(q(1,k,1)))
      wp = wp + dble(at1)
  120 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 130 j = 2, nxh
      at1 = real(ffc(j,1,1))*aimag(ffc(j,1,1))
      at2 = dnx*real(j - 1)*at1
      fxyz(1,j,1,1) = at2*cmplx(aimag(q(j,1,1)),-real(q(j,1,1)))
      fxyz(2,j,1,1) = zero
      fxyz(3,j,1,1) = zero
      fxyz(1,j,k1,1) = zero
      fxyz(2,j,k1,1) = zero
      fxyz(3,j,k1,1) = zero
      fxyz(1,j,1,l1) = zero
      fxyz(2,j,1,l1) = zero
      fxyz(3,j,1,l1) = zero
      fxyz(1,j,k1,l1) = zero
      fxyz(2,j,k1,l1) = zero
      fxyz(3,j,k1,l1) = zero
      at1 = at1*(q(j,1,1)*conjg(q(j,1,1)))
      wp = wp + dble(at1)
  130 continue
      fxyz(1,1,1,1) = zero
      fxyz(2,1,1,1) = zero
      fxyz(3,1,1,1) = zero
      fxyz(1,1,k1,1) = zero
      fxyz(2,1,k1,1) = zero
      fxyz(3,1,k1,1) = zero
      fxyz(1,1,1,l1) = zero
      fxyz(2,1,1,l1) = zero
      fxyz(3,1,1,l1) = zero
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      we = real(nx)*real(ny)*real(nz)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine CUPERP3(cu,nx,ny,nz,nxvh,nyv,nzv)
c this subroutine calculates the transverse current in fourier space
c input: all, output: cu
c approximate flop count is:
c 100*nxc*nyc*nzc + 36*(nxc*nyc + nxc*nzc + nyc*nzc)
c and (nx/2)*nyc*nzc divides
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c the transverse current is calculated using the equation:
c cux(kx,ky,kz) = cux(kx,ky,kz) - kx*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c cuy(kx,ky,kz) = cuy(kx,ky,kz) - ky*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c cuz(kx,ky,kz) = cuz(kx,ky,kz) - kz*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c cux(kx=pi) = cuy(kx=pi) = cuz(kx=pi) = 0,
c cux(ky=pi) = cuy(ky=pi) = cux(ky=pi) = 0,
c cux(kz=pi) = cuy(kz=pi) = cuz(kz=pi) = 0,
c cux(kx=0,ky=0,kz=0) = cuy(kx=0,ky=0,kz=0) = cuz(kx=0,ky=0,kz=0) = 0.
c cu(i,j,k,l) = complex current density for fourier mode (j-1,k-1,l-1)
c nx/ny/nz = system length in x/y/z direction
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
c nzv = fourth dimension of field arrays, must be >= nz
      implicit none
      integer nx, ny, nz, nxvh, nyv, nzv
      complex cu
      dimension cu(4,nxvh,nyv,nzv)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real dnx, dny, dnz, dkx, dky, dkz, dky2, dkz2, dkyz2, at1
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
c calculate transverse part of current
      do 50 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
      dkz2 = dkz*dkz
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dkyz2 = dky*dky + dkz2
!dir$ ivdep
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx + dkyz2)
      zt1 = at1*(dkx*cu(1,j,k,l) + dky*cu(2,j,k,l) + dkz*cu(3,j,k,l))
      cu(1,j,k,l) = cu(1,j,k,l) - dkx*zt1
      cu(2,j,k,l) = cu(2,j,k,l) - dky*zt1
      cu(3,j,k,l) = cu(3,j,k,l) - dkz*zt1
      zt1 = at1*(dkx*cu(1,j,k1,l) - dky*cu(2,j,k1,l) + dkz*cu(3,j,k1,l))
      cu(1,j,k1,l) = cu(1,j,k1,l) - dkx*zt1
      cu(2,j,k1,l) = cu(2,j,k1,l) + dky*zt1
      cu(3,j,k1,l) = cu(3,j,k1,l) - dkz*zt1
      zt1 = at1*(dkx*cu(1,j,k,l1) + dky*cu(2,j,k,l1) - dkz*cu(3,j,k,l1))
      cu(1,j,k,l1) = cu(1,j,k,l1) - dkx*zt1
      cu(2,j,k,l1) = cu(2,j,k,l1) - dky*zt1
      cu(3,j,k,l1) = cu(3,j,k,l1) + dkz*zt1
      zt1 = at1*(dkx*cu(1,j,k1,l1) - dky*cu(2,j,k1,l1)
     1    - dkz*cu(3,j,k1,l1))
      cu(1,j,k1,l1) = cu(1,j,k1,l1) - dkx*zt1
      cu(2,j,k1,l1) = cu(2,j,k1,l1) + dky*zt1
      cu(3,j,k1,l1) = cu(3,j,k1,l1) + dkz*zt1
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
!dir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      at1 = 1.0/(dky*dky + dkz2)
      zt1 = at1*(dky*cu(2,1,k,l) + dkz*cu(3,1,k,l))
      cu(2,1,k,l) = cu(2,1,k,l) - dky*zt1
      cu(3,1,k,l) = cu(3,1,k,l) - dkz*zt1
      cu(1,1,k1,l) = zero
      cu(2,1,k1,l) = zero
      cu(3,1,k1,l) = zero
      zt1 = at1*(dky*cu(2,1,k,l1) - dkz*cu(3,1,k,l1))
      cu(2,1,k,l1) = cu(2,1,k,l1) - dky*zt1
      cu(3,1,k,l1) = cu(3,1,k,l1) + dkz*zt1
      cu(1,1,k1,l1) = zero
      cu(2,1,k1,l1) = zero
      cu(3,1,k1,l1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx + dkz2)
      zt1 = at1*(dkx*cu(1,j,1,l) + dkz*cu(3,j,1,l))
      cu(1,j,1,l) = cu(1,j,1,l) - dkx*zt1
      cu(3,j,1,l) = cu(3,j,1,l) - dkz*zt1
      cu(1,j,k1,l) = zero
      cu(2,j,k1,l) = zero
      cu(3,j,k1,l) = zero
      zt1 = at1*(dkx*cu(1,j,1,l1) - dkz*cu(3,j,1,l1))
      cu(1,j,1,l1) = cu(1,j,1,l1) - dkx*zt1
      cu(3,j,1,l1) = cu(3,j,1,l1) + dkz*zt1
      cu(1,j,k1,l1) = zero
      cu(2,j,k1,l1) = zero
      cu(3,j,k1,l1) = zero
   40 continue
c mode numbers kx = 0, nx/2
      cu(3,1,1,l) = zero
      cu(1,1,k1,l) = zero
      cu(2,1,k1,l) = zero
      cu(3,1,k1,l) = zero
      cu(1,1,1,l1) = zero
      cu(2,1,1,l1) = zero
      cu(3,1,1,l1) = zero
      cu(1,1,k1,l1) = zero
      cu(2,1,k1,l1) = zero
      cu(3,1,k1,l1) = zero
   50 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
!dir$ ivdep
      do 60 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      zt1 = at1*(dkx*cu(1,j,k,1) + dky*cu(2,j,k,1))
      cu(1,j,k,1) = cu(1,j,k,1) - dkx*zt1
      cu(2,j,k,1) = cu(2,j,k,1) - dky*zt1
      zt1 = at1*(dkx*cu(1,j,k1,1) - dky*cu(2,j,k1,1))
      cu(1,j,k1,1) = cu(1,j,k1,1) - dkx*zt1
      cu(2,j,k1,1) = cu(2,j,k1,1) + dky*zt1
      cu(1,j,k,l1) = zero
      cu(2,j,k,l1) = zero
      cu(3,j,k,l1) = zero
      cu(1,j,k1,l1) = zero
      cu(2,j,k1,l1) = zero
      cu(3,j,k1,l1) = zero
   60 continue
   70 continue
c mode numbers kx = 0, nx/2
!dir$ ivdep
      do 80 k = 2, nyh
      k1 = ny2 - k
      cu(2,1,k,1) = zero
      cu(1,1,k1,1) = zero
      cu(2,1,k1,1) = zero
      cu(3,1,k1,1) = zero
      cu(1,1,k,l1) = zero
      cu(2,1,k,l1) = zero
      cu(3,1,k,l1) = zero
      cu(1,1,k1,l1) = zero
      cu(2,1,k1,l1) = zero
      cu(3,1,k1,l1) = zero
   80 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 90 j = 2, nxh
      cu(1,j,1,1) = zero
      cu(1,j,k1,1) = zero
      cu(2,j,k1,1) = zero
      cu(3,j,k1,1) = zero
      cu(1,j,1,l1) = zero
      cu(2,j,1,l1) = zero
      cu(3,j,1,l1) = zero
      cu(1,j,k1,l1) = zero
      cu(2,j,k1,l1) = zero
      cu(3,j,k1,l1) = zero
   90 continue
      cu(1,1,1,1) = zero
      cu(2,1,1,1) = zero
      cu(3,1,1,1) = zero
      cu(1,1,k1,1) = zero
      cu(2,1,k1,1) = zero
      cu(3,1,k1,1) = zero
      cu(1,1,1,l1) = zero
      cu(2,1,1,l1) = zero
      cu(3,1,1,l1) = zero
      cu(1,1,k1,l1) = zero
      cu(2,1,k1,l1) = zero
      cu(3,1,k1,l1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine VIBPOIS33(cu,bxyz,ffc,ci,wm,nx,ny,nz,nxvh,nyv,nzv,nxhd,
     1nyhd,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c magnetic field with periodic boundary conditions.
c input: cu,ffc,ci,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
c output: bxyz, wm
c approximate flop count is:
c 193*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c the magnetic field is calculated using the equations:
c bx(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (ky*cuz(kx,ky,kz)-kz*cuy(kx,ky,kz)),
c by(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kz*cux(kx,ky,kz)-kx*cuz(kx,ky,kz)),
c bz(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kx*cuy(kx,ky,kz)-ky*cux(kx,ky,kz)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bx(ky=pi) = 0,
c bx(kz=pi) = by(kz=pi) = bz(kz=pi) = 0,
c bx(kx=0,ky=0,kz=0) = by(kx=0,ky=0,kz=0) = bz(kx=0,ky=0,kz=0) = 0.
c cu(i,j,k,l) = complex current density for fourier mode (j-1,k-1,l-1)
c bxyz(i,j,k,l) = i component of complex magnetic field
c all for fourier mode (j-1,k-1,l-1)
c aimag(ffc(j,k,l)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1,l-1)
c real(ffc(j,k,l)) = potential green's function g
c for fourier mode (j-1,k-1,l-1)
c ci = reciprocal of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny/nz = system length in x/y/z direction
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
c nzv = fourth dimension of field arrays, must be >= nz
c nxhd = dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
c nzhd = third dimension of form factor array, must be >= nzh
c vectorizable version
      implicit none
      integer nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
      real ci, wm
      complex cu, bxyz, ffc
      dimension cu(4,nxvh,nyv,nzv), bxyz(4,nxvh,nyv,nzv)
      dimension ffc(nxhd,nyhd,nzhd)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real dnx, dny, dnz, dky, dkz, ci2, at1, at2, at3, at4
      complex zero, zt1, zt2, zt3
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 50 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
!dir$ ivdep
      do 10 j = 2, nxh
      at1 = ci2*real(ffc(j,k,l))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      at4 = dkz*at1
      at1 = at1*aimag(ffc(j,k,l))
      zt1 = cmplx(-aimag(cu(3,j,k,l)),real(cu(3,j,k,l)))
      zt2 = cmplx(-aimag(cu(2,j,k,l)),real(cu(2,j,k,l)))
      zt3 = cmplx(-aimag(cu(1,j,k,l)),real(cu(1,j,k,l)))
      bxyz(1,j,k,l) = at3*zt1 - at4*zt2
      bxyz(2,j,k,l) = at4*zt3 - at2*zt1
      bxyz(3,j,k,l) = at2*zt2 - at3*zt3
      zt1 = cmplx(-aimag(cu(3,j,k1,l)),real(cu(3,j,k1,l)))
      zt2 = cmplx(-aimag(cu(2,j,k1,l)),real(cu(2,j,k1,l)))
      zt3 = cmplx(-aimag(cu(1,j,k1,l)),real(cu(1,j,k1,l)))
      bxyz(1,j,k1,l) = -at3*zt1 - at4*zt2
      bxyz(2,j,k1,l) = at4*zt3 - at2*zt1
      bxyz(3,j,k1,l) = at2*zt2 + at3*zt3
      zt1 = cmplx(-aimag(cu(3,j,k,l1)),real(cu(3,j,k,l1)))
      zt2 = cmplx(-aimag(cu(2,j,k,l1)),real(cu(2,j,k,l1)))
      zt3 = cmplx(-aimag(cu(1,j,k,l1)),real(cu(1,j,k,l1)))
      bxyz(1,j,k,l1) = at3*zt1 + at4*zt2
      bxyz(2,j,k,l1) = -at4*zt3 - at2*zt1
      bxyz(3,j,k,l1) = at2*zt2 - at3*zt3
      zt1 = cmplx(-aimag(cu(3,j,k1,l1)),real(cu(3,j,k1,l1)))
      zt2 = cmplx(-aimag(cu(2,j,k1,l1)),real(cu(2,j,k1,l1)))
      zt3 = cmplx(-aimag(cu(1,j,k1,l1)),real(cu(1,j,k1,l1)))
      bxyz(1,j,k1,l1) = -at3*zt1 + at4*zt2
      bxyz(2,j,k1,l1) = -at4*zt3 - at2*zt1
      bxyz(3,j,k1,l1) = at2*zt2 + at3*zt3
      at1 = at1*(cu(1,j,k,l)*conjg(cu(1,j,k,l))
     1   + cu(2,j,k,l)*conjg(cu(2,j,k,l))
     2   + cu(3,j,k,l)*conjg(cu(3,j,k,l))
     3   + cu(1,j,k1,l)*conjg(cu(1,j,k1,l))
     4   + cu(2,j,k1,l)*conjg(cu(2,j,k1,l))
     5   + cu(3,j,k1,l)*conjg(cu(3,j,k1,l))
     6   + cu(1,j,k,l1)*conjg(cu(1,j,k,l1))
     7   + cu(2,j,k,l1)*conjg(cu(2,j,k,l1))
     8   + cu(3,j,k,l1)*conjg(cu(3,j,k,l1))
     9   + cu(1,j,k1,l1)*conjg(cu(1,j,k1,l1))
     a   + cu(2,j,k1,l1)*conjg(cu(2,j,k1,l1))
     b   + cu(3,j,k1,l1)*conjg(cu(3,j,k1,l1)))
      wp = wp + dble(at1)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
!dir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k,l))
      at3 = dny*real(k - 1)*at1
      at4 = dkz*at1
      at1 = at1*aimag(ffc(1,k,l))
      zt1 = cmplx(-aimag(cu(3,1,k,l)),real(cu(3,1,k,l)))
      zt2 = cmplx(-aimag(cu(2,1,k,l)),real(cu(2,1,k,l)))
      zt3 = cmplx(-aimag(cu(1,1,k,l)),real(cu(1,1,k,l)))
      bxyz(1,1,k,l) = at3*zt1 - at4*zt2
      bxyz(2,1,k,l) = at4*zt3
      bxyz(3,1,k,l) = -at3*zt3
      bxyz(1,1,k1,l) = zero
      bxyz(2,1,k1,l) = zero
      bxyz(3,1,k1,l) = zero
      zt1 = cmplx(-aimag(cu(3,1,k,l1)),real(cu(3,1,k,l1)))
      zt2 = cmplx(-aimag(cu(2,1,k,l1)),real(cu(2,1,k,l1)))
      zt3 = cmplx(-aimag(cu(1,1,k,l1)),real(cu(1,1,k,l1)))
      bxyz(1,1,k,l1) = at3*zt1 + at4*zt2
      bxyz(2,1,k,l1) = -at4*zt3
      bxyz(3,1,k,l1) = -at3*zt3
      bxyz(1,1,k1,l1) = zero
      bxyz(2,1,k1,l1) = zero
      bxyz(3,1,k1,l1) = zero
      at1 = at1*(cu(1,1,k,l)*conjg(cu(1,1,k,l))
     1   + cu(2,1,k,l)*conjg(cu(2,1,k,l))
     2   + cu(3,1,k,l)*conjg(cu(3,1,k,l))
     3   + cu(1,1,k,l1)*conjg(cu(1,1,k,l1))
     4   + cu(2,1,k,l1)*conjg(cu(2,1,k,l1))
     5   + cu(3,1,k,l1)*conjg(cu(3,1,k,l1)))
      wp = wp + dble(at1)
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 40 j = 2, nxh
      at1 = ci2*real(ffc(j,1,l))
      at2 = dnx*real(j - 1)*at1
      at4 = dkz*at1
      at1 = at1*aimag(ffc(j,1,l))
      zt1 = cmplx(-aimag(cu(3,j,1,l)),real(cu(3,j,1,l)))
      zt2 = cmplx(-aimag(cu(2,j,1,l)),real(cu(2,j,1,l)))
      zt3 = cmplx(-aimag(cu(1,j,1,l)),real(cu(1,j,1,l)))
      bxyz(1,j,1,l) = -at4*zt2
      bxyz(2,j,1,l) = at4*zt3 - at2*zt1
      bxyz(3,j,1,l) = at2*zt2
      bxyz(1,j,k1,l) = zero
      bxyz(2,j,k1,l) = zero
      bxyz(3,j,k1,l) = zero
      zt1 = cmplx(-aimag(cu(3,j,1,l1)),real(cu(3,j,1,l1)))
      zt2 = cmplx(-aimag(cu(2,j,1,l1)),real(cu(2,j,1,l1)))
      zt3 = cmplx(-aimag(cu(1,j,1,l1)),real(cu(1,j,1,l1)))
      bxyz(1,j,1,l1) = at4*zt2
      bxyz(2,j,1,l1) = -at4*zt3 - at2*zt1
      bxyz(3,j,1,l1) = at2*zt2
      bxyz(1,j,k1,l1) = zero
      bxyz(2,j,k1,l1) = zero
      bxyz(3,j,k1,l1) = zero
      at1 = at1*(cu(1,j,1,l)*conjg(cu(1,j,1,l))
     1   + cu(2,j,1,l)*conjg(cu(2,j,1,l))
     2   + cu(3,j,1,l)*conjg(cu(3,j,1,l))
     3   + cu(1,j,1,l1)*conjg(cu(1,j,1,l1))
     4   + cu(2,j,1,l1)*conjg(cu(2,j,1,l1))
     5   + cu(3,j,1,l1)*conjg(cu(3,j,1,l1)))
      wp = wp + dble(at1)
   40 continue
c mode numbers kx = 0, nx/2
      at1 = ci2*real(ffc(1,1,l))
      at4 = dkz*at1
      at1 = at1*aimag(ffc(1,1,l))
      zt2 = cmplx(-aimag(cu(2,1,1,l)),real(cu(2,1,1,l)))
      zt3 = cmplx(-aimag(cu(1,1,1,l)),real(cu(1,1,1,l)))
      bxyz(1,1,1,l) = -at4*zt2
      bxyz(2,1,1,l) = at4*zt3
      bxyz(3,1,1,l) = zero
      bxyz(1,1,k1,l) = zero
      bxyz(2,1,k1,l) = zero
      bxyz(3,1,k1,l) = zero
      bxyz(1,1,1,l1) = zero
      bxyz(2,1,1,l1) = zero
      bxyz(3,1,1,l1) = zero
      bxyz(1,1,k1,l1) = zero
      bxyz(2,1,k1,l1) = zero
      bxyz(3,1,k1,l1) = zero
      at1 = at1*(cu(1,1,1,l)*conjg(cu(1,1,1,l))
     1   + cu(2,1,1,l)*conjg(cu(2,1,1,l))
     2   + cu(3,1,1,l)*conjg(cu(3,1,1,l)))
      wp = wp + dble(at1)
   50 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
!dir$ ivdep
      do 60 j = 2, nxh
      at1 = ci2*real(ffc(j,k,1))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      at1 = at1*aimag(ffc(j,k,1))
      zt1 = cmplx(-aimag(cu(3,j,k,1)),real(cu(3,j,k,1)))
      zt2 = cmplx(-aimag(cu(2,j,k,1)),real(cu(2,j,k,1)))
      zt3 = cmplx(-aimag(cu(1,j,k,1)),real(cu(1,j,k,1)))
      bxyz(1,j,k,1) = at3*zt1
      bxyz(2,j,k,1) = -at2*zt1
      bxyz(3,j,k,1) = at2*zt2 - at3*zt3
      zt1 = cmplx(-aimag(cu(3,j,k1,1)),real(cu(3,j,k1,1)))
      zt2 = cmplx(-aimag(cu(2,j,k1,1)),real(cu(2,j,k1,1)))
      zt3 = cmplx(-aimag(cu(1,j,k1,1)),real(cu(1,j,k1,1)))
      bxyz(1,j,k1,1) = -at3*zt1
      bxyz(2,j,k1,1) = -at2*zt1
      bxyz(3,j,k1,1) = at2*zt2 + at3*zt3
      bxyz(1,j,k,l1) = zero
      bxyz(2,j,k,l1) = zero
      bxyz(3,j,k,l1) = zero
      bxyz(1,j,k1,l1) = zero
      bxyz(2,j,k1,l1) = zero
      bxyz(3,j,k1,l1) = zero
      at1 = at1*(cu(1,j,k,1)*conjg(cu(1,j,k,1))
     1   + cu(2,j,k,1)*conjg(cu(2,j,k,1))
     2   + cu(3,j,k,1)*conjg(cu(3,j,k,1))
     3   + cu(1,j,k1,1)*conjg(cu(1,j,k1,1))
     4   + cu(2,j,k1,1)*conjg(cu(2,j,k1,1))
     5   + cu(3,j,k1,1)*conjg(cu(3,j,k1,1)))
      wp = wp + dble(at1)
   60 continue
   70 continue
c mode numbers kx = 0, nx/2
!dir$ ivdep
      do 80 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k,1))
      at3 = dny*real(k - 1)*at1
      at1 = at1*aimag(ffc(1,k,1))
      zt1 = cmplx(-aimag(cu(3,1,k,1)),real(cu(3,1,k,1)))
      zt3 = cmplx(-aimag(cu(1,1,k,1)),real(cu(1,1,k,1)))
      bxyz(1,1,k,1) = at3*zt1
      bxyz(2,1,k,1) = zero
      bxyz(3,1,k,1) = -at3*zt3
      bxyz(1,1,k1,1) = zero
      bxyz(2,1,k1,1) = zero
      bxyz(3,1,k1,1) = zero
      bxyz(1,1,k,l1) = zero
      bxyz(2,1,k,l1) = zero
      bxyz(3,1,k,l1) = zero
      bxyz(1,1,k1,l1) = zero
      bxyz(2,1,k1,l1) = zero
      bxyz(3,1,k1,l1) = zero
      at1 = at1*(cu(1,1,k,1)*conjg(cu(1,1,k,1))
     1   + cu(2,1,k,1)*conjg(cu(2,1,k,1))
     2   + cu(3,1,k,1)*conjg(cu(3,1,k,1)))
      wp = wp + dble(at1)
   80 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 90 j = 2, nxh
      at1 = ci2*real(ffc(j,1,1))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffc(j,1,1))
      zt1 = cmplx(-aimag(cu(3,j,1,1)),real(cu(3,j,1,1)))
      zt2 = cmplx(-aimag(cu(2,j,1,1)),real(cu(2,j,1,1)))
      bxyz(1,j,1,1) = zero
      bxyz(2,j,1,1) = -at2*zt1
      bxyz(3,j,1,1) = at2*zt2
      bxyz(1,j,k1,1) = zero
      bxyz(2,j,k1,1) = zero
      bxyz(3,j,k1,1) = zero
      bxyz(1,j,1,l1) = zero
      bxyz(2,j,1,l1) = zero
      bxyz(3,j,1,l1) = zero
      bxyz(1,j,k1,l1) = zero
      bxyz(2,j,k1,l1) = zero
      bxyz(3,j,k1,l1) = zero
      at1 = at1*(cu(1,j,1,1)*conjg(cu(1,j,1,1))
     1   + cu(2,j,1,1)*conjg(cu(2,j,1,1))
     2   + cu(3,j,1,1)*conjg(cu(3,j,1,1)))
      wp = wp + dble(at1)
   90 continue
      bxyz(1,1,1,1) = zero
      bxyz(2,1,1,1) = zero
      bxyz(3,1,1,1) = zero
      bxyz(1,1,k1,1) = zero
      bxyz(2,1,k1,1) = zero
      bxyz(3,1,k1,1) = zero
      bxyz(1,1,1,l1) = zero
      bxyz(2,1,1,l1) = zero
      bxyz(3,1,1,l1) = zero
      bxyz(1,1,k1,l1) = zero
      bxyz(2,1,k1,l1) = zero
      bxyz(3,1,k1,l1) = zero
      wm = real(nx)*real(ny)*real(nz)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine VMAXWEL3(exyz,bxyz,cu,ffc,ci,dt,wf,wm,nx,ny,nz,nxvh,nyv
     1,nzv,nxhd,nyhd,nzhd)
c this subroutine solves 3d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions.
c input: all, output: wf, wm, exyz, bxyz
c approximate flop count is:
c 680*nxc*nyc*nzc + 149*(nxc*nyc + nxc*nzc + nyc*nzc)
c plus nxc*nyc*nzc divides
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky,kz) = bx(kx,ky,kz) - .5*dt*sqrt(-1)*
c                (ky*ez(kx,ky,kz)-kz*ey(kx,ky,kz))
c by(kx,ky,kz) = by(kx,ky,kz) - .5*dt*sqrt(-1)*
c               (kz*ex(kx,ky,kz)-kx*ez(kx,ky,kz))
c bz(kx,ky,kz) = bz(kx,ky,kz) - .5*dt*sqrt(-1)*
c               (kx*ey(kx,ky,kz)-ky*ex(kx,ky,kz))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky,kz) = ex(kx,ky,kz) + c2*dt*sqrt(-1)*
c  (ky*bz(kx,ky,kz)-kz*by(kx,ky,kz)) - affp*dt*cux(kx,ky,kz)*s(kx,ky,kz)
c ey(kx,ky,kz) = ey(kx,ky,kz) + c2*dt*sqrt(-1)*
c  (kz*bx(kx,ky,kz)-kx*bz(kx,ky,kz)) - affp*dt*cuy(kx,ky,kz)*s(kx,ky,kz)
c ez(kx,ky,kz) = ez(kx,ky,kz) + c2*dt*sqrt(-1)*
c  (kx*by(kx,ky,kz)-ky*bx(kx,ky,kz)) - affp*dt*cuz(kx,ky,kz)*s(kx,ky,kz)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, c2 = 1./(ci*ci)
c and s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)
c j,k,l = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
c ex(kz=pi) = ey(kz=pi) = ez(kz=pi) = 0,
c ex(kx=0,ky=0,kz=0) = ey(kx=0,ky=0,kz=0) = ez(kx=0,ky=0,kz=0) = 0.
c and similarly for bx, by, bz.
c cu(i,j,k,l) = complex current density
c exyz(i,j,k,l) = complex transverse electric field
c bxyz(i,j,k,l) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1,l-1)
c real(ffc(1,1,1)) = affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c aimag(ffc(j,k,l)) = finite-size particle shape factor s,
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
c for fourier mode (j-1,k-1,l-1)
c ci = reciprocal of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
c nzv = fourth dimension of field arrays, must be >= nz
c nxhd = second dimension of form factor array, must be >= nxh
c nyhd = third dimension of form factor array, must be >= nyh
c nzhd = fourth dimension of form factor array, must be >= nzh
c vectorizable version
      implicit none
      integer nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
      real ci, dt, wf, wm
      complex exyz, bxyz, cu, ffc
      dimension exyz(4,nxvh,nyv,nzv), bxyz(4,nxvh,nyv,nzv)
      dimension cu(4,nxvh,nyv,nzv)
      dimension ffc(nxhd,nyhd,nzhd)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real dnx, dny, dnz, dth, c2, cdt, affp, anorm, dkx, dky, dkz
      real adt, afdt
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      real at1
      double precision wp, ws
      if (ci.le.0.0) return
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      dth = .5*dt
      c2 = 1.0/(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1,1,1))
      adt = affp*dt
      zero = cmplx(0.0,0.0)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
      do 50 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
!dir$ ivdep
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j,k,l))
c update magnetic field half time step, ky > 0, kz > 0
      zt1 = cmplx(-aimag(exyz(3,j,k,l)),real(exyz(3,j,k,l)))
      zt2 = cmplx(-aimag(exyz(2,j,k,l)),real(exyz(2,j,k,l)))
      zt3 = cmplx(-aimag(exyz(1,j,k,l)),real(exyz(1,j,k,l)))
      zt4 = bxyz(1,j,k,l) - dth*(dky*zt1 - dkz*zt2)
      zt5 = bxyz(2,j,k,l) - dth*(dkz*zt3 - dkx*zt1)
      zt6 = bxyz(3,j,k,l) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,j,k,l) + cdt*(dky*zt1 - dkz*zt2) - afdt*cu(1,j,k,l)
      zt8 = exyz(2,j,k,l) + cdt*(dkz*zt3 - dkx*zt1) - afdt*cu(2,j,k,l)
      zt9 = exyz(3,j,k,l) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,j,k,l)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,j,k,l) = zt7
      exyz(2,j,k,l) = zt8
      exyz(3,j,k,l) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 - dth*(dky*zt1 - dkz*zt2)
      zt5 = zt5 - dth*(dkz*zt3 - dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
      bxyz(1,j,k,l) = zt4
      bxyz(2,j,k,l) = zt5
      bxyz(3,j,k,l) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
c update magnetic field half time step, ky < 0, kz > 0
      zt1 = cmplx(-aimag(exyz(3,j,k1,l)),real(exyz(3,j,k1,l)))
      zt2 = cmplx(-aimag(exyz(2,j,k1,l)),real(exyz(2,j,k1,l)))
      zt3 = cmplx(-aimag(exyz(1,j,k1,l)),real(exyz(1,j,k1,l)))
      zt4 = bxyz(1,j,k1,l) + dth*(dky*zt1 + dkz*zt2)
      zt5 = bxyz(2,j,k1,l) - dth*(dkz*zt3 - dkx*zt1)
      zt6 = bxyz(3,j,k1,l) - dth*(dkx*zt2 + dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,j,k1,l) - cdt*(dky*zt1 + dkz*zt2) - afdt*cu(1,j,k1,l)
      zt8 = exyz(2,j,k1,l) + cdt*(dkz*zt3 - dkx*zt1) - afdt*cu(2,j,k1,l)
      zt9 = exyz(3,j,k1,l) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,j,k1,l)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,j,k1,l) = zt7
      exyz(2,j,k1,l) = zt8
      exyz(3,j,k1,l) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 + dth*(dky*zt1 + dkz*zt2)
      zt5 = zt5 - dth*(dkz*zt3 - dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
      bxyz(1,j,k1,l) = zt4
      bxyz(2,j,k1,l) = zt5
      bxyz(3,j,k1,l) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
c update magnetic field half time step, ky > 0, kz < 0
      zt1 = cmplx(-aimag(exyz(3,j,k,l1)),real(exyz(3,j,k,l1)))
      zt2 = cmplx(-aimag(exyz(2,j,k,l1)),real(exyz(2,j,k,l1)))
      zt3 = cmplx(-aimag(exyz(1,j,k,l1)),real(exyz(1,j,k,l1)))
      zt4 = bxyz(1,j,k,l1) - dth*(dky*zt1 + dkz*zt2)
      zt5 = bxyz(2,j,k,l1) + dth*(dkz*zt3 + dkx*zt1)
      zt6 = bxyz(3,j,k,l1) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,j,k,l1) + cdt*(dky*zt1 + dkz*zt2) - afdt*cu(1,j,k,l1)
      zt8 = exyz(2,j,k,l1) - cdt*(dkz*zt3 + dkx*zt1) - afdt*cu(2,j,k,l1)
      zt9 = exyz(3,j,k,l1) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,j,k,l1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,j,k,l1) = zt7
      exyz(2,j,k,l1) = zt8
      exyz(3,j,k,l1) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 - dth*(dky*zt1 + dkz*zt2)
      zt5 = zt5 + dth*(dkz*zt3 + dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
      bxyz(1,j,k,l1) = zt4
      bxyz(2,j,k,l1) = zt5
      bxyz(3,j,k,l1) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
c update magnetic field half time step, ky < 0, kz < 0
      zt1 = cmplx(-aimag(exyz(3,j,k1,l1)),real(exyz(3,j,k1,l1)))
      zt2 = cmplx(-aimag(exyz(2,j,k1,l1)),real(exyz(2,j,k1,l1)))
      zt3 = cmplx(-aimag(exyz(1,j,k1,l1)),real(exyz(1,j,k1,l1)))
      zt4 = bxyz(1,j,k1,l1) + dth*(dky*zt1 - dkz*zt2)
      zt5 = bxyz(2,j,k1,l1) + dth*(dkz*zt3 + dkx*zt1)
      zt6 = bxyz(3,j,k1,l1) - dth*(dkx*zt2 + dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,j,k1,l1) - cdt*(dky*zt1 - dkz*zt2)
     1    - afdt*cu(1,j,k1,l1)
      zt8 = exyz(2,j,k1,l1) - cdt*(dkz*zt3 + dkx*zt1)
     1    - afdt*cu(2,j,k1,l1)
      zt9 = exyz(3,j,k1,l1) + cdt*(dkx*zt2 + dky*zt3)
     1    - afdt*cu(3,j,k1,l1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,j,k1,l1) = zt7
      exyz(2,j,k1,l1) = zt8
      exyz(3,j,k1,l1) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 + dth*(dky*zt1 - dkz*zt2)
      zt5 = zt5 + dth*(dkz*zt3 + dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
      bxyz(1,j,k1,l1) = zt4
      bxyz(2,j,k1,l1) = zt5
      bxyz(3,j,k1,l1) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
!dir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      afdt = adt*aimag(ffc(1,k,l))
c update magnetic field half time step, kz > 0
      zt1 = cmplx(-aimag(exyz(3,1,k,l)),real(exyz(3,1,k,l)))
      zt2 = cmplx(-aimag(exyz(2,1,k,l)),real(exyz(2,1,k,l)))
      zt3 = cmplx(-aimag(exyz(1,1,k,l)),real(exyz(1,1,k,l)))
      zt4 = bxyz(1,1,k,l) - dth*(dky*zt1 - dkz*zt2)
      zt5 = bxyz(2,1,k,l) - dth*(dkz*zt3)
      zt6 = bxyz(3,1,k,l) + dth*(dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,1,k,l) + cdt*(dky*zt1 - dkz*zt2) - afdt*cu(1,1,k,l)
      zt8 = exyz(2,1,k,l) + cdt*(dkz*zt3) - afdt*cu(2,1,k,l)
      zt9 = exyz(3,1,k,l) - cdt*(dky*zt3) - afdt*cu(3,1,k,l)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,1,k,l) = zt7
      exyz(2,1,k,l) = zt8
      exyz(3,1,k,l) = zt9  
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 - dth*(dky*zt1 - dkz*zt2)
      zt5 = zt5 - dth*(dkz*zt3)
      zt6 = zt6 + dth*(dky*zt3) 
      bxyz(1,1,k,l) = zt4
      bxyz(2,1,k,l) = zt5
      bxyz(3,1,k,l) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
      bxyz(1,1,k1,l) = zero
      bxyz(2,1,k1,l) = zero
      bxyz(3,1,k1,l) = zero
      exyz(1,1,k1,l) = zero
      exyz(2,1,k1,l) = zero
      exyz(3,1,k1,l) = zero
c update magnetic field half time step, kz < 0
      zt1 = cmplx(-aimag(exyz(3,1,k,l1)),real(exyz(3,1,k,l1)))
      zt2 = cmplx(-aimag(exyz(2,1,k,l1)),real(exyz(2,1,k,l1)))
      zt3 = cmplx(-aimag(exyz(1,1,k,l1)),real(exyz(1,1,k,l1)))
      zt4 = bxyz(1,1,k,l1) - dth*(dky*zt1 + dkz*zt2)
      zt5 = bxyz(2,1,k,l1) + dth*(dkz*zt3)
      zt6 = bxyz(3,1,k,l1) + dth*(dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,1,k,l1) + cdt*(dky*zt1 + dkz*zt2) - afdt*cu(1,1,k,l1)
      zt8 = exyz(2,1,k,l1) - cdt*(dkz*zt3) - afdt*cu(2,1,k,l1)
      zt9 = exyz(3,1,k,l1) - cdt*(dky*zt3) - afdt*cu(3,1,k,l1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,1,k,l1) = zt7
      exyz(2,1,k,l1) = zt8
      exyz(3,1,k,l1) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 - dth*(dky*zt1 + dkz*zt2)
      zt5 = zt5 + dth*(dkz*zt3)
      zt6 = zt6 + dth*(dky*zt3)
      bxyz(1,1,k,l1) = zt4
      bxyz(2,1,k,l1) = zt5
      bxyz(3,1,k,l1) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
      bxyz(1,1,k1,l1) = zero
      bxyz(2,1,k1,l1) = zero
      bxyz(3,1,k1,l1) = zero
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j,1,l))
c update magnetic field half time step, kz > 0
      zt1 = cmplx(-aimag(exyz(3,j,1,l)),real(exyz(3,j,1,l)))
      zt2 = cmplx(-aimag(exyz(2,j,1,l)),real(exyz(2,j,1,l)))
      zt3 = cmplx(-aimag(exyz(1,j,1,l)),real(exyz(1,j,1,l)))
      zt4 = bxyz(1,j,1,l) + dth*(dkz*zt2)
      zt5 = bxyz(2,j,1,l) - dth*(dkz*zt3 - dkx*zt1)
      zt6 = bxyz(3,j,1,l) - dth*(dkx*zt2)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,j,1,l) - cdt*(dkz*zt2) - afdt*cu(1,j,1,l)
      zt8 = exyz(2,j,1,l) + cdt*(dkz*zt3 - dkx*zt1) - afdt*cu(2,j,1,l)
      zt9 = exyz(3,j,1,l) + cdt*(dkx*zt2) - afdt*cu(3,j,1,l)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,j,1,l) = zt7
      exyz(2,j,1,l) = zt8
      exyz(3,j,1,l) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 + dth*(dkz*zt2)
      zt5 = zt5 - dth*(dkz*zt3 - dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      bxyz(1,j,1,l) = zt4
      bxyz(2,j,1,l) = zt5
      bxyz(3,j,1,l) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
      bxyz(1,j,k1,l) = zero
      bxyz(2,j,k1,l) = zero
      bxyz(3,j,k1,l) = zero
      exyz(1,j,k1,l) = zero
      exyz(2,j,k1,l) = zero
      exyz(3,j,k1,l) = zero
c update magnetic field half time step, kz > 0
      zt1 = cmplx(-aimag(exyz(3,j,1,l1)),real(exyz(3,j,1,l1)))
      zt2 = cmplx(-aimag(exyz(2,j,1,l1)),real(exyz(2,j,1,l1)))
      zt3 = cmplx(-aimag(exyz(1,j,1,l1)),real(exyz(1,j,1,l1)))
      zt4 = bxyz(1,j,1,l1) - dth*(dkz*zt2)
      zt5 = bxyz(2,j,1,l1) + dth*(dkz*zt3 + dkx*zt1)
      zt6 = bxyz(3,j,1,l1) - dth*(dkx*zt2)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,j,1,l1) + cdt*(dkz*zt2) - afdt*cu(1,j,1,l1)
      zt8 = exyz(2,j,1,l1) - cdt*(dkz*zt3 + dkx*zt1) - afdt*cu(2,j,1,l1)
      zt9 = exyz(3,j,1,l1) + cdt*(dkx*zt2) - afdt*cu(3,j,1,l1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,j,1,l1) = zt7
      exyz(2,j,1,l1) = zt8
      exyz(3,j,1,l1) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 - dth*(dkz*zt2)
      zt5 = zt5 + dth*(dkz*zt3 + dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      bxyz(1,j,1,l1) = zt4
      bxyz(2,j,1,l1) = zt5
      bxyz(3,j,1,l1) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
      bxyz(1,j,k1,l1) = zero
      bxyz(2,j,k1,l1) = zero
      bxyz(3,j,k1,l1) = zero
      exyz(1,j,k1,l1) = zero
      exyz(2,j,k1,l1) = zero
      exyz(3,j,k1,l1) = zero
   40 continue
c mode numbers kx = 0, nx/2
      afdt = adt*aimag(ffc(1,1,l))
c update magnetic field half time step
      zt2 = cmplx(-aimag(exyz(2,1,1,l)),real(exyz(2,1,1,l)))
      zt3 = cmplx(-aimag(exyz(1,1,1,l)),real(exyz(1,1,1,l)))
      zt4 = bxyz(1,1,1,l) + dth*(dkz*zt2)
      zt5 = bxyz(2,1,1,l) - dth*(dkz*zt3)
c update electric field whole time step
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,1,1,l) - cdt*(dkz*zt2) - afdt*cu(1,1,1,l)
      zt8 = exyz(2,1,1,l) + cdt*(dkz*zt3) - afdt*cu(2,1,1,l)
c update magnetic field half time step and store electric field
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,1,1,l) = zt7
      exyz(2,1,1,l) = zt8
      exyz(3,1,1,l) = zero
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8))
      ws = ws + dble(at1)
      zt4 = zt4 + dth*(dkz*zt2)
      zt5 = zt5 - dth*(dkz*zt3)
      bxyz(1,1,1,l) = zt4
      bxyz(2,1,1,l) = zt5
      bxyz(3,1,1,l) = zero
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5))
      wp = wp + dble(at1)
      bxyz(1,1,k1,l) = zero
      bxyz(2,1,k1,l) = zero
      bxyz(3,1,k1,l) = zero
      exyz(1,1,k1,l) = zero
      exyz(2,1,k1,l) = zero
      exyz(3,1,k1,l) = zero
      bxyz(1,1,1,l1) = zero
      bxyz(2,1,1,l1) = zero
      bxyz(3,1,1,l1) = zero
      exyz(1,1,1,l1) = zero
      exyz(2,1,1,l1) = zero
      exyz(3,1,1,l1) = zero
      bxyz(1,1,k1,l1) = zero
      bxyz(2,1,k1,l1) = zero
      bxyz(3,1,k1,l1) = zero
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
   50 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
!dir$ ivdep
      do 60 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j,k,1))
c update magnetic field half time step, ky > 0
      zt1 = cmplx(-aimag(exyz(3,j,k,1)),real(exyz(3,j,k,1)))
      zt2 = cmplx(-aimag(exyz(2,j,k,1)),real(exyz(2,j,k,1)))
      zt3 = cmplx(-aimag(exyz(1,j,k,1)),real(exyz(1,j,k,1)))
      zt4 = bxyz(1,j,k,1) - dth*(dky*zt1)
      zt5 = bxyz(2,j,k,1) + dth*(dkx*zt1)
      zt6 = bxyz(3,j,k,1) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,j,k,1) + cdt*(dky*zt1) - afdt*cu(1,j,k,1)
      zt8 = exyz(2,j,k,1) - cdt*(dkx*zt1) - afdt*cu(2,j,k,1)
      zt9 = exyz(3,j,k,1) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,j,k,1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,j,k,1) = zt7
      exyz(2,j,k,1) = zt8
      exyz(3,j,k,1) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 - dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
      bxyz(1,j,k,1) = zt4
      bxyz(2,j,k,1) = zt5
      bxyz(3,j,k,1) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
c update magnetic field half time step, ky < 0
      zt1 = cmplx(-aimag(exyz(3,j,k1,1)),real(exyz(3,j,k1,1)))
      zt2 = cmplx(-aimag(exyz(2,j,k1,1)),real(exyz(2,j,k1,1)))
      zt3 = cmplx(-aimag(exyz(1,j,k1,1)),real(exyz(1,j,k1,1)))
      zt4 = bxyz(1,j,k1,1) + dth*(dky*zt1)
      zt5 = bxyz(2,j,k1,1) + dth*(dkx*zt1)
      zt6 = bxyz(3,j,k1,1) - dth*(dkx*zt2 + dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,j,k1,1) - cdt*(dky*zt1) - afdt*cu(1,j,k1,1)
      zt8 = exyz(2,j,k1,1) - cdt*(dkx*zt1) - afdt*cu(2,j,k1,1)
      zt9 = exyz(3,j,k1,1) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,j,k1,1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,j,k1,1) = zt7
      exyz(2,j,k1,1) = zt8
      exyz(3,j,k1,1) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 + dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
      bxyz(1,j,k1,1) = zt4
      bxyz(2,j,k1,1) = zt5
      bxyz(3,j,k1,1) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
      bxyz(1,j,k,l1) = zero
      bxyz(2,j,k,l1) = zero
      bxyz(3,j,k,l1) = zero
      exyz(1,j,k,l1) = zero
      exyz(2,j,k,l1) = zero
      exyz(3,j,k,l1) = zero
      bxyz(1,j,k1,l1) = zero
      bxyz(2,j,k1,l1) = zero
      bxyz(3,j,k1,l1) = zero
      exyz(1,j,k1,l1) = zero
      exyz(2,j,k1,l1) = zero
      exyz(3,j,k1,l1) = zero
   60 continue
   70 continue
c mode numbers kx = 0, nx/2
!dir$ ivdep
      do 80 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      afdt = adt*aimag(ffc(1,k,1))
c update magnetic field half time step
      zt1 = cmplx(-aimag(exyz(3,1,k,1)),real(exyz(3,1,k,1)))
      zt3 = cmplx(-aimag(exyz(1,1,k,1)),real(exyz(1,1,k,1)))
      zt4 = bxyz(1,1,k,1) - dth*(dky*zt1)
      zt6 = bxyz(3,1,k,1) + dth*(dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyz(1,1,k,1) + cdt*(dky*zt1) - afdt*cu(1,1,k,1)
      zt9 = exyz(3,1,k,1) - cdt*(dky*zt3) - afdt*cu(3,1,k,1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyz(1,1,k,1) = zt7
      exyz(2,1,k,1) = zero
      exyz(3,1,k,1) = zt9
      at1 = anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt4 = zt4 - dth*(dky*zt1)
      zt6 = zt6 + dth*(dky*zt3)
      bxyz(1,1,k,1) = zt4
      bxyz(2,1,k,1) = zero
      bxyz(3,1,k,1) = zt6
      at1 = anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
      wp = wp + dble(at1)
      bxyz(1,1,k1,1) = zero
      bxyz(2,1,k1,1) = zero
      bxyz(3,1,k1,1) = zero
      exyz(1,1,k1,1) = zero
      exyz(2,1,k1,1) = zero
      exyz(3,1,k1,1) = zero
      bxyz(1,1,k,l1) = zero
      bxyz(2,1,k,l1) = zero
      bxyz(3,1,k,l1) = zero
      exyz(1,1,k,l1) = zero
      exyz(2,1,k,l1) = zero
      exyz(3,1,k,l1) = zero
      bxyz(1,1,k1,l1) = zero
      bxyz(2,1,k1,l1) = zero
      bxyz(3,1,k1,l1) = zero
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
   80 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 90 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j,1,1))
c update magnetic field half time step
      zt1 = cmplx(-aimag(exyz(3,j,1,1)),real(exyz(3,j,1,1)))
      zt2 = cmplx(-aimag(exyz(2,j,1,1)),real(exyz(2,j,1,1)))
      zt5 = bxyz(2,j,1,1) + dth*(dkx*zt1)
      zt6 = bxyz(3,j,1,1) - dth*(dkx*zt2)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt8 = exyz(2,j,1,1) - cdt*(dkx*zt1) - afdt*cu(2,j,1,1)
      zt9 = exyz(3,j,1,1) + cdt*(dkx*zt2) - afdt*cu(3,j,1,1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      exyz(1,j,1,1) = zero
      exyz(2,j,1,1) = zt8
      exyz(3,j,1,1) = zt9
      at1 = anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      ws = ws + dble(at1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      bxyz(1,j,1,1) = zero
      bxyz(2,j,1,1) = zt5
      bxyz(3,j,1,1) = zt6
      at1 = anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
      wp = wp + dble(at1)
      bxyz(1,j,k1,1) = zero
      bxyz(2,j,k1,1) = zero
      bxyz(3,j,k1,1) = zero
      exyz(1,j,k1,1) = zero
      exyz(2,j,k1,1) = zero
      exyz(3,j,k1,1) = zero
      bxyz(1,j,1,l1) = zero
      bxyz(2,j,1,l1) = zero
      bxyz(3,j,1,l1) = zero
      exyz(1,j,1,l1) = zero
      exyz(2,j,1,l1) = zero
      exyz(3,j,1,l1) = zero
      bxyz(1,j,k1,l1) = zero
      bxyz(2,j,k1,l1) = zero
      bxyz(3,j,k1,l1) = zero
      exyz(1,j,k1,l1) = zero
      exyz(2,j,k1,l1) = zero
      exyz(3,j,k1,l1) = zero
   90 continue
      bxyz(1,1,1,1) = zero
      bxyz(2,1,1,1) = zero
      bxyz(3,1,1,1) = zero
      exyz(1,1,1,1) = zero
      exyz(2,1,1,1) = zero
      exyz(3,1,1,1) = zero
      bxyz(1,1,k1,1) = zero
      bxyz(2,1,k1,1) = zero
      bxyz(3,1,k1,1) = zero
      exyz(1,1,k1,1) = zero
      exyz(2,1,k1,1) = zero
      exyz(3,1,k1,1) = zero
      bxyz(1,1,1,l1) = zero
      bxyz(2,1,1,l1) = zero
      bxyz(3,1,1,l1) = zero
      exyz(1,1,1,l1) = zero
      exyz(2,1,1,l1) = zero
      exyz(3,1,1,l1) = zero
      bxyz(1,1,k1,l1) = zero
      bxyz(2,1,k1,l1) = zero
      bxyz(3,1,k1,l1) = zero
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
      wf = real(nx)*real(ny)*real(nz)*ws
      wm = real(nx)*real(ny)*real(nz)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine VEMFIELD3(fxyz,exyz,ffc,isign,nx,ny,nz,nxvh,nyv,nzv,   
     1nxhd,nyhd,nzhd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign < 0
c includes additional smoothing
      implicit none
      integer isign, nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
      complex fxyz, exyz, ffc
      dimension fxyz(4,nxvh,nyv,nzv), exyz(4,nxvh,nyv,nzv)
      dimension ffc(nxhd,nyhd,nzhd)
c local data
      integer j, k, l, nxh, nyh, nzh, ny2, nz2, k1, l1
      real at1
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
c add the fields
      if (isign.gt.0) then
         do 40 l = 2, nzh
         l1 = nz2 - l
         do 20 k = 2, nyh
         k1 = ny2 - k
!dir$ ivdep
         do 10 j = 1, nxh
         at1 = aimag(ffc(j,k,l))
         fxyz(1,j,k,l) = fxyz(1,j,k,l) + exyz(1,j,k,l)*at1
         fxyz(2,j,k,l) = fxyz(2,j,k,l) + exyz(2,j,k,l)*at1
         fxyz(3,j,k,l) = fxyz(3,j,k,l) + exyz(3,j,k,l)*at1
         fxyz(1,j,k1,l) = fxyz(1,j,k1,l) + exyz(1,j,k1,l)*at1
         fxyz(2,j,k1,l) = fxyz(2,j,k1,l) + exyz(2,j,k1,l)*at1
         fxyz(3,j,k1,l) = fxyz(3,j,k1,l) + exyz(3,j,k1,l)*at1
         fxyz(1,j,k,l1) = fxyz(1,j,k,l1) + exyz(1,j,k,l1)*at1
         fxyz(2,j,k,l1) = fxyz(2,j,k,l1) + exyz(2,j,k,l1)*at1
         fxyz(3,j,k,l1) = fxyz(3,j,k,l1) + exyz(3,j,k,l1)*at1
         fxyz(1,j,k1,l1) = fxyz(1,j,k1,l1) + exyz(1,j,k1,l1)*at1
         fxyz(2,j,k1,l1) = fxyz(2,j,k1,l1) + exyz(2,j,k1,l1)*at1
         fxyz(3,j,k1,l1) = fxyz(3,j,k1,l1) + exyz(3,j,k1,l1)*at1
   10    continue
   20    continue
         k1 = nyh + 1
!dir$ ivdep
         do 30 j = 1, nxh
         at1 = aimag(ffc(j,1,l))
         fxyz(1,j,1,l) = fxyz(1,j,1,l) + exyz(1,j,1,l)*at1
         fxyz(2,j,1,l) = fxyz(2,j,1,l) + exyz(2,j,1,l)*at1
         fxyz(3,j,1,l) = fxyz(3,j,1,l) + exyz(3,j,1,l)*at1
         fxyz(1,j,k1,l) = fxyz(1,j,k1,l) + exyz(1,j,k1,l)*at1
         fxyz(2,j,k1,l) = fxyz(2,j,k1,l) + exyz(2,j,k1,l)*at1
         fxyz(3,j,k1,l) = fxyz(3,j,k1,l) + exyz(3,j,k1,l)*at1
         fxyz(1,j,1,l1) = fxyz(1,j,1,l1) + exyz(1,j,1,l1)*at1
         fxyz(2,j,1,l1) = fxyz(2,j,1,l1) + exyz(2,j,1,l1)*at1
         fxyz(3,j,1,l1) = fxyz(3,j,1,l1) + exyz(3,j,1,l1)*at1
         fxyz(1,j,k1,l1) = fxyz(1,j,k1,l1) + exyz(1,j,k1,l1)*at1
         fxyz(2,j,k1,l1) = fxyz(2,j,k1,l1) + exyz(2,j,k1,l1)*at1
         fxyz(3,j,k1,l1) = fxyz(3,j,k1,l1) + exyz(3,j,k1,l1)*at1
   30    continue
   40    continue
         l1 = nzh + 1
         do 60 k = 2, nyh
         k1 = ny2 - k
!dir$ ivdep
         do 50 j = 1, nxh
         at1 = aimag(ffc(j,k,1))
         fxyz(1,j,k,1) = fxyz(1,j,k,1) + exyz(1,j,k,1)*at1
         fxyz(2,j,k,1) = fxyz(2,j,k,1) + exyz(2,j,k,1)*at1
         fxyz(3,j,k,1) = fxyz(3,j,k,1) + exyz(3,j,k,1)*at1
         fxyz(1,j,k1,1) = fxyz(1,j,k1,1) + exyz(1,j,k1,1)*at1
         fxyz(2,j,k1,1) = fxyz(2,j,k1,1) + exyz(2,j,k1,1)*at1
         fxyz(3,j,k1,1) = fxyz(3,j,k1,1) + exyz(3,j,k1,1)*at1
         fxyz(1,j,k,l1) = fxyz(1,j,k,l1) + exyz(1,j,k,l1)*at1
         fxyz(2,j,k,l1) = fxyz(2,j,k,l1) + exyz(2,j,k,l1)*at1
         fxyz(3,j,k,l1) = fxyz(3,j,k,l1) + exyz(3,j,k,l1)*at1
         fxyz(1,j,k1,l1) = fxyz(1,j,k1,l1) + exyz(1,j,k1,l1)*at1
         fxyz(2,j,k1,l1) = fxyz(2,j,k1,l1) + exyz(2,j,k1,l1)*at1
         fxyz(3,j,k1,l1) = fxyz(3,j,k1,l1) + exyz(3,j,k1,l1)*at1
   50    continue
   60    continue
         k1 = nyh + 1
!dir$ ivdep
         do 70 j = 1, nxh
         at1 = aimag(ffc(j,1,1))
         fxyz(1,j,1,1) = fxyz(1,j,1,1) + exyz(1,j,1,1)*at1
         fxyz(2,j,1,1) = fxyz(2,j,1,1) + exyz(2,j,1,1)*at1
         fxyz(3,j,1,1) = fxyz(3,j,1,1) + exyz(3,j,1,1)*at1
         fxyz(1,j,k1,1) = fxyz(1,j,k1,1) + exyz(1,j,k1,1)*at1
         fxyz(2,j,k1,1) = fxyz(2,j,k1,1) + exyz(2,j,k1,1)*at1
         fxyz(3,j,k1,1) = fxyz(3,j,k1,1) + exyz(3,j,k1,1)*at1
         fxyz(1,j,1,l1) = fxyz(1,j,1,l1) + exyz(1,j,1,l1)*at1
         fxyz(2,j,1,l1) = fxyz(2,j,1,l1) + exyz(2,j,1,l1)*at1
         fxyz(3,j,1,l1) = fxyz(3,j,1,l1) + exyz(3,j,1,l1)*at1
         fxyz(1,j,k1,l1) = fxyz(1,j,k1,l1) + exyz(1,j,k1,l1)*at1
         fxyz(2,j,k1,l1) = fxyz(2,j,k1,l1) + exyz(2,j,k1,l1)*at1
         fxyz(3,j,k1,l1) = fxyz(3,j,k1,l1) + exyz(3,j,k1,l1)*at1
   70    continue
c copy the fields
      else if (isign.lt.0) then
         do 110 l = 2, nzh
         l1 = nz2 - l
         do 90 k = 2, nyh
         k1 = ny2 - k
!dir$ ivdep
         do 80 j = 1, nxh
         at1 = aimag(ffc(j,k,l))
         fxyz(1,j,k,l) = exyz(1,j,k,l)*at1
         fxyz(2,j,k,l) = exyz(2,j,k,l)*at1
         fxyz(3,j,k,l) = exyz(3,j,k,l)*at1
         fxyz(1,j,k1,l) = exyz(1,j,k1,l)*at1
         fxyz(2,j,k1,l) = exyz(2,j,k1,l)*at1
         fxyz(3,j,k1,l) = exyz(3,j,k1,l)*at1
         fxyz(1,j,k,l1) = exyz(1,j,k,l1)*at1
         fxyz(2,j,k,l1) = exyz(2,j,k,l1)*at1
         fxyz(3,j,k,l1) = exyz(3,j,k,l1)*at1
         fxyz(1,j,k1,l1) = exyz(1,j,k1,l1)*at1
         fxyz(2,j,k1,l1) = exyz(2,j,k1,l1)*at1
         fxyz(3,j,k1,l1) = exyz(3,j,k1,l1)*at1
   80    continue
   90    continue
         k1 = nyh + 1
!dir$ ivdep
         do 100 j = 1, nxh
         at1 = aimag(ffc(j,1,l))
         fxyz(1,j,1,l) = exyz(1,j,1,l)*at1
         fxyz(2,j,1,l) = exyz(2,j,1,l)*at1
         fxyz(3,j,1,l) = exyz(3,j,1,l)*at1
         fxyz(1,j,k1,l) = exyz(1,j,k1,l)*at1
         fxyz(2,j,k1,l) = exyz(2,j,k1,l)*at1
         fxyz(3,j,k1,l) = exyz(3,j,k1,l)*at1
         fxyz(1,j,1,l1) = exyz(1,j,1,l1)*at1
         fxyz(2,j,1,l1) = exyz(2,j,1,l1)*at1
         fxyz(3,j,1,l1) = exyz(3,j,1,l1)*at1
         fxyz(1,j,k1,l1) = exyz(1,j,k1,l1)*at1
         fxyz(2,j,k1,l1) = exyz(2,j,k1,l1)*at1
         fxyz(3,j,k1,l1) = exyz(3,j,k1,l1)*at1
  100    continue
  110    continue
         l1 = nzh + 1
         do 130 k = 2, nyh
         k1 = ny2 - k
!dir$ ivdep
         do 120 j = 1, nxh
         at1 = aimag(ffc(j,k,1))
         fxyz(1,j,k,1) = exyz(1,j,k,1)*at1
         fxyz(2,j,k,1) = exyz(2,j,k,1)*at1
         fxyz(3,j,k,1) = exyz(3,j,k,1)*at1
         fxyz(1,j,k1,1) = exyz(1,j,k1,1)*at1
         fxyz(2,j,k1,1) = exyz(2,j,k1,1)*at1
         fxyz(3,j,k1,1) = exyz(3,j,k1,1)*at1
         fxyz(1,j,k,l1) = exyz(1,j,k,l1)*at1
         fxyz(2,j,k,l1) = exyz(2,j,k,l1)*at1
         fxyz(3,j,k,l1) = exyz(3,j,k,l1)*at1
         fxyz(1,j,k1,l1) = exyz(1,j,k1,l1)*at1
         fxyz(2,j,k1,l1) = exyz(2,j,k1,l1)*at1
         fxyz(3,j,k1,l1) = exyz(3,j,k1,l1)*at1
  120    continue
  130    continue
         k1 = nyh + 1
!dir$ ivdep
         do 140 j = 1, nxh
         at1 = aimag(ffc(j,1,1))
         fxyz(1,j,1,1) = exyz(1,j,1,1)*at1
         fxyz(2,j,1,1) = exyz(2,j,1,1)*at1
         fxyz(3,j,1,1) = exyz(3,j,1,1)*at1
         fxyz(1,j,k1,1) = exyz(1,j,k1,1)*at1
         fxyz(2,j,k1,1) = exyz(2,j,k1,1)*at1
         fxyz(3,j,k1,1) = exyz(3,j,k1,1)*at1
         fxyz(1,j,1,l1) = exyz(1,j,1,l1)*at1
         fxyz(2,j,1,l1) = exyz(2,j,1,l1)*at1
         fxyz(3,j,1,l1) = exyz(3,j,1,l1)*at1
         fxyz(1,j,k1,l1) = exyz(1,j,k1,l1)*at1
         fxyz(2,j,k1,l1) = exyz(2,j,k1,l1)*at1
         fxyz(3,j,k1,l1) = exyz(3,j,k1,l1)*at1
  140    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
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
      subroutine WFFT3RVX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd,
     1nxhyzd,nxyzhd)
c wrapper function for real to complex fft
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, indz, nxhd, nyd, nzd, nxhyzd, nxyzhd
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer ny, nz, nyi, nzi
      data nyi, nzi /1,1/
c calculate range of indices
      ny = 2**indy
      nz = 2**indz
c inverse fourier transform
      if (isign.lt.0) then
c perform xy fft
         call FFT3RVXY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
     1nzd,nxhyzd,nxyzhd)
c perform z fft
         call FFT3RXZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd, 
     1nzd,nxhyzd,nxyzhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform z fft
         call FFT3RXZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd, 
     1nzd,nxhyzd,nxyzhd)
c perform xy fft
         call FFT3RVXY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
     1nzd,nxhyzd,nxyzhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT3RV3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd,
     1nxhyzd,nxyzhd)
c wrapper function for 3 2d real to complex ffts
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, indz, nxhd, nyd, nzd, nxhyzd, nxyzhd
      dimension f(4,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer ny, nz, nyi, nzi
      data nyi, nzi /1,1/
c calculate range of indices
      ny = 2**indy
      nz = 2**indz
c inverse fourier transform
      if (isign.lt.0) then
c perform xy fft
         call FFT3RV3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd
     1,nzd,nxhyzd,nxyzhd)
c perform z fft
         call FFT3RV3Z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,
     1nzd,nxhyzd,nxyzhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform z fft
         call FFT3RV3Z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,
     1nzd,nxhyzd,nxyzhd)
c perform xy fft
         call FFT3RV3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd
     1,nzd,nxhyzd,nxyzhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RVXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,nxhd,
     1nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the x-y part of a three dimensional real to
c complex fast fourier transform and its inverse, for a subset of z,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c f(n,m,i) = (1/nx*ny*nz)*sum(f(j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)*
c       exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k,l) = sum(f(n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nzi = initial z index used
c nzp = number of z indices used
c nxhd = first dimension of f
c nyd,nzd = second and third dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1,1,1)) = real part of mode nx/2,0,0
c imag(f(1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nxyz, nxhyz, nzt, nrx, nry
      integer i, j, k, l, n, j1, j2, k1, k2, ns, ns2, km, kmr
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nzt = nzi + nzp - 1
      if (isign.gt.0) go to 180
c inverse fourier transform
      do 170 n = nzi, nzt
c bit-reverse array elements in x
      nrx = nxhyz/nxh
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 i = 1, ny
      t1 = f(j1,i,n)
      f(j1,i,n) = f(j,i,n)
      f(j,i,n) = t1
   10 continue
   20 continue
c first transform in x
      nrx = nxyz/nxh
      do 60 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 i = 1, ny
      do 30 j = 1, ns
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(j+k2,i,n)
      f(j+k2,i,n) = f(j+k1,i,n) - t2
      f(j+k1,i,n) = f(j+k1,i,n) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      do 80 k = 1, ny
      do 70 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      t2 = conjg(f(nxh2-j,k,n))
      t1 = f(j,k,n) + t2
      t2 = (f(j,k,n) - t2)*t3
      f(j,k,n) = ani*(t1 + t2)
      f(nxh2-j,k,n) = ani*conjg(t1 - t2)
   70 continue
   80 continue
      ani = 2.0*ani
      do 90 k = 1, ny
      f(nxhh+1,k,n) = ani*conjg(f(nxhh+1,k,n))
      f(1,k,n) = ani*cmplx(real(f(1,k,n)) + aimag(f(1,k,n)),
     1                     real(f(1,k,n)) - aimag(f(1,k,n)))
   90 continue
c bit-reverse array elements in y
      nry = nxhyz/ny
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 110
      do 100 i = 1, nxh
      t1 = f(i,k1,n)
      f(i,k1,n) = f(i,k,n)
      f(i,k,n) = t1
  100 continue
  110 continue
c then transform in y
      nry = nxyz/ny
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
      t1 = sct(1+kmr*(j-1))
      do 120 i = 1, nxh
      t2 = t1*f(i,j2,n)
      f(i,j2,n) = f(i,j1,n) - t2
      f(i,j1,n) = f(i,j1,n) + t2
  120 continue
  130 continue
  140 continue
  150 continue
c unscramble modes kx = 0, nx/2
      do 160 k = 2, nyh
      t1 = f(1,ny2-k,n)
      f(1,ny2-k,n) = 0.5*cmplx(aimag(f(1,k,n) + t1),real(f(1,k,n) - t1))
      f(1,k,n) = 0.5*cmplx(real(f(1,k,n) + t1),aimag(f(1,k,n) - t1))
  160 continue
  170 continue
      return
c forward fourier transform
  180 do 350 n = nzi, nzt
c scramble modes kx = 0, nx/2
      do 190 k = 2, nyh
      t1 = cmplx(aimag(f(1,ny2-k,n)),real(f(1,ny2-k,n)))
      f(1,ny2-k,n) = conjg(f(1,k,n) - t1)
      f(1,k,n) = f(1,k,n) + t1
  190 continue
c bit-reverse array elements in y
      nry = nxhyz/ny
      do 210 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 210
      do 200 i = 1, nxh
      t1 = f(i,k1,n)
      f(i,k1,n) = f(i,k,n)
      f(i,k,n) = t1
  200 continue
  210 continue
c then transform in y
      nry = nxyz/ny
      do 250 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 220 i = 1, nxh
      t2 = t1*f(i,j2,n)
      f(i,j2,n) = f(i,j1,n) - t2
      f(i,j1,n) = f(i,j1,n) + t2
  220 continue
  230 continue
  240 continue
  250 continue
c scramble coefficients
      kmr = nxyz/nx
      do 270 k = 1, ny
      do 260 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      t2 = conjg(f(nxh2-j,k,n))
      t1 = f(j,k,n) + t2
      t2 = (f(j,k,n) - t2)*t3
      f(j,k,n) = t1 + t2
      f(nxh2-j,k,n) = conjg(t1 - t2)
  260 continue
  270 continue
      do 280 k = 1, ny
      f(nxhh+1,k,n) = 2.0*conjg(f(nxhh+1,k,n))
      f(1,k,n) = cmplx(real(f(1,k,n)) + aimag(f(1,k,n)),
     1                 real(f(1,k,n)) - aimag(f(1,k,n)))
  280 continue
c bit-reverse array elements in x
      nrx = nxhyz/nxh
      do 300 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 300
      do 290 i = 1, ny
      t1 = f(j1,i,n)
      f(j1,i,n) = f(j,i,n)
      f(j,i,n) = t1
  290 continue
  300 continue
c finally transform in x
      nrx = nxyz/nxh
      do 340 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 330 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 320 i = 1, ny
      do 310 j = 1, ns
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(j+k2,i,n)
      f(j+k2,i,n) = f(j+k1,i,n) - t2
      f(j+k1,i,n) = f(j+k1,i,n) + t2
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RXZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,nxhd, 
     1nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the z part of a three dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c f(j,k,l) = sum(f(j,k,i)*exp(-sqrt(-1)*2pi*l*i/nz))
c if isign = 1, a forward fourier transform is performed
c f(n,m,i) = sum(f(n,m,l)*exp(sqrt(-1)*2pi*l*i/nz))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f
c nyd,nzd = second and third dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(j,k,l)= real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1,1,1)) = real part of mode nx/2,0,0
c imag(f(1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh
      integer nz, nzh, nz2, nxyz, nxhyz, nyt, nrz
      integer i, j, k, l, n, j1, j2, k1, k2, l1, ns, ns2, km, kmr
      complex t1, t2
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 100
c inverse fourier transform
c bit-reverse array elements in z
      nrz = nxhyz/nz
      do 30 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 30
      do 20 n = nyi, nyt
      do 10 i = 1, nxh
      t1 = f(i,n,l1)
      f(i,n,l1) = f(i,n,l)
      f(i,n,l) = t1
   10 continue
   20 continue
   30 continue
c finally transform in z
      nrz = nxyz/nz
      do 80 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 n = nyi, nyt
      do 40 i = 1, nxh
      t2 = t1*f(i,n,j2)
      f(i,n,j2) = f(i,n,j1) - t2
      f(i,n,j1) = f(i,n,j1) + t2
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 90 n = 2, nzh
      if (nyi.eq.1) then
         t1 = f(1,1,nz2-n)
         f(1,1,nz2-n) = 0.5*cmplx(aimag(f(1,1,n) + t1),
     1                            real(f(1,1,n) - t1))
         f(1,1,n) = 0.5*cmplx(real(f(1,1,n) + t1),aimag(f(1,1,n) - t1))
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         t1 = f(1,nyh+1,nz2-n)
         f(1,nyh+1,nz2-n) = 0.5*cmplx(aimag(f(1,nyh+1,n) + t1),
     1                                real(f(1,nyh+1,n) - t1))
         f(1,nyh+1,n) = 0.5*cmplx(real(f(1,nyh+1,n) + t1),
     1                            aimag(f(1,nyh+1,n) - t1))
      endif
   90 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  100 do 110 n = 2, nzh
      if (nyi.eq.1) then
         t1 = cmplx(aimag(f(1,1,nz2-n)),real(f(1,1,nz2-n)))
         f(1,1,nz2-n) = conjg(f(1,1,n) - t1)
         f(1,1,n) = f(1,1,n) + t1
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         t1 = cmplx(aimag(f(1,nyh+1,nz2-n)),real(f(1,nyh+1,nz2-n)))
         f(1,nyh+1,nz2-n) = conjg(f(1,nyh+1,n) - t1)
         f(1,nyh+1,n) = f(1,nyh+1,n) + t1
      endif
  110 continue
c bit-reverse array elements in z
      nrz = nxhyz/nz
      do 140 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 140
      do 130 n = nyi, nyt
      do 120 i = 1, nxh
      t1 = f(i,n,l1)
      f(i,n,l1) = f(i,n,l)
      f(i,n,l) = t1
  120 continue
  130 continue
  140 continue
c first transform in z
      nrz = nxyz/nz
      do 190 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 160 n = nyi, nyt
      do 150 i = 1, nxh
      t2 = t1*f(i,n,j2)
      f(i,n,j2) = f(i,n,j1) - t2
      f(i,n,j1) = f(i,n,j1) + t2
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RV3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,nxhd
     1,nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the x-y part of 3 three dimensional complex
c to real fast fourier transforms and their inverses, for a subset of z,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, three inverse fourier transforms is performed
c f(1:3,n,m,i) = (1/nx*ny*nz)*sum(f(1:3,j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
c       *exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, three forward fourier transforms are performed
c f(1:3,j,k,l) = sum(f(1:3,n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nzi = initial z index used
c nzp = number of z indices used
c nxhd = second dimension of f
c nyd,nzd = third and fourth dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(1:3,j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1:3,1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1:3,1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1:3,1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1:3,1,1,1)) = real part of mode nx/2,0,0
c imag(f(1:3,1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1:3,1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1:3,1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
      integer nxhyzd,nxyzhd
      complex f, sct
      integer mixup
      dimension f(4,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nxyz, nxhyz, nzt, nrx, nry
      integer i, j, k, l, n, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      real at1, at2, ani
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nzt = nzi + nzp - 1
      if (isign.gt.0) go to 230
c inverse fourier transform
      do 220 n = nzi, nzt
c swap complex components
      do 20 i = 1, ny
      do 10 j = 1, nxh
      at1 = aimag(f(3,j,i,n))
      at2 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(real(f(2,j,i,n)),real(f(4,j,i,n)))
      f(2,j,i,n) = cmplx(aimag(f(1,j,i,n)),at1)
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
   10 continue
   20 continue
c bit-reverse array elements in x
      nrx = nxhyz/nxh
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 i = 1, ny
      t1 = f(1,j1,i,n)
      t2 = f(2,j1,i,n)
      t3 = f(3,j1,i,n)
      f(1,j1,i,n) = f(1,j,i,n)
      f(2,j1,i,n) = f(2,j,i,n)
      f(3,j1,i,n) = f(3,j,i,n)
      f(1,j,i,n) = t1
      f(2,j,i,n) = t2
      f(3,j,i,n) = t3
   30 continue
   40 continue
c first transform in x
      nrx = nxyz/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 i = 1, ny
      do 50 j = 1, ns
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(1,j+k2,i,n)
      t3 = t1*f(2,j+k2,i,n)
      t4 = t1*f(3,j+k2,i,n)
      f(1,j+k2,i,n) = f(1,j+k1,i,n) - t2
      f(2,j+k2,i,n) = f(2,j+k1,i,n) - t3
      f(3,j+k2,i,n) = f(3,j+k1,i,n) - t4
      f(1,j+k1,i,n) = f(1,j+k1,i,n) + t2
      f(2,j+k1,i,n) = f(2,j+k1,i,n) + t3
      f(3,j+k1,i,n) = f(3,j+k1,i,n) + t4
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      do 110 k = 1, ny
      do 100 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 90 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k,n))
      t1 = f(jj,j,k,n) + t2
      t2 = (f(jj,j,k,n) - t2)*t3
      f(jj,j,k,n) = ani*(t1 + t2)
      f(jj,nxh2-j,k,n) = ani*conjg(t1 - t2)
   90 continue
  100 continue
  110 continue
      ani = 2.0*ani
      do 130 k = 1, ny
      do 120 jj = 1, 3
      f(jj,nxhh+1,k,n) = ani*conjg(f(jj,nxhh+1,k,n))
      f(jj,1,k,n) = ani*cmplx(real(f(jj,1,k,n)) + aimag(f(jj,1,k,n)),
     1                        real(f(jj,1,k,n)) - aimag(f(jj,1,k,n)))
  120 continue
  130 continue
c bit-reverse array elements in y
      nry = nxhyz/ny
      do 150 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 150
      do 140 i = 1, nxh
      t1 = f(1,i,k1,n)
      t2 = f(2,i,k1,n)
      t3 = f(3,i,k1,n)
      f(1,i,k1,n) = f(1,i,k,n)
      f(2,i,k1,n) = f(2,i,k,n)
      f(3,i,k1,n) = f(3,i,k,n)
      f(1,i,k,n) = t1
      f(2,i,k,n) = t2
      f(3,i,k,n) = t3
  140 continue
  150 continue
c then transform in y
      nry = nxyz/ny
      do 190 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 160 i = 1, nxh
      t2 = t1*f(1,i,j2,n)
      t3 = t1*f(2,i,j2,n)
      t4 = t1*f(3,i,j2,n)
      f(1,i,j2,n) = f(1,i,j1,n) - t2
      f(2,i,j2,n) = f(2,i,j1,n) - t3
      f(3,i,j2,n) = f(3,i,j1,n) - t4
      f(1,i,j1,n) = f(1,i,j1,n) + t2
      f(2,i,j1,n) = f(2,i,j1,n) + t3
      f(3,i,j1,n) = f(3,i,j1,n) + t4
  160 continue
  170 continue
  180 continue
  190 continue
c unscramble modes kx = 0, nx/2
      do 210 k = 2, nyh
      do 200 jj = 1, 3
      t1 = f(jj,1,ny2-k,n)
      f(jj,1,ny2-k,n) = 0.5*cmplx(aimag(f(jj,1,k,n) + t1),
     1                            real(f(jj,1,k,n) - t1))
      f(jj,1,k,n) = 0.5*cmplx(real(f(jj,1,k,n) + t1),
     1                        aimag(f(jj,1,k,n) - t1))
  200 continue
  210 continue
  220 continue
      return
c forward fourier transform
  230 do 450 n = nzi, nzt
c scramble modes kx = 0, nx/2
      do 250 k = 2, nyh
      do 240 jj = 1, 3
      t1 = cmplx(aimag(f(jj,1,ny2-k,n)),real(f(jj,1,ny2-k,n)))
      f(jj,1,ny2-k,n) = conjg(f(jj,1,k,n) - t1)
      f(jj,1,k,n) = f(jj,1,k,n) + t1
  240 continue
  250 continue
c bit-reverse array elements in y
      nry = nxhyz/ny
      do 270 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 270
      do 260 i = 1, nxh
      t1 = f(1,i,k1,n)
      t2 = f(2,i,k1,n)
      t3 = f(3,i,k1,n)
      f(1,i,k1,n) = f(1,i,k,n)
      f(2,i,k1,n) = f(2,i,k,n)
      f(3,i,k1,n) = f(3,i,k,n)
      f(1,i,k,n) = t1
      f(2,i,k,n) = t2
      f(3,i,k,n) = t3
  260 continue
  270 continue
c then transform in y
      nry = nxyz/ny
      do 310 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 300 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 290 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 280 i = 1, nxh
      t2 = t1*f(1,i,j2,n)
      t3 = t1*f(2,i,j2,n)
      t4 = t1*f(3,i,j2,n)
      f(1,i,j2,n) = f(1,i,j1,n) - t2
      f(2,i,j2,n) = f(2,i,j1,n) - t3
      f(3,i,j2,n) = f(3,i,j1,n) - t4
      f(1,i,j1,n) = f(1,i,j1,n) + t2
      f(2,i,j1,n) = f(2,i,j1,n) + t3
      f(3,i,j1,n) = f(3,i,j1,n) + t4
  280 continue
  290 continue
  300 continue
  310 continue
c scramble coefficients
      kmr = nxyz/nx
      do 340 k = 1, ny
      do 330 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 320 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k,n))
      t1 = f(jj,j,k,n) + t2
      t2 = (f(jj,j,k,n) - t2)*t3
      f(jj,j,k,n) = t1 + t2
      f(jj,nxh2-j,k,n) = conjg(t1 - t2)
  320 continue
  330 continue
  340 continue
      do 360 k = 1, ny
      do 350 jj = 1, 3
      f(jj,nxhh+1,k,n) = 2.0*conjg(f(jj,nxhh+1,k,n))
      f(jj,1,k,n) = cmplx(real(f(jj,1,k,n)) + aimag(f(jj,1,k,n)),
     1                    real(f(jj,1,k,n)) - aimag(f(jj,1,k,n)))
  350 continue
  360 continue
c bit-reverse array elements in x
      nrx = nxhyz/nxh
      do 380 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 380
      do 370 i = 1, ny
      t1 = f(1,j1,i,n)
      t2 = f(2,j1,i,n)
      t3 = f(3,j1,i,n)
      f(1,j1,i,n) = f(1,j,i,n)
      f(2,j1,i,n) = f(2,j,i,n)
      f(3,j1,i,n) = f(3,j,i,n)
      f(1,j,i,n) = t1
      f(2,j,i,n) = t2
      f(3,j,i,n) = t3
  370 continue
  380 continue
c finally transform in x
      nrx = nxyz/nxh
      do 420 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 410 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 400 i = 1, ny
      do 390 j = 1, ns
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(1,j+k2,i,n)
      t3 = t1*f(2,j+k2,i,n)
      t4 = t1*f(3,j+k2,i,n)
      f(1,j+k2,i,n) = f(1,j+k1,i,n) - t2
      f(2,j+k2,i,n) = f(2,j+k1,i,n) - t3
      f(3,j+k2,i,n) = f(3,j+k1,i,n) - t4
      f(1,j+k1,i,n) = f(1,j+k1,i,n) + t2
      f(2,j+k1,i,n) = f(2,j+k1,i,n) + t3
      f(3,j+k1,i,n) = f(3,j+k1,i,n) + t4
  390 continue
  400 continue
  410 continue
  420 continue
c swap complex components
      do 440 i = 1, ny
      do 430 j = 1, nxh
      f(4,j,i,n) = cmplx(aimag(f(3,j,i,n)),aimag(f(4,j,i,n)))
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(aimag(f(1,j,i,n)),aimag(f(2,j,i,n)))
      at2 = real(f(2,j,i,n))
      f(2,j,i,n) = cmplx(at1,0.0)
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
  430 continue
  440 continue
  450 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RV3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,nxhd,
     1nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the z part of 3 three dimensional complex to
c real fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, three inverse fourier transforms is performed
c f(1:3,j,k,l) = sum(f(1:3,j,k,i)*exp(-sqrt(-1)*2pi*l*i/nz))
c if isign = 1, three forward fourier transforms are performed
c f(1:3,n,m,i) = sum(f(1:3,n,m,l)*exp(sqrt(-1)*2pi*l*i/nz))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = second dimension of f
c nyd,nzd = third and fourth dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(1:3,j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1:3,1,k,l), = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1:3,1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1:3,1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1:3,1,1,1)) = real part of mode nx/2,0,0
c imag(f(1:3,1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1:3,1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1:3,1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(4,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh
      integer nz, nzh, nz2, nxyz, nxhyz, nyt, nrz
      integer i, j, k, l, n, jj, j1, j2, k1, k2, l1, ns, ns2, km, kmr
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 120
c inverse fourier transform
c bit-reverse array elements in z
      nrz = nxhyz/nz
      do 30 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 30
      do 20 n = nyi, nyt
      do 10 i = 1, nxh
      t1 = f(1,i,n,l1)
      t2 = f(2,i,n,l1)
      t3 = f(3,i,n,l1)
      f(1,i,n,l1) = f(1,i,n,l)
      f(2,i,n,l1) = f(2,i,n,l)
      f(3,i,n,l1) = f(3,i,n,l)
      f(1,i,n,l) = t1
      f(2,i,n,l) = t2
      f(3,i,n,l) = t3
   10 continue
   20 continue
   30 continue
c finally transform in z
      nrz = nxyz/nz
      do 80 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 n = nyi, nyt
      do 40 i = 1, nxh
      t2 = t1*f(1,i,n,j2)
      t3 = t1*f(2,i,n,j2)
      t4 = t1*f(3,i,n,j2)
      f(1,i,n,j2) = f(1,i,n,j1) - t2
      f(2,i,n,j2) = f(2,i,n,j1) - t3
      f(3,i,n,j2) = f(3,i,n,j1) - t4
      f(1,i,n,j1) = f(1,i,n,j1) + t2
      f(2,i,n,j1) = f(2,i,n,j1) + t3
      f(3,i,n,j1) = f(3,i,n,j1) + t4
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 110 n = 2, nzh
      if (nyi.eq.1) then
         do 90 jj = 1, 3
         t1 = f(jj,1,1,nz2-n)
         f(jj,1,1,nz2-n) = 0.5*cmplx(aimag(f(jj,1,1,n) + t1),
     1                               real(f(jj,1,1,n) - t1))
         f(jj,1,1,n) = 0.5*cmplx(real(f(jj,1,1,n) + t1),
     1                           aimag(f(jj,1,1,n) - t1))
   90    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         do 100 jj = 1, 3
         t1 = f(jj,1,nyh+1,nz2-n)
         f(jj,1,nyh+1,nz2-n) = 0.5*cmplx(aimag(f(jj,1,nyh+1,n) + t1),
     1                                  real(f(jj,1,nyh+1,n) - t1))
         f(jj,1,nyh+1,n) = 0.5*cmplx(real(f(jj,1,nyh+1,n) + t1),
     1                              aimag(f(jj,1,nyh+1,n) - t1))
  100    continue
      endif
  110 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  120 do 150 n = 2, nzh
      if (nyi.eq.1) then
         do 130 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,1,nz2-n)),real(f(jj,1,1,nz2-n)))
         f(jj,1,1,nz2-n) = conjg(f(jj,1,1,n) - t1)
         f(jj,1,1,n) = f(jj,1,1,n) + t1
  130    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         do 140 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,nyh+1,nz2-n)),
     1              real(f(jj,1,nyh+1,nz2-n)))
         f(jj,1,nyh+1,nz2-n) = conjg(f(jj,1,nyh+1,n) - t1)
         f(jj,1,nyh+1,n) = f(jj,1,nyh+1,n) + t1
  140    continue
      endif
  150 continue
c bit-reverse array elements in z
      nrz = nxhyz/nz
      do 180 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 180
      do 170 n = nyi, nyt
      do 160 i = 1, nxh
      t1 = f(1,i,n,l1)
      t2 = f(2,i,n,l1)
      t3 = f(3,i,n,l1)
      f(1,i,n,l1) = f(1,i,n,l)
      f(2,i,n,l1) = f(2,i,n,l)
      f(3,i,n,l1) = f(3,i,n,l)
      f(1,i,n,l) = t1
      f(2,i,n,l) = t2
      f(3,i,n,l) = t3
  160 continue
  170 continue
  180 continue
c first transform in z
      nrz = nxyz/nz
      do 230 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 220 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 210 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 200 n = nyi, nyt
      do 190 i = 1, nxh
      t2 = t1*f(1,i,n,j2)
      t3 = t1*f(2,i,n,j2)
      t4 = t1*f(3,i,n,j2)
      f(1,i,n,j2) = f(1,i,n,j1) - t2
      f(2,i,n,j2) = f(2,i,n,j1) - t3
      f(3,i,n,j2) = f(3,i,n,j1) - t4
      f(1,i,n,j1) = f(1,i,n,j1) + t2
      f(2,i,n,j1) = f(2,i,n,j1) + t3
      f(3,i,n,j1) = f(3,i,n,j1) + t4
  190 continue
  200 continue
  210 continue
  220 continue
  230 continue
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
