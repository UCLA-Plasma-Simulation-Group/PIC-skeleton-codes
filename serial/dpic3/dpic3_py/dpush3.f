c Fortran Library for Skeleton 3D Darwin PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR3(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,idimp, 
     1nop,nx,ny,nz,ipbc)
c for 3d code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c part(4,n) = velocity vx of particle n
c part(5,n) = velocity vy of particle n
c part(6,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny/nz = system length in x/y/z direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, npz, idimp, nop, nx, ny, nz, ipbc
      real vtx, vty, vtz, vdx, vdy, vdz
      real part
      dimension part(idimp,nop)
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
      part(1,j+k1) = edgelx + at1*(real(j) - .5)
      part(2,j+k1) = at4
      part(3,j+k1) = at5
   10 continue
   20 continue
   30 continue
c maxwellian velocity distribution
      do 40 j = 1, npxyz
      part(4,j) = vtx*ranorm()
      part(5,j) = vty*ranorm()
      part(6,j) = vtz*ranorm()
   40 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 50 j = 1, npxyz
      dsum1 = dsum1 + part(4,j)
      dsum2 = dsum2 + part(5,j)
      dsum3 = dsum3 + part(6,j)
   50 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1.0/real(npxyz)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 60 j = 1, npxyz
      part(4,j) = part(4,j) - sum1
      part(5,j) = part(5,j) - sum2
      part(6,j) = part(6,j) - sum3
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GBPUSH3L(part,fxyz,bxyz,qbm,dt,dtc,ek,idimp,nop,nx,ny, 
     1nz,nxv,nyv,nzv,ipbc)
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c part(4,n) = velocity vx of particle n
c part(5,n) = velocity vy of particle n
c part(6,n) = velocity vz of particle n
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
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nop, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qbm, dt, dtc, ek
      real part, fxyz, bxyz
      dimension part(idimp,nop)
      dimension fxyz(3,nxv,nyv,nzv), bxyz(3,nxv,nyv,nzv)
c local data
      integer j, nn, mm, ll, np, mp, lp
      real qtmh, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
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
      nn = part(1,j)
      mm = part(2,j)
      ll = part(3,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      dzp = part(3,j) - real(ll)
      nn = nn + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      np = nn + 1
      mm = mm + 1
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1.0 - dzp
      ll = ll + 1
      amy = dxp*amy
      lp = ll + 1
c find electric field
      dx = amz*(amx*fxyz(1,nn,mm,ll) + amy*fxyz(1,np,mm,ll)
     1        + dyp*fxyz(1,nn,mp,ll) + dx1*fxyz(1,np,mp,ll))
     2   + dzp*(amx*fxyz(1,nn,mm,lp) + amy*fxyz(1,np,mm,lp) 
     3        + dyp*fxyz(1,nn,mp,lp) + dx1*fxyz(1,np,mp,lp))
      dy = amz*(amx*fxyz(2,nn,mm,ll) + amy*fxyz(2,np,mm,ll)
     1        + dyp*fxyz(2,nn,mp,ll) + dx1*fxyz(2,np,mp,ll))
     2   + dzp*(amx*fxyz(2,nn,mm,lp) + amy*fxyz(2,np,mm,lp)
     3        + dyp*fxyz(2,nn,mp,lp) + dx1*fxyz(2,np,mp,lp))
      dz = amz*(amx*fxyz(3,nn,mm,ll) + amy*fxyz(3,np,mm,ll)
     1        + dyp*fxyz(3,nn,mp,ll) + dx1*fxyz(3,np,mp,ll))
     2   + dzp*(amx*fxyz(3,nn,mm,lp) + amy*fxyz(3,np,mm,lp)
     3        + dyp*fxyz(3,nn,mp,lp) + dx1*fxyz(3,np,mp,lp))
c find magnetic field
      ox = amz*(amx*bxyz(1,nn,mm,ll) + amy*bxyz(1,np,mm,ll)
     1        + dyp*bxyz(1,nn,mp,ll) + dx1*bxyz(1,np,mp,ll))
     2   + dzp*(amx*bxyz(1,nn,mm,lp) + amy*bxyz(1,np,mm,lp)
     3        + dyp*bxyz(1,nn,mp,lp) + dx1*bxyz(1,np,mp,lp))
      oy = amz*(amx*bxyz(2,nn,mm,ll) + amy*bxyz(2,np,mm,ll)
     1        + dyp*bxyz(2,nn,mp,ll) + dx1*bxyz(2,np,mp,ll))
     2   + dzp*(amx*bxyz(2,nn,mm,lp) + amy*bxyz(2,np,mm,lp)
     3        + dyp*bxyz(2,nn,mp,lp) + dx1*bxyz(2,np,mp,lp))
      oz = amz*(amx*bxyz(3,nn,mm,ll) + amy*bxyz(3,np,mm,ll)
     1        + dyp*bxyz(3,nn,mp,ll) + dx1*bxyz(3,np,mp,ll))
     2   + dzp*(amx*bxyz(3,nn,mm,lp) + amy*bxyz(3,np,mm,lp)
     3        + dyp*bxyz(3,nn,mp,lp) + dx1*bxyz(3,np,mp,lp))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j) + dx
      acy = part(5,j) + dy
      acz = part(6,j) + dz
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
      part(4,j) = dx
      part(5,j) = dy
      part(6,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
      dz = part(3,j) + dz*dtc
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
            dx = part(1,j)
            part(4,j) = -part(4,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(5,j) = -part(5,j)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j)
            part(6,j) = -part(6,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(4,j) = -part(4,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(5,j) = -part(5,j)
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
      part(3,j) = dz
   10 continue
c normalize kinetic energy
      ek = ek + 0.5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine GPOST3L(part,q,qm,nop,idimp,nxv,nyv,nzv)
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c q(j,k,l) = charge density at grid point j,k,l
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 6
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
c nzv = third dimension of charge array, must be >= nz+1
      implicit none
      integer nop, idimp, nxv, nyv, nzv
      real qm
      real part, q
      dimension part(idimp,nop), q(nxv,nyv,nzv)
c local data
      integer j, nn, mm, ll, np, mp, lp
      real dx1, dxp, dyp, dzp, amx, amy, amz
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      ll = part(3,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      dzp = part(3,j) - real(ll)
      nn = nn + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      np = nn + 1
      mm = mm + 1
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1.0 - dzp
      ll = ll + 1
      amy = dxp*amy
      lp = ll + 1
c deposit charge
      q(nn,mm,ll) = q(nn,mm,ll) + amx*amz
      q(np,mm,ll) = q(np,mm,ll) + amy*amz
      q(nn,mp,ll) = q(nn,mp,ll) + dyp*amz
      q(np,mp,ll) = q(np,mp,ll) + dx1*amz
      q(nn,mm,lp) = q(nn,mm,lp) + amx*dzp
      q(np,mm,lp) = q(np,mm,lp) + amy*dzp
      q(nn,mp,lp) = q(nn,mp,lp) + dyp*dzp
      q(np,mp,lp) = q(np,mp,lp) + dx1*dzp
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GJPOST3L(part,cu,qm,dt,nop,idimp,nx,ny,nz,nxv,nyv,nzv, 
     1ipbc)
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c part(4,n) = x velocity of particle n
c part(5,n) = y velocity of particle n
c part(6,n) = z velocity of particle n
c cu(i,j,k,l) = ith component of current density at grid point j,k,l
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 6
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c nzv = fourth dimension of current array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer nop, idimp, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qm, dt
      real part, cu
      dimension part(idimp,nop), cu(3,nxv,nyv,nzv)
c local data
      integer j, nn, mm, ll, np, mp, lp
      real edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz
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
      nn = part(1,j)
      mm = part(2,j)
      ll = part(3,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      dzp = part(3,j) - real(ll)
      nn = nn + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      np = nn + 1
      mm = mm + 1
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1.0 - dzp
      ll = ll + 1
      amy = dxp*amy
      lp = ll + 1
c deposit current
      dx = amx*amz
      dy = amy*amz
      vx = part(4,j)
      vy = part(5,j)
      vz = part(6,j)
      cu(1,nn,mm,ll) = cu(1,nn,mm,ll) + vx*dx
      cu(2,nn,mm,ll) = cu(2,nn,mm,ll) + vy*dx
      cu(3,nn,mm,ll) = cu(3,nn,mm,ll) + vz*dx
      dx = dyp*amz
      cu(1,np,mm,ll) = cu(1,np,mm,ll) + vx*dy
      cu(2,np,mm,ll) = cu(2,np,mm,ll) + vy*dy
      cu(3,np,mm,ll) = cu(3,np,mm,ll) + vz*dy
      dy = dx1*amz
      cu(1,nn,mp,ll) = cu(1,nn,mp,ll) + vx*dx
      cu(2,nn,mp,ll) = cu(2,nn,mp,ll) + vy*dx
      cu(3,nn,mp,ll) = cu(3,nn,mp,ll) + vz*dx
      dx = amx*dzp
      cu(1,np,mp,ll) = cu(1,np,mp,ll) + vx*dy
      cu(2,np,mp,ll) = cu(2,np,mp,ll) + vy*dy
      cu(3,np,mp,ll) = cu(3,np,mp,ll) + vz*dy
      dy = amy*dzp
      cu(1,nn,mm,lp) = cu(1,nn,mm,lp) + vx*dx
      cu(2,nn,mm,lp) = cu(2,nn,mm,lp) + vy*dx
      cu(3,nn,mm,lp) = cu(3,nn,mm,lp) + vz*dx
      dx = dyp*dzp
      cu(1,np,mm,lp) = cu(1,np,mm,lp) + vx*dy
      cu(2,np,mm,lp) = cu(2,np,mm,lp) + vy*dy
      cu(3,np,mm,lp) = cu(3,np,mm,lp) + vz*dy
      dy = dx1*dzp
      cu(1,nn,mp,lp) = cu(1,nn,mp,lp) + vx*dx
      cu(2,nn,mp,lp) = cu(2,nn,mp,lp) + vy*dx
      cu(3,nn,mp,lp) = cu(3,nn,mp,lp) + vz*dx
      cu(1,np,mp,lp) = cu(1,np,mp,lp) + vx*dy
      cu(2,np,mp,lp) = cu(2,np,mp,lp) + vy*dy
      cu(3,np,mp,lp) = cu(3,np,mp,lp) + vz*dy
c advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
      dz = part(3,j) + vz*dt
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
            dx = part(1,j)
            part(4,j) = -part(4,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(5,j) = -part(5,j)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j)
            part(6,j) = -part(6,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(4,j) = -part(4,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(5,j) = -part(5,j)
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
      part(3,j) = dz
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GMJPOST3L(part,amu,qm,nop,idimp,nxv,nyv,nzv)
c for 3d code, this subroutine calculates particle momentum flux
c using first-order linear interpolation
c scalar version using guard cells
c 121 flops/particle, 52 loads, 48 stores
c input: all, output: part, amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vj*vk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = position z of particle n at t
c part(4,n) = x velocity of particle n at t - dt/2
c part(5,n) = y velocity of particle n at t - dt/2
c part(6,n) = z velocity of particle n at t - dt/2
c amu(i,j,k,l) = ith component of momentum flux at grid point j,k,l
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 6
c nxv = second dimension of flux array, must be >= nx+1
c nyv = third dimension of flux array, must be >= ny+1
c nzv = fourth dimension of flux array, must be >= nz+1
      implicit none
      integer nop, idimp, nxv, nyv, nzv
      real qm
      real part, amu
      dimension part(idimp,nop), amu(6,nxv,nyv,nzv)
c local data
      integer j, nn, mm, ll, np, mp, lp
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, vx, vy, vz
      real v1, v2, v3, v4, v5, v6
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j)
      mm = part(2,j)
      ll = part(3,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      dzp = part(3,j) - real(ll)
      nn = nn + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      np = nn + 1
      mm = mm + 1
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1.0 - dzp
      ll = ll + 1
      amy = dxp*amy
      lp = ll + 1
c deposit current
      dx = amx*amz
      dy = amy*amz
      vx = part(4,j)
      vy = part(5,j)
      vz = part(6,j)
      v1 = vx*vx
      v2 = vx*vy
      v3 = vx*vz
      v4 = vy*vy
      v5 = vy*vz
      v6 = vz*vz
      amu(1,nn,mm,ll) = amu(1,nn,mm,ll) + v1*dx
      amu(2,nn,mm,ll) = amu(2,nn,mm,ll) + v2*dx
      amu(3,nn,mm,ll) = amu(3,nn,mm,ll) + v3*dx
      amu(4,nn,mm,ll) = amu(4,nn,mm,ll) + v4*dx
      amu(5,nn,mm,ll) = amu(5,nn,mm,ll) + v5*dx
      amu(6,nn,mm,ll) = amu(6,nn,mm,ll) + v6*dx
      dx = dyp*amz
      amu(1,np,mm,ll) = amu(1,np,mm,ll) + v1*dy
      amu(2,np,mm,ll) = amu(2,np,mm,ll) + v2*dy
      amu(3,np,mm,ll) = amu(3,np,mm,ll) + v3*dy
      amu(4,np,mm,ll) = amu(4,np,mm,ll) + v4*dy
      amu(5,np,mm,ll) = amu(5,np,mm,ll) + v5*dy
      amu(6,np,mm,ll) = amu(6,np,mm,ll) + v6*dy
      dy = dx1*amz
      amu(1,nn,mp,ll) = amu(1,nn,mp,ll) + v1*dx
      amu(2,nn,mp,ll) = amu(2,nn,mp,ll) + v2*dx
      amu(3,nn,mp,ll) = amu(3,nn,mp,ll) + v3*dx
      amu(4,nn,mp,ll) = amu(4,nn,mp,ll) + v4*dx
      amu(5,nn,mp,ll) = amu(5,nn,mp,ll) + v5*dx
      amu(6,nn,mp,ll) = amu(6,nn,mp,ll) + v6*dx
      dx = amx*dzp
      amu(1,np,mp,ll) = amu(1,np,mp,ll) + v1*dy
      amu(2,np,mp,ll) = amu(2,np,mp,ll) + v2*dy
      amu(3,np,mp,ll) = amu(3,np,mp,ll) + v3*dy
      amu(4,np,mp,ll) = amu(4,np,mp,ll) + v4*dy
      amu(5,np,mp,ll) = amu(5,np,mp,ll) + v5*dy
      amu(6,np,mp,ll) = amu(6,np,mp,ll) + v6*dy
      dy = amy*dzp
      amu(1,nn,mm,lp) = amu(1,nn,mm,lp) + v1*dx
      amu(2,nn,mm,lp) = amu(2,nn,mm,lp) + v2*dx
      amu(3,nn,mm,lp) = amu(3,nn,mm,lp) + v3*dx
      amu(4,nn,mm,lp) = amu(4,nn,mm,lp) + v4*dx
      amu(5,nn,mm,lp) = amu(5,nn,mm,lp) + v5*dx
      amu(6,nn,mm,lp) = amu(6,nn,mm,lp) + v6*dx
      dx = dyp*dzp
      amu(1,np,mm,lp) = amu(1,np,mm,lp) + v1*dy
      amu(2,np,mm,lp) = amu(2,np,mm,lp) + v2*dy
      amu(3,np,mm,lp) = amu(3,np,mm,lp) + v3*dy
      amu(4,np,mm,lp) = amu(4,np,mm,lp) + v4*dy
      amu(5,np,mm,lp) = amu(5,np,mm,lp) + v5*dy
      amu(6,np,mm,lp) = amu(6,np,mm,lp) + v6*dy
      dy = dx1*dzp
      amu(1,nn,mp,lp) = amu(1,nn,mp,lp) + v1*dx
      amu(2,nn,mp,lp) = amu(2,nn,mp,lp) + v2*dx
      amu(3,nn,mp,lp) = amu(3,nn,mp,lp) + v3*dx
      amu(4,nn,mp,lp) = amu(4,nn,mp,lp) + v4*dx
      amu(5,nn,mp,lp) = amu(5,nn,mp,lp) + v5*dx
      amu(6,nn,mp,lp) = amu(6,nn,mp,lp) + v6*dx
      amu(1,np,mp,lp) = amu(1,np,mp,lp) + v1*dy
      amu(2,np,mp,lp) = amu(2,np,mp,lp) + v2*dy
      amu(3,np,mp,lp) = amu(3,np,mp,lp) + v3*dy
      amu(4,np,mp,lp) = amu(4,np,mp,lp) + v4*dy
      amu(5,np,mp,lp) = amu(5,np,mp,lp) + v5*dy
      amu(6,np,mp,lp) = amu(6,np,mp,lp) + v6*dy
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GDJPOST3L(part,fxyz,bxyz,dcu,amu,qm,qbm,dt,idimp,nop,  
     1nxv,nyv,nzv)
c for 3d code, this subroutine calculates particle momentum flux
c and acceleration density using first-order spline interpolation.
c scalar version using guard cells
c 350 flops/particle, 1 divide, 126 loads, 72 stores
c input: all, output: dcu, amu
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c dcu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c dcu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c dcu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c dcu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c dcu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c dcu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c dcu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c and qci = qm*vj*vk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c velocity equations at t=t+dt/2 are calculated from:
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
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = position z of particle n at t
c part(4,n) = velocity vx of particle n at t - dt/2
c part(5,n) = velocity vy of particle n at t - dt/2
c part(6,n) = velocity vz of particle n at t - dt/2
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c bxyz(1,j,k,l) = x component of magnetic field at grid (j,k,l)
c bxyz(2,j,k,l) = y component of magnetic field at grid (j,k,l)
c bxyz(3,j,k,l) = z component of magnetic field at grid (j,k,l)
c that is, the convolution of magnetic field over particle shape
c dcu(i,j,k,l) = ith component of acceleration density
c at grid point j,k for i = 1, 3
c amu(i,j,k,l) = ith component of momentum flux
c at grid point j,k,l for i = 1, 6
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive force calculations
c idimp = size of phase space = 6
c nop = number of particles
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
      implicit none
      integer idimp, nop, nxv, nyv, nzv
      real qm, qbm, dt
      real part, fxyz, bxyz, dcu, amu
      dimension part(idimp,nop)
      dimension fxyz(3,nxv,nyv,nzv), bxyz(3,nxv,nyv,nzv)
      dimension dcu(3,nxv,nyv,nzv), amu(6,nxv,nyv,nzv)
c local data
      integer j, nn, mm, ll, np, mp, lp
      real qtmh, dti, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz
      real ox, oy, oz, dx1, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, v3, v4, v5, v6
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      ll = part(3,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      dzp = part(3,j) - real(ll)
      nn = nn + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      np = nn + 1
      mm = mm + 1
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1.0 - dzp
      ll = ll + 1
      amy = dxp*amy
      lp = ll + 1
c find electric field
      dx = amz*(amx*fxyz(1,nn,mm,ll) + amy*fxyz(1,np,mm,ll)
     1        + dyp*fxyz(1,nn,mp,ll) + dx1*fxyz(1,np,mp,ll))
     2   + dzp*(amx*fxyz(1,nn,mm,lp) + amy*fxyz(1,np,mm,lp) 
     3        + dyp*fxyz(1,nn,mp,lp) + dx1*fxyz(1,np,mp,lp))
      dy = amz*(amx*fxyz(2,nn,mm,ll) + amy*fxyz(2,np,mm,ll)
     1        + dyp*fxyz(2,nn,mp,ll) + dx1*fxyz(2,np,mp,ll))
     2   + dzp*(amx*fxyz(2,nn,mm,lp) + amy*fxyz(2,np,mm,lp)
     3        + dyp*fxyz(2,nn,mp,lp) + dx1*fxyz(2,np,mp,lp))
      dz = amz*(amx*fxyz(3,nn,mm,ll) + amy*fxyz(3,np,mm,ll)
     1        + dyp*fxyz(3,nn,mp,ll) + dx1*fxyz(3,np,mp,ll))
     2   + dzp*(amx*fxyz(3,nn,mm,lp) + amy*fxyz(3,np,mm,lp)
     3        + dyp*fxyz(3,nn,mp,lp) + dx1*fxyz(3,np,mp,lp))
c find magnetic field
      ox = amz*(amx*bxyz(1,nn,mm,ll) + amy*bxyz(1,np,mm,ll)
     1        + dyp*bxyz(1,nn,mp,ll) + dx1*bxyz(1,np,mp,ll))
     2   + dzp*(amx*bxyz(1,nn,mm,lp) + amy*bxyz(1,np,mm,lp)
     3        + dyp*bxyz(1,nn,mp,lp) + dx1*bxyz(1,np,mp,lp))
      oy = amz*(amx*bxyz(2,nn,mm,ll) + amy*bxyz(2,np,mm,ll)
     1        + dyp*bxyz(2,nn,mp,ll) + dx1*bxyz(2,np,mp,ll))
     2   + dzp*(amx*bxyz(2,nn,mm,lp) + amy*bxyz(2,np,mm,lp)
     3        + dyp*bxyz(2,nn,mp,lp) + dx1*bxyz(2,np,mp,lp))
      oz = amz*(amx*bxyz(3,nn,mm,ll) + amy*bxyz(3,np,mm,ll)
     1        + dyp*bxyz(3,nn,mp,ll) + dx1*bxyz(3,np,mp,ll))
     2   + dzp*(amx*bxyz(3,nn,mm,lp) + amy*bxyz(3,np,mm,lp)
     3        + dyp*bxyz(3,nn,mp,lp) + dx1*bxyz(3,np,mp,lp))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(4,j)
      vy = part(5,j)
      vz = part(6,j)
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
c deposit momentum flux and acceleration density,
      amz = qm*amz
      dzp = qm*dzp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = amx*amz
      dy = amy*amz
      v1 = ox*ox
      v2 = ox*oy
      v3 = ox*oz
      v4 = oy*oy
      v5 = oy*oz
      v6 = oz*oz
      amu(1,nn,mm,ll) = amu(1,nn,mm,ll) + v1*dx
      amu(2,nn,mm,ll) = amu(2,nn,mm,ll) + v2*dx
      amu(3,nn,mm,ll) = amu(3,nn,mm,ll) + v3*dx
      amu(4,nn,mm,ll) = amu(4,nn,mm,ll) + v4*dx
      amu(5,nn,mm,ll) = amu(5,nn,mm,ll) + v5*dx
      amu(6,nn,mm,ll) = amu(6,nn,mm,ll) + v6*dx
      dcu(1,nn,mm,ll) = dcu(1,nn,mm,ll) + vx*dx
      dcu(2,nn,mm,ll) = dcu(2,nn,mm,ll) + vy*dx
      dcu(3,nn,mm,ll) = dcu(3,nn,mm,ll) + vz*dx
      dx = dyp*amz
      amu(1,np,mm,ll) = amu(1,np,mm,ll) + v1*dy
      amu(2,np,mm,ll) = amu(2,np,mm,ll) + v2*dy
      amu(3,np,mm,ll) = amu(3,np,mm,ll) + v3*dy
      amu(4,np,mm,ll) = amu(4,np,mm,ll) + v4*dy
      amu(5,np,mm,ll) = amu(5,np,mm,ll) + v5*dy
      amu(6,np,mm,ll) = amu(6,np,mm,ll) + v6*dy
      dcu(1,np,mm,ll) = dcu(1,np,mm,ll) + vx*dy
      dcu(2,np,mm,ll) = dcu(2,np,mm,ll) + vy*dy
      dcu(3,np,mm,ll) = dcu(3,np,mm,ll) + vz*dy
      dy = dx1*amz
      amu(1,nn,mp,ll) = amu(1,nn,mp,ll) + v1*dx
      amu(2,nn,mp,ll) = amu(2,nn,mp,ll) + v2*dx
      amu(3,nn,mp,ll) = amu(3,nn,mp,ll) + v3*dx
      amu(4,nn,mp,ll) = amu(4,nn,mp,ll) + v4*dx
      amu(5,nn,mp,ll) = amu(5,nn,mp,ll) + v5*dx
      amu(6,nn,mp,ll) = amu(6,nn,mp,ll) + v6*dx
      dcu(1,nn,mp,ll) = dcu(1,nn,mp,ll) + vx*dx
      dcu(2,nn,mp,ll) = dcu(2,nn,mp,ll) + vy*dx
      dcu(3,nn,mp,ll) = dcu(3,nn,mp,ll) + vz*dx
      dx = amx*dzp
      amu(1,np,mp,ll) = amu(1,np,mp,ll) + v1*dy
      amu(2,np,mp,ll) = amu(2,np,mp,ll) + v2*dy
      amu(3,np,mp,ll) = amu(3,np,mp,ll) + v3*dy
      amu(4,np,mp,ll) = amu(4,np,mp,ll) + v4*dy
      amu(5,np,mp,ll) = amu(5,np,mp,ll) + v5*dy
      amu(6,np,mp,ll) = amu(6,np,mp,ll) + v6*dy
      dcu(1,np,mp,ll) = dcu(1,np,mp,ll) + vx*dy
      dcu(2,np,mp,ll) = dcu(2,np,mp,ll) + vy*dy
      dcu(3,np,mp,ll) = dcu(3,np,mp,ll) + vz*dy
      dy = amy*dzp
      amu(1,nn,mm,lp) = amu(1,nn,mm,lp) + v1*dx
      amu(2,nn,mm,lp) = amu(2,nn,mm,lp) + v2*dx
      amu(3,nn,mm,lp) = amu(3,nn,mm,lp) + v3*dx
      amu(4,nn,mm,lp) = amu(4,nn,mm,lp) + v4*dx
      amu(5,nn,mm,lp) = amu(5,nn,mm,lp) + v5*dx
      amu(6,nn,mm,lp) = amu(6,nn,mm,lp) + v6*dx
      dcu(1,nn,mm,lp) = dcu(1,nn,mm,lp) + vx*dx
      dcu(2,nn,mm,lp) = dcu(2,nn,mm,lp) + vy*dx
      dcu(3,nn,mm,lp) = dcu(3,nn,mm,lp) + vz*dx
      dx = dyp*dzp
      amu(1,np,mm,lp) = amu(1,np,mm,lp) + v1*dy
      amu(2,np,mm,lp) = amu(2,np,mm,lp) + v2*dy
      amu(3,np,mm,lp) = amu(3,np,mm,lp) + v3*dy
      amu(4,np,mm,lp) = amu(4,np,mm,lp) + v4*dy
      amu(5,np,mm,lp) = amu(5,np,mm,lp) + v5*dy
      amu(6,np,mm,lp) = amu(6,np,mm,lp) + v6*dy
      dcu(1,np,mm,lp) = dcu(1,np,mm,lp) + vx*dy
      dcu(2,np,mm,lp) = dcu(2,np,mm,lp) + vy*dy
      dcu(3,np,mm,lp) = dcu(3,np,mm,lp) + vz*dy
      dy = dx1*dzp
      amu(1,nn,mp,lp) = amu(1,nn,mp,lp) + v1*dx
      amu(2,nn,mp,lp) = amu(2,nn,mp,lp) + v2*dx
      amu(3,nn,mp,lp) = amu(3,nn,mp,lp) + v3*dx
      amu(4,nn,mp,lp) = amu(4,nn,mp,lp) + v4*dx
      amu(5,nn,mp,lp) = amu(5,nn,mp,lp) + v5*dx
      amu(6,nn,mp,lp) = amu(6,nn,mp,lp) + v6*dx
      dcu(1,nn,mp,lp) = dcu(1,nn,mp,lp) + vx*dx
      dcu(2,nn,mp,lp) = dcu(2,nn,mp,lp) + vy*dx
      dcu(3,nn,mp,lp) = dcu(3,nn,mp,lp) + vz*dx
      amu(1,np,mp,lp) = amu(1,np,mp,lp) + v1*dy
      amu(2,np,mp,lp) = amu(2,np,mp,lp) + v2*dy
      amu(3,np,mp,lp) = amu(3,np,mp,lp) + v3*dy
      amu(4,np,mp,lp) = amu(4,np,mp,lp) + v4*dy
      amu(5,np,mp,lp) = amu(5,np,mp,lp) + v5*dy
      amu(6,np,mp,lp) = amu(6,np,mp,lp) + v6*dy
      dcu(1,np,mp,lp) = dcu(1,np,mp,lp) + vx*dy
      dcu(2,np,mp,lp) = dcu(2,np,mp,lp) + vy*dy
      dcu(3,np,mp,lp) = dcu(3,np,mp,lp) + vz*dy
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GDCJPOST3L(part,fxyz,bxyz,cu,dcu,amu,qm,qbm,dt,idimp,  
     1nop,nxv,nyv,nzv)
c for 3d code, this subroutine calculates particle momentum flux,
c acceleration density and current density using first-order spline
c interpolation.
c scalar version using guard cells
c 398 flops/particle, 1 divide, 150 loads, 96 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c dcu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c dcu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c dcu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c dcu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c dcu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c dcu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c dcu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c and qci = qm*vj*vk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c velocity equations at t=t+dt/2 are calculated from:
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
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = position z of particle n at t
c part(4,n) = velocity vx of particle n at t - dt/2
c part(5,n) = velocity vy of particle n at t - dt/2
c part(6,n) = velocity vz of particle n at t - dt/2
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c bxyz(1,j,k,l) = x component of magnetic field at grid (j,k,l)
c bxyz(2,j,k,l) = y component of magnetic field at grid (j,k,l)
c bxyz(3,j,k,l) = z component of magnetic field at grid (j,k,l)
c that is, the convolution of magnetic field over particle shape
c cu(i,j,k,l) = ith component of current density
c at grid point j,k for i = 1, 3
c dcu(i,j,k,l) = ith component of acceleration density
c at grid point j,k for i = 1, 3
c amu(i,j,k,l) = ith component of momentum flux
c at grid point j,k,l for i = 1, 6
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive force calculations
c idimp = size of phase space = 6
c nop = number of particles
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
      implicit none
      integer idimp, nop, nxv, nyv, nzv
      real qm, qbm, dt
      real part, fxyz, bxyz, cu, dcu, amu
      dimension part(idimp,nop)
      dimension fxyz(3,nxv,nyv,nzv), bxyz(3,nxv,nyv,nzv)
      dimension cu(3,nxv,nyv,nzv), dcu(3,nxv,nyv,nzv)
      dimension amu(6,nxv,nyv,nzv)
c local data
      integer j, nn, mm, ll, np, mp, lp
      real qtmh, dti, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz
      real ox, oy, oz, dx1, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, v3, v4, v5, v6
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      ll = part(3,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      dzp = part(3,j) - real(ll)
      nn = nn + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      np = nn + 1
      mm = mm + 1
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1.0 - dzp
      ll = ll + 1
      amy = dxp*amy
      lp = ll + 1
c find electric field
      dx = amz*(amx*fxyz(1,nn,mm,ll) + amy*fxyz(1,np,mm,ll)
     1        + dyp*fxyz(1,nn,mp,ll) + dx1*fxyz(1,np,mp,ll))
     2   + dzp*(amx*fxyz(1,nn,mm,lp) + amy*fxyz(1,np,mm,lp) 
     3        + dyp*fxyz(1,nn,mp,lp) + dx1*fxyz(1,np,mp,lp))
      dy = amz*(amx*fxyz(2,nn,mm,ll) + amy*fxyz(2,np,mm,ll)
     1        + dyp*fxyz(2,nn,mp,ll) + dx1*fxyz(2,np,mp,ll))
     2   + dzp*(amx*fxyz(2,nn,mm,lp) + amy*fxyz(2,np,mm,lp)
     3        + dyp*fxyz(2,nn,mp,lp) + dx1*fxyz(2,np,mp,lp))
      dz = amz*(amx*fxyz(3,nn,mm,ll) + amy*fxyz(3,np,mm,ll)
     1        + dyp*fxyz(3,nn,mp,ll) + dx1*fxyz(3,np,mp,ll))
     2   + dzp*(amx*fxyz(3,nn,mm,lp) + amy*fxyz(3,np,mm,lp)
     3        + dyp*fxyz(3,nn,mp,lp) + dx1*fxyz(3,np,mp,lp))
c find magnetic field
      ox = amz*(amx*bxyz(1,nn,mm,ll) + amy*bxyz(1,np,mm,ll)
     1        + dyp*bxyz(1,nn,mp,ll) + dx1*bxyz(1,np,mp,ll))
     2   + dzp*(amx*bxyz(1,nn,mm,lp) + amy*bxyz(1,np,mm,lp)
     3        + dyp*bxyz(1,nn,mp,lp) + dx1*bxyz(1,np,mp,lp))
      oy = amz*(amx*bxyz(2,nn,mm,ll) + amy*bxyz(2,np,mm,ll)
     1        + dyp*bxyz(2,nn,mp,ll) + dx1*bxyz(2,np,mp,ll))
     2   + dzp*(amx*bxyz(2,nn,mm,lp) + amy*bxyz(2,np,mm,lp)
     3        + dyp*bxyz(2,nn,mp,lp) + dx1*bxyz(2,np,mp,lp))
      oz = amz*(amx*bxyz(3,nn,mm,ll) + amy*bxyz(3,np,mm,ll)
     1        + dyp*bxyz(3,nn,mp,ll) + dx1*bxyz(3,np,mp,ll))
     2   + dzp*(amx*bxyz(3,nn,mm,lp) + amy*bxyz(3,np,mm,lp)
     3        + dyp*bxyz(3,nn,mp,lp) + dx1*bxyz(3,np,mp,lp))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(4,j)
      vy = part(5,j)
      vz = part(6,j)
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
      amz = qm*amz
      dzp = qm*dzp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = amx*amz
      dy = amy*amz
      v1 = ox*ox
      v2 = ox*oy
      v3 = ox*oz
      v4 = oy*oy
      v5 = oy*oz
      v6 = oz*oz
      amu(1,nn,mm,ll) = amu(1,nn,mm,ll) + v1*dx
      amu(2,nn,mm,ll) = amu(2,nn,mm,ll) + v2*dx
      amu(3,nn,mm,ll) = amu(3,nn,mm,ll) + v3*dx
      amu(4,nn,mm,ll) = amu(4,nn,mm,ll) + v4*dx
      amu(5,nn,mm,ll) = amu(5,nn,mm,ll) + v5*dx
      amu(6,nn,mm,ll) = amu(6,nn,mm,ll) + v6*dx
      dcu(1,nn,mm,ll) = dcu(1,nn,mm,ll) + vx*dx
      dcu(2,nn,mm,ll) = dcu(2,nn,mm,ll) + vy*dx
      dcu(3,nn,mm,ll) = dcu(3,nn,mm,ll) + vz*dx
      cu(1,nn,mm,ll) = cu(1,nn,mm,ll) + ox*dx
      cu(2,nn,mm,ll) = cu(2,nn,mm,ll) + oy*dx
      cu(3,nn,mm,ll) = cu(3,nn,mm,ll) + oz*dx
      dx = dyp*amz
      amu(1,np,mm,ll) = amu(1,np,mm,ll) + v1*dy
      amu(2,np,mm,ll) = amu(2,np,mm,ll) + v2*dy
      amu(3,np,mm,ll) = amu(3,np,mm,ll) + v3*dy
      amu(4,np,mm,ll) = amu(4,np,mm,ll) + v4*dy
      amu(5,np,mm,ll) = amu(5,np,mm,ll) + v5*dy
      amu(6,np,mm,ll) = amu(6,np,mm,ll) + v6*dy
      dcu(1,np,mm,ll) = dcu(1,np,mm,ll) + vx*dy
      dcu(2,np,mm,ll) = dcu(2,np,mm,ll) + vy*dy
      dcu(3,np,mm,ll) = dcu(3,np,mm,ll) + vz*dy
      cu(1,np,mm,ll) = cu(1,np,mm,ll) + ox*dy
      cu(2,np,mm,ll) = cu(2,np,mm,ll) + oy*dy
      cu(3,np,mm,ll) = cu(3,np,mm,ll) + oz*dy
      dy = dx1*amz
      amu(1,nn,mp,ll) = amu(1,nn,mp,ll) + v1*dx
      amu(2,nn,mp,ll) = amu(2,nn,mp,ll) + v2*dx
      amu(3,nn,mp,ll) = amu(3,nn,mp,ll) + v3*dx
      amu(4,nn,mp,ll) = amu(4,nn,mp,ll) + v4*dx
      amu(5,nn,mp,ll) = amu(5,nn,mp,ll) + v5*dx
      amu(6,nn,mp,ll) = amu(6,nn,mp,ll) + v6*dx
      dcu(1,nn,mp,ll) = dcu(1,nn,mp,ll) + vx*dx
      dcu(2,nn,mp,ll) = dcu(2,nn,mp,ll) + vy*dx
      dcu(3,nn,mp,ll) = dcu(3,nn,mp,ll) + vz*dx
      cu(1,nn,mp,ll) = cu(1,nn,mp,ll) + ox*dx
      cu(2,nn,mp,ll) = cu(2,nn,mp,ll) + oy*dx
      cu(3,nn,mp,ll) = cu(3,nn,mp,ll) + oz*dx
      dx = amx*dzp
      amu(1,np,mp,ll) = amu(1,np,mp,ll) + v1*dy
      amu(2,np,mp,ll) = amu(2,np,mp,ll) + v2*dy
      amu(3,np,mp,ll) = amu(3,np,mp,ll) + v3*dy
      amu(4,np,mp,ll) = amu(4,np,mp,ll) + v4*dy
      amu(5,np,mp,ll) = amu(5,np,mp,ll) + v5*dy
      amu(6,np,mp,ll) = amu(6,np,mp,ll) + v6*dy
      dcu(1,np,mp,ll) = dcu(1,np,mp,ll) + vx*dy
      dcu(2,np,mp,ll) = dcu(2,np,mp,ll) + vy*dy
      dcu(3,np,mp,ll) = dcu(3,np,mp,ll) + vz*dy
      cu(1,np,mp,ll) = cu(1,np,mp,ll) + ox*dy
      cu(2,np,mp,ll) = cu(2,np,mp,ll) + oy*dy
      cu(3,np,mp,ll) = cu(3,np,mp,ll) + oz*dy
      dy = amy*dzp
      amu(1,nn,mm,lp) = amu(1,nn,mm,lp) + v1*dx
      amu(2,nn,mm,lp) = amu(2,nn,mm,lp) + v2*dx
      amu(3,nn,mm,lp) = amu(3,nn,mm,lp) + v3*dx
      amu(4,nn,mm,lp) = amu(4,nn,mm,lp) + v4*dx
      amu(5,nn,mm,lp) = amu(5,nn,mm,lp) + v5*dx
      amu(6,nn,mm,lp) = amu(6,nn,mm,lp) + v6*dx
      dcu(1,nn,mm,lp) = dcu(1,nn,mm,lp) + vx*dx
      dcu(2,nn,mm,lp) = dcu(2,nn,mm,lp) + vy*dx
      dcu(3,nn,mm,lp) = dcu(3,nn,mm,lp) + vz*dx
      cu(1,nn,mm,lp) = cu(1,nn,mm,lp) + ox*dx
      cu(2,nn,mm,lp) = cu(2,nn,mm,lp) + oy*dx
      cu(3,nn,mm,lp) = cu(3,nn,mm,lp) + oz*dx
      dx = dyp*dzp
      amu(1,np,mm,lp) = amu(1,np,mm,lp) + v1*dy
      amu(2,np,mm,lp) = amu(2,np,mm,lp) + v2*dy
      amu(3,np,mm,lp) = amu(3,np,mm,lp) + v3*dy
      amu(4,np,mm,lp) = amu(4,np,mm,lp) + v4*dy
      amu(5,np,mm,lp) = amu(5,np,mm,lp) + v5*dy
      amu(6,np,mm,lp) = amu(6,np,mm,lp) + v6*dy
      dcu(1,np,mm,lp) = dcu(1,np,mm,lp) + vx*dy
      dcu(2,np,mm,lp) = dcu(2,np,mm,lp) + vy*dy
      dcu(3,np,mm,lp) = dcu(3,np,mm,lp) + vz*dy
      cu(1,np,mm,lp) = cu(1,np,mm,lp) + ox*dy
      cu(2,np,mm,lp) = cu(2,np,mm,lp) + oy*dy
      cu(3,np,mm,lp) = cu(3,np,mm,lp) + oz*dy
      dy = dx1*dzp
      amu(1,nn,mp,lp) = amu(1,nn,mp,lp) + v1*dx
      amu(2,nn,mp,lp) = amu(2,nn,mp,lp) + v2*dx
      amu(3,nn,mp,lp) = amu(3,nn,mp,lp) + v3*dx
      amu(4,nn,mp,lp) = amu(4,nn,mp,lp) + v4*dx
      amu(5,nn,mp,lp) = amu(5,nn,mp,lp) + v5*dx
      amu(6,nn,mp,lp) = amu(6,nn,mp,lp) + v6*dx
      dcu(1,nn,mp,lp) = dcu(1,nn,mp,lp) + vx*dx
      dcu(2,nn,mp,lp) = dcu(2,nn,mp,lp) + vy*dx
      dcu(3,nn,mp,lp) = dcu(3,nn,mp,lp) + vz*dx
      cu(1,nn,mp,lp) = cu(1,nn,mp,lp) + ox*dx
      cu(2,nn,mp,lp) = cu(2,nn,mp,lp) + oy*dx
      cu(3,nn,mp,lp) = cu(3,nn,mp,lp) + oz*dx
      amu(1,np,mp,lp) = amu(1,np,mp,lp) + v1*dy
      amu(2,np,mp,lp) = amu(2,np,mp,lp) + v2*dy
      amu(3,np,mp,lp) = amu(3,np,mp,lp) + v3*dy
      amu(4,np,mp,lp) = amu(4,np,mp,lp) + v4*dy
      amu(5,np,mp,lp) = amu(5,np,mp,lp) + v5*dy
      amu(6,np,mp,lp) = amu(6,np,mp,lp) + v6*dy
      dcu(1,np,mp,lp) = dcu(1,np,mp,lp) + vx*dy
      dcu(2,np,mp,lp) = dcu(2,np,mp,lp) + vy*dy
      dcu(3,np,mp,lp) = dcu(3,np,mp,lp) + vz*dy
      cu(1,np,mp,lp) = cu(1,np,mp,lp) + ox*dy
      cu(2,np,mp,lp) = cu(2,np,mp,lp) + oy*dy
      cu(3,np,mp,lp) = cu(3,np,mp,lp) + oz*dy
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DSORTP3YZL(parta,partb,npic,idimp,nop,ny1,nyz1)
c this subroutine sorts particles by y,z grid
c linear interpolation
c part = particle array
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c npic = address offset for reordering particles
c idimp = size of phase space = 6
c nop = number of particles
c ny1 = system length in y direction + 1
c nyz1 = ny1*nz1, where nz1 = system length in z direction + 1
      implicit none
      integer npic, idimp, nop, ny1, nyz1
      real parta, partb
      dimension parta(idimp,nop), partb(idimp,nop), npic(nyz1)
c local data
      integer i, j, k, m, l, isum, ist, ip
c clear counter array
      do 10 k = 1, nyz1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      m = parta(2,j)
      l = parta(3,j)
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
      m = parta(2,j)
      l = parta(3,j)
      l = m + ny1*l + 1
      ip = npic(l) + 1
      do 40 i = 1, idimp
      partb(i,ip) = parta(i,j)
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
      dimension fxyz(3,nxe,nye,nze)
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
      dimension cu(3,nxe,nye,nze)
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
      subroutine AMCGUARD3L(amu,nx,ny,nz,nxe,nye,nze,ndim)
c accumulate extended periodic tensor field amu
c linear interpolation
c nx/ny/nz = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
c nze = third dimension of field arrays, must be >= nz+1
c ndim = number of elements in tensor
      implicit none
      integer nx, ny, nz, nxe, nye, nze, ndim
      real amu
      dimension amu(ndim,nxe,nye,nze)
c local data
      integer i, j, k, l
c accumulate edges of extended field
      do 60 l = 1, nz
      do 20 k = 1, ny
      do 10 i = 1, ndim
      amu(i,1,k,l) = amu(i,1,k,l) + amu(i,nx+1,k,l)
      amu(i,nx+1,k,l) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, ndim
      amu(i,j,1,l) = amu(i,j,1,l) + amu(i,j,ny+1,l)
      amu(i,j,ny+1,l) = 0.0
   30 continue
   40 continue
      do 50 i = 1, ndim
      amu(i,1,1,l) = amu(i,1,1,l) + amu(i,nx+1,ny+1,l)
      amu(i,nx+1,ny+1,l) = 0.0
   50 continue
   60 continue
      do 100 k = 1, ny
      do 80 j = 1, nx
      do 70 i = 1, ndim
      amu(i,j,k,1) = amu(i,j,k,1) + amu(i,j,k,nz+1)
      amu(i,j,k,nz+1) = 0.0
   70 continue
   80 continue
      do 90 i = 1, ndim
      amu(i,1,k,1) = amu(i,1,k,1) + amu(i,nx+1,k,nz+1)
      amu(i,nx+1,k,nz+1) = 0.0
   90 continue
  100 continue
      do 120 j = 1, nx
      do 110 i = 1, ndim
      amu(i,j,1,1) = amu(i,j,1,1) + amu(i,j,ny+1,nz+1)
      amu(i,j,ny+1,nz+1) = 0.0
  110 continue
  120 continue
      do 130 i = 1, ndim
      amu(i,1,1,1) = amu(i,1,1,1) + amu(i,nx+1,ny+1,nz+1)
      amu(i,nx+1,ny+1,nz+1) = 0.0
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine ASCFGUARD3L(dcu,cus,q2m0,nx,ny,nz,nxe,nye,nze)
c add scaled field to extended periodic field
      implicit none
      real dcu, cus, q2m0
      integer nx, ny, nz, nxe, nye, nze
      dimension dcu(3,nxe,nye,nze), cus(3,nxe,nye,nze)
c local data
      integer i, j, k, l
      do 40 l = 1, nz
      do 30 k = 1, ny
      do 20 j = 1, nx
      do 10 i = 1, 3
      dcu(i,j,k,l) = dcu(i,j,k,l) - q2m0*cus(i,j,k,l)
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FWPMINMX3(qe,qbme,wpmax,wpmin,nx,ny,nz,nxe,nye,nze)
c calculates maximum and minimum plasma frequency. assumes guard cells
c have already been added
c qe = charge density for electrons
c qbme = charge/mass ratio for electrons
c wpmax/wpmin = maximum/minimum plasma frequency
c nx/ny/nz = system length in x/y/z direction
c nxe = first dimension of charge arrays, nxe must be >= nx
c nye = second dimension of charge arrays, nye must be >= ny
c nze = third dimension of charge arrays, nze must be >= nz
      implicit none
      real qe, qbme, wpmax, wpmin
      integer nx, ny, nz, nxe, nye, nze
      dimension qe(nxe,nye,nze)
c local data
      integer j, k, l
      real at1
      wpmax = qbme*qe(1,1,1)
      wpmin = wpmax
      do 30 l = 1, nz
      do 20 k = 1, ny
      do 10 j = 1, nx
      at1 = qbme*qe(j,k,l)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine POIS33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxvh,
     1nyv,nzv,nxhd,nyhd,nzhd)
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
      implicit none
      integer isign, nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
      real ax, ay, az, affp, we
      complex q, fxyz, ffc
      dimension q(nxvh,nyv,nzv), fxyz(3,nxvh,nyv,nzv)
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
      do 90 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 60 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
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
      wp = wp + at1*(q(j,k,l)*conjg(q(j,k,l))                           
     1   + q(j,k1,l)*conjg(q(j,k1,l)) + q(j,k,l1)*conjg(q(j,k,l1))      
     2   + q(j,k1,l1)*conjg(q(j,k1,l1)))
   50 continue
   60 continue
c mode numbers kx = 0, nx/2
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
      wp = wp + at1*(q(1,k,l)*conjg(q(1,k,l))                           
     1   + q(1,k,l1)*conjg(q(1,k,l1)))
   70 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
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
      wp = wp + at1*(q(j,1,l)*conjg(q(j,1,l))                           
     1   + q(j,1,l1)*conjg(q(j,1,l1)))
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
      wp = wp + at1*(q(1,1,l)*conjg(q(1,1,l)))
   90 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 110 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
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
      wp = wp + at1*(q(j,k,1)*conjg(q(j,k,1))
     1   + q(j,k1,1)*conjg(q(j,k1,1)))
  100 continue
  110 continue
c mode numbers kx = 0, nx/2
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
      wp = wp + at1*(q(1,k,1)*conjg(q(1,k,1)))
  120 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
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
      wp = wp + at1*(q(j,1,1)*conjg(q(j,1,1)))
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
c nx/ny/nz = system length in x/y/z direction
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
c nzv = fourth dimension of field arrays, must be >= nz
      implicit none
      integer nx, ny, nz, nxvh, nyv, nzv
      complex cu
      dimension cu(3,nxvh,nyv,nzv)
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
      subroutine BBPOIS33(cu,bxyz,ffc,ci,wm,nx,ny,nz,nxvh,nyv,nzv,nxhd, 
     1nyhd,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c with periodic boundary conditions.
c input: cu,ffc,ci,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
c output: bxyz, wm
c approximate flop count is:
c 193*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c the magnetic field is calculated using the equations:
c bx(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (ky*cuz(kx,ky,kz)-kz*cuy(kx,ky,kz))*s(kx,ky,kz),
c by(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kz*cux(kx,ky,kz)-kx*cuz(kx,ky,kz))*s(kx,ky,kz),
c bz(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kx*cuy(kx,ky,kz)-ky*cux(kx,ky,kz))*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bx(ky=pi) = 0,
c bx(kz=pi) = by(kz=pi) = bz(kz=pi) = 0,
c bx(kx=0,ky=0,kz=0) = by(kx=0,ky=0,kz=0) = bz(kx=0,ky=0,kz=0) = 0.
c cu(i,j,k,l) = complex current density for fourier mode (j-1,k-1,l-1)
c bxy(1,j,k,l) = x component of complex magnetic field
c bxy(2,j,k,l) = y component of complex magnetic field
c bxy(3,j,k,l) = z component of complex magnetic field
c all for fourier mode (j-1,k-1,l-1)
c aimag(ffc(j,k,l)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1,l-1)
c real(ffc(j,k,l)) = potential green's function g
c for fourier mode (j-1,k-1,l-1)
c ci = reciprocal of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny/nz = system length in x/y/z direction
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
c nzv = fourth dimension of field arrays, must be >= nz
c nxhd = dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
c nzhd = third dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
      real ci, wm
      complex cu, bxyz, ffc
      dimension cu(3,nxvh,nyv,nzv), bxyz(3,nxvh,nyv,nzv)
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
c calculate smoothed magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 50 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 10 j = 2, nxh
      at1 = ci2*real(ffc(j,k,l))*aimag(ffc(j,k,l))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      at4 = dkz*at1
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
      wp = wp + at1*(cu(1,j,k,l)*conjg(cu(1,j,k,l))
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
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k,l))*aimag(ffc(1,k,l))
      at3 = dny*real(k - 1)*at1
      at4 = dkz*at1
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
      wp = wp + at1*(cu(1,1,k,l)*conjg(cu(1,1,k,l))
     1   + cu(2,1,k,l)*conjg(cu(2,1,k,l))
     2   + cu(3,1,k,l)*conjg(cu(3,1,k,l))
     3   + cu(1,1,k,l1)*conjg(cu(1,1,k,l1))
     4   + cu(2,1,k,l1)*conjg(cu(2,1,k,l1))
     5   + cu(3,1,k,l1)*conjg(cu(3,1,k,l1)))
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      at1 = ci2*real(ffc(j,1,l))*aimag(ffc(j,1,l))
      at2 = dnx*real(j - 1)*at1
      at4 = dkz*at1
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
      wp = wp + at1*(cu(1,j,1,l)*conjg(cu(1,j,1,l))
     1   + cu(2,j,1,l)*conjg(cu(2,j,1,l))
     2   + cu(3,j,1,l)*conjg(cu(3,j,1,l))
     3   + cu(1,j,1,l1)*conjg(cu(1,j,1,l1))
     4   + cu(2,j,1,l1)*conjg(cu(2,j,1,l1))
     5   + cu(3,j,1,l1)*conjg(cu(3,j,1,l1)))
   40 continue
c mode numbers kx = 0, nx/2
      at1 = ci2*real(ffc(1,1,l))*aimag(ffc(1,1,l))
      at4 = dkz*at1
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
      wp = wp + at1*(cu(1,1,1,l)*conjg(cu(1,1,1,l))
     1   + cu(2,1,1,l)*conjg(cu(2,1,1,l))
     2   + cu(3,1,1,l)*conjg(cu(3,1,1,l)))
   50 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 60 j = 2, nxh
      at1 = ci2*real(ffc(j,k,1))*aimag(ffc(j,k,1))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
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
      wp = wp + at1*(cu(1,j,k,1)*conjg(cu(1,j,k,1))
     1   + cu(2,j,k,1)*conjg(cu(2,j,k,1))
     2   + cu(3,j,k,1)*conjg(cu(3,j,k,1))
     3   + cu(1,j,k1,1)*conjg(cu(1,j,k1,1))
     4   + cu(2,j,k1,1)*conjg(cu(2,j,k1,1))
     5   + cu(3,j,k1,1)*conjg(cu(3,j,k1,1)))
   60 continue
   70 continue
c mode numbers kx = 0, nx/2
      do 80 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k,1))*aimag(ffc(1,k,1))
      at3 = dny*real(k - 1)*at1
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
      wp = wp + at1*(cu(1,1,k,1)*conjg(cu(1,1,k,1))
     1   + cu(2,1,k,1)*conjg(cu(2,1,k,1))
     2   + cu(3,1,k,1)*conjg(cu(3,1,k,1)))
   80 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 90 j = 2, nxh
      at1 = ci2*real(ffc(j,1,1))*aimag(ffc(j,1,1))
      at2 = dnx*real(j - 1)*at1
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
      wp = wp + at1*(cu(1,j,1,1)*conjg(cu(1,j,1,1))
     1   + cu(2,j,1,1)*conjg(cu(2,j,1,1))
     2   + cu(3,j,1,1)*conjg(cu(3,j,1,1)))
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
      subroutine BADDEXT3(bxyz,omx,omy,omz,nx,ny,nz,nxe,nye,nze)
c adds constant to magnetic field for 3d code
c bxyz = magnetic field
c omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
c nx/ny/nz = system length in x/y/z direction
c nxe = second dimension of magnetic field array, nxe must be >= nx
c nye = third dimension of magnetic field array, nye must be >= ny
c nze = fourth dimension of magnetic field array, nze must be >= nz
      implicit none
      real bxyz, omx, omy, omz
      integer nx, ny, nz, nxe, nye, nze
      dimension bxyz(3,nxe,nye,nze)
c local data
      integer j, k, l
      do 30 l = 1, nz
      do 20 k = 1, ny
      do 10 j = 1, nx
      bxyz(1,j,k,l) = bxyz(1,j,k,l) + omx
      bxyz(2,j,k,l) = bxyz(2,j,k,l) + omy
      bxyz(3,j,k,l) = bxyz(3,j,k,l) + omz
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DCUPERP3(dcu,amu,nx,ny,nz,nxvh,nyv,nzv)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux
c in 3d with periodic boundary conditions.
c input: all except dcu, output: dcu
c approximate flop count is:
c 220*nxc*nyc*nzc + 70*(nxc*nyc + nxc*nzc + nyc*nzc)
c and (nx/2)*nyc*nzc divides
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky,kz) = -sqrt(-1)*(kx*vx*vx+ky*vx*vy+kz*vx*vz)
c dcu(2,kx,ky,kz) = -sqrt(-1)*(kx*vx*vy+ky*vy*vy+kz*vy*vz)
c dcu(3,kx,ky,kz) = -sqrt(-1)*(kx*vx*vz+ky*vy*vz+kz*vz*vz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c dcux(kx=pi) = dcuy(kx=pi) = dcuz(kx=pi) = 0,
c dcux(ky=pi) = dcuy(ky=pi) = dcux(ky=pi) = 0,
c dcux(kz=pi) = dcuy(kz=pi) = dcuz(kz=pi) = 0,
c dcux(kx=0,ky=0,kz=0) = dcuy(kx=0,ky=0,kz=0) = dcuz(kx=0,ky=0,kz=0) = 0
c the transverse part is calculated using the equation:
c dcu(1,kx,ky,kz) = dcu(1,kx,ky,kz)
c                 - kx*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
c                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c dcu(2,kx,ky,kz) = dcu(2,kx,ky,kz)
c                 - ky*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
c                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c dcu(3,kx,ky,kz) = dcu(3,kx,ky,kz)
c                 - kz*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
c                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c on output:
c dcu(i,j,k,l) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1,l-1)
c amu(1,j,k,l) = xx component of complex momentum flux
c amu(2,j,k,l) = xy component of complex momentum flux
c amu(3,j,k,l) = xz component of complex momentum flux
c amu(4,j,k,l) = yy component of complex momentum flux
c amu(5,j,k,l) = yz component of complex momentum flux
c amu(6,j,k,l) = zz component of complex momentum flux
c all for fourier mode (j-1,k-1,l-1)
c nx/ny/nz = system length in x/y/z direction
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
c nzv = fourth dimension of field arrays, must be >= nz
      implicit none
      integer nx, ny, nz, nxvh, nyv, nzv
      complex dcu, amu
      dimension dcu(3,nxvh,nyv,nzv), amu(6,nxvh,nyv,nzv)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real dnx, dny, dnz, dkx, dky, dkz, dky2, dkz2, dkyz2, at1
      complex zero, zt1, zt2, zt3, zt4, zt5
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
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 50 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
      dkz2 = dkz*dkz
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dkyz2 = dky*dky + dkz2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx + dkyz2)
      zt1 = cmplx(aimag(amu(1,j,k,l)),-real(amu(1,j,k,l)))
      zt2 = cmplx(aimag(amu(2,j,k,l)),-real(amu(2,j,k,l)))
      zt3 = cmplx(aimag(amu(3,j,k,l)),-real(amu(3,j,k,l)))
      zt1 = dkx*zt1 + dky*zt2 + dkz*zt3
      zt4 = cmplx(aimag(amu(4,j,k,l)),-real(amu(4,j,k,l)))
      zt5 = cmplx(aimag(amu(5,j,k,l)),-real(amu(5,j,k,l)))
      zt2 = dkx*zt2 + dky*zt4 + dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,k,l)),-real(amu(6,j,k,l)))
      zt3 = dkx*zt3 + dky*zt5 + dkz*zt4
      zt4 = at1*(dkx*zt1 + dky*zt2 + dkz*zt3)
      dcu(1,j,k,l) = zt1 - dkx*zt4
      dcu(2,j,k,l) = zt2 - dky*zt4
      dcu(3,j,k,l) = zt3 - dkz*zt4
      zt1 = cmplx(aimag(amu(1,j,k1,l)),-real(amu(1,j,k1,l)))
      zt2 = cmplx(aimag(amu(2,j,k1,l)),-real(amu(2,j,k1,l)))
      zt3 = cmplx(aimag(amu(3,j,k1,l)),-real(amu(3,j,k1,l)))
      zt1 = dkx*zt1 - dky*zt2 + dkz*zt3
      zt4 = cmplx(aimag(amu(4,j,k1,l)),-real(amu(4,j,k1,l)))
      zt5 = cmplx(aimag(amu(5,j,k1,l)),-real(amu(5,j,k1,l)))
      zt2 = dkx*zt2 - dky*zt4 + dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,k1,l)),-real(amu(6,j,k1,l)))
      zt3 = dkx*zt3 - dky*zt5 + dkz*zt4
      zt4 = at1*(dkx*zt1 - dky*zt2 + dkz*zt3)
      dcu(1,j,k1,l) = zt1 - dkx*zt4
      dcu(2,j,k1,l) = zt2 + dky*zt4
      dcu(3,j,k1,l) = zt3 - dkz*zt4
      zt1 = cmplx(aimag(amu(1,j,k,l1)),-real(amu(1,j,k,l1)))
      zt2 = cmplx(aimag(amu(2,j,k,l1)),-real(amu(2,j,k,l1)))
      zt3 = cmplx(aimag(amu(3,j,k,l1)),-real(amu(3,j,k,l1)))
      zt1 = dkx*zt1 + dky*zt2 - dkz*zt3
      zt4 = cmplx(aimag(amu(4,j,k,l1)),-real(amu(4,j,k,l1)))
      zt5 = cmplx(aimag(amu(5,j,k,l1)),-real(amu(5,j,k,l1)))
      zt2 = dkx*zt2 + dky*zt4 - dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,k,l1)),-real(amu(6,j,k,l1)))
      zt3 = dkx*zt3 + dky*zt5 - dkz*zt4
      zt4 = at1*(dkx*zt1 + dky*zt2 - dkz*zt3)
      dcu(1,j,k,l1) = zt1 - dkx*zt4
      dcu(2,j,k,l1) = zt2 - dky*zt4
      dcu(3,j,k,l1) = zt3 + dkz*zt4
      zt1 = cmplx(aimag(amu(1,j,k1,l1)),-real(amu(1,j,k1,l1)))
      zt2 = cmplx(aimag(amu(2,j,k1,l1)),-real(amu(2,j,k1,l1)))
      zt3 = cmplx(aimag(amu(3,j,k1,l1)),-real(amu(3,j,k1,l1)))
      zt1 = dkx*zt1 - dky*zt2 - dkz*zt3
      zt4 = cmplx(aimag(amu(4,j,k1,l1)),-real(amu(4,j,k1,l1)))
      zt5 = cmplx(aimag(amu(5,j,k1,l1)),-real(amu(5,j,k1,l1)))
      zt2 = dkx*zt2 - dky*zt4 - dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,k1,l1)),-real(amu(6,j,k1,l1)))
      zt3 = dkx*zt3 - dky*zt5 - dkz*zt4
      zt4 = at1*(dkx*zt1 - dky*zt2 - dkz*zt3)
      dcu(1,j,k1,l1) = zt1 - dkx*zt4
      dcu(2,j,k1,l1) = zt2 + dky*zt4
      dcu(3,j,k1,l1) = zt3 + dkz*zt4
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      at1 = 1.0/(dky*dky + dkz2)
      zt2 = cmplx(aimag(amu(2,1,k,l)),-real(amu(2,1,k,l)))
      zt3 = cmplx(aimag(amu(3,1,k,l)),-real(amu(3,1,k,l)))
      zt1 = dky*zt2 + dkz*zt3
      zt4 = cmplx(aimag(amu(4,1,k,l)),-real(amu(4,1,k,l)))
      zt5 = cmplx(aimag(amu(5,1,k,l)),-real(amu(5,1,k,l)))
      zt2 = dky*zt4 + dkz*zt5
      zt4 = cmplx(aimag(amu(6,1,k,l)),-real(amu(6,1,k,l)))
      zt3 = dky*zt5 + dkz*zt4
      zt4 = at1*(dky*zt2 + dkz*zt3)
      dcu(1,1,k,l) = zt1
      dcu(2,1,k,l) = zt2 - dky*zt4
      dcu(3,1,k,l) = zt3 - dkz*zt4
      dcu(1,1,k1,l) = zero
      dcu(2,1,k1,l) = zero
      dcu(3,1,k1,l) = zero
      zt2 = cmplx(aimag(amu(2,1,k,l1)),-real(amu(2,1,k,l1)))
      zt3 = cmplx(aimag(amu(3,1,k,l1)),-real(amu(3,1,k,l1)))
      zt1 = dky*zt2 - dkz*zt3
      zt4 = cmplx(aimag(amu(4,1,k,l1)),-real(amu(4,1,k,l1)))
      zt5 = cmplx(aimag(amu(5,1,k,l1)),-real(amu(5,1,k,l1)))
      zt2 = dky*zt4 - dkz*zt5
      zt4 = cmplx(aimag(amu(6,1,k,l1)),-real(amu(6,1,k,l1)))
      zt3 = dky*zt5 - dkz*zt4
      zt4 = at1*(dky*zt2 - dkz*zt3)
      dcu(1,1,k,l1) = zt1
      dcu(2,1,k,l1) = zt2 - dky*zt4
      dcu(3,1,k,l1) = zt3 + dkz*zt4
      dcu(1,1,k1,l1) = zero
      dcu(2,1,k1,l1) = zero
      dcu(3,1,k1,l1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx + dkz2)
      zt1 = cmplx(aimag(amu(1,j,1,l)),-real(amu(1,j,1,l)))
      zt2 = cmplx(aimag(amu(2,j,1,l)),-real(amu(2,j,1,l)))
      zt3 = cmplx(aimag(amu(3,j,1,l)),-real(amu(3,j,1,l)))
      zt1 = dkx*zt1 + dkz*zt3
      zt5 = cmplx(aimag(amu(5,j,1,l)),-real(amu(5,j,1,l)))
      zt2 = dkx*zt2 + dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,1,l)),-real(amu(6,j,1,l)))
      zt3 = dkx*zt3 + dkz*zt4
      zt4 = at1*(dkx*zt1 + dkz*zt3)
      dcu(1,j,1,l) = zt1 - dkx*zt4
      dcu(2,j,1,l) = zt2
      dcu(3,j,1,l) = zt3 - dkz*zt4
      dcu(1,j,k1,l) = zero
      dcu(2,j,k1,l) = zero
      dcu(3,j,k1,l) = zero
      zt1 = cmplx(aimag(amu(1,j,1,l1)),-real(amu(1,j,1,l1)))
      zt2 = cmplx(aimag(amu(2,j,1,l1)),-real(amu(2,j,1,l1)))
      zt3 = cmplx(aimag(amu(3,j,1,l1)),-real(amu(3,j,1,l1)))
      zt1 = dkx*zt1 - dkz*zt3
      zt5 = cmplx(aimag(amu(5,j,1,l1)),-real(amu(5,j,1,l1)))
      zt2 = dkx*zt2 - dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,1,l1)),-real(amu(6,j,1,l1)))
      zt3 = dkx*zt3 - dkz*zt4
      zt4 = at1*(dkx*zt1 - dkz*zt3)
      dcu(1,j,1,l1) = zt1 - dkx*zt4
      dcu(2,j,1,l1) = zt2
      dcu(3,j,1,l1) = zt3 + dkz*zt4
      dcu(1,j,k1,l1) = zero
      dcu(2,j,k1,l1) = zero
      dcu(3,j,k1,l1) = zero
   40 continue
c mode numbers kx = 0, nx/2
      zt3 = cmplx(aimag(amu(3,1,1,l)),-real(amu(3,1,1,l)))
      zt1 = dkz*zt3
      zt5 = cmplx(aimag(amu(5,1,1,l)),-real(amu(5,1,1,l)))
      zt2 = dkz*zt5
      dcu(1,1,1,l) = zt1
      dcu(2,1,1,l) = zt2
      dcu(3,1,1,l) = zero
      dcu(1,1,k1,l) = zero
      dcu(2,1,k1,l) = zero
      dcu(3,1,k1,l) = zero
      dcu(1,1,1,l1) = zero
      dcu(2,1,1,l1) = zero
      dcu(3,1,1,l1) = zero
      dcu(1,1,k1,l1) = zero
      dcu(2,1,k1,l1) = zero
      dcu(3,1,k1,l1) = zero
   50 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      do 60 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx + dky2)
      zt1 = cmplx(aimag(amu(1,j,k,1)),-real(amu(1,j,k,1)))
      zt2 = cmplx(aimag(amu(2,j,k,1)),-real(amu(2,j,k,1)))
      zt3 = cmplx(aimag(amu(3,j,k,1)),-real(amu(3,j,k,1)))
      zt1 = dkx*zt1 + dky*zt2
      zt4 = cmplx(aimag(amu(4,j,k,1)),-real(amu(4,j,k,1)))
      zt5 = cmplx(aimag(amu(5,j,k,1)),-real(amu(5,j,k,1)))
      zt2 = dkx*zt2 + dky*zt4
      zt3 = dkx*zt3 + dky*zt5
      zt4 = at1*(dkx*zt1 + dky*zt2)
      dcu(1,j,k,1) = zt1 - dkx*zt4
      dcu(2,j,k,1) = zt2 - dky*zt4
      dcu(3,j,k,1) = zt3
      dcu(1,j,k,l1) = zero
      dcu(2,j,k,l1) = zero
      dcu(3,j,k,l1) = zero
      zt1 = cmplx(aimag(amu(1,j,k1,1)),-real(amu(1,j,k1,1)))
      zt2 = cmplx(aimag(amu(2,j,k1,1)),-real(amu(2,j,k1,1)))
      zt3 = cmplx(aimag(amu(3,j,k1,1)),-real(amu(3,j,k1,1)))
      zt1 = dkx*zt1 - dky*zt2
      zt4 = cmplx(aimag(amu(4,j,k1,1)),-real(amu(4,j,k1,1)))
      zt5 = cmplx(aimag(amu(5,j,k1,1)),-real(amu(5,j,k1,1)))
      zt2 = dkx*zt2 - dky*zt4
      zt3 = dkx*zt3 - dky*zt5
      zt4 = at1*(dkx*zt1 - dky*zt2)
      dcu(1,j,k1,1) = zt1 - dkx*zt4
      dcu(2,j,k1,1) = zt2 + dky*zt4
      dcu(3,j,k1,1) = zt3
      dcu(1,j,k1,l1) = zero
      dcu(2,j,k1,l1) = zero
      dcu(3,j,k1,l1) = zero
   60 continue
   70 continue
c mode numbers kx = 0, nx/2
      do 80 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      zt2 = cmplx(aimag(amu(2,1,k,1)),-real(amu(2,1,k,1)))
      zt1 = dky*zt2
      zt5 = cmplx(aimag(amu(5,1,k,1)),-real(amu(5,1,k,1)))
      zt3 = dky*zt5
      dcu(1,1,k,1) = zt1
      dcu(2,1,k,1) = zero
      dcu(3,1,k,1) = zt3
      dcu(1,1,k1,1) = zero
      dcu(2,1,k1,1) = zero
      dcu(3,1,k1,1) = zero
      dcu(1,1,k,l1) = zero
      dcu(2,1,k,l1) = zero
      dcu(3,1,k,l1) = zero
      dcu(1,1,k1,l1) = zero
      dcu(2,1,k1,l1) = zero
      dcu(3,1,k1,l1) = zero
   80 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 90 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1,1)),-real(amu(2,j,1,1)))
      zt3 = cmplx(aimag(amu(3,j,1,1)),-real(amu(3,j,1,1)))
      zt2 = dkx*zt2
      zt3 = dkx*zt3
      dcu(1,j,1,1) = zero
      dcu(2,j,1,1) = zt2
      dcu(3,j,1,1) = zt3
      dcu(1,j,k1,1) = zero
      dcu(2,j,k1,1) = zero
      dcu(3,j,k1,1) = zero
      dcu(1,j,1,l1) = zero
      dcu(2,j,1,l1) = zero
      dcu(3,j,1,l1) = zero
      dcu(1,j,k1,l1) = zero
      dcu(2,j,k1,l1) = zero
      dcu(3,j,k1,l1) = zero
   90 continue
      dcu(1,1,1,1) = zero
      dcu(2,1,1,1) = zero
      dcu(3,1,1,1) = zero
      dcu(1,1,k1,1) = zero
      dcu(2,1,k1,1) = zero
      dcu(3,1,k1,1) = zero
      dcu(1,1,1,l1) = zero
      dcu(2,1,1,l1) = zero
      dcu(3,1,1,l1) = zero
      dcu(1,1,k1,l1) = zero
      dcu(2,1,k1,l1) = zero
      dcu(3,1,k1,l1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine ADCUPERP3(dcu,amu,nx,ny,nz,nxvh,nyv,nzv)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux and acceleration density
c in 3d with periodic boundary conditions.
c input: all, output: dcu
c approximate flop count is:
c 244*nxc*nyc*nzc + 82*(nxc*nyc + nxc*nzc + nyc*nzc)
c and (nx/2)*nyc*nzc divides
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky,kz) = dcu(1,kx,ky,kz)
c                -sqrt(-1)*(kx*vx*vx+ky*vx*vy+kz*vx*vz)
c dcu(2,kx,ky,kz) = dcu(2,kx,ky,kz)
c                -sqrt(-1)*(kx*vx*vy+ky*vy*vy+kz*vy*vz)
c dcu(3,kx,ky,kz) = dcu(3,kx,ky,kz)
c                 -sqrt(-1)*(kx*vx*vz+ky*vy*vz+kz*vz*vz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c dcux(kx=pi) = dcuy(kx=pi) = dcuz(kx=pi) = 0,
c dcux(ky=pi) = dcuy(ky=pi) = dcux(ky=pi) = 0,
c dcux(kz=pi) = dcuy(kz=pi) = dcuz(kz=pi) = 0,
c dcux(kx=0,ky=0,kz=0) = dcuy(kx=0,ky=0,kz=0) = dcuz(kx=0,ky=0,kz=0) = 0
c the transverse part is calculated using the equation:
c dcu(1,kx,ky,kz) = dcu(1,kx,ky,kz)
c                 - kx*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
c                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c dcu(2,kx,ky,kz) = dcu(2,kx,ky,kz)
c                 - ky*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
c                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c dcu(3,kx,ky,kz) = dcu(3,kx,ky,kz)
c                 - kz*(kx*dcu(1,kx,ky,kz)+ky*dcu(2,kx,ky,kz)
c                      +kz*dcu(3,kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c on input:
c dcu(i,j,k,l) = complex acceleration density for
c fourier mode (j-1,k-1,l-1)
c on output:
c dcu(i,j,k,l) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1,l-1)
c amu(1,j,k,l) = xx component of complex momentum flux
c amu(2,j,k,l) = xy component of complex momentum flux
c amu(3,j,k,l) = xz component of complex momentum flux
c amu(4,j,k,l) = yy component of complex momentum flux
c amu(5,j,k,l) = yz component of complex momentum flux
c amu(6,j,k,l) = zz component of complex momentum flux
c all for fourier mode (j-1,k-1,l-1)
c nx/ny/nz = system length in x/y/z direction
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
c nzv = fourth dimension of field arrays, must be >= nz
      implicit none
      integer nx, ny, nz, nxvh, nyv, nzv
      complex dcu, amu
      dimension dcu(3,nxvh,nyv,nzv), amu(6,nxvh,nyv,nzv)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real dnx, dny, dnz, dkx, dky, dkz, dky2, dkz2, dkyz2, at1
      complex zero, zt1, zt2, zt3, zt4, zt5
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
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 50 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
      dkz2 = dkz*dkz
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dkyz2 = dky*dky + dkz2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx + dkyz2)
      zt1 = cmplx(aimag(amu(1,j,k,l)),-real(amu(1,j,k,l)))
      zt2 = cmplx(aimag(amu(2,j,k,l)),-real(amu(2,j,k,l)))
      zt3 = cmplx(aimag(amu(3,j,k,l)),-real(amu(3,j,k,l)))
      zt1 = dcu(1,j,k,l) + dkx*zt1 + dky*zt2 + dkz*zt3
      zt4 = cmplx(aimag(amu(4,j,k,l)),-real(amu(4,j,k,l)))
      zt5 = cmplx(aimag(amu(5,j,k,l)),-real(amu(5,j,k,l)))
      zt2 = dcu(2,j,k,l) + dkx*zt2 + dky*zt4 + dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,k,l)),-real(amu(6,j,k,l)))
      zt3 = dcu(3,j,k,l) + dkx*zt3 + dky*zt5 + dkz*zt4
      zt4 = at1*(dkx*zt1 + dky*zt2 + dkz*zt3)
      dcu(1,j,k,l) = zt1 - dkx*zt4
      dcu(2,j,k,l) = zt2 - dky*zt4
      dcu(3,j,k,l) = zt3 - dkz*zt4
      zt1 = cmplx(aimag(amu(1,j,k1,l)),-real(amu(1,j,k1,l)))
      zt2 = cmplx(aimag(amu(2,j,k1,l)),-real(amu(2,j,k1,l)))
      zt3 = cmplx(aimag(amu(3,j,k1,l)),-real(amu(3,j,k1,l)))
      zt1 = dcu(1,j,k1,l) + dkx*zt1 - dky*zt2 + dkz*zt3
      zt4 = cmplx(aimag(amu(4,j,k1,l)),-real(amu(4,j,k1,l)))
      zt5 = cmplx(aimag(amu(5,j,k1,l)),-real(amu(5,j,k1,l)))
      zt2 = dcu(2,j,k1,l) + dkx*zt2 - dky*zt4 + dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,k1,l)),-real(amu(6,j,k1,l)))
      zt3 = dcu(3,j,k1,l) + dkx*zt3 - dky*zt5 + dkz*zt4
      zt4 = at1*(dkx*zt1 - dky*zt2 + dkz*zt3)
      dcu(1,j,k1,l) = zt1 - dkx*zt4
      dcu(2,j,k1,l) = zt2 + dky*zt4
      dcu(3,j,k1,l) = zt3 - dkz*zt4
      zt1 = cmplx(aimag(amu(1,j,k,l1)),-real(amu(1,j,k,l1)))
      zt2 = cmplx(aimag(amu(2,j,k,l1)),-real(amu(2,j,k,l1)))
      zt3 = cmplx(aimag(amu(3,j,k,l1)),-real(amu(3,j,k,l1)))
      zt1 = dcu(1,j,k,l1) + dkx*zt1 + dky*zt2 - dkz*zt3
      zt4 = cmplx(aimag(amu(4,j,k,l1)),-real(amu(4,j,k,l1)))
      zt5 = cmplx(aimag(amu(5,j,k,l1)),-real(amu(5,j,k,l1)))
      zt2 = dcu(2,j,k,l1) + dkx*zt2 + dky*zt4 - dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,k,l1)),-real(amu(6,j,k,l1)))
      zt3 = dcu(3,j,k,l1) + dkx*zt3 + dky*zt5 - dkz*zt4
      zt4 = at1*(dkx*zt1 + dky*zt2 - dkz*zt3)
      dcu(1,j,k,l1) = zt1 - dkx*zt4
      dcu(2,j,k,l1) = zt2 - dky*zt4
      dcu(3,j,k,l1) = zt3 + dkz*zt4
      zt1 = cmplx(aimag(amu(1,j,k1,l1)),-real(amu(1,j,k1,l1)))
      zt2 = cmplx(aimag(amu(2,j,k1,l1)),-real(amu(2,j,k1,l1)))
      zt3 = cmplx(aimag(amu(3,j,k1,l1)),-real(amu(3,j,k1,l1)))
      zt1 = dcu(1,j,k1,l1) + dkx*zt1 - dky*zt2 - dkz*zt3
      zt4 = cmplx(aimag(amu(4,j,k1,l1)),-real(amu(4,j,k1,l1)))
      zt5 = cmplx(aimag(amu(5,j,k1,l1)),-real(amu(5,j,k1,l1)))
      zt2 = dcu(2,j,k1,l1) + dkx*zt2 - dky*zt4 - dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,k1,l1)),-real(amu(6,j,k1,l1)))
      zt3 = dcu(3,j,k1,l1) + dkx*zt3 - dky*zt5 - dkz*zt4
      zt4 = at1*(dkx*zt1 - dky*zt2 - dkz*zt3)
      dcu(1,j,k1,l1) = zt1 - dkx*zt4
      dcu(2,j,k1,l1) = zt2 + dky*zt4
      dcu(3,j,k1,l1) = zt3 + dkz*zt4
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      at1 = 1.0/(dky*dky + dkz2)
      zt2 = cmplx(aimag(amu(2,1,k,l)),-real(amu(2,1,k,l)))
      zt3 = cmplx(aimag(amu(3,1,k,l)),-real(amu(3,1,k,l)))
      zt1 = dcu(1,1,k,l) + dky*zt2 + dkz*zt3
      zt4 = cmplx(aimag(amu(4,1,k,l)),-real(amu(4,1,k,l)))
      zt5 = cmplx(aimag(amu(5,1,k,l)),-real(amu(5,1,k,l)))
      zt2 = dcu(2,1,k,l) + dky*zt4 + dkz*zt5
      zt4 = cmplx(aimag(amu(6,1,k,l)),-real(amu(6,1,k,l)))
      zt3 = dcu(3,1,k,l) + dky*zt5 + dkz*zt4
      zt4 = at1*(dky*zt2 + dkz*zt3)
      dcu(1,1,k,l) = zt1
      dcu(2,1,k,l) = zt2 - dky*zt4
      dcu(3,1,k,l) = zt3 - dkz*zt4
      dcu(1,1,k1,l) = zero
      dcu(2,1,k1,l) = zero
      dcu(3,1,k1,l) = zero
      zt2 = cmplx(aimag(amu(2,1,k,l1)),-real(amu(2,1,k,l1)))
      zt3 = cmplx(aimag(amu(3,1,k,l1)),-real(amu(3,1,k,l1)))
      zt1 = dcu(1,1,k,l1) + dky*zt2 - dkz*zt3
      zt4 = cmplx(aimag(amu(4,1,k,l1)),-real(amu(4,1,k,l1)))
      zt5 = cmplx(aimag(amu(5,1,k,l1)),-real(amu(5,1,k,l1)))
      zt2 = dcu(2,1,k,l1) + dky*zt4 - dkz*zt5
      zt4 = cmplx(aimag(amu(6,1,k,l1)),-real(amu(6,1,k,l1)))
      zt3 = dcu(3,1,k,l1) + dky*zt5 - dkz*zt4
      zt4 = at1*(dky*zt2 - dkz*zt3)
      dcu(1,1,k,l1) = zt1
      dcu(2,1,k,l1) = zt2 - dky*zt4
      dcu(3,1,k,l1) = zt3 + dkz*zt4
      dcu(1,1,k1,l1) = zero
      dcu(2,1,k1,l1) = zero
      dcu(3,1,k1,l1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx + dkz2)
      zt1 = cmplx(aimag(amu(1,j,1,l)),-real(amu(1,j,1,l)))
      zt2 = cmplx(aimag(amu(2,j,1,l)),-real(amu(2,j,1,l)))
      zt3 = cmplx(aimag(amu(3,j,1,l)),-real(amu(3,j,1,l)))
      zt1 = dcu(1,j,1,l) + dkx*zt1 + dkz*zt3
      zt5 = cmplx(aimag(amu(5,j,1,l)),-real(amu(5,j,1,l)))
      zt2 = dcu(2,j,1,l) + dkx*zt2 + dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,1,l)),-real(amu(6,j,1,l)))
      zt3 = dcu(3,j,1,l) + dkx*zt3 + dkz*zt4
      zt4 = at1*(dkx*zt1 + dkz*zt3)
      dcu(1,j,1,l) = zt1 - dkx*zt4
      dcu(2,j,1,l) = zt2
      dcu(3,j,1,l) = zt3 - dkz*zt4
      dcu(1,j,k1,l) = zero
      dcu(2,j,k1,l) = zero
      dcu(3,j,k1,l) = zero
      zt1 = cmplx(aimag(amu(1,j,1,l1)),-real(amu(1,j,1,l1)))
      zt2 = cmplx(aimag(amu(2,j,1,l1)),-real(amu(2,j,1,l1)))
      zt3 = cmplx(aimag(amu(3,j,1,l1)),-real(amu(3,j,1,l1)))
      zt1 = dcu(1,j,1,l1) + dkx*zt1 - dkz*zt3
      zt5 = cmplx(aimag(amu(5,j,1,l1)),-real(amu(5,j,1,l1)))
      zt2 = dcu(2,j,1,l1) + dkx*zt2 - dkz*zt5
      zt4 = cmplx(aimag(amu(6,j,1,l1)),-real(amu(6,j,1,l1)))
      zt3 = dcu(3,j,1,l1) + dkx*zt3 - dkz*zt4
      zt4 = at1*(dkx*zt1 - dkz*zt3)
      dcu(1,j,1,l1) = zt1 - dkx*zt4
      dcu(2,j,1,l1) = zt2
      dcu(3,j,1,l1) = zt3 + dkz*zt4
      dcu(1,j,k1,l1) = zero
      dcu(2,j,k1,l1) = zero
      dcu(3,j,k1,l1) = zero
   40 continue
c mode numbers kx = 0, nx/2
      zt3 = cmplx(aimag(amu(3,1,1,l)),-real(amu(3,1,1,l)))
      zt1 = dcu(1,1,1,l) + dkz*zt3
      zt5 = cmplx(aimag(amu(5,1,1,l)),-real(amu(5,1,1,l)))
      zt2 = dcu(2,1,1,l) + dkz*zt5
      dcu(1,1,1,l) = zt1
      dcu(2,1,1,l) = zt2
      dcu(3,1,1,l) = zero
      dcu(1,1,k1,l) = zero
      dcu(2,1,k1,l) = zero
      dcu(3,1,k1,l) = zero
      dcu(1,1,1,l1) = zero
      dcu(2,1,1,l1) = zero
      dcu(3,1,1,l1) = zero
      dcu(1,1,k1,l1) = zero
      dcu(2,1,k1,l1) = zero
      dcu(3,1,k1,l1) = zero
   50 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      do 60 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx + dky2)
      zt1 = cmplx(aimag(amu(1,j,k,1)),-real(amu(1,j,k,1)))
      zt2 = cmplx(aimag(amu(2,j,k,1)),-real(amu(2,j,k,1)))
      zt3 = cmplx(aimag(amu(3,j,k,1)),-real(amu(3,j,k,1)))
      zt1 = dcu(1,j,k,1) + dkx*zt1 + dky*zt2
      zt4 = cmplx(aimag(amu(4,j,k,1)),-real(amu(4,j,k,1)))
      zt5 = cmplx(aimag(amu(5,j,k,1)),-real(amu(5,j,k,1)))
      zt2 = dcu(2,j,k,1) + dkx*zt2 + dky*zt4
      zt3 = dcu(3,j,k,1) + dkx*zt3 + dky*zt5
      zt4 = at1*(dkx*zt1 + dky*zt2)
      dcu(1,j,k,1) = zt1 - dkx*zt4
      dcu(2,j,k,1) = zt2 - dky*zt4
      dcu(3,j,k,1) = zt3
      dcu(1,j,k,l1) = zero
      dcu(2,j,k,l1) = zero
      dcu(3,j,k,l1) = zero
      zt1 = cmplx(aimag(amu(1,j,k1,1)),-real(amu(1,j,k1,1)))
      zt2 = cmplx(aimag(amu(2,j,k1,1)),-real(amu(2,j,k1,1)))
      zt3 = cmplx(aimag(amu(3,j,k1,1)),-real(amu(3,j,k1,1)))
      zt1 = dcu(1,j,k1,1) + dkx*zt1 - dky*zt2
      zt4 = cmplx(aimag(amu(4,j,k1,1)),-real(amu(4,j,k1,1)))
      zt5 = cmplx(aimag(amu(5,j,k1,1)),-real(amu(5,j,k1,1)))
      zt2 = dcu(2,j,k1,1) + dkx*zt2 - dky*zt4
      zt3 = dcu(3,j,k1,1) + dkx*zt3 - dky*zt5
      zt4 = at1*(dkx*zt1 - dky*zt2)
      dcu(1,j,k1,1) = zt1 - dkx*zt4
      dcu(2,j,k1,1) = zt2 + dky*zt4
      dcu(3,j,k1,1) = zt3
      dcu(1,j,k1,l1) = zero
      dcu(2,j,k1,l1) = zero
      dcu(3,j,k1,l1) = zero
   60 continue
   70 continue
c mode numbers kx = 0, nx/2
      do 80 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      zt2 = cmplx(aimag(amu(2,1,k,1)),-real(amu(2,1,k,1)))
      zt1 = dcu(1,1,k,1) + dky*zt2
      zt5 = cmplx(aimag(amu(5,1,k,1)),-real(amu(5,1,k,1)))
      zt3 = dcu(3,1,k,1) + dky*zt5
      dcu(1,1,k,1) = zt1
      dcu(2,1,k,1) = zero
      dcu(3,1,k,1) = zt3
      dcu(1,1,k1,1) = zero
      dcu(2,1,k1,1) = zero
      dcu(3,1,k1,1) = zero
      dcu(1,1,k,l1) = zero
      dcu(2,1,k,l1) = zero
      dcu(3,1,k,l1) = zero
      dcu(1,1,k1,l1) = zero
      dcu(2,1,k1,l1) = zero
      dcu(3,1,k1,l1) = zero
   80 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 90 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1,1)),-real(amu(2,j,1,1)))
      zt3 = cmplx(aimag(amu(3,j,1,1)),-real(amu(3,j,1,1)))
      zt2 = dcu(2,j,1,1) + dkx*zt2
      zt3 = dcu(3,j,1,1) + dkx*zt3
      dcu(1,j,1,1) = zero
      dcu(2,j,1,1) = zt2
      dcu(3,j,1,1) = zt3
      dcu(1,j,k1,1) = zero
      dcu(2,j,k1,1) = zero
      dcu(3,j,k1,1) = zero
      dcu(1,j,1,l1) = zero
      dcu(2,j,1,l1) = zero
      dcu(3,j,1,l1) = zero
      dcu(1,j,k1,l1) = zero
      dcu(2,j,k1,l1) = zero
      dcu(3,j,k1,l1) = zero
   90 continue
      dcu(1,1,1,1) = zero
      dcu(2,1,1,1) = zero
      dcu(3,1,1,1) = zero
      dcu(1,1,k1,1) = zero
      dcu(2,1,k1,1) = zero
      dcu(3,1,k1,1) = zero
      dcu(1,1,1,l1) = zero
      dcu(2,1,1,l1) = zero
      dcu(3,1,1,l1) = zero
      dcu(1,1,k1,l1) = zero
      dcu(2,1,k1,l1) = zero
      dcu(3,1,k1,l1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine EPOIS33(dcu,exyz,isign,ffe,ax,ay,az,affp,wp0,ci,wf,nx, 
     1ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c transverse electric field (or convolution of transverse electric field
c over particle shape), with periodic boundary conditions.
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c for isign = 0, output: ffe
c input: isign,ax,ay,az,affp,wp0,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
c for isign /= 0, output: exyz, wf
c input: dcu,ffe,isign,ci,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
c approximate flop count is:
c 128*nxc*nyc*nzc + 66*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c if isign = -1, smoothed transverse electric field is calculated
c using the equation:
c ex(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcux(kx,ky,kz)*s(kx,ky,kz)
c ey(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcuy(kx,ky,kz)*s(kx,ky,kz)
c ez(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcuz(kx,ky,kz)*s(kx,ky,kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
c ex(kz=pi) = ey(kz=pi) = ez(kz=pi) = 0,
c ex(kx=0,ky=0,kz=0) = ey(kx=0,ky=0,kz=0) = ez(kx=0,ky=0,kz=0) = 0.
c if isign = 1, unsmoothed transverse electric field is calculated
c using the equation:
c ex(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcux(kx,ky,kz)
c ey(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcuy(kx,ky,kz)
c ez(kx,ky,kz) = -ci*ci*g(kx,ky,kz)*dcuz(kx,ky,kz)
c dcu(i,j,k,l) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1,l-1)
c exyz(1,j,k,l) = x component of complex transverse electric field
c exyz(2,j,k,l) = y component of complex transverse electric field
c exyz(3,j,k,l) = z component of complex transverse electric field
c all for fourier mode (j-1,k-1,l-1)
c aimag(ffe(j,k,l)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1,l-1)
c real(ffe(j,k,l)) = potential green's function g
c for fourier mode (j-1,k-1,l-1)
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c wp0 = normalized total plasma frequency squared
c where np=number of particles
c ci = reciprocal of velocity of light
c transverse electric field energy is also calculated, using
c wf = nx*ny*nz*sum((affp/((kx**2+ky**2+kz**2)*ci*ci)**2)
c    |dcu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the derivative of current is
c divergence-free
c nx/ny/nz = system length in x/y/z direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nzv = third dimension of field arrays, must be >= nz
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
c nzhd = third dimension of form factor array, must be >= nzh
      implicit none
      integer isign, nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
      real ax, ay, az, affp, wp0, ci, wf
      complex dcu, exyz, ffe
      dimension dcu(3,nxvh,nyv,nzv), exyz(3,nxvh,nyv,nzv)
      dimension ffe(nxhd,nyhd,nzhd)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real dnx, dny, dnz, dkx, dky, dkz, ci2, wpc
      real at1, at2, at3, at4, at5, at6
      complex zero
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
      if (isign.ne.0) go to 40
      wpc = wp0*ci2
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
      if (at5.eq.0.0) then
         ffe(j,k,l) = cmplx(affp,1.0)
      else
         ffe(j,k,l) = cmplx(affp*at6/(at5 + wpc*at6*at6),at6)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 140
c calculate smoothed transverse electric field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 90 l = 2, nzh
      l1 = nz2 - l
      do 60 k = 2, nyh
      k1 = ny2 - k
      do 50 j = 2, nxh
      at2 = -ci2*real(ffe(j,k,l))
      at1 = at2*aimag(ffe(j,k,l))
      at2 = at2*at2
      exyz(1,j,k,l) = at1*dcu(1,j,k,l)
      exyz(2,j,k,l) = at1*dcu(2,j,k,l)
      exyz(3,j,k,l) = at1*dcu(3,j,k,l)
      exyz(1,j,k1,l) = at1*dcu(1,j,k1,l)
      exyz(2,j,k1,l) = at1*dcu(2,j,k1,l)
      exyz(3,j,k1,l) = at1*dcu(3,j,k1,l)
      exyz(1,j,k,l1) = at1*dcu(1,j,k,l1)
      exyz(2,j,k,l1) = at1*dcu(2,j,k,l1)
      exyz(3,j,k,l1) = at1*dcu(3,j,k,l1)
      exyz(1,j,k1,l1) = at1*dcu(1,j,k1,l1)
      exyz(2,j,k1,l1) = at1*dcu(2,j,k1,l1)
      exyz(3,j,k1,l1) = at1*dcu(3,j,k1,l1)
      wp = wp + at2*(dcu(1,j,k,l)*conjg(dcu(1,j,k,l))                   
     1   + dcu(2,j,k,l)*conjg(dcu(2,j,k,l))                             
     2   + dcu(3,j,k,l)*conjg(dcu(3,j,k,l))                             
     3   + dcu(1,j,k1,l)*conjg(dcu(1,j,k1,l))                           
     4   + dcu(2,j,k1,l)*conjg(dcu(2,j,k1,l))                           
     5   + dcu(3,j,k1,l)*conjg(dcu(3,j,k1,l))                           
     6   + dcu(1,j,k,l1)*conjg(dcu(1,j,k,l1))                           
     7   + dcu(2,j,k,l1)*conjg(dcu(2,j,k,l1))                           
     8   + dcu(3,j,k,l1)*conjg(dcu(3,j,k,l1))                           
     9   + dcu(1,j,k1,l1)*conjg(dcu(1,j,k1,l1))                         
     a   + dcu(2,j,k1,l1)*conjg(dcu(2,j,k1,l1))                         
     b   + dcu(3,j,k1,l1)*conjg(dcu(3,j,k1,l1)))
   50 continue
   60 continue
c mode numbers kx = 0, nx/2
      do 70 k = 2, nyh
      k1 = ny2 - k
      at2 = -ci2*real(ffe(1,k,l))
      at1 = at2*aimag(ffe(1,k,l))
      at2 = at2*at2
      exyz(1,1,k,l) = at1*dcu(1,1,k,l)
      exyz(2,1,k,l) = at1*dcu(2,1,k,l)
      exyz(3,1,k,l) = at1*dcu(3,1,k,l)
      exyz(1,1,k1,l) = zero
      exyz(2,1,k1,l) = zero
      exyz(3,1,k1,l) = zero
      exyz(1,1,k,l1) = at1*dcu(1,1,k,l1)
      exyz(2,1,k,l1) = at1*dcu(2,1,k,l1)
      exyz(3,1,k,l1) = at1*dcu(3,1,k,l1)
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
      wp = wp + at2*(dcu(1,1,k,l)*conjg(dcu(1,1,k,l))                   
     1   + dcu(2,1,k,l)*conjg(dcu(2,1,k,l))                             
     2   + dcu(3,1,k,l)*conjg(dcu(3,1,k,l))                             
     3   + dcu(1,1,k,l1)*conjg(dcu(1,1,k,l1))                           
     4   + dcu(2,1,k,l1)*conjg(dcu(2,1,k,l1))                           
     5   + dcu(3,1,k,l1)*conjg(dcu(3,1,k,l1)))
   70 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 80 j = 2, nxh
      at2 = -ci2*real(ffe(j,1,l))
      at1 = at2*aimag(ffe(j,1,l))
      at2 = at2*at2
      exyz(1,j,1,l) = at1*dcu(1,j,1,l)
      exyz(2,j,1,l) = at1*dcu(2,j,1,l)
      exyz(3,j,1,l) = at1*dcu(3,j,1,l)
      exyz(1,j,k1,l) = zero
      exyz(2,j,k1,l) = zero
      exyz(3,j,k1,l) = zero
      exyz(1,j,1,l1) = at1*dcu(1,j,1,l1)
      exyz(2,j,1,l1) = at1*dcu(2,j,1,l1)
      exyz(3,j,1,l1) = at1*dcu(3,j,1,l1)
      exyz(1,j,k1,l1) = zero
      exyz(2,j,k1,l1) = zero
      exyz(3,j,k1,l1) = zero
      wp = wp + at2*(dcu(1,j,1,l)*conjg(dcu(1,j,1,l))                   
     1   + dcu(2,j,1,l)*conjg(dcu(2,j,1,l))                             
     2   + dcu(3,j,1,l)*conjg(dcu(3,j,1,l))                             
     3   + dcu(1,j,1,l1)*conjg(dcu(1,j,1,l1))                           
     4   + dcu(2,j,1,l1)*conjg(dcu(2,j,1,l1))                           
     5   + dcu(3,j,1,l1)*conjg(dcu(3,j,1,l1)))
   80 continue
c mode numbers kx = 0, nx/2
      at2 = -ci2*real(ffe(1,1,l))
      at1 = at2*aimag(ffe(1,1,l))
      at2 = at2*at2
      exyz(1,1,1,l) = at1*dcu(1,1,1,l)
      exyz(2,1,1,l) = at1*dcu(2,1,1,l)
      exyz(3,1,1,l) = at1*dcu(3,1,1,l)
      exyz(1,1,k1,l) = zero
      exyz(2,1,k1,l) = zero
      exyz(3,1,k1,l) = zero
      exyz(1,1,1,l1) = zero
      exyz(2,1,1,l1) = zero
      exyz(3,1,1,l1) = zero
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
      wp = wp + at2*(dcu(1,1,1,l)*conjg(dcu(1,1,1,l))                   
     1   + dcu(2,1,1,l)*conjg(dcu(2,1,1,l))                             
     2   + dcu(3,1,1,l)*conjg(dcu(3,1,1,l)))
   90 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 110 k = 2, nyh
      k1 = ny2 - k
      do 100 j = 2, nxh
      at2 = -ci2*real(ffe(j,k,1))
      at1 = at2*aimag(ffe(j,k,1))
      at2 = at2*at2
      exyz(1,j,k,1) = at1*dcu(1,j,k,1)
      exyz(2,j,k,1) = at1*dcu(2,j,k,1)
      exyz(3,j,k,1) = at1*dcu(3,j,k,1)
      exyz(1,j,k1,1) = at1*dcu(1,j,k1,1)
      exyz(2,j,k1,1) = at1*dcu(2,j,k1,1)
      exyz(3,j,k1,1) = at1*dcu(3,j,k1,1)
      exyz(1,j,k,l1) = zero
      exyz(2,j,k,l1) = zero
      exyz(3,j,k,l1) = zero
      exyz(1,j,k1,l1) = zero
      exyz(2,j,k1,l1) = zero
      exyz(3,j,k1,l1) = zero
      wp = wp + at2*(dcu(1,j,k,1)*conjg(dcu(1,j,k,1))                   
     1   + dcu(2,j,k,1)*conjg(dcu(2,j,k,1))                             
     2   + dcu(3,j,k,1)*conjg(dcu(3,j,k,1))                             
     3   + dcu(1,j,k1,1)*conjg(dcu(1,j,k1,1))                           
     4   + dcu(2,j,k1,1)*conjg(dcu(2,j,k1,1))                           
     5   + dcu(3,j,k1,1)*conjg(dcu(3,j,k1,1)))
  100 continue
  110 continue
c mode numbers kx = 0, nx/2
      do 120 k = 2, nyh
      k1 = ny2 - k
      at2 = -ci2*real(ffe(1,k,1))
      at1 = at2*aimag(ffe(1,k,1))
      at2 = at2*at2
      exyz(1,1,k,1) = at1*dcu(1,1,k,1)
      exyz(2,1,k,1) = at1*dcu(2,1,k,1)
      exyz(3,1,k,1) = at1*dcu(3,1,k,1)
      exyz(1,1,k1,1) = zero
      exyz(2,1,k1,1) = zero
      exyz(3,1,k1,1) = zero
      exyz(1,1,k,l1) = zero
      exyz(2,1,k,l1) = zero
      exyz(3,1,k,l1) = zero
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
      wp = wp + at2*(dcu(1,1,k,1)*conjg(dcu(1,1,k,1))                   
     1   + dcu(2,1,k,1)*conjg(dcu(2,1,k,1))                             
     2   + dcu(3,1,k,1)*conjg(dcu(3,1,k,1)))
  120 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 130 j = 2, nxh
      at2 = -ci2*real(ffe(j,1,1))
      at1 = at2*aimag(ffe(j,1,1))
      at2 = at2*at2
      exyz(1,j,1,1) = at1*dcu(1,j,1,1)
      exyz(2,j,1,1) = at1*dcu(2,j,1,1)
      exyz(3,j,1,1) = at1*dcu(3,j,1,1)
      exyz(1,j,k1,1) = zero
      exyz(2,j,k1,1) = zero
      exyz(3,j,k1,1) = zero
      exyz(1,j,1,l1) = zero
      exyz(2,j,1,l1) = zero
      exyz(3,j,1,l1) = zero
      exyz(1,j,k1,l1) = zero
      exyz(2,j,k1,l1) = zero
      exyz(3,j,k1,l1) = zero
      wp = wp + at2*(dcu(1,j,1,1)*conjg(dcu(1,j,1,1))                   
     1   + dcu(2,j,1,1)*conjg(dcu(2,j,1,1))                             
     2   + dcu(3,j,1,1)*conjg(dcu(3,j,1,1)))
  130 continue
      exyz(1,1,1,1) = zero
      exyz(2,1,1,1) = zero
      exyz(3,1,1,1) = zero
      exyz(1,1,k1,1) = zero
      exyz(2,1,k1,1) = zero
      exyz(3,1,k1,1) = zero
      exyz(1,1,1,l1) = zero
      exyz(2,1,1,l1) = zero
      exyz(3,1,1,l1) = zero
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
      wf = real(nx*ny*nz)*wp/real(ffe(1,1,1))
      return
c calculate unsmoothed transverse electric field and sum field energy
  140 wp = 0.0d0
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 190 l = 2, nzh
      l1 = nz2 - l
      do 160 k = 2, nyh
      k1 = ny2 - k
      do 150 j = 2, nxh
      at2 = -ci2*real(ffe(j,k,l))
      at1 = at2*at2
      exyz(1,j,k,l) = at2*dcu(1,j,k,l)
      exyz(2,j,k,l) = at2*dcu(2,j,k,l)
      exyz(3,j,k,l) = at2*dcu(3,j,k,l)
      exyz(1,j,k1,l) = at2*dcu(1,j,k1,l)
      exyz(2,j,k1,l) = at2*dcu(2,j,k1,l)
      exyz(3,j,k1,l) = at2*dcu(3,j,k1,l)
      exyz(1,j,k,l1) = at2*dcu(1,j,k,l1)
      exyz(2,j,k,l1) = at2*dcu(2,j,k,l1)
      exyz(3,j,k,l1) = at2*dcu(3,j,k,l1)
      exyz(1,j,k1,l1) = at2*dcu(1,j,k1,l1)
      exyz(2,j,k1,l1) = at2*dcu(2,j,k1,l1)
      exyz(3,j,k1,l1) = at2*dcu(3,j,k1,l1)
      wp = wp + at1*(dcu(1,j,k,l)*conjg(dcu(1,j,k,l))                   
     1   + dcu(2,j,k,l)*conjg(dcu(2,j,k,l))                             
     2   + dcu(3,j,k,l)*conjg(dcu(3,j,k,l))                             
     3   + dcu(1,j,k1,l)*conjg(dcu(1,j,k1,l))                           
     4   + dcu(2,j,k1,l)*conjg(dcu(2,j,k1,l))                           
     5   + dcu(3,j,k1,l)*conjg(dcu(3,j,k1,l))                           
     6   + dcu(1,j,k,l1)*conjg(dcu(1,j,k,l1))                           
     7   + dcu(2,j,k,l1)*conjg(dcu(2,j,k,l1))                           
     8   + dcu(3,j,k,l1)*conjg(dcu(3,j,k,l1))                           
     9   + dcu(1,j,k1,l1)*conjg(dcu(1,j,k1,l1))                         
     a   + dcu(2,j,k1,l1)*conjg(dcu(2,j,k1,l1))                         
     b   + dcu(3,j,k1,l1)*conjg(dcu(3,j,k1,l1)))
  150 continue
  160 continue
c mode numbers kx = 0, nx/2
      do 170 k = 2, nyh
      k1 = ny2 - k
      at2 = -ci2*real(ffe(1,k,l))
      at1 = at2*at2
      exyz(1,1,k,l) = at2*dcu(1,1,k,l)
      exyz(2,1,k,l) = at2*dcu(2,1,k,l)
      exyz(3,1,k,l) = at2*dcu(3,1,k,l)
      exyz(1,1,k1,l) = zero
      exyz(2,1,k1,l) = zero
      exyz(3,1,k1,l) = zero
      exyz(1,1,k,l1) = at2*dcu(1,1,k,l1)
      exyz(2,1,k,l1) = at2*dcu(2,1,k,l1)
      exyz(3,1,k,l1) = at2*dcu(3,1,k,l1)
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
      wp = wp + at1*(dcu(1,1,k,l)*conjg(dcu(1,1,k,l))                   
     1   + dcu(2,1,k,l)*conjg(dcu(2,1,k,l))                             
     2   + dcu(3,1,k,l)*conjg(dcu(3,1,k,l))                             
     3   + dcu(1,1,k,l1)*conjg(dcu(1,1,k,l1))                           
     4   + dcu(2,1,k,l1)*conjg(dcu(2,1,k,l1))                           
     5   + dcu(3,1,k,l1)*conjg(dcu(3,1,k,l1)))
  170 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 180 j = 2, nxh
      at2 = -ci2*real(ffe(j,1,l))
      at1 = at2*at2
      exyz(1,j,1,l) = at2*dcu(1,j,1,l)
      exyz(2,j,1,l) = at2*dcu(2,j,1,l)
      exyz(3,j,1,l) = at2*dcu(3,j,1,l)
      exyz(1,j,k1,l) = zero
      exyz(2,j,k1,l) = zero
      exyz(3,j,k1,l) = zero
      exyz(1,j,1,l1) = at2*dcu(1,j,1,l1)
      exyz(2,j,1,l1) = at2*dcu(2,j,1,l1)
      exyz(3,j,1,l1) = at2*dcu(3,j,1,l1)
      exyz(1,j,k1,l1) = zero
      exyz(2,j,k1,l1) = zero
      exyz(3,j,k1,l1) = zero
      wp = wp + at1*(dcu(1,j,1,l)*conjg(dcu(1,j,1,l))                   
     1   + dcu(2,j,1,l)*conjg(dcu(2,j,1,l))                             
     2   + dcu(3,j,1,l)*conjg(dcu(3,j,1,l))                             
     3   + dcu(1,j,1,l1)*conjg(dcu(1,j,1,l1))                           
     4   + dcu(2,j,1,l1)*conjg(dcu(2,j,1,l1))                           
     5   + dcu(3,j,1,l1)*conjg(dcu(3,j,1,l1)))
  180 continue
c mode numbers kx = 0, nx/2
      at2 = -ci2*real(ffe(1,1,l))
      at1 = at2*at2
      exyz(1,1,1,l) = at2*dcu(1,1,1,l)
      exyz(2,1,1,l) = at2*dcu(2,1,1,l)
      exyz(3,1,1,l) = at2*dcu(3,1,1,l)
      exyz(1,1,k1,l) = zero
      exyz(2,1,k1,l) = zero
      exyz(3,1,k1,l) = zero
      exyz(1,1,1,l1) = zero
      exyz(2,1,1,l1) = zero
      exyz(3,1,1,l1) = zero
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
      wp = wp + at1*(dcu(1,1,1,l)*conjg(dcu(1,1,1,l))                   
     1   + dcu(2,1,1,l)*conjg(dcu(2,1,1,l))                             
     2   + dcu(3,1,1,l)*conjg(dcu(3,1,1,l)))
  190 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 210 k = 2, nyh
      k1 = ny2 - k
      do 200 j = 2, nxh
      at2 = -ci2*real(ffe(j,k,1))
      at1 = at2*at2
      exyz(1,j,k,1) = at2*dcu(1,j,k,1)
      exyz(2,j,k,1) = at2*dcu(2,j,k,1)
      exyz(3,j,k,1) = at2*dcu(3,j,k,1)
      exyz(1,j,k1,1) = at2*dcu(1,j,k1,1)
      exyz(2,j,k1,1) = at2*dcu(2,j,k1,1)
      exyz(3,j,k1,1) = at2*dcu(3,j,k1,1)
      exyz(1,j,k,l1) = zero
      exyz(2,j,k,l1) = zero
      exyz(3,j,k,l1) = zero
      exyz(1,j,k1,l1) = zero
      exyz(2,j,k1,l1) = zero
      exyz(3,j,k1,l1) = zero
      wp = wp + at1*(dcu(1,j,k,1)*conjg(dcu(1,j,k,1))                   
     1   + dcu(2,j,k,1)*conjg(dcu(2,j,k,1))                             
     2   + dcu(3,j,k,1)*conjg(dcu(3,j,k,1))                             
     3   + dcu(1,j,k1,1)*conjg(dcu(1,j,k1,1))                           
     4   + dcu(2,j,k1,1)*conjg(dcu(2,j,k1,1))                           
     5   + dcu(3,j,k1,1)*conjg(dcu(3,j,k1,1)))
  200 continue
  210 continue
c mode numbers kx = 0, nx/2
      do 220 k = 2, nyh
      k1 = ny2 - k
      at2 = -ci2*real(ffe(1,k,1))
      at1 = at2*at2
      exyz(1,1,k,1) = at2*dcu(1,1,k,1)
      exyz(2,1,k,1) = at2*dcu(2,1,k,1)
      exyz(3,1,k,1) = at2*dcu(3,1,k,1)
      exyz(1,1,k1,1) = zero
      exyz(2,1,k1,1) = zero
      exyz(3,1,k1,1) = zero
      exyz(1,1,k,l1) = zero
      exyz(2,1,k,l1) = zero
      exyz(3,1,k,l1) = zero
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
      wp = wp + at1*(dcu(1,1,k,1)*conjg(dcu(1,1,k,1))                   
     1   + dcu(2,1,k,1)*conjg(dcu(2,1,k,1))                             
     2   + dcu(3,1,k,1)*conjg(dcu(3,1,k,1)))
  220 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 230 j = 2, nxh
      at2 = -ci2*real(ffe(j,1,1))
      at1 = at2*at2
      exyz(1,j,1,1) = at2*dcu(1,j,1,1)
      exyz(2,j,1,1) = at2*dcu(2,j,1,1)
      exyz(3,j,1,1) = at2*dcu(3,j,1,1)
      exyz(1,j,k1,1) = zero
      exyz(2,j,k1,1) = zero
      exyz(3,j,k1,1) = zero
      exyz(1,j,1,l1) = zero
      exyz(2,j,1,l1) = zero
      exyz(3,j,1,l1) = zero
      exyz(1,j,k1,l1) = zero
      exyz(2,j,k1,l1) = zero
      exyz(3,j,k1,l1) = zero
      wp = wp + at1*(dcu(1,j,1,1)*conjg(dcu(1,j,1,1))                   
     1   + dcu(2,j,1,1)*conjg(dcu(2,j,1,1))                             
     2   + dcu(3,j,1,1)*conjg(dcu(3,j,1,1)))
  230 continue
      exyz(1,1,1,1) = zero
      exyz(2,1,1,1) = zero
      exyz(3,1,1,1) = zero
      exyz(1,1,k1,1) = zero
      exyz(2,1,k1,1) = zero
      exyz(3,1,k1,1) = zero
      exyz(1,1,1,l1) = zero
      exyz(2,1,1,l1) = zero
      exyz(3,1,1,l1) = zero
      exyz(1,1,k1,l1) = zero
      exyz(2,1,k1,l1) = zero
      exyz(3,1,k1,l1) = zero
      wf = real(nx*ny*nz)*wp/real(ffe(1,1,1))
      return
      end
c-----------------------------------------------------------------------
      subroutine ADDVRFIELD3(a,b,c,ndim,nxe,nye,nze)
c this subroutine calculates a = b + c for real vector fields
      implicit none
      integer ndim, nxe, nye, nze
      real a, b, c
      dimension a(ndim,nxe,nye,nze), b(ndim,nxe,nye,nze)
      dimension c(ndim,nxe,nye,nze)
c local data
      integer i, j, k, l
      do 40 l = 1, nze
      do 30 k = 1, nye
      do 20 j = 1, nxe
      do 10 i = 1, ndim
      a(i,j,k,l) = b(i,j,k,l) + c(i,j,k,l)
   10 continue
   20 continue
   30 continue
   40 continue
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
      subroutine WFFT3RX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd, 
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
         call FFT3RXY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
c perform z fft
         call FFT3RXZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform z fft
         call FFT3RXZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
c perform xy fft
         call FFT3RXY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT3R3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd, 
     1nxhyzd,nxyzhd)
c wrapper function for 3 2d real to complex ffts
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, indz, nxhd, nyd, nzd, nxhyzd, nxyzhd
      dimension f(3,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer ny, nz, nyi, nzi
      data nyi, nzi /1,1/
c calculate range of indices
      ny = 2**indy
      nz = 2**indz
c inverse fourier transform
      if (isign.lt.0) then
c perform xy fft
         call FFT3R3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
     1nzd,nxhyzd,nxyzhd)
c perform z fft
         call FFT3R3Z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform z fft
         call FFT3R3Z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
c perform xy fft
         call FFT3R3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
     1nzd,nxhyzd,nxyzhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT3RN(f,ss,isign,mixup,sct,indx,indy,indz,nxhd,nyd,  
     1nzd,ndim,nxhyzd,nxyzhd)
c wrapper function for multiple 3d real to complex ffts
      implicit none
      complex f, ss, sct
      integer mixup
      integer isign, indx, indy, indz, nxhd, nyd, nzd, ndim
      integer nxhyzd, nxyzhd
      dimension f(ndim,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
      dimension ss(ndim,nxhd)
c local data
      integer ny, nz, nyi, nzi
      data nyi, nzi /1,1/
c calculate range of indices
      ny = 2**indy
      nz = 2**indz
c inverse fourier transform
      if (isign.lt.0) then
c perform xy fft
         call FFT3RNXY(f,ss,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd, 
     1nyd,nzd,ndim,nxhyzd,nxyzhd)
c perform z fft
         call FFT3RNZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd, 
     1nzd,ndim,nxhyzd,nxyzhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform z fft
         call FFT3RNZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd, 
     1nzd,ndim,nxhyzd,nxyzhd)
c perform xy fft
         call FFT3RNXY(f,ss,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd, 
     1nyd,nzd,ndim,nxhyzd,nxyzhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,nxhd, 
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
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = 1, ny
      t2 = t1*f(j2,i,n)
      f(j2,i,n) = f(j1,i,n) - t2
      f(j1,i,n) = f(j1,i,n) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = 1, ny
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
      do 270 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 260 k = 1, ny
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
      do 320 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 310 i = 1, ny
      t2 = t1*f(j2,i,n)
      f(j2,i,n) = f(j1,i,n) - t2
      f(j1,i,n) = f(j1,i,n) + t2
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
      subroutine FFT3R3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,nxhd,
     1nyd,nzd,nxhyzd,nxyzhd)
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
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
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
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(real(f(2,j,i,n)),aimag(f(3,j,i,n)))
      at2 = aimag(f(2,j,i,n))
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
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = 1, ny
      t2 = t1*f(1,j2,i,n)
      t3 = t1*f(2,j2,i,n)
      t4 = t1*f(3,j2,i,n)
      f(1,j2,i,n) = f(1,j1,i,n) - t2
      f(2,j2,i,n) = f(2,j1,i,n) - t3
      f(3,j2,i,n) = f(3,j1,i,n) - t4
      f(1,j1,i,n) = f(1,j1,i,n) + t2
      f(2,j1,i,n) = f(2,j1,i,n) + t3
      f(3,j1,i,n) = f(3,j1,i,n) + t4
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = 1, ny
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
      do 340 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 330 k = 1, ny
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
      do 400 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 390 i = 1, ny
      t2 = t1*f(1,j2,i,n)
      t3 = t1*f(2,j2,i,n)
      t4 = t1*f(3,j2,i,n)
      f(1,j2,i,n) = f(1,j1,i,n) - t2
      f(2,j2,i,n) = f(2,j1,i,n) - t3
      f(3,j2,i,n) = f(3,j1,i,n) - t4
      f(1,j1,i,n) = f(1,j1,i,n) + t2
      f(2,j1,i,n) = f(2,j1,i,n) + t3
      f(3,j1,i,n) = f(3,j1,i,n) + t4
  390 continue
  400 continue
  410 continue
  420 continue
c swap complex components
      do 440 i = 1, ny
      do 430 j = 1, nxh
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(aimag(f(2,j,i,n)),aimag(f(3,j,i,n)))
      at2 = real(f(2,j,i,n))
      f(2,j,i,n) = cmplx(at1,aimag(f(1,j,i,n)))
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
  430 continue
  440 continue
  450 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3R3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,nxhd, 
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
      dimension f(3,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
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
      subroutine FFT3RNXY(f,ss,isign,mixup,sct,indx,indy,indz,nzi,nzp,  
     1nxhd,nyd,nzd,ndim,nxhyzd,nxyzhd)
c this subroutine performs the x-y part of N three dimensional complex
c to real fast fourier transforms and their inverses, for a subset of z,
c using complex arithmetic, where N = ndim
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)
c for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)
c where M = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, three inverse fourier transforms is performed
c f(1:N,n,m,i) = (1/nx*ny*nz)*sum(f(1:N,j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
c       *exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, three forward fourier transforms are performed
c f(1:N,j,k,l) = sum(f(1:N,n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c ss = scratch array
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nzi = initial z index used
c nzp = number of z indices used
c nxhd = second dimension of f
c nyd,nzd = third and fourth dimensions of f
c ndim = leading dimension of array f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(1:N,j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1:N,1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1:N,1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1:N,1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1:N,1,1,1)) = real part of mode nx/2,0,0
c imag(f(1:N,1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1:N,1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1:N,1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
      integer ndim, nxhyzd, nxyzhd
      complex f, ss, sct
      dimension ss(ndim,nxhd)
      dimension f(ndim,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nxyz, nxhyz, nzt, nrx, nry, ns, ns2, km, kmr
      integer i, j, k, l, n, k1, k2, j1, j2, jj
      complex t1, t2, t3
      real ani
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
      if (isign.gt.0) go to 250
c inverse fourier transform
c swap complex components
      call SWAP3CN(f,ss,isign,nxh,ny,nzi,nzt,nxhd,nyd,nzd,ndim)
      do 240 n = nzi, nzt
c bit-reverse array elements in x
      nrx = nxhyz/nxh
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 30
      do 20 i = 1, ny
      do 10 jj = 1, ndim
      t1 = f(jj,j1,i,n)
      f(jj,j1,i,n) = f(jj,j,i,n)
      f(jj,j,i,n) = t1
   10 continue
   20 continue
   30 continue
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
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = 1, ny
      do 40 jj = 1, ndim
      t2 = t1*f(jj,j2,i,n)
      f(jj,j2,i,n) = f(jj,j1,i,n) - t2
      f(jj,j1,i,n) = f(jj,j1,i,n) + t2
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = 1, ny
      do 90 jj = 1, ndim
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
      do 120 jj = 1, ndim
      f(jj,nxhh+1,k,n) = ani*conjg(f(jj,nxhh+1,k,n))
      f(jj,1,k,n) = ani*cmplx(real(f(jj,1,k,n)) + aimag(f(jj,1,k,n)),
     1                        real(f(jj,1,k,n)) - aimag(f(jj,1,k,n)))
  120 continue
  130 continue
c bit-reverse array elements in y
      nry = nxhyz/ny
      do 160 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 160
      do 150 i = 1, nxh
      do 140 jj = 1, ndim
      t1 = f(jj,i,k1,n)
      f(jj,i,k1,n) = f(jj,i,k,n)
      f(jj,i,k,n) = t1
  140 continue
  150 continue
  160 continue
c then transform in y
      nry = nxyz/ny
      do 210 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 200 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 190 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 180 i = 1, nxh
      do 170 jj = 1, ndim
      t2 = t1*f(jj,i,j2,n)
      f(jj,i,j2,n) = f(jj,i,j1,n) - t2
      f(jj,i,j1,n) = f(jj,i,j1,n) + t2
  170 continue
  180 continue
  190 continue
  200 continue
  210 continue
c unscramble modes kx = 0, nx/2
cdir$ ivdep
      do 230 k = 2, nyh
      do 220 jj = 1, ndim
      t1 = f(jj,1,ny2-k,n)
      f(jj,1,ny2-k,n) = 0.5*cmplx(aimag(f(jj,1,k,n) + t1),
     1                            real(f(jj,1,k,n) - t1))
      f(jj,1,k,n) = 0.5*cmplx(real(f(jj,1,k,n) + t1),
     1                        aimag(f(jj,1,k,n) - t1))
  220 continue
  230 continue
  240 continue
      return
c forward fourier transform
  250 do 490 n = nzi, nzt
c scramble modes kx = 0, nx/2
cdir$ ivdep
      do 270 k = 2, nyh
      do 260 jj = 1, ndim
      t1 = cmplx(aimag(f(jj,1,ny2-k,n)),real(f(jj,1,ny2-k,n)))
      f(jj,1,ny2-k,n) = conjg(f(jj,1,k,n) - t1)
      f(jj,1,k,n) = f(jj,1,k,n) + t1
  260 continue
  270 continue
c bit-reverse array elements in y
      nry = nxhyz/ny
      do 300 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 300
      do 290 i = 1, nxh
      do 280 jj = 1, ndim
      t1 = f(jj,i,k1,n)
      f(jj,i,k1,n) = f(jj,i,k,n)
      f(jj,i,k,n) = t1
  280 continue
  290 continue
  300 continue
c then transform in y
      nry = nxyz/ny
      do 350 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 340 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 330 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 320 i = 1, nxh
      do 310 jj = 1, ndim
      t2 = t1*f(jj,i,j2,n)
      f(jj,i,j2,n) = f(jj,i,j1,n) - t2
      f(jj,i,j1,n) = f(jj,i,j1,n) + t2
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
c scramble coefficients
      kmr = nxyz/nx
      do 380 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 370 k = 1, ny
      do 360 jj = 1, ndim
      t2 = conjg(f(jj,nxh2-j,k,n))
      t1 = f(jj,j,k,n) + t2
      t2 = (f(jj,j,k,n) - t2)*t3
      f(jj,j,k,n) = t1 + t2
      f(jj,nxh2-j,k,n) = conjg(t1 - t2)
  360 continue
  370 continue
  380 continue
      do 400 k = 1, ny
      do 390 jj = 1, ndim
      f(jj,nxhh+1,k,n) = 2.0*conjg(f(jj,nxhh+1,k,n))
      f(jj,1,k,n) = cmplx(real(f(jj,1,k,n)) + aimag(f(jj,1,k,n)),
     1                    real(f(jj,1,k,n)) - aimag(f(jj,1,k,n)))
  390 continue
  400 continue
c bit-reverse array elements in x
      nrx = nxhyz/nxh
      do 430 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 430
      do 420 i = 1, ny
      do 410 jj = 1, ndim
      t1 = f(jj,j1,i,n)
      f(jj,j1,i,n) = f(jj,j,i,n)
      f(jj,j,i,n) = t1
  410 continue
  420 continue
  430 continue
c finally transform in x
      nrx = nxyz/nxh
      do 480 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 470 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 460 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 450 i = 1, ny
      do 440 jj = 1, ndim
      t2 = t1*f(jj,j2,i,n)
      f(jj,j2,i,n) = f(jj,j1,i,n) - t2
      f(jj,j1,i,n) = f(jj,j1,i,n) + t2
  440 continue
  450 continue
  460 continue
  470 continue
  480 continue
  490 continue
c swap complex components
      call SWAP3CN(f,ss,isign,nxh,ny,nzi,nzt,nxhd,nyd,nzd,ndim)
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RNZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,nxhd, 
     1nyd,nzd,ndim,nxhyzd,nxyzhd)
c this subroutine performs the z part of N three dimensional complex to
c real fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic, where N = ndim
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)
c for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)
c where M = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, three inverse fourier transforms is performed
c f(1:N,j,k,l) = sum(f(1:N,j,k,i)*exp(-sqrt(-1)*2pi*l*i/nz))
c if isign = 1, three forward fourier transforms are performed
c f(1:N,n,m,i) = sum(f(1:N,n,m,l)*exp(sqrt(-1)*2pi*l*i/nz))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = second dimension of f
c nyd,nzd = third and fourth dimensions of f
c ndim = leading dimension of array f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(1:N,j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1:N,1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1:N,1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1:N,1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1:N,1,1,1)) = real part of mode nx/2,0,0
c imag(f(1:N,1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1:N,1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1:N,1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, indz, nyi, nyp, nxhd, nyd
      integer nzd, ndim, nxhyzd, nxyzhd
      complex f, sct
      dimension f(ndim,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh
      integer nz, nzh, nz2, nxyz, nxhyz, nyt, nrz, ns, ns2, km, kmr
      integer i, j, k, l, n, k1, k2, j1, j2, l1, jj
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
      if (isign.gt.0) go to 140
c inverse fourier transform
c bit-reverse array elements in z
      nrz = nxhyz/nz
      do 40 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 40
      do 30 n = nyi, nyt
      do 20 i = 1, nxh
      do 10 jj = 1, ndim
      t1 = f(jj,i,n,l1)
      f(jj,i,n,l1) = f(jj,i,n,l)
      f(jj,i,n,l) = t1
   10 continue
   20 continue
   30 continue
   40 continue
c finally transform in z
      nrz = nxyz/nz
      do 100 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 70 n = nyi, nyt
      do 60 i = 1, nxh
      do 50 jj = 1, ndim
      t2 = t1*f(jj,i,n,j2)
      f(jj,i,n,j2) = f(jj,i,n,j1) - t2
      f(jj,i,n,j1) = f(jj,i,n,j1) + t2
   50 continue
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble modes kx = 0, nx/2
cdir$ ivdep
      do 130 n = 2, nzh
      if (nyi.eq.1) then
         do 110 jj = 1, ndim
         t1 = f(jj,1,1,nz2-n)
         f(jj,1,1,nz2-n) = 0.5*cmplx(aimag(f(jj,1,1,n) + t1),
     1                               real(f(jj,1,1,n) - t1))
         f(jj,1,1,n) = 0.5*cmplx(real(f(jj,1,1,n) + t1),
     1                           aimag(f(jj,1,1,n) - t1))
  110    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         do 120 jj = 1, ndim
         t1 = f(jj,1,nyh+1,nz2-n)
         f(jj,1,nyh+1,nz2-n) = 0.5*cmplx(aimag(f(jj,1,nyh+1,n) + t1),
     1                                   real(f(jj,1,nyh+1,n) - t1))
         f(jj,1,nyh+1,n) = 0.5*cmplx(real(f(jj,1,nyh+1,n) + t1),
     1                               aimag(f(jj,1,nyh+1,n) - t1))
  120    continue
      endif
  130 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
cdir$ ivdep
  140 do 170 n = 2, nzh
      if (nyi.eq.1) then
         do 150 jj = 1, ndim
         t1 = cmplx(aimag(f(jj,1,1,nz2-n)),real(f(jj,1,1,nz2-n)))
         f(jj,1,1,nz2-n) = conjg(f(jj,1,1,n) - t1)
         f(jj,1,1,n) = f(jj,1,1,n) + t1
  150    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         do 160 jj = 1, ndim
         t1 = cmplx(aimag(f(jj,1,nyh+1,nz2-n)),
     1              real(f(jj,1,nyh+1,nz2-n)))
         f(jj,1,nyh+1,nz2-n) = conjg(f(jj,1,nyh+1,n) - t1)
         f(jj,1,nyh+1,n) = f(jj,1,nyh+1,n) + t1
  160    continue
      endif
  170 continue
c bit-reverse array elements in z
      nrz = nxhyz/nz
      do 210 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 210
      do 200 n = nyi, nyt
      do 190 i = 1, nxh
      do 180 jj = 1, ndim
      t1 = f(jj,i,n,l1)
      f(jj,i,n,l1) = f(jj,i,n,l)
      f(jj,i,n,l) = t1
  180 continue
  190 continue
  200 continue
  210 continue
c first transform in z
      nrz = nxyz/nz
      do 270 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 260 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 250 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 240 n = nyi, nyt
      do 230 i = 1, nxh
      do 220 jj = 1, ndim
      t2 = t1*f(jj,i,n,j2)
      f(jj,i,n,j2) = f(jj,i,n,j1) - t2
      f(jj,i,n,j1) = f(jj,i,n,j1) + t2
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine SWAP3CN(f,s,isign,nxh,ny,nzi,nzt,nxhd,nyd,nzd,ndim)
c this subroutine swaps components for multiple ffts
c f = input  array
c s = scratch array
c isign = (-1,1) = swap (real-to-complex,complex-to-real)
c nxh = complex dimension in x direction
c ny = complex dimension in y direction
c nzi/nzt = initial/final z index used
c nxhd = half of the second dimension of f
c nyd,nzd = third and fourth dimension of f
c ndim = leading dimension of array f
      implicit none
      integer isign, nxh, ny, nzi, nzt, nxhd, nyd, nzd, ndim
      real f, s
      dimension f(ndim,2*nxhd,nyd,nzd), s(2*ndim*nxhd)
c local data
      integer i, j, k, l, ioff
c swap complex components
c real to complex
      if (isign.lt.0) then
         do 70 l = nzi, nzt
         do 60 k = 1, ny
         do 20 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 10 i = 1, ndim
         s(2*i+ioff-1) = f(i,2*j-1,k,l)
         s(2*i+ioff) = f(i,2*j,k,l)
   10    continue
   20    continue
         do 50 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 30 i = 1, ndim
         f(i,2*j-1,k,l) = s(i+ioff)
   30    continue
         ioff = ioff + ndim
         do 40 i = 1, ndim
         f(i,2*j,k,l) = s(i+ioff)
   40    continue
   50    continue
   60    continue
   70    continue
c complex to real
      else if (isign.gt.0) then
         do 140 l = nzi, nzt
         do 130 k = 1, ny
         do 100 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 80 i = 1, ndim
         s(i+ioff) = f(i,2*j-1,k,l)
   80    continue
         ioff = ioff + ndim
         do 90 i = 1, ndim
         s(i+ioff) = f(i,2*j,k,l)
   90    continue
  100    continue
         do 120 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 110 i = 1, ndim
         f(i,2*j-1,k,l) = s(2*i+ioff-1)
         f(i,2*j,k,l) = s(2*i+ioff)
  110    continue
  120    continue
  130    continue
  140    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine GSJPOST3L(part,cu,qm,dt,nop,idimp,nx,ny,nz,nxv,nyv,    
     1nxyzv,ipbc)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c part(4,n) = x velocity of particle n
c part(5,n) = y velocity of particle n
c part(6,n) = z velocity of particle n
c cu(i,n) = ith component of current density at grid point j,k,l
c where n = j + nxv*(k-1) + nxv*nyv*(l-1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 6
c nx/ny/nz = system length in x/y/z direction
c nxv = second virtual dimension of current array, must be >= nx+1
c nyv = third virtual dimension of current array, must be >= ny+1
c nxyzv = actual dimension of current array, must be >= nxv*nyv*(nz+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer nop,idimp,nx,ny,nz,nxv,nyv,nxyzv,ipbc
      real qm, dt
      real part, cu
      dimension part(idimp,nop), cu(3,nxyzv)
c local data
      integer j, nnn, mmm, lll, nn, mm, ll, mp, lp, nxyv
      real dxn, dyn, dzn, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, vx, vy, vz
      real dx1, dx2, dx3
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmm = part(2,1)
      lll = part(3,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmm)
      dzn = part(3,1) - real(lll)
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
      nxyv = nxv*nyv
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      mm = nxv*mmm + nxyv*lll
      nnn = part(1,j)
      mmm = part(2,j)
      lll = part(3,j)
      dxp = qm*dxn
      dyp = dyn
      dzp = dzn
      dxn = part(1,j) - real(nnn)
      dyn = part(2,j) - real(mmm)
      dzn = part(3,j) - real(lll)
      amx = qm - dxp
      amy = 1.0 - dyp
      mm = mm + nn
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + nxv
      amx = amx*amy
      amz = 1.0 - dzp
      ll = mm + nxyv
      amy = dxp*amy
      lp = mp + nxyv
c deposit current
      dx = amx*amz
      dz = amy*amz
      vx = part(4,j-1)
      vy = part(5,j-1)
      vz = part(6,j-1)
      dxp = cu(1,mm) + vx*dx
      dx2 = cu(2,mm) + vy*dx
      dx3 = cu(3,mm) + vz*dx
      dx = cu(1,mm+1) + vx*dz
      dy = cu(2,mm+1) + vy*dz
      dz = cu(3,mm+1) + vz*dz
      cu(1,mm) = dxp
      cu(2,mm) = dx2
      cu(3,mm) = dx3
      cu(1,mm+1) = dx
      cu(2,mm+1) = dy
      cu(3,mm+1) = dz
      dx = dyp*amz
      dz = dx1*amz
      dxp = cu(1,mp) + vx*dx
      dx2 = cu(2,mp) + vy*dx
      dx3 = cu(3,mp) + vz*dx
      dx = cu(1,mp+1) + vx*dz
      dy = cu(2,mp+1) + vy*dz
      dz = cu(3,mp+1) + vz*dz
      cu(1,mp) = dxp
      cu(2,mp) = dx2
      cu(3,mp) = dx3
      cu(1,mp+1) = dx
      cu(2,mp+1) = dy
      cu(3,mp+1) = dz
      dx = amx*dzp
      dz = amy*dzp
      dxp = cu(1,ll) + vx*dx
      dx2 = cu(2,ll) + vy*dx
      dx3 = cu(3,ll) + vz*dx
      dx = cu(1,ll+1) + vx*dz
      dy = cu(2,ll+1) + vy*dz
      dz = cu(3,ll+1) + vz*dz
      cu(1,ll) = dxp
      cu(2,ll) = dx2
      cu(3,ll) = dx3
      cu(1,ll+1) = dx
      cu(2,ll+1) = dy
      cu(3,ll+1) = dz
      dx = dyp*dzp
      dz = dx1*dzp
      dxp = cu(1,lp) + vx*dx
      dx2 = cu(2,lp) + vy*dx
      dx3 = cu(3,lp) + vz*dx
      dx = cu(1,lp+1) + vx*dz
      dy = cu(2,lp+1) + vy*dz
      dz = cu(3,lp+1) + vz*dz
      cu(1,lp) = dxp
      cu(2,lp) = dx2
      cu(3,lp) = dx3
      cu(1,lp+1) = dx
      cu(2,lp+1) = dy
      cu(3,lp+1) = dz
c advance position half a time-step
      dx = part(1,j-1) + vx*dt
      dy = part(2,j-1) + vy*dt
      dz = part(3,j-1) + vz*dt
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
            dx = part(1,j-1)
            part(4,j-1) = -part(4,j-1)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1)
            part(5,j-1) = -part(5,j-1)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j-1)
            part(6,j-1) = -part(6,j-1)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(4,j-1) = -part(4,j-1)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1)
            part(5,j-1) = -part(5,j-1)
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(1,j-1) = dx
      part(2,j-1) = dy
      part(3,j-1) = dz
   10 continue
c deposit current for last particle
      nn = nnn + 1
      mm = nxv*mmm + nxyv*lll
      dxp = qm*dxn
      amx = qm - dxp
      amy = 1.0 - dyn
      mm = mm + nn
      dx1 = dxp*dyn
      dyp = amx*dyn
      mp = mm + nxv
      amx = amx*amy
      amz = 1.0 - dzn
      ll = mm + nxyv
      amy = dxp*amy
      lp = mp + nxyv
c deposit current
      dx = amx*amz
      dy = amy*amz
      vx = part(4,nop)
      vy = part(5,nop)
      vz = part(6,nop)
      cu(1,mm) = cu(1,mm) + vx*dx
      cu(2,mm) = cu(2,mm) + vy*dx
      cu(3,mm) = cu(3,mm) + vz*dx
      cu(1,mm+1) = cu(1,mm+1) + vx*dy
      cu(2,mm+1) = cu(2,mm+1) + vy*dy
      cu(3,mm+1) = cu(3,mm+1) + vz*dy
      dx = dyp*amz
      dy = dx1*amz
      cu(1,mp) = cu(1,mp) + vx*dx
      cu(2,mp) = cu(2,mp) + vy*dx
      cu(3,mp) = cu(3,mp) + vz*dx
      cu(1,mp+1) = cu(1,mp+1) + vx*dy
      cu(2,mp+1) = cu(2,mp+1) + vy*dy
      cu(3,mp+1) = cu(3,mp+1) + vz*dy
      dx = amx*dzn
      dy = amy*dzn
      cu(1,ll) = cu(1,ll) + vx*dx
      cu(2,ll) = cu(2,ll) + vy*dx
      cu(3,ll) = cu(3,ll) + vz*dx
      cu(1,ll+1) = cu(1,ll+1) + vx*dy
      cu(2,ll+1) = cu(2,ll+1) + vy*dy
      cu(3,ll+1) = cu(3,ll+1) + vz*dy
      dx = dyp*dzn
      dy = dx1*dzn
      cu(1,lp) = cu(1,lp) + vx*dx
      cu(2,lp) = cu(2,lp) + vy*dx
      cu(3,lp) = cu(3,lp) + vz*dx
      cu(1,lp+1) = cu(1,lp+1) + vx*dy
      cu(2,lp+1) = cu(2,lp+1) + vy*dy
      cu(3,lp+1) = cu(3,lp+1) + vz*dy
c advance position half a time-step
      dx = part(1,nop) + vx*dt
      dy = part(2,nop) + vy*dt
      dz = part(3,nop) + vz*dt
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
            dx = part(1,nop)
            part(4,nop) = -part(4,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(5,nop) = -part(5,nop)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,nop)
            part(6,nop) = -part(6,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(4,nop) = -part(4,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(5,nop) = -part(5,nop)
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
      part(3,nop) = dz
      return
      end
c-----------------------------------------------------------------------
      subroutine GSPOST3L(part,q,qm,nop,idimp,nxv,nyv,nxyzv)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c q(j,k,l) = charge density at grid point j,k,l
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 6
c nxv = first virtual dimension of charge array, must be >= nx+1
c nyv = second virtual dimension of charge array, must be >= ny+1
c nxyzv = actual dimension of charge array, must be >= nxv*nyv*(nz+1)
      implicit none
      integer nop, idimp, nxv, nyv, nxyzv
      real qm
      real part, q
      dimension part(idimp,nop), q(nxyzv)
c local data
      integer j, nnn, mmm, lll, nxyv, nn, mm, ll, mp, lp
      real dxn, dyn, dzn, dxp, dyp, dzp, amx, amy, amz, dx1, dx2, dx3
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmm = part(2,1)
      lll = part(3,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmm)
      dzn = part(3,1) - real(lll)
      nxyv = nxv*nyv
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      mm = nxv*mmm + nxyv*lll
      nnn = part(1,j)
      mmm = part(2,j)
      lll = part(3,j)
      dxp = qm*dxn
      dyp = dyn
      dzp = dzn
      dxn = part(1,j) - real(nnn)
      dyn = part(2,j) - real(mmm)
      dzn = part(3,j) - real(lll)
      amx = qm - dxp
      amy = 1.0 - dyp
      mm = mm + nn
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + nxv
      amx = amx*amy
      amz = 1.0 - dzp
      ll = mm + nxyv
      amy = dxp*amy
      lp = mp + nxyv
c deposit charge
      dxp = q(mm) + amx*amz
      dx2 = q(mm+1) + amy*amz
      dx3 = q(mp) + dyp*amz
      amz = q(mp+1) + dx1*amz
      amx = q(ll) + amx*dzp
      amy = q(ll+1) + amy*dzp
      dyp = q(lp) + dyp*dzp
      dzp = q(lp+1) + dx1*dzp
      q(mm) = dxp
      q(mm+1) = dx2
      q(mp) = dx3
      q(mp+1) = amz
      q(ll) = amx
      q(ll+1) = amy
      q(lp) = dyp
      q(lp+1) = dzp
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      mm = nxv*mmm + nxyv*lll
      dxp = qm*dxn
      amx = qm - dxp
      amy = 1.0 - dyn
      mm = mm + nn
      dx1 = dxp*dyn
      dyp = amx*dyn
      mp = mm + nxv
      amx = amx*amy
      amz = 1.0 - dzn
      ll = mm + nxyv
      amy = dxp*amy
      lp = mp + nxyv
      q(mm) = q(mm) + amx*amz
      q(mm+1) = q(mm+1) + amy*amz
      q(mp) = q(mp) + dyp*amz
      q(mp+1) = q(mp+1) + dx1*amz
      q(ll) = q(ll) + amx*dzn
      q(ll+1) = q(ll+1) + amy*dzn
      q(lp) = q(lp) + dyp*dzn
      q(lp+1) = q(lp+1) + dx1*dzn
      return
      end
c-----------------------------------------------------------------------
      subroutine GSBPUSH3L(part,fxyz,bxyz,qbm,dt,dtc,ek,idimp,nop,nx,ny,
     1nz,nxv,nyv,nxyzv,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field.  Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c part(4,n) = velocity vx of particle n
c part(5,n) = velocity vy of particle n
c part(6,n) = velocity vz of particle n
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
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c nxyzv = actual dimension of field array, must be >= nxv*nyv*(nz+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nop, nx, ny, nz, nxv, nyv, nxyzv, ipbc
      real qbm, dt, dtc, ek
      real part, fxyz, bxyz
      dimension part(idimp,nop)
      dimension fxyz(3,nxyzv), bxyz(3,nxyzv)
c local data
      integer nnn, mmm, lll, nop1, j, nn, mm, ll, mp, lp, nxyv
      real qtmh, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxn, dyn, dzn, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, dx1
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmm = part(2,1)
      lll = part(3,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmm)
      dzn = part(3,1) - real(lll)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
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
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmm + nxyv*lll
      nnn = part(1,j+1)
      mmm = part(2,j+1)
      lll = part(3,j+1)
      dxp = dxn
      dyp = dyn
      dzp = dzn
      dxn = part(1,j+1) - real(nnn)
      dyn = part(2,j+1) - real(mmm)
      dzn = part(3,j+1) - real(lll)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      mm = mm + nn
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + nxv
      amx = amx*amy
      amz = 1.0 - dzp
      ll = mm + nxyv
      amy = dxp*amy
      lp = mp + nxyv
c find electric field
      dx = amz*(amx*fxyz(1,mm) + amy*fxyz(1,mm+1) + dyp*fxyz(1,mp)
     1        + dx1*fxyz(1,mp+1))
     2   + dzp*(amx*fxyz(1,ll) + amy*fxyz(1,ll+1) + dyp*fxyz(1,lp)
     3        + dx1*fxyz(1,lp+1))
      dy = amz*(amx*fxyz(2,mm) + amy*fxyz(2,mm+1) + dyp*fxyz(2,mp)
     1        + dx1*fxyz(2,mp+1))
     2   + dzp*(amx*fxyz(2,ll) + amy*fxyz(2,ll+1) + dyp*fxyz(2,lp)
     3        + dx1*fxyz(2,lp+1))
      dz = amz*(amx*fxyz(3,mm) + amy*fxyz(3,mm+1) + dyp*fxyz(3,mp)
     1        + dx1*fxyz(3,mp+1))
     2   + dzp*(amx*fxyz(3,ll) + amy*fxyz(3,ll+1) + dyp*fxyz(3,lp)
     3        + dx1*fxyz(3,lp+1))
c find magnetic field
      ox = amz*(amx*bxyz(1,mm) + amy*bxyz(1,mm+1) + dyp*bxyz(1,mp)
     1        + dx1*bxyz(1,mp+1))
     2   + dzp*(amx*bxyz(1,ll) + amy*bxyz(1,ll+1) + dyp*bxyz(1,lp)
     3        + dx1*bxyz(1,lp+1))
      oy = amz*(amx*bxyz(2,mm) + amy*bxyz(2,mm+1) + dyp*bxyz(2,mp) 
     1        + dx1*bxyz(2,mp+1))
     2   + dzp*(amx*bxyz(2,ll) + amy*bxyz(2,ll+1) + dyp*bxyz(2,lp)
     3        + dx1*bxyz(2,lp+1))
      oz = amz*(amx*bxyz(3,mm) + amy*bxyz(3,mm+1) + dyp*bxyz(3,mp)
     1        + dx1*bxyz(3,mp+1))
     2   + dzp*(amx*bxyz(3,ll) + amy*bxyz(3,ll+1) + dyp*bxyz(3,lp)
     3        + dx1*bxyz(3,lp+1))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j) + dx
      acy = part(5,j) + dy
      acz = part(6,j) + dz
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
      part(4,j) = dx
      part(5,j) = dy
      part(6,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
      dz = part(3,j) + dz*dtc
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
            dx = part(1,j)
            part(4,j) = -part(4,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(5,j) = -part(5,j)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j)
            part(6,j) = -part(6,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(4,j) = -part(4,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(5,j) = -part(5,j)
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
      part(3,j) = dz
   10 continue
c push last particle
      nn = nnn + 1
      mm = nxv*mmm + nxyv*lll
      amx = 1.0 - dxn
      amy = 1.0 - dyn
      mm = mm + nn
      dx1 = dxn*dyn
      dyp = amx*dyn
      mp = mm + nxv
      amx = amx*amy
      amz = 1.0 - dzn
      ll = mm + nxyv
      amy = dxn*amy
      lp = mp + nxyv
c find electric field
      dx = amz*(amx*fxyz(1,mm) + amy*fxyz(1,mm+1) + dyp*fxyz(1,mp)
     1        + dx1*fxyz(1,mp+1))
     2   + dzn*(amx*fxyz(1,ll) + amy*fxyz(1,ll+1) + dyp*fxyz(1,lp)
     3        + dx1*fxyz(1,lp+1))
      dy = amz*(amx*fxyz(2,mm) + amy*fxyz(2,mm+1) + dyp*fxyz(2,mp)
     1        + dx1*fxyz(2,mp+1))
     2   + dzn*(amx*fxyz(2,ll) + amy*fxyz(2,ll+1) + dyp*fxyz(2,lp)
     3        + dx1*fxyz(2,lp+1))
      dz = amz*(amx*fxyz(3,mm) + amy*fxyz(3,mm+1) + dyp*fxyz(3,mp)
     1        + dx1*fxyz(3,mp+1))
     2   + dzn*(amx*fxyz(3,ll) + amy*fxyz(3,ll+1) + dyp*fxyz(3,lp)
     3        + dx1*fxyz(3,lp+1))
c find magnetic field
      ox = amz*(amx*bxyz(1,mm) + amy*bxyz(1,mm+1) + dyp*bxyz(1,mp)
     1        + dx1*bxyz(1,mp+1))
     2   + dzn*(amx*bxyz(1,ll) + amy*bxyz(1,ll+1) + dyp*bxyz(1,lp)
     3        + dx1*bxyz(1,lp+1))
      oy = amz*(amx*bxyz(2,mm) + amy*bxyz(2,mm+1) + dyp*bxyz(2,mp)
     1        + dx1*bxyz(2,mp+1))
     2   + dzn*(amx*bxyz(2,ll) + amy*bxyz(2,ll+1) + dyp*bxyz(2,lp)
     3        + dx1*bxyz(2,lp+1))
      oz = amz*(amx*bxyz(3,mm) + amy*bxyz(3,mm+1) + dyp*bxyz(3,mp)
     1        + dx1*bxyz(3,mp+1)) 
     2   + dzn*(amx*bxyz(3,ll) + amy*bxyz(3,ll+1) + dyp*bxyz(3,lp)
     3        + dx1*bxyz(3,lp+1))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,nop) + dx
      acy = part(5,nop) + dy
      acz = part(6,nop) + dz
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
      part(4,nop) = dx
      part(5,nop) = dy
      part(6,nop) = dz
c new position
      dx = part(1,nop) + dx*dtc
      dy = part(2,nop) + dy*dtc
      dz = part(3,nop) + dz*dtc
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
            dx = part(1,nop)
            part(4,nop) = -part(4,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(5,nop) = -part(5,nop)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,nop)
            part(6,nop) = -part(6,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(4,nop) = -part(4,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(5,nop) = -part(5,nop)
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
      part(3,nop) = dz
c normalize kinetic energy
      ek = ek + .5*sum1
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
