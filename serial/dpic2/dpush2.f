c Fortran Library for Skeleton 2-1/2D Darwin PIC Code
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
      subroutine GBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,  
     1nxv,nyv,ipbc)
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
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
c nx/ny = system length in x/y direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real qbm, dt, dtc, ek
      real part, fxy, bxy
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
c local data
      integer j, nn, mm, np, mp
      real qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
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
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1.0 - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp))                    
     1   + amy*(dxp*fxy(1,np,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp))                    
     1   + amy*(dxp*fxy(2,np,mm) + amx*fxy(2,nn,mm))
      dz = dyp*(dxp*fxy(3,np,mp) + amx*fxy(3,nn,mp))                    
     1   + amy*(dxp*fxy(3,np,mm) + amx*fxy(3,nn,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp))                    
     1   + amy*(dxp*bxy(1,np,mm) + amx*bxy(1,nn,mm))
      oy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp))                    
     1   + amy*(dxp*bxy(2,np,mm) + amx*bxy(2,nn,mm))
      oz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp))                    
     1   + amy*(dxp*bxy(3,np,mm) + amx*bxy(3,nn,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
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
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c normalize kinetic energy
      ek = ek + 0.5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine GPOST2L(part,q,qm,nop,idimp,nxv,nyv)
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(j,k) = charge density at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
      implicit none
      integer nop, idimp, nxv, nyv
      real qm
      real part, q
      dimension part(idimp,nop), q(nxv,nyv)
c local data
      integer j, nn, mm, np, mp
      real dxp, dyp, amx, amy
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c deposit charge
      q(np,mp) = q(np,mp) + dxp*dyp
      q(nn,mp) = q(nn,mp) + amx*dyp
      q(np,mm) = q(np,mm) + dxp*amy
      q(nn,mm) = q(nn,mm) + amx*amy
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,j,k) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = second dimension of current array, must be >= nx+1
c nyv = third dimension of current array, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxv, nyv, ipbc
      real qm, dt
      real part, cu
      dimension part(idimp,nop), cu(3,nxv,nyv)
c local data
      integer j, nn, mm, np, mp
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, vx, vy, vz
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
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
c deposit current
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      cu(1,np,mp) = cu(1,np,mp) + vx*dx
      cu(2,np,mp) = cu(2,np,mp) + vy*dx
      cu(3,np,mp) = cu(3,np,mp) + vz*dx
      dx = dxp*amy
      cu(1,nn,mp) = cu(1,nn,mp) + vx*dy
      cu(2,nn,mp) = cu(2,nn,mp) + vy*dy
      cu(3,nn,mp) = cu(3,nn,mp) + vz*dy
      dy = amx*amy
      cu(1,np,mm) = cu(1,np,mm) + vx*dx
      cu(2,np,mm) = cu(2,np,mm) + vy*dx
      cu(3,np,mm) = cu(3,np,mm) + vz*dx
      cu(1,nn,mm) = cu(1,nn,mm) + vx*dy
      cu(2,nn,mm) = cu(2,nn,mm) + vy*dy
      cu(3,nn,mm) = cu(3,nn,mm) + vz*dy
c advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GMJPOST2L(part,amu,qm,nop,idimp,nxv,nyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using first-order spline interpolation
c scalar version using guard cells
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
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x velocity of particle n at t - dt/2
c part(4,n) = y velocity of particle n at t - dt/2
c part(5,n) = z velocity of particle n at t - dt/2
c amu(i,j,k) = ith component of momentum flux at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 5
c nxv = second dimension of flux array, must be >= nx+1
c nyv = third dimension of flux array, must be >= ny+1
      implicit none
      integer nop, idimp, nxv, nyv
      real part, amu, qm
      dimension part(idimp,nop), amu(4,nxv,nyv)
c local data
      integer j, nn, mm, np, mp
      real dxp, dyp, amx, amy
      real dx, dy, vx, vy, vz, v1, v2, v3, v4
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
c deposit momentum flux
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      amu(1,np,mp) = amu(1,np,mp) + v1*dx
      amu(2,np,mp) = amu(2,np,mp) + v2*dx
      amu(3,np,mp) = amu(3,np,mp) + v3*dx
      amu(4,np,mp) = amu(4,np,mp) + v4*dx
      dx = dxp*amy
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      amu(3,nn,mp) = amu(3,nn,mp) + v3*dy
      amu(4,nn,mp) = amu(4,nn,mp) + v4*dy
      dy = amx*amy
      amu(1,np,mm) = amu(1,np,mm) + v1*dx
      amu(2,np,mm) = amu(2,np,mm) + v2*dx
      amu(3,np,mm) = amu(3,np,mm) + v3*dx
      amu(4,np,mm) = amu(4,np,mm) + v4*dx
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      amu(3,nn,mm) = amu(3,nn,mm) + v3*dy
      amu(4,nn,mm) = amu(4,nn,mm) + v4*dy
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GDJPOST2L(part,fxy,bxy,dcu,amu,qm,qbm,dt,idimp,nop,nxv,
     1nyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c and acceleration density using first-order spline interpolation.
c scalar version using guard cells
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
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
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
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 5
c nop = number of particles
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, bxy, dcu, amu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension dcu(3,nxv,nyv), amu(4,nxv,nyv)
c local data
      integer j, nn, mm, np, mp
      real qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, v3, v4
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp))                    
     1   + amy*(dxp*fxy(1,np,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp))                    
     1   + amy*(dxp*fxy(2,np,mm) + amx*fxy(2,nn,mm))
      dz = dyp*(dxp*fxy(3,np,mp) + amx*fxy(3,nn,mp))                    
     1   + amy*(dxp*fxy(3,np,mm) + amx*fxy(3,nn,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp))                    
     1   + amy*(dxp*bxy(1,np,mm) + amx*bxy(1,nn,mm))
      oy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp))                    
     1   + amy*(dxp*bxy(2,np,mm) + amx*bxy(2,nn,mm))
      oz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp))                    
     1   + amy*(dxp*bxy(3,np,mm) + amx*bxy(3,nn,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
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
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
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
      dx = dxp*dyp
      dy = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      amu(1,np,mp) = amu(1,np,mp) + v1*dx
      amu(2,np,mp) = amu(2,np,mp) + v2*dx
      amu(3,np,mp) = amu(3,np,mp) + v3*dx
      amu(4,np,mp) = amu(4,np,mp) + v4*dx
      dcu(1,np,mp) = dcu(1,np,mp) + vx*dx
      dcu(2,np,mp) = dcu(2,np,mp) + vy*dx
      dcu(3,np,mp) = dcu(3,np,mp) + vz*dx
      dx = dxp*amy
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      amu(3,nn,mp) = amu(3,nn,mp) + v3*dy
      amu(4,nn,mp) = amu(4,nn,mp) + v4*dy
      dcu(1,nn,mp) = dcu(1,nn,mp) + vx*dy
      dcu(2,nn,mp) = dcu(2,nn,mp) + vy*dy
      dcu(3,nn,mp) = dcu(3,nn,mp) + vz*dy
      dy = amx*amy
      amu(1,np,mm) = amu(1,np,mm) + v1*dx
      amu(2,np,mm) = amu(2,np,mm) + v2*dx
      amu(3,np,mm) = amu(3,np,mm) + v3*dx
      amu(4,np,mm) = amu(4,np,mm) + v4*dx
      dcu(1,np,mm) = dcu(1,np,mm) + vx*dx
      dcu(2,np,mm) = dcu(2,np,mm) + vy*dx
      dcu(3,np,mm) = dcu(3,np,mm) + vz*dx
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      amu(3,nn,mm) = amu(3,nn,mm) + v3*dy
      amu(4,nn,mm) = amu(4,nn,mm) + v4*dy
      dcu(1,nn,mm) = dcu(1,nn,mm) + vx*dy
      dcu(2,nn,mm) = dcu(2,nn,mm) + vy*dy
      dcu(3,nn,mm) = dcu(3,nn,mm) + vz*dy
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,nop,
     1nxv,nyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density and current density using first-order spline
c interpolation.
c scalar version using guard cells
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
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
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
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 5
c nop = number of particles
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension cu(3,nxv,nyv), dcu(3,nxv,nyv), amu(4,nxv,nyv)
c local data
      integer j, nn, mm, np, mp
      real qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, v3, v4
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1.0 - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp))                    
     1   + amy*(dxp*fxy(1,np,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp))                    
     1   + amy*(dxp*fxy(2,np,mm) + amx*fxy(2,nn,mm))
      dz = dyp*(dxp*fxy(3,np,mp) + amx*fxy(3,nn,mp))                    
     1   + amy*(dxp*fxy(3,np,mm) + amx*fxy(3,nn,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp))                    
     1   + amy*(dxp*bxy(1,np,mm) + amx*bxy(1,nn,mm))
      oy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp))                    
     1   + amy*(dxp*bxy(2,np,mm) + amx*bxy(2,nn,mm))
      oz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp))                    
     1   + amy*(dxp*bxy(3,np,mm) + amx*bxy(3,nn,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
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
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
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
      dx = dxp*dyp
      dy = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      amu(1,np,mp) = amu(1,np,mp) + v1*dx
      amu(2,np,mp) = amu(2,np,mp) + v2*dx
      amu(3,np,mp) = amu(3,np,mp) + v3*dx
      amu(4,np,mp) = amu(4,np,mp) + v4*dx
      dcu(1,np,mp) = dcu(1,np,mp) + vx*dx
      dcu(2,np,mp) = dcu(2,np,mp) + vy*dx
      dcu(3,np,mp) = dcu(3,np,mp) + vz*dx
      cu(1,np,mp) = cu(1,np,mp) + ox*dx
      cu(2,np,mp) = cu(2,np,mp) + oy*dx
      cu(3,np,mp) = cu(3,np,mp) + oz*dx
      dx = dxp*amy
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      amu(3,nn,mp) = amu(3,nn,mp) + v3*dy
      amu(4,nn,mp) = amu(4,nn,mp) + v4*dy
      dcu(1,nn,mp) = dcu(1,nn,mp) + vx*dy
      dcu(2,nn,mp) = dcu(2,nn,mp) + vy*dy
      dcu(3,nn,mp) = dcu(3,nn,mp) + vz*dy
      cu(1,nn,mp) = cu(1,nn,mp) + ox*dy
      cu(2,nn,mp) = cu(2,nn,mp) + oy*dy
      cu(3,nn,mp) = cu(3,nn,mp) + oz*dy
      dy = amx*amy
      amu(1,np,mm) = amu(1,np,mm) + v1*dx
      amu(2,np,mm) = amu(2,np,mm) + v2*dx
      amu(3,np,mm) = amu(3,np,mm) + v3*dx
      amu(4,np,mm) = amu(4,np,mm) + v4*dx
      dcu(1,np,mm) = dcu(1,np,mm) + vx*dx
      dcu(2,np,mm) = dcu(2,np,mm) + vy*dx
      dcu(3,np,mm) = dcu(3,np,mm) + vz*dx
      cu(1,np,mm) = cu(1,np,mm) + ox*dx
      cu(2,np,mm) = cu(2,np,mm) + oy*dx
      cu(3,np,mm) = cu(3,np,mm) + oz*dx
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      amu(3,nn,mm) = amu(3,nn,mm) + v3*dy
      amu(4,nn,mm) = amu(4,nn,mm) + v4*dy
      dcu(1,nn,mm) = dcu(1,nn,mm) + vx*dy
      dcu(2,nn,mm) = dcu(2,nn,mm) + vy*dy
      dcu(3,nn,mm) = dcu(3,nn,mm) + vz*dy
      cu(1,nn,mm) = cu(1,nn,mm) + ox*dy
      cu(2,nn,mm) = cu(2,nn,mm) + oy*dy
      cu(3,nn,mm) = cu(3,nn,mm) + oz*dy
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DSORTP2YL(parta,partb,npic,idimp,nop,ny1)
c this subroutine sorts particles by y grid
c linear interpolation
c parta/partb = input/output particle arrays
c parta(2,n) = position y of particle n
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c ny1 = system length in y direction + 1
      implicit none
      integer npic, idimp, nop, ny1
      real parta, partb
      dimension parta(idimp,nop), partb(idimp,nop), npic(ny1)
c local data
      integer i, j, k, m, isum, ist, ip
c clear counter array
      do 10 k = 1, ny1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      m = parta(2,j)
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
      m = parta(2,j)
      m = m + 1
      ip = npic(m) + 1
      do 40 i = 1, idimp
      partb(i,ip) = parta(i,j)
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
c nye = second dimension of field arrays, must be >= ny+1
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
      do 30 k = 1, ny
      do 20 j = 1, nx
      do 10 i = 1, 3
      dcu(i,j,k) = dcu(i,j,k) - q2m0*cus(i,j,k)
   10 continue
   20 continue
   30 continue
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
      do 20 k = 1, ny
      do 10 j = 1, nx
      at1 = qbme*qe(j,k)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine POIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,   
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
   50 continue
c mode numbers kx = 0, nx/2
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
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
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
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
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
   20 continue
c mode numbers kx = 0, nx/2
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
      subroutine BBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c with periodic boundary conditions.
c input: cu,ffc,ci,nx,ny,nxvh,nxhd,nyhd, output: bxy,wm
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
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
c calculate smoothed magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
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
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
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
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
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
      do 20 k = 1, ny
      do 10 j = 1, nx
      bxy(1,j,k) = bxy(1,j,k) + omx
      bxy(2,j,k) = bxy(2,j,k) + omy
      bxy(3,j,k) = bxy(3,j,k) + omz
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DCUPERP23(dcu,amu,nx,ny,nxvh,nyv)
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
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      zt2 = cmplx(aimag(amu(2,1,k)),-real(amu(2,1,k)))
      dcu(1,1,k) = dky*zt2
      dcu(2,1,k) = zero
      zt2 = cmplx(aimag(amu(4,1,k)),-real(amu(4,1,k)))
      dcu(3,1,k) = dky*zt2
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1)),-real(amu(2,j,1)))
      dcu(1,j,1) = zero
      dcu(2,j,1) = dkx*zt2
      zt1 = cmplx(aimag(amu(3,j,1)),-real(amu(3,j,1)))
      dcu(3,j,1) = dkx*zt1
      dcu(1,j,k1) = zero
      dcu(2,j,k1) = zero
      dcu(3,j,k1) = zero
   40 continue
      dcu(1,1,1) = zero
      dcu(2,1,1) = zero
      dcu(3,1,1) = zero
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine ADCUPERP23(dcu,amu,nx,ny,nxvh,nyv)
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
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      zt2 = cmplx(aimag(amu(2,1,k)),-real(amu(2,1,k)))
      dcu(1,1,k) = dcu(1,1,k) + dky*zt2
      dcu(2,1,k) = zero
      zt2 = cmplx(aimag(amu(4,1,k)),-real(amu(4,1,k)))
      dcu(3,1,k) = dcu(3,1,k) + dky*zt2
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1)),-real(amu(2,j,1)))
      dcu(1,j,1) = zero
      dcu(2,j,1) = dcu(2,j,1) + dkx*zt2
      zt1 = cmplx(aimag(amu(3,j,1)),-real(amu(3,j,1)))
      dcu(3,j,1) = dcu(3,j,1) + dkx*zt1
      dcu(1,j,k1) = zero
      dcu(2,j,k1) = zero
      dcu(3,j,k1) = zero
   40 continue
      dcu(1,1,1) = zero
      dcu(2,1,1) = zero
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine EPOIS23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,  
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
      double precision wp
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
   30 if (isign.gt.0) go to 80
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
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
   50 continue
c mode numbers kx = 0, nx/2
      do 60 k = 2, nyh
      k1 = ny2 - k
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
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
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
   70 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = real(nx*ny)*wp/real(ffe(1,1))
      return
c calculate unsmoothed transverse electric field and sum field energy
   80 wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 2, nxh
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
   90 continue
  100 continue
c mode numbers kx = 0, nx/2
      do 110 k = 2, nyh
      k1 = ny2 - k
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
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nxh
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
  120 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = real(nx*ny)*wp/real(ffe(1,1))
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
      do 30 k = 1, nye
      do 20 j = 1, nxe
      do 10 i = 1, ndim
      a(i,j,k) = b(i,j,k) + c(i,j,k)
   10 continue
   20 continue
   30 continue
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
      subroutine WFFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,    
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
         call FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
c perform y fft
         call FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c perform x fft
         call FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT2R3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,    
     1nxyhd)
c wrapper function for 3 2d real to complex ffts
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
         call FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
c perform y fft
         call FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c perform x fft
         call FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT2RN(f,ss,isign,mixup,sct,indx,indy,nxhd,nyd,ndim,  
     1nxhyd,nxyhd)
c wrapper function for multiple 2d real to complex ffts
      implicit none
      complex f, ss, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, ndim, nxhyd, nxyhd
      dimension f(ndim,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      dimension ss(ndim,nxhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2RNX(f,ss,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,   
     1ndim,nxhyd,nxyhd)
c perform y fft
         call FFT2RNY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,ndim,
     1nxhyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RNY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,ndim,
     1nxhyd,nxyhd)
c perform x fft
         call FFT2RNX(f,ss,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,   
     1ndim,nxhyd,nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd,  
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
c f(j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
c except for f(1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c aimag(f(1,1)) = real part of mode nx/2,0 and
c aimag(f(1,ny/2+1)) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr
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
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = nyi, nyt
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1.0/real(2*nx*ny)
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = nyi, nyt
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
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 k = nyi, nyt
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
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 160 i = nyi, nyt
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
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
c f(j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
c except for f(1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c aimag(f(1,1)) = real part of mode nx/2,0 and
c aimag(f(1,ny/2+1)) = real part of mode nx/2,ny/2
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
      subroutine FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd,  
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
c f(1:3,j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
c except for f(1:3,1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c aimag(f(1:3,1,1)) = real part of mode nx/2,0 and
c aimag(f(1:3,1,ny/2+1)) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
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
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(real(f(2,j,i)),aimag(f(3,j,i)))
      at2 = aimag(f(2,j,i))
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
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = nyi, nyt
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1.0/real(2*nx*ny)
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = nyi, nyt
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
      do 170 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 160 k = nyi, nyt
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
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 220 i = nyi, nyt
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
  220 continue
  230 continue
  240 continue
  250 continue
c swap complex components
      do 270 i = nyi, nyt
      do 260 j = 1, nxh
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(aimag(f(2,j,i)),aimag(f(3,j,i)))
      at2 = real(f(2,j,i))
      f(2,j,i) = cmplx(at1,aimag(f(1,j,i)))
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd,  
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
c f(1:3,j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
c except for f(1:3,1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c aimag(f(1:3,1,1)) = real part of mode nx/2,0 and
c aimag(f(1:3,1,ny/2+1)) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
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
      subroutine FFT2RNX(f,ss,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd
     1,ndim,nxhyd,nxyhd)
c this subroutine performs the x part of N two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c y, using complex arithmetic, where N = ndim
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
c f(1:N,j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
c except for f(1:N,1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c aimag(f(1:N,1,1)) = real part of mode nx/2,0 and
c aimag(f(1:N,1,ny/2+1)) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, ndim, nxhyd, nxyhd
      complex f, ss, sct
      integer mixup
      dimension f(ndim,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      dimension ss(ndim,nxhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, ns, ns2, km, kmr, i, j, k, l, k1, k2, j1, j2, jj
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
      if (isign.gt.0) go to 140
c inverse fourier transform
c swap complex components
      call SWAPC2N(f,ss,isign,nxh,nyi,nyt,nxhd,nyd,ndim)
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 30
      do 20 k = nyi, nyt
      do 10 jj = 1, ndim
      t1 = f(jj,j1,k)
      f(jj,j1,k) = f(jj,j,k)
      f(jj,j,k) = t1
   10 continue
   20 continue
   30 continue
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
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = nyi, nyt
      do 40 jj = 1, ndim
      t2 = t1*f(jj,j2,i)
      f(jj,j2,i) = f(jj,j1,i) - t2
      f(jj,j1,i) = f(jj,j1,i) + t2
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1.0/real(2*nx*ny)
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = nyi, nyt
      do 90 jj = 1, ndim
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,nxh2-j,k) = ani*conjg(t1 - t2)
   90 continue
  100 continue
  110 continue
      ani = 2.0*ani
      do 130 k = nyi, nyt
      do 120 jj = 1, ndim
      f(jj,nxhh+1,k) = ani*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = ani*cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),
     1                      real(f(jj,1,k)) - aimag(f(jj,1,k)))
  120 continue
  130 continue
      return
c forward fourier transform
c scramble coefficients
  140 kmr = nxy/nx
      do 170 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 160 k = nyi, nyt
      do 150 jj = 1, ndim
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,nxh2-j,k) = conjg(t1 - t2)
  150 continue
  160 continue
  170 continue
      do 190 k = nyi, nyt
      do 180 jj = 1, ndim
      f(jj,nxhh+1,k) = 2.0*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),
     1                  real(f(jj,1,k)) - aimag(f(jj,1,k)))
  180 continue
  190 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 220 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 220
      do 210 k = nyi, nyt
      do 200 jj = 1, ndim
      t1 = f(jj,j1,k)
      f(jj,j1,k) = f(jj,j,k)
      f(jj,j,k) = t1
  200 continue
  210 continue
  220 continue
c then transform in x
      nrx = nxy/nxh
      do 270 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 260 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 250 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 240 i = nyi, nyt
      do 230 jj = 1, ndim
      t2 = t1*f(jj,j2,i)
      f(jj,j2,i) = f(jj,j1,i) - t2
      f(jj,j1,i) = f(jj,j1,i) + t2
  230 continue
  240 continue
  250 continue
  260 continue
  270 continue
c swap complex components
      call SWAPC2N(f,ss,isign,nxh,nyi,nyt,nxhd,nyd,ndim)
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RNY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd,  
     1ndim,nxhyd,nxyhd)
c this subroutine performs the y part of N two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c x, using complex arithmetic, where N = ndim
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
c f(1:N,j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
c except for f(1:N,1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c aimag(f(1:N,1,1)) = real part of mode nx/2,0 and
c aimag(f(1:N,1,ny/2+1)) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, ndim, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(ndim,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, ns, ns2, km, kmr, i, j, k, l, k1, k2, j1, j2, jj
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
      if (isign.gt.0) go to 110
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 30 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 30
      do 20 j = nxi, nxt
      do 10 jj = 1, ndim
      t1 = f(jj,j,k1)
      f(jj,j,k1) = f(jj,j,k)
      f(jj,j,k) = t1
   10 continue
   20 continue
   30 continue
c then transform in y
      nry = nxy/ny
      do 80 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = nxi, nxt
      do 40 jj = 1, ndim
      t2 = t1*f(jj,i,j2)
      f(jj,i,j2) = f(jj,i,j1) - t2
      f(jj,i,j1) = f(jj,i,j1) + t2
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 100 k = 2, nyh
      if (nxi.eq.1) then
         do 90 jj = 1, ndim
         t1 = f(jj,1,ny2-k)
         f(jj,1,ny2-k) = 0.5*cmplx(aimag(f(jj,1,k) + t1),
     1                             real(f(jj,1,k) - t1))
         f(jj,1,k) = 0.5*cmplx(real(f(jj,1,k) + t1),
     1                         aimag(f(jj,1,k) - t1))
   90    continue
      endif
  100 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  110 do 130 k = 2, nyh
      if (nxi.eq.1) then
         do 120 jj = 1, ndim
         t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
         f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
         f(jj,1,k) = f(jj,1,k) + t1
  120    continue
      endif
  130 continue
c bit-reverse array elements in y
      nry = nxhy/ny
      do 160 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 160
      do 150 j = nxi, nxt
      do 140 jj = 1, ndim
      t1 = f(jj,j,k1)
      f(jj,j,k1) = f(jj,j,k)
      f(jj,j,k) = t1
  140 continue
  150 continue
  160 continue
c first transform in y
      nry = nxy/ny
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
      t1 = conjg(sct(1+kmr*(j-1)))
      do 180 i = nxi, nxt
      do 170 jj = 1, ndim
      t2 = t1*f(jj,i,j2)
      f(jj,i,j2) = f(jj,i,j1) - t2
      f(jj,i,j1) = f(jj,i,j1) + t2
  170 continue
  180 continue
  190 continue
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine SWAPC2N(f,s,isign,nxh,nyi,nyt,nxhd,nyd,ndim)
c this subroutine swaps components for multiple ffts
c f = input  array
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
      dimension f(ndim,2*nxhd,nyd), s(2*ndim*nxhd)
c local data
      integer i, j, k, ioff
c swap complex components
c real to complex
      if (isign.lt.0) then
         do 60 k = nyi, nyt
         do 20 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 10 i = 1, ndim
         s(2*i+ioff-1) = f(i,2*j-1,k)
         s(2*i+ioff) = f(i,2*j,k)
   10    continue
   20    continue
         do 50 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 30 i = 1, ndim
         f(i,2*j-1,k) = s(i+ioff)
   30    continue
         ioff = ioff + ndim
         do 40 i = 1, ndim
         f(i,2*j,k) = s(i+ioff)
   40    continue
   50    continue
   60    continue
c complex to real
      else if (isign.gt.0) then
         do 120 k = nyi, nyt
         do 90 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 70 i = 1, ndim
         s(i+ioff) = f(i,2*j-1,k)
   70    continue
         ioff = ioff + ndim
         do 80 i = 1, ndim
         s(i+ioff) = f(i,2*j,k)
   80    continue
   90    continue
         do 110 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 100 i = 1, ndim
         f(i,2*j-1,k) = s(2*i+ioff-1)
         f(i,2*j,k) = s(2*i+ioff)
  100    continue
  110    continue
  120    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine GSJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 37 flops/particle, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of current array, must be >= nx+1
c nxyv = dimension of current array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc
      real qm, dt
      real part, cu
      dimension part(idimp,nop), cu(3,nxyv)
c local data
      integer j, nnn, mmn, nn, mm, mp
      real dxn, dyn, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, vx, vy, vz, dx1, dy1
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
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
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j)
      mmn = part(2,j)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j) - real(nnn)
      dyn = part(2,j) - real(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyp
c deposit current
      dx = dxp*dyp
      dz = amx*dyp
      vx = part(3,j-1)
      vy = part(4,j-1)
      vz = part(5,j-1)
      dx1 = cu(1,mp+1) + vx*dx
      dy1 = cu(2,mp+1) + vy*dx
      dyp = cu(3,mp+1) + vz*dx
      dx = cu(1,mp) + vx*dz
      dy = cu(2,mp) + vy*dz
      dz = cu(3,mp) + vz*dz
      cu(1,mp+1) = dx1
      cu(2,mp+1) = dy1
      cu(3,mp+1) = dyp
      cu(1,mp) = dx
      cu(2,mp) = dy
      cu(3,mp) = dz
      dx = dxp*amy
      dz = amx*amy
      dxp = cu(1,mm+1) + vx*dx
      amx = cu(2,mm+1) + vy*dx
      dyp = cu(3,mm+1) + vz*dx
      dx = cu(1,mm) + vx*dz
      dy = cu(2,mm) + vy*dz
      dz = cu(3,mm) + vz*dz
      cu(1,mm+1) = dxp
      cu(2,mm+1) = amx
      cu(3,mm+1) = dyp
      cu(1,mm) = dx
      cu(2,mm) = dy
      cu(3,mm) = dz
c advance position half a time-step
      dx = part(1,j-1) + vx*dt
      dy = part(2,j-1) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1)
            part(4,j-1) = -part(4,j-1)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j-1) = dx
      part(2,j-1) = dy
   10 continue
c deposit current for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyn
c deposit current
      dx = dxp*dyn
      dy = amx*dyn
      vx = part(3,nop)
      vy = part(4,nop)
      vz = part(5,nop)
      cu(1,mp+1) = cu(1,mp+1) + vx*dx
      cu(2,mp+1) = cu(2,mp+1) + vy*dx
      cu(3,mp+1) = cu(3,mp+1) + vz*dx
      cu(1,mp) = cu(1,mp) + vx*dy
      cu(2,mp) = cu(2,mp) + vy*dy
      cu(3,mp) = cu(3,mp) + vz*dy
      dx = dxp*amy
      dy = amx*amy
      cu(1,mm+1) = cu(1,mm+1) + vx*dx
      cu(2,mm+1) = cu(2,mm+1) + vy*dx
      cu(3,mm+1) = cu(3,mm+1) + vz*dx
      cu(1,mm) = cu(1,mm) + vx*dy
      cu(2,mm) = cu(2,mm) + vy*dy
      cu(3,mm) = cu(3,mm) + vz*dy
c advance position half a time-step
      dx = part(1,nop) + vx*dt
      dy = part(2,nop) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
      return
      end
c-----------------------------------------------------------------------
      subroutine GSMJPOST2L(part,amu,qm,nop,idimp,nxv,nxyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using first-order linear interpolation
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
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
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x velocity of particle n at t - dt/2
c part(4,n) = y velocity of particle n at t - dt/2
c part(5,n) = z velocity of particle n at t - dt/2
c amu(i,n) = ith component of momentum flux at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 5
c nxv = first virtual dimension of current array, must be >= nx+1
c nxyv = dimension of current array, must be >= nxv*(ny+1)
      implicit none
      integer nop, idimp, nxv, nxyv
      real part, amu, qm
      dimension part(idimp,nop), amu(4,nxyv)
      integer nnn, mmn, nn, mm, mp, j
      real dxn, dyn, dxp, dyp, amx, amy, dx, dy, dz, vx, vy, vz
      real v1, v2, v3, v4, dx1, dy1
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j)
      mmn = part(2,j)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j) - real(nnn)
      dyn = part(2,j) - real(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyp
c deposit momentum flux
      dx = dxp*dyp
      dz = amx*dyp
      vx = part(3,j-1)
      vy = part(4,j-1)
      vz = part(5,j-1)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      dx1 = amu(1,mp+1) + v1*dx
      dy1 = amu(2,mp+1) + v2*dx
      dyp = amu(3,mp+1) + v3*dx
      vx = amu(4,mp+1) + v4*dx
      dx = amu(1,mp) + v1*dz
      dy = amu(2,mp) + v2*dz
      vy = amu(3,mp) + v3*dz
      dz = amu(4,mp) + v4*dz
      amu(1,mp+1) = dx1
      amu(2,mp+1) = dy1
      amu(3,mp+1) = dyp
      amu(4,mp+1) = vx
      amu(1,mp) = dx
      amu(2,mp) = dy
      amu(3,mp) = vy
      amu(4,mp) = dz
      dx = dxp*amy
      dz = amx*amy
      dxp = amu(1,mm+1) + v1*dx
      amx = amu(2,mm+1) + v2*dx
      dyp = amu(3,mm+1) + v3*dx
      vx = amu(4,mm+1) + v4*dx
      dx = amu(1,mm) + v1*dz
      dy = amu(2,mm) + v2*dz
      vy = amu(3,mm) + v3*dz
      dz = amu(4,mm) + v4*dz
      amu(1,mm+1) = dxp
      amu(2,mm+1) = amx
      amu(3,mm+1) = dyp
      amu(4,mm+1) = vx
      amu(1,mm) = dx
      amu(2,mm) = dy
      amu(3,mm) = vy
      amu(4,mm) = dz
   10 continue
c deposit momentum flux for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyn
c deposit momentum flux
      dx = dxp*dyn
      dy = amx*dyn
      vx = part(3,nop)
      vy = part(4,nop)
      vz = part(5,nop)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      amu(1,mp+1) = amu(1,mp+1) + v1*dx
      amu(2,mp+1) = amu(2,mp+1) + v2*dx
      amu(3,mp+1) = amu(3,mp+1) + v3*dx
      amu(4,mp+1) = amu(4,mp+1) + v4*dx
      amu(1,mp) = amu(1,mp) + v1*dy
      amu(2,mp) = amu(2,mp) + v2*dy
      amu(3,mp) = amu(3,mp) + v3*dy
      amu(4,mp) = amu(4,mp) + v4*dy
      dx = dxp*amy
      dy = amx*amy
      amu(1,mm+1) = amu(1,mm+1) + v1*dx
      amu(2,mm+1) = amu(2,mm+1) + v2*dx
      amu(3,mm+1) = amu(3,mm+1) + v3*dx
      amu(4,mm+1) = amu(4,mm+1) + v4*dx
      amu(1,mm) = amu(1,mm) + v1*dy
      amu(2,mm) = amu(2,mm) + v2*dy
      amu(3,mm) = amu(3,mm) + v3*dy
      amu(4,mm) = amu(4,mm) + v4*dy
      return
      end
c-----------------------------------------------------------------------
      subroutine GSDJPOST2L(part,fxy,bxy,dcu,amu,qm,qbm,dt,idimp,nop,nxv
     1,nxyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c and acceleration density using first-order spline interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 194 flops/particle, 1 divide, 57 loads, 28 stores
c input: all, output: dcu, amu
c acceleration density is approximated by values at the grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c dcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c where n,m = leftmost grid points and dx = x-n, dy = y-m
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
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c dcu(i,n) = ith component of acceleration density at grid point j,k
c where n = j + nxv*(k-1)
c amu(i,n) = ith component of momentum at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 5
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
      implicit none
      integer idimp, nop, nxv, nxyv
      real part, fxy, bxy, dcu, amu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), bxy(3,nxyv)
      dimension dcu(3,nxyv), amu(4,nxyv)
      integer nnn, mmn, nop1, j, nn, mm, mp
      real qtmh, dti, dxn, dyn, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, v3, v4, dx1, dy1, dx2, dy2, dx3, dy3
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1)
      mmn = part(2,j+1)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1) - real(nnn)
      dyn = part(2,j+1) - real(mmn)
      mm = mm + nn
      amx = 1. - dxp
      mp = mm + nxv
      amy = 1. - dyp
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyp*(dxp*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxp*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxp*bxy(1,mm+1)
     1+ amx*bxy(1,mm))
      oy = dyp*(dxp*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxp*bxy(2,mm+1)
     1+ amx*bxy(2,mm))
      oz = dyp*(dxp*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxp*bxy(3,mm+1)
     1+ amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
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
      dx = dxp*dyp
      dz = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      dx1 = amu(1,mp+1) + v1*dx
      dy1 = amu(2,mp+1) + v2*dx
      dyp = amu(3,mp+1) + v3*dx
      dx2 = amu(4,mp+1) + v4*dx
      dy2 = amu(1,mp) + v1*dz
      dx3 = amu(2,mp) + v2*dz
      dy3 = amu(3,mp) + v3*dz
      dy = amu(4,mp) + v4*dz
      amu(1,mp+1) = dx1
      amu(2,mp+1) = dy1
      amu(3,mp+1) = dyp
      amu(4,mp+1) = dx2
      amu(1,mp) = dy2
      amu(2,mp) = dx3
      amu(3,mp) = dy3
      amu(4,mp) = dy
      dx1 = dcu(1,mp+1) + vx*dx
      dy1 = dcu(2,mp+1) + vy*dx
      dyp = dcu(3,mp+1) + vz*dx
      dx2 = dcu(1,mp) + vx*dz
      dy2 = dcu(2,mp) + vy*dz
      dy = dcu(3,mp) + vz*dz
      dcu(1,mp+1) = dx1
      dcu(2,mp+1) = dy1
      dcu(3,mp+1) = dyp
      dcu(1,mp) = dx2
      dcu(2,mp) = dy2
      dcu(3,mp) = dy
      dx = dxp*amy
      dz = amx*amy
      dxp = amu(1,mm+1) + v1*dx
      amx = amu(2,mm+1) + v2*dx
      dyp = amu(3,mm+1) + v3*dx
      dx1 = amu(4,mm+1) + v4*dx
      dy1 = amu(1,mm) + v1*dz
      dx2 = amu(2,mm) + v2*dz
      dy2 = amu(3,mm) + v3*dz
      dy = amu(4,mm) + v4*dz
      amu(1,mm+1) = dxp
      amu(2,mm+1) = amx
      amu(3,mm+1) = dyp
      amu(4,mm+1) = dx1
      amu(1,mm) = dy1
      amu(2,mm) = dx2
      amu(3,mm) = dy2
      amu(4,mm) = dy
      dxp = dcu(1,mm+1) + vx*dx
      amx = dcu(2,mm+1) + vy*dx
      dyp = dcu(3,mm+1) + vz*dx
      dx1 = dcu(1,mm) + vx*dz
      dy1 = dcu(2,mm) + vy*dz
      dy = dcu(3,mm) + vz*dz
      dcu(1,mm+1) = dxp
      dcu(2,mm+1) = amx
      dcu(3,mm+1) = dyp
      dcu(1,mm) = dx1
      dcu(2,mm) = dy1
      dcu(3,mm) = dy
   10 continue
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1. - dxn
      mp = mm + nxv
      amy = 1. - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyn*(dxn*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxn*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxn*bxy(1,mm+1)
     1+ amx*bxy(1,mm))
      oy = dyn*(dxn*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxn*bxy(2,mm+1)
     1+ amx*bxy(2,mm))
      oz = dyn*(dxn*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxn*bxy(3,mm+1)
     1+ amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,nop)
      vy = part(4,nop)
      vz = part(5,nop)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxn
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = dxp*dyn
      dy = amx*dyn
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      amu(1,mp+1) = amu(1,mp+1) + v1*dx
      amu(2,mp+1) = amu(2,mp+1) + v2*dx
      amu(3,mp+1) = amu(3,mp+1) + v3*dx
      amu(4,mp+1) = amu(4,mp+1) + v4*dx
      amu(1,mp) = amu(1,mp) + v1*dy
      amu(2,mp) = amu(2,mp) + v2*dy
      amu(3,mp) = amu(3,mp) + v3*dy
      amu(4,mp) = amu(4,mp) + v4*dy
      dcu(1,mp+1) = dcu(1,mp+1) + vx*dx
      dcu(2,mp+1) = dcu(2,mp+1) + vy*dx
      dcu(3,mp+1) = dcu(3,mp+1) + vz*dx
      dcu(1,mp) = dcu(1,mp) + vx*dy
      dcu(2,mp) = dcu(2,mp) + vy*dy
      dcu(3,mp) = dcu(3,mp) + vz*dy
      dx = dxp*amy
      dy = amx*amy
      amu(1,mm+1) = amu(1,mm+1) + v1*dx
      amu(2,mm+1) = amu(2,mm+1) + v2*dx
      amu(3,mm+1) = amu(3,mm+1) + v3*dx
      amu(4,mm+1) = amu(4,mm+1) + v4*dx
      amu(1,mm) = amu(1,mm) + v1*dy
      amu(2,mm) = amu(2,mm) + v2*dy
      amu(3,mm) = amu(3,mm) + v3*dy
      amu(4,mm) = amu(4,mm) + v4*dy
      dcu(1,mm+1) = dcu(1,mm+1) + vx*dx
      dcu(2,mm+1) = dcu(2,mm+1) + vy*dx
      dcu(3,mm+1) = dcu(3,mm+1) + vz*dx
      dcu(1,mm) = dcu(1,mm) + vx*dy
      dcu(2,mm) = dcu(2,mm) + vy*dy
      dcu(3,mm) = dcu(3,mm) + vz*dy
      return
      end
c-----------------------------------------------------------------------
      subroutine GSDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,nop
     1,nxv,nxyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using first-order spline
c interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 218 flops/particle, 1 divide, 69 loads, 40 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c dcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c where n,m = leftmost grid points and dx = x-n, dy = y-m
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
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c dcu(i,n) = ith component of acceleration density at grid point j,k
c where n = j + nxv*(k-1)
c amu(i,n) = ith component of momentum at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 5
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
      implicit none
      integer idimp, nop, nxv, nxyv
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), bxy(3,nxyv)
      dimension cu(3,nxyv), dcu(3,nxyv), amu(4,nxyv)
      integer nnn, mmn, nop1, j, nn, mm, mp
      real qtmh, dti, dxn, dyn, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, v3, v4, dx1, dy1, dx2, dy2, dx3, dy3
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1)
      mmn = part(2,j+1)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1) - real(nnn)
      dyn = part(2,j+1) - real(mmn)
      mm = mm + nn
      amx = 1. - dxp
      mp = mm + nxv
      amy = 1. - dyp
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyp*(dxp*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxp*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxp*bxy(1,mm+1)
     1+ amx*bxy(1,mm))
      oy = dyp*(dxp*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxp*bxy(2,mm+1)
     1+ amx*bxy(2,mm))
      oz = dyp*(dxp*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxp*bxy(3,mm+1)
     1+ amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
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
      dx = dxp*dyp
      dz = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      dx1 = amu(1,mp+1) + v1*dx
      dy1 = amu(2,mp+1) + v2*dx
      dyp = amu(3,mp+1) + v3*dx
      dx2 = amu(4,mp+1) + v4*dx
      dy2 = amu(1,mp) + v1*dz
      dx3 = amu(2,mp) + v2*dz
      dy3 = amu(3,mp) + v3*dz
      dy = amu(4,mp) + v4*dz
      amu(1,mp+1) = dx1
      amu(2,mp+1) = dy1
      amu(3,mp+1) = dyp
      amu(4,mp+1) = dx2
      amu(1,mp) = dy2
      amu(2,mp) = dx3
      amu(3,mp) = dy3
      amu(4,mp) = dy
      dx1 = dcu(1,mp+1) + vx*dx
      dy1 = dcu(2,mp+1) + vy*dx
      dyp = dcu(3,mp+1) + vz*dx
      dx2 = dcu(1,mp) + vx*dz
      dy2 = dcu(2,mp) + vy*dz
      dy = dcu(3,mp) + vz*dz
      dcu(1,mp+1) = dx1
      dcu(2,mp+1) = dy1
      dcu(3,mp+1) = dyp
      dcu(1,mp) = dx2
      dcu(2,mp) = dy2
      dcu(3,mp) = dy
      dx1 = cu(1,mp+1) + ox*dx
      dy1 = cu(2,mp+1) + oy*dx
      dyp = cu(3,mp+1) + oz*dx
      dx2 = cu(1,mp) + ox*dz
      dy2 = cu(2,mp) + oy*dz
      dy = cu(3,mp) + oz*dz
      cu(1,mp+1) = dx1
      cu(2,mp+1) = dy1
      cu(3,mp+1) = dyp
      cu(1,mp) = dx2
      cu(2,mp) = dy2
      cu(3,mp) = dy
      dx = dxp*amy
      dz = amx*amy
      dxp = amu(1,mm+1) + v1*dx
      amx = amu(2,mm+1) + v2*dx
      dyp = amu(3,mm+1) + v3*dx
      dx1 = amu(4,mm+1) + v4*dx
      dy1 = amu(1,mm) + v1*dz
      dx2 = amu(2,mm) + v2*dz
      dy2 = amu(3,mm) + v3*dz
      dy = amu(4,mm) + v4*dz
      amu(1,mm+1) = dxp
      amu(2,mm+1) = amx
      amu(3,mm+1) = dyp
      amu(4,mm+1) = dx1
      amu(1,mm) = dy1
      amu(2,mm) = dx2
      amu(3,mm) = dy2
      amu(4,mm) = dy
      dxp = dcu(1,mm+1) + vx*dx
      amx = dcu(2,mm+1) + vy*dx
      dyp = dcu(3,mm+1) + vz*dx
      dx1 = dcu(1,mm) + vx*dz
      dy1 = dcu(2,mm) + vy*dz
      dy = dcu(3,mm) + vz*dz
      dcu(1,mm+1) = dxp
      dcu(2,mm+1) = amx
      dcu(3,mm+1) = dyp
      dcu(1,mm) = dx1
      dcu(2,mm) = dy1
      dcu(3,mm) = dy
      dxp = cu(1,mm+1) + ox*dx
      amx = cu(2,mm+1) + oy*dx
      dyp = cu(3,mm+1) + oz*dx
      dx1 = cu(1,mm) + ox*dz
      dy1 = cu(2,mm) + oy*dz
      dy = cu(3,mm) + oz*dz
      cu(1,mm+1) = dxp
      cu(2,mm+1) = amx
      cu(3,mm+1) = dyp
      cu(1,mm) = dx1
      cu(2,mm) = dy1
      cu(3,mm) = dy
   10 continue
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1. - dxn
      mp = mm + nxv
      amy = 1. - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyn*(dxn*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxn*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxn*bxy(1,mm+1)
     1+ amx*bxy(1,mm))
      oy = dyn*(dxn*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxn*bxy(2,mm+1)
     1+ amx*bxy(2,mm))
      oz = dyn*(dxn*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxn*bxy(3,mm+1)
     1+ amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,nop)
      vy = part(4,nop)
      vz = part(5,nop)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxn
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = dxp*dyn
      dy = amx*dyn
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      amu(1,mp+1) = amu(1,mp+1) + v1*dx
      amu(2,mp+1) = amu(2,mp+1) + v2*dx
      amu(3,mp+1) = amu(3,mp+1) + v3*dx
      amu(4,mp+1) = amu(4,mp+1) + v4*dx
      amu(1,mp) = amu(1,mp) + v1*dy
      amu(2,mp) = amu(2,mp) + v2*dy
      amu(3,mp) = amu(3,mp) + v3*dy
      amu(4,mp) = amu(4,mp) + v4*dy
      dcu(1,mp+1) = dcu(1,mp+1) + vx*dx
      dcu(2,mp+1) = dcu(2,mp+1) + vy*dx
      dcu(3,mp+1) = dcu(3,mp+1) + vz*dx
      dcu(1,mp) = dcu(1,mp) + vx*dy
      dcu(2,mp) = dcu(2,mp) + vy*dy
      dcu(3,mp) = dcu(3,mp) + vz*dy
      cu(1,mp+1) = cu(1,mp+1) + ox*dx
      cu(2,mp+1) = cu(2,mp+1) + oy*dx
      cu(3,mp+1) = cu(3,mp+1) + oz*dx
      cu(1,mp) = cu(1,mp) + ox*dy
      cu(2,mp) = cu(2,mp) + oy*dy
      cu(3,mp) = cu(3,mp) + oz*dy
      dx = dxp*amy
      dy = amx*amy
      amu(1,mm+1) = amu(1,mm+1) + v1*dx
      amu(2,mm+1) = amu(2,mm+1) + v2*dx
      amu(3,mm+1) = amu(3,mm+1) + v3*dx
      amu(4,mm+1) = amu(4,mm+1) + v4*dx
      amu(1,mm) = amu(1,mm) + v1*dy
      amu(2,mm) = amu(2,mm) + v2*dy
      amu(3,mm) = amu(3,mm) + v3*dy
      amu(4,mm) = amu(4,mm) + v4*dy
      dcu(1,mm+1) = dcu(1,mm+1) + vx*dx
      dcu(2,mm+1) = dcu(2,mm+1) + vy*dx
      dcu(3,mm+1) = dcu(3,mm+1) + vz*dx
      dcu(1,mm) = dcu(1,mm) + vx*dy
      dcu(2,mm) = dcu(2,mm) + vy*dy
      dcu(3,mm) = dcu(3,mm) + vz*dy
      cu(1,mm+1) = cu(1,mm+1) + ox*dx
      cu(2,mm+1) = cu(2,mm+1) + oy*dx
      cu(3,mm+1) = cu(3,mm+1) + oz*dx
      cu(1,mm) = cu(1,mm) + ox*dy
      cu(2,mm) = cu(2,mm) + oy*dy
      cu(3,mm) = cu(3,mm) + oz*dy
      return
      end
c-----------------------------------------------------------------------
      subroutine GSPOST2L(part,q,qm,nop,idimp,nxv,nxyv)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(j,k) = charge density at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first virtual dimension of charge array, must be >= nx+1
c nxyv = dimension of charge array, must be >= nxv*(ny+1)
      implicit none
      integer nop, idimp, nxv, nxyv
      real qm
      real part, q
      dimension part(idimp,nop), q(nxyv)
c local data
      integer j, nnn, mmn, nn, mm, mp
      real dxn, dyn, dxp, dyp, amx, amy, dx1
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j)
      mmn = part(2,j)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j) - real(nnn)
      dyn = part(2,j) - real(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyp
c deposit charge
      dx1 = q(mp+1) + dxp*dyp
      dyp = q(mp) + amx*dyp
      dxp = q(mm+1) + dxp*amy
      amy = q(mm) + amx*amy
      q(mp+1) = dx1
      q(mp) = dyp
      q(mm+1) = dxp
      q(mm) = amy
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyn
c deposit charge
      q(mp+1) = q(mp+1) + dxp*dyn
      q(mp) = q(mp) + amx*dyn
      q(mm+1) = q(mm+1) + dxp*amy
      q(mm) = q(mm) + amx*amy
      return
      end
c-----------------------------------------------------------------------
      subroutine GSBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny, 
     1nxv,nxyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
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
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real qbm, dt, dtc, ek
      real part, fxy, bxy
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), bxy(3,nxyv)
c local data
      integer j, nnn, mmn, nop1, nn, mm, mp
      real dxn, dyn, qtmh, edgelx, edgely, edgerx, edgery
      real dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz, acx, acy, acz
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
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
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1)
      mmn = part(2,j+1)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1) - real(nnn)
      dyn = part(2,j+1) - real(mmn)
      mm = mm + nn
      amx = 1.0 - dxp
      mp = mm + nxv
      amy = 1.0 - dyp
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp))                        
     1   + amy*(dxp*fxy(1,mm+1) + amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp))                        
     1   + amy*(dxp*fxy(2,mm+1) + amx*fxy(2,mm))
      dz = dyp*(dxp*fxy(3,mp+1) + amx*fxy(3,mp))                        
     1   + amy*(dxp*fxy(3,mm+1) + amx*fxy(3,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1) + amx*bxy(1,mp))                        
     1   + amy*(dxp*bxy(1,mm+1) + amx*bxy(1,mm))
      oy = dyp*(dxp*bxy(2,mp+1) + amx*bxy(2,mp))                        
     1   + amy*(dxp*bxy(2,mm+1) + amx*bxy(2,mm))
      oz = dyp*(dxp*bxy(3,mp+1) + amx*bxy(3,mp))                        
     1   + amy*(dxp*bxy(3,mm+1) + amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
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
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1.0 - dxn
      mp = mm + nxv
      amy = 1.0 - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp))                        
     1   + amy*(dxn*fxy(1,mm+1) + amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp))                        
     1   + amy*(dxn*fxy(2,mm+1) + amx*fxy(2,mm))
      dz = dyn*(dxn*fxy(3,mp+1) + amx*fxy(3,mp))                        
     1   + amy*(dxn*fxy(3,mm+1) + amx*fxy(3,mm))
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1) + amx*bxy(1,mp))                        
     1   + amy*(dxn*bxy(1,mm+1) + amx*bxy(1,mm))
      oy = dyn*(dxn*bxy(2,mp+1) + amx*bxy(2,mp))                        
     1   + amy*(dxn*bxy(2,mm+1) + amx*bxy(2,mm))
      oz = dyn*(dxn*bxy(3,mp+1) + amx*bxy(3,mp))                        
     1   + amy*(dxn*bxy(3,mm+1) + amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
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
      part(3,nop) = dx
      part(4,nop) = dy
      part(5,nop) = dz
c new position
      dx = part(1,nop) + dx*dtc
      dy = part(2,nop) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c normalize kinetic energy
   20 ek = ek + 0.5*sum1
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
