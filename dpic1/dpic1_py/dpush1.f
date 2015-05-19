c Fortran Library for Skeleton 1-2/2D Darwin PIC Code
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
      subroutine GBPUSH13L(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,nop,nx,
     1nxv,ipbc)
c for 1-2/2d code, this subroutine updates particle co-ordinate and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells
c 78 flops/particle, 1 divide, 14 loads, 4 stores
c input: all, output: part, ek
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
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c fxyz(1,j) = x component of force/charge at grid (j)
c fxyz(2,j) = y component of force/charge at grid (j)
c fxyz(3,j) = z component of force/charge at grid (j)
c that is, convolution of electric field over particle shape
c byz(1,j) = y component of magnetic field at grid (j)
c byz(2,j) = z component of magnetic field at grid (j)
c that is, the convolution of magnetic field over particle shape
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxv = second dimension of field arrays, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fxyz, byz, omx, qbm, dt, dtc, ek
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
c local data
      integer j, nn
      real qtmh, edgelx, edgerx, dxp, amx, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      dxp = part(1,j) - real(nn)
      nn = nn + 1
      amx = 1.0 - dxp
c find electric field
      dx = amx*fxyz(1,nn) + dxp*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + dxp*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + dxp*fxyz(3,nn+1)
c find magnetic field
      ox = omx
      oy = amx*byz(1,nn) + dxp*byz(1,nn+1)
      oz = amx*byz(2,nn) + dxp*byz(2,nn+1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(2,j) + dx
      acy = part(3,j) + dy
      acz = part(4,j) + dz
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
      part(2,j) = dx
      part(3,j) = dy
      part(4,j) = dz
c new position
      dx = part(1,j) + dx*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(2,j) = -part(2,j)
         endif
      endif
c set new position
      part(1,j) = dx
   10 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine GJPOST1L(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
c for 1-2/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 14 flops/particle, 8 loads, 5 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vi, where i = y,z
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c cu(i,j) = ith component of current density at grid point j
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of current array, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer nop, idimp, nx, nxv, ipbc
      real part, cu, qm, dt
      dimension part(idimp,nop), cu(2,nxv)
c local data
      integer j, nn
      real edgelx, edgerx, dxp, amx, vy, vz, dx
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) 
      dxp = qm*(part(1,j) - real(nn))
      nn = nn + 1
      amx = qm - dxp
c deposit current
      vy = part(3,j)
      vz = part(4,j)
      cu(1,nn) = cu(1,nn) + vy*amx
      cu(2,nn) = cu(2,nn) + vz*amx
      cu(1,nn+1) = cu(1,nn+1) + vy*dxp
      cu(2,nn+1) = cu(2,nn+1) + vz*dxp
c advance position half a time-step
      dx = part(1,j) + part(2,j)*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(2,j) = -part(2,j)
         endif
      endif
c set new position
      part(1,j) = dx
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GPOST1L(part,q,qm,nop,idimp,nxv)
c for 1d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells
c 7 flops/particle, 3 loads, 3 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n)=qm*(1.-dx) and q(n+1)=qm*dx
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c q(j) = charge density at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 2
c nxv = first dimension of charge array, must be >= nx+1
      implicit none
      integer nop, idimp, nxv
      real part, q, qm
      dimension part(idimp,nop), q(nxv)
c local data
      integer j, nn
      real dx
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) 
      dx = qm*(part(1,j) - real(nn))
      nn = nn + 1
c deposit charge
      q(nn) = q(nn) + (qm - dx)
      q(nn+1) = q(nn+1) + dx
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GMJPOST1L(part,amu,qm,nop,idimp,nxv)
c for 1-2/2d code, this subroutine calculates particle momentum flux
c using first-order spline interpolation
c scalar version using guard cells
c 14 flops/particle, 8 loads, 4 stores
c input: all, output: part, amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c amu(i,j) = ith component of momentum flux at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of flux array, must be >= nx+1
      implicit none
      integer nop, idimp, nxv
      real part, amu, qm
      dimension part(idimp,nop), amu(2,nxv)
c local data
      integer j, nn
      real dxp, amx, vx, v1, v2
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) 
      dxp = qm*(part(1,j) - real(nn))
      nn = nn + 1
      amx = qm - dxp
c deposit momentum flux
      vx = part(2,j)
      v1 = vx*part(3,j)
      v2 = vx*part(4,j)
      amu(1,nn) = amu(1,nn) + v1*amx
      amu(2,nn) = amu(2,nn) + v2*amx
      amu(1,nn+1) = amu(1,nn+1) + v1*dxp
      amu(2,nn+1) = amu(2,nn+1) + v2*dxp
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GDJPOST1L(part,fxyz,byz,dcu,amu,omx,qm,qbm,dt,idimp,nop
     1,nxv)
c for 1-2/2d code, this subroutine calculates particle momentum flux
c and acceleration density using first-order spline interpolation.
c scalar version using guard cells
c 100 flops/particle, 1 divide, 12 loads, 8 stores
c input: all, output: dcu, amu
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n)=qci*(1.-dx) and dcu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*dvj/dt, where j = y,z, for i = 1, 2
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c velocity equations at t=t+dt/2 are calculated from:
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
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
c where n = nearest grid point and dx = x-n
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c fxyz(1,j) = x component of force/charge at grid (j)
c fxyz(2,j) = y component of force/charge at grid (j)
c fxyz(3,j) = z component of force/charge at grid (j)
c that is, convolution of electric field over particle shape
c byz(1,j) = y component of magnetic field at grid (j)
c byz(2,j) = z component of magnetic field at grid (j)
c that is, the convolution of magnetic field over particle shape
c dcu(i,j) = ith component of acceleration density at grid point j
c for i = 1, 2
c amu(i,j) = ith component of momentum flux at grid point j
c for i = 1, 2
c omx = magnetic field electron cyclotron frequency in x
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c nop = number of particles
c nxv = second dimension of field arrays, must be >= nx+1
      implicit none
      integer idimp, nop, nxv
      real part, fxyz, byz, dcu, amu, omx, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension dcu(2,nxv), amu(2,nxv)
c local data
      integer j, nn
      real qtmh, dti, dxp, amx, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      dxp = part(1,j) - real(nn)
      nn = nn + 1
      amx = 1.0 - dxp
c find electric field
      dx = amx*fxyz(1,nn) + dxp*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + dxp*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + dxp*fxyz(3,nn+1)
c find magnetic field
      ox = omx
      oy = amx*byz(1,nn) + dxp*byz(1,nn+1)
      oz = amx*byz(2,nn) + dxp*byz(2,nn+1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(2,j)
      vy = part(3,j)
      vz = part(4,j)
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
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      v1 = ox*oy
      v2 = ox*oz
      amu(1,nn) = amu(1,nn) + v1*amx
      amu(2,nn) = amu(2,nn) + v2*amx
      amu(1,nn+1) = amu(1,nn+1) + v1*dxp
      amu(2,nn+1) = amu(2,nn+1) + v2*dxp
      dcu(1,nn) = dcu(1,nn) + vy*amx
      dcu(2,nn) = dcu(2,nn) + vz*amx
      dcu(1,nn+1) = dcu(1,nn+1) + vy*dxp
      dcu(2,nn+1) = dcu(2,nn+1) + vz*dxp
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GDCJPOST1L(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt,idimp
     1,nop,nxv)
c for 1-2/2d code, this subroutine calculates particle momentum flux,
c acceleration density and current density using first-order spline
c interpolation.
c scalar version using guard cells
c 108 flops/particle, 1 divide, 16 loads, 12 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vj, where j = y,z, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n)=qci*(1.-dx) and dcu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*dvj/dt, where j = y,z, for i = 1, 2
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c velocity equations at t=t+dt/2 are calculated from:
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
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
c where n = nearest grid point and dx = x-n
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c fxyz(1,j) = x component of force/charge at grid (j)
c fxyz(2,j) = y component of force/charge at grid (j)
c fxyz(3,j) = z component of force/charge at grid (j)
c that is, convolution of electric field over particle shape
c byz(1,j) = y component of magnetic field at grid (j)
c byz(2,j) = z component of magnetic field at grid (j)
c that is, the convolution of magnetic field over particle shape
c cu(i,j) = ith component of current density at grid point j
c for i = 1, 2
c dcu(i,j) = ith component of acceleration density at grid point j
c for i = 1, 2
c amu(i,j) = ith component of momentum flux at grid point j
c for i = 1, 2
c omx = magnetic field electron cyclotron frequency in x
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c nop = number of particles
c nxv = second dimension of field arrays, must be >= nx+1
      implicit none
      integer idimp, nop, nxv
      real part, fxyz, byz, cu, dcu, amu, omx, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension cu(2,nxv), dcu(2,nxv), amu(2,nxv)
c local data
      integer j, nn
      real qtmh, dti, dxp, amx, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      dxp = part(1,j) - real(nn)
      nn = nn + 1
      amx = 1.0 - dxp
c find electric field
      dx = amx*fxyz(1,nn) + dxp*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + dxp*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + dxp*fxyz(3,nn+1)
c find magnetic field
      ox = omx
      oy = amx*byz(1,nn) + dxp*byz(1,nn+1)
      oz = amx*byz(2,nn) + dxp*byz(2,nn+1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(2,j)
      vy = part(3,j)
      vz = part(4,j)
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
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      v1 = ox*oy
      v2 = ox*oz
      amu(1,nn) = amu(1,nn) + v1*amx
      amu(2,nn) = amu(2,nn) + v2*amx
      amu(1,nn+1) = amu(1,nn+1) + v1*dxp
      amu(2,nn+1) = amu(2,nn+1) + v2*dxp
      dcu(1,nn) = dcu(1,nn) + vy*amx
      dcu(2,nn) = dcu(2,nn) + vz*amx
      dcu(1,nn+1) = dcu(1,nn+1) + vy*dxp
      dcu(2,nn+1) = dcu(2,nn+1) + vz*dxp
      cu(1,nn) = cu(1,nn) + oy*amx
      cu(2,nn) = cu(2,nn) + oz*amx
      cu(1,nn+1) = cu(1,nn+1) + oy*dxp
      cu(2,nn+1) = cu(2,nn+1) + oz*dxp
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DSORTP1XL(parta,partb,npic,idimp,nop,nx1)
c this subroutine sorts particles by grid
c linear interpolation
c parta/partb = input/output particle arrays
c parta(1,n) = position x of particle n
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c nx1 = system length in x direction + 1
      implicit none
      integer npic, idimp, nop, nx1
      real parta, partb
      dimension parta(idimp,nop), partb(idimp,nop), npic(nx1)
c local data
      integer i, j, k, n, isum, ist, ip
c clear counter array
      do 10 k = 1, nx1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      n = parta(1,j)
      n = n + 1
      npic(n) = npic(n) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nx1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, nop
      n = parta(1,j)
      n = n + 1
      ip = npic(n) + 1
      do 40 i = 1, idimp
      partb(i,ip) = parta(i,j)
   40 continue
      npic(n) = ip
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DGUARD1L(fx,nx,nxe)
c replicate extended periodic scalar field fx
c linear interpolation
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      real fx
      integer nx, nxe
      dimension fx(nxe)
      fx(nx+1) = fx(1)
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
      subroutine ASCFGUARD1L(dcu,cus,q2m0,nx,nxe)
c add scaled field to extended periodic field
c linear interpolation
      implicit none
      real dcu, cus, q2m0
      integer nx, nxe
      dimension dcu(2,nxe), cus(2,nxe)
c local data
      integer i, j
      do 20 j = 1, nx
      do 10 i = 1, 2
      dcu(i,j) = dcu(i,j) - q2m0*cus(i,j)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FWPMINMX1(qe,qbme,wpmax,wpmin,nx,nxe)
c calculates maximum and minimum plasma frequency.  assumes guard cells
c have already been added
c qe = charge density for electrons
c qbme = charge/mass ratio for electrons
c wpmax/wpmin = maximum/minimum plasma frequency
c nx = system length in x direction
c nxe = first dimension of charge arrays, nxe must be >= nx
      implicit none
      real qe, qbme, wpmax, wpmin
      integer nx, nxe
      dimension qe(nxe)
c local data
      integer j
      real at1
      wpmax = qbme*qe(1)
      wpmin = wpmax
      do 10 j = 1, nx
      at1 = qbme*qe(j)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
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
      subroutine BBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
c this subroutine solves 1-2/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c with periodic boundary conditions.
c input: cu,ffc,ci,nx,nxvh,nxhd output: byz,wm
c approximate flop count is: 30*nxc
c where nxc = nx/2 - 1
c the magnetic field is calculated using the equations:
c by(kx) = -ci*ci*sqrt(-1)*g(kx)*kx*cuz(kx)*s(kx),
c bz(kx) = ci*ci*sqrt(-1)*g(kx)*kx*cuy(kx)*s(kx),
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/kx**2)*s(kx),
c s(kx) = exp(-((kx*ax)**2+)/2), except for
c by(kx=pi) = bz(kx=pi) = 0, and by(kx=0) = bz(kx=0) = 0.
c cu(i,j) = complex current density for fourier mode (j-1)
c byz(1,j) = y component of complex magnetic field
c byz(2,j) = z component of complex magnetic field
c all for fourier mode (j-1)
c aimag(ffc(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffc(j)) = potential green's function g
c for fourier mode (j-1)
c ci = reciprocal of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2)
c affp = normalization constant = nx/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
c nxhd = dimension of form factor array, must be >= nxh
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
      at1 = ci2*real(ffc(j))*aimag(ffc(j))
      at2 = dnx*real(j - 1)*at1
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
      subroutine BADDEXT1(byz,omy,omz,nx,nxe)
c adds constant to magnetic field for 1-2/2d code
c byz = magnetic field
c omy/omz = magnetic field electron cyclotron frequency in y/z 
c nx = system length in x direction
c nxe = second dimension of magnetic field array, nxe must be >= nx
      implicit none
      real byz, omy, omz
      integer nx, nxe
      dimension byz(2,nxe)
c local data
      integer j
      do 10 j = 1, nx
      byz(1,j) = byz(1,j) + omy
      byz(2,j) = byz(2,j) + omz
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DCUPERP13(dcu,amu,nx,nxvh)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux
c in 1-2/2d with periodic boundary conditions.
c the transverse part of the derivative of the current is calculated
c using the equations:
c dcu(1,kx) = -sqrt(-1)*kx*vx*vy
c dcu(2,kx) = -sqrt(-1)*kx*vx*vz
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,kx=0) = 0.
c amu(1,j) = xy component of complex momentum flux
c amu(2,j) = xz component of complex momentum flux
c all for fourier mode (j-1)
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
      implicit none
      integer nx, nxvh
      complex dcu, amu
      dimension dcu(2,nxvh), amu(2,nxvh)
c local data
      integer nxh, j
      real dnx, dkx
      complex zero, zt1, zt2
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(1,j)),-real(amu(1,j)))
      dcu(1,j) = dkx*zt2
      zt1 = cmplx(aimag(amu(2,j)),-real(amu(2,j)))
      dcu(2,j) = dkx*zt1
   10 continue
      dcu(1,1) = zero
      dcu(2,1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine ADCUPERP13(dcu,amu,nx,nxvh)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux and acceleration density
c in 1-2/2d with periodic boundary conditions.
c the derivative of the current is calculated using the equations:
c dcu(1,kx) = dcu(1,kx)-sqrt(-1)*kx*vx*vy
c dcu(2,kx) = dcu(2,kx)-sqrt(-1)*kx*vx*vz
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,kx=0) = 0.
c on input:
c dcu(i,j) = complex acceleration density for fourier mode (j-1)
c on output:
c dcu(i,j) = transverse part of complex derivative of current for
c fourier mode (j-1)
c amu(1,j) = xy component of complex momentum flux
c amu(2,j) = xz component of complex momentum flux
c all for fourier mode (j-1)
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
      implicit none
      integer nx, nxvh
      complex dcu, amu
      dimension dcu(2,nxvh), amu(2,nxvh)
c local data
      integer nxh, j
      real dnx, dkx
      complex zero, zt1, zt2
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(1,j)),-real(amu(1,j)))
      dcu(1,j) = dcu(1,j) + dkx*zt2
      zt1 = cmplx(aimag(amu(2,j)),-real(amu(2,j)))
      dcu(2,j) = dcu(2,j) + dkx*zt1
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine EPOIS13(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,nxvh,   
     1nxhd)
c this subroutine solves 1-2/2d poisson's equation in fourier space for
c transverse electric field (or convolution of transverse electric field
c over particle shape), with periodic boundary conditions.
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c for isign = 0, input: isign,ax,affp,wp0,nx,nxvh, output:ffe
c for isign =/ 0, input: dcu,ffe,isign,ci,nx,nxvh,nxhd, output: eyz,wf
c approximate flop count is: 25*nxc
c where nxc = nx/2 - 1
c if isign = 0, form factor array is prepared
c if isign = -1, smoothed transverse electric field is calculated
c using the equation:
c ey(kx) = -ci*ci*g(kx)*dcuy(kx)*s(kx)
c ez(kx) = -ci*ci*g(kx)*dcuz(kx)*s(kx)
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/(kx**2+wp0*ci2*s(kx)**2))*s(kx),
c s(kx) = exp(-((kx*ax)**2+)/2), except for
c ey(kx=pi) = ez(kx=pi) = 0, and ey(kx=0) = ez(kx=0,) = 0.
c if isign = 1, unsmoothed transverse electric field is calculated
c using the equation:
c ey(kx) = -ci*ci*g(kx)*dcuy(kx)
c ez(kx) = -ci*ci*g(kx)*dcuz(kx)
c dcu(i,j) = transverse part of complex derivative of current for
c fourier mode (j-1)
c eyz(1,j) = y component of complex transverse electric field
c eyz(2,j) = z component of complex transverse electric field
c all for fourier mode (j-1)
c aimag(ffe(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffe(j)) = potential green's function g
c for fourier mode (j-1)
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np=number of particles
c wp0 = normalized total plasma frequency squared
c ci = reciprocal of velocity of light
c transverse electric field energy is also calculated, using
c wf = nx*sum((affp/(kx**2*ci*ci)**2)*|dcu(kx)*s(kx)|**2)
c this expression is valid only if the derivative of current is
c divergence-free
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
c nxhd = second dimension of form factor array, must be >= nxh
      implicit none
      integer isign, nx, nxvh, nxhd
      real ax, affp, wp0, ci, wf
      complex dcu, eyz, ffe
      dimension dcu(2,nxvh), eyz(2,nxvh)
      dimension ffe(nxhd)
c local data
      integer nxh, j
      real dnx, ci2, wpc, dkx, at1, at2
      complex zero
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
      if (isign.ne.0) go to 20
      wpc = wp0*ci2
c prepare form factor array
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at2 = exp(-.5*(dkx*ax)**2)
      ffe(j) = cmplx(affp*at2/(dkx*dkx+ wpc*at2*at2),at2)
   10 continue
      ffe(1) = cmplx(affp,1.0)
      return
c calculate smoothed transverse electric field and sum field energy
   20 if (isign.gt.0) go to 40
      wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 30 j = 2, nxh
      at2 = -ci2*real(ffe(j))
      at1 = at2*aimag(ffe(j))
      at2 = at2*at2
      eyz(1,j) = at1*dcu(1,j)
      eyz(2,j) = at1*dcu(2,j)
      wp = wp + at2*(dcu(1,j)*conjg(dcu(1,j)) + dcu(2,j)*conjg(dcu(2,j))
     1)
   30 continue
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*wp/real(ffe(1))
      return
c calculate unsmoothed transverse electric field and sum field energy
   40 wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 50 j = 2, nxh
      at2 = -ci2*real(ffe(j))
      at1 = at2*at2
      eyz(1,j) = at2*dcu(1,j)
      eyz(2,j) = at2*dcu(2,j)
      wp = wp + at1*(dcu(1,j)*conjg(dcu(1,j)) + dcu(2,j)*conjg(dcu(2,j))
     1)
   50 continue
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*wp/real(ffe(1))
      return
      end
c-----------------------------------------------------------------------
      subroutine ADDVRFIELD13(fxyze,eyze,fxe,nxe)
c this subroutine merges longitudinal and transverse electric fields
c fxyze(1,:) = fxe
c fxyze(2:3,:) = eyze
      implicit none
      integer nxe
      real fxyze, fxe, eyze
      dimension fxyze(3,nxe), eyze(2,nxe), fxe(nxe)
c local data
      integer j
      do 10 j = 1, nxe
      fxyze(1,j) = fxe(j)
      fxyze(2,j) = eyze(1,j)
      fxyze(3,j) = eyze(2,j)
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
