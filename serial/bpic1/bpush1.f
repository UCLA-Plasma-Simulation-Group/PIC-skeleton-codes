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
      subroutine GRBPUSH13L(part,fxyz,byz,omx,qbm,dt,dtc,ci,ek,idimp,nop
     1,nx,nxv,ipbc)
c for 1-2/2d code, this subroutine updates particle co-ordinate and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells
c 90 flops/particle, 4 divides, 2 sqrts, 14 loads, 4 stores
c input: all, output: part, ek
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
c part(1,n) = position x of particle n
c part(2,n) = momentum px of particle n
c part(3,n) = momentum py of particle n
c part(4,n) = momentum pz of particle n
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxv = second dimension of field arrays, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fxyz, byz, omx, qbm, dt, dtc, ci, ek
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
c local data
      integer j, nn
      real qtmh, ci2, edgelx, edgerx, dxp, amx, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real p2, gami, qtmg, dtg
      double precision sum1
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(2,j) + dx
      acy = part(3,j) + dy
      acz = part(4,j) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = omx
      oy = amx*byz(1,nn) + dxp*byz(1,nn+1)
      oz = amx*byz(2,nn) + dxp*byz(2,nn+1)
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
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
      ek = ek + sum1
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
      subroutine GRJPOST1L(part,cu,qm,dt,ci,nop,idimp,nx,nxv,ipbc)
c for 1-2/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 24 flops/particle, 1 divide, 1 sqrt, 9 loads, 5 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*pi*gami, where i = y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = x momentum of particle n
c part(3,n) = y momentum of particle n
c part(4,n) = z momentum of particle n
c cu(i,j) = ith component of current density at grid point j
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of current array, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer nop, idimp, nx, nxv, ipbc
      real part, cu, qm, dt, ci
      dimension part(idimp,nop), cu(2,nxv)
c local data
      integer j, nn
      real edgelx, ci2, edgerx, dxp, amx, vx, vy, vz, dx, p2, gami
      ci2 = ci*ci
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
c find inverse gamma
      vx = part(2,j)
      vy = part(3,j)
      vz = part(4,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c deposit current
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      cu(1,nn) = cu(1,nn) + vy*amx
      cu(2,nn) = cu(2,nn) + vz*amx
      cu(1,nn+1) = cu(1,nn+1) + vy*dxp
      cu(2,nn+1) = cu(2,nn+1) + vz*dxp
c advance position half a time-step
      dx = part(1,j) + vx*dt
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
