c Fortran Library for Skeleton 2-1/2D Electromagnetic MPI PIC Code
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
      subroutine PPGBPUSH23L(part,fxy,bxy,edges,npp,noff,ihole,qbm,dt,  
     1dtc,ek,nx,ny,idimp,npmax,nxv,nypmx,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space with magnetic field. Using the Boris Mover.
c scalar version using guard cells, for distributed data
c also determines list of particles which are leaving this processor
c 117 flops/particle, 1 divide, 25 loads, 5 stores
c input: all except ihole, output: part, ihole, ek
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
c fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = velocity vx of particle n in partition
c part(4,n) = velocity vy of particle n in partition
c part(5,n) = velocity vz of particle n in partition
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c fxy(3,j,k) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff - 1
c bxy(1,j,k) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c where kk = k + noff - 1
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = second dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c idps = number of partition boundaries
c ntmax = size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nypmx
      integer ipbc
      real qbm, dt, dtc, ek
      real part, fxy, bxy, edges
      integer ihole
      dimension part(idimp,npmax), fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, j, nn, mm, np, mp, ih, nh
      real qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = 0.5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      mnoff = noff - 1
      ih = 0
      nh = 0
      do 10 j = 1, npp
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm - mnoff
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
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c set end of file flag
      if (nh.gt.0) ih = -ih
      ihole(1) = ih
c normalize kinetic energy
      ek = ek + 0.5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGRBPUSH23L(part,fxy,bxy,edges,npp,noff,ihole,qbm,dt, 
     1dtc,ci,ek,nx,ny,idimp,npmax,nxv,nypmx,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells, for distributed data
c also determines list of particles which are leaving this processor
c 129 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
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
c fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = momentum px of particle n in partition
c part(4,n) = momentum py of particle n in partition
c part(5,n) = momentum pz of particle n in partition
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c fxy(3,j,k) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff - 1
c bxy(1,j,k) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c where kk = k + noff - 1
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = second dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c idps = number of partition boundaries
c ntmax = size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nypmx
      integer ipbc
      real qbm, dt, dtc, ci, ek
      real part, fxy, bxy, edges
      integer ihole
      dimension part(idimp,npmax), fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, j, nn, mm, np, mp, ih, nh
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      mnoff = noff - 1
      ih = 0
      nh = 0
      do 10 j = 1, npp
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm - mnoff
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp))                    
     1   + amy*(dxp*bxy(1,np,mm) + amx*bxy(1,nn,mm))
      oy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp))                    
     1   + amy*(dxp*bxy(2,np,mm) + amx*bxy(2,nn,mm))
      oz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp))                    
     1   + amy*(dxp*bxy(3,np,mm) + amx*bxy(3,nn,mm))
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
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c set end of file flag
      if (nh.gt.0) ih = -ih
      ihole(1) = ih
c normalize kinetic energy
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGPOST2L(part,q,npp,noff,qm,idimp,npmax,nxv,nypmx)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells, for distributed data
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c q(j,k) = charge density at grid point (j,kk),
c where kk = k + noff - 1
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nxv = first dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nxv, nypmx
      real qm
      real part, q
      dimension part(idimp,npmax), q(nxv,nypmx)
c local data
      integer mnoff, j, nn, np, mm, mp
      real dxp, dyp, amx, amy
      mnoff = noff - 1
      do 10 j = 1, npp
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = qm - dxp
      mp = mm + 1
      amy = 1.0 - dyp
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
      subroutine PPGJPOST2L(part,cu,edges,npp,noff,ihole,qm,dt,nx,ny,   
     1idimp,npmax,nxv,nypmx,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation, and distributed data.
c in addition, particle positions are advanced a half time-step
c also determines list of particles which are leaving this processor
c scalar version using guard cells, for distributed data
c 35 flops/particle, 17 loads, 14 stores
c input: all except ihole, output: part, ihole, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = x velocity of particle n in partition
c part(4,n) = y velocity of particle n in partition
c part(5,n) = z velocity of particle n in partition
c cu(i,j,k) = ith component of current density at grid point (j,kk),
c where kk = k + noff - 1
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = second dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nypmx
      integer ipbc
      real qm, dt
      real part, cu, edges
      integer ihole
      dimension part(idimp,npmax), cu(3,nxv,nypmx)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, j, nn, mm, np, mp, ih, nh
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, vx, vy, vz
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      mnoff = noff - 1
      ih = 0
      nh = 0
      do 10 j = 1, npp
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm - mnoff
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
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c set end of file flag
      if (nh.gt.0) ih = -ih
      ihole(1) = ih
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGRJPOST2L(part,cu,edges,npp,noff,ihole,qm,dt,ci,nx,ny
     1,idimp,npmax,nxv,nypmx,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles,
c and distributed data.
c in addition, particle positions are advanced a half time-step
c also determines list of particles which are leaving this processor
c scalar version using guard cells, for distributed data
c 45 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
c input: all except ihole, output: part, ihole, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = x momentum of particle n in partition
c part(4,n) = y momentum of particle n in partition
c part(5,n) = z momentum of particle n in partition
c cu(i,j,k) = ith component of current density at grid point (j,kk),
c where kk = k + noff - 1
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = second dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nypmx
      integer ipbc
      real qm, dt, ci
      real part, cu, edges
      integer ihole
      dimension part(idimp,npmax), cu(3,nxv,nypmx)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, j, nn, mm, np, mp, ih, nh
      real ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, vx, vy, vz, p2, gami
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
      mnoff = noff - 1
      ih = 0
      nh = 0
      do 10 j = 1, npp
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn + 1
      mm = mm - mnoff
      amx = qm - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
c deposit current
      dx = dxp*dyp
      dy = amx*dyp
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
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
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c set end of file flag
      if (nh.gt.0) ih = -ih
      ihole(1) = ih
      return
      end
c-----------------------------------------------------------------------
      subroutine PPDSORTP2YL(parta,partb,npic,npp,noff,nyp,idimp,npmax, 
     1nypm1)
c this subroutine sorts particles by y grid
c linear interpolation, spatial decomposition in y direction
c parta/partb = input/output particle array
c part(2,n) = position y of particle n in partition
c npic = address offset for reordering particles
c npp = number of particles in partition
c noff = backmost global gridpoint in particle partition
c nyp = number of primary gridpoints in particle partition
c idimp = size of phase space
c npmax = maximum number of particles in each partition
c nypm1 = maximum size of particle partition plus one
      implicit none
      integer idimp, npmax, nypm1
      integer npic, npp, noff, nyp
      real parta, partb
      dimension parta(idimp,npmax), partb(idimp,npmax)
      dimension npic(nypm1)
c local data
      integer i, j, k, m, mnoff, nyp1, isum, ist, ip
      mnoff = noff - 1
      nyp1 = nyp + 1
c clear counter array
      do 10 k = 1, nyp1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, npp
      m = parta(2,j)
      m = m - mnoff
      npic(m) = npic(m) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyp1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, npp
      m = parta(2,j)
      m = m - mnoff
      ip = npic(m) + 1
      do 40 i = 1, idimp
      partb(i,ip) = parta(i,j)
   40 continue
      npic(m) = ip
   50 continue
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
      subroutine PPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv, 
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
      double precision wp
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
   30 wp = 0.0d0
      if (kstrt.gt.nxh) go to 70
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 j = 1, kxps
      dkx = dnx*real(j + joff)
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
   50 continue
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
   70 continue
      we = real(nx)*real(ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPCUPERP2(cu,nx,ny,kstrt,nyv,kxp)
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
c nyv = seconsecond dimension of field arrays, must be >= ny
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
      subroutine IPPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field with periodic boundary conditions for distributed data.
c input: cu,ffc,ci,nx,ny,kstrt,nyv,kxp,nyhd, output: bxy,wm
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
c ci = reciprocal of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
c    |cu(kx,ky)*s(kx,ky)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = second dimension of field arrays, must be >= ny
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
      double precision wp
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
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 40
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
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
   20 continue
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
   40 continue
      wm = real(nx)*real(ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,kstrt, 
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
c ci = reciprocal of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
c nx/ny = system length in x/y direction
c kxp = number of data values per block
c kstrt = starting data block number
c nyv = second dimension of field arrays, must be >= ny
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
      double precision wp, ws
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
      ws = 0.0d0
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 40
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
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
   20 continue
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
   40 continue
      wf = real(nx)*real(ny)*ws
      wm = real(nx)*real(ny)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,nyhd)
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
c copy the fields
      else if (isign.lt.0) then
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
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
      subroutine WPPFFT2R(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy
     1,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
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
         call PPFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh
     1,kypd,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp,  
     1kypd)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PPFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,
     1kxp,nxhyd,nxyhd)
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
         call PPFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,
     1kxp,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,kypd, 
     1kxp)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PPFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh
     1,kypd,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPPFFT2R3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,   
     1indy,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
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
         call PPFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,   
     1nxvh,kypd,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,3,nxvh,nyv,kxp
     1,kypd)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PPFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv
     1,kxp,nxhyd,nxyhd)
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
         call PPFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv
     1,kxp,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,3,nyv,nxvh,   
     1kypd,kxp)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PPFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,   
     1nxvh,kypd,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp, 
     1nxvh,kypd,nxhyd,nxyhd)
c this subroutine performs the x part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks
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
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
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
      if (isign.gt.0) go to 100
c inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 k = kypi, kypt
      t = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t
   10 continue
   20 continue
c first transform in x
      nrx = nxy/nxh
      do 60 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = kypi, kypt
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 80 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = kypi, kypt
      t = conjg(f(nxh2-j,k))
      s = f(j,k) + t
      t = (f(j,k) - t)*t1
      f(j,k) = ani*(s + t)
      f(nxh2-j,k) = ani*conjg(s - t)
   70 continue
   80 continue
      do 90 k = kypi, kypt
      f(1,k) = 2.0*ani*cmplx(real(f(1,k)) + aimag(f(1,k)),              
     1                       real(f(1,k)) - aimag(f(1,k)))
      if (nxhh.gt.0) f(nxhh+1,k) = 2.0*ani*conjg(f(nxhh+1,k))
   90 continue
      return
c forward fourier transform
  100 kmr = nxy/nx
c scramble coefficients
      do 120 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 k = kypi, kypt
      t = conjg(f(nxh2-j,k))
      s = f(j,k) + t
      t = (f(j,k) - t)*t1
      f(j,k) = s + t
      f(nxh2-j,k) = conjg(s - t)
  110 continue
  120 continue
      do 130 k = kypi, kypt
      f(1,k) = cmplx(real(f(1,k)) + aimag(f(1,k)),                      
     1               real(f(1,k)) - aimag(f(1,k)))
      if (nxhh.gt.0) f(nxhh+1,k) = 2.0*conjg(f(nxhh+1,k))
  130 continue
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 150 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 150
      do 140 k = kypi, kypt
      t = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t
  140 continue
  150 continue
c then transform in x
      nrx = nxy/nxh
      do 190 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 160 i = kypi, kypt
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp, 
     1nyv,kxp,nxhyd,nxyhd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks
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
      integer isign, mixup, indx, indy, kstrt, kxpi, kxpp, nyv
      integer kxp, nxhyd, nxyhd
      complex g, sct
      dimension g(nyv,kxp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
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
      if (isign.gt.0) go to 80
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = kxpi, kxpt
      t = g(k1,j)
      g(k1,j) = g(k,j)
      g(k,j) = t
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
      do 60 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = kxpi, kxpt
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      if (ks.gt.0) return
      do 70 k = 2, nyh
      if (kxpi.eq.1) then
         s = g(ny2-k,1)
         g(ny2-k,1) = 0.5*cmplx(aimag(g(k,1) + s),real(g(k,1) - s))
         g(k,1) = 0.5*cmplx(real(g(k,1) + s),aimag(g(k,1) - s))
      endif
   70 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
   80 nry = nxhy/ny
      if (ks.gt.0) go to 100
      do 90 k = 2, nyh
      if (kxpi.eq.1) then
         s = cmplx(aimag(g(ny2-k,1)),real(g(ny2-k,1)))
         g(ny2-k,1) = conjg(g(k,1) - s)
         g(k,1) = g(k,1) + s
      endif
   90 continue
c bit-reverse array elements in y
  100 do 120 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 120
      do 110 j = kxpi, kxpt
      t = g(k1,j)
      g(k1,j) = g(k,j)
      g(k,j) = t
  110 continue
  120 continue
c first transform in y
      nry = nxy/ny
      do 160 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 150 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 140 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 130 i = kxpi, kxpt
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
  130 continue
  140 continue
  150 continue
  160 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,
     1nxvh,kypd,nxhyd,nxyhd)
c this subroutine performs the x part of 3 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks
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
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
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
      if (isign.gt.0) go to 140
c inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrx = nxhy/nxh
c swap complex components
      do 20 k = kypi, kypt
      do 10 j = 1, nxh
      at1 = real(f(3,j,k))
      f(3,j,k) = cmplx(real(f(2,j,k)),aimag(f(3,j,k)))
      at2 = aimag(f(2,j,k))
      f(2,j,k) = cmplx(aimag(f(1,j,k)),at1)
      f(1,j,k) = cmplx(real(f(1,j,k)),at2)
   10 continue
   20 continue
c bit-reverse array elements in x
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 k = kypi, kypt
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
      do 80 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 50 i = kypi, kypt
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      t3 = s*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(3,j2,i) = f(3,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
      f(3,j1,i) = f(3,j1,i) + t3
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 110 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = kypi, kypt
      do 90 i = 1, 3
      t = conjg(f(i,nxh2-j,k))
      s = f(i,j,k) + t
      t = (f(i,j,k) - t)*t1
      f(i,j,k) = ani*(s + t)
      f(i,nxh2-j,k) = ani*conjg(s - t)
   90 continue
  100 continue
  110 continue
      do 130 k = kypi, kypt
      do 120 i = 1, 3
      f(i,1,k) = 2.0*ani*cmplx(real(f(i,1,k)) + aimag(f(i,1,k)),        
     1                         real(f(i,1,k)) - aimag(f(i,1,k)))
      if (nxhh.gt.0) f(i,nxhh+1,k) = 2.0*ani*conjg(f(i,nxhh+1,k))
  120 continue
  130 continue
      return
c forward fourier transform
  140 kmr = nxy/nx
c scramble coefficients
      do 170 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 160 k = kypi, kypt
      do 150 i = 1, 3
      t = conjg(f(i,nxh2-j,k))
      s = f(i,j,k) + t
      t = (f(i,j,k) - t)*t1
      f(i,j,k) = s + t
      f(i,nxh2-j,k) = conjg(s - t)
  150 continue
  160 continue
  170 continue
      do 190 k = kypi, kypt
      do 180 i = 1, 3
      f(i,1,k) = cmplx(real(f(i,1,k)) + aimag(f(i,1,k)),                
     1                 real(f(i,1,k)) - aimag(f(i,1,k)))
      if (nxhh.gt.0) f(i,nxhh+1,k) = 2.0*conjg(f(i,nxhh+1,k))
  180 continue
  190 continue
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 210 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 210
      do 200 k = kypi, kypt
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
      do 250 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 220 i = kypi, kypt
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      t3 = s*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(3,j2,i) = f(3,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
      f(3,j1,i) = f(3,j1,i) + t3
  220 continue
  230 continue
  240 continue
  250 continue
c swap complex components
      do 270 k = kypi, kypt
      do 260 j = 1, nxh
      at1 = real(f(3,j,k))
      f(3,j,k) = cmplx(aimag(f(2,j,k)),aimag(f(3,j,k)))
      at2 = real(f(2,j,k))
      f(2,j,k) = cmplx(at1,aimag(f(1,j,k)))
      f(1,j,k) = cmplx(real(f(1,j,k)),at2)
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,
     1nyv,kxp,nxhyd,nxyhd)
c this subroutine performs the y part of 3 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks
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
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
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
      if (isign.gt.0) go to 90
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = kxpi, kxpt
      t1 = g(1,k1,j)
      t2 = g(2,k1,j)
      t3 = g(3,k1,j)
      g(1,k1,j) = g(1,k,j)
      g(2,k1,j) = g(2,k,j)
      g(3,k1,j) = g(3,k,j)
      g(1,k,j) = t1
      g(2,k,j) = t2
      g(3,k,j) = t3
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
      do 60 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = kxpi, kxpt
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      t3 = s*g(3,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(3,j2,i) = g(3,j1,i) - t3
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
      g(3,j1,i) = g(3,j1,i) + t3
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      if (ks.gt.0) return
      do 80 k = 2, nyh
      if (kxpi.eq.1) then
         do 70 i = 1, 3
         s = g(i,ny2-k,1)
         g(i,ny2-k,1) = 0.5*cmplx(aimag(g(i,k,1) + s),                  
     1                            real(g(i,k,1) - s))
         g(i,k,1) = 0.5*cmplx(real(g(i,k,1) + s),aimag(g(i,k,1) - s))
   70    continue
      endif
   80 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
   90 nry = nxhy/ny
      if (ks.gt.0) go to 120
      do 110 k = 2, nyh
      if (kxpi.eq.1) then
         do 100 i = 1, 3
         s = cmplx(aimag(g(i,ny2-k,1)),real(g(i,ny2-k,1)))
         g(i,ny2-k,1) = conjg(g(i,k,1) - s)
         g(i,k,1) = g(i,k,1) + s
  100    continue
      endif
  110 continue
c bit-reverse array elements in y
  120 do 140 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 140
      do 130 j = kxpi, kxpt
      t1 = g(1,k1,j)
      t2 = g(2,k1,j)
      t3 = g(3,k1,j)
      g(1,k1,j) = g(1,k,j)
      g(2,k1,j) = g(2,k,j)
      g(3,k1,j) = g(3,k,j)
      g(1,k,j) = t1
      g(2,k,j) = t2
      g(3,k,j) = t3
  130 continue
  140 continue
c first transform in y
      nry = nxy/ny
      do 180 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 150 i = kxpi, kxpt
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      t3 = s*g(3,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(3,j2,i) = g(3,j1,i) - t3
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
      g(3,j1,i) = g(3,j1,i) + t3
  150 continue
  160 continue
  170 continue
  180 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGSPOST2L(part,q,npp,noff,qm,idimp,npmax,nxv,nxyp)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c q(j,k) = charge density at grid point (j,kk),
c where kk = k + noff - 1
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nxv = first dimension of charge array, must be >= nx+1
c nxyp = actual first dimension of charge array, must be >= nxv*nypmx
      implicit none
      integer npp, noff, idimp, npmax, nxv, nxyp
      real qm
      real part, q
      dimension part(idimp,npmax), q(nxyp)
c local data
      integer mnoff, j, nnn, mmn, nn, mm, mp
      real dxn, dyn, dxp, dyp, amx, amy, dx1
      if (npp.lt.1) return
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp
c find interpolation weights
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
      mmn = mmn - mnoff
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
      subroutine PPGSJPOST2L(part,cu,edges,npp,noff,ihole,qm,dt,nx,ny,  
     1idimp,npmax,nxv,nxyp,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation, and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c also determines list of particles which are leaving this processor
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 35 flops/particle, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qm*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qm*dx*(1.-dy)
c cu(i,n,m+1)=qm*(1.-dx)*dy
c cu(i,n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = x velocity of particle n in partition
c part(4,n) = y velocity of particle n in partition
c part(5,n) = z velocity of particle n in partition
c cu(i,j,k) = ith component of current density at grid point (j,kk),
c where kk = k + noff - 1
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = first dimension of current array, must be >= nx+1
c nxyp = actual first dimension of current array, must be >= nxv*nypmx
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nxyp
      integer ipbc
      real qm, dt
      real part, cu, edges
      integer ihole
      dimension part(idimp,npmax), cu(3,nxyp)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, mmn, nnn, nop, j, nn, mm, mp, ih, nh
      real edgelx, edgely, edgerx, edgery, dxn, dyn, dxp, dyp, amx, amy
      real dx, dy, dz, vx, vy, vz, dx1, dy1
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      ih = 0
      nh = 0
      if (npp.lt.1) go to 20
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp
c find interpolation weights
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
      mmn = mmn - mnoff
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
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j - 1
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j-1) = dx
      part(2,j-1) = dy
   10 continue
      nop = npp
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
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = nop
         else
            nh = 1
         endif
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c set end of file flag
   20 if (nh.gt.0) ih = -ih
      ihole(1) = ih
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGSRJPOST2L(part,cu,edges,npp,noff,ihole,qm,dt,ci,nx, 
     1ny,idimp,npmax,nxv,nxyp,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles,
c and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c also determines list of particles which are leaving this processor
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 45 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qm*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qm*dx*(1.-dy)
c cu(i,n,m+1)=qm*(1.-dx)*dy
c cu(i,n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = x momentum of particle n in partition
c part(4,n) = y momentum of particle n in partition
c part(5,n) = z momentum of particle n in partition
c cu(i,j,k) = ith component of current density at grid point (j,kk),
c where kk = k + noff - 1
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = first dimension of current array, must be >= nx+1
c nxyp = actual first dimension of current array, must be >= nxv*nypmx
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nxyp
      integer ipbc
      real qm, dt, ci
      real part, cu, edges
      integer ihole
      dimension part(idimp,npmax), cu(3,nxyp)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, mmn, nnn, nop, j, nn, mm, mp, ih, nh
      real ci2, edgelx, edgely, edgerx, edgery, dxn, dyn, dxp, dyp
      real amx, amy, dx, dy, dz, vx, vy, vz, dx1, dy1, p2, gami
      real vxn, vyn, vzn
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
      ih = 0
      nh = 0
      if (npp.lt.1) go to 20
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
c find inverse gamma
      vxn = part(3,1)
      vyn = part(4,1)
      vzn = part(5,1)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
      gami = 1.0/sqrt(1.0 + p2*ci2)
      mmn = mmn - mnoff
      do 10 j = 2, npp
c find interpolation weights
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
      mmn = mmn - mnoff
c calculate weights
      dx = dxp*dyp
      dz = amx*dyp
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
c get momentum for next particle
      vxn = part(3,j)
      vyn = part(4,j)
      vzn = part(5,j)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
c deposit current
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
c find inverse gamma for next particle
      gami = 1.0/sqrt(1.0 + p2*ci2)
c advance position half a time-step
      dx = part(1,j-1) + vx*dt
      dy = part(2,j-1) + vy*dt
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j - 1
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j-1) = dx
      part(2,j-1) = dy
   10 continue
      nop = npp
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
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
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
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = nop
         else
            nh = 1
         endif
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c set end of file flag
   20 if (nh.gt.0) ih = -ih
      ihole(1) = ih
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGSBPUSH23L(part,fxy,bxy,edges,npp,noff,ihole,qbm,dt, 
     1dtc,ek,nx,ny,idimp,npmax,nxv,nxyp,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c also determines list of particles which are leaving this processor
c 117 flops/particle, 1 divide, 25 loads, 5 stores
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
c fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = velocity vx of particle n in partition
c part(4,n) = velocity vy of particle n in partition
c part(5,n) = velocity vz of particle n in partition
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c fxy(3,j,k) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff - 1
c bxy(1,j,k) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c edges(1:2) = lower:upper boundary of particle partition
c qbm = particle charge/mass ratio
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = first dimension of field arrays, must be >= nx+1
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nxyp
      integer ipbc
      real qbm, dt, dtc, ek
      real part, fxy, bxy, edges
      integer ihole
      dimension part(idimp,npmax), fxy(3,nxyp), bxy(3,nxyp)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, mmn, nnn, nop, nop1, j, nn, mm, mp, ih, nh
      real qtmh, edgelx, edgely, edgerx, edgery, dxn, dyn
      real dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      ih = 0
      nh = 0
      if (npp.lt.1) go to 20
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      mmn = mmn - mnoff
      nop1 = npp - 1
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
      mmn = mmn - mnoff
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
      anorm = 2.0/(1. + omt)
      omt = 0.5*(1. - omt)
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
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      nop = npp
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
      anorm = 2.0/(1. + omt)
      omt = 0.5*(1. - omt)
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
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = nop
         else
            nh = 1
         endif
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c set end of file flag
   20 if (nh.gt.0) ih = -ih
      ihole(1) = ih
c normalize kinetic energy
      ek = ek + 0.5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGSRBPUSH23L(part,fxy,bxy,edges,npp,noff,ihole,qbm,dt,
     1dtc,ci,ek,nx,ny,idimp,npmax,nxv,nxyp,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c also determines list of particles which are leaving this processor
c 129 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
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
c fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = momentum px of particle n in partition
c part(4,n) = momentum py of particle n in partition
c part(5,n) = momentum pz of particle n in partition
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c fxy(3,j,k) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff - 1
c bxy(1,j,k) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = first dimension of field arrays, must be >= nx+1
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nxyp
      integer ipbc
      real qbm, dt, dtc, ci, ek
      real part, fxy, bxy, edges
      integer ihole
      dimension part(idimp,npmax), fxy(3,nxyp), bxy(3,nxyp)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, mmn, nnn, nop, nop1, j, nn, mm, mp, ih, nh
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxn, dyn
      real dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = .5*qbm*dt
      ci2 = ci*ci
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      ih = 0
      nh = 0
      if (npp.lt.1) go to 20
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      mmn = mmn - mnoff
      nop1 = npp - 1
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
      mmn = mmn - mnoff
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp))                        
     1   + amy*(dxp*fxy(1,mm+1) + amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp))                        
     1   + amy*(dxp*fxy(2,mm+1) + amx*fxy(2,mm))
      dz = dyp*(dxp*fxy(3,mp+1) + amx*fxy(3,mp))                        
     1   + amy*(dxp*fxy(3,mm+1) + amx*fxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1) + amx*bxy(1,mp))                        
     1   + amy*(dxp*bxy(1,mm+1) + amx*bxy(1,mm))
      oy = dyp*(dxp*bxy(2,mp+1) + amx*bxy(2,mp))                        
     1   + amy*(dxp*bxy(2,mm+1) + amx*bxy(2,mm))
      oz = dyp*(dxp*bxy(3,mp+1) + amx*bxy(3,mp))                        
     1   + amy*(dxp*bxy(3,mm+1) + amx*bxy(3,mm))
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
      anorm = 2.0/(1. + omt)
      omt = 0.5*(1. - omt)
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      nop = npp
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1) + amx*bxy(1,mp))                        
     1   + amy*(dxn*bxy(1,mm+1) + amx*bxy(1,mm))
      oy = dyn*(dxn*bxy(2,mp+1) + amx*bxy(2,mp))                        
     1   + amy*(dxn*bxy(2,mm+1) + amx*bxy(2,mm))
      oz = dyn*(dxn*bxy(3,mp+1) + amx*bxy(3,mp))                        
     1   + amy*(dxn*bxy(3,mm+1) + amx*bxy(3,mm))
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
      anorm = 2.0/(1. + omt)
      omt = 0.5*(1. - omt)
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop) + dx*dtg
      dy = part(2,nop) + dy*dtg
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
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
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = nop
         else
            nh = 1
         endif
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c set end of file flag
   20 if (nh.gt.0) ih = -ih
      ihole(1) = ih
c normalize kinetic energy
      ek = ek + sum1
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
