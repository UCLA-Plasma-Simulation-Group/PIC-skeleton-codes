c Fortran Library for Skeleton 2D Electrostatic PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR2(part,vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx,ny,   
     1ipbc)
c for 2d code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c vtx/vty = thermal velocity of electrons in x/y direction
c vdx/vdy = drift velocity of beam electrons in x/y direction
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, idimp, nop, nx, ny, ipbc
      real vtx, vty, vdx, vdy
      real part
      dimension part(idimp,nop)
c local data
      integer j, k, k1, npxy
      real edgelx, edgely, at1, at2, at3, sum1, sum2
      double precision dsum1, dsum2
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
   30 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      do 40 j = 1, npxy
      dsum1 = dsum1 + part(3,j)
      dsum2 = dsum2 + part(4,j)
   40 continue
      sum1 = dsum1
      sum2 = dsum2
      at1 = 1.0/real(npxy)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      do 50 j = 1, npxy
      part(3,j) = part(3,j) - sum1
      part(4,j) = part(4,j) - sum2
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,ipbc
     1)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions.
c scalar version using guard cells
c 44 flops/particle, 12 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = second dimension of field arrays, must be >= nx+1
c nyv = third dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real qbm, dt, ek
      real part, fxy
      dimension part(idimp,nop), fxy(2,nxv,nyv)
c local data
      integer j, nn, mm, np, mp
      real qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy
      double precision sum1
      qtm = qbm*dt
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
c find acceleration
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp))                    
     1   + amy*(dxp*fxy(1,np,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp))                    
     1   + amy*(dxp*fxy(2,np,mm) + amx*fxy(2,nn,mm))
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
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
      ek = ek + 0.125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine DGPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,   
     1ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions.
c scalar version using guard cells
c 44 flops/particle, 12 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      double precision qbm, dt, ek
      double precision part, fxy
      dimension part(idimp,nop), fxy(2,nxv,nyv)
c local data
      integer j, nn, mm, np, mp
      double precision qtm, edgelx, edgely, edgerx, edgery
      double precision dxp, dyp, amx, amy, dx, dy
      double precision sum1
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0d0
      edgely = 0.0d0
      edgerx = dble(nx)
      edgery = dble(ny)
      if (ipbc.eq.2) then
         edgelx = 1.0d0
         edgely = 1.0d0
         edgerx = dble(nx-1)
         edgery = dble(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0d0
         edgerx = dble(nx-1)
      endif
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - dble(nn)
      dyp = part(2,j) - dble(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1.0d0 - dxp
      mp = mm + 1
      amy = 1.0d0 - dyp
      np = nn + 1
c find acceleration
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp))                    
     1   + amy*(dxp*fxy(1,np,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp))                    
     1   + amy*(dxp*fxy(2,np,mm) + amx*fxy(2,nn,mm))
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
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
      ek = ek + 0.125d0*sum1
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
c deposit charge
      q(np,mp) = q(np,mp) + dxp*dyp
      q(nn,mp) = q(nn,mp) + amx*dyp
      q(np,mm) = q(np,mm) + dxp*amy
      q(nn,mm) = q(nn,mm) + amx*amy
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
      subroutine CGUARD2L(fxy,nx,ny,nxe,nye)
c replicate extended periodic vector field fxy
c linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
      implicit none
      real fxy
      integer nx, ny, nxe, nye
      dimension fxy(2,nxe,nye)
c local data
      integer j, k
c copy edges of extended field
      do 10 k = 1, ny
      fxy(1,nx+1,k) = fxy(1,1,k)
      fxy(2,nx+1,k) = fxy(2,1,k)
   10 continue
      do 20 j = 1, nx
      fxy(1,j,ny+1) = fxy(1,j,1)
      fxy(2,j,ny+1) = fxy(2,j,1)
   20 continue
      fxy(1,nx+1,ny+1) = fxy(1,1,1)
      fxy(2,nx+1,ny+1) = fxy(2,1,1)
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
      subroutine POIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,   
     1nxhd,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,nxvh,nyhd, output: fxy,we
c approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
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
      dimension q(nxvh,nyv), fxy(2,nxvh,nyv)
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
      if (at3.eq.0.) then
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
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
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
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
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
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   70 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      we = real(nx*ny)*wp
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
      subroutine WFFT2R2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,    
     1nxyhd)
c wrapper function for 2 2d real to complex ffts
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2R2X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
c perform y fft
         call FFT2R2Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2R2Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c perform x fft
         call FFT2R2X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
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
      subroutine FFT2R2X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd,  
     1nxhyd,nxyhd)
c this subroutine performs the x part of 2 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c y, using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, two inverse fourier transforms are performed
c f(1:2,n,m) = (1/nx*ny)*sum(f(1:2,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
c f(1:2,j,k) = sum(f(1:2,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = second dimension of f
c nyd = third dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(1:2,j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
c except for f(1:2,1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c aimag(f(1:2,1,1)) = real part of mode nx/2,0 and
c aimag(f(1:2,1,ny/2+1)) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      real at1, ani
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
      do 20 k = nyi, nyt
      do 10 j = 1, nxh
      at1 = aimag(f(1,j,k))
      f(1,j,k) = cmplx(real(f(1,j,k)),real(f(2,j,k)))
      f(2,j,k) = cmplx(at1,aimag(f(2,j,k)))
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
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
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
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
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
      do 90 jj = 1, 2
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
      do 120 jj = 1, 2
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
      do 150 jj = 1, 2
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,nxh2-j,k) = conjg(t1 - t2)
  150 continue
  160 continue
  170 continue
      do 190 k = nyi, nyt
      do 180 jj = 1, 2
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
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
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
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
  220 continue
  230 continue
  240 continue
  250 continue
c swap complex components
      do 270 k = nyi, nyt
      do 260 j = 1, nxh
      at1 = aimag(f(1,j,k))
      f(1,j,k) = cmplx(real(f(1,j,k)),real(f(2,j,k)))
      f(2,j,k) = cmplx(at1,aimag(f(2,j,k)))
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2R2Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd,  
     1nxhyd,nxyhd)
c this subroutine performs the y part of 2 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c x, using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, two inverse fourier transforms are performed
c f(1:2,n,m) = (1/nx*ny)*sum(f(1:2,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
c f(1:2,j,k) = sum(f(1:2,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = second dimension of f
c nyd = third dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(1:2,j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
c except for f(1:2,1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c aimag(f(1:2,1,1)) = real part of mode nx/2,0 and
c aimag(f(1:2,1,ny/2+1)) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      complex t1, t2, t3
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
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
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
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      do 80 k = 2, nyh
      if (nxi.eq.1) then
         do 70 jj = 1, 2
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
         do 100 jj = 1, 2
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
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
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
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
  140 continue
  150 continue
  160 continue
  170 continue
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
      subroutine GSPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxyv,  
     1ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 44 flops/particle, 12 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1.m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1.m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = second virtual dimension of field arrays, must be >= nx+1
c nxyv = dimension of field arrays, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real qbm, dt, ek
      real part, fxy
      dimension part(idimp,nop), fxy(2,nxyv)
c local data
      integer j, nnn, mmn, nop1, nn, mm, mp
      real dxn, dyn, qtm, edgelx, edgely, edgerx, edgery
      real dxp, dyp, amx, amy, dx, dy
      double precision sum1
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      nop1 = nop - 1
      qtm = qbm*dt
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
c find acceleration
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp))                        
     1   + amy*(dxp*fxy(1,mm+1) + amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp))                        
     1   + amy*(dxp*fxy(2,mm+1) + amx*fxy(2,mm))
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
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
c find acceleration
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp))                        
     1   + amy*(dxn*fxy(1,mm+1) + amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp))                        
     1   + amy*(dxn*fxy(2,mm+1) + amx*fxy(2,mm))
c new velocity
      dx = part(3,nop) + qtm*dx
      dy = part(4,nop) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,nop))**2 + (dy + part(4,nop))**2
      part(3,nop) = dx
      part(4,nop) = dy
c new position
      dx = part(1,nop) + dx*dt
      dy = part(2,nop) + dy*dt
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
   20 ek = ek + 0.125*sum1
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
