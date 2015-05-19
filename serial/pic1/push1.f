c Fortran Library for Skeleton 1D Electrostatic PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR1(part,vtx,vdx,npx,idimp,nop,nx,ipbc)
c for 1d code, this subroutine calculates initial particle co-ordinate
c and velocity, with uniform density and maxwellian velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c vtx = thermal velocity of particles in x direction
c vdx = drift velocity of particles x direction
c npx = number of particles distributed in x direction
c idimp = size of phase space = 2
c nop = number of particles
c nx = system length in x direction
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, idimp, nop, nx, ipbc
      real part, vtx, vdx
      dimension part(idimp,nop)
c local data
      integer j
      real edgelx, at1, sum1
      double precision dsum1
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
   20 continue
c add correct drift
      dsum1 = 0.0d0
      do 30 j = 1, npx
      dsum1 = dsum1 + part(2,j)
   30 continue
      sum1 = dsum1
      sum1 = sum1/real(npx) - vdx
      do 40 j = 1, npx
      part(2,j) = part(2,j) - sum1
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
c for 1d code, this subroutine updates particle co-ordinate and velocity
c using leap-frog scheme in time and first-order linear interpolation
c in space, with various boundary conditions.
c scalar version using guard cells
c 16 flops/particle, 4 loads, 2 stores
c input: all, output: part, ek
c equations used are:
c v(t+dt/2) = v(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
c and x(t+dt) = x(t) + v(t+dt/2)*dt
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c fx(j) = force/charge at grid point j, that is convolution of electric
c field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((v(t+dt/2)+v(t-dt/2))**2)
c idimp = size of phase space = 2
c nop = number of particles
c nx = system length in x direction
c nxv = first dimension of charge array, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fx, qbm, dt, ek
      dimension part(idimp,nop), fx(nxv)
c local data
      integer j, nn
      real qtm, edgelx, edgerx, dx
      double precision sum1
      qtm = qbm*dt
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
      dx = part(1,j) - real(nn)
      nn = nn + 1
c find acceleration
      dx = (1.0 - dx)*fx(nn) + dx*fx(nn+1)
c new velocity
      dx = part(2,j) + qtm*dx
c average kinetic energy
      sum1 = sum1 + (part(2,j) + dx)**2
      part(2,j) = dx
c new position
      dx = part(1,j) + dx*dt
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
      ek = ek + .125*sum1
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
      subroutine CGUARD1L(fx,nx,nxe)
c replicate extended periodic field fx
c linear interpolation
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      real fx
      integer nx, nxe
      dimension fx(nxe)
c copy edge of extended field
      fx(nx+1) = fx(1)
      return
      end
c-----------------------------------------------------------------------
      subroutine AGUARD1L(q,nx,nxe)
c accumulate extended periodic field
c linear interpolation
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      real q
      integer nx, nxe
      dimension q(nxe)
c accumulate edge of extended field
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
      ani = 1./real(2*nx)
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
