c Fortran Library for Skeleton 2-1/2D Electromagnetic PIC Code
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
      integer mxy1, noff, moff, npp, j, k, ist, nn, mm, ierr
      real edgelx, edgely, edgerx, edgery, dx, dy
      mxy1 = mx1*my1
      ierr = 0
c sanity check
      noff = 0
      moff = 0
c loop over tiles
      do 20 k = 1, mxy1
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
      if (ist.gt.0) ierr = k
   10 continue
      noff = noff + mx
      if (noff.ge.(mx*mx1)) then
         noff = 0
         moff = moff + my
      endif
   20 continue
      if (ierr.ne.0) irc = ierr
      return
      end
c-----------------------------------------------------------------------
      subroutine POIS23T(qt,fxyt,isign,ffct,ax,ay,affp,we,nx,ny,nxvh,nyv
     1,nxhd,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions, without packed data.
c Zeros out z component.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffct
c for isign /= 0, input: qt,ffct,isign,nx,ny,nxvh,nyhd, output: fxyt,we
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
c qt(k,j) = complex charge density for fourier mode (k-1,j-1)
c fxyt(k,1,j) = x component of complex force/charge,
c fxyt(k,2,j) = y component of complex force/charge,
c fxyt(k,3,j) = zero,
c all for fourier mode (k-1,j-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffct(k,j)) = finite-size particle shape factor s
c for fourier mode (k-1,j-1)
c real(ffct(k,j)) = potential green's function g
c for fourier mode (k-1,j-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = second dimension of field arrays, must be >= nxh+1
c nyv = first dimension of field arrays, must be >= ny
c nxhd = second dimension of form factor array, must be >= nxh
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, we
      complex qt, fxyt, ffct
      dimension qt(nyv,nxvh), fxyt(nyv,3,nxvh)
      dimension ffct(nyhd,nxhd)
c local data
      integer nxh, nyh, ny2, nxh1, j, k, k1
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
      do 20 j = 1, nxh
      dkx = dnx*real(j - 1)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffct(k,j) = cmplx(affp,1.0)
      else
         ffct(k,j) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 j = 2, nxh
      dkx = dnx*real(j - 1)
      do 40 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffct(k,j))*aimag(ffct(k,j))
      at2 = dkx*at1
      at3 = dny*real(k - 1)*at1
      zt1 = cmplx(aimag(qt(k,j)),-real(qt(k,j)))
      zt2 = cmplx(aimag(qt(k1,j)),-real(qt(k1,j)))
      fxyt(k,1,j) = at2*zt1
      fxyt(k,2,j) = at3*zt1
      fxyt(k,3,j) = zero
      fxyt(k1,1,j) = at2*zt2
      fxyt(k1,2,j) = -at3*zt2
      fxyt(k1,3,j) = zero
      wp = wp + at1*(qt(k,j)*conjg(qt(k,j)) + qt(k1,j)*conjg(qt(k1,j)))
   40 continue
   50 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nxh
      at1 = real(ffct(1,j))*aimag(ffct(1,j))
      at2 = dnx*real(j - 1)*at1
      zt1 = cmplx(aimag(qt(1,j)),-real(qt(1,j)))
      fxyt(1,1,j) = at2*zt1
      fxyt(1,2,j) = zero
      fxyt(1,3,j) = zero
      fxyt(k1,1,j) = zero
      fxyt(k1,2,j) = zero
      fxyt(k1,3,j) = zero
      wp = wp + at1*(qt(1,j)*conjg(qt(1,j)))
   60 continue
c mode numbers kx = 0, nx/2
      nxh1 = nxh + 1
cdir$ ivdep
      do 70 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffct(k,1))*aimag(ffct(k,1))
      at3 = dny*real(k - 1)*at1
      zt1 = cmplx(aimag(qt(k,1)),-real(qt(k,1)))
      fxyt(k,1,1) = zero
      fxyt(k,2,1) = at3*zt1
      fxyt(k,3,1) = zero
      fxyt(k1,1,1) = zero
      fxyt(k1,2,1) = at3*conjg(zt1)
      fxyt(k1,3,1) = zero
      fxyt(k,1,nxh1) = zero
      fxyt(k,2,nxh1) = zero
      fxyt(k,3,nxh1) = zero
      fxyt(k1,1,nxh1) = zero
      fxyt(k1,2,nxh1) = zero
      fxyt(k1,3,nxh1) = zero
      wp = wp + at1*(qt(k,1)*conjg(qt(k,1)))
   70 continue
      k1 = nyh + 1
      fxyt(1,1,1) = zero
      fxyt(1,2,1) = zero
      fxyt(1,3,1) = zero
      fxyt(k1,1,1) = zero
      fxyt(k1,2,1) = zero
      fxyt(k1,3,1) = zero
      fxyt(1,1,nxh1) = zero
      fxyt(1,2,nxh1) = zero
      fxyt(1,3,nxh1) = zero
      fxyt(k1,1,nxh1) = zero
      fxyt(k1,2,nxh1) = zero
      fxyt(k1,3,nxh1) = zero
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
