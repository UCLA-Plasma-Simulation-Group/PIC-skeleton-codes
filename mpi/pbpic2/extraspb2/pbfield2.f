c Fortran library for Skeleton 2-1/2D Electromagnetic MPI PIC Code field
c diagnostics
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c-----------------------------------------------------------------------
      subroutine PPOTP2(q,pot,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c potential, with periodic boundary conditions for distributed data.
c input: q,ffc,nx,ny,kstrt,nyv,kxp,nyhd, output: pot,we
c approximate flop count is: 21*nxc*nyc + 11*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c potential is calculated using the equation:
c pot(kx,ky) = g(kx,ky)*q(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c pot(kx=pi) = pot(ky=pi) = 0, and pot(kx=0,ky=0) = 0.
c q(k,j) = complex charge density for fourier mode (jj-1,k-1)
c pot(k,j) = x component of complex potential
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffc(k,j)) = finite-size particle shape factor s
c real(ffc(k,j)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c where affp = normalization constant = nx*ny/np,
c where np=number of particles
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real we
      complex q, pot, ffc
      dimension q(nyv,kxp), pot(nyv,kxp)
      dimension ffc(nyhd,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real at1, at2
      complex zero
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      zero = cmplx(0.0,0.0)
c calculate potential and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 40
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 1, kxps
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at2 = real(ffc(k,j))
         at1 = at2*aimag(ffc(k,j))
         pot(k,j) = at2*q(k,j)
         pot(k1,j) = at2*q(k1,j)
         wp = wp + at1*(q(k,j)*conjg(q(k,j)) + q(k1,j)*conjg(q(k1,j)))
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = real(ffc(1,j))
         at1 = at2*aimag(ffc(1,j))
         pot(1,j) = at2*q(1,j)
         pot(k1,j) = zero
         wp = wp + at1*(q(1,j)*conjg(q(1,j)))
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         at2 = real(ffc(k,1))
         at1 = at2*aimag(ffc(k,1))
         pot(k,1) = at2*q(k,1)
         pot(k1,1) = zero
         wp = wp + at1*(q(k,1)*conjg(q(k,1)))
   30    continue
         k1 = nyh + 1
         pot(1,1) = zero
         pot(k1,1) = zero
      endif
   40 continue
      we = real(nx)*real(ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp)
c this subroutine calculates the divergence in fourier space
c input: all except df, output: df
c approximate flop count is: 16*nxc*nyc + 5*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the divergence is calculated using the equation:
c df(kx,ky) = sqrt(-1)*(kx*fx(kx,ky)+ky*fy(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for df(kx=pi) = df(ky=pi) = df(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp = number of data values per block
      implicit none
      integer nx, ny, kstrt, ndim, nyv, kxp
      complex f, df
      dimension f(ndim,nyv,kxp), df(nyv,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky
      complex zero, zt1
      if (ndim.lt.2) return
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
c calculate the divergence
      if (kstrt.gt.nxh) return
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 1, kxps
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         zt1 = dkx*f(1,k,j) + dky*f(2,k,j)
         df(k,j) = cmplx(-aimag(zt1),real(zt1))
         zt1 = dkx*f(1,k1,j) - dky*f(2,k1,j)
         df(k1,j) = cmplx(-aimag(zt1),real(zt1))
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         df(1,j) = dkx*cmplx(-aimag(f(1,1,j)),real(f(1,1,j)))
         df(k1,j) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         df(k,1) = dky*cmplx(-aimag(f(2,k,1)),real(f(2,k,1)))
         df(k1,1) = zero
   30    continue
         k1 = nyh + 1
         df(1,1) = zero
         df(k1,1) = zero
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp)
c this subroutine calculates the gradient in fourier space
c input: all except f, output: f
c approximate flop count is: 12*nxc*nyc + 4*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the gradient is calculated using the equations:
c fx(kx,ky) = sqrt(-1)*kx*df(kx,ky)
c fy(kx,ky) = sqrt(-1)*ky*df(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for fx(kx=pi) = fy(kx=pi) = 0, fx(ky=pi) = fy(ky=pi) = 0,
c and fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp = number of data values per block
      implicit none
      integer nx, ny, kstrt, ndim, nyv, kxp
      complex df, f
      dimension df(nyv,kxp), f(ndim,nyv,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky
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
c calculate the gradient
      if (kstrt.gt.nxh) return
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 1, kxps
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         zt1 = cmplx(-aimag(df(k,j)),real(df(k,j)))
         f(1,k,j) = dkx*zt1
         f(2,k,j) = dky*zt1
         zt1 = cmplx(-aimag(df(k1,j)),real(df(k1,j)))
         f(1,k1,j) = dkx*zt1
         f(2,k1,j) = -dky*zt1
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         f(1,1,j) = dkx*cmplx(-aimag(df(1,j)),real(df(1,j)))
         f(2,1,j) = zero
         f(1,k1,j) = zero
         f(2,k1,j) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         f(1,k,1) = zero
         f(2,k,1) = dky*cmplx(-aimag(df(k,1)),real(df(k,1)))
         f(1,k1,1) = zero
         f(2,k1,1) = zero
   30    continue
         k1 = nyh + 1
         f(1,1,1) = zero
         f(2,1,1) = zero
         f(1,k1,1) = zero
         f(2,k1,1) = zero
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPCURLF2(f,g,nx,ny,kstrt,nyv,kxp)
c this subroutine calculates the curl in fourier space
c input: all except g, output: g
c approximate flop count is: 32*nxc*nyc + 10*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the curl is calculated using the equations:
c gx(kx,ky) = sqrt(-1)*ky*fz(kx,ky)
c gy(kx,ky) = -sqrt(-1)*kx*fz(kx,ky)
c gz(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
c and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = second dimension of field arrays, must be >= ny
c kxp = number of data values per block
      implicit none
      integer nx, ny, kstrt, nyv, kxp
      complex f, g
      dimension f(3,nyv,kxp), g(3,nyv,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky
      complex zero, zt1, zt2, zt3
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
c calculate the curl
      if (kstrt.gt.nxh) return
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 1, kxps
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         zt1 = cmplx(-aimag(f(3,k,j)),real(f(3,k,j)))
         zt2 = cmplx(-aimag(f(2,k,j)),real(f(2,k,j)))
         zt3 = cmplx(-aimag(f(1,k,j)),real(f(1,k,j)))
         g(1,k,j) = dky*zt1
         g(2,k,j) = -dkx*zt1
         g(3,k,j) = dkx*zt2 - dky*zt3
         zt1 = cmplx(-aimag(f(3,k1,j)),real(f(3,k1,j)))
         zt2 = cmplx(-aimag(f(2,k1,j)),real(f(2,k1,j)))
         zt3 = cmplx(-aimag(f(1,k1,j)),real(f(1,k1,j)))
         g(1,k1,j) = -dky*zt1
         g(2,k1,j) = -dkx*zt1
         g(3,k1,j) = dkx*zt2 + dky*zt3
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         zt1 = cmplx(-aimag(f(3,1,j)),real(f(3,1,j)))
         zt2 = cmplx(-aimag(f(2,1,j)),real(f(2,1,j)))
         g(1,1,j) = zero
         g(2,1,j) = -dkx*zt1
         g(3,1,j) = dkx*zt2
         g(1,k1,j) = zero
         g(2,k1,j) = zero
         g(3,k1,j) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         zt1 = cmplx(-aimag(f(3,k,1)),real(f(3,k,1)))
         zt2 = cmplx(-aimag(f(1,k,1)),real(f(1,k,1)))
         g(1,k,1) = dky*zt1
         g(2,k,1) = zero
         g(3,k,1) = -dky*zt2
         g(1,k1,1) = zero
         g(2,k1,1) = zero
         g(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         g(1,1,1) = zero
         g(2,1,1) = zero
         g(3,1,1) = zero
         g(1,k1,1) = zero
         g(2,k1,1) = zero
         g(3,k1,1) = zero
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPAVPOT23(bxy,axy,nx,ny,kstrt,nyv,kxp)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with periodic boundary conditions
c for distributed data.
c input: bxy,nx,ny,kstrt,nyv,kxp, output: axy
c approximate flop count is: 38*nxc*nyc + 10*(nxc + nyc)
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0,
c ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c bxy(i,k,j) = i-th component of complex magnetic field,
c axy(i,k,j) = i-th component of complex vector potential,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c kxp = number of data values per block
c kstrt = starting data block number
c nx/ny = system length in x/y direction
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp
      complex bxy, axy
      dimension bxy(3,nyv,kxp), axy(3,nyv,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky, dkx2, at1, at2, at3
      complex zero, zt1, zt2, zt3
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
c calculate vector potential
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
         at2 = dky*at1
         at3 = dkx*at1
         zt1 = cmplx(-aimag(bxy(3,k,j)),real(bxy(3,k,j)))
         zt2 = cmplx(-aimag(bxy(2,k,j)),real(bxy(2,k,j)))
         zt3 = cmplx(-aimag(bxy(1,k,j)),real(bxy(1,k,j)))
         axy(1,k,j) = at2*zt1
         axy(2,k,j) = -at3*zt1
         axy(3,k,j) = at3*zt2 - at2*zt3
         zt1 = cmplx(-aimag(bxy(3,k1,j)),real(bxy(3,k1,j)))
         zt2 = cmplx(-aimag(bxy(2,k1,j)),real(bxy(2,k1,j)))
         zt3 = cmplx(-aimag(bxy(1,k1,j)),real(bxy(1,k1,j)))
         axy(1,k1,j) = -at2*zt1
         axy(2,k1,j) = -at3*zt1
         axy(3,k1,j) = at3*zt2 + at2*zt3
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = 1.0/dkx
         zt1 = cmplx(-aimag(bxy(3,1,j)),real(bxy(3,1,j)))
         zt2 = cmplx(-aimag(bxy(2,1,j)),real(bxy(2,1,j)))
         axy(1,1,j) = zero
         axy(2,1,j) = -at2*zt1
         axy(3,1,j) = at2*zt2
         axy(1,k1,j) = zero
         axy(2,k1,j) = zero
         axy(3,k1,j) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at2 = 1.0/dky
         zt1 = cmplx(-aimag(bxy(3,k,1)),real(bxy(3,k,1)))
         zt2 = cmplx(-aimag(bxy(1,k,1)),real(bxy(1,k,1)))
         axy(1,k,1) = at2*zt1
         axy(2,k,1) = zero
         axy(3,k,1) = -at2*zt2
         axy(1,k1,1) = zero
         axy(2,k1,1) = zero
         axy(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         axy(1,1,1) = zero
         axy(2,1,1) = zero
         axy(3,1,1) = zero
         axy(1,k1,1) = zero
         axy(2,k1,1) = zero
         axy(3,k1,1) = zero
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPAVRPOT23(axy,bxy,ffc,affp,ci,nx,ny,kstrt,nyv,kxp,    
     1nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c the radiative part of the vector potential
c with periodic boundary conditions, for distributed data.
c input: all, output: axy
c approximate flop count is: 68*nxc*nyc + 20*(nxc + nyc)
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the radiative vector potential is updated using the equations:
c ax(kx,ky) = (sqrt(-1)*ky*bz(kx,ky)
c                       - affp*ci2*cux(kx,ky)*s(kx,ky)/(kx*kx+ky*ky)
c ay(kx,ky) = -(sqrt(-1)*kx*bz(kx,ky)
c                       + affp*ci2*cuy(kx,ky)*s(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = (sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*ci2*cuz(kx,ky)*s(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, ci2 = ci*ci
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0,
c ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c axy(i,k,j) = on entry, i-th component of complex current density cu,
c axy(i,k,j) = on exit, i-th component of complex radiative vector
c potential,
c bxy(i,k,j) = i-th component of complex magnetic field,
c aimag(ffc(k,j)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = second dimension of field arrays, must be >= ny
c kxp = number of data values per block
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real affp, ci
      complex axy, bxy, ffc
      dimension axy(3,nyv,kxp), bxy(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, afc2, dkx, dkx2, dky, at1, at2
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      afc2 = affp*ci*ci
      zero = cmplx(0.,0.)
c calculate the radiative vector potential
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
         at2 = afc2*aimag(ffc(k,j))
c update radiative vector potential, ky > 0
         zt1 = cmplx(-aimag(bxy(3,k,j)),real(bxy(3,k,j)))
         zt2 = cmplx(-aimag(bxy(2,k,j)),real(bxy(2,k,j)))
         zt3 = cmplx(-aimag(bxy(1,k,j)),real(bxy(1,k,j)))
         axy(1,k,j) = at1*(dky*zt1 - at2*axy(1,k,j))
         axy(2,k,j) = -at1*(dkx*zt1 + at2*axy(2,k,j))
         axy(3,k,j) = at1*((dkx*zt2 - dky*zt3) - at2*axy(3,k,j))
c update radiative vector potential, ky < 0
         zt1 = cmplx(-aimag(bxy(3,k1,j)),real(bxy(3,k1,j)))
         zt2 = cmplx(-aimag(bxy(2,k1,j)),real(bxy(2,k1,j)))
         zt3 = cmplx(-aimag(bxy(1,k1,j)),real(bxy(1,k1,j)))
         axy(1,k1,j) = -at1*(dky*zt1 + at2*axy(1,k1,j))
         axy(2,k1,j) = -at1*(dkx*zt1 + at2*axy(2,k1,j))
         axy(3,k1,j) = at1*((dkx*zt2 + dky*zt3) - at2*axy(3,k1,j))
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = 1.0/dkx2
         at2 = afc2*aimag(ffc(1,j))
c update radiative vector potential
         zt1 = cmplx(-aimag(bxy(3,1,j)),real(bxy(3,1,j)))
         zt2 = cmplx(-aimag(bxy(2,1,j)),real(bxy(2,1,j)))
         axy(1,1,j) = zero
         axy(2,1,j) = -at1*(dkx*zt1 + at2*axy(2,1,j))
         axy(3,1,j) = at1*(dkx*zt2 - at2*axy(3,1,j))
         axy(1,k1,j) = zero
         axy(2,k1,j) = zero
         axy(3,k1,j) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*real(k - 1)
         at1 = 1.0/(dky*dky)
         at2 = afc2*aimag(ffc(k,1))
c update radiative vector potential
         zt1 = cmplx(-aimag(bxy(3,k,1)),real(bxy(3,k,1)))
         zt3 = cmplx(-aimag(bxy(1,k,1)),real(bxy(1,k,1)))
         axy(1,k,1) = at1*(dky*zt1 - at2*axy(1,k,1))
         axy(2,k,1) = zero
         axy(3,k,1) = -at1*(dky*zt3 + at2*axy(3,k,1))
         axy(1,k1,1) = zero
         axy(2,k1,1) = zero
         axy(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         axy(1,1,1) = zero
         axy(2,1,1) = zero
         axy(3,1,1) = zero
         axy(1,k1,1) = zero
         axy(2,k1,1) = zero
         axy(3,k1,1) = zero
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPSMOOTH2(q,qs,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
c this subroutine provides a 2d scalar smoothing function
c in fourier space, with periodic boundary conditions
c for distributed data.
c input: q,ffc,nx,ny,kstrt,nyv,kxp,nyhd, output: qs
c approximate flop count is: 4*nxc*nyc + 2*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c smoothing is calculated using the equation:
c qs(kx,ky) = q(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c qs(kx=pi) = qs(ky=pi) = 0, and qs(kx=0,ky=0) = 0.
c q(k,j) = complex charge density
c qs(k,j) = complex smoothed charge density,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffc(k,j)) = finite-size particle shape factor s
c real(ffc(k,j)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      complex q, qs, ffc
      dimension q(nyv,kxp), qs(nyv,kxp)
      dimension ffc(nyhd,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real at1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      zero = cmplx(0.0,0.0)
c calculate smoothing
      if (kstrt.gt.nxh) return
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 1, kxps
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j))
         qs(k,j) = at1*q(k,j)
         qs(k1,j) = at1*q(k1,j)
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j))
         qs(1,j) = at1*q(1,j)
         qs(k1,j) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,1))
         qs(k,1) = at1*q(k,1)
         qs(k1,1) = zero
   30    continue
         k1 = nyh + 1
         qs(1,1) = cmplx(aimag(ffc(1,1))*real(q(1,1)),0.)
         qs(k1,1) = zero
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPSMOOTH23(cu,cus,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
c this subroutine provides a 2d vector smoothing function
c in fourier space, with periodic boundary conditions.
c for distributed data.
c input: cu,ffc,nx,ny,kstrt,nyv,kxp,nyhd, output: cus
c approximate flop count is: 12*nxc*nyc + 6*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c smoothing is calculated using the equation:
c cusx(kx,ky) = cux(kx,ky)*s(kx,ky)
c cusy(kx,ky) = cuy(kx,ky)*s(kx,ky)
c cusz(kx,ky) = cuz(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c cusx(kx=pi) = cusy(kx=pi) = cusz(kx=pi) = 0,
c cusx(ky=pi) = cusy(ky=pi) = cusz(ky=pi) = 0,
c cusx(kx=0,ky=0) = cusy(kx=0,ky=0) = cusz(kx=0,ky=0) = 0.
c cu(i,k,j) = i-th component of complex current density and
c cus(i,k,j) = i-th component of complex smoothed current density
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffc(k,j)) = finite-size particle shape factor s
c real(ffc(k,j)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c nx/ny = system length in x/y direction
c nyv = second dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      complex cu, cus, ffc
      dimension cu(3,nyv,kxp), cus(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real at1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      zero = cmplx(0.0,0.0)
c calculate smoothing
      if (kstrt.gt.nxh) return
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 1, kxps
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j))
         cus(1,k,j) = at1*cu(1,k,j)
         cus(2,k,j) = at1*cu(2,k,j)
         cus(3,k,j) = at1*cu(3,k,j)
         cus(1,k1,j) = at1*cu(1,k1,j)
         cus(2,k1,j) = at1*cu(2,k1,j)
         cus(3,k1,j) = at1*cu(3,k1,j)
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j))
         cus(1,1,j) = at1*cu(1,1,j)
         cus(2,1,j) = at1*cu(2,1,j)
         cus(3,1,j) = at1*cu(3,1,j)
         cus(1,k1,j) = zero
         cus(2,k1,j) = zero
         cus(3,k1,j) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,1))
         cus(1,k,1) = at1*cu(1,k,1)
         cus(2,k,1) = at1*cu(2,k,1)
         cus(3,k,1) = at1*cu(3,k,1)
         cus(1,k1,1) = zero
         cus(2,k1,1) = zero
         cus(3,k1,1) = zero
   30    continue
         k1 = nyh + 1
         at1 = aimag(ffc(1,1))
         cus(1,1,1) = cmplx(at1*real(cu(1,1,1)),0.)
         cus(2,1,1) = cmplx(at1*real(cu(2,1,1)),0.)
         cus(3,1,1) = cmplx(at1*real(cu(3,1,1)),0.)
         cus(1,k1,1) = zero
         cus(2,k1,1) = zero
         cus(3,k1,1) = zero
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPRDMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,kxp, 
     1modesxpd,modesyd)
c this subroutine extracts lowest order modes from packed complex array
c pot and stores them into a location in an unpacked complex array pott
c modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
c and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c except kx = NX/2 is stored at location kxp+1 when idproc=0.
c nx/ny = system length in x/y direction
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c kstrt = starting data block number
c nyv = first dimension of input array pot, nyv >= ny
c kxp = number of data values per block
c modesyd = second dimension of array pott,
c where modesyd >= min(2*modesy-1,ny)
c modesxpd = third dimension of array pott, modesxpd >= min(modesx,kxp)
c unless modesx = nx/2+1, in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, modesx, modesy, kstrt, nyv, kxp
      integer modesxpd, modesyd
      complex pot, pott
      dimension pot(nyv,kxp), pott(modesyd,modesxpd)
c local data
      integer nxh, nyh, jmax, kmax, ny2, j, k, j1, k1, ks, joff
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 1
      if (kstrt.gt.nxh) return
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*ks
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
      do 20 j = 1, jmax
      if ((j+joff).gt.1) then
         do 10 k = 2, kmax
         k1 = ny2 - k
         pott(2*k-2,j) = pot(k,j)
         pott(2*k-1,j) = pot(k1,j)
   10    continue
c mode numbers ky = 0, ny/2
         pott(1,j) = pot(1,j)
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            pott(ny,j) = pot(k1,j)
         endif
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 30 k = 2, kmax
         k1 = ny2 - k
         pott(2*k-2,1) = pot(k,1)
         pott(2*k-1,1) = conjg(pot(k,1))
         if (modesx.gt.nxh) then
            pott(2*k-2,j1) = conjg(pot(k1,1))
            pott(2*k-1,j1) = pot(k1,1)
         endif
   30    continue
         pott(1,1) = cmplx(real(pot(1,1)),0.0)
         if (modesx.gt.nxh) then
            pott(1,j1) = cmplx(aimag(pot(1,1)),0.0)
         endif
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            pott(ny,1) = cmplx(real(pot(k1,1)),0.0)
            if (modesx.gt.nxh) then
               pott(ny,j1) = cmplx(aimag(pot(k1,1)),0.0)
            endif
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPWRMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,kxp, 
     1modesxpd,modesyd)
c this subroutine extracts lowest order modes from a location in an
c unpacked complex array pott and stores them into a packed complex
c array pot
c modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
c and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c except kx = NX/2 is stored at location kxp+1 when idproc=0.
c nx/ny = system length in x/y direction
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c kstrt = starting data block number
c nyv = first dimension of input array pot, nyv >= ny
c kxp = number of data values per block
c modesyd = second dimension of array pott,
c where modesyd  >= min(2*modesy-1,ny)
c modesxpd = third dimension of array pott, modesxpd >= min(modesx,kxp)
c unless modesx = nx/2+1, in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, modesx, modesy, kstrt, nyv, kxp
      integer modesxpd, modesyd
      complex pot, pott
      dimension pot(nyv,kxp), pott(modesyd,modesxpd)
c local data
      integer nxh, nyh, jmax, kmax, ny2, j, k, j1, k1, ks, joff
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 1
      zero = cmplx(0.0,0.0)
      if (kstrt.gt.nxh) return
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*ks
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
      do 30 j = 1, jmax
      if ((j+joff).gt.1) then
         do 10 k = 2, kmax
         k1 = ny2 - k
         pot(k,j) = pott(2*k-2,j)
         pot(k1,j) = pott(2*k-1,j)
   10    continue
         do 20 k = kmax+1, nyh
         k1 = ny2 - k
         pot(k,j) = zero
         pot(k1,j) = zero
   20    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         pot(1,j) = pott(1,j)
         pot(k1,j) = zero
         if (modesy.gt.nyh) then
            pot(k1,j) = pott(ny,j)
         endif
      endif
   30 continue
      do 50 j = jmax+1, kxp
      if ((j+joff).gt.1) then
         do 40 k = 2, nyh
         k1 = ny2 - k
         pot(k,j) = zero
         pot(k1,j) = zero
   40    continue
         k1 = nyh + 1
c mode numbers ky = 0, ny/2
         pot(1,j) = zero
         pot(k1,j) = zero
      endif
   50 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 60 k = 2, kmax
         k1 = ny2 - k
         pot(k,1) = pott(2*k-2,1)
         pot(k1,1) = zero
         if (modesx.gt.nxh) then
            pot(k1,1) = conjg(pott(2*k-2,j1))
         endif
   60    continue
         do 70 k = kmax+1, nyh
         k1 = ny2 - k
         pot(k,1) = zero
         pot(k1,1) = zero
   70    continue
         k1 = nyh + 1
         pot(1,1) = cmplx(real(pott(1,1)),0.0)
         pot(k1,1) = zero
         if (modesx.gt.nxh) then
            pot(1,1) = cmplx(real(pot(1,1)),real(pott(1,j1)))
         endif
         if (modesy.gt.nyh) then
            pot(k1,1) = cmplx(real(pott(ny,1)),0.0)
            if (modesx.gt.nxh) then
               pot(k1,1) = cmplx(real(pot(k1,1)),real(pott(ny,j1)))
            endif
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPRDVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,kstrt, 
     1nyv,kxp,modesxpd,modesyd)
c this subroutine extracts lowest order modes from packed complex vector
c array vpot and stores them into a location in an unpacked complex
c vector array vpott
c modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
c and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c except kx = NX/2 is stored at location kxp+1 when idproc=0.
c nx/ny = system length in x/y direction
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c ndim = number of field arrays, must be >= 1
c kstrt = starting data block number
c nyv = second dimension of input array vpot, nyv >= ny
c kxp = number of data values per block
c modesyd = third dimension of array vpott,
c where modesyd >= min(2*modesy-1,ny)
c modesxpd = fourth dimension of array vpott,
c modesxpd >= min(modesx,kxp),  unless modesx = nx/2+1,
c in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, modesx, modesy, ndim, kstrt, nyv, kxp
      integer modesxpd, modesyd
      complex vpot, vpott
      dimension vpot(ndim,nyv,kxp)
      dimension vpott(ndim,modesyd,modesxpd)
c local data
      integer nxh, nyh, jmax, kmax, ny2, i, j, k, j1, k1, ks, joff
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 1
      if (kstrt.gt.nxh) return
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*ks
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
      do 40 j = 1, jmax
      if ((j+joff).gt.1) then
         do 20 k = 2, kmax
         k1 = ny2 - k
         do 10 i = 1, ndim
         vpott(i,2*k-2,j) = vpot(i,k,j)
         vpott(i,2*k-1,j) = vpot(i,k1,j)
   10    continue
   20    continue
c mode numbers ky = 0, ny/2
         do 30 i = 1, ndim
         vpott(i,1,j) = vpot(i,1,j)
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            vpott(i,ny,j) = vpot(i,k1,j)
         endif
   30    continue
      endif
   40 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 60 k = 2, kmax
         k1 = ny2 - k
         do 50 i = 1, ndim
         vpott(i,2*k-2,1) = vpot(i,k,1)
         vpott(i,2*k-1,1) = conjg(vpot(i,k,1))
         if (modesx.gt.nxh) then
            vpott(i,2*k-2,j1) = conjg(vpot(i,k1,1))
            vpott(i,2*k-1,j1) = vpot(i,k1,1)
         endif
   50    continue
   60    continue
         do 70 i = 1, ndim
         vpott(i,1,1) = cmplx(real(vpot(i,1,1)),0.0)
         if (modesx.gt.nxh) then
            vpott(i,1,j1) = cmplx(aimag(vpot(i,1,1)),0.0)
         endif
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            vpott(i,ny,1) = cmplx(real(vpot(i,k1,1)),0.0)
            if (modesx.gt.nxh) then
               vpott(i,ny,j1) = cmplx(aimag(vpot(i,k1,1)),0.0)
            endif
         endif
   70    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPWRVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,kstrt, 
     1nyv,kxp,modesxpd,modesyd)
c this subroutine extracts lowest order modes from a location in an
c unpacked complex vector array vpott and stores them into a packed
c complex vector array vpot
c modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
c and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c except kx = NX/2 is stored at location kxp+1 when idproc=0.
c nx/ny = system length in x/y direction
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c ndim = number of field arrays, must be >= 1
c kstrt = starting data block number
c nyv = second dimension of input array vpot, nyv >= ny
c kxp = number of data values per block
c modesyd = third dimension of array vpott,
c where modesyd  >= min(2*modesy-1,ny)
c modesxpd = fourth dimension of array vpott,
c modesxpd >= min(modesx,kxp) unless modesx = nx/2+1,
c in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, modesx, modesy, ndim, kstrt, nyv, kxp
      integer modesxpd, modesyd
      complex vpot, vpott
      dimension vpot(ndim,nyv,kxp)
      dimension vpott(ndim,modesyd,modesxpd)
c local data
      integer nxh, nyh, jmax, kmax, ny2, i, j, k, j1, k1, ks, joff
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 1
      zero = cmplx(0.0,0.0)
      if (kstrt.gt.nxh) return
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*ks
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
      do 60 j = 1, jmax
      if ((j+joff).gt.1) then
         do 20 k = 2, kmax
         do 10 i = 1, ndim
         k1 = ny2 - k
         vpot(i,k,j) = vpott(i,2*k-2,j)
         vpot(i,k1,j) = vpott(i,2*k-1,j)
   10    continue
   20    continue
         do 40 k = kmax+1, nyh
         k1 = ny2 - k
         do 30 i = 1, ndim
         vpot(i,k,j) = zero
         vpot(i,k1,j) = zero
   30    continue
   40    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         do 50 i = 1, ndim
         vpot(i,1,j) = vpott(i,1,j)
         vpot(i,k1,j) = zero
         if (modesy.gt.nyh) then
            vpot(i,k1,j) = vpott(i,ny,j)
         endif
   50    continue
      endif
   60 continue
      do 100 j = jmax+1, kxp
      if ((j+joff).gt.1) then
         do 80 k = 2, nyh
         k1 = ny2 - k
         do 70 i = 1, ndim
         vpot(i,k,j) = zero
         vpot(i,k1,j) = zero
   70    continue
   80    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         do 90 i = 1, ndim
         vpot(i,1,j) = zero
         vpot(i,k1,j) = zero
   90    continue
      endif
  100 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 120 k = 2, kmax
         k1 = ny2 - k
         do 110 i = 1, ndim
         vpot(i,k,1) = vpott(i,2*k-2,1)
         vpot(i,k1,1) = zero
         if (modesx.gt.nxh) then
            vpot(i,k1,1) = conjg(vpott(i,2*k-2,j1))
         endif
  110    continue
  120    continue
         do 140 k = kmax+1, nyh
         k1 = ny2 - k
         do 130 i = 1, ndim
         vpot(i,k,1) = zero
         vpot(i,k1,1) = zero
  130    continue
  140    continue
         k1 = nyh + 1
         do 150 i = 1, ndim
         vpot(i,1,1) = cmplx(real(vpott(i,1,1)),0.0)
         vpot(i,k1,1) = zero
         if (modesx.gt.nxh) then
            vpot(i,1,1) = cmplx(real(vpot(i,1,1)),real(vpott(i,1,j1)))
         endif
         if (modesy.gt.nyh) then
            vpot(i,k1,1) = cmplx(real(vpott(i,ny,1)),0.0)
            if (modesx.gt.nxh) then
               vpot(i,k1,1) = cmplx(real(vpot(i,k1,1)),
     1                              real(vpott(i,ny,j1)))
            endif
         endif
  150    continue
      endif
      return
      end
