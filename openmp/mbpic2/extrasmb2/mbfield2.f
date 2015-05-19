c Fortran library for Skeleton 2-1/2D Electromagnetic PIC Code field
c diagnostics
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c-----------------------------------------------------------------------
      subroutine MPOTP2(q,pot,ffc,we,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c potential with periodic boundary conditions.
c input: q,ffc,nx,ny,nxvh,nyv,nxhd,nyhd, output: pot,we
c approximate flop count is: 14*nxc*nyc + 8*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c pot(kx,ky) = g(kx,ky)*q(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c pot(kx=pi) = 0, pot(ky=pi) = 0, and pot(kx=0,ky=0) = 0.
c q(j,k)) = complex charge density for fourier mode (j-1,k-1)
c pot(j,k) = complex potential ffor fourier mode (j-1,k-1)
c real(ffc(j,k) = finite-size particle shape factor s
c imag(ffc(j,k) = potential green's function g
c for fourier mode (j-1,k-1)
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c where affp = normalization constant = nx*ny/np,
c where np=number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer  nx, ny, nxvh, nyv, nxhd, nyhd
      real we, sum1
      complex q, pot, ffc
      dimension q(nxvh,nyv), pot(nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real at1, at2
      complex zero
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      zero = cmplx(0.0,0.0)
c calculate potential and sum field energy
      sum1 = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1,at2,wp) REDUCTION(+:sum1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      wp = 0.0d0
      do 10 j = 2, nxh
      at2 = real(ffc(j,k))
      at1 = at2*aimag(ffc(j,k))
      pot(j,k) = at2*q(j,k)
      pot(j,k1) = at2*q(j,k1)
      wp = wp + at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
   10 continue
c mode numbers kx = 0, nx/2
      at2 = real(ffc(1,k))
      at1 = at2*aimag(ffc(1,k))
      pot(1,k) = at2*q(1,k)
      pot(1,k1) = zero
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      wp = 0.0d0
      k1 = nyh + 1
      do 30 j = 2, nxh
      at2 = real(ffc(j,1))
      at1 = at2*aimag(ffc(j,1))
      pot(j,1) = at2*q(j,1)
      pot(j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   30 continue
      pot(1,1) = zero
      pot(1,k1) = zero
      sum1 = sum1 + wp
      we = real(nx*ny)*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine MDIVF2(f,df,nx,ny,ndim,nxvh,nyv)
c this subroutine calculates the divergence in fourier space
c input: all except df, output: df
c approximate flop count is: 16*nxc*nyc + 5*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the divergence is calculated using the equation:
c df(kx,ky) = sqrt(-1)*(kx*fx(kx,ky)+ky*fy(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for df(kx=pi) = df(ky=pi) = df(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, ndim, nxvh, nyv
      complex f, df
      dimension f(ndim,nxvh,nyv), df(nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dkx, dky
      complex zero, zt1
      if (ndim.lt.2) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c calculate the divergence
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,dkx,zt1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = dkx*f(1,j,k) + dky*f(2,j,k)
      df(j,k) = cmplx(-aimag(zt1),real(zt1))
      zt1 = dkx*f(1,j,k1) - dky*f(2,j,k1)
      df(j,k1) = cmplx(-aimag(zt1),real(zt1))
   10 continue
c mode numbers kx = 0, nx/2
      df(1,k) = dky*cmplx(-aimag(f(2,1,k)),real(f(2,1,k)))
      df(1,k1) = zero
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      dkx = dnx*real(j - 1)
      df(j,1) = dkx*cmplx(-aimag(f(1,j,1)),real(f(1,j,1)))
      df(j,k1) = zero
   30 continue
      df(1,1) = zero
      df(1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRADF2(df,f,nx,ny,ndim,nxvh,nyv)
c this subroutine calculates the gradient in fourier space
c input: all except f, output: f
c approximate flop count is: 12*nxc*nyc + 4*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the gradient is calculated using the equations:
c fx(kx,ky) = sqrt(-1)*kx*df(kx,ky)
c fy(kx,ky) = sqrt(-1)*ky*df(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for fx(kx=pi) = fy(kx=pi) = 0, fx(ky=pi) = fy(ky=pi) = 0,
c and fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, ndim, nxvh, nyv
      complex df, f
      dimension df(nxvh,nyv), f(ndim,nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dkx, dky
      complex zero, zt1
      if (ndim.lt.2) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c calculate the gradient
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,dkx,zt1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = cmplx(-aimag(df(j,k)),real(df(j,k)))
      f(1,j,k) = dkx*zt1
      f(2,j,k) = dky*zt1
      zt1 = cmplx(-aimag(df(j,k1)),real(df(j,k1)))
      f(1,j,k1) = dkx*zt1
      f(2,j,k1) = -dky*zt1
   10 continue
c mode numbers kx = 0, nx/2
      f(1,1,k) = zero
      f(2,1,k) = dky*cmplx(-aimag(df(1,k)),real(df(1,k)))
      f(1,1,k1) = zero
      f(2,1,k1) = zero
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      dkx = dnx*real(j - 1)
      f(1,j,1) = dkx*cmplx(-aimag(df(j,1)),real(df(j,1)))
      f(2,j,1) = zero
      f(1,j,k1) = zero
      f(2,j,k1) = zero
   30 continue
      f(1,1,1) = zero
      f(2,1,1) = zero
      f(1,1,k1) = zero
      f(2,1,k1) = zero
      if (ndim.eq.2) return
c handle case of ndim = 3
      do 50 k = 2, nyh
      k1 = ny2 - k
      do 40 j = 2, nxh
      f(3,j,k) = zero
      f(3,j,k1) = zero
   40 continue
   50 continue
      k1 = nyh + 1
      do 60 j = 2, nxh
      f(3,j,1) = zero
      f(3,j,k1) = zero
   60 continue
      f(3,1,1) = zero
      f(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MCURLF2(f,g,nx,ny,nxvh,nyv)
c this subroutine calculates the curl in fourier space
c input: all except g, output: g
c approximate flop count is: 32*nxc*nyc + 10*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the curl is calculated using the equations:
c gx(kx,ky) = sqrt(-1)*ky*fz(kx,ky)
c gy(kx,ky) = -sqrt(-1)*kx*fz(kx,ky)
c gz(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
c and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex f, g
      dimension f(3,nxvh,nyv), g(3,nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dkx, dky
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c calculate the curl
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,dkx,zt1,zt2,zt3)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = cmplx(-aimag(f(3,j,k)),real(f(3,j,k)))
      zt2 = cmplx(-aimag(f(2,j,k)),real(f(2,j,k)))
      zt3 = cmplx(-aimag(f(1,j,k)),real(f(1,j,k)))
      g(1,j,k) = dky*zt1
      g(2,j,k) = -dkx*zt1
      g(3,j,k) = dkx*zt2 - dky*zt3
      zt1 = cmplx(-aimag(f(3,j,k1)),real(f(3,j,k1)))
      zt2 = cmplx(-aimag(f(2,j,k1)),real(f(2,j,k1)))
      zt3 = cmplx(-aimag(f(1,j,k1)),real(f(1,j,k1)))
      g(1,j,k1) = -dky*zt1
      g(2,j,k1) = -dkx*zt1
      g(3,j,k1) = dkx*zt2 + dky*zt3
   10 continue
c mode numbers kx = 0, nx/2
      zt1 = cmplx(-aimag(f(3,1,k)),real(f(3,1,k)))
      zt3 = cmplx(-aimag(f(1,1,k)),real(f(1,1,k)))
      g(1,1,k) = dky*zt1
      g(2,1,k) = zero
      g(3,1,k) = -dky*zt3
      g(1,1,k1) = zero
      g(2,1,k1) = zero
      g(3,1,k1) = zero
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = cmplx(-aimag(f(3,j,1)),real(f(3,j,1)))
      zt2 = cmplx(-aimag(f(2,j,1)),real(f(2,j,1)))
      g(1,j,1) = zero
      g(2,j,1) = -dkx*zt1
      g(3,j,1) = dkx*zt2
      g(1,j,k1) = zero
      g(2,j,k1) = zero
      g(3,j,k1) = zero
   30 continue
      g(1,1,1) = zero
      g(2,1,1) = zero
      g(3,1,1) = zero
      g(1,1,k1) = zero
      g(2,1,k1) = zero
      g(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MAVPOT23(bxy,axy,nx,ny,nxvh,nyv)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with periodic boundary conditions.
c input: bxy, nx, ny, nxvh, nyv, output: axy
c approximate flop count is: 38*nxc*nyc + 10*(nxc + nyc)
c and nxc*nyc divides, where nxc = nx/2 - 1, nyc = ny/2 - 1
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = ax(ky=pi) = ay(ky=pi) = az(ky=pi) 
c = 0, and ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c bxy(i,j,k) = i component of complex magnetic field
c axy(i,j,k) = i component of complex vector potential
c all for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex bxy, axy
      dimension bxy(3,nxvh,nyv), axy(3,nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dkx, dky, dky2
      real at1, at2, at3
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c calculate vector potential
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,dky2,dkx,at1,at2,at3,zt1,zt2,zt3)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx + dky2)
      at2 = dkx*at1
      at3 = dky*at1
      zt1 = cmplx(-aimag(bxy(3,j,k)),real(bxy(3,j,k)))
      zt2 = cmplx(-aimag(bxy(2,j,k)),real(bxy(2,j,k)))
      zt3 = cmplx(-aimag(bxy(1,j,k)),real(bxy(1,j,k)))
      axy(1,j,k) = at3*zt1
      axy(2,j,k) = -at2*zt1
      axy(3,j,k) = at2*zt2 - at3*zt3
      zt1 = cmplx(-aimag(bxy(3,j,k1)),real(bxy(3,j,k1)))
      zt2 = cmplx(-aimag(bxy(2,j,k1)),real(bxy(2,j,k1)))
      zt3 = cmplx(-aimag(bxy(1,j,k1)),real(bxy(1,j,k1)))
      axy(1,j,k1) = -at3*zt1
      axy(2,j,k1) = -at2*zt1
      axy(3,j,k1) = at2*zt2 + at3*zt3
   10 continue
c mode numbers kx = 0, nx/2
      at3 = 1.0/dky
      zt1 = cmplx(-aimag(bxy(3,1,k)),real(bxy(3,1,k)))
      zt3 = cmplx(-aimag(bxy(1,1,k)),real(bxy(1,1,k)))
      axy(1,1,k) = at3*zt1
      axy(2,1,k) = zero
      axy(3,1,k) = -at3*zt3
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      dkx = dnx*real(j - 1)
      at2 = 1.0/dkx
      zt1 = cmplx(-aimag(bxy(3,j,1)),real(bxy(3,j,1)))
      zt2 = cmplx(-aimag(bxy(2,j,1)),real(bxy(2,j,1)))
      axy(1,j,1) = zero
      axy(2,j,1) = -at2*zt1
      axy(3,j,1) = at2*zt2
      axy(1,j,k1) = zero
      axy(2,j,k1) = zero
      axy(3,j,k1) = zero
   30 continue
      axy(1,1,1) = zero
      axy(2,1,1) = zero
      axy(3,1,1) = zero
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MAVRPOT23(axy,bxy,ffc,ci,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c the radiative part of the vector potential
c with periodic boundary conditions.
c input: all, output: axy
c approximate flop count is: 68*nxc*nyc + 20*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
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
c ax(ky=pi) = ay(ky=pi) = ax(ky=pi) = 0,
c ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c axy(i,j,k) = on entry, complex current density cu
c axy(i,j,k) = on exit, complex radiative vector potential
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffc(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real ci
      complex axy, bxy, ffc
      dimension axy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, afc2, dkx, dky, dky2, at1, at2
      complex zero, zt1, zt2, zt3
      if (ci.le.0.) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      afc2 = real(ffc(1,1))*ci*ci
      zero = cmplx(0.,0.)
c calculate the radiative vector potential
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,dky2,dkx,at1,at2,zt1,zt2,zt3)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx + dky2)
      at2 = afc2*aimag(ffc(j,k))
c update radiative vector potential, ky > 0
      zt1 = cmplx(-aimag(bxy(3,j,k)),real(bxy(3,j,k)))
      zt2 = cmplx(-aimag(bxy(2,j,k)),real(bxy(2,j,k)))
      zt3 = cmplx(-aimag(bxy(1,j,k)),real(bxy(1,j,k)))
      axy(1,j,k) = at1*(dky*zt1 - at2*axy(1,j,k))
      axy(2,j,k) = -at1*(dkx*zt1 + at2*axy(2,j,k))
      axy(3,j,k) = at1*((dkx*zt2 - dky*zt3) - at2*axy(3,j,k))
c update radiative vector potential, ky < 0
      zt1 = cmplx(-aimag(bxy(3,j,k1)),real(bxy(3,j,k1)))
      zt2 = cmplx(-aimag(bxy(2,j,k1)),real(bxy(2,j,k1)))
      zt3 = cmplx(-aimag(bxy(1,j,k1)),real(bxy(1,j,k1)))
      axy(1,j,k1) = -at1*(dky*zt1 + at2*axy(1,j,k1))
      axy(2,j,k1) = -at1*(dkx*zt1 + at2*axy(2,j,k1))
      axy(3,j,k1) = at1*((dkx*zt2 + dky*zt3) - at2*axy(3,j,k1))
   10 continue
c mode numbers kx = 0, nx/2
      at1 = 1.0/(dky*dky)
      at2 = afc2*aimag(ffc(1,k))
c update radiative vector potential
      zt1 = cmplx(-aimag(bxy(3,1,k)),real(bxy(3,1,k)))
      zt3 = cmplx(-aimag(bxy(1,1,k)),real(bxy(1,1,k)))
      axy(1,1,k) = at1*(dky*zt1 - at2*axy(1,1,k))
      axy(2,1,k) = zero
      axy(3,1,k) = -at1*(dky*zt3 + at2*axy(3,1,k))
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx)
      at2 = afc2*aimag(ffc(j,1))
c update radiative vector potential
      zt1 = cmplx(-aimag(bxy(3,j,1)),real(bxy(3,j,1)))
      zt2 = cmplx(-aimag(bxy(2,j,1)),real(bxy(2,j,1)))
      axy(1,j,1) = zero
      axy(2,j,1) = -at1*(dkx*zt1 + at2*axy(2,j,1))
      axy(3,j,1) = at1*(dkx*zt2 - at2*axy(3,j,1))
      axy(1,j,k1) = zero
      axy(2,j,k1) = zero
      axy(3,j,k1) = zero
   30 continue
      axy(1,1,1) = zero
      axy(2,1,1) = zero
      axy(3,1,1) = zero
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MSMOOTH2(q,qs,ffc,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine provides a 2d scalar smoothing function
c in fourier space, with periodic boundary conditions.
c input: q,ffc,nx,ny,nxvh,nyv,nxhd,nyhd, output: qs
c approximate flop count is: 4*nxc*nyc + 2*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c smoothing is calculated using the equation:
c qs(kx,ky) = q(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),h
c where affp = normalization constant = nx*ny/np,
c and np=number of particles
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c qs(kx=pi) = qs(kx=pi) = 0, and qs(kx=0,ky=0) = 0.
c q(j,k) = complex charge density
c qs(j,k) = complex smoothed charge density
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = first dimension of scalar field arrays, must be >= nx/2
c nyv = second dimension of scalar field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nx/2
c nyhd = second dimension of form factor array, must be >= ny/2
      implicit none
      integer  nx, ny, nxvh, nyv, nxhd, nyhd
      complex q, qs, ffc
      dimension q(nxvh,nyv), qs(nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real at1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      zero = cmplx(0.0,0.0)
c calculate smoothing
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nxh
      at1 = aimag(ffc(j,k))
      qs(j,k) = at1*q(j,k)
      qs(j,k1) = at1*q(j,k1)
   10 continue
c mode numbers kx = 0, nx/2
      at1 = aimag(ffc(1,k))
      qs(1,k) = at1*q(1,k)
      qs(1,k1) = zero
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      at1 = aimag(ffc(j,1))
      qs(j,1) = at1*q(j,1)
      qs(j,k1) = zero
   30 continue
      qs(1,1) = cmplx(aimag(ffc(1,1))*real(q(1,1)),0.0)
      qs(1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MSMOOTH23(cu,cus,ffc,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine provides a 2d vector smoothing function
c in fourier space, with periodic boundary conditions.
c input: cu,ffc,nx,ny,nxvh,nyv,nxhd,nyhd, output: cus
c approximate flop count is: 12*nxc*nyc + 6*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c smoothing is calculated using the equation:
c cusx(kx,ky) = cux(kx,ky)*s(kx,ky)
c cusy(kx,ky) = cuy(kx,ky)*s(kx,ky)
c cusy(kx,ky) = cuz(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c where affp = normalization constant = nx*ny/np,
c and np=number of particles
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c cusx(kx=pi) = cusy(kx=pi) = cusz(kx=pi) = 0 
c cusx(ky=pi) = cusy(ky=pi) = cusz(ky=pi) = 0
c and cusx(kx=0,ky=0) = cusy(kx=0,ky=0) = cusz(kx=0,ky=0) = 0.
c cu(i,j,k) = ith-component of complex current density
c cus(j,k) = ith-component of complex smoothed charge density
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      complex cu, cus, ffc
      dimension cu(3,nxvh,nyv), cus(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real at1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      zero = cmplx(0.0,0.0)
c calculate smoothing
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nxh
      at1 = aimag(ffc(j,k))
      cus(1,j,k) = at1*cu(1,j,k)
      cus(2,j,k) = at1*cu(2,j,k)
      cus(3,j,k) = at1*cu(3,j,k)
      cus(1,j,k1) = at1*cu(1,j,k1)
      cus(2,j,k1) = at1*cu(2,j,k1)
      cus(3,j,k1) = at1*cu(3,j,k1)
   10 continue
c mode numbers kx = 0, nx/2
      at1 = aimag(ffc(1,k))
      cus(1,1,k) = at1*cu(1,1,k)
      cus(2,1,k) = at1*cu(2,1,k)
      cus(3,1,k) = at1*cu(3,1,k)
      cus(1,1,k1) = zero
      cus(2,1,k1) = zero
      cus(3,1,k1) = zero
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      at1 = aimag(ffc(j,1))
      cus(1,j,1) = at1*cu(1,j,1)
      cus(2,j,1) = at1*cu(2,j,1)
      cus(3,j,1) = at1*cu(3,j,1)
      cus(1,j,k1) = zero
      cus(2,j,k1) = zero
      cus(3,j,k1) = zero
   30 continue
      at1 = aimag(ffc(1,1))
      cus(1,1,1) = cmplx(at1*real(cu(1,1,1)),0.)
      cus(2,1,1) = cmplx(at1*real(cu(2,1,1)),0.)
      cus(3,1,1) = cmplx(at1*real(cu(3,1,1)),0.)
      cus(1,1,k1) = zero
      cus(2,1,k1) = zero
      cus(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine RDMODES2(pot,pott,nx,ny,modesx,modesy,nxvh,nyv,modesxd,
     1modesyd)
c this subroutine extracts lowest order modes from packed complex array
c pot and stores them into a location in an unpacked complex array pott
c modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c nx/ny = system length in x/y direction
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c nxvh = first dimension of input array pot, nxvh >= nx/2
c nyv = second dimension of input array pot, nyv >= ny
c modesxd = first dimension of output array pott, modesxd >= modesx
c modesyd = second dimension of output array pott, modesyd  = 2*modesy
c where modesyd  >= min(2*modesy-1,ny)
      implicit none
      integer nx, ny, modesx, modesy, nxvh, nyv, modesxd, modesyd
      complex pot, pott
      dimension pot(nxvh,nyv), pott(modesxd,modesyd)
c local data
      integer nxh, nyh, kmax, jmax, ny2, j, k, j1, k1
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1)
      do 20 k = 2, kmax
      k1 = ny2 - k
      do 10 j = 2, jmax
      pott(j,2*k-2) = pot(j,k)
      pott(j,2*k-1) = pot(j,k1)
   10 continue
c mode numbers kx = 0, nx/2
      pott(1,2*k-2) = pot(1,k)
      pott(1,2*k-1) = conjg(pot(1,k))
      if (modesx.gt.nxh) then
         pott(j1,2*k-2) = conjg(pot(1,k1))
         pott(j1,2*k-1) = pot(1,k1)
      endif
   20 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      do 30 j = 2, jmax
      pott(j,1) = pot(j,1)
   30 continue
      pott(1,1) = cmplx(real(pot(1,1)),0.0)
      if (modesx.gt.nxh) then
         pott(j1,1) = cmplx(aimag(pot(1,1)),0.0)
      endif
      if (modesy.gt.nyh) then
         k1 = nyh + 1
         do 40 j = 2, jmax
         pott(j,ny) = pot(j,k1)
   40    continue
         pott(1,ny) = cmplx(real(pot(1,k1)),0.0)
         if (modesx.gt.nxh) then
            pott(j1,ny) = cmplx(aimag(pot(1,k1)),0.0)
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WRMODES2(pot,pott,nx,ny,modesx,modesy,nxvh,nyv,modesxd,
     1modesyd)
c this subroutine extracts lowest order modes from a location in an
c unpacked complex array pott and stores them into a packed complex
c array pot
c modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2),
c nx/ny = system length in x/y direction
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c nxvh = first dimension of input array pot, nxvh >= nx/2
c nyv = second dimension of input array pot, nyv >= ny
c modesxd = first dimension of output array pott, modesxd >= modesx
c modesyd = second dimension of output array pott,
c where modesyd  >= min(2*modesy-1,ny)
      implicit none
      integer nx, ny, modesx, modesy, nxvh, nyv, modesxd, modesyd
      complex pot, pott
      dimension pot(nxvh,nyv), pott(modesxd,modesyd)
c local data
      integer nxh, nyh, kmax, jmax, ny2, j, k, j1, k1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      j1 = nxh + 1
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1)
      do 30 k = 2, kmax
      k1 = ny2 - k
      do 10 j = 2, jmax
      pot(j,k) = pott(j,2*k-2)
      pot(j,k1) = pott(j,2*k-1)
   10 continue
      do 20 j = jmax+1, nxh
      pot(j,k) = zero
      pot(j,k1) = zero
   20 continue
c mode numbers kx = 0, nx/2
      pot(1,k) = pott(1,2*k-2)
      pot(1,k1) = zero
      if (modesx.gt.nxh) then
         pot(1,k1) = conjg(pott(j1,2*k-2))
      endif
   30 continue
!$OMP END PARALLEL DO
      do 50 k = kmax+1, nyh
      k1 = ny2 - k
      do 40 j = 1, nxh
      pot(j,k) = zero
      pot(j,k1) = zero
   40 continue
   50 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, jmax
      pot(j,1) = pott(j,1)
      pot(j,k1) = zero
   60 continue
      do 70 j = jmax+1, nxh
      pot(j,1) = zero
      pot(j,k1) = zero
   70 continue
      pot(1,1) = cmplx(real(pott(1,1)),0.0)
      pot(1,k1) = zero
      if (modesx.gt.nxh) then
         pot(1,1) = cmplx(real(pot(1,1)),real(pott(j1,1)))
      endif
      if (modesy.gt.nyh) then
         k1 = nyh + 1
         do 80 j = 2, jmax
         pot(j,k1) = pott(j,ny)
   80    continue
         pot(1,k1) = cmplx(real(pott(1,ny)),0.0)
         if (modesx.gt.nxh) then
            pot(1,k1) = cmplx(real(pot(1,k1)),real(pott(j1,ny)))
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine RDVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,nxvh,nyv,
     1modesxd,modesyd)
c this subroutine extracts lowest order modes from packed complex vector
c array vpot and stores them into a location in an unpacked complex
c vector array vpott
c modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c nx/ny = system length in x/y direction
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c ndim = number of field arrays, must be >= 1
c nxvh = second dimension of input array vpot, nxvh >= nx/2
c nyv = third dimension of input array vpot, nyv >= ny
c modesxd = second dimension of output array vpott, modesxd >= modesx
c modesyd = third dimension of output array vpott,
c where modesyd  >= min(2*modesy-1,ny)
      implicit none
      integer nx, ny, modesx, modesy, ndim, nxvh, nyv, modesxd, modesyd
      complex vpot, vpott
      dimension vpot(ndim,nxvh,nyv), vpott(ndim,modesxd,modesyd)
c local data
      integer nxh, nyh, kmax, jmax, ny2, i, j, k, j1, k1
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(i,j,k,k1)
      do 40 k = 2, kmax
      k1 = ny2 - k
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpott(i,j,2*k-2) = vpot(i,j,k)
      vpott(i,j,2*k-1) = vpot(i,j,k1)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 i = 1, ndim
      vpott(i,1,2*k-2) = vpot(i,1,k)
      vpott(i,1,2*k-1) = conjg(vpot(i,1,k))
      if (modesx.gt.nxh) then
         vpott(i,j1,2*k-2) = conjg(vpot(i,1,k1))
         vpott(i,j1,2*k-1) = vpot(i,1,k1)
      endif
   30 continue
   40 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      do 60 j = 2, jmax
      do 50 i = 1, ndim
      vpott(i,j,1) = vpot(i,j,1)
   50 continue
   60 continue
      do 70 i = 1, ndim
      vpott(i,1,1) = cmplx(real(vpot(i,1,1)),0.0)
      if (modesx.gt.nxh) then
         vpott(i,j1,1) = cmplx(aimag(vpot(i,1,1)),0.0)
      endif
   70 continue
      if (modesy.gt.nyh) then
         k1 = nyh + 1
         do 90 j = 2, jmax
         do 80 i = 1, ndim
         vpott(i,j,ny) = vpot(i,j,k1)
   80    continue
   90    continue
         do 100 i = 1, ndim
         vpott(i,1,ny) = cmplx(real(vpot(i,1,k1)),0.0)
         if (modesx.gt.nxh) then
            vpott(i,j1,ny) = cmplx(aimag(vpot(i,1,k1)),0.0)
         endif
  100    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WRVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,nxvh,nyv,
     1modesxd,modesyd)
c this subroutine extracts lowest order modes from a location in an
c unpacked complex vector array vpott and stores them into a packed
c complex vector array vpot
c modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2),
c nx/ny = system length in x/y direction
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c ndim = number of field arrays, must be >= 1
c nxvh = second dimension of input array vpot, nxvh >= nx/2
c nyv = third dimension of input array vpot, nyv >= ny
c modesxd = second dimension of output array vpott, modesxd >= modesx
c modesyd = third dimension of output array vpott,
c where modesyd  >= min(2*modesy-1,ny)
      implicit none
      integer nx, ny, modesx, modesy, ndim, nxvh, nyv, modesxd, modesyd
      complex vpot, vpott
      dimension vpot(ndim,nxvh,nyv), vpott(ndim,modesxd,modesyd)
c local data
      integer nxh, nyh, kmax, jmax, ny2, i, j, k, j1, k1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      j1 = nxh + 1
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(i,j,k,k1)
      do 60 k = 2, kmax
      k1 = ny2 - k
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpot(i,j,k) = vpott(i,j,2*k-2)
      vpot(i,j,k1) = vpott(i,j,2*k-1)
   10 continue
   20 continue
      do 40 j = jmax+1, nxh
      do 30 i = 1, ndim
      vpot(i,j,k) = zero
      vpot(i,j,k1) = zero
   30 continue
   40 continue
c mode numbers kx = 0, nx/2
      do 50 i = 1, ndim
      vpot(i,1,k) = vpott(i,1,2*k-2)
      vpot(i,1,k1) = zero
      if (modesx.gt.nxh) then
         vpot(i,1,k1) = conjg(vpott(i,j1,2*k-2))
      endif
   50 continue
   60 continue
!$OMP END PARALLEL DO
      do 90 k = kmax+1, nyh
      k1 = ny2 - k
      do 80 j = 1, nxh
      do 70 i = 1, ndim
      vpot(i,j,k) = zero
      vpot(i,j,k1) = zero
   70 continue
   80 continue
   90 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 110 j = 2, jmax
      do 100 i = 1, ndim
      vpot(i,j,1) = vpott(i,j,1)
      vpot(i,j,k1) = zero
  100 continue
  110 continue
      do 130 j = jmax+1, nxh
      do 120 i = 1, ndim
      vpot(i,j,1) = zero
      vpot(i,j,k1) = zero
  120 continue
  130 continue
      do 140 i = 1, ndim
      vpot(i,1,1) = cmplx(real(vpott(i,1,1)),0.0)
      vpot(i,1,k1) = zero
      if (modesx.gt.nxh) then
         vpot(i,1,1) = cmplx(real(vpot(i,1,1)),real(vpott(i,j1,1)))
      endif
  140 continue
      if (modesy.gt.nyh) then
         k1 = nyh + 1
         do 160 j = 2, jmax
         do 150 i = 1, ndim
         vpot(i,j,k1) = vpott(i,j,ny)
  150    continue
  160    continue
         do 170 i = 1, ndim
         vpot(i,1,k1) = cmplx(real(vpott(i,1,ny)),0.0)
         if (modesx.gt.nxh) then
            vpot(i,1,k1) = cmplx(real(vpot(i,1,k1)),
     1                           real(vpott(i,j1,ny)))
         endif
  170    continue
      endif
      return
      end
