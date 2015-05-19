c Fortran library for Skeleton 2D Electrostatic PIC Code field
c diagnostics
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c-----------------------------------------------------------------------
      subroutine POTP2(q,pot,ffc,we,nx,ny,nxvh,nyv,nxhd,nyhd)
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
c if isign = 0, form factor array is prepared
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
      real we
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
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nxh
      at2 = real(ffc(j,k))
      at1 = at2*aimag(ffc(j,k))
      pot(j,k) = at2*q(j,k)
      pot(j,k1) = at2*q(j,k1)
      wp = wp + at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      at2 = real(ffc(1,k))
      at1 = at2*aimag(ffc(1,k))
      pot(1,k) = at2*q(1,k)
      pot(1,k1) = zero
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      at2 = real(ffc(j,1))
      at1 = at2*aimag(ffc(j,1))
      pot(j,1) = at2*q(j,1)
      pot(j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   40 continue
      pot(1,1) = zero
      pot(1,k1) = zero
      we = real(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine DIVF2(f,df,nx,ny,ndim,nxvh,nyv)
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
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      df(1,k) = dky*cmplx(-aimag(f(2,1,k)),real(f(2,1,k)))
      df(1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      df(j,1) = dkx*cmplx(-aimag(f(1,j,1)),real(f(1,j,1)))
      df(j,k1) = zero
   40 continue
      df(1,1) = zero
      df(1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine GRADF2(df,f,nx,ny,ndim,nxvh,nyv)
c this subroutine calculates the gradient in fourier space
c input: all except f, output: f
c approximate flop count is: 12*nxc*nyc + 4*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the gradient is calculated using the equations:
c fx(kx,ky) = sqrt(-1)*kx*df(kx,ky)
c fy(kx,ky) = sqrt(-1)*ky*df(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k,l = fourier mode numbers,
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
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      f(1,1,k) = zero
      f(2,1,k) = dky*cmplx(-aimag(df(1,k)),real(df(1,k)))
      f(1,1,k1) = zero
      f(2,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      f(1,j,1) = dkx*cmplx(-aimag(df(j,1)),real(df(j,1)))
      f(2,j,1) = zero
      f(1,j,k1) = zero
      f(2,j,k1) = zero
   40 continue
      f(1,1,1) = zero
      f(2,1,1) = zero
      f(1,1,k1) = zero
      f(2,1,k1) = zero
      if (ndim.eq.2) return
c handle case of ndim = 3
      do 60 k = 2, nyh
      k1 = ny2 - k
      do 50 j = 2, nxh
      f(3,j,k) = zero
      f(3,j,k1) = zero
   50 continue
   60 continue
      k1 = nyh + 1
      do 70 j = 2, nxh
      f(3,j,1) = zero
      f(3,j,k1) = zero
   70 continue
      f(3,1,1) = zero
      f(3,1,k1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine SMOOTH2(q,qs,ffc,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine provides a 2d scalar smoothing function
c in fourier space, with periodic boundary conditions.
c input: q,ffc,nx,ny,nxv,nyhd, output: qs
c approximate flop count is: 4*nxc*nyc + 2*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c smoothing is calculated using the equation:
c qs(kx,ky) = q(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c qs(kx=pi) = qs(kx=pi) = 0, and qs(kx=0,ky=0) = 0.
c q(j,k) = complex charge density
c qs(j,k) = complex smoothed charge density
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
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
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nxh
      at1 = aimag(ffc(j,k))
      qs(j,k) = at1*q(j,k)
      qs(j,k1) = at1*q(j,k1)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      at1 = aimag(ffc(1,k))
      qs(1,k) = at1*q(1,k)
      qs(1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      at1 = aimag(ffc(j,1))
      qs(j,1) = at1*q(j,1)
      qs(j,k1) = zero
   40 continue
      qs(1,1) = cmplx(aimag(ffc(1,1))*real(q(1,1)),0.0)
      qs(1,k1) = zero
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
