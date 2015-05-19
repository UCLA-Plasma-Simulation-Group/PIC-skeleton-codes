c Fortran library for Skeleton 1D Electrostatic PIC Code field
c diagnostics
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c-----------------------------------------------------------------------
      subroutine POTP1(q,pot,ffc,we,nx,nxvh,nxhd)
c this subroutine solves 1d poisson's equation in fourier space for
c potential with periodic boundary conditions.
c input: q,ffc,nx,nxvh,nxhd, output: pot,we
c approximate flop count is: 3*nx
c potential is calculated using the equation:
c pot(k) = g(k)*q(k), where k = 2pi*j/nx, j=fourier mode,
c g(k) = (affp/k**2)*s(k), and s(k) = exp(-(k*ax)**2/2), except for
c pot(k=0) = pot(k=pi) = 0.
c q(j) = complex charge density for fourier mode j-1
c pot(j) = complex potential for fourier mode j-1
c real(ffc(j)) = potential green's function g for fourier mode j-1
c imag(ffc(j)) = finite-size particle shape factor s for fourier mode j-1
c electric field energy is also calculated, using
c we = nx*sum((affp/k**2)*|q(k)*s(k)|**2)
c where affp = normalization constant = nx/np,
c where np = number of particles
c nx = system length in x direction
c nxvh = first dimension of field arrays, must be >= nxh
c nxhd = first dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real we
      complex q, pot, ffc
      dimension q(nxvh), pot(nxvh), ffc(nxhd)
c local data
      integer j, nxh
      real at1, at2
      double precision wp
      nxh = nx/2
c calculate potential and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at2 = real(ffc(j))
      at1 = at2*aimag(ffc(j))
      pot(j) = at2*q(j)
      wp = wp + at1*(q(j)*conjg(q(j)))
   10 continue
      pot(1) = cmplx(0.0,0.0)
      we = real(nx)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine DIVF1(f,df,nx,ndim,nxvh)
c this subroutine calculates the divergence in fourier space
c input: all except df, output: df
c approximate flop count is: 15*nxc
c where nxc = nx/2 - 1
c the divergence is calculated using the equation:
c df(kx) = sqrt(-1)*kx*fx(kx)
c where kx = 2pi*j/nx, and j = fourier mode number,
c except for df(kx=pi) = 0.
c nx = system length in x direction
c ndim = number of field arrays
c nxvh = first dimension of field arrays, must be >= nxh
      implicit none
      integer nx, ndim, nxvh
      complex f, df
      dimension f(ndim,nxvh), df(nxvh)
c local data
      integer nxh, j
      real dnx, dkx
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
c calculate the divergence
c mode numbers 0 < kx < nx/2
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      df(j) = dkx*cmplx(-aimag(f(1,j)),real(f(1,j)))
   40 continue
      df(1) = cmplx(0.0,0.0)
      return
      end
c-----------------------------------------------------------------------
      subroutine GRADF1(df,f,nx,ndim,nxvh)
c this subroutine calculates the gradient in fourier space
c input: all except f, output: f
c approximate flop count is: 5*nxc
c where nxc = nx/2 - 1
c the gradient is calculated using the equations:
c fx(kx) = sqrt(-1)*kx*df(kx)
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c except for fx(kx=pi) = 0,
c nx = system length in x direction
c ndim = number of field arrays
c nxvh = first dimension of field arrays, must be >= nxh
      implicit none
      integer nx, ndim, nxvh
      complex df, f
      dimension df(nxvh), f(ndim,nxvh)
c local data
      integer nxh, i, j
      real dnx, dkx
      complex zero
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
c calculate the gradient
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      f(1,j) = dkx*cmplx(-aimag(df(j)),real(df(j)))
   10 continue
      f(1,1) = zero
      if (ndim.eq.1) return
c handle case of ndim > 1
      do 30 j = 1, nxh
      do 20 i = 2, ndim
      f(i,j) = zero
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine SMOOTH1(q,qs,ffc,nx,nxvh,nxhd)
c this subroutine provides a 1d scalar smoothing function
c in fourier space, with periodic boundary conditions.
c input: q,ffc,nx,nxvh,nxhd, output: qs, flop count = nx
c smoothing is calculated using the equation:
c qs(k) = q(k)*s(k), where k = 2pi*j/nx, j=fourier mode,
c g(k) = (affp/k**2)*s(k), and s(k) = exp(-(k*ax)**2/2),
c except for qs(k=pi) = 0.
c q(j) = complex charge density for fourier mode j-1
c qs(j) = complex smoothed charge density for fourier mode j-1
c real(ffc(j,k)) = potential green's function g
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode j-1
c nx = system length in x direction
c nxvh = first dimension of scalar field arrays, must be >= nx/2
c nxhd = first dimension of form factor array, must be >= nx/2
      implicit none
      integer nx, nxvh, nxhd
      complex q, qs, ffc
      dimension q(nxvh), qs(nxvh), ffc(nxhd)
c local data
      integer j, nxh
      real at1
      nxh = nx/2
c calculate smoothing
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at1 = aimag(ffc(j))
      qs(j) = at1*q(j)
   10 continue
      qs(1) = cmplx(aimag(ffc(1))*real(q(1)),0.0)
      return
      end
c-----------------------------------------------------------------------
      subroutine RDMODES1(pot,pott,nx,modesx,nxvh,modesxd)
c this subroutine extracts lowest order modes from packed complex array
c pot and stores them into a location in an unpacked complex array pott
c vector array pott
c modes stored: kx=(0,1,...,NX/2)
c nx = system length in x direction
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c nxvh = second dimension of input array vpot, nxvh >= nx/2
c modesxd = second dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, modesx, nxvh, modesxd
      complex pot, pott
      dimension pot(nxvh), pott(modesxd)
c local data
      integer nxh, jmax, j, j1
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2
      do 10 j = 2, jmax
      pott(j) = pot(j)
   10 continue
c mode numbers kx = 0, nx/2
      pott(1) = cmplx(real(pot(1)),0.0)
      if (modesx.gt.nxh) then
         pott(j1) = cmplx(aimag(pot(1)),0.0)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WRMODES1(pot,pott,nx,modesx,nxvh,modesxd)
c this subroutine extracts lowest order modes from a location in an
c unpacked complex array pott and stores them into a packed complex
c array pot
c modes stored: kx=(0,1,...,NX/2)
c nx = system length in x direction
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c nxvh = second dimension of input array vpot, nxvh >= nx/2
c modesxd = second dimension of output array vpott, modesxd >= modesx
      implicit none
      integer nx, modesx, nxvh, modesxd
      complex pot, pott
      dimension pot(nxvh), pott(modesxd)
c local data
      integer nxh, jmax, j, j1
      complex zero
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nx/2
      do 10 j = 2, jmax
      pot(j) = pott(j)
   10 continue
      do 20 j = jmax+1, nxh
      pot(j) = zero
   20 continue
c mode numbers kx = 0, nx/2
      pot(1) = cmplx(real(pott(1)),0.0)
      if (modesx.gt.nxh) then
         pot(1) = cmplx(real(pot(1)),real(pott(j1)))
      endif
      return
      end
