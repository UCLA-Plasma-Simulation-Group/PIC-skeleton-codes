c Fortran library for Skeleton 1-2/2D Darwin PIC Code field diagnostics
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
      subroutine CURLF1(f,g,nx,nxvh)
c this subroutine calculates the curl in fourier space
c input: all except g, output: g
c approximate flop count is: 10*nxc
c where nxc = nx/2 - 1
c the curl is calculated using the equations:
c gy(kx) = -sqrt(-1)*kx*fz(kx)
c gz(kx) = sqrt(-1)*kx*fy(kx)
c where kx = 2pi*j/nx, and j = fourier mode number,
c except for gy(kx=pi) = gz(kx=pi) = 0,
c nx = system length in x direction
c nxvh = first dimension of field arrays, must be >= nxh
      implicit none
      integer nx, nxvh
      complex f, g
      dimension f(2,nxvh), g(2,nxvh)
c local data
      integer nxh, j
      real dnx, dkx
      complex zero, zt1, zt2
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
c calculate the curl
c mode numbers 0 < kx < nx/2
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = cmplx(-aimag(f(2,j)),real(f(2,j)))
      zt2 = cmplx(-aimag(f(1,j)),real(f(1,j)))
      g(1,j) = -dkx*zt1
      g(2,j) = dkx*zt2
   40 continue
      g(1,1) = zero
      g(2,1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine APOTP13(cu,ayz,ffc,ci,wm,nx,nxvh,nxhd)
c this subroutine solves 1-2/2d poisson's equation in fourier space for
c vector potential with periodic boundary conditions.
c input: cu,ffc,ci,nx,nxvh,nxhd, output: ayz,wm
c approximate flop count is: 23*nxc
c where nxc = nx/2 - 1
c vector potential is calculated using the equation:
c by(kx) = ci*ci*g(kx)*cuy(kx)
c bz(kx) = ci*ci*g(kx)*cuz(kx)
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/kx**2)*s(kx),
c s(kx) = exp(-((kx*ax)**2+)/2), except for
c by(kx=pi) = bz(kx=pi) = 0, and by(kx=0) = bz(kx=0) = 0.
c cu(i,j) = complex current density for fourier mode (j-1)
c ayz(1,j) = y component of complex vector potential
c ayz(2,j) = z component of complex vector potential
c all for fourier mode (j-1)
c aimag(ffc(j)) = finite-size particle shape factor s
c real(ffc(j)) = potential green's function g
c for fourier mode (j-1)
c ci = reciprocal of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2)
c where affp = normalization constant = nx/np,
c where np=number of particles
c this expression is valid only if the current is divergence-free
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
c nxhd = dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, wm
      complex cu, ayz, ffc
      dimension cu(2,nxvh), ayz(2,nxvh), ffc(nxhd)
c local data
      integer j, nxh
      real ci2, at1, at2
      complex zero
      double precision wp
      nxh = nx/2
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
c calculate vector potential and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at2 = ci2*real(ffc(j))
      at1 = at2*aimag(ffc(j))
      ayz(1,j) = at2*cu(1,j)
      ayz(2,j) = at2*cu(2,j)
      wp = wp + at1*(cu(1,j)*conjg(cu(1,j)) + cu(2,j)*conjg(cu(2,j)))
   10 continue
      ayz(1,1) = zero
      ayz(2,1) = zero
      wm = real(nx)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine ETFIELD13(dcu,eyz,ffe,ci,wf,nx,nxvh,nxhd)
c this subroutine solves 1-2/2d poisson's equation in fourier space for
c unsmoothed transverse electric field, with periodic boundary
c conditions.
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c input: dcu,ffe,isign,ci,nx,nxvh,nxhd, output: eyz,wf
c approximate flop count is: 25*nxc
c where nxc = nx/2 - 1
c unsmoothed transverse electric field is calculated using the equation:
c ey(kx) = -ci*ci*g(kx)*dcuy(kx)
c ez(kx) = -ci*ci*g(kx)*dcuz(kx)
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/(kx**2+wp0*ci2*s(kx)**2))*s(kx),
c s(kx) = exp(-((kx*ax)**2+)/2), except for
c ey(kx=pi) = ez(kx=pi) = 0, and ey(kx=0) = ez(kx=0,) = 0.
c dcu(i,j) = transverse part of complex derivative of current for
c fourier mode (j-1)
c eyz(1,j) = y component of complex transverse electric field
c eyz(2,j) = z component of complex transverse electric field
c all for fourier mode (j-1)
c aimag(ffe(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffe(j)) = potential green's function g
c for fourier mode (j-1)
c ci = reciprocal of velocity of light
c transverse electric field energy is also calculated, using
c wf = nx*sum((affp/(kx**2*ci*ci)**2)*|dcu(kx)*s(kx)|**2)
c where affp = normalization constant = nx/np, and where
c np=number of particles
c this expression is valid only if the derivative of current is
c divergence-free
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
c nxhd = second dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, wf
      complex dcu, eyz, ffe
      dimension dcu(2,nxvh), eyz(2,nxvh)
      dimension ffe(nxhd)
c local data
      integer nxh, j
      real ci2, at1, at2
      complex zero
      double precision wp
      nxh = nx/2
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
c calculate unsmoothed transverse electric field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at2 = -ci2*real(ffe(j))
      at1 = at2*at2
      eyz(1,j) = at2*dcu(1,j)
      eyz(2,j) = at2*dcu(2,j)
      wp = wp + at1*(dcu(1,j)*conjg(dcu(1,j)) + dcu(2,j)*conjg(dcu(2,j))
     1)
   10 continue
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*wp/real(ffe(1))
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
      subroutine SMOOTH13(cu,cus,ffc,nx,nxvh,nxhd)
c this subroutine provides a 1d vector smoothing function
c in fourier space, with periodic boundary conditions.
c input: cu,ffc,nx,nxvh, output: cus
c approximate flop count is: 4*nxc, where nxc = nx/2 - 1
c smoothing is calculated using the equation:
c cusy(kx) = cuy(kx)*s(kx)
c cusz(kx) = cuz(kx)*s(kx)
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/kx**2)*s(kx),
c where affp = normalization constant = nx/np,
c and np=number of particles
c s(kx) = exp(-((kx*ax)**2+)/2), except for
c cusy(kx=pi) = cusz(kx=pi) = 0
c cu(i,j) = ith-component of complex current density
c cus(i,j) = ith-component of complex smoothed charge density
c all for fourier mode (j-1)
c aimag(ffc(j)) = finite-size particle shape factor s
c real(ffc(j)) = potential green's function g
c for fourier mode (j-1)
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
c nxhd = dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      complex cu, cus, ffc
      dimension cu(2,nxvh), cus(2,nxvh), ffc(nxhd)
c local data
      integer j, nxh
      real at1
      nxh = nx/2
c calculate smoothing
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at1 = aimag(ffc(j))
      cus(1,j) = at1*cu(1,j)
      cus(2,j) = at1*cu(2,j)
   10 continue
      at1 = aimag(ffc(1))
      cus(1,1) = cmplx(at1*real(cu(1,1)),0.0)
      cus(2,1) = cmplx(at1*real(cu(2,1)),0.0)
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
c-----------------------------------------------------------------------
      subroutine RDVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
c this subroutine extracts lowest order modes from packed complex vector
c array vpot and stores them into a location in an unpacked complex
c vector array vpott
c modes stored: kx=(0,1,...,NX/2)
c nx = system length in x direction
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c ndim = number of field arrays, must be >= 1
c nxvh = second dimension of input array vpot, nxvh >= nx/2
c modesxd = second dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, modesx, ndim, nxvh, modesxd
      complex vpot, vpott
      dimension vpot(ndim,nxvh), vpott(ndim,modesxd)
c local data
      integer nxh, jmax, i, j, j1
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpott(i,j) = vpot(i,j)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 i = 1, ndim
      vpott(i,1) = cmplx(real(vpot(i,1)),0.0)
      if (modesx.gt.nxh) then
         vpott(i,j1) = cmplx(aimag(vpot(i,1)),0.0)
      endif
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WRVMODES1(vpot,vpott,nx,modesx,ndim,nxvh,modesxd)
c this subroutine extracts lowest order modes from a location in an
c unpacked complex vector array vpott and stores them into a packed
c complex vector array vpot
c modes stored: kx=(0,1,...,NX/2)
c nx = system length in x direction
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c ndim = number of field arrays, must be >= 1
c nxvh = second dimension of input array vpot, nxvh >= nx/2
c modesxd = second dimension of output array vpott, modesxd >= modesx
      implicit none
      integer nx, modesx, ndim, nxvh, modesxd
      complex vpot, vpott
      dimension vpot(ndim,nxvh), vpott(ndim,modesxd)
c local data
      integer nxh, jmax, i, j, j1
      complex zero
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nx/2
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpot(i,j) = vpott(i,j)
   10 continue
   20 continue
      do 40 j = jmax+1, nxh
      do 30 i = 1, ndim
      vpot(i,j) = zero
   30 continue
   40 continue
c mode numbers kx = 0, nx/2
      do 50 i = 1, ndim
      vpot(i,1) = cmplx(real(vpott(i,1)),0.0)
      if (modesx.gt.nxh) then
         vpot(i,1) = cmplx(real(vpot(i,1)),real(vpott(i,j1)))
      endif
   50 continue
      return
      end
