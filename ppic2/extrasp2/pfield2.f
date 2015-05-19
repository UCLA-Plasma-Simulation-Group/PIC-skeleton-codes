c Fortran library for Skeleton 2D Electrostatic MPI PIC Code field
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
