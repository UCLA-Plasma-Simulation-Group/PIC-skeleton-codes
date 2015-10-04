c Fortran library for Skeleton 3D Electrostatic MPI PIC Code field
c diagnostics
c written by viktor k. decyk, ucla
c copyright 1994-2015, regents of the university of california
c-----------------------------------------------------------------------
      subroutine PPOTP32(q,pot,ffc,we,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,
     1kyzp,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c potential, with periodic boundary conditions for distributed data,
c with 2D spatial decomposition
c input: q,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,nzhd
c output: pot,we
c approximate flop count is:
c 41*nxc*nyc*nzc + 21*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c potential is calculated using the equation:
c pot(kx,ky,kz) = g(kx,ky,kz)*q(kx,ky,kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c pot(kx=pi) = 0, pot(ky=pi) = 0, pot(kz=pi) = 0, and
c pot(kx=0,ky=0,kz=0) = 0.
c q(l,j,k) = complex charge density for fourier mode jj-1,kk-1,l-1
c pot(l,j,k) = complex potential
c aimag(ffc(l,j,k)) = finite-size particle shape factor s
c real(ffc(l,j,k)) = potential green's function g
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
c kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c where affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      real we
      complex q, pot, ffc
      dimension q(nzv,kxyp,kyzp), pot(nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
c local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1
      real at1, at2
      complex zero
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      zero = cmplx(0.0,0.0)
c find processor id and offsets in y/z
c js/ks = processor co-ordinates in x/y => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
c calculate potential and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 120
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         do 20 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at2 = real(ffc(l,j,k))
            at1 = at2*aimag(ffc(l,j,k))
            pot(l,j,k) = at2*q(l,j,k)
            pot(l1,j,k) = at2*q(l1,j,k)
            wp = wp + at1*(q(l,j,k)*conjg(q(l,j,k))
     1              + q(l1,j,k)*conjg(q(l1,j,k)))
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = real(ffc(1,j,k))
            at1 = at2*aimag(ffc(1,j,k))
            pot(1,j,k) = at2*q(1,j,k)
            pot(l1,j,k) = zero
            wp = wp + at1*(q(1,j,k)*conjg(q(1,j,k)))
         endif
   20    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at2 = real(ffc(l,1,k))
               at1 = at2*aimag(ffc(l,1,k))
               pot(l,1,k) = at2*q(l,1,k)
               pot(l1,1,k) = at2*q(l1,1,k)
               wp = wp + at1*(q(l,1,k)*conjg(q(l,1,k))
     1                 + q(l1,1,k)*conjg(q(l1,1,k)))
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at2 = real(ffc(1,1,k))
               at1 = at2*aimag(ffc(1,1,k))
               pot(1,1,k) = at2*q(1,1,k)
               pot(l1,1,k) = zero
               wp = wp + at1*(q(1,1,k)*conjg(q(1,1,k)))
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               pot(l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
c keep ky = 0
      if (ks.eq.0) then
         do 70 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at2 = real(ffc(l,j,1))
            at1 = at2*aimag(ffc(l,j,1))
            pot(l,j,1) = at2*q(l,j,1)
            pot(l1,j,1) = at2*q(l1,j,1)
            wp = wp + at1*(q(l,j,1)*conjg(q(l,j,1))
     1              + q(l1,j,1)*conjg(q(l1,j,1)))
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = real(ffc(1,j,1))
            at1 = at2*aimag(ffc(1,j,1))
            pot(1,j,1) = at2*q(1,j,1)
            pot(l1,j,1) = zero
            wp = wp + at1*(q(1,j,1)*conjg(q(1,j,1)))
         endif
   70    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at2 = real(ffc(l,1,1))
            at1 = at2*aimag(ffc(l,1,1))
            pot(l,1,1) = at2*q(l,1,1)
            pot(l1,1,1) = zero
            wp = wp + at1*(q(l,1,1)*conjg(q(l,1,1)))
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            pot(1,1,1) = zero
            pot(l1,1,1) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            pot(l,j,k1) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            pot(l,1,k1) = zero
  110       continue
         endif
      endif
  120 continue
      we = real(nx)*real(ny)*real(nz)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPDIVF32(f,df,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
c this subroutine calculates the divergence in fourier space
c for distributed data with 2D spatial decomposition
c input: all except df, output: df
c approximate flop count is:
c 35*nxc*nyc*nzc + 16*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the divergence is calculated using the equations:
c df(kx,ky,kz) = sqrt(-1)*(kx*fx(kx,ky,kz)+ky*fy(kx,ky,kz)
c                       +kz*fz(kx,ky,kz))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c df(kx=pi) = 0, df(ky=pi) = 0, df(kz=pi) = 0
c and df(kx=0,ky=0,kz=0) = 0.
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      complex f, df
      dimension f(3,nzv,kxyp,kyzp), df(nzv,kxyp,kyzp)
c local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1
      real dnx, dny, dnz, dkx, dky, dkz
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
c find processor id and offsets in y/z
c js/ks = processor co-ordinates in x/y => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
c calculate the divergence
      if (kstrt.gt.(nvpy*nvpz)) return
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         do 20 j = 1, kxyps
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = dkx*f(1,l,j,k) + dky*f(2,l,j,k) + dkz*f(3,l,j,k)
            df(l,j,k) = cmplx(-aimag(zt1),real(zt1))
            zt1 = dkx*f(1,l1,j,k) + dky*f(2,l1,j,k) - dkz*f(3,l1,j,k)
            df(l1,j,k) = cmplx(-aimag(zt1),real(zt1))
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = dkx*f(1,1,j,k) + dky*f(2,1,j,k)
            df(1,j,k) = cmplx(-aimag(zt1),real(zt1))
            df(l1,j,k) = zero
         endif
   20    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               zt1 = dky*f(2,l,1,k) + dkz*f(3,l,1,k)
               df(l,1,k) = cmplx(-aimag(zt1),real(zt1))
               zt1 = dky*f(2,l1,1,k) - dkz*f(3,l1,1,k)
               df(l1,1,k) = cmplx(-aimag(zt1),real(zt1))
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               zt1 = dky*f(2,1,1,k)
               df(1,1,k) = cmplx(-aimag(zt1),real(zt1))
               df(l1,1,k) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               df(l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
c keep ky = 0
      if (ks.eq.0) then
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = dkx*f(1,l,j,1) + dkz*f(3,l,j,1)
            df(l,j,1) = cmplx(-aimag(zt1),real(zt1))
            zt1 = dkx*f(1,l1,j,1) - dkz*f(3,l1,j,1)
            df(l1,j,1) = cmplx(-aimag(zt1),real(zt1))
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = dkx*f(1,1,j,1)
            df(1,j,1) = cmplx(-aimag(zt1),real(zt1))
            df(l1,j,1) = zero
         endif
   70    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = dkz*f(3,l,1,1)
            df(l,1,1) = cmplx(-aimag(zt1),real(zt1))
            df(l1,1,1) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            df(1,1,1) = zero
            df(l1,1,1) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            df(l,j,k1) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            df(l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGRADF32(df,f,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
c this subroutine calculates the gradient in fourier space
c for distributed data with 2D spatial decomposition
c input: all except f, output: f
c approximate flop count is:
c 30*nxc*nyc*nzc + 12*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the gradient is calculated using the equations:
c fx(kx,ky,kz) = sqrt(-1)*kx*df(kx,ky,kz)
c fy(kx,ky,kz) = sqrt(-1)*ky*df(kx,ky,kz)
c fz(kx,ky,kz) = sqrt(-1)*kz*df(kx,ky,kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
c fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
c fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
c fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      complex df, f
      dimension df(nzv,kxyp,kyzp), f(3,nzv,kxyp,kyzp)
c local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1
      real dnx, dny, dnz, dkx, dky, dkz
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
c find processor id and offsets in y/z
c js/ks = processor co-ordinates in x/y => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
c calculate the gradient
      if (kstrt.gt.(nvpy*nvpz)) return
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         do 20 j = 1, kxyps
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = cmplx(-aimag(df(l,j,k)),real(df(l,j,k)))
            f(1,l,j,k) = dkx*zt1
            f(2,l,j,k) = dky*zt1
            f(3,l,j,k) = dkz*zt1
            zt1 = cmplx(-aimag(df(l1,j,k)),real(df(l1,j,k)))
            f(1,l1,j,k) = dkx*zt1
            f(2,l1,j,k) = dky*zt1
            f(3,l1,j,k) = -dkz*zt1
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(df(1,j,k)),real(df(1,j,k)))
            f(1,1,j,k) = dkx*zt1
            f(2,1,j,k) = dky*zt1
            f(3,1,j,k) = zero
            f(1,l1,j,k) = zero
            f(2,l1,j,k) = zero
            f(3,l1,j,k) = zero
         endif
   20    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               zt1 = cmplx(-aimag(df(l,1,k)),real(df(l,1,k)))
               f(1,l,1,k) = zero
               f(2,l,1,k) = dky*zt1
               f(3,l,1,k) = dkz*zt1
               zt1 = cmplx(-aimag(df(l1,1,k)),real(df(l1,1,k)))
               f(1,l1,1,k) = zero
               f(2,l1,1,k) = dky*zt1
               f(3,l1,1,k) = -dkz*zt1
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               zt1 = cmplx(-aimag(df(1,1,k)),real(df(1,1,k)))
               f(1,1,1,k) = zero
               f(2,1,1,k) = dky*zt1
               f(3,1,1,k) = zero
               f(1,l1,1,k) = zero
               f(2,l1,1,k) = zero
               f(3,l1,1,k) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               f(1,l,1,k) = zero
               f(2,l,1,k) = zero
               f(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
c keep ky = 0
      if (ks.eq.0) then
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = cmplx(-aimag(df(l,j,1)),real(df(l,j,1)))
            f(1,l,j,1) = dkx*zt1
            f(2,l,j,1) = zero
            f(3,l,j,1) = dkz*zt1
            zt1 = cmplx(-aimag(df(l1,j,1)),real(df(l1,j,1)))
            f(1,l1,j,1) = dkx*zt1
            f(2,l1,j,1) = zero
            f(3,l1,j,1) = -dkz*zt1
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(df(1,j,1)),real(df(1,j,1)))
            f(1,1,j,1) = dkx*zt1
            f(2,1,j,1) = zero
            f(3,1,j,1) = zero
            f(1,l1,j,1) = zero
            f(2,l1,j,1) = zero
            f(3,l1,j,1) = zero
         endif
   70    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = cmplx(-aimag(df(l,1,1)),real(df(l,1,1)))
            f(1,l,1,1) = zero
            f(2,l,1,1) = zero
            f(3,l,1,1) = dkz*zt1
            f(1,l1,1,1) = zero
            f(2,l1,1,1) = zero
            f(3,l1,1,1) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            f(1,1,1,1) = zero
            f(2,1,1,1) = zero
            f(3,1,1,1) = zero
            f(1,l1,1,1) = zero
            f(2,l1,1,1) = zero
            f(3,l1,1,1) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            f(1,l,j,k1) = zero
            f(2,l,j,k1) = zero
            f(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            f(1,l,1,k1) = zero
            f(2,l,1,k1) = zero
            f(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPSMOOTH32(q,qs,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp, 
     1kyzp,nzhd)
c this subroutine provides a 3d smoothing function,
c in fourier space, with periodic boundary conditions
c for distributed data, with 2D spatial decomposition
c input: q,ffc,nx,ny,nz,kstrt,nzv,kxyp,kyzp,nzhd, output: qs
c approximate flop count is:
c 8*nxc*nyc*nzc + 4*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c smoothing is calculated using the equation:
c qs(kx,ky,kz) = q(kx,ky,kz)*s(kx,ky,kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c qs(kx=pi) = 0, qs(ky=pi) = 0, qs(kz=pi) = 0, and
c qs(kx=0,ky=0,kz=0) = 0.
c q(l,j,k) = complex charge density for fourier mode jj-1,kk-1,l-1
c qs(l,j,k) = complex smoothed charge density,
c aimag(ffc(l,j,k)) = finite-size particle shape factor s
c real(ffc(l,j,k)) = potential green's function g
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
c kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      complex q, qs, ffc
      dimension q(nzv,kxyp,kyzp), qs(nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
c local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1
      real at1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      zero = cmplx(0.0,0.0)
c find processor id and offsets in y/z
c js/ks = processor co-ordinates in x/y => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      joff = joff - 1
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
c calculate smoothing
      if (kstrt.gt.(nvpy*nvpz)) return
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         do 20 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,k))
            qs(l,j,k) = at1*q(l,j,k)
            qs(l1,j,k) = at1*q(l1,j,k)
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,k))
            qs(1,j,k) = at1*q(1,j,k)
            qs(l1,j,k) = zero
         endif
   20    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at1 = aimag(ffc(l,1,k))
               qs(l,1,k) = at1*q(l,1,k)
               qs(l1,1,k) = at1*q(l1,1,k)
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = aimag(ffc(1,1,k))
               qs(1,1,k) = at1*q(1,1,k)
               qs(l1,1,k) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               qs(l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
c keep ky = 0
      if (ks.eq.0) then
         do 70 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,1))
            qs(l,j,1) = at1*q(l,j,1)
            qs(l1,j,1) = at1*q(l1,j,1)
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,1))
            qs(1,j,1) = at1*q(1,j,1)
            qs(l1,j,1) = zero
         endif
   70    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,1,1))
            qs(l,1,1) = at1*q(l,1,1)
            qs(l1,1,1) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,1,1))
            qs(1,1,1) = cmplx(at1*real(q(1,1,1)),0.0)
            qs(l1,1,1) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            qs(l,j,k1) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            qs(l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPRDMODES32(pot,pott,nx,ny,nz,modesx,modesy,modesz,    
     1kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
c this subroutine extracts lowest order modes from packed complex array
c pot and stores them into a location in an unpacked complex array pott
c modes stored: kx = (kxyp*js+(0,1,...kxyp-1)),
c and ky = (kyzp*ks+(0,1,...kyzp-1)), when ks < nvpy/2
c and ky = (kyzp*(ks-nvpy+1)-(kyzp-1,...,1,0)-1), when ks >= nvpy/2
c where js/ks = processor co-ordinates in x/y => idproc = js + nvpy*ks
c and kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
c except kx = NX/2 is stored at location kxyp+1 when js=0,
c and ky = NY/2 is stored at location 1 when ks=nvp/2.
c nx/ny/nz = system length in x/y/z direction
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c modeszd = first dimension of array pott,
c where modeszd  = min(2*modesz-1,nz)
c modesypd = third dimension of array pott, modesypd >= kyzp
c modesxpd = second dimension of array pott,
c modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
c in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, kstrt, nvpy, nvpz, nzv
      integer kxyp, kyzp, modesxpd, modesypd, modeszd
      complex pot, pott
      dimension pot(nzv,kxyp,kyzp), pott(modeszd,modesxpd,modesypd)
c local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, j1, k1, l1, jmax, kmax, kmin, lmax
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      nz2 = nz + 2
      kmax = min0(modesy,nyh)
      kmin = ny - kmax
      lmax = min0(modesz,nzh)
      j1 = kxyp + 1
c find processor id and offsets in y/z
c js/ks = processor co-ordinates in x/y => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      if (kstrt.gt.(nvpy*nvpz)) return
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      jmax = modesx - joff
      if (jmax.gt.kxyps) then
         jmax = kxyps
      else if (jmax.le.0) then
         jmax = 0
      endif
      do 50 k = 1, kyzps
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))
     1 then
         if (k1.gt.nyh) k1 = k1 - ny
         do 20 j = 1, jmax
         if ((j+joff).gt.1) then
            do 10 l = 2, lmax
            l1 = nz2 - l
            pott(2*l-2,j,k) = pot(l,j,k)
            pott(2*l-1,j,k) = pot(l1,j,k)
   10       continue
c mode numbers kz = 0, nz/2
            pott(1,j,k) = pot(1,j,k)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               pott(nz,j,k) = pot(l1,j,k)
            endif
         endif
   20    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c kx = 0
            if (k1.gt.0) then
               do 30 l = 2, lmax
               l1 = nz2 - l
               pott(2*l-2,1,k) = pot(l,1,k)
               pott(2*l-1,1,k) = pot(l1,1,k)
   30          continue
c mode numbers kz = 0, nz/2
               pott(1,1,k) = pot(1,1,k)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(nz,1,k) = pot(l1,1,k)
               endif
c kx = nx/2
            else
               if (modesx.gt.nxh) then
                  do 40 l = 2, lmax
                  l1 = nz2 - l
                  pott(2*l-2,1,k) = pot(l,1,k)
                  pott(2*l-1,1,k) = pot(l1,1,k)
   40             continue
c mode numbers kz = 0, nz/2
                  pott(1,1,k) = pot(1,1,k)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pott(nz,1,k) = pot(l1,1,k)
                  endif
               endif
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
c ky = 0
      if (ks.eq.0) then
         do 70 j = 1, jmax
         if ((j+joff).gt.1) then
            do 60 l = 2, lmax
            l1 = nz2 - l
            pott(2*l-2,j,1) = pot(l,j,1)
            pott(2*l-1,j,1) = pot(l1,j,1)
   60       continue
c mode numbers kz = 0, nz/2
            pott(1,j,1) = pot(1,j,1)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               pott(nz,j,1) = pot(l1,j,1)
            endif
         endif
   70    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c kx = 0
            do 80 l = 2, lmax
            pott(2*l-2,1,1) = pot(l,1,1)
            pott(2*l-1,1,1) = conjg(pot(l,1,1))
   80       continue
c mode numbers kz = 0, nz/2
            pott(1,1,1) = cmplx(real(pot(1,1,1)),0.0)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               pott(nz,1,1) = cmplx(real(pot(l1,1,1)),0.0)
            endif
c kx = nx/2
            if (modesx.gt.nxh) then
               do 90 l = 2, lmax
               l1 = nz2 - l
               pott(2*l-2,j1,1) = conjg(pot(l1,1,1))
               pott(2*l-1,j1,1) = pot(l1,1,1)
   90          continue
c mode numbers kz = 0, nz/2
               pott(1,j1,1) = cmplx(aimag(pot(1,1,1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(nz,j1,1) = cmplx(aimag(pot(l1,1,1)),0.0)
               endif
            endif
         endif
      endif
c ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         if (modesy.gt.nyh) then
            do 110 j = 1, jmax
            if ((j+joff).gt.1) then
               do 100 l = 2, lmax
               l1 = nz2 - l
               pott(2*l-2,j,k1) = pot(l,j,k1)
               pott(2*l-1,j,k1) = pot(l1,j,k1)
  100          continue
c mode numbers kz = 0, nz/2
               pott(1,j,k1) = pot(1,j,k1)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(nz,j,k1) = pot(l1,j,k1)
               endif
            endif
  110       continue
c mode numbers kx = 0, nx/2
            if (js.eq.0) then
c kx = 0
               do 120 l = 2, lmax
               pott(2*l-2,1,k1) = pot(l,1,k1)
               pott(2*l-1,1,k1) = conjg(pot(l,1,k1))
  120          continue
c mode numbers kz = 0, nz/2
               pott(1,1,k1) = cmplx(real(pot(1,1,k1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(nz,1,k1) = cmplx(real(pot(l1,1,k1)),0.0)
               endif
c kx  = nx/2
               if (modesx.gt.nxh) then
                  do 130 l = 2, lmax
                  l1 = nz2 - l
                  pott(2*l-2,j1,k1) = conjg(pot(l1,1,k1))
                  pott(2*l-1,j1,k1) = pot(l1,1,k1)
  130             continue
c mode numbers kz = 0, nz/2
                  pott(1,j1,k1) = cmplx(aimag(pot(1,1,k1)),0.0)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pott(nz,j1,k1) = cmplx(aimag(pot(l1,1,k1)),0.0)
                  endif
               endif
            endif
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPWRMODES32(pot,pott,nx,ny,nz,modesx,modesy,modesz,
     1kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
c this subroutine extracts lowest order modes from a location in an
c unpacked complex array pott and stores them into a packed complex
c array pot
c modes stored: kx = (kxyp*js+(0,1,...kxyp-1)),
c and ky = (kyzp*ks+(0,1,...kyzp-1)), when ks < nvpy/2
c and ky = (kyzp*(ks-nvpy+1)-(kyzp-1,...,1,0)-1), when ks >= nvpy/2
c where js/ks = processor co-ordinates in x/y => idproc = js + nvpy*ks
c and kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
c except kx = NX/2 is stored at location kxyp+1 when js=0,
c and ky = NY/2 is stored at location 1 when ks=nvp/2.
c nx/ny/nz = system length in x/y/z direction
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c modeszd = first dimension of array pott,
c where modeszd  = min(2*modesz-1,nz)
c modesypd = third dimension of array pott, modesypd >= kyzp
c modesxpd = second dimension of array pott,
c modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
c in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, kstrt, nvpy, nvpz, nzv
      integer kxyp, kyzp, modesxpd, modesypd, modeszd
      complex pot, pott
      dimension pot(nzv,kxyp,kyzp), pott(modeszd,modesxpd,modesypd)
c local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, j1, k1, l1, jmax, kmax, kmin, lmax
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      nz2 = nz + 2
      kmax = min0(modesy,nyh)
      kmin = ny - kmax
      lmax = min0(modesz,nzh)
      j1 = kxyp + 1
c find processor id and offsets in y/z
c js/ks = processor co-ordinates in x/y => idproc = js + nvpy*ks
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      joff = kxyp*js
      kxyps = min(kxyp,max(0,nxh-joff))
      koff = kyzp*ks
      kyzps = min(kyzp,max(0,ny-koff))
      koff = koff - 1
      zero = cmplx(0.0,0.0)
      if (kstrt.gt.(nvpy*nvpz)) return
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      jmax = modesx - joff
      if (jmax.gt.kxyps) then
         jmax = kxyps
      else if (jmax.le.0) then
         jmax = 0
      endif
      do 100 k = 1, kyzps
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))
     1 then
         if (k1.gt.nyh) k1 = k1 - ny
         do 30 j = 1, jmax
         if ((j+joff).gt.1) then
            do 10 l = 2, lmax
            l1 = nz2 - l
            pot(l,j,k) = pott(2*l-2,j,k)
            pot(l1,j,k) = pott(2*l-1,j,k)
   10       continue
            do 20 l = lmax+1, nzh
            l1 = nz2 - l
            pot(l,j,k) = zero
            pot(l1,j,k) = zero
   20       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1 
            pot(1,j,k) = pott(1,j,k)
            pot(l1,j,k) = zero
            if (modesz.gt.nzh) then
               pot(l1,j,k) = pott(nz,j,k)
            endif
         endif
   30    continue
         do 50 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 40 l = 1, nz
            pot(l,j,k) = zero
   40       continue
         endif
   50    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c kx = 0
            if (k1.gt.0) then
               do 60 l = 2, lmax
               l1 = nz2 - l
               pot(l,1,k) = pott(2*l-2,1,k)
               pot(l1,1,k) = pott(2*l-1,1,k)
   60          continue
               do 70 l = lmax+1, nzh
               l1 = nz2 - l
               pot(l,1,k) = zero
               pot(l1,1,k) = zero
   70          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               pot(1,1,k) = pott(1,1,k)
               pot(l1,1,k) = zero
               if (modesz.gt.nzh) then
                  pot(l1,1,k) = pott(nz,1,k)
               endif
c kx = nx/2
            else
               do 80 l = 1, nz
               pot(l,1,k) = zero
   80          continue
               if (modesx.gt.nxh) then
                  do 90 l = 2, lmax
                  l1 = nz2 - l
                  pot(l,1,k) = pott(2*l-2,1,k)
                  pot(l1,1,k) = pott(2*l-1,1,k)
   90             continue
c mode numbers kz = 0, nz/2
                  pot(1,1,k) = pott(1,1,k)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pot(l1,1,k) = pott(nz,1,k)
                  endif
               endif
            endif
         endif
      endif
  100 continue
      do 140 k = 1, kyzps
      k1 = k + koff
      if ((k1.ge.kmax).and.(k1.le.kmin).and.(k1.ne.nyh)) then
         do 120 j = 1, kxyps
         if ((j+joff).gt.1) then
            do 110 l = 1, nz
            pot(l,j,k) = zero
  110       continue
         endif
  120    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 130 l = 1, nz
            pot(l,1,k) = zero
  130       continue
         endif
      endif
  140 continue
c mode numbers ky = 0, ny/2
c ky = 0
      if (ks.eq.0) then
         do 170 j = 1, jmax
         if ((j+joff).gt.1) then
            do 150 l = 2, lmax
            l1 = nz2 - l
            pot(l,j,1) = pott(2*l-2,j,1)
            pot(l1,j,1) = pott(2*l-1,j,1)
  150       continue
            do 160 l = lmax+1, nzh
            l1 = nz2 - l
            pot(l,j,1) = zero
            pot(l1,j,1) = zero
  160       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            pot(1,j,1) = pott(1,j,1)
            pot(l1,j,1) = zero
            if (modesz.gt.nzh) then
               pot(l1,j,1) = pott(nz,j,1)
            endif
         endif
  170    continue
         do 190 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 180 l = 1, nz
            pot(l,j,1) = zero
  180       continue
         endif
  190    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c kx = 0
            do 200 l = 2, lmax
            l1 = nz2 - l
            pot(l,1,1) = pott(2*l-2,1,1)
            pot(l1,1,1) = zero
  200       continue
            do 210 l = lmax+1, nzh
            l1 = nz2 - l
            pot(l,1,1) = zero
            pot(l1,1,1) = zero
  210       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            pot(1,1,1) = cmplx(real(pott(1,1,1)),0.0)
            pot(l1,1,1) = zero
            if (modesz.gt.nzh) then
               pot(l1,1,1) = cmplx(real(pott(nz,1,1)),0.0)
            endif
c kx = nx/2
            if (modesx.gt.nxh) then
               do 220 l = 2, lmax
               l1 = nz2 - l
               pot(l1,1,1) = conjg(pott(2*l-2,j1,1))
  220          continue
c mode numbers kz = 0, nz/2
               pot(1,1,1) = cmplx(real(pot(1,1,1)),real(pott(1,j1,1)))
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pot(l1,1,1) = cmplx(real(pot(l1,1,1)),
     1                                real(pott(nz,j1,1)))
               endif
            endif
         endif
      endif
c ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 240 j = 1, jmax
         if ((j+joff).gt.1) then
            do 230 l = 1, nz
            pot(l,j,k1) = zero
  230       continue
         endif
  240    continue
         do 260 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 250 l = 1, nz
            pot(l,j,k1) = zero
  250       continue
         endif
  260    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 270 l = 1, nz
            pot(l,1,k1) = zero
  270       continue
         endif
         if (modesy.gt.nyh) then
            do 290 j = 1, jmax
            if ((j+joff).gt.1) then
               do 280 l = 2, lmax
               l1 = nz2 - l
               pot(l,j,k1) = pott(2*l-2,j,k1)
               pot(l1,j,k1) = pott(2*l-1,j,k1)
  280          continue
c mode numbers kz = 0, nz/2
               pot(1,j,k1) = pott(1,j,k1)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pot(l1,j,k1) = pott(nz,j,k1)
               endif
            endif
  290       continue
c mode numbers kx = 0, nx/2
            if (js.eq.0) then
c kx = 0
               do 300 l = 2, lmax
               pot(l,1,k1) = pott(2*l-2,1,k1)
  300          continue
c mode numbers kz = 0, nz/2
               pot(1,1,k1) = cmplx(real(pott(1,1,k1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pot(l1,1,k1) = cmplx(real(pott(nz,1,k1)),0.0)
               endif
c kx  = nx/2
               if (modesx.gt.nxh) then
                  do 310 l = 2, lmax
                  l1 = nz2 - l
                  pot(l1,1,k1) = conjg(pott(2*l-2,j1,k1))
  310             continue
c mode numbers kz = 0, nz/2
                  pot(1,1,k1) = cmplx(real(pot(1,1,k1)),
     1                                real(pott(1,j1,k1)))
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pot(l1,1,k1) = cmplx(real(pot(l1,1,k1)),
     1                                    real(pott(nz,j1,k1)))
                  endif
               endif
            endif
         endif
      endif
      return
      end
