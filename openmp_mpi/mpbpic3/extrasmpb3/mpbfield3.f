c Fortran library for Skeleton 3D Electromagnetic MPI/OpenMP PIC Code
c field diagnostics
c written by viktor k. decyk, ucla
c copyright 1994-2015, regents of the university of california
c-----------------------------------------------------------------------
      subroutine MPPOTP32(q,pot,ffc,we,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp
     1,kyzp,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c potential, with periodic boundary conditions for distributed data,
c with 2D spatial decomposition and OpenMP
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
      integer kxyps, kyzps, k1, l1, kk
      real at1, at2
      complex zero
      double precision wp, sum1, sum2
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
      sum1 = 0.0d0
      if (kstrt.gt.(nvpy*nvpz)) go to 120
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,at1,at2,wp) REDUCTION(+:sum1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      wp = 0.0d0
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
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
      endif
      sum1 = sum1 + wp
   20 continue
!$OMP END DO
!$OMP END PARALLEL
c mode numbers kx = 0, nx/2
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,at1,at2,wp) REDUCTION(+:sum2)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         wp = 0.0d0
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
         sum2 = sum2 + wp
      endif
   50 continue
!$OMP END PARALLEL DO
      sum1 = sum1 + sum2
c mode numbers ky = 0, ny/2
      sum2 = 0.0d0
c keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,at1,at2,wp) REDUCTION(+:sum2)
         do 70 j = 1, kxyps
         wp = 0.0d0
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
         sum2 = sum2 + wp
   70    continue
!$OMP END PARALLEL DO
         wp = 0.0d0
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
         sum2 = sum2 + wp
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
      sum1 = sum1 + sum2
  120 continue
      we = real(nx)*real(ny)*real(nz)*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPDIVF32(f,df,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
c this subroutine calculates the divergence in fourier space
c for distributed data with 2D spatial decomposition and OpenMP
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
      integer kxyps, kyzps, k1, l1, kk
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
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,zt1)
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
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
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
c mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dkz,zt1)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
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
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
c keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkz,zt1)
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
!$OMP END PARALLEL DO
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
      subroutine MPPGRADF32(df,f,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
c this subroutine calculates the gradient in fourier space
c for distributed data with 2D spatial decomposition and OpenMP
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
      integer kxyps, kyzps, k1, l1, kk
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
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,zt1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
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
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
c mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dkz,zt1)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
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
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
c keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkz,zt1)
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
!$OMP END PARALLEL DO
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
      subroutine MPPCURLF32(f,g,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp)
c this subroutine calculates the curl in fourier space
c for distributed data with 2D spatial decomposition and OpenMP
c input: all except g, output: g
c approximate flop count is:
c 86*nxc*nyc*nzc + 32*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the curl is calculated using the equations:
c gx(kx,ky,kz) = sqrt(-1)*(ky*fz(kx,ky,kz)-kz*fy(kx,ky,kz))
c gy(kx,ky,kz) = sqrt(-1)*(kz*fx(kx,ky,kz)-kx*fz(kx,ky,kz))
c gz(kx,ky,kz) = sqrt(-1)*(kx*fy(kx,ky,kz)-ky*fx(kx,ky,kz))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c gx(kx=pi) = gy(kx=pi) = gz(kx=pi) = 0,
c gx(ky=pi) = gy(ky=pi) = gx(ky=pi) = 0,
c gx(kz=pi) = gy(kz=pi) = gz(kz=pi) = 0,
c gx(kx=0,ky=0,kz=0) = gy(kx=0,ky=0,kz=0) = gz(kx=0,ky=0,kz=0) = 0.
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      complex f, g
      dimension f(3,nzv,kxyp,kyzp), g(3,nzv,kxyp,kyzp)
c local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz
      complex zero, zt1, zt2, zt3
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
c calculate the curl
      if (kstrt.gt.(nvpy*nvpz)) return
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,zt1,zt2,zt3)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = cmplx(-aimag(f(3,l,j,k)),real(f(3,l,j,k)))
            zt2 = cmplx(-aimag(f(2,l,j,k)),real(f(2,l,j,k)))
            zt3 = cmplx(-aimag(f(1,l,j,k)),real(f(1,l,j,k)))
            g(1,l,j,k) = dky*zt1 - dkz*zt2
            g(2,l,j,k) = dkz*zt3 - dkx*zt1
            g(3,l,j,k) = dkx*zt2 - dky*zt3
            zt1 = cmplx(-aimag(f(3,l1,j,k)),real(f(3,l1,j,k)))
            zt2 = cmplx(-aimag(f(2,l1,j,k)),real(f(2,l1,j,k)))
            zt3 = cmplx(-aimag(f(1,l1,j,k)),real(f(1,l1,j,k)))
            g(1,l1,j,k) = dky*zt1 + dkz*zt2
            g(2,l1,j,k) = -dkz*zt3 - dkx*zt1
            g(3,l1,j,k) = dkx*zt2 - dky*zt3
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(f(3,1,j,k)),real(f(3,1,j,k)))
            zt2 = cmplx(-aimag(f(2,1,j,k)),real(f(2,1,j,k)))
            zt3 = cmplx(-aimag(f(1,1,j,k)),real(f(1,1,j,k)))
            g(1,1,j,k) = dky*zt1
            g(2,1,j,k) = -dkx*zt1
            g(3,1,j,k) = dkx*zt2 - dky*zt3
            g(1,l1,j,k) = zero
            g(2,l1,j,k) = zero
            g(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
c mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dkz,zt1,zt2,zt3)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         if (js.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               zt1 = cmplx(-aimag(f(3,l,1,k)),real(f(3,l,1,k)))
               zt2 = cmplx(-aimag(f(2,l,1,k)),real(f(2,l,1,k)))
               zt3 = cmplx(-aimag(f(1,l,1,k)),real(f(1,l,1,k)))
               g(1,l,1,k) = dky*zt1 - dkz*zt2
               g(2,l,1,k) = dkz*zt3
               g(3,l,1,k) = -dky*zt3
               zt1 = cmplx(-aimag(f(3,l1,1,k)),real(f(3,l1,1,k)))
               zt2 = cmplx(-aimag(f(2,l1,1,k)),real(f(2,l1,1,k)))
               zt3 = cmplx(-aimag(f(1,l1,1,k)),real(f(1,l1,1,k)))
               g(1,l1,1,k) = dky*zt1 + dkz*zt2
               g(2,l1,1,k) = -dkz*zt3
               g(3,l1,1,k) = -dky*zt3
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               zt1 = cmplx(-aimag(f(3,1,1,k)),real(f(3,1,1,k)))
               zt3 = cmplx(-aimag(f(1,1,1,k)),real(f(1,1,1,k)))
               g(1,1,1,k) = dky*zt1
               g(2,1,1,k) = zero
               g(3,1,1,k) = -dky*zt3
               g(1,l1,1,k) = zero
               g(2,l1,1,k) = zero
               g(3,l1,1,k) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               g(1,l,1,k) = zero
               g(2,l,1,k) = zero
               g(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
c keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkz,zt1,zt2,zt3)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt1 = cmplx(-aimag(f(3,l,j,1)),real(f(3,l,j,1)))
            zt2 = cmplx(-aimag(f(2,l,j,1)),real(f(2,l,j,1)))
            zt3 = cmplx(-aimag(f(1,l,j,1)),real(f(1,l,j,1)))
            g(1,l,j,1) = -dkz*zt2
            g(2,l,j,1) = dkz*zt3 - dkx*zt1
            g(3,l,j,1) = dkx*zt2
            zt1 = cmplx(-aimag(f(3,l1,j,1)),real(f(3,l1,j,1)))
            zt2 = cmplx(-aimag(f(2,l1,j,1)),real(f(2,l1,j,1)))
            zt3 = cmplx(-aimag(f(1,l1,j,1)),real(f(1,l1,j,1)))
            g(1,l1,j,1) = dkz*zt2
            g(2,l1,j,1) = -dkz*zt3 - dkx*zt1
            g(3,l1,j,1) = dkx*zt2
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(f(3,1,j,1)),real(f(3,1,j,1)))
            zt2 = cmplx(-aimag(f(2,1,j,1)),real(f(2,1,j,1)))
            g(1,1,j,1) = zero
            g(2,1,j,1) = -dkx*zt1
            g(3,1,j,1) = dkx*zt2
            g(1,l1,j,1) = zero
            g(2,l1,j,1) = zero
            g(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            zt2 = cmplx(-aimag(f(2,l,1,1)),real(f(2,l,1,1)))
            zt3 = cmplx(-aimag(f(1,l,1,1)),real(f(1,l,1,1)))
            g(1,l,1,1) = -dkz*zt2
            g(2,l,1,1) = dkz*zt3
            g(3,l,1,1) = zero
            g(1,l1,1,1) = zero
            g(2,l1,1,1) = zero
            g(3,l1,1,1) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            g(1,1,1,1) = zero
            g(2,1,1,1) = zero
            g(3,1,1,1) = zero
            g(1,l1,1,1) = zero
            g(2,l1,1,1) = zero
            g(3,l1,1,1) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            g(1,l,j,k1) = zero
            g(2,l,j,k1) = zero
            g(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            g(1,l,1,k1) = zero
            g(2,l,1,k1) = zero
            g(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPAVPOT332(bxyz,axyz,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp
     1,kyzp)
c this subroutine calculates 3d vector potential from magnetic field
c in fourier space with periodic boundary conditions,
c for distributed data, with 2D spatial decomposition and OpenMP
c input: bxyz,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp
c output: axyz
c approximate flop count is:
c 99*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the vector potential is calculated using the equations:
c ax(kx,ky,kz) = sqrt(-1)*
c                (ky*bz(kx,ky,kz)-kz*by(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c ay(kx,ky,kz) = sqrt(-1)*
c                (kz*bx(kx,ky,kz)-kx*bz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c az(kx,ky,kz) = sqrt(-1)*
c                (kx*by(kx,ky,kz)-ky*bx(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = ax(ky=pi) = 0,
c ax(kz=pi) = ay(kz=pi) = az(kz=pi) = 0,
c ax(kx=0,ky=0,kz=0) = ay(kx=0,ky=0,kz=0) = az(kx=0,ky=0,kz=0) = 0.
c bxyz(i,l,j,k) = i component of complex magnetic field
c axyz(i,l,j,k) = i component of complex vector potential
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
c kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp
      complex bxyz, axyz
      dimension bxyz(3,nzv,kxyp,kyzp), axyz(3,nzv,kxyp,kyzp)
c local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, dkx, dky, dkz, dkx2, dky2, dkxy2
      real at1, at2, at3, at4
      complex zero, zt1, zt2, zt3
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
c calculate vector potential
      if (kstrt.gt.(nvpy*nvpz)) return
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,dky2,dkxy2,at1,at2,at3,at4, 
!$OMP& zt1,zt2,zt3)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         dkx = dnx*real(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkxy2)
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dkz*at1
            zt1 = cmplx(-aimag(bxyz(3,l,j,k)),real(bxyz(3,l,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,l,j,k)),real(bxyz(2,l,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,l,j,k)),real(bxyz(1,l,j,k)))
            axyz(1,l,j,k) = at3*zt1 - at4*zt2
            axyz(2,l,j,k) = at4*zt3 - at2*zt1
            axyz(3,l,j,k) = at2*zt2 - at3*zt3
            zt1 = cmplx(-aimag(bxyz(3,l1,j,k)),real(bxyz(3,l1,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,l1,j,k)),real(bxyz(2,l1,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,l1,j,k)),real(bxyz(1,l1,j,k)))
            axyz(1,l1,j,k) = at3*zt1 + at4*zt2
            axyz(2,l1,j,k) = -at4*zt3 - at2*zt1
            axyz(3,l1,j,k) = at2*zt2 - at3*zt3
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = 1.0/dkxy2
            at2 = dkx*at1
            at3 = dky*at1
            zt1 = cmplx(-aimag(bxyz(3,1,j,k)),real(bxyz(3,1,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,1,j,k)),real(bxyz(2,1,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,1,j,k)),real(bxyz(1,1,j,k)))
            axyz(1,1,j,k) = at3*zt1
            axyz(2,1,j,k) = -at2*zt1
            axyz(3,1,j,k) = at2*zt2 - at3*zt3
            axyz(1,l1,j,k) = zero
            axyz(2,l1,j,k) = zero
            axyz(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
c mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dky2,dkz,at1,at3,at4,zt1,zt2,zt3
!$OMP& )
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         if (js.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               at1 = 1.0/(dkz*dkz + dky2)
               at3 = dky*at1
               at4 = dkz*at1
               zt1 = cmplx(-aimag(bxyz(3,l,1,k)),real(bxyz(3,l,1,k)))
               zt2 = cmplx(-aimag(bxyz(2,l,1,k)),real(bxyz(2,l,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,l,1,k)),real(bxyz(1,l,1,k)))
               axyz(1,l,1,k) = at3*zt1 - at4*zt2
               axyz(2,l,1,k) = at4*zt3
               axyz(3,l,1,k) = -at3*zt3
               zt1 = cmplx(-aimag(bxyz(3,l1,1,k)),real(bxyz(3,l1,1,k)))
               zt2 = cmplx(-aimag(bxyz(2,l1,1,k)),real(bxyz(2,l1,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,l1,1,k)),real(bxyz(1,l1,1,k)))
               axyz(1,l1,1,k) = at3*zt1 + at4*zt2
               axyz(2,l1,1,k) = -at4*zt3
               axyz(3,l1,1,k) = -at3*zt3
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at3 = 1.0/dky
               zt1 = cmplx(-aimag(bxyz(3,1,1,k)),real(bxyz(3,1,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,1,1,k)),real(bxyz(1,1,1,k)))
               axyz(1,1,1,k) = at3*zt1
               axyz(2,1,1,k) = zero
               axyz(3,1,1,k) = -at3*zt3
               axyz(1,l1,1,k) = zero
               axyz(2,l1,1,k) = zero
               axyz(3,l1,1,k) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               axyz(1,l,1,k) = zero
               axyz(2,l,1,k) = zero
               axyz(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
c keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkx2,dkz,at1,at2,at4,zt1,zt2,zt3)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         dkx2 = dkx*dkx
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkx2)
            at2 = dkx*at1
            at4 = dkz*at1
            zt1 = cmplx(-aimag(bxyz(3,l,j,1)),real(bxyz(3,l,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,l,j,1)),real(bxyz(2,l,j,1)))
            zt3 = cmplx(-aimag(bxyz(1,l,j,1)),real(bxyz(1,l,j,1)))
            axyz(1,l,j,1) = -at4*zt2
            axyz(2,l,j,1) = at4*zt3 - at2*zt1
            axyz(3,l,j,1) = at2*zt2
            zt1 = cmplx(-aimag(bxyz(3,l1,j,1)),real(bxyz(3,l1,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,l1,j,1)),real(bxyz(2,l1,j,1)))
            zt3 = cmplx(-aimag(bxyz(1,l1,j,1)),real(bxyz(1,l1,j,1)))
            axyz(1,l1,j,1) = at4*zt2
            axyz(2,l1,j,1) = -at4*zt3 - at2*zt1
            axyz(3,l1,j,1) = at2*zt2
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = 1.0/dkx
            zt1 = cmplx(-aimag(bxyz(3,1,j,1)),real(bxyz(3,1,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,1,j,1)),real(bxyz(2,1,j,1)))
            axyz(1,1,j,1) = zero
            axyz(2,1,j,1) = -at2*zt1
            axyz(3,1,j,1) = at2*zt2
            axyz(1,l1,j,1) = zero
            axyz(2,l1,j,1) = zero
            axyz(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at4 = 1.0/dkz
            zt2 = cmplx(-aimag(bxyz(2,l,1,1)),real(bxyz(2,l,1,1)))
            zt3 = cmplx(-aimag(bxyz(1,l,1,1)),real(bxyz(1,l,1,1)))
            axyz(1,l,1,1) = -at4*zt2
            axyz(2,l,1,1) = at4*zt3
            axyz(3,l,1,1) = zero
            axyz(1,l1,1,1) = zero
            axyz(2,l1,1,1) = zero
            axyz(3,l1,1,1) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            axyz(1,1,1,1) = zero
            axyz(2,1,1,1) = zero
            axyz(3,1,1,1) = zero
            axyz(1,l1,1,1) = zero
            axyz(2,l1,1,1) = zero
            axyz(3,l1,1,1) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            axyz(1,l,j,k1) = zero
            axyz(2,l,j,k1) = zero
            axyz(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            axyz(1,l,1,k1) = zero
            axyz(2,l,1,k1) = zero
            axyz(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPAVRPOT332(axyz,bxyz,ffc,affp,ci,nx,ny,nz,kstrt,nvpy,
     1nvpz,nzv,kxyp,kyzp,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for the
c radiative part of the vector potential
c with periodic boundary conditions,
c for distributed data, with 2D spatial decomposition and OpenMP
c input: axyz,bxyz,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,kyzp
c output: axyz
c approximate flop count is:
c 105*nxc*nyc*nzc + 42*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the vector potential is calculated using the equations:
c ax(kx,ky,kz) = sqrt(-1)*
c (ky*bz(kx,ky,kz)-kz*by(kx,ky,kz) - affp*ci2*cux(kx,ky,kz)*s(kx,ky,kz))
c /(kx*kx+ky*ky+kz*kz)
c ay(kx,ky,kz) = sqrt(-1)*
c (kz*bx(kx,ky,kz)-kx*bz(kx,ky,kz) + affp*ci2*cuy(kx,ky,kz)*s(kx,ky,kz))
c /(kx*kx+ky*ky+kz*kz)
c az(kx,ky,kz) = sqrt(-1)*
c (kx*by(kx,ky,kz)-ky*bx(kx,ky,kz) - affp*ci2*cuz(kx,ky,kz)*s(kx,ky,kz))
c /(kx*kx+ky*ky+kz*kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = ax(ky=pi) = 0,
c ax(kz=pi) = ay(kz=pi) = az(kz=pi) = 0,
c ax(kx=0,ky=0,kz=0) = ay(kx=0,ky=0,kz=0) = az(kx=0,ky=0,kz=0) = 0.
c axyz(i,l,j,k) = on entry, complex current density cu
c axyz(i,l,j,k) = on exit, complex current radiative vector potential
c bxyz(i,l,j,k) = complex magnetic field
c aimag(ffc(l,j,k)) = finite-size particle shape factor s
c real(ffc(l,j,k)) = potential green's function g
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
c kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
c kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c ci = reciprical of velocity of light
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      real affp, ci
      complex axyz, bxyz, ffc
      dimension axyz(3,nzv,kxyp,kyzp), bxyz(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
c local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
      real dnx, dny, dnz, afc2, dkx, dky, dkz, dkx2, dky2, dkxy2
      real at1, at2
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      afc2 = affp*ci*ci
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
c calculate the radiative vector potential
      if (kstrt.gt.(nvpy*nvpz)) return
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,dkx,dky,dkz,dky2,dkxy2,at1,at2,zt1,zt2,
!$OMP& zt3)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         dkx = dnx*real(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkxy2)
            at2 = afc2*aimag(ffc(l,j,k))
            zt1 = cmplx(-aimag(bxyz(3,l,j,k)),real(bxyz(3,l,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,l,j,k)),real(bxyz(2,l,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,l,j,k)),real(bxyz(1,l,j,k)))
            axyz(1,l,j,k) = at1*(dky*zt1 - dkz*zt2 - at2*axyz(1,l,j,k))
            axyz(2,l,j,k) = at1*(dkz*zt3 - dkx*zt1 - at2*axyz(2,l,j,k))
            axyz(3,l,j,k) = at1*(dkx*zt2 - dky*zt3 - at2*axyz(3,l,j,k))
            zt1 = cmplx(-aimag(bxyz(3,l1,j,k)),real(bxyz(3,l1,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,l1,j,k)),real(bxyz(2,l1,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,l1,j,k)),real(bxyz(1,l1,j,k)))
            axyz(1,l1,j,k) = at1*(dky*zt1 + dkz*zt2
     1                     - at2*axyz(1,l1,j,k))
            axyz(2,l1,j,k) = -at1*(dkz*zt3 + dkx*zt1
     1                     + at2*axyz(2,l1,j,k))
            axyz(3,l1,j,k) = at1*(dkx*zt2 - dky*zt3
     1                     - at2*axyz(3,l1,j,k))
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = 1.0/dkxy2
            at2 = afc2*aimag(ffc(1,j,k))
            zt1 = cmplx(-aimag(bxyz(3,1,j,k)),real(bxyz(3,1,j,k)))
            zt2 = cmplx(-aimag(bxyz(2,1,j,k)),real(bxyz(2,1,j,k)))
            zt3 = cmplx(-aimag(bxyz(1,1,j,k)),real(bxyz(1,1,j,k)))
            axyz(1,1,j,k) = at1*(dky*zt1 - at2*axyz(1,1,j,k))
            axyz(2,1,j,k) = -at1*(dkx*zt1 + at2*axyz(2,1,j,k))
            axyz(3,1,j,k) = at1*(dkx*zt2 - dky*zt3 - at2*axyz(3,1,j,k))
            axyz(1,l1,j,k) = zero
            axyz(2,l1,j,k) = zero
            axyz(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
c mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,dky,dkz,dky2,at1,at2,zt1,zt2,zt3)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*real(k1)
         dky2 = dky*dky
         if (js.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*real(l - 1)
               at1 = 1.0/(dkz*dkz + dky2)
               at2 = afc2*aimag(ffc(l,1,k))
               zt1 = cmplx(-aimag(bxyz(3,l,1,k)),real(bxyz(3,l,1,k)))
               zt2 = cmplx(-aimag(bxyz(2,l,1,k)),real(bxyz(2,l,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,l,1,k)),real(bxyz(1,l,1,k)))
               axyz(1,l,1,k) = at1*(dky*zt1 - dkz*zt2 
     1                       - at2*axyz(1,l,1,k))
               axyz(2,l,1,k) = at1*(dkz*zt3 - at2*axyz(2,l,1,k))
               axyz(3,l,1,k) = -at1*(dky*zt3 + at2*axyz(3,l,1,k))
               zt1 = cmplx(-aimag(bxyz(3,l1,1,k)),real(bxyz(3,l1,1,k)))
               zt2 = cmplx(-aimag(bxyz(2,l1,1,k)),real(bxyz(2,l1,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,l1,1,k)),real(bxyz(1,l1,1,k)))
               axyz(1,l1,1,k) = at1*(dky*zt1 + dkz*zt2
     1                        - at2*axyz(1,l1,1,k))
               axyz(2,l1,1,k) = -at1*(dkz*zt3 + at2*axyz(2,l1,1,k))
               axyz(3,l1,1,k) = -at1*(dky*zt3 + at2*axyz(3,l1,1,k))
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = 1.0/(dky*dky)
               at2 = afc2*aimag(ffc(1,1,k))
               zt1 = cmplx(-aimag(bxyz(3,1,1,k)),real(bxyz(3,1,1,k)))
               zt3 = cmplx(-aimag(bxyz(1,1,1,k)),real(bxyz(1,1,1,k)))
               axyz(1,1,1,k) = at1*(dky*zt1 - at2*axyz(1,1,1,k))
               axyz(2,1,1,k) = zero
               axyz(3,1,1,k) = -at1*(dky*zt3 + at2*axyz(3,1,1,k))
               axyz(1,l1,1,k) = zero
               axyz(2,l1,1,k) = zero
               axyz(3,l1,1,k) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               axyz(1,l,1,k) = zero
               axyz(2,l,1,k) = zero
               axyz(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
c keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,dkx,dkz,dkx2,at1,at2,zt1,zt2,zt3)
         do 70 j = 1, kxyps
         dkx = dnx*real(j + joff)
         dkx2 = dkx*dkx
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz + dkx2)
            at2 = afc2*aimag(ffc(l,j,1))
            zt1 = cmplx(-aimag(bxyz(3,l,j,1)),real(bxyz(3,l,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,l,j,1)),real(bxyz(2,l,j,1)))
            zt3 = cmplx(-aimag(bxyz(1,l,j,1)),real(bxyz(1,l,j,1)))
            axyz(1,l,j,1) = -at1*(dkz*zt2 + at2*axyz(1,l,j,1))
            axyz(2,l,j,1) = at1*(dkz*zt3 - dkx*zt1 - at2*axyz(2,l,j,1))
            axyz(3,l,j,1) = at1*(dkx*zt2 - at2*axyz(3,l,j,1))
            zt1 = cmplx(-aimag(bxyz(3,l1,j,1)),real(bxyz(3,l1,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,l1,j,1)),real(bxyz(2,l1,j,1)))
            zt3 = cmplx(-aimag(bxyz(1,l1,j,1)),real(bxyz(1,l1,j,1)))
            axyz(1,l1,j,1) = at1*(dkz*zt2 - at2*axyz(1,l1,j,1))
            axyz(2,l1,j,1) = -at1*(dkz*zt3 + dkx*zt1 
     1                     + at2*axyz(2,l1,j,1))
            axyz(3,l1,j,1) = at1*(dkx*zt2 - at2*axyz(3,l1,j,1))
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = 1.0/(dkx*dkx)
            at2 = afc2*aimag(ffc(1,j,1))
            zt1 = cmplx(-aimag(bxyz(3,1,j,1)),real(bxyz(3,1,j,1)))
            zt2 = cmplx(-aimag(bxyz(2,1,j,1)),real(bxyz(2,1,j,1)))
            axyz(1,1,j,1) = zero
            axyz(2,1,j,1) = -at1*(dkx*zt1 + at2*axyz(2,1,j,1))
            axyz(3,1,j,1) = at1*(dkx*zt2 - at2*axyz(3,1,j,1))
            axyz(1,l1,j,1) = zero
            axyz(2,l1,j,1) = zero
            axyz(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*real(l - 1)
            at1 = 1.0/(dkz*dkz)
            at2 = afc2*aimag(ffc(l,1,1))
            zt2 = cmplx(-aimag(bxyz(2,l,1,1)),real(bxyz(2,l,1,1)))
            zt3 = cmplx(-aimag(bxyz(1,l,1,1)),real(bxyz(1,l,1,1)))
            axyz(1,l,1,1) = -at1*(dkz*zt2 + at2*axyz(1,l,1,1))
            axyz(2,l,1,1) = at1*(dkz*zt3 - at2*axyz(2,l,1,1))
            axyz(3,l,1,1) = zero
            axyz(1,l1,1,1) = zero
            axyz(2,l1,1,1) = zero
            axyz(3,l1,1,1) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            axyz(1,1,1,1) = zero
            axyz(2,1,1,1) = zero
            axyz(3,1,1,1) = zero
            axyz(1,l1,1,1) = zero
            axyz(2,l1,1,1) = zero
            axyz(3,l1,1,1) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            axyz(1,l,j,k1) = zero
            axyz(2,l,j,k1) = zero
            axyz(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            axyz(1,l,1,k1) = zero
            axyz(2,l,1,k1) = zero
            axyz(3,l,1,k1) = zero
  110       continue
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPSMOOTH32(q,qs,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,kxyp,
     1kyzp,nzhd)
c this subroutine provides a 3d smoothing function,
c in fourier space, with periodic boundary conditions
c for distributed data, with 2D spatial decomposition and OpenMP
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
      integer kxyps, kyzps, k1, l1, kk
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
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,at1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
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
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
c mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,at1)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
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
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
c keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,at1)
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
!$OMP END PARALLEL DO
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
      subroutine MPPSMOOTH332(cu,cus,ffc,nx,ny,nz,kstrt,nvpy,nvpz,nzv,  
     1kxyp,kyzp,nzhd)
c this subroutine provides a 3d vector smoothing function
c in fourier space, with periodic boundary conditions
c or distributed data, with 2D spatial decomposition and OpenMP
c input: cu,ffc,isign,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd
c output: cus
c approximate flop count is:
c 24*nxc*nyc*nzc + 12*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c smoothing is calculated using the equation:
c cusx(kx,ky,kz) = cux(kx,ky,kz)*s(kx,ky,kz)
c cusy(kx,ky,kz) = cuy(kx,ky,kz)*s(kx,ky,kz)
c cusz(kx,ky,kz) = cuz(kx,ky,kz)*s(kx,ky,kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c cusx(kx=pi) = cusy(kx=pi) = cusz(kx=pi) = 0,
c cusx(ky=pi) = cusy(ky=pi) = cusx(ky=pi) = 0,
c cusx(kz=pi) = cusy(kz=pi) = cusz(kz=pi) = 0,
c cusx(kx=0,ky=0,kz=0) = cusy(kx=0,ky=0,kz=0) = cusz(kx=0,ky=0,kz=0) = 0
c cu(l,j,k,m) = complex current density for fourier mode jj-1,kk-1,l-1
c cus(1,l,j,k,m) = x component of complex smoothed current density
c cus(2,l,j,k,m) = y component of complex smoothed current density
c cus(3,l,j,k,m) = z component of complex smoothed current density
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*js and
c kk = k + kyzp*ks, and MPI rank idproc = js + nvpy*ks
c aimag(ffc(l,j,k)) = finite-size particle shape factor s
c real(ffc(l,j,k)) = potential green's function g
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, nvpy, nvpz, nzv, kxyp, kyzp, nzhd
      complex cu, cus, ffc
      dimension cu(3,nzv,kxyp,kyzp), cus(3,nzv,kxyp,kyzp)
      dimension ffc(nzhd,kxyp,kyzp)
c local data
      integer j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
      integer kxyps, kyzps, k1, l1, kk
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
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,kk,k1,l1,at1)
      do 20 kk = 1, kyzps*kxyps
      k = (kk - 1)/kxyps
      j = kk - kxyps*k
      k = k + 1
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,k))
            cus(1,l,j,k) = at1*cu(1,l,j,k)
            cus(2,l,j,k) = at1*cu(2,l,j,k)
            cus(3,l,j,k) = at1*cu(3,l,j,k)
            cus(1,l1,j,k) = at1*cu(1,l1,j,k)
            cus(2,l1,j,k) = at1*cu(2,l1,j,k)
            cus(3,l1,j,k) = at1*cu(3,l1,j,k)
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,k))
            cus(1,1,j,k) = at1*cu(1,1,j,k)
            cus(2,1,j,k) = at1*cu(2,1,j,k)
            cus(3,1,j,k) = at1*cu(3,1,j,k)
            cus(1,l1,j,k) = zero
            cus(2,l1,j,k) = zero
            cus(3,l1,j,k) = zero
         endif
      endif
   20 continue
!$OMP END DO
!$OMP END PARALLEL
c mode numbers kx = 0, nx/2
!$OMP PARALLEL DO PRIVATE(k,l,k1,l1,at1)
      do 50 k = 1, kyzps
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         if (js.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at1 = aimag(ffc(l,1,k))
               cus(1,l,1,k) = at1*cu(1,l,1,k)
               cus(2,l,1,k) = at1*cu(2,l,1,k)
               cus(3,l,1,k) = at1*cu(3,l,1,k)
               cus(1,l1,1,k) = at1*cu(1,l1,1,k)
               cus(2,l1,1,k) = at1*cu(2,l1,1,k)
               cus(3,l1,1,k) = at1*cu(3,l1,1,k)
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = aimag(ffc(1,1,k))
               cus(1,1,1,k) = at1*cu(1,1,1,k)
               cus(2,1,1,k) = at1*cu(2,1,1,k)
               cus(3,1,1,k) = at1*cu(3,1,1,k)
               cus(1,l1,1,k) = zero
               cus(2,l1,1,k) = zero
               cus(3,l1,1,k) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               cus(1,l,1,k) = zero
               cus(2,l,1,k) = zero
               cus(3,l,1,k) = zero
   40          continue
            endif
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
c keep ky = 0
      if (ks.eq.0) then
!$OMP PARALLEL DO PRIVATE(j,l,l1,at1)
         do 70 j = 1, kxyps
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,1))
            cus(1,l,j,1) = at1*cu(1,l,j,1)
            cus(2,l,j,1) = at1*cu(2,l,j,1)
            cus(3,l,j,1) = at1*cu(3,l,j,1)
            cus(1,l1,j,1) = at1*cu(1,l1,j,1)
            cus(2,l1,j,1) = at1*cu(2,l1,j,1)
            cus(3,l1,j,1) = at1*cu(3,l1,j,1)
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,1))
            cus(1,1,j,1) = at1*cu(1,1,j,1)
            cus(2,1,j,1) = at1*cu(2,1,j,1)
            cus(3,1,j,1) = at1*cu(3,1,j,1)
            cus(1,l1,j,1) = zero
            cus(2,l1,j,1) = zero
            cus(3,l1,j,1) = zero
         endif
   70    continue
!$OMP END PARALLEL DO
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,1,1))
            cus(1,l,1,1) = at1*cu(1,l,1,1)
            cus(2,l,1,1) = at1*cu(2,l,1,1)
            cus(3,l,1,1) = at1*cu(3,l,1,1)
            cus(1,l1,1,1) = zero
            cus(2,l1,1,1) = zero
            cus(3,l1,1,1) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,1,1))
            cus(1,1,1,1) = cmplx(at1*real(cu(1,1,1,1)),0.)
            cus(2,1,1,1) = cmplx(at1*real(cu(2,1,1,1)),0.)
            cus(3,1,1,1) = cmplx(at1*real(cu(3,1,1,1)),0.)
            cus(1,l1,1,1) = zero
            cus(2,l1,1,1) = zero
            cus(3,l1,1,1) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 100 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            cus(1,l,j,k1) = zero
            cus(2,l,j,k1) = zero
            cus(3,l,j,k1) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 110 l = 1, nz
            cus(1,l,1,k1) = zero
            cus(2,l,1,k1) = zero
            cus(3,l,1,k1) = zero
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
!$OMP PARALLEL DO PRIVATE(j,k,l,k1,l1)
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
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(j,k,l,k1,l1)
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
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(j,k,l,k1,l1)
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
!$OMP END PARALLEL DO
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
c-----------------------------------------------------------------------
      subroutine PPRDVMODES32(vpot,vpott,nx,ny,nz,modesx,modesy,modesz, 
     1ndim,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
c this subroutine extracts lowest order modes from packed complex vector
c array vpot and stores them into a location in an unpacked complex
c vector array vpott
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
c ndim = number of field arrays, must be >= 1
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c modeszd = second dimension of array vpott,
c where modeszd  = min(2*modesz-1,nz)
c modesypd = fourth dimension of array vpott, modesypd >= kyzp
c modesxpd = thid dimension of array vpott,
c modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
c in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, ndim, kstrt
      integer nvpy, nvpz, nzv, kxyp, kyzp, modesxpd, modesypd, modeszd
      complex vpot, vpott
      dimension vpot(ndim,nzv,kxyp,kyzp)
      dimension vpott(ndim,modeszd,modesxpd,modesypd)
c local data
      integer i, j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
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
!$OMP PARALLEL DO PRIVATE(i,j,k,l,k1,l1)
      do 110 k = 1, kyzps
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))
     1 then
         if (k1.gt.nyh) k1 = k1 - ny
         do 40 j = 1, jmax
         if ((j+joff).gt.1) then
            do 20 l = 2, lmax
            l1 = nz2 - l
            do 10 i = 1, ndim
            vpott(i,2*l-2,j,k) = vpot(i,l,j,k)
            vpott(i,2*l-1,j,k) = vpot(i,l1,j,k)
   10       continue
   20       continue
c mode numbers kz = 0, nz/2
            do 30 i = 1, ndim
            vpott(i,1,j,k) = vpot(i,1,j,k)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               vpott(i,nz,j,k) = vpot(i,l1,j,k)
            endif
   30       continue
         endif
   40    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c kx = 0
            if (k1.gt.0) then
               do 60 l = 2, lmax
               l1 = nz2 - l
               do 50 i = 1, ndim
               vpott(i,2*l-2,1,k) = vpot(i,l,1,k)
               vpott(i,2*l-1,1,k) = vpot(i,l1,1,k)
   50          continue
   60          continue
c mode numbers kz = 0, nz/2
               do 70 i = 1, ndim
               vpott(i,1,1,k) = vpot(i,1,1,k)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(i,nz,1,k) = vpot(i,l1,1,k)
               endif
   70          continue
c kx = nx/2
            else
               if (modesx.gt.nxh) then
                  do 90 l = 2, lmax
                  l1 = nz2 - l
                  do 80 i = 1, ndim
                  vpott(i,2*l-2,1,k) = vpot(i,l,1,k)
                  vpott(i,2*l-1,1,k) = vpot(i,l1,1,k)
   80             continue
   90             continue
c mode numbers kz = 0, nz/2
                  do 100 i = 1, ndim
                  vpott(i,1,1,k) = vpot(i,1,1,k)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpott(i,nz,1,k) = vpot(i,l1,1,k)
                  endif
  100             continue
               endif
            endif
         endif
      endif
  110 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
c ky = 0
      if (ks.eq.0) then
         do 150 j = 1, jmax
         if ((j+joff).gt.1) then
            do 130 l = 2, lmax
            l1 = nz2 - l
            do 120 i = 1, ndim
            vpott(i,2*l-2,j,1) = vpot(i,l,j,1)
            vpott(i,2*l-1,j,1) = vpot(i,l1,j,1)
  120       continue
  130       continue
c mode numbers kz = 0, nz/2
            do 140 i = 1, ndim
            vpott(i,1,j,1) = vpot(i,1,j,1)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               vpott(i,nz,j,1) = vpot(i,l1,j,1)
            endif
  140       continue
         endif
  150    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c kx = 0
            do 170 l = 2, lmax
            do 160 i = 1, ndim
            vpott(i,2*l-2,1,1) = vpot(i,l,1,1)
            vpott(i,2*l-1,1,1) = conjg(vpot(i,l,1,1))
  160       continue
  170       continue
c mode numbers kz = 0, nz/2
            do 180 i = 1, ndim
            vpott(i,1,1,1) = cmplx(real(vpot(i,1,1,1)),0.0)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               vpott(i,nz,1,1) = cmplx(real(vpot(i,l1,1,1)),0.0)
            endif
  180       continue
c kx = nx/2
            if (modesx.gt.nxh) then
               do 200 l = 2, lmax
               l1 = nz2 - l
               do 190 i = 1, ndim
               vpott(i,2*l-2,j1,1) = conjg(vpot(i,l1,1,1))
               vpott(i,2*l-1,j1,1) = vpot(i,l1,1,1)
  190          continue
  200          continue
c mode numbers kz = 0, nz/2
               do 210 i = 1, ndim
               vpott(i,1,j1,1) = cmplx(aimag(vpot(i,1,1,1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(i,nz,j1,1) = cmplx(aimag(vpot(i,l1,1,1)),
     10.)
               endif
  210          continue
            endif
         endif
      endif
c ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         if (modesy.gt.nyh) then
            do 250 j = 1, jmax
            if ((j+joff).gt.1) then
               do 230 l = 2, lmax
               l1 = nz2 - l
               do 220 i = 1, ndim
               vpott(i,2*l-2,j,k1) = vpot(i,l,j,k1)
               vpott(i,2*l-1,j,k1) = vpot(i,l1,j,k1)
  220          continue
  230          continue
c mode numbers kz = 0, nz/2
               do 240 i = 1, ndim
               vpott(i,1,j,k1) = vpot(i,1,j,k1)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(i,nz,j,k1) = vpot(i,l1,j,k1)
               endif
  240          continue
            endif
  250       continue
c mode numbers kx = 0, nx/2
            if (js.eq.0) then
c kx = 0
               do 270 l = 2, lmax
               do 260 i = 1, ndim
               vpott(i,2*l-2,1,k1) = vpot(i,l,1,k1)
               vpott(i,2*l-1,1,k1) = conjg(vpot(i,l,1,k1))
  260          continue
  270          continue
c mode numbers kz = 0, nz/2
               do 280 i = 1, ndim
               vpott(i,1,1,k1) = cmplx(real(vpot(i,1,1,k1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(i,nz,1,k1) = cmplx(real(vpot(i,l1,1,k1)),0.0)
               endif
  280          continue
c kx  = nx/2
               if (modesx.gt.nxh) then
                  do 300 l = 2, lmax
                  l1 = nz2 - l
                  do 290 i = 1, ndim
                  vpott(i,2*l-2,j1,k1) = conjg(vpot(i,l1,1,k1))
                  vpott(i,2*l-1,j1,k1) = vpot(i,l1,1,k1)
  290             continue
  300             continue
c mode numbers kz = 0, nz/2
                  do 310 i = 1, ndim
                  vpott(i,1,j1,k1) = cmplx(aimag(vpot(i,1,1,k1)),0.0)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpott(i,nz,j1,k1) = cmplx(aimag(vpot(i,l1,1,k1)),
     1                                         0.0)
                  endif
  310             continue
               endif
            endif
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPWRVMODES32(vpot,vpott,nx,ny,nz,modesx,modesy,modesz, 
     1ndim,kstrt,nvpy,nvpz,nzv,kxyp,kyzp,modesxpd,modesypd,modeszd)
c this subroutine extracts lowest order modes from a location in an
c unpacked complex vector array vpott and stores them into packed
c complex vector array vpot
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
c ndim = number of field arrays, must be >= 1
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c modeszd = second dimension of array vpott,
c where modeszd  = min(2*modesz-1,nz)
c modesypd = fourth dimension of array vpott, modesypd >= kyzp
c modesxpd = thid dimension of array vpott,
c modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
c in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, ndim, kstrt
      integer nvpy, nvpz, nzv, kxyp, kyzp, modesxpd, modesypd, modeszd
      complex vpot, vpott
      dimension vpot(ndim,nzv,kxyp,kyzp)
      dimension vpott(ndim,modeszd,modesxpd,modesypd)
c local data
      integer i, j, k, l, nxh, nyh, nzh, nz2, js, ks, joff, koff
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
!$OMP PARALLEL DO PRIVATE(i,j,k,l,k1,l1)
      do 200 k = 1, kyzps
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))
     1 then
         if (k1.gt.nyh) k1 = k1 - ny
         do 60 j = 1, jmax
         if ((j+joff).gt.1) then
            do 20 l = 2, lmax
            l1 = nz2 - l
            do 10 i = 1, ndim
            vpot(i,l,j,k) = vpott(i,2*l-2,j,k)
            vpot(i,l1,j,k) = vpott(i,2*l-1,j,k)
   10       continue
   20       continue
            do 40 l = lmax+1, nzh
            l1 = nz2 - l
            do 30 i = 1, ndim
            vpot(i,l,j,k) = zero
            vpot(i,l1,j,k) = zero
   30       continue
   40       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1 
            do 50 i = 1, ndim
            vpot(i,1,j,k) = vpott(i,1,j,k)
            vpot(i,l1,j,k) = zero
            if (modesz.gt.nzh) then
               vpot(i,l1,j,k) = vpott(i,nz,j,k)
            endif
   50       continue
         endif
   60    continue
         do 90 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 80 l = 1, nz
            do 70 i = 1, ndim
            vpot(i,l,j,k) = zero
   70       continue
   80       continue
         endif
   90    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c kx = 0
            if (k1.gt.0) then
               do 110 l = 2, lmax
               l1 = nz2 - l
               do 100 i = 1, ndim
               vpot(i,l,1,k) = vpott(i,2*l-2,1,k)
               vpot(i,l1,1,k) = vpott(i,2*l-1,1,k)
  100          continue
  110          continue
               do 130 l = lmax+1, nzh
               l1 = nz2 - l
               do 120 i = 1, ndim
               vpot(i,l,1,k) = zero
               vpot(i,l1,1,k) = zero
  120          continue
  130          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               do 140 i = 1, ndim
               vpot(i,1,1,k) = vpott(i,1,1,k)
               vpot(i,l1,1,k) = zero
               if (modesz.gt.nzh) then
                  vpot(i,l1,1,k) = vpott(i,nz,1,k)
               endif
  140          continue
c kx = nx/2
            else
               do 160 l = 1, nz
               do 150 i = 1, ndim
               vpot(i,l,1,k) = zero
  150          continue
  160          continue
               if (modesx.gt.nxh) then
                  do 180 l = 2, lmax
                  l1 = nz2 - l
                  do 170 i = 1, ndim
                  vpot(i,l,1,k) = vpott(i,2*l-2,1,k)
                  vpot(i,l1,1,k) = vpott(i,2*l-1,1,k)
  170             continue
  180             continue
c mode numbers kz = 0, nz/2
                  do 190 i = 1, ndim
                  vpot(i,1,1,k) = vpott(i,1,1,k)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpot(i,l1,1,k) = vpott(i,nz,1,k)
                  endif
  190             continue
               endif
            endif
         endif
      endif
  200 continue
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i,j,k,l,k1,l1)
      do 260 k = 1, kyzps
      k1 = k + koff
      if ((k1.ge.kmax).and.(k1.le.kmin).and.(k1.ne.nyh)) then
         do 230 j = 1, kxyps
         if ((j+joff).gt.1) then
            do 220 l = 1, nz
            do 210 i = 1, ndim
            vpot(i,l,j,k) = zero
  210       continue
  220       continue
         endif
  230    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 250 l = 1, nz
            do 240 i = 1, ndim
            vpot(i,l,1,k) = zero
  240       continue
  250       continue
         endif
      endif
  260 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
c ky = 0
      if (ks.eq.0) then
         do 320 j = 1, jmax
         if ((j+joff).gt.1) then
            do 280 l = 2, lmax
            l1 = nz2 - l
            do 270 i = 1, ndim
            vpot(i,l,j,1) = vpott(i,2*l-2,j,1)
            vpot(i,l1,j,1) = vpott(i,2*l-1,j,1)
  270       continue
  280       continue
            do 300 l = lmax+1, nzh
            l1 = nz2 - l
            do 290 i = 1, ndim
            vpot(i,l,j,1) = zero
            vpot(i,l1,j,1) = zero
  290       continue
  300       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            do 310 i = 1, ndim
            vpot(i,1,j,1) = vpott(i,1,j,1)
            vpot(i,l1,j,1) = zero
            if (modesz.gt.nzh) then
               vpot(i,l1,j,1) = vpott(i,nz,j,1)
            endif
  310       continue
         endif
  320    continue
         do 350 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 340 l = 1, nz
            do 330 i = 1, ndim
            vpot(i,l,j,1) = zero
  330       continue
  340       continue
         endif
  350    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
c kx = 0
            do 370 l = 2, lmax
            l1 = nz2 - l
            do 360 i = 1, ndim
            vpot(i,l,1,1) = vpott(i,2*l-2,1,1)
            vpot(i,l1,1,1) = zero
  360       continue
  370       continue
            do 390 l = lmax+1, nzh
            l1 = nz2 - l
            do 380 i = 1, ndim
            vpot(i,l,1,1) = zero
            vpot(i,l1,1,1) = zero
  380       continue
  390       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            do 400 i = 1, ndim
            vpot(i,1,1,1) = cmplx(real(vpott(i,1,1,1)),0.0)
            vpot(i,l1,1,1) = zero
            if (modesz.gt.nzh) then
               vpot(i,l1,1,1) = cmplx(real(vpott(i,nz,1,1)),0.0)
            endif
  400       continue
c kx = nx/2
            if (modesx.gt.nxh) then
               do 420 l = 2, lmax
               l1 = nz2 - l
               do 410 i = 1, ndim
               vpot(i,l1,1,1) = conjg(vpott(i,2*l-2,j1,1))
  410          continue
  420          continue
c mode numbers kz = 0, nz/2
               do 430 i = 1, ndim
               vpot(i,1,1,1) = cmplx(real(vpot(i,1,1,1)),
     1                               real(vpott(i,1,j1,1)))
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpot(i,l1,1,1) = cmplx(real(vpot(i,l1,1,1)),
     1                                   real(vpott(i,nz,j1,1)))
               endif
  430          continue
            endif
         endif
      endif
c ky = ny/2
      k1 = nyh/kyzp
      if (ks.eq.k1) then
         k1 = nyh - kyzp*k1 + 1
         do 460 j = 1, jmax
         if ((j+joff).gt.1) then
            do 450 l = 1, nz
            do 440 i = 1, ndim
            vpot(i,l,j,k1) = zero
  440       continue
  450       continue
         endif
  460    continue
         do 490 j = jmax+1, kxyps
         if ((j+joff).gt.1) then
            do 480 l = 1, nz
            do 470 i = 1, ndim
            vpot(i,l,j,k1) = zero
  470       continue
  480       continue
         endif
  490    continue
c mode numbers kx = 0, nx/2
         if (js.eq.0) then
            do 510 l = 1, nz
            do 500 i = 1, ndim
            vpot(i,l,1,k1) = zero
  500       continue
  510       continue
         endif
         if (modesy.gt.nyh) then
            do 550 j = 1, jmax
            if ((j+joff).gt.1) then
               do 530 l = 2, lmax
               l1 = nz2 - l
               do 520 i = 1, ndim
               vpot(i,l,j,k1) = vpott(i,2*l-2,j,k1)
               vpot(i,l1,j,k1) = vpott(i,2*l-1,j,k1)
  520          continue
  530          continue
c mode numbers kz = 0, nz/2
               do 540 i = 1, ndim
               vpot(i,1,j,k1) = vpott(i,1,j,k1)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpot(i,l1,j,k1) = vpott(i,nz,j,k1)
               endif
  540          continue
            endif
  550       continue
c mode numbers kx = 0, nx/2
            if (js.eq.0) then
c kx = 0
               do 570 l = 2, lmax
               do 560 i = 1, ndim
               vpot(i,l,1,k1) = vpott(i,2*l-2,1,k1)
  560          continue
  570          continue
c mode numbers kz = 0, nz/2
               do 580 i = 1, ndim
               vpot(i,1,1,k1) = cmplx(real(vpott(i,1,1,k1)),0.0)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpot(i,l1,1,k1) = cmplx(real(vpott(i,nz,1,k1)),0.0)
               endif
  580          continue
c kx  = nx/2
               if (modesx.gt.nxh) then
                  do 600 l = 2, lmax
                  l1 = nz2 - l
                  do 590 i = 1, ndim
                  vpot(i,l1,1,k1) = conjg(vpott(i,2*l-2,j1,k1))
  590             continue
  600             continue
c mode numbers kz = 0, nz/2
                  do 610 i = 1, ndim
                  vpot(i,1,1,k1) = cmplx(real(vpot(i,1,1,k1)),
     1                                   real(vpott(i,1,j1,k1)))

                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpot(i,l1,1,k1) = cmplx(real(vpot(i,l1,1,k1)),
     1                                       real(vpott(i,nz,j1,k1)))
                  endif
  610             continue
               endif
            endif
         endif
      endif
      return
      end
