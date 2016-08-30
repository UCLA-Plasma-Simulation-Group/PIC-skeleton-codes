c-----------------------------------------------------------------------
c Fortran library for Skeleton 3D Electrostatic OpenMP PIC Code field
c diagnostics
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c-----------------------------------------------------------------------
      subroutine MPOTP3(q,pot,ffc,we,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,   
     1nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c for potential with periodic boundary conditions.
c input: q,ffc,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd, output: pot,we
c approximate flop count is:
c 26*nxc*nyc*nzc + 14*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c potential is calculated using the equation:
c pot(kx,ky,kz) = g(kx,ky,kz)*q(kx,ky,kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c pot(kx=pi) = 0, pot(ky=pi) = 0, pot(kz=pi) = 0, and
c pot(kx=0,ky=0,kz=0) = 0
c q(j,k,l) = complex charge density for fourier mode (j-1,k-1,l-1)
c pot(j,k,l) = complex potential for fourier mode (j-1,k-1,l-1)
c real(ffc(j,k,l)) = potential green's function g
c imag(ffc(j,k,l)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1,l-1)
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c where affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c nx/ny/nz = system length in x/y/z direction
c nxvh = first dimension of scalar field arrays, must be >= nx/2
c nyv = second dimension of scalar field arrays, must be >= ny
c nzv = third dimension of scalar field arrays, must be >= nz
c nxhd = first dimension of form factor array, must be >= nx/2
c nyhd = second dimension of form factor array, must be >= ny/2
c nzhd = third dimension of form factor array, must be >= nz/2
      implicit none
      integer  nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
      real we
      complex q, pot, ffc
      dimension q(nxvh,nyv,nzv), pot(nxvh,nyv,nzv)
      dimension ffc(nxhd,nyhd,nzhd)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real at1, at2
      complex zero
      double precision wp, sum1, sum2
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
      zero = cmplx(0.0,0.0)
c calculate potential and sum field energy
      sum1 = 0.0d0
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,k1,l1,at1,at2,wp) REDUCTION(+:sum1)
      do 50 l = 2, nzh
      l1 = nz2 - l
      wp = 0.0d0
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nxh
      at2 = real(ffc(j,k,l))
      at1 = at2*aimag(ffc(j,k,l))
      pot(j,k,l) = at2*q(j,k,l)
      pot(j,k1,l) = at2*q(j,k1,l)
      pot(j,k,l1) = at2*q(j,k,l1)
      pot(j,k1,l1) = at2*q(j,k1,l1)
      wp = wp + at1*(q(j,k,l)*conjg(q(j,k,l))                           
     1   + q(j,k1,l)*conjg(q(j,k1,l)) + q(j,k,l1)*conjg(q(j,k,l1))      
     2   + q(j,k1,l1)*conjg(q(j,k1,l1)))
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      at2 = real(ffc(1,k,l))
      at1 = at2*aimag(ffc(1,k,l))
      pot(1,k,l) = at2*q(1,k,l)
      pot(1,k1,l) = zero
      pot(1,k,l1) = at2*q(1,k,l1)
      pot(1,k1,l1) = zero
      wp = wp + at1*(q(1,k,l)*conjg(q(1,k,l))                           
     1   + q(1,k,l1)*conjg(q(1,k,l1)))
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      at2 = real(ffc(j,1,l))
      at1 = at2*aimag(ffc(j,1,l))
      pot(j,1,l) = at2*q(j,1,l)
      pot(j,k1,l) = zero
      pot(j,1,l1) = at2*q(j,1,l1)
      pot(j,k1,l1) = zero
      wp = wp + at1*(q(j,1,l)*conjg(q(j,1,l))                           
     1   + q(j,1,l1)*conjg(q(j,1,l1)))
   40 continue
c mode numbers kx = 0, nx/2
      at2 = real(ffc(1,1,l))
      at1 = at2*aimag(ffc(1,1,l))
      pot(1,1,l) = at2*q(1,1,l)
      pot(1,k1,l) = zero
      pot(1,1,l1) = zero
      pot(1,k1,l1) = zero
      wp = wp + at1*(q(1,1,l)*conjg(q(1,1,l)))
      sum1 = sum1 + wp
   50 continue
!$OMP END DO NOWAIT
!$OMP END PARALLEL
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      sum2 = 0.0d0
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1,at2,wp) REDUCTION(+:sum2)
      do 70 k = 2, nyh
      k1 = ny2 - k
      wp = 0.0d0
      do 60 j = 2, nxh
      at2 = real(ffc(j,k,1))
      at1 = at2*aimag(ffc(j,k,1))
      pot(j,k,1) = at2*q(j,k,1)
      pot(j,k1,1) = at2*q(j,k1,1)
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
      wp = wp + at1*(q(j,k,1)*conjg(q(j,k,1))
     1   + q(j,k1,1)*conjg(q(j,k1,1)))
   60 continue
c mode numbers kx = 0, nx/2
      at2 = real(ffc(1,k,1))
      at1 = at2*aimag(ffc(1,k,1))
      pot(1,k,1) = at2*q(1,k,1)
      pot(1,k1,1) = zero
      pot(1,k,l1) = zero
      pot(1,k1,l1) = zero
      wp = wp + at1*(q(1,k,1)*conjg(q(1,k,1)))
      sum2 = sum2 + wp
   70 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 80 j = 2, nxh
      at2 = real(ffc(j,1,1))
      at1 = at2*aimag(ffc(j,1,1))
      pot(j,1,1) = at2*q(j,1,1)
      pot(j,k1,1) = zero
      pot(j,1,l1) = zero
      pot(j,k1,l1) = zero
      wp = wp + at1*(q(j,1,1)*conjg(q(j,1,1)))
   80 continue
      pot(1,1,1) = zero
      pot(1,k1,1) = zero
      pot(1,1,l1) = zero
      pot(1,k1,l1) = zero
      we = real(nx)*real(ny)*real(nz)*(sum1 + sum2 + wp)
      return
      end
c-----------------------------------------------------------------------
      subroutine MDIVF3(f,df,nx,ny,nz,nxvh,nyv,nzv)
c this subroutine calculates the divergence in fourier space
c input: all except df, output: df
c approximate flop count is:
c 35*nxc*nyc*nzc + 16*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c the divergence is calculated using the equation:
c df(kx,ky,kz) = sqrt(-1)*(kx*fx(kx,ky,kz)+ky*fy(kx,ky,kz)
c                       +kz*fz(kx,ky,kz))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c df(kx=pi) = 0, df(ky=pi) = 0, df(kz=pi) = 0
c and df(kx=0,ky=0,kz=0) = 0.
c nx/ny/nz = system length in x/y/z direction
c nxvh = first dimension of scalar field arrays, must be >= nx/2
c nyv = second dimension of scalar field arrays, must be >= ny
c nzv = third dimension of scalar field arrays, must be >= nz
      implicit none
      integer nx, ny, nz, nxvh, nyv, nzv
      complex f, df
      dimension f(3,nxvh,nyv,nzv), df(nxvh,nyv,nzv)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real dnx, dny, dnz, dkx, dky, dkz
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
c calculate the divergence
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,k1,l1,dkx,dky,dkz,zt1)
      do 50 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = dkx*f(1,j,k,l) + dky*f(2,j,k,l) + dkz*f(3,j,k,l)
      df(j,k,l) = cmplx(-aimag(zt1),real(zt1))
      zt1 = dkx*f(1,j,k1,l) - dky*f(2,j,k1,l) + dkz*f(3,j,k1,l)
      df(j,k1,l) = cmplx(-aimag(zt1),real(zt1))
      zt1 = dkx*f(1,j,k,l1) + dky*f(2,j,k,l1) - dkz*f(3,j,k,l1)
      df(j,k,l1) = cmplx(-aimag(zt1),real(zt1))
      zt1 = dkx*f(1,j,k1,l1) - dky*f(2,j,k1,l1) - dkz*f(3,j,k1,l1)
      df(j,k1,l1) = cmplx(-aimag(zt1),real(zt1))
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      zt1 = dky*f(2,1,k,l) + dkz*f(3,1,k,l)
      df(1,k,l) = cmplx(-aimag(zt1),real(zt1))
      df(1,k1,l) = zero
      zt1 = dky*f(2,1,k,l1) - dkz*f(3,1,k,l1)
      df(1,k,l1) = cmplx(-aimag(zt1),real(zt1))
      df(1,k1,l1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = dkx*f(1,j,1,l) + dkz*f(3,j,1,l)
      df(j,1,l) = cmplx(-aimag(zt1),real(zt1))
      df(j,k1,l) = zero
      zt1 = dkx*f(1,j,1,l1) - dkz*f(3,j,1,l1)
      df(j,1,l1) = cmplx(-aimag(zt1),real(zt1))
      df(j,k1,l1) = zero
   40 continue
c mode numbers kx = 0, nx/2
      df(1,1,l) = dkz*cmplx(-aimag(f(3,1,1,l)),real(f(3,1,1,l)))
      df(1,k1,l) = zero
      df(1,1,l1) = zero
      df(1,k1,l1) = zero
   50 continue
!$OMP END DO NOWAIT
!$OMP END PARALLEL
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,zt1)
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 60 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = dkx*f(1,j,k,1) + dky*f(2,j,k,1)
      df(j,k,1) = cmplx(-aimag(zt1),real(zt1))
      zt1 = dkx*f(1,j,k1,1) - dky*f(2,j,k1,1)
      df(j,k1,1) = cmplx(-aimag(zt1),real(zt1))
      df(j,k,l1) = zero
      df(j,k1,l1) = zero
   60 continue
c mode numbers kx = 0, nx/2
      df(1,k,1) = dky*cmplx(-aimag(f(2,1,k,1)),real(f(2,1,k,1)))
      df(1,k1,1) = zero
      df(1,k,l1) = zero
      df(1,k1,l1) = zero
   70 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 80 j = 2, nxh
      dkx = dnx*real(j - 1)
      df(j,1,1) = dkx*cmplx(-aimag(f(1,j,1,1)),real(f(1,j,1,1)))
      df(j,k1,1) = zero
      df(j,1,l1) = zero
      df(j,k1,l1) = zero
   80 continue
      df(1,1,1) = zero
      df(1,k1,1) = zero
      df(1,1,l1) = zero
      df(1,k1,l1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRADF3(df,f,nx,ny,nz,nxvh,nyv,nzv)
c this subroutine calculates the gradient in fourier space
c input: all except f, output: f
c approximate flop count is:
c 30*nxc*nyc*nzc + 12*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
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
c nxvh = first dimension of scalar field arrays, must be >= nxh
c nyv = second dimension of scalar field arrays, must be >= ny
c nzv = third dimension of scalar field arrays, must be >= nz
      implicit none
      integer nx, ny, nz, nxvh, nyv, nzv
      complex df, f
      dimension df(nxvh,nyv,nzv), f(3,nxvh,nyv,nzv)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real dnx, dny, dnz, dkx, dky, dkz
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
c calculate the gradient
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,k1,l1,dkx,dky,dkz,zt1)
      do 50 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = cmplx(-aimag(df(j,k,l)),real(df(j,k,l)))
      f(1,j,k,l) = dkx*zt1
      f(2,j,k,l) = dky*zt1
      f(3,j,k,l) = dkz*zt1
      zt1 = cmplx(-aimag(df(j,k1,l)),real(df(j,k1,l)))
      f(1,j,k1,l) = dkx*zt1
      f(2,j,k1,l) = -dky*zt1
      f(3,j,k1,l) = dkz*zt1
      zt1 = cmplx(-aimag(df(j,k,l1)),real(df(j,k,l1)))
      f(1,j,k,l1) = dkx*zt1
      f(2,j,k,l1) = dky*zt1
      f(3,j,k,l1) = -dkz*zt1
      zt1 = cmplx(-aimag(df(j,k1,l1)),real(df(j,k1,l1)))
      f(1,j,k1,l1) = dkx*zt1
      f(2,j,k1,l1) = -dky*zt1
      f(3,j,k1,l1) = -dkz*zt1
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      zt1 = cmplx(-aimag(df(1,k,l)),real(df(1,k,l)))
      f(1,1,k,l) = 0.
      f(2,1,k,l) = dky*zt1
      f(3,1,k,l) = dkz*zt1
      f(1,1,k1,l) = zero
      f(2,1,k1,l) = zero
      f(3,1,k1,l) = zero
      zt1 = cmplx(-aimag(df(1,k,l1)),real(df(1,k,l1)))
      f(1,1,k,l1) = 0.
      f(2,1,k,l1) = dky*zt1
      f(3,1,k,l1) = -dkz*zt1
      f(1,1,k1,l1) = zero
      f(2,1,k1,l1) = zero
      f(3,1,k1,l1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = cmplx(-aimag(df(j,1,l)),real(df(j,1,l)))
      f(1,j,1,l) = dkx*zt1
      f(2,j,1,l) = zero
      f(3,j,1,l) = dkz*zt1
      f(1,j,k1,l) = zero
      f(2,j,k1,l) = zero
      f(3,j,k1,l) = zero
      zt1 = cmplx(-aimag(df(j,1,l1)),real(df(j,1,l1)))
      f(1,j,1,l1) = dkx*zt1
      f(2,j,1,l1) = zero
      f(3,j,1,l1) = -dkz*zt1
      f(1,j,k1,l1) = zero
      f(2,j,k1,l1) = zero
      f(3,j,k1,l1) = zero
   40 continue
c mode numbers kx = 0, nx/2
      f(1,1,1,l) = zero
      f(2,1,1,l) = zero
      f(3,1,1,l) = dkz*cmplx(-aimag(df(1,1,l)),real(df(1,1,l)))
      f(1,1,k1,l) = zero
      f(2,1,k1,l) = zero
      f(3,1,k1,l) = zero
      f(1,1,1,l1) = zero
      f(2,1,1,l1) = zero
      f(3,1,1,l1) = zero
      f(1,1,k1,l1) = zero
      f(2,1,k1,l1) = zero
      f(3,1,k1,l1) = zero
   50 continue
!$OMP END DO NOWAIT
!$OMP END PARALLEL
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,zt1)
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 60 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt1 = cmplx(-aimag(df(j,k,1)),real(df(j,k,1)))
      f(1,j,k,1) = dkx*zt1
      f(2,j,k,1) = dky*zt1
      f(3,j,k,1) = zero
      zt1 = cmplx(-aimag(df(j,k1,1)),real(df(j,k1,1)))
      f(1,j,k1,1) = dkx*zt1
      f(2,j,k1,1) = -dky*zt1
      f(3,j,k1,1) = zero
      f(1,j,k,l1) = zero
      f(2,j,k,l1) = zero
      f(3,j,k,l1) = zero
      f(1,j,k1,l1) = zero
      f(2,j,k1,l1) = zero
      f(3,j,k1,l1) = zero
   60 continue
c mode numbers kx = 0, nx/2
      f(1,1,k,1) = zero
      f(2,1,k,1) = dky*cmplx(-aimag(df(1,k,1)),real(df(1,k,1)))
      f(3,1,k,1) = zero
      f(1,1,k1,1) = zero
      f(2,1,k1,1) = zero
      f(3,1,k1,1) = zero
      f(1,1,k,l1) = zero
      f(2,1,k,l1) = zero
      f(3,1,k,l1) = zero
      f(1,1,k1,l1) = zero
      f(2,1,k1,l1) = zero
      f(3,1,k1,l1) = zero
   70 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 80 j = 2, nxh
      dkx = dnx*real(j - 1)
      f(1,j,1,1) = dkx*cmplx(-aimag(df(j,1,1)),real(df(j,1,1)))
      f(2,j,1,1) = zero
      f(3,j,1,1) = zero
      f(1,j,k1,1) = zero
      f(2,j,k1,1) = zero
      f(3,j,k1,1) = zero
      f(1,j,1,l1) = zero
      f(2,j,1,l1) = zero
      f(3,j,1,l1) = zero
      f(1,j,k1,l1) = zero
      f(2,j,k1,l1) = zero
      f(3,j,k1,l1) = zero
   80 continue
      f(1,1,1,1) = zero
      f(2,1,1,1) = zero
      f(3,1,1,1) = zero
      f(1,1,k1,1) = zero
      f(2,1,k1,1) = zero
      f(3,1,k1,1) = zero
      f(1,1,1,l1) = zero
      f(2,1,1,l1) = zero
      f(3,1,1,l1) = zero
      f(1,1,k1,l1) = zero
      f(2,1,k1,l1) = zero
      f(3,1,k1,l1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine MSMOOTH3(q,qs,ffc,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd)
c this subroutine provides a 3d scalar smoothing function
c in fourier space, with periodic boundary conditions.
c input: q,ffc,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd, output: qs
c approximate flop count is:
c 8*nxc*nyc*nzc + 4*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c smoothing is calculated using the equation:
c qs(kx,ky,kz) = q(kx,ky,kz)*s(kx,ky,kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c where affp = normalization constant = nx*ny*nz/np,
c and np=number of particles
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c qs(kx=pi) = qs(ky=pi) = qs(kz=pi) = 0, and qs(kx=0,ky=0,kz=0) = 0.
c q(j,k,l) = complex charge density
c qs(j,k,l) = complex smoothed charge density
c for fourier mode (j-1,k-1,l-1)
c real(ffc(j,k,l)) = potential green's function g
c aimag(ffc(j,k,l)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1,l-1)
c nx/ny/nz = system length in x/y/z direction
c nxvh = first dimension of scalar field arrays, must be >= nx/2
c nyv = second dimension of field arrays, must be >= ny
c nzv = third dimension of field arrays, must be >= nz
c nxhd = first dimension of form factor array, must be >= nx/2
c nyhd = second dimension of form factor array, must be >= ny/2
c nzhd = third dimension of form factor array, must be >= nz/2
      implicit none
      integer  nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
      complex q, qs, ffc
      dimension q(nxvh,nyv,nzv), qs(nxvh,nyv,nzv)
      dimension ffc(nxhd,nyhd,nzhd)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real at1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
      zero = cmplx(0.0,0.0)
c calculate smoothing
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,k1,l1,at1)
      do 50 l = 2, nzh
      l1 = nz2 - l
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nxh
      at1 = aimag(ffc(j,k,l))
      qs(j,k,l) = at1*q(j,k,l)
      qs(j,k1,l) = at1*q(j,k1,l)
      qs(j,k,l1) = at1*q(j,k,l1)
      qs(j,k1,l1) = at1*q(j,k1,l1)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      at1 = aimag(ffc(1,k,l))
      qs(1,k,l) = at1*q(1,k,l)
      qs(1,k1,l) = zero
      qs(1,k,l1) = at1*q(1,k,l1)
      qs(1,k1,l1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      at1 = aimag(ffc(j,1,l))
      qs(j,1,l) = at1*q(j,1,l)
      qs(j,k1,l) = zero
      qs(j,1,l1) = at1*q(j,1,l1)
      qs(j,k1,l1) = zero
   40 continue
c mode numbers kx = 0, nx/2
      at1 = aimag(ffc(1,1,l))
      qs(1,1,l) = at1*q(1,1,l)
      qs(1,k1,l) = zero
      qs(1,1,l1) = zero
      qs(1,k1,l1) = zero
   50 continue
!$OMP END DO NOWAIT
!$OMP END PARALLEL
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
!$OMP PARALLEL DO PRIVATE(j,k,k1,at1)
      do 70 k = 2, nyh
      k1 = ny2 - k
      do 60 j = 2, nxh
      at1 = aimag(ffc(j,k,1))
      qs(j,k,1) = at1*q(j,k,1)
      qs(j,k1,1) = at1*q(j,k1,1)
      qs(j,k,l1) = zero
      qs(j,k1,l1) = zero
   60 continue
c mode numbers kx = 0, nx/2
      at1 = aimag(ffc(1,k,1))
      qs(1,k,1) = at1*q(1,k,1)
      qs(1,k1,1) = zero
      qs(1,k,l1) = zero
      qs(1,k1,l1) = zero
   70 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 80 j = 2, nxh
      at1 = aimag(ffc(j,1,1))
      qs(j,1,1) = at1*q(j,1,1)
      qs(j,k1,1) = zero
      qs(j,1,l1) = zero
      qs(j,k1,l1) = zero
   80 continue
      qs(1,1,1) = cmplx(aimag(ffc(1,1,1))*real(q(1,1,1)),0.0)
      qs(1,k1,1) = zero
      qs(1,1,l1) = zero
      qs(1,k1,l1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine RDMODES3(pot,pott,nx,ny,nz,modesx,modesy,modesz,nxvh,  
     1nyv,nzv,modesxd,modesyd,modeszd)
c this subroutine extracts lowest order modes from packed complex array
c pot and stores them into a location in an unpacked complex array pott
c modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2),
c kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
c nx/ny/nz = system length in x/y/z direction
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
c nxvh = first dimension of input array pot, nxvh >= nx/2
c nyv = second dimension of input array pot, nyv >= ny
c nzv = third dimension of input array pot, nzv >= nz
c modesxd = first dimension of output array pott, modesxd >= modesx
c modesyd = second dimension of output array pott,
c where modesyd  >= min(2*modesy-1,ny)
c modeszd = third dimension of output array pott,
c where modeszd  >= min(2*modesz-1,nz)
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, nxvh, nyv, nzv
      integer modesxd, modesyd, modeszd
      complex pot, pott
      dimension pot(nxvh,nyv,nzv), pott(modesxd,modesyd,modeszd)
c local data
      integer nxh, nyh, nzh, lmax, kmax, jmax, ny2, nz2, j, k, l
      integer j1, k1, l1
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      ny2 = ny + 2
      nz2 = nz + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      lmax = min0(modesz,nzh)
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL DO PRIVATE(j,k,l,k1,l1)
      do 50 l = 2, lmax
      l1 = nz2 - l
      do 20 k = 2, kmax
      k1 = ny2 - k
      do 10 j = 2, jmax
      pott(j,2*k-2,2*l-2) = pot(j,k,l)
      pott(j,2*k-1,2*l-2) = pot(j,k1,l)
      pott(j,2*k-2,2*l-1) = pot(j,k,l1)
      pott(j,2*k-1,2*l-1) = pot(j,k1,l1)
   10 continue
c mode numbers kx = 0, nx/2
      pott(1,2*k-2,2*l-2) = pot(1,k,l)
      pott(1,2*k-1,2*l-2) = conjg(pot(1,k,l1))
      pott(1,2*k-2,2*l-1) = pot(1,k,l1)
      pott(1,2*k-1,2*l-1) = conjg(pot(1,k,l))
      if (modesx.gt.nxh) then
         pott(j1,2*k-2,2*l-2) = conjg(pot(1,k1,l1))
         pott(j1,2*k-1,2*l-2) = pot(1,k1,l)
         pott(j1,2*k-2,2*l-1) = conjg(pot(1,k1,l))
         pott(j1,2*k-1,2*l-1) = pot(1,k1,l1)
      endif
   20 continue
c mode numbers ky = 0
      do 30 j = 2, jmax
      pott(j,1,2*l-2) = pot(j,1,l)
      pott(j,1,2*l-1) = pot(j,1,l1)
   30 continue
c mode numbers kx = 0, nx/2
      pott(1,1,2*l-2) = pot(1,1,l)
      pott(1,1,2*l-1) = conjg(pot(1,1,l))
      if (modesx.gt.nxh) then
         pott(j1,1,2*l-2) = conjg(pot(1,1,l1))
         pott(j1,1,2*l-1) = pot(1,1,l1)
      endif
c mode numbers ky = ny/2
      if (modesy.gt.nyh) then
         k1 = nyh + 1
         do 40 j = 2, jmax
         pott(j,ny,2*l-2) = pot(j,k1,l)
         pott(j,ny,2*l-1) = pot(j,k1,l1)
   40    continue
         pott(1,ny,2*l-2) = pot(1,k1,l)
         pott(1,ny,2*l-1) = conjg(pot(1,k1,l))
         if (modesx.gt.nxh) then
            pott(j1,ny,2*l-2) = conjg(pot(1,k1,l1))
            pott(j1,ny,2*l-1) = pot(1,k1,l1)
         endif
      endif
   50 continue
!$OMP END PARALLEL DO
c mode numbers kz = 0
      do 70 k = 2, kmax
      k1 = ny2 - k
      do 60 j = 2, jmax
      pott(j,2*k-2,1) = pot(j,k,1)
      pott(j,2*k-1,1) = pot(j,k1,1)
   60 continue
c mode numbers kx = 0, nx/2
      pott(1,2*k-2,1) = pot(1,k,1)
      pott(1,2*k-1,1) = conjg(pot(1,k,1))
      if (modesx.gt.nxh) then
         pott(j1,2*k-2,1) = conjg(pot(1,k1,1))
         pott(j1,2*k-1,1) = pot(1,k1,1)
      endif
   70 continue
c mode numbers ky = 0
      do 80 j = 2, jmax
      pott(j,1,1) = pot(j,1,1)
   80 continue
      pott(1,1,1) = cmplx(real(pot(1,1,1)),0.0)
      if (modesx.gt.nxh) then
         pott(j1,1,1) = cmplx(aimag(pot(1,1,1)),0.0)
      endif
      if (modesy.gt.nyh) then
         k1 = nyh + 1
c mode numbers ky = ny/2
         do 90 j = 2, jmax
         pott(j,ny,1) = pot(j,k1,1)
   90    continue
         pott(1,ny,1) = cmplx(real(pot(1,k1,1)),0.0)
         if (modesx.gt.nxh) then
            pott(j1,ny,1) = cmplx(aimag(pot(1,k1,1)),0.0)
         endif
      endif
c mode numbers kz = nz/2
      if (modesz.gt.nzh) then
         l1 = nzh + 1
         do 110 k = 2, kmax
         k1 = ny2 - k
         do 100 j = 2, jmax
         pott(j,2*k-2,nz) = pot(j,k,l1)
         pott(j,2*k-1,nz) = pot(j,k1,l1)
  100    continue
c mode numbers kx = 0, nx/2
         pott(1,2*k-2,nz) = pot(1,k,l1)
         pott(1,2*k-1,nz) = conjg(pot(1,k,l1))
         if (modesx.gt.nxh) then
            pott(j1,2*k-2,nz) = conjg(pot(1,k1,l1))
            pott(j1,2*k-1,nz) = pot(1,k1,l1)
         endif
  110    continue
c mode numbers ky = 0
         do 120 j = 2, jmax
         pott(j,1,nz) = pot(j,1,l1)
  120    continue
         pott(1,1,nz) = cmplx(real(pot(1,1,l1)),0.0)
         if (modesx.gt.nxh) then
            pott(j1,1,nz) = cmplx(aimag(pot(1,1,l1)),0.0)
         endif
         if (modesy.gt.nyh) then
            k1 = nyh + 1
c mode numbers ky = ny/2
            do 130 j = 2, jmax
            pott(j,ny,nz) = pot(j,k1,l1)
  130       continue
            pott(1,ny,nz) = cmplx(real(pot(1,k1,l1)),0.0)
            if (modesx.gt.nxh) then
               pott(j1,ny,nz) = cmplx(aimag(pot(1,k1,l1)),0.0)
            endif
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WRMODES3(pot,pott,nx,ny,nz,modesx,modesy,modesz,nxvh,  
     1nyv,nzv,modesxd,modesyd,modeszd)
c this subroutine extracts lowest order modes from a location in an
c unpacked complex array pott and stores them into a packed complex
c array pot
c modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2),
c kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
c nx/ny/nz = system length in x/y/z direction
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
c nxvh = first dimension of output array pot, nxvh >= nx/2
c nyv = second dimension of output array pot, nyv >= ny
c nzv = third dimension of output array pot, nzv >= nz
c modesxd = first dimension of input array pott, modesxd >= modesx
c modesyd = second dimension of output array pott,
c where modesyd  >= min(2*modesy-1,ny)
c modeszd = third dimension of output array pott,
c where modeszd  >= min(2*modesz-1,nz)
      implicit none
      integer nx, ny, nz, modesx, modesy, modesz, nxvh, nyv, nzv
      integer modesxd, modesyd, modeszd
      complex pot, pott
      dimension pot(nxvh,nyv,nzv), pott(modesxd,modesyd,modeszd)
c local data
      integer nxh, nyh, nzh, lmax, kmax, jmax, ny2, nz2, j, k, l
      integer j1, k1, l1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      ny2 = ny + 2
      nz2 = nz + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      lmax = min0(modesz,nzh)
      j1 = nxh + 1
      zero = cmplx(0.0,0.0)
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL DO PRIVATE(j,k,l,k1,l1)
      do 90 l = 2, lmax
      l1 = nz2 - l
      do 30 k = 2, kmax
      k1 = ny2 - k
      do 10 j = 2, jmax
      pot(j,k,l) = pott(j,2*k-2,2*l-2)
      pot(j,k1,l) = pott(j,2*k-1,2*l-2)
      pot(j,k,l1) = pott(j,2*k-2,2*l-1)
      pot(j,k1,l1) = pott(j,2*k-1,2*l-1)
   10 continue
      do 20 j = jmax+1, nxh
      pot(j,k,l) = zero
      pot(j,k1,l) = zero
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
   20 continue
c mode numbers kx = 0, nx/2
      pot(1,k,l) = pott(1,2*k-2,2*l-2)
      pot(1,k1,l) = zero
      pot(1,k,l1) = pott(1,2*k-2,2*l-1)
      pot(1,k1,l1) = zero
      if (modesx.gt.nxh) then
         pot(1,k1,l) = conjg(pott(j1,2*k-2,2*l-1))
         pot(1,k1,l1) = conjg(pott(j1,2*k-2,2*l-2))
      endif
   30 continue
      do 50 k = kmax+1, nyh
      k1 = ny2 - k
      do 40 j = 1, nxh
      pot(j,k,l) = zero
      pot(j,k1,l) = zero
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
   40 continue
   50 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, jmax
      pot(j,1,l) = pott(j,1,2*l-2)
      pot(j,k1,l) = zero
      pot(j,1,l1) = pott(j,1,2*l-1)
      pot(j,k1,l1) = zero
   60 continue
c mode numbers kx = 0, nx/2
      pot(1,1,l) = pott(1,1,2*l-2)
      pot(1,k1,l) = zero
      pot(1,1,l1) = zero
      pot(1,k1,l1) = zero
      do 70 j = jmax+1, nxh
      pot(j,1,l) = zero
      pot(j,k1,l) = zero
      pot(j,1,l1) = zero
      pot(j,k1,l1) = zero
   70 continue
      if (modesx.gt.nxh) then
         pot(1,1,l1) = conjg(pott(j1,1,2*l-2))
      endif
c mode numbers ky = ny/2
      if (modesy.gt.nyh) then
         do 80 j = 2, jmax
         pot(j,k1,l) = pott(j,ny,2*l-2)
         pot(j,k1,l1) = pott(j,ny,2*l-1)
   80    continue
         pot(1,k1,l) = pott(1,ny,2*l-2)
         if (modesx.gt.nxh) then
            pot(1,k1,l1) = conjg(pott(j1,ny,2*l-2))
         endif
      endif
   90 continue
!$OMP END PARALLEL DO
      do 130 l = modesz+1, nzh
      l1 = nz2 - l
      do 110 k = 2, nyh
      k1 = ny2 - k
      do 100 j = 1, nxh
      pot(j,k,l) = zero
      pot(j,k1,l) = zero
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
  100 continue
  110 continue
      k1 = nyh + 1
      do 120 j = 1, nxh
      pot(j,1,l) = zero
      pot(j,k1,l) = zero
      pot(j,1,l1) = zero
      pot(j,k1,l1) = zero
  120 continue
  130 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 160 k = 2, kmax
      k1 = ny2 - k
      do 140 j = 2, jmax
      pot(j,k,1) = pott(j,2*k-2,1)
      pot(j,k1,1) = pott(j,2*k-1,1)
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
  140 continue
      do 150 j = jmax+1, nxh
      pot(j,k,1) = zero
      pot(j,k1,1) = zero
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
  150 continue
c mode numbers kx = 0, nx/2
      pot(1,k,1) = pott(1,2*k-2,1)
      pot(1,k1,1) = zero
      pot(1,k,l1) = zero
      pot(1,k1,l1) = zero
      if (modesx.gt.nxh) then
         pot(1,k1,1) = conjg(pott(j1,2*k-2,1))
      endif
  160 continue
      do 180 k = modesy+1, nyh
      k1 = ny2 - k
      do 170 j = 1, nxh
      pot(j,k,1) = zero
      pot(j,k1,1) = zero
      pot(j,k,l1) = zero
      pot(j,k1,l1) = zero
  170 continue
  180 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 190 j = 2, jmax
      pot(j,1,1) = pott(j,1,1)
      pot(j,k1,1) = zero
      pot(j,1,l1) = zero
      pot(j,k1,l1) = zero
  190 continue
      do 200 j = jmax+1, nxh
      pot(j,1,1) = zero
      pot(j,k1,1) = zero
      pot(j,1,l1) = zero
      pot(j,k1,l1) = zero
  200 continue
      pot(1,1,1) = cmplx(real(pott(1,1,1)),0.0)
      pot(1,k1,1) = zero
      pot(1,1,l1) = zero
      pot(1,k1,l1) = zero
      if (modesx.gt.nxh) then
         pot(1,1,1) = cmplx(real(pot(1,1,1)),real(pott(j1,1,1)))
      endif
      if (modesy.gt.nyh) then
         do 210 j = 2, jmax
         pot(j,k1,1) = pott(j,ny,1)
  210    continue
         pot(1,k1,1) = cmplx(real(pott(1,ny,1)),0.0)
         if (modesx.gt.nxh) then
            pot(1,k1,1) = cmplx(real(pot(1,k1,1)),real(pott(j1,ny,1)))
         endif
      endif
c mode numbers kz = nz/2
      if (modesz.gt.nzh) then
         do 230 k = 2, kmax
         k1 = ny2 - k
         do 220 j = 2, jmax
         pot(j,k,l1) = pott(j,2*k-2,nz)
         pot(j,k1,l1) = pott(j,2*k-1,nz)
  220    continue
c mode numbers kx = 0, nx/2
         pot(1,k,l1) = pott(1,2*k-2,nz)
         if (modesx.gt.nxh) then
            pot(1,k1,l1) = conjg(pott(j1,2*k-2,nz))
         endif
  230 continue
c mode numbers ky = 0
         do 240 j = 2, jmax
         pot(j,1,l1) = pott(j,1,nz)
  240    continue
         pot(1,1,l1) = cmplx(real(pott(1,1,nz)),0.0)
         if (modesx.gt.nxh) then
            pot(1,1,l1) = cmplx(real(pot(1,1,l1)),real(pott(j1,1,nz)))
         endif
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            do 250 j = 2, jmax
            pot(j,k1,l1) = pott(j,ny,nz)
  250       continue
            pot(1,k1,l1) = cmplx(real(pott(1,ny,nz)),0.0)
            if (modesx.gt.nxh) then
               pot(1,k1,l1) = cmplx(real(pot(1,k1,l1)),
     1                              real(pott(j1,ny,nz)))
            endif
         endif
      endif
      return
      end
