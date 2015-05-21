c Fortran Library for Skeleton 2D Electrostatic MPI PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine PDICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,idps)
c this subroutine determines spatial boundaries for uniform particle
c decomposition, calculates number of grid points in each spatial
c region, and the offset of these grid points from the global address
c nvp must be < ny.  some combinations of ny and nvp result in a zero
c value of nyp.  this is not supported.
c integer boundaries are set.
c input: ny, kstrt, nvp, idps, output: edges, nyp, noff, nypmx, nypmn
c edges(1) = lower boundary of particle partition
c edges(2) = upper boundary of particle partition
c nyp = number of primary (complete) gridpoints in particle partition
c noff = lowermost global gridpoint in particle partition
c nypmx = maximum size of particle partition, including guard cells
c nypmn = minimum value of nyp
c ny = system length in y direction
c kstrt = starting data block number (processor id + 1)
c nvp = number of real or virtual processors
c idps = number of partition boundaries
      implicit none
      integer nyp, noff, nypmx, nypmn, ny, kstrt, nvp, idps
      real edges
      dimension edges(idps)
c local data
      integer kb, kyp
      real at1, any
      integer mypm, iwork2
      dimension mypm(2), iwork2(2)
      any = real(ny)
c determine decomposition
      kb = kstrt - 1
      kyp = (ny - 1)/nvp + 1
      at1 = real(kyp)
      edges(1) = at1*real(kb)
      if (edges(1).gt.any) edges(1) = any
      noff = edges(1)
      edges(2) = at1*real(kb + 1)
      if (edges(2).gt.any) edges(2) = any
      kb = edges(2)
      nyp = kb - noff
c find maximum/minimum partition size
      mypm(1) = nyp
      mypm(2) = -nyp
      call PPIMAX(mypm,iwork2,2)
      nypmx = mypm(1) + 1
      nypmn = -mypm(2)
      return
      end
c-----------------------------------------------------------------------
      subroutine PDISTR2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx, 
     1ny,idimp,npmax,idps,ipbc,ierr)
c for 2d code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c for distributed data.
c input: all except part, npp, ierr, output: part, npp, ierr
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = velocity vx of particle n in partition
c part(4,n) = velocity vy of particle n in partition
c edges(1) = lower boundary of particle partition
c edges(2) = upper boundary of particle partition
c npp = number of particles in partition
c nps = starting address of particles in partition
c vtx/vty = thermal velocity of electrons in x/y direction
c vdx/vdy = drift velocity of beam electrons in x/y direction
c npx/npy = initial number of particles distributed in x/y direction
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c idps = number of partition boundaries
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ierr = (0,1) = (no,yes) error condition exists
c ranorm = gaussian random number with zero mean and unit variance
c with spatial decomposition
      implicit none
      integer npp, nps, npx, npy, nx, ny, idimp, npmax, idps, ipbc, ierr
      real vtx, vty, vdx, vdy
      real part, edges
      dimension part(idimp,npmax), edges(idps)
c local data
      integer j, k, npt, npxyp
      real edgelx, edgely, at1, at2, xt, yt, vxt, vyt
      double precision dnpx, dnpxy, dt1
      integer ierr1, iwork1
      double precision sum3, work3
      dimension ierr1(1), iwork1(1), sum3(3), work3(3)
      double precision ranorm
      ierr = 0
c particle distribution constant
      dnpx = dble(npx)
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
      yt = edgely + at2*(real(k) - 0.5)
      do 10 j = 1, npx
      xt = edgelx + at1*(real(j) - 0.5)
c maxwellian velocity distribution
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      if ((yt.ge.edges(1)).and.(yt.lt.edges(2))) then
         npt = npp + 1
         if (npt.le.npmax) then
            part(1,npt) = xt
            part(2,npt) = yt
            part(3,npt) = vxt
            part(4,npt) = vyt
            npp = npt
         else
            ierr = ierr + 1
         endif
      endif
   10 continue
   20 continue
      npxyp = 0
c add correct drift
      sum3(1) = 0.0d0
      sum3(2) = 0.0d0
      do 30 j = nps, npp
      npxyp = npxyp + 1
      sum3(1) = sum3(1) + part(3,j)
      sum3(2) = sum3(2) + part(4,j)
   30 continue
      sum3(3) = npxyp
      call PPDSUM(sum3,work3,3)
      dnpxy = sum3(3)
      ierr1(1) = ierr
      call PPIMAX(ierr1,iwork1,1)
      ierr = ierr1(1)
      dt1 = 1.0d0/dnpxy
      sum3(1) = dt1*sum3(1) - vdx
      sum3(2) = dt1*sum3(2) - vdy
      do 40 j = nps, npp
      part(3,j) = part(3,j) - sum3(1)
      part(4,j) = part(4,j) - sum3(2)
   40 continue
c process errors
      dnpxy = dnpxy - dnpx*dble(npy)
      if (dnpxy.ne.0.0d0) ierr = dnpxy
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGPUSH2L(part,fxy,edges,npp,noff,ihole,qbm,dt,ek,nx,ny
     1,idimp,npmax,nxv,nypmx,idps,ntmax,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions
c also determines list of particles which are leaving this processor
c scalar version using guard cells, for distributed data
c 42 flops/particle, 12 loads, 4 stores
c input: all except ihole, output: part, ihole, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = velocity vx of particle n in partition
c part(4,n) = velocity vy of particle n in partition
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c in other words, fxy are the convolutions of the electric field
c over the particle shape, where kk = k + noff - 1
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nxv = first dimension of field array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c idps = number of partition boundaries
c ntmax = size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nypmx
      integer ipbc
      real qbm, dt, ek
      real part, fxy, edges
      integer ihole
      dimension part(idimp,npmax), fxy(2,nxv,nypmx)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, j, nn, mm, np, mp, ih, nh
      real qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy
      double precision sum1
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      mnoff = noff - 1
      ih = 0
      nh = 0
      do 10 j = 1, npp
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
c find acceleration
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp))                    
     1   + amy*(dxp*fxy(1,np,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp))                    
     1   + amy*(dxp*fxy(2,np,mm) + amx*fxy(2,nn,mm))
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c set end of file flag
      if (nh.gt.0) ih = -ih
      ihole(1) = ih
c normalize kinetic energy
      ek = ek + 0.125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGPOST2L(part,q,npp,noff,qm,idimp,npmax,nxv,nypmx)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells, for distributed data
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c q(j,k) = charge density at grid point (j,kk),
c where kk = k + noff - 1
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nxv = first dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nxv, nypmx
      real qm
      real part, q
      dimension part(idimp,npmax), q(nxv,nypmx)
c local data
      integer mnoff, j, nn, np, mm, mp
      real dxp, dyp, amx, amy
      mnoff = noff - 1
      do 10 j = 1, npp
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = qm - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
c deposit charge
      q(np,mp) = q(np,mp) + dxp*dyp
      q(nn,mp) = q(nn,mp) + amx*dyp
      q(np,mm) = q(np,mm) + dxp*amy
      q(nn,mm) = q(nn,mm) + amx*amy
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPDSORTP2YL(parta,partb,npic,npp,noff,nyp,idimp,npmax, 
     1nypm1)
c this subroutine sorts particles by y grid
c linear interpolation, spatial decomposition in y direction
c parta/partb = input/output particle array
c part(2,n) = position y of particle n in partition
c npic = address offset for reordering particles
c npp = number of particles in partition
c noff = backmost global gridpoint in particle partition
c nyp = number of primary gridpoints in particle partition
c idimp = size of phase space
c npmax = maximum number of particles in each partition
c nypm1 = maximum size of particle partition plus one
      implicit none
      integer idimp, npmax, nypm1
      integer npic, npp, noff, nyp
      real parta, partb
      dimension parta(idimp,npmax), partb(idimp,npmax)
      dimension npic(nypm1)
c local data
      integer i, j, k, m, mnoff, nyp1, isum, ist, ip
      mnoff = noff - 1
      nyp1 = nyp + 1
c clear counter array
      do 10 k = 1, nyp1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, npp
      m = parta(2,j)
      m = m - mnoff
      npic(m) = npic(m) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyp1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, npp
      m = parta(2,j)
      m = m - mnoff
      ip = npic(m) + 1
      do 40 i = 1, idimp
      partb(i,ip) = parta(i,j)
   40 continue
      npic(m) = ip
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPCGUARD2XL(fxy,nyp,nx,ndim,nxe,nypmx)
c replicate extended periodic vector field in x direction
c linear interpolation, for distributed data
c nyp = number of primary (complete) gridpoints in particle partition
c nx = system length in x direction
c ndim = leading dimension of array fxy
c nxe = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells
      implicit none
      integer nyp, nx, ndim, nxe, nypmx
      real fxy
      dimension fxy(ndim,nxe,nypmx)
c local data
      integer i, k, myp1
c replicate edges of extended field
      myp1 = nyp + 1
      do 20 k = 1, myp1
      do 10 i = 1, ndim
      fxy(i,nx+1,k) = fxy(i,1,k)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPAGUARD2XL(q,nyp,nx,nxe,nypmx)
c accumulate extended periodic scalar field in x direction
c linear interpolation, for distributed data
c nyp = number of primary (complete) gridpoints in particle partition
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells
      implicit none
      integer nyp, nx, nxe, nypmx
      real q
      dimension q(nxe,nypmx)
c local data
      integer k, myp1
c accumulate edges of extended field
      myp1 = nyp + 1
      do 10 k = 1, myp1
      q(1,k) = q(1,k) + q(nx+1,k)
      q(nx+1,k) = 0.0
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv, 
     1kxp,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions, for distributed data.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp,nyhd,
c output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,kstrt,nyv,kxp,nyhd,
c output: fxy,we
c approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c q(k,j) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j) = x component of complex force/charge,
c fxy(2,k,j) = y component of complex force/charge,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c kxp = number of data values per block
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated.
c aimag(ffc(k,j)) = finite-size particle shape factor s
c real(ffc(k,j)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp, nyhd
      real ax, ay, affp, we
      complex q, fxy, ffc
      dimension q(nyv,kxp), fxy(2,nyv,kxp)
      dimension ffc(nyhd,kxp)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2
      double precision wp
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
      if (isign.ne.0) go to 30
      if (kstrt.gt.nxh) return
c prepare form factor array
      do 20 j = 1, kxps
      dkx = dnx*real(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.0) then
         ffc(k,j) = cmplx(affp,1.0)
      else
         ffc(k,j) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
      if (kstrt.gt.nxh) go to 70
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 j = 1, kxps
      dkx = dnx*real(j + joff)
      if ((j+joff).gt.0) then
         do 40 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j))*aimag(ffc(k,j))
         at2 = dkx*at1
         at3 = dny*real(k - 1)*at1
         zt1 = cmplx(aimag(q(k,j)),-real(q(k,j)))
         zt2 = cmplx(aimag(q(k1,j)),-real(q(k1,j)))
         fxy(1,k,j) = at2*zt1
         fxy(2,k,j) = at3*zt1
         fxy(1,k1,j) = at2*zt2
         fxy(2,k1,j) = -at3*zt2
         wp = wp + at1*(q(k,j)*conjg(q(k,j)) + q(k1,j)*conjg(q(k1,j)))
   40    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffc(1,j))*aimag(ffc(1,j))
         at3 = dkx*at1
         zt1 = cmplx(aimag(q(1,j)),-real(q(1,j)))
         fxy(1,1,j) = at3*zt1
         fxy(2,1,j) = zero
         fxy(1,k1,j) = zero
         fxy(2,k1,j) = zero
         wp = wp + at1*(q(1,j)*conjg(q(1,j)))
      endif
   50 continue
c mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         do 60 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,1))*aimag(ffc(k,1))
         at2 = dny*real(k - 1)*at1
         zt1 = cmplx(aimag(q(k,1)),-real(q(k,1)))
         fxy(1,k,1) = zero
         fxy(2,k,1) = at2*zt1
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
         wp = wp + at1*(q(k,1)*conjg(q(k,1)))
   60    continue
         k1 = nyh + 1
         fxy(1,1,1) = zero
         fxy(2,1,1) = zero
         fxy(1,k1,1) = zero
         fxy(2,k1,1) = zero
      endif
   70 continue
      we = real(nx)*real(ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
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
      subroutine WPPFFT2R(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy
     1,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv, kxp, kyp
      integer kypd, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(nxvh,kypd), g(nyv,kxp)
      dimension bs(kxp,kyp), br(kxp,kyp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi, ks, kxpp, kypp
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      ks = kstrt - 1
      kxpp = min(kxp,max(0,nxh-kxp*ks))
      kypp = min(kyp,max(0,ny-kyp*ks))
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PPFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh
     1,kypd,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp,  
     1kypd)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PPFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,
     1kxp,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,   
     1kypd,kxp)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp
     1,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PPFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,
     1kxp,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,kypd, 
     1kxp)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PPFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh
     1,kypd,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPPFFT2R2(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,ind
     1y,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv, kxp, kyp
      integer kypd, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(2,nxvh,kypd), g(2,nyv,kxp)
      dimension bs(2,kxp,kyp), br(2,kxp,kyp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi, ks, kxpp, kypp
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      ks = kstrt - 1
      kxpp = min(kxp,max(0,nxh-kxp*ks))
      kypp = min(kyp,max(0,ny-kyp*ks))
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PPFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,   
     1nxvh,kypd,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,2,nxvh,nyv,kxp
     1,kypd)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PPFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv
     1,kxp,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,2,nyv,nxvh,
     1kypd,kxp)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,2,nxvh,nyv,
     1kxp,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PPFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv
     1,kxp,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,2,nyv,nxvh,   
     1kypd,kxp)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PPFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,   
     1nxvh,kypd,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp, 
     1nxvh,kypd,nxhyd,nxyhd)
c this subroutine performs the x part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kypd = second dimension of f
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k) = mode j-1,kk-1, where kk = k + kyp*(kstrt - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1) = real part of mode nx/2,0 on mode kstrt=0
c imaginary part of f(1,1) = real part of mode nx/2,ny/2
c on mode kstrt=(ny/2)/kyp
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, kstrt, nxvh, kypi, kypp, kypd
      integer nxhyd, nxyhd, mixup
      complex f, sct
      dimension f(nxvh,kypd)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, kypt, j, k, nrx
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      real ani
      complex s, t, t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 100
c inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 k = kypi, kypt
      t = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t
   10 continue
   20 continue
c first transform in x
      nrx = nxy/nxh
      do 60 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = kypi, kypt
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 80 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = kypi, kypt
      t = conjg(f(nxh2-j,k))
      s = f(j,k) + t
      t = (f(j,k) - t)*t1
      f(j,k) = ani*(s + t)
      f(nxh2-j,k) = ani*conjg(s - t)
   70 continue
   80 continue
      do 90 k = kypi, kypt
      f(1,k) = 2.0*ani*cmplx(real(f(1,k)) + aimag(f(1,k)),              
     1                       real(f(1,k)) - aimag(f(1,k)))
      if (nxhh.gt.0) f(nxhh+1,k) = 2.0*ani*conjg(f(nxhh+1,k))
   90 continue
      return
c forward fourier transform
  100 kmr = nxy/nx
c scramble coefficients
      do 120 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 k = kypi, kypt
      t = conjg(f(nxh2-j,k))
      s = f(j,k) + t
      t = (f(j,k) - t)*t1
      f(j,k) = s + t
      f(nxh2-j,k) = conjg(s - t)
  110 continue
  120 continue
      do 130 k = kypi, kypt
      f(1,k) = cmplx(real(f(1,k)) + aimag(f(1,k)),                      
     1               real(f(1,k)) - aimag(f(1,k)))
      if (nxhh.gt.0) f(nxhh+1,k) = 2.0*conjg(f(nxhh+1,k))
  130 continue
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 150 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 150
      do 140 k = kypi, kypt
      t = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t
  140 continue
  150 continue
c then transform in x
      nrx = nxy/nxh
      do 190 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 160 i = kypi, kypt
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp, 
     1nyv,kxp,nxhyd,nxyhd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c g(m,n) = sum(g(k,j)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(k,j) = sum(g(m,n)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c kxp = number of x indices per block
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g
c kxp = number of data values per block in x
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(k,j) = mode jj-1,k-1, where jj = j + kxp*(kstrt - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1) = real part of mode nx/2,ny/2
c on node kstrt=0
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxp
      integer nxhyd, nxyhd, mixup
      complex g, sct
      dimension g(nyv,kxp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      complex s, t
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 80
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = kxpi, kxpt
      t = g(k1,j)
      g(k1,j) = g(k,j)
      g(k,j) = t
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
      do 60 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = kxpi, kxpt
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      if (ks.gt.0) return
      do 70 k = 2, nyh
      if (kxpi.eq.1) then
         s = g(ny2-k,1)
         g(ny2-k,1) = 0.5*cmplx(aimag(g(k,1) + s),real(g(k,1) - s))
         g(k,1) = 0.5*cmplx(real(g(k,1) + s),aimag(g(k,1) - s))
      endif
   70 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
   80 nry = nxhy/ny
      if (ks.gt.0) go to 100
      do 90 k = 2, nyh
      if (kxpi.eq.1) then
         s = cmplx(aimag(g(ny2-k,1)),real(g(ny2-k,1)))
         g(ny2-k,1) = conjg(g(k,1) - s)
         g(k,1) = g(k,1) + s
      endif
   90 continue
c bit-reverse array elements in y
  100 do 120 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 120
      do 110 j = kxpi, kxpt
      t = g(k1,j)
      g(k1,j) = g(k,j)
      g(k,j) = t
  110 continue
  120 continue
c first transform in y
      nry = nxy/ny
      do 160 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 150 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 140 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 130 i = kxpi, kxpt
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
  130 continue
  140 continue
  150 continue
  160 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,
     1nxvh,kypd,nxhyd,nxyhd)
c this subroutine performs the x part of 2 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(1:2,n,m) = (1/nx*ny)*sum(f(1:2,j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(1:2,j,k) = sum(f(1:2,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kypd = second dimension of f
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(1:2,j,k) = mode j-1,kk-1, where kk = k + kyp*(kstrt - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1:2,1,k) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1:2,1,1) = real part of mode nx/2,0
c on mode kstrt=0
c imaginary part of f(1:2,1,1) = real part of mode nx/2,ny/2
c on mode kstrt=(ny/2)/kyp
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, kstrt, nxvh, kypi, kypp, kypd
      integer nxhyd, nxyhd, mixup
      complex f, sct
      dimension f(2,nxvh,kypd)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, kypt, j, k, nrx
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      real ani, at1
      complex s, t, t1, t2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 140
c inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrx = nxhy/nxh
c swap complex components
      do 20 k = kypi, kypt
      do 10 j = 1, nxh
      at1 = aimag(f(1,j,k))
      f(1,j,k) = cmplx(real(f(1,j,k)),real(f(2,j,k)))
      f(2,j,k) = cmplx(at1,aimag(f(2,j,k)))
   10 continue
   20 continue
c bit-reverse array elements in x
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 k = kypi, kypt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
   30 continue
   40 continue
c first transform in x
      nrx = nxy/nxh
      do 80 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 50 i = kypi, kypt
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 110 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = kypi, kypt
      do 90 i = 1, 2
      t = conjg(f(i,nxh2-j,k))
      s = f(i,j,k) + t
      t = (f(i,j,k) - t)*t1
      f(i,j,k) = ani*(s + t)
      f(i,nxh2-j,k) = ani*conjg(s - t)
   90 continue
  100 continue
  110 continue
      do 130 k = kypi, kypt
      do 120 i = 1, 2
      f(i,1,k) = 2.0*ani*cmplx(real(f(i,1,k)) + aimag(f(i,1,k)),        
     1                         real(f(i,1,k)) - aimag(f(i,1,k)))
      if (nxhh.gt.0) f(i,nxhh+1,k) = 2.0*ani*conjg(f(i,nxhh+1,k))
  120 continue
  130 continue
      return
c forward fourier transform
  140 kmr = nxy/nx
c scramble coefficients
      do 170 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 160 k = kypi, kypt
      do 150 i = 1, 2
      t = conjg(f(i,nxh2-j,k))
      s = f(i,j,k) + t
      t = (f(i,j,k) - t)*t1
      f(i,j,k) = s + t
      f(i,nxh2-j,k) = conjg(s - t)
  150 continue
  160 continue
  170 continue
      do 190 k = kypi, kypt
      do 180 i = 1, 2
      f(i,1,k) = cmplx(real(f(i,1,k)) + aimag(f(i,1,k)),                
     1                 real(f(i,1,k)) - aimag(f(i,1,k)))
      if (nxhh.gt.0) f(i,nxhh+1,k) = 2.0*conjg(f(i,nxhh+1,k))
  180 continue
  190 continue
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 210 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 210
      do 200 k = kypi, kypt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
  200 continue
  210 continue
c then transform in x
      nrx = nxy/nxh
      do 250 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 220 i = kypi, kypt
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
  220 continue
  230 continue
  240 continue
  250 continue
c swap complex components
      do 270 k = kypi, kypt
      do 260 j = 1, nxh
      at1 = aimag(f(1,j,k))
      f(1,j,k) = cmplx(real(f(1,j,k)),real(f(2,j,k)))
      f(2,j,k) = cmplx(at1,aimag(f(2,j,k)))
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,
     1nyv,kxp,nxhyd,nxyhd)
c this subroutine performs the y part of 2 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c g(1:2,m,n) = sum(g(1:2,k,j)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(1:2,k,j) = sum(g(1:2,m,n)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g
c kxp = number of data values per block in x
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(1:2,k,j) = mode jj-1,k-1, where jj = j + kxp*(kstrt - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(1:2,k,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1:2,1,1) = real part of mode nx/2,0 and
c imaginary part of g(1:2,ny/2+1,1) = real part of mode nx/2,ny/2
c on node kstrt=0
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxp
      integer nxhyd, nxyhd, mixup
      complex g, sct
      dimension g(2,nyv,kxp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      complex s, t1, t2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 90
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = kxpi, kxpt
      t1 = g(1,k1,j)
      t2 = g(2,k1,j)
      g(1,k1,j) = g(1,k,j)
      g(2,k1,j) = g(2,k,j)
      g(1,k,j) = t1
      g(2,k,j) = t2
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
      do 60 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = kxpi, kxpt
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      if (ks.gt.0) return
      do 80 k = 2, nyh
      if (kxpi.eq.1) then
         do 70 i = 1, 2
         s = g(i,ny2-k,1)
         g(i,ny2-k,1) = 0.5*cmplx(aimag(g(i,k,1) + s),                  
     1                            real(g(i,k,1) - s))
         g(i,k,1) = 0.5*cmplx(real(g(i,k,1) + s),aimag(g(i,k,1) - s))
   70    continue
      endif
   80 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
   90 nry = nxhy/ny
      if (ks.gt.0) go to 120
      do 110 k = 2, nyh
      if (kxpi.eq.1) then
         do 100 i = 1, 2
         s = cmplx(aimag(g(i,ny2-k,1)),real(g(i,ny2-k,1)))
         g(i,ny2-k,1) = conjg(g(i,k,1) - s)
         g(i,k,1) = g(i,k,1) + s
  100    continue
      endif
  110 continue
c bit-reverse array elements in y
  120 do 140 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 140
      do 130 j = kxpi, kxpt
      t1 = g(1,k1,j)
      t2 = g(2,k1,j)
      g(1,k1,j) = g(1,k,j)
      g(2,k1,j) = g(2,k,j)
      g(1,k,j) = t1
      g(2,k,j) = t2
  130 continue
  140 continue
c first transform in y
      nry = nxy/ny
      do 180 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 150 i = kxpi, kxpt
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
  150 continue
  160 continue
  170 continue
  180 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGSPOST2L(part,q,npp,noff,qm,idimp,npmax,nxv,nxyp)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c q(j,k) = charge density at grid point (j,kk),
c where kk = k + noff - 1
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nxv = first dimension of charge array, must be >= nx+1
c nxyp = actual first dimension of charge array, must be >= nxv*nypmx
      implicit none
      integer npp, noff, idimp, npmax, nxv, nxyp
      real qm
      real part, q
      dimension part(idimp,npmax), q(nxyp)
c local data
      integer mnoff, j, nnn, mmn, nn, mm, mp
      real dxn, dyn, dxp, dyp, amx, amy, dx1
      if (npp.lt.1) return
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j)
      mmn = part(2,j)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j) - real(nnn)
      dyn = part(2,j) - real(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyp
      mmn = mmn - mnoff
c deposit charge
      dx1 = q(mp+1) + dxp*dyp
      dyp = q(mp) + amx*dyp
      dxp = q(mm+1) + dxp*amy
      amy = q(mm) + amx*amy
      q(mp+1) = dx1
      q(mp) = dyp
      q(mm+1) = dxp
      q(mm) = amy
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyn
c deposit charge
      q(mp+1) = q(mp+1) + dxp*dyn
      q(mp) = q(mp) + amx*dyn
      q(mm+1) = q(mm+1) + dxp*amy
      q(mm) = q(mm) + amx*amy
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGSPUSH2L(part,fxy,edges,npp,noff,ihole,qbm,dt,ek,nx, 
     1ny,idimp,npmax,nxv,nxyp,idps,ntmax,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions,
c also determines list of particles which are leaving this processor
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996)
c 42 flops/particle, 12 loads, 4 stores
c input: all except ihole, output: part, ihole, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = velocity vx of particle n in partition
c part(4,n) = velocity vy of particle n in partition
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c in other words, fxy are the convolutions of the electric field
c over the particle shape, where kk = k + noff - 1
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nxv = second virtual dimension of field array, must be >= nx+1
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c idps = number of partition boundaries
c ntmax = size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, nxv, nxyp, idps, ntmax
      integer ipbc
      real qbm, dt, ek
      real part, fxy, edges
      integer ihole
      dimension part(idimp,npmax), fxy(2,nxyp)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, j, nnn, mmn, nn, mm, mp, ih, nh, nop, nop1
      real qtm, edgelx, edgely, edgerx, edgery, dxn, dyn, dxp, dyp
      real amx, amy, dx, dy
      double precision sum1
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      ih = 0
      nh = 0
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      mmn = mmn - mnoff
      nop1 = npp - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1)
      mmn = part(2,j+1)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1) - real(nnn)
      dyn = part(2,j+1) - real(mmn)
      mm = mm + nn
      amx = 1.0 - dxp
      mp = mm + nxv
      amy = 1.0 - dyp
      mmn = mmn - mnoff
c find acceleration
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp))                        
     1   + amy*(dxp*fxy(1,mm+1) + amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp))                        
     1   + amy*(dxp*fxy(2,mm+1) + amx*fxy(2,mm))
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      nop = npp
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1.0 - dxn
      mp = mm + nxv
      amy = 1.0 - dyn
c find acceleration
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp))                        
     1   + amy*(dxn*fxy(1,mm+1) + amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp))                        
     1   + amy*(dxn*fxy(2,mm+1) + amx*fxy(2,mm))
c new velocity
      dx = part(3,nop) + qtm*dx
      dy = part(4,nop) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,nop))**2 + (dy + part(4,nop))**2
      part(3,nop) = dx
      part(4,nop) = dy
c new position
      dx = part(1,nop) + dx*dt
      dy = part(2,nop) + dy*dt
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c set end of file flag
      if (nh.gt.0) ih = -ih
      ihole(1) = ih
c normalize kinetic energy
      ek = ek + 0.125*sum1
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
