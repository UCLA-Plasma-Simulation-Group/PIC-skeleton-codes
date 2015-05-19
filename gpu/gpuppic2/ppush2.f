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
      subroutine PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,  
     1mx1,mxyp1,irc)
c this subroutine finds the maximum number of particles in each tile of
c mx, my to calculate size of segmented particle array ppart
c linear interpolation, spatial decomposition in y direction
c input: all except kpic, nppmx, output: kpic, nppmx
c part = input particle array
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c kpic = output number of particles per tile
c nppmx = return maximum number of particles in tile
c npp = number of particles in partition
c noff = backmost global gridpoint in particle partition
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c mx/my = number of grids in sorting cell in x and y
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
      integer kpic, npp, noff
      real part
      dimension part(idimp,npmax)
      dimension kpic(mxyp1)
c local data
      integer j, k, n, m, mnoff, isum, ist, npx, ierr
      mnoff = noff
      ierr = 0
c clear counter array
      do 10 k = 1, mxyp1
      kpic(k) = 0
   10 continue
c find how many particles in each tile
      do 20 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      m = n + mx1*m
      if (m.le.mxyp1) then
         kpic(m) = kpic(m) + 1
      else
         ierr = max(ierr,m-mxyp1)
      endif
   20 continue
c find maximum
      isum = 0
      npx = 0
      do 30 k = 1, mxyp1
      ist = kpic(k)
      npx = max(npx,ist)
      isum = isum + ist
   30 continue
      nppmx = npx
c check for errors
      if (ierr.gt.0) then
         irc = ierr
      else if (isum.ne.npp) then
         irc = -1
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPMOVIN2LT(part,ppart,kpic,npp,noff,nppmx,idimp,npmax,
     1mx,my,mx1,mxyp1,irc)
c this subroutine sorts particles by x,y grid in tiles of
c mx, my and copies to segmented array ppart
c linear interpolation, spatial decomposition in y direction
c input: all except ppart, kpic, output: ppart, kpic
c part/ppart = input/output particle arrays
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c kpic = output number of particles per tile
c nppmx = maximum number of particles in tile
c npp = number of particles in partition
c noff = backmost global gridpoint in particle partition
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c mx/my = number of grids in sorting cell in x and y
c mx1 = (system length in x direction - 1)/mx + 1
c mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
      integer kpic, npp, noff
      real part, ppart
      dimension part(idimp,npmax), ppart(nppmx,idimp,mxyp1)
      dimension kpic(mxyp1)
c local data
      integer i, j, k, n, m, mnoff, ip, ierr
      mnoff = noff
      ierr = 0
c clear counter array
      do 10 k = 1, mxyp1
      kpic(k) = 0
   10 continue
c find addresses of particles at each tile and reorder particles
      do 30 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      m = n + mx1*m
      ip = kpic(m) + 1
      if (ip.le.nppmx) then
         do 20 i = 1, idimp
         ppart(ip,i,m) = part(i,j)
   20    continue
      else
         ierr = max(ierr,ip-nppmx)
      endif
      kpic(m) = ip
   30 continue
      if (ierr.gt.0) irc = ierr
      return
      end
c-----------------------------------------------------------------------
      subroutine PPPCHECK2LT(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,  
     1mx1,myp1,irc)
c this subroutine performs a sanity check to make sure particles sorted
c by x,y grid in tiles of mx, my, are all within bounds.
c tiles are assumed to be arranged in 2D linear memory, and transposed
c input: all except irc
c output: irc
c ppart(n,1,k) = position x of particle n in tile k
c ppart(n,2,k) = position y of particle n in tile k
c kpic(k) = number of reordered output particles in tile k
c noff = lowermost global gridpoint in particle partition.
c nyp = number of primary (complete) gridpoints in particle partition
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx = system length in x direction
c mx/my = number of grids in sorting cell in x/y
c mx1 = (system length in x direction - 1)/mx + 1
c myp1 = (partition length in y direction - 1)/my + 1
c irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, mx1, myp1, irc
      real ppart
      integer kpic
      dimension ppart(nppmx,idimp,mx1*myp1)
      dimension kpic(mx1*myp1)
c local data
      integer mxyp1, noffp, moffp, nppp, j, k, ist, nn, mm, ierr
      real edgelx, edgely, edgerx, edgery, dx, dy
      mxyp1 = mx1*myp1
      ierr = 0
c sanity check
      noffp = 0
      moffp = 0
c loop over tiles
      do 20 k = 1, mxyp1
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
c loop over particles in tile
      do 10 j = 1, nppp
      dx = ppart(j,1,k)
      dy = ppart(j,2,k)
c find particles going out of bounds
      ist = 0
      if (dx.lt.edgelx) ist = 1
      if (dx.ge.edgerx) ist = 2
      if (dy.lt.edgely) ist = ist + 3
      if (dy.ge.edgery) ist = ist + 6
      if (ist.gt.0) ierr = k
   10 continue
      noffp = noffp + mx
      if (noffp.ge.(mx*mx1)) then
         noffp = 0
         moffp = moffp + my
      endif
   20 continue
      if (ierr.ne.0) irc = ierr
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOIS22T(qt,fxyt,isign,ffct,ax,ay,affp,we,nx,ny,kstrt, 
     1nyv,kxp1,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions, for distributed data.
c vector length is second dimension.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp,nyhd,
c output: ffct
c for isign /= 0, input: qt,ffct,isign,nx,ny,kstrt,nyv,kxp,nyhd,
c output: fxyt,we
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
c qt(k,j) = complex charge density for fourier mode (jj-1,k-1)
c fxyt(k,1,j) = x component of complex force/charge,
c fxyt(k,2,j) = y component of complex force/charge,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c kxp1 = number of data values per block for unpacked field data
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated.
c aimag(ffct(k,j)) = finite-size particle shape factor s
c real(ffct(k,j)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp1, nyhd
      real ax, ay, affp, we
      complex qt, fxyt, ffct
      dimension qt(nyv,kxp1), fxyt(nyv,2,kxp1)
      dimension ffct(nyhd,kxp1)
c local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, j0, j1, k1
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp1*ks
      j1 = nxh + 1
      kxps = min(kxp1,max(0,j1-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
      if (kstrt.gt.j1) return
c prepare form factor array
      do 20 j = 1, kxps
      j0 = j + joff
      dkx = dnx*real(j0)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      if ((j0.ge.0).and.(j0.lt.nxh)) then
         do 10 k = 1, nyh
         dky = dny*real(k - 1)
         at3 = dky*dky + at1
         at4 = exp(-.5*((dky*ay)**2 + at2))
         if (at3.eq.0.0) then
            ffct(k,j) = cmplx(affp,1.0)
         else
            ffct(k,j) = cmplx(affp*at4/at3,at4)
         endif
   10    continue
      endif
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
      if (kstrt.gt.j1) go to 80
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 j = 1, kxps
      j0 = j + joff
      dkx = dnx*real(j0)
      if ((j0.gt.0).and.(j0.lt.nxh)) then
         do 40 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffct(k,j))*aimag(ffct(k,j))
         at2 = dkx*at1
         at3 = dny*real(k - 1)*at1
         zt1 = cmplx(aimag(qt(k,j)),-real(qt(k,j)))
         zt2 = cmplx(aimag(qt(k1,j)),-real(qt(k1,j)))
         fxyt(k,1,j) = at2*zt1
         fxyt(k1,1,j) = at2*zt2
         fxyt(k,2,j) = at3*zt1
         fxyt(k1,2,j) = -at3*zt2
         wp = wp + at1*(qt(k,j)*conjg(qt(k,j))                          
     1                + qt(k1,j)*conjg(qt(k1,j)))
   40    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffct(1,j))*aimag(ffct(1,j))
         at3 = dkx*at1
         zt1 = cmplx(aimag(qt(1,j)),-real(qt(1,j)))
         fxyt(1,1,j) = at3*zt1
         fxyt(k1,1,j) = zero
         fxyt(1,2,j) = zero
         fxyt(k1,2,j) = zero
         wp = wp + at1*(qt(1,j)*conjg(qt(1,j)))
      endif
   50 continue
c mode numbers kx = 0
      if (ks.eq.0) then
         do 60 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffct(k,1))*aimag(ffct(k,1))
         at2 = dny*real(k - 1)*at1
         zt1 = cmplx(aimag(qt(k,1)),-real(qt(k,1)))
         fxyt(k,1,1) = zero
         fxyt(k1,1,1) = zero
         fxyt(k,2,1) = at2*zt1
         fxyt(k1,2,1) = at2*conjg(zt1)
         wp = wp + at1*(qt(k,1)*conjg(qt(k,1)))
   60    continue
         k1 = nyh + 1
         fxyt(1,1,1) = zero
         fxyt(k1,1,1) = zero
         fxyt(1,2,1) = zero
         fxyt(k1,2,1) = zero
      endif
c mode numbers kx = nx/2
      if (ks.eq.(nxh/kxp1)) then
         do 70 k = 1, ny
         fxyt(k,1,kxps) = zero
         fxyt(k,2,kxps) = zero
   70    continue
      endif
   80 continue
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
