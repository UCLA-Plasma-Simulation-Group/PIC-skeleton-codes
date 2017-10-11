c Fortran Library for Skeleton 3D Electrostatic OpenMP/Vector PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR3(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,idimp, 
     1nop,nx,ny,nz,ipbc)
c for 3d code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c part(4,n) = velocity vx of particle n
c part(5,n) = velocity vy of particle n
c part(6,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny/nz = system length in x/y/z direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, npz, idimp, nop, nx, ny, nz, ipbc
      real vtx, vty, vtz, vdx, vdy, vdz
      real part
      dimension part(idimp,nop)
c local data
      integer j, k, l, k1, l1, npxy, npxyz
      real edgelx, edgely, edgelz, at1, at2, at3, at4, at5
      real sum1, sum2, sum3
      double precision dsum1, dsum2, dsum3
      double precision ranorm
      npxy = npx*npy
      npxyz = npxy*npz
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      at3 = real(nz)/real(npz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
         at3 = real(nz-2)/real(npz)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 0.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
      endif
c uniform density profile
      do 30 l = 1, npz
      l1 = npxy*(l - 1)
      at5 = edgelz + at3*(real(l) - .5)
      do 20 k = 1, npy
      k1 = npx*(k - 1) + l1
      at4 = edgely + at2*(real(k) - .5)
      do 10 j = 1, npx
      part(1,j+k1) = edgelx + at1*(real(j) - .5)
      part(2,j+k1) = at4
      part(3,j+k1) = at5
   10 continue
   20 continue
   30 continue
c maxwellian velocity distribution
      do 40 j = 1, npxyz
      part(4,j) = vtx*ranorm()
      part(5,j) = vty*ranorm()
      part(6,j) = vtz*ranorm()
   40 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 50 j = 1, npxyz
      dsum1 = dsum1 + part(4,j)
      dsum2 = dsum2 + part(5,j)
      dsum3 = dsum3 + part(6,j)
   50 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1.0/real(npxyz)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 60 j = 1, npxyz
      part(4,j) = part(4,j) - sum1
      part(5,j) = part(5,j) - sum2
      part(6,j) = part(6,j) - sum3
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DBLKP3L(part,kpic,nppmx,idimp,nop,mx,my,mz,mx1,my1,    
     1mxyz1,irc)
c this subroutine finds the maximum number of particles in each tile of
c mx, my, mz to calculate size of segmented particle array ppart
c linear interpolation
c input: all except kpic, nppmx, output: kpic, nppmx
c part = input particle array
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c kpic = output number of particles per tile
c nppmx = return maximum number of particles in tile
c idimp = size of phase space = 6
c nop = number of particles
c mx/my/mz = number of grids in sorting cell in x, y and z
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kpic, nppmx, idimp, nop, mx, my, mz, mx1, my1, mxyz1, irc
      real part
      dimension part(idimp,nop), kpic(mxyz1)
c local data
      integer j, k, n, m, l, mxy1, isum, ist, npx, ierr
      ierr = 0
      mxy1 = mx1*my1
c clear counter array
      do 10 k = 1, mxyz1
      kpic(k) = 0
   10 continue
c find how many particles in each tile
      do 20 j = 1, nop
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = m/my
      l = part(3,j)
      l = l/mz
      m = n + mx1*m + mxy1*l
      if (m.le.mxyz1) then
         kpic(m) = kpic(m) + 1
      else
         ierr = max(ierr,m-mxyz1)
      endif
   20 continue
c find maximum
      isum = 0
      npx = 0
      do 30 k = 1, mxyz1
      ist = kpic(k)
      npx = max(npx,ist)
      isum = isum + ist
   30 continue
      nppmx = npx
c check for errors
      if (ierr.gt.0) then
         irc = ierr
      else if (isum.ne.nop) then
         irc = -1
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPMOVIN3LT(part,ppart,kpic,nppmx,idimp,nop,mx,my,mz,mx1
     1,my1,mxyz1,irc)
c this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
c and copies to segmented array ppart
c linear interpolation
c input: all except ppart, kpic, output: ppart, kpic
c part/ppart = input/output particle arrays
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppart(n,4,m) = velocity vx of particle n in tile m
c ppart(n,5,m) = velocity vy of particle n in tile m
c ppart(n,6,m) = velocity vz of particle n in tile m
c kpic = output number of particles per tile
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 6
c nop = number of particles
c mx/my/mz = number of grids in sorting cell in x, y and z
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kpic, nppmx, idimp, nop, mx, my, mz, mx1, my1, mxyz1, irc
      real part, ppart
      dimension part(idimp,nop), ppart(nppmx,idimp,mxyz1)
      dimension kpic(mxyz1)
c local data
      integer i, j, k, n, m, l, mxy1, ip, ierr
      ierr = 0
      mxy1 = mx1*my1
c clear counter array
      do 10 k = 1, mxyz1
      kpic(k) = 0
   10 continue
c find addresses of particles at each tile and reorder particles
      do 30 j = 1, nop
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = m/my
      l = part(3,j)
      l = l/mz
      m = n + mx1*m + mxy1*l
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
      subroutine PPMOVIN3LTP(part,ppart,kpic,kp,nppmx,idimp,nop,mx,my,mz
     1,mx1,my1,mxyz1,irc)
c this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
c and copies to segmented array ppart
c designed for NUMA architectures, where memory is associated with the
c processor which first writes a memory location.
c linear interpolation
c input: all except ppart, kpic, output: ppart, kpic
c part/ppart = input/output particle arrays
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppart(n,4,m) = velocity vx of particle n in tile m
c ppart(n,5,m) = velocity vy of particle n in tile m
c ppart(n,6,m) = velocity vz of particle n in tile m
c kpic = output number of particles per tile
c kp = original location of reordered particle
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 6
c nop = number of particles
c mx/my/mz = number of grids in sorting cell in x, y and z
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kpic, kp, nppmx, idimp, nop, mx, my, mz, mx1, my1, mxyz1
      integer irc
      real part, ppart
      dimension part(idimp,nop), ppart(nppmx,idimp,mxyz1)
      dimension kpic(mxyz1), kp(nppmx,mxyz1)
c local data
      integer i, j, k, n, m, l, mxy1, ip, npp, ierr
      ierr = 0
      mxy1 = mx1*my1
c clear counter array
      do 10 k = 1, mxyz1
      kpic(k) = 0
   10 continue
c find addresses of particles at each tile to reorder particles
      do 20 j = 1, nop
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = m/my
      l = part(3,j)
      l = l/mz
      m = n + mx1*m + mxy1*l
      ip = kpic(m) + 1
      if (ip.le.nppmx) then
         kp(ip,m) = j
      else
         ierr = max(ierr,ip-nppmx)
      endif
      kpic(m) = ip
   20 continue
c check for overflow
      if (ierr.gt.0) then
         irc = ierr
         return
      endif
c copy reordered particles
!$OMP PARALLEL DO PRIVATE(i,j,k,m,npp)
      do 50 k = 1, mxyz1
      npp = kpic(k)
      do 40 j = 1, npp
      m = kp(j,k)
      do 30 i = 1, idimp
      ppart(j,i,k) = part(i,m)
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPCHECK3LT(ppart,kpic,idimp,nppmx,nx,ny,nz,mx,my,mz,mx1
     1,my1,mz1,irc)
c this subroutine performs a sanity check to make sure particles sorted
c by x,y,z grid in tiles of mx, my, mz, are all within bounds.
c tiles are assumed to be arranged in 3D linear memory, and transposed
c input: all except irc
c output: irc
c ppart(n,1,l) = position x of particle n in tile l
c ppart(n,2,l) = position y of particle n in tile l
c ppart(n,3,l) = position a of particle n in tile l
c kpic(l) = number of reordered output particles in tile l
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = number of grids in sorting cell in x/y/z
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mz1 = (system length in z direction - 1)/mz + 1
c irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, my1, mz1, irc
      real ppart
      integer kpic
      dimension ppart(nppmx,idimp,mx1*my1*mz1)
      dimension kpic(mx1*my1*mz1)
c local data
      integer mxy1, mxyz1, noff, moff, loff, npp, j, k, l, nn, mm, ll
      integer ist
      real edgelx, edgely, edgelz, edgerx, edgery, edgerz, dx, dy, dz
      mxy1 = mx1*my1
      mxyz1 = mxy1*mz1
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,l,noff,moff,loff,npp,nn,mm,ll,ist,edgelx,edgely,     
!$OMP& edgelz,edgerx,edgery,edgerz,dx,dy,dz)
      do 20 l = 1, mxyz1
      loff = (l - 1)/mxy1
      k = l - mxy1*loff
      loff = mz*loff
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(l)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      ll = min(mz,nz-loff)
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
      edgelz = loff
      edgerz = loff + ll
c loop over particles in tile
      do 10 j = 1, npp
      dx = ppart(j,1,l)
      dy = ppart(j,2,l)
      dz = ppart(j,3,l)
c find particles going out of bounds
      ist = 0
      if (dx.lt.edgelx) ist = 1
      if (dx.ge.edgerx) ist = 2
      if (dy.lt.edgely) ist = ist + 3
      if (dy.ge.edgery) ist = ist + 6
      if (dz.lt.edgelz) ist = ist + 9
      if (dz.ge.edgerz) ist = ist + 18
      if (ist.gt.0) irc = l
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GPPUSH3LT(ppart,fxyz,kpic,qbm,dt,ek,idimp,nppmx,nx,ny, 
     1nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions.
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 94 flops/particle, 30 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppart(n,4,m) = velocity vx of particle n in tile m
c ppart(n,5,m) = velocity vy of particle n in tile m
c ppart(n,6,m) = velocity vz of particle n in tile m
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c kpic = number of particles per tile
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = system length in x/y/z direction
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = second dimension of field array, must be >= nx+1
c nyv = third dimension of field array, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nyv, nzv
      integer mx1, my1, mxyz1, ipbc
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic
      dimension ppart(nppmx,idimp,mxyz1), fxyz(4,nxv*nyv*nzv)
      dimension kpic(mxyz1)
c local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxy1, noff, moff, loff, npp
      integer i, j, k, l, nn, mm, ll, lxv, lxyv, nxyv
      real qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real sfxyz
      dimension sfxyz(4,MXV*MYV*MZV)
c     dimension sfxyz(4,(mx+1)*(my+1)*(mz+1))
      double precision sum1, sum2
      mxy1 = mx1*my1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nyv
      qtm = qbm*dt
      sum2 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      edgerz = real(nz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
         edgerz = real(nz-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      endif
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noff,moff,loff,npp,nn,mm,ll,x,y,z,dxp,dyp,dzp,amx
!$OMP& ,amy,amz,dx1,dx,dy,dz,vx,vy,vz,sum1,sfxyz)
!$OMP& REDUCTION(+:sum2)
      do 50 l = 1, mxyz1
      loff = (l - 1)/mxy1
      k = l - mxy1*loff
      loff = mz*loff
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(l)
c load local fields from global array
      do 30 k = 1, min(mz,nz-loff)+1
      do 20 j = 1, min(my,ny-moff)+1
      do 10 i = 1, min(mx,nx-noff)+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(1,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(2,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(3,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
   10 continue
   20 continue
   30 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 40 j = 1, npp
c find interpolation weights
      x = ppart(j,1,l)
      y = ppart(j,2,l)
      z = ppart(j,3,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noff + 1 + lxv*(mm - moff) + lxyv*(ll - loff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find acceleration
      dx = amx*sfxyz(1,nn) + amy*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + amy*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + amy*sfxyz(3,nn+1)
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
c new velocity
      dxp = ppart(j,4,l)
      dyp = ppart(j,5,l)
      dzp = ppart(j,6,l)
      vx = dxp + qtm*dx
      vy = dyp + qtm*dy
      vz = dzp + qtm*dz
c average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
c new position
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
      endif
c set new position
      ppart(j,1,l) = dx
      ppart(j,2,l) = dy
      ppart(j,3,l) = dz
c set new velocity
      ppart(j,4,l) = vx
      ppart(j,5,l) = vy
      ppart(j,6,l) = vz
   40 continue
      sum2 = sum2 + sum1
   50 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine GPPUSHF3LT(ppart,fxyz,kpic,ncl,ihole,qbm,dt,ek,idimp,  
     1nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax,irc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 94 flops/particle, 30 loads, 6 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppart(n,4,m) = velocity vx of particle n in tile m
c ppart(n,5,m) = velocity vy of particle n in tile m
c ppart(n,6,m) = velocity vz of particle n in tile m
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c kpic(l) = number of particles in tile l
c ncl(i,l) = number of particles going to destination i, tile l
c ihole(1,:,l) = location of hole in array left by departing particle
c ihole(2,:,l) = direction destination of particle leaving hole
c all for tile l
c ihole(1,1,l) = ih, number of holes left (error, if negative)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = system length in x/y/z direction
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = second dimension of field array, must be >= nx+1
c nyv = third dimension of field array, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nyv, nzv
      integer mx1, my1, mxyz1, ntmax, irc
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mxyz1), fxyz(4,nxv*nyv*nzv)
      dimension kpic(mxyz1), ncl(26,mxyz1)
      dimension ihole(2,ntmax+1,mxyz1)
c local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxy1, noff, moff, loff, npp
      integer i, j, k, l, ih, nh, nn, mm, ll, lxv, lxyv, nxyv
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real qtm, x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real sfxyz
      dimension sfxyz(4,MXV*MYV*MZV)
c     dimension sfxyz(4,(mx+1)*(my+1)*(mz+1))
      double precision sum1, sum2
      mxy1 = mx1*my1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nyv
      qtm = qbm*dt
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
      sum2 = 0.0d0
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noff,moff,loff,npp,nn,mm,ll,ih,nh,x,y,z,dxp,dyp, 
!$OMP& dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,edgelx,edgely,edgelz,edgerx
!$OMP& ,edgery,edgerz,sum1,sfxyz)
!$OMP& REDUCTION(+:sum2)
      do 60 l = 1, mxyz1
      loff = (l - 1)/mxy1
      k = l - mxy1*loff
      loff = mz*loff
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(l)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      ll = min(mz,nz-loff)
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
      edgelz = loff
      edgerz = loff + ll
      ih = 0
      nh = 0
c load local fields from global array
      do 30 k = 1, ll+1
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(1,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(2,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(3,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
   10 continue
   20 continue
   30 continue
c clear counters
      do 40 j = 1, 26
      ncl(j,l) = 0
   40 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 50 j = 1, npp
c find interpolation weights
      x = ppart(j,1,l)
      y = ppart(j,2,l)
      z = ppart(j,3,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noff + 1 + lxv*(mm - moff) + lxyv*(ll - loff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find acceleration
      dx = amx*sfxyz(1,nn) + amy*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + amy*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + amy*sfxyz(3,nn+1)
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
c new velocity
      dxp = ppart(j,4,l)
      dyp = ppart(j,5,l)
      dzp = ppart(j,6,l)
      vx = dxp + qtm*dx
      vy = dyp + qtm*dy
      vz = dzp + qtm*dz
c average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
c new position
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
c find particles going out of bounds
      mm = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) dz = dz - anz
         mm = mm + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               mm = mm + 9
            else
               dz = 0.0
            endif
         else
            mm = mm + 9
         endif
      endif
c set new position
      ppart(j,1,l) = dx
      ppart(j,2,l) = dy
      ppart(j,3,l) = dz
c set new velocity
      ppart(j,4,l) = vx
      ppart(j,5,l) = vy
      ppart(j,6,l) = vz
c increment counters
      if (mm.gt.0) then
         ncl(mm,l) = ncl(mm,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j
            ihole(2,ih+1,l) = mm
         else
            nh = 1
         endif
      endif
   50 continue
      sum2 = sum2 + sum1
c set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
   60 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine VGPPUSH3LT(ppart,fxyz,kpic,qbm,dt,ek,idimp,nppmx,nx,ny,
     1nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions.
c vectorizable/OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 94 flops/particle, 30 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppart(n,4,m) = velocity vx of particle n in tile m
c ppart(n,5,m) = velocity vy of particle n in tile m
c ppart(n,6,m) = velocity vz of particle n in tile m
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c kpic = number of particles per tile
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = system length in x/y/z direction
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = second dimension of field array, must be >= nx+1
c nyv = third dimension of field array, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nyv, nzv
      integer mx1, my1, mxyz1, ipbc
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic
      dimension ppart(nppmx,idimp,mxyz1), fxyz(4,nxv*nyv*nzv)
      dimension kpic(mxyz1)
c local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxy1, noff, moff, loff, npp, ipp, joff, nps
      integer i, j, k, l, m, nn, mm, ll, lxv, lxyv, nxyv
      real qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real sfxyz
      dimension sfxyz(4,MXV*MYV*MZV)
c     dimension sfxyz(4,(mx+1)*(my+1)*(mz+1))
c scratch arrays
      integer n
      real s, t
      dimension n(npblk), s(npblk,lvect), t(npblk,3)
!dir$ attributes align: 64:: n, s, t
      double precision sum1, sum2
      mxy1 = mx1*my1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nyv
      qtm = qbm*dt
      sum2 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      edgerz = real(nz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
         edgerz = real(nz-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      endif
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,m,noff,moff,loff,npp,ipp,joff,nps,nn,mm,ll,x,y,z,
!$OMP& dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,sum1,sfxyz,n,s,t)
!$OMP& REDUCTION(+:sum2)
      do 110 l = 1, mxyz1
      loff = (l - 1)/mxy1
      k = l - mxy1*loff
      loff = mz*loff
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(l)
c load local fields from global array
      do 30 k = 1, min(mz,nz-loff)+1
      do 20 j = 1, min(my,ny-moff)+1
!dir$ ivdep
      do 10 i = 1, min(mx,nx-noff)+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(1,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(2,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(3,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
   10 continue
   20 continue
   30 continue
      sum1 = 0.0d0
c loop over particles in tile
      ipp = npp/npblk
c outer loop over number of full blocks
      do 90 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
!dir$ vector aligned
      do 40 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,l)
      y = ppart(j+joff,2,l)
      z = ppart(j+joff,3,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn - noff + lxv*(mm - moff) + lxyv*(ll - loff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   40 continue
c find acceleration
      do 60 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      ll = nn + lxyv - 4
      k = ll + lxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      do 50 i = 1, lvect
      if (i.gt.6) then
         nn = k
      else if (i.gt.4) then
         nn = ll
      else if (i.gt.2) then
         nn = mm
      endif
      dx = dx + sfxyz(1,i+nn)*s(j,i)
      dy = dy + sfxyz(2,i+nn)*s(j,i)
      dz = dz + sfxyz(3,i+nn)*s(j,i)
   50 continue
      s(j,1) = dx
      s(j,2) = dy
      s(j,3) = dz
   60 continue
c new velocity
!dir$ vector aligned
      do 70 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
      dxp = ppart(j+joff,4,l)
      dyp = ppart(j+joff,5,l)
      dzp = ppart(j+joff,6,l)
      vx = dxp + qtm*s(j,1)
      vy = dyp + qtm*s(j,2)
      vz = dzp + qtm*s(j,3)
c average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
c new position
      s(j,1) = x + vx*dt
      s(j,2) = y + vy*dt
      s(j,3) = z + vz*dt
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
   70 continue
c check boundary conditions
!dir$ vector aligned
      do 80 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      dz = s(j,3)
      vx = s(j,4)
      vy = s(j,5)
      vz = s(j,6)
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = t(j,3)
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
      endif
c set new position
      ppart(j+joff,1,l) = dx
      ppart(j+joff,2,l) = dy
      ppart(j+joff,3,l) = dz
c set new velocity
      ppart(j+joff,4,l) = vx
      ppart(j+joff,5,l) = vy
      ppart(j+joff,6,l) = vz
   80 continue
   90 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 100 j = nps, npp
c find interpolation weights
      x = ppart(j,1,l)
      y = ppart(j,2,l)
      z = ppart(j,3,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noff + 1 + lxv*(mm - moff) + lxyv*(ll - loff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find acceleration
      dx = amx*sfxyz(1,nn) + amy*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + amy*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + amy*sfxyz(3,nn+1)
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
c new velocity
      dxp = ppart(j,4,l)
      dyp = ppart(j,5,l)
      dzp = ppart(j,6,l)
      vx = dxp + qtm*dx
      vy = dyp + qtm*dy
      vz = dzp + qtm*dz
c average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
c new position
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
      endif
c set new position
      ppart(j,1,l) = dx
      ppart(j,2,l) = dy
      ppart(j,3,l) = dz
c set new velocity
      ppart(j,4,l) = vx
      ppart(j,5,l) = vy
      ppart(j,6,l) = vz
  100 continue
      sum2 = sum2 + sum1
  110 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine VGPPUSHF3LT(ppart,fxyz,kpic,ncl,ihole,qbm,dt,ek,idimp, 
     1nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax,irc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c vectorizable/OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 94 flops/particle, 30 loads, 6 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppart(n,4,m) = velocity vx of particle n in tile m
c ppart(n,5,m) = velocity vy of particle n in tile m
c ppart(n,6,m) = velocity vz of particle n in tile m
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c kpic(l) = number of particles in tile l
c ncl(i,l) = number of particles going to destination i, tile l
c ihole(1,:,l) = location of hole in array left by departing particle
c ihole(2,:,l) = direction destination of particle leaving hole
c all for tile l
c ihole(1,1,l) = ih, number of holes left (error, if negative)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = system length in x/y/z direction
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = second dimension of field array, must be >= nx+1
c nyv = third dimension of field array, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nyv, nzv
      integer mx1, my1, mxyz1, ntmax, irc
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mxyz1), fxyz(4,nxv*nyv*nzv)
      dimension kpic(mxyz1), ncl(26,mxyz1)
      dimension ihole(2,ntmax+1,mxyz1)
c local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxy1, noff, moff, loff, npp, ipp, joff, nps
      integer i, j, k, l, m, ih, nh, nn, mm, ll, lxv, lxyv, nxyv
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real qtm, x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real sfxyz
      dimension sfxyz(4,MXV*MYV*MZV)
c     dimension sfxyz(4,(mx+1)*(my+1)*(mz+1))
c scratch arrays
      integer n
      real s, t
      dimension n(npblk), s(npblk,lvect), t(npblk,3)
!dir$ attributes align: 64:: n, s, t
      double precision sum1, sum2
      mxy1 = mx1*my1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nyv
      qtm = qbm*dt
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
      sum2 = 0.0d0
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,m,noff,moff,loff,npp,ipp,joff,nps,nn,mm,ll,ih,nh,
!$OMP& x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,edgelx,edgely
!$OMP& ,edgelz,edgerx,edgery,edgerz,sum1,sfxyz,n,s,t)
!$OMP& REDUCTION(+:sum2)
      do 130 l = 1, mxyz1
      loff = (l - 1)/mxy1
      k = l - mxy1*loff
      loff = mz*loff
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(l)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      ll = min(mz,nz-loff)
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
      edgelz = loff
      edgerz = loff + ll
      ih = 0
      nh = 0
c load local fields from global array
      do 30 k = 1, ll+1
      do 20 j = 1, mm+1
!dir$ ivdep
      do 10 i = 1, nn+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(1,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(2,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(3,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
   10 continue
   20 continue
   30 continue
c clear counters
      do 40 j = 1, 26
      ncl(j,l) = 0
   40 continue
      sum1 = 0.0d0
c loop over particles in tile
      ipp = npp/npblk
c outer loop over number of full blocks
      do 110 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
!dir$ vector aligned
      do 50 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,l)
      y = ppart(j+joff,2,l)
      z = ppart(j+joff,3,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn - noff + lxv*(mm - moff) + lxyv*(ll - loff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   50 continue
c find acceleration
      do 70 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      ll = nn + lxyv - 4
      k = ll + lxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      do 60 i = 1, lvect
      if (i.gt.6) then
         nn = k
      else if (i.gt.4) then
         nn = ll
      else if (i.gt.2) then
         nn = mm
      endif
      dx = dx + sfxyz(1,i+nn)*s(j,i)
      dy = dy + sfxyz(2,i+nn)*s(j,i)
      dz = dz + sfxyz(3,i+nn)*s(j,i)
   60 continue
      s(j,1) = dx
      s(j,2) = dy
      s(j,3) = dz
   70 continue
c new velocity
!dir$ vector aligned
      do 80 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
      dxp = ppart(j+joff,4,l)
      dyp = ppart(j+joff,5,l)
      dzp = ppart(j+joff,6,l)
      vx = dxp + qtm*s(j,1)
      vy = dyp + qtm*s(j,2)
      vz = dzp + qtm*s(j,3)
c average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
c new position
      s(j,1) = x + vx*dt
      s(j,2) = y + vy*dt
      s(j,3) = z + vz*dt
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
   80 continue
c check boundary conditions
!dir$ vector aligned
      do 90 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      dz = s(j,3)
c find particles going out of bounds
      mm = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) dz = dz - anz
         mm = mm + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               mm = mm + 9
            else
               dz = 0.0
            endif
         else
            mm = mm + 9
         endif
      endif
c set new position
      ppart(j+joff,1,l) = dx
      ppart(j+joff,2,l) = dy
      ppart(j+joff,3,l) = dz
c set new velocity
      ppart(j+joff,4,l) = s(j,4)
      ppart(j+joff,5,l) = s(j,5)
      ppart(j+joff,6,l) = s(j,6)
      n(j) = mm
   90 continue
c increment counters
      do 100 j = 1, npblk
      mm = n(j)
      if (mm.gt.0) then
         ncl(mm,l) = ncl(mm,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j + joff
            ihole(2,ih+1,l) = mm
         else
            nh = 1
         endif
      endif
  100 continue
  110 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 120 j = nps, npp
c find interpolation weights
      x = ppart(j,1,l)
      y = ppart(j,2,l)
      z = ppart(j,3,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noff + 1 + lxv*(mm - moff) + lxyv*(ll - loff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find acceleration
      dx = amx*sfxyz(1,nn) + amy*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + amy*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + amy*sfxyz(3,nn+1)
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
c new velocity
      dxp = ppart(j,4,l)
      dyp = ppart(j,5,l)
      dzp = ppart(j,6,l)
      vx = dxp + qtm*dx
      vy = dyp + qtm*dy
      vz = dzp + qtm*dz
c average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
c new position
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
c find particles going out of bounds
      mm = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) dz = dz - anz
         mm = mm + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               mm = mm + 9
            else
               dz = 0.0
            endif
         else
            mm = mm + 9
         endif
      endif
c set new position
      ppart(j,1,l) = dx
      ppart(j,2,l) = dy
      ppart(j,3,l) = dz
c set new velocity
      ppart(j,4,l) = vx
      ppart(j,5,l) = vy
      ppart(j,6,l) = vz
c increment counters
      if (mm.gt.0) then
         ncl(mm,l) = ncl(mm,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j
            ihole(2,ih+1,l) = mm
         else
            nh = 1
         endif
      endif
  120 continue
      sum2 = sum2 + sum1
c set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
  130 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine V2GPPUSH3LT(ppart,fxyz,kpic,qbm,dt,ek,idimp,nppmx,nx,ny
     1,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions.
c vectorizable/OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 94 flops/particle, 30 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppart(n,4,m) = velocity vx of particle n in tile m
c ppart(n,5,m) = velocity vy of particle n in tile m
c ppart(n,6,m) = velocity vz of particle n in tile m
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c kpic = number of particles per tile
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = system length in x/y/z direction
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = second dimension of field array, must be >= nx+1
c nyv = third dimension of field array, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nyv, nzv
      integer mx1, my1, mxyz1, ipbc
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic
      dimension ppart(nppmx,idimp,mxyz1), fxyz(4,nxv*nyv*nzv)
      dimension kpic(mxyz1)
c local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxy1, noff, moff, loff, npp, ipp, joff, nps
      integer i, j, k, l, m, nn, mm, ll, lxv, lxyv, nxyv
      real qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real sfxyz
      dimension sfxyz(4,MXV*MYV*MZV)
c     dimension sfxyz(4,(mx+1)*(my+1)*(mz+1))
c scratch arrays
      integer n, mn
      real s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,3)
!dir$ attributes align: 64:: n, mn, s, t
      double precision sum1, sum2
      mxy1 = mx1*my1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nyv
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      mn(5) = lxyv
      mn(6) = lxyv + 1
      mn(7) = lxyv + lxv
      mn(8) = lxyv + lxv + 1
      qtm = qbm*dt
      sum2 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      edgerz = real(nz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
         edgerz = real(nz-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      endif
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,m,noff,moff,loff,npp,ipp,joff,nps,nn,mm,ll,x,y,z,
!$OMP& dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,sum1,sfxyz,n,s,t)
!$OMP& REDUCTION(+:sum2)
      do 110 l = 1, mxyz1
      loff = (l - 1)/mxy1
      k = l - mxy1*loff
      loff = mz*loff
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(l)
c load local fields from global array
      do 30 k = 1, min(mz,nz-loff)+1
      do 20 j = 1, min(my,ny-moff)+1
!dir$ ivdep
      do 10 i = 1, min(mx,nx-noff)+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(1,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(2,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(3,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
   10 continue
   20 continue
   30 continue
      sum1 = 0.0d0
c loop over particles in tile
      ipp = npp/npblk
c outer loop over number of full blocks
      do 90 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
!dir$ vector aligned
      do 40 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,l)
      y = ppart(j+joff,2,l)
      z = ppart(j+joff,3,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn - noff + lxv*(mm - moff) + lxyv*(ll - loff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   40 continue
c find acceleration
      do 60 j = 1, npblk
      dx = 0.0
      dy = 0.0
      dz = 0.0
!dir$ ivdep
      do 50 i = 1, lvect
      dx = dx + sfxyz(1,n(j)+mn(i))*s(j,i)
      dy = dy + sfxyz(2,n(j)+mn(i))*s(j,i)
      dz = dz + sfxyz(3,n(j)+mn(i))*s(j,i)
   50 continue
      s(j,1) = dx
      s(j,2) = dy
      s(j,3) = dz
   60 continue
c new velocity
!dir$ vector aligned
      do 70 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
      dxp = ppart(j+joff,4,l)
      dyp = ppart(j+joff,5,l)
      dzp = ppart(j+joff,6,l)
      vx = dxp + qtm*s(j,1)
      vy = dyp + qtm*s(j,2)
      vz = dzp + qtm*s(j,3)
c average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
c new position
      s(j,1) = x + vx*dt
      s(j,2) = y + vy*dt
      s(j,3) = z + vz*dt
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
   70 continue
c check boundary conditions
!dir$ vector aligned
      do 80 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      dz = s(j,3)
      vx = s(j,4)
      vy = s(j,5)
      vz = s(j,6)
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = t(j,3)
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = t(j,1)
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = t(j,2)
            vy = -vy
         endif
      endif
c set new position
      ppart(j+joff,1,l) = dx
      ppart(j+joff,2,l) = dy
      ppart(j+joff,3,l) = dz
c set new velocity
      ppart(j+joff,4,l) = vx
      ppart(j+joff,5,l) = vy
      ppart(j+joff,6,l) = vz
   80 continue
   90 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 100 j = nps, npp
c find interpolation weights
      x = ppart(j,1,l)
      y = ppart(j,2,l)
      z = ppart(j,3,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noff + 1 + lxv*(mm - moff) + lxyv*(ll - loff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find acceleration
      dx = amx*sfxyz(1,nn) + amy*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + amy*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + amy*sfxyz(3,nn+1)
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
c new velocity
      dxp = ppart(j,4,l)
      dyp = ppart(j,5,l)
      dzp = ppart(j,6,l)
      vx = dxp + qtm*dx
      vy = dyp + qtm*dy
      vz = dzp + qtm*dz
c average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
c new position
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = z
            vz = -vz
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            vx = -vx
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            vy = -vy
         endif
      endif
c set new position
      ppart(j,1,l) = dx
      ppart(j,2,l) = dy
      ppart(j,3,l) = dz
c set new velocity
      ppart(j,4,l) = vx
      ppart(j,5,l) = vy
      ppart(j,6,l) = vz
  100 continue
      sum2 = sum2 + sum1
  110 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine V2GPPUSHF3LT(ppart,fxyz,kpic,ncl,ihole,qbm,dt,ek,idimp,
     1nppmx,nx,ny,nz,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1,ntmax,irc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c vectorizable/OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 94 flops/particle, 30 loads, 6 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppart(n,4,m) = velocity vx of particle n in tile m
c ppart(n,5,m) = velocity vy of particle n in tile m
c ppart(n,6,m) = velocity vz of particle n in tile m
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c kpic(l) = number of particles in tile l
c ncl(i,l) = number of particles going to destination i, tile l
c ihole(1,:,l) = location of hole in array left by departing particle
c ihole(2,:,l) = direction destination of particle leaving hole
c all for tile l
c ihole(1,1,l) = ih, number of holes left (error, if negative)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = system length in x/y/z direction
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = second dimension of field array, must be >= nx+1
c nyv = third dimension of field array, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, nxv, nyv, nzv
      integer mx1, my1, mxyz1, ntmax, irc
      real qbm, dt, ek
      real ppart, fxyz
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mxyz1), fxyz(4,nxv*nyv*nzv)
      dimension kpic(mxyz1), ncl(26,mxyz1)
      dimension ihole(2,ntmax+1,mxyz1)
c local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxy1, noff, moff, loff, npp, ipp, joff, nps
      integer i, j, k, l, m, ih, nh, nn, mm, ll, lxv, lxyv, nxyv
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real qtm, x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      real vx, vy, vz
      real sfxyz
      dimension sfxyz(4,MXV*MYV*MZV)
c     dimension sfxyz(4,(mx+1)*(my+1)*(mz+1))
c scratch arrays
      integer n, mn
      real s, t
      dimension n(npblk), mn(lvect), s(npblk,lvect), t(npblk,3)
!dir$ attributes align: 64:: n, mn, s, t
      double precision sum1, sum2
      mxy1 = mx1*my1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nyv
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      mn(5) = lxyv
      mn(6) = lxyv + 1
      mn(7) = lxyv + lxv
      mn(8) = lxyv + lxv + 1
      qtm = qbm*dt
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
      sum2 = 0.0d0
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,m,noff,moff,loff,npp,ipp,joff,nps,nn,mm,ll,ih,nh,
!$OMP& x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,edgelx,edgely
!$OMP& ,edgelz,edgerx,edgery,edgerz,sum1,sfxyz,n,s,t)
!$OMP& REDUCTION(+:sum2)
      do 130 l = 1, mxyz1
      loff = (l - 1)/mxy1
      k = l - mxy1*loff
      loff = mz*loff
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(l)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      ll = min(mz,nz-loff)
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
      edgelz = loff
      edgerz = loff + ll
      ih = 0
      nh = 0
c load local fields from global array
      do 30 k = 1, ll+1
      do 20 j = 1, mm+1
!dir$ ivdep
      do 10 i = 1, nn+1
      sfxyz(1,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(1,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(2,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(2,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
      sfxyz(3,i+lxv*(j-1)+lxyv*(k-1)) = 
     1 fxyz(3,i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))
   10 continue
   20 continue
   30 continue
c clear counters
      do 40 j = 1, 26
      ncl(j,l) = 0
   40 continue
      sum1 = 0.0d0
c loop over particles in tile
      ipp = npp/npblk
c outer loop over number of full blocks
      do 110 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
!dir$ vector aligned
      do 50 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,l)
      y = ppart(j+joff,2,l)
      z = ppart(j+joff,3,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn - noff + lxv*(mm - moff) + lxyv*(ll - loff) + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
      t(j,1) = x
      t(j,2) = y
      t(j,3) = z
   50 continue
c find acceleration
      do 70 j = 1, npblk
      dx = 0.0
      dy = 0.0
      dz = 0.0
!dir$ ivdep
      do 60 i = 1, lvect
      dx = dx + sfxyz(1,n(j)+mn(i))*s(j,i)
      dy = dy + sfxyz(2,n(j)+mn(i))*s(j,i)
      dz = dz + sfxyz(3,n(j)+mn(i))*s(j,i)
   60 continue
      s(j,1) = dx
      s(j,2) = dy
      s(j,3) = dz
   70 continue
c new velocity
!dir$ vector aligned
      do 80 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      z = t(j,3)
      dxp = ppart(j+joff,4,l)
      dyp = ppart(j+joff,5,l)
      dzp = ppart(j+joff,6,l)
      vx = dxp + qtm*s(j,1)
      vy = dyp + qtm*s(j,2)
      vz = dzp + qtm*s(j,3)
c average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
c new position
      s(j,1) = x + vx*dt
      s(j,2) = y + vy*dt
      s(j,3) = z + vz*dt
      s(j,4) = vx
      s(j,5) = vy
      s(j,6) = vz
   80 continue
c check boundary conditions
!dir$ vector aligned
      do 90 j = 1, npblk
      dx = s(j,1)
      dy = s(j,2)
      dz = s(j,3)
c find particles going out of bounds
      mm = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) dz = dz - anz
         mm = mm + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               mm = mm + 9
            else
               dz = 0.0
            endif
         else
            mm = mm + 9
         endif
      endif
c set new position
      ppart(j+joff,1,l) = dx
      ppart(j+joff,2,l) = dy
      ppart(j+joff,3,l) = dz
c set new velocity
      ppart(j+joff,4,l) = s(j,4)
      ppart(j+joff,5,l) = s(j,5)
      ppart(j+joff,6,l) = s(j,6)
      n(j) = mm
   90 continue
c increment counters
      do 100 j = 1, npblk
      mm = n(j)
      if (mm.gt.0) then
         ncl(mm,l) = ncl(mm,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j + joff
            ihole(2,ih+1,l) = mm
         else
            nh = 1
         endif
      endif
  100 continue
  110 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 120 j = nps, npp
c find interpolation weights
      x = ppart(j,1,l)
      y = ppart(j,2,l)
      z = ppart(j,3,l)
      nn = x
      mm = y
      ll = z
      dxp = x - real(nn)
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noff + 1 + lxv*(mm - moff) + lxyv*(ll - loff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c find acceleration
      dx = amx*sfxyz(1,nn) + amy*sfxyz(1,nn+1)
      dy = amx*sfxyz(2,nn) + amy*sfxyz(2,nn+1)
      dz = amx*sfxyz(3,nn) + amy*sfxyz(3,nn+1)
      dx = amz*(dx + dyp*sfxyz(1,nn+lxv) + dx1*sfxyz(1,nn+1+lxv))
      dy = amz*(dy + dyp*sfxyz(2,nn+lxv) + dx1*sfxyz(2,nn+1+lxv))
      dz = amz*(dz + dyp*sfxyz(3,nn+lxv) + dx1*sfxyz(3,nn+1+lxv))
      mm = nn + lxyv
      vx = amx*sfxyz(1,mm) + amy*sfxyz(1,mm+1)
      vy = amx*sfxyz(2,mm) + amy*sfxyz(2,mm+1)
      vz = amx*sfxyz(3,mm) + amy*sfxyz(3,mm+1)
      dx = dx + dzp*(vx + dyp*sfxyz(1,mm+lxv) + dx1*sfxyz(1,mm+1+lxv))
      dy = dy + dzp*(vy + dyp*sfxyz(2,mm+lxv) + dx1*sfxyz(2,mm+1+lxv))
      dz = dz + dzp*(vz + dyp*sfxyz(3,mm+lxv) + dx1*sfxyz(3,mm+1+lxv))
c new velocity
      dxp = ppart(j,4,l)
      dyp = ppart(j,5,l)
      dzp = ppart(j,6,l)
      vx = dxp + qtm*dx
      vy = dyp + qtm*dy
      vz = dzp + qtm*dz
c average kinetic energy
      dxp = dxp + vx
      dyp = dyp + vy
      dzp = dzp + vz
      sum1 = sum1 + (dxp*dxp + dyp*dyp + dzp*dzp)
c new position
      dx = x + vx*dt
      dy = y + vy*dt
      dz = z + vz*dt
c find particles going out of bounds
      mm = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) dz = dz - anz
         mm = mm + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               mm = mm + 9
            else
               dz = 0.0
            endif
         else
            mm = mm + 9
         endif
      endif
c set new position
      ppart(j,1,l) = dx
      ppart(j,2,l) = dy
      ppart(j,3,l) = dz
c set new velocity
      ppart(j,4,l) = vx
      ppart(j,5,l) = vy
      ppart(j,6,l) = vz
c increment counters
      if (mm.gt.0) then
         ncl(mm,l) = ncl(mm,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j
            ihole(2,ih+1,l) = mm
         else
            nh = 1
         endif
      endif
  120 continue
      sum2 = sum2 + sum1
c set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
  130 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine GPPOST3LT(ppart,q,kpic,qm,nppmx,idimp,mx,my,mz,nxv,nyv,
     1nzv,mx1,my1,mxyz1)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c OpenMP version using guard cells
c data deposited in tiles
c particles stored segmented array
c 33 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c q(j,k,l) = charge density at grid point j,k,l
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 6
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
c nzv = third dimension of charge array, must be >= nz+1
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
      implicit none
      integer nppmx, idimp, mx, my, mz, nxv, nyv, nzv, mx1, my1, mxyz1
      real qm
      real ppart, q
      integer kpic
      dimension ppart(nppmx,idimp,mxyz1), q(nxv*nyv*nzv)
      dimension kpic(mxyz1)
c local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer mxy1, noff, moff, loff, npp
      integer i, j, k, l, nn, mm, ll, nm, lm, lxv, lxyv, nxyv
      real x, y, z, w, dxp, dyp, dzp, amx, amy, amz, dx1
      real sq
c     dimension sq(MXV*MYV*MZV)
      dimension sq((mx+1)*(my+1)*(mz+1))
      mxy1 = mx1*my1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nyv
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noff,moff,loff,npp,nn,mm,ll,nm,lm,x,y,z,w,dxp,
!$OMP& dyp,dzp,amx,amy,amz,dx1,sq)
      do 130 l = 1, mxyz1
      loff = (l - 1)/mxy1
      k = l - mxy1*loff
      loff = mz*loff
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(l)
c zero out local accumulator
      do 10 j = 1, (mx+1)*(my+1)*(mz+1)
      sq(j) = 0.0
   10 continue
c loop over particles in tile
      do 20 j = 1, npp
c find interpolation weights
      x = ppart(j,1,l)
      y = ppart(j,2,l)
      z = ppart(j,3,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noff + 1 + lxv*(mm - moff) + lxyv*(ll - loff)
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c deposit charge within tile to local accumulator
      x = sq(nn) + amx*amz
      y = sq(nn+1) + amy*amz
      z = sq(nn+lxv) + dyp*amz
      w = sq(nn+1+lxv) + dx1*amz
      sq(nn) = x
      sq(nn+1) = y
      sq(nn+lxv) = z
      sq(nn+1+lxv) = w
      mm = nn + lxyv
      x = sq(mm) + amx*dzp
      y = sq(mm+1) + amy*dzp
      z = sq(mm+lxv) + dyp*dzp
      w = sq(mm+1+lxv) + dx1*dzp
      sq(mm) = x
      sq(mm+1) = y
      sq(mm+lxv) = z
      sq(mm+1+lxv) = w
   20 continue
c deposit charge to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      ll = min(mz,nzv-loff)
      do 50 k = 2, ll
      do 40 j = 2, mm
!dir$ ivdep
      do 30 i = 2, nn
      q(i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1)) =                        
     1q(i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1)) +                        
     2 sq(i+lxv*(j-1)+lxyv*(k-1))
   30 continue
   40 continue
   50 continue
c deposit charge to edge points in global array
      lm = min(mz+1,nzv-loff)
      do 70 j = 2, mm
      do 60 i = 2, nn
!$OMP ATOMIC
      q(i+noff+nxv*(j+moff-1)+nxyv*loff) =                              
     1q(i+noff+nxv*(j+moff-1)+nxyv*loff) + sq(i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         q(i+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1)) =                    
     1   q(i+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1))
     1   + sq(i+lxv*(j-1)+lxyv*(lm-1))
      endif
   60 continue
   70 continue
      nm = min(mx+1,nxv-noff)
      mm = min(my+1,nyv-moff)
      do 100 k = 1, ll
      do 80 i = 2, nn
!$OMP ATOMIC
      q(i+noff+nxv*moff+nxyv*(k+loff-1)) =                              
     1q(i+noff+nxv*moff+nxyv*(k+loff-1)) + sq(i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         q(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff-1)) =      
     1   q(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff-1))                      
     2   + sq(i+lxv*(mm-1)+lxyv*(k-1))
      endif
   80 continue
      do 90 j = 1, mm
!$OMP ATOMIC
      q(1+noff+nxv*(j+moff-1)+nxyv*(k+loff-1)) =                        
     1q(1+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))                          
     2+ sq(1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         q(nm+noff+nxv*(j+moff-1)+nxyv*(k+loff-1)) =                    
     1   q(nm+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))                      
     1   + sq(nm+lxv*(j-1)+lxyv*(k-1))
      endif
   90 continue
  100 continue
      if (lm > mz) then
         do 110 i = 2, nn
!$OMP ATOMIC
         q(i+noff+nxv*moff+nxyv*(lm+loff-1)) =                          
     1   q(i+noff+nxv*moff+nxyv*(lm+loff-1)) + sq(i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            q(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1)) =                
     1      q(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))                  
     1      + sq(i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  110    continue
         do 120 j = 1, mm
!$OMP ATOMIC
         q(1+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1)) =                    
     1   q(1+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1))                      
     1   + sq(1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            q(nm+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1)) =                
     1      q(nm+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1))                  
     1      + sq(nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  120    continue
      endif
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine VGPPOST3LT(ppart,q,kpic,qm,nppmx,idimp,mx,my,mz,nxv,nyv
     1,nzv,mx1,my1,mxyz1)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c vectorizable/OpenMP version using guard cells
c data deposited in tiles
c particles stored segmented array
c 33 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c q(j,k,l) = charge density at grid point j,k,l
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 6
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
c nzv = third dimension of charge array, must be >= nz+1
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
      implicit none
      integer nppmx, idimp, mx, my, mz, nxv, nyv, nzv, mx1, my1, mxyz1
      real qm
      real ppart, q
      integer kpic
      dimension ppart(nppmx,idimp,mxyz1), q(nxv*nyv*nzv)
      dimension kpic(mxyz1)
c local data
      integer MXV, MYV, MZV
      parameter(MXV=17,MYV=17,MZV=17)
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer mxy1, noff, moff, loff, npp, ipp, joff, nps
      integer i, j, k, l, m, nn, mm, ll, nm, lm, lxv, lxyv, nxyv
      real x, y, z, w, dxp, dyp, dzp, amx, amy, amz, dx1
      real sq
c     dimension sq(MXV*MYV*MZV)
      dimension sq((mx+1)*(my+1)*(mz+1))
c scratch arrays
      integer n, mn
      real s
      dimension n(npblk), mn(lvect), s(npblk,lvect)
!dir$ attributes align: 64:: n, mn, s
      mxy1 = mx1*my1
      lxv = mx + 1
      lxyv = lxv*(my + 1)
      nxyv = nxv*nyv
      mn(1) = 0
      mn(2) = 1
      mn(3) = lxv
      mn(4) = lxv + 1
      mn(5) = lxyv
      mn(6) = lxyv + 1
      mn(7) = lxyv + lxv
      mn(8) = lxyv + lxv + 1
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV).or.(mz.ge.MZV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,m,noff,moff,loff,npp,ipp,joff,nps,nn,mm,ll,nm,lm,
!$OMP& x,y,z,w,dxp,dyp,dzp,amx,amy,amz,dx1,sq,n,s)
      do 170 l = 1, mxyz1
      loff = (l - 1)/mxy1
      k = l - mxy1*loff
      loff = mz*loff
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(l)
c zero out local accumulator
      do 10 j = 1, (mx+1)*(my+1)*(mz+1)
      sq(j) = 0.0
   10 continue
c loop over particles in tile
      ipp = npp/npblk
c outer loop over number of full blocks
      do 50 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
!dir$ vector aligned
      do 20 j = 1, npblk
c find interpolation weights
      x = ppart(j+joff,1,l)
      y = ppart(j+joff,2,l)
      z = ppart(j+joff,3,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn - noff + lxv*(mm - moff) + lxyv*(ll - loff) + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
   20 continue
c deposit charge within tile to local accumulator
      do 40 j = 1, npblk
!dir$ ivdep
      do 30 i = 1, lvect
      sq(n(j)+mn(i)) = sq(n(j)+mn(i)) + s(j,i)
   30 continue
   40 continue
   50 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 60 j = nps, npp
c find interpolation weights
      x = ppart(j,1,l)
      y = ppart(j,2,l)
      z = ppart(j,3,l)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn - noff + 1 + lxv*(mm - moff) + lxyv*(ll - loff)
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c deposit charge within tile to local accumulator
      x = sq(nn) + amx*amz
      y = sq(nn+1) + amy*amz
      z = sq(nn+lxv) + dyp*amz
      w = sq(nn+1+lxv) + dx1*amz
      sq(nn) = x
      sq(nn+1) = y
      sq(nn+lxv) = z
      sq(nn+1+lxv) = w
      mm = nn + lxyv
      x = sq(mm) + amx*dzp
      y = sq(mm+1) + amy*dzp
      z = sq(mm+lxv) + dyp*dzp
      w = sq(mm+1+lxv) + dx1*dzp
      sq(mm) = x
      sq(mm+1) = y
      sq(mm+lxv) = z
      sq(mm+1+lxv) = w
   60 continue
c deposit charge to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      ll = min(mz,nzv-loff)
      do 90 k = 2, ll
      do 80 j = 2, mm
!dir$ ivdep
      do 70 i = 2, nn
      q(i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1)) =                        
     1q(i+noff+nxv*(j+moff-1)+nxyv*(k+loff-1)) +                        
     2 sq(i+lxv*(j-1)+lxyv*(k-1))
   70 continue
   80 continue
   90 continue
c deposit charge to edge points in global array
      lm = min(mz+1,nzv-loff)
      do 110 j = 2, mm
      do 100 i = 2, nn
!$OMP ATOMIC
      q(i+noff+nxv*(j+moff-1)+nxyv*loff) =                              
     1q(i+noff+nxv*(j+moff-1)+nxyv*loff) + sq(i+lxv*(j-1))
      if (lm > mz) then
!$OMP ATOMIC
         q(i+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1)) =                    
     1   q(i+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1))
     1   + sq(i+lxv*(j-1)+lxyv*(lm-1))
      endif
  100 continue
  110 continue
      nm = min(mx+1,nxv-noff)
      mm = min(my+1,nyv-moff)
      do 140 k = 1, ll
      do 120 i = 2, nn
!$OMP ATOMIC
      q(i+noff+nxv*moff+nxyv*(k+loff-1)) =                              
     1q(i+noff+nxv*moff+nxyv*(k+loff-1)) + sq(i+lxyv*(k-1))
      if (mm > my) then
!$OMP ATOMIC
         q(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff-1)) =      
     1   q(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff-1))                      
     2   + sq(i+lxv*(mm-1)+lxyv*(k-1))
      endif
  120 continue
      do 130 j = 1, mm
!$OMP ATOMIC
      q(1+noff+nxv*(j+moff-1)+nxyv*(k+loff-1)) =                        
     1q(1+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))                          
     2+ sq(1+lxv*(j-1)+lxyv*(k-1))
      if (nm > mx) then
!$OMP ATOMIC
         q(nm+noff+nxv*(j+moff-1)+nxyv*(k+loff-1)) =                    
     1   q(nm+noff+nxv*(j+moff-1)+nxyv*(k+loff-1))                      
     1   + sq(nm+lxv*(j-1)+lxyv*(k-1))
      endif
  130 continue
  140 continue
      if (lm > mz) then
         do 150 i = 2, nn
!$OMP ATOMIC
         q(i+noff+nxv*moff+nxyv*(lm+loff-1)) =                          
     1   q(i+noff+nxv*moff+nxyv*(lm+loff-1)) + sq(i+lxyv*(lm-1))
         if (mm > my) then
!$OMP ATOMIC
            q(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1)) =                
     1      q(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))                  
     1      + sq(i+lxv*(mm-1)+lxyv*(lm-1))
         endif
  150    continue
         do 160 j = 1, mm
!$OMP ATOMIC
         q(1+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1)) =                    
     1   q(1+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1))                      
     1   + sq(1+lxv*(j-1)+lxyv*(lm-1))
         if (nm > mx) then
!$OMP ATOMIC
            q(nm+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1)) =                
     1      q(nm+noff+nxv*(j+moff-1)+nxyv*(lm+loff-1))                  
     1      + sq(nm+lxv*(j-1)+lxyv*(lm-1))
         endif
  160    continue
      endif
  170 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine viscan2(isdata,mb,nths)
c performs vectorizable prefix reduction of integer data
c using binary tree method.
      implicit none
      integer nths
      integer isdata, mb
      dimension isdata(nths), mb(nths/2)
c local data
      integer j, kxs, lb, ns
      ns = nths/2
      do 10 j = 1, ns
      mb(j) = j - 1
   10 continue
      kxs = 1
   20 if (kxs.lt.nths) then
!dir$ ivdep
         do 30 j = 1, ns
         lb = kxs*mb(j)
         if ((j+lb+kxs).le.nths) then
            isdata(j+lb+kxs) = isdata(j+lb+kxs) + isdata(2*lb+kxs)
         endif
         mb(j) = mb(j)/2
   30    continue
         kxs = kxs + kxs
         go to 20
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPORDER3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx, 
     1ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
c this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 3D linear memory
c algorithm has 3 steps.  first, one finds particles leaving tile and
c stores their number in each directon, location, and destination in ncl
c and ihole.  second, a prefix scan of ncl is performed and departing
c particles are buffered in ppbuff in direction order.  finally, we copy
c the incoming particles from other tiles into ppart.
c input: all except ppbuff, ncl, ihole, irc
c output: ppart, ppbuff, kpic, ncl, ihole, irc
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppbuff(n,i,l) = i co-ordinate of particle n in tile l
c kpic(l) = number of particles in tile l
c ncl(i,l) = number of particles going to destination i, tile l
c ihole(1,:,l) = location of hole in array left by departing particle
c ihole(2,:,l) = direction destination of particle leaving hole
c all for tile l
c ihole(1,1,l) = ih, number of holes left (error, if negative)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = system length in x/y/z direction
c mx/my/mz = number of grids in sorting cell in x/y/z
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mz1 = (system length in z direction - 1)/mz + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, my1, mz1
      integer npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mx1*my1*mz1)
      dimension ppbuff(npbmx,idimp,mx1*my1*mz1)
      dimension kpic(mx1*my1*mz1), ncl(26,mx1*my1*mz1)
      dimension ihole(2,ntmax+1,mx1*my1*mz1)
c local data
      integer mxy1, mxyz1, noff, moff, loff, npp, ncoff
      integer i, j, k, l, ii, kx, ky, kz, ih, nh, ist, nn, mm, ll, isum
      integer ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dx, dy, dz
      integer ks
      dimension ks(26)
      mxy1 = mx1*my1
      mxyz1 = mxy1*mz1
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
c find and count particles leaving tiles and determine destination
c update ppart, ihole, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,l,noff,moff,loff,npp,nn,mm,ll,ih,nh,ist,dx,dy,dz,    
!$OMP& edgelx,edgely,edgelz,edgerx,edgery,edgerz)
      do 30 l = 1, mxyz1
      loff = (l - 1)/mxy1
      k = l - mxy1*loff
      loff = mz*loff
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(l)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      ll = min(mz,nz-loff)
      ih = 0
      nh = 0
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
      edgelz = loff
      edgerz = loff + ll
c clear counters
      do 10 j = 1, 26
      ncl(j,l) = 0
   10 continue
c loop over particles in tile
      do 20 j = 1, npp
      dx = ppart(j,1,l)
      dy = ppart(j,2,l)
      dz = ppart(j,3,l)
c find particles going out of bounds
      ist = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(j,1,l) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(j,1,l) = dx
         else
            ist = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(j,2,l) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(j,2,l) = dy
         else
            ist = ist + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) ppart(j,3,l) = dz - anz
         ist = ist + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               ist = ist + 9
            else
               dz = 0.0
            endif
            ppart(j,3,l) = dz
         else
            ist = ist + 9
         endif
      endif
      if (ist.gt.0) then
         ncl(ist,l) = ncl(ist,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j
            ihole(2,ih+1,l) = ist
         else
            nh = 1
         endif
      endif
   20 continue
c set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
   30 continue
!$OMP END PARALLEL DO
c ihole overflow
      if (irc.gt.0) return
c
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,isum,ist,nh,ip,j1,ii)
      do 70 l = 1, mxyz1
c find address offset for ordered ppbuff array
      isum = 0
      do 40 j = 1, 26
      ist = ncl(j,l)
      ncl(j,l) = isum
      isum = isum + ist
   40 continue
      nh = ihole(1,1,l)
      ip = 0
c loop over particles leaving tile
      do 60 j = 1, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,l)
      ist = ihole(2,j+1,l)
      ii = ncl(ist,l) + 1
      if (ii.le.npbmx) then
         do 50 i = 1, idimp
         ppbuff(ii,i,l) = ppart(j1,i,l)
   50    continue
      else
         ip = 1
      endif
      ncl(ist,l) = ii
   60 continue
c set error
      if (ip.gt.0) irc = ncl(26,l)
   70 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ii,kk,npp,kx,ky,kz,kl,kr,kxl,kxr,lk,ll,lr,ih,nh, 
!$OMP& nn,ncoff,ist,j1,j2,ip,ks)
      do 130 l = 1, mxyz1
      npp = kpic(l)
      kz = (l - 1)/mxy1
      k = l - mxy1*kz
      kz = kz + 1
c loop over tiles in z, assume periodic boundary conditions
      lk = (kz - 1)*mxy1
c find tile behind
      ll = kz - 1 
      if (ll.lt.1) ll = ll + mz1
      ll = (ll - 1)*mxy1
c find tile in front
      lr = kz + 1
      if (lr.gt.mz1) lr = lr - mz1
      lr = (lr - 1)*mxy1
      ky = (k - 1)/mx1 + 1
c loop over tiles in y, assume periodic boundary conditions
      kk = (ky - 1)*mx1
c find tile above
      kl = ky - 1 
      if (kl.lt.1) kl = kl + my1
      kl = (kl - 1)*mx1
c find tile below
      kr = ky + 1
      if (kr.gt.my1) kr = kr - my1
      kr = (kr - 1)*mx1
c loop over tiles in x, assume periodic boundary conditions
      kx = k - (ky - 1)*mx1
      kxl = kx - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = kx + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
c find tile number for different directions
      ks(1) = kxr + kk + lk
      ks(2) = kxl + kk + lk
      ks(3) = kx + kr + lk
      ks(4) = kxr + kr + lk
      ks(5) = kxl + kr + lk
      ks(6) = kx + kl + lk
      ks(7) = kxr + kl + lk
      ks(8) = kxl + kl + lk
      ks(9) = kx + kk + lr
      ks(10) = kxr + kk + lr
      ks(11) = kxl + kk + lr
      ks(12) = kx + kr + lr
      ks(13) = kxr + kr + lr
      ks(14) = kxl + kr + lr
      ks(15) = kx + kl + lr
      ks(16) = kxr + kl + lr
      ks(17) = kxl + kl + lr
      ks(18) = kx + kk + ll
      ks(19) = kxr + kk + ll
      ks(20) = kxl + kk + ll
      ks(21) = kx + kr + ll
      ks(22) = kxr + kr + ll
      ks(23) = kxl + kr + ll
      ks(24) = kx + kl + ll
      ks(25) = kxr + kl + ll
      ks(26) = kxl + kl + ll
c loop over directions
      nh = ihole(1,1,l)
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 100 ii = 1, 26
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
c ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
      do 90 j = 1, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,l)
c place overflow at end of array
      else
         j1 = npp + 1
         npp = j1
      endif
      if (j1.le.nppmx) then
         do 80 i = 1, idimp
         ppart(j1,i,l) = ppbuff(j+ncoff,i,ks(ii))
   80    continue
      else
         ist = 1
      endif
   90 continue
  100 continue
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
c holes with locations great than npp-ip do not need to be filled
      if (ih.lt.nh) then
         ip = nh - ih
         ii = nh + 1
         nn = ihole(1,ii,l)
         ih = ih + 2
         j2 = ihole(1,ih,l)
c move particles from end into remaining holes
c holes are processed in increasing order
         do 120 j = 1, ip
         j1 = npp - j + 1
         if (j1.eq.nn) then
            ii = ii - 1
            nn = ihole(1,ii,l)
         else
            do 110 i = 1, idimp
            ppart(j2,i,l) = ppart(j1,i,l)
  110       continue
            ih = ih + 1
            j2 = ihole(1,ih,l)
         endif
  120    continue
         npp = npp - ip
      endif
      kpic(l) = npp
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPORDERF3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,mx1
     1,my1,mz1,npbmx,ntmax,irc)
c this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 3D linear memory
c the algorithm has 2 steps.  first, a prefix scan of ncl is performed
c and departing particles are buffered in ppbuff in direction order.
c then we copy the incoming particles from other tiles into ppart.
c it assumes that the number, location, and destination of particles 
c leaving a tile have been previously stored in ncl and ihole by the
c GPPUSHF3LT subroutine.
c input: all except ppbuff, irc
c output: ppart, ppbuff, kpic, ncl, irc
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppbuff(n,i,l) = i co-ordinate of particle n in tile l
c kpic(l) = number of particles in tile l
c ncl(i,l) = number of particles going to destination i, tile l
c ihole(1,:,l) = location of hole in array left by departing particle
c ihole(2,:,l) = direction destination of particle leaving hole
c all for tile l
c ihole(1,1,l) = ih, number of holes left (error, if negative)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mz1 = (system length in z direction - 1)/mz + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, my1, mz1, npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mx1*my1*mz1)
      dimension ppbuff(npbmx,idimp,mx1*my1*mz1)
      dimension kpic(mx1*my1*mz1), ncl(26,mx1*my1*mz1)
      dimension ihole(2,ntmax+1,mx1*my1*mz1)
c local data
      integer mxy1, mxyz1, npp, ncoff
      integer i, j, k, l, ii, kx, ky, kz, ih, nh, ist, nn, ll, isum
      integer ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr
      integer ks
      dimension ks(26)
      mxy1 = mx1*my1
      mxyz1 = mxy1*mz1
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,isum,ist,nh,ip,j1,ii)
      do 40 l = 1, mxyz1
c find address offset for ordered ppbuff array
      isum = 0
      do 10 j = 1, 26
      ist = ncl(j,l)
      ncl(j,l) = isum
      isum = isum + ist
   10 continue
      nh = ihole(1,1,l)
      ip = 0
c loop over particles leaving tile
      do 30 j = 1, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,l)
      ist = ihole(2,j+1,l)
      ii = ncl(ist,l) + 1
      if (ii.le.npbmx) then
         do 20 i = 1, idimp
         ppbuff(ii,i,l) = ppart(j1,i,l)
   20    continue
      else
         ip = 1
      endif
      ncl(ist,l) = ii
   30 continue
c set error
      if (ip.gt.0) irc = ncl(26,l)
   40 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ii,kk,npp,kx,ky,kz,kl,kr,kxl,kxr,lk,ll,lr,ih,nh, 
!$OMP& nn,ncoff,ist,j1,j2,ip,ks)
      do 100 l = 1, mxyz1
      npp = kpic(l)
      kz = (l - 1)/mxy1
      k = l - mxy1*kz
      kz = kz + 1
c loop over tiles in z, assume periodic boundary conditions
      lk = (kz - 1)*mxy1
c find tile behind
      ll = kz - 1 
      if (ll.lt.1) ll = ll + mz1
      ll = (ll - 1)*mxy1
c find tile in front
      lr = kz + 1
      if (lr.gt.mz1) lr = lr - mz1
      lr = (lr - 1)*mxy1
      ky = (k - 1)/mx1 + 1
c loop over tiles in y, assume periodic boundary conditions
      kk = (ky - 1)*mx1
c find tile above
      kl = ky - 1 
      if (kl.lt.1) kl = kl + my1
      kl = (kl - 1)*mx1
c find tile below
      kr = ky + 1
      if (kr.gt.my1) kr = kr - my1
      kr = (kr - 1)*mx1
c loop over tiles in x, assume periodic boundary conditions
      kx = k - (ky - 1)*mx1
      kxl = kx - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = kx + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
c find tile number for different directions
      ks(1) = kxr + kk + lk
      ks(2) = kxl + kk + lk
      ks(3) = kx + kr + lk
      ks(4) = kxr + kr + lk
      ks(5) = kxl + kr + lk
      ks(6) = kx + kl + lk
      ks(7) = kxr + kl + lk
      ks(8) = kxl + kl + lk
      ks(9) = kx + kk + lr
      ks(10) = kxr + kk + lr
      ks(11) = kxl + kk + lr
      ks(12) = kx + kr + lr
      ks(13) = kxr + kr + lr
      ks(14) = kxl + kr + lr
      ks(15) = kx + kl + lr
      ks(16) = kxr + kl + lr
      ks(17) = kxl + kl + lr
      ks(18) = kx + kk + ll
      ks(19) = kxr + kk + ll
      ks(20) = kxl + kk + ll
      ks(21) = kx + kr + ll
      ks(22) = kxr + kr + ll
      ks(23) = kxl + kr + ll
      ks(24) = kx + kl + ll
      ks(25) = kxr + kl + ll
      ks(26) = kxl + kl + ll
c loop over directions
      nh = ihole(1,1,l)
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 70 ii = 1, 26
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
c ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
      do 60 j = 1, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,l)
c place overflow at end of array
      else
         j1 = npp + 1
         npp = j1
      endif
      if (j1.le.nppmx) then
         do 50 i = 1, idimp
         ppart(j1,i,l) = ppbuff(j+ncoff,i,ks(ii))
   50    continue
      else
         ist = 1
      endif
   60 continue
   70 continue
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
c holes with locations great than npp-ip do not need to be filled
      if (ih.lt.nh) then
         ip = nh - ih
         ii = nh + 1
         nn = ihole(1,ii,l)
         ih = ih + 2
         j2 = ihole(1,ih,l)
c move particles from end into remaining holes
c holes are processed in increasing order
         do 90 j = 1, ip
         j1 = npp - j + 1
         if (j1.eq.nn) then
            ii = ii - 1
            nn = ihole(1,ii,l)
         else
            do 80 i = 1, idimp
            ppart(j2,i,l) = ppart(j1,i,l)
   80       continue
            ih = ih + 1
            j2 = ihole(1,ih,l)
         endif
   90    continue
         npp = npp - ip
      endif
      kpic(l) = npp
  100 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine VPPORDER3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx,
     1ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,irc)
c this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 3D linear memory
c algorithm has 3 steps.  first, one finds particles leaving tile and
c stores their number in each directon, location, and destination in ncl
c and ihole.  second, a prefix scan of ncl is performed and departing
c particles are buffered in ppbuff in direction order.  finally, we copy
c the incoming particles from other tiles into ppart.
c input: all except ppbuff, ncl, ihole, irc
c output: ppart, ppbuff, kpic, ncl, ihole, irc
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppbuff(n,i,l) = i co-ordinate of particle n in tile l
c kpic(l) = number of particles in tile l
c ncl(i,l) = number of particles going to destination i, tile l
c ihole(1,:,l) = location of hole in array left by departing particle
c ihole(2,:,l) = direction destination of particle leaving hole
c all for tile l
c ihole(1,1,l) = ih, number of holes left (error, if negative)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c nx/ny/nz = system length in x/y/z direction
c mx/my/mz = number of grids in sorting cell in x/y/z
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mz1 = (system length in z direction - 1)/mz + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, ny, nz, mx, my, mz, mx1, my1, mz1
      integer npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mx1*my1*mz1)
      dimension ppbuff(npbmx,idimp,mx1*my1*mz1)
      dimension kpic(mx1*my1*mz1), ncl(26,mx1*my1*mz1)
      dimension ihole(2,ntmax+1,mx1*my1*mz1)
c local data
      integer npblk
      parameter(npblk=16)
      integer mxy1, mxyz1, noff, moff, loff, npp, ipp, joff, nps, ncoff
      integer i, j, k, l, m, ii, kx, ky, kz, ih, nh, ist, nn, mm, ll, in
      integer ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr, lb, kxs
      real anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dx, dy, dz
      integer sncl, ks
      dimension sncl(26), ks(26)
!dir$ attributes align: 64:: sncl, ks
c scratch arrays
      integer n
      dimension n(npblk,3)
!dir$ attributes align: 64:: n
      mxy1 = mx1*my1
      mxyz1 = mxy1*mz1
      anx = real(nx)
      any = real(ny)
      anz = real(nz)
c find and count particles leaving tiles and determine destination
c update ppart, ihole, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,l,noff,moff,loff,npp,ipp,joff,nps,nn,mm,ll,ih,nh,ist,
!$OMP& dx,dy,dz,edgelx,edgely,edgelz,edgerx,edgery,edgerz,n)
      do 60 l = 1, mxyz1
      loff = (l - 1)/mxy1
      k = l - mxy1*loff
      loff = mz*loff
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(l)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      ll = min(mz,nz-loff)
      ih = 0
      nh = 0
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
      edgelz = loff
      edgerz = loff + ll
c clear counters
      do 10 j = 1, 26
      ncl(j,l) = 0
   10 continue
c loop over particles in tile
      ipp = npp/npblk
c outer loop over number of full blocks
      do 40 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
!dir$ vector aligned
      do 20 j = 1, npblk
      dx = ppart(j+joff,1,l)
      dy = ppart(j+joff,2,l)
      dz = ppart(j+joff,3,l)
c find particles going out of bounds
      ist = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(j+joff,1,l) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(j+joff,1,l) = dx
         else
            ist = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(j+joff,2,l) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(j+joff,2,l) = dy
         else
            ist = ist + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) ppart(j+joff,3,l) = dz - anz
         ist = ist + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               ist = ist + 9
            else
               dz = 0.0
            endif
            ppart(j+joff,3,l) = dz
         else
            ist = ist + 9
         endif
      endif
      n(j,1) = ist
   20 continue
! store outgoing particle address and destination
      do 30 j = 1, npblk
      ist = n(j,1)
      if (ist.gt.0) then
         ncl(ist,l) = ncl(ist,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j + joff
            ihole(2,ih+1,l) = ist
         else
            nh = 1
         endif
      endif
   30 continue
   40 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 50 j = nps, npp
      dx = ppart(j,1,l)
      dy = ppart(j,2,l)
      dz = ppart(j,3,l)
c find particles going out of bounds
      ist = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(j,1,l) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(j,1,l) = dx
         else
            ist = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(j,2,l) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(j,2,l) = dy
         else
            ist = ist + 3
         endif
      endif
      if (dz.ge.edgerz) then
         if (dz.ge.anz) ppart(j,3,l) = dz - anz
         ist = ist + 18
      else if (dz.lt.edgelz) then
         if (dz.lt.0.0) then
            dz = dz + anz
            if (dz.lt.anz) then
               ist = ist + 9
            else
               dz = 0.0
            endif
            ppart(j,3,l) = dz
         else
            ist = ist + 9
         endif
      endif
      if (ist.gt.0) then
         ncl(ist,l) = ncl(ist,l) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,l) = j
            ihole(2,ih+1,l) = ist
         else
            nh = 1
         endif
      endif
   50 continue
c set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,l) = ih
   60 continue
!$OMP END PARALLEL DO
c ihole overflow
      if (irc.gt.0) return
c
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,m,kxs,lb,ist,nh,ip,ipp,nps,joff,j1,ii,sncl,ks,n)
      do 190 l = 1, mxyz1
c find address offset for ordered ppbuff array
      do 70 j = 1, 26
      sncl(j) = ncl(j,l)
      ks(j) = j - 1
   70 continue
      kxs = 1
   80 if (kxs.lt.26) then
!dir$ ivdep
         do 90 j = 1, 13
         lb = kxs*ks(j)
         if ((j+lb+kxs).le.26) then
            sncl(j+lb+kxs) = sncl(j+lb+kxs) + sncl(2*lb+kxs)
         endif
         ks(j) = ks(j)/2
   90    continue
         kxs = kxs + kxs
         go to 80
      endif
      do 100 j = 1, 26
      sncl(j) = sncl(j) - ncl(j,l)
  100 continue
      nh = ihole(1,1,l)
      ip = 0
c buffer particles that are leaving tile, in direction order
c loop over particles leaving tile
      ipp = nh/npblk
c outer loop over number of full blocks
      do 150 m = 1, ipp
      joff = npblk*(m - 1) + 1
c inner loop over particles in block
      do 110 j = 1, npblk
      n(j,1) = ihole(1,j+joff,l)
      n(j,2) = ihole(2,j+joff,l)
  110 continue
c calculate offsets
      do 120 j = 1, npblk
      ist = n(j,2)
      ii = sncl(ist) + 1
      n(j,2) = ii
      sncl(ist) = ii
  120 continue
c buffer particles that are leaving tile, in direction order
      do 140 i = 1, idimp
      do 130 j = 1, npblk
      j1 = n(j,1)
      ii = n(j,2)
      if (ii.le.npbmx) then
         ppbuff(ii,i,l) = ppart(j1,i,l)
      else
         ip = 1
      endif
  130 continue
  140 continue
  150 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 170 j = nps, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,l)
      ist = ihole(2,j+1,l)
      ii = sncl(ist) + 1
      if (ii.le.npbmx) then
         do 160 i = 1, idimp
         ppbuff(ii,i,l) = ppart(j1,i,l)
  160    continue
      else
         ip = 1
      endif
      sncl(ist) = ii
  170 continue
      do 180 j = 1, 26
      ncl(j,l) = sncl(j)
  180 continue
c set error
      if (ip.gt.0) irc = ncl(26,l)
  190 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,m,ii,kk,in,npp,ipp,joff,nps,kx,ky,kz,kl,kr,kxl,  
!$OMP& kxr,lk,ll,lr,ih,nh,nn,mm,ncoff,ist,j1,j2,ip,ks,n)
      do 340 l = 1, mxyz1
      npp = kpic(l)
      kz = (l - 1)/mxy1
      k = l - mxy1*kz
      kz = kz + 1
c loop over tiles in z, assume periodic boundary conditions
      lk = (kz - 1)*mxy1
c find tile behind
      ll = kz - 1 
      if (ll.lt.1) ll = ll + mz1
      ll = (ll - 1)*mxy1
c find tile in front
      lr = kz + 1
      if (lr.gt.mz1) lr = lr - mz1
      lr = (lr - 1)*mxy1
      ky = (k - 1)/mx1 + 1
c loop over tiles in y, assume periodic boundary conditions
      kk = (ky - 1)*mx1
c find tile above
      kl = ky - 1 
      if (kl.lt.1) kl = kl + my1
      kl = (kl - 1)*mx1
c find tile below
      kr = ky + 1
      if (kr.gt.my1) kr = kr - my1
      kr = (kr - 1)*mx1
c loop over tiles in x, assume periodic boundary conditions
      kx = k - (ky - 1)*mx1
      kxl = kx - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = kx + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
c find tile number for different directions
      ks(1) = kxr + kk + lk
      ks(2) = kxl + kk + lk
      ks(3) = kx + kr + lk
      ks(4) = kxr + kr + lk
      ks(5) = kxl + kr + lk
      ks(6) = kx + kl + lk
      ks(7) = kxr + kl + lk
      ks(8) = kxl + kl + lk
      ks(9) = kx + kk + lr
      ks(10) = kxr + kk + lr
      ks(11) = kxl + kk + lr
      ks(12) = kx + kr + lr
      ks(13) = kxr + kr + lr
      ks(14) = kxl + kr + lr
      ks(15) = kx + kl + lr
      ks(16) = kxr + kl + lr
      ks(17) = kxl + kl + lr
      ks(18) = kx + kk + ll
      ks(19) = kxr + kk + ll
      ks(20) = kxl + kk + ll
      ks(21) = kx + kr + ll
      ks(22) = kxr + kr + ll
      ks(23) = kxl + kr + ll
      ks(24) = kx + kl + ll
      ks(25) = kxr + kl + ll
      ks(26) = kxl + kl + ll
c loop over directions
      nh = ihole(1,1,l)
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 260 ii = 1, 26
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
c ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
c loop over particles coming from direction ii
      ipp = ip/npblk
c outer loop over number of full blocks
      do 230 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 200 j = 1, npblk
c insert incoming particles into holes
      if ((j+ih).le.nh) then
         j1 = ihole(1,j+ih+1,l)
c place overflow at end of array
      else
         j1 = npp + j + ih - nh
      endif
      n(j,1) = j1
  200 continue
      do 220 i = 1, idimp
      do 210 j = 1, npblk
      j1 = n(j,1)
      if (j1.le.nppmx) then
         ppart(j1,i,l) = ppbuff(j+joff+ncoff,i,ks(ii))
      else
         ist = 1
      endif
  210 continue
  220 continue
      ih = ih + npblk
  230 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 250 j = nps, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,l)
c place overflow at end of array
      else
         j1 = npp + ih - nh
      endif
      if (j1.le.nppmx) then
         do 240 i = 1, idimp
         ppart(j1,i,l) = ppbuff(j+ncoff,i,ks(ii))
  240    continue
      else
         ist = 1
      endif
  250 continue
  260 continue
      if (ih > nh) npp = npp + ih - nh
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
c holes with locations great than npp-ip do not need to be filled
      if (ih.lt.nh) then
         ip = nh - ih
c move particles from end into remaining holes
c holes are processed in increasing order
         ii = nh + 1
         ipp = ip/npblk
c outer loop over number of full blocks
         do 310 m = 1, ipp
         joff = npblk*(m - 1)
c inner loop over particles in block
         do 270 j = 1, npblk
         n(j,2) = ihole(1,ih+j+1,l)
         n(j,3) = ihole(1,ii-j+1,l)
  270    continue
         in = 1
         mm = 1
         nn = n(in,3)
         do 280 j = 1, npblk
         j1 = npp - j - joff + 1
         n(j,1) = n(mm,2)
         if (j1.eq.nn) then
            in = in + 1
            nn = n(in,3)
            n(j,1) = -1
         else
            mm = mm + 1
         endif
  280    continue
         do 300 i = 1, idimp
!dir$ ivdep
         do 290 j = 1, npblk
         j1 = npp - j - joff + 1
         j2 = n(j,1)
         if (j2.gt.0) ppart(j2,i,l) = ppart(j1,i,l)
  290    continue
  300    continue
         ii = ii - in + 1
         ih = ih + mm - 1
  310    continue
         nps = npblk*ipp + 1
         nn = ihole(1,ii,l)
         ih = ih + 2
         j2 = ihole(1,ih,l)
c loop over remaining particles
         do 330 j = nps, ip
         j1 = npp - j + 1
         if (j1.eq.nn) then
            ii = ii - 1
            nn = ihole(1,ii,l)
         else
            do 320 i = 1, idimp
            ppart(j2,i,l) = ppart(j1,i,l)
  320       continue
            ih = ih + 1
            j2 = ihole(1,ih,l)
         endif
  330    continue
         npp = npp - ip
      endif
      kpic(l) = npp
  340 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine VPPORDERF3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,  
     1mx1,my1,mz1,npbmx,ntmax,irc)
c this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 3D linear memory
c the algorithm has 2 steps.  first, a prefix scan of ncl is performed
c and departing particles are buffered in ppbuff in direction order.
c then we copy the incoming particles from other tiles into ppart.
c it assumes that the number, location, and destination of particles 
c leaving a tile have been previously stored in ncl and ihole by the
c GPPUSHF3LT subroutine.
c input: all except ppbuff, irc
c output: ppart, ppbuff, kpic, ncl, irc
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppbuff(n,i,l) = i co-ordinate of particle n in tile l
c kpic(l) = number of particles in tile l
c ncl(i,l) = number of particles going to destination i, tile l
c ihole(1,:,l) = location of hole in array left by departing particle
c ihole(2,:,l) = direction destination of particle leaving hole
c all for tile l
c ihole(1,1,l) = ih, number of holes left (error, if negative)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mz1 = (system length in z direction - 1)/mz + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, my1, mz1, npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mx1*my1*mz1)
      dimension ppbuff(npbmx,idimp,mx1*my1*mz1)
      dimension kpic(mx1*my1*mz1), ncl(26,mx1*my1*mz1)
      dimension ihole(2,ntmax+1,mx1*my1*mz1)
c local data
      integer npblk
      parameter(npblk=16)
      integer mxy1, mxyz1, npp, ncoff
      integer i, j, k, l, ii, kx, ky, kz, ih, nh, ist, nn, ll, mm, in
      integer ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr
      integer lb, kxs, m, ipp, nps, joff
      integer sncl, ks
      dimension sncl(26), ks(26)
!dir$ attributes align: 64:: sncl, ks
c scratch arrays
      integer n
      dimension n(npblk,3)
!dir$ attributes align: 64:: n
      mxy1 = mx1*my1
      mxyz1 = mxy1*mz1
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,m,kxs,lb,ist,nh,ip,ipp,nps,joff,j1,ii,sncl,ks,n)
      do 130 l = 1, mxyz1
c find address offset for ordered ppbuff array
      do 10 j = 1, 26
      sncl(j) = ncl(j,l)
      ks(j) = j - 1
   10 continue
      kxs = 1
   20 if (kxs.lt.26) then
!dir$ ivdep
         do 30 j = 1, 13
         lb = kxs*ks(j)
         if ((j+lb+kxs).le.26) then
            sncl(j+lb+kxs) = sncl(j+lb+kxs) + sncl(2*lb+kxs)
         endif
         ks(j) = ks(j)/2
   30    continue
         kxs = kxs + kxs
         go to 20
      endif
      do 40 j = 1, 26
      sncl(j) = sncl(j) - ncl(j,l)
   40 continue
      nh = ihole(1,1,l)
      ip = 0
c buffer particles that are leaving tile, in direction order
c loop over particles leaving tile
      ipp = nh/npblk
c outer loop over number of full blocks
      do 90 m = 1, ipp
      joff = npblk*(m - 1) + 1
c inner loop over particles in block
      do 50 j = 1, npblk
      n(j,1) = ihole(1,j+joff,l)
      n(j,2) = ihole(2,j+joff,l)
   50 continue
c calculate offsets
      do 60 j = 1, npblk
      ist = n(j,2)
      ii = sncl(ist) + 1
      n(j,2) = ii
      sncl(ist) = ii
   60 continue
c buffer particles that are leaving tile, in direction order
      do 80 i = 1, idimp
      do 70 j = 1, npblk
      j1 = n(j,1)
      ii = n(j,2)
      if (ii.le.npbmx) then
         ppbuff(ii,i,l) = ppart(j1,i,l)
      else
         ip = 1
      endif
   70 continue
   80 continue
   90 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 110 j = nps, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,l)
      ist = ihole(2,j+1,l)
      ii = sncl(ist) + 1
      if (ii.le.npbmx) then
         do 100 i = 1, idimp
         ppbuff(ii,i,l) = ppart(j1,i,l)
  100    continue
      else
         ip = 1
      endif
      sncl(ist) = ii
  110 continue
      do 120 j = 1, 26
      ncl(j,l) = sncl(j)
  120 continue
c set error
      if (ip.gt.0) irc = ncl(26,l)
  130 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,m,ii,kk,in,npp,ipp,joff,nps,kx,ky,kz,kl,kr,kxl,  
!$OMP& kxr,lk,ll,lr,ih,nh,nn,mm,ncoff,ist,j1,j2,ip,ks,n)
      do 280 l = 1, mxyz1
      npp = kpic(l)
      kz = (l - 1)/mxy1
      k = l - mxy1*kz
      kz = kz + 1
c loop over tiles in z, assume periodic boundary conditions
      lk = (kz - 1)*mxy1
c find tile behind
      ll = kz - 1 
      if (ll.lt.1) ll = ll + mz1
      ll = (ll - 1)*mxy1
c find tile in front
      lr = kz + 1
      if (lr.gt.mz1) lr = lr - mz1
      lr = (lr - 1)*mxy1
      ky = (k - 1)/mx1 + 1
c loop over tiles in y, assume periodic boundary conditions
      kk = (ky - 1)*mx1
c find tile above
      kl = ky - 1 
      if (kl.lt.1) kl = kl + my1
      kl = (kl - 1)*mx1
c find tile below
      kr = ky + 1
      if (kr.gt.my1) kr = kr - my1
      kr = (kr - 1)*mx1
c loop over tiles in x, assume periodic boundary conditions
      kx = k - (ky - 1)*mx1
      kxl = kx - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = kx + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
c find tile number for different directions
      ks(1) = kxr + kk + lk
      ks(2) = kxl + kk + lk
      ks(3) = kx + kr + lk
      ks(4) = kxr + kr + lk
      ks(5) = kxl + kr + lk
      ks(6) = kx + kl + lk
      ks(7) = kxr + kl + lk
      ks(8) = kxl + kl + lk
      ks(9) = kx + kk + lr
      ks(10) = kxr + kk + lr
      ks(11) = kxl + kk + lr
      ks(12) = kx + kr + lr
      ks(13) = kxr + kr + lr
      ks(14) = kxl + kr + lr
      ks(15) = kx + kl + lr
      ks(16) = kxr + kl + lr
      ks(17) = kxl + kl + lr
      ks(18) = kx + kk + ll
      ks(19) = kxr + kk + ll
      ks(20) = kxl + kk + ll
      ks(21) = kx + kr + ll
      ks(22) = kxr + kr + ll
      ks(23) = kxl + kr + ll
      ks(24) = kx + kl + ll
      ks(25) = kxr + kl + ll
      ks(26) = kxl + kl + ll
c loop over directions
      nh = ihole(1,1,l)
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 200 ii = 1, 26
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
c ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
c loop over particles coming from direction ii
      ipp = ip/npblk
c outer loop over number of full blocks
      do 170 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 140 j = 1, npblk
c insert incoming particles into holes
      if ((j+ih).le.nh) then
         j1 = ihole(1,j+ih+1,l)
c place overflow at end of array
      else
         j1 = npp + j + ih - nh
      endif
      n(j,1) = j1
  140 continue
      do 160 i = 1, idimp
      do 150 j = 1, npblk
      j1 = n(j,1)
      if (j1.le.nppmx) then
         ppart(j1,i,l) = ppbuff(j+joff+ncoff,i,ks(ii))
      else
         ist = 1
      endif
  150 continue
  160 continue
      ih = ih + npblk
  170 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 190 j = nps, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,l)
c place overflow at end of array
      else
         j1 = npp + ih - nh
      endif
      if (j1.le.nppmx) then
         do 180 i = 1, idimp
         ppart(j1,i,l) = ppbuff(j+ncoff,i,ks(ii))
  180    continue
      else
         ist = 1
      endif
  190 continue
  200 continue
      if (ih > nh) npp = npp + ih - nh
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
c holes with locations great than npp-ip do not need to be filled
      if (ih.lt.nh) then
         ip = nh - ih
c move particles from end into remaining holes
c holes are processed in increasing order
         ii = nh + 1
         ipp = ip/npblk
c outer loop over number of full blocks
         do 250 m = 1, ipp
         joff = npblk*(m - 1)
c inner loop over particles in block
         do 210 j = 1, npblk
         n(j,2) = ihole(1,ih+j+1,l)
         n(j,3) = ihole(1,ii-j+1,l)
  210    continue
         in = 1
         mm = 1
         nn = n(in,3)
         do 220 j = 1, npblk
         j1 = npp - j - joff + 1
         n(j,1) = n(mm,2)
         if (j1.eq.nn) then
            in = in + 1
            nn = n(in,3)
            n(j,1) = -1
         else
            mm = mm + 1
         endif
  220    continue
         do 240 i = 1, idimp
!dir$ ivdep
         do 230 j = 1, npblk
         j1 = npp - j - joff + 1
         j2 = n(j,1)
         if (j2.gt.0) ppart(j2,i,l) = ppart(j1,i,l)
  230    continue
  240    continue
         ii = ii - in + 1
         ih = ih + mm - 1
  250    continue
         nps = npblk*ipp + 1
         nn = ihole(1,ii,l)
         ih = ih + 2
         j2 = ihole(1,ih,l)
c loop over remaining particles
         do 270 j = nps, ip
         j1 = npp - j + 1
         if (j1.eq.nn) then
            ii = ii - 1
            nn = ihole(1,ii,l)
         else
            do 260 i = 1, idimp
            ppart(j2,i,l) = ppart(j1,i,l)
  260       continue
            ih = ih + 1
            j2 = ihole(1,ih,l)
         endif
  270    continue
         npp = npp - ip
      endif
      kpic(l) = npp
  280 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine V2PPORDERF3LT(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, 
     1mx1,my1,mz1,npbmx,ntmax,irc)
c this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 3D linear memory
c the algorithm has 2 steps.  first, a prefix scan of ncl is performed
c and departing particles are buffered in ppbuff in direction order.
c then we copy the incoming particles from other tiles into ppart.
c it assumes that the number, location, and destination of particles 
c leaving a tile have been previously stored in ncl and ihole by the
c GPPUSHF3LT subroutine.
c input: all except ppbuff, irc
c output: ppart, ppbuff, kpic, ncl, irc
c ppart(n,1,m) = position x of particle n in tile m
c ppart(n,2,m) = position y of particle n in tile m
c ppart(n,3,m) = position z of particle n in tile m
c ppbuff(i,n,l) = i co-ordinate of particle n in tile l
c kpic(l) = number of particles in tile l
c ncl(i,l) = number of particles going to destination i, tile l
c ihole(1,:,l) = location of hole in array left by departing particle
c ihole(2,:,l) = direction destination of particle leaving hole
c all for tile l
c ihole(1,1,l) = ih, number of holes left (error, if negative)
c idimp = size of phase space = 6
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mz1 = (system length in z direction - 1)/mz + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, my1, mz1, npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mx1*my1*mz1)
      dimension ppbuff(idimp,npbmx,mx1*my1*mz1)
      dimension kpic(mx1*my1*mz1), ncl(26,mx1*my1*mz1)
      dimension ihole(2,ntmax+1,mx1*my1*mz1)
c local data
      integer npblk
      parameter(npblk=16)
      integer mxy1, mxyz1, npp, ncoff
      integer i, j, k, l, ii, kx, ky, kz, ih, nh, ist, nn, ll, mm, in
      integer ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr
      integer lb, kxs, m, ipp, nps, joff
      integer sncl, ks
      dimension sncl(26), ks(26)
!dir$ attributes align: 64:: sncl, ks
c scratch arrays
      integer n
      dimension n(npblk,3)
!dir$ attributes align: 64:: n
      mxy1 = mx1*my1
      mxyz1 = mxy1*mz1
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,l,m,kxs,lb,ist,nh,ip,ipp,nps,joff,j1,ii,sncl,ks,n)
      do 130 l = 1, mxyz1
c find address offset for ordered ppbuff array
      do 10 j = 1, 26
      sncl(j) = ncl(j,l)
      ks(j) = j - 1
   10 continue
      kxs = 1
   20 if (kxs.lt.26) then
!dir$ ivdep
         do 30 j = 1, 13
         lb = kxs*ks(j)
         if ((j+lb+kxs).le.26) then
            sncl(j+lb+kxs) = sncl(j+lb+kxs) + sncl(2*lb+kxs)
         endif
         ks(j) = ks(j)/2
   30    continue
         kxs = kxs + kxs
         go to 20
      endif
      do 40 j = 1, 26
      sncl(j) = sncl(j) - ncl(j,l)
   40 continue
      nh = ihole(1,1,l)
      ip = 0
c buffer particles that are leaving tile, in direction order
c loop over particles leaving tile
      ipp = nh/npblk
c outer loop over number of full blocks
      do 90 m = 1, ipp
      joff = npblk*(m - 1) + 1
c inner loop over particles in block
      do 50 j = 1, npblk
      n(j,1) = ihole(1,j+joff,l)
      n(j,2) = ihole(2,j+joff,l)
   50 continue
c calculate offsets
      do 60 j = 1, npblk
      ist = n(j,2)
      ii = sncl(ist) + 1
      n(j,2) = ii
      sncl(ist) = ii
   60 continue
c buffer particles that are leaving tile, in direction order
      do 80 j = 1, npblk
      j1 = n(j,1)
      ii = n(j,2)
      if (ii.le.npbmx) then
         do 70 i = 1, idimp
         ppbuff(i,ii,l) = ppart(j1,i,l)
   70    continue
      else
         ip = 1
      endif
   80 continue
   90 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 110 j = nps, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,l)
      ist = ihole(2,j+1,l)
      ii = sncl(ist) + 1
      if (ii.le.npbmx) then
         do 100 i = 1, idimp
         ppbuff(i,ii,l) = ppart(j1,i,l)
  100    continue
      else
         ip = 1
      endif
      sncl(ist) = ii
  110 continue
      do 120 j = 1, 26
      ncl(j,l) = sncl(j)
  120 continue
c set error
      if (ip.gt.0) irc = ncl(26,l)
  130 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,m,ii,kk,in,npp,ipp,joff,nps,kx,ky,kz,kl,kr,kxl,  
!$OMP& kxr,lk,ll,lr,ih,nh,nn,mm,ncoff,ist,j1,j2,ip,ks,n)
      do 280 l = 1, mxyz1
      npp = kpic(l)
      kz = (l - 1)/mxy1
      k = l - mxy1*kz
      kz = kz + 1
c loop over tiles in z, assume periodic boundary conditions
      lk = (kz - 1)*mxy1
c find tile behind
      ll = kz - 1 
      if (ll.lt.1) ll = ll + mz1
      ll = (ll - 1)*mxy1
c find tile in front
      lr = kz + 1
      if (lr.gt.mz1) lr = lr - mz1
      lr = (lr - 1)*mxy1
      ky = (k - 1)/mx1 + 1
c loop over tiles in y, assume periodic boundary conditions
      kk = (ky - 1)*mx1
c find tile above
      kl = ky - 1 
      if (kl.lt.1) kl = kl + my1
      kl = (kl - 1)*mx1
c find tile below
      kr = ky + 1
      if (kr.gt.my1) kr = kr - my1
      kr = (kr - 1)*mx1
c loop over tiles in x, assume periodic boundary conditions
      kx = k - (ky - 1)*mx1
      kxl = kx - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = kx + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
c find tile number for different directions
      ks(1) = kxr + kk + lk
      ks(2) = kxl + kk + lk
      ks(3) = kx + kr + lk
      ks(4) = kxr + kr + lk
      ks(5) = kxl + kr + lk
      ks(6) = kx + kl + lk
      ks(7) = kxr + kl + lk
      ks(8) = kxl + kl + lk
      ks(9) = kx + kk + lr
      ks(10) = kxr + kk + lr
      ks(11) = kxl + kk + lr
      ks(12) = kx + kr + lr
      ks(13) = kxr + kr + lr
      ks(14) = kxl + kr + lr
      ks(15) = kx + kl + lr
      ks(16) = kxr + kl + lr
      ks(17) = kxl + kl + lr
      ks(18) = kx + kk + ll
      ks(19) = kxr + kk + ll
      ks(20) = kxl + kk + ll
      ks(21) = kx + kr + ll
      ks(22) = kxr + kr + ll
      ks(23) = kxl + kr + ll
      ks(24) = kx + kl + ll
      ks(25) = kxr + kl + ll
      ks(26) = kxl + kl + ll
c loop over directions
      nh = ihole(1,1,l)
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 200 ii = 1, 26
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
c ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
c loop over particles coming from direction ii
      ipp = ip/npblk
c outer loop over number of full blocks
      do 170 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 140 j = 1, npblk
c insert incoming particles into holes
      if ((j+ih).le.nh) then
         j1 = ihole(1,j+ih+1,l)
c place overflow at end of array
      else
         j1 = npp + j + ih - nh
      endif
      n(j,1) = j1
  140 continue
      do 160 j = 1, npblk
      j1 = n(j,1)
      if (j1.le.nppmx) then
         do 150 i = 1, idimp
         ppart(j1,i,l) = ppbuff(i,j+joff+ncoff,ks(ii))
  150    continue
      else
         ist = 1
      endif
  160 continue
      ih = ih + npblk
  170 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 190 j = nps, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,l)
c place overflow at end of array
      else
         j1 = npp + ih - nh
      endif
      if (j1.le.nppmx) then
         do 180 i = 1, idimp
         ppart(j1,i,l) = ppbuff(i,j+ncoff,ks(ii))
  180    continue
      else
         ist = 1
      endif
  190 continue
  200 continue
      if (ih > nh) npp = npp + ih - nh
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
c holes with locations great than npp-ip do not need to be filled
      if (ih.lt.nh) then
         ip = nh - ih
c move particles from end into remaining holes
c holes are processed in increasing order
         ii = nh + 1
         ipp = ip/npblk
c outer loop over number of full blocks
         do 250 m = 1, ipp
         joff = npblk*(m - 1)
c inner loop over particles in block
         do 210 j = 1, npblk
         n(j,2) = ihole(1,ih+j+1,l)
         n(j,3) = ihole(1,ii-j+1,l)
  210    continue
         in = 1
         mm = 1
         nn = n(in,3)
         do 220 j = 1, npblk
         j1 = npp - j - joff + 1
         n(j,1) = n(mm,2)
         if (j1.eq.nn) then
            in = in + 1
            nn = n(in,3)
            n(j,1) = -1
         else
            mm = mm + 1
         endif
  220    continue
         do 240 i = 1, idimp
!dir$ ivdep
         do 230 j = 1, npblk
         j1 = npp - j - joff + 1
         j2 = n(j,1)
         if (j2.gt.0) ppart(j2,i,l) = ppart(j1,i,l)
  230    continue
  240    continue
         ii = ii - in + 1
         ih = ih + mm - 1
  250    continue
         nps = npblk*ipp + 1
         nn = ihole(1,ii,l)
         ih = ih + 2
         j2 = ihole(1,ih,l)
c loop over remaining particles
         do 270 j = nps, ip
         j1 = npp - j + 1
         if (j1.eq.nn) then
            ii = ii - 1
            nn = ihole(1,ii,l)
         else
            do 260 i = 1, idimp
            ppart(j2,i,l) = ppart(j1,i,l)
  260       continue
            ih = ih + 1
            j2 = ihole(1,ih,l)
         endif
  270    continue
         npp = npp - ip
      endif
      kpic(l) = npp
  280 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine CGUARD3L(fxyz,nx,ny,nz,nxe,nye,nze)
c replicate extended periodic vector field fxyz
c linear interpolation
c nx/ny/nz = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
c nze = third dimension of field arrays, must be >= nz+1
      implicit none
      real fxyz
      integer nx, ny, nz, nxe, nye, nze
      dimension fxyz(4,nxe,nye,nze)
c local data
      integer j, k, l
c copy edges of extended field
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l)
      do 30 l = 1, nz
      do 10 k = 1, ny
      fxyz(1,nx+1,k,l) = fxyz(1,1,k,l)
      fxyz(2,nx+1,k,l) = fxyz(2,1,k,l)
      fxyz(3,nx+1,k,l) = fxyz(3,1,k,l)
   10 continue
      do 20 j = 1, nx
      fxyz(1,j,ny+1,l) = fxyz(1,j,1,l)
      fxyz(2,j,ny+1,l) = fxyz(2,j,1,l)
      fxyz(3,j,ny+1,l) = fxyz(3,j,1,l)
   20 continue
      fxyz(1,nx+1,ny+1,l) = fxyz(1,1,1,l)
      fxyz(2,nx+1,ny+1,l) = fxyz(2,1,1,l)
      fxyz(3,nx+1,ny+1,l) = fxyz(3,1,1,l)
   30 continue
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(j,k)
      do 50 k = 1, ny
      do 40 j = 1, nx
      fxyz(1,j,k,nz+1) = fxyz(1,j,k,1)
      fxyz(2,j,k,nz+1) = fxyz(2,j,k,1)
      fxyz(3,j,k,nz+1) = fxyz(3,j,k,1)
   40 continue
      fxyz(1,nx+1,k,nz+1) = fxyz(1,1,k,1)
      fxyz(2,nx+1,k,nz+1) = fxyz(2,1,k,1)
      fxyz(3,nx+1,k,nz+1) = fxyz(3,1,k,1)
   50 continue
!$OMP END DO
!$OMP END PARALLEL
      do 60 j = 1, nx
      fxyz(1,j,ny+1,nz+1) = fxyz(1,j,1,1)
      fxyz(2,j,ny+1,nz+1) = fxyz(2,j,1,1)
      fxyz(3,j,ny+1,nz+1) = fxyz(3,j,1,1)
   60 continue
      fxyz(1,nx+1,ny+1,nz+1) = fxyz(1,1,1,1)
      fxyz(2,nx+1,ny+1,nz+1) = fxyz(2,1,1,1)
      fxyz(3,nx+1,ny+1,nz+1) = fxyz(3,1,1,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine AGUARD3L(q,nx,ny,nz,nxe,nye,nze)
c accumulate extended periodic scalar field q
c linear interpolation
c nx/ny/nz = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
c nze = third dimension of field arrays, must be >= nz+1
      implicit none
      real q
      integer nx, ny, nz, nxe, nye, nze
      dimension q(nxe,nye,nze)
c local data
      integer j, k, l
c accumulate edges of extended field
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l)
      do 30 l = 1, nz
      do 10 k = 1, ny
      q(1,k,l) = q(1,k,l) + q(nx+1,k,l)
      q(nx+1,k,l) = 0.0
   10 continue
      do 20 j = 1, nx
      q(j,1,l) = q(j,1,l) + q(j,ny+1,l)
      q(j,ny+1,l) = 0.0
   20 continue
      q(1,1,l) = q(1,1,l) + q(nx+1,ny+1,l)
      q(nx+1,ny+1,l) = 0.0
   30 continue
!$OMP END DO
!$OMP DO PRIVATE(j,k)
      do 50 k = 1, ny
      do 40 j = 1, nx
      q(j,k,1) = q(j,k,1) + q(j,k,nz+1)
      q(j,k,nz+1) = 0.0
   40 continue
      q(1,k,1) = q(1,k,1) + q(nx+1,k,nz+1)
      q(nx+1,k,nz+1) = 0.0
   50 continue
!$OMP END DO
!$OMP END PARALLEL
      do 60 j = 1, nx
      q(j,1,1) = q(j,1,1) + q(j,ny+1,nz+1)
      q(j,ny+1,nz+1) = 0.0
   60 continue
      q(1,1,1) = q(1,1,1) + q(nx+1,ny+1,nz+1)
      q(nx+1,ny+1,nz+1) = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine VMPOIS33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,   
     1nxvh,nyv,nzv,nxhd,nyhd,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.
c for isign = 0, output: ffc
c input: isign,ax,ay,az,affp,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
c for isign = -1, output: fxyz, we
c input: q,ffc,isign,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
c approximate flop count is:
c 59*nxc*nyc*nzc + 26*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c equation used is:
c fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*s(kx,ky,kz),
c fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*s(kx,ky,kz),
c fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
c fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
c fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
c fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
c q(j,k,l) = complex charge density for fourier mode (j-1,k-1,l-1)
c fxyz(1,j,k,l) = x component of complex force/charge
c fxyz(2,j,k,l) = y component of complex force/charge
c fxyz(3,j,k,l) = z component of complex force/charge
c all for fourier mode (j-1,k-1,l-1)
c aimag(ffc(j,k,l)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1,l-1)
c real(ffc(j,k,l)) = potential green's function g
c for fourier mode (j-1,k-1,l-1)
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nzv = third dimension of field arrays, must be >= nz
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
c nzhd = third dimension of form factor array, must be >= nzh
c vectorizable version
      implicit none
      integer isign, nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
      real ax, ay, az, affp, we
      complex q, fxyz, ffc
      dimension q(nxvh,nyv,nzv), fxyz(4,nxvh,nyv,nzv)
      dimension ffc(nxhd,nyhd,nzhd)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real dnx, dny, dnz, dkx, dky, dkz, at1, at2, at3, at4, at5, at6
      complex zero, zt1, zt2
      double precision wp, sum1, sum2
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 40
c prepare form factor array
      do 30 l = 1, nzh
      dkz = dnz*real(l - 1)
      at1 = dkz*dkz
      at2 = (dkz*az)**2
      do 20 k = 1, nyh
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = (dky*ay)**2 + at2
      do 10 j = 1, nxh
      dkx = dnx*real(j - 1)
      at5 = dkx*dkx + at3
      at6 = exp(-.5*((dkx*ax)**2 + at4))
      if (at5.eq.0.) then
         ffc(j,k,l) = cmplx(affp,1.0)
      else
         ffc(j,k,l) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 sum1 = 0.0d0
!$OMP PARALLEL
!$OMP DO PRIVATE(j,k,l,k1,l1,dky,dkz,at1,at2,at3,at4,zt1,zt2,wp)
!$OMP& REDUCTION(+:sum1)
      do 90 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
      wp = 0.0d0
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 60 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
!dir$ ivdep
      do 50 j = 2, nxh
      zt1 = ffc(j,k,l)
      at1 = real(zt1)*aimag(zt1)
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      at4 = dkz*at1
      zt1 = q(j,k,l)
      zt1 = cmplx(aimag(zt1),-real(zt1))
      zt2 = q(j,k1,l)
      zt2 = cmplx(aimag(zt2),-real(zt2))
      fxyz(1,j,k,l) = at2*zt1
      fxyz(2,j,k,l) = at3*zt1
      fxyz(3,j,k,l) = at4*zt1
      fxyz(1,j,k1,l) = at2*zt2
      fxyz(2,j,k1,l) = -at3*zt2
      fxyz(3,j,k1,l) = at4*zt2
      zt1 = q(j,k,l1)
      zt1 = cmplx(aimag(zt1),-real(zt1))
      zt2 = q(j,k1,l1)
      zt2 = cmplx(aimag(zt2),-real(zt2))
      fxyz(1,j,k,l1) = at2*zt1
      fxyz(2,j,k,l1) = at3*zt1
      fxyz(3,j,k,l1) = -at4*zt1
      fxyz(1,j,k1,l1) = at2*zt2
      fxyz(2,j,k1,l1) = -at3*zt2
      fxyz(3,j,k1,l1) = -at4*zt2
      at1 = at1*(q(j,k,l)*conjg(q(j,k,l)) + q(j,k1,l)*conjg(q(j,k1,l))  
     1    + q(j,k,l1)*conjg(q(j,k,l1)) + q(j,k1,l1)*conjg(q(j,k1,l1)))
      wp = wp + dble(at1)
   50 continue
   60 continue
c mode numbers kx = 0, nx/2
!dir$ ivdep
      do 70 k = 2, nyh
      k1 = ny2 - k
      zt1 = ffc(1,k,l)
      at1 = real(zt1)*aimag(zt1)
      at3 = dny*real(k - 1)*at1
      at4 = dkz*at1
      zt1 = q(1,k,l)
      zt1 = cmplx(aimag(zt1),-real(zt1))
      zt2 = q(1,k,l1)
      zt2 = cmplx(aimag(zt2),-real(zt2))
      fxyz(1,1,k,l) = zero
      fxyz(2,1,k,l) = at3*zt1
      fxyz(3,1,k,l) = at4*zt1
      fxyz(1,1,k1,l) = zero
      fxyz(2,1,k1,l) = zero
      fxyz(3,1,k1,l) = zero
      fxyz(1,1,k,l1) = zero
      fxyz(2,1,k,l1) = at3*zt2
      fxyz(3,1,k,l1) = -at4*zt2
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      at1 = at1*(q(1,k,l)*conjg(q(1,k,l)) + q(1,k,l1)*conjg(q(1,k,l1)))
      wp = wp + dble(at1)
   70 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 80 j = 2, nxh
      zt1 = ffc(j,1,l)
      at1 = real(zt1)*aimag(zt1)
      at2 = dnx*real(j - 1)*at1
      at4 = dkz*at1
      zt1 = q(j,1,l)
      zt1 = cmplx(aimag(zt1),-real(zt1))
      zt2 = q(j,1,l1)
      zt2 = cmplx(aimag(zt2),-real(zt2))
      fxyz(1,j,1,l) = at2*zt1
      fxyz(2,j,1,l) = zero
      fxyz(3,j,1,l) = at4*zt1
      fxyz(1,j,k1,l) = zero
      fxyz(2,j,k1,l) = zero
      fxyz(3,j,k1,l) = zero
      fxyz(1,j,1,l1) = at2*zt2
      fxyz(2,j,1,l1) = zero
      fxyz(3,j,1,l1) = -at4*zt2
      fxyz(1,j,k1,l1) = zero
      fxyz(2,j,k1,l1) = zero
      fxyz(3,j,k1,l1) = zero
      at1 = at1*(q(j,1,l)*conjg(q(j,1,l)) + q(j,1,l1)*conjg(q(j,1,l1)))
      wp = wp + dble(at1)
   80 continue
c mode numbers kx = 0, nx/2
      zt1 = ffc(1,1,l)
      at1 = real(zt1)*aimag(zt1)
      at4 = dkz*at1
      fxyz(1,1,1,l) = zero
      fxyz(2,1,1,l) = zero
      zt1 = q(1,1,l)
      fxyz(3,1,1,l) = at4*cmplx(aimag(zt1),-real(zt1))
      fxyz(1,1,k1,l) = zero
      fxyz(2,1,k1,l) = zero
      fxyz(3,1,k1,l) = zero
      fxyz(1,1,1,l1) = zero
      fxyz(2,1,1,l1) = zero
      fxyz(3,1,1,l1) = zero
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      at1 = at1*(q(1,1,l)*conjg(q(1,1,l)))
      wp = wp + dble(at1)
      sum1 = sum1 + wp
   90 continue
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      l1 = nzh + 1
      sum2 = 0.0d0
c mode numbers kz = 0, nz/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,at1,at2,at3,zt1,zt2,wp)
!$OMP& REDUCTION(+:sum2)
      do 110 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      wp = 0.0d0
!dir$ ivdep
      do 100 j = 2, nxh
      zt1 = ffc(j,k,1)
      at1 = real(zt1)*aimag(zt1)
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      zt1 = q(j,k,1)
      zt1 = cmplx(aimag(zt1),-real(zt1))
      zt2 = q(j,k1,1)
      zt2 = cmplx(aimag(zt2),-real(zt2))
      fxyz(1,j,k,1) = at2*zt1
      fxyz(2,j,k,1) = at3*zt1
      fxyz(3,j,k,1) = zero
      fxyz(1,j,k1,1) = at2*zt2
      fxyz(2,j,k1,1) = -at3*zt2
      fxyz(3,j,k1,1) = zero
      fxyz(1,j,k,l1) = zero
      fxyz(2,j,k,l1) = zero
      fxyz(3,j,k,l1) = zero
      fxyz(1,j,k1,l1) = zero
      fxyz(2,j,k1,l1) = zero
      fxyz(3,j,k1,l1) = zero
      at1 = at1*(q(j,k,1)*conjg(q(j,k,1)) + q(j,k1,1)*conjg(q(j,k1,1)))
      wp = wp + dble(at1)
  100 continue
c mode numbers kx = 0, nx/2
      zt1 = ffc(1,k,1)
      at1 = real(zt1)*aimag(zt1)
      at3 = dny*real(k - 1)*at1
      fxyz(1,1,k,1) = zero
      zt1 = q(1,k,1)
      fxyz(2,1,k,1) = at3*cmplx(aimag(zt1),-real(zt1))
      fxyz(3,1,k,1) = zero
      fxyz(1,1,k1,1) = zero
      fxyz(2,1,k1,1) = zero
      fxyz(3,1,k1,1) = zero
      fxyz(1,1,k,l1) = zero
      fxyz(2,1,k,l1) = zero
      fxyz(3,1,k,l1) = zero
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      at1 = at1*(q(1,k,1)*conjg(q(1,k,1)))
      wp = wp + dble(at1)
      sum2 = sum2 + wp
  110 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
!dir$ ivdep
      do 120 j = 2, nxh
      zt1 = ffc(j,1,1)
      at1 = real(zt1)*aimag(zt1)
      at2 = dnx*real(j - 1)*at1
      zt1 = q(j,1,1)
      fxyz(1,j,1,1) = at2*cmplx(aimag(zt1),-real(zt1))
      fxyz(2,j,1,1) = zero
      fxyz(3,j,1,1) = zero
      fxyz(1,j,k1,1) = zero
      fxyz(2,j,k1,1) = zero
      fxyz(3,j,k1,1) = zero
      fxyz(1,j,1,l1) = zero
      fxyz(2,j,1,l1) = zero
      fxyz(3,j,1,l1) = zero
      fxyz(1,j,k1,l1) = zero
      fxyz(2,j,k1,l1) = zero
      fxyz(3,j,k1,l1) = zero
      at1 = at1*(q(j,1,1)*conjg(q(j,1,1)))
      wp = wp + dble(at1)
  120 continue
      fxyz(1,1,1,1) = zero
      fxyz(2,1,1,1) = zero
      fxyz(3,1,1,1) = zero
      fxyz(1,1,k1,1) = zero
      fxyz(2,1,k1,1) = zero
      fxyz(3,1,k1,1) = zero
      fxyz(1,1,1,l1) = zero
      fxyz(2,1,1,l1) = zero
      fxyz(3,1,1,l1) = zero
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      we = real(nx)*real(ny)*real(nz)*(sum1 + sum2 + wp)
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
c this subroutine calculates tables needed by a three dimensional
c real to complex fast fourier transform and its inverse.
c input: indx, indy, indz, nxhyzd, nxyzhd
c output: mixup, sct
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, indz, nxhyzd, nxyzhd
      integer mixup
      complex sct
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, ny, nz, nxyz, nxhyz, nxyzh
      integer j, k, lb, ll, jb, it
      real dnxyz, arg
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhyz
      lb = j - 1
      ll = 0
      do 10 k = 1, ndx1yz
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles 2*n*pi/nxyz
      nxyzh = nxyz/2
      dnxyz = 6.28318530717959/real(nxyz)
      do 30 j = 1, nxyzh
      arg = dnxyz*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT3RVMX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd
     1,nxhyzd,nxyzhd)
c wrapper function for real to complex fft, with packed data
c parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, indz, nxhd, nyd, nzd, nxhyzd, nxyzhd
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer ny, nz, nyi, nzi
      data nyi, nzi /1,1/
c calculate range of indices
      ny = 2**indy
      nz = 2**indz
c inverse fourier transform
      if (isign.lt.0) then
c perform xy fft
         call FFT3RVMXY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd
     1,nzd,nxhyzd,nxyzhd)
c perform z fft
         call FFT3RVMXZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd
     1,nzd,nxhyzd,nxyzhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform z fft
         call FFT3RVMXZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd
     1,nzd,nxhyzd,nxyzhd)
c perform xy fft
         call FFT3RVMXY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd
     1,nzd,nxhyzd,nxyzhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT3RVM3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd
     1,nxhyzd,nxyzhd)
c wrapper function for 3 2d real to complex ffts, with packed data
c parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, indz, nxhd, nyd, nzd, nxhyzd, nxyzhd
      dimension f(4,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer ny, nz, nyi, nzi
      data nyi, nzi /1,1/
c calculate range of indices
      ny = 2**indy
      nz = 2**indz
c inverse fourier transform
      if (isign.lt.0) then
c perform xy fft
         call FFT3RVM3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,  
     1nyd,nzd,nxhyzd,nxyzhd)
c perform z fft
         call FFT3RVM3Z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd
     1,nzd,nxhyzd,nxyzhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform z fft
         call FFT3RVM3Z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd
     1,nzd,nxhyzd,nxyzhd)
c perform xy fft
         call FFT3RVM3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,  
     1nyd,nzd,nxhyzd,nxyzhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RVMXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,nxhd
     1,nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the x-y part of a three dimensional real to
c complex fast fourier transform and its inverse, for a subset of z,
c using complex arithmetic, with Vector/OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform in x and y is performed
c f(n,m,i) = (1/nx*ny*nz)*sum(f(j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)*
c       exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform in x and y is performed
c f(j,k,l) = sum(f(n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nzi = initial z index used
c nzp = number of z indices used
c nxhd = first dimension of f
c nyd,nzd = second and third dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1,1,1)) = real part of mode nx/2,0,0
c imag(f(1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nxyz, nxhyz, nzt, nrx, nry, nrxb, nryb
      integer i, j, k, l, n, j1, k1, k2, ns, ns2, km, kmr
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nzt = nzi + nzp - 1
      if (isign.gt.0) go to 180
c inverse fourier transform
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
      nryb = nxhyz/ny
      nry = nxyz/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,ani,t1,t2,t3)
      do 170 n = nzi, nzt
c bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 10 i = 1, ny
         t1 = f(j1,i,n)
         f(j1,i,n) = f(j,i,n)
         f(j,i,n) = t1
   10    continue
      endif
   20 continue
c first transform in x
      do 60 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      t1 = sct(1+kmr*(j-1))
      do 30 i = 1, ny
      t2 = t1*f(j+k2,i,n)
      f(j+k2,i,n) = f(j+k1,i,n) - t2
      f(j+k1,i,n) = f(j+k1,i,n) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = 1, ny
      t2 = conjg(f(nxh2-j,k,n))
      t1 = f(j,k,n) + t2
      t2 = (f(j,k,n) - t2)*t3
      f(j,k,n) = ani*(t1 + t2)
      f(nxh2-j,k,n) = ani*conjg(t1 - t2)
   70 continue
   80 continue
      ani = 2.0*ani
      do 90 k = 1, ny
      f(nxhh+1,k,n) = ani*conjg(f(nxhh+1,k,n))
      f(1,k,n) = ani*cmplx(real(f(1,k,n)) + aimag(f(1,k,n)),
     1                     real(f(1,k,n)) - aimag(f(1,k,n)))
   90 continue
c bit-reverse array elements in y
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 100 i = 1, nxh
         t1 = f(i,k1,n)
         f(i,k1,n) = f(i,k,n)
         f(i,k,n) = t1
  100    continue
      endif
  110 continue
c then transform in y
      do 150 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      t1 = sct(1+kmr*(j-1))
      do 120 i = 1, nxh
      t2 = t1*f(i,j+k2,n)
      f(i,j+k2,n) = f(i,j+k1,n) - t2
      f(i,j+k1,n) = f(i,j+k1,n) + t2
  120 continue
  130 continue
  140 continue
  150 continue
c unscramble modes kx = 0, nx/2
!dir$ ivdep
      do 160 k = 2, nyh
      t1 = f(1,ny2-k,n)
      f(1,ny2-k,n) = 0.5*cmplx(aimag(f(1,k,n) + t1),real(f(1,k,n) - t1))
      f(1,k,n) = 0.5*cmplx(real(f(1,k,n) + t1),aimag(f(1,k,n) - t1))
  160 continue
  170 continue
!$OMP END PARALLEL DO
      return
c forward fourier transform
  180 nryb = nxhyz/ny
      nry = nxyz/ny
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,t1,t2,t3)
      do 350 n = nzi, nzt
c scramble modes kx = 0, nx/2
!dir$ ivdep
      do 190 k = 2, nyh
      t1 = cmplx(aimag(f(1,ny2-k,n)),real(f(1,ny2-k,n)))
      f(1,ny2-k,n) = conjg(f(1,k,n) - t1)
      f(1,k,n) = f(1,k,n) + t1
  190 continue
c bit-reverse array elements in y
      do 210 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 200 i = 1, nxh
         t1 = f(i,k1,n)
         f(i,k1,n) = f(i,k,n)
         f(i,k,n) = t1
  200    continue
      endif
  210 continue
c then transform in y
      do 250 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      t1 = conjg(sct(1+kmr*(j-1)))
      do 220 i = 1, nxh
      t2 = t1*f(i,j+k2,n)
      f(i,j+k2,n) = f(i,j+k1,n) - t2
      f(i,j+k1,n) = f(i,j+k1,n) + t2
  220 continue
  230 continue
  240 continue
  250 continue
c scramble coefficients
      kmr = nxyz/nx
      do 270 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 260 k = 1, ny
      t2 = conjg(f(nxh2-j,k,n))
      t1 = f(j,k,n) + t2
      t2 = (f(j,k,n) - t2)*t3
      f(j,k,n) = t1 + t2
      f(nxh2-j,k,n) = conjg(t1 - t2)
  260 continue
  270 continue
      do 280 k = 1, ny
      f(nxhh+1,k,n) = 2.0*conjg(f(nxhh+1,k,n))
      f(1,k,n) = cmplx(real(f(1,k,n)) + aimag(f(1,k,n)),
     1                 real(f(1,k,n)) - aimag(f(1,k,n)))
  280 continue
c bit-reverse array elements in x
      do 300 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 290 i = 1, ny
         t1 = f(j1,i,n)
         f(j1,i,n) = f(j,i,n)
         f(j,i,n) = t1
  290    continue
      endif
  300 continue
c finally transform in x
      do 340 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 330 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 320 j = 1, ns
      t1 = conjg(sct(1+kmr*(j-1)))
      do 310 i = 1, ny
      t2 = t1*f(j+k2,i,n)
      f(j+k2,i,n) = f(j+k1,i,n) - t2
      f(j+k1,i,n) = f(j+k1,i,n) + t2
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RVMXZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,nxhd
     1,nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the z part of a three dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic, with Vector/OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform in z is performed
c f(j,k,l) = sum(f(j,k,i)*exp(-sqrt(-1)*2pi*l*i/nz))
c if isign = 1, a forward fourier transform in z is performed
c f(n,m,i) = sum(f(n,m,l)*exp(sqrt(-1)*2pi*l*i/nz))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f
c nyd,nzd = second and third dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1,1,1)) = real part of mode nx/2,0,0
c imag(f(1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh
      integer nz, nzh, nz2, nxyz, nxhyz, nyt, nrz, nrzb
      integer i, j, k, l, n, k1, k2, l1, ns, ns2, km, kmr
      complex t1, t2
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 100
c inverse fourier transform
      nrzb = nxhyz/nz
      nrz = nxyz/nz
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,l1,t1,t2)
      do 70 n = nyi, nyt
c bit-reverse array elements in z
      do 20 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         do 10 i = 1, nxh
         t1 = f(i,n,l1)
         f(i,n,l1) = f(i,n,l)
         f(i,n,l) = t1
   10    continue
      endif
   20 continue
c finally transform in z
      do 60 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      t1 = sct(1+kmr*(j-1))
      do 30 i = 1, nxh
      t2 = t1*f(i,n,j+k2)
      f(i,n,j+k2) = f(i,n,j+k1) - t2
      f(i,n,j+k1) = f(i,n,j+k1) + t2
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if (nyi.eq.1) then
!dir$ ivdep
         do 80 n = 2, nzh
         t1 = f(1,1,nz2-n)
         f(1,1,nz2-n) = 0.5*cmplx(aimag(f(1,1,n) + t1),
     1                            real(f(1,1,n) - t1))
         f(1,1,n) = 0.5*cmplx(real(f(1,1,n) + t1),aimag(f(1,1,n) - t1))
   80    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
!dir$ ivdep
         do 90 n = 2, nzh
         t1 = f(1,nyh+1,nz2-n)
         f(1,nyh+1,nz2-n) = 0.5*cmplx(aimag(f(1,nyh+1,n) + t1),
     1                                real(f(1,nyh+1,n) - t1))
         f(1,nyh+1,n) = 0.5*cmplx(real(f(1,nyh+1,n) + t1),
     1                            aimag(f(1,nyh+1,n) - t1))
   90    continue
      endif
      return
c forward fourier transform
  100 nrzb = nxhyz/nz
      nrz = nxyz/nz
c scramble modes kx = 0, nx/2
      if (nyi.eq.1) then
!dir$ ivdep
         do 110 n = 2, nzh
         t1 = cmplx(aimag(f(1,1,nz2-n)),real(f(1,1,nz2-n)))
         f(1,1,nz2-n) = conjg(f(1,1,n) - t1)
         f(1,1,n) = f(1,1,n) + t1
  110    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
!dir$ ivdep
         do 120 n = 2, nzh
         t1 = cmplx(aimag(f(1,nyh+1,nz2-n)),real(f(1,nyh+1,nz2-n)))
         f(1,nyh+1,nz2-n) = conjg(f(1,nyh+1,n) - t1)
         f(1,nyh+1,n) = f(1,nyh+1,n) + t1
  120    continue
      endif
c bit-reverse array elements in z
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,l1,t1,t2)
      do 190 n = nyi, nyt
      do 140 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         do 130 i = 1, nxh
         t1 = f(i,n,l1)
         f(i,n,l1) = f(i,n,l)
         f(i,n,l) = t1
  130    continue
      endif
  140 continue
c first transform in z
      do 180 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      t1 = conjg(sct(1+kmr*(j-1)))
      do 150 i = 1, nxh
      t2 = t1*f(i,n,j+k2)
      f(i,n,j+k2) = f(i,n,j+k1) - t2
      f(i,n,j+k1) = f(i,n,j+k1) + t2
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RVM3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,   
     1nxhd,nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the x-y part of 3 three dimensional complex
c to real fast fourier transforms and their inverses, for a subset of z,
c using complex arithmetic, with Vector/OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, three inverse fourier transforms in x and y are
c performed
c f(1:3,n,m,i) = (1/nx*ny*nz)*sum(f(1:3,j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
c       *exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, three forward fourier transforms in x and y are
c performed
c f(1:3,j,k,l) = sum(f(1:3,n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nzi = initial z index used
c nzp = number of z indices used
c nxhd = second dimension of f
c nyd,nzd = third and fourth dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(1:3,j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1:3,1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1:3,1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1:3,1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1:3,1,1,1)) = real part of mode nx/2,0,0
c imag(f(1:3,1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1:3,1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1:3,1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
      integer nxhyzd,nxyzhd
      complex f, sct
      integer mixup
      dimension f(4,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nxyz, nxhyz, nzt, nrx, nry, nrxb, nryb
      integer i, j, k, l, n, jj, j1, k1, k2, ns, ns2, km, kmr
      real at1, at2, ani
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nzt = nzi + nzp - 1
      if (isign.gt.0) go to 230
c inverse fourier transform
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
      nryb = nxhyz/ny
      nry = nxyz/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,jj,j1,at1,at2,ani,t1,t2,t3,
!$OMP& t4)
      do 220 n = nzi, nzt
c swap complex components
      do 20 i = 1, ny
      do 10 j = 1, nxh
      at1 = aimag(f(3,j,i,n))
      at2 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(real(f(2,j,i,n)),real(f(4,j,i,n)))
      f(2,j,i,n) = cmplx(aimag(f(1,j,i,n)),at1)
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
   10 continue
   20 continue
c bit-reverse array elements in x
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
      do 30 i = 1, ny
         t1 = f(1,j1,i,n)
         t2 = f(2,j1,i,n)
         t3 = f(3,j1,i,n)
         f(1,j1,i,n) = f(1,j,i,n)
         f(2,j1,i,n) = f(2,j,i,n)
         f(3,j1,i,n) = f(3,j,i,n)
         f(1,j,i,n) = t1
         f(2,j,i,n) = t2
         f(3,j,i,n) = t3
   30    continue
      endif
   40 continue
c first transform in x
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      t1 = sct(1+kmr*(j-1))
      do 50 i = 1, ny
      t2 = t1*f(1,j+k2,i,n)
      t3 = t1*f(2,j+k2,i,n)
      t4 = t1*f(3,j+k2,i,n)
      f(1,j+k2,i,n) = f(1,j+k1,i,n) - t2
      f(2,j+k2,i,n) = f(2,j+k1,i,n) - t3
      f(3,j+k2,i,n) = f(3,j+k1,i,n) - t4
      f(1,j+k1,i,n) = f(1,j+k1,i,n) + t2
      f(2,j+k1,i,n) = f(2,j+k1,i,n) + t3
      f(3,j+k1,i,n) = f(3,j+k1,i,n) + t4
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = 1, ny
      do 90 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k,n))
      t1 = f(jj,j,k,n) + t2
      t2 = (f(jj,j,k,n) - t2)*t3
      f(jj,j,k,n) = ani*(t1 + t2)
      f(jj,nxh2-j,k,n) = ani*conjg(t1 - t2)
   90 continue
  100 continue
  110 continue
      ani = 2.0*ani
      do 130 k = 1, ny
      do 120 jj = 1, 3
      f(jj,nxhh+1,k,n) = ani*conjg(f(jj,nxhh+1,k,n))
      f(jj,1,k,n) = ani*cmplx(real(f(jj,1,k,n)) + aimag(f(jj,1,k,n)),
     1                        real(f(jj,1,k,n)) - aimag(f(jj,1,k,n)))
  120 continue
  130 continue
c bit-reverse array elements in y
      do 150 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 140 i = 1, nxh
         t1 = f(1,i,k1,n)
         t2 = f(2,i,k1,n)
         t3 = f(3,i,k1,n)
         f(1,i,k1,n) = f(1,i,k,n)
         f(2,i,k1,n) = f(2,i,k,n)
         f(3,i,k1,n) = f(3,i,k,n)
         f(1,i,k,n) = t1
         f(2,i,k,n) = t2
         f(3,i,k,n) = t3
  140    continue
      endif
  150 continue
c then transform in y
      do 190 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      t1 = sct(1+kmr*(j-1))
      do 160 i = 1, nxh
      t2 = t1*f(1,i,j+k2,n)
      t3 = t1*f(2,i,j+k2,n)
      t4 = t1*f(3,i,j+k2,n)
      f(1,i,j+k2,n) = f(1,i,j+k1,n) - t2
      f(2,i,j+k2,n) = f(2,i,j+k1,n) - t3
      f(3,i,j+k2,n) = f(3,i,j+k1,n) - t4
      f(1,i,j+k1,n) = f(1,i,j+k1,n) + t2
      f(2,i,j+k1,n) = f(2,i,j+k1,n) + t3
      f(3,i,j+k1,n) = f(3,i,j+k1,n) + t4
  160 continue
  170 continue
  180 continue
  190 continue
c unscramble modes kx = 0, nx/2
      do 210 k = 2, nyh
      do 200 jj = 1, 3
      t1 = f(jj,1,ny2-k,n)
      f(jj,1,ny2-k,n) = 0.5*cmplx(aimag(f(jj,1,k,n) + t1),
     1                            real(f(jj,1,k,n) - t1))
      f(jj,1,k,n) = 0.5*cmplx(real(f(jj,1,k,n) + t1),
     1                        aimag(f(jj,1,k,n) - t1))
  200 continue
  210 continue
  220 continue
!$OMP END PARALLEL DO
      return
c forward fourier transform
  230 nryb = nxhyz/ny
      nry = nxyz/ny
      nrxb = nxhyz/nxh
      nrx = nxyz/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,jj,j1,at1,at2,t1,t2,t3,t4)
      do 450 n = nzi, nzt
c scramble modes kx = 0, nx/2
      do 250 k = 2, nyh
      do 240 jj = 1, 3
      t1 = cmplx(aimag(f(jj,1,ny2-k,n)),real(f(jj,1,ny2-k,n)))
      f(jj,1,ny2-k,n) = conjg(f(jj,1,k,n) - t1)
      f(jj,1,k,n) = f(jj,1,k,n) + t1
  240 continue
  250 continue
c bit-reverse array elements in y
      do 270 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 260 i = 1, nxh
         t1 = f(1,i,k1,n)
         t2 = f(2,i,k1,n)
         t3 = f(3,i,k1,n)
         f(1,i,k1,n) = f(1,i,k,n)
         f(2,i,k1,n) = f(2,i,k,n)
         f(3,i,k1,n) = f(3,i,k,n)
         f(1,i,k,n) = t1
         f(2,i,k,n) = t2
         f(3,i,k,n) = t3
  260    continue
      endif
  270 continue
c then transform in y
      do 310 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 300 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 290 j = 1, ns
      t1 = conjg(sct(1+kmr*(j-1)))
      do 280 i = 1, nxh
      t2 = t1*f(1,i,j+k2,n)
      t3 = t1*f(2,i,j+k2,n)
      t4 = t1*f(3,i,j+k2,n)
      f(1,i,j+k2,n) = f(1,i,j+k1,n) - t2
      f(2,i,j+k2,n) = f(2,i,j+k1,n) - t3
      f(3,i,j+k2,n) = f(3,i,j+k1,n) - t4
      f(1,i,j+k1,n) = f(1,i,j+k1,n) + t2
      f(2,i,j+k1,n) = f(2,i,j+k1,n) + t3
      f(3,i,j+k1,n) = f(3,i,j+k1,n) + t4
  280 continue
  290 continue
  300 continue
  310 continue
c scramble coefficients
      kmr = nxyz/nx
      do 340 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 330 k = 1, ny
      do 320 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k,n))
      t1 = f(jj,j,k,n) + t2
      t2 = (f(jj,j,k,n) - t2)*t3
      f(jj,j,k,n) = t1 + t2
      f(jj,nxh2-j,k,n) = conjg(t1 - t2)
  320 continue
  330 continue
  340 continue
      do 360 k = 1, ny
      do 350 jj = 1, 3
      f(jj,nxhh+1,k,n) = 2.0*conjg(f(jj,nxhh+1,k,n))
      f(jj,1,k,n) = cmplx(real(f(jj,1,k,n)) + aimag(f(jj,1,k,n)),
     1                    real(f(jj,1,k,n)) - aimag(f(jj,1,k,n)))
  350 continue
  360 continue
c bit-reverse array elements in x
      do 380 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 370 i = 1, ny
         t1 = f(1,j1,i,n)
         t2 = f(2,j1,i,n)
         t3 = f(3,j1,i,n)
         f(1,j1,i,n) = f(1,j,i,n)
         f(2,j1,i,n) = f(2,j,i,n)
         f(3,j1,i,n) = f(3,j,i,n)
         f(1,j,i,n) = t1
         f(2,j,i,n) = t2
         f(3,j,i,n) = t3
  370    continue
      endif
  380 continue
c finally transform in x
      do 420 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 410 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 400 j = 1, ns
      t1 = conjg(sct(1+kmr*(j-1)))
      do 390 i = 1, ny
      t2 = t1*f(1,j+k2,i,n)
      t3 = t1*f(2,j+k2,i,n)
      t4 = t1*f(3,j+k2,i,n)
      f(1,j+k2,i,n) = f(1,j+k1,i,n) - t2
      f(2,j+k2,i,n) = f(2,j+k1,i,n) - t3
      f(3,j+k2,i,n) = f(3,j+k1,i,n) - t4
      f(1,j+k1,i,n) = f(1,j+k1,i,n) + t2
      f(2,j+k1,i,n) = f(2,j+k1,i,n) + t3
      f(3,j+k1,i,n) = f(3,j+k1,i,n) + t4
  390 continue
  400 continue
  410 continue
  420 continue
c swap complex components
      do 440 i = 1, ny
      do 430 j = 1, nxh
      f(4,j,i,n) = cmplx(aimag(f(3,j,i,n)),aimag(f(4,j,i,n)))
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(aimag(f(1,j,i,n)),aimag(f(2,j,i,n)))
      at2 = real(f(2,j,i,n))
      f(2,j,i,n) = cmplx(at1,0.0)
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
  430 continue
  440 continue
  450 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RVM3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,nxhd
     1,nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the z part of 3 three dimensional complex to
c real fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic, with Vector/OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, three inverse fourier transforms in z are performed
c f(1:3,j,k,l) = sum(f(1:3,j,k,i)*exp(-sqrt(-1)*2pi*l*i/nz))
c if isign = 1, three forward fourier transforms in z are performed
c f(1:3,n,m,i) = sum(f(1:3,n,m,l)*exp(sqrt(-1)*2pi*l*i/nz))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = second dimension of f
c nyd,nzd = third and fourth dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(1:3,2*j-1,k,l),f(2*j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1:3,1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1:3,1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1:3,1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1:3,1,1,1)) = real part of mode nx/2,0,0
c imag(f(1:3,1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1:3,1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1:3,1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(4,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh
      integer nz, nzh, nz2, nxyz, nxhyz, nyt, nrz, nrzb
      integer i, j, k, l, n, jj, k1, k2, l1, ns, ns2, km, kmr
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 120
c inverse fourier transform
      nrzb = nxhyz/nz
      nrz = nxyz/nz
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,l1,t1,t2,t3,t4)
      do 70 n = nyi, nyt
c bit-reverse array elements in z
      do 20 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
      do 10 i = 1, nxh
         t1 = f(1,i,n,l1)
         t2 = f(2,i,n,l1)
         t3 = f(3,i,n,l1)
         f(1,i,n,l1) = f(1,i,n,l)
         f(2,i,n,l1) = f(2,i,n,l)
         f(3,i,n,l1) = f(3,i,n,l)
         f(1,i,n,l) = t1
         f(2,i,n,l) = t2
         f(3,i,n,l) = t3
   10    continue
      endif
   20 continue
c finally transform in z
      do 60 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      t1 = sct(1+kmr*(j-1))
      do 30 i = 1, nxh
      t2 = t1*f(1,i,n,j+k2)
      t3 = t1*f(2,i,n,j+k2)
      t4 = t1*f(3,i,n,j+k2)
      f(1,i,n,j+k2) = f(1,i,n,j+k1) - t2
      f(2,i,n,j+k2) = f(2,i,n,j+k1) - t3
      f(3,i,n,j+k2) = f(3,i,n,j+k1) - t4
      f(1,i,n,j+k1) = f(1,i,n,j+k1) + t2
      f(2,i,n,j+k1) = f(2,i,n,j+k1) + t3
      f(3,i,n,j+k1) = f(3,i,n,j+k1) + t4
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if (nyi.eq.1) then
         do 90 n = 2, nzh
         do 80 jj = 1, 3
         t1 = f(jj,1,1,nz2-n)
         f(jj,1,1,nz2-n) = 0.5*cmplx(aimag(f(jj,1,1,n) + t1),
     1                               real(f(jj,1,1,n) - t1))
         f(jj,1,1,n) = 0.5*cmplx(real(f(jj,1,1,n) + t1),
     1                           aimag(f(jj,1,1,n) - t1))
   80    continue
   90    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         do 110 n = 2, nzh
         do 100 jj = 1, 3
         t1 = f(jj,1,nyh+1,nz2-n)
         f(jj,1,nyh+1,nz2-n) = 0.5*cmplx(aimag(f(jj,1,nyh+1,n) + t1),
     1                                  real(f(jj,1,nyh+1,n) - t1))
         f(jj,1,nyh+1,n) = 0.5*cmplx(real(f(jj,1,nyh+1,n) + t1),
     1                              aimag(f(jj,1,nyh+1,n) - t1))
  100    continue
  110    continue
      endif
      return
c forward fourier transform
  120 nrzb = nxhyz/nz
      nrz = nxyz/nz
c scramble modes kx = 0, nx/2
      if (nyi.eq.1) then
         do 140 n = 2, nzh
         do 130 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,1,nz2-n)),real(f(jj,1,1,nz2-n)))
         f(jj,1,1,nz2-n) = conjg(f(jj,1,1,n) - t1)
         f(jj,1,1,n) = f(jj,1,1,n) + t1
  130    continue
  140    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         do 160 n = 2, nzh
         do 150 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,nyh+1,nz2-n)),
     1              real(f(jj,1,nyh+1,nz2-n)))
         f(jj,1,nyh+1,nz2-n) = conjg(f(jj,1,nyh+1,n) - t1)
         f(jj,1,nyh+1,n) = f(jj,1,nyh+1,n) + t1
  150    continue
  160    continue
      endif
c bit-reverse array elements in z
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,l1,t1,t2,t3,t4)
      do 230 n = nyi, nyt
      do 180 l = 1, nz
      l1 = (mixup(l) - 1)/nrzb + 1
      if (l.lt.l1) then
         do 170 i = 1, nxh
         t1 = f(1,i,n,l1)
         t2 = f(2,i,n,l1)
         t3 = f(3,i,n,l1)
         f(1,i,n,l1) = f(1,i,n,l)
         f(2,i,n,l1) = f(2,i,n,l)
         f(3,i,n,l1) = f(3,i,n,l)
         f(1,i,n,l) = t1
         f(2,i,n,l) = t2
         f(3,i,n,l) = t3
  170    continue
      endif
  180 continue
c first transform in z
      do 220 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 210 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 200 j = 1, ns
      t1 = conjg(sct(1+kmr*(j-1)))
      do 190 i = 1, nxh
      t2 = t1*f(1,i,n,j+k2)
      t3 = t1*f(2,i,n,j+k2)
      t4 = t1*f(3,i,n,j+k2)
      f(1,i,n,j+k2) = f(1,i,n,j+k1) - t2
      f(2,i,n,j+k2) = f(2,i,n,j+k1) - t3
      f(3,i,n,j+k2) = f(3,i,n,j+k1) - t4
      f(1,i,n,j+k1) = f(1,i,n,j+k1) + t2
      f(2,i,n,j+k1) = f(2,i,n,j+k1) + t3
      f(3,i,n,j+k1) = f(3,i,n,j+k1) + t4
  190 continue
  200 continue
  210 continue
  220 continue
  230 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine SET_SZERO3(q,mx,my,mz,nxv,nyv,nzv,mx1,my1,mxyz1)
c for 3d code, this subroutine zeros out charge density array.
c for Intel NUMA architecture with first touch policy, this associates
c array segments with appropriate threads
c OpenMP version
c input: all, output: q
c q(j,k,l) = charge density at grid point j,k,l
c mx/my/mz = number of grids in sorting cell in x/y/z
c nxv = first dimension of charge array, must be >= nx+ng
c nyv = second dimension of charge array, must be >= ny+ng
c nzv = third dimension of charge array, must be >= nz+ng
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
      implicit none
      integer mx, my, mz, nxv, nyv, nzv, mx1, my1, mxyz1
      real q
      dimension q(nxv,nyv,nzv)
c local data
      integer mxy1, mz1, noff, moff, loff
      integer i, j, k, l, nn, mm, ll
      mxy1 = mx1*my1
      mz1 = mxyz1/mxy1
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,noff,moff,loff,nn,mm,ll)
      do 40 l = 1, mxyz1
      i = (l - 1)/mxy1
      k = l - mxy1*i
      loff = mz*i
      ll = mz
      if ((i+1).eq.mz1) ll = nzv - loff
      j = (k - 1)/mx1
      moff = my*j
      mm = my
      if ((j+1).eq.my1) mm = nyv - moff
      k = k - mx1*j
      noff = mx*(k - 1)
      nn = mx
      if (k.eq.mx1) nn = nxv - noff
c zero charge in global array
      do 30 k = 1, ll
      do 20 j = 1, mm
!dir$ ivdep
      do 10 i = 1, nn
      q(i+noff,j+moff,k+loff) = 0.0
   10 continue
   20 continue
   30 continue
   40 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine SET_VZERO3(cu,mx,my,mz,ndim,nxv,nyv,nzv,mx1,my1,mxyz1)
c for 3d code, this subroutine zeros out current density array.
c for Intel NUMA architecture with first touch policy, this associates
c array segments with appropriate threads
c OpenMP version
c input: all, output: cu
c cu(m,j,k,l) = charge density at grid point m,j,k,l
c mx/my/mz = number of grids in sorting cell in x/y/z
c ndim = first dimension of current array
c nxv = second dimension of current array, must be >= nx+ng
c nyv = third dimension of current array, must be >= ny+ng
c nzv = fourth dimension of current array, must be >= nz+ng
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c mxyz1 = mx1*my1*mz1,
c where mz1 = (system length in z direction - 1)/mz + 1
      implicit none
      integer mx, my, mz, ndim, nxv, nyv, nzv, mx1, my1, mxyz1
      real cu
      dimension cu(ndim,nxv,nyv,nzv)
c local data
      integer mxy1, mz1, noff, moff, loff
      integer i, j, k, l, m, nn, mm, ll
      mxy1 = mx1*my1
      mz1 = mxyz1/mxy1
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,m,noff,moff,loff,nn,mm,ll)
      do 50 l = 1, mxyz1
      i = (l - 1)/mxy1
      k = l - mxy1*i
      loff = mz*i
      ll = mz
      if ((i+1).eq.mz1) ll = nzv - loff
      j = (k - 1)/mx1
      moff = my*j
      mm = my
      if ((j+1).eq.my1) mm = nyv - moff
      k = k - mx1*j
      noff = mx*(k - 1)
      nn = mx
      if (k.eq.mx1) nn = nxv - noff
c zero current in global array
      do 40 k = 1, ll
      do 30 j = 1, mm
      do 20 i = 1, nn
      do 10 m = 1, ndim
      cu(m,i+noff,j+moff,k+loff) = 0.0
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine SET_CVZERO3(exyz,nx,ny,nz,ndim,nxvh,nyv,nzv)
c for 3d code, this subroutine zeros out transverse field array.
c for Intel NUMA architecture with first touch policy, this associates
c array segments with appropriate threads
c OpenMP version
c input: all, output: exyz
c exyz(i,j,k,l) = complex transverse electric field
c nx/ny/nz = system length in x/y/z direction
c ndim = first dimension of field array
c nxvh = second dimension of field array, must be >= nxh
c nyv = third dimension of field array, must be >= ny
c nzv = fourth dimension of field array, must be >= nz
      implicit none
      integer nx, ny, nz, ndim, nxvh, nyv, nzv
      complex exyz
      dimension exyz(ndim,nxvh,nyv,nzv)
c local data
      integer nxh, nyh, nzh, ny2, nz2, i, j, k, l, k1, l1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
      zero = cmplx(0.0,0.0)
c loop over mode numbers
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
!$OMP PARALLEL
!$OMP DO PRIVATE(i,j,k,l,k1,l1)
      do 90 l = 2, nzh
      l1 = nz2 - l
      do 30 k = 2, nyh
      k1 = ny2 - k
      do 20 j = 2, nxh
      do 10 i = 1, ndim
      exyz(i,j,k,l) = zero
      exyz(i,j,k1,l) = zero
      exyz(i,j,k,l1) = zero
      exyz(i,j,k1,l1) = zero
   10 continue
   20 continue
   30 continue
c mode numbers kx = 0, nx/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      do 40 i = 1, ndim
      exyz(i,1,k,l) = zero
      exyz(i,1,k1,l) = zero
      exyz(i,1,k,l1) = zero
      exyz(i,1,k1,l1) = zero
   40 continue
   50 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      do 60 i = 1, ndim
      exyz(i,j,1,l) = zero
      exyz(i,j,k1,l) = zero
      exyz(i,j,1,l1) = zero
      exyz(i,j,k1,l1) = zero
   60 continue
   70 continue
c mode numbers kx = 0, nx/2
      do 80 i = 1, ndim
      exyz(i,1,1,l) = zero
      exyz(i,1,k1,l) = zero
      exyz(i,1,1,l1) = zero
      exyz(i,1,k1,l1) = zero
   80 continue
   90 continue
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      l1 = nzh + 1
c mode numbers kz = 0, nz/2
!$OMP PARALLEL DO PRIVATE(i,j,k,k1)
      do 130 k = 2, nyh
      k1 = ny2 - k
      do 110 j = 2, nxh
      do 100 i = 1, ndim
      exyz(i,j,k,1) = zero
      exyz(i,j,k1,1) = zero
      exyz(i,j,k,l1) = zero
      exyz(i,j,k1,l1) = zero
  100 continue
  110 continue
c mode numbers kx = 0, nx/2
      do 120 i = 1, ndim
      exyz(i,1,k,1) = zero
      exyz(i,1,k1,1) = zero
      exyz(i,1,k,l1) = zero
      exyz(i,1,k1,l1) = zero
  120 continue
  130 continue
!$OMP END PARALLEL DO
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 150 j = 2, nxh
      do 140 i = 1, ndim
      exyz(i,j,1,1) = zero
      exyz(i,j,k1,1) = zero
      exyz(i,j,1,l1) = zero
      exyz(i,j,k1,l1) = zero
  140 continue
  150 continue
      do 160 i = 1, ndim
      exyz(i,1,1,1) = zero
      exyz(i,1,k1,1) = zero
      exyz(i,1,1,l1) = zero
      exyz(i,1,k1,l1) = zero
  160 continue
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
