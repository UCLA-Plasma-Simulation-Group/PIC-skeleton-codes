c Fortran Library for Skeleton 2D Electrostatic OpenMP PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR2(part,vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx,ny,   
     1ipbc)
c for 2d code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c vtx/vty = thermal velocity of electrons in x/y direction
c vdx/vdy = drift velocity of beam electrons in x/y direction
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, idimp, nop, nx, ny, ipbc
      real vtx, vty, vdx, vdy
      real part
      dimension part(idimp,nop)
c local data
      integer j, k, k1, npxy
      real edgelx, edgely, at1, at2, at3, sum1, sum2
      double precision dsum1, dsum2
      double precision ranorm
      npxy = npx*npy
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
      k1 = npx*(k - 1)
      at3 = edgely + at2*(real(k) - 0.5)
      do 10 j = 1, npx
      part(1,j+k1) = edgelx + at1*(real(j) - 0.5)
      part(2,j+k1) = at3
   10 continue
   20 continue
c maxwellian velocity distribution
      do 30 j = 1, npxy
      part(3,j) = vtx*ranorm()
      part(4,j) = vty*ranorm()
   30 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      do 40 j = 1, npxy
      dsum1 = dsum1 + part(3,j)
      dsum2 = dsum2 + part(4,j)
   40 continue
      sum1 = dsum1
      sum2 = dsum2
      at1 = 1.0/real(npxy)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      do 50 j = 1, npxy
      part(3,j) = part(3,j) - sum1
      part(4,j) = part(4,j) - sum2
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DBLKP2L(part,kpic,nppmx,idimp,nop,mx,my,mx1,mxy1,irc)
c this subroutine finds the maximum number of particles in each tile of
c mx, my to calculate size of segmented particle array ppart
c linear interpolation
c input: all except kpic, nppmx, output: kpic, nppmx
c part = input particle array
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c kpic = output number of particles per tile
c nppmx = return maximum number of particles in tile
c idimp = size of phase space = 4
c nop = number of particles
c mx/my = number of grids in sorting cell in x and y
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kpic, nppmx, idimp, nop, mx, my, mx1, mxy1, irc
      real part
      dimension part(idimp,nop), kpic(mxy1)
c local data
      integer j, k, n, m, isum, ist, npx, ierr
      ierr = 0
c clear counter array
      do 10 k = 1, mxy1
      kpic(k) = 0
   10 continue
c find how many particles in each tile
      do 20 j = 1, nop
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = m/my
      m = n + mx1*m
      if (m.le.mxy1) then
         kpic(m) = kpic(m) + 1
      else
         ierr = max(ierr,m-mxy1)
      endif
   20 continue
c find maximum
      isum = 0
      npx = 0
      do 30 k = 1, mxy1
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
      subroutine PPMOVIN2L(part,ppart,kpic,nppmx,idimp,nop,mx,my,mx1,   
     1mxy1,irc)
c this subroutine sorts particles by x,y grid in tiles of mx, my
c and copies to segmented array ppart
c linear interpolation
c input: all except ppart, kpic, output: ppart, kpic
c part/ppart = input/output particle arrays
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c ppart(3,n,m) = velocity vx of particle n in tile m
c ppart(4,n,m) = velocity vy of particle n in tile m
c kpic = output number of particles per tile
c nppmx = rmaximum number of particles in tile
c idimp = size of phase space = 4
c nop = number of particles
c mx/my = number of grids in sorting cell in x and y
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer kpic, nppmx, idimp, nop, mx, my, mx1, mxy1, irc
      real part, ppart
      dimension part(idimp,nop), ppart(idimp,nppmx,mxy1)
      dimension kpic(mxy1)
c local data
      integer i, j, k, n, m, ip, ierr
      ierr = 0
c clear counter array
      do 10 k = 1, mxy1
      kpic(k) = 0
   10 continue
c find addresses of particles at each tile and reorder particles
      do 30 j = 1, nop
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = m/my
      m = n + mx1*m
      ip = kpic(m) + 1
      if (ip.le.nppmx) then
         do 20 i = 1, idimp
         ppart(i,ip,m) = part(i,j)
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
      subroutine PPCHECK2L(ppart,kpic,idimp,nppmx,nx,ny,mx,my,mx1,my1,  
     1irc)
c this subroutine performs a sanity check to make sure particles sorted
c by x,y grid in tiles of mx, my, are all within bounds.
c tiles are assumed to be arranged in 2D linear memory
c input: all except irc
c output: irc
c ppart(1,n,k) = position x of particle n in tile k
c ppart(2,n,k) = position y of particle n in tile k
c kpic(k) = number of reordered output particles in tile k
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, ny, mx, my, mx1, my1, irc
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mx1*my1)
      dimension kpic(mx1*my1)
c local data
      integer mxy1, noff, moff, npp, j, k, ist, nn, mm
      real edgelx, edgely, edgerx, edgery, dx, dy
      mxy1 = mx1*my1
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,moff,npp,nn,mm,ist,edgelx,edgely,edgerx,edgery,
!$OMP& dx,dy)
      do 20 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
c loop over particles in tile
      do 10 j = 1, npp
      dx = ppart(1,j,k)
      dy = ppart(2,j,k)
c find particles going out of bounds
      ist = 0
      if (dx.lt.edgelx) ist = 1
      if (dx.ge.edgerx) ist = 2
      if (dy.lt.edgely) ist = ist + 3
      if (dy.ge.edgery) ist = ist + 6
      if (ist.gt.0) irc = k
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine GPPUSH2L(ppart,fxy,kpic,qbm,dt,ek,idimp,nppmx,nx,ny,mx,
     1my,nxv,nyv,mx1,mxy1,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions.
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 44 flops/particle, 12 loads, 4 stores
c input: all, output: ppart, ek
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
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c ppart(3,n,m) = velocity vx of particle n in tile m
c ppart(4,n,m) = velocity vy of particle n in tile m
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c kpic = number of particles per tile
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ipbc
      real qbm, dt, ek
      real ppart, fxy
      integer kpic
      dimension ppart(idimp,nppmx,mxy1), fxy(2,nxv,nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm
      real qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy
      real sfxy
      dimension sfxy(2,MXV,MYV)
c     dimension sfxy(2,mx+1,my+1)
      double precision sum1, sum2
      qtm = qbm*dt
      sum2 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,x,y,dxp,dyp,amx,amy,dx,dy,vx,vy
!$OMP& ,sum1,sfxy)
!$OMP& REDUCTION(+:sum2)
      do 40 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c load local fields from global array
      do 20 j = 1, min(my,ny-moff)+1
      do 10 i = 1, min(mx,nx-noff)+1
      sfxy(1,i,j) = fxy(1,i+noff,j+moff)
      sfxy(2,i,j) = fxy(2,i+noff,j+moff)
   10 continue
   20 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 30 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1
      mm = mm - moff + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find acceleration
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      vx = amx*sfxy(1,nn,mm+1)
      vy = amx*sfxy(2,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + vx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + vy)
c new velocity
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      dx = vx + qtm*dx
      dy = vy + qtm*dy
c average kinetic energy
      vx = vx + dx
      vy = vy + dy
      sum1 = sum1 + (vx*vx + vy*vy)
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
c new position
      dx = x + dx*dt
      dy = y + dy*dt
c reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
      endif
c set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   30 continue
      sum2 = sum2 + sum1
   40 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine GPPUSHF2L(ppart,fxy,kpic,ncl,ihole,qbm,dt,ek,idimp,    
     1nppmx,nx,ny,mx,my,nxv,nyv,mx1,mxy1,ntmax,irc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with periodic boundary conditions.
c also determines list of particles which are leaving this tile
c OpenMP version using guard cells
c data read in tiles
c particles stored segmented array
c 44 flops/particle, 12 loads, 4 stores
c input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
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
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c ppart(3,n,m) = velocity vx of particle n in tile m
c ppart(4,n,m) = velocity vy of particle n in tile m
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = destination of particle leaving hole
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
c optimized version
      implicit none
      integer idimp, nppmx, nx, ny, mx, my, nxv, nyv, mx1, mxy1, ntmax
      integer irc
      real qbm, dt, ek
      real ppart, fxy
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxy1), fxy(2,nxv,nyv)
      dimension kpic(mxy1), ncl(8,mxy1)
      dimension ihole(2,ntmax+1,mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, ih, nh, nn, mm
      real qtm, dxp, dyp, amx, amy
      real x, y, dx, dy, vx, vy
      real anx, any, edgelx, edgely, edgerx, edgery
      real sfxy
      dimension sfxy(2,MXV,MYV)
c     dimension sfxy(2,mx+1,my+1)
      double precision sum1, sum2
      qtm = qbm*dt
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,ih,nh,x,y,dxp,dyp,amx,amy,dx,dy
!$OMP& ,vx,vy,edgelx,edgely,edgerx,edgery,sum1,sfxy)
!$OMP& REDUCTION(+:sum2)
      do 50 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
      ih = 0
      nh = 0
c load local fields from global array
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxy(1,i,j) = fxy(1,i+noff,j+moff)
      sfxy(2,i,j) = fxy(2,i+noff,j+moff)
   10 continue
   20 continue
c clear counters
      do 30 j = 1, 8
      ncl(j,k) = 0
   30 continue
      sum1 = 0.0d0
c loop over particles in tile
      do 40 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noff + 1
      mm = mm - moff + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find acceleration
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      vx = amx*sfxy(1,nn,mm+1)
      vy = amx*sfxy(2,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + vx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + vy)
c new velocity
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      dx = vx + qtm*dx
      dy = vy + qtm*dy
c average kinetic energy
      vx = vx + dx
      vy = vy + dy
      sum1 = sum1 + (vx*vx + vy*vy)
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
c new position
      dx = x + dx*dt
      dy = y + dy*dt
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
c set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
c increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
   40 continue
      sum2 = sum2 + sum1
c set error and end of file flag
c ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   50 continue
!$OMP END PARALLEL DO
c normalize kinetic energy
      ek = ek + 0.125*sum2
      return
      end
c-----------------------------------------------------------------------
      subroutine GPPOST2L(ppart,q,kpic,qm,nppmx,idimp,mx,my,nxv,nyv,mx1,
     1mxy1)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c OpenMP version using guard cells
c data deposited in tiles
c particles stored segmented array
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c ppart(1,n,m) = position x of particle n in tile m
c ppart(2,n,m) = position y of particle n in tile m
c q(j,k) = charge density at grid point j,k
c kpic = number of particles per tile
c qm = charge on particle, in units of e
c nppmx = maximum number of particles in tile
c idimp = size of phase space = 4
c mx/my = number of grids in sorting cell in x/y
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
c mx1 = (system length in x direction - 1)/mx + 1
c mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
      implicit none
      integer nppmx, idimp, mx, my, nxv, nyv, mx1, mxy1
      real qm
      real ppart, q
      integer kpic
      dimension ppart(idimp,nppmx,mxy1), q(nxv,nyv)
      dimension kpic(mxy1)
c local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noff, moff, npp
      integer i, j, k, nn, mm
      real x, y, dxp, dyp, amx, amy
      real sq
c     dimension sq(MXV,MYV)
      dimension sq(mx+1,my+1)
c error if local array is too small
c     if ((mx.ge.MXV).or.(my.ge.MYV)) return
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noff,moff,npp,nn,mm,x,y,dxp,dyp,amx,amy,sq)
      do 80 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
c zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      sq(i,j) = 0.0
   10 continue
   20 continue
c loop over particles in tile
      do 30 j = 1, npp
c find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noff + 1
      mm = mm - moff + 1
      amx = qm - dxp
      amy = 1.0 - dyp
c deposit charge within tile to local accumulator
      x = sq(nn,mm) + amx*amy
      y = sq(nn+1,mm) + dxp*amy
      sq(nn,mm) = x
      sq(nn+1,mm) = y
      x = sq(nn,mm+1) + amx*dyp
      y = sq(nn+1,mm+1) + dxp*dyp
      sq(nn,mm+1) = x
      sq(nn+1,mm+1) = y
   30 continue
c deposit charge to interior points in global array
      nn = min(mx,nxv-noff)
      mm = min(my,nyv-moff)
      do 50 j = 2, mm
      do 40 i = 2, nn
      q(i+noff,j+moff) = q(i+noff,j+moff) + sq(i,j)
   40 continue
   50 continue
c deposit charge to edge points in global array
      mm = min(my+1,nyv-moff)
      do 60 i = 2, nn
!$OMP ATOMIC
      q(i+noff,1+moff) = q(i+noff,1+moff) + sq(i,1)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noff,mm+moff) = q(i+noff,mm+moff) + sq(i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noff)
      do 70 j = 1, mm
!$OMP ATOMIC
      q(1+noff,j+moff) = q(1+noff,j+moff) + sq(1,j)
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noff,j+moff) = q(nn+noff,j+moff) + sq(nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPORDER2L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx,ny
     1,mx,my,mx1,my1,npbmx,ntmax,irc)
c this subroutine sorts particles by x,y grid in tiles of mx, my
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 2D linear memory
c algorithm has 3 steps.  first, one finds particles leaving tile and
c stores their number in each directon, location, and destination in ncl
c and ihole.  second, a prefix scan of ncl is performed and departing
c particles are buffered in ppbuff in direction order.  finally, we copy
c the incoming particles from other tiles into ppart.
c input: all except ppbuff, ncl, ihole, irc
c output: ppart, ppbuff, kpic, ncl, ihole, irc
c ppart(1,n,k) = position x of particle n in tile k
c ppart(2,n,k) = position y of particle n in tile k
c ppbuff(i,n,k) = i co-ordinate of particle n in tile k
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = direction destination of particle leaving hole
c all for tile k
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c nx/ny = system length in x/y direction
c mx/my = number of grids in sorting cell in x/y
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, nx, ny, mx, my, mx1, my1, npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1*my1), ppbuff(idimp,npbmx,mx1*my1)
      dimension kpic(mx1*my1), ncl(8,mx1*my1)
      dimension ihole(2,ntmax+1,mx1*my1)
c local data
      integer mxy1, noff, moff, npp, ncoff
      integer i, j, k, ii, kx, ky, ih, nh, ist, nn, mm, isum
      integer ip, j1, j2, kxl, kxr, kk, kl, kr
      real anx, any, edgelx, edgely, edgerx, edgery, dx, dy
      integer ks
      dimension ks(8)
      mxy1 = mx1*my1
      anx = real(nx)
      any = real(ny)
c find and count particles leaving tiles and determine destination
c update ppart, ihole, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noff,moff,npp,nn,mm,ih,nh,ist,dx,dy,edgelx,edgely,
!$OMP& edgerx,edgery)
      do 30 k = 1, mxy1
      noff = (k - 1)/mx1
      moff = my*noff
      noff = mx*(k - mx1*noff - 1)
      npp = kpic(k)
      nn = min(mx,nx-noff)
      mm = min(my,ny-moff)
      ih = 0
      nh = 0
      edgelx = noff
      edgerx = noff + nn
      edgely = moff
      edgery = moff + mm
c clear counters
      do 10 j = 1, 8
      ncl(j,k) = 0
   10 continue
c loop over particles in tile
      do 20 j = 1, npp
      dx = ppart(1,j,k)
      dy = ppart(2,j,k)
c find particles going out of bounds
      ist = 0
c count how many particles are going in each direction in ncl
c save their address and destination in ihole
c use periodic boundary conditions and check for roundoff error
c ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(1,j,k) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(1,j,k) = dx
         else
            ist = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(2,j,k) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(2,j,k) = dy
         else
            ist = ist + 3
         endif
      endif
      if (ist.gt.0) then
         ncl(ist,k) = ncl(ist,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = ist
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
      ihole(1,1,k) = ih
   30 continue
!$OMP END PARALLEL DO
c ihole overflow
      if (irc.gt.0) return
c
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 70 k = 1, mxy1
c find address offset for ordered ppbuff array
      isum = 0
      do 40 j = 1, 8
      ist = ncl(j,k)
      ncl(j,k) = isum
      isum = isum + ist
   40 continue
      nh = ihole(1,1,k)
      ip = 0
c loop over particles leaving tile
      do 60 j = 1, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = ncl(ist,k) + 1
      if (ii.le.npbmx) then
         do 50 i = 1, idimp
         ppbuff(i,ii,k) = ppart(i,j1,k)
   50    continue
      else
         ip = 1
      endif
      ncl(ist,k) = ii
   60 continue
c set error
      if (ip.gt.0) irc = ncl(8,k)
   70 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,kk,npp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,ist,j1,j2
!$OMP& ,ip,ks)
      do 130 k = 1, mxy1
      npp = kpic(k)
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
      ks(1) = kxr + kk
      ks(2) = kxl + kk
      ks(3) = kx + kr
      ks(4) = kxr + kr
      ks(5) = kxl + kr
      ks(6) = kx + kl
      ks(7) = kxr + kl
      ks(8) = kxl + kl
c loop over directions
      nh = ihole(1,1,k)
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 100 ii = 1, 8
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
c ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
      do 90 j = 1, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,k)
c place overflow at end of array
      else
         j1 = npp + 1
         npp = j1
      endif
      if (j1.le.nppmx) then
         do 80 i = 1, idimp
         ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
   80    continue
      else
         ist = 1
      endif
   90 continue
  100 continue
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
      if (ih.lt.nh) then
         ip = nh - ih
         do 120 j = 1, ip
         j1 = npp - j + 1
         j2 = ihole(1,nh-j+2,k)
         if (j1.gt.j2) then
c move particle only if it is below current hole
            do 110 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
  110       continue
         endif
  120    continue
         npp = npp - ip
      endif
      kpic(k) = npp
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine PPORDERF2L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,mx1,
     1my1,npbmx,ntmax,irc)
c this subroutine sorts particles by x,y grid in tiles of mx, my
c linear interpolation, with periodic boundary conditions
c tiles are assumed to be arranged in 2D linear memory.
c the algorithm has 2 steps.  first, a prefix scan of ncl is performed
c and departing particles are buffered in ppbuff in direction order.
c then we copy the incoming particles from other tiles into ppart.
c it assumes that the number, location, and destination of particles 
c leaving a tile have been previously stored in ncl and ihole by the
c GPPUSHF2L subroutine.
c input: all except ppbuff, irc
c output: ppart, ppbuff, kpic, ncl, irc
c ppart(1,n,k) = position x of particle n in tile k
c ppart(2,n,k) = position y of particle n in tile k
c ppbuff(i,n,k) = i co-ordinate of particle n in tile k
c kpic(k) = number of particles in tile k
c ncl(i,k) = number of particles going to destination i, tile k
c ihole(1,:,k) = location of hole in array left by departing particle
c ihole(2,:,k) = direction destination of particle leaving hole
c all for tile k
c ihole(1,1,k) = ih, number of holes left (error, if negative)
c idimp = size of phase space = 4
c nppmx = maximum number of particles in tile
c mx1 = (system length in x direction - 1)/mx + 1
c my1 = (system length in y direction - 1)/my + 1
c npbmx = size of buffer array ppbuff
c ntmax = size of hole array for particles leaving tiles
c irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, my1, npbmx, ntmax, irc
      real ppart, ppbuff
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mx1*my1), ppbuff(idimp,npbmx,mx1*my1)
      dimension kpic(mx1*my1), ncl(8,mx1*my1)
      dimension ihole(2,ntmax+1,mx1*my1)
c local data
      integer mxy1, npp, ncoff
      integer i, j, k, ii, kx, ky, ih, nh, ist, isum
      integer ip, j1, j2, kxl, kxr, kk, kl, kr
      integer ks
      dimension ks(8)
      mxy1 = mx1*my1
c buffer particles that are leaving tile: update ppbuff, ncl
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 40 k = 1, mxy1
c find address offset for ordered ppbuff array
      isum = 0
      do 10 j = 1, 8
      ist = ncl(j,k)
      ncl(j,k) = isum
      isum = isum + ist
   10 continue
      nh = ihole(1,1,k)
      ip = 0
c loop over particles leaving tile
      do 30 j = 1, nh
c buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = ncl(ist,k) + 1
      if (ii.le.npbmx) then
         do 20 i = 1, idimp
         ppbuff(i,ii,k) = ppart(i,j1,k)
   20    continue
      else
         ip = 1
      endif
      ncl(ist,k) = ii
   30 continue
c set error
      if (ip.gt.0) irc = ncl(8,k)
   40 continue
!$OMP END PARALLEL DO
c ppbuff overflow
      if (irc.gt.0) return
c
c copy incoming particles from buffer into ppart: update ppart, kpic
c loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,kk,npp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,ist,j1,j2
!$OMP& ,ip,ks)
      do 100 k = 1, mxy1
      npp = kpic(k)
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
      ks(1) = kxr + kk
      ks(2) = kxl + kk
      ks(3) = kx + kr
      ks(4) = kxr + kr
      ks(5) = kxl + kr
      ks(6) = kx + kl
      ks(7) = kxr + kl
      ks(8) = kxl + kl
c loop over directions
      nh = ihole(1,1,k)
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 70 ii = 1, 8
      if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
c ip = number of particles coming from direction ii
      ip = ncl(ii,ks(ii)) - ncoff
      do 60 j = 1, ip
      ih = ih + 1
c insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,k)
c place overflow at end of array
      else
         j1 = npp + 1
         npp = j1
      endif
      if (j1.le.nppmx) then
         do 50 i = 1, idimp
         ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
   50    continue
      else
         ist = 1
      endif
   60 continue
   70 continue
c set error
      if (ist.gt.0) irc = j1
c fill up remaining holes in particle array with particles from bottom
      if (ih.lt.nh) then
         ip = nh - ih
         do 90 j = 1, ip
         j1 = npp - j + 1
         j2 = ihole(1,nh-j+2,k)
         if (j1.gt.j2) then
c move particle only if it is below current hole
            do 80 i = 1, idimp
            ppart(i,j2,k) = ppart(i,j1,k)
   80       continue
         endif
   90    continue
         npp = npp - ip
      endif
      kpic(k) = npp
  100 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine CGUARD2L(fxy,nx,ny,nxe,nye)
c replicate extended periodic vector field fxy
c linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nxe = second dimension of field arrays, must be >= ny+1
      implicit none
      real fxy
      integer nx, ny, nxe, nye
      dimension fxy(2,nxe,nye)
c local data
      integer j, k
c copy edges of extended field
      do 10 k = 1, ny
      fxy(1,nx+1,k) = fxy(1,1,k)
      fxy(2,nx+1,k) = fxy(2,1,k)
   10 continue
      do 20 j = 1, nx
      fxy(1,j,ny+1) = fxy(1,j,1)
      fxy(2,j,ny+1) = fxy(2,j,1)
   20 continue
      fxy(1,nx+1,ny+1) = fxy(1,1,1)
      fxy(2,nx+1,ny+1) = fxy(2,1,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine AGUARD2L(q,nx,ny,nxe,nye)
c accumulate extended periodic scalar field q
c linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nxe = second dimension of field arrays, must be >= ny+1
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k
c accumulate edges of extended field
      do 10 k = 1, ny
      q(1,k) = q(1,k) + q(nx+1,k)
      q(nx+1,k) = 0.0
   10 continue
      do 20 j = 1, nx
      q(j,1) = q(j,1) + q(j,ny+1)
      q(j,ny+1) = 0.0
   20 continue
      q(1,1) = q(1,1) + q(nx+1,ny+1)
      q(nx+1,ny+1) = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine MPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,  
     1nxhd,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,nxvh,nyhd, output: fxy,we
c approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, we
      complex q, fxy, ffc
      dimension q(nxvh,nyv), fxy(2,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2
      double precision wp, sum1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*real(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*real(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffc(j,k) = cmplx(affp,1.0)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 sum1 = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dky,at1,at2,at3,zt1,zt2,wp)
!$OMP& REDUCTION(+:sum1)
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      wp = 0.0d0
      do 40 j = 2, nxh
      at1 = real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      wp = wp + at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
   40 continue
c mode numbers kx = 0, nx/2
      at1 = real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dky*at1
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at3*zt1
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nxh
      at1 = real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*real(j - 1)*at1
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      fxy(1,j,1) = at2*zt1
      fxy(2,j,1) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   60 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      sum1 = sum1 + wp
      we = real(nx*ny)*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
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
      subroutine WFFT2RMX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,   
     1nxyhd)
c wrapper function for real to complex fft, with packed data
c parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
c perform y fft
         call FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c perform x fft
         call FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT2RM2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,   
     1nxyhd)
c wrapper function for 2 2d real to complex ffts, with packed data
c parallelized with OpenMP
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2RM2X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
c perform y fft
         call FFT2RM2Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RM2Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,    
     1nxhyd,nxyhd)
c perform x fft
         call FFT2RM2X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd
     1,nxyhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RMXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd, 
     1nxhyd,nxyhd)
c this subroutine performs the x part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic, with OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform in x is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, a forward fourier transform in x is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2
c nyd = second dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c imag(f(1,1)) = real part of mode nx/2,0 and
c imag(f(1,ny/2+1) ) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, nrxb
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 70
c inverse fourier transform
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,ani,t1,t2,t3)
      do 60 i = nyi, nyt
c bit-reverse array elements in x
      do 10 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(j1,i)
         f(j1,i) = f(j,i)
         f(j,i) = t1
      endif
   10 continue
c then transform in x
      do 40 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
   20 continue
   30 continue
   40 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 0.5/(real(nx)*real(ny))
      do 50 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      t2 = conjg(f(nxh2-j,i))
      t1 = f(j,i) + t2
      t2 = (f(j,i) - t2)*t3
      f(j,i) = ani*(t1 + t2)
      f(nxh2-j,i) = ani*conjg(t1 - t2)
   50 continue
      ani = 2.0*ani
      f(nxhh+1,i) = ani*conjg(f(nxhh+1,i))
      f(1,i) = ani*cmplx(real(f(1,i)) + aimag(f(1,i)),                  
     1                   real(f(1,i)) - aimag(f(1,i)))
   60 continue
!$OMP END PARALLEL DO
      return
c forward fourier transform
   70 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,t1,t2,t3)
      do 130 i = nyi, nyt
c scramble coefficients
      kmr = nxy/nx
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      t2 = conjg(f(nxh2-j,i))
      t1 = f(j,i) + t2
      t2 = (f(j,i) - t2)*t3
      f(j,i) = t1 + t2
      f(nxh2-j,i) = conjg(t1 - t2)
   80 continue
      f(nxhh+1,i) = 2.0*conjg(f(nxhh+1,i))
      f(1,i) = cmplx(real(f(1,i)) + aimag(f(1,i)),                      
     1               real(f(1,i)) - aimag(f(1,i)))
c bit-reverse array elements in x
      do 90 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(j1,i)
         f(j1,i) = f(j,i)
         f(j,i) = t1
      endif
   90 continue
c then transform in x
      do 120 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
  100 continue
  110 continue
  120 continue
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RMXY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd, 
     1nxhyd,nxyhd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic, with OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform in y is performed
c f(n,m) = sum(f(j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform in y is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = first dimension of f >= nx/2
c nyd = second dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c imag(f(1,1)) = real part of mode nx/2,0 and
c imag(f(1,ny/2+1) ) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, nryb
      complex t1, t2
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 70
c inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,t1,t2)
      do 50 i = nxi, nxt
c bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(i,k1)
         f(i,k1) = f(i,k)
         f(i,k) = t1
      endif
   10 continue
c then transform in y
      do 40 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 60 k = 2, nyh
         t1 = f(1,ny2-k)
         f(1,ny2-k) = 0.5*cmplx(aimag(f(1,k) + t1),real(f(1,k) - t1))
         f(1,k) = 0.5*cmplx(real(f(1,k) + t1),aimag(f(1,k) - t1))
   60    continue
      endif
      return
c forward fourier transform
   70 nryb = nxhy/ny
      nry = nxy/ny
c scramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 80 k = 2, nyh
         t1 = cmplx(aimag(f(1,ny2-k)),real(f(1,ny2-k)))
         f(1,ny2-k) = conjg(f(1,k) - t1)
         f(1,k) = f(1,k) + t1
   80    continue
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,t1,t2)
      do 130 i = nxi, nxt
c bit-reverse array elements in y
      do 90 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(i,k1)
         f(i,k1) = f(i,k)
         f(i,k) = t1
      endif
   90 continue
c then transform in y
      do 120 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
  100 continue
  110 continue
  120 continue
  130 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RM2X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd, 
     1nxhyd,nxyhd)
c this subroutine performs the x part of 2 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c y, using complex arithmetic, with OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, two inverse fourier transforms in x are performed
c f(1:2,n,m) = (1/nx*ny)*sum(f(1:2,j,k)*exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, two forward fourier transforms in x are performed
c f(1:2,j,k) = sum(f(1:2,n,m)*exp(sqrt(-1)*2pi*n*j/nx))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = second dimension of f >= nx/2
c nyd = third dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(1:2,j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1:2,1,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c imag(f(1:2,1,1)) = real part of mode nx/2,0 and
c imag(f(1:2,1,ny/2+1) ) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      integer nrxb
      real at1, ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 100
c inverse fourier transform
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,ani,t1,t2,t3)
      do 90 i = nyi, nyt
c swap complex components
      do 10 j = 1, nxh
      at1 = aimag(f(1,j,i))
      f(1,j,i) = cmplx(real(f(1,j,i)),real(f(2,j,i)))
      f(2,j,i) = cmplx(at1,aimag(f(2,j,i)))
   10 continue
c bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i)
         t2 = f(2,j1,i)
         f(1,j1,i) = f(1,j,i)
         f(2,j1,i) = f(2,j,i)
         f(1,j,i) = t1
         f(2,j,i) = t2
      endif
   20 continue
c then transform in x
      do 50 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
   30 continue
   40 continue
   50 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 0.5/(real(nx)*real(ny))
      do 70 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 60 jj = 1, 2
      t2 = conjg(f(jj,nxh2-j,i))
      t1 = f(jj,j,i) + t2
      t2 = (f(jj,j,i) - t2)*t3
      f(jj,j,i) = ani*(t1 + t2)
      f(jj,nxh2-j,i) = ani*conjg(t1 - t2)
   60 continue
   70 continue
      ani = 2.0*ani
      do 80 jj = 1, 2
      f(jj,nxhh+1,i) = ani*conjg(f(jj,nxhh+1,i))
      f(jj,1,i) = ani*cmplx(real(f(jj,1,i)) + aimag(f(jj,1,i)),         
     1                      real(f(jj,1,i)) - aimag(f(jj,1,i)))
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
c forward fourier transform
  100 nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,at1,t1,t2,t3)
      do 190 i = nyi, nyt
c scramble coefficients
      kmr = nxy/nx
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 jj = 1, 2
      t2 = conjg(f(jj,nxh2-j,i))
      t1 = f(jj,j,i) + t2
      t2 = (f(jj,j,i) - t2)*t3
      f(jj,j,i) = t1 + t2
      f(jj,nxh2-j,i) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 jj = 1, 2
      f(jj,nxhh+1,i) = 2.0*conjg(f(jj,nxhh+1,i))
      f(jj,1,i) = cmplx(real(f(jj,1,i)) + aimag(f(jj,1,i)),             
     1                  real(f(jj,1,i)) - aimag(f(jj,1,i)))
  130 continue
c bit-reverse array elements in x
      do 140 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t1 = f(1,j1,i)
         t2 = f(2,j1,i)
         f(1,j1,i) = f(1,j,i)
         f(2,j1,i) = f(2,j,i)
         f(1,j,i) = t1
         f(2,j,i) = t2
      endif
  140 continue
c then transform in x
      do 170 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 160 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 150 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
  150 continue
  160 continue
  170 continue
c swap complex components
      do 180 j = 1, nxh
      at1 = aimag(f(1,j,i))
      f(1,j,i) = cmplx(real(f(1,j,i)),real(f(2,j,i)))
      f(2,j,i) = cmplx(at1,aimag(f(2,j,i)))
  180 continue
  190 continue
!$OMP END PARALLEL DO
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT2RM2Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd, 
     1nxhyd,nxyhd)
c this subroutine performs the y part of 2 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c x, using complex arithmetic, with OpenMP
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, two inverse fourier transforms in y are performed
c f(1:2,n,m) = sum(f(1:2,j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms in y are performed
c f(1:2,j,k) = sum(f(1:2,n,m)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = second dimension of f >= nx/2
c nyd = third dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(1:2,j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1:2,1,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c imag(f(1:2,1,1)) = real part of mode nx/2,0 and
c imag(f(1:2,1,ny/2+1) ) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      integer nryb
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 80
c inverse fourier transform
      nryb = nxhy/ny
      nry = nxy/ny
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,t1,t2,t3)
      do 50 i = nxi, nxt
c bit-reverse array elements in y
      do 10 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(1,i,k1)
         t2 = f(2,i,k1)
         f(1,i,k1) = f(1,i,k)
         f(2,i,k1) = f(2,i,k)
         f(1,i,k) = t1
         f(2,i,k) = t2
      endif
   10 continue
c then transform in y
      do 40 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 30 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 20 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
   20 continue
   30 continue
   40 continue
   50 continue
!$OMP END PARALLEL DO
c unscramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 70 k = 2, nyh
         do 60 jj = 1, 2
         t1 = f(jj,1,ny2-k)
         f(jj,1,ny2-k) = 0.5*cmplx(aimag(f(jj,1,k) + t1),               
     1                             real(f(jj,1,k) - t1))
         f(jj,1,k) = 0.5*cmplx(real(f(jj,1,k) + t1),                    
     1                   aimag(f(jj,1,k) - t1))
   60    continue
   70    continue
      endif
      return
c forward fourier transform
   80 nryb = nxhy/ny
      nry = nxy/ny
c scramble modes kx = 0, nx/2
      if (nxi.eq.1) then
         do 100 k = 2, nyh
         do 90 jj = 1, 2
         t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
         f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
         f(jj,1,k) = f(jj,1,k) + t1
   90    continue
  100    continue
      endif
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,t1,t2,t3)
      do 150 i = nxi, nxt
c bit-reverse array elements in y
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t1 = f(1,i,k1)
         t2 = f(2,i,k1)
         f(1,i,k1) = f(1,i,k)
         f(2,i,k1) = f(2,i,k)
         f(1,i,k) = t1
         f(2,i,k) = t2
      endif
  110 continue
c then transform in y
      do 140 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
  120 continue
  130 continue
  140 continue
  150 continue
!$OMP END PARALLEL DO
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
