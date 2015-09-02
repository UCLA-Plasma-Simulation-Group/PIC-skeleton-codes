/* C Library for Skeleton 2D Electrostatic MPI PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "ppush2.h"
#include "pplib2.h"

/*--------------------------------------------------------------------*/
double ranorm() {
/* this program calculates a random number y from a gaussian distribution
   with zero mean and unit variance, according to the method of
   mueller and box:
      y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
      y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
   where x is a random number uniformly distributed on (0,1).
   written for the ibm by viktor k. decyk, ucla
local data                                                              */
   static int r1 = 885098780, r2 = 1824280461;
   static int r4 = 1396483093, r5 = 55318673;
   static int iflg = 0;
   static double h1l = 65531.0, h1u = 32767.0, h2l = 65525.0;
   static double r0 = 0.0;
   int isc, i1;
   double ranorm, r3, asc, bsc, temp;
   if (iflg==1) {
      ranorm = r0;
      r0 = 0.0;
      iflg = 0;
      return ranorm;
   }
   isc = 65536;
   asc = (double) isc;
   bsc = asc*asc;
   i1 = r1 - (r1/isc)*isc;
   r3 = h1l*(double) r1 + asc*h1u*(double) i1;
   i1 = r3/bsc;
   r3 -= ((double) i1)*bsc;
   bsc = 0.5*bsc;
   i1 = r2/isc;
   isc = r2 - i1*isc;
   r0 = h1l*(double) r2 + asc*h1u*(double) isc;
   asc = 1.0/bsc;
   isc = r0*asc;
   r2 = r0 - ((double) isc)*bsc;
   r3 += (double) isc + 2.0*h1u*(double) i1;
   isc = r3*asc;
   r1 = r3 - ((double) isc)*bsc;
   temp = sqrt(-2.0*log((((double) r1) + ((double) r2)*asc)*asc));
   isc = 65536;
   asc = (double) isc;
   bsc = asc*asc;
   i1 = r4 - (r4/isc)*isc;
   r3 = h2l*(double) r4 + asc*h1u*(double) i1;
   i1 = r3/bsc;
   r3 -= ((double) i1)*bsc;
   bsc = 0.5*bsc;
   i1 = r5/isc;
   isc = r5 - i1*isc;
   r0 = h2l*(double) r5 + asc*h1u*(double) isc;
   asc = 1.0/bsc;
   isc = r0*asc;
   r5 = r0 - ((double) isc)*bsc;
   r3 += (double) isc + 2.0*h1u*(double) i1;
   isc = r3*asc;
   r4 = r3 - ((double) isc)*bsc;
   r0 = 6.28318530717959*((((double) r4) + ((double) r5)*asc)*asc);
   ranorm = temp*sin(r0);
   r0 = temp*cos(r0);
   iflg = 1;
   return ranorm;
}

/*--------------------------------------------------------------------*/
void cpdicomp2l(float edges[], int *nyp, int *noff, int *nypmx,
                int *nypmn, int ny, int kstrt, int nvp, int idps) {
/* this subroutine determines spatial boundaries for uniform particle
   decomposition, calculates number of grid points in each spatial
   region, and the offset of these grid points from the global address
   nvp must be < ny.  some combinations of ny and nvp result in a zero
   value of nyp.  this is not supported.
   integer boundaries are set.
   input: ny, kstrt, nvp, idps, output: edges, nyp, noff, nypmx, nypmn
   edges[0] = lower boundary of particle partition
   edges[1] = upper boundary of particle partition
   nyp = number of primary (complete) gridpoints in particle partition
   noff = lowermost global gridpoint in particle partition
   nypmx = maximum size of particle partition, including guard cells
   nypmn = minimum value of nyp
   ny = system length in y direction
   kstrt = starting data block number (processor id + 1)
   nvp = number of real or virtual processors
   idps = number of partition boundaries
local data                                                            */
   int kb, kyp;
   float at1, any;
   int mypm[2], iwork2[2];
   any = (float) ny;
/* determine decomposition */
   kb = kstrt - 1;
   kyp = (ny - 1)/nvp + 1;
   at1 = (float) kyp;
   edges[0] = at1*(float) kb;
   if (edges[0] > any)
      edges[0] = any;
   *noff = edges[0];
   edges[1] = at1*(float) (kb + 1);
   if (edges[1] > any)
      edges[1] = any;
   kb = edges[1];
   *nyp = kb - *noff;
/* find maximum/minimum partition size */
   mypm[0] = *nyp;
   mypm[1] = -(*nyp);
   cppimax(mypm,iwork2,2);
   *nypmx = mypm[0] + 1;
   *nypmn = -mypm[1];
   return;
}

/*--------------------------------------------------------------------*/
void cpdistr2(float part[], float edges[], int *npp, int nps, float vtx,
              float vty, float vdx, float vdy, int npx, int npy, int nx,
              int ny, int idimp, int npmax, int idps, int ipbc, int *ierr) {
/* for 2d code, this subroutine calculates initial particle co-ordinates
   and velocities with uniform density and maxwellian velocity with drift
   for distributed data.
   input: all except part, npp, ierr, output: part, npp, ierr
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   part[n][2] = velocity vx of particle n in partition
   part[n][3] = velocity vy of particle n in partition
   edges[0] = lower boundary of particle partition
   edges[1] = upper boundary of particle partition
   npp = number of particles in partition
   nps = starting address of particles in partition
   vtx/vty = thermal velocity of electrons in x/y direction
   vdx/vdy = drift velocity of beam electrons in x/y direction
   npx/npy = initial number of particles distributed in x/y direction
   nx/ny = system length in x/y direction
   idimp = size of phase space = 4
   npmax = maximum number of particles in each partition
   idps = number of partition boundaries
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   ierr = (0,1) = (no,yes) error condition exists
   ranorm = gaussian random number with zero mean and unit variance
   with spatial decomposition
local data                                                            */
   int j, k, npt, k1, npxyp;
   float edgelx, edgely, at1, at2, xt, yt, vxt, vyt;
   double dnpx, dnpxy, dt1;
   int ierr1[1], iwork1[1];
   double sum3[3], work3[3];
   *ierr = 0;
/* particle distribution constant */
   dnpx = (double) npx;
/* set boundary values */
   edgelx = 0.0;
   edgely = 0.0;
   at1 = (float) nx/(float) npx;
   at2 = (float) ny/(float) npy;
   if (ipbc==2) {
      edgelx = 1.0;
      edgely = 1.0;
      at1 = (float) (nx-2)/(float) npx;
      at2 = (float) (ny-2)/(float) npy;
   }
   else if (ipbc==3) {
      edgelx = 1.0;
      at1 = (float) (nx-2)/(float) npx;
   }
   npt = *npp;
/* uniform density profile */
   for (k = 0; k < npy; k++) {
      yt = edgely + at2*(((float) k) + 0.5);
      for (j = 0; j < npx; j++) {
         xt = edgelx + at1*(((float) j) + 0.5);
/* maxwellian velocity distribution */
         vxt = vtx*ranorm();
         vyt = vty*ranorm();
         if ((yt >= edges[0]) && (yt < edges[1])) {
            if (npt < npmax) {
               k1 = idimp*npt;
               part[k1] = xt;
               part[1+k1] = yt;
               part[2+k1] = vxt;
               part[3+k1] = vyt;
               npt += 1;
            }
            else
               *ierr += 1;
         }
      }
   }
   npxyp = 0;
/* add correct drift */
   sum3[0] = 0.0;
   sum3[1] = 0.0;
   for (j = nps-1; j < npt; j++) {
      npxyp += 1;
      sum3[0] += part[2+idimp*j];
      sum3[1] += part[3+idimp*j];
   }
   sum3[2] = npxyp;
   cppdsum(sum3,work3,3);
   dnpxy = sum3[2];
   ierr1[0] = *ierr;
   cppimax(ierr1,iwork1,1);
   *ierr = ierr1[0];
   dt1 = 1.0/dnpxy;
   sum3[0] = dt1*sum3[0] - vdx;
   sum3[1] = dt1*sum3[1] - vdy;
   for (j = nps-1; j < npt; j++) {
      part[2+idimp*j] -= sum3[0];
      part[3+idimp*j] -= sum3[1];
   }
/* process errors */
   dnpxy -= dnpx*(double) npy;
   if (dnpxy != 0.0)
      *ierr = dnpxy;
   *npp = npt;
   return;
}

/*--------------------------------------------------------------------*/
void cppgpush2l(float part[], float fxy[], float edges[], int npp,
                int noff, int ihole[], float qbm, float dt, float *ek,
                int nx, int ny, int idimp, int npmax, int nxv,
                int nypmx, int idps, int ntmax, int ipbc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions
   also determines list of particles which are leaving this processor
   scalar version using guard cells, for distributed data
   42 flops/particle, 12 loads, 4 stores
   input: all except ihole, output: part, ihole, ek
   equations used are:
   vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
   vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
   where q/m is charge/mass, and
   x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
   fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
   the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
      + dx*fy(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   part[n][2] = velocity vx of particle n in partition
   part[n][3] = velocity vy of particle n in partition
   fxy[k][j][0] = x component of force/charge at grid (j,kk)
   fxy[k][j][1] = y component of force/charge at grid (j,kk)
   in other words, fxy are the convolutions of the electric field
   over the particle shape, where kk = k + noff
   edges[0:1] = lower:upper boundary of particle partition
   npp = number of particles in partition
   noff = lowermost global gridpoint in particle partition.
   ihole = location of hole left in particle arrays
   ihole[0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   nx/ny = system length in x/y direction
   idimp = size of phase space = 4
   npmax = maximum number of particles in each partition
   nxv = first dimension of field array, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   idps = number of partition boundaries
   ntmax = size of hole array for particles leaving processors
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int mnoff, j, nn, mm, np, mp, ih, nh, nxv2;
   float qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, vx, vy;
   double sum1;
   nxv2 = 2*nxv;
   qtm = qbm*dt;
   sum1 = 0.0;
/* set boundary values */
   edgelx = 0.0;
   edgely = 1.0;
   edgerx = (float) nx;
   edgery = (float) (ny-1);
   if ((ipbc==2) || (ipbc==3)) {
      edgelx = 1.0;
      edgerx = (float) (nx-1);
   }
   mnoff = noff;
   ih = 0;
   nh = 0;
   for (j = 0; j < npp; j++) {
/* find interpolation weights */
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = part[idimp*j] - (float) nn;
      dyp = part[1+idimp*j] - (float) mm;
      nn = 2*nn;
      mm = nxv2*(mm - mnoff);
      amx = 1.0 - dxp;
      mp = mm + nxv2;
      amy = 1.0 - dyp;
      np = nn + 2;
/* find acceleration */
      dx = dyp*(dxp*fxy[np+mp] + amx*fxy[nn+mp])
         + amy*(dxp*fxy[np+mm] + amx*fxy[nn+mm]);
      dy = dyp*(dxp*fxy[1+np+mp] + amx*fxy[1+nn+mp])
         + amy*(dxp*fxy[1+np+mm] + amx*fxy[1+nn+mm]);
/* new velocity */
      vx = part[2+idimp*j];
      vy = part[3+idimp*j];
      dx = vx + qtm*dx;
      dy = vy + qtm*dy;
/* average kinetic energy */
      vx += dx;
      vy += dy;
      sum1 += vx*vx + vy*vy;
      part[2+idimp*j] = dx;
      part[3+idimp*j] = dy;
/* new position */
      dx = part[idimp*j] + dx*dt;
      dy = part[1+idimp*j] + dy*dt;
/* periodic boundary conditions in x */
      if (ipbc==1) {
         if (dx < edgelx) dx += edgerx;
         if (dx >= edgerx) dx -= edgerx;
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = part[idimp*j];
            part[2+idimp*j] = -part[2+idimp*j];
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = part[1+idimp*j];
            part[3+idimp*j] = -part[3+idimp*j];
         }
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = part[idimp*j];
            part[2+idimp*j] = -part[2+idimp*j];
         }
      }
/* find particles out of bounds */
      if ((dy < edges[0]) || (dy >= edges[1])) {
         if (ih < ntmax)
            ihole[ih+1] = j + 1;
         else
            nh = 1;
         ih += 1;
      }
/* set new position */
      part[idimp*j] = dx;
      part[1+idimp*j] = dy;
   }
/* set end of file flag */
      if (nh > 0)
         ih = -ih;
      ihole[0] = ih;
/* normalize kinetic energy */
   *ek += 0.125*sum1;
   return;
}

/*--------------------------------------------------------------------*/
void cppgpost2l(float part[], float q[], int npp, int noff, float qm,
                int idimp, int npmax, int nxv, int nypmx) {
/* for 2d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   scalar version using guard cells, for distributed data
   17 flops/particle, 6 loads, 4 stores
   input: all, output: q
   charge density is approximated by values at the nearest grid points
   q(n,m)=qm*(1.-dx)*(1.-dy)
   q(n+1,m)=qm*dx*(1.-dy)
   q(n,m+1)=qm*(1.-dx)*dy
   q(n+1,m+1)=qm*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   q[k][j] = charge density at grid point (j,kk),
   where kk = k + noff
   npp = number of particles in partition
   noff = lowermost global gridpoint in particle partition.
   qm = charge on particle, in units of e
   idimp = size of phase space = 4
   npmax = maximum number of particles in each partition
   nxv = first dimension of charge array, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
local data                                                            */
   int  mnoff, j, nn, np, mm, mp;
   float dxp, dyp, amx, amy;
   mnoff = noff;
   for (j = 0; j < npp; j++) {
/* find interpolation weights */
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = qm*(part[idimp*j] - (float) nn);
      dyp = part[1+idimp*j] - (float) mm;
      mm = nxv*(mm - mnoff);
      amx = qm - dxp;
      mp = mm + nxv;
      amy = 1.0 - dyp;
      np = nn + 1;
/* deposit charge */
      q[np+mp] += dxp*dyp;
      q[nn+mp] += amx*dyp;
      q[np+mm] += dxp*amy;
      q[nn+mm] += amx*amy;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppdsortp2yl(float parta[], float partb[], int npic[], int npp,
                  int noff, int nyp, int idimp, int npmax, int nypm1) {
/* this subroutine sorts particles by y grid
   linear interpolation, spatial decomposition in y direction
   parta/partb = input/output particle array
   part[n][1] = position y of particle n in partition
   npic = address offset for reordering particles
   npp = number of particles in partition
   noff = backmost global gridpoint in particle partition
   nyp = number of primary gridpoints in particle partition
   idimp = size of phase space
   npmax = maximum number of particles in each partition
   nypm1 = maximum size of particle partition plus one
local data                                                            */
   int i, j, k, m, mnoff, nyp1, isum, ist, ip;
   mnoff = noff;
   nyp1 = nyp + 1;
/* clear counter array */
   for (k = 0; k < nyp1; k++) {
      npic[k] = 0;
   }
/* find how many particles in each grid */
   for (j = 0; j < npp; j++) {
      m = parta[1+idimp*j];
      m -= mnoff;
      npic[m] += 1;
   }
/* find address offset */
   isum = 0;
   for (k = 0; k < nyp1; k++) {
      ist = npic[k];
      npic[k] = isum;
      isum += ist;
   }
/* find addresses of particles at each grid and reorder particles */
   for (j = 0; j < npp; j++) {
      m = parta[1+idimp*j];
      m -= mnoff;
      ip = npic[m];
      for (i = 0; i < idimp; i++) {
         partb[i+idimp*ip] = parta[i+idimp*j];
      }
      npic[m] = ip + 1;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppcguard2xl(float fxy[], int nyp, int nx, int ndim, int nxe,
                  int nypmx) {
/* replicate extended periodic vector field in x direction
   linear interpolation, for distributed data
   nyp = number of primary (complete) gridpoints in particle partition
   nx = system length in x direction
   ndim = leading dimension of array fxy
   nxe = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells
local data                                                 */
   int i, k, kk, myp1;
/* replicate edges of extended field */
   myp1 = nyp + 1;
   for (k = 0; k < myp1; k++) {
      kk = ndim*nxe*k;
      for (i = 0; i < ndim; i++) {
         fxy[i+ndim*nx+kk] = fxy[i+kk];
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppaguard2xl(float q[], int nyp, int nx, int nxe, int nypmx) {
/* accumulate extended periodic scalar field in x direction
   linear interpolation, for distributed data
   nyp = number of primary (complete) gridpoints in particle partition
   nx = system length in x direction
   nxe = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells
local data                                                 */
   int k, myp1;
/* accumulate edges of extended field */
   myp1 = nyp + 1;
   for (k = 0; k < myp1; k++) {
      q[nxe*k] += q[nx+nxe*k];
      q[nx+nxe*k] = 0.0;
   }
   return;
}
      
/*--------------------------------------------------------------------*/
void cppois22(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int kstrt, int nyv, int kxp,
              int nyhd) {
/* this subroutine solves 2d poisson's equation in fourier space for
   force/charge (or convolution of electric field over particle shape)
   with periodic boundary conditions, for distributed data.
   for isign = 0, input: isign,ax,ay,affp,nx,ny,jblok,nyv,kxp,nyhd,
   output: ffc
   for isign /= 0, input: q,ffc,isign,nx,ny,nyv,kxp,jblok,nyhd,
   output: fxy,we
   approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
   where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
   the equation used is:
   fx[ky][kx] = -sqrt(-1)*kx*g[ky][kx]*s[ky][kx]*q[ky][kx],
   fy[ky][kx] = -sqrt(-1)*ky*g[ky][kx]*s[ky][kx]*q[ky][kx],
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   g[ky][kx] = (affp/(kx**2+ky**2))*s[ky][kx],
   s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
   fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
   fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
   q[k][j] = complex charge density for fourier mode (jj,k)
   fxy[k][j][0] = x component of complex force/charge,
   fxy[k][j][1] = y component of complex force/charge,
   for fourier mode (jj,k), where jj = j + kxp*(kstrt - 1)
   kxp = number of data values per block
   kstrt = starting data block number
   if isign = 0, form factor array is prepared
   if isign is not equal to 0, force/charge is calculated.
   aimag(ffc[k][j]) = finite-size particle shape factor s
   real(ffc[k][j])) = potential green's function g
   for fourier mode (jj,k), where jj = j + kxp*(kstrt - 1)
   ax/ay = half-width of particle in x/y direction
   affp = normalization constant = nx*ny/np, where np=number of particles
   electric field energy is also calculated, using
   we = nx*ny*sum((affp/(kx**2+ky**2))*|q[ky][kx]*s[ky][kx]|**2)
   nx/ny = system length in x/y direction
   nyv = first dimension of field arrays, must be >= ny
   nyhd = first dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, ks, joff, kxps, j, jj, jk, k, k1;
   float dnx, dny, dkx, dky, at1, at2, at3, at4;
   float complex zero, zt1, zt2;
   double wp;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   ks = kstrt - 1;
   joff = kxp*ks;
   kxps = nxh - joff;
   kxps = 0 > kxps ? 0 : kxps;
   kxps = kxp < kxps ? kxp : kxps;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
   if (isign != 0)
      goto L30;
   if (kstrt > nxh) return;
/* prepare form factor array */
   for (j = 0; j < kxps; j++) {
      dkx = dnx*(float) (j + joff);
      jj = nyhd*j;
      at1 = dkx*dkx;
      at2 = pow((dkx*ax),2);
      for (k = 0; k < nyh; k++) {
         dky = dny*(float) k;
         at3 = dky*dky + at1;
         at4 = exp(-.5*(pow((dky*ay),2) + at2));
         if (at3==0.0) {
            ffc[k+jj] = affp + 1.0*_Complex_I;
         }
         else {
            ffc[k+jj] = (affp*at4/at3) + at4*_Complex_I;
         }
      }
   }
   return;
/* calculate force/charge and sum field energy */
L30: wp = 0.0;
   if (kstrt > nxh)
      goto L70;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (j = 0; j < kxps; j++) {
      dkx = dnx*(float) (j + joff);
      jj = nyhd*j;
      jk = nyv*j;
      if ((j+joff) > 0) {
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            at1 = crealf(ffc[k+jj])*cimagf(ffc[k+jj]);
            at2 = dkx*at1;
            at3 = dny*at1*(float) k;
            zt1 = cimagf(q[k+jk]) - crealf(q[k+jk])*_Complex_I;
            zt2 = cimagf(q[k1+jk]) - crealf(q[k1+jk])*_Complex_I;
            fxy[2*k+2*jk] = at2*zt1;
            fxy[1+2*k+2*jk] = at3*zt1;
            fxy[2*k1+2*jk] = at2*zt2;
            fxy[1+2*k1+2*jk] = -at3*zt2;
            wp += at1*(q[k+jk]*conjf(q[k+jk])
                  + q[k1+jk]*conjf(q[k1+jk]));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         at1 = crealf(ffc[jj])*cimagf(ffc[jj]);
         at3 = dkx*at1;
         zt1 = cimagf(q[jk]) - crealf(q[jk])*_Complex_I;
         fxy[2*jk] = at3*zt1;
         fxy[1+2*jk] = zero;
         fxy[2*k1+2*jk] = zero;
         fxy[1+2*k1+2*jk] = zero;
         wp += at1*(q[jk]*conjf(q[jk]));
      }
   }
/* mode numbers kx = 0, nx/2 */
   if (ks==0) {
      for (k = 1; k < nyh; k++) {
         k1 = ny - k;
         at1 = crealf(ffc[k])*cimagf(ffc[k]);
         at2 = dny*at1*(float) k;
         zt1 = cimagf(q[k]) - crealf(q[k])*_Complex_I;
         fxy[2*k] = zero;
         fxy[1+2*k] = at2*zt1;
         fxy[2*k1] = zero;
         fxy[1+2*k1] = zero;
         wp += at1*(q[k]*conjf(q[k]));
      }
      k1 = 2*nyh;
      fxy[0] = zero;
      fxy[1] = zero;
      fxy[k1] = zero;
      fxy[1+k1] = zero;
   }
L70:
   *we = wp*((float) nx)*((float) ny);
   return;
}

/*--------------------------------------------------------------------*/
void cwpfft2rinit(int mixup[], float complex sct[], int indx, int indy,
                  int nxhyd, int nxyhd) {
/* this subroutine calculates tables needed by a two dimensional
   real to complex fast fourier transform and its inverse.
   input: indx, indy, nxhyd, nxyhd
   output: mixup, sct
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   nxhyd = maximum of (nx/2,ny)
   nxyhd = one half of maximum of (nx,ny)
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, ny, nxy, nxhy, nxyh;
   int  j, k, lb, ll, jb, it;
   float dnxy, arg;
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   ny = 1L<<indy;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
/* bit-reverse index table: mixup[j] = 1 + reversed bits of j */
   for (j = 0; j < nxhy; j++) {
      lb = j;
      ll = 0;
      for (k = 0; k < indx1y; k++) {
         jb = lb/2;
         it = lb - 2*jb;
         lb = jb;
         ll = 2*ll + it;
      }
      mixup[j] = ll + 1;
   }
/* sine/cosine table for the angles 2*n*pi/nxy */
   nxyh = nxy/2;
   dnxy = 6.28318530717959/(float) nxy;
   for (j = 0; j < nxyh; j++) {
      arg = dnxy*(float) j;
      sct[j] = cosf(arg) - sinf(arg)*_Complex_I;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rxx(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int kstrt,
                int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                int nxyhd) {
/* this subroutine performs the x part of a two dimensional real to
   complex fast fourier transform and its inverse, for a subset of y,
   using complex arithmetic, for data which is distributed in blocks
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
   for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
   where N = (nx/2)*ny, and nvp = number of procs
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, an inverse fourier transform is performed
   f[m][n] = (1/nx*ny)*sum(f[k][j]*exp(-sqrt(-1)*2pi*n*j/nx)
   if isign = 1, a forward fourier transform is performed
   f[k][j] = sum(f[m][n]*exp(sqrt(-1)*2pi*n*j/nx)
   kstrt = starting data block number
   kypi = initial y index used
   kypp = number of y indices used
   nxvh = first dimension of f
   kypd = second dimension of f
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxhyd = maximum of (nx/2,ny)
   nxyhd = one half of maximum of (nx,ny)
   the real data is stored in a complex array of length nx/2, ny
   with the odd/even x points stored in the real/imaginary parts.
   in complex notation, fourier coefficients are stored as follows:
   f[k][j] = mode j,kk, where kk = k + kyp*(kstrt - 1)
   0 <= j < nx/2 and 0 <= kk < ny, except for
   f[k][0] = mode nx/2,kk, where ny/2+1 <= kk < ny, and
   imaginary part of f[0][0] = real part of mode nx/2,0 on mode kstrt=0
   imaginary part of f[0][0] = real part of mode nx/2,ny/2
   on mode kstrt=(ny/2)/kyp
   written by viktor k. decyk, ucla
   parallel, RISC optimized version
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny;
   int nxy, nxhy, kypt, j, k, nrx;
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, joff;
   float ani;
   float complex s, t, t1;
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   nxh = nx/2;
   nxhh = nx/4;
   ny = 1L<<indy;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   kypt = kypi + kypp - 1;
   if (kstrt > ny)
      return;
   if (isign > 0)
      goto L100;
/* inverse fourier transform */
   ani = 0.5/(((float) nx)*((float) ny));
   nrx = nxhy/nxh;
/* bit-reverse array elements in x */
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = kypi-1; k < kypt; k++) {
         joff = nxvh*k;
         t = f[j1+joff];
         f[j1+joff] = f[j+joff];
         f[j+joff] = t;
      }
   }
/* first transform in x */
   nrx = nxy/nxh;
   ns = 1;
   for (m = 0; m < indx1; m++) {
      ns2 = ns + ns;
      km = nxhh/ns;
      kmr = km*nrx;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = j + k1;
            j2 = j + k2;
            s = sct[kmr*j];
            for (i = kypi-1; i < kypt; i++) {
               joff = nxvh*i;
               t = s*f[j2+joff];
               f[j2+joff] = f[j1+joff] - t;
               f[j1+joff] += t;
            }
         }
      }
      ns = ns2;
   }
/* unscramble coefficients and normalize */
   kmr = nxy/nx;
   for (j = 1; j < nxhh; j++) {
      t1 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
      for (k = kypi-1; k < kypt; k++) {
         joff = nxvh*k;
         t = conjf(f[nxh-j+joff]);
         s = f[j+joff] + t;
         t = (f[j+joff] - t)*t1;
         f[j+joff] = ani*(s + t);
         f[nxh-j+joff] = ani*conjf(s - t);
      }
   }
   ani = 2.0*ani;
   for (k = kypi-1; k < kypt; k++) {
      joff = nxvh*k;
      f[joff] = ani*((crealf(f[joff]) + cimagf(f[joff]))
                + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I);
      if (nxhh > 0)
         f[nxhh+joff] = ani*conjf(f[nxhh+joff]);
   }
   return;
/* forward fourier transform */
L100: kmr = nxy/nx;
/* scramble coefficients */
   for (j = 1; j < nxhh; j++) {
      t1 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
      for (k = kypi-1; k < kypt; k++) {
         joff = nxvh*k;
         t = conjf(f[nxh-j+joff]);
         s = f[j+joff] + t;
         t = (f[j+joff] - t)*t1;
         f[j+joff] = s + t;
         f[nxh-j+joff] = conjf(s - t);
      }
   }
   for (k = kypi-1; k < kypt; k++) {
      joff = nxvh*k;
      f[joff] = (crealf(f[joff]) + cimagf(f[joff]))
                + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I;
      if (nxhh > 0)
         f[nxhh+joff] = 2.0*conjf(f[nxhh+joff]);
   }
   nrx = nxhy/nxh;
/* bit-reverse array elements in x */
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = kypi-1; k < kypt; k++) {
         joff = nxvh*k;
         t = f[j1+joff];
         f[j1+joff] = f[j+joff];
         f[j+joff] = t;
      }
   }
/* then transform in x */
   nrx = nxy/nxh;
   ns = 1;
   for (m = 0; m < indx1; m++) {
      ns2 = ns + ns;
      km = nxhh/ns;
      kmr = km*nrx;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = j + k1;
            j2 = j + k2;
            s = conjf(sct[kmr*j]);
            for (i = kypi-1; i < kypt; i++) {
               joff = nxvh*i;
               t = s*f[j2+joff];
               f[j2+joff] = f[j1+joff] - t;
               f[j1+joff] += t;
            }
         }
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rxy(float complex g[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int kstrt,
                int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                int nxyhd) {
/* this subroutine performs the y part of a two dimensional real to
   complex fast fourier transform and its inverse, for a subset of x,
   using complex arithmetic, for data which is distributed in blocks
   for isign = (-1,1), input: all, output: g
   for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
   for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
   where N = (nx/2)*ny, and nvp = number of procs
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, an inverse fourier transform is performed
   g[m][n] = sum(g[k][j]*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, a forward fourier transform is performed
   g[k][j] = sum(g[m][n]*exp(sqrt(-1)*2pi*m*k/ny))
   kstrt = starting data block number
   kxp = number of x indices per block
   kxpi = initial x index used
   kxpp = number of x indices used
   nyv = first dimension of g
   kxp = number of data values per block in x
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxhyd = maximum of (nx/2,ny)
   nxyhd = one half of maximum of (nx,ny)
   the real data is stored in a complex array of length nx/2, ny
   with the odd/even x points stored in the real/imaginary parts.
   in complex notation, fourier coefficients are stored as follows:
   g[k][j] = mode jj,k, where jj = j + kxp*(kstrt - 1)
   0 <= jj < nx/2 and 0 <= k < ny, except for
   g[0][k] = mode nx/2,k, where ny/2+1 <= k < ny, and
   imaginary part of g[0][0] = real part of mode nx/2,0 and
   imaginary part of g[1][ny/2] = real part of mode nx/2,ny/2
   on node kstrt=0
   written by viktor k. decyk, ucla
   parallel, RISC optimized version
local data                                                            */
   int indx1, indx1y, nx, nxh, ny, nyh;
   int nxy, nxhy, ks, kxpt, j, k, nry;
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, koff;
   float complex s, t;
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   nxh = nx/2;
   ny = 1L<<indy;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   ks = kstrt - 1;
   kxpt = kxpi + kxpp - 1;
   if (kstrt > nxh)
      return;
   if (isign > 0)
      goto L80;
/* inverse fourier transform */
   nry = nxhy/ny;
/* bit-reverse array elements in y */
   for (k = 0; k < ny; k++) {
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      for (j = kxpi-1; j < kxpt; j++) {
         koff = nyv*j;
         t = g[k1+koff];
         g[k1+koff] = g[k+koff];
         g[k+koff] = t;
      }
   }
/* then transform in y */
   nry = nxy/ny;
   ns = 1;
   for (m = 0; m < indy; m++) {
      ns2 = ns + ns;
      km = nyh/ns;
      kmr = km*nry;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = j + k1;
            j2 = j + k2;
            s = sct[kmr*j];
            for (i = kxpi-1; i < kxpt; i++) {
               koff = nyv*i;
               t = s*g[j2+koff];
               g[j2+koff] = g[j1+koff] - t;
               g[j1+koff] += t;
            }
         }
      }
      ns = ns2;
   }
/* unscramble modes kx = 0, nx/2 */
   if (ks > 0)
      return;
   for (k = 1; k < nyh; k++) {
      if (kxpi==1) {
         s = g[ny-k];
         g[ny-k] = 0.5*(cimagf(g[k] + s) + crealf(g[k] - s)*_Complex_I);
         g[k] = 0.5*(crealf(g[k] + s) + cimagf(g[k] - s)*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
/* scramble modes kx = 0, nx/2 */
L80: if (ks==0) {
      for (k = 1; k < nyh; k++) {
         if (kxpi==1) {
            s = cimagf(g[ny-k]) + crealf(g[ny-k])*_Complex_I;
            g[ny-k] = conjf(g[k] - s);
            g[k] += s;
         }
      }
   }
   nry = nxhy/ny;
/* bit-reverse array elements in y */
   for (k = 0; k < ny; k++) {
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      for (j = kxpi-1; j < kxpt; j++) {
         koff = nyv*j;
         t = g[k1+koff];
         g[k1+koff] = g[k+koff];
         g[k+koff] = t;
      }
   }
/* first transform in y */
   nry = nxy/ny;
   ns = 1;
   for (m = 0; m < indy; m++) {
      ns2 = ns + ns;
      km = nyh/ns;
      kmr = km*nry;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = j + k1;
            j2 = j + k2;
            s = conjf(sct[kmr*j]);
            for (i = kxpi-1; i < kxpt; i++) {
               koff = nyv*i;
               t = s*g[j2+koff];
               g[j2+koff] = g[j1+koff] - t;
               g[j1+koff] += t;
            }
         }
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2r2xx(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                 int nxyhd) {
/* this subroutine performs the x part of 2 two dimensional real to
   complex fast fourier transforms and their inverses, for a subset of y,
   using complex arithmetic, for data which is distributed in blocks
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
   for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
   where N = (nx/2)*ny, and nvp = number of procs
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, an inverse fourier transform is performed
   f[m][n][0:1] = (1/nx*ny)*sum(f[k][j][0:1]*exp(-sqrt(-1)*2pi*n*j/nx)
   if isign = 1, a forward fourier transform is performed
   f[k][j][0:1] = sum(f[m][n][0:1]*exp(sqrt(-1)*2pi*n*j/nx)*
   kstrt = starting data block number
   kypi = initial y index used
   kypp = number of y indices used
   nxvh = first dimension of f
   kypd = second dimension of f
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxhyd = maximum of (nx/2,ny)
   nxyhd = one half of maximum of (nx,ny)
   the real data is stored in a complex array of length nx/2, ny
   with the odd/even x points stored in the real/imaginary parts.
   in complex notation, fourier coefficients are stored as follows:
   f[k][j][0:1] = mode j,kk, where kk = k + kyp*(kstrt - 1)
   0 <= j < nx/2 and 0 <= kk < ny, except for
   f[k][0][0:1] = mode nx/2,kk, where ny/2+1 <= kk < ny, and
   imaginary part of f[0][0][0:1] = real part of mode nx/2,0
   on mode kstrt=0
   imaginary part of f[0][0][0:1] = real part of mode nx/2,ny/2
   on mode kstrt=(ny/2)/kyp
   written by viktor k. decyk, ucla
   parallel, RISC optimized version
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny;
   int nxy, nxhy, kypt, j, k, nrx;
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, joff;
   float ani, at1;
   float complex s, t, t1, t2;
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   nxh = nx/2;
   nxhh = nx/4;
   ny = 1L<<indy;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   kypt = kypi + kypp - 1;
   if (kstrt > ny)
      return;
   if (isign > 0)
      goto L140;
/* inverse fourier transform */
   ani = 0.5/(((float) nx)*((float) ny));
   nrx = nxhy/nxh;
/* swap complex components */
   for (k = kypi-1; k < kypt; k++) {
      for (j = 0; j < nxh; j++) {
         joff = 2*nxvh*k;
         at1 = cimagf(f[2*j+joff]);
         f[2*j+joff] = crealf(f[2*j+joff])
                       + crealf(f[1+2*j+joff])*_Complex_I;
         f[1+2*j+joff] = at1 + cimagf(f[1+2*j+joff])*_Complex_I;
       }
   }
/* bit-reverse array elements in x */
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = kypi-1; k < kypt; k++) {
         joff = 2*nxvh*k;
         t1 = f[2*j1+joff];
         t2 = f[1+2*j1+joff];
         f[2*j1+joff] = f[2*j+joff];
         f[1+2*j1+joff] = f[1+2*j+joff];
         f[2*j+joff] = t1;
         f[1+2*j+joff] = t2;
      }
   }
/* first transform in x */
   nrx = nxy/nxh;
   ns = 1;
   for (m = 0; m < indx1; m++) {
      ns2 = ns + ns;
      km = nxhh/ns;
      kmr = km*nrx;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = j + k1;
            j2 = j + k2;
            s = sct[kmr*j];
            for (i = kypi-1; i < kypt; i++) {
               joff = 2*nxvh*i;
               t1 = s*f[2*j2+joff];
               t2 = s*f[1+2*j2+joff];
               f[2*j2+joff] = f[2*j1+joff] - t1;
               f[1+2*j2+joff] = f[1+2*j1+joff] - t2;
               f[2*j1+joff] += t1;
               f[1+2*j1+joff] += t2;
            }
         }
      }
      ns = ns2;
   }
/* unscramble coefficients and normalize */
   kmr = nxy/nx;
   for (j = 1; j < nxhh; j++) {
      t1 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
      for (k = kypi-1; k < kypt; k++) {
         joff = 2*nxvh*k;
         for (i = 0; i < 2; i++) {
            t = conjf(f[i+2*(nxh-j)+joff]);
            s = f[i+2*j+joff] + t;
            t = (f[i+2*j+joff] - t)*t1;
            f[i+2*j+joff] = ani*(s + t);
            f[i+2*(nxh-j)+joff] = ani*conjf(s - t);
         }
      }
   }
   ani = 2.0*ani;
   for (k = kypi-1; k < kypt; k++) {
      joff = 2*nxvh*k;
      for (i = 0; i < 2; i++) {
         f[i+joff] = ani*((crealf(f[i+joff]) + cimagf(f[i+joff]))
                      + (crealf(f[i+joff]) - cimagf(f[i+joff]))*_Complex_I);
         if (nxhh > 0)
            f[i+2*nxhh+joff] = ani*conjf(f[i+2*nxhh+joff]);
      }
   }
   return;
/* forward fourier transform */
L140: kmr = nxy/nx;
/* scramble coefficients */
   for (j = 1; j < nxhh; j++) {
      t1 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
      for (k = kypi-1; k < kypt; k++) {
         joff = 2*nxvh*k;
         for (i = 0; i < 2; i++) {
            t = conjf(f[i+2*(nxh-j)+joff]);
            s = f[i+2*j+joff] + t;
            t = (f[i+2*j+joff] - t)*t1;
            f[i+2*j+joff] = s + t;
            f[i+2*(nxh-j)+joff] = conjf(s - t);
         }
      }
   }
   for (k = kypi-1; k < kypt; k++) {
      joff = 2*nxvh*k;
      for (i = 0; i < 2; i++) {
         f[i+joff] = (crealf(f[i+joff]) + cimagf(f[i+joff]))
                      + (crealf(f[i+joff]) - cimagf(f[i+joff]))*_Complex_I;
         if (nxhh > 0)
            f[i+2*nxhh+joff] = 2.0*conjf(f[i+2*nxhh+joff]);
      }
   }
   nrx = nxhy/nxh;
/* bit-reverse array elements in x */
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = kypi-1; k < kypt; k++) {
         joff = 2*nxvh*k;
         t1 = f[2*j1+joff];
         t2 = f[1+2*j1+joff];
         f[2*j1+joff] = f[2*j+joff];
         f[1+2*j1+joff] = f[1+2*j+joff];
         f[2*j+joff] = t1;
         f[1+2*j+joff] = t2;
      }
   }
/* then transform in x */
   nrx = nxy/nxh;
   ns = 1;
   for (m = 0; m < indx1; m++) {
      ns2 = ns + ns;
      km = nxhh/ns;
      kmr = km*nrx;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = j + k1;
            j2 = j + k2;
            s = conjf(sct[kmr*j]);
            for (i = kypi-1; i < kypt; i++) {
               joff = 2*nxvh*i;
               t1 = s*f[2*j2+joff];
               t2 = s*f[1+2*j2+joff];
               f[2*j2+joff] = f[2*j1+joff] - t1;
               f[1+2*j2+joff] = f[1+2*j1+joff] - t2;
               f[2*j1+joff] += t1;
               f[1+2*j1+joff] +=  t2;
            }
         }
      }
      ns = ns2;
   }
/* swap complex components */
   for (k = kypi-1; k < kypt; k++) {
      joff = 2*nxvh*k;
      for (j = 0; j < nxh; j++) {
         at1 = cimagf(f[2*j+joff]);
         f[2*j+joff] = crealf(f[2*j+joff])
                       + crealf(f[1+2*j+joff])*_Complex_I;
         f[1+2*j+joff] = at1 + cimagf(f[1+2*j+joff])*_Complex_I;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2r2xy(float complex g[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                 int nxyhd) {
/* this subroutine performs the y part of 2 two dimensional real to
   complex fast fourier transforms and their inverses, for a subset of x,
   using complex arithmetic, for data which is distributed in blocks
   for isign = (-1,1), input: all, output: g
   for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
   for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
   where N = (nx/2)*ny, and nvp = number of procs
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, an inverse fourier transform is performed
   g[n][m][0:1] = sum(g[j][k][0:1]*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, a forward fourier transform is performed
   g[j][k][0:1] = sum(g[n][m][0:1]*exp(sqrt(-1)*2pi*m*k/ny))
   kstrt = starting data block number
   kxpi = initial x index used
   kxpp = number of x indices used
   nyv = first dimension of g
   kxp = number of data values per block in x
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxhyd = maximum of (nx/2,ny)
   nxyhd = one half of maximum of (nx,ny)
   the real data is stored in a complex array of length nx/2, ny
   with the odd/even x points stored in the real/imaginary parts.
   in complex notation, fourier coefficients are stored as follows:
   g[j][k][0:1] = mode jj,k, where jj = j + kxp*(kstrt - 1)
   0 <= jj < nx/2 and 0 <= k < ny, except for
   g[0][k][0:1] = mode nx/2,k, where ny/2+1 <= k < ny, and
   imaginary part of g[0][0][0:1] = real part of mode nx/2,0 and
   imaginary part of g[0][ny/2][0:1] = real part of mode nx/2,ny/2
   on node kstrt=0
   written by viktor k. decyk, ucla
   parallel, RISC optimized version
local data                                                            */
   int indx1, indx1y, nx, nxh, ny, nyh;
   int nxy, nxhy, ks, kxpt, j, k, nry;
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, koff;
   float complex s, t1, t2;
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   nxh = nx/2;
   ny = 1L<<indy;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   ks = kstrt - 1;
   kxpt = kxpi + kxpp - 1;
   if (kstrt > nxh)
      return;
   if (isign > 0)
      goto L90;
/* inverse fourier transform */
   nry = nxhy/ny;
/* bit-reverse array elements in y */
   for (k = 0; k < ny; k++) {
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      for (j = kxpi-1; j < kxpt; j++) {
         koff = 2*nyv*j;
         t1 = g[2*k1+koff];
         t2 = g[1+2*k1+koff];
         g[2*k1+koff] = g[2*k+koff];
         g[1+2*k1+koff] = g[1+2*k+koff];
         g[2*k+koff] = t1;
         g[1+2*k+koff] = t2;
      }
   }
/* then transform in y */
   nry = nxy/ny;
   ns = 1;
   for (m = 0; m < indy; m++) {
      ns2 = ns + ns;
      km = nyh/ns;
      kmr = km*nry;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = j + k1;
            j2 = j + k2;
            s = sct[kmr*j];
            for (i = kxpi-1; i < kxpt; i++) {
               koff = 2*nyv*i;
               t1 = s*g[2*j2+koff];
               t2 = s*g[1+2*j2+koff];
               g[2*j2+koff] = g[2*j1+koff] - t1;
               g[1+2*j2+koff] = g[1+2*j1+koff] - t2;
               g[2*j1+koff] += t1;
               g[1+2*j1+koff] += t2;
            }
         }
      }
      ns = ns2;
   }
/* unscramble modes kx = 0, nx/2 */
   if (ks > 0)
      return;
   for (k = 1; k < nyh; k++) {
      if (kxpi==1) {
         for (i = 0; i < 2; i++) {
            s = g[i+2*(ny-k)];
            g[i+2*(ny-k)] = 0.5*(cimagf(g[i+2*k] + s)
                             + crealf(g[i+2*k] - s)*_Complex_I);
            g[i+2*k] = 0.5*(crealf(g[i+2*k] + s)
                             + cimagf(g[i+2*k] - s)*_Complex_I);
         }
      }
   }
   return;
/* forward fourier transform */
/* scramble modes kx = 0, nx/2 */
L90: if (ks==0) {
      for (k = 1; k < nyh; k++) {
         if (kxpi==1) {
            for (i = 0; i < 2; i++) {
               s = cimagf(g[i+2*(ny-k)])
                   + crealf(g[i+2*(ny-k)])*_Complex_I;
               g[i+2*(ny-k)] = conjf(g[i+2*k] - s);
               g[i+2*k] += s;
            }
         }
      }
   }
   nry = nxhy/ny;
/* bit-reverse array elements in y */
   for (k = 0; k < ny; k++) {
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      for (j = kxpi-1; j < kxpt; j++) {
         koff = 2*nyv*j;
         t1 = g[2*k1+koff];
         t2 = g[1+2*k1+koff];
         g[2*k1+koff] = g[2*k+koff];
         g[1+2*k1+koff] = g[1+2*k+koff];
         g[2*k+koff] = t1;
         g[1+2*k+koff] = t2;
      }
   }
/* first transform in y */
   nry = nxy/ny;
   ns = 1;
   for (m = 0; m < indy; m++) {
      ns2 = ns + ns;
      km = nyh/ns;
      kmr = km*nry;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = j + k1;
            j2 = j + k2;
            s = conjf(sct[kmr*j]);
            for (i = kxpi-1; i < kxpt; i++) {
               koff = 2*nyv*i;
               t1 = s*g[2*j2+koff];
               t2 = s*g[1+2*j2+koff];
               g[2*j2+koff] = g[2*j1+koff] - t1;
               g[1+2*j2+koff] = g[1+2*j1+koff] - t2;
               g[2*j1+koff] += t1;
               g[1+2*j1+koff] += t2;
            }
         }
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2r(float complex f[], float complex g[], float complex bs[],
               float complex br[], int isign, int ntpose, int mixup[],
               float complex sct[], float *ttp, int indx, int indy,
               int kstrt, int nvp, int nxvh, int nyv, int kxp, int kyp,
               int kypd, int nxhyd, int nxyhd) {
/* wrapper function for 2d real to complex fft, with packed data */
/* parallelized with MPI */
/* local data */
   int nxh, ny, ks, kxpp, kypp;
   static int kxpi = 1, kypi = 1;
   float tf;
   double dtime;
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   ny = 1L<<indy;
   ks = kstrt - 1;
   kxpp = nxh - kxp*ks;
   kxpp = 0 > kxpp ? 0 : kxpp;
   kxpp = kxp < kxpp ? kxp : kxpp;
   kypp = ny - kyp*ks;
   kypp = 0 > kypp ? 0 : kypp;
   kypp = kyp < kypp ? kyp : kypp;
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cppfft2rxx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,kypd,
                 nxhyd,nxyhd);
/* transpose f array to g */
      cpwtimera(-1,ttp,&dtime);
      cpptpose(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp,kypd);
      cpwtimera(1,ttp,&dtime);
/* perform y fft */
      cppfft2rxy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
                 nxhyd,nxyhd);
/* transpose g array to f */
      if (ntpose==0) {
         cpwtimera(-1,&tf,&dtime);
         cpptpose(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,kypd,kxp);
         cpwtimera(1,&tf,&dtime);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* transpose f array to g */
      if (ntpose==0) {
         cpwtimera(-1,&tf,&dtime);
         cpptpose(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp,kypd);
         cpwtimera(1,&tf,&dtime);
      }
/* perform y fft */
      cppfft2rxy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
                 nxhyd,nxyhd);
/* transpose g array to f */
      cpwtimera(-1,ttp,&dtime);
      cpptpose(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,kypd,kxp);
      cpwtimera(1,ttp,&dtime);
/* perform x fft */
      cppfft2rxx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,kypd,
                 nxhyd,nxyhd);
   }
   if (ntpose==0)
      *ttp += tf;
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2r2(float complex f[], float complex g[], float complex bs[],
                float complex br[], int isign, int ntpose, int mixup[],
                float complex sct[], float *ttp, int indx, int indy,
                int kstrt, int nvp, int nxvh, int nyv, int kxp, int kyp,
                int kypd, int nxhyd, int nxyhd) {
/* wrapper function for 2 2d real to complex ffts, with packed data */
/* parallelized with MPI */
/* local data */
   int nxh, ny, ks, kxpp, kypp;
   static int kxpi = 1, kypi = 1;
   float tf;
   double dtime;
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   ny = 1L<<indy;
   ks = kstrt - 1;
   kxpp = nxh - kxp*ks;
   kxpp = 0 > kxpp ? 0 : kxpp;
   kxpp = kxp < kxpp ? kxp : kxpp;
   kypp = ny - kyp*ks;
   kypp = 0 > kypp ? 0 : kypp;
   kypp = kyp < kypp ? kyp : kypp;
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cppfft2r2xx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,kypd,
                  nxhyd,nxyhd);
/* transpose f array to g */
      cpwtimera(-1,ttp,&dtime);
      cppntpose(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,2,nxvh,nyv,kxp,kypd);
      cpwtimera(1,ttp,&dtime);
/* perform y fft */
      cppfft2r2xy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
                  nxhyd,nxyhd);
/* transpose g array to f */
      if (ntpose==0) {
         cpwtimera(-1,&tf,&dtime);
         cppntpose(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,2,nyv,nxvh,kypd,
                   kxp);
         cpwtimera(1,&tf,&dtime);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* transpose f array to g */
      if (ntpose==0) {
         cpwtimera(-1,&tf,&dtime);
         cppntpose(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,2,nxvh,nyv,kxp,
                   kypd);
         cpwtimera(1,&tf,&dtime);
      }
/* perform y fft */
      cppfft2r2xy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
                  nxhyd,nxyhd);
/* transpose g array to f */
      cpwtimera(-1,ttp,&dtime);
      cppntpose(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,2,nyv,nxvh,kypd,kxp);
      cpwtimera(1,ttp,&dtime);
/* perform x fft */
      cppfft2r2xx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,kypd,
                  nxhyd,nxyhd);
   }
   if (ntpose==0)
      *ttp += tf;
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void cpdicomp2l_(float *edges, int *nyp, int *noff, int *nypmx,
                 int *nypmn, int *ny, int *kstrt, int *nvp, int *idps) {
   cpdicomp2l(edges,nyp,noff,nypmx,nypmn,*ny,*kstrt,*nvp,*idps);
   return;
}

/*--------------------------------------------------------------------*/
void cpdistr2_(float *part, float *edges, int *npp, int *nps, float *vtx,
               float *vty, float *vdx, float *vdy, int *npx, int *npy,
               int *nx, int *ny, int *idimp, int *npmax, int *idps,
               int *ipbc, int *ierr) {
   cpdistr2(part,edges,npp,*nps,*vtx,*vty,*vdx,*vdy,*npx,*npy,*nx,*ny,
            *idimp,*npmax,*idps,*ipbc,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void cppgpush2l_(float *part, float *fxy, float *edges, int *npp,
                 int *noff, int *ihole, float *qbm, float *dt, float *ek,
                 int *nx, int *ny, int *idimp, int *npmax, int *nxv,
                 int *nypmx, int *idps, int *ntmax, int *ipbc) {
   cppgpush2l(part,fxy,edges,*npp,*noff,ihole,*qbm,*dt,ek,*nx,*ny,*idimp,
              *npmax,*nxv,*nypmx,*idps,*ntmax,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgpost2l_(float *part, float *q, int *npp, int *noff, float *qm,
                 int *idimp, int *npmax, int *nxv, int *nypmx) {
   cppgpost2l(part,q,*npp,*noff,*qm,*idimp,*npmax,*nxv,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppdsortp2yl_(float *parta, float *partb, int *npic, int *npp,
                   int *noff, int *nyp, int *idimp, int *npmax,
                   int *nypm1) {
   cppdsortp2yl(parta,partb,npic,*npp,*noff,*nyp,*idimp,*npmax,*nypm1);
   return;
}

/*--------------------------------------------------------------------*/
void cppcguard2xl_(float *fxy, int *nyp, int *nx, int *ndim, int *nxe,
                   int *nypmx) {
   cppcguard2xl(fxy,*nyp,*nx,*ndim,*nxe,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppaguard2xl_(float *q, int *nyp, int *nx, int *nxe, int *nypmx) {
   cppaguard2xl(q,*nyp,*nx,*nxe,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppois22_(float complex *q, float complex *fxy, int *isign,
               float complex *ffc, float *ax, float *ay, float *affp,
               float *we, int *nx, int *ny, int *kstrt, int *nyv,
               int *kxp, int *nyhd) {
   cppois22(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*kstrt,*nyv,*kxp,
            *nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwpfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                   int *nxhyd, int *nxyhd) {
   cwpfft2rinit(mixup,sct,*indx,*indy,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rxx_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *kstrt,
                 int *kypi, int *kypp, int *nxvh, int *kypd, int *nxhyd,
                 int *nxyhd) {
   cppfft2rxx(f,*isign,mixup,sct,*indx,*indy,*kstrt,*kypi,*kypp,*nxvh,
              *kypd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rxy_(float complex *g, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *kstrt,
                 int *kxpi, int *kxpp, int *nyv, int *kxp, int *nxhyd,
                 int *nxyhd) {
   cppfft2rxy(g,*isign,mixup,sct,*indx,*indy,*kstrt,*kxpi,*kxpp,*nyv,
              *kxp,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2r2xx_(float complex *f, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *kstrt,
                  int *kypi, int *kypp, int *nxvh, int *kypd, int *nxhyd,
                  int *nxyhd) {
   cppfft2r2xx(f,*isign,mixup,sct,*indx,*indy,*kstrt,*kypi,*kypp,*nxvh,
               *kypd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2r2xy_(float complex *g, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *kstrt,
                  int *kxpi, int *kxpp, int *nyv, int *kxp, int *nxhyd,
                  int *nxyhd) {
   cppfft2r2xy(g,*isign,mixup,sct,*indx,*indy,*kstrt,*kxpi,*kxpp,*nyv,
               *kxp,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2r_(float complex *f, float complex *g, float complex *bs,
                float complex *br, int *isign, int *ntpose, int *mixup,
                float complex *sct, float *ttp, int *indx, int *indy,
                int *kstrt, int *nvp, int *nxvh, int *nyv, int *kxp,
                int *kyp, int *kypd, int *nxhyd, int *nxyhd) {
   cwppfft2r(f,g,bs,br,*isign,*ntpose,mixup,sct,ttp,*indx,*indy,*kstrt,
             *nvp,*nxvh,*nyv,*kxp,*kyp,*kypd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2r2_(float complex *f, float complex *g, float complex *bs,
                 float complex *br, int *isign, int *ntpose, int *mixup,
                 float complex *sct, float *ttp, int *indx, int *indy,
                 int *kstrt, int *nvp, int *nxvh, int *nyv, int *kxp,
                 int *kyp, int *kypd, int *nxhyd, int *nxyhd) {
   cwppfft2r2(f,g,bs,br,*isign,*ntpose,mixup,sct,ttp,*indx,*indy,*kstrt,
              *nvp,*nxvh,*nyv,*kxp,*kyp,*kypd,*nxhyd,*nxyhd);
   return;
}
