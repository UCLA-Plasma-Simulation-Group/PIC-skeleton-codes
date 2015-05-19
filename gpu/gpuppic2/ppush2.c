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
void cppdblkp2l(float part[], int kpic[], int npp, int noff, int *nppmx,
                int idimp, int npmax, int mx, int my, int mx1,
                int mxyp1, int *irc) {
/* this subroutine finds the maximum number of particles in each tile of
   mx, my to calculate size of segmented particle array ppart
   linear interpolation, spatial decomposition in y direction
   input: all except kpic, nppmx, output: kpic, nppmx
   part = input particle array
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   kpic = output number of particles per tile
   nppmx = return maximum number of particles in tile
   npp = number of particles in partition
   noff = backmost global gridpoint in particle partition
   idimp = size of phase space = 4
   npmax = maximum number of particles in each partition
   mx/my = number of grids in sorting cell in x and y
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int j, k, n, m, mnoff, isum, ist, npx, ierr;
   mnoff = noff;
   ierr = 0;
/* clear counter array */
   for (k = 0; k < mxyp1; k++) {
      kpic[k] = 0;
   }
/* find how many particles in each tile */
   for (j = 0; j < npp; j++) {
      n = part[idimp*j];
      m = part[1+idimp*j];
      n = n/mx;
      m = (m - mnoff)/my;
      m = n + mx1*m;
      if (m < mxyp1) {
         kpic[m] += 1;
      }
      else {
         ierr = ierr > m-mxyp1+1 ? ierr : m-mxyp1+1;
      }
   }
/* find maximum */
   isum = 0;
   npx = 0;
   for (k = 0; k < mxyp1; k++) {
      ist = kpic[k];
      npx = npx > ist ? npx : ist;
      isum += ist;
   }
   *nppmx = npx;
/* check for errors */
   if (ierr > 0) {
      *irc = ierr;
   }
   else if (isum != npp) {
      *irc = -1;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cpppmovin2lt(float part[], float ppart[], int kpic[], int npp,
                  int noff, int nppmx, int idimp, int npmax, int mx,
                  int my, int mx1, int mxyp1, int *irc) {
/* this subroutine sorts particles by x,y grid in tiles of
   mx, my and copies to segmented array ppart
   linear interpolation, spatial decomposition in y direction
   input: all except ppart, kpic, output: ppart, kpic
   part/ppart = input/output particle arrays
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   kpic = output number of particles per tile
   nppmx = maximum number of particles in tile
   npp = number of particles in partition
   noff = backmost global gridpoint in particle partition
   idimp = size of phase space = 4
   npmax = maximum number of particles in each partition
   mx/my = number of grids in sorting cell in x and y
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int i, j, k, n, m, mnoff, ip, ierr;
   mnoff = noff;
   ierr = 0;
/* clear counter array */
   for (k = 0; k < mxyp1; k++) {
      kpic[k] = 0;
   }
/* find addresses of particles at each tile and reorder particles */
   for (j = 0; j < npp; j++) {
      n = part[idimp*j];
      m = part[1+idimp*j];
      n = n/mx;
      m = (m - mnoff)/my;
      m = n + mx1*m;
      ip = kpic[m];
      if (ip < nppmx) {
         for (i = 0; i < idimp; i++) {
            ppart[ip+nppmx*(i+idimp*m)] = part[i+idimp*j];
         }
      }
      else {
         ierr = ierr > ip-nppmx+1 ? ierr : ip-nppmx+1;
      }
      kpic[m] = ip + 1;
   }
   if (ierr > 0)
      *irc = ierr;
   return;
}

/*--------------------------------------------------------------------*/
void cpppcheck2lt(float ppart[], int kpic[], int noff, int nyp,
                  int idimp, int nppmx, int nx, int mx, int my, int mx1,
                  int myp1, int *irc) {
/* this subroutine performs a sanity check to make sure particles sorted
   by x,y grid in tiles of mx, my, are all within bounds.
   tiles are assumed to be arranged in 2D linear memory, and transposed
   input: all except irc
   output: irc
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k
   kpic[k] = number of reordered output particles in tile k
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx = system length in x direction
   mx/my = number of grids in sorting cell in x/y
   mx1 = (system length in x direction - 1)/mx + 1
   myp1 = (partition length in y direction - 1)/my + 1
   irc = particle error, returned only if error occurs, when irc > 0
local data                                                            */
   int mxyp1, noffp, moffp, nppp, j, k, ist, nn, mm, ierr;
   float edgelx, edgely, edgerx, edgery, dx, dy;
   mxyp1 = mx1*myp1;
   ierr = 0;
/* sanity check */
   noffp = 0;
   moffp = 0;
/* loop over tiles */
   for (k = 0; k < mxyp1; k++) {
      nppp = kpic[k];
      nn = nx - noffp;
      nn = mx < nn ? mx : nn;
      mm = nyp - moffp;
      mm = my < mm ? my : mm;
      edgelx = noffp;
      edgerx = noffp + nn;
      edgely = noff + moffp;
      edgery = noff + moffp + mm;
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
         dx = ppart[j+nppmx*(idimp*k)];
         dy = ppart[j+nppmx*(1+idimp*k)];
/* find particles going out of bounds */
         ist = 0;
         if (dx < edgelx)
            ist = 1;
         if (dx >= edgerx)
            ist = 2;
         if (dy < edgely)
            ist += 3;
         if (dy >= edgery)
            ist += 6;
         if (ist > 0)
            ierr = k + 1;
      }
      noffp += mx;
      if (noffp >= (mx*mx1)) {
         noffp = 0;
         moffp += my;
      }
   }
   if (ierr != 0)
      *irc = ierr;
   return;
}

/*--------------------------------------------------------------------*/
void cppois22t(float complex qt[], float complex fxyt[], int isign,
               float complex ffct[], float ax, float ay, float affp,
               float *we, int nx, int ny, int kstrt, int nyv, int kxp1,
               int nyhd) {
/* this subroutine solves 2d poisson's equation in fourier space for
   force/charge (or convolution of electric field over particle shape)
   with periodic boundary conditions, for distributed data.
   vector length is second dimension.
   for isign = 0, input: isign,ax,ay,affp,nx,ny,nyv,kxp1,nyhd,
   output: ffct
   for isign /= 0, input: qt,ffct,isign,nx,ny,nyv,kxp1,nyhd,
   output: fxyt,we
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
   qt[k][j] = complex charge density for fourier mode (jj,k)
   fxyt[j][0][k] = x component of complex force/charge,
   fxyt[j][1][k] = y component of complex force/charge,
   for fourier mode (jj,k), where jj = j + kxp1*(kstrt - 1)
   kxp1 = number of data values per block for unpacked field data
   kstrt = starting data block number
   if isign = 0, form factor array is prepared
   if isign is not equal to 0, force/charge is calculated.
   aimag(ffct[j][k]) = finite-size particle shape factor s
   real(ffct[j][k])) = potential green's function g
   for fourier mode (jj,k), where jj = j + kxp1*(kstrt - 1)
   ax/ay = half-width of particle in x/y direction
   affp = normalization constant = nx*ny/np, where np=number of particles
   electric field energy is also calculated, using
   we = nx*ny*sum((affp/(kx**2+ky**2))*|q[ky][kx]*s[ky][kx]|**2)
   nx/ny = system length in x/y direction
   nyv = first dimension of field arrays, must be >= ny
   nyhd = first dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, ks, joff, kxps, j, jj, jk, k, j0, j1, k1;
   float dnx, dny, dkx, dky, at1, at2, at3, at4;
   float complex zero, zt1, zt2;
   double wp;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   ks = kstrt - 1;
   joff = kxp1*ks;
   j1 = nxh + 1;
   kxps = j1 - joff;
   kxps = 0 > kxps ? 0 : kxps;
   kxps = kxp1 < kxps ? kxp1 : kxps;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
   if (isign != 0)
      goto L30;
   if (kstrt > j1) return;
/* prepare form factor array */
   for (j = 0; j < kxps; j++) {
      j0 = j + joff;
      dkx = dnx*(float) j0;
      jj = nyhd*j;
      if ((j0 >= 0) && (j0 < nxh)) {
         at1 = dkx*dkx;
         at2 = pow((dkx*ax),2);
         for (k = 0; k < nyh; k++) {
            dky = dny*(float) k;
            at3 = dky*dky + at1;
            at4 = exp(-.5*(pow((dky*ay),2) + at2));
            if (at3==0.0) {
               ffct[k+jj] = affp + 1.0*_Complex_I;
            }
            else {
               ffct[k+jj] = (affp*at4/at3) + at4*_Complex_I;
            }
         }
      }
   }
   return;
/* calculate force/charge and sum field energy */
L30: wp = 0.0;
   if (kstrt > j1)
      goto L80;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (j = 0; j < kxps; j++) {
      j0 = j + joff;
      dkx = dnx*(float) j0;
      jj = nyhd*j;
      jk = nyv*j;
      if ((j0 > 0) && (j0 < nxh)) {
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            at1 = crealf(ffct[k+jj])*cimagf(ffct[k+jj]);
            at2 = dkx*at1;
            at3 = dny*at1*(float) k;
            zt1 = cimagf(qt[k+jk]) - crealf(qt[k+jk])*_Complex_I;
            zt2 = cimagf(qt[k1+jk]) - crealf(qt[k1+jk])*_Complex_I;
            fxyt[k+2*jk] = at2*zt1;
            fxyt[k1+2*jk] = at2*zt2;
            fxyt[k+nyv+2*jk] = at3*zt1;
            fxyt[k1+nyv+2*jk] = -at3*zt2;
            wp += at1*(qt[k+jk]*conjf(qt[k+jk])
                  + qt[k1+jk]*conjf(qt[k1+jk]));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         at1 = crealf(ffct[jj])*cimagf(ffct[jj]);
         at3 = dkx*at1;
         zt1 = cimagf(qt[jk]) - crealf(qt[jk])*_Complex_I;
         fxyt[2*jk] = at3*zt1;
         fxyt[k1+2*jk] = zero;
         fxyt[nyv+2*jk] = zero;
         fxyt[k1+nyv+2*jk] = zero;
         wp += at1*(qt[jk]*conjf(qt[jk]));
      }
   }
/* mode numbers kx = 0 */
   if (ks==0) {
      for (k = 1; k < nyh; k++) {
         k1 = ny - k;
         at1 = crealf(ffct[k])*cimagf(ffct[k]);
         at2 = dny*at1*(float) k;
         zt1 = cimagf(qt[k]) - crealf(qt[k])*_Complex_I;
         fxyt[k] = zero;
         fxyt[k1] = zero;
         fxyt[k+nyv] = at2*zt1;
         fxyt[k1+nyv] = at2*conjf(zt1);
         wp += at1*(qt[k]*conjf(qt[k]));
      }
      k1 = nyh;
      fxyt[0] = zero;
      fxyt[nyv] = zero;
      fxyt[k1] = zero;
      fxyt[k1+nyv] = zero;
   }
/* mode numbers kx = nx/2 */
   if (ks==(nxh/kxp1)) {
      jk = 2*nyv*(kxps-1);
      for (k = 0; k < ny; k++) {
         fxyt[k+jk] = zero;
         fxyt[k+nyv+jk] = zero;
      }
   }
L80:
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
void cppdblkp2l_(float *part, int *kpic, int *npp, int *noff, 
                 int *nppmx, int *idimp, int *npmax, int *mx, int *my,
                 int *mx1,int *mxyp1, int *irc) {
   cppdblkp2l(part,kpic,*npp,*noff,nppmx,*idimp,*npmax,*mx,*my,*mx1,
              *mxyp1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpppmovin2lt_(float *part, float *ppart, int *kpic, int *npp,
                   int *noff, int *nppmx, int *idimp, int *npmax,
                   int *mx, int *my, int *mx1, int *mxyp1, int *irc) {
   cpppmovin2lt(part,ppart,kpic,*npp,*noff,*nppmx,*idimp,*npmax,*mx,*my,
                *mx1,*mxyp1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpppcheck2lt_(float *ppart, int *kpic, int *noff, int *nyp,
                   int *idimp, int *nppmx, int *nx, int *mx, int *my,
                   int *mx1, int *myp1, int *irc) {
   cpppcheck2lt(ppart,kpic,*noff,*nyp,*idimp,*nppmx,*nx,*mx,*my,*mx1,
                *myp1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppois22t_(float complex *qt, float complex *fxyt, int *isign,
                float complex *ffct, float *ax, float *ay, float *affp,
                float *we, int *nx, int *ny, int *kstrt, int *nyv,
                int *kxp1, int *nyhd) {
   cppois22t(qt,fxyt,*isign,ffct,*ax,*ay,*affp,we,*nx,*ny,*kstrt,*nyv,
             *kxp1,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwpfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                   int *nxhyd, int *nxyhd) {
   cwpfft2rinit(mixup,sct,*indx,*indy,*nxhyd,*nxyhd);
   return;
}

