/* C Library for Skeleton 2D Electrostatic PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "push2.h"

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
void cdistr2(float part[], float vtx, float vty, float vdx, float vdy,
             int npx, int npy, int idimp, int nop, int nx, int ny,
             int ipbc) {
/* for 2d code, this subroutine calculates initial particle co-ordinates
   and velocities with uniform density and maxwellian velocity with drift
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = velocity vx of particle n
   part[n][3] = velocity vy of particle n
   vtx/vty = thermal velocity of electrons in x/y direction
   vdx/vdy = drift velocity of beam electrons in x/y direction
   npx/npy = initial number of particles distributed in x/y direction
   idimp = size of phase space = 4
   nop = number of particles
   nx/ny = system length in x/y direction
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   ranorm = gaussian random number with zero mean and unit variance
local data                                                            */
   int j, k, k1, npxy;
   float edgelx, edgely, at1, at2, at3, sum1, sum2;
   double dsum1, dsum2;
   npxy = npx*npy;
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
/* uniform density profile */
   for (k = 0; k < npy; k++) {
      k1 = idimp*npx*k;
      at3 = edgely + at2*(((float) k) + 0.5);
      for (j = 0; j < npx; j++) {
         part[idimp*j+k1] = edgelx + at1*(((float) j) + 0.5);
         part[1+idimp*j+k1] = at3;
      }
   }
/* maxwellian velocity distribution */
   for (j = 0; j < npxy; j++) {
      part[2+idimp*j] = vtx*ranorm();
      part[3+idimp*j] = vty*ranorm();
   }
/* add correct drift */
   dsum1 = 0.0;
   dsum2 = 0.0;
   for (j = 0; j < npxy; j++) {
      dsum1 += part[2+idimp*j];
      dsum2 += part[3+idimp*j];
   }
   sum1 = dsum1;
   sum2 = dsum2;
   at1 = 1.0/(float) npxy;
   sum1 = at1*sum1 - vdx;
   sum2 = at1*sum2 - vdy;
   for (j = 0; j < npxy; j++) {
      part[2+idimp*j] -= sum1;
      part[3+idimp*j] -= sum2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cdblkp2l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int my, int mx1, int mxy1, int *irc) {
/* this subroutine finds the maximum number of particles in each tile of
   mx, my to calculate size of segmented particle array ppart
   linear interpolation
   part = input particle array
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   kpic = output number of particles per tile
   nppmx = return maximum number of particles in tile
   idimp = size of phase space = 4
   nop = number of particles
   mx/my = number of grids in sorting cell in x and y
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int j, k, n, m, isum, ist, npx, ierr;
   ierr = 0;
/* clear counter array */
   for (k = 0; k < mxy1; k++) {
      kpic[k] = 0;
   }
/* find how many particles in each tile */
   for (j = 0; j < nop; j++) {
      n = part[idimp*j];
      m = part[1+idimp*j];
      n = n/mx;
      m = m/my;
      m = n + mx1*m;
      if (m < mxy1) {
         kpic[m] += 1;
      }
      else {
         ierr = ierr > (m - mxy1 + 1) ? ierr : (m - mxy1 + 1);
      }
   }
/* find maximum */
   isum = 0;
   npx = 0;
   for (k = 0; k < mxy1; k++) {
      ist = kpic[k];
      npx = npx > ist ? npx : ist;
      isum += ist;
   }
   *nppmx = npx;
/* check for errors */
      if (ierr > 0) {
         *irc = ierr;
      }
      else if (isum != nop) {
         *irc = -1;
      }
   return;
}

/*--------------------------------------------------------------------*/
void cppmovin2lt(float part[], float ppart[], int kpic[], int nppmx,
                 int idimp, int nop, int mx, int my, int mx1, int mxy1,
                 int *irc) {
/* this subroutine sorts particles by x,y grid in tiles of mx, my
   and copies to segmented array ppart
   linear interpolation
   input: all except ppart, kpic, output: ppart, kpic
   part/ppart = input/output particle arrays
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   kpic = output number of particles per tile
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 4
   nop = number of particles
   mx/my = number of grids in sorting cell in x and y
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int i, j, k, n, m, ip, ierr;
   ierr = 0;
/* clear counter array */
   for (k = 0; k < mxy1; k++) {
      kpic[k] = 0;
   }
/* find addresses of particles at each tile and reorder particles */
   for (j = 0; j < nop; j++) {
      n = part[idimp*j];
      m = part[1+idimp*j];
      n = n/mx;
      m = m/my;
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
void cppcheck2lt(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                 int ny, int mx, int my, int mx1, int my1, 
                 int *irc) {
/* this subroutine performs a sanity check to make sure particles sorted
   by x,y grid in tiles of mx, my, are all within bounds.
   tiles are assumed to be arranged in 2D linear memory, and transposed
   input: all except irc
   output: irc
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k
   kpic[k] = number of reordered output particles in tile k
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   irc = particle error, returned only if error occurs, when irc > 0
local data                                                            */
   int mxy1, noff, moff, npp, j, k, ist, nn, mm, ierr;
   float edgelx, edgely, edgerx, edgery, dx, dy;
   mxy1 = mx1*my1;
   ierr = 0;
/* sanity check */
   noff = 0;
   moff = 0;
/* loop over tiles */
   for (k = 0; k < mxy1; k++) {
      npp = kpic[k];
      nn = nx - noff;
      nn = mx < nn ? mx : nn;
      mm = ny - moff;
      mm = my < mm ? my : mm;
      edgelx = noff;
      edgerx = noff + nn;
      edgely = moff;
      edgery = moff + mm;
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
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
      noff += mx;
      if (noff >= (mx*mx1)) {
         noff = 0;
         moff += my;
      }
   }
   if (ierr != 0)
      *irc = ierr;
   return;
}

/*--------------------------------------------------------------------*/
void cpois22t(float complex qt[], float complex fxyt[], int isign,
              float complex ffct[], float ax, float ay, float affp,
              float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
              int nyhd) {
/* this subroutine solves 2d poisson's equation in fourier space for
   force/charge (or convolution of electric field over particle shape)
   with periodic boundary conditions, without packed data.
   for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffct
   for isign /= 0, input: qt,ffct,isign,nx,ny,nxvh,nyhd, output: fxyt,we
   approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1
   equation used is:
   fx[ky][kx] = -sqrt(-1)*kx*g[ky][kx]*s[ky][kx]*q[ky][kx],
   fy[ky][kx] = -sqrt(-1)*ky*g[ky][kx]*s[ky][kx]*q[ky][kx],
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   g[ky][kx] = (affp/(kx**2+ky**2))*s(kx,ky),
   s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
   fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
   fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
   qt[j][k] = complex charge density for fourier mode (k,j)
   fxyt[j][0][k] = x component of complex force/charge,
   fxyt[j][1][k] = y component of complex force/charge,
   all for fourier mode (k,j)
   if isign = 0, form factor array is prepared
   if isign is not equal to 0, force/charge is calculated
   aimag(ffct[j][k]) = finite-size particle shape factor s
   for fourier mode (k,j)
   real(ffct([j][k])) = potential green's function g
   for fourier mode (k,j)
   ax/ay = half-width of particle in x/y direction
   affp = normalization constant = nx*ny/np, where np=number of particles
   electric field energy is also calculated, using
   we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
   nx/ny = system length in x/y direction
   nxvh = second dimension of field arrays, must be >= nxh+1
   nyv = first dimension of field arrays, must be >= ny
   nxhd = second dimension of form factor array, must be >= nxh
   nyhd = first dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, nxh1, j, k, k1, jj, jk, jk2;
   float dnx, dny, dkx, dky, at1, at2, at3, at4;
   float complex zero, zt1, zt2;
   double wp;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
   if (isign != 0)
      goto L30;
/* prepare form factor array */
   for (j = 0; j < nxh; j++) {
      dkx = dnx*(float) j;
      jj = nyhd*j;
      at1 = dkx*dkx;
      at2 = pow((dkx*ax),2);
      for (k = 0; k < nyh; k++) {
         dky = dny*(float) k;
         at3 = dky*dky + at1;
         at4 = exp(-0.5*(pow((dky*ay),2) + at2));
         if (at3==0.0) {
            ffct[k+jj] = affp + 1.0*_Complex_I;
         }
         else {
            ffct[k+jj] = (affp*at4/at3) + at4*_Complex_I;
         }
      }
   }
   return;
/* calculate force/charge and sum field energy */
L30: wp = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (j = 1; j < nxh; j++) {
      dkx = dnx*(float) j;
      jj = nyhd*j;
      jk = nyv*j;
      jk2 = 2*jk;
      for (k = 1; k < nyh; k++) {
         k1 = ny - k;
         at1 = crealf(ffct[k+jj])*cimagf(ffct[k+jj]);
         at2 = at1*dkx;
         at3 = at1*dny*(float) k;
         zt1 = cimagf(qt[k+jk]) - crealf(qt[k+jk])*_Complex_I;
         zt2 = cimagf(qt[k1+jk]) - crealf(qt[k1+jk])*_Complex_I;
         fxyt[k+jk2] = at2*zt1;
         fxyt[k+nyv+jk2] = at3*zt1;
         fxyt[k1+jk2] = at2*zt2;
         fxyt[k1+nyv+jk2] = -at3*zt2;
         wp += at1*(qt[k+jk]*conjf(qt[k+jk])
             + qt[k1+jk]*conjf(qt[k1+jk]));
      }
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nyh;
   for (j = 1; j < nxh; j++) {
      jj = nyhd*j;
      jk = nyv*j;
      jk2 = 2*jk;
      at1 = crealf(ffct[jj])*cimagf(ffct[jj]);
      at2 = at1*dnx*(float) j;
      zt1 = cimagf(qt[jk]) - crealf(qt[jk])*_Complex_I;
      fxyt[jk2] = at2*zt1;
      fxyt[nyv+jk2] = zero;
      fxyt[k1+jk2] = zero;
      fxyt[k1+nyv+jk2] = zero;
      wp += at1*(qt[jk]*conjf(qt[jk]));
   }
/* mode numbers kx = 0, nx/2 */
   nxh1 = 2*nyv*nxh;
   for (k = 1; k < nyh; k++) {
      k1 = ny - k;
      at1 = crealf(ffct[k])*cimagf(ffct[k]);
      at3 = at1*dny*(float) k;
      zt1 = cimagf(qt[k]) - crealf(qt[k])*_Complex_I;
      fxyt[k] = zero;
      fxyt[k+nyv] = at3*zt1;
      fxyt[k1] = zero;
      fxyt[k1+nyv] = at3*conjf(zt1);
      fxyt[k+nxh1] = zero;
      fxyt[k+nyv+nxh1] = zero;
      fxyt[k1+nxh1] = zero;
      fxyt[k1+nyv+nxh1] = zero;
      wp += at1*(qt[k]*conjf(qt[k]));
   }
   k1 = nyh;
   fxyt[0] = zero;
   fxyt[nyv] = zero;
   fxyt[k1] = zero;
   fxyt[k1+nyv] = zero;
   fxyt[nxh1] = zero;
   fxyt[nyv+nxh1] = zero;
   fxyt[k1+nxh1] = zero;
   fxyt[k1+nyv+nxh1] = zero;
   *we = wp*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rinit(int mixup[], float complex sct[], int indx, int indy, 
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
   int j, k, lb, ll, jb, it;
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
void cdistr2_(float *part, float *vtx, float *vty, float *vdx, float *vdy,
              int *npx, int *npy, int *idimp, int *nop, int *nx, int *ny,
              int *ipbc) {
   cdistr2(part,*vtx,*vty,*vdx,*vdy,*npx,*npy,*idimp,*nop,*nx,*ny,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdblkp2l_(float *part, int *kpic, int *nppmx, int *idimp, int *nop,
               int *mx, int *my, int *mx1, int *mxy1, int *irc) {
   cdblkp2l(part,kpic,nppmx,*idimp,*nop,*mx,*my,*mx1,*mxy1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppmovin2lt_(float *part, float *ppart, int *kpic, int *nppmx,
                  int *idimp, int *nop, int *mx, int *my, int *mx1,
                  int *mxy1, int *irc) {
   cppmovin2lt(part,ppart,kpic,*nppmx,*idimp,*nop,*mx,*my,*mx1,*mxy1,
               irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcheck2lt_(float *ppart, int *kpic, int *idimp, int *nppmx,
                  int *nx, int *ny, int *mx, int *my, int *mx1,
                  int *my1, int *irc) {
   cppcheck2lt(ppart,kpic,*idimp,*nppmx,*nx,*ny,*mx,*my,*mx1,*my1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpois22t_(float complex *qt, float complex *fxyt, int *isign,
               float complex *ffct, float *ax, float *ay, float *affp,
               float *we, int *nx, int *ny, int *nxvh, int *nyv,
               int *nxhd, int *nyhd) {
   cpois22t(qt,fxyt,*isign,ffct,*ax,*ay,*affp,we,*nx,*ny,*nxvh,*nyv,
            *nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                  int *nxhyd, int *nxyhd) {
   cwfft2rinit(mixup,sct,*indx,*indy,*nxhyd,*nxyhd);
   return;
}

