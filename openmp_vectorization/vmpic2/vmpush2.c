/* C Library for Skeleton 2D Electrostatic OpenMP/Vector PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "vmpush2.h"

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
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k
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
void cppmovin2ltp(float part[], float ppart[], int kpic[], int kp[],
                  int nppmx, int idimp, int nop, int mx, int my,
                  int mx1, int mxy1, int *irc) {
/* this subroutine sorts particles by x,y grid in tiles of mx, my
   and copies to segmented array ppart
   designed for NUMA architectures, where memory is associated with the
   processor which first writes a memory location.
   linear interpolation
   input: all except ppart, kpic, output: ppart, kpic
   part/ppart = input/output particle arrays
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k
   kpic = output number of particles per tile
   kp = original location of reordered particle
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 4
   nop = number of particles
   mx/my = number of grids in sorting cell in x and y
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int i, j, k, n, m, ip, npp, ierr;
   ierr = 0;
/* clear counter array */
   for (k = 0; k < mxy1; k++) {
      kpic[k] = 0;
   }
/* find addresses of particles at each tile to reorder particles */
   for (j = 0; j < nop; j++) {
      n = part[idimp*j];
      m = part[1+idimp*j];
      n = n/mx;
      m = m/my;
      m = n + mx1*m;
      ip = kpic[m];
      if (ip < nppmx) {
         kp[ip+nppmx*m] = j;
      }
      else {
         ierr = ierr > ip-nppmx+1 ? ierr : ip-nppmx+1;
      }
      kpic[m] = ip + 1;
   }
/* check for overflow */
   if (ierr > 0) {
      *irc = ierr;
      return;
   }
/* copy reordered particles */
#pragma omp parallel for private(i,j,k,m,npp)
   for (k = 0; k < mxy1; k++) {
      npp = kpic[k];
      for (j = 0; j < npp; j++) {
         m = kp[j+nppmx*k];
         for (i = 0; i < idimp; i++) {
            ppart[j+nppmx*(i+idimp*k)] = part[i+idimp*m];
         }
      }
   }
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
   int mxy1, noff, moff, npp, j, k, ist, nn, mm;
   float edgelx, edgely, edgerx, edgery, dx, dy;
   mxy1 = mx1*my1;
/* loop over tiles */
#pragma omp parallel for \
private(j,k,noff,moff,npp,nn,mm,ist,edgelx,edgely,edgerx,edgery,dx,dy)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
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
            *irc = k + 1;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cgppush2lt(float ppart[], float fxy[], int kpic[], float qbm,
                float dt, float *ek, int idimp, int nppmx, int nx,
                int ny, int mx, int my, int nxv, int nyv, int mx1,
                int mxy1, int ipbc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions.
   OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   44 flops/particle, 12 loads, 4 stores
   input: all, output: ppart, ek
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = velocity vx of particle n in tile m
   ppart[m][3][n] = velocity vy of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp;
   int i, j, k, nn, mm, mxv;
   float qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
   float sfxy[2*MXV*MYV];
/* float sfxy[2*(mx+1)*(my+1)]; */
   double sum1, sum2;
   mxv = mx + 1;
   qtm = qbm*dt;
   sum2 = 0.0;
/* set boundary values */
   edgelx = 0.0f;
   edgely = 0.0f;
   edgerx = (float) nx;
   edgery = (float) ny;
   if (ipbc==2) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0f;
      edgerx = (float) (nx-1);
   }
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,x,y,dxp,dyp,amx,amy,dx,dy,vx, \
vy,sum1,sfxy) \
reduction(+:sum2)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[2*(i+mxv*j)] = fxy[2*(i+noff+nxv*(j+moff))];
            sfxy[1+2*(i+mxv*j)] = fxy[1+2*(i+noff+nxv*(j+moff))];
         }
      }
      sum1 = 0.0;
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nn = 2*(nn - noff + mxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find acceleration */
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dx = amy*(dxp*sfxy[nn+2] + dx);
         dy = amy*(dxp*sfxy[nn+3] + dy);
         nn += 2*mxv;
         vx = amx*sfxy[nn];
         vy = amx*sfxy[nn+1];
         dx += dyp*(dxp*sfxy[nn+2] + vx);
         dy += dyp*(dxp*sfxy[nn+3] + vy);
/* new velocity */
         dxp = ppart[j+2*nppmx+npoff];
         dyp = ppart[j+3*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         sum1 += dxp*dxp + dyp*dyp;
/* new position */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               vx = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               vy = -vy;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               vx = -vx;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
/* set new velocity */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
      }
      sum2 += sum1;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum2;
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgppushf2lt(float ppart[], float fxy[], int kpic[], int ncl[],
                 int ihole[], float qbm, float dt, float *ek, int idimp,
                 int nppmx, int nx, int ny, int mx, int my, int nxv,
                 int nyv, int mx1, int mxy1, int ntmax, int *irc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   44 flops/particle, 12 loads, 4 stores
   input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = velocity vx of particle n in tile m
   ppart[m][3][n] = velocity vy of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp;
   int i, j, k, ih, nh, nn, mm, mxv;
   float qtm, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float sfxy[2*MXV*MYV];
/* float sfxy[2*(mx+1)*(my+1)]; */
   double sum1, sum2;
   mxv = mx + 1;
   qtm = qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   sum2 = 0.0;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,ih,nh,x,y,dxp,dyp,amx,amy, \
dx,dy,vx,vy,edgelx,edgely,edgerx,edgery,sum1,sfxy) \
reduction(+:sum2)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
      nn = nx - noff;
      nn = mx < nn ? mx : nn;
      mm = ny - moff;
      mm = my < mm ? my : mm;
      edgelx = noff;
      edgerx = noff + nn;
      edgely = moff;
      edgery = moff + mm;
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
/* load local fields from global array */
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[2*(i+mxv*j)] = fxy[2*(i+noff+nxv*(j+moff))];
            sfxy[1+2*(i+mxv*j)] = fxy[1+2*(i+noff+nxv*(j+moff))];
         }
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
      sum1 = 0.0;
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nn = 2*(nn - noff + mxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find acceleration */
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dx = amy*(dxp*sfxy[nn+2] + dx);
         dy = amy*(dxp*sfxy[nn+3] + dy);
         nn += 2*mxv;
         vx = amx*sfxy[nn];
         vy = amx*sfxy[nn+1];
         dx += dyp*(dxp*sfxy[nn+2] + vx);
         dy += dyp*(dxp*sfxy[nn+3] + vy);
/* new velocity */
         dxp = ppart[j+2*nppmx+npoff];
         dyp = ppart[j+3*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         sum1 += dxp*dxp + dyp*dyp;
/* new position */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* find particles going out of bounds */
         mm = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                              */
         if (dx >= edgerx) {
            if (dx >= anx)
               dx -= anx;
            mm = 2;
         }
         else if (dx < edgelx) {
            if (dx < 0.0f) {
               dx += anx;
               if (dx < anx)
                  mm = 1;
               else
                  dx = 0.0;
            }
            else {
               mm = 1;
            }
         }
         if (dy >= edgery) {
            if (dy >= any)
               dy -= any;
            mm += 6;
         }
         else if (dy < edgely) {
            if (dy < 0.0) {
               dy += any;
               if (dy < any)
                  mm += 3;
               else
                  dy = 0.0;
            }
            else {
               mm += 3;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
/* set new velocity */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
/* increment counters */
         if (mm > 0) {
            ncl[mm+8*k-1] += 1;
            ih += 1;
            if (ih <= ntmax) {
               ihole[2*(ih+(ntmax+1)*k)] = j + 1;
               ihole[1+2*(ih+(ntmax+1)*k)] = mm;
            }
            else {
               nh = 1;
            }
         }
      }
      sum2 += sum1;
/* set error and end of file flag */
/* ihole overflow */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*k] = ih;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum2;
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cvgppush2lt(float ppart[], float fxy[], int kpic[], float qbm,
                 float dt, float *ek, int idimp, int nppmx, int nx,
                 int ny, int mx, int my, int nxv, int nyv, int mx1,
                 int mxy1, int ipbc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions.
   vectorizable/OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   44 flops/particle, 12 loads, 4 stores
   input: all, output: ppart, ek
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = velocity vx of particle n in tile m
   ppart[m][3][n] = velocity vy of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
#define NPBLK             32
#define LVECT             4
   int noff, moff, npoff, npp, ipp, joff, nps;
   int i, j, k, m, nn, mm, lxv;
   float qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
   float sfxy[2*MXV*MYV];
/* float sfxy[2*(mx+1)*(my+1)]; */
/* scratch arrays */
   int n[NPBLK];
   float s[NPBLK*LVECT], t[NPBLK*2];
   double sum1, sum2;
   lxv = mx + 1;
   qtm = qbm*dt;
   sum2 = 0.0;
/* set boundary values */
   edgelx = 0.0f;
   edgely = 0.0f;
   edgerx = (float) nx;
   edgery = (float) ny;
   if (ipbc==2) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0f;
      edgerx = (float) (nx-1);
   }
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,noff,moff,npp,npoff,ipp,joff,nps,nn,mm,x,y,dxp,dyp, \
amx,amy,dx,dy,vx,vy,sum1,sfxy,n,s,t) \
reduction(+:sum2)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[2*(i+lxv*j)] = fxy[2*(i+noff+nxv*(j+moff))];
            sfxy[1+2*(i+lxv*j)] = fxy[1+2*(i+noff+nxv*(j+moff))];
         }
      }
      sum1 = 0.0;
/* loop over particles in tile */
      ipp = npp/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m;
/* inner loop over particles in block */
         for (j = 0; j < NPBLK; j++) {
/* find interpolation weights */
            x = ppart[j+joff+npoff];
            y = ppart[j+joff+nppmx+npoff];
            nn = x;
            mm = y;
            dxp = x - (float) nn;
            dyp = y - (float) mm;
            n[j] = nn - noff + lxv*(mm - moff);
            amx = 1.0f - dxp;
            amy = 1.0f - dyp;
            s[j] = amx*amy;
            s[j+NPBLK] = dxp*amy;
            s[j+2*NPBLK] = amx*dyp;
            s[j+3*NPBLK] = dxp*dyp;
            t[j] = x;
            t[j+NPBLK] = y;
         }
/* find acceleration */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + lxv - 2;
            dx = 0.0f;
            dy = 0.0f;
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               if (i > 1)
                  nn = mm;
               dx += sfxy[2*(i+nn)]*s[j+NPBLK*i];
               dy += sfxy[1+2*(i+nn)]*s[j+NPBLK*i];
            }
            s[j] = dx;
            s[j+NPBLK] = dy;
         }
/* new velocity */
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
            dxp = ppart[j+joff+2*nppmx+npoff];
            dyp = ppart[j+joff+3*nppmx+npoff];
            vx = dxp + qtm*s[j];
            vy = dyp + qtm*s[j+NPBLK];
/* average kinetic energy */
            dxp += vx;
            dyp += vy;
            sum1 += dxp*dxp + dyp*dyp;
/* new position */
            s[j] = x + vx*dt;
            s[j+NPBLK] = y + vy*dt;
            s[j+2*NPBLK] = vx;
            s[j+3*NPBLK] = vy;
         }
/* check boundary conditions */
#pragma novector
         for (j = 0; j < NPBLK; j++) {
            dx = s[j];
            dy = s[j+NPBLK];
            vx = s[j+2*NPBLK];
            vy = s[j+3*NPBLK];
/* reflecting boundary conditions */
            if (ipbc==2) {
               if ((dx < edgelx) || (dx >= edgerx)) {
                  dx = t[j];
                  vx = -vx;
               }
               if ((dy < edgely) || (dy >= edgery)) {
                  dy = t[j+NPBLK];
                  vy = -vy;
               }
            }
/* mixed reflecting/periodic boundary conditions */
            else if (ipbc==3) {
               if ((dx < edgelx) || (dx >= edgerx)) {
                  dx = t[j];
                  vx = -vx;
               }
            }
/* set new position */
            ppart[j+joff+npoff] = dx;
            ppart[j+joff+nppmx+npoff] = dy;
/* set new velocity */
            ppart[j+joff+2*nppmx+npoff] = vx;
            ppart[j+joff+3*nppmx+npoff] = vy;
         }
      }
      nps = NPBLK*ipp;
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nn = 2*(nn - noff + lxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find acceleration */
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dx = amy*(dxp*sfxy[nn+2] + dx);
         dy = amy*(dxp*sfxy[nn+3] + dy);
         nn += 2*lxv;
         vx = amx*sfxy[nn];
         vy = amx*sfxy[nn+1];
         dx += dyp*(dxp*sfxy[nn+2] + vx);
         dy += dyp*(dxp*sfxy[nn+3] + vy);
/* new velocity */
         dxp = ppart[j+2*nppmx+npoff];
         dyp = ppart[j+3*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         sum1 += dxp*dxp + dyp*dyp;
/* new position */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               vx = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               vy = -vy;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               vx = -vx;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
/* set new velocity */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
      }
      sum2 += sum1;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum2;
   return;
#undef LVECT
#undef NPBLK
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cvgppushf2lt(float ppart[], float fxy[], int kpic[], int ncl[],
                  int ihole[], float qbm, float dt, float *ek,
                  int idimp, int nppmx, int nx, int ny, int mx, int my,
                  int nxv, int nyv, int mx1, int mxy1, int ntmax,
                  int *irc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   vectorizable/OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   44 flops/particle, 12 loads, 4 stores
   input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = velocity vx of particle n in tile m
   ppart[m][3][n] = velocity vy of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             33
#define MYV             33
#define NPBLK             32
#define LVECT             4
   int noff, moff, npoff, npp, ipp, joff, nps;
   int i, j, k, m, ih, nh, nn, mm, lxv;
   float qtm, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float sfxy[2*MXV*MYV];
/* float sfxy[2*(mx+1)*(my+1)]; */
/* scratch arrays */
   int n[NPBLK];
   float s[NPBLK*LVECT], t[NPBLK*2];
   double sum1, sum2;
   lxv = mx + 1;
   qtm = qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   sum2 = 0.0;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,noff,moff,npp,npoff,ipp,joff,nps,nn,mm,ih,nh,x,y,dxp, \
dyp,amx,amy,dx,dy,vx,vy,edgelx,edgely,edgerx,edgery,sum1,sfxy,n,s,t) \
reduction(+:sum2)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
      nn = nx - noff;
      nn = mx < nn ? mx : nn;
      mm = ny - moff;
      mm = my < mm ? my : mm;
      edgelx = noff;
      edgerx = noff + nn;
      edgely = moff;
      edgery = moff + mm;
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
/* load local fields from global array */
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[2*(i+lxv*j)] = fxy[2*(i+noff+nxv*(j+moff))];
            sfxy[1+2*(i+lxv*j)] = fxy[1+2*(i+noff+nxv*(j+moff))];
         }
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
      sum1 = 0.0;
/* loop over particles in tile */
      ipp = npp/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m;
/* inner loop over particles in block */
         for (j = 0; j < NPBLK; j++) {
/* find interpolation weights */
            x = ppart[j+joff+npoff];
            y = ppart[j+joff+nppmx+npoff];
            nn = x;
            mm = y;
            dxp = x - (float) nn;
            dyp = y - (float) mm;
            n[j] = nn - noff + lxv*(mm - moff);
            amx = 1.0f - dxp;
            amy = 1.0f - dyp;
            s[j] = amx*amy;
            s[j+NPBLK] = dxp*amy;
            s[j+2*NPBLK] = amx*dyp;
            s[j+3*NPBLK] = dxp*dyp;
            t[j] = x;
            t[j+NPBLK] = y;
         }
/* find acceleration */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + lxv - 2;
            dx = 0.0f;
            dy = 0.0f;
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               if (i > 1)
                  nn = mm;
               dx += sfxy[2*(i+nn)]*s[j+NPBLK*i];
               dy += sfxy[1+2*(i+nn)]*s[j+NPBLK*i];
            }
            s[j] = dx;
            s[j+NPBLK] = dy;
         }
/* new velocity */
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
            dxp = ppart[j+joff+2*nppmx+npoff];
            dyp = ppart[j+joff+3*nppmx+npoff];
            vx = dxp + qtm*s[j];
            vy = dyp + qtm*s[j+NPBLK];
/* average kinetic energy */
            dxp += vx;
            dyp += vy;
            sum1 += dxp*dxp + dyp*dyp;
/* new position */
            s[j] = x + vx*dt;
            s[j+NPBLK] = y + vy*dt;
            s[j+2*NPBLK] = vx;
            s[j+3*NPBLK] = vy;
         }
/* check boundary conditions */
#pragma novector
         for (j = 0; j < NPBLK; j++) {
            dx = s[j];
            dy = s[j+NPBLK];
/* find particles going out of bounds */
            mm = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                              */
            if (dx >= edgerx) {
               if (dx >= anx)
                  dx -= anx;
               mm = 2;
            }
            else if (dx < edgelx) {
               if (dx < 0.0f) {
                  dx += anx;
                  if (dx < anx)
                     mm = 1;
                  else
                     dx = 0.0;
               }
               else {
                  mm = 1;
               }
            }
            if (dy >= edgery) {
               if (dy >= any)
                  dy -= any;
               mm += 6;
            }
            else if (dy < edgely) {
               if (dy < 0.0) {
                  dy += any;
                  if (dy < any)
                     mm += 3;
                  else
                     dy = 0.0;
               }
               else {
                  mm += 3;
               }
            }
/* set new position */
            ppart[j+joff+npoff] = dx;
            ppart[j+joff+nppmx+npoff] = dy;
/* set new velocity */
            ppart[j+joff+2*nppmx+npoff] = s[j+2*NPBLK];
            ppart[j+joff+3*nppmx+npoff] = s[j+3*NPBLK];
/* increment counters */
            if (mm > 0) {
               ncl[mm+8*k-1] += 1;
               ih += 1;
               if (ih <= ntmax) {
                  ihole[2*(ih+(ntmax+1)*k)] = j + joff + 1;
                  ihole[1+2*(ih+(ntmax+1)*k)] = mm;
               }
               else {
                  nh = 1;
               }
            }
         }
      }
      nps = NPBLK*ipp;
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nn = 2*(nn - noff + lxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find acceleration */
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dx = amy*(dxp*sfxy[nn+2] + dx);
         dy = amy*(dxp*sfxy[nn+3] + dy);
         nn += 2*lxv;
         vx = amx*sfxy[nn];
         vy = amx*sfxy[nn+1];
         dx += dyp*(dxp*sfxy[nn+2] + vx);
         dy += dyp*(dxp*sfxy[nn+3] + vy);
/* new velocity */
         dxp = ppart[j+2*nppmx+npoff];
         dyp = ppart[j+3*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         sum1 += dxp*dxp + dyp*dyp;
/* new position */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* find particles going out of bounds */
         mm = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                              */
         if (dx >= edgerx) {
            if (dx >= anx)
               dx -= anx;
            mm = 2;
         }
         else if (dx < edgelx) {
            if (dx < 0.0f) {
               dx += anx;
               if (dx < anx)
                  mm = 1;
               else
                  dx = 0.0;
            }
            else {
               mm = 1;
            }
         }
         if (dy >= edgery) {
            if (dy >= any)
               dy -= any;
            mm += 6;
         }
         else if (dy < edgely) {
            if (dy < 0.0) {
               dy += any;
               if (dy < any)
                  mm += 3;
               else
                  dy = 0.0;
            }
            else {
               mm += 3;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
/* set new velocity */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
/* increment counters */
         if (mm > 0) {
            ncl[mm+8*k-1] += 1;
            ih += 1;
            if (ih <= ntmax) {
               ihole[2*(ih+(ntmax+1)*k)] = j + 1;
               ihole[1+2*(ih+(ntmax+1)*k)] = mm;
            }
            else {
               nh = 1;
            }
         }
      }
      sum2 += sum1;
/* set error and end of file flag */
/* ihole overflow */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*k] = ih;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum2;
   return;
#undef LVECT
#undef NPBLK
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgppost2lt(float ppart[], float q[], int kpic[], float qm,
                int nppmx, int idimp, int mx, int my, int nxv, int nyv,
                int mx1, int mxy1) {
/* for 2d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   17 flops/particle, 6 loads, 4 stores
   input: all, output: q
   charge density is approximated by values at the nearest grid points
   q(n,m)=qm*(1.-dx)*(1.-dy)
   q(n+1,m)=qm*dx*(1.-dy)
   q(n,m+1)=qm*(1.-dx)*dy
   q(n+1,m+1)=qm*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   q[k][j] = charge density at grid point j,k
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 4
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of charge array, must be >= nx+1
   nyv = second dimension of charge array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nn, mm;
   float x, y, dxp, dyp, amx, amy;
   float sq[MXV*MYV];
/* float sq[(mx+1)*(my+1)]; */
   mxv = mx + 1;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,x,y,dxp,dyp,amx,amy,sq)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < mxv*(my+1); j++) {
         sq[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         nn = nn - noff + mxv*(mm - moff);
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit charge within tile to local accumulator */
         x = sq[nn] + amx*amy;
         y = sq[nn+1] + dxp*amy;
         sq[nn] = x;
         sq[nn+1] = y;
         nn += mxv;
         x = sq[nn] + amx*dyp;
         y = sq[nn+1] + dxp*dyp;
         sq[nn] = x;
         sq[nn+1] = y;
      }
/* deposit charge to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            q[i+noff+nxv*(j+moff)] += sq[i+mxv*j];
         }
      }
/* deposit charge to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         q[i+noff+nxv*moff] += sq[i];
         if (mm > my) {
#pragma omp atomic
            q[i+noff+nxv*(mm+moff-1)] += sq[i+mxv*(mm-1)];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         q[noff+nxv*(j+moff)] += sq[mxv*j];
         if (nn > mx) {
#pragma omp atomic
            q[nn+noff-1+nxv*(j+moff)] += sq[nn-1+mxv*j];
         }
      }
   }
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cvgppost2lt(float ppart[], float q[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int my, int nxv, int nyv,
                 int mx1, int mxy1) {
/* for 2d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   vectorizable/OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   17 flops/particle, 6 loads, 4 stores
   input: all, output: q
   charge density is approximated by values at the nearest grid points
   q(n,m)=qm*(1.-dx)*(1.-dy)
   q(n+1,m)=qm*dx*(1.-dy)
   q(n,m+1)=qm*(1.-dx)*dy
   q(n+1,m+1)=qm*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   q[k][j] = charge density at grid point j,k
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 4
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of charge array, must be >= nx+1
   nyv = second dimension of charge array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
local data                                                            */
#define MXV             33
#define MYV             33
#define NPBLK             32
#define LVECT             4
   int noff, moff, npoff, npp, ipp, joff, nps;
   int i, j, k, m, nn, mm, lxv;
   float x, y, dxp, dyp, amx, amy;
   float sq[MXV*MYV];
/* float sq[(mx+1)*(my+1)]; */
/* scratch arrays */
   int n[NPBLK];
   float s[NPBLK*LVECT];
   lxv = mx + 1;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,noff,moff,npp,npoff,ipp,joff,nps,nn,mm,x,y,dxp,dyp, \
amx,amy,sq,n,s)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < lxv*(my+1); j++) {
         sq[j] = 0.0f;
      }
/* loop over particles in tile */
      ipp = npp/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m;
/* inner loop over particles in block */
         for (j = 0; j < NPBLK; j++) {
/* find interpolation weights */
            x = ppart[j+joff+npoff];
            y = ppart[j+joff+nppmx+npoff];
            nn = x;
            mm = y;
            dxp = qm*(x - (float) nn);
            dyp = y - (float) mm;
            n[j] = nn - noff + lxv*(mm - moff);
            amx = qm - dxp;
            amy = 1.0f - dyp;
            s[j] = amx*amy;
            s[j+NPBLK] = dxp*amy;
            s[j+2*NPBLK] = amx*dyp;
            s[j+3*NPBLK] = dxp*dyp;
        }
/* deposit charge within tile to local accumulator */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + lxv - 2;
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               if (i > 1)
                  nn = mm;
               sq[i+nn] += s[j+NPBLK*i];
            }
         }
      }
      nps = NPBLK*ipp;
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         nn = nn - noff + lxv*(mm - moff);
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit charge within tile to local accumulator */
         x = sq[nn] + amx*amy;
         y = sq[nn+1] + dxp*amy;
         sq[nn] = x;
         sq[nn+1] = y;
         nn += lxv;
         x = sq[nn] + amx*dyp;
         y = sq[nn+1] + dxp*dyp;
         sq[nn] = x;
         sq[nn+1] = y;
      }
/* deposit charge to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
#pragma ivdep
         for (i = 1; i < nn; i++) {
            q[i+noff+nxv*(j+moff)] += sq[i+lxv*j];
         }
      }
/* deposit charge to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         q[i+noff+nxv*moff] += sq[i];
         if (mm > my) {
#pragma omp atomic
            q[i+noff+nxv*(mm+moff-1)] += sq[i+lxv*(mm-1)];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         q[noff+nxv*(j+moff)] += sq[lxv*j];
         if (nn > mx) {
#pragma omp atomic
            q[nn+noff-1+nxv*(j+moff)] += sq[nn-1+lxv*j];
         }
      }
   }
   return;
#undef LVECT
#undef NPBLK
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cviscan2(int *isdata, int *mb, int nths) {
/* performs vectorizable prefix reduction of integer data */
/* using binary tree method. */
/* local data */
   int j, kxs, lb, ns;
   ns = nths/2;
   for (j = 0; j < ns; j++) {
      mb[j] = j;
   }
   kxs = 1;
   while (kxs < nths) {
#pragma ivdep
      for (j = 0; j < ns; j++) {
         lb = kxs*mb[j];
         if ((j+lb+kxs) < nths) {
            isdata[j+lb+kxs] += isdata[2*lb+kxs-1];
         }
         mb[j] >>= 1;
      }
      kxs <<= 1;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cpporder2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int nx, int ny,
                 int mx, int my, int mx1, int my1, int npbmx, int ntmax,
                 int *irc) {
/* this subroutine sorts particles by x,y grid in tiles of mx, my
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 2D linear memory
   algorithm has 3 steps.  first, one finds particles leaving tile and
   stores their number in each directon, location, and destination in ncl
   and ihole.  second, a prefix scan of ncl is performed and departing
   particles are buffered in ppbuff in direction order.  finally, we copy
   the incoming particles from other tiles into ppart.
   input: all except ppbuff, ncl, ihole, irc
   output: ppart, ppbuff, kpic, ncl, ihole, irc
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k 
   ppbuff[k][i][n] = i co-ordinate of particle n in tile k
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = direction destination of particle leaving hole
   all for tile k
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int mxy1, noff, moff, npoff, npp, nboff, ncoff;
   int i, j, k, ii, kx, ky, ih, nh, ist, nn, mm, isum;
   int ip, j1, j2, kxl, kxr, kk, kl, kr;
   float anx, any, edgelx, edgely, edgerx, edgery, dx, dy;
   int ks[8];
   mxy1 = mx1*my1;
   anx = (float) nx;
   any = (float) ny;
/* find and count particles leaving tiles and determine destination */
/* update ppart, ihole, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(j,k,noff,moff,npp,npoff,nn,mm,ih,nh,ist,dx,dy,edgelx,edgely, \
edgerx,edgery)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
      nn = nx - noff;
      nn = mx < nn ? mx : nn;
      mm = ny - moff;
      mm = my < mm ? my : mm;
      ih = 0;
      nh = 0;
      edgelx = noff;
      edgerx = noff + nn;
      edgely = moff;
      edgery = moff + mm;
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
         dx = ppart[j+npoff];
         dy = ppart[j+nppmx+npoff];
/* find particles going out of bounds */
         ist = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* ist = direction particle is going                             */
         if (dx >= edgerx) {
            if (dx >= anx)
               ppart[j+npoff] = dx - anx;
            ist = 2;
         }
         else if (dx < edgelx) {
            if (dx < 0.0) {
               dx += anx;
               if (dx < anx)
                  ist = 1;
               else
                  dx = 0.0;
               ppart[j+npoff] = dx;
            }
            else {
               ist = 1;
            }
         }
         if (dy >= edgery) {
            if (dy >= any)
               ppart[j+nppmx+npoff] = dy - any;
            ist += 6;
         }
         else if (dy < edgely) {
            if (dy < 0.0) {
               dy += any;
               if (dy < any)
                  ist += 3;
               else
                  dy = 0.0;
               ppart[j+nppmx+npoff] = dy;
            }
            else {
               ist += 3;
            }
         }
         if (ist > 0) {
            ncl[ist+8*k-1] += 1;
            ih += 1;
            if (ih <= ntmax) {
               ihole[2*(ih+(ntmax+1)*k)] = j + 1;
               ihole[1+2*(ih+(ntmax+1)*k)] = ist;
            }
            else {
               nh = 1;
            }
         }
      }
/* set error and end of file flag */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*k] = ih;
   }
/* ihole overflow */
   if (*irc > 0)
      return;

/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,npoff,nboff,isum,ist,nh,ip,j1,ii)
   for (k = 0; k < mxy1; k++) {
      npoff = idimp*nppmx*k;
      nboff = idimp*npbmx*k;
/* find address offset for ordered ppbuff array */
      isum = 0;
      for (j = 0; j < 8; j++) {
         ist = ncl[j+8*k];
         ncl[j+8*k] = isum;
         isum += ist;
      }
      nh = ihole[2*(ntmax+1)*k];
      ip = 0;
/* loop over particles leaving tile */
      for (j = 0; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+(ntmax+1)*k)] - 1;
         ist = ihole[1+2*(j+1+(ntmax+1)*k)];
         ii = ncl[ist+8*k-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*i+nboff]
               = ppart[j1+nppmx*i+npoff];
            }
         }
         else {
            ip = 1;
         }
         ncl[ist+8*k-1] = ii + 1;
      }
/* set error */
      if (ip > 0)
         *irc = ncl[7+8*k];
   }
/* ppbuff overflow */
   if (*irc > 0)
      return;

/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,ii,kk,npp,npoff,nboff,kx,ky,kl,kr,kxl,kxr,ih,nh,nn, \
ncoff,ist,j1,j2,ip,ks)
   for (k = 0; k < mxy1; k++) {
      npp = kpic[k];
      npoff = idimp*nppmx*k;
      ky = k/mx1;
/* loop over tiles in y, assume periodic boundary conditions */
      kk = ky*mx1;
/* find tile above */
      kl = ky - 1;
      if (kl < 0)
         kl += my1;
      kl = kl*mx1;
/* find tile below */
      kr = ky + 1;
      if (kr >= my1)
          kr -= my1;
      kr = kr*mx1;
/* loop over tiles in x, assume periodic boundary conditions */
      kx = k - ky*mx1;
      kxl = kx - 1;
      if (kxl < 0)
         kxl += mx1;
      kxr = kx + 1;
      if (kxr >= mx1)
         kxr -= mx1;
/* find tile number for different directions */
      ks[0] = kxr + kk;
      ks[1] = kxl + kk;
      ks[2] = kx + kr;
      ks[3] = kxr + kr;
      ks[4] = kxl + kr;
      ks[5] = kx + kl;
      ks[6] = kxr + kl;
      ks[7] = kxl + kl;
/* loop over directions */
      nh = ihole[2*(ntmax+1)*k];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      for (ii = 0; ii < 8; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+8*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+8*ks[ii]] - ncoff;
         for (j = 0; j < ip; j++) {
            ih += 1;
/* insert incoming particles into holes */
            if (ih <= nh) {
               j1 = ihole[2*(ih+(ntmax+1)*k)] - 1;
            }
/* place overflow at end of array */
            else {
               j1 = npp;
               npp += 1;
            }
            if (j1 < nppmx) {
               for (i = 0; i < idimp; i++) {
                  ppart[j1+nppmx*i+npoff]
                  = ppbuff[j+ncoff+npbmx*i+nboff];
                }
            }
            else {
               ist = 1;
            }
         }
      }
/* set error */
      if (ist > 0)
         *irc = j1+1;
/* fill up remaining holes in particle array with particles from bottom */
/* holes with locations great than npp-ip do not need to be filled      */
      if (ih < nh) {
         ip = nh - ih;
         ii = nh;
         nn = ihole[2*(ii+(ntmax+1)*k)] - 1;
         ih += 1;
         j2 = ihole[2*(ih+(ntmax+1)*k)] - 1;
/* move particles from end into remaining holes */
/* holes are processed in increasing order      */
         for (j = 0; j < ip; j++) {
            j1 = npp - j - 1;
            if (j1==nn) {
               ii -= 1;
               nn = ihole[2*(ii+(ntmax+1)*k)] - 1;
            }
            else {
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*i+npoff]
                  = ppart[j1+nppmx*i+npoff];
               }
               ih += 1;
               j2 = ihole[2*(ih+(ntmax+1)*k)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[k] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int mx1, int my1,
                  int npbmx, int ntmax, int *irc) {
/* this subroutine sorts particles by x,y grid in tiles of mx, my
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 2D linear memory
   the algorithm has 2 steps.  first, a prefix scan of ncl is performed
   and departing particles are buffered in ppbuff in direction order.
   then we copy the incoming particles from other tiles into ppart.
   it assumes that the number, location, and destination of particles 
   leaving a tile have been previously stored in ncl and ihole by the
   cgppushf2lt procedure.
   input: all except ppbuff, irc
   output: ppart, ppbuff, kpic, ncl, irc
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k 
   ppbuff[k][i][n] = i co-ordinate of particle n in tile k
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = direction destination of particle leaving hole
   all for tile k
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int mxy1, npoff, npp, nboff, ncoff;
   int i, j, k, ii, kx, ky, ih, nh, ist, nn, isum;
   int ip, j1, j2, kxl, kxr, kk, kl, kr;
   int ks[8];
   mxy1 = mx1*my1;
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,npoff,nboff,isum,ist,nh,ip,j1,ii)
   for (k = 0; k < mxy1; k++) {
      npoff = idimp*nppmx*k;
      nboff = idimp*npbmx*k;
/* find address offset for ordered ppbuff array */
      isum = 0;
      for (j = 0; j < 8; j++) {
         ist = ncl[j+8*k];
         ncl[j+8*k] = isum;
         isum += ist;
      }
      nh = ihole[2*(ntmax+1)*k];
      ip = 0;
/* loop over particles leaving tile */
      for (j = 0; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+(ntmax+1)*k)] - 1;
         ist = ihole[1+2*(j+1+(ntmax+1)*k)];
         ii = ncl[ist+8*k-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*i+nboff]
               = ppart[j1+nppmx*i+npoff];
            }
         }
         else {
            ip = 1;
         }
         ncl[ist+8*k-1] = ii + 1;
      }
/* set error */
      if (ip > 0)
         *irc = ncl[7+8*k];
   }
/* ppbuff overflow */
   if (*irc > 0)
      return;

/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,ii,kk,npp,npoff,nboff,kx,ky,kl,kr,kxl,kxr,ih,nh,nn, \
ncoff,ist,j1,j2,ip,ks)
   for (k = 0; k < mxy1; k++) {
      npp = kpic[k];
      npoff = idimp*nppmx*k;
      ky = k/mx1;
/* loop over tiles in y, assume periodic boundary conditions */
      kk = ky*mx1;
/* find tile above */
      kl = ky - 1;
      if (kl < 0)
         kl += my1;
      kl = kl*mx1;
/* find tile below */
      kr = ky + 1;
      if (kr >= my1)
          kr -= my1;
      kr = kr*mx1;
/* loop over tiles in x, assume periodic boundary conditions */
      kx = k - ky*mx1;
      kxl = kx - 1;
      if (kxl < 0)
         kxl += mx1;
      kxr = kx + 1;
      if (kxr >= mx1)
         kxr -= mx1;
/* find tile number for different directions */
      ks[0] = kxr + kk;
      ks[1] = kxl + kk;
      ks[2] = kx + kr;
      ks[3] = kxr + kr;
      ks[4] = kxl + kr;
      ks[5] = kx + kl;
      ks[6] = kxr + kl;
      ks[7] = kxl + kl;
/* loop over directions */
      nh = ihole[2*(ntmax+1)*k];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      for (ii = 0; ii < 8; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+8*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+8*ks[ii]] - ncoff;
         for (j = 0; j < ip; j++) {
            ih += 1;
/* insert incoming particles into holes */
            if (ih <= nh) {
               j1 = ihole[2*(ih+(ntmax+1)*k)] - 1;
            }
/* place overflow at end of array */
            else {
               j1 = npp;
               npp += 1;
            }
            if (j1 < nppmx) {
               for (i = 0; i < idimp; i++) {
                  ppart[j1+nppmx*i+npoff]
                  = ppbuff[j+ncoff+npbmx*i+nboff];
                }
            }
            else {
               ist = 1;
            }
         }
      }
/* set error */
      if (ist > 0)
         *irc = j1+1;
/* fill up remaining holes in particle array with particles from bottom */
/* holes with locations great than npp-ip do not need to be filled      */
      if (ih < nh) {
         ip = nh - ih;
         ii = nh;
         nn = ihole[2*(ii+(ntmax+1)*k)] - 1;
         ih += 1;
         j2 = ihole[2*(ih+(ntmax+1)*k)] - 1;
/* move particles from end into remaining holes */
/* holes are processed in increasing order      */
         for (j = 0; j < ip; j++) {
            j1 = npp - j - 1;
            if (j1==nn) {
               ii -= 1;
               nn = ihole[2*(ii+(ntmax+1)*k)] - 1;
            }
            else {
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*i+npoff]
                  = ppart[j1+nppmx*i+npoff];
               }
               ih += 1;
               j2 = ihole[2*(ih+(ntmax+1)*k)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[k] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cvpporder2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int nx, int ny,
                  int mx, int my, int mx1, int my1, int npbmx,
                  int ntmax, int *irc) {
/* this subroutine sorts particles by x,y grid in tiles of mx, my
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 2D linear memory
   algorithm has 3 steps.  first, one finds particles leaving tile and
   stores their number in each directon, location, and destination in ncl
   and ihole.  second, a prefix scan of ncl is performed and departing
   particles are buffered in ppbuff in direction order.  finally, we copy
   the incoming particles from other tiles into ppart.
   input: all except ppbuff, ncl, ihole, irc
   output: ppart, ppbuff, kpic, ncl, ihole, irc
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k 
   ppbuff[k][i][n] = i co-ordinate of particle n in tile k
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = direction destination of particle leaving hole
   all for tile k
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
#define NPBLK             16
   int mxy1, noff, moff, npoff, npp, ipp, joff, nps, nboff, ncoff;
   int i, j, k, m, ii, kx, ky, ih, nh, ist, nn, mm, in;
   int ip, j1, j2, kxl, kxr, kk, kl, kr, lb, kxs;
   float anx, any, edgelx, edgely, edgerx, edgery, dx, dy;
   int sncl[8], ks[8];
/* scratch arrays */
   int n[NPBLK*3];
   mxy1 = mx1*my1;
   anx = (float) nx;
   any = (float) ny;
/* find and count particles leaving tiles and determine destination */
/* update ppart, ihole, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(j,k,noff,moff,npp,npoff,nn,mm,ih,nh,ist,dx,dy,edgelx,edgely, \
edgerx,edgery)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
      nn = nx - noff;
      nn = mx < nn ? mx : nn;
      mm = ny - moff;
      mm = my < mm ? my : mm;
      ih = 0;
      nh = 0;
      edgelx = noff;
      edgerx = noff + nn;
      edgely = moff;
      edgery = moff + mm;
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
         dx = ppart[j+npoff];
         dy = ppart[j+nppmx+npoff];
/* find particles going out of bounds */
         ist = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* ist = direction particle is going                             */
         if (dx >= edgerx) {
            if (dx >= anx)
               ppart[j+npoff] = dx - anx;
            ist = 2;
         }
         else if (dx < edgelx) {
            if (dx < 0.0) {
               dx += anx;
               if (dx < anx)
                  ist = 1;
               else
                  dx = 0.0;
               ppart[j+npoff] = dx;
            }
            else {
               ist = 1;
            }
         }
         if (dy >= edgery) {
            if (dy >= any)
               ppart[j+nppmx+npoff] = dy - any;
            ist += 6;
         }
         else if (dy < edgely) {
            if (dy < 0.0) {
               dy += any;
               if (dy < any)
                  ist += 3;
               else
                  dy = 0.0;
               ppart[j+nppmx+npoff] = dy;
            }
            else {
               ist += 3;
            }
         }
         if (ist > 0) {
            ncl[ist+8*k-1] += 1;
            ih += 1;
            if (ih <= ntmax) {
               ihole[2*(ih+(ntmax+1)*k)] = j + 1;
               ihole[1+2*(ih+(ntmax+1)*k)] = ist;
            }
            else {
               nh = 1;
            }
         }
      }
/* set error and end of file flag */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*k] = ih;
   }
/* ihole overflow */
   if (*irc > 0)
      return;

/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,kxs,lb,npoff,nboff,ist,nh,ip,ipp,nps,joff,j1,ii,sncl, \
ks,n)
   for (k = 0; k < mxy1; k++) {
      npoff = idimp*nppmx*k;
      nboff = idimp*npbmx*k;
/* find address offset for ordered ppbuff array */
      for (j = 0; j < 8; j++) {
         sncl[j] = ncl[j+8*k];
         ks[j] = j;
      }
      kxs = 1;
      while (kxs < 8) {
#pragma ivdep
         for (j = 0; j < 4; j++) {
            lb = kxs*ks[j];
            sncl[j+lb+kxs] += sncl[2*lb+kxs-1];
            ks[j] >>= 1;
         }     
         kxs <<= 1;
      }
      for (j = 0; j < 8; j++) {
         sncl[j] -= ncl[j+8*k];
      }
      nh = ihole[2*(ntmax+1)*k];
      ip = 0;
/* buffer particles that are leaving tile, in direction order */
/* loop over particles leaving tile */
      ipp = nh/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m + 1;
/* inner loop over particles in block */
         for (j = 0; j < NPBLK; j++) {
            n[j] = ihole[2*(j+joff+(ntmax+1)*k)] - 1;
            n[j+NPBLK] = ihole[1+2*(j+joff+(ntmax+1)*k)];
         }
/* calculate offsets */
         for (j = 0; j < NPBLK; j++) {
            ist = n[j+NPBLK];
            ii = sncl[ist-1];
            n[j+NPBLK] = ii;
            sncl[ist-1] = ii + 1;
         }
/* buffer particles that are leaving tile, in direction order */
         for (i = 0; i < idimp; i++) {
            for (j = 0; j < NPBLK; j++) {
               j1 = n[j];
               ii = n[j+NPBLK];
               if (ii < npbmx) {
                 ppbuff[ii+npbmx*i+nboff]
                  = ppart[j1+nppmx*i+npoff];
               }
               else {
                  ip = 1;
               }
            }
         }
      }
      nps = NPBLK*ipp;
/* loop over remaining particles */
      for (j = nps; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+(ntmax+1)*k)] - 1;
         ist = ihole[1+2*(j+1+(ntmax+1)*k)];
         ii = sncl[ist-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*i+nboff]
               = ppart[j1+nppmx*i+npoff];
            }
         }
         else {
            ip = 1;
         }
         sncl[ist-1] = ii + 1;
      }
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = sncl[j];
      }
/* set error */
      if (ip > 0)
         *irc = ncl[7+8*k];
   }
/* ppbuff overflow */
   if (*irc > 0)
      return;

/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,ii,kk,in,npp,npoff,nboff,ipp,joff,nps,kx,ky,kl,kr,kxl, \
kxr,ih,nh,nn,mm,ncoff,ist,j1,j2,ip,ks,n)
   for (k = 0; k < mxy1; k++) {
      npp = kpic[k];
      npoff = idimp*nppmx*k;
      ky = k/mx1;
/* loop over tiles in y, assume periodic boundary conditions */
      kk = ky*mx1;
/* find tile above */
      kl = ky - 1;
      if (kl < 0)
         kl += my1;
      kl = kl*mx1;
/* find tile below */
      kr = ky + 1;
      if (kr >= my1)
          kr -= my1;
      kr = kr*mx1;
/* loop over tiles in x, assume periodic boundary conditions */
      kx = k - ky*mx1;
      kxl = kx - 1;
      if (kxl < 0)
         kxl += mx1;
      kxr = kx + 1;
      if (kxr >= mx1)
         kxr -= mx1;
/* find tile number for different directions */
      ks[0] = kxr + kk;
      ks[1] = kxl + kk;
      ks[2] = kx + kr;
      ks[3] = kxr + kr;
      ks[4] = kxl + kr;
      ks[5] = kx + kl;
      ks[6] = kxr + kl;
      ks[7] = kxl + kl;
/* loop over directions */
      nh = ihole[2*(ntmax+1)*k];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      for (ii = 0; ii < 8; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+8*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+8*ks[ii]] - ncoff;
/* loop over particles coming from direction ii */
         ipp = ip/NPBLK;
/* outer loop over number of full blocks */
         for (m = 0; m < ipp; m++) {
            joff = NPBLK*m;
/* inner loop over particles in block */
            for (j = 0; j < NPBLK; j++) {
/* insert incoming particles into holes */
               if ((j+ih) < nh) {
                  j1 = ihole[2*(j+ih+1+(ntmax+1)*k)] - 1;
               }
/* place overflow at end of array */
               else {
                  j1 = npp + j + ih - nh;
               }
               n[j] = j1;
            }
            for (i = 0; i < idimp; i++) {
               for (j = 0; j < NPBLK; j++) {
                  j1 = n[j];
                  if (j1 < nppmx) {
                     ppart[j1+nppmx*i+npoff]
                     = ppbuff[j+joff+ncoff+npbmx*i+nboff];
                  }
                  else {
                    ist = 1;
                  }
               }
            }
            ih += NPBLK;
         }
         nps = NPBLK*ipp;
/* loop over remaining particles */
         for (j = nps; j < ip; j++) {
            ih += 1;
/* insert incoming particles into holes */
            if (ih <= nh) {
               j1 = ihole[2*(ih+(ntmax+1)*k)] - 1;
            }
/* place overflow at end of array */
            else {
               j1 = npp + ih - nh - 1;
            }
            if (j1 < nppmx) {
               for (i = 0; i < idimp; i++) {
                  ppart[j1+nppmx*i+npoff]
                  = ppbuff[j+ncoff+npbmx*i+nboff];
                }
            }
            else {
               ist = 1;
            }
         }
      }
      if (ih > nh)
         npp = npp + ih - nh;
/* set error */
      if (ist > 0)
         *irc = j1+1;
/* fill up remaining holes in particle array with particles from bottom */
/* holes with locations great than npp-ip do not need to be filled      */
      if (ih < nh) {
         ip = nh - ih;
/* move particles from end into remaining holes */
/* holes are processed in increasing order      */
         ii = nh;
         ipp = ip/NPBLK;
/* outer loop over number of full blocks */
         for (m = 0; m < ipp; m++) {
            joff = NPBLK*m;
/* inner loop over particles in block */
            for (j = 0; j < NPBLK; j++) {
               n[j+NPBLK] = ihole[2*(ih+j+1+(ntmax+1)*k)] - 1;
               n[j+2*NPBLK] = ihole[2*(ii-j+(ntmax+1)*k)] - 1;
            }
            in = 0;
            mm = 0;
            nn = n[in+2*NPBLK];
            for (j = 0; j < NPBLK; j++) {
               j1 = npp - j - joff - 1;
               n[j] = n[mm+NPBLK];
               if (j1==nn) {
                  in += 1;
                  nn = n[in+2*NPBLK];
                  n[j] = -1;
               }
               else {
                  mm += 1;
               }
            }
            for (i = 0; i < idimp; i++) {
#pragma ivdep
               for (j = 0; j < NPBLK; j++) {
                  j1 = npp - j - joff - 1;
                  j2 = n[j];
                  if (j2 >= 0) {
                     ppart[j2+nppmx*i+npoff]
                     = ppart[j1+nppmx*i+npoff];
                  }
               }
            }
            ii -= in;
            ih += mm;
         }
         nps = NPBLK*ipp;
         nn = ihole[2*(ii+(ntmax+1)*k)] - 1;
         ih += 1;
         j2 = ihole[2*(ih+(ntmax+1)*k)] - 1;
/* loop over remaining particles */
         for (j = nps; j < ip; j++) {
            j1 = npp - j - 1;
            if (j1==nn) {
               ii -= 1;
               nn = ihole[2*(ii+(ntmax+1)*k)] - 1;
            }
            else {
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*i+npoff]
                  = ppart[j1+nppmx*i+npoff];
               }
               ih += 1;
               j2 = ihole[2*(ih+(ntmax+1)*k)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[k] = npp;
   }
   return;
#undef NPBLK
}

/*--------------------------------------------------------------------*/
void cvpporderf2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                   int ihole[], int idimp, int nppmx, int mx1, int my1,
                   int npbmx, int ntmax, int *irc) {
/* this subroutine sorts particles by x,y grid in tiles of mx, my
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 2D linear memory
   the algorithm has 2 steps.  first, a prefix scan of ncl is performed
   and departing particles are buffered in ppbuff in direction order.
   then we copy the incoming particles from other tiles into ppart.
   it assumes that the number, location, and destination of particles 
   leaving a tile have been previously stored in ncl and ihole by the
   cgppushf2lt procedure.
   input: all except ppbuff, irc
   output: ppart, ppbuff, kpic, ncl, irc
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k 
   ppbuff[k][i][n] = i co-ordinate of particle n in tile k
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = direction destination of particle leaving hole
   all for tile k
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
#define NPBLK             16
   int mxy1, npoff, npp, nboff, ncoff;
   int i, j, k, ii, kx, ky, ih, nh, ist, nn, mm, in;
   int ip, j1, j2, kxl, kxr, kk, kl, kr;
   int lb, kxs, m, ipp, nps, joff;
   int sncl[8], ks[8];
/* scratch arrays */
   int n[NPBLK*3];
   mxy1 = mx1*my1;
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,kxs,lb,npoff,nboff,ist,nh,ip,ipp,nps,joff,j1,ii,sncl, \
ks,n)
   for (k = 0; k < mxy1; k++) {
      npoff = idimp*nppmx*k;
      nboff = idimp*npbmx*k;
/* find address offset for ordered ppbuff array */
      for (j = 0; j < 8; j++) {
         sncl[j] = ncl[j+8*k];
         ks[j] = j;
      }
      kxs = 1;
      while (kxs < 8) {
#pragma ivdep
         for (j = 0; j < 4; j++) {
            lb = kxs*ks[j];
            sncl[j+lb+kxs] += sncl[2*lb+kxs-1];
            ks[j] >>= 1;
         }     
         kxs <<= 1;
      }
      for (j = 0; j < 8; j++) {
         sncl[j] -= ncl[j+8*k];
      }
      nh = ihole[2*(ntmax+1)*k];
      ip = 0;
/* buffer particles that are leaving tile, in direction order */
/* loop over particles leaving tile */
      ipp = nh/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m + 1;
/* inner loop over particles in block */
         for (j = 0; j < NPBLK; j++) {
            n[j] = ihole[2*(j+joff+(ntmax+1)*k)] - 1;
            n[j+NPBLK] = ihole[1+2*(j+joff+(ntmax+1)*k)];
         }
/* calculate offsets */
         for (j = 0; j < NPBLK; j++) {
            ist = n[j+NPBLK];
            ii = sncl[ist-1];
            n[j+NPBLK] = ii;
            sncl[ist-1] = ii + 1;
         }
/* buffer particles that are leaving tile, in direction order */
         for (i = 0; i < idimp; i++) {
            for (j = 0; j < NPBLK; j++) {
               j1 = n[j];
               ii = n[j+NPBLK];
               if (ii < npbmx) {
                 ppbuff[ii+npbmx*i+nboff]
                  = ppart[j1+nppmx*i+npoff];
               }
               else {
                  ip = 1;
               }
            }
         }
      }
      nps = NPBLK*ipp;
/* loop over remaining particles */
      for (j = nps; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+(ntmax+1)*k)] - 1;
         ist = ihole[1+2*(j+1+(ntmax+1)*k)];
         ii = sncl[ist-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*i+nboff]
               = ppart[j1+nppmx*i+npoff];
            }
         }
         else {
            ip = 1;
         }
         sncl[ist-1] = ii + 1;
      }
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = sncl[j];
      }
/* set error */
      if (ip > 0)
         *irc = ncl[7+8*k];
   }
/* ppbuff overflow */
   if (*irc > 0)
      return;

/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,ii,kk,in,npp,npoff,nboff,ipp,joff,nps,kx,ky,kl,kr,kxl, \
kxr,ih,nh,nn,mm,ncoff,ist,j1,j2,ip,ks,n)
   for (k = 0; k < mxy1; k++) {
      npp = kpic[k];
      npoff = idimp*nppmx*k;
      ky = k/mx1;
/* loop over tiles in y, assume periodic boundary conditions */
      kk = ky*mx1;
/* find tile above */
      kl = ky - 1;
      if (kl < 0)
         kl += my1;
      kl = kl*mx1;
/* find tile below */
      kr = ky + 1;
      if (kr >= my1)
          kr -= my1;
      kr = kr*mx1;
/* loop over tiles in x, assume periodic boundary conditions */
      kx = k - ky*mx1;
      kxl = kx - 1;
      if (kxl < 0)
         kxl += mx1;
      kxr = kx + 1;
      if (kxr >= mx1)
         kxr -= mx1;
/* find tile number for different directions */
      ks[0] = kxr + kk;
      ks[1] = kxl + kk;
      ks[2] = kx + kr;
      ks[3] = kxr + kr;
      ks[4] = kxl + kr;
      ks[5] = kx + kl;
      ks[6] = kxr + kl;
      ks[7] = kxl + kl;
/* loop over directions */
      nh = ihole[2*(ntmax+1)*k];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      for (ii = 0; ii < 8; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+8*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+8*ks[ii]] - ncoff;
/* loop over particles coming from direction ii */
         ipp = ip/NPBLK;
/* outer loop over number of full blocks */
         for (m = 0; m < ipp; m++) {
            joff = NPBLK*m;
/* inner loop over particles in block */
            for (j = 0; j < NPBLK; j++) {
/* insert incoming particles into holes */
               if ((j+ih) < nh) {
                  j1 = ihole[2*(j+ih+1+(ntmax+1)*k)] - 1;
               }
/* place overflow at end of array */
               else {
                  j1 = npp + j + ih - nh;
               }
               n[j] = j1;
            }
            for (i = 0; i < idimp; i++) {
               for (j = 0; j < NPBLK; j++) {
                  j1 = n[j];
                  if (j1 < nppmx) {
                     ppart[j1+nppmx*i+npoff]
                     = ppbuff[j+joff+ncoff+npbmx*i+nboff];
                  }
                  else {
                    ist = 1;
                  }
               }
            }
            ih += NPBLK;
         }
         nps = NPBLK*ipp;
/* loop over remaining particles */
         for (j = nps; j < ip; j++) {
            ih += 1;
/* insert incoming particles into holes */
            if (ih <= nh) {
               j1 = ihole[2*(ih+(ntmax+1)*k)] - 1;
            }
/* place overflow at end of array */
            else {
               j1 = npp + ih - nh - 1;
            }
            if (j1 < nppmx) {
               for (i = 0; i < idimp; i++) {
                  ppart[j1+nppmx*i+npoff]
                  = ppbuff[j+ncoff+npbmx*i+nboff];
                }
            }
            else {
               ist = 1;
            }
         }
      }
      if (ih > nh)
         npp = npp + ih - nh;
/* set error */
      if (ist > 0)
         *irc = j1+1;
/* fill up remaining holes in particle array with particles from bottom */
/* holes with locations great than npp-ip do not need to be filled      */
      if (ih < nh) {
         ip = nh - ih;
/* move particles from end into remaining holes */
/* holes are processed in increasing order      */
         ii = nh;
         ipp = ip/NPBLK;
/* outer loop over number of full blocks */
         for (m = 0; m < ipp; m++) {
            joff = NPBLK*m;
/* inner loop over particles in block */
            for (j = 0; j < NPBLK; j++) {
               n[j+NPBLK] = ihole[2*(ih+j+1+(ntmax+1)*k)] - 1;
               n[j+2*NPBLK] = ihole[2*(ii-j+(ntmax+1)*k)] - 1;
            }
            in = 0;
            mm = 0;
            nn = n[in+2*NPBLK];
            for (j = 0; j < NPBLK; j++) {
               j1 = npp - j - joff - 1;
               n[j] = n[mm+NPBLK];
               if (j1==nn) {
                  in += 1;
                  nn = n[in+2*NPBLK];
                  n[j] = -1;
               }
               else {
                  mm += 1;
               }
            }
            for (i = 0; i < idimp; i++) {
#pragma ivdep
               for (j = 0; j < NPBLK; j++) {
                  j1 = npp - j - joff - 1;
                  j2 = n[j];
                  if (j2 >= 0) {
                     ppart[j2+nppmx*i+npoff]
                     = ppart[j1+nppmx*i+npoff];
                  }
               }
            }
            ii -= in;
            ih += mm;
         }
         nps = NPBLK*ipp;
         nn = ihole[2*(ii+(ntmax+1)*k)] - 1;
         ih += 1;
         j2 = ihole[2*(ih+(ntmax+1)*k)] - 1;
/* loop over remaining particles */
         for (j = nps; j < ip; j++) {
            j1 = npp - j - 1;
            if (j1==nn) {
               ii -= 1;
               nn = ihole[2*(ii+(ntmax+1)*k)] - 1;
            }
            else {
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*i+npoff]
                  = ppart[j1+nppmx*i+npoff];
               }
               ih += 1;
               j2 = ihole[2*(ih+(ntmax+1)*k)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[k] = npp;
   }
   return;
#undef NPBLK
}

/*--------------------------------------------------------------------*/
void ccguard2l(float fxy[], int nx, int ny, int nxe, int nye) {
/* replicate extended periodic vector field fxy
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nxe = second dimension of field arrays, must be >= ny+1
local data                                                 */
   int j, k;
/* copy edges of extended field */
   for (k = 0; k < ny; k++) {
      fxy[2*nx+2*nxe*k] = fxy[2*nxe*k];
      fxy[1+2*nx+2*nxe*k] = fxy[1+2*nxe*k];
   }
   for (j = 0; j < nx; j++) {
      fxy[2*j+2*nxe*ny] = fxy[2*j];
      fxy[1+2*j+2*nxe*ny] = fxy[1+2*j];
   }
   fxy[2*nx+2*nxe*ny] = fxy[0];
   fxy[1+2*nx+2*nxe*ny] = fxy[1];
   return;
}

/*--------------------------------------------------------------------*/
void caguard2l(float q[], int nx, int ny, int nxe, int nye) {
/* accumulate extended periodic scalar field q
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nxe = second dimension of field arrays, must be >= ny+1
local data                                                 */
   int j, k;
/* accumulate edges of extended field */
   for (k = 0; k < ny; k++) {
      q[nxe*k] += q[nx+nxe*k];
      q[nx+nxe*k] = 0.0;
   }
   for (j = 0; j < nx; j++) {
      q[j] += q[j+nxe*ny];
      q[j+nxe*ny] = 0.0;
   }
   q[0] += q[nx+nxe*ny];
   q[nx+nxe*ny] = 0.0;
   return;
}

/*--------------------------------------------------------------------*/
void cvmpois22(float complex q[], float complex fxy[], int isign,
               float complex ffc[], float ax, float ay, float affp,
               float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
               int nyhd) {
/* this subroutine solves 2d poisson's equation in fourier space for
   force/charge (or convolution of electric field over particle shape)
   with periodic boundary conditions.
   for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
   for isign /= 0, input: q,ffc,isign,nx,ny,nxvh,nyhd, output: fxy,we
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
   q[k][j] = complex charge density for fourier mode (j,k)
   fxy[k][j][0] = x component of complex force/charge,
   fxy[k][j][1] = y component of complex force/charge,
   all for fourier mode (j,k)
   if isign = 0, form factor array is prepared
   if isign is not equal to 0, force/charge is calculated
   cimag(ffc[k][j]) = finite-size particle shape factor s
   for fourier mode (j,k)
   creal(ffc[k][j]) = potential green's function g
   for fourier mode (j,k)
   ax/ay = half-width of particle in x/y direction
   affp = normalization constant = nx*ny/np, where np=number of particles
   electric field energy is also calculated, using
   we = nx*ny*sum((affp/(kx**2+ky**2))*|q[ky][kx]*s[ky][kx]|**2)
   nx/ny = system length in x/y direction
   nxvh = first dimension of field arrays, must be >= nxh
   nyv = second dimension of field arrays, must be >= ny
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
   vectorizable version
local data                                                 */
   int nxh, nyh, j, k, k1, kk, kj;
   float dnx, dny, dkx, dky, at1, at2, at3, at4;
   float complex zero, zt1, zt2;
   double wp, sum1;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
   if (isign != 0)
      goto L30;
/* prepare form factor array */
   for (k = 0; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      at1 = dky*dky;
      at2 = pow((dky*ay),2);
      for (j = 0; j < nxh; j++) {
         dkx = dnx*(float) j;
         at3 = dkx*dkx + at1;
         at4 = exp(-0.5*(pow((dkx*ax),2) + at2));
         if (at3==0.0) {
            ffc[j+kk] = affp + 1.0*_Complex_I;
         }
         else {
            ffc[j+kk] = (affp*at4/at3) + at4*_Complex_I;
         }
      }
   }
   return;
/* calculate force/charge and sum field energy */
L30: sum1 = 0.0;
#pragma omp parallel for \
private(j,k,k1,kk,kj,dky,at1,at2,at3,zt1,zt2,wp) \
reduction(+:sum1)
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
      k1 = ny - k;
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      wp = 0.0;
#pragma ivdep
      for (j = 1; j < nxh; j++) {
         at1 = crealf(ffc[j+kk])*cimagf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         zt1 = cimagf(q[j+kj]) - crealf(q[j+kj])*_Complex_I;
         zt2 = cimagf(q[j+k1]) - crealf(q[j+k1])*_Complex_I;
         fxy[2*j+2*kj] = at2*zt1;
         fxy[1+2*j+2*kj] = at3*zt1;
         fxy[2*j+2*k1] = at2*zt2;
         fxy[1+2*j+2*k1] = -at3*zt2;
         at1 = at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1]));
         wp += (double) at1;
      }
      at1 = crealf(ffc[kk])*cimagf(ffc[kk]);
      at3 = at1*dny*(float) k;
      zt1 = cimagf(q[kj]) - crealf(q[kj])*_Complex_I;
      fxy[2*kj] = zero;
      fxy[1+2*kj] = at3*zt1;
      fxy[2*k1] = zero;
      fxy[1+2*k1] = zero;
      at1 = at1*(q[kj]*conjf(q[kj]));
      wp += (double) at1;
      sum1 += wp;
   }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = 2*nxvh*nyh;
#pragma ivdep
   for (j = 1; j < nxh; j++) {
      at1 = crealf(ffc[j])*cimagf(ffc[j]);
      at2 = at1*dnx*(float) j;  
      zt1 = cimagf(q[j]) - crealf(q[j])*_Complex_I;
      fxy[2*j] = at2*zt1;
      fxy[1+2*j] = zero;
      fxy[2*j+k1] = zero;
      fxy[1+2*j+k1] = zero;
      at1 = at1*(q[j]*conjf(q[j]));
      wp += (double) at1;
   }
   fxy[0] = zero;
   fxy[1] = zero;
   fxy[k1] = zero;
   fxy[1+k1] = zero;
   sum1 += wp;
   *we = sum1*(float) (nx*ny);
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

/*--------------------------------------------------------------------*/
void cfft2rvmxx(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nyi,
                int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the x part of a two dimensional real to
   complex fast fourier transform and its inverse, for a subset of y,
   using complex arithmetic, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, an inverse fourier transform in x is performed
   f[m][n] = (1/nx*ny)*sum(f[k][j]*exp(-sqrt(-1)*2pi*n*j/nx))
   if isign = 1, a forward fourier transform in x is performed
   f[k][j] = sum(f[m][n]*exp(sqrt(-1)*2pi*n*j/nx))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = first dimension of f >= nx/2
   nyd = second dimension of f >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0]) = real part of mode nx/2,0 and
   imag(f[0][ny/2]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, j1, k1, k2, ns, ns2, km, kmr, nrxb, joff;
   float ani;
   float complex t1, t2, t3;
   if (isign==0)
      return;
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   nxh = nx/2;
   nxhh = nx/4;
   ny = 1L<<indy;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   nyt = nyi + nyp - 1;
   if (isign > 0)
      goto L70;
/* inverse fourier transform */
   nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,joff,ani,t1,t2,t3)
   for (i = nyi-1; i < nyt; i++) {
      joff = nxhd*i;
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            t1 = f[j1+joff];
            f[j1+joff] = f[j+joff];
            f[j+joff] = t1;
         }
      }
/* then transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               t1 = sct[kmr*j];
               t2 = t1*f[j+k2+joff];
               f[j+k2+joff] = f[j+k1+joff] - t2;
               f[j+k1+joff] += t2;
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxy/nx;
      ani = 0.5/(((float) nx)*((float) ny));
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         t2 = conjf(f[nxh-j+joff]);
         t1 = f[j+joff] + t2;
         t2 = (f[j+joff] - t2)*t3;
         f[j+joff] = ani*(t1 + t2);
         f[nxh-j+joff] = ani*conjf(t1 - t2);
      }
      ani = 2.0*ani;
      f[nxhh+joff] = ani*conjf(f[nxhh+joff]);
      f[joff] = ani*((crealf(f[joff]) + cimagf(f[joff]))
                + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I);
   }
   return;
/* forward fourier transform */
L70: nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,joff,t1,t2,t3)
   for (i = nyi-1; i < nyt; i++) {
      joff = nxhd*i;
/* scramble coefficients */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         t2 = conjf(f[nxh-j+joff]);
         t1 = f[j+joff] + t2;
         t2 = (f[j+joff] - t2)*t3;
         f[j+joff] = t1 + t2;
         f[nxh-j+joff] = conjf(t1 - t2);
      }
      f[nxhh+joff] = 2.0*conjf(f[nxhh+joff]);
      f[joff] = (crealf(f[joff]) + cimagf(f[joff]))
                + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I;
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            t1 = f[j1+joff];
            f[j1+joff] = f[j+joff];
            f[j+joff] = t1;
         }
      }
/* then transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               t1 = conjf(sct[kmr*j]);
               t2 = t1*f[j+k2+joff];
               f[j+k2+joff] = f[j+k1+joff] - t2;
               f[j+k1+joff] += t2;
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the y part of a two dimensional real to
   complex fast fourier transform and its inverse, for a subset of x,
   using complex arithmetic, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, an inverse fourier transform in y is performed
   f[m][n] = sum(f[k][j]*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, a forward fourier transform in y is performed
   f[k][j] = sum(f[m][n]*exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxi = initial x index used
   nxp = number of x indices used
   nxhd = first dimension of f >= nx/2
   nyd = second dimension of f >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0]) = real part of mode nx/2,0 and
   imag(f[0][ny/2]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, nryb, koff;
   float complex t1, t2;
   if (isign==0)
      return;
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   ny = 1L<<indy;
   nyh = ny/2;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   nxt = nxi + nxp - 1;
   if (isign > 0)
      goto L70;
/* inverse fourier transform */
   nryb = nxhy/ny;
   nry = nxy/ny;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,koff,t1,t2)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nxhd*k1;
            t1 = f[i+k1];
            f[i+k1] = f[i+koff];
            f[i+koff] = t1;
         }
      }
/* then transform in y */
      ns = 1;
      for (l = 0; l < indy; l++) {
         ns2 = ns + ns;
         km = nyh/ns;
         kmr = km*nry;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = nxhd*(j + k1);
               j2 = nxhd*(j + k2);
               t1 = sct[kmr*j];
               t2 = t1*f[i+j2];
               f[i+j2] = f[i+j1] - t2;
               f[i+j1] += t2;
            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if (nxi==1) {
      for (k = 1; k < nyh; k++) {
         koff = nxhd*k;
         k1 = nxhd*ny - koff;
         t1 = f[k1];
         f[k1] = 0.5*(cimagf(f[koff] + t1)
                  + crealf(f[koff] - t1)*_Complex_I);
         f[koff] = 0.5*(crealf(f[koff] + t1)
                    + cimagf(f[koff] - t1)*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
L70: nryb = nxhy/ny;
   nry = nxy/ny;
/* scramble modes kx = 0, nx/2 */
   if (nxi==1) {
      for (k = 1; k < nyh; k++) {
         koff = nxhd*k;
         k1 = nxhd*ny - koff;
         t1 = cimagf(f[k1]) + crealf(f[k1])*_Complex_I;
         f[k1] = conjf(f[koff] - t1);
         f[koff] += t1;
      }
   }
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,koff,t1,t2)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nxhd*k1;
            t1 = f[i+k1];
            f[i+k1] = f[i+koff];
            f[i+koff] = t1;
         }
      }
/* then transform in y */
      ns = 1;
      for (l = 0; l < indy; l++) {
         ns2 = ns + ns;
         km = nyh/ns;
         kmr = km*nry;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = nxhd*(j + k1);
               j2 = nxhd*(j + k2);
               t1 = conjf(sct[kmr*j]);
               t2 = t1*f[i+j2];
               f[i+j2] = f[i+j1] - t2;
               f[i+j1] += t2;
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rvm2x(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nyi,
                int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the x part of 2 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   y, using complex arithmetic, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, two inverse fourier transforms in x are performed
   f[m][n][0:1] = (1/nx*ny)*sum(f[k][j][0:1]*exp(-sqrt(-1)*2pi*n*j/nx))
   if isign = 1, two forward fourier transforms in x are performed
   f[k][j][0:1] = sum(f[m][n][0:1]*exp(sqrt(-1)*2pi*n*j/nx))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = second dimension of f >= nx/2
   nyd = third dimension of f >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j][0:1] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][1][0:1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0][0:1]) = real part of mode nx/2,0 and
   imag(f[0][ny/2][0:1]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, jj, j1, k1, k2, ns, ns2, km, kmr, joff;
   int nrxb;
   float at1, ani;
   float complex t1, t2, t3;
   if (isign==0)
      return;
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   nxh = nx/2;
   nxhh = nx/4;
   ny = 1L<<indy;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   nyt = nyi + nyp - 1;
   if (isign > 0)
      goto L100;
/* inverse fourier transform */
   nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,joff,at1,ani,t1,t2,t3)
   for (i = nyi-1; i < nyt; i++) {
      joff = 2*nxhd*i;
/* swap complex components */
      for (j = 0; j < nxh; j++) {
         at1 = cimagf(f[2*j+joff]);
         f[2*j+joff] = crealf(f[2*j+joff])
                       + crealf(f[1+2*j+joff])*_Complex_I;
         f[1+2*j+joff] = at1 + cimagf(f[1+2*j+joff])*_Complex_I;
       }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            t1 = f[2*j1+joff];
            t2 = f[1+2*j1+joff];
            f[2*j1+joff] = f[2*j+joff];
            f[1+2*j1+joff] = f[1+2*j+joff];
            f[2*j+joff] = t1;
            f[1+2*j+joff] = t2;
         }
      }
/* then transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = 2*ns2*k;
            k2 = k1 + 2*ns;
            for (j = 0; j < ns; j++) {
               t1 = sct[kmr*j];
               t2 = t1*f[2*j+k2+joff];
               t3 = t1*f[1+2*j+k2+joff];
               f[2*j+k2+joff] = f[2*j+k1+joff] - t2;
               f[1+2*j+k2+joff] = f[1+2*j+k1+joff] - t3;
               f[2*j+k1+joff] += t2;
               f[1+2*j+k1+joff] += t3;
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxy/nx;
      ani = 0.5/(((float) nx)*((float) ny));
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         for (jj = 0; jj < 2; jj++) {
            t2 = conjf(f[jj+2*(nxh-j)+joff]);
            t1 = f[jj+2*j+joff] + t2;
            t2 = (f[jj+2*j+joff] - t2)*t3;
            f[jj+2*j+joff] = ani*(t1 + t2);
            f[jj+2*(nxh-j)+joff] = ani*conjf(t1 - t2);
         }
      }
      ani = 2.0*ani;
      for (jj = 0; jj < 2; jj++) {
         f[jj+2*nxhh+joff] = ani*conjf(f[jj+2*nxhh+joff]);
         f[jj+joff] = ani*((crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
L100: nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,joff,at1,t1,t2,t3)
   for (i = nyi-1; i < nyt; i++) {
      joff = 2*nxhd*i;
/* scramble coefficients */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         for (jj = 0; jj < 2; jj++) {
            t2 = conjf(f[jj+2*(nxh-j)+joff]);
            t1 = f[jj+2*j+joff] + t2;
            t2 = (f[jj+2*j+joff] - t2)*t3;
            f[jj+2*j+joff] = t1 + t2;
            f[jj+2*(nxh-j)+joff] = conjf(t1 - t2);
         }
      }
      for (jj = 0; jj < 2; jj++) {
         f[jj+2*nxhh+joff] = 2.0*conjf(f[jj+2*nxhh+joff]);
         f[jj+joff] = (crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I;
      }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            t1 = f[2*j1+joff];
            t2 = f[1+2*j1+joff];
            f[2*j1+joff] = f[2*j+joff];
            f[1+2*j1+joff] = f[1+2*j+joff];
            f[2*j+joff] = t1;
            f[1+2*j+joff] = t2;
         }
      }
/* then transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = 2*ns2*k;
            k2 = k1 + 2*ns;
            for (j = 0; j < ns; j++) {
               t1 = conjf(sct[kmr*j]);
               t2 = t1*f[2*j+k2+joff];
               t3 = t1*f[1+2*j+k2+joff];
               f[2*j+k2+joff] = f[2*j+k1+joff] - t2;
               f[1+2*j+k2+joff] = f[1+2*j+k1+joff] - t3;
               f[2*j+k1+joff] += t2;
               f[1+2*j+k1+joff] +=  t3;
            }
         }
         ns = ns2;
      }
/* swap complex components */
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
void cfft2rm2y(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the y part of 2 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   x, using complex arithmetic, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, two inverse fourier transforms in y are performed
   f[m][n][0:1] = *sum(f[k][j][0:1]*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, two forward fourier transforms in y are performed
   f[k][j][0:1] = sum(f[m][n][0:1]*exp(sqrt(-1)*2pi*n*j/nx))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxi = initial x index used
   nxp = number of x indices used
   nxhd = second dimension of f >= nx/2
   nyd = third dimension of f >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j][0:1] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][1][0:1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0][0:1]) = real part of mode nx/2,0 and
   imag(f[0][ny/2][0:1]) = real part of mode nx/2,ny/2
  written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, koff;
   int nryb;
   float complex t1, t2, t3;
   if (isign==0)
      return;
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   ny = 1L<<indy;
   nyh = ny/2;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   nxt = nxi + nxp - 1;
   if (isign > 0)
      goto L80;
/* inverse fourier transform */
   nryb = nxhy/ny;
   nry = nxy/ny;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,koff,t1,t2,t3)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = 2*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = 2*nxhd*k1;
            t1 = f[2*i+k1];
            t2 = f[1+2*i+k1];
            f[2*i+k1] = f[2*i+koff];
            f[1+2*i+k1] = f[1+2*i+koff];
            f[2*i+koff] = t1;
            f[1+2*i+koff] = t2;
         }
      }
/* then transform in y */
      ns = 1;
      for (l = 0; l < indy; l++) {
         ns2 = ns + ns;
         km = nyh/ns;
         kmr = km*nry;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = 2*nxhd*(j + k1);
               j2 = 2*nxhd*(j + k2);
               t1 = sct[kmr*j];
               t2 = t1*f[2*i+j2];
               t3 = t1*f[1+2*i+j2];
               f[2*i+j2] = f[2*i+j1] - t2;
               f[1+2*i+j2] = f[1+2*i+j1] - t3;
               f[2*i+j1] += t2;
               f[1+2*i+j1] += t3;
            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if (nxi==1) {
      for (k = 1; k < nyh; k++) {
         koff = 2*nxhd*k;
         k1 = 2*nxhd*ny - koff;
         for (jj = 0; jj < 2; jj++) {
            t1 = f[jj+k1];
            f[jj+k1] = 0.5*(cimagf(f[jj+koff] + t1)
                        + crealf(f[jj+koff] - t1)*_Complex_I);
            f[jj+koff] = 0.5*(crealf(f[jj+koff] + t1)
                          + cimagf(f[jj+koff] - t1)*_Complex_I);
         }
      }
   }
   return;
/* forward fourier transform */
L80: nryb = nxhy/ny;
   nry = nxy/ny;
/* scramble modes kx = 0, nx/2 */
   if (nxi==1) {
      for (k = 1; k < nyh; k++) {
         koff = 2*nxhd*k;
         k1 = 2*nxhd*ny - koff;
         for (jj = 0; jj < 2; jj++) {
            t1 = cimagf(f[jj+k1]) + crealf(f[jj+k1])*_Complex_I;
            f[jj+k1] = conjf(f[jj+koff] - t1);
            f[jj+koff] += t1;
         }
      }
   }
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,koff,t1,t2,t3)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = 2*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = 2*nxhd*k1;
            t1 = f[2*i+k1];
            t2 = f[1+2*i+k1];
            f[2*i+k1] = f[2*i+koff];
            f[1+2*i+k1] = f[1+2*i+koff];
            f[2*i+koff] = t1;
            f[1+2*i+koff] = t2;
         }
      }
/* then transform in y */
      ns = 1;
      for (l = 0; l < indy; l++) {
         ns2 = ns + ns;
         km = nyh/ns;
         kmr = km*nry;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = 2*nxhd*(j + k1);
               j2 = 2*nxhd*(j + k2);
               t1 = conjf(sct[kmr*j]);
               t2 = t1*f[2*i+j2];
               t3 = t1*f[1+2*i+j2];
               f[2*i+j2] = f[2*i+j1] - t2;
               f[1+2*i+j2] = f[1+2*i+j1] - t3;
               f[2*i+j1] += t2;
               f[1+2*i+j1] += t3;
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rvmx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd) {
/* wrapper function for real to complex fft, with packed data */
/* parallelized with OpenMP */
/* local data */
   int nxh, ny;
   static int nxi = 1, nyi = 1;
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   ny = 1L<<indy;
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cfft2rvmxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                 nxyhd);
/* perform y fft */
      cfft2rmxy(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfft2rmxy(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                nxyhd);
/* perform x fft */
      cfft2rvmxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                 nxyhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rvm2(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nxhd,
                int nyd, int nxhyd, int nxyhd) {
/* wrapper function for 2 2d real to complex ffts, with packed data */
/* parallelized with OpenMP */
/* local data */
   int nxh, ny;
   static int nxi = 1, nyi = 1;
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   ny = 1L<<indy;
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cfft2rvm2x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                 nxyhd);
/* perform y fft */
      cfft2rm2y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfft2rm2y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                nxyhd);
/* perform x fft */
      cfft2rvm2x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                 nxyhd);
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
void cppmovin2ltp_(float *part, float *ppart, int *kpic, int *kp,
                   int *nppmx, int *idimp, int *nop, int *mx, int *my,
                   int *mx1, int *mxy1, int *irc) {
   cppmovin2ltp(part,ppart,kpic,kp,*nppmx,*idimp,*nop,*mx,*my,*mx1,
                *mxy1,irc);
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
void cgppush2lt_(float *ppart, float *fxy, int *kpic, float *qbm,
                 float *dt, float *ek, int *idimp, int *nppmx, int *nx,
                 int *ny, int *mx, int *my, int *nxv, int *nyv,
                 int *mx1, int *mxy1, int *ipbc) {
   cgppush2lt(ppart,fxy,kpic,*qbm,*dt,ek,*idimp,*nppmx,*nx,*ny,*mx,*my,
              *nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppushf2lt_(float *ppart, float *fxy, int *kpic, int *ncl,
                  int *ihole, float *qbm, float *dt, float *ek,
                  int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                  int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                  int *ntmax, int *irc) {
   cgppushf2lt(ppart,fxy,kpic,ncl,ihole,*qbm,*dt,ek,*idimp,*nppmx,*nx,
               *ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppush2lt_(float *ppart, float *fxy, int *kpic, float *qbm,
                  float *dt, float *ek, int *idimp, int *nppmx, int *nx,
                  int *ny, int *mx, int *my, int *nxv, int *nyv,
                  int *mx1, int *mxy1, int *ipbc) {
   cvgppush2lt(ppart,fxy,kpic,*qbm,*dt,ek,*idimp,*nppmx,*nx,*ny,*mx,*my,
               *nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppushf2lt_(float *ppart, float *fxy, int *kpic, int *ncl,
                   int *ihole, float *qbm, float *dt, float *ek,
                   int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                   int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                   int *ntmax, int *irc) {
   cvgppushf2lt(ppart,fxy,kpic,ncl,ihole,*qbm,*dt,ek,*idimp,*nppmx,*nx,
                *ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppost2lt_(float *ppart, float *q, int *kpic, float *qm,
                 int *nppmx, int *idimp, int *mx, int *my, int *nxv,
                 int *nyv, int *mx1, int *mxy1) {
   cgppost2lt(ppart,q,kpic,*qm,*nppmx,*idimp,*mx,*my,*nxv,*nyv,*mx1,
              *mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppost2lt_(float *ppart, float *q, int *kpic, float *qm,
                  int *nppmx, int *idimp, int *mx, int *my, int *nxv,
                  int *nyv, int *mx1, int *mxy1) {
   cvgppost2lt(ppart,q,kpic,*qm,*nppmx,*idimp,*mx,*my,*nxv,*nyv,*mx1,
               *mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cviscan2_(int *isdata, int *mb, int *nths) {
   cviscan2(isdata,mb,*nths);
   return;
}

/*--------------------------------------------------------------------*/
void cpporder2lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                  int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                  int *mx, int *my, int *mx1, int *my1, int *npbmx,
                  int *ntmax, int *irc) {
   cpporder2lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*nx,*ny,*mx,
                *my,*mx1,*my1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf2lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                   int *ihole, int *idimp, int *nppmx, int *mx1,
                   int *my1, int *npbmx, int *ntmax, int *irc) {
   cpporderf2lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*mx1,*my1,
                 *npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvpporder2lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                   int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                   int *mx, int *my, int *mx1, int *my1, int *npbmx,
                   int *ntmax, int *irc) {
   cvpporder2lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*nx,*ny,*mx,
                *my,*mx1,*my1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvpporderf2lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                    int *ihole, int *idimp, int *nppmx, int *mx1,
                    int *my1, int *npbmx, int *ntmax, int *irc) {
   cvpporderf2lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*mx1,*my1,
                 *npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void ccguard2l_(float *fxy, int *nx, int *ny, int *nxe, int *nye) {
   ccguard2l(fxy,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void caguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye) {
   caguard2l(q,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void cvmpois22_(float complex *q, float complex *fxy, int *isign,
                float complex *ffc, float *ax, float *ay, float *affp,
                float *we, int *nx, int *ny, int *nxvh, int *nyv,
                int *nxhd, int *nyhd) {
   cvmpois22(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*nxvh,*nyv,*nxhd,
             *nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                  int *nxhyd, int *nxyhd) {
   cwfft2rinit(mixup,sct,*indx,*indy,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rvmxx_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *nyi,
                 int *nyp, int *nxhd, int *nyd, int *nxhyd,
                 int *nxyhd) {
   cfft2rvmxx(f,*isign,mixup,sct,*indx,*indy,*nyi,*nyp,*nxhd,*nyd,
              *nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmxy_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nxi,
                int *nxp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd) {
   cfft2rmxy(f,*isign,mixup,sct,*indx,*indy,*nxi,*nxp,*nxhd,*nyd,*nxhyd,
             *nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rvm2x_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *nyi,
                 int *nyp, int *nxhd, int *nyd, int *nxhyd,
                 int *nxyhd) {
   cfft2rvm2x(f,*isign,mixup,sct,*indx,*indy,*nyi,*nyp,*nxhd,*nyd,
              *nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rm2y_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nxi,
                int *nxp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd) {
   cfft2rm2y(f,*isign,mixup,sct,*indx,*indy,*nxi,*nxp,*nxhd,*nyd,*nxhyd,
             *nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rvmx_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *nxhd,
                 int *nyd, int *nxhyd, int *nxyhd) {
   cwfft2rvmx(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rvm2_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *nxhd,
                 int *nyd, int *nxhyd, int *nxyhd) {
   cwfft2rvm2(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,*nxyhd);
   return;
}

