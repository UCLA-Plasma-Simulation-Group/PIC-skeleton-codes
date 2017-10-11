/* C Library for Skeleton 3D Electrostatic OpenMP/Vector PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "vmpush3.h"

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
void cdistr3(float part[], float vtx, float vty, float vtz, float vdx,
             float vdy, float vdz, int npx, int npy, int npz, int idimp,
             int nop, int nx, int ny, int nz, int ipbc) {
/* for 3d code, this subroutine calculates initial particle co-ordinates
   and velocities with uniform density and maxwellian velocity with drift
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = position z of particle n
   part[n][3] = velocity vx of particle n
   part[n][4] = velocity vy of particle n
   part[n][5] = velocity vz of particle n
   vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
   vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
   npx/npy/npz = initial number of particles distributed in x/y/z
   direction
   idimp = size of phase space = 6
   nop = number of particles
   nx/ny/nz = system length in x/y/z direction
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
   ranorm = gaussian random number with zero mean and unit variance
local data                                                            */
   int j, k, l, k1, l1, npxy, npxyz;
   float edgelx, edgely, edgelz, at1, at2, at3, at4, at5;
   float sum1, sum2, sum3;
   double dsum1, dsum2, dsum3;
   npxy = npx*npy;
   npxyz = npxy*npz;
/* set boundary values */
   edgelx = 0.0;
   edgely = 0.0;
   edgelz = 0.0;
   at1 = (float) nx/(float) npx;
   at2 = (float) ny/(float) npy;
   at3 = (float) nz/(float) npz;
   if (ipbc==2) {
      edgelx = 1.0;
      edgely = 1.0;
      edgelz = 1.0;
      at1 = (float) (nx-2)/(float) npx;
      at2 = (float) (ny-2)/(float) npy;
      at3 = (float) (nz-2)/(float) npz;
   }
   else if (ipbc==3) {
      edgelx = 1.0;
      edgely = 1.0;
      edgelz = 0.0;
      at1 = (float) (nx-2)/(float) npx;
      at2 = (float) (ny-2)/(float) npy;
   }
/* uniform density profile */
   for (l = 0; l < npz; l++) {
      l1 = idimp*npxy*l;
      at5 = edgelz + at3*(((float) l) + 0.5);
      for (k = 0; k < npy; k++) {
         k1 = idimp*npx*k + l1;
         at4 = edgely + at2*(((float) k) + 0.5);
         for (j = 0; j < npx; j++) {
            part[idimp*j+k1] = edgelx + at1*(((float) j) + 0.5);
            part[1+idimp*j+k1] = at4;
            part[2+idimp*j+k1] = at5;
         }
      }
   }
/* maxwellian velocity distribution */
   for (j = 0; j < npxyz; j++) {
      part[3+idimp*j] = vtx*ranorm();
      part[4+idimp*j] = vty*ranorm();
      part[5+idimp*j] = vtz*ranorm();
   }
/* add correct drift */
   dsum1 = 0.0;
   dsum2 = 0.0;
   dsum3 = 0.0;
   for (j = 0; j < npxyz; j++) {
      dsum1 += part[3+idimp*j];
      dsum2 += part[4+idimp*j];
      dsum3 += part[5+idimp*j];
   }
   sum1 = dsum1;
   sum2 = dsum2;
   sum3 = dsum3;
   at1 = 1.0/(float) npxyz;
   sum1 = at1*sum1 - vdx;
   sum2 = at1*sum2 - vdy;
   sum3 = at1*sum3 - vdz;
   for (j = 0; j < npxyz; j++) {
      part[3+idimp*j] -= sum1;
      part[4+idimp*j] -= sum2;
      part[5+idimp*j] -= sum3;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cdblkp3l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int my, int mz, int mx1, int my1, int mxyz1,
              int *irc) {
/* this subroutine finds the maximum number of particles in each tile of
   mx, my, mz to calculate size of segmented particle array ppart
   linear interpolation
   input: all except kpic, nppmx, output: kpic, nppmx
   part = input particle array
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = position z of particle n
   kpic = output number of particles per tile
   nppmx = return maximum number of particles in tile
   idimp = size of phase space = 6
   nop = number of particles
   mx/my/mz = number of grids in sorting cell in x, y and z
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int j, k, n, m, l, mxy1, isum, ist, npx, ierr;
   ierr = 0;
   mxy1 = mx1*my1;
/* clear counter array */
   for (k = 0; k < mxyz1; k++) {
      kpic[k] = 0;
   }
/* find how many particles in each tile */
   for (j = 0; j < nop; j++) {
      n = part[idimp*j];
      n = n/mx;
      m = part[1+idimp*j];
      m = m/my;
      l = part[2+idimp*j];
      l = l/mz;
      m = n + mx1*m + mxy1*l;
      if (m < mxyz1) {
         kpic[m] += 1;
      }
      else {
         ierr = ierr > (m - mxyz1 + 1) ? ierr : (m - mxyz1 + 1);
      }
   }
/* find maximum */
   isum = 0;
   npx = 0;
   for (k = 0; k < mxyz1; k++) {
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
void cppmovin3lt(float part[], float ppart[], int kpic[], int nppmx,
                 int idimp, int nop, int mx, int my, int mz, int mx1,
                 int my1, int mxyz1, int *irc) {
/* this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
   and copies to segmented array ppart
   linear interpolation
   input: all except ppart, kpic, output: ppart, kpic
   part/ppart = input/output particle arrays
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = position z of particle n
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
   kpic = output number of particles per tile
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 6
   nop = number of particles
   mx/my/mz = number of grids in sorting cell in x, y and z
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int i, j, k, n, m, l, mxy1, ip, ierr;
   ierr = 0;
   mxy1 = mx1*my1;
/* clear counter array */
   for (k = 0; k < mxyz1; k++) {
      kpic[k] = 0;
   }
/* find addresses of particles at each tile and reorder particles */
   for (j = 0; j < nop; j++) {
      n = part[idimp*j];
      n = n/mx;
      m = part[1+idimp*j];
      m = m/my;
      l = part[2+idimp*j];
      l = l/mz;
      m = n + mx1*m + mxy1*l;
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
void cppmovin3ltp(float part[], float ppart[], int kpic[], int kp[],
                  int nppmx, int idimp, int nop, int mx, int my, int mz,
                  int mx1, int my1, int mxyz1, int *irc) {
/* this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
   and copies to segmented array ppart
   designed for NUMA architectures, where memory is associated with the
   processor which first writes a memory location.
   linear interpolation
   input: all except ppart, kpic, output: ppart, kpic
   part/ppart = input/output particle arrays
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = position z of particle n
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
   kpic = output number of particles per tile
   kp = original location of reordered particle
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 6
   nop = number of particles
   mx/my/mz = number of grids in sorting cell in x, y and z
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int i, j, k, n, m, l, mxy1, ip, npp, ierr;
   ierr = 0;
   mxy1 = mx1*my1;
/* clear counter array */
   for (k = 0; k < mxyz1; k++) {
      kpic[k] = 0;
   }
/* find addresses of particles at each tile to reorder particles */
   for (j = 0; j < nop; j++) {
      n = part[idimp*j];
      n = n/mx;
      m = part[1+idimp*j];
      m = m/my;
      l = part[2+idimp*j];
      l = l/mz;
      m = n + mx1*m + mxy1*l;
      ip = kpic[m];
      if (ip < nppmx) {
         for (i = 0; i < idimp; i++) {
            kp[ip+nppmx*m] = j;
         }
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
   for (k = 0; k < mxyz1; k++) {
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
void cppcheck3lt(float ppart[], int kpic[], int idimp, int nppmx,
                 int nx, int ny, int nz, int mx, int my, int mz,
                 int mx1, int my1, int mz1, int *irc) {
/* this subroutine performs a sanity check to make sure particles sorted
   by x,y,z grid in tiles of mx, my, mz, are all within bounds.
   tiles are assumed to be arranged in 3D linear memory, and transposed
   input: all except irc
   output: irc
   ppart[l][0][n] = position x of particle n in tile l
   ppart[l][1][n] = position y of particle n in tile l
   ppart[l][2][n] = position z of particle n in tile l
   kpic[l] = number of reordered output particles in tile l
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = number of grids in sorting cell in x/y/z
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mz1 = (system length in z direction - 1)/mz + 1
   irc = particle error, returned only if error occurs, when irc > 0
local data                                                            */
   int mxy1, mxyz1, noff, moff, loff, npp, j, k, l, nn, mm, ll, ist;
   float edgelx, edgely, edgelz, edgerx, edgery, edgerz, dx, dy, dz;
   mxy1 = mx1*my1;
   mxyz1 = mxy1*mz1;
/* loop over tiles */
#pragma omp parallel for \
private(j,k,l,noff,moff,loff,npp,nn,mm,ll,ist,edgelx,edgely,edgelz, \
edgerx,edgery,edgerz,dx,dy,dz)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      nn = nx - noff;
      nn = mx < nn ? mx : nn;
      mm = ny - moff;
      mm = my < mm ? my : mm;
      ll = nz - loff;
      ll = mz < ll ? mz : ll;
      edgelx = noff;
      edgerx = noff + nn;
      edgely = moff;
      edgery = moff + mm;
      edgelz = loff;
      edgerz = loff + ll;
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
         dx = ppart[j+nppmx*(idimp*l)];
         dy = ppart[j+nppmx*(1+idimp*l)];
         dz = ppart[j+nppmx*(2+idimp*l)];
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
         if (dz < edgelz)
            ist += 9;
         if (dz >= edgerz)
            ist += 18;
         if (ist > 0)
            *irc = l + 1;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cgppush3lt(float ppart[], float fxyz[], int kpic[], float qbm,
                float dt, float *ek, int idimp, int nppmx, int nx,
                int ny, int nz, int mx, int my, int mz, int nxv, int nyv,
                int nzv, int mx1, int my1, int mxyz1, int ipbc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions.
   OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   94 flops/particle, 30 loads, 6 stores
   input: all, output: part, ek
   equations used are:
   vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
   vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
   vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
   where q/m is charge/mass, and
   x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
   z(t+dt) = z(t) + vz(t+dt/2)*dt
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
                  + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
                  + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
   fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
                  + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
                  + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
   (vz(t+dt/2)+vz(t-dt/2))**2)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of field array, must be >= nx+1
   nyv = third dimension of field array, must be >= ny+1
   nzv = fourth dimension of field array, must be >= nz+1
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
#define N                4
   int mxy1, noff, moff, loff, npoff, npp;
   int i, j, k, l, nn, mm, ll, mxv, myv, mxyv, nxyv;
   float qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz;
   float vx, vy, vz;
   float sfxyz[N*MXV*MYV*MZV];
/* float sfxyz[N*(mx+1)*(my+1)*(mz+1)]; */
   double sum1, sum2;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   mxy1 = mx1*my1;
   qtm = qbm*dt;
   sum2 = 0.0;
/* set boundary values */
   edgelx = 0.0f;
   edgely = 0.0f;
   edgelz = 0.0f;
   edgerx = (float) nx;
   edgery = (float) ny;
   edgerz = (float) nz;
   if (ipbc==2) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgelz = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
      edgerz = (float) (nz-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,x,y,z,dxp,dyp,dzp, \
amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,sum1,sfxyz) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      ll = (mz < nz-loff ? mz : nz-loff) + 1;
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
            for (i = 0; i < nn; i++) {
               sfxyz[N*(i+mxv*j+mxyv*k)]
               = fxyz[N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+N*(i+mxv*j+mxyv*k)]
               = fxyz[1+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+N*(i+mxv*j+mxyv*k)]
               = fxyz[2+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
      sum1 = 0.0;
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = N*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find acceleration */
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+N];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+N];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+N];
         dx = amz*(dx + dyp*sfxyz[nn+N*mxv] + dx1*sfxyz[nn+N*mxv+N]);
         dy = amz*(dy + dyp*sfxyz[nn+N*mxv+1] + dx1*sfxyz[nn+N*mxv+1+N]);
         dz = amz*(dz + dyp*sfxyz[nn+N*mxv+2] + dx1*sfxyz[nn+N*mxv+2+N]);
         mm = nn + N*mxyv;
         vx = amx*sfxyz[mm] + amy*sfxyz[mm+N];
         vy = amx*sfxyz[mm+1] + amy*sfxyz[mm+1+N];
         vz = amx*sfxyz[mm+2] + amy*sfxyz[mm+2+N];
         dx = dx + dzp*(vx + dyp*sfxyz[mm+N*mxv] + dx1*sfxyz[mm+N*mxv+N]);
         dy = dy + dzp*(vy + dyp*sfxyz[mm+N*mxv+1] + dx1*sfxyz[mm+N*mxv+1+N]);
         dz = dz + dzp*(vz + dyp*sfxyz[mm+N*mxv+2] + dx1*sfxyz[mm+N*mxv+2+N]);
/* new velocity */
         dxp = ppart[j+3*nppmx+npoff];
         dyp = ppart[j+4*nppmx+npoff];
         dzp = ppart[j+5*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
         vz = dzp + qtm*dz;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         dzp += vz;
         sum1 += dxp*dxp + dyp*dyp+ dzp*dzp;
/* new position */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
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
            if ((dz < edgelz) || (dz >= edgerz)) {
               dz = z;
               vz = -vz;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               vx = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               vy = -vy;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
/* set new velocity */
         ppart[j+3*nppmx+npoff] = vx;
         ppart[j+4*nppmx+npoff] = vy;
         ppart[j+5*nppmx+npoff] = vz;
      }
      sum2 += sum1;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum2;
   return;
#undef N
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void cgppushf3lt(float ppart[], float fxyz[], int kpic[], int ncl[],
                 int ihole[], float qbm, float dt, float *ek, int idimp,   
                 int nppmx, int nx, int ny, int nz, int mx, int my,
                 int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                 int mxyz1, int ntmax, int *irc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   94 flops/particle, 30 loads, 6 stores
   input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
   equations used are:
   vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
   vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
   vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
   where q/m is charge/mass, and
   x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
   z(t+dt) = z(t) + vz(t+dt/2)*dt
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
                  + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
                  + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
   fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
                  + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
                  + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   kpic[l] = number of particles in tile l
   ncl[l][i] = number of particles going to destination i, tile l
   ihole[l][:][0] = location of hole in array left by departing particle
   ihole[l][:][1] = direction destination of particle leaving hole
   all for tile l
   ihole[l][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
   (vz(t+dt/2)+vz(t-dt/2))**2)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of field array, must be >= nx+1
   nyv = third dimension of field array, must be >= ny+1
   nzv = fourth dimension of field array, must be >= nz+1
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
#define N                4
   int mxy1, noff, moff, loff, npoff, npp;
   int i, j, k, l, ih, nh, nn, mm, ll, mxv, myv, mxyv, nxyv;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float qtm, x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz;
   float vx, vy, vz;
   float sfxyz[N*MXV*MYV*MZV];
/* float sfxyz[N*(mx+1)*(my+1)*(mz+1)]; */
   double sum1, sum2;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   mxy1 = mx1*my1;
   qtm = qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
   sum2 = 0.0;
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,ih,nh,x,y,z,dxp,dyp, \
dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,edgelx,edgely,edgelz,edgerx, \
edgery,edgerz,sum1,sfxyz) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      nn = nx - noff;
      nn = mx < nn ? mx : nn;
      mm = ny - moff;
      mm = my < mm ? my : mm;
      ll = nz - loff;
      ll = mz < ll ? mz : ll;
      edgelx = noff;
      edgerx = noff + nn;
      edgely = moff;
      edgery = moff + mm;
      edgelz = loff;
      edgerz = loff + ll;
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
      ll += 1;
/* load local fields from global array */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
            for (i = 0; i < nn; i++) {
               sfxyz[N*(i+mxv*j+mxyv*k)]
               = fxyz[N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+N*(i+mxv*j+mxyv*k)]
               = fxyz[1+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+N*(i+mxv*j+mxyv*k)]
               = fxyz[2+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
/* clear counters */
      for (j = 0; j < 26; j++) {
         ncl[j+26*l] = 0;
      }
      sum1 = 0.0;
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = N*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find acceleration */
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+N];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+N];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+N];
         dx = amz*(dx + dyp*sfxyz[nn+N*mxv] + dx1*sfxyz[nn+N*mxv+N]);
         dy = amz*(dy + dyp*sfxyz[nn+N*mxv+1] + dx1*sfxyz[nn+N*mxv+1+N]);
         dz = amz*(dz + dyp*sfxyz[nn+N*mxv+2] + dx1*sfxyz[nn+N*mxv+2+N]);
         mm = nn + N*mxyv;
         vx = amx*sfxyz[mm] + amy*sfxyz[mm+N];
         vy = amx*sfxyz[mm+1] + amy*sfxyz[mm+1+N];
         vz = amx*sfxyz[mm+2] + amy*sfxyz[mm+2+N];
         dx = dx + dzp*(vx + dyp*sfxyz[mm+N*mxv] + dx1*sfxyz[mm+N*mxv+N]);
         dy = dy + dzp*(vy + dyp*sfxyz[mm+N*mxv+1] + dx1*sfxyz[mm+N*mxv+1+N]);
         dz = dz + dzp*(vz + dyp*sfxyz[mm+N*mxv+2] + dx1*sfxyz[mm+N*mxv+2+N]);
/* new velocity */
         dxp = ppart[j+3*nppmx+npoff];
         dyp = ppart[j+4*nppmx+npoff];
         dzp = ppart[j+5*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
         vz = dzp + qtm*dz;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         dzp += vz;
         sum1 += dxp*dxp + dyp*dyp+ dzp*dzp;
/* new position */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
/* find particles going out of bounds */
         mm = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                             */
         if (dx >= edgerx) {
            if (dx >= anx)
               dx = dx - anx;
            mm = 2;
         }
         else if (dx < edgelx) {
            if (dx < 0.0f) {
               dx += anx;
               if (dx < anx)
                  mm = 1;
               else
                  dx = 0.0f;
            }
            else {
               mm = 1;
            }
         }
         if (dy >= edgery) {
            if (dy >= any)
               dy = dy - any;
            mm += 6;
         }
         else if (dy < edgely) {
            if (dy < 0.0f) {
               dy += any;
               if (dy < any)
                  mm += 3;
               else
                  dy = 0.0f;
            }
            else {
               mm += 3;
            }
         }
         if (dz >= edgerz) {
            if (dz >= anz)
               dz = dz - anz;
            mm += 18;
         }
         else if (dz < edgelz) {
            if (dz < 0.0f) {
               dz += anz;
               if (dz < anz)
                  mm += 9;
               else
                  dz = 0.0f;
            }
            else {
               mm += 9;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
/* set new velocity */
         ppart[j+3*nppmx+npoff] = vx;
         ppart[j+4*nppmx+npoff] = vy;
         ppart[j+5*nppmx+npoff] = vz;
/* increment counters */
         if (mm > 0) {
            ncl[mm+26*l-1] += 1;
            ih += 1;
            if (ih <= ntmax) {
               ihole[2*(ih+(ntmax+1)*l)] = j + 1;
               ihole[1+2*(ih+(ntmax+1)*l)] = mm;
            }
            else {
               nh = 1;
            }
         }
      }
      sum2 += sum1;
/* set error and end of file flag */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*l] = ih;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum2;
   return;
#undef N
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void cvgppush3lt(float ppart[], float fxyz[], int kpic[], float qbm,
                 float dt, float *ek, int idimp, int nppmx, int nx,
                 int ny, int nz, int mx, int my, int mz, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1,
                 int ipbc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions.
   vectorizable/OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   94 flops/particle, 30 loads, 6 stores
   input: all, output: ppart, ek
   equations used are:
   vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
   vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
   vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
   where q/m is charge/mass, and
   x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
   z(t+dt) = z(t) + vz(t+dt/2)*dt
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
                  + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
                  + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
   fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
                  + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
                  + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
   (vz(t+dt/2)+vz(t-dt/2))**2)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of field array, must be >= nx+1
   nyv = third dimension of field array, must be >= ny+1
   nzv = fourth dimension of field array, must be >= nz+1
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
#define NPBLK           32
#define LVECT            8
#define N                4
   int mxy1, noff, moff, loff, npoff, npp, ipp, joff, nps;
   int i, j, k, l, m, nn, mm, ll, lxv, lyv, lxyv, nxyv;
   float qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz;
   float vx, vy, vz;
   float sfxyz[N*MXV*MYV*MZV];
/* float sfxyz[N*(mx+1)*(my+1)*(mz+1)]; */
/* scratch arrays */
   __attribute__((aligned(64))) int n[NPBLK];
   __attribute__((aligned(64))) float s[NPBLK*LVECT], t[NPBLK*3];
   double sum1, sum2;
   mxy1 = mx1*my1;
/* lxv = MXV; */
/* lyv = MYV; */
   lxv = mx + 1;
   lyv = my + 1;
   lxyv = lxv*lyv;
   nxyv = nxv*nyv;
   qtm = qbm*dt;
   sum2 = 0.0;
/* set boundary values */
   edgelx = 0.0f;
   edgely = 0.0f;
   edgelz = 0.0f;
   edgerx = (float) nx;
   edgery = (float) ny;
   edgerz = (float) nz;
   if (ipbc==2) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgelz = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
      edgerz = (float) (nz-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,ipp,joff,nps,nn,mm,ll,x,y,z, \
dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,sum1,sfxyz,n,s,t) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      ll = (mz < nz-loff ? mz : nz-loff) + 1;
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
#pragma ivdep
            for (i = 0; i < nn; i++) {
               sfxyz[N*(i+lxv*j+lxyv*k)]
               = fxyz[N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+N*(i+lxv*j+lxyv*k)]
               = fxyz[1+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+N*(i+lxv*j+lxyv*k)]
               = fxyz[2+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
      sum1 = 0.0;
/* loop over particles in tile */
      ipp = npp/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m;
/* inner loop over particles in block */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
/* find interpolation weights */
            x = ppart[j+joff+npoff];
            y = ppart[j+joff+nppmx+npoff];
            z = ppart[j+joff+2*nppmx+npoff];
            nn = x;
            mm = y;
            ll = z;
            dxp = x - (float) nn;
            dyp = y - (float) mm;
            dzp = z - (float) ll;
            n[j] = nn - noff + lxv*(mm - moff) + lxyv*(ll - loff);
            amx = 1.0f - dxp;
            amy = 1.0f - dyp;
            dx1 = dxp*dyp;
            dyp = amx*dyp;
            amx = amx*amy;
            amz = 1.0f - dzp;
            amy = dxp*amy;
            s[j] = amx*amz;
            s[j+NPBLK] = amy*amz;
            s[j+2*NPBLK] = dyp*amz;
            s[j+3*NPBLK] = dx1*amz;
            s[j+4*NPBLK] = amx*dzp;
            s[j+5*NPBLK] = amy*dzp;
            s[j+6*NPBLK] = dyp*dzp;
            s[j+7*NPBLK] = dx1*dzp;
            t[j] = x;
            t[j+NPBLK] = y;
            t[j+2*NPBLK] = z;
         }
/* find acceleration */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + lxv - 2;
            ll = nn + lxyv - 4;
            k = ll + lxv - 2;
            dx = 0.0f;
            dy = 0.0f;
            dz = 0.0f;
            for (i = 0; i < LVECT; i++) {
               if (i > 5) {
                  nn = k;
               }
               else if (i > 3) {
                  nn = ll;
               }
               else if (i > 1) {
                  nn = mm;
               }
               dx += sfxyz[N*(i+nn)]*s[j+NPBLK*i];
               dy += sfxyz[1+N*(i+nn)]*s[j+NPBLK*i];
               dz += sfxyz[2+N*(i+nn)]*s[j+NPBLK*i];
            }
            s[j] = dx;
            s[j+NPBLK] = dy;
            s[j+2*NPBLK] = dz;
         }
/* new velocity */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
            z = t[j+2*NPBLK];
            dxp = ppart[j+joff+3*nppmx+npoff];
            dyp = ppart[j+joff+4*nppmx+npoff];
            dzp = ppart[j+joff+5*nppmx+npoff];
            vx = dxp + qtm*s[j];
            vy = dyp + qtm*s[j+NPBLK];
            vz = dzp + qtm*s[j+2*NPBLK];
/* average kinetic energy */
            dxp += vx;
            dyp += vy;
            dzp += vz;
            sum1 += dxp*dxp + dyp*dyp + dzp*dzp;
/* new position */
            s[j] = x + vx*dt;
            s[j+NPBLK] = y + vy*dt;
            s[j+2*NPBLK] = z + vz*dt;
            s[j+3*NPBLK] = vx;
            s[j+4*NPBLK] = vy;
            s[j+5*NPBLK] = vz;
         }
/* check boundary conditions */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
            dx = s[j];
            dy = s[j+NPBLK];
            dz = s[j+2*NPBLK];
            vx = s[j+3*NPBLK];
            vy = s[j+4*NPBLK];
            vz = s[j+5*NPBLK];
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
               if ((dz < edgelz) || (dz >= edgerz)) {
                  dz = t[j+2*NPBLK];
                  vz = -vz;
               }
            }
/* mixed reflecting/periodic boundary conditions */
            else if (ipbc==3) {
               if ((dx < edgelx) || (dx >= edgerx)) {
                  dx = t[j];
                  vx = -vx;
               }
               if ((dy < edgely) || (dy >= edgery)) {
                  dy = t[j+NPBLK];
                  vy = -vy;
               }
            }
/* set new position */
            ppart[j+joff+npoff] = dx;
            ppart[j+joff+nppmx+npoff] = dy;
            ppart[j+joff+2*nppmx+npoff] = dz;
/* set new velocity */
            ppart[j+joff+3*nppmx+npoff] = vx;
            ppart[j+joff+4*nppmx+npoff] = vy;
            ppart[j+joff+5*nppmx+npoff] = vz;
         }
      }
      nps = NPBLK*ipp;
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = N*(nn - noff + lxv*(mm - moff) + lxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find acceleration */
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+N];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+N];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+N];
         dx = amz*(dx + dyp*sfxyz[nn+N*lxv] + dx1*sfxyz[nn+N*lxv+N]);
         dy = amz*(dy + dyp*sfxyz[nn+N*lxv+1] + dx1*sfxyz[nn+N*lxv+1+N]);
         dz = amz*(dz + dyp*sfxyz[nn+N*lxv+2] + dx1*sfxyz[nn+N*lxv+2+N]);
         mm = nn + N*lxyv;
         vx = amx*sfxyz[mm] + amy*sfxyz[mm+N];
         vy = amx*sfxyz[mm+1] + amy*sfxyz[mm+1+N];
         vz = amx*sfxyz[mm+2] + amy*sfxyz[mm+2+N];
         dx = dx + dzp*(vx + dyp*sfxyz[mm+N*lxv] + dx1*sfxyz[mm+N*lxv+N]);
         dy = dy + dzp*(vy + dyp*sfxyz[mm+N*lxv+1] + dx1*sfxyz[mm+N*lxv+1+N]);
         dz = dz + dzp*(vz + dyp*sfxyz[mm+N*lxv+2] + dx1*sfxyz[mm+N*lxv+2+N]);
/* new velocity */
         dxp = ppart[j+3*nppmx+npoff];
         dyp = ppart[j+4*nppmx+npoff];
         dzp = ppart[j+5*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
         vz = dzp + qtm*dz;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         dzp += vz;
         sum1 += dxp*dxp + dyp*dyp+ dzp*dzp;
/* new position */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
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
            if ((dz < edgelz) || (dz >= edgerz)) {
               dz = z;
               vz = -vz;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               vx = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               vy = -vy;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
/* set new velocity */
         ppart[j+3*nppmx+npoff] = vx;
         ppart[j+4*nppmx+npoff] = vy;
         ppart[j+5*nppmx+npoff] = vz;
      }
      sum2 += sum1;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum2;
#undef LVECT
#undef NPBLK
#undef N
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void cvgppushf3lt(float ppart[], float fxyz[], int kpic[], int ncl[],
                  int ihole[], float qbm, float dt, float *ek,
                  int idimp, int nppmx, int nx, int ny, int nz, int mx,
                  int my, int mz, int nxv, int nyv, int nzv, int mx1,
                  int my1, int mxyz1, int ntmax, int *irc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   vectorizable/OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   94 flops/particle, 30 loads, 6 stores
   input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
   equations used are:
   vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
   vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
   vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
   where q/m is charge/mass, and
   x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
   z(t+dt) = z(t) + vz(t+dt/2)*dt
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
                  + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
                  + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
   fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
                  + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
                  + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   kpic[l] = number of particles in tile l
   ncl[l][i] = number of particles going to destination i, tile l
   ihole[l][:][0] = location of hole in array left by departing particle
   ihole[l][:][1] = direction destination of particle leaving hole
   all for tile l
   ihole[l][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
   (vz(t+dt/2)+vz(t-dt/2))**2)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of field array, must be >= nx+1
   nyv = third dimension of field array, must be >= ny+1
   nzv = fourth dimension of field array, must be >= nz+1
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
#define NPBLK           32
#define LVECT            8
#define N                4
   int mxy1, noff, moff, loff, npp, npoff, ipp, joff, nps;
   int i, j, k, l, m, ih, nh, nn, mm, ll, lxv, lyv, lxyv, nxyv;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float qtm, x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz;
   float vx, vy, vz;
   float sfxyz[N*MXV*MYV*MZV];
/* float sfxyz[N*(mx+1)*(my+1)*(mz+1)]; */
/* scratch arrays */
   __attribute__((aligned(64))) int n[NPBLK];
   __attribute__((aligned(64))) float s[NPBLK*LVECT], t[NPBLK*3];
   double sum1, sum2;
   mxy1 = mx1*my1;
/* lxv = MXV; */
/* lyv = MYV; */
   lxv = mx + 1;
   lyv = my + 1;
   lxyv = lxv*lyv;
   nxyv = nxv*nyv;
   qtm = qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
   sum2 = 0.0;
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,ipp,joff,nps,nn,mm,ll,ih,nh, \
x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,edgelx,edgely, \
edgelz,edgerx,edgery,edgerz,sum1,sfxyz,n,s,t) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      nn = nx - noff;
      nn = mx < nn ? mx : nn;
      mm = ny - moff;
      mm = my < mm ? my : mm;
      ll = nz - loff;
      ll = mz < ll ? mz : ll;
      edgelx = noff;
      edgerx = noff + nn;
      edgely = moff;
      edgery = moff + mm;
      edgelz = loff;
      edgerz = loff + ll;
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
      ll += 1;
/* load local fields from global array */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
#pragma ivdep
            for (i = 0; i < nn; i++) {
               sfxyz[N*(i+lxv*j+lxyv*k)]
               = fxyz[N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+N*(i+lxv*j+lxyv*k)]
               = fxyz[1+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+N*(i+lxv*j+lxyv*k)]
               = fxyz[2+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
/* clear counters */
      for (j = 0; j < 26; j++) {
         ncl[j+26*l] = 0;
      }
      sum1 = 0.0;
/* loop over particles in tile */
      ipp = npp/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m;
/* inner loop over particles in block */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
/* find interpolation weights */
            x = ppart[j+joff+npoff];
            y = ppart[j+joff+nppmx+npoff];
            z = ppart[j+joff+2*nppmx+npoff];
            nn = x;
            mm = y;
            ll = z;
            dxp = x - (float) nn;
            dyp = y - (float) mm;
            dzp = z - (float) ll;
            n[j] = nn - noff + lxv*(mm - moff) + lxyv*(ll - loff);
            amx = 1.0f - dxp;
            amy = 1.0f - dyp;
            dx1 = dxp*dyp;
            dyp = amx*dyp;
            amx = amx*amy;
            amz = 1.0f - dzp;
            amy = dxp*amy;
            s[j] = amx*amz;
            s[j+NPBLK] = amy*amz;
            s[j+2*NPBLK] = dyp*amz;
            s[j+3*NPBLK] = dx1*amz;
            s[j+4*NPBLK] = amx*dzp;
            s[j+5*NPBLK] = amy*dzp;
            s[j+6*NPBLK] = dyp*dzp;
            s[j+7*NPBLK] = dx1*dzp;
            t[j] = x;
            t[j+NPBLK] = y;
            t[j+2*NPBLK] = z;
         }
/* find acceleration */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + lxv - 2;
            ll = nn + lxyv - 4;
            k = ll + lxv - 2;
            dx = 0.0f;
            dy = 0.0f;
            dz = 0.0f;
            for (i = 0; i < LVECT; i++) {
               if (i > 5) {
                  nn = k;
               }
               else if (i > 3) {
                  nn = ll;
               }
               else if (i > 1) {
                  nn = mm;
               }
               dx += sfxyz[N*(i+nn)]*s[j+NPBLK*i];
               dy += sfxyz[1+N*(i+nn)]*s[j+NPBLK*i];
               dz += sfxyz[2+N*(i+nn)]*s[j+NPBLK*i];
            }
            s[j] = dx;
            s[j+NPBLK] = dy;
            s[j+2*NPBLK] = dz;
         }
/* new velocity */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
            z = t[j+2*NPBLK];
            dxp = ppart[j+joff+3*nppmx+npoff];
            dyp = ppart[j+joff+4*nppmx+npoff];
            dzp = ppart[j+joff+5*nppmx+npoff];
            vx = dxp + qtm*s[j];
            vy = dyp + qtm*s[j+NPBLK];
            vz = dzp + qtm*s[j+2*NPBLK];
/* average kinetic energy */
            dxp += vx;
            dyp += vy;
            dzp += vz;
            sum1 += dxp*dxp + dyp*dyp + dzp*dzp;
/* new position */
            s[j] = x + vx*dt;
            s[j+NPBLK] = y + vy*dt;
            s[j+2*NPBLK] = z + vz*dt;
            s[j+3*NPBLK] = vx;
            s[j+4*NPBLK] = vy;
            s[j+5*NPBLK] = vz;
         }
/* check boundary conditions */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
            dx = s[j];
            dy = s[j+NPBLK];
            dz = s[j+2*NPBLK];
/* find particles going out of bounds */
            mm = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                             */
            if (dx >= edgerx) {
               if (dx >= anx)
                  dx = dx - anx;
               mm = 2;
            }
            else if (dx < edgelx) {
               if (dx < 0.0f) {
                  dx += anx;
                  if (dx < anx)
                     mm = 1;
                  else
                     dx = 0.0f;
               }
               else {
                  mm = 1;
               }
            }
            if (dy >= edgery) {
               if (dy >= any)
                  dy = dy - any;
               mm += 6;
            }
            else if (dy < edgely) {
               if (dy < 0.0f) {
                  dy += any;
                  if (dy < any)
                     mm += 3;
                  else
                     dy = 0.0f;
               }
               else {
                  mm += 3;
               }
            }
            if (dz >= edgerz) {
               if (dz >= anz)
                  dz = dz - anz;
               mm += 18;
            }
            else if (dz < edgelz) {
               if (dz < 0.0f) {
                  dz += anz;
                  if (dz < anz)
                     mm += 9;
                  else
                     dz = 0.0f;
               }
               else {
                  mm += 9;
               }
            }
/* set new position */
            ppart[j+joff+npoff] = dx;
            ppart[j+joff+nppmx+npoff] = dy;
            ppart[j+joff+2*nppmx+npoff] = dz;
/* set new velocity */
            ppart[j+joff+3*nppmx+npoff] = s[j+3*NPBLK];
            ppart[j+joff+4*nppmx+npoff] = s[j+4*NPBLK];
            ppart[j+joff+5*nppmx+npoff] = s[j+5*NPBLK];
            n[j] = mm;
         }
         for (j = 0; j < NPBLK; j++) {
            mm = n[j];
/* increment counters */
            if (mm > 0) {
               ncl[mm+26*l-1] += 1;
               ih += 1;
               if (ih <= ntmax) {
                  ihole[2*(ih+(ntmax+1)*l)] = j + joff + 1;
                  ihole[1+2*(ih+(ntmax+1)*l)] = mm;
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
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = N*(nn - noff + lxv*(mm - moff) + lxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find acceleration */
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+N];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+N];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+N];
         dx = amz*(dx + dyp*sfxyz[nn+N*lxv] + dx1*sfxyz[nn+N*lxv+N]);
         dy = amz*(dy + dyp*sfxyz[nn+N*lxv+1] + dx1*sfxyz[nn+N*lxv+1+N]);
         dz = amz*(dz + dyp*sfxyz[nn+N*lxv+2] + dx1*sfxyz[nn+N*lxv+2+N]);
         mm = nn + N*lxyv;
         vx = amx*sfxyz[mm] + amy*sfxyz[mm+N];
         vy = amx*sfxyz[mm+1] + amy*sfxyz[mm+1+N];
         vz = amx*sfxyz[mm+2] + amy*sfxyz[mm+2+N];
         dx = dx + dzp*(vx + dyp*sfxyz[mm+N*lxv] + dx1*sfxyz[mm+N*lxv+N]);
         dy = dy + dzp*(vy + dyp*sfxyz[mm+N*lxv+1] + dx1*sfxyz[mm+N*lxv+1+N]);
         dz = dz + dzp*(vz + dyp*sfxyz[mm+N*lxv+2] + dx1*sfxyz[mm+N*lxv+2+N]);
/* new velocity */
         dxp = ppart[j+3*nppmx+npoff];
         dyp = ppart[j+4*nppmx+npoff];
         dzp = ppart[j+5*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
         vz = dzp + qtm*dz;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         dzp += vz;
         sum1 += dxp*dxp + dyp*dyp+ dzp*dzp;
/* new position */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
/* find particles going out of bounds */
         mm = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                             */
         if (dx >= edgerx) {
            if (dx >= anx)
               dx = dx - anx;
            mm = 2;
         }
         else if (dx < edgelx) {
            if (dx < 0.0f) {
               dx += anx;
               if (dx < anx)
                  mm = 1;
               else
                  dx = 0.0f;
            }
            else {
               mm = 1;
            }
         }
         if (dy >= edgery) {
            if (dy >= any)
               dy = dy - any;
            mm += 6;
         }
         else if (dy < edgely) {
            if (dy < 0.0f) {
               dy += any;
               if (dy < any)
                  mm += 3;
               else
                  dy = 0.0f;
            }
            else {
               mm += 3;
            }
         }
         if (dz >= edgerz) {
            if (dz >= anz)
               dz = dz - anz;
            mm += 18;
         }
         else if (dz < edgelz) {
            if (dz < 0.0f) {
               dz += anz;
               if (dz < anz)
                  mm += 9;
               else
                  dz = 0.0f;
            }
            else {
               mm += 9;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
/* set new velocity */
         ppart[j+3*nppmx+npoff] = vx;
         ppart[j+4*nppmx+npoff] = vy;
         ppart[j+5*nppmx+npoff] = vz;
/* increment counters */
         if (mm > 0) {
            ncl[mm+26*l-1] += 1;
            ih += 1;
            if (ih <= ntmax) {
               ihole[2*(ih+(ntmax+1)*l)] = j + 1;
               ihole[1+2*(ih+(ntmax+1)*l)] = mm;
            }
            else {
               nh = 1;
            }
         }
      }
      sum2 += sum1;
/* set error and end of file flag */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*l] = ih;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum2;
#undef LVECT
#undef NPBLK
#undef N
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void cv2gppush3lt(float ppart[], float fxyz[], int kpic[], float qbm,
                  float dt, float *ek, int idimp, int nppmx, int nx,
                  int ny, int nz, int mx, int my, int mz, int nxv,
                  int nyv, int nzv, int mx1, int my1, int mxyz1,
                  int ipbc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions.
   vectorizable/OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   94 flops/particle, 30 loads, 6 stores
   input: all, output: part, ek
   equations used are:
   vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
   vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
   vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
   where q/m is charge/mass, and
   x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
   z(t+dt) = z(t) + vz(t+dt/2)*dt
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
                  + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
                  + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
   fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
                  + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
                  + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
   (vz(t+dt/2)+vz(t-dt/2))**2)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of field array, must be >= nx+1
   nyv = third dimension of field array, must be >= ny+1
   nzv = fourth dimension of field array, must be >= nz+1
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
#define NPBLK           32
#define LVECT            8
#define N                4
   int mxy1, noff, moff, loff, npoff, npp, ipp, joff, nps;
   int i, j, k, l, m, nn, mm, ll, lxv, lyv, lxyv, nxyv;
   float qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz;
   float vx, vy, vz;
   float sfxyz[N*MXV*MYV*MZV];
/* float sfxyz[N*(mx+1)*(my+1)*(mz+1)]; */
/* scratch arrays */
   __attribute__((aligned(64))) int n[NPBLK], mn[LVECT];;
   __attribute__((aligned(64))) float s[NPBLK*LVECT], t[NPBLK*3];
   double sum1, sum2;
   mxy1 = mx1*my1;
/* lxv = MXV; */
/* lyv = MYV; */
   lxv = mx + 1;
   lyv = my + 1;
   lxyv = lxv*lyv;
   nxyv = nxv*nyv;
   mn[0] = 0;
   mn[1] = N;
   mn[2] = N*lxv;
   mn[3] = N*(lxv + 1);
   mn[4] = N*lxyv;
   mn[5] = N*(lxyv + 1);
   mn[6] = N*(lxyv + lxv);
   mn[7] = N*(lxyv + lxv + 1);
   qtm = qbm*dt;
   sum2 = 0.0;
/* set boundary values */
   edgelx = 0.0f;
   edgely = 0.0f;
   edgelz = 0.0f;
   edgerx = (float) nx;
   edgery = (float) ny;
   edgerz = (float) nz;
   if (ipbc==2) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgelz = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
      edgerz = (float) (nz-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,ipp,joff,nps,nn,mm,ll,x,y,z, \
dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,sum1,sfxyz,n,s,t) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      ll = (mz < nz-loff ? mz : nz-loff) + 1;
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
#pragma ivdep
            for (i = 0; i < nn; i++) {
               sfxyz[N*(i+lxv*j+lxyv*k)]
               = fxyz[N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+N*(i+lxv*j+lxyv*k)]
               = fxyz[1+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+N*(i+lxv*j+lxyv*k)]
               = fxyz[2+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
      sum1 = 0.0;
/* loop over particles in tile */
      ipp = npp/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m;
/* inner loop over particles in block */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
/* find interpolation weights */
            x = ppart[j+joff+npoff];
            y = ppart[j+joff+nppmx+npoff];
            z = ppart[j+joff+2*nppmx+npoff];
            nn = x;
            mm = y;
            ll = z;
            dxp = x - (float) nn;
            dyp = y - (float) mm;
            dzp = z - (float) ll;
            n[j] = N*(nn - noff + lxv*(mm - moff) + lxyv*(ll - loff));
            amx = 1.0f - dxp;
            amy = 1.0f - dyp;
            dx1 = dxp*dyp;
            dyp = amx*dyp;
            amx = amx*amy;
            amz = 1.0f - dzp;
            amy = dxp*amy;
            s[j] = amx*amz;
            s[j+NPBLK] = amy*amz;
            s[j+2*NPBLK] = dyp*amz;
            s[j+3*NPBLK] = dx1*amz;
            s[j+4*NPBLK] = amx*dzp;
            s[j+5*NPBLK] = amy*dzp;
            s[j+6*NPBLK] = dyp*dzp;
            s[j+7*NPBLK] = dx1*dzp;
            t[j] = x;
            t[j+NPBLK] = y;
            t[j+2*NPBLK] = z;
         }
/* find acceleration */
         for (j = 0; j < NPBLK; j++) {
            dx = 0.0f;
            dy = 0.0f;
            dz = 0.0f;
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               dx += sfxyz[n[j]+mn[i]]*s[j+NPBLK*i];
               dy += sfxyz[1+n[j]+mn[i]]*s[j+NPBLK*i];
               dz += sfxyz[2+n[j]+mn[i]]*s[j+NPBLK*i];
            }
            s[j] = dx;
            s[j+NPBLK] = dy;
            s[j+2*NPBLK] = dz;
         }
/* new velocity */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
            z = t[j+2*NPBLK];
            dxp = ppart[j+joff+3*nppmx+npoff];
            dyp = ppart[j+joff+4*nppmx+npoff];
            dzp = ppart[j+joff+5*nppmx+npoff];
            vx = dxp + qtm*s[j];
            vy = dyp + qtm*s[j+NPBLK];
            vz = dzp + qtm*s[j+2*NPBLK];
/* average kinetic energy */
            dxp += vx;
            dyp += vy;
            dzp += vz;
            sum1 += dxp*dxp + dyp*dyp + dzp*dzp;
/* new position */
            s[j] = x + vx*dt;
            s[j+NPBLK] = y + vy*dt;
            s[j+2*NPBLK] = z + vz*dt;
            s[j+3*NPBLK] = vx;
            s[j+4*NPBLK] = vy;
            s[j+5*NPBLK] = vz;
         }
/* check boundary conditions */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
            dx = s[j];
            dy = s[j+NPBLK];
            dz = s[j+2*NPBLK];
            vx = s[j+3*NPBLK];
            vy = s[j+4*NPBLK];
            vz = s[j+5*NPBLK];
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
               if ((dz < edgelz) || (dz >= edgerz)) {
                  dz = t[j+2*NPBLK];
                  vz = -vz;
               }
            }
/* mixed reflecting/periodic boundary conditions */
            else if (ipbc==3) {
               if ((dx < edgelx) || (dx >= edgerx)) {
                  dx = t[j];
                  vx = -vx;
               }
               if ((dy < edgely) || (dy >= edgery)) {
                  dy = t[j+NPBLK];
                  vy = -vy;
               }
            }
/* set new position */
            ppart[j+joff+npoff] = dx;
            ppart[j+joff+nppmx+npoff] = dy;
            ppart[j+joff+2*nppmx+npoff] = dz;
/* set new velocity */
            ppart[j+joff+3*nppmx+npoff] = vx;
            ppart[j+joff+4*nppmx+npoff] = vy;
            ppart[j+joff+5*nppmx+npoff] = vz;
         }
      }
      nps = NPBLK*ipp;
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = N*(nn - noff + lxv*(mm - moff) + lxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find acceleration */
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+N];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+N];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+N];
         dx = amz*(dx + dyp*sfxyz[nn+N*lxv] + dx1*sfxyz[nn+N*lxv+N]);
         dy = amz*(dy + dyp*sfxyz[nn+N*lxv+1] + dx1*sfxyz[nn+N*lxv+1+N]);
         dz = amz*(dz + dyp*sfxyz[nn+N*lxv+2] + dx1*sfxyz[nn+N*lxv+2+N]);
         mm = nn + N*lxyv;
         vx = amx*sfxyz[mm] + amy*sfxyz[mm+N];
         vy = amx*sfxyz[mm+1] + amy*sfxyz[mm+1+N];
         vz = amx*sfxyz[mm+2] + amy*sfxyz[mm+2+N];
         dx = dx + dzp*(vx + dyp*sfxyz[mm+N*lxv] + dx1*sfxyz[mm+N*lxv+N]);
         dy = dy + dzp*(vy + dyp*sfxyz[mm+N*lxv+1] + dx1*sfxyz[mm+N*lxv+1+N]);
         dz = dz + dzp*(vz + dyp*sfxyz[mm+N*lxv+2] + dx1*sfxyz[mm+N*lxv+2+N]);
/* new velocity */
         dxp = ppart[j+3*nppmx+npoff];
         dyp = ppart[j+4*nppmx+npoff];
         dzp = ppart[j+5*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
         vz = dzp + qtm*dz;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         dzp += vz;
         sum1 += dxp*dxp + dyp*dyp+ dzp*dzp;
/* new position */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
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
            if ((dz < edgelz) || (dz >= edgerz)) {
               dz = z;
               vz = -vz;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               vx = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               vy = -vy;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
/* set new velocity */
         ppart[j+3*nppmx+npoff] = vx;
         ppart[j+4*nppmx+npoff] = vy;
         ppart[j+5*nppmx+npoff] = vz;
      }
      sum2 += sum1;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum2;
#undef LVECT
#undef NPBLK
#undef N
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void cv2gppushf3lt(float ppart[], float fxyz[], int kpic[], int ncl[],
                   int ihole[], float qbm, float dt, float *ek,
                   int idimp, int nppmx, int nx, int ny, int nz, int mx,
                   int my, int mz, int nxv, int nyv, int nzv, int mx1,
                   int my1, int mxyz1, int ntmax, int *irc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   vectorizable/OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   94 flops/particle, 30 loads, 6 stores
   input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
   equations used are:
   vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
   vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
   vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
   where q/m is charge/mass, and
   x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
   z(t+dt) = z(t) + vz(t+dt/2)*dt
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
                  + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
                  + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
   fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
                  + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
                  + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   kpic[l] = number of particles in tile l
   ncl[l][i] = number of particles going to destination i, tile l
   ihole[l][:][0] = location of hole in array left by departing particle
   ihole[l][:][1] = direction destination of particle leaving hole
   all for tile l
   ihole[l][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
   (vz(t+dt/2)+vz(t-dt/2))**2)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of field array, must be >= nx+1
   nyv = third dimension of field array, must be >= ny+1
   nzv = fourth dimension of field array, must be >= nz+1
   ipbc = particle boundary condition = (0,1,2,3) =
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
#define NPBLK           32
#define LVECT            8
#define N                4
   int mxy1, noff, moff, loff, npp, npoff, ipp, joff, nps;
   int i, j, k, l, m, ih, nh, nn, mm, ll, lxv, lyv, lxyv, nxyv;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float qtm, x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz;
   float vx, vy, vz;
   float sfxyz[N*MXV*MYV*MZV];
/* float sfxyz[N*(mx+1)*(my+1)*(mz+1)]; */
/* scratch arrays */
   __attribute__((aligned(64))) int n[NPBLK], mn[LVECT];;
   __attribute__((aligned(64))) float s[NPBLK*LVECT], t[NPBLK*3];
   double sum1, sum2;
   mxy1 = mx1*my1;
/* lxv = MXV; */
/* lyv = MYV; */
   lxv = mx + 1;
   lyv = my + 1;
   lxyv = lxv*lyv;
   nxyv = nxv*nyv;
   mn[0] = 0;
   mn[1] = N;
   mn[2] = N*lxv;
   mn[3] = N*(lxv + 1);
   mn[4] = N*lxyv;
   mn[5] = N*(lxyv + 1);
   mn[6] = N*(lxyv + lxv);
   mn[7] = N*(lxyv + lxv + 1);
   qtm = qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
   sum2 = 0.0;
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,ipp,joff,nps,nn,mm,ll,ih,nh, \
x,y,z,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,edgelx,edgely, \
edgelz,edgerx,edgery,edgerz,sum1,sfxyz,n,s,t) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      nn = nx - noff;
      nn = mx < nn ? mx : nn;
      mm = ny - moff;
      mm = my < mm ? my : mm;
      ll = nz - loff;
      ll = mz < ll ? mz : ll;
      edgelx = noff;
      edgerx = noff + nn;
      edgely = moff;
      edgery = moff + mm;
      edgelz = loff;
      edgerz = loff + ll;
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
      ll += 1;
/* load local fields from global array */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
#pragma ivdep
            for (i = 0; i < nn; i++) {
               sfxyz[N*(i+lxv*j+lxyv*k)]
               = fxyz[N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+N*(i+lxv*j+lxyv*k)]
               = fxyz[1+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+N*(i+lxv*j+lxyv*k)]
               = fxyz[2+N*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
/* clear counters */
      for (j = 0; j < 26; j++) {
         ncl[j+26*l] = 0;
      }
      sum1 = 0.0;
/* loop over particles in tile */
      ipp = npp/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m;
/* inner loop over particles in block */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
/* find interpolation weights */
            x = ppart[j+joff+npoff];
            y = ppart[j+joff+nppmx+npoff];
            z = ppart[j+joff+2*nppmx+npoff];
            nn = x;
            mm = y;
            ll = z;
            dxp = x - (float) nn;
            dyp = y - (float) mm;
            dzp = z - (float) ll;
            n[j] = N*(nn - noff + lxv*(mm - moff) + lxyv*(ll - loff));
            amx = 1.0f - dxp;
            amy = 1.0f - dyp;
            dx1 = dxp*dyp;
            dyp = amx*dyp;
            amx = amx*amy;
            amz = 1.0f - dzp;
            amy = dxp*amy;
            s[j] = amx*amz;
            s[j+NPBLK] = amy*amz;
            s[j+2*NPBLK] = dyp*amz;
            s[j+3*NPBLK] = dx1*amz;
            s[j+4*NPBLK] = amx*dzp;
            s[j+5*NPBLK] = amy*dzp;
            s[j+6*NPBLK] = dyp*dzp;
            s[j+7*NPBLK] = dx1*dzp;
            t[j] = x;
            t[j+NPBLK] = y;
            t[j+2*NPBLK] = z;
         }
/* find acceleration */
         for (j = 0; j < NPBLK; j++) {
            dx = 0.0f;
            dy = 0.0f;
            dz = 0.0f;
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               dx += sfxyz[n[j]+mn[i]]*s[j+NPBLK*i];
               dy += sfxyz[1+n[j]+mn[i]]*s[j+NPBLK*i];
               dz += sfxyz[2+n[j]+mn[i]]*s[j+NPBLK*i];
            }
            s[j] = dx;
            s[j+NPBLK] = dy;
            s[j+2*NPBLK] = dz;
         }
/* new velocity */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
            z = t[j+2*NPBLK];
            dxp = ppart[j+joff+3*nppmx+npoff];
            dyp = ppart[j+joff+4*nppmx+npoff];
            dzp = ppart[j+joff+5*nppmx+npoff];
            vx = dxp + qtm*s[j];
            vy = dyp + qtm*s[j+NPBLK];
            vz = dzp + qtm*s[j+2*NPBLK];
/* average kinetic energy */
            dxp += vx;
            dyp += vy;
            dzp += vz;
            sum1 += dxp*dxp + dyp*dyp + dzp*dzp;
/* new position */
            s[j] = x + vx*dt;
            s[j+NPBLK] = y + vy*dt;
            s[j+2*NPBLK] = z + vz*dt;
            s[j+3*NPBLK] = vx;
            s[j+4*NPBLK] = vy;
            s[j+5*NPBLK] = vz;
         }
/* check boundary conditions */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
            dx = s[j];
            dy = s[j+NPBLK];
            dz = s[j+2*NPBLK];
/* find particles going out of bounds */
            mm = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                             */
            if (dx >= edgerx) {
               if (dx >= anx)
                  dx = dx - anx;
               mm = 2;
            }
            else if (dx < edgelx) {
               if (dx < 0.0f) {
                  dx += anx;
                  if (dx < anx)
                     mm = 1;
                  else
                     dx = 0.0f;
               }
               else {
                  mm = 1;
               }
            }
            if (dy >= edgery) {
               if (dy >= any)
                  dy = dy - any;
               mm += 6;
            }
            else if (dy < edgely) {
               if (dy < 0.0f) {
                  dy += any;
                  if (dy < any)
                     mm += 3;
                  else
                     dy = 0.0f;
               }
               else {
                  mm += 3;
               }
            }
            if (dz >= edgerz) {
               if (dz >= anz)
                  dz = dz - anz;
               mm += 18;
            }
            else if (dz < edgelz) {
               if (dz < 0.0f) {
                  dz += anz;
                  if (dz < anz)
                     mm += 9;
                  else
                     dz = 0.0f;
               }
               else {
                  mm += 9;
               }
            }
/* set new position */
            ppart[j+joff+npoff] = dx;
            ppart[j+joff+nppmx+npoff] = dy;
            ppart[j+joff+2*nppmx+npoff] = dz;
/* set new velocity */
            ppart[j+joff+3*nppmx+npoff] = s[j+3*NPBLK];
            ppart[j+joff+4*nppmx+npoff] = s[j+4*NPBLK];
            ppart[j+joff+5*nppmx+npoff] = s[j+5*NPBLK];
            n[j] = mm;
         }
         for (j = 0; j < NPBLK; j++) {
            mm = n[j];
/* increment counters */
            if (mm > 0) {
               ncl[mm+26*l-1] += 1;
               ih += 1;
               if (ih <= ntmax) {
                  ihole[2*(ih+(ntmax+1)*l)] = j + joff + 1;
                  ihole[1+2*(ih+(ntmax+1)*l)] = mm;
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
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = N*(nn - noff + lxv*(mm - moff) + lxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find acceleration */
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+N];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+N];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+N];
         dx = amz*(dx + dyp*sfxyz[nn+N*lxv] + dx1*sfxyz[nn+N*lxv+N]);
         dy = amz*(dy + dyp*sfxyz[nn+N*lxv+1] + dx1*sfxyz[nn+N*lxv+1+N]);
         dz = amz*(dz + dyp*sfxyz[nn+N*lxv+2] + dx1*sfxyz[nn+N*lxv+2+N]);
         mm = nn + N*lxyv;
         vx = amx*sfxyz[mm] + amy*sfxyz[mm+N];
         vy = amx*sfxyz[mm+1] + amy*sfxyz[mm+1+N];
         vz = amx*sfxyz[mm+2] + amy*sfxyz[mm+2+N];
         dx = dx + dzp*(vx + dyp*sfxyz[mm+N*lxv] + dx1*sfxyz[mm+N*lxv+N]);
         dy = dy + dzp*(vy + dyp*sfxyz[mm+N*lxv+1] + dx1*sfxyz[mm+N*lxv+1+N]);
         dz = dz + dzp*(vz + dyp*sfxyz[mm+N*lxv+2] + dx1*sfxyz[mm+N*lxv+2+N]);
/* new velocity */
         dxp = ppart[j+3*nppmx+npoff];
         dyp = ppart[j+4*nppmx+npoff];
         dzp = ppart[j+5*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
         vz = dzp + qtm*dz;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         dzp += vz;
         sum1 += dxp*dxp + dyp*dyp+ dzp*dzp;
/* new position */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
/* find particles going out of bounds */
         mm = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                             */
         if (dx >= edgerx) {
            if (dx >= anx)
               dx = dx - anx;
            mm = 2;
         }
         else if (dx < edgelx) {
            if (dx < 0.0f) {
               dx += anx;
               if (dx < anx)
                  mm = 1;
               else
                  dx = 0.0f;
            }
            else {
               mm = 1;
            }
         }
         if (dy >= edgery) {
            if (dy >= any)
               dy = dy - any;
            mm += 6;
         }
         else if (dy < edgely) {
            if (dy < 0.0f) {
               dy += any;
               if (dy < any)
                  mm += 3;
               else
                  dy = 0.0f;
            }
            else {
               mm += 3;
            }
         }
         if (dz >= edgerz) {
            if (dz >= anz)
               dz = dz - anz;
            mm += 18;
         }
         else if (dz < edgelz) {
            if (dz < 0.0f) {
               dz += anz;
               if (dz < anz)
                  mm += 9;
               else
                  dz = 0.0f;
            }
            else {
               mm += 9;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
/* set new velocity */
         ppart[j+3*nppmx+npoff] = vx;
         ppart[j+4*nppmx+npoff] = vy;
         ppart[j+5*nppmx+npoff] = vz;
/* increment counters */
         if (mm > 0) {
            ncl[mm+26*l-1] += 1;
            ih += 1;
            if (ih <= ntmax) {
               ihole[2*(ih+(ntmax+1)*l)] = j + 1;
               ihole[1+2*(ih+(ntmax+1)*l)] = mm;
            }
            else {
               nh = 1;
            }
         }
      }
      sum2 += sum1;
/* set error and end of file flag */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*l] = ih;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum2;
   return;
#undef LVECT
#undef NPBLK
#undef N
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void cgppost3lt(float ppart[], float q[], int kpic[], float qm,
                int nppmx, int idimp, int mx, int my, int mz, int nxv,
                int nyv, int nzv, int mx1, int my1, int mxyz1) {
/* for 3d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   33 flops/particle, 11 loads, 8 stores
   input: all, output: q
   charge density is approximated by values at the nearest grid points
   q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
   q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
   q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
   q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
   q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
   q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
   q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
   q(n+1,m+1,l+1)=qm*dx*dy*dz
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   q[l][k][j] = charge density at grid point j,k,l
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 6
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = first dimension of charge array, must be >= nx+1
   nyv = second dimension of charge array, must be >= ny+1
   nzv = third dimension of charge array, must be >= nz+1
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp;
   int i, j, k, l, nn, mm, ll, nm, lm, mxv, myv, mxyv, nxyv;
   float x, y, z, w, dxp, dyp, dzp, amx, amy, amz, dx1;
   float sq[MXV*MYV*MZV];
/* float sq[(mx+1)*(my+1)*(mz+1)]; */
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx + 1;
   myv = my + 1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   mxy1 = mx1*my1;
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
#pragma omp parallel for \
private(i,j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,nm,lm,x,y,z,w,dxp, \
dyp,dzp,amx,amy,amz,dx1,sq)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* zero out local accumulator */
      for (j = 0; j < mxyv*(mz+1); j++) {
         sq[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = nn - noff + mxv*(mm - moff) + mxyv*(ll - loff);
         amx = qm - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* deposit charge within tile to local accumulator */
         x = sq[nn] + amx*amz;
         y = sq[nn+1] + amy*amz;
         z = sq[nn+mxv] + dyp*amz;
         w = sq[nn+1+mxv] + dx1*amz;
         sq[nn] = x;
         sq[nn+1] = y;
         sq[nn+mxv] = z;
         sq[nn+1+mxv] = w;
         nn += mxyv;
         x = sq[nn] + amx*dzp;
         y = sq[nn+1] + amy*dzp;
         z = sq[nn+mxv] + dyp*dzp;
         w = sq[nn+1+mxv] + dx1*dzp;
         sq[nn] = x;
         sq[nn+1] = y;
         sq[nn+mxv] = z;
         sq[nn+1+mxv] = w;
      }
/* deposit charge to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
#pragma ivdep
            for (i = 1; i < nn; i++) {
               q[i+noff+nxv*(j+moff)+nxyv*(k+loff)]
               += sq[i+mxv*j+mxyv*k];
            }
         }
      }
/* deposit charge to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            q[i+noff+nxv*(j+moff)+nxyv*loff] += sq[i+mxv*j];
            if (lm > mz) {
#pragma omp atomic
               q[i+noff+nxv*(j+moff)+nxyv*(lm+loff-1)]
               += sq[i+mxv*j+mxyv*(lm-1)];
            }
         }
      }
      nm = nxv - noff;
      nm = mx+1 < nm ? mx+1 : nm;
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (k = 0; k < ll; k++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            q[i+noff+nxv*moff+nxyv*(k+loff)] += sq[i+mxyv*k];
            if (mm > my) {
#pragma omp atomic
               q[i+noff+nxv*(mm+moff-1)+nxyv*(k+loff)]
               += sq[i+mxv*(mm-1)+mxyv*k];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            q[noff+nxv*(j+moff)+nxyv*(k+loff)] += sq[mxv*j+mxyv*k];
            if (nm > mx) {
#pragma omp atomic
               q[nm+noff-1+nxv*(j+moff)+nxyv*(k+loff)]
               += sq[nm-1+mxv*j+mxyv*k];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            q[i+noff+nxv*moff+nxyv*(lm+loff-1)] += sq[i+mxyv*(lm-1)];
            if (mm > my) {
#pragma omp atomic
               q[i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1)]
               += sq[i+mxv*(mm-1)+mxyv*(lm-1)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            q[noff+nxv*(j+moff)+nxyv*(lm+loff-1)]
            += sq[mxv*j+mxyv*(lm-1)];
            if (nm > mx) {
#pragma omp atomic
               q[nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1)]
               += sq[nm-1+mxv*j+mxyv*(lm-1)];
            }
         }
      }
   }
   return;
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void cvgppost3lt(float ppart[], float q[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int my, int mz, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1) {
/* for 3d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   vectorizable/OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   33 flops/particle, 11 loads, 8 stores
   input: all, output: q
   charge density is approximated by values at the nearest grid points
   q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
   q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
   q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
   q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
   q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
   q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
   q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
   q(n+1,m+1,l+1)=qm*dx*dy*dz
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   q[l][k][j] = charge density at grid point j,k,l
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 6
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = first dimension of charge array, must be >= nx+1
   nyv = second dimension of charge array, must be >= ny+1
   nzv = third dimension of charge array, must be >= nz+1
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
#define NPBLK             32
#define LVECT             8
   int mxy1, noff, moff, loff, npoff, npp, ipp, joff, nps;
   int i, j, k, l, m, nn, mm, ll, nm, lm, lxv, lyv, lxyv, nxyv;
   float x, y, z, w, dxp, dyp, dzp, amx, amy, amz, dx1;
   float sq[MXV*MYV*MZV];
/* float sq[(mx+1)*(my+1)*(mz+1)]; */
/* scratch arrays */
   __attribute__((aligned(64))) int n[NPBLK], mn[LVECT];
   __attribute__((aligned(64))) float s[NPBLK*LVECT];
/* lxv = MXV; */
/* lyv = MYV; */
   mxy1 = mx1*my1;
   lxv = mx + 1;
   lyv = my + 1;
   lxyv = lxv*lyv;
   nxyv = nxv*nyv;
   mn[0] = 0;
   mn[1] = 1;
   mn[2] = lxv;
   mn[3] = lxv + 1;
   mn[4] = lxyv;
   mn[5] = lxyv + 1;
   mn[6] = lxyv + lxv;
   mn[7] = lxyv + lxv + 1;
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,ipp,joff,nps,nn,mm,ll,nm,lm, \
x,y,z,w,dxp,dyp,dzp,amx,amy,amz,dx1,sq,n,s)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* zero out local accumulator */
      for (j = 0; j < lxyv*(mz+1); j++) {
         sq[j] = 0.0f;
      }
/* loop over particles in tile */
      ipp = npp/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m;
/* inner loop over particles in block */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
/* find interpolation weights */
            x = ppart[j+joff+npoff];
            y = ppart[j+joff+nppmx+npoff];
            z = ppart[j+joff+2*nppmx+npoff];
            nn = x;
            mm = y;
            ll = z;
            dxp = qm*(x - (float) nn);
            dyp = y - (float) mm;
            dzp = z - (float) ll;
            n[j] = nn - noff + lxv*(mm - moff) + lxyv*(ll - loff);
            amx = qm - dxp;
            amy = 1.0f - dyp;
            dx1 = dxp*dyp;
            dyp = amx*dyp;
            amx = amx*amy;
            amz = 1.0f - dzp;
            amy = dxp*amy;
            s[j] = amx*amz;
            s[j+NPBLK] = amy*amz;
            s[j+2*NPBLK] = dyp*amz;
            s[j+3*NPBLK] = dx1*amz;
            s[j+4*NPBLK] = amx*dzp;
            s[j+5*NPBLK] = amy*dzp;
            s[j+6*NPBLK] = dyp*dzp;
            s[j+7*NPBLK] = dx1*dzp;
         }
/* deposit charge within tile to local accumulator */
         for (j = 0; j < NPBLK; j++) {
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               sq[n[j]+mn[i]] += s[j+NPBLK*i];
            }
         }
      }
      nps = NPBLK*ipp;
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = nn - noff + lxv*(mm - moff) + lxyv*(ll - loff);
         amx = qm - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* deposit charge within tile to local accumulator */
         x = sq[nn] + amx*amz;
         y = sq[nn+1] + amy*amz;
         z = sq[nn+lxv] + dyp*amz;
         w = sq[nn+1+lxv] + dx1*amz;
         sq[nn] = x;
         sq[nn+1] = y;
         sq[nn+lxv] = z;
         sq[nn+1+lxv] = w;
         nn += lxyv;
         x = sq[nn] + amx*dzp;
         y = sq[nn+1] + amy*dzp;
         z = sq[nn+lxv] + dyp*dzp;
         w = sq[nn+1+lxv] + dx1*dzp;
         sq[nn] = x;
         sq[nn+1] = y;
         sq[nn+lxv] = z;
         sq[nn+1+lxv] = w;
      }
/* deposit charge to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
#pragma ivdep
            for (i = 1; i < nn; i++) {
               q[i+noff+nxv*(j+moff)+nxyv*(k+loff)]
               += sq[i+lxv*j+lxyv*k];
            }
         }
      }
/* deposit charge to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            q[i+noff+nxv*(j+moff)+nxyv*loff] += sq[i+lxv*j];
            if (lm > mz) {
#pragma omp atomic
               q[i+noff+nxv*(j+moff)+nxyv*(lm+loff-1)]
               += sq[i+lxv*j+lxyv*(lm-1)];
            }
         }
      }
      nm = nxv - noff;
      nm = mx+1 < nm ? mx+1 : nm;
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (k = 0; k < ll; k++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            q[i+noff+nxv*moff+nxyv*(k+loff)] += sq[i+lxyv*k];
            if (mm > my) {
#pragma omp atomic
               q[i+noff+nxv*(mm+moff-1)+nxyv*(k+loff)]
               += sq[i+lxv*(mm-1)+lxyv*k];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            q[noff+nxv*(j+moff)+nxyv*(k+loff)] += sq[lxv*j+lxyv*k];
            if (nm > mx) {
#pragma omp atomic
               q[nm+noff-1+nxv*(j+moff)+nxyv*(k+loff)]
               += sq[nm-1+lxv*j+lxyv*k];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            q[i+noff+nxv*moff+nxyv*(lm+loff-1)] += sq[i+lxyv*(lm-1)];
            if (mm > my) {
#pragma omp atomic
               q[i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1)]
               += sq[i+lxv*(mm-1)+lxyv*(lm-1)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            q[noff+nxv*(j+moff)+nxyv*(lm+loff-1)]
            += sq[lxv*j+lxyv*(lm-1)];
            if (nm > mx) {
#pragma omp atomic
               q[nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1)]
               += sq[nm-1+lxv*j+lxyv*(lm-1)];
            }
         }
      }
   }
   return;
#undef LVECT
#undef NPBLK
#undef MXV
#undef MYV
#undef MZV
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
void cpporder3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int nx, int ny,
                 int nz, int mx, int my, int mz, int mx1, int my1,
                 int mz1, int npbmx, int ntmax, int *irc) {
/* this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 3D linear memory
   algorithm has 3 steps.  first, one finds particles leaving tile and
   stores their number in each directon, location, and destination in ncl
   and ihole.  second, a prefix scan of ncl is performed and departing
   particles are buffered in ppbuff in direction order.  finally, we copy
   the incoming particles from other tiles into ppart.
   input: all except ppbuff, ncl, ihole, irc
   output: ppart, ppbuff, kpic, ncl, ihole, irc
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppbuff[m][i][n] = i co-ordinate of particle n in tile m
   kpic[m] = number of particles in tile m
   ncl[m][i] = number of particles going to destination i, tile m
   ihole[m][:][0] = location of hole in array left by departing particle
   ihole[m][:][1] = direction destination of particle leaving hole
   all for tile m
   ihole[m][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mz1 = (system length in z direction - 1)/mz + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int mxy1, mxyz1, noff, moff, loff, npoff, npp, nboff, ncoff;
   int i, j, k, l, ii, kx, ky, kz, ih, nh, ist, nn, mm, ll, isum;
   int ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dx, dy, dz;
   int ks[26];
   mxy1 = mx1*my1;
   mxyz1 = mxy1*mz1;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
/* find and count particles leaving tiles and determine destination */
/* update ppart, ihole, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,ih,nh,ist,dx,dy,dz, \
edgelx,edgely,edgelz,edgerx,edgery,edgerz)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      nn = nx - noff;
      nn = mx < nn ? mx : nn;
      mm = ny - moff;
      mm = my < mm ? my : mm;
      ll = nz - loff;
      ll = mz < ll ? mz : ll;
      ih = 0;
      nh = 0;
      edgelx = noff;
      edgerx = noff + nn;
      edgely = moff;
      edgery = moff + mm;
      edgelz = loff;
      edgerz = loff + ll;
/* clear counters */
      for (j = 0; j < 26; j++) {
         ncl[j+26*l] = 0;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
         dx = ppart[j+npoff];
         dy = ppart[j+nppmx+npoff];
         dz = ppart[j+2*nppmx+npoff];
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
         if (dz >= edgerz) {
            if (dz >= anz)
               ppart[j+2*nppmx+npoff] = dz - anz;
            ist += 18;
         }
         else if (dz < edgelz) {
            if (dz < 0.0) {
               dz += anz;
               if (dz < anz)
                  ist += 9;
               else
                  dz = 0.0;
               ppart[j+2*nppmx+npoff] = dz;
            }
            else {
               ist += 9;
            }
         }
         if (ist > 0) {
            ncl[ist+26*l-1] += 1;
            ih += 1;
            if (ih <= ntmax) {
               ihole[2*(ih+(ntmax+1)*l)] = j + 1;
               ihole[1+2*(ih+(ntmax+1)*l)] = ist;
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
      ihole[2*(ntmax+1)*l] = ih;
   }
/* ihole overflow */
   if (*irc > 0)
      return;

/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,l,npoff,nboff,isum,ist,nh,ip,j1,ii)
   for (l = 0; l < mxyz1; l++) {
      npoff = idimp*nppmx*l;
      nboff = idimp*npbmx*l;
/* find address offset for ordered ppbuff array */
      isum = 0;
      for (j = 0; j < 26; j++) {
         ist = ncl[j+26*l];
         ncl[j+26*l] = isum;
         isum += ist;
      }
      nh = ihole[2*(ntmax+1)*l];
      ip = 0;
/* loop over particles leaving tile */
      for (j = 0; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+(ntmax+1)*l)] - 1;
         ist = ihole[1+2*(j+1+(ntmax+1)*l)];
         ii = ncl[ist+26*l-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*i+nboff]
               = ppart[j1+nppmx*i+npoff];
            }
         }
         else {
            ip = 1;
         }
         ncl[ist+26*l-1] = ii + 1;
      }
/* set error */
      if (ip > 0)
         *irc = ncl[25+26*l];
   }
/* ppbuff overflow */
   if (*irc > 0)
      return;

/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,ii,kk,npp,npoff,nboff,kx,ky,kz,kl,kr,kxl,kxr,lk,ll,lr, \
ih,nh,nn,ncoff,ist,j1,j2,ip,ks)
   for (l = 0; l < mxyz1; l++) {
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      kz = l/mxy1;
      k = l - mxy1*kz;
/* loop over tiles in z, assume periodic boundary conditions */
      lk = kz*mxy1;
/* find tile behind */
      ll = kz - 1;
      if (ll < 0)
         ll += mz1;
      ll = ll*mxy1;
/* find tile in front */
      lr = kz + 1;
      if (lr >= mz1)
         lr -= mz1;
      lr = lr*mxy1;
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
      kxl = kx - 1 ;
      if (kxl < 0)
         kxl += mx1;
      kxr = kx + 1;
      if (kxr >= mx1)
         kxr -= mx1;
/* find tile number for different directions */
      ks[0] = kxr + kk + lk;
      ks[1] = kxl + kk + lk;
      ks[2] = kx + kr + lk;
      ks[3] = kxr + kr + lk;
      ks[4] = kxl + kr + lk;
      ks[5] = kx + kl + lk;
      ks[6] = kxr + kl + lk;
      ks[7] = kxl + kl + lk;
      ks[8] = kx + kk + lr;
      ks[9] = kxr + kk + lr;
      ks[10] = kxl + kk + lr;
      ks[11] = kx + kr + lr;
      ks[12] = kxr + kr + lr;
      ks[13] = kxl + kr + lr;
      ks[14] = kx + kl + lr;
      ks[15] = kxr + kl + lr;
      ks[16] = kxl + kl + lr;
      ks[17] = kx + kk + ll;
      ks[18] = kxr + kk + ll;
      ks[19] = kxl + kk + ll;
      ks[20] = kx + kr + ll;
      ks[21] = kxr + kr + ll;
      ks[22] = kxl + kr + ll;
      ks[23] = kx + kl + ll;
      ks[24] = kxr + kl + ll;
      ks[25] = kxl + kl + ll;
/* loop over directions */
      nh = ihole[2*(ntmax+1)*l];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      for (ii = 0; ii < 26; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+26*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+26*ks[ii]] - ncoff;
         for (j = 0; j < ip; j++) {
            ih += 1;
/* insert incoming particles into holes */
            if (ih <= nh) {
               j1 = ihole[2*(ih+(ntmax+1)*l)] - 1;
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
/* holes with locations great than npp-ip do not need to be filled */
      if (ih < nh) {
         ip = nh - ih;
         ii = nh;
         nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
         ih += 1;
         j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
/* move particles from end into remaining holes */
/* holes are processed in increasing order */
         for (j = 0; j < ip; j++) {
            j1 = npp - j - 1;
            if (j1==nn) {
               ii -= 1;
               nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
            }
            else {
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*i+npoff]
                  = ppart[j1+nppmx*i+npoff];
               }
               ih += 1;
               j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[l] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int mx1, int my1,
                  int mz1, int npbmx, int ntmax, int *irc) {
/* this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 3D linear memory
   the algorithm has 2 steps.  first, a prefix scan of ncl is performed
   and departing particles are buffered in ppbuff in direction order.
   then we copy the incoming particles from other tiles into ppart.
   it assumes that the number, location, and destination of particles 
   leaving a tile have been previously stored in ncl and ihole by the
   cgppushf3lt subroutine.
   input: all except ppbuff, irc
   output: ppart, ppbuff, kpic, ncl, irc
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppbuff[m][i][n] = i co-ordinate of particle n in tile m
   kpic[m] = number of particles in tile m
   ncl[m][i] = number of particles going to destination i, tile m
   ihole[m][:][0] = location of hole in array left by departing particle
   ihole[m][:][1] = direction destination of particle leaving hole
   all for tile m
   ihole[m][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mz1 = (system length in z direction - 1)/mz + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int mxy1, mxyz1, npp, npoff, nboff, ncoff;
   int i, j, k, l, ii, kx, ky, kz, ih, nh, ist, nn, ll, isum;
   int ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr;
   int ks[26];
   mxy1 = mx1*my1;
   mxyz1 = mxy1*mz1;
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,l,npoff,nboff,isum,ist,nh,ip,j1,ii)
   for (l = 0; l < mxyz1; l++) {
      npoff = idimp*nppmx*l;
      nboff = idimp*npbmx*l;
/* find address offset for ordered ppbuff array */
      isum = 0;
      for (j = 0; j < 26; j++) {
         ist = ncl[j+26*l];
         ncl[j+26*l] = isum;
         isum += ist;
      }
      nh = ihole[2*(ntmax+1)*l];
      ip = 0;
/* loop over particles leaving tile */
      for (j = 0; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+(ntmax+1)*l)] - 1;
         ist = ihole[1+2*(j+1+(ntmax+1)*l)];
         ii = ncl[ist+26*l-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*i+nboff]
               = ppart[j1+nppmx*i+npoff];
            }
         }
         else {
            ip = 1;
         }
         ncl[ist+26*l-1] = ii + 1;
      }
/* set error */
      if (ip > 0)
         *irc = ncl[25+26*l];
   }
/* ppbuff overflow */
   if (*irc > 0)
      return;

/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,ii,kk,npp,npoff,nboff,kx,ky,kz,kl,kr,kxl,kxr,lk,ll,lr, \
ih,nh,nn,ncoff,ist,j1,j2,ip,ks)
   for (l = 0; l < mxyz1; l++) {
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      kz = l/mxy1;
      k = l - mxy1*kz;
/* loop over tiles in z, assume periodic boundary conditions */
      lk = kz*mxy1;
/* find tile behind */
      ll = kz - 1;
      if (ll < 0)
         ll += mz1;
      ll = ll*mxy1;
/* find tile in front */
      lr = kz + 1;
      if (lr >= mz1)
         lr -= mz1;
      lr = lr*mxy1;
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
      kxl = kx - 1 ;
      if (kxl < 0)
         kxl += mx1;
      kxr = kx + 1;
      if (kxr >= mx1)
         kxr -= mx1;
/* find tile number for different directions */
      ks[0] = kxr + kk + lk;
      ks[1] = kxl + kk + lk;
      ks[2] = kx + kr + lk;
      ks[3] = kxr + kr + lk;
      ks[4] = kxl + kr + lk;
      ks[5] = kx + kl + lk;
      ks[6] = kxr + kl + lk;
      ks[7] = kxl + kl + lk;
      ks[8] = kx + kk + lr;
      ks[9] = kxr + kk + lr;
      ks[10] = kxl + kk + lr;
      ks[11] = kx + kr + lr;
      ks[12] = kxr + kr + lr;
      ks[13] = kxl + kr + lr;
      ks[14] = kx + kl + lr;
      ks[15] = kxr + kl + lr;
      ks[16] = kxl + kl + lr;
      ks[17] = kx + kk + ll;
      ks[18] = kxr + kk + ll;
      ks[19] = kxl + kk + ll;
      ks[20] = kx + kr + ll;
      ks[21] = kxr + kr + ll;
      ks[22] = kxl + kr + ll;
      ks[23] = kx + kl + ll;
      ks[24] = kxr + kl + ll;
      ks[25] = kxl + kl + ll;
/* loop over directions */
      nh = ihole[2*(ntmax+1)*l];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      for (ii = 0; ii < 26; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+26*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+26*ks[ii]] - ncoff;
         for (j = 0; j < ip; j++) {
            ih += 1;
/* insert incoming particles into holes */
            if (ih <= nh) {
               j1 = ihole[2*(ih+(ntmax+1)*l)] - 1;
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
/* holes with locations great than npp-ip do not need to be filled */
      if (ih < nh) {
         ip = nh - ih;
         ii = nh;
         nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
         ih += 1;
         j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
/* move particles from end into remaining holes */
/* holes are processed in increasing order */
         for (j = 0; j < ip; j++) {
            j1 = npp - j - 1;
            if (j1==nn) {
               ii -= 1;
               nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
            }
            else {
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*i+npoff]
                  = ppart[j1+nppmx*i+npoff];
               }
               ih += 1;
               j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[l] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cvpporder3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int nx, int ny,
                  int nz, int mx, int my, int mz, int mx1, int my1,
                  int mz1, int npbmx, int ntmax, int *irc) {
/* this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 3D linear memory
   algorithm has 3 steps.  first, one finds particles leaving tile and
   stores their number in each directon, location, and destination in ncl
   and ihole.  second, a prefix scan of ncl is performed and departing
   particles are buffered in ppbuff in direction order.  finally, we copy
   the incoming particles from other tiles into ppart.
   input: all except ppbuff, ncl, ihole, irc
   output: ppart, ppbuff, kpic, ncl, ihole, irc
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppbuff[m][i][n] = i co-ordinate of particle n in tile m
   kpic[m] = number of particles in tile m
   ncl[m][i] = number of particles going to destination i, tile m
   ihole[m][:][0] = location of hole in array left by departing particle
   ihole[m][:][1] = direction destination of particle leaving hole
   all for tile m
   ihole[m][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mz1 = (system length in z direction - 1)/mz + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
#define NPBLK             16
   int mxy1, mxyz1, noff, moff, loff, npp, npoff, nboff, ncoff;
   int ipp, joff, nps;
   int i, j, k, l, m, ii, kx, ky, kz, ih, nh, ist, nn, mm, ll, in;
   int ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr, lb, kxs;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dx, dy, dz;
   __attribute__((aligned(64))) int sncl[26], ks[26];
/* scratch arrays */
   __attribute__((aligned(64))) int n[NPBLK*3];
   mxy1 = mx1*my1;
   mxyz1 = mxy1*mz1;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
/* find and count particles leaving tiles and determine destination */
/* update ppart, ihole, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,ih,nh,ist,dx,dy,dz, \
edgelx,edgely,edgelz,edgerx,edgery,edgerz,n)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      nn = nx - noff;
      nn = mx < nn ? mx : nn;
      mm = ny - moff;
      mm = my < mm ? my : mm;
      ll = nz - loff;
      ll = mz < ll ? mz : ll;
      ih = 0;
      nh = 0;
      edgelx = noff;
      edgerx = noff + nn;
      edgely = moff;
      edgery = moff + mm;
      edgelz = loff;
      edgerz = loff + ll;
/* clear counters */
      for (j = 0; j < 26; j++) {
         ncl[j+26*l] = 0;
      }
/* loop over particles in tile */
      ipp = npp/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m;
/* inner loop over particles in block */
#pragma vector aligned
         for (j = 0; j < NPBLK; j++) {
            dx = ppart[j+joff+npoff];
            dy = ppart[j+joff+nppmx+npoff];
            dz = ppart[j+joff+2*nppmx+npoff];
/* find particles going out of bounds */
            ist = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* ist = direction particle is going                             */
            if (dx >= edgerx) {
               if (dx >= anx)
                  ppart[j+joff+npoff] = dx - anx;
               ist = 2;
            }
            else if (dx < edgelx) {
               if (dx < 0.0) {
                  dx += anx;
                  if (dx < anx)
                     ist = 1;
                  else
                     dx = 0.0;
                  ppart[j+joff+npoff] = dx;
              }
               else {
                  ist = 1;
               }
            }
            if (dy >= edgery) {
               if (dy >= any)
                  ppart[j+joff+nppmx+npoff] = dy - any;
               ist += 6;
            }
            else if (dy < edgely) {
               if (dy < 0.0) {
                  dy += any;
                  if (dy < any)
                     ist += 3;
                  else
                     dy = 0.0;
                  ppart[j+joff+nppmx+npoff] = dy;
               }
               else {
                  ist += 3;
               }
            }
            if (dz >= edgerz) {
               if (dz >= anz)
                  ppart[j+joff+2*nppmx+npoff] = dz - anz;
               ist += 18;
            }
            else if (dz < edgelz) {
               if (dz < 0.0) {
                     dz += anz;
               if (dz < anz)
                     ist += 9;
                  else
                     dz = 0.0;
                  ppart[j+joff+2*nppmx+npoff] = dz;
               }
               else {
                  ist += 9;
               }
            }
            n[j] = ist;
         }
/* store outgoing particle address and destination */
         for (j = 0; j < NPBLK; j++) {
            ist = n[j];
            if (ist > 0) {
               ncl[ist+26*l-1] += 1;
               ih += 1;
               if (ih <= ntmax) {
                  ihole[2*(ih+(ntmax+1)*l)] = j + joff + 1;
                  ihole[1+2*(ih+(ntmax+1)*l)] = ist;
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
         dx = ppart[j+npoff];
         dy = ppart[j+nppmx+npoff];
         dz = ppart[j+2*nppmx+npoff];
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
         if (dz >= edgerz) {
            if (dz >= anz)
               ppart[j+2*nppmx+npoff] = dz - anz;
            ist += 18;
         }
         else if (dz < edgelz) {
            if (dz < 0.0) {
               dz += anz;
               if (dz < anz)
                  ist += 9;
               else
                  dz = 0.0;
               ppart[j+2*nppmx+npoff] = dz;
            }
            else {
               ist += 9;
            }
         }
         if (ist > 0) {
            ncl[ist+26*l-1] += 1;
            ih += 1;
            if (ih <= ntmax) {
               ihole[2*(ih+(ntmax+1)*l)] = j + 1;
               ihole[1+2*(ih+(ntmax+1)*l)] = ist;
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
      ihole[2*(ntmax+1)*l] = ih;
   }
/* ihole overflow */
   if (*irc > 0)
      return;

/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,l,m,npoff,nboff,kxs,lb,ist,nh,ip,ipp,nps,joff,j1,ii,sncl, \
ks,n)
   for (l = 0; l < mxyz1; l++) {
      npoff = idimp*nppmx*l;
      nboff = idimp*npbmx*l;
/* find address offset for ordered ppbuff array */
/* find address offset for ordered ppbuff array */
      for (j = 0; j < 26; j++) {
         sncl[j] = ncl[j+26*l];
         ks[j] = j;
      }
      kxs = 1;
      while (kxs < 26) {
#pragma ivdep
         for (j = 0; j < 13; j++) {
            lb = kxs*ks[j];
            if ((j+lb+kxs) < 26)
               sncl[j+lb+kxs] += sncl[2*lb+kxs-1];
            ks[j] >>= 1;
         }     
         kxs <<= 1;
      }
      for (j = 0; j < 26; j++) {
         sncl[j] -= ncl[j+26*l];
      }
      nh = ihole[2*(ntmax+1)*l];
      ip = 0;
/* buffer particles that are leaving tile, in direction order */
/* loop over particles leaving tile */
      ipp = nh/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m + 1;
/* inner loop over particles in block */
         for (j = 0; j < NPBLK; j++) {
            n[j] = ihole[2*(j+joff+(ntmax+1)*l)] - 1;
            n[j+NPBLK] = ihole[1+2*(j+joff+(ntmax+1)*l)];
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
         j1 = ihole[2*(j+1+(ntmax+1)*l)] - 1;
         ist = ihole[1+2*(j+1+(ntmax+1)*l)];
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
      for (j = 0; j < 26; j++) {
         ncl[j+26*l] = sncl[j];
      }
/* set error */
      if (ip > 0)
         *irc = ncl[25+26*l];
   }
/* ppbuff overflow */
   if (*irc > 0)
      return;

/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,ii,kk,in,npp,npoff,nboff,ipp,joff,nps,kx,ky,kz,kl,kr, \
kxl,kxr,lk,ll,lr,ih,nh,nn,mm,ncoff,ist,j1,j2,ip,ks,n)
   for (l = 0; l < mxyz1; l++) {
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      kz = l/mxy1;
      k = l - mxy1*kz;
/* loop over tiles in z, assume periodic boundary conditions */
      lk = kz*mxy1;
/* find tile behind */
      ll = kz - 1;
      if (ll < 0)
         ll += mz1;
      ll = ll*mxy1;
/* find tile in front */
      lr = kz + 1;
      if (lr >= mz1)
         lr -= mz1;
      lr = lr*mxy1;
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
      kxl = kx - 1 ;
      if (kxl < 0)
         kxl += mx1;
      kxr = kx + 1;
      if (kxr >= mx1)
         kxr -= mx1;
/* find tile number for different directions */
      ks[0] = kxr + kk + lk;
      ks[1] = kxl + kk + lk;
      ks[2] = kx + kr + lk;
      ks[3] = kxr + kr + lk;
      ks[4] = kxl + kr + lk;
      ks[5] = kx + kl + lk;
      ks[6] = kxr + kl + lk;
      ks[7] = kxl + kl + lk;
      ks[8] = kx + kk + lr;
      ks[9] = kxr + kk + lr;
      ks[10] = kxl + kk + lr;
      ks[11] = kx + kr + lr;
      ks[12] = kxr + kr + lr;
      ks[13] = kxl + kr + lr;
      ks[14] = kx + kl + lr;
      ks[15] = kxr + kl + lr;
      ks[16] = kxl + kl + lr;
      ks[17] = kx + kk + ll;
      ks[18] = kxr + kk + ll;
      ks[19] = kxl + kk + ll;
      ks[20] = kx + kr + ll;
      ks[21] = kxr + kr + ll;
      ks[22] = kxl + kr + ll;
      ks[23] = kx + kl + ll;
      ks[24] = kxr + kl + ll;
      ks[25] = kxl + kl + ll;
/* loop over directions */
      nh = ihole[2*(ntmax+1)*l];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      for (ii = 0; ii < 26; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+26*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+26*ks[ii]] - ncoff;
/* loop over particles coming from direction ii */
         ipp = ip/NPBLK;
/* outer loop over number of full blocks */
         for (m = 0; m < ipp; m++) {
            joff = NPBLK*m;
/* inner loop over particles in block */
            for (j = 0; j < NPBLK; j++) {
/* insert incoming particles into holes */
               if ((j+ih) < nh) {
                  j1 = ihole[2*(j+ih+1+(ntmax+1)*l)] - 1;
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
               j1 = ihole[2*(ih+(ntmax+1)*l)] - 1;
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
/* holes with locations great than npp-ip do not need to be filled */
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
               n[j+NPBLK] = ihole[2*(ih+j+1+(ntmax+1)*l)] - 1;
               n[j+2*NPBLK] = ihole[2*(ii-j+(ntmax+1)*l)] - 1;
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
         nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
         ih += 1;
         j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
/* loop over remaining particles */
         for (j = nps; j < ip; j++) {
            j1 = npp - j - 1;
            if (j1==nn) {
               ii -= 1;
               nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
            }
            else {
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*i+npoff]
                  = ppart[j1+nppmx*i+npoff];
               }
               ih += 1;
               j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[l] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cvpporderf3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                   int ihole[], int idimp, int nppmx, int mx1, int my1,
                   int mz1, int npbmx, int ntmax, int *irc) {
/* this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 3D linear memory
   the algorithm has 2 steps.  first, a prefix scan of ncl is performed
   and departing particles are buffered in ppbuff in direction order.
   then we copy the incoming particles from other tiles into ppart.
   it assumes that the number, location, and destination of particles 
   leaving a tile have been previously stored in ncl and ihole by the
   cvgppushf3lt subroutine.
   input: all except ppbuff, irc
   output: ppart, ppbuff, kpic, ncl, irc
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppbuff[m][i][n] = i co-ordinate of particle n in tile m
   kpic[m] = number of particles in tile m
   ncl[m][i] = number of particles going to destination i, tile m
   ihole[m][:][0] = location of hole in array left by departing particle
   ihole[m][:][1] = direction destination of particle leaving hole
   all for tile m
   ihole[m][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mz1 = (system length in z direction - 1)/mz + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
#define NPBLK             16
   int mxy1, mxyz1, npp, ncoff, npoff, nboff;
   int i, j, k, l, ii, kx, ky, kz, ih, nh, ist, nn, ll, mm, in;
   int ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr;
   int lb, kxs, m, ipp, nps, joff;
   __attribute__((aligned(64))) int sncl[26], ks[26];
/* scratch arrays */
   __attribute__((aligned(64))) int n[NPBLK*3];
   mxy1 = mx1*my1;
   mxyz1 = mxy1*mz1;
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,l,m,npoff,nboff,kxs,lb,ist,nh,ip,ipp,nps,joff,j1,ii,sncl, \
ks,n)
   for (l = 0; l < mxyz1; l++) {
      npoff = idimp*nppmx*l;
      nboff = idimp*npbmx*l;
/* find address offset for ordered ppbuff array */
/* find address offset for ordered ppbuff array */
      for (j = 0; j < 26; j++) {
         sncl[j] = ncl[j+26*l];
         ks[j] = j;
      }
      kxs = 1;
      while (kxs < 26) {
#pragma ivdep
         for (j = 0; j < 13; j++) {
            lb = kxs*ks[j];
            if ((j+lb+kxs) < 26)
               sncl[j+lb+kxs] += sncl[2*lb+kxs-1];
            ks[j] >>= 1;
         }     
         kxs <<= 1;
      }
      for (j = 0; j < 26; j++) {
         sncl[j] -= ncl[j+26*l];
      }
      nh = ihole[2*(ntmax+1)*l];
      ip = 0;
/* buffer particles that are leaving tile, in direction order */
/* loop over particles leaving tile */
      ipp = nh/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m + 1;
/* inner loop over particles in block */
         for (j = 0; j < NPBLK; j++) {
            n[j] = ihole[2*(j+joff+(ntmax+1)*l)] - 1;
            n[j+NPBLK] = ihole[1+2*(j+joff+(ntmax+1)*l)];
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
         j1 = ihole[2*(j+1+(ntmax+1)*l)] - 1;
         ist = ihole[1+2*(j+1+(ntmax+1)*l)];
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
      for (j = 0; j < 26; j++) {
         ncl[j+26*l] = sncl[j];
      }
/* set error */
      if (ip > 0)
         *irc = ncl[25+26*l];
   }
/* ppbuff overflow */
   if (*irc > 0)
      return;

/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,ii,kk,in,npp,npoff,nboff,ipp,joff,nps,kx,ky,kz,kl,kr, \
kxl,kxr,lk,ll,lr,ih,nh,nn,mm,ncoff,ist,j1,j2,ip,ks,n)
   for (l = 0; l < mxyz1; l++) {
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      kz = l/mxy1;
      k = l - mxy1*kz;
/* loop over tiles in z, assume periodic boundary conditions */
      lk = kz*mxy1;
/* find tile behind */
      ll = kz - 1;
      if (ll < 0)
         ll += mz1;
      ll = ll*mxy1;
/* find tile in front */
      lr = kz + 1;
      if (lr >= mz1)
         lr -= mz1;
      lr = lr*mxy1;
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
      kxl = kx - 1 ;
      if (kxl < 0)
         kxl += mx1;
      kxr = kx + 1;
      if (kxr >= mx1)
         kxr -= mx1;
/* find tile number for different directions */
      ks[0] = kxr + kk + lk;
      ks[1] = kxl + kk + lk;
      ks[2] = kx + kr + lk;
      ks[3] = kxr + kr + lk;
      ks[4] = kxl + kr + lk;
      ks[5] = kx + kl + lk;
      ks[6] = kxr + kl + lk;
      ks[7] = kxl + kl + lk;
      ks[8] = kx + kk + lr;
      ks[9] = kxr + kk + lr;
      ks[10] = kxl + kk + lr;
      ks[11] = kx + kr + lr;
      ks[12] = kxr + kr + lr;
      ks[13] = kxl + kr + lr;
      ks[14] = kx + kl + lr;
      ks[15] = kxr + kl + lr;
      ks[16] = kxl + kl + lr;
      ks[17] = kx + kk + ll;
      ks[18] = kxr + kk + ll;
      ks[19] = kxl + kk + ll;
      ks[20] = kx + kr + ll;
      ks[21] = kxr + kr + ll;
      ks[22] = kxl + kr + ll;
      ks[23] = kx + kl + ll;
      ks[24] = kxr + kl + ll;
      ks[25] = kxl + kl + ll;
/* loop over directions */
      nh = ihole[2*(ntmax+1)*l];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      for (ii = 0; ii < 26; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+26*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+26*ks[ii]] - ncoff;
/* loop over particles coming from direction ii */
         ipp = ip/NPBLK;
/* outer loop over number of full blocks */
         for (m = 0; m < ipp; m++) {
            joff = NPBLK*m;
/* inner loop over particles in block */
            for (j = 0; j < NPBLK; j++) {
/* insert incoming particles into holes */
               if ((j+ih) < nh) {
                  j1 = ihole[2*(j+ih+1+(ntmax+1)*l)] - 1;
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
               j1 = ihole[2*(ih+(ntmax+1)*l)] - 1;
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
/* holes with locations great than npp-ip do not need to be filled */
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
               n[j+NPBLK] = ihole[2*(ih+j+1+(ntmax+1)*l)] - 1;
               n[j+2*NPBLK] = ihole[2*(ii-j+(ntmax+1)*l)] - 1;
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
         nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
         ih += 1;
         j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
/* loop over remaining particles */
         for (j = nps; j < ip; j++) {
            j1 = npp - j - 1;
            if (j1==nn) {
               ii -= 1;
               nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
            }
            else {
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*i+npoff]
                  = ppart[j1+nppmx*i+npoff];
               }
               ih += 1;
               j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[l] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cv2pporderf3lt(float ppart[], float ppbuff[], int kpic[],
                    int ncl[], int ihole[], int idimp, int nppmx,
                    int mx1, int my1, int mz1, int npbmx, int ntmax,
                    int *irc) {
/* this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 3D linear memory
   the algorithm has 2 steps.  first, a prefix scan of ncl is performed
   and departing particles are buffered in ppbuff in direction order.
   then we copy the incoming particles from other tiles into ppart.
   it assumes that the number, location, and destination of particles 
   leaving a tile have been previously stored in ncl and ihole by the
   cvgppushf3lt subroutine.
   input: all except ppbuff, irc
   output: ppart, ppbuff, kpic, ncl, irc
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppbuff[m][n][i] = i co-ordinate of particle n in tile m
   kpic[m] = number of particles in tile m
   ncl[m][i] = number of particles going to destination i, tile m
   ihole[m][:][0] = location of hole in array left by departing particle
   ihole[m][:][1] = direction destination of particle leaving hole
   all for tile m
   ihole[m][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mz1 = (system length in z direction - 1)/mz + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
#define NPBLK             16
   int mxy1, mxyz1, npp, ncoff, npoff, nboff;
   int i, j, k, l, ii, kx, ky, kz, ih, nh, ist, nn, ll, mm, in;
   int ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr;
   int lb, kxs, m, ipp, nps, joff;
   __attribute__((aligned(64))) int sncl[26], ks[26];
/* scratch arrays */
   __attribute__((aligned(64))) int n[NPBLK*3];
   mxy1 = mx1*my1;
   mxyz1 = mxy1*mz1;
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,l,m,npoff,nboff,kxs,lb,ist,nh,ip,ipp,nps,joff,j1,ii,sncl, \
ks,n)
   for (l = 0; l < mxyz1; l++) {
      npoff = idimp*nppmx*l;
      nboff = idimp*npbmx*l;
/* find address offset for ordered ppbuff array */
/* find address offset for ordered ppbuff array */
      for (j = 0; j < 26; j++) {
         sncl[j] = ncl[j+26*l];
         ks[j] = j;
      }
      kxs = 1;
      while (kxs < 26) {
#pragma ivdep
         for (j = 0; j < 13; j++) {
            lb = kxs*ks[j];
            if ((j+lb+kxs) < 26)
               sncl[j+lb+kxs] += sncl[2*lb+kxs-1];
            ks[j] >>= 1;
         }     
         kxs <<= 1;
      }
      for (j = 0; j < 26; j++) {
         sncl[j] -= ncl[j+26*l];
      }
      nh = ihole[2*(ntmax+1)*l];
      ip = 0;
/* buffer particles that are leaving tile, in direction order */
/* loop over particles leaving tile */
      ipp = nh/NPBLK;
/* outer loop over number of full blocks */
      for (m = 0; m < ipp; m++) {
         joff = NPBLK*m + 1;
/* inner loop over particles in block */
         for (j = 0; j < NPBLK; j++) {
            n[j] = ihole[2*(j+joff+(ntmax+1)*l)] - 1;
            n[j+NPBLK] = ihole[1+2*(j+joff+(ntmax+1)*l)];
         }
/* calculate offsets */
         for (j = 0; j < NPBLK; j++) {
            ist = n[j+NPBLK];
            ii = sncl[ist-1];
            n[j+NPBLK] = ii;
            sncl[ist-1] = ii + 1;
         }
/* buffer particles that are leaving tile, in direction order */
         for (j = 0; j < NPBLK; j++) {
            j1 = n[j];
            ii = n[j+NPBLK];
            if (ii < npbmx) {
               for (i = 0; i < idimp; i++) {
                 ppbuff[i+idimp*ii+nboff]
                  = ppart[j1+nppmx*i+npoff];
               }
            }
            else {
               ip = 1;
            }
         }
      }
      nps = NPBLK*ipp;
/* loop over remaining particles */
      for (j = nps; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+(ntmax+1)*l)] - 1;
         ist = ihole[1+2*(j+1+(ntmax+1)*l)];
         ii = sncl[ist-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[i+idimp*ii+nboff]
               = ppart[j1+nppmx*i+npoff];
            }
         }
         else {
            ip = 1;
         }
         sncl[ist-1] = ii + 1;
      }
      for (j = 0; j < 26; j++) {
         ncl[j+26*l] = sncl[j];
      }
/* set error */
      if (ip > 0)
         *irc = ncl[25+26*l];
   }
/* ppbuff overflow */
   if (*irc > 0)
      return;

/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,ii,kk,in,npp,npoff,nboff,ipp,joff,nps,kx,ky,kz,kl,kr, \
kxl,kxr,lk,ll,lr,ih,nh,nn,mm,ncoff,ist,j1,j2,ip,ks,n)
   for (l = 0; l < mxyz1; l++) {
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      kz = l/mxy1;
      k = l - mxy1*kz;
/* loop over tiles in z, assume periodic boundary conditions */
      lk = kz*mxy1;
/* find tile behind */
      ll = kz - 1;
      if (ll < 0)
         ll += mz1;
      ll = ll*mxy1;
/* find tile in front */
      lr = kz + 1;
      if (lr >= mz1)
         lr -= mz1;
      lr = lr*mxy1;
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
      kxl = kx - 1 ;
      if (kxl < 0)
         kxl += mx1;
      kxr = kx + 1;
      if (kxr >= mx1)
         kxr -= mx1;
/* find tile number for different directions */
      ks[0] = kxr + kk + lk;
      ks[1] = kxl + kk + lk;
      ks[2] = kx + kr + lk;
      ks[3] = kxr + kr + lk;
      ks[4] = kxl + kr + lk;
      ks[5] = kx + kl + lk;
      ks[6] = kxr + kl + lk;
      ks[7] = kxl + kl + lk;
      ks[8] = kx + kk + lr;
      ks[9] = kxr + kk + lr;
      ks[10] = kxl + kk + lr;
      ks[11] = kx + kr + lr;
      ks[12] = kxr + kr + lr;
      ks[13] = kxl + kr + lr;
      ks[14] = kx + kl + lr;
      ks[15] = kxr + kl + lr;
      ks[16] = kxl + kl + lr;
      ks[17] = kx + kk + ll;
      ks[18] = kxr + kk + ll;
      ks[19] = kxl + kk + ll;
      ks[20] = kx + kr + ll;
      ks[21] = kxr + kr + ll;
      ks[22] = kxl + kr + ll;
      ks[23] = kx + kl + ll;
      ks[24] = kxr + kl + ll;
      ks[25] = kxl + kl + ll;
/* loop over directions */
      nh = ihole[2*(ntmax+1)*l];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      for (ii = 0; ii < 26; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+26*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+26*ks[ii]] - ncoff;
/* loop over particles coming from direction ii */
         ipp = ip/NPBLK;
/* outer loop over number of full blocks */
         for (m = 0; m < ipp; m++) {
            joff = NPBLK*m;
/* inner loop over particles in block */
            for (j = 0; j < NPBLK; j++) {
/* insert incoming particles into holes */
               if ((j+ih) < nh) {
                  j1 = ihole[2*(j+ih+1+(ntmax+1)*l)] - 1;
               }
/* place overflow at end of array */
               else {
                  j1 = npp + j + ih - nh;
               }
               n[j] = j1;
            }
            for (j = 0; j < NPBLK; j++) {
               j1 = n[j];
               if (j1 < nppmx) {
                  for (i = 0; i < idimp; i++) {
                     ppart[j1+nppmx*i+npoff]
                     = ppbuff[i+idimp*(j+joff+ncoff)+nboff];
                  }
               }
               else {
                  ist = 1;
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
               j1 = ihole[2*(ih+(ntmax+1)*l)] - 1;
            }
/* place overflow at end of array */
            else {
               j1 = npp + ih - nh - 1;
            }
            if (j1 < nppmx) {
               for (i = 0; i < idimp; i++) {
                  ppart[j1+nppmx*i+npoff]
                  = ppbuff[i+idimp*(j+ncoff)+nboff];
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
/* holes with locations great than npp-ip do not need to be filled */
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
               n[j+NPBLK] = ihole[2*(ih+j+1+(ntmax+1)*l)] - 1;
               n[j+2*NPBLK] = ihole[2*(ii-j+(ntmax+1)*l)] - 1;
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
         nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
         ih += 1;
         j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
/* loop over remaining particles */
         for (j = nps; j < ip; j++) {
            j1 = npp - j - 1;
            if (j1==nn) {
               ii -= 1;
               nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
            }
            else {
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*i+npoff]
                  = ppart[j1+nppmx*i+npoff];
               }
               ih += 1;
               j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[l] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ccguard3l(float fxyz[], int nx, int ny, int nz, int nxe, int nye,
               int nze) {
/* replicate extended periodic vector field fxyz
   linear interpolation
   nx/ny/nz = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   nze = third dimension of field arrays, must be >= nz+1
local data                                                            */
#define N 4
   int j, k, l, nnxye, ll;
   nnxye = N*nxe*nye;
/* copy edges of extended field */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,ll)
      for (l = 0; l < nz; l++) {
         ll = nnxye*l;
         for (k = 0; k < ny; k++) {
            fxyz[N*(nx+nxe*k)+ll] = fxyz[N*nxe*k+ll];
            fxyz[1+N*(nx+nxe*k)+ll] = fxyz[1+N*nxe*k+ll];
            fxyz[2+N*(nx+nxe*k)+ll] = fxyz[2+N*nxe*k+ll];
         }
         for (j = 0; j < nx; j++) {
            fxyz[N*(j+nxe*ny)+ll] = fxyz[N*j+ll];
            fxyz[1+N*(j+nxe*ny)+ll] = fxyz[1+N*j+ll];
            fxyz[2+N*(j+nxe*ny)+ll] = fxyz[2+N*j+ll];
         }
         fxyz[N*(nx+nxe*ny)+ll] = fxyz[ll];
         fxyz[1+N*(nx+nxe*ny)+ll] = fxyz[1+ll];
         fxyz[2+N*(nx+nxe*ny)+ll] = fxyz[2+ll];
      }
#pragma omp for \
private(j,k)
      for (k = 0; k < ny; k++) {
         for (j = 0; j < nx; j++) {
            fxyz[N*(j+nxe*k)+nnxye*nz] = fxyz[N*(j+nxe*k)];
            fxyz[1+N*(j+nxe*k)+nnxye*nz] = fxyz[1+N*(j+nxe*k)];
            fxyz[2+N*(j+nxe*k)+nnxye*nz] = fxyz[2+N*(j+nxe*k)];
         }
         fxyz[N*(nx+nxe*k)+nnxye*nz] = fxyz[N*nxe*k];
         fxyz[1+N*(nx+nxe*k)+nnxye*nz] = fxyz[1+N*nxe*k];
         fxyz[2+N*(nx+nxe*k)+nnxye*nz] = fxyz[2+N*nxe*k];
      }
   }
   for (j = 0; j < nx; j++) {
      fxyz[N*(j+nxe*ny)+nnxye*nz] = fxyz[N*j];
      fxyz[1+N*(j+nxe*ny)+nnxye*nz] = fxyz[1+N*j];
      fxyz[2+N*(j+nxe*ny)+nnxye*nz] = fxyz[2+N*j];
   }
   fxyz[N*(nx+nxe*ny)+nnxye*nz] = fxyz[0];
   fxyz[1+N*(nx+nxe*ny)+nnxye*nz] = fxyz[1];
   fxyz[2+N*(nx+nxe*ny)+nnxye*nz] = fxyz[2];
   return;
#undef N
}

/*--------------------------------------------------------------------*/
void caguard3l(float q[], int nx, int ny, int nz, int nxe, int nye,
               int nze) {
/* accumulate extended periodic scalar field q
   linear interpolation
   nx/ny/nz = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   nze = third dimension of field arrays, must be >= nz+1
local data                                                            */
   int j, k, l, nxye, ll;
   nxye = nxe*nye;
/* accumulate edges of extended field */
#pragma omp parallel
   {
#pragma omp for \
private(j,k,l,ll)
      for (l = 0; l < nz; l++) {
         ll = nxye*l;
         for (k = 0; k < ny; k++) {
            q[nxe*k+ll] += q[nx+nxe*k+ll];
            q[nx+nxe*k+ll] = 0.0;
         }
         for (j = 0; j < nx; j++) {
            q[j+ll] += q[j+nxe*ny+ll];
            q[j+nxe*ny+ll] = 0.0;
         }
         q[ll] += q[nx+nxe*ny+ll];
         q[nx+nxe*ny+ll] = 0.0;
      }
#pragma omp for \
private(j,k)
      for (k = 0; k < ny; k++) {
         for (j = 0; j < nx; j++) {
            q[j+nxe*k] += q[j+nxe*k+nxye*nz];
            q[j+nxe*k+nxye*nz] = 0.0;
         }
         q[nxe*k] += q[nx+nxe*k+nxye*nz];
         q[nx+nxe*k+nxye*nz] = 0.0;
      }
   }
   for (j = 0; j < nx; j++) {
      q[j] += q[j+nxe*ny+nxye*nz];
      q[j+nxe*ny+nxye*nz] = 0.0;
   }
   q[0] += q[nx+nxe*ny+nxye*nz];
   q[nx+nxe*ny+nxye*nz] = 0.0;
   return;
}

/*--------------------------------------------------------------------*/
void cvmpois33(float complex q[], float complex fxyz[], int isign,
               float complex ffc[], float ax, float ay, float az,
               float affp, float *we, int nx, int ny, int nz, int nxvh,
               int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
/* this subroutine solves 3d poisson's equation in fourier space for
   force/charge (or convolution of electric field over particle shape)
   with periodic boundary conditions.
   for isign = 0, output: ffc
   input: isign,ax,ay,az,affp,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
   for isign = -1, output: fxyz, we
   input: q,ffc,isign,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
   approximate flop count is:
   59*nxc*nyc*nzc + 26*(nxc*nyc + nxc*nzc + nyc*nzc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   if isign = 0, form factor array is prepared
   if isign is not equal to 0, force/charge is calculated
   equation used is:
   fx[kz][ky][kx] = -sqrt(-1)*kx*g[kz][ky][kx]*s[kz][ky][kx],
   fy[kz][ky][kx] = -sqrt(-1)*ky*g[kz][ky][kx]*s[kz][ky][kx],
   fz[kz][ky][kx] = -sqrt(-1)*kz*g[kz][ky][kx]*s[kz][ky][kx],
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
   j,k,l = fourier mode numbers,
   g[kz][ky][kx] = (affp/(kx**2+ky**2+kz**2))*s[kz][ky][kx],
   s[kz][ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
   fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
   fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
   fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
   fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
   q[l][k][j] = complex charge density for fourier mode (j,k,l)
   fxyz[l][k][j][0] = x component of complex force/charge
   fxyz[l][k][j][1] = y component of complex force/charge
   fxyz[l][k][j][2] = z component of complex force/charge
   all for fourier mode (j,k,l)
   cimag(ffc[l][k][j]) = finite-size particle shape factor s
   for fourier mode (j,k,l)
   creal(ffc[l][k][j]) = potential green's function g
   for fourier mode (j,k,l)
   ax/ay/az = half-width of particle in x/y/z direction
   affp = normalization constant = nx*ny*nz/np,
   where np=number of particles
   electric field energy is also calculated, using
   we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
      |q[kz][ky][kx]*s[kz][ky][kx]|**2)
   nx/ny/nz = system length in x/y/z direction
   nxvh = first dimension of field arrays, must be >= nxh
   nyv = second dimension of field arrays, must be >= ny
   nzv = third dimension of field arrays, must be >= nz
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
   nzhd = third dimension of form factor array, must be >= nzh
   vectorizable version
local data                                                            */
#define N 4
   int nxh, nyh, nzh, j, k, l, k1, l1, kk, kj, ll, lj, nxyhd, nxvyh;
   float dnx, dny, dnz, dkx, dky, dkz, at1, at2, at3, at4, at5, at6;
   float complex zero, zt1, zt2;
   double wp, sum1, sum2;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxyhd = nxhd*nyhd;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   zero = 0.0 + 0.0*_Complex_I;
   if (isign != 0)
      goto L40;
/* prepare form factor array */
   for (l = 0; l < nzh; l++) {
      dkz = dnz*(float) l;
      ll = nxyhd*l;
      at1 = dkz*dkz;
      at2 = pow((dkz*az),2);
      for (k = 0; k < nyh; k++) {
         dky = dny*(float) k;
         kk = nxhd*k;
         at3 = dky*dky + at1;
         at4 = pow((dky*ay),2) + at2;
         for (j = 0; j < nxh; j++) {
            dkx = dnx*(float) j;
            at5 = dkx*dkx + at3;
            at6 = exp(-0.5*(pow((dkx*ax),2) + at4));
            if (at5==0.0) {
               ffc[j+kk+ll] = affp + 1.0*_Complex_I;
            }
            else {
               ffc[j+kk+ll] = (affp*at6/at5) + at6*_Complex_I;
            }
         }
      }
   }
   return;
/* calculate force/charge and sum field energy */
L40: sum1 = 0.0;
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,k1,l1,ll,lj,kk,kj,dky,dkz,at1,at2,at3,at4,zt1,zt2,wp) \
reduction(+:sum1)
      for (l = 1; l < nzh; l++) {
         dkz = dnz*(float) l;
         ll = nxyhd*l;
         lj = nxvyh*l;
         l1 = nxvyh*nz - lj;
         wp = 0.0;
         for (k = 1; k < nyh; k++) {
            dky = dny*(float) k;
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
#pragma ivdep
            for (j = 1; j < nxh; j++) {
               at1 = crealf(ffc[j+kk+ll])*cimagf(ffc[j+kk+ll]);
               at2 = at1*dnx*(float) j;
               at3 = dky*at1;
               at4 = dkz*at1;
               zt1 = cimagf(q[j+kj+lj]) - crealf(q[j+kj+lj])*_Complex_I;
               zt2 = cimagf(q[j+k1+lj]) - crealf(q[j+k1+lj])*_Complex_I;
               fxyz[N*(j+kj+lj)] = at2*zt1;
               fxyz[1+N*(j+kj+lj)] = at3*zt1;
               fxyz[2+N*(j+kj+lj)] = at4*zt1;
               fxyz[N*(j+k1+lj)] = at2*zt2;
               fxyz[1+N*(j+k1+lj)] = -at3*zt2;
               fxyz[2+N*(j+k1+lj)] = at4*zt2;
               zt1 = cimagf(q[j+kj+l1]) - crealf(q[j+kj+l1])*_Complex_I;
               zt2 = cimagf(q[j+k1+l1]) - crealf(q[j+k1+l1])*_Complex_I;
               fxyz[N*(j+kj+l1)] = at2*zt1;
               fxyz[1+N*(j+kj+l1)] = at3*zt1;
               fxyz[2+N*(j+kj+l1)] = -at4*zt1;
               fxyz[N*(j+k1+l1)] = at2*zt2;
               fxyz[1+N*(j+k1+l1)] = -at3*zt2;
               fxyz[2+N*(j+k1+l1)] = -at4*zt2;
               at1 = at1*(q[j+kj+lj]*conjf(q[j+kj+lj])                 
                   + q[j+k1+lj]*conjf(q[j+k1+lj])
                   + q[j+kj+l1]*conjf(q[j+kj+l1])      
                   + q[j+k1+l1]*conjf(q[j+k1+l1]));
               wp += (double) at1;
            }
         }
/* mode numbers kx = 0, nx/2 */
#pragma ivdep
         for (k = 1; k < nyh; k++) {
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            at1 = crealf(ffc[kk+ll])*cimagf(ffc[kk+ll]);
            at3 = at1*dny*(float) k;
            at4 = dkz*at1;
            zt1 = cimagf(q[kj+lj]) - crealf(q[kj+lj])*_Complex_I;
            zt2 = cimagf(q[kj+l1]) - crealf(q[kj+l1])*_Complex_I;
            fxyz[N*(kj+lj)] = zero;
            fxyz[1+N*(kj+lj)] = at3*zt1;
            fxyz[2+N*(kj+lj)] = at4*zt1;
            fxyz[N*(k1+lj)] = zero;
            fxyz[1+N*(k1+lj)] = zero;
            fxyz[2+N*(k1+lj)] = zero;
            fxyz[N*(kj+l1)] = zero;
            fxyz[1+N*(kj+l1)] = at3*zt2;
            fxyz[2+N*(kj+l1)] = -at4*zt2;
            fxyz[N*(k1+l1)] = zero;
            fxyz[1+N*(k1+l1)] = zero;
            fxyz[2+N*(k1+l1)] = zero;
            at1 = at1*(q[kj+lj]*conjf(q[kj+lj])
                + q[kj+l1]*conjf(q[kj+l1]));
            wp += (double) at1;
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
#pragma ivdep
         for (j = 1; j < nxh; j++) {
            at1 = crealf(ffc[j+ll])*cimagf(ffc[j+ll]);
            at2 = at1*dnx*(float) j;  
            at4 = dkz*at1;
            zt1 = cimagf(q[j+lj]) - crealf(q[j+lj])*_Complex_I;
            zt2 = cimagf(q[j+l1]) - crealf(q[j+l1])*_Complex_I;
            fxyz[N*(j+lj)] = at2*zt1;
            fxyz[1+N*(j+lj)] = zero;
            fxyz[2+N*(j+lj)] = at4*zt1;
            fxyz[N*(j+k1+lj)] = zero;
            fxyz[1+N*(j+k1+lj)] = zero;
            fxyz[2+N*(j+k1+lj)] = zero;
            fxyz[N*(j+l1)] = at2*zt2;
            fxyz[1+N*(j+l1)] = zero;
            fxyz[2+N*(j+l1)] = -at4*zt2;
            fxyz[N*(j+k1+l1)] = zero;
            fxyz[1+N*(j+k1+l1)] = zero;
            fxyz[2+N*(j+k1+l1)] = zero;
            at1 = at1*(q[j+lj]*conjf(q[j+lj])                           
                + q[j+l1]*conjf(q[j+l1]));
            wp += (double) at1;
         }
/* mode numbers kx = 0, nx/2 */
         at1 = crealf(ffc[ll])*cimagf(ffc[ll]);
         at4 = dkz*at1;
         zt1 = cimagf(q[lj]) - crealf(q[lj])*_Complex_I;
         fxyz[N*lj] = zero;
         fxyz[1+N*lj] = zero;
         fxyz[2+N*lj] = at4*zt1;
         fxyz[N*(k1+lj)] = zero;
         fxyz[1+N*(k1+lj)] = zero;
         fxyz[2+N*(k1+lj)] = zero;
         fxyz[N*l1] = zero;
         fxyz[1+N*l1] = zero;
         fxyz[2+N*l1] = zero;
         fxyz[N*(k1+l1)] = zero;
         fxyz[1+N*(k1+l1)] = zero;
         fxyz[2+N*(k1+l1)] = zero;
         at1 = at1*(q[lj]*conjf(q[lj]));
         wp += (double) at1;
         sum1 += wp;
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   sum2 = 0.0;
#pragma omp parallel for \
private(j,k,k1,kk,kj,dky,at1,at2,at3,zt1,zt2,wp) \
reduction(+:sum2)
   for (k = 1; k < nyh; k++) {
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
         fxyz[N*(j+kj)] = at2*zt1;
         fxyz[1+N*(j+kj)] = at3*zt1;
         fxyz[2+N*(j+kj)] = zero;
         fxyz[N*(j+k1)] = at2*zt2;
         fxyz[1+N*(j+k1)] = -at3*zt2;
         fxyz[2+N*(j+k1)] = zero;
         fxyz[N*(j+kj+l1)] = zero;
         fxyz[1+N*(j+kj+l1)] = zero;
         fxyz[2+N*(j+kj+l1)] = zero;
         fxyz[N*(j+k1+l1)] = zero;
         fxyz[1+N*(j+k1+l1)] = zero;
         fxyz[2+N*(j+k1+l1)] = zero;
         at1 = at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1]));
         wp += (double) at1;
      }
/* mode numbers kx = 0, nx/2 */
      at1 = crealf(ffc[kk])*cimagf(ffc[kk]);
      at3 = at1*dny*(float) k;
      zt1 = cimagf(q[kj]) - crealf(q[kj])*_Complex_I;
      fxyz[N*kj] = zero;
      fxyz[1+N*kj] = at3*zt1;
      fxyz[2+N*kj] = zero;
      fxyz[N*k1] = zero;
      fxyz[1+N*k1] = zero;
      fxyz[2+N*k1] = zero;
      fxyz[N*(kj+l1)] = zero;
      fxyz[1+N*(kj+l1)] = zero;
      fxyz[2+N*(kj+l1)] = zero;
      fxyz[N*(k1+l1)] = zero;
      fxyz[1+N*(k1+l1)] = zero;
      fxyz[2+N*(k1+l1)] = zero;
      at1 = at1*(q[kj]*conjf(q[kj]));
      wp += (double) at1;
      sum2 += wp;
   }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
#pragma ivdep
   for (j = 1; j < nxh; j++) {
      at1 = crealf(ffc[j])*cimagf(ffc[j]);
      at2 = at1*dnx*(float) j;
      zt1 = cimagf(q[j]) - crealf(q[j])*_Complex_I;
      fxyz[N*j] = at2*zt1;
      fxyz[1+N*j] = zero;
      fxyz[2+N*j] = zero;
      fxyz[N*(j+k1)] = zero;
      fxyz[1+N*(j+k1)] = zero;
      fxyz[2+N*(j+k1)] = zero;
      fxyz[N*(j+l1)] = zero;
      fxyz[1+N*(j+l1)] = zero;
      fxyz[2+N*(j+l1)] = zero;
      fxyz[N*(j+k1+l1)] = zero;
      fxyz[1+N*(j+k1+l1)] = zero;
      fxyz[2+N*(j+k1+l1)] = zero;
      at1 = at1*(q[j]*conjf(q[j]));
      wp += (double) at1;
   }
   fxyz[0] = zero;
   fxyz[1] = zero;
   fxyz[2] = zero;
   fxyz[N*k1] = zero;
   fxyz[1+N*k1] = zero;
   fxyz[2+N*k1] = zero;
   fxyz[N*l1] = zero;
   fxyz[1+N*l1] = zero;
   fxyz[2+N*l1] = zero;
   fxyz[N*(k1+l1)] = zero;
   fxyz[1+N*(k1+l1)] = zero;
   fxyz[2+N*(k1+l1)] = zero;
   *we = (sum1 + sum2 + wp)*((float) nx)*((float) ny)*((float) nz);
   return;
#undef N
}

/*--------------------------------------------------------------------*/
void cwfft3rinit(int mixup[], float complex sct[], int indx, int indy,
                 int indz, int nxhyzd, int nxyzhd) {
/* this subroutine calculates tables needed by a three dimensional
   real to complex fast fourier transform and its inverse.
   input: indx, indy, indz, nxhyzd, nxyzhd
   output: mixup, sct
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   nxhyzd = maximum of (nx/2,ny,nz)
   nxyzhd = one half of maximum of (nx,ny,nz)
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, ny, nz, nxyz, nxhyz, nxyzh;
   int j, k, lb, ll, jb, it;
   float dnxyz, arg;
   indx1 = indx - 1;
   ndx1yz = indx1 > indy ? indx1 : indy;
   ndx1yz = ndx1yz > indz ? ndx1yz : indz;
   nx = 1L<<indx;
   ny = 1L<<indy;
   nz = 1L<<indz;
   nxyz = nx > ny ? nx : ny;
   nxyz = nxyz > nz ? nxyz : nz;
   nxhyz = 1L<<ndx1yz;
/* bit-reverse index table: mixup[j] = 1 + reversed bits of j */
   for (j = 0; j < nxhyz; j++) {
      lb = j;
      ll = 0;
      for (k = 0; k < ndx1yz; k++) {
         jb = lb/2;
         it = lb - 2*jb;
         lb = jb;
         ll = 2*ll + it;
      }
      mixup[j] = ll + 1;
   }
/* sine/cosine table for the angles 2*n*pi/nxyz */
   nxyzh = nxyz/2;
   dnxyz = 6.28318530717959/(float) nxyz;
   for (j = 0; j < nxyzh; j++) {
      arg = dnxyz*(float) j;
      sct[j] = cosf(arg) - sinf(arg)*_Complex_I;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rvmxy(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nzi, int nzp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd) {
/* this subroutine performs the x-y part of a three dimensional real to
   complex fast fourier transform and its inverse, for a subset of z,
   using complex arithmetic, with Vector/OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, an inverse fourier transform in x and y is performed
   f[i][m][n] = (1/nx*ny*nz)*sum(f[i][k][j]*exp(-sqrt(-1)*2pi*n*j/nx)*
         exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, a forward fourier transform in x and y is performed
   f[l][k][j] = sum(f[l][m][n]*exp(sqrt(-1)*2pi*n*j/nx)*
         exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nzi = initial z index used
   nzp = number of z indices used
   nxhd = first dimension of f
   nyd,nzd = second and third dimensions of f
   nxhyzd = maximum of (nx/2,ny,nz)
   nxyzhd = maximum of (nx,ny,nz)/2
   fourier coefficients are stored as follows:
   f[l][k][j] = real, imaginary part of mode j,k,l
   where 0 <= j < nx/2, 0 <= k < ny, 0 <= l < nz, except for
   f[l][k][0] = real, imaginary part of mode nx/2,k,l,
   where ny/2+1 <= k < ny and 0 <= l < nz, and
   f[l][0][0] = real, imaginary part of mode nx/2,0,l,
   f[l][ny/2][0] = real, imaginary part mode nx/2,ny/2,l,
   where nz/2+1 <= l < nz, and
   imag(f[0][0][0]) = real part of mode nx/2,0,0
   imag(f[0][ny/2][0]) = real part of mode nx/2,ny/2,0
   imag(f[nz/2][0][0]) = real part of mode nx/2,0,nz/2
   imag(f[nz/2][ny/2][0]) = real part of mode nx/2,ny/2,nz/2
   using jpl storage convention, as described in:
   E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
   Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
   Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
   December 1993.
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, nxhh, ny, nyh;
   int nz, nxyz, nxhyz, nzt, nrx, nry, nrxb, nryb, nxhyd;
   int i, j, k, l, n, nn, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
   float ani;
   float complex t1, t2, t3;
   if (isign==0)
      return;
   indx1 = indx - 1;
   ndx1yz = indx1 > indy ? indx1 : indy;
   ndx1yz = ndx1yz > indz ? ndx1yz : indz;
   nx = 1L<<indx;
   nxh = nx/2;
   nxhh = nx/4;
   ny = 1L<<indy;
   nyh = ny/2;
   nz = 1L<<indz;
   nxyz = nx > ny ? nx : ny;
   nxyz = nxyz > nz ? nxyz : nz;
   nxhyz = 1L<<ndx1yz;
   nzt = nzi + nzp - 1;
   nxhyd = nxhd*nyd;
   if (isign > 0)
      goto L180;
/* inverse fourier transform */
   nrxb = nxhyz/nxh;
   nrx = nxyz/nxh;
   nryb = nxhyz/ny;
   nry = nxyz/ny;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,nn,joff,ani,t1,t2,t3)
   for (n = nzi-1; n < nzt; n++) {
      nn = nxhyd*n;
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            for (i = 0; i < ny; i++) {
               joff = nxhd*i + nn;
               t1 = f[j1+joff];
               f[j1+joff] = f[j+joff];
               f[j+joff] = t1;
            }
         }
      }
/* first transform in x */
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
               for (i = 0; i < ny; i++) {
                  joff = nxhd*i + nn;
                  t2 = t1*f[j+k2+joff];
                  f[j+k2+joff] = f[j+k1+joff] - t2;
                  f[j+k1+joff] += t2;
               }
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxyz/nx;
      ani = 0.5/(((float) nx)*((float) ny)*((float) nz));
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < ny; k++) {
            joff = nxhd*k + nn;
            t2 = conjf(f[nxh-j+joff]);
            t1 = f[j+joff] + t2;
            t2 = (f[j+joff] - t2)*t3;
            f[j+joff] = ani*(t1 + t2);
            f[nxh-j+joff] = ani*conjf(t1 - t2);
         }
      }
      ani = 2.0*ani;
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
         f[nxhh+joff] = ani*conjf(f[nxhh+joff]);
         f[joff] = ani*((crealf(f[joff]) + cimagf(f[joff]))
                   + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I);
      }
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nxhd*k1 + nn;
            for (i = 0; i < nxh; i++) {
               t1 = f[i+k1];
               f[i+k1] = f[i+joff];
               f[i+joff] = t1;
            }
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
               j1 = nxhd*(j + k1) + nn;
               j2 = nxhd*(j + k2) + nn;
               t1 = sct[kmr*j];
               for (i = 0; i < nxh; i++) {
                  t2 = t1*f[i+j2];
                  f[i+j2] = f[i+j1] - t2;
                  f[i+j1] += t2;
               }
            }
         }
         ns = ns2;
      }
/* unscramble modes kx = 0, nx/2 */
#pragma ivdep
      for (k = 1; k < nyh; k++) {
         joff = nxhd*k;
         k1 = nxhd*ny - joff + nn;
         joff += nn;
         t1 = f[k1];
         f[k1] = 0.5*(cimagf(f[joff] + t1)
                  + crealf(f[joff] - t1)*_Complex_I);
         f[joff] = 0.5*(crealf(f[joff] + t1)
                    + cimagf(f[joff] - t1)*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
L180: nryb = nxhyz/ny;
   nry = nxyz/ny;
   nrxb = nxhyz/nxh;
   nrx = nxyz/nxh;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,nn,joff,t1,t2,t3)
   for (n = nzi-1; n < nzt; n++) {
      nn = nxhyd*n;
/* scramble modes kx = 0, nx/2 */
#pragma ivdep
      for (k = 1; k < nyh; k++) {
         joff = nxhd*k;
         k1 = nxhd*ny - joff + nn;
         joff += nn;
         t1 = cimagf(f[k1]) + crealf(f[k1])*_Complex_I;
         f[k1] = conjf(f[joff] - t1);
         f[joff] += t1;
      }
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nxhd*k1 + nn;
            for (i = 0; i < nxh; i++) {
               t1 = f[i+k1];
               f[i+k1] = f[i+joff];
               f[i+joff] = t1;
            }
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
               j1 = nxhd*(j + k1) + nn;
               j2 = nxhd*(j + k2) + nn;
               t1 = conjf(sct[kmr*j]);
               for (i = 0; i < nxh; i++) {
                  t2 = t1*f[i+j2];
                  f[i+j2] = f[i+j1] - t2;
                  f[i+j1] += t2;
               }
            }
         }
         ns = ns2;
      }
/* scramble coefficients */
      kmr = nxyz/nx;
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < ny; k++) {
            joff = nxhd*k + nn;
            t2 = conjf(f[nxh-j+joff]);
            t1 = f[j+joff] + t2;
            t2 = (f[j+joff] - t2)*t3;
            f[j+joff] = t1 + t2;
            f[nxh-j+joff] = conjf(t1 - t2);
         }
      }
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
         f[nxhh+joff] = 2.0*conjf(f[nxhh+joff]);
         f[joff] = (crealf(f[joff]) + cimagf(f[joff]))
                   + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I;
      }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            for (i = 0; i < ny; i++) {
               joff = nxhd*i + nn;
               t1 = f[j1+joff];
               f[j1+joff] = f[j+joff];
               f[j+joff] = t1;
            }
         }
      }
/* finally transform in x */
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
               for (i = 0; i < ny; i++) {
                  joff = nxhd*i + nn;
                  t2 = t1*f[j+k2+joff];
                  f[j+k2+joff] = f[j+k1+joff] - t2;
                  f[j+k1+joff] += t2;
               }
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rvmxz(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nyi, int nyp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd) {
/* this subroutine performs the z part of a three dimensional real to
   complex fast fourier transform and its inverse, for a subset of y,
   using complex arithmetic, with Vector/OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, an inverse fourier transform in z is performed
   f[l][k][j] = sum(f[i][k][j]*exp(-sqrt(-1)*2pi*l*i/nz))
   if isign = 1, a forward fourier transform in z is performed
   f[i][m][n] = sum(f[l][m][n]*exp(sqrt(-1)*2pi*l*i/nz))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = first dimension of f
   nyd,nzd = second and third dimensions of f
   nxhyzd = maximum of (nx/2,ny,nz)
   nxyzhd = maximum of (nx,ny,nz)/2
   fourier coefficients are stored as follows:
   f[l][k][j] = real, imaginary part of mode j,k,l
   where 0 <= j < nx/2, 0 <= k < ny, 0 <= l < nz, except for
   f[l][k][0] = real, imaginary part of mode nx/2,k,l,
   where ny/2+1 <= k < ny and 0 <= l < nz, and
   f[l][0][0] = real, imaginary part of mode nx/2,0,l,
   f[l][ny/2][0] = real, imaginary part mode nx/2,ny/2,l,
   where nz/2+1 <= l < nz, and
   imag(f[0][0][0]) = real part of mode nx/2,0,0
   imag(f[0][ny/2][0]) = real part of mode nx/2,ny/2,0
   imag(f[nz/2][0][0]) = real part of mode nx/2,0,nz/2
   imag(f[nz/2][ny/2][0]) = real part of mode nx/2,ny/2,nz/2
   using jpl storage convention, as described in:
   E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
   Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
   Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
   December 1993.
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, ny, nyh;
   int nz, nzh, nxyz, nxhyz, nyt, nrz, nrzb, nxhyd, ioff;
   int i, j, k, l, n, ll, j1, j2, k1, k2, l1, ns, ns2, km, kmr, i0, i1;
   float complex t1, t2;
   if (isign==0)
      return;
   indx1 = indx - 1;
   ndx1yz = indx1 > indy ? indx1 : indy;
   ndx1yz = ndx1yz > indz ? ndx1yz : indz;
   nx = 1L<<indx;
   nxh = nx/2;
   ny = 1L<<indy;
   nyh = ny/2;
   nz = 1L<<indz;
   nzh = nz/2;
   nxyz = nx > ny ? nx : ny;
   nxyz = nxyz > nz ? nxyz : nz;
   nxhyz = 1L<<ndx1yz;
   nyt = nyi + nyp - 1;
   nxhyd = nxhd*nyd;
   if (isign > 0)
      goto L100;
/* inverse fourier transform */
   nrzb = nxhyz/nz;
   nrz = nxyz/nz;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,ll,l1,i0,i1,ioff,t1,t2)
   for (n = nyi-1; n < nyt; n++) {
      ioff = nxhd*n;
/* bit-reverse array elements in z */
      for (l = 0; l < nz; l++) {
         ll = nxhyd*l;
         l1 = (mixup[l] - 1)/nrzb;
         if (l < l1) {
            l1 = nxhyd*l1;
            i0 = ioff + ll;
            i1 = ioff + l1;
            for (i = 0; i < nxh; i++) {
               t1 = f[i+i1];
               f[i+i1] = f[i+i0];
               f[i+i0] = t1;
            }
         }
      }
/* finally transform in z */
      ns = 1;
      for (l = 0; l < indz; l++) {
         ns2 = ns + ns;
         km = nzh/ns;
         kmr = km*nrz;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = nxhyd*(j + k1);
               j2 = nxhyd*(j + k2);
               t1 = sct[kmr*j];
               i0 = ioff + j1;
               i1 = ioff + j2;
               for (i = 0; i < nxh; i++) {
                  t2 = t1*f[i+i1];
                  f[i+i1] = f[i+i0] - t2;
                  f[i+i0] += t2;
               }
            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if (nyi==1) {
#pragma ivdep
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         t1 = f[l1];
         f[l1] = 0.5*(cimagf(f[ll] + t1)
                    + crealf(f[ll] - t1)*_Complex_I);
         f[ll] = 0.5*(crealf(f[ll] + t1)
                    + cimagf(f[ll] - t1)*_Complex_I);
      }
   }
   if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
#pragma ivdep
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         i1 = nxhd*nyh;
         i0 = i1 + ll;
         i1 += l1;
         t1 = f[i1];
         f[i1] = 0.5*(cimagf(f[i0] + t1)
                  +   crealf(f[i0] - t1)*_Complex_I);
         f[i0] = 0.5*(crealf(f[i0] + t1)
                    + cimagf(f[i0] - t1)*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
L100: nrzb = nxhyz/nz;
   nrz = nxyz/nz;
/* scramble modes kx = 0, nx/2 */
   if (nyi==1) {
#pragma ivdep
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         t1 = cimagf(f[l1]) + crealf(f[l1])*_Complex_I;
         f[l1] = conjf(f[ll] - t1);
         f[ll] += t1;
      }
   }
   if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
#pragma ivdep
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         i1 = nxhd*nyh;
         i0 = i1 + ll;
         i1 += l1;
         t1 = cimagf(f[i1]) + crealf(f[i1])*_Complex_I;
         f[i1] = conjf(f[i0] - t1);
         f[i0] += t1;
      }
   }
/* bit-reverse array elements in z */
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,ll,l1,i0,i1,ioff,t1,t2)
   for (n = nyi-1; n < nyt; n++) {
      ioff = nxhd*n;
      for (l = 0; l < nz; l++) {
         ll = nxhyd*l;
         l1 = (mixup[l] - 1)/nrzb;
         if (l < l1) {
            l1 = nxhyd*l1;
            i0 = ioff + ll;
            i1 = ioff + l1;
            for (i = 0; i < nxh; i++) {
               t1 = f[i+i1];
               f[i+i1] = f[i+i0];
               f[i+i0] = t1;
            }
         }
      }
/* first transform in z */
      ns = 1;
      for (l = 0; l < indz; l++) {
         ns2 = ns + ns;
         km = nzh/ns;
         kmr = km*nrz;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = nxhyd*(j + k1);
               j2 = nxhyd*(j + k2);
               t1 = conjf(sct[kmr*j]);
               i0 = ioff + j1;
               i1 = ioff + j2;
               for (i = 0; i < nxh; i++) {
                  t2 = t1*f[i+i1];
                  f[i+i1] = f[i+i0] - t2;
                  f[i+i0] += t2;
               }
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rvm3xy(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int nzi, int nzp, int nxhd, int nyd, int nzd,
                 int nxhyzd, int nxyzhd) {
/* this subroutine performs the x-y part of 3 three dimensional complex
   to real fast fourier transforms and their inverses, for a subset of z,
   using complex arithmetic, with Vector/OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, three inverse fourier transforms in x and y are
   performed
   f[i][m][n][0:2] = (1/nx*ny*nz)*sum(f[i][k][j][0:2]*exp(-sqrt(-1)*2pi*n*j/nx)
         *exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, three forward fourier transforms in x and y are
   performed
   f[l][k][j][0:2] = sum(f[l][m][n][0:2]*exp(sqrt(-1)*2pi*n*j/nx)*
         exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nzi = initial z index used
   nzp = number of z indices used
   nxhd = second dimension of f
   nyd,nzd = third and fourth dimensions of f
   nxhyzd = maximum of (nx/2,ny,nz)
   nxyzhd = maximum of (nx,ny,nz)/2
   fourier coefficients are stored as follows:
   f[l][k][j][0:2] = real, imaginary part of mode j,k,l
   where 0 <= j < nx/2, 0 <= k < ny, 0 <= l < nz, except for
   f[l][k][0][0:2] = real, imaginary part of mode nx/2,k,l,
   where ny/2+1 <= k < ny and 0 <= l < nz, and
   f[l][0][0][0:2] = real, imaginary part of mode nx/2,0,l,
   f[l][ny/2][0][0:2] = real, imaginary part mode nx/2,ny/2,l,
   where nz/2+1 <= l < nz, and
   imag(f[0][0][0][0:2]) = real part of mode nx/2,0,0
   imag(f[0][ny/2][0][0:2]) = real part of mode nx/2,ny/2,0
   imag(f[nz/2][0][0][0:2]) = real part of mode nx/2,0,nz/2
   imag(f[nz/2][ny/2][0][0:2]) = real part of mode nx/2,ny/2,nz/2
   using jpl storage convention, as described in:
   E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
   Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
   Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
   December 1993.
   written by viktor k. decyk, ucla
local data                                                            */
#define N 4
   int indx1, ndx1yz, nx, nxh, nxhh, ny, nyh;
   int nz, nxyz, nxhyz, nzt, nrx, nry, nrxb, nryb, nnxhd, nxhyd;
   int i, j, k, l, n, nn, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
   float at1, at2, ani;
   float complex t1, t2, t3, t4;
   if (isign==0)
      return;
   indx1 = indx - 1;
   ndx1yz = indx1 > indy ? indx1 : indy;
   ndx1yz = ndx1yz > indz ? ndx1yz : indz;
   nx = 1L<<indx;
   nxh = nx/2;
   nxhh = nx/4;
   ny = 1L<<indy;
   nyh = ny/2;
   nz = 1L<<indz;
   nxyz = nx > ny ? nx : ny;
   nxyz = nxyz > nz ? nxyz : nz;
   nxhyz = 1L<<ndx1yz;
   nzt = nzi + nzp - 1;
   nnxhd = N*nxhd;
   nxhyd = nnxhd*nyd;
   if (isign > 0)
      goto L230;
/* inverse fourier transform */
   nrxb = nxhyz/nxh;
   nrx = nxyz/nxh;
   nryb = nxhyz/ny;
   nry = nxyz/ny;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,jj,j1,j2,nn,joff,at1,at2,ani,t1, \
t2,t3,t4)
   for (n = nzi-1; n < nzt; n++) {
      nn = nxhyd*n;
/* swap complex components */
      for (i = 0; i < ny; i++) {
         joff = nnxhd*i + nn;
         for (j = 0; j < nxh; j++) {
            at1 = cimagf(f[2+N*j+joff]);
            at2 = crealf(f[2+N*j+joff]);
            f[2+N*j+joff] = crealf(f[1+N*j+joff])
                            + crealf(f[3+N*j+joff])*_Complex_I;
            f[1+N*j+joff] = cimagf(f[N*j+joff]) + at1*_Complex_I;
            f[N*j+joff] = crealf(f[N*j+joff]) + at2*_Complex_I;
         }
      }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            for (i = 0; i < ny; i++) {
               joff = nnxhd*i + nn;
               t1 = f[N*j1+joff];
               t2 = f[1+N*j1+joff];
               t3 = f[2+N*j1+joff];
               f[N*j1+joff] = f[N*j+joff];
               f[1+N*j1+joff] = f[1+N*j+joff];
               f[2+N*j1+joff] = f[2+N*j+joff];
               f[N*j+joff] = t1;
               f[1+N*j+joff] = t2;
               f[2+N*j+joff] = t3;
            }
         }
      }
/* first transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = N*ns2*k;
            k2 = k1 + 4*ns;
            for (j = 0; j < ns; j++) {
               t1 = sct[kmr*j];
               for (i = 0; i < ny; i++) {
                  joff = nnxhd*i + nn;
                  t2 = t1*f[N*j+k2+joff];
                  t3 = t1*f[1+N*j+k2+joff];
                  t4 = t1*f[2+N*j+k2+joff];
                  f[N*j+k2+joff] = f[N*j+k1+joff] - t2;
                  f[1+N*j+k2+joff] = f[1+N*j+k1+joff] - t3;
                  f[2+N*j+k2+joff] = f[2+N*j+k1+joff] - t4;
                  f[N*j+k1+joff] += t2;
                  f[1+N*j+k1+joff] += t3;
                  f[2+N*j+k1+joff] += t4;
               }
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxyz/nx;
      ani = 0.5/(((float) nx)*((float) ny)*((float) nz));
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < ny; k++) {
            joff = nnxhd*k + nn;
            for (jj = 0; jj < 3; jj++) {
               t2 = conjf(f[jj+N*(nxh-j)+joff]);
               t1 = f[jj+N*j+joff] + t2;
               t2 = (f[jj+N*j+joff] - t2)*t3;
               f[jj+N*j+joff] = ani*(t1 + t2);
               f[jj+N*(nxh-j)+joff] = ani*conjf(t1 - t2);
            }
         }
      }
      ani = 2.0*ani;
      for (k = 0; k < ny; k++) {
         joff = nnxhd*k + nn;
         for (jj = 0; jj < 3; jj++) {
            f[jj+N*nxhh+joff] = ani*conjf(f[jj+N*nxhh+joff]);
            f[jj+joff] = ani*((crealf(f[jj+joff]) 
                          + cimagf(f[jj+joff]))
                          + (crealf(f[jj+joff])
                          - cimagf(f[jj+joff]))*_Complex_I);
         }
      }
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         joff = nnxhd*k + nn;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nnxhd*k1 + nn;
            for (i = 0; i < nxh; i++) {
               t1 = f[N*i+k1];
               t2 = f[1+N*i+k1];
               t3 = f[2+N*i+k1];
               f[N*i+k1] = f[N*i+joff];
               f[1+N*i+k1] = f[1+N*i+joff];
               f[2+N*i+k1] = f[2+N*i+joff];
               f[N*i+joff] = t1;
               f[1+N*i+joff] = t2;
               f[2+N*i+joff] = t3;
            }
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
               j1 = nnxhd*(j + k1) + nn;
               j2 = nnxhd*(j + k2) + nn;
               t1 = sct[kmr*j];
               for (i = 0; i < nxh; i++) {
                  t2 = t1*f[N*i+j2];
                  t3 = t1*f[1+N*i+j2];
                  t4 = t1*f[2+N*i+j2];
                  f[N*i+j2] = f[N*i+j1] - t2;
                  f[1+N*i+j2] = f[1+N*i+j1] - t3;
                  f[2+N*i+j2] = f[2+N*i+j1] - t4;
                  f[N*i+j1] += t2;
                  f[1+N*i+j1] += t3;
                  f[2+N*i+j1] += t4;
               }
            }
         }
         ns = ns2;
      }
/* unscramble modes kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         joff = nnxhd*k;
         k1 = nnxhd*ny - joff + nn;
         joff += nn;
         for (jj = 0; jj < 3; jj++) {
            t1 = f[jj+k1];
            f[jj+k1] = 0.5*(cimagf(f[jj+joff] + t1)
                        + crealf(f[jj+joff] - t1)*_Complex_I);
            f[jj+joff] = 0.5*(crealf(f[jj+joff] + t1)
                          + cimagf(f[jj+joff] - t1)*_Complex_I);
         }
      }
   }
   return;
/* forward fourier transform */
L230: nryb = nxhyz/ny;
   nry = nxyz/ny;
   nrxb = nxhyz/nxh;
   nrx = nxyz/nxh;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,jj,j1,j2,nn,joff,at1,at2,t1,t2, \
t3,t4)
   for (n = nzi-1; n < nzt; n++) {
      nn = nxhyd*n;
/* scramble modes kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         joff = nnxhd*k;
         k1 = nnxhd*ny - joff + nn;
         joff += nn;
         for (jj = 0; jj < 3; jj++) {
            t1 = cimagf(f[jj+k1]) + crealf(f[jj+k1])*_Complex_I;
            f[jj+k1] = conjf(f[jj+joff] - t1);
            f[jj+joff] += t1;
         }
      }
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         joff = nnxhd*k + nn;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nnxhd*k1 + nn;
            for (i = 0; i < nxh; i++) {
               t1 = f[N*i+k1];
               t2 = f[1+N*i+k1];
               t3 = f[2+N*i+k1];
               f[N*i+k1] = f[N*i+joff];
               f[1+N*i+k1] = f[1+N*i+joff];
               f[2+N*i+k1] = f[2+N*i+joff];
               f[N*i+joff] = t1;
               f[1+N*i+joff] = t2;
               f[2+N*i+joff] = t3;
            }
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
               j1 = nnxhd*(j + k1) + nn;
               j2 = nnxhd*(j + k2) + nn;
               t1 = conjf(sct[kmr*j]);
               for (i = 0; i < nxh; i++) {
                  t2 = t1*f[N*i+j2];
                  t3 = t1*f[1+N*i+j2];
                  t4 = t1*f[2+N*i+j2];
                  f[N*i+j2] = f[N*i+j1] - t2;
                  f[1+N*i+j2] = f[1+N*i+j1] - t3;
                  f[2+N*i+j2] = f[2+N*i+j1] - t4;
                  f[N*i+j1] += t2;
                  f[1+N*i+j1] += t3;
                  f[2+N*i+j1] += t4;
               }
            }
         }
         ns = ns2;
      }
/* scramble coefficients */
      kmr = nxyz/nx;
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < ny; k++) {
            joff = nnxhd*k + nn;
            for (jj = 0; jj < 3; jj++) {
               t2 = conjf(f[jj+N*(nxh-j)+joff]);
               t1 = f[jj+N*j+joff] + t2;
               t2 = (f[jj+N*j+joff] - t2)*t3;
               f[jj+N*j+joff] = t1 + t2;
               f[jj+N*(nxh-j)+joff] = conjf(t1 - t2);
            }
         }
      }
      for (k = 0; k < ny; k++) {
         joff = nnxhd*k + nn;
         for (jj = 0; jj < 3; jj++) {
            f[jj+N*nxhh+joff] = 2.0*conjf(f[jj+N*nxhh+joff]);
            f[jj+joff] = (crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                       + (crealf(f[jj+joff])
                       - cimagf(f[jj+joff]))*_Complex_I;
         }
      }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            for (i = 0; i < ny; i++) {
               joff = nnxhd*i + nn;
               t1 = f[N*j1+joff];
               t2 = f[1+N*j1+joff];
               t3 = f[2+N*j1+joff];
               f[N*j1+joff] = f[N*j+joff];
               f[1+N*j1+joff] = f[1+N*j+joff];
               f[2+N*j1+joff] = f[2+N*j+joff];
               f[N*j+joff] = t1;
               f[1+N*j+joff] = t2;
               f[2+N*j+joff] = t3;
            }
         }
      }
/* finally transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = N*ns2*k;
            k2 = k1 + 4*ns;
            for (j = 0; j < ns; j++) {
               t1 = conjf(sct[kmr*j]);
               for (i = 0; i < ny; i++) {
                  joff = nnxhd*i + nn;
                  t2 = t1*f[N*j+k2+joff];
                  t3 = t1*f[1+N*j+k2+joff];
                  t4 = t1*f[2+N*j+k2+joff];
                  f[N*j+k2+joff] = f[N*j+k1+joff] - t2;
                  f[1+N*j+k2+joff] = f[1+N*j+k1+joff] - t3;
                  f[2+N*j+k2+joff] = f[2+N*j+k1+joff] - t4;
                  f[N*j+k1+joff] += t2;
                  f[1+N*j+k1+joff] += t3;
                  f[2+N*j+k1+joff] += t4;
               }
            }
         }
         ns = ns2;
      }
/* swap complex components */
      for (i = 0; i < ny; i++) {
         joff = nnxhd*i + nn;
         for (j = 0; j < nxh; j++) {
            f[3+N*j+joff] = cimagf(f[2+N*j+joff])
                            + cimagf(f[3+N*j+joff])*_Complex_I;
            at1 = crealf(f[2+N*j+joff]);
            f[2+N*j+joff] = cimagf(f[N*j+joff])
                            + cimagf(f[1+N*j+joff])*_Complex_I;
            at2 = crealf(f[1+N*j+joff]);
            f[1+N*j+joff] = at1 + 0.0*_Complex_I;
            f[N*j+joff] = crealf(f[N*j+joff]) + at2*_Complex_I;
         }
      }
   }
   return;
#undef N
}

/*--------------------------------------------------------------------*/
void cfft3rvm3z(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nyi, int nyp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd) {
/* this subroutine performs the z part of 3 three dimensional complex to
   real fast fourier transforms and their inverses, for a subset of y,
   using complex arithmetic, with Vector/OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, three inverse fourier transforms in z are performed
   f[l][k][j][0:2] = sum(f[i][k][j][0:2]*exp(-sqrt(-1)*2pi*l*i/nz))
   if isign = 1, three forward fourier transforms in z are performed
   f[i][m][n][0:2] = sum(f[l][m][n][0:2]*exp(sqrt(-1)*2pi*l*i/nz))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = second dimension of f
   nyd,nzd = third and fourth dimensions of f
   nxhyzd = maximum of (nx/2,ny,nz)
   nxyzhd = maximum of (nx,ny,nz)/2
   fourier coefficients are stored as follows:
   f[l][k][j][0:2] = real, imaginary part of mode j,k,l
   where 0 <= j < nx/2, 0 <= k < ny, 0 <= l < nz, except for
   f[l][k][0][0:2], = real, imaginary part of mode nx/2,k,l,
   where ny/2+1 <= k < ny and 0 <= l < nz, and
   f[l][0][0][0:2] = real, imaginary part of mode nx/2,0,l,
   f[l][ny/2][0][0:2] = real, imaginary part mode nx/2,ny/2,l,
   where nz/2+1 <= l < nz, and
   imag(f[0][0][0][0:2]) = real part of mode nx/2,0,0
   imag(f[0][ny/2][0][0:2]) = real part of mode nx/2,ny/2,0
   imag(f[nz/2][0][0][0:2]) = real part of mode nx/2,0,nz/2
   imag(f[nz/2][ny/2][0][0:2]) = real part of mode nx/2,ny/2,nz/2
   using jpl storage convention, as described in:
   E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
   Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
   Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
   December 1993.
   written by viktor k. decyk, ucla
local data                                                            */
#define N 4
   int indx1, ndx1yz, nx, nxh, ny, nyh;
   int nz, nzh, nxyz, nxhyz, nyt, nrz, nrzb, nnxhd, nxhyd, ioff;
   int i, j, k, l, n, ll, jj, j1, j2, k1, k2, l1, ns, ns2, km, kmr;
   int i0, i1;
   float complex t1, t2, t3, t4;
   if (isign==0)
      return;
   indx1 = indx - 1;
   ndx1yz = indx1 > indy ? indx1 : indy;
   ndx1yz = ndx1yz > indz ? ndx1yz : indz;
   nx = 1L<<indx;
   nxh = nx/2;
   ny = 1L<<indy;
   nyh = ny/2;
   nz = 1L<<indz;
   nzh = nz/2;
   nxyz = nx > ny ? nx : ny;
   nxyz = nxyz > nz ? nxyz : nz;
   nxhyz = 1L<<ndx1yz;
   nyt = nyi + nyp - 1;
   nnxhd = N*nxhd;
   nxhyd = nnxhd*nyd;
   if (isign > 0)
      goto L110;
/* inverse fourier transform */
   nrzb = nxhyz/nz;
   nrz = nxyz/nz;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,ll,l1,i0,i1,ioff,t1,t2,t3, \
t4)
   for (n = nyi-1; n < nyt; n++) {
      ioff = nnxhd*n;
      for (l = 0; l < nz; l++) {
         ll = nxhyd*l;
         l1 = (mixup[l] - 1)/nrzb;
         if (l < l1) {
            l1 = nxhyd*l1;
            i0 = ioff + ll;
            i1 = ioff + l1;
            for (i = 0; i < nxh; i++) {
               t1 = f[N*i+i1];
               t2 = f[1+N*i+i1];
               t3 = f[2+N*i+i1];
               f[N*i+i1] = f[N*i+i0];
               f[1+N*i+i1] = f[1+N*i+i0];
               f[2+N*i+i1] = f[2+N*i+i0];
               f[N*i+i0] = t1;
               f[1+N*i+i0] = t2;
               f[2+N*i+i0] = t3;
            }
         }
      }
/* finally transform in z */
      ns = 1;
      for (l = 0; l < indz; l++) {
         ns2 = ns + ns;
         km = nzh/ns;
         kmr = km*nrz;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = nxhyd*(j + k1);
               j2 = nxhyd*(j + k2);
               t1 = sct[kmr*j];
               i0 = ioff + j1;
               i1 = ioff + j2;
               for (i = 0; i < nxh; i++) {
                  t2 = t1*f[N*i+i1];
                  t3 = t1*f[1+N*i+i1];
                  t4 = t1*f[2+N*i+i1];
                  f[N*i+i1] = f[N*i+i0] - t2;
                  f[1+N*i+i1] = f[1+N*i+i0] - t3;
                  f[2+N*i+i1] = f[2+N*i+i0] - t4;
                  f[N*i+i0] += t2;
                  f[1+N*i+i0] += t3;
                  f[2+N*i+i0] += t4;
               }
            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if (nyi==1) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         for (jj = 0; jj < 3; jj++) {
            t1 = f[jj+l1];
            f[jj+l1] = 0.5*(cimagf(f[jj+ll] + t1)
                          + crealf(f[jj+ll] - t1)*_Complex_I);
            f[jj+ll] = 0.5*(crealf(f[jj+ll] + t1)
                          + cimagf(f[jj+ll] - t1)*_Complex_I);
         }
      }
   }
   if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         for (jj = 0; jj < 3; jj++) {
            i1 = nnxhd*nyh;
            i0 = i1 + ll;
            i1 += l1;
            t1 = f[jj+i1];
            f[jj+i1] = 0.5*(cimagf(f[jj+i0] + t1)
                        +   crealf(f[jj+i0] - t1)*_Complex_I);
            f[jj+i0] = 0.5*(crealf(f[jj+i0] + t1)
                          + cimagf(f[jj+i0] - t1)*_Complex_I);
         }
      }
   }
   return;
/* forward fourier transform */
L110: nrzb = nxhyz/nz;
   nrz = nxyz/nz;
/* scramble modes kx = 0, nx/2 */
   if (nyi==1) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         for (jj = 0; jj < 3; jj++) {
            t1 = cimagf(f[jj+l1]) + crealf(f[jj+l1])*_Complex_I;
            f[jj+l1] = conjf(f[jj+ll] - t1);
            f[jj+ll] += t1;
         }
      }
   }
   if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         for (jj = 0; jj < 3; jj++) {
            i1 = nnxhd*nyh;
            i0 = i1 + ll;
            i1 += l1;
            t1 = cimagf(f[jj+i1]) + crealf(f[jj+i1])*_Complex_I;
            f[jj+i1] = conjf(f[jj+i0] - t1);
            f[jj+i0] += t1;
         }
      }
   }
/* bit-reverse array elements in z */
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,ll,l1,i0,i1,ioff,t1,t2,t3, \
t4)
   for (n = nyi-1; n < nyt; n++) {
      ioff = nnxhd*n;
      for (l = 0; l < nz; l++) {
         ll = nxhyd*l;
         l1 = (mixup[l] - 1)/nrzb;
         if (l < l1) {
            l1 = nxhyd*l1;
            i0 = ioff + ll;
            i1 = ioff + l1;
            for (i = 0; i < nxh; i++) {
               t1 = f[N*i+i1];
               t2 = f[1+N*i+i1];
               t3 = f[2+N*i+i1];
               f[N*i+i1] = f[N*i+i0];
               f[1+N*i+i1] = f[1+N*i+i0];
               f[2+N*i+i1] = f[2+N*i+i0];
               f[N*i+i0] = t1;
               f[1+N*i+i0] = t2;
               f[2+N*i+i0] = t3;
            }
         }
      }
/* first transform in z */
      ns = 1;
      for (l = 0; l < indz; l++) {
         ns2 = ns + ns;
         km = nzh/ns;
         kmr = km*nrz;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = nxhyd*(j + k1);
               j2 = nxhyd*(j + k2);
               t1 = conjf(sct[kmr*j]);
               i0 = ioff + j1;
               i1 = ioff + j2;
               for (i = 0; i < nxh; i++) {
                  t2 = t1*f[N*i+i1];
                  t3 = t1*f[1+N*i+i1];
                  t4 = t1*f[2+N*i+i1];
                  f[N*i+i1] = f[N*i+i0] - t2;
                  f[1+N*i+i1] = f[1+N*i+i0] - t3;
                  f[2+N*i+i1] = f[2+N*i+i0] - t4;
                  f[N*i+i0] += t2;
                  f[1+N*i+i0] += t3;
                  f[2+N*i+i0] += t4;
               }
            }
         }
         ns = ns2;
      }
   }
   return;
#undef N
}

/*--------------------------------------------------------------------*/
void cwfft3rvmx(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
/* wrapper function for real to complex fft, with packed data */
/* parallelized with OpenMP*/
/* local data */
   int ny, nz;
   static int nyi = 1, nzi = 1;
/* calculate range of indices */
   ny = 1L<<indy;
   nz = 1L<<indz;
/* inverse fourier transform */
   if (isign < 0) {
/* perform xy fft */
      cfft3rvmxy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                 nxhyzd,nxyzhd);
/* perform z fft */
      cfft3rvmxz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                 nxhyzd,nxyzhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform z fft */
      cfft3rvmxz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                 nxhyzd,nxyzhd);
/* perform xy fft */
      cfft3rvmxy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                 nxhyzd,nxyzhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rvm3(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
/* wrapper function for 3 3d real to complex ffts, with packed data */
/* parallelized with OpenMP */
/* local data */
   int ny, nz;
   static int nyi = 1, nzi = 1;
/* calculate range of indices */
   ny = 1L<<indy;
   nz = 1L<<indz;
/* inverse fourier transform */
   if (isign < 0) {
/* perform xy fft */
      cfft3rvm3xy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                 nxhyzd,nxyzhd);
/* perform z fft */
      cfft3rvm3z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                nxhyzd,nxyzhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform z fft */
      cfft3rvm3z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                nxhyzd,nxyzhd);
/* perform xy fft */
      cfft3rvm3xy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                 nxhyzd,nxyzhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cset_szero3(float q[], int mx, int my, int mz, int nxv, int nyv,
                 int nzv, int mx1, int my1, int mxyz1) {
/* for 3d code, this subroutine zeros out charge density array.
   for Intel NUMA architecture with first touch policy, this associates
   array segments with appropriate threads
   OpenMP version
   input: all, output: q
   q[l][k][j] = charge density at grid point j,k,l
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = first dimension of charge array, must be >= nx+ng
   nyv = second dimension of charge array, must be >= ny+ng
   nzv = third dimension of charge array, must be >= nz+ng
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
local data                                                            */
   int mxy1, mz1, noff, moff, loff;
   int i, j, k, l, nn, mm, ll;
   mxy1 = mx1*my1;
   mz1 = mxyz1/mxy1;
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,noff,moff,loff,nn,mm,ll)
   for (l = 0; l < mxyz1; l++) {
      i = l/mxy1;
      k = l - mxy1*i;
      loff = mz*i;
      ll = mz;
      if ((i+1)==mz1)
         ll = nzv - loff;
      j = k/mx1;
      moff = my*j;
      mm = my;
      if ((j+1)==my1)
         mm = nyv - moff;
      k = k - mx1*j;
      noff = mx*k;
      nn = mx;
      if ((k+1)==mx1)
         nn = nxv - noff;
/* zero charge in global array */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
#pragma ivdep
            for (i = 0; i < nn; i++) {
               q[i+noff+nxv*(j+moff+nyv*(k+loff))] = 0.0f;
            }
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cset_vzero3(float cu[], int mx, int my, int mz, int ndim, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1) {
/* for 3d code, this subroutine zeros out current density array.
   for Intel NUMA architecture with first touch policy, this associates
   array segments with appropriate threads
   OpenMP version
   input: all, output: cu
   cu[l][k][j][m] = charge density at grid point m,j,k,l
   mx/my/mz = number of grids in sorting cell in x/y/z
   ndim = first dimension of current array
   nxv = second dimension of current array, must be >= nx+ng
   nyv = third dimension of current array, must be >= ny+ng
   nzv = fourth dimension of current array, must be >= nz+ng
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mxyz1 = mx1*my1*mz1,
   where mz1 = (system length in z direction - 1)/mz + 1
local data                                                            */
   int mxy1, mz1, noff, moff, loff;
   int i, j, k, l, m, nn, mm, ll;
   mxy1 = mx1*my1;
   mz1 = mxyz1/mxy1;
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,nn,mm,ll)
   for (l = 0; l < mxyz1; l++) {
      i = l/mxy1;
      k = l - mxy1*i;
      loff = mz*i;
      ll = mz;
      if ((i+1)==mz1)
         ll = nzv - loff;
      j = k/mx1;
      moff = my*j;
      mm = my;
      if ((j+1)==my1)
         mm = nyv - moff;
      k = k - mx1*j;
      noff = mx*k;
      nn = mx;
      if ((k+1)==mx1)
         nn = nxv - noff;
/* zero current in global array */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
#pragma ivdep
            for (i = 0; i < nn; i++) {
               for (m = 0; m < ndim; m++) {
                  cu[m+ndim*(i+noff+nxv*(j+moff+nyv*(k+loff)))] = 0.0f;
               }
            }
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cset_cvzero3(float complex exyz[], int nx, int ny, int nz,
                  int ndim, int nxvh, int nyv, int nzv) {
/* for 3d code, this subroutine zeros out transverse field array.
   for Intel NUMA architecture with first touch policy, this associates
   array segments with appropriate threads
   OpenMP version
   input: all, output: exyz
   exyz[l][k][j][i] = complex transverse electric field
   nx/ny/nz = system length in x/y/z direction
   ndim = first dimension of field array
   nxvh = second dimension of field array, must be >= nxh
   nyv = third dimension of field array, must be >= ny
   nzv = fourth dimension of field array, must be >= nz
local data                                                            */
   int nxh, nyh, nzh, i, j, k, l, k1, l1, kj, lj, nxvyh;
   float complex zero;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxvyh = nxvh*nyv;
   zero = 0.0 + 0.0*_Complex_I;
/* loop over mode numbers */
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
#pragma omp parallel
   {
#pragma omp for nowait \
private(i,j,k,l,k1,l1,lj,kj)
      for (l = 1; l < nzh; l++) {
         lj = nxvyh*l;
         l1 = nxvyh*nz - lj;
         for (k = 1; k < nyh; k++) {
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            for (j = 1; j < nxh; j++) {
               for (i = 0; i < ndim; i++) {
                  exyz[i+ndim*(j+kj+lj)] = zero;
                  exyz[i+ndim*(j+k1+lj)] = zero;
                  exyz[i+ndim*(j+kj+l1)] = zero;
                  exyz[i+ndim*(j+k1+l1)] = zero;
               }
            }
         }
/* mode numbers kx = 0, nx/2 */
         for (k = 1; k < nyh; k++) {
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            for (i = 0; i < ndim; i++) {
               exyz[i+ndim*(kj+lj)] = zero;
               exyz[i+ndim*(k1+lj)] = zero;
               exyz[i+ndim*(kj+l1)] = zero;
               exyz[i+ndim*(k1+l1)] = zero;
            }
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
         for (j = 1; j < nxh; j++) {
            for (i = 0; i < ndim; i++) {
               exyz[i+ndim*(j+lj)] = zero;
               exyz[i+ndim*(j+k1+lj)] = zero;
               exyz[i+ndim*(j+l1)] = zero;
               exyz[i+ndim*(j+k1+l1)] = zero;
            }
         }
/* mode numbers kx = 0, nx/2 */
         for (i = 0; i < ndim; i++) {
            exyz[i+ndim*lj] = zero;
            exyz[i+ndim*(k1+lj)] = zero;
            exyz[i+ndim*l1] = zero;
            exyz[i+ndim*(k1+l1)] = zero;
         }
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
#pragma omp parallel for private(i,j,k,k1,kj)
   for (k = 1; k < nyh; k++) {
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      for (j = 1; j < nxh; j++) {
         for (i = 0; i < ndim; i++) {
            exyz[i+ndim*(j+kj)] = zero;
            exyz[i+ndim*(j+k1)] = zero;
            exyz[i+ndim*(j+kj+l1)] = zero;
            exyz[i+ndim*(j+k1+l1)] = zero;
         }
      }
/* mode numbers kx = 0, nx/2 */
      for (i = 0; i < ndim; i++) {
         exyz[i+ndim*kj] = zero;
         exyz[i+ndim*k1] = zero;
         exyz[i+ndim*(kj+l1)] = zero;
         exyz[i+ndim*(k1+l1)] = zero;
      }
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      for (i = 0; i < ndim; i++) {
         exyz[i+ndim*j] = zero;
         exyz[i+ndim*(j+k1)] = zero;
         exyz[i+ndim*(j+l1)] = zero;
         exyz[i+ndim*(j+k1+l1)] = zero;
      }
   }
   for (i = 0; i < ndim; i++) {
      exyz[i] = zero;
      exyz[i+ndim*k1] = zero;
      exyz[i+ndim*l1] = zero;
      exyz[i+ndim*(k1+l1)] = zero;
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void cdistr3_(float *part, float *vtx, float *vty, float *vtz,
              float *vdx, float *vdy, float *vdz, int *npx, int *npy,
              int *npz, int *idimp, int *nop, int *nx, int *ny, int *nz, 
              int *ipbc) {
   cdistr3(part,*vtx,*vty,*vtz,*vdx,*vdy,*vdz,*npx,*npy,*npz,*idimp,
           *nop,*nx,*ny,*nz,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdblkp3l_(float *part, int *kpic, int *nppmx, int *idimp, int *nop,
               int *mx, int *my, int *mz, int *mx1, int *my1,
               int *mxyz1, int *irc) {
   cdblkp3l(part,kpic,nppmx,*idimp,*nop,*mx,*my,*mz,*mx1,*my1,*mxyz1,
            irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppmovin3lt_(float *part, float *ppart, int *kpic, int *nppmx,
                  int *idimp, int *nop, int *mx, int *my, int *mz,
                  int *mx1, int *my1, int *mxyz1, int *irc) {
   cppmovin3lt(part,ppart,kpic,*nppmx,*idimp,*nop,*mx,*my,*mz,*mx1,*my1,
               *mxyz1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppmovin3ltp_(float *part, float *ppart, int *kpic, int *kp,
                   int *nppmx, int *idimp, int *nop, int *mx, int *my,
                   int *mz, int *mx1, int *my1, int *mxyz1, int *irc) {
   cppmovin3ltp(part,ppart,kpic,kp,*nppmx,*idimp,*nop,*mx,*my,*mz,*mx1,
                *my1,*mxyz1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcheck3lt_(float *ppart, int *kpic, int *idimp, int *nppmx,
                  int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                  int *mx1, int *my1, int *mz1, int *irc) {
   cppcheck3lt(ppart,kpic,*idimp,*nppmx,*nx,*ny,*nz,*mx,*my,*mz,*mx1,
               *my1,*mz1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppush3lt_(float *ppart, float *fxyz, int *kpic, float *qbm,
                 float *dt, float *ek, int *idimp, int *nppmx, int *nx,
                 int *ny, int *nz, int *mx, int *my, int *mz, int *nxv,
                 int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1,
                 int *ipbc) {
   cgppush3lt(ppart,fxyz,kpic,*qbm,*dt,ek,*idimp,*nppmx,*nx,*ny,*nz,*mx,
              *my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppushf3lt_(float *ppart, float *fxyz, int *kpic, int *ncl,
                  int *ihole, float *qbm, float *dt, float *ek,
                  int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                  int *mx, int *my, int *mz, int *nxv, int *nyv,
                  int *nzv, int *mx1, int *my1, int *mxyz1, int *ntmax,
                  int *irc) {
   cgppushf3lt(ppart,fxyz,kpic,ncl,ihole,*qbm,*dt,ek,*idimp,*nppmx,*nx,
               *ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1,
               *ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppush3lt_(float *ppart, float *fxyz, int *kpic, float *qbm,
                  float *dt, float *ek, int *idimp, int *nppmx, int *nx,
                  int *ny, int *nz, int *mx, int *my, int *mz, int *nxv,
                  int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1,
                  int *ipbc) {
   cvgppush3lt(ppart,fxyz,kpic,*qbm,*dt,ek,*idimp,*nppmx,*nx,*ny,*nz,
               *mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppushf3lt_(float *ppart, float *fxyz, int *kpic, int *ncl,
                   int *ihole, float *qbm, float *dt, float *ek,
                   int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                   int *mx, int *my, int *mz, int *nxv, int *nyv,
                   int *nzv, int *mx1, int *my1, int *mxyz1, int *ntmax,
                   int *irc) {
   cvgppushf3lt(ppart,fxyz,kpic,ncl,ihole,*qbm,*dt,ek,*idimp,*nppmx,*nx,
                *ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1,
                *ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cv2gppush3lt_(float *ppart, float *fxyz, int *kpic, float *qbm,
                   float *dt, float *ek, int *idimp, int *nppmx,
                   int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                   int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                   int *mxyz1, int *ipbc) {
   cv2gppush3lt(ppart,fxyz,kpic,*qbm,*dt,ek,*idimp,*nppmx,*nx,*ny,*nz,
                *mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cv2gppushf3lt_(float *ppart, float *fxyz, int *kpic, int *ncl,
                    int *ihole, float *qbm, float *dt, float *ek,
                    int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                    int *mx, int *my, int *mz, int *nxv, int *nyv,
                    int *nzv, int *mx1, int *my1, int *mxyz1,
                    int *ntmax, int *irc) {
   cv2gppushf3lt(ppart,fxyz,kpic,ncl,ihole,*qbm,*dt,ek,*idimp,*nppmx,
                 *nx,*ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,
                 *mxyz1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppost3lt_(float *ppart, float *q, int *kpic, float *qm,
                 int *nppmx, int *idimp, int *mx, int *my, int *mz,
                 int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                 int *mxyz1) {
   cgppost3lt(ppart,q,kpic,*qm,*nppmx,*idimp,*mx,*my,*mz,*nxv,*nyv,*nzv,
              *mx1,*my1,*mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppost3lt_(float *ppart, float *q, int *kpic, float *qm,
                  int *nppmx, int *idimp, int *mx, int *my, int *mz,
                  int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                  int *mxyz1) {
   cvgppost3lt(ppart,q,kpic,*qm,*nppmx,*idimp,*mx,*my,*mz,*nxv,*nyv,
               *nzv,*mx1,*my1,*mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cviscan2_(int *isdata, int *mb, int *nths) {
   cviscan2(isdata,mb,*nths);
   return;
}

/*--------------------------------------------------------------------*/
void cpporder3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                  int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                  int *nz, int *mx, int *my, int *mz, int *mx1,
                  int *my1, int *mz1, int *npbmx, int *ntmax,
                  int *irc) {
   cpporder3lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*nx,*ny,*nz,
               *mx,*my,*mz,*mx1,*my1,*mz1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                   int *ihole, int *idimp, int *nppmx, int *mx1,
                   int *my1, int *mz1, int *npbmx, int *ntmax,
                   int *irc) {
   cpporderf3lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*mx1,*my1,
                *mz1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvpporder3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                   int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                   int *nz, int *mx, int *my, int *mz, int *mx1,
                   int *my1, int *mz1, int *npbmx, int *ntmax,
                   int *irc) {
   cvpporder3lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*nx,*ny,*nz,
                *mx,*my,*mz,*mx1,*my1,*mz1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvpporderf3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                    int *ihole, int *idimp, int *nppmx, int *mx1,
                    int *my1, int *mz1, int *npbmx, int *ntmax,
                    int *irc) {
   cvpporderf3lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*mx1,*my1,
                 *mz1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cv2pporderf3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                     int *ihole, int *idimp, int *nppmx, int *mx1, 
                     int *my1, int *mz1, int *npbmx, int *ntmax,
                     int *irc) {
   cv2pporderf3lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*mx1,*my1,
                  *mz1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void ccguard3l_(float *fxyz, int *nx, int *ny, int *nz, int *nxe,
                int *nye, int *nze) {
   ccguard3l(fxyz,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void caguard3l_(float *q, int *nx, int *ny, int *nz, int *nxe, int *nye,
                int *nze) {
   caguard3l(q,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void cvmpois33_(float complex *q, float complex *fxyz, int *isign,
                float complex *ffc, float *ax, float *ay, float *az,
                float *affp, float *we, int *nx, int *ny, int *nz,
                int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
                int *nzhd) {
 cvmpois33(q,fxyz,*isign,ffc,*ax,*ay,*az,*affp,we,*nx,*ny,*nz,*nxvh,
           *nyv,*nzv,*nxhd,*nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                  int *indz, int *nxhyzd, int *nxyzhd) {
   cwfft3rinit(mixup,sct,*indx,*indy,*indz,*nxhyzd,*nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rvmx_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *indz,
                 int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                 int *nxyzhd) {
   cwfft3rvmx(f,*isign,mixup,sct,*indx,*indy,*indz,*nxhd,*nyd,*nzd,
              *nxhyzd,*nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rvm3_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *indz,
                 int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                 int *nxyzhd) {
   cwfft3rvm3(f,*isign,mixup,sct,*indx,*indy,*indz,*nxhd,*nyd,*nzd,
              *nxhyzd,*nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cset_szero3_(float *q, int *mx, int *my, int *mz, int *nxv,
                  int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1) {
   cset_szero3(q,*mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cset_vzero3_(float *cu, int *mx, int *my, int *mz, int *ndim,
                  int *nxv, int *nyv, int *nzv, int *mx1, int *my1, 
                  int *mxyz1) {
   cset_vzero3(cu,*mx,*my,*mz,*ndim,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cset_cvzero3_(float complex *exyz, int *nx, int *ny, int *nz,
                   int *ndim, int *nxvh, int *nyv, int *nzv) {
   cset_cvzero3(exyz,*nx,*ny,*nz,*ndim,*nxvh,*nyv,*nzv);
   return;
}
