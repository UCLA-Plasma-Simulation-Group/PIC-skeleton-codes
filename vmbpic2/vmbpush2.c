/* C Library for Skeleton 2-1/2D Electromagnetic OpenMP/Vector PIC */
/* Code                                                            */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "vmbpush2.h"

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
void cdistr2h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int npy, int idimp, int nop,
              int nx, int ny, int ipbc) {
/* for 2-1/2d code, this subroutine calculates initial particle
   co-ordinates and velocities with uniform density and maxwellian
   velocity with drift
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = velocity vx of particle n
   part[n][3] = velocity vy of particle n
   part[n][4] = velocity vz of particle n
   vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
   vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
   npx/npy = initial number of particles distributed in x/y direction
   idimp = size of phase space = 5
   nop = number of particles
   nx/ny = system length in x/y direction
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   ranorm = gaussian random number with zero mean and unit variance
local data                                                            */
   int j, k, k1, npxy;
   float edgelx, edgely, at1, at2, at3, sum1, sum2, sum3;
   double dsum1, dsum2, dsum3;
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
      edgely = 0.0;
      at1 = (float) (nx-2)/(float) npx;
      at2 = (float) ny/(float) npy;
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
      part[4+idimp*j] = vtz*ranorm();
   }
/* add correct drift */
   dsum1 = 0.0;
   dsum2 = 0.0;
   dsum3 = 0.0;
   for (j = 0; j < npxy; j++) {
      dsum1 += part[2+idimp*j];
      dsum2 += part[3+idimp*j];
      dsum3 += part[4+idimp*j];
   }
   sum1 = dsum1;
   sum2 = dsum2;
   sum3 = dsum3;
   at1 = 1.0/(float) npxy;
   sum1 = at1*sum1 - vdx;
   sum2 = at1*sum2 - vdy;
   sum3 = at1*sum3 - vdz;
   for (j = 0; j < npxy; j++) {
      part[2+idimp*j] -= sum1;
      part[3+idimp*j] -= sum2;
      part[4+idimp*j] -= sum3;
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
void cgbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                  float qbm, float dt, float dtc, float *ek, int idimp,
                  int nppmx, int nx, int ny, int mx, int my, int nxv,
                  int nyv, int mx1, int mxy1, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field. Using the Boris Mover.
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   119 flops/particle, 1 divide, 29 loads, 5 stores
   input: all, output: ppart, ek
   velocity equations used are:
   vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fz(x(t),y(t))*dt)
   where q/m is charge/mass, and the rotation matrix is given by:
      rot[0] = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[1] = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[2] = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[3] = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[4] = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[5] = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[6] = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[7] = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[8] = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
   and om**2 = omx**2 + omy**2 + omz**2
   the rotation matrix is determined by:
   omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
   omz = (q/m)*bz(x(t),y(t)).
   position equations used are:
   x(t+dt)=x(t) + vx(t+dt/2)*dt
   y(t+dt)=y(t) + vy(t+dt/2)*dt
   fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
        (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
   idimp = size of phase space = 5
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
#define N 4
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nn, mm, nm;
   float qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt;
   float anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, vx, vy, vz;
   float sfxy[N*MXV*MYV], sbxy[N*MXV*MYV];
/* float sfxy[N*(mx+1)*(my+1)], sbxy[N*(mx+1)*(my+1)]; */
   double sum1, sum2;
   mxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
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
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,nm,x,y,vx,vy,vz,dxp,dyp,amx, \
amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2, \
rot3,rot4,rot5,rot6,rot7,rot8,rot9,sum1,sfxy,sbxy) \
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
            sfxy[N*(i+mxv*j)] = fxy[N*(i+noff+nxv*(j+moff))];
            sfxy[1+N*(i+mxv*j)] = fxy[1+N*(i+noff+nxv*(j+moff))];
            sfxy[2+N*(i+mxv*j)] = fxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[N*(i+mxv*j)] = bxy[N*(i+noff+nxv*(j+moff))];
            sbxy[1+N*(i+mxv*j)] = bxy[1+N*(i+noff+nxv*(j+moff))];
            sbxy[2+N*(i+mxv*j)] = bxy[2+N*(i+noff+nxv*(j+moff))];
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
         nm = N*(nn - noff + mxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + N;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += N*mxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + N;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + N;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += N*mxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + N;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
/* time-centered kinetic energy */
         sum1 += (acx*acx + acy*acy + acz*acz);
/* calculate cyclotron frequency */
         omxt = qtmh*ox;
         omyt = qtmh*oy;
         omzt = qtmh*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0f/(1.0f + omt);
         omt = 0.5f*(1.0f - omt);
         rot4 = omxt*omyt;
         rot7 = omxt*omzt;
         rot8 = omyt*omzt;
         rot1 = omt + omxt*omxt;
         rot5 = omt + omyt*omyt;
         rot9 = omt + omzt*omzt;
         rot2 = omzt + rot4;
         rot4 -= omzt;
         rot3 = -omyt + rot7;
         rot7 += omyt;
         rot6 = omxt + rot8;
         rot8 -= omxt;
/* new velocity */
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* new position */
         dx = x + vx*dtc;
         dy = y + vy*dtc;
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
         ppart[j+4*nppmx+npoff] = vz;
      }
      sum2 += sum1;
   }
/* normalize kinetic energy */
   *ek += 0.5*sum2;
   return;
#undef N
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgbppushf23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                   int ncl[], int ihole[], float qbm, float dt,
                   float dtc, float *ek, int idimp, int nppmx, int nx,
                   int ny, int mx, int my, int nxv, int nyv, int mx1,
                   int mxy1, int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field. Using the Boris Mover.
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   119 flops/particle, 1 divide, 29 loads, 5 stores
   input: all except ncl, ihole, irc, output: ppart, ncl, ihole, irc, ek
   velocity equations used are:
   vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fz(x(t),y(t))*dt)
   where q/m is charge/mass, and the rotation matrix is given by:
      rot[0] = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[1] = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[2] = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[3] = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[4] = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[5] = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[6] = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[7] = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[8] = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
   and om**2 = omx**2 + omy**2 + omz**2
   the rotation matrix is determined by:
   omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
   omz = (q/m)*bz(x(t),y(t)).
   position equations used are:
   x(t+dt)=x(t) + vx(t+dt/2)*dt
   y(t+dt)=y(t) + vy(t+dt/2)*dt
   fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
        (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
   idimp = size of phase space = 5
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
#define N 4
   int noff, moff, npoff, npp, mxv;
   int i, j, k, ih, nh, nn, mm, nm;
   float qtmh, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float x, y, vx, vy, vz;
   float sfxy[N*MXV*MYV], sbxy[N*MXV*MYV];
/* float sfxy[N*(mx+1)*(my+1)], sbxy[N*(mx+1)*(my+1)]; */
   double sum1, sum2;
   mxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   sum2 = 0.0;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,nm,ih,nh,x,y,vx,vy,vz,dxp,dyp, \
amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1, \
rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,edgery, \
sum1,sfxy,sbxy) \
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
            sfxy[N*(i+mxv*j)] = fxy[N*(i+noff+nxv*(j+moff))];
            sfxy[1+N*(i+mxv*j)] = fxy[1+N*(i+noff+nxv*(j+moff))];
            sfxy[2+N*(i+mxv*j)] = fxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[N*(i+mxv*j)] = bxy[N*(i+noff+nxv*(j+moff))];
            sbxy[1+N*(i+mxv*j)] = bxy[1+N*(i+noff+nxv*(j+moff))];
            sbxy[2+N*(i+mxv*j)] = bxy[2+N*(i+noff+nxv*(j+moff))];
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
         nm = N*(nn - noff + mxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + N;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += N*mxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + N;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + N;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += N*mxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + N;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
/* time-centered kinetic energy */
         sum1 += (acx*acx + acy*acy + acz*acz);
/* calculate cyclotron frequency */
         omxt = qtmh*ox;
         omyt = qtmh*oy;
         omzt = qtmh*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0f/(1.0f + omt);
         omt = 0.5f*(1.0f - omt);
         rot4 = omxt*omyt;
         rot7 = omxt*omzt;
         rot8 = omyt*omzt;
         rot1 = omt + omxt*omxt;
         rot5 = omt + omyt*omyt;
         rot9 = omt + omzt*omzt;
         rot2 = omzt + rot4;
         rot4 -= omzt;
         rot3 = -omyt + rot7;
         rot7 += omyt;
         rot6 = omxt + rot8;
         rot8 -= omxt;
/* new velocity */
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* new position */
         dx = x + vx*dtc;
         dy = y + vy*dtc;
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
         ppart[j+4*nppmx+npoff] = vz;
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
   *ek += 0.5*sum2;
   return;
#undef N
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgrbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                   float qbm, float dt, float dtc, float ci, float *ek,
                   int idimp, int nppmx, int nx, int ny, int mx, int my,
                   int nxv, int nyv, int mx1, int mxy1, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
   input: all, output: ppart, ek
   momentum equations used are:
   px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fz(x(t),y(t))*dt)
   where q/m is charge/mass, and the rotation matrix is given by:
      rot[0] = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[1] = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[2] = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[3] = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[4] = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[5] = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[6] = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[7] = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[8] = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
   and om**2 = omx**2 + omy**2 + omz**2
   the rotation matrix is determined by:
   omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
   omz = (q/m)*bz(x(t),y(t))*gami,
   where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
   position equations used are:
   x(t+dt) = x(t) + px(t+dt/2)*dtg
   y(t+dt) = y(t) + py(t+dt/2)*dtg
   where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
   pz(t+dt/2)*pz(t+dt/2))*ci*ci)
   fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   idimp = size of phase space = 5
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
#define N 4
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nn, mm, nm;
   float qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg;
   float omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, vx, vy, vz;
   float sfxy[N*MXV*MYV], sbxy[N*MXV*MYV];
/* float sfxy[N*(mx+1)*(my+1)], sbxy[N*(mx+1)*(my+1)]; */
   double sum1, sum2;
   mxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
   ci2 = ci*ci;
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
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,nm,x,y,vx,vy,vz,dxp,dyp,amx, \
amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2, \
rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,dtg,sum1,sfxy,sbxy) \
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
            sfxy[N*(i+mxv*j)] = fxy[N*(i+noff+nxv*(j+moff))];
            sfxy[1+N*(i+mxv*j)] = fxy[1+N*(i+noff+nxv*(j+moff))];
            sfxy[2+N*(i+mxv*j)] = fxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[N*(i+mxv*j)] = bxy[N*(i+noff+nxv*(j+moff))];
            sbxy[1+N*(i+mxv*j)] = bxy[1+N*(i+noff+nxv*(j+moff))];
            sbxy[2+N*(i+mxv*j)] = bxy[2+N*(i+noff+nxv*(j+moff))];
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
         nm = N*(nn - noff + mxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + N;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += N*mxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + N;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + N;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += N*mxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + N;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
/* find inverse gamma */
         p2 = acx*acx + acy*acy + acz*acz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* renormalize magnetic field */
         qtmg = qtmh*gami;
/* time-centered kinetic energy */
         sum1 += gami*p2/(1.0f + gami);
/* calculate cyclotron frequency */
         omxt = qtmg*ox;
         omyt = qtmg*oy;
         omzt = qtmg*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0f/(1.0f + omt);
         omt = 0.5f*(1.0f - omt);
         rot4 = omxt*omyt;
         rot7 = omxt*omzt;
         rot8 = omyt*omzt;
         rot1 = omt + omxt*omxt;
         rot5 = omt + omyt*omyt;
         rot9 = omt + omzt*omzt;
         rot2 = omzt + rot4;
         rot4 -= omzt;
         rot3 = -omyt + rot7;
         rot7 += omyt;
         rot6 = omxt + rot8;
         rot8 -= omxt;
/* new velocity */
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* update inverse gamma */
         p2 = vx*vx + vy*vy + vz*vz;
         dtg = dtc/sqrtf(1.0f + p2*ci2);
/* new position */
         dx = x + vx*dtg;
         dy = y + vy*dtg;
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
/* set new momentum */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
         ppart[j+4*nppmx+npoff] = vz;
      }
      sum2 += sum1;
   }
/* normalize kinetic energy */
   *ek += sum2;
   return;
#undef N
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgrbppushf23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                    int ncl[], int ihole[], float qbm, float dt,
                    float dtc, float ci, float *ek, int idimp,
                    int nppmx, int nx, int ny, int mx, int my, int nxv,
                    int nyv, int mx1, int mxy1, int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   with periodic boundary conditions.
   Using the Boris Mover.
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
   input: all except ncl, ihole, irc, output: ppart, ncl, ihole, irc, ek
   momentum equations used are:
   px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fz(x(t),y(t))*dt)
   where q/m is charge/mass, and the rotation matrix is given by:
      rot[0] = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[1] = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[2] = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[3] = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[4] = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[5] = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[6] = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[7] = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[8] = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
   and om**2 = omx**2 + omy**2 + omz**2
   the rotation matrix is determined by:
   omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
   omz = (q/m)*bz(x(t),y(t))*gami,
   where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
   position equations used are:
   x(t+dt) = x(t) + px(t+dt/2)*dtg
   y(t+dt) = y(t) + py(t+dt/2)*dtg
   where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
   pz(t+dt/2)*pz(t+dt/2))*ci*ci)
   fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   idimp = size of phase space = 5
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
#define N 4
   int noff, moff, npoff, npp, mxv;
   int i, j, k, ih, nh, nn, mm, nm;
   float qtmh, ci2, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz;
   float acx, acy, acz, p2, gami, qtmg, dtg, omxt, omyt, omzt, omt;
   float anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float x, y, vx, vy, vz;
   float sfxy[N*MXV*MYV], sbxy[N*MXV*MYV];
/* float sfxy[N*(mx+1)*(my+1)], sbxy[N*(mx+1)*(my+1)]; */
   double sum1, sum2;
   mxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
   ci2 = ci*ci;
   anx = (float) nx;
   any = (float) ny;
   sum2 = 0.0;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,nm,ih,nh,x,y,vx,vy,vz,dxp,dyp, \
amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1, \
rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,edgery,p2, \
gami,qtmg,dtg,sum1,sfxy,sbxy) \
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
            sfxy[N*(i+mxv*j)] = fxy[N*(i+noff+nxv*(j+moff))];
            sfxy[1+N*(i+mxv*j)] = fxy[1+N*(i+noff+nxv*(j+moff))];
            sfxy[2+N*(i+mxv*j)] = fxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[N*(i+mxv*j)] = bxy[N*(i+noff+nxv*(j+moff))];
            sbxy[1+N*(i+mxv*j)] = bxy[1+N*(i+noff+nxv*(j+moff))];
            sbxy[2+N*(i+mxv*j)] = bxy[2+N*(i+noff+nxv*(j+moff))];
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
         nm = N*(nn - noff + mxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + N;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += N*mxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + N;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + N;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += N*mxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + N;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
/* find inverse gamma */
         p2 = acx*acx + acy*acy + acz*acz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* renormalize magnetic field */
         qtmg = qtmh*gami;
/* time-centered kinetic energy */
         sum1 += gami*p2/(1.0f + gami);
/* calculate cyclotron frequency */
         omxt = qtmg*ox;
         omyt = qtmg*oy;
         omzt = qtmg*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0f/(1.0f + omt);
         omt = 0.5f*(1.0f - omt);
         rot4 = omxt*omyt;
         rot7 = omxt*omzt;
         rot8 = omyt*omzt;
         rot1 = omt + omxt*omxt;
         rot5 = omt + omyt*omyt;
         rot9 = omt + omzt*omzt;
         rot2 = omzt + rot4;
         rot4 -= omzt;
         rot3 = -omyt + rot7;
         rot7 += omyt;
         rot6 = omxt + rot8;
         rot8 -= omxt;
/* new momentum */
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* update inverse gamma */
         p2 = vx*vx + vy*vy + vz*vz;
         dtg = dtc/sqrtf(1.0f + p2*ci2);
/* new position */
         dx = x + vx*dtg;
         dy = y + vy*dtg;
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
/* set new momentum */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
         ppart[j+4*nppmx+npoff] = vz;
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
   *ek += sum2;
   return;
#undef N
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cvgbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                   float qbm, float dt, float dtc, float *ek, int idimp,
                   int nppmx, int nx, int ny, int mx, int my, int nxv,
                   int nyv, int mx1, int mxy1, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field. Using the Boris Mover.
   vectorizable/OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   119 flops/particle, 1 divide, 29 loads, 5 stores
   input: all, output: ppart, ek
   velocity equations used are:
   vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fz(x(t),y(t))*dt)
   where q/m is charge/mass, and the rotation matrix is given by:
      rot[0] = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[1] = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[2] = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[3] = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[4] = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[5] = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[6] = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[7] = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[8] = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
   and om**2 = omx**2 + omy**2 + omz**2
   the rotation matrix is determined by:
   omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
   omz = (q/m)*bz(x(t),y(t)).
   position equations used are:
   x(t+dt)=x(t) + vx(t+dt/2)*dt
   y(t+dt)=y(t) + vy(t+dt/2)*dt
   fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
        (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
   idimp = size of phase space = 5
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
#define N 4
   int noff, moff, npoff, npp, ipp, joff, nps;
   int i, j, k, m, nn, mm, nm, lxv;
   float qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt;
   float anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, vx, vy, vz;
   float sfxy[N*MXV*MYV], sbxy[N*MXV*MYV];
/* float sfxy[N*(mx+1)*(my+1)], sbxy[N*(mx+1)*(my+1)]; */
/* scratch arrays */
   int n[NPBLK];
   float s1[NPBLK*LVECT], s2[NPBLK*LVECT], t[NPBLK*2];
   double sum1, sum2;
   lxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
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
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,noff,moff,npp,npoff,ipp,joff,nps,nn,mm,nm,x,y,vx,vy,vz, \
dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm, \
rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sum1,sfxy,sbxy,n,s1,s2,t) \
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
            sfxy[N*(i+lxv*j)] = fxy[N*(i+noff+nxv*(j+moff))];
            sfxy[1+N*(i+lxv*j)] = fxy[1+N*(i+noff+nxv*(j+moff))];
            sfxy[2+N*(i+lxv*j)] = fxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[N*(i+lxv*j)] = bxy[N*(i+noff+nxv*(j+moff))];
            sbxy[1+N*(i+lxv*j)] = bxy[1+N*(i+noff+nxv*(j+moff))];
            sbxy[2+N*(i+lxv*j)] = bxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
      sum1 = 0.0;
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
            n[j] = N*(nn - noff + lxv*(mm - moff));
            amx = 1.0f - dxp;
            amy = 1.0f - dyp;
            s1[j] = amx*amy;
            s1[j+NPBLK] = dxp*amy;
            s1[j+2*NPBLK] = amx*dyp;
            s1[j+3*NPBLK] = dxp*dyp;
            t[j] = x;
            t[j+NPBLK] = y;
         }
/* find acceleration */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + N*(lxv - 2);
            dx = 0.0f;
            dy = 0.0f;
            dz = 0.0f;
            ox = 0.0f;
            oy = 0.0f;
            oz = 0.0f;
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               if (i > 1)
                  nn = mm;
               dx += sfxy[N*i+nn]*s1[j+NPBLK*i];
               dy += sfxy[1+N*i+nn]*s1[j+NPBLK*i];
               dz += sfxy[2+N*i+nn]*s1[j+NPBLK*i];
               ox += sbxy[N*i+nn]*s1[j+NPBLK*i];
               oy += sbxy[1+N*i+nn]*s1[j+NPBLK*i];
               oz += sbxy[2+N*i+nn]*s1[j+NPBLK*i];
            }
            s1[j] = dx;
            s1[j+NPBLK] = dy;
            s1[j+2*NPBLK] = dz;
            s2[j] = ox;
            s2[j+NPBLK] = oy;
            s2[j+2*NPBLK] = oz;
         }
/* new velocity */
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
/* calculate half impulse */
            dx = qtmh*s1[j];
            dy = qtmh*s1[j+NPBLK];
            dz = qtmh*s1[j+2*NPBLK];
/* half acceleration */
            acx = ppart[j+joff+2*nppmx+npoff] + dx;
            acy = ppart[j+joff+3*nppmx+npoff] + dy;
            acz = ppart[j+joff+4*nppmx+npoff] + dz;
/* time-centered kinetic energy */
            sum1 += (acx*acx + acy*acy + acz*acz);
/* calculate cyclotron frequency */
            omxt = qtmh*s2[j];
            omyt = qtmh*s2[j+NPBLK];
            omzt = qtmh*s2[j+2*NPBLK];
/* calculate rotation matrix */
            omt = omxt*omxt + omyt*omyt + omzt*omzt;
            anorm = 2.0f/(1.0f + omt);
            omt = 0.5f*(1.0f - omt);
            rot4 = omxt*omyt;
            rot7 = omxt*omzt;
            rot8 = omyt*omzt;
            rot1 = omt + omxt*omxt;
            rot5 = omt + omyt*omyt;
            rot9 = omt + omzt*omzt;
            rot2 = omzt + rot4;
            rot4 -= omzt;
            rot3 = -omyt + rot7;
            rot7 += omyt;
            rot6 = omxt + rot8;
            rot8 -= omxt;
/* new velocity */
            vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
            vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
            vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* new position */
            s1[j] = x + vx*dtc;
            s1[j+NPBLK] = y + vy*dtc;
            s2[j] = vx;
            s2[j+NPBLK] = vy;
            s2[j+2*NPBLK] = vz;
         }
/* check boundary conditions */
#pragma novector
         for (j = 0; j < NPBLK; j++) {
            dx = s1[j];
            dy = s1[j+NPBLK];
            vx = s2[j];
            vy = s2[j+NPBLK];
            vz = s2[j+2*NPBLK];
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
            ppart[j+joff+4*nppmx+npoff] = vz;
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
         nm = N*(nn - noff + lxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + N;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += N*lxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + N;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + N;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += N*lxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + N;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
/* time-centered kinetic energy */
         sum1 += (acx*acx + acy*acy + acz*acz);
/* calculate cyclotron frequency */
         omxt = qtmh*ox;
         omyt = qtmh*oy;
         omzt = qtmh*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0f/(1.0f + omt);
         omt = 0.5f*(1.0f - omt);
         rot4 = omxt*omyt;
         rot7 = omxt*omzt;
         rot8 = omyt*omzt;
         rot1 = omt + omxt*omxt;
         rot5 = omt + omyt*omyt;
         rot9 = omt + omzt*omzt;
         rot2 = omzt + rot4;
         rot4 -= omzt;
         rot3 = -omyt + rot7;
         rot7 += omyt;
         rot6 = omxt + rot8;
         rot8 -= omxt;
/* new velocity */
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* new position */
         dx = x + vx*dtc;
         dy = y + vy*dtc;
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
         ppart[j+4*nppmx+npoff] = vz;
      }
      sum2 += sum1;
   }
/* normalize kinetic energy */
   *ek += 0.5*sum2;
   return;
#undef N
#undef LVECT
#undef NPBLK
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cvgbppushf23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                    int ncl[], int ihole[], float qbm, float dt,
                    float dtc, float *ek, int idimp, int nppmx, int nx,
                    int ny, int mx, int my, int nxv, int nyv, int mx1,
                    int mxy1, int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field. Using the Boris Mover.
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   vectorizable/OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   119 flops/particle, 1 divide, 29 loads, 5 stores
   input: all except ncl, ihole, irc, output: ppart, ncl, ihole, irc, ek
   velocity equations used are:
   vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fz(x(t),y(t))*dt)
   where q/m is charge/mass, and the rotation matrix is given by:
      rot[0] = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[1] = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[2] = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[3] = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[4] = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[5] = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[6] = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[7] = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[8] = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
   and om**2 = omx**2 + omy**2 + omz**2
   the rotation matrix is determined by:
   omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
   omz = (q/m)*bz(x(t),y(t)).
   position equations used are:
   x(t+dt)=x(t) + vx(t+dt/2)*dt
   y(t+dt)=y(t) + vy(t+dt/2)*dt
   fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
        (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
   idimp = size of phase space = 5
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
#define N 4
   int noff, moff, npoff, npp, ipp, joff, nps;
   int i, j, k, m, ih, nh, nn, mm, nm, lxv;
   float qtmh, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float x, y, vx, vy, vz;
   float sfxy[N*MXV*MYV], sbxy[N*MXV*MYV];
/* float sfxy[N*(mx+1)*(my+1)], sbxy[N*(mx+1)*(my+1)]; */
/* scratch arrays */
   int n[NPBLK];
   float s1[NPBLK*LVECT], s2[NPBLK*LVECT], t[NPBLK*2];
   double sum1, sum2;
   lxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   sum2 = 0.0;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,noff,moff,npp,npoff,ipp,joff,nps,nn,mm,nm,ih,nh,x,y, \
vx,vy,vz,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt, \
omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely, \
edgerx,edgery,sum1,sfxy,sbxy,n,s1,s2,t) \
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
            sfxy[N*(i+lxv*j)] = fxy[N*(i+noff+nxv*(j+moff))];
            sfxy[1+N*(i+lxv*j)] = fxy[1+N*(i+noff+nxv*(j+moff))];
            sfxy[2+N*(i+lxv*j)] = fxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[N*(i+lxv*j)] = bxy[N*(i+noff+nxv*(j+moff))];
            sbxy[1+N*(i+lxv*j)] = bxy[1+N*(i+noff+nxv*(j+moff))];
            sbxy[2+N*(i+lxv*j)] = bxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
      sum1 = 0.0;
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
            n[j] = N*(nn - noff + lxv*(mm - moff));
            amx = 1.0f - dxp;
            amy = 1.0f - dyp;
            s1[j] = amx*amy;
            s1[j+NPBLK] = dxp*amy;
            s1[j+2*NPBLK] = amx*dyp;
            s1[j+3*NPBLK] = dxp*dyp;
            t[j] = x;
            t[j+NPBLK] = y;
         }
/* find acceleration */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + N*(lxv - 2);
            dx = 0.0f;
            dy = 0.0f;
            dz = 0.0f;
            ox = 0.0f;
            oy = 0.0f;
            oz = 0.0f;
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               if (i > 1)
                  nn = mm;
               dx += sfxy[N*i+nn]*s1[j+NPBLK*i];
               dy += sfxy[1+N*i+nn]*s1[j+NPBLK*i];
               dz += sfxy[2+N*i+nn]*s1[j+NPBLK*i];
               ox += sbxy[N*i+nn]*s1[j+NPBLK*i];
               oy += sbxy[1+N*i+nn]*s1[j+NPBLK*i];
               oz += sbxy[2+N*i+nn]*s1[j+NPBLK*i];
            }
            s1[j] = dx;
            s1[j+NPBLK] = dy;
            s1[j+2*NPBLK] = dz;
            s2[j] = ox;
            s2[j+NPBLK] = oy;
            s2[j+2*NPBLK] = oz;
         }
/* new velocity */
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
/* calculate half impulse */
            dx = qtmh*s1[j];
            dy = qtmh*s1[j+NPBLK];
            dz = qtmh*s1[j+2*NPBLK];
/* half acceleration */
            acx = ppart[j+joff+2*nppmx+npoff] + dx;
            acy = ppart[j+joff+3*nppmx+npoff] + dy;
            acz = ppart[j+joff+4*nppmx+npoff] + dz;
/* time-centered kinetic energy */
            sum1 += (acx*acx + acy*acy + acz*acz);
/* calculate cyclotron frequency */
            omxt = qtmh*s2[j];
            omyt = qtmh*s2[j+NPBLK];
            omzt = qtmh*s2[j+2*NPBLK];
/* calculate rotation matrix */
            omt = omxt*omxt + omyt*omyt + omzt*omzt;
            anorm = 2.0f/(1.0f + omt);
            omt = 0.5f*(1.0f - omt);
            rot4 = omxt*omyt;
            rot7 = omxt*omzt;
            rot8 = omyt*omzt;
            rot1 = omt + omxt*omxt;
            rot5 = omt + omyt*omyt;
            rot9 = omt + omzt*omzt;
            rot2 = omzt + rot4;
            rot4 -= omzt;
            rot3 = -omyt + rot7;
            rot7 += omyt;
            rot6 = omxt + rot8;
            rot8 -= omxt;
/* new velocity */
            vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
            vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
            vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* new position */
            s1[j] = x + vx*dtc;
            s1[j+NPBLK] = y + vy*dtc;
            s2[j] = vx;
            s2[j+NPBLK] = vy;
            s2[j+2*NPBLK] = vz;
         }
/* check boundary conditions */
#pragma novector
         for (j = 0; j < NPBLK; j++) {
            dx = s1[j];
            dy = s1[j+NPBLK];
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
            ppart[j+joff+2*nppmx+npoff] = s2[j];
            ppart[j+joff+3*nppmx+npoff] = s2[j+NPBLK];
            ppart[j+joff+4*nppmx+npoff] = s2[j+2*NPBLK];
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
         nm = N*(nn - noff + lxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + N;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += N*lxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + N;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + N;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += N*lxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + N;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
/* time-centered kinetic energy */
         sum1 += (acx*acx + acy*acy + acz*acz);
/* calculate cyclotron frequency */
         omxt = qtmh*ox;
         omyt = qtmh*oy;
         omzt = qtmh*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0f/(1.0f + omt);
         omt = 0.5f*(1.0f - omt);
         rot4 = omxt*omyt;
         rot7 = omxt*omzt;
         rot8 = omyt*omzt;
         rot1 = omt + omxt*omxt;
         rot5 = omt + omyt*omyt;
         rot9 = omt + omzt*omzt;
         rot2 = omzt + rot4;
         rot4 -= omzt;
         rot3 = -omyt + rot7;
         rot7 += omyt;
         rot6 = omxt + rot8;
         rot8 -= omxt;
/* new velocity */
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* new position */
         dx = x + vx*dtc;
         dy = y + vy*dtc;
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
         ppart[j+4*nppmx+npoff] = vz;
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
   *ek += 0.5*sum2;
   return;
#undef N
#undef LVECT
#undef NPBLK
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cvgrbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                    float qbm, float dt, float dtc, float ci, float *ek,
                    int idimp, int nppmx, int nx, int ny, int mx,
                    int my, int nxv, int nyv, int mx1, int mxy1,
                    int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   vectorizable/OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
   input: all, output: ppart, ek
   momentum equations used are:
   px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fz(x(t),y(t))*dt)
   where q/m is charge/mass, and the rotation matrix is given by:
      rot[0] = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[1] = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[2] = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[3] = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[4] = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[5] = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[6] = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[7] = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[8] = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
   and om**2 = omx**2 + omy**2 + omz**2
   the rotation matrix is determined by:
   omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
   omz = (q/m)*bz(x(t),y(t))*gami,
   where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
   position equations used are:
   x(t+dt) = x(t) + px(t+dt/2)*dtg
   y(t+dt) = y(t) + py(t+dt/2)*dtg
   where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
   pz(t+dt/2)*pz(t+dt/2))*ci*ci)
   fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   idimp = size of phase space = 5
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
#define N 4
   int noff, moff, npoff, npp, ipp, joff, nps;
   int i, j, k, m, nn, mm, nm, lxv;
   float qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg;
   float omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, vx, vy, vz;
   float sfxy[N*MXV*MYV], sbxy[N*MXV*MYV];
/* float sfxy[N*(mx+1)*(my+1)], sbxy[N*(mx+1)*(my+1)]; */
/* scratch arrays */
   int n[NPBLK];
   float s1[NPBLK*LVECT], s2[NPBLK*LVECT], t[NPBLK*2];
   double sum1, sum2;
   lxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
   ci2 = ci*ci;
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
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,noff,moff,npp,npoff,ipp,joff,nps,nn,mm,nm,x,y,vx,vy,vz, \
dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm, \
rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,dtg,sum1, \
sfxy,sbxy,n,s1,s2,t) \
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
            sfxy[N*(i+lxv*j)] = fxy[N*(i+noff+nxv*(j+moff))];
            sfxy[1+N*(i+lxv*j)] = fxy[1+N*(i+noff+nxv*(j+moff))];
            sfxy[2+N*(i+lxv*j)] = fxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[N*(i+lxv*j)] = bxy[N*(i+noff+nxv*(j+moff))];
            sbxy[1+N*(i+lxv*j)] = bxy[1+N*(i+noff+nxv*(j+moff))];
            sbxy[2+N*(i+lxv*j)] = bxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
      sum1 = 0.0;
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
            n[j] = N*(nn - noff + lxv*(mm - moff));
            amx = 1.0f - dxp;
            amy = 1.0f - dyp;
            s1[j] = amx*amy;
            s1[j+NPBLK] = dxp*amy;
            s1[j+2*NPBLK] = amx*dyp;
            s1[j+3*NPBLK] = dxp*dyp;
            t[j] = x;
            t[j+NPBLK] = y;
         }
/* find acceleration */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + N*(lxv - 2);
            dx = 0.0f;
            dy = 0.0f;
            dz = 0.0f;
            ox = 0.0f;
            oy = 0.0f;
            oz = 0.0f;
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               if (i > 1)
                  nn = mm;
               dx += sfxy[N*i+nn]*s1[j+NPBLK*i];
               dy += sfxy[1+N*i+nn]*s1[j+NPBLK*i];
               dz += sfxy[2+N*i+nn]*s1[j+NPBLK*i];
               ox += sbxy[N*i+nn]*s1[j+NPBLK*i];
               oy += sbxy[1+N*i+nn]*s1[j+NPBLK*i];
               oz += sbxy[2+N*i+nn]*s1[j+NPBLK*i];
            }
            s1[j] = dx;
            s1[j+NPBLK] = dy;
            s1[j+2*NPBLK] = dz;
            s2[j] = ox;
            s2[j+NPBLK] = oy;
            s2[j+2*NPBLK] = oz;
         }
/* new momentum */
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
/* calculate half impulse */
            dx = qtmh*s1[j];
            dy = qtmh*s1[j+NPBLK];
            dz = qtmh*s1[j+2*NPBLK];
/* half acceleration */
            acx = ppart[j+joff+2*nppmx+npoff] + dx;
            acy = ppart[j+joff+3*nppmx+npoff] + dy;
            acz = ppart[j+joff+4*nppmx+npoff] + dz;
/* find inverse gamma */
            p2 = acx*acx + acy*acy + acz*acz;
            gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* renormalize magnetic field */
            qtmg = qtmh*gami;
/* time-centered kinetic energy */
            sum1 += gami*p2/(1.0f + gami);
/* calculate cyclotron frequency */
            omxt = qtmg*s2[j];
            omyt = qtmg*s2[j+NPBLK];
            omzt = qtmg*s2[j+2*NPBLK];
/* calculate rotation matrix */
            omt = omxt*omxt + omyt*omyt + omzt*omzt;
            anorm = 2.0f/(1.0f + omt);
            omt = 0.5f*(1.0f - omt);
            rot4 = omxt*omyt;
            rot7 = omxt*omzt;
            rot8 = omyt*omzt;
            rot1 = omt + omxt*omxt;
            rot5 = omt + omyt*omyt;
            rot9 = omt + omzt*omzt;
            rot2 = omzt + rot4;
            rot4 -= omzt;
            rot3 = -omyt + rot7;
            rot7 += omyt;
            rot6 = omxt + rot8;
            rot8 -= omxt;
/* new velocity */
            vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
            vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
            vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* update inverse gamma */
            p2 = vx*vx + vy*vy + vz*vz;
            dtg = dtc/sqrtf(1.0f + p2*ci2);
/* new position */
            s1[j] = x + vx*dtg;
            s1[j+NPBLK] = y + vy*dtg;
            s2[j] = vx;
            s2[j+NPBLK] = vy;
            s2[j+2*NPBLK] = vz;
         }
/* check boundary conditions */
#pragma novector
         for (j = 0; j < NPBLK; j++) {
            dx = s1[j];
            dy = s1[j+NPBLK];
            vx = s2[j];
            vy = s2[j+NPBLK];
            vz = s2[j+2*NPBLK];
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
/* set new momentum */
            ppart[j+joff+2*nppmx+npoff] = vx;
            ppart[j+joff+3*nppmx+npoff] = vy;
            ppart[j+joff+4*nppmx+npoff] = vz;
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
         nm = N*(nn - noff + lxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + N;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += N*lxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + N;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + N;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += N*lxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + N;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
/* find inverse gamma */
         p2 = acx*acx + acy*acy + acz*acz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* renormalize magnetic field */
         qtmg = qtmh*gami;
/* time-centered kinetic energy */
         sum1 += gami*p2/(1.0f + gami);
/* calculate cyclotron frequency */
         omxt = qtmg*ox;
         omyt = qtmg*oy;
         omzt = qtmg*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0f/(1.0f + omt);
         omt = 0.5f*(1.0f - omt);
         rot4 = omxt*omyt;
         rot7 = omxt*omzt;
         rot8 = omyt*omzt;
         rot1 = omt + omxt*omxt;
         rot5 = omt + omyt*omyt;
         rot9 = omt + omzt*omzt;
         rot2 = omzt + rot4;
         rot4 -= omzt;
         rot3 = -omyt + rot7;
         rot7 += omyt;
         rot6 = omxt + rot8;
         rot8 -= omxt;
/* new velocity */
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* update inverse gamma */
         p2 = vx*vx + vy*vy + vz*vz;
         dtg = dtc/sqrtf(1.0f + p2*ci2);
/* new position */
         dx = x + vx*dtg;
         dy = y + vy*dtg;
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
/* set new momentum */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
         ppart[j+4*nppmx+npoff] = vz;
      }
      sum2 += sum1;
   }
/* normalize kinetic energy */
   *ek += sum2;
   return;
#undef N
#undef LVECT
#undef NPBLK
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cvgrbppushf23lt(float ppart[], float fxy[], float bxy[],
                     int kpic[], int ncl[], int ihole[], float qbm,
                     float dt, float dtc, float ci, float *ek, 
                     int idimp, int nppmx, int nx, int ny, int mx,
                     int my, int nxv, int nyv, int mx1, int mxy1,
                     int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   with periodic boundary conditions.
   Using the Boris Mover.
   also determines list of particles which are leaving this tile
   vectorizable/OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
   input: all except ncl, ihole, irc, output: ppart, ncl, ihole, irc, ek
   momentum equations used are:
   px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fz(x(t),y(t))*dt)
   where q/m is charge/mass, and the rotation matrix is given by:
      rot[0] = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[1] = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[2] = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[3] = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[4] = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[5] = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[6] = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[7] = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[8] = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
   and om**2 = omx**2 + omy**2 + omz**2
   the rotation matrix is determined by:
   omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
   omz = (q/m)*bz(x(t),y(t))*gami,
   where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
   position equations used are:
   x(t+dt) = x(t) + px(t+dt/2)*dtg
   y(t+dt) = y(t) + py(t+dt/2)*dtg
   where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
   pz(t+dt/2)*pz(t+dt/2))*ci*ci)
   fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   idimp = size of phase space = 5
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
#define N 4
   int noff, moff, npoff, npp, ipp, joff, nps;
   int i, j, k, m, ih, nh, nn, mm, nm, lxv;
   float qtmh, ci2, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz;
   float acx, acy, acz, p2, gami, qtmg, dtg, omxt, omyt, omzt, omt;
   float anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float x, y, vx, vy, vz;
   float sfxy[N*MXV*MYV], sbxy[N*MXV*MYV];
/* float sfxy[N*(mx+1)*(my+1)], sbxy[N*(mx+1)*(my+1)]; */
/* scratch arrays */
   int n[NPBLK];
   float s1[NPBLK*LVECT], s2[NPBLK*LVECT], t[NPBLK*2];
   double sum1, sum2;
   lxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
   ci2 = ci*ci;
   anx = (float) nx;
   any = (float) ny;
   sum2 = 0.0;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,noff,moff,npp,npoff,ipp,joff,nps,nn,mm,nm,ih,nh,x,y,vx, \
vy,vz,dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt, \
anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely, \
edgerx,edgery,p2,gami,qtmg,dtg,sum1,sfxy,sbxy,n,s1,s2,t) \
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
            sfxy[N*(i+lxv*j)] = fxy[N*(i+noff+nxv*(j+moff))];
            sfxy[1+N*(i+lxv*j)] = fxy[1+N*(i+noff+nxv*(j+moff))];
            sfxy[2+N*(i+lxv*j)] = fxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[N*(i+lxv*j)] = bxy[N*(i+noff+nxv*(j+moff))];
            sbxy[1+N*(i+lxv*j)] = bxy[1+N*(i+noff+nxv*(j+moff))];
            sbxy[2+N*(i+lxv*j)] = bxy[2+N*(i+noff+nxv*(j+moff))];
         }
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
      sum1 = 0.0;
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
            n[j] = N*(nn - noff + lxv*(mm - moff));
            amx = 1.0f - dxp;
            amy = 1.0f - dyp;
            s1[j] = amx*amy;
            s1[j+NPBLK] = dxp*amy;
            s1[j+2*NPBLK] = amx*dyp;
            s1[j+3*NPBLK] = dxp*dyp;
            t[j] = x;
            t[j+NPBLK] = y;
         }
/* find acceleration */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + N*(lxv - 2);
            dx = 0.0f;
            dy = 0.0f;
            dz = 0.0f;
            ox = 0.0f;
            oy = 0.0f;
            oz = 0.0f;
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               if (i > 1)
                  nn = mm;
               dx += sfxy[N*i+nn]*s1[j+NPBLK*i];
               dy += sfxy[1+N*i+nn]*s1[j+NPBLK*i];
               dz += sfxy[2+N*i+nn]*s1[j+NPBLK*i];
               ox += sbxy[N*i+nn]*s1[j+NPBLK*i];
               oy += sbxy[1+N*i+nn]*s1[j+NPBLK*i];
               oz += sbxy[2+N*i+nn]*s1[j+NPBLK*i];
            }
            s1[j] = dx;
            s1[j+NPBLK] = dy;
            s1[j+2*NPBLK] = dz;
            s2[j] = ox;
            s2[j+NPBLK] = oy;
            s2[j+2*NPBLK] = oz;
         }
/* new momentum */
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
/* calculate half impulse */
            dx = qtmh*s1[j];
            dy = qtmh*s1[j+NPBLK];
            dz = qtmh*s1[j+2*NPBLK];
/* half acceleration */
            acx = ppart[j+joff+2*nppmx+npoff] + dx;
            acy = ppart[j+joff+3*nppmx+npoff] + dy;
            acz = ppart[j+joff+4*nppmx+npoff] + dz;
/* find inverse gamma */
            p2 = acx*acx + acy*acy + acz*acz;
            gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* renormalize magnetic field */
            qtmg = qtmh*gami;
/* time-centered kinetic energy */
            sum1 += gami*p2/(1.0f + gami);
/* calculate cyclotron frequency */
            omxt = qtmg*s2[j];
            omyt = qtmg*s2[j+NPBLK];
            omzt = qtmg*s2[j+2*NPBLK];
/* calculate rotation matrix */
            omt = omxt*omxt + omyt*omyt + omzt*omzt;
            anorm = 2.0f/(1.0f + omt);
            omt = 0.5f*(1.0f - omt);
            rot4 = omxt*omyt;
            rot7 = omxt*omzt;
            rot8 = omyt*omzt;
            rot1 = omt + omxt*omxt;
            rot5 = omt + omyt*omyt;
            rot9 = omt + omzt*omzt;
            rot2 = omzt + rot4;
            rot4 -= omzt;
            rot3 = -omyt + rot7;
            rot7 += omyt;
            rot6 = omxt + rot8;
            rot8 -= omxt;
/* new momentum */
            vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
            vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
            vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* update inverse gamma */
            p2 = vx*vx + vy*vy + vz*vz;
            dtg = dtc/sqrtf(1.0f + p2*ci2);
/* new position */
            s1[j] = x + vx*dtg;
            s1[j+NPBLK] = y + vy*dtg;
            s2[j] = vx;
            s2[j+NPBLK] = vy;
            s2[j+2*NPBLK] = vz;
         }
/* check boundary conditions */
#pragma novector
         for (j = 0; j < NPBLK; j++) {
            dx = s1[j];
            dy = s1[j+NPBLK];
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
/* set new momentum */
            ppart[j+joff+2*nppmx+npoff] = s2[j];
            ppart[j+joff+3*nppmx+npoff] = s2[j+NPBLK];
            ppart[j+joff+4*nppmx+npoff] = s2[j+2*NPBLK];
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
         nm = N*(nn - noff + lxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + N;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += N*lxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + N;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + N;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += N*lxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + N;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
/* find inverse gamma */
         p2 = acx*acx + acy*acy + acz*acz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* renormalize magnetic field */
         qtmg = qtmh*gami;
/* time-centered kinetic energy */
         sum1 += gami*p2/(1.0f + gami);
/* calculate cyclotron frequency */
         omxt = qtmg*ox;
         omyt = qtmg*oy;
         omzt = qtmg*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0f/(1.0f + omt);
         omt = 0.5f*(1.0f - omt);
         rot4 = omxt*omyt;
         rot7 = omxt*omzt;
         rot8 = omyt*omzt;
         rot1 = omt + omxt*omxt;
         rot5 = omt + omyt*omyt;
         rot9 = omt + omzt*omzt;
         rot2 = omzt + rot4;
         rot4 -= omzt;
         rot3 = -omyt + rot7;
         rot7 += omyt;
         rot6 = omxt + rot8;
         rot8 -= omxt;
/* new momentum */
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* update inverse gamma */
         p2 = vx*vx + vy*vy + vz*vz;
         dtg = dtc/sqrtf(1.0f + p2*ci2);
/* new position */
         dx = x + vx*dtg;
         dy = y + vy*dtg;
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
/* set new momentum */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
         ppart[j+4*nppmx+npoff] = vz;
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
   *ek += sum2;
   return;
#undef N
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
void cgjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                 float dt, int nppmx, int idimp, int nx, int ny, int mx,
                 int my, int nxv, int nyv, int mx1, int mxy1,
                 int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   41 flops/particle, 17 loads, 14 stores
   input: all, output: ppart, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*vi, where i = x,y,z
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
#define N 4
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nn, mm;
   float edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz;
   float scu[N*MXV*MYV];
/* float scu[N*(mx+1)*(my+1)]; */
   mxv = mx + 1;
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
vy,vz,scu)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < N*mxv*(my+1); j++) {
         scu[j] = 0.0f;
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
         nn = N*(nn - noff + mxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ppart[j+2*nppmx+npoff];
         vy = ppart[j+3*nppmx+npoff];
         vz = ppart[j+4*nppmx+npoff];
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += N*mxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+3*nppmx+npoff] = -vy;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -vx;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[N*(i+noff+nxv*(j+moff))] += scu[N*(i+mxv*j)];
            cu[1+N*(i+noff+nxv*(j+moff))] += scu[1+N*(i+mxv*j)];
            cu[2+N*(i+noff+nxv*(j+moff))] += scu[2+N*(i+mxv*j)];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[N*(i+noff+nxv*moff)] += scu[N*i];
#pragma omp atomic
         cu[1+N*(i+noff+nxv*moff)] += scu[1+N*i];
#pragma omp atomic
         cu[2+N*(i+noff+nxv*moff)] += scu[2+N*i];
         if (mm > my) {
#pragma omp atomic
            cu[N*(i+noff+nxv*(mm+moff-1))] += scu[N*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[1+N*(i+noff+nxv*(mm+moff-1))] += scu[1+N*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[2+N*(i+noff+nxv*(mm+moff-1))] += scu[2+N*(i+mxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[N*(noff+nxv*(j+moff))] += scu[N*mxv*j];
#pragma omp atomic
         cu[1+N*(noff+nxv*(j+moff))] += scu[1+N*mxv*j];
#pragma omp atomic
         cu[2+N*(noff+nxv*(j+moff))] += scu[2+N*mxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[N*(nn+noff-1+nxv*(j+moff))] += scu[N*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[1+N*(nn+noff-1+nxv*(j+moff))] += scu[1+N*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[2+N*(nn+noff-1+nxv*(j+moff))] += scu[2+N*((nn-1)+mxv*j)];
         }
      }
   }
   return;
#undef N
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                  int ihole[], float qm, float dt, int nppmx, int idimp,
                  int nx, int ny, int mx, int my, int nxv, int nyv,
                  int mx1, int mxy1, int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   41 flops/particle, 17 loads, 14 stores
   input: all except ncl, ihole, irc,
   output: ppart, cu, ncl, ihole, irc
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*vi, where i = x,y,z
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             33
#define MYV             33
#define N 4
   int noff, moff, npoff, npp;
   int i, j, k, ih, nh, nn, mm, mxv;
   float dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float scu[N*MXV*MYV];
/* float scu[N*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   anx = (float) nx;
   any = (float) ny;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,ih,nh,x,y,dxp,dyp,amx,amy,dx, \
dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,scu)
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
/* zero out local accumulator */
      for (j = 0; j < N*mxv*(my+1); j++) {
         scu[j] = 0.0f;
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
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
         nn = N*(nn - noff + mxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ppart[j+2*nppmx+npoff];
         vy = ppart[j+3*nppmx+npoff];
         vz = ppart[j+4*nppmx+npoff];
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += N*mxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
/* advance position half a time-step */
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
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[N*(i+noff+nxv*(j+moff))] += scu[N*(i+mxv*j)];
            cu[1+N*(i+noff+nxv*(j+moff))] += scu[1+N*(i+mxv*j)];
            cu[2+N*(i+noff+nxv*(j+moff))] += scu[2+N*(i+mxv*j)];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[N*(i+noff+nxv*moff)] += scu[N*i];
#pragma omp atomic
         cu[1+N*(i+noff+nxv*moff)] += scu[1+N*i];
#pragma omp atomic
         cu[2+N*(i+noff+nxv*moff)] += scu[2+N*i];
         if (mm > my) {
#pragma omp atomic
            cu[N*(i+noff+nxv*(mm+moff-1))] += scu[N*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[1+N*(i+noff+nxv*(mm+moff-1))] += scu[1+N*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[2+N*(i+noff+nxv*(mm+moff-1))] += scu[2+N*(i+mxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[N*(noff+nxv*(j+moff))] += scu[N*mxv*j];
#pragma omp atomic
         cu[1+N*(noff+nxv*(j+moff))] += scu[1+N*mxv*j];
#pragma omp atomic
         cu[2+N*(noff+nxv*(j+moff))] += scu[2+N*mxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[N*(nn+noff-1+nxv*(j+moff))] += scu[N*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[1+N*(nn+noff-1+nxv*(j+moff))] += scu[1+N*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[2+N*(nn+noff-1+nxv*(j+moff))] += scu[2+N*((nn-1)+mxv*j)];
         }
      }
/* set error and end of file flag */
/* ihole overflow */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*k] = ih;
   }
   return;
#undef N
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgrjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                  float dt, float ci, int nppmx, int idimp, int nx,
                  int ny, int mx, int my, int nxv, int nyv, int mx1,
                  int mxy1, int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
   input: all, output: ppart, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*pi*gami, where i = x,y,z
   where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
#define N 4
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nn, mm;
   float ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami;
   float scu[N*MXV*MYV];
/* float scu[N*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   ci2 = ci*ci;
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
vy,vz,ux,uy,uz,p2,gami,scu)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < N*mxv*(my+1); j++) {
         scu[j] = 0.0f;
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
/* find inverse gamma */
         ux = ppart[j+2*nppmx+npoff];
         uy = ppart[j+3*nppmx+npoff];
         uz = ppart[j+4*nppmx+npoff];
         p2 = ux*ux + uy*uy + uz*uz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* calculate weights */
         nn = N*(nn - noff + mxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ux*gami;
         vy = uy*gami;
         vz = uz*gami;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += N*mxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -ux;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+3*nppmx+npoff] = -uy;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -ux;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[N*(i+noff+nxv*(j+moff))] += scu[N*(i+mxv*j)];
            cu[1+N*(i+noff+nxv*(j+moff))] += scu[1+N*(i+mxv*j)];
            cu[2+N*(i+noff+nxv*(j+moff))] += scu[2+N*(i+mxv*j)];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[N*(i+noff+nxv*moff)] += scu[N*i];
#pragma omp atomic
         cu[1+N*(i+noff+nxv*moff)] += scu[1+N*i];
#pragma omp atomic
         cu[2+N*(i+noff+nxv*moff)] += scu[2+N*i];
         if (mm > my) {
#pragma omp atomic
            cu[N*(i+noff+nxv*(mm+moff-1))] += scu[N*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[1+N*(i+noff+nxv*(mm+moff-1))] += scu[1+N*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[2+N*(i+noff+nxv*(mm+moff-1))] += scu[2+N*(i+mxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[N*(noff+nxv*(j+moff))] += scu[N*mxv*j];
#pragma omp atomic
         cu[1+N*(noff+nxv*(j+moff))] += scu[1+N*mxv*j];
#pragma omp atomic
         cu[2+N*(noff+nxv*(j+moff))] += scu[2+N*mxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[N*(nn+noff-1+nxv*(j+moff))] += scu[N*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[1+N*(nn+noff-1+nxv*(j+moff))] += scu[1+N*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[2+N*(nn+noff-1+nxv*(j+moff))] += scu[2+N*((nn-1)+mxv*j)];
         }
      }
   }
   return;
#undef N
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgrjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                   int ihole[], float qm, float dt, float ci, int nppmx,
                   int idimp, int nx, int ny, int mx, int my, int nxv,
                   int nyv, int mx1, int mxy1, int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation for relativistic particles
   in addition, particle positions are advanced a half time-step
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
   input: all except ncl, ihole, irc,
   output: ppart, cu, ncl, ihole, irc
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*pi*gami, where i = x,y,z
   where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             33
#define MYV             33
#define N 4
   int noff, moff, npoff, npp;
   int i, j, k, ih, nh, nn, mm, mxv;
   float ci2, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float scu[N*MXV*MYV];
/* float scu[N*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   ci2 = ci*ci;
   anx = (float) nx;
   any = (float) ny;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,ih,nh,x,y,dxp,dyp,amx,amy,dx, \
dy,vx,vy,vz,ux,uy,uz,edgelx,edgely,edgerx,edgery,p2,gami,scu)
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
/* zero out local accumulator */
      for (j = 0; j < N*mxv*(my+1); j++) {
         scu[j] = 0.0f;
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
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
/* find inverse gamma */
         ux = ppart[j+2*nppmx+npoff];
         uy = ppart[j+3*nppmx+npoff];
         uz = ppart[j+4*nppmx+npoff];
         p2 = ux*ux + uy*uy + uz*uz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* calculate weights */
         nn = N*(nn - noff + mxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ux*gami;
         vy = uy*gami;
         vz = uz*gami;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += N*mxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
/* advance position half a time-step */
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
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[N*(i+noff+nxv*(j+moff))] += scu[N*(i+mxv*j)];
            cu[1+N*(i+noff+nxv*(j+moff))] += scu[1+N*(i+mxv*j)];
            cu[2+N*(i+noff+nxv*(j+moff))] += scu[2+N*(i+mxv*j)];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[N*(i+noff+nxv*moff)] += scu[N*i];
#pragma omp atomic
         cu[1+N*(i+noff+nxv*moff)] += scu[1+N*i];
#pragma omp atomic
         cu[2+N*(i+noff+nxv*moff)] += scu[2+N*i];
         if (mm > my) {
#pragma omp atomic
            cu[N*(i+noff+nxv*(mm+moff-1))] += scu[N*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[1+N*(i+noff+nxv*(mm+moff-1))] += scu[1+N*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[2+N*(i+noff+nxv*(mm+moff-1))] += scu[2+N*(i+mxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[N*(noff+nxv*(j+moff))] += scu[N*mxv*j];
#pragma omp atomic
         cu[1+N*(noff+nxv*(j+moff))] += scu[1+N*mxv*j];
#pragma omp atomic
         cu[2+N*(noff+nxv*(j+moff))] += scu[2+N*mxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[N*(nn+noff-1+nxv*(j+moff))] += scu[N*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[1+N*(nn+noff-1+nxv*(j+moff))] += scu[1+N*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[2+N*(nn+noff-1+nxv*(j+moff))] += scu[2+N*((nn-1)+mxv*j)];
         }
      }
/* set error and end of file flag */
/* ihole overflow */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*k] = ih;
   }
   return;
#undef N
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cvgjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                  float dt, int nppmx, int idimp, int nx, int ny,
                  int mx, int my, int nxv, int nyv, int mx1, int mxy1,
                  int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   vectorizable/OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   41 flops/particle, 17 loads, 14 stores
   input: all, output: ppart, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*vi, where i = x,y,z
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
#define NPBLK             32
#define LVECT             4
#define N 4
   int noff, moff, npoff, npp, lxv;
   int i, j, k, m, ipp, joff, nps, nn, mm;
   float edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz;
   float scu[N*MXV*MYV];
/* float scu[N*(mx+1)*(my+1)]; */
/* scratch arrays */
   int n[NPBLK];
   float s1[NPBLK*LVECT], s2[NPBLK*LVECT], t[NPBLK*2];
   lxv = mx + 1;
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
private(i,j,k,m,noff,moff,npp,npoff,ipp,joff,nps,nn,mm,x,y,dxp,dyp,amx, \
amy,dx,dy,vx,vy,vz,scu,n,s1,s2,t)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < N*lxv*(my+1); j++) {
         scu[j] = 0.0f;
      }
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
            n[j] = N*(nn - noff + lxv*(mm - moff));
            amx = qm - dxp;
            amy = 1.0f - dyp;
            s1[j] = amx*amy;
            s1[j+NPBLK] = dxp*amy;
            s1[j+2*NPBLK] = amx*dyp;
            s1[j+3*NPBLK] = dxp*dyp;
            t[j] = x;
            t[j+NPBLK] = y;
            s2[j] = ppart[j+joff+2*nppmx+npoff];
            s2[j+NPBLK] = ppart[j+joff+3*nppmx+npoff];
            s2[j+2*NPBLK] = ppart[j+joff+4*nppmx+npoff];
         }
/* deposit current */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + N*(lxv - 2);
            vx = s2[j];
            vy = s2[j+NPBLK];
            vz = s2[j+2*NPBLK];
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               if (i > 1)
                  nn = mm;
               scu[N*i+nn] += vx*s1[j+NPBLK*i];
               scu[1+N*i+nn] += vy*s1[j+NPBLK*i];
               scu[2+N*i+nn] += vz*s1[j+NPBLK*i];
            }
         }
/* advance position half a time-step */
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
            vx = s2[j];
            vy = s2[j+NPBLK];
            dx = x + vx*dt;
            dy = y + vy*dt;
/* reflecting boundary conditions */
            if (ipbc==2) {
               if ((dx < edgelx) || (dx >= edgerx)) {
                  dx = x;
                  ppart[j+joff+2*nppmx+npoff] = -vx;
               }
               if ((dy < edgely) || (dy >= edgery)) {
                  dy = y;
                  ppart[j+joff+3*nppmx+npoff] = -vy;
               }
            }
/* mixed reflecting/periodic boundary conditions */
            else if (ipbc==3) {
               if ((dx < edgelx) || (dx >= edgerx)) {
                  dx = x;
                  ppart[j+joff+2*nppmx+npoff] = -vx;
               }
            }
/* set new position */
            ppart[j+joff+npoff] = dx;
            ppart[j+joff+nppmx+npoff] = dy;
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
         nn = N*(nn - noff + lxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ppart[j+2*nppmx+npoff];
         vy = ppart[j+3*nppmx+npoff];
         vz = ppart[j+4*nppmx+npoff];
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += N*lxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+3*nppmx+npoff] = -vy;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -vx;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[N*(i+noff+nxv*(j+moff))] += scu[N*(i+lxv*j)];
            cu[1+N*(i+noff+nxv*(j+moff))] += scu[1+N*(i+lxv*j)];
            cu[2+N*(i+noff+nxv*(j+moff))] += scu[2+N*(i+lxv*j)];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[N*(i+noff+nxv*moff)] += scu[N*i];
#pragma omp atomic
         cu[1+N*(i+noff+nxv*moff)] += scu[1+N*i];
#pragma omp atomic
         cu[2+N*(i+noff+nxv*moff)] += scu[2+N*i];
         if (mm > my) {
#pragma omp atomic
            cu[N*(i+noff+nxv*(mm+moff-1))] += scu[N*(i+lxv*(mm-1))];
#pragma omp atomic
            cu[1+N*(i+noff+nxv*(mm+moff-1))] += scu[1+N*(i+lxv*(mm-1))];
#pragma omp atomic
            cu[2+N*(i+noff+nxv*(mm+moff-1))] += scu[2+N*(i+lxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[N*(noff+nxv*(j+moff))] += scu[N*lxv*j];
#pragma omp atomic
         cu[1+N*(noff+nxv*(j+moff))] += scu[1+N*lxv*j];
#pragma omp atomic
         cu[2+N*(noff+nxv*(j+moff))] += scu[2+N*lxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[N*(nn+noff-1+nxv*(j+moff))] += scu[N*((nn-1)+lxv*j)];
#pragma omp atomic
            cu[1+N*(nn+noff-1+nxv*(j+moff))] += scu[1+N*((nn-1)+lxv*j)];
#pragma omp atomic
            cu[2+N*(nn+noff-1+nxv*(j+moff))] += scu[2+N*((nn-1)+lxv*j)];
         }
      }
   }
   return;
#undef N
#undef LVECT
#undef NPBLK
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cvgjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                   int ihole[], float qm, float dt, int nppmx,
                   int idimp, int nx, int ny, int mx, int my, int nxv,
                   int nyv, int mx1, int mxy1, int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   vectorizable/OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   41 flops/particle, 17 loads, 14 stores
   input: all except ncl, ihole, irc,
   output: ppart, cu, ncl, ihole, irc
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*vi, where i = x,y,z
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
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
#define N 4
   int noff, moff, npoff, npp, lxv;
   int i, j, k, m, ih, nh, ipp, joff, nps, nn, mm;
   float dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float scu[N*MXV*MYV];
/* float scu[N*(mx+1)*(my+1)]; */
/* scratch arrays */
   int n[NPBLK];
   float s1[NPBLK*LVECT], s2[NPBLK*LVECT], t[NPBLK*2];
   lxv = mx + 1;
   anx = (float) nx;
   any = (float) ny;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,noff,moff,npp,npoff,ipp,joff,nps,nn,mm,ih,nh,x,y,dxp, \
dyp,amx,amy,dx,dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,scu,n,s1,s2,t)
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
/* zero out local accumulator */
      for (j = 0; j < N*lxv*(my+1); j++) {
         scu[j] = 0.0f;
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
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
            n[j] = N*(nn - noff + lxv*(mm - moff));
            amx = qm - dxp;
            amy = 1.0f - dyp;
            s1[j] = amx*amy;
            s1[j+NPBLK] = dxp*amy;
            s1[j+2*NPBLK] = amx*dyp;
            s1[j+3*NPBLK] = dxp*dyp;
            t[j] = x;
            t[j+NPBLK] = y;
            s2[j] = ppart[j+joff+2*nppmx+npoff];
            s2[j+NPBLK] = ppart[j+joff+3*nppmx+npoff];
            s2[j+2*NPBLK] = ppart[j+joff+4*nppmx+npoff];
         }
/* deposit current */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + N*(lxv - 2);
            vx = s2[j];
            vy = s2[j+NPBLK];
            vz = s2[j+2*NPBLK];
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               if (i > 1)
                  nn = mm;
               scu[N*i+nn] += vx*s1[j+NPBLK*i];
               scu[1+N*i+nn] += vy*s1[j+NPBLK*i];
               scu[2+N*i+nn] += vz*s1[j+NPBLK*i];
            }
         }
/* advance position half a time-step */
         for (j = 0; j < NPBLK; j++) {
            dx = t[j] + s2[j]*dt;
            dy = t[j+NPBLK] + s2[j+NPBLK]*dt;
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
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         nn = N*(nn - noff + lxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ppart[j+2*nppmx+npoff];
         vy = ppart[j+3*nppmx+npoff];
         vz = ppart[j+4*nppmx+npoff];
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += N*lxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
/* advance position half a time-step */
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
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[N*(i+noff+nxv*(j+moff))] += scu[N*(i+lxv*j)];
            cu[1+N*(i+noff+nxv*(j+moff))] += scu[1+N*(i+lxv*j)];
            cu[2+N*(i+noff+nxv*(j+moff))] += scu[2+N*(i+lxv*j)];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[N*(i+noff+nxv*moff)] += scu[N*i];
#pragma omp atomic
         cu[1+N*(i+noff+nxv*moff)] += scu[1+N*i];
#pragma omp atomic
         cu[2+N*(i+noff+nxv*moff)] += scu[2+N*i];
         if (mm > my) {
#pragma omp atomic
            cu[N*(i+noff+nxv*(mm+moff-1))] += scu[N*(i+lxv*(mm-1))];
#pragma omp atomic
            cu[1+N*(i+noff+nxv*(mm+moff-1))] += scu[1+N*(i+lxv*(mm-1))];
#pragma omp atomic
            cu[2+N*(i+noff+nxv*(mm+moff-1))] += scu[2+N*(i+lxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[N*(noff+nxv*(j+moff))] += scu[N*lxv*j];
#pragma omp atomic
         cu[1+N*(noff+nxv*(j+moff))] += scu[1+N*lxv*j];
#pragma omp atomic
         cu[2+N*(noff+nxv*(j+moff))] += scu[2+N*lxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[N*(nn+noff-1+nxv*(j+moff))] += scu[N*((nn-1)+lxv*j)];
#pragma omp atomic
            cu[1+N*(nn+noff-1+nxv*(j+moff))] += scu[1+N*((nn-1)+lxv*j)];
#pragma omp atomic
            cu[2+N*(nn+noff-1+nxv*(j+moff))] += scu[2+N*((nn-1)+lxv*j)];
         }
      }
/* set error and end of file flag */
/* ihole overflow */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*k] = ih;
   }
   return;
#undef N
#undef LVECT
#undef NPBLK
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cvgrjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                   float dt, float ci, int nppmx, int idimp, int nx,
                   int ny, int mx, int my, int nxv, int nyv, int mx1,
                   int mxy1, int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   vectorizable/OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
   input: all, output: ppart, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*pi*gami, where i = x,y,z
   where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
#define NPBLK             32
#define LVECT             4
#define N 4
   int noff, moff, npoff, npp, lxv;
   int i, j, k, m, ipp, joff, nps, nn, mm;
   float ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami;
   float scu[N*MXV*MYV];
/* float scu[N*(mx+1)*(my+1)]; */
/* scratch arrays */
   int n[NPBLK];
   float s1[NPBLK*LVECT], s2[NPBLK*LVECT], t[NPBLK*4];
   lxv = mx + 1;
   ci2 = ci*ci;
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
private(i,j,k,m,noff,moff,npp,npoff,ipp,joff,nps,nn,mm,x,y,dxp,dyp,amx, \
amy,dx,dy,vx,vy,vz,ux,uy,uz,p2,gami,scu,n,s1,s2,t)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < N*lxv*(my+1); j++) {
         scu[j] = 0.0f;
      }
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
            n[j] = N*(nn - noff + lxv*(mm - moff));
            amx = qm - dxp;
            amy = 1.0f - dyp;
            s1[j] = amx*amy;
            s1[j+NPBLK] = dxp*amy;
            s1[j+2*NPBLK] = amx*dyp;
            s1[j+3*NPBLK] = dxp*dyp;
            t[j] = x;
            t[j+NPBLK] = y;
/* find inverse gamma */
            ux = ppart[j+joff+2*nppmx+npoff];
            uy = ppart[j+joff+3*nppmx+npoff];
            uz = ppart[j+joff+4*nppmx+npoff];
            p2 = ux*ux + uy*uy + uz*uz;
            gami = 1.0f/sqrtf(1.0f + p2*ci2);
            s2[j] = ux*gami;
            s2[j+NPBLK] = uy*gami;
            s2[j+2*NPBLK] = uz*gami;
            t[j+2*NPBLK] = ux;
            t[j+3*NPBLK] = uy;
         }
/* deposit current */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + N*(lxv - 2);
            vx = s2[j];
            vy = s2[j+NPBLK];
            vz = s2[j+2*NPBLK];
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               if (i > 1)
                  nn = mm;
               scu[N*i+nn] += vx*s1[j+NPBLK*i];
               scu[1+N*i+nn] += vy*s1[j+NPBLK*i];
               scu[2+N*i+nn] += vz*s1[j+NPBLK*i];
            }
         }
/* advance position half a time-step */
         for (j = 0; j < NPBLK; j++) {
            x = t[j];
            y = t[j+NPBLK];
            vx = s2[j];
            vy = s2[j+NPBLK];
            ux = t[j+2*NPBLK];
            uy = t[j+3*NPBLK];
            dx = x + vx*dt;
            dy = y + vy*dt;
/* reflecting boundary conditions */
            if (ipbc==2) {
               if ((dx < edgelx) || (dx >= edgerx)) {
                  dx = x;
                  ppart[j+joff+2*nppmx+npoff] = -ux;
               }
               if ((dy < edgely) || (dy >= edgery)) {
                  dy = y;
                  ppart[j+joff+3*nppmx+npoff] = -uy;
               }
            }
/* mixed reflecting/periodic boundary conditions */
            else if (ipbc==3) {
               if ((dx < edgelx) || (dx >= edgerx)) {
                  dx = x;
                  ppart[j+joff+2*nppmx+npoff] = -ux;
               }
            }
/* set new position */
            ppart[j+joff+npoff] = dx;
            ppart[j+joff+nppmx+npoff] = dy;
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
/* find inverse gamma */
         ux = ppart[j+2*nppmx+npoff];
         uy = ppart[j+3*nppmx+npoff];
         uz = ppart[j+4*nppmx+npoff];
         p2 = ux*ux + uy*uy + uz*uz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* calculate weights */
         nn = N*(nn - noff + lxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ux*gami;
         vy = uy*gami;
         vz = uz*gami;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += N*lxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -ux;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+3*nppmx+npoff] = -uy;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -ux;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[N*(i+noff+nxv*(j+moff))] += scu[N*(i+lxv*j)];
            cu[1+N*(i+noff+nxv*(j+moff))] += scu[1+N*(i+lxv*j)];
            cu[2+N*(i+noff+nxv*(j+moff))] += scu[2+N*(i+lxv*j)];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[N*(i+noff+nxv*moff)] += scu[N*i];
#pragma omp atomic
         cu[1+N*(i+noff+nxv*moff)] += scu[1+N*i];
#pragma omp atomic
         cu[2+N*(i+noff+nxv*moff)] += scu[2+N*i];
         if (mm > my) {
#pragma omp atomic
            cu[N*(i+noff+nxv*(mm+moff-1))] += scu[N*(i+lxv*(mm-1))];
#pragma omp atomic
            cu[1+N*(i+noff+nxv*(mm+moff-1))] += scu[1+N*(i+lxv*(mm-1))];
#pragma omp atomic
            cu[2+N*(i+noff+nxv*(mm+moff-1))] += scu[2+N*(i+lxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[N*(noff+nxv*(j+moff))] += scu[N*lxv*j];
#pragma omp atomic
         cu[1+N*(noff+nxv*(j+moff))] += scu[1+N*lxv*j];
#pragma omp atomic
         cu[2+N*(noff+nxv*(j+moff))] += scu[2+N*lxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[N*(nn+noff-1+nxv*(j+moff))] += scu[N*((nn-1)+lxv*j)];
#pragma omp atomic
            cu[1+N*(nn+noff-1+nxv*(j+moff))] += scu[1+N*((nn-1)+lxv*j)];
#pragma omp atomic
            cu[2+N*(nn+noff-1+nxv*(j+moff))] += scu[2+N*((nn-1)+lxv*j)];
         }
      }
   }
   return;
#undef N
#undef LVECT
#undef NPBLK
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cvgrjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                    int ihole[], float qm, float dt, float ci,
                    int nppmx, int idimp, int nx, int ny, int mx,
                    int my, int nxv, int nyv, int mx1, int mxy1,
                    int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation for relativistic particles
   in addition, particle positions are advanced a half time-step
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   vectorizable/OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
   input: all except ncl, ihole, irc,
   output: ppart, cu, ncl, ihole, irc
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*pi*gami, where i = x,y,z
   where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
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
#define N 4
   int noff, moff, npoff, npp, lxv;
   int i, j, k, m, ih, nh, ipp, joff, nps, nn, mm;
   float ci2, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz, p2, gami;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float scu[N*MXV*MYV];
/* float scu[N*(mx+1)*(my+1)]; */
/* scratch arrays */
   int n[NPBLK];
   float s1[NPBLK*LVECT], s2[NPBLK*LVECT], t[NPBLK*2];
   lxv = mx + 1;
   ci2 = ci*ci;
   anx = (float) nx;
   any = (float) ny;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,m,noff,moff,npp,npoff,ipp,joff,nps,nn,mm,ih,nh,x,y,dxp, \
dyp,amx,amy,dx,dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,p2,gami,scu,n, \
s1,s2,t)
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
/* zero out local accumulator */
      for (j = 0; j < N*lxv*(my+1); j++) {
         scu[j] = 0.0f;
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
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
            n[j] = N*(nn - noff + lxv*(mm - moff));
            amx = qm - dxp;
            amy = 1.0f - dyp;
            s1[j] = amx*amy;
            s1[j+NPBLK] = dxp*amy;
            s1[j+2*NPBLK] = amx*dyp;
            s1[j+3*NPBLK] = dxp*dyp;
            t[j] = x;
            t[j+NPBLK] = y;
/* find inverse gamma */
            vx = ppart[j+joff+2*nppmx+npoff];
            vy = ppart[j+joff+3*nppmx+npoff];
            vz = ppart[j+joff+4*nppmx+npoff];
            p2 = vx*vx + vy*vy + vz*vz;
            gami = 1.0f/sqrtf(1.0f + p2*ci2);
            s2[j] = vx*gami;
            s2[j+NPBLK] = vy*gami;
            s2[j+2*NPBLK] = vz*gami;
         }
/* deposit current */
         for (j = 0; j < NPBLK; j++) {
            nn = n[j];
            mm = nn + N*(lxv - 2);
            vx = s2[j];
            vy = s2[j+NPBLK];
            vz = s2[j+2*NPBLK];
#pragma ivdep
            for (i = 0; i < LVECT; i++) {
               if (i > 1)
                  nn = mm;
               scu[N*i+nn] += vx*s1[j+NPBLK*i];
               scu[1+N*i+nn] += vy*s1[j+NPBLK*i];
               scu[2+N*i+nn] += vz*s1[j+NPBLK*i];
            }
         }
/* advance position half a time-step */
         for (j = 0; j < NPBLK; j++) {
            dx = t[j] + s2[j]*dt;
            dy = t[j+NPBLK] + s2[j+NPBLK]*dt;
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
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
/* find inverse gamma */
         vx = ppart[j+2*nppmx+npoff];
         vy = ppart[j+3*nppmx+npoff];
         vz = ppart[j+4*nppmx+npoff];
         p2 = vx*vx + vy*vy + vz*vz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* calculate weights */
         nn = N*(nn - noff + lxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx *= gami;
         vy *= gami;
         vz *= gami;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += N*lxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + N;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
/* advance position half a time-step */
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
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[N*(i+noff+nxv*(j+moff))] += scu[N*(i+lxv*j)];
            cu[1+N*(i+noff+nxv*(j+moff))] += scu[1+N*(i+lxv*j)];
            cu[2+N*(i+noff+nxv*(j+moff))] += scu[2+N*(i+lxv*j)];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[N*(i+noff+nxv*moff)] += scu[N*i];
#pragma omp atomic
         cu[1+N*(i+noff+nxv*moff)] += scu[1+N*i];
#pragma omp atomic
         cu[2+N*(i+noff+nxv*moff)] += scu[2+N*i];
         if (mm > my) {
#pragma omp atomic
            cu[N*(i+noff+nxv*(mm+moff-1))] += scu[N*(i+lxv*(mm-1))];
#pragma omp atomic
            cu[1+N*(i+noff+nxv*(mm+moff-1))] += scu[1+N*(i+lxv*(mm-1))];
#pragma omp atomic
            cu[2+N*(i+noff+nxv*(mm+moff-1))] += scu[2+N*(i+lxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[N*(noff+nxv*(j+moff))] += scu[N*lxv*j];
#pragma omp atomic
         cu[1+N*(noff+nxv*(j+moff))] += scu[1+N*lxv*j];
#pragma omp atomic
         cu[2+N*(noff+nxv*(j+moff))] += scu[2+N*lxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[N*(nn+noff-1+nxv*(j+moff))] += scu[N*((nn-1)+lxv*j)];
#pragma omp atomic
            cu[1+N*(nn+noff-1+nxv*(j+moff))] += scu[1+N*((nn-1)+lxv*j)];
#pragma omp atomic
            cu[2+N*(nn+noff-1+nxv*(j+moff))] += scu[2+N*((nn-1)+lxv*j)];
         }
      }
/* set error and end of file flag */
/* ihole overflow */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*k] = ih;
   }
   return;
#undef N
#undef LVECT
#undef NPBLK
#undef MXV
#undef MYV
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
void cbguard2l(float bxy[], int nx, int ny, int nxe, int nye) {
/* replicate extended periodic vector field bxy
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nxe = second dimension of field arrays, must be >= ny+1
local data                                                 */
#define N 4
   int j, k, kk;
/* copy edges of extended field */
   for (k = 0; k < ny; k++) {
      kk = N*nxe*k;
      bxy[N*nx+kk] = bxy[kk];
      bxy[1+N*nx+kk] = bxy[1+kk];
      bxy[2+N*nx+kk] = bxy[2+kk];
   }
   kk = N*nxe*ny;
   for (j = 0; j < nx; j++) {
      bxy[N*j+kk] = bxy[N*j];
      bxy[1+N*j+kk] = bxy[1+N*j];
      bxy[2+N*j+kk] = bxy[2+N*j];
   }
   bxy[N*nx+kk] = bxy[0];
   bxy[1+N*nx+kk] = bxy[1];
   bxy[2+N*nx+kk] = bxy[2];
   return;
#undef N
}

/*--------------------------------------------------------------------*/
void cacguard2l(float cu[], int nx, int ny, int nxe, int nye) {
/* accumulate extended periodic vector field cu
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nxe = second dimension of field arrays, must be >= ny+1
local data                                                 */
#define N 4
   int j, k, kk;
/* accumulate edges of extended field */
   for (k = 0; k < ny; k++) {
      kk = N*nxe*k;
      cu[kk] += cu[N*nx+kk];
      cu[1+kk] += cu[1+N*nx+kk];
      cu[2+kk] += cu[2+N*nx+kk];
      cu[N*nx+kk] = 0.0;
      cu[1+N*nx+kk] = 0.0;
      cu[2+N*nx+kk] = 0.0;
   }
   kk = N*nxe*ny;
   for (j = 0; j < nx; j++) {
      cu[N*j] += cu[N*j+kk];
      cu[1+N*j] += cu[1+N*j+kk];
      cu[2+N*j] += cu[2+N*j+kk];
      cu[N*j+kk] = 0.0;
      cu[1+N*j+kk] = 0.0;
      cu[2+N*j+kk] = 0.0;
   }
   cu[0] += cu[N*nx+kk];
   cu[1] += cu[1+N*nx+kk];
   cu[2] += cu[2+N*nx+kk];
   cu[N*nx+kk] = 0.0;
   cu[1+N*nx+kk] = 0.0;
   cu[2+N*nx+kk] = 0.0;
   return;
#undef N
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
void cvmpois23(float complex q[], float complex fxy[], int isign,
               float complex ffc[], float ax, float ay, float affp,
               float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
               int nyhd) {
/* this subroutine solves 2-1/2d poisson's equation in fourier space for
   force/charge (or convolution of electric field over particle shape)
   with periodic boundary conditions.  Zeros out z component.
   for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
   for isign /= 0, input: q,ffc,isign,nx,ny,nxvh,nyhd, output: fxy,we
   approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1
   equation used is:
   fx[ky][kx] = -sqrt(-1)*kx*g[ky][kx]*s[ky][kx]*q[ky][kx],
   fy[ky][kx] = -sqrt(-1)*ky*g[ky][kx]*s[ky][kx]*q[ky][kx],
   fz[ky][kx] = zero,
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   g[ky][kx] = (affp/(kx**2+ky**2))*s[ky][kx],
   s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
   fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
   fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
   q[k][j] = complex charge density for fourier mode (j,k)
   fxy[k][j][0] = x component of complex force/charge,
   fxy[k][j][1] = y component of complex force/charge,
   fxy[k][j][2] = zero,
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
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
local data                                                 */
#define N 4
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
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,kk,kj,dky,at1,at2,at3,zt1,zt2,wp) \
reduction(+:sum1)
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
         fxy[N*(j+kj)] = at2*zt1;
         fxy[1+N*(j+kj)] = at3*zt1;
         fxy[2+N*(j+kj)] = zero;
         fxy[N*(j+k1)] = at2*zt2;
         fxy[1+N*(j+k1)] = -at3*zt2;
         fxy[2+N*(j+k1)] = zero;
         at1 = at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1]));
         wp += (double) at1;
      }
/* mode numbers kx = 0, nx/2 */
      at1 = crealf(ffc[kk])*cimagf(ffc[kk]);
      at3 = at1*dny*(float) k;
      zt1 = cimagf(q[kj]) - crealf(q[kj])*_Complex_I;
      fxy[N*kj] = zero;
      fxy[1+N*kj] = at3*zt1;
      fxy[2+N*kj] = zero;
      fxy[N*k1] = zero;
      fxy[1+N*k1] = zero;
      fxy[2+N*k1] = zero;
      at1 = at1*(q[kj]*conjf(q[kj]));
      wp += (double) at1;
      sum1 += wp;
   }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = N*nxvh*nyh;
#pragma ivdep
   for (j = 1; j < nxh; j++) {
      at1 = crealf(ffc[j])*cimagf(ffc[j]);
      at2 = at1*dnx*(float) j;  
      zt1 = cimagf(q[j]) - crealf(q[j])*_Complex_I;
      fxy[N*j] = at2*zt1;
      fxy[1+N*j] = zero;
      fxy[2+N*j] = zero;
      fxy[N*j+k1] = zero;
      fxy[1+N*j+k1] = zero;
      fxy[2+N*j+k1] = zero;
      at1 = at1*(q[j]*conjf(q[j]));
      wp += (double) at1;
   }
   fxy[0] = zero;
   fxy[1] = zero;
   fxy[2] = zero;
   fxy[k1] = zero;
   fxy[1+k1] = zero;
   fxy[2+k1] = zero;
   sum1 += wp;
   *we = sum1*(float) (nx*ny);
   return;
#undef N
}

/*--------------------------------------------------------------------*/
void cmcuperp2(float complex cu[], int nx, int ny, int nxvh, int nyv) {
/* this subroutine calculates the transverse current in fourier space
   input: all, output: cu
   approximate flop count is: 36*nxc*nyc
   and nxc*nyc divides
   where nxc = nx/2 - 1, nyc = ny/2 - 1
   the transverse current is calculated using the equation:
   cux[ky][kx] = cux[ky][kx]
                 -kx*(kx*cux[ky][kx]+ky*cuy[ky][kx])/(kx*kx+ky*ky)
   cuy[ky][kx] = cuy[ky][kx]
                 -ky*(kx*cux[ky][kx]+ky*cuy[ky][kx])/(kx*kx+ky*ky)
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
   and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
   cu[k][j][i] = complex current density for fourier mode (j,k)
   nx/ny = system length in x/y direction
   nxvh = second dimension of current array, must be >= nxh
   nyv = third dimension of current array, must be >= ny
local data                                                 */
#define N 4
   int nxh, nyh, j, k, k1, kj;
   float dnx, dny, dkx, dky, dky2, at1;
   float complex zero, zt1;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
/* calculate transverse part of current */
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for private(j,k,k1,kj,dky,dky2,dkx,at1,zt1)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      dky2 = dky*dky;
      kj = N*nxvh*k;
      k1 = N*nxvh*ny - kj;
#pragma ivdep
      for (j = 1; j < nxh; j++) {
         dkx = dnx*(float) j;
         at1 = 1./(dkx*dkx + dky2);
         zt1 = at1*(dkx*cu[N*j+kj] + dky*cu[1+N*j+kj]);
         cu[N*j+kj] -= dkx*zt1;
         cu[1+N*j+kj] -= dky*zt1;
         zt1 = at1*(dkx*cu[N*j+k1] - dky*cu[1+N*j+k1]);
         cu[N*j+k1] -= dkx*zt1;
         cu[1+N*j+k1] += dky*zt1;
      }
/* mode numbers kx = 0, nx/2 */
      cu[1+kj] = zero;
      cu[k1] = zero;
      cu[1+k1] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = N*nxvh*nyh;
#pragma ivdep
   for (j = 1; j < nxh; j++) {
      cu[N*j] = zero;
      cu[N*j+k1] = zero;
      cu[1+N*j+k1] = zero;
   }
   cu[0] = zero;
   cu[1] = zero;
   cu[k1] = zero;
   cu[1+k1] = zero;
   return;
#undef N
}

/*--------------------------------------------------------------------*/
void cvmibpois23(float complex cu[], float complex bxy[],
                 float complex ffc[], float ci, float *wm, int nx,
                 int ny, int nxvh, int nyv, int nxhd, int nyhd) {
/* this subroutine solves 2-1/2d poisson's equation in fourier space for
   magnetic field, with periodic boundary conditions.
   input: cu,ffc,ci,nx,ny,nxv,nyhd, output: bxy,wm
   approximate flop count is: 90*nxc*nyc + 40*(nxc + nyc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1
   the magnetic field is calculated using the equations:
   bx[ky][kx] = ci*ci*sqrt(-1)*g[ky][kx]*ky*cuz[ky][kx],
   by[ky][kx] = -ci*ci*sqrt(-1)*g[ky][kx]*kx*cuz[ky][kx],
   bz[ky][kx] = ci*ci*sqrt(-1)*g[ky][kx]*(kx*cuy[ky][kx]-ky*cux[ky][kx]),
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   g[ky][kx] = (affp/(kx**2+ky**2))*s[ky][kx],
   s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
   bx(kx=pi) = by(kx=pi) = bz(kx=pi) = bx(ky=pi) = by(ky=pi) = bz(ky=pi) 
   = 0, and bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
   cu[k][j][i] = complex current density for fourier mode (j,k)
   bxy[k][j][i] = i component of complex magnetic field
   all for fourier mode (j,k)
   cimag(ffc[k][j]) = finite-size particle shape factor s
   for fourier mode (j,k)
   creal(ffc[k][j]) = potential green's function g
   for fourier mode (j,k)
   ci = reciprocal of velocity of light
   magnetic field energy is also calculated, using
   wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
      |cu[ky][kx]*s[ky][kx]|**2), where
   affp = normalization constant = nx*ny/np, where np=number of particles
   this expression is valid only if the current is divergence-free
   nx/ny = system length in x/y direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
local data                                                 */
#define N 4
   int nxh, nyh, j, k, k1, kk, kj;
   float dnx, dny, dky, ci2, at1, at2, at3;
   float complex zero, zt1, zt2, zt3;
   double wp, sum1;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
   ci2 = ci*ci;
/* calculate magnetic field and sum field energy */
   sum1 = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,kk,kj,dky,at1,at2,at3,zt1,zt2,zt3,wp) \
reduction(+:sum1)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = N*nxvh*k;
      k1 = N*nxvh*ny - kj;
      wp = 0.0;
#pragma ivdep
      for (j = 1; j < nxh; j++) {
         at1 = ci2*crealf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         at1 = at1*cimagf(ffc[j+kk]);
         zt1 = -cimagf(cu[2+N*j+kj])
               + crealf(cu[2+N*j+kj])*_Complex_I;
         zt2 = -cimagf(cu[1+N*j+kj])
             + crealf(cu[1+N*j+kj])*_Complex_I;
         zt3 = -cimagf(cu[N*j+kj]) + crealf(cu[N*j+kj])*_Complex_I;
         bxy[N*j+kj] = at3*zt1;
         bxy[1+N*j+kj] = -at2*zt1;
         bxy[2+N*j+kj] = at2*zt2 - at3*zt3;
         zt1 = -cimagf(cu[2+N*j+k1])
               + crealf(cu[2+N*j+k1])*_Complex_I;
         zt2 = -cimagf(cu[1+N*j+k1])
             + crealf(cu[1+N*j+k1])*_Complex_I;
         zt3 = -cimagf(cu[N*j+k1]) + crealf(cu[N*j+k1])*_Complex_I;
         bxy[N*j+k1] = -at3*zt1;
         bxy[1+N*j+k1] = -at2*zt1;
         bxy[2+N*j+k1] = at2*zt2 + at3*zt3;
         at1 = at1*(cu[N*j+kj]*conjf(cu[N*j+kj])
               + cu[1+N*j+kj]*conjf(cu[1+N*j+kj])
               + cu[2+N*j+kj]*conjf(cu[2+N*j+kj])
               + cu[N*j+k1]*conjf(cu[N*j+k1])
               + cu[1+N*j+k1]*conjf(cu[1+N*j+k1])
               + cu[2+N*j+k1]*conjf(cu[2+N*j+k1]));
         wp += (double) at1;
      }
/* mode numbers kx = 0, nx/2 */
      at1 = ci2*crealf(ffc[kk]);
      at3 = at1*dny*(float) k;
      at1 = at1*cimagf(ffc[kk]);
      zt1 = -cimagf(cu[2+kj]) + crealf(cu[2+kj])*_Complex_I;
      zt3 = -cimagf(cu[kj]) + crealf(cu[kj])*_Complex_I;
      bxy[kj] = at3*zt1;
      bxy[1+kj] = zero;
      bxy[2+kj] = -at3*zt3;
      bxy[k1] = zero;
      bxy[1+k1] = zero;
      bxy[2+k1] = zero;
      at1 = at1*(cu[kj]*conjf(cu[kj]) + cu[1+kj]*conjf(cu[1+kj])
            + cu[2+kj]*conjf(cu[2+kj]));
      wp += (double) at1;
      sum1 += wp;
   }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = N*nxvh*nyh;
#pragma ivdep
   for (j = 1; j < nxh; j++) {
      at1 = ci2*crealf(ffc[j]);
      at2 = at1*dnx*(float) j; 
      at1 = at1*cimagf(ffc[j]);
      zt1 = -cimagf(cu[2+N*j]) + crealf(cu[2+N*j])*_Complex_I;
      zt2 = -cimagf(cu[1+N*j]) + crealf(cu[1+N*j])*_Complex_I;
      bxy[N*j] = zero;
      bxy[1+N*j] = -at2*zt1;
      bxy[2+N*j] = at2*zt2;
      bxy[N*j+k1] = zero;
      bxy[1+N*j+k1] = zero;
      bxy[2+N*j+k1] = zero;
      at1 = at1*(cu[N*j]*conjf(cu[N*j]) + cu[1+N*j]*conjf(cu[1+N*j])
            + cu[2+N*j]*conjf(cu[2+N*j]));
      wp += (double) at1;
   }
   bxy[0] = zero;
   bxy[1] = zero;
   bxy[2] = zero;
   bxy[k1] = zero;
   bxy[1+k1] = zero;
   bxy[2+k1] = zero;
   sum1 += wp;
   *wm = sum1*(float) (nx*ny);
   return;
#undef N
}

/*--------------------------------------------------------------------*/
void cvmmaxwel2(float complex exy[], float complex bxy[],
                float complex cu[], float complex ffc[], float ci,
                float dt, float *wf, float *wm, int nx, int ny,
                int nxvh, int nyv, int nxhd, int nyhd) {
/* this subroutine solves 2-1/2d maxwell's equation in fourier space for
   transverse electric and magnetic fields with periodic boundary
   conditions
   input: all, output: wf, wm, exy, bxy
   approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1
   the magnetic field is first updated half a step using the equations:
   bx[ky][kx] = bx[ky][kx] - .5*dt*sqrt(-1)*ky*ez[ky][kx]
   by[ky][kx] = by[ky][kx] + .5*dt*sqrt(-1)*kx*ez[ky][kx]
   bz[ky][kx] = bz[ky][kx] - .5*dt*sqrt(-1)*(kx*ey[ky][kx]-ky*ex[ky][kx])
   the electric field is then updated a whole step using the equations:
   ex[ky][kx] = ex[ky][kx] + c2*dt*sqrt(-1)*ky*bz[ky][kx]
                         - affp*dt*cux[ky][kx]*s[ky][kx]
   ey[ky][kx] = ey[ky][kx] - c2*dt*sqrt(-1)*kx*bz[ky][kx]
                         - affp*dt*cuy[ky][kx]*s[ky][kx]
   ez[ky][kx] = ez[ky][kx] + c2*dt*sqrt(-1)*(kx*by[ky][kx]-ky*bx[ky][kx])
                         - affp*dt*cuz[ky][kx]*s[ky][kx]
   the magnetic field is finally updated the remaining half step with
   the new electric field and the previous magnetic field equations.
   where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
   and s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)
   j,k = fourier mode numbers, except for
   ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
   ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
   ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
   and similarly for bx, by, bz.
   cu[k][j][i] = complex current density
   exy[k][j][i] = complex transverse electric field
   bxy[k][j][i] = complex magnetic field
   for component i, all for fourier mode (j,k)
   creal(ffc[0][0]) = affp = normalization constant = nx*ny/np,
   where np=number of particles
   cimag(ffc[k][j]) = finite-size particle shape factor s.
   s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)/2)
   for fourier mode (j-1,k-1)
   ci = reciprocal of velocity of light
   dt = time interval between successive calculations
   transverse electric field energy is also calculated, using
   wf = nx*ny**sum((1/affp)*|exy[ky][kx]|**2)
   magnetic field energy is also calculated, using
   wm = nx*ny**sum((c2/affp)*|bxy[ky][kx]|**2)
   nx/ny = system length in x/y direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
local data                                                 */
#define N 4
   int nxh, nyh, j, k, k1, kk, kj;
   float dnx, dny, dth, c2, cdt, affp, anorm, dkx, dky, afdt, adt;
   float at1;
   float complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9;
   double wp, ws, sum1, sum2;
   if (ci <= 0.0)
      return;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dth = 0.5*dt;
   c2 = 1.0/(ci*ci);
   cdt = c2*dt;
   affp = creal(ffc[0]);
   adt = affp*dt;
   zero = 0.0 + 0.0*_Complex_I;
   anorm = 1.0/affp;
/* update electromagnetic field and sum field energies */
   sum1 = 0.0;
   sum2 = 0.0;
/* calculate the electromagnetic fields */
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,kk,kj,dky,dkx,afdt,at1,zt1,zt2,zt3,zt4,zt5,zt6,zt7,zt8, \
zt9,ws,wp) \
reduction(+:sum1,sum2)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = N*nxvh*k;
      k1 = N*nxvh*ny - kj;
      ws = 0.0;
      wp = 0.0;
#pragma ivdep
      for (j = 1; j < nxh; j++) {
         dkx = dnx*(float) j;
         afdt = adt*cimagf(ffc[j+kk]);
/* update magnetic field half time step, ky > 0 */
         zt1 = -cimagf(exy[2+N*j+kj])
               + crealf(exy[2+N*j+kj])*_Complex_I;
         zt2 = -cimagf(exy[1+N*j+kj])
             + crealf(exy[1+N*j+kj])*_Complex_I;
         zt3 = -cimagf(exy[N*j+kj]) + crealf(exy[N*j+kj])*_Complex_I;
         zt4 = bxy[N*j+kj] - dth*(dky*zt1);
         zt5 = bxy[1+N*j+kj] + dth*(dkx*zt1);
         zt6 = bxy[2+N*j+kj] - dth*(dkx*zt2 - dky*zt3);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exy[N*j+kj] + cdt*(dky*zt1) - afdt*cu[N*j+kj];
         zt8 = exy[1+N*j+kj] - cdt*(dkx*zt1) - afdt*cu[1+N*j+kj];
         zt9 = exy[2+N*j+kj] + cdt*(dkx*zt2 - dky*zt3)
               - afdt*cu[2+N*j+kj];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exy[N*j+kj] = zt7;
         exy[1+N*j+kj] = zt8;
         exy[2+N*j+kj] = zt9;
         at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         ws += (double) at1;
         zt4 -= dth*(dky*zt1);
         zt5 += dth*(dkx*zt1);
         zt6 -= dth*(dkx*zt2 - dky*zt3);
         bxy[N*j+kj] = zt4;
         bxy[1+N*j+kj] = zt5;
         bxy[2+N*j+kj] = zt6;
         at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
         wp += (double) at1;
/* update magnetic field half time step, ky < 0 */
         zt1 = -cimagf(exy[2+N*j+k1])
               + crealf(exy[2+N*j+k1])*_Complex_I;
         zt2 = -cimagf(exy[1+N*j+k1])
               + crealf(exy[1+N*j+k1])*_Complex_I;
         zt3 = -cimagf(exy[N*j+k1]) + crealf(exy[N*j+k1])*_Complex_I;
         zt4 = bxy[N*j+k1] + dth*(dky*zt1);
         zt5 = bxy[1+N*j+k1] + dth*(dkx*zt1);
         zt6 = bxy[2+N*j+k1] - dth*(dkx*zt2 + dky*zt3);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exy[N*j+k1] - cdt*(dky*zt1) - afdt*cu[N*j+k1];
         zt8 = exy[1+N*j+k1] - cdt*(dkx*zt1) - afdt*cu[1+N*j+k1];
         zt9 = exy[2+N*j+k1] + cdt*(dkx*zt2 + dky*zt3)
               - afdt*cu[2+N*j+k1];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exy[N*j+k1] = zt7;
         exy[1+N*j+k1] = zt8;
         exy[2+N*j+k1] = zt9;
         at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         ws += (double) at1;
         zt4 += dth*(dky*zt1);
         zt5 += dth*(dkx*zt1);
         zt6 -= dth*(dkx*zt2 + dky*zt3);
         bxy[N*j+k1] = zt4;
         bxy[1+N*j+k1] = zt5;
         bxy[2+N*j+k1] = zt6;
         at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
         wp += (double) at1;
      }
/* mode numbers kx = 0, nx/2 */
      afdt = adt*cimagf(ffc[kk]);
/* update magnetic field half time step */
      zt1 = -cimagf(exy[2+kj]) + crealf(exy[2+kj])*_Complex_I;
      zt3 = -cimagf(exy[kj]) + crealf(exy[kj])*_Complex_I;
      zt4 = bxy[kj] - dth*(dky*zt1);
      zt6 = bxy[2+kj] + dth*(dky*zt3);
/* update electric field whole time step */
      zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
      zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
      zt7 = exy[kj] + cdt*(dky*zt1) - afdt*cu[kj];
      zt9 = exy[2+kj] - cdt*(dky*zt3) - afdt*cu[2+kj];
/* update magnetic field half time step and store electric field */
      zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
      zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
      exy[kj] = zt7;
      exy[1+kj] = zero;
      exy[2+kj] = zt9;
      at1 = anorm*(zt7*conjf(zt7) + zt9*conjf(zt9));
      ws += (double) at1;
      zt4 -= dth*(dky*zt1);
      zt6 += dth*(dky*zt3);
      bxy[kj] = zt4;
      bxy[1+kj] = zero;
      bxy[2+kj] = zt6;
      at1 = anorm*(zt4*conjf(zt4) + zt6*conjf(zt6));
      wp += (double) at1;
      bxy[k1] = zero;
      bxy[1+k1] = zero;
      bxy[2+k1] = zero;
      exy[k1] = zero;
      exy[1+k1] = zero;
      exy[2+k1] = zero;
      sum1 += ws;
      sum2 += wp;
   }
   ws = 0.0;
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = N*nxvh*nyh;
#pragma ivdep
   for (j = 1; j < nxh; j++) {
      dkx = dnx*(float) j; 
      afdt = adt*cimagf(ffc[j]);
/* update magnetic field half time step */
      zt1 = -cimagf(exy[2+N*j]) + crealf(exy[2+N*j])*_Complex_I;
      zt2 = -cimagf(exy[1+N*j]) + crealf(exy[1+N*j])*_Complex_I;
      zt5 = bxy[1+N*j] + dth*(dkx*zt1);
      zt6 = bxy[2+N*j] - dth*(dkx*zt2);
/* update electric field whole time step */
      zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
      zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
      zt8 = exy[1+N*j] - cdt*(dkx*zt1) - afdt*cu[1+N*j];
      zt9 = exy[2+N*j] + cdt*(dkx*zt2) - afdt*cu[2+N*j];
/* update magnetic field half time step and store electric field */
      zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
      zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
      exy[N*j] = zero;
      exy[1+N*j] = zt8;
      exy[2+N*j] = zt9;
      at1 = anorm*(zt8*conjf(zt8) + zt9*conjf(zt9));
      ws += (double) at1;
      zt5 += dth*(dkx*zt1);
      zt6 -= dth*(dkx*zt2);
      bxy[N*j] = zero;
      bxy[1+N*j] = zt5;
      bxy[2+N*j] = zt6;
      at1 = anorm*(zt5*conjf(zt5) + zt6*conjf(zt6));
      wp += (double) at1;
      bxy[N*j+k1] = zero;
      bxy[1+N*j+k1] = zero;
      bxy[2+N*j+k1] = zero;
      exy[N*j+k1] = zero;
      exy[1+N*j+k1] = zero;
      exy[2+N*j+k1] = zero;
   }
   bxy[0] = zero;
   bxy[1] = zero;
   bxy[2] = zero;
   exy[0] = zero;
   exy[1] = zero;
   exy[2] = zero;
   bxy[k1] = zero;
   bxy[1+k1] = zero;
   bxy[2+k1] = zero;
   exy[k1] = zero;
   exy[1+k1] = zero;
   exy[2+k1] = zero;
   sum1 += ws;
   sum2 += wp;
   *wf = sum1*(float) (nx*ny);
   *wm = sum2*c2*(float) (nx*ny);
   return;
#undef N
}

/*--------------------------------------------------------------------*/
void cvmemfield2(float complex fxy[], float complex exy[],
                 float complex ffc[], int isign, int nx, int ny,
                 int nxvh, int nyv, int nxhd, int nyhd) {
/* this subroutine either adds complex vector fields if isign > 0
   or copies complex vector fields if isign < 0
   includes additional smoothing
local data                                                 */
#define N 4
   int j, k, nxh, nyh, k1, kk, kj;
   float at1;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
/* add the fields */
   if (isign > 0) {
#pragma omp parallel for private(j,k,k1,kk,kj,at1)
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = N*nxvh*k;
         k1 = N*nxvh*ny - kj;
#pragma ivdep
         for (j = 0; j < nxh; j++) {
            at1 = cimagf(ffc[j+kk]);
            fxy[N*j+kj] += exy[N*j+kj]*at1;
            fxy[1+N*j+kj] += exy[1+N*j+kj]*at1;
            fxy[2+N*j+kj] += exy[2+N*j+kj]*at1;
            fxy[N*j+k1] += exy[N*j+k1]*at1;
            fxy[1+N*j+k1] += exy[1+N*j+k1]*at1;
            fxy[2+N*j+k1] += exy[2+N*j+k1]*at1;
         }
      }
      k1 = N*nxvh*nyh;
#pragma ivdep
      for (j = 0; j < nxh; j++) {
         at1 = cimagf(ffc[j]);
         fxy[N*j] += exy[N*j]*at1;
         fxy[1+N*j] += exy[1+N*j]*at1;
         fxy[2+N*j] += exy[2+N*j]*at1;
         fxy[N*j+k1] += exy[N*j+k1]*at1;
         fxy[1+N*j+k1] += exy[1+N*j+k1]*at1;
         fxy[2+N*j+k1] += exy[2+N*j+k1]*at1;
      }
   }
/* copy the fields */
   else if (isign < 0) {
#pragma omp parallel for private(j,k,k1,kk,kj,at1)
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = N*nxvh*k;
         k1 = N*nxvh*ny - kj;
#pragma ivdep
         for (j = 0; j < nxh; j++) {
            at1 = cimagf(ffc[j+kk]);
            fxy[N*j+kj] = exy[N*j+kj]*at1;
            fxy[1+N*j+kj] = exy[1+N*j+kj]*at1;
            fxy[2+N*j+kj] = exy[2+N*j+kj]*at1;
            fxy[N*j+k1] = exy[N*j+k1]*at1;
            fxy[1+N*j+k1] = exy[1+N*j+k1]*at1;
            fxy[2+N*j+k1] = exy[2+N*j+k1]*at1;
         }
      }
      k1 = N*nxvh*nyh;
      for (j = 0; j < nxh; j++) {
         at1 = cimagf(ffc[j]);
         fxy[N*j] = exy[N*j]*at1;
         fxy[1+N*j] = exy[1+N*j]*at1;
         fxy[2+N*j] = exy[2+N*j]*at1;
         fxy[N*j+k1] = exy[N*j+k1]*at1;
         fxy[1+N*j+k1] = exy[1+N*j+k1]*at1;
         fxy[2+N*j+k1] = exy[2+N*j+k1]*at1;
      }
   }
   return;
#undef N
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
void cfft2rvm3x(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nyi,
                int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the x part of 3 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   y, using complex arithmetic, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, two inverse fourier transforms are performed
   f[m][n][0:2] = (1/nx*ny)*sum(f[k][j][0:2]*
         exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, two forward fourier transforms are performed
   f[k][j][0:2] = sum(f[m][n][0:2]*exp(sqrt(-1)*2pi*n*j/nx)*
         exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = second dimension of f >= nx/2
   nyd = third dimension of f >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j][0:2] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][1][0:2] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0][0:2]) = real part of mode nx/2,0 and
   imag(f[0][ny/2][0:2]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, jj, j1, k1, k2, ns, ns2, km, kmr, joff;
   int nrxb;
   float at1, at2, ani;
   float complex t1, t2, t3, t4;
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
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,joff,at1,at2,ani,t1,t2,t3,t4)
   for (i = nyi-1; i < nyt; i++) {
      joff = 4*nxhd*i;
/* swap complex components */
      for (j = 0; j < nxh; j++) {
         at1 = cimagf(f[2+4*j+joff]);
         at2 = crealf(f[2+4*j+joff]);
         f[2+4*j+joff] = crealf(f[1+4*j+joff])
                         + crealf(f[3+4*j+joff])*_Complex_I;
         f[1+4*j+joff] = cimagf(f[4*j+joff]) + at1*_Complex_I;
         f[4*j+joff] = crealf(f[4*j+joff]) + at2*_Complex_I;
       }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            t1 = f[4*j1+joff];
            t2 = f[1+4*j1+joff];
            t3 = f[2+4*j1+joff];
            f[4*j1+joff] = f[4*j+joff];
            f[1+4*j1+joff] = f[1+4*j+joff];
            f[2+4*j1+joff] = f[2+4*j+joff];
            f[4*j+joff] = t1;
            f[1+4*j+joff] = t2;
            f[2+4*j+joff] = t3;
         }
      }
/* then transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = 4*ns2*k;
            k2 = k1 + 4*ns;
            for (j = 0; j < ns; j++) {
               t1 = sct[kmr*j];
               t2 = t1*f[4*j+k2+joff];
               t3 = t1*f[1+4*j+k2+joff];
               t4 = t1*f[2+4*j+k2+joff];
               f[4*j+k2+joff] = f[4*j+k1+joff] - t2;
               f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3;
               f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4;
               f[4*j+k1+joff] += t2;
               f[1+4*j+k1+joff] += t3;
               f[2+4*j+k1+joff] += t4;
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxy/nx;
      ani = 0.5/(((float) nx)*((float) ny));
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         for (jj = 0; jj < 3; jj++) {
            t2 = conjf(f[jj+4*(nxh-j)+joff]);
            t1 = f[jj+4*j+joff] + t2;
            t2 = (f[jj+4*j+joff] - t2)*t3;
            f[jj+4*j+joff] = ani*(t1 + t2);
            f[jj+4*(nxh-j)+joff] = ani*conjf(t1 - t2);
         }
      }
      ani = 2.0*ani;
      for (jj = 0; jj < 3; jj++) {
         f[jj+4*nxhh+joff] = ani*conjf(f[jj+4*nxhh+joff]);
         f[jj+joff] = ani*((crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
L100: nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,joff,at1,at2,t1,t2,t3,t4)
   for (i = nyi-1; i < nyt; i++) {
      joff = 4*nxhd*i;
/* scramble coefficients */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         for (jj = 0; jj < 3; jj++) {
            t2 = conjf(f[jj+4*(nxh-j)+joff]);
            t1 = f[jj+4*j+joff] + t2;
            t2 = (f[jj+4*j+joff] - t2)*t3;
            f[jj+4*j+joff] = t1 + t2;
            f[jj+4*(nxh-j)+joff] = conjf(t1 - t2);
         }
      }
      for (jj = 0; jj < 3; jj++) {
         f[jj+4*nxhh+joff] = 2.0*conjf(f[jj+4*nxhh+joff]);
         f[jj+joff] = (crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I;
      }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            t1 = f[4*j1+joff];
            t2 = f[1+4*j1+joff];
            t3 = f[2+4*j1+joff];
            f[4*j1+joff] = f[4*j+joff];
            f[1+4*j1+joff] = f[1+4*j+joff];
            f[2+4*j1+joff] = f[2+4*j+joff];
            f[4*j+joff] = t1;
            f[1+4*j+joff] = t2;
            f[2+4*j+joff] = t3;
         }
      }
/* then transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = 4*ns2*k;
            k2 = k1 + 4*ns;
            for (j = 0; j < ns; j++) {
               t1 = conjf(sct[kmr*j]);
               t2 = t1*f[4*j+k2+joff];
               t3 = t1*f[1+4*j+k2+joff];
               t4 = t1*f[2+4*j+k2+joff];
               f[4*j+k2+joff] = f[4*j+k1+joff] - t2;
               f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3;
               f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4;
               f[4*j+k1+joff] += t2;
               f[1+4*j+k1+joff] += t3;
               f[2+4*j+k1+joff] += t4;
            }
         }
         ns = ns2;
      }
/* swap complex components */
      for (j = 0; j < nxh; j++) {
         f[3+4*j+joff] = cimagf(f[2+4*j+joff])
                         + cimagf(f[3+4*j+joff])*_Complex_I;
         at1 = crealf(f[2+4*j+joff]);
         f[2+4*j+joff] = cimagf(f[4*j+joff])
                         + cimagf(f[1+4*j+joff])*_Complex_I;
         at2 = crealf(f[1+4*j+joff]);
         f[1+4*j+joff] = at1 + 0.0*_Complex_I;
         f[4*j+joff] = crealf(f[4*j+joff]) + at2*_Complex_I;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rvm3y(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nxi,
                int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the y part of 3 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   x, using complex arithmetic, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, two inverse fourier transforms are performed
   f[m][n][0:2] = (1/nx*ny)*sum(f[k][j][0:2] *
         exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, two forward fourier transforms are performed
   f[k][j][0:2] = sum(f[m][n][0:2]*exp(sqrt(-1)*2pi*n*j/nx)*
         exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxi = initial x index used
   nxp = number of x indices used
   nxhd = second dimension of f >= nx/2
   nyd = third dimension of f >= ny
  nxhyd = maximum of (nx/2,ny)
  nxyhd = maximum of (nx,ny)/2
  fourier coefficients are stored as follows:
   f[k][j][0:2] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][1][0:2] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0][0:2]) = real part of mode nx/2,0 and
   imag(f[0][ny/2][0:2]) = real part of mode nx/2,ny/2
  written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, koff;
   int nryb;
   float complex t1, t2, t3, t4;
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
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,koff,t1,t2,t3,t4)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = 4*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = 4*nxhd*k1;
            t1 = f[4*i+k1];
            t2 = f[1+4*i+k1];
            t3 = f[2+4*i+k1];
            f[4*i+k1] = f[4*i+koff];
            f[1+4*i+k1] = f[1+4*i+koff];
            f[2+4*i+k1] = f[2+4*i+koff];
            f[4*i+koff] = t1;
            f[1+4*i+koff] = t2;
            f[2+4*i+koff] = t3;
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
               j1 = 4*nxhd*(j + k1);
               j2 = 4*nxhd*(j + k2);
               t1 = sct[kmr*j];
               t2 = t1*f[4*i+j2];
               t3 = t1*f[1+4*i+j2];
               t4 = t1*f[2+4*i+j2];
               f[4*i+j2] = f[4*i+j1] - t2;
               f[1+4*i+j2] = f[1+4*i+j1] - t3;
               f[2+4*i+j2] = f[2+4*i+j1] - t4;
               f[4*i+j1] += t2;
               f[1+4*i+j1] += t3;
               f[2+4*i+j1] += t4;
            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if (nxi==1) {
      for (k = 1; k < nyh; k++) {
         koff = 4*nxhd*k;
         k1 = 4*nxhd*ny - koff;
         for (jj = 0; jj < 3; jj++) {
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
         koff = 4*nxhd*k;
         k1 = 4*nxhd*ny - koff;
         for (jj = 0; jj < 3; jj++) {
            t1 = cimagf(f[jj+k1]) + crealf(f[jj+k1])*_Complex_I;
            f[jj+k1] = conjf(f[jj+koff] - t1);
            f[jj+koff] += t1;
         }
      }
   }
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,koff,t1,t2,t3,t4)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = 4*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = 4*nxhd*k1;
            t1 = f[4*i+k1];
            t2 = f[1+4*i+k1];
            t3 = f[2+4*i+k1];
            f[4*i+k1] = f[4*i+koff];
            f[1+4*i+k1] = f[1+4*i+koff];
            f[2+4*i+k1] = f[2+4*i+koff];
            f[4*i+koff] = t1;
            f[1+4*i+koff] = t2;
            f[2+4*i+koff] = t3;
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
               j1 = 4*nxhd*(j + k1);
               j2 = 4*nxhd*(j + k2);
               t1 = conjf(sct[kmr*j]);
               t2 = t1*f[4*i+j2];
               t3 = t1*f[1+4*i+j2];
               t4 = t1*f[2+4*i+j2];
               f[4*i+j2] = f[4*i+j1] - t2;
               f[1+4*i+j2] = f[1+4*i+j1] - t3;
               f[2+4*i+j2] = f[2+4*i+j1] - t4;
               f[4*i+j1] += t2;
               f[1+4*i+j1] += t3;
               f[2+4*i+j1] += t4;
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
void cwfft2rvm3(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd) {
/* wrapper function for 3 2d real to complex ffts */
/* local data */
   int nxh, ny;
   static int nxi = 1, nyi = 1;
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   ny = 1L<<indy;
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cfft2rvm3x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                 nxyhd);
/* perform y fft */
      cfft2rvm3y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfft2rvm3y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                 nxyhd);
/* perform x fft */
      cfft2rvm3x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                nxyhd);
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void cdistr2h_(float *part, float *vtx, float *vty, float *vtz,
               float *vdx, float *vdy, float *vdz, int *npx, int *npy,
               int *idimp, int *nop, int *nx, int *ny, int *ipbc) {
   cdistr2h(part,*vtx,*vty,*vtz,*vdx,*vdy,*vdz,*npx,*npy,*idimp,*nop,
            *nx,*ny,*ipbc);
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
void cgbppush23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                   float *qbm, float *dt, float *dtc, float *ek,
                   int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                   int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                   int *ipbc) {
   cgbppush23lt(ppart,fxy,bxy,kpic,*qbm,*dt,*dtc,ek,*idimp,*nppmx,*nx,
                *ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppushf23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                    int *ncl, int *ihole, float *qbm, float *dt,
                    float *dtc, float *ek, int *idimp, int *nppmx,
                    int *nx, int *ny, int *mx, int *my, int *nxv,
                    int *nyv, int *mx1, int *mxy1, int *ntmax, 
                    int *irc) {
   cgbppushf23lt(ppart,fxy,bxy,kpic,ncl,ihole,*qbm,*dt,*dtc,ek,*idimp,
                 *nppmx,*nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbppush23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                    float *qbm, float *dt, float *dtc, float *ci,
                    float *ek, int *idimp, int *nppmx, int *nx, int *ny,
                    int *mx, int *my, int *nxv, int *nyv, int *mx1,
                    int *mxy1, int *ipbc) {
   cgrbppush23lt(ppart,fxy,bxy,kpic,*qbm,*dt,*dtc,*ci,ek,*idimp,*nppmx,
                 *nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbppushf23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                     int *ncl, int *ihole, float *qbm, float *dt,
                     float *dtc, float *ci, float *ek, int *idimp,
                     int *nppmx, int *nx, int *ny, int *mx, int *my,
                     int *nxv, int *nyv, int *mx1, int *mxy1,
                     int *ntmax, int *irc) {
   cgrbppushf23lt(ppart,fxy,bxy,kpic,ncl,ihole,*qbm,*dt,*dtc,*ci,ek,
                  *idimp,*nppmx,*nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,
                  *ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgbppush23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                    float *qbm, float *dt, float *dtc, float *ek,
                    int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                    int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                    int *ipbc) {
   cvgbppush23lt(ppart,fxy,bxy,kpic,*qbm,*dt,*dtc,ek,*idimp,*nppmx,*nx,
                 *ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgbppushf23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                     int *ncl, int *ihole, float *qbm, float *dt,
                     float *dtc, float *ek, int *idimp, int *nppmx,
                     int *nx, int *ny, int *mx, int *my, int *nxv,
                     int *nyv, int *mx1, int *mxy1, int *ntmax, 
                     int *irc) {
   cvgbppushf23lt(ppart,fxy,bxy,kpic,ncl,ihole,*qbm,*dt,*dtc,ek,*idimp,
                  *nppmx,*nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,
                  irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrbppush23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                     float *qbm, float *dt, float *dtc, float *ci,
                     float *ek, int *idimp, int *nppmx, int *nx,
                     int *ny, int *mx, int *my, int *nxv, int *nyv,
                     int *mx1, int *mxy1, int *ipbc) {
   cvgrbppush23lt(ppart,fxy,bxy,kpic,*qbm,*dt,*dtc,*ci,ek,*idimp,*nppmx,
                  *nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrbppushf23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                      int *ncl, int *ihole, float *qbm, float *dt,
                      float *dtc, float *ci, float *ek, int *idimp,
                      int *nppmx, int *nx, int *ny, int *mx, int *my,
                      int *nxv, int *nyv, int *mx1, int *mxy1,
                      int *ntmax, int *irc) {
   cvgrbppushf23lt(ppart,fxy,bxy,kpic,ncl,ihole,*qbm,*dt,*dtc,*ci,ek,
                   *idimp,*nppmx,*nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,
                   *ntmax,irc);
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
void cgjppost2lt_(float *ppart, float *cu, int *kpic, float *qm,
                  float *dt, int *nppmx, int *idimp, int *nx, int *ny, 
                  int *mx, int *my, int *nxv, int *nyv, int *mx1,
                  int *mxy1, int *ipbc) {
   cgjppost2lt(ppart,cu,kpic,*qm,*dt,*nppmx,*idimp,*nx,*ny,*mx,*my,*nxv,
               *nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgjppostf2lt_(float *ppart, float *cu, int *kpic, int *ncl,
                   int *ihole, float *qm, float *dt, int *nppmx, 
                   int *idimp, int *nx, int *ny, int *mx, int *my,
                   int *nxv, int *nyv, int *mx1, int *mxy1, int *ntmax,
                   int *irc) {
   cgjppostf2lt(ppart,cu,kpic,ncl,ihole,*qm,*dt,*nppmx,*idimp,*nx,*ny,
                *mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjppost2lt_(float *ppart, float *cu, int *kpic, float *qm,
                   float *dt, float *ci, int *nppmx, int *idimp,
                   int *nx, int *ny, int *mx, int *my, int *nxv,
                   int *nyv, int *mx1, int *mxy1, int *ipbc) {
   cgrjppost2lt(ppart,cu,kpic,*qm,*dt,*ci,*nppmx,*idimp,*nx,*ny,*mx,*my,
                *nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjppostf2lt_(float *ppart, float *cu, int *kpic, int *ncl,
                    int *ihole, float *qm, float *dt, float *ci,
                    int *nppmx, int *idimp, int *nx, int *ny, int *mx,
                    int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                    int *ntmax, int *irc) {
   cgrjppostf2lt(ppart,cu,kpic,ncl,ihole,*qm,*dt,*ci,*nppmx,*idimp,*nx,
                 *ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgjppost2lt_(float *ppart, float *cu, int *kpic, float *qm,
                   float *dt, int *nppmx, int *idimp, int *nx, int *ny, 
                   int *mx, int *my, int *nxv, int *nyv, int *mx1,
                   int *mxy1, int *ipbc) {
   cvgjppost2lt(ppart,cu,kpic,*qm,*dt,*nppmx,*idimp,*nx,*ny,*mx,*my,
                *nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgjppostf2lt_(float *ppart, float *cu, int *kpic, int *ncl,
                    int *ihole, float *qm, float *dt, int *nppmx, 
                    int *idimp, int *nx, int *ny, int *mx, int *my,
                    int *nxv, int *nyv, int *mx1, int *mxy1, int *ntmax,
                    int *irc) {
   cvgjppostf2lt(ppart,cu,kpic,ncl,ihole,*qm,*dt,*nppmx,*idimp,*nx,*ny,
                 *mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrjppost2lt_(float *ppart, float *cu, int *kpic, float *qm,
                    float *dt, float *ci, int *nppmx, int *idimp,
                    int *nx, int *ny, int *mx, int *my, int *nxv,
                    int *nyv, int *mx1, int *mxy1, int *ipbc) {
   cvgrjppost2lt(ppart,cu,kpic,*qm,*dt,*ci,*nppmx,*idimp,*nx,*ny,*mx,
                 *my,*nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrjppostf2lt_(float *ppart, float *cu, int *kpic, int *ncl,
                     int *ihole, float *qm, float *dt, float *ci,
                     int *nppmx, int *idimp, int *nx, int *ny, int *mx,
                     int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                     int *ntmax, int *irc) {
   cvgrjppostf2lt(ppart,cu,kpic,ncl,ihole,*qm,*dt,*ci,*nppmx,*idimp,*nx,
                  *ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
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
void cbguard2l_(float *bxy, int *nx, int *ny, int *nxe, int *nye) {
   cbguard2l(bxy,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void cacguard2l_(float *cu, int *nx, int *ny, int *nxe, int *nye) {
   cacguard2l(cu,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void caguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye) {
   caguard2l(q,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void cvmpois23_(float complex *q, float complex *fxy, int *isign,
                float complex *ffc, float *ax, float *ay, float *affp,
                float *we, int *nx, int *ny, int *nxvh, int *nyv,
                int *nxhd, int *nyhd) {
   cvmpois23(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*nxvh,*nyv,
             *nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmcuperp2_(float complex *cu, int *nx, int *ny, int *nxvh,
                int *nyv) {
   cmcuperp2(cu,*nx,*ny,*nxvh,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cvmibpois23_(float complex *cu, float complex *bxy,
                  float complex *ffc, float *ci, float *wm, int *nx,
                  int *ny, int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   cvmibpois23(cu,bxy,ffc,*ci,wm,*nx,*ny,*nxvh,*nyv,*nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cvmmaxwel2_(float complex *exy, float complex *bxy,
                 float complex *cu, float complex *ffc, float *ci,
                 float *dt, float *wf, float *wm, int *nx, int *ny,
                 int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   cvmmaxwel2(exy,bxy,cu,ffc,*ci,*dt,wf,wm,*nx,*ny,*nxvh,*nyv,*nxhd,
              *nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cvmemfield2_(float complex *fxy, float complex *exy,
                  float complex *ffc, int *isign, int *nx, int *ny,
                  int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   cvmemfield2(fxy,exy,ffc,*isign,*nx,*ny,*nxvh,*nyv,*nxhd,*nyhd);
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
void cfft2rvm3x_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *nyi,
                 int *nyp, int *nxhd, int *nyd, int *nxhyd,
                 int *nxyhd) {
   cfft2rvm3x(f,*isign,mixup,sct,*indx,*indy,*nyi,*nyp,*nxhd,*nyd,
              *nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rvm3y_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *nxi,
                 int *nxp, int *nxhd, int *nyd, int *nxhyd,
                 int *nxyhd) {
   cfft2rvm3y(f,*isign,mixup,sct,*indx,*indy,*nxi,*nxp,*nxhd,*nyd,
              *nxhyd,*nxyhd);
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
void cwfft2rvm3_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *nxhd,
                 int *nyd, int *nxhyd, int *nxyhd) {
   cwfft2rvm3(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,*nxyhd);
   return;
}