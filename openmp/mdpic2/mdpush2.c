/* C Library for Skeleton 2-1/2D Darwin OpenMP PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "mdpush2.h"

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
   input: all except kpic, nppmx, irc, output: kpic, nppmx, irc
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
void cppmovin2l(float part[], float ppart[], int kpic[], int nppmx,
                int idimp, int nop, int mx, int my, int mx1, int mxy1,
                int *irc) {
/* this subroutine sorts particles by x,y grid in tiles of mx, my
   and copies to segmented array ppart
   linear interpolation
   input: all except ppart, kpic, output: ppart, kpic
   part/ppart = input/output particle arrays
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
   ppart[m][n][2] = velocity vx of particle n in tile m
   ppart[m][n][3] = velocity vy of particle n in tile m
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
            ppart[i+idimp*(ip+nppmx*m)] = part[i+idimp*j];
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
void cppcheck2l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int ny, int mx, int my, int mx1, int my1, 
                int *irc) {
/* this subroutine performs a sanity check to make sure particles sorted
   by x,y grid in tiles of mx, my, are all within bounds.
   tiles are assumed to be arranged in 2D linear memory
   input: all except irc
   output: irc
   ppart[k][n][0] = position x of particle n in tile k
   ppart[k][n][1] = position y of particle n in tile k
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
         dx = ppart[idimp*(j+nppmx*k)];
         dy = ppart[1+idimp*(j+nppmx*k)];
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
void cgbppush23l(float ppart[], float fxy[], float bxy[], int kpic[],
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
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
   ppart[m][n][2] = x velocity of particle n in tile m
   ppart[m][n][3] = y velocity of particle n in tile m
   ppart[m][n][4] = z velocity of particle n in tile m
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
   int noff, moff, npoff, npp, mxv3;
   int i, j, k, nn, mm, nm;
   float qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt;
   float anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y;
   float sfxy[3*MXV*MYV], sbxy[3*MXV*MYV];
/* float sfxy[3*(mx+1)*(my+1)], sbxy[3*(mx+1)*(my+1)]; */
   double sum1, sum2;
   mxv3 = 3*(mx + 1);
   qtmh = 0.5*qbm*dt;
   sum2 = 0.0;
/* set boundary values */
   edgelx = 0.0;
   edgely = 0.0;
   edgerx = (float) nx;
   edgery = (float) ny;
   if (ipbc==2) {
      edgelx = 1.0;
      edgely = 1.0;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0;
      edgerx = (float) (nx-1);
   }
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,nm,x,y,dxp,dyp,amx,amy,dx,dy, \
dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4, \
rot5,rot6,rot7,rot8,rot9,sum1,sfxy,sbxy) \
reduction(+:sum2)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = nppmx*k;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[3*i+mxv3*j] = fxy[3*(i+noff+nxv*(j+moff))];
            sfxy[1+3*i+mxv3*j] = fxy[1+3*(i+noff+nxv*(j+moff))];
            sfxy[2+3*i+mxv3*j] = fxy[2+3*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[3*i+mxv3*j] = bxy[3*(i+noff+nxv*(j+moff))];
            sbxy[1+3*i+mxv3*j] = bxy[1+3*(i+noff+nxv*(j+moff))];
            sbxy[2+3*i+mxv3*j] = bxy[2+3*(i+noff+nxv*(j+moff))];
         }
      }
      sum1 = 0.0;
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 3*(nn - noff) + mxv3*(mm - moff);
         amx = 1.0 - dxp;
         amy = 1.0 - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + 3;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += mxv3;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + 3;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + 3;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += mxv3;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + 3;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[2+idimp*(j+npoff)] + dx;
         acy = ppart[3+idimp*(j+npoff)] + dy;
         acz = ppart[4+idimp*(j+npoff)] + dz;
/* time-centered kinetic energy */
         sum1 += (acx*acx + acy*acy + acz*acz);
/* calculate cyclotron frequency */
         omxt = qtmh*ox;
         omyt = qtmh*oy;
         omzt = qtmh*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0/(1.0 + omt);
         omt = 0.5*(1.0 - omt);
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
         dx += (rot1*acx + rot2*acy + rot3*acz)*anorm;
         dy += (rot4*acx + rot5*acy + rot6*acz)*anorm;
         dz += (rot7*acx + rot8*acy + rot9*acz)*anorm;
         ppart[2+idimp*(j+npoff)] = dx;
         ppart[3+idimp*(j+npoff)] = dy;
         ppart[4+idimp*(j+npoff)] = dz;
/* new position */
         dx = x + dx*dtc;
         dy = y + dy*dtc;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = ppart[idimp*(j+npoff)];
               ppart[2+idimp*(j+npoff)] = -ppart[2+idimp*(j+npoff)];
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = ppart[1+idimp*(j+npoff)];
               ppart[3+idimp*(j+npoff)] = -ppart[3+idimp*(j+npoff)];
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = ppart[idimp*(j+npoff)];
               ppart[2+idimp*(j+npoff)] = -ppart[2+idimp*(j+npoff)];
            }
         }
/* set new position */
         ppart[idimp*(j+npoff)] = dx;
         ppart[1+idimp*(j+npoff)] = dy;
      }
      sum2 += sum1;
   }
/* normalize kinetic energy */
   *ek += 0.5*sum2;
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgbppushf23l(float ppart[], float fxy[], float bxy[], int kpic[],
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
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
   ppart[m][n][2] = x velocity of particle n in tile m
   ppart[m][n][3] = y velocity of particle n in tile m
   ppart[m][n][4] = z velocity of particle n in tile m
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
   int noff, moff, npoff, npp, mxv3;
   int i, j, k, ih, nh, nn, mm, nm;
   float qtmh, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float x, y;
   float sfxy[3*MXV*MYV], sbxy[3*MXV*MYV];
/* float sfxy[3*(mx+1)*(my+1)], sbxy[3*(mx+1)*(my+1)]; */
   double sum1, sum2;
   mxv3 = 3*(mx + 1);
   qtmh = 0.5*qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   sum2 = 0.0;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,nm,ih,nh,x,y,dxp,dyp,amx,amy, \
dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3, \
rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,edgery,sum1,sfxy, \
sbxy) \
reduction(+:sum2)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = nppmx*k;
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
            sfxy[3*i+mxv3*j] = fxy[3*(i+noff+nxv*(j+moff))];
            sfxy[1+3*i+mxv3*j] = fxy[1+3*(i+noff+nxv*(j+moff))];
            sfxy[2+3*i+mxv3*j] = fxy[2+3*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[3*i+mxv3*j] = bxy[3*(i+noff+nxv*(j+moff))];
            sbxy[1+3*i+mxv3*j] = bxy[1+3*(i+noff+nxv*(j+moff))];
            sbxy[2+3*i+mxv3*j] = bxy[2+3*(i+noff+nxv*(j+moff))];
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
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 3*(nn - noff) + mxv3*(mm - moff);
         amx = 1.0 - dxp;
         amy = 1.0 - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + 3;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += mxv3;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + 3;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + 3;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += mxv3;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + 3;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[2+idimp*(j+npoff)] + dx;
         acy = ppart[3+idimp*(j+npoff)] + dy;
         acz = ppart[4+idimp*(j+npoff)] + dz;
/* time-centered kinetic energy */
         sum1 += (acx*acx + acy*acy + acz*acz);
/* calculate cyclotron frequency */
         omxt = qtmh*ox;
         omyt = qtmh*oy;
         omzt = qtmh*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0/(1.0 + omt);
         omt = 0.5*(1.0 - omt);
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
         dx += (rot1*acx + rot2*acy + rot3*acz)*anorm;
         dy += (rot4*acx + rot5*acy + rot6*acz)*anorm;
         dz += (rot7*acx + rot8*acy + rot9*acz)*anorm;
         ppart[2+idimp*(j+npoff)] = dx;
         ppart[3+idimp*(j+npoff)] = dy;
         ppart[4+idimp*(j+npoff)] = dz;
/* new position */
         dx = x + dx*dtc;
         dy = y + dy*dtc;
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
         ppart[idimp*(j+npoff)] = dx;
         ppart[1+idimp*(j+npoff)] = dy;
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
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgppost2l(float ppart[], float q[], int kpic[], float qm,
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
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
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
      npoff = nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < mxv*(my+1); j++) {
         sq[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
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
void cgjppost2l(float ppart[], float cu[], int kpic[], float qm,
                float dt, int nppmx, int idimp, int nx, int ny, int mx,
                int my, int nxv, int nyv, int mx1, int mxy1, int ipbc) {
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
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
   ppart[m][n][2] = x velocity of particle n in tile m
   ppart[m][n][3] = y velocity of particle n in tile m
   ppart[m][n][4] = z velocity of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of current array, must be >= nx+1
   nyv = second dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv3;
   int i, j, k, nn, mm;
   float edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz;
   float scu[3*MXV*MYV];
/* float scu[3*(mx+1)*(my+1)]; */
   mxv3 = 3*(mx + 1);
/* set boundary values */
   edgelx = 0.0;
   edgely = 0.0;
   edgerx = (float) nx;
   edgery = (float) ny;
   if (ipbc==2) {
      edgelx = 1.0;
      edgely = 1.0;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0;
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
      npoff = nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < mxv3*(my+1); j++) {
         scu[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         nn = 3*(nn - noff) + mxv3*(mm - moff);
         amx = qm - dxp;
         amy = 1.0 - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ppart[2+idimp*(j+npoff)];
         vy = ppart[3+idimp*(j+npoff)];
         vz = ppart[4+idimp*(j+npoff)];
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + 3;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += mxv3;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + 3;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = ppart[idimp*(j+npoff)];
               ppart[2+idimp*(j+npoff)] = -ppart[2+idimp*(j+npoff)];
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = ppart[1+idimp*(j+npoff)];
               ppart[3+idimp*(j+npoff)] = -ppart[3+idimp*(j+npoff)];
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = ppart[idimp*(j+npoff)];
               ppart[2+idimp*(j+npoff)] = -ppart[2+idimp*(j+npoff)];
            }
         }
/* set new position */
         ppart[idimp*(j+npoff)] = dx;
         ppart[1+idimp*(j+npoff)] = dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[3*(i+noff+nxv*(j+moff))] += scu[3*i+mxv3*j];
            cu[1+3*(i+noff+nxv*(j+moff))] += scu[1+3*i+mxv3*j];
            cu[2+3*(i+noff+nxv*(j+moff))] += scu[2+3*i+mxv3*j];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[3*(i+noff+nxv*moff)] += scu[3*i];
#pragma omp atomic
         cu[1+3*(i+noff+nxv*moff)] += scu[1+3*i];
#pragma omp atomic
         cu[2+3*(i+noff+nxv*moff)] += scu[2+3*i];
         if (mm > my) {
#pragma omp atomic
            cu[3*(i+noff+nxv*(mm+moff-1))] += scu[3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[1+3*(i+noff+nxv*(mm+moff-1))] += scu[1+3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[2+3*(i+noff+nxv*(mm+moff-1))] += scu[2+3*i+mxv3*(mm-1)];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[3*(noff+nxv*(j+moff))] += scu[mxv3*j];
#pragma omp atomic
         cu[1+3*(noff+nxv*(j+moff))] += scu[1+mxv3*j];
#pragma omp atomic
         cu[2+3*(noff+nxv*(j+moff))] += scu[2+mxv3*j];
         if (nn > mx) {
#pragma omp atomic
            cu[3*(nn+noff-1+nxv*(j+moff))] += scu[3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[1+3*(nn+noff-1+nxv*(j+moff))] += scu[1+3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[2+3*(nn+noff-1+nxv*(j+moff))] += scu[2+3*(nn-1)+mxv3*j];
         }
      }
   }
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgmjppost2l(float ppart[], float amu[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int my, int nxv, int nyv,
                 int mx1, int mxy1) {
/* for 2-1/2d code, this subroutine calculates particle momentum flux
   using first-order spline interpolation
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   51 flops/particle, 21 loads, 16 stores
   input: all, output: ppart, amu
   momentum flux is approximated by values at the nearest grid points
   amu(i,n,m)=qci*(1.-dx)*(1.-dy)
   amu(i,n+1,m)=qci*dx*(1.-dy)
   amu(i,n,m+1)=qci*(1.-dx)*dy
   amu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
   where vj = vj(t-dt/2) and vk = vk(t-dt/2)
   ppart[m][n][0] = position x of particle n in tile m at t
   ppart[m][n][1] = position y of particle n in tile m at t
   ppart[m][n][2] = x velocity of particle n in tile m at t - dt/2
   ppart[m][n][3] = y velocity of particle n in tile m at t - dt/2
   ppart[m][n][4] = z velocity of particle n in tile m at t - dt/2
   amu[k][j][i] = ith component of momentum flux at grid point j,k
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv4;
   int i, j, k, nn, mm;
   float dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz, v1, v2, v3, v4;
   float samu[4*MXV*MYV];
/* float samu[4*(mx+1)*(my+1)]; */
   mxv4 = 4*(mx + 1);
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,x,y,dxp,dyp,amx,amy,dx,dy,vx, \
vy,vz,v1,v2,v3,v4,samu)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < mxv4*(my+1); j++) {
         samu[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         nn = 4*(nn - noff) + mxv4*(mm - moff);
         amx = qm - dxp;
         amy = 1.0 - dyp;
/* deposit momentum flux */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ppart[2+idimp*(j+npoff)];
         vy = ppart[3+idimp*(j+npoff)];
         vz = ppart[4+idimp*(j+npoff)];
         v1 = vx*vx - vy*vy;
         v2 = vx*vy;
         v3 = vz*vx;
         v4 = vz*vy;
         samu[nn] += v1*dx;
         samu[nn+1] += v2*dx;
         samu[nn+2] += v3*dx;
         samu[nn+3] += v4*dx;
         dx = amx*dyp;
         mm = nn + 4;
         samu[mm] += v1*dy;
         samu[mm+1] += v2*dy;
         samu[mm+2] += v3*dy;
         samu[mm+3] += v4*dy;
         dy = dxp*dyp;
         nn += mxv4;
         samu[nn] += v1*dx;
         samu[nn+1] += v2*dx;
         samu[nn+2] += v3*dx;
         samu[nn+3] += v4*dx;
         mm = nn + 4;
         samu[mm] += v1*dy;
         samu[mm+1] += v2*dy;
         samu[mm+2] += v3*dy;
         samu[mm+3] += v4*dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            amu[4*(i+noff+nxv*(j+moff))] += samu[4*i+mxv4*j];
            amu[1+4*(i+noff+nxv*(j+moff))] += samu[1+4*i+mxv4*j];
            amu[2+4*(i+noff+nxv*(j+moff))] += samu[2+4*i+mxv4*j];
            amu[3+4*(i+noff+nxv*(j+moff))] += samu[3+4*i+mxv4*j];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         amu[4*(i+noff+nxv*moff)] += samu[4*i];
#pragma omp atomic
         amu[1+4*(i+noff+nxv*moff)] += samu[1+4*i];
#pragma omp atomic
         amu[2+4*(i+noff+nxv*moff)] += samu[2+4*i];
#pragma omp atomic
         amu[3+4*(i+noff+nxv*moff)] += samu[3+4*i];
         if (mm > my) {
#pragma omp atomic
            amu[4*(i+noff+nxv*(mm+moff-1))] += samu[4*i+mxv4*(mm-1)];
#pragma omp atomic
            amu[1+4*(i+noff+nxv*(mm+moff-1))] += samu[1+4*i+mxv4*(mm-1)];
#pragma omp atomic
            amu[2+4*(i+noff+nxv*(mm+moff-1))] += samu[2+4*i+mxv4*(mm-1)];
#pragma omp atomic
            amu[3+4*(i+noff+nxv*(mm+moff-1))] += samu[3+4*i+mxv4*(mm-1)];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         amu[4*(noff+nxv*(j+moff))] += samu[mxv4*j];
#pragma omp atomic
         amu[1+4*(noff+nxv*(j+moff))] += samu[1+mxv4*j];
#pragma omp atomic
         amu[2+4*(noff+nxv*(j+moff))] += samu[2+mxv4*j];
#pragma omp atomic
         amu[3+4*(noff+nxv*(j+moff))] += samu[3+mxv4*j];
         if (nn > mx) {
#pragma omp atomic
            amu[4*(nn+noff-1+nxv*(j+moff))] += samu[4*(nn-1)+mxv4*j];
#pragma omp atomic
            amu[1+4*(nn+noff-1+nxv*(j+moff))] += samu[1+4*(nn-1)+mxv4*j];
#pragma omp atomic
            amu[2+4*(nn+noff-1+nxv*(j+moff))] += samu[2+4*(nn-1)+mxv4*j];
#pragma omp atomic
            amu[3+4*(nn+noff-1+nxv*(j+moff))] += samu[3+4*(nn-1)+mxv4*j];
         }
      }
   }
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgdjppost2l(float ppart[], float fxy[], float bxy[], float dcu[],
                 float amu[], int kpic[], float qm, float qbm, float dt,
                 int idimp, int nppmx, int nx, int ny, int mx, int my,
                 int nxv, int nyv, int mx1, int mxy1) {
/* for 2-1/2d code, this subroutine calculates particle momentum flux
   and acceleration density using first-order spline interpolation.
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   194 flops/particle, 1 divide, 57 loads, 28 stores
   input: all, output: dcu, amu
   acceleration density is approximated by values at the nearest grid
   points
   dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
   dcu(i,n+1,m)=qci*dx*(1.-dy)
   dcu(i,n,m+1)=qci*(1.-dx)*dy
   dcu(i,n+1,m+1)=qci*dx*dy
   and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
   where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
   momentum flux is approximated by values at the nearest grid points
   amu(i,n,m)=qci*(1.-dx)*(1.-dy)
   amu(i,n+1,m)=qci*dx*(1.-dy)
   amu(i,n,m+1)=qci*(1.-dx)*dy
   amu(i,n+1,m+1)=qci*dx*dy
   and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
   where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
   and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
   where n,m = nearest grid points and dx = x-n, dy = y-m
   velocity equations at t=t+dt/2 are calculated from:
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
   fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   ppart[m][n][0] = position x of particle n in tile m at t
   ppart[m][n][1] = position y of particle n in tile m at t
   ppart[m][n][2] = x velocity of particle n in tile m at t - dt/2
   ppart[m][n][3] = y velocity of particle n in tile m at t - dt/2
   ppart[m][n][4] = z velocity of particle n in tile m at t - dt/2
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   dcu[k][j][i] = ith component of acceleration density
   at grid point j,k for i = 0, 2
   amu[k][j][i] = ith component of momentum flux
   at grid point j,k for i = 0, 3
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   idimp = size of phase space = 5
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv3, mxv4;
   int i, j, k, nn, mm, nm, mn;
   float qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, vx, vy, vz, v1, v2, v3, v4;
   float sfxy[3*MXV*MYV], sbxy[3*MXV*MYV];
   float sdcu[3*MXV*MYV], samu[4*MXV*MYV];
/* float sfxy[3*(mx+1)*(my+1)], sbxy[3*(mx+1)*(my+1)]; */
/* float sdcu[3*(mx+1)*(my+1)] samu[4*(mx+1)*(my+1)];  */
   mxv3 = 3*(mx + 1);
   mxv4 = 4*(mx + 1);
   qtmh = 0.5*qbm*dt;
   dti = 1.0/dt;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,nm,mn,x,y,vx,vy,vz,v1,v2,v3,v4, \
dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm, \
rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sfxy,sbxy,sdcu,samu)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = nppmx*k;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[3*i+mxv3*j] = fxy[3*(i+noff+nxv*(j+moff))];
            sfxy[1+3*i+mxv3*j] = fxy[1+3*(i+noff+nxv*(j+moff))];
            sfxy[2+3*i+mxv3*j] = fxy[2+3*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[3*i+mxv3*j] = bxy[3*(i+noff+nxv*(j+moff))];
            sbxy[1+3*i+mxv3*j] = bxy[1+3*(i+noff+nxv*(j+moff))];
            sbxy[2+3*i+mxv3*j] = bxy[2+3*(i+noff+nxv*(j+moff))];
         }
      }
/* zero out local accumulators */
      for (j = 0; j < mxv3*(my+1); j++) {
         sdcu[j] = 0.0f;
      }
      for (j = 0; j < mxv4*(my+1); j++) {
         samu[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 3*(nn - noff) + mxv3*(mm - moff);
         mn = 4*(nn - noff) + mxv4*(mm - moff);
         amx = 1.0 - dxp;
         amy = 1.0 - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + 3;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += mxv3;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + 3;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + 3;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += mxv3;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + 3;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         vx = ppart[2+idimp*(j+npoff)];
         vy = ppart[3+idimp*(j+npoff)];
         vz = ppart[4+idimp*(j+npoff)];
         acx = vx + dx;
         acy = vy + dy;
         acz = vz + dz;
/* calculate cyclotron frequency */
         omxt = qtmh*ox;
         omyt = qtmh*oy;
         omzt = qtmh*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0/(1.0 + omt);
         omt = 0.5*(1.0 - omt);
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
         dx += (rot1*acx + rot2*acy + rot3*acz)*anorm;
         dy += (rot4*acx + rot5*acy + rot6*acz)*anorm;
         dz += (rot7*acx + rot8*acy + rot9*acz)*anorm;
/* deposit momentum flux and acceleration density */
         amx = qm*amx;
         dxp = qm*dxp;
         ox = 0.5*(dx + vx);
         oy = 0.5*(dy + vy);
         oz = 0.5*(dz + vz);
         vx = dti*(dx - vx);
         vy = dti*(dy - vy);
         vz = dti*(dz - vz);
         dx = amx*amy;
         dy = dxp*amy;
         v1 = ox*ox - oy*oy;
         v2 = ox*oy;
         v3 = oz*ox;
         v4 = oz*oy;
         nn = mn;
         samu[nn] += v1*dx;
         samu[nn+1] += v2*dx;
         samu[nn+2] += v3*dx;
         samu[nn+3] += v4*dx;
         dx = amx*dyp;
         mm = nn + 4;
         samu[mm] += v1*dy;
         samu[mm+1] += v2*dy;
         samu[mm+2] += v3*dy;
         samu[mm+3] += v4*dy;
         dy = dxp*dyp;
         nn += mxv4;
         samu[nn] += v1*dx;
         samu[nn+1] += v2*dx;
         samu[nn+2] += v3*dx;
         samu[nn+3] += v4*dx;
         mm = nn + 4;
         samu[mm] += v1*dy;
         samu[mm+1] += v2*dy;
         samu[mm+2] += v3*dy;
         samu[mm+3] += v4*dy;
         dx = amx*amy;
         dy = dxp*amy;
         nn = nm;
         sdcu[nn] += vx*dx;
         sdcu[nn+1] += vy*dx;
         sdcu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + 3;
         sdcu[mm] += vx*dy;
         sdcu[mm+1] += vy*dy;
         sdcu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += mxv3;
         sdcu[nn] += vx*dx;
         sdcu[nn+1] += vy*dx;
         sdcu[nn+2] += vz*dx;
         mm = nn + 3;
         sdcu[mm] += vx*dy;
         sdcu[mm+1] += vy*dy;
         sdcu[mm+2] += vz*dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            amu[4*(i+noff+nxv*(j+moff))] += samu[4*i+mxv4*j];
            amu[1+4*(i+noff+nxv*(j+moff))] += samu[1+4*i+mxv4*j];
            amu[2+4*(i+noff+nxv*(j+moff))] += samu[2+4*i+mxv4*j];
            amu[3+4*(i+noff+nxv*(j+moff))] += samu[3+4*i+mxv4*j];
            dcu[3*(i+noff+nxv*(j+moff))] += sdcu[3*i+mxv3*j];
            dcu[1+3*(i+noff+nxv*(j+moff))] += sdcu[1+3*i+mxv3*j];
            dcu[2+3*(i+noff+nxv*(j+moff))] += sdcu[2+3*i+mxv3*j];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         amu[4*(i+noff+nxv*moff)] += samu[4*i];
#pragma omp atomic
         amu[1+4*(i+noff+nxv*moff)] += samu[1+4*i];
#pragma omp atomic
         amu[2+4*(i+noff+nxv*moff)] += samu[2+4*i];
#pragma omp atomic
         amu[3+4*(i+noff+nxv*moff)] += samu[3+4*i];
#pragma omp atomic
         dcu[3*(i+noff+nxv*moff)] += sdcu[3*i];
#pragma omp atomic
         dcu[1+3*(i+noff+nxv*moff)] += sdcu[1+3*i];
#pragma omp atomic
         dcu[2+3*(i+noff+nxv*moff)] += sdcu[2+3*i];
         if (mm > my) {
#pragma omp atomic
            amu[4*(i+noff+nxv*(mm+moff-1))] += samu[4*i+mxv4*(mm-1)];
#pragma omp atomic
            amu[1+4*(i+noff+nxv*(mm+moff-1))] += samu[1+4*i+mxv4*(mm-1)];
#pragma omp atomic
            amu[2+4*(i+noff+nxv*(mm+moff-1))] += samu[2+4*i+mxv4*(mm-1)];
#pragma omp atomic
            amu[3+4*(i+noff+nxv*(mm+moff-1))] += samu[3+4*i+mxv4*(mm-1)];
#pragma omp atomic
            dcu[3*(i+noff+nxv*(mm+moff-1))] += sdcu[3*i+mxv3*(mm-1)];
#pragma omp atomic
            dcu[1+3*(i+noff+nxv*(mm+moff-1))] += sdcu[1+3*i+mxv3*(mm-1)];
#pragma omp atomic
            dcu[2+3*(i+noff+nxv*(mm+moff-1))] += sdcu[2+3*i+mxv3*(mm-1)];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         amu[4*(noff+nxv*(j+moff))] += samu[mxv4*j];
#pragma omp atomic
         amu[1+4*(noff+nxv*(j+moff))] += samu[1+mxv4*j];
#pragma omp atomic
         amu[2+4*(noff+nxv*(j+moff))] += samu[2+mxv4*j];
#pragma omp atomic
         amu[3+4*(noff+nxv*(j+moff))] += samu[3+mxv4*j];
#pragma omp atomic
         dcu[3*(noff+nxv*(j+moff))] += sdcu[mxv3*j];
#pragma omp atomic
         dcu[1+3*(noff+nxv*(j+moff))] += sdcu[1+mxv3*j];
#pragma omp atomic
         dcu[2+3*(noff+nxv*(j+moff))] += sdcu[2+mxv3*j];
         if (nn > mx) {
#pragma omp atomic
            amu[4*(nn+noff-1+nxv*(j+moff))] += samu[4*(nn-1)+mxv4*j];
#pragma omp atomic
            amu[1+4*(nn+noff-1+nxv*(j+moff))] += samu[1+4*(nn-1)+mxv4*j];
#pragma omp atomic
            amu[2+4*(nn+noff-1+nxv*(j+moff))] += samu[2+4*(nn-1)+mxv4*j];
#pragma omp atomic
            amu[3+4*(nn+noff-1+nxv*(j+moff))] += samu[3+4*(nn-1)+mxv4*j];
#pragma omp atomic
            dcu[3*(nn+noff-1+nxv*(j+moff))] += sdcu[3*(nn-1)+mxv3*j];
#pragma omp atomic
            dcu[1+3*(nn+noff-1+nxv*(j+moff))] += sdcu[1+3*(nn-1)+mxv3*j];
#pragma omp atomic
            dcu[2+3*(nn+noff-1+nxv*(j+moff))] += sdcu[2+3*(nn-1)+mxv3*j];
         }
      }
   }
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cgdcjppost2l(float ppart[], float fxy[], float bxy[], float cu[],
                  float dcu[], float amu[], int kpic[], float qm, 
                  float qbm, float dt, int idimp, int nppmx, int nx,
                  int ny, int mx, int my, int nxv, int nyv, int mx1,
                  int mxy1) {
/* for 2-1/2d code, this subroutine calculates particle momentum flux,
   acceleration density and current density using first-order spline
   interpolation.
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   218 flops/particle, 1 divide, 69 loads, 40 stores
   input: all, output: cu, dcu, amu
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   and qci = qm*vj, where j = x,y,z, for i = 1, 3
   where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
   acceleration density is approximated by values at the nearest grid
   points
   dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
   dcu(i,n+1,m)=qci*dx*(1.-dy)
   dcu(i,n,m+1)=qci*(1.-dx)*dy
   dcu(i,n+1,m+1)=qci*dx*dy
   and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
   where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
   momentum flux is approximated by values at the nearest grid points
   amu(i,n,m)=qci*(1.-dx)*(1.-dy)
   amu(i,n+1,m)=qci*dx*(1.-dy)
   amu(i,n,m+1)=qci*(1.-dx)*dy
   amu(i,n+1,m+1)=qci*dx*dy
   and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
   where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
   and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
   where n,m = nearest grid points and dx = x-n, dy = y-m
   velocity equations at t=t+dt/2 are calculated from:
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
   fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   ppart[m][n][0] = position x of particle n in tile m at t
   ppart[m][n][1] = position y of particle n in tile m at t
   ppart[m][n][2] = x velocity of particle n in tile m at t - dt/2
   ppart[m][n][3] = y velocity of particle n in tile m at t - dt/2
   ppart[m][n][4] = z velocity of particle n in tile m at t - dt/2
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   cu[k][j][i] = ith component of current density at grid point j,k
   at grid point j,k for i = 0, 2
   dcu[k][j][i] = ith component of acceleration density
   at grid point j,k for i = 0, 2
   amu[k][j][i] = ith component of momentum flux
   at grid point j,k for i = 0, 3
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   idimp = size of phase space = 5
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv3, mxv4;
   int i, j, k, nn, mm, nm, mn;
   float qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, vx, vy, vz, v1, v2, v3, v4;
   float sfxy[3*MXV*MYV], sbxy[3*MXV*MYV];
   float scu[3*MXV*MYV], sdcu[3*MXV*MYV], samu[4*MXV*MYV];
/* float sfxy[3*(mx+1)*(my+1)], sbxy[3*(mx+1)*(my+1)]; */
/* float scu[3*(mx+1)*(my+1)], sdcu[3*(mx+1)*(my+1)];  */
/* float samu[4*(mx+1)*(my+1)];                        */
   mxv3 = 3*(mx + 1);
   mxv4 = 4*(mx + 1);
   qtmh = 0.5*qbm*dt;
   dti = 1.0/dt;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nn,mm,nm,mn,x,y,vx,vy,vz,v1,v2,v3,v4, \
dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm, \
rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sfxy,sbxy,scu,sdcu,samu)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = nppmx*k;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[3*i+mxv3*j] = fxy[3*(i+noff+nxv*(j+moff))];
            sfxy[1+3*i+mxv3*j] = fxy[1+3*(i+noff+nxv*(j+moff))];
            sfxy[2+3*i+mxv3*j] = fxy[2+3*(i+noff+nxv*(j+moff))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[3*i+mxv3*j] = bxy[3*(i+noff+nxv*(j+moff))];
            sbxy[1+3*i+mxv3*j] = bxy[1+3*(i+noff+nxv*(j+moff))];
            sbxy[2+3*i+mxv3*j] = bxy[2+3*(i+noff+nxv*(j+moff))];
         }
      }
/* zero out local accumulators */
      for (j = 0; j < mxv3*(my+1); j++) {
         scu[j] = 0.0f;
         sdcu[j] = 0.0f;
      }
      for (j = 0; j < mxv4*(my+1); j++) {
         samu[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 3*(nn - noff) + mxv3*(mm - moff);
         mn = 4*(nn - noff) + mxv4*(mm - moff);
         amx = 1.0 - dxp;
         amy = 1.0 - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + 3;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += mxv3;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + 3;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + 3;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += mxv3;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + 3;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         vx = ppart[2+idimp*(j+npoff)];
         vy = ppart[3+idimp*(j+npoff)];
         vz = ppart[4+idimp*(j+npoff)];
         acx = vx + dx;
         acy = vy + dy;
         acz = vz + dz;
/* calculate cyclotron frequency */
         omxt = qtmh*ox;
         omyt = qtmh*oy;
         omzt = qtmh*oz;
/* calculate rotation matrix */
         omt = omxt*omxt + omyt*omyt + omzt*omzt;
         anorm = 2.0/(1.0 + omt);
         omt = 0.5*(1.0 - omt);
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
         dx += (rot1*acx + rot2*acy + rot3*acz)*anorm;
         dy += (rot4*acx + rot5*acy + rot6*acz)*anorm;
         dz += (rot7*acx + rot8*acy + rot9*acz)*anorm;
/* deposit momentum flux, acceleration density, and current density */
         amx = qm*amx;
         dxp = qm*dxp;
         ox = 0.5*(dx + vx);
         oy = 0.5*(dy + vy);
         oz = 0.5*(dz + vz);
         vx = dti*(dx - vx);
         vy = dti*(dy - vy);
         vz = dti*(dz - vz);
         dx = amx*amy;
         dy = dxp*amy;
         v1 = ox*ox - oy*oy;
         v2 = ox*oy;
         v3 = oz*ox;
         v4 = oz*oy;
         nn = mn;
         samu[nn] += v1*dx;
         samu[nn+1] += v2*dx;
         samu[nn+2] += v3*dx;
         samu[nn+3] += v4*dx;
         dx = amx*dyp;
         mm = nn + 4;
         samu[mm] += v1*dy;
         samu[mm+1] += v2*dy;
         samu[mm+2] += v3*dy;
         samu[mm+3] += v4*dy;
         dy = dxp*dyp;
         nn += mxv4;
         samu[nn] += v1*dx;
         samu[nn+1] += v2*dx;
         samu[nn+2] += v3*dx;
         samu[nn+3] += v4*dx;
         mm = nn + 4;
         samu[mm] += v1*dy;
         samu[mm+1] += v2*dy;
         samu[mm+2] += v3*dy;
         samu[mm+3] += v4*dy;
         dx = amx*amy;
         dy = dxp*amy;
         nn = nm;
         sdcu[nn] += vx*dx;
         sdcu[nn+1] += vy*dx;
         sdcu[nn+2] += vz*dx;
         scu[nn] += ox*dx;
         scu[nn+1] += oy*dx;
         scu[nn+2] += oz*dx;
         dx = amx*dyp;
         mm = nn + 3;
         sdcu[mm] += vx*dy;
         sdcu[mm+1] += vy*dy;
         sdcu[mm+2] += vz*dy;
         scu[mm] += ox*dy;
         scu[mm+1] += oy*dy;
         scu[mm+2] += oz*dy;
         dy = dxp*dyp;
         nn += mxv3;
         sdcu[nn] += vx*dx;
         sdcu[nn+1] += vy*dx;
         sdcu[nn+2] += vz*dx;
         scu[nn] += ox*dx;
         scu[nn+1] += oy*dx;
         scu[nn+2] += oz*dx;
         mm = nn + 3;
         sdcu[mm] += vx*dy;
         sdcu[mm+1] += vy*dy;
         sdcu[mm+2] += vz*dy;
         scu[mm] += ox*dy;
         scu[mm+1] += oy*dy;
         scu[mm+2] += oz*dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            amu[4*(i+noff+nxv*(j+moff))] += samu[4*i+mxv4*j];
            amu[1+4*(i+noff+nxv*(j+moff))] += samu[1+4*i+mxv4*j];
            amu[2+4*(i+noff+nxv*(j+moff))] += samu[2+4*i+mxv4*j];
            amu[3+4*(i+noff+nxv*(j+moff))] += samu[3+4*i+mxv4*j];
            dcu[3*(i+noff+nxv*(j+moff))] += sdcu[3*i+mxv3*j];
            dcu[1+3*(i+noff+nxv*(j+moff))] += sdcu[1+3*i+mxv3*j];
            dcu[2+3*(i+noff+nxv*(j+moff))] += sdcu[2+3*i+mxv3*j];
            cu[3*(i+noff+nxv*(j+moff))] += scu[3*i+mxv3*j];
            cu[1+3*(i+noff+nxv*(j+moff))] += scu[1+3*i+mxv3*j];
            cu[2+3*(i+noff+nxv*(j+moff))] += scu[2+3*i+mxv3*j];
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         amu[4*(i+noff+nxv*moff)] += samu[4*i];
#pragma omp atomic
         amu[1+4*(i+noff+nxv*moff)] += samu[1+4*i];
#pragma omp atomic
         amu[2+4*(i+noff+nxv*moff)] += samu[2+4*i];
#pragma omp atomic
         amu[3+4*(i+noff+nxv*moff)] += samu[3+4*i];
#pragma omp atomic
         dcu[3*(i+noff+nxv*moff)] += sdcu[3*i];
#pragma omp atomic
         dcu[1+3*(i+noff+nxv*moff)] += sdcu[1+3*i];
#pragma omp atomic
         dcu[2+3*(i+noff+nxv*moff)] += sdcu[2+3*i];
#pragma omp atomic
         cu[3*(i+noff+nxv*moff)] += scu[3*i];
#pragma omp atomic
         cu[1+3*(i+noff+nxv*moff)] += scu[1+3*i];
#pragma omp atomic
         cu[2+3*(i+noff+nxv*moff)] += scu[2+3*i];
         if (mm > my) {
#pragma omp atomic
            amu[4*(i+noff+nxv*(mm+moff-1))] += samu[4*i+mxv4*(mm-1)];
#pragma omp atomic
            amu[1+4*(i+noff+nxv*(mm+moff-1))] += samu[1+4*i+mxv4*(mm-1)];
#pragma omp atomic
            amu[2+4*(i+noff+nxv*(mm+moff-1))] += samu[2+4*i+mxv4*(mm-1)];
#pragma omp atomic
            amu[3+4*(i+noff+nxv*(mm+moff-1))] += samu[3+4*i+mxv4*(mm-1)];
#pragma omp atomic
            dcu[3*(i+noff+nxv*(mm+moff-1))] += sdcu[3*i+mxv3*(mm-1)];
#pragma omp atomic
            dcu[1+3*(i+noff+nxv*(mm+moff-1))] += sdcu[1+3*i+mxv3*(mm-1)];
#pragma omp atomic
            dcu[2+3*(i+noff+nxv*(mm+moff-1))] += sdcu[2+3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[3*(i+noff+nxv*(mm+moff-1))] += scu[3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[1+3*(i+noff+nxv*(mm+moff-1))] += scu[1+3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[2+3*(i+noff+nxv*(mm+moff-1))] += scu[2+3*i+mxv3*(mm-1)];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         amu[4*(noff+nxv*(j+moff))] += samu[mxv4*j];
#pragma omp atomic
         amu[1+4*(noff+nxv*(j+moff))] += samu[1+mxv4*j];
#pragma omp atomic
         amu[2+4*(noff+nxv*(j+moff))] += samu[2+mxv4*j];
#pragma omp atomic
         amu[3+4*(noff+nxv*(j+moff))] += samu[3+mxv4*j];
#pragma omp atomic
         dcu[3*(noff+nxv*(j+moff))] += sdcu[mxv3*j];
#pragma omp atomic
         dcu[1+3*(noff+nxv*(j+moff))] += sdcu[1+mxv3*j];
#pragma omp atomic
         dcu[2+3*(noff+nxv*(j+moff))] += sdcu[2+mxv3*j];
#pragma omp atomic
         cu[3*(noff+nxv*(j+moff))] += scu[mxv3*j];
#pragma omp atomic
         cu[1+3*(noff+nxv*(j+moff))] += scu[1+mxv3*j];
#pragma omp atomic
         cu[2+3*(noff+nxv*(j+moff))] += scu[2+mxv3*j];
         if (nn > mx) {
#pragma omp atomic
            amu[4*(nn+noff-1+nxv*(j+moff))] += samu[4*(nn-1)+mxv4*j];
#pragma omp atomic
            amu[1+4*(nn+noff-1+nxv*(j+moff))] += samu[1+4*(nn-1)+mxv4*j];
#pragma omp atomic
            amu[2+4*(nn+noff-1+nxv*(j+moff))] += samu[2+4*(nn-1)+mxv4*j];
#pragma omp atomic
            amu[3+4*(nn+noff-1+nxv*(j+moff))] += samu[3+4*(nn-1)+mxv4*j];
#pragma omp atomic
            dcu[3*(nn+noff-1+nxv*(j+moff))] += sdcu[3*(nn-1)+mxv3*j];
#pragma omp atomic
            dcu[1+3*(nn+noff-1+nxv*(j+moff))] += sdcu[1+3*(nn-1)+mxv3*j];
#pragma omp atomic
            dcu[2+3*(nn+noff-1+nxv*(j+moff))] += sdcu[2+3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[3*(nn+noff-1+nxv*(j+moff))] += scu[3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[1+3*(nn+noff-1+nxv*(j+moff))] += scu[1+3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[2+3*(nn+noff-1+nxv*(j+moff))] += scu[2+3*(nn-1)+mxv3*j];
         }
      }
   }
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cpporder2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
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
   ppart[k][n][0] = position x of particle n in tile k
   ppart[k][n][1] = position y of particle n in tile k 
   ppbuff[k][n][i] = i co-ordinate of particle n in tile k
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
   int mxy1, noff, moff, npp, ncoff;
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
private(j,k,noff,moff,npp,nn,mm,ih,nh,ist,dx,dy,edgelx,edgely,edgerx, \
edgery)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
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
         dx = ppart[idimp*(j+nppmx*k)];
         dy = ppart[1+idimp*(j+nppmx*k)];
/* find particles going out of bounds */
         ist = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* ist = direction particle is going                             */
         if (dx >= edgerx) {
            if (dx >= anx)
               ppart[idimp*(j+nppmx*k)] = dx - anx;
            ist = 2;
         }
         else if (dx < edgelx) {
            if (dx < 0.0) {
               dx += anx;
               if (dx < anx)
                  ist = 1;
               else
                  dx = 0.0;
               ppart[idimp*(j+nppmx*k)] = dx;
            }
            else {
               ist = 1;
            }
         }
         if (dy >= edgery) {
            if (dy >= any)
               ppart[1+idimp*(j+nppmx*k)] = dy - any;
            ist += 6;
         }
         else if (dy < edgely) {
            if (dy < 0.0) {
               dy += any;
               if (dy < any)
                  ist += 3;
               else
                  dy = 0.0;
               ppart[1+idimp*(j+nppmx*k)] = dy;
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
private(i,j,k,isum,ist,nh,ip,j1,ii)
   for (k = 0; k < mxy1; k++) {
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
               ppbuff[i+idimp*(ii+npbmx*k)]
               = ppart[i+idimp*(j1+nppmx*k)];
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
private(i,j,k,ii,kk,npp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,ist,j1,j2,ip,ks)
   for (k = 0; k < mxy1; k++) {
      npp = kpic[k];
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
                  ppart[i+idimp*(j1+nppmx*k)]
                  = ppbuff[i+idimp*(j+ncoff+npbmx*ks[ii])];
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
      if (ih < nh) {
         ip = nh - ih;
         for (j = 0; j < ip; j++) {
            j1 = npp - j - 1;
            j2 = ihole[2*(nh-j+(ntmax+1)*k)] - 1;
            if (j1 > j2) {
/* move particle only if it is below current hole */
               for (i = 0; i < idimp; i++) {
                  ppart[i+idimp*(j2+nppmx*k)]
                  = ppart[i+idimp*(j1+nppmx*k)];
               }
            }
         }
         npp -= ip;
      }
      kpic[k] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
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
   cgppushf2l procedure.
   input: all except ppbuff, irc
   output: ppart, ppbuff, kpic, ncl, irc
   ppart[k][n][0] = position x of particle n in tile k
   ppart[k][n][1] = position y of particle n in tile k 
   ppbuff[k][n][i] = i co-ordinate of particle n in tile k
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
   int mxy1, npp, ncoff;
   int i, j, k, ii, kx, ky, ih, nh, ist, isum;
   int ip, j1, j2, kxl, kxr, kk, kl, kr;
   int ks[8];
   mxy1 = mx1*my1;
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,isum,ist,nh,ip,j1,ii)
   for (k = 0; k < mxy1; k++) {
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
               ppbuff[i+idimp*(ii+npbmx*k)]
               = ppart[i+idimp*(j1+nppmx*k)];
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
private(i,j,k,ii,kk,npp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,ist,j1,j2,ip,ks)
   for (k = 0; k < mxy1; k++) {
      npp = kpic[k];
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
                  ppart[i+idimp*(j1+nppmx*k)]
                  = ppbuff[i+idimp*(j+ncoff+npbmx*ks[ii])];
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
      if (ih < nh) {
         ip = nh - ih;
         for (j = 0; j < ip; j++) {
            j1 = npp - j - 1;
            j2 = ihole[2*(nh-j+(ntmax+1)*k)] - 1;
            if (j1 > j2) {
/* move particle only if it is below current hole */
               for (i = 0; i < idimp; i++) {
                  ppart[i+idimp*(j2+nppmx*k)]
                  = ppart[i+idimp*(j1+nppmx*k)];
               }
            }
         }
         npp -= ip;
      }
      kpic[k] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cbguard2l(float bxy[], int nx, int ny, int nxe, int nye) {
/* replicate extended periodic vector field bxy
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nxe = second dimension of field arrays, must be >= ny+1
local data                                                 */
   int j, k;
/* copy edges of extended field */
   for (k = 0; k < ny; k++) {
      bxy[3*nx+3*nxe*k] = bxy[3*nxe*k];
      bxy[1+3*nx+3*nxe*k] = bxy[1+3*nxe*k];
      bxy[2+3*nx+3*nxe*k] = bxy[2+3*nxe*k];
   }
   for (j = 0; j < nx; j++) {
      bxy[3*j+3*nxe*ny] = bxy[3*j];
      bxy[1+3*j+3*nxe*ny] = bxy[1+3*j];
      bxy[2+3*j+3*nxe*ny] = bxy[2+3*j];
   }
   bxy[3*nx+3*nxe*ny] = bxy[0];
   bxy[1+3*nx+3*nxe*ny] = bxy[1];
   bxy[2+3*nx+3*nxe*ny] = bxy[2];
   return;
}

/*--------------------------------------------------------------------*/
void cacguard2l(float cu[], int nx, int ny, int nxe, int nye) {
/* accumulate extended periodic vector field cu
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nxe = second dimension of field arrays, must be >= ny+1
local data                                                 */
   int j, k;
/* accumulate edges of extended field */
   for (k = 0; k < ny; k++) {
      cu[3*nxe*k] += cu[3*nx+3*nxe*k];
      cu[1+3*nxe*k] += cu[1+3*nx+3*nxe*k];
      cu[2+3*nxe*k] += cu[2+3*nx+3*nxe*k];
      cu[3*nx+3*nxe*k] = 0.0;
      cu[1+3*nx+3*nxe*k] = 0.0;
      cu[2+3*nx+3*nxe*k] = 0.0;
   }
   for (j = 0; j < nx; j++) {
      cu[3*j] += cu[3*j+3*nxe*ny];
      cu[1+3*j] += cu[1+3*j+3*nxe*ny];
      cu[2+3*j] += cu[2+3*j+3*nxe*ny];
      cu[3*j+3*nxe*ny] = 0.0;
      cu[1+3*j+3*nxe*ny] = 0.0;
      cu[2+3*j+3*nxe*ny] = 0.0;
   }
   cu[0] += cu[3*nx+3*nxe*ny];
   cu[1] += cu[1+3*nx+3*nxe*ny];
   cu[2] += cu[2+3*nx+3*nxe*ny];
   cu[3*nx+3*nxe*ny] = 0.0;
   cu[1+3*nx+3*nxe*ny] = 0.0;
   cu[2+3*nx+3*nxe*ny] = 0.0;
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
void camcguard2l(float amu[], int nx, int ny, int nxe, int nye,
                 int ndim) {
/* accumulate extended periodic tensor field
   linear interpolation
local data                                                 */
   int i, j, k, nnxe;
   nnxe = ndim*nxe;
/* accumulate edges of extended field */
   for (k = 0; k < ny; k++) {
      for (i = 0; i < ndim; i++) {
         amu[i+nnxe*k] += amu[i+ndim*nx+nnxe*k];
         amu[i+ndim*nx+nnxe*k] = 0.0;
      }
   }
   for (j = 0; j < nx; j++) {
      for (i = 0; i < ndim; i++) {
         amu[i+ndim*j] += amu[i+ndim*j+nnxe*ny];
         amu[i+ndim*j+nnxe*ny] = 0.0;
      }
   }
   for (i = 0; i < ndim; i++) {
      amu[i] += amu[i+ndim*nx+nnxe*ny];
      amu[i+ndim*nx+nnxe*ny] = 0.0;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cascfguard2l(float dcu[], float cus[], float q2m0, int nx, int ny,
                  int nxe, int nye) {
/* add scaled field to extended periodic field
local data                                                 */
   int i, j, k;
#pragma omp parallel for private(i,j,k)
   for (k = 0; k < ny; k++) {
      for (j = 0; j < nx; j++) {
         for (i = 0; i < 3; i++) {
            dcu[i+3*j+3*nxe*k] -= q2m0*cus[i+3*j+3*nxe*k];
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfwpminmx2(float qe[], float qbme, float *wpmax, float *wpmin,
                int nx, int ny, int nxe, int nye) {

/* calculates maximum and minimum plasma frequency.  assumes guard cells
   have already been added
   qe = charge density for electrons
   qbme = charge/mass ratio for electrons
   wpmax/wpmin = maximum/minimum plasma frequency
   nx/ny = system length in x/y direction
   nxe = first dimension of charge arrays, nxe must be >= nx
   nye = second dimension of charge arrays, nye must be >= ny
local data                                                 */
   int j, k;
   float tpmax, tpmin, at1;
   tpmax = qbme*qe[0];
   tpmin = tpmax;
#pragma omp parallel for private(j,k)
   for (k = 0; k < ny; k++) {
      for (j = 0; j < nx; j++) {
         at1 = qbme*qe[j+nxe*k];
#pragma omp critical
         tpmax = at1 > tpmax ? at1 : tpmax;
#pragma omp critical
         tpmin = at1 < tpmin ? at1 : tpmin;
      }
   }
   *wpmax = tpmax;
   *wpmin = tpmin;
   return;
}

/*--------------------------------------------------------------------*/
void cmpois23(float complex q[], float complex fxy[], int isign,
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
   nxvh = first dimension of field arrays, must be >= nxh
   nyv = second dimension of field arrays, must be >= ny
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
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
      for (j = 1; j < nxh; j++) {
         at1 = crealf(ffc[j+kk])*cimagf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         zt1 = cimagf(q[j+kj]) - crealf(q[j+kj])*_Complex_I;
         zt2 = cimagf(q[j+k1]) - crealf(q[j+k1])*_Complex_I;
         fxy[3*j+3*kj] = at2*zt1;
         fxy[1+3*j+3*kj] = at3*zt1;
         fxy[2+3*j+3*kj] = zero;
         fxy[3*j+3*k1] = at2*zt2;
         fxy[1+3*j+3*k1] = -at3*zt2;
         fxy[2+3*j+3*k1] = zero;
         wp += at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1]));
      }
/* mode numbers kx = 0, nx/2 */
      at1 = crealf(ffc[kk])*cimagf(ffc[kk]);
      at3 = at1*dny*(float) k;
      zt1 = cimagf(q[kj]) - crealf(q[kj])*_Complex_I;
      fxy[3*kj] = zero;
      fxy[1+3*kj] = at3*zt1;
      fxy[2+3*kj] = zero;
      fxy[3*k1] = zero;
      fxy[1+3*k1] = zero;
      fxy[2+3*k1] = zero;
      wp += at1*(q[kj]*conjf(q[kj]));
      sum1 += wp;
   }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = 3*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      at1 = crealf(ffc[j])*cimagf(ffc[j]);
      at2 = at1*dnx*(float) j;  
      zt1 = cimagf(q[j]) - crealf(q[j])*_Complex_I;
      fxy[3*j] = at2*zt1;
      fxy[1+3*j] = zero;
      fxy[2+3*j] = zero;
      fxy[3*j+k1] = zero;
      fxy[1+3*j+k1] = zero;
      fxy[2+3*j+k1] = zero;
      wp += at1*(q[j]*conjf(q[j]));
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
   nxvh = first dimension of current array, must be >= nxh
   nyv = second dimension of current array, must be >= ny
local data                                                 */
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
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      for (j = 1; j < nxh; j++) {
         dkx = dnx*(float) j;
         at1 = 1./(dkx*dkx + dky2);
         zt1 = at1*(dkx*cu[3*j+3*kj] + dky*cu[1+3*j+3*kj]);
         cu[3*j+3*kj] -= dkx*zt1;
         cu[1+3*j+3*kj] -= dky*zt1;
         zt1 = at1*(dkx*cu[3*j+3*k1] - dky*cu[1+3*j+3*k1]);
         cu[3*j+3*k1] -= dkx*zt1;
         cu[1+3*j+3*k1] += dky*zt1;
      }
/* mode numbers kx = 0, nx/2 */
      cu[1+3*kj] = zero;
      cu[3*k1] = zero;
      cu[1+3*k1] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = 3*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      cu[3*j] = zero;
      cu[3*j+k1] = zero;
      cu[1+3*j+k1] = zero;
   }
   cu[0] = zero;
   cu[1] = zero;
   cu[k1] = zero;
   cu[1+k1] = zero;
   return;
}

/*--------------------------------------------------------------------*/
void cmbbpois23(float complex cu[], float complex bxy[],
                float complex ffc[], float ci, float *wm, int nx,
                int ny, int nxvh, int nyv, int nxhd, int nyhd) {
/* this subroutine solves 2-1/2d poisson's equation in fourier space for
   magnetic field (or convolution of magnetic field over particle shape)
   with periodic boundary conditions.
   input: cu,ffc,ci,nx,ny,nxv,nyhd, output: bxy,wm
   approximate flop count is: 90*nxc*nyc + 40*(nxc + nyc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1
   the magnetic field is calculated using the equations:
   bx[ky][kx] = ci*ci*sqrt(-1)*g[ky][kx]*ky*cuz[ky][kx]*s[ky][kx],
   by[ky][kx] = -ci*ci*sqrt(-1)*g[ky][kx]*kx*cuz[ky][kx]*s[ky][kx],
   bz[ky][kx] = ci*ci*sqrt(-1)*g[ky][kx]*(kx*cuy[ky][kx]-ky*cux[ky][kx])*
               s[ky][kx],
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   g[ky][kx] = (affp/(kx**2+ky**2))*s[ky][kx],
   s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
   bx(kx=pi) = by(kx=pi) = bz(kx=pi) = bx(ky=pi) = by(ky=pi) = bz(ky=pi) 
   = 0, and bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
   cu[k][j][i] = complex current density for fourier mode (j,k)
   bxy[k][j][0] = x component of complex magnetic field
   bxy[k][j][1] = y component of complex magnetic field
   bxy[k][j][2] = z component of complex magnetic field
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
   nxvh = first dimension of field arrays, must be >= nxh
   nyv = second dimension of field arrays, must be >= ny
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
local data                                                 */
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
/* calculate smoothed magnetic field and sum field energy */
   sum1 = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,kk,kj,dky,at1,at2,at3,zt1,zt2,zt3,wp) \
reduction(+:sum1)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      wp = 0.0;
      for (j = 1; j < nxh; j++) {
         at1 = ci2*crealf(ffc[j+kk])*cimagf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         zt1 = -cimagf(cu[2+3*j+3*kj])
               + crealf(cu[2+3*j+3*kj])*_Complex_I;
         zt2 = -cimagf(cu[1+3*j+3*kj])
             + crealf(cu[1+3*j+3*kj])*_Complex_I;
         zt3 = -cimagf(cu[3*j+3*kj]) + crealf(cu[3*j+3*kj])*_Complex_I;
         bxy[3*j+3*kj] = at3*zt1;
         bxy[1+3*j+3*kj] = -at2*zt1;
         bxy[2+3*j+3*kj] = at2*zt2 - at3*zt3;
         zt1 = -cimagf(cu[2+3*j+3*k1])
               + crealf(cu[2+3*j+3*k1])*_Complex_I;
         zt2 = -cimagf(cu[1+3*j+3*k1])
             + crealf(cu[1+3*j+3*k1])*_Complex_I;
         zt3 = -cimagf(cu[3*j+3*k1]) + crealf(cu[3*j+3*k1])*_Complex_I;
         bxy[3*j+3*k1] = -at3*zt1;
         bxy[1+3*j+3*k1] = -at2*zt1;
         bxy[2+3*j+3*k1] = at2*zt2 + at3*zt3;
         wp += at1*(cu[3*j+3*kj]*conjf(cu[3*j+3*kj])
               + cu[1+3*j+3*kj]*conjf(cu[1+3*j+3*kj])
               + cu[2+3*j+3*kj]*conjf(cu[2+3*j+3*kj])
               + cu[3*j+3*k1]*conjf(cu[3*j+3*k1])
               + cu[1+3*j+3*k1]*conjf(cu[1+3*j+3*k1])
               + cu[2+3*j+3*k1]*conjf(cu[2+3*j+3*k1]));
      }
/* mode numbers kx = 0, nx/2 */
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      at1 = ci2*crealf(ffc[kk])*cimagf(ffc[kk]);
      at3 = at1*dny*(float) k;
      zt1 = -cimagf(cu[2+3*kj]) + crealf(cu[2+3*kj])*_Complex_I;
      zt3 = -cimagf(cu[3*kj]) + crealf(cu[3*kj])*_Complex_I;
      bxy[3*kj] = at3*zt1;
      bxy[1+3*kj] = zero;
      bxy[2+3*kj] = -at3*zt3;
      bxy[3*k1] = zero;
      bxy[1+3*k1] = zero;
      bxy[2+3*k1] = zero;
      wp += at1*(cu[3*kj]*conjf(cu[3*kj]) + cu[1+3*kj]*conjf(cu[1+3*kj])
            + cu[2+3*kj]*conjf(cu[2+3*kj]));
      sum1 += wp;
   }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = 3*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      at1 = ci2*crealf(ffc[j])*cimagf(ffc[j]);
      at2 = at1*dnx*(float) j; 
      zt1 = -cimagf(cu[2+3*j]) + crealf(cu[2+3*j])*_Complex_I;
      zt2 = -cimagf(cu[1+3*j]) + crealf(cu[1+3*j])*_Complex_I;
      bxy[3*j] = zero;
      bxy[1+3*j] = -at2*zt1;
      bxy[2+3*j] = at2*zt2;
      bxy[3*j+k1] = zero;
      bxy[1+3*j+k1] = zero;
      bxy[2+3*j+k1] = zero;
      wp += at1*(cu[3*j]*conjf(cu[3*j]) + cu[1+3*j]*conjf(cu[1+3*j])
            + cu[2+3*j]*conjf(cu[2+3*j]));
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
}

/*--------------------------------------------------------------------*/
void cbaddext2(float bxy[], float omx, float omy, float omz, int nx,
               int ny, int nxe, int nye) {
/* adds constant to magnetic field for 2-1/2d code
   bxy = magnetic field
   omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
   nx/ny = system length in x/y direction
   nxe = second dimension of magnetic field array, nxe must be >= nx
   nye = third dimension of magnetic field array, nye must be >= ny
local data                                                 */
   int j, k;
#pragma omp parallel for private(j,k)
   for (k = 0; k < ny; k++) {
      for (j = 0; j < nx; j++) {
         bxy[3*j+3*nxe*k] += omx;
         bxy[1+3*j+3*nxe*k] += omy;
         bxy[2+3*j+3*nxe*k] += omz;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cmdcuperp23(float complex dcu[], float complex amu[], int nx,
                 int ny, int nxvh, int nyv) {
/* this subroutine calculates transverse part of the derivative of
   the current density from the momentum flux
   in 2-1/2d with periodic boundary conditions.
   the derivative of the current is calculated using the equations:
   dcu[ky][kx][0] = -sqrt(-1)*(kx*vx*vx+ky*vx*vy)
   dcu[ky][kx][1] = -sqrt(-1)*(kx*vx*vy+ky*vy*vy)
   dcu[ky][kx][2] = -sqrt(-1)*(kx*vx*vz+ky*vy*vz)
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   except for dcu(i,kx=pi) = dcu(i,ky=pi) =  dcu(i,kx=0,ky=0) = 0.
   the transverse part is calculated using the equation:
   dcu[ky][kx][0] = dcu[ky][kx][0]-kx*(kx*dcu[ky][kx][0]
                  + ky*dcu[ky][kx][1])/(kx*kx+ky*ky)
   dcu[ky][kx][1] = dcu[ky][kx][1]-ky*(kx*dcu[ky][kx][1]
                  + ky*dcu[ky][kx][1])/(kx*kx+ky*ky)
   on output:
   dcu[k][j][i] = transverse part of complex derivative of current for
   fourier mode (j-1,k-1)
   amu[k][j][0] = xx-yy component of complex momentum flux
   amu[k][j][1] = xy component of complex momentum flux
   amu[k][j][2] = zx component of complex momentum flux
   amu[k][j][3] = zy component of complex momentum flux
   all for fourier mode (j-1,k-1)
   nx/ny = system length in x/y direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
local data                                                 */
   int nxh, nyh, j, k, k1, kj;
   float dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2, at1;
   float complex zero, zt1, zt2, zt3;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,kj,dky,dky2,dkx,dkx2,dkxy,dkxy2,at1,zt1,zt2,zt3)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      dky2 = dky*dky;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      for (j = 1; j < nxh; j++) {
         dkx = dnx*(float) j;
         dkx2 = dkx*dkx;
         dkxy = dkx*dky;
         dkxy2 = dky2 - dkx2;
         at1 = 1.0/(dkx2 + dky2);
         zt1 = cimagf(amu[4*j+4*kj]) - crealf(amu[4*j+4*kj])*_Complex_I;
         zt2 = cimagf(amu[1+4*j+4*kj])
             - crealf(amu[1+4*j+4*kj])*_Complex_I;
         zt3 = at1*(dkxy*zt1 + dkxy2*zt2);
         dcu[3*j+3*kj] = dky*zt3;
         dcu[1+3*j+3*kj] = -dkx*zt3;
         zt1 = cimagf(amu[2+4*j+4*kj])
             - crealf(amu[2+4*j+4*kj])*_Complex_I;
         zt2 = cimagf(amu[3+4*j+4*kj])
             - crealf(amu[3+4*j+4*kj])*_Complex_I;
         dcu[2+3*j+3*kj] = dkx*zt1 + dky*zt2;
         zt1 = cimagf(amu[4*j+4*k1]) - crealf(amu[4*j+4*k1])*_Complex_I;
         zt2 = cimagf(amu[1+4*j+4*k1])
             - crealf(amu[1+4*j+4*k1])*_Complex_I;
         zt3 = at1*(dkxy*zt1 - dkxy2*zt2);
         dcu[3*j+3*k1] = dky*zt3;
         dcu[1+3*j+3*k1] = dkx*zt3;
         zt1 = cimagf(amu[2+4*j+4*k1])
             - crealf(amu[2+4*j+4*k1])*_Complex_I;
         zt2 = cimagf(amu[3+4*j+4*k1])
             - crealf(amu[3+4*j+4*k1])*_Complex_I;
         dcu[2+3*j+3*k1] = dkx*zt1 - dky*zt2;
      }
/* mode numbers kx = 0, nx/2 */
      zt2 = cimagf(amu[1+4*kj]) - crealf(amu[1+4*kj])*_Complex_I;
      dcu[3*kj] = dky*zt2;
      dcu[1+3*kj] = zero;
      zt2 = cimagf(amu[3+4*kj]) - crealf(amu[3+4*kj])*_Complex_I;
      dcu[2+3*kj] = dky*zt2;
      dcu[3*k1] = zero;
      dcu[1+3*k1] = zero;
      dcu[2+3*k1] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = 3*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      dkx = dnx*(float) j;
      zt2 = cimagf(amu[1+4*j]) - crealf(amu[1+4*j])*_Complex_I;
      dcu[3*j] = zero;
      dcu[1+3*j] = dkx*zt2;
      zt1 = cimagf(amu[2+4*j]) - crealf(amu[2+4*j])*_Complex_I;
      dcu[2+3*j] = dkx*zt1;
      dcu[3*j+k1] = zero;
      dcu[1+3*j+k1] = zero;
      dcu[2+3*j+k1] = zero;
   }
   dcu[0] = zero;
   dcu[1] = zero;
   dcu[2] = zero;
   dcu[k1] = zero;
   dcu[1+k1] = zero;
   dcu[2+k1] = zero;
   return;
}

/*--------------------------------------------------------------------*/
void cmadcuperp23(float complex dcu[], float complex amu[], int nx,
                  int ny, int nxvh, int nyv) {
/* this subroutine calculates transverse part of the derivative of
   the current density from the momentum flux and acceleration density
   in 2-1/2d with periodic boundary conditions.
   the derivative of the current is calculated using the equations:
   dcu[ky][kx][0] = dcu[ky][kx][0]-sqrt(-1)*(kx*vx*vx+ky*vx*vy)
   dcu[ky][kx][1] = dcu[ky][kx][1]-sqrt(-1)*(kx*vx*vy+ky*vy*vy)
   dcu[ky][kx][2] = dcu[ky][kx][2]-sqrt(-1)*(kx*vx*vz+ky*vy*vz)
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   except for dcu(i,kx=pi) = dcu(i,ky=pi) =  dcu(i,kx=0,ky=0) = 0.
   the transverse part is calculated using the equation:
   dcu[ky][kx][0] = dcu[ky][kx][0]-kx*(kx*dcu[ky][kx][0]
                  + ky*dcu[ky][kx][1])/(kx*kx+ky*ky)
   dcu[ky][kx][1] = dcu[ky][kx][1]-ky*(kx*dcu[ky][kx][0]
                  + ky*dcu[ky][kx][1])/(kx*kx+ky*ky)
   on input:
   dcu[k][j][i] = complex acceleration density for fourier mode (j-1,k-1)
   on output:
   dcu[k][j][i] = transverse part of complex derivative of current for
   fourier mode (j-1,k-1)
   amu[k][j][0] = xx-yy component of complex momentum flux
   amu[k][j][1] = xy component of complex momentum flux
   amu[k][j][2] = zx component of complex momentum flux
   amu[k][j][3] = zy component of complex momentum flux
   all for fourier mode (j-1,k-1)
   nx/ny = system length in x/y direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
local data                                                 */
   int nxh, nyh, j, k, k1, kj;
   float dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2, at1;
   float complex zero, zt1, zt2, zt3;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,kj,dky,dky2,dkx,dkx2,dkxy,dkxy2,at1,zt1,zt2,zt3)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      dky2 = dky*dky;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      for (j = 1; j < nxh; j++) {
         dkx = dnx*(float) j;
         dkx2 = dkx*dkx;
         dkxy = dkx*dky;
         dkxy2 = dky2 - dkx2;
         at1 = 1.0/(dkx2 + dky2);
         zt1 = cimagf(amu[4*j+4*kj]) - crealf(amu[4*j+4*kj])*_Complex_I;
         zt2 = cimagf(amu[1+4*j+4*kj])
             - crealf(amu[1+4*j+4*kj])*_Complex_I;
         zt3 = at1*(dky*dcu[3*j+3*kj] - dkx*dcu[1+3*j+3*kj] + dkxy*zt1  
             + dkxy2*zt2);
         dcu[3*j+3*kj] = dky*zt3;
         dcu[1+3*j+3*kj] = -dkx*zt3;
         zt1 = cimagf(amu[2+4*j+4*kj])
             - crealf(amu[2+4*j+4*kj])*_Complex_I;
         zt2 = cimagf(amu[3+4*j+4*kj])
             - crealf(amu[3+4*j+4*kj])*_Complex_I;
         dcu[2+3*j+3*kj] += dkx*zt1 + dky*zt2;
         zt1 = cimagf(amu[4*j+4*k1]) - crealf(amu[4*j+4*k1])*_Complex_I;
         zt2 = cimagf(amu[1+4*j+4*k1])
             - crealf(amu[1+4*j+4*k1])*_Complex_I;
         zt3 = at1*(dky*dcu[3*j+3*k1] + dkx*dcu[1+3*j+3*k1] + dkxy*zt1  
             - dkxy2*zt2);
         dcu[3*j+3*k1] = dky*zt3;
         dcu[1+3*j+3*k1] = dkx*zt3;
         zt1 = cimagf(amu[2+4*j+4*k1])
             - crealf(amu[2+4*j+4*k1])*_Complex_I;
         zt2 = cimagf(amu[3+4*j+4*k1])
             - crealf(amu[3+4*j+4*k1])*_Complex_I;
         dcu[2+3*j+3*k1] += dkx*zt1 - dky*zt2;
      }
/* mode numbers kx = 0, nx/2 */
      zt2 = cimagf(amu[1+4*kj]) - crealf(amu[1+4*kj])*_Complex_I;
      dcu[3*kj] += dky*zt2;
      dcu[1+3*kj] = zero;
      zt2 = cimagf(amu[3+4*kj]) - crealf(amu[3+4*kj])*_Complex_I;
      dcu[2+3*kj] += dky*zt2;
      dcu[3*k1] = zero;
      dcu[1+3*k1] = zero;
      dcu[2+3*k1] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = 3*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      dkx = dnx*(float) j;
      zt2 = cimagf(amu[1+4*j]) - crealf(amu[1+4*j])*_Complex_I;
      dcu[3*j] = zero;
      dcu[1+3*j] += dkx*zt2;
      zt1 = cimagf(amu[2+4*j]) - crealf(amu[2+4*j])*_Complex_I;
      dcu[2+3*j] += dkx*zt1;
      dcu[3*j+k1] = zero;
      dcu[1+3*j+k1] = zero;
      dcu[2+3*j+k1] = zero;
   }
   dcu[0] = zero;
   dcu[1] = zero;
   dcu[k1] = zero;
   dcu[1+k1] = zero;
   dcu[2+k1] = zero;
   return;
}

/*--------------------------------------------------------------------*/
void cmepois23(float complex dcu[], float complex exy[], int isign,
               float complex ffe[], float ax, float ay, float affp,
               float wp0, float ci, float *wf, int nx, int ny, int nxvh,
               int nyv, int nxhd, int nyhd) {
/* this subroutine solves 2-1/2d poisson's equation in fourier space for
   transverse electric field (or convolution of transverse electric field
   over particle shape), with periodic boundary conditions.
   using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
   A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
   for isign = 0, input: isign,ax,ay,affp,wp0,nx,ny,nxvh,nyhd, output:ffe
   for isign /= 0, input: dcu,ffe,isign,ci,nx,ny,nxvh,nyv,nxhd,nyhd,
   output: exy,wf
   approximate flop count is: 68*nxc*nyc + 33*(nxc + nyc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1
   if isign = 0, form factor array is prepared
   if isign = -1, smoothed transverse electric field is calculated
   using the equation:
   ex[ky][kx] = -ci*ci*g[ky][kx]*dcux[ky][kx]*s[ky][kx]
   ey[ky][kx] = -ci*ci*g[ky][kx]*dcuy[ky][kx]*s[ky][kx]
   ez[ky][kx] = -ci*ci*g[ky][kx]*dcuz[ky][kx]*s[ky][kx]
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   g[ky][kx] = (affp/(kx**2+ky**2+wp0*ci2*s[ky][kx]**2))*s[ky][kx],
   s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
   ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = ex(ky=pi) = ey(ky=pi) = ez(ky=pi) 
   = 0, and ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
   if isign = 1, unsmoothed transverse electric field is calculated
   using the equation:
   ex[ky][kx] = -ci*ci*g[ky][kx]*dcux([ky][kx]
   ey[ky][kx] = -ci*ci*g[ky][kx]*dcuy[ky][kx]
   ez[ky][kx] = -ci*ci*g[ky][kx]*dcuz[ky][kx]
   dcu(i,j,k) = transverse part of complex derivative of current for
   fourier mode (j-1,k-1)
   exy[ky][kx][0] = x component of complex transverse electric field
   exy[ky][kx][1] = y component of complex transverse electric field
   exy[ky][kx][2] = z component of complex transverse electric field
   all for fourier mode (j-1,k-1)
   cimagf(ffe[ky][kx]) = finite-size particle shape factor s
   for fourier mode (j-1,k-1)
   crealf(ffe[ky][kx]) = potential green's function g
   for fourier mode (j-1,k-1)
   ax/ay = half-width of particle in x/y direction
   affp = normalization constant = nx*ny/np, where np=number of particles
   wp0 = normalized total plasma frequency squared
   ci = reciprocal of velocity of light
   transverse electric field energy is also calculated, using
   wf = nx*ny*sum((affp/((kx**2+ky**2)*ci*ci)**2)
      |dcu(kx,ky)*s[ky][kx]|**2)
   this expression is valid only if the derivative of current is
   divergence-free
   nx/ny = system length in x/y direction
   nxvh = first dimension of field arrays, must be >= nxh
   nyv = second dimension of field arrays, must be >= ny
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, j, k, k1, kk, kj;
   float dnx, dny, ci2, wpc, dkx, dky, at1, at2, at3, at4;
   float complex zero;
   double wp, sum1;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
   ci2 = ci*ci;
   if (isign != 0)
      goto L30;
   wpc = wp0*ci2;
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
            ffe[j+kk] = affp + 1.0*_Complex_I;
         }
         else {
            ffe[j+kk] = (affp*at4/(at3 + wpc*at4*at4)) + at4*_Complex_I;
         }
      }
   }
   return;
/* calculate smoothed transverse electric field and sum field energy */
L30: if (isign > 0)
        goto L70;
   sum1 = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for private(j,k,k1,kk,kj,at1,at2,wp) \
reduction(+:sum1)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      wp = 0.0;
      for (j = 1; j < nxh; j++) {
         at2 = -ci2*crealf(ffe[j+kk]);
         at1 = at2*cimagf(ffe[j+kk]);
         at2 = at2*at2;
         exy[3*j+3*kj] = at1*dcu[3*j+3*kj];
         exy[1+3*j+3*kj] = at1*dcu[1+3*j+3*kj];
         exy[2+3*j+3*kj] = at1*dcu[2+3*j+3*kj];
         exy[3*j+3*k1] = at1*dcu[3*j+3*k1];
         exy[1+3*j+3*k1] = at1*dcu[1+3*j+3*k1];
         exy[2+3*j+3*k1] = at1*dcu[2+3*j+3*k1];
         wp += at2*(dcu[3*j+3*kj]*conjf(dcu[3*j+3*kj])
               + dcu[1+3*j+3*kj]*conjf(dcu[1+3*j+3*kj])
               + dcu[2+3*j+3*kj]*conjf(dcu[2+3*j+3*kj])
               + dcu[3*j+3*k1]*conjf(dcu[3*j+3*k1])
               + dcu[1+3*j+3*k1]*conjf(dcu[1+3*j+3*k1])
               + dcu[2+3*j+3*k1]*conjf(dcu[2+3*j+3*k1]));
      }
/* mode numbers kx = 0, nx/2 */
      at2 = -ci2*crealf(ffe[kk]);
      at1 = at2*cimagf(ffe[kk]);
      at2 = at2*at2;
      exy[3*kj] = at1*dcu[3*kj];
      exy[1+3*kj] = at1*dcu[1+3*kj];
      exy[2+3*kj] = at1*dcu[2+3*kj];
      exy[3*k1] = zero;
      exy[1+3*k1] = zero;
      exy[2+3*k1] = zero;
      wp += at2*(dcu[3*kj]*conjf(dcu[3*kj])
         + dcu[1+3*kj]*conjf(dcu[1+3*kj])
         + dcu[2+3*kj]*conjf(dcu[2+3*kj]));
      sum1 += wp;
   }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = 3*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      at2 = -ci2*crealf(ffe[j]);
      at1 = at2*cimagf(ffe[j]);
      at2 = at2*at2;
      exy[3*j] = at1*dcu[3*j];
      exy[1+3*j] = at1*dcu[1+3*j];
      exy[2+3*j] = at1*dcu[2+3*j];
      exy[3*j+k1] = zero;
      exy[1+3*j+k1] = zero;
      exy[2+3*j+k1] = zero;
      wp += at2*(dcu[3*j]*conjf(dcu[3*j]) + dcu[1+3*j]*conjf(dcu[1+3*j])
         + dcu[2+3*j]*conjf(dcu[2+3*j]));
   }
   exy[0] = zero;
   exy[1] = zero;
   exy[2] = zero;
   exy[k1] = zero;
   exy[1+k1] = zero;
   exy[2+k1] = zero;
   sum1 += wp;
   *wf = (sum1*(float) (nx*ny))/crealf(ffe[0]);
   return;
/* calculate unsmoothed transverse electric field and sum field energy */
L70: sum1 = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for private(j,k,k1,kk,kj,at1,at2,wp) \
reduction(+:sum1)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      wp = 0.0;
      for (j = 1; j < nxh; j++) {
         at2 = -ci2*crealf(ffe[j+kk]);
         at1 = at2*at2;
         exy[3*j+3*kj] = at2*dcu[3*j+3*kj];
         exy[1+3*j+3*kj] = at2*dcu[1+3*j+3*kj];
         exy[2+3*j+3*kj] = at2*dcu[2+3*j+3*kj];
         exy[3*j+3*k1] = at2*dcu[3*j+3*k1];
         exy[1+3*j+3*k1] = at2*dcu[1+3*j+3*k1];
         exy[2+3*j+3*k1] = at2*dcu[2+3*j+3*k1];
         wp += at1*(dcu[3*j+3*kj]*conjf(dcu[3*j+3*kj])
               + dcu[1+3*j+3*kj]*conjf(dcu[1+3*j+3*kj])
               + dcu[2+3*j+3*kj]*conjf(dcu[2+3*j+3*kj])
               + dcu[3*j+3*k1]*conjf(dcu[3*j+3*k1])
               + dcu[1+3*j+3*k1]*conjf(dcu[1+3*j+3*k1])
               + dcu[2+3*j+3*k1]*conjf(dcu[2+3*j+3*k1]));
     }
/* mode numbers kx = 0, nx/2 */
      at2 = -ci2*crealf(ffe[kk]);
      at1 = at2*at2;
      exy[3*kj] = at2*dcu[3*kj];
      exy[1+3*kj] = at2*dcu[1+3*kj];
      exy[2+3*kj] = at2*dcu[2+3*kj];
      exy[3*k1] = zero;
      exy[1+3*k1] = zero;
      exy[2+3*k1] = zero;
      wp += at1*(dcu[3*kj]*conjf(dcu[3*kj])
         + dcu[1+3*kj]*conjf(dcu[1+3*kj])
         + dcu[2+3*kj]*conjf(dcu[2+3*kj]));
      sum1 += wp;
  }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = 3*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      at2 = -ci2*crealf(ffe[j]);
      at1 = at2*at2;
      exy[3*j] = at2*dcu[3*j];
      exy[1+3*j] = at2*dcu[1+3*j];
      exy[2+3*j] = at2*dcu[2+3*j];
      exy[3*j+k1] = zero;
      exy[1+3*j+k1] = zero;
      exy[2+3*j+k1] = zero;
      wp += at1*(dcu[3*j]*conjf(dcu[3*j]) + dcu[1+3*j]*conjf(dcu[1+3*j])
         + dcu[2+3*j]*conjf(dcu[2+3*j]));
   }
   exy[0] = zero;
   exy[1] = zero;
   exy[2] = zero;
   exy[k1] = zero;
   exy[1+k1] = zero;
   exy[2+k1] = zero;
   sum1 += wp;
   *wf = (sum1*(float) (nx*ny))/crealf(ffe[0]);
   return;
}

/*--------------------------------------------------------------------*/
void caddvrfield2(float a[], float b[], float c[], int ndim, int nxe,
                  int nye) {
/* this subroutine calculates a = b + c for real vector fields
local data                                                 */
   int i, j, k, nnxe;
   nnxe = ndim*nxe;
#pragma omp parallel for private(i,j,k)
   for (k = 0; k < nye; k++) {
      for (j = 0; j < nxe; j++) {
         for (i = 0; i < ndim; i++) {
            a[i+ndim*j+nnxe*k] = b[i+ndim*j+nnxe*k]
                               + c[i+ndim*j+nnxe*k];
         }
      }
   }
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
void cfft2rmxx(float complex f[], int isign, int mixup[],
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
   int nrx, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, nrxb, joff;
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
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,joff,ani,t1,t2,t3)
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
               j1 = j + k1;
               j2 = j + k2;
               t1 = sct[kmr*j];
               t2 = t1*f[j2+joff];
               f[j2+joff] = f[j1+joff] - t2;
               f[j1+joff] += t2;
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
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,joff,t1,t2,t3)
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
               j1 = j + k1;
               j2 = j + k2;
               t1 = conjf(sct[kmr*j]);
               t2 = t1*f[j2+joff];
               f[j2+joff] = f[j1+joff] - t2;
               f[j1+joff] += t2;
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
void cfft2rm3x(float complex f[], int isign, int mixup[],
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
   int nrx, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
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
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,joff,at1,at2,ani,t1,t2,t3 \
,t4)
   for (i = nyi-1; i < nyt; i++) {
      joff = 3*nxhd*i;
/* swap complex components */
      for (j = 0; j < nxh; j++) {
         at1 = crealf(f[2+3*j+joff]);
         f[2+3*j+joff] = crealf(f[1+3*j+joff])
                         + cimagf(f[2+3*j+joff])*_Complex_I;
         at2 = cimagf(f[1+3*j+joff]);
         f[1+3*j+joff] = cimagf(f[3*j+joff]) + at1*_Complex_I;
         f[3*j+joff] = crealf(f[3*j+joff]) + at2*_Complex_I;
       }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            t1 = f[3*j1+joff];
            t2 = f[1+3*j1+joff];
            t3 = f[2+3*j1+joff];
            f[3*j1+joff] = f[3*j+joff];
            f[1+3*j1+joff] = f[1+3*j+joff];
            f[2+3*j1+joff] = f[2+3*j+joff];
            f[3*j+joff] = t1;
            f[1+3*j+joff] = t2;
            f[2+3*j+joff] = t3;
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
               j1 = j + k1;
               j2 = j + k2;
               t1 = sct[kmr*j];
               t2 = t1*f[3*j2+joff];
               t3 = t1*f[1+3*j2+joff];
               t4 = t1*f[2+3*j2+joff];
               f[3*j2+joff] = f[3*j1+joff] - t2;
               f[1+3*j2+joff] = f[1+3*j1+joff] - t3;
               f[2+3*j2+joff] = f[2+3*j1+joff] - t4;
               f[3*j1+joff] += t2;
               f[1+3*j1+joff] += t3;
               f[2+3*j1+joff] += t4;
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
            t2 = conjf(f[jj+3*(nxh-j)+joff]);
            t1 = f[jj+3*j+joff] + t2;
            t2 = (f[jj+3*j+joff] - t2)*t3;
            f[jj+3*j+joff] = ani*(t1 + t2);
            f[jj+3*(nxh-j)+joff] = ani*conjf(t1 - t2);
         }
      }
      ani = 2.0*ani;
      for (jj = 0; jj < 3; jj++) {
         f[jj+3*nxhh+joff] = ani*conjf(f[jj+3*nxhh+joff]);
         f[jj+joff] = ani*((crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
L100: nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,joff,at1,at2,t1,t2,t3,t4)
   for (i = nyi-1; i < nyt; i++) {
      joff = 3*nxhd*i;
/* scramble coefficients */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         for (jj = 0; jj < 3; jj++) {
            t2 = conjf(f[jj+3*(nxh-j)+joff]);
            t1 = f[jj+3*j+joff] + t2;
            t2 = (f[jj+3*j+joff] - t2)*t3;
            f[jj+3*j+joff] = t1 + t2;
            f[jj+3*(nxh-j)+joff] = conjf(t1 - t2);
         }
      }
      for (jj = 0; jj < 3; jj++) {
         f[jj+3*nxhh+joff] = 2.0*conjf(f[jj+3*nxhh+joff]);
         f[jj+joff] = (crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I;
      }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            t1 = f[3*j1+joff];
            t2 = f[1+3*j1+joff];
            t3 = f[2+3*j1+joff];
            f[3*j1+joff] = f[3*j+joff];
            f[1+3*j1+joff] = f[1+3*j+joff];
            f[2+3*j1+joff] = f[2+3*j+joff];
            f[3*j+joff] = t1;
            f[1+3*j+joff] = t2;
            f[2+3*j+joff] = t3;
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
               j1 = j + k1;
               j2 = j + k2;
               t1 = conjf(sct[kmr*j]);
               t2 = t1*f[3*j2+joff];
               t3 = t1*f[1+3*j2+joff];
               t4 = t1*f[2+3*j2+joff];
               f[3*j2+joff] = f[3*j1+joff] - t2;
               f[1+3*j2+joff] = f[1+3*j1+joff] - t3;
               f[2+3*j2+joff] = f[2+3*j1+joff] - t4;
               f[3*j1+joff] += t2;
               f[1+3*j1+joff] += t3;
               f[2+3*j1+joff] += t4;
            }
         }
         ns = ns2;
      }
/* swap complex components */
      for (j = 0; j < nxh; j++) {
         at1 = crealf(f[2+3*j+joff]);
         f[2+3*j+joff] = cimagf(f[1+3*j+joff])
                         + cimagf(f[2+3*j+joff])*_Complex_I;
         at2 = crealf(f[1+3*j+joff]);
         f[1+3*j+joff] = at1 + cimagf(f[3*j+joff])*_Complex_I;
         f[3*j+joff] = crealf(f[3*j+joff]) + at2*_Complex_I;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rm3y(float complex f[], int isign, int mixup[],
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
         koff = 3*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = 3*nxhd*k1;
            t1 = f[3*i+k1];
            t2 = f[1+3*i+k1];
            t3 = f[2+3*i+k1];
            f[3*i+k1] = f[3*i+koff];
            f[1+3*i+k1] = f[1+3*i+koff];
            f[2+3*i+k1] = f[2+3*i+koff];
            f[3*i+koff] = t1;
            f[1+3*i+koff] = t2;
            f[2+3*i+koff] = t3;
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
               j1 = 3*nxhd*(j + k1);
               j2 = 3*nxhd*(j + k2);
               t1 = sct[kmr*j];
               t2 = t1*f[3*i+j2];
               t3 = t1*f[1+3*i+j2];
               t4 = t1*f[2+3*i+j2];
               f[3*i+j2] = f[3*i+j1] - t2;
               f[1+3*i+j2] = f[1+3*i+j1] - t3;
               f[2+3*i+j2] = f[2+3*i+j1] - t4;
               f[3*i+j1] += t2;
               f[1+3*i+j1] += t3;
               f[2+3*i+j1] += t4;
            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if (nxi==1) {
      for (k = 1; k < nyh; k++) {
         koff = 3*nxhd*k;
         k1 = 3*nxhd*ny - koff;
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
         koff = 3*nxhd*k;
         k1 = 3*nxhd*ny - koff;
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
         koff = 3*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = 3*nxhd*k1;
            t1 = f[3*i+k1];
            t2 = f[1+3*i+k1];
            t3 = f[2+3*i+k1];
            f[3*i+k1] = f[3*i+koff];
            f[1+3*i+k1] = f[1+3*i+koff];
            f[2+3*i+k1] = f[2+3*i+koff];
            f[3*i+koff] = t1;
            f[1+3*i+koff] = t2;
            f[2+3*i+koff] = t3;
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
               j1 = 3*nxhd*(j + k1);
               j2 = 3*nxhd*(j + k2);
               t1 = conjf(sct[kmr*j]);
               t2 = t1*f[3*i+j2];
               t3 = t1*f[1+3*i+j2];
               t4 = t1*f[2+3*i+j2];
               f[3*i+j2] = f[3*i+j1] - t2;
               f[1+3*i+j2] = f[1+3*i+j1] - t3;
               f[2+3*i+j2] = f[2+3*i+j1] - t4;
               f[3*i+j1] += t2;
               f[1+3*i+j1] += t3;
               f[2+3*i+j1] += t4;
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cmswapc2n(float f[], float s[], int isign, int nxh, int nyi,
               int nyt, int nxhd, int nyd, int ndim) {
/* this subroutine swaps components for multiple ffts
   f = input  array
   s = scratch array
   isign = (-1,1) = swap (real-to-complex,complex-to-real)
   nxh = complex dimension in x direction
   nyi/nyt = initial/final y index used
   nxhd = half of the second dimension of f
   nyd = third dimension of f
   ndim = leading dimension of array f
local data                                                            */
   int i, j, k, ioff, nk;
/* swap complex components */
/* real to complex */
   if (isign < 0){
#pragma omp parallel for private(i,j,k,ioff)
      for (k = nyi-1; k < nyt; k++) {
         nk = 2*ndim*nxhd*k;
         for (j = 0; j < nxh; j++) {
            ioff = 2*ndim*j;
            for (i = 0; i < ndim; i++) {
               s[2*i+ioff+nk] = f[i+ndim*(2*j)+nk];
               s[2*i+ioff+1+nk] = f[i+ndim*(2*j+1)+nk];
            }
         }
         for (j = 0; j < nxh; j++) {
            ioff = 2*ndim*j;
            for (i = 0; i < ndim; i++) {
               f[i+ndim*(2*j)+nk] = s[i+ioff+nk];
            }
            ioff += ndim;
            for (i = 0; i < ndim; i++) {
               f[i+ndim*(2*j+1)+nk] = s[i+ioff+nk];
            }
         }
      }
   }
/* complex to real */
   else if (isign > 0) {
#pragma omp parallel for private(i,j,k,ioff)
      for (k = nyi-1; k < nyt; k++) {
         nk = 2*ndim*nxhd*k;
         for (j = 0; j < nxh; j++) {
            ioff = 2*ndim*j;
            for (i = 0; i < ndim; i++) {
               s[i+ioff+nk] = f[i+ndim*(2*j)+nk];
            }
            ioff += ndim;
            for (i = 0; i < ndim; i++) {
               s[i+ioff+nk] = f[i+ndim*(2*j+1)+nk];
            }
         }
         for (j = 0; j < nxh; j++) {
            ioff = 2*ndim*j;
            for (i = 0; i < ndim; i++) {
               f[i+ndim*(2*j)+nk] = s[2*i+ioff+nk];
               f[i+ndim*(2*j+1)+nk] = s[2*i+ioff+1+nk];
            }
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmnx(float complex f[], float complex ss[], int isign,
               int mixup[], float complex sct[], int indx, int indy,
               int nyi, int nyp, int nxhd, int nyd, int ndim, int nxhyd,
               int nxyhd) {
/* this subroutine performs the x part of N two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   y, using complex arithmetic, where N = ndim, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)
   for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)
   where M = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, two inverse fourier transforms are performed
   f[m][n][0:N-1] = (1/nx*ny)*sum(f[k][j][0:N-1]*
         exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, two forward fourier transforms are performed
   f[k][j][0:N-1] = sum(f[m][n][0:N-1]*exp(sqrt(-1)*2pi*n*j/nx)*
         exp(sqrt(-1)*2pi*m*k/ny))
   ss = scratch array
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = second dimension of f
   nyd = third dimension of f
   ndim = leading dimension of array f
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j][0:N-1], = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][0][0:N-1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0][0:N-1]) = real part of mode nx/2,0 and
   imagf(f[ny/2][0][0:N-1]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
   int nrxb;
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
      goto L110;
/* inverse fourier transform */
   nrxb = nxhy/nxh;
   nrx = nxy/nxh;
/* swap complex components */
   cmswapc2n((float *)f,(float *)ss,isign,nxh,nyi,nyt,nxhd,nyd,ndim);
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,joff,ani,t1,t2,t3)
   for (i = nyi-1; i < nyt; i++) {
      joff = ndim*nxhd*i;
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            for (jj = 0; jj < ndim; jj++) {
               t1 = f[jj+ndim*j1+joff];
               f[jj+ndim*j1+joff] = f[jj+ndim*j+joff];
               f[jj+ndim*j+joff] = t1;
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
               j1 = j + k1;
               j2 = j + k2;
               t1 = sct[kmr*j];
               for (jj = 0; jj < ndim; jj++) {
                  t2 = t1*f[jj+ndim*j2+joff];
                  f[jj+ndim*j2+joff] = f[jj+ndim*j1+joff] - t2;
                  f[jj+ndim*j1+joff] += t2;
               }

            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxy/nx;
      ani = 1.0/(float) (2*nx*ny);
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         for (jj = 0; jj < ndim; jj++) {
            t2 = conjf(f[jj+ndim*(nxh-j)+joff]);
            t1 = f[jj+ndim*j+joff] + t2;
            t2 = (f[jj+ndim*j+joff] - t2)*t3;
            f[jj+ndim*j+joff] = ani*(t1 + t2);
            f[jj+ndim*(nxh-j)+joff] = ani*conjf(t1 - t2);
         }
      }
      ani = 2.0*ani;
      for (jj = 0; jj < ndim; jj++) {
         f[jj+ndim*nxhh+joff] = ani*conjf(f[jj+ndim*nxhh+joff]);
         f[jj+joff] = ani*((crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                    + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
L110: nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,joff,t1,t2,t3)
   for (i = nyi-1; i < nyt; i++) {
      joff = ndim*nxhd*i;
/* scramble coefficients */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         for (jj = 0; jj < ndim; jj++) {
            t2 = conjf(f[jj+ndim*(nxh-j)+joff]);
            t1 = f[jj+ndim*j+joff] + t2;
            t2 = (f[jj+ndim*j+joff] - t2)*t3;
            f[jj+ndim*j+joff] = t1 + t2;
            f[jj+ndim*(nxh-j)+joff] = conjf(t1 - t2);
         }
      }
      for (jj = 0; jj < ndim; jj++) {
         f[jj+ndim*nxhh+joff] = 2.0*conjf(f[jj+ndim*nxhh+joff]);
         f[jj+joff] = (crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                    + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I;
      }

/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            for (jj = 0; jj < ndim; jj++) {
               t1 = f[jj+ndim*j1+joff];
               f[jj+ndim*j1+joff] = f[jj+ndim*j+joff];
               f[jj+ndim*j+joff] = t1;
            }
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
               j1 = j + k1;
               j2 = j + k2;
               t1 = conjf(sct[kmr*j]);
               for (jj = 0; jj < ndim; jj++) {
                  t2 = t1*f[jj+ndim*j2+joff];
                  f[jj+ndim*j2+joff] = f[jj+ndim*j1+joff] - t2;
                  f[jj+ndim*j1+joff] += t2;
               }

            }
         }
         ns = ns2;
      }
   }
/* swap complex components */
   cmswapc2n((float *)f,(float *)ss,isign,nxh,nyi,nyt,nxhd,nyd,ndim);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmny(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int ndim, int nxhyd,
               int nxyhd) {
/* this subroutine performs the y part of N two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   x, using complex arithmetic, where N = ndim, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)
   for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)
   where M = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, two inverse fourier transforms are performed
   f[m][n][0:N-1] = (1/nx*ny)*sum(f[k][j][0:N-1] *
         exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, two forward fourier transforms are performed
   f[k][j][0:N-1] = sum(f[m][n][0:N-1]*exp(sqrt(-1)*2pi*n*j/nx)*
         exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxi = initial x index used
   nxp = number of x indices used
   nxhd = second dimension of f
   nyd = third dimension of f
   ndim = leading dimension of array f
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j][0:N-1] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][0][0:N-1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0][0:N-1]) = real part of mode nx/2,0 and
   imag(f[ny/2][0][0:N-1]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
   int nryb;
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
      goto L100;
/* inverse fourier transform */
   nryb = nxhy/ny;
   nry = nxy/ny;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,joff,t1,t2)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         joff = ndim*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = ndim*nxhd*k1;
            for (jj = 0; jj < ndim; jj++) {
               t1 = f[jj+ndim*i+k1];
               f[jj+ndim*i+k1] = f[jj+ndim*i+joff];
               f[jj+ndim*i+joff] = t1;
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
               j1 = ndim*nxhd*(j + k1);
               j2 = ndim*nxhd*(j + k2);
               t1 = sct[kmr*j];

               for (jj = 0; jj < ndim; jj++) {
                  t2 = t1*f[jj+ndim*i+j2];
                  f[jj+ndim*i+j2] = f[jj+ndim*i+j1] - t2;
                  f[jj+ndim*i+j1] += t2;
               }
            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      if (nxi==1) {
         joff = ndim*nxhd*k;
         k1 = ndim*nxhd*ny - joff;
         for (jj = 0; jj < ndim; jj++) {
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
L100: nryb = nxhy/ny;
   nry = nxy/ny;
/* scramble modes kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      if (nxi==1) {
         joff = ndim*nxhd*k;
         k1 = ndim*nxhd*ny - joff;
         for (jj = 0; jj < ndim; jj++) {
            t1 = cimagf(f[jj+k1]) + crealf(f[jj+k1])*_Complex_I;
            f[jj+k1] = conjf(f[jj+joff] - t1);
            f[jj+joff] += t1;
         }
      }
   }
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,joff,t1,t2)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         joff = ndim*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = ndim*nxhd*k1;
            for (jj = 0; jj < ndim; jj++) {
               t1 = f[jj+ndim*i+k1];
               f[jj+ndim*i+k1] = f[jj+ndim*i+joff];
               f[jj+ndim*i+joff] = t1;
            }
         }
      }
/* first transform in y */
      ns = 1;
      for (l = 0; l < indy; l++) {
         ns2 = ns + ns;
         km = nyh/ns;
         kmr = km*nry;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = ndim*nxhd*(j + k1);
               j2 = ndim*nxhd*(j + k2);
               t1 = conjf(sct[kmr*j]);
               for (jj = 0; jj < ndim; jj++) {
                  t2 = t1*f[jj+ndim*i+j2];
                  f[jj+ndim*i+j2] = f[jj+ndim*i+j1] - t2;
                  f[jj+ndim*i+j1] += t2;
               }
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rmx(float complex f[], int isign, int mixup[],
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
      cfft2rmxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
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
      cfft2rmxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                nxyhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rm3(float complex f[], int isign, int mixup[],
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
      cfft2rm3x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                nxyhd);
/* perform y fft */
      cfft2rm3y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
               nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfft2rm3y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                nxyhd);
/* perform x fft */
      cfft2rm3x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                nxyhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rmn(float complex f[], float complex ss[], int isign,
               int mixup[], float complex sct[], int indx, int indy,
               int nxhd, int nyd, int ndim, int nxhyd, int nxyhd) {
/* wrapper function for multiple 2d real to complex ffts */
/* local data */
   int nxh, ny;
   static int nxi = 1, nyi = 1;
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   ny = 1L<<indy;
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cfft2rmnx(f,ss,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,ndim,
               nxhyd,nxyhd);
/* perform y fft */
      cfft2rmny(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,ndim,nxhyd,
               nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfft2rmny(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,ndim,nxhyd,
               nxyhd);
/* perform x fft */
      cfft2rmnx(f,ss,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,ndim,
               nxhyd,nxyhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppcopyout(float part[], float ppart[], int kpic[], int nop,
                int nppmx, int idimp, int mxy1, int *irc) {
/* for 2d code, this subroutine copies segmented particle data ppart to
   the array part with original tiled layout
   input: all except part, output: part
   part[j][i] = i-th coordinate for particle j
   ppart[k][j][i] = i-th coordinate for particle j in tile k
   kpic = number of particles per tilees
   nop = number of particles
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   mxy1 = total number of tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int i, j, k, npoff, npp, ne, ierr;
   npoff = 0;
   ierr = 0;
/* loop over tiles */
   for (k = 0; k < mxy1; k++) {
      npp = kpic[k];
      ne = npp + npoff;
      if (ne > nop)
         ierr = ierr > ne-nop ? ierr : ne-nop;
      if (ierr > 0)
         npp = 0;
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
         for (i = 0; i < idimp; i++) {
            part[i+idimp*(j+npoff)] = ppart[i+idimp*(j+nppmx*k)];
         }
      }
      npoff += npp;
   }
   if (ierr > 0)
      *irc = ierr;
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
void cppmovin2l_(float *part, float *ppart, int *kpic, int *nppmx,
                 int *idimp, int *nop, int *mx, int *my, int *mx1,
                 int *mxy1, int *irc) {
   cppmovin2l(part,ppart,kpic,*nppmx,*idimp,*nop,*mx,*my,*mx1,*mxy1,
              irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcheck2l_(float *ppart, int *kpic, int *idimp, int *nppmx,
                 int *nx, int *ny, int *mx, int *my, int *mx1, int *my1, 
                 int *irc) {
   cppcheck2l(ppart,kpic,*idimp,*nppmx,*nx,*ny,*mx,*my,*mx1,*my1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppush23l_(float *ppart, float *fxy, float *bxy, int *kpic,
                  float *qbm, float *dt, float *dtc, float *ek,
                  int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                  int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                  int *ipbc) {
   cgbppush23l(ppart,fxy,bxy,kpic,*qbm,*dt,*dtc,ek,*idimp,*nppmx,*nx,
               *ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppushf23l_(float *ppart, float *fxy, float *bxy, int *kpic,
                   int *ncl, int *ihole, float *qbm, float *dt,
                   float *dtc, float *ek, int *idimp, int *nppmx,
                   int *nx, int *ny, int *mx, int *my, int *nxv,
                   int *nyv, int *mx1, int *mxy1, int *ntmax, 
                   int *irc) {
   cgbppushf23l(ppart,fxy,bxy,kpic,ncl,ihole,*qbm,*dt,*dtc,ek,*idimp,
                *nppmx,*nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppost2l_(float *ppart, float *q, int *kpic, float *qm,
                int *nppmx, int *idimp, int *mx, int *my, int *nxv,
                int *nyv, int *mx1, int *mxy1) {
   cgppost2l(ppart,q,kpic,*qm,*nppmx,*idimp,*mx,*my,*nxv,*nyv,*mx1,
             *mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cgjppost2l_(float *ppart, float *cu, int *kpic, float *qm,
                 float *dt, int *nppmx, int *idimp, int *nx, int *ny, 
                 int *mx, int *my, int *nxv, int *nyv, int *mx1,
                 int *mxy1, int *ipbc) {
   cgjppost2l(ppart,cu,kpic,*qm,*dt,*nppmx,*idimp,*nx,*ny,*mx,*my,*nxv,
              *nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgmjppost2l_(float *ppart, float *amu, int *kpic, float *qm,
                  int *nppmx, int *idimp, int *mx, int *my, int *nxv,
                  int *nyv, int *mx1, int *mxy1) {
   cgmjppost2l(ppart,amu,kpic,*qm,*nppmx,*idimp,*mx,*my,*nxv,*nyv,*mx1,
               *mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cgdjppost2l_(float *ppart, float *fxy, float *bxy, float *dcu,
                  float *amu, int *kpic, float *qm, float *qbm,
                  float *dt, int *idimp, int *nppmx, int *nx, int *ny,
                  int *mx, int *my, int *nxv, int *nyv, int *mx1,
                  int *mxy1) {
   cgdjppost2l(ppart,fxy,bxy,dcu,amu,kpic,*qm,*qbm,*dt,*idimp,*nppmx,
               *nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cgdcjppost2l_(float *ppart, float *fxy, float *bxy, float *cu,
                   float *dcu, float *amu, int *kpic, float *qm,
                   float *qbm, float *dt, int *idimp, int *nppmx,
                   int *nx, int *ny, int *mx, int *my, int *nxv,
                   int *nyv, int *mx1, int *mxy1) {
   cgdcjppost2l(ppart,fxy,bxy,cu,dcu,amu,kpic,*qm,*qbm,*dt,*idimp,
                *nppmx,*nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cpporder2l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                 int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                 int *mx, int *my, int *mx1, int *my1, int *npbmx,
                 int *ntmax, int *irc) {
   cpporder2l(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*nx,*ny,*mx,*my,
              *mx1,*my1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf2l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                  int *ihole, int *idimp, int *nppmx, int *mx1,
                  int *my1, int *npbmx, int *ntmax, int *irc) {
   cpporderf2l(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*mx1,*my1,
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
void camcguard2l_(float *amu, int *nx, int *ny, int *nxe, int *nye,
                  int *ndim) {
   camcguard2l(amu,*nx,*ny,*nxe,*nye,*ndim);
   return;
}

/*--------------------------------------------------------------------*/
void cascfguard2l_(float *dcu, float *cus, float *q2m0, int *nx,
                   int *ny, int *nxe, int *nye) {
   cascfguard2l(dcu,cus,*q2m0,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void cfwpminmx2_(float *qe, float *qbme, float *wpmax, float *wpmin,
                 int *nx, int *ny, int *nxe, int *nye) {
   cfwpminmx2(qe,*qbme,wpmax,wpmin,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void cmpois23_(float complex *q, float complex *fxy, int *isign,
               float complex *ffc, float *ax, float *ay, float *affp,
               float *we, int *nx, int *ny, int *nxvh, int *nyv,
               int *nxhd, int *nyhd) {
   cmpois23(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*nxvh,*nyv,*nxhd,
            *nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmcuperp2_(float complex *cu, int *nx, int *ny, int *nxvh,
                int *nyv) {
   cmcuperp2(cu,*nx,*ny,*nxvh,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmbbpois23_(float complex *cu, float complex *bxy,
                 float complex *ffc, float *ci, float *wm, int *nx,
                 int *ny, int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   cmbbpois23(cu,bxy,ffc,*ci,wm,*nx,*ny,*nxvh,*nyv,*nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cbaddext2_(float *bxy, float *omx, float *omy, float *omz, int *nx,
                int *ny, int *nxe, int *nye) {
 cbaddext2(bxy,*omx,*omy,*omz,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void cmdcuperp23_(float complex *dcu, float complex *amu, int *nx,
                  int *ny, int *nxvh, int *nyv) {
 cmdcuperp23(dcu,amu,*nx,*ny,*nxvh,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmadcuperp23_(float complex *dcu, float complex *amu, int *nx,
                   int *ny, int *nxvh, int *nyv) {
   cmadcuperp23(dcu,amu,*nx,*ny,*nxvh,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmepois23_(float complex *dcu, float complex *exy, int *isign,
                float complex *ffe, float *ax, float *ay, float *affp,
                float *wp0, float *ci, float *wf, int *nx, int *ny,
                int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   cmepois23(dcu,exy,*isign,ffe,*ax,*ay,*affp,*wp0,*ci,wf,*nx,*ny,*nxvh,
             *nyv,*nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void caddvrfield2_(float *a, float *b, float *c, int *ndim, int *nxe,
                   int *nye) {
   caddvrfield2(a,b,c,*ndim,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                  int *nxhyd, int *nxyhd) {
   cwfft2rinit(mixup,sct,*indx,*indy,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmxx_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nyi,
                int *nyp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd) {
   cfft2rmxx(f,*isign,mixup,sct,*indx,*indy,*nyi,*nyp,*nxhd,*nyd,*nxhyd,
             *nxyhd);
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
void cfft2rm3x_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nyi,
                int *nyp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd) {
   cfft2rm3x(f,*isign,mixup,sct,*indx,*indy,*nyi,*nyp,*nxhd,*nyd,*nxhyd,
             *nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rm3y_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nxi,
                int *nxp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd) {
   cfft2rm3y(f,*isign,mixup,sct,*indx,*indy,*nxi,*nxp,*nxhd,*nyd,*nxhyd,
             *nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmnx_(float complex *f, float complex *ss, int *isign,
               int *mixup, float complex *sct, int *indx, int *indy,
               int *nyi, int *nyp, int *nxhd, int *nyd, int *ndim,
               int *nxhyd, int *nxyhd) {
   cfft2rmnx(f,ss,*isign,mixup,sct,*indx,*indy,*nyi,*nyp,*nxhd,*nyd,
             *ndim,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmny_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nxi,
                int *nxp, int *nxhd, int *nyd, int *ndim, int *nxhyd,
                int *nxyhd) {
   cfft2rmny(f,*isign,mixup,sct,*indx,*indy,*nxi,*nxp,*nxhd,*nyd,*ndim,
             *nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rmx_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nxhd,
                int *nyd, int *nxhyd, int *nxyhd) {
   cwfft2rmx(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rm3_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nxhd,
                int *nyd, int *nxhyd, int *nxyhd) {
   cwfft2rm3(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,*nxyhd);
   return;
}


/*--------------------------------------------------------------------*/
void cwfft2rmn_(float complex *f, float complex *ss, int *isign,
               int *mixup, float complex *sct, int *indx, int *indy,
               int *nxhd, int *nyd, int *ndim, int *nxhyd, int *nxyhd) {
   cwfft2rmn(f,ss,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*ndim,*nxhyd,
             *nxyhd);
   return;
}
