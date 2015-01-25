/* C Library for Skeleton 2D Electrostatic Vector PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <xmmintrin.h>
#include "vpush2.h"

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
double randum() {
/* this is a version of the random number generator dprandom due to
   c. bingham and the yale computer center, producing numbers
   in the interval (0,1).  written for the sun by viktor k. decyk, ucla
local data                                                              */
   static int r1 = 1271199957, r2 = 1013501921;
   static double h1l = 65533.0, h1u = 32767.0;
   int isc, i1;
   double randum, r0, r3, asc, bsc;
   isc = 65536;
   asc = (double) isc;
   bsc = asc*asc;
   i1 = r1 - (r1/isc)*isc;
   r3 = h1l*(double) r1 + asc*h1u*(double) i1;
   i1 = r3/bsc;
   r3 = r3 - ((double) i1)*bsc;
   bsc = 0.5*bsc;
   i1 = r2/isc;
   isc = r2 - i1*isc;
   r0 = h1l*(double) r2 + asc*h1u*(double) isc;
   asc = 1.0/bsc;
   isc = r0*asc;
   r2 = r0 - ((double) isc)*bsc;
   r3 = r3 + ((double) isc + 2.0*h1u*(double) i1);
   isc = r3*asc;
   r1 = r3 - ((double) isc)*bsc;
   randum = ((double) r1 + ((double) r2)*asc)*asc;
   return randum;
}

/*--------------------------------------------------------------------*/
void cdistr2t(float part[], float vtx, float vty, float vdx, float vdy,
              int npx, int npy, int idimp, int npe, int nx, int ny,
              int ipbc) {
/* for 2d code, this subroutine calculates initial particle co-ordinates
   and velocities with uniform density and maxwellian velocity with drift
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = velocity vx of particle n
   part[3][n] = velocity vy of particle n
   vtx/vty = thermal velocity of electrons in x/y direction
   vdx/vdy = drift velocity of beam electrons in x/y direction
   npx/npy = initial number of particles distributed in x/y direction
   idimp = size of phase space = 4
   npe = first dimension of particle array
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
      k1 = npx*k;
      at3 = edgely + at2*(((float) k) + 0.5);
      for (j = 0; j < npx; j++) {
         part[j+k1] = edgelx + at1*(((float) j) + 0.5);
         part[j+k1+npe] = at3;
      }
   }
/* maxwellian velocity distribution */
   for (j = 0; j < npxy; j++) {
      part[j+2*npe] = vtx*ranorm();
      part[j+3*npe] = vty*ranorm();
   }
/* add correct drift */
   dsum1 = 0.0;
   dsum2 = 0.0;
   for (j = 0; j < npxy; j++) {
      dsum1 += part[j+2*npe];
      dsum2 += part[j+3*npe];
   }
   sum1 = dsum1;
   sum2 = dsum2;
   at1 = 1.0/(float) npxy;
   sum1 = at1*sum1 - vdx;
   sum2 = at1*sum2 - vdy;
   for (j = 0; j < npxy; j++) {
      part[j+2*npe] -= sum1;
      part[j+3*npe] -= sum2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cgpush2lt(float part[], float fxy[], float qbm, float dt,
               float *ek, int idimp, int nop, int npe, int nx, int ny,
               int nxv, int nyv, int ipbc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions.
   scalar version using guard cells
   44 flops/particle, 12 loads, 4 stores
   input: all, output: part, ek
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
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = velocity vx of particle n
   part[3][n] = velocity vy of particle n
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   idimp = size of phase space = 4
   nop = number of particles
   npe = first dimension of particle array
   nx/ny = system length in x/y direction
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int j, nn, mm, nxv2;
   float qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
   double sum1;
   nxv2 = 2*nxv;
   qtm = qbm*dt;
   sum1 = 0.0;
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
#pragma ivdep
   for (j = 0; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      nn = x;
      mm = y;
      dxp = x - (float) nn;
      dyp = y - (float) mm;
      nn = 2*nn + nxv2*mm;
      amx = 1.0f - dxp;
      amy = 1.0f - dyp;
/* find acceleration */
      dx = amx*fxy[nn];
      dy = amx*fxy[nn+1];
      dx = amy*(dxp*fxy[nn+2] + dx);
      dy = amy*(dxp*fxy[nn+3] + dy);
      nn += nxv2;
      vx = amx*fxy[nn];
      vy = amx*fxy[nn+1];
      dx += dyp*(dxp*fxy[nn+2] + vx);
      dy += dyp*(dxp*fxy[nn+3] + vy);
/* new velocity */
      dxp = part[j+2*npe];
      dyp = part[j+3*npe];
      vx = dxp + qtm*dx;
      vy = dyp + qtm*dy;
/* average kinetic energy */
      dxp += vx;
      dyp += vy;
      sum1 += dxp*dxp + dyp*dyp;
/* new position */
      dx = x + vx*dt;
      dy = y + vy*dt;
/* periodic boundary conditions */
      if (ipbc==1) {
         if (dx < edgelx) dx += edgerx;
         if (dx >= edgerx) dx -= edgerx;
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
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
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* set new position */
      part[j] = dx;
      part[j+npe] = dy;
/* set new velocity */
      part[j+2*npe] = vx;
      part[j+3*npe] = vy;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum1;
   return;
}

/*--------------------------------------------------------------------*/
void cvgpush2lt(float part[], float fxy[], float qbm, float dt,
                float *ek, int idimp, int nop, int npe, int nx, int ny,
                int nxv, int nyv, int ipbc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions.
   vectorizable version using guard cells
   44 flops/particle, 12 loads, 4 stores
   input: all, output: part, ek
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
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = velocity vx of particle n
   part[3][n] = velocity vy of particle n
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   idimp = size of phase space = 4
   nop = number of particles
   npe = first dimension of particle array
   nx/ny = system length in x/y direction
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define NPBLK             32
#define LVECT             4
   int i, j, k, ipp, joff, nps, nn, mm;
   float qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
/* scratch arrays */
   int n[NPBLK];
   float s[NPBLK*LVECT], t[NPBLK*2];
   double sum1;
   qtm = qbm*dt;
   sum1 = 0.0;
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
   ipp = nop/NPBLK;
/* outer loop over number of full blocks */
   for (k = 0; k < ipp; k++) {
      joff = NPBLK*k;
/* inner loop over particles in block */
      for (j = 0; j < NPBLK; j++) {
/* find interpolation weights */
         x = part[j+joff];
         y = part[j+joff+npe];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         n[j] = nn + nxv*mm;
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
         mm = nn + nxv - 2;
         dx = 0.0f;
         dy = 0.0f;
#pragma ivdep
         for (i = 0; i < LVECT; i++) {
            if (i > 1)
               nn = mm;
            dx += fxy[2*(i+nn)]*s[j+NPBLK*i];
            dy += fxy[1+2*(i+nn)]*s[j+NPBLK*i];
         }
         s[j] = dx;
         s[j+NPBLK] = dy;
      }
/* new velocity */
      for (j = 0; j < NPBLK; j++) {
         x = t[j];
         y = t[j+NPBLK];
         dxp = part[j+joff+2*npe];
         dyp = part[j+joff+3*npe];
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
/* periodic boundary conditions */
         if (ipbc==1) {
            if (dx < edgelx) dx += edgerx;
            if (dx >= edgerx) dx -= edgerx;
            if (dy < edgely) dy += edgery;
            if (dy >= edgery) dy -= edgery;
         }
/* reflecting boundary conditions */
         else if (ipbc==2) {
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
            if (dy < edgely) dy += edgery;
            if (dy >= edgery) dy -= edgery;
         }
/* set new position */
         part[j+joff] = dx;
         part[j+joff+npe] = dy;
/* set new velocity */
         part[j+joff+2*npe] = vx;
         part[j+joff+3*npe] = vy;
      }
   }
   nps = NPBLK*ipp;
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      nn = x;
      mm = y;
      dxp = x - (float) nn;
      dyp = y - (float) mm;
      nn = 2*(nn + nxv*mm);
      amx = 1.0f - dxp;
      amy = 1.0f - dyp;
/* find acceleration */
      dx = amx*fxy[nn];
      dy = amx*fxy[nn+1];
      dx = amy*(dxp*fxy[nn+2] + dx);
      dy = amy*(dxp*fxy[nn+3] + dy);
      nn += 2*nxv;
      vx = amx*fxy[nn];
      vy = amx*fxy[nn+1];
      dx += dyp*(dxp*fxy[nn+2] + vx);
      dy += dyp*(dxp*fxy[nn+3] + vy);
/* new velocity */
      dxp = part[j+2*npe];
      dyp = part[j+3*npe];
      vx = dxp + qtm*dx;
      vy = dyp + qtm*dy;
/* average kinetic energy */
      dxp += vx;
      dyp += vy;
      sum1 += dxp*dxp + dyp*dyp;
/* new position */
      dx = x + vx*dt;
      dy = y + vy*dt;
/* periodic boundary conditions */
      if (ipbc==1) {
         if (dx < edgelx) dx += edgerx;
         if (dx >= edgerx) dx -= edgerx;
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
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
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* set new position */
      part[j] = dx;
      part[j+npe] = dy;
/* set new velocity */
      part[j+2*npe] = vx;
      part[j+3*npe] = vy;
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum1;
   return;
#undef LVECT
#undef NPBLK
}

/*--------------------------------------------------------------------*/
void cvgpost2lt(float part[], float q[], float qm, int nop, int npe,
                int idimp, int nxv, int nyv) {
/* for 2d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   vectorizable version using guard cells
   17 flops/particle, 6 loads, 4 stores
   input: all, output: q
   charge density is approximated by values at the nearest grid points
   q(n,m)=qm*(1.-dx)*(1.-dy)
   q(n+1,m)=qm*dx*(1.-dy)
   q(n,m+1)=qm*(1.-dx)*dy
   q(n+1,m+1)=qm*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   q[k][j] = charge density at grid point j,k
   qm = charge on particle, in units of e
   nop = number of particles
   npe = first dimension of particle array
   idimp = size of phase space = 4
   nxv = first dimension of charge array, must be >= nx+1
   nyv = second dimension of charge array, must be >= ny+1
local data                                                            */
#define NPBLK             32
#define LVECT             4
   int i, j, k, ipp, joff, nps, nn, mm;
   float x, y, dxp, dyp, amx, amy;
/* scratch arrays */
   int n[NPBLK];
   float s[NPBLK*LVECT];
   ipp = nop/NPBLK;
/* outer loop over number of full blocks */
   for (k = 0; k < ipp; k++) {
      joff = NPBLK*k;
/* inner loop over particles in block */
      for (j = 0; j < NPBLK; j++) {
/* find interpolation weights */
         x = part[j+joff];
         y = part[j+joff+npe];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         n[j] = nn + nxv*mm;
         amx = qm - dxp;
         amy = 1.0f - dyp;
         s[j] = amx*amy;
         s[j+NPBLK] = dxp*amy;
         s[j+2*NPBLK] = amx*dyp;
         s[j+3*NPBLK] = dxp*dyp;
     }
/* deposit charge */
      for (j = 0; j < NPBLK; j++) {
         nn = n[j];
         mm = nn + nxv - 2;
#pragma ivdep
         for (i = 0; i < LVECT; i++) {
            if (i > 1)
               nn = mm;
            q[i+nn] += s[j+NPBLK*i];
         }
      }
   }
   nps = NPBLK*ipp;
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      nn = x;
      mm = y;
      dxp = qm*(x - (float) nn);
      dyp = y - (float) mm;
      nn = nn + nxv*mm;
      amx = qm - dxp;
      amy = 1.0f - dyp;
/* deposit charge */
      x = q[nn] + amx*amy;
      y = q[nn+1] + dxp*amy;
      q[nn] = x;
      q[nn+1] = y;
      nn += nxv;
      x = q[nn] + amx*dyp;
      y = q[nn+1] + dxp*dyp;
      q[nn] = x;
      q[nn+1] = y;
   }
   return;
#undef LVECT
#undef NPBLK
}

/*--------------------------------------------------------------------*/
void cgpost2lt(float part[], float q[], float qm, int nop, int npe,
               int idimp, int nxv, int nyv) {
/* for 2d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   scalar version using guard cells
   17 flops/particle, 6 loads, 4 stores
   input: all, output: q
   charge density is approximated by values at the nearest grid points
   q(n,m)=qm*(1.-dx)*(1.-dy)
   q(n+1,m)=qm*dx*(1.-dy)
   q(n,m+1)=qm*(1.-dx)*dy
   q(n+1,m+1)=qm*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   q[k][j] = charge density at grid point j,k
   qm = charge on particle, in units of e
   nop = number of particles
   npe = first dimension of particle array
   idimp = size of phase space = 4
   nxv = first dimension of charge array, must be >= nx+1
   nyv = second dimension of charge array, must be >= ny+1
local data                                                            */
   int j, nn, mm;
   float x, y, dxp, dyp, amx, amy;
   for (j = 0; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      nn = x;
      mm = y;
      dxp = qm*(x - (float) nn);
      dyp = y - (float) mm;
      nn = nn + nxv*mm;
      amx = qm - dxp;
      amy = 1.0f - dyp;
/* deposit charge */
      x = q[nn] + amx*amy;
      y = q[nn+1] + dxp*amy;
      q[nn] = x;
      q[nn+1] = y;
      nn += nxv;
      x = q[nn] + amx*dyp;
      y = q[nn+1] + dxp*dyp;
      q[nn] = x;
      q[nn+1] = y;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cdsortp2ylt(float parta[], float partb[], int npic[], int idimp,
                 int nop, int npe, int ny1) {
/* this subroutine sorts particles by y grid
   linear interpolation
   parta/partb = input/output particle arrays
   parta[1][n] = position y of particle n
   npic = address offset for reordering particles
   idimp = size of phase space = 4
   nop = number of particles
   npe = first dimension of particle array
   ny1 = system length in y direction + 1
local data                                                            */
   int i, j, k, m, isum, ist, ip;
/* clear counter array */
   for (k = 0; k < ny1; k++) {
      npic[k] = 0;
   }
/* find how many particles in each grid */
   for (j = 0; j < nop; j++) {
      m = parta[j+npe];
      npic[m] += 1;
   }
/* find address offset */
   isum = 0;
   for (k = 0; k < ny1; k++) {
      ist = npic[k];
      npic[k] = isum;
      isum += ist;
   }
/* find addresses of particles at each grid and reorder particles */
   for (j = 0; j < nop; j++) {
      m = parta[j+npe];
      ip = npic[m];
      npic[m] = ip + 1;
      for (i = 0; i < idimp; i++) {
         partb[ip+npe*i] = parta[j+npe*i];
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void ccguard2l(float fxy[], int nx, int ny, int nxe, int nye) {
/* replicate extended periodic vector field fxy
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
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
   nye = second dimension of field arrays, must be >= ny+1
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
void cvpois22(float complex q[], float complex fxy[], int isign,
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
   g[ky][kx] = (affp/(kx**2+ky**2))*s[ky][kx],
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
   double wp;
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
L30: wp = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
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
   }
/* mode numbers kx = 0, nx/2 */
#pragma ivdep
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      at1 = crealf(ffc[kk])*cimagf(ffc[kk]);
      at3 = at1*dny*(float) k;
      zt1 = cimagf(q[kj]) - crealf(q[kj])*_Complex_I;
      fxy[2*kj] = zero;
      fxy[1+2*kj] = at3*zt1;
      fxy[2*k1] = zero;
      fxy[1+2*k1] = zero;
      at1 = at1*(q[kj]*conjf(q[kj]));
      wp += (double) at1;
   }
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

/*--------------------------------------------------------------------*/
void cfft2rvxx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the x part of a two dimensional real to
   complex fast fourier transform and its inverse, for a subset of y,
   using complex arithmetic.
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, an inverse fourier transform is performed
   f[m][n] = (1/nx*ny)*sum(f[k][j]*
         exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, a forward fourier transform is performed
   f[k][j] = sum(f[m][n]*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = first dimension of f >= nx/2
   nyd = second dimension of f >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][2*j],f[k][2*j+1] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k <= ny, except for
   f[k][0],f[k][1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   f[1][0] = real part of mode nx/2,0 and
   f[1][ny/2] = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, j1, k1, k2, ns, ns2, km, kmr, joff;
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
      goto L100;
/* inverse fourier transform */
/* bit-reverse array elements in x */
   nrx = nxhy/nxh;
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = nyi-1; k < nyt; k++) {
         joff = nxhd*k;
         t1 = f[j1+joff];
         f[j1+joff] = f[j+joff];
         f[j+joff] = t1;
      }
   }
/* first transform in x */
   nrx = nxy/nxh;
   ns = 1;
   for (l = 0; l < indx1; l++) {
      ns2 = ns + ns;
      km = nxhh/ns;
      kmr = km*nrx;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (i = nyi-1; i < nyt; i++) {
            joff = nxhd*i;
            for (j = 0; j < ns; j++) {
               t1 = sct[kmr*j];
               t2 = t1*f[j+k2+joff];
               f[j+k2+joff] = f[j+k1+joff] - t2;
               f[j+k1+joff] += t2;
            }
         }
      }
      ns = ns2;
   }
/* unscramble coefficients and normalize */
   kmr = nxy/nx;
   ani = 1.0/(float) (2*nx*ny);
   for (k = nyi-1; k < nyt; k++) {
      joff = nxhd*k;
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         t2 = conjf(f[nxh-j+joff]);
         t1 = f[j+joff] + t2;
         t2 = (f[j+joff] - t2)*t3;
         f[j+joff] = ani*(t1 + t2);
         f[nxh-j+joff] = ani*conjf(t1 - t2);
      }
   }
   ani = 2.0*ani;
   for (k = nyi-1; k < nyt; k++) {
      joff = nxhd*k;
      f[nxhh+joff] = ani*conjf(f[nxhh+joff]);
      f[joff] = ani*((crealf(f[joff]) + cimagf(f[joff]))
                + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I);
   }
   return;
/* forward fourier transform */
/* scramble coefficients */
L100: kmr = nxy/nx;
   for (k = nyi-1; k < nyt; k++) {
      joff = nxhd*k;
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         t2 = conjf(f[nxh-j+joff]);
         t1 = f[j+joff] + t2;
         t2 = (f[j+joff] - t2)*t3;
         f[j+joff] = t1 + t2;
         f[nxh-j+joff] = conjf(t1 - t2);
      }
   }
   for (k = nyi-1; k < nyt; k++) {
      joff = nxhd*k;
      f[nxhh+joff] = 2.0*conjf(f[nxhh+joff]);
      f[joff] = (crealf(f[joff]) + cimagf(f[joff]))
                + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I;
   }
/* bit-reverse array elements in x */
   nrx = nxhy/nxh;
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = nyi-1; k < nyt; k++) {
         joff = nxhd*k;
         t1 = f[j1+joff];
         f[j1+joff] = f[j+joff];
         f[j+joff] = t1;
      }
   }
/* then transform in x */
   nrx = nxy/nxh;
   ns = 1;
   for (l = 0; l < indx1; l++) {
      ns2 = ns + ns;
      km = nxhh/ns;
      kmr = km*nrx;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (i = nyi-1; i < nyt; i++) {
            joff = nxhd*i;
            for (j = 0; j < ns; j++) {
               t1 = conjf(sct[kmr*j]);
               t2 = t1*f[j+k2+joff];
               f[j+k2+joff] = f[j+k1+joff] - t2;
               f[j+k1+joff] += t2;
            }
         }
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rxy(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the y part of a two dimensional real to
   complex fast fourier transform and its inverse, for a subset of x,
   using complex arithmetic
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, an inverse fourier transform is performed
   f[m][n] = (1/nx*ny)*sum(f[k][j]*
         exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, a forward fourier transform is performed
   f[k][j] = sum(f[m][n]*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxi = initial x index used
   nxp = number of x indices used
   nxhd = first dimension of f >= nx/2
   nyd = second dimension of f >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][2*j],f[k][2*j+1] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][0],f[k][1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   f[1][0] = real part of mode nx/2,0 and
   f[1][ny/2] = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
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
      goto L80;
/* inverse fourier transform */
   nry = nxhy/ny;
/* bit-reverse array elements in y */
   for (k = 0; k < ny; k++) {
      joff = nxhd*k;
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      k1 = nxhd*k1;
      for (j = nxi-1; j < nxt; j++) {
         t1 = f[j+k1];
         f[j+k1] = f[j+joff];
         f[j+joff] = t1;
      }
   }
/* then transform in y */
   nry = nxy/ny;
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
            for (i = nxi-1; i < nxt; i++) {
               t2 = t1*f[i+j2];
               f[i+j2] = f[i+j1] - t2;
               f[i+j1] += t2;
            }
         }
      }
      ns = ns2;
   }
/* unscramble modes kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      if (nxi==1) {
         joff = nxhd*k;
         k1 = nxhd*ny - joff;
         t1 = f[k1];
         f[k1] = 0.5*(cimagf(f[joff] + t1)
                  + crealf(f[joff] - t1)*_Complex_I);
         f[joff] = 0.5*(crealf(f[joff] + t1)
                    + cimagf(f[joff] - t1)*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
/* scramble modes kx = 0, nx/2 */
L80: for (k = 1; k < nyh; k++) {
      if (nxi==1) {
         joff = nxhd*k;
         k1 = nxhd*ny - joff;
         t1 = cimagf(f[k1]) + crealf(f[k1])*_Complex_I;
         f[k1] = conjf(f[joff] - t1);
         f[joff] += t1;
      }
   }
/* bit-reverse array elements in y */
   nry = nxhy/ny;
   for (k = 0; k < ny; k++) {
      joff = nxhd*k;
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      k1 = nxhd*k1;
      for (j = nxi-1; j < nxt; j++) {
         t1 = f[j+k1];
         f[j+k1] = f[j+joff];
         f[j+joff] = t1;
      }
   }
/* first transform in y */
   nry = nxy/ny;
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
            for (i = nxi-1; i < nxt; i++) {
               t2 = t1*f[i+j2];
               f[i+j2] = f[i+j1] - t2;
               f[i+j1] += t2;
            }
         }
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rv2x(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the x part of 2 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   y, using complex arithmetic
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, two inverse fourier transforms are performed
   f[m][n][0:1] = (1/nx*ny)*sum(f[k][j][0:1]*
         exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, two forward fourier transforms are performed
   f[k][j][0:1] = sum(f[m][n][0:1]*exp(sqrt(-1)*2pi*n*j/nx)*
         exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = second dimension of f
   nyd = third dimension of f
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][2*j],f[k][2*j+1] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][0],f[k][1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   f[1][0] = real part of mode nx/2,0 and
   f[1][ny/2] = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, jj, j1, k1, k2, ns, ns2, km, kmr, joff;
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
      goto L140;
/* inverse fourier transform */
/* swap complex components */
   for (k = nyi-1; k < nyt; k++) {
      for (j = 0; j < nxh; j++) {
         joff = 2*nxhd*k;
         at1 = cimagf(f[2*j+joff]);
         f[2*j+joff] = crealf(f[2*j+joff])
                       + crealf(f[1+2*j+joff])*_Complex_I;
         f[1+2*j+joff] = at1 + cimagf(f[1+2*j+joff])*_Complex_I;
       }
   }
/* bit-reverse array elements in x */
   nrx = nxhy/nxh;
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = nyi-1; k < nyt; k++) {
         joff = 2*nxhd*k;
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
   for (l = 0; l < indx1; l++) {
      ns2 = ns + ns;
      km = nxhh/ns;
      kmr = km*nrx;
      for (k = 0; k < km; k++) {
         k1 = 2*ns2*k;
         k2 = k1 + 2*ns;
         for (i = nyi-1; i < nyt; i++) {
            joff = 2*nxhd*i;
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
      }
      ns = ns2;
   }
/* unscramble coefficients and normalize */
   kmr = nxy/nx;
   ani = 1.0/(float) (2*nx*ny);
   for (k = nyi-1; k < nyt; k++) {
      joff = 2*nxhd*k;
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
   }
   ani = 2.0*ani;
   for (k = nyi-1; k < nyt; k++) {
      joff = 2*nxhd*k;
      for (jj = 0; jj < 2; jj++) {
         f[jj+2*nxhh+joff] = ani*conjf(f[jj+2*nxhh+joff]);
         f[jj+joff] = ani*((crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
/* scramble coefficients */
L140: kmr = nxy/nx;
   for (k = nyi-1; k < nyt; k++) {
      joff = 2*nxhd*k;
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
   }
   for (k = nyi-1; k < nyt; k++) {
      joff = 2*nxhd*k;
      for (jj = 0; jj < 2; jj++) {
         f[jj+2*nxhh+joff] = 2.0*conjf(f[jj+2*nxhh+joff]);
         f[jj+joff] = (crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I;
      }
   }
/* bit-reverse array elements in x */
   nrx = nxhy/nxh;
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = nyi-1; k < nyt; k++) {
         joff = 2*nxhd*k;
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
   for (l = 0; l < indx1; l++) {
      ns2 = ns + ns;
      km = nxhh/ns;
      kmr = km*nrx;
      for (k = 0; k < km; k++) {
         k1 = 2*ns2*k;
         k2 = k1 + 2*ns;
         for (i = nyi-1; i < nyt; i++) {
            joff = 2*nxhd*i;
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
      }
      ns = ns2;
   }
/* swap complex components */
   for (k = nyi-1; k < nyt; k++) {
      joff = 2*nxhd*k;
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
void cfft2r2y(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the y part of 2 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   x, using complex arithmetic
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, two inverse fourier transforms are performed
   f[m][n][0:1] = (1/nx*ny)*sum(f[k][j][0:1] *
         exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, two forward fourier transforms are performed
   f[k][j][0:1] = sum(f[m][n][0:1]*exp(sqrt(-1)*2pi*n*j/nx)*
         exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxi = initial x index used
   nxp = number of x indices used
   nxhd = second dimension of f
   nyd = third dimension of f
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][2*j],f[k][2*j+1] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][0],f[k][1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   f[1][0] = real part of mode nx/2,0 and
   f[1][ny/2] = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
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
      goto L90;
/* inverse fourier transform */
   nry = nxhy/ny;
/* bit-reverse array elements in y */
   for (k = 0; k < ny; k++) {
      joff = 2*nxhd*k;
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      k1 = 2*nxhd*k1;
      for (j = nxi-1; j < nxt; j++) {
         t1 = f[2*j+k1];
         t2 = f[1+2*j+k1];
         f[2*j+k1] = f[2*j+joff];
         f[1+2*j+k1] = f[1+2*j+joff];
         f[2*j+joff] = t1;
         f[1+2*j+joff] = t2;
      }
   }
/* then transform in y */
   nry = nxy/ny;
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
            for (i = nxi-1; i < nxt; i++) {
               t2 = t1*f[2*i+j2];
               t3 = t1*f[1+2*i+j2];
               f[2*i+j2] = f[2*i+j1] - t2;
               f[1+2*i+j2] = f[1+2*i+j1] - t3;
               f[2*i+j1] += t2;
               f[1+2*i+j1] += t3;
            }
         }
      }
      ns = ns2;
   }
/* unscramble modes kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      if (nxi==1) {
         joff = 2*nxhd*k;
         k1 = 2*nxhd*ny - joff;
         for (jj = 0; jj < 2; jj++) {
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
/* scramble modes kx = 0, nx/2 */
L90: for (k = 1; k < nyh; k++) {
      if (nxi==1) {
         joff = 2*nxhd*k;
         k1 = 2*nxhd*ny - joff;
         for (jj = 0; jj < 2; jj++) {
            t1 = cimagf(f[jj+k1]) + crealf(f[jj+k1])*_Complex_I;
            f[jj+k1] = conjf(f[jj+joff] - t1);
            f[jj+joff] += t1;
         }
      }
   }
/* bit-reverse array elements in y */
   nry = nxhy/ny;
   for (k = 0; k < ny; k++) {
      joff = 2*nxhd*k;
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      k1 = 2*nxhd*k1;
      for (j = nxi-1; j < nxt; j++) {
         t1 = f[2*j+k1];
         t2 = f[1+2*j+k1];
         f[2*j+k1] = f[2*j+joff];
         f[1+2*j+k1] = f[1+2*j+joff];
         f[2*j+joff] = t1;
         f[1+2*j+joff] = t2;
      }
   }
/* first transform in y */
   nry = nxy/ny;
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
            for (i = nxi-1; i < nxt; i++) {
               t2 = t1*f[2*i+j2];
               t3 = t1*f[1+2*i+j2];
               f[2*i+j2] = f[2*i+j1] - t2;
               f[1+2*i+j2] = f[1+2*i+j1] - t3;
               f[2*i+j1] += t2;
               f[1+2*i+j1] += t3;
            }
         }
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rvx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd, int nyd,
               int nxhyd, int nxyhd) {
/* wrapper function for real to complex fft, with packed data */
/* local data */
   int nxh, ny;
   static int nxi = 1, nyi = 1;
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   ny = 1L<<indy;
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cfft2rvxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                nxyhd);
/* perform y fft */
      cfft2rxy(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfft2rxy(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,nxyhd);
/* perform x fft */
      cfft2rvxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                nxyhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rv2(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd) {
/* wrapper function for 2 2d real to complex ffts */
/* local data */
   int nxh, ny;
   static int nxi = 1, nyi = 1;
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   ny = 1L<<indy;
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cfft2rv2x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                nxyhd);
/* perform y fft */
      cfft2r2y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfft2r2y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                nxyhd);
/* perform x fft */
      cfft2rv2x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                nxyhd);
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void cdistr2t_(float *part, float *vtx, float *vty, float *vdx,
               float *vdy, int *npx, int *npy, int *idimp, int *npe,
               int *nx, int *ny, int *ipbc) {
   cdistr2t(part,*vtx,*vty,*vdx,*vdy,*npx,*npy,*idimp,*npe,*nx,*ny,
            *ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpush2lt_(float *part, float *fxy, float *qbm, float *dt,
                float *ek, int *idimp, int *nop, int *npe, int *nx,
                int *ny, int *nxv, int *nyv, int *ipbc) {
   cgpush2lt(part,fxy,*qbm,*dt,ek,*idimp,*nop,*npe,*nx,*ny,*nxv,*nyv,
             *ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgpush2lt_(float *part, float *fxy, float *qbm, float *dt,
                 float *ek, int *idimp, int *nop, int *npe, int *nx,
                 int *ny, int *nxv, int *nyv, int *ipbc) {
   cvgpush2lt(part,fxy,*qbm,*dt,ek,*idimp,*nop,*npe,*nx,*ny,*nxv,*nyv,
              *ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpost2lt_(float *part, float *q, float *qm, int *nop, int *npe,
                int *idimp, int *nxv, int *nyv) {
   cgpost2lt(part,q,*qm,*nop,*npe,*idimp,*nxv,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cvgpost2lt_(float *part, float *q, float *qm, int *nop, int *npe,
                 int *idimp, int *nxv, int *nyv) {
   cvgpost2lt(part,q,*qm,*nop,*npe,*idimp,*nxv,*nyv);
   return;
} 

/*--------------------------------------------------------------------*/
void cdsortp2ylt_(float *parta, float *partb, int *npic, int *idimp,
                  int *nop, int *npe, int *ny1) {
   cdsortp2ylt(parta,partb,npic,*idimp,*nop,*npe,*ny1);
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
void cvpois22_(float complex *q, float complex *fxy, int *isign,
               float complex *ffc, float *ax, float *ay, float *affp,
               float *we, int *nx, int *ny, int *nxvh, int *nyv,
               int *nxhd, int *nyhd) {
   cvpois22(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*nxvh,*nyv,*nxhd,
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
void cwfft2rvx_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nxhd,
                int *nyd, int *nxhyd, int *nxyhd) {
   cwfft2rvx(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rv2_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nxhd,
                int *nyd, int *nxhyd, int *nxyhd) {
   cwfft2rv2(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,*nxyhd);
   return;
}
