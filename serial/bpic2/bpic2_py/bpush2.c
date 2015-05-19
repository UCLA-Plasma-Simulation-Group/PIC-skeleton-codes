/* C Library for Skeleton 2-1/2D Electromagnetic PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include "bpush2.h"

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
void cgbpush23l(float part[], float fxy[], float bxy[], float qbm,
                float dt, float dtc, float *ek, int idimp, int nop,
                int nx, int ny, int nxv, int nyv, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field. Using the Boris Mover.
   scalar version using guard cells
   119 flops/particle, 1 divide, 29 loads, 5 stores
   input: all, output: part, ek
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
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = velocity vx of particle n
   part[n][3] = velocity vy of particle n
   part[n][4] = velocity vz of particle n
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
        (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
   idimp = size of phase space = 5
   nop = number of particles
   nx/ny = system length in x/y direction
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int j, nn, mm, np, mp, nxv3;
   float qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt;
   float anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   double sum1;
   nxv3 = 3*nxv;
   qtmh = 0.5*qbm*dt;
   sum1 = 0.0;
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
   for (j = 0; j < nop; j++) {
/* find interpolation weights */
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = part[idimp*j] - (float) nn;
      dyp = part[1+idimp*j] - (float) mm;
      nn = 3*nn;
      mm = nxv3*mm;
      amx = 1.0 - dxp;
      mp = mm + nxv3;
      amy = 1.0 - dyp;
      np = nn + 3;
/* find electric field */
      dx = dyp*(dxp*fxy[np+mp] + amx*fxy[nn+mp])
         + amy*(dxp*fxy[np+mm] + amx*fxy[nn+mm]);
      dy = dyp*(dxp*fxy[1+np+mp] + amx*fxy[1+nn+mp])
         + amy*(dxp*fxy[1+np+mm] + amx*fxy[1+nn+mm]);
      dz = dyp*(dxp*fxy[2+np+mp] + amx*fxy[2+nn+mp])
         + amy*(dxp*fxy[2+np+mm] + amx*fxy[2+nn+mm]);
/* find magnetic field */
      ox = dyp*(dxp*bxy[np+mp] + amx*bxy[nn+mp])
         + amy*(dxp*bxy[np+mm] + amx*bxy[nn+mm]);
      oy = dyp*(dxp*bxy[1+np+mp] + amx*bxy[1+nn+mp])
         + amy*(dxp*bxy[1+np+mm] + amx*bxy[1+nn+mm]);
      oz = dyp*(dxp*bxy[2+np+mp] + amx*bxy[2+nn+mp])
         + amy*(dxp*bxy[2+np+mm] + amx*bxy[2+nn+mm]);
/* calculate half impulse */
      dx *= qtmh;
      dy *= qtmh;
      dz *= qtmh;
/* half acceleration */
      acx = part[2+idimp*j] + dx;
      acy = part[3+idimp*j] + dy;
      acz = part[4+idimp*j] + dz;
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
      part[2+idimp*j] = dx;
      part[3+idimp*j] = dy;
      part[4+idimp*j] = dz;
/* new position */
      dx = part[idimp*j] + dx*dtc;
      dy = part[1+idimp*j] + dy*dtc;
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
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* set new position */
      part[idimp*j] = dx;
      part[1+idimp*j] = dy;
   }
/* normalize kinetic energy */
   *ek += 0.5*sum1;
   return;
}


/*--------------------------------------------------------------------*/
void cdgbpush23l(double part[], double fxy[], double bxy[], double qbm,
                 double dt, double dtc, double *ek, int idimp, int nop,
                 int nx, int ny, int nxv, int nyv, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field. Using the Boris Mover.
   scalar version using guard cells
   119 flops/particle, 1 divide, 29 loads, 5 stores
   input: all, output: part, ek
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
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = velocity vx of particle n
   part[n][3] = velocity vy of particle n
   part[n][4] = velocity vz of particle n
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
        (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
   idimp = size of phase space = 5
   nop = number of particles
   nx/ny = system length in x/y direction
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int j, nn, mm, np, mp, nxv3;
   double qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   double dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt;
   double anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   double sum1;
   nxv3 = 3*nxv;
   qtmh = 0.5*qbm*dt;
   sum1 = 0.0;
/* set boundary values */
   edgelx = 0.0;
   edgely = 0.0;
   edgerx = (double) nx;
   edgery = (double) ny;
   if (ipbc==2) {
      edgelx = 1.0;
      edgely = 1.0;
      edgerx = (double) (nx-1);
      edgery = (double) (ny-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0;
      edgerx = (double) (nx-1);
   }
   for (j = 0; j < nop; j++) {
/* find interpolation weights */
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = part[idimp*j] - (double) nn;
      dyp = part[1+idimp*j] - (double) mm;
      nn = 3*nn;
      mm = nxv3*mm;
      amx = 1.0 - dxp;
      mp = mm + nxv3;
      amy = 1.0 - dyp;
      np = nn + 3;
/* find electric field */
      dx = dyp*(dxp*fxy[np+mp] + amx*fxy[nn+mp])
         + amy*(dxp*fxy[np+mm] + amx*fxy[nn+mm]);
      dy = dyp*(dxp*fxy[1+np+mp] + amx*fxy[1+nn+mp])
         + amy*(dxp*fxy[1+np+mm] + amx*fxy[1+nn+mm]);
      dz = dyp*(dxp*fxy[2+np+mp] + amx*fxy[2+nn+mp])
         + amy*(dxp*fxy[2+np+mm] + amx*fxy[2+nn+mm]);
/* find magnetic field */
      ox = dyp*(dxp*bxy[np+mp] + amx*bxy[nn+mp])
         + amy*(dxp*bxy[np+mm] + amx*bxy[nn+mm]);
      oy = dyp*(dxp*bxy[1+np+mp] + amx*bxy[1+nn+mp])
         + amy*(dxp*bxy[1+np+mm] + amx*bxy[1+nn+mm]);
      oz = dyp*(dxp*bxy[2+np+mp] + amx*bxy[2+nn+mp])
         + amy*(dxp*bxy[2+np+mm] + amx*bxy[2+nn+mm]);
/* calculate half impulse */
      dx *= qtmh;
      dy *= qtmh;
      dz *= qtmh;
/* half acceleration */
      acx = part[2+idimp*j] + dx;
      acy = part[3+idimp*j] + dy;
      acz = part[4+idimp*j] + dz;
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
      part[2+idimp*j] = dx;
      part[3+idimp*j] = dy;
      part[4+idimp*j] = dz;
/* new position */
      dx = part[idimp*j] + dx*dtc;
      dy = part[1+idimp*j] + dy*dtc;
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
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* set new position */
      part[idimp*j] = dx;
      part[1+idimp*j] = dy;
   }
/* normalize kinetic energy */
   *ek += 0.5*sum1;
   return;
}

/*--------------------------------------------------------------------*/
void cgrbpush23l(float part[], float fxy[], float bxy[], float qbm,
                 float dt, float dtc, float ci, float *ek, int idimp,
                 int nop, int nx, int ny, int nxv, int nyv, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   scalar version using guard cells
   131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
   input: all, output: part, ek
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
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = momentum px of particle n
   part[n][3] = momentum py of particle n
   part[n][4] = momentum pz of particle n
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   idimp = size of phase space = 5
   nop = number of particles
   nx/ny = system length in x/y direction
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int j, nn, mm, np, mp, nxv3;
   float qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg;
   float omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   double sum1;
   nxv3 = 3*nxv;
   qtmh = 0.5*qbm*dt;
   ci2 = ci*ci;
   sum1 = 0.0;
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
   for (j = 0; j < nop; j++) {
/* find interpolation weights */
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = part[idimp*j] - (float) nn;
      dyp = part[1+idimp*j] - (float) mm;
      nn = 3*nn;
      mm = nxv3*mm;
      amx = 1.0 - dxp;
      mp = mm + nxv3;
      amy = 1.0 - dyp;
      np = nn + 3;
/* find electric field */
      dx = dyp*(dxp*fxy[np+mp] + amx*fxy[nn+mp])
         + amy*(dxp*fxy[np+mm] + amx*fxy[nn+mm]);
      dy = dyp*(dxp*fxy[1+np+mp] + amx*fxy[1+nn+mp])
         + amy*(dxp*fxy[1+np+mm] + amx*fxy[1+nn+mm]);
      dz = dyp*(dxp*fxy[2+np+mp] + amx*fxy[2+nn+mp])
         + amy*(dxp*fxy[2+np+mm] + amx*fxy[2+nn+mm]);
/* calculate half impulse */
      dx *= qtmh;
      dy *= qtmh;
      dz *= qtmh;
/* half acceleration */
      acx = part[2+idimp*j] + dx;
      acy = part[3+idimp*j] + dy;
      acz = part[4+idimp*j] + dz;
/* find inverse gamma */
      p2 = acx*acx + acy*acy + acz*acz;
      gami = 1.0/sqrtf(1.0 + p2*ci2);
/* find magnetic field */
      ox = dyp*(dxp*bxy[np+mp] + amx*bxy[nn+mp])
         + amy*(dxp*bxy[np+mm] + amx*bxy[nn+mm]);
      oy = dyp*(dxp*bxy[1+np+mp] + amx*bxy[1+nn+mp])
         + amy*(dxp*bxy[1+np+mm] + amx*bxy[1+nn+mm]);
      oz = dyp*(dxp*bxy[2+np+mp] + amx*bxy[2+nn+mp])
         + amy*(dxp*bxy[2+np+mm] + amx*bxy[2+nn+mm]);
/* renormalize magnetic field */
      qtmg = qtmh*gami;
/* time-centered kinetic energy */
      sum1 += gami*p2/(1.0 + gami);
/* calculate cyclotron frequency */
      omxt = qtmg*ox;
      omyt = qtmg*oy;
      omzt = qtmg*oz;
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
/* new momentum */
      dx += (rot1*acx + rot2*acy + rot3*acz)*anorm;
      dy += (rot4*acx + rot5*acy + rot6*acz)*anorm;
      dz += (rot7*acx + rot8*acy + rot9*acz)*anorm;
      part[2+idimp*j] = dx;
      part[3+idimp*j] = dy;
      part[4+idimp*j] = dz;
/* update inverse gamma */
      p2 = dx*dx + dy*dy + dz*dz;
      dtg = dtc/sqrtf(1.0 + p2*ci2);
/* new position */
      dx = part[idimp*j] + dx*dtg;
      dy = part[1+idimp*j] + dy*dtg;
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
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* set new position */
      part[idimp*j] = dx;
      part[1+idimp*j] = dy;
   }
/* normalize kinetic energy */
   *ek += sum1;
   return;
}

/*--------------------------------------------------------------------*/
void cdgrbpush23l(double part[], double fxy[], double bxy[], double qbm,
                  double dt, double dtc, double ci, double *ek, int idimp,
                  int nop, int nx, int ny, int nxv, int nyv, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   scalar version using guard cells
   131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
   input: all, output: part, ek
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
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = momentum px of particle n
   part[n][3] = momentum py of particle n
   part[n][4] = momentum pz of particle n
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   idimp = size of phase space = 5
   nop = number of particles
   nx/ny = system length in x/y direction
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int j, nn, mm, np, mp, nxv3;
   double qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   double dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg;
   double omxt, omyt, omzt, omt, anorm;
   double rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   double sum1;
   nxv3 = 3*nxv;
   qtmh = 0.5*qbm*dt;
   ci2 = ci*ci;
   sum1 = 0.0;
/* set boundary values */
   edgelx = 0.0;
   edgely = 0.0;
   edgerx = (double) nx;
   edgery = (double) ny;
   if (ipbc==2) {
      edgelx = 1.0;
      edgely = 1.0;
      edgerx = (double) (nx-1);
      edgery = (double) (ny-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0;
      edgerx = (double) (nx-1);
   }
   for (j = 0; j < nop; j++) {
/* find interpolation weights */
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = part[idimp*j] - (double) nn;
      dyp = part[1+idimp*j] - (double) mm;
      nn = 3*nn;
      mm = nxv3*mm;
      amx = 1.0 - dxp;
      mp = mm + nxv3;
      amy = 1.0 - dyp;
      np = nn + 3;
/* find electric field */
      dx = dyp*(dxp*fxy[np+mp] + amx*fxy[nn+mp])
         + amy*(dxp*fxy[np+mm] + amx*fxy[nn+mm]);
      dy = dyp*(dxp*fxy[1+np+mp] + amx*fxy[1+nn+mp])
         + amy*(dxp*fxy[1+np+mm] + amx*fxy[1+nn+mm]);
      dz = dyp*(dxp*fxy[2+np+mp] + amx*fxy[2+nn+mp])
         + amy*(dxp*fxy[2+np+mm] + amx*fxy[2+nn+mm]);
/* calculate half impulse */
      dx *= qtmh;
      dy *= qtmh;
      dz *= qtmh;
/* half acceleration */
      acx = part[2+idimp*j] + dx;
      acy = part[3+idimp*j] + dy;
      acz = part[4+idimp*j] + dz;
/* find inverse gamma */
      p2 = acx*acx + acy*acy + acz*acz;
      gami = 1.0/sqrt(1.0 + p2*ci2);
/* find magnetic field */
      ox = dyp*(dxp*bxy[np+mp] + amx*bxy[nn+mp])
         + amy*(dxp*bxy[np+mm] + amx*bxy[nn+mm]);
      oy = dyp*(dxp*bxy[1+np+mp] + amx*bxy[1+nn+mp])
         + amy*(dxp*bxy[1+np+mm] + amx*bxy[1+nn+mm]);
      oz = dyp*(dxp*bxy[2+np+mp] + amx*bxy[2+nn+mp])
         + amy*(dxp*bxy[2+np+mm] + amx*bxy[2+nn+mm]);
/* renormalize magnetic field */
      qtmg = qtmh*gami;
/* time-centered kinetic energy */
      sum1 += gami*p2/(1.0 + gami);
/* calculate cyclotron frequency */
      omxt = qtmg*ox;
      omyt = qtmg*oy;
      omzt = qtmg*oz;
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
/* new momentum */
      dx += (rot1*acx + rot2*acy + rot3*acz)*anorm;
      dy += (rot4*acx + rot5*acy + rot6*acz)*anorm;
      dz += (rot7*acx + rot8*acy + rot9*acz)*anorm;
      part[2+idimp*j] = dx;
      part[3+idimp*j] = dy;
      part[4+idimp*j] = dz;
/* update inverse gamma */
      p2 = dx*dx + dy*dy + dz*dz;
      dtg = dtc/sqrt(1.0 + p2*ci2);
/* new position */
      dx = part[idimp*j] + dx*dtg;
      dy = part[1+idimp*j] + dy*dtg;
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
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* set new position */
      part[idimp*j] = dx;
      part[1+idimp*j] = dy;
   }
/* normalize kinetic energy */
   *ek += sum1;
   return;
}

/*--------------------------------------------------------------------*/
void cgpost2l(float part[], float q[], float qm, int nop, int idimp,
              int nxv, int nyv) {
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
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   q[k][j] = charge density at grid point j,k
   qm = charge on particle, in units of e
   nop = number of particles
   idimp = size of phase space = 4
   nxv = first dimension of charge array, must be >= nx+1
   nyv = second dimension of charge array, must be >= ny+1
local data                                                            */
   int j, nn, mm, np, mp;
   float dxp, dyp, amx, amy;
/* find interpolation weights */
   for (j = 0; j < nop; j++) {
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = qm*(part[idimp*j] - (float) nn);
      dyp = part[1+idimp*j] - (float) mm;
      mm = nxv*mm;
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
void cgjpost2l(float part[], float cu[], float qm, float dt, int nop,
               int idimp, int nx, int ny, int nxv, int nyv, int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   scalar version using guard cells
   41 flops/particle, 17 loads, 14 stores
   input: all, output: part, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*vi, where i = x,y,z
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = x velocity of particle n
   part[n][3] = y velocity of particle n
   part[n][4] = z velocity of particle n
   cu[k][j][i] = ith component of current density at grid point j,k
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nop = number of particles
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   nxv = first dimension of current array, must be >= nx+1
   nyv = second dimension of current array, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int j, nn, mm, np, mp, nxv3;
   float edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, vx, vy, vz;
   nxv3 = 3*nxv;
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
   for (j = 0; j < nop; j++) {
/* find interpolation weights */
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = qm*(part[idimp*j] - (float) nn);
      dyp = part[1+idimp*j] - (float) mm;
      nn = 3*nn;
      mm = nxv3*mm;
      amx = qm - dxp;
      mp = mm + nxv3;
      amy = 1.0 - dyp;
      np = nn + 3;
/* deposit current */
      dx = dxp*dyp;
      dy = amx*dyp;
      vx = part[2+idimp*j];
      vy = part[3+idimp*j];
      vz = part[4+idimp*j];
      cu[np+mp] += vx*dx;
      cu[1+np+mp] += vy*dx;
      cu[2+np+mp] += vz*dx;
      dx = dxp*amy;
      cu[nn+mp] += vx*dy;
      cu[1+nn+mp] += vy*dy;
      cu[2+nn+mp] += vz*dy;
      dy = amx*amy;
      cu[np+mm] += vx*dx;
      cu[1+np+mm] += vy*dx;
      cu[2+np+mm] += vz*dx;
      cu[nn+mm] += vx*dy;
      cu[1+nn+mm] += vy*dy;
      cu[2+nn+mm] += vz*dy;
/* advance position half a time-step */
      dx = part[idimp*j] + vx*dt;
      dy = part[1+idimp*j] + vy*dt;
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
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* set new position */
      part[idimp*j] = dx;
      part[1+idimp*j] = dy;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cgrjpost2l(float part[], float cu[], float qm, float dt, float ci,
                int nop, int idimp, int nx, int ny, int nxv, int nyv,
                int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation for relativistic particles
   in addition, particle positions are advanced a half time-step
   scalar version using guard cells
   47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
   input: all, output: part, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*pi*gami, where i = x,y,z
   where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
   part[n][0] = position x of particle n
   part[n][1] = position y of particle n
   part[n][2] = x momentum of particle n
   part[n][3] = y momentum of particle n
   part[n][4] = z momentum of particle n
   cu[k][j][i] = ith component of current density at grid point j,k
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
   nop = number of particles
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   nxv = first dimension of current array, must be >= nx+1
   nyv = second dimension of current array, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int j, nn, mm, np, mp, nxv3;
   float ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, vx, vy, vz, p2, gami;
   nxv3 = 3*nxv;
   ci2 = ci*ci;
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
   for (j = 0; j < nop; j++) {
/* find interpolation weights */
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = qm*(part[idimp*j] - (float) nn);
      dyp = part[1+idimp*j] - (float) mm;
/* find inverse gamma */
      vx = part[2+idimp*j];
      vy = part[3+idimp*j];
      vz = part[4+idimp*j];
      p2 = vx*vx + vy*vy + vz*vz;
      gami = 1.0/sqrtf(1.0 + p2*ci2);
/* calculate weights */
      nn = 3*nn;
      mm = nxv3*mm;
      amx = qm - dxp;
      mp = mm + nxv3;
      amy = 1.0 - dyp;
      np = nn + 3;
/* deposit current */
      dx = dxp*dyp;
      dy = amx*dyp;
      vx *= gami;
      vy *= gami;
      vz *= gami;
      cu[np+mp] += vx*dx;
      cu[1+np+mp] += vy*dx;
      cu[2+np+mp] += vz*dx;
      dx = dxp*amy;
      cu[nn+mp] += vx*dy;
      cu[1+nn+mp] += vy*dy;
      cu[2+nn+mp] += vz*dy;
      dy = amx*amy;
      cu[np+mm] += vx*dx;
      cu[1+np+mm] += vy*dx;
      cu[2+np+mm] += vz*dx;
      cu[nn+mm] += vx*dy;
      cu[1+nn+mm] += vy*dy;
      cu[2+nn+mm] += vz*dy;
/* advance position half a time-step */
      dx = part[idimp*j] + vx*dt;
      dy = part[1+idimp*j] + vy*dt;
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
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* set new position */
      part[idimp*j] = dx;
      part[1+idimp*j] = dy;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cdsortp2yl(float parta[], float partb[], int npic[], int idimp,
                int nop, int ny1) {
/* this subroutine sorts particles by y grid
   linear interpolation
   parta/partb = input/output particle arrays
   parta[n][1] = position y of particle n
   npic = address offset for reordering particles
   idimp = size of phase space = 4
   nop = number of particles
   ny1 = system length in y direction + 1
      dimension parta(idimp,nop), partb(idimp,nop), npic(ny1)
local data                                                            */
   int i, j, k, m, isum, ist, ip;
/* clear counter array */
   for (k = 0; k < ny1; k++) {
      npic[k] = 0;
   }
/* find how many particles in each grid */
   for (j = 0; j < nop; j++) {
      m = parta[1+idimp*j];
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
      m = parta[1+idimp*j];
      ip = npic[m];
      for (i = 0; i < idimp; i++) {
         partb[i+idimp*ip] = parta[i+idimp*j];
      }
      npic[m] = ip + 1;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cbguard2l(float bxy[], int nx, int ny, int nxe, int nye) {
/* replicate extended periodic vector field bxy
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
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
   nye = second dimension of field arrays, must be >= ny+1
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
void cpois23(float complex q[], float complex fxy[], int isign,
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
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
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
   }
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
   *we = wp*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void ccuperp2(float complex cu[], int nx, int ny, int nxvh, int nyv) {
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
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
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
void cibpois23(float complex cu[], float complex bxy[],
               float complex ffc[], float ci, float *wm, int nx, int ny,
               int nxvh, int nyv, int nxhd, int nyhd) {
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
   nxvh = first dimension of field arrays, must be >= nxh
   nyv = second dimension of field arrays, must be >= ny
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, j, k, k1, kk, kj;
   float dnx, dny, dky, ci2, at1, at2, at3;
   float complex zero, zt1, zt2, zt3;
   double wp;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
   ci2 = ci*ci;
/* calculate magnetic field and sum field energy */
   wp = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      for (j = 1; j < nxh; j++) {
         at1 = ci2*crealf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         at1 = at1*cimagf(ffc[j+kk]);
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
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      at1 = ci2*crealf(ffc[kk]);
      at3 = at1*dny*(float) k;
      at1 = at1*cimagf(ffc[kk]);
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
   }
/* mode numbers ky = 0, ny/2 */
   k1 = 3*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      at1 = ci2*crealf(ffc[j]);
      at2 = at1*dnx*(float) j; 
      at1 = at1*cimagf(ffc[j]);
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
   *wm = wp*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void cmaxwel2(float complex exy[], float complex bxy[],
              float complex cu[], float complex ffc[], float ci,
              float dt, float *wf, float *wm, int nx, int ny, int nxvh,
              int nyv, int nxhd, int nyhd) {
/* this subroutine solves 2-1/2d maxwell's equation in fourier space for
   transverse electric and magnetic fields with periodic boundary
   conditions.
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
   nxvh = first dimension of field arrays, must be >= nxh
   nyv = second dimension of field arrays, must be >= ny
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, j, k, k1, kk, kj;
   float dnx, dny, dth, c2, cdt, affp, anorm, dkx, dky, afdt, adt;
   float complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9;
   double wp, ws;
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
   ws = 0.0;
   wp = 0.0;
/* calculate the electromagnetic fields */
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      for (j = 1; j < nxh; j++) {
         dkx = dnx*(float) j;
         afdt = adt*cimagf(ffc[j+kk]);
/* update magnetic field half time step, ky > 0 */
         zt1 = -cimagf(exy[2+3*j+3*kj])
               + crealf(exy[2+3*j+3*kj])*_Complex_I;
         zt2 = -cimagf(exy[1+3*j+3*kj])
             + crealf(exy[1+3*j+3*kj])*_Complex_I;
         zt3 = -cimagf(exy[3*j+3*kj]) + crealf(exy[3*j+3*kj])*_Complex_I;
         zt4 = bxy[3*j+3*kj] - dth*(dky*zt1);
         zt5 = bxy[1+3*j+3*kj] + dth*(dkx*zt1);
         zt6 = bxy[2+3*j+3*kj] - dth*(dkx*zt2 - dky*zt3);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exy[3*j+3*kj] + cdt*(dky*zt1) - afdt*cu[3*j+3*kj];
         zt8 = exy[1+3*j+3*kj] - cdt*(dkx*zt1) - afdt*cu[1+3*j+3*kj];
         zt9 = exy[2+3*j+3*kj] + cdt*(dkx*zt2 - dky*zt3)
               - afdt*cu[2+3*j+3*kj];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exy[3*j+3*kj] = zt7;
         exy[1+3*j+3*kj] = zt8;
         exy[2+3*j+3*kj] = zt9;
         ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         zt4 -= dth*(dky*zt1);
         zt5 += dth*(dkx*zt1);
         zt6 -= dth*(dkx*zt2 - dky*zt3);
         bxy[3*j+3*kj] = zt4;
         bxy[1+3*j+3*kj] = zt5;
         bxy[2+3*j+3*kj] = zt6;
         wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
/* update magnetic field half time step, ky < 0 */
         zt1 = -cimagf(exy[2+3*j+3*k1])
               + crealf(exy[2+3*j+3*k1])*_Complex_I;
         zt2 = -cimagf(exy[1+3*j+3*k1])
               + crealf(exy[1+3*j+3*k1])*_Complex_I;
         zt3 = -cimagf(exy[3*j+3*k1]) + crealf(exy[3*j+3*k1])*_Complex_I;
         zt4 = bxy[3*j+3*k1] + dth*(dky*zt1);
         zt5 = bxy[1+3*j+3*k1] + dth*(dkx*zt1);
         zt6 = bxy[2+3*j+3*k1] - dth*(dkx*zt2 + dky*zt3);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exy[3*j+3*k1] - cdt*(dky*zt1) - afdt*cu[3*j+3*k1];
         zt8 = exy[1+3*j+3*k1] - cdt*(dkx*zt1) - afdt*cu[1+3*j+3*k1];
         zt9 = exy[2+3*j+3*k1] + cdt*(dkx*zt2 + dky*zt3)
               - afdt*cu[2+3*j+3*k1];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exy[3*j+3*k1] = zt7;
         exy[1+3*j+3*k1] = zt8;
         exy[2+3*j+3*k1] = zt9;
         ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         zt4 += dth*(dky*zt1);
         zt5 += dth*(dkx*zt1);
         zt6 -= dth*(dkx*zt2 + dky*zt3);
         bxy[3*j+3*k1] = zt4;
         bxy[1+3*j+3*k1] = zt5;
         bxy[2+3*j+3*k1] = zt6;
         wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
      }
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      dky = dny*(float) k;
      afdt = adt*cimagf(ffc[kk]);
/* update magnetic field half time step */
      zt1 = -cimagf(exy[2+3*kj]) + crealf(exy[2+3*kj])*_Complex_I;
      zt3 = -cimagf(exy[3*kj]) + crealf(exy[3*kj])*_Complex_I;
      zt4 = bxy[3*kj] - dth*(dky*zt1);
      zt6 = bxy[2+3*kj] + dth*(dky*zt3);
/* update electric field whole time step */
      zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
      zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
      zt7 = exy[3*kj] + cdt*(dky*zt1) - afdt*cu[3*kj];
      zt9 = exy[2+3*kj] - cdt*(dky*zt3) - afdt*cu[2+3*kj];
/* update magnetic field half time step and store electric field */
      zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
      zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
      exy[3*kj] = zt7;
      exy[1+3*kj] = zero;
      exy[2+3*kj] = zt9;
      ws += anorm*(zt7*conjf(zt7) + zt9*conjf(zt9));
      zt4 -= dth*(dky*zt1);
      zt6 += dth*(dky*zt3);
      bxy[3*kj] = zt4;
      bxy[1+3*kj] = zero;
      bxy[2+3*kj] = zt6;
      wp += anorm*(zt4*conjf(zt4) + zt6*conjf(zt6));
      bxy[3*k1] = zero;
      bxy[1+3*k1] = zero;
      bxy[2+3*k1] = zero;
      exy[3*k1] = zero;
      exy[1+3*k1] = zero;
      exy[2+3*k1] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = 3*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      dkx = dnx*(float) j; 
      afdt = adt*cimagf(ffc[j]);
/* update magnetic field half time step */
      zt1 = -cimagf(exy[2+3*j]) + crealf(exy[2+3*j])*_Complex_I;
      zt2 = -cimagf(exy[1+3*j]) + crealf(exy[1+3*j])*_Complex_I;
      zt5 = bxy[1+3*j] + dth*(dkx*zt1);
      zt6 = bxy[2+3*j] - dth*(dkx*zt2);
/* update electric field whole time step */
      zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
      zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
      zt8 = exy[1+3*j] - cdt*(dkx*zt1) - afdt*cu[1+3*j];
      zt9 = exy[2+3*j] + cdt*(dkx*zt2) - afdt*cu[2+3*j];
/* update magnetic field half time step and store electric field */
      zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
      zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
      exy[3*j] = zero;
      exy[1+3*j] = zt8;
      exy[2+3*j] = zt9;
      ws += anorm*(zt8*conjf(zt8) + zt9*conjf(zt9));
      zt5 += dth*(dkx*zt1);
      zt6 -= dth*(dkx*zt2);
      bxy[3*j] = zero;
      bxy[1+3*j] = zt5;
      bxy[2+3*j] = zt6;
      wp += anorm*(zt5*conjf(zt5) + zt6*conjf(zt6));
      bxy[3*j+k1] = zero;
      bxy[1+3*j+k1] = zero;
      bxy[2+3*j+k1] = zero;
      exy[3*j+k1] = zero;
      exy[1+3*j+k1] = zero;
      exy[2+3*j+k1] = zero;
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
   *wf = ws*(float) (nx*ny);
   *wm = c2*wp*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void cemfield2(float complex fxy[], float complex exy[],
               float complex ffc[], int isign, int nx, int ny, int nxvh,
               int nyv, int nxhd, int nyhd) {
/* this subroutine either adds complex vector fields if isign > 0
   or copies complex vector fields if isign < 0
   includes additional smoothing
local data                                                 */
   int i, j, k, nxh, nyh, k1, kk, kj;
   float at1;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
/* add the fields */
   if (isign > 0) {
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
         for (j = 0; j < nxh; j++) {
            at1 = cimagf(ffc[j+kk]);
            for (i = 0; i < 3; i++) {
              fxy[i+3*j+3*kj] += exy[i+3*j+3*kj]*at1;
              fxy[i+3*j+3*k1] += exy[i+3*j+3*k1]*at1;
            }
         }
      }
      k1 = 3*nxvh*nyh;
      for (j = 0; j < nxh; j++) {
         at1 = cimagf(ffc[j]);
         for (i = 0; i < 3; i++) {
            fxy[i+3*j] += exy[i+3*j]*at1;
            fxy[i+3*j+k1] += exy[i+3*j+k1]*at1;
         }
      }
   }
/* copy the fields */
   else if (isign < 0) {
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
         for (j = 0; j < nxh; j++) {
            at1 = cimagf(ffc[j+kk]);
            for (i = 0; i < 3; i++) {
               fxy[i+3*j+3*kj] = exy[i+3*j+3*kj]*at1;
               fxy[i+3*j+3*k1] = exy[i+3*j+3*k1]*at1;
            }
         }
      }
      k1 = 3*nxvh*nyh;
      for (j = 0; j < nxh; j++) {
         at1 = cimagf(ffc[j]);
         for (i = 0; i < 3; i++) {
            fxy[i+3*j] = exy[i+3*j]*at1;
            fxy[i+3*j+k1] = exy[i+3*j+k1]*at1;
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
void cfft2rxx(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nyi, int nyp,
              int nxhd, int nyd, int nxhyd, int nxyhd) {
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
   f[k][j] = mode j,k, where 0 <= j < nx/2 and 0 <= k < ny,
   except for f[k][0] =  mode nx/2,k-1, where ny/2+1 <= k < ny, and
   imag(f[0][0]) = real part of mode nx/2,0 and
   imag(f[ny/2][0]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
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
         for (j = 0; j < ns; j++) {
            j1 = j + k1;
            j2 = j + k2;
            t1 = sct[kmr*j];
            for (i = nyi-1; i < nyt; i++) {
               joff = nxhd*i;
               t2 = t1*f[j2+joff];
               f[j2+joff] = f[j1+joff] - t2;
               f[j1+joff] += t2;
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
      for (k = nyi-1; k < nyt; k++) {
         joff = nxhd*k;
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
   for (j = 1; j < nxhh; j++) {
      t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
      for (k = nyi-1; k < nyt; k++) {
         joff = nxhd*k;
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
         for (j = 0; j < ns; j++) {
            j1 = j + k1;
            j2 = j + k2;
            t1 = conjf(sct[kmr*j]);
            for (i = nyi-1; i < nyt; i++) {
               joff = nxhd*i;
               t2 = t1*f[j2+joff];
               f[j2+joff] = f[j1+joff] - t2;
               f[j1+joff] += t2;
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
   f[k][j] = mode j,k, where 0 <= j < nx/2 and 0 <= k < ny,
   except for f[k][0] =  mode nx/2,k-1, where ny/2+1 <= k < ny, and
   imag(f[0][0]) = real part of mode nx/2,0 and
   imag(f[ny/2][0]) = real part of mode nx/2,ny/2
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
void cfft2r3x(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nyi, int nyp,
              int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the x part of 3 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   y, using complex arithmetic
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
   nxhd = second dimension of f
   nyd = third dimension of f
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j][0:2] = mode j,k, where 0 <= j < nx/2 and 0 <= k < ny,
   except for f[k][0][0:2] =  mode nx/2,k-1, where ny/2+1 <= k < ny, and
   imag(f[0][0][0:2]) = real part of mode nx/2,0 and
   imag(f[ny/2][0][0:2]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
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
      goto L140;
/* inverse fourier transform */
/* swap complex components */
   for (k = nyi-1; k < nyt; k++) {
      for (j = 0; j < nxh; j++) {
         joff = 3*nxhd*k;
         at1 = crealf(f[2+3*j+joff]);
         f[2+3*j+joff] = crealf(f[1+3*j+joff])
                         + cimagf(f[2+3*j+joff])*_Complex_I;
         at2 = cimagf(f[1+3*j+joff]);
         f[1+3*j+joff] = cimagf(f[3*j+joff]) + at1*_Complex_I;
         f[3*j+joff] = crealf(f[3*j+joff]) + at2*_Complex_I;
       }
   }
/* bit-reverse array elements in x */
   nrx = nxhy/nxh;
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = nyi-1; k < nyt; k++) {
         joff = 3*nxhd*k;
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
         for (j = 0; j < ns; j++) {
            j1 = j + k1;
            j2 = j + k2;
            t1 = sct[kmr*j];
            for (i = nyi-1; i < nyt; i++) {
               joff = 3*nxhd*i;
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
      }
      ns = ns2;
   }
/* unscramble coefficients and normalize */
   kmr = nxy/nx;
   ani = 1.0/(float) (2*nx*ny);
   for (j = 1; j < nxhh; j++) {
      t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
      for (k = nyi-1; k < nyt; k++) {
         joff = 3*nxhd*k;
         for (jj = 0; jj < 3; jj++) {
            t2 = conjf(f[jj+3*(nxh-j)+joff]);
            t1 = f[jj+3*j+joff] + t2;
            t2 = (f[jj+3*j+joff] - t2)*t3;
            f[jj+3*j+joff] = ani*(t1 + t2);
            f[jj+3*(nxh-j)+joff] = ani*conjf(t1 - t2);
         }
      }
   }
   ani = 2.0*ani;
   for (k = nyi-1; k < nyt; k++) {
      joff = 3*nxhd*k;
      for (jj = 0; jj < 3; jj++) {
         f[jj+3*nxhh+joff] = ani*conjf(f[jj+3*nxhh+joff]);
         f[jj+joff] = ani*((crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
/* scramble coefficients */
L140: kmr = nxy/nx;
   for (j = 1; j < nxhh; j++) {
      t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
      for (k = nyi-1; k < nyt; k++) {
         joff = 3*nxhd*k;
         for (jj = 0; jj < 3; jj++) {
            t2 = conjf(f[jj+3*(nxh-j)+joff]);
            t1 = f[jj+3*j+joff] + t2;
            t2 = (f[jj+3*j+joff] - t2)*t3;
            f[jj+3*j+joff] = t1 + t2;
            f[jj+3*(nxh-j)+joff] = conjf(t1 - t2);
         }
      }
   }
   for (k = nyi-1; k < nyt; k++) {
      joff = 3*nxhd*k;
      for (jj = 0; jj < 3; jj++) {
         f[jj+3*nxhh+joff] = 2.0*conjf(f[jj+3*nxhh+joff]);
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
         joff = 3*nxhd*k;
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
   nrx = nxy/nxh;
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
            for (i = nyi-1; i < nyt; i++) {
               joff = 3*nxhd*i;
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
      }
      ns = ns2;
   }
/* swap complex components */
   for (k = nyi-1; k < nyt; k++) {
      joff = 3*nxhd*k;
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
void cfft2r3y(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the y part of 3 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   x, using complex arithmetic
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
   nxhd = second dimension of f
   nyd = third dimension of f
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j][0:2] = mode j,k, where 0 <= j < nx/2 and 0 <= k < ny,
   except for f[k][0][0:2] =  mode nx/2,k-1, where ny/2+1 <= k < ny, and
   imag(f[0][0][0:2]) = real part of mode nx/2,0 and
   imag(f[ny/2][0][0:2]) = real part of mode nx/2,ny/2
  written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
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
      goto L90;
/* inverse fourier transform */
   nry = nxhy/ny;
/* bit-reverse array elements in y */
   for (k = 0; k < ny; k++) {
      joff = 3*nxhd*k;
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      k1 = 3*nxhd*k1;
      for (j = nxi-1; j < nxt; j++) {
         t1 = f[3*j+k1];
         t2 = f[1+3*j+k1];
         t3 = f[2+3*j+k1];
         f[3*j+k1] = f[3*j+joff];
         f[1+3*j+k1] = f[1+3*j+joff];
         f[2+3*j+k1] = f[2+3*j+joff];
         f[3*j+joff] = t1;
         f[1+3*j+joff] = t2;
         f[2+3*j+joff] = t3;
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
            j1 = 3*nxhd*(j + k1);
            j2 = 3*nxhd*(j + k2);
            t1 = sct[kmr*j];
            for (i = nxi-1; i < nxt; i++) {
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
      }
      ns = ns2;
   }
/* unscramble modes kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      if (nxi==1) {
         joff = 3*nxhd*k;
         k1 = 3*nxhd*ny - joff;
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
/* scramble modes kx = 0, nx/2 */
L90: for (k = 1; k < nyh; k++) {
      if (nxi==1) {
         joff = 3*nxhd*k;
         k1 = 3*nxhd*ny - joff;
         for (jj = 0; jj < 3; jj++) {
            t1 = cimagf(f[jj+k1]) + crealf(f[jj+k1])*_Complex_I;
            f[jj+k1] = conjf(f[jj+joff] - t1);
            f[jj+joff] += t1;
         }
      }
   }
/* bit-reverse array elements in y */
   nry = nxhy/ny;
   for (k = 0; k < ny; k++) {
      joff = 3*nxhd*k;
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      k1 = 3*nxhd*k1;
      for (j = nxi-1; j < nxt; j++) {
         t1 = f[3*j+k1];
         t2 = f[1+3*j+k1];
         t3 = f[2+3*j+k1];
         f[3*j+k1] = f[3*j+joff];
         f[1+3*j+k1] = f[1+3*j+joff];
         f[2+3*j+k1] = f[2+3*j+joff];
         f[3*j+joff] = t1;
         f[1+3*j+joff] = t2;
         f[2+3*j+joff] = t3;
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
            j1 = 3*nxhd*(j + k1);
            j2 = 3*nxhd*(j + k2);
            t1 = conjf(sct[kmr*j]);
            for (i = nxi-1; i < nxt; i++) {
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
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rx(float complex f[], int isign, int mixup[],
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
      cfft2rxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,nxyhd);
/* perform y fft */
      cfft2rxy(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfft2rxy(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,nxyhd);
/* perform x fft */
      cfft2rxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,nxyhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2r3(float complex f[],int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxhd, int nyd,
              int nxhyd, int nxyhd) {
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
      cfft2r3x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,nxyhd);
/* perform y fft */
      cfft2r3y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfft2r3y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,nxyhd);
/* perform x fft */
      cfft2r3x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,nxyhd);
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
void cdgbpush23l_(double *part, double *fxy, double *bxy, double *qbm,
                  double *dt, double *dtc, double *ek, int *idimp,
                  int *nop, int *nx, int *ny, int *nxv, int *nyv,
                  int *ipbc) {
   cdgbpush23l(part,fxy,bxy,*qbm,*dt,*dtc,ek,*idimp,*nop,*nx,*ny,*nxv,
               *nyv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbpush23l_(float *part, float *fxy, float *bxy, float *qbm,
                 float *dt, float *dtc, float *ek, int *idimp, int *nop,
                 int *nx, int *ny, int *nxv, int *nyv, int *ipbc) {
   cgbpush23l(part,fxy,bxy,*qbm,*dt,*dtc,ek,*idimp,*nop,*nx,*ny,*nxv,
              *nyv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdgrbpush23l_(double *part, double *fxy, double *bxy, double *qbm,
                   double *dt, double *dtc, double *ci, double *ek,
                   int *idimp, int *nop, int *nx, int *ny, int *nxv,
                   int *nyv, int *ipbc) {
   cdgrbpush23l(part,fxy,bxy,*qbm,*dt,*dtc,*ci,ek,*idimp,*nop,*nx,*ny,
                *nxv,*nyv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbpush23l_(float *part, float *fxy, float *bxy, float *qbm,
                  float *dt, float *dtc, float *ci, float *ek,
                  int *idimp, int *nop, int *nx, int *ny, int *nxv,
                  int *nyv, int *ipbc) {
   cgrbpush23l(part,fxy,bxy,*qbm,*dt,*dtc,*ci,ek,*idimp,*nop,*nx,*ny,
               *nxv,*nyv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpost2l_(float *part, float *q, float *qm, int *nop, int *idimp,
               int *nxv, int *nyv) {
   cgpost2l(part,q,*qm,*nop,*idimp,*nxv,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cgjpost2l_(float *part, float *cu, float *qm, float *dt, int *nop,
                int *idimp, int *nx, int *ny, int *nxv, int *nyv,
                int *ipbc) {
   cgjpost2l(part,cu,*qm,*dt,*nop,*idimp,*nx,*ny,*nxv,*nyv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjpost2l_(float *part, float *cu, float *qm, float *dt, float *ci,
                 int *nop, int *idimp, int *nx, int *ny, int *nxv,
                 int *nyv, int *ipbc) {
   cgrjpost2l(part,cu,*qm,*dt,*ci,*nop,*idimp,*nx,*ny,*nxv,*nyv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdsortp2yl_(float *parta, float *partb, int *npic, int *idimp,
                 int *nop, int *ny1) {
   cdsortp2yl(parta,partb,npic,*idimp,*nop,*ny1);
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
void cpois23_(float complex *q, float complex *fxy, int *isign,
              float complex *ffc, float *ax, float *ay, float *affp,
              float *we, int *nx, int *ny, int *nxvh, int *nyv,
              int *nxhd, int *nyhd) {
   cpois23(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*nxvh,*nyv,*nxhd,
           *nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void ccuperp2_(float complex *cu, int *nx, int *ny, int *nxvh, int *nyv) {
   ccuperp2(cu,*nx,*ny,*nxvh,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cibpois23_(float complex *cu, float complex *bxy,
                float complex *ffc, float *ci, float *wm, int *nx,
                int *ny, int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   cibpois23(cu,bxy,ffc,*ci,wm,*nx,*ny,*nxvh,*nyv,*nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmaxwel2_(float complex *exy, float complex *bxy,
               float complex *cu, float complex *ffc, float *ci,
               float *dt, float *wf, float *wm, int *nx, int *ny,
               int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   cmaxwel2(exy,bxy,cu,ffc,*ci,*dt,wf,wm,*nx,*ny,*nxvh,*nyv,*nxhd,
            *nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cemfield2_(float complex *fxy, float complex *exy,
                float complex *ffc, int *isign, int *nx, int *ny,
                int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   cemfield2(fxy,exy,ffc,*isign,*nx,*ny,*nxvh,*nyv,*nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                  int *nxhyd, int *nxyhd) {
   cwfft2rinit(mixup,sct,*indx,*indy,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rx_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxhd,
               int *nyd, int *nxhyd, int *nxyhd) {
   cwfft2rx(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2r3_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxhd,
               int *nyd, int *nxhyd, int *nxyhd) {
   cwfft2r3(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,*nxyhd);
   return;
}
