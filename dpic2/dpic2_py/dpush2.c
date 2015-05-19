/* C Library for Skeleton 2-1/2D Darwin PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include "dpush2.h"

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
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
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
void cgmjpost2l(float part[], float amu[], float qm, int nop, int idimp,
                int nxv, int nyv) {
/* for 2-1/2d code, this subroutine calculates particle momentum flux
   using first-order spline interpolation
   scalar version using guard cells
   51 flops/particle, 21 loads, 16 stores
   input: all, output: amu
   momentum flux is approximated by values at the nearest grid points
   amu(i,n,m)=qci*(1.-dx)*(1.-dy)
   amu(i,n+1,m)=qci*dx*(1.-dy)
   amu(i,n,m+1)=qci*(1.-dx)*dy
   amu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
   where vj = vj(t-dt/2) and vk = vk(t-dt/2)
   part[n][0] = position x of particle n at t
   part[n][1] = position y of particle n at t
   part[n][2] = x velocity of particle n at t - dt/2
   part[n][3] = y velocity of particle n at t - dt/2
   part[n][4] = z velocity of particle n at t - dt/2
   amu[k][j][i] = ith component of momentum flux at grid point j,k
   qm = charge on particle, in units of e
   nop = number of particles
   idimp = size of phase space = 5
   nxv = second dimension of flux array, must be >= nx+1
   nyv = third dimension of flux array, must be >= ny+1
local data                                                            */
   int j, nn, mm, np, mp, nxv4;
   float dxp, dyp, amx, amy;
   float dx, dy, vx, vy, vz, v1, v2, v3, v4;
   nxv4 = 4*nxv;
   for (j = 0; j < nop; j++) {
/* find interpolation weights */
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = qm*(part[idimp*j] - (float) nn);
      dyp = part[1+idimp*j] - (float) mm;
      nn = 4*nn;
      mm = nxv4*mm;
      amx = qm - dxp;
      mp = mm + nxv4;
      amy = 1.0 - dyp;
      np = nn + 4;
/* deposit momentum flux */
      dx = dxp*dyp;
      dy = amx*dyp;
      vx = part[2+idimp*j];
      vy = part[3+idimp*j];
      vz = part[4+idimp*j];
      v1 = vx*vx - vy*vy;
      v2 = vx*vy;
      v3 = vz*vx;
      v4 = vz*vy;
      amu[np+mp] += v1*dx;
      amu[1+np+mp] += v2*dx;
      amu[2+np+mp] += v3*dx;
      amu[3+np+mp] += v4*dx;
      dx = dxp*amy;
      amu[nn+mp] += v1*dy;
      amu[1+nn+mp] += v2*dy;
      amu[2+nn+mp] += v3*dy;
      amu[3+nn+mp] += v4*dy;
      dy = amx*amy;
      amu[np+mm] += v1*dx;
      amu[1+np+mm] += v2*dx;
      amu[2+np+mm] += v3*dx;
      amu[3+np+mm] += v4*dx;
      amu[nn+mm] += v1*dy;
      amu[1+nn+mm] += v2*dy;
      amu[2+nn+mm] += v3*dy;
      amu[3+nn+mm] += v4*dy;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cgdjpost2l(float part[], float fxy[], float bxy[], float dcu[],
                float amu[], float qm, float qbm, float dt, int idimp,
                int nop, int nxv, int nyv) {
/* for 2-1/2d code, this subroutine calculates particle momentum flux,
   and acceleration density using first-order spline interpolation.
   scalar version using guard cells
   194 flops/particle, 1 divide, 67 loads, 28 stores
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
      rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
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
   dcu[k][j][i] = ith component of acceleration density
   at grid point j,k for i = 0, 2
   amu[k][j][i] = ith component of momentum flux
   at grid point j,k for i = 0, 3
   qm = charge on particle, in units of e
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   idimp = size of phase space = 5
   nop = number of particles
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
local data                                                            */
   int j, nn, mm, np, mp, nnn, mmm, nnp, mmp, nxv3, nxv4;
   float qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float vx, vy, vz, v1, v2, v3, v4;
   nxv3 = 3*nxv;
   nxv4 = 4*nxv;
   qtmh = 0.5*qbm*dt;
   dti = 1.0/dt;
   for (j = 0; j < nop; j++) {
/* find interpolation weights */
      nnn = part[idimp*j];
      mmm = part[1+idimp*j];
      dxp = part[idimp*j] - (float) nnn;
      dyp = part[1+idimp*j] - (float) mmm;
      nn = 3*nnn;
      mm = nxv3*mmm;
      amx = 1.0 - dxp;
      mp = mm + nxv3;
      amy = 1.0 - dyp;
      np = nn + 3;
      nnn = 4*nnn;
      mmm = nxv4*mmm;
      mmp = mmm + nxv4;
      nnp = nnn + 4;
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
      vx = part[2+idimp*j];
      vy = part[3+idimp*j];
      vz = part[4+idimp*j];
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
      dx = dxp*dyp;
      dy = amx*dyp;
      v1 = ox*ox - oy*oy;
      v2 = ox*oy;
      v3 = oz*ox;
      v4 = oz*oy;
      amu[nnp+mmp] += v1*dx;
      amu[1+nnp+mmp] += v2*dx;
      amu[2+nnp+mmp] += v3*dx;
      amu[3+nnp+mmp] += v4*dx;
      dcu[np+mp] += vx*dx;
      dcu[1+np+mp] += vy*dx;
      dcu[2+np+mp] += vz*dx;
      dx = dxp*amy;
      amu[nnn+mmp] += v1*dy;
      amu[1+nnn+mmp] += v2*dy;
      amu[2+nnn+mmp] += v3*dy;
      amu[3+nnn+mmp] += v4*dy;
      dcu[nn+mp] += vx*dy;
      dcu[1+nn+mp] += vy*dy;
      dcu[2+nn+mp] += vz*dy;
      dy = amx*amy;
      amu[nnp+mmm] += v1*dx;
      amu[1+nnp+mmm] += v2*dx;
      amu[2+nnp+mmm] += v3*dx;
      amu[3+nnp+mmm] += v4*dx;
      dcu[np+mm] += vx*dx;
      dcu[1+np+mm] += vy*dx;
      dcu[2+np+mm] += vz*dx;
      amu[nnn+mmm] += v1*dy;
      amu[1+nnn+mmm] += v2*dy;
      amu[2+nnn+mmm] += v3*dy;
      amu[3+nnn+mmm] += v4*dy;
      dcu[nn+mm] += vx*dy;
      dcu[1+nn+mm] += vy*dy;
      dcu[2+nn+mm] += vz*dy;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cgdcjpost2l(float part[], float fxy[], float bxy[], float cu[],
                 float dcu[], float amu[], float qm, float qbm,
                 float dt, int idimp, int nop, int nxv, int nyv) {
/* for 2-1/2d code, this subroutine calculates particle momentum flux,
   acceleration density and current density using first-order spline
   interpolation.
   scalar version using guard cells
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
      rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
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
   part[n][0] = position x of particle n at t
   part[n][1] = position y of particle n at t
   part[n][2] = velocity vx of particle n at t - dt/2
   part[n][3] = velocity vy of particle n at t - dt/2
   part[n][4] = velocity vz of particle n at t - dt/2
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
   qm = charge on particle, in units of e
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   idimp = size of phase space = 5
   nop = number of particles
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
local data                                                            */
   int j, nn, mm, np, mp, nnn, mmm, nnp, mmp, nxv3, nxv4;
   float qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float vx, vy, vz, v1, v2, v3, v4;
   nxv3 = 3*nxv;
   nxv4 = 4*nxv;
   qtmh = 0.5*qbm*dt;
   dti = 1.0/dt;
   for (j = 0; j < nop; j++) {
/* find interpolation weights */
      nnn = part[idimp*j];
      mmm = part[1+idimp*j];
      dxp = part[idimp*j] - (float) nnn;
      dyp = part[1+idimp*j] - (float) mmm;
      nn = 3*nnn;
      mm = nxv3*mmm;
      amx = 1.0 - dxp;
      mp = mm + nxv3;
      amy = 1.0 - dyp;
      np = nn + 3;
      nnn = 4*nnn;
      mmm = nxv4*mmm;
      mmp = mmm + nxv4;
      nnp = nnn + 4;
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
      vx = part[2+idimp*j];
      vy = part[3+idimp*j];
      vz = part[4+idimp*j];
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
      dx = dxp*dyp;
      dy = amx*dyp;
      v1 = ox*ox - oy*oy;
      v2 = ox*oy;
      v3 = oz*ox;
      v4 = oz*oy;
      amu[nnp+mmp] += v1*dx;
      amu[1+nnp+mmp] += v2*dx;
      amu[2+nnp+mmp] += v3*dx;
      amu[3+nnp+mmp] += v4*dx;
      dcu[np+mp] += vx*dx;
      dcu[1+np+mp] += vy*dx;
      dcu[2+np+mp] += vz*dx;
      cu[np+mp] += ox*dx;
      cu[1+np+mp] += oy*dx;
      cu[2+np+mp] += oz*dx;
      dx = dxp*amy;
      amu[nnn+mmp] += v1*dy;
      amu[1+nnn+mmp] += v2*dy;
      amu[2+nnn+mmp] += v3*dy;
      amu[3+nnn+mmp] += v4*dy;
      dcu[nn+mp] += vx*dy;
      dcu[1+nn+mp] += vy*dy;
      dcu[2+nn+mp] += vz*dy;
      cu[nn+mp] += ox*dy;
      cu[1+nn+mp] += oy*dy;
      cu[2+nn+mp] += oz*dy;
      dy = amx*amy;
      amu[nnp+mmm] += v1*dx;
      amu[1+nnp+mmm] += v2*dx;
      amu[2+nnp+mmm] += v3*dx;
      amu[3+nnp+mmm] += v4*dx;
      dcu[np+mm] += vx*dx;
      dcu[1+np+mm] += vy*dx;
      dcu[2+np+mm] += vz*dx;
      cu[np+mm] += ox*dx;
      cu[1+np+mm] += oy*dx;
      cu[2+np+mm] += oz*dx;
      amu[nnn+mmm] += v1*dy;
      amu[1+nnn+mmm] += v2*dy;
      amu[2+nnn+mmm] += v3*dy;
      amu[3+nnn+mmm] += v4*dy;
      dcu[nn+mm] += vx*dy;
      dcu[1+nn+mm] += vy*dy;
      dcu[2+nn+mm] += vz*dy;
      cu[nn+mm] += ox*dy;
      cu[1+nn+mm] += oy*dy;
      cu[2+nn+mm] += oz*dy;
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
   for (k = 0; k < ny; k++) {
      for (j = 0; j < nx; j++) {
         at1 = qbme*qe[j+nxe*k];
         tpmax = at1 > tpmax ? at1 : tpmax;
         tpmin = at1 < tpmin ? at1 : tpmin;
      }
   }
   *wpmax = tpmax;
   *wpmin = tpmin;
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
         at1 = 1.0/(dkx*dkx + dky2);
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
void cbbpois23(float complex cu[], float complex bxy[],
               float complex ffc[], float ci, float *wm, int nx, int ny,
               int nxvh, int nyv, int nxhd, int nyhd) {
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
   double wp;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
   ci2 = ci*ci;
/* calculate smoothed magnetic field and sum field energy */
   wp = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
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
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
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
   }
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
   *wm = wp*(float) (nx*ny);
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
void cdcuperp23(float complex dcu[], float complex amu[], int nx,
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
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
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
void cadcuperp23(float complex dcu[], float complex amu[], int nx,
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
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
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
void cepois23(float complex dcu[], float complex exy[], int isign,
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
   double wp;
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
        goto L80;
   wp = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
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
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
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
   }
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
   *wf = (wp*(float) (nx*ny))/crealf(ffe[0]);
   return;
/* calculate unsmoothed transverse electric field and sum field energy */
L80: wp = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
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
  }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
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
  }
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
   *wf = (wp*(float) (nx*ny))/crealf(ffe[0]);
   return;
}

/*--------------------------------------------------------------------*/
void caddvrfield2(float a[], float b[], float c[], int ndim, int nxe,
                  int nye) {
/* this subroutine calculates a = b + c for real vector fields
local data                                                 */
   int i, j, k, nnxe;
   nnxe = ndim*nxe;
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
void cswapc2n(float f[], float s[], int isign, int nxh, int nyi,
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

      dimension f(ndim,2*nxhd,nyd), s(2*ndim*nxhd)

   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = nyi-1; k < nyt; k++) {

local data                                                            */
   int i, j, k, ioff, nk;
/* swap complex components */
/* real to complex */
   if (isign < 0){
      for (k = nyi-1; k < nyt; k++) {
         nk = 2*ndim*nxhd*k;
         for (j = 0; j < nxh; j++) {
            ioff = 2*ndim*j;
            for (i = 0; i < ndim; i++) {
               s[2*i+ioff] = f[i+ndim*(2*j)+nk];
               s[2*i+ioff+1] = f[i+ndim*(2*j+1)+nk];
            }
         }
         for (j = 0; j < nxh; j++) {
            ioff = 2*ndim*j;
            for (i = 0; i < ndim; i++) {
               f[i+ndim*(2*j)+nk] = s[i+ioff];
            }
            ioff += ndim;
            for (i = 0; i < ndim; i++) {
               f[i+ndim*(2*j+1)+nk] = s[i+ioff];
            }
         }
      }
   }
/* complex to real */
   else if (isign > 0) {
      for (k = nyi-1; k < nyt; k++) {
         nk = 2*ndim*nxhd*k;
         for (j = 0; j < nxh; j++) {
            ioff = 2*ndim*j;
            for (i = 0; i < ndim; i++) {
               s[i+ioff] = f[i+ndim*(2*j)+nk];
            }
            ioff += ndim;
            for (i = 0; i < ndim; i++) {
               s[i+ioff] = f[i+ndim*(2*j+1)+nk];
            }
         }
         for (j = 0; j < nxh; j++) {
            ioff = 2*ndim*j;
            for (i = 0; i < ndim; i++) {
               f[i+ndim*(2*j)+nk] = s[2*i+ioff];
               f[i+ndim*(2*j+1)+nk] = s[2*i+ioff+1];
            }
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rnx(float complex f[], float complex ss[], int isign,
              int mixup[], float complex sct[], int indx, int indy,
              int nyi, int nyp, int nxhd, int nyd, int ndim, int nxhyd,
              int nxyhd) {
/* this subroutine performs the x part of N two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   y, using complex arithmetic, where N = ndim
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
   f[k][j][0:N-1] = mode j,k, where 0 <= j < nx/2 and 0 <= k < ny,
   except for f[k][0][0:N-1] =  mode nx/2,k-1, where ny/2+1 <= k < ny, and
   imag(f[0][0][0:N-1]) = real part of mode nx/2,0 and
   imag(f[ny/2][0][0:N-1]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
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
      goto L140;
/* inverse fourier transform */
/* swap complex components */
   cswapc2n((float *)f,(float *)ss,isign,nxh,nyi,nyt,nxhd,nyd,ndim);
/* bit-reverse array elements in x */
   nrx = nxhy/nxh;
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = nyi-1; k < nyt; k++) {
         joff = ndim*nxhd*k;
         for (jj = 0; jj < ndim; jj++) {
            t1 = f[jj+ndim*j1+joff];
            f[jj+ndim*j1+joff] = f[jj+ndim*j+joff];
            f[jj+ndim*j+joff] = t1;
         }
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
               joff = ndim*nxhd*i;
               for (jj = 0; jj < ndim; jj++) {
                  t2 = t1*f[jj+ndim*j2+joff];
                  f[jj+ndim*j2+joff] = f[jj+ndim*j1+joff] - t2;
                  f[jj+ndim*j1+joff] += t2;
               }
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
         joff = ndim*nxhd*k;
         for (jj = 0; jj < ndim; jj++) {
            t2 = conjf(f[jj+ndim*(nxh-j)+joff]);
            t1 = f[jj+ndim*j+joff] + t2;
            t2 = (f[jj+ndim*j+joff] - t2)*t3;
            f[jj+ndim*j+joff] = ani*(t1 + t2);
            f[jj+ndim*(nxh-j)+joff] = ani*conjf(t1 - t2);
         }
      }
   }
   ani = 2.0*ani;
   for (k = nyi-1; k < nyt; k++) {
      joff = ndim*nxhd*k;
      for (jj = 0; jj < ndim; jj++) {
         f[jj+ndim*nxhh+joff] = ani*conjf(f[jj+ndim*nxhh+joff]);
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
         joff = ndim*nxhd*k;
         for (jj = 0; jj < ndim; jj++) {
            t2 = conjf(f[jj+ndim*(nxh-j)+joff]);
            t1 = f[jj+ndim*j+joff] + t2;
            t2 = (f[jj+ndim*j+joff] - t2)*t3;
            f[jj+ndim*j+joff] = t1 + t2;
            f[jj+ndim*(nxh-j)+joff] = conjf(t1 - t2);
         }
      }
   }
   for (k = nyi-1; k < nyt; k++) {
      joff = ndim*nxhd*k;
      for (jj = 0; jj < ndim; jj++) {
         f[jj+ndim*nxhh+joff] = 2.0*conjf(f[jj+ndim*nxhh+joff]);
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
         joff = ndim*nxhd*k;
         for (jj = 0; jj < ndim; jj++) {
            t1 = f[jj+ndim*j1+joff];
            f[jj+ndim*j1+joff] = f[jj+ndim*j+joff];
            f[jj+ndim*j+joff] = t1;
         }
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
               joff = ndim*nxhd*i;
               for (jj = 0; jj < ndim; jj++) {
                  t2 = t1*f[jj+ndim*j2+joff];
                  f[jj+ndim*j2+joff] = f[jj+ndim*j1+joff] - t2;
                  f[jj+ndim*j1+joff] += t2;
               }
            }
         }
      }
      ns = ns2;
   }
/* swap complex components */
   cswapc2n((float *)f,(float *)ss,isign,nxh,nyi,nyt,nxhd,nyd,ndim);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rny(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int ndim, int nxhyd, int nxyhd) {
/* this subroutine performs the y part of N two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   x, using complex arithmetic, where N = ndim
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
   f[k][j][0:N-1] = mode j,k, where 0 <= j < nx/2 and 0 <= k < ny,
   except for f[k][0][0:N-1] =  mode nx/2,k-1, where ny/2+1 <= k < ny, and
   imag(f[0][0][0:N-1]) = real part of mode nx/2,0 and
   imag(f[ny/2][0][0:N-1]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
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
      goto L90;
/* inverse fourier transform */
   nry = nxhy/ny;
/* bit-reverse array elements in y */
   for (k = 0; k < ny; k++) {
      joff = ndim*nxhd*k;
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      k1 = ndim*nxhd*k1;
      for (j = nxi-1; j < nxt; j++) {
         for (jj = 0; jj < ndim; jj++) {
            t1 = f[jj+ndim*j+k1];
            f[jj+ndim*j+k1] = f[jj+ndim*j+joff];
            f[jj+ndim*j+joff] = t1;
         }
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
            j1 = ndim*nxhd*(j + k1);
            j2 = ndim*nxhd*(j + k2);
            t1 = sct[kmr*j];
            for (i = nxi-1; i < nxt; i++) {
               for (jj = 0; jj < ndim; jj++) {
                  t2 = t1*f[jj+ndim*i+j2];
                  f[jj+ndim*i+j2] = f[jj+ndim*i+j1] - t2;
                  f[jj+ndim*i+j1] += t2;
               }
            }
         }
      }
      ns = ns2;
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
/* scramble modes kx = 0, nx/2 */
L90: for (k = 1; k < nyh; k++) {
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
/* bit-reverse array elements in y */
   nry = nxhy/ny;
   for (k = 0; k < ny; k++) {
      joff = ndim*nxhd*k;
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      k1 = ndim*nxhd*k1;
      for (j = nxi-1; j < nxt; j++) {
         for (jj = 0; jj < ndim; jj++) {
            t1 = f[jj+ndim*j+k1];
            f[jj+ndim*j+k1] = f[jj+ndim*j+joff];
            f[jj+ndim*j+joff] = t1;
         }
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
            j1 = ndim*nxhd*(j + k1);
            j2 = ndim*nxhd*(j + k2);
            t1 = conjf(sct[kmr*j]);
            for (i = nxi-1; i < nxt; i++) {
               for (jj = 0; jj < ndim; jj++) {
                  t2 = t1*f[jj+ndim*i+j2];
                  f[jj+ndim*i+j2] = f[jj+ndim*i+j1] - t2;
                  f[jj+ndim*i+j1] += t2;
               }
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
void cwfft2r3(float complex f[], int isign, int mixup[],
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

/*--------------------------------------------------------------------*/
void cwfft2rn(float complex f[], float complex ss[], int isign,
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
      cfft2rnx(f,ss,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,ndim,
               nxhyd,nxyhd);
/* perform y fft */
      cfft2rny(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,ndim,nxhyd,
               nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfft2rny(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,ndim,nxhyd,
               nxyhd);
/* perform x fft */
      cfft2rnx(f,ss,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,ndim,
               nxhyd,nxyhd);
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
void cgbpush23l_(float *part, float *fxy, float *bxy, float *qbm,
                 float *dt, float *dtc, float *ek, int *idimp, int *nop,
                 int *nx, int *ny, int *nxv, int *nyv, int *ipbc) {
   cgbpush23l(part,fxy,bxy,*qbm,*dt,*dtc,ek,*idimp,*nop,*nx,*ny,*nxv,
              *nyv,*ipbc);
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
void cgmjpost2l_(float *part, float *amu, float *qm, int *nop,
                 int *idimp, int *nxv, int *nyv) {
   cgmjpost2l(part,amu,*qm,*nop,*idimp,*nxv,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cgdjpost2l_(float *part, float *fxy, float *bxy, float *dcu,
                 float *amu, float *qm, float *qbm, float *dt,
                 int *idimp, int *nop, int *nxv, int *nyv) {
   cgdjpost2l(part,fxy,bxy,dcu,amu,*qm,*qbm,*dt,*idimp,*nop,*nxv,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cgdcjpost2l_(float *part, float *fxy, float *bxy, float *cu,
                  float *dcu, float *amu, float *qm, float *qbm,
                  float *dt, int *idimp, int *nop, int *nxv, int *nyv) {
   cgdcjpost2l(part,fxy,bxy,cu,dcu,amu,*qm,*qbm,*dt,*idimp,*nop,*nxv,
               *nyv);
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
void cbbpois23_(float complex *cu, float complex *bxy,
                float complex *ffc, float *ci, float *wm, int *nx,
                int *ny, int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   cbbpois23(cu,bxy,ffc,*ci,wm,*nx,*ny,*nxvh,*nyv,*nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cbaddext2_(float *bxy, float *omx, float *omy, float *omz, int *nx,
                int *ny, int *nxe, int *nye) {
 cbaddext2(bxy,*omx,*omy,*omz,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void cdcuperp23_(float complex *dcu, float complex *amu, int *nx,
                 int *ny, int *nxvh, int *nyv) {
 cdcuperp23(dcu,amu,*nx,*ny,*nxvh,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cadcuperp23_(float complex *dcu, float complex *amu, int *nx,
                  int *ny, int *nxvh, int *nyv) {
   cadcuperp23(dcu,amu,*nx,*ny,*nxvh,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cepois23_(float complex *dcu, float complex *exy, int *isign,
               float complex *ffe, float *ax, float *ay, float *affp,
               float *wp0, float *ci, float *wf, int *nx, int *ny,
               int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   cepois23(dcu,exy,*isign,ffe,*ax,*ay,*affp,*wp0,*ci,wf,*nx,*ny,*nxvh,
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

/*--------------------------------------------------------------------*/
void cwfft2rn_(float complex *f, float complex *ss, int *isign,
              int *mixup, float complex *sct, int *indx, int *indy,
              int *nxhd, int *nyd, int *ndim, int *nxhyd, int *nxyhd) {
   cwfft2rn(f,ss,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*ndim,*nxhyd,
            *nxyhd);
   return;
}
