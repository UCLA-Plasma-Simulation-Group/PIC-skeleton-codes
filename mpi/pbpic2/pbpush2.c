/* C Library for Skeleton 2-1/2D Electromagnetic MPI PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "pbpush2.h"
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
void cpdistr2h(float part[], float edges[], int *npp, int nps,
               float vtx, float vty, float vtz, float vdx, float vdy,
               float vdz, int npx, int npy, int nx, int ny, int idimp,
               int npmax, int idps, int ipbc, int *ierr) {
/* for 2-1/2d code, this subroutine calculates initial particle
   co-ordinates and velocities with uniform density and maxwellian
   velocity with drift for distributed data.
   input: all except part, ierr, output: part, npp, ierr
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   part[n][2] = velocity vx of particle n in partition
   part[n][3] = velocity vy of particle n in partition
   part[n][4] = velocity vz of particle n in partition
   edges[0] = lower boundary of particle partition
   edges[1] = upper boundary of particle partition
   npp = number of particles in partition
   nps = starting address of particles in partition
   vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
   vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
   npx/npy = initial number of particles distributed in x/y direction
   nx/ny = system length in x/y direction
   idimp = size of phase space = 5
   npmax = maximum number of particles in each partition
   idps = number of partition boundaries
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   ierr = (0,1) = (no,yes) error condition exists
   ranorm = gaussian random number with zero mean and unit variance
   with spatial decomposition
local data                                                            */
   int j, k, npt, k1, npxyp;
   float edgelx, edgely, at1, at2, xt, yt, vxt, vyt, vzt;
   double dnpx, dnpxy, dt1;
   int ierr1[1], iwork1[1];
   double sum4[4], work4[4];
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
         vzt = vtz*ranorm();
         if ((yt >= edges[0]) && (yt < edges[1])) {
            if (npt < npmax) {
               k1 = idimp*npt;
               part[k1] = xt;
               part[1+k1] = yt;
               part[2+k1] = vxt;
               part[3+k1] = vyt;
               part[4+k1] = vzt;
               npt += 1;
            }
            else
               *ierr += 1;
         }
      }
   }
   npxyp = 0;
/* add correct drift */
   sum4[0] = 0.0;
   sum4[1] = 0.0;
   sum4[2] = 0.0;
   for (j = nps-1; j < npt; j++) {
      npxyp += 1;
      sum4[0] += part[2+idimp*j];
      sum4[1] += part[3+idimp*j];
      sum4[2] += part[4+idimp*j];
   }
   sum4[3] = npxyp;
   cppdsum(sum4,work4,4);
   dnpxy = sum4[3];
   ierr1[0] = *ierr;
   cppimax(ierr1,iwork1,1);
   *ierr = ierr1[0];
   dt1 = 1.0/dnpxy;
   sum4[0] = dt1*sum4[0] - vdx;
   sum4[1] = dt1*sum4[1] - vdy;
   sum4[2] = dt1*sum4[2] - vdz;
   for (j = nps-1; j < npt; j++) {
      part[2+idimp*j] -= sum4[0];
      part[3+idimp*j] -= sum4[1];
      part[4+idimp*j] -= sum4[2];
   }
/* process errors */
   dnpxy -= dnpx*(double) npy;
   if (dnpxy != 0.0)
      *ierr = dnpxy;
   *npp = npt;
   return;
}

/*--------------------------------------------------------------------*/
void cppgbpush23l(float part[], float fxy[], float bxy[], float edges[],
                  int npp, int noff, int ihole[], float qbm, float dt, 
                  float dtc, float *ek, int nx, int ny, int idimp,
                  int npmax, int nxv, int nypmx, int idps, int ntmax,
                  int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space with magnetic field. Using the Boris Mover.
   scalar version using guard cells, for distributed data
   also determines list of particles which are leaving this processor
   117 flops/particle, 1 divide, 25 loads, 5 stores
   input: all except ihole, output: part, ihole, ek
   velocity equations used are:
   vx(t+dt/2) = rot[0]*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[1]*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot[2]*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   vy(t+dt/2) = rot[3]*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[4]*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot[5]*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   vz(t+dt/2) = rot[6]*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[7]*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
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
   fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   part[n][2] = velocity vx of particle n in partition
   part[n][3] = velocity vy of particle n in partition
   part[n][4] = velocity vz of particle n in partition
   fxy[k][j][0] = x component of force/charge at grid (j,kk)
   fxy[k][j][1] = y component of force/charge at grid (j,kk)
   fxy[k][j][2] = z component of force/charge at grid (j,kk)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,kk)
   bxy[k][j][1] = y component of magnetic field at grid (j,kk)
   bxy[k][j][2] = z component of magnetic field at grid (j,kk)
   that is, the convolution of magnetic field over particle shape
   where kk = k + noff
   edges[0:1] = lower:upper boundary of particle partition
   npp = number of particles in partition
   noff = lowermost global gridpoint in particle partition.
   ihole = location of hole left in particle arrays
   ihole[0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
        .25*(vz(t+dt/2) + vz(t-dt/2))**2)
   nx/ny = system length in x/y direction
   idimp = size of phase space = 5
   npmax = maximum number of particles in each partition
   nxv = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   idps = number of partition boundaries
   ntmax = size of hole array for particles leaving processors
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int mnoff, j, nn, mm, np, mp, ih, nh, nxv3;
   float qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt;
   float anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   double sum1;
   nxv3 = 3*nxv;
   qtmh = 0.5*qbm*dt;
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
      nn = 3*nn;
      mm = nxv3*(mm - mnoff);
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
   *ek += 0.5*sum1;
   return;
}

/*--------------------------------------------------------------------*/
void cdppgbpush23l(double part[], double fxy[], double bxy[], 
                   double edges[], int npp, int noff, int ihole[],
                   double qbm, double dt, double dtc, double *ek,
                   int nx, int ny, int idimp, int npmax, int nxv,
                   int nypmx, int idps, int ntmax, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space with magnetic field. Using the Boris Mover.
   scalar version using guard cells, for distributed data
   also determines list of particles which are leaving this processor
   117 flops/particle, 1 divide, 25 loads, 5 stores
   input: all except ihole, output: part, ihole, ek
   velocity equations used are:
   vx(t+dt/2) = rot[0]*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[1]*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot[2]*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   vy(t+dt/2) = rot[3]*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[4]*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot[5]*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   vz(t+dt/2) = rot[6]*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[7]*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
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
   fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   part[n][2] = velocity vx of particle n in partition
   part[n][3] = velocity vy of particle n in partition
   part[n][4] = velocity vz of particle n in partition
   fxy[k][j][0] = x component of force/charge at grid (j,kk)
   fxy[k][j][1] = y component of force/charge at grid (j,kk)
   fxy[k][j][2] = z component of force/charge at grid (j,kk)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,kk)
   bxy[k][j][1] = y component of magnetic field at grid (j,kk)
   bxy[k][j][2] = z component of magnetic field at grid (j,kk)
   that is, the convolution of magnetic field over particle shape
   where kk = k + noff
   edges[0:1] = lower:upper boundary of particle partition
   npp = number of particles in partition
   noff = lowermost global gridpoint in particle partition.
   ihole = location of hole left in particle arrays
   ihole[0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
        .25*(vz(t+dt/2) + vz(t-dt/2))**2)
   nx/ny = system length in x/y direction
   idimp = size of phase space = 5
   npmax = maximum number of particles in each partition
   nxv = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   idps = number of partition boundaries
   ntmax = size of hole array for particles leaving processors
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int mnoff, j, nn, mm, np, mp, ih, nh, nxv3;
   double qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   double dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt;
   double anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   double sum1;
   nxv3 = 3*nxv;
   qtmh = 0.5*qbm*dt;
   sum1 = 0.0;
/* set boundary values */
   edgelx = 0.0;
   edgely = 1.0;
   edgerx = (double) nx;
   edgery = (double) (ny-1);
   if ((ipbc==2) || (ipbc==3)) {
      edgelx = 1.0;
      edgerx = (double) (nx-1);
   }
   mnoff = noff;
   ih = 0;
   nh = 0;
   for (j = 0; j < npp; j++) {
/* find interpolation weights */
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = part[idimp*j] - (double) nn;
      dyp = part[1+idimp*j] - (double) mm;
      nn = 3*nn;
      mm = nxv3*(mm - mnoff);
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
   *ek += 0.5*sum1;
   return;
}

/*--------------------------------------------------------------------*/
void cppgrbpush23l(float part[], float fxy[], float bxy[],
                   float edges[], int npp, int noff, int ihole[],
                   float qbm, float dt, float dtc, float ci, float *ek,
                   int nx, int ny, int idimp, int npmax, int nxv,
                   int nypmx, int idps, int ntmax, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   scalar version using guard cells, for distributed data
   also determines list of particles which are leaving this processor
   129 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
   input: all except ihole, output: part, ihole, ek
   momentum equations used are:
   px(t+dt/2) = rot[0]*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[1]*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot[2]*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   py(t+dt/2) = rot[3]*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[4]*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot[5]*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   pz(t+dt/2) = rot[6]*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[7]*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot[8]*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
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
   fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   part[n][2] = momentum px of particle n in partition
   part[n][3] = momentum py of particle n in partition
   part[n][4] = momentum pz of particle n in partition
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   edges[0:1] = lower:upper boundary of particle partition
   npp = number of particles in partition
   noff = lowermost global gridpoint in particle partition.
   ihole = location of hole left in particle arrays
   ihole(1) = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   nx/ny = system length in x/y direction
   idimp = size of phase space = 5
   npmax = maximum number of particles in each partition
   nxv = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   idps = number of partition boundaries
   ntmax = size of hole array for particles leaving processors
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int mnoff, j, nn, mm, np, mp, ih, nh, nxv3;
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
      nn = 3*nn;
      mm = nxv3*(mm - mnoff);
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
   *ek += sum1;
   return;
}

/*--------------------------------------------------------------------*/
void cdppgrbpush23l(double part[], double fxy[], double bxy[],
                    double edges[], int npp, int noff, int ihole[],
                    double qbm, double dt, double dtc, double ci, 
                    double *ek, int nx, int ny, int idimp, int npmax,
                    int nxv, int nypmx, int idps, int ntmax, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   scalar version using guard cells, for distributed data
   also determines list of particles which are leaving this processor
   129 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
   input: all, output: part, ek
   momentum equations used are:
   px(t+dt/2) = rot[0]*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[1]*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot[2]*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fx(x(t),y(t))*dt)
   py(t+dt/2) = rot[3]*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[4]*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot[5]*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
      .5*(q/m)*fy(x(t),y(t))*dt)
   pz(t+dt/2) = rot[6]*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
      rot[7]*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
      rot[8]*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
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
   fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
   bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   part[n][2] = momentum px of particle n in partition
   part[n][3] = momentum py of particle n in partition
   part[n][4] = momentum pz of particle n in partition
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   edges[0:1] = lower:upper boundary of particle partition
   npp = number of particles in partition
   noff = lowermost global gridpoint in particle partition.
   ihole = location of hole left in particle arrays
   ihole(1) = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   nx/ny = system length in x/y direction
   idimp = size of phase space = 5
   npmax = maximum number of particles in each partition
   nxv = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   idps = number of partition boundaries
   ntmax = size of hole array for particles leaving processors
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int mnoff, j, nn, mm, np, mp, ih, nh, nxv3;
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
   edgely = 1.0;
   edgerx = (double) nx;
   edgery = (double) (ny-1);
   if ((ipbc==2) || (ipbc==3)) {
      edgelx = 1.0;
      edgerx = (double) (nx-1);
   }
   mnoff = noff;
   ih = 0;
   nh = 0;
   for (j = 0; j < npp; j++) {
/* find interpolation weights */
      nn = part[idimp*j];
      mm = part[1+idimp*j];
      dxp = part[idimp*j] - (double) nn;
      dyp = part[1+idimp*j] - (double) mm;
      nn = 3*nn;
      mm = nxv3*(mm - mnoff);
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
   *ek += sum1;
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
void cppgjpost2l(float part[], float cu[], float edges[], int npp,
                 int noff, int ihole[], float qm, float dt, int nx,
                 int ny, int idimp, int npmax, int nxv, int nypmx,
                 int idps, int ntmax, int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation, and distributed data.
   in addition, particle positions are advanced a half time-step
   also determines list of particles which are leaving this processor
   scalar version using guard cells, for distributed data
   35 flops/particle, 17 loads, 14 stores
   input: all except ihole, output: part, ihole, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*vi, where i = x,y,z
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   part[n][2] = velocity vx of particle n in partition
   part[n][3] = velocity vy of particle n in partition
   part[n][4] = velocity vz of particle n in partition
   cu[k][j][i] = ith component of current density at grid point j,kk,
   where kk = k + noff
   edges[0:1] = lower:upper boundary of particle partition
   npp = number of particles in partition
   noff = lowermost global gridpoint in particle partition.
   ihole = location of hole left in particle arrays
   ihole[0] = ih, number of holes left (error, if negative)
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nx/ny = system length in x/y direction
   idimp = size of phase space = 5
   npmax = maximum number of particles in each partition
   nxv = first dimension of current array, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   idps = number of partition boundaries
   ntmax =  size of hole array for particles leaving processors
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int mnoff, j, nn, mm, np, mp, ih, nh, nxv3;
   float edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, vx, vy, vz;
   nxv3 = 3*nxv;
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
      dxp = qm*(part[idimp*j] - (float) nn);
      dyp = part[1+idimp*j] - (float) mm;
      nn = 3*nn;
      mm = nxv3*(mm - mnoff);
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
   return;
}

/*--------------------------------------------------------------------*/
void cppgrjpost2l(float part[], float cu[], float edges[], int npp,
                  int noff, int ihole[], float qm, float dt, float ci,
                  int nx, int ny, int idimp, int npmax, int nxv,
                  int nypmx, int idps, int ntmax, int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation for relativistic particles,
   in addition, particle positions are advanced a half time-step
   also determines list of particles which are leaving this processor
   scalar version using guard cells, for distributed data
   45 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
   input: all except ihole, output: part, ihole, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*pi*gami, where i = x,y,z
   where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   part[n][2] = momentum px of particle n in partition
   part[n][3] = momentum py of particle n in partition
   part[n][4] = momentum pz of particle n in partition
   cu[k][j][i] = ith component of current density at grid point j,kk,
   where kk = k + noff
   edges[0:1] = lower:upper boundary of particle partition
   npp = number of particles in partition
   noff = lowermost global gridpoint in particle partition.
   ihole = location of hole left in particle arrays
   ihole[0] = ih, number of holes left (error, if negative)
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
   nx/ny = system length in x/y direction
   idimp = size of phase space = 5
   npmax = maximum number of particles in each partition
   nxv = first dimension of current array, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   idps = number of partition boundaries
   ntmax =  size of hole array for particles leaving processors
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int mnoff, j, nn, mm, np, mp, ih, nh, nxv3;
   float ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, vx, vy, vz, p2, gami;
   nxv3 = 3*nxv;
   ci2 = ci*ci;
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
      mm = nxv3*(mm - mnoff);
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
void cppcguard2xl(float fxy[], int myp, int nx, int ndim, int nxe,
                  int nypmx) {
/* replicate extended periodic vector field in x direction
   linear interpolation, for distributed data
   myp = number of full or partial grids in particle partition
   nx = system length in x direction
   ndim = leading dimension of array fxy
   nxe = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells
local data                                                 */
   int i, k, kk, myp1;
/* replicate edges of extended field */
   myp1 = myp + 1;
   for (k = 0; k < myp1; k++) {
      kk = ndim*nxe*k;
      for (i = 0; i < ndim; i++) {
         fxy[i+ndim*nx+kk] = fxy[i+kk];
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppaguard2xl(float q[], int myp, int nx, int nxe, int nypmx) {
/* accumulate extended periodic scalar field in x direction
   linear interpolation, for distributed data
   myp = number of full or partial grids in particle partition
   nx = system length in x direction
   nxe = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells
local data                                                 */
   int k, myp1;
/* accumulate edges of extended field */
   myp1 = myp + 1;
   for (k = 0; k < myp1; k++) {
      q[nxe*k] += q[nx+nxe*k];
      q[nx+nxe*k] = 0.0;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppacguard2xl(float cu[], int myp, int nx, int ndim, int nxe,
                   int nypmx) {
/* accumulate extended periodic vector field in x direction
   linear interpolation, for distributed data
   myp = number of full or partial grids in particle partition
   nx = system length in x direction
   ndim = leading dimension of array fxy
   nxe = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells
      implicit none
      real cu
      integer myp, nx, ndim, nxe, nypmx
      dimension cu(ndim,nxe,nypmx)
local data                                                 */
   int i, k, kk, myp1;
/* accumulate edges of extended field */
   myp1 = myp + 1;
   for (k = 0; k < myp1; k++) {
      kk = ndim*nxe*k;
      for (i = 0; i < ndim; i++) {
         cu[i+kk] += cu[i+ndim*nx+kk];
         cu[i+ndim*nx+kk] = 0.0;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppois23(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int kstrt, int nyv, int kxp,
              int nyhd) {
/* this subroutine solves 2d poisson's equation in fourier space for
   force/charge (or convolution of electric field over particle shape)
   with periodic boundary conditions.  Zeros out z component.
   for distributed data.
   for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp,nyhd,
   output: ffc
   for isign /= 0, input: q,ffc,isign,nx,ny,kstrt,nyv,kxp,nyhd,
   output: fxy,we
   approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
   where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
   the equation used is:
   fx[ky][kx] = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
   fy[ky][kx] = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
   fz[ky][kx] = zero,
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   g[ky][kx] = (affp/(kx**2+ky**2))*s(kx,ky),
   s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
   fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
   fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
   q[k][j] = complex charge density for fourier mode (jj-1,k-1)
   fxy[k][j][0] = x component of complex force/charge,
   fxy[k][j][1] = y component of complex force/charge,
   fxy[k][j][2] = zero,
   for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
   kxp = number of data values per block
   kstrt = starting data block number
   if isign = 0, form factor array is prepared
   if isign is not equal to 0, force/charge is calculated.
   aimag(ffc[k][j]) = finite-size particle shape factor s
   real(ffc[k][j])) = potential green's function g
   for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
   ax/ay = half-width of particle in x/y direction
   affp = normalization constant = nx*ny/np, where np=number of particles
   electric field energy is also calculated, using
   we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
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
            fxy[3*k+3*jk] = at2*zt1;
            fxy[1+3*k+3*jk] = at3*zt1;
            fxy[2+3*k+3*jk] = zero;
            fxy[3*k1+3*jk] = at2*zt2;
            fxy[1+3*k1+3*jk] = -at3*zt2;
            fxy[2+3*k1+3*jk] = zero;
            wp += at1*(q[k+jk]*conjf(q[k+jk])
                  + q[k1+jk]*conjf(q[k1+jk]));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         at1 = crealf(ffc[jj])*cimagf(ffc[jj]);
         at3 = dkx*at1;
         zt1 = cimagf(q[jk]) - crealf(q[jk])*_Complex_I;
         fxy[3*jk] = at3*zt1;
         fxy[1+3*jk] = zero;
         fxy[2+3*jk] = zero;
         fxy[3*k1+3*jk] = zero;
         fxy[1+3*k1+3*jk] = zero;
         fxy[2+3*k1+3*jk] = zero;
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
         fxy[3*k] = zero;
         fxy[1+3*k] = at2*zt1;
         fxy[2+3*k] = zero;
         fxy[3*k1] = zero;
         fxy[1+3*k1] = zero;
         fxy[2+3*k1] = zero;
         wp += at1*(q[k]*conjf(q[k]));
      }
      k1 = 3*nyh;
      fxy[0] = zero;
      fxy[1] = zero;
      fxy[2] = zero;
      fxy[k1] = zero;
      fxy[1+k1] = zero;
      fxy[2+k1] = zero;
   }
L70:
   *we = wp*((float) nx)*((float) ny);
   return;
}

/*--------------------------------------------------------------------*/
void cppcuperp2(float complex cu[], int nx, int ny, int kstrt, int nyv,
                int kxp) {
/* this subroutine calculates the transverse current in fourier space
   input: all, output: cu
   approximate flop count is: 36*nxc*nyc
   and nxc*nyc divides
   where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
   the transverse current is calculated using the equation:
   cux[ky][kx] = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
   cuy[ky][kx] = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
   and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
   cu[j][k][i] = i-th component of complex current density and
   for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
   nx/ny = system length in x/y direction
   kstrt = starting data block number
   nyv = first dimension of field arrays, must be >= ny
   kxp = number of data values per block
local data                                                          */
   int nxh, nyh, ks, joff, kxps, j, jk, k, k1;
   float dnx, dny, dkx, dky, dkx2, at1;
   float complex zero, zt1;
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
/* calculate transverse part of current */
   if (kstrt > nxh)
      return;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (j = 0; j < kxps; j++) {
      dkx = dnx*(float) (j + joff);
      dkx2 = dkx*dkx;
      jk = nyv*j;
      if ((j+joff) > 0) {
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            dky = dny*(float) k;
            at1 = 1.0/(dky*dky + dkx2);
            zt1 = at1*(dkx*cu[3*k+3*jk] + dky*cu[1+3*k+3*jk]);
            cu[3*k+3*jk] -= dkx*zt1;
            cu[1+3*k+3*jk] -= dky*zt1;
            zt1 = at1*(dkx*cu[3*k1+3*jk] - dky*cu[1+3*k1+3*jk]);
            cu[3*k1+3*jk] -= dkx*zt1;
            cu[1+3*k1+3*jk] += dky*zt1;
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         cu[3*jk] = zero;
         cu[3*k1+3*jk] = zero;
         cu[1+3*k1+3*jk] = zero;
      }

   }
/* mode numbers kx = 0, nx/2 */
   if (ks==0) {
      for (k = 1; k < nyh; k++) {
         k1 = ny - k;
         cu[1+3*k] = zero;
         cu[3*k1] = zero;
         cu[1+3*k1] = zero;
      }
      k1 = 3*nyh;
      cu[0] = zero;
      cu[1] = zero;
      cu[k1] = zero;
      cu[1+k1] = zero;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cippbpoisp23(float complex cu[], float complex bxy[],
                  float complex ffc[], float ci, float *wm, int nx,
                  int ny, int kstrt, int nyv, int kxp, int nyhd) {
/* this subroutine solves 2-1/2d poisson's equation in fourier space for
   magnetic field with periodic boundary conditions for distributed data.
   input: cu,ffc,ci,nx,ny,kstrt,nyv,kxp,jblok,nyhd, output: bxy,wm
   approximate flop count is: 85*nxc*nyc + 36*(nxc + nyc)
   where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
   magnetic field is calculated using the equations:
   bx[ky][kx] = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
   by[ky][kx] = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
   bz[ky][kx] = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   g[ky][kx] = (affp/(kx**2+ky**2))*s(kx,ky),
   s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
   bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
   bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0,
   bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
   cu[j][k][i] = i-th component of complex current density and
   bxy[j][k][i] = i-th component of complex magnetic field,
   for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
   kxp = number of data values per block
   kstrt = starting data block number
   imag(ffc[j][k]) = finite-size particle shape factor s
   real(ffc[j][k]) = potential green's function g
   for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
   ci = reciprocal of velocity of light
   magnetic field energy is also calculated, using
   wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
      |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
   affp = normalization constant = nx*ny/np, where np=number of particles
   this expression is valid only if the current is divergence-free
   nx/ny = system length in x/y direction
   nyv = second dimension of field arrays, must be >= ny
   nyhd = first dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, ks, joff, kxps, j, jj, jk, k, k1;
   float ci2, dnx, dny, dkx, dky, at1, at2, at3;
   float complex zero, zt1, zt2, zt3;
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
   ci2 = ci*ci;
/* calculate magnetic field and sum field energy */
   wp = 0.0;
   if (kstrt > nxh)
      goto L40;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (j = 0; j < kxps; j++) {
      dkx = dnx*(float) (j + joff);
      jj = nyhd*j;
      jk = nyv*j;
      if ((j+joff) > 0) {
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            dky = dny*(float) k;
            at1 = ci2*crealf(ffc[k+jj]);
            at2 = dky*at1;
            at3 = dkx*at1;
            at1 = at1*cimagf(ffc[k+jj]);
            zt1 = -cimagf(cu[2+3*k+3*jk])
                + crealf(cu[2+3*k+3*jk])*_Complex_I;
            zt2 = -cimagf(cu[1+3*k+3*jk])
                + crealf(cu[1+3*k+3*jk])*_Complex_I;
            zt3 = -cimagf(cu[3*k+3*jk])
                + crealf(cu[3*k+3*jk])*_Complex_I;
            bxy[3*k+3*jk] = at2*zt1;
            bxy[1+3*k+3*jk] = -at3*zt1;
            bxy[2+3*k+3*jk] = at3*zt2 - at2*zt3;
            zt1 = -cimagf(cu[2+3*k1+3*jk])
                + crealf(cu[2+3*k1+3*jk])*_Complex_I;
            zt2 = -cimagf(cu[1+3*k1+3*jk])
                + crealf(cu[1+3*k1+3*jk])*_Complex_I;
            zt3 = -cimagf(cu[3*k1+3*jk])
                + crealf(cu[3*k1+3*jk])*_Complex_I;
            bxy[3*k1+3*jk] = -at2*zt1;
            bxy[1+3*k1+3*jk] = -at3*zt1;
            bxy[2+3*k1+3*jk] = at3*zt2 + at2*zt3;
            wp += at1*(cu[3*k+3*jk]*conjf(cu[3*k+3*jk])
               + cu[1+3*k+3*jk]*conjf(cu[1+3*k+3*jk])
               + cu[2+3*k+3*jk]*conjf(cu[2+3*k+3*jk])
               + cu[3*k1+3*jk]*conjf(cu[3*k1+3*jk])
               + cu[1+3*k1+3*jk]*conjf(cu[1+3*k1+3*jk])
               + cu[2+3*k1+3*jk]*conjf(cu[2+3*k1+3*jk]));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         at1 = ci2*crealf(ffc[jj]);
         at2 = dkx*at1;
         at1 = at1*cimagf(ffc[jj]);
         zt1 = -cimagf(cu[2+3*jk])
             + crealf(cu[2+3*jk])*_Complex_I;
         zt2 = -cimagf(cu[1+3*jk])
             + crealf(cu[1+3*jk])*_Complex_I;
         bxy[3*jk] = zero;
         bxy[1+3*jk] = -at2*zt1;
         bxy[2+3*jk] = at2*zt2;
         bxy[3*k1+3*jk] = zero;
         bxy[1+3*k1+3*jk] = zero;
         bxy[2+3*k1+3*jk] = zero;
         wp += at1*(cu[3*jk]*conjf(cu[3*jk])
            + cu[1+3*jk]*conjf(cu[1+3*jk])
            + cu[2+3*jk]*conjf(cu[2+3*jk]));
      }
   }
/* mode numbers kx = 0, nx/2 */
   if (ks==0) {
      for (k = 1; k < nyh; k++) {
         k1 = ny - k;
         dky = dny*(float) k;
         at1 = ci2*crealf(ffc[k]);
         at2 = dky*at1;
         at1 = at1*cimagf(ffc[k]);
         zt1 = -cimagf(cu[2+3*k]) + crealf(cu[2+3*k])*_Complex_I;
         zt2 = -cimagf(cu[3*k]) + crealf(cu[3*k])*_Complex_I;
         bxy[3*k] = at2*zt1;
         bxy[1+3*k] = zero;
         bxy[2+3*k] = -at2*zt2;
         bxy[3*k1] = zero;
         bxy[1+3*k1] = zero;
         bxy[2+3*k1] = zero;
         wp += at1*(cu[3*k]*conjf(cu[3*k]) + cu[1+3*k]*conjf(cu[1+3*k])
            + cu[2+3*k]*conjf(cu[2+3*k]));
      }
      k1 = 3*nyh;
      bxy[0] = zero;
      bxy[1] = zero;
      bxy[2] = zero;
      bxy[k1] = zero;
      bxy[1+k1] = zero;
      bxy[2+k1] = zero;
   }
L40:
   *wm = wp*((float) nx)*((float) ny);
   return;
}

/*--------------------------------------------------------------------*/
void cppmaxwel2(float complex exy[], float complex bxy[],
                float complex cu[], float complex ffc[], float affp,
                float ci, float dt, float *wf, float *wm, int nx,
                int ny, int kstrt, int nyv, int kxp, int nyhd) {
/* this subroutine solves 2d maxwell's equation in fourier space for
   transverse electric and magnetic fields with periodic boundary
   conditions.
   input: all, output: wf, wm, exy, bxy
   approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
   where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
   the magnetic field is first updated half a step using the equations:
   bx[ky][kx] = bx[ky][kx] - .5*dt*sqrt(-1)*ky*ez(kx,ky)
   by[ky][kx] = by[ky][kx] + .5*dt*sqrt(-1)*kx*ez(kx,ky)
   bz[ky][kx] = bz[ky][kx] - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
   the electric field is then updated a whole step using the equations:
   ex[ky][kx] = ex[ky][kx] + c2*dt*sqrt(-1)*ky*bz(kx,ky)
                           - affp*dt*cux(kx,ky)*s(kx,ky)
   ey[ky][kx] = ey[ky][kx] - c2*dt*sqrt(-1)*kx*bz(kx,ky)
                           - affp*dt*cuy(kx,ky)*s(kx,ky)
   ez[ky][kx] = ez[ky][kx] + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
                           - affp*dt*cuz(kx,ky)*s(kx,ky)
   the magnetic field is finally updated the remaining half step with
   the new electric field and the previous magnetic field equations.
   where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
   and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
   j,k = fourier mode numbers, except for
   ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
   ex(ky=pi) = ey(ky=pi) = ez(ky=pi) = 0,
   ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
   and similarly for bx, by, bz.
   cu[j][k][i] = i-th component of complex current density and
   exy[j][k][i] = i-th component of complex electric field,
   bxy[j][k][i] = i-th component of complex magnetic field,
   for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
   aimag(ffc[k][j]) = finite-size particle shape factor s
   s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)
   for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
   affp = normalization constant = nx*ny/np, where np=number of particles
   ci = reciprocal of velocity of light
   dt = time interval between successive calculations
   transverse electric field energy is also calculated, using
   wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
   magnetic field energy is also calculated, using
   wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
   nx/ny = system length in x/y direction
   kxp = number of data values per block
   kstrt = starting data block number
   nyv = first dimension of field arrays, must be >= ny
   nyhd = first dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, ks, joff, kxps, j, jj, jk, k, k1;
   float dnx, dny, dth, c2, cdt, adt, anorm, dkx, dky, afdt;
   float complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9;
   double wp, ws;
   if (ci <= 0.0)
      return;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   ks = kstrt - 1;
   joff = kxp*ks;
   kxps = nxh - joff;
   kxps = 0 > kxps ? 0 : kxps;
   kxps = kxp < kxps ? kxp : kxps;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dth = 0.5*dt;
   c2 = 1.0/(ci*ci);
   cdt = c2*dt;
   adt = affp*dt;
   zero = 0.0 + 0.0*_Complex_I;
   anorm = 1.0/affp;
/* calculate magnetic field and sum field energy */
   ws = 0.0;
   wp = 0.0;
   if (kstrt > nxh)
      goto L40;
/* calculate the electromagnetic fields */
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (j = 0; j < kxps; j++) {
      dkx = dnx*(float) (j + joff);
      jj = nyhd*j;
      jk = nyv*j;
      if ((j+joff) > 0) {
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            dky = dny*(float) k;
            afdt = adt*cimagf(ffc[k+jj]);
/* update magnetic field half time step, ky > 0 */
            zt1 = -cimagf(exy[2+3*k+3*jk])
                + crealf(exy[2+3*k+3*jk])*_Complex_I;
            zt2 = -cimagf(exy[1+3*k+3*jk])
                + crealf(exy[1+3*k+3*jk])*_Complex_I;
            zt3 = -cimagf(exy[3*k+3*jk])
                + crealf(exy[3*k+3*jk])*_Complex_I;
            zt4 = bxy[3*k+3*jk] - dth*(dky*zt1);
            zt5 = bxy[1+3*k+3*jk] + dth*(dkx*zt1);
            zt6 = bxy[2+3*k+3*jk] - dth*(dkx*zt2 - dky*zt3);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exy[3*k+3*jk] + cdt*(dky*zt1) - afdt*cu[3*k+3*jk];
            zt8 = exy[1+3*k+3*jk] - cdt*(dkx*zt1) - afdt*cu[1+3*k+3*jk];
            zt9 = exy[2+3*k+3*jk] + cdt*(dkx*zt2 - dky*zt3) 
                - afdt*cu[2+3*k+3*jk];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exy[3*k+3*jk] = zt7;
            exy[1+3*k+3*jk] = zt8;
            exy[2+3*k+3*jk] = zt9;
            ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) 
               + zt9*conjf(zt9));
            zt4 -= dth*(dky*zt1);
            zt5 += dth*(dkx*zt1);
            zt6 -= dth*(dkx*zt2 - dky*zt3);
            bxy[3*k+3*jk] = zt4;
            bxy[1+3*k+3*jk] = zt5;
            bxy[2+3*k+3*jk] = zt6;
            wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
               + zt6*conjf(zt6));
/* update magnetic field half time step, ky < 0 */
            zt1 = -cimagf(exy[2+3*k1+3*jk])
                + crealf(exy[2+3*k1+3*jk])*_Complex_I;
            zt2 = -cimagf(exy[1+3*k1+3*jk])
                + crealf(exy[1+3*k1+3*jk])*_Complex_I;
            zt3 = -cimagf(exy[3*k1+3*jk])
                + crealf(exy[3*k1+3*jk])*_Complex_I;
            zt4 = bxy[3*k1+3*jk] + dth*(dky*zt1);
            zt5 = bxy[1+3*k1+3*jk] + dth*(dkx*zt1);
            zt6 = bxy[2+3*k1+3*jk] - dth*(dkx*zt2 + dky*zt3);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exy[3*k1+3*jk] - cdt*(dky*zt1) - afdt*cu[3*k1+3*jk];
            zt8 = exy[1+3*k1+3*jk] - cdt*(dkx*zt1)
                - afdt*cu[1+3*k1+3*jk];
            zt9 = exy[2+3*k1+3*jk] + cdt*(dkx*zt2 + dky*zt3)
                - afdt*cu[2+3*k1+3*jk];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exy[3*k1+3*jk] = zt7;
            exy[1+3*k1+3*jk] = zt8;
            exy[2+3*k1+3*jk] = zt9;
            ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) 
               + zt9*conjf(zt9));
            zt4 += dth*(dky*zt1);
            zt5 += dth*(dkx*zt1);
            zt6 -= dth*(dkx*zt2 + dky*zt3);
            bxy[3*k1+3*jk] = zt4;
            bxy[1+3*k1+3*jk] = zt5;
            bxy[2+3*k1+3*jk] = zt6;
            wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
               + zt6*conjf(zt6));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         afdt = adt*cimagf(ffc[jj]);
/* update magnetic field half time step */
         zt1 = -cimagf(exy[2+3*jk])
             + crealf(exy[2+3*jk])*_Complex_I;
         zt2 = -cimagf(exy[1+3*jk])
             + crealf(exy[1+3*jk])*_Complex_I;
         zt5 = bxy[1+3*jk] + dth*(dkx*zt1);
         zt6 = bxy[2+3*jk] - dth*(dkx*zt2);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt8 = exy[1+3*jk] - cdt*(dkx*zt1) - afdt*cu[1+3*jk];
         zt9 = exy[2+3*jk] + cdt*(dkx*zt2) - afdt*cu[2+3*jk];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         exy[3*jk] = zero;
         exy[1+3*jk] = zt8;
         exy[2+3*jk] = zt9;
         ws += anorm*(zt8*conjf(zt8) + zt9*conjf(zt9));
         zt5 = zt5 + dth*(dkx*zt1);
         zt6 = zt6 - dth*(dkx*zt2);
         bxy[3*jk] = zero;
         bxy[1+3*jk] = zt5;
         bxy[2+3*jk] = zt6;
         wp += anorm*(zt5*conjf(zt5) + zt6*conjf(zt6));
         bxy[3*k1+3*jk] = zero;
         bxy[1+3*k1+3*jk] = zero;
         bxy[2+3*k1+3*jk] = zero;
         exy[3*k1+3*jk] = zero;
         exy[1+3*k1+3*jk] = zero;
         exy[2+3*k1+3*jk] = zero;
      }
   }
/* mode numbers kx = 0, nx/2 */
   if (ks==0) {
      for (k = 1; k < nyh; k++) {
         k1 = ny - k;
         dky = dny*(float) k;
         afdt = adt*cimagf(ffc[k]);
/* update magnetic field half time step */
         zt1 = -cimagf(exy[2+3*k]) + crealf(exy[2+3*k])*_Complex_I;
         zt3 = -cimagf(exy[3*k]) + crealf(exy[3*k])*_Complex_I;
         zt4 = bxy[3*k] - dth*(dky*zt1);
         zt6 = bxy[2+3*k] + dth*(dky*zt3);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exy[3*k] + cdt*(dky*zt1) - afdt*cu[3*k];
         zt9 = exy[2+3*k] - cdt*(dky*zt3) - afdt*cu[2+3*k];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exy[3*k] = zt7;
         exy[1+3*k] = zero;
         exy[2+3*k] = zt9;
         ws += anorm*(zt7*conjf(zt7) + zt9*conjf(zt9));
         zt4 -= dth*(dky*zt1);
         zt6 += dth*(dky*zt3);
         bxy[3*k] = zt4;
         bxy[1+3*k] = zero;
         bxy[2+3*k] = zt6;
         wp += anorm*(zt4*conjf(zt4) + zt6*conjf(zt6));
         bxy[3*k1] = zero;
         bxy[1+3*k1] = zero;
         bxy[2+3*k1] = zero;
         exy[3*k1] = zero;
         exy[1+3*k1] = zero;
         exy[2+3*k1] = zero;
      }
      k1 = 3*nyh;
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
   }
L40:
   *wf = ws*((float) nx)*((float) ny);
   *wm = c2*wp*((float) nx)*((float) ny);
   return;
}


/*--------------------------------------------------------------------*/
void cppemfield2(float complex fxy[], float complex exy[],
                 float complex ffc[], int isign, int nx, int ny,
                 int kstrt, int nyv, int kxp, int nyhd) {
/* this subroutine either adds complex vector fields if isign > 0
   or copies complex vector fields if isign < 0
   includes additional smoothing
local data                                                 */
   int i, nxh, nyh, ks, joff, kxps, j, jj, jk, k, k1;
   float at1;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   ks = kstrt - 1;
   joff = kxp*ks;
   kxps = nxh - joff;
   kxps = 0 > kxps ? 0 : kxps;
   kxps = kxp < kxps ? kxp : kxps;
   if (kstrt > nxh) return;
/* add the fields */
   if (isign > 0) {
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
      for (j = 0; j < kxps; j++) {
         jj = nyhd*j;
         jk = nyv*j;
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            at1 = cimagf(ffc[k+jj]);
            for (i = 0; i < 3; i++) {
               fxy[i+3*k+3*jk] += exy[i+3*k+3*jk]*at1;
               fxy[i+3*k1+3*jk] += exy[i+3*k1+3*jk]*at1;
            }
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         at1 = cimagf(ffc[jj]);
         for (i = 0; i < 3; i++) {
            fxy[i+3*jk] += exy[i+3*jk]*at1;
            fxy[i+3*k1+3*jk] += exy[i+3*k1+3*jk]*at1;
         }
      }
   }
/* copy the fields */
   else if (isign < 0) {
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
      for (j = 0; j < kxps; j++) {
         jj = nyhd*j;
         jk = nyv*j;
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            at1 = cimagf(ffc[k+jj]);
            for (i = 0; i < 3; i++) {
               fxy[i+3*k+3*jk] = exy[i+3*k+3*jk]*at1;
               fxy[i+3*k1+3*jk] = exy[i+3*k1+3*jk]*at1;
            }
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         at1 = cimagf(ffc[jj]);
         for (i = 0; i < 3; i++) {
            fxy[i+3*jk] = exy[i+3*jk]*at1;
            fxy[i+3*k1+3*jk] = exy[i+3*k1+3*jk]*at1;
         }
      }
   }
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
void cppfft2r3xx(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                 int nxyhd) {
/* this subroutine performs the x part of 3 two dimensional real to
   complex fast fourier transforms and their inverses, for a subset of y,
   using complex arithmetic, for data which is distributed in blocks
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
   for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
   where N = (nx/2)*ny, and nvp = number of procs
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, an inverse fourier transform is performed
   f[m][n][0:2] = (1/nx*ny)*sum(f[k][j][0:2]*exp(-sqrt(-1)*2pi*n*j/nx)
   if isign = 1, a forward fourier transform is performed
   f[k][j][0:2] = sum(f[m][n][0:2]*exp(sqrt(-1)*2pi*n*j/nx)*
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
   f[k][j][0:2] = mode j,kk, where kk = k + kyp*(kstrt - 1)
   0 <= j < nx/2 and 0 <= kk < ny, except for
   f[k][0][0:2] = mode nx/2,kk, where ny/2+1 <= kk < ny, and
   imaginary part of f[0][0][0:2] = real part of mode nx/2,0
   on mode kstrt=0
   imaginary part of f[0][0][0:2] = real part of mode nx/2,ny/2
   on mode kstrt=(ny/2)/kyp
   written by viktor k. decyk, ucla
   parallel, RISC optimized version
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny;
   int nxy, nxhy, kypt, j, k, nrx;
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, joff;
   float ani, at1, at2;
   float complex s, t, t1, t2, t3;
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
      joff = 3*nxvh*k;
      for (j = 0; j < nxh; j++) {
         at1 = crealf(f[2+3*j+joff]);
         f[2+3*j+joff] = crealf(f[1+3*j+joff])
                       + cimagf(f[2+3*j+joff])*_Complex_I;
         at2 = cimagf(f[1+3*j+joff]);
         f[1+3*j+joff] = cimagf(f[3*j+joff]) + at1*_Complex_I;
         f[3*j+joff] = crealf(f[3*j+joff]) + at2*_Complex_I;
       }
   }
/* bit-reverse array elements in x */
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = kypi-1; k < kypt; k++) {
         joff = 3*nxvh*k;
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
               joff = 3*nxvh*i;
               t1 = s*f[3*j2+joff];
               t2 = s*f[1+3*j2+joff];
               t3 = s*f[2+3*j2+joff];
               f[3*j2+joff] = f[3*j1+joff] - t1;
               f[1+3*j2+joff] = f[1+3*j1+joff] - t2;
               f[2+3*j2+joff] = f[2+3*j1+joff] - t3;
               f[3*j1+joff] += t1;
               f[1+3*j1+joff] += t2;
               f[2+3*j1+joff] += t3;
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
         joff = 3*nxvh*k;
         for (i = 0; i < 3; i++) {
            t = conjf(f[i+3*(nxh-j)+joff]);
            s = f[i+3*j+joff] + t;
            t = (f[i+3*j+joff] - t)*t1;
            f[i+3*j+joff] = ani*(s + t);
            f[i+3*(nxh-j)+joff] = ani*conjf(s - t);
         }
      }
   }
   ani = 2.0*ani;
   for (k = kypi-1; k < kypt; k++) {
      joff = 3*nxvh*k;
      for (i = 0; i < 3; i++) {
         f[i+joff] = ani*((crealf(f[i+joff]) + cimagf(f[i+joff]))
                      + (crealf(f[i+joff]) - cimagf(f[i+joff]))*_Complex_I);
         if (nxhh > 0)
            f[i+3*nxhh+joff] = ani*conjf(f[i+3*nxhh+joff]);
      }
   }
   return;
/* forward fourier transform */
L140: kmr = nxy/nx;
/* scramble coefficients */
   for (j = 1; j < nxhh; j++) {
      t1 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
      for (k = kypi-1; k < kypt; k++) {
         joff = 3*nxvh*k;
         for (i = 0; i < 3; i++) {
            t = conjf(f[i+3*(nxh-j)+joff]);
            s = f[i+3*j+joff] + t;
            t = (f[i+3*j+joff] - t)*t1;
            f[i+3*j+joff] = s + t;
            f[i+3*(nxh-j)+joff] = conjf(s - t);
         }
      }
   }
   for (k = kypi-1; k < kypt; k++) {
      joff = 3*nxvh*k;
      for (i = 0; i < 3; i++) {
         f[i+joff] = (crealf(f[i+joff]) + cimagf(f[i+joff]))
                      + (crealf(f[i+joff]) - cimagf(f[i+joff]))*_Complex_I;
         if (nxhh > 0)
            f[i+3*nxhh+joff] = 2.0*conjf(f[i+3*nxhh+joff]);
      }
   }
   nrx = nxhy/nxh;
/* bit-reverse array elements in x */
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = kypi-1; k < kypt; k++) {
         joff = 3*nxvh*k;
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
               joff = 3*nxvh*i;
               t1 = s*f[3*j2+joff];
               t2 = s*f[1+3*j2+joff];
               t3 = s*f[2+3*j2+joff];
               f[3*j2+joff] = f[3*j1+joff] - t1;
               f[1+3*j2+joff] = f[1+3*j1+joff] - t2;
               f[2+3*j2+joff] = f[2+3*j1+joff] - t3;
               f[3*j1+joff] += t1;
               f[1+3*j1+joff] +=  t2;
               f[2+3*j1+joff] +=  t3;
            }
         }
      }
      ns = ns2;
   }
/* swap complex components */
   for (k = kypi-1; k < kypt; k++) {
      joff = 3*nxvh*k;
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
void cppfft2r3xy(float complex g[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                 int nxyhd) {
/* this subroutine performs the y part of 3 two dimensional real to
   complex fast fourier transforms and their inverses, for a subset of x,
   using complex arithmetic, for data which is distributed in blocks
   for isign = (-1,1), input: all, output: g
   for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
   for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
   where N = (nx/2)*ny, and nvp = number of procs
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, an inverse fourier transform is performed
   g[n][m][0:2] = sum(g[j][k][0:2]*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, a forward fourier transform is performed
   g[j][k][0:2] = sum(g[n][m][0:2]*exp(sqrt(-1)*2pi*m*k/ny))
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
   g[j][k][0:2] = mode jj,k, where jj = j + kxp*(kstrt - 1)
   0 <= jj < nx/2 and 0 <= k < ny, except for
   g[0][k][0:2] = mode nx/2,k, where ny/2+1 <= k < ny, and
   imaginary part of g[0][0][0:2] = real part of mode nx/2,0 and
   imaginary part of g[0][ny/2][0:2] = real part of mode nx/2,ny/2
   on node kstrt=0
   written by viktor k. decyk, ucla
   parallel, RISC optimized version
local data                                                            */
   int indx1, indx1y, nx, nxh, ny, nyh;
   int nxy, nxhy, ks, kxpt, j, k, nry;
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, koff;
   float complex s, t1, t2, t3;
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
         koff = 3*nyv*j;
         t1 = g[3*k1+koff];
         t2 = g[1+3*k1+koff];
         t3 = g[2+3*k1+koff];
         g[3*k1+koff] = g[3*k+koff];
         g[1+3*k1+koff] = g[1+3*k+koff];
         g[2+3*k1+koff] = g[2+3*k+koff];
         g[3*k+koff] = t1;
         g[1+3*k+koff] = t2;
         g[2+3*k+koff] = t3;
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
               koff = 3*nyv*i;
               t1 = s*g[3*j2+koff];
               t2 = s*g[1+3*j2+koff];
               t3 = s*g[2+3*j2+koff];
               g[3*j2+koff] = g[3*j1+koff] - t1;
               g[1+3*j2+koff] = g[1+3*j1+koff] - t2;
               g[2+3*j2+koff] = g[2+3*j1+koff] - t3;
               g[3*j1+koff] += t1;
               g[1+3*j1+koff] += t2;
               g[2+3*j1+koff] += t3;
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
         for (i = 0; i < 3; i++) {
            s = g[i+3*(ny-k)];
            g[i+3*(ny-k)] = 0.5*(cimagf(g[i+3*k] + s)
                             + crealf(g[i+3*k] - s)*_Complex_I);
            g[i+3*k] = 0.5*(crealf(g[i+3*k] + s)
                             + cimagf(g[i+3*k] - s)*_Complex_I);
         }
      }
   }
   return;
/* forward fourier transform */
/* scramble modes kx = 0, nx/2 */
L90: if (ks==0) {
      for (k = 1; k < nyh; k++) {
         if (kxpi==1) {
            for (i = 0; i < 3; i++) {
               s = cimagf(g[i+3*(ny-k)])
                   + crealf(g[i+3*(ny-k)])*_Complex_I;
               g[i+3*(ny-k)] = conjf(g[i+3*k] - s);
               g[i+3*k] += s;
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
         koff = 3*nyv*j;
         t1 = g[3*k1+koff];
         t2 = g[1+3*k1+koff];
         t3 = g[2+3*k1+koff];
         g[3*k1+koff] = g[3*k+koff];
         g[1+3*k1+koff] = g[1+3*k+koff];
         g[2+3*k1+koff] = g[2+3*k+koff];
         g[3*k+koff] = t1;
         g[1+3*k+koff] = t2;
         g[2+3*k+koff] = t3;
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
               koff = 3*nyv*i;
               t1 = s*g[3*j2+koff];
               t2 = s*g[1+3*j2+koff];
               t3 = s*g[2+3*j2+koff];
               g[3*j2+koff] = g[3*j1+koff] - t1;
               g[1+3*j2+koff] = g[1+3*j1+koff] - t2;
               g[2+3*j2+koff] = g[2+3*j1+koff] - t3;
               g[3*j1+koff] += t1;
               g[1+3*j1+koff] += t2;
               g[2+3*j1+koff] += t3;
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
void cwppfft2r3(float complex f[], float complex g[], float complex bs[],
                float complex br[], int isign, int ntpose, int mixup[],
                float complex sct[], float *ttp, int indx, int indy,
                int kstrt, int nvp, int nxvh, int nyv, int kxp, int kyp,
                int kypd, int nxhyd, int nxyhd) {
/* wrapper function for 3 2d real to complex ffts, with packed data */
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
      cppfft2r3xx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,kypd,
                  nxhyd,nxyhd);
/* transpose f array to g */
      cpwtimera(-1,ttp,&dtime);
      cppntpose(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,3,nxvh,nyv,kxp,kypd);
      cpwtimera(1,ttp,&dtime);
/* perform y fft */
      cppfft2r3xy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
                  nxhyd,nxyhd);
/* transpose g array to f */
      if (ntpose==0) {
         cpwtimera(-1,&tf,&dtime);
         cppntpose(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,3,nyv,nxvh,kypd,
                   kxp);
         cpwtimera(1,&tf,&dtime);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* transpose f array to g */
      if (ntpose==0) {
         cpwtimera(-1,&tf,&dtime);
         cppntpose(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,3,nxvh,nyv,kxp,
                   kypd);
         cpwtimera(1,&tf,&dtime);
      }
/* perform y fft */
      cppfft2r3xy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
                  nxhyd,nxyhd);
/* transpose g array to f */
      cpwtimera(-1,ttp,&dtime);
      cppntpose(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,3,nyv,nxvh,kypd,kxp);
      cpwtimera(1,ttp,&dtime);
/* perform x fft */
      cppfft2r3xx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,kypd,
                  nxhyd,nxyhd);
   }
   if (ntpose==0)
      *ttp += tf;
   return;
}

/*
c-----------------------------------------------------------------------
      subroutine WPPFFT2R(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy
     1,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv, kxp, kyp
      integer kypd, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(nxvh,kypd), g(nyv,kxp)
      dimension bs(kxp,kyp), br(kxp,kyp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi, ks, kxpp, kypp
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      ks = kstrt - 1
      kxpp = min(kxp,max(0,nxh-kxp*ks))
      kypp = min(kyp,max(0,ny-kyp*ks))
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PPFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh
     1,kypd,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp,ky
     1pd)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PPFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,
     1kxp,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,kyp
     1d,kxp)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp
     1,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PPFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,
     1kxp,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,kypd,k
     1xp)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PPFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh
     1,kypd,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPPFFT2R3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,ind
     1y,kstrt,nvp,nxvh,nyv,kxp,kyp,kypd,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv, kxp, kyp
      integer kypd, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(3,nxvh,kypd), g(3,nyv,kxp)
      dimension bs(3,kxp,kyp), br(3,kxp,kyp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi, ks, kxpp, kypp
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      ks = kstrt - 1
      kxpp = min(kxp,max(0,nxh-kxp*ks))
      kypp = min(kyp,max(0,ny-kyp*ks))
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PPFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxv
     1h,kypd,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,3,nxvh,nyv,kxp
     1,kypd)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PPFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv
     1,kxp,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,3,nyv,nxvh,
     1kypd,kxp)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPNTPOSE(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,3,nxvh,nyv,
     1kxp,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PPFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv
     1,kxp,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPNTPOSE(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,3,nyv,nxvh,kyp
     1d,kxp)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PPFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxv
     1h,kypd,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,n
     1xvh,kypd,nxhyd,nxyhd)
c this subroutine performs the x part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kypd = second dimension of f
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k) = mode j-1,kk-1, where kk = k + kyp*(kstrt - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1) = real part of mode nx/2,0 on mode kstrt=0
c imaginary part of f(1,1) = real part of mode nx/2,ny/2
c on mode kstrt=(ny/2)/kyp
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, kstrt, nxvh, kypi, kypp, kypd
      integer nxhyd, nxyhd, mixup
      complex f, sct
      dimension f(nxvh,kypd)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, kypt, j, k, nrx
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      real ani
      complex s, t, t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 100
c inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 k = kypi, kypt
      t = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t
   10 continue
   20 continue
c first transform in x
      nrx = nxy/nxh
      do 60 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = kypi, kypt
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 80 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = kypi, kypt
      t = conjg(f(nxh2-j,k))
      s = f(j,k) + t
      t = (f(j,k) - t)*t1
      f(j,k) = ani*(s + t)
      f(nxh2-j,k) = ani*conjg(s - t)
   70 continue
   80 continue
      do 90 k = kypi, kypt
      f(1,k) = 2.0*ani*cmplx(real(f(1,k)) + aimag(f(1,k)),real(f(1,k)) -
     1 aimag(f(1,k)))
      if (nxhh.gt.0) f(nxhh+1,k) = 2.0*ani*conjg(f(nxhh+1,k))
   90 continue
      return
c forward fourier transform
  100 kmr = nxy/nx
c scramble coefficients
      do 120 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 k = kypi, kypt
      t = conjg(f(nxh2-j,k))
      s = f(j,k) + t
      t = (f(j,k) - t)*t1
      f(j,k) = s + t
      f(nxh2-j,k) = conjg(s - t)
  110 continue
  120 continue
      do 130 k = kypi, kypt
      f(1,k) = cmplx(real(f(1,k)) + aimag(f(1,k)),real(f(1,k)) - aimag(f
     1(1,k)))
      if (nxhh.gt.0) f(nxhh+1,k) = 2.0*conjg(f(nxhh+1,k))
  130 continue
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 150 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 150
      do 140 k = kypi, kypt
      t = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t
  140 continue
  150 continue
c then transform in x
      nrx = nxy/nxh
      do 190 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 160 i = kypi, kypt
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,n
     1yv,kxp,nxhyd,nxyhd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c g(m,n) = sum(g(k,j)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(k,j) = sum(g(m,n)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c kxp = number of x indices per block
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g
c kxp = number of data values per block in x
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(k,j) = mode jj-1,k-1, where jj = j + kxp*(kstrt - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1) = real part of mode nx/2,ny/2
c on node kstrt=0
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxpi, kxpp, nyv
      integer kxp, nxhyd, nxyhd
      complex g, sct
      dimension g(nyv,kxp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      complex s, t
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 80
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = kxpi, kxpt
      t = g(k1,j)
      g(k1,j) = g(k,j)
      g(k,j) = t
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
      do 60 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = kxpi, kxpt
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      if (ks.gt.0) return
      do 70 k = 2, nyh
      if (kxpi.eq.1) then
         s = g(ny2-k,1)
         g(ny2-k,1) = 0.5*cmplx(aimag(g(k,1) + s),real(g(k,1) - s))
         g(k,1) = 0.5*cmplx(real(g(k,1) + s),aimag(g(k,1) - s))
      endif
   70 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
   80 nry = nxhy/ny
      if (ks.gt.0) go to 100
      do 90 k = 2, nyh
      if (kxpi.eq.1) then
         s = cmplx(aimag(g(ny2-k,1)),real(g(ny2-k,1)))
         g(ny2-k,1) = conjg(g(k,1) - s)
         g(k,1) = g(k,1) + s
      endif
   90 continue
c bit-reverse array elements in y
  100 do 120 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 120
      do 110 j = kxpi, kxpt
      t = g(k1,j)
      g(k1,j) = g(k,j)
      g(k,j) = t
  110 continue
  120 continue
c first transform in y
      nry = nxy/ny
      do 160 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 150 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 140 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 130 i = kxpi, kxpt
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
  130 continue
  140 continue
  150 continue
  160 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,
     1nxvh,kypd,nxhyd,nxyhd)
c this subroutine performs the x part of 3 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kypd = second dimension of f
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(1:3,j,k) = mode j-1,kk-1, where kk = k + kyp*(kstrt - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1:3,1,k) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1:3,1,1) = real part of mode nx/2,0
c on mode kstrt=0
c imaginary part of f(1:3,1,1) = real part of mode nx/2,ny/2
c on mode kstrt=(ny/2)/kyp
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, nxvh, kypi, kypp
      integer kypd, nxhyd, nxyhd
      complex f, sct
      dimension f(3,nxvh,kypd)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, kypt, j, k, nrx
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      real ani, at1, at2
      complex s, t, t1, t2, t3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 140
c inverse fourier transform
      ani = 0.5/(real(nx)*real(ny))
      nrx = nxhy/nxh
c swap complex components
      do 20 k = kypi, kypt
      do 10 j = 1, nxh
      at1 = real(f(3,j,k))
      f(3,j,k) = cmplx(real(f(2,j,k)),aimag(f(3,j,k)))
      at2 = aimag(f(2,j,k))
      f(2,j,k) = cmplx(aimag(f(1,j,k)),at1)
      f(1,j,k) = cmplx(real(f(1,j,k)),at2)
   10 continue
   20 continue
c bit-reverse array elements in x
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 k = kypi, kypt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      t3 = f(3,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(3,j1,k) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
   30 continue
   40 continue
c first transform in x
      nrx = nxy/nxh
      do 80 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 50 i = kypi, kypt
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      t3 = s*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(3,j2,i) = f(3,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
      f(3,j1,i) = f(3,j1,i) + t3
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 110 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = kypi, kypt
      do 90 i = 1, 3
      t = conjg(f(i,nxh2-j,k))
      s = f(i,j,k) + t
      t = (f(i,j,k) - t)*t1
      f(i,j,k) = ani*(s + t)
      f(i,nxh2-j,k) = ani*conjg(s - t)
   90 continue
  100 continue
  110 continue
      do 130 k = kypi, kypt
      do 120 i = 1, 3
      f(i,1,k) = 2.0*ani*cmplx(real(f(i,1,k)) + aimag(f(i,1,k)),real(f(i
     1,1,k)) - aimag(f(i,1,k)))
      if (nxhh.gt.0) f(i,nxhh+1,k) = 2.0*ani*conjg(f(i,nxhh+1,k))
  120 continue
  130 continue
      return
c forward fourier transform
  140 kmr = nxy/nx
c scramble coefficients
      do 170 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 160 k = kypi, kypt
      do 150 i = 1, 3
      t = conjg(f(i,nxh2-j,k))
      s = f(i,j,k) + t
      t = (f(i,j,k) - t)*t1
      f(i,j,k) = s + t
      f(i,nxh2-j,k) = conjg(s - t)
  150 continue
  160 continue
  170 continue
      do 190 k = kypi, kypt
      do 180 i = 1, 3
      f(i,1,k) = cmplx(real(f(i,1,k)) + aimag(f(i,1,k)),real(f(i,1,k)) -
     1aimag(f(i,1,k)))
      if (nxhh.gt.0) f(i,nxhh+1,k) = 2.0*conjg(f(i,nxhh+1,k))
  180 continue
  190 continue
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 210 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 210
      do 200 k = kypi, kypt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      t3 = f(3,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(3,j1,k) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  200 continue
  210 continue
c then transform in x
      nrx = nxy/nxh
      do 250 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 220 i = kypi, kypt
      t1 = s*f(1,j2,i)
      t2 = s*f(2,j2,i)
      t3 = s*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t1
      f(2,j2,i) = f(2,j1,i) - t2
      f(3,j2,i) = f(3,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t1
      f(2,j1,i) = f(2,j1,i) + t2
      f(3,j1,i) = f(3,j1,i) + t3
  220 continue
  230 continue
  240 continue
  250 continue
c swap complex components
      do 270 k = kypi, kypt
      do 260 j = 1, nxh
      at1 = real(f(3,j,k))
      f(3,j,k) = cmplx(aimag(f(2,j,k)),aimag(f(3,j,k)))
      at2 = real(f(2,j,k))
      f(2,j,k) = cmplx(at1,aimag(f(1,j,k)))
      f(1,j,k) = cmplx(real(f(1,j,k)),at2)
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,
     1nyv,kxp,nxhyd,nxyhd)
c this subroutine performs the y part of 3 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c g(1:3,m,n) = sum(g(1:3,k,j)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(1:3,k,j) = sum(g(1:3,m,n)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g
c kxp = number of data values per block in x
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(1:3,k,j) = mode jj-1,k-1, where jj = j + kxp*(kstrt - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(1:3,k,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1:3,1,1) = real part of mode nx/2,0 and
c imaginary part of g(1:3,ny/2+1,1) = real part of mode nx/2,ny/2
c on node kstrt=0
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxpi, kxpp, nyv
      integer kxp, nxhyd, nxyhd
      complex g, sct
      dimension g(3,nyv,kxp)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      complex s, t1, t2, t3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 90
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = kxpi, kxpt
      t1 = g(1,k1,j)
      t2 = g(2,k1,j)
      t3 = g(3,k1,j)
      g(1,k1,j) = g(1,k,j)
      g(2,k1,j) = g(2,k,j)
      g(3,k1,j) = g(3,k,j)
      g(1,k,j) = t1
      g(2,k,j) = t2
      g(3,k,j) = t3
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
      do 60 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = kxpi, kxpt
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      t3 = s*g(3,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(3,j2,i) = g(3,j1,i) - t3
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
      g(3,j1,i) = g(3,j1,i) + t3
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      if (ks.gt.0) return
      do 80 k = 2, nyh
      if (kxpi.eq.1) then
         do 70 i = 1, 3
         s = g(i,ny2-k,1)
         g(i,ny2-k,1) = 0.5*cmplx(aimag(g(i,k,1) + s),real(g(i,k,1) - s)
     1)
         g(i,k,1) = 0.5*cmplx(real(g(i,k,1) + s),aimag(g(i,k,1) - s))
   70    continue
      endif
   80 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
   90 nry = nxhy/ny
      if (ks.gt.0) go to 120
      do 110 k = 2, nyh
      if (kxpi.eq.1) then
         do 100 i = 1, 3
         s = cmplx(aimag(g(i,ny2-k,1)),real(g(i,ny2-k,1)))
         g(i,ny2-k,1) = conjg(g(i,k,1) - s)
         g(i,k,1) = g(i,k,1) + s
  100    continue
      endif
  110 continue
c bit-reverse array elements in y
  120 do 140 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 140
      do 130 j = kxpi, kxpt
      t1 = g(1,k1,j)
      t2 = g(2,k1,j)
      t3 = g(3,k1,j)
      g(1,k1,j) = g(1,k,j)
      g(2,k1,j) = g(2,k,j)
      g(3,k1,j) = g(3,k,j)
      g(1,k,j) = t1
      g(2,k,j) = t2
      g(3,k,j) = t3
  130 continue
  140 continue
c first transform in y
      nry = nxy/ny
      do 180 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 150 i = kxpi, kxpt
      t1 = s*g(1,j2,i)
      t2 = s*g(2,j2,i)
      t3 = s*g(3,j2,i)
      g(1,j2,i) = g(1,j1,i) - t1
      g(2,j2,i) = g(2,j1,i) - t2
      g(3,j2,i) = g(3,j1,i) - t3
      g(1,j1,i) = g(1,j1,i) + t1
      g(2,j1,i) = g(2,j1,i) + t2
      g(3,j1,i) = g(3,j1,i) + t3
  150 continue
  160 continue
  170 continue
  180 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGSPOST2L(part,q,npp,noff,qm,idimp,npmax,nxv,nxyp)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c q(j,k) = charge density at grid point (j,kk),
c where kk = k + noff - 1
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nxv = first dimension of charge array, must be >= nx+1
c nxyp = actual first dimension of charge array, must be >= nxv*nypmx
      implicit none
      integer npp, noff, idimp, npmax, nxv, nxyp
      real qm
      real part, q
      dimension part(idimp,npmax), q(nxyp)
c local data
      integer mnoff, j, nnn, mmn, nn, mm, mp
      real dxn, dyn, dxp, dyp, amx, amy, dx1
      if (npp.lt.1) return
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j)
      mmn = part(2,j)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j) - real(nnn)
      dyn = part(2,j) - real(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyp
      mmn = mmn - mnoff
c deposit charge
      dx1 = q(mp+1) + dxp*dyp
      dyp = q(mp) + amx*dyp
      dxp = q(mm+1) + dxp*amy
      amy = q(mm) + amx*amy
      q(mp+1) = dx1
      q(mp) = dyp
      q(mm+1) = dxp
      q(mm) = amy
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyn
c deposit charge
      q(mp+1) = q(mp+1) + dxp*dyn
      q(mp) = q(mp) + amx*dyn
      q(mm+1) = q(mm+1) + dxp*amy
      q(mm) = q(mm) + amx*amy
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGSJPOST2L(part,cu,edges,npp,noff,ihole,qm,dt,nx,ny,id
     1imp,npmax,nxv,nxyp,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation, and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c also determines list of particles which are leaving this processor
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 35 flops/particle, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qm*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qm*dx*(1.-dy)
c cu(i,n,m+1)=qm*(1.-dx)*dy
c cu(i,n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = x velocity of particle n in partition
c part(4,n) = y velocity of particle n in partition
c part(5,n) = z velocity of particle n in partition
c cu(i,j,k) = ith component of current density at grid point (j,kk),
c where kk = k + noff - 1
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = first dimension of current array, must be >= nx+1
c nxyp = actual first dimension of current array, must be >= nxv*nypmx
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nxyp
      integer ipbc
      real qm, dt
      real part, cu, edges
      integer ihole
      dimension part(idimp,npmax), cu(3,nxyp)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, mmn, nnn, nop, j, nn, mm, mp, ih, nh
      real edgelx, edgely, edgerx, edgery, dxn, dyn, dxp, dyp, amx, amy
      real dx, dy, dz, vx, vy, vz, dx1, dy1
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      ih = 0
      nh = 0
      if (npp.lt.1) go to 20
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j)
      mmn = part(2,j)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j) - real(nnn)
      dyn = part(2,j) - real(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyp
      mmn = mmn - mnoff
c deposit current
      dx = dxp*dyp
      dz = amx*dyp
      vx = part(3,j-1)
      vy = part(4,j-1)
      vz = part(5,j-1)
      dx1 = cu(1,mp+1) + vx*dx
      dy1 = cu(2,mp+1) + vy*dx
      dyp = cu(3,mp+1) + vz*dx
      dx = cu(1,mp) + vx*dz
      dy = cu(2,mp) + vy*dz
      dz = cu(3,mp) + vz*dz
      cu(1,mp+1) = dx1
      cu(2,mp+1) = dy1
      cu(3,mp+1) = dyp
      cu(1,mp) = dx
      cu(2,mp) = dy
      cu(3,mp) = dz
      dx = dxp*amy
      dz = amx*amy
      dxp = cu(1,mm+1) + vx*dx
      amx = cu(2,mm+1) + vy*dx
      dyp = cu(3,mm+1) + vz*dx
      dx = cu(1,mm) + vx*dz
      dy = cu(2,mm) + vy*dz
      dz = cu(3,mm) + vz*dz
      cu(1,mm+1) = dxp
      cu(2,mm+1) = amx
      cu(3,mm+1) = dyp
      cu(1,mm) = dx
      cu(2,mm) = dy
      cu(3,mm) = dz
c advance position half a time-step
      dx = part(1,j-1) + vx*dt
      dy = part(2,j-1) + vy*dt
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1)
            part(4,j-1) = -part(4,j-1)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j - 1
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j-1) = dx
      part(2,j-1) = dy
   10 continue
      nop = npp
c deposit current for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyn
c deposit current
      dx = dxp*dyn
      dy = amx*dyn
      vx = part(3,nop)
      vy = part(4,nop)
      vz = part(5,nop)
      cu(1,mp+1) = cu(1,mp+1) + vx*dx
      cu(2,mp+1) = cu(2,mp+1) + vy*dx
      cu(3,mp+1) = cu(3,mp+1) + vz*dx
      cu(1,mp) = cu(1,mp) + vx*dy
      cu(2,mp) = cu(2,mp) + vy*dy
      cu(3,mp) = cu(3,mp) + vz*dy
      dx = dxp*amy
      dy = amx*amy
      cu(1,mm+1) = cu(1,mm+1) + vx*dx
      cu(2,mm+1) = cu(2,mm+1) + vy*dx
      cu(3,mm+1) = cu(3,mm+1) + vz*dx
      cu(1,mm) = cu(1,mm) + vx*dy
      cu(2,mm) = cu(2,mm) + vy*dy
      cu(3,mm) = cu(3,mm) + vz*dy
c advance position half a time-step
      dx = part(1,nop) + vx*dt
      dy = part(2,nop) + vy*dt
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = nop
         else
            nh = 1
         endif
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c set end of file flag
   20 if (nh.gt.0) ih = -ih
      ihole(1) = ih
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGSRJPOST2L(part,cu,edges,npp,noff,ihole,qm,dt,ci,nx,n
     1y,idimp,npmax,nxv,nxyp,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles,
c and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c also determines list of particles which are leaving this processor
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 45 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qm*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qm*dx*(1.-dy)
c cu(i,n,m+1)=qm*(1.-dx)*dy
c cu(i,n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = x momentum of particle n in partition
c part(4,n) = y momentum of particle n in partition
c part(5,n) = z momentum of particle n in partition
c cu(i,j,k) = ith component of current density at grid point (j,kk),
c where kk = k + noff - 1
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = first dimension of current array, must be >= nx+1
c nxyp = actual first dimension of current array, must be >= nxv*nypmx
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nxyp
      integer ipbc
      real qm, dt, ci
      real part, cu, edges
      integer ihole
      dimension part(idimp,npmax), cu(3,nxyp)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, mmn, nnn, nop, j, nn, mm, mp, ih, nh
      real ci2, edgelx, edgely, edgerx, edgery, dxn, dyn, dxp, dyp
      real amx, amy, dx, dy, dz, vx, vy, vz, dx1, dy1, p2, gami
      real vxn, vyn, vzn
      ci2 = ci*ci
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      ih = 0
      nh = 0
      if (npp.lt.1) go to 20
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
c find inverse gamma
      vxn = part(3,1)
      vyn = part(4,1)
      vzn = part(5,1)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
      gami = 1.0/sqrt(1.0 + p2*ci2)
      mmn = mmn - mnoff
      do 10 j = 2, npp
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j)
      mmn = part(2,j)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j) - real(nnn)
      dyn = part(2,j) - real(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyp
      mmn = mmn - mnoff
c calculate weights
      dx = dxp*dyp
      dz = amx*dyp
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
c get momentum for next particle
      vxn = part(3,j)
      vyn = part(4,j)
      vzn = part(5,j)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
c deposit current
      dx1 = cu(1,mp+1) + vx*dx
      dy1 = cu(2,mp+1) + vy*dx
      dyp = cu(3,mp+1) + vz*dx
      dx = cu(1,mp) + vx*dz
      dy = cu(2,mp) + vy*dz
      dz = cu(3,mp) + vz*dz
      cu(1,mp+1) = dx1
      cu(2,mp+1) = dy1
      cu(3,mp+1) = dyp
      cu(1,mp) = dx
      cu(2,mp) = dy
      cu(3,mp) = dz
      dx = dxp*amy
      dz = amx*amy
      dxp = cu(1,mm+1) + vx*dx
      amx = cu(2,mm+1) + vy*dx
      dyp = cu(3,mm+1) + vz*dx
      dx = cu(1,mm) + vx*dz
      dy = cu(2,mm) + vy*dz
      dz = cu(3,mm) + vz*dz
      cu(1,mm+1) = dxp
      cu(2,mm+1) = amx
      cu(3,mm+1) = dyp
      cu(1,mm) = dx
      cu(2,mm) = dy
      cu(3,mm) = dz
c find inverse gamma for next particle
      gami = 1.0/sqrt(1.0 + p2*ci2)
c advance position half a time-step
      dx = part(1,j-1) + vx*dt
      dy = part(2,j-1) + vy*dt
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1)
            part(4,j-1) = -part(4,j-1)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j - 1
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j-1) = dx
      part(2,j-1) = dy
   10 continue
      nop = npp
c deposit current for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1.0 - dyn
c deposit current
      dx = dxp*dyn
      dy = amx*dyn
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
      cu(1,mp+1) = cu(1,mp+1) + vx*dx
      cu(2,mp+1) = cu(2,mp+1) + vy*dx
      cu(3,mp+1) = cu(3,mp+1) + vz*dx
      cu(1,mp) = cu(1,mp) + vx*dy
      cu(2,mp) = cu(2,mp) + vy*dy
      cu(3,mp) = cu(3,mp) + vz*dy
      dx = dxp*amy
      dy = amx*amy
      cu(1,mm+1) = cu(1,mm+1) + vx*dx
      cu(2,mm+1) = cu(2,mm+1) + vy*dx
      cu(3,mm+1) = cu(3,mm+1) + vz*dx
      cu(1,mm) = cu(1,mm) + vx*dy
      cu(2,mm) = cu(2,mm) + vy*dy
      cu(3,mm) = cu(3,mm) + vz*dy
c advance position half a time-step
      dx = part(1,nop) + vx*dt
      dy = part(2,nop) + vy*dt
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = nop
         else
            nh = 1
         endif
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c set end of file flag
   20 if (nh.gt.0) ih = -ih
      ihole(1) = ih
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGSBPUSH23L(part,fxy,bxy,edges,npp,noff,ihole,qbm,dt,d
     1tc,ek,nx,ny,idimp,npmax,nxv,nxyp,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c also determines list of particles which are leaving this processor
c 117 flops/particle, 1 divide, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = velocity vx of particle n in partition
c part(4,n) = velocity vy of particle n in partition
c part(5,n) = velocity vz of particle n in partition
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c fxy(3,j,k) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff - 1
c bxy(1,j,k) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c edges(1:2) = lower:upper boundary of particle partition
c qbm = particle charge/mass ratio
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = first dimension of field arrays, must be >= nx+1
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nxyp
      integer ipbc
      real qbm, dt, dtc, ek
      real part, fxy, bxy, edges
      integer ihole
      dimension part(idimp,npmax), fxy(3,nxyp), bxy(3,nxyp)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, mmn, nnn, nop, nop1, j, nn, mm, mp, ih, nh
      real qtmh, edgelx, edgely, edgerx, edgery, dxn, dyn
      real dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      ih = 0
      nh = 0
      if (npp.lt.1) go to 20
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      mmn = mmn - mnoff
      nop1 = npp - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1)
      mmn = part(2,j+1)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1) - real(nnn)
      dyn = part(2,j+1) - real(mmn)
      mm = mm + nn
      amx = 1.0 - dxp
      mp = mm + nxv
      amy = 1.0 - dyp
      mmn = mmn - mnoff
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1) 
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1) 
     1+ amx*fxy(2,mm))
      dz = dyp*(dxp*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxp*fxy(3,mm+1) 
     1+ amx*fxy(3,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxp*bxy(1,mm+1) 
     1+ amx*bxy(1,mm))
      oy = dyp*(dxp*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxp*bxy(2,mm+1) 
     1+ amx*bxy(2,mm))
      oz = dyp*(dxp*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxp*bxy(3,mm+1) 
     1+ amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1. + omt)
      omt = 0.5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      nop = npp
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1.0 - dxn
      mp = mm + nxv
      amy = 1.0 - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1) 
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1) 
     1+ amx*fxy(2,mm))
      dz = dyn*(dxn*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxn*fxy(3,mm+1) 
     1+ amx*fxy(3,mm))
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxn*bxy(1,mm+1) 
     1+ amx*bxy(1,mm))
      oy = dyn*(dxn*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxn*bxy(2,mm+1) 
     1+ amx*bxy(2,mm))
      oz = dyn*(dxn*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxn*bxy(3,mm+1) 
     1+ amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1. + omt)
      omt = 0.5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,nop) = dx
      part(4,nop) = dy
      part(5,nop) = dz
c new position
      dx = part(1,nop) + dx*dtc
      dy = part(2,nop) + dy*dtc
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = nop
         else
            nh = 1
         endif
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c set end of file flag
   20 if (nh.gt.0) ih = -ih
      ihole(1) = ih
c normalize kinetic energy
      ek = ek + 0.5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PPGSRBPUSH23L(part,fxy,bxy,edges,npp,noff,ihole,qbm,dt,
     1dtc,ci,ek,nx,ny,idimp,npmax,nxv,nxyp,idps,ntmax,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c also determines list of particles which are leaving this processor
c 129 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n in partition
c part(2,n) = position y of particle n in partition
c part(3,n) = momentum px of particle n in partition
c part(4,n) = momentum py of particle n in partition
c part(5,n) = momentum pz of particle n in partition
c fxy(1,j,k) = x component of force/charge at grid (j,kk)
c fxy(2,j,k) = y component of force/charge at grid (j,kk)
c fxy(3,j,k) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff - 1
c bxy(1,j,k) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c edges(1:2) = lower:upper boundary of particle partition
c npp = number of particles in partition
c noff = lowermost global gridpoint in particle partition.
c ihole = location of hole left in particle arrays
c ihole(1) = ih, number of holes left (error, if negative)
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nxv = first dimension of field arrays, must be >= nx+1
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, noff, nx, ny, idimp, npmax, idps, ntmax, nxv, nxyp
      integer ipbc
      real qbm, dt, dtc, ci, ek
      real part, fxy, bxy, edges
      integer ihole
      dimension part(idimp,npmax), fxy(3,nxyp), bxy(3,nxyp)
      dimension edges(idps), ihole(ntmax+1)
c local data
      integer mnoff, mmn, nnn, nop, nop1, j, nn, mm, mp, ih, nh
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxn, dyn
      real dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = .5*qbm*dt
      ci2 = ci*ci
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 1.0
      edgerx = real(nx)
      edgery = real(ny-1)
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      ih = 0
      nh = 0
      if (npp.lt.1) go to 20
      mnoff = noff
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmn)
      mmn = mmn - mnoff
      nop1 = npp - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1)
      mmn = part(2,j+1)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1) - real(nnn)
      dyn = part(2,j+1) - real(mmn)
      mm = mm + nn
      amx = 1.0 - dxp
      mp = mm + nxv
      amy = 1.0 - dyp
      mmn = mmn - mnoff
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1) 
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1) 
     1+ amx*fxy(2,mm))
      dz = dyp*(dxp*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxp*fxy(3,mm+1) 
     1+ amx*fxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxp*bxy(1,mm+1) 
     1+ amx*bxy(1,mm))
      oy = dyp*(dxp*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxp*bxy(2,mm+1) 
     1+ amx*bxy(2,mm))
      oz = dyp*(dxp*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxp*bxy(3,mm+1) 
     1+ amx*bxy(3,mm))
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1. + omt)
      omt = 0.5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = j
         else
            nh = 1
         endif
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      nop = npp
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1.0 - dxn
      mp = mm + nxv
      amy = 1.0 - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1) 
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1) 
     1+ amx*fxy(2,mm))
      dz = dyn*(dxn*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxn*fxy(3,mm+1) 
     1+ amx*fxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxn*bxy(1,mm+1) 
     1+ amx*bxy(1,mm))
      oy = dyn*(dxn*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxn*bxy(2,mm+1) 
     1+ amx*bxy(2,mm))
      oz = dyn*(dxn*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxn*bxy(3,mm+1) 
     1+ amx*bxy(3,mm))
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1. + omt)
      omt = 0.5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,nop) = dx
      part(4,nop) = dy
      part(5,nop) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop) + dx*dtg
      dy = part(2,nop) + dy*dtg
c periodic boundary conditions in x
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
      endif
c find particles out of bounds
      if ((dy.lt.edges(1)).or.(dy.ge.edges(2))) then
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(ih+1) = nop
         else
            nh = 1
         endif
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c set end of file flag
   20 if (nh.gt.0) ih = -ih
      ihole(1) = ih
c normalize kinetic energy
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      function ranorm()
c this program calculates a random number y from a gaussian distribution
c with zero mean and unit variance, according to the method of
c mueller and box:
c    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
c    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
c where x is a random number uniformly distributed on (0,1).
c written for the ibm by viktor k. decyk, ucla
      implicit none
      integer iflg,isc,i1,r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end
c-----------------------------------------------------------------------
      function randum()
c this is a version of the random number generator dprandom due to
c c. bingham and the yale computer center, producing numbers
c in the interval (0,1).  written for the sun by viktor k. decyk, ucla
      implicit none
      integer isc,i1,r1,r2
      double precision randum,h1l,h1u,r0,r3,asc,bsc
      save r1,r2,h1l,h1u
      data r1,r2 /1271199957,1013501921/
      data h1l,h1u /65533.0d0,32767.0d0/
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      randum = (dble(r1) + dble(r2)*asc)*asc
      return
      end
*/

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void cpdicomp2l_(float *edges, int *nyp, int *noff, int *nypmx,
                 int *nypmn, int *ny, int *kstrt, int *nvp, int *idps) {
   cpdicomp2l(edges,nyp,noff,nypmx,nypmn,*ny,*kstrt,*nvp,*idps);
   return;
}

/*--------------------------------------------------------------------*/
void cpdistr2h_(float *part, float *edges, int *npp, int *nps,
                float *vtx, float *vty, float *vtz, float *vdx,
                float *vdy, float *vdz, int *npx, int *npy, int *nx,
                int *ny, int *idimp, int *npmax, int *idps, int *ipbc,
                int *ierr) {
   cpdistr2h(part,edges,npp,*nps,*vtx,*vty,*vtz,*vdx,*vdy,*vdz,*npx,
             *npy,*nx,*ny,*idimp,*npmax,*idps,*ipbc,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void cppgbpush23l_(float *part, float *fxy, float *bxy, float *edges,
                   int *npp, int *noff, int *ihole, float *qbm,
                   float *dt, float *dtc, float *ek, int *nx, int *ny,
                   int *idimp, int *npmax, int *nxv, int *nypmx,
                   int *idps, int *ntmax, int *ipbc) {
   cppgbpush23l(part,fxy,bxy,edges,*npp,*noff,ihole,*qbm,*dt,*dtc,ek,
                *nx,*ny,*idimp,*npmax,*nxv,*nypmx,*idps,*ntmax,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdppgbpush23l_(double *part, double *fxy, double *bxy,
                    double *edges, int *npp, int *noff, int *ihole,
                    double *qbm, double *dt, double *dtc, double *ek,
                    int *nx, int *ny, int *idimp, int *npmax, int *nxv,
                    int *nypmx,int *idps, int *ntmax, int *ipbc) {
   cdppgbpush23l(part,fxy,bxy,edges,*npp,*noff,ihole,*qbm,*dt,*dtc,ek,
                 *nx,*ny,*idimp,*npmax,*nxv,*nypmx,*idps,*ntmax,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrbpush23l_(float *part, float *fxy, float *bxy, float *edges,
                    int *npp, int *noff, int *ihole, float *qbm,
                    float *dt, float *dtc, float *ci, float *ek,
                    int *nx, int *ny, int *idimp, int *npmax, int *nxv,
                    int *nypmx, int *idps, int *ntmax, int *ipbc) {
   cppgrbpush23l(part,fxy,bxy,edges,*npp,*noff,ihole,*qbm,*dt,*dtc,*ci,
                 ek,*nx,*ny,*idimp,*npmax,*nxv,*nypmx,*idps,*ntmax,
                 *ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdppgrbpush23l_(double *part, double *fxy, double *bxy,
                     double *edges, int *npp, int *noff, int *ihole,
                     double *qbm, double *dt, double *dtc, double *ci, 
                     double *ek, int *nx, int *ny, int *idimp,
                     int *npmax, int *nxv, int *nypmx, int *idps,
                     int *ntmax, int *ipbc) {
   cdppgrbpush23l(part,fxy,bxy,edges,*npp,*noff,ihole,*qbm,*dt,*dtc,*ci,
                  ek,*nx,*ny,*idimp,*npmax,*nxv,*nypmx,*idps,*ntmax,
                  *ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgpost2l_(float *part, float *q, int *npp, int *noff, float *qm,
                 int *idimp, int *npmax, int *nxv, int *nypmx) {
   cppgpost2l(part,q,*npp,*noff,*qm,*idimp,*npmax,*nxv,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppgjpost2l_(float *part, float *cu, float *edges, int *npp,
                  int *noff, int *ihole, float *qm, float *dt, int *nx,
                  int *ny, int *idimp, int *npmax, int *nxv, int *nypmx,
                  int *idps, int *ntmax, int *ipbc) {
   cppgjpost2l(part,cu,edges,*npp,*noff,ihole,*qm,*dt,*nx,*ny,*idimp,
               *npmax,*nxv,*nypmx,*idps,*ntmax,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrjpost2l_(float *part, float *cu, float *edges, int *npp,
                   int *noff, int *ihole, float *qm, float *dt,
                   float *ci, int *nx, int *ny, int *idimp, int *npmax,
                   int *nxv, int *nypmx, int *idps, int *ntmax,
                   int *ipbc) {
   cppgrjpost2l(part,cu,edges,*npp,*noff,ihole,*qm,*dt,*ci,*nx,*ny,
                *idimp,*npmax,*nxv,*nypmx,*idps,*ntmax,*ipbc);
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
void cppcguard2xl_(float *fxy, int *myp, int *nx, int *ndim, int *nxe,
                   int *nypmx) {
   cppcguard2xl(fxy,*myp,*nx,*ndim,*nxe,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppaguard2xl_(float *q, int *myp, int *nx, int *nxe, int *nypmx) {
   cppaguard2xl(q,*myp,*nx,*nxe,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppacguard2xl_(float *cu, int *myp, int *nx, int *ndim, int *nxe,
                    int *nypmx) {
   cppacguard2xl(cu,*myp,*nx,*ndim,*nxe,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppois23_(float complex *q, float complex *fxy, int *isign,
               float complex *ffc, float *ax, float *ay, float *affp,
               float *we, int *nx, int *ny, int *kstrt, int *nyv,
               int *kxp, int *nyhd) {
   cppois23(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*kstrt,*nyv,*kxp,
            *nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppcuperp2_(float complex *cu, int *nx, int *ny, int *kstrt,
                 int *nyv, int *kxp) {
   cppcuperp2(cu,*nx,*ny,*kstrt,*nyv,*kxp);
   return;
}

/*--------------------------------------------------------------------*/
void cippbpoisp23_(float complex *cu, float complex *bxy,
                   float complex *ffc, float *ci, float *wm, int *nx,
                   int *ny, int *kstrt, int *nyv, int *kxp, int *nyhd) {
   cippbpoisp23(cu,bxy,ffc,*ci,wm,*nx,*ny,*kstrt,*nyv,*kxp,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppmaxwel2_(float complex *exy, float complex *bxy,
                 float complex *cu, float complex *ffc, float *affp,
                 float *ci, float *dt, float *wf, float *wm, int *nx,
                 int *ny, int *kstrt, int *nyv, int *kxp, int *nyhd) {
   cppmaxwel2(exy,bxy,cu,ffc,*affp,*ci,*dt,wf,wm,*nx,*ny,*kstrt,*nyv,
              *kxp,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppemfield2_(float complex *fxy, float complex *exy,
                  float complex *ffc, int *isign, int *nx, int *ny,
                  int *kstrt, int *nyv, int *kxp, int *nyhd) {
   cppemfield2(fxy,exy,ffc,*isign,*nx,*ny,*kstrt,*nyv,*kxp,*nyhd);
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
void cppfft2r3xx_(float complex *f, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *kstrt,
                  int *kypi, int *kypp, int *nxvh, int *kypd, int *nxhyd,
                  int *nxyhd) {
   cppfft2r3xx(f,*isign,mixup,sct,*indx,*indy,*kstrt,*kypi,*kypp,*nxvh,
               *kypd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2r3xy_(float complex *g, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *kstrt,
                  int *kxpi, int *kxpp, int *nyv, int *kxp, int *nxhyd,
                  int *nxyhd) {
   cppfft2r3xy(g,*isign,mixup,sct,*indx,*indy,*kstrt,*kxpi,*kxpp,*nyv,
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
void cwppfft2r3_(float complex *f, float complex *g, float complex *bs,
                 float complex *br, int *isign, int *ntpose, int *mixup,
                 float complex *sct, float *ttp, int *indx, int *indy,
                 int *kstrt, int *nvp, int *nxvh, int *nyv, int *kxp,
                 int *kyp, int *kypd, int *nxhyd, int *nxyhd) {
   cwppfft2r3(f,g,bs,br,*isign,*ntpose,mixup,sct,ttp,*indx,*indy,*kstrt,
              *nvp,*nxvh,*nyv,*kxp,*kyp,*kypd,*nxhyd,*nxyhd);
   return;
}
