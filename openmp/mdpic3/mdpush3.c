/* C Library for Skeleton 3D Darwin OpenMP PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "mdpush3.h"

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
void cppmovin3l(float part[], float ppart[], int kpic[], int nppmx,
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
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
   ppart[m][n][2] = position z of particle n in tile m
   ppart[m][n][3] = velocity vx of particle n in tile m
   ppart[m][n][4] = velocity vy of particle n in tile m
   ppart[m][n][5] = velocity vz of particle n in tile m
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
void cppcheck3l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int ny, int nz, int mx, int my, int mz, int mx1,
                int my1, int mz1, int *irc) {
/* this subroutine performs a sanity check to make sure particles sorted
   by x,y,z grid in tiles of mx, my, mz, are all within bounds.
   tiles are assumed to be arranged in 3D linear memory
   input: all except irc
   output: irc
   ppart[l][n][0] = position x of particle n in tile l
   ppart[l][n][1] = position y of particle n in tile l
   ppart[l][n][2] = position a of particle n in tile l
   kpic(l) = number of reordered output particles in tile l
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
         dx = ppart[idimp*(j+nppmx*l)];
         dy = ppart[1+idimp*(j+nppmx*l)];
         dz = ppart[2+idimp*(j+nppmx*l)];
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
void cgbppush3l(float ppart[], float fxyz[], float bxyz[], int kpic[],
                float qbm, float dt, float dtc, float *ek, int idimp,
                int nppmx, int nx, int ny, int nz, int mx, int my,
                int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                int mxyz1, int ipbc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field.  Using the Boris Mover.
   OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   190 flops/particle, 1 divide, 54 loads, 6 stores
   input: all, output: ppart, ek
   velocity equations used are:
   vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fx(x(t),y(t),z(t))*dt)
   vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fy(x(t),y(t),z(t))*dt)
   vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
   omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
   omz = (q/m)*bz(x(t),y(t),z(t)).
   position equations used are:
   x(t+dt)=x(t) + vx(t+dt/2)*dt
   y(t+dt)=y(t) + vy(t+dt/2)*dt
   z(t+dt)=z(t) + vz(t+dt/2)*dt
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
   bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
   ppart[m][n][2] = position z of particle n in tile m
   ppart[m][n][3] = velocity vx of particle n in tile m
   ppart[m][n][4] = velocity vy of particle n in tile m
   ppart[m][n][5] = velocity vz of particle n in tile m
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   bxyz[l][k][j][0] = x component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][1] = y component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][2] = z component of magnetic field at grid (j,k,l)
   that is, the convolution of magnetic field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass ratio
   dt = time interval between successive force calculations
   dtc = time interval between successive co-ordinate calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
        .25*(vz(t+dt/2) + vz(t-dt/2))**2)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
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
   int mxy1, noff, moff, loff, npoff, npp;
   int i, j, k, l, nn, mm, ll, nm, mxv, myv, mxyv, nxyv;
   float qtmh, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, z;
   float sfxyz[3*MXV*MYV*MZV], sbxyz[3*MXV*MYV*MZV];
/* float sfxyz[3*(mx+1)*(my+1)*(mz+1)]; */
/* float sbxyz[3*(mx+1)*(my+1)*(mz+1)]; */
   double sum1, sum2;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   mxy1 = mx1*my1;
   qtmh = 0.5f*qbm*dt;
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
private(i,j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,nm,x,y,z,dxp,dyp,dzp, \
amx,amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm, \
rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sum1,sfxyz,sbxyz) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = nppmx*l;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      ll = (mz < nz-loff ? mz : nz-loff) + 1;
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
            for (i = 0; i < nn; i++) {
               sfxyz[3*(i+mxv*j+mxyv*k)]
               = fxyz[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+3*(i+mxv*j+mxyv*k)]
               = fxyz[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+3*(i+mxv*j+mxyv*k)]
               = fxyz[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
            for (i = 0; i < nn; i++) {
               sbxyz[3*(i+mxv*j+mxyv*k)]
               = bxyz[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[1+3*(i+mxv*j+mxyv*k)]
               = bxyz[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[2+3*(i+mxv*j+mxyv*k)]
               = bxyz[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
      sum1 = 0.0;
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         z = ppart[2+idimp*(j+npoff)];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nm = 3*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find electric field */
         nn = nm;
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+3];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+3];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+3];
         mm = nn + 3*mxv;
         dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+3]);
         dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+3]);
         dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+3]);
         nn += 3*mxyv;
         acx = amx*sfxyz[nn] + amy*sfxyz[nn+3];
         acy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+3];
         acz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+3];
         mm = nn + 3*mxv;
         dx = dx + dzp*(acx + dyp*sfxyz[mm] + dx1*sfxyz[mm+3]);
         dy = dy + dzp*(acy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+3]);
         dz = dz + dzp*(acz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+3]);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxyz[nn] + amy*sbxyz[nn+3];
         oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+3];
         oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+3];
         mm = nn + 3*mxv;
         ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+3]);
         oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+3]);
         oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+3]);
         nn += 3*mxyv;
         acx = amx*sbxyz[nn] + amy*sbxyz[nn+3];
         acy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+3];
         acz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+3];
         mm = nn + 3*mxv;
         ox = ox + dzp*(acx + dyp*sbxyz[mm] + dx1*sbxyz[mm+3]);
         oy = oy + dzp*(acy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+3]);
         oz = oz + dzp*(acz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+3]);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[3+idimp*(j+npoff)] + dx;
         acy = ppart[4+idimp*(j+npoff)] + dy;
         acz = ppart[5+idimp*(j+npoff)] + dz;
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
         dx += (rot1*acx + rot2*acy + rot3*acz)*anorm;
         dy += (rot4*acx + rot5*acy + rot6*acz)*anorm;
         dz += (rot7*acx + rot8*acy + rot9*acz)*anorm;
         ppart[3+idimp*(j+npoff)] = dx;
         ppart[4+idimp*(j+npoff)] = dy;
         ppart[5+idimp*(j+npoff)] = dz;
/* new position */
         dx = x + dx*dtc;
         dy = y + dy*dtc;
         dz = z + dz*dtc;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[3+idimp*(j+npoff)] = -ppart[3+idimp*(j+npoff)];
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[4+idimp*(j+npoff)] = -ppart[4+idimp*(j+npoff)];
            }
            if ((dz < edgelz) || (dz >= edgerz)) {
               dz = z;
               ppart[5+idimp*(j+npoff)] = -ppart[5+idimp*(j+npoff)];
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[3+idimp*(j+npoff)] = -ppart[3+idimp*(j+npoff)];
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[4+idimp*(j+npoff)] = -ppart[4+idimp*(j+npoff)];
            }
         }
/* set new position */
         ppart[idimp*(j+npoff)] = dx;
         ppart[1+idimp*(j+npoff)] = dy;
         ppart[2+idimp*(j+npoff)] = dz;
      }
      sum2 += sum1;
   }
/* normalize kinetic energy */
   *ek += 0.5f*sum2;
   return;
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void cgbppushf3l(float ppart[], float fxyz[], float bxyz[], int kpic[],
                 int ncl[], int ihole[], float qbm, float dt, float dtc,
                 float *ek, int idimp, int nppmx, int nx, int ny,
                 int nz, int mx, int my, int mz, int nxv, int nyv,
                 int nzv, int mx1, int my1, int mxyz1, int ntmax,
                 int *irc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field.  Using the Boris Mover.
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells
   data read in tiles
   particles stored segmented array
   190 flops/particle, 1 divide, 54 loads, 6 stores
   input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
   velocity equations used are:
   vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fx(x(t),y(t),z(t))*dt)
   vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fy(x(t),y(t),z(t))*dt)
   vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
     .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
   omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
   omz = (q/m)*bz(x(t),y(t),z(t)).
   position equations used are:
   x(t+dt)=x(t) + vx(t+dt/2)*dt
   y(t+dt)=y(t) + vy(t+dt/2)*dt
   z(t+dt)=z(t) + vz(t+dt/2)*dt
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
   bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
   ppart[m][n][2] = position z of particle n in tile m
   ppart[m][n][3] = velocity vx of particle n in tile m
   ppart[m][n][4] = velocity vy of particle n in tile m
   ppart[m][n][5] = velocity vz of particle n in tile m
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   bxyz[l][k][j][0] = x component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][1] = y component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][2] = z component of magnetic field at grid (j,k,l)
   that is, the convolution of magnetic field over particle shape
   kpic[l] = number of particles in tile l
   ncl[l][i] = number of particles going to destination i, tile l
   ihole[l][:][0] = location of hole in array left by departing particle
   ihole[l][:][1] = direction destination of particle leaving hole
   all for tile l
   ihole[l][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive force calculations
   dtc = time interval between successive co-ordinate calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
        .25*(vz(t+dt/2) + vz(t-dt/2))**2)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
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
   int mxy1, noff, moff, loff, npoff, npp;
   int i, j, k, l, ih, nh, nn, mm, ll, nm, mxv, myv, mxyv, nxyv;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1;
   float qtmh, acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, z;
   float sfxyz[3*MXV*MYV*MZV], sbxyz[3*MXV*MYV*MZV];
/* float sfxyz[3*(mx+1)*(my+1)*(mz+1)]; */
/* float sbxyz[3*(mx+1)*(my+1)*(mz+1)]; */
   double sum1, sum2;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   mxy1 = mx1*my1;
   qtmh = 0.5f*qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
   sum2 = 0.0;
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,nm,ih,nh,x,y,z,dxp, \
dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt, \
omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely, \
edgelz,edgerx,edgery,edgerz,sum1,sfxyz,sbxyz) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = nppmx*l;
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
               sfxyz[3*(i+mxv*j+mxyv*k)]
               = fxyz[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+3*(i+mxv*j+mxyv*k)]
               = fxyz[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+3*(i+mxv*j+mxyv*k)]
               = fxyz[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
            for (i = 0; i < nn; i++) {
               sbxyz[3*(i+mxv*j+mxyv*k)]
               = bxyz[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[1+3*(i+mxv*j+mxyv*k)]
               = bxyz[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[2+3*(i+mxv*j+mxyv*k)]
               = bxyz[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
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
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         z = ppart[2+idimp*(j+npoff)];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nm = 3*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find electric field */
         nn = nm;
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+3];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+3];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+3];
         mm = nn + 3*mxv;
         dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+3]);
         dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+3]);
         dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+3]);
         nn += 3*mxyv;
         acx = amx*sfxyz[nn] + amy*sfxyz[nn+3];
         acy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+3];
         acz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+3];
         mm = nn + 3*mxv;
         dx = dx + dzp*(acx + dyp*sfxyz[mm] + dx1*sfxyz[mm+3]);
         dy = dy + dzp*(acy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+3]);
         dz = dz + dzp*(acz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+3]);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxyz[nn] + amy*sbxyz[nn+3];
         oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+3];
         oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+3];
         mm = nn + 3*mxv;
         ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+3]);
         oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+3]);
         oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+3]);
         nn += 3*mxyv;
         acx = amx*sbxyz[nn] + amy*sbxyz[nn+3];
         acy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+3];
         acz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+3];
         mm = nn + 3*mxv;
         ox = ox + dzp*(acx + dyp*sbxyz[mm] + dx1*sbxyz[mm+3]);
         oy = oy + dzp*(acy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+3]);
         oz = oz + dzp*(acz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+3]);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[3+idimp*(j+npoff)] + dx;
         acy = ppart[4+idimp*(j+npoff)] + dy;
         acz = ppart[5+idimp*(j+npoff)] + dz;
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
         dx += (rot1*acx + rot2*acy + rot3*acz)*anorm;
         dy += (rot4*acx + rot5*acy + rot6*acz)*anorm;
         dz += (rot7*acx + rot8*acy + rot9*acz)*anorm;
         ppart[3+idimp*(j+npoff)] = dx;
         ppart[4+idimp*(j+npoff)] = dy;
         ppart[5+idimp*(j+npoff)] = dz;
/* new position */
         dx = x + dx*dtc;
         dy = y + dy*dtc;
         dz = z + dz*dtc;
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
         ppart[idimp*(j+npoff)] = dx;
         ppart[1+idimp*(j+npoff)] = dy;
         ppart[2+idimp*(j+npoff)] = dz;
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
   *ek += 0.5f*sum2;
   return;
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void cgppost3l(float ppart[], float q[], int kpic[], float qm,
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
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
   ppart[m][n][2] = position z of particle n in tile m
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
   float x, y, z, dxp, dyp, dzp, amx, amy, amz, dx1;
   float sq[MXV*MYV*MZV];
/* float sq[(mx+1)*(my+1)*(mz+1)]; */
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   mxy1 = mx1*my1;
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
#pragma omp parallel for \
private(i,j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,nm,lm,x,y,z,dxp,dyp, \
dzp,amx,amy,amz,dx1,sq)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = nppmx*l;
/* zero out local accumulator */
      for (j = 0; j < mxyv*(mz+1); j++) {
         sq[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         z = ppart[2+idimp*(j+npoff)];
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
         sq[nn] = x;
         sq[nn+1] = y;
         mm = nn + mxv;
         x = sq[mm] + dyp*amz;
         y = sq[mm+1] + dx1*amz;
         sq[mm] = x;
         sq[mm+1] = y;
         nn += mxyv;
         x = sq[nn] + amx*dzp;
         y = sq[nn+1] + amy*dzp;
         sq[nn] = x;
         sq[nn+1] = y;
         mm = nn + mxv;
         x = sq[mm] + dyp*dzp;
         y = sq[mm+1] + dx1*dzp;
         sq[mm] = x;
         sq[mm+1] = y;
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
void cgjppost3l(float ppart[], float cu[], int kpic[], float qm,
                float dt, int nppmx, int idimp, int nx, int ny, int nz,
                int mx, int my, int mz, int nxv, int nyv, int nzv,
                int mx1, int my1, int mxyz1, int ipbc) {
/* for 3d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   69 flops/particle, 30 loads, 27 stores
   input: all, output: ppart, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
   cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
   cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
   cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
   cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
   cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
   cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
   cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   and qci = qm*vi, where i = x,y,z
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
   ppart[m][n][2] = position z of particle n in tile m
   ppart[m][n][3] = velocity vx of particle n in tile m
   ppart[m][n][4] = velocity vy of particle n in tile m
   ppart[m][n][5] = velocity vz of particle n in tile m
   cu[l][k][j][i] = ith component of current density at grid point j,k,l
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 6
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
   nzv = fourth dimension of current array, must be >= nz+1
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
   int mxy1, noff, moff, loff, npoff, npp;
   int i, j, k, l, nn, mm, ll, nm, lm, mxv, myv, mxyv, nxyv;
   float edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz;
   float x, y, z;
   float scu[3*MXV*MYV*MZV];
/* float scu[3*(mx+1)*(my+1)*(mz+1)]; */
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   mxy1 = mx1*my1;
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
private(i,j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,nm,lm,x,y,z,dxp,dyp, \
dzp,amx,amy,amz,dx1,dx,dy,dz,vx,vy,vz,scu)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = nppmx*l;
/* zero out local accumulator */
      for (j = 0; j < 3*mxyv*(mz+1); j++) {
         scu[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         z = ppart[2+idimp*(j+npoff)];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = 3*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* deposit current within tile to local accumulator */
         dx = amx*amz;
         dy = amy*amz;
         vx = ppart[3+idimp*(j+npoff)];
         vy = ppart[4+idimp*(j+npoff)];
         vz = ppart[5+idimp*(j+npoff)];
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*amz;
         scu[nn+3] += vx*dy;
         scu[nn+1+3] += vy*dy;
         scu[nn+2+3] += vz*dy;
         dy = dx1*amz;
         mm = nn + 3*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         dx = amx*dzp;
         scu[mm+3] += vx*dy;
         scu[mm+1+3] += vy*dy;
         scu[mm+2+3] += vz*dy;
         dy = amy*dzp;
         nn += 3*mxyv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*dzp;
         scu[nn+3] += vx*dy;
         scu[nn+1+3] += vy*dy;
         scu[nn+2+3] += vz*dy;
         dy = dx1*dzp;
         mm = nn + 3*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         scu[mm+3] += vx*dy;
         scu[mm+1+3] += vy*dy;
         scu[mm+2+3] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[3+idimp*(j+npoff)] = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[4+idimp*(j+npoff)] = -vy;
            }
            if ((dz < edgelz) || (dz >= edgerz)) {
               dz = z;
               ppart[5+idimp*(j+npoff)] = -vz;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[3+idimp*(j+npoff)] = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[4+idimp*(j+npoff)] = -vy;
            }
         }
/* set new position */
         ppart[idimp*(j+npoff)] = dx;
         ppart[1+idimp*(j+npoff)] = dy;
         ppart[2+idimp*(j+npoff)] = dz;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
            for (i = 1; i < nn; i++) {
               cu[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[3*(i+mxv*j+mxyv*k)];
               cu[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+3*(i+mxv*j+mxyv*k)];
               cu[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+3*(i+mxv*j+mxyv*k)];
            }
         }
      }
/* deposit current to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[3*(i+noff+nxv*(j+moff)+nxyv*loff)] += scu[3*(i+mxv*j)];
#pragma omp atomic
            cu[1+3*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[1+3*(i+mxv*j)];
#pragma omp atomic
            cu[2+3*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[2+3*(i+mxv*j)];
            if (lm > mz) {
#pragma omp atomic
               cu[3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[3*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+3*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+3*(i+mxv*j+mxyv*(lm-1))];
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
            cu[3*(i+noff+nxv*moff+nxyv*(k+loff))] += scu[3*(i+mxyv*k)];
#pragma omp atomic
            cu[1+3*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[1+3*(i+mxyv*k)];
#pragma omp atomic
            cu[2+3*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[2+3*(i+mxyv*k)];
            if (mm > my) {
#pragma omp atomic
               cu[3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[3*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[1+3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[1+3*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[2+3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[2+3*(i+mxv*(mm-1)+mxyv*k)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[3*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[1+3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[1+3*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[2+3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[2+3*(mxv*j+mxyv*k)];
            if (nm > mx) {
#pragma omp atomic
               cu[3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[3*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[1+3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+3*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[2+3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+3*(nm-1+mxv*j+mxyv*k)];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[3*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[1+3*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[2+3*(i+mxyv*(lm-1))];
            if (mm > my) {
#pragma omp atomic
               cu[3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[3*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[1+3*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[2+3*(i+mxv*(mm-1)+mxyv*(lm-1))];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[3*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[1+3*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[2+3*(mxv*j+mxyv*(lm-1))];
            if (nm > mx) {
#pragma omp atomic
               cu[3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[3*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+3*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+3*(nm-1+mxv*j+mxyv*(lm-1))];
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
void cgmjppost3l(float ppart[], float amu[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int my, int mz, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1) {
/* for 3d code, this subroutine calculates particle momentum flux
   using first-order linear interpolation
   OpenMP version using guard cells
   data deposited in tiles
   particles stored segmented array
   121 flops/particle, 52 loads, 48 stores
   input: all, output: ppart, amu
   momentum flux is approximated by values at the nearest grid points
   amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
   amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
   amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
   amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
   amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
   amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
   amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
   amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   and qci = qm*vj*vk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
   where vj = vj(t-dt/2) and vk = vk(t-dt/2)
   ppart[m][n][0] = position x of particle n in tile m at t
   ppart[m][n][1] = position y of particle n in tile m at t
   ppart[m][n][2] = position z of particle n in tile m at t
   ppart[m][n][3] = velocity vx of particle n in tile m at t - dt/2
   ppart[m][n][4] = velocity vy of particle n in tile m at t - dt/2
   ppart[m][n][5] = velocity vz of particle n in tile m at t - dt/2
   amu[l][k][j][i] = ith component of momentum flux at grid point j,k,l
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 6
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
   nzv = fourth dimension of current array, must be >= nz+1
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
   float dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, vx, vy, vz;
   float x, y, z, v1, v2, v3, v4, v5, v6;
   float samu[6*MXV*MYV*MZV];
/* float samu[6*(mx+1)*(my+1)*(mz+1)]; */
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   mxy1 = mx1*my1;
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,nm,lm,x,y,z,dxp,dyp, \
dzp,amx,amy,amz,dx1,dx,dy,vx,vy,vz,v1,v2,v3,v4,v5,v6,samu)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = nppmx*l;
/* zero out local accumulator */
      for (j = 0; j < 6*mxyv*(mz+1); j++) {
         samu[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         z = ppart[2+idimp*(j+npoff)];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = 6*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* deposit momentum flux within tile to local accumulator */
         dx = amx*amz;
         dy = amy*amz;
         vx = ppart[3+idimp*(j+npoff)];
         vy = ppart[4+idimp*(j+npoff)];
         vz = ppart[5+idimp*(j+npoff)];
         v1 = vx*vx;
         v2 = vx*vy;
         v3 = vx*vz;
         v4 = vy*vy;
         v5 = vy*vz;
         v6 = vz*vz;
         samu[nn] += v1*dx;
         samu[nn+1] += v2*dx;
         samu[nn+2] += v3*dx;
         samu[nn+3] += v4*dx;
         samu[nn+4] += v5*dx;
         samu[nn+5] += v6*dx;
         dx = dyp*amz;
         samu[nn+6] += v1*dy;
         samu[nn+1+6] += v2*dy;
         samu[nn+2+6] += v3*dy;
         samu[nn+3+6] += v4*dy;
         samu[nn+4+6] += v5*dy;
         samu[nn+5+6] += v6*dy;
         dy = dx1*amz;
         mm = nn + 6*mxv;
         samu[mm] += v1*dx;
         samu[mm+1] += v2*dx;
         samu[mm+2] += v3*dx;
         samu[mm+3] += v4*dx;
         samu[mm+4] += v5*dx;
         samu[mm+5] += v6*dx;
         dx = amx*dzp;
         samu[mm+6] += v1*dy;
         samu[mm+1+6] += v2*dy;
         samu[mm+2+6] += v3*dy;
         samu[mm+3+6] += v4*dy;
         samu[mm+4+6] += v5*dy;
         samu[mm+5+6] += v6*dy;
         dy = amy*dzp;
         nn += 6*mxyv;
         samu[nn] += v1*dx;
         samu[nn+1] += v2*dx;
         samu[nn+2] += v3*dx;
         samu[nn+3] += v4*dx;
         samu[nn+4] += v5*dx;
         samu[nn+5] += v6*dx;
         dx = dyp*dzp;
         samu[nn+6] += v1*dy;
         samu[nn+1+6] += v2*dy;
         samu[nn+2+6] += v3*dy;
         samu[nn+3+6] += v4*dy;
         samu[nn+4+6] += v5*dy;
         samu[nn+5+6] += v6*dy;
         dy = dx1*dzp;
         mm = nn + 6*mxv;
         samu[mm] += v1*dx;
         samu[mm+1] += v2*dx;
         samu[mm+2] += v3*dx;
         samu[mm+3] += v4*dx;
         samu[mm+4] += v5*dx;
         samu[mm+5] += v6*dx;
         samu[mm+6] += v1*dy;
         samu[mm+1+6] += v2*dy;
         samu[mm+2+6] += v3*dy;
         samu[mm+3+6] += v4*dy;
         samu[mm+4+6] += v5*dy;
         samu[mm+5+6] += v6*dy;
      }
/* deposit momentum flux to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
            for (i = 1; i < nn; i++) {
               amu[6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[6*(i+mxv*j+mxyv*k)];
               amu[1+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[1+6*(i+mxv*j+mxyv*k)];
               amu[2+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[2+6*(i+mxv*j+mxyv*k)];
               amu[3+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[3+6*(i+mxv*j+mxyv*k)];
               amu[4+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[4+6*(i+mxv*j+mxyv*k)];
               amu[5+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[5+6*(i+mxv*j+mxyv*k)];
            }
         }
      }
/* deposit momentum flux to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            amu[6*(i+noff+nxv*(j+moff)+nxyv*loff)] += samu[6*(i+mxv*j)];
#pragma omp atomic
            amu[1+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[1+6*(i+mxv*j)];
#pragma omp atomic
            amu[2+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[2+6*(i+mxv*j)];
#pragma omp atomic
            amu[3+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[3+6*(i+mxv*j)];
#pragma omp atomic
            amu[4+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[4+6*(i+mxv*j)];
#pragma omp atomic
            amu[5+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[5+6*(i+mxv*j)];
            if (lm > mz) {
#pragma omp atomic
               amu[6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[1+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[1+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[2+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[2+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[3+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[3+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[4+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[4+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[5+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[5+6*(i+mxv*j+mxyv*(lm-1))];
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
            amu[6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[6*(i+mxyv*k)];
#pragma omp atomic
            amu[1+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[1+6*(i+mxyv*k)];
#pragma omp atomic
            amu[2+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[2+6*(i+mxyv*k)];
#pragma omp atomic
            amu[3+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[3+6*(i+mxyv*k)];
#pragma omp atomic
            amu[4+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[4+6*(i+mxyv*k)];
#pragma omp atomic
            amu[5+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[5+6*(i+mxyv*k)];
            if (mm > my) {
#pragma omp atomic
               amu[6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[1+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[1+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[2+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[2+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[3+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[3+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[4+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[4+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[5+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[5+6*(i+mxv*(mm-1)+mxyv*k)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            amu[6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[1+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[1+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[2+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[2+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[3+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[3+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[4+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[4+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[5+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[5+6*(mxv*j+mxyv*k)];
            if (nm > mx) {
#pragma omp atomic
               amu[6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[1+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[1+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[2+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[2+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[3+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[3+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[4+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[4+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[5+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[5+6*(nm-1+mxv*j+mxyv*k)];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            amu[6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[1+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[1+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[2+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[2+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[3+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[3+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[4+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[4+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[5+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[5+6*(i+mxyv*(lm-1))];
            if (mm > my) {
#pragma omp atomic
               amu[6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[1+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[1+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[2+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[2+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[3+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[3+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[4+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[4+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[5+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[5+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            amu[6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[1+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[1+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[2+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[2+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[3+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[3+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[4+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[4+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[5+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[5+6*(mxv*j+mxyv*(lm-1))];
            if (nm > mx) {
#pragma omp atomic
               amu[6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[1+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[1+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[2+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[2+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[3+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[3+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[4+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[4+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[5+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[5+6*(nm-1+mxv*j+mxyv*(lm-1))];
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
void cgdjppost3l(float ppart[], float fxyz[], float bxyz[], int kpic[],
                 float dcu[], float amu[], float qm, float qbm,
                 float dt, int idimp, int nppmx, int nx, int ny, int nz,
                 int mx, int my, int mz, int nxv, int nyv, int nzv,
                 int mx1, int my1, int mxyz1) {
/* for 3 code, this subroutine calculates particle momentum flux,
   and acceleration density using first-order spline interpolation.
   OpenMP version using guard cells
   data read/written in tiles
   particles stored segmented array
   350 flops/particle, 1 divide, 126 loads, 72 stores
   acceleration density is approximated by values at the nearest grid
   points
   dcu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
   dcu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
   dcu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
   dcu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
   dcu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
   dcu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
   dcu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
   dcu(i,n+1,m+1,l+1)=qci*dx*dy*dz
   and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
   where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
   momentum flux is approximated by values at the nearest grid points
   amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
   amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
   amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
   amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
   amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
   amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
   amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
   amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
   and qci = qm*vj*vk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
   where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
   and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   velocity equations at t=t+dt/2 are calculated from:
   vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fx(x(t),y(t),z(t))*dt)
   vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fy(x(t),y(t),z(t))*dt)
   vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
   omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
   omz = (q/m)*bz(x(t),y(t),z(t)).
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
   bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
   ppart[m][n][0] = position x of particle n in tile m at t
   ppart[m][n][1] = position y of particle n in tile m at t
   ppart[m][n][2] = position z of particle n in tile m at t
   ppart[m][n][3] = velocity vx of particle n in tile m at t - dt/2
   ppart[m][n][4] = velocity vy of particle n in tile m at t - dt/2
   ppart[m][n][5] = velocity vz of particle n in tile m at t - dt/2
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   bxyz[l][k][j][0] = x component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][1] = y component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][2] = z component of magnetic field at grid (j,k,l)
   that is, the convolution of magnetic field over particle shape
   dcu[l][k][j][i] = ith component of acceleration density
   at grid point j,k for i = 1, 3
   amu[l][k][j][i] = ith component of momentum flux
   at grid point j,k,l for i = 1, 6
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   qbm = particle charge/mass ratio
   dt = time interval between successive force calculations
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   nzv = fourth dimension of field array, must be >= nz+1
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
   float qtmh, dti, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz;
   float ox, oy, oz, dx1, acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, z, vx, vy, vz, v1, v2, v3, v4, v5, v6;
   float sfxyz[3*MXV*MYV*MZV], sbxyz[3*MXV*MYV*MZV];
   float sdcu[3*MXV*MYV*MZV], samu[6*MXV*MYV*MZV];
/* float sfxyz[3*(mx+1)*(my+1)*(mz+1)]; */
/* float sbxyz[3*(mx+1)*(my+1)*(mz+1)]; */
/* float sdcu[3*(mx+1)*(my+1)*(mz+1)]; */
/* float samu[6*(mx+1)*(my+1)*(mz+1)]; */
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   mxy1 = mx1*my1;
   qtmh = 0.5f*qbm*dt;
   dti = 1.0f/dt;
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,nm,lm,x,y,z,vx,vy,vz, \
v1,v2,v3,v4,v5,v6,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx, \
acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7, \
rot8,rot9,sfxyz,sbxyz,sdcu,samu)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = nppmx*l;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      ll = (mz < nz-loff ? mz : nz-loff) + 1;
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
            for (i = 0; i < nn; i++) {
               sfxyz[3*(i+mxv*j+mxyv*k)]
               = fxyz[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+3*(i+mxv*j+mxyv*k)]
               = fxyz[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+3*(i+mxv*j+mxyv*k)]
               = fxyz[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
            for (i = 0; i < nn; i++) {
               sbxyz[3*(i+mxv*j+mxyv*k)]
               = bxyz[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[1+3*(i+mxv*j+mxyv*k)]
               = bxyz[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[2+3*(i+mxv*j+mxyv*k)]
               = bxyz[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
/* zero out local accumulators */
      for (j = 0; j < 3*mxyv*(mz+1); j++) {
         sdcu[j] = 0.0f;
      }
      for (j = 0; j < 6*mxyv*(mz+1); j++) {
         samu[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         z = ppart[2+idimp*(j+npoff)];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nm = 3*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find electric field */
         nn = nm;
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+3];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+3];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+3];
         mm = nn + 3*mxv;
         dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+3]);
         dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+3]);
         dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+3]);
         nn += 3*mxyv;
         acx = amx*sfxyz[nn] + amy*sfxyz[nn+3];
         acy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+3];
         acz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+3];
         mm = nn + 3*mxv;
         dx = dx + dzp*(acx + dyp*sfxyz[mm] + dx1*sfxyz[mm+3]);
         dy = dy + dzp*(acy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+3]);
         dz = dz + dzp*(acz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+3]);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxyz[nn] + amy*sbxyz[nn+3];
         oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+3];
         oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+3];
         mm = nn + 3*mxv;
         ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+3]);
         oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+3]);
         oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+3]);
         nn += 3*mxyv;
         acx = amx*sbxyz[nn] + amy*sbxyz[nn+3];
         acy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+3];
         acz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+3];
         mm = nn + 3*mxv;
         ox = ox + dzp*(acx + dyp*sbxyz[mm] + dx1*sbxyz[mm+3]);
         oy = oy + dzp*(acy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+3]);
         oz = oz + dzp*(acz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+3]);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         vx = ppart[3+idimp*(j+npoff)];
         vy = ppart[4+idimp*(j+npoff)];
         vz = ppart[5+idimp*(j+npoff)];
         acx = vx + dx;
         acy = vy + dy;
         acz = vz + dz;
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
         dx += (rot1*acx + rot2*acy + rot3*acz)*anorm;
         dy += (rot4*acx + rot5*acy + rot6*acz)*anorm;
         dz += (rot7*acx + rot8*acy + rot9*acz)*anorm;
/* deposit momentum flux and acceleration density */
         nn = nm;
         amz = qm*amz;
         dzp = qm*dzp;
         ox = 0.5*(dx + vx);
         oy = 0.5*(dy + vy);
         oz = 0.5*(dz + vz);
         vx = dti*(dx - vx);
         vy = dti*(dy - vy);
         vz = dti*(dz - vz);
         dx = amx*amz;
         dy = amy*amz;
         v1 = ox*ox;
         v2 = ox*oy;
         v3 = ox*oz;
         v4 = oy*oy;
         v5 = oy*oz;
         v6 = oz*oz;
         samu[2*nn] += v1*dx;
         samu[2*nn+1] += v2*dx;
         samu[2*nn+2] += v3*dx;
         samu[2*nn+3] += v4*dx;
         samu[2*nn+4] += v5*dx;
         samu[2*nn+5] += v6*dx;
         sdcu[nn] += vx*dx;
         sdcu[nn+1] += vy*dx;
         sdcu[nn+2] += vz*dx;
         dx = dyp*amz;
         samu[2*nn+6] += v1*dy;
         samu[2*nn+1+6] += v2*dy;
         samu[2*nn+2+6] += v3*dy;
         samu[2*nn+3+6] += v4*dy;
         samu[2*nn+4+6] += v5*dy;
         samu[2*nn+5+6] += v6*dy;
         sdcu[nn+3] += vx*dy;
         sdcu[nn+1+3] += vy*dy;
         sdcu[nn+2+3] += vz*dy;
         dy = dx1*amz;
         mm = nn + 3*mxv;
         samu[2*mm] += v1*dx;
         samu[2*mm+1] += v2*dx;
         samu[2*mm+2] += v3*dx;
         samu[2*mm+3] += v4*dx;
         samu[2*mm+4] += v5*dx;
         samu[2*mm+5] += v6*dx;
         sdcu[mm] += vx*dx;
         sdcu[mm+1] += vy*dx;
         sdcu[mm+2] += vz*dx;
         dx = amx*dzp;
         samu[2*mm+6] += v1*dy;
         samu[2*mm+1+6] += v2*dy;
         samu[2*mm+2+6] += v3*dy;
         samu[2*mm+3+6] += v4*dy;
         samu[2*mm+4+6] += v5*dy;
         samu[2*mm+5+6] += v6*dy;
         sdcu[mm+3] += vx*dy;
         sdcu[mm+1+3] += vy*dy;
         sdcu[mm+2+3] += vz*dy;
         dy = amy*dzp;
         nn += 3*mxyv;
         samu[2*nn] += v1*dx;
         samu[2*nn+1] += v2*dx;
         samu[2*nn+2] += v3*dx;
         samu[2*nn+3] += v4*dx;
         samu[2*nn+4] += v5*dx;
         samu[2*nn+5] += v6*dx;
         sdcu[nn] += vx*dx;
         sdcu[nn+1] += vy*dx;
         sdcu[nn+2] += vz*dx;
         dx = dyp*dzp;
         samu[2*nn+6] += v1*dy;
         samu[2*nn+1+6] += v2*dy;
         samu[2*nn+2+6] += v3*dy;
         samu[2*nn+3+6] += v4*dy;
         samu[2*nn+4+6] += v5*dy;
         samu[2*nn+5+6] += v6*dy;
         sdcu[nn+3] += vx*dy;
         sdcu[nn+1+3] += vy*dy;
         sdcu[nn+2+3] += vz*dy;
         dy = dx1*dzp;
         mm = nn + 3*mxv;
         samu[2*mm] += v1*dx;
         samu[2*mm+1] += v2*dx;
         samu[2*mm+2] += v3*dx;
         samu[2*mm+3] += v4*dx;
         samu[2*mm+4] += v5*dx;
         samu[2*mm+5] += v6*dx;
         sdcu[mm] += vx*dx;
         sdcu[mm+1] += vy*dx;
         sdcu[mm+2] += vz*dx;
         samu[2*mm+6] += v1*dy;
         samu[2*mm+1+6] += v2*dy;
         samu[2*mm+2+6] += v3*dy;
         samu[2*mm+3+6] += v4*dy;
         samu[2*mm+4+6] += v5*dy;
         samu[2*mm+5+6] += v6*dy;
         sdcu[mm+3] += vx*dy;
         sdcu[mm+1+3] += vy*dy;
         sdcu[mm+2+3] += vz*dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
            for (i = 1; i < nn; i++) {
               amu[6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[6*(i+mxv*j+mxyv*k)];
               amu[1+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[1+6*(i+mxv*j+mxyv*k)];
               amu[2+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[2+6*(i+mxv*j+mxyv*k)];
               amu[3+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[3+6*(i+mxv*j+mxyv*k)];
               amu[4+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[4+6*(i+mxv*j+mxyv*k)];
               amu[5+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[5+6*(i+mxv*j+mxyv*k)];
               dcu[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[3*(i+mxv*j+mxyv*k)];
               dcu[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[1+3*(i+mxv*j+mxyv*k)];
               dcu[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[2+3*(i+mxv*j+mxyv*k)];
            }
         }
      }
/* deposit current to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            amu[6*(i+noff+nxv*(j+moff)+nxyv*loff)] += samu[6*(i+mxv*j)];
#pragma omp atomic
            amu[1+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[1+6*(i+mxv*j)];
#pragma omp atomic
            amu[2+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[2+6*(i+mxv*j)];
#pragma omp atomic
            amu[3+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[3+6*(i+mxv*j)];
#pragma omp atomic
            amu[4+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[4+6*(i+mxv*j)];
#pragma omp atomic
            amu[5+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[5+6*(i+mxv*j)];
#pragma omp atomic
            dcu[3*(i+noff+nxv*(j+moff)+nxyv*loff)] += sdcu[3*(i+mxv*j)];
#pragma omp atomic
            dcu[1+3*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += sdcu[1+3*(i+mxv*j)];
#pragma omp atomic
            dcu[2+3*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += sdcu[2+3*(i+mxv*j)];
            if (lm > mz) {
#pragma omp atomic
               amu[6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[1+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[1+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[2+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[2+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[3+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[3+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[4+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[4+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[5+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[5+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[3*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[1+3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[1+3*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[2+3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[2+3*(i+mxv*j+mxyv*(lm-1))];
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
            amu[6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[6*(i+mxyv*k)];
#pragma omp atomic
            amu[1+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[1+6*(i+mxyv*k)];
#pragma omp atomic
            amu[2+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[2+6*(i+mxyv*k)];
#pragma omp atomic
            amu[3+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[3+6*(i+mxyv*k)];
#pragma omp atomic
            amu[4+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[4+6*(i+mxyv*k)];
#pragma omp atomic
            amu[5+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[5+6*(i+mxyv*k)];
#pragma omp atomic
            dcu[3*(i+noff+nxv*moff+nxyv*(k+loff))]
            += sdcu[3*(i+mxyv*k)];
#pragma omp atomic
            dcu[1+3*(i+noff+nxv*moff+nxyv*(k+loff))]
            += sdcu[1+3*(i+mxyv*k)];
#pragma omp atomic
            dcu[2+3*(i+noff+nxv*moff+nxyv*(k+loff))]
            += sdcu[2+3*(i+mxyv*k)];
            if (mm > my) {
#pragma omp atomic
               amu[6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[1+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[1+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[2+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[2+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[3+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[3+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[4+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[4+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[5+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[5+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               dcu[3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += sdcu[3*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               dcu[1+3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += sdcu[1+3*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               dcu[2+3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += sdcu[2+3*(i+mxv*(mm-1)+mxyv*k)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            amu[6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[1+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[1+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[2+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[2+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[3+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[3+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[4+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[4+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[5+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[5+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            dcu[3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += sdcu[3*(mxv*j+mxyv*k)];
#pragma omp atomic
            dcu[1+3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += sdcu[1+3*(mxv*j+mxyv*k)];
#pragma omp atomic
            dcu[2+3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += sdcu[2+3*(mxv*j+mxyv*k)];
            if (nm > mx) {
#pragma omp atomic
               amu[6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[1+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[1+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[2+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[2+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[3+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[3+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[4+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[4+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[5+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[5+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               dcu[3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[3*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               dcu[1+3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[1+3*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               dcu[2+3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[2+3*(nm-1+mxv*j+mxyv*k)];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            amu[6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[1+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[1+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[2+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[2+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[3+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[3+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[4+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[4+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[5+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[5+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            dcu[3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += sdcu[3*(i+mxyv*(lm-1))];
#pragma omp atomic
            dcu[1+3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += sdcu[1+3*(i+mxyv*(lm-1))];
#pragma omp atomic
            dcu[2+3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += sdcu[2+3*(i+mxyv*(lm-1))];
            if (mm > my) {
#pragma omp atomic
               amu[6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[1+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[1+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[2+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[2+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[3+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[3+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[4+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[4+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[5+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[5+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               dcu[3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += sdcu[3*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               dcu[1+3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += sdcu[1+3*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               dcu[2+3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += sdcu[2+3*(i+mxv*(mm-1)+mxyv*(lm-1))];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            amu[6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[1+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[1+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[2+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[2+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[3+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[3+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[4+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[4+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[5+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[5+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            dcu[3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += sdcu[3*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            dcu[1+3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += sdcu[1+3*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            dcu[2+3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += sdcu[2+3*(mxv*j+mxyv*(lm-1))];
            if (nm > mx) {
#pragma omp atomic
               amu[6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[1+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[1+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[2+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[2+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[3+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[3+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[4+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[4+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[5+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[5+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[3*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[1+3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[1+3*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[2+3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[2+3*(nm-1+mxv*j+mxyv*(lm-1))];
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
void cgdcjppost3l(float ppart[], float fxyz[], float bxyz[], int kpic[],
                  float cu[], float dcu[], float amu[], float qm,
                  float qbm, float dt, int idimp, int nppmx, int nx,
                  int ny, int nz, int mx, int my, int mz, int nxv,
                  int nyv, int nzv, int mx1, int my1, int mxyz1) {
/* for 3 code, this subroutine calculates particle momentum flux,
   acceleration density and current density using first-order spline
   interpolation.
   OpenMP version using guard cells
   data read/written in tiles
   particles stored segmented array
   398 flops/particle, 1 divide, 150 loads, 96 stores
   input: all, output: cu, dcu, amu
   current density is approximated by values at the nearest grid points
   cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
   cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
   cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
   cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
   cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
   cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
   cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
   cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
   and qci = qm*vj, where j = x,y,z, for i = 1, 3
   where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
   acceleration density is approximated by values at the nearest grid
   points
   dcu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
   dcu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
   dcu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
   dcu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
   dcu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
   dcu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
   dcu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
   dcu(i,n+1,m+1,l+1)=qci*dx*dy*dz
   and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
   where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
   momentum flux is approximated by values at the nearest grid points
   amu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
   amu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
   amu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
   amu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
   amu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
   amu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
   amu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
   amu(i,n+1,m+1,l+1)=qci*dx*dy*dz
   and qci = qm*vj*vk, where jk = xx,xy,xz,yy,yz,zz, for i = 1, 6
   where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
   and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   velocity equations at t=t+dt/2 are calculated from:
   vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fx(x(t),y(t),z(t))*dt)
   vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fy(x(t),y(t),z(t))*dt)
   vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
   omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
   omz = (q/m)*bz(x(t),y(t),z(t)).
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
   bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
   ppart[m][n][0] = position x of particle n in tile m at t
   ppart[m][n][1] = position y of particle n in tile m at t
   ppart[m][n][2] = position z of particle n in tile m at t
   ppart[m][n][3] = velocity vx of particle n in tile m at t - dt/2
   ppart[m][n][4] = velocity vy of particle n in tile m at t - dt/2
   ppart[m][n][5] = velocity vz of particle n in tile m at t - dt/2
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   bxyz[l][k][j][0] = x component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][1] = y component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][2] = z component of magnetic field at grid (j,k,l)
   that is, the convolution of magnetic field over particle shape
   cu[l][k][j][i] = ith component of current density
   at grid point j,k for i = 1, 3
   dcu[l][k][j][i] = ith component of acceleration density
   at grid point j,k for i = 1, 3
   amu[l][k][j][i] = ith component of momentum flux
   at grid point j,k,l for i = 1, 6
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   qbm = particle charge/mass ratio
   dt = time interval between successive force calculations
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   nx/ny/nz = system length in x/y/z direction
   mx/my/mz = number of grids in sorting cell in x/y/z
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   nzv = fourth dimension of field array, must be >= nz+1
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
   float qtmh, dti, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz;
   float ox, oy, oz, dx1, acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, z, vx, vy, vz, v1, v2, v3, v4, v5, v6;
   float sfxyz[3*MXV*MYV*MZV], sbxyz[3*MXV*MYV*MZV];
   float scu[3*MXV*MYV*MZV], sdcu[3*MXV*MYV*MZV];
   float samu[6*MXV*MYV*MZV];
/* float sfxyz[3*(mx+1)*(my+1)*(mz+1)]; */
/* float sbxyz[3*(mx+1)*(my+1)*(mz+1)]; */
/* float scu[3*(mx+1)*(my+1)*(mz+1)]; */
/* float sdcu[3*(mx+1)*(my+1)*(mz+1)]; */
/* float samu[6*(mx+1)*(my+1)*(mz+1)]; */
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   mxy1 = mx1*my1;
   qtmh = 0.5f*qbm*dt;
   dti = 1.0f/dt;
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,noff,moff,loff,npp,npoff,nn,mm,ll,nm,lm,x,y,z,vx,vy,vz, \
v1,v2,v3,v4,v5,v6,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx, \
acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7, \
rot8,rot9,sfxyz,sbxyz,scu,sdcu,samu)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = nppmx*l;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      ll = (mz < nz-loff ? mz : nz-loff) + 1;
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
            for (i = 0; i < nn; i++) {
               sfxyz[3*(i+mxv*j+mxyv*k)]
               = fxyz[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+3*(i+mxv*j+mxyv*k)]
               = fxyz[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+3*(i+mxv*j+mxyv*k)]
               = fxyz[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
            for (i = 0; i < nn; i++) {
               sbxyz[3*(i+mxv*j+mxyv*k)]
               = bxyz[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[1+3*(i+mxv*j+mxyv*k)]
               = bxyz[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[2+3*(i+mxv*j+mxyv*k)]
               = bxyz[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
/* zero out local accumulators */
      for (j = 0; j < 3*mxyv*(mz+1); j++) {
         scu[j] = 0.0f;
         sdcu[j] = 0.0f;
      }
      for (j = 0; j < 6*mxyv*(mz+1); j++) {
         samu[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < npp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         z = ppart[2+idimp*(j+npoff)];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nm = 3*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find electric field */
         nn = nm;
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+3];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+3];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+3];
         mm = nn + 3*mxv;
         dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+3]);
         dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+3]);
         dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+3]);
         nn += 3*mxyv;
         acx = amx*sfxyz[nn] + amy*sfxyz[nn+3];
         acy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+3];
         acz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+3];
         mm = nn + 3*mxv;
         dx = dx + dzp*(acx + dyp*sfxyz[mm] + dx1*sfxyz[mm+3]);
         dy = dy + dzp*(acy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+3]);
         dz = dz + dzp*(acz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+3]);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxyz[nn] + amy*sbxyz[nn+3];
         oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+3];
         oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+3];
         mm = nn + 3*mxv;
         ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+3]);
         oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+3]);
         oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+3]);
         nn += 3*mxyv;
         acx = amx*sbxyz[nn] + amy*sbxyz[nn+3];
         acy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+3];
         acz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+3];
         mm = nn + 3*mxv;
         ox = ox + dzp*(acx + dyp*sbxyz[mm] + dx1*sbxyz[mm+3]);
         oy = oy + dzp*(acy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+3]);
         oz = oz + dzp*(acz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+3]);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         vx = ppart[3+idimp*(j+npoff)];
         vy = ppart[4+idimp*(j+npoff)];
         vz = ppart[5+idimp*(j+npoff)];
         acx = vx + dx;
         acy = vy + dy;
         acz = vz + dz;
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
         dx += (rot1*acx + rot2*acy + rot3*acz)*anorm;
         dy += (rot4*acx + rot5*acy + rot6*acz)*anorm;
         dz += (rot7*acx + rot8*acy + rot9*acz)*anorm;
/* deposit momentum flux, acceleration density, and current density */
         nn = nm;
         amz = qm*amz;
         dzp = qm*dzp;
         ox = 0.5*(dx + vx);
         oy = 0.5*(dy + vy);
         oz = 0.5*(dz + vz);
         vx = dti*(dx - vx);
         vy = dti*(dy - vy);
         vz = dti*(dz - vz);
         dx = amx*amz;
         dy = amy*amz;
         v1 = ox*ox;
         v2 = ox*oy;
         v3 = ox*oz;
         v4 = oy*oy;
         v5 = oy*oz;
         v6 = oz*oz;
         samu[2*nn] += v1*dx;
         samu[2*nn+1] += v2*dx;
         samu[2*nn+2] += v3*dx;
         samu[2*nn+3] += v4*dx;
         samu[2*nn+4] += v5*dx;
         samu[2*nn+5] += v6*dx;
         sdcu[nn] += vx*dx;
         sdcu[nn+1] += vy*dx;
         sdcu[nn+2] += vz*dx;
         scu[nn] += ox*dx;
         scu[nn+1] += oy*dx;
         scu[nn+2] += oz*dx;
         dx = dyp*amz;
         samu[2*nn+6] += v1*dy;
         samu[2*nn+1+6] += v2*dy;
         samu[2*nn+2+6] += v3*dy;
         samu[2*nn+3+6] += v4*dy;
         samu[2*nn+4+6] += v5*dy;
         samu[2*nn+5+6] += v6*dy;
         sdcu[nn+3] += vx*dy;
         sdcu[nn+1+3] += vy*dy;
         sdcu[nn+2+3] += vz*dy;
         scu[nn+3] += ox*dy;
         scu[nn+1+3] += oy*dy;
         scu[nn+2+3] += oz*dy;
         dy = dx1*amz;
         mm = nn + 3*mxv;
         samu[2*mm] += v1*dx;
         samu[2*mm+1] += v2*dx;
         samu[2*mm+2] += v3*dx;
         samu[2*mm+3] += v4*dx;
         samu[2*mm+4] += v5*dx;
         samu[2*mm+5] += v6*dx;
         sdcu[mm] += vx*dx;
         sdcu[mm+1] += vy*dx;
         sdcu[mm+2] += vz*dx;
         scu[mm] += ox*dx;
         scu[mm+1] += oy*dx;
         scu[mm+2] += oz*dx;
         dx = amx*dzp;
         samu[2*mm+6] += v1*dy;
         samu[2*mm+1+6] += v2*dy;
         samu[2*mm+2+6] += v3*dy;
         samu[2*mm+3+6] += v4*dy;
         samu[2*mm+4+6] += v5*dy;
         samu[2*mm+5+6] += v6*dy;
         sdcu[mm+3] += vx*dy;
         sdcu[mm+1+3] += vy*dy;
         sdcu[mm+2+3] += vz*dy;
         scu[mm+3] += ox*dy;
         scu[mm+1+3] += oy*dy;
         scu[mm+2+3] += oz*dy;
         dy = amy*dzp;
         nn += 3*mxyv;
         samu[2*nn] += v1*dx;
         samu[2*nn+1] += v2*dx;
         samu[2*nn+2] += v3*dx;
         samu[2*nn+3] += v4*dx;
         samu[2*nn+4] += v5*dx;
         samu[2*nn+5] += v6*dx;
         sdcu[nn] += vx*dx;
         sdcu[nn+1] += vy*dx;
         sdcu[nn+2] += vz*dx;
         scu[nn] += ox*dx;
         scu[nn+1] += oy*dx;
         scu[nn+2] += oz*dx;
         dx = dyp*dzp;
         samu[2*nn+6] += v1*dy;
         samu[2*nn+1+6] += v2*dy;
         samu[2*nn+2+6] += v3*dy;
         samu[2*nn+3+6] += v4*dy;
         samu[2*nn+4+6] += v5*dy;
         samu[2*nn+5+6] += v6*dy;
         sdcu[nn+3] += vx*dy;
         sdcu[nn+1+3] += vy*dy;
         sdcu[nn+2+3] += vz*dy;
         scu[nn+3] += ox*dy;
         scu[nn+1+3] += oy*dy;
         scu[nn+2+3] += oz*dy;
         dy = dx1*dzp;
         mm = nn + 3*mxv;
         samu[2*mm] += v1*dx;
         samu[2*mm+1] += v2*dx;
         samu[2*mm+2] += v3*dx;
         samu[2*mm+3] += v4*dx;
         samu[2*mm+4] += v5*dx;
         samu[2*mm+5] += v6*dx;
         sdcu[mm] += vx*dx;
         sdcu[mm+1] += vy*dx;
         sdcu[mm+2] += vz*dx;
         scu[mm] += ox*dx;
         scu[mm+1] += oy*dx;
         scu[mm+2] += oz*dx;
         samu[2*mm+6] += v1*dy;
         samu[2*mm+1+6] += v2*dy;
         samu[2*mm+2+6] += v3*dy;
         samu[2*mm+3+6] += v4*dy;
         samu[2*mm+4+6] += v5*dy;
         samu[2*mm+5+6] += v6*dy;
         sdcu[mm+3] += vx*dy;
         sdcu[mm+1+3] += vy*dy;
         sdcu[mm+2+3] += vz*dy;
         scu[mm+3] += ox*dy;
         scu[mm+1+3] += oy*dy;
         scu[mm+2+3] += oz*dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
            for (i = 1; i < nn; i++) {
               amu[6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[6*(i+mxv*j+mxyv*k)];
               amu[1+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[1+6*(i+mxv*j+mxyv*k)];
               amu[2+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[2+6*(i+mxv*j+mxyv*k)];
               amu[3+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[3+6*(i+mxv*j+mxyv*k)];
               amu[4+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[4+6*(i+mxv*j+mxyv*k)];
               amu[5+6*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[5+6*(i+mxv*j+mxyv*k)];
               dcu[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[3*(i+mxv*j+mxyv*k)];
               dcu[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[1+3*(i+mxv*j+mxyv*k)];
               dcu[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[2+3*(i+mxv*j+mxyv*k)];
               cu[3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[3*(i+mxv*j+mxyv*k)];
               cu[1+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+3*(i+mxv*j+mxyv*k)];
               cu[2+3*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+3*(i+mxv*j+mxyv*k)];
            }
         }
      }
/* deposit current to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            amu[6*(i+noff+nxv*(j+moff)+nxyv*loff)] += samu[6*(i+mxv*j)];
#pragma omp atomic
            amu[1+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[1+6*(i+mxv*j)];
#pragma omp atomic
            amu[2+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[2+6*(i+mxv*j)];
#pragma omp atomic
            amu[3+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[3+6*(i+mxv*j)];
#pragma omp atomic
            amu[4+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[4+6*(i+mxv*j)];
#pragma omp atomic
            amu[5+6*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += samu[5+6*(i+mxv*j)];
#pragma omp atomic
            dcu[3*(i+noff+nxv*(j+moff)+nxyv*loff)] += sdcu[3*(i+mxv*j)];
#pragma omp atomic
            dcu[1+3*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += sdcu[1+3*(i+mxv*j)];
#pragma omp atomic
            dcu[2+3*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += sdcu[2+3*(i+mxv*j)];
#pragma omp atomic
            cu[3*(i+noff+nxv*(j+moff)+nxyv*loff)] += scu[3*(i+mxv*j)];
#pragma omp atomic
            cu[1+3*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[1+3*(i+mxv*j)];
#pragma omp atomic
            cu[2+3*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[2+3*(i+mxv*j)];
            if (lm > mz) {
#pragma omp atomic
               amu[6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[1+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[1+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[2+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[2+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[3+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[3+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[4+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[4+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[5+6*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[5+6*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[3*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[1+3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[1+3*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[2+3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[2+3*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[3*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+3*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+3*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+3*(i+mxv*j+mxyv*(lm-1))];
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
            amu[6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[6*(i+mxyv*k)];
#pragma omp atomic
            amu[1+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[1+6*(i+mxyv*k)];
#pragma omp atomic
            amu[2+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[2+6*(i+mxyv*k)];
#pragma omp atomic
            amu[3+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[3+6*(i+mxyv*k)];
#pragma omp atomic
            amu[4+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[4+6*(i+mxyv*k)];
#pragma omp atomic
            amu[5+6*(i+noff+nxv*moff+nxyv*(k+loff))]
            += samu[5+6*(i+mxyv*k)];
#pragma omp atomic
            dcu[3*(i+noff+nxv*moff+nxyv*(k+loff))]
            += sdcu[3*(i+mxyv*k)];
#pragma omp atomic
            dcu[1+3*(i+noff+nxv*moff+nxyv*(k+loff))]
            += sdcu[1+3*(i+mxyv*k)];
#pragma omp atomic
            dcu[2+3*(i+noff+nxv*moff+nxyv*(k+loff))]
            += sdcu[2+3*(i+mxyv*k)];
#pragma omp atomic
            cu[3*(i+noff+nxv*moff+nxyv*(k+loff))] += scu[3*(i+mxyv*k)];
#pragma omp atomic
            cu[1+3*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[1+3*(i+mxyv*k)];
#pragma omp atomic
            cu[2+3*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[2+3*(i+mxyv*k)];
            if (mm > my) {
#pragma omp atomic
               amu[6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[1+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[1+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[2+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[2+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[3+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[3+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[4+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[4+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               amu[5+6*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += samu[5+6*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               dcu[3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += sdcu[3*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               dcu[1+3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += sdcu[1+3*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               dcu[2+3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += sdcu[2+3*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[3*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[1+3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[1+3*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[2+3*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[2+3*(i+mxv*(mm-1)+mxyv*k)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            amu[6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[1+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[1+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[2+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[2+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[3+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[3+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[4+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[4+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            amu[5+6*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += samu[5+6*(mxv*j+mxyv*k)];
#pragma omp atomic
            dcu[3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += sdcu[3*(mxv*j+mxyv*k)];
#pragma omp atomic
            dcu[1+3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += sdcu[1+3*(mxv*j+mxyv*k)];
#pragma omp atomic
            dcu[2+3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += sdcu[2+3*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[3*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[1+3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[1+3*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[2+3*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[2+3*(mxv*j+mxyv*k)];
            if (nm > mx) {
#pragma omp atomic
               amu[6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[1+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[1+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[2+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[2+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[3+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[3+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[4+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[4+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               amu[5+6*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += samu[5+6*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               dcu[3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[3*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               dcu[1+3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[1+3*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               dcu[2+3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += sdcu[2+3*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[3*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[1+3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+3*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[2+3*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+3*(nm-1+mxv*j+mxyv*k)];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            amu[6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[1+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[1+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[2+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[2+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[3+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[3+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[4+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[4+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            amu[5+6*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += samu[5+6*(i+mxyv*(lm-1))];
#pragma omp atomic
            dcu[3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += sdcu[3*(i+mxyv*(lm-1))];
#pragma omp atomic
            dcu[1+3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += sdcu[1+3*(i+mxyv*(lm-1))];
#pragma omp atomic
            dcu[2+3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += sdcu[2+3*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[3*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[1+3*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+3*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[2+3*(i+mxyv*(lm-1))];
            if (mm > my) {
#pragma omp atomic
               amu[6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[1+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[1+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[2+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[2+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[3+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[3+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[4+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[4+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               amu[5+6*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += samu[5+6*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               dcu[3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += sdcu[3*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               dcu[1+3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += sdcu[1+3*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               dcu[2+3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += sdcu[2+3*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[3*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[1+3*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+3*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[2+3*(i+mxv*(mm-1)+mxyv*(lm-1))];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            amu[6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[1+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[1+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[2+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[2+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[3+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[3+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[4+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[4+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            amu[5+6*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += samu[5+6*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            dcu[3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += sdcu[3*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            dcu[1+3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += sdcu[1+3*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            dcu[2+3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += sdcu[2+3*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[3*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[1+3*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+3*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[2+3*(mxv*j+mxyv*(lm-1))];
            if (nm > mx) {
#pragma omp atomic
               amu[6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[1+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[1+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[2+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[2+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[3+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[3+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[4+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[4+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               amu[5+6*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += samu[5+6*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[3*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[1+3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[1+3*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               dcu[2+3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += sdcu[2+3*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[3*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+3*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+3*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+3*(nm-1+mxv*j+mxyv*(lm-1))];
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
void cpporder3l(float ppart[], float ppbuff[], int kpic[], int ncl[],
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
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
   ppart[m][n][2] = position z of particle n in tile m
   ppbuff[l][n][i] = i co-ordinate of particle n in tile l
   kpic[l] = number of particles in tile l
   ncl[l][i] = number of particles going to destination i, tile l
   ihole[l][:][0] = location of hole in array left by departing particle
   ihole[l][:][1] = direction destination of particle leaving hole
   all for tile l
   ihole[l][0][0] = ih, number of holes left (error, if negative)
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
   int mxy1, mxyz1, noff, moff, loff, npp, ncoff;
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
private(j,k,l,noff,moff,loff,npp,nn,mm,ll,ih,nh,ist,dx,dy,dz,edgelx, \
edgely,edgelz,edgerx,edgery,edgerz)
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
         dx = ppart[idimp*(j+nppmx*l)];
         dy = ppart[1+idimp*(j+nppmx*l)];
         dz = ppart[2+idimp*(j+nppmx*l)];
/* find particles going out of bounds */
         ist = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* ist = direction particle is going                             */
         if (dx >= edgerx) {
            if (dx >= anx)
               ppart[idimp*(j+nppmx*l)] = dx - anx;
            ist = 2;
         }
         else if (dx < edgelx) {
            if (dx < 0.0) {
               dx += anx;
               if (dx < anx)
                  ist = 1;
               else
                  dx = 0.0;
               ppart[idimp*(j+nppmx*l)] = dx;
            }
            else {
               ist = 1;
            }
         }
         if (dy >= edgery) {
            if (dy >= any)
               ppart[1+idimp*(j+nppmx*l)] = dy - any;
            ist += 6;
         }
         else if (dy < edgely) {
            if (dy < 0.0) {
               dy += any;
               if (dy < any)
                  ist += 3;
               else
                  dy = 0.0;
               ppart[1+idimp*(j+nppmx*l)] = dy;
            }
            else {
               ist += 3;
            }
         }
         if (dz >= edgerz) {
            if (dz >= anz)
               ppart[2+idimp*(j+nppmx*l)] = dz - anz;
            ist += 18;
         }
         else if (dz < edgelz) {
            if (dz < 0.0) {
               dz += anz;
               if (dz < anz)
                  ist += 9;
               else
                  dz = 0.0;
               ppart[2+idimp*(j+nppmx*l)] = dz;
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
private(i,j,l,isum,ist,nh,ip,j1,ii)
   for (l = 0; l < mxyz1; l++) {
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
               ppbuff[i+idimp*(ii+npbmx*l)]
               = ppart[i+idimp*(j1+nppmx*l)];
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
private(i,j,k,l,ii,kk,npp,kx,ky,kz,kl,kr,kxl,kxr,lk,ll,lr,ih,nh,ncoff, \
ist,j1,j2,ip,ks)
   for (l = 0; l < mxyz1; l++) {
      npp = kpic[l];
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
                  ppart[i+idimp*(j1+nppmx*l)]
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
            j2 = ihole[2*(nh-j+(ntmax+1)*l)] - 1;
            if (j1 > j2) {
/* move particle only if it is below current hole */
               for (i = 0; i < idimp; i++) {
                  ppart[i+idimp*(j2+nppmx*l)]
                  = ppart[i+idimp*(j1+nppmx*l)];
               }
            }
         }
         npp -= ip;
      }
      kpic[l] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf3l(float ppart[], float ppbuff[], int kpic[], int ncl[],
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
   cgppushf3l procedure.
   input: all except ppbuff, irc
   output: ppart, ppbuff, kpic, ncl, irc
   ppart[m][n][0] = position x of particle n in tile m
   ppart[m][n][1] = position y of particle n in tile m
   ppart[m][n][2] = position z of particle n in tile m
   ppbuff[l][n][i] = i co-ordinate of particle n in tile l
   kpic[l] = number of particles in tile l
   ncl[l][i] = number of particles going to destination i, tile l
   ihole[l][:][0] = location of hole in array left by departing particle
   ihole[l][:][1] = direction destination of particle leaving hole
   all for tile l
   ihole[l][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mz1 = (system length in z direction - 1)/mz + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int mxy1, mxyz1, npp, ncoff;
   int i, j, k, l, ii, kx, ky, kz, ih, nh, ist, ll, isum;
   int ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr;
   int ks[26];
   mxy1 = mx1*my1;
   mxyz1 = mxy1*mz1;
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,l,isum,ist,nh,ip,j1,ii)
   for (l = 0; l < mxyz1; l++) {
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
               ppbuff[i+idimp*(ii+npbmx*l)]
               = ppart[i+idimp*(j1+nppmx*l)];
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
private(i,j,k,l,ii,kk,npp,kx,ky,kz,kl,kr,kxl,kxr,lk,ll,lr,ih,nh,ncoff, \
ist,j1,j2,ip,ks)
   for (l = 0; l < mxyz1; l++) {
      npp = kpic[l];
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
                  ppart[i+idimp*(j1+nppmx*l)]
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
            j2 = ihole[2*(nh-j+(ntmax+1)*l)] - 1;
            if (j1 > j2) {
/* move particle only if it is below current hole */
               for (i = 0; i < idimp; i++) {
                  ppart[i+idimp*(j2+nppmx*l)]
                  = ppart[i+idimp*(j1+nppmx*l)];
               }
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
   int j, k, l, nxye3, ll;
   nxye3 = 3*nxe*nye;
/* copy edges of extended field */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,ll)
      for (l = 0; l < nz; l++) {
         ll = nxye3*l;
         for (k = 0; k < ny; k++) {
            fxyz[3*(nx+nxe*k)+ll] = fxyz[3*nxe*k+ll];
            fxyz[1+3*(nx+nxe*k)+ll] = fxyz[1+3*nxe*k+ll];
            fxyz[2+3*(nx+nxe*k)+ll] = fxyz[2+3*nxe*k+ll];
         }
         for (j = 0; j < nx; j++) {
            fxyz[3*(j+nxe*ny)+ll] = fxyz[3*j+ll];
            fxyz[1+3*(j+nxe*ny)+ll] = fxyz[1+3*j+ll];
            fxyz[2+3*(j+nxe*ny)+ll] = fxyz[2+3*j+ll];
         }
         fxyz[3*(nx+nxe*ny)+ll] = fxyz[ll];
         fxyz[1+3*(nx+nxe*ny)+ll] = fxyz[1+ll];
         fxyz[2+3*(nx+nxe*ny)+ll] = fxyz[2+ll];
      }
#pragma omp for private(j,k)
      for (k = 0; k < ny; k++) {
         for (j = 0; j < nx; j++) {
            fxyz[3*(j+nxe*k)+nxye3*nz] = fxyz[3*(j+nxe*k)];
            fxyz[1+3*(j+nxe*k)+nxye3*nz] = fxyz[1+3*(j+nxe*k)];
            fxyz[2+3*(j+nxe*k)+nxye3*nz] = fxyz[2+3*(j+nxe*k)];
         }
         fxyz[3*(nx+nxe*k)+nxye3*nz] = fxyz[3*nxe*k];
         fxyz[1+3*(nx+nxe*k)+nxye3*nz] = fxyz[1+3*nxe*k];
         fxyz[2+3*(nx+nxe*k)+nxye3*nz] = fxyz[2+3*nxe*k];
      }
   }
   for (j = 0; j < nx; j++) {
      fxyz[3*(j+nxe*ny)+nxye3*nz] = fxyz[3*j];
      fxyz[1+3*(j+nxe*ny)+nxye3*nz] = fxyz[1+3*j];
      fxyz[2+3*(j+nxe*ny)+nxye3*nz] = fxyz[2+3*j];
   }
   fxyz[3*(nx+nxe*ny)+nxye3*nz] = fxyz[0];
   fxyz[1+3*(nx+nxe*ny)+nxye3*nz] = fxyz[1];
   fxyz[2+3*(nx+nxe*ny)+nxye3*nz] = fxyz[2];
   return;
}

/*--------------------------------------------------------------------*/
void cacguard3l(float cu[], int nx, int ny, int nz, int nxe, int nye,
                int nze) {
/* accumulate extended periodic vector field cu
   linear interpolation
   nx/ny/nz = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   nze = third dimension of field arrays, must be >= nz+1
local data                                                            */
   int j, k, l, nxye3, ll;
   nxye3 = 3*nxe*nye;
/* accumulate edges of extended field */
#pragma omp parallel
   {
#pragma omp for \
private(j,k,l,ll)
      for (l = 0; l < nz; l++) {
         ll = nxye3*l;
         for (k = 0; k < ny; k++) {
            cu[3*nxe*k+ll] += cu[3*(nx+nxe*k)+ll];
            cu[1+3*nxe*k+ll] += cu[1+3*(nx+nxe*k)+ll];
            cu[2+3*nxe*k+ll] += cu[2+3*(nx+nxe*k)+ll];
            cu[3*(nx+nxe*k)+ll] = 0.0;
            cu[1+3*(nx+nxe*k)+ll] = 0.0;
            cu[2+3*(nx+nxe*k)+ll] = 0.0;
         }
         for (j = 0; j < nx; j++) {
            cu[3*j+ll] += cu[3*(j+nxe*ny)+ll];
            cu[1+3*j+ll] += cu[1+3*(j+nxe*ny)+ll];
            cu[2+3*j+ll] += cu[2+3*(j+nxe*ny)+ll];
            cu[3*(j+nxe*ny)+ll] = 0.0;
            cu[1+3*(j+nxe*ny)+ll] = 0.0;
            cu[2+3*(j+nxe*ny)+ll] = 0.0;
         }
         cu[ll] += cu[3*(nx+nxe*ny)+ll];
         cu[1+ll] += cu[1+3*(nx+nxe*ny)+ll];
         cu[2+ll] += cu[2+3*(nx+nxe*ny)+ll];
         cu[3*(nx+nxe*ny)+ll] = 0.0;
         cu[1+3*(nx+nxe*ny)+ll] = 0.0;
         cu[2+3*(nx+nxe*ny)+ll] = 0.0;
      }
#pragma omp for private(j,k)
      for (k = 0; k < ny; k++) {
         for (j = 0; j < nx; j++) {
            cu[3*(j+nxe*k)] += cu[3*(j+nxe*k)+nxye3*nz];
            cu[1+3*(j+nxe*k)] += cu[1+3*(j+nxe*k)+nxye3*nz];
            cu[2+3*(j+nxe*k)] += cu[2+3*(j+nxe*k)+nxye3*nz];
            cu[3*(j+nxe*k)+nxye3*nz] = 0.0;
            cu[1+3*(j+nxe*k)+nxye3*nz] = 0.0;
            cu[2+3*(j+nxe*k)+nxye3*nz] = 0.0;
         }
         cu[3*nxe*k] += cu[3*(nx+nxe*k)+nxye3*nz];
         cu[1+3*nxe*k] += cu[1+3*(nx+nxe*k)+nxye3*nz];
         cu[2+3*nxe*k] += cu[2+3*(nx+nxe*k)+nxye3*nz];
         cu[3*(nx+nxe*k)+nxye3*nz] = 0.0;
         cu[1+3*(nx+nxe*k)+nxye3*nz] = 0.0;
         cu[2+3*(nx+nxe*k)+nxye3*nz] = 0.0;
      }
   }
   for (j = 0; j < nx; j++) {
      cu[3*j] += cu[3*(j+nxe*ny)+nxye3*nz];
      cu[1+3*j] += cu[1+3*(j+nxe*ny)+nxye3*nz];
      cu[2+3*j] += cu[2+3*(j+nxe*ny)+nxye3*nz];
      cu[3*(j+nxe*ny)+nxye3*nz] = 0.0;
      cu[1+3*(j+nxe*ny)+nxye3*nz] = 0.0;
      cu[2+3*(j+nxe*ny)+nxye3*nz] = 0.0;
   }
   cu[0] += cu[3*(nx+nxe*ny)+nxye3*nz];
   cu[1] += cu[1+3*(nx+nxe*ny)+nxye3*nz];
   cu[2] += cu[2+3*(nx+nxe*ny)+nxye3*nz];
   cu[3*(nx+nxe*ny)+nxye3*nz] = 0.0;
   cu[1+3*(nx+nxe*ny)+nxye3*nz] = 0.0;
   cu[2+3*(nx+nxe*ny)+nxye3*nz] = 0.0;
   return;
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
#pragma omp for private(j,k)
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
void camcguard3l(float amu[], int nx, int ny, int nz, int nxe, int nye,
                 int nze, int ndim) {
/* accumulate extended periodic tensor field amu
   linear interpolation
   nx/ny/nz = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   nze = third dimension of field arrays, must be >= nz+1
   ndim = number of elements in tensor
local data                                                            */
   int i, j, k, l, nnxye, ll;
   nnxye = ndim*nxe*nye;
/* accumulate edges of extended field */
#pragma omp parallel
   {
#pragma omp for \
private(i,j,k,l,ll)
      for (l = 0; l < nz; l++) {
         ll = nnxye*l;
         for (k = 0; k < ny; k++) {
            for (i = 0; i < ndim; i++) {
               amu[i+ndim*nxe*k+ll] += amu[i+ndim*(nx+nxe*k)+ll];
               amu[i+ndim*(nx+nxe*k)+ll] = 0.0;
            }
         }
         for (j = 0; j < nx; j++) {
            for (i = 0; i < ndim; i++) {
               amu[i+ndim*j+ll] += amu[i+ndim*(j+nxe*ny)+ll];
               amu[i+ndim*(j+nxe*ny)+ll] = 0.0;
            }
         }
         for (i = 0; i < ndim; i++) {
            amu[i+ll] += amu[i+ndim*(nx+nxe*ny)+ll];
            amu[i+ndim*(nx+nxe*ny)+ll] = 0.0;
         }
      }
#pragma omp for private(i,j,k)
      for (k = 0; k < ny; k++) {
         for (j = 0; j < nx; j++) {
            for (i = 0; i < ndim; i++) {
               amu[i+ndim*(j+nxe*k)] += amu[i+ndim*(j+nxe*k)+nnxye*nz];
               amu[i+ndim*(j+nxe*k)+nnxye*nz] = 0.0;
            }
         }
         for (i = 0; i < ndim; i++) {
            amu[i+ndim*nxe*k] += amu[i+ndim*(nx+nxe*k)+nnxye*nz];
            amu[i+ndim*(nx+nxe*k)+nnxye*nz] = 0.0;
         }
      }
   }
   for (j = 0; j < nx; j++) {
      for (i = 0; i < ndim; i++) {
         amu[i+ndim*j] += amu[i+ndim*(j+nxe*ny)+nnxye*nz];
         amu[i+ndim*(j+nxe*ny)+nnxye*nz] = 0.0;
      }
   }
   for (i = 0; i < ndim; i++) {
      amu[i] += amu[i+ndim*(nx+nxe*ny)+nnxye*nz];
      amu[i+ndim*(nx+nxe*ny)+nnxye*nz] = 0.0;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cascfguard3l(float dcu[], float cus[], float q2m0, int nx, int ny,
                  int nz, int nxe, int nye, int nze) {
/* add scaled field to extended periodic field */
/* local data */
   int i, j, k, l, nxye3, ll;
   nxye3 = 3*nxe*nye;
#pragma omp for private(i,j,k,l,ll)
   for (l = 0; l < nz; l++) {
      ll = nxye3*l;
      for (k = 0; k < ny; k++) {
         for (j = 0; j < nx; j++) {
           for (i = 0; i < 3; i++) {
               dcu[i+3*(j+nxe*k)+ll] -= q2m0*cus[i+3*(j+nxe*k)+ll];
            }
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfwpminmx3(float qe[], float qbme, float *wpmax, float *wpmin,
                int nx, int ny, int nz, int nxe, int nye, int nze) {
/* calculates maximum and minimum plasma frequency. assumes guard cells
   have already been added
   qe = charge density for electrons
   qbme = charge/mass ratio for electrons
   wpmax/wpmin = maximum/minimum plasma frequency
   nx/ny/nz = system length in x/y/z direction
   nxe = first dimension of charge arrays, nxe must be >= nx
   nye = second dimension of charge arrays, nye must be >= ny
   nze = third dimension of charge arrays, nze must be >= nz
local data                                                 */
   int j, k, l, nxye, ll;
   float tpmax, tpmin, at1;
   nxye = nxe*nye;
   tpmax = qbme*qe[0];
   tpmin = tpmax;
#pragma omp for private(j,k,l,ll)
   for (l = 0; l < nz; l++) {
      ll = nxye*l;
      for (k = 0; k < ny; k++) {
         for (j = 0; j < nx; j++) {
            at1 = qbme*qe[j+nxe*k+ll];
#pragma omp critical
            tpmax = at1 > tpmax ? at1 : tpmax;
#pragma omp critical
            tpmin = at1 < tpmin ? at1 : tpmin;
         }
      }
   }
   *wpmax = tpmax;
   *wpmin = tpmin;
   return;
}

/*--------------------------------------------------------------------*/
void cmpois33(float complex q[], float complex fxyz[], int isign,
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
local data                                                            */
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
            for (j = 1; j < nxh; j++) {
               at1 = crealf(ffc[j+kk+ll])*cimagf(ffc[j+kk+ll]);
               at2 = at1*dnx*(float) j;
               at3 = dky*at1;
               at4 = dkz*at1;
               zt1 = cimagf(q[j+kj+lj]) - crealf(q[j+kj+lj])*_Complex_I;
               zt2 = cimagf(q[j+k1+lj]) - crealf(q[j+k1+lj])*_Complex_I;
               fxyz[3*(j+kj+lj)] = at2*zt1;
               fxyz[1+3*(j+kj+lj)] = at3*zt1;
               fxyz[2+3*(j+kj+lj)] = at4*zt1;
               fxyz[3*(j+k1+lj)] = at2*zt2;
               fxyz[1+3*(j+k1+lj)] = -at3*zt2;
               fxyz[2+3*(j+k1+lj)] = at4*zt2;
               zt1 = cimagf(q[j+kj+l1]) - crealf(q[j+kj+l1])*_Complex_I;
               zt2 = cimagf(q[j+k1+l1]) - crealf(q[j+k1+l1])*_Complex_I;
               fxyz[3*(j+kj+l1)] = at2*zt1;
               fxyz[1+3*(j+kj+l1)] = at3*zt1;
               fxyz[2+3*(j+kj+l1)] = -at4*zt1;
               fxyz[3*(j+k1+l1)] = at2*zt2;
               fxyz[1+3*(j+k1+l1)] = -at3*zt2;
               fxyz[2+3*(j+k1+l1)] = -at4*zt2;
               wp += at1*(q[j+kj+lj]*conjf(q[j+kj+lj])                 
                  + q[j+k1+lj]*conjf(q[j+k1+lj])
                  + q[j+kj+l1]*conjf(q[j+kj+l1])      
                  + q[j+k1+l1]*conjf(q[j+k1+l1]));
            }
         }
/* mode numbers kx = 0, nx/2 */
         for (k = 1; k < nyh; k++) {
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            at1 = crealf(ffc[kk+ll])*cimagf(ffc[kk+ll]);
            at3 = at1*dny*(float) k;
            at4 = dkz*at1;
            zt1 = cimagf(q[kj+lj]) - crealf(q[kj+lj])*_Complex_I;
            zt2 = cimagf(q[kj+l1]) - crealf(q[kj+l1])*_Complex_I;
            fxyz[3*(kj+lj)] = zero;
            fxyz[1+3*(kj+lj)] = at3*zt1;
            fxyz[2+3*(kj+lj)] = at4*zt1;
            fxyz[3*(k1+lj)] = zero;
            fxyz[1+3*(k1+lj)] = zero;
            fxyz[2+3*(k1+lj)] = zero;
            fxyz[3*(kj+l1)] = zero;
            fxyz[1+3*(kj+l1)] = at3*zt2;
            fxyz[2+3*(kj+l1)] = -at4*zt2;
            fxyz[3*(k1+l1)] = zero;
            fxyz[1+3*(k1+l1)] = zero;
            fxyz[2+3*(k1+l1)] = zero;
            wp += at1*(q[kj+lj]*conjf(q[kj+lj])
               + q[kj+l1]*conjf(q[kj+l1]));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
         for (j = 1; j < nxh; j++) {
            at1 = crealf(ffc[j+ll])*cimagf(ffc[j+ll]);
            at2 = at1*dnx*(float) j;  
            at4 = dkz*at1;
            zt1 = cimagf(q[j+lj]) - crealf(q[j+lj])*_Complex_I;
            zt2 = cimagf(q[j+l1]) - crealf(q[j+l1])*_Complex_I;
            fxyz[3*(j+lj)] = at2*zt1;
            fxyz[1+3*(j+lj)] = zero;
            fxyz[2+3*(j+lj)] = at4*zt1;
            fxyz[3*(j+k1+lj)] = zero;
            fxyz[1+3*(j+k1+lj)] = zero;
            fxyz[2+3*(j+k1+lj)] = zero;
            fxyz[3*(j+l1)] = at2*zt2;
            fxyz[1+3*(j+l1)] = zero;
            fxyz[2+3*(j+l1)] = -at4*zt2;
            fxyz[3*(j+k1+l1)] = zero;
            fxyz[1+3*(j+k1+l1)] = zero;
            fxyz[2+3*(j+k1+l1)] = zero;
            wp += at1*(q[j+lj]*conjf(q[j+lj])                           
               + q[j+l1]*conjf(q[j+l1]));
         }
/* mode numbers kx = 0, nx/2 */
         at1 = crealf(ffc[ll])*cimagf(ffc[ll]);
         at4 = dkz*at1;
         zt1 = cimagf(q[lj]) - crealf(q[lj])*_Complex_I;
         fxyz[3*lj] = zero;
         fxyz[1+3*lj] = zero;
         fxyz[2+3*lj] = at4*zt1;
         fxyz[3*(k1+lj)] = zero;
         fxyz[1+3*(k1+lj)] = zero;
         fxyz[2+3*(k1+lj)] = zero;
         fxyz[3*l1] = zero;
         fxyz[1+3*l1] = zero;
         fxyz[2+3*l1] = zero;
         fxyz[3*(k1+l1)] = zero;
         fxyz[1+3*(k1+l1)] = zero;
         fxyz[2+3*(k1+l1)] = zero;
         wp += at1*(q[lj]*conjf(q[lj]));
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
      for (j = 1; j < nxh; j++) {
         at1 = crealf(ffc[j+kk])*cimagf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         zt1 = cimagf(q[j+kj]) - crealf(q[j+kj])*_Complex_I;
         zt2 = cimagf(q[j+k1]) - crealf(q[j+k1])*_Complex_I;
         fxyz[3*(j+kj)] = at2*zt1;
         fxyz[1+3*(j+kj)] = at3*zt1;
         fxyz[2+3*(j+kj)] = zero;
         fxyz[3*(j+k1)] = at2*zt2;
         fxyz[1+3*(j+k1)] = -at3*zt2;
         fxyz[2+3*(j+k1)] = zero;
         fxyz[3*(j+kj+l1)] = zero;
         fxyz[1+3*(j+kj+l1)] = zero;
         fxyz[2+3*(j+kj+l1)] = zero;
         fxyz[3*(j+k1+l1)] = zero;
         fxyz[1+3*(j+k1+l1)] = zero;
         fxyz[2+3*(j+k1+l1)] = zero;
         wp += at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1]));
      }
/* mode numbers kx = 0, nx/2 */
      at1 = crealf(ffc[kk])*cimagf(ffc[kk]);
      at3 = at1*dny*(float) k;
      zt1 = cimagf(q[kj]) - crealf(q[kj])*_Complex_I;
      fxyz[3*kj] = zero;
      fxyz[1+3*kj] = at3*zt1;
      fxyz[2+3*kj] = zero;
      fxyz[3*k1] = zero;
      fxyz[1+3*k1] = zero;
      fxyz[2+3*k1] = zero;
      fxyz[3*(kj+l1)] = zero;
      fxyz[1+3*(kj+l1)] = zero;
      fxyz[2+3*(kj+l1)] = zero;
      fxyz[3*(k1+l1)] = zero;
      fxyz[1+3*(k1+l1)] = zero;
      fxyz[2+3*(k1+l1)] = zero;
      wp += at1*(q[kj]*conjf(q[kj]));
      sum2 += wp;
   }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      at1 = crealf(ffc[j])*cimagf(ffc[j]);
      at2 = at1*dnx*(float) j;
      zt1 = cimagf(q[j]) - crealf(q[j])*_Complex_I;
      fxyz[3*j] = at2*zt1;
      fxyz[1+3*j] = zero;
      fxyz[2+3*j] = zero;
      fxyz[3*(j+k1)] = zero;
      fxyz[1+3*(j+k1)] = zero;
      fxyz[2+3*(j+k1)] = zero;
      fxyz[3*(j+l1)] = zero;
      fxyz[1+3*(j+l1)] = zero;
      fxyz[2+3*(j+l1)] = zero;
      fxyz[3*(j+k1+l1)] = zero;
      fxyz[1+3*(j+k1+l1)] = zero;
      fxyz[2+3*(j+k1+l1)] = zero;
      wp += at1*(q[j]*conjf(q[j]));
   }
   fxyz[0] = zero;
   fxyz[1] = zero;
   fxyz[2] = zero;
   fxyz[3*k1] = zero;
   fxyz[1+3*k1] = zero;
   fxyz[2+3*k1] = zero;
   fxyz[3*l1] = zero;
   fxyz[1+3*l1] = zero;
   fxyz[2+3*l1] = zero;
   fxyz[3*(k1+l1)] = zero;
   fxyz[1+3*(k1+l1)] = zero;
   fxyz[2+3*(k1+l1)] = zero;
   *we = (sum1 + sum2 + wp)*((float) nx)*((float) ny)*((float) nz);
   return;
}

/*--------------------------------------------------------------------*/
void cmcuperp3(float complex cu[], int nx, int ny, int nz, int nxvh,
               int nyv, int nzv) {
/* this subroutine calculates the transverse current in fourier space
   input: all, output: cu
   approximate flop count is:
   100*nxc*nyc*nzc + 36*(nxc*nyc + nxc*nzc + nyc*nzc)
   and (nx/2)*nyc*nzc divides
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   the transverse current is calculated using the equation:
   cux[kz][ky][kx] = cux[kz][ky][kx]
                   - kx*(kx*cux[kz][ky][kx]+ky*cuy[kz][ky][kx]
                   +     kz*cuz[kz][ky][kx])/(kx*kx+ky*ky+kz*kz)
   cuy[kz][ky][kx] = cuy[kz][ky][kx]
                   - ky*(kx*cux([kz][ky][kx]+ky*cuy[kz][ky][kx]
                   +     kz*cuz[kz][ky][kx])/(kx*kx+ky*ky+kz*kz)
   cuz[kz][ky][kx] = cuz[kz][ky][kx]
                   - kz*(kx*cux[kz][ky][kx]+ky*cuy[kz][ky][kx]
                   +     kz*cuz[kz][ky][kx])/(kx*kx+ky*ky+kz*kz)
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
   j,k,l = fourier mode numbers, except for
   cux(kx=pi) = cuy(kx=pi) = cuz(kx=pi) = 0,
   cux(ky=pi) = cuy(ky=pi) = cux(ky=pi) = 0,
   cux(kz=pi) = cuy(kz=pi) = cuz(kz=pi) = 0,
   cux(kx=0,ky=0,kz=0) = cuy(kx=0,ky=0,kz=0) = cuz(kx=0,ky=0,kz=0) = 0.
   cu[l][k][j][i] = complex current density for fourier mode (j,k,l)
   nx/ny/nz = system length in x/y/z direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nzv = fourth dimension of field arrays, must be >= nz
local data                                                 */
   int nxh, nyh, nzh, j, k, l, k1, l1, kj, lj, nxvyh;
   float dnx, dny, dnz, dkx, dky, dkz, dky2, dkz2, dkyz2, at1;
   float complex zero, zt1;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   zero = 0.0 + 0.0*_Complex_I;
/* calculate transverse part of current */
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,k1,l1,lj,kj,dkx,dky,dkz,dkz2,dkyz2,at1,zt1)
      for (l = 1; l < nzh; l++) {
         dkz = dnz*(float) l;
         lj = nxvyh*l;
         l1 = nxvyh*nz - lj;
         dkz2 = dkz*dkz;
         for (k = 1; k < nyh; k++) {
            dky = dny*(float) k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            dkyz2 = dky*dky + dkz2;
            for (j = 1; j < nxh; j++) {
               dkx = dnx*(float) j;
               at1 = 1.0/(dkx*dkx + dkyz2);
               zt1 = at1*(dkx*cu[3*(j+kj+lj)] + dky*cu[1+3*(j+kj+lj)]
                        + dkz*cu[2+3*(j+kj+lj)]);
               cu[3*(j+kj+lj)] -= dkx*zt1;
               cu[1+3*(j+kj+lj)] -= dky*zt1;
               cu[2+3*(j+kj+lj)] -= dkz*zt1;
               zt1 = at1*(dkx*cu[3*(j+k1+lj)] - dky*cu[1+3*(j+k1+lj)]
                        + dkz*cu[2+3*(j+k1+lj)]);
               cu[3*(j+k1+lj)] -= dkx*zt1;
               cu[1+3*(j+k1+lj)] += dky*zt1;
               cu[2+3*(j+k1+lj)] -= dkz*zt1;
               zt1 = at1*(dkx*cu[3*(j+kj+l1)] + dky*cu[1+3*(j+kj+l1)]
                        - dkz*cu[2+3*(j+kj+l1)]);
               cu[3*(j+kj+l1)] -= dkx*zt1;
               cu[1+3*(j+kj+l1)] -= dky*zt1;
               cu[2+3*(j+kj+l1)] += dkz*zt1;
               zt1 = at1*(dkx*cu[3*(j+k1+l1)] - dky*cu[1+3*(j+k1+l1)]
                        - dkz*cu[2+3*(j+k1+l1)]);
               cu[3*(j+k1+l1)] -= dkx*zt1;
               cu[1+3*(j+k1+l1)] += dky*zt1;
               cu[2+3*(j+k1+l1)] += dkz*zt1;
            }
         }
/* mode numbers kx = 0, nx/2 */
         for (k = 1; k < nyh; k++) {
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            dky = dny*(float) k;
            at1 = 1.0/(dky*dky + dkz2);
            zt1 = at1*(dky*cu[1+3*(kj+lj)] + dkz*cu[2+3*(kj+lj)]);
            cu[1+3*(kj+lj)] -= dky*zt1;
            cu[2+3*(kj+lj)] -= dkz*zt1;
            cu[3*(k1+lj)] = zero;
            cu[1+3*(k1+lj)] = zero;
            cu[2+3*(k1+lj)] = zero;
            zt1 = at1*(dky*cu[1+3*(kj+l1)] - dkz*cu[2+3*(kj+l1)]);
            cu[1+3*(kj+l1)] -= dky*zt1;
            cu[2+3*(kj+l1)] += dkz*zt1;
            cu[3*(k1+l1)] = zero;
            cu[1+3*(k1+l1)] = zero;
            cu[2+3*(k1+l1)] = zero;
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
         for (j = 1; j < nxh; j++) {
            dkx = dnx*(float) j;
            at1 = 1.0/(dkx*dkx + dkz2);
            zt1 = at1*(dkx*cu[3*(j+lj)] + dkz*cu[2+3*(j+lj)]);
            cu[3*(j+lj)] -= dkx*zt1;
            cu[2+3*(j+lj)] -= dkz*zt1;
            cu[3*(j+k1+lj)] = zero;
            cu[1+3*(j+k1+lj)] = zero;
            cu[2+3*(j+k1+lj)] = zero;
            zt1 = at1*(dkx*cu[3*(j+l1)] - dkz*cu[2+3*(j+l1)]);
            cu[3*(j+l1)] -= dkx*zt1;
            cu[2+3*(j+l1)] += dkz*zt1;
            cu[3*(j+k1+l1)] = zero;
            cu[1+3*(j+k1+l1)] = zero;
            cu[2+3*(j+k1+l1)] = zero;
         }
/* mode numbers kx = 0, nx/2 */
         cu[2+3*lj] = zero;
         cu[3*(k1+lj)] = zero;
         cu[1+3*(k1+lj)] = zero;
         cu[2+3*(k1+lj)] = zero;
         cu[3*l1] = zero;
         cu[1+3*l1] = zero;
         cu[2+3*l1] = zero;
         cu[3*(k1+l1)] = zero;
         cu[1+3*(k1+l1)] = zero;
         cu[2+3*(k1+l1)] = zero;
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
#pragma omp parallel for \
private(j,k,k1,kj,dky,dky2,dkx,at1,zt1)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      dky2 = dky*dky;
      for (j = 1; j < nxh; j++) {
         dkx = dnx*(float) j;
         at1 = 1.0/(dkx*dkx + dky2);
         zt1 = at1*(dkx*cu[3*(j+kj)] + dky*cu[1+3*(j+kj)]);
         cu[3*(j+kj)] -= dkx*zt1;
         cu[1+3*(j+kj)] -= dky*zt1;
         zt1 = at1*(dkx*cu[3*(j+k1)]- dky*cu[1+3*(j+k1)]);
         cu[3*(j+k1)] -= dkx*zt1;
         cu[1+3*(j+k1)] += dky*zt1;
         cu[3*(j+kj+l1)] = zero;
         cu[1+3*(j+kj+l1)] = zero;
         cu[2+3*(j+kj+l1)] = zero;
         cu[3*(j+k1+l1)] = zero;
         cu[1+3*(j+k1+l1)] = zero;
         cu[2+3*(j+k1+l1)] = zero;
      }
/* mode numbers kx = 0, nx/2 */
      cu[1+3*kj] = zero;
      cu[3*k1] = zero;
      cu[1+3*k1] = zero;
      cu[2+3*k1] = zero;
      cu[3*(kj+l1)] = zero;
      cu[1+3*(kj+l1)] = zero;
      cu[2+3*(kj+l1)] = zero;
      cu[3*(k1+l1)] = zero;
      cu[1+3*(k1+l1)] = zero;
      cu[2+3*(k1+l1)] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      cu[3*j] = zero;
      cu[3*(j+k1)] = zero;
      cu[1+3*(j+k1)] = zero;
      cu[2+3*(j+k1)] = zero;
      cu[3*(j+l1)] = zero;
      cu[1+3*(j+l1)] = zero;
      cu[2+3*(j+l1)] = zero;
      cu[3*(j+k1+l1)] = zero;
      cu[1+3*(j+k1+l1)] = zero;
      cu[2+3*(j+k1+l1)] = zero;
   }
   cu[0] = zero;
   cu[1] = zero;
   cu[2] = zero;
   cu[3*k1] = zero;
   cu[1+3*k1] = zero;
   cu[2+3*k1] = zero;
   cu[3*l1] = zero;
   cu[1+3*l1] = zero;
   cu[2+3*l1] = zero;
   cu[3*(k1+l1)] = zero;
   cu[1+3*(k1+l1)] = zero;
   cu[2+3*(k1+l1)] = zero;
   return;
}

/*--------------------------------------------------------------------*/
void cmbbpois33(float complex cu[], float complex bxyz[],
                float complex ffc[], float ci, float *wm, int nx,
                int ny, int nz, int nxvh, int nyv, int nzv, int nxhd,
                int nyhd, int nzhd) {
/* this subroutine solves 3d poisson's equation in fourier space for
   magnetic field (or convolution of magnetic field over particle shape)
   with periodic boundary conditions.
   input: cu,ffc,ci,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
   output: bxyz, wm
   approximate flop count is:
   193*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   the magnetic field is calculated using the equations:
   bx[kz][ky][kx] = ci*ci*sqrt(-1)*g[kz][ky][kx]*
                  (ky*cuz[kz][ky][kx]-kz*cuy[kz][ky][kx])*s[kz][ky][kx],
   by[kz][ky][kx] = ci*ci*sqrt(-1)*g[kz][ky][kx]*
                  (kz*cux[kz][ky][kx]-kx*cuz[kz][ky][kx])*s[kz][ky][kx],
   bz[kz][ky][kx] = ci*ci*sqrt(-1)*g[kz][ky][kx]*
                  (kx*cuy[kz][ky][kx]-ky*cux[kz][ky][kx])*s[kz][ky][kx],
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
   j,k,l = fourier mode numbers,
   g[kz][ky][kx] = (affp/(kx**2+ky**2+kz**2))*s[kz][ky][kx],
   s[kz][ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
   bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
   bx(ky=pi) = by(ky=pi) = bx(ky=pi) = 0,
   bx(kz=pi) = by(kz=pi) = bz(kz=pi) = 0,
   bx(kx=0,ky=0,kz=0) = by(kx=0,ky=0,kz=0) = bz(kx=0,ky=0,kz=0) = 0.
   cu[l][k][j][i] = complex current density for fourier mode (j,k,l)
   bxy[l][k][j][0] = x component of complex magnetic field
   bxy[l][k][j][1] = y component of complex magnetic field
   bxy[l][k][j][2] = z component of complex magnetic field
   all for fourier mode (j,k,l)
   aimag(ffc[l][k][j]) = finite-size particle shape factor s
   for fourier mode (j,k,l)
   real(ffc[l][k][j]) = potential green's function g
   for fourier mode (j,k,l)
   ci = reciprocal of velocity of light
   magnetic field energy is also calculated, using
   wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
      |cu[kz][ky][kx]*s[kz][ky][kx]|**2), where
   affp = normalization constant = nx*ny/np, where np=number of particles
   this expression is valid only if the current is divergence-free
   nx/ny/nz = system length in x/y/z direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nzv = fourth dimension of field arrays, must be >= nz
   nxhd = dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
   nzhd = third dimension of form factor array, must be >= nzh
local data                                                 */
   int nxh, nyh, nzh, j, k, l, k1, l1, kk, kj, ll, lj, nxyhd, nxvyh;
   float dnx, dny, dnz, dky, dkz, ci2, at1, at2, at3, at4;
   float complex zero, zt1, zt2, zt3;
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
   ci2 = ci*ci;
/* calculate smoothed magnetic field and sum field energy */
   sum1 = 0.0;
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,k1,l1,ll,lj,kk,kj,dky,dkz,at1,at2,at3,at4,zt1,zt2,zt3,wp) \
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
            for (j = 1; j < nxh; j++) {
               at1 = ci2*crealf(ffc[j+kk+ll])*cimagf(ffc[j+kk+ll]);
               at2 = at1*dnx*(float) j;
               at3 = dky*at1;
               at4 = dkz*at1;
               zt1 = -cimagf(cu[2+3*(j+kj+lj)])
                    + crealf(cu[2+3*(j+kj+lj)])*_Complex_I;
               zt2 = -cimagf(cu[1+3*(j+kj+lj)])
                    + crealf(cu[1+3*(j+kj+lj)])*_Complex_I;
               zt3 = -cimagf(cu[3*(j+kj+lj)])
                    + crealf(cu[3*(j+kj+lj)])*_Complex_I;
               bxyz[3*(j+kj+lj)] = at3*zt1 - at4*zt2;
               bxyz[1+3*(j+kj+lj)] = at4*zt3 - at2*zt1;
               bxyz[2+3*(j+kj+lj)] = at2*zt2 - at3*zt3;
               zt1 = -cimagf(cu[2+3*(j+k1+lj)])
                    + crealf(cu[2+3*(j+k1+lj)])*_Complex_I;
               zt2 = -cimagf(cu[1+3*(j+k1+lj)])
                    + crealf(cu[1+3*(j+k1+lj)])*_Complex_I;
               zt3 = -cimagf(cu[3*(j+k1+lj)])
                    + crealf(cu[3*(j+k1+lj)])*_Complex_I;
               bxyz[3*(j+k1+lj)] = -at3*zt1 - at4*zt2;
               bxyz[1+3*(j+k1+lj)] = at4*zt3 - at2*zt1;
               bxyz[2+3*(j+k1+lj)] = at2*zt2 + at3*zt3;
               zt1 = -cimagf(cu[2+3*(j+kj+l1)])
                    + crealf(cu[2+3*(j+kj+l1)])*_Complex_I;
               zt2 = -cimagf(cu[1+3*(j+kj+l1)])
                    + crealf(cu[1+3*(j+kj+l1)])*_Complex_I;
               zt3 = -cimagf(cu[3*(j+kj+l1)])
                    + crealf(cu[3*(j+kj+l1)])*_Complex_I;
               bxyz[3*(j+kj+l1)] = at3*zt1 + at4*zt2;
               bxyz[1+3*(j+kj+l1)] = -at4*zt3 - at2*zt1;
               bxyz[2+3*(j+kj+l1)] = at2*zt2 - at3*zt3;
               zt1 = -cimagf(cu[2+3*(j+k1+l1)])
                    + crealf(cu[2+3*(j+k1+l1)])*_Complex_I;
               zt2 = -cimagf(cu[1+3*(j+k1+l1)])
                    + crealf(cu[1+3*(j+k1+l1)])*_Complex_I;
               zt3 = -cimagf(cu[3*(j+k1+l1)])
                    + crealf(cu[3*(j+k1+l1)])*_Complex_I;
               bxyz[3*(j+k1+l1)] = -at3*zt1 + at4*zt2;
               bxyz[1+3*(j+k1+l1)] = -at4*zt3 - at2*zt1;
               bxyz[2+3*(j+k1+l1)] = at2*zt2 + at3*zt3;
               wp += at1*(cu[3*(j+kj+lj)]*conjf(cu[3*(j+kj+lj)])
                  + cu[1+3*(j+kj+lj)]*conjf(cu[1+3*(j+kj+lj)])
                  + cu[2+3*(j+kj+lj)]*conjf(cu[2+3*(j+kj+lj)])
                  + cu[3*(j+k1+lj)]*conjf(cu[3*(j+k1+lj)])
                  + cu[1+3*(j+k1+lj)]*conjf(cu[1+3*(j+k1+lj)])
                  + cu[2+3*(j+k1+lj)]*conjf(cu[2+3*(j+k1+lj)])
                  + cu[3*(j+kj+l1)]*conjf(cu[3*(j+kj+l1)])
                  + cu[1+3*(j+kj+l1)]*conjf(cu[1+3*(j+kj+l1)])
                  + cu[2+3*(j+kj+l1)]*conjf(cu[2+3*(j+kj+l1)])
                  + cu[3*(j+k1+l1)]*conjf(cu[3*(j+k1+l1)])
                  + cu[1+3*(j+k1+l1)]*conjf(cu[1+3*(j+k1+l1)])
                  + cu[2+3*(j+k1+l1)]*conjf(cu[2+3*(j+k1+l1)]));
            }
         }
/* mode numbers kx = 0, nx/2 */
         for (k = 1; k < nyh; k++) {
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            at1 = ci2*crealf(ffc[kk+ll])*cimagf(ffc[kk+ll]);
            at3 = at1*dny*(float) k;
            at4 = dkz*at1;
            zt1 = -cimagf(cu[2+3*(kj+lj)])
                 + crealf(cu[2+3*(kj+lj)])*_Complex_I;
            zt2 = -cimagf(cu[1+3*(kj+lj)])
                 + crealf(cu[1+3*(kj+lj)])*_Complex_I;
            zt3 = -cimagf(cu[3*(kj+lj)])
                 + crealf(cu[3*(kj+lj)])*_Complex_I;
            bxyz[3*(kj+lj)] = at3*zt1 - at4*zt2;
            bxyz[1+3*(kj+lj)] = at4*zt3;
            bxyz[2+3*(kj+lj)] = -at3*zt3;
            bxyz[3*(k1+lj)] = zero;
            bxyz[1+3*(k1+lj)] = zero;
            bxyz[2+3*(k1+lj)] = zero;
            zt1 = -cimagf(cu[2+3*(kj+l1)])
                 + crealf(cu[2+3*(kj+l1)])*_Complex_I;
            zt2 = -cimagf(cu[1+3*(kj+l1)])
                 + crealf(cu[1+3*(kj+l1)])*_Complex_I;
            zt3 = -cimagf(cu[3*(kj+l1)])
                 + crealf(cu[3*(kj+l1)])*_Complex_I;
            bxyz[3*(kj+l1)] = at3*zt1 + at4*zt2;
            bxyz[1+3*(kj+l1)] = -at4*zt3;
            bxyz[2+3*(kj+l1)] = -at3*zt3;
            bxyz[3*(k1+l1)] = zero;
            bxyz[1+3*(k1+l1)]  = zero;
            bxyz[2+3*(k1+l1)] = zero;
            wp += at1*(cu[3*(kj+lj)]*conjf(cu[3*(kj+lj)])
               + cu[1+3*(kj+lj)]*conjf(cu[1+3*(kj+lj)])
               + cu[2+3*(kj+lj)]*conjf(cu[2+3*(kj+lj)])
               + cu[3*(kj+l1)]*conjf(cu[3*(kj+l1)])
               + cu[1+3*(kj+l1)]*conjf(cu[1+3*(kj+l1)])
               + cu[2+3*(kj+l1)]*conjf(cu[2+3*(kj+l1)]));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
         for (j = 1; j < nxh; j++) {
            at1 = ci2*crealf(ffc[j+ll])*cimagf(ffc[j+ll]);
            at2 = at1*dnx*(float) j;  
            at4 = dkz*at1;
            zt1 = -cimagf(cu[2+3*(j+lj)])
                 + crealf(cu[2+3*(j+lj)])*_Complex_I;
            zt2 = -cimagf(cu[1+3*(j+lj)])
                 + crealf(cu[1+3*(j+lj)])*_Complex_I;
            zt3 = -cimagf(cu[3*(j+lj)])
                 + crealf(cu[3*(j+lj)])*_Complex_I;
             bxyz[3*(j+lj)] = -at4*zt2;
            bxyz[1+3*(j+lj)] = at4*zt3 - at2*zt1;
            bxyz[2+3*(j+lj)] = at2*zt2;
            bxyz[3*(j+k1+lj)] = zero;
            bxyz[1+3*(j+k1+lj)] = zero;
            bxyz[2+3*(j+k1+lj)] = zero;
            zt1 = -cimagf(cu[2+3*(j+l1)])
                 + crealf(cu[2+3*(j+l1)])*_Complex_I;
            zt2 = -cimagf(cu[1+3*(j+l1)])
                 + crealf(cu[1+3*(j+l1)])*_Complex_I;
            zt3 = -cimagf(cu[3*(j+l1)])
                 + crealf(cu[3*(j+l1)])*_Complex_I;
            bxyz[3*(j+l1)] = at4*zt2;
            bxyz[1+3*(j+l1)] = -at4*zt3 - at2*zt1;
            bxyz[2+3*(j+l1)] = at2*zt2;
            bxyz[3*(j+k1+l1)] = zero;
            bxyz[1+3*(j+k1+l1)] = zero;
            bxyz[2+3*(j+k1+l1)] = zero;
            wp += at1*(cu[3*(j+lj)]*conjf(cu[3*(j+lj)])
               + cu[1+3*(j+lj)]*conjf(cu[1+3*(j+lj)])
               + cu[2+3*(j+lj)]*conjf(cu[2+3*(j+lj)])
               + cu[3*(j+l1)]*conjf(cu[3*(j+l1)])
               + cu[1+3*(j+l1)]*conjf(cu[1+3*(j+l1)])
               + cu[2+3*(j+l1)]*conjf(cu[2+3*(j+l1)]));
         }
/* mode numbers kx = 0, nx/2 */
         at1 = ci2*crealf(ffc[ll])*cimagf(ffc[ll]);
         at4 = dkz*at1;
         zt2 = -cimagf(cu[1+3*(lj)]) + crealf(cu[1+3*(lj)])*_Complex_I;
         zt3 = -cimagf(cu[3*(lj)]) + crealf(cu[3*(lj)])*_Complex_I;
         bxyz[3*lj] = -at4*zt2;
         bxyz[1+3*lj] = at4*zt3;
         bxyz[2+3*lj] = zero;
         bxyz[3*(k1+lj)] = zero;
         bxyz[1+3*(k1+lj)] = zero;
         bxyz[2+3*(k1+lj)] = zero;
         bxyz[3*l1] = zero;
         bxyz[1+3*l1] = zero;
         bxyz[2+3*l1] = zero;
         bxyz[3*(k1+l1)] = zero;
         bxyz[1+3*(k1+l1)] = zero;
         bxyz[2+3*(k1+l1)] = zero;
         wp += at1*(cu[3*lj]*conjf(cu[3*lj])
            + cu[1+3*lj]*conjf(cu[1+3*lj])
            + cu[2+3*lj]*conjf(cu[2+3*lj]));
         sum1 += wp;
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   sum2 = 0.0;
#pragma omp parallel for \
private(j,k,k1,kk,kj,dky,at1,at2,at3,zt1,zt2,zt3,wp) \
reduction(+:sum2)
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
         zt1 = -cimagf(cu[2+3*(j+kj)])
              + crealf(cu[2+3*(j+kj)])*_Complex_I;
         zt2 = -cimagf(cu[1+3*(j+kj)])
              + crealf(cu[1+3*(j+kj)])*_Complex_I;
         zt3 = -cimagf(cu[3*(j+kj)])
              + crealf(cu[3*(j+kj)])*_Complex_I;
         bxyz[3*(j+kj)] = at3*zt1;
         bxyz[1+3*(j+kj)] = -at2*zt1;
         bxyz[2+3*(j+kj)] = at2*zt2 - at3*zt3;
         zt1 = -cimagf(cu[2+3*(j+k1)])
              + crealf(cu[2+3*(j+k1)])*_Complex_I;
         zt2 = -cimagf(cu[1+3*(j+k1)])
              + crealf(cu[1+3*(j+k1)])*_Complex_I;
         zt3 = -cimagf(cu[3*(j+k1)])
              + crealf(cu[3*(j+k1)])*_Complex_I;
         bxyz[3*(j+k1)] = -at3*zt1;
         bxyz[1+3*(j+k1)] = -at2*zt1;
         bxyz[2+3*(j+k1)] = at2*zt2 + at3*zt3;
         bxyz[3*(j+kj+l1)] = zero;
         bxyz[1+3*(j+kj+l1)] = zero;
         bxyz[2+3*(j+kj+l1)] = zero;
         bxyz[3*(j+k1+l1)] = zero;
         bxyz[1+3*(j+k1+l1)] = zero;
         bxyz[2+3*(j+k1+l1)] = zero;
         wp += at1*(cu[3*(j+kj)]*conjf(cu[3*(j+kj)])
            + cu[1+3*(j+kj)]*conjf(cu[1+3*(j+kj)])
            + cu[2+3*(j+kj)]*conjf(cu[2+3*(j+kj)])
            + cu[3*(j+k1)]*conjf(cu[3*(j+k1)])
            + cu[1+3*(j+k1)]*conjf(cu[1+3*(j+k1)])
            + cu[2+3*(j+k1)]*conjf(cu[2+3*(j+k1)]));
      }
/* mode numbers kx = 0, nx/2 */
      at1 = ci2*crealf(ffc[kk])*cimagf(ffc[kk]);
      at3 = at1*dny*(float) k;
      zt1 = -cimagf(cu[2+3*(kj)]) + crealf(cu[2+3*(kj)])*_Complex_I;
      zt3 = -cimagf(cu[3*(kj)]) + crealf(cu[3*(kj)])*_Complex_I;
      bxyz[3*kj] = at3*zt1;
      bxyz[1+3*kj] = zero;
      bxyz[2+3*kj] = -at3*zt3;
      bxyz[3*k1] = zero;
      bxyz[1+3*k1] = zero;
      bxyz[2+3*k1] = zero;
      bxyz[3*(kj+l1)] = zero;
      bxyz[1+3*(kj+l1)] = zero;
      bxyz[2+3*(kj+l1)] = zero;
      bxyz[3*(k1+l1)] = zero;
      bxyz[1+3*(k1+l1)] = zero;
      bxyz[2+3*(k1+l1)] = zero;
      wp += at1*(cu[3*kj]*conjf(cu[3*kj])
         + cu[1+3*kj]*conjf(cu[1+3*kj])
         + cu[2+3*kj]*conjf(cu[2+3*kj]));
      sum2 += wp;
   }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      at1 = ci2*crealf(ffc[j])*cimagf(ffc[j]);
      at2 = at1*dnx*(float) j;
      zt1 = -cimagf(cu[2+3*j]) + crealf(cu[2+3*j])*_Complex_I;
      zt2 = -cimagf(cu[1+3*j]) + crealf(cu[1+3*j])*_Complex_I;
      bxyz[3*j] = zero;
      bxyz[1+3*j] = -at2*zt1;
      bxyz[2+3*j] = at2*zt2;
      bxyz[3*(j+k1)] = zero;
      bxyz[1+3*(j+k1)] = zero;
      bxyz[2+3*(j+k1)] = zero;
      bxyz[3*(j+l1)] = zero;
      bxyz[1+3*(j+l1)] = zero;
      bxyz[2+3*(j+l1)] = zero;
      bxyz[3*(j+k1+l1)] = zero;
      bxyz[1+3*(j+k1+l1)] = zero;
      bxyz[2+3*(j+k1+l1)] = zero;
      wp += at1*(cu[3*j]*conjf(cu[3*j])
         + cu[1+3*j]*conjf(cu[1+3*j])
         + cu[2+3*j]*conjf(cu[2+3*j]));
   }
   bxyz[0] = zero;
   bxyz[1] = zero;
   bxyz[2] = zero;
   bxyz[3*k1] = zero;
   bxyz[1+3*k1] = zero;
   bxyz[2+3*k1] = zero;
   bxyz[3*l1] = zero;
   bxyz[1+3*l1] = zero;
   bxyz[2+3*l1] = zero;
   bxyz[3*(k1+l1)] = zero;
   bxyz[1+3*(k1+l1)] = zero;
   bxyz[2+3*(k1+l1)] = zero;
   *wm = (sum1 + sum2 + wp)*((float) nx)*((float) ny)*((float) nz);
   return;
}

/*--------------------------------------------------------------------*/
void cbaddext3(float bxyz[], float omx, float omy, float omz, int nx,
               int ny, int nz, int nxe, int nye, int nze) {
/* adds constant to magnetic field for 3d code
   bxy = magnetic field
   omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
   nx/ny/nz = system length in x/y/z direction
   nxe = second dimension of magnetic field array, nxe must be >= nx
   nye = third dimension of magnetic field array, nye must be >= ny
   nze = fourth dimension of magnetic field array, nze must be >= nz
local data                                                 */
   int j, k, l, nxye3, ll;
   nxye3 = 3*nxe*nye;
#pragma omp parallel for private(j,k,l,ll)
   for (l = 0; l < nz; l++) {
      ll = nxye3*l;
      for (k = 0; k < ny; k++) {
         for (j = 0; j < nx; j++) {
            bxyz[3*j+3*nxe*k+ll] += omx;
            bxyz[1+3*j+3*nxe*k+ll] += omy;
            bxyz[2+3*j+3*nxe*k+ll] += omz;
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cmdcuperp3(float complex dcu[], float complex amu[], int nx,
                int ny, int nz, int nxvh, int nyv, int nzv) {
/* this subroutine calculates transverse part of the derivative of
   the current density from the momentum flux
   in 3d with periodic boundary conditions.
   input: all, output: dcu
   approximate flop count is:
   244*nxc*nyc*nzc + 82*(nxc*nyc + nxc*nzc + nyc*nzc)
   and (nx/2)*nyc*nzc divides
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   the derivative of the current is calculated using the equations:
   dcu[kz][ky][kx][0] = -sqrt(-1)*(kx*vx*vx+ky*vx*vy+kz*vx*vz)
   dcu[kz][ky][kx][1] = -sqrt(-1)*(kx*vx*vy+ky*vy*vy+kz*vy*vz)
   dcu[kz][ky][kx][2] = -sqrt(-1)*(kx*vx*vz+ky*vy*vz+kz*vz*vz)
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
   j,k,l = fourier mode numbers, except for
   dcux(kx=pi) = dcuy(kx=pi) = dcuz(kx=pi) = 0,
   dcux(ky=pi) = dcuy(ky=pi) = dcux(ky=pi) = 0,
   dcux(kz=pi) = dcuy(kz=pi) = dcuz(kz=pi) = 0,
   dcux(kx=0,ky=0,kz=0) = dcuy(kx=0,ky=0,kz=0) = dcuz(kx=0,ky=0,kz=0) = 0
   the transverse part is calculated using the equation:
   dcu[kz][ky][kx][0] = dcu[kz][ky][kx][0]
                   - kx*(kx*dcu([kz][ky][kx][0]+ky*dcu[kz][ky][kx][1]
                   + kz*dcu[kz][ky][kx][2)/(kx*kx+ky*ky+kz*kz)
   dcu[kz][ky][kx][1] = dcu(2,kx,ky,kz)
                   - ky*(kx*dcu[kz][ky][kx][0]+ky*dcu[kz][ky][kx][1]
                   + kz*dcu[kz][ky][kx][2)/(kx*kx+ky*ky+kz*kz)
   dcu[kz][ky][kx][2] = dcu(3,kx,ky,kz)
                   - kz*(kx*dcu[kz][ky][kx][0]+ky*dcu[kz][ky][kx][1]
                   + kz*dcu[kz][ky][kx][2)/(kx*kx+ky*ky+kz*kz)
   on output:
   dcu[l][k][j][i] = transverse part of complex derivative of current for
   fourier mode (j,k,l)
   amu[l][k][j][0] = xx component of complex momentum flux
   amu[l][k][j][1] = xy component of complex momentum flux
   amu[l][k][j][2] = xz component of complex momentum flux
   amu[l][k][j][3] = yy component of complex momentum flux
   amu[l][k][j][4] = yz component of complex momentum flux
   amu[l][k][j][5] = zz component of complex momentum flux
   all for fourier mode (j,k,l)
   nx/ny/nz = system length in x/y/z direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nzv = fourth dimension of field arrays, must be >= nz
local data                                                 */
   int nxh, nyh, nzh, j, k, l, k1, l1, kj, lj, nxvyh;
   float dnx, dny, dnz, dkx, dky, dkz, dky2, dkz2, dkyz2, at1;
   float complex zero, zt1, zt2, zt3, zt4, zt5;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   zero = 0.0 + 0.0*_Complex_I;
/* calculate transverse part of current */
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,k1,l1,lj,kj,dkx,dky,dkz,dkz2,dkyz2,at1,zt1,zt2,zt3,zt4, \
zt5)
      for (l = 1; l < nzh; l++) {
         dkz = dnz*(float) l;
         lj = nxvyh*l;
         l1 = nxvyh*nz - lj;
         dkz2 = dkz*dkz;
         for (k = 1; k < nyh; k++) {
            dky = dny*(float) k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            dkyz2 = dky*dky + dkz2;
            for (j = 1; j < nxh; j++) {
               dkx = dnx*(float) j;
               at1 = 1.0/(dkx*dkx + dkyz2);
               zt1 = cimagf(amu[6*(j+kj+lj)])
                   - crealf(amu[6*(j+kj+lj)])*_Complex_I;
               zt2 = cimagf(amu[1+6*(j+kj+lj)])
                   - crealf(amu[1+6*(j+kj+lj)])*_Complex_I;
               zt3 = cimagf(amu[2+6*(j+kj+lj)])
                   - crealf(amu[2+6*(j+kj+lj)])*_Complex_I;
               zt1 = dkx*zt1 + dky*zt2 + dkz*zt3;
               zt4 = cimagf(amu[3+6*(j+kj+lj)])
                   - crealf(amu[3+6*(j+kj+lj)])*_Complex_I;
               zt5 = cimagf(amu[4+6*(j+kj+lj)])
                   - crealf(amu[4+6*(j+kj+lj)])*_Complex_I;
               zt2 = dkx*zt2 + dky*zt4 + dkz*zt5;
               zt4 = cimagf(amu[5+6*(j+kj+lj)])
                   - crealf(amu[5+6*(j+kj+lj)])*_Complex_I;
               zt3 = dkx*zt3 + dky*zt5 + dkz*zt4;
               zt4 = at1*(dkx*zt1 + dky*zt2 + dkz*zt3);
               dcu[3*(j+kj+lj)] = zt1 - dkx*zt4;
               dcu[1+3*(j+kj+lj)] = zt2 - dky*zt4;
               dcu[2+3*(j+kj+lj)] = zt3 - dkz*zt4;
               zt1 = cimagf(amu[6*(j+k1+lj)])
                   - crealf(amu[6*(j+k1+lj)])*_Complex_I;
               zt2 = cimagf(amu[1+6*(j+k1+lj)])
                   - crealf(amu[1+6*(j+k1+lj)])*_Complex_I;
               zt3 = cimagf(amu[2+6*(j+k1+lj)])
                   - crealf(amu[2+6*(j+k1+lj)])*_Complex_I;
               zt1 = dkx*zt1 - dky*zt2 + dkz*zt3;
               zt4 = cimagf(amu[3+6*(j+k1+lj)])
                   - crealf(amu[3+6*(j+k1+lj)])*_Complex_I;
               zt5 = cimagf(amu[4+6*(j+k1+lj)])
                   - crealf(amu[4+6*(j+k1+lj)])*_Complex_I;
               zt2 = dkx*zt2 - dky*zt4 + dkz*zt5;
               zt4 = cimagf(amu[5+6*(j+k1+lj)])
                   - crealf(amu[5+6*(j+k1+lj)])*_Complex_I;
               zt3 = dkx*zt3 - dky*zt5 + dkz*zt4;
               zt4 = at1*(dkx*zt1 - dky*zt2 + dkz*zt3);
               dcu[3*(j+k1+lj)] = zt1 - dkx*zt4;
               dcu[1+3*(j+k1+lj)] = zt2 + dky*zt4;
               dcu[2+3*(j+k1+lj)] = zt3 - dkz*zt4;
               zt1 = cimagf(amu[6*(j+kj+l1)])
                   - crealf(amu[6*(j+kj+l1)])*_Complex_I;
               zt2 = cimagf(amu[1+6*(j+kj+l1)])
                   - crealf(amu[1+6*(j+kj+l1)])*_Complex_I;
               zt3 = cimagf(amu[2+6*(j+kj+l1)])
                   - crealf(amu[2+6*(j+kj+l1)])*_Complex_I;
               zt1 = dkx*zt1 + dky*zt2 - dkz*zt3;
               zt4 = cimagf(amu[3+6*(j+kj+l1)])
                   - crealf(amu[3+6*(j+kj+l1)])*_Complex_I;
               zt5 = cimagf(amu[4+6*(j+kj+l1)])
                   - crealf(amu[4+6*(j+kj+l1)])*_Complex_I;
               zt2 = dkx*zt2 + dky*zt4 - dkz*zt5;
               zt4 = cimagf(amu[5+6*(j+kj+l1)])
                   - crealf(amu[5+6*(j+kj+l1)])*_Complex_I;
               zt3 = dkx*zt3 + dky*zt5 - dkz*zt4;
               zt4 = at1*(dkx*zt1 + dky*zt2 - dkz*zt3);
               dcu[3*(j+kj+l1)] = zt1 - dkx*zt4;
               dcu[1+3*(j+kj+l1)] = zt2 - dky*zt4;
               dcu[2+3*(j+kj+l1)] = zt3 + dkz*zt4;
               zt1 = cimagf(amu[6*(j+k1+l1)])
                  - crealf(amu[6*(j+k1+l1)])*_Complex_I;
               zt2 = cimagf(amu[1+6*(j+k1+l1)])
                   - crealf(amu[1+6*(j+k1+l1)])*_Complex_I;
               zt3 = cimagf(amu[2+6*(j+k1+l1)])
                   - crealf(amu[2+6*(j+k1+l1)])*_Complex_I;
               zt1 = dkx*zt1 - dky*zt2 - dkz*zt3;
               zt4 = cimagf(amu[3+6*(j+k1+l1)])
                   - crealf(amu[3+6*(j+k1+l1)])*_Complex_I;
               zt5 = cimagf(amu[4+6*(j+k1+l1)])
                   - crealf(amu[4+6*(j+k1+l1)])*_Complex_I;
               zt2 = dkx*zt2 - dky*zt4 - dkz*zt5;
               zt4 = cimagf(amu[5+6*(j+k1+l1)])
                   - crealf(amu[5+6*(j+k1+l1)])*_Complex_I;
               zt3 = dkx*zt3 - dky*zt5 - dkz*zt4;
               zt4 = at1*(dkx*zt1 - dky*zt2 - dkz*zt3);
               dcu[3*(j+k1+l1)] = zt1 - dkx*zt4;
               dcu[1+3*(j+k1+l1)] = zt2 + dky*zt4;
               dcu[2+3*(j+k1+l1)] = zt3 + dkz*zt4;
            }
         }
/* mode numbers kx = 0, nx/2 */
         for (k = 1; k < nyh; k++) {
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            dky = dny*(float) k;
            at1 = 1.0/(dky*dky + dkz2);
            zt2 = cimagf(amu[1+6*(kj+lj)])
                - crealf(amu[1+6*(kj+lj)])*_Complex_I;
            zt3 = cimagf(amu[2+6*(kj+lj)])
                - crealf(amu[2+6*(kj+lj)])*_Complex_I;
            zt1 = dky*zt2 + dkz*zt3;
            zt4 = cimagf(amu[3+6*(kj+lj)])
                - crealf(amu[3+6*(kj+lj)])*_Complex_I;
            zt5 = cimagf(amu[4+6*(kj+lj)])
                - crealf(amu[4+6*(kj+lj)])*_Complex_I;
            zt2 = dky*zt4 + dkz*zt5;
            zt4 = cimagf(amu[5+6*(kj+lj)])
                - crealf(amu[5+6*(kj+lj)])*_Complex_I;
            zt3 = dky*zt5 + dkz*zt4;
            zt4 = at1*(dky*zt2 + dkz*zt3);
            dcu[3*(kj+lj)] = zt1;
            dcu[1+3*(kj+lj)] = zt2 - dky*zt4;
            dcu[2+3*(kj+lj)] = zt3 - dkz*zt4;
            dcu[3*(k1+lj)] = zero;
            dcu[1+3*(k1+lj)] = zero;
            dcu[2+3*(k1+lj)] = zero;
            zt2 = cimagf(amu[1+6*(kj+l1)])
                - crealf(amu[1+6*(kj+l1)])*_Complex_I;
            zt3 = cimagf(amu[2+6*(kj+l1)])
                - crealf(amu[2+6*(kj+l1)])*_Complex_I;
            zt1 = dky*zt2 - dkz*zt3;
            zt4 = cimagf(amu[3+6*(kj+l1)])
                - crealf(amu[3+6*(kj+l1)])*_Complex_I;
            zt5 = cimagf(amu[4+6*(kj+l1)])
                - crealf(amu[4+6*(kj+l1)])*_Complex_I;
            zt2 = dky*zt4 - dkz*zt5;
            zt4 = cimagf(amu[5+6*(kj+l1)])
                - crealf(amu[5+6*(kj+l1)])*_Complex_I;
            zt3 = dky*zt5 - dkz*zt4;
            zt4 = at1*(dky*zt2 - dkz*zt3);
            dcu[3*(kj+l1)] = zt1;
            dcu[1+3*(kj+l1)] = zt2 - dky*zt4;
            dcu[2+3*(kj+l1)] = zt3 + dkz*zt4;
            dcu[3*(k1+l1)] = zero;
            dcu[1+3*(k1+l1)] = zero;
            dcu[2+3*(k1+l1)] = zero;
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
         for (j = 1; j < nxh; j++) {
            dkx = dnx*(float) j;
            at1 = 1.0/(dkx*dkx + dkz2);
            zt1 = cimagf(amu[6*(j+lj)])
                - crealf(amu[6*(j+lj)])*_Complex_I;
            zt2 = cimagf(amu[1+6*(j+lj)])
               - crealf(amu[1+6*(j+lj)])*_Complex_I;
            zt3 = cimagf(amu[2+6*(j+lj)])
                - crealf(amu[2+6*(j+lj)])*_Complex_I;
            zt1 = dkx*zt1 + dkz*zt3;
            zt5 = cimagf(amu[4+6*(j+lj)])
                - crealf(amu[4+6*(j+lj)])*_Complex_I;
            zt2 = dkx*zt2 + dkz*zt5;
            zt4 = cimagf(amu[5+6*(j+lj)])
                - crealf(amu[5+6*(j+lj)])*_Complex_I;
            zt3 = dkx*zt3 + dkz*zt4;
            zt4 = at1*(dkx*zt1 + dkz*zt3);
            dcu[3*(j+lj)] = zt1 - dkx*zt4;
            dcu[1+3*(j+lj)] = zt2;
            dcu[2+3*(j+lj)] = zt3 - dkz*zt4;
            dcu[3*(j+k1+lj)] = zero;
            dcu[1+3*(j+k1+lj)] = zero;
            dcu[2+3*(j+k1+lj)] = zero;
            zt1 = cimagf(amu[6*(j+l1)])
                - crealf(amu[6*(j+l1)])*_Complex_I;
            zt2 = cimagf(amu[1+6*(j+l1)])
                - crealf(amu[1+6*(j+l1)])*_Complex_I;
            zt3 = cimagf(amu[2+6*(j+l1)])
                - crealf(amu[2+6*(j+l1)])*_Complex_I;
            zt1 = dkx*zt1 - dkz*zt3;
            zt5 = cimagf(amu[4+6*(j+l1)])
                - crealf(amu[4+6*(j+l1)])*_Complex_I;
            zt2 = dkx*zt2 - dkz*zt5;
            zt4 = cimagf(amu[5+6*(j+l1)])
                - crealf(amu[5+6*(j+l1)])*_Complex_I;
            zt3 = dkx*zt3 - dkz*zt4;
            zt4 = at1*(dkx*zt1 - dkz*zt3);
            dcu[3*(j+l1)] = zt1 - dkx*zt4;
            dcu[1+3*(j+l1)] = zt2;
            dcu[2+3*(j+l1)] = zt3 + dkz*zt4;
            dcu[3*(j+k1+l1)] = zero;
            dcu[1+3*(j+k1+l1)] = zero;
            dcu[2+3*(j+k1+l1)] = zero;
         }
/* mode numbers kx = 0, nx/2 */
         zt3 = cimagf(amu[2+6*(lj)]) - crealf(amu[2+6*(lj)])*_Complex_I;
         zt1 = dkz*zt3;
         zt5 = cimagf(amu[4+6*(lj)]) - crealf(amu[4+6*(lj)])*_Complex_I;
         zt2 = dkz*zt5;
         dcu[3*lj] = zt1;
         dcu[1+3*lj]= zt2;
         dcu[2+3*lj] = zero;
         dcu[3*(k1+lj)] = zero;
         dcu[1+3*(k1+lj)] = zero;
         dcu[2+3*(k1+lj)] = zero;
         dcu[3*l1] = zero;
         dcu[1+3*l1] = zero;
         dcu[2+3*l1] = zero;
         dcu[3*(k1+l1)] = zero;
         dcu[1+3*(k1+l1)] = zero;
         dcu[2+3*(k1+l1)] = zero;
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
#pragma omp parallel for \
private(j,k,k1,kj,dky,dky2,dkx,at1,zt1,zt2,zt3,zt4,zt5)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      dky2 = dky*dky;
      for (j = 1; j < nxh; j++) {
         dkx = dnx*(float) j;
         at1 = 1.0/(dkx*dkx + dky2);
         zt1 = cimagf(amu[6*(j+kj)])
             - crealf(amu[6*(j+kj)])*_Complex_I;
         zt2 = cimagf(amu[1+6*(j+kj)])
             - crealf(amu[1+6*(j+kj)])*_Complex_I;
         zt3 = cimagf(amu[2+6*(j+kj)])
             - crealf(amu[2+6*(j+kj)])*_Complex_I;
         zt1 = dkx*zt1 + dky*zt2;
         zt4 = cimagf(amu[3+6*(j+kj)])
             - crealf(amu[3+6*(j+kj)])*_Complex_I;
         zt5 = cimagf(amu[4+6*(j+kj)])
             - crealf(amu[4+6*(j+kj)])*_Complex_I;
         zt2 = dkx*zt2 + dky*zt4;
         zt3 = dkx*zt3 + dky*zt5;
         zt4 = at1*(dkx*zt1 + dky*zt2);
         dcu[3*(j+kj)] = zt1 - dkx*zt4;
         dcu[1+3*(j+kj)] = zt2 - dky*zt4;
         dcu[2+3*(j+kj)] = zt3;
         dcu[3*(j+kj+l1)] = zero;
         dcu[1+3*(j+kj+l1)] = zero;
         dcu[2+3*(j+kj+l1)] = zero;
         zt1 = cimagf(amu[6*(j+k1)])
             - crealf(amu[6*(j+k1)])*_Complex_I;
         zt2 = cimagf(amu[1+6*(j+k1)])
             - crealf(amu[1+6*(j+k1)])*_Complex_I;
         zt3 = cimagf(amu[2+6*(j+k1)])
             - crealf(amu[2+6*(j+k1)])*_Complex_I;
         zt1 = dkx*zt1 - dky*zt2;
         zt4 = cimagf(amu[3+6*(j+k1)])
             - crealf(amu[3+6*(j+k1)])*_Complex_I;
         zt5 = cimagf(amu[4+6*(j+k1)])
             - crealf(amu[4+6*(j+k1)])*_Complex_I;
         zt2 = dkx*zt2 - dky*zt4;
         zt3 = dkx*zt3 - dky*zt5;
         zt4 = at1*(dkx*zt1 - dky*zt2);
         dcu[3*(j+k1)] = zt1 - dkx*zt4;
         dcu[1+3*(j+k1)] = zt2 + dky*zt4;
         dcu[2+3*(j+k1)] = zt3;
         dcu[3*(j+k1+l1)] = zero;
         dcu[1+3*(j+k1+l1)] = zero;
         dcu[2+3*(j+k1+l1)]= zero;
      }
/* mode numbers kx = 0, nx/2 */
      zt2 = cimagf(amu[1+6*(kj)]) - crealf(amu[1+6*(kj)])*_Complex_I;
      zt1 = dky*zt2;
      zt5 = cimagf(amu[4+6*(kj)]) - crealf(amu[4+6*(kj)])*_Complex_I;
      zt3 = dky*zt5;
      dcu[3*kj] = zt1;
      dcu[1+3*kj] = zero;
      dcu[2+3*kj] = zt3;
      dcu[3*k1] = zero;
      dcu[1+3*k1] = zero;
      dcu[2+3*k1] = zero;
      dcu[3*(kj+l1)] = zero;
      dcu[1+3*(kj+l1)] = zero;
      dcu[2+3*(kj+l1)] = zero;
      dcu[3*(k1+l1)] = zero;
      dcu[1+3*(k1+l1)] = zero;
      dcu[2+3*(k1+l1)] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      dkx = dnx*(float) j;
      zt2 = cimagf(amu[1+6*j]) - crealf(amu[1+6*j])*_Complex_I;
      zt3 = cimagf(amu[2+6*j]) - crealf(amu[2+6*j])*_Complex_I;
      zt2 = dkx*zt2;
      zt3 = dkx*zt3;
      dcu[3*j] = zero;
      dcu[1+3*j] = zt2;
      dcu[2+3*j] = zt3;
      dcu[3*(j+k1)] = zero;
      dcu[1+3*(j+k1)] = zero;
      dcu[2+3*(j+k1)] = zero;
      dcu[3*(j+l1)] = zero;
      dcu[1+3*(j+l1)] = zero;
      dcu[2+3*(j+l1)] = zero;
      dcu[3*(j+k1+l1)] = zero;
      dcu[1+3*(j+k1+l1)] = zero;
      dcu[2+3*(j+k1+l1)] = zero;
   }
   dcu[0] = zero;
   dcu[1] = zero;
   dcu[2] = zero;
   dcu[3*k1] = zero;
   dcu[1+3*k1] = zero;
   dcu[2+3*k1] = zero;
   dcu[3*l1] = zero;
   dcu[1+3*l1] = zero;
   dcu[2+3*l1] = zero;
   dcu[3*(k1+l1)] = zero;
   dcu[1+3*(k1+l1)] = zero;
   dcu[2+3*(k1+l1)] = zero;
   return;
}

/*--------------------------------------------------------------------*/
void cmadcuperp3(float complex dcu[], float complex amu[], int nx,
                 int ny, int nz, int nxvh, int nyv, int nzv) {
/* this subroutine calculates transverse part of the derivative of
   the current density from the momentum flux and acceleration density
   in 3d with periodic boundary conditions.
   input: all, output: dcu
   approximate flop count is:
   244*nxc*nyc*nzc + 82*(nxc*nyc + nxc*nzc + nyc*nzc)
   and (nx/2)*nyc*nzc divides
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   the derivative of the current is calculated using the equations:
   dcu[kz][ky][kx][0] = dcu[kz][ky][kx][0]
                  -sqrt(-1)*(kx*vx*vx+ky*vx*vy+kz*vx*vz)
   dcu[kz][ky][kx][1] = dcu[kz][ky][kx][1]
                  -sqrt(-1)*(kx*vx*vy+ky*vy*vy+kz*vy*vz)
   dcu[kz][ky][kx][2] = dcu[kz][ky][kx][2]
                   -sqrt(-1)*(kx*vx*vz+ky*vy*vz+kz*vz*vz)
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
   j,k,l = fourier mode numbers, except for
   dcux(kx=pi) = dcuy(kx=pi) = dcuz(kx=pi) = 0,
   dcux(ky=pi) = dcuy(ky=pi) = dcux(ky=pi) = 0,
   dcux(kz=pi) = dcuy(kz=pi) = dcuz(kz=pi) = 0,
   dcux(kx=0,ky=0,kz=0) = dcuy(kx=0,ky=0,kz=0) = dcuz(kx=0,ky=0,kz=0) = 0
   the transverse part is calculated using the equation:
   dcu[kz][ky][kx][0] = dcu[kz][ky][kx][0]
                   - kx*(kx*dcu[kz][ky][kx][0]+ky*dcu[kz][ky][kx][1]
                   + kz*dcu[kz][ky][kx][2])/(kx*kx+ky*ky+kz*kz)
   dcu[kz][ky][kx][1] = dcu(2,kx,ky,kz)
                   - ky*(kx*dcu[kz][ky][kx][0]+ky*dcu[kz][ky][kx][1]
                   + kz*dcu[kz][ky][kx][2])/(kx*kx+ky*ky+kz*kz)
   dcu[kz][ky][kx][2] = dcu(3,kx,ky,kz)
                   - kz*(kx*dcu[kz][ky][kx][0]+ky*dcu[kz][ky][kx][1]
                   + kz*dcu[kz][ky][kx][2])/(kx*kx+ky*ky+kz*kz)
   on input:
   dcu[l][k][j][i] = complex acceleration density for fourier mode (j-1,k-1)
   on output:
   dcu[l][k][j][i] = transverse part of complex derivative of current for
   fourier mode (j,k,l)
   amu[l][k][j][0] = xx component of complex momentum flux
   amu[l][k][j][1] = xy component of complex momentum flux
   amu[l][k][j][2] = xz component of complex momentum flux
   amu[l][k][j][3] = yy component of complex momentum flux
   amu[l][k][j][4] = yz component of complex momentum flux
   amu[l][k][j][5] = zz component of complex momentum flux
   all for fourier mode (j,k,l)
   nx/ny/nz = system length in x/y/z direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nzv = fourth dimension of field arrays, must be >= nz
local data                                                 */
   int nxh, nyh, nzh, j, k, l, k1, l1, kj, lj, nxvyh;
   float dnx, dny, dnz, dkx, dky, dkz, dky2, dkz2, dkyz2, at1;
   float complex zero, zt1, zt2, zt3, zt4, zt5;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   zero = 0.0 + 0.0*_Complex_I;
/* calculate transverse part of current */
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,k1,l1,lj,kj,dkx,dky,dkz,dkz2,dkyz2,at1,zt1,zt2,zt3,zt4, \
zt5)
      for (l = 1; l < nzh; l++) {
         dkz = dnz*(float) l;
         lj = nxvyh*l;
         l1 = nxvyh*nz - lj;
         dkz2 = dkz*dkz;
         for (k = 1; k < nyh; k++) {
            dky = dny*(float) k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            dkyz2 = dky*dky + dkz2;
            for (j = 1; j < nxh; j++) {
               dkx = dnx*(float) j;
               at1 = 1.0/(dkx*dkx + dkyz2);
               zt1 = cimagf(amu[6*(j+kj+lj)])
                   - crealf(amu[6*(j+kj+lj)])*_Complex_I;
               zt2 = cimagf(amu[1+6*(j+kj+lj)])
                   - crealf(amu[1+6*(j+kj+lj)])*_Complex_I;
               zt3 = cimagf(amu[2+6*(j+kj+lj)])
                   - crealf(amu[2+6*(j+kj+lj)])*_Complex_I;
               zt1 = dcu[3*(j+kj+lj)] + dkx*zt1 + dky*zt2 + dkz*zt3;
               zt4 = cimagf(amu[3+6*(j+kj+lj)])
                   - crealf(amu[3+6*(j+kj+lj)])*_Complex_I;
               zt5 = cimagf(amu[4+6*(j+kj+lj)])
                   - crealf(amu[4+6*(j+kj+lj)])*_Complex_I;
               zt2 = dcu[1+3*(j+kj+lj)] + dkx*zt2 + dky*zt4 + dkz*zt5;
               zt4 = cimagf(amu[5+6*(j+kj+lj)])
                   - crealf(amu[5+6*(j+kj+lj)])*_Complex_I;
               zt3 = dcu[2+3*(j+kj+lj)] + dkx*zt3 + dky*zt5 + dkz*zt4;
               zt4 = at1*(dkx*zt1 + dky*zt2 + dkz*zt3);
               dcu[3*(j+kj+lj)] = zt1 - dkx*zt4;
               dcu[1+3*(j+kj+lj)] = zt2 - dky*zt4;
               dcu[2+3*(j+kj+lj)] = zt3 - dkz*zt4;
               zt1 = cimagf(amu[6*(j+k1+lj)])
                   - crealf(amu[6*(j+k1+lj)])*_Complex_I;
               zt2 = cimagf(amu[1+6*(j+k1+lj)])
                   - crealf(amu[1+6*(j+k1+lj)])*_Complex_I;
               zt3 = cimagf(amu[2+6*(j+k1+lj)])
                   - crealf(amu[2+6*(j+k1+lj)])*_Complex_I;
               zt1 = dcu[3*(j+k1+lj)] + dkx*zt1 - dky*zt2 + dkz*zt3;
               zt4 = cimagf(amu[3+6*(j+k1+lj)])
                   - crealf(amu[3+6*(j+k1+lj)])*_Complex_I;
               zt5 = cimagf(amu[4+6*(j+k1+lj)])
                   - crealf(amu[4+6*(j+k1+lj)])*_Complex_I;
               zt2 = dcu[1+3*(j+k1+lj)] + dkx*zt2 - dky*zt4 + dkz*zt5;
               zt4 = cimagf(amu[5+6*(j+k1+lj)])
                   - crealf(amu[5+6*(j+k1+lj)])*_Complex_I;
               zt3 = dcu[2+3*(j+k1+lj)] + dkx*zt3 - dky*zt5 + dkz*zt4;
               zt4 = at1*(dkx*zt1 - dky*zt2 + dkz*zt3);
               dcu[3*(j+k1+lj)] = zt1 - dkx*zt4;
               dcu[1+3*(j+k1+lj)] = zt2 + dky*zt4;
               dcu[2+3*(j+k1+lj)] = zt3 - dkz*zt4;
               zt1 = cimagf(amu[6*(j+kj+l1)])
                   - crealf(amu[6*(j+kj+l1)])*_Complex_I;
               zt2 = cimagf(amu[1+6*(j+kj+l1)])
                   - crealf(amu[1+6*(j+kj+l1)])*_Complex_I;
               zt3 = cimagf(amu[2+6*(j+kj+l1)])
                   - crealf(amu[2+6*(j+kj+l1)])*_Complex_I;
               zt1 = dcu[3*(j+kj+l1)] + dkx*zt1 + dky*zt2 - dkz*zt3;
               zt4 = cimagf(amu[3+6*(j+kj+l1)])
                   - crealf(amu[3+6*(j+kj+l1)])*_Complex_I;
               zt5 = cimagf(amu[4+6*(j+kj+l1)])
                   - crealf(amu[4+6*(j+kj+l1)])*_Complex_I;
               zt2 = dcu[1+3*(j+kj+l1)] + dkx*zt2 + dky*zt4 - dkz*zt5;
               zt4 = cimagf(amu[5+6*(j+kj+l1)])
                   - crealf(amu[5+6*(j+kj+l1)])*_Complex_I;
               zt3 = dcu[2+3*(j+kj+l1)] + dkx*zt3 + dky*zt5 - dkz*zt4;
               zt4 = at1*(dkx*zt1 + dky*zt2 - dkz*zt3);
               dcu[3*(j+kj+l1)] = zt1 - dkx*zt4;
               dcu[1+3*(j+kj+l1)] = zt2 - dky*zt4;
               dcu[2+3*(j+kj+l1)] = zt3 + dkz*zt4;
               zt1 = cimagf(amu[6*(j+k1+l1)])
                   - crealf(amu[6*(j+k1+l1)])*_Complex_I;
               zt2 = cimagf(amu[1+6*(j+k1+l1)])
                   - crealf(amu[1+6*(j+k1+l1)])*_Complex_I;
               zt3 = cimagf(amu[2+6*(j+k1+l1)])
                   - crealf(amu[2+6*(j+k1+l1)])*_Complex_I;
               zt1 = dcu[3*(j+k1+l1)] + dkx*zt1 - dky*zt2 - dkz*zt3;
               zt4 = cimagf(amu[3+6*(j+k1+l1)])
                   - crealf(amu[3+6*(j+k1+l1)])*_Complex_I;
               zt5 = cimagf(amu[4+6*(j+k1+l1)])
                   - crealf(amu[4+6*(j+k1+l1)])*_Complex_I;
               zt2 = dcu[1+3*(j+k1+l1)] + dkx*zt2 - dky*zt4 - dkz*zt5;
               zt4 = cimagf(amu[5+6*(j+k1+l1)])
                   - crealf(amu[5+6*(j+k1+l1)])*_Complex_I;
               zt3 = dcu[2+3*(j+k1+l1)] + dkx*zt3 - dky*zt5 - dkz*zt4;
               zt4 = at1*(dkx*zt1 - dky*zt2 - dkz*zt3);
               dcu[3*(j+k1+l1)] = zt1 - dkx*zt4;
               dcu[1+3*(j+k1+l1)] = zt2 + dky*zt4;
               dcu[2+3*(j+k1+l1)] = zt3 + dkz*zt4;
            }
         }
/* mode numbers kx = 0, nx/2 */
         for (k = 1; k < nyh; k++) {
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            dky = dny*(float) k;
            at1 = 1.0/(dky*dky + dkz2);
            zt2 = cimagf(amu[1+6*(kj+lj)])
                - crealf(amu[1+6*(kj+lj)])*_Complex_I;
            zt3 = cimagf(amu[2+6*(kj+lj)])
                - crealf(amu[2+6*(kj+lj)])*_Complex_I;
            zt1 = dcu[3*(kj+lj)] + dky*zt2 + dkz*zt3;
            zt4 = cimagf(amu[3+6*(kj+lj)])
                - crealf(amu[3+6*(kj+lj)])*_Complex_I;
            zt5 = cimagf(amu[4+6*(kj+lj)])
                - crealf(amu[4+6*(kj+lj)])*_Complex_I;
            zt2 = dcu[1+3*(kj+lj)] + dky*zt4 + dkz*zt5;
            zt4 = cimagf(amu[5+6*(kj+lj)])
                - crealf(amu[5+6*(kj+lj)])*_Complex_I;
            zt3 = dcu[2+3*(kj+lj)] + dky*zt5 + dkz*zt4;
            zt4 = at1*(dky*zt2 + dkz*zt3);
            dcu[3*(kj+lj)] = zt1;
            dcu[1+3*(kj+lj)] = zt2 - dky*zt4;
            dcu[2+3*(kj+lj)] = zt3 - dkz*zt4;
            dcu[3*(k1+lj)] = zero;
            dcu[1+3*(k1+lj)] = zero;
            dcu[2+3*(k1+lj)] = zero;
            zt2 = cimagf(amu[1+6*(kj+l1)])
                - crealf(amu[1+6*(kj+l1)])*_Complex_I;
            zt3 = cimagf(amu[2+6*(kj+l1)])
                - crealf(amu[2+6*(kj+l1)])*_Complex_I;
            zt1 = dcu[3*(kj+l1)] + dky*zt2 - dkz*zt3;
            zt4 = cimagf(amu[3+6*(kj+l1)])
                - crealf(amu[3+6*(kj+l1)])*_Complex_I;
            zt5 = cimagf(amu[4+6*(kj+l1)])
                - crealf(amu[4+6*(kj+l1)])*_Complex_I;
            zt2 = dcu[1+3*(kj+l1)] + dky*zt4 - dkz*zt5;
            zt4 = cimagf(amu[5+6*(kj+l1)])
                - crealf(amu[5+6*(kj+l1)])*_Complex_I;
            zt3 = dcu[2+3*(kj+l1)] + dky*zt5 - dkz*zt4;
            zt4 = at1*(dky*zt2 - dkz*zt3);
            dcu[3*(kj+l1)] = zt1;
            dcu[1+3*(kj+l1)] = zt2 - dky*zt4;
            dcu[2+3*(kj+l1)] = zt3 + dkz*zt4;
            dcu[3*(k1+l1)] = zero;
            dcu[1+3*(k1+l1)] = zero;
            dcu[2+3*(k1+l1)] = zero;
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
         for (j = 1; j < nxh; j++) {
            dkx = dnx*(float) j;
            at1 = 1.0/(dkx*dkx + dkz2);
            zt1 = cimagf(amu[6*(j+lj)])
                - crealf(amu[6*(j+lj)])*_Complex_I;
            zt2 = cimagf(amu[1+6*(j+lj)])
                - crealf(amu[1+6*(j+lj)])*_Complex_I;
            zt3 = cimagf(amu[2+6*(j+lj)])
                - crealf(amu[2+6*(j+lj)])*_Complex_I;
            zt1 = dcu[3*(j+lj)] + dkx*zt1 + dkz*zt3;
            zt5 = cimagf(amu[4+6*(j+lj)])
                - crealf(amu[4+6*(j+lj)])*_Complex_I;
            zt2 = dcu[1+3*(j+lj)] + dkx*zt2 + dkz*zt5;
            zt4 = cimagf(amu[5+6*(j+lj)])
                - crealf(amu[5+6*(j+lj)])*_Complex_I;
            zt3 = dcu[2+3*(j+lj)] + dkx*zt3 + dkz*zt4;
            zt4 = at1*(dkx*zt1 + dkz*zt3);
            dcu[3*(j+lj)] = zt1 - dkx*zt4;
            dcu[1+3*(j+lj)] = zt2;
            dcu[2+3*(j+lj)] = zt3 - dkz*zt4;
            dcu[3*(j+k1+lj)] = zero;
            dcu[1+3*(j+k1+lj)] = zero;
            dcu[2+3*(j+k1+lj)] = zero;
            zt1 = cimagf(amu[6*(j+l1)])
                - crealf(amu[6*(j+l1)])*_Complex_I;
            zt2 = cimagf(amu[1+6*(j+l1)])
                - crealf(amu[1+6*(j+l1)])*_Complex_I;
            zt3 = cimagf(amu[2+6*(j+l1)])
                - crealf(amu[2+6*(j+l1)])*_Complex_I;
            zt1 = dcu[3*(j+l1)] + dkx*zt1 - dkz*zt3;
            zt5 = cimagf(amu[4+6*(j+l1)])
                - crealf(amu[4+6*(j+l1)])*_Complex_I;
            zt2 = dcu[1+3*(j+l1)] + dkx*zt2 - dkz*zt5;
            zt4 = cimagf(amu[5+6*(j+l1)])
                - crealf(amu[5+6*(j+l1)])*_Complex_I;
            zt3 = dcu[2+3*(j+l1)] + dkx*zt3 - dkz*zt4;
            zt4 = at1*(dkx*zt1 - dkz*zt3);
            dcu[3*(j+l1)] = zt1 - dkx*zt4;
            dcu[1+3*(j+l1)] = zt2;
            dcu[2+3*(j+l1)] = zt3 + dkz*zt4;
            dcu[3*(j+k1+l1)] = zero;
            dcu[1+3*(j+k1+l1)] = zero;
            dcu[2+3*(j+k1+l1)] = zero;
         }
/* mode numbers kx = 0, nx/2 */
         zt3 = cimagf(amu[2+6*(lj)]) - crealf(amu[2+6*(lj)])*_Complex_I;
         zt1 = dcu[3*lj] + dkz*zt3;
         zt5 = cimagf(amu[4+6*(lj)]) - crealf(amu[4+6*(lj)])*_Complex_I;
         zt2 = dcu[1+3*lj] + dkz*zt5;
         dcu[3*lj] = zt1;
         dcu[1+3*lj]= zt2;
         dcu[2+3*lj] = zero;
         dcu[3*(k1+lj)] = zero;
         dcu[1+3*(k1+lj)] = zero;
         dcu[2+3*(k1+lj)] = zero;
         dcu[3*l1] = zero;
         dcu[1+3*l1] = zero;
         dcu[2+3*l1] = zero;
         dcu[3*(k1+l1)] = zero;
         dcu[1+3*(k1+l1)] = zero;
         dcu[2+3*(k1+l1)] = zero;
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
#pragma omp parallel for \
private(j,k,k1,kj,dky,dky2,dkx,at1,zt1,zt2,zt3,zt4,zt5)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      dky2 = dky*dky;
      for (j = 1; j < nxh; j++) {
         dkx = dnx*(float) j;
         at1 = 1.0/(dkx*dkx + dky2);
         zt1 = cimagf(amu[6*(j+kj)])
             - crealf(amu[6*(j+kj)])*_Complex_I;
         zt2 = cimagf(amu[1+6*(j+kj)])
             - crealf(amu[1+6*(j+kj)])*_Complex_I;
         zt3 = cimagf(amu[2+6*(j+kj)])
             - crealf(amu[2+6*(j+kj)])*_Complex_I;
         zt1 = dcu[3*(j+kj)] + dkx*zt1 + dky*zt2;
         zt4 = cimagf(amu[3+6*(j+kj)])
             - crealf(amu[3+6*(j+kj)])*_Complex_I;
         zt5 = cimagf(amu[4+6*(j+kj)])
             - crealf(amu[4+6*(j+kj)])*_Complex_I;
         zt2 = dcu[1+3*(j+kj)] + dkx*zt2 + dky*zt4;
         zt3 = dcu[2+3*(j+kj)] + dkx*zt3 + dky*zt5;
         zt4 = at1*(dkx*zt1 + dky*zt2);
         dcu[3*(j+kj)] = zt1 - dkx*zt4;
         dcu[1+3*(j+kj)] = zt2 - dky*zt4;
         dcu[2+3*(j+kj)] = zt3;
         dcu[3*(j+kj+l1)] = zero;
         dcu[1+3*(j+kj+l1)] = zero;
         dcu[2+3*(j+kj+l1)] = zero;
         zt1 = cimagf(amu[6*(j+k1)])
             - crealf(amu[6*(j+k1)])*_Complex_I;
         zt2 = cimagf(amu[1+6*(j+k1)])
             - crealf(amu[1+6*(j+k1)])*_Complex_I;
         zt3 = cimagf(amu[2+6*(j+k1)])
             - crealf(amu[2+6*(j+k1)])*_Complex_I;
         zt1 = dcu[3*(j+k1)] + dkx*zt1 - dky*zt2;
         zt4 = cimagf(amu[3+6*(j+k1)])
             - crealf(amu[3+6*(j+k1)])*_Complex_I;
         zt5 = cimagf(amu[4+6*(j+k1)])
             - crealf(amu[4+6*(j+k1)])*_Complex_I;
         zt2 = dcu[1+3*(j+k1)] + dkx*zt2 - dky*zt4;
         zt3 = dcu[2+3*(j+k1)] + dkx*zt3 - dky*zt5;
         zt4 = at1*(dkx*zt1 - dky*zt2);
         dcu[3*(j+k1)] = zt1 - dkx*zt4;
         dcu[1+3*(j+k1)] = zt2 + dky*zt4;
         dcu[2+3*(j+k1)] = zt3;
         dcu[3*(j+k1+l1)] = zero;
         dcu[1+3*(j+k1+l1)] = zero;
         dcu[2+3*(j+k1+l1)]= zero;
      }
      zt2 = cimagf(amu[1+6*(kj)]) - crealf(amu[1+6*(kj)])*_Complex_I;
      zt1 = dcu[3*kj] + dky*zt2;
      zt5 = cimagf(amu[4+6*(kj)]) - crealf(amu[4+6*(kj)])*_Complex_I;
      zt3 = dcu[2+3*kj] + dky*zt5;
      dcu[3*kj] = zt1;
      dcu[1+3*kj] = zero;
      dcu[2+3*kj] = zt3;
      dcu[3*k1] = zero;
      dcu[1+3*k1] = zero;
      dcu[2+3*k1] = zero;
      dcu[3*(kj+l1)] = zero;
      dcu[1+3*(kj+l1)] = zero;
      dcu[2+3*(kj+l1)] = zero;
      dcu[3*(k1+l1)] = zero;
      dcu[1+3*(k1+l1)] = zero;
      dcu[2+3*(k1+l1)] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      dkx = dnx*(float) j;
      zt2 = cimagf(amu[1+6*j]) - crealf(amu[1+6*j])*_Complex_I;
      zt3 = cimagf(amu[2+6*j]) - crealf(amu[2+6*j])*_Complex_I;
      zt2 = dcu[1+3*j] + dkx*zt2;
      zt3 = dcu[2+3*j] + dkx*zt3;
      dcu[3*j] = zero;
      dcu[1+3*j] = zt2;
      dcu[2+3*j] = zt3;
      dcu[3*(j+k1)] = zero;
      dcu[1+3*(j+k1)] = zero;
      dcu[2+3*(j+k1)] = zero;
      dcu[3*(j+l1)] = zero;
      dcu[1+3*(j+l1)] = zero;
      dcu[2+3*(j+l1)] = zero;
      dcu[3*(j+k1+l1)] = zero;
      dcu[1+3*(j+k1+l1)] = zero;
      dcu[2+3*(j+k1+l1)] = zero;
   }
   dcu[0] = zero;
   dcu[1] = zero;
   dcu[2] = zero;
   dcu[3*k1] = zero;
   dcu[1+3*k1] = zero;
   dcu[2+3*k1] = zero;
   dcu[3*l1] = zero;
   dcu[1+3*l1] = zero;
   dcu[2+3*l1] = zero;
   dcu[3*(k1+l1)] = zero;
   dcu[1+3*(k1+l1)] = zero;
   dcu[2+3*(k1+l1)] = zero;
   return;
}

/*--------------------------------------------------------------------*/
void cmepois33(float complex dcu[], float complex exyz[], int isign,
               float complex ffe[], float ax, float ay, float az,
               float affp, float wp0, float ci, float *wf, int nx,
               int ny, int nz, int nxvh, int nyv, int nzv, int nxhd,
               int nyhd, int nzhd) {
/* this subroutine solves 3d poisson's equation in fourier space for
   transverse electric field (or convolution of transverse electric field
   over particle shape), with periodic boundary conditions.
   using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
   A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
   for isign = 0, output: ffe
   input: isign,ax,ay,az,affp,wp0,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
   for isign /= 0, output: exyz, wf
   input: dcu,ffe,isign,ci,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
   approximate flop count is:
   128*nxc*nyc*nzc + 66*(nxc*nyc + nxc*nzc + nyc*nzc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   if isign = -1, smoothed transverse electric field is calculated
   using the equation:
   ex[kz][ky][kx] = -ci*ci*g[kz][ky][kx]*dcux[kz][ky][kx]*s[kz][ky][kx]
   ey[kz][ky][kx] = -ci*ci*g[kz][ky][kx]*dcuy[kz][ky][kx]*s[kz][ky][kx]
   ez[kz][ky][kx] = -ci*ci*g[kz][ky][kx]*dcuz[kz][ky][kx]*s[kz][ky][kx]
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
   j,k,l = fourier mode numbers,
   g[kz][ky][kx] = (affp/(kx**2+ky**2+kz**2))*s[kz][ky][kx],
   s[kz][ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
   ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
   ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
   ex(kz=pi) = ey(kz=pi) = ez(kz=pi) = 0,
   ex(kx=0,ky=0,kz=0) = ey(kx=0,ky=0,kz=0) = ez(kx=0,ky=0,kz=0) = 0.
   if isign = 1, unsmoothed transverse electric field is calculated
   using the equation:
   ex[kz][ky][kx] = -ci*ci*g[kz][ky][kx]*dcux[kz][ky][kx]
   ey[kz][ky][kx] = -ci*ci*g[kz][ky][kx])*dcuy[kz][ky][kx]
   ez[kz][ky][kx] = -ci*ci*g[kz][ky][kx]*dcuz[kz][ky][kx]
   dcu[l][k][j][i] = transverse part of complex derivative of current for
   fourier mode (j,k,l)
   exyz[l][k][j][0] = x component of complex transverse electric field
   exyz[l][k][j][1] = y component of complex transverse electric field
   exyz[l][k][j][2] = z component of complex transverse electric field
   all for fourier mode (j,k,l)
   cimag(ffe[l][k][j]) = finite-size particle shape factor s
   for fourier mode (j,k,l)
   creal(ffe[l][k][j]) = potential green's function g
   for fourier mode (j,k,l)
   ax/ay/az = half-width of particle in x/y/z direction
   affp = normalization constant = nx*ny*nz/np,
   where np=number of particles
   wp0 = normalized total plasma frequency squared
   where np=number of particles
   ci = reciprocal of velocity of light
   transverse electric field energy is also calculated, using
   wf = nx*ny*nz*sum((affp/((kx**2+ky**2+kz**2)*ci*ci)**2)
      |dcu[kz][ky][kx]*s[kz][ky][kx]|**2)
   this expression is valid only if the derivative of current is
   divergence-free
   nx/ny/nz = system length in x/y/z direction
   nxvh = first dimension of field arrays, must be >= nxh
   nyv = second dimension of field arrays, must be >= ny
   nzv = third dimension of field arrays, must be >= nz
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
   nzhd = third dimension of form factor array, must be >= nzh
local data                                                 */
   int nxh, nyh, nzh, j, k, l, k1, l1, kk, kj, ll, lj, nxyhd, nxvyh;
   float dnx, dny, dnz, dkx, dky, dkz, ci2, wpc;
   float at1, at2, at3, at4, at5, at6;
   float complex zero;
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
   ci2 = ci*ci;
   if (isign != 0)
      goto L40;
   wpc = wp0*ci2;
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
               ffe[j+kk+ll] = affp + 1.0*_Complex_I;
            }
            else {
               ffe[j+kk+ll] = (affp*at6/(at5 + wpc*at6*at6))
                            + at6*_Complex_I;
            }
         }
      }
   }
   return;
L40: if (isign > 0)
      goto L130;
/* calculate smoothed transverse electric field and sum field energy */
   sum1 = 0.0;
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,k1,l1,ll,lj,kk,kj,at1,at2,wp) \
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
            for (j = 1; j < nxh; j++) {
               at2 = -ci2*crealf(ffe[j+kk+ll]);
               at1 = at2*cimagf(ffe[j+kk+ll]);
               at2 = at2*at2;
               exyz[3*(j+kj+lj)] = at1*dcu[3*(j+kj+lj)];
               exyz[1+3*(j+kj+lj)] = at1*dcu[1+3*(j+kj+lj)];
               exyz[2+3*(j+kj+lj)] = at1*dcu[2+3*(j+kj+lj)];
               exyz[3*(j+k1+lj)] = at1*dcu[3*(j+k1+lj)];
               exyz[1+3*(j+k1+lj)] = at1*dcu[1+3*(j+k1+lj)];
               exyz[2+3*(j+k1+lj)] = at1*dcu[2+3*(j+k1+lj)];
               exyz[3*(j+kj+l1)] = at1*dcu[3*(j+kj+l1)];
               exyz[1+3*(j+kj+l1)] = at1*dcu[1+3*(j+kj+l1)];
               exyz[2+3*(j+kj+l1)] = at1*dcu[2+3*(j+kj+l1)];
               exyz[3*(j+k1+l1)] = at1*dcu[3*(j+k1+l1)];
               exyz[1+3*(j+k1+l1)] = at1*dcu[1+3*(j+k1+l1)];
               exyz[2+3*(j+k1+l1)] = at1*dcu[2+3*(j+k1+l1)];
               wp += at2*(dcu[3*(j+kj+lj)]*conjf(dcu[3*(j+kj+lj)])
                  + dcu[1+3*(j+kj+lj)]*conjf(dcu[1+3*(j+kj+lj)])
                  + dcu[2+3*(j+kj+lj)]*conjf(dcu[2+3*(j+kj+lj)])
                  + dcu[3*(j+k1+lj)]*conjf(dcu[3*(j+k1+lj)])
                  + dcu[1+3*(j+k1+lj)]*conjf(dcu[1+3*(j+k1+lj)])
                  + dcu[2+3*(j+k1+lj)]*conjf(dcu[2+3*(j+k1+lj)])
                  + dcu[3*(j+kj+l1)]*conjf(dcu[3*(j+kj+l1)])
                  + dcu[1+3*(j+kj+l1)]*conjf(dcu[1+3*(j+kj+l1)])
                  + dcu[2+3*(j+kj+l1)]*conjf(dcu[2+3*(j+kj+l1)])
                  + dcu[3*(j+k1+l1)]*conjf(dcu[3*(j+k1+l1)])
                  + dcu[1+3*(j+k1+l1)]*conjf(dcu[1+3*(j+k1+l1)])
                  + dcu[2+3*(j+k1+l1)]*conjf(dcu[2+3*(j+k1+l1)]));
            }
         }
/* mode numbers kx = 0, nx/2 */
         for (k = 1; k < nyh; k++) {
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            at2 = -ci2*crealf(ffe[kk+ll]);
            at1 = at2*cimagf(ffe[kk+ll]);
            at2 = at2*at2;
            exyz[3*(kj+lj)] = at1*dcu[3*(kj+lj)];
            exyz[1+3*(kj+lj)] = at1*dcu[1+3*(kj+lj)];
            exyz[2+3*(kj+lj)] = at1*dcu[2+3*(kj+lj)];
            exyz[3*(k1+lj)] = zero;
            exyz[1+3*(k1+lj)] = zero;
            exyz[2+3*(k1+lj)] = zero;
            exyz[3*(kj+l1)] = at1*dcu[3*(kj+l1)];
            exyz[1+3*(kj+l1)] = at1*dcu[1+3*(kj+l1)];
            exyz[2+3*(kj+l1)] = at1*dcu[2+3*(kj+l1)];
            exyz[3*(k1+l1)] = zero;
            exyz[1+3*(k1+l1)] = zero;
            exyz[2+3*(k1+l1)] = zero;
            wp += at2*(dcu[3*(kj+lj)]*conjf(dcu[3*(kj+lj)])
               + dcu[1+3*(kj+lj)]*conjf(dcu[1+3*(kj+lj)])
               + dcu[2+3*(kj+lj)]*conjf(dcu[2+3*(kj+lj)])
               + dcu[3*(kj+l1)]*conjf(dcu[3*(kj+l1)])
               + dcu[1+3*(kj+l1)]*conjf(dcu[1+3*(kj+l1)])
               + dcu[2+3*(kj+l1)]*conjf(dcu[2+3*(kj+l1)]));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
         for (j = 1; j < nxh; j++) {
            at2 = -ci2*crealf(ffe[j+ll]);
            at1 = at2*cimagf(ffe[j+ll]);
            at2 = at2*at2;
            exyz[3*(j+lj)] = at1*dcu[3*(j+lj)];
            exyz[1+3*(j+lj)] = at1*dcu[1+3*(j+lj)];
            exyz[2+3*(j+lj)] = at1*dcu[2+3*(j+lj)];
            exyz[3*(j+k1+lj)] = zero;
            exyz[1+3*(j+k1+lj)] = zero;
            exyz[2+3*(j+k1+lj)] = zero;
            exyz[3*(j+l1)] = at1*dcu[3*(j+l1)];
            exyz[1+3*(j+l1)] = at1*dcu[1+3*(j+l1)];
            exyz[2+3*(j+l1)] = at1*dcu[2+3*(j+l1)];
            exyz[3*(j+k1+l1)] = zero;
            exyz[1+3*(j+k1+l1)] = zero;
            exyz[2+3*(j+k1+l1)] = zero;
            wp += at2*(dcu[3*(j+lj)]*conjf(dcu[3*(j+lj)])
               + dcu[1+3*(j+lj)]*conjf(dcu[1+3*(j+lj)])
               + dcu[2+3*(j+lj)]*conjf(dcu[2+3*(j+lj)])
               + dcu[3*(j+l1)]*conjf(dcu[3*(j+l1)])
               + dcu[1+3*(j+l1)]*conjf(dcu[1+3*(j+l1)])
               + dcu[2+3*(j+l1)]*conjf(dcu[2+3*(j+l1)]));
         }
/* mode numbers kx = 0, nx/2 */
         at2 = -ci2*crealf(ffe[ll]);
         at1 = at2*cimagf(ffe[ll]);
         at2 = at2*at2;
         exyz[3*lj] = at1*dcu[3*lj];
         exyz[1+3*lj] = at1*dcu[1+3*lj];
         exyz[2+3*lj] = at1*dcu[2+3*lj];
         exyz[3*(k1+lj)] = zero;
         exyz[1+3*(k1+lj)] = zero;
         exyz[2+3*(k1+lj)] = zero;
         exyz[3*l1] = zero;
         exyz[1+3*l1] = zero;
         exyz[2+3*l1] = zero;
         exyz[3*(k1+l1)] = zero;
         exyz[1+3*(k1+l1)] = zero;
         exyz[2+3*(k1+l1)] = zero;
         wp += at2*(dcu[3*lj]*conjf(dcu[3*lj])
            + dcu[1+3*lj]*conjf(dcu[1+3*lj])
            + dcu[2+3*lj]*conjf(dcu[2+3*lj]));
         sum1 += wp;
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   sum2 = 0.0;
#pragma omp parallel for \
private(j,k,k1,kk,kj,at1,at2,wp) \
reduction(+:sum2)
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
         exyz[3*(j+kj)] = at1*dcu[3*(j+kj)];
         exyz[1+3*(j+kj)] = at1*dcu[1+3*(j+kj)];
         exyz[2+3*(j+kj)] = at1*dcu[2+3*(j+kj)];
         exyz[3*(j+k1)] = at1*dcu[3*(j+k1)];
         exyz[1+3*(j+k1)] = at1*dcu[1+3*(j+k1)];
         exyz[2+3*(j+k1)] = at1*dcu[2+3*(j+k1)];
         exyz[3*(j+kj+l1)] = zero;
         exyz[1+3*(j+kj+l1)] = zero;
         exyz[2+3*(j+kj+l1)] = zero;
         exyz[3*(j+k1+l1)] = zero;
         exyz[1+3*(j+k1+l1)] = zero;
         exyz[2+3*(j+k1+l1)] = zero;
         wp += at2*(dcu[3*(j+kj)]*conjf(dcu[3*(j+kj)])
            + dcu[1+3*(j+kj)]*conjf(dcu[1+3*(j+kj)])
            + dcu[2+3*(j+kj)]*conjf(dcu[2+3*(j+kj)])
            + dcu[3*(j+k1)]*conjf(dcu[3*(j+k1)])
            + dcu[1+3*(j+k1)]*conjf(dcu[1+3*(j+k1)])
            + dcu[2+3*(j+k1)]*conjf(dcu[2+3*(j+k1)]));
      }
/* mode numbers kx = 0, nx/2 */
      at2 = -ci2*crealf(ffe[kk]);
      at1 = at2*cimagf(ffe[kk]);
      at2 = at2*at2;
      exyz[3*kj] = at1*dcu[3*kj];
      exyz[1+3*kj] = at1*dcu[1+3*kj];
      exyz[2+3*kj] = at1*dcu[2+3*kj];
      exyz[3*k1] = zero;
      exyz[1+3*k1] = zero;
      exyz[2+3*k1] = zero;
      exyz[3*(kj+l1)] = zero;
      exyz[1+3*(kj+l1)] = zero;
      exyz[2+3*(kj+l1)] = zero;
      exyz[3*(k1+l1)] = zero;
      exyz[1+3*(k1+l1)] = zero;
      exyz[2+3*(k1+l1)] = zero;
      wp += at2*(dcu[3*kj]*conjf(dcu[3*kj])
         + dcu[1+3*kj]*conjf(dcu[1+3*kj])
         + dcu[2+3*kj]*conjf(dcu[2+3*kj]));
      sum2 += wp;
   }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      at2 = -ci2*crealf(ffe[j]);
      at1 = at2*cimagf(ffe[j]);
      at2 = at2*at2;
      exyz[3*j] = at1*dcu[3*j];
      exyz[1+3*j] = at1*dcu[1+3*j];
      exyz[2+3*j] = at1*dcu[2+3*j];
      exyz[3*(j+k1)]= zero;
      exyz[1+3*(j+k1)] = zero;
      exyz[2+3*(j+k1)] = zero;
      exyz[3*(j+l1)] = zero;
      exyz[1+3*(j+l1)] = zero;
      exyz[2+3*(j+l1)] = zero;
      exyz[3*(j+k1+l1)] = zero;
      exyz[1+3*(j+k1+l1)] = zero;
      exyz[2+3*(j+k1+l1)] = zero;
      wp += at2*(dcu[3*j]*conjf(dcu[3*j])
         + dcu[1+3*j]*conjf(dcu[1+3*j])
         + dcu[2+3*j]*conjf(dcu[2+3*j]));
   }
   exyz[0] = zero;
   exyz[1] = zero;
   exyz[2] = zero;
   exyz[3*k1] = zero;
   exyz[1+3*k1] = zero;
   exyz[2+3*k1] = zero;
   exyz[3*l1] = zero;
   exyz[1+3*l1] = zero;
   exyz[2+3*l1] = zero;
   exyz[3*(k1+l1)] = zero;
   exyz[1+3*(k1+l1)] = zero;
   exyz[2+3*(k1+l1)] = zero;
   *wf = (sum1 + sum2 + wp)*((float) nx)*((float) ny)
       *((float) nz)/crealf(ffe[0]);
   return;
/* calculate unsmoothed transverse electric field and sum field energy */
L130: sum1 = 0.0;
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,k1,l1,ll,lj,kk,kj,at1,at2,wp) \
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
            for (j = 1; j < nxh; j++) {
               at2 = -ci2*crealf(ffe[j+kk+ll]);
               at1 = at2*at2;
               exyz[3*(j+kj+lj)] = at2*dcu[3*(j+kj+lj)];
               exyz[1+3*(j+kj+lj)] = at2*dcu[1+3*(j+kj+lj)];
               exyz[2+3*(j+kj+lj)] = at2*dcu[2+3*(j+kj+lj)];
               exyz[3*(j+k1+lj)] = at2*dcu[3*(j+k1+lj)];
               exyz[1+3*(j+k1+lj)] = at2*dcu[1+3*(j+k1+lj)];
               exyz[2+3*(j+k1+lj)] = at2*dcu[2+3*(j+k1+lj)];
               exyz[3*(j+kj+l1)] = at2*dcu[3*(j+kj+l1)];
               exyz[1+3*(j+kj+l1)] = at2*dcu[1+3*(j+kj+l1)];
               exyz[2+3*(j+kj+l1)] = at2*dcu[2+3*(j+kj+l1)];
               exyz[3*(j+k1+l1)] = at2*dcu[3*(j+k1+l1)];
               exyz[1+3*(j+k1+l1)] = at2*dcu[1+3*(j+k1+l1)];
               exyz[2+3*(j+k1+l1)] = at2*dcu[2+3*(j+k1+l1)];
               wp += at1*(dcu[3*(j+kj+lj)]*conjf(dcu[3*(j+kj+lj)])
                  + dcu[1+3*(j+kj+lj)]*conjf(dcu[1+3*(j+kj+lj)])
                  + dcu[2+3*(j+kj+lj)]*conjf(dcu[2+3*(j+kj+lj)])
                  + dcu[3*(j+k1+lj)]*conjf(dcu[3*(j+k1+lj)])
                  + dcu[1+3*(j+k1+lj)]*conjf(dcu[1+3*(j+k1+lj)])
                  + dcu[2+3*(j+k1+lj)]*conjf(dcu[2+3*(j+k1+lj)])
                  + dcu[3*(j+kj+l1)]*conjf(dcu[3*(j+kj+l1)])
                  + dcu[1+3*(j+kj+l1)]*conjf(dcu[1+3*(j+kj+l1)])
                  + dcu[2+3*(j+kj+l1)]*conjf(dcu[2+3*(j+kj+l1)])
                  + dcu[3*(j+k1+l1)]*conjf(dcu[3*(j+k1+l1)])
                  + dcu[1+3*(j+k1+l1)]*conjf(dcu[1+3*(j+k1+l1)])
                  + dcu[2+3*(j+k1+l1)]*conjf(dcu[2+3*(j+k1+l1)]));
            }
         }
/* mode numbers kx = 0, nx/2 */
         for (k = 1; k < nyh; k++) {
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            at2 = -ci2*crealf(ffe[kk+ll]);
            at1 = at2*at2;
            exyz[3*(kj+lj)] = at2*dcu[3*(kj+lj)];
            exyz[1+3*(kj+lj)] = at2*dcu[1+3*(kj+lj)];
            exyz[2+3*(kj+lj)] = at2*dcu[2+3*(kj+lj)];
            exyz[3*(k1+lj)] = zero;
            exyz[1+3*(k1+lj)] = zero;
            exyz[2+3*(k1+lj)] = zero;
            exyz[3*(kj+l1)] = at2*dcu[3*(kj+l1)];
            exyz[1+3*(kj+l1)] = at2*dcu[1+3*(kj+l1)];
            exyz[2+3*(kj+l1)] = at2*dcu[2+3*(kj+l1)];
            exyz[3*(k1+l1)] = zero;
            exyz[1+3*(k1+l1)] = zero;
            exyz[2+3*(k1+l1)] = zero;
            wp += at1*(dcu[3*(kj+lj)]*conjf(dcu[3*(kj+lj)])
               + dcu[1+3*(kj+lj)]*conjf(dcu[1+3*(kj+lj)])
               + dcu[2+3*(kj+lj)]*conjf(dcu[2+3*(kj+lj)])
               + dcu[3*(kj+l1)]*conjf(dcu[3*(kj+l1)])
               + dcu[1+3*(kj+l1)]*conjf(dcu[1+3*(kj+l1)])
               + dcu[2+3*(kj+l1)]*conjf(dcu[2+3*(kj+l1)]));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
         for (j = 1; j < nxh; j++) {
            at2 = -ci2*crealf(ffe[j+ll]);
            at1 = at2*at2;
            exyz[3*(j+lj)] = at2*dcu[3*(j+lj)];
            exyz[1+3*(j+lj)] = at2*dcu[1+3*(j+lj)];
            exyz[2+3*(j+lj)] = at2*dcu[2+3*(j+lj)];
            exyz[3*(j+k1+lj)] = zero;
            exyz[1+3*(j+k1+lj)] = zero;
            exyz[2+3*(j+k1+lj)] = zero;
            exyz[3*(j+l1)] = at2*dcu[3*(j+l1)];
            exyz[1+3*(j+l1)] = at2*dcu[1+3*(j+l1)];
            exyz[2+3*(j+l1)] = at2*dcu[2+3*(j+l1)];
            exyz[3*(j+k1+l1)] = zero;
            exyz[1+3*(j+k1+l1)] = zero;
            exyz[2+3*(j+k1+l1)] = zero;
            wp += at1*(dcu[3*(j+lj)]*conjf(dcu[3*(j+lj)])
               + dcu[1+3*(j+lj)]*conjf(dcu[1+3*(j+lj)])
               + dcu[2+3*(j+lj)]*conjf(dcu[2+3*(j+lj)])
               + dcu[3*(j+l1)]*conjf(dcu[3*(j+l1)])
               + dcu[1+3*(j+l1)]*conjf(dcu[1+3*(j+l1)])
              + dcu[2+3*(j+l1)]*conjf(dcu[2+3*(j+l1)]));
         }
/* mode numbers kx = 0, nx/2 */
         at2 = -ci2*crealf(ffe[ll]);
         at1 = at2*at2;
         exyz[3*lj] = at2*dcu[3*lj];
         exyz[1+3*lj] = at2*dcu[1+3*lj];
         exyz[2+3*lj] = at2*dcu[2+3*lj];
         exyz[3*(k1+lj)] = zero;
         exyz[1+3*(k1+lj)] = zero;
         exyz[2+3*(k1+lj)] = zero;
         exyz[3*l1] = zero;
         exyz[1+3*l1] = zero;
         exyz[2+3*l1] = zero;
         exyz[3*(k1+l1)] = zero;
         exyz[1+3*(k1+l1)] = zero;
         exyz[2+3*(k1+l1)] = zero;
         wp += at1*(dcu[3*lj]*conjf(dcu[3*lj])
            + dcu[1+3*lj]*conjf(dcu[1+3*lj])
            + dcu[2+3*lj]*conjf(dcu[2+3*lj]));
         sum1 += wp;
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   sum2 = 0.0;
#pragma omp parallel for \
private(j,k,k1,kk,kj,at1,at2,wp) \
reduction(+:sum2)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      wp = 0.0;
      for (j = 1; j < nxh; j++) {
         at2 = -ci2*crealf(ffe[j+kk]);
         at1 = at2*at2;
         exyz[3*(j+kj)] = at2*dcu[3*(j+kj)];
         exyz[1+3*(j+kj)] = at2*dcu[1+3*(j+kj)];
         exyz[2+3*(j+kj)] = at2*dcu[2+3*(j+kj)];
         exyz[3*(j+k1)] = at2*dcu[3*(j+k1)];
         exyz[1+3*(j+k1)] = at2*dcu[1+3*(j+k1)];
         exyz[2+3*(j+k1)] = at2*dcu[2+3*(j+k1)];
         exyz[3*(j+kj+l1)] = zero;
         exyz[1+3*(j+kj+l1)] = zero;
         exyz[2+3*(j+kj+l1)] = zero;
         exyz[3*(j+k1+l1)] = zero;
         exyz[1+3*(j+k1+l1)] = zero;
         exyz[2+3*(j+k1+l1)] = zero;
         wp += at1*(dcu[3*(j+kj)]*conjf(dcu[3*(j+kj)])
            + dcu[1+3*(j+kj)]*conjf(dcu[1+3*(j+kj)])
            + dcu[2+3*(j+kj)]*conjf(dcu[2+3*(j+kj)])
            + dcu[3*(j+k1)]*conjf(dcu[3*(j+k1)])
            + dcu[1+3*(j+k1)]*conjf(dcu[1+3*(j+k1)])
            + dcu[2+3*(j+k1)]*conjf(dcu[2+3*(j+k1)]));
      }
/* mode numbers kx = 0, nx/2 */
      at2 = -ci2*crealf(ffe[kk]);
      at1 = at2*at2;
      exyz[3*kj] = at2*dcu[3*kj];
      exyz[1+3*kj] = at2*dcu[1+3*kj];
      exyz[2+3*kj] = at2*dcu[2+3*kj];
      exyz[3*k1]= zero;
      exyz[1+3*k1] = zero;
      exyz[2+3*k1] = zero;
      exyz[3*(kj+l1)] = zero;
      exyz[1+3*(kj+l1)] = zero;
      exyz[2+3*(kj+l1)] = zero;
      exyz[3*(k1+l1)] = zero;
      exyz[1+3*(k1+l1)] = zero;
      exyz[2+3*(k1+l1)] = zero;
      wp += at1*(dcu[3*kj]*conjf(dcu[3*kj])
         + dcu[1+3*kj]*conjf(dcu[1+3*kj])
         + dcu[2+3*kj]*conjf(dcu[2+3*kj]));
      sum2 += wp;
   }
   wp = 0.0;
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      at2 = -ci2*crealf(ffe[j]);
      at1 = at2*at2;
      exyz[3*j] = at2*dcu[3*j];
      exyz[1+3*j] = at2*dcu[1+3*j];
      exyz[2+3*j] = at2*dcu[2+3*j];
      exyz[3*(j+k1)] = zero;
      exyz[1+3*(j+k1)] = zero;
      exyz[2+3*(j+k1)] = zero;
      exyz[3*(j+l1)] = zero;
      exyz[1+3*(j+l1)] = zero;
      exyz[2+3*(j+l1)] = zero;
      exyz[3*(j+k1+l1)] = zero;
      exyz[1+3*(j+k1+l1)] = zero;
      exyz[2+3*(j+k1+l1)] = zero;
      wp += at1*(dcu[3*j]*conjf(dcu[3*j])
         + dcu[1+3*j]*conjf(dcu[1+3*j])
         + dcu[2+3*j]*conjf(dcu[2+3*j]));
   }
   exyz[0] = zero;
   exyz[1] = zero;
   exyz[2] = zero;
   exyz[3*k1] = zero;
   exyz[1+3*k1] = zero;
   exyz[2+3*k1] = zero;
   exyz[3*l1] = zero;
   exyz[1+3*l1] = zero;
   exyz[2+3*l1] = zero;
   exyz[3*(k1+l1)] = zero;
   exyz[1+3*(k1+l1)] = zero;
   exyz[2+3*(k1+l1)] = zero;
   *wf = (sum1+sum2+wp)*((float) nx)*((float) ny)
       *((float) nz)/crealf(ffe[0]);
   return;
}

/*--------------------------------------------------------------------*/
void caddvrfield3(float a[], float b[], float c[], int ndim, int nxe,
                  int nye, int nze) {
/* this subroutine calculates a = b + c for real vector fields
local data                                                  */
   int i, j, k, l, nnxye, ll;
   nnxye = ndim*nxe*nye;
#pragma omp parallel for private(i,j,k,l,ll)
   for (l = 0; l < nze; l++) {
      ll = nnxye*l;
      for (k = 0; k < nye; k++) {
         for (j = 0; j < nxe; j++) {
            for (i = 0; i < ndim; i++) {
               a[i+ndim*(j+nxe*k)+ll] = b[i+ndim*(j+nxe*k)+ll]
                                      + c[i+ndim*(j+nxe*k)+ll];
            }
         }
      }
   }
   return;
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
void cfft3rmxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nzi, int nzp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd) {
/* this subroutine performs the x-y part of a three dimensional real to
   complex fast fourier transform and its inverse, for a subset of z,
   using complex arithmetic, with OpenMP
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
               j1 = j + k1;
               j2 = j + k2;
               t1 = sct[kmr*j];
               for (i = 0; i < ny; i++) {
                  joff = nxhd*i + nn;
                  t2 = t1*f[j2+joff];
                  f[j2+joff] = f[j1+joff] - t2;
                  f[j1+joff] += t2;
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
               j1 = j + k1;
               j2 = j + k2;
               t1 = conjf(sct[kmr*j]);
               for (i = 0; i < ny; i++) {
                  joff = nxhd*i + nn;
                  t2 = t1*f[j2+joff];
                  f[j2+joff] = f[j1+joff] - t2;
                  f[j1+joff] += t2;
               }
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rmxz(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd) {
/* this subroutine performs the z part of a three dimensional real to
   complex fast fourier transform and its inverse, for a subset of y,
   using complex arithmetic, with OpenMP
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
      goto L90;
/* inverse fourier transform */
   nrzb = nxhyz/nz;
   nrz = nxyz/nz;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,l1,i0,i1,ioff,t1,t2)
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
   for (n = 1; n < nzh; n++) {
      ll = nxhyd*n;
      l1 = nxhyd*nz - ll;
      if (nyi==1) {
         t1 = f[l1];
         f[l1] = 0.5*(cimagf(f[ll] + t1)
                    + crealf(f[ll] - t1)*_Complex_I);
         f[ll] = 0.5*(crealf(f[ll] + t1)
                    + cimagf(f[ll] - t1)*_Complex_I);
      }
      if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
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
L90: nrzb = nxhyz/nz;
   nrz = nxyz/nz;
/* scramble modes kx = 0, nx/2 */
   for (n = 1; n < nzh; n++) {
      ll = nxhyd*n;
      l1 = nxhyd*nz - ll;
      if (nyi==1) {
         t1 = cimagf(f[l1]) + crealf(f[l1])*_Complex_I;
         f[l1] = conjf(f[ll] - t1);
         f[ll] += t1;
      }
      if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
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
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,l1,i0,i1,ioff,t1,t2)
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
void cfft3rm3xy(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nzi, int nzp, int nxhd, int nyd, int nzd, 
                int nxhyzd, int nxyzhd) {
/* this subroutine performs the x-y part of 3 three dimensional complex
   to real fast fourier transforms and their inverses, for a subset of z,
   using complex arithmetic, with OpenMP
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
   int indx1, ndx1yz, nx, nxh, nxhh, ny, nyh;
   int nz, nxyz, nxhyz, nzt, nrx, nry, nrxb, nryb, nxhd3, nxhyd;
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
   nxhd3 = 3*nxhd;
   nxhyd = nxhd3*nyd;
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
         joff = nxhd3*i + nn;
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
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            for (i = 0; i < ny; i++) {
               joff = nxhd3*i + nn;
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
               for (i = 0; i < ny; i++) {
                  joff = nxhd3*i + nn;
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
      kmr = nxyz/nx;
      ani = 0.5/(((float) nx)*((float) ny)*((float) nz));
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < ny; k++) {
            joff = nxhd3*k + nn;
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
      for (k = 0; k < ny; k++) {
         joff = nxhd3*k + nn;
         for (jj = 0; jj < 3; jj++) {
            f[jj+3*nxhh+joff] = ani*conjf(f[jj+3*nxhh+joff]);
            f[jj+joff] = ani*((crealf(f[jj+joff]) 
                          + cimagf(f[jj+joff]))
                          + (crealf(f[jj+joff])
                          - cimagf(f[jj+joff]))*_Complex_I);
         }
      }
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         joff = nxhd3*k + nn;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nxhd3*k1 + nn;
            for (i = 0; i < nxh; i++) {
               t1 = f[3*i+k1];
               t2 = f[1+3*i+k1];
               t3 = f[2+3*i+k1];
               f[3*i+k1] = f[3*i+joff];
               f[1+3*i+k1] = f[1+3*i+joff];
               f[2+3*i+k1] = f[2+3*i+joff];
               f[3*i+joff] = t1;
               f[1+3*i+joff] = t2;
               f[2+3*i+joff] = t3;
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
               j1 = nxhd3*(j + k1) + nn;
               j2 = nxhd3*(j + k2) + nn;
               t1 = sct[kmr*j];
               for (i = 0; i < nxh; i++) {
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
         joff = nxhd3*k;
         k1 = nxhd3*ny - joff + nn;
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
         joff = nxhd3*k;
         k1 = nxhd3*ny - joff + nn;
         joff += nn;
         for (jj = 0; jj < 3; jj++) {
            t1 = cimagf(f[jj+k1]) + crealf(f[jj+k1])*_Complex_I;
            f[jj+k1] = conjf(f[jj+joff] - t1);
            f[jj+joff] += t1;
         }
      }
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         joff = nxhd3*k + nn;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nxhd3*k1 + nn;
            for (i = 0; i < nxh; i++) {
               t1 = f[3*i+k1];
               t2 = f[1+3*i+k1];
               t3 = f[2+3*i+k1];
               f[3*i+k1] = f[3*i+joff];
               f[1+3*i+k1] = f[1+3*i+joff];
               f[2+3*i+k1] = f[2+3*i+joff];
               f[3*i+joff] = t1;
               f[1+3*i+joff] = t2;
               f[2+3*i+joff] = t3;
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
               j1 = nxhd3*(j + k1) + nn;
               j2 = nxhd3*(j + k2) + nn;
               t1 = conjf(sct[kmr*j]);
               for (i = 0; i < nxh; i++) {
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
/* scramble coefficients */
      kmr = nxyz/nx;
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < ny; k++) {
            joff = nxhd3*k + nn;
            for (jj = 0; jj < 3; jj++) {
               t2 = conjf(f[jj+3*(nxh-j)+joff]);
               t1 = f[jj+3*j+joff] + t2;
               t2 = (f[jj+3*j+joff] - t2)*t3;
               f[jj+3*j+joff] = t1 + t2;
               f[jj+3*(nxh-j)+joff] = conjf(t1 - t2);
            }
         }
      }
      for (k = 0; k < ny; k++) {
         joff = nxhd3*k + nn;
         for (jj = 0; jj < 3; jj++) {
            f[jj+3*nxhh+joff] = 2.0*conjf(f[jj+3*nxhh+joff]);
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
               joff = nxhd3*i + nn;
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
               j1 = j + k1;
               j2 = j + k2;
               t1 = conjf(sct[kmr*j]);
               for (i = 0; i < ny; i++) {
                  joff = nxhd3*i + nn;
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
      for (i = 0; i < ny; i++) {
         joff = nxhd3*i + nn;
         for (j = 0; j < nxh; j++) {
            at1 = crealf(f[2+3*j+joff]);
            f[2+3*j+joff] = cimagf(f[1+3*j+joff])
                            + cimagf(f[2+3*j+joff])*_Complex_I;
            at2 = crealf(f[1+3*j+joff]);
            f[1+3*j+joff] = at1 + cimagf(f[3*j+joff])*_Complex_I;
            f[3*j+joff] = crealf(f[3*j+joff]) + at2*_Complex_I;
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rm3z(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd, 
               int nxyzhd) {
/* this subroutine performs the z part of 3 three dimensional complex to
   real fast fourier transforms and their inverses, for a subset of y,
   using complex arithmetic, with OpenMP
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
   int indx1, ndx1yz, nx, nxh, ny, nyh;
   int nz, nzh, nxyz, nxhyz, nyt, nrz, nrzb, nxhd3, nxhyd, ioff;
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
   nxhd3 = 3*nxhd;
   nxhyd = nxhd3*nyd;
   if (isign > 0)
      goto L110;
/* inverse fourier transform */
   nrzb = nxhyz/nz;
   nrz = nxyz/nz;
/* bit-reverse array elements in z */
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,l1,i0,i1,ioff,t1,t2,t3,t4)
   for (n = nyi-1; n < nyt; n++) {
      ioff = nxhd3*n;
      for (l = 0; l < nz; l++) {
         ll = nxhyd*l;
         l1 = (mixup[l] - 1)/nrzb;
         if (l < l1) {
            l1 = nxhyd*l1;
            i0 = ioff + ll;
            i1 = ioff + l1;
            for (i = 0; i < nxh; i++) {
               t1 = f[3*i+i1];
               t2 = f[1+3*i+i1];
               t3 = f[2+3*i+i1];
               f[3*i+i1] = f[3*i+i0];
               f[1+3*i+i1] = f[1+3*i+i0];
               f[2+3*i+i1] = f[2+3*i+i0];
               f[3*i+i0] = t1;
               f[1+3*i+i0] = t2;
               f[2+3*i+i0] = t3;
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
                  t2 = t1*f[3*i+i1];
                  t3 = t1*f[1+3*i+i1];
                  t4 = t1*f[2+3*i+i1];
                  f[3*i+i1] = f[3*i+i0] - t2;
                  f[1+3*i+i1] = f[1+3*i+i0] - t3;
                  f[2+3*i+i1] = f[2+3*i+i0] - t4;
                  f[3*i+i0] += t2;
                  f[1+3*i+i0] += t3;
                  f[2+3*i+i0] += t4;
               }

            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   for (n = 1; n < nzh; n++) {
      ll = nxhyd*n;
      l1 = nxhyd*nz - ll;
      if (nyi==1) {
         for (jj = 0; jj < 3; jj++) {
            t1 = f[jj+l1];
            f[jj+l1] = 0.5*(cimagf(f[jj+ll] + t1)
                          + crealf(f[jj+ll] - t1)*_Complex_I);
            f[jj+ll] = 0.5*(crealf(f[jj+ll] + t1)
                          + cimagf(f[jj+ll] - t1)*_Complex_I);
         }
      }
      if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
         for (jj = 0; jj < 3; jj++) {
            i1 = nxhd3*nyh;
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
   for (n = 1; n < nzh; n++) {
      ll = nxhyd*n;
      l1 = nxhyd*nz - ll;
      if (nyi==1) {
         for (jj = 0; jj < 3; jj++) {
            t1 = cimagf(f[jj+l1]) + crealf(f[jj+l1])*_Complex_I;
            f[jj+l1] = conjf(f[jj+ll] - t1);
            f[jj+ll] += t1;
         }
      }
      if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
         for (jj = 0; jj < 3; jj++) {
            i1 = nxhd3*nyh;
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
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,l1,i0,i1,ioff,t1,t2,t3,t4)
   for (n = nyi-1; n < nyt; n++) {
      ioff = nxhd3*n;
      for (l = 0; l < nz; l++) {
         ll = nxhyd*l;
         l1 = (mixup[l] - 1)/nrzb;
         if (l < l1) {
            l1 = nxhyd*l1;
            i0 = ioff + ll;
            i1 = ioff+ l1;
            for (i = 0; i < nxh; i++) {
               t1 = f[3*i+i1];
               t2 = f[1+3*i+i1];
               t3 = f[2+3*i+i1];
               f[3*i+i1] = f[3*i+i0];
               f[1+3*i+i1] = f[1+3*i+i0];
               f[2+3*i+i1] = f[2+3*i+i0];
               f[3*i+i0] = t1;
               f[1+3*i+i0] = t2;
               f[2+3*i+i0] = t3;
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
                  t2 = t1*f[3*i+i1];
                  t3 = t1*f[1+3*i+i1];
                  t4 = t1*f[2+3*i+i1];
                  f[3*i+i1] = f[3*i+i0] - t2;
                  f[1+3*i+i1] = f[1+3*i+i0] - t3;
                  f[2+3*i+i1] = f[2+3*i+i0] - t4;
                  f[3*i+i0] += t2;
                  f[1+3*i+i0] += t3;
                  f[2+3*i+i0] += t4;
               }
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cmswap3cn(float f[], float s[], int isign, int nxh, int ny,
               int nzi, int nzt, int nxhd, int nyd, int nzd, int ndim) {
/* this subroutine swaps components for multiple ffts
   f = input array
   s = scratch array
   isign = (-1,1) = swap (real-to-complex,complex-to-real)
   nxh = complex dimension in x direction
   ny = complex dimension in y direction
   nzi/nzt = initial/final z index used
   nxhd = half of the second dimension of f
   nyd,nzd = third and fourth dimension of f
   ndim = leading dimension of array f
local data                                                            */
   int i, j, k, l, ll, ioff, nnxhd2, nnxhyd, nl, nk;
/* swap complex components */
      nnxhd2 = 2*ndim*nxhd;
      nnxhyd = nnxhd2*nyd;
/* real to complex */
   if (isign < 0) {
#pragma omp for private(i,j,k,l,ll,nl,nk,ioff)
      for (l = nzi-1; l < nzt; l++) {
         ll = nnxhd2*l;
         nl = nnxhyd*l;
         for (k = 0; k < ny; k++) {
            nk = nnxhd2*k + nl;
            for (j = 0; j < nxh; j++) {
               ioff = 2*ndim*j + ll;
               for (i = 0; i < ndim; i++) {
                  s[2*i+ioff] = f[i+ndim*(2*j)+nk];
                  s[2*i+ioff+1] = f[i+ndim*(2*j+1)+nk];
               }
            }
            for (j = 0; j < nxh; j++) {
               ioff = 2*ndim*j + ll;
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
   }
/* complex to real */
   else if (isign > 0) {
#pragma omp for private(i,j,k,l,ll,nl,nk,ioff)
      for (l = nzi-1; l < nzt; l++) {
         ll = nnxhd2*l;
         nl = nnxhyd*l;
         for (k = 0; k < ny; k++) {
            nk = nnxhd2*k + nl;
            for (j = 0; j < nxh; j++) {
               ioff = 2*ndim*j + ll;
               for (i = 0; i < ndim; i++) {
                  s[i+ioff] = f[i+ndim*(2*j)+nk];
               }
               ioff += ndim;
               for (i = 0; i < ndim; i++) {
                  s[i+ioff] = f[i+ndim*(2*j+1)+nk];
               }
            }
            for (j = 0; j < nxh; j++) {
               ioff = 2*ndim*j + ll;
               for (i = 0; i < ndim; i++) {
                  f[i+ndim*(2*j)+nk] = s[2*i+ioff];
                  f[i+ndim*(2*j+1)+nk] = s[2*i+ioff+1];
               }
            }
         }
      }
   }
   return;
}


/*--------------------------------------------------------------------*/
void cfft3rmnxy(float complex f[], float complex ss[], int isign,
                int mixup[], float complex sct[], int indx, int indy,
                int indz, int nzi, int nzp, int nxhd, int nyd, int nzd,
                int ndim, int nxhyzd, int nxyzhd) {
/* this subroutine performs the x-y part of N three dimensional complex
   to real fast fourier transforms and their inverses, for a subset of z,
   using complex arithmetic, where N = ndim, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)
   for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)
   where M = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, N inverse fourier transforms in x and y are performed
   f[i][m][n][0:N-1] = (1/nx*ny*nz)*sum(f[i][k][j][0:N-1]*
                     exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, N forward fourier transforms in x and y are performed
   f[l][k][j][0:N-1] = sum([l][m][n][0:N-1]*exp(sqrt(-1)*2pi*n*j/nx)*
                       exp(sqrt(-1)*2pi*m*k/ny))
   ss = scratch array
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nzi = initial z index used
   nzp = number of z indices used
   nxhd = second dimension of f
   nyd,nzd = third and fourth dimensions of f
   ndim = leading dimension of array f
   nxhyzd = maximum of (nx/2,ny,nz)
   nxyzhd = maximum of (nx,ny,nz)/2
   fourier coefficients are stored as follows:
   f[l][k][j][0:N-1] = real, imaginary part of mode j,k,l
   where 0 <= j < nx/2, 0 <= k < ny, 0 <= l < nz, except for
   f[l][k][0][0:N-1] = real, imaginary part of mode nx/2,k,l,
   where ny/2+1 <= k < ny and 0 <= l < nz, and
   f[l][0][0][0:N-1] = real, imaginary part of mode nx/2,0,l,
   f[l][ny/2][0][0:N-1] = real, imaginary part mode nx/2,ny/2,l,
   where nz/2+1 <= l < nz, and
   imag(f[0][0][0][0:N-1]) = real part of mode nx/2,0,0
   imag(f[0][ny/2][0][0:N-1]) = real part of mode nx/2,ny/2,0
   imag(f[nz/2][0][0][0:N-1]) = real part of mode nx/2,0,nz/2
   imag(f[nz/2][ny/2][0][0:N-1]) = real part of mode nx/2,ny/2,nz/2
   using jpl storage convention, as described in:
   E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
   Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
   Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
   December 1993.
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, nxhh, ny, nyh, nnxhd, nnxhyd;
   int nz, nxyz, nxhyz, nzt, nrx, nry, ns, ns2, km, kmr, joff;
   int i, j, k, l, n, nn, k1, k2, j1, j2, jj, nrxb, nryb;
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
   nnxhd = ndim*nxhd;
   nnxhyd = nnxhd*nyd;
   if (isign > 0)
      goto L250;
/* inverse fourier transform */
   nrxb = nxhyz/nxh;
   nrx = nxyz/nxh;
   nryb = nxhyz/ny;
   nry = nxyz/ny;
/* swap complex components */
   cmswap3cn((float *)f,(float *)ss,isign,nxh,ny,nzi,nzt,nxhd,nyd,nzd,
             ndim);
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,jj,j1,j2,nn,joff,ani,t1,t2,t3)
   for (n = nzi-1; n < nzt; n++) {
      nn = nnxhyd*n;
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            for (i = 0; i < ny; i++) {
               joff = nnxhd*i + nn;
               for (jj = 0; jj < ndim; jj++) {
                  t1 = f[jj+ndim*j1+joff];
                  f[jj+ndim*j1+joff] = f[jj+ndim*j+joff];
                  f[jj+ndim*j+joff] = t1;
               }
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
               for (i = 0; i < ny; i++) {
                  joff = nnxhd*i + nn;
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
      kmr = nxyz/nx;
      ani = 0.5/(((float) nx)*((float) ny)*((float) nz));
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < ny; k++) {
            joff = nnxhd*k + nn;
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
      for (k = 0; k < ny; k++) {
         joff = nnxhd*k + nn;
         for (jj = 0; jj < ndim; jj++) {
            f[jj+ndim*nxhh+joff] = ani*conjf(f[jj+ndim*nxhh+joff]);
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
               for (jj = 0; jj < ndim; jj++) {
                  t1 = f[jj+ndim*i+k1];
                  f[jj+ndim*i+k1] = f[jj+ndim*i+joff];
                  f[jj+ndim*i+joff] = t1;
               }
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
         joff = nnxhd*k;
         k1 = nnxhd*ny - joff + nn;
         joff += nn;
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
L250: nryb = nxhyz/ny;
   nry = nxyz/ny;
   nrxb = nxhyz/nxh;
   nrx = nxyz/nxh;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,jj,j1,j2,nn,joff,ani,t1,t2,t3)
   for (n = nzi-1; n < nzt; n++) {
      nn = nnxhyd*n;
/* scramble modes kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         joff = nnxhd*k;
         k1 = nnxhd*ny - joff + nn;
         joff += nn;
         for (jj = 0; jj < ndim; jj++) {
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
               for (jj = 0; jj < ndim; jj++) {
                  t1 = f[jj+ndim*i+k1];
                  f[jj+ndim*i+k1] = f[jj+ndim*i+joff];
                  f[jj+ndim*i+joff] = t1;
               }
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
/* scramble coefficients */
      kmr = nxyz/nx;
      for (j = 1; j < nxhh; j++) {
         t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < ny; k++) {
            joff = nnxhd*k + nn;
            for (jj = 0; jj < ndim; jj++) {
               t2 = conjf(f[jj+ndim*(nxh-j)+joff]);
               t1 = f[jj+ndim*j+joff] + t2;
               t2 = (f[jj+ndim*j+joff] - t2)*t3;
               f[jj+ndim*j+joff] = t1 + t2;
               f[jj+ndim*(nxh-j)+joff] = conjf(t1 - t2);
            }
         }
      }
      for (k = 0; k < ny; k++) {
         joff = nnxhd*k + nn;
         for (jj = 0; jj < ndim; jj++) {
            f[jj+ndim*nxhh+joff] = 2.0*conjf(f[jj+ndim*nxhh+joff]);
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
               for (jj = 0; jj < ndim; jj++) {
                  t1 = f[jj+ndim*j1+joff];
                  f[jj+ndim*j1+joff] = f[jj+ndim*j+joff];
                  f[jj+ndim*j+joff] = t1;
               }
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
               j1 = j + k1;
               j2 = j + k2;
               t1 = conjf(sct[kmr*j]);
               for (i = 0; i < ny; i++) {
                  joff = nnxhd*i + nn;
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
   }
/* swap complex components */
   cmswap3cn((float *)f,(float *)ss,isign,nxh,ny,nzi,nzt,nxhd,nyd,nzd,
             ndim);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rmnz(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nyi, int nyp, int nxhd, int nyd, int nzd, int ndim,
               int nxhyzd, int nxyzhd) {
/* this subroutine performs the z part of N three dimensional complex to
   real fast fourier transforms and their inverses, for a subset of y,
   using complex arithmetic, where N = ndim, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)
   for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)
   where M = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, three inverse fourier transforms is performed
   f[l][k][j][0:N-1] = sum(f[i][k][j][0:N-1]*exp(-sqrt(-1)*2pi*l*i/nz))
   if isign = 1, three forward fourier transforms are performed
   f[i][m][n][0:N-1] = sum(f[l][m][n][0:N-1]*exp(sqrt(-1)*2pi*l*i/nz))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = second dimension of f
   nyd,nzd = third and fourth dimensions of f
   ndim = leading dimension of array f
   nxhyzd = maximum of (nx/2,ny,nz)
   nxyzhd = maximum of (nx,ny,nz)/2
   fourier coefficients are stored as follows:
   f[l][k][j][0:N-1] = real, imaginary part of mode j,k,l
   where 0 <= j < nx/2, 0 <= k < ny, 0 <= l < nz, except for
   f[l][k][0][0:N-1], = real, imaginary part of mode nx/2,k,l,
   where ny/2+1 <= k < ny and 0 <= l < nz, and
   f[l][0][0][0:N-1] = real, imaginary part of mode nx/2,0,l,
   f[l][ny/2][0][0:N-1] = real, imaginary part mode nx/2,ny/2,l,
   where nz/2+1 <= l < nz, and
   imag(f[0][0][0][0:N-1]) = real part of mode nx/2,0,0
   imag(f[0][ny/2][0][0:N-1]) = real part of mode nx/2,ny/2,0
   imag(f[nz/2][0][0][0:N-1]) = real part of mode nx/2,0,nz/2
   imag(f[nz/2][ny/2][0][0:N-1]) = real part of mode nx/2,ny/2,nz/2
   using jpl storage convention, as described in:
   E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
   Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
   Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
   December 1993.
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, ny, nyh, nnxhd, nnxhyd;
   int nz, nzh, nxyz, nxhyz, nyt, nrz, ns, ns2, km, kmr;
   int i, j, k, l, n, nn, ll, k1, k2, j1, j2, l1, jj, i0, i1, nrzb;
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
   nnxhd = ndim*nxhd;
   nnxhyd = nnxhd*nyd;
   if (isign > 0)
      goto L130;
/* inverse fourier transform */
   nrzb = nxhyz/nz;
   nrz = nxyz/nz;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,jj,j1,j2,nn,ll,l1,i0,i1,t1,t2)
   for (n = nyi-1; n < nyt; n++) {
      nn = nnxhd*n;
/* bit-reverse array elements in z */
      for (l = 0; l < nz; l++) {
         ll = nnxhyd*l;
         l1 = (mixup[l] - 1)/nrzb;
         if (l < l1) {
            l1 = nnxhyd*l1;
            i1 = nn;
            i0 = i1 + ll;
            i1 += l1;
            for (i = 0; i < nxh; i++) {
               for (jj = 0; jj < ndim; jj++) {
                  t1 = f[jj+ndim*i+i1];
                  f[jj+ndim*i+i1] = f[jj+ndim*i+i0];
                  f[jj+ndim*i+i0] = t1;
               }
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
               j1 = nnxhyd*(j + k1);
               j2 = nnxhyd*(j + k2);
               t1 = sct[kmr*j];
               i1 = nn;
               i0 = i1 + j1;
               i1 += j2;
               for (i = 0; i < nxh; i++) {
                  for (jj = 0; jj < ndim; jj++) {
                     t2 = t1*f[jj+ndim*i+i1];
                     f[jj+ndim*i+i1] = f[jj+ndim*i+i0] - t2;
                     f[jj+ndim*i+i0] += t2;
                  }
               }

            }
         }
      ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   for (n = 1; n < nzh; n++) {
      ll = nnxhyd*n;
      l1 = nnxhyd*nz - ll;
      if (nyi==1) {
         for (jj = 0; jj < ndim; jj++) {
            t1 = f[jj+l1];
            f[jj+l1] = 0.5*(cimagf(f[jj+ll] + t1)
                          + crealf(f[jj+ll] - t1)*_Complex_I);
            f[jj+ll] = 0.5*(crealf(f[jj+ll] + t1)
                          + cimagf(f[jj+ll] - t1)*_Complex_I);
         }
      }
      if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
         for (jj = 0; jj < ndim; jj++) {
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
L130: nrzb = nxhyz/nz;
   nrz = nxyz/nz;
/* scramble modes kx = 0, nx/2 */
   for (n = 1; n < nzh; n++) {
      ll = nnxhyd*n;
      l1 = nnxhyd*nz - ll;
      if (nyi==1) {
         for (jj = 0; jj < ndim; jj++) {
            t1 = cimagf(f[jj+l1]) + crealf(f[jj+l1])*_Complex_I;
            f[jj+l1] = conjf(f[jj+ll] - t1);
            f[jj+ll] += t1;
         }
      }
      if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
         for (jj = 0; jj < ndim; jj++) {
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
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,jj,j1,j2,nn,ll,l1,i0,i1,t1,t2)
   for (n = nyi-1; n < nyt; n++) {
      nn = nnxhd*n;
      for (l = 0; l < nz; l++) {
         ll = nnxhyd*l;
         l1 = (mixup[l] - 1)/nrzb;
         if (l < l1) {
            l1 = nnxhyd*l1;
            i1 = nn;
            i0 = i1 + ll;
            i1 += l1;
            for (i = 0; i < nxh; i++) {
               for (jj = 0; jj < ndim; jj++) {
                  t1 = f[jj+ndim*i+i1];
                  f[jj+ndim*i+i1] = f[jj+ndim*i+i0];
                  f[jj+ndim*i+i0] = t1;
               }
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
               j1 = nnxhyd*(j + k1);
               j2 = nnxhyd*(j + k2);
               t1 = conjf(sct[kmr*j]);
               i1 = nn;
               i0 = i1 + j1;
               i1 += j2;
               for (i = 0; i < nxh; i++) {
                  for (jj = 0; jj < ndim; jj++) {
                     t2 = t1*f[jj+ndim*i+i1];
                     f[jj+ndim*i+i1] = f[jj+ndim*i+i0] - t2;
                     f[jj+ndim*i+i0] += t2;
                  }
               }
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rmx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
/* wrapper function for real to complex fft, with packed data */
/* local data */
   int ny, nz;
   static int nyi = 1, nzi = 1;
/* calculate range of indices */
   ny = 1L<<indy;
   nz = 1L<<indz;
/* inverse fourier transform */
   if (isign < 0) {
/* perform xy fft */
      cfft3rmxy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                nxhyzd,nxyzhd);
/* perform z fft */
      cfft3rmxz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                nxhyzd,nxyzhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform z fft */
      cfft3rmxz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                nxhyzd,nxyzhd);
/* perform xy fft */
      cfft3rmxy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                nxhyzd,nxyzhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rm3(float complex f[], int isign, int mixup[],
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
      cfft3rm3xy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                 nxhyzd,nxyzhd);
/* perform z fft */
      cfft3rm3z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                nxhyzd,nxyzhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform z fft */
      cfft3rm3z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                nxhyzd,nxyzhd);
/* perform xy fft */
      cfft3rm3xy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                 nxhyzd,nxyzhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rmn(float complex f[], float complex ss[], int isign,
               int mixup[], float complex sct[], int indx, int indy,
               int indz, int nxhd, int nyd, int nzd, int ndim, 
               int nxhyzd, int nxyzhd) {
/* wrapper function for multiple 3d real to complex ffts, packed data */
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
      cfft3rmnxy(f,ss,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
                 nzd,ndim,nxhyzd,nxyzhd);
/* perform z fft */
      cfft3rmnz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                ndim,nxhyzd,nxyzhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform z fft */
      cfft3rmnz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                ndim,nxhyzd,nxyzhd);
/* perform xy fft */
      cfft3rmnxy(f,ss,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
                 nzd,ndim,nxhyzd,nxyzhd);
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
void cppmovin3l_(float *part, float *ppart, int *kpic, int *nppmx,
                 int *idimp, int *nop, int *mx, int *my, int *mz,
                 int *mx1, int *my1, int *mxyz1, int *irc) {
   cppmovin3l(part,ppart,kpic,*nppmx,*idimp,*nop,*mx,*my,*mz,*mx1,*my1,
              *mxyz1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcheck3l_(float *ppart, int *kpic, int *idimp, int *nppmx,
                 int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                 int *mx1, int *my1, int *mz1, int *irc) {
   cppcheck3l(ppart,kpic,*idimp,*nppmx,*nx,*ny,*nz,*mx,*my,*mz,*mx1,
              *my1,*mz1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppush3l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                 float *qbm, float *dt, float *dtc, float *ek,
                 int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                 int *mx, int *my, int *mz, int *nxv, int *nyv,
                 int *nzv, int *mx1, int *my1, int *mxyz1, int *ipbc) {
   cgbppush3l(ppart,fxyz,bxyz,kpic,*qbm,*dt,*dtc,ek,*idimp,*nppmx,*nx,
              *ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,* mx1,*my1,*mxyz1,
              *ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppushf3l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                  int *ncl, int *ihole, float *qbm, float *dt,
                  float *dtc, float *ek, int *idimp, int *nppmx,
                  int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                  int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                  int *mxyz1, int *ntmax, int *irc) {
   cgbppushf3l(ppart,fxyz,bxyz,kpic,ncl,ihole,*qbm,*dt,*dtc,ek,*idimp,
               *nppmx,*nx,*ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,
               *mxyz1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppost3l_(float *ppart, float *q, int *kpic, float *qm,
                int *nppmx, int *idimp, int *mx, int *my, int *mz,
                int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                int *mxyz1) {
   cgppost3l(ppart,q,kpic,*qm,*nppmx,*idimp,*mx,*my,*mz,*nxv,*nyv,*nzv,
             *mx1,*my1,*mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cgjppost3l_(float *ppart, float *cu, int *kpic, float *qm,
                 float *dt, int *nppmx, int *idimp, int *nx, int *ny,
                 int *nz, int *mx, int *my, int *mz, int *nxv, int *nyv,
                 int *nzv, int *mx1, int *my1, int *mxyz1, int *ipbc) {
   cgjppost3l(ppart,cu,kpic,*qm,*dt,*nppmx,*idimp,*nx,*ny,*nz,*mx,*my,
              *mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgmjppost3l_(float *ppart, float *amu, int *kpic, float *qm,
                  int *nppmx, int *idimp, int *mx, int *my, int *mz,
                  int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                  int *mxyz1) {
   cgmjppost3l(ppart,amu,kpic,*qm,*nppmx,*idimp,*mx,*my,*mz,*nxv,*nyv,
               *nzv,*mx1,*my1,*mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cgdjppost3l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                  float *dcu, float *amu, float *qm, float *qbm,
                  float *dt, int *idimp, int *nppmx, int *nx, int *ny,
                  int *nz, int *mx, int *my, int *mz, int *nxv,
                  int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1) {
   cgdjppost3l(ppart,fxyz,bxyz,kpic,dcu,amu,*qm,*qbm,*dt,*idimp,*nppmx,
               *nx,*ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cgdcjppost3l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                  float *cu, float *dcu, float *amu, float *qm,
                  float *qbm, float *dt, int *idimp, int *nppmx,
                  int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                  int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                  int *mxyz1) {
   cgdcjppost3l(ppart,fxyz,bxyz,kpic,cu,dcu,amu,*qm,*qbm,*dt,*idimp,
                *nppmx,*nx,*ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,
                *mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cpporder3l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                 int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                 int *nz, int *mx, int *my, int *mz, int *mx1,
                 int *my1, int *mz1, int *npbmx, int *ntmax, int *irc) {
   cpporder3l(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*nx,*ny,*nz,*mx,
              *my,*mz,*mx1,*my1,*mz1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf3l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                  int *ihole, int *idimp, int *nppmx, int *mx1, 
                  int *my1, int *mz1, int *npbmx, int *ntmax,
                  int *irc) {
   cpporderf3l(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*mx1,*my1,*mz1,
               *npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void ccguard3l_(float *fxyz, int *nx, int *ny, int *nz, int *nxe,
                int *nye, int *nze) {
   ccguard3l(fxyz,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void cacguard3l_(float *cu, int *nx, int *ny, int *nz, int *nxe, 
                int *nye, int *nze) {
   cacguard3l(cu,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void caguard3l_(float *q, int *nx, int *ny, int *nz, int *nxe, int *nye,
                int *nze) {
   caguard3l(q,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void camcguard3l_(float *amu, int *nx, int *ny, int *nz, int *nxe,
                  int *nye, int *nze, int *ndim) {
   camcguard3l(amu,*nx,*ny,*nz,*nxe,*nye,*nze,*ndim);
   return;
}

/*--------------------------------------------------------------------*/
void cascfguard3l_(float *dcu, float *cus, float *q2m0, int *nx,
                   int *ny, int *nz, int *nxe, int *nye, int *nze) {
   cascfguard3l(dcu,cus,*q2m0,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void cfwpminmx3_(float *qe, float *qbme, float *wpmax, float *wpmin,
                 int *nx, int *ny, int *nz, int *nxe, int *nye,
                 int *nze) {
   cfwpminmx3(qe,*qbme,wpmax,wpmin,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void cmpois33_(float complex *q, float complex *fxyz, int *isign,
               float complex *ffc, float *ax, float *ay, float *az,
               float *affp, float *we, int *nx, int *ny, int *nz,
               int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
               int *nzhd) {
   cmpois33(q,fxyz,*isign,ffc,*ax,*ay,*az,*affp,we,*nx,*ny,*nz,*nxvh,
            *nyv,*nzv,*nxhd,*nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmcuperp3_(float complex *cu, int *nx, int *ny, int *nz, int *nxvh,
                int *nyv, int *nzv) {
   cmcuperp3(cu,*nx,*ny,*nz,*nxvh,*nyv,*nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cmbbpois33_(float complex *cu, float complex *bxyz,
                 float complex *ffc, float *ci, float *wm, int *nx,
                 int *ny, int *nz, int *nxvh, int *nyv, int *nzv,
                 int *nxhd, int *nyhd, int *nzhd) {
   cmbbpois33(cu,bxyz,ffc,*ci,wm,*nx,*ny,*nz,*nxvh,*nyv,*nzv,*nxhd,
              *nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cbaddext3_(float *bxyz, float *omx, float *omy, float *omz,
                int *nx, int *ny, int *nz, int *nxe, int *nye,
                int *nze) {
   cbaddext3(bxyz,*omx,*omy,*omz,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void cmdcuperp3_(float complex *dcu, float complex *amu, int *nx,
                 int *ny, int *nz, int *nxvh, int *nyv, int *nzv) {
   cmdcuperp3(dcu,amu,*nx,*ny,*nz,*nxvh,*nyv,*nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cmadcuperp3_(float complex *dcu, float complex *amu, int *nx,
                  int *ny, int *nz, int *nxvh, int *nyv, int *nzv) {
   cmadcuperp3(dcu,amu,*nx,*ny,*nz,*nxvh,*nyv,*nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cmepois33_(float complex *dcu, float complex *exyz, int *isign,
                float complex *ffe, float *ax, float *ay, float *az,
                float *affp, float *wp0, float *ci, float *wf, int *nx,
                int *ny, int *nz, int *nxvh, int *nyv, int *nzv, 
                int *nxhd, int *nyhd, int *nzhd) {
   cmepois33(dcu,exyz,*isign,ffe,*ax,*ay,*az,*affp,*wp0,*ci,wf,*nx,*ny,
             *nz,*nxvh,*nyv,*nzv,*nxhd,*nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void caddvrfield3_(float *a, float *b, float *c, int *ndim, int *nxe,
                   int *nye, int *nze) {
   caddvrfield3(a,b,c,*ndim,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                  int *indz, int *nxhyzd, int *nxyzhd) {
   cwfft3rinit(mixup,sct,*indx,*indy,*indz,*nxhyzd,*nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rmx_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                int *nxyzhd) {
   cwfft3rmx(f,*isign,mixup,sct,*indx,*indy,*indz,*nxhd,*nyd,*nzd,
             *nxhyzd,*nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rm3_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                int *nxyzhd) {
   cwfft3rm3(f,*isign,mixup,sct,*indx,*indy,*indz,*nxhd,*nyd,*nzd,
             *nxhyzd,*nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rmn_(float complex *f, float complex *ss, int *isign,
               int *mixup, float complex *sct, int *indx, int *indy,
               int *indz, int *nxhd, int *nyd, int *nzd, int *ndim, 
               int *nxhyzd, int *nxyzhd) {
   cwfft3rmn(f,ss,*isign,mixup,sct,*indx,*indy,*indz,*nxhd,*nyd,*nzd,
             *ndim,*nxhyzd,*nxyzhd);
   return;
}
