/* C Library for Skeleton 2-1/2D Electromagnetic MPI/OpenMP PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "mpbpush2.h"
#include "mpplib2.h"

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
void cppdblkp2l(float part[], int kpic[], int npp, int noff, int *nppmx,
                int idimp, int npmax, int mx, int my, int mx1,
                int mxyp1, int *irc) {
/* this subroutine finds the maximum number of particles in each tile of
   mx, my to calculate size of segmented particle array ppart
   linear interpolation, spatial decomposition in y direction
   input: all except kpic, nppmx, output: kpic, nppmx
   part = input particle array
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   kpic = output number of particles per tile
   nppmx = return maximum number of particles in tile
   npp = number of particles in partition
   noff = backmost global gridpoint in particle partition
   idimp = size of phase space = 4
   npmax = maximum number of particles in each partition
   mx/my = number of grids in sorting cell in x and y
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int j, k, n, m, mnoff, isum, ist, npx, ierr;
   mnoff = noff;
   ierr = 0;
/* clear counter array */
   for (k = 0; k < mxyp1; k++) {
      kpic[k] = 0;
   }
/* find how many particles in each tile */
   for (j = 0; j < npp; j++) {
      n = part[idimp*j];
      m = part[1+idimp*j];
      n = n/mx;
      m = (m - mnoff)/my;
      m = n + mx1*m;
      if (m < mxyp1) {
         kpic[m] += 1;
      }
      else {
         ierr = ierr > m-mxyp1+1 ? ierr : m-mxyp1+1;
      }
   }
/* find maximum */
   isum = 0;
   npx = 0;
   for (k = 0; k < mxyp1; k++) {
      ist = kpic[k];
      npx = npx > ist ? npx : ist;
      isum += ist;
   }
   *nppmx = npx;
/* check for errors */
   if (ierr > 0) {
      *irc = ierr;
   }
   else if (isum != npp) {
      *irc = -1;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cpppmovin2l(float part[], float ppart[], int kpic[], int npp,
                 int noff, int nppmx, int idimp, int npmax, int mx,
                 int my, int mx1, int mxyp1, int *irc) {
/* this subroutine sorts particles by x,y grid in tiles of
   mx, my and copies to segmented array ppart
   linear interpolation, spatial decomposition in y direction
   input: all except ppart, kpic, output: ppart, kpic
   part/ppart = input/output particle arrays
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   kpic = output number of particles per tile
   nppmx = maximum number of particles in tile
   npp = number of particles in partition
   noff = backmost global gridpoint in particle partition
   idimp = size of phase space = 4
   npmax = maximum number of particles in each partition
   mx/my = number of grids in sorting cell in x and y
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int i, j, k, n, m, mnoff, ip, ierr;
   mnoff = noff;
   ierr = 0;
/* clear counter array */
   for (k = 0; k < mxyp1; k++) {
      kpic[k] = 0;
   }
/* find addresses of particles at each tile and reorder particles */
   for (j = 0; j < npp; j++) {
      n = part[idimp*j];
      m = part[1+idimp*j];
      n = n/mx;
      m = (m - mnoff)/my;
      m = n + mx1*m;
      ip = kpic[m];
      if (ip < nppmx) {
         for (i = 0; i < idimp; i++) {
            ppart[i+idimp*(ip+nppmx*m)] = part[i+idimp*j];         }
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
void cpppcheck2l(float ppart[], int kpic[], int noff, int nyp,
                 int idimp, int nppmx, int nx, int mx, int my, int mx1,
                 int myp1, int *irc) {
/* this subroutine performs a sanity check to make sure particles sorted
   by x,y grid in tiles of mx, my, are all within bounds.
   tiles are assumed to be arranged in 2D linear memory
   input: all except irc
   output: irc
   ppart[k][n][0] = position x of particle n in tile k
   ppart[k][n][1] = position y of particle n in tile k
   kpic[k] = number of reordered output particles in tile k
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx = system length in x direction
   mx/my = number of grids in sorting cell in x/y
   mx1 = (system length in x direction - 1)/mx + 1
   myp1 = (partition length in y direction - 1)/my + 1
   irc = particle error, returned only if error occurs, when irc > 0
local data                                                            */
   int mxyp1, noffp, moffp, nppp, j, k, ist, nn, mm;
   float edgelx, edgely, edgerx, edgery, dx, dy;
   mxyp1 = mx1*myp1;
/* loop over tiles */
#pragma omp parallel for \
private(j,k,noffp,moffp,nppp,nn,mm,ist,edgelx,edgely,edgerx,edgery,dx, \
dy)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      nn = nx - noffp;
      nn = mx < nn ? mx : nn;
      mm = nyp - moffp;
      mm = my < mm ? my : mm;
      edgelx = noffp;
      edgerx = noffp + nn;
      edgely = noff + moffp;
      edgery = noff + moffp + mm;
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
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
void cppgbppush23l(float ppart[], float fxy[], float bxy[], int kpic[],
                   int noff, int nyp, float qbm, float dt, float dtc,
                   float *ek, int idimp, int nppmx, int nx, int ny,
                   int mx, int my, int nxv, int nypmx, int mx1,
                   int mxyp1, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field. Using the Boris Mover.
   OpenMP version using guard cells, for distributed data
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
   ppart[m][n][0] = position x of particle n in partition in tile m
   ppart[m][n][1] = position y of particle n in partition in tile m
   ppart[m][n][2] = x velocity of particle n in partition in tile m
   ppart[m][n][3] = y velocity of particle n in partition in tile m
   ppart[m][n][4] = z velocity of particle n in partition in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,kk)
   fxy[k][j][1] = y component of force/charge at grid (j,kk)
   fxy[k][j][2] = z component of force/charge at grid (j,kk)
   that is, convolution of electric field over particle shape,
   where kk = k + noff
   bxy[k][j][0] = x component of magnetic field at grid (j,kk)
   bxy[k][j][1] = y component of magnetic field at grid (j,kk)
   bxy[k][j][2] = z component of magnetic field at grid (j,kk)
   that is, the convolution of magnetic field over particle shape,
   where kk = k + noff
   kpic = number of particles per tile
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
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
   nxv = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
   int noffp, moffp, npoff, nppp, mxv3;
   int mnoff, i, j, k, nn, mm, nm;
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
   edgelx = 0.0f;
   edgely = 1.0f;
   edgerx = (float) (nx);
   edgery = (float) (ny-1);
   if ((ipbc==2) || (ipbc==3)) {
      edgelx = 1.0f;
      edgerx = (float) (nx-1);
   }
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noffp,moffp,nppp,npoff,nn,mm,nm,mnoff,x,y,dxp,dyp,amx, \
amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2, \
rot3,rot4,rot5,rot6,rot7,rot8,rot9,sum1,sfxy,sbxy) \
reduction(+:sum2)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      mnoff = moffp + noff;
      npoff = nppmx*k;
/* load local fields from global array */
      nn = (mx < nx-noffp ? mx : nx-noffp) + 1;
      mm = (my < nyp-moffp ? my : nyp-moffp) + 1;
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[3*i+mxv3*j] = fxy[3*(i+noffp+nxv*(j+moffp))];
            sfxy[1+3*i+mxv3*j] = fxy[1+3*(i+noffp+nxv*(j+moffp))];
            sfxy[2+3*i+mxv3*j] = fxy[2+3*(i+noffp+nxv*(j+moffp))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[3*i+mxv3*j] = bxy[3*(i+noffp+nxv*(j+moffp))];
            sbxy[1+3*i+mxv3*j] = bxy[1+3*(i+noffp+nxv*(j+moffp))];
            sbxy[2+3*i+mxv3*j] = bxy[2+3*(i+noffp+nxv*(j+moffp))];
         }
      }
      sum1 = 0.0;
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 3*(nn - noffp) + mxv3*(mm - mnoff);
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
void cppgbppushf23l(float ppart[], float fxy[], float bxy[], int kpic[],
                    int ncl[], int ihole[], int noff, int nyp,
                    float qbm, float dt, float dtc, float *ek,
                    int idimp, int nppmx, int nx, int ny,
                    int mx, int my, int nxv, int nypmx, int mx1,
                    int mxyp1, int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field. Using the Boris Mover.
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells, for distributed data
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
   ppart[m][n][0] = position x of particle n in partition in tile m
   ppart[m][n][1] = position y of particle n in partition in tile m
   ppart[m][n][2] = x velocity of particle n in partition in tile m
   ppart[m][n][3] = y velocity of particle n in partition in tile m
   ppart[m][n][4] = z velocity of particle n in partition in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,kk)
   fxy[k][j][1] = y component of force/charge at grid (j,kk)
   fxy[k][j][2] = z component of force/charge at grid (j,kk)
   that is, convolution of electric field over particle shape,
   where kk = k + noff
   bxy[k][j][0] = x component of magnetic field at grid (j,kk)
   bxy[k][j][1] = y component of magnetic field at grid (j,kk)
   bxy[k][j][2] = z component of magnetic field at grid (j,kk)
   that is, the convolution of magnetic field over particle shape,
   where kk = k + noff
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
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
   nxv = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             33
#define MYV             33
   int noffp, moffp, npoff, nppp, mxv3;
   int mnoff, i, j, k, ih, nh, nn, mm, nm;
   float qtmh, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt;
   float anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
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
private(i,j,k,noffp,moffp,nppp,npoff,nn,mm,nm,ih,nh,mnoff,x,y,dxp,dyp, \
amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1, \
rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,edgery, \
sum1,sfxy,sbxy) \
reduction(+:sum2)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      nn = nx - noffp;
      nn = mx < nn ? mx : nn;
      mm = nyp - moffp;
      mm = my < mm ? my : mm;
      edgelx = noffp;
      edgerx = noffp + nn;
      edgely = noff + moffp;
      edgery = noff + moffp + mm;
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
      mnoff = moffp + noff;
      npoff = nppmx*k;
/* load local fields from global array */
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[3*i+mxv3*j] = fxy[3*(i+noffp+nxv*(j+moffp))];
            sfxy[1+3*i+mxv3*j] = fxy[1+3*(i+noffp+nxv*(j+moffp))];
            sfxy[2+3*i+mxv3*j] = fxy[2+3*(i+noffp+nxv*(j+moffp))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[3*i+mxv3*j] = bxy[3*(i+noffp+nxv*(j+moffp))];
            sbxy[1+3*i+mxv3*j] = bxy[1+3*(i+noffp+nxv*(j+moffp))];
            sbxy[2+3*i+mxv3*j] = bxy[2+3*(i+noffp+nxv*(j+moffp))];
         }
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
      sum1 = 0.0;
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 3*(nn - noffp) + mxv3*(mm - mnoff);
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
void cppgrbppush23l(float ppart[], float fxy[], float bxy[], int kpic[],
                    int noff, int nyp, float qbm, float dt, float dtc,
                    float ci, float *ek, int idimp, int nppmx, int nx,
                    int ny, int mx, int my, int nxv, int nypmx, int mx1,
                    int mxyp1, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   OpenMP version using guard cells, for distributed data
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
   ppart[m][n][0] = position x of particle n in partition in tile m
   ppart[m][n][1] = position y of particle n in partition in tile m
   ppart[m][n][2] = x momentum of particle n in partition in tile m
   ppart[m][n][3] = y momentum of particle n in partition in tile m
   ppart[m][n][4] = z momentum of particle n in partition in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,kk)
   fxy[k][j][1] = y component of force/charge at grid (j,kk)
   fxy[k][j][2] = z component of force/charge at grid (j,kk)
   that is, convolution of electric field over particle shape,
   where kk = k + noff
   bxy[k][j][0] = x component of magnetic field at grid (j,kk)
   bxy[k][j][1] = y component of magnetic field at grid (j,kk)
   bxy[k][j][2] = z component of magnetic field at grid (j,kk)
   that is, the convolution of magnetic field over particle shape,
   where kk = k + noff
   kpic = number of particles per tile
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprical of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   idimp = size of phase space = 5
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
   int noffp, moffp, npoff, nppp, mxv3;
   int mnoff, i, j, k, nn, mm, nm;
   float qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg;
   float omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y;
   float sfxy[3*MXV*MYV], sbxy[3*MXV*MYV];
/* float sfxy[3*(mx+1)*(my+1)], sbxy[3*(mx+1)*(my+1)]; */
   double sum1, sum2;
   mxv3 = 3*(mx + 1);
   qtmh = 0.5*qbm*dt;
   ci2 = ci*ci;
   sum2 = 0.0;
/* set boundary values */
   edgelx = 0.0f;
   edgely = 1.0f;
   edgerx = (float) (nx);
   edgery = (float) (ny-1);
   if ((ipbc==2) || (ipbc==3)) {
      edgelx = 1.0f;
      edgerx = (float) (nx-1);
   }
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noffp,moffp,nppp,npoff,nn,mm,nm,mnoff,x,y,dxp,dyp,amx, \
amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2, \
rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,dtg,sum1,sfxy,sbxy) \
reduction(+:sum2)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      mnoff = moffp + noff;
      npoff = nppmx*k;
/* load local fields from global array */
      nn = (mx < nx-noffp ? mx : nx-noffp) + 1;
      mm = (my < nyp-moffp ? my : nyp-moffp) + 1;
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[3*i+mxv3*j] = fxy[3*(i+noffp+nxv*(j+moffp))];
            sfxy[1+3*i+mxv3*j] = fxy[1+3*(i+noffp+nxv*(j+moffp))];
            sfxy[2+3*i+mxv3*j] = fxy[2+3*(i+noffp+nxv*(j+moffp))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[3*i+mxv3*j] = bxy[3*(i+noffp+nxv*(j+moffp))];
            sbxy[1+3*i+mxv3*j] = bxy[1+3*(i+noffp+nxv*(j+moffp))];
            sbxy[2+3*i+mxv3*j] = bxy[2+3*(i+noffp+nxv*(j+moffp))];
         }
      }
      sum1 = 0.0;
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 3*(nn - noffp) + mxv3*(mm - mnoff);
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
/* find inverse gamma */
         p2 = acx*acx + acy*acy + acz*acz;
         gami = 1.0/sqrtf(1.0 + p2*ci2);
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
         ppart[2+idimp*(j+npoff)] = dx;
         ppart[3+idimp*(j+npoff)] = dy;
         ppart[4+idimp*(j+npoff)] = dz;
/* update inverse gamma */
         p2 = dx*dx + dy*dy + dz*dz;
         dtg = dtc/sqrtf(1.0 + p2*ci2);
/* new position */
         dx = x + dx*dtg;
         dy = y + dy*dtg;
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
   *ek += sum2;
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cppgrbppushf23l(float ppart[], float fxy[], float bxy[],
                     int kpic[], int ncl[], int ihole[], int noff,
                     int nyp, float qbm, float dt, float dtc, float ci,
                     float *ek, int idimp, int nppmx, int nx, int ny,
                     int mx, int my, int nxv, int nypmx, int mx1,
                     int mxyp1, int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells, for distributed data
   data deposited in tiles
   particles stored segmented array
   131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
   119 flops/particle, 1 divide, 29 loads, 5 stores
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
   ppart[m][n][0] = position x of particle n in partition in tile m
   ppart[m][n][1] = position y of particle n in partition in tile m
   ppart[m][n][2] = x momentum of particle n in partition in tile m
   ppart[m][n][3] = y momentum of particle n in partition in tile m
   ppart[m][n][4] = z momentum of particle n in partition in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,kk)
   fxy[k][j][1] = y component of force/charge at grid (j,kk)
   fxy[k][j][2] = z component of force/charge at grid (j,kk)
   that is, convolution of electric field over particle shape,
   where kk = k + noff
   bxy[k][j][0] = x component of magnetic field at grid (j,kk)
   bxy[k][j][1] = y component of magnetic field at grid (j,kk)
   bxy[k][j][2] = z component of magnetic field at grid (j,kk)
   that is, the convolution of magnetic field over particle shape,
   where kk = k + noff
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprical of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   idimp = size of phase space = 5
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             33
#define MYV             33
   int noffp, moffp, npoff, nppp, mxv3;
   int mnoff, i, j, k, ih, nh, nn, mm, nm;
   float qtmh, ci2, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg;
   float omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float x, y;
   float sfxy[3*MXV*MYV], sbxy[3*MXV*MYV];
/* float sfxy[3*(mx+1)*(my+1)], sbxy[3*(mx+1)*(my+1)]; */
   double sum1, sum2;
   mxv3 = 3*(mx + 1);
   qtmh = 0.5*qbm*dt;
   ci2 = ci*ci;
   anx = (float) nx;
   any = (float) ny;
   sum2 = 0.0;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noffp,moffp,nppp,npoff,nn,mm,nm,ih,nh,mnoff,x,y,dxp,dyp, \
amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1, \
rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,edgery, \
p2,gami,qtmg,dtg,sum1,sfxy,sbxy) \
reduction(+:sum2)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      nn = nx - noffp;
      nn = mx < nn ? mx : nn;
      mm = nyp - moffp;
      mm = my < mm ? my : mm;
      edgelx = noffp;
      edgerx = noffp + nn;
      edgely = noff + moffp;
      edgery = noff + moffp + mm;
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
      mnoff = moffp + noff;
      npoff = nppmx*k;
/* load local fields from global array */
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[3*i+mxv3*j] = fxy[3*(i+noffp+nxv*(j+moffp))];
            sfxy[1+3*i+mxv3*j] = fxy[1+3*(i+noffp+nxv*(j+moffp))];
            sfxy[2+3*i+mxv3*j] = fxy[2+3*(i+noffp+nxv*(j+moffp))];
         }
      }
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sbxy[3*i+mxv3*j] = bxy[3*(i+noffp+nxv*(j+moffp))];
            sbxy[1+3*i+mxv3*j] = bxy[1+3*(i+noffp+nxv*(j+moffp))];
            sbxy[2+3*i+mxv3*j] = bxy[2+3*(i+noffp+nxv*(j+moffp))];
         }
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
      sum1 = 0.0;
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 3*(nn - noffp) + mxv3*(mm - mnoff);
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
/* find inverse gamma */
         p2 = acx*acx + acy*acy + acz*acz;
         gami = 1.0/sqrtf(1.0 + p2*ci2);
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
         ppart[2+idimp*(j+npoff)] = dx;
         ppart[3+idimp*(j+npoff)] = dy;
         ppart[4+idimp*(j+npoff)] = dz;
/* update inverse gamma */
         p2 = dx*dx + dy*dy + dz*dz;
         dtg = dtc/sqrtf(1.0 + p2*ci2);
/* new position */
         dx = x + dx*dtg;
         dy = y + dy*dtg;
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
   *ek += sum2;
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cppgppost2l(float ppart[], float q[], int kpic[], int noff, 
                 float qm, int idimp, int nppmx, int mx, int my,
                 int nxv, int nypmx, int mx1, int mxyp1) {
/* for 2d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   OpenMP version using guard cells, for distributed data
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
   ppart[m][n][0] = position x of particle n in partition in tile m
   ppart[m][n][1] = position y of particle n in partition in tile m
   q[k][j] = charge density at grid point (j,kk),
   where kk = k + noff
   kpic = number of particles per tile
   noff = lowermost global gridpoint in particle partition.
   qm = charge on particle, in units of e
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of charge array, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
local data                                                            */
#define MXV             33
#define MYV             33
   int noffp, moffp, npoff, nppp, mxv;
   int mnoff, i, j, k, nn, mm;
   float x, y, dxp, dyp, amx, amy;
   float sq[MXV*MYV];
/* float sq[(mx+1)*(my+1)]; */
   mxv = mx + 1;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noffp,moffp,nppp,npoff,mnoff,nn,mm,x,y,dxp,dyp,amx,amy, \
sq)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      npoff = nppmx*k;
      mnoff = moffp + noff;
/* zero out local accumulator */
      for (j = 0; j < my+1; j++) {
         for (i = 0; i < mx+1; i++) {
            sq[i+mxv*j] = 0.0f;
         }
      }
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         nn = nn - noffp + mxv*(mm - mnoff);
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
      nn = nxv - noffp;
      mm = nypmx - moffp;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            q[i+noffp+nxv*(j+moffp)] += sq[i+mxv*j];
         }
      }
/* deposit charge to edge points in global array */
      mm = nypmx - moffp;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         q[i+noffp+nxv*moffp] += sq[i];
         if (mm > my) {
#pragma omp atomic
            q[i+noffp+nxv*(mm+moffp-1)] += sq[i+mxv*(mm-1)];
         }
      }
      nn = nxv - noffp;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         q[noffp+nxv*(j+moffp)] += sq[mxv*j];
         if (nn > mx) {
#pragma omp atomic
            q[nn+noffp-1+nxv*(j+moffp)] += sq[nn-1+mxv*j];
         }
      }
   }
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cppgjppost2l(float ppart[], float cu[], int kpic[], int noff,
                  float qm, float dt, int nppmx, int idimp, int nx,
                  int ny, int mx, int my, int nxv, int nypmx, int mx1,
                  int mxyp1, int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   OpenMP version using guard cells, for distributed data
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
   ppart[m][n][0] = position x of particle n in partition in tile m
   ppart[m][n][1] = position y of particle n in partition in tile m
   ppart[m][n][2] = x velocity of particle n in partition in tile m
   ppart[m][n][3] = y velocity of particle n in partition in tile m
   ppart[m][n][4] = z velocity of particle n in partition in tile m
   cu[k][j][i] = ith component of current density at grid point (j,kk),
   where kk = k + noff
   kpic = number of particles per tile
   noff = lowermost global gridpoint in particle partition.
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of current array, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
   int noffp, moffp, npoff, nppp, mxv3;
   int mnoff, i, j, k, nn, mm;
   float edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz;
   float scu[3*MXV*MYV];
/* float scu[3*(mx+1)*(my+1)]; */
   mxv3 = 3*(mx + 1);
/* set boundary values */
   edgelx = 0.0f;
   edgely = 1.0f;
   edgerx = (float) (nx);
   edgery = (float) (ny-1);
   if ((ipbc==2) || (ipbc==3)) {
      edgelx = 1.0f;
      edgerx = (float) (nx-1);
   }
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noffp,moffp,nppp,npoff,nn,mm,mnoff,x,y,dxp,dyp,amx,amy, \
dx,dy,vx,vy,vz,scu)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      mnoff = moffp + noff;
      npoff = nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < mxv3*(my+1); j++) {
         scu[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         nn = 3*(nn - noffp) + mxv3*(mm - mnoff);
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
      nn = nxv - noffp;
      mm = nypmx - moffp;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[3*(i+noffp+nxv*(j+moffp))] += scu[3*i+mxv3*j];
            cu[1+3*(i+noffp+nxv*(j+moffp))] += scu[1+3*i+mxv3*j];
            cu[2+3*(i+noffp+nxv*(j+moffp))] += scu[2+3*i+mxv3*j];
         }
      }
/* deposit current to edge points in global array */
      mm = nypmx - moffp;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[3*(i+noffp+nxv*moffp)] += scu[3*i];
#pragma omp atomic
         cu[1+3*(i+noffp+nxv*moffp)] += scu[1+3*i];
#pragma omp atomic
         cu[2+3*(i+noffp+nxv*moffp)] += scu[2+3*i];
         if (mm > my) {
#pragma omp atomic
            cu[3*(i+noffp+nxv*(mm+moffp-1))] += scu[3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[1+3*(i+noffp+nxv*(mm+moffp-1))] += scu[1+3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[2+3*(i+noffp+nxv*(mm+moffp-1))] += scu[2+3*i+mxv3*(mm-1)];
         }
      }
      nn = nxv - noffp;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[3*(noffp+nxv*(j+moffp))] += scu[mxv3*j];
#pragma omp atomic
         cu[1+3*(noffp+nxv*(j+moffp))] += scu[1+mxv3*j];
#pragma omp atomic
         cu[2+3*(noffp+nxv*(j+moffp))] += scu[2+mxv3*j];
         if (nn > mx) {
#pragma omp atomic
            cu[3*(nn+noffp-1+nxv*(j+moffp))] += scu[3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[1+3*(nn+noffp-1+nxv*(j+moffp))] += scu[1+3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[2+3*(nn+noffp-1+nxv*(j+moffp))] += scu[2+3*(nn-1)+mxv3*j];
         }
      }
   }
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cppgjppostf2l(float ppart[], float cu[], int kpic[], int ncl[],
                   int ihole[], int noff, int nyp, float qm, float dt,
                   int nppmx, int idimp, int nx, int ny, int mx, int my,
                   int nxv, int nypmx, int mx1, int mxyp1, int ntmax,
                   int *irc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells, for distributed data
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
   ppart[m][n][0] = position x of particle n in partition in tile m
   ppart[m][n][1] = position y of particle n in partition in tile m
   ppart[m][n][2] = x velocity of particle n in partition in tile m
   ppart[m][n][3] = y velocity of particle n in partition in tile m
   ppart[m][n][4] = z velocity of particle n in partition in tile m
   cu[k][j][i] = ith component of current density at grid point (j,kk),
   where kk = k + noff
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of current array, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             33
#define MYV             33
   int noffp, moffp, npoff, nppp, mxv3;
   int mnoff, i, j, k, ih, nh, nn, mm;
   float dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float scu[3*MXV*MYV];
/* float scu[3*(mx+1)*(my+1)]; */
   mxv3 = 3*(mx + 1);
   anx = (float) nx;
   any = (float) ny;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noffp,moffp,nppp,npoff,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx, \
amy,dx,dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,scu)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      nn = nx - noffp;
      nn = mx < nn ? mx : nn;
      mm = nyp - moffp;
      mm = my < mm ? my : mm;
      edgelx = noffp;
      edgerx = noffp + nn;
      edgely = noff + moffp;
      edgery = noff + moffp + mm;
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
      mnoff = moffp + noff;
      npoff = nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < mxv3*(my+1); j++) {
         scu[j] = 0.0f;
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         nn = 3*(nn - noffp) + mxv3*(mm - mnoff);
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
/* deposit current to interior points in global array */
      nn = nxv - noffp;
      mm = nypmx - moffp;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[3*(i+noffp+nxv*(j+moffp))] += scu[3*i+mxv3*j];
            cu[1+3*(i+noffp+nxv*(j+moffp))] += scu[1+3*i+mxv3*j];
            cu[2+3*(i+noffp+nxv*(j+moffp))] += scu[2+3*i+mxv3*j];
         }
      }
/* deposit current to edge points in global array */
      mm = nypmx - moffp;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[3*(i+noffp+nxv*moffp)] += scu[3*i];
#pragma omp atomic
         cu[1+3*(i+noffp+nxv*moffp)] += scu[1+3*i];
#pragma omp atomic
         cu[2+3*(i+noffp+nxv*moffp)] += scu[2+3*i];
         if (mm > my) {
#pragma omp atomic
            cu[3*(i+noffp+nxv*(mm+moffp-1))] += scu[3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[1+3*(i+noffp+nxv*(mm+moffp-1))] += scu[1+3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[2+3*(i+noffp+nxv*(mm+moffp-1))] += scu[2+3*i+mxv3*(mm-1)];
         }
      }
      nn = nxv - noffp;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[3*(noffp+nxv*(j+moffp))] += scu[mxv3*j];
#pragma omp atomic
         cu[1+3*(noffp+nxv*(j+moffp))] += scu[1+mxv3*j];
#pragma omp atomic
         cu[2+3*(noffp+nxv*(j+moffp))] += scu[2+mxv3*j];
         if (nn > mx) {
#pragma omp atomic
            cu[3*(nn+noffp-1+nxv*(j+moffp))] += scu[3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[1+3*(nn+noffp-1+nxv*(j+moffp))] += scu[1+3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[2+3*(nn+noffp-1+nxv*(j+moffp))] += scu[2+3*(nn-1)+mxv3*j];
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
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cppgrjppost2l(float ppart[], float cu[], int kpic[], int noff,
                   float qm, float dt, float ci, int nppmx, int idimp,
                   int nx, int ny, int mx, int my, int nxv, int nypmx,
                   int mx1, int mxyp1, int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation for relativistic particles
   in addition, particle positions are advanced a half time-step
   OpenMP version using guard cells, for distributed data
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
   ppart[m][n][0] = position x of particle n in partition in tile m
   ppart[m][n][1] = position y of particle n in partition in tile m
   ppart[m][n][2] = x momentum of particle n in partition in tile m
   ppart[m][n][3] = y momentum of particle n in partition in tile m
   ppart[m][n][4] = z momentum of particle n in partition in tile m
   cu[k][j][i] = ith component of current density at grid point (j,kk),
   where kk = k + noff
   kpic = number of particles per tile
   noff = lowermost global gridpoint in particle partition.
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprical of velocity of light
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of current array, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
   int noffp, moffp, npoff, nppp, mxv3;
   int mnoff, i, j, k, nn, mm;
   float ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz, p2, gami;
   float scu[3*MXV*MYV];
/* float scu[3*(mx+1)*(my+1)]; */
   mxv3 = 3*(mx + 1);
   ci2 = ci*ci;
/* set boundary values */
   edgelx = 0.0f;
   edgely = 1.0f;
   edgerx = (float) (nx);
   edgery = (float) (ny-1);
   if ((ipbc==2) || (ipbc==3)) {
      edgelx = 1.0f;
      edgerx = (float) (nx-1);
   }
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noffp,moffp,nppp,npoff,nn,mm,mnoff,x,y,dxp,dyp,amx,amy, \
dx,dy,vx,vy,vz,p2,gami,scu)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      mnoff = moffp + noff;
      npoff = nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < mxv3*(my+1); j++) {
         scu[j] = 0.0f;
      }
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
/* find inverse gamma */
         vx = ppart[2+idimp*(j+npoff)];
         vy = ppart[3+idimp*(j+npoff)];
         vz = ppart[4+idimp*(j+npoff)];
         p2 = vx*vx + vy*vy + vz*vz;
         gami = 1.0/sqrtf(1.0 + p2*ci2);
/* calculate weights */
         nn = 3*(nn - noffp) + mxv3*(mm - mnoff);
         amx = qm - dxp;
         amy = 1.0 - dyp;
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
      nn = nxv - noffp;
      mm = nypmx - moffp;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[3*(i+noffp+nxv*(j+moffp))] += scu[3*i+mxv3*j];
            cu[1+3*(i+noffp+nxv*(j+moffp))] += scu[1+3*i+mxv3*j];
            cu[2+3*(i+noffp+nxv*(j+moffp))] += scu[2+3*i+mxv3*j];
         }
      }
/* deposit current to edge points in global array */
      mm = nypmx - moffp;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[3*(i+noffp+nxv*moffp)] += scu[3*i];
#pragma omp atomic
         cu[1+3*(i+noffp+nxv*moffp)] += scu[1+3*i];
#pragma omp atomic
         cu[2+3*(i+noffp+nxv*moffp)] += scu[2+3*i];
         if (mm > my) {
#pragma omp atomic
            cu[3*(i+noffp+nxv*(mm+moffp-1))] += scu[3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[1+3*(i+noffp+nxv*(mm+moffp-1))] += scu[1+3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[2+3*(i+noffp+nxv*(mm+moffp-1))] += scu[2+3*i+mxv3*(mm-1)];
         }
      }
      nn = nxv - noffp;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[3*(noffp+nxv*(j+moffp))] += scu[mxv3*j];
#pragma omp atomic
         cu[1+3*(noffp+nxv*(j+moffp))] += scu[1+mxv3*j];
#pragma omp atomic
         cu[2+3*(noffp+nxv*(j+moffp))] += scu[2+mxv3*j];
         if (nn > mx) {
#pragma omp atomic
            cu[3*(nn+noffp-1+nxv*(j+moffp))] += scu[3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[1+3*(nn+noffp-1+nxv*(j+moffp))] += scu[1+3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[2+3*(nn+noffp-1+nxv*(j+moffp))] += scu[2+3*(nn-1)+mxv3*j];
         }
      }
   }
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cppgrjppostf2l(float ppart[], float cu[], int kpic[], int ncl[],
                    int ihole[], int noff, int nyp, float qm, float dt,
                    float ci, int nppmx, int idimp, int nx, int ny,
                    int mx, int my, int nxv, int nypmx, int mx1,
                    int mxyp1, int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation for relativistic particles
   in addition, particle positions are advanced a half time-step
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells, for distributed data
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
   ppart[m][n][0] = position x of particle n in partition in tile m
   ppart[m][n][1] = position y of particle n in partition in tile m
   ppart[m][n][2] = x momentum of particle n in partition in tile m
   ppart[m][n][3] = y momentum of particle n in partition in tile m
   ppart[m][n][4] = z momentum of particle n in partition in tile m
   cu[k][j][i] = ith component of current density at grid point (j,kk),
   where kk = k + noff
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprical of velocity of light
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of current array, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             33
#define MYV             33
   int noffp, moffp, npoff, nppp, mxv3;
   int mnoff, i, j, k, ih, nh, nn, mm;
   float ci2, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz, p2, gami;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float scu[3*MXV*MYV];
/* float scu[3*(mx+1)*(my+1)]; */
   mxv3 = 3*(mx + 1);
   ci2 = ci*ci;
   anx = (float) nx;
   any = (float) ny;
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noffp,moffp,nppp,npoff,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx, \
amy,dx,dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,p2,gami,scu)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      nn = nx - noffp;
      nn = mx < nn ? mx : nn;
      mm = nyp - moffp;
      mm = my < mm ? my : mm;
      edgelx = noffp;
      edgerx = noffp + nn;
      edgely = noff + moffp;
      edgery = noff + moffp + mm;
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
      mnoff = moffp + noff;
      npoff = nppmx*k;
/* zero out local accumulator */
      for (j = 0; j < mxv3*(my+1); j++) {
         scu[j] = 0.0f;
      }
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
/* find interpolation weights */
         x = ppart[idimp*(j+npoff)];
         y = ppart[1+idimp*(j+npoff)];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
/* find inverse gamma */
         vx = ppart[2+idimp*(j+npoff)];
         vy = ppart[3+idimp*(j+npoff)];
         vz = ppart[4+idimp*(j+npoff)];
         p2 = vx*vx + vy*vy + vz*vz;
         gami = 1.0/sqrtf(1.0 + p2*ci2);
/* calculate weights */
         nn = 3*(nn - noffp) + mxv3*(mm - mnoff);
         amx = qm - dxp;
         amy = 1.0 - dyp;
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
/* deposit current to interior points in global array */
      nn = nxv - noffp;
      mm = nypmx - moffp;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
            cu[3*(i+noffp+nxv*(j+moffp))] += scu[3*i+mxv3*j];
            cu[1+3*(i+noffp+nxv*(j+moffp))] += scu[1+3*i+mxv3*j];
            cu[2+3*(i+noffp+nxv*(j+moffp))] += scu[2+3*i+mxv3*j];
         }
      }
/* deposit current to edge points in global array */
      mm = nypmx - moffp;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[3*(i+noffp+nxv*moffp)] += scu[3*i];
#pragma omp atomic
         cu[1+3*(i+noffp+nxv*moffp)] += scu[1+3*i];
#pragma omp atomic
         cu[2+3*(i+noffp+nxv*moffp)] += scu[2+3*i];
         if (mm > my) {
#pragma omp atomic
            cu[3*(i+noffp+nxv*(mm+moffp-1))] += scu[3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[1+3*(i+noffp+nxv*(mm+moffp-1))] += scu[1+3*i+mxv3*(mm-1)];
#pragma omp atomic
            cu[2+3*(i+noffp+nxv*(mm+moffp-1))] += scu[2+3*i+mxv3*(mm-1)];
         }
      }
      nn = nxv - noffp;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[3*(noffp+nxv*(j+moffp))] += scu[mxv3*j];
#pragma omp atomic
         cu[1+3*(noffp+nxv*(j+moffp))] += scu[1+mxv3*j];
#pragma omp atomic
         cu[2+3*(noffp+nxv*(j+moffp))] += scu[2+mxv3*j];
         if (nn > mx) {
#pragma omp atomic
            cu[3*(nn+noffp-1+nxv*(j+moffp))] += scu[3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[1+3*(nn+noffp-1+nxv*(j+moffp))] += scu[1+3*(nn-1)+mxv3*j];
#pragma omp atomic
            cu[2+3*(nn+noffp-1+nxv*(j+moffp))] += scu[2+3*(nn-1)+mxv3*j];
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
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cppporder2la(float ppart[], float ppbuff[], float sbufl[],
                  float sbufr[], int kpic[], int ncl[], int ihole[],
                  int ncll[], int nclr[], int noff, int nyp, int idimp,
                  int nppmx, int nx, int ny, int mx, int my, int mx1,
                  int myp1, int npbmx, int ntmax, int nbmax, int *irc) {
/* this subroutine performs first part of a particle sort by x,y grid
   in tiles of mx, my
   linear interpolation, with periodic boundary conditions
   for distributed data, with 1d domain decomposition in y.
   tiles are assumed to be arranged in 2D linear memory
   this part of the algorithm has 3 steps.  first, one finds particles
   leaving tile and stores their number in each directon, location, and
   destination in ncl and ihole.  then, a prefix scan of ncl is performed
   and departing particles are buffered in ppbuff in direction order.
   finally, we buffer particles leaving the processor in sbufl and sbufr,
   and store particle number offsets in ncll and nclr.
   input: all except ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
   output: ppart, ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
   ppart[k][n][0] = position x of particle n in tile k
   ppart[k][n][1] = position y of particle n in tile k 
   ppbuff[k][n][i] = i co-ordinate of particle n in tile k
   sbufl = buffer for particles being sent to lower processor
   sbufr = buffer for particles being sent to upper processor
   kpic[k] = number of particles in tile k
   ncl(i,k) = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = direction destination of particle leaving hole
   all for tile k
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   ncll = number offset being sent to lower processor
   nclr = number offset being sent to upper processor
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   mx1 = (system length in x direction - 1)/mx + 1
   myp1 = (partition length in y direction - 1)/my + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   nbmax =  size of buffers for passing particles between processors
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int mxyp1, noffp, moffp, nppp;
   int i, j, k, ii, jj, ih, nh, ist, nn, mm, isum, ip, j1, kk;
   float anx, any, edgelx, edgely, edgerx, edgery, dx, dy;
   mxyp1 = mx1*myp1;
   anx = (float) nx;
   any = (float) ny;
/* find and count particles leaving tiles and determine destination */
/* update ppart, ihole, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(j,k,noffp,moffp,nppp,nn,mm,ih,nh,ist,dx,dy,edgelx,edgely, \
edgerx,edgery)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      nn = nx - noffp;
      nn = mx < nn ? mx : nn;
      mm = nyp - moffp;
      mm = my < mm ? my : mm;
      ih = 0;
      nh = 0;
      edgelx = noffp;
      edgerx = noffp + nn;
      edgely = noff + moffp;
      edgery = noff + moffp + mm;
/* clear counters */
      for (j = 0; j < 8; j++) {
         ncl[j+8*k] = 0;
      }
/* loop over particles in tile */
      for (j = 0; j < nppp; j++) {
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
   for (k = 0; k < mxyp1; k++) {
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

/* buffer particles and their number leaving the node: */
/* update sbufl, sbufr, ncll, nclr */
   kk = mx1*(myp1 - 1);
#pragma omp parallel for private(k)
   for (k = 0; k < mx1; k++) {
      ncll[3*k] = ncl[4+8*k] - ncl[1+8*k];
      nclr[3*k] = ncl[7+8*(k+kk)] - ncl[4+8*(k+kk)];
   }
/* perform prefix scan */
   kk = 1;
L90: if (kk >= mx1)
      goto L110;
#pragma omp parallel for private(k,ii,nn,mm)
   for (k = 0; k < mx1; k++) {
      ii = k/kk;
      nn = kk*ii;
      mm = 2*nn + kk - 1;
      nn += k + kk;
      if (nn < mx1) {
         ncll[3*nn] += ncll[3*mm];
         nclr[3*nn] += nclr[3*mm];
      }
   }
   kk += kk;
   goto L90;
L110: kk = mx1*(myp1 - 1);
#pragma omp parallel for private(i,j,k,ii,nn,mm)
   for (k = 0; k < mx1; k++) {
      ii = ncl[4+8*k] - ncl[1+8*k];
      nn = ncll[3*k] - ii;
      jj = nbmax - nn;
      jj = ii < jj ? ii : jj;
      for (j = 0; j < jj; j++) {
         for (i = 0; i < idimp; i++) {
            sbufl[i+idimp*(j+nn)]
            = ppbuff[i+idimp*(j+ncl[1+8*k]+npbmx*k)];
         }
      }
      for (i = 0; i < 3; i++) {
         ncll[i+3*k] = ncl[i+2+8*k] - ncl[1+8*k] + nn;
      }
      ii = ncl[7+8*(k+kk)] - ncl[4+8*(k+kk)];
      mm = nclr[3*k] - ii;
      jj = nbmax - mm;
      jj = ii < jj ? ii : jj;
      for (j = 0; j < jj; j++) {
         for (i = 0; i < idimp; i++) {
            sbufr[i+idimp*(j+mm)]
            = ppbuff[i+idimp*(j+ncl[4+8*(k+kk)]+npbmx*(k+kk))];
         }
      }
      for (i = 0; i < 3; i++) {
         nclr[i+3*k] = ncl[i+5+8*(k+kk)] - ncl[4+8*(k+kk)] + mm;
      }
   }
/* sbufl or sbufr overflow */
   nn = ncll[3*mx1-1];
   mm = nclr[3*mx1-1];
   ii = nn > mm ? nn : mm;
   if (ii > nbmax)
      *irc = ii;
   return;
}

/*--------------------------------------------------------------------*/
void cppporderf2la(float ppart[], float ppbuff[], float sbufl[],
                   float sbufr[], int ncl[], int ihole[], int ncll[],
                   int nclr[], int idimp, int nppmx, int mx1, int myp1,
                   int npbmx, int ntmax, int nbmax, int *irc) {
/* this subroutine performs first part of a particle sort by x,y grid
   in tiles of mx, my
   linear interpolation, with periodic boundary conditions
   for distributed data, with 1d domain decomposition in y.
   tiles are assumed to be arranged in 2D linear memory
   this part of the algorithm has 2 steps.  first, a prefix scan of ncl
   is performed and departing particles are buffered in ppbuff in
   direction order. then, we buffer particles leaving the processor in
   sbufl and sbufr, and store particle number offsets in ncll and nclr.
   it assumes that the number, location, and destination of particles 
   leaving a tile have been previously stored in ncl and ihole by the
   cppgppushf2l procedure.
   input: all except ppbuff, sbufl, sbufr, ncll, nclr, irc
   output: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
   ppart[k][n][0] = position x of particle n in tile k
   ppart[k][n][1] = position y of particle n in tile k 
   ppbuff[k][n][i] = i co-ordinate of particle n in tile k
   sbufl = buffer for particles being sent to lower processor
   sbufr = buffer for particles being sent to upper processor
   ncl(i,k) = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = direction destination of particle leaving hole
   all for tile k
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   ncll = number offset being sent to lower processor
   nclr = number offset being sent to upper processor
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   mx1 = (system length in x direction - 1)/mx + 1
   myp1 = (partition length in y direction - 1)/my + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   nbmax =  size of buffers for passing particles between processors
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int mxyp1;
   int i, j, k, ii, jj, nh, ist, nn, mm, isum, ip, j1, kk;
   mxyp1 = mx1*myp1;
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,isum,ist,nh,ip,j1,ii)
   for (k = 0; k < mxyp1; k++) {
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
 
/* buffer particles and their number leaving the node: */
/* update sbufl, sbufr, ncll, nclr */
   kk = mx1*(myp1 - 1);
#pragma omp parallel for private(k)
   for (k = 0; k < mx1; k++) {
      ncll[3*k] = ncl[4+8*k] - ncl[1+8*k];
      nclr[3*k] = ncl[7+8*(k+kk)] - ncl[4+8*(k+kk)];
   }
/* perform prefix scan */
   kk = 1;
L90: if (kk >= mx1)
      goto L110;
#pragma omp parallel for private(k,ii,nn,mm)
   for (k = 0; k < mx1; k++) {
      ii = k/kk;
      nn = kk*ii;
      mm = 2*nn + kk - 1;
      nn += k + kk;
      if (nn < mx1) {
         ncll[3*nn] += ncll[3*mm];
         nclr[3*nn] += nclr[3*mm];
      }
   }
   kk += kk;
   goto L90;
L110: kk = mx1*(myp1 - 1);
#pragma omp parallel for private(i,j,k,ii,nn,mm)
   for (k = 0; k < mx1; k++) {
      ii = ncl[4+8*k] - ncl[1+8*k];
      nn = ncll[3*k] - ii;
      jj = nbmax - nn;
      jj = ii < jj ? ii : jj;
      for (j = 0; j < jj; j++) {
         for (i = 0; i < idimp; i++) {
            sbufl[i+idimp*(j+nn)]
            = ppbuff[i+idimp*(j+ncl[1+8*k]+npbmx*k)];
         }
      }
      for (i = 0; i < 3; i++) {
         ncll[i+3*k] = ncl[i+2+8*k] - ncl[1+8*k] + nn;
      }
      ii = ncl[7+8*(k+kk)] - ncl[4+8*(k+kk)];
      mm = nclr[3*k] - ii;
      jj = nbmax - mm;
      jj = ii < jj ? ii : jj;
      for (j = 0; j < jj; j++) {
         for (i = 0; i < idimp; i++) {
            sbufr[i+idimp*(j+mm)]
            = ppbuff[i+idimp*(j+ncl[4+8*(k+kk)]+npbmx*(k+kk))];
         }
      }
      for (i = 0; i < 3; i++) {
         nclr[i+3*k] = ncl[i+5+8*(k+kk)] - ncl[4+8*(k+kk)] + mm;
      }
   }
/* sbufl or sbufr overflow */
   nn = ncll[3*mx1-1];
   mm = nclr[3*mx1-1];
   ii = nn > mm ? nn : mm;
   if (ii > nbmax)
      *irc = ii;
   return;
}

/*--------------------------------------------------------------------*/
void cppporder2lb(float ppart[], float ppbuff[], float rbufl[],
                  float rbufr[], int kpic[], int ncl[], int ihole[],
                  int mcll[], int mclr[], int idimp, int nppmx, int mx1,
                  int myp1, int npbmx, int ntmax, int nbmax, int *irc) {
/* this subroutine performs second part of a particle sort by x,y grid
   in tiles of mx, my
   linear interpolation, with periodic boundary conditions
   for distributed data, with 1d domain decomposition in y.
   tiles are assumed to be arranged in 2D linear memory
   incoming particles from other tiles are copied from ppbuff, rbufl, and
   rbufr into ppart
   input: all except ppart, kpic, irc
   output: ppart, kpic, irc
   ppart[k][n][0] = position x of particle n in tile k
   ppart[k][n][1] = position y of particle n in tile k 
   ppbuff[k][n][i] = i co-ordinate of particle n in tile k
   rbufl = buffer for particles being received from lower processor
   rbufr = buffer for particles being received from upper processor
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = direction destination of particle leaving hole
   all for tile k
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   mcll = number offset being received from lower processor
   mclr = number offset being received from upper processor
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   mx1 = (system length in x direction - 1)/mx + 1
   myp1 = (partition length in y direction - 1)/my + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   nbmax =  size of buffers for passing particles between processors
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int mxyp1, nppp, ncoff, noff, moff;
   int i, j, k, ii, kx, ky, ih, nh, ist;
   int ip, j1, j2, kxl, kxr, kk, kl, kr;
   int ks[8];
   mxyp1 = mx1*myp1;
/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,ii,kk,nppp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,noff,moff, \
ist,j1,j2,ip,ks)
   for (k = 0; k < mxyp1; k++) {
      nppp = kpic[k];
      ky = k/mx1;
/* loop over tiles in y */
      kk = ky*mx1;
/* find tile above */
      kl = (ky - 1)*mx1;
/* find tile below */
      kr = (ky + 1)*mx1;
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
      noff = 0;
      moff = 0;
      if (ky==0) {
         if (kx > 0)
            noff = mcll[2+3*(kx-1)];
      }
      if (ky==(myp1-1)) {
         if (kx > 0)
            moff = mclr[2+3*(kx-1)];
      }
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      for (ii = 0; ii < 8; ii++) {
/* ip = number of particles coming from direction ii */
         if (ks[ii] < 0) {
            if (ii > 5)
               noff = mcll[ii-6+3*(ks[ii]+mx1)];
            ip = mcll[ii-5+3*(ks[ii]+mx1)] - noff;
         }
         else if (ks[ii] >= mxyp1) {
            if (ii > 2)
               moff = mclr[ii-3+3*(ks[ii]-mxyp1)];
            ip = mclr[ii-2+3*(ks[ii]-mxyp1)] - moff;
         }
         else {
            if (ii > 0)
               ncoff = ncl[ii-1+8*ks[ii]];
            ip = ncl[ii+8*ks[ii]] - ncoff;
         }
         for (j = 0; j < ip; j++) {
            ih += 1;
/* insert incoming particles into holes */
            if (ih <= nh) {
               j1 = ihole[2*(ih+(ntmax+1)*k)] - 1;
            }
/* place overflow at end of array */
            else {
               j1 = nppp;
               nppp += 1;
            }
            if (j1 < nppmx) {
               if (ks[ii] < 0) {
                  for (i = 0; i < idimp; i++) {
                     ppart[i+idimp*(j1+nppmx*k)]
                     = rbufl[i+idimp*(j+noff)];
                  }
               }
               else if (ks[ii] >= mxyp1) {
                  for (i = 0; i < idimp; i++) {
                     ppart[i+idimp*(j1+nppmx*k)]
                     = rbufr[i+idimp*(j+moff)];
                  }
               }
               else {
                  for (i = 0; i < idimp; i++) {
                     ppart[i+idimp*(j1+nppmx*k)]
                     = ppbuff[i+idimp*(j+ncoff+npbmx*ks[ii])];
                  }
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
            j1 = nppp - j - 1;
            j2 = ihole[2*(nh-j+(ntmax+1)*k)] - 1;
            if (j1 > j2) {
/* move particle only if it is below current hole */
               for (i = 0; i < idimp; i++) {
                  ppart[i+idimp*(j2+nppmx*k)]
                  = ppart[i+idimp*(j1+nppmx*k)];
               }
            }
         }
         nppp -= ip;
      }
      kpic[k] = nppp;
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
void cmppois23(float complex q[], float complex fxy[], int isign,
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
   int nxh, nyh, ks, joff, kxps, j, jj, jk, jk3, k, k1;
   float dnx, dny, dkx, dky, at1, at2, at3, at4;
   float complex zero, zt1, zt2;
   double wp, sum1;
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
L30: sum1 = 0.0;
   if (kstrt > nxh)
      goto L70;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,jj,jk,jk3,dkx,at1,at2,at3,zt1,zt2,wp) \
reduction(+:sum1)
   for (j = 0; j < kxps; j++) {
      dkx = dnx*(float) (j + joff);
      jj = nyhd*j;
      jk = nyv*j;
      jk3 = 3*jk;
      wp = 0.0;
      if ((j+joff) > 0) {
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            at1 = crealf(ffc[k+jj])*cimagf(ffc[k+jj]);
            at2 = dkx*at1;
            at3 = dny*at1*(float) k;
            zt1 = cimagf(q[k+jk]) - crealf(q[k+jk])*_Complex_I;
            zt2 = cimagf(q[k1+jk]) - crealf(q[k1+jk])*_Complex_I;
            fxy[3*k+jk3] = at2*zt1;
            fxy[1+3*k+jk3] = at3*zt1;
            fxy[2+3*k+jk3] = zero;
            fxy[3*k1+jk3] = at2*zt2;
            fxy[1+3*k1+jk3] = -at3*zt2;
            fxy[2+3*k1+jk3] = zero;
            wp += at1*(q[k+jk]*conjf(q[k+jk])
                  + q[k1+jk]*conjf(q[k1+jk]));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         at1 = crealf(ffc[jj])*cimagf(ffc[jj]);
         at3 = dkx*at1;
         zt1 = cimagf(q[jk]) - crealf(q[jk])*_Complex_I;
         fxy[jk3] = at3*zt1;
         fxy[1+jk3] = zero;
         fxy[2+jk3] = zero;
         fxy[3*k1+jk3] = zero;
         fxy[1+3*k1+jk3] = zero;
         fxy[2+3*k1+jk3] = zero;
         wp += at1*(q[jk]*conjf(q[jk]));
      }
      sum1 += wp;
   }
   wp = 0.0;
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
   sum1 += wp;
L70:
   *we = sum1*((float) nx)*((float) ny);
   return;
}

/*--------------------------------------------------------------------*/
void cmppcuperp2(float complex cu[], int nx, int ny, int kstrt, int nyv,
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
   int nxh, nyh, ks, joff, kxps, j, jk3, k, k1;
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
#pragma omp parallel for private(j,k,k1,jk3,dkx,dkx2,dky,at1,zt1)
   for (j = 0; j < kxps; j++) {
      dkx = dnx*(float) (j + joff);
      dkx2 = dkx*dkx;
      jk3 = 3*nyv*j;
      if ((j+joff) > 0) {
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            dky = dny*(float) k;
            at1 = 1.0/(dky*dky + dkx2);
            zt1 = at1*(dkx*cu[3*k+jk3] + dky*cu[1+3*k+jk3]);
            cu[3*k+jk3] -= dkx*zt1;
            cu[1+3*k+jk3] -= dky*zt1;
            zt1 = at1*(dkx*cu[3*k1+jk3] - dky*cu[1+3*k1+jk3]);
            cu[3*k1+jk3] -= dkx*zt1;
            cu[1+3*k1+jk3] += dky*zt1;
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         cu[jk3] = zero;
         cu[3*k1+jk3] = zero;
         cu[1+3*k1+jk3] = zero;
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
void cmippbpoisp23(float complex cu[], float complex bxy[],
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
   aimag(ffc[k][j]) = finite-size particle shape factor s
   real(ffc[k][j])) = potential green's function g
   for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
   ci = reciprical of velocity of light
   magnetic field energy is also calculated, using
   wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
      |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
   affp = normalization constant = nx*ny/np, where np=number of particles
   this expression is valid only if the current is divergence-free
   nx/ny = system length in x/y direction
   nyv = first dimension of field arrays, must be >= ny
   nyhd = first dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, ks, joff, kxps, j, jj, jk3, k, k1;
   float ci2, dnx, dny, dkx, dky, at1, at2, at3;
   float complex zero, zt1, zt2, zt3;
   double wp, sum1;
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
   sum1 = 0.0;
   if (kstrt > nxh)
      goto L40;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,jj,jk3,dkx,dky,at1,at2,at3,zt1,zt2,zt3,wp) \
reduction(+:sum1)
   for (j = 0; j < kxps; j++) {
      dkx = dnx*(float) (j + joff);
      jj = nyhd*j;
      jk3 = 3*nyv*j;
      wp = 0.0;
      if ((j+joff) > 0) {
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            dky = dny*(float) k;
            at1 = ci2*crealf(ffc[k+jj]);
            at2 = dky*at1;
            at3 = dkx*at1;
            at1 = at1*cimagf(ffc[k+jj]);
            zt1 = -cimagf(cu[2+3*k+jk3])
                + crealf(cu[2+3*k+jk3])*_Complex_I;
            zt2 = -cimagf(cu[1+3*k+jk3])
                + crealf(cu[1+3*k+jk3])*_Complex_I;
            zt3 = -cimagf(cu[3*k+jk3])
                + crealf(cu[3*k+jk3])*_Complex_I;
            bxy[3*k+jk3] = at2*zt1;
            bxy[1+3*k+jk3] = -at3*zt1;
            bxy[2+3*k+jk3] = at3*zt2 - at2*zt3;
            zt1 = -cimagf(cu[2+3*k1+jk3])
                + crealf(cu[2+3*k1+jk3])*_Complex_I;
            zt2 = -cimagf(cu[1+3*k1+jk3])
                + crealf(cu[1+3*k1+jk3])*_Complex_I;
            zt3 = -cimagf(cu[3*k1+jk3])
                + crealf(cu[3*k1+jk3])*_Complex_I;
            bxy[3*k1+jk3] = -at2*zt1;
            bxy[1+3*k1+jk3] = -at3*zt1;
            bxy[2+3*k1+jk3] = at3*zt2 + at2*zt3;
            wp += at1*(cu[3*k+jk3]*conjf(cu[3*k+jk3])
               + cu[1+3*k+jk3]*conjf(cu[1+3*k+jk3])
               + cu[2+3*k+jk3]*conjf(cu[2+3*k+jk3])
               + cu[3*k1+jk3]*conjf(cu[3*k1+jk3])
               + cu[1+3*k1+jk3]*conjf(cu[1+3*k1+jk3])
               + cu[2+3*k1+jk3]*conjf(cu[2+3*k1+jk3]));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         at1 = ci2*crealf(ffc[jj]);
         at2 = dkx*at1;
         at1 = at1*cimagf(ffc[jj]);
         zt1 = -cimagf(cu[2+jk3])
             + crealf(cu[2+jk3])*_Complex_I;
         zt2 = -cimagf(cu[1+jk3])
             + crealf(cu[1+jk3])*_Complex_I;
         bxy[jk3] = zero;
         bxy[1+jk3] = -at2*zt1;
         bxy[2+jk3] = at2*zt2;
         bxy[3*k1+jk3] = zero;
         bxy[1+3*k1+jk3] = zero;
         bxy[2+3*k1+jk3] = zero;
         wp += at1*(cu[jk3]*conjf(cu[jk3])
            + cu[1+jk3]*conjf(cu[1+jk3])
            + cu[2+jk3]*conjf(cu[2+jk3]));
      }
      sum1 += wp;
   }
   wp = 0.0;
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
   sum1 += wp;
L40:
   *wm = sum1*((float) nx)*((float) ny);
   return;
}

/*--------------------------------------------------------------------*/
void cmppmaxwel2(float complex exy[], float complex bxy[],
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
   ci = reciprical of velocity of light
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
   int nxh, nyh, ks, joff, kxps, j, jj, jk3, k, k1;
   float dnx, dny, dth, c2, cdt, adt, anorm, dkx, dky, afdt;
   float complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9;
   double wp, ws, sum1, sum2;
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
   sum1 = 0.0;
   sum2 = 0.0;
   if (kstrt > nxh)
      goto L40;
/* calculate the electromagnetic fields */
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,jj,jk3,dkx,dky,afdt,zt1,zt2,zt3,zt4,zt5,zt6,zt7,zt8,zt9, \
ws,wp) \
reduction(+:sum1,sum2)
   for (j = 0; j < kxps; j++) {
      dkx = dnx*(float) (j + joff);
      jj = nyhd*j;
      jk3 = 3*nyv*j;
      ws = 0.0;
      wp = 0.0;
      if ((j+joff) > 0) {
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            dky = dny*(float) k;
            afdt = adt*cimagf(ffc[k+jj]);
/* update magnetic field half time step, ky > 0 */
            zt1 = -cimagf(exy[2+3*k+jk3])
                + crealf(exy[2+3*k+jk3])*_Complex_I;
            zt2 = -cimagf(exy[1+3*k+jk3])
                + crealf(exy[1+3*k+jk3])*_Complex_I;
            zt3 = -cimagf(exy[3*k+jk3])
                + crealf(exy[3*k+jk3])*_Complex_I;
            zt4 = bxy[3*k+jk3] - dth*(dky*zt1);
            zt5 = bxy[1+3*k+jk3] + dth*(dkx*zt1);
            zt6 = bxy[2+3*k+jk3] - dth*(dkx*zt2 - dky*zt3);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exy[3*k+jk3] + cdt*(dky*zt1) - afdt*cu[3*k+jk3];
            zt8 = exy[1+3*k+jk3] - cdt*(dkx*zt1) - afdt*cu[1+3*k+jk3];
            zt9 = exy[2+3*k+jk3] + cdt*(dkx*zt2 - dky*zt3) 
                - afdt*cu[2+3*k+jk3];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exy[3*k+jk3] = zt7;
            exy[1+3*k+jk3] = zt8;
            exy[2+3*k+jk3] = zt9;
            ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) 
               + zt9*conjf(zt9));
            zt4 -= dth*(dky*zt1);
            zt5 += dth*(dkx*zt1);
            zt6 -= dth*(dkx*zt2 - dky*zt3);
            bxy[3*k+jk3] = zt4;
            bxy[1+3*k+jk3] = zt5;
            bxy[2+3*k+jk3] = zt6;
            wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
               + zt6*conjf(zt6));
/* update magnetic field half time step, ky < 0 */
            zt1 = -cimagf(exy[2+3*k1+jk3])
                + crealf(exy[2+3*k1+jk3])*_Complex_I;
            zt2 = -cimagf(exy[1+3*k1+jk3])
                + crealf(exy[1+3*k1+jk3])*_Complex_I;
            zt3 = -cimagf(exy[3*k1+jk3])
                + crealf(exy[3*k1+jk3])*_Complex_I;
            zt4 = bxy[3*k1+jk3] + dth*(dky*zt1);
            zt5 = bxy[1+3*k1+jk3] + dth*(dkx*zt1);
            zt6 = bxy[2+3*k1+jk3] - dth*(dkx*zt2 + dky*zt3);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exy[3*k1+jk3] - cdt*(dky*zt1) - afdt*cu[3*k1+jk3];
            zt8 = exy[1+3*k1+jk3] - cdt*(dkx*zt1)
                - afdt*cu[1+3*k1+jk3];
            zt9 = exy[2+3*k1+jk3] + cdt*(dkx*zt2 + dky*zt3)
                - afdt*cu[2+3*k1+jk3];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exy[3*k1+jk3] = zt7;
            exy[1+3*k1+jk3] = zt8;
            exy[2+3*k1+jk3] = zt9;
            ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) 
               + zt9*conjf(zt9));
            zt4 += dth*(dky*zt1);
            zt5 += dth*(dkx*zt1);
            zt6 -= dth*(dkx*zt2 + dky*zt3);
            bxy[3*k1+jk3] = zt4;
            bxy[1+3*k1+jk3] = zt5;
            bxy[2+3*k1+jk3] = zt6;
            wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
               + zt6*conjf(zt6));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         afdt = adt*cimagf(ffc[jj]);
/* update magnetic field half time step */
         zt1 = -cimagf(exy[2+jk3])
             + crealf(exy[2+jk3])*_Complex_I;
         zt2 = -cimagf(exy[1+jk3])
             + crealf(exy[1+jk3])*_Complex_I;
         zt5 = bxy[1+jk3] + dth*(dkx*zt1);
         zt6 = bxy[2+jk3] - dth*(dkx*zt2);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt8 = exy[1+jk3] - cdt*(dkx*zt1) - afdt*cu[1+jk3];
         zt9 = exy[2+jk3] + cdt*(dkx*zt2) - afdt*cu[2+jk3];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         exy[jk3] = zero;
         exy[1+jk3] = zt8;
         exy[2+jk3] = zt9;
         ws += anorm*(zt8*conjf(zt8) + zt9*conjf(zt9));
         zt5 = zt5 + dth*(dkx*zt1);
         zt6 = zt6 - dth*(dkx*zt2);
         bxy[jk3] = zero;
         bxy[1+jk3] = zt5;
         bxy[2+jk3] = zt6;
         wp += anorm*(zt5*conjf(zt5) + zt6*conjf(zt6));
         bxy[3*k1+jk3] = zero;
         bxy[1+3*k1+jk3] = zero;
         bxy[2+3*k1+jk3] = zero;
         exy[3*k1+jk3] = zero;
         exy[1+3*k1+jk3] = zero;
         exy[2+3*k1+jk3] = zero;
      }
      sum1 += ws;
      sum2 += wp;
   }
   ws = 0.0;
   wp = 0.0;
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
   sum1 += ws;
   sum2 += wp;
L40:
   *wf = sum1*((float) nx)*((float) ny);
   *wm = c2*sum2*((float) nx)*((float) ny);
   return;
}

/*--------------------------------------------------------------------*/
void cmppemfield2(float complex fxy[], float complex exy[],
                  float complex ffc[], int isign, int nx, int ny,
                  int kstrt, int nyv, int kxp, int nyhd) {
/* this subroutine either adds complex vector fields if isign > 0
   or copies complex vector fields if isign < 0
   includes additional smoothing
local data                                                 */
   int i, nxh, nyh, ks, joff, kxps, j, jj, jk3, k, k1;
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
#pragma omp parallel for private(i,j,k,k1,jj,jk3,at1)
      for (j = 0; j < kxps; j++) {
         jj = nyhd*j;
         jk3 = 3*nyv*j;
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            at1 = cimagf(ffc[k+jj]);
            for (i = 0; i < 3; i++) {
               fxy[i+3*k+jk3] += exy[i+3*k+jk3]*at1;
               fxy[i+3*k1+jk3] += exy[i+3*k1+jk3]*at1;
            }
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         at1 = cimagf(ffc[jj]);
         for (i = 0; i < 3; i++) {
            fxy[i+jk3] += exy[i+jk3]*at1;
            fxy[i+3*k1+jk3] += exy[i+3*k1+jk3]*at1;
         }
      }
   }
/* copy the fields */
   else if (isign < 0) {
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for private(i,j,k,k1,jj,jk3,at1)
      for (j = 0; j < kxps; j++) {
         jj = nyhd*j;
         jk3 = 3*nyv*j;
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            at1 = cimagf(ffc[k+jj]);
            for (i = 0; i < 3; i++) {
               fxy[i+3*k+jk3] = exy[i+3*k+jk3]*at1;
               fxy[i+3*k1+jk3] = exy[i+3*k1+jk3]*at1;
            }
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         at1 = cimagf(ffc[jj]);
         for (i = 0; i < 3; i++) {
            fxy[i+jk3] = exy[i+jk3]*at1;
            fxy[i+3*k1+jk3] = exy[i+3*k1+jk3]*at1;
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
void cppfft2rmxx(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                 int nxyhd) {
/* this subroutine performs the x part of a two dimensional real to
   complex fast fourier transform and its inverse, for a subset of y,
   using complex arithmetic, with OpenMP,
   for data which is distributed in blocks
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
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nrxb, joff;
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
      goto L70;
/* inverse fourier transform */
   ani = 0.5/(((float) nx)*((float) ny));
   nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,joff,s,t,t1)
   for (i = kypi-1; i < kypt; i++) {
      joff = nxvh*i;
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            t = f[j1+joff];
            f[j1+joff] = f[j+joff];
            f[j+joff] = t;
         }
      }
/* then transform in x */
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
               t = s*f[j2+joff];
               f[j2+joff] = f[j1+joff] - t;
               f[j1+joff] += t;
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
         t1 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         t = conjf(f[nxh-j+joff]);
         s = f[j+joff] + t;
         t = (f[j+joff] - t)*t1;
         f[j+joff] = ani*(s + t);
         f[nxh-j+joff] = ani*conjf(s - t);
      }
      f[joff] = 2.0*ani*((crealf(f[joff]) + cimagf(f[joff]))
                + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I);
      if (nxhh > 0)
         f[nxhh+joff] = 2.0*ani*conjf(f[nxhh+joff]);
   }
   return;
/* forward fourier transform */
L70: nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,joff,s,t,t1)
   for (i = kypi-1; i < kypt; i++) {
      joff = nxvh*i;
/* scramble coefficients */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
         t1 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         t = conjf(f[nxh-j+joff]);
         s = f[j+joff] + t;
         t = (f[j+joff] - t)*t1;
         f[j+joff] = s + t;
         f[nxh-j+joff] = conjf(s - t);
      }
      f[joff] = (crealf(f[joff]) + cimagf(f[joff]))
                + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I;
      if (nxhh > 0)
         f[nxhh+joff] = 2.0*conjf(f[nxhh+joff]);
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            t = f[j1+joff];
            f[j1+joff] = f[j+joff];
            f[j+joff] = t;
         }
      }
/* then transform in x */
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
               t = s*f[j2+joff];
               f[j2+joff] = f[j1+joff] - t;
               f[j1+joff] += t;
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rmxy(float complex g[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                 int nxyhd) {
/* this subroutine performs the y part of a two dimensional real to
   complex fast fourier transform and its inverse, for a subset of x,
   using complex arithmetic, with OpenMP,
   for data which is distributed in blocks
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
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nryb, koff;
   float complex s, t;
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   nxh = nx/2;
   ny = 1L<<indy;
   nyh = ny/2;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   ks = kstrt - 1;
   kxpt = kxpi + kxpp - 1;
   if (kstrt > nxh)
      return;
   if (isign > 0)
      goto L70;
/* inverse fourier transform */
   nryb = nxhy/ny;
   nry = nxy/ny;
#pragma omp parallel for \
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,koff,s,t)
   for (i = kxpi-1; i < kxpt; i++) {
      koff = nyv*i;
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            t = g[k1+koff];
            g[k1+koff] = g[k+koff];
            g[k+koff] = t;
         }
      }
/* then transform in y */
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
               t = s*g[j2+koff];
               g[j2+koff] = g[j1+koff] - t;
               g[j1+koff] += t;
            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if ((ks==0) && (kxpi==1)) {
      for (k = 1; k < nyh; k++) {
         s = g[ny-k];
         g[ny-k] = 0.5*(cimagf(g[k] + s) + crealf(g[k] - s)*_Complex_I);
         g[k] = 0.5*(crealf(g[k] + s) + cimagf(g[k] - s)*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
L70: nryb = nxhy/ny;
   nry = nxy/ny;
/* scramble modes kx = 0, nx/2 */
   if ((ks==0) && (kxpi==1)) {
      for (k = 1; k < nyh; k++) {
         s = cimagf(g[ny-k]) + crealf(g[ny-k])*_Complex_I;
         g[ny-k] = conjf(g[k] - s);
         g[k] += s;
      }
   }
#pragma omp parallel for \
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,koff,s,t)
   for (i = kxpi-1; i < kxpt; i++) {
      koff = nyv*i;
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            t = g[k1+koff];
            g[k1+koff] = g[k+koff];
            g[k+koff] = t;
         }
      }
/* then transform in y */
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
               t = s*g[j2+koff];
               g[j2+koff] = g[j1+koff] - t;
               g[j1+koff] += t;
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rm3xx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int kstrt,
                  int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                  int nxyhd) {
/* this subroutine performs the x part of 3 two dimensional real to
   complex fast fourier transforms and their inverses, for a subset of y,
   using complex arithmetic, with OpenMP,
   for data which is distributed in blocks
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
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nrxb, joff;
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
      goto L100;
/* inverse fourier transform */
   ani = 0.5/(((float) nx)*((float) ny));
   nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,joff,at1,at2,s,t,t1,t2,t3)
   for (i = kypi-1; i < kypt; i++) {
      joff = 3*nxvh*i;
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
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
         t1 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < 3; k++) {
            t = conjf(f[k+3*(nxh-j)+joff]);
            s = f[k+3*j+joff] + t;
            t = (f[k+3*j+joff] - t)*t1;
            f[k+3*j+joff] = ani*(s + t);
            f[k+3*(nxh-j)+joff] = ani*conjf(s - t);
         }
      }
      for (k = 0; k < 3; k++) {
         f[k+joff] = 2.0*ani*((crealf(f[k+joff]) + cimagf(f[k+joff]))
                     + (crealf(f[k+joff]) - cimagf(f[k+joff]))*_Complex_I);
         if (nxhh > 0)
            f[k+3*nxhh+joff] = 2.0*ani*conjf(f[k+3*nxhh+joff]);
      }
   }
   return;
/* forward fourier transform */
L100: nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,joff,at1,at2,s,t,t1,t2,t3)
   for (i = kypi-1; i < kypt; i++) {
      joff = 3*nxvh*i;
/* scramble coefficients */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
         t1 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < 3; k++) {
            t = conjf(f[k+3*(nxh-j)+joff]);
            s = f[k+3*j+joff] + t;
            t = (f[k+3*j+joff] - t)*t1;
            f[k+3*j+joff] = s + t;
            f[k+3*(nxh-j)+joff] = conjf(s - t);
         }
      }
      for (k = 0; k < 3; k++) {
         f[k+joff] = (crealf(f[k+joff]) + cimagf(f[k+joff]))
                      + (crealf(f[k+joff]) - cimagf(f[k+joff]))*_Complex_I;
         if (nxhh > 0)
            f[k+3*nxhh+joff] = 2.0*conjf(f[k+3*nxhh+joff]);
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
void cppfft2rm3xy(float complex g[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int kstrt,
                  int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                  int nxyhd) {
/* this subroutine performs the y part of 3 two dimensional real to
   complex fast fourier transforms and their inverses, for a subset of x,
   using complex arithmetic, with OpenMP,
   for data which is distributed in blocks
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
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nryb, koff;
   float complex s, t1, t2, t3;
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   nxh = nx/2;
   ny = 1L<<indy;
   nyh = ny/2;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   ks = kstrt - 1;
   kxpt = kxpi + kxpp - 1;
   if (kstrt > nxh)
      return;
   if (isign > 0)
      goto L80;
/* inverse fourier transform */
   nryb = nxhy/ny;
   nry = nxy/ny;
#pragma omp parallel for \
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,koff,s,t1,t2,t3)
   for (i = kxpi-1; i < kxpt; i++) {
      koff = 3*nyv*i;
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
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
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if ((ks==0) && (kxpi==1)) {
      for (k = 1; k < nyh; k++) {
         for (j = 0; j < 3; j++) {
            s = g[j+3*(ny-k)];
            g[j+3*(ny-k)] = 0.5*(cimagf(g[j+3*k] + s)
                             + crealf(g[j+3*k] - s)*_Complex_I);
            g[j+3*k] = 0.5*(crealf(g[j+3*k] + s)
                        + cimagf(g[j+3*k] - s)*_Complex_I);
         }
      }
   }
   return;
/* forward fourier transform */
L80: nryb = nxhy/ny;
   nry = nxy/ny;
/* scramble modes kx = 0, nx/2 */
   if ((ks==0) && (kxpi==1)) {
      for (k = 1; k < nyh; k++) {
         for (j = 0; j < 3; j++) {
             s = cimagf(g[j+3*(ny-k)])
                 + crealf(g[j+3*(ny-k)])*_Complex_I;
             g[j+3*(ny-k)] = conjf(g[j+3*k] - s);
             g[j+3*k] += s;
         }
      }
   }
#pragma omp parallel for \
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,koff,s,t1,t2,t3)
   for (i = kxpi-1; i < kxpt; i++) {
      koff = 3*nyv*i;
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
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
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2rm(float complex f[], float complex g[],
                float complex bs[], float complex br[], int isign,
                int ntpose, int mixup[], float complex sct[],
                float *ttp, int indx, int indy, int kstrt, int nvp,
                int nxvh, int nyv, int kxp, int kyp, int kypd,
                int nxhyd, int nxyhd) {
/* wrapper function for parallel real to complex fft */
/* parallelized with OpenMP */
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
      cppfft2rmxx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,kypd,
                 nxhyd,nxyhd);
/* transpose f array to g */
      cpwtimera(-1,ttp,&dtime);
      cpptpose(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,nxvh,nyv,kxp,kypd);
      cpwtimera(1,ttp,&dtime);
/* perform y fft */
      cppfft2rmxy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
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
      cppfft2rmxy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
                 nxhyd,nxyhd);
/* transpose g array to f */
      cpwtimera(-1,ttp,&dtime);
      cpptpose(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,nyv,nxvh,kypd,kxp);
      cpwtimera(1,ttp,&dtime);
/* perform x fft */
      cppfft2rmxx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,kypd,
                 nxhyd,nxyhd);
   }
   if (ntpose==0)
      *ttp += tf;
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2rm3(float complex f[], float complex g[],
                 float complex bs[], float complex br[], int isign,
                 int ntpose, int mixup[], float complex sct[],
                 float *ttp, int indx, int indy, int kstrt, int nvp,
                 int nxvh, int nyv, int kxp, int kyp, int kypd,
                 int nxhyd, int nxyhd) {
/* wrapper function for parallel real to complex fft */
/* parallelized with OpenMP */
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
      cppfft2rm3xx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,
                   kypd,nxhyd,nxyhd);
/* transpose f array to g */
      cpwtimera(-1,ttp,&dtime);
      cppntpose(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,3,nxvh,nyv,kxp,kypd);
      cpwtimera(1,ttp,&dtime);
/* perform y fft */
      cppfft2rm3xy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
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
      cppfft2rm3xy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
                   nxhyd,nxyhd);
/* transpose g array to f */
      cpwtimera(-1,ttp,&dtime);
      cppntpose(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,3,nyv,nxvh,kypd,kxp);
      cpwtimera(1,ttp,&dtime);
/* perform x fft */
      cppfft2rm3xx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,
                   kypd,nxhyd,nxyhd);
   }
   if (ntpose==0)
      *ttp += tf;
   return;
}

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
void cppdblkp2l_(float *part, int *kpic, int *npp, int *noff, 
                 int *nppmx, int *idimp, int *npmax, int *mx, int *my,
                 int *mx1,int *mxyp1, int *irc) {
   cppdblkp2l(part,kpic,*npp,*noff,nppmx,*idimp,*npmax,*mx,*my,*mx1,
              *mxyp1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpppmovin2l_(float *part, float *ppart, int *kpic, int *npp,
                  int *noff, int *nppmx, int *idimp, int *npmax,
                  int *mx, int *my, int *mx1, int *mxyp1, int *irc) {
   cpppmovin2l(part,ppart,kpic,*npp,*noff,*nppmx,*idimp,*npmax,*mx,*my,
               *mx1,*mxyp1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpppcheck2l_(float *ppart, int *kpic, int *noff, int *nyp,
                  int *idimp, int *nppmx, int *nx, int *mx, int *my,
                  int *mx1, int *myp1, int *irc) {
   cpppcheck2l(ppart,kpic,*noff,*nyp,*idimp,*nppmx,*nx,*mx,*my,*mx1,
               *myp1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgbppush23l_(float *ppart, float *fxy, float *bxy, int *kpic,
                    int *noff, int *nyp, float *qbm, float *dt,
                    float *dtc, float *ek, int *idimp, int *nppmx, 
                    int *nx, int *ny, int *mx, int *my, int *nxv,
                    int *nypmx, int *mx1, int *mxyp1, int *ipbc) {
   cppgbppush23l(ppart,fxy,bxy,kpic,*noff,*nyp,*qbm,*dt,*dtc,ek,*idimp,
                 *nppmx,*nx,*ny,*mx,*my,*nxv,*nypmx,*mx1,*mxyp1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgbppushf23l_(float *ppart, float *fxy, float *bxy, int *kpic,
                     int *ncl, int *ihole, int *noff, int *nyp,
                     float *qbm, float *dt, float *dtc, float *ek,
                     int *idimp, int *nppmx, int *nx, int *ny,
                     int *mx, int *my, int *nxv, int *nypmx, int *mx1,
                     int *mxyp1, int *ntmax, int *irc) {
   cppgbppushf23l(ppart,fxy,bxy,kpic,ncl,ihole,*noff,*nyp,*qbm,*dt,*dtc,
                  ek,*idimp,*nppmx,*nx,*ny,*mx,*my,*nxv,*nypmx,*mx1,
                  *mxyp1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrbppush23l_(float *ppart, float *fxy, float *bxy, int *kpic,
                     int *noff, int *nyp, float *qbm, float *dt, float *dtc, float *ci,
                     float *ek, int *idimp, int *nppmx, int *nx, int *ny,
                     int *mx, int *my, int *nxv, int *nypmx, int *mx1,
                     int *mxyp1, int *ipbc) {
   cppgrbppush23l(ppart,fxy,bxy,kpic,*noff,*nyp,*qbm,*dt,*dtc,*ci,ek,
                  *idimp,*nppmx,*nx,*ny,*mx,*my,*nxv,*nypmx,*mx1,*mxyp1,
                  *ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrbppushf23l_(float *ppart, float *fxy, float *bxy, int *kpic,
                      int *ncl, int *ihole, int *noff, int *nyp,
                      float *qbm, float *dt, float *dtc, float *ci,
                      float *ek, int *idimp, int *nppmx, int *nx,
                      int *ny, int *mx, int *my, int *nxv, int *nypmx,
                      int *mx1, int *mxyp1, int *ntmax, int *irc) {
   cppgrbppushf23l(ppart,fxy,bxy,kpic,ncl,ihole,*noff,*nyp,*qbm,*dt,
                   *dtc,*ci,ek,*idimp,*nppmx,*nx,*ny,*mx,*my,*nxv,
                   *nypmx,*mx1,*mxyp1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgppost2l_(float *ppart, float *q, int *kpic, int *noff,
                  float *qm, int *idimp, int *nppmx, int *mx, int *my,
                  int *nxv, int *nypmx, int *mx1, int *mxyp1) {
   cppgppost2l(ppart,q,kpic,*noff, *qm,*idimp,*nppmx,*mx,*my,*nxv,
               *nypmx,*mx1,*mxyp1);
   return;
}

/*--------------------------------------------------------------------*/
void cppgjppost2l_(float *ppart, float *cu, int *kpic, int *noff,
                   float *qm, float *dt, int *nppmx, int *idimp,
                   int *nx, int *ny, int *mx, int *my, int *nxv,
                   int *nypmx, int *mx1, int *mxyp1, int *ipbc) {
   cppgjppost2l(ppart,cu,kpic,*noff,*qm,*dt,*nppmx,*idimp,*nx,*ny,*mx,
                *my,*nxv,*nypmx,*mx1,*mxyp1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgjppostf2l_(float *ppart, float *cu, int *kpic, int *ncl,
                    int *ihole, int *noff, int *nyp, float *qm,
                    float *dt, int *nppmx, int *idimp, int *nx, int *ny,
                    int *mx, int *my, int *nxv, int *nypmx, int *mx1,
                    int *mxyp1, int *ntmax, int *irc) {
   cppgjppostf2l(ppart,cu,kpic,ncl,ihole,*noff,*nyp,*qm,*dt,*nppmx,
                 *idimp,*nx,*ny,*mx,*my,*nxv,*nypmx,*mx1,*mxyp1,*ntmax,
                 irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrjppost2l_(float *ppart, float *cu, int *kpic, int *noff,
                    float *qm, float *dt, float *ci, int *nppmx,
                    int *idimp, int *nx, int *ny, int *mx, int *my,
                    int *nxv, int *nypmx, int *mx1, int *mxyp1,
                    int *ipbc) {
   cppgrjppost2l(ppart,cu,kpic,*noff,*qm,*dt,*ci,*nppmx,*idimp,*nx,*ny,
                 *mx,*my,*nxv,*nypmx,*mx1,*mxyp1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrjppostf2l_(float *ppart, float *cu, int *kpic, int *ncl,
                     int *ihole, int *noff, int *nyp, float *qm,
                     float *dt, float *ci, int *nppmx, int *idimp,
                     int *nx, int *ny, int *mx, int *my, int *nxv, 
                     int *nypmx, int *mx1, int *mxyp1, int *ntmax, 
                     int *irc) {
   cppgrjppostf2l(ppart,cu,kpic,ncl,ihole,*noff,*nyp,*qm,*dt,*ci,*nppmx,
                  *idimp,*nx,*ny,*mx,*my,*nxv,*nypmx,*mx1,*mxyp1,*ntmax,
                  irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppporder2la_(float *ppart, float *ppbuff, float *sbufl,
                   float *sbufr, int *kpic, int *ncl, int *ihole,
                   int *ncll, int *nclr, int *noff, int *nyp,
                   int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                   int *my, int *mx1, int *myp1, int *npbmx, int *ntmax,
                   int *nbmax, int *irc) {
   cppporder2la(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,nclr,*noff,
                *nyp,*idimp,*nppmx,*nx,*ny,*mx,*my,*mx1,*myp1,*npbmx,
                *ntmax,*nbmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppporderf2la_(float *ppart, float *ppbuff, float *sbufl,
                    float *sbufr, int *ncl, int *ihole, int *ncll,
                    int *nclr, int *idimp, int *nppmx, int *mx1,
                    int *myp1, int *npbmx, int *ntmax, int *nbmax,
                    int *irc) {
   cppporderf2la(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr,*idimp,
                 *nppmx,*mx1,*myp1,*npbmx,*ntmax,*nbmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppporder2lb_(float *ppart, float *ppbuff, float *rbufl,
                   float *rbufr, int *kpic, int *ncl, int *ihole,
                   int *mcll, int *mclr, int *idimp, int *nppmx,
                   int *mx1, int *myp1, int *npbmx, int *ntmax,
                   int *nbmax, int *irc) {
   cppporder2lb(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,mclr,
                *idimp,*nppmx,*mx1,*myp1,*npbmx,*ntmax,*nbmax,irc);
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
void cmppois23_(float complex *q, float complex *fxy, int *isign,
                float complex *ffc, float *ax, float *ay, float *affp,
                float *we, int *nx, int *ny, int *kstrt, int *nyv,
                int *kxp, int *nyhd) {
   cmppois23(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*kstrt,*nyv,*kxp,
             *nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppcuperp2_(float complex *cu, int *nx, int *ny, int *kstrt,
                  int *nyv, int *kxp) {
   cmppcuperp2(cu,*nx,*ny,*kstrt,*nyv,*kxp);
   return;
}

/*--------------------------------------------------------------------*/
void cmippbpoisp23_(float complex *cu, float complex *bxy,
                    float complex *ffc, float *ci, float *wm, int *nx,
                    int *ny, int *kstrt, int *nyv, int *kxp,
                    int *nyhd) {
   cmippbpoisp23(cu,bxy,ffc,*ci,wm,*nx,*ny,*kstrt,*nyv,*kxp,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppmaxwel2_(float complex *exy, float complex *bxy,
                  float complex *cu, float complex *ffc, float *affp,
                  float *ci, float *dt, float *wf, float *wm, int *nx,
                  int *ny, int *kstrt, int *nyv, int *kxp, int *nyhd) {
   cmppmaxwel2(exy,bxy,cu,ffc,*affp,*ci,*dt,wf,wm,*nx,*ny,*kstrt,*nyv,
               *kxp,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppemfield2_(float complex *fxy, float complex *exy,
                   float complex *ffc, int *isign, int *nx, int *ny,
                   int *kstrt, int *nyv, int *kxp, int *nyhd) {
   cmppemfield2(fxy,exy,ffc,*isign,*nx,*ny,*kstrt,*nyv,*kxp,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwpfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                   int *nxhyd, int *nxyhd) {
   cwpfft2rinit(mixup,sct,*indx,*indy,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rmxx_(float complex *f, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *kstrt,
                  int *kypi, int *kypp, int *nxvh, int *kypd,
                  int *nxhyd, int *nxyhd) {
   cppfft2rmxx(f,*isign,mixup,sct,*indx,*indy,*kstrt,*kypi,*kypp,*nxvh,
                *kypd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rmxy_(float complex *g, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *kstrt,
                  int *kxpi, int *kxpp, int *nyv, int *kxp, int *nxhyd,
                  int *nxyhd) {
   cppfft2rmxy(g,*isign,mixup,sct,*indx,*indy,*kstrt,*kxpi,*kxpp,*nyv,
               *kxp,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rm3xx_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *kstrt,
                   int *kypi, int *kypp, int *nxvh, int *kypd,
                   int *nxhyd, int *nxyhd) {
   cppfft2rm3xx(f,*isign,mixup,sct,*indx,*indy,*kstrt,*kypi,*kypp,*nxvh,
                *kypd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rm3xy_(float complex *g, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *kstrt,
                   int *kxpi, int *kxpp, int *nyv, int *kxp, int *nxhyd,
                   int *nxyhd) {
   cppfft2rm3xy(g,*isign,mixup,sct,*indx,*indy,*kstrt,*kxpi,*kxpp,*nyv,
                *kxp,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2rm_(float complex *f, float complex *g, float complex *bs,
                 float complex *br, int *isign, int *ntpose, int *mixup,
                 float complex *sct, float *ttp, int *indx, int *indy,
                 int *kstrt, int *nvp, int *nxvh, int *nyv, int *kxp,
                 int *kyp, int *kypd, int *nxhyd, int *nxyhd) {
   cwppfft2rm(f,g,bs,br,*isign,*ntpose,mixup,sct,ttp,*indx,*indy,*kstrt,
              *nvp,*nxvh,*nyv,*kxp,*kyp,*kypd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2rm3_(float complex *f, float complex *g, float complex *bs,
                  float complex *br, int *isign, int *ntpose,
                  int *mixup, float complex *sct, float *ttp, int *indx,
                  int *indy, int *kstrt, int *nvp, int *nxvh, int *nyv,
                  int *kxp, int *kyp, int *kypd, int *nxhyd,
                  int *nxyhd) {
   cwppfft2rm3(f,g,bs,br,*isign,*ntpose,mixup,sct,ttp,*indx,*indy,
               *kstrt,*nvp,*nxvh,*nyv,*kxp,*kyp,*kypd,*nxhyd,*nxyhd);
   return;
}
