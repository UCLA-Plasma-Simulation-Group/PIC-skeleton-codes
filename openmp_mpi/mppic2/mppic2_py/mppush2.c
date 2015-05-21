/* C Library for Skeleton 2D Electrostatic MPI/OpenMP PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "mppush2.h"
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
void cpdistr2(float part[], float edges[], int *npp, int nps, float vtx,
              float vty, float vdx, float vdy, int npx, int npy, int nx,
              int ny, int idimp, int npmax, int idps, int ipbc, int *ierr) {
/* for 2d code, this subroutine calculates initial particle co-ordinates
   and velocities with uniform density and maxwellian velocity with drift
   for distributed data.
   input: all except part, npp, ierr, output: part, npp, ierr
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   part[n][2] = velocity vx of particle n in partition
   part[n][3] = velocity vy of particle n in partition
   edges[0] = lower boundary of particle partition
   edges[1] = upper boundary of particle partition
   npp = number of particles in partition
   nps = starting address of particles in partition
   vtx/vty = thermal velocity of electrons in x/y direction
   vdx/vdy = drift velocity of beam electrons in x/y direction
   npx/npy = initial number of particles distributed in x/y direction
   nx/ny = system length in x/y direction
   idimp = size of phase space = 4
   npmax = maximum number of particles in each partition
   idps = number of partition boundaries
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   ierr = (0,1) = (no,yes) error condition exists
   ranorm = gaussian random number with zero mean and unit variance
   with spatial decomposition
local data                                                            */
   int j, k, npt, k1, npxyp;
   float edgelx, edgely, at1, at2, xt, yt, vxt, vyt;
   double dnpx, dnpxy, dt1;
   int ierr1[1], iwork1[1];
   double sum3[3], work3[3];
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
         if ((yt >= edges[0]) && (yt < edges[1])) {
            if (npt < npmax) {
               k1 = idimp*npt;
               part[k1] = xt;
               part[1+k1] = yt;
               part[2+k1] = vxt;
               part[3+k1] = vyt;
               npt += 1;
            }
            else
               *ierr += 1;
         }
      }
   }
   npxyp = 0;
/* add correct drift */
   sum3[0] = 0.0;
   sum3[1] = 0.0;
   for (j = nps-1; j < npt; j++) {
      npxyp += 1;
      sum3[0] += part[2+idimp*j];
      sum3[1] += part[3+idimp*j];
   }
   sum3[2] = npxyp;
   cppdsum(sum3,work3,3);
   dnpxy = sum3[2];
   ierr1[0] = *ierr;
   cppimax(ierr1,iwork1,1);
   *ierr = ierr1[0];
   dt1 = 1.0/dnpxy;
   sum3[0] = dt1*sum3[0] - vdx;
   sum3[1] = dt1*sum3[1] - vdy;
   for (j = nps-1; j < npt; j++) {
      part[2+idimp*j] -= sum3[0];
      part[3+idimp*j] -= sum3[1];
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
void cppgppush2l(float ppart[], float fxy[], int kpic[], int noff,
                 int nyp, float qbm, float dt, float *ek, int nx,
                 int ny, int mx, int my, int idimp, int nppmx, int nxv,
                 int nypmx, int mx1, int mxyp1, int ipbc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions
   OpenMP version using guard cells, for distributed data
   data read in tiles
   particles stored segmented array
   42 flops/particle, 12 loads, 4 stores
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
   ppart[m][n][0] = position x of particle n in partition in tile m
   ppart[m][n][1] = position y of particle n in partition in tile m
   ppart[m][n][2] = velocity vx of particle n in partition in tile m
   ppart[m][n][3] = velocity vy of particle n in partition in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,kk)
   fxy[k][j][1] = y component of force/charge at grid (j,kk)
   in other words, fxy are the convolutions of the electric field
   over the particle shape, where kk = k + noff
   kpic = number of particles per tile
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nxv = first dimension of field array, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
#define MXV             33
#define MYV             33
   int noffp, moffp, npoff, nppp;
   int mnoff, i, j, k, nn, mm, mxv;
   float qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
   float sfxy[2*MXV*MYV];
/* float sfxy[2*(mx+1)*(my+1)]; */
   double sum1, sum2;
/* mxv2 = 2*MXV; */
   mxv = mx + 1;
   qtm = qbm*dt;
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
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noffp,moffp,nppp,npoff,mnoff,nn,mm,x,y,dxp,dyp,amx,amy, \
dx,dy,vx,vy,sum1,sfxy) \
reduction(+:sum2)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      npoff = nppmx*k;
      mnoff = moffp + noff;
/* load local fields from global array */
      nn = (mx < nx-noffp ? mx : nx-noffp) + 1;
      mm = (my < nyp-moffp ? my : nyp-moffp) + 1;
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[2*(i+mxv*j)] = fxy[2*(i+noffp+nxv*(j+moffp))];
            sfxy[1+2*(i+mxv*j)] = fxy[1+2*(i+noffp+nxv*(j+moffp))];
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
         nn = 2*(nn - noffp) + 2*mxv*(mm - mnoff);
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
         vx = ppart[2+idimp*(j+npoff)];
         vy = ppart[3+idimp*(j+npoff)];
         dx = vx + qtm*dx;
         dy = vy + qtm*dy;
/* average kinetic energy */
         vx += dx;
         vy += dy;
         sum1 += vx*vx + vy*vy;
         ppart[2+idimp*(j+npoff)] = dx;
         ppart[3+idimp*(j+npoff)] = dy;
/* new position */
         dx = x + dx*dt;
         dy = y + dy*dt;
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
   *ek += 0.125f*sum2;
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void cppgppushf2l(float ppart[], float fxy[], int kpic[], int ncl[],
                  int ihole[], int noff, int nyp, float qbm, float dt,
                  float *ek, int nx, int ny, int mx, int my, int idimp,
                  int nppmx, int nxv, int nypmx, int mx1, int mxyp1,
                  int ntmax, int *irc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with periodic boundary conditions
   also determines list of particles which are leaving this tile
   OpenMP version using guard cells, for distributed data
   data read in tiles
   particles stored segmented array
   42 flops/particle, 12 loads, 4 stores
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
   ppart[m][n][0] = position x of particle n in partition in tile m
   ppart[m][n][1] = position y of particle n in partition in tile m
   ppart[m][n][2] = velocity vx of particle n in partition in tile m
   ppart[m][n][3] = velocity vy of particle n in partition in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,kk)
   fxy[k][j][1] = y component of force/charge at grid (j,kk)
   in other words, fxy are the convolutions of the electric field
   over the particle shape, where kk = k + noff
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nxv = first dimension of field array, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells.
   mx1 = (system length in x direction - 1)/mx + 1
   mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
#define MXV             33
#define MYV             33
   int noffp, moffp, npoff, nppp;
   int mnoff, i, j, k, ih, nh, nn, mm, mxv;
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
private(i,j,k,noffp,moffp,nppp,npoff,nn,mm,ih,nh,mnoff,x,y,dxp,dyp, \
amx,amy,dx,dy,vx,vy,edgelx,edgely,edgerx,edgery,sum1,sfxy) \
reduction(+:sum2)
   for (k = 0; k < mxyp1; k++) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      npoff = nppmx*k;
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
/* load local fields from global array */
      for (j = 0; j < mm; j++) {
         for (i = 0; i < nn; i++) {
            sfxy[2*(i+mxv*j)] = fxy[2*(i+noffp+nxv*(j+moffp))];
            sfxy[1+2*(i+mxv*j)] = fxy[1+2*(i+noffp+nxv*(j+moffp))];
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
         nn = 2*(nn - noffp) + 2*mxv*(mm - mnoff);
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
         vx = ppart[2+idimp*(j+npoff)];
         vy = ppart[3+idimp*(j+npoff)];
         dx = vx + qtm*dx;
         dy = vy + qtm*dy;
/* average kinetic energy */
         vx += dx;
         vy += dy;
         sum1 += vx*vx + vy*vy;
         ppart[2+idimp*(j+npoff)] = dx;
         ppart[3+idimp*(j+npoff)] = dy;
/* new position */
         dx = x + dx*dt;
         dy = y + dy*dt;
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
   *ek += 0.125f*sum2;
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
void cppcguard2xl(float fxy[], int nyp, int nx, int ndim, int nxe,
                  int nypmx) {
/* replicate extended periodic vector field in x direction
   linear interpolation, for distributed data
   nyp = number of primary (complete) gridpoints in particle partition
   nx = system length in x direction
   ndim = leading dimension of array fxy
   nxe = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells
local data                                                 */
   int i, k, kk, myp1;
/* replicate edges of extended field */
   myp1 = nyp + 1;
   for (k = 0; k < myp1; k++) {
      kk = ndim*nxe*k;
      for (i = 0; i < ndim; i++) {
         fxy[i+ndim*nx+kk] = fxy[i+kk];
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppaguard2xl(float q[], int nyp, int nx, int nxe, int nypmx) {
/* accumulate extended periodic scalar field in x direction
   linear interpolation, for distributed data
   nyp = number of primary (complete) gridpoints in particle partition
   nx = system length in x direction
   nxe = first dimension of field arrays, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells
local data                                                 */
   int k, myp1;
/* accumulate edges of extended field */
   myp1 = nyp + 1;
   for (k = 0; k < myp1; k++) {
      q[nxe*k] += q[nx+nxe*k];
      q[nx+nxe*k] = 0.0;
   }
   return;
}
      
/*--------------------------------------------------------------------*/
void cmppois22(float complex q[], float complex fxy[], int isign,
               float complex ffc[], float ax, float ay, float affp,
               float *we, int nx, int ny, int kstrt, int nyv, int kxp,
               int nyhd) {
/* this subroutine solves 2d poisson's equation in fourier space for
   force/charge (or convolution of electric field over particle shape)
   with periodic boundary conditions, for distributed data.
   for isign = 0, input: isign,ax,ay,affp,nx,ny,jblok,nyv,kxp,nyhd,
   output: ffc
   for isign /= 0, input: q,ffc,isign,nx,ny,nyv,kxp,jblok,nyhd,
   output: fxy,we
   approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
   where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
   the equation used is:
   fx[ky][kx] = -sqrt(-1)*kx*g[ky][kx]*s[ky][kx]*q[ky][kx],
   fy[ky][kx] = -sqrt(-1)*ky*g[ky][kx]*s[ky][kx]*q[ky][kx],
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   g[ky][kx] = (affp/(kx**2+ky**2))*s[ky][kx],
   s[ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
   fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
   fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
   q[k][j] = complex charge density for fourier mode (jj,k)
   fxy[k][j][0] = x component of complex force/charge,
   fxy[k][j][1] = y component of complex force/charge,
   for fourier mode (jj,k), where jj = j + kxp*(kstrt - 1)
   kxp = number of data values per block
   kstrt = starting data block number
   if isign = 0, form factor array is prepared
   if isign is not equal to 0, force/charge is calculated.
   aimag(ffc[k][j]) = finite-size particle shape factor s
   real(ffc[k][j])) = potential green's function g
   for fourier mode (jj,k), where jj = j + kxp*(kstrt - 1)
   ax/ay = half-width of particle in x/y direction
   affp = normalization constant = nx*ny/np, where np=number of particles
   electric field energy is also calculated, using
   we = nx*ny*sum((affp/(kx**2+ky**2))*|q[ky][kx]*s[ky][kx]|**2)
   nx/ny = system length in x/y direction
   nyv = first dimension of field arrays, must be >= ny
   nyhd = first dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, ks, joff, kxps, j, jj, jk, jk2, k, k1;
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
private(j,k,k1,jj,jk,jk2,dkx,at1,at2,at3,zt1,zt2,wp) \
reduction(+:sum1)
   for (j = 0; j < kxps; j++) {
      dkx = dnx*(float) (j + joff);
      jj = nyhd*j;
      jk = nyv*j;
      jk2 = 2*jk;
      wp = 0.0;
      if ((j+joff) > 0) {
         for (k = 1; k < nyh; k++) {
            k1 = ny - k;
            at1 = crealf(ffc[k+jj])*cimagf(ffc[k+jj]);
            at2 = dkx*at1;
            at3 = dny*at1*(float) k;
            zt1 = cimagf(q[k+jk]) - crealf(q[k+jk])*_Complex_I;
            zt2 = cimagf(q[k1+jk]) - crealf(q[k1+jk])*_Complex_I;
            fxy[2*k+jk2] = at2*zt1;
            fxy[1+2*k+jk2] = at3*zt1;
            fxy[2*k1+jk2] = at2*zt2;
            fxy[1+2*k1+jk2] = -at3*zt2;
            wp += at1*(q[k+jk]*conjf(q[k+jk])
                  + q[k1+jk]*conjf(q[k1+jk]));
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nyh;
         at1 = crealf(ffc[jj])*cimagf(ffc[jj]);
         at3 = dkx*at1;
         zt1 = cimagf(q[jk]) - crealf(q[jk])*_Complex_I;
         fxy[jk2] = at3*zt1;
         fxy[1+jk2] = zero;
         fxy[2*k1+jk2] = zero;
         fxy[1+2*k1+jk2] = zero;
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
         fxy[2*k] = zero;
         fxy[1+2*k] = at2*zt1;
         fxy[2*k1] = zero;
         fxy[1+2*k1] = zero;
         wp += at1*(q[k]*conjf(q[k]));
      }
      k1 = 2*nyh;
      fxy[0] = zero;
      fxy[1] = zero;
      fxy[k1] = zero;
      fxy[1+k1] = zero;
   }
   sum1 += wp;
L70:
   *we = sum1*((float) nx)*((float) ny);
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
void cppfft2rm2xx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int kstrt,
                  int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                  int nxyhd) {
/* this subroutine performs the x part of 2 two dimensional real to
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
   f[m][n][0:1] = (1/nx*ny)*sum(f[k][j][0:1]*exp(-sqrt(-1)*2pi*n*j/nx)
   if isign = 1, a forward fourier transform is performed
   f[k][j][0:1] = sum(f[m][n][0:1]*exp(sqrt(-1)*2pi*n*j/nx)*
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
   f[k][j][0:1] = mode j,kk, where kk = k + kyp*(kstrt - 1)
   0 <= j < nx/2 and 0 <= kk < ny, except for
   f[k][0][0:1] = mode nx/2,kk, where ny/2+1 <= kk < ny, and
   imaginary part of f[0][0][0:1] = real part of mode nx/2,0
   on mode kstrt=0
   imaginary part of f[0][0][0:1] = real part of mode nx/2,ny/2
   on mode kstrt=(ny/2)/kyp
   written by viktor k. decyk, ucla
   parallel, RISC optimized version
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny;
   int nxy, nxhy, kypt, j, k, nrx;
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nrxb, joff;
   float ani, at1;
   float complex s, t, t1, t2;
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
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,joff,at1,s,t,t1,t2)
   for (i = kypi-1; i < kypt; i++) {
      joff = 2*nxvh*i;
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
               t1 = s*f[2*j2+joff];
               t2 = s*f[1+2*j2+joff];
               f[2*j2+joff] = f[2*j1+joff] - t1;
               f[1+2*j2+joff] = f[1+2*j1+joff] - t2;
               f[2*j1+joff] += t1;
               f[1+2*j1+joff] += t2;
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
         t1 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < 2; k++) {
            t = conjf(f[k+2*(nxh-j)+joff]);
            s = f[k+2*j+joff] + t;
            t = (f[k+2*j+joff] - t)*t1;
            f[k+2*j+joff] = ani*(s + t);
            f[k+2*(nxh-j)+joff] = ani*conjf(s - t);
         }
      }
      for (k = 0; k < 2; k++) {
         f[k+joff] = 2.0*ani*((crealf(f[k+joff]) + cimagf(f[k+joff]))
                     + (crealf(f[k+joff]) - cimagf(f[k+joff]))*_Complex_I);
         if (nxhh > 0)
            f[k+2*nxhh+joff] = 2.0*ani*conjf(f[k+2*nxhh+joff]);
      }
   }
   return;
/* forward fourier transform */
L100: nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,joff,at1,s,t,t1,t2)
   for (i = kypi-1; i < kypt; i++) {
      joff = 2*nxvh*i;
/* scramble coefficients */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
         t1 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
         for (k = 0; k < 2; k++) {
            t = conjf(f[k+2*(nxh-j)+joff]);
            s = f[k+2*j+joff] + t;
            t = (f[k+2*j+joff] - t)*t1;
            f[k+2*j+joff] = s + t;
            f[k+2*(nxh-j)+joff] = conjf(s - t);
         }
      }
      for (k = 0; k < 2; k++) {
         f[k+joff] = (crealf(f[k+joff]) + cimagf(f[k+joff]))
                      + (crealf(f[k+joff]) - cimagf(f[k+joff]))*_Complex_I;
         if (nxhh > 0)
            f[k+2*nxhh+joff] = 2.0*conjf(f[k+2*nxhh+joff]);
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
               t1 = s*f[2*j2+joff];
               t2 = s*f[1+2*j2+joff];
               f[2*j2+joff] = f[2*j1+joff] - t1;
               f[1+2*j2+joff] = f[1+2*j1+joff] - t2;
               f[2*j1+joff] += t1;
               f[1+2*j1+joff] +=  t2;
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
void cppfft2rm2xy(float complex g[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int kstrt,
                  int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                  int nxyhd) {
/* this subroutine performs the y part of 2 two dimensional real to
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
   g[n][m][0:1] = sum(g[j][k][0:1]*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, a forward fourier transform is performed
   g[j][k][0:1] = sum(g[n][m][0:1]*exp(sqrt(-1)*2pi*m*k/ny))
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
   g[j][k][0:1] = mode jj,k, where jj = j + kxp*(kstrt - 1)
   0 <= jj < nx/2 and 0 <= k < ny, except for
   g[0][k][0:1] = mode nx/2,k, where ny/2+1 <= k < ny, and
   imaginary part of g[0][0][0:1] = real part of mode nx/2,0 and
   imaginary part of g[0][ny/2][0:1] = real part of mode nx/2,ny/2
   on node kstrt=0
   written by viktor k. decyk, ucla
   parallel, RISC optimized version
local data                                                            */
   int indx1, indx1y, nx, nxh, ny, nyh;
   int nxy, nxhy, ks, kxpt, j, k, nry;
   int i, m, ns, ns2, km, kmr, k1, k2, j1, j2, nryb, koff;
   float complex s, t1, t2;
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
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,koff,s,t1,t2)
   for (i = kxpi-1; i < kxpt; i++) {
      koff = 2*nyv*i;
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            t1 = g[2*k1+koff];
            t2 = g[1+2*k1+koff];
            g[2*k1+koff] = g[2*k+koff];
            g[1+2*k1+koff] = g[1+2*k+koff];
            g[2*k+koff] = t1;
            g[1+2*k+koff] = t2;
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
               t1 = s*g[2*j2+koff];
               t2 = s*g[1+2*j2+koff];
               g[2*j2+koff] = g[2*j1+koff] - t1;
               g[1+2*j2+koff] = g[1+2*j1+koff] - t2;
               g[2*j1+koff] += t1;
               g[1+2*j1+koff] += t2;
            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if ((ks==0) && (kxpi==1)) {
      for (k = 1; k < nyh; k++) {
         for (j = 0; j < 2; j++) {
            s = g[j+2*(ny-k)];
            g[j+2*(ny-k)] = 0.5*(cimagf(g[j+2*k] + s)
                             + crealf(g[j+2*k] - s)*_Complex_I);
            g[j+2*k] = 0.5*(crealf(g[j+2*k] + s)
                        + cimagf(g[j+2*k] - s)*_Complex_I);
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
         for (j = 0; j < 2; j++) {
            s = cimagf(g[j+2*(ny-k)])
                + crealf(g[j+2*(ny-k)])*_Complex_I;
            g[j+2*(ny-k)] = conjf(g[j+2*k] - s);
            g[j+2*k] += s;
         }
      }
   }
#pragma omp parallel for \
private(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,koff,s,t1,t2)
   for (i = kxpi-1; i < kxpt; i++) {
      koff = 2*nyv*i;
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            t1 = g[2*k1+koff];
            t2 = g[1+2*k1+koff];
            g[2*k1+koff] = g[2*k+koff];
            g[1+2*k1+koff] = g[1+2*k+koff];
            g[2*k+koff] = t1;
            g[1+2*k+koff] = t2;
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
               t1 = s*g[2*j2+koff];
               t2 = s*g[1+2*j2+koff];
               g[2*j2+koff] = g[2*j1+koff] - t1;
               g[1+2*j2+koff] = g[1+2*j1+koff] - t2;
               g[2*j1+koff] += t1;
               g[1+2*j1+koff] += t2;
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
void cwppfft2rm2(float complex f[], float complex g[],
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
      cppfft2rm2xx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,
                   kypd,nxhyd,nxyhd);
/* transpose f array to g */
      cpwtimera(-1,ttp,&dtime);
      cppntpose(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,2,nxvh,nyv,kxp,kypd);
      cpwtimera(1,ttp,&dtime);
/* perform y fft */
      cppfft2rm2xy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
                  nxhyd,nxyhd);
/* transpose g array to f */
      if (ntpose==0) {
         cpwtimera(-1,&tf,&dtime);
         cppntpose(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,2,nyv,nxvh,kypd,
                   kxp);
         cpwtimera(1,&tf,&dtime);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* transpose f array to g */
      if (ntpose==0) {
         cpwtimera(-1,&tf,&dtime);
         cppntpose(f,g,bs,br,nxh,ny,kxp,kyp,kstrt,nvp,2,nxvh,nyv,kxp,
                   kypd);
         cpwtimera(1,&tf,&dtime);
      }
/* perform y fft */
      cppfft2rm2xy(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,nyv,kxp,
                  nxhyd,nxyhd);
/* transpose g array to f */
      cpwtimera(-1,ttp,&dtime);
      cppntpose(g,f,br,bs,ny,nxh,kyp,kxp,kstrt,nvp,2,nyv,nxvh,kypd,kxp);
      cpwtimera(1,ttp,&dtime);
/* perform x fft */
      cppfft2rm2xx(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nxvh,
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
void cpdistr2_(float *part, float *edges, int *npp, int *nps, float *vtx,
               float *vty, float *vdx, float *vdy, int *npx, int *npy,
               int *nx, int *ny, int *idimp, int *npmax, int *idps,
               int *ipbc, int *ierr) {
   cpdistr2(part,edges,npp,*nps,*vtx,*vty,*vdx,*vdy,*npx,*npy,*nx,*ny,
            *idimp,*npmax,*idps,*ipbc,ierr);
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
void cppgppush2l_(float *ppart, float *fxy, int *kpic, int *noff,
                  int *nyp, float *qbm, float *dt, float *ek, int *nx,
                  int *ny, int *mx, int *my, int *idimp, int *nppmx,
                  int *nxv, int *nypmx, int *mx1, int *mxyp1, 
                  int *ipbc) {
   cppgppush2l(ppart,fxy,kpic,*noff,*nyp,*qbm,*dt,ek,*nx,*ny,*mx,*my,
               *idimp,*nppmx,*nxv,*nypmx,*mx1,*mxyp1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgppushf2l_(float *ppart, float *fxy, int *kpic, int *ncl,
                   int *ihole, int *noff, int *nyp, float *qbm,
                   float *dt, float *ek, int *nx, int *ny, int *mx,
                   int *my, int *idimp, int *nppmx, int *nxv,
                   int *nypmx, int *mx1, int *mxyp1, int *ntmax,
                   int *irc) {
   cppgppushf2l(ppart,fxy,kpic,ncl,ihole,*noff,*nyp,*qbm,*dt,ek,*nx,*ny,
                *mx,*my,*idimp,*nppmx,*nxv,*nypmx,*mx1,*mxyp1,*ntmax,
                irc);
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
void cppcguard2xl_(float *fxy, int *nyp, int *nx, int *ndim, int *nxe,
                   int *nypmx) {
   cppcguard2xl(fxy,*nyp,*nx,*ndim,*nxe,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppaguard2xl_(float *q, int *nyp, int *nx, int *nxe, int *nypmx) {
   cppaguard2xl(q,*nyp,*nx,*nxe,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cmppois22_(float complex *q, float complex *fxy, int *isign,
                float complex *ffc, float *ax, float *ay, float *affp,
                float *we, int *nx, int *ny, int *kstrt, int *nyv,
                int *kxp, int *nyhd) {
   cmppois22(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*kstrt,*nyv,*kxp,
             *nyhd);
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
void cppfft2rm2xx_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *kstrt,
                   int *kypi, int *kypp, int *nxvh, int *kypd,
                   int *nxhyd, int *nxyhd) {
   cppfft2rm2xx(f,*isign,mixup,sct,*indx,*indy,*kstrt,*kypi,*kypp,*nxvh,
                *kypd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rm2xy_(float complex *g, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *kstrt,
                   int *kxpi, int *kxpp, int *nyv, int *kxp, int *nxhyd,
                   int *nxyhd) {
   cppfft2rm2xy(g,*isign,mixup,sct,*indx,*indy,*kstrt,*kxpi,*kxpp,*nyv,
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
void cwppfft2rm2_(float complex *f, float complex *g, float complex *bs,
                  float complex *br, int *isign, int *ntpose,
                  int *mixup, float complex *sct, float *ttp, int *indx,
                  int *indy, int *kstrt, int *nvp, int *nxvh, int *nyv,
                  int *kxp, int *kyp, int *kypd, int *nxhyd,
                  int *nxyhd) {
   cwppfft2rm2(f,g,bs,br,*isign,*ntpose,mixup,sct,ttp,*indx,*indy,
               *kstrt,*nvp,*nxvh,*nyv,*kxp,*kyp,*kypd,*nxhyd,*nxyhd);
   return;
}
