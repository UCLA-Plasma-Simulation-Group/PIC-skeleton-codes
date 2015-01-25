/* CUDA Library for Skeleton 2D Electrostatic GPU PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include "cuda.h"

extern int nblock_size;
extern int maxgsx;

static cudaError_t crc;

/*--------------------------------------------------------------------*/
__device__ void liscan2(int *isdata, int nths) {
/* performs local prefix reduction of integer data shared by threads */
/* using binary tree method. */
/* local data */
   int l, mb, kxs, lb, kb;
   l = threadIdx.x;
   mb = l;
   kxs = 1;
   while (kxs < nths) {
      lb = kxs*mb;
      kb = 2*lb + kxs - 1;
      lb += l + kxs;
      if (lb < nths) {
         isdata[lb] += isdata[kb];
      }
      __syncthreads();
      mb >>= 1;
      kxs <<= 1;
   }
   return;
}

/*--------------------------------------------------------------------*/
__device__ void lsum2(float *sdata, int n) {
/* finds local sum of nths data items shared by threads */
/* using binary tree method. input is modified. */
/* local data */
   int l, k;
   float s;
   l = threadIdx.x;
   k = blockDim.x >> 1;
   s = 0.0f;
   if (l < n) s = sdata[l];
   while (k > 0) {
      if (l < k) {
         if ((l+k) < n) {
            s += sdata[l+k];
            sdata[l] = s;
         }
      }
      __syncthreads();
      k >>= 1;
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppush2l(float ppart[], float fxy[], int kpic[],
                           float qbm, float dt, float *ek, int idimp,
                           int nppmx, int nx, int ny, int mx, int my,
                           int nxv, int nyv, int mx1, int mxy1,
                           int ipbc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions.
   threaded version using guard cells
   data read in tiles
   particles stored segmented array
   44 flops/particle, 12 loads, 4 stores
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = velocity vx of particle n in tile m
   ppart[m][3][n] = velocity vy of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of field arrays, must be >= nx+1
   nyv = second dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
local data                                                            */
   int noff, moff, npoff, npp, mxv;
   int i, j, k, ii, nn, mm;
   float qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
/* The sizes of the shared memory arrays are as follows: */
/* float sfxy[2*(mx+1)*(my+1)], sek[blockDim.x];         */
/* to conserve memory, sek overlaps with sfxy            */
/* and the name sfxy is used instead of sek              */
   extern __shared__ float sfxy[];
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
   mxv = mx + 1;
/* k = tile number */
   k = blockIdx.x + gridDim.x*blockIdx.y;
/* loop over tiles */
   if (k < mxy1) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      ii = threadIdx.x;
      while (ii < mxv*(my+1)) {
         j = ii/mxv;
         i = ii - mxv*j;
         if ((i < nn) && (j < mm)) {
            sfxy[2*ii] = fxy[2*(i+noff+nxv*(j+moff))];
            sfxy[1+2*ii] = fxy[1+2*(i+noff+nxv*(j+moff))];
         }
         ii += blockDim.x;
      }
/* synchronize threads */
      __syncthreads();
/* loop over particles in tile */
      j = threadIdx.x;
      while (j < npp) {
/* find interpolation weights */
         x = ppart[j+npoff];
         nn = x;
         y = ppart[j+npoff+nppmx];
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nn = 2*(nn - noff) + 2*mxv*(mm - moff);
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find acceleration */
         dx = amx*sfxy[nn];
         dy = amx*sfxy[1+nn];
         dx = amy*(dxp*sfxy[2+nn] + dx);
         dy = amy*(dxp*sfxy[3+nn] + dy);
         nn += 2*mxv;
         vx = amx*sfxy[nn];
         vy = amx*sfxy[1+nn];
         dx += dyp*(dxp*sfxy[2+nn] + vx);
         dy += dyp*(dxp*sfxy[3+nn] + vy);
/* new velocity */
         vx = ppart[j+npoff+nppmx*2];
         vy = ppart[j+npoff+nppmx*3];
         dx = vx + qtm*dx;
         dy = vy + qtm*dy;
/* average kinetic energy */
         vx += dx;
         vy += dy;
         sum1 += (double) (vx*vx + vy*vy);
         ppart[j+npoff+nppmx*2] = dx;
         ppart[j+npoff+nppmx*3] = dy;
/* new position */
         dx = x + dx*dt;
         dy = y + dy*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = ppart[j+npoff];
               ppart[j+npoff+nppmx*2] = -ppart[j+npoff+nppmx*2];
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = ppart[j+npoff+nppmx];
               ppart[j+npoff+nppmx*3] = -ppart[j+npoff+nppmx*3];
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = ppart[j+npoff];
               ppart[j+npoff+nppmx*2] = -ppart[j+npoff+nppmx*2];
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+npoff+nppmx] = dy;
         j += blockDim.x;
      }
/* synchronize threads */
      __syncthreads();
/* add kinetic energies in tile */
      sfxy[threadIdx.x] = (float) sum1;
/* synchronize threads */
      __syncthreads();
      lsum2(sfxy,blockDim.x);
/* normalize kinetic energy of tile */
      if (threadIdx.x==0) {
         ek[k] = 0.125f*sfxy[0];
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppushf2l(float ppart[], float fxy[], int kpic[],
                            int ncl[], int ihole[], float qbm, float dt,
                            float *ek, int idimp, int nppmx, int nx,
                            int ny, int mx, int my, int nxv, int nyv,
                            int mx1, int mxy1, int ntmax, int *irc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   threaded version using guard cells
   data read in tiles
   particles stored segmented array
   44 flops/particle, 12 loads, 4 stores
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = velocity vx of particle n in tile m
   ppart[m][3][n] = velocity vy of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of field arrays, must be >= nx+1
   nyv = second dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
local data                                                            */
   int noff, moff, npoff, nhoff, mhoff, npp, mxv;
   int i, j, k, ii, ih, nn, mm;
   float qtm, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
   float anx, any, edgelx, edgely, edgerx, edgery;
/* The sizes of the shared memory arrays are as follows: */
/* float sfxy[2*(mx+1)*(my+1)], sek[blockDim.x];         */
/* int sih[blockDim.x], sncl[8], nh[1];                  */
/* to conserve memory, sek overlaps with sfxy and sih    */
/* and the name sfxy is used instead of sek              */
   float *sfxy;
   int *sncl, *sih, *nh;
   extern __shared__ int shm[];
   sfxy = (float *)&shm[0];
   sih = (int *)&sfxy[2*(mx+1)*(my+1)];
   sncl = (int *)&sih[blockDim.x];
   nh = (int *)&sfxy[blockDim.x];
   sncl = sncl > nh ? sncl : nh;
   nh = (int *)&sncl[8];
   double sum1;
   qtm = qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   sum1 = 0.0;
   mxv = mx + 1;
/* k = tile number */
   k = blockIdx.x + gridDim.x*blockIdx.y;
/* loop over tiles */
   if (k < mxy1) {
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
/* load local fields from global array */
      ii = threadIdx.x;
      while (ii < mxv*(my+1)) {
         j = ii/mxv;
         i = ii - mxv*j;
         if ((i < nn+1) && (j < mm+1)) {
            sfxy[2*ii] = fxy[2*(i+noff+nxv*(j+moff))];
            sfxy[1+2*ii] = fxy[1+2*(i+noff+nxv*(j+moff))];
         }
         ii += blockDim.x;
      }
/* clear counters */
      j = threadIdx.x;
      while (j < 8) {
         sncl[j] = 0;
         j += blockDim.x;
      }
      if (threadIdx.x==0) {
         nh[0] = 0;
      }
/* synchronize threads */
      __syncthreads();
/* loop over particles in tile */
      ii = (npp - 1)/(int) blockDim.x + 1;
      nhoff = 0;
      for (i = 0; i < ii; i++) {
         j = threadIdx.x + blockDim.x*i;
         sih[threadIdx.x] = 0;
         if (j < npp) {
/* find interpolation weights */
            x = ppart[j+npoff];
            nn = x;
            y = ppart[j+npoff+nppmx];
            mm = y;
            dxp = x - (float) nn;
            dyp = y - (float) mm;
            nn = 2*(nn - noff) + 2*mxv*(mm - moff);
            amx = 1.0f - dxp;
            amy = 1.0f - dyp;
/* find acceleration */
            dx = amx*sfxy[nn];
            dy = amx*sfxy[1+nn];
            dx = amy*(dxp*sfxy[2+nn] + dx);
            dy = amy*(dxp*sfxy[3+nn] + dy);
            nn += 2*mxv;
            vx = amx*sfxy[nn];
            vy = amx*sfxy[1+nn];
            dx += dyp*(dxp*sfxy[2+nn] + vx);
            dy += dyp*(dxp*sfxy[3+nn] + vy);
/* new velocity */
            vx = ppart[j+npoff+nppmx*2];
            vy = ppart[j+npoff+nppmx*3];
            dx = vx + qtm*dx;
            dy = vy + qtm*dy;
/* average kinetic energy */
            vx += dx;
            vy += dy;
            sum1 += (double) (vx*vx + vy*vy);
            ppart[j+npoff+nppmx*2] = dx;
            ppart[j+npoff+nppmx*3] = dy;
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
                     dx = 0.0f;
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
/* using prefix scan for ih to keep holes ordered */
            if (mm > 0) {
               atomicAdd(&sncl[mm-1],1);
               sih[threadIdx.x] = 1;
            }
         }
/* synchronize threads */
         __syncthreads();
         nn = npp - blockDim.x*i;
         if (nn > blockDim.x)
            nn = blockDim.x;
/* perform local prefix reduction */
         liscan2(sih,nn);
         if (j < npp) {
/* set new position */
            ppart[j+npoff] = dx;
            ppart[j+npoff+nppmx] = dy;
/* write out location and direction of departing particles */
            ih = sih[threadIdx.x];
            mhoff = 0;
            if (threadIdx.x > 0)
               mhoff = sih[threadIdx.x-1];
/* this thread has a hole present */
            if (ih > mhoff) {
               ih += nhoff;
               if (ih <= ntmax) {
                  ihole[2*(ih+(ntmax+1)*k)] = j + 1;
                  ihole[1+2*(ih+(ntmax+1)*k)] = mm;
               }
               else {
                  nh[0] = 1;
               }
            }
         }
/* update number of holes in this iteration */
         if (nn > 0)
            nhoff += sih[nn-1];
/* synchronize threads */
         __syncthreads();
      }
/* add kinetic energies in tile */
      sfxy[threadIdx.x] = (float) sum1;
/* synchronize threads */
      __syncthreads();
      lsum2(sfxy,blockDim.x);
/* write out counters */
      j = threadIdx.x;
      while (j < 8) {
         ncl[j+8*k] = sncl[j];
         j += blockDim.x;
      }
/* set error and end of file flag */
      if (threadIdx.x==0) {
/* ihole overflow */
         ih  = nhoff;
         if (nh[0] > 0) {
            *irc = ih;
            ih = -ih;
         }
         ihole[2*(ntmax+1)*k] = ih;
/* normalize kinetic energy of tile */
         ek[k] = 0.125f*sfxy[0];
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpu2ppost2l(float ppart[], float q[], int kpic[],
                            float qm, int nppmx, int idimp, int mx,
                            int my, int nxv, int nyv, int mx1,
                            int mxy1) {
/* for 2d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   threaded version using guard cells
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
   int noff, moff, npoff, npp, mxv;
   int i, j, k, ii, nn, mm, np, mp;
   float dxp, dyp, amx, amy;
/* The size of the shared memory array is as follows: */
/* float sq[(mx+1)*(my+1)]                            */
   extern __shared__ float sq[];
   mxv = mx + 1;
/* k = tile number */
   k = blockIdx.x + gridDim.x*blockIdx.y;
/* loop over tiles */
   if (k < mxy1) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
/* zero out local accumulator */
      i = threadIdx.x;
      while (i < mxv*(my+1)) {
         sq[i] = 0.0f;
         i += blockDim.x;
      }
/* synchronize threads */
       __syncthreads();
/* loop over particles in tile */
      j = threadIdx.x;
      while (j < npp) {
/* find interpolation weights */
         dxp = ppart[j+npoff];
         nn = dxp;
         dyp = ppart[j+npoff+nppmx];
         mm = dyp;
         dxp = qm*(dxp - (float) nn);
         dyp = dyp - (float) mm;
         nn = nn - noff;
         mm = mxv*(mm - moff);
         amx = qm - dxp;
         mp = mm + mxv;
         amy = 1.0f - dyp;
         np = nn + 1;
/* deposit charge within tile to local accumulator */
/* original deposit charge, has data hazard on GPU */
/*       sq[np+mp] += dxp*dyp; */
/*       sq[nn+mp] += amx*dyp; */
/*       sq[np+mm] += dxp*amy; */
/*       sq[nn+mm] += amx*amy; */
/* for devices with compute capability 2.x */
         atomicAdd(&sq[np+mp],dxp*dyp);
         atomicAdd(&sq[nn+mp],amx*dyp);
         atomicAdd(&sq[np+mm],dxp*amy);
         atomicAdd(&sq[nn+mm],amx*amy);
         j += blockDim.x;
      }
/* synchronize threads */
      __syncthreads();
/* deposit charge to global array */
      nn = mxv < nxv-noff ? mxv : nxv-noff;
      mm = my+1 < nyv-moff ? my+1 : nyv-moff;
      ii = threadIdx.x;
      while (ii < mxv*(my+1)) {
         j = ii/mxv;
         i = ii - mxv*j;
         if ((i < nn) && (j < mm)) {
/* original deposit charge, has data hazard on GPU */
/*          q[i+noff+nxv*(j+moff)] += sq[ii]; */
/* for devices with compute capability 2.x */
            atomicAdd(&q[i+noff+nxv*(j+moff)],sq[ii]);
         }
         ii += blockDim.x;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpucaguard2l(float2 qc[], float q[], int nx, int ny,
                             int nxe, int nye, int nxvh, int nyv) {
/* copy and accumulate extended periodic scalar field q 
   into complex output field qc
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of input field array q, must be >= nx+1
   nye = second dimension of input field array q, must be >= ny+1
   nxvh = first dimension of output field array qc, must be >= nx/2+1
   nyv = second dimension of output field array qc, must be >= ny     */
/* local data */
   int j, k, nxh;
   float at1, at2;
   float2 a;
   nxh = nx/2;
   k = blockIdx.x;
   if (k < ny) {
      j = threadIdx.x;
      at2 = 0.0f;
      while (j < nxh) {
         if (k==0) {
            at1 = q[2*j+nxe*ny];
            at2 = q[2*j+1+nxe*ny];
            if (j==0) {
               at1 += q[nx] + q[nx+nxe*ny];
            }
         }
         if (k > 0) {
            at1 = 0.0f;
            if (j==0) {
               at1 = q[nx+nxe*k];
            }
         }
         a.x = q[2*j+nxe*k] + at1;
         a.y = q[2*j+1+nxe*k] + at2;
         qc[j+nxvh*k] = a;
         j += blockDim.x;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuccguard2l(float2 fxyc[], float fxy[], int nx, int ny,
                             int nxe, int nye, int nxvh, int nyv) {
/* copy and replicate complex input 2d vector field fxyc
   into extended periodic field fxy
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = second dimension of input field array fxy, must be >= nx+1
   nye = third dimension of input field array fxy, must be >= ny+1
   nxvh = first dimension of input field array fxyc, must be >= nx/2+1
   nyv = third dimension of input field array fxyc, must be >= ny     */
/* local data */
   int j, k, nxh;
   float2 a, b;
   nxh = nx/2;
   k = blockIdx.x;
/* copy interior points */
   if (k < ny) {
      j = threadIdx.x;
      while (j < nxh) {
         a = fxyc[j+nxvh*2*k];
         b = fxyc[j+nxvh*(1+2*k)];
         fxy[2*(2*j+nxe*k)] = a.x;
         fxy[1+2*(2*j+nxe*k)] = b.x;
         fxy[2*(2*j+1+nxe*k)] = a.y;
         fxy[1+2*(2*j+1+nxe*k)] = b.y;
         j += blockDim.x;
      }
   }
/* accumulate edges of extended field */
   if (blockIdx.x==0) {
      k = threadIdx.x;
      while (k < ny) {
         a = fxyc[nxvh*2*k];
         b = fxyc[nxvh*(1+2*k)];
         fxy[2*(nx+nxe*k)] = a.x;
         fxy[1+2*(nx+nxe*k)] = b.x;
         k += blockDim.x;
      }
      j = threadIdx.x;
      while (j < nxh) {
         a = fxyc[j];
         b = fxyc[j+nxvh];
         fxy[2*(2*j+nxe*ny)] = a.x;
         fxy[1+2*(2*j+nxe*ny)] = b.x;
         fxy[2*(2*j+1+nxe*ny)] = a.y;
         fxy[1+2*(2*j+1+nxe*ny)] = b.y;
         j += blockDim.x;
      }
      if (threadIdx.x==0) {
         a = fxyc[0];
         b = fxyc[nxvh*1];
         fxy[2*(nx+nxe*ny)] = a.x;
         fxy[1+2*(nx+nxe*ny)] = b.x;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppfnd2l(float ppart[], int kpic[], int ncl[],
                           int ihole[], int idimp, int nppmx, int nx,
                           int ny, int mx, int my, int mx1, int my1,
                           int ntmax, int *irc) {
/* this subroutine performs first step of a particle sort by x,y grid
   in tiles of mx, my, where one finds the particles leaving tile and
   stores their number, location, and destination in ncl and ihole.
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 2D linear memory
   input: all except ncl, ihole, irc
   output: ppart, ncl, ihole, irc
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k 
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int mxy1, noff, moff, npp, j, k, ih, ist, nn, mm, nths;
   float anx, any, edgelx, edgely, edgerx, edgery, dx, dy;
/* The sizes of the shared memory arrays are as follows: */
/* int sncl[8], sih[blockDim.x], nh[1];                  */
   int *sncl, *sih, *nh;
   extern __shared__ int shm[];
   sncl = (int *)&shm[0];
   sih = (int *)&shm[8];
   nh = (int *)&shm[8+blockDim.x];
   mxy1 = mx1*my1;
   anx = (float) nx;
   any = (float) ny;
/* k = tile number */
   k = blockIdx.x + gridDim.x*blockIdx.y;
/* find and count particles leaving tiles and determine destination */
/* update ppart, ihole, ncl */
/* loop over tiles */
   if (k < mxy1) {
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
/* clear counters */
      j = threadIdx.x;
      while (j < 8) {
         sncl[j] = 0;
         j += blockDim.x;
      }
      if (threadIdx.x==0) {
         nh[0] = 0;
      }
/* synchronize threads */
      __syncthreads();
/* loop over particles in tile */
      mm = (npp - 1)/(int) blockDim.x + 1;
      noff = 0;
      for (nn = 0; nn < mm; nn++) {
         j = threadIdx.x + blockDim.x*nn;
         sih[threadIdx.x] = 0;
         if (j < npp) {
            dx = ppart[j+nppmx*(idimp*k)];
            dy = ppart[j+nppmx*(1+idimp*k)];
/* find particles going out of bounds */
            ist = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* ist = direction particle is going                             */
            if (dx >= edgerx) {
               if (dx >= anx)
                  ppart[j+nppmx*(idimp*k)] = dx - anx;
               ist = 2;
            }
            else if (dx < edgelx) {
               if (dx < 0.0f) {
                  dx += anx;
                  if (dx < anx)
                     ist = 1;
                  else
                     dx = 0.0f;
                  ppart[j+nppmx*(idimp*k)] = dx;
               }
               else {
                  ist = 1;
               }
            }
            if (dy >= edgery) {
               if (dy >= any)
                  ppart[j+nppmx*(1+idimp*k)] = dy - any;
               ist += 6;
            }
            else if (dy < edgely) {
               if (dy < 0.0f) {
                  dy += any;
                  if (dy < any)
                     ist += 3;
                  else
                     dy = 0.0f;
                  ppart[j+nppmx*(1+idimp*k)] = dy;
               }
               else {
                  ist += 3;
               }
            }
/* using prefix scan for ih to keep holes ordered */
            if (ist > 0) {
               atomicAdd(&sncl[ist-1],1);
               sih[threadIdx.x] = 1;
            }
         }
/* synchronize threads */
         __syncthreads();
         nths = npp - blockDim.x*nn;
         if (nths > blockDim.x)
            nths = blockDim.x;
/* perform local prefix reduction */
         liscan2(sih,nths);
         if (j < npp) {
            ih = sih[threadIdx.x];
            moff = 0;
            if (threadIdx.x > 0)
               moff = sih[threadIdx.x-1];
/* this thread has a hole present */
            if (ih > moff) {
               ih += noff;
               if (ih <= ntmax) {
                  ihole[2*(ih+(ntmax+1)*k)] = j + 1;
                  ihole[1+2*(ih+(ntmax+1)*k)] = ist;
               }
               else {
                  nh[0] = 1;
               }
            }
         }
/* update number of holes in this iteration */
         if (nths > 0)
            noff += sih[nths-1];
/* synchronize threads */
         __syncthreads();
      }
/* write out counters */
      j = threadIdx.x;
      while (j < 8) {
         ncl[j+8*k] = sncl[j];
         j += blockDim.x;
      }
/* set error and end of file flag */
      if (threadIdx.x==0) {
/* ihole overflow */
         ih  = noff;
         if (nh[0] > 0) {
            *irc = ih;
            ih = -ih;
         }
         ihole[2*(ntmax+1)*k] = ih;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppmov2l(float ppart[], float ppbuff[], int ncl[],
                           int ihole[], int idimp, int nppmx, int mx1,
                           int my1, int npbmx, int ntmax, int *irc) {
/* this subroutine performs second step of a particle sort by x,y grid
   in tiles of mx, my, where prefix scan of ncl is performed and
   departing particles are buffered in ppbuff in direction order.
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 2D linear memory
   input: all except ppbuff, irc
   output: ppbuff, ncl, irc
   ppart[k][i][n] = i co-ordinate of particle n in tile k 
   ppbuff[k][i][n] = i co-ordinate of particle n in tile k
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
   int mxy1, i, j, k, ii, nh, ist, j1, ierr;
/* The sizes of the shared memory arrays are as follows: */
/* int sncl[8], ip[1];                                   */
/* blockDim.x should be >= 8                             */
   int *sncl, *ip;
   extern __shared__ int shm[];
   sncl = (int *)&shm[0];
   ip = (int *)&shm[8];
   mxy1 = mx1*my1;
   ierr = 0;
/* k = tile number */
   k = blockIdx.x + gridDim.x*blockIdx.y;
   j = threadIdx.x;
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
   if (k < mxy1) {
/* find address offset for ordered ppbuff array */
      if (j < 8) {
         ist = ncl[j+8*k];
         sncl[j] = ist;
      }
      if (threadIdx.x==0)
         ip[0] = 0;
/* synchronize threads */
      __syncthreads();
/* perform local prefix reduction */
      liscan2(sncl,8);
      if (j < 8)
         sncl[j] -= ist;
/* synchronize threads */
      __syncthreads();
      nh = ihole[2*(ntmax+1)*k];
/* loop over particles leaving tile */
      while (j < nh) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+(ntmax+1)*k)] - 1;
         ist = ihole[1+2*(j+1+(ntmax+1)*k)];
         ii = atomicAdd(&sncl[ist-1],1);
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*(i+idimp*k)]
               = ppart[j1+nppmx*(i+idimp*k)];
            }
         }
         else {
            ip[0] = 1;
         }
         j += blockDim.x;
      }
/* synchronize threads */
      __syncthreads();
/* write out counters */
      j = threadIdx.x;
      if (j < 8) {
         ncl[j+8*k] = sncl[j];
      }
/* set error */
      if (threadIdx.x==0) {
         if (ip[0] > 0)
            ierr = ierr > sncl[7] ? ierr : sncl[7];
      }
   }
/* ppbuff overflow */
   if (ierr > 0)
      *irc = ierr;
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppord2l(float ppart[], float ppbuff[], int kpic[],
                           int ncl[], int ihole[], int idimp, int nppmx,
                           int mx1, int my1, int npbmx, int ntmax,
                           int *irc) {
/* this subroutine performs third step of a particle sort by x,y grid
   in tiles of mx, my, where incoming particles from other tiles are
   copied into ppart.
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 2D linear memory
   input: all except irc
   output: ppart, kpic, irc
   ppart[k][i][n] = i co-ordinate of particle n in tile k 
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
   int mxy1, npp, ncoff, i, j, k, ii, jj, kx, ky, ni, nh;
   int nn, mm, ll, ip, j1, j2, kxl, kxr, kk, kl, kr;
   int nths;
/* The sizes of the shared memory arrays are as follows: */
/* int ks[8], sip[8], sj[blockDim.x], sj1[1], ist[1];    */
   int *ks, *sip, *sj, *sj1, *ist;
   extern __shared__ int shm[];
   ks = (int *)&shm[0];
   sip = (int *)&shm[8];
   sj = (int *)&shm[16];
   sj1 = (int *)&shm[16+blockDim.x];
   ist = (int *)&shm[17+blockDim.x];
   mxy1 = mx1*my1;
/* k = tile number */
   k = blockIdx.x + gridDim.x*blockIdx.y;
/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
   if (k < mxy1) {
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
      if (threadIdx.x==0) {
         ks[0] = kxr + kk;
         ks[1] = kxl + kk;
         ks[2] = kx + kr;
         ks[3] = kxr + kr;
         ks[4] = kxl + kr;
         ks[5] = kx + kl;
         ks[6] = kxr + kl;
         ks[7] = kxl + kl;
         sj1[0] = 0;
         ist[0] = 0;
      }
/* synchronize threads */
      __syncthreads();
/* find number of incoming particles */
      kk = 0;
      ncoff = 0;
      ip = 0;
      ii = threadIdx.x;
      if (ii < 8) {
         kk = ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+8*kk];
         ip = ncl[ii+8*kk] - ncoff;
         kk = ncoff + idimp*npbmx*kk;
         sip[ii] = ip;
      }
/* synchronize threads */
      __syncthreads();
/* perform local prefix reduction */
      liscan2(sip,8);
      ni = sip[7];
/* loop over directions */
      nh = ihole[2*(ntmax+1)*k];
      j1 = 0;
      mm = (ni - 1)/(int) blockDim.x + 1;
      for (nn = 0; nn < mm; nn++) {
         j = threadIdx.x + blockDim.x*nn;
         sj[threadIdx.x] = 0;
         if (threadIdx.x==0)
            sj[0] = sj1[0];
/* synchronize threads */
         __syncthreads();
/* calculate offset for reading from particle buffer */
         if (ii < 8) {
/* mark next location where direction ii changes */
            jj = sip[ii] - blockDim.x*nn;
            if ((jj >= 0) && (jj < blockDim.x)) {
               if (ip > 0)
                  sj[jj] -= kk + ip;
            }
         }
/* synchronize threads */
         __syncthreads();
/* calculate offset for reading from particle buffer */
         if (ii < 8) {
/* mark location where direction ii starts */
            jj -= ip;
            if ((jj >= 0) && (jj < blockDim.x)) {
               if (ip > 0)
                  sj[jj] += kk;
            }
         }
         nths = ni - blockDim.x*nn;
         if (nths > blockDim.x)
            nths = blockDim.x;
/* synchronize threads */
         __syncthreads();
/* perform local prefix reduction */
         liscan2(sj,nths);
/* save last value for next time */
         if (threadIdx.x==0) {
            jj = 0;
            if (nths > 0)
               jj = sj[nths-1];
            sj1[0] = jj;
         }
         if (j < ni) {
/* insert incoming particles into holes */
            if (j < nh) {
               j1 = ihole[2*(j+1+(ntmax+1)*k)] - 1;
            }
/* place overflow at end of array */
            else {
               j1 = npp + (j - nh);
            }
            if (j1 < nppmx) {
               jj = sj[threadIdx.x];
               for (i = 0; i < idimp; i++) {
                  ppart[j1+nppmx*(i+idimp*k)]
                  = ppbuff[j+jj+npbmx*i];
                }
            }
            else {
               ist[0] = 1;
            }
         }
/* synchronize threads */
         __syncthreads();
      }
/* update particle number if all holes have been filled */
      jj = ni - nh;
      if (jj > 0)
         npp += jj;
/* fill up remaining holes in particle array with particles from end */
      ip = nh - ni;
      if (ip > 0) {
         mm = (ip - 1)/(int) blockDim.x + 1;
         kk = 0;
         ll = 0;
/* loop over holes */
         for (nn = 0; nn < mm; nn++) {
            j = threadIdx.x + blockDim.x*nn;
/* j1 = locations of particles to fill holes, in decreasing order */
            j1 = 0;
            if (j < ip) {
               j1 = npp - j - 1;
            }
/* j2 = locations of holes at the end, in decreasing order */
            j2 = 0;
            jj = nh - ll - threadIdx.x;
            if (jj > 0) {
               j2 = ihole[2*(jj+(ntmax+1)*k)] - 1;
            }
/* holes with locations greater than npp-ip do not need to be filled */
/* identify such holes */
            sj[threadIdx.x] = 1;
/* synchronize threads */
            __syncthreads();
/* omit particles at end that are holes */
            ii = npp - (j2 + blockDim.x*nn) - 1;
            if ((ii >= 0) && (ii < blockDim.x))
               sj[ii] = 0;
            nths = ip - blockDim.x*nn;
            if (nths > blockDim.x)
               nths = blockDim.x;
/* synchronize threads */
            __syncthreads();
/* perform local prefix reduction */
            liscan2(sj,nths);
/* ii = number particles at end to be moved */
            ii = 0;
            if (nths > 0)
               ii = sj[nths-1];
/* identify which particles at end to be moved */
            if (ii < nths) {
               ncoff = 0;
               if (j < ip) {
                  if (threadIdx.x > 0)
                    ncoff = sj[threadIdx.x-1];
                  jj = sj[threadIdx.x];
               }
/* synchronize threads */
               __syncthreads();
               if (j < ip) {
                  if (jj > ncoff) {
                     sj[jj-1] = j1;
                  }
               }
/* synchronize threads */
               __syncthreads();
            }
/* j2 = locations of holes to be filled in increasing order */
            j2 = 0;
            if (j < ip) {
               j1 = npp - j - 1;
               jj = threadIdx.x + ni + kk + 1;
               if (jj <= nh)
                  j2 = ihole[2*(jj+(ntmax+1)*k)] - 1;
            }
/* move particles from end into remaining holes */
            if (j < (ii+blockDim.x*nn)) {
               if (ii < nths)
                  j1 = sj[threadIdx.x];
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*(i+idimp*k)]
                  = ppart[j1+nppmx*(i+idimp*k)];
               }
            }
/* accumulate number of holes filled */
            kk += ii;
/* accumulate number of holes skipped over */
            ii = nths - ii;
            ll += ii;
         }
/* update number of particles */
         npp -= ip;
      }
/* set error and update particle */
      if (threadIdx.x==0) {
/* ppart overflow */
         if (ist[0] > 0)
            *irc = npp;
         kpic[k] = npp;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpupois22t(float2 qt[], float2 fxyt[], float2 ffct[],
                           float *we, int nx, int ny, int nxvh, int nyv,
                           int nxhd, int nyhd) {
/* this subroutine solves 2d poisson's equation in fourier space for
   force/charge (or convolution of electric field over particle shape)
   with periodic boundary conditions, without packed data.
   vector length is second dimension
   input: qt,ffct,nx,ny,nxvh,nyv,nxhd,nyhd, output: fxyt,we
   approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1
   equation used is:
   fx[kx][ky] = -sqrt(-1)*kx*g[kx][ky]*s[kx][ky]*q[kx][ky],
   fy[kx][ky] = -sqrt(-1)*ky*g[kx][ky]*s[kx][ky]*q[kx][ky],
   where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
   g[kx][ky] = (affp/(kx**2+ky**2))*s[kx,ky],
   s[kx][ky] = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
   fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
   fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
   qt[j][k] = complex charge density for fourier mode (k,j)
   fxyt[j][0][k] = x component of complex force/charge,
   fxyt[j][1][k] = y component of complex force/charge,
   all for fourier mode (k,j)
   caimag(ffct[j][k]) = finite-size particle shape factor s
   creal(ffct([j][k])) = potential green's function g
   for fourier mode (k,j)
   electric field energy is also calculated, using
   we = nx*ny*sum((affp/(kx**2+ky**2))*|q[kx][ky]*s[kx][ky]|**2)
   nx/ny = system length in x/y direction
   nxvh = second dimension of field arrays, must be >= nxh+1
   nyv = first dimension of field arrays, must be >= ny
   nxhd = second dimension of form factor array, must be >= nxh
   nyhd = first dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, nxh1, j, k, k1, jj, jk, jk2;
   float dnx, dny, dkx, at1, at2, at3, at4;
   float2 zero, zt1, zt2, zt3;
/* The size of the shared memory array is as follows: */
/* float ss[blockDim.x];                              */
   extern __shared__ float ss[];
   double wp;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959f/(float) nx;
   dny = 6.28318530717959f/(float) ny;
   zero.x = 0.0f;
   zero.y = 0.0f;
/* calculate force/charge and sum field energy */
   wp = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
/* for (j = 1; j < nxh; j++) { */
   j = blockIdx.x;
   if ((j > 0) && (j < nxh)) {
      dkx = dnx*(float) j;
      jj = nyhd*j;
      jk = nyv*j;
      jk2 = 2*jk;
/*    for (k = 1; k < nyh; k++) { */
      k = threadIdx.x;
      while (k < nyh) {
         if (k > 0) {
            k1 = ny - k;
            zt1 = ffct[k+jj];
            at1 = zt1.x*zt1.y;
            at2 = at1*dkx;
            at3 = at1*dny*(float) k;
            zt1 = qt[k+jk];
            at4 = zt1.x;
            zt1.x = zt1.y;
            zt1.y = -at4;
            zt2 = qt[k1+jk];
            at4 = zt2.x;
            zt2.x = zt2.y;
            zt2.y = -at4;
            zt3.x = at2*zt1.x;
            zt3.y = at2*zt1.y;
            fxyt[k+jk2] = zt3;
            zt3.x = at3*zt1.x;
            zt3.y = at3*zt1.y;
            fxyt[k+nyv+jk2] = zt3;
            zt3.x = at2*zt2.x;
            zt3.y = at2*zt2.y;
            fxyt[k1+jk2] = zt3;
            zt3.x = -at3*zt2.x;
            zt3.y = -at3*zt2.y;
            fxyt[k1+nyv+jk2] = zt3;
            wp += (double) (at1*(zt1.x*zt1.x + zt1.y*zt1.y
                + zt2.x*zt2.x + zt2.y*zt2.y));
         }
         k += blockDim.x;
      }
   }
/* mode numbers ky = 0, ny/2 */
   if (blockIdx.x==0) {
      k1 = nyh;
/*    for (j = 1; j < nxh; j++) { */
      j = threadIdx.x;
      while (j < nxh) {
         if (j > 0) {
            jj = nyhd*j;
            jk = nyv*j;
            jk2 = 2*jk;
            zt1 = ffct[jj];
            at1 = zt1.x*zt1.y;
            at2 = at1*dnx*(float) j;
            zt1 = qt[jk];
            at4 = zt1.x;
            zt3.x = at2*zt1.y;
            zt3.y = -at2*at4;
            fxyt[jk2] = zt3;
            fxyt[nyv+jk2] = zero;
            fxyt[k1+jk2] = zero;
            fxyt[k1+nyv+jk2] = zero;
            wp += (double) (at1*(zt1.x*zt1.x + zt1.y*zt1.y));
         }
         j += blockDim.x;
      }
/* mode numbers kx = 0, nx/2 */
      nxh1 = 2*nyv*nxh;
/*    for (k = 1; k < nyh; k++) { */
      k = threadIdx.x;
      while (k < nyh) {
         if (k > 0) {
            k1 = ny - k;
            zt1 = ffct[k];
            at1 = zt1.x*zt1.y;
            at3 = at1*dny*(float) k;
            zt1 = qt[k];
            at4 = zt1.x;
            zt3.x = at3*zt1.y;
            zt3.y = -at3*at4;
            fxyt[k] = zero;
            fxyt[k+nyv] = zt3;
            fxyt[k1] = zero;
            zt3.y = -zt3.y;
            fxyt[k1+nyv] = zt3;
            fxyt[k+nxh1] = zero;
            fxyt[k+nyv+nxh1] = zero;
            fxyt[k1+nxh1] = zero;
            fxyt[k1+nyv+nxh1] = zero;
            wp += (double) (at1*(zt1.x*zt1.x + zt1.y*zt1.y));
         }
         k += blockDim.x;
      }
      if (threadIdx.x==0) {
         k1 = nyh;
         fxyt[0] = zero;
         fxyt[nyv] = zero;
         fxyt[k1] = zero;
         fxyt[k1+nyv] = zero;
         fxyt[nxh1] = zero;
         fxyt[nxh1+nyv] = zero;
         fxyt[k1+nxh1] = zero;
         fxyt[k1+nyv+nxh1] = zero;
      }
   }
   j = blockIdx.x;
   if (j <= nxh) {
/* sum potential energies for each x co-ordinate */
      ss[threadIdx.x] = (float) wp;
/* synchronize threads */
      __syncthreads();
      lsum2(ss,blockDim.x);
/* normalize potential energy for each x co-ordinate */
      if (threadIdx.x==0)
         we[j] = ss[0]*((float) (nx*ny));
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuctpose4(float2 f[], float2 g[], int nx, int ny,
                           int nxv, int nyv) {
/* complex transpose using blocking algorithm with gaps */
/* local data */
   int j, k, js, ks, joff, koff, mx, mxv;
/* The size of the shared memory array is as follows: */
/* float2 shm2[(mx + 1)*mx];                          */
   extern __shared__ float2 shm2[];
   mx = blockDim.x;
   mxv = mx + 1;
   joff = mx*blockIdx.x;
   koff = mx*blockIdx.y;
   js = threadIdx.x;
   ks = threadIdx.y;
/* copy into block */
   j = js + joff;
   k = ks + koff;
   if ((j < nx) && (k < ny)) {
      shm2[js+mxv*ks] = f[j+nxv*k];
   }
   __syncthreads();
/* copy out from block */
   j = ks + joff;
   k = js + koff;
   if ((j < nx) && (k < ny)) {
      g[k+nyv*j] = shm2[ks+mxv*js];
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuctpose4n(float2 fn[], float2 gn[], int nx, int ny,
                            int ndim, int nxv, int nyv) {
/* complex vector transpose using blocking algorithm with gaps */
/* ndim = vector dimension                                     */
/* local data */
   int i, j, k, js, ks, joff, koff, mx, mxv, nmxv, nnxv, nnyv, jj, kk;
/* The size of the shared memory array is as follows: */
/* float2 shmn2[ndim*(mx + 1)*mx];                    */
   extern __shared__ float2 shmn2[];
   mx = blockDim.x;
   mxv = mx + 1;
   joff = mx*blockIdx.x;
   koff = mx*blockIdx.y;
   js = threadIdx.x;
   ks = threadIdx.y;
   nmxv = ndim*mxv;
   nnxv = ndim*nxv;
   nnyv = ndim*nyv;
/* copy into block */
   j = js + joff;
   k = ks + koff;
   if ((j < nx) && (k < ny)) {
      jj = j + nnxv*k;
      kk = js + nmxv*ks;
      for (i = 0; i < ndim; i++) {
         shmn2[kk+mxv*i] = fn[jj+nxv*i];
      }
   }
   __syncthreads();
/* copy out from block */
   j = ks + joff;
   k = js + koff;
   if ((j < nx) && (k < ny)) {
      kk = k + nnyv*j;
      jj = ks + nmxv*js;
      for (i = 0; i < ndim; i++) {
         gn[kk+nyv*i] = shmn2[jj+mxv*i];
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpufft2rcxs(float2 f[], int isign, int mixup[],
                            float2 sct[], int indx, int indy, int nyi,
                            int nyp, int nxhd, int nyd, int nxhyd,
                            int nxyhd, int nsize) {
/* this subroutine performs the x part of a two dimensional real to
   complex fast fourier transform and its inverse, for a subset of y,
   using complex arithmetic, with data not packed
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
   nxhd = first dimension of f >= nx/2+1
   nyd = second dimension of f >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   nsize = amount of scratch complex memory used
   fourier coefficients are stored as follows:
   f[k][j].x, f[k][j].y = real, imaginary part of mode j,k, where
   0 <= j < nx/2+1 and 0 <= k < ny
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, jj, kk;
   int n, nn, in, nt, nh;
   float ani, at1, at2;
   float2 t1, t2, t3;
/* The size of the shared memory array is as follows: */
/* float2 s[nsize];                                   */
   extern __shared__ float2 s[];
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   nxh = nx/2;
   nxhh = nx/4;
   ny = 1L<<indy;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   nyt = nyi + nyp - 1;
/* calculate extent of shared memory usage: */
/* nn = size of shared memory in x          */
   nn = nxh;
   in = 0;
   while (nn > nsize) {
      nn = nn/2;
      in += 1;
   }
/* nt = number of iterations in x           */
   nt = 1L<<in;
   in = indx1 - in;
   nh = nn/2;
/* inverse fourier transform */
   if (isign < 0) {
/* bit-reverse array elements in x */
      nrx = nxhy/nxh;
/* for (k = nyi-1; k < nyt; k++) { */
      k = blockIdx.x + nyi - 1;
      if (k < nyt) {
         jj = nxhd*k;
/*       for (j = 0; j < nxh; j++) { */
         j = threadIdx.x;
         while (j < nxh) {
            j1 = (mixup[j] - 1)/nrx;
            if (j < j1) {
               t1 = f[j1+jj];
               f[j1+jj] = f[j+jj];
               f[j+jj] = t1;
            }
            j += blockDim.x;
         }
/* synchronize threads */
         __syncthreads();
      }
/* copy data to local memory */
      nrx = nxy/nxh;
/*    for (i = nyi-1; i < nyt; i++) { */
      i = blockIdx.x + nyi - 1;
      if (i < nyt) {
         jj = nxhd*i;
         for (n = 0; n < nt; n++) {
/*          for (kk = 0; kk < nn; kk++) { */
            kk = threadIdx.x;
            while (kk < nn) {
               s[kk] = f[kk+nn*n+jj];
               kk += blockDim.x;
            }
/* synchronize threads */
            __syncthreads();
/* transform using local data in x */
            ns = 1;
            for (l = 0; l < in; l++) {
               ns2 = ns + ns;
               km = nxhh/ns;
               kmr = km*nrx;
/*             for (kk = 0; kk < nh; kk++) { */
               kk = threadIdx.x;
               while (kk < nh) {
                  k = kk/ns;
                  j = kk - ns*k;
                  k1 = ns2*k;
                  k2 = k1 + ns;
                  j1 = j + k1;
                  j2 = j + k2;
                  t1 = sct[kmr*j];
                  t2 = s[j2];
                  at1 = t1.x*t2.x - t1.y*t2.y;
                  at2 = t1.x*t2.y + t1.y*t2.x;
                  t2 = s[j1];
                  t3.x = t2.x - at1;
                  t3.y = t2.y - at2;
                  s[j2] = t3;
                  t3.x = t2.x + at1;
                  t3.y = t2.y + at2;
                  s[j1] = t3;
                  kk += blockDim.x;
               }
               ns = ns2;
/* synchronize threads */
               __syncthreads();
            }
/* copy data to global memory */
/*          for (kk = 0; kk < nn; kk++) { */
            kk = threadIdx.x;
            while (kk < nn) {
               f[kk+nn*n+jj] = s[kk];
               kk += blockDim.x;
            }
/* synchronize threads */
            __syncthreads();
         }
/* transform using global data in x */
         ns = 1L<<in;
         for (l = in; l < indx1; l++) {
            ns2 = ns + ns;
            km = nxhh/ns;
            kmr = km*nrx;
/*          for (kk = 0; kk < nxhh; kk++) { */
            kk = threadIdx.x;
            while (kk < nxhh) {
               k = kk/ns;
               j = kk - ns*k;
               k1 = ns2*k;
               k2 = k1 + ns;
               j1 = j + k1;
               j2 = j + k2;
               t1 = sct[kmr*j];
               t2 = f[j2+jj];
               at1 = t1.x*t2.x - t1.y*t2.y;
               at2 = t1.x*t2.y + t1.y*t2.x;
               t2 = f[j1+jj];
               t3.x = t2.x - at1;
               t3.y = t2.y - at2;
               f[j2+jj] = t3;
               t3.x = t2.x + at1;
               t3.y = t2.y + at2;
               f[j1+jj] = t3;
               kk += blockDim.x;
            }
            ns = ns2;
/* synchronize threads */
            __syncthreads();
         }
      }
/* unscramble coefficients and normalize */
      kmr = nxy/nx;
      ani = 0.5f/(((float) nx)*((float) ny));
/* for (k = nyi-1; k < nyt; k++) */
      k = blockIdx.x + nyi - 1;
      if (k < nyt) {
         jj = nxhd*k;
/*       for (j = 1; j < nxhh; j++) { */
         j = threadIdx.x;
         while (j < nxhh) {
            if (j > 0) {
               t3 = sct[kmr*j];
               at1 = t3.y;
               at2 = -t3.x;
               t2 = f[nxh-j+jj];
               t2.y = -t2.y;
               t3 = f[j+jj];
               t1.x = t3.x + t2.x;
               t1.y = t3.y + t2.y;
               t3.x -= t2.x;
               t3.y -= t2.y;
               t2.x = t3.x*at1 - t3.y*at2;
               t2.y = t3.x*at2 + t3.y*at1;
               t3.x = ani*(t1.x + t2.x);
               t3.y = ani*(t1.y + t2.y);
               f[j+jj] = t3;
               t3.x = ani*(t1.x - t2.x);
               t3.y = ani*(t2.y - t1.y);
               f[nxh-j+jj] = t3;
            }
            j += blockDim.x;
         }
         if (threadIdx.x==0) {
            ani = 2.0f*ani;
            t3 = f[nxhh+jj];
            t3.x = ani*t3.x;
            t3.y = -ani*t3.y;
            f[nxhh+jj] = t3;
            t3 = f[jj];
            at1 = t3.x;
            at2 = t3.y;
            t3.x = ani*(at1 - at2);
            t3.y = 0.0f;
            f[nxh+jj] = t3;
            t3.x = ani*(at1 + at2);
            f[jj] = t3;
         }
/* synchronize threads */
         __syncthreads();
      }
   }
/* forward fourier transform */
   if (isign > 0) {
/* scramble coefficients */
      kmr = nxy/nx;
/*    for (k = nyi-1; k < nyt; k++) { */
      k = blockIdx.x + nyi - 1;
      if (k < nyt) {
         jj = nxhd*k;
/*       for (j = 1; j < nxhh; j++) { */
         j = threadIdx.x;
         while (j < nxhh) {
            if (j > 0) {
               t3 = sct[kmr*j];
               at1 = t3.y;
               at2 = t3.x;
               t2 = f[nxh-j+jj];
               t2.y = -t2.y;
               t3 = f[j+jj];
               t1.x = t3.x + t2.x;
               t1.y = t3.y + t2.y;
               t3.x -= t2.x;
               t3.y -= t2.y;
               t2.x = t3.x*at1 - t3.y*at2;
               t2.y = t3.x*at2 + t3.y*at1;
               t3.x = t1.x + t2.x;
               t3.y = t1.y + t2.y;
               f[j+jj] = t3;
               t3.x = t1.x - t2.x;
               t3.y = t2.y - t1.y;
               f[nxh-j+jj] = t3;
            }
            j += blockDim.x;
         }
         if (threadIdx.x==0) {
            t3 = f[nxhh+jj];
            t3.x = 2.0f*t3.x;
            t3.y = -2.0f*t3.y;
            f[nxhh+jj] = t3;
            t3 = f[jj];
            at1 = t3.x;
            t3 = f[nxh+jj];
            at2 = t3.x;
            t3.x = at1 + at2;
            t3.y = at1 - at2;
            f[jj] = t3;
         }
/* synchronize threads */
         __syncthreads();
      }
/* bit-reverse array elements in x */
      nrx = nxhy/nxh;
/* for (k = nyi-1; k < nyt; k++) { */
      k = blockIdx.x + nyi - 1;
      if (k < nyt) {
         jj = nxhd*k;
/*       for (j = 0; j < nxh; j++) { */
         j = threadIdx.x;
         while (j < nxh) {
            j1 = (mixup[j] - 1)/nrx;
            if (j < j1) {
               t1 = f[j1+jj];
               f[j1+jj] = f[j+jj];
               f[j+jj] = t1;
            }
            j += blockDim.x;
         }
/* synchronize threads */
         __syncthreads();
      }
/* copy data to local memory */
      nrx = nxy/nxh;
/*    for (i = nyi-1; i < nyt; i++) { */
      i = blockIdx.x + nyi - 1;
      if (i < nyt) {
         jj = nxhd*i;
         for (n = 0; n < nt; n++) {
/*          for (kk = 0; kk < nn; kk++) { */
            kk = threadIdx.x;
            while (kk < nn) {
               s[kk] = f[kk+nn*n+jj];
               kk += blockDim.x;
            }
/* synchronize threads */
            __syncthreads();
/* transform using local data in x */
            ns = 1;
            for (l = 0; l < in; l++) {
               ns2 = ns + ns;
               km = nxhh/ns;
               kmr = km*nrx;
/*             for (kk = 0; kk < nh; kk++) { */
               kk = threadIdx.x;
               while (kk < nh) {
                  k = kk/ns;
                  j = kk - ns*k;
                  k1 = ns2*k;
                  k2 = k1 + ns;
                  j1 = j + k1;
                  j2 = j + k2;
                  t1 = sct[kmr*j];
                  t1.y = -t1.y;
                  t2 = s[j2];
                  at1 = t1.x*t2.x - t1.y*t2.y;
                  at2 = t1.x*t2.y + t1.y*t2.x;
                  t2 = s[j1];
                  t3.x = t2.x - at1;
                  t3.y = t2.y - at2;
                  s[j2] = t3;
                  t3.x = t2.x + at1;
                  t3.y = t2.y + at2;
                  s[j1] = t3;
                  kk += blockDim.x;
               }
               ns = ns2;
/* synchronize threads */
               __syncthreads();
            }
/* copy data to global memory */
/*          for (kk = 0; kk < nn; kk++) { */
            kk = threadIdx.x;
            while (kk < nn) {
               f[kk+nn*n+jj] = s[kk];
               kk += blockDim.x;
            }
/* synchronize threads */
            __syncthreads();
         }
/* transform using global data in x */
         ns = 1L<<in;
         for (l = in; l < indx1; l++) {
            ns2 = ns + ns;
            km = nxhh/ns;
            kmr = km*nrx;
/*          for (kk = 0; kk < nxhh; kk++) { */
            kk = threadIdx.x;
            while (kk < nxhh) {
               k = kk/ns;
               j = kk - ns*k;
               k1 = ns2*k;
               k2 = k1 + ns;
               j1 = j + k1;
               j2 = j + k2;
               t1 = sct[kmr*j];
               t1.y = -t1.y;
               t2 = f[j2+jj];
               at1 = t1.x*t2.x - t1.y*t2.y;
               at2 = t1.x*t2.y + t1.y*t2.x;
               t2 = f[j1+jj];
               t3.x = t2.x - at1;
               t3.y = t2.y - at2;
               f[j2+jj] = t3;
               t3.x = t2.x + at1;
               t3.y = t2.y + at2;
               f[j1+jj] = t3;
               kk += blockDim.x;
            }
            ns = ns2;
/* synchronize threads */
            __syncthreads();
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpufft2rcys(float2 g[], int isign, int mixup[],
                            float2 sct[], int indx, int indy, int nxi,
                            int nxp, int nxhd, int nyd, int nxhyd,
                            int nxyhd, int nsize) {
/* this subroutine performs the y part of a two dimensional real to
   complex fast fourier transform and its inverse, for a subset of x,
   using complex arithmetic, with data not packed
   for isign = (-1,1), input: all, output: g
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, an inverse fourier transform in y is performed
   g[n][m] = sum(g[j][k]*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, a forward fourier transform in y is performed
   g[j][k] = sum(g[n][m]*exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxi = initial x index used
   nxp = number of x indices used
   nxhd = second dimension of g >= nx/2+1
   nyd = first dimension of g >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   nsize = amount of scratch complex memory used
   fourier coefficients are stored as follows:
   g[j][k] = real, imaginary part of mode j,k, where
   0 <= j < nx/2+1 and 0 <= k < ny
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, koff, kk;
   int n, nn, in, nt, nh;
   float at1, at2;
   float2 t1, t2, t3;
/* The size of the shared memory array is as follows: */
/* float2 s[nsize];                                   */
   extern __shared__ float2 s[];
   indx1 = indx - 1;
   indx1y = indx1 > indy ? indx1 : indy;
   nx = 1L<<indx;
   ny = 1L<<indy;
   nyh = ny/2;
   nxy = nx > ny ? nx : ny;
   nxhy = 1L<<indx1y;
   nxt = nxi + nxp - 1;
/* calculate extent of shared memory usage: */
/* nn = size of shared memory in y          */
   nn = ny;
   in = 0;
   while (nn > nsize) {
      nn = nn/2;
      in += 1;
   }
/* nt = number of iterations in y           */
   nt = 1L<<in;
   in = indy - in;
   nh = nn/2;
/* bit-reverse array elements in y */
   nry = nxhy/ny;
/* for (j = nxi-1; j < nxt; j++) { */
   j = blockIdx.x + nxi - 1;
   if (j < nxt) {
      kk = nyd*j;
/*    for (k = 0; k < ny; k++) { */
      k = threadIdx.x;
      while (k < ny) {
         k1 = (mixup[k] - 1)/nry;
         if (k < k1) {
            t1 = g[k1+kk];
            g[k1+kk] = g[k+kk];
            g[k+kk] = t1;
         }
         k += blockDim.x;
      }
/* synchronize threads */
      __syncthreads();
   }
   nry = nxy/ny;
/* inverse fourier transform in y */
   if (isign < 0) {
/* copy data to local memory */
/*    for (i = nxi-1; i < nxt; i++) { */
      i = blockIdx.x + nxi - 1;
      if (i < nxt) {
         koff = nyd*i;
         for (n = 0; n < nt; n++) {
/*          for (kk = 0; kk < nn; kk++) { */
            kk = threadIdx.x;
            while (kk < nn) {
               s[kk] = g[kk+nn*n+koff];
               kk += blockDim.x;
            }
/* synchronize threads */
            __syncthreads();
/* transform using local data in y */
            ns = 1;
            for (l = 0; l < in; l++) {
               ns2 = ns + ns;
               km = nyh/ns;
               kmr = km*nry;
/*             for (kk = 0; kk < nh; kk++) { */
               kk = threadIdx.x;
               while (kk < nh) {
                  k = kk/ns;
                  j = kk - ns*k;
                  k1 = ns2*k;
                  k2 = k1 + ns;
                  j1 = j + k1;
                  j2 = j + k2;
                  t1 = sct[kmr*j];
                  t2 = s[j2];
                  at1 = t1.x*t2.x - t1.y*t2.y;
                  at2 = t1.x*t2.y + t1.y*t2.x;
                  t3.x = t2.x - at1;
                  t3.y = t2.y - at2;
                  t2 = s[j1];
                  t3.x = t2.x - at1;
                  t3.y = t2.y - at2;
                  s[j2] = t3;
                  t3.x = t2.x + at1;
                  t3.y = t2.y + at2;
                  s[j1] = t3;
                  kk += blockDim.x;
               }
               ns = ns2;
/* synchronize threads */
               __syncthreads();
            }
/* copy data to global memory */
/*          for (kk = 0; kk < nn; kk++) { */
            kk = threadIdx.x;
            while (kk < nn) {
               g[kk+nn*n+koff] = s[kk];
               kk += blockDim.x;
            }
/* synchronize threads */
            __syncthreads();
         }
/* transform using global data in y */
         ns = 1L<<in;
         for (l = in; l < indy; l++) {
            ns2 = ns + ns;
            km = nyh/ns;
            kmr = km*nry;
/*          for (kk = 0; kk < nyh; kk++) { */
            kk = threadIdx.x;
            while (kk < nyh) {
               k = kk/ns;
               j = kk - ns*k;
               k1 = ns2*k;
               k2 = k1 + ns;
               j1 = j + k1;
               j2 = j + k2;
               t1 = sct[kmr*j];
               t2 = g[j2+koff];
               at1 = t1.x*t2.x - t1.y*t2.y;
               at2 = t1.x*t2.y + t1.y*t2.x;
               t3.x = t2.x - at1;
               t3.y = t2.y - at2;
               t2 = g[j1+koff];
               t3.x = t2.x - at1;
               t3.y = t2.y - at2;
               g[j2+koff] = t3;
               t3.x = t2.x + at1;
               t3.y = t2.y + at2;
               g[j1+koff] = t3;
               kk += blockDim.x;
            }
            ns = ns2;
/* synchronize threads */
            __syncthreads();
         }
      }
   }
/* forward fourier transform in y */
   if (isign > 0) {
/* copy data to local memory */
/*    for (i = nxi-1; i < nxt; i++) { */
      i = blockIdx.x + nxi - 1;
      if (i < nxt) {
         koff = nyd*i;
         for (n = 0; n < nt; n++) {
/*          for (kk = 0; kk < nn; kk++) { */
            kk = threadIdx.x;
            while (kk < nn) {
               s[kk] = g[kk+nn*n+koff];
               kk += blockDim.x;
            }
/* synchronize threads */
            __syncthreads();
/* transform using local data in y */
            ns = 1;
            for (l = 0; l < in; l++) {
               ns2 = ns + ns;
               km = nyh/ns;
               kmr = km*nry;
/*             for (kk = 0; kk < nh; kk++) { */
               kk = threadIdx.x;
               while (kk < nh) {
                  k = kk/ns;
                  j = kk - ns*k;
                  k1 = ns2*k;
                  k2 = k1 + ns;
                  j1 = j + k1;
                  j2 = j + k2;
                  t1 = sct[kmr*j];
                  t1.y = -t1.y;
                  t2 = s[j2];
                  at1 = t1.x*t2.x - t1.y*t2.y;
                  at2 = t1.x*t2.y + t1.y*t2.x;
                  t3.x = t2.x - at1;
                  t3.y = t2.y - at2;
                  t2 = s[j1];
                  t3.x = t2.x - at1;
                  t3.y = t2.y - at2;
                  s[j2] = t3;
                  t3.x = t2.x + at1;
                  t3.y = t2.y + at2;
                  s[j1] = t3;
                  kk += blockDim.x;
               }
                ns = ns2;
/* synchronize threads */
               __syncthreads();
            }
/* copy data to global memory */
/*          for (kk = 0; kk < nn; kk++) { */
            kk = threadIdx.x;
            while (kk < nn) {
               g[kk+nn*n+koff] = s[kk];
               kk += blockDim.x;
            }
/* synchronize threads */
            __syncthreads();
         }
/* transform using global data in y */
         ns = 1L<<in;
         for (l = in; l < indy; l++) {
            ns2 = ns + ns;
            km = nyh/ns;
            kmr = km*nry;
/*          for (kk = 0; kk < nyh; kk++) { */
            kk = threadIdx.x;
            while (kk < nyh) {
               k = kk/ns;
               j = kk - ns*k;
               k1 = ns2*k;
               k2 = k1 + ns;
               j1 = j + k1;
               j2 = j + k2;
               t1 = sct[kmr*j];
               t1.y = -t1.y;
               t2 = g[j2+koff];
               at1 = t1.x*t2.x - t1.y*t2.y;
               at2 = t1.x*t2.y + t1.y*t2.x;
               t3.x = t2.x - at1;
               t3.y = t2.y - at2;
               t2 = g[j1+koff];
               t3.x = t2.x - at1;
               t3.y = t2.y - at2;
               g[j2+koff] = t3;
               t3.x = t2.x + at1;
               t3.y = t2.y + at2;
               g[j1+koff] = t3;
               kk += blockDim.x;
            }
            ns = ns2;
/* synchronize threads */
            __syncthreads();
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpusum1(float a[], float *sa, int nx) {
/* 1d serial sum reduction */
/* nx = length of data     */
/* sa = sum(a)             */
/* local data */
   int j, js, jb, mx, joff, mxm;
   float t;
/* The size of the shared memory array is as follows: */
/* ss[blockDim.x];                                    */
   extern __shared__ float ss[];
   mx = blockDim.x;
   js = threadIdx.x;
   jb = blockIdx.x;
   joff = mx*jb;
   j = js + joff;
/* copy global data to shared memory */
   if (j < nx) ss[js] = a[j];
/* synchronize to make sure each thread in block has the data */
   __syncthreads();
   if (js==0) {
      mxm = nx - joff;
      if (mxm > mx) mxm = mx;
/* perform serial local sum reduction: result in t */
      t = 0.0f;
      for (j = 0; j < mxm; j++) {
         t += ss[j];
      }
/* accumulate results to global memory for each block */
/* for devices with compute capability 2.x            */
      atomicAdd(&sa[0],t);
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpusum2(float a[], float d[], int nx) {
/* segmented 1d sum reductions, each of length mx = blockDim.x */
/* nx = length of data                                         */
/* forall (j = 1:nbx); d(j) = sum(a(1+mx*(j-1):min(nx,mx*j)))  */
/* local data */
   int j, js, jb, mx, joff, mxm;
/* The size of the shared memory array is as follows: */
/* ss[blockDim.x];                                    */
   extern __shared__ float ss[];
   mx = blockDim.x;
   js = threadIdx.x;
   jb = blockIdx.x;
   joff = mx*jb;
   j = js + joff;
/* copy global data to shared memory */
   if (j < nx) ss[js] = a[j];
/* synchronize to make sure each thread in block has the data */
   __syncthreads();
   mxm = nx - joff;
   if (mxm > mx) mxm = mx;
/* perform parallel local sum reduction: result in ss[0] */
   lsum2(ss,mxm);
/* write out result to global memory for each block */
   if (js==0) d[jb] = ss[0];
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppush2l(float *ppart, float *fxy, int *kpic,
                            float qbm, float dt, float *ek, int idimp,
                            int nppmx, int nx, int ny, int mx, int my,
                            int nxv, int nyv, int mx1, int mxy1,
                            int ipbc) {
/* Push Interface for C */
   int n, m, ns;
   dim3 dimBlock(nblock_size);
   n = mxy1;
   m = (n - 1)/maxgsx + 1;
   n = n < maxgsx ? n : maxgsx;
   dim3 dimGrid(n,m);
   ns = 2*(mx + 1)*(my + 1)*sizeof(float);
   n = nblock_size*sizeof(float);
   ns = ns > n ? ns : n;
   crc = cudaGetLastError();
   gpuppush2l<<<dimGrid,dimBlock,ns>>>(ppart,fxy,kpic,qbm,dt,ek,idimp,
                                       nppmx,nx,ny,mx,my,nxv,nyv,mx1,
                                       mxy1,ipbc);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppush2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppushf2l(float *ppart, float *fxy, int *kpic,
                             int *ncl, int *ihole, float qbm, float dt,
                             float *ek, int idimp, int nppmx, int nx,
                             int ny, int mx, int my, int nxv, int nyv,
                             int mx1, int mxy1, int ntmax, int *irc) {
/* Push Interface for C */
   int n, m, ns;
   dim3 dimBlock(nblock_size);
   n = mxy1;
   m = (n - 1)/maxgsx + 1;
   n = n < maxgsx ? n : maxgsx;
   dim3 dimGrid(n,m);
   ns = 2*(mx + 1)*(my + 1)*sizeof(float) + nblock_size*sizeof(int);
   n = nblock_size*sizeof(float);
   ns = ns > n ? ns : n;
   ns += 9*sizeof(int);
   crc = cudaGetLastError();
   gpuppushf2l<<<dimGrid,dimBlock,ns>>>(ppart,fxy,kpic,ncl,ihole,qbm,dt,
                                        ek,idimp,nppmx,nx,ny,mx,my,nxv,
                                        nyv,mx1,mxy1,ntmax,irc);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppushf2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpu2ppost2l(float *ppart, float *q, int *kpic,
                             float qm, int nppmx, int idimp, int mx,
                             int my, int nxv, int nyv, int mx1,
                             int mxy1) {
/* Deposit Interface for C */
   int n, m, ns;
   dim3 dimBlock(nblock_size);
   n = mxy1;
   m = (n - 1)/maxgsx + 1;
   n = n < maxgsx ? n : maxgsx;
   dim3 dimGrid(n,m);
   ns = (mx + 1)*(my + 1)*sizeof(float);
   crc = cudaGetLastError();
   gpu2ppost2l<<<dimGrid,dimBlock,ns>>>(ppart,q,kpic,qm,nppmx,idimp,mx,
                                        my,nxv,nyv,mx1,mxy1);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpu2ppost2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpucaguard2l(float2 *qc, float *q, int nx, int ny,
                              int nxe, int nye, int nxvh, int nyv) {
/* Guard Cell Interface for C */
   dim3 dimBlock(nblock_size);
   dim3 dimGrid(ny);
   crc = cudaGetLastError();
   gpucaguard2l<<<dimGrid,dimBlock>>>(qc,q,nx,ny,nxe,nye,nxvh,nyv);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpucaguard2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuccguard2l(float2 *fxyc, float *fxy, int nx, int ny,
                              int nxe, int nye, int nxvh, int nyv) {
/* Guard Cell Interface for C */
   dim3 dimBlock(nblock_size);
   dim3 dimGrid(ny);
   crc = cudaGetLastError();
   gpuccguard2l<<<dimGrid,dimBlock>>>(fxyc,fxy,nx,ny,nxe,nye,nxvh,nyv);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuccguard2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppord2l(float *ppart, float *ppbuff, int *kpic,
                            int *ncl, int *ihole, int idimp, int nppmx,
                            int nx, int ny, int mx, int my, int mx1,
                            int my1, int npbmx, int ntmax, int *irc) {
/* Sort Interface for C */
   int mxy1, n, m, ns;
   dim3 dimBlock(nblock_size);
   mxy1 = mx1*my1;
   m = (mxy1 - 1)/maxgsx + 1;
   n = mxy1 < maxgsx ? mxy1 : maxgsx;
   dim3 dimGrid(n,m);
/* find which particles are leaving tile */
   ns = (nblock_size+9)*sizeof(int);
   crc = cudaGetLastError();
   gpuppfnd2l<<<dimGrid,dimBlock,ns>>>(ppart,kpic,ncl,ihole,idimp,nppmx,
                                       nx,ny,mx,my,mx1,my1,ntmax,irc);
/* cudaThreadSynchronize(); */
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppfnd2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
/* buffer particles that are leaving tile and sum ncl */
   ns = 9*sizeof(int);
   crc = cudaGetLastError();
   gpuppmov2l<<<dimGrid,dimBlock,ns>>>(ppart,ppbuff,ncl,ihole,idimp,
                                       nppmx,mx1,my1,npbmx,ntmax,irc);
/* cudaThreadSynchronize(); */
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppmov2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
/* copy incoming particles from ppbuff into ppart, update kpic */
   ns = (nblock_size+18)*sizeof(int);
   crc = cudaGetLastError();
   gpuppord2l<<<dimGrid,dimBlock,ns>>>(ppart,ppbuff,kpic,ncl,ihole,
                                       idimp,nppmx,mx1,my1,npbmx,ntmax,
                                       irc);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppord2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppordf2l(float *ppart, float *ppbuff, int *kpic,
                             int *ncl, int *ihole, int idimp, int nppmx,
                             int mx1, int my1, int npbmx, int ntmax,
                             int *irc) {
/* Sort Interface for C */
   int mxy1, n, m, ns;
   dim3 dimBlock(nblock_size);
   mxy1 = mx1*my1;
   m = (mxy1 - 1)/maxgsx + 1;
   n = mxy1 < maxgsx ? mxy1 : maxgsx;
   dim3 dimGrid(n,m);
/* buffer particles that are leaving tile and sum ncl */
   ns = 9*sizeof(int);
   crc = cudaGetLastError();
   gpuppmov2l<<<dimGrid,dimBlock,ns>>>(ppart,ppbuff,ncl,ihole,idimp,
                                       nppmx,mx1,my1,npbmx,ntmax,irc);
/* cudaThreadSynchronize(); */
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppmov2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
/* copy incoming particles from ppbuff into ppart, update kpic */
   ns = (nblock_size+18)*sizeof(int);
   crc = cudaGetLastError();
   gpuppord2l<<<dimGrid,dimBlock,ns>>>(ppart,ppbuff,kpic,ncl,ihole,
                                       idimp,nppmx,mx1,my1,npbmx,ntmax,
                                       irc);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppord2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C"  void cgpupois22t(float2 *qt, float2 *fxyt, float2 *ffct,
                             float *we, int nx, int ny, int nxvh,
                             int nyv, int nxhd, int nyhd) {
/* Poisson Solver Interface for C */
   int nxh1, ns;
   dim3 dimBlock(nblock_size);
   nxh1 = nx/2 + 1;
   dim3 dimGrid(nxh1);
   ns = nblock_size*sizeof(float);
   crc = cudaGetLastError();
   gpupois22t<<<dimGrid,dimBlock,ns>>>(qt,fxyt,ffct,we,nx,ny,nxvh,nyv,
                                       nxhd,nyhd);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpupois22t error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwfft2rcs(float2 *f, float2 *g, int isign,
                             int *mixup, float2 *sct, int indx,
                             int indy, int nxhd, int nyd, int nxhyd,
                             int nxyhd) {
/* wrapper function for real to complex fft, without packed data */
/* if isign = -1, f = input, g = output                          */
/* if isign = 1, g = input, f = output                           */
/* nxhd must be >= nx/2 + 1                                      */
/* local data */
   int nxh, nxh1, ny, nsize, ns;
   int nxi = 1, nyi = 1, mx = 16;
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   nxh1 = nxh + 1;
   ny = 1L<<indy;
   dim3 dimGridx(nxh1);
   dim3 dimGridy(ny);
   dim3 dimGridtx((nxh1-1)/mx+1,(ny-1)/mx+1);
   dim3 dimGridty((ny-1)/mx+1,(nxh1-1)/mx+1);
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      nsize = nxh < 1024 ? nxh : 1024;
      ns = nsize*sizeof(float2);
      crc = cudaGetLastError();
      gpufft2rcxs<<<dimGridy,dimBlock,ns>>>(f,isign,mixup,sct,indx,indy,
                                            nyi,ny,nxhd,nyd,nxhyd,nxyhd,
                                            nsize);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpufft2rcxs error=%d:%s\n",
                crc,cudaGetErrorString(crc));
         exit(1);
      }
/* transpose f to g */
      ns = (mx+1)*mx*sizeof(float2);
      crc = cudaGetLastError();
      gpuctpose4<<<dimGridtx,dimBlockt,ns>>>(f,g,nxh1,ny,nxhd,nyd);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuctpose4 error=%d:%s\n",crc,cudaGetErrorString(crc));
         exit(1);
      }
/* perform y fft */
      nsize = ny < 1024 ? ny : 1024;
      ns = nsize*sizeof(float2);
      crc = cudaGetLastError();
      gpufft2rcys<<<dimGridx,dimBlock,ns>>>(g,isign,mixup,sct,indx,indy,
                                            nxi,nxh1,nxhd,nyd,nxhyd,
                                            nxyhd,nsize);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpufft2rcys error=%d:%s\n",
                crc,cudaGetErrorString(crc));
         exit(1);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      nsize = ny < 1024 ? ny : 1024;
      ns = nsize*sizeof(float2);
      crc = cudaGetLastError();
      gpufft2rcys<<<dimGridx,dimBlock,ns>>>(g,isign,mixup,sct,indx,indy,
                                            nxi,nxh1,nxhd,nyd,nxhyd,
                                            nxyhd,nsize);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpufft2rcys error=%d:%s\n",
                crc,cudaGetErrorString(crc));
         exit(1);
      }
/* transpose g to f */
      ns = (mx+1)*mx*sizeof(float2);
      crc = cudaGetLastError();
      gpuctpose4<<<dimGridty,dimBlockt,ns>>>(g,f,ny,nxh1,nyd,nxhd);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuctpose4 error=%d:%s\n",crc,cudaGetErrorString(crc));
         exit(1);
      }
/* perform x fft */
      nsize = nxh < 1024 ? nxh : 1024;
      ns = nsize*sizeof(float2);
      crc = cudaGetLastError();
      gpufft2rcxs<<<dimGridy,dimBlock,ns>>>(f,isign,mixup,sct,indx,indy,
                                            nyi,ny,nxhd,nyd,nxhyd,nxyhd,
                                            nsize);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpufft2rcxs error=%d:%s\n",
                crc,cudaGetErrorString(crc));
         exit(1);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwfft2rcsn(float2 *fn, float2 *gn, int isign,
                              int *mixup, float2 *sct, int indx,
                              int indy, int ndim, int nxhd, int nyd,
                              int nxhyd, int nxyhd) {
/* wrapper function for multiple real to complex ffts, */
/* without packed data                                 */
/* if isign = -1, fn = input, gn = output              */
/* if isign = 1, gn = input, fn = output               */
/* ndim = vector dimension                             */
/* nxhd must be >= nx/2 + 1                            */
/* local data */
   int nxh, nxh1, ny, nxp, nyp, nnxd, nnyd, nsize, ns;
   int nxi = 1, nyi = 1, mx = 16;
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   nxh1 = nxh + 1;
   ny = 1L<<indy;
   nxp = ndim*nxh1;
   nyp = ndim*ny;
   nnxd = ndim*nxhd;
   nnyd = ndim*nyd;
   dim3 dimGridx(nxp);
   dim3 dimGridy(nyp);
   dim3 dimGridtx((nxh1-1)/mx+1,(ny-1)/mx+1);
   dim3 dimGridty((ny-1)/mx+1,(nxh1-1)/mx+1);
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      nsize = nxh < 1024 ? nxh : 1024;
      ns = nsize*sizeof(float2);
      crc = cudaGetLastError();
      gpufft2rcxs<<<dimGridy,dimBlock,ns>>>(fn,isign,mixup,sct,indx,
                                            indy,nyi,nyp,nxhd,nnyd,
                                            nxhyd,nxyhd,nsize);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpufft2rcxs error=%d:%s\n",
                crc,cudaGetErrorString(crc));
         exit(1);
      }
/* transpose f to g */
      ns = ndim*(mx+1)*mx*sizeof(float2);
      crc = cudaGetLastError();
      gpuctpose4n<<<dimGridtx,dimBlockt,ns>>>(fn,gn,nxh1,ny,ndim,nxhd,
                                              nyd);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuctpose4n error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform y fft */
      nsize = ny < 1024 ? ny : 1024;
      ns = nsize*sizeof(float2);
      crc = cudaGetLastError();
      gpufft2rcys<<<dimGridx,dimBlock,ns>>>(gn,isign,mixup,sct,indx,
                                            indy,nxi,nxp,nnxd,nyd,nxhyd,
                                            nxyhd,nsize);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpufft2rcys error=%d:%s\n",
                crc,cudaGetErrorString(crc));
         exit(1);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      nsize = ny < 1024 ? ny : 1024;
      ns = nsize*sizeof(float2);
      crc = cudaGetLastError();
      gpufft2rcys<<<dimGridx,dimBlock,ns>>>(gn,isign,mixup,sct,indx,
                                            indy,nxi,nxp,nnxd,nyd,nxhyd,
                                            nxyhd,nsize);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpufft2rcys error=%d:%s\n",
                crc,cudaGetErrorString(crc));
         exit(1);
      }
/* transpose g to f */
      ns = ndim*(mx+1)*mx*sizeof(float2);
      crc = cudaGetLastError();
      gpuctpose4n<<<dimGridty,dimBlockt,ns>>>(gn,fn,ny,nxh1,ndim,nyd,
                                              nxhd);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuctpose4n error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform x fft */
      nsize = nxh < 1024 ? nxh : 1024;
      ns = nsize*sizeof(float2);
      crc = cudaGetLastError();
      gpufft2rcxs<<<dimGridy,dimBlock,ns>>>(fn,isign,mixup,sct,indx,
                                            indy,nyi,nyp,nxhd,nnyd,
                                            nxhyd,nxyhd,nsize);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpufft2rcxs error=%d:%s\n",
                crc,cudaGetErrorString(crc));
         exit(1);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpusum2(float *a, float *sa, int nx) {
/* segmented 1d parallel sum reduction of input array a, of length nx */
/* first reduce individual blocks in parallel, writing result to scr  */
/* then reduce scr serially, result is written to sa                  */
/* local data */
   int nbx, nbs, ns;
   void *gptr;
   static int len = 0;
   static float *scr = NULL;
   nbx = (nx - 1)/nblock_size + 1;
   dim3 dimBlock(nblock_size);
   dim3 dimGrid(nbx);
   nbs = (nbx - 1)/nblock_size + 1;
   dim3 dimGrid1(nbs);
/* create scratch array */
   if (len < nbx) {
      if (len > 0)
         crc = cudaFree((void *)scr);
      crc = cudaMalloc(&gptr,sizeof(float)*nbx);
      if (crc) {
         printf("cudaMalloc cgpusum2 float Error=%d:%s,l=%d\n",crc,
                 cudaGetErrorString(crc),nbx);
         exit(1);
      }
      scr = (float *)gptr;
      len = nbx;
   }
/* reduce individual blocks in parallel */
   ns = nblock_size*sizeof(float);
   crc = cudaGetLastError();
   gpusum2<<<dimGrid,dimBlock,ns>>>(a,scr,nx);
/* cudaThreadSynchronize(); */
   crc = cudaGetLastError();
   if (crc) {
      printf("gpusum2 error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
/* 1d serial reduction */
   crc = cudaGetLastError();
   gpusum1<<<dimGrid1,dimBlock,ns>>>(scr,sa,nbx);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpusum1 error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
extern "C" void cgpuppush2l_(unsigned long *gp_ppart,
                             unsigned long *gp_fxy,
                             unsigned long *gp_kpic,
                             float *qbm, float *dt, 
                             unsigned long *gp_ek, int *idimp,
                             int *nppmx, int *nx, int *ny, int *mx,
                             int *my, int *nxv, int *nyv, int *mx1,
                             int *mxy1, int *ipbc) {
   float *ppart, *fxy, *ek;
   int *kpic;
   ppart = (float *)*gp_ppart;
   fxy = (float *)*gp_fxy;
   kpic = (int *)*gp_kpic;
   ek = (float *)*gp_ek;
   cgpuppush2l(ppart,fxy,kpic,*qbm,*dt,ek,*idimp,*nppmx,*nx,*ny,*mx,*my,
               *nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppushf2l_(unsigned long *gp_ppart,
                              unsigned long *gp_fxy,
                              unsigned long *gp_kpic,
                              unsigned long *gp_ncl,
                              unsigned long *gp_ihole, float *qbm,
                              float *dt, unsigned long *gp_ek,
                              int *idimp, int *nppmx, int *nx, int *ny,
                              int *mx, int *my, int *nxv, int *nyv,
                              int *mx1, int *mxy1, int *ntmax,
                              unsigned long *gp_irc) {
   float *ppart, *fxy, *ek;
   int *kpic, *ncl, *ihole, *irc;
   ppart = (float *)*gp_ppart;
   fxy = (float *)*gp_fxy;
   kpic = (int *)*gp_kpic;
   ncl = (int *)*gp_ncl;
   ihole = (int *)*gp_ihole;
   ek = (float *)*gp_ek;
   irc = (int *)*gp_irc;
   cgpuppushf2l(ppart,fxy,kpic,ncl,ihole,*qbm,*dt,ek,*idimp,*nppmx,*nx,
                *ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpu2ppost2l_(unsigned long *gp_ppart,
                              unsigned long *gp_q,
                              unsigned long *gp_kpic, float *qm,
                              int *nppmx, int *idimp, int *mx, int *my,
                              int *nxv, int *nyv, int *mx1, int *mxy1) {
   float *ppart, *q;
   int *kpic;
   ppart = (float *)*gp_ppart;
   q = (float *)*gp_q;
   kpic = (int *)*gp_kpic;
   cgpu2ppost2l(ppart,q,kpic,*qm,*nppmx,*idimp,*mx,*my,*nxv,*nyv,*mx1,
                *mxy1);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpucaguard2l_(unsigned long *gp_qc,
                               unsigned long *gp_q, int *nx, int *ny,
                               int *nxe, int *nye, int *nxvh,
                               int *nyv) {
   float2 *qc;
   float *q;
   qc = (float2 *)*gp_qc;
   q = (float *)*gp_q;
   cgpucaguard2l(qc,q,*nx,*ny,*nxe,*nye,*nxvh,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuccguard2l_(unsigned long *gp_fxyc,
                               unsigned long *gp_fxy, int *nx, int *ny,
                               int *nxe, int *nye, int *nxvh,
                               int *nyv) {
   float2 *fxyc;
   float *fxy;
   fxyc = (float2 *)*gp_fxyc;
   fxy = (float *)*gp_fxy;
   cgpuccguard2l(fxyc,fxy,*nx,*ny,*nxe,*nye,*nxvh,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void  cgpuppord2l_(unsigned long *gp_ppart,
                              unsigned long *gp_ppbuff,
                              unsigned long *gp_kpic,
                              unsigned long *gp_ncl,
                              unsigned long *gp_ihole, int *idimp,
                              int *nppmx, int *nx, int *ny, int *mx,
                              int *my, int *mx1, int *my1, int *npbmx,
                              int *ntmax, unsigned long *gp_irc) {
   float *ppart, *ppbuff;
   int *kpic, *ncl, *ihole, *irc;
   ppart = (float *)*gp_ppart;
   ppbuff = (float *)*gp_ppbuff;
   kpic = (int *)*gp_kpic;
   ncl = (int *)*gp_ncl;
   ihole = (int *)*gp_ihole;
   irc = (int *)*gp_irc;
   cgpuppord2l(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*nx,*ny,*mx,
               *my,*mx1,*my1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppordf2l_(unsigned long *gp_ppart,
                              unsigned long *gp_ppbuff,
                              unsigned long *gp_kpic,
                              unsigned long *gp_ncl,
                              unsigned long *gp_ihole, int *idimp,
                              int *nppmx, int *mx1, int *my1, int *npbmx,
                              int *ntmax, unsigned long *gp_irc) {
   float *ppart, *ppbuff;
   int *kpic, *ncl, *ihole, *irc;
   ppart = (float *)*gp_ppart;
   ppbuff = (float *)*gp_ppbuff;
   kpic = (int *)*gp_kpic;
   ncl = (int *)*gp_ncl;
   ihole = (int *)*gp_ihole;
   irc = (int *)*gp_irc;
   cgpuppordf2l(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*mx1,*my1,
                *npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
extern "C"  void cgpupois22t_(unsigned long *gp_qt, 
                              unsigned long *gp_fxyt,
                              unsigned long *gp_ffct,
                              unsigned long *gp_we, int *nx, int *ny,
                              int *nxvh, int *nyv, int *nxhd,
                              int *nyhd) {
   float2 *qt, *fxyt, *ffct;
   float *we;
   qt = (float2 *)*gp_qt;
   fxyt = (float2 *)*gp_fxyt;
   ffct = (float2 *)*gp_ffct;
   we = (float *)*gp_we;
   cgpupois22t(qt,fxyt,ffct,we,*nx,*ny,*nxvh,*nyv,*nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwfft2rcs_(unsigned long *gp_f, unsigned long *gp_g,
                              int *isign, unsigned long *gp_mixup,
                              unsigned long *gp_sct, int *indx,
                              int *indy, int *nxhd, int *nyd,
                              int *nxhyd, int *nxyhd) {
   float2 *f, *g, *sct;
   int *mixup;
   f = (float2 *)*gp_f;
   g = (float2 *)*gp_g;
   mixup = (int *)*gp_mixup;
   sct = (float2 *)*gp_sct;
   cgpuwfft2rcs(f,g,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,
                *nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwfft2rcsn_(unsigned long *gp_fn,
                               unsigned long *gp_gn, int *isign,
                               unsigned long *gp_mixup,
                               unsigned long *gp_sct, int *indx,
                               int *indy, int *ndim, int *nxhd,
                               int *nyd, int *nxhyd, int *nxyhd) {
   float2 *fn, *gn, *sct;
   int *mixup;
   fn = (float2 *)*gp_fn;
   gn = (float2 *)*gp_gn;
   mixup = (int *)*gp_mixup;
   sct = (float2 *)*gp_sct;
   cgpuwfft2rcsn(fn,gn,*isign,mixup,sct,*indx,*indy,*ndim,*nxhd,*nyd,
                 *nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpusum2_(unsigned long *gp_a, unsigned long *gp_sa,
                          int *nx) {
   float *a, *sa;
   a = (float *)*gp_a;
   sa = (float *)*gp_sa;
   cgpusum2(a,sa,*nx);
   return;
}
