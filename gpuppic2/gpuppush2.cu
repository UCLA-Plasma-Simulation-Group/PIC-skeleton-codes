/* CUDA Library for Skeleton 2D Electrostatic GPU-MPI PIC Code */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include "cuda.h"

extern int nblock_size;
extern int maxgsx;

static cudaError_t crc;

extern "C" void gpu_deallocate(void *g_d, int *irc);

extern "C" void gpu_iallocate(int **g_i, int nsize, int *irc);

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
__global__ void gpuppgppush2l(float ppart[], float fxy[], int kpic[],
                              int noff, int nyp, float qbm, float dt,
                              float *ek, int nx, int ny, int mx, int my,
                              int idimp, int nppmx, int nxv, int nypmx,
                              int mx1, int mxyp1, int ipbc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions
   threaded version using guard cells, for distributed data
   data read in tiles
   particles stored segmented array
   42 flops/particle, 12 loads, 4 stores
   input: all except ihole, output: ppart, ek
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
   ppart[m][0][n] = position x of particle n in partition in tile m
   ppart[m][1][n] = position y of particle n in partition in tile m
   ppart[m][2][n] = velocity vx of particle n in partition in tile m
   ppart[m][3][n] = velocity vy of particle n in partition in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,kk)
   fxy[k][j][1] = y component of force/charge at grid (j,kk)
   in other words, fxy are the convolutions of the electric field
   over the particle shape, where kk = k + noff
   kpic = number of particles per tile
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
   qbm = particle charge/mass ratio
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
   int noffp, moffp, npoff, nppp, mxv;
   int mnoff, i, j, k, ii, nn, mm;
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
   edgely = 1.0f;
   edgerx = (float) (nx);
   edgery = (float) (ny-1);
   if ((ipbc==2) || (ipbc==3)) {
      edgelx = 1.0f;
      edgerx = (float) (nx-1);
   }
   mxv = mx + 1;
/* k = tile number */
   k = blockIdx.x + gridDim.x*blockIdx.y;
/* loop over tiles */
   if (k < mxyp1) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      npoff = idimp*nppmx*k;
      mnoff = moffp + noff;
/* load local fields from global array */
      nn = (mx < nx-noffp ? mx : nx-noffp) + 1;
      mm = (my < nyp-moffp ? my : nyp-moffp) + 1;
      ii = threadIdx.x;
      while (ii < mxv*(my+1)) {
         j = ii/mxv;
         i = ii - mxv*j;
         if ((i < nn) && (j < mm)) {
            sfxy[2*ii] = fxy[2*(i+noffp+nxv*(j+moffp))];
            sfxy[1+2*ii] = fxy[1+2*(i+noffp+nxv*(j+moffp))];
         }
         ii += blockDim.x;
      }
/* synchronize threads */
      __syncthreads();
/* loop over particles in tile */
      j = threadIdx.x;
      while (j < nppp) {
/* find interpolation weights */
         x = ppart[j+npoff];
         nn = x;
         y = ppart[j+npoff+nppmx];
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nn = 2*(nn - noffp) + 2*mxv*(mm - mnoff);
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
__global__ void gpu2ppgppost2l(float ppart[], float q[], int kpic[],
                               int noff, float qm, int idimp, int nppmx,
                               int mx, int my, int nxv, int nypmx, 
                               int mx1, int mxyp1) {
/* for 2d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   threaded version using guard cells, for distributed data
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
   int noffp, moffp, npoff, nppp, mxv;
   int mnoff, i, j, k, ii, nn, np, mm, mp;
   float dxp, dyp, amx, amy;
/* The size of the shared memory array is as follows: */
/* float sq[(mx+1)*(my+1)]                            */
   extern __shared__ float sq[];
   mxv = mx + 1;
/* k = tile number */
   k = blockIdx.x + gridDim.x*blockIdx.y;
/* loop over tiles */
   if (k < mxyp1) {
      noffp = k/mx1;
      moffp = my*noffp;
      noffp = mx*(k - mx1*noffp);
      nppp = kpic[k];
      npoff = idimp*nppmx*k;
      mnoff = moffp + noff;
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
      while (j < nppp) {
/* find interpolation weights */
         dxp = ppart[j+npoff];
         nn = dxp;
         dyp = ppart[j+npoff+nppmx];
         mm = dyp;
         dxp = qm*(dxp - (float) nn);
         dyp = dyp - (float) mm;
         nn = nn - noffp;
         mm = mxv*(mm - mnoff);
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
      nn = mxv < nxv-noffp ? mxv : nxv-noffp;
      mm = my+1 < nypmx-moffp ? my+1 : nypmx-moffp;
      ii = threadIdx.x;
      while (ii < mxv*(my+1)) {
         j = ii/mxv;
         i = ii - mxv*j;
         if ((i < nn) && (j < mm)) {
/* original deposit charge, has data hazard on GPU */
/*          q[i+noffp+nxv*(j+moffp)] += sq[ii]; */
/* for devices with compute capability 2.x */
            atomicAdd(&q[i+noffp+nxv*(j+moffp)],sq[ii]);
         }
         ii += blockDim.x;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppcaguard2xl(float2 qc[], float2 scs[], float q[],
                                int nyp, int nx, int nxe, int nypmx,
                                int nxvh, int kypd) {
/* copy and accumulate extended periodic scalar field q in x direction
   into complex output fields qc, scs
   linear interpolation, for distributed data
   scs[j] = data to send to another processor
   nyp = number of primary (complete) gridpoints in particle partition
   nx = system length in x direction
   nxe = first dimension of input field array q, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells
   nxvh = first dimension of output field array qc, must be >= nx/2+1
   kypd = second dimension of output field array qc, must be >= nyp */
/* local data */
   int j, k, nxh;
   float at1;
   float2 a;
   nxh = nx/2;
   k = blockIdx.x;
/* copy interior points */
   if (k < nyp) {
      j = threadIdx.x;
      while (j < nxh) {
         at1 = 0.0f;
         if (j==0) {
            at1 = q[nx+nxe*k];
         }
         a.x = q[2*j+nxe*k] + at1;
         a.y = q[2*j+1+nxe*k];
         qc[j+nxvh*k] = a;
         j += blockDim.x;
      }
   }
/* copy exterior points */
   if (k==0) {
      j = threadIdx.x;
      while (j < nxh) {
         at1 = 0.0f;
         if (j==0) {
            at1 = q[nx+nxe*nyp];
         }
         a.x = q[2*j+nxe*nyp] + at1;
         a.y = q[2*j+1+nxe*nyp];
         scs[j] = a;
         j += blockDim.x;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppcaguard2yl(float2 fc[], float2 scr[], int nx,
                                int nxvh, int kypd) {
/* this subroutine adds data from guard cells from remote processors
   fc[k][j] = complex data for grid j,k in particle partition.
   output: fc
   scr[j] = complex input array for arriving data
   nx = system length in x direction
   nxvh = first dimension of fc, must be >= nx/2+1
   kypd = maximum size of field partition, including guard cells.
   linear interpolation, for distributed data
local data */
   int j, nxh;
   float2 a, b;
   nxh = nx/2;
   j = threadIdx.x+blockDim.x*blockIdx.x;
/* add up the guard cells from remote processors */
   if (j < nxh) {
      a = fc[j];
      b = scr[j];
      a.x += b.x;
      a.y += b.y;
      fc[j] = a;
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppccguard2xl(float2 fxyc[], float2 scs[],
                                float fxy[], int nyp, int nx, int nxe,
                                int nypmx, int nxvh, int kypd) {
/* copy and replicate complex input 2d vector field fxyc in x direction
   into extended periodic fields fxy, scs
   linear interpolation, for distributed data
   scs[2][j] = data to send to another processor
   nyp = number of primary (complete) gridpoints in particle partition
   nx = system length in x direction
   nxe = second dimension of input field array fxy, must be >= nx+1
   nypmx = maximum size of particle partition, including guard cells
   nxvh = first dimension of input field array fxyc, must be >= nx/2+1
   kypd = third dimension of input field array fxyc, must be >= nyp+1 */
/* local data */
   int j, k, nxh;
   float2 a, b;
   nxh = nx/2;
   k = blockIdx.x;
/* copy interior points */
   if (k < nyp) {
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
/* copy edges of extended field */
   if (blockIdx.x==0) {
      k = threadIdx.x;
      while (k < nyp) {
         a = fxyc[nxvh*2*k];
         b = fxyc[nxvh*(1+2*k)];
         fxy[2*(nx+nxe*k)] = a.x;
         fxy[1+2*(nx+nxe*k)] = b.x;
         k += blockDim.x;
      }
/* copy exterior points */
      j = threadIdx.x;
      while (j < nxh) {
         scs[j] = fxyc[j];
         scs[j+nxvh] = fxyc[j+nxvh];
         j += blockDim.x;
      }
/* copy edges of extended field */
      if (threadIdx.x==0) {
         a = fxyc[0];
         b = fxyc[nxvh];
         a.y = 0.0f;
         b.y = 0.0f;
         scs[nxh] = a;
         scs[nxh+nxvh] = b;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppccguard2yl(float fxy[], float2 scr[], int nyp,
                                int nx, int nxe, int nxvh, int nypmx) {
/* this subroutine copies data to guard cells from remote processors
   fxy[k][j] = real data for grid j,k in particle partition.
   the grid is non-uniform and includes one extra guard cell.
   output: fxy
   scr[2][j] = complex input array for arriving data
   nyp = number of primary gridpoints in field partition
   it is assumed the nyp > 0.
   nx = system length in x direction
   nxe = second dimension of input field array fxy, must be >= nx+1
   nxvh = first dimension of scr, must be >= nx/2+1
   nypmx = maximum size of field partition, including guard cell.
   linear interpolation, for distributed data
local data */
   int j, nxh;
   float2 a, b;
   nxh = nx/2;
   j = threadIdx.x+blockDim.x*blockIdx.x;
/* copy to guard cells */
   if (j < nxh) {
      a = scr[j];
      b = scr[j+nxvh];
      fxy[2*(2*j+nxe*nyp)] = a.x;
      fxy[1+2*(2*j+nxe*nyp)] = b.x;
      fxy[2*(2*j+1+nxe*nyp)] = a.y;
      fxy[1+2*(2*j+1+nxe*nyp)] = b.y;
   }
   if (j==0) {
      a = scr[nxh];
      b = scr[nxh+nxvh];
      fxy[2*(nx+nxe*nyp)] = a.x;
      fxy[1+2*(nx+nxe*nyp)] = b.x;
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpupppfnd2l(float ppart[], int kpic[], int ncl[],
                            int ihole[], int noff, int nyp, int idimp,
                            int nppmx, int nx, int ny, int mx, int my,
                            int mx1, int myp1, int ntmax, int *irc) {
/* this subroutine performs first step of a particle sort by x,y grid
   in tiles of mx, my, where one finds the particles leaving tile and
   stores their number, location, and destination in ncl and ihole.
   linear interpolation, with periodic boundary conditions
   for distributed data, with 1d domain decomposition in y.
   tiles are assumed to be arranged in 2D linear memory
   input: all except ncl, ihole, irc
   output: ncl, ihole, irc
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k 
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   noff = lowermost global gridpoint in particle partition.
   nyp = number of primary (complete) gridpoints in particle partition
   idimp = size of phase space = 4
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   mx1 = (system length in x direction - 1)/mx + 1
   myp1 = (partition length in y direction - 1)/my + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int mxyp1, noffp, moffp, nppp, j, k, ih, ist, nn, mm, nths;
   float anx, any, edgelx, edgely, edgerx, edgery, dx, dy;
/* The sizes of the shared memory arrays are as follows: */
/* int sncl[8], sih[blockDim.x], nh[1];                  */
   int *sncl, *sih, *nh;
   extern __shared__ int shm[];
   sncl = (int *)&shm[0];
   sih = (int *)&shm[8];
   nh = (int *)&shm[8+blockDim.x];
   mxyp1 = mx1*myp1;
   anx = (float) nx;
   any = (float) ny;
/* k = tile number */
   k = blockIdx.x + gridDim.x*blockIdx.y;
/* find and count particles leaving tiles and determine destination */
/* update ppart, ihole, ncl */
/* loop over tiles */
   if (k < mxyp1) {
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
      mm = (nppp - 1)/(int) blockDim.x + 1;
      noffp = 0;
      for (nn = 0; nn < mm; nn++) {
         j = threadIdx.x + blockDim.x*nn;
         sih[threadIdx.x] = 0;
         if (j < nppp) {
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
         nths = nppp - blockDim.x*nn;
         if (nths > blockDim.x)
            nths = blockDim.x;
/* perform local prefix reduction */
         liscan2(sih,nths);
         if (j < nppp) {
            ih = sih[threadIdx.x];
            moffp = 0;
            if (threadIdx.x > 0)
               moffp = sih[threadIdx.x-1];
/* this thread has a hole present */
            if (ih > moffp) {
               ih += noffp;
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
            noffp += sih[nths-1];
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
         ih  = noffp;
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
__global__ void gpupppmov2l(float ppart[], float ppbuff[], int ncl[],
                            int ihole[], int idimp, int nppmx, int mx1,
                            int myp1, int npbmx, int ntmax, int *irc) {
/* this subroutine performs second step of a particle sort by x,y grid
   in tiles of mx, my, where prefix scan of ncl is performed and
   departing particles are buffered in ppbuff in direction order.
   linear interpolation, with periodic boundary conditions
   for distributed data, with 1d domain decomposition in y.
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
   myp1 = (partition length in y direction - 1)/my + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int mxyp1, i, j, k, ii, nh, ist, j1, ierr;
/* The sizes of the shared memory arrays are as follows: */
/* int sncl[8], ip[1];                                   */
/* blockDim.x should be >= 8                             */
   int *sncl, *ip;
   extern __shared__ int shm[];
   sncl = (int *)&shm[0];
   ip = (int *)&shm[8];
   mxyp1 = mx1*myp1;
   ierr = 0;
/* k = tile number */
   k = blockIdx.x + gridDim.x*blockIdx.y;
   j = threadIdx.x;
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
   if (k < mxyp1) {
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
__global__ void nciscan2(int ncl[], int ncll[], int nclr[], int nscr[], 
                         int mx1, int myp1) {
/* calculate number number offsets for particles leaving processor
   by performing prefix scan (running sum) for each block.  the last
   value of each block is written to nscr to add to each block later.
   ncl[k][i] = number of particles going to destination i, tile k
   ncll = number offset being sent to lower processor
   nclr = number offset being sent to upper processor
   nscr = scratch integer array, of size 2*((mx1-1)/blockDim.x+1)
   mx1 = (system length in x direction - 1)/mx + 1
   myp1 = (partition length in y direction - 1)/my + 1
local data                                                            */
   int l, koff, ii, jj, kk, nths;
/* The size of the shared memory array are as follows: */
/* int isdata[2*blockDim.x];                           */
   int *isdata;
   extern __shared__ int shm[];
   isdata = (int *)&shm[0];
   l = threadIdx.x;
   koff = blockDim.x*blockIdx.x;
   kk = mx1*(myp1 - 1) + koff;
   nths = mx1 - koff;
   if (nths > blockDim.x)
      nths = blockDim.x;
   if ((l+koff) < mx1) {
      ii = ncl[4+8*(l+koff)] - ncl[1+8*(l+koff)];
      isdata[l] = ii;
      jj = ncl[7+8*(l+kk)] - ncl[4+8*(l+kk)];
      isdata[l+blockDim.x] = jj;
   }
/* synchronize threads */
   __syncthreads();
/* perform local prefix reductions */
   liscan2(isdata,nths);
   liscan2(&isdata[blockDim.x],nths);
   if ((l+koff) < mx1) {
      ncll[3*(l+koff)] = ii;
      ncll[3*(l+koff)+1] = isdata[l] - ii;
      nclr[3*(l+koff)] = jj;
      nclr[3*(l+koff)+1] = isdata[l+blockDim.x] - jj;
   }
   if (l==0) {
      ii = 0;
      jj = 0;
      if (nths > 0) {
         ii = isdata[nths-1];
         jj = isdata[nths+blockDim.x-1];
      }
      nscr[blockIdx.x] = ii;
      nscr[blockIdx.x+gridDim.x] = jj;
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpupppbuf2l(float ppbuff[], float sbufl[],
                            float sbufr[], int ncl[], int ncll[],
                            int nclr[], int nscr[], int idimp, int mx1,
                            int myp1, int npbmx, int nbmax, int *irc) {
/* this subroutine performs third step of a particle sort by x,y grid
   in tiles of mx, my, where particles leaving the processor are
   buffered in sbufl and sbufr, and particle number offsets are stored
   in ncll and nclr.
   gpupppbuf2l and nciscan2 should use the same blocksize.
   linear interpolation, with periodic boundary conditions
   for distributed data, with 1d domain decomposition in y.
   tiles are assumed to be arranged in 2D linear memory
   input: all except sbufl, sbufr, ncll, nclr, irc
   output: sbufl, sbufr, ncll, nclr, irc
   ppbuff[k][i][n] = i co-ordinate of particle n in tile k
   sbufl = buffer for particles being sent to lower processor
   sbufr = buffer for particles being sent to upper processor
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ncll = number offset being sent to lower processor
   nclr = number offset being sent to upper processor
   nscr = scratch integer array, of size 2*((mx1-1)/blockDim.x+1)
   idimp = size of phase space = 4
   mx1 = (system length in x direction - 1)/mx + 1
   myp1 = (partition length in y direction - 1)/my + 1
   npbmx = size of buffer array ppbuff
   nbmax =  size of buffers for passing particles between processors
   irc = maximum overflow, returned only if error occurs, when irc > 0
local data                                                            */
   int i, j, k, ii, jj, nl, nr, nn, mm, nbl, kk, ll, im;
   k = blockIdx.x;
/* buffer particles and their number leaving the node: */
/* update sbufl, sbufr, ncll, nclr */
   nbl = (mx1 - 1)/blockDim.x + 1;
   nl = 0;
   nr = 0;
   nn = 0;
   mm = 0;
   kk = mx1*(myp1 - 1);
/* loop over row of tiles */
   if (k < mx1) {
      j = k/blockDim.x + 1;
/* find how many particles must be buffered */
      for (i = 1; i < nbl; i++) {
         nl += nscr[i-1];
         nr += nscr[i+nbl-1];
/* save offsets */
         if (i==(j-1)) {
            nn = nl;
            mm = nr;
         }
      }
      nl += nscr[nbl-1];
      nr += nscr[2*nbl-1];
      ii = ncll[3*k];
      nn += ncll[3*k+1];
      im = nclr[3*k];
      mm += nclr[3*k+1];
/* synchronize threads */
      __syncthreads();
      ll = ncl[1+8*k];
      jj = nl - nn;
      jj = ii < jj ? ii : jj;
      j = threadIdx.x;
      while (j < jj) {
         for (i = 0; i < idimp; i++) {
            sbufl[j+nn+nl*i]
            = ppbuff[j+ll+npbmx*(i+idimp*k)];

         }
         j += blockDim.x;
      }
      ll = nn - ll;
      if (threadIdx.x < 3) {
         ncll[threadIdx.x+3*k] = ncl[threadIdx.x+2+8*k] + ll;
      }
      nn += ii;
      ii = im;
      ll = ncl[4+8*(k+kk)];
      jj = nr - mm;
      jj = ii < jj ? ii : jj;
      j = threadIdx.x;
      while (j < jj) {
         for (i = 0; i < idimp; i++) {
            sbufr[j+mm+nr*i]
            = ppbuff[j+ll+npbmx*(i+idimp*(k+kk))];
         }
         j += blockDim.x;
      }
      ll = mm - ll;
      if (threadIdx.x < 3) {
         nclr[threadIdx.x+3*k] = ncl[threadIdx.x+5+8*(k+kk)] + ll;
      }
      mm += ii;
   }
/* sbufl or sbufr overflow */
   ii = nn > mm ? nn : mm;
   if (ii > nbmax)
      *irc = ii;
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpupppord2l(float ppart[], float ppbuff[],
                            float rbufl[], float rbufr[], int kpic[],
                            int ncl[], int ihole[], int mcll[],
                            int mclr[], int idimp, int nppmx, int mx1,
                            int myp1, int npbmx, int ntmax, int nbmax,
                            int *irc) {
/* this subroutine performs third step of a particle sort by x,y grid
   in tiles of mx, my, where incoming particles from other tiles are
   copied into ppart from ppbuff, rbufl, and rbufr
   linear interpolation, with periodic boundary conditions
   for distributed data, with 1d domain decomposition in y.
   tiles are assumed to be arranged in 2D linear memory
   input: all except irc
   output: ppart, kpic, irc
   ppart[k][i][n] = i co-ordinate of particle n in tile k 
   ppbuff[k][i][n] = i co-ordinate of particle n in tile k
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
   int mxyp1, nppp, ncoff, noff, moff, i, j, k, ii, jj, kx, ky, ni, nh;
   int nn, mm, ll, ip, j1, j2, kxl, kxr, kk, kl, kr, nr, mr;
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
   mxyp1 = mx1*myp1;
   noff = 0;
   moff = 0;
/* k = tile number */
   k = blockIdx.x + gridDim.x*blockIdx.y;
/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
   if (k < mxyp1) {
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
      if (ky==0) {
         nr = mcll[3*mx1-1];
         if (kx > 0)
            noff = mcll[3*kx-1];
      }
      if (ky==(myp1-1)) {
         mr = mclr[3*mx1-1];
         if (kx > 0)
            moff = mclr[3*kx-1];
      }
      ncoff = 0;
      ip = 0;
      ii = threadIdx.x;
      if (ii < 8) {
         kk = ks[ii];
/* ip = number of particles coming from direction ii */
         if (kk < 0) {
            if (ii > 5)
               noff = mcll[ii-6+3*(kk+mx1)];
            ip = mcll[ii-5+3*(kk+mx1)] - noff;
            kk = noff;
         }
         else if (kk >= mxyp1) {
            if (ii > 2)
               moff = mclr[ii-3+3*(kk-mxyp1)];
            ip = mclr[ii-2+3*(kk-mxyp1)] - moff;
            kk = moff;
         }
         else {
            if (ii > 0)
               ncoff = ncl[ii-1+8*kk];
            ip = ncl[ii+8*kk] - ncoff;
            kk = ncoff + idimp*npbmx*kk;
         }
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
               j1 = nppp + (j - nh);
            }
            if (j1 < nppmx) {
               jj = sj[threadIdx.x];
               if ((ky==0) && (j>=sip[4]) && (j<sip[7])) {
                  for (i = 0; i < idimp; i++) {
                     ppart[j1+nppmx*(i+idimp*k)]
                     = rbufl[j+jj+nr*i];
                  }
               }
               else if ((ky==(myp1-1)) && (j>=sip[1]) && (j<sip[4])) {
                  for (i = 0; i < idimp; i++) {
                     ppart[j1+nppmx*(i+idimp*k)]
                     = rbufr[j+jj+mr*i];
                  }
               }
               else {
                  for (i = 0; i < idimp; i++) {
                     ppart[j1+nppmx*(i+idimp*k)]
                     = ppbuff[j+jj+npbmx*i];
                  }
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
         nppp += jj;
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
               j1 = nppp - j - 1;
            }
/* j2 = locations of holes at the end, in decreasing order */
            j2 = 0;
            jj = nh - ll - threadIdx.x;
            if (jj > 0) {
               j2 = ihole[2*(jj+(ntmax+1)*k)] - 1;
            }
/* holes with locations greater than nppp-ip do not need to be filled */
/* identify such holes */
            sj[threadIdx.x] = 1;
/* synchronize threads */
            __syncthreads();
/* omit particles at end that are holes */
            ii = nppp - (j2 + blockDim.x*nn) - 1;
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
               j1 = nppp - j - 1;
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
         nppp -= ip;
      }
/* set error and update particle */
      if (threadIdx.x==0) {
/* ppart overflow */
         if (ist[0] > 0)
            *irc = nppp;
         kpic[k] = nppp;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppois22t(float2 qt[], float2 fxyt[], float2 ffct[],
                            float *we, int nx, int ny, int kstrt,
                            int nyv, int kxp1, int nyhd) {
/* this subroutine solves 2d poisson's equation in fourier space for
   force/charge (or convolution of electric field over particle shape)
   with periodic boundary conditions, for distributed data.
   vector length is second dimension.
   input: qt,ffct,nx,ny,nyv,kxp1,nyhd, output: fxyt,we
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
   qt[j][k] = complex charge density for fourier mode (jj,k)
   fxyt[j][0][k] = x component of complex force/charge,
   fxyt[j][1][k] = y component of complex force/charge,
   for fourier mode (jj,k), where jj = j + kxp1*(kstrt - 1)
   kxp1 = number of data values per block for unpacked field data
   kstrt = starting data block number
   aimag(ffct[j][k]) = finite-size particle shape factor s
   real(ffct[j][k])) = potential green's function g
   for fourier mode (jj,k), where jj = j + kxp1*(kstrt - 1)
   electric field energy is also calculated, using
   we = nx*ny*sum((affp/(kx**2+ky**2))*|q[ky][kx]*s[ky][kx]|**2)
   where affp = nx*ny/np, where np=number of particles
   nx/ny = system length in x/y direction
   nyv = first dimension of field arrays, must be >= ny
   nyhd = first dimension of form factor array, must be >= nyh
local data                                                 */
   int nxh, nyh, ks, joff, kxps, j, jj, jk, k, j0, j1, k1, jk2;
   float dnx, dny, dkx, at1, at2, at3, at4;
   float2 zero, zt1, zt2, zt3;
/* The size of the shared memory array is as follows: */
/* float ss[blockDim.x];                              */
   extern __shared__ float ss[];
   double wp;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   ks = kstrt - 1;
   joff = kxp1*ks;
   j1 = nxh + 1;
   kxps = j1 - joff;
   kxps = 0 > kxps ? 0 : kxps;
   kxps = kxp1 < kxps ? kxp1 : kxps;
   dnx = 6.28318530717959f/(float) nx;
   dny = 6.28318530717959f/(float) ny;
   zero.x = 0.0f;
   zero.y = 0.0f;
/* calculate force/charge and sum field energy */
   wp = 0.0;
   if (kstrt <= j1) {
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
/*    for (j = 0; j < kxps; j++) { */
      j = blockIdx.x;
      j0 = j + joff;
      dkx = dnx*(float) j0;
      jj = nyhd*j;
      jk = nyv*j;
      jk2 = 2*jk;
      if ((j0 > 0) && (j0 < nxh)) {
/*       for (k = 1; k < nyh; k++) { */
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
               zt3.x = at2*zt2.x;
               zt3.y = at2*zt2.y;
               fxyt[k1+jk2] = zt3;
               zt3.x = at3*zt1.x;
               zt3.y = at3*zt1.y;
               fxyt[k+nyv+jk2] = zt3;
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
/*       for (j = 1; j < kxps; j++) { */
         j = threadIdx.x;
         while (j < kxps) {
            j0 = j + joff;
            dkx = dnx*(float) j0;
            jj = nyhd*j;
            jk = nyv*j;
            jk2 = 2*jk;
            if ((j0 > 0) && (j0 < nxh)) {
               zt1 = ffct[jj];
               at1 = zt1.x*zt1.y;
               at2 = dkx*at1;
               zt1 = qt[jk];
               at4 = zt1.x;
               zt3.x = at2*zt1.y;
               zt3.y = -at2*at4;
               fxyt[jk2] = zt3;
               fxyt[k1+jk2] = zero;
               fxyt[nyv+jk2] = zero;
               fxyt[k1+nyv+jk2] = zero;
               wp += (double) (at1*(zt1.x*zt1.x + zt1.y*zt1.y));
            }
            j += blockDim.x;
         }
/* mode numbers kx = 0 */
         if (ks==0) {
/*          for (k = 1; k < nyh; k++) { */
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
                  fxyt[k1] = zero;
                  fxyt[k+nyv] = zt3;
                  zt3.y = -zt3.y;
                  fxyt[k1+nyv] = zt3;
                  wp += (double) (at1*(zt1.x*zt1.x + zt1.y*zt1.y));
               }
               k += blockDim.x;
            }
            if (threadIdx.x==0) {
               k1 = nyh;
               fxyt[0] = zero;
               fxyt[k1] = zero;
               fxyt[nyv] = zero;
               fxyt[k1+nyv] = zero;
            }
         }
/* mode numbers kx = nx/2 */
         if (ks==(nxh/kxp1)) {
            jk = 2*nyv*(kxps-1);
/*          for (k = 0; k < ny; k++) { */
            k = threadIdx.x;
            while (k < ny) {
               fxyt[k+jk] = zero;
               fxyt[k+nyv+jk] = zero;
               k += blockDim.x;
            }
         }
      }
   }
   j = blockIdx.x;
   if (j < kxps) {
/* sum potential energies for each x co-ordinate */
      ss[threadIdx.x] = (float) wp;
/* synchronize threads */
      __syncthreads();
      lsum2(ss,blockDim.x);
/* normalize potential energy for each x co-ordinate */
      if (threadIdx.x==0)
         we[j] = ss[0]*((float) nx)*((float) ny);
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
   using complex arithmetic.
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
__global__ void gpuppmtposes(float2 f[], float2 sm[], int nx, int kxp,
                             int kyps, int kstrt, int nvp, int kxyp,
                             int nxv, int kypd) {
/* extract data to send */
/* local data */
   int ks, j, k, n, nn, id, joff, ld;
   ks = kstrt - 1;
/* for (n = 0; n < nvp; n++) { */
   n = blockIdx.y;
   if (n < nvp) {
      id = n - ks;
      if (id < 0)
         id += nvp;
/* find which node sends to itself */
      nn = 2*ks;
      if (nn >= nvp)
         nn -= nvp;
/* adjust counter */
      if (n > nn)
         n -= 1;
/* do not send local data */
      if (id != ks) {
         joff = kxp*id;
         ld = nx - joff;
         ld = 0 > ld ? 0 : ld;
         ld = kxp < ld ? kxp : ld;
/*       for (k = 0; k < kyps; k++) { */
         k = blockIdx.x;
         if (k < kyps) {
/*          for (j = 0; j < ld; j++) { */
            j = threadIdx.x;
            while (j < ld) {
               sm[j+ld*k+kxyp*n] = f[j+joff+nxv*k];
               j += blockDim.x;
            }
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppmtposer(float2 g[], float2 tm[], int ny, int kyp,
                             int kxps, int kstrt, int nvp, int kxyp,
                             int nyv, int kxpd) {
/* transpose data received */
/* local data */
   int kt, mxv, j, k, n, nn, id, koff, ld, js, ks, jj, kk;
/* The size of the shared memory array is as follows: */
/* float2 s2[(mx + 1)*mx];                            */
   extern __shared__ float2 s2[];
   kt = kstrt - 1;
   mxv = blockDim.x + 1;
/* for (n = 0; n < nvp; n++) { */
   n = blockIdx.z;
   if (n < nvp) {
      id = n - kt;
      if (id < 0)
         id += nvp;
/* find which node sends to itself */
      nn = 2*kt;
      if (nn >= nvp)
         nn -= nvp;
/* adjust counter */
      if (n > nn)
         n -= 1;
/* do not transpose local data */
      if (id != kt) {
         koff = kyp*id;
         ld = ny - koff;
         ld = 0 > ld ? 0 : ld;
         ld = kyp < ld ? kyp : ld;
         js = threadIdx.x;
         ks = threadIdx.y;
         jj = blockDim.x*blockIdx.x;
         kk = blockDim.y*blockIdx.y;
         j = js + jj;
         k = ks + kk;
         if ((j < kxps) && (k < ld)) {
            s2[js+mxv*ks] = tm[j+kxps*k+kxyp*n];
         }
/* synchronize threads */
         __syncthreads();
         j = ks + jj;
         k = js + kk;
         if ((j < kxps) && (k < ld)) {
            g[k+koff+nyv*j] = s2[ks+mxv*js];
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppmtposesn(float2 fn[], float2 sm[], int nx, int kxp,
                              int kyps, int kstrt, int nvp, int ndim,
                              int kxyp, int nxv, int kypd) {
/* extract vector data to send */
/* local data */
   int ks, i, j, k, n, nn, id, joff, ld, nnxv, nkxyp;
   ks = kstrt - 1;
   nnxv = ndim*nxv;
   nkxyp = ndim*kxyp;
/* for (n = 0; n < nvp; n++) { */
   n = blockIdx.y;
   if (n < nvp) {
      id = n - ks;
      if (id < 0)
         id += nvp;
/* find which node sends to itself */
      nn = 2*ks;
      if (nn >= nvp)
         nn -= nvp;
/* adjust counter */
      if (n > nn)
         n -= 1;
/* do not send local data */
      if (id != ks) {
         joff = kxp*id;
         ld = nx - joff;
         ld = 0 > ld ? 0 : ld;
         ld = kxp < ld ? kxp : ld;
/*       for (k = 0; k < kyps; k++) { */
         k = blockIdx.x;
         if (k < kyps) {
/*          for (j = 0; j < ld; j++) { */
            j = threadIdx.x;
            while (j < ld) {
               for (i = 0; i < ndim; i++) {
                  sm[j+ld*(i+ndim*k)+nkxyp*n] = fn[j+joff+nxv*i+nnxv*k];
               }
               j += blockDim.x;
            }
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppmtposern(float2 gn[], float2 tm[], int ny, int kyp,
                              int kxps, int kstrt, int nvp, int ndim,
                              int kxyp, int nyv, int kxpd) {
/* transpose vector data received */
/* local data */
   int kt, mxv, i, j, k, n, nn, id, koff, ld, js, ks, jj, kk;
   int nnyv, nkxyp;
/* The size of the shared memory array is as follows: */
/* float2 s2n[ndim*(mx + 1)*mx];                      */
   extern __shared__ float2 s2n[];
   kt = kstrt - 1;
   mxv = blockDim.x + 1;
   nnyv = ndim*nyv;
   nkxyp = ndim*kxyp;
/* for (n = 0; n < nvp; n++) { */
   n = blockIdx.z;
   if (n < nvp) {
      id = n - kt;
      if (id < 0)
         id += nvp;
/* find which node sends to itself */
      nn = 2*kt;
      if (nn >= nvp)
         nn -= nvp;
/* adjust counter */
      if (n > nn)
         n -= 1;
/* do not transpose local data */
      if (id != kt) {
         koff = kyp*id;
         ld = ny - koff;
         ld = 0 > ld ? 0 : ld;
         ld = kyp < ld ? kyp : ld;
         js = threadIdx.x;
         ks = threadIdx.y;
         jj = blockDim.x*blockIdx.x;
         kk = blockDim.y*blockIdx.y;
         j = js + jj;
         k = ks + kk;
         if ((j < kxps) && (k < ld)) {
            for (i = 0; i < ndim; i++) {
               s2n[js+mxv*(i+ndim*ks)] = tm[j+kxps*(i+ndim*k)+nkxyp*n];
            }
         }
/* synchronize threads */
         __syncthreads();
         j = ks + jj;
         k = js + kk;
         if ((j < kxps) && (k < ld)) {
            for (i = 0; i < ndim; i++) {
               gn[k+koff+nyv*i+nnyv*j] = s2n[ks+mxv*(i+ndim*js)];
            }
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppltpose(float2 f[], float2 g[], int nx, int ny,
                            int kxp, int kyp, int kstrt, int nxv,
                            int nyv) {
/* transpose local data */
/* local data */
   int mxv, j, k, ks, kxps, kyps, joff, koff, js, jj, kk;
/* The size of the shared memory array is as follows: */
/* float2 s2[(mx + 1)*mx];                            */
   extern __shared__ float2 s2[];
   mxv = blockDim.x + 1;
   ks = kstrt - 1;
   joff = kxp*ks;
   koff = kyp*ks;
   kxps = nx - joff;
   kxps = 0 > kxps ? 0 : kxps;
   kxps = kxp < kxps ? kxp : kxps;
   kyps = ny - koff;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = kyp < kyps ? kyp : kyps;
   js = threadIdx.x;
   ks = threadIdx.y;
   jj = blockDim.x*blockIdx.x;
   kk = blockDim.y*blockIdx.y;
   j = js + jj;
   k = ks + kk;
   if ((j < kxps) && (k < kyps)) {
      s2[js+mxv*ks] = f[j+joff+nxv*k];
   }
/* synchronize threads */
   __syncthreads();
   j = ks + jj;
   k = js + kk;
   if ((j < kxps) && (k < kyps)) {
      g[k+koff+nyv*j] = s2[ks+mxv*js];
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppltposen(float2 fn[], float2 gn[], int nx, int ny,
                             int kxp, int kyp, int kstrt, int ndim,
                             int nxv, int nyv) {
/* transpose local vector data */
/* local data */
   int mxv, i, j, k, ks, kxps, kyps, joff, koff, js, jj, kk;
   int nnxv, nnyv;
/* The size of the shared memory array is as follows: */
/* float2 s2n[ndim*(mx + 1)*mx];                      */
   extern __shared__ float2 s2n[];
   mxv = blockDim.x + 1;
   ks = kstrt - 1;
   nnxv = ndim*nxv;
   nnyv = ndim*nyv;
   joff = kxp*ks;
   koff = kyp*ks;
   kxps = nx - joff;
   kxps = 0 > kxps ? 0 : kxps;
   kxps = kxp < kxps ? kxp : kxps;
   kyps = ny - koff;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = kyp < kyps ? kyp : kyps;
   js = threadIdx.x;
   ks = threadIdx.y;
   jj = blockDim.x*blockIdx.x;
   kk = blockDim.y*blockIdx.y;
   j = js + jj;
   k = ks + kk;
   if ((j < kxps) && (k < kyps)) {
      for (i = 0; i < ndim; i++) {
         s2n[js+mxv*(i+ndim*ks)] = fn[j+joff+nxv*i+nnxv*k];
      }
   }
/* synchronize threads */
   __syncthreads();
   j = ks + jj;
   k = js + kk;
   if ((j < kxps) && (k < kyps)) {
      for (i = 0; i < ndim; i++) {
         gn[k+koff+nyv*i+nnyv*j] = s2n[ks+mxv*(i+ndim*js)];
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
extern "C" void cgpuppgppush2l(float *ppart, float *fxy, int *kpic,
                               int noff, int nyp, float qbm, float dt,
                               float *ek, int nx, int ny, int mx,
                               int my, int idimp, int nppmx, int nxv,
                               int nypmx, int mx1, int mxyp1,
                               int ipbc) {
/* Push Interface for C */
   int n, m, ns;
   dim3 dimBlock(nblock_size);
   n = mxyp1;
   m = (n - 1)/maxgsx + 1;
   n = n < maxgsx ? n : maxgsx;
   dim3 dimGrid(n,m);
   ns = 2*(mx + 1)*(my + 1)*sizeof(float);
   n = nblock_size*sizeof(float);
   ns = ns > n ? ns : n;
   crc = cudaGetLastError();
   gpuppgppush2l<<<dimGrid,dimBlock,ns>>>(ppart,fxy,kpic,noff,nyp,qbm,
                                          dt,ek,nx,ny,mx,my,idimp,nppmx,
                                          nxv,nypmx,mx1,mxyp1,ipbc);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppgppush2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpu2ppgppost2l(float *ppart, float *q, int *kpic,
                                int noff, float qm, int idimp,
                                int nppmx, int mx, int my, int nxv,
                                int nypmx, int mx1, int mxyp1) {
/* Deposit Interface for C */
   int n, m, ns;
   dim3 dimBlock(nblock_size);
   n = mxyp1;
   m = (n - 1)/maxgsx + 1;
   n = n < maxgsx ? n : maxgsx;
   dim3 dimGrid(n,m);
   ns = (mx + 1)*(my + 1)*sizeof(float);
   crc = cudaGetLastError();
   gpu2ppgppost2l<<<dimGrid,dimBlock,ns>>>(ppart,q,kpic,noff,qm,idimp,
                                           nppmx,mx,my,nxv,nypmx,mx1,
                                           mxyp1);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpu2ppgppost2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppcaguard2xl(float2 *qc, float2 *scs, float *q,
                                 int nyp, int nx, int nxe, int nypmx,
                                 int nxvh, int kypd) {
/* Guard Cell Interface for C */
   dim3 dimBlock(nblock_size);
   dim3 dimGrid(nyp);
   crc = cudaGetLastError();
   gpuppcaguard2xl<<<dimGrid,dimBlock>>>(qc,scs,q,nyp,nx,nxe,nypmx,nxvh,
                                         kypd);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppcaguard2xl error=%d:%s\n",crc,
             cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppcaguard2yl(float2 *qc, float2 *scr, int nx,
                                 int nxvh, int kypd) {
/* Guard Cell Interface for C */
   int nxh;
   nxh = nx/2;
   dim3 dimBlock(nblock_size);
   dim3 dimGrid((nxh - 1)/nblock_size+1);
   crc = cudaGetLastError();
   gpuppcaguard2yl<<<dimGrid,dimBlock>>>(qc,scr,nx,nxvh,kypd);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppcaguard2yl error=%d:%s\n",crc,
             cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppccguard2xl(float2 *fxyc, float2 *scs, float *fxy,
                                 int nyp, int nx, int nxe, int nypmx,
                                 int nxvh, int kypd) {
/* Guard Cell Interface for C */
   dim3 dimBlock(nblock_size);
   dim3 dimGrid(nyp);
   crc = cudaGetLastError();
   gpuppccguard2xl<<<dimGrid,dimBlock>>>(fxyc,scs,fxy,nyp,nx,nxe,nypmx,
                                         nxvh,kypd);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppccguard2xl error=%d:%s\n",crc,
             cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppccguard2yl(float *fxy, float2 *scr, int nyp,
                                 int nx, int nxe, int nxvh, int nypmx) {
/* Guard Cell Interface for C */
   int nxh;
   nxh = nx/2;
   dim3 dimBlock(nblock_size);
   dim3 dimGrid((nxh - 1)/nblock_size+1);
   crc = cudaGetLastError();
   gpuppccguard2yl<<<dimGrid,dimBlock>>>(fxy,scr,nyp,nx,nxe,nxvh,nypmx);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppccguard2yl error=%d:%s\n",crc,
             cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpupppord2la(float *ppart, float *ppbuff, float *sbufl,
                              float *sbufr, int *kpic, int *ncl,
                              int *ihole, int *ncll, int *nclr,
                              int noff, int nyp, int idimp, int nppmx,
                              int nx, int ny, int mx, int my, int mx1,
                              int myp1, int npbmx, int ntmax, int nbmax,
                              int *irc) {
/* Sort Interface for C */
   int mxyp1, n, m, ns, nbl, ierr;
   static int *g_block = NULL;
   static int lg_block = 0;
   dim3 dimBlock(nblock_size);
   mxyp1 = mx1*myp1;
   m = (mxyp1 - 1)/maxgsx + 1;
   n = mxyp1 < maxgsx ? mxyp1 : maxgsx;
   dim3 dimGrid(n,m);
/* find which particles are leaving tile */
   ns = (nblock_size+9)*sizeof(int);
   crc = cudaGetLastError();
   gpupppfnd2l<<<dimGrid,dimBlock,ns>>>(ppart,kpic,ncl,ihole,noff,nyp,
                                        idimp,nppmx,nx,ny,mx,my,mx1,
                                        myp1,ntmax,irc);
/* cudaThreadSynchronize(); */
   crc = cudaGetLastError();
   if (crc) {
      printf("gpupppfnd2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
/* buffer particles that are leaving tile and sum ncl */
   ns = 9*sizeof(int);
   crc = cudaGetLastError();
   gpupppmov2l<<<dimGrid,dimBlock,ns>>>(ppart,ppbuff,ncl,ihole,idimp,
                                        nppmx,mx1,myp1,npbmx,ntmax,irc);
/* cudaThreadSynchronize(); */
   crc = cudaGetLastError();
   if (crc) {
      printf("gpupppmov2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
/* find address offsets */
   nbl = (mx1 - 1)/nblock_size + 1;
   dim3 dimGrids(nbl);
/* allocate scratch memory needed by prefix scan */
   ierr = 0;
   if (lg_block < nbl) {
      if (lg_block > 0) 
         gpu_deallocate((void *)g_block,&ierr);
      lg_block = 0;
      gpu_iallocate(&g_block,2*nbl,&ierr);
      if (ierr != 0) {
         printf("GPU g_block allocate error!\n");
         exit(1);
      }
      lg_block = nbl;
   }
   ns = 2*nblock_size*sizeof(int);
   crc = cudaGetLastError();
   nciscan2<<<dimGrids,dimBlock,ns>>>(ncl,ncll,nclr,g_block,mx1,myp1);
/* cudaThreadSynchronize(); */
   crc = cudaGetLastError();
   if (crc) {
      printf("nciscan2 error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
/* copy particles and offsets leaving processor */
/* gpupppbuf2l and nciscan2 should use the same blocksize */
   dim3 dimGridg(mx1);
   crc = cudaGetLastError();
   gpupppbuf2l<<<dimGridg,dimBlock>>>(ppbuff,sbufl,sbufr,ncl,ncll,nclr,
                                      g_block,idimp,mx1,myp1,npbmx,
                                      nbmax,irc);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpupppbuf2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpupppord2lb(float *ppart, float *ppbuff, float *rbufl,
                              float *rbufr, int *kpic, int *ncl,
                              int *ihole, int *mcll, int *mclr,
                              int idimp, int nppmx, int mx1, int myp1,
                              int npbmx, int ntmax, int nbmax,
                              int *irc) {
/* Sort Interface for C */
   int mxyp1, n, m, ns;
   dim3 dimBlock(nblock_size);
   mxyp1 = mx1*myp1;
   m = (mxyp1 - 1)/maxgsx + 1;
   n = mxyp1 < maxgsx ? mxyp1 : maxgsx;
   dim3 dimGrid(n,m);
/* copy incoming particles from ppbuff, rbufl, and rbufr into ppart */
/* and update kpic */
   ns = (nblock_size+18)*sizeof(int);
   crc = cudaGetLastError();
   gpupppord2l<<<dimGrid,dimBlock,ns>>>(ppart,ppbuff,rbufl,rbufr,kpic,
                                        ncl,ihole,mcll,mclr,idimp,nppmx,
                                        mx1,myp1,npbmx,ntmax,nbmax,irc);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpupppord2l error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C"  void cgpuppois22t(float2 *qt, float2 *fxyt, float2 *ffct,
                              float *we, int nx, int ny, int kstrt,
                              int nyv, int kxp1, int nyhd) {
/* Poisson Solver Interface for C */
   int nxh1, ks, kxpp, ns;
   dim3 dimBlock(nblock_size);
   nxh1 = nx/2 + 1;
   ks = kstrt - 1;
   kxpp = nxh1 - kxp1*ks;
   kxpp = 0 > kxpp ? 0 : kxpp;
   kxpp = kxp1 < kxpp ? kxp1 : kxpp;
   if (kxpp <= 0)
      return;
   dim3 dimGrid(kxpp);
   ns = nblock_size*sizeof(float);
   crc = cudaGetLastError();
   gpuppois22t<<<dimGrid,dimBlock,ns>>>(qt,fxyt,ffct,we,nx,ny,kstrt,nyv,
                                        kxp1,nyhd);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppois22t error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwppfft2rcsx(float2 *f, float2 *bsm, int isign,
                                int *mixup, float2 *sct, int indx,
                                int indy, int kstrt, int nvp, int kxp1,
                                int kyp, int nxhd, int kypd, int nxhyd,
                                int nxyhd) {
/* wrapper function for parallel real to complex fft in x, */
/* without packed data                                     */
/* nxhd must be >= nx/2 + 1                                */
/* local data */
   int nxh, nxh1, ny, ks, kypp, kxyp, nsize, ns;
   static int kypi = 1, mx = 16;
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   nxh1 = nxh + 1;
   ny = 1L<<indy;
   ks = kstrt - 1;
   kypp = ny - kyp*ks;
   kypp = 0 > kypp ? 0 : kypp;
   kypp = kyp < kypp ? kyp : kypp;
   if (kypp <= 0)
      return;
   kxyp = kxp1*kyp;
   dim3 dimGridy(kypp);
   dim3 dimGrids(kypp,nvp);
   dim3 dimGridty((kyp-1)/mx+1,(kxp1-1)/mx+1,nvp);
/* inverse fourier transform */
   if (isign < 0) {
      nsize = nxh < 1024 ? nxh : 1024;
      ns = nsize*sizeof(float2);
/* perform x fft */
      if (kstrt <= ny) {
         crc = cudaGetLastError();
         gpufft2rcxs<<<dimGridy,dimBlock,ns>>>(f,isign,mixup,sct,indx,
                                               indy,kypi,kypp,nxhd,kypd,
                                               nxhyd,nxyhd,nsize);
/*       cudaThreadSynchronize(); */
         crc = cudaGetLastError();
         if (crc) {
            printf("gpufft2rcxs error=%d:%s\n",
                   crc,cudaGetErrorString(crc));
            exit(1);
         }
      }
/* extract data to send */
      crc = cudaGetLastError();
      gpuppmtposes<<<dimGrids,dimBlock>>>(f,bsm,nxh1,kxp1,kypp,kstrt,
                                          nvp,kxyp,nxhd,kypd);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposes error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
      ns = (mx+1)*mx*sizeof(float2);
/* transpose data received */
      crc = cudaGetLastError();
      gpuppmtposer<<<dimGridty,dimBlockt,ns>>>(f,bsm,nxh1,kxp1,kypp,
                                               kstrt,nvp,kxyp,nxhd,
                                               kypd);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposer error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform x fft */
      nsize = nxh < 1024 ? nxh : 1024;
      ns = nsize*sizeof(float2);
      if (kstrt <= ny) {
         gpufft2rcxs<<<dimGridy,dimBlock,ns>>>(f,isign,mixup,sct,indx,
                                               indy,kypi,kypp,nxhd,kypd,
                                               nxhyd,nxyhd,nsize);
         cudaThreadSynchronize();
         crc = cudaGetLastError();
         if (crc) {
            printf("gpufft2rcxs error=%d:%s\n",
                   crc,cudaGetErrorString(crc));
            exit(1);
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwppfft2rcsy(float2 *g, float2 *brm, int isign,
                                int *mixup, float2 *sct, int indx,
                                int indy, int kstrt, int nvp, int kxp1,
                                int kyp, int nyd, int nxhyd,
                                int nxyhd) {
/* wrapper function for parallel real to complex fft in y, */
/* without packed data                                     */
/* nyd must be >= ny                                       */
/* local data */
   int nxh, nxh1, ny, ks, kxpp, kxyp, nsize, ns;
   static int kxpi = 1, mx = 16;
/* double dtime; */
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   nxh1 = nxh + 1;
   ny = 1L<<indy;
   ks = kstrt - 1;
   kxpp = nxh1 - kxp1*ks;
   kxpp = 0 > kxpp ? 0 : kxpp;
   kxpp = kxp1 < kxpp ? kxp1 : kxpp;
   if (kxpp <= 0)
      return;
   kxyp = kxp1*kyp;
   dim3 dimGridx(kxpp);
   dim3 dimGrids(kxpp,nvp);
   dim3 dimGridtx((kxp1-1)/mx+1,(kyp-1)/mx+1,nvp);
/* inverse fourier transform */
   if (isign < 0) {
      ns = (mx+1)*mx*sizeof(float2);
/* transpose data received */
      crc = cudaGetLastError();
      gpuppmtposer<<<dimGridtx,dimBlockt,ns>>>(g,brm,ny,kyp,kxpp,kstrt,
                                               nvp,kxyp,nyd,kxp1);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposer error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform y fft */
      nsize = ny < 1024 ? ny : 1024;
      ns = nsize*sizeof(float2);
      if (kstrt <= nxh1) {
         crc = cudaGetLastError();
         gpufft2rcys<<<dimGridx,dimBlock,ns>>>(g,isign,mixup,sct,indx,
                                               indy,kxpi,kxpp,kxp1,nyd,
                                               nxhyd,nxyhd,nsize);
         cudaThreadSynchronize();
         crc = cudaGetLastError();
         if (crc) {
            printf("gpufft2rcys error=%d:%s\n",
                   crc,cudaGetErrorString(crc));
            exit(1);
         }
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
      nsize = ny < 1024 ? ny : 1024;
      ns = nsize*sizeof(float2);
/* perform y fft */
      if (kstrt <= nxh1) {
         crc = cudaGetLastError();
         gpufft2rcys<<<dimGridx,dimBlock,ns>>>(g,isign,mixup,sct,indx,
                                               indy,kxpi,kxpp,kxp1,nyd,
                                               nxhyd,nxyhd,nsize);
/*       cudaThreadSynchronize(); */
         crc = cudaGetLastError();
         if (crc) {
            printf("gpufft2rcys error=%d:%s\n",
                   crc,cudaGetErrorString(crc));
            exit(1);
         }
      }
/* extract data to send */
      crc = cudaGetLastError();
      gpuppmtposes<<<dimGrids,dimBlock>>>(g,brm,ny,kyp,kxpp,kstrt,
                                          nvp,kxyp,nyd,kxp1);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposes error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }


   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwppfft2rcsxn(float2 *fn, float2 *bsm, int isign,
                                 int *mixup, float2 *sct, int indx,
                                 int indy, int ndim, int kstrt, int nvp,
                                 int kxp1, int kyp, int nxhd, int kypd,
                                 int nxhyd, int nxyhd) {
/* wrapper function for multiple parallel real to complex ffts in x, */
/* without packed data                                               */
/* ndim = vector dimension                                           */
/* nxhd must be >= nx/2 + 1                                          */
/* local data */
   int nxh, nxh1, ny, ks, kypp, kxyp, nkypd, nkypp, nsize, ns;
   static int kypi = 1, mx = 16;
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   nxh1 = nxh + 1;
   ny = 1L<<indy;
   ks = kstrt - 1;
   kypp = ny - kyp*ks;
   kypp = 0 > kypp ? 0 : kypp;
   kypp = kyp < kypp ? kyp : kypp;
   if (kypp <= 0)
      return;
   kxyp = kxp1*kyp;
   nkypd = ndim*kypd;
   nkypp = ndim*kypp;
   dim3 dimGridy(nkypp);
   dim3 dimGrids(kypp,nvp);
   dim3 dimGridty((kyp-1)/mx+1,(kxp1-1)/mx+1,nvp);
/* inverse fourier transform */
   if (isign < 0) {
      nsize = nxh < 1024 ? nxh : 1024;
      ns = nsize*sizeof(float2);
/* perform x fft */
      if (kstrt <= ny) {
         crc = cudaGetLastError();
         gpufft2rcxs<<<dimGridy,dimBlock,ns>>>(fn,isign,mixup,sct,indx,
                                               indy,kypi,nkypp,nxhd,
                                               nkypd,nxhyd,nxyhd,nsize);
/*       cudaThreadSynchronize(); */
         crc = cudaGetLastError();
         if (crc) {
            printf("gpufft2rcxs error=%d:%s\n",
                   crc,cudaGetErrorString(crc));
            exit(1);
         }
      }
/* extract data to send */
      crc = cudaGetLastError();
      gpuppmtposesn<<<dimGrids,dimBlock>>>(fn,bsm,nxh1,kxp1,kypp,kstrt,
                                           nvp,ndim,kxyp,nxhd,kypd);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposesn error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
      ns = ndim*(mx+1)*mx*sizeof(float2);
/* transpose data received */
      crc = cudaGetLastError();
      gpuppmtposern<<<dimGridty,dimBlockt,ns>>>(fn,bsm,nxh1,kxp1,kypp,
                                                kstrt,nvp,ndim,kxyp,
                                                nxhd,kypd);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposern error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform x fft */
      nsize = nxh < 1024 ? nxh : 1024;
      ns = nsize*sizeof(float2);
      if (kstrt <= ny) {
         crc = cudaGetLastError();
         gpufft2rcxs<<<dimGridy,dimBlock,ns>>>(fn,isign,mixup,sct,indx,
                                               indy,kypi,nkypp,nxhd,
                                               nkypd,nxhyd,nxyhd,nsize);
         cudaThreadSynchronize();
         crc = cudaGetLastError();
         if (crc) {
            printf("gpufft2rcxs error=%d:%s\n",
                   crc,cudaGetErrorString(crc));
            exit(1);
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwppfft2rcsyn(float2 *gn, float2 *brm, int isign,
                                 int *mixup, float2 *sct, int indx,
                                 int indy, int ndim, int kstrt, int nvp,
                                 int kxp1, int kyp, int nyd, int nxhyd,
                                 int nxyhd) {
/* wrapper function for multiple parallel real to complex ffts in y, */
/* without packed data                                               */
/* ndim = vector dimension                                           */
/* nyd must be >= ny                                                 */
/* local data */
   int nxh, nxh1, ny, ks, kxpp, kxyp, nkxp1, nkxpp, nsize, ns;
   static int kxpi = 1, mx = 16;
/* double dtime; */
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nxh = 1L<<(indx - 1);
   nxh1 = nxh + 1;
   ny = 1L<<indy;
   ks = kstrt - 1;
   kxpp = nxh1 - kxp1*ks;
   kxpp = 0 > kxpp ? 0 : kxpp;
   kxpp = kxp1 < kxpp ? kxp1 : kxpp;
   if (kxpp <= 0)
      return;
   kxyp = kxp1*kyp;
   nkxp1 = ndim*kxp1;
   nkxpp = ndim*kxpp;
   dim3 dimGridx(nkxpp);
   dim3 dimGrids(kxpp,nvp);
   dim3 dimGridtx((kxp1-1)/mx+1,(kyp-1)/mx+1,nvp);
/* inverse fourier transform */
   if (isign < 0) {
      ns = ndim*(mx+1)*mx*sizeof(float2);
/* transpose data received */
      crc = cudaGetLastError();
      gpuppmtposern<<<dimGridtx,dimBlockt,ns>>>(gn,brm,ny,kyp,kxpp,
                                                kstrt,nvp,ndim,kxyp,nyd,
                                                kxp1);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposern error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform y fft */
      nsize = ny < 1024 ? ny : 1024;
      ns = nsize*sizeof(float2);
      if (kstrt <= nxh1) {
         crc = cudaGetLastError();
         gpufft2rcys<<<dimGridx,dimBlock,ns>>>(gn,isign,mixup,sct,indx,
                                               indy,kxpi,nkxpp,nkxp1,
                                               nyd,nxhyd,nxyhd,nsize);
         cudaThreadSynchronize();
         crc = cudaGetLastError();
         if (crc) {
            printf("gpufft2rcys error=%d:%s\n",
                   crc,cudaGetErrorString(crc));
            exit(1);
         }
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      nsize = ny < 1024 ? ny : 1024;
      ns = nsize*sizeof(float2);
      if (kstrt <= nxh1) {
         crc = cudaGetLastError();
         gpufft2rcys<<<dimGridx,dimBlock,ns>>>(gn,isign,mixup,sct,indx,
                                               indy,kxpi,nkxpp,nkxp1,
                                               nyd,nxhyd,nxyhd,nsize);
/*       cudaThreadSynchronize(); */
         crc = cudaGetLastError();
         if (crc) {
            printf("gpufft2rcys error=%d:%s\n",
                   crc,cudaGetErrorString(crc));
            exit(1);
         }
      }
/* extract data to send */
      crc = cudaGetLastError();
      gpuppmtposesn<<<dimGrids,dimBlock>>>(gn,brm,ny,kyp,kxpp,kstrt,
                                           nvp,ndim,kxyp,nyd,kxp1);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposesn error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppltpose(float2 *f, float2 *g, int nx, int ny,
                             int kxp, int kyp, int kstrt, int nxv,
                             int nyv) {
/* local complex transpose using blocking algorithm with gaps */
/* input = f, output = g                                      */
/* local data */
   int ns;
   static int mx = 16;
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   dim3 dimGridtx((kxp-1)/mx+1,(kyp-1)/mx+1);
   ns = (mx+1)*mx*sizeof(float2);
/* local transpose f to g */
   crc = cudaGetLastError();
   gpuppltpose<<<dimGridtx,dimBlockt,ns>>>(f,g,nx,ny,kxp,kyp,kstrt,nxv,
                                           nyv);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppltpose error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppltposen(float2 *fn, float2 *gn, int nx, int ny,
                              int kxp, int kyp, int kstrt, int ndim,
                              int nxv, int nyv) {
/* local complex vector transpose */
/* input = fn, output = gn        */
/* local data */
   int ns;
   static int mx = 16;
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   dim3 dimGridtx((kxp-1)/mx+1,(kyp-1)/mx+1);
   ns = ndim*(mx+1)*mx*sizeof(float2);
/* local transpose f to g */
   crc = cudaGetLastError();
   gpuppltposen<<<dimGridtx,dimBlockt,ns>>>(fn,gn,nx,ny,kxp,kyp,kstrt,
                                            ndim,nxv,nyv);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppltposen error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
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
extern "C" void cgpuppgppush2l_(unsigned long *gp_ppart,
                                unsigned long *gp_fxy,
                                unsigned long *gp_kpic,
                                int *noff, int *nyp, float *qbm,
                                float *dt, unsigned long *gp_ek,
                                int *nx, int *ny, int *mx, int *my,
                                int *idimp, int *nppmx, int *nxv,
                                int *nypmx, int *mx1, int *mxyp1,
                                int *ipbc) {
   float *ppart, *fxy, *ek;
   int *kpic;
   ppart = (float *)*gp_ppart;
   fxy = (float *)*gp_fxy;
   kpic = (int *)*gp_kpic;
   ek = (float *)*gp_ek;
   cgpuppgppush2l(ppart,fxy,kpic,*noff,*nyp,*qbm,*dt,ek,*nx,*ny,*mx,*my,
                  *idimp,*nppmx,*nxv,*nypmx,*mx1,*mxyp1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpu2ppgppost2l_(unsigned long *gp_ppart,
                                 unsigned long *gp_q,
                                 unsigned long *gp_kpic,
                                 int *noff, float *qm, int *idimp,
                                 int *nppmx, int *mx, int *my, int *nxv,
                                 int *nypmx, int *mx1, int *mxyp1) {
   float *ppart, *q;
   int *kpic;
   ppart = (float *)*gp_ppart;
   q = (float *)*gp_q;
   kpic = (int *)*gp_kpic;
   cgpu2ppgppost2l(ppart,q,kpic,*noff,*qm,*idimp,*nppmx,*mx,*my,*nxv,
                   *nypmx,*mx1,*mxyp1);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppcaguard2xl_(unsigned long *gp_qc,
                                  unsigned long *gp_scs,
                                  unsigned long *gp_q, int *nyp,
                                  int *nx, int *nxe, int *nypmx,
                                  int *nxvh, int *kypd) {
   float2 *qc, *scs;
   float *q;
   qc = (float2 *)*gp_qc;
   scs = (float2 *)*gp_scs;
   q = (float *)*gp_q;
   cgpuppcaguard2xl(qc,scs,q,*nyp,*nx,*nxe,*nypmx,*nxvh,*kypd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppcaguard2yl_(unsigned long *gp_fc,
                                 unsigned long *gp_scr, int *nx,
                                 int *nxvh, int *kypd) {
   float2 *fc, *scr;
   fc = (float2 *)*gp_fc;
   scr = (float2 *)*gp_scr;
   cgpuppcaguard2yl(fc,scr,*nx,*nxvh,*kypd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppccguard2xl_(unsigned long *gp_fxyc,
                                  unsigned long *gp_scs,
                                  unsigned long *gp_fxy, int *nyp,
                                  int *nx, int *nxe, int *nypmx,
                                  int *nxvh, int *kypd) {
   float2 *fxyc, *scs;
   float *fxy;
   fxyc = (float2 *)*gp_fxyc;
   scs = (float2 *)*gp_scs;
   fxy = (float *)*gp_fxy;
   cgpuppccguard2xl(fxyc,scs,fxy,*nyp,*nx,*nxe,*nypmx,*nxvh,*kypd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppccguard2yl_(unsigned long *gp_fxy,
                                  unsigned long *gp_scr, int *nyp,
                                  int *nx, int *nxe, int *nxvh,
                                  int *nypmx) {
   float *fxy;
   float2 *scr;
   fxy = (float *)*gp_fxy;
   scr = (float2 *)*gp_scr;
   cgpuppccguard2yl(fxy,scr,*nyp,*nx,*nxe,*nxvh,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpupppord2la_(unsigned long *gp_ppart,
                               unsigned long *gp_ppbuff,
                               unsigned long *gp_sbufl,
                               unsigned long *gp_sbufr,
                               unsigned long *gp_kpic,
                               unsigned long *gp_ncl,
                               unsigned long *gp_ihole,
                               unsigned long *gp_ncll,
                               unsigned long *gp_nclr,
                               int *noff, int *nyp, int *idimp,
                               int *nppmx, int *nx, int *ny, int *mx,
                               int *my, int *mx1, int *myp1, int *npbmx,
                               int *ntmax, int *nbmax,
                               unsigned long *gp_irc) {
   float *ppart, *ppbuff, *sbufl, *sbufr;
   int *kpic, *ncl, *ihole, *ncll, *nclr, *irc;
   ppart = (float *)*gp_ppart;
   ppbuff = (float *)*gp_ppbuff;
   sbufl = (float *)*gp_sbufl;
   sbufr = (float *)*gp_sbufr;
   kpic = (int *)*gp_kpic;
   ncl = (int *)*gp_ncl;
   ihole = (int *)*gp_ihole;
   ncll = (int *)*gp_ncll;
   nclr = (int *)*gp_nclr;
   irc = (int *)*gp_irc;
   cgpupppord2la(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,nclr,
                 *noff,*nyp,*idimp,*nppmx,*nx,*ny,*mx,*my,*mx1,*myp1,
                 *npbmx,*ntmax,*nbmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpupppord2lb_(unsigned long *gp_ppart,
                               unsigned long *gp_ppbuff,
                               unsigned long *gp_rbufl,
                               unsigned long *gp_rbufr,
                               unsigned long *gp_kpic,
                               unsigned long *gp_ncl,
                               unsigned long *gp_ihole,
                               unsigned long *gp_mcll,
                               unsigned long *gp_mclr,
                               int *idimp, int *nppmx, int *mx1,
                               int *myp1, int *npbmx, int *ntmax,
                               int *nbmax, unsigned long *gp_irc) {
   float *ppart, *ppbuff, *rbufl, *rbufr;
   int *kpic, *ncl, *ihole, *mcll, *mclr, *irc;
   ppart = (float *)*gp_ppart;
   ppbuff = (float *)*gp_ppbuff;
   rbufl = (float *)*gp_rbufl;
   rbufr = (float *)*gp_rbufr;
   kpic = (int *)*gp_kpic;
   ncl = (int *)*gp_ncl;
   ihole = (int *)*gp_ihole;
   mcll = (int *)*gp_mcll;
   mclr = (int *)*gp_mclr;
   irc = (int *)*gp_irc;
   cgpupppord2lb(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,mclr,
                 *idimp,*nppmx,*mx1,*myp1,*npbmx,*ntmax,*nbmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
extern "C"  void cgpuppois22t_(unsigned long *gp_qt, 
                               unsigned long *gp_fxyt,
                               unsigned long *gp_ffct,
                               unsigned long *gp_we, int *nx, int *ny,
                               int *kstrt, int *nyv, int *kxp1,
                               int *nyhd) {
   float2 *qt, *fxyt, *ffct;
   float *we;
   qt = (float2 *)*gp_qt;
   fxyt = (float2 *)*gp_fxyt;
   ffct = (float2 *)*gp_ffct;
   we = (float *)*gp_we;
   cgpuppois22t(qt,fxyt,ffct,we,*nx,*ny,*kstrt,*nyv,*kxp1,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwppfft2rcsx_(unsigned long *gp_f,
                                 unsigned long *gp_bsm, int *isign,
                                 unsigned long *gp_mixup,
                                 unsigned long *gp_sct, int *indx,
                                 int *indy, int *kstrt, int *nvp,
                                 int *kxp1, int *kyp, int *nxhd,
                                 int *kypd, int *nxhyd, int *nxyhd) {
   float2 *f, *bsm, *sct;
   int *mixup;
   f = (float2 *)*gp_f;
   bsm = (float2 *)*gp_bsm;
   mixup = (int *)*gp_mixup;
   sct = (float2 *)*gp_sct;
   cgpuwppfft2rcsx(f,bsm,*isign,mixup,sct,*indx,*indy,*kstrt,*nvp,*kxp1,
                   *kyp,*nxhd,*kypd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwppfft2rcsy_(unsigned long *gp_g,
                                 unsigned long *gp_brm, int *isign,
                                 unsigned long *gp_mixup,
                                 unsigned long *gp_sct, int *indx,
                                 int *indy, int *kstrt, int *nvp,
                                 int *kxp1, int *kyp, int *nyd,
                                 int *nxhyd, int *nxyhd) {
   float2 *g, *brm, *sct;
   int *mixup;
   g = (float2 *)*gp_g;
   brm = (float2 *)*gp_brm;
   mixup = (int *)*gp_mixup;
   sct = (float2 *)*gp_sct;
   cgpuwppfft2rcsy(g,brm,*isign,mixup,sct,*indx,*indy,*kstrt,*nvp,*kxp1,
                   *kyp,*nyd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwppfft2rcsxn_(unsigned long *gp_fn,
                                  unsigned long *gp_bsm, int *isign,
                                  unsigned long *gp_mixup,
                                  unsigned long *gp_sct, int *indx,
                                  int *indy, int *ndim, int *kstrt,
                                  int *nvp, int *kxp1, int *kyp,
                                  int *nxhd, int *kypd, int *nxhyd,
                                  int *nxyhd) {
   float2 *fn, *bsm, *sct;
   int *mixup;
   fn = (float2 *)*gp_fn;
   bsm = (float2 *)*gp_bsm;
   mixup = (int *)*gp_mixup;
   sct = (float2 *)*gp_sct;
   cgpuwppfft2rcsxn(fn,bsm,*isign,mixup,sct,*indx,*indy,*ndim,*kstrt,
                    *nvp,*kxp1,*kyp,*nxhd,*kypd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuwppfft2rcsyn_(unsigned long *gp_gn,
                                  unsigned long *gp_brm, int *isign,
                                  unsigned long *gp_mixup,
                                  unsigned long *gp_sct, int *indx,
                                  int *indy, int *ndim, int *kstrt,
                                  int *nvp, int *kxp1, int *kyp,
                                  int *nyd, int *nxhyd, int *nxyhd) {
   float2 *gn, *brm, *sct;
   int *mixup;
   gn = (float2 *)*gp_gn;
   brm = (float2 *)*gp_brm;
   mixup = (int *)*gp_mixup;
   sct = (float2 *)*gp_sct;
   cgpuwppfft2rcsyn(gn,brm,*isign,mixup,sct,*indx,*indy,*ndim,*kstrt,
                    *nvp,*kxp1,*kyp,*nyd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppltpose_(unsigned long *gp_f, unsigned long *gp_g,
                              int *nx, int *ny, int *kxp, int *kyp,
                              int *kstrt, int *nxv, int *nyv) {
   float2 *f, *g;
   f = (float2 *)*gp_f;
   g = (float2 *)*gp_g;
   cgpuppltpose(f,g,*nx,*ny,*kxp,*kyp,*kstrt,*nxv,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void cgpuppltposen_(unsigned long *gp_fn,
                               unsigned long *gp_gn, int *nx, int *ny,
                               int *kxp, int *kyp, int *kstrt, int *ndim,
                               int *nxv, int *nyv) {
   float2 *fn, *gn;
   fn = (float2 *)*gp_fn;
   gn = (float2 *)*gp_gn;
   cgpuppltposen(fn,gn,*nx,*ny,*kxp,*kyp,*kstrt,*ndim,*nxv,*nyv);
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
