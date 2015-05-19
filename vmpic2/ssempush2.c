/* SSE2 C Library for Skeleton 2D Electrostatic OpenMP/Vector PIC Code */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <xmmintrin.h>
#include "ssempush2.h"

void cfft2r2x(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nyi, int nyp,
              int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2r2y(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd);

/*--------------------------------------------------------------------*/
void csse2gppush2lt(float ppart[], float fxy[], int kpic[], float qbm,
                    float dt, float *ek, int idimp, int nppmx, int nx,
                    int ny, int mx, int my, int nxv, int nyv, int mx1,
                    int mxy1, int ipbc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions.
   OpenMP/vector version using guard cells
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
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   requires SSE2, ppart needs to be 16 byte aligned
   nppmx needs to be a multiple of 4
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, nps;
   int i, j, k, nn, mm, mxv;
   float qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
   double sum1, sum2;
   __m128i v_noff, v_moff, v_mxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qtm, v_dt, v_one;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_vx, v_vy;
   __m128 v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 a, b, c, d;
   __m128d v_sum1, v_d;
   __attribute__((aligned(16))) unsigned int ll[4];
   __attribute__((aligned(16))) double dd[2];
   __attribute__((aligned(16))) float sfxy[2*MXV*MYV];
/* __attribute__((aligned(16))) float sfxy[2*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   qtm = qbm*dt;
   sum2 = 0.0;
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
   v_mxv = _mm_set1_epi32(mxv);
   v_qtm = _mm_set1_ps(qtm);
   v_one = _mm_set1_ps(1.0f);
   v_dt = _mm_set1_ps(dt);
   v_edgelx = _mm_set1_ps(edgelx);
   v_edgely = _mm_set1_ps(edgely);
   v_edgerx = _mm_set1_ps(edgerx);
   v_edgery = _mm_set1_ps(edgery);
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,nps,npoff,nn,mm,x,y,dxp,dyp,amx,amy,dx,dy, \
vx,vy,sum1,v_noff,v_moff,v_nn,v_mm,v_it,v_x,v_y,v_dxp,v_dyp,v_amx,v_amy, \
v_dx,v_dy,v_vx,v_vy,v_at,v_d,v_sum1,a,b,c,d,ll,dd,sfxy)  \
reduction(+:sum2)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm_set1_epi32(noff);
      v_moff = _mm_set1_epi32(moff);
      npp = kpic[k];
      npoff = idimp*nppmx*k;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      nps = 4*((2*nn)/4);
      for (j = 0; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*       for (i = 0; i < nn; i++) {                               */
/*          sfxy[2*(i+mxv*j)] = fxy[2*(i+noff+nxv*(j+moff))];     */
/*          sfxy[1+2*(i+mxv*j)] = fxy[1+2*(i+noff+nxv*(j+moff))]; */
/*       }                                                        */
         for (i = 0; i < nps; i+=4) {
            v_at = _mm_loadu_ps(&fxy[i+2*(noff+nxv*(j+moff))]);
           _mm_storeu_ps(&sfxy[i+2*(mxv*j)],v_at);
         }
/* loop over remaining elements */
         for (i = nps; i < 2*nn; i++) {
            sfxy[i+2*(mxv*j)] = fxy[i+2*(noff+nxv*(j+moff))];
         }
      }
      nps = 4*(npp/4);
      sum1 = 0.0;
      v_sum1 = _mm_set1_pd(0.0);
/* loop over particles in tile in groups of 4 */
      for (j = 0; j < nps; j+=4) {
/* find interpolation weights */
/*       x = ppart[j+npoff];       */
/*       y = ppart[j+nppmx+npoff]; */
         v_x = _mm_load_ps(&ppart[j+npoff]);
         v_y = _mm_load_ps(&ppart[j+nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
         v_nn = _mm_cvttps_epi32(v_x);
         v_mm = _mm_cvttps_epi32(v_y);
/*       dxp = x - (float) nn; */
         v_dxp = _mm_sub_ps(v_x,_mm_cvtepi32_ps(v_nn));
/*       dyp = y - (float) mm; */
         v_dyp = _mm_sub_ps(v_y,_mm_cvtepi32_ps(v_mm));
/*       nn = 2*(nn - noff + mxv*(mm - moff)); */
         v_nn = _mm_sub_epi32(v_nn,v_noff);
         v_mm = _mm_sub_epi32(v_mm,v_moff);
         v_it = _mm_mul_epu32(v_mxv,_mm_srli_si128(v_mm,4));
         v_mm = _mm_mul_epu32(v_mm,v_mxv);
         v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
         v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),1);
/*       amx = 1.0f - dxp; */
/*       amy = 1.0f - dyp; */
         v_amx = _mm_sub_ps(v_one,v_dxp);
         v_amy = _mm_sub_ps(v_one,v_dyp);
/* find acceleration */
/* load fields, for lower left/right */
         _mm_store_si128((__m128i *)ll,v_nn);
         a = _mm_loadu_ps(&sfxy[ll[0]]);
         b = _mm_loadu_ps(&sfxy[ll[1]]);
         c = _mm_loadu_ps(&sfxy[ll[2]]);
         d = _mm_loadu_ps(&sfxy[ll[3]]);
/* transpose so a,b,c,d contain first 4 fields for each of 4 particles */
         _MM_TRANSPOSE4_PS(a,b,c,d);
/*       dx = amx*sfxy[nn];   */
/*       dy = amx*sfxy[nn+1]; */
         v_dx = _mm_mul_ps(v_amx,a);
         v_dy = _mm_mul_ps(v_amx,b);
/*       dx = amy*(dxp*sfxy[nn+2] + dx); */
/*       dy = amy*(dxp*sfxy[nn+3] + dy); */
         v_dx = _mm_mul_ps(v_amy,_mm_add_ps(_mm_mul_ps(v_dxp,c),v_dx));
         v_dy = _mm_mul_ps(v_amy,_mm_add_ps(_mm_mul_ps(v_dxp,d),v_dy));
/*       nn += 2*mxv; */
/* load fields, for upper left/right */
         a = _mm_loadu_ps(&sfxy[ll[0]+2*mxv]);
         b = _mm_loadu_ps(&sfxy[ll[1]+2*mxv]);
         c = _mm_loadu_ps(&sfxy[ll[2]+2*mxv]);
         d = _mm_loadu_ps(&sfxy[ll[3]+2*mxv]);
/* transpose so a,b,c,d contain next 4 fields for each of 4 particles */
         _MM_TRANSPOSE4_PS(a,b,c,d);
/*       vx = amx*sfxy[nn];   */
/*       vy = amx*sfxy[nn+1]; */
         a = _mm_mul_ps(v_amx,a);
         b = _mm_mul_ps(v_amx,b);
/*       dx += dyp*(dxp*sfxy[nn+2] + vx); */
/*       dy += dyp*(dxp*sfxy[nn+3] + vy); */
         a = _mm_mul_ps(v_dyp,_mm_add_ps(_mm_mul_ps(v_dxp,c),a));
         b = _mm_mul_ps(v_dyp,_mm_add_ps(_mm_mul_ps(v_dxp,d),b));
         v_dx = _mm_add_ps(v_dx,a);
         v_dy = _mm_add_ps(v_dy,b);
/* new velocity */
/*       dxp = ppart[j+2*nppmx+npoff]; */
/*       dyp = ppart[j+3*nppmx+npoff]; */
         v_dxp = _mm_load_ps(&ppart[j+2*nppmx+npoff]);
         v_dyp = _mm_load_ps(&ppart[j+3*nppmx+npoff]);
/*       vx = dxp + qtm*dx; */
/*       vy = dyp + qtm*dy; */
         v_vx = _mm_add_ps(v_dxp,_mm_mul_ps(v_qtm,v_dx));
         v_vy = _mm_add_ps(v_dyp,_mm_mul_ps(v_qtm,v_dy));
/* average kinetic energy */
/*       dxp += vx; */
/*       dyp += vy; */
         v_dxp = _mm_add_ps(v_dxp,v_vx);
         v_dyp = _mm_add_ps(v_dyp,v_vy);
/*       sum1 += dxp*dxp + dyp*dyp; */
         v_at = _mm_mul_ps(v_dxp,v_dxp);
         v_at = _mm_add_ps(v_at,_mm_mul_ps(v_dyp,v_dyp));
/* convert to double precision before accumulating */
         v_d = _mm_cvtps_pd(v_at);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
         v_it = _mm_srli_si128((__m128i)v_at,8);
         v_d = _mm_cvtps_pd((__m128)v_it);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
/* new position */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
         v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_dt));
         v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_dt));
/* reflecting boundary conditions */
         if (ipbc==2) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             vx = -vx;                           */
/*          }                                      */
            v_at = _mm_cmplt_ps(v_dx,v_edgelx);
            v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dx,v_edgerx));
            v_x = _mm_and_ps(v_at,v_x);
            v_dx = _mm_add_ps(_mm_andnot_ps(v_at,v_dx),v_x);
            v_dxp = _mm_and_ps(v_at,v_vx);
            v_vx = _mm_sub_ps(_mm_andnot_ps(v_at,v_vx),v_dxp);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             vy = -vy;                           */
/*          }                                      */
            v_at = _mm_cmplt_ps(v_dy,v_edgely);
            v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dy,v_edgery));
            v_y = _mm_and_ps(v_at,v_y);
            v_dy = _mm_add_ps(_mm_andnot_ps(v_at,v_dy),v_y);
            v_dyp = _mm_and_ps(v_at,v_vy);
            v_vy = _mm_sub_ps(_mm_andnot_ps(v_at,v_vy),v_dyp);
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             vx = -vx;                           */
/*          }                                      */
            v_at = _mm_cmplt_ps(v_dx,v_edgelx);
            v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dx,v_edgerx));
            v_x = _mm_and_ps(v_at,v_x);
            v_dx = _mm_add_ps(_mm_andnot_ps(v_at,v_dx),v_x);
            v_dxp = _mm_and_ps(v_at,v_vx);
            v_vx = _mm_sub_ps(_mm_andnot_ps(v_at,v_vx),v_dxp);
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy; */
         _mm_store_ps(&ppart[j+npoff],v_dx);
         _mm_store_ps(&ppart[j+nppmx+npoff],v_dy);
/* set new velocity */
/*       ppart[j+2*nppmx+npoff] = vx; */
/*       ppart[j+3*nppmx+npoff] = vy; */
         _mm_store_ps(&ppart[j+2*nppmx+npoff],v_vx);
         _mm_store_ps(&ppart[j+3*nppmx+npoff],v_vy);
      }
/* loop over remaining particles in tile */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nn = 2*(nn - noff + mxv*(mm - moff));
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
         dxp = ppart[j+2*nppmx+npoff];
         dyp = ppart[j+3*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         sum1 += dxp*dxp + dyp*dyp;
/* new position */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
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
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
/* set new velocity */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
      }
/*    sum2 += sum1; */
      _mm_store_pd(&dd[0],v_sum1);
      for (j = 1; j < 2; j++) {
         dd[0] += dd[j];
      }
      sum2 += (sum1 + dd[0]);
   }
/* normalize kinetic energy */
   *ek += 0.125f*sum2;
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void csse2gppushf2lt(float ppart[], float fxy[], int kpic[], int ncl[],
                     int ihole[], float qbm, float dt, float *ek, 
                     int idimp, int nppmx, int nx, int ny, int mx,
                     int my, int nxv, int nyv, int mx1, int mxy1,
                     int ntmax, int *irc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   OpenMP/vector version using guard cells
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
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   requires SSE2, ppart needs to be 16 byte aligned
   nppmx needs to be a multiple of 4
   optimized version
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, nps;
   int i, j, k, ih, nh, nn, mm, kk, mxv;
   float qtm, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
   float anx, any, edgelx, edgely, edgerx, edgery;
   double sum1, sum2;
   __m128i v_noff, v_moff, v_mxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qtm, v_dt, v_one;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_st, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_vx, v_vy;
   __m128 v_anx, v_any, v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 v_zero, v_two, v_three, v_six;
   __m128 a, b, c, d;
   __m128d v_sum1, v_d;
   __attribute__((aligned(16))) unsigned int ll[4], lm[8];
   __attribute__((aligned(16))) unsigned long jj[1];
   __attribute__((aligned(16))) double dd[2];
   __attribute__((aligned(16))) float sfxy[2*MXV*MYV];
/* __attribute__((aligned(16))) float sfxy[2*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   qtm = qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   sum2 = 0.0;
   v_mxv = _mm_set1_epi32(mxv);
   v_qtm = _mm_set1_ps(qtm);
   v_one = _mm_set1_ps(1.0f);
   v_dt = _mm_set1_ps(dt);
   v_anx = _mm_set1_ps(anx);
   v_any = _mm_set1_ps(any);
   v_zero = _mm_setzero_ps();
   v_two = _mm_set1_ps(2.0f);
   v_three = _mm_set1_ps(3.0f);
   v_six = _mm_set1_ps(6.0f);
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,nps,npoff,nn,mm,kk,ih,nh,x,y,dxp,dyp,amx, \
amy,dx,dy,vx,vy,edgelx,edgely,edgerx,edgery,sum1,v_noff,v_moff,v_nn, \
v_mm,v_it,v_x,v_y,v_dxp,v_dyp,v_amx,v_amy,v_dx,v_dy,v_vx,v_vy,v_st, \
v_at,v_edgelx,v_edgely,v_edgerx,v_edgery,v_d,v_sum1,a,b,c,d,jj,ll,lm, \
dd,sfxy)  \
reduction(+:sum2)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm_set1_epi32(noff);
      v_moff = _mm_set1_epi32(moff);
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
      v_edgelx = _mm_set1_ps(edgelx);
      v_edgely = _mm_set1_ps(edgely);
      v_edgerx = _mm_set1_ps(edgerx);
      v_edgery = _mm_set1_ps(edgery);
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      nps = 4*((2*nn)/4);
      for (j = 0; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*       for (i = 0; i < nn; i++) {                               */
/*          sfxy[2*(i+mxv*j)] = fxy[2*(i+noff+nxv*(j+moff))];     */
/*          sfxy[1+2*(i+mxv*j)] = fxy[1+2*(i+noff+nxv*(j+moff))]; */
/*       }                                                        */
         for (i = 0; i < nps; i+=4) {
            v_at = _mm_loadu_ps(&fxy[i+2*(noff+nxv*(j+moff))]);
           _mm_storeu_ps(&sfxy[i+2*(mxv*j)],v_at);
         }
/* loop over remaining elements */
         for (i = nps; i < 2*nn; i++) {
            sfxy[i+2*(mxv*j)] = fxy[i+2*(noff+nxv*(j+moff))];
         }
      }
/* clear counters */
/*    for (j = 0; j < 8; j++) { */
/*       ncl[j+8*k] = 0;        */
/*    }                         */
      memset((void*)&ncl[8*k],0,8*sizeof(int));
      nps = 4*(npp/4);
      sum1 = 0.0;
      v_sum1 = _mm_set1_pd(0.0);
/* loop over particles in tile in groups of 4 */
      for (j = 0; j < nps; j+=4) {
/* find interpolation weights */
/*       x = ppart[j+npoff];       */
/*       y = ppart[j+nppmx+npoff]; */
         v_x = _mm_load_ps(&ppart[j+npoff]);
         v_y = _mm_load_ps(&ppart[j+nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
         v_nn = _mm_cvttps_epi32(v_x);
         v_mm = _mm_cvttps_epi32(v_y);
/*       dxp = x - (float) nn; */
         v_dxp = _mm_sub_ps(v_x,_mm_cvtepi32_ps(v_nn));
/*       dyp = y - (float) mm; */
         v_dyp = _mm_sub_ps(v_y,_mm_cvtepi32_ps(v_mm));
/*       nn = 2*(nn - noff + mxv*(mm - moff)); */
         v_nn = _mm_sub_epi32(v_nn,v_noff);
         v_mm = _mm_sub_epi32(v_mm,v_moff);
         v_it = _mm_mul_epu32(v_mxv,_mm_srli_si128(v_mm,4));
         v_mm = _mm_mul_epu32(v_mm,v_mxv);
         v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
         v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),1);
/*       amx = 1.0f - dxp; */
/*       amy = 1.0f - dyp; */
         v_amx = _mm_sub_ps(v_one,v_dxp);
         v_amy = _mm_sub_ps(v_one,v_dyp);
/* find acceleration */
/* load fields, for lower left/right */
         _mm_store_si128((__m128i *)ll,v_nn);
         a = _mm_loadu_ps(&sfxy[ll[0]]);
         b = _mm_loadu_ps(&sfxy[ll[1]]);
         c = _mm_loadu_ps(&sfxy[ll[2]]);
         d = _mm_loadu_ps(&sfxy[ll[3]]);
/* transpose so a,b,c,d contain first 4 fields for each of 4 particles */
         _MM_TRANSPOSE4_PS(a,b,c,d);
/*       dx = amx*sfxy[nn];   */
/*       dy = amx*sfxy[nn+1]; */
         v_dx = _mm_mul_ps(v_amx,a);
         v_dy = _mm_mul_ps(v_amx,b);
/*       dx = amy*(dxp*sfxy[nn+2] + dx); */
/*       dy = amy*(dxp*sfxy[nn+3] + dy); */
         v_dx = _mm_mul_ps(v_amy,_mm_add_ps(_mm_mul_ps(v_dxp,c),v_dx));
         v_dy = _mm_mul_ps(v_amy,_mm_add_ps(_mm_mul_ps(v_dxp,d),v_dy));
/*       nn += 2*mxv; */
/* load fields, for upper left/right */
         a = _mm_loadu_ps(&sfxy[ll[0]+2*mxv]);
         b = _mm_loadu_ps(&sfxy[ll[1]+2*mxv]);
         c = _mm_loadu_ps(&sfxy[ll[2]+2*mxv]);
         d = _mm_loadu_ps(&sfxy[ll[3]+2*mxv]);
/* transpose so a,b,c,d contain next 4 fields for each of 4 particles */
         _MM_TRANSPOSE4_PS(a,b,c,d);
/*       vx = amx*sfxy[nn];   */
/*       vy = amx*sfxy[nn+1]; */
         a = _mm_mul_ps(v_amx,a);
         b = _mm_mul_ps(v_amx,b);
/*       dx += dyp*(dxp*sfxy[nn+2] + vx); */
/*       dy += dyp*(dxp*sfxy[nn+3] + vy); */
         a = _mm_mul_ps(v_dyp,_mm_add_ps(_mm_mul_ps(v_dxp,c),a));
         b = _mm_mul_ps(v_dyp,_mm_add_ps(_mm_mul_ps(v_dxp,d),b));
         v_dx = _mm_add_ps(v_dx,a);
         v_dy = _mm_add_ps(v_dy,b);
/* new velocity */
/*       dxp = ppart[j+2*nppmx+npoff]; */
/*       dyp = ppart[j+3*nppmx+npoff]; */
         v_dxp = _mm_load_ps(&ppart[j+2*nppmx+npoff]);
         v_dyp = _mm_load_ps(&ppart[j+3*nppmx+npoff]);
/*       vx = dxp + qtm*dx; */
/*       vy = dyp + qtm*dy; */
         v_vx = _mm_add_ps(v_dxp,_mm_mul_ps(v_qtm,v_dx));
         v_vy = _mm_add_ps(v_dyp,_mm_mul_ps(v_qtm,v_dy));
/* average kinetic energy */
/*       dxp += vx; */
/*       dyp += vy; */
         v_dxp = _mm_add_ps(v_dxp,v_vx);
         v_dyp = _mm_add_ps(v_dyp,v_vy);
/*       sum1 += dxp*dxp + dyp*dyp; */
         v_at = _mm_mul_ps(v_dxp,v_dxp);
         v_at = _mm_add_ps(v_at,_mm_mul_ps(v_dyp,v_dyp));
/* convert to double precision before accumulating */
         v_d = _mm_cvtps_pd(v_at);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
         v_it = _mm_srli_si128((__m128i)v_at,8);
         v_d = _mm_cvtps_pd((__m128)v_it);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
/* new position */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
         v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_dt));
         v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_dt));
/* find particles going out of bounds */
/*       mm = 0; */
         v_st = v_zero;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                              */
/*       if (dx >= edgerx) { */
/*          if (dx >= anx)   */
/*             dx -= anx;    */
/*          mm = 2;          */
/*       }                   */
         v_x = _mm_cmpge_ps(v_dx,v_edgerx);
         v_y = _mm_cmplt_ps(v_dx,v_edgelx);
         v_at = _mm_or_ps(v_x,v_y);
         v_it = _mm_srli_si128((__m128i)v_at,8);
         v_it = _mm_add_epi64((__m128i)v_at,v_it);
         _mm_storel_epi64((__m128i *)&jj[0],v_it);
/* execute if either test result is true for any particle */
         if (jj[0] != 0) {
            v_st = _mm_and_ps(v_two,v_x);
            v_x = _mm_and_ps(v_x,_mm_cmpge_ps(v_dx,v_anx));
            v_dx = _mm_sub_ps(v_dx,_mm_and_ps(v_anx,v_x));
/*          if (dx < edgelx) {  */
/*             if (dx < 0.0f) { */
/*                dx += anx;    */
/*                if (dx < anx) */
/*                   mm = 1;    */
/*                else          */
/*                   dx = 0.0;  */
/*             }                */
/*            else {            */
/*                mm = 1;       */
/*            }                 */
/*          }                   */
            v_at = _mm_and_ps(v_one,v_y);
            v_x = _mm_and_ps(v_y,_mm_cmplt_ps(v_dx,v_zero));
            v_dx = _mm_add_ps(v_dx,_mm_and_ps(v_anx,v_x));
            v_y = _mm_cmplt_ps(v_dx,v_anx);
            v_dx = _mm_and_ps(v_dx,v_y);
            v_st = _mm_add_ps(v_st,_mm_and_ps(v_at,v_y));
         }
/*       if (dy >= edgery) { */
/*          if (dy >= any)   */
/*             dy -= any;    */
/*          mm += 6;         */
/*       }                   */
         v_y = _mm_cmpge_ps(v_dy,v_edgery);
         v_x = _mm_cmplt_ps(v_dy,v_edgely);
         v_at = _mm_or_ps(v_x,v_y);
         v_it = _mm_srli_si128((__m128i)v_at,8);
         v_it = _mm_add_epi64((__m128i)v_at,v_it);
         _mm_storel_epi64((__m128i *)&jj[0],v_it);
/* execute if either test result is true for any particle */
         if (jj[0] != 0) {
            v_st = _mm_add_ps(v_st,_mm_and_ps(v_six,v_y));
            v_y = _mm_and_ps(v_y,_mm_cmpge_ps(v_dy,v_any));
            v_dy = _mm_sub_ps(v_dy,_mm_and_ps(v_any,v_y));
/*          if (dy < edgely) {  */
/*             if (dy < 0.0) {  */
/*                dy += any;    */
/*                if (dy < any) */
/*                   mm += 3;   */
/*                else          */
/*                   dy = 0.0;  */
/*             }                */
/*             else {           */
/*                mm += 3;      */
/*             }                */
/*          }                   */
            v_at = _mm_and_ps(v_three,v_x);
            v_y = _mm_and_ps(v_x,_mm_cmplt_ps(v_dy,v_zero));
            v_dy = _mm_add_ps(v_dy,_mm_and_ps(v_any,v_y));
            v_x = _mm_cmplt_ps(v_dy,v_any);
            v_dy = _mm_and_ps(v_dy,v_x);
            v_st = _mm_add_ps(v_st,_mm_and_ps(v_at,v_x));
         }
/* set new position */
/*       ppart[j+npoff] = dx;      */
/*       ppart[j+nppmx+npoff] = dy; */
         _mm_store_ps(&ppart[j+npoff],v_dx);
         _mm_store_ps(&ppart[j+nppmx+npoff],v_dy);
/* set new velocity */
/*       ppart[j+2*nppmx+npoff] = vx; */
/*       ppart[j+3*nppmx+npoff] = vy; */
         _mm_store_ps(&ppart[j+2*nppmx+npoff],v_vx);
         _mm_store_ps(&ppart[j+3*nppmx+npoff],v_vy);
/* increment counters */
/*       if (mm > 0) {                            */
/*          ncl[mm+8*k-1] += 1;                   */
/*          ih += 1;                              */
/*          if (ih <= ntmax) {                    */
/*             ihole[2*(ih+(ntmax+1)*k)] = j + 1; */
/*             ihole[1+2*(ih+(ntmax+1)*k)] = mm;  */
/*          }                                     */
/*          else {                                */
/*             nh = 1;                            */
/*          }                                     */
/*       }                                        */
         _mm_store_si128((__m128i *)ll,_mm_cvttps_epi32(v_st));
/* remove zero ist values and left shift data */
         kk = 0;
         memset((void*)lm,0,8*sizeof(int));
         for (i = 0; i < 4; i++) {
            mm = ll[i];
            if (mm > 0) {
               lm[2*kk] = j + i + 1;
               lm[1+2*kk] = mm;
               ncl[mm+8*k-1] += 1;
               kk += 1;
            }
         }
         if (kk > 0) {
            if ((ih+kk) > ntmax) {
               nh = 1;
            }
            else {
            v_it = _mm_load_si128((__m128i *)lm);
               _mm_storeu_si128((__m128i *)&ihole[2*(ih+1+(ntmax+1)*k)],v_it);
               if (kk > 2) {
                  v_it = _mm_load_si128((__m128i *)&lm[4]);
                  _mm_storeu_si128((__m128i *)&ihole[2*(ih+3+(ntmax+1)*k)],v_it);
               }
            }
            ih += kk;
         }
      }
/* loop over remaining particles in tile */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nn = 2*(nn - noff + mxv*(mm - moff));
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
         dxp = ppart[j+2*nppmx+npoff];
         dyp = ppart[j+3*nppmx+npoff];
         vx = dxp + qtm*dx;
         vy = dyp + qtm*dy;
/* average kinetic energy */
         dxp += vx;
         dyp += vy;
         sum1 += dxp*dxp + dyp*dyp;
/* new position */
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
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
/* set new velocity */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
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
/*    sum2 += sum1; */
      _mm_store_pd(&dd[0],v_sum1);
      for (j = 1; j < 2; j++) {
         dd[0] += dd[j];
      }
      sum2 += (sum1 + dd[0]);
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
void csse2gppost2lt(float ppart[], float q[], int kpic[], float qm,
                    int nppmx, int idimp, int mx, int my, int nxv,
                    int nyv, int mx1, int mxy1) {
/* for 2d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   OpenMP/vector version using guard cells
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
   requires SSE2, ppart needs to be 16 byte aligned
   nppmx needs to be a multiple of 4
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, nps, mxv;
   int i, j, k, nn, mm, it;
   float x, y, dxp, dyp, amx, amy;
   __m128i v_noff, v_moff, v_mxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qm, v_one, v_m;
   __m128 v_x, v_y, v_dxp, v_dyp, v_amx, v_amy;
   __m128 a, b, c, d;
   __attribute__((aligned(16))) unsigned int ll[4];
   __attribute__((aligned(16))) float sq[MXV*MYV];
/* __attribute__((aligned(16))) float sq[(mx+1)*(my+1)]; */
   mxv = mx + 1;
   v_mxv = _mm_set1_epi32(mxv);
   v_qm = _mm_set1_ps(qm);
   v_one = _mm_set1_ps(1.0f);
   v_m = _mm_castsi128_ps(_mm_set_epi32(-1,-1,-1,0));
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,nps,npoff,nn,mm,it,x,y,dxp,dyp,amx,amy, \
v_noff,v_moff,v_nn,v_mm,v_it,v_x,v_y,v_dxp,v_dyp,v_amx,v_amy,a,b,c,d, \
ll,sq)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm_set1_epi32(noff);
      v_moff = _mm_set1_epi32(moff);
      npp = kpic[k];
      nps = 4*(npp/4);
      npoff = idimp*nppmx*k;
/* zero out local accumulator */
/*    for (j = 0; j < mxv*(my+1); j++) { */
/*       sq[j] = 0.0f;                   */
/*    }                                  */
      memset((void*)sq,0,mxv*(my+1)*sizeof(float));
/* loop over particles in tile in groups of 4 */
      for (j = 0; j < nps; j+=4) {
/* find interpolation weights */
/*       x = ppart[j+npoff];       */
/*       y = ppart[j+nppmx+npoff]; */
         v_x = _mm_load_ps(&ppart[j+npoff]);
         v_y = _mm_load_ps(&ppart[j+nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
         v_nn = _mm_cvttps_epi32(v_x);
         v_mm = _mm_cvttps_epi32(v_y);
/*       dxp = qm*(x - (float) nn); */
         v_dxp = _mm_sub_ps(v_x,_mm_cvtepi32_ps(v_nn));
         v_dxp = _mm_mul_ps(v_dxp,v_qm);
/*       dyp = y - (float) mm; */
         v_dyp = _mm_sub_ps(v_y,_mm_cvtepi32_ps(v_mm));
/*       nn = nn - noff + mxv*(mm - moff); */
         v_nn = _mm_sub_epi32(v_nn,v_noff);
         v_mm = _mm_sub_epi32(v_mm,v_moff);
         v_it = _mm_mul_epu32(v_mxv,_mm_srli_si128(v_mm,4));
         v_mm = _mm_mul_epu32(v_mm,v_mxv);
         v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
         v_nn = _mm_add_epi32(v_nn,v_mm);
/*       amx = qm - dxp;   */
/*       amy = 1.0f - dyp; */
         v_amx = _mm_sub_ps(v_qm,v_dxp);
         v_amy = _mm_sub_ps(v_one,v_dyp);
/* calculate weights, for lower left/right, upper left/right */
         a = _mm_mul_ps(v_amx,v_amy);
         b = _mm_mul_ps(v_dxp,v_amy);
         c = _mm_mul_ps(v_amx,v_dyp);
         d = _mm_mul_ps(v_dxp,v_dyp);
         _mm_store_si128((__m128i *)ll,v_nn);
/* transpose so a,b,c,d contain the 4 weights for each of 4 particles */
         _MM_TRANSPOSE4_PS(a,b,c,d);
/* deposit charge within tile to local accumulator */
/*       x = q[nn] + amx*amy;   */
/*       y = q[nn+1] + dxp*amy; */
/*       q[nn] = x;             */
/*       q[nn+1] = y;           */
/*       nn += nxv;             */
/*       x = q[nn] + amx*dyp;   */
/*       y = q[nn+1] + dxp*dyp; */
/*       q[nn] = x;             */
/*       q[nn+1] = y;           */
/* deposit for first particle */
         mm = ll[0];
         v_x = _mm_loadl_pi(v_x,(__m64 *)&sq[mm]);
         v_x = _mm_loadh_pi(v_x,(__m64 *)&sq[mm+mxv]);
         v_x = _mm_add_ps(v_x,a);
         _mm_storel_pi((__m64 *)&sq[mm],v_x);
         _mm_storeh_pi((__m64 *)&sq[mm+mxv],v_x);
/* deposit for second particle */
         mm = ll[1];
         v_y = _mm_loadl_pi(v_y,(__m64 *)&sq[mm]);
         v_y = _mm_loadh_pi(v_y,(__m64 *)&sq[mm+mxv]);
         v_y = _mm_add_ps(v_y,b);
         _mm_storel_pi((__m64 *)&sq[mm],v_y);
         _mm_storeh_pi((__m64 *)&sq[mm+mxv],v_y);
/* deposit for third particle */
         mm = ll[2];
         v_x = _mm_loadl_pi(v_x,(__m64 *)&sq[mm]);
         v_x = _mm_loadh_pi(v_x,(__m64 *)&sq[mm+mxv]);
         v_x = _mm_add_ps(v_x,c);
         _mm_storel_pi((__m64 *)&sq[mm],v_x);
         _mm_storeh_pi((__m64 *)&sq[mm+mxv],v_x);
/* deposit for fourth particle */
         mm = ll[3];
         v_y = _mm_loadl_pi(v_y,(__m64 *)&sq[mm]);
         v_y = _mm_loadh_pi(v_y,(__m64 *)&sq[mm+mxv]);
         v_y = _mm_add_ps(v_y,d);
         _mm_storel_pi((__m64 *)&sq[mm],v_y);
         _mm_storeh_pi((__m64 *)&sq[mm+mxv],v_y);
      }
/* loop over remaining particles in tile */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
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
      nps = 4*(nn/4);
      for (j = 1; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*       for (i = 1; i < nn; i++) {                */
/*          q[i+noff+nxv*(j+moff)] += sq[i+mxv*j]; */
/*       }                                         */
         for (i = 0; i < nps; i+=4) {
            v_x = _mm_loadu_ps(&q[i+noff+nxv*(j+moff)]);
            v_y = _mm_loadu_ps(&sq[i+mxv*j]);
/* zero out first element for i = 0 */
            if (i==0)
               v_y = _mm_and_ps(v_y,v_m);
            v_x = _mm_add_ps(v_x,v_y);
            _mm_storeu_ps(&q[i+noff+nxv*(j+moff)],v_x);
         }
/* loop over remaining elements */
         it = 1 > nps ? 1 : nps;
         for (i = it; i < nn; i++) {
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
void csse2pporder2lt(float ppart[], float ppbuff[], int kpic[],
                     int ncl[], int ihole[], int idimp, int nppmx,
                     int nx, int ny, int mx, int my, int mx1, int my1,
                     int npbmx, int ntmax, int *irc) {
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
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k 
   ppbuff[k][i][n] = i co-ordinate of particle n in tile k
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
   requires SSE2, ppart, ppbuff need to be 16 byte aligned
   nppmx, npbmx need to be a multiple of 4
local data                                                            */
   int mxy1, noff, moff, npoff, npp, nps, nboff, ncoff;
   int i, j, k, ii, kx, ky, ih, nh, ist, nn, mm, isum;
   int ip, in, j1, j2, kxl, kxr, kk, kl, kr;
   float anx, any, edgelx, edgely, edgerx, edgery, dx, dy;
   __m128i v_it, v_is, v_in, v_m1, v_m2;
   __m128 v_dx, v_dy, v_st, v_at, v_x, v_y;
   __m128 v_anx, v_any, v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 v_zero, v_one, v_two, v_three, v_six;
   __attribute__((aligned(16))) unsigned int ll[8], lm[8];
   __attribute__((aligned(16))) unsigned long jj[1];
   int ks[8];
   mxy1 = mx1*my1;
   anx = (float) nx;
   any = (float) ny;
/* find and count particles leaving tiles and determine destination */
/* update ppart, ihole, ncl */
   v_anx = _mm_set1_ps(anx);
   v_any = _mm_set1_ps(any);
   v_zero = _mm_setzero_ps();
   v_one = _mm_set1_ps(1.0f);
   v_two = _mm_set1_ps(2.0f);
   v_three = _mm_set1_ps(3.0f);
   v_six = _mm_set1_ps(6.0f);
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,nps,npoff,nn,mm,ih,nh,ist,kk,dx,dy, \
edgelx,edgely,edgerx,edgery,v_it,v_edgelx,v_edgely,v_edgerx,v_edgery,\
v_dx,v_dy,v_st,v_at,v_x,v_y,jj,ll,lm)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[k];
/*    nps = 4*(npp/4); */
      nps = (npp >> 2) << 2;
      npoff = idimp*nppmx*k;
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
      noff = (ntmax+1)*k;
      v_edgelx = _mm_set1_ps(edgelx);
      v_edgely = _mm_set1_ps(edgely);
      v_edgerx = _mm_set1_ps(edgerx);
      v_edgery = _mm_set1_ps(edgery);
/* clear counters */
/*    for (j = 0; j < 8; j++) { */
/*       ncl[j+8*k] = 0;        */
/*    }                         */
      memset((void*)&ncl[8*k],0,8*sizeof(int));
/* loop over particles in tile in groups of 4 */
      for (j = 0; j < nps; j+=4) {
/*       dx = ppart[j+npoff];       */
/*       dy = ppart[j+nppmx+npoff]; */
         v_dx = _mm_load_ps(&ppart[j+npoff]);
         v_dy = _mm_load_ps(&ppart[j+nppmx+npoff]);
/* find particles going out of bounds */
/*       ist = 0; */
         v_st = v_zero;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* ist = direction particle is going                             */
/*       if (dx >= edgerx) {              */
/*          if (dx >= anx)                */
/*             ppart[j+npoff] = dx - anx; */
/*          ist = 2;                      */
/*       }                                */
         v_x = _mm_cmpge_ps(v_dx,v_edgerx);
         v_y = _mm_cmplt_ps(v_dx,v_edgelx);
         v_at = _mm_or_ps(v_x,v_y);
         v_it = _mm_srli_si128((__m128i)v_at,8);
         v_it = _mm_add_epi64((__m128i)v_at,v_it);
         _mm_storel_epi64((__m128i *)&jj[0],v_it);
/* execute if either test result is true for any particle */
         if (jj[0] != 0) {
            v_st = _mm_and_ps(v_two,v_x);
            v_x = _mm_and_ps(v_x,_mm_cmpge_ps(v_dx,v_anx));
/* write output if test result is true for any particle */
            v_it = _mm_srli_si128((__m128i)v_x,8);
            v_it = _mm_add_epi64((__m128i)v_x,v_it);
            _mm_storel_epi64((__m128i *)&jj[0],v_it);
            if (jj[0] != 0) {
               v_x = _mm_sub_ps(v_dx,_mm_and_ps(v_anx,v_x));
               _mm_store_ps(&ppart[j+npoff],v_x);
            }
/*          if (dx < edgelx) {         */
/*             if (dx < 0.0) {         */
/*                dx += anx;           */
/*                if (dx < anx)        */
/*                   ist += 1;         */
/*                else                 */
/*                   dx = 0.0;         */
/*                ppart[j+npoff] = dx; */
/*             }                       */
/*             else {                  */
/*                ist += 1;            */
/*             }                       */
/*          }                          */
            v_at = _mm_and_ps(v_one,v_y);
            v_x = _mm_and_ps(v_y,_mm_cmplt_ps(v_dx,v_zero));
/* write output if test result is true for any particle */
            v_it = _mm_srli_si128((__m128i)v_x,8);
            v_it = _mm_add_epi64((__m128i)v_x,v_it);
            _mm_storel_epi64((__m128i *)&jj[0],v_it);
            if (jj[0] != 0) {
               v_x = _mm_add_ps(v_dx,_mm_and_ps(v_anx,v_x));
               v_y = _mm_cmplt_ps(v_x,v_anx);
               v_at = _mm_and_ps(v_at,v_y);
               v_x = _mm_and_ps(v_x,v_y);
               _mm_store_ps(&ppart[j+npoff],v_x);
            }
            v_st = _mm_add_ps(v_st,v_at);
         }
/*       if (dy >= edgery) {                    */
/*          if (dy >= any)                      */
/*             ppart[j+nppmx+npoff] = dy - any; */
/*          ist += 6;                           */
/*       }                                      */
         v_y = _mm_cmpge_ps(v_dy,v_edgery);
         v_x = _mm_cmplt_ps(v_dy,v_edgely);
         v_at = _mm_or_ps(v_x,v_y);
         v_it = _mm_srli_si128((__m128i)v_at,8);
         v_it = _mm_add_epi64((__m128i)v_at,v_it);
         _mm_storel_epi64((__m128i *)&jj[0],v_it);
/* execute if either test result is true for any particle */
         if (jj[0] != 0) {
            v_st = _mm_add_ps(v_st,_mm_and_ps(v_six,v_y));
            v_y = _mm_and_ps(v_y,_mm_cmpge_ps(v_dy,v_any));
/* write output if test result is true for any particle */
            v_it = _mm_srli_si128((__m128i)v_y,8);
            v_it = _mm_add_epi64((__m128i)v_y,v_it);
            _mm_storel_epi64((__m128i *)&jj[0],v_it);
            if (jj[0] != 0) {
               v_y = _mm_sub_ps(v_dy,_mm_and_ps(v_any,v_y));
               _mm_store_ps(&ppart[j+nppmx+npoff],v_y);
            }
/*          if (dy < edgely) {               */
/*             if (dy < 0.0) {               */
/*                dy += any;                 */
/*                if (dy < any)              */
/*                   ist += 3;               */
/*                else                       */
/*                   dy = 0.0;               */
/*                ppart[j+nppmx+npoff] = dy; */
/*             }                             */
/*             else {                        */
/*                ist += 3;                  */
/*             }                             */
/*          }                                */
            v_at = _mm_and_ps(v_three,v_x);
            v_y = _mm_and_ps(v_x,_mm_cmplt_ps(v_dy,v_zero));
/* write output if test result is true for any particle */
            v_it = _mm_srli_si128((__m128i)v_y,8);
            v_it = _mm_add_epi64((__m128i)v_y,v_it);
            _mm_storel_epi64((__m128i *)&jj[0],v_it);
            if (jj[0] != 0) {
               v_y = _mm_add_ps(v_dy,_mm_and_ps(v_any,v_y));
               v_x = _mm_cmplt_ps(v_y,v_any);
               v_at = _mm_and_ps(v_at,v_x);
               v_y = _mm_and_ps(v_y,v_x);
               _mm_store_ps(&ppart[j+nppmx+npoff],v_y);
            }
            v_st = _mm_add_ps(v_st,v_at);
         }
/* increment counters */
/*       if (ist > 0) {                           */
/*          ncl[ist+8*k-1] += 1;                  */
/*          ih += 1;                              */
/*          if (ih <= ntmax) {                    */
/*             ihole[2*(ih+(ntmax+1)*k)] = j + 1; */
/*             ihole[1+2*(ih+(ntmax+1)*k)] = ist; */
/*          }                                     */
/*          else {                                */
/*             nh = 1;                            */
/*          }                                     */
/*       }                                        */
         _mm_store_si128((__m128i *)ll,_mm_cvttps_epi32(v_st));
/* remove zero ist values and left shift data */
         kk = 0;
         memset((void*)lm,0,8*sizeof(int));
         for (i = 0; i < 4; i++) {
            ist = ll[i];
            if (ist > 0) {
               lm[2*kk] = j + i + 1;
               lm[1+2*kk] = ist;
               ncl[ist+8*k-1] += 1;
               kk += 1;
            }
         }
         if (kk > 0) {
            if ((ih+kk) > ntmax) {
               nh = 1;
            }
            else {
            v_it = _mm_load_si128((__m128i *)lm);
               _mm_storeu_si128((__m128i *)&ihole[2*(ih+1+noff)],v_it);
               if (kk > 2) {
                  v_it = _mm_load_si128((__m128i *)&lm[4]);
                  _mm_storeu_si128((__m128i *)&ihole[2*(ih+3+noff)],v_it);
               }
            }
            ih += kk;
         }
      }
/* loop over remaining particles in tile */
      for (j = nps; j < npp; j++) {
         dx = ppart[j+npoff];
         dy = ppart[j+nppmx+npoff];
/* find particles going out of bounds */
         ist = 0;
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* ist = direction particle is going                             */
         if (dx >= edgerx) {
            if (dx >= anx)
               ppart[j+npoff] = dx - anx;
            ist = 2;
         }
         else if (dx < edgelx) {
            if (dx < 0.0) {
               dx += anx;
               if (dx < anx)
                  ist = 1;
               else
                  dx = 0.0;
               ppart[j+npoff] = dx;
            }
            else {
               ist = 1;
            }
         }
         if (dy >= edgery) {
            if (dy >= any)
               ppart[j+nppmx+npoff] = dy - any;
            ist += 6;
         }
         else if (dy < edgely) {
            if (dy < 0.0) {
               dy += any;
               if (dy < any)
                  ist += 3;
               else
                  dy = 0.0;
               ppart[j+nppmx+npoff] = dy;
            }
            else {
               ist += 3;
            }
         }
         if (ist > 0) {
            ncl[ist+8*k-1] += 1;
            ih += 1;
            if (ih <= ntmax) {
               ihole[2*(ih+noff)] = j + 1;
               ihole[1+2*(ih+noff)] = ist;
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
      ihole[2*noff] = ih;
   }
/* ihole overflow */
   if (*irc > 0)
      return;

/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
   v_m1 = _mm_set_epi32(0,-1,0,-1);
   v_m2 = _mm_set_epi32(0,-1,-1,0);
#pragma omp parallel for \
private(i,j,k,noff,npoff,nboff,isum,ist,nh,nps,ip,j1,ii,kk, \
v_it,v_is,v_in,lm)
   for (k = 0; k < mxy1; k++) {
      npoff = idimp*nppmx*k;
      nboff = idimp*npbmx*k;
      noff = (ntmax+1)*k;
/* find address offset for ordered ppbuff array */
      isum = 0;
/*    for (j = 0; j < 8; j++) { */
/*       ist = ncl[j+8*k];      */
/*       ncl[j+8*k] = isum;     */
/*       isum += ist;           */
/*    }                         */
/* perform exclusive prefix scan */
      v_is = _mm_setzero_si128();
      for (i = 0; i < 8; i+=4) {
         v_it = _mm_load_si128((__m128i *)&ncl[i+8*k]);
/* save last entry */
         v_in = _mm_srli_si128(v_it,12);
/* shift and add last entry from previous read */
         v_it = _mm_add_epi32(v_is,_mm_slli_si128(v_it,4));
/* first pass */
         v_is = _mm_slli_si128(_mm_and_si128(v_it,v_m1),4);
         v_it = _mm_add_epi32(v_is,v_it);
/* second pass */
         v_is = _mm_shuffle_epi32(v_it,212);
         v_is = _mm_slli_si128(_mm_and_si128(v_is,v_m2),4);
         v_it = _mm_add_epi32(v_is,v_it);
/* add last sum to next entry */
         v_is = _mm_add_epi32(v_in,_mm_srli_si128(v_it,12));
        _mm_store_si128((__m128i *)&ncl[i+8*k],v_it);
      }
      nh = ihole[2*noff];
/*    nps = 4*(nh/4); */
      nps = (nh >> 2) << 2;
      ip = 0;
/* loop over particles leaving tile in groups of 4 */
      for (j = 0; j < nps; j+=4) {
/* buffer particles that are leaving tile, in direction order */
/*       j1 = ihole[2*(j+1+noff)] - 1; */
/*       ist = ihole[1+2*(j+1+noff)];  */
         v_it = _mm_loadu_si128((__m128i *)&ihole[2*(j+1+noff)]);
         _mm_store_si128((__m128i *)lm,v_it);
         v_it = _mm_loadu_si128((__m128i *)&ihole[2*(j+3+noff)]);
         _mm_store_si128((__m128i *)&lm[4],v_it);
         for (kk = 0; kk < 4; kk++) {
            j1 = lm[2*kk] - 1;
            ist = lm[1+2*kk];
            ii = ncl[ist+8*k-1];
            if (ii < npbmx) {
               for (i = 0; i < idimp; i++) {
                  ppbuff[ii+npbmx*i+nboff]
                  = ppart[j1+nppmx*i+npoff];
               }
            }
            else {
               ip = 1;
            }
            ncl[ist+8*k-1] = ii + 1;
         }
      }
/* loop over remaining particles leaving tile */
      for (j = nps; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+noff)] - 1;
         ist = ihole[1+2*(j+1+noff)];
         ii = ncl[ist+8*k-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*i+nboff]
               = ppart[j1+nppmx*i+npoff];
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
private(i,j,k,ii,kk,npp,nps,noff,npoff,nboff,kx,ky,kl,kr,kxl,kxr,ih,nh, \
nn,mm,ncoff,ist,j1,j2,ip,in,v_it,v_is,v_in,v_x,ks,ll,lm)
   for (k = 0; k < mxy1; k++) {
      npp = kpic[k];
      npoff = idimp*nppmx*k;
      noff = (ntmax+1)*k;
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
      nh = ihole[2*noff];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      v_in = _mm_set1_epi32(1);
      for (ii = 0; ii < 8; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+8*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+8*ks[ii]] - ncoff;
/*       nps = 4*(ip/4); */
         nps = (ip >> 2) << 2;
/* loop over particles in this direction in groups of 4 */
         for (j = 0; j < nps; j+=4) {
/* insert incoming particles into holes */
/*          ih += 1;                        */
/*          if (ih <= nh) {                 */
/*             j1 = ihole[2*(ih+noff)] - 1; */
/*          }                               */
            if (ih < nh) {
               v_it = _mm_loadu_si128((__m128i *)&ihole[2*(ih+1+noff)]);
               _mm_store_si128((__m128i *)lm,_mm_sub_epi32(v_it,v_in));
            }
            if ((ih+2) < nh) {
               v_is = _mm_loadu_si128((__m128i *)&ihole[2*(ih+3+noff)]);
               _mm_store_si128((__m128i *)&lm[4],_mm_sub_epi32(v_is,v_in));
            }
/* place overflow at end of array */
/*          else {       */
/*             j1 = npp; */
/*             npp += 1; */
/*          }            */
            ih += 4;
            nn = ih - nh;
            if (nn >= 4) {
               for (kk = 0; kk < 4; kk++) {
                  lm[2*kk] = npp + kk;
               }
               npp += 4;
            }
            else if (nn > 0) {
               nn = nn < 4 ? nn : 4;
               for (kk = 4-nn; kk < 4; kk++) {
                  lm[2*kk] = npp;
                  npp += 1;
               }
            }
            for (i = 0; i < idimp; i++) {
/*             if (j1 < nppmx)                     */
/*                ppart[j1+nppmx*i+npoff]          */
/*                = ppbuff[j+ncoff+npbmx*i+nboff]; */
               v_x = _mm_loadu_ps(&ppbuff[j+ncoff+npbmx*i+nboff]);
               for (kk = 0; kk < 4; kk++) {
                  j1 = lm[2*kk];
                  if (j1 < nppmx) {
                     _mm_store_ss(&ppart[j1+nppmx*i+npoff],v_x);
                     v_x = (__m128)_mm_srli_si128((__m128i)v_x,4);
                  }
                  else {
                     ist = 1;
                  }
               }
            }
         }
/* loop over remaining particles in this direction */
         for (j = nps; j < ip; j++) {
            ih += 1;
/* insert incoming particles into holes */
            if (ih <= nh) {
               j1 = ihole[2*(ih+noff)] - 1;
            }
/* place overflow at end of array */
            else {
               j1 = npp;
               npp += 1;
            }
            if (j1 < nppmx) {
               for (i = 0; i < idimp; i++) {
                  ppart[j1+nppmx*i+npoff]
                  = ppbuff[j+ncoff+npbmx*i+nboff];
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
         ii = nh;
         ih += 1;
/* move particles from end into remaining holes */
/* holes are processed in increasing order      */
/*       nps = 4*(ip/4); */
         nps = (ip >> 2) << 2;
/* loop over particles in groups of 4 */
         for (j = 0; j < nps; j+=4) {
/*          nn = ihole[2*(ii+noff)] - 1; */
            v_it = _mm_loadu_si128((__m128i *)&ihole[2*(ii-3+noff)]);
            _mm_store_si128((__m128i *)ll,_mm_sub_epi32(v_it,v_in));
            v_is = _mm_loadu_si128((__m128i *)&ihole[2*(ii-1+noff)]);
            _mm_store_si128((__m128i *)&ll[4],_mm_sub_epi32(v_is,v_in));
/*          j2 = ihole[2*(ih+(ntmax+1)*k)] - 1; */
            v_it = _mm_loadu_si128((__m128i *)&ihole[2*(ih+noff)]);
            _mm_store_si128((__m128i *)lm,_mm_sub_epi32(v_it,v_in));
            v_is = _mm_loadu_si128((__m128i *)&ihole[2*(ih+2+noff)]);
            _mm_store_si128((__m128i *)&lm[4],_mm_sub_epi32(v_is,v_in));
/* holes with locations great than npp-ip do not need to be filled */
            in = 0;
            mm = 0;
            nn = ll[6];
            j2 = lm[0];
            for (kk = 0; kk < 4; kk++) {
               j1 = npp - (j + kk) - 1;
               ll[2*kk+1] = nn;
               lm[2*kk+1] = j2;
               if (j1==nn) {
                  in += 1; 
                  if (in < 4)
                     nn = ll[6-2*in];
               }
               else {
                  mm += 1;
                  if (mm < 4)
                     j2 = lm[2*mm];
               }
            }
            ii -= in;
            ih += mm;
/* fill holes */
            for (i = 0; i < idimp; i++) {
               ist = npp - j - 1;
               v_x = _mm_loadu_ps(&ppart[ist-3+nppmx*i+npoff]);
               v_x = _mm_shuffle_ps(v_x,v_x,27);
/*             j1 = npp - j - 1;               */
/*             if (j1==nn) {                   */
/*                ii -= 1;                     */
/*                nn = ihole[2*(ii+noff)] - 1; */
/*             }                               */
               for (kk = 0; kk < 4; kk++) {
                  j1 = ist - kk;
                  nn = ll[2*kk+1];
                  if (j1 != nn) {
/*                   ppart[j2+nppmx*i+npoff]    */
/*                   = ppart[j1+nppmx*i+npoff]; */
                     _mm_store_ss(&ppart[lm[2*kk+1]+nppmx*i+npoff],v_x);
                  }
                  v_x = (__m128)_mm_srli_si128((__m128i)v_x,4);
               }
            }
         }
/* loop over remaining particles */
         if (nps < ip) {
            nn = ihole[2*(ii+noff)] - 1;
            j2 = ihole[2*(ih+noff)] - 1;
         }
         for (j = nps; j < ip; j++) {
            j1 = npp - j - 1;
            if (j1==nn) {
               ii -= 1;
               nn = ihole[2*(ii+noff)] - 1;
            }
            else {
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*i+npoff]
                  = ppart[j1+nppmx*i+npoff];
               }
               ih += 1;
               j2 = ihole[2*(ih+(ntmax+1)*k)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[k] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2pporderf2lt(float ppart[], float ppbuff[], int kpic[],
                      int ncl[], int ihole[], int idimp, int nppmx,
                      int mx1, int my1, int npbmx, int ntmax,
                      int *irc) {
/* this subroutine sorts particles by x,y grid in tiles of mx, my
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 2D linear memory
   the algorithm has 2 steps.  first, a prefix scan of ncl is performed
   and departing particles are buffered in ppbuff in direction order.
   then we copy the incoming particles from other tiles into ppart.
   it assumes that the number, location, and destination of particles 
   leaving a tile have been previously stored in ncl and ihole by the
   cgppushf2lt procedure.
   input: all except ppbuff, irc
   output: ppart, ppbuff, kpic, ncl, irc
   ppart[k][0][n] = position x of particle n in tile k
   ppart[k][1][n] = position y of particle n in tile k 
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
   int mxy1, noff, npoff, npp, nps, nboff, ncoff;
   int i, j, k, ii, kx, ky, ih, nh, ist, nn, mm, isum;
   int ip, in, j1, j2, kxl, kxr, kk, kl, kr;
   __m128i v_it, v_is, v_in, v_m1, v_m2;
   __m128 v_x;
   __attribute__((aligned(16))) unsigned int ll[8], lm[8];
   int ks[8];
   mxy1 = mx1*my1;
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
   v_m1 = _mm_set_epi32(0,-1,0,-1);
   v_m2 = _mm_set_epi32(0,-1,-1,0);
#pragma omp parallel for \
private(i,j,k,noff,npoff,nboff,isum,ist,nh,nps,ip,j1,ii,kk, \
v_it,v_is,v_in,lm)
   for (k = 0; k < mxy1; k++) {
      npoff = idimp*nppmx*k;
      nboff = idimp*npbmx*k;
      noff = (ntmax+1)*k;
/* find address offset for ordered ppbuff array */
      isum = 0;
/*    for (j = 0; j < 8; j++) { */
/*       ist = ncl[j+8*k];      */
/*       ncl[j+8*k] = isum;     */
/*       isum += ist;           */
/*    }                         */
/* perform exclusive prefix scan */
      v_is = _mm_setzero_si128();
      for (i = 0; i < 8; i+=4) {
         v_it = _mm_load_si128((__m128i *)&ncl[i+8*k]);
/* save last entry */
         v_in = _mm_srli_si128(v_it,12);
/* shift and add last entry from previous read */
         v_it = _mm_add_epi32(v_is,_mm_slli_si128(v_it,4));
/* first pass */
         v_is = _mm_slli_si128(_mm_and_si128(v_it,v_m1),4);
         v_it = _mm_add_epi32(v_is,v_it);
/* second pass */
         v_is = _mm_shuffle_epi32(v_it,212);
         v_is = _mm_slli_si128(_mm_and_si128(v_is,v_m2),4);
         v_it = _mm_add_epi32(v_is,v_it);
/* add last sum to next entry */
         v_is = _mm_add_epi32(v_in,_mm_srli_si128(v_it,12));
        _mm_store_si128((__m128i *)&ncl[i+8*k],v_it);
      }
      nh = ihole[2*noff];
/*    nps = 4*(nh/4); */
      nps = (nh >> 2) << 2;
      ip = 0;
/* loop over particles leaving tile in groups of 4 */
      for (j = 0; j < nps; j+=4) {
/* buffer particles that are leaving tile, in direction order */
/*       j1 = ihole[2*(j+1+noff)] - 1; */
/*       ist = ihole[1+2*(j+1+noff)];  */
         v_it = _mm_loadu_si128((__m128i *)&ihole[2*(j+1+noff)]);
         _mm_store_si128((__m128i *)lm,v_it);
         v_it = _mm_loadu_si128((__m128i *)&ihole[2*(j+3+noff)]);
         _mm_store_si128((__m128i *)&lm[4],v_it);
         for (kk = 0; kk < 4; kk++) {
            j1 = lm[2*kk] - 1;
            ist = lm[1+2*kk];
            ii = ncl[ist+8*k-1];
            if (ii < npbmx) {
               for (i = 0; i < idimp; i++) {
                  ppbuff[ii+npbmx*i+nboff]
                  = ppart[j1+nppmx*i+npoff];
               }
            }
            else {
               ip = 1;
            }
            ncl[ist+8*k-1] = ii + 1;
         }
      }
/* loop over remaining particles leaving tile */
      for (j = nps; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+noff)] - 1;
         ist = ihole[1+2*(j+1+noff)];
         ii = ncl[ist+8*k-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*i+nboff]
               = ppart[j1+nppmx*i+npoff];
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
private(i,j,k,ii,kk,npp,nps,noff,npoff,nboff,kx,ky,kl,kr,kxl,kxr,ih,nh, \
nn,mm,ncoff,ist,j1,j2,ip,in,v_it,v_is,v_in,v_x,ks,ll,lm)
   for (k = 0; k < mxy1; k++) {
      npp = kpic[k];
      npoff = idimp*nppmx*k;
      noff = (ntmax+1)*k;
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
      nh = ihole[2*noff];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      v_in = _mm_set1_epi32(1);
      for (ii = 0; ii < 8; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+8*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+8*ks[ii]] - ncoff;
/*       nps = 4*(ip/4); */
         nps = (ip >> 2) << 2;
/* loop over particles in this direction in groups of 4 */
         for (j = 0; j < nps; j+=4) {
/* insert incoming particles into holes */
/*          ih += 1;                        */
/*          if (ih <= nh) {                 */
/*             j1 = ihole[2*(ih+noff)] - 1; */
/*          }                               */
            if (ih < nh) {
               v_it = _mm_loadu_si128((__m128i *)&ihole[2*(ih+1+noff)]);
               _mm_store_si128((__m128i *)lm,_mm_sub_epi32(v_it,v_in));
            }
            if ((ih+2) < nh) {
               v_is = _mm_loadu_si128((__m128i *)&ihole[2*(ih+3+noff)]);
               _mm_store_si128((__m128i *)&lm[4],_mm_sub_epi32(v_is,v_in));
            }
/* place overflow at end of array */
/*          else {       */
/*             j1 = npp; */
/*             npp += 1; */
/*          }            */
            ih += 4;
            nn = ih - nh;
            if (nn >= 4) {
               for (kk = 0; kk < 4; kk++) {
                  lm[2*kk] = npp + kk;
               }
               npp += 4;
            }
            else if (nn > 0) {
               nn = nn < 4 ? nn : 4;
               for (kk = 4-nn; kk < 4; kk++) {
                  lm[2*kk] = npp;
                  npp += 1;
               }
            }
            for (i = 0; i < idimp; i++) {
/*             if (j1 < nppmx)                     */
/*                ppart[j1+nppmx*i+npoff]          */
/*                = ppbuff[j+ncoff+npbmx*i+nboff]; */
               v_x = _mm_loadu_ps(&ppbuff[j+ncoff+npbmx*i+nboff]);
               for (kk = 0; kk < 4; kk++) {
                  j1 = lm[2*kk];
                  if (j1 < nppmx) {
                     _mm_store_ss(&ppart[j1+nppmx*i+npoff],v_x);
                     v_x = (__m128)_mm_srli_si128((__m128i)v_x,4);
                  }
                  else {
                     ist = 1;
                  }
               }
            }
         }
/* loop over remaining particles in this direction */
         for (j = nps; j < ip; j++) {
            ih += 1;
/* insert incoming particles into holes */
            if (ih <= nh) {
               j1 = ihole[2*(ih+noff)] - 1;
            }
/* place overflow at end of array */
            else {
               j1 = npp;
               npp += 1;
            }
            if (j1 < nppmx) {
               for (i = 0; i < idimp; i++) {
                  ppart[j1+nppmx*i+npoff]
                  = ppbuff[j+ncoff+npbmx*i+nboff];
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
         ii = nh;
         ih += 1;
/* move particles from end into remaining holes */
/* holes are processed in increasing order      */
/*       nps = 4*(ip/4); */
         nps = (ip >> 2) << 2;
/* loop over particles in groups of 4 */
         for (j = 0; j < nps; j+=4) {
/*          nn = ihole[2*(ii+noff)] - 1; */
            v_it = _mm_loadu_si128((__m128i *)&ihole[2*(ii-3+noff)]);
            _mm_store_si128((__m128i *)ll,_mm_sub_epi32(v_it,v_in));
            v_is = _mm_loadu_si128((__m128i *)&ihole[2*(ii-1+noff)]);
            _mm_store_si128((__m128i *)&ll[4],_mm_sub_epi32(v_is,v_in));
/*          j2 = ihole[2*(ih+(ntmax+1)*k)] - 1; */
            v_it = _mm_loadu_si128((__m128i *)&ihole[2*(ih+noff)]);
            _mm_store_si128((__m128i *)lm,_mm_sub_epi32(v_it,v_in));
            v_is = _mm_loadu_si128((__m128i *)&ihole[2*(ih+2+noff)]);
            _mm_store_si128((__m128i *)&lm[4],_mm_sub_epi32(v_is,v_in));
/* holes with locations great than npp-ip do not need to be filled */
            in = 0;
            mm = 0;
            nn = ll[6];
            j2 = lm[0];
            for (kk = 0; kk < 4; kk++) {
               j1 = npp - (j + kk) - 1;
               ll[2*kk+1] = nn;
               lm[2*kk+1] = j2;
               if (j1==nn) {
                  in += 1; 
                  if (in < 4)
                     nn = ll[6-2*in];
               }
               else {
                  mm += 1;
                  if (mm < 4)
                     j2 = lm[2*mm];
               }
            }
            ii -= in;
            ih += mm;
/* fill holes */
            for (i = 0; i < idimp; i++) {
               ist = npp - j - 1;
               v_x = _mm_loadu_ps(&ppart[ist-3+nppmx*i+npoff]);
               v_x = _mm_shuffle_ps(v_x,v_x,27);
/*             j1 = npp - j - 1;               */
/*             if (j1==nn) {                   */
/*                ii -= 1;                     */
/*                nn = ihole[2*(ii+noff)] - 1; */
/*             }                               */
               for (kk = 0; kk < 4; kk++) {
                  j1 = ist - kk;
                  nn = ll[2*kk+1];
                  if (j1 != nn) {
/*                   ppart[j2+nppmx*i+npoff]    */
/*                   = ppart[j1+nppmx*i+npoff]; */
                     _mm_store_ss(&ppart[lm[2*kk+1]+nppmx*i+npoff],v_x);
                  }
                  v_x = (__m128)_mm_srli_si128((__m128i)v_x,4);
               }
            }
         }
/* loop over remaining particles */
         if (nps < ip) {
            nn = ihole[2*(ii+noff)] - 1;
            j2 = ihole[2*(ih+noff)] - 1;
         }
         for (j = nps; j < ip; j++) {
            j1 = npp - j - 1;
            if (j1==nn) {
               ii -= 1;
               nn = ihole[2*(ii+noff)] - 1;
            }
            else {
               for (i = 0; i < idimp; i++) {
                  ppart[j2+nppmx*i+npoff]
                  = ppart[j1+nppmx*i+npoff];
               }
               ih += 1;
               j2 = ihole[2*(ih+(ntmax+1)*k)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[k] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2cguard2l(float fxy[], int nx, int ny, int nxe, int nye) {
/* replicate extended periodic vector field fxy
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   requires SSE2, fxy needs to be 16 byte aligned
   nxe*ny needs to be a multiple of 2
local data                                                 */
   int j, k, nxs;
   nxs = 2*(nx/2);
/* copy edges of extended field */
   for (k = 0; k < ny; k++) {
      fxy[2*nx+2*nxe*k] = fxy[2*nxe*k];
      fxy[1+2*nx+2*nxe*k] = fxy[1+2*nxe*k];
   }
/* vector loop over elements in blocks of 2 */
   for (j = 0; j < nxs; j+=2) {
      _mm_store_ps(&fxy[2*j+2*nxe*ny],_mm_load_ps(&fxy[2*j]));
   }
/* loop over remaining elements */
   for (j = nxs; j < nx; j++) {
      fxy[2*j+2*nxe*ny] = fxy[2*j];
      fxy[1+2*j+2*nxe*ny] = fxy[1+2*j];
   }
   fxy[2*nx+2*nxe*ny] = fxy[0];
   fxy[1+2*nx+2*nxe*ny] = fxy[1];
   return;
}

/*--------------------------------------------------------------------*/
void csse2aguard2l(float q[], int nx, int ny, int nxe, int nye) {
/* accumulate extended periodic scalar field q
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   requires SSE2, q needs to be 16 byte aligned
   nxe*ny needs to be a multiple of 4
local data                                                 */
   int j, k, nxs;
   __m128 v_q;
   nxs = 4*(nx/4);
/* accumulate edges of extended field */
   for (k = 0; k < ny; k++) {
      q[nxe*k] += q[nx+nxe*k];
      q[nx+nxe*k] = 0.0;
   }
/* vector loop over elements in blocks of 4 */
   for (j = 0; j < nxs; j+=4) {
      v_q = _mm_add_ps(_mm_load_ps(&q[j]),_mm_load_ps(&q[j+nxe*ny]));
      _mm_store_ps(&q[j],v_q);
      _mm_store_ps(&q[j+nxe*ny],_mm_setzero_ps());
   }
/* loop over remaining elements */
   for (j = nxs; j < nx; j++) {
      q[j] += q[j+nxe*ny];
      q[j+nxe*ny] = 0.0;
   }
   q[0] += q[nx+nxe*ny];
   q[nx+nxe*ny] = 0.0;
   return;
}

/*--------------------------------------------------------------------*/
void csse2mpois22(float complex q[], float complex fxy[], int isign,
                  float complex ffc[], float ax, float ay, float affp,
                  float *we, int nx, int ny, int nxvh, int nyv,
                  int nxhd, int nyhd) {
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
   g[ky][kx] = (affp/(kx**2+ky**2))*s(kx,ky),
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
   requires SSE2, q, fxy, ffc need to be 16 byte aligned
   nxhd, nxvh need to be a multiple of 2
local data                                                 */
   int nxh, nyh, nxhs, j, k, k1, kk, kj, it;
   float dnx, dny, dkx, dky, at1, at2, at3, at4;
   float complex zero, zt1, zt2;
   double wp, sum1;
   __m128i v_j, v_it;
   __m128 v_dnx, v_dny, v_dky, v_at1, v_at2, v_at3, v_at4;
   __m128 v_zero, v_m, v_zt1, v_zt2, v_zt3, v_zt4;
   __m128d v_wp, v_d;
   __attribute__((aligned(16))) double dd[2];
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nxhs = 2*(nxh/2);
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
   v_j = _mm_set_epi32(1,1,0,0);
   v_dnx = _mm_set1_ps(dnx);
   v_dny = _mm_set1_ps(dny);
   v_zero = _mm_set1_ps(0.0f);
   v_m = _mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);
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
#pragma omp parallel for \
private(j,k,k1,kk,kj,dky,at1,at2,at3,zt1,zt2,wp,v_it,v_dky,v_at1, \
v_at2,v_at3,v_at4,v_zt1,v_zt2,v_zt3,v_zt4,v_wp,v_d,dd) \
reduction(+:sum1)
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
      k1 = ny - k;
      dky = dny*(float) k;
      v_dky = _mm_mul_ps(v_dny,_mm_cvtepi32_ps(_mm_set1_epi32(k)));
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      wp = 0.0;
      v_wp = _mm_set1_pd(0.0);
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       at1 = crealf(ffc[j+kk])*cimagf(ffc[j+kk]); */
         v_at1 = _mm_load_ps((float *)&ffc[j+kk]);
         v_at1 = _mm_mul_ps(v_at1,_mm_shuffle_ps(v_at1,v_at1,177));
/*       at2 = at1*dnx*(float) j; */
         v_it = _mm_add_epi32(_mm_set1_epi32(j),v_j);
         v_at2 = _mm_mul_ps(v_dnx,_mm_cvtepi32_ps(v_it));
         v_at2 = _mm_mul_ps(v_at1,v_at2);
/*       at3 = dky*at1; */
         v_at3 = _mm_mul_ps(v_dky,v_at1);
/*       zt1 = cimagf(q[j+kj]) - crealf(q[j+kj])*_Complex_I; */
         v_zt1 = _mm_load_ps((float *)&q[j+kj]);
         v_zt1 = _mm_mul_ps(v_zt1,v_m);
         v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*       zt2 = cimagf(q[j+k1]) - crealf(q[j+k1])*_Complex_I; */
         v_zt2 = _mm_load_ps((float *)&q[j+k1]);
         v_zt2 = _mm_mul_ps(v_zt2,v_m);
         v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
/* zero out kx = 0 mode */
         if (j==0) {
            v_at4 = _mm_castsi128_ps(_mm_set_epi32(-1,-1,0,0));
            v_zt1 = _mm_and_ps(v_zt1,v_at4);
            v_zt2 = _mm_and_ps(v_zt2,v_at4);
         }
/*       fxy[2*j+2*kj] = at2*zt1;   */
/*       fxy[1+2*j+2*kj] = at3*zt1; */
         v_at4 = _mm_mul_ps(v_at2,v_zt1);
         v_zt4 = _mm_mul_ps(v_at3,v_zt1);
         v_at4 = _mm_mul_ps(v_at2,v_zt1);
         v_zt4 = _mm_mul_ps(v_at3,v_zt1);
/* reorder write */
         v_zt3 = _mm_shuffle_ps(v_at4,v_zt4,68);
         v_zt4 = _mm_shuffle_ps(v_at4,v_zt4,238);
         _mm_store_ps((float *)&fxy[2*(j+kj)],v_zt3);
         _mm_store_ps((float *)&fxy[2*(j+1+kj)],v_zt4);
/*       fxy[2*j+2*k1] = at2*zt2;    */
/*       fxy[1+2*j+2*k1] = -at3*zt2; */
         v_at4 = _mm_mul_ps(v_at2,v_zt2);
         v_zt4 = _mm_sub_ps(v_zero,_mm_mul_ps(v_at3,v_zt2));
/* reorder write */
         v_zt3 = _mm_shuffle_ps(v_at4,v_zt4,68);
         v_zt4 = _mm_shuffle_ps(v_at4,v_zt4,238);
         _mm_store_ps((float *)&fxy[2*(j+k1)],v_zt3);
         _mm_store_ps((float *)&fxy[2*(j+1+k1)],v_zt4);
/*       wp += at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1])); */
         v_at4 = _mm_mul_ps(v_zt1,v_zt1);
         v_at4 = _mm_add_ps(v_at4,_mm_mul_ps(v_zt2,v_zt2));
         v_at4 = _mm_mul_ps(v_at1,v_at4);
/* convert to double precision before accumulating */
         v_d = _mm_cvtps_pd(v_at4);
         v_wp = _mm_add_pd(v_wp,v_d);
         v_it = _mm_srli_si128((__m128i)v_at4,8);
         v_d = _mm_cvtps_pd((__m128)v_it);
         v_wp = _mm_add_pd(v_wp,v_d);
      }
/* loop over remaining elements */
      it = 1 > nxhs ? 1 : nxhs;
      for (j = it; j < nxh; j++) {
         at1 = crealf(ffc[j+kk])*cimagf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         zt1 = cimagf(q[j+kj]) - crealf(q[j+kj])*_Complex_I;
         zt2 = cimagf(q[j+k1]) - crealf(q[j+k1])*_Complex_I;
         fxy[2*j+2*kj] = at2*zt1;
         fxy[1+2*j+2*kj] = at3*zt1;
         fxy[2*j+2*k1] = at2*zt2;
         fxy[1+2*j+2*k1] = -at3*zt2;
         wp += at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1]));
      }
/*    sum1 += wp; */
      _mm_store_pd(&dd[0],v_wp);
      for (j = 1; j < 2; j++) {
         dd[0] += dd[j];
      }
      sum1 += (wp + dd[0]);
   }
   wp = 0.0;
   v_wp = _mm_set1_pd(0.0);
/* mode numbers kx = 0, nx/2 */
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
      wp += at1*(q[kj]*conjf(q[kj]));
   }
/* mode numbers ky = 0, ny/2 */
   k1 = 2*nxvh*nyh;
/* vector loop over elements in blocks of 2 */
   for (j = 0; j < nxhs; j+=2) {
/*    at1 = crealf(ffc[j])*cimagf(ffc[j]); */
      v_at1 = _mm_load_ps((float *)&ffc[j]);
      v_at1 = _mm_mul_ps(v_at1,_mm_shuffle_ps(v_at1,v_at1,177));
/*    at2 = at1*dnx*(float) j; */
      v_it = _mm_add_epi32(_mm_set1_epi32(j),v_j);
      v_at2 = _mm_mul_ps(v_dnx,_mm_cvtepi32_ps(v_it));
      v_at2 = _mm_mul_ps(v_at1,v_at2);
/*    zt1 = cimagf(q[j]) - crealf(q[j])*_Complex_I; */
      v_zt1 = _mm_load_ps((float *)&q[j]);
      v_zt1 = _mm_mul_ps(v_zt1,v_m);
      v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/* zero out kx = 0 mode */
      if (j==0) {
         v_at4 = _mm_castsi128_ps(_mm_set_epi32(-1,-1,0,0));
         v_zt1 = _mm_and_ps(v_zt1,v_at4);
      }
/*    fxy[2*j] = at2*zt1; */
/*    fxy[1+2*j] = zero;  */
      v_at4 = _mm_mul_ps(v_at2,v_zt1);
/* reorder write */
      v_zt3 = _mm_shuffle_ps(v_at4,v_zero,68);
      v_zt4 = _mm_shuffle_ps(v_at4,v_zero,238);
      _mm_store_ps((float *)&fxy[2*j],v_zt3);
      _mm_store_ps((float *)&fxy[2*j+2],v_zt4);
/*    fxy[2*j+k1] = zero;   */
/*    fxy[1+2*j+k1] = zero; */
      _mm_store_ps((float *)&fxy[2*j+k1],v_zero);
      _mm_store_ps((float *)&fxy[2*j+2+k1],v_zero);
/*    wp += at1*(q[j]*conjf(q[j])); */
      v_at4 = _mm_mul_ps(v_at1,_mm_mul_ps(v_zt1,v_zt1));
/* convert to double precision before accumulating */
      v_d = _mm_cvtps_pd(v_at4);
      v_wp = _mm_add_pd(v_wp,v_d);
      v_it = _mm_srli_si128((__m128i)v_at4,8);
      v_d = _mm_cvtps_pd((__m128)v_it);
      v_wp = _mm_add_pd(v_wp,v_d);
   }
/* loop over remaining elements */
   it = 1 > nxhs ? 1 : nxhs;
   for (j = it; j < nxh; j++) {
      at1 = crealf(ffc[j])*cimagf(ffc[j]);
      at2 = at1*dnx*(float) j;  
      zt1 = cimagf(q[j]) - crealf(q[j])*_Complex_I;
      fxy[2*j] = at2*zt1;
      fxy[1+2*j] = zero;
      fxy[2*j+k1] = zero;
      fxy[1+2*j+k1] = zero;
      wp += at1*(q[j]*conjf(q[j]));
   }
   fxy[0] = zero;
   fxy[1] = zero;
   fxy[k1] = zero;
   fxy[1+k1] = zero;
   sum1 += wp;
/* *we = sum1*(float) (nx*ny); */
   _mm_store_pd(&dd[0],v_wp);
   for (j = 1; j < 2; j++) {
      dd[0] += dd[j];
   }
   *we = (sum1 + dd[0])*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void csse2fft2rmxx(float complex f[], int isign, int mixup[],
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
   requires SSE2, f needs to be 16 byte aligned
   nxhd need to be a multiple of 2
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, j1, k1, k2, ns, ns2, km, kmr, nrxb, joff;
   int nss, nxhhs, it;
   float ani;
   float complex t1, t2, t3;
   __m128 v_m, v_n, v_t1, v_t2, v_t3, v_t4, v_ani;
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
   nxhhs = 2*(nxhh/2);
   v_m = _mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);
   v_n = _mm_set_ps(-1.0f,1.0f,-1.0f,1.0f);
   v_t1 = _mm_setzero_ps();
   v_t2 = _mm_setzero_ps();
   v_t3 = _mm_setzero_ps();
   if (isign > 0)
      goto L70;
/* inverse fourier transform */
   nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,nss,km,kmr,k1,k2,j1,joff,it,ani,t1,t2,t3,v_t1, \
v_t2,v_t3,v_t4)
   for (i = nyi-1; i < nyt; i++) {
      joff = nxhd*i;
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
/*          t1 = f[j1+joff]; */
            v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&f[j1+joff]);
/*          f[j1+joff] = f[j+joff]; */
            v_t2 = _mm_loadl_pi(v_t2,(__m64 *)&f[j+joff]);
            _mm_storel_pi((__m64 *)&f[j1+joff],v_t2);
/*          f[j+joff] = t1; */
            _mm_storel_pi((__m64 *)&f[j+joff],v_t1);
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
            nss = 2*(ns/2);
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nss; j+=2) {
/*             t1 = sct[kmr*j]; */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_loadh_pi(v_t1,(__m64 *)&sct[kmr*j+kmr]);
/*             t2 = t1*f[j+k2+joff]; */
               v_t2 = _mm_load_ps((float *)&f[j+k2+joff]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             f[j+k2+joff] = f[j+k1+joff] - t2; */
               v_t3 = _mm_load_ps((float *)&f[j+k1+joff]);
               _mm_store_ps((float *)&f[j+k2+joff],_mm_sub_ps(v_t3,v_t2));
/*             f[j+k1+joff] += t2; */
               _mm_store_ps((float *)&f[j+k1+joff],_mm_add_ps(v_t3,v_t2));
            }
/* loop over remaining elements */
            for (j = nss; j < ns; j++) {
               t1 = sct[kmr*j];
               t2 = t1*f[j+k2+joff];
               f[j+k2+joff] = f[j+k1+joff] - t2;
               f[j+k1+joff] += t2;
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxy/nx;
      ani = 0.5/(((float) nx)*((float) ny));
      v_ani = _mm_set1_ps(ani);
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhhs; j+=2) {
/*       t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I; */
         v_t3 = _mm_loadl_pi(v_t3,(__m64 *)&sct[kmr*j]);
         v_t3 = _mm_loadh_pi(v_t3,(__m64 *)&sct[kmr*j+kmr]);
         v_t3 = _mm_mul_ps(v_t3,v_m);
         v_t3 = _mm_shuffle_ps(v_t3,v_t3,177);
/*       t2 = conjf(f[nxh-j+joff]); */
         if (j==0) {
            v_t2 = _mm_setzero_ps();
         }
         else {
            v_t2 = _mm_loadl_pi(v_t2,(__m64 *)&f[nxh-j+joff]);
         }
         v_t2 = _mm_loadh_pi(v_t2,(__m64 *)&f[nxh-j-1+joff]);
         v_t2 = _mm_mul_ps(v_t2,v_n);
/*       t1 = f[j+joff] + t2; */
         v_t4 = _mm_load_ps((float *)&f[j+joff]);
         v_t1 = _mm_add_ps(v_t4,v_t2);
/*       t2 = (f[j+joff] - t2)*t3; */
         v_t2 = _mm_sub_ps(v_t4,v_t2);
         v_t4 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,160));
         v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
         v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,245));
         v_t2 = _mm_add_ps(v_t4,_mm_mul_ps(v_t2,v_m));
/*       f[j+joff] = ani*(t1 + t2);          */
/*       f[nxh-j+joff] = ani*conjf(t1 - t2); */
         v_t3 = _mm_mul_ps(v_ani,_mm_add_ps(v_t1,v_t2));
         v_t4 = _mm_mul_ps(v_ani,_mm_mul_ps(_mm_sub_ps(v_t1,v_t2),v_n));
         if (j==0) {
            _mm_storeh_pi((__m64 *)&f[joff+1],v_t3);
            _mm_storeh_pi((__m64 *)&f[nxh-1+joff],v_t4);
         }
         else {
            _mm_store_ps((float *)&f[j+joff],v_t3);
            _mm_storel_pi((__m64 *)&f[nxh-j+joff],v_t4);
            _mm_storeh_pi((__m64 *)&f[nxh-j-1+joff],v_t4);
         }
      }
/* loop over remaining elements */
      it = 1 > nxhhs ? 1 : nxhhs;
      for (j = it; j < nxhh; j++) {
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
private(i,j,k,l,ns,ns2,nss,km,kmr,k1,k2,j1,joff,it,t1,t2,t3,v_t1,v_t2, \
v_t3,v_t4)
   for (i = nyi-1; i < nyt; i++) {
      joff = nxhd*i;
/* scramble coefficients */
      kmr = nxy/nx;
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhhs; j+=2) {
/*       t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I; */
         v_t3 = _mm_loadl_pi(v_t3,(__m64 *)&sct[kmr*j]);
         v_t3 = _mm_loadh_pi(v_t3,(__m64 *)&sct[kmr*j+kmr]);
         v_t3 = _mm_shuffle_ps(v_t3,v_t3,177);
/*       t2 = conjf(f[nxh-j+joff]); */
         if (j==0) {
            v_t2 = _mm_setzero_ps();
         }
         else {
            v_t2 = _mm_loadl_pi(v_t2,(__m64 *)&f[nxh-j+joff]);
         }
         v_t2 = _mm_loadh_pi(v_t2,(__m64 *)&f[nxh-j-1+joff]);
         v_t2 = _mm_mul_ps(v_t2,v_n);
/*       t1 = f[j+joff] + t2; */
         v_t4 = _mm_load_ps((float *)&f[j+joff]);
         v_t1 = _mm_add_ps(v_t4,v_t2);
/*       t2 = (f[j+joff] - t2)*t3; */
         v_t2 = _mm_sub_ps(v_t4,v_t2);
         v_t4 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,160));
         v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
         v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,245));
         v_t2 = _mm_add_ps(v_t4,_mm_mul_ps(v_t2,v_m));
/*       f[j+joff] = t1 + t2;            */
/*       f[nxh-j+joff] = conjf(t1 - t2); */
         v_t3 = _mm_add_ps(v_t1,v_t2);
         v_t4 = _mm_mul_ps(_mm_sub_ps(v_t1,v_t2),v_n);
         if (j==0) {
            _mm_storeh_pi((__m64 *)&f[joff+1],v_t3);
            _mm_storeh_pi((__m64 *)&f[nxh-1+joff],v_t4);
         }
         else {
            _mm_store_ps((float *)&f[j+joff],v_t3);
            _mm_storel_pi((__m64 *)&f[nxh-j+joff],v_t4);
            _mm_storeh_pi((__m64 *)&f[nxh-j-1+joff],v_t4);
         }
      }
/* loop over remaining elements */
      it = 1 > nxhhs ? 1 : nxhhs;
      for (j = it; j < nxhh; j++) {
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
/*          t1 = f[j1+joff]; */
            v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&f[j1+joff]);
/*          f[j1+joff] = f[j+joff]; */
            v_t2 = _mm_loadl_pi(v_t2,(__m64 *)&f[j+joff]);
            _mm_storel_pi((__m64 *)&f[j1+joff],v_t2);
/*          f[j+joff] = t1; */
            _mm_storel_pi((__m64 *)&f[j+joff],v_t1);
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
            nss = 2*(ns/2);
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nss; j+=2) {
/*             t1 = conjf(sct[kmr*j]); */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_loadh_pi(v_t1,(__m64 *)&sct[kmr*j+kmr]);
               v_t1 = _mm_mul_ps(v_t1,v_n);
/*             t2 = t1*f[j+k2+joff]; */
               v_t2 = _mm_load_ps((float *)&f[j+k2+joff]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             f[j+k2+joff] = f[j+k1+joff] - t2; */
               v_t3 = _mm_load_ps((float *)&f[j+k1+joff]);
               _mm_store_ps((float *)&f[j+k2+joff],_mm_sub_ps(v_t3,v_t2));
/*             f[j+k1+joff] += t2; */
               _mm_store_ps((float *)&f[j+k1+joff],_mm_add_ps(v_t3,v_t2));
            }
/* loop over remaining elements */
            for (j = nss; j < ns; j++) {
               t1 = conjf(sct[kmr*j]);
               t2 = t1*f[j+k2+joff];
               f[j+k2+joff] = f[j+k1+joff] - t2;
               f[j+k1+joff] += t2;
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2fft2rmxy(float complex f[], int isign, int mixup[],
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
   requires SSE2, f needs to be 16 byte aligned
   nxhd needs to be a multiple of 2, and nxi needs to be odd
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, nryb, koff;
   int nss;
   float complex t1, t2;
   __m128 v_m, v_n, v_t1, v_t2, v_t3;
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
   v_m = _mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);
   v_n = _mm_set_ps(-1.0f,1.0f,-1.0f,1.0f);
   if (isign > 0)
      goto L70;
/* inverse fourier transform */
   nryb = nxhy/ny;
   nry = nxy/ny;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,nss,km,kmr,k1,k2,j1,j2,koff,t1,t2,v_t1,v_t2,v_t3)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nxhd*k1;
/*          t1 = f[i+k1]; */
            v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&f[i+k1]);
/*          f[i+k1] = f[i+koff]; */
            v_t2 = _mm_loadl_pi(v_t2,(__m64 *)&f[i+koff]);
            _mm_storel_pi((__m64 *)&f[i+k1],v_t2);
/*          f[i+koff] = t1; */
            _mm_storel_pi((__m64 *)&f[i+koff],v_t1);
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
            nss = 2*(ns/2);
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nss; j+=2) {
               j1 = nxhd*(j + k1);
               j2 = nxhd*(j + k2);
/*             t1 = sct[kmr*j]; */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_loadh_pi(v_t1,(__m64 *)&sct[kmr*j+kmr]);
/*             t2 = t1*f[i+j2]; */
               v_t2 = _mm_loadl_pi(v_t2,(__m64 *)&f[i+j2]);
               v_t2 = _mm_loadh_pi(v_t2,(__m64 *)&f[i+j2+nxhd]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             f[i+j2] = f[i+j1] - t2; */
               v_t3 = _mm_loadl_pi(v_t3,(__m64 *)&f[i+j1]);
               v_t3 = _mm_loadh_pi(v_t3,(__m64 *)&f[i+j1+nxhd]);
               v_t1 = _mm_sub_ps(v_t3,v_t2);
               _mm_storel_pi((__m64 *)&f[i+j2],v_t1);
               _mm_storeh_pi((__m64 *)&f[i+j2+nxhd],v_t1);
/*             f[i+j1] += t2; */
               v_t1 = _mm_add_ps(v_t3,v_t2);
               _mm_storel_pi((__m64 *)&f[i+j1],v_t1);
               _mm_storeh_pi((__m64 *)&f[i+j1+nxhd],v_t1);
            }
/* loop over remaining elements */
            for (j = nss; j < ns; j++) {
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
private(i,j,k,l,ns,ns2,nss,km,kmr,k1,k2,j1,j2,koff,t1,t2,v_t1,v_t2,v_t3)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nxhd*k1;
/*          t1 = f[i+k1]; */
            v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&f[i+k1]);
/*          f[i+k1] = f[i+koff]; */
            v_t2 = _mm_loadl_pi(v_t2,(__m64 *)&f[i+koff]);
            _mm_storel_pi((__m64 *)&f[i+k1],v_t2);
/*          f[i+koff] = t1; */
            _mm_storel_pi((__m64 *)&f[i+koff],v_t1);
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
            nss = 2*(ns/2);
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nss; j+=2) {
               j1 = nxhd*(j + k1);
               j2 = nxhd*(j + k2);
/*             t1 = conjf(sct[kmr*j]); */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_loadh_pi(v_t1,(__m64 *)&sct[kmr*j+kmr]);
               v_t1 = _mm_mul_ps(v_t1,v_n);
/*             t2 = t1*f[i+j2]; */
               v_t2 = _mm_loadl_pi(v_t2,(__m64 *)&f[i+j2]);
               v_t2 = _mm_loadh_pi(v_t2,(__m64 *)&f[i+j2+nxhd]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             f[i+j2] = f[i+j1] - t2; */
               v_t3 = _mm_loadl_pi(v_t3,(__m64 *)&f[i+j1]);
               v_t3 = _mm_loadh_pi(v_t3,(__m64 *)&f[i+j1+nxhd]);
               v_t1 = _mm_sub_ps(v_t3,v_t2);
               _mm_storel_pi((__m64 *)&f[i+j2],v_t1);
               _mm_storeh_pi((__m64 *)&f[i+j2+nxhd],v_t1);
/*             f[i+j1] += t2; */
               v_t1 = _mm_add_ps(v_t3,v_t2);
               _mm_storel_pi((__m64 *)&f[i+j1],v_t1);
               _mm_storeh_pi((__m64 *)&f[i+j1+nxhd],v_t1);
            }
/* loop over remaining elements */
            for (j = nss; j < ns; j++) {
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
void csse2fft2rm2x(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nyi,
                   int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the x part of 2 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   y, using complex arithmetic, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, two inverse fourier transforms in x are performed
   f[m][n][0:1] = (1/nx*ny)*sum(f[k][j][0:1]*exp(-sqrt(-1)*2pi*n*j/nx))
   if isign = 1, two forward fourier transforms in x are performed
   f[k][j][0:1] = sum(f[m][n][0:1]*exp(sqrt(-1)*2pi*n*j/nx))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = second dimension of f >= nx/2
   nyd = third dimension of f >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j][0:1] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][1][0:1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0][0:1]) = real part of mode nx/2,0 and
   imag(f[0][ny/2][0:1]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
   requires SSE2, f needs to be 16 byte aligned
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, jj, j1, k1, k2, ns, ns2, km, kmr, joff;
   int nrxb;
   float ani;
/* float complex t1, t2, t3; */
   __m128 v_m, v_n, v_t1, v_t2, v_t3, v_t4, v_ani;
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
   v_m = _mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);
   v_n = _mm_set_ps(-1.0f,1.0f,-1.0f,1.0f);
   v_t1 = _mm_setzero_ps();
   v_t3 = _mm_setzero_ps();
   if (isign > 0)
      goto L100;
/* inverse fourier transform */
   nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,joff,ani,v_t1,v_t2,v_t3,v_t4, \
v_ani)
   for (i = nyi-1; i < nyt; i++) {
      joff = 2*nxhd*i;
/* swap complex components */
      for (j = 0; j < nxh; j++) {
         v_t1 = _mm_load_ps((float *)&f[2*j+joff]);
         v_t1 = _mm_shuffle_ps(v_t1,v_t1,216);
         _mm_store_ps((float *)&f[2*j+joff],v_t1);
       }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
/*          t1 = f[2*j1+joff];   */
/*          t2 = f[1+2*j1+joff]; */
            v_t1 = _mm_load_ps((float *)&f[2*j1+joff]);
/*          f[2*j1+joff] = f[2*j+joff];     */
/*          f[1+2*j1+joff] = f[1+2*j+joff]; */
            v_t2 = _mm_load_ps((float *)&f[2*j+joff]);
            _mm_store_ps((float *)&f[2*j1+joff],v_t2);
/*          f[2*j+joff] = t1;   */
/*          f[1+2*j+joff] = t2; */
            _mm_store_ps((float *)&f[2*j+joff],v_t1);
         }
      }
/* then transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = 2*ns2*k;
            k2 = k1 + 2*ns;
            for (j = 0; j < ns; j++) {
/*             t1 = sct[kmr*j]; */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_movelh_ps(v_t1,v_t1);
/*             t2 = t1*f[2*j+k2+joff];   */
/*             t3 = t1*f[1+2*j+k2+joff]; */
               v_t2 = _mm_load_ps((float *)&f[2*j+k2+joff]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             f[2*j+k2+joff] = f[2*j+k1+joff] - t2;     */
/*             f[1+2*j+k2+joff] = f[1+2*j+k1+joff] - t3; */
               v_t3 = _mm_load_ps((float *)&f[2*j+k1+joff]);
               v_t4 = _mm_sub_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[2*j+k2+joff],v_t4);
/*             f[2*j+k1+joff] += t2;   */
/*             f[1+2*j+k1+joff] += t3; */
               v_t2 = _mm_add_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[2*j+k1+joff],v_t2);
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxy/nx;
      ani = 0.5/(((float) nx)*((float) ny));
      v_ani = _mm_set1_ps(ani);
      for (j = 1; j < nxhh; j++) {
/*       t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I; */
         v_t3 = _mm_loadl_pi(v_t3,(__m64 *)&sct[kmr*j]);
         v_t3 = _mm_movelh_ps(v_t3,v_t3);
         v_t3 = _mm_mul_ps(v_t3,v_m);
         v_t3 = _mm_shuffle_ps(v_t3,v_t3,177);
/*       t2 = conjf(f[jj+2*(nxh-j)+joff]); */
         v_t2 = _mm_load_ps((float *)&f[2*(nxh-j)+joff]);
         v_t2 = _mm_mul_ps(v_t2,v_n);
/*       t1 = f[jj+2*j+joff] + t2; */
         v_t4 = _mm_load_ps((float *)&f[2*j+joff]);
         v_t1 = _mm_add_ps(v_t4,v_t2);
/*       t2 = (f[jj+2*j+joff] - t2)*t3; */
         v_t2 = _mm_sub_ps(v_t4,v_t2);
         v_t4 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,160));
         v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
         v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,245));
         v_t2 = _mm_add_ps(v_t4,_mm_mul_ps(v_t2,v_m));
/*       f[jj+2*j+joff] = ani*(t1 + t2); */
         v_t3 = _mm_mul_ps(v_ani,_mm_add_ps(v_t1,v_t2));
         _mm_store_ps((float *)&f[2*j+joff],v_t3);
/*       f[jj+2*(nxh-j)+joff] = ani*conjf(t1 - t2); */
         v_t4 = _mm_mul_ps(v_ani,_mm_mul_ps(_mm_sub_ps(v_t1,v_t2),v_n));
         _mm_store_ps((float *)&f[2*(nxh-j)+joff],v_t4);
      }
      ani = 2.0*ani;
      for (jj = 0; jj < 2; jj++) {
         f[jj+2*nxhh+joff] = ani*conjf(f[jj+2*nxhh+joff]);
         f[jj+joff] = ani*((crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
L100: nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,joff,v_t1,v_t2,v_t3,v_t4)
   for (i = nyi-1; i < nyt; i++) {
      joff = 2*nxhd*i;
/* scramble coefficients */
      kmr = nxy/nx;
      for (j = 1; j < nxhh; j++) {
/*       t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I; */
         v_t3 = _mm_loadl_pi(v_t3,(__m64 *)&sct[kmr*j]);
         v_t3 = _mm_movelh_ps(v_t3,v_t3);
         v_t3 = _mm_shuffle_ps(v_t3,v_t3,177);
/*       t2 = conjf(f[jj+2*(nxh-j)+joff]); */
         v_t2 = _mm_load_ps((float *)&f[2*(nxh-j)+joff]);
         v_t2 = _mm_mul_ps(v_t2,v_n);
/*       t1 = f[jj+2*j+joff] + t2; */
         v_t4 = _mm_load_ps((float *)&f[2*j+joff]);
         v_t1 = _mm_add_ps(v_t4,v_t2);
/*       t2 = (f[jj+2*j+joff] - t2)*t3; */
         v_t2 = _mm_sub_ps(v_t4,v_t2);
         v_t4 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,160));
         v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
         v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,245));
         v_t2 = _mm_add_ps(v_t4,_mm_mul_ps(v_t2,v_m));
/*       f[jj+2*j+joff] = t1 + t2; */
         v_t3 = _mm_add_ps(v_t1,v_t2);
          _mm_store_ps((float *)&f[2*j+joff],v_t3);
/*       f[jj+2*(nxh-j)+joff] = conjf(t1 - t2); */
         v_t4 = _mm_mul_ps(_mm_sub_ps(v_t1,v_t2),v_n);
         _mm_store_ps((float *)&f[2*(nxh-j)+joff],v_t4);
      }
      for (jj = 0; jj < 2; jj++) {
         f[jj+2*nxhh+joff] = 2.0*conjf(f[jj+2*nxhh+joff]);
         f[jj+joff] = (crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I;
      }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
/*          t1 = f[2*j1+joff];   */
/*          t2 = f[1+2*j1+joff]; */
            v_t1 = _mm_load_ps((float *)&f[2*j1+joff]);
/*          f[2*j1+joff] = f[2*j+joff];     */
/*          f[1+2*j1+joff] = f[1+2*j+joff]; */
            v_t2 = _mm_load_ps((float *)&f[2*j+joff]);
            _mm_store_ps((float *)&f[2*j1+joff],v_t2);
/*          f[2*j+joff] = t1;   */
/*          f[1+2*j+joff] = t2; */
            _mm_store_ps((float *)&f[2*j+joff],v_t1);
         }
      }
/* then transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = 2*ns2*k;
            k2 = k1 + 2*ns;
            for (j = 0; j < ns; j++) {
/*             t1 = conjf(sct[kmr*j]); */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_movelh_ps(v_t1,v_t1);
               v_t1 = _mm_mul_ps(v_t1,v_n);
/*             t2 = t1*f[2*j+k2+joff];   */
/*             t3 = t1*f[1+2*j+k2+joff]; */
               v_t2 = _mm_load_ps((float *)&f[2*j+k2+joff]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             f[2*j+k2+joff] = f[2*j+k1+joff] - t2;    */
/*             f[1+2*j+k2+joff] = f[1+2*j+k1+joff] - t3; */
               v_t3 = _mm_load_ps((float *)&f[2*j+k1+joff]);
               v_t4 = _mm_sub_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[2*j+k2+joff],v_t4);
/*             f[2*j+k1+joff] += t2;    */
/*             f[1+2*j+k1+joff] +=  t3; */
               v_t2 = _mm_add_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[2*j+k1+joff],v_t2);
            }
         }
         ns = ns2;
      }
/* swap complex components */
      for (j = 0; j < nxh; j++) {
         v_t1 = _mm_load_ps((float *)&f[2*j+joff]);
         v_t1 = _mm_shuffle_ps(v_t1,v_t1,216);
         _mm_store_ps((float *)&f[2*j+joff],v_t1);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2fft2rm2y(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxi,
                   int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the y part of 2 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   x, using complex arithmetic, with OpenMP
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny
   indx/indy = exponent which determines length in x/y direction,
   where nx=2**indx, ny=2**indy
   if isign = -1, two inverse fourier transforms in y are performed
   f[m][n][0:1] = *sum(f[k][j][0:1]*exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, two forward fourier transforms in y are performed
   f[k][j][0:1] = sum(f[m][n][0:1]*exp(sqrt(-1)*2pi*n*j/nx))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nxi = initial x index used
   nxp = number of x indices used
   nxhd = second dimension of f >= nx/2
   nyd = third dimension of f >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j][0:1] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][1][0:1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0][0:1]) = real part of mode nx/2,0 and
   imag(f[0][ny/2][0:1]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
   requires SSE2, f needs to be 16 byte aligned
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, koff;
   int nryb;
   float complex t1;
   __m128 v_m, v_n, v_t1, v_t2, v_t3, v_t4;
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
   v_m = _mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);
   v_n = _mm_set_ps(-1.0f,1.0f,-1.0f,1.0f);
   v_t1 = _mm_setzero_ps();
   if (isign > 0)
      goto L80;
/* inverse fourier transform */
   nryb = nxhy/ny;
   nry = nxy/ny;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,koff,v_t1,v_t2,v_t3,v_t4)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = 2*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = 2*nxhd*k1;
/*          t1 = f[2*i+k1];   */
/*          t2 = f[1+2*i+k1]; */
            v_t1 = _mm_load_ps((float *)&f[2*i+k1]);
/*          f[2*i+k1] = f[2*i+koff];     */
/*          f[1+2*i+k1] = f[1+2*i+koff]; */
            v_t2 = _mm_load_ps((float *)&f[2*i+koff]);
            _mm_store_ps((float *)&f[2*i+k1],v_t2);
/*          f[2*i+koff] = t1;   */
/*          f[1+2*i+koff] = t2; */
            _mm_store_ps((float *)&f[2*i+koff],v_t1);
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
               j1 = 2*nxhd*(j + k1);
               j2 = 2*nxhd*(j + k2);
/*             t1 = sct[kmr*j]; */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_movelh_ps(v_t1,v_t1);
/*             t2 = t1*f[2*i+j2];   */
/*             t3 = t1*f[1+2*i+j2]; */
               v_t2 = _mm_load_ps((float *)&f[2*i+j2]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             f[2*i+j2] = f[2*i+j1] - t2;     */
/*             f[1+2*i+j2] = f[1+2*i+j1] - t3; */
               v_t3 = _mm_load_ps((float *)&f[2*i+j1]);
               v_t4 = _mm_sub_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[2*i+j2],v_t4);
/*             f[2*i+j1] += t2;   */
/*             f[1+2*i+j1] += t3; */
               v_t2 = _mm_add_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[2*i+j1],v_t2);
            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if (nxi==1) {
      for (k = 1; k < nyh; k++) {
         koff = 2*nxhd*k;
         k1 = 2*nxhd*ny - koff;
         for (jj = 0; jj < 2; jj++) {
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
         koff = 2*nxhd*k;
         k1 = 2*nxhd*ny - koff;
         for (jj = 0; jj < 2; jj++) {
            t1 = cimagf(f[jj+k1]) + crealf(f[jj+k1])*_Complex_I;
            f[jj+k1] = conjf(f[jj+koff] - t1);
            f[jj+koff] += t1;
         }
      }
   }
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,j2,koff,v_t1,v_t2,v_t3,v_t4)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = 2*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = 2*nxhd*k1;
/*          t1 = f[2*i+k1];   */
/*          t2 = f[1+2*i+k1]; */
            v_t1 = _mm_load_ps((float *)&f[2*i+k1]);
/*          f[2*i+k1] = f[2*i+koff];     */
/*          f[1+2*i+k1] = f[1+2*i+koff]; */
            v_t2 = _mm_load_ps((float *)&f[2*i+koff]);
            _mm_store_ps((float *)&f[2*i+k1],v_t2);
/*          f[2*i+koff] = t1;   */
/*          f[1+2*i+koff] = t2; */
            _mm_store_ps((float *)&f[2*i+koff],v_t1);
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
               j1 = 2*nxhd*(j + k1);
               j2 = 2*nxhd*(j + k2);
/*             t1 = conjf(sct[kmr*j]); */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_movelh_ps(v_t1,v_t1);
               v_t1 = _mm_mul_ps(v_t1,v_n);
/*             t2 = t1*f[2*i+j2];   */
/*             t3 = t1*f[1+2*i+j2]; */
               v_t2 = _mm_load_ps((float *)&f[2*i+j2]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             f[2*i+j2] = f[2*i+j1] - t2;     */
/*             f[1+2*i+j2] = f[1+2*i+j1] - t3; */
               v_t3 = _mm_load_ps((float *)&f[2*i+j1]);
               v_t4 = _mm_sub_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[2*i+j2],v_t4);
/*             f[2*i+j1] += t2;   */
/*             f[1+2*i+j1] += t3; */
               v_t2 = _mm_add_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[2*i+j1],v_t2);
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2wfft2rmx(float complex f[], int isign, int mixup[],
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
      csse2fft2rmxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                    nxyhd);
/* perform y fft */
      csse2fft2rmxy(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                    nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      csse2fft2rmxy(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                    nxyhd);
/* perform x fft */
      csse2fft2rmxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                    nxyhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2wfft2rm2(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxhd,
                   int nyd, int nxhyd, int nxyhd) {
/* wrapper function for 2 2d real to complex ffts, with packed data */
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
      csse2fft2rm2x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                    nxyhd);
/* perform y fft */
      csse2fft2rm2y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                    nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      csse2fft2rm2y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                    nxyhd);
/* perform x fft */
      csse2fft2rm2x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                    nxyhd);
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void csse2gppush2lt_(float *ppart, float *fxy, int *kpic, float *qbm,
                     float *dt, float *ek, int *idimp, int *nppmx,
                     int *nx, int *ny, int *mx, int *my, int *nxv,
                     int *nyv, int *mx1, int *mxy1, int *ipbc) {
   csse2gppush2lt(ppart,fxy,kpic,*qbm,*dt,ek,*idimp,*nppmx,*nx,*ny,*mx,
                  *my,*nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2gppushf2lt_(float *ppart, float *fxy, int *kpic, int *ncl,
                      int *ihole, float *qbm, float *dt, float *ek, 
                      int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                      int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                      int *ntmax, int *irc) {
   csse2gppushf2lt(ppart,fxy,kpic,ncl,ihole,*qbm,*dt,ek,*idimp,*nppmx,
                   *nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2gppost2lt_(float *ppart, float *q, int *kpic, float *qm,
                     int *nppmx, int *idimp, int *mx, int *my, int *nxv,
                     int *nyv, int *mx1, int *mxy1) {
   csse2gppost2lt(ppart,q,kpic,*qm,*nppmx,*idimp,*mx,*my,*nxv,*nyv,*mx1,
                  *mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void csse2pporder2lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                      int *ihole, int *idimp, int *nppmx, int *nx,
                      int *ny, int *mx, int *my, int *mx1, int *my1,
                      int *npbmx, int *ntmax, int *irc) {
   csse2pporder2lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*nx,*ny,
                   *mx,*my,*mx1,*my1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2pporderf2lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                       int *ihole, int *idimp, int *nppmx, int *mx1,
                       int *my1, int *npbmx, int *ntmax, int *irc) {
   csse2pporderf2lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*mx1,*my1,
                    *npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2cguard2l_(float *fxy, int *nx, int *ny, int *nxe, int *nye) {
   csse2cguard2l(fxy,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void csse2aguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye) {
   csse2aguard2l(q,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void csse2mpois22_(float complex *q, float complex *fxy, int *isign,
                   float complex *ffc, float *ax, float *ay,
                   float *affp, float *we, int *nx, int *ny, int *nxvh,
                   int *nyv, int *nxhd, int *nyhd) {
   csse2mpois22(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*nxvh,*nyv,
                *nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csse2wfft2rmx_(float complex *f, int *isign, int *mixup,
                    float complex *sct, int *indx, int *indy, int *nxhd,
                    int *nyd, int *nxhyd, int *nxyhd) {
   csse2wfft2rmx(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,
                 *nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csse2wfft2rm2_(float complex *f, int *isign, int *mixup,
                    float complex *sct, int *indx, int *indy, int *nxhd,
                    int *nyd, int *nxhyd, int *nxyhd) {
   csse2wfft2rm2(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,
                 *nxyhd);
   return;
}

