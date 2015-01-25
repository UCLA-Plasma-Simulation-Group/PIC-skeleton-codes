/* SSE2 C Library for Skeleton 2D Electrostatic Vector PIC Code */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <xmmintrin.h>
#include "ssepush2.h"

/*--------------------------------------------------------------------*/
void csse2xiscan2(int *isdata, int nths) {
/* performs local prefix reduction of integer data shared by threads */
/* using binary tree method, exclusive scan. */
/* requires SSE2, isdata needs to be 16 byte aligned */
/* local data */
   int j, ns, isum, ist;
   __m128i v_m1, v_m2, v_it0, v_it, v_is, v_ioff;
   __attribute__((aligned(16))) int ll[4];
   ns = 4*(nths/4);
   v_m1 = _mm_set_epi32(0,-1,0,-1);
   v_m2 = _mm_set_epi32(0,-1,-1,0);
   isum = 0;
   v_ioff = _mm_set1_epi32(isum);
/* vector loop over elements in blocks of 4 */
   for (j = 0; j < ns; j+=4) {
/* load data */
      v_it0 = _mm_load_si128((__m128i *)&isdata[j]);
/* first pass */
      v_is = _mm_slli_si128(_mm_and_si128(v_it0,v_m1),4);
      v_it = _mm_add_epi32(v_is,v_it0);
/* second pass */
      v_is = _mm_shuffle_epi32(v_it,212);
      v_is = _mm_slli_si128(_mm_and_si128(v_is,v_m2),4);
      v_it = _mm_add_epi32(v_is,v_it);
/* add offset */
      v_it = _mm_add_epi32(v_it,v_ioff);
/* next offset */
      v_ioff = _mm_shuffle_epi32(v_it,255);
      v_it = _mm_sub_epi32(v_it,v_it0);
/* write data */
      _mm_store_si128((__m128i *)&isdata[j],v_it);
   }
   if (ns > 0) {
      _mm_store_si128((__m128i *)ll,v_ioff);
      isum = ll[0];
   }
/* loop over remaining elements */
   for (j = ns; j < nths; j++) {
      ist = isdata[j];
      isdata[j] = isum;
      isum += ist;
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2gpush2lt(float part[], float fxy[], float qbm, float dt,
                   float *ek, int idimp, int nop, int npe, int nx,
                   int ny, int nxv, int nyv, int ipbc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions.
   vector version using guard cells
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
   requires SSE2, part needs to be 16 byte aligned
   npe needs to be a multiple of 4
local data                                                            */
   int j, nps, nn, mm;
   float qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
   double sum1;
   __m128i v_nxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qtm, v_dt, v_one;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_vx, v_vy;
   __m128 v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 a, b, c, d;
   __m128d v_sum1, v_d;
   __attribute__((aligned(16))) unsigned int ll[4];
   __attribute__((aligned(16))) double dd[2];
   qtm = qbm*dt;
   sum1 = 0.0;
   nps = 4*(nop/4);
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
   v_nxv = _mm_set1_epi32(nxv);
   v_qtm = _mm_set1_ps(qtm);
   v_one = _mm_set1_ps(1.0f);
   v_dt = _mm_set1_ps(dt);
   v_edgelx = _mm_set1_ps(edgelx);
   v_edgely = _mm_set1_ps(edgely);
   v_edgerx = _mm_set1_ps(edgerx);
   v_edgery = _mm_set1_ps(edgery);
   v_sum1 = _mm_set1_pd(0.0);
/* vector loop over particles in blocks of 4 */
   for (j = 0; j < nps; j+=4) {
/* find interpolation weights */
/*    x = part[j];     */
/*    y = part[j+npe]; */
      v_x = _mm_load_ps(&part[j]);
      v_y = _mm_load_ps(&part[j+npe]);
/*    nn = x; */
/*    mm = y; */
      v_nn = _mm_cvttps_epi32(v_x);
      v_mm = _mm_cvttps_epi32(v_y);
/*    dxp = x - (float) nn; */
      v_dxp = _mm_sub_ps(v_x,_mm_cvtepi32_ps(v_nn));
/*    dyp = y - (float) mm; */
      v_dyp = _mm_sub_ps(v_y,_mm_cvtepi32_ps(v_mm));
/*    nn = 2*(nn + nxv*mm); */
      v_it = _mm_mul_epu32(v_nxv,_mm_srli_si128(v_mm,4));
      v_mm = _mm_mul_epu32(v_mm,v_nxv);
      v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
      v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),1);
/*    amx = 1.0f - dxp; */
/*    amy = 1.0f - dyp; */
      v_amx = _mm_sub_ps(v_one,v_dxp);
      v_amy = _mm_sub_ps(v_one,v_dyp);
/* find acceleration */
/* load fields, for lower left/right */
      _mm_store_si128((__m128i *)ll,v_nn);
      a = _mm_loadu_ps(&fxy[ll[0]]);
      b = _mm_loadu_ps(&fxy[ll[1]]);
      c = _mm_loadu_ps(&fxy[ll[2]]);
      d = _mm_loadu_ps(&fxy[ll[3]]);
/* transpose so a,b,c,d contain first 4 fields for each of 4 particles */
      _MM_TRANSPOSE4_PS(a,b,c,d);
/*    dx = amx*fxy[nn];   */
/*    dy = amx*fxy[nn+1]; */
      v_dx = _mm_mul_ps(v_amx,a);
      v_dy = _mm_mul_ps(v_amx,b);
/*    dx = amy*(dxp*fxy[nn+2] + dx); */
/*    dy = amy*(dxp*fxy[nn+3] + dy); */
      v_dx = _mm_mul_ps(v_amy,_mm_add_ps(_mm_mul_ps(v_dxp,c),v_dx));
      v_dy = _mm_mul_ps(v_amy,_mm_add_ps(_mm_mul_ps(v_dxp,d),v_dy));
/*    nn += 2*nxv; */
/* load fields, for upper left/right */
      a = _mm_loadu_ps(&fxy[ll[0]+2*nxv]);
      b = _mm_loadu_ps(&fxy[ll[1]+2*nxv]);
      c = _mm_loadu_ps(&fxy[ll[2]+2*nxv]);
      d = _mm_loadu_ps(&fxy[ll[3]+2*nxv]);
/* transpose so a,b,c,d contain next 4 fields for each of 4 particles */
      _MM_TRANSPOSE4_PS(a,b,c,d);
/*    vx = amx*fxy[nn];   */
/*    vy = amx*fxy[nn+1]; */
      a = _mm_mul_ps(v_amx,a);
      b = _mm_mul_ps(v_amx,b);
/*    dx += dyp*(dxp*fxy[nn+2] + vx); */
/*    dy += dyp*(dxp*fxy[nn+3] + vy); */
      a = _mm_mul_ps(v_dyp,_mm_add_ps(_mm_mul_ps(v_dxp,c),a));
      b = _mm_mul_ps(v_dyp,_mm_add_ps(_mm_mul_ps(v_dxp,d),b));
      v_dx = _mm_add_ps(v_dx,a);
      v_dy = _mm_add_ps(v_dy,b);
/* new velocity */
/*    dxp = part[j+2*npe]; */
/*    dyp = part[j+3*npe]; */
      v_dxp = _mm_load_ps(&part[j+2*npe]);
      v_dyp = _mm_load_ps(&part[j+3*npe]);
/*    vx = dxp + qtm*dx; */
/*    vy = dyp + qtm*dy; */
      v_vx = _mm_add_ps(v_dxp,_mm_mul_ps(v_qtm,v_dx));
      v_vy = _mm_add_ps(v_dyp,_mm_mul_ps(v_qtm,v_dy));
/* average kinetic energy */
/*    dxp += vx; */
/*    dyp += vy; */
      v_dxp = _mm_add_ps(v_dxp,v_vx);
      v_dyp = _mm_add_ps(v_dyp,v_vy);
/*    sum1 += dxp*dxp + dyp*dyp; */
      v_at = _mm_mul_ps(v_dxp,v_dxp);
      v_at = _mm_add_ps(v_at,_mm_mul_ps(v_dyp,v_dyp));
/* convert to double precision before accumulating */
      v_d = _mm_cvtps_pd(v_at);
      v_sum1 = _mm_add_pd(v_sum1,v_d);
      v_it = _mm_srli_si128((__m128i)v_at,8);
      v_d = _mm_cvtps_pd((__m128)v_it);
      v_sum1 = _mm_add_pd(v_sum1,v_d);
/* new position */
/*    dx = x + vx*dt; */
/*    dy = y + vy*dt; */
      v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_dt));
      v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_dt));
/* periodic boundary conditions */
      if (ipbc==1) {
/*       if (dx < edgelx) dx += edgerx; */
         v_dxp = _mm_and_ps(v_edgerx,_mm_cmplt_ps(v_dx,v_edgelx));
         v_dx = _mm_add_ps(v_dx,v_dxp);
/*       if (dx >= edgerx) dx -= edgerx; */
         v_dyp = _mm_and_ps(v_edgerx,_mm_cmpge_ps(v_dx,v_edgerx));
         v_dx = _mm_sub_ps(v_dx,v_dyp);
/*       if (dy < edgely) dy += edgery; */
         v_dxp = _mm_and_ps(v_edgery,_mm_cmplt_ps(v_dy,v_edgely));
         v_dy = _mm_add_ps(v_dy,v_dxp);
/*       if (dy >= edgery) dy -= edgery; */
         v_dyp = _mm_and_ps(v_edgery,_mm_cmpge_ps(v_dy,v_edgery));
         v_dy = _mm_sub_ps(v_dy,v_dyp);
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          vx = -vx;                           */
/*       }                                      */
         v_at = _mm_cmplt_ps(v_dx,v_edgelx);
         v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dx,v_edgerx));
         v_x = _mm_and_ps(v_at,v_x);
         v_dx = _mm_add_ps(_mm_andnot_ps(v_at,v_dx),v_x);
         v_dxp = _mm_and_ps(v_at,v_vx);
         v_vx = _mm_sub_ps(_mm_andnot_ps(v_at,v_vx),v_dxp);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          vy = -vy;                           */
/*       }                                      */
         v_at = _mm_cmplt_ps(v_dy,v_edgely);
         v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dy,v_edgery));
         v_y = _mm_and_ps(v_at,v_y);
         v_dy = _mm_add_ps(_mm_andnot_ps(v_at,v_dy),v_y);
         v_dyp = _mm_and_ps(v_at,v_vy);
         v_vy = _mm_sub_ps(_mm_andnot_ps(v_at,v_vy),v_dyp);
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          vx = -vx;                           */
/*       }                                      */
         v_at = _mm_cmplt_ps(v_dx,v_edgelx);
         v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dx,v_edgerx));
         v_x = _mm_and_ps(v_at,v_x);
         v_dx = _mm_add_ps(_mm_andnot_ps(v_at,v_dx),v_x);
         v_dxp = _mm_and_ps(v_at,v_vx);
         v_vx = _mm_sub_ps(_mm_andnot_ps(v_at,v_vx),v_dxp);
/*       if (dy < edgely) dy += edgery; */
         v_dxp = _mm_and_ps(v_edgery,_mm_cmplt_ps(v_dy,v_edgely));
         v_dy = _mm_add_ps(v_dy,v_dxp);
/*       if (dy >= edgery) dy -= edgery; */
         v_dyp = _mm_and_ps(v_edgery,_mm_cmpge_ps(v_dy,v_edgery));
         v_dy = _mm_sub_ps(v_dy,v_dyp);
      }
/* set new position */
/*    part[j] = dx;     */
/*    part[j+npe] = dy; */
      _mm_store_ps(&part[j],v_dx);
      _mm_store_ps(&part[j+npe],v_dy);
/* set new velocity */
/*    part[j+2*npe] = vx; */
/*    part[j+3*npe] = vy; */
      _mm_store_ps(&part[j+2*npe],v_vx);
      _mm_store_ps(&part[j+3*npe],v_vy);
   }
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
/* *ek += 0.125f*sum1; */
   _mm_store_pd(&dd[0],v_sum1);
   for (j = 1; j < 2; j++) {
      dd[0] += dd[j];
   }
   *ek += 0.125f*(sum1 + dd[0]);
   return;
}

/*--------------------------------------------------------------------*/
void csse2gpost2lt(float part[], float q[], float qm, int nop, int npe,
                   int idimp, int nxv, int nyv) {
/* for 2d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   vector version using guard cells
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
   requires SSE2, part needs to be 16 byte aligned
   npe needs to be a multiple of 4
local data                                                            */
   int j, nps, nn, mm;
   float x, y, dxp, dyp, amx, amy;
   __m128i v_nxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qm, v_one;
   __m128 v_x, v_y, v_dxp, v_dyp, v_amx, v_amy;
   __m128 a, b, c, d;
   __attribute__((aligned(16))) unsigned int ll[4];
   nps = 4*(nop/4);
   v_nxv = _mm_set1_epi32(nxv);
   v_qm = _mm_set1_ps(qm);
   v_one = _mm_set1_ps(1.0f);
/* vector loop over particles in blocks of 4 */
   for (j = 0; j < nps; j+=4) {
/* find interpolation weights */
/*    x = part[j];     */
/*    y = part[j+npe]; */
      v_x = _mm_load_ps(&part[j]);
      v_y = _mm_load_ps(&part[j+npe]);
/*    nn = x; */
/*    mm = y; */
      v_nn = _mm_cvttps_epi32(v_x);
      v_mm = _mm_cvttps_epi32(v_y);
/*    dxp = qm*(x - (float) nn); */
      v_dxp = _mm_sub_ps(v_x,_mm_cvtepi32_ps(v_nn));
      v_dxp = _mm_mul_ps(v_dxp,v_qm);
/*    dyp = y - (float) mm; */
      v_dyp = _mm_sub_ps(v_y,_mm_cvtepi32_ps(v_mm));
/*    nn = nn + nxv*mm; */
      v_it = _mm_mul_epu32(v_nxv,_mm_srli_si128(v_mm,4));
      v_mm = _mm_mul_epu32(v_mm,v_nxv);
      v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
      v_nn = _mm_add_epi32(v_nn,v_mm);
/*    amx = qm - dxp;   */
/*    amy = 1.0f - dyp; */
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
/* deposit charge */
/*    x = q[nn] + amx*amy;   */
/*    y = q[nn+1] + dxp*amy; */
/*    q[nn] = x;             */
/*    q[nn+1] = y;           */
/*    nn += nxv;             */
/*    x = q[nn] + amx*dyp;   */
/*    y = q[nn+1] + dxp*dyp; */
/*    q[nn] = x;             */
/*    q[nn+1] = y;           */
/* deposit for first particle */
      mm = ll[0];
      v_x = _mm_loadl_pi(v_x,(__m64 *)&q[mm]);
      v_x = _mm_loadh_pi(v_x,(__m64 *)&q[mm+nxv]);
      v_x = _mm_add_ps(v_x,a);
      _mm_storel_pi((__m64 *)&q[mm],v_x);
      _mm_storeh_pi((__m64 *)&q[mm+nxv],v_x);
/* deposit for second particle */
      mm = ll[1];
      v_y = _mm_loadl_pi(v_y,(__m64 *)&q[mm]);
      v_y = _mm_loadh_pi(v_y,(__m64 *)&q[mm+nxv]);
      v_y = _mm_add_ps(v_y,b);
      _mm_storel_pi((__m64 *)&q[mm],v_y);
      _mm_storeh_pi((__m64 *)&q[mm+nxv],v_y);
/* deposit for third particle */
      mm = ll[2];
      v_x = _mm_loadl_pi(v_x,(__m64 *)&q[mm]);
      v_x = _mm_loadh_pi(v_x,(__m64 *)&q[mm+nxv]);
      v_x = _mm_add_ps(v_x,c);
      _mm_storel_pi((__m64 *)&q[mm],v_x);
      _mm_storeh_pi((__m64 *)&q[mm+nxv],v_x);
/* deposit for fourth particle */
      mm = ll[3];
      v_y = _mm_loadl_pi(v_y,(__m64 *)&q[mm]);
      v_y = _mm_loadh_pi(v_y,(__m64 *)&q[mm+nxv]);
      v_y = _mm_add_ps(v_y,d);
      _mm_storel_pi((__m64 *)&q[mm],v_y);
      _mm_storeh_pi((__m64 *)&q[mm+nxv],v_y);
   }
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
}

/*--------------------------------------------------------------------*/
void csse2dsortp2ylt(float parta[], float partb[], int npic[],
                     int idimp, int nop, int npe, int ny1) {
/* this subroutine sorts particles by y grid
   linear interpolation
   parta/partb = input/output particle arrays
   parta[1][n] = position y of particle n
   npic = address offset for reordering particles
   idimp = size of phase space = 4
   nop = number of particles
   npe = first dimension of particle array
   ny1 = system length in y direction + 1
   requires SSE2, part needs to be 16 byte aligned
   npe needs to be a multiple of 4
local data                                                            */
   int i, j, k, m, nps, ip;
/* int isum, ist; */
   __m128i v_m, v_it, v_one;
   __m128 v_at;
   __attribute__((aligned(16))) unsigned int ll[4], pp[4];
   nps = 4*(nop/4);
   v_one = _mm_set1_epi32(1);
/* clear counter array */
/* for (k = 0; k < ny1; k++) { */
/*    npic[k] = 0;             */
/* }                           */
   memset((void *)npic,0,ny1*sizeof(int));
/* find how many particles in each grid */
/* vector loop over particles in blocks of 4 */
   for (j = 0; j < nps; j+=4) {
/*    m = parta[j+npe]; */
      v_m = _mm_cvttps_epi32(_mm_load_ps(&parta[j+npe]));
      _mm_store_si128((__m128i *)ll,v_m);
/*    npic[m] += 1; */
      for (k = 0; k < 4; k++) {
         m = ll[k];
         v_it = (__m128i)_mm_load_ss((float *)&npic[m]);
         v_it = _mm_add_epi32(v_it,v_one);
         _mm_store_ss((float *)&npic[m],(__m128)v_it);
      }
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
      m = parta[j+npe];
      npic[m] += 1;
   }
/* find address offset */
/* isum = 0;                   */
/* for (k = 0; k < ny1; k++) { */
/*    ist = npic[k];           */
/*    npic[k] = isum;          */
/*    isum += ist;             */
/* }                           */
   csse2xiscan2(npic,ny1);
/* find addresses of particles at each grid and reorder particles */
/* vector loop over particles in blocks of 4 */
   for (j = 0; j < nps; j+=4) {
/*    m = parta[j+npe]; */
      v_m = _mm_cvttps_epi32(_mm_load_ps(&parta[j+npe]));
      _mm_store_si128((__m128i *)ll,v_m);
/*    ip = npic[m]; */
/*    npic[m] = ip + 1; */
      for (k = 0; k < 4; k++) {
         m = ll[k];
         ip = npic[m];
         npic[m] = ip + 1;
         pp[k] = ip;
      }
      for (i = 0; i < idimp; i++) {
/*       partb[ip+npe*i] = parta[j+npe*i]; */
         v_at = _mm_load_ps(&parta[j+npe*i]);
         for (k = 0; k < 4; k++) {
            _mm_store_ss(&partb[pp[k]+npe*i],v_at);
            v_at = (__m128)_mm_srli_si128((__m128i)v_at,4);
         }
      }
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
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
void csse2pois22(float complex q[], float complex fxy[], int isign,
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
   double wp;
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
L30: wp = 0.0;
   v_wp = _mm_set1_pd(0.0);
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      v_dky = _mm_mul_ps(v_dny,_mm_cvtepi32_ps(_mm_set1_epi32(k)));
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
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
/*       fxy[2*j+2*kj] = at2*zt1; */
/*       fxy[1+2*j+2*kj] = at3*zt1; */
         v_at4 = _mm_mul_ps(v_at2,v_zt1);
         v_zt4 = _mm_mul_ps(v_at3,v_zt1);
/* reorder write */
         v_zt3 = _mm_shuffle_ps(v_at4,v_zt4,68);
         v_zt4 = _mm_shuffle_ps(v_at4,v_zt4,238);
         _mm_store_ps((float *)&fxy[2*(j+kj)],v_zt3);
         _mm_store_ps((float *)&fxy[2*(j+1+kj)],v_zt4);
/*       fxy[2*j+2*k1] = at2*zt2; */
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
   }
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
/*    fxy[1+2*j] = zero; */
      v_at4 = _mm_mul_ps(v_at2,v_zt1);
/* reorder write */
      v_zt3 = _mm_shuffle_ps(v_at4,v_zero,68);
      v_zt4 = _mm_shuffle_ps(v_at4,v_zero,238);
      _mm_store_ps((float *)&fxy[2*j],v_zt3);
      _mm_store_ps((float *)&fxy[2*j+2],v_zt4);
/*    fxy[2*j+k1] = zero; */
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
/* *we = wp*(float) (nx*ny); */
   _mm_store_pd(&dd[0],v_wp);
   for (j = 1; j < 2; j++) {
      dd[0] += dd[j];
   }
   *we = (wp + dd[0])*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void csse2fft2rxx(float complex f[], int isign, int mixup[],
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
   requires SSE2, f needs to be 16 byte aligned
   nxhd need to be a multiple of 2
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, j1, k1, k2, ns, ns2, km, kmr, joff;
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
/*       t1 = f[j1+joff]; */
         v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&f[j1+joff]);
/*       f[j1+joff] = f[j+joff]; */
         v_t2 = _mm_loadl_pi(v_t2,(__m64 *)&f[j+joff]);
         _mm_storel_pi((__m64 *)&f[j1+joff],v_t2);
/*       f[j+joff] = t1; */
         _mm_storel_pi((__m64 *)&f[j+joff],v_t1);
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
         nss = 2*(ns/2);
         for (i = nyi-1; i < nyt; i++) {
            joff = nxhd*i;
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
      }
      ns = ns2;
   }
/* unscramble coefficients and normalize */
   kmr = nxy/nx;
   ani = 1.0/(float) (2*nx*ny);
   v_ani = _mm_set1_ps(ani);
   for (k = nyi-1; k < nyt; k++) {
      joff = nxhd*k;
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
/*       t1 = f[j1+joff]; */
         v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&f[j1+joff]);
/*       f[j1+joff] = f[j+joff]; */
         v_t2 = _mm_loadl_pi(v_t2,(__m64 *)&f[j+joff]);
         _mm_storel_pi((__m64 *)&f[j1+joff],v_t2);
/*       f[j+joff] = t1; */
         _mm_storel_pi((__m64 *)&f[j+joff],v_t1);
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
         nss = 2*(ns/2);
         for (i = nyi-1; i < nyt; i++) {
            joff = nxhd*i;
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nss; j+=2) {
/*             t1 = conjf(sct[kmr*j]); */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_loadh_pi(v_t1,(__m64 *)&sct[kmr*j+kmr]);
               v_t1 = _mm_mul_ps(v_t1,v_n);
/*             t2 = t1*f[j2+joff]; */
               v_t2 = _mm_load_ps((float *)&f[j+k2+joff]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             f[j2+joff] = f[j1+joff] - t2; */
               v_t3 = _mm_load_ps((float *)&f[j+k1+joff]);
               _mm_store_ps((float *)&f[j+k2+joff],_mm_sub_ps(v_t3,v_t2));
/*             f[j1+joff] += t2; */
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
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2fft2rxy(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int nxi, 
                 int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
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
   requires SSE2, f needs to be 16 byte aligned
   nxhd needs to be a multiple of 2, and nxi needs to be odd
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
   int nxts;
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
   nxts = nxi + 2*(nxp/2) - 1;
   v_m = _mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);
   v_n = _mm_set_ps(-1.0f,1.0f,-1.0f,1.0f);
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
/* vector loop over elements in blocks of 2 */
      for (j = nxi-1; j < nxts; j+=2) {
/*       t1 = f[j+k1]; */
         v_t1 = _mm_load_ps((float *)&f[j+k1]);
/*       f[j+k1] = f[j+joff]; */
         v_t2 = _mm_load_ps((float *)&f[j+joff]);
         _mm_store_ps((float *)&f[j+k1],v_t2);
/*       f[j+joff] = t1; */
         _mm_store_ps((float *)&f[j+joff],v_t1);
      }
/* loop over remaining elements */
      for (j = nxts; j < nxt; j++) {
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
            v_t1 = _mm_set_ps(cimagf(t1),crealf(t1),cimagf(t1),crealf(t1));
/* vector loop over elements in blocks of 2 */
            for (i = nxi-1; i < nxts; i+=2) {
/*             t2 = t1*f[i+j2]; */
               v_t2 = _mm_load_ps((float *)&f[i+j2]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             f[i+j2] = f[i+j1] - t2; */
               v_t3 = _mm_load_ps((float *)&f[i+j1]);
               _mm_store_ps((float *)&f[i+j2],_mm_sub_ps(v_t3,v_t2));
/*             f[i+j1] += t2; */
               _mm_store_ps((float *)&f[i+j1],_mm_add_ps(v_t3,v_t2));
            }
/* loop over remaining elements */
            for (i = nxts; i < nxt; i++) {
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
/* vector loop over elements in blocks of 2 */
      for (j = nxi-1; j < nxts; j+=2) {
/*       t1 = f[j+k1]; */
         v_t1 = _mm_load_ps((float *)&f[j+k1]);
/*       f[j+k1] = f[j+joff]; */
         v_t2 = _mm_load_ps((float *)&f[j+joff]);
         _mm_store_ps((float *)&f[j+k1],v_t2);
/*       f[j+joff] = t1; */
         _mm_store_ps((float *)&f[j+joff],v_t1);
      }
/* loop over remaining elements */
      for (j = nxts; j < nxt; j++) {
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
            v_t1 = _mm_set_ps(cimagf(t1),crealf(t1),cimagf(t1),crealf(t1));
/* vector loop over elements in blocks of 2 */
            for (i = nxi-1; i < nxts; i+=2) {
/*             t2 = t1*f[i+j2]; */
               v_t2 = _mm_load_ps((float *)&f[i+j2]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             f[i+j2] = f[i+j1] - t2; */
               v_t3 = _mm_load_ps((float *)&f[i+j1]);
               _mm_store_ps((float *)&f[i+j2],_mm_sub_ps(v_t3,v_t2));
/*             f[i+j1] += t2; */
               _mm_store_ps((float *)&f[i+j1],_mm_add_ps(v_t3,v_t2));
            }
/* loop over remaining elements */
            for (i = nxts; i < nxt; i++) {
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
void csse2fft2r2x(float complex f[], int isign, int mixup[],
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
   requires SSE2, f needs to be 16 byte aligned
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, jj, j1, k1, k2, ns, ns2, km, kmr, joff;
   float ani;
/* float complex t1, t2; */
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
      goto L140;
/* inverse fourier transform */
/* swap complex components */
   for (k = nyi-1; k < nyt; k++) {
      joff = 2*nxhd*k;
      for (j = 0; j < nxh; j++) {
         v_t1 = _mm_load_ps((float *)&f[2*j+joff]);
         v_t1 = _mm_shuffle_ps(v_t1,v_t1,216);
         _mm_store_ps((float *)&f[2*j+joff],v_t1);
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
/*       t1 = f[2*j1+joff];   */
/*       t2 = f[1+2*j1+joff]; */
         v_t1 = _mm_load_ps((float *)&f[2*j1+joff]);
/*       f[2*j1+joff] = f[2*j+joff];     */
/*       f[1+2*j1+joff] = f[1+2*j+joff]; */
         v_t2 = _mm_load_ps((float *)&f[2*j+joff]);
         _mm_store_ps((float *)&f[2*j1+joff],v_t2);
/*       f[2*j+joff] = t1;   */
/*       f[1+2*j+joff] = t2; */
         _mm_store_ps((float *)&f[2*j+joff],v_t1);
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
      }
      ns = ns2;
   }
/* unscramble coefficients and normalize */
   kmr = nxy/nx;
   ani = 1.0/(float) (2*nx*ny);
   v_ani = _mm_set1_ps(ani);
   for (k = nyi-1; k < nyt; k++) {
      joff = 2*nxhd*k;
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
/*       t1 = f[2*j1+joff];   */
/*       t2 = f[1+2*j1+joff]; */
         v_t1 = _mm_load_ps((float *)&f[2*j1+joff]);
/*       f[2*j1+joff] = f[2*j+joff];     */
/*       f[1+2*j1+joff] = f[1+2*j+joff]; */
         v_t2 = _mm_load_ps((float *)&f[2*j+joff]);
         _mm_store_ps((float *)&f[2*j1+joff],v_t2);
/*       f[2*j+joff] = t1;   */
/*       f[1+2*j+joff] = t2; */
         _mm_store_ps((float *)&f[2*j+joff],v_t1);
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
/*             f[2*j+k2+joff] = f[2*j+k1+joff] - t2;     */
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
      }
      ns = ns2;
   }
/* swap complex components */
   for (k = nyi-1; k < nyt; k++) {
      joff = 2*nxhd*k;
      for (j = 0; j < nxh; j++) {
         v_t1 = _mm_load_ps((float *)&f[2*j+joff]);
         v_t1 = _mm_shuffle_ps(v_t1,v_t1,216);
         _mm_store_ps((float *)&f[2*j+joff],v_t1);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2fft2r2y(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxi,
                  int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
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
   requires SSE2, f needs to be 16 byte aligned
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
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
/*       t1 = f[2*j+k1];   */
/*       t2 = f[1+2*j+k1]; */
         v_t1 = _mm_load_ps((float *)&f[2*j+k1]);
/*       f[2*j+k1] = f[2*j+joff];     */
/*       f[1+2*j+k1] = f[1+2*j+joff]; */
         v_t2 = _mm_load_ps((float *)&f[2*j+joff]);
         _mm_store_ps((float *)&f[2*j+k1],v_t2);
/*       f[2*j+joff] = t1;   */
/*       f[1+2*j+joff] = t2; */
         _mm_store_ps((float *)&f[2*j+joff],v_t1);
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
/*          t1 = sct[kmr*j]; */
            v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
            v_t1 = _mm_movelh_ps(v_t1,v_t1);
            for (i = nxi-1; i < nxt; i++) {
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
/*       t1 = f[2*j+k1];   */
/*       t2 = f[1+2*j+k1]; */
         v_t1 = _mm_load_ps((float *)&f[2*j+k1]);
/*       f[2*j+k1] = f[2*j+joff];     */
/*       f[1+2*j+k1] = f[1+2*j+joff]; */
         v_t2 = _mm_load_ps((float *)&f[2*j+joff]);
         _mm_store_ps((float *)&f[2*j+k1],v_t2);
/*       f[2*j+joff] = t1;   */
/*       f[1+2*j+joff] = t2; */
         _mm_store_ps((float *)&f[2*j+joff],v_t1);
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
/*          t1 = conjf(sct[kmr*j]); */
            v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
            v_t1 = _mm_movelh_ps(v_t1,v_t1);
            v_t1 = _mm_mul_ps(v_t1,v_n);
            for (i = nxi-1; i < nxt; i++) {
/*             t2 = t1*f[2*i+j2];    */
/*             t3 = t1*f[1+2*i+j2];  */
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
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2wfft2rx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxhd,
                  int nyd, int nxhyd, int nxyhd) {
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
      csse2fft2rxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                   nxyhd);
/* perform y fft */
      csse2fft2rxy(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                   nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      csse2fft2rxy(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                   nxyhd);
/* perform x fft */
      csse2fft2rxx(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                   nxyhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2wfft2r2(float complex f[], int isign, int mixup[],
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
      csse2fft2r2x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                   nxyhd);
/* perform y fft */
      csse2fft2r2y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                   nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      csse2fft2r2y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                   nxyhd);
/* perform x fft */
      csse2fft2r2x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                   nxyhd);
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void csse2xiscan2_(int *isdata, int *nths) {
   csse2xiscan2(isdata,*nths);
   return;
}

/*--------------------------------------------------------------------*/
void csse2gpush2lt_(float *part, float *fxy, float *qbm, float *dt,
                    float *ek, int *idimp, int *nop, int *npe, int *nx,
                    int *ny, int *nxv, int *nyv, int *ipbc) {
   csse2gpush2lt(part,fxy,*qbm,*dt,ek,*idimp,*nop,*npe,*nx,*ny,*nxv,
                 *nyv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2gpost2lt_(float *part, float *q, float *qm, int *nop, int *npe,
                    int *idimp, int *nxv, int *nyv) {
   csse2gpost2lt(part,q,*qm,*nop,*npe,*idimp,*nxv,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void csse2dsortp2ylt_(float *parta, float *partb, int *npic, int *idimp,
                      int *nop, int *npe, int *ny1) {
   csse2dsortp2ylt(parta,partb,npic,*idimp,*nop,*npe,*ny1);
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
void csse2pois22_(float complex *q, float complex *fxy, int *isign,
                  float complex *ffc, float *ax, float *ay, float *affp,
                  float *we, int *nx, int *ny, int *nxvh, int *nyv,
                  int *nxhd, int *nyhd) {
   csse2pois22(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*nxvh,*nyv,*nxhd,
               *nyhd);
   return;
}

void csse2wfft2rx_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *nxhd,
                   int *nyd, int *nxhyd, int *nxyhd) {
   csse2wfft2rx(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,
                *nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csse2wfft2r2_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *nxhd,
                   int *nyd, int *nxhyd, int *nxyhd) {
   csse2wfft2r2(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,
                *nxyhd);
   return;
}
