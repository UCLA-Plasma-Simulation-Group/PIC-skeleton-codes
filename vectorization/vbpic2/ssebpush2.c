/* SE2 C Library for Skeleton 2-1/2D Electromagnetic Vector PIC Code */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <xmmintrin.h>
#include "ssebpush2.h"

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
void csse2gbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                     float dt, float dtc, float *ek, int idimp, int nop,
                     int npe, int nx, int ny, int nxv, int nyv,
                     int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field. Using the Boris Mover.
   vector version using guard cells
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
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = velocity vx of particle n
   part[3][n] = velocity vy of particle n
   part[4][n] = velocity vz of particle n
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
   npe = first dimension of particle array
   nx/ny = system length in x/y direction
   nxv = first dimension of field arrays, must be >= nx+1
   nyv = second dimension of field arrays, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   requires SSE2, part, fxy, and bxy need to be 16 byte aligned
   npe needs to be a multiple of 4, fxy, bxy need to have 4 components
local data                                                            */
   int j, nps, nn, mm, nm;
   float qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt;
   float anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, vx, vy, vz;
   double sum1;
   __m128i v_nxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qtmh, v_dt, v_dtc, v_one, v_two, v_half;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m128 v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 a, b, c, d, e, f, g, h;
   __m128d v_sum1, v_d;
   __attribute__((aligned(16))) unsigned int ll[4];
   __attribute__((aligned(16))) double dd[2];
   qtmh = 0.5f*qbm*dt;
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
   v_qtmh = _mm_set1_ps(qtmh);
   v_dt = _mm_set1_ps(dt);
   v_dtc = _mm_set1_ps(dtc);
   v_one = _mm_set1_ps(1.0f);
   v_two = _mm_set1_ps(2.0f);
   v_half = _mm_set1_ps(0.5f);
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
/*    nm = 4*(nn + nxv*mm); */
      v_it = _mm_mul_epu32(v_nxv,_mm_srli_si128(v_mm,4));
      v_mm = _mm_mul_epu32(v_mm,v_nxv);
      v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
      v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
/*    amx = 1.0f - dxp; */
/*    amy = 1.0f - dyp; */
      v_amx = _mm_sub_ps(v_one,v_dxp);
      v_amy = _mm_sub_ps(v_one,v_dyp);
      _mm_store_si128((__m128i *)ll,v_nn);
/* find electric field */
/*    nn = nm;                         */
/*    dx = amx*fxy[nn];                */
/*    dy = amx*fxy[nn+1];              */
/*    dz = amx*fxy[nn+2];              */
/*    mm = nn + 4;                     */
/*    dx = amy*(dxp*fxy[mm] + dx);     */
/*    dy = amy*(dxp*fxy[mm+1] + dy);   */
/*    dz = amy*(dxp*fxy[mm+2] + dz);   */
/*    nn += 4*nxv;                     */
/*    acx = amx*fxy[nn];               */
/*    acy = amx*fxy[nn+1];             */
/*    acz = amx*fxy[nn+2];             */
/*    mm = nn + 4;                     */
/*    dx += dyp*(dxp*fxy[mm] + acx);   */
/*    dy += dyp*(dxp*fxy[mm+1] + acy); */
/*    dz += dyp*(dxp*fxy[mm+2] + acz); */
/* find magnetic field */
/*    nn = nm;                           */
/*    ox = amx*bxy[nn];                  */
/*    oy = amx*bxy[nn+1];                */
/*    oz = amx*bxy[nn+2];                */
/*    mm = nn + 4;                       */
/*    ox = amy*(dxp*bxy[mm] + ox);       */
/*    oy = amy*(dxp*bxy[mm+1] + oy);     */
/*    oz = amy*(dxp*bxy[mm+2] + oz);     */
/*    nn += 4*nxv;                       */
/*    acx = amx*bxy[nn];                 */
/*    acy = amx*bxy[nn+1];               */
/*    acz = amx*bxy[nn+2];               */
/*    mm = nn + 4;                       */
/*    ox += dyp*(dxp*bxy[mm] + acx);     */
/*    oy += dyp*(dxp*bxy[mm+1] + acy);   */
/*    oz += dyp*(dxp*bxy[mm+2] + acz);   */
/* interpolate electric and magnetic fields for first particle */
      nn = ll[0];
      v_at = _mm_shuffle_ps(v_amx,v_amx,0);
      a = _mm_mul_ps(v_at,_mm_load_ps(&fxy[nn]));
      e = _mm_mul_ps(v_at,_mm_load_ps(&bxy[nn]));
      mm = nn + 4*nxv;
      v_dx = _mm_mul_ps(v_at,_mm_load_ps(&fxy[mm]));
      v_dy = _mm_mul_ps(v_at,_mm_load_ps(&bxy[mm]));
      v_at = _mm_shuffle_ps(v_dxp,v_dxp,0);
      nn = nn + 4;
      a = _mm_add_ps(a,_mm_mul_ps(v_at,_mm_load_ps(&fxy[nn])));
      e = _mm_add_ps(e,_mm_mul_ps(v_at,_mm_load_ps(&bxy[nn])));
      mm = mm + 4;
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&fxy[mm])));
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&bxy[mm])));
      v_at = _mm_shuffle_ps(v_amy,v_amy,0);
      a = _mm_mul_ps(a,v_at);
      e = _mm_mul_ps(e,v_at);
      v_at = _mm_shuffle_ps(v_dyp,v_dyp,0);
      a = _mm_add_ps(a,_mm_mul_ps(v_dx,v_at));
      e = _mm_add_ps(e,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for second particle */
      nn = ll[1];
      v_at = _mm_shuffle_ps(v_amx,v_amx,85);
      b = _mm_mul_ps(v_at,_mm_load_ps(&fxy[nn]));
      f = _mm_mul_ps(v_at,_mm_load_ps(&bxy[nn]));
      mm = nn + 4*nxv;
      v_dx = _mm_mul_ps(v_at,_mm_load_ps(&fxy[mm]));
      v_dy = _mm_mul_ps(v_at,_mm_load_ps(&bxy[mm]));
      v_at = _mm_shuffle_ps(v_dxp,v_dxp,85);
      nn = nn + 4;
      b = _mm_add_ps(b,_mm_mul_ps(v_at,_mm_load_ps(&fxy[nn])));
      f = _mm_add_ps(f,_mm_mul_ps(v_at,_mm_load_ps(&bxy[nn])));
      mm = mm + 4;
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&fxy[mm])));
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&bxy[mm])));
      v_at = _mm_shuffle_ps(v_amy,v_amy,85);
      b = _mm_mul_ps(b,v_at);
      f = _mm_mul_ps(f,v_at);
      v_at = _mm_shuffle_ps(v_dyp,v_dyp,85);
      b = _mm_add_ps(b,_mm_mul_ps(v_dx,v_at));
      f = _mm_add_ps(f,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for third particle */
      nn = ll[2];
      v_at = _mm_shuffle_ps(v_amx,v_amx,170);
      c = _mm_mul_ps(v_at,_mm_load_ps(&fxy[nn]));
      g = _mm_mul_ps(v_at,_mm_load_ps(&bxy[nn]));
      mm = nn + 4*nxv;
      v_dx = _mm_mul_ps(v_at,_mm_load_ps(&fxy[mm]));
      v_dy = _mm_mul_ps(v_at,_mm_load_ps(&bxy[mm]));
      v_at = _mm_shuffle_ps(v_dxp,v_dxp,170);
      nn = nn + 4;
      c = _mm_add_ps(c,_mm_mul_ps(v_at,_mm_load_ps(&fxy[nn])));
      g = _mm_add_ps(g,_mm_mul_ps(v_at,_mm_load_ps(&bxy[nn])));
      mm = mm + 4;
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&fxy[mm])));
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&bxy[mm])));
      v_at = _mm_shuffle_ps(v_amy,v_amy,170);
      c = _mm_mul_ps(c,v_at);
      g = _mm_mul_ps(g,v_at);
      v_at = _mm_shuffle_ps(v_dyp,v_dyp,170);
      c = _mm_add_ps(c,_mm_mul_ps(v_dx,v_at));
      g = _mm_add_ps(g,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for fourth particle */
      nn = ll[3];
      v_at = _mm_shuffle_ps(v_amx,v_amx,255);
      d = _mm_mul_ps(v_at,_mm_load_ps(&fxy[nn]));
      h = _mm_mul_ps(v_at,_mm_load_ps(&bxy[nn]));
      mm = nn + 4*nxv;
      v_dx = _mm_mul_ps(v_at,_mm_load_ps(&fxy[mm]));
      v_dy = _mm_mul_ps(v_at,_mm_load_ps(&bxy[mm]));
      v_at = _mm_shuffle_ps(v_dxp,v_dxp,255);
      nn = nn + 4;
      d = _mm_add_ps(d,_mm_mul_ps(v_at,_mm_load_ps(&fxy[nn])));
      h = _mm_add_ps(h,_mm_mul_ps(v_at,_mm_load_ps(&bxy[nn])));
      mm = mm + 4;
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&fxy[mm])));
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&bxy[mm])));
      v_at = _mm_shuffle_ps(v_amy,v_amy,255);
      d = _mm_mul_ps(d,v_at);
      h = _mm_mul_ps(h,v_at);
      v_at = _mm_shuffle_ps(v_dyp,v_dyp,255);
      d = _mm_add_ps(d,_mm_mul_ps(v_dx,v_at));
      h = _mm_add_ps(h,_mm_mul_ps(v_dy,v_at));
/* transpose so a,b,c,d contain electric fields for each particle */
      _MM_TRANSPOSE4_PS(a,b,c,d);
/* transpose so e,f,g,h contain magnetic fields for each particle */
      _MM_TRANSPOSE4_PS(e,f,g,h);
/* calculate half impulse */
/*    dx *= qtmh; */
/*    dy *= qtmh; */
/*    dz *= qtmh; */
      v_dx = _mm_mul_ps(a,v_qtmh);
      v_dy = _mm_mul_ps(b,v_qtmh);
      v_dz = _mm_mul_ps(c,v_qtmh);
/* half acceleration */
/*    acx = part[j+2*npe] + dx; */
/*    acy = part[j+3*npe] + dy; */
/*    acz = part[j+4*npe] + dz; */
      a = _mm_add_ps(v_dx,_mm_load_ps(&part[j+2*npe]));
      b = _mm_add_ps(v_dy,_mm_load_ps(&part[j+3*npe]));
      c = _mm_add_ps(v_dz,_mm_load_ps(&part[j+4*npe]));
/* time-centered kinetic energy */
/*    sum1 += (acx*acx + acy*acy + acz*acz); */
      v_at = _mm_add_ps(_mm_mul_ps(a,a),_mm_mul_ps(b,b));
      v_at = _mm_add_ps(v_at,_mm_mul_ps(c,c));
/* convert to double precision before accumulating */
      v_d = _mm_cvtps_pd(v_at);
      v_sum1 = _mm_add_pd(v_sum1,v_d);
      v_it = _mm_srli_si128((__m128i)v_at,8);
      v_d = _mm_cvtps_pd((__m128)v_it);
      v_sum1 = _mm_add_pd(v_sum1,v_d);
/* calculate cyclotron frequency */
/*    omxt = qtmh*ox; */
/*    omyt = qtmh*oy; */
/*    omzt = qtmh*oz; */
      e = _mm_mul_ps(v_qtmh,e);
      f = _mm_mul_ps(v_qtmh,f);
      g = _mm_mul_ps(v_qtmh,g);
/* calculate rotation matrix */
/*    vx = omxt*omxt; */
      v_vx = _mm_mul_ps(e,e);
/*    vy = omyt*omyt; */
      v_vy = _mm_mul_ps(f,f);
/*    vz = omzt*omzt; */
      v_vz = _mm_mul_ps(g,g);
/*    omt = omxt*omxt + omyt*omyt + omzt*omzt; */
      v_at = _mm_add_ps(_mm_add_ps(v_vx,v_vy),v_vz);
/*    anorm = 2.0f/(1.0f + omt); */
      d = _mm_div_ps(v_two,_mm_add_ps(v_one,v_at));
/*    omt = 0.5f*(1.0f - omt); */
      h = _mm_mul_ps(v_half,_mm_sub_ps(v_one,v_at));
/*    vx = (omt + vx)*acx; */
      v_vx = _mm_mul_ps(_mm_add_ps(h,v_vx),a);
/*    vy = (omt + vy)*acy; */
      v_vy = _mm_mul_ps(_mm_add_ps(h,v_vy),b);
/*    vz = (omt + vz)*acz; */
      v_vz = _mm_mul_ps(_mm_add_ps(h,v_vz),c);
/*    omt = omxt*omyt; */
      h = _mm_mul_ps(e,f);
/*    vx = vx + (omzt + omt)*acy; */
      v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_add_ps(h,g),b));
/*    vy = vy + (omt - omzt)*acx; */
      v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_sub_ps(h,g),a));
/*    omt = omxt*omzt;  */
      h = _mm_mul_ps(e,g);
/*    vx = vx + (omt - omyt)*acz; */
      v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_sub_ps(h,f),c));
/*    vz = vz + (omt + omyt)*acx; */
      v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_add_ps(h,f),a));
/*    omt = omyt*omzt; */
      h = _mm_mul_ps(f,g);
/*    vy = vy + (omt + omxt)*acz; */
      v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_add_ps(h,e),c));
/*    vz = vz + (omt - omxt)*acy; */
      v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_sub_ps(h,e),b));
/* new momentum */
/*    vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm; */
/*    vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm; */
/*    vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm; */
      v_vx = _mm_add_ps(v_dx,_mm_mul_ps(v_vx,d));
      v_vy = _mm_add_ps(v_dy,_mm_mul_ps(v_vy,d));
      v_vz = _mm_add_ps(v_dz,_mm_mul_ps(v_vz,d));
/* new position */
/*    dx = x + vx*dtc; */
/*    dy = y + vy*dtc; */
      v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_dtc));
      v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_dtc));
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
/* set new momemtum */
/*    part[j+2*npe] = vx; */
/*    part[j+3*npe] = vy; */
/*    part[j+4*npe] = vz; */
      _mm_store_ps(&part[j+2*npe],v_vx);
      _mm_store_ps(&part[j+3*npe],v_vy);
      _mm_store_ps(&part[j+4*npe],v_vz);
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
      nm = 4*(nn + nxv*mm);
      amx = 1.0f - dxp;
      amy = 1.0f - dyp;
/* find electric field */
      nn = nm;
      dx = amx*fxy[nn];
      dy = amx*fxy[nn+1];
      dz = amx*fxy[nn+2];
      mm = nn + 4;
      dx = amy*(dxp*fxy[mm] + dx);
      dy = amy*(dxp*fxy[mm+1] + dy);
      dz = amy*(dxp*fxy[mm+2] + dz);
      nn += 4*nxv;
      acx = amx*fxy[nn];
      acy = amx*fxy[nn+1];
      acz = amx*fxy[nn+2];
      mm = nn + 4;
      dx += dyp*(dxp*fxy[mm] + acx);
      dy += dyp*(dxp*fxy[mm+1] + acy);
      dz += dyp*(dxp*fxy[mm+2] + acz);
/* find magnetic field */
      nn = nm;
      ox = amx*bxy[nn];
      oy = amx*bxy[nn+1];
      oz = amx*bxy[nn+2];
      mm = nn + 4;
      ox = amy*(dxp*bxy[mm] + ox);
      oy = amy*(dxp*bxy[mm+1] + oy);
      oz = amy*(dxp*bxy[mm+2] + oz);
      nn += 4*nxv;
      acx = amx*bxy[nn];
      acy = amx*bxy[nn+1];
      acz = amx*bxy[nn+2];
      mm = nn + 4;
      ox += dyp*(dxp*bxy[mm] + acx);
      oy += dyp*(dxp*bxy[mm+1] + acy);
      oz += dyp*(dxp*bxy[mm+2] + acz);
/* calculate half impulse */
      dx *= qtmh;
      dy *= qtmh;
      dz *= qtmh;
/* half acceleration */
      acx = part[j+2*npe] + dx;
      acy = part[j+3*npe] + dy;
      acz = part[j+4*npe] + dz;
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
      vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm;
      vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm;
      vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm;
/* new position */
      dx = x + vx*dtc;
      dy = y + vy*dtc;
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
      part[j+4*npe] = vz;
   }

/* normalize kinetic energy */
/* *ek += 0.5f*sum1; */
   _mm_store_pd(&dd[0],v_sum1);
   for (j = 1; j < 2; j++) {
      dd[0] += dd[j];
   }
   *ek += 0.5f*(sum1 + dd[0]);
   return;
}

/*--------------------------------------------------------------------*/
void csse2grbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                      float dt, float dtc, float ci, float *ek,
                      int idimp, int nop, int npe, int nx, int ny,
                      int nxv, int nyv, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   vector version using guard cells
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
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = momentum px of particle n
   part[3][n] = momentum py of particle n
   part[4][n] = momentum pz of particle n
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
   npe = first dimension of particle array
   nx/ny = system length in x/y direction
   nxv = first dimension of field arrays, must be >= nx+1
   nyv = second dimension of field arrays, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   requires SSE2, part, fxy, and bxy need to be 16 byte aligned
   npe needs to be a multiple of 4, fxy, bxy need to have 4 components
local data                                                            */
   int j, nps, nn, mm, nm;
   float qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg;
   float omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, vx, vy, vz;
   double sum1;
   __m128i v_nxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qtmh, v_ci2, v_dt, v_dtc, v_one, v_two, v_half;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_gami, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m128 v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 a, b, c, d, e, f, g, h;
   __m128d v_sum1, v_d;
   __attribute__((aligned(16))) unsigned int ll[4];
   __attribute__((aligned(16))) double dd[2];
   qtmh = 0.5f*qbm*dt;
   ci2 = ci*ci;
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
   v_qtmh = _mm_set1_ps(qtmh);
   v_ci2 = _mm_set1_ps(ci2);
   v_dt = _mm_set1_ps(dt);
   v_dtc = _mm_set1_ps(dtc);
   v_one = _mm_set1_ps(1.0f);
   v_two = _mm_set1_ps(2.0f);
   v_half = _mm_set1_ps(0.5f);
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
/*    nm = 4*(nn + nxv*mm); */
      v_it = _mm_mul_epu32(v_nxv,_mm_srli_si128(v_mm,4));
      v_mm = _mm_mul_epu32(v_mm,v_nxv);
      v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
      v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
/*    amx = 1.0f - dxp; */
/*    amy = 1.0f - dyp; */
      v_amx = _mm_sub_ps(v_one,v_dxp);
      v_amy = _mm_sub_ps(v_one,v_dyp);
      _mm_store_si128((__m128i *)ll,v_nn);
/* find electric field */
/*    nn = nm;                         */
/*    dx = amx*fxy[nn];                */
/*    dy = amx*fxy[nn+1];              */
/*    dz = amx*fxy[nn+2];              */
/*    mm = nn + 4;                     */
/*    dx = amy*(dxp*fxy[mm] + dx);     */
/*    dy = amy*(dxp*fxy[mm+1] + dy);   */
/*    dz = amy*(dxp*fxy[mm+2] + dz);   */
/*    nn += 4*nxv;                     */
/*    acx = amx*fxy[nn];               */
/*    acy = amx*fxy[nn+1];             */
/*    acz = amx*fxy[nn+2];             */
/*    mm = nn + 4;                     */
/*    dx += dyp*(dxp*fxy[mm] + acx);   */
/*    dy += dyp*(dxp*fxy[mm+1] + acy); */
/*    dz += dyp*(dxp*fxy[mm+2] + acz); */
/* find magnetic field */
/*    nn = nm;                           */
/*    ox = amx*bxy[nn];                  */
/*    oy = amx*bxy[nn+1];                */
/*    oz = amx*bxy[nn+2];                */
/*    mm = nn + 4;                       */
/*    ox = amy*(dxp*bxy[mm] + ox);       */
/*    oy = amy*(dxp*bxy[mm+1] + oy);     */
/*    oz = amy*(dxp*bxy[mm+2] + oz);     */
/*    nn += 4*nxv;                       */
/*    acx = amx*bxy[nn];                 */
/*    acy = amx*bxy[nn+1];               */
/*    acz = amx*bxy[nn+2];               */
/*    mm = nn + 4;                       */
/*    ox += dyp*(dxp*bxy[mm] + acx);     */
/*    oy += dyp*(dxp*bxy[mm+1] + acy);   */
/*    oz += dyp*(dxp*bxy[mm+2] + acz);   */
/* interpolate electric and magnetic fields for first particle */
      nn = ll[0];
      v_at = _mm_shuffle_ps(v_amx,v_amx,0);
      a = _mm_mul_ps(v_at,_mm_load_ps(&fxy[nn]));
      e = _mm_mul_ps(v_at,_mm_load_ps(&bxy[nn]));
      mm = nn + 4*nxv;
      v_dx = _mm_mul_ps(v_at,_mm_load_ps(&fxy[mm]));
      v_dy = _mm_mul_ps(v_at,_mm_load_ps(&bxy[mm]));
      v_at = _mm_shuffle_ps(v_dxp,v_dxp,0);
      nn = nn + 4;
      a = _mm_add_ps(a,_mm_mul_ps(v_at,_mm_load_ps(&fxy[nn])));
      e = _mm_add_ps(e,_mm_mul_ps(v_at,_mm_load_ps(&bxy[nn])));
      mm = mm + 4;
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&fxy[mm])));
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&bxy[mm])));
      v_at = _mm_shuffle_ps(v_amy,v_amy,0);
      a = _mm_mul_ps(a,v_at);
      e = _mm_mul_ps(e,v_at);
      v_at = _mm_shuffle_ps(v_dyp,v_dyp,0);
      a = _mm_add_ps(a,_mm_mul_ps(v_dx,v_at));
      e = _mm_add_ps(e,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for second particle */
      nn = ll[1];
      v_at = _mm_shuffle_ps(v_amx,v_amx,85);
      b = _mm_mul_ps(v_at,_mm_load_ps(&fxy[nn]));
      f = _mm_mul_ps(v_at,_mm_load_ps(&bxy[nn]));
      mm = nn + 4*nxv;
      v_dx = _mm_mul_ps(v_at,_mm_load_ps(&fxy[mm]));
      v_dy = _mm_mul_ps(v_at,_mm_load_ps(&bxy[mm]));
      v_at = _mm_shuffle_ps(v_dxp,v_dxp,85);
      nn = nn + 4;
      b = _mm_add_ps(b,_mm_mul_ps(v_at,_mm_load_ps(&fxy[nn])));
      f = _mm_add_ps(f,_mm_mul_ps(v_at,_mm_load_ps(&bxy[nn])));
      mm = mm + 4;
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&fxy[mm])));
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&bxy[mm])));
      v_at = _mm_shuffle_ps(v_amy,v_amy,85);
      b = _mm_mul_ps(b,v_at);
      f = _mm_mul_ps(f,v_at);
      v_at = _mm_shuffle_ps(v_dyp,v_dyp,85);
      b = _mm_add_ps(b,_mm_mul_ps(v_dx,v_at));
      f = _mm_add_ps(f,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for third particle */
      nn = ll[2];
      v_at = _mm_shuffle_ps(v_amx,v_amx,170);
      c = _mm_mul_ps(v_at,_mm_load_ps(&fxy[nn]));
      g = _mm_mul_ps(v_at,_mm_load_ps(&bxy[nn]));
      mm = nn + 4*nxv;
      v_dx = _mm_mul_ps(v_at,_mm_load_ps(&fxy[mm]));
      v_dy = _mm_mul_ps(v_at,_mm_load_ps(&bxy[mm]));
      v_at = _mm_shuffle_ps(v_dxp,v_dxp,170);
      nn = nn + 4;
      c = _mm_add_ps(c,_mm_mul_ps(v_at,_mm_load_ps(&fxy[nn])));
      g = _mm_add_ps(g,_mm_mul_ps(v_at,_mm_load_ps(&bxy[nn])));
      mm = mm + 4;
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&fxy[mm])));
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&bxy[mm])));
      v_at = _mm_shuffle_ps(v_amy,v_amy,170);
      c = _mm_mul_ps(c,v_at);
      g = _mm_mul_ps(g,v_at);
      v_at = _mm_shuffle_ps(v_dyp,v_dyp,170);
      c = _mm_add_ps(c,_mm_mul_ps(v_dx,v_at));
      g = _mm_add_ps(g,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for fourth particle */
      nn = ll[3];
      v_at = _mm_shuffle_ps(v_amx,v_amx,255);
      d = _mm_mul_ps(v_at,_mm_load_ps(&fxy[nn]));
      h = _mm_mul_ps(v_at,_mm_load_ps(&bxy[nn]));
      mm = nn + 4*nxv;
      v_dx = _mm_mul_ps(v_at,_mm_load_ps(&fxy[mm]));
      v_dy = _mm_mul_ps(v_at,_mm_load_ps(&bxy[mm]));
      v_at = _mm_shuffle_ps(v_dxp,v_dxp,255);
      nn = nn + 4;
      d = _mm_add_ps(d,_mm_mul_ps(v_at,_mm_load_ps(&fxy[nn])));
      h = _mm_add_ps(h,_mm_mul_ps(v_at,_mm_load_ps(&bxy[nn])));
      mm = mm + 4;
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&fxy[mm])));
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&bxy[mm])));
      v_at = _mm_shuffle_ps(v_amy,v_amy,255);
      d = _mm_mul_ps(d,v_at);
      h = _mm_mul_ps(h,v_at);
      v_at = _mm_shuffle_ps(v_dyp,v_dyp,255);
      d = _mm_add_ps(d,_mm_mul_ps(v_dx,v_at));
      h = _mm_add_ps(h,_mm_mul_ps(v_dy,v_at));
/* transpose so a,b,c,d contain electric fields for each particle */
      _MM_TRANSPOSE4_PS(a,b,c,d);
/* transpose so e,f,g,h contain magnetic fields for each particle */
      _MM_TRANSPOSE4_PS(e,f,g,h);
/* calculate half impulse */
/*    dx *= qtmh; */
/*    dy *= qtmh; */
/*    dz *= qtmh; */
      v_dx = _mm_mul_ps(a,v_qtmh);
      v_dy = _mm_mul_ps(b,v_qtmh);
      v_dz = _mm_mul_ps(c,v_qtmh);
/* half acceleration */
/*    acx = part[j+2*npe] + dx; */
/*    acy = part[j+3*npe] + dy; */
/*    acz = part[j+4*npe] + dz; */
      a = _mm_add_ps(v_dx,_mm_load_ps(&part[j+2*npe]));
      b = _mm_add_ps(v_dy,_mm_load_ps(&part[j+3*npe]));
      c = _mm_add_ps(v_dz,_mm_load_ps(&part[j+4*npe]));
/* find inverse gamma */
/*    p2 = acx*acx + acy*acy + acz*acz; */
      v_at = _mm_add_ps(_mm_mul_ps(a,a),_mm_mul_ps(b,b));
      v_at = _mm_add_ps(v_at,_mm_mul_ps(c,c));
/*    gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*    v_gami = _mm_rsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* full accuracy calculation */
      v_gami = _mm_sqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2)));
      v_gami = _mm_div_ps(v_one,v_gami);
/* full accuracy calculation with SVML */
/*    v_gami = _mm_invsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* time-centered kinetic energy */
/*    sum1 += gami*p2/(1.0f + gami); */
      v_at = _mm_mul_ps(v_gami,v_at);
      v_at = _mm_div_ps(v_at,_mm_add_ps(v_one,v_gami));
/* convert to double precision before accumulating */
      v_d = _mm_cvtps_pd(v_at);
      v_sum1 = _mm_add_pd(v_sum1,v_d);
      v_it = _mm_srli_si128((__m128i)v_at,8);
      v_d = _mm_cvtps_pd((__m128)v_it);
      v_sum1 = _mm_add_pd(v_sum1,v_d);
/* renormalize magnetic field */
/*    qtmg = qtmh*gami; */
      v_at = _mm_mul_ps(v_qtmh,v_gami);
/* calculate cyclotron frequency */
/*    omxt = qtmg*ox; */
/*    omyt = qtmg*oy; */
/*    omzt = qtmg*oz; */
      e = _mm_mul_ps(v_at,e);
      f = _mm_mul_ps(v_at,f);
      g = _mm_mul_ps(v_at,g);
/* calculate rotation matrix */
/*    vx = omxt*omxt; */
      v_vx = _mm_mul_ps(e,e);
/*    vy = omyt*omyt; */
      v_vy = _mm_mul_ps(f,f);
/*    vz = omzt*omzt; */
      v_vz = _mm_mul_ps(g,g);
/*    omt = omxt*omxt + omyt*omyt + omzt*omzt; */
      v_at = _mm_add_ps(_mm_add_ps(v_vx,v_vy),v_vz);
/*    anorm = 2.0f/(1.0f + omt); */
      d = _mm_div_ps(v_two,_mm_add_ps(v_one,v_at));
/*    omt = 0.5f*(1.0f - omt); */
      h = _mm_mul_ps(v_half,_mm_sub_ps(v_one,v_at));
/*    vx = (omt + vx)*acx; */
      v_vx = _mm_mul_ps(_mm_add_ps(h,v_vx),a);
/*    vy = (omt + vy)*acy; */
      v_vy = _mm_mul_ps(_mm_add_ps(h,v_vy),b);
/*    vz = (omt + vz)*acz; */
      v_vz = _mm_mul_ps(_mm_add_ps(h,v_vz),c);
/*    omt = omxt*omyt; */
      h = _mm_mul_ps(e,f);
/*    vx = vx + (omzt + omt)*acy; */
      v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_add_ps(h,g),b));
/*    vy = vy + (omt - omzt)*acx; */
      v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_sub_ps(h,g),a));
/*    omt = omxt*omzt;  */
      h = _mm_mul_ps(e,g);
/*    vx = vx + (omt - omyt)*acz; */
      v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_sub_ps(h,f),c));
/*    vz = vz + (omt + omyt)*acx; */
      v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_add_ps(h,f),a));
/*    omt = omyt*omzt; */
      h = _mm_mul_ps(f,g);
/*    vy = vy + (omt + omxt)*acz; */
      v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_add_ps(h,e),c));
/*    vz = vz + (omt - omxt)*acy; */
      v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_sub_ps(h,e),b));
/* new momentum */
/*    vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm; */
/*    vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm; */
/*    vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm; */
      v_vx = _mm_add_ps(v_dx,_mm_mul_ps(v_vx,d));
      v_vy = _mm_add_ps(v_dy,_mm_mul_ps(v_vy,d));
      v_vz = _mm_add_ps(v_dz,_mm_mul_ps(v_vz,d));
/* update inverse gamma */
/*    p2 = vx*vx + vy*vy + vz*vz; */
      v_at = _mm_mul_ps(v_vx,v_vx);
      v_at = _mm_add_ps(v_at,_mm_mul_ps(v_vy,v_vy));
      v_at = _mm_add_ps(v_at,_mm_mul_ps(v_vz,v_vz));
/*    dtg = dtc/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*    v_at = _mm_rsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/*    v_at = _mm_mul_ps(v_dtc,v_at); */
/* full accuracy calculation */
      v_at = _mm_sqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2)));
      v_at = _mm_div_ps(v_dtc,v_at);
/* full accuracy calculation with SVML */
/*    v_at = _mm_invsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/*    v_at = _mm_mul_ps(v_dtc,v_at); */
/* new position */
/*    dx = x + vx*dtg; */
/*    dy = y + vy*dtg; */
      v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_at));
      v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_at));
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
/*       } */
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
/*    part[j+4*npe] = vz; */
      _mm_store_ps(&part[j+2*npe],v_vx);
      _mm_store_ps(&part[j+3*npe],v_vy);
      _mm_store_ps(&part[j+4*npe],v_vz);
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
      nm = 4*(nn + nxv*mm);
      amx = 1.0f - dxp;
      amy = 1.0f - dyp;
/* find electric field */
      nn = nm;
      dx = amx*fxy[nn];
      dy = amx*fxy[nn+1];
      dz = amx*fxy[nn+2];
      mm = nn + 4;
      dx = amy*(dxp*fxy[mm] + dx);
      dy = amy*(dxp*fxy[mm+1] + dy);
      dz = amy*(dxp*fxy[mm+2] + dz);
      nn += 4*nxv;
      acx = amx*fxy[nn];
      acy = amx*fxy[nn+1];
      acz = amx*fxy[nn+2];
      mm = nn + 4;
      dx += dyp*(dxp*fxy[mm] + acx);
      dy += dyp*(dxp*fxy[mm+1] + acy);
      dz += dyp*(dxp*fxy[mm+2] + acz);
/* find magnetic field */
      nn = nm;
      ox = amx*bxy[nn];
      oy = amx*bxy[nn+1];
      oz = amx*bxy[nn+2];
      mm = nn + 4;
      ox = amy*(dxp*bxy[mm] + ox);
      oy = amy*(dxp*bxy[mm+1] + oy);
      oz = amy*(dxp*bxy[mm+2] + oz);
      nn += 4*nxv;
      acx = amx*bxy[nn];
      acy = amx*bxy[nn+1];
      acz = amx*bxy[nn+2];
      mm = nn + 4;
      ox += dyp*(dxp*bxy[mm] + acx);
      oy += dyp*(dxp*bxy[mm+1] + acy);
      oz += dyp*(dxp*bxy[mm+2] + acz);
/* calculate half impulse */
      dx *= qtmh;
      dy *= qtmh;
      dz *= qtmh;
/* half acceleration */
      acx = part[j+2*npe] + dx;
      acy = part[j+3*npe] + dy;
      acz = part[j+4*npe] + dz;
/* find inverse gamma */
      p2 = acx*acx + acy*acy + acz*acz;
      gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* renormalize magnetic field */
      qtmg = qtmh*gami;
/* time-centered kinetic energy */
      sum1 += gami*p2/(1.0f + gami);
/* calculate cyclotron frequency */
      omxt = qtmg*ox;
      omyt = qtmg*oy;
      omzt = qtmg*oz;
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
/* new momentum */
      vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm;
      vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm;
      vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm;
/* update inverse gamma */
      p2 = vx*vx + vy*vy + vz*vz;
      dtg = dtc/sqrtf(1.0f + p2*ci2);
/* new position */
      dx = x + vx*dtg;
      dy = y + vy*dtg;
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
      part[j+4*npe] = vz;
   }
/* normalize kinetic energy */
/* *ek += sum1; */
   _mm_store_pd(&dd[0],v_sum1);
   for (j = 1; j < 2; j++) {
      dd[0] += dd[j];
   }
   *ek += sum1 + dd[0];
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
void csse2gjpost2lt(float part[], float cu[], float qm, float dt,
                    int nop, int npe, int idimp, int nx, int ny, 
                    int nxv, int nyv, int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   vector version using guard cells
   41 flops/particle, 17 loads, 14 stores
   input: all, output: part, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m)=qci*(1.-dx)*(1.-dy)
   cu(i,n+1,m)=qci*dx*(1.-dy)
   cu(i,n,m+1)=qci*(1.-dx)*dy
   cu(i,n+1,m+1)=qci*dx*dy
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   and qci = qm*vi, where i = x,y,z
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = x velocity of particle n
   part[3][n] = y velocity of particle n
   part[4][n] = z velocity of particle n
   cu[k][j][i] = ith component of current density at grid point j,k
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nop = number of particles
   npe = first dimension of particle array
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   nxv = first dimension of current array, must be >= nx+1
   nyv = second dimension of current array, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   requires SSE2, part and cu need to be 16 byte aligned
   npe needs to be a multiple of 4, cu needs to have 4 components
local data                                                            */
   int j, nps, nn, mm;
   float edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz;
   __m128i v_nxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qm, v_dt, v_one;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_vx, v_vy;
   __m128 v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 a, b, c, d, va, vb, vc, vd;
   __attribute__((aligned(16))) unsigned int ll[4];
   __attribute__((aligned(16))) unsigned long kk[1];
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
   v_qm = _mm_set1_ps(qm);
   v_one = _mm_set1_ps(1.0f);
   v_dt = _mm_set1_ps(dt);
   v_edgelx = _mm_set1_ps(edgelx);
   v_edgely = _mm_set1_ps(edgely);
   v_edgerx = _mm_set1_ps(edgerx);
   v_edgery = _mm_set1_ps(edgery);
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
/*    dyp = y - (float) mm;      */
      v_dyp = _mm_sub_ps(v_y,_mm_cvtepi32_ps(v_mm));
/*    nn = 4*(nn + nxv*mm); */
      v_it = _mm_mul_epu32(v_nxv,_mm_srli_si128(v_mm,4));
      v_mm = _mm_mul_epu32(v_mm,v_nxv);
      v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
      v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
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
/* deposit current */
/*    vx = part[j+2*npe]; */
/*    vy = part[j+3*npe]; */
/*    vz = part[j+4*npe]; */
      v_vx = _mm_load_ps(&part[j+2*npe]);
      v_vy = _mm_load_ps(&part[j+3*npe]);
      va = v_vx;
      vb = v_vy;
      vc = _mm_load_ps(&part[j+4*npe]);
      vd = _mm_setzero_ps();
/* transpose so va,vb,vc,vd contain the 3 velocities plus zero */
/* for each of 4 particles                                     */
      _MM_TRANSPOSE4_PS(va,vb,vc,vd);
/*    dx = amx*amy;      */
/*    cu[nn] += vx*dx;   */
/*    cu[nn+1] += vy*dx; */
/*    cu[nn+2] += vz*dx; */
/*    dy = dxp*amy;      */
/*    mm = nn + 4;       */
/*    cu[mm] += vx*dy;   */
/*    cu[mm+1] += vy*dy; */
/*    cu[mm+2] += vz*dy; */
/*    dx = amx*dyp;      */
/*    nn += 4*nxv;       */
/*    cu[nn] += vx*dx;   */
/*    cu[nn+1] += vy*dx; */
/*    cu[nn+2] += vz*dx; */
/*    dy = dxp*dyp;      */
/*    mm = nn + 4;       */
/*    cu[mm] += vx*dy;   */
/*    cu[mm+1] += vy*dy; */
/*    cu[mm+2] += vz*dy; */
/* deposit for first particle */
      mm = ll[0];
      v_dx = _mm_load_ps(&cu[mm]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(a,a,0)));
      _mm_store_ps(&cu[mm],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(b,b,0)));
      _mm_store_ps(&cu[mm+4],v_dy);
      v_dx = _mm_load_ps(&cu[mm+4*nxv]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(c,c,0)));
      _mm_store_ps(&cu[mm+4*nxv],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4+4*nxv]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(d,d,0)));
      _mm_store_ps(&cu[mm+4+4*nxv],v_dy);
/* deposit for second particle */
      mm = ll[1];
      v_dx = _mm_load_ps(&cu[mm]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(a,a,85)));
      _mm_store_ps(&cu[mm],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(b,b,85)));
      _mm_store_ps(&cu[mm+4],v_dy);
      v_dx = _mm_load_ps(&cu[mm+4*nxv]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(c,c,85)));
      _mm_store_ps(&cu[mm+4*nxv],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4+4*nxv]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(d,d,85)));
      _mm_store_ps(&cu[mm+4+4*nxv],v_dy);
/* deposit for third particle */
      mm = ll[2];
      v_dx = _mm_load_ps(&cu[mm]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(a,a,170)));
      _mm_store_ps(&cu[mm],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(b,b,170)));
      _mm_store_ps(&cu[mm+4],v_dy);
      v_dx = _mm_load_ps(&cu[mm+4*nxv]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(c,c,170)));
      _mm_store_ps(&cu[mm+4*nxv],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4+4*nxv]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(d,d,170)));
      _mm_store_ps(&cu[mm+4+4*nxv],v_dy);
/* deposit for fourth particle */
      mm = ll[3];
      v_dx = _mm_load_ps(&cu[mm]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(a,a,255)));
      _mm_store_ps(&cu[mm],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(b,b,255)));
      _mm_store_ps(&cu[mm+4],v_dy);
      v_dx = _mm_load_ps(&cu[mm+4*nxv]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(c,c,255)));
      _mm_store_ps(&cu[mm+4*nxv],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4+4*nxv]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(d,d,255)));
      _mm_store_ps(&cu[mm+4+4*nxv],v_dy);
/* advance position half a time-step */
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
/*          part[j+2*npe] = -vx;                */
/*       }                                      */
         v_at = _mm_cmplt_ps(v_dx,v_edgelx);
         v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dx,v_edgerx));
         v_x = _mm_and_ps(v_at,v_x);
         v_dx = _mm_add_ps(_mm_andnot_ps(v_at,v_dx),v_x);
         v_dxp = _mm_and_ps(v_at,v_vx);
         v_vx = _mm_sub_ps(_mm_andnot_ps(v_at,v_vx),v_dxp);
/* write output if test result is true for any particle */
         v_mm = _mm_srli_si128((__m128i)v_at,8);
         v_mm = _mm_add_epi64((__m128i)v_at,v_mm);
         _mm_storel_epi64((__m128i *)&kk[0],v_mm);
         if (kk[0] != 0)
            _mm_store_ps(&part[j+2*npe],v_vx);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          part[j+3*npe] = -vy;                */
/*       }                                      */
         v_at = _mm_cmplt_ps(v_dy,v_edgely);
         v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dy,v_edgery));
         v_y = _mm_and_ps(v_at,v_y);
         v_dy = _mm_add_ps(_mm_andnot_ps(v_at,v_dy),v_y);
         v_dyp = _mm_and_ps(v_at,v_vy);
         v_vy = _mm_sub_ps(_mm_andnot_ps(v_at,v_vy),v_dyp);
/* write output if test result is true for any particle */
         v_mm = _mm_srli_si128((__m128i)v_at,8);
         v_mm = _mm_add_epi64((__m128i)v_at,v_mm);
         _mm_storel_epi64((__m128i *)&kk[0],v_mm);
         if (kk[0] != 0)
            _mm_store_ps(&part[j+3*npe],v_vy);
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          part[j+2*npe] = -vx;                */
/*       }                                      */
         v_at = _mm_cmplt_ps(v_dx,v_edgelx);
         v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dx,v_edgerx));
         v_x = _mm_and_ps(v_at,v_x);
         v_dx = _mm_add_ps(_mm_andnot_ps(v_at,v_dx),v_x);
         v_dxp = _mm_and_ps(v_at,v_vx);
         v_vx = _mm_sub_ps(_mm_andnot_ps(v_at,v_vx),v_dxp);
/* write output if test result is true for any particle */
         v_mm = _mm_srli_si128((__m128i)v_at,8);
         v_mm = _mm_add_epi64((__m128i)v_at,v_mm);
         _mm_storel_epi64((__m128i *)&kk[0],v_mm);
         if (kk[0] != 0)
            _mm_store_ps(&part[j+2*npe],v_vx);
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
      nn = 4*(nn + nxv*mm);
      amx = qm - dxp;
      amy = 1.0f - dyp;
/* deposit current */
      dx = amx*amy;
      dy = dxp*amy;
      vx = part[j+2*npe];
      vy = part[j+3*npe];
      vz = part[j+4*npe];
      cu[nn] += vx*dx;
      cu[nn+1] += vy*dx;
      cu[nn+2] += vz*dx;
      dx = amx*dyp;
      mm = nn + 4;
      cu[mm] += vx*dy;
      cu[mm+1] += vy*dy;
      cu[mm+2] += vz*dy;
      dy = dxp*dyp;
      nn += 4*nxv;
      cu[nn] += vx*dx;
      cu[nn+1] += vy*dx;
      cu[nn+2] += vz*dx;
      mm = nn + 4;
      cu[mm] += vx*dy;
      cu[mm+1] += vy*dy;
      cu[mm+2] += vz*dy;
/* advance position half a time-step */
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
            part[j+2*npe] = -vx;
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = y;
            part[j+3*npe] = -vy;
         }
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            part[j+2*npe] = -vx;
         }
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* set new position */
      part[j] = dx;
      part[j+npe] = dy;
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2grjpost2lt(float part[], float cu[], float qm, float dt,
                     float ci, int nop, int npe, int idimp, int nx,
                     int ny, int nxv, int nyv, int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation for relativistic particles
   in addition, particle positions are advanced a half time-step
   vector version using guard cells
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
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = x momentum of particle n
   part[3][n] = y momentum of particle n
   part[4][n] = z momentum of particle n
   cu[k][j][i] = ith component of current density at grid point j,k
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
   nop = number of particles
   npe = first dimension of particle array
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   nxv = first dimension of current array, must be >= nx+1
   nyv = second dimension of current array, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   requires SSE2, part and cu need to be 16 byte aligned
   npe needs to be a multiple of 4, cu needs to have 4 components
local data                                                            */
   int j, nps, nn, mm;
   float ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami;
   __m128i v_nxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qm, v_dt, v_one, v_ci2;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_vx, v_vy, v_ux, v_uy;
   __m128 v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 a, b, c, d, va, vb, vc, vd;
   __attribute__((aligned(16))) unsigned int ll[4];
   __attribute__((aligned(16))) unsigned long kk[1];
   ci2 = ci*ci;
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
   v_qm = _mm_set1_ps(qm);
   v_one = _mm_set1_ps(1.0f);
   v_dt = _mm_set1_ps(dt);
   v_ci2 = _mm_set1_ps(ci2);
   v_edgelx = _mm_set1_ps(edgelx);
   v_edgely = _mm_set1_ps(edgely);
   v_edgerx = _mm_set1_ps(edgerx);
   v_edgery = _mm_set1_ps(edgery);
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
/*    dyp = y - (float) mm;      */
      v_dyp = _mm_sub_ps(v_y,_mm_cvtepi32_ps(v_mm));
/* find inverse gamma */
/*    ux = part[j+2*npe]; */
/*    uy = part[j+3*npe]; */
/*    uz = part[j+4*npe]; */
      v_ux = _mm_load_ps(&part[j+2*npe]);
      v_uy = _mm_load_ps(&part[j+3*npe]);
      vc = _mm_load_ps(&part[j+4*npe]);
/*    p2 = ux*ux + uy*uy + uz*uz; */
      v_at = _mm_mul_ps(v_ux,v_ux);
      v_at = _mm_add_ps(v_at,_mm_mul_ps(v_uy,v_uy));
      v_at = _mm_add_ps(v_at,_mm_mul_ps(vc,vc));
/*    gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*    v_at = _mm_rsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* full accuracy calculation */
      v_at = _mm_sqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2)));
      v_at = _mm_div_ps(v_one,v_at);
/* full accuracy calculation with SVML */
/*    v_at = _mm_invsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* calculate weights */
/*    nn = 4*(nn + nxv*mm); */
      v_it = _mm_mul_epu32(v_nxv,_mm_srli_si128(v_mm,4));
      v_mm = _mm_mul_epu32(v_mm,v_nxv);
      v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
      v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
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
/* deposit current */
/*    vx = ux*gami; */
/*    vy = uy*gami; */
/*    vz = uz*gami; */
      v_vx = _mm_mul_ps(v_ux,v_at);
      v_vy = _mm_mul_ps(v_uy,v_at);
      va = v_vx;
      vb = v_vy;
      vc = _mm_mul_ps(vc,v_at);
      vd = _mm_setzero_ps();
/* transpose so va,vb,vc,vd contain the 3 velocities plus zero */
/* for each of 4 particles                                     */
      _MM_TRANSPOSE4_PS(va,vb,vc,vd);
/*    dx = amx*amy; */
/*    cu[nn] += vx*dx;   */
/*    cu[nn+1] += vy*dx; */
/*    cu[nn+2] += vz*dx; */
/*    dy = dxp*amy;      */
/*    mm = nn + 4;       */
/*    cu[mm] += vx*dy;   */
/*    cu[mm+1] += vy*dy; */
/*    cu[mm+2] += vz*dy; */
/*    dx = amx*dyp;      */
/*    nn += 4*nxv;       */
/*    cu[nn] += vx*dx;   */
/*    cu[nn+1] += vy*dx; */
/*    cu[nn+2] += vz*dx; */
/*    dy = dxp*dyp;      */
/*    mm = nn + 4;       */
/*    cu[mm] += vx*dy;   */
/*    cu[mm+1] += vy*dy; */
/*    cu[mm+2] += vz*dy; */
/* deposit for first particle */
      mm = ll[0];
      v_dx = _mm_load_ps(&cu[mm]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(a,a,0)));
      _mm_store_ps(&cu[mm],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(b,b,0)));
      _mm_store_ps(&cu[mm+4],v_dy);
      v_dx = _mm_load_ps(&cu[mm+4*nxv]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(c,c,0)));
      _mm_store_ps(&cu[mm+4*nxv],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4+4*nxv]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(d,d,0)));
      _mm_store_ps(&cu[mm+4+4*nxv],v_dy);
/* deposit for second particle */
      mm = ll[1];
      v_dx = _mm_load_ps(&cu[mm]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(a,a,85)));
      _mm_store_ps(&cu[mm],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(b,b,85)));
      _mm_store_ps(&cu[mm+4],v_dy);
      v_dx = _mm_load_ps(&cu[mm+4*nxv]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(c,c,85)));
      _mm_store_ps(&cu[mm+4*nxv],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4+4*nxv]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(d,d,85)));
      _mm_store_ps(&cu[mm+4+4*nxv],v_dy);
/* deposit for third particle */
      mm = ll[2];
      v_dx = _mm_load_ps(&cu[mm]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(a,a,170)));
      _mm_store_ps(&cu[mm],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(b,b,170)));
      _mm_store_ps(&cu[mm+4],v_dy);
      v_dx = _mm_load_ps(&cu[mm+4*nxv]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(c,c,170)));
      _mm_store_ps(&cu[mm+4*nxv],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4+4*nxv]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(d,d,170)));
      _mm_store_ps(&cu[mm+4+4*nxv],v_dy);
/* deposit for fourth particle */
      mm = ll[3];
      v_dx = _mm_load_ps(&cu[mm]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(a,a,255)));
      _mm_store_ps(&cu[mm],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(b,b,255)));
      _mm_store_ps(&cu[mm+4],v_dy);
      v_dx = _mm_load_ps(&cu[mm+4*nxv]);
      v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(c,c,255)));
      _mm_store_ps(&cu[mm+4*nxv],v_dx);
      v_dy = _mm_load_ps(&cu[mm+4+4*nxv]);
      v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(d,d,255)));
      _mm_store_ps(&cu[mm+4+4*nxv],v_dy);
/* advance position half a time-step */
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
/*          part[j+2*npe] = -ux;                */
/*       }                                      */
         v_at = _mm_cmplt_ps(v_dx,v_edgelx);
         v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dx,v_edgerx));
         v_x = _mm_and_ps(v_at,v_x);
         v_dx = _mm_add_ps(_mm_andnot_ps(v_at,v_dx),v_x);
         v_dxp = _mm_and_ps(v_at,v_ux);
         v_ux = _mm_sub_ps(_mm_andnot_ps(v_at,v_ux),v_dxp);
/* write output if test result is true for any particle */
         v_mm = _mm_srli_si128((__m128i)v_at,8);
         v_mm = _mm_add_epi64((__m128i)v_at,v_mm);
         _mm_storel_epi64((__m128i *)&kk[0],v_mm);
         if (kk[0] != 0)
            _mm_store_ps(&part[j+2*npe],v_ux);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          part[j+3*npe] = -uy;                */
/*       }                                      */
         v_at = _mm_cmplt_ps(v_dy,v_edgely);
         v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dy,v_edgery));
         v_y = _mm_and_ps(v_at,v_y);
         v_dy = _mm_add_ps(_mm_andnot_ps(v_at,v_dy),v_y);
         v_dyp = _mm_and_ps(v_at,v_uy);
         v_uy = _mm_sub_ps(_mm_andnot_ps(v_at,v_uy),v_dyp);
/* write output if test result is true for any particle */
         v_mm = _mm_srli_si128((__m128i)v_at,8);
         v_mm = _mm_add_epi64((__m128i)v_at,v_mm);
         _mm_storel_epi64((__m128i *)&kk[0],v_mm);
         if (kk[0] != 0)
            _mm_store_ps(&part[j+3*npe],v_uy);
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          part[j+2*npe] = -ux;                */
/*       }                                      */
         v_at = _mm_cmplt_ps(v_dx,v_edgelx);
         v_at = _mm_or_ps(v_at,_mm_cmpge_ps(v_dx,v_edgerx));
         v_x = _mm_and_ps(v_at,v_x);
         v_dx = _mm_add_ps(_mm_andnot_ps(v_at,v_dx),v_x);
         v_dxp = _mm_and_ps(v_at,v_ux);
         v_ux = _mm_sub_ps(_mm_andnot_ps(v_at,v_ux),v_dxp);
/* write output if test result is true for any particle */
         v_mm = _mm_srli_si128((__m128i)v_at,8);
         v_mm = _mm_add_epi64((__m128i)v_at,v_mm);
         _mm_storel_epi64((__m128i *)&kk[0],v_mm);
         if (kk[0] != 0)
            _mm_store_ps(&part[j+2*npe],v_ux);
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
/* find inverse gamma */
      ux = part[j+2*npe];
      uy = part[j+3*npe];
      uz = part[j+4*npe];
      p2 = ux*ux + uy*uy + uz*uz;
      gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* calculate weights */
      nn = 4*(nn + nxv*mm);
      amx = qm - dxp;
      amy = 1.0f - dyp;
/* deposit current */
      dx = amx*amy;
      dy = dxp*amy;
      vx = ux*gami;
      vy = uy*gami;
      vz = uz*gami;
      cu[nn] += vx*dx;
      cu[nn+1] += vy*dx;
      cu[nn+2] += vz*dx;
      dx = amx*dyp;
      mm = nn + 4;
      cu[mm] += vx*dy;
      cu[mm+1] += vy*dy;
      cu[mm+2] += vz*dy;
      dy = dxp*dyp;
      nn += 4*nxv;
      cu[nn] += vx*dx;
      cu[nn+1] += vy*dx;
      cu[nn+2] += vz*dx;
      mm = nn + 4;
      cu[mm] += vx*dy;
      cu[mm+1] += vy*dy;
      cu[mm+2] += vz*dy;
/* advance position half a time-step */
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
            part[j+2*npe] = -ux;
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = y;
            part[j+3*npe] = -uy;
         }
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            part[j+2*npe] = -ux;
         }
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* set new position */
      part[j] = dx;
      part[j+npe] = dy;
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
void csse2bguard2l(float bxy[], int nx, int ny, int nxe, int nye) {
/* replicate extended periodic vector field bxy
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   requires SSE2, bxy needs to be 16 byte aligned
   bxy needs to have 4 components
local data                                                 */
   int j, k, kk;
/* copy edges of extended field */
   for (k = 0; k < ny; k++) {
      kk = 4*nxe*k;
      bxy[4*nx+kk] = bxy[kk];
      bxy[1+4*nx+kk] = bxy[1+kk];
      bxy[2+4*nx+kk] = bxy[2+kk];
   }
   kk = 4*nxe*ny;
   for (j = 0; j < nx; j++) {
/*    bxy[4*j+kk] = bxy[4*j];     */
/*    bxy[1+4*j+kk] = bxy[1+4*j]; */
/*    bxy[2+4*j+kk] = bxy[2+4*j]; */
      _mm_store_ps(&bxy[4*j+kk],_mm_load_ps(&bxy[4*j]));
   }
   bxy[4*nx+kk] = bxy[0];
   bxy[1+4*nx+kk] = bxy[1];
   bxy[2+4*nx+kk] = bxy[2];
   return;
}


/*--------------------------------------------------------------------*/
void csse2acguard2l(float cu[], int nx, int ny, int nxe, int nye) {
/* accumulate extended periodic vector field cu
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   requires SSE2, cu needs to be 16 byte aligned
   cu needs to have 4 components
local data                                                 */
   int j, k, kk;
   __m128 v_cu, v_zero;
   v_zero = _mm_set1_ps(0.0f);
/* accumulate edges of extended field */
   for (k = 0; k < ny; k++) {
      kk = 4*nxe*k;
      cu[kk] += cu[4*nx+kk];
      cu[1+kk] += cu[1+4*nx+kk];
      cu[2+kk] += cu[2+4*nx+kk];
      cu[4*nx+kk] = 0.0;
      cu[1+4*nx+kk] = 0.0;
      cu[2+4*nx+kk] = 0.0;
   }
   kk = 4*nxe*ny;
   for (j = 0; j < nx; j++) {
/*    cu[4*j] += cu[4*j+kk];     */
/*    cu[1+4*j] += cu[1+4*j+kk]; */
/*    cu[2+4*j] += cu[2+4*j+kk]; */
/*    cu[4*j+kk] = 0.0;          */
/*    cu[1+4*j+kk] = 0.0;        */
/*    cu[2+4*j+kk] = 0.0;        */
      v_cu = _mm_add_ps(_mm_load_ps(&cu[4*j]),_mm_load_ps(&cu[4*j+kk]));
      _mm_store_ps(&cu[4*j],v_cu);
      _mm_store_ps(&cu[4*j+kk],v_zero);
   }
   cu[0] += cu[4*nx+kk];
   cu[1] += cu[1+4*nx+kk];
   cu[2] += cu[2+4*nx+kk];
   cu[4*nx+kk] = 0.0;
   cu[1+4*nx+kk] = 0.0;
   cu[2+4*nx+kk] = 0.0;
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
void csse2pois23(float complex q[], float complex fxy[], int isign,
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
   requires SSE2, q, fxy, ffc need to be 16 byte aligned
   nxhd, nxvh need to be a multiple of 2
   fxy needs to have 4 components
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
/*       fxy[4*j+4*kj] = at2*zt1;   */
/*       fxy[1+4*j+4*kj] = at3*zt1; */
/*       fxy[2+4*j+4*kj] = zero;    */
         v_at4 = _mm_mul_ps(v_at2,v_zt1);
         v_zt4 = _mm_mul_ps(v_at3,v_zt1);
/* reorder write */
         v_zt3 = _mm_shuffle_ps(v_at4,v_zt4,68);
         v_zt4 = _mm_shuffle_ps(v_at4,v_zt4,238);
         _mm_store_ps((float *)&fxy[4*(j+kj)],v_zt3);
         _mm_store_ps((float *)&fxy[2+4*(j+kj)],v_zero);
         _mm_store_ps((float *)&fxy[4*(j+1+kj)],v_zt4);
         _mm_store_ps((float *)&fxy[2+4*(j+1+kj)],v_zero);
/*       fxy[4*j+4*k1] = at2*zt2;    */
/*       fxy[1+4*j+4*k1] = -at3*zt2; */
/*       fxy[2+4*j+4*k1] = zero;     */
         v_at4 = _mm_mul_ps(v_at2,v_zt2);
         v_zt4 = _mm_sub_ps(v_zero,_mm_mul_ps(v_at3,v_zt2));
/* reorder write */
         v_zt3 = _mm_shuffle_ps(v_at4,v_zt4,68);
         v_zt4 = _mm_shuffle_ps(v_at4,v_zt4,238);
         _mm_store_ps((float *)&fxy[4*(j+k1)],v_zt3);
         _mm_store_ps((float *)&fxy[2+4*(j+k1)],v_zero);
         _mm_store_ps((float *)&fxy[4*(j+1+k1)],v_zt4);
         _mm_store_ps((float *)&fxy[2+4*(j+1+k1)],v_zero);
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
         fxy[4*j+4*kj] = at2*zt1;
         fxy[1+4*j+4*kj] = at3*zt1;
         fxy[2+4*j+4*kj] = zero;
         fxy[4*j+4*k1] = at2*zt2;
         fxy[1+4*j+4*k1] = -at3*zt2;
         fxy[2+4*j+4*k1] = zero;
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
      fxy[4*kj] = zero;
      fxy[1+4*kj] = at3*zt1;
      fxy[2+4*kj] = zero;
      fxy[4*k1] = zero;
      fxy[1+4*k1] = zero;
      fxy[2+4*k1] = zero;
      wp += at1*(q[kj]*conjf(q[kj]));
   }
/* mode numbers ky = 0, ny/2 */
   k1 = 4*nxvh*nyh;
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
/*    fxy[4*j] = at2*zt1; */
/*    fxy[1+4*j] = zero;  */
/*    fxy[2+4*j] = zero;  */
      v_at4 = _mm_mul_ps(v_at2,v_zt1);
/* reorder write */
      v_zt3 = _mm_shuffle_ps(v_at4,v_zero,68);
      v_zt4 = _mm_shuffle_ps(v_at4,v_zero,238);
      _mm_store_ps((float *)&fxy[4*j],v_zt3);
      _mm_store_ps((float *)&fxy[4*j+2],v_zero);
      _mm_store_ps((float *)&fxy[4*j+4],v_zt4);
      _mm_store_ps((float *)&fxy[4*j+6],v_zero);
/*    fxy[4*j+k1] = zero;   */
/*    fxy[1+4*j+k1] = zero; */
/*    fxy[2+4*j+k1] = zero; */
      _mm_store_ps((float *)&fxy[4*j+k1],v_zero);
      _mm_store_ps((float *)&fxy[4*j+2+k1],v_zero);
      _mm_store_ps((float *)&fxy[4*j+4+k1],v_zero);
      _mm_store_ps((float *)&fxy[4*j+6+k1],v_zero);
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
      fxy[4*j] = at2*zt1;
      fxy[1+4*j] = zero;
      fxy[2+4*j] = zero;
      fxy[4*j+k1] = zero;
      fxy[1+4*j+k1] = zero;
      fxy[2+4*j+k1] = zero;
      wp += at1*(q[j]*conjf(q[j]));
   }
   fxy[0] = zero;
   fxy[1] = zero;
   fxy[2] = zero;
   fxy[k1] = zero;
   fxy[1+k1] = zero;
   fxy[2+k1] = zero;
/* *we = wp*(float) (nx*ny); */
   _mm_store_pd(&dd[0],v_wp);
   for (j = 1; j < 2; j++) {
      dd[0] += dd[j];
   }
   *we = (wp + dd[0])*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void csse2cuperp2(float complex cu[], int nx, int ny, int nxvh,
                  int nyv) {
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
   requires SSE2, cu needs to be 16 byte aligned
   cu needs to have 4 components
local data                                                 */
   int nxh, nyh, j, k, k1, kj;
   float dnx, dny;
   float complex zero;
   __m128 v_dnx, v_dny, v_dkx, v_dky, v_dky2, v_at1, v_zt1, v_zt2;
   __m128 v_zero, v_one, v_n, v_at;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
   v_dnx = _mm_set1_ps(dnx);
   v_dny = _mm_set1_ps(dny);
   v_zero = _mm_set1_ps(0.0f);
   v_one = _mm_set1_ps(1.0f);
   v_n = _mm_set_ps(-1.0f,-1.0f,1.0f,1.0f);
/* calculate transverse part of current */
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
/*    dky = dny*(float) k; */
      v_dky = _mm_mul_ps(v_dny,_mm_cvtepi32_ps(_mm_set1_epi32(k)));
/*    dky2 = dky*dky; */
      v_dky2 = _mm_mul_ps(v_dky,v_dky);
      kj = 4*nxvh*k;
      k1 = 4*nxvh*ny - kj;
      for (j = 1; j < nxh; j++) {
/*       dkx = dnx*(float) j; */
         v_dkx = _mm_mul_ps(v_dnx,_mm_cvtepi32_ps(_mm_set1_epi32(j)));
/*       at1 = 1./(dkx*dkx + dky2); */
         v_at1 = _mm_add_ps(_mm_mul_ps(v_dkx,v_dkx),v_dky2);
         v_at1 = _mm_div_ps(v_one,v_at1);
/*       zt1 = at1*(dkx*cu[4*j+kj] + dky*cu[1+4*j+kj]); */
         v_dkx = _mm_movelh_ps(v_dkx,v_dky);
         v_zt2 = _mm_load_ps((float *)&cu[4*j+kj]);
         v_zt1 = _mm_mul_ps(v_at1,_mm_mul_ps(v_dkx,v_zt2));
         v_at = _mm_movelh_ps(v_zt1,v_zt1);
         v_zt1 = _mm_add_ps(_mm_movehl_ps(v_zt1,v_zt1),v_at);
/*       cu[4*j+kj] -= dkx*zt1;   */
/*       cu[1+4*j+kj] -= dky*zt1; */
         v_zt2 = _mm_sub_ps(v_zt2,_mm_mul_ps(v_dkx,v_zt1));
         _mm_store_ps((float *)&cu[4*j+kj],v_zt2);
/*       zt1 = at1*(dkx*cu[4*j+k1] - dky*cu[1+4*j+k1]); */
         v_dkx = _mm_mul_ps(v_dkx,v_n);
         v_zt2 = _mm_load_ps((float *)&cu[4*j+k1]);
         v_zt1 = _mm_mul_ps(v_at1,_mm_mul_ps(v_dkx,v_zt2));
         v_at = _mm_movelh_ps(v_zt1,v_zt1);
         v_zt1 = _mm_add_ps(_mm_movehl_ps(v_zt1,v_zt1),v_at);
/*       cu[4*j+k1] -= dkx*zt1;   */
/*       cu[1+4*j+k1] += dky*zt1; */
         v_zt2 = _mm_sub_ps(v_zt2,_mm_mul_ps(v_dkx,v_zt1));
         _mm_store_ps((float *)&cu[4*j+k1],v_zt2);
      }
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kj = 4*nxvh*k;
      k1 = 4*nxvh*ny - kj;
      cu[1+kj] = zero;
      cu[k1] = zero;
      cu[1+k1] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = 4*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
      cu[4*j] = zero;
/*    cu[4*j+k1] = zero;   */
/*    cu[1+4*j+k1] = zero; */
      _mm_store_ps((float *)&cu[4*j+k1],v_zero);
   }
/* cu[0] = zero; */
/* cu[1] = zero; */
   _mm_store_ps((float *)&cu[0],v_zero);
/* cu[k1] = zero; */
/* cu[1+k1] = zero; */
   _mm_store_ps((float *)&cu[k1],v_zero);
   return;
}

/*--------------------------------------------------------------------*/
void csse2ibpois23(float complex cu[], float complex bxy[],
                   float complex ffc[], float ci, float *wm, int nx,
                   int ny, int nxvh, int nyv, int nxhd, int nyhd) {
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
   requires SSE2, cu, bxy, ffc need to be 16 byte aligned
   nxhd, nxvh need to be a multiple of 2
   cu, bxy need to have 4 components
local data                                                 */
   int nxh, nyh, j, k, k1, kk, kj;
   float dnx, dny, ci2, at1, at3;
   float complex zero, zt1, zt3;
   double wp;
   __m128i v_it;
   __m128 v_dnx, v_dny, v_dky, v_ci2, v_at1, v_at2, v_at3, v_at4;
   __m128 v_zero, v_n, v_m, v_zt1, v_zt2, v_zt3, v_zt4;
   __m128d v_wp, v_d;
   __attribute__((aligned(16))) double dd[2];
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   zero = 0.0 + 0.0*_Complex_I;
   ci2 = ci*ci;
   v_dnx = _mm_set1_ps(dnx);
   v_dny = _mm_set1_ps(dny);
   v_zero = _mm_set1_ps(0.0f);
   v_ci2 = _mm_set1_ps(ci2);
   v_n = _mm_set_ps(-1.0f,1.0f,-1.0f,1.0f);
   v_m = _mm_set_ps(1.0f,1.0f,-1.0f,-1.0f);
/* calculate magnetic field and sum field energy */
   wp = 0.0;
   v_wp = _mm_set1_pd(0.0);
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
/*    dky = dny*(float) k; */
      v_dky = _mm_mul_ps(v_dny,_mm_cvtepi32_ps(_mm_set1_epi32(k)));
      kk = nxhd*k;
      kj = 4*nxvh*k;
      k1 = 4*nxvh*ny - kj;
      for (j = 1; j < nxh; j++) {
/*       at1 = ci2*crealf(ffc[j+kk]); */
         v_at3 = _mm_loadl_pi(v_zero,(__m64 *)&ffc[j+kk]);
         v_at1 = _mm_mul_ps(v_ci2,_mm_shuffle_ps(v_at3,v_at3,0));
/*       at2 = at1*dnx*(float) j; */
/*       at3 = dky*at1;           */
         v_at2 = _mm_mul_ps(v_dnx,_mm_cvtepi32_ps(_mm_set1_epi32(j)));
         v_at2 = _mm_sub_ps(v_zero,v_at2);
         v_at2 = _mm_mul_ps(v_at1,_mm_movelh_ps(v_dky,v_at2));
/*       at1 = at1*cimagf(ffc[j+kk]); */
         v_at3 = _mm_movelh_ps(v_at3,v_at3);
         v_at1 = _mm_mul_ps(v_at1,_mm_shuffle_ps(v_at3,v_at3,245));
/*       zt1 = -cimagf(cu[2+4*j+kj])                                */
/*             + crealf(cu[2+4*j+kj])*_Complex_I;                   */
/*       zt2 = -cimagf(cu[1+4*j+kj])                                */
/*           + crealf(cu[1+4*j+kj])*_Complex_I;                     */
/*       zt3 = -cimagf(cu[4*j+kj]) + crealf(cu[4*j+kj])*_Complex_I; */
         v_zt3 = _mm_load_ps((float *)&cu[4*j+kj]);
         v_zt2 = _mm_mul_ps(v_zt3,v_n);
         v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
         v_zt4 = _mm_load_ps((float *)&cu[2+4*j+kj]);
         v_zt1 = _mm_mul_ps(v_zt4,v_n);
         v_zt1 = _mm_movelh_ps(v_zt1,v_zt1);
         v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*       bxy[4*j+kj] = at3*zt1;    */
/*       bxy[1+4*j+kj] = -at2*zt1; */
         v_zt1 = _mm_mul_ps(v_at2,v_zt1);
         _mm_store_ps((float *)&bxy[4*j+kj],v_zt1);
/*       bxy[2+4*j+kj] = at2*zt2 - at3*zt3; */
         v_zt2 = _mm_sub_ps(v_zero,_mm_mul_ps(v_at2,v_zt2));
         v_at3 = _mm_movelh_ps(v_zt2,v_zt2);
         v_zt2 = _mm_add_ps(_mm_movehl_ps(v_zt2,v_zt2),v_at3);
         v_zt2 = _mm_movelh_ps(v_zt2,v_zero);
         _mm_store_ps((float *)&bxy[2+4*j+kj],v_zt2);
/*       wp += at1*(cu[4*j+kj]*conjf(cu[4*j+kj])    */
/*             + cu[1+4*j+kj]*conjf(cu[1+4*j+kj])   */
/*             + cu[2+4*j+kj]*conjf(cu[2+4*j+kj])); */
         v_zt3 = _mm_mul_ps(v_zt3,v_zt3);
         v_zt4 = _mm_movelh_ps(v_zt4,v_zero);
         v_at4 = _mm_add_ps(v_zt3,_mm_mul_ps(v_zt4,v_zt4));
/*       zt1 = -cimagf(cu[2+4*j+k1])                                */
/*             + crealf(cu[2+4*j+k1])*_Complex_I;                   */
/*       zt2 = -cimagf(cu[1+4*j+k1])                                */
/*           + crealf(cu[1+4*j+k1])*_Complex_I;                     */
/*       zt3 = -cimagf(cu[4*j+k1]) + crealf(cu[4*j+k1])*_Complex_I; */
         v_zt3 = _mm_load_ps((float *)&cu[4*j+k1]);
         v_zt2 = _mm_mul_ps(v_zt3,v_n);
         v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
         v_zt4 = _mm_load_ps((float *)&cu[2+4*j+k1]);
         v_zt1 = _mm_movelh_ps(v_zt4,v_zt4);
         v_zt1 = _mm_mul_ps(v_zt1,v_n);
         v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*       bxy[4*j+k1] = -at3*zt1;   */
/*       bxy[1+4*j+k1] = -at2*zt1; */
         v_at2 = _mm_mul_ps(v_at2,v_m);
         v_zt1 = _mm_mul_ps(v_at2,v_zt1);
         _mm_store_ps((float *)&bxy[4*j+k1],v_zt1);
/*       bxy[2+4*j+k1] = at2*zt2 + at3*zt3; */
         v_zt2 = _mm_sub_ps(v_zero,_mm_mul_ps(v_at2,v_zt2));
         v_at3 = _mm_movelh_ps(v_zt2,v_zt2);
         v_zt2 = _mm_add_ps(_mm_movehl_ps(v_zt2,v_zt2),v_at3);
         v_zt2 = _mm_movelh_ps(v_zt2,v_zero);
         _mm_store_ps((float *)&bxy[2+4*j+k1],v_zt2);
/*       wp += at1*(cu[4*j+k1]*conjf(cu[4*j+k1])    */
/*             + cu[1+4*j+k1]*conjf(cu[1+4*j+k1])   */
/*             + cu[2+4*j+k1]*conjf(cu[2+4*j+k1])); */
         v_zt3 = _mm_mul_ps(v_zt3,v_zt3);
         v_zt4 = _mm_movelh_ps(v_zt4,v_zero);
         v_zt3 = _mm_add_ps(v_zt3,_mm_mul_ps(v_zt4,v_zt4));
         v_at4 = _mm_mul_ps(v_at1,_mm_add_ps(v_at4,v_zt3));
/* convert to double precision before accumulating */
         v_d = _mm_cvtps_pd(v_at4);
         v_wp = _mm_add_pd(v_wp,v_d);
         v_it = _mm_srli_si128((__m128i)v_at4,8);
         v_d = _mm_cvtps_pd((__m128)v_it);
         v_wp = _mm_add_pd(v_wp,v_d);
      }
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = 4*nxvh*k;
      k1 = 4*nxvh*ny - kj;
      at1 = ci2*crealf(ffc[kk]);
      at3 = at1*dny*(float) k;
      at1 = at1*cimagf(ffc[kk]);
      zt1 = -cimagf(cu[2+kj]) + crealf(cu[2+kj])*_Complex_I;
      zt3 = -cimagf(cu[kj]) + crealf(cu[kj])*_Complex_I;
      bxy[kj] = at3*zt1;
      bxy[1+kj] = zero;
      bxy[2+kj] = -at3*zt3;
      bxy[3+kj] = zero;
      bxy[k1] = zero;
      bxy[1+k1] = zero;
      bxy[2+k1] = zero;
      bxy[3+k1] = zero;
      wp += at1*(cu[kj]*conjf(cu[kj]) + cu[1+kj]*conjf(cu[1+kj])
            + cu[2+kj]*conjf(cu[2+kj]));
   }
/* mode numbers ky = 0, ny/2 */
   k1 = 4*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
/*    at1 = ci2*crealf(ffc[j]); */
      v_at3 = _mm_loadl_pi(v_zero,(__m64 *)&ffc[j]);
      v_at1 = _mm_mul_ps(v_ci2,_mm_shuffle_ps(v_at3,v_at3,0));
/*    at2 = at1*dnx*(float) j; */
      v_at2 = _mm_mul_ps(v_dnx,_mm_cvtepi32_ps(_mm_set1_epi32(j)));
      v_at2 = _mm_sub_ps(v_zero,v_at2);
      v_at2 = _mm_mul_ps(v_at1,_mm_movelh_ps(v_zero,v_at2));
/*    at1 = at1*cimagf(ffc[j]); */
      v_at3 = _mm_movelh_ps(v_at3,v_at3);
      v_at1 = _mm_mul_ps(v_at1,_mm_shuffle_ps(v_at3,v_at3,245));
/*    zt1 = -cimagf(cu[2+4*j]) + crealf(cu[2+4*j])*_Complex_I; */
/*    zt2 = -cimagf(cu[1+4*j]) + crealf(cu[1+4*j])*_Complex_I; */
      v_zt3 = _mm_load_ps((float *)&cu[4*j]);
      v_zt2 = _mm_mul_ps(v_zt3,v_n);
      v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
      v_zt4 = _mm_load_ps((float *)&cu[2+4*j]);
      v_zt1 = _mm_mul_ps(v_zt4,v_n);
      v_zt1 = _mm_movelh_ps(v_zt1,v_zt1);
      v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*    bxy[4*j] = zero;       */
/*    bxy[1+4*j] = -at2*zt1; */
      v_zt1 = _mm_mul_ps(v_at2,v_zt1);
      _mm_store_ps((float *)&bxy[4*j],v_zt1);
/*    bxy[2+4*j] = at2*zt2; */
      v_zt2 = _mm_sub_ps(v_zero,_mm_mul_ps(v_at2,v_zt2));
      v_zt2 = _mm_movehl_ps(v_zt2,v_zt2);
      v_zt2 = _mm_movelh_ps(v_zt2,v_zero);
      _mm_store_ps((float *)&bxy[2+4*j],v_zt2);
/*    bxy[4*j+k1] = zero;   */
/*    bxy[1+4*j+k1] = zero; */
/*    bxy[2+4*j+k1] = zero; */
      _mm_store_ps((float *)&bxy[4*j+k1],v_zero);
      _mm_store_ps((float *)&bxy[2+4*j+k1],v_zero);
/*    wp += at1*(cu[4*j]*conjf(cu[4*j]) + cu[1+4*j]*conjf(cu[1+4*j]) */
/*          + cu[2+4*j]*conjf(cu[2+4*j]));                           */
      v_zt3 = _mm_mul_ps(v_zt3,v_zt3);
      v_zt4 = _mm_movelh_ps(v_zt4,v_zero);
      v_at4 = _mm_add_ps(v_zt3,_mm_mul_ps(v_zt4,v_zt4));
      v_at4 = _mm_mul_ps(v_at1,v_at4);
/* convert to double precision before accumulating */
      v_d = _mm_cvtps_pd(v_at4);
      v_wp = _mm_add_pd(v_wp,v_d);
      v_it = _mm_srli_si128((__m128i)v_at4,8);
      v_d = _mm_cvtps_pd((__m128)v_it);
      v_wp = _mm_add_pd(v_wp,v_d);
   }
/* bxy[0] = zero; */
/* bxy[1] = zero; */
/* bxy[2] = zero; */
   _mm_store_ps((float *)&bxy[0],v_zero);
   _mm_store_ps((float *)&bxy[2],v_zero);
/* bxy[k1] = zero;   */
/* bxy[1+k1] = zero; */
/* bxy[2+k1] = zero; */
   _mm_store_ps((float *)&bxy[k1],v_zero);
    _mm_store_ps((float *)&bxy[2+k1],v_zero);
/* *wm = wp*(float) (nx*ny); */
   _mm_store_pd(&dd[0],v_wp);
   for (j = 1; j < 2; j++) {
      dd[0] += dd[j];
   }
   *wm = (wp + dd[0])*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void csse2maxwel2(float complex exy[], float complex bxy[],
                  float complex cu[], float complex ffc[], float ci,
                  float dt, float *wf, float *wm, int nx, int ny,
                  int nxvh, int nyv, int nxhd, int nyhd) {
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
   requires SSE2, cu, bxy, ffc need to be 16 byte aligned
   nxhd, nxvh need to be a multiple of 2
   cu, exy, bxy need to have 4 components
local data                                                 */
   int nxh, nyh, j, k, k1, kk, kj;
   float dnx, dny, dth, c2, cdt, affp, anorm, dky, afdt, adt;
   float complex zero, zt1, zt3, zt4, zt6, zt7, zt9;
   double wp, ws;
   __m128i v_it;
   __m128 v_dnx, v_dkx, v_dny, v_dky, v_cdt, v_adt, v_afdt, v_dth;
   __m128 v_anorm, v_at1, v_at2, v_at3, v_at4;
   __m128 v_zero, v_n, v_m, v_zt1, v_zt2, v_zt4, v_zt6, v_zt7, v_zt9;
   __m128d v_wp, v_ws, v_d;
   __attribute__((aligned(16))) double dd[2];
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
   v_dnx = _mm_set1_ps(dnx);
   v_dny = _mm_set1_ps(dny);
   v_zero = _mm_set1_ps(0.0f);
   v_cdt = _mm_set1_ps(cdt);
   v_adt = _mm_set1_ps(adt);
   v_dth = _mm_set1_ps(dth);
   v_anorm = _mm_set1_ps(anorm);
   v_n = _mm_set_ps(-1.0f,1.0f,-1.0f,1.0f);
   v_m = _mm_set_ps(1.0f,1.0f,-1.0f,-1.0f);
/* update electromagnetic field and sum field energies */
   ws = 0.0;
   wp = 0.0;
   v_wp = _mm_set1_pd(0.0);
   v_ws = _mm_set1_pd(0.0);
/* calculate the electromagnetic fields */
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
   for (k = 1; k < nyh; k++) {
/*    dky = dny*(float) k; */
      v_dky = _mm_mul_ps(v_dny,_mm_cvtepi32_ps(_mm_set1_epi32(k)));
      kk = nxhd*k;
      kj = 4*nxvh*k;
      k1 = 4*nxvh*ny - kj;
      for (j = 1; j < nxh; j++) {
/*       dkx = dnx*(float) j; */
         v_dkx = _mm_mul_ps(v_dnx,_mm_cvtepi32_ps(_mm_set1_epi32(j)));
         v_dkx = _mm_movelh_ps(_mm_sub_ps(v_zero,v_dky),v_dkx);
/*       afdt = adt*cimagf(ffc[j+kk]); */
         v_afdt = _mm_loadl_pi(v_zero,(__m64 *)&ffc[j+kk]);
         v_afdt = _mm_movelh_ps(v_afdt,v_afdt);
         v_afdt = _mm_mul_ps(v_adt,_mm_shuffle_ps(v_afdt,v_afdt,245));
/* update magnetic field half time step, ky > 0 */
/*       zt1 = -cimagf(exy[2+4*j+kj])                                 */
/*             + crealf(exy[2+4*j+kj])*_Complex_I;                    */
/*       zt2 = -cimagf(exy[1+4*j+kj])                                 */
/*           + crealf(exy[1+4*j+kj])*_Complex_I;                      */
/*       zt3 = -cimagf(exy[4*j+kj]) + crealf(exy[4*j+kj])*_Complex_I; */
         v_zt7 = _mm_load_ps((float *)&exy[4*j+kj]);
         v_zt2 = _mm_mul_ps(v_zt7,v_n);
         v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
         v_zt9 = _mm_load_ps((float *)&exy[2+4*j+kj]);
         v_zt1 = _mm_mul_ps(v_zt9,v_n);
         v_zt1 = _mm_movelh_ps(v_zt1,v_zt1);
         v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*       zt4 = bxy[4*j+kj] - dth*(dky*zt1);   */
/*       zt5 = bxy[1+4*j+kj] + dth*(dkx*zt1); */
         v_zt4 = _mm_load_ps((float *)&bxy[4*j+kj]);
         v_zt1 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt1));
         v_zt4 = _mm_add_ps(v_zt4,v_zt1);
/*       zt6 = bxy[2+4*j+kj] - dth*(dkx*zt2 - dky*zt3); */
         v_zt6 = _mm_load_ps((float *)&bxy[2+4*j+kj]);
         v_zt2 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt2));
         v_at1 = _mm_movelh_ps(v_zt2,v_zt2);
         v_zt2 = _mm_add_ps(_mm_movehl_ps(v_zt2,v_zt2),v_at1);
         v_zt6 = _mm_movelh_ps(_mm_sub_ps(v_zt6,v_zt2),v_zero);
/* update electric field whole time step */
/*       zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*       zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*       zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
         v_zt2 = _mm_mul_ps(v_zt4,v_n);
         v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
         v_zt1 = _mm_mul_ps(v_zt6,v_n);
         v_zt1 = _mm_movelh_ps(v_zt1,v_zt1);
         v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*       zt7 = exy[4*j+kj] + cdt*(dky*zt1) - afdt*cu[4*j+kj];     */
/*       zt8 = exy[1+4*j+kj] - cdt*(dkx*zt1) - afdt*cu[1+4*j+kj]; */
         v_at2 = _mm_load_ps((float *)&cu[4*j+kj]);
         v_zt1 = _mm_mul_ps(v_cdt,_mm_mul_ps(v_dkx,v_zt1));
         v_zt1 = _mm_add_ps(_mm_mul_ps(v_afdt,v_at2),v_zt1);
         v_zt7 = _mm_sub_ps(v_zt7,v_zt1);
/*       zt9 = exy[2+4*j+kj] + cdt*(dkx*zt2 - dky*zt3) */
/*             - afdt*cu[2+4*j+kj];                    */
         v_at2 = _mm_load_ps((float *)&cu[2+4*j+kj]);
         v_zt2 = _mm_mul_ps(v_cdt,_mm_mul_ps(v_dkx,v_zt2));
         v_at1 = _mm_movelh_ps(v_zt2,v_zt2);
         v_zt2 = _mm_add_ps(_mm_movehl_ps(v_zt2,v_zt2),v_at1);
         v_zt2 = _mm_sub_ps(v_zt2,_mm_mul_ps(v_afdt,v_at2));
         v_zt9 = _mm_movelh_ps(_mm_add_ps(v_zt9,v_zt2),v_zero);
/* update magnetic field half time step and store electric field */
/*       zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*       zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*       zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
         v_zt2 = _mm_mul_ps(v_zt7,v_n);
         v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
         v_zt1 = _mm_mul_ps(v_zt9,v_n);
         v_zt1 = _mm_movelh_ps(v_zt1,v_zt1);
         v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*       exy[4*j+kj] = zt7;   */
/*       exy[1+4*j+kj] = zt8; */
/*       exy[2+4*j+kj] = zt9; */
         _mm_store_ps((float *)&exy[4*j+kj],v_zt7);
         _mm_store_ps((float *)&exy[2+4*j+kj],v_zt9);
/*       ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
         v_zt7 = _mm_mul_ps(v_zt7,v_zt7);
         v_at3 = _mm_add_ps(v_zt7,_mm_mul_ps(v_zt9,v_zt9));
/*       zt4 -= dth*(dky*zt1); */
/*       zt5 += dth*(dkx*zt1); */
         v_zt1 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt1));
         v_zt4 = _mm_add_ps(v_zt4,v_zt1);
/*       zt6 -= dth*(dkx*zt2 - dky*zt3); */
         v_zt2 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt2));
         v_at1 = _mm_movelh_ps(v_zt2,v_zt2);
         v_zt2 = _mm_add_ps(_mm_movehl_ps(v_zt2,v_zt2),v_at1);
         v_zt6 = _mm_movelh_ps(_mm_sub_ps(v_zt6,v_zt2),v_zero);
/*       bxy[4*j+kj] = zt4;   */
/*       bxy[1+4*j+kj] = zt5; */
/*       bxy[2+4*j+kj] = zt6; */
         _mm_store_ps((float *)&bxy[4*j+kj],v_zt4);
         _mm_store_ps((float *)&bxy[2+4*j+kj],v_zt6);
/*       wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
         v_zt4 = _mm_mul_ps(v_zt4,v_zt4);
         v_at4 = _mm_add_ps(v_zt4,_mm_mul_ps(v_zt6,v_zt6));
/* update magnetic field half time step, ky < 0 */
         v_dkx = _mm_mul_ps(v_dkx,v_m);
/*       zt1 = -cimagf(exy[2+4*j+k1])                                 */
/*             + crealf(exy[2+4*j+k1])*_Complex_I;                    */
/*       zt2 = -cimagf(exy[1+4*j+k1])                                 */
/*             + crealf(exy[1+4*j+k1])*_Complex_I;                    */
/*       zt3 = -cimagf(exy[4*j+k1]) + crealf(exy[4*j+k1])*_Complex_I; */
         v_zt7 = _mm_load_ps((float *)&exy[4*j+k1]);
         v_zt2 = _mm_mul_ps(v_zt7,v_n);
         v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
         v_zt9 = _mm_load_ps((float *)&exy[2+4*j+k1]);
         v_zt1 = _mm_mul_ps(v_zt9,v_n);
         v_zt1 = _mm_movelh_ps(v_zt1,v_zt1);
         v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*       zt4 = bxy[4*j+k1] + dth*(dky*zt1);  */
/*       zt5 = bxy[1+4*j+k1] + dth*(dkx*zt1); */
         v_zt4 = _mm_load_ps((float *)&bxy[4*j+k1]);
         v_zt1 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt1));
         v_zt4 = _mm_add_ps(v_zt4,v_zt1);
/*       zt6 = bxy[2+4*j+k1] - dth*(dkx*zt2 + dky*zt3); */
         v_zt6 = _mm_load_ps((float *)&bxy[2+4*j+k1]);
         v_zt2 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt2));
         v_at1 = _mm_movelh_ps(v_zt2,v_zt2);
         v_zt2 = _mm_add_ps(_mm_movehl_ps(v_zt2,v_zt2),v_at1);
         v_zt6 = _mm_movelh_ps(_mm_sub_ps(v_zt6,v_zt2),v_zero);
/* update electric field whole time step */
/*       zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*       zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*       zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
         v_zt2 = _mm_mul_ps(v_zt4,v_n);
         v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
         v_zt1 = _mm_mul_ps(v_zt6,v_n);
         v_zt1 = _mm_movelh_ps(v_zt1,v_zt1);
         v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*       zt7 = exy[4*j+k1] - cdt*(dky*zt1) - afdt*cu[4*j+k1];     */
/*       zt8 = exy[1+4*j+k1] - cdt*(dkx*zt1) - afdt*cu[1+4*j+k1]; */
         v_at2 = _mm_load_ps((float *)&cu[4*j+k1]);
         v_zt1 = _mm_mul_ps(v_cdt,_mm_mul_ps(v_dkx,v_zt1));
         v_zt1 = _mm_add_ps(_mm_mul_ps(v_afdt,v_at2),v_zt1);
         v_zt7 = _mm_sub_ps(v_zt7,v_zt1);
/*       zt9 = exy[2+4*j+k1] + cdt*(dkx*zt2 + dky*zt3) */
/*             - afdt*cu[2+4*j+k1];                    */
         v_at2 = _mm_load_ps((float *)&cu[2+4*j+k1]);
         v_zt2 = _mm_mul_ps(v_cdt,_mm_mul_ps(v_dkx,v_zt2));
         v_at1 = _mm_movelh_ps(v_zt2,v_zt2);
         v_zt2 = _mm_add_ps(_mm_movehl_ps(v_zt2,v_zt2),v_at1);
         v_zt2 = _mm_sub_ps(v_zt2,_mm_mul_ps(v_afdt,v_at2));
         v_zt9 = _mm_movelh_ps(_mm_add_ps(v_zt9,v_zt2),v_zero);
/* update magnetic field half time step and store electric field */
/*       zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*       zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*       zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
         v_zt2 = _mm_mul_ps(v_zt7,v_n);
         v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
         v_zt1 = _mm_mul_ps(v_zt9,v_n);
         v_zt1 = _mm_movelh_ps(v_zt1,v_zt1);
         v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*       exy[4*j+k1] = zt7;   */
/*       exy[1+4*j+k1] = zt8; */
/*       exy[2+4*j+k1] = zt9; */
         _mm_store_ps((float *)&exy[4*j+k1],v_zt7);
         _mm_store_ps((float *)&exy[2+4*j+k1],v_zt9);
/*       ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
         v_zt7 = _mm_mul_ps(v_zt7,v_zt7);
         v_zt7 = _mm_add_ps(v_zt7,_mm_mul_ps(v_zt9,v_zt9));
         v_at3 = _mm_mul_ps(v_anorm,_mm_add_ps(v_at3,v_zt7));
/* convert to double precision before accumulating */
         v_d = _mm_cvtps_pd(v_at3);
         v_ws = _mm_add_pd(v_ws,v_d);
         v_it = _mm_srli_si128((__m128i)v_at3,8);
         v_d = _mm_cvtps_pd((__m128)v_it);
         v_ws = _mm_add_pd(v_ws,v_d);
/*       zt4 += dth*(dky*zt1); */
/*       zt5 += dth*(dkx*zt1); */
         v_zt1 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt1));
         v_zt4 = _mm_add_ps(v_zt4,v_zt1);
/*       zt6 -= dth*(dkx*zt2 + dky*zt3); */
         v_zt2 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt2));
         v_at1 = _mm_movelh_ps(v_zt2,v_zt2);
         v_zt2 = _mm_add_ps(_mm_movehl_ps(v_zt2,v_zt2),v_at1);
         v_zt6 = _mm_movelh_ps(_mm_sub_ps(v_zt6,v_zt2),v_zero);
/*       bxy[4*j+k1] = zt4;   */
/*       bxy[1+4*j+k1] = zt5; */
/*       bxy[2+4*j+k1] = zt6; */
         _mm_store_ps((float *)&bxy[4*j+k1],v_zt4);
         _mm_store_ps((float *)&bxy[2+4*j+k1],v_zt6);
/*       wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
         v_zt4 = _mm_mul_ps(v_zt4,v_zt4);
         v_zt4 = _mm_add_ps(v_zt4,_mm_mul_ps(v_zt6,v_zt6));
         v_at4 = _mm_mul_ps(v_anorm,_mm_add_ps(v_at4,v_zt4));
/* convert to double precision before accumulating */
         v_d = _mm_cvtps_pd(v_at4);
         v_wp = _mm_add_pd(v_wp,v_d);
         v_it = _mm_srli_si128((__m128i)v_at4,8);
         v_d = _mm_cvtps_pd((__m128)v_it);
         v_wp = _mm_add_pd(v_wp,v_d);
      }
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = 4*nxvh*k;
      k1 = 4*nxvh*ny - kj;
      dky = dny*(float) k;
      afdt = adt*cimagf(ffc[kk]);
/* update magnetic field half time step */
      zt1 = -cimagf(exy[2+kj]) + crealf(exy[2+kj])*_Complex_I;
      zt3 = -cimagf(exy[kj]) + crealf(exy[kj])*_Complex_I;
      zt4 = bxy[kj] - dth*(dky*zt1);
      zt6 = bxy[2+kj] + dth*(dky*zt3);
/* update electric field whole time step */
      zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
      zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
      zt7 = exy[kj] + cdt*(dky*zt1) - afdt*cu[kj];
      zt9 = exy[2+kj] - cdt*(dky*zt3) - afdt*cu[2+kj];
/* update magnetic field half time step and store electric field */
      zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
      zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
      exy[kj] = zt7;
      exy[1+kj] = zero;
      exy[2+kj] = zt9;
      exy[3+kj] = zero;
      ws += anorm*(zt7*conjf(zt7) + zt9*conjf(zt9));
      zt4 -= dth*(dky*zt1);
      zt6 += dth*(dky*zt3);
      bxy[kj] = zt4;
      bxy[1+kj] = zero;
      bxy[2+kj] = zt6;
      bxy[3+kj] = zero;
      wp += anorm*(zt4*conjf(zt4) + zt6*conjf(zt6));
      bxy[k1] = zero;
      bxy[1+k1] = zero;
      bxy[2+k1] = zero;
      bxy[3+k1] = zero;
      exy[k1] = zero;
      exy[1+k1] = zero;
      exy[2+k1] = zero;
      exy[3+k1] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = 4*nxvh*nyh;
   for (j = 1; j < nxh; j++) {
/*    dkx = dnx*(float) j; */
      v_dkx = _mm_mul_ps(v_dnx,_mm_cvtepi32_ps(_mm_set1_epi32(j)));
      v_dkx = _mm_movelh_ps(v_zero,v_dkx);
/*    afdt = adt*cimagf(ffc[j]); */
      v_afdt = _mm_loadl_pi(v_zero,(__m64 *)&ffc[j]);
      v_afdt = _mm_movelh_ps(v_afdt,v_afdt);
      v_afdt = _mm_mul_ps(v_adt,_mm_shuffle_ps(v_afdt,v_afdt,245));
/* update magnetic field half time step */
/*    zt1 = -cimagf(exy[2+4*j]) + crealf(exy[2+4*j])*_Complex_I; */
/*    zt2 = -cimagf(exy[1+4*j]) + crealf(exy[1+4*j])*_Complex_I; */
      v_zt7 = _mm_load_ps((float *)&exy[4*j]);
      v_zt2 = _mm_mul_ps(v_zt7,v_n);
      v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
      v_zt9 = _mm_load_ps((float *)&exy[2+4*j]);
      v_zt1 = _mm_mul_ps(v_zt9,v_n);
      v_zt1 = _mm_movelh_ps(v_zt1,v_zt1);
      v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*    zt5 = bxy[1+4*j] + dth*(dkx*zt1); */
      v_zt4 = _mm_load_ps((float *)&bxy[4*j]);
      v_zt1 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt1));
      v_zt4 = _mm_add_ps(v_zt4,v_zt1);
/*    zt6 = bxy[2+4*j] - dth*(dkx*zt2); */
      v_zt6 = _mm_load_ps((float *)&bxy[2+4*j]);
      v_zt2 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt2));
      v_zt2 = _mm_movehl_ps(v_zero,v_zt2);
      v_zt6 = _mm_movelh_ps(_mm_sub_ps(v_zt6,v_zt2),v_zero);
/* update electric field whole time step */
/*    zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*    zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
      v_zt2 = _mm_mul_ps(v_zt4,v_n);
      v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
      v_zt1 = _mm_mul_ps(v_zt6,v_n);
      v_zt1 = _mm_movelh_ps(v_zt1,v_zt1);
      v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*    zt8 = exy[1+4*j] - cdt*(dkx*zt1) - afdt*cu[1+4*j]; */
      v_at2 = _mm_load_ps((float *)&cu[4*j]);
      v_zt1 = _mm_mul_ps(v_cdt,_mm_mul_ps(v_dkx,v_zt1));
      v_zt1 = _mm_add_ps(_mm_mul_ps(v_afdt,v_at2),v_zt1);
      v_zt7 = _mm_sub_ps(v_zt7,v_zt1);
/*    zt9 = exy[2+4*j] + cdt*(dkx*zt2) - afdt*cu[2+4*j]; */
      v_at2 = _mm_load_ps((float *)&cu[2+4*j]);
      v_zt2 = _mm_mul_ps(v_cdt,_mm_mul_ps(v_dkx,v_zt2));
      v_zt2 = _mm_movehl_ps(v_zero,v_zt2);
      v_zt2 = _mm_sub_ps(v_zt2,_mm_mul_ps(v_afdt,v_at2));
      v_zt9 = _mm_movelh_ps(_mm_add_ps(v_zt9,v_zt2),v_zero);
/* update magnetic field half time step and store electric field */
/*    zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*    zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
      v_zt2 = _mm_mul_ps(v_zt7,v_n);
      v_zt2 = _mm_shuffle_ps(v_zt2,v_zt2,177);
      v_zt1 = _mm_mul_ps(v_zt9,v_n);
      v_zt1 = _mm_movelh_ps(v_zt1,v_zt1);
      v_zt1 = _mm_shuffle_ps(v_zt1,v_zt1,177);
/*    exy[4*j] = zero;  */
/*    exy[1+4*j] = zt8; */
/*    exy[2+4*j] = zt9; */
      _mm_store_ps((float *)&exy[4*j],v_zt7);
      _mm_store_ps((float *)&exy[2+4*j],v_zt9);
/*    ws += anorm*(zt8*conjf(zt8) + zt9*conjf(zt9)); */
      v_zt7 = _mm_mul_ps(v_zt7,v_zt7);
      v_at3 = _mm_add_ps(v_zt7,_mm_mul_ps(v_zt9,v_zt9));
      v_at3 = _mm_mul_ps(v_anorm,v_at3);
/* convert to double precision before accumulating */
      v_d = _mm_cvtps_pd(v_at3);
      v_ws = _mm_add_pd(v_ws,v_d);
      v_it = _mm_srli_si128((__m128i)v_at3,8);
      v_d = _mm_cvtps_pd((__m128)v_it);
      v_ws = _mm_add_pd(v_ws,v_d);
/*    zt5 += dth*(dkx*zt1); */
      v_zt1 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt1));
      v_zt4 = _mm_add_ps(v_zt4,v_zt1);
/*    zt6 -= dth*(dkx*zt2); */
      v_zt2 = _mm_mul_ps(v_dth,_mm_mul_ps(v_dkx,v_zt2));
      v_zt2 = _mm_movehl_ps(v_zero,v_zt2);
      v_zt6 = _mm_movelh_ps(_mm_sub_ps(v_zt6,v_zt2),v_zero);
/*    bxy[4*j] = zero;  */
/*    bxy[1+4*j] = zt5; */
/*    bxy[2+4*j] = zt6; */
      _mm_store_ps((float *)&bxy[4*j],v_zt4);
      _mm_store_ps((float *)&bxy[2+4*j],v_zt6);
/*    wp += anorm*(zt5*conjf(zt5) + zt6*conjf(zt6)); */
      v_zt4 = _mm_mul_ps(v_zt4,v_zt4);
      v_at4 = _mm_add_ps(v_zt4,_mm_mul_ps(v_zt6,v_zt6));
      v_at4 = _mm_mul_ps(v_anorm,v_at4);
/* convert to double precision before accumulating */
      v_d = _mm_cvtps_pd(v_at4);
      v_wp = _mm_add_pd(v_wp,v_d);
      v_it = _mm_srli_si128((__m128i)v_at4,8);
      v_d = _mm_cvtps_pd((__m128)v_it);
      v_wp = _mm_add_pd(v_wp,v_d);
/*    bxy[4*j+k1] = zero;   */ 
/*    bxy[1+4*j+k1] = zero; */
/*    bxy[2+4*j+k1] = zero; */
      _mm_store_ps((float *)&bxy[4*j+k1],v_zero);
      _mm_store_ps((float *)&bxy[2+4*j+k1],v_zero);
/*    exy[4*j+k1] = zero;   */
/*    exy[1+4*j+k1] = zero; */
/*    exy[2+4*j+k1] = zero; */
      _mm_store_ps((float *)&exy[4*j+k1],v_zero);
      _mm_store_ps((float *)&exy[2+4*j+k1],v_zero);
   }
/* bxy[0] = zero; */
/* bxy[1] = zero; */
/* bxy[2] = zero; */
   _mm_store_ps((float *)&bxy[0],v_zero);
   _mm_store_ps((float *)&bxy[2],v_zero);
/* exy[0] = zero; */
/* exy[1] = zero; */
/* exy[2] = zero; */
   _mm_store_ps((float *)&exy[0],v_zero);
   _mm_store_ps((float *)&exy[2],v_zero);
/* bxy[k1] = zero;   */
/* bxy[1+k1] = zero; */
/* bxy[2+k1] = zero; */
   _mm_store_ps((float *)&bxy[k1],v_zero);
   _mm_store_ps((float *)&bxy[2+k1],v_zero);
/* exy[k1] = zero;   */
/* exy[1+k1] = zero; */
/* exy[2+k1] = zero; */
   _mm_store_ps((float *)&exy[k1],v_zero);
   _mm_store_ps((float *)&exy[2+k1],v_zero);
/* *wf = ws*(float) (nx*ny); */
   _mm_store_pd(&dd[0],v_ws);
   for (j = 1; j < 2; j++) {
      dd[0] += dd[j];
   }
   *wf = (ws + dd[0])*(float) (nx*ny);
/* *wm = c2*wp*(float) (nx*ny); */
   _mm_store_pd(&dd[0],v_wp);
   for (j = 1; j < 2; j++) {
      dd[0] += dd[j];
   }
   *wm = c2*(wp + dd[0])*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void csse2emfield2(float complex fxy[], float complex exy[],
                   float complex ffc[], int isign, int nx, int ny,
                   int nxvh, int nyv, int nxhd, int nyhd) {
/* this subroutine either adds complex vector fields if isign > 0
   or copies complex vector fields if isign < 0
   includes additional smoothing
   requires SSE2, fxy, exy, ffc need to be 16 byte aligned
   nxhd, nxvh need to be a multiple of 2
   fxy, exy, need to have 4 components
local data                                                 */
   int j, k, nxh, nyh, k1, kk, kj;
   __m128 v_at1, v_zero, v_zt1, v_zt2;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   v_zero = _mm_set1_ps(0.0f);
/* add the fields */
   if (isign > 0) {
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = 4*nxvh*k;
         k1 = 4*nxvh*ny - kj;
         for (j = 0; j < nxh; j++) {
/*          at1 = cimagf(ffc[j+kk]); */
            v_at1 = _mm_loadl_pi(v_zero,(__m64 *)&ffc[j+kk]);
            v_at1 = _mm_movelh_ps(v_at1,v_at1);
            v_at1 = _mm_shuffle_ps(v_at1,v_at1,245);
/*          for (i = 0; i < 3; i++) {             */
/*            fxy[i+4*j+kj] += exy[i+4*j+kj]*at1; */
/*            fxy[i+4*j+k1] += exy[i+4*j+k1]*at1; */
/*          }                                     */
            v_zt1 = _mm_load_ps((float *)&exy[4*j+kj]);
            v_zt2 = _mm_load_ps((float *)&fxy[4*j+kj]);
            v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
            _mm_store_ps((float *)&fxy[4*j+kj],v_zt2);
            v_zt1 = _mm_load_ps((float *)&exy[2+4*j+kj]);
            v_zt2 = _mm_load_ps((float *)&fxy[2+4*j+kj]);
            v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
            _mm_store_ps((float *)&fxy[2+4*j+kj],v_zt2);
            v_zt1 = _mm_load_ps((float *)&exy[4*j+k1]);
            v_zt2 = _mm_load_ps((float *)&fxy[4*j+k1]);
            v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
            _mm_store_ps((float *)&fxy[4*j+k1],v_zt2);
            v_zt1 = _mm_load_ps((float *)&exy[2+4*j+k1]);
            v_zt2 = _mm_load_ps((float *)&fxy[2+4*j+k1]);
            v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
            _mm_store_ps((float *)&fxy[2+4*j+k1],v_zt2);
         }
      }
      k1 = 4*nxvh*nyh;
      for (j = 0; j < nxh; j++) {
/*       at1 = cimagf(ffc[j]); */
         v_at1 = _mm_loadl_pi(v_zero,(__m64 *)&ffc[j]);
         v_at1 = _mm_movelh_ps(v_at1,v_at1);
         v_at1 = _mm_shuffle_ps(v_at1,v_at1,245);
/*       for (i = 0; i < 3; i++) {              */
/*          fxy[i+4*j] += exy[i+4*j]*at1;       */
/*          fxy[i+4*j+k1] += exy[i+4*j+k1]*at1; */
/*       }                                      */
         v_zt1 = _mm_load_ps((float *)&exy[4*j]);
         v_zt2 = _mm_load_ps((float *)&fxy[4*j]);
         v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
         _mm_store_ps((float *)&fxy[4*j],v_zt2);
         v_zt1 = _mm_load_ps((float *)&exy[2+4*j]);
         v_zt2 = _mm_load_ps((float *)&fxy[2+4*j]);
         v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
         _mm_store_ps((float *)&fxy[2+4*j],v_zt2);
         v_zt1 = _mm_load_ps((float *)&exy[4*j+k1]);
         v_zt2 = _mm_load_ps((float *)&fxy[4*j+k1]);
         v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
         _mm_store_ps((float *)&fxy[4*j+k1],v_zt2);
         v_zt1 = _mm_load_ps((float *)&exy[2+4*j+k1]);
         v_zt2 = _mm_load_ps((float *)&fxy[2+4*j+k1]);
         v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
         _mm_store_ps((float *)&fxy[2+4*j+k1],v_zt2);
      }
   }
/* copy the fields */
   else if (isign < 0) {
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = 4*nxvh*k;
         k1 = 4*nxvh*ny - kj;
         for (j = 0; j < nxh; j++) {
/*          at1 = cimagf(ffc[j+kk]); */
            v_at1 = _mm_loadl_pi(v_zero,(__m64 *)&ffc[j+kk]);
            v_at1 = _mm_movelh_ps(v_at1,v_at1);
            v_at1 = _mm_shuffle_ps(v_at1,v_at1,245);
/*          for (i = 0; i < 3; i++) {             */
/*             fxy[i+4*j+kj] = exy[i+4*j+kj]*at1; */
/*             fxy[i+4*j+k1] = exy[i+4*j+k1]*at1; */
/*          }                                     */
            v_zt1 = _mm_load_ps((float *)&exy[4*j+kj]);
            v_zt2 = _mm_mul_ps(v_zt1,v_at1);
            _mm_store_ps((float *)&fxy[4*j+kj],v_zt2);
            v_zt1 = _mm_load_ps((float *)&exy[2+4*j+kj]);
            v_zt2 = _mm_mul_ps(v_zt1,v_at1);
            _mm_store_ps((float *)&fxy[2+4*j+kj],v_zt2);
            v_zt1 = _mm_load_ps((float *)&exy[4*j+k1]);
            v_zt2 = _mm_mul_ps(v_zt1,v_at1);
            _mm_store_ps((float *)&fxy[4*j+k1],v_zt2);
            v_zt1 = _mm_load_ps((float *)&exy[2+4*j+k1]);
            v_zt2 = _mm_mul_ps(v_zt1,v_at1);
            _mm_store_ps((float *)&fxy[2+4*j+k1],v_zt2);
         }
      }
      k1 = 4*nxvh*nyh;
      for (j = 0; j < nxh; j++) {
/*       at1 = cimagf(ffc[j]); */
         v_at1 = _mm_loadl_pi(v_zero,(__m64 *)&ffc[j]);
         v_at1 = _mm_movelh_ps(v_at1,v_at1);
         v_at1 = _mm_shuffle_ps(v_at1,v_at1,245);
/*       for (i = 0; i < 3; i++) {             */
/*          fxy[i+4*j] = exy[i+4*j]*at1;       */
/*          fxy[i+4*j+k1] = exy[i+4*j+k1]*at1; */
/*       }                                     */
         v_zt1 = _mm_load_ps((float *)&exy[4*j]);
         v_zt2 = _mm_mul_ps(v_zt1,v_at1);
         _mm_store_ps((float *)&fxy[4*j],v_zt2);
         v_zt1 = _mm_load_ps((float *)&exy[2+4*j]);
         v_zt2 = _mm_mul_ps(v_zt1,v_at1);
         _mm_store_ps((float *)&fxy[2+4*j],v_zt2);
         v_zt1 = _mm_load_ps((float *)&exy[4*j+k1]);
         v_zt2 = _mm_mul_ps(v_zt1,v_at1);
         _mm_store_ps((float *)&fxy[4*j+k1],v_zt2);
         v_zt1 = _mm_load_ps((float *)&exy[2+4*j+k1]);
         v_zt2 = _mm_mul_ps(v_zt1,v_at1);
         _mm_store_ps((float *)&fxy[2+4*j+k1],v_zt2);
      }
   }
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
void csse2fft2r3x(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nyi,
                  int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
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
   f[k][2*j],f[k][2*j+1] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][0],f[k][1] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   f[1][0] = real part of mode nx/2,0 and
   f[1][ny/2] = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
   requires SSE2, f needs to be 16 byte aligned
   f needs to have 4 components
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, jj, j1, k1, k2, ns, ns2, km, kmr, joff;
   float ani;
/* float complex t1, t2, t3; */
   __m128 v_m, v_n, v_t1, v_t2, v_t3, v_t4, v_t5, v_ani;
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
      joff = 4*nxhd*k;
      for (j = 0; j < nxh; j++) {
         v_t1 = _mm_load_ps((float *)&f[4*j+joff]);
         v_t2 = _mm_load_ps((float *)&f[2+4*j+joff]);
         v_t3 = _mm_movelh_ps(v_t1,v_t2);
         v_t3 = _mm_shuffle_ps(v_t3,v_t3,216);
         _mm_store_ps((float *)&f[4*j+joff],v_t3);
         v_t3 = _mm_movehl_ps(v_t2,v_t1);
         v_t3 = _mm_shuffle_ps(v_t3,v_t3,216);
         _mm_store_ps((float *)&f[2+4*j+joff],v_t3);
      } 
   }
/* bit-reverse array elements in x */
   nrx = nxhy/nxh;
   for (j = 0; j < nxh; j++) {
      j1 = (mixup[j] - 1)/nrx;
      if (j >= j1)
         continue;
      for (k = nyi-1; k < nyt; k++) {
         joff = 4*nxhd*k;
/*       t1 = f[4*j1+joff];   */
/*       t2 = f[1+4*j1+joff]; */
/*       t3 = f[2+4*j1+joff]; */
         v_t1 = _mm_load_ps((float *)&f[4*j1+joff]);
         v_t3 = _mm_load_ps((float *)&f[2+4*j1+joff]);
/*       f[4*j1+joff] = f[4*j+joff];     */
/*       f[1+4*j1+joff] = f[1+4*j+joff]; */
/*       f[2+4*j1+joff] = f[2+4*j+joff]; */
         v_t2 = _mm_load_ps((float *)&f[4*j+joff]);
         _mm_store_ps((float *)&f[4*j1+joff],v_t2);
         v_t2 = _mm_load_ps((float *)&f[2+4*j+joff]);
         _mm_store_ps((float *)&f[2+4*j1+joff],v_t2);
/*       f[4*j+joff] = t1;   */
/*       f[1+4*j+joff] = t2; */
/*       f[2+4*j+joff] = t3; */
         _mm_store_ps((float *)&f[4*j+joff],v_t1);
         _mm_store_ps((float *)&f[2+4*j+joff],v_t3);
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
         k1 = 4*ns2*k;
         k2 = k1 + 4*ns;
         for (i = nyi-1; i < nyt; i++) {
            joff = 4*nxhd*i;
            for (j = 0; j < ns; j++) {
/*             t1 = sct[kmr*j]; */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_movelh_ps(v_t1,v_t1);
/*             t2 = t1*f[4*j+k2+joff];   */
/*             t3 = t1*f[1+4*j+k2+joff]; */
               v_t2 = _mm_load_ps((float *)&f[4*j+k2+joff]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             t4 = t1*f[2+4*j+k2+joff]; */
               v_t4 = _mm_load_ps((float *)&f[2+4*j+k2+joff]);
               v_t3 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t4 = _mm_shuffle_ps(v_t4,v_t4,177);
               v_t4 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t4 = _mm_add_ps(v_t3,_mm_mul_ps(v_t4,v_m));
/*             f[4*j+k2+joff] = f[4*j+k1+joff] - t2;     */
/*             f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3; */
               v_t3 = _mm_load_ps((float *)&f[4*j+k1+joff]);
               v_t5 = _mm_sub_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[4*j+k2+joff],v_t5);
/*             f[4*j+k1+joff] += t2;   */
/*             f[1+4*j+k1+joff] += t3; */
               v_t2 = _mm_add_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[4*j+k1+joff],v_t2);
/*             f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4; */
               v_t3 = _mm_load_ps((float *)&f[2+4*j+k1+joff]);
               v_t5 = _mm_sub_ps(v_t3,v_t4);
               _mm_store_ps((float *)&f[2+4*j+k2+joff],v_t5);
/*             f[2+4*j+k1+joff] += t4; */
               v_t4 = _mm_add_ps(v_t3,v_t4);
               _mm_store_ps((float *)&f[2+4*j+k1+joff],v_t4);
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
      joff = 4*nxhd*k;
      for (j = 1; j < nxhh; j++) {
/*       t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I; */
         v_t3 = _mm_loadl_pi(v_t3,(__m64 *)&sct[kmr*j]);
         v_t3 = _mm_movelh_ps(v_t3,v_t3);
         v_t3 = _mm_mul_ps(v_t3,v_m);
         v_t3 = _mm_shuffle_ps(v_t3,v_t3,177);
/*       t2 = conjf(f[jj+4*(nxh-j)+joff]); */
         v_t2 = _mm_load_ps((float *)&f[4*(nxh-j)+joff]);
         v_t2 = _mm_mul_ps(v_t2,v_n);
         v_t4 = _mm_load_ps((float *)&f[2+4*(nxh-j)+joff]);
         v_t4 = _mm_mul_ps(v_t4,v_n);
/* first block, jj=1:2 */
/*       t1 = f[jj+4*j+joff] + t2; */
         v_t5 = _mm_load_ps((float *)&f[4*j+joff]);
         v_t1 = _mm_add_ps(v_t5,v_t2);
/*       t2 = (f[jj+4*j+joff] - t2)*t3; */
         v_t2 = _mm_sub_ps(v_t5,v_t2);
         v_t5 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,160));
         v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
         v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,245));
         v_t2 = _mm_add_ps(v_t5,_mm_mul_ps(v_t2,v_m));
/*       f[jj+4*j+joff] = ani*(t1 + t2); */
         v_t5 = _mm_mul_ps(v_ani,_mm_add_ps(v_t1,v_t2));
         _mm_store_ps((float *)&f[4*j+joff],v_t5);
/*       f[jj+4*(nxh-j)+joff] = ani*conjf(t1 - t2); */
         v_t5 = _mm_mul_ps(v_ani,_mm_mul_ps(_mm_sub_ps(v_t1,v_t2),v_n));
         _mm_store_ps((float *)&f[4*(nxh-j)+joff],v_t5);
/* second block, jj=3:4 */
/*       t1 = f[jj+4*j+joff] + t2; */
         v_t5 = _mm_load_ps((float *)&f[2+4*j+joff]);
         v_t1 = _mm_add_ps(v_t5,v_t4);
/*       t2 = (f[jj+4*j+joff] - t2)*t3; */
         v_t4 = _mm_sub_ps(v_t5,v_t4);
         v_t5 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t3,v_t3,160));
         v_t4 = _mm_shuffle_ps(v_t4,v_t4,177);
         v_t4 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t3,v_t3,245));
         v_t4 = _mm_add_ps(v_t5,_mm_mul_ps(v_t4,v_m));
/*       f[jj+4*j+joff] = ani*(t1 + t2); */
         v_t5 = _mm_mul_ps(v_ani,_mm_add_ps(v_t1,v_t4));
         _mm_store_ps((float *)&f[2+4*j+joff],v_t5);
/*       f[jj+4*(nxh-j)+joff] = ani*conjf(t1 - t2); */
         v_t5 = _mm_mul_ps(v_ani,_mm_mul_ps(_mm_sub_ps(v_t1,v_t4),v_n));
         _mm_store_ps((float *)&f[2+4*(nxh-j)+joff],v_t5);
      }
   }
   ani = 2.0*ani;
   for (k = nyi-1; k < nyt; k++) {
      joff = 4*nxhd*k;
      for (jj = 0; jj < 3; jj++) {
         f[jj+4*nxhh+joff] = ani*conjf(f[jj+4*nxhh+joff]);
         f[jj+joff] = ani*((crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
/* scramble coefficients */
L140: kmr = nxy/nx;
   for (k = nyi-1; k < nyt; k++) {
      joff = 4*nxhd*k;
      for (j = 1; j < nxhh; j++) {
/*       t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I; */
         v_t3 = _mm_loadl_pi(v_t3,(__m64 *)&sct[kmr*j]);
         v_t3 = _mm_movelh_ps(v_t3,v_t3);
         v_t3 = _mm_shuffle_ps(v_t3,v_t3,177);
/*       t2 = conjf(f[jj+4*(nxh-j)+joff]); */
         v_t2 = _mm_load_ps((float *)&f[4*(nxh-j)+joff]);
         v_t2 = _mm_mul_ps(v_t2,v_n);
         v_t4 = _mm_load_ps((float *)&f[2+4*(nxh-j)+joff]);
         v_t4 = _mm_mul_ps(v_t4,v_n);
/* first block, jj=1:2 */
/*       t1 = f[jj+4*j+joff] + t2; */
         v_t5 = _mm_load_ps((float *)&f[4*j+joff]);
         v_t1 = _mm_add_ps(v_t5,v_t2);
/*       t2 = (f[jj+4*j+joff] - t2)*t3; */
         v_t2 = _mm_sub_ps(v_t5,v_t2);
         v_t5 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,160));
         v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
         v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t3,v_t3,245));
         v_t2 = _mm_add_ps(v_t5,_mm_mul_ps(v_t2,v_m));
/*       f[jj+4*j+joff] = t1 + t2; */
         v_t5 = _mm_add_ps(v_t1,v_t2);
         _mm_store_ps((float *)&f[4*j+joff],v_t5);
/*       f[jj+4*(nxh-j)+joff] = conjf(t1 - t2); */
         v_t5 = _mm_mul_ps(_mm_sub_ps(v_t1,v_t2),v_n);
         _mm_store_ps((float *)&f[4*(nxh-j)+joff],v_t5);
/* second block, jj=3:4 */
/*       t1 = f[jj+4*j+joff] + t2; */
         v_t5 = _mm_load_ps((float *)&f[2+4*j+joff]);
         v_t1 = _mm_add_ps(v_t5,v_t4);
/*       t2 = (f[jj+4*j+joff] - t2)*t3; */
         v_t4 = _mm_sub_ps(v_t5,v_t4);
         v_t5 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t3,v_t3,160));
         v_t4 = _mm_shuffle_ps(v_t4,v_t4,177);
         v_t4 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t3,v_t3,245));
         v_t4 = _mm_add_ps(v_t5,_mm_mul_ps(v_t4,v_m));
/*       f[jj+4*j+joff] = t1 + t2; */
         v_t5 = _mm_add_ps(v_t1,v_t4);
          _mm_store_ps((float *)&f[2+4*j+joff],v_t5);
/*       f[jj+4*(nxh-j)+joff] = conjf(t1 - t2); */
         v_t5 = _mm_mul_ps(_mm_sub_ps(v_t1,v_t4),v_n);
         _mm_store_ps((float *)&f[2+4*(nxh-j)+joff],v_t5);
      }
   }
   for (k = nyi-1; k < nyt; k++) {
      joff = 4*nxhd*k;
      for (jj = 0; jj < 3; jj++) {
         f[jj+4*nxhh+joff] = 2.0*conjf(f[jj+4*nxhh+joff]);
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
         joff = 4*nxhd*k;
/*       t1 = f[4*j1+joff];   */
/*       t2 = f[1+4*j1+joff]; */
/*       t3 = f[2+4*j1+joff]; */
         v_t1 = _mm_load_ps((float *)&f[4*j1+joff]);
         v_t3 = _mm_load_ps((float *)&f[2+4*j1+joff]);
/*       f[4*j1+joff] = f[4*j+joff];     */
/*       f[1+4*j1+joff] = f[1+4*j+joff]; */
/*       f[2+4*j1+joff] = f[2+4*j+joff]; */
         v_t2 = _mm_load_ps((float *)&f[4*j+joff]);
         _mm_store_ps((float *)&f[4*j1+joff],v_t2);
         v_t2 = _mm_load_ps((float *)&f[2+4*j+joff]);
         _mm_store_ps((float *)&f[2+4*j1+joff],v_t2);
/*       f[4*j+joff] = t1;   */
/*       f[1+4*j+joff] = t2; */
/*       f[2+4*j+joff] = t3; */
         _mm_store_ps((float *)&f[4*j+joff],v_t1);
         _mm_store_ps((float *)&f[2+4*j+joff],v_t3);
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
         k1 = 4*ns2*k;
         k2 = k1 + 4*ns;
         for (i = nyi-1; i < nyt; i++) {
            joff = 4*nxhd*i;
            for (j = 0; j < ns; j++) {
/*             t1 = conjf(sct[kmr*j]); */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_movelh_ps(v_t1,v_t1);
               v_t1 = _mm_mul_ps(v_t1,v_n);
/*             t2 = t1*f[4*j+k2+joff];   */
/*             t3 = t1*f[1+4*j+k2+joff]; */
               v_t2 = _mm_load_ps((float *)&f[4*j+k2+joff]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             t4 = t1*f[2+4*j+k2+joff]; */
               v_t4 = _mm_load_ps((float *)&f[2+4*j+k2+joff]);
               v_t3 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t4 = _mm_shuffle_ps(v_t4,v_t4,177);
               v_t4 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t4 = _mm_add_ps(v_t3,_mm_mul_ps(v_t4,v_m));
/*             f[4*j+k2+joff] = f[4*j+k1+joff] - t2;     */
/*             f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3; */
               v_t3 = _mm_load_ps((float *)&f[4*j+k1+joff]);
               v_t5 = _mm_sub_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[4*j+k2+joff],v_t5);
/*             f[4*j+k1+joff] += t2;   */
/*             f[1+4*j+k1+joff] += t3; */
               v_t2 = _mm_add_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[4*j+k1+joff],v_t2);
/*             f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4; */
               v_t3 = _mm_load_ps((float *)&f[2+4*j+k1+joff]);
               v_t5 = _mm_sub_ps(v_t3,v_t4);
               _mm_store_ps((float *)&f[2+4*j+k2+joff],v_t5);
/*             f[2+4*j+k1+joff] += t4; */
               v_t4 = _mm_add_ps(v_t3,v_t4);
               _mm_store_ps((float *)&f[2+4*j+k1+joff],v_t4);
            }
         }
      }
      ns = ns2;
   }
/* swap complex components */
   for (k = nyi-1; k < nyt; k++) {
      joff = 4*nxhd*k;
      for (j = 0; j < nxh; j++) {
         v_t1 = _mm_load_ps((float *)&f[4*j+joff]);
         v_t2 = _mm_load_ps((float *)&f[2+4*j+joff]);
         v_t1 = _mm_shuffle_ps(v_t1,v_t1,216);
         v_t2 = _mm_shuffle_ps(v_t2,v_t2,216);
         v_t3 = _mm_movelh_ps(v_t1,v_t2);
         _mm_store_ps((float *)&f[4*j+joff],v_t3);
         v_t3 = _mm_movehl_ps(v_t2,v_t1);
         _mm_store_ps((float *)&f[2+4*j+joff],v_t3);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void csse2fft2r3y(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxi,
                  int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
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
  f[k][2*j],f[k][2*j+1] = real, imaginary part of mode j,k, where
  0 <= j < nx/2 and 0 <= k < ny, except for
  f[k][0],f[k][1] = real, imaginary part of mode nx/2,k, where
  ny/2+1 <= k < ny, and
  f[1][0] = real part of mode nx/2,0 and
  f[1][ny/2] = real part of mode nx/2,ny/2
  written by viktor k. decyk, ucla
  requires SSE2, f needs to be 16 byte aligned
  f needs to have 4 components
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
   float complex t1;
   __m128 v_m, v_n, v_t1, v_t2, v_t3, v_t4, v_t5;
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
      joff = 4*nxhd*k;
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      k1 = 4*nxhd*k1;
      for (j = nxi-1; j < nxt; j++) {
/*       t1 = f[4*j+k1];   */
/*       t2 = f[1+4*j+k1]; */
/*       t3 = f[2+4*j+k1]; */
         v_t1 = _mm_load_ps((float *)&f[4*j+k1]);
         v_t3 = _mm_load_ps((float *)&f[2+4*j+k1]);
/*       f[4*j+k1] = f[4*j+joff];     */
/*       f[1+4*j+k1] = f[1+4*j+joff]; */
/*       f[2+4*j+k1] = f[2+4*j+joff]; */
         v_t2 = _mm_load_ps((float *)&f[4*j+joff]);
         _mm_store_ps((float *)&f[4*j+k1],v_t2);
         v_t2 = _mm_load_ps((float *)&f[2+4*j+joff]);
         _mm_store_ps((float *)&f[2+4*j+k1],v_t2);
/*       f[4*j+joff] = t1;   */
/*       f[1+4*j+joff] = t2; */
/*       f[2+4*j+joff] = t3; */
         _mm_store_ps((float *)&f[4*j+joff],v_t1);
         _mm_store_ps((float *)&f[2+4*j+joff],v_t3);
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
            j1 = 4*nxhd*(j + k1);
            j2 = 4*nxhd*(j + k2);
/*          t1 = sct[kmr*j]; */
            v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
            v_t1 = _mm_movelh_ps(v_t1,v_t1);
            for (i = nxi-1; i < nxt; i++) {
/*             t2 = t1*f[4*i+j2];   */
/*             t3 = t1*f[1+4*i+j2]; */
               v_t2 = _mm_load_ps((float *)&f[4*i+j2]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             t4 = t1*f[2+4*i+j2]; */
               v_t4 = _mm_load_ps((float *)&f[2+4*i+j2]);
               v_t3 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t4 = _mm_shuffle_ps(v_t4,v_t4,177);
               v_t4 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t4 = _mm_add_ps(v_t3,_mm_mul_ps(v_t4,v_m));
/*             f[4*i+j2] = f[4*i+j1] - t2;     */
/*             f[1+4*i+j2] = f[1+4*i+j1] - t3; */
               v_t3 = _mm_load_ps((float *)&f[4*i+j1]);
               v_t5 = _mm_sub_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[4*i+j2],v_t5);
/*             f[4*i+j1] += t2;   */
/*             f[1+4*i+j1] += t3; */
               v_t2 = _mm_add_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[4*i+j1],v_t2);
/*             f[2+4*i+j2] = f[2+4*i+j1] - t4; */
               v_t3 = _mm_load_ps((float *)&f[2+4*i+j1]);
               v_t5 = _mm_sub_ps(v_t3,v_t4);
               _mm_store_ps((float *)&f[2+4*i+j2],v_t5);
/*             f[2+4*i+j1] += t4; */
               v_t4 = _mm_add_ps(v_t3,v_t4);
               _mm_store_ps((float *)&f[2+4*i+j1],v_t4);
            }
         }
      }
      ns = ns2;
   }
/* unscramble modes kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      if (nxi==1) {
         joff = 4*nxhd*k;
         k1 = 4*nxhd*ny - joff;
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
         joff = 4*nxhd*k;
         k1 = 4*nxhd*ny - joff;
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
      joff = 4*nxhd*k;
      k1 = (mixup[k] - 1)/nry;
      if (k >= k1)
         continue;
      k1 = 4*nxhd*k1;
      for (j = nxi-1; j < nxt; j++) {
/*       t1 = f[4*j+k1];   */
/*       t2 = f[1+4*j+k1]; */
/*       t3 = f[2+4*j+k1]; */
         v_t1 = _mm_load_ps((float *)&f[4*j+k1]);
         v_t3 = _mm_load_ps((float *)&f[2+4*j+k1]);
/*       f[4*j+k1] = f[4*j+joff];     */
/*       f[1+4*j+k1] = f[1+4*j+joff]; */
/*       f[2+4*j+k1] = f[2+4*j+joff]; */
         v_t2 = _mm_load_ps((float *)&f[4*j+joff]);
         _mm_store_ps((float *)&f[4*j+k1],v_t2);
         v_t2 = _mm_load_ps((float *)&f[2+4*j+joff]);
         _mm_store_ps((float *)&f[2+4*j+k1],v_t2);
/*       f[4*j+joff] = t1;   */
/*       f[1+4*j+joff] = t2; */
/*       f[2+4*j+joff] = t3; */
         _mm_store_ps((float *)&f[4*j+joff],v_t1);
         _mm_store_ps((float *)&f[2+4*j+joff],v_t3);
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
            j1 = 4*nxhd*(j + k1);
            j2 = 4*nxhd*(j + k2);
/*          t1 = conjf(sct[kmr*j]); */
            v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
            v_t1 = _mm_movelh_ps(v_t1,v_t1);
            v_t1 = _mm_mul_ps(v_t1,v_n);
            for (i = nxi-1; i < nxt; i++) {
/*             t2 = t1*f[4*i+j2];   */
/*             t3 = t1*f[1+4*i+j2]; */
               v_t2 = _mm_load_ps((float *)&f[4*i+j2]);
               v_t3 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t2 = _mm_shuffle_ps(v_t2,v_t2,177);
               v_t2 = _mm_mul_ps(v_t2,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t2 = _mm_add_ps(v_t3,_mm_mul_ps(v_t2,v_m));
/*             t4 = t1*f[2+4*i+j2]; */
               v_t4 = _mm_load_ps((float *)&f[2+4*i+j2]);
               v_t3 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t1,v_t1,160));
               v_t4 = _mm_shuffle_ps(v_t4,v_t4,177);
               v_t4 = _mm_mul_ps(v_t4,_mm_shuffle_ps(v_t1,v_t1,245));
               v_t4 = _mm_add_ps(v_t3,_mm_mul_ps(v_t4,v_m));
/*             f[4*i+j2] = f[4*i+j1] - t2;    */
/*             f[1+4*i+j2] = f[1+4*i+j1] - t3; */
               v_t3 = _mm_load_ps((float *)&f[4*i+j1]);
               v_t5 = _mm_sub_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[4*i+j2],v_t5);
/*             f[4*i+j1] += t2;   */
/*             f[1+4*i+j1] += t3; */
               v_t2 = _mm_add_ps(v_t3,v_t2);
               _mm_store_ps((float *)&f[4*i+j1],v_t2);
/*             f[2+4*i+j2] = f[2+4*i+j1] - t4; */
               v_t3 = _mm_load_ps((float *)&f[2+4*i+j1]);
               v_t5 = _mm_sub_ps(v_t3,v_t4);
               _mm_store_ps((float *)&f[2+4*i+j2],v_t5);
/*             f[2+4*i+j1] += t4; */
               v_t4 = _mm_add_ps(v_t3,v_t4);
               _mm_store_ps((float *)&f[2+4*i+j1],v_t4);
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
void csse2wfft2r3(float complex f[],int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxhd,
                  int nyd, int nxhyd, int nxyhd) {
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
      csse2fft2r3x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                   nxyhd);
/* perform y fft */
      csse2fft2r3y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                   nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      csse2fft2r3y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                   nxyhd);
/* perform x fft */
      csse2fft2r3x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
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
void csse2gbpush23lt_(float *part, float *fxy, float *bxy, float *qbm,
                      float *dt, float *dtc, float *ek, int *idimp,
                      int *nop, int *npe, int *nx, int *ny, int *nxv,
                      int *nyv, int *ipbc) {
   csse2gbpush23lt(part,fxy,bxy,*qbm,*dt,*dtc,ek,*idimp,*nop,*npe,*nx,
                   *ny,*nxv,*nyv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2grbpush23lt_(float *part, float *fxy, float *bxy, float *qbm,
                       float *dt, float *dtc, float *ci, float *ek,
                       int *idimp, int *nop, int *npe, int *nx, int *ny,
                       int *nxv, int *nyv, int *ipbc) {
   csse2grbpush23lt(part,fxy,bxy,*qbm,*dt,*dtc,*ci,ek,*idimp,*nop,*npe,
                    *nx,*ny,*nxv,*nyv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2gpost2lt_(float *part, float *q, float *qm, int *nop,
                    int *npe, int *idimp, int *nxv, int *nyv) {
   csse2gpost2lt(part,q,*qm,*nop,*npe,*idimp,*nxv,*nyv);
   return;
}


/*--------------------------------------------------------------------*/
void csse2gjpost2lt_(float *part, float *cu, float *qm, float *dt,
                     int *nop, int *npe, int *idimp, int *nx, int *ny, 
                     int *nxv, int *nyv, int *ipbc) {
   csse2gjpost2lt(part,cu,*qm,*dt,*nop,*npe,*idimp,*nx,*ny,*nxv,*nyv,
                  *ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2grjpost2lt_(float *part, float *cu, float *qm, float *dt,
                      float *ci, int *nop, int *npe, int *idimp,
                      int *nx, int *ny, int *nxv, int *nyv, int *ipbc) {
   csse2grjpost2lt(part,cu,*qm,*dt,*ci,*nop,*npe,*idimp,*nx,*ny,*nxv,
                   *nyv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2dsortp2ylt_(float *parta, float *partb, int *npic, int *idimp,
                      int *nop, int *npe, int *ny1) {
   csse2dsortp2ylt(parta,partb,npic,*idimp,*nop,*npe,*ny1);
   return;
}

/*--------------------------------------------------------------------*/
void csse2bguard2l_(float *bxy, int *nx, int *ny, int *nxe, int *nye) {
   csse2bguard2l(bxy,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void csse2acguard2l_(float *cu, int *nx, int *ny, int *nxe, int *nye) {
   csse2acguard2l(cu,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void csse2aguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye) {
   csse2aguard2l(q,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void csse2pois23_(float complex *q, float complex *fxy, int *isign,
                  float complex *ffc, float *ax, float *ay, float *affp,
                  float *we, int *nx, int *ny, int *nxvh, int *nyv,
                  int *nxhd, int *nyhd) {
   csse2pois23(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*nxvh,*nyv,*nxhd,
               *nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csse2cuperp2_(float complex *cu, int *nx, int *ny, int *nxvh,
                   int *nyv) {
   csse2cuperp2(cu,*nx,*ny,*nxvh,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void csse2ibpois23_(float complex *cu, float complex *bxy,
                    float complex *ffc, float *ci, float *wm, int *nx,
                    int *ny, int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   csse2ibpois23(cu,bxy,ffc,*ci,wm,*nx,*ny,*nxvh,*nyv,*nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csse2maxwel2_(float complex *exy, float complex *bxy,
                   float complex *cu, float complex *ffc, float *ci,
                   float *dt, float *wf, float *wm, int *nx, int *ny,
                   int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   csse2maxwel2(exy,bxy,cu,ffc,*ci,*dt,wf,wm,*nx,*ny,*nxvh,*nyv,*nxhd,
                *nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csse2emfield2_(float complex *fxy, float complex *exy,
                    float complex *ffc, int *isign, int *nx, int *ny,
                    int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   csse2emfield2(fxy,exy,ffc,*isign,*nx,*ny,*nxvh,*nyv,*nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csse2wfft2rx_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *nxhd,
                   int *nyd, int *nxhyd, int *nxyhd) {
   csse2wfft2rx(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,
                *nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csse2wfft2r3_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *nxhd,
                   int *nyd, int *nxhyd, int *nxyhd) {
   csse2wfft2r3(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,
                *nxyhd);
   return;
}
