/* SSE2 C Library for Skeleton 2-1/2D Electromagnetic OpenMP/Vector */
/* PIC Code                                                         */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <xmmintrin.h>
#include "ssembpush2.h"

/*--------------------------------------------------------------------*/
void csse2gbppush23lt(float ppart[], float fxy[], float bxy[],
                      int kpic[], float qbm, float dt, float dtc,
                      float *ek, int idimp, int nppmx, int nx, int ny,
                      int mx, int my, int nxv, int nyv, int mx1,
                      int mxy1, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field. Using the Boris Mover.
   vector/OpenMP version using guard cells
   particles stored in segmented array
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic = number of particles per tile
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
   nyv = second dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   requires SSE2, ppart, fxy, and bxy need to be 16 byte aligned
   nppmx needs to be a multiple of 4, fxy, bxy need to have 4 components
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nps, nn, mm, nm;
   float qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt;
   float anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, vx, vy, vz;
   double sum1, sum2;
   __m128i v_noff, v_moff, v_mxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qtmh, v_dt, v_dtc, v_one, v_two, v_half;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m128 v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 a, b, c, d, e, f, g, h;
   __m128d v_sum1, v_d;
   __attribute__((aligned(16))) unsigned int ll[4];
   __attribute__((aligned(16))) double dd[2];
   __attribute__((aligned(16))) float sfxy[4*MXV*MYV], sbxy[4*MXV*MYV];
/* __attribute__((aligned(16))) float sfxy[4*(mx+1)*(my+1)]; */
/* __attribute__((aligned(16))) float sbxy[4*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
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
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nps,nn,mm,nm,x,y,vx,vy,vz,dxp,dyp,amx, \
amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2, \
rot3,rot4,rot5,rot6,rot7,rot8,rot9,sum1,v_noff,v_moff,v_nn,v_mm,v_it, \
v_x,v_y,v_vx,v_vy,v_vz,v_dxp,v_dyp,v_amx,v_amy,v_dx,v_dy,v_dz,v_at, \
v_d,v_sum1,a,b,c,d,e,f,g,h,ll,dd,sfxy,sbxy) \
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
      for (j = 0; j < mm; j++) {
/*       for (i = 0; i < nn; i++) {                               */
/*          sfxy[4*(i+mxv*j)] = fxy[4*(i+noff+nxv*(j+moff))];     */
/*          sfxy[1+4*(i+mxv*j)] = fxy[1+4*(i+noff+nxv*(j+moff))]; */
/*          sfxy[2+4*(i+mxv*j)] = fxy[2+4*(i+noff+nxv*(j+moff))]; */
/*       }                                                        */
         for (i = 0; i < nn; i++) {
            v_at = _mm_loadu_ps(&fxy[4*(i+noff+nxv*(j+moff))]);
           _mm_storeu_ps(&sfxy[4*(i+mxv*j)],v_at);
         }
      }
      for (j = 0; j < mm; j++) {
/*       for (i = 0; i < nn; i++) {                               */
/*          sbxy[4*(i+mxv*j)] = bxy[4*(i+noff+nxv*(j+moff))];     */
/*          sbxy[1+4*(i+mxv*j)] = bxy[1+4*(i+noff+nxv*(j+moff))]; */
/*          sbxy[2+4*(i+mxv*j)] = bxy[2+4*(i+noff+nxv*(j+moff))]; */
/*       }                                                        */
         for (i = 0; i < nn; i++) {
            v_at = _mm_loadu_ps(&bxy[4*(i+noff+nxv*(j+moff))]);
           _mm_storeu_ps(&sbxy[4*(i+mxv*j)],v_at);
         }
      }
      nps = 4*(npp/4);
      sum1 = 0.0;
      v_sum1 = _mm_set1_pd(0.0);
/* vector loop over particles in blocks of 4 */
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
/*       nm = 4*(nn - noff + mxv*(mm - moff)); */
         v_nn = _mm_sub_epi32(v_nn,v_noff);
         v_mm = _mm_sub_epi32(v_mm,v_moff);
         v_it = _mm_mul_epu32(v_mxv,_mm_srli_si128(v_mm,4));
         v_mm = _mm_mul_epu32(v_mm,v_mxv);
         v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
         v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
/*       amx = 1.0f - dxp; */
/*       amy = 1.0f - dyp; */
         v_amx = _mm_sub_ps(v_one,v_dxp);
         v_amy = _mm_sub_ps(v_one,v_dyp);
        _mm_store_si128((__m128i *)ll,v_nn);
/* find electric field */
/*       nn = nm;                          */
/*       dx = amx*sfxy[nn];                */
/*       dy = amx*sfxy[nn+1];              */
/*       dz = amx*sfxy[nn+2];              */
/*       mm = nn + 4;                      */
/*       dx = amy*(dxp*sfxy[mm] + dx);     */
/*       dy = amy*(dxp*sfxy[mm+1] + dy);   */
/*       dz = amy*(dxp*sfxy[mm+2] + dz);   */
/*       nn += 4*mxv;                      */
/*       acx = amx*sfxy[nn];               */
/*       acy = amx*sfxy[nn+1];             */
/*       acz = amx*sfxy[nn+2];             */
/*       mm = nn + 4;                      */
/*       dx += dyp*(dxp*sfxy[mm] + acx);   */
/*       dy += dyp*(dxp*sfxy[mm+1] + acy); */
/*       dz += dyp*(dxp*sfxy[mm+2] + acz); */
/* find magnetic field */
/*       nn = nm;                          */
/*       ox = amx*sbxy[nn];                */
/*       oy = amx*sbxy[nn+1];              */
/*       oz = amx*sbxy[nn+2];              */
/*       mm = nn + 4;                      */
/*       ox = amy*(dxp*sbxy[mm] + ox);     */
/*       oy = amy*(dxp*sbxy[mm+1] + oy);   */
/*       oz = amy*(dxp*sbxy[mm+2] + oz);   */
/*       nn += 4*mxv;                      */
/*       acx = amx*sbxy[nn];               */
/*       acy = amx*sbxy[nn+1];             */
/*       acz = amx*sbxy[nn+2];             */
/*       mm = nn + 4;                      */
/*       ox += dyp*(dxp*sbxy[mm] + acx);   */
/*       oy += dyp*(dxp*sbxy[mm+1] + acy); */
/*       oz += dyp*(dxp*sbxy[mm+2] + acz); */
/* interpolate electric and magnetic fields for first particle */
         nn = ll[0];
         v_at = _mm_shuffle_ps(v_amx,v_amx,0);
         a = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         e = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,0);
         nn += 4;
         a = _mm_add_ps(a,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         e = _mm_add_ps(e,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,0);
         a = _mm_mul_ps(a,v_at);
         e = _mm_mul_ps(e,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,0);
         a = _mm_add_ps(a,_mm_mul_ps(v_dx,v_at));
         e = _mm_add_ps(e,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for second particle */
         nn = ll[1];
         v_at = _mm_shuffle_ps(v_amx,v_amx,85);
         b = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         f = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,85);
         nn += 4;
         b = _mm_add_ps(b,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         f = _mm_add_ps(f,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,85);
         b = _mm_mul_ps(b,v_at);
         f = _mm_mul_ps(f,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,85);
         b = _mm_add_ps(b,_mm_mul_ps(v_dx,v_at));
         f = _mm_add_ps(f,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for third particle */
         nn = ll[2];
         v_at = _mm_shuffle_ps(v_amx,v_amx,170);
         c = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         g = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,170);
         nn += 4;
         c = _mm_add_ps(c,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         g = _mm_add_ps(g,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,170);
         c = _mm_mul_ps(c,v_at);
         g = _mm_mul_ps(g,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,170);
         c = _mm_add_ps(c,_mm_mul_ps(v_dx,v_at));
         g = _mm_add_ps(g,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for fourth particle */
         nn = ll[3];
         v_at = _mm_shuffle_ps(v_amx,v_amx,255);
         d = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         h = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,255);
         nn += 4;
         d = _mm_add_ps(d,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         h = _mm_add_ps(h,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
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
/*       dx *= qtmh; */
/*       dy *= qtmh; */
/*       dz *= qtmh; */
         v_dx = _mm_mul_ps(a,v_qtmh);
         v_dy = _mm_mul_ps(b,v_qtmh);
         v_dz = _mm_mul_ps(c,v_qtmh);
/* half acceleration */
/*       acx = ppart[j+2*nppmx+npoff] + dx; */
/*       acy = ppart[j+3*nppmx+npoff] + dy; */
/*       acz = ppart[j+4*nppmx+npoff] + dz; */
         a = _mm_add_ps(v_dx,_mm_load_ps(&ppart[j+2*nppmx+npoff]));
         b = _mm_add_ps(v_dy,_mm_load_ps(&ppart[j+3*nppmx+npoff]));
         c = _mm_add_ps(v_dz,_mm_load_ps(&ppart[j+4*nppmx+npoff]));
/* time-centered kinetic energy */
/*       sum1 += (acx*acx + acy*acy + acz*acz); */
         v_at = _mm_add_ps(_mm_mul_ps(a,a),_mm_mul_ps(b,b));
         v_at = _mm_add_ps(v_at,_mm_mul_ps(c,c));
/* convert to double precision before accumulating */
         v_d = _mm_cvtps_pd(v_at);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
         v_it = _mm_srli_si128((__m128i)v_at,8);
         v_d = _mm_cvtps_pd((__m128)v_it);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
/* calculate cyclotron frequency */
/*       omxt = qtmh*ox; */
/*       omyt = qtmh*oy; */
/*       omzt = qtmh*oz; */
         e = _mm_mul_ps(v_qtmh,e);
         f = _mm_mul_ps(v_qtmh,f);
         g = _mm_mul_ps(v_qtmh,g);
/* calculate rotation matrix */
/*       vx = omxt*omxt; */
         v_vx = _mm_mul_ps(e,e);
/*       vy = omyt*omyt; */
         v_vy = _mm_mul_ps(f,f);
/*       vz = omzt*omzt; */
         v_vz = _mm_mul_ps(g,g);
/*       omt = omxt*omxt + omyt*omyt + omzt*omzt; */
         v_at = _mm_add_ps(_mm_add_ps(v_vx,v_vy),v_vz);
/*       anorm = 2.0f/(1.0f + omt); */
         d = _mm_div_ps(v_two,_mm_add_ps(v_one,v_at));
/*       omt = 0.5f*(1.0f - omt); */
         h = _mm_mul_ps(v_half,_mm_sub_ps(v_one,v_at));
/*       vx = (omt + vx)*acx; */
         v_vx = _mm_mul_ps(_mm_add_ps(h,v_vx),a);
/*       vy = (omt + vy)*acy; */
         v_vy = _mm_mul_ps(_mm_add_ps(h,v_vy),b);
/*       vz = (omt + vz)*acz; */
         v_vz = _mm_mul_ps(_mm_add_ps(h,v_vz),c);
/*       omt = omxt*omyt; */
         h = _mm_mul_ps(e,f);
/*       vx = vx + (omzt + omt)*acy; */
         v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_add_ps(h,g),b));
/*       vy = vy + (omt - omzt)*acx; */
         v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_sub_ps(h,g),a));
/*       omt = omxt*omzt;  */
         h = _mm_mul_ps(e,g);
/*       vx = vx + (omt - omyt)*acz; */
         v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_sub_ps(h,f),c));
/*       vz = vz + (omt + omyt)*acx; */
         v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_add_ps(h,f),a));
/*       omt = omyt*omzt; */
         h = _mm_mul_ps(f,g);
/*       vy = vy + (omt + omxt)*acz; */
         v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_add_ps(h,e),c));
/*       vz = vz + (omt - omxt)*acy; */
         v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_sub_ps(h,e),b));
/* new velocity */
/*       vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx; */
/*       vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy; */
/*       vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz; */
         v_vx = _mm_add_ps(v_dx,_mm_mul_ps(v_vx,d));
         v_vy = _mm_add_ps(v_dy,_mm_mul_ps(v_vy,d));
         v_vz = _mm_add_ps(v_dz,_mm_mul_ps(v_vz,d));
/* new position */
/*       dx = x + vx*dtc; */
/*       dy = y + vy*dtc; */
         v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_dtc));
         v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_dtc));
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
/*       ppart[j+4*nppmx+npoff] = vz; */
         _mm_store_ps(&ppart[j+2*nppmx+npoff],v_vx);
         _mm_store_ps(&ppart[j+3*nppmx+npoff],v_vy);
         _mm_store_ps(&ppart[j+4*nppmx+npoff],v_vz);
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 4*(nn - noff + mxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + 4;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += 4*mxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + 4;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + 4;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += 4*mxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + 4;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
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
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* new position */
         dx = x + vx*dtc;
         dy = y + vy*dtc;
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
         ppart[j+4*nppmx+npoff] = vz;
      }
/*    sum2 += sum1; */
      _mm_store_pd(&dd[0],v_sum1);
      for (j = 1; j < 2; j++) {
         dd[0] += dd[j];
      }
      sum2 += (sum1 + dd[0]);
   }
/* normalize kinetic energy */
   *ek += 0.5*sum2;
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void csse2gbppushf23lt(float ppart[], float fxy[], float bxy[],
                       int kpic[], int ncl[], int ihole[], float qbm,
                       float dt, float dtc, float *ek, int idimp,
                       int nppmx, int nx, int ny, int mx, int my,
                       int nxv, int nyv, int mx1, int mxy1, int ntmax,
                       int *irc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field. Using the Boris Mover.
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   vector/OpenMP version using guard cells
   particles stored in segmented array
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
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
   nyv = second dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
   requires SSE2, ppart, fxy, and bxy need to be 16 byte aligned
   nppmx needs to be a multiple of 4, fxy, bxy need to have 4 components
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nps, ih, nh, nn, mm, nm, kk;
   float qtmh, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float x, y, vx, vy, vz;
   double sum1, sum2;
   __m128i v_noff, v_moff, v_mxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qtmh, v_dt, v_dtc, v_one, v_two, v_half;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_st, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m128 v_anx, v_any, v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 v_zero, v_three, v_six;
   __m128 a, b, c, d, e, f, g, h;
   __m128d v_sum1, v_d;
   __attribute__((aligned(16))) unsigned int ll[4], lm[8];
   __attribute__((aligned(16))) unsigned long jj[1];
   __attribute__((aligned(16))) double dd[2];
   __attribute__((aligned(16))) float sfxy[4*MXV*MYV], sbxy[4*MXV*MYV];
/* __attribute__((aligned(16))) float sfxy[4*(mx+1)*(my+1)]; */
/* __attribute__((aligned(16))) float sbxy[4*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   sum2 = 0.0;
   v_mxv = _mm_set1_epi32(mxv);
   v_qtmh = _mm_set1_ps(qtmh);
   v_dt = _mm_set1_ps(dt);
   v_dtc = _mm_set1_ps(dtc);
   v_anx = _mm_set1_ps(anx);
   v_any = _mm_set1_ps(any);
   v_zero = _mm_setzero_ps();
   v_one = _mm_set1_ps(1.0f);
   v_two = _mm_set1_ps(2.0f);
   v_half = _mm_set1_ps(0.5f);
   v_three = _mm_set1_ps(3.0f);
   v_six = _mm_set1_ps(6.0f);
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nps,nn,mm,nm,kk,ih,nh,x,y,vx,vy,vz, \
dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm, \
rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx, \
edgery,sum1,v_noff,v_moff,v_nn,v_mm,v_it,v_x,v_y,v_vx,v_vy,v_vz,v_dxp, \
v_dyp,v_amx,v_amy,v_dx,v_dy,v_dz,v_st,v_at,v_edgelx,v_edgely,v_edgerx, \
v_edgery,v_d,v_sum1,a,b,c,d,e,f,g,h,jj,ll,lm,dd,sfxy,sbxy) \
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
      for (j = 0; j < mm; j++) {
/*       for (i = 0; i < nn; i++) {                               */
/*          sfxy[4*(i+mxv*j)] = fxy[4*(i+noff+nxv*(j+moff))];     */
/*          sfxy[1+4*(i+mxv*j)] = fxy[1+4*(i+noff+nxv*(j+moff))]; */
/*          sfxy[2+4*(i+mxv*j)] = fxy[2+4*(i+noff+nxv*(j+moff))]; */
/*       }                                                        */
         for (i = 0; i < nn; i++) {
            v_at = _mm_loadu_ps(&fxy[4*(i+noff+nxv*(j+moff))]);
           _mm_storeu_ps(&sfxy[4*(i+mxv*j)],v_at);
         }
      }
      for (j = 0; j < mm; j++) {
/*       for (i = 0; i < nn; i++) {                               */
/*          sbxy[4*(i+mxv*j)] = bxy[4*(i+noff+nxv*(j+moff))];     */
/*          sbxy[1+4*(i+mxv*j)] = bxy[1+4*(i+noff+nxv*(j+moff))]; */
/*          sbxy[2+4*(i+mxv*j)] = bxy[2+4*(i+noff+nxv*(j+moff))]; */
/*       }                                                        */
         for (i = 0; i < nn; i++) {
            v_at = _mm_loadu_ps(&bxy[4*(i+noff+nxv*(j+moff))]);
           _mm_storeu_ps(&sbxy[4*(i+mxv*j)],v_at);
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
/* vector loop over particles in blocks of 4 */
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
/*       nm = 4*(nn - noff + mxv*(mm - moff)); */
         v_nn = _mm_sub_epi32(v_nn,v_noff);
         v_mm = _mm_sub_epi32(v_mm,v_moff);
         v_it = _mm_mul_epu32(v_mxv,_mm_srli_si128(v_mm,4));
         v_mm = _mm_mul_epu32(v_mm,v_mxv);
         v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
         v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
/*       amx = 1.0f - dxp; */
/*       amy = 1.0f - dyp; */
         v_amx = _mm_sub_ps(v_one,v_dxp);
         v_amy = _mm_sub_ps(v_one,v_dyp);
        _mm_store_si128((__m128i *)ll,v_nn);
/* find electric field */
/*       nn = nm;                          */
/*       dx = amx*sfxy[nn];                */
/*       dy = amx*sfxy[nn+1];              */
/*       dz = amx*sfxy[nn+2];              */
/*       mm = nn + 4;                      */
/*       dx = amy*(dxp*sfxy[mm] + dx);     */
/*       dy = amy*(dxp*sfxy[mm+1] + dy);   */
/*       dz = amy*(dxp*sfxy[mm+2] + dz);   */
/*       nn += 4*mxv;                      */
/*       acx = amx*sfxy[nn];               */
/*       acy = amx*sfxy[nn+1];             */
/*       acz = amx*sfxy[nn+2];             */
/*       mm = nn + 4;                      */
/*       dx += dyp*(dxp*sfxy[mm] + acx);   */
/*       dy += dyp*(dxp*sfxy[mm+1] + acy); */
/*       dz += dyp*(dxp*sfxy[mm+2] + acz); */
/* find magnetic field */
/*       nn = nm;                          */
/*       ox = amx*sbxy[nn];                */
/*       oy = amx*sbxy[nn+1];              */
/*       oz = amx*sbxy[nn+2];              */
/*       mm = nn + 4;                      */
/*       ox = amy*(dxp*sbxy[mm] + ox);     */
/*       oy = amy*(dxp*sbxy[mm+1] + oy);   */
/*       oz = amy*(dxp*sbxy[mm+2] + oz);   */
/*       nn += 4*mxv;                      */
/*       acx = amx*sbxy[nn];               */
/*       acy = amx*sbxy[nn+1];             */
/*       acz = amx*sbxy[nn+2];             */
/*       mm = nn + 4;                      */
/*       ox += dyp*(dxp*sbxy[mm] + acx);   */
/*       oy += dyp*(dxp*sbxy[mm+1] + acy); */
/*       oz += dyp*(dxp*sbxy[mm+2] + acz); */
/* interpolate electric and magnetic fields for first particle */
         nn = ll[0];
         v_at = _mm_shuffle_ps(v_amx,v_amx,0);
         a = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         e = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,0);
         nn += 4;
         a = _mm_add_ps(a,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         e = _mm_add_ps(e,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,0);
         a = _mm_mul_ps(a,v_at);
         e = _mm_mul_ps(e,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,0);
         a = _mm_add_ps(a,_mm_mul_ps(v_dx,v_at));
         e = _mm_add_ps(e,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for second particle */
         nn = ll[1];
         v_at = _mm_shuffle_ps(v_amx,v_amx,85);
         b = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         f = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,85);
         nn += 4;
         b = _mm_add_ps(b,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         f = _mm_add_ps(f,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,85);
         b = _mm_mul_ps(b,v_at);
         f = _mm_mul_ps(f,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,85);
         b = _mm_add_ps(b,_mm_mul_ps(v_dx,v_at));
         f = _mm_add_ps(f,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for third particle */
         nn = ll[2];
         v_at = _mm_shuffle_ps(v_amx,v_amx,170);
         c = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         g = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,170);
         nn += 4;
         c = _mm_add_ps(c,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         g = _mm_add_ps(g,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,170);
         c = _mm_mul_ps(c,v_at);
         g = _mm_mul_ps(g,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,170);
         c = _mm_add_ps(c,_mm_mul_ps(v_dx,v_at));
         g = _mm_add_ps(g,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for fourth particle */
         nn = ll[3];
         v_at = _mm_shuffle_ps(v_amx,v_amx,255);
         d = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         h = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,255);
         nn += 4;
         d = _mm_add_ps(d,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         h = _mm_add_ps(h,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
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
/*       dx *= qtmh; */
/*       dy *= qtmh; */
/*       dz *= qtmh; */
         v_dx = _mm_mul_ps(a,v_qtmh);
         v_dy = _mm_mul_ps(b,v_qtmh);
         v_dz = _mm_mul_ps(c,v_qtmh);
/* half acceleration */
/*       acx = ppart[j+2*nppmx+npoff] + dx; */
/*       acy = ppart[j+3*nppmx+npoff] + dy; */
/*       acz = ppart[j+4*nppmx+npoff] + dz; */
         a = _mm_add_ps(v_dx,_mm_load_ps(&ppart[j+2*nppmx+npoff]));
         b = _mm_add_ps(v_dy,_mm_load_ps(&ppart[j+3*nppmx+npoff]));
         c = _mm_add_ps(v_dz,_mm_load_ps(&ppart[j+4*nppmx+npoff]));
/* time-centered kinetic energy */
/*       sum1 += (acx*acx + acy*acy + acz*acz); */
         v_at = _mm_add_ps(_mm_mul_ps(a,a),_mm_mul_ps(b,b));
         v_at = _mm_add_ps(v_at,_mm_mul_ps(c,c));
/* convert to double precision before accumulating */
         v_d = _mm_cvtps_pd(v_at);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
         v_it = _mm_srli_si128((__m128i)v_at,8);
         v_d = _mm_cvtps_pd((__m128)v_it);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
/* calculate cyclotron frequency */
/*       omxt = qtmh*ox; */
/*       omyt = qtmh*oy; */
/*       omzt = qtmh*oz; */
         e = _mm_mul_ps(v_qtmh,e);
         f = _mm_mul_ps(v_qtmh,f);
         g = _mm_mul_ps(v_qtmh,g);
/* calculate rotation matrix */
/*       vx = omxt*omxt; */
         v_vx = _mm_mul_ps(e,e);
/*       vy = omyt*omyt; */
         v_vy = _mm_mul_ps(f,f);
/*       vz = omzt*omzt; */
         v_vz = _mm_mul_ps(g,g);
/*       omt = omxt*omxt + omyt*omyt + omzt*omzt; */
         v_at = _mm_add_ps(_mm_add_ps(v_vx,v_vy),v_vz);
/*       anorm = 2.0f/(1.0f + omt); */
         d = _mm_div_ps(v_two,_mm_add_ps(v_one,v_at));
/*       omt = 0.5f*(1.0f - omt); */
         h = _mm_mul_ps(v_half,_mm_sub_ps(v_one,v_at));
/*       vx = (omt + vx)*acx; */
         v_vx = _mm_mul_ps(_mm_add_ps(h,v_vx),a);
/*       vy = (omt + vy)*acy; */
         v_vy = _mm_mul_ps(_mm_add_ps(h,v_vy),b);
/*       vz = (omt + vz)*acz; */
         v_vz = _mm_mul_ps(_mm_add_ps(h,v_vz),c);
/*       omt = omxt*omyt; */
         h = _mm_mul_ps(e,f);
/*       vx = vx + (omzt + omt)*acy; */
         v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_add_ps(h,g),b));
/*       vy = vy + (omt - omzt)*acx; */
         v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_sub_ps(h,g),a));
/*       omt = omxt*omzt;  */
         h = _mm_mul_ps(e,g);
/*       vx = vx + (omt - omyt)*acz; */
         v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_sub_ps(h,f),c));
/*       vz = vz + (omt + omyt)*acx; */
         v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_add_ps(h,f),a));
/*       omt = omyt*omzt; */
         h = _mm_mul_ps(f,g);
/*       vy = vy + (omt + omxt)*acz; */
         v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_add_ps(h,e),c));
/*       vz = vz + (omt - omxt)*acy; */
         v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_sub_ps(h,e),b));
/* new velocity */
/*       vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx; */
/*       vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy; */
/*       vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz; */
         v_vx = _mm_add_ps(v_dx,_mm_mul_ps(v_vx,d));
         v_vy = _mm_add_ps(v_dy,_mm_mul_ps(v_vy,d));
         v_vz = _mm_add_ps(v_dz,_mm_mul_ps(v_vz,d));
/* new position */
/*       dx = x + vx*dtc; */
/*       dy = y + vy*dtc; */
         v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_dtc));
         v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_dtc));
/* find particles going out of bounds */
         mm = 0;
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
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy; */
         _mm_store_ps(&ppart[j+npoff],v_dx);
         _mm_store_ps(&ppart[j+nppmx+npoff],v_dy);
/* set new velocity */
/*       ppart[j+2*nppmx+npoff] = vx; */
/*       ppart[j+3*nppmx+npoff] = vy; */
/*       ppart[j+4*nppmx+npoff] = vz; */
         _mm_store_ps(&ppart[j+2*nppmx+npoff],v_vx);
         _mm_store_ps(&ppart[j+3*nppmx+npoff],v_vy);
         _mm_store_ps(&ppart[j+4*nppmx+npoff],v_vz);
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
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 4*(nn - noff + mxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + 4;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += 4*mxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + 4;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + 4;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += 4*mxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + 4;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
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
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* new position */
         dx = x + vx*dtc;
         dy = y + vy*dtc;
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
         ppart[j+4*nppmx+npoff] = vz;
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
   *ek += 0.5*sum2;
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void csse2grbppush23lt(float ppart[], float fxy[], float bxy[],
                       int kpic[], float qbm, float dt, float dtc, 
                       float ci, float *ek, int idimp, int nppmx,
                       int nx, int ny, int mx, int my, int nxv, int nyv,
                       int mx1, int mxy1, int ipbc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   vector/OpenMP version using guard cells
   particles stored in segmented array
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic = number of particles per tile
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   idimp = size of phase space = 5
   nppmx = maximum number of particles in tile
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of field arrays, must be >= nx+1
   nyv = second dimension of field arrays, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   requires SSE2, ppart, fxy, and bxy need to be 16 byte aligned
   nppmx needs to be a multiple of 4, fxy, bxy need to have 4 components
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nps, nn, mm, nm;
   float qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg;
   float omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, vx, vy, vz;
   double sum1, sum2;
   __m128i v_noff, v_moff, v_mxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qtmh, v_ci2, v_dt, v_dtc, v_one, v_two, v_half;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_gami, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m128 v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 a, b, c, d, e, f, g, h;
   __m128d v_sum1, v_d;
   __attribute__((aligned(16))) unsigned int ll[4];
   __attribute__((aligned(16))) double dd[2];
   __attribute__((aligned(16))) float sfxy[4*MXV*MYV], sbxy[4*MXV*MYV];
/* __attribute__((aligned(16))) float sfxy[4*(mx+1)*(my+1)]; */
/* __attribute__((aligned(16))) float sbxy[4*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
   ci2 = ci*ci;
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
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nps,nn,mm,nm,x,y,vx,vy,vz,dxp,dyp, \
amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1, \
rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,dtg,sum1,v_noff, \
v_moff,v_nn,v_mm,v_it,v_x,v_y,v_vx,v_vy,v_vz,v_dxp,v_dyp,v_amx,v_amy, \
v_dx,v_dy,v_dz,v_gami,v_at,v_d,v_sum1,a,b,c,d,e,f,g,h,ll,dd,sfxy,sbxy) \
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
      for (j = 0; j < mm; j++) {
/*       for (i = 0; i < nn; i++) {                               */
/*          sfxy[4*(i+mxv*j)] = fxy[4*(i+noff+nxv*(j+moff))];     */
/*          sfxy[1+4*(i+mxv*j)] = fxy[1+4*(i+noff+nxv*(j+moff))]; */
/*          sfxy[2+4*(i+mxv*j)] = fxy[2+4*(i+noff+nxv*(j+moff))]; */
/*       }                                                        */
         for (i = 0; i < nn; i++) {
            v_at = _mm_loadu_ps(&fxy[4*(i+noff+nxv*(j+moff))]);
           _mm_storeu_ps(&sfxy[4*(i+mxv*j)],v_at);
         }
      }
      for (j = 0; j < mm; j++) {
/*       for (i = 0; i < nn; i++) {                               */
/*          sbxy[4*(i+mxv*j)] = bxy[4*(i+noff+nxv*(j+moff))];     */
/*          sbxy[1+4*(i+mxv*j)] = bxy[1+4*(i+noff+nxv*(j+moff))]; */
/*          sbxy[2+4*(i+mxv*j)] = bxy[2+4*(i+noff+nxv*(j+moff))]; */
/*       }                                                        */
         for (i = 0; i < nn; i++) {
            v_at = _mm_loadu_ps(&bxy[4*(i+noff+nxv*(j+moff))]);
           _mm_storeu_ps(&sbxy[4*(i+mxv*j)],v_at);
         }
      }
      nps = 4*(npp/4);
      sum1 = 0.0;
      v_sum1 = _mm_set1_pd(0.0);
/* vector loop over particles in blocks of 4 */
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
/*       nm = 4*(nn - noff + mxv*(mm - moff)); */
         v_nn = _mm_sub_epi32(v_nn,v_noff);
         v_mm = _mm_sub_epi32(v_mm,v_moff);
         v_it = _mm_mul_epu32(v_mxv,_mm_srli_si128(v_mm,4));
         v_mm = _mm_mul_epu32(v_mm,v_mxv);
         v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
         v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
/*       amx = 1.0f - dxp; */
/*       amy = 1.0f - dyp; */
         v_amx = _mm_sub_ps(v_one,v_dxp);
         v_amy = _mm_sub_ps(v_one,v_dyp);
        _mm_store_si128((__m128i *)ll,v_nn);
/* find electric field */
/*       nn = nm;                          */
/*       dx = amx*sfxy[nn];                */
/*       dy = amx*sfxy[nn+1];              */
/*       dz = amx*sfxy[nn+2];              */
/*       mm = nn + 4;                      */
/*       dx = amy*(dxp*sfxy[mm] + dx);     */
/*       dy = amy*(dxp*sfxy[mm+1] + dy);   */
/*       dz = amy*(dxp*sfxy[mm+2] + dz);   */
/*       nn += 4*mxv;                      */
/*       acx = amx*sfxy[nn];               */
/*       acy = amx*sfxy[nn+1];             */
/*       acz = amx*sfxy[nn+2];             */
/*       mm = nn + 4;                      */
/*       dx += dyp*(dxp*sfxy[mm] + acx);   */
/*       dy += dyp*(dxp*sfxy[mm+1] + acy); */
/*       dz += dyp*(dxp*sfxy[mm+2] + acz); */
/* find magnetic field */
/*       nn = nm;                          */
/*       ox = amx*sbxy[nn];                */
/*       oy = amx*sbxy[nn+1];              */
/*       oz = amx*sbxy[nn+2];              */
/*       mm = nn + 4;                      */
/*       ox = amy*(dxp*sbxy[mm] + ox);     */
/*       oy = amy*(dxp*sbxy[mm+1] + oy);   */
/*       oz = amy*(dxp*sbxy[mm+2] + oz);   */
/*       nn += 4*mxv;                      */
/*       acx = amx*sbxy[nn];               */
/*       acy = amx*sbxy[nn+1];             */
/*       acz = amx*sbxy[nn+2];             */
/*       mm = nn + 4;                      */
/*       ox += dyp*(dxp*sbxy[mm] + acx);   */
/*       oy += dyp*(dxp*sbxy[mm+1] + acy); */
/*       oz += dyp*(dxp*sbxy[mm+2] + acz); */
/* interpolate electric and magnetic fields for first particle */
         nn = ll[0];
         v_at = _mm_shuffle_ps(v_amx,v_amx,0);
         a = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         e = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,0);
         nn += 4;
         a = _mm_add_ps(a,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         e = _mm_add_ps(e,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,0);
         a = _mm_mul_ps(a,v_at);
         e = _mm_mul_ps(e,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,0);
         a = _mm_add_ps(a,_mm_mul_ps(v_dx,v_at));
         e = _mm_add_ps(e,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for second particle */
         nn = ll[1];
         v_at = _mm_shuffle_ps(v_amx,v_amx,85);
         b = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         f = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,85);
         nn += 4;
         b = _mm_add_ps(b,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         f = _mm_add_ps(f,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,85);
         b = _mm_mul_ps(b,v_at);
         f = _mm_mul_ps(f,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,85);
         b = _mm_add_ps(b,_mm_mul_ps(v_dx,v_at));
         f = _mm_add_ps(f,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for third particle */
         nn = ll[2];
         v_at = _mm_shuffle_ps(v_amx,v_amx,170);
         c = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         g = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,170);
         nn += 4;
         c = _mm_add_ps(c,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         g = _mm_add_ps(g,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,170);
         c = _mm_mul_ps(c,v_at);
         g = _mm_mul_ps(g,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,170);
         c = _mm_add_ps(c,_mm_mul_ps(v_dx,v_at));
         g = _mm_add_ps(g,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for fourth particle */
         nn = ll[3];
         v_at = _mm_shuffle_ps(v_amx,v_amx,255);
         d = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         h = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,255);
         nn += 4;
         d = _mm_add_ps(d,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         h = _mm_add_ps(h,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
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
/*       dx *= qtmh; */
/*       dy *= qtmh; */
/*       dz *= qtmh; */
         v_dx = _mm_mul_ps(a,v_qtmh);
         v_dy = _mm_mul_ps(b,v_qtmh);
         v_dz = _mm_mul_ps(c,v_qtmh);
/* half acceleration */
/*       acx = ppart[j+2*nppmx+npoff] + dx; */
/*       acy = ppart[j+3*nppmx+npoff] + dy; */
/*       acz = ppart[j+4*nppmx+npoff] + dz; */
         a = _mm_add_ps(v_dx,_mm_load_ps(&ppart[j+2*nppmx+npoff]));
         b = _mm_add_ps(v_dy,_mm_load_ps(&ppart[j+3*nppmx+npoff]));
         c = _mm_add_ps(v_dz,_mm_load_ps(&ppart[j+4*nppmx+npoff]));
/* find inverse gamma */
/*       p2 = acx*acx + acy*acy + acz*acz; */
         v_at = _mm_add_ps(_mm_mul_ps(a,a),_mm_mul_ps(b,b));
         v_at = _mm_add_ps(v_at,_mm_mul_ps(c,c));
/*       gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_gami = _mm_rsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* full accuracy calculation */
         v_gami = _mm_sqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2)));
         v_gami = _mm_div_ps(v_one,v_gami);
/* full accuracy calculation with SVML */
/*       v_gami = _mm_invsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* time-centered kinetic energy */
/*       sum1 += gami*p2/(1.0f + gami); */
         v_at = _mm_mul_ps(v_gami,v_at);
         v_at = _mm_div_ps(v_at,_mm_add_ps(v_one,v_gami));
/* convert to double precision before accumulating */
         v_d = _mm_cvtps_pd(v_at);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
         v_it = _mm_srli_si128((__m128i)v_at,8);
         v_d = _mm_cvtps_pd((__m128)v_it);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
/* renormalize magnetic field */
/*       qtmg = qtmh*gami; */
         v_at = _mm_mul_ps(v_qtmh,v_gami);
/* calculate cyclotron frequency */
/*       omxt = qtmg*ox; */
/*       omyt = qtmg*oy; */
/*       omzt = qtmg*oz; */
         e = _mm_mul_ps(v_at,e);
         f = _mm_mul_ps(v_at,f);
         g = _mm_mul_ps(v_at,g);
/* calculate rotation matrix */
/*       vx = omxt*omxt; */
         v_vx = _mm_mul_ps(e,e);
/*       vy = omyt*omyt; */
         v_vy = _mm_mul_ps(f,f);
/*       vz = omzt*omzt; */
         v_vz = _mm_mul_ps(g,g);
/*       omt = omxt*omxt + omyt*omyt + omzt*omzt; */
         v_at = _mm_add_ps(_mm_add_ps(v_vx,v_vy),v_vz);
/*       anorm = 2.0f/(1.0f + omt); */
         d = _mm_div_ps(v_two,_mm_add_ps(v_one,v_at));
/*       omt = 0.5f*(1.0f - omt); */
         h = _mm_mul_ps(v_half,_mm_sub_ps(v_one,v_at));
/*       vx = (omt + vx)*acx; */
         v_vx = _mm_mul_ps(_mm_add_ps(h,v_vx),a);
/*       vy = (omt + vy)*acy; */
         v_vy = _mm_mul_ps(_mm_add_ps(h,v_vy),b);
/*       vz = (omt + vz)*acz; */
         v_vz = _mm_mul_ps(_mm_add_ps(h,v_vz),c);
/*       omt = omxt*omyt; */
         h = _mm_mul_ps(e,f);
/*       vx = vx + (omzt + omt)*acy; */
         v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_add_ps(h,g),b));
/*       vy = vy + (omt - omzt)*acx; */
         v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_sub_ps(h,g),a));
/*       omt = omxt*omzt;  */
         h = _mm_mul_ps(e,g);
/*       vx = vx + (omt - omyt)*acz; */
         v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_sub_ps(h,f),c));
/*       vz = vz + (omt + omyt)*acx; */
         v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_add_ps(h,f),a));
/*       omt = omyt*omzt; */
         h = _mm_mul_ps(f,g);
/*       vy = vy + (omt + omxt)*acz; */
         v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_add_ps(h,e),c));
/*       vz = vz + (omt - omxt)*acy; */
         v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_sub_ps(h,e),b));
/* new momentum */
/*       vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm; */
/*       vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm; */
/*       vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm; */
         v_vx = _mm_add_ps(v_dx,_mm_mul_ps(v_vx,d));
         v_vy = _mm_add_ps(v_dy,_mm_mul_ps(v_vy,d));
         v_vz = _mm_add_ps(v_dz,_mm_mul_ps(v_vz,d));
/* update inverse gamma */
/*       p2 = vx*vx + vy*vy + vz*vz; */
         v_at = _mm_mul_ps(v_vx,v_vx);
         v_at = _mm_add_ps(v_at,_mm_mul_ps(v_vy,v_vy));
         v_at = _mm_add_ps(v_at,_mm_mul_ps(v_vz,v_vz));
/*       dtg = dtc/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_at = _mm_rsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/*       v_at = _mm_mul_ps(v_dtc,v_at); */
/* full accuracy calculation */
         v_at = _mm_sqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2)));
         v_at = _mm_div_ps(v_dtc,v_at);
/* full accuracy calculation with SVML */
/*       v_at = _mm_invsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/*       v_at = _mm_mul_ps(v_dtc,v_at); */
/* new position */
/*       dx = x + vx*dtg; */
/*       dy = y + vy*dtg; */
         v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_at));
         v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_at));
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
/* set new momentum */
/*       ppart[j+2*nppmx+npoff] = vx; */
/*       ppart[j+3*nppmx+npoff] = vy; */
/*       ppart[j+4*nppmx+npoff] = vz; */
         _mm_store_ps(&ppart[j+2*nppmx+npoff],v_vx);
         _mm_store_ps(&ppart[j+3*nppmx+npoff],v_vy);
         _mm_store_ps(&ppart[j+4*nppmx+npoff],v_vz);
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 4*(nn - noff + mxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + 4;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += 4*mxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + 4;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + 4;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += 4*mxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + 4;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
/* find inverse gamma */
         p2 = acx*acx + acy*acy + acz*acz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
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
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* update inverse gamma */
         p2 = vx*vx + vy*vy + vz*vz;
         dtg = dtc/sqrtf(1.0f + p2*ci2);
/* new position */
         dx = x + vx*dtg;
         dy = y + vy*dtg;
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
/* set new momentum */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
         ppart[j+4*nppmx+npoff] = vz;
      }
/*    sum2 += sum1; */
      _mm_store_pd(&dd[0],v_sum1);
      for (j = 1; j < 2; j++) {
         dd[0] += dd[j];
      }
      sum2 += (sum1 + dd[0]);
   }
/* normalize kinetic energy */
   *ek += sum2;
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void csse2grbppushf23lt(float ppart[], float fxy[], float bxy[],
                        int kpic[], int ncl[], int ihole[], float qbm,
                        float dt, float dtc, float ci, float *ek,
                        int idimp, int nppmx, int nx, int ny, int mx,
                        int my, int nxv, int nyv, int mx1, int mxy1,
                        int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   with periodic boundary conditions.
   Using the Boris Mover.
   also determines list of particles which are leaving this tile
   vector/OpenMP version using guard cells
   particles stored in segmented array
   131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   fxy[k][j][2] = z component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   bxy[k][j][0] = x component of magnetic field at grid (j,k)
   bxy[k][j][1] = y component of magnetic field at grid (j,k)
   bxy[k][j][2] = z component of magnetic field at grid (j,k)
   that is, the convolution of magnetic field over particle shape
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   idimp = size of phase space = 5
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
   requires SSE2, ppart, fxy, and bxy need to be 16 byte aligned
   nppmx needs to be a multiple of 4, fxy, bxy need to have 4 components
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nps, ih, nh, nn, mm, nm, kk;
   float qtmh, ci2, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz;
   float acx, acy, acz, p2, gami, qtmg, dtg, omxt, omyt, omzt, omt;
   float anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float anx, any, edgelx, edgely, edgerx, edgery;
   float x, y, vx, vy, vz;
   double sum1, sum2;
   __m128i v_noff, v_moff, v_mxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qtmh, v_ci2, v_dt, v_dtc, v_one, v_two, v_half;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_gami, v_st, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m128 v_anx, v_any, v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 v_zero, v_three, v_six;
   __m128 a, b, c, d, e, f, g, h;
   __m128d v_sum1, v_d;
   __attribute__((aligned(16))) unsigned int ll[4], lm[8];
   __attribute__((aligned(16))) unsigned long jj[1];
   __attribute__((aligned(16))) double dd[2];
   __attribute__((aligned(16))) float sfxy[4*MXV*MYV], sbxy[4*MXV*MYV];
/* __attribute__((aligned(16))) float sfxy[4*(mx+1)*(my+1)]; */
/* __attribute__((aligned(16))) float sbxy[4*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   qtmh = 0.5f*qbm*dt;
   ci2 = ci*ci;
   anx = (float) nx;
   any = (float) ny;
   sum2 = 0.0;
   v_mxv = _mm_set1_epi32(mxv);
   v_qtmh = _mm_set1_ps(qtmh);
   v_ci2 = _mm_set1_ps(ci2);
   v_dt = _mm_set1_ps(dt);
   v_dtc = _mm_set1_ps(dtc);
   v_anx = _mm_set1_ps(anx);
   v_any = _mm_set1_ps(any);
   v_zero = _mm_setzero_ps();
   v_one = _mm_set1_ps(1.0f);
   v_two = _mm_set1_ps(2.0f);
   v_half = _mm_set1_ps(0.5f);
   v_three = _mm_set1_ps(3.0f);
   v_six = _mm_set1_ps(6.0f);
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV)) */
/*    return;                      */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nps,nn,mm,nm,kk,ih,nh,x,y,vx,vy,vz, \
dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm, \
rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx, \
edgery,p2,gami,qtmg,dtg,sum1,v_noff,v_moff,v_nn,v_mm,v_it,v_x,v_y,v_vx, \
v_vy,v_vz,v_dxp,v_dyp,v_amx,v_amy,v_dx,v_dy,v_dz,v_gami,v_st,v_at, \
v_edgelx,v_edgely,v_edgerx,v_edgery,v_d,v_sum1,a,b,c,d,e,f,g,h,jj,ll, \
lm,dd,sfxy,sbxy) \
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
      for (j = 0; j < mm; j++) {
/*       for (i = 0; i < nn; i++) {                               */
/*          sfxy[4*(i+mxv*j)] = fxy[4*(i+noff+nxv*(j+moff))];     */
/*          sfxy[1+4*(i+mxv*j)] = fxy[1+4*(i+noff+nxv*(j+moff))]; */
/*          sfxy[2+4*(i+mxv*j)] = fxy[2+4*(i+noff+nxv*(j+moff))]; */
/*       }                                                        */
         for (i = 0; i < nn; i++) {
            v_at = _mm_loadu_ps(&fxy[4*(i+noff+nxv*(j+moff))]);
           _mm_storeu_ps(&sfxy[4*(i+mxv*j)],v_at);
         }
      }
      for (j = 0; j < mm; j++) {
/*       for (i = 0; i < nn; i++) {                               */
/*          sbxy[4*(i+mxv*j)] = bxy[4*(i+noff+nxv*(j+moff))];     */
/*          sbxy[1+4*(i+mxv*j)] = bxy[1+4*(i+noff+nxv*(j+moff))]; */
/*          sbxy[2+4*(i+mxv*j)] = bxy[2+4*(i+noff+nxv*(j+moff))]; */
/*       }                                                        */
         for (i = 0; i < nn; i++) {
            v_at = _mm_loadu_ps(&bxy[4*(i+noff+nxv*(j+moff))]);
           _mm_storeu_ps(&sbxy[4*(i+mxv*j)],v_at);
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
/* vector loop over particles in blocks of 4 */
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
/*       nm = 4*(nn - noff + mxv*(mm - moff)); */
         v_nn = _mm_sub_epi32(v_nn,v_noff);
         v_mm = _mm_sub_epi32(v_mm,v_moff);
         v_it = _mm_mul_epu32(v_mxv,_mm_srli_si128(v_mm,4));
         v_mm = _mm_mul_epu32(v_mm,v_mxv);
         v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
         v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
/*       amx = 1.0f - dxp; */
/*       amy = 1.0f - dyp; */
         v_amx = _mm_sub_ps(v_one,v_dxp);
         v_amy = _mm_sub_ps(v_one,v_dyp);
        _mm_store_si128((__m128i *)ll,v_nn);
/* find electric field */
/*       nn = nm;                          */
/*       dx = amx*sfxy[nn];                */
/*       dy = amx*sfxy[nn+1];              */
/*       dz = amx*sfxy[nn+2];              */
/*       mm = nn + 4;                      */
/*       dx = amy*(dxp*sfxy[mm] + dx);     */
/*       dy = amy*(dxp*sfxy[mm+1] + dy);   */
/*       dz = amy*(dxp*sfxy[mm+2] + dz);   */
/*       nn += 4*mxv;                      */
/*       acx = amx*sfxy[nn];               */
/*       acy = amx*sfxy[nn+1];             */
/*       acz = amx*sfxy[nn+2];             */
/*       mm = nn + 4;                      */
/*       dx += dyp*(dxp*sfxy[mm] + acx);   */
/*       dy += dyp*(dxp*sfxy[mm+1] + acy); */
/*       dz += dyp*(dxp*sfxy[mm+2] + acz); */
/* find magnetic field */
/*       nn = nm;                          */
/*       ox = amx*sbxy[nn];                */
/*       oy = amx*sbxy[nn+1];              */
/*       oz = amx*sbxy[nn+2];              */
/*       mm = nn + 4;                      */
/*       ox = amy*(dxp*sbxy[mm] + ox);     */
/*       oy = amy*(dxp*sbxy[mm+1] + oy);   */
/*       oz = amy*(dxp*sbxy[mm+2] + oz);   */
/*       nn += 4*mxv;                      */
/*       acx = amx*sbxy[nn];               */
/*       acy = amx*sbxy[nn+1];             */
/*       acz = amx*sbxy[nn+2];             */
/*       mm = nn + 4;                      */
/*       ox += dyp*(dxp*sbxy[mm] + acx);   */
/*       oy += dyp*(dxp*sbxy[mm+1] + acy); */
/*       oz += dyp*(dxp*sbxy[mm+2] + acz); */
/* interpolate electric and magnetic fields for first particle */
         nn = ll[0];
         v_at = _mm_shuffle_ps(v_amx,v_amx,0);
         a = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         e = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,0);
         nn += 4;
         a = _mm_add_ps(a,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         e = _mm_add_ps(e,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,0);
         a = _mm_mul_ps(a,v_at);
         e = _mm_mul_ps(e,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,0);
         a = _mm_add_ps(a,_mm_mul_ps(v_dx,v_at));
         e = _mm_add_ps(e,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for second particle */
         nn = ll[1];
         v_at = _mm_shuffle_ps(v_amx,v_amx,85);
         b = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         f = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,85);
         nn += 4;
         b = _mm_add_ps(b,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         f = _mm_add_ps(f,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,85);
         b = _mm_mul_ps(b,v_at);
         f = _mm_mul_ps(f,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,85);
         b = _mm_add_ps(b,_mm_mul_ps(v_dx,v_at));
         f = _mm_add_ps(f,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for third particle */
         nn = ll[2];
         v_at = _mm_shuffle_ps(v_amx,v_amx,170);
         c = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         g = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,170);
         nn += 4;
         c = _mm_add_ps(c,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         g = _mm_add_ps(g,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
         v_at = _mm_shuffle_ps(v_amy,v_amy,170);
         c = _mm_mul_ps(c,v_at);
         g = _mm_mul_ps(g,v_at);
         v_at = _mm_shuffle_ps(v_dyp,v_dyp,170);
         c = _mm_add_ps(c,_mm_mul_ps(v_dx,v_at));
         g = _mm_add_ps(g,_mm_mul_ps(v_dy,v_at));
/* interpolate electric and magnetic fields for fourth particle */
         nn = ll[3];
         v_at = _mm_shuffle_ps(v_amx,v_amx,255);
         d = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn]));
         h = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn]));
         mm = nn + 4*mxv;
         v_dx = _mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm]));
         v_dy = _mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm]));
         v_at = _mm_shuffle_ps(v_dxp,v_dxp,255);
         nn += 4;
         d = _mm_add_ps(d,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[nn])));
         h = _mm_add_ps(h,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[nn])));
         mm += 4;
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(v_at,_mm_load_ps(&sfxy[mm])));
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(v_at,_mm_load_ps(&sbxy[mm])));
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
/*       dx *= qtmh; */
/*       dy *= qtmh; */
/*       dz *= qtmh; */
         v_dx = _mm_mul_ps(a,v_qtmh);
         v_dy = _mm_mul_ps(b,v_qtmh);
         v_dz = _mm_mul_ps(c,v_qtmh);
/* half acceleration */
/*       acx = ppart[j+2*nppmx+npoff] + dx; */
/*       acy = ppart[j+3*nppmx+npoff] + dy; */
/*       acz = ppart[j+4*nppmx+npoff] + dz; */
         a = _mm_add_ps(v_dx,_mm_load_ps(&ppart[j+2*nppmx+npoff]));
         b = _mm_add_ps(v_dy,_mm_load_ps(&ppart[j+3*nppmx+npoff]));
         c = _mm_add_ps(v_dz,_mm_load_ps(&ppart[j+4*nppmx+npoff]));
/* find inverse gamma */
/*       p2 = acx*acx + acy*acy + acz*acz; */
         v_at = _mm_add_ps(_mm_mul_ps(a,a),_mm_mul_ps(b,b));
         v_at = _mm_add_ps(v_at,_mm_mul_ps(c,c));
/*       gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_gami = _mm_rsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* full accuracy calculation */
         v_gami = _mm_sqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2)));
         v_gami = _mm_div_ps(v_one,v_gami);
/* full accuracy calculation with SVML */
/*       v_gami = _mm_invsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* time-centered kinetic energy */
/*       sum1 += gami*p2/(1.0f + gami); */
         v_at = _mm_mul_ps(v_gami,v_at);
         v_at = _mm_div_ps(v_at,_mm_add_ps(v_one,v_gami));
/* convert to double precision before accumulating */
         v_d = _mm_cvtps_pd(v_at);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
         v_it = _mm_srli_si128((__m128i)v_at,8);
         v_d = _mm_cvtps_pd((__m128)v_it);
         v_sum1 = _mm_add_pd(v_sum1,v_d);
/* renormalize magnetic field */
/*       qtmg = qtmh*gami; */
         v_at = _mm_mul_ps(v_qtmh,v_gami);
/* calculate cyclotron frequency */
/*       omxt = qtmg*ox; */
/*       omyt = qtmg*oy; */
/*       omzt = qtmg*oz; */
         e = _mm_mul_ps(v_at,e);
         f = _mm_mul_ps(v_at,f);
         g = _mm_mul_ps(v_at,g);
/* calculate rotation matrix */
/*       vx = omxt*omxt; */
         v_vx = _mm_mul_ps(e,e);
/*       vy = omyt*omyt; */
         v_vy = _mm_mul_ps(f,f);
/*       vz = omzt*omzt; */
         v_vz = _mm_mul_ps(g,g);
/*       omt = omxt*omxt + omyt*omyt + omzt*omzt; */
         v_at = _mm_add_ps(_mm_add_ps(v_vx,v_vy),v_vz);
/*       anorm = 2.0f/(1.0f + omt); */
         d = _mm_div_ps(v_two,_mm_add_ps(v_one,v_at));
/*       omt = 0.5f*(1.0f - omt); */
         h = _mm_mul_ps(v_half,_mm_sub_ps(v_one,v_at));
/*       vx = (omt + vx)*acx; */
         v_vx = _mm_mul_ps(_mm_add_ps(h,v_vx),a);
/*       vy = (omt + vy)*acy; */
         v_vy = _mm_mul_ps(_mm_add_ps(h,v_vy),b);
/*       vz = (omt + vz)*acz; */
         v_vz = _mm_mul_ps(_mm_add_ps(h,v_vz),c);
/*       omt = omxt*omyt; */
         h = _mm_mul_ps(e,f);
/*       vx = vx + (omzt + omt)*acy; */
         v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_add_ps(h,g),b));
/*       vy = vy + (omt - omzt)*acx; */
         v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_sub_ps(h,g),a));
/*       omt = omxt*omzt;  */
         h = _mm_mul_ps(e,g);
/*       vx = vx + (omt - omyt)*acz; */
         v_vx = _mm_add_ps(v_vx,_mm_mul_ps(_mm_sub_ps(h,f),c));
/*       vz = vz + (omt + omyt)*acx; */
         v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_add_ps(h,f),a));
/*       omt = omyt*omzt; */
         h = _mm_mul_ps(f,g);
/*       vy = vy + (omt + omxt)*acz; */
         v_vy = _mm_add_ps(v_vy,_mm_mul_ps(_mm_add_ps(h,e),c));
/*       vz = vz + (omt - omxt)*acy; */
         v_vz = _mm_add_ps(v_vz,_mm_mul_ps(_mm_sub_ps(h,e),b));
/* new momentum */
/*       vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm; */
/*       vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm; */
/*       vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm; */
         v_vx = _mm_add_ps(v_dx,_mm_mul_ps(v_vx,d));
         v_vy = _mm_add_ps(v_dy,_mm_mul_ps(v_vy,d));
         v_vz = _mm_add_ps(v_dz,_mm_mul_ps(v_vz,d));
/* update inverse gamma */
/*       p2 = vx*vx + vy*vy + vz*vz; */
         v_at = _mm_mul_ps(v_vx,v_vx);
         v_at = _mm_add_ps(v_at,_mm_mul_ps(v_vy,v_vy));
         v_at = _mm_add_ps(v_at,_mm_mul_ps(v_vz,v_vz));
/*       dtg = dtc/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_at = _mm_rsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/*       v_at = _mm_mul_ps(v_dtc,v_at); */
/* full accuracy calculation */
         v_at = _mm_sqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2)));
         v_at = _mm_div_ps(v_dtc,v_at);
/* full accuracy calculation with SVML */
/*       v_at = _mm_invsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/*       v_at = _mm_mul_ps(v_dtc,v_at); */
/* new position */
/*       dx = x + vx*dtg; */
/*       dy = y + vy*dtg; */
         v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_at));
         v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_at));
/* find particles going out of bounds */
         mm = 0;
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
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy; */
         _mm_store_ps(&ppart[j+npoff],v_dx);
         _mm_store_ps(&ppart[j+nppmx+npoff],v_dy);
/* set new momentum */
/*       ppart[j+2*nppmx+npoff] = vx; */
/*       ppart[j+3*nppmx+npoff] = vy; */
/*       ppart[j+4*nppmx+npoff] = vz; */
         _mm_store_ps(&ppart[j+2*nppmx+npoff],v_vx);
         _mm_store_ps(&ppart[j+3*nppmx+npoff],v_vy);
         _mm_store_ps(&ppart[j+4*nppmx+npoff],v_vz);
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
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         nm = 4*(nn - noff + mxv*(mm - moff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
/* find electric field */
         nn = nm;
         dx = amx*sfxy[nn];
         dy = amx*sfxy[nn+1];
         dz = amx*sfxy[nn+2];
         mm = nn + 4;
         dx = amy*(dxp*sfxy[mm] + dx);
         dy = amy*(dxp*sfxy[mm+1] + dy);
         dz = amy*(dxp*sfxy[mm+2] + dz);
         nn += 4*mxv;
         acx = amx*sfxy[nn];
         acy = amx*sfxy[nn+1];
         acz = amx*sfxy[nn+2];
         mm = nn + 4;
         dx += dyp*(dxp*sfxy[mm] + acx);
         dy += dyp*(dxp*sfxy[mm+1] + acy);
         dz += dyp*(dxp*sfxy[mm+2] + acz);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxy[nn];
         oy = amx*sbxy[nn+1];
         oz = amx*sbxy[nn+2];
         mm = nn + 4;
         ox = amy*(dxp*sbxy[mm] + ox);
         oy = amy*(dxp*sbxy[mm+1] + oy);
         oz = amy*(dxp*sbxy[mm+2] + oz);
         nn += 4*mxv;
         acx = amx*sbxy[nn];
         acy = amx*sbxy[nn+1];
         acz = amx*sbxy[nn+2];
         mm = nn + 4;
         ox += dyp*(dxp*sbxy[mm] + acx);
         oy += dyp*(dxp*sbxy[mm+1] + acy);
         oz += dyp*(dxp*sbxy[mm+2] + acz);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+2*nppmx+npoff] + dx;
         acy = ppart[j+3*nppmx+npoff] + dy;
         acz = ppart[j+4*nppmx+npoff] + dz;
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
         vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx;
         vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy;
         vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz;
/* update inverse gamma */
         p2 = vx*vx + vy*vy + vz*vz;
         dtg = dtc/sqrtf(1.0f + p2*ci2);
/* new position */
         dx = x + vx*dtg;
         dy = y + vy*dtg;
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
/* set new momentum */
         ppart[j+2*nppmx+npoff] = vx;
         ppart[j+3*nppmx+npoff] = vy;
         ppart[j+4*nppmx+npoff] = vz;
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
   *ek += sum2;
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
void csse2gjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                    float dt, int nppmx, int idimp, int nx, int ny,
                    int mx, int my, int nxv, int nyv, int mx1, int mxy1,
                    int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   vector/OpenMP version using guard cells
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of current array, must be >= nx+1
   nyv = second dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   requires SSE2, ppart and cu need to be 16 byte aligned
   nppmx needs to be a multiple of 4, cu needs to have 4 components
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nps, nn, mm;
   float edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz;
   __m128i v_noff, v_moff, v_mxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qm, v_dt, v_one;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_vx, v_vy;
   __m128 v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 a, b, c, d, va, vb, vc, vd;
   __attribute__((aligned(16))) unsigned int ll[4];
   __attribute__((aligned(16))) unsigned long kk[1];
   __attribute__((aligned(16))) float scu[4*MXV*MYV];
/* __attribute__((aligned(16))) float scu[4*(mx+1)*(my+1)]; */
   mxv = mx + 1;
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
   v_qm = _mm_set1_ps(qm);
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
private(i,j,k,noff,moff,npp,npoff,nps,nn,mm,x,y,dxp,dyp,amx,amy,dx,dy, \
vx,vy,vz,v_noff,v_moff,v_nn,v_mm,v_it,v_x,v_y,v_vx,v_vy,v_dxp,v_dyp, \
v_amx,v_amy,v_dx,v_dy,v_at,a,b,c,d,va,vb,vc,vd,ll,kk,scu)
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
/*    for (j = 0; j < 4*mxv*(my+1); j++) { */
/*       scu[j] = 0.0f;                    */
/*    }                                    */
      memset((void*)scu,0,4*mxv*(my+1)*sizeof(float));
/* vector loop over particles in blocks of 4 */
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
/*       dyp = y - (float) mm;      */
         v_dyp = _mm_sub_ps(v_y,_mm_cvtepi32_ps(v_mm));
/*       nm = 4*(nn - noff + mxv*(mm - moff)); */
         v_nn = _mm_sub_epi32(v_nn,v_noff);
         v_mm = _mm_sub_epi32(v_mm,v_moff);
         v_it = _mm_mul_epu32(v_mxv,_mm_srli_si128(v_mm,4));
         v_mm = _mm_mul_epu32(v_mm,v_mxv);
         v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
         v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
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
/* deposit current */
/*       vx = ppart[j+2*nppmx+npoff]; */
/*       vy = ppart[j+3*nppmx+npoff]; */
/*       vz = ppart[j+4*nppmx+npoff]; */
         v_vx = _mm_load_ps(&ppart[j+2*nppmx+npoff]);
         v_vy = _mm_load_ps(&ppart[j+3*nppmx+npoff]);
         va = v_vx;
         vb = v_vy;
         vc = _mm_load_ps(&ppart[j+4*nppmx+npoff]);
         vd = _mm_setzero_ps();
/* transpose so va,vb,vc,vd contain the 3 velocities plus zero */
/* for each of 4 particles                                     */
         _MM_TRANSPOSE4_PS(va,vb,vc,vd);
/*       dx = amx*amy;       */
/*       scu[nn] += vx*dx;   */
/*       scu[nn+1] += vy*dx; */
/*       scu[nn+2] += vz*dx; */
/*       dy = dxp*amy;       */
/*       mm = nn + 4;        */
/*       scu[mm] += vx*dy;   */
/*       scu[mm+1] += vy*dy; */
/*       scu[mm+2] += vz*dy; */
/*       dx = amx*dyp;       */
/*       nn += 4*nxv;        */
/*       scu[nn] += vx*dx;   */
/*       scu[nn+1] += vy*dx; */
/*       scu[nn+2] += vz*dx; */
/*       dy = dxp*dyp;       */
/*       mm = nn + 4;        */
/*       scu[mm] += vx*dy;   */
/*       scu[mm+1] += vy*dy; */
/*       scu[mm+2] += vz*dy; */
/* deposit for first particle */
         mm = ll[0];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(a,a,0)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(b,b,0)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(c,c,0)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(d,d,0)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for second particle */
         mm = ll[1];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(a,a,85)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(b,b,85)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(c,c,85)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(d,d,85)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for third particle */
         mm = ll[2];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(a,a,170)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(b,b,170)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(c,c,170)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(d,d,170)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for fourth particle */
         mm = ll[3];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(a,a,255)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(b,b,255)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(c,c,255)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(d,d,255)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* advance position half a time-step */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
         v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_dt));
         v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_dt));
/* reflecting boundary conditions */
         if (ipbc==2) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+2*nppmx+npoff] = -vx;       */
/*          }                                      */
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
               _mm_store_ps(&ppart[j+2*nppmx+npoff],v_vx);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             ppart[j+3*nppmx+npoff] = -vy;       */
/*          }                                      */
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
               _mm_store_ps(&ppart[j+3*nppmx+npoff],v_vy);
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+2*nppmx+npoff] = -vx;       */
/*          }                                      */
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
               _mm_store_ps(&ppart[j+2*nppmx+npoff],v_vx);
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy; */
         _mm_store_ps(&ppart[j+npoff],v_dx);
         _mm_store_ps(&ppart[j+nppmx+npoff],v_dy);
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         nn = 4*(nn - noff + mxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ppart[j+2*nppmx+npoff];
         vy = ppart[j+3*nppmx+npoff];
         vz = ppart[j+4*nppmx+npoff];
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + 4;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += 4*mxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + 4;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+3*nppmx+npoff] = -vy;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -vx;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
/*          cu[4*(i+noff+nxv*(j+moff))] += scu[4*(i+mxv*j)]; */
/*          cu[1+4*(i+noff+nxv*(j+moff))] += scu[1+4*(i+mxv*j)]; */
/*          cu[2+4*(i+noff+nxv*(j+moff))] += scu[2+4*(i+mxv*j)]; */
            v_x = _mm_loadu_ps(&cu[4*(i+noff+nxv*(j+moff))]);
            v_y = _mm_loadu_ps(&scu[4*(i+mxv*j)]);
            v_x = _mm_add_ps(v_x,v_y);
            _mm_storeu_ps(&cu[4*(i+noff+nxv*(j+moff))],v_x);
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[4*(i+noff+nxv*moff)] += scu[4*i];
#pragma omp atomic
         cu[1+4*(i+noff+nxv*moff)] += scu[1+4*i];
#pragma omp atomic
         cu[2+4*(i+noff+nxv*moff)] += scu[2+4*i];
         if (mm > my) {
#pragma omp atomic
            cu[4*(i+noff+nxv*(mm+moff-1))] += scu[4*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*(mm+moff-1))] += scu[1+4*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*(mm+moff-1))] += scu[2+4*(i+mxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[4*(noff+nxv*(j+moff))] += scu[4*mxv*j];
#pragma omp atomic
         cu[1+4*(noff+nxv*(j+moff))] += scu[1+4*mxv*j];
#pragma omp atomic
         cu[2+4*(noff+nxv*(j+moff))] += scu[2+4*mxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[4*(nn+noff-1+nxv*(j+moff))] += scu[4*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[1+4*(nn+noff-1+nxv*(j+moff))] += scu[1+4*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[2+4*(nn+noff-1+nxv*(j+moff))] += scu[2+4*((nn-1)+mxv*j)];
         }
      }
   }
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void csse2gjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                      int ihole[], float qm, float dt, int nppmx,
                      int idimp, int nx, int ny, int mx, int my,
                      int nxv, int nyv, int mx1, int mxy1, int ntmax,
                      int *irc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   vector/OpenMP version using guard cells
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x velocity of particle n in tile m
   ppart[m][3][n] = y velocity of particle n in tile m
   ppart[m][4][n] = z velocity of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of current array, must be >= nx+1
   nyv = second dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
   requires SSE2, ppart and cu need to be 16 byte aligned
   nppmx needs to be a multiple of 4, cu needs to have 4 components
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nps, ih, nh, nn, mm, kk;
   float dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz;
   float anx, any, edgelx, edgely, edgerx, edgery;
   __m128i v_noff, v_moff, v_mxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qm, v_dt, v_one;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_st, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_vx, v_vy;
   __m128 v_anx, v_any, v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 v_zero, v_two, v_three, v_six;
   __m128 a, b, c, d, va, vb, vc, vd;
   __attribute__((aligned(16))) unsigned int ll[4], lm[8];
   __attribute__((aligned(16))) unsigned long jj[1];
   __attribute__((aligned(16))) float scu[4*MXV*MYV];
/* __attribute__((aligned(16))) float scu[4*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   anx = (float) nx;
   any = (float) ny;
   v_mxv = _mm_set1_epi32(mxv);
   v_qm = _mm_set1_ps(qm);
   v_anx = _mm_set1_ps(anx);
   v_any = _mm_set1_ps(any);
   v_zero = _mm_setzero_ps();
   v_one = _mm_set1_ps(1.0f);
   v_dt = _mm_set1_ps(dt);
   v_two = _mm_set1_ps(2.0f);
   v_three = _mm_set1_ps(3.0f);
   v_six = _mm_set1_ps(6.0f);
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nps,nn,mm,kk,ih,nh,x,y,dxp,dyp,amx, \
amy,dx,dy,vx,vy,vz,edgelx,edgely,edgerx,edgery,v_noff,v_moff,v_nn,v_mm, \
v_it,v_x,v_y,v_vx,v_vy,v_dxp,v_dyp,v_amx,v_amy,v_dx,v_dy,v_at,v_st, \
v_edgelx,v_edgely,v_edgerx,v_edgery,a,b,c,d,va,vb,vc,vd,jj,ll,lm,scu)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm_set1_epi32(noff);
      v_moff = _mm_set1_epi32(moff);
      npp = kpic[k];
      nps = 4*(npp/4);
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
/* zero out local accumulator */
/*    for (j = 0; j < 4*mxv*(my+1); j++) { */
/*       scu[j] = 0.0f;                    */
/*    }                                    */
      memset((void*)scu,0,4*mxv*(my+1)*sizeof(float));
/* clear counters */
/*    for (j = 0; j < 8; j++) { */
/*       ncl[j+8*k] = 0;        */
/*    }                         */
      memset((void*)&ncl[8*k],0,8*sizeof(int));
/* vector loop over particles in blocks of 4 */
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
/*       dyp = y - (float) mm;      */
         v_dyp = _mm_sub_ps(v_y,_mm_cvtepi32_ps(v_mm));
/*       nm = 4*(nn - noff + mxv*(mm - moff)); */
         v_nn = _mm_sub_epi32(v_nn,v_noff);
         v_mm = _mm_sub_epi32(v_mm,v_moff);
         v_it = _mm_mul_epu32(v_mxv,_mm_srli_si128(v_mm,4));
         v_mm = _mm_mul_epu32(v_mm,v_mxv);
         v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
         v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
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
/* deposit current */
/*       vx = ppart[j+2*nppmx+npoff]; */
/*       vy = ppart[j+3*nppmx+npoff]; */
/*       vz = ppart[j+4*nppmx+npoff]; */
         v_vx = _mm_load_ps(&ppart[j+2*nppmx+npoff]);
         v_vy = _mm_load_ps(&ppart[j+3*nppmx+npoff]);
         va = v_vx;
         vb = v_vy;
         vc = _mm_load_ps(&ppart[j+4*nppmx+npoff]);
         vd = _mm_setzero_ps();
/* transpose so va,vb,vc,vd contain the 3 velocities plus zero */
/* for each of 4 particles                                     */
         _MM_TRANSPOSE4_PS(va,vb,vc,vd);
/*       dx = amx*amy;       */
/*       scu[nn] += vx*dx;   */
/*       scu[nn+1] += vy*dx; */
/*       scu[nn+2] += vz*dx; */
/*       dy = dxp*amy;       */
/*       mm = nn + 4;        */
/*       scu[mm] += vx*dy;   */
/*       scu[mm+1] += vy*dy; */
/*       scu[mm+2] += vz*dy; */
/*       dx = amx*dyp;       */
/*       nn += 4*nxv;        */
/*       scu[nn] += vx*dx;   */
/*       scu[nn+1] += vy*dx; */
/*       scu[nn+2] += vz*dx; */
/*       dy = dxp*dyp;       */
/*       mm = nn + 4;        */
/*       scu[mm] += vx*dy;   */
/*       scu[mm+1] += vy*dy; */
/*       scu[mm+2] += vz*dy; */
/* deposit for first particle */
         mm = ll[0];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(a,a,0)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(b,b,0)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(c,c,0)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(d,d,0)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for second particle */
         mm = ll[1];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(a,a,85)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(b,b,85)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(c,c,85)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(d,d,85)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for third particle */
         mm = ll[2];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(a,a,170)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(b,b,170)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(c,c,170)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(d,d,170)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for fourth particle */
         mm = ll[3];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(a,a,255)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(b,b,255)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(c,c,255)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(d,d,255)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* advance position half a time-step */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
         v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_dt));
         v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_dt));
/* find particles going out of bounds */
         mm = 0;
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
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy; */
         _mm_store_ps(&ppart[j+npoff],v_dx);
         _mm_store_ps(&ppart[j+nppmx+npoff],v_dy);
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
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         nn = 4*(nn - noff + mxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ppart[j+2*nppmx+npoff];
         vy = ppart[j+3*nppmx+npoff];
         vz = ppart[j+4*nppmx+npoff];
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + 4;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += 4*mxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + 4;
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
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
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
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
/*          cu[4*(i+noff+nxv*(j+moff))] += scu[4*(i+mxv*j)]; */
/*          cu[1+4*(i+noff+nxv*(j+moff))] += scu[1+4*(i+mxv*j)]; */
/*          cu[2+4*(i+noff+nxv*(j+moff))] += scu[2+4*(i+mxv*j)]; */
            v_x = _mm_loadu_ps(&cu[4*(i+noff+nxv*(j+moff))]);
            v_y = _mm_loadu_ps(&scu[4*(i+mxv*j)]);
            v_x = _mm_add_ps(v_x,v_y);
            _mm_storeu_ps(&cu[4*(i+noff+nxv*(j+moff))],v_x);
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[4*(i+noff+nxv*moff)] += scu[4*i];
#pragma omp atomic
         cu[1+4*(i+noff+nxv*moff)] += scu[1+4*i];
#pragma omp atomic
         cu[2+4*(i+noff+nxv*moff)] += scu[2+4*i];
         if (mm > my) {
#pragma omp atomic
            cu[4*(i+noff+nxv*(mm+moff-1))] += scu[4*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*(mm+moff-1))] += scu[1+4*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*(mm+moff-1))] += scu[2+4*(i+mxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[4*(noff+nxv*(j+moff))] += scu[4*mxv*j];
#pragma omp atomic
         cu[1+4*(noff+nxv*(j+moff))] += scu[1+4*mxv*j];
#pragma omp atomic
         cu[2+4*(noff+nxv*(j+moff))] += scu[2+4*mxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[4*(nn+noff-1+nxv*(j+moff))] += scu[4*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[1+4*(nn+noff-1+nxv*(j+moff))] += scu[1+4*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[2+4*(nn+noff-1+nxv*(j+moff))] += scu[2+4*((nn-1)+mxv*j)];
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
void csse2grjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                      float dt, float ci, int nppmx, int idimp, int nx,
                      int ny, int mx, int my, int nxv, int nyv, int mx1,
                      int mxy1, int ipbc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   vecgtor/OpenMP version using guard cells
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of current array, must be >= nx+1
   nyv = second dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   requires SSE2, ppart and cu need to be 16 byte aligned
   nppmx needs to be a multiple of 4, cu needs to have 4 components
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nps, nn, mm;
   float ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami;
   __m128i v_noff, v_moff, v_mxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qm, v_dt, v_one, v_ci2;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_vx, v_vy, v_ux, v_uy;
   __m128 v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 a, b, c, d, va, vb, vc, vd;
   __attribute__((aligned(16))) unsigned int ll[4];
   __attribute__((aligned(16))) unsigned long kk[1];
   __attribute__((aligned(16))) float scu[4*MXV*MYV];
/* __attribute__((aligned(16))) float scu[4*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   ci2 = ci*ci;
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
   v_qm = _mm_set1_ps(qm);
   v_one = _mm_set1_ps(1.0f);
   v_dt = _mm_set1_ps(dt);
   v_ci2 = _mm_set1_ps(ci2);
   v_edgelx = _mm_set1_ps(edgelx);
   v_edgely = _mm_set1_ps(edgely);
   v_edgerx = _mm_set1_ps(edgerx);
   v_edgery = _mm_set1_ps(edgery);
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nps,nn,mm,x,y,dxp,dyp,amx,amy,dx,dy, \
vx,vy,vz,ux,uy,uz,p2,gami,v_noff,v_moff,v_nn,v_mm,v_it,v_x,v_y,v_vx, \
v_vy,v_ux,v_uy,v_dxp,v_dyp,v_amx,v_amy,v_dx,v_dy,v_at,a,b,c,d,va,vb, \
vc,vd,ll,kk,scu)
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
/*    for (j = 0; j < 4*mxv*(my+1); j++) { */
/*       scu[j] = 0.0f;                    */
/*    }                                    */
      memset((void*)scu,0,4*mxv*(my+1)*sizeof(float));
/* vector loop over particles in blocks of 4 */
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
/*       dyp = y - (float) mm;      */
         v_dyp = _mm_sub_ps(v_y,_mm_cvtepi32_ps(v_mm));
/* find inverse gamma */
/*       ux = ppart[j+2*nppmx+npoff]; */
/*       uy = ppart[j+3*nppmx+npoff]; */
/*       uz = ppart[j+4*nppmx+npoff]; */
         v_ux = _mm_load_ps(&ppart[j+2*nppmx+npoff]);
         v_uy = _mm_load_ps(&ppart[j+3*nppmx+npoff]);
         vc = _mm_load_ps(&ppart[j+4*nppmx+npoff]);
/*       p2 = ux*ux + uy*uy + uz*uz; */
         v_at = _mm_mul_ps(v_ux,v_ux);
         v_at = _mm_add_ps(v_at,_mm_mul_ps(v_uy,v_uy));
         v_at = _mm_add_ps(v_at,_mm_mul_ps(vc,vc));
/*       gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_at = _mm_rsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* full accuracy calculation */
         v_at = _mm_sqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2)));
         v_at = _mm_div_ps(v_one,v_at);
/* full accuracy calculation with SVML */
/*       v_at = _mm_invsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* calculate weights */
/*       nm = 4*(nn - noff + mxv*(mm - moff)); */
         v_nn = _mm_sub_epi32(v_nn,v_noff);
         v_mm = _mm_sub_epi32(v_mm,v_moff);
         v_it = _mm_mul_epu32(v_mxv,_mm_srli_si128(v_mm,4));
         v_mm = _mm_mul_epu32(v_mm,v_mxv);
         v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
         v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
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
/* deposit current */
/*       vx = ux*gami; */
/*       vy = uy*gami; */
/*       vz = uz*gami; */
         v_vx = _mm_mul_ps(v_ux,v_at);
         v_vy = _mm_mul_ps(v_uy,v_at);
         va = v_vx;
         vb = v_vy;
         vc = _mm_mul_ps(vc,v_at);
         vd = _mm_setzero_ps();
/* transpose so va,vb,vc,vd contain the 3 velocities plus zero */
/* for each of 4 particles                                     */
         _MM_TRANSPOSE4_PS(va,vb,vc,vd);
/*       dx = amx*amy;       */
/*       scu[nn] += vx*dx;   */
/*       scu[nn+1] += vy*dx; */
/*       scu[nn+2] += vz*dx; */
/*       dy = dxp*amy;       */
/*       mm = nn + 4;        */
/*       scu[mm] += vx*dy;   */
/*       scu[mm+1] += vy*dy; */
/*       scu[mm+2] += vz*dy; */
/*       dx = amx*dyp;       */
/*       nn += 4*nxv;        */
/*       scu[nn] += vx*dx;   */
/*       scu[nn+1] += vy*dx; */
/*       scu[nn+2] += vz*dx; */
/*       dy = dxp*dyp;       */
/*       mm = nn + 4;        */
/*       scu[mm] += vx*dy;   */
/*       scu[mm+1] += vy*dy; */
/*       scu[mm+2] += vz*dy; */
/* deposit for first particle */
         mm = ll[0];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(a,a,0)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(b,b,0)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(c,c,0)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(d,d,0)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for second particle */
         mm = ll[1];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(a,a,85)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(b,b,85)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(c,c,85)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(d,d,85)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for third particle */
         mm = ll[2];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(a,a,170)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(b,b,170)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(c,c,170)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(d,d,170)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for fourth particle */
         mm = ll[3];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(a,a,255)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(b,b,255)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(c,c,255)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(d,d,255)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* advance position half a time-step */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
         v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_dt));
         v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_dt));
/* reflecting boundary conditions */
         if (ipbc==2) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+2*nppmx+npoff] = -ux;       */
/*          }                                      */
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
               _mm_store_ps(&ppart[j+2*nppmx+npoff],v_ux);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             ppart[j+3*nppmx+npoff] = -uy;       */
/*          }                                      */
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
               _mm_store_ps(&ppart[j+3*nppmx+npoff],v_uy);
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+2*nppmx+npoff] = -ux;       */
/*          }                                      */
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
               _mm_store_ps(&ppart[j+2*nppmx+npoff],v_ux);
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy; */
         _mm_store_ps(&ppart[j+npoff],v_dx);
         _mm_store_ps(&ppart[j+nppmx+npoff],v_dy);
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
/* find inverse gamma */
         ux = ppart[j+2*nppmx+npoff];
         uy = ppart[j+3*nppmx+npoff];
         uz = ppart[j+4*nppmx+npoff];
         p2 = ux*ux + uy*uy + uz*uz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* calculate weights */
         nn = 4*(nn - noff + mxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ux*gami;
         vy = uy*gami;
         vz = uz*gami;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + 4;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += 4*mxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + 4;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -ux;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+3*nppmx+npoff] = -uy;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+2*nppmx+npoff] = -ux;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
/*          cu[4*(i+noff+nxv*(j+moff))] += scu[4*(i+mxv*j)]; */
/*          cu[1+4*(i+noff+nxv*(j+moff))] += scu[1+4*(i+mxv*j)]; */
/*          cu[2+4*(i+noff+nxv*(j+moff))] += scu[2+4*(i+mxv*j)]; */
            v_x = _mm_loadu_ps(&cu[4*(i+noff+nxv*(j+moff))]);
            v_y = _mm_loadu_ps(&scu[4*(i+mxv*j)]);
            v_x = _mm_add_ps(v_x,v_y);
            _mm_storeu_ps(&cu[4*(i+noff+nxv*(j+moff))],v_x);
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[4*(i+noff+nxv*moff)] += scu[4*i];
#pragma omp atomic
         cu[1+4*(i+noff+nxv*moff)] += scu[1+4*i];
#pragma omp atomic
         cu[2+4*(i+noff+nxv*moff)] += scu[2+4*i];
         if (mm > my) {
#pragma omp atomic
            cu[4*(i+noff+nxv*(mm+moff-1))] += scu[4*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*(mm+moff-1))] += scu[1+4*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*(mm+moff-1))] += scu[2+4*(i+mxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[4*(noff+nxv*(j+moff))] += scu[4*mxv*j];
#pragma omp atomic
         cu[1+4*(noff+nxv*(j+moff))] += scu[1+4*mxv*j];
#pragma omp atomic
         cu[2+4*(noff+nxv*(j+moff))] += scu[2+4*mxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[4*(nn+noff-1+nxv*(j+moff))] += scu[4*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[1+4*(nn+noff-1+nxv*(j+moff))] += scu[1+4*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[2+4*(nn+noff-1+nxv*(j+moff))] += scu[2+4*((nn-1)+mxv*j)];
         }
      }
   }
   return;
#undef MXV
#undef MYV
}

/*--------------------------------------------------------------------*/
void csse2grjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                       int ihole[], float qm, float dt, float ci,
                       int nppmx, int idimp, int nx, int ny, int mx,
                       int my, int nxv, int nyv, int mx1, int mxy1,
                       int ntmax, int *irc) {
/* for 2-1/2d code, this subroutine calculates particle current density
   using first-order linear interpolation for relativistic particles
   in addition, particle positions are advanced a half time-step
   with periodic boundary conditions.
   also determines list of particles which are leaving this tile
   vector/OpenMP version using guard cells
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = x momentum of particle n in tile m
   ppart[m][3][n] = y momentum of particle n in tile m
   ppart[m][4][n] = z momentum of particle n in tile m
   cu[k][j][i] = ith component of current density at grid point j,k
   kpic[k] = number of particles in tile k
   ncl[k][i] = number of particles going to destination i, tile k
   ihole[k][:][0] = location of hole in array left by departing particle
   ihole[k][:][1] = destination of particle leaving hole
   ihole[k][0][0] = ih, number of holes left (error, if negative)
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
   nppmx = maximum number of particles in tile
   idimp = size of phase space = 5
   nx/ny = system length in x/y direction
   mx/my = number of grids in sorting cell in x/y
   nxv = first dimension of current array, must be >= nx+1
   nyv = second dimension of current array, must be >= ny+1
   mx1 = (system length in x direction - 1)/mx + 1
   mxy1 = mx1*my1, where my1 = (system length in y direction - 1)/my + 1
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   optimized version
   requires SSE2, ppart and cu need to be 16 byte aligned
   nppmx needs to be a multiple of 4, cu needs to have 4 components
local data                                                            */
#define MXV             33
#define MYV             33
   int noff, moff, npoff, npp, mxv;
   int i, j, k, nps, ih, nh, nn, mm, kk;
   float ci2, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy, vz, ux, uy, uz, p2, gami;
   float anx, any, edgelx, edgely, edgerx, edgery;
   __m128i v_noff, v_moff, v_mxv;
   __m128i v_nn, v_mm, v_it;
   __m128 v_qm, v_dt, v_one, v_ci2;
   __m128 v_dxp, v_dyp, v_amx, v_amy, v_st, v_at;
   __m128 v_x, v_y, v_dx, v_dy, v_vx, v_vy, v_ux, v_uy;
   __m128 v_anx, v_any, v_edgelx, v_edgely, v_edgerx, v_edgery;
   __m128 v_zero, v_two, v_three, v_six;
   __m128 a, b, c, d, va, vb, vc, vd;
   __attribute__((aligned(16))) unsigned int ll[4], lm[8];
   __attribute__((aligned(16))) unsigned long jj[1];
   __attribute__((aligned(16))) float scu[4*MXV*MYV];
/* __attribute__((aligned(16))) float scu[4*(mx+1)*(my+1)]; */
   mxv = mx + 1;
   ci2 = ci*ci;
   anx = (float) nx;
   any = (float) ny;
   v_mxv = _mm_set1_epi32(mxv);
   v_qm = _mm_set1_ps(qm);
   v_anx = _mm_set1_ps(anx);
   v_any = _mm_set1_ps(any);
   v_zero = _mm_setzero_ps();
   v_one = _mm_set1_ps(1.0f);
   v_dt = _mm_set1_ps(dt);
   v_ci2 = _mm_set1_ps(ci2);
   v_two = _mm_set1_ps(2.0f);
   v_three = _mm_set1_ps(3.0f);
   v_six = _mm_set1_ps(6.0f);
/* error if local array is too small */
/* if ((mx >= MXV) || (my >= MYV))   */
/*    return;                        */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,noff,moff,npp,npoff,nps,nn,mm,kk,ih,nh,x,y,dxp,dyp,amx, \
amy,dx,dy,vx,vy,vz,ux,uy,uz,edgelx,edgely,edgerx,edgery,p2,gami,v_noff, \
v_moff,v_nn,v_mm,v_it,v_x,v_y,v_vx,v_vy,v_ux,v_uy,v_dxp,v_dyp,v_amx, \
v_amy,v_dx,v_dy,v_at,v_st,v_edgelx,v_edgely,v_edgerx,v_edgery,a,b,c,d, \
va,vb,vc,vd,jj,ll,lm,scu)
   for (k = 0; k < mxy1; k++) {
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm_set1_epi32(noff);
      v_moff = _mm_set1_epi32(moff);
      npp = kpic[k];
      nps = 4*(npp/4);
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
/* zero out local accumulator */
/*    for (j = 0; j < 4*mxv*(my+1); j++) { */
/*       scu[j] = 0.0f;                    */
/*    }                                    */
      memset((void*)scu,0,4*mxv*(my+1)*sizeof(float));
/* clear counters */
/*    for (j = 0; j < 8; j++) { */
/*       ncl[j+8*k] = 0;        */
/*    }                         */
      memset((void*)&ncl[8*k],0,8*sizeof(int));
/* vector loop over particles in blocks of 4 */
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
/*       dyp = y - (float) mm;      */
         v_dyp = _mm_sub_ps(v_y,_mm_cvtepi32_ps(v_mm));
/* find inverse gamma */
/*       ux = ppart[j+2*nppmx+npoff]; */
/*       uy = ppart[j+3*nppmx+npoff]; */
/*       uz = ppart[j+4*nppmx+npoff]; */
         v_ux = _mm_load_ps(&ppart[j+2*nppmx+npoff]);
         v_uy = _mm_load_ps(&ppart[j+3*nppmx+npoff]);
         vc = _mm_load_ps(&ppart[j+4*nppmx+npoff]);
/*       p2 = ux*ux + uy*uy + uz*uz; */
         v_at = _mm_mul_ps(v_ux,v_ux);
         v_at = _mm_add_ps(v_at,_mm_mul_ps(v_uy,v_uy));
         v_at = _mm_add_ps(v_at,_mm_mul_ps(vc,vc));
/*       gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_at = _mm_rsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* full accuracy calculation */
         v_at = _mm_sqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2)));
         v_at = _mm_div_ps(v_one,v_at);
/* full accuracy calculation with SVML */
/*       v_at = _mm_invsqrt_ps(_mm_add_ps(v_one,_mm_mul_ps(v_at,v_ci2))); */
/* calculate weights */
/*       nm = 4*(nn - noff + mxv*(mm - moff)); */
         v_nn = _mm_sub_epi32(v_nn,v_noff);
         v_mm = _mm_sub_epi32(v_mm,v_moff);
         v_it = _mm_mul_epu32(v_mxv,_mm_srli_si128(v_mm,4));
         v_mm = _mm_mul_epu32(v_mm,v_mxv);
         v_mm = _mm_add_epi32(v_mm,_mm_slli_si128(v_it,4));
         v_nn = _mm_slli_epi32(_mm_add_epi32(v_nn,v_mm),2);
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
/* deposit current */
/*       vx = ux*gami; */
/*       vy = uy*gami; */
/*       vz = uz*gami; */
         v_vx = _mm_mul_ps(v_ux,v_at);
         v_vy = _mm_mul_ps(v_uy,v_at);
         va = v_vx;
         vb = v_vy;
         vc = _mm_mul_ps(vc,v_at);
         vd = _mm_setzero_ps();
/* transpose so va,vb,vc,vd contain the 3 velocities plus zero */
/* for each of 4 particles                                     */
         _MM_TRANSPOSE4_PS(va,vb,vc,vd);
/*       dx = amx*amy;       */
/*       scu[nn] += vx*dx;   */
/*       scu[nn+1] += vy*dx; */
/*       scu[nn+2] += vz*dx; */
/*       dy = dxp*amy;       */
/*       mm = nn + 4;        */
/*       scu[mm] += vx*dy;   */
/*       scu[mm+1] += vy*dy; */
/*       scu[mm+2] += vz*dy; */
/*       dx = amx*dyp;       */
/*       nn += 4*nxv;        */
/*       scu[nn] += vx*dx;   */
/*       scu[nn+1] += vy*dx; */
/*       scu[nn+2] += vz*dx; */
/*       dy = dxp*dyp;       */
/*       mm = nn + 4;        */
/*       scu[mm] += vx*dy;   */
/*       scu[mm+1] += vy*dy; */
/*       scu[mm+2] += vz*dy; */
/* deposit for first particle */
         mm = ll[0];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(a,a,0)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(b,b,0)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(va,_mm_shuffle_ps(c,c,0)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(va,_mm_shuffle_ps(d,d,0)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for second particle */
         mm = ll[1];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(a,a,85)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(b,b,85)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vb,_mm_shuffle_ps(c,c,85)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vb,_mm_shuffle_ps(d,d,85)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for third particle */
         mm = ll[2];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(a,a,170)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(b,b,170)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vc,_mm_shuffle_ps(c,c,170)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vc,_mm_shuffle_ps(d,d,170)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* deposit for fourth particle */
         mm = ll[3];
         v_dx = _mm_load_ps(&scu[mm]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(a,a,255)));
         _mm_store_ps(&scu[mm],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(b,b,255)));
         _mm_store_ps(&scu[mm+4],v_dy);
         v_dx = _mm_load_ps(&scu[mm+4*mxv]);
         v_dx = _mm_add_ps(v_dx,_mm_mul_ps(vd,_mm_shuffle_ps(c,c,255)));
         _mm_store_ps(&scu[mm+4*mxv],v_dx);
         v_dy = _mm_load_ps(&scu[mm+4+4*mxv]);
         v_dy = _mm_add_ps(v_dy,_mm_mul_ps(vd,_mm_shuffle_ps(d,d,255)));
         _mm_store_ps(&scu[mm+4+4*mxv],v_dy);
/* advance position half a time-step */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
         v_dx = _mm_add_ps(v_x,_mm_mul_ps(v_vx,v_dt));
         v_dy = _mm_add_ps(v_y,_mm_mul_ps(v_vy,v_dt));
/* find particles going out of bounds */
         mm = 0;
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
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy; */
         _mm_store_ps(&ppart[j+npoff],v_dx);
         _mm_store_ps(&ppart[j+nppmx+npoff],v_dy);
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
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         nn = x;
         mm = y;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
/* find inverse gamma */
         ux = ppart[j+2*nppmx+npoff];
         uy = ppart[j+3*nppmx+npoff];
         uz = ppart[j+4*nppmx+npoff];
         p2 = ux*ux + uy*uy + uz*uz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* calculate weights */
         nn = 4*(nn - noff + mxv*(mm - moff));
         amx = qm - dxp;
         amy = 1.0f - dyp;
/* deposit current */
         dx = amx*amy;
         dy = dxp*amy;
         vx = ux*gami;
         vy = uy*gami;
         vz = uz*gami;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = amx*dyp;
         mm = nn + 4;
         scu[mm] += vx*dy;
         scu[mm+1] += vy*dy;
         scu[mm+2] += vz*dy;
         dy = dxp*dyp;
         nn += 4*mxv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         mm = nn + 4;
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
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
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
      nn = nxv - noff;
      mm = nyv - moff;
      nn = mx < nn ? mx : nn;
      mm = my < mm ? my : mm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
/*          cu[4*(i+noff+nxv*(j+moff))] += scu[4*(i+mxv*j)]; */
/*          cu[1+4*(i+noff+nxv*(j+moff))] += scu[1+4*(i+mxv*j)]; */
/*          cu[2+4*(i+noff+nxv*(j+moff))] += scu[2+4*(i+mxv*j)]; */
            v_x = _mm_loadu_ps(&cu[4*(i+noff+nxv*(j+moff))]);
            v_y = _mm_loadu_ps(&scu[4*(i+mxv*j)]);
            v_x = _mm_add_ps(v_x,v_y);
            _mm_storeu_ps(&cu[4*(i+noff+nxv*(j+moff))],v_x);
         }
      }
/* deposit current to edge points in global array */
      mm = nyv - moff;
      mm = my+1 < mm ? my+1 : mm;
      for (i = 1; i < nn; i++) {
#pragma omp atomic
         cu[4*(i+noff+nxv*moff)] += scu[4*i];
#pragma omp atomic
         cu[1+4*(i+noff+nxv*moff)] += scu[1+4*i];
#pragma omp atomic
         cu[2+4*(i+noff+nxv*moff)] += scu[2+4*i];
         if (mm > my) {
#pragma omp atomic
            cu[4*(i+noff+nxv*(mm+moff-1))] += scu[4*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*(mm+moff-1))] += scu[1+4*(i+mxv*(mm-1))];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*(mm+moff-1))] += scu[2+4*(i+mxv*(mm-1))];
         }
      }
      nn = nxv - noff;
      nn = mx+1 < nn ? mx+1 : nn;
      for (j = 0; j < mm; j++) {
#pragma omp atomic
         cu[4*(noff+nxv*(j+moff))] += scu[4*mxv*j];
#pragma omp atomic
         cu[1+4*(noff+nxv*(j+moff))] += scu[1+4*mxv*j];
#pragma omp atomic
         cu[2+4*(noff+nxv*(j+moff))] += scu[2+4*mxv*j];
         if (nn > mx) {
#pragma omp atomic
            cu[4*(nn+noff-1+nxv*(j+moff))] += scu[4*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[1+4*(nn+noff-1+nxv*(j+moff))] += scu[1+4*((nn-1)+mxv*j)];
#pragma omp atomic
            cu[2+4*(nn+noff-1+nxv*(j+moff))] += scu[2+4*((nn-1)+mxv*j)];
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
void csse2bguard2l(float bxy[], int nx, int ny, int nxe, int nye) {
/* replicate extended periodic vector field bxy
   linear interpolation
   nx/ny = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nxe = second dimension of field arrays, must be >= ny+1
   requires SSE2, bxy needs to be 16 byte aligned
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
   nxe = second dimension of field arrays, must be >= ny+1
   requires SSE2, bxy needs to be 16 byte aligned
local data                                                 */
   int j, k, kk;
   __m128 v_cu;
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
      v_cu = _mm_add_ps(_mm_load_ps(&cu[4*j]),_mm_load_ps(&cu[4*j+kk]));
      _mm_store_ps(&cu[4*j],v_cu);
/*    cu[4*j+kk] = 0.0;   */
/*    cu[1+4*j+kk] = 0.0; */
/*    cu[2+4*j+kk] = 0.0; */
      _mm_store_ps(&cu[4*j+kk],_mm_setzero_ps());
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
void csse2mpois23(float complex q[], float complex fxy[], int isign,
                  float complex ffc[], float ax, float ay, float affp,
                  float *we, int nx, int ny, int nxvh, int nyv,
                  int nxhd, int nyhd) {
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
   requires SSE2, q, fxy, ffc need to be 16 byte aligned
   nxhd, nxvh need to be a multiple of 2
   fxy needs to have 4 components
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
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,kk,kj,at1,at2,at3,zt1,zt2,wp,v_it,v_dky,v_at1,v_at2, \
v_at3,v_at4,v_zt1,v_zt2,v_zt3,v_zt4,v_wp,v_d,dd) \
reduction(+:sum1)
   for (k = 1; k < nyh; k++) {
/*    dky = dny*(float) k; */
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
/*       fxy[4*(j+kj)] = at2*zt1;   */
/*       fxy[1+4*(j+kj)] = at3*zt1; */
/*       fxy[2+4*(j+kj)] = zero;    */
         v_at4 = _mm_mul_ps(v_at2,v_zt1);
         v_zt4 = _mm_mul_ps(v_at3,v_zt1);
/* reorder write */
         v_zt3 = _mm_shuffle_ps(v_at4,v_zt4,68);
         v_zt4 = _mm_shuffle_ps(v_at4,v_zt4,238);
         _mm_store_ps((float *)&fxy[4*(j+kj)],v_zt3);
         _mm_store_ps((float *)&fxy[2+4*(j+kj)],v_zero);
         _mm_store_ps((float *)&fxy[4*(j+1+kj)],v_zt4);
         _mm_store_ps((float *)&fxy[2+4*(j+1+kj)],v_zero);
/*       fxy[4*(j+k1)] = at2*zt2;    */
/*       fxy[1+4*(j+k1)] = -at3*zt2; */
/*       fxy[2+4*(j+k1)] = zero;     */
         v_at4 = _mm_mul_ps(v_at2,v_zt2);
         v_zt4 = _mm_sub_ps(v_zero,_mm_mul_ps(v_at3,v_zt2));
/* reorder write */
         v_zt3 = _mm_shuffle_ps(v_at4,v_zt4,68);
         v_zt4 = _mm_shuffle_ps(v_at4,v_zt4,238);
         _mm_store_ps((float *)&fxy[4*(j+k1)],v_zt3);
         _mm_store_ps((float *)&fxy[2+4*(j+k1)],v_zero);
         _mm_store_ps((float *)&fxy[4*(j+1+k1)],v_zt4);
         _mm_store_ps((float *)&fxy[2+4*(j+1+k1)],v_zero);
/*       at1 = at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1])); */
/*       wp += (double) at1;                                          */
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
#pragma ivdep
      for (j = it; j < nxh; j++) {
         at1 = crealf(ffc[j+kk])*cimagf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         zt1 = cimagf(q[j+kj]) - crealf(q[j+kj])*_Complex_I;
         zt2 = cimagf(q[j+k1]) - crealf(q[j+k1])*_Complex_I;
         fxy[4*(j+kj)] = at2*zt1;
         fxy[1+4*(j+kj)] = at3*zt1;
         fxy[2+4*(j+kj)] = zero;
         fxy[4*(j+k1)] = at2*zt2;
         fxy[1+4*(j+k1)] = -at3*zt2;
         fxy[2+4*(j+k1)] = zero;
         at1 = at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1]));
         wp += (double) at1;
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
#pragma ivdep
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
      at1 = at1*(q[kj]*conjf(q[kj]));
      wp += (double) at1;
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
/*    at1 = at1*(q[j]*conjf(q[j])); */
/*    wp += (double) at1;           */
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
#pragma ivdep
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
      at1 = at1*(q[j]*conjf(q[j]));
      wp += (double) at1;
   }
   fxy[0] = zero;
   fxy[1] = zero;
   fxy[2] = zero;
   fxy[k1] = zero;
   fxy[1+k1] = zero;
   fxy[2+k1] = zero;
   sum1 += wp;
/* *we = wp*(float) (nx*ny); */
   _mm_store_pd(&dd[0],v_wp);
   for (j = 1; j < 2; j++) {
      dd[0] += dd[j];
   }
   *we = (sum1 + dd[0])*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void csse2mcuperp2(float complex cu[], int nx, int ny, int nxvh,
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
#pragma omp parallel for private(j,k,k1,kj,v_dkx,v_dky,v_dky2,v_at1, \
v_zt1,v_zt2,v_at)
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
/* mode numbers kx = 0, nx/2 */
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
/* cu[k1] = zero;   */
/* cu[1+k1] = zero; */
   _mm_store_ps((float *)&cu[k1],v_zero);
   return;
}

/*--------------------------------------------------------------------*/
void csse2mibpois23(float complex cu[], float complex bxy[],
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
   double wp, sum1;
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
   sum1 = 0.0;
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,kk,kj,v_it,v_dky,v_at1,v_at2,v_at3,v_at4,v_zt1,v_zt2, \
v_zt3,v_zt4,v_wp,v_d,dd) \
reduction(+:sum1)
   for (k = 1; k < nyh; k++) {
/*    dky = dny*(float) k; */
      v_dky = _mm_mul_ps(v_dny,_mm_cvtepi32_ps(_mm_set1_epi32(k)));
      kk = nxhd*k;
      kj = 4*nxvh*k;
      k1 = 4*nxvh*ny - kj;
      v_wp = _mm_set1_pd(0.0);
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
/*       at1 = at1*(cu[4*j+kj]*conjf(cu[4*j+kj])    */
/*             + cu[1+4*j+kj]*conjf(cu[1+4*j+kj])   */
/*             + cu[2+4*j+kj]*conjf(cu[2+4*j+kj])); */
/*       wp += (double) at1;                        */
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
/*       at1 = at1*(cu[4*j+k1]*conjf(cu[4*j+k1])    */
/*             + cu[1+4*j+k1]*conjf(cu[1+4*j+k1])   */
/*             + cu[2+4*j+k1]*conjf(cu[2+4*j+k1])); */
/*       wp += (double) at1;                        */
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
/*    sum1 += wp; */
      _mm_store_pd(&dd[0],v_wp);
      for (j = 1; j < 2; j++) {
         dd[0] += dd[j];
      }
      sum1 += dd[0];
   }
   wp = 0.0;
   v_wp = _mm_set1_pd(0.0);
/* mode numbers kx = 0, nx/2 */
#pragma ivdep
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
      bxy[k1] = zero;
      bxy[1+k1] = zero;
      bxy[2+k1] = zero;
      at1 = at1*(cu[kj]*conjf(cu[kj]) + cu[1+kj]*conjf(cu[1+kj])
            + cu[2+kj]*conjf(cu[2+kj]));
      wp += (double) at1;
   }
   sum1 += wp;
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
/*    at1 = at1*(cu[4*j]*conjf(cu[4*j]) + cu[1+4*j]*conjf(cu[1+4*j]) */
/*          + cu[2+4*j]*conjf(cu[2+4*j]));                           */
/*    wp += (double) at1;                                            */
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
   *wm = (sum1 + dd[0])*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void csse2mmaxwel2(float complex exy[], float complex bxy[],
                   float complex cu[], float complex ffc[], float ci,
                   float dt, float *wf, float *wm, int nx, int ny,
                   int nxvh, int nyv, int nxhd, int nyhd) {
/* this subroutine solves 2-1/2d maxwell's equation in fourier space for
   transverse electric and magnetic fields with periodic boundary
   conditions
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
   float at1;
   float complex zero, zt1, zt3, zt4, zt6, zt7, zt9;
   double wp, ws, sum1, sum2;
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
   sum1 = 0.0;
   sum2 = 0.0;
/* calculate the electromagnetic fields */
/* mode numbers 0 < kx < nx/2 and 0 < ky < ny/2 */
#pragma omp parallel for \
private(j,k,k1,kk,kj,v_it,v_dkx,v_dky,v_afdt,v_at1,v_at2,v_at3,v_at4, \
v_zt1,v_zt2,v_zt4,v_zt6,v_zt7,v_zt9,v_ws,v_wp,v_d,dd) \
reduction(+:sum1,sum2)
   for (k = 1; k < nyh; k++) {
/*    dky = dny*(float) k; */
      v_dky = _mm_mul_ps(v_dny,_mm_cvtepi32_ps(_mm_set1_epi32(k)));
      kk = nxhd*k;
      kj = 4*nxvh*k;
      k1 = 4*nxvh*ny - kj;
      v_wp = _mm_set1_pd(0.0);
      v_ws = _mm_set1_pd(0.0);
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
/*       at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
/*       ws += (double) at1;                                             */
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
/*       at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
/*       wp += (double) at1;                                             */
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
/*       zt4 = bxy[4*j+k1] + dth*(dky*zt1);   */
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
/*       at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
/*       ws += (double) at1;                                             */
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
/*       at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
/*       wp += (double) at1;                                             */
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
/*    sum1 += ws; */
      _mm_store_pd(&dd[0],v_ws);
      for (j = 1; j < 2; j++) {
         dd[0] += dd[j];
      }
      sum1 += dd[0];
/*    sum2 += wp; */
      _mm_store_pd(&dd[0],v_wp);
      for (j = 1; j < 2; j++) {
         dd[0] += dd[j];
      }
      sum2 += dd[0];
   }
   ws = 0.0;
   wp = 0.0;
   v_wp = _mm_set1_pd(0.0);
   v_ws = _mm_set1_pd(0.0);
/* mode numbers kx = 0, nx/2 */
#pragma ivdep
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
      at1 = anorm*(zt7*conjf(zt7) + zt9*conjf(zt9));
      ws += (double) at1;
      zt4 -= dth*(dky*zt1);
      zt6 += dth*(dky*zt3);
      bxy[kj] = zt4;
      bxy[1+kj] = zero;
      bxy[2+kj] = zt6;
      at1 = anorm*(zt4*conjf(zt4) + zt6*conjf(zt6));
      wp += (double) at1;
      bxy[k1] = zero;
      bxy[1+k1] = zero;
      bxy[2+k1] = zero;
      exy[k1] = zero;
      exy[1+k1] = zero;
      exy[2+k1] = zero;
   }
   sum1 += ws;
   sum2 += wp;
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
/*    at1 = anorm*(zt8*conjf(zt8) + zt9*conjf(zt9)); */
/*    ws += (double) at1;                            */
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
/*    at1 = anorm*(zt5*conjf(zt5) + zt6*conjf(zt6)); */
/*    wp += (double) at1;                            */
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
/* *wf = sum1*(float) (nx*ny);    */
   _mm_store_pd(&dd[0],v_ws);
   for (j = 1; j < 2; j++) {
      dd[0] += dd[j];
   }
   *wf = (sum1 + dd[0])*(float) (nx*ny);
/* *wm = sum2*c2*(float) (nx*ny); */
      _mm_store_pd(&dd[0],v_wp);
      for (j = 1; j < 2; j++) {
         dd[0] += dd[j];
      }
   *wm = (sum2 + dd[0])*c2*(float) (nx*ny);
   return;
}

/*--------------------------------------------------------------------*/
void csse2memfield2(float complex fxy[], float complex exy[],
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
#pragma omp parallel for private(j,k,k1,kk,kj,v_at1,v_zt1,v_zt2)
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = 4*nxvh*k;
         k1 = 4*nxvh*ny - kj;
         for (j = 0; j < nxh; j++) {
/*          at1 = cimagf(ffc[j+kk]); */
            v_at1 = _mm_loadl_pi(v_zero,(__m64 *)&ffc[j+kk]);
            v_at1 = _mm_movelh_ps(v_at1,v_at1);
            v_at1 = _mm_shuffle_ps(v_at1,v_at1,245);
/*          fxy[4*j+kj] += exy[4*j+kj]*at1;     */
/*          fxy[1+4*j+kj] += exy[1+4*j+kj]*at1; */
/*          fxy[2+4*j+kj] += exy[2+4*j+kj]*at1; */
            v_zt1 = _mm_load_ps((float *)&exy[4*j+kj]);
            v_zt2 = _mm_load_ps((float *)&fxy[4*j+kj]);
            v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
            _mm_store_ps((float *)&fxy[4*j+kj],v_zt2);
            v_zt1 = _mm_load_ps((float *)&exy[2+4*j+kj]);
            v_zt2 = _mm_load_ps((float *)&fxy[2+4*j+kj]);
            v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
            _mm_store_ps((float *)&fxy[2+4*j+kj],v_zt2);
/*          fxy[4*j+k1] += exy[4*j+k1]*at1;     */
/*          fxy[1+4*j+k1] += exy[1+4*j+k1]*at1; */
/*          fxy[2+4*j+k1] += exy[2+4*j+k1]*at1; */
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
/*       fxy[4*j] += exy[4*j]*at1;     */
/*       fxy[1+4*j] += exy[1+4*j]*at1; */
/*       fxy[2+4*j] += exy[2+4*j]*at1; */
         v_zt1 = _mm_load_ps((float *)&exy[4*j]);
         v_zt2 = _mm_load_ps((float *)&fxy[4*j]);
         v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
         _mm_store_ps((float *)&fxy[4*j],v_zt2);
         v_zt1 = _mm_load_ps((float *)&exy[2+4*j]);
         v_zt2 = _mm_load_ps((float *)&fxy[2+4*j]);
         v_zt2 = _mm_add_ps(v_zt2,_mm_mul_ps(v_zt1,v_at1));
         _mm_store_ps((float *)&fxy[2+4*j],v_zt2);
/*       fxy[4*j+k1] += exy[4*j+k1]*at1;     */
/*       fxy[1+4*j+k1] += exy[1+4*j+k1]*at1; */
/*       fxy[2+4*j+k1] += exy[2+4*j+k1]*at1; */
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
#pragma omp parallel for private(j,k,k1,kk,kj,v_at1,v_zt1,v_zt2)
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = 4*nxvh*k;
         k1 = 4*nxvh*ny - kj;
         for (j = 0; j < nxh; j++) {
/*          at1 = cimagf(ffc[j+kk]); */
            v_at1 = _mm_loadl_pi(v_zero,(__m64 *)&ffc[j+kk]);
            v_at1 = _mm_movelh_ps(v_at1,v_at1);
            v_at1 = _mm_shuffle_ps(v_at1,v_at1,245);
/*          fxy[4*j+kj] = exy[4*j+kj]*at1;     */
/*          fxy[1+4*j+kj] = exy[1+4*j+kj]*at1; */
/*          fxy[2+4*j+kj] = exy[2+4*j+kj]*at1; */
            v_zt1 = _mm_load_ps((float *)&exy[4*j+kj]);
            v_zt2 = _mm_mul_ps(v_zt1,v_at1);
            _mm_store_ps((float *)&fxy[4*j+kj],v_zt2);
            v_zt1 = _mm_load_ps((float *)&exy[2+4*j+kj]);
            v_zt2 = _mm_mul_ps(v_zt1,v_at1);
            _mm_store_ps((float *)&fxy[2+4*j+kj],v_zt2);
/*          fxy[4*j+k1] = exy[4*j+k1]*at1;     */
/*          fxy[1+4*j+k1] = exy[1+4*j+k1]*at1; */
/*          fxy[2+4*j+k1] = exy[2+4*j+k1]*at1; */
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
/*       fxy[4*j] = exy[4*j]*at1;     */
/*       fxy[1+4*j] = exy[1+4*j]*at1; */
/*       fxy[2+4*j] = exy[2+4*j]*at1; */
         v_zt1 = _mm_load_ps((float *)&exy[4*j]);
         v_zt2 = _mm_mul_ps(v_zt1,v_at1);
         _mm_store_ps((float *)&fxy[4*j],v_zt2);
         v_zt1 = _mm_load_ps((float *)&exy[2+4*j]);
         v_zt2 = _mm_mul_ps(v_zt1,v_at1);
         _mm_store_ps((float *)&fxy[2+4*j],v_zt2);
/*       fxy[4*j+k1] = exy[4*j+k1]*at1;     */
/*       fxy[1+4*j+k1] = exy[1+4*j+k1]*at1; */
/*       fxy[2+4*j+k1] = exy[2+4*j+k1]*at1; */
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
void csse2fft2rm3x(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nyi,
                   int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the x part of 3 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   y, using complex arithmetic, with OpenMP
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
   nxhd = second dimension of f >= nx/2
   nyd = third dimension of f >= ny
   nxhyd = maximum of (nx/2,ny)
   nxyhd = maximum of (nx,ny)/2
   fourier coefficients are stored as follows:
   f[k][j][0:2] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][1][0:2] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0][0:2]) = real part of mode nx/2,0 and
   imag(f[0][ny/2][0:2]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
   requires SSE2, f needs to be 16 byte aligned
   f needs to have 4 components
local data                                                            */
   int indx1, indx1y, nx, nxh, nxhh, ny, nxy, nxhy, nyt;
   int nrx, i, j, k, l, jj, j1, k1, k2, ns, ns2, km, kmr, joff;
   int nrxb;
   float ani;
/* float complex t1, t2, t3, t4; */
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
      goto L100;
/* inverse fourier transform */
   nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,joff,ani,v_t1,v_t2,v_t3, \
v_t4,v_t5,v_ani)
   for (i = nyi-1; i < nyt; i++) {
      joff = 4*nxhd*i;
/* swap complex components */
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
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
/*          t1 = f[4*j1+joff];   */
/*          t2 = f[1+4*j1+joff]; */
/*          t3 = f[2+4*j1+joff]; */
            v_t1 = _mm_load_ps((float *)&f[4*j1+joff]);
            v_t3 = _mm_load_ps((float *)&f[2+4*j1+joff]);
/*          f[4*j1+joff] = f[4*j+joff];     */
/*          f[1+4*j1+joff] = f[1+4*j+joff]; */
/*          f[2+4*j1+joff] = f[2+4*j+joff]; */
            v_t2 = _mm_load_ps((float *)&f[4*j+joff]);
            _mm_store_ps((float *)&f[4*j1+joff],v_t2);
            v_t2 = _mm_load_ps((float *)&f[2+4*j+joff]);
            _mm_store_ps((float *)&f[2+4*j1+joff],v_t2);
/*          f[4*j+joff] = t1;   */
/*          f[1+4*j+joff] = t2; */
/*          f[2+4*j+joff] = t3; */
            _mm_store_ps((float *)&f[4*j+joff],v_t1);
            _mm_store_ps((float *)&f[2+4*j+joff],v_t3);
         }
      }
/* then transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = 4*ns2*k;
            k2 = k1 + 4*ns;
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
      ani = 2.0*ani;
      for (jj = 0; jj < 3; jj++) {
         f[jj+4*nxhh+joff] = ani*conjf(f[jj+4*nxhh+joff]);
         f[jj+joff] = ani*((crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
L100: nrxb = nxhy/nxh;
   nrx = nxy/nxh;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,jj,j1,joff,v_t1,v_t2,v_t3,v_t4, \
v_t5,v_ani)
   for (i = nyi-1; i < nyt; i++) {
      joff = 4*nxhd*i;
/* scramble coefficients */
      kmr = nxy/nx;
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
      for (jj = 0; jj < 3; jj++) {
         f[jj+4*nxhh+joff] = 2.0*conjf(f[jj+4*nxhh+joff]);
         f[jj+joff] = (crealf(f[jj+joff]) + cimagf(f[jj+joff]))
                      + (crealf(f[jj+joff]) - cimagf(f[jj+joff]))*_Complex_I;
      }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
/*          t1 = f[4*j1+joff];   */
/*          t2 = f[1+4*j1+joff]; */
/*          t3 = f[2+4*j1+joff]; */
            v_t1 = _mm_load_ps((float *)&f[4*j1+joff]);
            v_t3 = _mm_load_ps((float *)&f[2+4*j1+joff]);
/*          f[4*j1+joff] = f[4*j+joff];     */
/*          f[1+4*j1+joff] = f[1+4*j+joff]; */
/*          f[2+4*j1+joff] = f[2+4*j+joff]; */
            v_t2 = _mm_load_ps((float *)&f[4*j+joff]);
            _mm_store_ps((float *)&f[4*j1+joff],v_t2);
            v_t2 = _mm_load_ps((float *)&f[2+4*j+joff]);
            _mm_store_ps((float *)&f[2+4*j1+joff],v_t2);
/*          f[4*j+joff] = t1;   */
/*          f[1+4*j+joff] = t2; */
/*          f[2+4*j+joff] = t3; */
            _mm_store_ps((float *)&f[4*j+joff],v_t1);
            _mm_store_ps((float *)&f[2+4*j+joff],v_t3);
         }
      }
/* then transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         for (k = 0; k < km; k++) {
            k1 = 4*ns2*k;
            k2 = k1 + 4*ns;
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
         ns = ns2;
      }
/* swap complex components */
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
void csse2fft2rm3y(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxi,
                   int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
/* this subroutine performs the y part of 3 two dimensional real to
   complex fast fourier transforms, and their inverses, for a subset of
   x, using complex arithmetic, with OpenMP
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
   nxhd = second dimension of f >= nx/2
   nyd = third dimension of f >= ny
  nxhyd = maximum of (nx/2,ny)
  nxyhd = maximum of (nx,ny)/2
  fourier coefficients are stored as follows:
   f[k][j][0:2] = real, imaginary part of mode j,k, where
   0 <= j < nx/2 and 0 <= k < ny, except for
   f[k][1][0:2] = real, imaginary part of mode nx/2,k, where
   ny/2+1 <= k < ny, and
   imag(f[0][0][0:2]) = real part of mode nx/2,0 and
   imag(f[0][ny/2][0:2]) = real part of mode nx/2,ny/2
   written by viktor k. decyk, ucla
   requires SSE2, f needs to be 16 byte aligned
   f needs to have 4 components
local data                                                            */
   int indx1, indx1y, nx, ny, nyh, nxy, nxhy, nxt;
   int nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr, koff;
   int nryb;
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
      goto L80;
/* inverse fourier transform */
   nryb = nxhy/ny;
   nry = nxy/ny;
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,koff,v_t1,v_t2,v_t3,v_t4,v_t5)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = 4*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = 4*nxhd*k1;
/*          t1 = f[4*i+k1];   */
/*          t2 = f[1+4*i+k1]; */
/*          t3 = f[2+4*i+k1]; */
            v_t1 = _mm_load_ps((float *)&f[4*i+k1]);
            v_t3 = _mm_load_ps((float *)&f[2+4*i+k1]);
/*          f[4*i+k1] = f[4*i+koff];     */
/*          f[1+4*i+k1] = f[1+4*i+koff]; */
/*          f[2+4*i+k1] = f[2+4*i+koff]; */
            v_t2 = _mm_load_ps((float *)&f[4*i+koff]);
            _mm_store_ps((float *)&f[4*i+k1],v_t2);
            v_t2 = _mm_load_ps((float *)&f[2+4*i+koff]);
            _mm_store_ps((float *)&f[2+4*i+k1],v_t2);
/*          f[4*i+koff] = t1;   */
/*          f[1+4*i+koff] = t2; */
/*          f[2+4*i+koff] = t3; */
            _mm_store_ps((float *)&f[4*i+koff],v_t1);
            _mm_store_ps((float *)&f[2+4*i+koff],v_t3);
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
               j1 = 4*nxhd*(j + k1);
               j2 = 4*nxhd*(j + k2);
/*             t1 = sct[kmr*j]; */
               v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
               v_t1 = _mm_movelh_ps(v_t1,v_t1);
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
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if (nxi==1) {
      for (k = 1; k < nyh; k++) {
         koff = 4*nxhd*k;
         k1 = 4*nxhd*ny - koff;
         for (jj = 0; jj < 3; jj++) {
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
         koff = 4*nxhd*k;
         k1 = 4*nxhd*ny - koff;
         for (jj = 0; jj < 3; jj++) {
            t1 = cimagf(f[jj+k1]) + crealf(f[jj+k1])*_Complex_I;
            f[jj+k1] = conjf(f[jj+koff] - t1);
            f[jj+koff] += t1;
         }
      }
   }
#pragma omp parallel for \
private(i,j,k,l,ns,ns2,km,kmr,k1,k2,j1,j2,koff,v_t1,v_t2,v_t3,v_t4,v_t5)
   for (i = nxi-1; i < nxt; i++) {
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         koff = 4*nxhd*k;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = 4*nxhd*k1;
/*          t1 = f[4*i+k1];   */
/*          t2 = f[1+4*i+k1]; */
/*          t3 = f[2+4*i+k1]; */
            v_t1 = _mm_load_ps((float *)&f[4*i+k1]);
            v_t3 = _mm_load_ps((float *)&f[2+4*i+k1]);
/*          f[4*i+k1] = f[4*i+koff];     */
/*          f[1+4*i+k1] = f[1+4*i+koff]; */
/*          f[2+4*i+k1] = f[2+4*i+koff]; */
            v_t2 = _mm_load_ps((float *)&f[4*i+koff]);
            _mm_store_ps((float *)&f[4*i+k1],v_t2);
            v_t2 = _mm_load_ps((float *)&f[2+4*i+koff]);
            _mm_store_ps((float *)&f[2+4*i+k1],v_t2);
/*          f[4*i+koff] = t1;   */
/*          f[1+4*i+koff] = t2; */
/*          f[2+4*i+koff] = t3; */
            _mm_store_ps((float *)&f[4*i+koff],v_t1);
            _mm_store_ps((float *)&f[2+4*i+koff],v_t3);
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
               j1 = 4*nxhd*(j + k1);
               j2 = 4*nxhd*(j + k2);
/*             t1 = conjf(sct[kmr*j]); */
            v_t1 = _mm_loadl_pi(v_t1,(__m64 *)&sct[kmr*j]);
            v_t1 = _mm_movelh_ps(v_t1,v_t1);
            v_t1 = _mm_mul_ps(v_t1,v_n);
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
void csse2wfft2rm3(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxhd,
                   int nyd, int nxhyd, int nxyhd) {
/* wrapper function for 3 2d real to complex ffts, with packed data */
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
      csse2fft2rm3x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                    nxyhd);
/* perform y fft */
      csse2fft2rm3y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                    nxyhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      csse2fft2rm3y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd,
                    nxyhd);
/* perform x fft */
      csse2fft2rm3x(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
                    nxyhd);
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void csse2gbppush23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                       float *qbm, float *dt, float *dtc, float *ek,
                       int *idimp, int *nppmx, int *nx, int *ny,
                       int *mx, int *my, int *nxv, int *nyv, int *mx1,
                       int *mxy1, int *ipbc) {
   csse2gbppush23lt(ppart,fxy,bxy,kpic,*qbm,*dt,*dtc,ek,*idimp,*nppmx,
                    *nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2gbppushf23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                        int *ncl, int *ihole, float *qbm, float *dt,
                        float *dtc, float *ek, int *idimp, int *nppmx,
                        int *nx, int *ny, int *mx, int *my, int *nxv,
                        int *nyv, int *mx1, int *mxy1, int *ntmax,
                        int *irc) {
   csse2gbppushf23lt(ppart,fxy,bxy,kpic,ncl,ihole,*qbm,*dt,*dtc,ek,
                     *idimp,*nppmx,*nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,
                     *ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2grbppush23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                        float *qbm, float *dt, float *dtc, float *ci,
                        float *ek, int *idimp, int *nppmx, int *nx,
                        int *ny, int *mx, int *my, int *nxv, int *nyv,
                        int *mx1, int *mxy1, int *ipbc) {
   csse2grbppush23lt(ppart,fxy,bxy,kpic,*qbm,*dt,*dtc,*ci,ek,*idimp,
                     *nppmx,*nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2grbppushf23lt_(float *ppart, float *fxy, float *bxy,
                         int *kpic, int *ncl, int *ihole, float *qbm,
                         float *dt, float *dtc, float *ci, float *ek,
                         int *idimp, int *nppmx, int *nx, int *ny,
                         int *mx, int *my, int *nxv, int *nyv, int *mx1,
                         int *mxy1, int *ntmax, int *irc) {
   csse2grbppushf23lt(ppart,fxy,bxy,kpic,ncl,ihole,*qbm,*dt,*dtc,*ci,ek,
                      *idimp,*nppmx,*nx,*ny,*mx,*my,*nxv,*nyv,*mx1,
                      *mxy1,*ntmax,irc);
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
void csse2gjppost2lt_(float *ppart, float *cu, int *kpic, float *qm,
                      float *dt, int *nppmx, int *idimp, int *nx,
                      int *ny, int *mx, int *my, int *nxv, int *nyv,
                      int *mx1, int *mxy1, int *ipbc) {
   csse2gjppost2lt(ppart,cu,kpic,*qm,*dt,*nppmx,*idimp,*nx,*ny,*mx,*my,
                   *nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2gjppostf2lt_(float *ppart, float *cu, int *kpic, int *ncl,
                       int *ihole, float *qm, float *dt, int *nppmx,
                       int *idimp, int *nx, int *ny, int *mx, int *my,
                       int *nxv, int *nyv, int *mx1, int *mxy1,
                       int *ntmax, int *irc) {
   csse2gjppostf2lt(ppart,cu,kpic,ncl,ihole,*qm,*dt,*nppmx,*idimp,*nx,
                    *ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2grjppost2lt_(float *ppart, float *cu, int *kpic, float *qm,
                       float *dt, float *ci, int *nppmx, int *idimp,
                       int *nx, int *ny, int *mx, int *my, int *nxv,
                       int *nyv, int *mx1, int *mxy1, int *ipbc) {
   csse2grjppost2lt(ppart,cu,kpic,*qm,*dt,*ci,*nppmx,*idimp,*nx,*ny,*mx,
                    *my,*nxv,*nyv,*mx1,*mxy1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void csse2grjppostf2lt_(float *ppart, float *cu, int *kpic, int *ncl,
                        int *ihole, float *qm, float *dt, float *ci,
                        int *nppmx, int *idimp, int *nx, int *ny,
                        int *mx, int *my, int *nxv, int *nyv, int *mx1,
                        int *mxy1, int *ntmax, int *irc) {
   csse2grjppostf2lt(ppart,cu,kpic,ncl,ihole,*qm,*dt,*ci,*nppmx,*idimp,
                     *nx,*ny,*mx,*my,*nxv,*nyv,*mx1,*mxy1,*ntmax,irc);
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
void csse2mpois23_(float complex *q, float complex *fxy, int *isign,
                   float complex *ffc, float *ax, float *ay,
                   float *affp, float *we, int *nx, int *ny, int *nxvh,
                   int *nyv, int *nxhd, int *nyhd) {
   csse2mpois23(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*nxvh,*nyv,
                *nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csse2mcuperp2_(float complex *cu, int *nx, int *ny, int *nxvh,
                    int *nyv) {
   csse2mcuperp2(cu,*nx,*ny,*nxvh,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void csse2mibpois23_(float complex *cu, float complex *bxy,
                     float complex *ffc, float *ci, float *wm, int *nx,
                     int *ny, int *nxvh, int *nyv, int *nxhd,
                     int *nyhd) {
   csse2mibpois23(cu,bxy,ffc,*ci,wm,*nx,*ny,*nxvh,*nyv,*nxhd,*nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csse2mmaxwel2_(float complex *exy, float complex *bxy,
                    float complex *cu, float complex *ffc, float *ci,
                    float *dt, float *wf, float *wm, int *nx, int *ny,
                    int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   csse2mmaxwel2(exy,bxy,cu,ffc,*ci,*dt,wf,wm,*nx,*ny,*nxvh,*nyv,*nxhd,
                 *nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csse2memfield2_(float complex *fxy, float complex *exy,
                     float complex *ffc, int *isign, int *nx, int *ny,
                     int *nxvh, int *nyv, int *nxhd, int *nyhd) {
   csse2memfield2(fxy,exy,ffc,*isign,*nx,*ny,*nxvh,*nyv,*nxhd,*nyhd);
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
void csse2wfft2rm3_(float complex *f, int *isign, int *mixup,
                    float complex *sct, int *indx, int *indy, int *nxhd,
                    int *nyd, int *nxhyd, int *nxyhd) {
   csse2wfft2rm3(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,
                 *nxyhd);
   return;
}

