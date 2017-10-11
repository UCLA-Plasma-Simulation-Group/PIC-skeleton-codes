/* KNC C Library for Skeleton 3D Electromagnetic Vector PIC Code */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <immintrin.h>
#include "kncmbpush3.h"

/*--------------------------------------------------------------------*/
void ckncgbppush3lt(float ppart[], float fxyz[], float bxyz[],
                    int kpic[], float qbm, float dt, float dtc,
                    float *ek, int idimp, int nppmx, int nx, int ny,
                    int nz, int mx, int my, int mz, int nxv, int nyv,
                    int nzv, int mx1, int my1, int mxyz1,int ipbc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field.  Using the Boris Mover.
   OpenMP/vector version using guard cells
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
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
   requires KNC, ppart needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
   fxyz needs to have 4 components, although one is not used
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, nn, mm, ll, nm, mxv, myv, mxyv, nxyv;
   float qtmh, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, z, vx, vy, vz;
   double sum1, sum2;
   __m512i v_noff, v_moff, v_loff, v_mxv4, v_mxyv4;
   __m512i v_nn, v_mm, v_ll, v_nm, v_it, v_perm;
   __m512 v_qtmh, v_dt, v_dtc, v_one, v_zero;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_at, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s;
   __m512 v_two, v_half, v_ox, v_oy, v_oz;
   __m512d v_sum1, v_d;
   __mmask16 msk;
   __attribute__((aligned(64))) unsigned int kk[16];
   __attribute__((aligned(64))) double dd[8];
   __attribute__((aligned(64))) float sfxyz[4*MXV*MYV*MZV];
   __attribute__((aligned(64))) float sbxyz[4*MXV*MYV*MZV];
/* __attribute__((aligned(64))) float sfxyz[4*(mx+1)*(my+1)*(mz+1)]; */
/* __attribute__((aligned(64))) float sbxyz[4*(mx+1)*(my+1)*(mz+1)]; */
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
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
   v_mxv4 = _mm512_set1_epi32(4*mxv);
   v_mxyv4 = _mm512_set1_epi32(4*mxyv);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qtmh = _mm512_set1_ps(qtmh);
   v_dt = _mm512_set1_ps(dt);
   v_dtc = _mm512_set1_ps(dtc);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_two = _mm512_set1_ps(2.0f);
   v_half = _mm512_set1_ps(0.5f);
   v_edgelx = _mm512_set1_ps(edgelx);
   v_edgely = _mm512_set1_ps(edgely);
   v_edgelz = _mm512_set1_ps(edgelz);
   v_edgerx = _mm512_set1_ps(edgerx);
   v_edgery = _mm512_set1_ps(edgery);
   v_edgerz = _mm512_set1_ps(edgerz);
   v_sum1 = _mm512_set1_pd(0.0);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,nps,nn,mm,ll,nm,x,y,z,vx, \
vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt, \
omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,sum1, \
v_noff,v_moff,v_loff,v_nn,v_mm,v_ll,v_nm,v_it,v_x,v_y,v_z,v_dxp,v_dyp, \
v_dzp,v_amx,v_amy,v_amz,v_dx1,v_dx,v_dy,v_dz,v_vx,v_vy,v_vz,v_ox,v_oy, \
v_oz,v_at,v_d,v_sum1,a,b,c,d,e,f,g,h,p,q,r,s,msk,kk,dd,sfxyz,sbxyz) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      ll = (mz < nz-loff ? mz : nz-loff) + 1;
      nps = 4*(nn/4);
/* load electric field */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 0; i < nn; i++) {                          */
/*             sfxyz[4*(i+mxv*j+mxyv*k)]                        */
/*             = fxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];   */
/*             sfxyz[1+4*(i+mxv*j+mxyv*k)]                      */
/*             = fxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*             sfxyz[2+4*(i+mxv*j+mxyv*k)]                      */
/*             = fxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*          }                                                   */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&fxyz[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&fxyz[m+16]);
               m = 4*(i + mxv*j + mxyv*k);
               _mm512_packstorelo_ps(&sfxyz[m],v_at);
               _mm512_packstorehi_ps(&sfxyz[m+16],v_at);
            }
/* loop over remaining elements */
            for (i = nps; i < nn; i++) {
               sfxyz[4*(i+mxv*j+mxyv*k)]
               = fxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+4*(i+mxv*j+mxyv*k)]
               = fxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+4*(i+mxv*j+mxyv*k)]
               = fxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[3+4*(i+mxv*j+mxyv*k)]
               = fxyz[3+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
/* load magnetic field */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 0; i < nn; i++) {                          */
/*             sbxyz[4*(i+mxv*j+mxyv*k)]                        */
/*             = bxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];   */
/*             sbxyz[1+4*(i+mxv*j+mxyv*k)]                      */
/*             = bxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*             sbxyz[2+4*(i+mxv*j+mxyv*k)]                      */
/*             = bxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*          }                                                   */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&bxyz[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&bxyz[m+16]);
               m = 4*(i + mxv*j + mxyv*k);
               _mm512_packstorelo_ps(&sbxyz[m],v_at);
               _mm512_packstorehi_ps(&sbxyz[m+16],v_at);
            }
/* loop over remaining elements */
            for (i = nps; i < nn; i++) {
               sbxyz[4*(i+mxv*j+mxyv*k)]
               = bxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[1+4*(i+mxv*j+mxyv*k)]
               = bxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[2+4*(i+mxv*j+mxyv*k)]
               = bxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[3+4*(i+mxv*j+mxyv*k)]
               = bxyz[3+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
      nps = 16*(npp/16);
      sum1 = 0.0;
      v_sum1 = _mm512_set1_pd(0.0);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = x - (float) nn; */
/*       dyp = y - (float) mm; */
/*       dzp = z - (float) ll; */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_sub_ps(v_x,v_dxp);
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*       nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff)); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv4,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv4,v_mm));
         v_nm = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*       amx = 1.0f - dxp; */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_one,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/* find electric field */
/*       nn = nm; */
         _mm512_store_epi32(kk,v_nm);
/* load sfxyz[nn:nn+3] and sfxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[nn:nn+3] field components */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[nn+4:nn+7] field components */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of electric field */
/*       dx = amx*sfxyz[nn] + amy*sfxyz[nn+4];                     */
         v_dx = _mm512_mul_ps(v_amx,a);
         v_dx = _mm512_fmadd_ps(v_amy,p,v_dx);
/*       dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];                 */
         v_dy = _mm512_mul_ps(v_amx,b);
         v_dy = _mm512_fmadd_ps(v_amy,q,v_dy);
/*       dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];                 */
         v_dz = _mm512_mul_ps(v_amx,c);
         v_dz = _mm512_fmadd_ps(v_amy,r,v_dz);
/*       mm = nn + 4*mxv;                                     */
/* load sfxyz[mm:mm+3] and sfxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxv;                                     */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxv;                                       */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of electric field */
/*       dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);          */
         v_dx = _mm512_fmadd_ps(v_dyp,a,v_dx);
         v_dx = _mm512_fmadd_ps(v_dx1,p,v_dx);
         v_dx = _mm512_mul_ps(v_amz,v_dx);
/*       dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);      */
         v_dy = _mm512_fmadd_ps(v_dyp,b,v_dy);
         v_dy = _mm512_fmadd_ps(v_dx1,q,v_dy);
         v_dy = _mm512_mul_ps(v_amz,v_dy);
/*       dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);      */
         v_dz = _mm512_fmadd_ps(v_dyp,c,v_dz);
         v_dz = _mm512_fmadd_ps(v_dx1,r,v_dz);
         v_dz = _mm512_mul_ps(v_amz,v_dz);
/*       nn += 4*mxyv;                                           */
         v_nn = _mm512_add_epi32(v_nm,v_mxyv4);
         _mm512_store_epi32(kk,v_nn);
/* load sfxyz[nn:nn+3] and sfxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[nn:nn+3] field components */
/* where nn = nn + 4*mxyv;                                    */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*mxyv;                                      */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of electric field */
/*       vx = amx*sfxyz[nn] + amy*sfxyz[nn+4];                     */
         v_vx = _mm512_mul_ps(v_amx,a);
         v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*       vy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];                 */
         v_vy = _mm512_mul_ps(v_amx,b);
         v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*       vz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];                 */
         v_vz = _mm512_mul_ps(v_amx,c);
         v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*       mm = nn + 4*mxv;                                        */
/* load sfxyz[mm:mm+3] and sfxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
              &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                            */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                              */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of electric field */
/*       dx = dx + dzp*(vx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);     */
         v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
         v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
         v_dx = _mm512_fmadd_ps(v_dzp,v_vx,v_dx);
/*       dy = dy + dzp*(vy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]); */
         v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
         v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
         v_dy = _mm512_fmadd_ps(v_dzp,v_vy,v_dy);
/*       dz = dz + dzp*(vz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]); */
         v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
         v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
         v_dz = _mm512_fmadd_ps(v_dzp,v_vz,v_dz);
/* find magnetic field */
/*       nn = nm; */
         _mm512_store_epi32(kk,v_nm);
/* load sbxyz[nn:nn+3] and sbxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[nn:nn+3] field components */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[nn+4:nn+7] field components */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of magnetic field */
/*       ox = amx*sbxyz[nn] + amy*sbxyz[nn+4];                     */
         v_ox = _mm512_mul_ps(v_amx,a);
         v_ox = _mm512_fmadd_ps(v_amy,p,v_ox);
/*       oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];                 */
         v_oy = _mm512_mul_ps(v_amx,b);
         v_oy = _mm512_fmadd_ps(v_amy,q,v_oy);
/*       oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];                 */
         v_oz = _mm512_mul_ps(v_amx,c);
         v_oz = _mm512_fmadd_ps(v_amy,r,v_oz);
/*       mm = nn + 4*mxv;                                     */
/* load sbxyz[mm:mm+3] and sbxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxv;                                     */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxv;                                       */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of magnetic field */
/*       ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);          */
         v_ox = _mm512_fmadd_ps(v_dyp,a,v_ox);
         v_ox = _mm512_fmadd_ps(v_dx1,p,v_ox);
         v_ox = _mm512_mul_ps(v_amz,v_ox);
/*       oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);      */
         v_oy = _mm512_fmadd_ps(v_dyp,b,v_oy);
         v_oy = _mm512_fmadd_ps(v_dx1,q,v_oy);
         v_oy = _mm512_mul_ps(v_amz,v_oy);
/*       oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);      */
         v_oz = _mm512_fmadd_ps(v_dyp,c,v_oz);
         v_oz = _mm512_fmadd_ps(v_dx1,r,v_oz);
         v_oz = _mm512_mul_ps(v_amz,v_oz);
/*       nn += 4*mxyv;                                           */
         v_nn = _mm512_add_epi32(v_nm,v_mxyv4);
         _mm512_store_epi32(kk,v_nn);
/* load sbxyz[nn:nn+3] and sbxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[nn:nn+3] field components */
/* where nn = nn + 4*mxyv;                                    */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*mxyv;                                     */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of magnetic field */
/*       vx = amx*sbxyz[nn] + amy*sbxyz[nn+4];                     */
         v_vx = _mm512_mul_ps(v_amx,a);
         v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*       vy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];                 */
         v_vy = _mm512_mul_ps(v_amx,b);
         v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*       vz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];                 */
         v_vz = _mm512_mul_ps(v_amx,c);
         v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*       mm = nn + 4*mxv;                                        */
/* load sbxyz[mm:mm+3] and sbxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                            */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                           */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of magnetic field */
/*       ox = ox + dzp*(vx + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);     */
         v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
         v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
         v_ox = _mm512_fmadd_ps(v_dzp,v_vx,v_ox);
/*       oy = oy + dzp*(vy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]); */
         v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
         v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
         v_oy = _mm512_fmadd_ps(v_dzp,v_vy,v_oy);
/*       oz = oz + dzp*(vz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]); */
         v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
         v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
         v_oz = _mm512_fmadd_ps(v_dzp,v_vz,v_oz);
/* calculate half impulse */
/*       dx *= qtmh; */
/*       dy *= qtmh; */
/*       dz *= qtmh; */
         v_dx = _mm512_mul_ps(v_dx,v_qtmh);
         v_dy = _mm512_mul_ps(v_dy,v_qtmh);
         v_dz = _mm512_mul_ps(v_dz,v_qtmh);
/* half acceleration */
/*       acx = ppart[j+3*nppmx+npoff] + dx; */
/*       acy = ppart[j+4*nppmx+npoff] + dy; */
/*       acz = ppart[j+5*nppmx+npoff] + dz; */
         a = _mm512_add_ps(v_dx,_mm512_load_ps(&ppart[j+3*nppmx+npoff]));
         b = _mm512_add_ps(v_dy,_mm512_load_ps(&ppart[j+4*nppmx+npoff]));
         c = _mm512_add_ps(v_dz,_mm512_load_ps(&ppart[j+5*nppmx+npoff]));
/* time-centered kinetic energy */
/*       sum1 += (acx*acx + acy*acy + acz*acz); */
         v_at = _mm512_fmadd_ps(b,b,_mm512_mul_ps(a,a));
         v_at = _mm512_fmadd_ps(c,c,v_at);
/* convert to double precision before accumulating */
         v_sum1 = _mm512_add_pd(v_sum1,_mm512_cvtpslo_pd(v_at));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_at,78));
         v_sum1 = _mm512_add_pd(v_sum1,v_d);
/* calculate cyclotron frequency */
/*       omxt = qtmh*ox; */
/*       omyt = qtmh*oy; */
/*       omzt = qtmh*oz; */
         e = _mm512_mul_ps(v_qtmh,v_ox);
         f = _mm512_mul_ps(v_qtmh,v_oy);
         g = _mm512_mul_ps(v_qtmh,v_oz);
/* calculate rotation matrix */
/*       vx = omxt*omxt; */
         v_vx = _mm512_mul_ps(e,e);
/*       vy = omyt*omyt; */
         v_vy = _mm512_mul_ps(f,f);
/*       vz = omzt*omzt; */
         v_vz = _mm512_mul_ps(g,g);
/*       omt = omxt*omxt + omyt*omyt + omzt*omzt; */
         v_at = _mm512_add_ps(_mm512_add_ps(v_vx,v_vy),v_vz);
/*       anorm = 2.0f/(1.0f + omt); */
         d = _mm512_div_ps(v_two,_mm512_add_ps(v_one,v_at));
/*       omt = 0.5f*(1.0f - omt); */
         h = _mm512_mul_ps(v_half,_mm512_sub_ps(v_one,v_at));
/*       vx = (omt + vx)*acx; */
         v_vx = _mm512_mul_ps(_mm512_add_ps(h,v_vx),a);
/*       vy = (omt + vy)*acy; */
         v_vy = _mm512_mul_ps(_mm512_add_ps(h,v_vy),b);
/*       vz = (omt + vz)*acz; */
         v_vz = _mm512_mul_ps(_mm512_add_ps(h,v_vz),c);
/*       omt = omxt*omyt; */
         h = _mm512_mul_ps(e,f);
/*       vx = vx + (omzt + omt)*acy; */
         v_vx = _mm512_fmadd_ps(_mm512_add_ps(h,g),b,v_vx);
/*       vy = vy + (omt - omzt)*acx; */
         v_vy = _mm512_fmadd_ps(_mm512_sub_ps(h,g),a,v_vy);
/*       omt = omxt*omzt;  */
         h = _mm512_mul_ps(e,g);
/*       vx = vx + (omt - omyt)*acz; */
         v_vx = _mm512_fmadd_ps(_mm512_sub_ps(h,f),c,v_vx);
/*       vz = vz + (omt + omyt)*acx; */
         v_vz = _mm512_fmadd_ps(_mm512_add_ps(h,f),a,v_vz);
/*       omt = omyt*omzt; */
         h = _mm512_mul_ps(f,g);
/*       vy = vy + (omt + omxt)*acz; */
         v_vy = _mm512_fmadd_ps(_mm512_add_ps(h,e),c,v_vy);
/*       vz = vz + (omt - omxt)*acy; */
         v_vz = _mm512_fmadd_ps(_mm512_sub_ps(h,e),b,v_vz);
/* new velocity */
/*       vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm; */
/*       vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm; */
/*       vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm; */
         v_vx = _mm512_fmadd_ps(v_vx,d,v_dx);
         v_vy = _mm512_fmadd_ps(v_vy,d,v_dy);
         v_vz = _mm512_fmadd_ps(v_vz,d,v_dz);
/* new position */
/*       dx = x + vx*dtc; */
/*       dy = y + vy*dtc; */
/*       dz = z + vz*dtc; */
         v_dx = _mm512_fmadd_ps(v_vx,v_dtc,v_x);
         v_dy = _mm512_fmadd_ps(v_vy,v_dtc,v_y);
         v_dz = _mm512_fmadd_ps(v_vz,v_dtc,v_z);
/* reflecting boundary conditions */
         if (ipbc==2) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             vx = -vx;                           */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             vy = -vy;                           */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/*          if ((dz < edgelz) || (dz >= edgerz)) { */
/*             dz = z;                             */
/*             vz = -vz;                           */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dz,v_edgerz,
                  _MM_CMPINT_GE));
            v_dz = _mm512_mask_blend_ps(msk,v_dz,v_z);
            v_vz = _mm512_mask_sub_ps(v_vz,msk,v_zero,v_vz);
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             vx = -vx;                           */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             vy = -vy;                           */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy;   */
/*       ppart[j+2*nppmx+npoff] = dz; */
         _mm512_store_ps(&ppart[j+npoff],v_dx);
         _mm512_store_ps(&ppart[j+nppmx+npoff],v_dy);
         _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_dz);
/* set new velocity */
/*       ppart[j+3*nppmx+npoff] = vx; */
/*       ppart[j+4*nppmx+npoff] = vy; */
/*       ppart[j+5*nppmx+npoff] = vz; */
         _mm512_store_ps(&ppart[j+3*nppmx+npoff],v_vx);
         _mm512_store_ps(&ppart[j+4*nppmx+npoff],v_vy);
         _mm512_store_ps(&ppart[j+5*nppmx+npoff],v_vz);
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nm = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find electric field */
         nn = nm;
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+4];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];
         mm = nn + 4*mxv;
         dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);
         dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);
         dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);
         nn += 4*mxyv;
         acx = amx*sfxyz[nn] + amy*sfxyz[nn+4];
         acy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];
         acz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];
         mm = nn + 4*mxv;
         dx = dx + dzp*(acx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);
         dy = dy + dzp*(acy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);
         dz = dz + dzp*(acz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxyz[nn] + amy*sbxyz[nn+4];
         oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];
         oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];
         mm = nn + 4*mxv;
         ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);
         oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);
         oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);
         nn += 4*mxyv;
         acx = amx*sbxyz[nn] + amy*sbxyz[nn+4];
         acy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];
         acz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];
         mm = nn + 4*mxv;
         ox = ox + dzp*(acx + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);
         oy = oy + dzp*(acy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);
         oz = oz + dzp*(acz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+3*nppmx+npoff] + dx;
         acy = ppart[j+4*nppmx+npoff] + dy;
         acz = ppart[j+5*nppmx+npoff] + dz;
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
         dz = z + vz*dtc;
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
            if ((dz < edgelz) || (dz >= edgerz)) {
               dz = z;
               vz = -vz;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               vx = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               vy = -vy;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
/* set new velocity */
         ppart[j+3*nppmx+npoff] = vx;
         ppart[j+4*nppmx+npoff] = vy;
         ppart[j+5*nppmx+npoff] = vz;
      }
/*    sum2 += sum1; */
      _mm512_store_pd(&dd[0],v_sum1);
      for (j = 1; j < 8; j++) {
         dd[0] += dd[j];
      }
      sum2 += (sum1 + dd[0]);
   }
/* normalize kinetic energy */
   *ek += 0.5f*sum2;
   return;
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void ckncgbppushf3lt(float ppart[], float fxyz[], float bxyz[],
                     int kpic[], int ncl[], int ihole[], float qbm,
                     float dt, float dtc, float *ek, int idimp,
                     int nppmx, int nx, int ny, int nz, int mx, int my,
                     int mz, int nxv, int nyv, int nzv, int mx1,
                     int my1, int mxyz1, int ntmax, int *irc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field.  Using the Boris Mover.
   also determines list of particles which are leaving this tile
   OpenMP/vector version using guard cells
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
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
   requires KNC, ppart needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
   fxyz needs to have 4 components, although one is not used
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, ii, ih, nh, nn, mm, ll, nm, mxv, myv, mxyv, nxyv;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float qtmh, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, z, vx, vy, vz;
   double sum1, sum2;
   __m512i v_noff, v_moff, v_loff, v_mxv4, v_mxyv4;
   __m512i v_nn, v_mm, v_ll, v_nm, v_it, v_0, v_1, v_3, v_9, v_perm;
   __m512 v_qtmh, v_dt, v_dtc, v_one, v_zero, v_anx, v_any, v_anz;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_at, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s;
   __m512 v_two, v_half, v_ox, v_oy, v_oz;
   __m512d v_sum1, v_d;
   __mmask16 msk1, msk2;
   __attribute__((aligned(64))) unsigned int kk[16];
   __attribute__((aligned(64))) double dd[8];
   __attribute__((aligned(64))) float sfxyz[4*MXV*MYV*MZV];
   __attribute__((aligned(64))) float sbxyz[4*MXV*MYV*MZV];
/* __attribute__((aligned(64))) float sfxyz[4*(mx+1)*(my+1)*(mz+1)]; */
/* __attribute__((aligned(64))) float sbxyz[4*(mx+1)*(my+1)*(mz+1)]; */
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   qtmh = 0.5f*qbm*dt;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
   sum2 = 0.0;
/* set boundary values */
   v_mxv4 = _mm512_set1_epi32(4*mxv);
   v_mxyv4 = _mm512_set1_epi32(4*mxyv);
   v_0 = _mm512_set1_epi32(0);
   v_1 = _mm512_set1_epi32(1);
   v_3 = _mm512_set1_epi32(3);
   v_9 = _mm512_set1_epi32(9);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qtmh = _mm512_set1_ps(qtmh);
   v_dt = _mm512_set1_ps(dt);
   v_dtc = _mm512_set1_ps(dtc);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_two = _mm512_set1_ps(2.0f);
   v_half = _mm512_set1_ps(0.5f);
   v_anx = _mm512_set1_ps(anx);
   v_any = _mm512_set1_ps(any);
   v_anz = _mm512_set1_ps(anz);
   v_sum1 = _mm512_set1_pd(0.0);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,ii,noff,moff,loff,npp,npoff,nps,nn,mm,ll,nm,ih,nh,x, \
y,z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx,acy,acz, \
omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9, \
edgelx,edgely,edgelz,edgerx,edgery,edgerz,sum1,v_noff,v_moff,v_loff, \
v_nn,v_mm,v_ll,v_nm,v_it,v_x,v_y,v_z,v_dxp,v_dyp,v_dzp,v_amx,v_amy, \
v_amz,v_dx1,v_dx,v_dy,v_dz,v_vx,v_vy,v_vz,v_ox,v_oy,v_oz,v_at,v_edgelx, \
v_edgely,v_edgelz,v_edgerx,v_edgery,v_edgerz,v_d,v_sum1,a,b,c,d,e,f,g, \
h,p,q,r,s,msk1,msk2,kk,dd,sfxyz,sbxyz) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
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
      v_edgelx = _mm512_set1_ps(edgelx);
      v_edgely = _mm512_set1_ps(edgely);
      v_edgelz = _mm512_set1_ps(edgelz);
      v_edgerx = _mm512_set1_ps(edgerx);
      v_edgery = _mm512_set1_ps(edgery);
      v_edgerz = _mm512_set1_ps(edgerz);
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
      ll += 1;
/* load local fields from global array */
      nps = 4*(nn/4);
/* load electric field */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 0; i < nn; i++) {                          */
/*             sfxyz[4*(i+mxv*j+mxyv*k)]                        */
/*             = fxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];   */
/*             sfxyz[1+4*(i+mxv*j+mxyv*k)]                      */
/*             = fxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*             sfxyz[2+4*(i+mxv*j+mxyv*k)]                      */
/*             = fxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*          }                                                   */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&fxyz[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&fxyz[m+16]);
               m = 4*(i + mxv*j + mxyv*k);
               _mm512_packstorelo_ps(&sfxyz[m],v_at);
               _mm512_packstorehi_ps(&sfxyz[m+16],v_at);
            }
/* loop over remaining elements */
            for (i = nps; i < nn; i++) {
               sfxyz[4*(i+mxv*j+mxyv*k)]
               = fxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+4*(i+mxv*j+mxyv*k)]
               = fxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+4*(i+mxv*j+mxyv*k)]
               = fxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[3+4*(i+mxv*j+mxyv*k)]
               = fxyz[3+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
/* load magnetic field */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 0; i < nn; i++) {                          */
/*             sbxyz[4*(i+mxv*j+mxyv*k)]                        */
/*             = bxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];   */
/*             sbxyz[1+4*(i+mxv*j+mxyv*k)]                      */
/*             = bxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*             sbxyz[2+4*(i+mxv*j+mxyv*k)]                      */
/*             = bxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*          }                                                   */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&bxyz[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&bxyz[m+16]);
               m = 4*(i + mxv*j + mxyv*k);
               _mm512_packstorelo_ps(&sbxyz[m],v_at);
               _mm512_packstorehi_ps(&sbxyz[m+16],v_at);
            }
/* loop over remaining elements */
            for (i = nps; i < nn; i++) {
               sbxyz[4*(i+mxv*j+mxyv*k)]
               = bxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[1+4*(i+mxv*j+mxyv*k)]
               = bxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[2+4*(i+mxv*j+mxyv*k)]
               = bxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[3+4*(i+mxv*j+mxyv*k)]
               = bxyz[3+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
/* clear counters */
/*    for (j = 0; j < 26; j++) { */
/*       ncl[j+26*l] = 0;        */
/*    }                         */
      memset((void*)&ncl[26*l],0,26*sizeof(int));
      nps = 16*(npp/16);
      sum1 = 0.0;
      v_sum1 = _mm512_set1_pd(0.0);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = x - (float) nn; */
/*       dyp = y - (float) mm; */
/*       dzp = z - (float) ll; */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_sub_ps(v_x,v_dxp);
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*       nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff)); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv4,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv4,v_mm));
         v_nm = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*       amx = 1.0f - dxp; */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_one,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/* find electric field */
/*       nn = nm; */
         _mm512_store_epi32(kk,v_nm);
/* load sfxyz[nn:nn+3] and sfxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[nn:nn+3] field components */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[nn+4:nn+7] field components */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of electric field */
/*       dx = amx*sfxyz[nn] + amy*sfxyz[nn+4];                     */
         v_dx = _mm512_mul_ps(v_amx,a);
         v_dx = _mm512_fmadd_ps(v_amy,p,v_dx);
/*       dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];                 */
         v_dy = _mm512_mul_ps(v_amx,b);
         v_dy = _mm512_fmadd_ps(v_amy,q,v_dy);
/*       dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];                 */
         v_dz = _mm512_mul_ps(v_amx,c);
         v_dz = _mm512_fmadd_ps(v_amy,r,v_dz);
/*       mm = nn + 4*mxv;                                     */
/* load sfxyz[mm:mm+3] and sfxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxv;                                     */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxv;                                       */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of electric field */
/*       dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);          */
         v_dx = _mm512_fmadd_ps(v_dyp,a,v_dx);
         v_dx = _mm512_fmadd_ps(v_dx1,p,v_dx);
         v_dx = _mm512_mul_ps(v_amz,v_dx);
/*       dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);      */
         v_dy = _mm512_fmadd_ps(v_dyp,b,v_dy);
         v_dy = _mm512_fmadd_ps(v_dx1,q,v_dy);
         v_dy = _mm512_mul_ps(v_amz,v_dy);
/*       dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);      */
         v_dz = _mm512_fmadd_ps(v_dyp,c,v_dz);
         v_dz = _mm512_fmadd_ps(v_dx1,r,v_dz);
         v_dz = _mm512_mul_ps(v_amz,v_dz);
/*       nn += 4*mxyv;                                           */
         v_nn = _mm512_add_epi32(v_nm,v_mxyv4);
         _mm512_store_epi32(kk,v_nn);
/* load sfxyz[nn:nn+3] and sfxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[nn:nn+3] field components */
/* where nn = nn + 4*mxyv;                                    */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*mxyv;                                      */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of electric field */
/*       vx = amx*sfxyz[nn] + amy*sfxyz[nn+4];                     */
         v_vx = _mm512_mul_ps(v_amx,a);
         v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*       vy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];                 */
         v_vy = _mm512_mul_ps(v_amx,b);
         v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*       vz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];                 */
         v_vz = _mm512_mul_ps(v_amx,c);
         v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*       mm = nn + 4*mxv;                                        */
/* load sfxyz[mm:mm+3] and sfxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
              &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                            */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                              */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of electric field */
/*       dx = dx + dzp*(vx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);     */
         v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
         v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
         v_dx = _mm512_fmadd_ps(v_dzp,v_vx,v_dx);
/*       dy = dy + dzp*(vy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]); */
         v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
         v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
         v_dy = _mm512_fmadd_ps(v_dzp,v_vy,v_dy);
/*       dz = dz + dzp*(vz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]); */
         v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
         v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
         v_dz = _mm512_fmadd_ps(v_dzp,v_vz,v_dz);
/* find magnetic field */
/*       nn = nm; */
         _mm512_store_epi32(kk,v_nm);
/* load sbxyz[nn:nn+3] and sbxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[nn:nn+3] field components */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[nn+4:nn+7] field components */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of magnetic field */
/*       ox = amx*sbxyz[nn] + amy*sbxyz[nn+4];                     */
         v_ox = _mm512_mul_ps(v_amx,a);
         v_ox = _mm512_fmadd_ps(v_amy,p,v_ox);
/*       oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];                 */
         v_oy = _mm512_mul_ps(v_amx,b);
         v_oy = _mm512_fmadd_ps(v_amy,q,v_oy);
/*       oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];                 */
         v_oz = _mm512_mul_ps(v_amx,c);
         v_oz = _mm512_fmadd_ps(v_amy,r,v_oz);
/*       mm = nn + 4*mxv;                                     */
/* load sbxyz[mm:mm+3] and sbxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxv;                                     */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxv;                                       */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of magnetic field */
/*       ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);          */
         v_ox = _mm512_fmadd_ps(v_dyp,a,v_ox);
         v_ox = _mm512_fmadd_ps(v_dx1,p,v_ox);
         v_ox = _mm512_mul_ps(v_amz,v_ox);
/*       oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);      */
         v_oy = _mm512_fmadd_ps(v_dyp,b,v_oy);
         v_oy = _mm512_fmadd_ps(v_dx1,q,v_oy);
         v_oy = _mm512_mul_ps(v_amz,v_oy);
/*       oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);      */
         v_oz = _mm512_fmadd_ps(v_dyp,c,v_oz);
         v_oz = _mm512_fmadd_ps(v_dx1,r,v_oz);
         v_oz = _mm512_mul_ps(v_amz,v_oz);
/*       nn += 4*mxyv;                                           */
         v_nn = _mm512_add_epi32(v_nm,v_mxyv4);
         _mm512_store_epi32(kk,v_nn);
/* load sbxyz[nn:nn+3] and sbxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[nn:nn+3] field components */
/* where nn = nn + 4*mxyv;                                    */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*mxyv;                                     */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of magnetic field */
/*       vx = amx*sbxyz[nn] + amy*sbxyz[nn+4];                     */
         v_vx = _mm512_mul_ps(v_amx,a);
         v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*       vy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];                 */
         v_vy = _mm512_mul_ps(v_amx,b);
         v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*       vz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];                 */
         v_vz = _mm512_mul_ps(v_amx,c);
         v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*       mm = nn + 4*mxv;                                        */
/* load sbxyz[mm:mm+3] and sbxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                            */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                           */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of magnetic field */
/*       ox = ox + dzp*(vx + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);     */
         v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
         v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
         v_ox = _mm512_fmadd_ps(v_dzp,v_vx,v_ox);
/*       oy = oy + dzp*(vy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]); */
         v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
         v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
         v_oy = _mm512_fmadd_ps(v_dzp,v_vy,v_oy);
/*       oz = oz + dzp*(vz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]); */
         v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
         v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
         v_oz = _mm512_fmadd_ps(v_dzp,v_vz,v_oz);
/* calculate half impulse */
/*       dx *= qtmh; */
/*       dy *= qtmh; */
/*       dz *= qtmh; */
         v_dx = _mm512_mul_ps(v_dx,v_qtmh);
         v_dy = _mm512_mul_ps(v_dy,v_qtmh);
         v_dz = _mm512_mul_ps(v_dz,v_qtmh);
/* half acceleration */
/*       acx = ppart[j+3*nppmx+npoff] + dx; */
/*       acy = ppart[j+4*nppmx+npoff] + dy; */
/*       acz = ppart[j+5*nppmx+npoff] + dz; */
         a = _mm512_add_ps(v_dx,_mm512_load_ps(&ppart[j+3*nppmx+npoff]));
         b = _mm512_add_ps(v_dy,_mm512_load_ps(&ppart[j+4*nppmx+npoff]));
         c = _mm512_add_ps(v_dz,_mm512_load_ps(&ppart[j+5*nppmx+npoff]));
/* time-centered kinetic energy */
/*       sum1 += (acx*acx + acy*acy + acz*acz); */
         v_at = _mm512_fmadd_ps(b,b,_mm512_mul_ps(a,a));
         v_at = _mm512_fmadd_ps(c,c,v_at);
/* convert to double precision before accumulating */
         v_sum1 = _mm512_add_pd(v_sum1,_mm512_cvtpslo_pd(v_at));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_at,78));
         v_sum1 = _mm512_add_pd(v_sum1,v_d);
/* calculate cyclotron frequency */
/*       omxt = qtmh*ox; */
/*       omyt = qtmh*oy; */
/*       omzt = qtmh*oz; */
         e = _mm512_mul_ps(v_qtmh,v_ox);
         f = _mm512_mul_ps(v_qtmh,v_oy);
         g = _mm512_mul_ps(v_qtmh,v_oz);
/* calculate rotation matrix */
/*       vx = omxt*omxt; */
         v_vx = _mm512_mul_ps(e,e);
/*       vy = omyt*omyt; */
         v_vy = _mm512_mul_ps(f,f);
/*       vz = omzt*omzt; */
         v_vz = _mm512_mul_ps(g,g);
/*       omt = omxt*omxt + omyt*omyt + omzt*omzt; */
         v_at = _mm512_add_ps(_mm512_add_ps(v_vx,v_vy),v_vz);
/*       anorm = 2.0f/(1.0f + omt); */
         d = _mm512_div_ps(v_two,_mm512_add_ps(v_one,v_at));
/*       omt = 0.5f*(1.0f - omt); */
         h = _mm512_mul_ps(v_half,_mm512_sub_ps(v_one,v_at));
/*       vx = (omt + vx)*acx; */
         v_vx = _mm512_mul_ps(_mm512_add_ps(h,v_vx),a);
/*       vy = (omt + vy)*acy; */
         v_vy = _mm512_mul_ps(_mm512_add_ps(h,v_vy),b);
/*       vz = (omt + vz)*acz; */
         v_vz = _mm512_mul_ps(_mm512_add_ps(h,v_vz),c);
/*       omt = omxt*omyt; */
         h = _mm512_mul_ps(e,f);
/*       vx = vx + (omzt + omt)*acy; */
         v_vx = _mm512_fmadd_ps(_mm512_add_ps(h,g),b,v_vx);
/*       vy = vy + (omt - omzt)*acx; */
         v_vy = _mm512_fmadd_ps(_mm512_sub_ps(h,g),a,v_vy);
/*       omt = omxt*omzt;  */
         h = _mm512_mul_ps(e,g);
/*       vx = vx + (omt - omyt)*acz; */
         v_vx = _mm512_fmadd_ps(_mm512_sub_ps(h,f),c,v_vx);
/*       vz = vz + (omt + omyt)*acx; */
         v_vz = _mm512_fmadd_ps(_mm512_add_ps(h,f),a,v_vz);
/*       omt = omyt*omzt; */
         h = _mm512_mul_ps(f,g);
/*       vy = vy + (omt + omxt)*acz; */
         v_vy = _mm512_fmadd_ps(_mm512_add_ps(h,e),c,v_vy);
/*       vz = vz + (omt - omxt)*acy; */
         v_vz = _mm512_fmadd_ps(_mm512_sub_ps(h,e),b,v_vz);
/* new velocity */
/*       vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm; */
/*       vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm; */
/*       vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm; */
         v_vx = _mm512_fmadd_ps(v_vx,d,v_dx);
         v_vy = _mm512_fmadd_ps(v_vy,d,v_dy);
         v_vz = _mm512_fmadd_ps(v_vz,d,v_dz);
/* new position */
/*       dx = x + vx*dtc; */
/*       dy = y + vy*dtc; */
/*       dz = z + vz*dtc; */
         v_dx = _mm512_fmadd_ps(v_vx,v_dtc,v_x);
         v_dy = _mm512_fmadd_ps(v_vy,v_dtc,v_y);
         v_dz = _mm512_fmadd_ps(v_vz,v_dtc,v_z);
/* find particles going out of bounds */
/*       mm = 0; */
         v_mm = _mm512_setzero_epi32();
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                              */
/*       if (dx >= edgerx) {              */
/*          if (dx >= anx)                */
/*             ppart[j+npoff] = dx - anx; */
/*          mm = 2;                       */
/*       }                                */
         msk1 = _mm512_cmp_ps_mask(v_dx,v_edgerx,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dx;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_1,v_1);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dx,v_anx,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dx,v_anx);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dx = v_x;
            }
/*          if (dx < edgelx) {         */
/*             if (dx < 0.0) {         */
/*                dx += anx;           */
/*                if (dx < anx)        */
/*                   mm = 1;           */
/*                else                 */
/*                   dx = 0.0;         */
/*                ppart[j+npoff] = dx; */
/*             }                       */
/*             else {                  */
/*                mm = 1;              */
/*             }                       */
/*          }                          */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_1);
               msk2 = _mm512_cmp_ps_mask(v_dx,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dx,v_anx);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anx,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dx = v_x;
            }
         }
/*       if (dy >= edgery) {                    */
/*          if (dy >= any)                      */
/*             ppart[j+nppmx+npoff] = dy - any; */
/*          mm += 6;                            */
/*       }                                      */
         msk1 = _mm512_cmp_ps_mask(v_dy,v_edgery,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dy;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_3,v_3);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dy,v_any,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dy,v_any);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dy = v_x;
            }
/*          if (dy < edgely) {               */
/*             if (dy < 0.0) {               */
/*                dy += any;                 */
/*                if (dy < any)              */
/*                   mm += 3;                */
/*                else                       */
/*                   dy = 0.0;               */
/*                ppart[j+nppmx+npoff] = dy; */
/*             }                             */
/*             else {                        */
/*                mm += 3;                   */
/*             }                             */
/*          }                                */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_3);
               msk2 = _mm512_cmp_ps_mask(v_dy,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dy,v_any);
               msk1 = _mm512_cmp_ps_mask(v_x,v_any,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dy = v_x;
            }
         }
/*       if (dz >= edgerz) {                      */
/*          if (dz >= anz)                        */
/*             ppart[j+2*nppmx+npoff] = dz - anz; */
/*          mm += 18;                             */
/*       }                                        */
         msk1 = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dz;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_9,v_9);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dz,v_anz,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dz,v_anz);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dz = v_x;
            }
/*          if (dz < edgelz) {                 */
/*             if (dz < 0.0) {                 */
/*                dz += anz;                   */
/*                if (dz < anz)                */
/*                   mm += 9;                  */
/*                else                         */
/*                   dz = 0.0;                 */
/*                ppart[j+2*nppmx+npoff] = dz; */
/*             }                               */
/*             else {                          */
/*                mm += 9;                     */
/*             }                               */
/*          }                                  */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_9);
               msk2 = _mm512_cmp_ps_mask(v_dz,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dz,v_anz);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anz,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dz = v_x;
            }
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy;   */
/*       ppart[j+2*nppmx+npoff] = dz; */
         _mm512_store_ps(&ppart[j+npoff],v_dx);
         _mm512_store_ps(&ppart[j+nppmx+npoff],v_dy);
         _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_dz);
/* set new velocity */
/*       ppart[j+3*nppmx+npoff] = vx; */
/*       ppart[j+4*nppmx+npoff] = vy; */
/*       ppart[j+5*nppmx+npoff] = vz; */
         _mm512_store_ps(&ppart[j+3*nppmx+npoff],v_vx);
         _mm512_store_ps(&ppart[j+4*nppmx+npoff],v_vy);
         _mm512_store_ps(&ppart[j+5*nppmx+npoff],v_vz);
/* increment counters */
/*       if (mm > 0) {                                */
/*          ncl[mm+26*l-1] += 1;                      */
/*          ih += 1;                                  */
/*          if (ih <= ntmax) {                        */
/*             ihole[2*(ih+(ntmax+1)*l)] = j + i + 1; */
/*             ihole[1+2*(ih+(ntmax+1)*l)] = mm;      */
/*          }                                         */
/*          else {                                    */
/*             nh = 1;                                */
/*          }                                         */
/*       }                                            */
         _mm512_store_epi32(kk,v_mm);
         for (i = 0; i < 16; i++) {
            mm = kk[i];
            if (mm > 0) {
               ncl[mm+26*l-1] += 1;
               ih += 1; 
               if (ih <= ntmax) { 
                  ihole[2*(ih+(ntmax+1)*l)] = j + i + 1;
                  ihole[1+2*(ih+(ntmax+1)*l)] = mm;
               } 
               else {
                  nh = 1;
               }
            }
         }
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nm = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find electric field */
         nn = nm;
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+4];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];
         mm = nn + 4*mxv;
         dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);
         dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);
         dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);
         nn += 4*mxyv;
         acx = amx*sfxyz[nn] + amy*sfxyz[nn+4];
         acy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];
         acz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];
         mm = nn + 4*mxv;
         dx = dx + dzp*(acx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);
         dy = dy + dzp*(acy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);
         dz = dz + dzp*(acz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxyz[nn] + amy*sbxyz[nn+4];
         oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];
         oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];
         mm = nn + 4*mxv;
         ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);
         oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);
         oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);
         nn += 4*mxyv;
         acx = amx*sbxyz[nn] + amy*sbxyz[nn+4];
         acy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];
         acz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];
         mm = nn + 4*mxv;
         ox = ox + dzp*(acx + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);
         oy = oy + dzp*(acy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);
         oz = oz + dzp*(acz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+3*nppmx+npoff] + dx;
         acy = ppart[j+4*nppmx+npoff] + dy;
         acz = ppart[j+5*nppmx+npoff] + dz;
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
         dz = z + vz*dtc;
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
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
/* set new velocity */
         ppart[j+3*nppmx+npoff] = vx;
         ppart[j+4*nppmx+npoff] = vy;
         ppart[j+5*nppmx+npoff] = vz;
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
/*    sum2 += sum1; */
      _mm512_store_pd(&dd[0],v_sum1);
      for (j = 1; j < 8; j++) {
         dd[0] += dd[j];
      }
      sum2 += (sum1 + dd[0]);
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
void ckncgrbppush3lt(float ppart[], float fxyz[], float bxyz[],
                     int kpic[], float qbm, float dt, float dtc,
                     float ci, float *ek, int idimp, int nppmx, int nx,
                     int ny, int nz, int mx, int my, int mz, int nxv,
                     int nyv, int nzv, int mx1, int my1, int mxyz1,
                     int ipbc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   OpenMP/vector version using guard cells
   data read in tiles
   particles stored segmented array
   202 flops/particle, 4 divides, 2 sqrts, 54 loads, 6 stores
   input: all, output: ppart, ek
   momentum equations used are:
   px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fx(x(t),y(t),z(t))*dt)
   py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fy(x(t),y(t),z(t))*dt)
   pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
   omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
   omy = (q/m)*by(x(t),y(t),z(t))*gami,
   omz = (q/m)*bz(x(t),y(t),z(t))*gami,
   where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
   position equations used are:
   x(t+dt) = x(t) + px(t+dt/2)*dtg
   y(t+dt) = y(t) + py(t+dt/2)*dtg
   z(t+dt) = z(t) + pz(t+dt/2)*dtg
   where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
   pz(t+dt/2)*pz(t+dt/2))*ci*ci)
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
   bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = momentum px of particle n in tile m
   ppart[m][4][n] = momentum py of particle n in tile m
   ppart[m][5][n] = momentum pz of particle n in tile m
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
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
       (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
   requires KNC, ppart needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
   fxyz needs to have 4 components, although one is not used
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, nn, mm, ll, nm, mxv, myv, mxyv, nxyv;
   float qtmh, ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1;
   float acx, acy, acz, omxt, p2, gami, qtmg, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9, dtg;
   float x, y, z, vx, vy, vz;
   double sum1, sum2;
   __m512i v_noff, v_moff, v_loff, v_mxv4, v_mxyv4;
   __m512i v_nn, v_mm, v_ll, v_nm, v_it, v_perm;
   __m512 v_qtmh, v_ci2, v_dt, v_dtc, v_one, v_zero;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_gami, v_at, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s;
   __m512 v_two, v_half, v_ox, v_oy, v_oz;
   __m512d v_sum1, v_d;
   __mmask16 msk;
   __attribute__((aligned(64))) unsigned int kk[16];
   __attribute__((aligned(64))) double dd[8];
   __attribute__((aligned(64))) float sfxyz[4*MXV*MYV*MZV];
   __attribute__((aligned(64))) float sbxyz[4*MXV*MYV*MZV];
/* __attribute__((aligned(64))) float sfxyz[4*(mx+1)*(my+1)*(mz+1)]; */
/* __attribute__((aligned(64))) float sbxyz[4*(mx+1)*(my+1)*(mz+1)]; */
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   qtmh = 0.5f*qbm*dt;
   ci2 = ci*ci;
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
   v_mxv4 = _mm512_set1_epi32(4*mxv);
   v_mxyv4 = _mm512_set1_epi32(4*mxyv);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qtmh = _mm512_set1_ps(qtmh);
   v_ci2 = _mm512_set1_ps(ci2);
   v_dt = _mm512_set1_ps(dt);
   v_dtc = _mm512_set1_ps(dtc);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_two = _mm512_set1_ps(2.0f);
   v_half = _mm512_set1_ps(0.5f);
   v_edgelx = _mm512_set1_ps(edgelx);
   v_edgely = _mm512_set1_ps(edgely);
   v_edgelz = _mm512_set1_ps(edgelz);
   v_edgerx = _mm512_set1_ps(edgerx);
   v_edgery = _mm512_set1_ps(edgery);
   v_edgerz = _mm512_set1_ps(edgerz);
   v_sum1 = _mm512_set1_pd(0.0);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,nps,nn,mm,ll,nm,x,y,z,vx, \
vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt, \
omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2, \
gami,qtmg,dtg,sum1,v_noff,v_moff,v_loff,v_nn,v_mm,v_ll,v_nm,v_it,v_x, \
v_y,v_z,v_dxp,v_dyp,v_dzp,v_amx,v_amy,v_amz,v_dx1,v_dx,v_dy,v_dz,v_vx, \
v_vy,v_vz,v_ox,v_oy,v_oz,v_gami,v_at,v_d,v_sum1,a,b,c,d,e,f,g,h,p,q,r, \
s,msk,kk,dd,sfxyz,sbxyz) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* load local fields from global array */
      nn = (mx < nx-noff ? mx : nx-noff) + 1;
      mm = (my < ny-moff ? my : ny-moff) + 1;
      ll = (mz < nz-loff ? mz : nz-loff) + 1;
      nps = 4*(nn/4);
/* load electric field */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 0; i < nn; i++) {                          */
/*             sfxyz[4*(i+mxv*j+mxyv*k)]                        */
/*             = fxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];   */
/*             sfxyz[1+4*(i+mxv*j+mxyv*k)]                      */
/*             = fxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*             sfxyz[2+4*(i+mxv*j+mxyv*k)]                      */
/*             = fxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*          }                                                   */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&fxyz[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&fxyz[m+16]);
               m = 4*(i + mxv*j + mxyv*k);
               _mm512_packstorelo_ps(&sfxyz[m],v_at);
               _mm512_packstorehi_ps(&sfxyz[m+16],v_at);
            }
/* loop over remaining elements */
            for (i = nps; i < nn; i++) {
               sfxyz[4*(i+mxv*j+mxyv*k)]
               = fxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+4*(i+mxv*j+mxyv*k)]
               = fxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+4*(i+mxv*j+mxyv*k)]
               = fxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[3+4*(i+mxv*j+mxyv*k)]
               = fxyz[3+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
/* load magnetic field */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 0; i < nn; i++) {                          */
/*             sbxyz[4*(i+mxv*j+mxyv*k)]                        */
/*             = bxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];   */
/*             sbxyz[1+4*(i+mxv*j+mxyv*k)]                      */
/*             = bxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*             sbxyz[2+4*(i+mxv*j+mxyv*k)]                      */
/*             = bxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*          }                                                   */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&bxyz[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&bxyz[m+16]);
               m = 4*(i + mxv*j + mxyv*k);
               _mm512_packstorelo_ps(&sbxyz[m],v_at);
               _mm512_packstorehi_ps(&sbxyz[m+16],v_at);
            }
/* loop over remaining elements */
            for (i = nps; i < nn; i++) {
               sbxyz[4*(i+mxv*j+mxyv*k)]
               = bxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[1+4*(i+mxv*j+mxyv*k)]
               = bxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[2+4*(i+mxv*j+mxyv*k)]
               = bxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[3+4*(i+mxv*j+mxyv*k)]
               = bxyz[3+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
      nps = 16*(npp/16);
      sum1 = 0.0;
      v_sum1 = _mm512_set1_pd(0.0);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = x - (float) nn; */
/*       dyp = y - (float) mm; */
/*       dzp = z - (float) ll; */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_sub_ps(v_x,v_dxp);
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*       nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff)); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv4,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv4,v_mm));
         v_nm = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*       amx = 1.0f - dxp; */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_one,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/* find electric field */
/*       nn = nm; */
         _mm512_store_epi32(kk,v_nm);
/* load sfxyz[nn:nn+3] and sfxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[nn:nn+3] field components */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[nn+4:nn+7] field components */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of electric field */
/*       dx = amx*sfxyz[nn] + amy*sfxyz[nn+4];                     */
         v_dx = _mm512_mul_ps(v_amx,a);
         v_dx = _mm512_fmadd_ps(v_amy,p,v_dx);
/*       dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];                 */
         v_dy = _mm512_mul_ps(v_amx,b);
         v_dy = _mm512_fmadd_ps(v_amy,q,v_dy);
/*       dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];                 */
         v_dz = _mm512_mul_ps(v_amx,c);
         v_dz = _mm512_fmadd_ps(v_amy,r,v_dz);
/*       mm = nn + 4*mxv;                                     */
/* load sfxyz[mm:mm+3] and sfxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxv;                                     */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxv;                                       */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of electric field */
/*       dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);          */
         v_dx = _mm512_fmadd_ps(v_dyp,a,v_dx);
         v_dx = _mm512_fmadd_ps(v_dx1,p,v_dx);
         v_dx = _mm512_mul_ps(v_amz,v_dx);
/*       dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);      */
         v_dy = _mm512_fmadd_ps(v_dyp,b,v_dy);
         v_dy = _mm512_fmadd_ps(v_dx1,q,v_dy);
         v_dy = _mm512_mul_ps(v_amz,v_dy);
/*       dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);      */
         v_dz = _mm512_fmadd_ps(v_dyp,c,v_dz);
         v_dz = _mm512_fmadd_ps(v_dx1,r,v_dz);
         v_dz = _mm512_mul_ps(v_amz,v_dz);
/*       nn += 4*mxyv;                                           */
         v_nn = _mm512_add_epi32(v_nm,v_mxyv4);
         _mm512_store_epi32(kk,v_nn);
/* load sfxyz[nn:nn+3] and sfxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[nn:nn+3] field components */
/* where nn = nn + 4*mxyv;                                    */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*mxyv;                                      */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of electric field */
/*       vx = amx*sfxyz[nn] + amy*sfxyz[nn+4];                     */
         v_vx = _mm512_mul_ps(v_amx,a);
         v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*       vy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];                 */
         v_vy = _mm512_mul_ps(v_amx,b);
         v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*       vz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];                 */
         v_vz = _mm512_mul_ps(v_amx,c);
         v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*       mm = nn + 4*mxv;                                        */
/* load sfxyz[mm:mm+3] and sfxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
              &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                            */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                              */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of electric field */
/*       dx = dx + dzp*(vx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);     */
         v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
         v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
         v_dx = _mm512_fmadd_ps(v_dzp,v_vx,v_dx);
/*       dy = dy + dzp*(vy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]); */
         v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
         v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
         v_dy = _mm512_fmadd_ps(v_dzp,v_vy,v_dy);
/*       dz = dz + dzp*(vz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]); */
         v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
         v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
         v_dz = _mm512_fmadd_ps(v_dzp,v_vz,v_dz);
/* find magnetic field */
/*       nn = nm; */
         _mm512_store_epi32(kk,v_nm);
/* load sbxyz[nn:nn+3] and sbxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[nn:nn+3] field components */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[nn+4:nn+7] field components */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of magnetic field */
/*       ox = amx*sbxyz[nn] + amy*sbxyz[nn+4];                     */
         v_ox = _mm512_mul_ps(v_amx,a);
         v_ox = _mm512_fmadd_ps(v_amy,p,v_ox);
/*       oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];                 */
         v_oy = _mm512_mul_ps(v_amx,b);
         v_oy = _mm512_fmadd_ps(v_amy,q,v_oy);
/*       oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];                 */
         v_oz = _mm512_mul_ps(v_amx,c);
         v_oz = _mm512_fmadd_ps(v_amy,r,v_oz);
/*       mm = nn + 4*mxv;                                     */
/* load sbxyz[mm:mm+3] and sbxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxv;                                     */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxv;                                       */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of magnetic field */
/*       ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);          */
         v_ox = _mm512_fmadd_ps(v_dyp,a,v_ox);
         v_ox = _mm512_fmadd_ps(v_dx1,p,v_ox);
         v_ox = _mm512_mul_ps(v_amz,v_ox);
/*       oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);      */
         v_oy = _mm512_fmadd_ps(v_dyp,b,v_oy);
         v_oy = _mm512_fmadd_ps(v_dx1,q,v_oy);
         v_oy = _mm512_mul_ps(v_amz,v_oy);
/*       oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);      */
         v_oz = _mm512_fmadd_ps(v_dyp,c,v_oz);
         v_oz = _mm512_fmadd_ps(v_dx1,r,v_oz);
         v_oz = _mm512_mul_ps(v_amz,v_oz);
/*       nn += 4*mxyv;                                           */
         v_nn = _mm512_add_epi32(v_nm,v_mxyv4);
         _mm512_store_epi32(kk,v_nn);
/* load sbxyz[nn:nn+3] and sbxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[nn:nn+3] field components */
/* where nn = nn + 4*mxyv;                                    */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*mxyv;                                     */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of magnetic field */
/*       vx = amx*sbxyz[nn] + amy*sbxyz[nn+4];                     */
         v_vx = _mm512_mul_ps(v_amx,a);
         v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*       vy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];                 */
         v_vy = _mm512_mul_ps(v_amx,b);
         v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*       vz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];                 */
         v_vz = _mm512_mul_ps(v_amx,c);
         v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*       mm = nn + 4*mxv;                                        */
/* load sbxyz[mm:mm+3] and sbxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                            */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                           */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of magnetic field */
/*       ox = ox + dzp*(vx + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);     */
         v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
         v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
         v_ox = _mm512_fmadd_ps(v_dzp,v_vx,v_ox);
/*       oy = oy + dzp*(vy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]); */
         v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
         v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
         v_oy = _mm512_fmadd_ps(v_dzp,v_vy,v_oy);
/*       oz = oz + dzp*(vz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]); */
         v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
         v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
         v_oz = _mm512_fmadd_ps(v_dzp,v_vz,v_oz);
/* calculate half impulse */
/*       dx *= qtmh; */
/*       dy *= qtmh; */
/*       dz *= qtmh; */
         v_dx = _mm512_mul_ps(v_dx,v_qtmh);
         v_dy = _mm512_mul_ps(v_dy,v_qtmh);
         v_dz = _mm512_mul_ps(v_dz,v_qtmh);
/* half acceleration */
/*       acx = ppart[j+3*nppmx+npoff] + dx; */
/*       acy = ppart[j+4*nppmx+npoff] + dy; */
/*       acz = ppart[j+5*nppmx+npoff] + dz; */
         a = _mm512_add_ps(v_dx,_mm512_load_ps(&ppart[j+3*nppmx+npoff]));
         b = _mm512_add_ps(v_dy,_mm512_load_ps(&ppart[j+4*nppmx+npoff]));
         c = _mm512_add_ps(v_dz,_mm512_load_ps(&ppart[j+5*nppmx+npoff]));
/* find inverse gamma */
/*       p2 = acx*acx + acy*acy + acz*acz; */
         v_at = _mm512_fmadd_ps(b,b,_mm512_mul_ps(a,a));
         v_at = _mm512_fmadd_ps(c,c,v_at);
/*       gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_gami = _mm512_rsqrt23_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* full accuracy calculation */
         v_gami = _mm512_sqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one));
         v_gami = _mm512_div_ps(v_one,v_gami);
/* full accuracy calculation with SVML */
/*       v_gami = _mm512_invsqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* time-centered kinetic energy */
/*       sum1 += gami*p2/(1.0f + gami); */
         v_at = _mm512_mul_ps(v_gami,v_at);
         v_at = _mm512_div_ps(v_at,_mm512_add_ps(v_one,v_gami));
/* convert to double precision before accumulating */
         v_sum1 = _mm512_add_pd(v_sum1,_mm512_cvtpslo_pd(v_at));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_at,78));
         v_sum1 = _mm512_add_pd(v_sum1,v_d);
/* renormalize magnetic field */
/*       qtmg = qtmh*gami; */
         v_at = _mm512_mul_ps(v_qtmh,v_gami);
/* calculate cyclotron frequency */
/*       omxt = qtmg*ox; */
/*       omyt = qtmg*oy; */
/*       omzt = qtmg*oz; */
         e = _mm512_mul_ps(v_at,v_ox);
         f = _mm512_mul_ps(v_at,v_oy);
         g = _mm512_mul_ps(v_at,v_oz);
/* calculate rotation matrix */
/*       vx = omxt*omxt; */
         v_vx = _mm512_mul_ps(e,e);
/*       vy = omyt*omyt; */
         v_vy = _mm512_mul_ps(f,f);
/*       vz = omzt*omzt; */
         v_vz = _mm512_mul_ps(g,g);
/*       omt = omxt*omxt + omyt*omyt + omzt*omzt; */
         v_at = _mm512_add_ps(_mm512_add_ps(v_vx,v_vy),v_vz);
/*       anorm = 2.0f/(1.0f + omt); */
         d = _mm512_div_ps(v_two,_mm512_add_ps(v_one,v_at));
/*       omt = 0.5f*(1.0f - omt); */
         h = _mm512_mul_ps(v_half,_mm512_sub_ps(v_one,v_at));
/*       vx = (omt + vx)*acx; */
         v_vx = _mm512_mul_ps(_mm512_add_ps(h,v_vx),a);
/*       vy = (omt + vy)*acy; */
         v_vy = _mm512_mul_ps(_mm512_add_ps(h,v_vy),b);
/*       vz = (omt + vz)*acz; */
         v_vz = _mm512_mul_ps(_mm512_add_ps(h,v_vz),c);
/*       omt = omxt*omyt; */
         h = _mm512_mul_ps(e,f);
/*       vx = vx + (omzt + omt)*acy; */
         v_vx = _mm512_fmadd_ps(_mm512_add_ps(h,g),b,v_vx);
/*       vy = vy + (omt - omzt)*acx; */
         v_vy = _mm512_fmadd_ps(_mm512_sub_ps(h,g),a,v_vy);
/*       omt = omxt*omzt;  */
         h = _mm512_mul_ps(e,g);
/*       vx = vx + (omt - omyt)*acz; */
         v_vx = _mm512_fmadd_ps(_mm512_sub_ps(h,f),c,v_vx);
/*       vz = vz + (omt + omyt)*acx; */
         v_vz = _mm512_fmadd_ps(_mm512_add_ps(h,f),a,v_vz);
/*       omt = omyt*omzt; */
         h = _mm512_mul_ps(f,g);
/*       vy = vy + (omt + omxt)*acz; */
         v_vy = _mm512_fmadd_ps(_mm512_add_ps(h,e),c,v_vy);
/*       vz = vz + (omt - omxt)*acy; */
         v_vz = _mm512_fmadd_ps(_mm512_sub_ps(h,e),b,v_vz);
/* new momentum */
/*       vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm; */
/*       vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm; */
/*       vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm; */
         v_vx = _mm512_fmadd_ps(v_vx,d,v_dx);
         v_vy = _mm512_fmadd_ps(v_vy,d,v_dy);
         v_vz = _mm512_fmadd_ps(v_vz,d,v_dz);
/* update inverse gamma */
/*       p2 = vx*vx + vy*vy + vz*vz; */
         v_at = _mm512_fmadd_ps(v_vy,v_vy,_mm512_mul_ps(v_vx,v_vx));
         v_at = _mm512_fmadd_ps(v_vz,v_vz,v_at);
/*       dtg = dtc/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_at = _mm512_rsqrt23_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/*       v_at = _mm512_mul_ps(v_dtc,v_at);                            */
/* full accuracy calculation */
         v_at = _mm512_sqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one));
         v_at = _mm512_div_ps(v_dtc,v_at);
/* full accuracy calculation with SVML */
/*       v_gami = _mm512_invsqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/*       v_at = _mm512_div_ps(v_dtc,v_at);                              */
/* new position */
/*       dx = x + vx*dtg; */
/*       dy = y + vy*dtg; */
/*       dz = z + vz*dtg; */
         v_dx = _mm512_fmadd_ps(v_vx,v_at,v_x);
         v_dy = _mm512_fmadd_ps(v_vy,v_at,v_y);
         v_dz = _mm512_fmadd_ps(v_vz,v_at,v_z);
/* reflecting boundary conditions */
         if (ipbc==2) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             vx = -vx;                           */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             vy = -vy;                           */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/*          if ((dz < edgelz) || (dz >= edgerz)) { */
/*             dz = z;                             */
/*             vz = -vz;                           */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dz,v_edgerz,
                  _MM_CMPINT_GE));
            v_dz = _mm512_mask_blend_ps(msk,v_dz,v_z);
            v_vz = _mm512_mask_sub_ps(v_vz,msk,v_zero,v_vz);
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             vx = -vx;                           */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             vy = -vy;                           */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy;   */
/*       ppart[j+2*nppmx+npoff] = dz; */
         _mm512_store_ps(&ppart[j+npoff],v_dx);
         _mm512_store_ps(&ppart[j+nppmx+npoff],v_dy);
         _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_dz);
/* set new momentum */
/*       ppart[j+3*nppmx+npoff] = vx; */
/*       ppart[j+4*nppmx+npoff] = vy; */
/*       ppart[j+5*nppmx+npoff] = vz; */
         _mm512_store_ps(&ppart[j+3*nppmx+npoff],v_vx);
         _mm512_store_ps(&ppart[j+4*nppmx+npoff],v_vy);
         _mm512_store_ps(&ppart[j+5*nppmx+npoff],v_vz);
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nm = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find electric field */
         nn = nm;
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+4];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];
         mm = nn + 4*mxv;
         dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);
         dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);
         dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);
         nn += 4*mxyv;
         acx = amx*sfxyz[nn] + amy*sfxyz[nn+4];
         acy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];
         acz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];
         mm = nn + 4*mxv;
         dx = dx + dzp*(acx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);
         dy = dy + dzp*(acy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);
         dz = dz + dzp*(acz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxyz[nn] + amy*sbxyz[nn+4];
         oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];
         oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];
         mm = nn + 4*mxv;
         ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);
         oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);
         oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);
         nn += 4*mxyv;
         acx = amx*sbxyz[nn] + amy*sbxyz[nn+4];
         acy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];
         acz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];
         mm = nn + 4*mxv;
         ox = ox + dzp*(acx + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);
         oy = oy + dzp*(acy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);
         oz = oz + dzp*(acz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+3*nppmx+npoff] + dx;
         acy = ppart[j+4*nppmx+npoff] + dy;
         acz = ppart[j+5*nppmx+npoff] + dz;
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
         dz = z + vz*dtg;
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
            if ((dz < edgelz) || (dz >= edgerz)) {
               dz = z;
               vz = -vz;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               vx = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               vy = -vy;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
/* set new momentum */
         ppart[j+3*nppmx+npoff] = vx;
         ppart[j+4*nppmx+npoff] = vy;
         ppart[j+5*nppmx+npoff] = vz;
      }
/*    sum2 += sum1; */
      _mm512_store_pd(&dd[0],v_sum1);
      for (j = 1; j < 8; j++) {
         dd[0] += dd[j];
      }
      sum2 += (sum1 + dd[0]);
   }
/* normalize kinetic energy */
   *ek += sum2;
   return;
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void ckncgrbppushf3lt(float ppart[], float fxyz[], float bxyz[],
                      int kpic[], int ncl[], int ihole[], float qbm,
                      float dt, float dtc, float ci, float *ek,
                      int idimp, int nppmx, int nx, int ny, int nz,
                      int mx, int my, int mz, int nxv, int nyv, int nzv,
                      int mx1, int my1, int mxyz1, int ntmax,
                      int *irc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   also determines list of particles which are leaving this tile
   OpenMP/vector version using guard cells
   data read in tiles
   particles stored segmented array
   202 flops/particle, 4 divides, 2 sqrts, 54 loads, 6 stores
   input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
   momentum equations used are:
   px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fx(x(t),y(t),z(t))*dt)
   py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fy(x(t),y(t),z(t))*dt)
   pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fz(x(t),y(t),z(t))*dt)
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
   omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
   omy = (q/m)*by(x(t),y(t),z(t))*gami,
   omz = (q/m)*bz(x(t),y(t),z(t))*gami,
   where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
   position equations used are:
   x(t+dt) = x(t) + px(t+dt/2)*dtg
   y(t+dt) = y(t) + py(t+dt/2)*dtg
   z(t+dt) = z(t) + pz(t+dt/2)*dtg
   where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
   pz(t+dt/2)*pz(t+dt/2))*ci*ci)
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
   bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = momentum px of particle n in tile m
   ppart[m][4][n] = momentum py of particle n in tile m
   ppart[m][5][n] = momentum pz of particle n in tile m
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
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
       (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
   requires KNC, ppart needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
   fxyz needs to have 4 components, although one is not used
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, ii, ih, nh, nn, mm, ll, nm, mxv, myv, mxyv, nxyv;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1;
   float acx, acy, acz, omxt, p2, gami, qtmg, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9, dtg;
   float qtmh, ci2, x, y, z, vx, vy, vz;
   double sum1, sum2;
   __m512i v_noff, v_moff, v_loff, v_mxv4, v_mxyv4;
   __m512i v_nn, v_mm, v_ll, v_nm, v_it, v_0, v_1, v_3, v_9, v_perm;
   __m512 v_dt, v_dtc, v_one, v_zero, v_anx, v_any, v_anz;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_gami, v_at, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s;
   __m512 v_qtmh, v_ci2, v_two, v_half, v_ox, v_oy, v_oz;
   __m512d v_sum1, v_d;
   __mmask16 msk1, msk2;
   __attribute__((aligned(64))) unsigned int kk[16];
   __attribute__((aligned(64))) double dd[8];
   __attribute__((aligned(64))) float sfxyz[4*MXV*MYV*MZV];
   __attribute__((aligned(64))) float sbxyz[4*MXV*MYV*MZV];
/* __attribute__((aligned(64))) float sfxyz[4*(mx+1)*(my+1)*(mz+1)]; */
/* __attribute__((aligned(64))) float sbxyz[4*(mx+1)*(my+1)*(mz+1)]; */
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   qtmh = 0.5f*qbm*dt;
   ci2 = ci*ci;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
   sum2 = 0.0;
/* set boundary values */
   v_mxv4 = _mm512_set1_epi32(4*mxv);
   v_mxyv4 = _mm512_set1_epi32(4*mxyv);
   v_0 = _mm512_set1_epi32(0);
   v_1 = _mm512_set1_epi32(1);
   v_3 = _mm512_set1_epi32(3);
   v_9 = _mm512_set1_epi32(9);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qtmh = _mm512_set1_ps(qtmh);
   v_ci2 = _mm512_set1_ps(ci2);
   v_dt = _mm512_set1_ps(dt);
   v_dtc = _mm512_set1_ps(dtc);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_two = _mm512_set1_ps(2.0f);
   v_half = _mm512_set1_ps(0.5f);
   v_anx = _mm512_set1_ps(anx);
   v_any = _mm512_set1_ps(any);
   v_anz = _mm512_set1_ps(anz);
   v_sum1 = _mm512_set1_pd(0.0);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,ii,noff,moff,loff,npp,npoff,nps,nn,mm,ll,nm,ih,nh,x, \
y,z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,ox,oy,oz,acx,acy,acz, \
omxt,omyt,omzt,omt,anorm,rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9, \
edgelx,edgely,edgelz,edgerx,edgery,edgerz,p2,gami,qtmg,dtg,sum1,v_noff, \
v_moff,v_loff,v_nn,v_mm,v_ll,v_nm,v_it,v_x,v_y,v_z,v_dxp,v_dyp,v_dzp, \
v_amx,v_amy,v_amz,v_dx1,v_dx,v_dy,v_dz,v_vx,v_vy,v_vz,v_ox,v_oy,v_oz, \
v_gami,v_at,v_edgelx,v_edgely,v_edgelz,v_edgerx,v_edgery,v_edgerz,v_d, \
v_sum1,a,b,c,d,e,f,g,h,p,q,r,s,msk1,msk2,kk,dd,sfxyz,sbxyz) \
reduction(+:sum2)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
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
      v_edgelx = _mm512_set1_ps(edgelx);
      v_edgely = _mm512_set1_ps(edgely);
      v_edgelz = _mm512_set1_ps(edgelz);
      v_edgerx = _mm512_set1_ps(edgerx);
      v_edgery = _mm512_set1_ps(edgery);
      v_edgerz = _mm512_set1_ps(edgerz);
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
      ll += 1;
/* load local fields from global array */
      nps = 4*(nn/4);
/* load electric field */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 0; i < nn; i++) {                          */
/*             sfxyz[4*(i+mxv*j+mxyv*k)]                        */
/*             = fxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];   */
/*             sfxyz[1+4*(i+mxv*j+mxyv*k)]                      */
/*             = fxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*             sfxyz[2+4*(i+mxv*j+mxyv*k)]                      */
/*             = fxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*          }                                                   */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&fxyz[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&fxyz[m+16]);
               m = 4*(i + mxv*j + mxyv*k);
               _mm512_packstorelo_ps(&sfxyz[m],v_at);
               _mm512_packstorehi_ps(&sfxyz[m+16],v_at);
            }
/* loop over remaining elements */
            for (i = nps; i < nn; i++) {
               sfxyz[4*(i+mxv*j+mxyv*k)]
               = fxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[1+4*(i+mxv*j+mxyv*k)]
               = fxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[2+4*(i+mxv*j+mxyv*k)]
               = fxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sfxyz[3+4*(i+mxv*j+mxyv*k)]
               = fxyz[3+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
/* load magnetic field */
      for (k = 0; k < ll; k++) {
         for (j = 0; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 0; i < nn; i++) {                          */
/*             sbxyz[4*(i+mxv*j+mxyv*k)]                        */
/*             = bxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];   */
/*             sbxyz[1+4*(i+mxv*j+mxyv*k)]                      */
/*             = bxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*             sbxyz[2+4*(i+mxv*j+mxyv*k)]                      */
/*             = bxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]; */
/*          }                                                   */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&bxyz[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&bxyz[m+16]);
               m = 4*(i + mxv*j + mxyv*k);
               _mm512_packstorelo_ps(&sbxyz[m],v_at);
               _mm512_packstorehi_ps(&sbxyz[m+16],v_at);
            }
/* loop over remaining elements */
            for (i = nps; i < nn; i++) {
               sbxyz[4*(i+mxv*j+mxyv*k)]
               = bxyz[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[1+4*(i+mxv*j+mxyv*k)]
               = bxyz[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[2+4*(i+mxv*j+mxyv*k)]
               = bxyz[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
               sbxyz[3+4*(i+mxv*j+mxyv*k)]
               = bxyz[3+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))];
            }
         }
      }
/* clear counters */
/*    for (j = 0; j < 26; j++) { */
/*       ncl[j+26*l] = 0;        */
/*    }                         */
      memset((void*)&ncl[26*l],0,26*sizeof(int));
      nps = 16*(npp/16);
      sum1 = 0.0;
      v_sum1 = _mm512_set1_pd(0.0);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = x - (float) nn; */
/*       dyp = y - (float) mm; */
/*       dzp = z - (float) ll; */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_sub_ps(v_x,v_dxp);
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*       nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff)); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv4,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv4,v_mm));
         v_nm = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*       amx = 1.0f - dxp; */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_one,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/* find electric field */
/*       nn = nm; */
         _mm512_store_epi32(kk,v_nm);
/* load sfxyz[nn:nn+3] and sfxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[nn:nn+3] field components */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[nn+4:nn+7] field components */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of electric field */
/*       dx = amx*sfxyz[nn] + amy*sfxyz[nn+4];                     */
         v_dx = _mm512_mul_ps(v_amx,a);
         v_dx = _mm512_fmadd_ps(v_amy,p,v_dx);
/*       dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];                 */
         v_dy = _mm512_mul_ps(v_amx,b);
         v_dy = _mm512_fmadd_ps(v_amy,q,v_dy);
/*       dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];                 */
         v_dz = _mm512_mul_ps(v_amx,c);
         v_dz = _mm512_fmadd_ps(v_amy,r,v_dz);
/*       mm = nn + 4*mxv;                                     */
/* load sfxyz[mm:mm+3] and sfxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxv;                                     */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxv;                                       */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of electric field */
/*       dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);          */
         v_dx = _mm512_fmadd_ps(v_dyp,a,v_dx);
         v_dx = _mm512_fmadd_ps(v_dx1,p,v_dx);
         v_dx = _mm512_mul_ps(v_amz,v_dx);
/*       dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);      */
         v_dy = _mm512_fmadd_ps(v_dyp,b,v_dy);
         v_dy = _mm512_fmadd_ps(v_dx1,q,v_dy);
         v_dy = _mm512_mul_ps(v_amz,v_dy);
/*       dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);      */
         v_dz = _mm512_fmadd_ps(v_dyp,c,v_dz);
         v_dz = _mm512_fmadd_ps(v_dx1,r,v_dz);
         v_dz = _mm512_mul_ps(v_amz,v_dz);
/*       nn += 4*mxyv;                                           */
         v_nn = _mm512_add_epi32(v_nm,v_mxyv4);
         _mm512_store_epi32(kk,v_nn);
/* load sfxyz[nn:nn+3] and sfxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[nn:nn+3] field components */
/* where nn = nn + 4*mxyv;                                    */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*mxyv;                                      */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of electric field */
/*       vx = amx*sfxyz[nn] + amy*sfxyz[nn+4];                     */
         v_vx = _mm512_mul_ps(v_amx,a);
         v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*       vy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];                 */
         v_vy = _mm512_mul_ps(v_amx,b);
         v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*       vz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];                 */
         v_vz = _mm512_mul_ps(v_amx,c);
         v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*       mm = nn + 4*mxv;                                        */
/* load sfxyz[mm:mm+3] and sfxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
              &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sfxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sfxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sfxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                            */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sfxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                              */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of electric field */
/*       dx = dx + dzp*(vx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);     */
         v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
         v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
         v_dx = _mm512_fmadd_ps(v_dzp,v_vx,v_dx);
/*       dy = dy + dzp*(vy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]); */
         v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
         v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
         v_dy = _mm512_fmadd_ps(v_dzp,v_vy,v_dy);
/*       dz = dz + dzp*(vz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]); */
         v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
         v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
         v_dz = _mm512_fmadd_ps(v_dzp,v_vz,v_dz);
/* find magnetic field */
/*       nn = nm; */
         _mm512_store_epi32(kk,v_nm);
/* load sbxyz[nn:nn+3] and sbxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[nn:nn+3] field components */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[nn+4:nn+7] field components */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of magnetic field */
/*       ox = amx*sbxyz[nn] + amy*sbxyz[nn+4];                     */
         v_ox = _mm512_mul_ps(v_amx,a);
         v_ox = _mm512_fmadd_ps(v_amy,p,v_ox);
/*       oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];                 */
         v_oy = _mm512_mul_ps(v_amx,b);
         v_oy = _mm512_fmadd_ps(v_amy,q,v_oy);
/*       oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];                 */
         v_oz = _mm512_mul_ps(v_amx,c);
         v_oz = _mm512_fmadd_ps(v_amy,r,v_oz);
/*       mm = nn + 4*mxv;                                     */
/* load sbxyz[mm:mm+3] and sbxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxv;                                     */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxv;                                       */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of magnetic field */
/*       ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);          */
         v_ox = _mm512_fmadd_ps(v_dyp,a,v_ox);
         v_ox = _mm512_fmadd_ps(v_dx1,p,v_ox);
         v_ox = _mm512_mul_ps(v_amz,v_ox);
/*       oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);      */
         v_oy = _mm512_fmadd_ps(v_dyp,b,v_oy);
         v_oy = _mm512_fmadd_ps(v_dx1,q,v_oy);
         v_oy = _mm512_mul_ps(v_amz,v_oy);
/*       oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);      */
         v_oz = _mm512_fmadd_ps(v_dyp,c,v_oz);
         v_oz = _mm512_fmadd_ps(v_dx1,r,v_oz);
         v_oz = _mm512_mul_ps(v_amz,v_oz);
/*       nn += 4*mxyv;                                           */
         v_nn = _mm512_add_epi32(v_nm,v_mxyv4);
         _mm512_store_epi32(kk,v_nn);
/* load sbxyz[nn:nn+3] and sbxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
         mm = kk[0];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14];
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15];
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[nn:nn+3] field components */
/* where nn = nn + 4*mxyv;                                    */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*mxyv;                                     */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of magnetic field */
/*       vx = amx*sbxyz[nn] + amy*sbxyz[nn+4];                     */
         v_vx = _mm512_mul_ps(v_amx,a);
         v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*       vy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];                 */
         v_vy = _mm512_mul_ps(v_amx,b);
         v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*       vz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];                 */
         v_vz = _mm512_mul_ps(v_amx,c);
         v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*       mm = nn + 4*mxv;                                        */
/* load sbxyz[mm:mm+3] and sbxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
         mm = kk[0] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[1] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[2] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[3] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
         mm = kk[4] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[5] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[6] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[7] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
         mm = kk[8] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[9] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[10] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[11] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
         mm = kk[12] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[13] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
             &sbxyz[mm+16]);
         mm = kk[14] + 4*mxv;
         e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm]);
         e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         mm = kk[15] + 4*mxv;
         f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm]);
         f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
             &sbxyz[mm+16]);
         d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
         s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for sbxyz[mm:mm+3] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                            */
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
         f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
         g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
         b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for sbxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*mxyv + 4*mxv;                           */
         p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
         q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
         r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
         s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
         e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
         f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
         g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
         q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
         p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
         r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
         q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of magnetic field */
/*       ox = ox + dzp*(vx + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);     */
         v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
         v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
         v_ox = _mm512_fmadd_ps(v_dzp,v_vx,v_ox);
/*       oy = oy + dzp*(vy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]); */
         v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
         v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
         v_oy = _mm512_fmadd_ps(v_dzp,v_vy,v_oy);
/*       oz = oz + dzp*(vz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]); */
         v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
         v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
         v_oz = _mm512_fmadd_ps(v_dzp,v_vz,v_oz);
/* calculate half impulse */
/*       dx *= qtmh; */
/*       dy *= qtmh; */
/*       dz *= qtmh; */
         v_dx = _mm512_mul_ps(v_dx,v_qtmh);
         v_dy = _mm512_mul_ps(v_dy,v_qtmh);
         v_dz = _mm512_mul_ps(v_dz,v_qtmh);
/* half acceleration */
/*       acx = ppart[j+3*nppmx+npoff] + dx; */
/*       acy = ppart[j+4*nppmx+npoff] + dy; */
/*       acz = ppart[j+5*nppmx+npoff] + dz; */
         a = _mm512_add_ps(v_dx,_mm512_load_ps(&ppart[j+3*nppmx+npoff]));
         b = _mm512_add_ps(v_dy,_mm512_load_ps(&ppart[j+4*nppmx+npoff]));
         c = _mm512_add_ps(v_dz,_mm512_load_ps(&ppart[j+5*nppmx+npoff]));
/* find inverse gamma */
/*       p2 = acx*acx + acy*acy + acz*acz; */
         v_at = _mm512_fmadd_ps(b,b,_mm512_mul_ps(a,a));
         v_at = _mm512_fmadd_ps(c,c,v_at);
/*       gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_gami = _mm512_rsqrt23_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* full accuracy calculation */
         v_gami = _mm512_sqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one));
         v_gami = _mm512_div_ps(v_one,v_gami);
/* full accuracy calculation with SVML */
/*       v_gami = _mm512_invsqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* time-centered kinetic energy */
/*       sum1 += gami*p2/(1.0f + gami); */
         v_at = _mm512_mul_ps(v_gami,v_at);
         v_at = _mm512_div_ps(v_at,_mm512_add_ps(v_one,v_gami));
/* convert to double precision before accumulating */
         v_sum1 = _mm512_add_pd(v_sum1,_mm512_cvtpslo_pd(v_at));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_at,78));
         v_sum1 = _mm512_add_pd(v_sum1,v_d);
/* renormalize magnetic field */
/*       qtmg = qtmh*gami; */
         v_at = _mm512_mul_ps(v_qtmh,v_gami);
/* calculate cyclotron frequency */
/*       omxt = qtmg*ox; */
/*       omyt = qtmg*oy; */
/*       omzt = qtmg*oz; */
         e = _mm512_mul_ps(v_at,v_ox);
         f = _mm512_mul_ps(v_at,v_oy);
         g = _mm512_mul_ps(v_at,v_oz);
/* calculate rotation matrix */
/*       vx = omxt*omxt; */
         v_vx = _mm512_mul_ps(e,e);
/*       vy = omyt*omyt; */
         v_vy = _mm512_mul_ps(f,f);
/*       vz = omzt*omzt; */
         v_vz = _mm512_mul_ps(g,g);
/*       omt = omxt*omxt + omyt*omyt + omzt*omzt; */
         v_at = _mm512_add_ps(_mm512_add_ps(v_vx,v_vy),v_vz);
/*       anorm = 2.0f/(1.0f + omt); */
         d = _mm512_div_ps(v_two,_mm512_add_ps(v_one,v_at));
/*       omt = 0.5f*(1.0f - omt); */
         h = _mm512_mul_ps(v_half,_mm512_sub_ps(v_one,v_at));
/*       vx = (omt + vx)*acx; */
         v_vx = _mm512_mul_ps(_mm512_add_ps(h,v_vx),a);
/*       vy = (omt + vy)*acy; */
         v_vy = _mm512_mul_ps(_mm512_add_ps(h,v_vy),b);
/*       vz = (omt + vz)*acz; */
         v_vz = _mm512_mul_ps(_mm512_add_ps(h,v_vz),c);
/*       omt = omxt*omyt; */
         h = _mm512_mul_ps(e,f);
/*       vx = vx + (omzt + omt)*acy; */
         v_vx = _mm512_fmadd_ps(_mm512_add_ps(h,g),b,v_vx);
/*       vy = vy + (omt - omzt)*acx; */
         v_vy = _mm512_fmadd_ps(_mm512_sub_ps(h,g),a,v_vy);
/*       omt = omxt*omzt;  */
         h = _mm512_mul_ps(e,g);
/*       vx = vx + (omt - omyt)*acz; */
         v_vx = _mm512_fmadd_ps(_mm512_sub_ps(h,f),c,v_vx);
/*       vz = vz + (omt + omyt)*acx; */
         v_vz = _mm512_fmadd_ps(_mm512_add_ps(h,f),a,v_vz);
/*       omt = omyt*omzt; */
         h = _mm512_mul_ps(f,g);
/*       vy = vy + (omt + omxt)*acz; */
         v_vy = _mm512_fmadd_ps(_mm512_add_ps(h,e),c,v_vy);
/*       vz = vz + (omt - omxt)*acy; */
         v_vz = _mm512_fmadd_ps(_mm512_sub_ps(h,e),b,v_vz);
/* new momentum */
/*       vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm; */
/*       vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm; */
/*       vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm; */
         v_vx = _mm512_fmadd_ps(v_vx,d,v_dx);
         v_vy = _mm512_fmadd_ps(v_vy,d,v_dy);
         v_vz = _mm512_fmadd_ps(v_vz,d,v_dz);
/* update inverse gamma */
/*       p2 = vx*vx + vy*vy + vz*vz; */
         v_at = _mm512_fmadd_ps(v_vy,v_vy,_mm512_mul_ps(v_vx,v_vx));
         v_at = _mm512_fmadd_ps(v_vz,v_vz,v_at);
/*       dtg = dtc/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_at = _mm512_rsqrt23_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/*       v_at = _mm512_mul_ps(v_dtc,v_at);                            */
/* full accuracy calculation */
         v_at = _mm512_sqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one));
         v_at = _mm512_div_ps(v_dtc,v_at);
/* full accuracy calculation with SVML */
/*       v_gami = _mm512_invsqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/*       v_at = _mm512_div_ps(v_dtc,v_at);                              */
/* new position */
/*       dx = x + vx*dtg; */
/*       dy = y + vy*dtg; */
/*       dz = z + vz*dtg; */
         v_dx = _mm512_fmadd_ps(v_vx,v_at,v_x);
         v_dy = _mm512_fmadd_ps(v_vy,v_at,v_y);
         v_dz = _mm512_fmadd_ps(v_vz,v_at,v_z);
/* find particles going out of bounds */
/*       mm = 0; */
         v_mm = _mm512_setzero_epi32();
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                              */
/*       if (dx >= edgerx) {              */
/*          if (dx >= anx)                */
/*             ppart[j+npoff] = dx - anx; */
/*          mm = 2;                       */
/*       }                                */
         msk1 = _mm512_cmp_ps_mask(v_dx,v_edgerx,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dx;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_1,v_1);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dx,v_anx,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dx,v_anx);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dx = v_x;
            }
/*          if (dx < edgelx) {         */
/*             if (dx < 0.0) {         */
/*                dx += anx;           */
/*                if (dx < anx)        */
/*                   mm = 1;           */
/*                else                 */
/*                   dx = 0.0;         */
/*                ppart[j+npoff] = dx; */
/*             }                       */
/*             else {                  */
/*                mm = 1;              */
/*             }                       */
/*          }                          */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_1);
               msk2 = _mm512_cmp_ps_mask(v_dx,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dx,v_anx);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anx,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dx = v_x;
            }
         }
/*       if (dy >= edgery) {                    */
/*          if (dy >= any)                      */
/*             ppart[j+nppmx+npoff] = dy - any; */
/*          mm += 6;                            */
/*       }                                      */
         msk1 = _mm512_cmp_ps_mask(v_dy,v_edgery,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dy;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_3,v_3);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dy,v_any,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dy,v_any);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dy = v_x;
            }
/*          if (dy < edgely) {               */
/*             if (dy < 0.0) {               */
/*                dy += any;                 */
/*                if (dy < any)              */
/*                   mm += 3;                */
/*                else                       */
/*                   dy = 0.0;               */
/*                ppart[j+nppmx+npoff] = dy; */
/*             }                             */
/*             else {                        */
/*                mm += 3;                   */
/*             }                             */
/*          }                                */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_3);
               msk2 = _mm512_cmp_ps_mask(v_dy,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dy,v_any);
               msk1 = _mm512_cmp_ps_mask(v_x,v_any,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dy = v_x;
            }
         }
/*       if (dz >= edgerz) {                      */
/*          if (dz >= anz)                        */
/*             ppart[j+2*nppmx+npoff] = dz - anz; */
/*          mm += 18;                             */
/*       }                                        */
         msk1 = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dz;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_9,v_9);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dz,v_anz,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dz,v_anz);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dz = v_x;
            }
/*          if (dz < edgelz) {                 */
/*             if (dz < 0.0) {                 */
/*                dz += anz;                   */
/*                if (dz < anz)                */
/*                   mm += 9;                  */
/*                else                         */
/*                   dz = 0.0;                 */
/*                ppart[j+2*nppmx+npoff] = dz; */
/*             }                               */
/*             else {                          */
/*                mm += 9;                     */
/*             }                               */
/*          }                                  */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_9);
               msk2 = _mm512_cmp_ps_mask(v_dz,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dz,v_anz);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anz,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dz = v_x;
            }
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy;   */
/*       ppart[j+2*nppmx+npoff] = dz; */
         _mm512_store_ps(&ppart[j+npoff],v_dx);
         _mm512_store_ps(&ppart[j+nppmx+npoff],v_dy);
         _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_dz);
/* set new momentum */
/*       ppart[j+3*nppmx+npoff] = vx; */
/*       ppart[j+4*nppmx+npoff] = vy; */
/*       ppart[j+5*nppmx+npoff] = vz; */
         _mm512_store_ps(&ppart[j+3*nppmx+npoff],v_vx);
         _mm512_store_ps(&ppart[j+4*nppmx+npoff],v_vy);
         _mm512_store_ps(&ppart[j+5*nppmx+npoff],v_vz);
/* increment counters */
/*       if (mm > 0) {                                */
/*          ncl[mm+26*l-1] += 1;                      */
/*          ih += 1;                                  */
/*          if (ih <= ntmax) {                        */
/*             ihole[2*(ih+(ntmax+1)*l)] = j + i + 1; */
/*             ihole[1+2*(ih+(ntmax+1)*l)] = mm;      */
/*          }                                         */
/*          else {                                    */
/*             nh = 1;                                */
/*          }                                         */
/*       }                                            */
         _mm512_store_epi32(kk,v_mm);
         for (i = 0; i < 16; i++) {
            mm = kk[i];
            if (mm > 0) {
               ncl[mm+26*l-1] += 1;
               ih += 1; 
               if (ih <= ntmax) { 
                  ihole[2*(ih+(ntmax+1)*l)] = j + i + 1;
                  ihole[1+2*(ih+(ntmax+1)*l)] = mm;
               } 
               else {
                  nh = 1;
               }
            }
         }
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = x - (float) nn;
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nm = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
         amx = 1.0f - dxp;
         amy = 1.0f - dyp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amz = 1.0f - dzp;
         amy = dxp*amy;
/* find electric field */
         nn = nm;
         dx = amx*sfxyz[nn] + amy*sfxyz[nn+4];
         dy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];
         dz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];
         mm = nn + 4*mxv;
         dx = amz*(dx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);
         dy = amz*(dy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);
         dz = amz*(dz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);
         nn += 4*mxyv;
         acx = amx*sfxyz[nn] + amy*sfxyz[nn+4];
         acy = amx*sfxyz[nn+1] + amy*sfxyz[nn+1+4];
         acz = amx*sfxyz[nn+2] + amy*sfxyz[nn+2+4];
         mm = nn + 4*mxv;
         dx = dx + dzp*(acx + dyp*sfxyz[mm] + dx1*sfxyz[mm+4]);
         dy = dy + dzp*(acy + dyp*sfxyz[mm+1] + dx1*sfxyz[mm+1+4]);
         dz = dz + dzp*(acz + dyp*sfxyz[mm+2] + dx1*sfxyz[mm+2+4]);
/* find magnetic field */
         nn = nm;
         ox = amx*sbxyz[nn] + amy*sbxyz[nn+4];
         oy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];
         oz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];
         mm = nn + 4*mxv;
         ox = amz*(ox + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);
         oy = amz*(oy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);
         oz = amz*(oz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);
         nn += 4*mxyv;
         acx = amx*sbxyz[nn] + amy*sbxyz[nn+4];
         acy = amx*sbxyz[nn+1] + amy*sbxyz[nn+1+4];
         acz = amx*sbxyz[nn+2] + amy*sbxyz[nn+2+4];
         mm = nn + 4*mxv;
         ox = ox + dzp*(acx + dyp*sbxyz[mm] + dx1*sbxyz[mm+4]);
         oy = oy + dzp*(acy + dyp*sbxyz[mm+1] + dx1*sbxyz[mm+1+4]);
         oz = oz + dzp*(acz + dyp*sbxyz[mm+2] + dx1*sbxyz[mm+2+4]);
/* calculate half impulse */
         dx *= qtmh;
         dy *= qtmh;
         dz *= qtmh;
/* half acceleration */
         acx = ppart[j+3*nppmx+npoff] + dx;
         acy = ppart[j+4*nppmx+npoff] + dy;
         acz = ppart[j+5*nppmx+npoff] + dz;
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
         dz = z + vz*dtg;
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
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
/* set new momentum */
         ppart[j+3*nppmx+npoff] = vx;
         ppart[j+4*nppmx+npoff] = vy;
         ppart[j+5*nppmx+npoff] = vz;
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
/*    sum2 += sum1; */
      _mm512_store_pd(&dd[0],v_sum1);
      for (j = 1; j < 8; j++) {
         dd[0] += dd[j];
      }
      sum2 += (sum1 + dd[0]);
/* set error and end of file flag */
      if (nh > 0) {
         *irc = ih;
         ih = -ih;
      }
      ihole[2*(ntmax+1)*l] = ih;
   }
/* normalize kinetic energy */
   *ek += sum2;
   return;
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void ckncgppost3lt(float ppart[], float q[], int kpic[], float qm,
                   int nppmx, int idimp, int mx, int my, int mz,
                   int nxv, int nyv, int nzv, int mx1, int my1,
                   int mxyz1) {
/* for 3d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   OpenMP/vector version using guard cells
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
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
   requires KNC, ppart needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, nn, mm, ll, nm, lm, mxv, myv, mxyv, nxyv;
   float x, y, z, w, dx1, dxp, dyp, dzp, amx, amy, amz;
   __m512i v_noff, v_moff, v_loff, v_mxv, v_mxyv;
   __m512i v_nn, v_mm, v_ll, v_it;
   __m512 v_qm, v_one;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_as, v_at;
   __m512 a, b, c, d, e, f, g, h, qp, qr;
   __mmask16 msk, msks, v_m;
   __attribute__((aligned(64))) unsigned int kk[16];
   __attribute__((aligned(64))) float sq[MXV*MYV*MZV];
/* __attribute__((aligned(64))) float sq[(mx+1)*(my+1)*(mz+1)]; */
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx + 1;
   myv = my + 1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   v_mxv = _mm512_set1_epi32(mxv);
   v_mxyv = _mm512_set1_epi32(mxyv);
   v_qm = _mm512_set1_ps(qm);
   v_one = _mm512_set1_ps(1.0f);
   v_at = _mm512_set_ps(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        1.);
   v_m = _mm512_cmp_ps_mask(v_at,v_one,_MM_CMPINT_LT);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,nps,nn,mm,ll,nm,lm,x,y,z,w, \
dxp,dyp,dzp,amx,amy,amz,dx1,v_noff,v_moff,v_loff,v_nn,v_mm,v_ll,v_it, \
v_x,v_y,v_z,v_dxp,v_dyp,v_dzp,v_amx,v_amy,v_amz,v_dx1,v_at,v_as,a,b,c, \
d,e,f,g,h,qp,qr,msk,msks,kk,sq)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* zero out local accumulator */
/*    for (j = 0; j < mxyv*(mz+1); j++) { */
/*       sq[j] = 0.0f;                    */
/*    }                                   */
      memset((void*)sq,0,mxyv*(mz+1)*sizeof(float));
      nps = 16*(npp/16);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = qm*(x - (float) nn); */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
/*       dyp = y - (float) mm; */
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
/*       dzp = z - (float) ll; */
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*       nn = nn - noff + mxv*(mm - moff) + mxyv*(ll - loff); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv,v_mm));
         v_nn = _mm512_add_epi32(v_nn,v_it);
/*       amx = qm - dxp;   */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_qm,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*       a = amx*amz; */
/*       b = amy*amz; */
/*       d = dyp*amz; */
/*       d = dx1*amz; */
         a = _mm512_mul_ps(v_amx,v_amz);
         b = _mm512_mul_ps(v_amy,v_amz);
         c = _mm512_mul_ps(v_dyp,v_amz);
         d = _mm512_mul_ps(v_dx1,v_amz);
/*       e = amx*dzp; */
/*       f = amy*dzp; */
/*       g = dyp*dzp; */
/*       h = dx1*dzp; */
         e = _mm512_mul_ps(v_amx,v_dzp);
         f = _mm512_mul_ps(v_amy,v_dzp);
         g = _mm512_mul_ps(v_dyp,v_dzp);
         h = _mm512_mul_ps(v_dx1,v_dzp);
         _mm512_store_epi32(kk,v_nn);
/* deposit charge */
/*       x = sq[nn] + amx*amz;       */
/*       y = sq[nn+1] + amy*amz;     */
/*       z = sq[nn+mxv] + dyp*amz;   */
/*       w = sq[nn+1+mxv] + dx1*amz; */
/*       sq[nn] = x;                 */
/*       sq[nn+1] = y;               */
/*       sq[nn+mxv] = z;             */
/*       sq[nn+1+mxv] = w;           */
/*       mm = nn + mxyv;             */
/*       x = sq[mm] + amx*dzp;       */
/*       y = sq[mm+1] + amy*dzp;     */
/*       z = sq[mm+mxv] + dyp*dzp;   */
/*       w = sq[mm+1+mxv] + dx1*dzp; */
/*       sq[mm] = x;                 */
/*       sq[mm+1] = y;               */
/*       sq[mm+mxv] = z;             */
/*       sq[mm+1+mxv] = w;           */
/* deposit charge for two particles at a time */
         for (i = 0; i < 8; i++) {
/* first particle */
            mm = kk[2*i];
            msk = _mm512_int2mask(3<<(2*i));
            msks = _mm512_int2mask(2<<(2*i));
            qp = _mm512_mask_loadunpacklo_ps(qp,msk,&sq[mm]);
            qp = _mm512_mask_loadunpackhi_ps(qp,msk,&sq[mm+16]);
            v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)a,msks,
                   (__m512i)b,177);
            qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
            _mm512_mask_packstorelo_ps(&sq[mm],msk,qp);
            _mm512_mask_packstorehi_ps(&sq[mm+16],msk,qp);
            ll = mm + mxv;
            qr = _mm512_mask_loadunpacklo_ps(qr,msk,&sq[ll]);
            qr = _mm512_mask_loadunpackhi_ps(qr,msk,&sq[ll+16]);
            v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)c,msks,
                   (__m512i)d,177);
            qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
            _mm512_mask_packstorelo_ps(&sq[ll],msk,qr);
            _mm512_mask_packstorehi_ps(&sq[ll+16],msk,qr);
            mm = mm + mxyv;
            qp = _mm512_mask_loadunpacklo_ps(qp,msk,&sq[mm]);
            qp = _mm512_mask_loadunpackhi_ps(qp,msk,&sq[mm+16]);
            v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)e,msks,
                   (__m512i)f,177);
            qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
            _mm512_mask_packstorelo_ps(&sq[mm],msk,qp);
            _mm512_mask_packstorehi_ps(&sq[mm+16],msk,qp);
            ll = mm + mxv;
            qr = _mm512_mask_loadunpacklo_ps(qr,msk,&sq[ll]);
            qr = _mm512_mask_loadunpackhi_ps(qr,msk,&sq[ll+16]);
            v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)g,msks,
                   (__m512i)h,177);
            qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
            _mm512_mask_packstorelo_ps(&sq[ll],msk,qr);
            _mm512_mask_packstorehi_ps(&sq[ll+16],msk,qr);
/* second particle */
            mm = kk[2*i+1];
            msks = _mm512_int2mask(1<<(2*i));
            qp = _mm512_mask_loadunpacklo_ps(qp,msk,&sq[mm]);
            qp = _mm512_mask_loadunpackhi_ps(qp,msk,&sq[mm+16]);
            v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)b,msks,
                   (__m512i)a,177);
            qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
            _mm512_mask_packstorelo_ps(&sq[mm],msk,qp);
            _mm512_mask_packstorehi_ps(&sq[mm+16],msk,qp);
            ll = mm + mxv;
            qr = _mm512_mask_loadunpacklo_ps(qr,msk,&sq[ll]);
            qr = _mm512_mask_loadunpackhi_ps(qr,msk,&sq[ll+16]);
            v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)d,msks,
                  (__m512i)c,177);
            qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
            _mm512_mask_packstorelo_ps(&sq[ll],msk,qr);
            _mm512_mask_packstorehi_ps(&sq[ll+16],msk,qr);
            mm = mm + mxyv;
            qp = _mm512_mask_loadunpacklo_ps(qp,msk,&sq[mm]);
            qp = _mm512_mask_loadunpackhi_ps(qp,msk,&sq[mm+16]);
            v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)f,msks,
                  (__m512i)e,177);
            qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
            _mm512_mask_packstorelo_ps(&sq[mm],msk,qp);
            _mm512_mask_packstorehi_ps(&sq[mm+16],msk,qp);
            ll = mm + mxv;
            qr = _mm512_mask_loadunpacklo_ps(qr,msk,&sq[ll]);
            qr = _mm512_mask_loadunpackhi_ps(qr,msk,&sq[ll+16]);
            v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)h,msks,
                   (__m512i)g,177);
            qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
            _mm512_mask_packstorelo_ps(&sq[ll],msk,qr);
            _mm512_mask_packstorehi_ps(&sq[ll+16],msk,qr);
         }
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = nn - noff + mxv*(mm - moff) + mxyv*(ll - loff);
         amx = qm - dxp;
         amy = 1.0f - dyp;
         amz = 1.0f - dzp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amy = dxp*amy;
/* deposit charge */
         x = sq[nn] + amx*amz;
         y = sq[nn+1] + amy*amz;
         z = sq[nn+mxv] + dyp*amz;
         w = sq[nn+1+mxv] + dx1*amz;
         sq[nn] = x;
         sq[nn+1] = y;
         sq[nn+mxv] = z;
         sq[nn+1+mxv] = w;
         mm = nn + mxyv;
         x = sq[mm] + amx*dzp;
         y = sq[mm+1] + amy*dzp;
         z = sq[mm+mxv] + dyp*dzp;
         w = sq[mm+1+mxv] + dx1*dzp;
         sq[mm] = x;
         sq[mm+1] = y;
         sq[mm+mxv] = z;
         sq[mm+1+mxv] = w;
      }
/* deposit charge to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      nps = 16*(nn/16);
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 1; i < nn; i++) {              */
/*             q[i+noff+nxv*(j+moff)+nxyv*(k+loff)] */
/*             += sq[i+mxv*j+mxyv*k];               */
/*          }                                       */
            for (i = 0; i < nps; i+=16) {
               m = i + mxv*j + mxyv*k;
               v_as = _mm512_loadunpacklo_ps(v_as,&sq[m]);
               v_as = _mm512_loadunpackhi_ps(v_as,&sq[m+16]);
               m = i + noff + nxv*(j + moff) + nxyv*(k + loff);
               v_at = _mm512_loadunpacklo_ps(v_at,&q[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&q[m+16]);
/* skip add for first element for i = 0 */
               if (i==0)
                  v_at = _mm512_mask_add_ps(v_at,v_m,v_at,v_as);
               else
                  v_at = _mm512_add_ps(v_at,v_as);
               _mm512_packstorelo_ps(&q[m],v_at);
               _mm512_packstorehi_ps(&q[m+16],v_at);
            }
/* loop over remaining elements */
            m = 1 > nps ? 1 : nps;
            for (i = m ; i < nn; i++) {
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
void cknc2gppost3lt(float ppart[], float q[], int kpic[], float qm,
                    int nppmx, int idimp, int mx, int my, int mz,
                    int nxv, int nyv, int nzv, int mx1, int my1,
                    int mxyz1) {
/* for 3d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   OpenMP/vector version using guard cells
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
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
   requires KNC, ppart needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, nn, mm, ll, nm, lm, mxv, myv, mxyv, nxyv;
   float x, y, z, w, dx1, dxp, dyp, dzp, amx, amy, amz;
   __m512i v_noff, v_moff, v_loff, v_mxv, v_mxyv;
   __m512i v_nn, v_mm, v_ll, v_it;
   __m512 v_qm, v_one;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_as, v_at;
   __mmask16 v_m;
   __attribute__((aligned(64))) unsigned int kk[16];
   typedef union vfloat {float v[16]; __m512 v16;} vf;
   __attribute__((aligned(64))) float sq[MXV*MYV*MZV];
/* __attribute__((aligned(64))) float sq[(mx+1)*(my+1)*(mz+1)]; */
   vf vv[8];
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx + 1;
   myv = my + 1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   v_mxv = _mm512_set1_epi32(mxv);
   v_mxyv = _mm512_set1_epi32(mxyv);
   v_qm = _mm512_set1_ps(qm);
   v_one = _mm512_set1_ps(1.0f);
   v_at = _mm512_set_ps(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        1.);
   v_m = _mm512_cmp_ps_mask(v_at,v_one,_MM_CMPINT_LT);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,nps,nn,mm,ll,nm,lm,x,y,z,w, \
dxp,dyp,dzp,amx,amy,amz,dx1,v_noff,v_moff,v_loff,v_nn,v_mm,v_ll,v_it, \
v_x,v_y,v_z,v_dxp,v_dyp,v_dzp,v_amx,v_amy,v_amz,v_dx1,v_at,v_as,kk,sq,vv)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* zero out local accumulator */
/*    for (j = 0; j < mxyv*(mz+1); j++) { */
/*       sq[j] = 0.0f;                    */
/*    }                                   */
      memset((void*)sq,0,mxyv*(mz+1)*sizeof(float));
      nps = 16*(npp/16);
/* vector loop over particles in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = qm*(x - (float) nn); */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
/*       dyp = y - (float) mm; */
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
/*       dzp = z - (float) ll; */
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*       nn = nn - noff + mxv*(mm - moff) + mxyv*(ll - loff); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv,v_mm));
         v_nn = _mm512_add_epi32(v_nn,v_it);
/*       amx = qm - dxp;   */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_qm,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*       x = amx*amz; */
/*       y = amy*amz; */
/*       z = dyp*amz; */
/*       w = dx1*amz; */
         vv[0].v16 = _mm512_mul_ps(v_amx,v_amz);
         vv[1].v16 = _mm512_mul_ps(v_amy,v_amz);
         vv[2].v16 = _mm512_mul_ps(v_dyp,v_amz);
         vv[3].v16 = _mm512_mul_ps(v_dx1,v_amz);
         vv[4].v16 = _mm512_mul_ps(v_amx,v_dzp);
         vv[5].v16 = _mm512_mul_ps(v_amy,v_dzp);
         vv[6].v16 = _mm512_mul_ps(v_dyp,v_dzp);
         vv[7].v16 = _mm512_mul_ps(v_dx1,v_dzp);
         _mm512_store_epi32(kk,v_nn);
/* deposit charge */
/*       x = sq[nn] + amx*amz;       */
/*       y = sq[nn+1] + amy*amz;     */
/*       z = sq[nn+mxv] + dyp*amz;   */
/*       w = sq[nn+1+mxv] + dx1*amz; */
/*       sq[nn] = x;                 */
/*       sq[nn+1] = y;               */
/*       sq[nn+mxv] = z;             */
/*       sq[nn+1+mxv] = w;           */
/*       mm = nn + mxyv;            */
/*       x = sq[mm] + amx*dzp;       */
/*       y = sq[mm+1] + amy*dzp;     */
/*       z = sq[mm+mxv] + dyp*dzp;   */
/*       w = sq[mm+1+mxv] + dx1*dzp; */
/*       sq[mm] = x;                 */
/*       sq[mm+1] = y;               */
/*       sq[mm+mxv] = z;             */
/*       sq[mm+1+mxv] = w;           */
         for (i = 0; i < 16; i++) {
            nn = kk[i];
            x = sq[nn] + vv[0].v[i];
            y = sq[nn+1] + vv[1].v[i];
            z = sq[nn+mxv] + vv[2].v[i];
            w = sq[nn+1+mxv] + vv[3].v[i];
            sq[nn] = x;
            sq[nn+1] = y;
            sq[nn+mxv] = z;
            sq[nn+1+mxv] = w;
            mm = nn + mxyv;
            x = sq[mm] + vv[4].v[i];
            y = sq[mm+1] + vv[5].v[i];
            z = sq[mm+mxv] + vv[6].v[i];
            w = sq[mm+1+mxv] + vv[7].v[i];
            sq[mm] = x;
            sq[mm+1] = y;
            sq[mm+mxv] = z;
            sq[mm+1+mxv] = w;
         }
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = nn - noff + mxv*(mm - moff) + mxyv*(ll - loff);
         amx = qm - dxp;
         amy = 1.0f - dyp;
         amz = 1.0f - dzp;
         dx1 = dxp*dyp;
         dyp = amx*dyp;
         amx = amx*amy;
         amy = dxp*amy;
/* deposit charge */
         x = sq[nn] + amx*amz;
         y = sq[nn+1] + amy*amz;
         z = sq[nn+mxv] + dyp*amz;
         w = sq[nn+1+mxv] + dx1*amz;
         sq[nn] = x;
         sq[nn+1] = y;
         sq[nn+mxv] = z;
         sq[nn+1+mxv] = w;
         mm = nn + mxyv;
         x = sq[mm] + amx*dzp;
         y = sq[mm+1] + amy*dzp;
         z = sq[mm+mxv] + dyp*dzp;
         w = sq[mm+1+mxv] + dx1*dzp;
         sq[mm] = x;
         sq[mm+1] = y;
         sq[mm+mxv] = z;
         sq[mm+1+mxv] = w;
     }
/* deposit charge to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      nps = 16*(nn/16);
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 1; i < nn; i++) {              */
/*             q[i+noff+nxv*(j+moff)+nxyv*(k+loff)] */
/*             += sq[i+mxv*j+mxyv*k];               */
/*          }                                       */
            for (i = 0; i < nps; i+=16) {
               m = i + mxv*j + mxyv*k;
               v_as = _mm512_loadunpacklo_ps(v_as,&sq[m]);
               v_as = _mm512_loadunpackhi_ps(v_as,&sq[m+16]);
               m = i + noff + nxv*(j + moff) + nxyv*(k + loff);
               v_at = _mm512_loadunpacklo_ps(v_at,&q[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&q[m+16]);
/* skip add for first element for i = 0 */
               if (i==0)
                  v_at = _mm512_mask_add_ps(v_at,v_m,v_at,v_as);
               else
                  v_at = _mm512_add_ps(v_at,v_as);
               _mm512_packstorelo_ps(&q[m],v_at);
               _mm512_packstorehi_ps(&q[m+16],v_at);
            }
/* loop over remaining elements */
            m = 1 > nps ? 1 : nps;
            for (i = m ; i < nn; i++) {
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
void ckncgjppost3lt(float ppart[], float cu[], int kpic[], float qm,
                    float dt, int nppmx, int idimp, int nx, int ny,
                    int nz, int mx, int my, int mz, int nxv, int nyv,
                    int nzv, int mx1, int my1, int mxyz1, int ipbc) {
/* for 3d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   OpenMP/vector version using guard cells
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
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
   requires KNC, part needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
   cu needs to have 4 components, although one is not used
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, nn, mm, ll, ii, nm, lm, mxv, myv, mxyv, nxyv;
   float edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz;
   float x, y, z;
   __m512i v_noff, v_moff, v_loff, v_mxv4, v_mxyv4;
   __m512i v_nn, v_mm, v_ll, v_it, v_perm;
   __m512 v_qm, v_dt, v_one, v_zero;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_at, v_as, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s, t, u, v, ws, wt, wu, wv;
   __m512 cp, cr;
   __mmask16 msk, v_m;
   __attribute__((aligned(64))) unsigned int kk[16];
   __attribute__((aligned(64))) float scu[4*MXV*MYV*MZV];
/* __attribute__((aligned(64))) float scu[4*(mx+1)*(my+1)*(mz+1)]; */
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
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
   v_mxv4 = _mm512_set1_epi32(4*mxv);
   v_mxyv4 = _mm512_set1_epi32(4*mxyv);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qm = _mm512_set1_ps(qm);
   v_dt = _mm512_set1_ps(dt);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_edgelx = _mm512_set1_ps(edgelx);
   v_edgely = _mm512_set1_ps(edgely);
   v_edgelz = _mm512_set1_ps(edgelz);
   v_edgerx = _mm512_set1_ps(edgerx);
   v_edgery = _mm512_set1_ps(edgery);
   v_edgerz = _mm512_set1_ps(edgerz);
   v_at = _mm512_set_ps(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,
                        1.);
   v_m = _mm512_cmp_ps_mask(v_at,v_one,_MM_CMPINT_LT);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,nps,nn,mm,ll,ii,nm,lm,x,y, \
z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,v_noff,v_moff,v_loff, \
v_nn,v_mm,v_ll,v_it,v_x,v_y,v_z,v_dxp,v_dyp,v_dzp,v_amx,v_amy,v_amz, \
v_dx1,v_dx,v_dy,v_dz,v_vx,v_vy,v_vz,v_at,v_as,a,b,c,d,e,f,g,h,p,q,r, \
s,t,u,v,ws,wt,wu,wv,cp,cr,msk,kk,scu)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* zero out local accumulator */
/*    for (j = 0; j < 4*mxyv*(mz+1); j++) { */
/*       scu[j] = 0.0f;                     */
/*    }                                     */
      memset((void*)scu,0,4*mxyv*(mz+1)*sizeof(float));
      nps = 16*(npp/16);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = qm*(x - (float) nn); */
/*       dyp = y - (float) mm;      */
/*       dzp = z - (float) ll;      */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*       nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff)); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv4,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv4,v_mm));
         v_nn = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*       amx = qm - dxp;   */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_qm,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*       a = amx*amz; */
/*       b = amy*amz; */
/*       c = dyp*amz; */
/*       d = dx1*amz; */
         a = _mm512_mul_ps(v_amx,v_amz);
         b = _mm512_mul_ps(v_amy,v_amz);
         c = _mm512_mul_ps(v_dyp,v_amz);
         d = _mm512_mul_ps(v_dx1,v_amz);
/*       e = amx*dzp; */
/*       f = amy*dzp; */
/*       g = dyp*dzp; */
/*       h = dx1*dzp; */
         e = _mm512_mul_ps(v_amx,v_dzp);
         f = _mm512_mul_ps(v_amy,v_dzp);
         g = _mm512_mul_ps(v_dyp,v_dzp);
         h = _mm512_mul_ps(v_dx1,v_dzp);
/* deposit current */
/*       vx = ppart[j+3*nppmx+npoff]; */
/*       vy = ppart[j+4*nppmx+npoff]; */
/*       vz = ppart[j+5*nppmx+npoff]; */
         v_vx = _mm512_load_ps(&ppart[j+3*nppmx+npoff]);
         v_vy = _mm512_load_ps(&ppart[j+4*nppmx+npoff]);
         v_vz = _mm512_load_ps(&ppart[j+5*nppmx+npoff]);
         v_ll = _mm512_add_epi32(v_nn,v_mxyv4);
/* deposit charge for one particle at a time */
         for (i = 0; i < 16; i++) {
            ii = i >> 2;
            if (i==(ii<<2)) {
               switch (ii)
               {
                  case 0:
/* replicate velocities of first group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,0);
                     q = _mm512_permute4f128_ps(v_vy,0);
                     r = _mm512_permute4f128_ps(v_vz,0);
/* regroup weights for first group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(a,
                         _mm512_int2mask(61680),b,177);
                     t = _mm512_mask_permute4f128_ps(c,
                         _mm512_int2mask(61680),d,177);
                     u = _mm512_mask_permute4f128_ps(e,
                         _mm512_int2mask(61680),f,177);
                     v = _mm512_mask_permute4f128_ps(g,
                         _mm512_int2mask(61680),h,177);
                     break;
                  case 1:
/* replicate velocities of second group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,85);
                     q = _mm512_permute4f128_ps(v_vy,85);
                     r = _mm512_permute4f128_ps(v_vz,85);
/* regroup weights for second group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(b,
                         _mm512_int2mask(3855),a,177);
                     t = _mm512_mask_permute4f128_ps(d,
                         _mm512_int2mask(3855),c,177);
                     u = _mm512_mask_permute4f128_ps(f,
                         _mm512_int2mask(3855),e,177);
                     v = _mm512_mask_permute4f128_ps(h,
                         _mm512_int2mask(3855),g,177);
                     break;
                  case 2:
/* replicate velocities of third group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,170);
                     q = _mm512_permute4f128_ps(v_vy,170);
                     r = _mm512_permute4f128_ps(v_vz,170);
/* regroup weights for third group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(a,
                         _mm512_int2mask(61680),b,177);
                     s = _mm512_permute4f128_ps(s,78);
                     t = _mm512_mask_permute4f128_ps(c,
                         _mm512_int2mask(61680),d,177);
                     t = _mm512_permute4f128_ps(t,78);
                     u = _mm512_mask_permute4f128_ps(e,
                         _mm512_int2mask(61680),f,177);
                     u = _mm512_permute4f128_ps(u,78);
                     v = _mm512_mask_permute4f128_ps(g,
                         _mm512_int2mask(61680),h,177);
                     v = _mm512_permute4f128_ps(v,78);
                     break;
                  case 3:
/* replicate velocities of fourth group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,255);
                     q = _mm512_permute4f128_ps(v_vy,255);
                     r = _mm512_permute4f128_ps(v_vz,255);
/* regroup weights for fourth group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(b,
                         _mm512_int2mask(3855),a,177);
                     s = _mm512_permute4f128_ps(s,78);
                     t = _mm512_mask_permute4f128_ps(d,
                         _mm512_int2mask(3855),c,177);
                     t = _mm512_permute4f128_ps(t,78);
                     u = _mm512_mask_permute4f128_ps(f,
                         _mm512_int2mask(3855),e,177);
                     u = _mm512_permute4f128_ps(u,78);
                     v = _mm512_mask_permute4f128_ps(h,
                         _mm512_int2mask(3855),g,177);
                     v = _mm512_permute4f128_ps(v,78);
                     break;
               }
            }
            v_it = _mm512_setzero_epi32();
            switch (i-(ii<<2))
            {
/* first particle */
               case 0:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                         _mm512_int2mask(170),(__m512i)q,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,78);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,0);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,0);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,0);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,0);
                  break;
/* second particle */
               case 1:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                         _mm512_int2mask(85),(__m512i)p,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,24);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,85);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,85);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,85);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,85);
                  break;
/* third particle */
               case 2:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                         _mm512_int2mask(170),(__m512i)q,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)r,
                         _mm512_int2mask(51),(__m512i)v_at,78);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                        _mm512_int2mask(119),(__m512i)s,170);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,170);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,170);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,170);
                   break;
/* fourth particle */
               case 3:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                         _mm512_int2mask(85),(__m512i)p,177);
                  v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_at,78);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,177);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,255);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,255);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,255);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,255);
                  break;
            }
            _mm512_store_epi32(kk,v_nn);
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dx = amx*amz;        */
/*          dy = amy*amz;        */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          dx = dyp*amz;        */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            mm = kk[i];
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),ws,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cp);
/*          mm = nn + 4*mxv;                                 */
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dx = dyp*amz;        */
/*          dy = dx1*amz;        */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            mm = kk[i] + 4*mxv;
            cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
                 &scu[mm]);
            cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
                 &scu[mm+16]);
            cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wt,cr);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cr);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cr);
            _mm512_store_epi32(kk,v_ll);
/*          nn += 4*mxyv;                                    */
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dx = amx*dzp;        */
/*          dy = amy*dzp;        */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            mm = kk[i];
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wu,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cp);
/*          mm = nn + 4*mxv;                                 */
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dx = dyp*dzp;        */
/*          dy = dx1*dzp;        */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            mm = kk[i] + 4*mxv;
            cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
                 &scu[mm]);
            cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
                 &scu[mm+16]);
            cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wv,cr);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cr);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cr);
         }
/* advance position half a time-step */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
/*       dz = z + vz*dt; */
         v_dx = _mm512_fmadd_ps(v_vx,v_dt,v_x);
         v_dy = _mm512_fmadd_ps(v_vy,v_dt,v_y);
         v_dz = _mm512_fmadd_ps(v_vz,v_dt,v_z);
/* reflecting boundary conditions */
        if (ipbc==2) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+3*nppmx+npoff] = -vx;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+3*nppmx+npoff],v_vx);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             ppart[j+4*nppmx+npoff] = -vy;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+4*nppmx+npoff],v_vy);
/*          if ((dz < edgelz) || (dz >= edgerz)) { */
/*             dz = z;                             */
/*             ppart[j+5*nppmx+npoff] = -vz;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dz,v_edgerz,
                  _MM_CMPINT_GE));
            v_dz = _mm512_mask_blend_ps(msk,v_dz,v_z);
            v_vz = _mm512_mask_sub_ps(v_vz,msk,v_zero,v_vz);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+5*nppmx+npoff],v_vz);
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+3*nppmx+npoff] = -vx;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/* write output if test result is true for any particle */
           if (msk)
               _mm512_store_ps(&ppart[j+3*nppmx+npoff],v_vx);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             ppart[j+4*nppmx+npoff] = -vy;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+4*nppmx+npoff],v_vy);
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy;   */
/*       ppart[j+2*nppmx+npoff] = dz; */
         _mm512_store_ps(&ppart[j+npoff],v_dx);
         _mm512_store_ps(&ppart[j+nppmx+npoff],v_dy);
         _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_dz);
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
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
         vx = ppart[j+3*nppmx+npoff];
         vy = ppart[j+4*nppmx+npoff];
         vz = ppart[j+5*nppmx+npoff];
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*amz;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*amz;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         dx = amx*dzp;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
         dy = amy*dzp;
         nn += 4*mxyv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*dzp;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*dzp;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+3*nppmx+npoff] = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+4*nppmx+npoff] = -vy;
            }
            if ((dz < edgelz) || (dz >= edgerz)) {
               dz = z;
               ppart[j+5*nppmx+npoff] = -vz;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+3*nppmx+npoff] = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+4*nppmx+npoff] = -vy;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      nps = 4*(nn/4);
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 1; i < nn; i++) {                     */
/*             cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]   */
/*             += scu[4*(i+mxv*j+mxyv*k)];                 */
/*             cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[1+4*(i+mxv*j+mxyv*k)];               */
/*             cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[2+4*(i+mxv*j+mxyv*k)];               */
/*          }                                              */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + mxv*j + mxyv*k);
               v_as = _mm512_loadunpacklo_ps(v_as,&scu[m]);
               v_as = _mm512_loadunpackhi_ps(v_as,&scu[m+16]);
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&cu[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&cu[m+16]);
/* skip add for first elements for i = 0 */
               if (i==0)
                  v_at = _mm512_mask_add_ps(v_at,v_m,v_at,v_as);
               else
                  v_at = _mm512_add_ps(v_at,v_as);
               _mm512_packstorelo_ps(&cu[m],v_at);
               _mm512_packstorehi_ps(&cu[m+16],v_at);
            }
/* loop over remaining elements */
            m = 1 > nps ? 1 : nps;
            for (i = m; i < nn; i++) {
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(i+mxv*j+mxyv*k)];
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*j+mxyv*k)];
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*j+mxyv*k)];
            }
         }
      }
/* deposit current to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*(j+moff)+nxyv*loff)] += scu[4*(i+mxv*j)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[1+4*(i+mxv*j)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[2+4*(i+mxv*j)];
            if (lm > mz) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*j+mxyv*(lm-1))];
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
            cu[4*(i+noff+nxv*moff+nxyv*(k+loff))] += scu[4*(i+mxyv*k)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[1+4*(i+mxyv*k)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[2+4*(i+mxyv*k)];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*k)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[1+4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[2+4*(mxv*j+mxyv*k)];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(nm-1+mxv*j+mxyv*k)];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[1+4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[2+4*(i+mxyv*(lm-1))];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[1+4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[2+4*(mxv*j+mxyv*(lm-1))];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(nm-1+mxv*j+mxyv*(lm-1))];
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
void ckncgjppostf3lt(float ppart[], float cu[], int kpic[], int ncl[],
                     int ihole[], float qm, float dt, int nppmx,
                     int idimp, int nx, int ny, int nz, int mx, int my,
                     int mz, int nxv, int nyv, int nzv, int mx1,
                     int my1, int mxyz1, int ntmax, int *irc) {
/* for 3d code, this subroutine calculates particle current density
   using first-order linear interpolation, with periodic boundary
   conditions.
   in addition, particle positions are advanced a half time-step
   also determines list of particles which are leaving this tile
   OpenMP/vector version using guard cells
   data deposited in tiles
   particles stored segmented array
   69 flops/particle, 30 loads, 27 stores
   input: all except ncl, ihole, irc
   output: ppart, cu, ncl, ihole, ek, irc
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
   cu[l][k][j][i] = ith component of current density at grid point j,k,l
   kpic[l] = number of particles in tile l
   ncl[l][i] = number of particles going to destination i, tile l
   ihole[l][:][0] = location of hole in array left by departing particle
   ihole[l][:][1] = direction destination of particle leaving hole
   all for tile l
   ihole[l][0][0] = ih, number of holes left (error, if negative)
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
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   requires KNC, part needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
   cu needs to have 4 components, although one is not used
   optimized version
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, ih, nh, nn, mm, ll, ii, nm, lm, mxv, myv, mxyv;
   int nxyv;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz;
   float x, y, z;
   __m512i v_noff, v_moff, v_loff, v_mxv4, v_mxyv4;
   __m512i v_nn, v_mm, v_ll, v_it, v_0, v_1, v_3, v_9, v_perm;
   __m512 v_qm, v_dt, v_one, v_zero, v_anx, v_any, v_anz;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_at, v_as, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s, t, u, v, ws, wt, wu, wv;
   __m512 cp, cr;
   __mmask16 msk, msk1, msk2, v_m;
   __attribute__((aligned(64))) unsigned int kk[16];
   __attribute__((aligned(64))) float scu[4*MXV*MYV*MZV];
/* __attribute__((aligned(64))) float scu[4*(mx+1)*(my+1)*(mz+1)]; */
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
/* set boundary values */
   v_mxv4 = _mm512_set1_epi32(4*mxv);
   v_mxyv4 = _mm512_set1_epi32(4*mxyv);
   v_0 = _mm512_set1_epi32(0);
   v_1 = _mm512_set1_epi32(1);
   v_3 = _mm512_set1_epi32(3);
   v_9 = _mm512_set1_epi32(9);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qm = _mm512_set1_ps(qm);
   v_dt = _mm512_set1_ps(dt);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_anx = _mm512_set1_ps(anx);
   v_any = _mm512_set1_ps(any);
   v_anz = _mm512_set1_ps(anz);
   v_at = _mm512_set_ps(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,
                        1.);
   v_m = _mm512_cmp_ps_mask(v_at,v_one,_MM_CMPINT_LT);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,nps,nn,mm,ll,ih,nh,ii,nm,lm, \
x,y,z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,edgelx,edgely, \
edgelz,edgerx,edgery,edgerz,v_noff,v_moff,v_loff,v_nn,v_mm,v_ll,v_it, \
v_x,v_y,v_z,v_dxp,v_dyp,v_dzp,v_amx,v_amy,v_amz,v_dx1,v_dx,v_dy,v_dz, \
v_vx,v_vy,v_vz,v_edgelx,v_edgely,v_edgelz,v_edgerx,v_edgery,v_edgerz, \
v_at,v_as,a,b,c,d,e,f,g,h,p,q,r,s,t,u,v,ws,wt,wu,wv,cp,cr,msk,msk1, \
msk2,kk,scu)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
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
      v_edgelx = _mm512_set1_ps(edgelx);
      v_edgely = _mm512_set1_ps(edgely);
      v_edgelz = _mm512_set1_ps(edgelz);
      v_edgerx = _mm512_set1_ps(edgerx);
      v_edgery = _mm512_set1_ps(edgery);
      v_edgerz = _mm512_set1_ps(edgerz);
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
      ll += 1;
/* zero out local accumulator */
/*    for (j = 0; j < 4*mxyv*(mz+1); j++) { */
/*       scu[j] = 0.0f;                     */
/*    }                                     */
      memset((void*)scu,0,4*mxyv*(mz+1)*sizeof(float));
/* clear counters */
/*    for (j = 0; j < 26; j++) { */
/*       ncl[j+26*l] = 0;        */
/*    }                         */
      memset((void*)&ncl[26*l],0,26*sizeof(int));
      nps = 16*(npp/16);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = qm*(x - (float) nn); */
/*       dyp = y - (float) mm;      */
/*       dzp = z - (float) ll;      */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*       nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff)); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv4,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv4,v_mm));
         v_nn = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*       amx = qm - dxp;   */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_qm,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*       a = amx*amz; */
/*       b = amy*amz; */
/*       c = dyp*amz; */
/*       d = dx1*amz; */
         a = _mm512_mul_ps(v_amx,v_amz);
         b = _mm512_mul_ps(v_amy,v_amz);
         c = _mm512_mul_ps(v_dyp,v_amz);
         d = _mm512_mul_ps(v_dx1,v_amz);
/*       e = amx*dzp; */
/*       f = amy*dzp; */
/*       g = dyp*dzp; */
/*       h = dx1*dzp; */
         e = _mm512_mul_ps(v_amx,v_dzp);
         f = _mm512_mul_ps(v_amy,v_dzp);
         g = _mm512_mul_ps(v_dyp,v_dzp);
         h = _mm512_mul_ps(v_dx1,v_dzp);
/* deposit current */
/*       vx = ppart[j+3*nppmx+npoff]; */
/*       vy = ppart[j+4*nppmx+npoff]; */
/*       vz = ppart[j+5*nppmx+npoff]; */
         v_vx = _mm512_load_ps(&ppart[j+3*nppmx+npoff]);
         v_vy = _mm512_load_ps(&ppart[j+4*nppmx+npoff]);
         v_vz = _mm512_load_ps(&ppart[j+5*nppmx+npoff]);
         v_ll = _mm512_add_epi32(v_nn,v_mxyv4);
/* deposit charge for one particle at a time */
         for (i = 0; i < 16; i++) {
            ii = i >> 2;
            if (i==(ii<<2)) {
               switch (ii)
               {
                  case 0:
/* replicate velocities of first group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,0);
                     q = _mm512_permute4f128_ps(v_vy,0);
                     r = _mm512_permute4f128_ps(v_vz,0);
/* regroup weights for first group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(a,
                         _mm512_int2mask(61680),b,177);
                     t = _mm512_mask_permute4f128_ps(c,
                         _mm512_int2mask(61680),d,177);
                     u = _mm512_mask_permute4f128_ps(e,
                         _mm512_int2mask(61680),f,177);
                     v = _mm512_mask_permute4f128_ps(g,
                         _mm512_int2mask(61680),h,177);
                     break;
                  case 1:
/* replicate velocities of second group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,85);
                     q = _mm512_permute4f128_ps(v_vy,85);
                     r = _mm512_permute4f128_ps(v_vz,85);
/* regroup weights for second group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(b,
                         _mm512_int2mask(3855),a,177);
                     t = _mm512_mask_permute4f128_ps(d,
                         _mm512_int2mask(3855),c,177);
                     u = _mm512_mask_permute4f128_ps(f,
                         _mm512_int2mask(3855),e,177);
                     v = _mm512_mask_permute4f128_ps(h,
                         _mm512_int2mask(3855),g,177);
                     break;
                  case 2:
/* replicate velocities of third group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,170);
                     q = _mm512_permute4f128_ps(v_vy,170);
                     r = _mm512_permute4f128_ps(v_vz,170);
/* regroup weights for third group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(a,
                         _mm512_int2mask(61680),b,177);
                     s = _mm512_permute4f128_ps(s,78);
                     t = _mm512_mask_permute4f128_ps(c,
                         _mm512_int2mask(61680),d,177);
                     t = _mm512_permute4f128_ps(t,78);
                     u = _mm512_mask_permute4f128_ps(e,
                         _mm512_int2mask(61680),f,177);
                     u = _mm512_permute4f128_ps(u,78);
                     v = _mm512_mask_permute4f128_ps(g,
                         _mm512_int2mask(61680),h,177);
                     v = _mm512_permute4f128_ps(v,78);
                     break;
                  case 3:
/* replicate velocities of fourth group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,255);
                     q = _mm512_permute4f128_ps(v_vy,255);
                     r = _mm512_permute4f128_ps(v_vz,255);
/* regroup weights for fourth group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(b,
                         _mm512_int2mask(3855),a,177);
                     s = _mm512_permute4f128_ps(s,78);
                     t = _mm512_mask_permute4f128_ps(d,
                         _mm512_int2mask(3855),c,177);
                     t = _mm512_permute4f128_ps(t,78);
                     u = _mm512_mask_permute4f128_ps(f,
                         _mm512_int2mask(3855),e,177);
                     u = _mm512_permute4f128_ps(u,78);
                     v = _mm512_mask_permute4f128_ps(h,
                         _mm512_int2mask(3855),g,177);
                     v = _mm512_permute4f128_ps(v,78);
                     break;
               }
            }
            v_it = _mm512_setzero_epi32();
            switch (i-(ii<<2))
            {
/* first particle */
               case 0:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                         _mm512_int2mask(170),(__m512i)q,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,78);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,0);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,0);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,0);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,0);
                  break;
/* second particle */
               case 1:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                         _mm512_int2mask(85),(__m512i)p,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,24);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,85);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,85);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,85);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,85);
                  break;
/* third particle */
               case 2:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                         _mm512_int2mask(170),(__m512i)q,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)r,
                         _mm512_int2mask(51),(__m512i)v_at,78);
/* reorder weights */
                 ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                        _mm512_int2mask(119),(__m512i)s,170);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,170);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,170);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,170);
                   break;
/* fourth particle */
               case 3:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                         _mm512_int2mask(85),(__m512i)p,177);
                  v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_at,78);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,177);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,255);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,255);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,255);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,255);
                  break;
            }
            _mm512_store_epi32(kk,v_nn);
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dx = amx*amz;        */
/*          dy = amy*amz;        */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          dx = dyp*amz;        */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            mm = kk[i];
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),ws,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cp);
/*          mm = nn + 4*mxv;                                 */
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dx = dyp*amz;        */
/*          dy = dx1*amz;        */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            mm = kk[i] + 4*mxv;
            cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
                 &scu[mm]);
            cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
                 &scu[mm+16]);
            cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wt,cr);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cr);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cr);
            _mm512_store_epi32(kk,v_ll);
/*          nn += 4*mxyv;                                    */
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dx = amx*dzp;        */
/*          dy = amy*dzp;        */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            mm = kk[i];
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wu,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cp);
/*          mm = nn + 4*mxv;                                 */
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dx = dyp*dzp;        */
/*          dy = dx1*dzp;        */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            mm = kk[i] + 4*mxv;
            cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
                 &scu[mm]);
            cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
                 &scu[mm+16]);
            cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wv,cr);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cr);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cr);
         }
/* advance position half a time-step */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
/*       dz = z + vz*dt; */
         v_dx = _mm512_fmadd_ps(v_vx,v_dt,v_x);
         v_dy = _mm512_fmadd_ps(v_vy,v_dt,v_y);
         v_dz = _mm512_fmadd_ps(v_vz,v_dt,v_z);
/* find particles going out of bounds */
/*       mm = 0; */
         v_mm = _mm512_setzero_epi32();
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                              */
/*       if (dx >= edgerx) {              */
/*          if (dx >= anx)                */
/*             ppart[j+npoff] = dx - anx; */
/*          mm = 2;                       */
/*       }                                */
         msk1 = _mm512_cmp_ps_mask(v_dx,v_edgerx,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dx;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_1,v_1);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dx,v_anx,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dx,v_anx);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dx = v_x;
            }
/*          if (dx < edgelx) {         */
/*             if (dx < 0.0) {         */
/*                dx += anx;           */
/*                if (dx < anx)        */
/*                   mm = 1;           */
/*                else                 */
/*                   dx = 0.0;         */
/*                ppart[j+npoff] = dx; */
/*             }                       */
/*             else {                  */
/*                mm = 1;              */
/*             }                       */
/*          }                          */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_1);
               msk2 = _mm512_cmp_ps_mask(v_dx,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dx,v_anx);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anx,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dx = v_x;
            }
         }
/*       if (dy >= edgery) {                    */
/*          if (dy >= any)                      */
/*             ppart[j+nppmx+npoff] = dy - any; */
/*          mm += 6;                            */
/*       }                                      */
         msk1 = _mm512_cmp_ps_mask(v_dy,v_edgery,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dy;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_3,v_3);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dy,v_any,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dy,v_any);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dy = v_x;
            }
/*          if (dy < edgely) {               */
/*             if (dy < 0.0) {               */
/*                dy += any;                 */
/*                if (dy < any)              */
/*                   mm += 3;                */
/*                else                       */
/*                   dy = 0.0;               */
/*                ppart[j+nppmx+npoff] = dy; */
/*             }                             */
/*             else {                        */
/*                mm += 3;                   */
/*             }                             */
/*          }                                */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_3);
               msk2 = _mm512_cmp_ps_mask(v_dy,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dy,v_any);
               msk1 = _mm512_cmp_ps_mask(v_x,v_any,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dy = v_x;
            }
         }
/*       if (dz >= edgerz) {                      */
/*          if (dz >= anz)                        */
/*             ppart[j+2*nppmx+npoff] = dz - anz; */
/*          mm += 18;                             */
/*       }                                        */
         msk1 = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dz;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_9,v_9);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dz,v_anz,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dz,v_anz);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dz = v_x;
            }
/*          if (dz < edgelz) {                 */
/*             if (dz < 0.0) {                 */
/*                dz += anz;                   */
/*                if (dz < anz)                */
/*                   mm += 9;                  */
/*                else                         */
/*                   dz = 0.0;                 */
/*                ppart[j+2*nppmx+npoff] = dz; */
/*             }                               */
/*             else {                          */
/*                mm += 9;                     */
/*             }                               */
/*          }                                  */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_9);
               msk2 = _mm512_cmp_ps_mask(v_dz,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dz,v_anz);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anz,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dz = v_x;
            }
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy;   */
/*       ppart[j+2*nppmx+npoff] = dz; */
         _mm512_store_ps(&ppart[j+npoff],v_dx);
         _mm512_store_ps(&ppart[j+nppmx+npoff],v_dy);
         _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_dz);
/* increment counters */
/*       if (mm > 0) {                                */
/*          ncl[mm+26*l-1] += 1;                      */
/*          ih += 1;                                  */
/*          if (ih <= ntmax) {                        */
/*             ihole[2*(ih+(ntmax+1)*l)] = j + i + 1; */
/*             ihole[1+2*(ih+(ntmax+1)*l)] = mm;      */
/*          }                                         */
/*          else {                                    */
/*             nh = 1;                                */
/*          }                                         */
/*       }                                            */
         _mm512_store_epi32(kk,v_mm);
         for (i = 0; i < 16; i++) {
            mm = kk[i];
            if (mm > 0) {
               ncl[mm+26*l-1] += 1;
               ih += 1; 
               if (ih <= ntmax) { 
                  ihole[2*(ih+(ntmax+1)*l)] = j + i + 1;
                  ihole[1+2*(ih+(ntmax+1)*l)] = mm;
               } 
               else {
                  nh = 1;
               }
            }
         }
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
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
         vx = ppart[j+3*nppmx+npoff];
         vy = ppart[j+4*nppmx+npoff];
         vz = ppart[j+5*nppmx+npoff];
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*amz;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*amz;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         dx = amx*dzp;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
         dy = amy*dzp;
         nn += 4*mxyv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*dzp;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*dzp;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
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
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
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
/* deposit current to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      nps = 4*(nn/4);
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 1; i < nn; i++) {                     */
/*             cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]   */
/*             += scu[4*(i+mxv*j+mxyv*k)];                 */
/*             cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[1+4*(i+mxv*j+mxyv*k)];               */
/*             cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[2+4*(i+mxv*j+mxyv*k)];               */
/*          }                                              */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + mxv*j + mxyv*k);
               v_as = _mm512_loadunpacklo_ps(v_as,&scu[m]);
               v_as = _mm512_loadunpackhi_ps(v_as,&scu[m+16]);
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&cu[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&cu[m+16]);
/* skip add for first elements for i = 0 */
               if (i==0)
                  v_at = _mm512_mask_add_ps(v_at,v_m,v_at,v_as);
               else
                  v_at = _mm512_add_ps(v_at,v_as);
               _mm512_packstorelo_ps(&cu[m],v_at);
               _mm512_packstorehi_ps(&cu[m+16],v_at);
            }
/* loop over remaining elements */
            m = 1 > nps ? 1 : nps;
            for (i = m; i < nn; i++) {
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(i+mxv*j+mxyv*k)];
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*j+mxyv*k)];
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*j+mxyv*k)];
            }
         }
      }
/* deposit current to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*(j+moff)+nxyv*loff)] += scu[4*(i+mxv*j)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[1+4*(i+mxv*j)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[2+4*(i+mxv*j)];
            if (lm > mz) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*j+mxyv*(lm-1))];
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
            cu[4*(i+noff+nxv*moff+nxyv*(k+loff))] += scu[4*(i+mxyv*k)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[1+4*(i+mxyv*k)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[2+4*(i+mxyv*k)];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*k)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[1+4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[2+4*(mxv*j+mxyv*k)];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(nm-1+mxv*j+mxyv*k)];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[1+4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[2+4*(i+mxyv*(lm-1))];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[1+4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[2+4*(mxv*j+mxyv*(lm-1))];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(nm-1+mxv*j+mxyv*(lm-1))];
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
   return;
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void ckncgrjppost3lt(float ppart[], float cu[], int kpic[], float qm,
                     float dt, float ci, int nppmx, int idimp, int nx,
                     int ny, int nz, int mx, int my, int mz, int nxv,
                     int nyv, int nzv, int mx1, int my1, int mxyz1,
                     int ipbc) {
/* for 3d code, this subroutine calculates particle current density
   using first-order linear interpolation for relativistic particles
   in addition, particle positions are advanced a half time-step
   OpenMP/vector version using guard cells
   data deposited in tiles
   particles stored segmented array
   79 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
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
   and qci = qm*pi*gami, where i = x,y,z
   where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = momentum vx of particle n in tile m
   ppart[m][4][n] = momentum vy of particle n in tile m
   ppart[m][5][n] = momentum vz of particle n in tile m
   cu[l][k][j][i] = ith component of current density at grid point j,k,l
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
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
   requires KNC, part needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
   cu needs to have 4 components, although one is not used
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, nn, mm, ll, ii, nm, lm, mxv, myv, mxyv, nxyv;
   float ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz;
   float p2, gami;
   float x, y, z, ux, uy, uz;
   __m512i v_noff, v_moff, v_loff, v_mxv4, v_mxyv4;
   __m512i v_nn, v_mm, v_ll, v_it, v_perm;
   __m512 v_qm, v_ci2, v_dt, v_one, v_zero;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_gami, v_at, v_as, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_ux, v_uy, v_uz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s, t, u, v, ws, wt, wu, wv;
   __m512 cp, cr;
   __mmask16 msk, v_m;
   __attribute__((aligned(64))) unsigned int kk[16];
   __attribute__((aligned(64))) float scu[4*MXV*MYV*MZV];
/* __attribute__((aligned(64))) float scu[4*(mx+1)*(my+1)*(mz+1)]; */
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   ci2 = ci*ci;
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
   v_mxv4 = _mm512_set1_epi32(4*mxv);
   v_mxyv4 = _mm512_set1_epi32(4*mxyv);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qm = _mm512_set1_ps(qm);
   v_ci2 = _mm512_set1_ps(ci2);
   v_dt = _mm512_set1_ps(dt);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_edgelx = _mm512_set1_ps(edgelx);
   v_edgely = _mm512_set1_ps(edgely);
   v_edgelz = _mm512_set1_ps(edgelz);
   v_edgerx = _mm512_set1_ps(edgerx);
   v_edgery = _mm512_set1_ps(edgery);
   v_edgerz = _mm512_set1_ps(edgerz);
   v_at = _mm512_set_ps(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,
                        1.);
   v_m = _mm512_cmp_ps_mask(v_at,v_one,_MM_CMPINT_LT);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,nps,nn,mm,ll,ii,nm,lm,x,y,z, \
vx,vy,vz,ux,uy,uz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,p2,gami,v_noff, \
v_moff,v_loff,v_nn,v_mm,v_ll,v_it,v_x,v_y,v_z,v_dxp,v_dyp,v_dzp,v_amx, \
v_amy,v_amz,v_dx1,v_dx,v_dy,v_dz,v_vx,v_vy,v_vz,v_ux,v_uy,v_uz,v_gami, \
v_at,v_as,a,b,c,d,e,f,g,h,p,q,r,s,t,u,v,ws,wt,wu,wv,cp,cr,msk,kk,scu)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* zero out local accumulator */
/*    for (j = 0; j < 4*mxyv*(mz+1); j++) { */
/*       scu[j] = 0.0f;                     */
/*    }                                     */
      memset((void*)scu,0,4*mxyv*(mz+1)*sizeof(float));
      nps = 16*(npp/16);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = qm*(x - (float) nn); */
/*       dyp = y - (float) mm;      */
/*       dzp = z - (float) ll;      */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/* find inverse gamma */
/*       ux = ppart[j+3*nppmx+npoff]; */
/*       uy = ppart[j+4*nppmx+npoff]; */
/*       uz = ppart[j+5*nppmx+npoff]; */
         v_ux = _mm512_load_ps(&ppart[j+3*nppmx+npoff]);
         v_uy = _mm512_load_ps(&ppart[j+4*nppmx+npoff]);
         v_uz = _mm512_load_ps(&ppart[j+5*nppmx+npoff]);
/*       p2 = ux*ux + uy*uy + uz*uz;       */
         v_at = _mm512_fmadd_ps(v_uy,v_uy,_mm512_mul_ps(v_ux,v_ux));
         v_at = _mm512_fmadd_ps(v_uz,v_uz,v_at);
/*       gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_gami = _mm512_rsqrt23_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* full accuracy calculation */
         v_gami = _mm512_sqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one));
         v_gami = _mm512_div_ps(v_one,v_gami);
/* full accuracy calculation with SVML */
/*       v_gami = _mm512_invsqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* calculate weights */
/*       nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff)); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv4,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv4,v_mm));
         v_nn = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*       amx = qm - dxp;   */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_qm,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*       a = amx*amz; */
/*       b = amy*amz; */
/*       c = dyp*amz; */
/*       d = dx1*amz; */
         a = _mm512_mul_ps(v_amx,v_amz);
         b = _mm512_mul_ps(v_amy,v_amz);
         c = _mm512_mul_ps(v_dyp,v_amz);
         d = _mm512_mul_ps(v_dx1,v_amz);
/*       e = amx*dzp; */
/*       f = amy*dzp; */
/*       g = dyp*dzp; */
/*       h = dx1*dzp; */
         e = _mm512_mul_ps(v_amx,v_dzp);
         f = _mm512_mul_ps(v_amy,v_dzp);
         g = _mm512_mul_ps(v_dyp,v_dzp);
         h = _mm512_mul_ps(v_dx1,v_dzp);
/* deposit current */
/*       vx = ux*gami; */
/*       vy = uy*gami; */
/*       vz = uz*gami; */
         v_vx = _mm512_mul_ps(v_ux,v_gami);
         v_vy = _mm512_mul_ps(v_uy,v_gami);
         v_vz = _mm512_mul_ps(v_uz,v_gami);
         v_ll = _mm512_add_epi32(v_nn,v_mxyv4);
/* deposit charge for one particle at a time */
         for (i = 0; i < 16; i++) {
            ii = i >> 2;
            if (i==(ii<<2)) {
               switch (ii)
               {
                  case 0:
/* replicate velocities of first group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,0);
                     q = _mm512_permute4f128_ps(v_vy,0);
                     r = _mm512_permute4f128_ps(v_vz,0);
/* regroup weights for first group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(a,
                         _mm512_int2mask(61680),b,177);
                     t = _mm512_mask_permute4f128_ps(c,
                         _mm512_int2mask(61680),d,177);
                     u = _mm512_mask_permute4f128_ps(e,
                         _mm512_int2mask(61680),f,177);
                     v = _mm512_mask_permute4f128_ps(g,
                         _mm512_int2mask(61680),h,177);
                     break;
                  case 1:
/* replicate velocities of second group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,85);
                     q = _mm512_permute4f128_ps(v_vy,85);
                     r = _mm512_permute4f128_ps(v_vz,85);
/* regroup weights for second group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(b,
                         _mm512_int2mask(3855),a,177);
                     t = _mm512_mask_permute4f128_ps(d,
                         _mm512_int2mask(3855),c,177);
                     u = _mm512_mask_permute4f128_ps(f,
                         _mm512_int2mask(3855),e,177);
                     v = _mm512_mask_permute4f128_ps(h,
                         _mm512_int2mask(3855),g,177);
                     break;
                  case 2:
/* replicate velocities of third group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,170);
                     q = _mm512_permute4f128_ps(v_vy,170);
                     r = _mm512_permute4f128_ps(v_vz,170);
/* regroup weights for third group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(a,
                         _mm512_int2mask(61680),b,177);
                     s = _mm512_permute4f128_ps(s,78);
                     t = _mm512_mask_permute4f128_ps(c,
                         _mm512_int2mask(61680),d,177);
                     t = _mm512_permute4f128_ps(t,78);
                     u = _mm512_mask_permute4f128_ps(e,
                         _mm512_int2mask(61680),f,177);
                     u = _mm512_permute4f128_ps(u,78);
                     v = _mm512_mask_permute4f128_ps(g,
                         _mm512_int2mask(61680),h,177);
                     v = _mm512_permute4f128_ps(v,78);
                     break;
                  case 3:
/* replicate velocities of fourth group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,255);
                     q = _mm512_permute4f128_ps(v_vy,255);
                     r = _mm512_permute4f128_ps(v_vz,255);
/* regroup weights for fourth group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(b,
                         _mm512_int2mask(3855),a,177);
                     s = _mm512_permute4f128_ps(s,78);
                     t = _mm512_mask_permute4f128_ps(d,
                         _mm512_int2mask(3855),c,177);
                     t = _mm512_permute4f128_ps(t,78);
                     u = _mm512_mask_permute4f128_ps(f,
                         _mm512_int2mask(3855),e,177);
                     u = _mm512_permute4f128_ps(u,78);
                     v = _mm512_mask_permute4f128_ps(h,
                         _mm512_int2mask(3855),g,177);
                     v = _mm512_permute4f128_ps(v,78);
                     break;
               }
            }
            v_it = _mm512_setzero_epi32();
            switch (i-(ii<<2))
            {
/* first particle */
               case 0:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                         _mm512_int2mask(170),(__m512i)q,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,78);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,0);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,0);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,0);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,0);
                  break;
/* second particle */
               case 1:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                         _mm512_int2mask(85),(__m512i)p,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,24);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,85);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,85);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,85);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,85);
                  break;
/* third particle */
               case 2:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                         _mm512_int2mask(170),(__m512i)q,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)r,
                         _mm512_int2mask(51),(__m512i)v_at,78);
/* reorder weights */
                 ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                        _mm512_int2mask(119),(__m512i)s,170);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,170);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,170);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,170);
                   break;
/* fourth particle */
               case 3:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                         _mm512_int2mask(85),(__m512i)p,177);
                  v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_at,78);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,177);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,255);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,255);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,255);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,255);
                  break;
            }
            _mm512_store_epi32(kk,v_nn);
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dx = amx*amz;        */
/*          dy = amy*amz;        */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          dx = dyp*amz;        */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            mm = kk[i];
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),ws,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cp);
/*          mm = nn + 4*mxv;                                 */
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dx = dyp*amz;        */
/*          dy = dx1*amz;        */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            mm = kk[i] + 4*mxv;
            cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
                 &scu[mm]);
            cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
                 &scu[mm+16]);
            cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wt,cr);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cr);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cr);
            _mm512_store_epi32(kk,v_ll);
/*          nn += 4*mxyv;                                    */
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dx = amx*dzp;        */
/*          dy = amy*dzp;        */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            mm = kk[i];
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wu,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cp);
/*          mm = nn + 4*mxv;                                 */
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dx = dyp*dzp;        */
/*          dy = dx1*dzp;        */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            mm = kk[i] + 4*mxv;
            cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
                 &scu[mm]);
            cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
                 &scu[mm+16]);
            cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wv,cr);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cr);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cr);
         }
/* advance position half a time-step */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
/*       dz = z + vz*dt; */
         v_dx = _mm512_fmadd_ps(v_vx,v_dt,v_x);
         v_dy = _mm512_fmadd_ps(v_vy,v_dt,v_y);
         v_dz = _mm512_fmadd_ps(v_vz,v_dt,v_z);
/* reflecting boundary conditions */
        if (ipbc==2) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+3*nppmx+npoff] = -ux;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            v_ux = _mm512_mask_sub_ps(v_ux,msk,v_zero,v_ux);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+3*nppmx+npoff],v_ux);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             ppart[j+4*nppmx+npoff] = -uy;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            v_uy = _mm512_mask_sub_ps(v_uy,msk,v_zero,v_uy);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+4*nppmx+npoff],v_uy);
/*          if ((dz < edgelz) || (dz >= edgerz)) { */
/*             dz = z;                             */
/*             ppart[j+5*nppmx+npoff] = -uz;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dz,v_edgerz,
                  _MM_CMPINT_GE));
            v_dz = _mm512_mask_blend_ps(msk,v_dz,v_z);
            v_uz = _mm512_mask_sub_ps(v_uz,msk,v_zero,v_uz);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+5*nppmx+npoff],v_uz);
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+3*nppmx+npoff] = -ux;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            v_ux = _mm512_mask_sub_ps(v_ux,msk,v_zero,v_ux);
/* write output if test result is true for any particle */
           if (msk)
               _mm512_store_ps(&ppart[j+3*nppmx+npoff],v_ux);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             ppart[j+4*nppmx+npoff] = -uy;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            v_uy = _mm512_mask_sub_ps(v_uy,msk,v_zero,v_uy);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+4*nppmx+npoff],v_uy);
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy;   */
/*       ppart[j+2*nppmx+npoff] = dz; */
         _mm512_store_ps(&ppart[j+npoff],v_dx);
         _mm512_store_ps(&ppart[j+nppmx+npoff],v_dy);
         _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_dz);
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
/* find inverse gamma */
         ux = ppart[j+3*nppmx+npoff];
         uy = ppart[j+4*nppmx+npoff];
         uz = ppart[j+5*nppmx+npoff];
         p2 = ux*ux + uy*uy + uz*uz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* calculate weights */
         nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
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
         vx = ux*gami;
         vy = uy*gami;
         vz = uz*gami;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*amz;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*amz;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         dx = amx*dzp;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
         dy = amy*dzp;
         nn += 4*mxyv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*dzp;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*dzp;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+3*nppmx+npoff] = -ux;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+4*nppmx+npoff] = -uy;
            }
            if ((dz < edgelz) || (dz >= edgerz)) {
               dz = z;
               ppart[j+5*nppmx+npoff] = -uz;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+3*nppmx+npoff] = -ux;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+4*nppmx+npoff] = -uy;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      nps = 4*(nn/4);
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 1; i < nn; i++) {                     */
/*             cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]   */
/*             += scu[4*(i+mxv*j+mxyv*k)];                 */
/*             cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[1+4*(i+mxv*j+mxyv*k)];               */
/*             cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[2+4*(i+mxv*j+mxyv*k)];               */
/*          }                                              */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + mxv*j + mxyv*k);
               v_as = _mm512_loadunpacklo_ps(v_as,&scu[m]);
               v_as = _mm512_loadunpackhi_ps(v_as,&scu[m+16]);
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&cu[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&cu[m+16]);
/* skip add for first elements for i = 0 */
               if (i==0)
                  v_at = _mm512_mask_add_ps(v_at,v_m,v_at,v_as);
               else
                  v_at = _mm512_add_ps(v_at,v_as);
               _mm512_packstorelo_ps(&cu[m],v_at);
               _mm512_packstorehi_ps(&cu[m+16],v_at);
            }
/* loop over remaining elements */
            m = 1 > nps ? 1 : nps;
            for (i = m; i < nn; i++) {
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(i+mxv*j+mxyv*k)];
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*j+mxyv*k)];
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*j+mxyv*k)];
            }
         }
      }
/* deposit current to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*(j+moff)+nxyv*loff)] += scu[4*(i+mxv*j)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[1+4*(i+mxv*j)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[2+4*(i+mxv*j)];
            if (lm > mz) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*j+mxyv*(lm-1))];
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
            cu[4*(i+noff+nxv*moff+nxyv*(k+loff))] += scu[4*(i+mxyv*k)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[1+4*(i+mxyv*k)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[2+4*(i+mxyv*k)];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*k)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[1+4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[2+4*(mxv*j+mxyv*k)];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(nm-1+mxv*j+mxyv*k)];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[1+4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[2+4*(i+mxyv*(lm-1))];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[1+4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[2+4*(mxv*j+mxyv*(lm-1))];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(nm-1+mxv*j+mxyv*(lm-1))];
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
void ckncgrjppostf3lt(float ppart[], float cu[], int kpic[], int ncl[],
                      int ihole[], float qm, float dt, float ci,
                      int nppmx, int idimp, int nx, int ny, int nz,
                      int mx, int my, int mz, int nxv, int nyv, int nzv,
                      int mx1, int my1, int mxyz1, int ntmax,
                      int *irc) {
/* for 3d code, this subroutine calculates particle current density
   using first-order linear interpolation with periodic boundary
   conditions for relativistic particles.
   in addition, particle positions are advanced a half time-step
   also determines list of particles which are leaving this tile
   OpenMP/vector version using guard cells
   data deposited in tiles
   particles stored segmented array
   79 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
   input: all except ncl, ihole, irc
   output: ppart, cu, ncl, ihole, ek, irc
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
   and qci = qm*pi*gami, where i = x,y,z
   where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = momentum vx of particle n in tile m
   ppart[m][4][n] = momentum vy of particle n in tile m
   ppart[m][5][n] = momentum vz of particle n in tile m
   cu[l][k][j][i] = ith component of current density at grid point j,k,l
   kpic[l] = number of particles in tile l
   ncl[l][i] = number of particles going to destination i, tile l
   ihole[l][:][0] = location of hole in array left by departing particle
   ihole[l][:][1] = direction destination of particle leaving hole
   all for tile l
   ihole[l][0][0] = ih, number of holes left (error, if negative)
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
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
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   requires KNC, part needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
   cu needs to have 4 components, although one is not used
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, ih, nh, nn, mm, ll, ii, nm, lm, mxv, myv, mxyv;
   int nxyv;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz;
   float ci2, p2, gami;
   float x, y, z, ux, uy, uz;
   __m512i v_noff, v_moff, v_loff, v_mxv4, v_mxyv4;
   __m512i v_nn, v_mm, v_ll, v_it, v_0, v_1, v_3, v_9, v_perm;
   __m512 v_qm, v_ci2, v_dt, v_one, v_zero, v_anx, v_any, v_anz;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_gami, v_at, v_as, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_ux, v_uy, v_uz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s, t, u, v, ws, wt, wu, wv;
   __m512 cp, cr;
   __mmask16 msk, msk1, msk2, v_m;
   __attribute__((aligned(64))) unsigned int kk[16];
   __attribute__((aligned(64))) float scu[4*MXV*MYV*MZV];
/* __attribute__((aligned(64))) float scu[4*(mx+1)*(my+1)*(mz+1)]; */
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   ci2 = ci*ci;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
/* set boundary values */
   v_mxv4 = _mm512_set1_epi32(4*mxv);
   v_mxyv4 = _mm512_set1_epi32(4*mxyv);
   v_0 = _mm512_set1_epi32(0);
   v_1 = _mm512_set1_epi32(1);
   v_3 = _mm512_set1_epi32(3);
   v_9 = _mm512_set1_epi32(9);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qm = _mm512_set1_ps(qm);
   v_ci2 = _mm512_set1_ps(ci2);
   v_dt = _mm512_set1_ps(dt);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_anx = _mm512_set1_ps(anx);
   v_any = _mm512_set1_ps(any);
   v_anz = _mm512_set1_ps(anz);
   v_at = _mm512_set_ps(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,
                        1.);
   v_m = _mm512_cmp_ps_mask(v_at,v_one,_MM_CMPINT_LT);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,nps,nn,mm,ll,ii,nm,lm,ih,nh, \
x,y,z,vx,vy,vz,ux,uy,uz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,edgelx, \
edgely,edgelz,edgerx,edgery,edgerz,p2,gami,v_noff,v_moff,v_loff,v_nn, \
v_mm,v_ll,v_it,v_x,v_y,v_z,v_dxp,v_dyp,v_dzp,v_amx,v_amy,v_amz,v_dx1, \
v_dx,v_dy,v_dz,v_vx,v_vy,v_vz,v_ux,v_uy,v_uz,v_edgelx,v_edgely, \
v_edgelz,v_edgerx,v_edgery,v_edgerz,v_gami,v_at,v_as,a,b,c,d,e,f,g,h, \
p,q,r,s,t,u,v,ws,wt,wu,wv,cp,cr,msk,msk1,msk2,kk,scu)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
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
      v_edgelx = _mm512_set1_ps(edgelx);
      v_edgely = _mm512_set1_ps(edgely);
      v_edgelz = _mm512_set1_ps(edgelz);
      v_edgerx = _mm512_set1_ps(edgerx);
      v_edgery = _mm512_set1_ps(edgery);
      v_edgerz = _mm512_set1_ps(edgerz);
      ih = 0;
      nh = 0;
      nn += 1;
      mm += 1;
      ll += 1;
/* zero out local accumulator */
/*    for (j = 0; j < 4*mxyv*(mz+1); j++) { */
/*       scu[j] = 0.0f;                     */
/*    }                                     */
      memset((void*)scu,0,4*mxyv*(mz+1)*sizeof(float));
/* clear counters */
/*    for (j = 0; j < 26; j++) { */
/*       ncl[j+26*l] = 0;        */
/*    }                         */
      memset((void*)&ncl[26*l],0,26*sizeof(int));
      nps = 16*(npp/16);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = qm*(x - (float) nn); */
/*       dyp = y - (float) mm;      */
/*       dzp = z - (float) ll;      */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/* find inverse gamma */
/*       ux = ppart[j+3*nppmx+npoff]; */
/*       uy = ppart[j+4*nppmx+npoff]; */
/*       uz = ppart[j+5*nppmx+npoff]; */
         v_ux = _mm512_load_ps(&ppart[j+3*nppmx+npoff]);
         v_uy = _mm512_load_ps(&ppart[j+4*nppmx+npoff]);
         v_uz = _mm512_load_ps(&ppart[j+5*nppmx+npoff]);
/*       p2 = ux*ux + uy*uy + uz*uz;       */
         v_at = _mm512_fmadd_ps(v_uy,v_uy,_mm512_mul_ps(v_ux,v_ux));
         v_at = _mm512_fmadd_ps(v_uz,v_uz,v_at);
/*       gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_gami = _mm512_rsqrt23_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* full accuracy calculation */
         v_gami = _mm512_sqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one));
         v_gami = _mm512_div_ps(v_one,v_gami);
/* full accuracy calculation with SVML */
/*       v_gami = _mm512_invsqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* calculate weights */
/*       nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff)); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv4,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv4,v_mm));
         v_nn = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*       amx = qm - dxp;   */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_qm,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*       a = amx*amz; */
/*       b = amy*amz; */
/*       c = dyp*amz; */
/*       d = dx1*amz; */
         a = _mm512_mul_ps(v_amx,v_amz);
         b = _mm512_mul_ps(v_amy,v_amz);
         c = _mm512_mul_ps(v_dyp,v_amz);
         d = _mm512_mul_ps(v_dx1,v_amz);
/*       e = amx*dzp; */
/*       f = amy*dzp; */
/*       g = dyp*dzp; */
/*       h = dx1*dzp; */
         e = _mm512_mul_ps(v_amx,v_dzp);
         f = _mm512_mul_ps(v_amy,v_dzp);
         g = _mm512_mul_ps(v_dyp,v_dzp);
         h = _mm512_mul_ps(v_dx1,v_dzp);
/* deposit current */
/*       vx = ux*gami; */
/*       vy = uy*gami; */
/*       vz = uz*gami; */
         v_vx = _mm512_mul_ps(v_ux,v_gami);
         v_vy = _mm512_mul_ps(v_uy,v_gami);
         v_vz = _mm512_mul_ps(v_uz,v_gami);
         v_ll = _mm512_add_epi32(v_nn,v_mxyv4);
/* deposit charge for one particle at a time */
         for (i = 0; i < 16; i++) {
            ii = i >> 2;
            if (i==(ii<<2)) {
               switch (ii)
               {
                  case 0:
/* replicate velocities of first group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,0);
                     q = _mm512_permute4f128_ps(v_vy,0);
                     r = _mm512_permute4f128_ps(v_vz,0);
/* regroup weights for first group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(a,
                         _mm512_int2mask(61680),b,177);
                     t = _mm512_mask_permute4f128_ps(c,
                         _mm512_int2mask(61680),d,177);
                     u = _mm512_mask_permute4f128_ps(e,
                         _mm512_int2mask(61680),f,177);
                     v = _mm512_mask_permute4f128_ps(g,
                         _mm512_int2mask(61680),h,177);
                     break;
                  case 1:
/* replicate velocities of second group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,85);
                     q = _mm512_permute4f128_ps(v_vy,85);
                     r = _mm512_permute4f128_ps(v_vz,85);
/* regroup weights for second group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(b,
                         _mm512_int2mask(3855),a,177);
                     t = _mm512_mask_permute4f128_ps(d,
                         _mm512_int2mask(3855),c,177);
                     u = _mm512_mask_permute4f128_ps(f,
                         _mm512_int2mask(3855),e,177);
                     v = _mm512_mask_permute4f128_ps(h,
                         _mm512_int2mask(3855),g,177);
                     break;
                  case 2:
/* replicate velocities of third group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,170);
                     q = _mm512_permute4f128_ps(v_vy,170);
                     r = _mm512_permute4f128_ps(v_vz,170);
/* regroup weights for third group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(a,
                         _mm512_int2mask(61680),b,177);
                     s = _mm512_permute4f128_ps(s,78);
                     t = _mm512_mask_permute4f128_ps(c,
                         _mm512_int2mask(61680),d,177);
                     t = _mm512_permute4f128_ps(t,78);
                     u = _mm512_mask_permute4f128_ps(e,
                         _mm512_int2mask(61680),f,177);
                     u = _mm512_permute4f128_ps(u,78);
                     v = _mm512_mask_permute4f128_ps(g,
                         _mm512_int2mask(61680),h,177);
                     v = _mm512_permute4f128_ps(v,78);
                     break;
                  case 3:
/* replicate velocities of fourth group of 4 particles */
                     p = _mm512_permute4f128_ps(v_vx,255);
                     q = _mm512_permute4f128_ps(v_vy,255);
                     r = _mm512_permute4f128_ps(v_vz,255);
/* regroup weights for fourth group of 4 particles */
                     s = _mm512_mask_permute4f128_ps(b,
                         _mm512_int2mask(3855),a,177);
                     s = _mm512_permute4f128_ps(s,78);
                     t = _mm512_mask_permute4f128_ps(d,
                         _mm512_int2mask(3855),c,177);
                     t = _mm512_permute4f128_ps(t,78);
                     u = _mm512_mask_permute4f128_ps(f,
                         _mm512_int2mask(3855),e,177);
                     u = _mm512_permute4f128_ps(u,78);
                     v = _mm512_mask_permute4f128_ps(h,
                         _mm512_int2mask(3855),g,177);
                     v = _mm512_permute4f128_ps(v,78);
                     break;
               }
            }
            v_it = _mm512_setzero_epi32();
            switch (i-(ii<<2))
            {
/* first particle */
               case 0:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                         _mm512_int2mask(170),(__m512i)q,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,78);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,0);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,0);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,0);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,0);
                  break;
/* second particle */
               case 1:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                         _mm512_int2mask(85),(__m512i)p,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,24);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,85);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,85);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,85);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,85);
                  break;
/* third particle */
               case 2:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                         _mm512_int2mask(170),(__m512i)q,177);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)r,
                         _mm512_int2mask(51),(__m512i)v_at,78);
/* reorder weights */
                 ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                        _mm512_int2mask(119),(__m512i)s,170);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,170);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,170);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,170);
                   break;
/* fourth particle */
               case 3:
/* reorder velocity components */
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                         _mm512_int2mask(85),(__m512i)p,177);
                  v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_at,78);
                  v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                         _mm512_int2mask(68),(__m512i)r,177);
/* reorder weights */
                  ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)s,255);
                  wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)t,255);
                  wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)u,255);
                  wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                       _mm512_int2mask(119),(__m512i)v,255);
                  break;
            }
            _mm512_store_epi32(kk,v_nn);
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dx = amx*amz;        */
/*          dy = amy*amz;        */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          dx = dyp*amz;        */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            mm = kk[i];
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),ws,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cp);
/*          mm = nn + 4*mxv;                                 */
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dx = dyp*amz;        */
/*          dy = dx1*amz;        */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            mm = kk[i] + 4*mxv;
            cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
                 &scu[mm]);
            cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
                 &scu[mm+16]);
            cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wt,cr);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cr);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cr);
            _mm512_store_epi32(kk,v_ll);
/*          nn += 4*mxyv;                                    */
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dx = amx*dzp;        */
/*          dy = amy*dzp;        */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            mm = kk[i];
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wu,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cp);
/*          mm = nn + 4*mxv;                                 */
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dx = dyp*dzp;        */
/*          dy = dx1*dzp;        */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            mm = kk[i] + 4*mxv;
            cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
                 &scu[mm]);
            cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
                 &scu[mm+16]);
            cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wv,cr);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),cr);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),cr);
         }
/* advance position half a time-step */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
/*       dz = z + vz*dt; */
         v_dx = _mm512_fmadd_ps(v_vx,v_dt,v_x);
         v_dy = _mm512_fmadd_ps(v_vy,v_dt,v_y);
         v_dz = _mm512_fmadd_ps(v_vz,v_dt,v_z);
/* find particles going out of bounds */
/*       mm = 0; */
         v_mm = _mm512_setzero_epi32();
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* mm = direction particle is going                              */
/*       if (dx >= edgerx) {              */
/*          if (dx >= anx)                */
/*             ppart[j+npoff] = dx - anx; */
/*          mm = 2;                       */
/*       }                                */
         msk1 = _mm512_cmp_ps_mask(v_dx,v_edgerx,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dx;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_1,v_1);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dx,v_anx,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dx,v_anx);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dx = v_x;
            }
/*          if (dx < edgelx) {         */
/*             if (dx < 0.0) {         */
/*                dx += anx;           */
/*                if (dx < anx)        */
/*                   mm = 1;           */
/*                else                 */
/*                   dx = 0.0;         */
/*                ppart[j+npoff] = dx; */
/*             }                       */
/*             else {                  */
/*                mm = 1;              */
/*             }                       */
/*          }                          */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_1);
               msk2 = _mm512_cmp_ps_mask(v_dx,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dx,v_anx);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anx,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dx = v_x;
            }
         }
/*       if (dy >= edgery) {                    */
/*          if (dy >= any)                      */
/*             ppart[j+nppmx+npoff] = dy - any; */
/*          mm += 6;                            */
/*       }                                      */
         msk1 = _mm512_cmp_ps_mask(v_dy,v_edgery,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dy;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_3,v_3);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dy,v_any,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dy,v_any);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dy = v_x;
            }
/*          if (dy < edgely) {               */
/*             if (dy < 0.0) {               */
/*                dy += any;                 */
/*                if (dy < any)              */
/*                   mm += 3;                */
/*                else                       */
/*                   dy = 0.0;               */
/*                ppart[j+nppmx+npoff] = dy; */
/*             }                             */
/*             else {                        */
/*                mm += 3;                   */
/*             }                             */
/*          }                                */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_3);
               msk2 = _mm512_cmp_ps_mask(v_dy,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dy,v_any);
               msk1 = _mm512_cmp_ps_mask(v_x,v_any,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dy = v_x;
            }
         }
/*       if (dz >= edgerz) {                      */
/*          if (dz >= anz)                        */
/*             ppart[j+2*nppmx+npoff] = dz - anz; */
/*          mm += 18;                             */
/*       }                                        */
         msk1 = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dz;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_9,v_9);
               v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dz,v_anz,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dz,v_anz);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  v_dz = v_x;
            }
/*          if (dz < edgelz) {                 */
/*             if (dz < 0.0) {                 */
/*                dz += anz;                   */
/*                if (dz < anz)                */
/*                   mm += 9;                  */
/*                else                         */
/*                   dz = 0.0;                 */
/*                ppart[j+2*nppmx+npoff] = dz; */
/*             }                               */
/*             else {                          */
/*                mm += 9;                     */
/*             }                               */
/*          }                                  */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_9);
               msk2 = _mm512_cmp_ps_mask(v_dz,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dz,v_anz);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anz,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_mm = _mm512_add_epi32(v_mm,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  v_dz = v_x;
            }
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy;   */
/*       ppart[j+2*nppmx+npoff] = dz; */
         _mm512_store_ps(&ppart[j+npoff],v_dx);
         _mm512_store_ps(&ppart[j+nppmx+npoff],v_dy);
         _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_dz);
/* increment counters */
/*       if (mm > 0) {                                */
/*          ncl[mm+26*l-1] += 1;                      */
/*          ih += 1;                                  */
/*          if (ih <= ntmax) {                        */
/*             ihole[2*(ih+(ntmax+1)*l)] = j + i + 1; */
/*             ihole[1+2*(ih+(ntmax+1)*l)] = mm;      */
/*          }                                         */
/*          else {                                    */
/*             nh = 1;                                */
/*          }                                         */
/*       }                                            */
         _mm512_store_epi32(kk,v_mm);
         for (i = 0; i < 16; i++) {
            mm = kk[i];
            if (mm > 0) {
               ncl[mm+26*l-1] += 1;
               ih += 1; 
               if (ih <= ntmax) { 
                  ihole[2*(ih+(ntmax+1)*l)] = j + i + 1;
                  ihole[1+2*(ih+(ntmax+1)*l)] = mm;
               } 
               else {
                  nh = 1;
               }
            }
         }
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
/* find inverse gamma */
         ux = ppart[j+3*nppmx+npoff];
         uy = ppart[j+4*nppmx+npoff];
         uz = ppart[j+5*nppmx+npoff];
         p2 = ux*ux + uy*uy + uz*uz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* calculate weights */
         nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
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
         vx = ux*gami;
         vy = uy*gami;
         vz = uz*gami;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*amz;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*amz;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         dx = amx*dzp;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
         dy = amy*dzp;
         nn += 4*mxyv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*dzp;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*dzp;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
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
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
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
/* deposit current to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      nps = 4*(nn/4);
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 1; i < nn; i++) {                     */
/*             cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]   */
/*             += scu[4*(i+mxv*j+mxyv*k)];                 */
/*             cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[1+4*(i+mxv*j+mxyv*k)];               */
/*             cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[2+4*(i+mxv*j+mxyv*k)];               */
/*          }                                              */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + mxv*j + mxyv*k);
               v_as = _mm512_loadunpacklo_ps(v_as,&scu[m]);
               v_as = _mm512_loadunpackhi_ps(v_as,&scu[m+16]);
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&cu[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&cu[m+16]);
/* skip add for first elements for i = 0 */
               if (i==0)
                  v_at = _mm512_mask_add_ps(v_at,v_m,v_at,v_as);
               else
                  v_at = _mm512_add_ps(v_at,v_as);
               _mm512_packstorelo_ps(&cu[m],v_at);
               _mm512_packstorehi_ps(&cu[m+16],v_at);
            }
/* loop over remaining elements */
            m = 1 > nps ? 1 : nps;
            for (i = m; i < nn; i++) {
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(i+mxv*j+mxyv*k)];
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*j+mxyv*k)];
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*j+mxyv*k)];
            }
         }
      }
/* deposit current to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*(j+moff)+nxyv*loff)] += scu[4*(i+mxv*j)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[1+4*(i+mxv*j)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[2+4*(i+mxv*j)];
            if (lm > mz) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*j+mxyv*(lm-1))];
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
            cu[4*(i+noff+nxv*moff+nxyv*(k+loff))] += scu[4*(i+mxyv*k)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[1+4*(i+mxyv*k)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[2+4*(i+mxyv*k)];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*k)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[1+4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[2+4*(mxv*j+mxyv*k)];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(nm-1+mxv*j+mxyv*k)];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[1+4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[2+4*(i+mxyv*(lm-1))];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[1+4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[2+4*(mxv*j+mxyv*(lm-1))];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(nm-1+mxv*j+mxyv*(lm-1))];
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
   return;
#undef MXV
#undef MYV
#undef MZV
}

/*--------------------------------------------------------------------*/
void cknc2gjppost3lt(float ppart[], float cu[], int kpic[], float qm,
                     float dt, int nppmx, int idimp, int nx, int ny,
                     int nz, int mx, int my, int mz, int nxv, int nyv,
                     int nzv, int mx1, int my1, int mxyz1, int ipbc) {
/* for 3d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   OpenMP/vector version using guard cells
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = velocity vx of particle n in tile m
   ppart[m][4][n] = velocity vy of particle n in tile m
   ppart[m][5][n] = velocity vz of particle n in tile m
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
   requires KNC, part needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
   cu needs to have 4 components, although one is not used
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, nn, mm, ll, ii, nm, lm, mxv, myv, mxyv, nxyv;
   float edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz;
   float x, y, z;
   __m512i v_noff, v_moff, v_loff, v_mxv4, v_mxyv4;
   __m512i v_nn, v_mm, v_ll, v_it;
   __m512 v_qm, v_dt, v_one, v_zero;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_at, v_as, v_dx, v_dy, v_dz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 cp;
   __mmask16 msk, v_m;
   __attribute__((aligned(64))) unsigned int kk[16];
   typedef union vfloat {float v[16]; __m512 v16;} vf;
   __attribute__((aligned(64))) float scu[4*MXV*MYV*MZV];
/* __attribute__((aligned(64))) float scu[4*(mx+1)*(my+1)*(mz+1)]; */
   vf vv[8], vc[3], vu;
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
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
   v_mxv4 = _mm512_set1_epi32(4*mxv);
   v_mxyv4 = _mm512_set1_epi32(4*mxyv);
   v_qm = _mm512_set1_ps(qm);
   v_dt = _mm512_set1_ps(dt);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_edgelx = _mm512_set1_ps(edgelx);
   v_edgely = _mm512_set1_ps(edgely);
   v_edgelz = _mm512_set1_ps(edgelz);
   v_edgerx = _mm512_set1_ps(edgerx);
   v_edgery = _mm512_set1_ps(edgery);
   v_edgerz = _mm512_set1_ps(edgerz);
   v_at = _mm512_set_ps(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,
                        1.);
   v_m = _mm512_cmp_ps_mask(v_at,v_one,_MM_CMPINT_LT);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,nps,nn,mm,ll,ii,nm,lm,x,y, \
z,vx,vy,vz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,v_noff,v_moff,v_loff, \
v_nn,v_mm,v_ll,v_it,v_x,v_y,v_z,v_dxp,v_dyp,v_dzp,v_amx,v_amy,v_amz, \
v_dx1,v_dx,v_dy,v_dz,v_at,v_as,cp,msk,kk,scu,vv,vc,vu)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* zero out local accumulator */
/*    for (j = 0; j < 4*mxyv*(mz+1); j++) { */
/*       scu[j] = 0.0f;                     */
/*    }                                     */
      memset((void*)scu,0,4*mxyv*(mz+1)*sizeof(float));
      nps = 16*(npp/16);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = qm*(x - (float) nn); */
/*       dyp = y - (float) mm;      */
/*       dzp = z - (float) ll;      */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*       nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff)); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv4,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv4,v_mm));
         v_nn = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*       amx = qm - dxp;   */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_qm,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*       a = amx*amz; */
/*       b = amy*amz; */
/*       c = dyp*amz; */
/*       d = dx1*amz; */
         vv[0].v16 = _mm512_mul_ps(v_amx,v_amz);
         vv[1].v16 = _mm512_mul_ps(v_amy,v_amz);
         vv[2].v16 = _mm512_mul_ps(v_dyp,v_amz);
         vv[3].v16 = _mm512_mul_ps(v_dx1,v_amz);
/*       e = amx*dzp; */
/*       f = amy*dzp; */
/*       g = dyp*dzp; */
/*       h = dx1*dzp; */
         vv[4].v16 = _mm512_mul_ps(v_amx,v_dzp);
         vv[5].v16 = _mm512_mul_ps(v_amy,v_dzp);
         vv[6].v16 = _mm512_mul_ps(v_dyp,v_dzp);
         vv[7].v16 = _mm512_mul_ps(v_dx1,v_dzp);
         _mm512_store_epi32(kk,v_nn);
/* deposit current */
/*       vx = ppart[j+3*nppmx+npoff]; */
/*       vy = ppart[j+4*nppmx+npoff]; */
/*       vz = ppart[j+5*nppmx+npoff]; */
         vc[0].v16 = _mm512_load_ps(&ppart[j+3*nppmx+npoff]);
         vc[1].v16 = _mm512_load_ps(&ppart[j+4*nppmx+npoff]);
         vc[2].v16 = _mm512_load_ps(&ppart[j+5*nppmx+npoff]);
/* deposit charge for one particle at a time */
         for (i = 0; i < 16; i++) {
            nn = kk[i];
            vu.v16 = _mm512_setzero_ps();
            vu.v[0] = vc[0].v[i];
            vu.v[1] = vc[1].v[i];
            vu.v[2] = vc[2].v[i];
            v_at = _mm512_permute4f128_ps(vu.v16,0);
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dx = amx*amz;         */
/*          dy = amy*amz;         */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          dx = dyp*amz;         */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            vu.v[0] = vv[0].v[i];
            vu.v[4] = vv[1].v[i];
            v_as =  (__m512)_mm512_shuffle_epi32 ((__m512i)vu.v16,0);
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[nn]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[nn+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),v_as,cp);
            _mm512_mask_packstorelo_ps(&scu[nn],_mm512_int2mask(255),
                                       cp);
            _mm512_mask_packstorehi_ps(&scu[nn+16],_mm512_int2mask(255),
                                       cp);
            mm = nn + 4*mxv;
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dy = dx1*amz;         */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          dx = amx*dzp;         */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
/*          dx = amx*dzp;         */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            vu.v[0] = vv[2].v[i];
            vu.v[4] = vv[3].v[i];
            v_as =  (__m512)_mm512_shuffle_epi32 ((__m512i)vu.v16,0);
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),v_as,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),
                                       cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),
                                       cp);
            nn += 4*mxyv;
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dy = amy*dzp;         */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          dx = dyp*dzp;         */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            vu.v[0] = vv[4].v[i];
            vu.v[4] = vv[5].v[i];
            v_as =  (__m512)_mm512_shuffle_epi32 ((__m512i)vu.v16,0);
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[nn]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[nn+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),v_as,cp);
            _mm512_mask_packstorelo_ps(&scu[nn],_mm512_int2mask(255),
                                       cp);
            _mm512_mask_packstorehi_ps(&scu[nn+16],_mm512_int2mask(255),
                                       cp);
            mm = nn + 4*mxv;
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dy = dx1*dzp;         */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            vu.v[0] = vv[6].v[i];
            vu.v[4] = vv[7].v[i];
            v_as =  (__m512)_mm512_shuffle_epi32 ((__m512i)vu.v16,0);
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),v_as,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),
                                       cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),
                                       cp);
         }
/* advance position half a time-step */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
/*       dz = z + vz*dt; */
         v_dx = _mm512_fmadd_ps(vc[0].v16,v_dt,v_x);
         v_dy = _mm512_fmadd_ps(vc[1].v16,v_dt,v_y);
         v_dz = _mm512_fmadd_ps(vc[2].v16,v_dt,v_z);
/* reflecting boundary conditions */
        if (ipbc==2) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+3*nppmx+npoff] = -vx;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            vc[0].v16 = _mm512_mask_sub_ps(vc[0].v16,msk,v_zero,vc[0].v16);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+3*nppmx+npoff],vc[0].v16);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             ppart[j+4*nppmx+npoff] = -vy;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            vc[1].v16 = _mm512_mask_sub_ps(vc[1].v16,msk,v_zero,vc[1].v16);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+4*nppmx+npoff],vc[1].v16);
/*          if ((dz < edgelz) || (dz >= edgerz)) { */
/*             dz = z;                             */
/*             ppart[j+5*nppmx+npoff] = -vz;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dz,v_edgerz,
                  _MM_CMPINT_GE));
            v_dz = _mm512_mask_blend_ps(msk,v_dz,v_z);
            vc[2].v16 = _mm512_mask_sub_ps(vc[2].v16,msk,v_zero,vc[2].v16);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+5*nppmx+npoff],vc[2].v16);
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+3*nppmx+npoff] = -vx;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            vc[0].v16 = _mm512_mask_sub_ps(vc[0].v16,msk,v_zero,vc[0].v16);
/* write output if test result is true for any particle */
           if (msk)
               _mm512_store_ps(&ppart[j+3*nppmx+npoff],vc[0].v16);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             ppart[j+4*nppmx+npoff] = -vy;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            vc[1].v16 = _mm512_mask_sub_ps(vc[1].v16,msk,v_zero,vc[1].v16);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+4*nppmx+npoff],vc[1].v16);
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy;   */
/*       ppart[j+2*nppmx+npoff] = dz; */
         _mm512_store_ps(&ppart[j+npoff],v_dx);
         _mm512_store_ps(&ppart[j+nppmx+npoff],v_dy);
         _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_dz);
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
         nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
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
         vx = ppart[j+3*nppmx+npoff];
         vy = ppart[j+4*nppmx+npoff];
         vz = ppart[j+5*nppmx+npoff];
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*amz;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*amz;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         dx = amx*dzp;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
         dy = amy*dzp;
         nn += 4*mxyv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*dzp;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*dzp;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+3*nppmx+npoff] = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+4*nppmx+npoff] = -vy;
            }
            if ((dz < edgelz) || (dz >= edgerz)) {
               dz = z;
               ppart[j+5*nppmx+npoff] = -vz;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+3*nppmx+npoff] = -vx;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+4*nppmx+npoff] = -vy;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      nps = 4*(nn/4);
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 1; i < nn; i++) {                     */
/*             cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]   */
/*             += scu[4*(i+mxv*j+mxyv*k)];                 */
/*             cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[1+4*(i+mxv*j+mxyv*k)];               */
/*             cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[2+4*(i+mxv*j+mxyv*k)];               */
/*          }                                              */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + mxv*j + mxyv*k);
               v_as = _mm512_loadunpacklo_ps(v_as,&scu[m]);
               v_as = _mm512_loadunpackhi_ps(v_as,&scu[m+16]);
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&cu[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&cu[m+16]);
/* skip add for first elements for i = 0 */
               if (i==0)
                  v_at = _mm512_mask_add_ps(v_at,v_m,v_at,v_as);
               else
                  v_at = _mm512_add_ps(v_at,v_as);
               _mm512_packstorelo_ps(&cu[m],v_at);
               _mm512_packstorehi_ps(&cu[m+16],v_at);
            }
/* loop over remaining elements */
            m = 1 > nps ? 1 : nps;
            for (i = m; i < nn; i++) {
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(i+mxv*j+mxyv*k)];
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*j+mxyv*k)];
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*j+mxyv*k)];
            }
         }
      }
/* deposit current to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*(j+moff)+nxyv*loff)] += scu[4*(i+mxv*j)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[1+4*(i+mxv*j)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[2+4*(i+mxv*j)];
            if (lm > mz) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*j+mxyv*(lm-1))];
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
            cu[4*(i+noff+nxv*moff+nxyv*(k+loff))] += scu[4*(i+mxyv*k)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[1+4*(i+mxyv*k)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[2+4*(i+mxyv*k)];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*k)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[1+4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[2+4*(mxv*j+mxyv*k)];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(nm-1+mxv*j+mxyv*k)];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[1+4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[2+4*(i+mxyv*(lm-1))];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[1+4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[2+4*(mxv*j+mxyv*(lm-1))];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(nm-1+mxv*j+mxyv*(lm-1))];
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
void cknc2grjppost3lt(float ppart[], float cu[], int kpic[], float qm,
                      float dt, float ci, int nppmx, int idimp, int nx,
                      int ny, int nz, int mx, int my, int mz, int nxv,
                      int nyv, int nzv, int mx1, int my1, int mxyz1,
                      int ipbc) {
/* for 3d code, this subroutine calculates particle current density
   using first-order linear interpolation for relativistic particles
   in addition, particle positions are advanced a half time-step
   OpenMP/vector version using guard cells
   data deposited in tiles
   particles stored segmented array
   79 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
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
   and qci = qm*pi*gami, where i = x,y,z
   where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppart[m][3][n] = momentum vx of particle n in tile m
   ppart[m][4][n] = momentum vy of particle n in tile m
   ppart[m][5][n] = momentum vz of particle n in tile m
   cu[l][k][j][i] = ith component of current density at grid point j,k,l
   kpic = number of particles per tile
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
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
   requires KNC, part needs to be 64 byte aligned
   nppmx needs to be a multiple of 16
   cu needs to have 4 components, although one is not used
local data                                                            */
#define MXV             17
#define MYV             17
#define MZV             17
   int mxy1, noff, moff, loff, npoff, npp, nps;
   int i, j, k, l, m, nn, mm, ll, ii, nm, lm, mxv, myv, mxyv, nxyv;
   float ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz;
   float p2, gami;
   float x, y, z, ux, uy, uz;
   __m512i v_noff, v_moff, v_loff, v_mxv4, v_mxyv4;
   __m512i v_nn, v_mm, v_ll, v_it;
   __m512 v_qm, v_ci2, v_dt, v_one, v_zero;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_gami, v_at, v_as, v_dx, v_dy, v_dz, v_ux, v_uy, v_uz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 cp;
   __mmask16 msk, v_m;
   __attribute__((aligned(64))) unsigned int kk[16];
   typedef union vfloat {float v[16]; __m512 v16;} vf;
   __attribute__((aligned(64))) float scu[4*MXV*MYV*MZV];
/* __attribute__((aligned(64))) float scu[4*(mx+1)*(my+1)*(mz+1)]; */
   vf vv[8], vc[3], vu;
   mxy1 = mx1*my1;
/* mxv = MXV; */
/* myv = MYV; */
   mxv = mx+1;
   myv = my+1;
   mxyv = mxv*myv;
   nxyv = nxv*nyv;
   ci2 = ci*ci;
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
   v_mxv4 = _mm512_set1_epi32(4*mxv);
   v_mxyv4 = _mm512_set1_epi32(4*mxyv);
   v_qm = _mm512_set1_ps(qm);
   v_ci2 = _mm512_set1_ps(ci2);
   v_dt = _mm512_set1_ps(dt);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_edgelx = _mm512_set1_ps(edgelx);
   v_edgely = _mm512_set1_ps(edgely);
   v_edgelz = _mm512_set1_ps(edgelz);
   v_edgerx = _mm512_set1_ps(edgerx);
   v_edgery = _mm512_set1_ps(edgery);
   v_edgerz = _mm512_set1_ps(edgerz);
   v_at = _mm512_set_ps(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,
                        1.);
   v_m = _mm512_cmp_ps_mask(v_at,v_one,_MM_CMPINT_LT);
/* error if local array is too small                */
/* if ((mx >= MXV) || (my >= MYV) || (mz >= MZV))   */
/*    return;                                       */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,m,noff,moff,loff,npp,npoff,nps,nn,mm,ll,ii,nm,lm,x,y,z, \
vx,vy,vz,ux,uy,uz,dxp,dyp,dzp,amx,amy,amz,dx1,dx,dy,dz,p2,gami,v_noff, \
v_moff,v_loff,v_nn,v_mm,v_ll,v_it,v_x,v_y,v_z,v_dxp,v_dyp,v_dzp,v_amx, \
v_amy,v_amz,v_dx1,v_dx,v_dy,v_dz,v_ux,v_uy,v_uz,v_gami,v_at,v_as,cp, \
msk,kk,scu,vv,vc,vu)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      v_noff = _mm512_set1_epi32(noff);
      v_moff = _mm512_set1_epi32(moff);
      v_loff = _mm512_set1_epi32(loff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
/* zero out local accumulator */
/*    for (j = 0; j < 4*mxyv*(mz+1); j++) { */
/*       scu[j] = 0.0f;                     */
/*    }                                     */
      memset((void*)scu,0,4*mxyv*(mz+1)*sizeof(float));
      nps = 16*(npp/16);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*       x = ppart[j+npoff];         */
/*       y = ppart[j+nppmx+npoff];   */
/*       z = ppart[j+2*nppmx+npoff]; */
         v_x = _mm512_load_ps(&ppart[j+npoff]);
         v_y = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_z = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/*       nn = x; */
/*       mm = y; */
/*       ll = z; */
         v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
         v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
                _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*       dxp = qm*(x - (float) nn); */
/*       dyp = y - (float) mm;      */
/*       dzp = z - (float) ll;      */
         v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
         v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dyp = _mm512_sub_ps(v_y,v_dyp);
         v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dzp = _mm512_sub_ps(v_z,v_dzp);
/* find inverse gamma */
/*       ux = ppart[j+3*nppmx+npoff]; */
/*       uy = ppart[j+4*nppmx+npoff]; */
/*       uz = ppart[j+5*nppmx+npoff]; */
         v_ux = _mm512_load_ps(&ppart[j+3*nppmx+npoff]);
         v_uy = _mm512_load_ps(&ppart[j+4*nppmx+npoff]);
         v_uz = _mm512_load_ps(&ppart[j+5*nppmx+npoff]);
/*       p2 = ux*ux + uy*uy + uz*uz;       */
         v_at = _mm512_fmadd_ps(v_uy,v_uy,_mm512_mul_ps(v_ux,v_ux));
         v_at = _mm512_fmadd_ps(v_uz,v_uz,v_at);
/*       gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*       v_gami = _mm512_rsqrt23_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* full accuracy calculation */
         v_gami = _mm512_sqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one));
         v_gami = _mm512_div_ps(v_one,v_gami);
/* full accuracy calculation with SVML */
/*       v_gami = _mm512_invsqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* calculate weights */
/*       nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff)); */
         v_nn = _mm512_sub_epi32(v_nn,v_noff);
         v_mm = _mm512_sub_epi32(v_mm,v_moff);
         v_ll = _mm512_sub_epi32(v_ll,v_loff);
         v_it = _mm512_mullo_epi32(v_mxyv4,v_ll);
         v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_mxv4,v_mm));
         v_nn = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*       amx = qm - dxp;   */
/*       amy = 1.0f - dyp; */
/*       amz = 1.0f - dzp; */
         v_amx = _mm512_sub_ps(v_qm,v_dxp);
         v_amy = _mm512_sub_ps(v_one,v_dyp);
         v_amz = _mm512_sub_ps(v_one,v_dzp);
/*       dx1 = dxp*dyp; */
/*       dyp = amx*dyp; */
/*       amx = amx*amy; */
/*       amy = dxp*amy; */
         v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
         v_dyp = _mm512_mul_ps(v_amx,v_dyp);
         v_amx = _mm512_mul_ps(v_amx,v_amy);
         v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*       a = amx*amz; */
/*       b = amy*amz; */
/*       c = dyp*amz; */
/*       d = dx1*amz; */
         vv[0].v16 = _mm512_mul_ps(v_amx,v_amz);
         vv[1].v16 = _mm512_mul_ps(v_amy,v_amz);
         vv[2].v16 = _mm512_mul_ps(v_dyp,v_amz);
         vv[3].v16 = _mm512_mul_ps(v_dx1,v_amz);
/*       e = amx*dzp; */
/*       f = amy*dzp; */
/*       g = dyp*dzp; */
/*       h = dx1*dzp; */
         vv[4].v16 = _mm512_mul_ps(v_amx,v_dzp);
         vv[5].v16 = _mm512_mul_ps(v_amy,v_dzp);
         vv[6].v16 = _mm512_mul_ps(v_dyp,v_dzp);
         vv[7].v16 = _mm512_mul_ps(v_dx1,v_dzp);
         _mm512_store_epi32(kk,v_nn);
/* deposit current */
/*       vx = ux*gami; */
/*       vy = uy*gami; */
/*       vz = uz*gami; */
         vc[0].v16 = _mm512_mul_ps(v_ux,v_gami);
         vc[1].v16 = _mm512_mul_ps(v_uy,v_gami);
         vc[2].v16 = _mm512_mul_ps(v_uz,v_gami);
         v_ll = _mm512_add_epi32(v_nn,v_mxyv4);
/* deposit charge for one particle at a time */
         for (i = 0; i < 16; i++) {
            nn = kk[i];
            vu.v16 = _mm512_setzero_ps();
            vu.v[0] = vc[0].v[i];
            vu.v[1] = vc[1].v[i];
            vu.v[2] = vc[2].v[i];
            v_at = _mm512_permute4f128_ps(vu.v16,0);
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dx = amx*amz;         */
/*          dy = amy*amz;         */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          dx = dyp*amz;         */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            vu.v[0] = vv[0].v[i];
            vu.v[4] = vv[1].v[i];
            v_as =  (__m512)_mm512_shuffle_epi32 ((__m512i)vu.v16,0);
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[nn]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[nn+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),v_as,cp);
            _mm512_mask_packstorelo_ps(&scu[nn],_mm512_int2mask(255),
                                       cp);
            _mm512_mask_packstorehi_ps(&scu[nn+16],_mm512_int2mask(255),
                                       cp);
            mm = nn + 4*mxv;
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dy = dx1*amz;         */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          dx = amx*dzp;         */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
/*          dx = amx*dzp;         */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            vu.v[0] = vv[2].v[i];
            vu.v[4] = vv[3].v[i];
            v_as =  (__m512)_mm512_shuffle_epi32 ((__m512i)vu.v16,0);
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),v_as,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),
                                       cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),
                                       cp);
            nn += 4*mxyv;
/* load scu[nn:nn+3] and scu[nn+4:nn+7] field components */
/*          dy = amy*dzp;         */
/*          scu[nn] += vx*dx;     */
/*          scu[nn+1] += vy*dx;   */
/*          scu[nn+2] += vz*dx;   */
/*          dx = dyp*dzp;         */
/*          scu[nn+4] += vx*dy;   */
/*          scu[nn+1+4] += vy*dy; */
/*          scu[nn+2+4] += vz*dy; */
            vu.v[0] = vv[4].v[i];
            vu.v[4] = vv[5].v[i];
            v_as =  (__m512)_mm512_shuffle_epi32 ((__m512i)vu.v16,0);
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[nn]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[nn+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),v_as,cp);
            _mm512_mask_packstorelo_ps(&scu[nn],_mm512_int2mask(255),
                                       cp);
            _mm512_mask_packstorehi_ps(&scu[nn+16],_mm512_int2mask(255),
                                       cp);
            mm = nn + 4*mxv;
/* load scu[mm:mm+3] and scu[mm+4:mm+7] field components */
/*          dy = dx1*dzp;         */
/*          scu[mm] += vx*dx;     */
/*          scu[mm+1] += vy*dx;   */
/*          scu[mm+2] += vz*dx;   */
/*          scu[mm+4] += vx*dy;   */
/*          scu[mm+1+4] += vy*dy; */
/*          scu[mm+2+4] += vz*dy; */
            vu.v[0] = vv[6].v[i];
            vu.v[4] = vv[7].v[i];
            v_as =  (__m512)_mm512_shuffle_epi32 ((__m512i)vu.v16,0);
            cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
                 &scu[mm]);
            cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
                 &scu[mm+16]);
            cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),v_as,cp);
            _mm512_mask_packstorelo_ps(&scu[mm],_mm512_int2mask(255),
                                       cp);
            _mm512_mask_packstorehi_ps(&scu[mm+16],_mm512_int2mask(255),
                                       cp);
         }
/* advance position half a time-step */
/*       dx = x + vx*dt; */
/*       dy = y + vy*dt; */
/*       dz = z + vz*dt; */
         v_dx = _mm512_fmadd_ps(vc[0].v16,v_dt,v_x);
         v_dy = _mm512_fmadd_ps(vc[1].v16,v_dt,v_y);
         v_dz = _mm512_fmadd_ps(vc[2].v16,v_dt,v_z);
/* reflecting boundary conditions */
        if (ipbc==2) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+3*nppmx+npoff] = -ux;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            v_ux = _mm512_mask_sub_ps(v_ux,msk,v_zero,v_ux);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+3*nppmx+npoff],v_ux);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             ppart[j+4*nppmx+npoff] = -uy;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            v_uy = _mm512_mask_sub_ps(v_uy,msk,v_zero,v_uy);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+4*nppmx+npoff],v_uy);
/*          if ((dz < edgelz) || (dz >= edgerz)) { */
/*             dz = z;                             */
/*             ppart[j+5*nppmx+npoff] = -uz;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dz,v_edgerz,
                  _MM_CMPINT_GE));
            v_dz = _mm512_mask_blend_ps(msk,v_dz,v_z);
            v_uz = _mm512_mask_sub_ps(v_uz,msk,v_zero,v_uz);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+5*nppmx+npoff],v_uz);
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
/*          if ((dx < edgelx) || (dx >= edgerx)) { */
/*             dx = x;                             */
/*             ppart[j+3*nppmx+npoff] = -ux;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
                  _MM_CMPINT_GE));
            v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
            v_ux = _mm512_mask_sub_ps(v_ux,msk,v_zero,v_ux);
/* write output if test result is true for any particle */
           if (msk)
               _mm512_store_ps(&ppart[j+3*nppmx+npoff],v_ux);
/*          if ((dy < edgely) || (dy >= edgery)) { */
/*             dy = y;                             */
/*             ppart[j+4*nppmx+npoff] = -uy;       */
/*          }                                      */
            msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
            msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
                  _MM_CMPINT_GE));
            v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
            v_uy = _mm512_mask_sub_ps(v_uy,msk,v_zero,v_uy);
/* write output if test result is true for any particle */
            if (msk)
               _mm512_store_ps(&ppart[j+4*nppmx+npoff],v_uy);
         }
/* set new position */
/*       ppart[j+npoff] = dx;       */
/*       ppart[j+nppmx+npoff] = dy;   */
/*       ppart[j+2*nppmx+npoff] = dz; */
         _mm512_store_ps(&ppart[j+npoff],v_dx);
         _mm512_store_ps(&ppart[j+nppmx+npoff],v_dy);
         _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_dz);
      }
/* loop over remaining particles */
      for (j = nps; j < npp; j++) {
/* find interpolation weights */
         x = ppart[j+npoff];
         y = ppart[j+nppmx+npoff];
         z = ppart[j+2*nppmx+npoff];
         nn = x;
         mm = y;
         ll = z;
         dxp = qm*(x - (float) nn);
         dyp = y - (float) mm;
         dzp = z - (float) ll;
/* find inverse gamma */
         ux = ppart[j+3*nppmx+npoff];
         uy = ppart[j+4*nppmx+npoff];
         uz = ppart[j+5*nppmx+npoff];
         p2 = ux*ux + uy*uy + uz*uz;
         gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* calculate weights */
         nn = 4*(nn - noff + mxv*(mm - moff) + mxyv*(ll - loff));
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
         vx = ux*gami;
         vy = uy*gami;
         vz = uz*gami;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*amz;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*amz;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         dx = amx*dzp;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
         dy = amy*dzp;
         nn += 4*mxyv;
         scu[nn] += vx*dx;
         scu[nn+1] += vy*dx;
         scu[nn+2] += vz*dx;
         dx = dyp*dzp;
         scu[nn+4] += vx*dy;
         scu[nn+1+4] += vy*dy;
         scu[nn+2+4] += vz*dy;
         dy = dx1*dzp;
         mm = nn + 4*mxv;
         scu[mm] += vx*dx;
         scu[mm+1] += vy*dx;
         scu[mm+2] += vz*dx;
         scu[mm+4] += vx*dy;
         scu[mm+1+4] += vy*dy;
         scu[mm+2+4] += vz*dy;
/* advance position half a time-step */
         dx = x + vx*dt;
         dy = y + vy*dt;
         dz = z + vz*dt;
/* reflecting boundary conditions */
         if (ipbc==2) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+3*nppmx+npoff] = -ux;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+4*nppmx+npoff] = -uy;
            }
            if ((dz < edgelz) || (dz >= edgerz)) {
               dz = z;
               ppart[j+5*nppmx+npoff] = -uz;
            }
         }
/* mixed reflecting/periodic boundary conditions */
         else if (ipbc==3) {
            if ((dx < edgelx) || (dx >= edgerx)) {
               dx = x;
               ppart[j+3*nppmx+npoff] = -ux;
            }
            if ((dy < edgely) || (dy >= edgery)) {
               dy = y;
               ppart[j+4*nppmx+npoff] = -uy;
            }
         }
/* set new position */
         ppart[j+npoff] = dx;
         ppart[j+nppmx+npoff] = dy;
         ppart[j+2*nppmx+npoff] = dz;
      }
/* deposit current to interior points in global array */
      nn = nxv - noff;
      nn = mx < nn ? mx : nn;
      mm = nyv - moff;
      mm = my < mm ? my : mm;
      ll = nzv - loff;
      ll = mz < ll ? mz : ll;
      nps = 4*(nn/4);
      for (k = 1; k < ll; k++) {
         for (j = 1; j < mm; j++) {
/* vector loop over elements in blocks of 4 */
/*          for (i = 1; i < nn; i++) {                     */
/*             cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]   */
/*             += scu[4*(i+mxv*j+mxyv*k)];                 */
/*             cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[1+4*(i+mxv*j+mxyv*k)];               */
/*             cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))] */
/*             += scu[2+4*(i+mxv*j+mxyv*k)];               */
/*          }                                              */
            for (i = 0; i < nps; i+=4) {
               m = 4*(i + mxv*j + mxyv*k);
               v_as = _mm512_loadunpacklo_ps(v_as,&scu[m]);
               v_as = _mm512_loadunpackhi_ps(v_as,&scu[m+16]);
               m = 4*(i + noff + nxv*(j + moff) + nxyv*(k + loff));
               v_at = _mm512_loadunpacklo_ps(v_at,&cu[m]);
               v_at = _mm512_loadunpackhi_ps(v_at,&cu[m+16]);
/* skip add for first elements for i = 0 */
               if (i==0)
                  v_at = _mm512_mask_add_ps(v_at,v_m,v_at,v_as);
               else
                  v_at = _mm512_add_ps(v_at,v_as);
               _mm512_packstorelo_ps(&cu[m],v_at);
               _mm512_packstorehi_ps(&cu[m+16],v_at);
            }
/* loop over remaining elements */
            m = 1 > nps ? 1 : nps;
            for (i = m; i < nn; i++) {
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(i+mxv*j+mxyv*k)];
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*j+mxyv*k)];
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*j+mxyv*k)];
            }
         }
      }
/* deposit current to edge points in global array */
      lm = nzv - loff;
      lm = mz+1 < lm ? mz+1 : lm;
      for (j = 1; j < mm; j++) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*(j+moff)+nxyv*loff)] += scu[4*(i+mxv*j)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[1+4*(i+mxv*j)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*(j+moff)+nxyv*loff)]
            += scu[2+4*(i+mxv*j)];
            if (lm > mz) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*j+mxyv*(lm-1))];
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
            cu[4*(i+noff+nxv*moff+nxyv*(k+loff))] += scu[4*(i+mxyv*k)];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[1+4*(i+mxyv*k)];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(k+loff))]
            += scu[2+4*(i+mxyv*k)];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*k)];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(k+loff))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*k)];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[1+4*(mxv*j+mxyv*k)];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(k+loff))]
            += scu[2+4*(mxv*j+mxyv*k)];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[1+4*(nm-1+mxv*j+mxyv*k)];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(k+loff))]
               += scu[2+4*(nm-1+mxv*j+mxyv*k)];
            }
         }
      }
      if (lm > mz) {
         for (i = 1; i < nn; i++) {
#pragma omp atomic
            cu[4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[1+4*(i+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(i+noff+nxv*moff+nxyv*(lm+loff-1))]
            += scu[2+4*(i+mxyv*(lm-1))];
            if (mm > my) {
#pragma omp atomic
               cu[4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[1+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(i+noff+nxv*(mm+moff-1)+nxyv*(lm+loff-1))]
               += scu[2+4*(i+mxv*(mm-1)+mxyv*(lm-1))];
            }
         }
         for (j = 0; j < mm; j++) {
#pragma omp atomic
            cu[4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[1+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[1+4*(mxv*j+mxyv*(lm-1))];
#pragma omp atomic
            cu[2+4*(noff+nxv*(j+moff)+nxyv*(lm+loff-1))]
            += scu[2+4*(mxv*j+mxyv*(lm-1))];
            if (nm > mx) {
#pragma omp atomic
               cu[4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[1+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[1+4*(nm-1+mxv*j+mxyv*(lm-1))];
#pragma omp atomic
               cu[2+4*(nm+noff-1+nxv*(j+moff)+nxyv*(lm+loff-1))]
               += scu[2+4*(nm-1+mxv*j+mxyv*(lm-1))];
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
void ckncpporder3lt(float ppart[], float ppbuff[], int kpic[],
                    int ncl[], int ihole[], int idimp, int nppmx, 
                    int nx, int ny, int nz, int mx, int my, int mz,
                    int mx1, int my1, int mz1, int npbmx, int ntmax,
                    int *irc) {
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppbuff[m][i][n] = i co-ordinate of particle n in tile m
   kpic[m] = number of particles in tile m
   ncl[m][i] = number of particles going to destination i, tile m
   ihole[m][:][0] = location of hole in array left by departing particle
   ihole[m][:][1] = direction destination of particle leaving hole
   all for tile m
   ihole[m][0][0] = ih, number of holes left (error, if negative)
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
   requires KNC, ppart, ppbuff need to be 64 byte aligned
   nppmx, npbmx need to be a multiple of 16
local data                                                            */
   int mxy1, mxyz1, noff, moff, loff, npoff, npp, nps, nboff, ncoff;
   int i, j, k, l, ii, kx, ky, kz, ih, nh, ist, nn, mm, ll;
   int ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dx, dy, dz;
   int ks[26];
   __m512i v_ist, v_it, v_0, v_1, v_3, v_9;
   __m512i v_m1, v_m2, v_m3, v_npp, v_mm, v_is, v_it0, v_ioff;
   __m512 v_anx, v_any, v_anz;
   __m512 v_dx, v_dy, v_dz, v_x;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 v_zero;
   __mmask16 msk1, msk2;
   __attribute__((aligned(64))) unsigned int ls[32], lm[32];
   mxy1 = mx1*my1;
   mxyz1 = mxy1*mz1;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
   v_0 = _mm512_set1_epi32(0);
   v_1 = _mm512_set1_epi32(1);
   v_3 = _mm512_set1_epi32(3);
   v_9 = _mm512_set1_epi32(9);
   v_anx = _mm512_set1_ps(anx);
   v_any = _mm512_set1_ps(any);
   v_anz = _mm512_set1_ps(anz);
   v_zero = _mm512_setzero_ps();
/* find and count particles leaving tiles and determine destination */
/* update ppart, ihole, ncl */
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,ii,noff,moff,loff,npp,npoff,nps,nn,mm,ll,ih,nh,ist,dx, \
dy,dz,edgelx,edgely,edgelz,edgerx,edgery,edgerz,v_it,v_ist,v_edgelx, \
v_edgely,v_edgelz,v_edgerx,v_edgery,v_edgerz,v_dx,v_dy,v_dz,v_x,msk1, \
msk2,ls)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
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
      noff = (ntmax+1)*l;
      v_edgelx = _mm512_set1_ps(edgelx);
      v_edgely = _mm512_set1_ps(edgely);
      v_edgelz = _mm512_set1_ps(edgelz);
      v_edgerx = _mm512_set1_ps(edgerx);
      v_edgery = _mm512_set1_ps(edgery);
      v_edgerz = _mm512_set1_ps(edgerz);
/* clear counters */
/*    for (j = 0; j < 26; j++) { */
/*       ncl[j+26*l] = 0;        */
/*    }                          */
      memset((void*)&ncl[26*l],0,26*sizeof(int));
      nps = 16*(npp/16);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/*       dx = ppart[j+npoff];         */
/*       dy = ppart[j+nppmx+npoff];   */
/*       dz = ppart[j+2*nppmx+npoff]; */
         v_dx = _mm512_load_ps(&ppart[j+npoff]);
         v_dy = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_dz = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/* find particles going out of bounds */
/*       ist = 0; */
         v_ist = _mm512_setzero_epi32();
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* ist = direction particle is going                             */
/*       if (dx >= edgerx) {              */
/*          if (dx >= anx)                */
/*             ppart[j+npoff] = dx - anx; */
/*          ist = 2;                      */
/*       }                                */
         msk1 = _mm512_cmp_ps_mask(v_dx,v_edgerx,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dx;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_1,v_1);
               v_ist = _mm512_mask_add_epi32(v_ist,msk1,v_ist,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dx,v_anx,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dx,v_anx);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+npoff],v_x);
            }
/*          if (dx < edgelx) {         */
/*             if (dx < 0.0) {         */
/*                dx += anx;           */
/*                if (dx < anx)        */
/*                   ist = 1;          */
/*                else                 */
/*                   dx = 0.0;         */
/*                ppart[j+npoff] = dx; */
/*             }                       */
/*             else {                  */
/*                ist = 1;             */
/*             }                       */
/*          }                          */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_1);
               msk2 = _mm512_cmp_ps_mask(v_dx,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dx,v_anx);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anx,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_ist = _mm512_add_epi32(v_ist,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+npoff],v_x);
            }
         }
/*       if (dy >= edgery) {                    */
/*          if (dy >= any)                      */
/*             ppart[j+nppmx+npoff] = dy - any; */
/*          ist += 6;                           */
/*       }                                      */
         msk1 = _mm512_cmp_ps_mask(v_dy,v_edgery,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dy;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_3,v_3);
               v_ist = _mm512_mask_add_epi32(v_ist,msk1,v_ist,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dy,v_any,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dy,v_any);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+nppmx+npoff],v_x);
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
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_3);
               msk2 = _mm512_cmp_ps_mask(v_dy,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dy,v_any);
               msk1 = _mm512_cmp_ps_mask(v_x,v_any,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_ist = _mm512_add_epi32(v_ist,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+nppmx+npoff],v_x);
            }
         }
/*       if (dz >= edgerz) {                      */
/*          if (dz >= anz)                        */
/*             ppart[j+2*nppmx+npoff] = dz - anz; */
/*          ist += 18;                            */
/*       }                                        */
         msk1 = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dz;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_9,v_9);
               v_ist = _mm512_mask_add_epi32(v_ist,msk1,v_ist,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dz,v_anz,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dz,v_anz);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_x);
            }
/*          if (dz < edgelz) {                 */
/*             if (dz < 0.0) {                 */
/*                dz += anz;                   */
/*                if (dz < anz)                */
/*                   ist += 9;                 */
/*                else                         */
/*                   dz = 0.0;                 */
/*                ppart[j+2*nppmx+npoff] = dz; */
/*             }                               */
/*             else {                          */
/*                ist += 9;                    */
/*             }                               */
/*          }                                  */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_9);
               msk2 = _mm512_cmp_ps_mask(v_dz,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dz,v_anz);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anz,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_ist = _mm512_add_epi32(v_ist,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_x);
            }
         }
/* increment counters */
/*       if (ist > 0) {                               */
/*          ncl[ist+26*l-1] += 1;                     */
/*          ih += 1;                                  */
/*          if (ih <= ntmax) {                        */
/*             ihole[2*(ih+(ntmax+1)*l)] = j + i + 1; */
/*             ihole[1+2*(ih+(ntmax+1)*l)] = ist;     */
/*          }                                         */
/*          else {                                    */
/*             nh = 1;                                */
/*          }                                         */
/*       }                                            */
         _mm512_store_epi32(ls,v_ist);
         for (i = 0; i < 16; i++) {
            ist = ls[i];
            if (ist > 0) {
               ncl[ist+26*l-1] += 1;
               ih += 1; 
               if (ih <= ntmax) { 
                  ihole[2*(ih+noff)] = j + i + 1;
                  ihole[1+2*(ih+noff)] = ist;
               } 
               else {
                  nh = 1;
               }
            }
         }
      }
/* loop over remaining particles in tile */
      for (j = nps; j < npp; j++) {
         dx = ppart[j+npoff];
         dy = ppart[j+nppmx+npoff];
         dz = ppart[j+2*nppmx+npoff];
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
         if (dz >= edgerz) {
            if (dz >= anz)
               ppart[j+2*nppmx+npoff] = dz - anz;
            ist += 18;
         }
         else if (dz < edgelz) {
            if (dz < 0.0) {
               dz += anz;
               if (dz < anz)
                  ist += 9;
               else
                  dz = 0.0;
               ppart[j+2*nppmx+npoff] = dz;
            }
            else {
               ist += 9;
            }
         }
         if (ist > 0) {
            ncl[ist+26*l-1] += 1;
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
   msk1 = _mm512_int2mask(1023);
   v_m1 = _mm512_set_epi32(11,11,11,11,11,10,9,8,3,3,3,3,3,2,1,0);
   v_m2 = _mm512_set_epi32(7,7,7,7,7,7,7,7,7,6,5,4,3,2,1,0);
#pragma omp parallel for \
private(i,j,l,npoff,nboff,noff,nps,mm,ii,ll,j1,ist,nh,ip,v_it,v_is, \
v_it0,v_ioff,ls,lm)
   for (l = 0; l < mxyz1; l++) {
      npoff = idimp*nppmx*l;
      nboff = idimp*npbmx*l;
      noff = (ntmax+1)*l;
/* find address offset for ordered ppbuff array */
/*    isum = 0;                  */
/*    for (j = 0; j < 26; j++) { */
/*       ist = ncl[j+26*l];      */
/*       ncl[j+26*l] = isum;     */
/*       isum += ist;            */
/*    }                          */
/* perform exclusive prefix scan */
/* load 26 data elements into 32 length vector with zero padding */
      mm = 26*l;
      v_it = _mm512_loadunpacklo_epi32(v_0,&ncl[mm]);
      v_it = _mm512_loadunpackhi_epi32(v_it,&ncl[mm+16]);
      _mm512_store_epi32(ls,v_it);
      v_is = _mm512_mask_loadunpacklo_epi32(v_0,msk1,&ncl[mm+16]);
      v_is = _mm512_mask_loadunpackhi_epi32(v_is,msk1,&ncl[mm+32]);
      _mm512_store_epi32(&ls[16],v_is);
      v_ioff = _mm512_setzero_epi32();
/* vector loop over elements in blocks of 16 */
      for (j = 0; j < 32; j+=16) {
/* load data */
         v_it0 = _mm512_load_epi32(&ls[j]);
/* first pass */
         v_is = _mm512_shuffle_epi32(v_it0,177);
         v_it = _mm512_mask_add_epi32(v_it0,_mm512_int2mask(43690),
                v_it0,v_is);
/* second pass */
         v_is = _mm512_shuffle_epi32(v_it,80);
         v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(52428),v_it,
                v_is);
/* third pass */
         v_is = _mm512_permutevar_epi32(v_m1,v_it);
         v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(61680),v_it,
                v_is);
/* fourth pass */
         v_is = _mm512_permutevar_epi32(v_m2,v_it);
         v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(65280),v_it,
                v_is);
/* add offset */
         v_it = _mm512_add_epi32(v_it,v_ioff);
/* next offset */
         if (j==0) {
            v_ioff = _mm512_shuffle_epi32(v_it,255);
            v_ioff = _mm512_permute4f128_epi32(v_ioff,255);
         }
/* subtract for exclusive scan */
         v_it = _mm512_sub_epi32(v_it,v_it0);
/* write data */
         _mm512_store_epi32(&ls[j],v_it);
      }
      nh = ihole[2*noff];
      nps = 16*(nh/16);
/*    nps = (nh >> 4) << 4; */
      ip = 0;
/* loop over particles leaving tile in groups of 16 */
      for (j = 0; j < nps; j+=16) {
/*       j1 = ihole[2*(j+1+(ntmax+1)*l)] - 1; */
/*       ist = ihole[1+2*(j+1+(ntmax+1)*l)];  */
         mm = 2*(j+1+noff);
         v_it = _mm512_loadunpacklo_epi32(v_0,&ihole[mm]);
         v_it = _mm512_loadunpackhi_epi32(v_it,&ihole[mm+16]);
         _mm512_store_epi32(lm,v_it);
         mm += 16;
         v_is = _mm512_loadunpacklo_epi32(v_0,&ihole[mm]);
         v_is = _mm512_loadunpackhi_epi32(v_is,&ihole[mm+16]);
         _mm512_store_epi32(&lm[16],v_is);
/* buffer particles that are leaving tile, in direction order */
         for (ll = 0; ll < 16; ll++) {
            j1 = lm[2*ll] - 1;
            ist = lm[1+2*ll];
            ii = ls[ist-1];
            if (ii < npbmx) {
               for (i = 0; i < idimp; i++) {
                  ppbuff[ii+npbmx*i+nboff]
                  = ppart[j1+nppmx*i+npoff];
               }
            }
            else {
               ip = 1;
            }
            ls[ist-1] = ii + 1;
         }
      }
/* loop over remaining particles leaving tile */
      for (j = nps; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+noff)] - 1;
         ist = ihole[1+2*(j+1+noff)];
         ii = ls[ist-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*i+nboff]
               = ppart[j1+nppmx*i+npoff];
            }
         }
         else {
            ip = 1;
         }
         ls[ist-1] = ii + 1;
      }
/* store 26 data elements into ncl */
      mm = 26*l;
      v_it = _mm512_load_epi32(ls);
      v_is = _mm512_load_epi32(&ls[16]);
     _mm512_packstorelo_epi32(&ncl[mm],v_it);
     _mm512_packstorehi_epi32(&ncl[mm+16],v_it);
     _mm512_mask_packstorelo_epi32(&ncl[mm+16],msk1,v_is);
     _mm512_mask_packstorehi_epi32(&ncl[mm+32],msk1,v_is);
/* set error */
      if (ip > 0)
         *irc = ncl[25+26*l];
   }
/* ppbuff overflow */
   if (*irc > 0)
      return;

/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
   v_ioff = _mm512_set_epi32(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
   v_m1 = _mm512_set1_epi32(nppmx);
#pragma omp parallel for \
private(i,j,k,l,ii,kk,npp,nps,npoff,noff,nboff,kx,ky,kz,kl,kr,kxl,kxr, \
lk,ll,lr,ih,nh,nn,mm,ncoff,ist,j1,j2,ip,v_m2,v_m3,v_it,v_is,v_it0,v_mm, \
v_npp,v_x,msk1,ks,ls)
   for (l = 0; l < mxyz1; l++) {
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      noff = (ntmax+1)*l;
      v_m2 = _mm512_set1_epi32(noff+1);
      v_m3 = _mm512_set1_epi32(npoff);
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
      nh = ihole[2*noff];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      v_it0 = _mm512_set1_epi32(nh);
      v_is = _mm512_add_epi32(v_m2,v_it0);
      v_it0 = _mm512_sub_epi32(v_ioff,v_it0);
      v_npp = _mm512_set1_epi32(npp);
      for (ii = 0; ii < 26; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+26*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+26*ks[ii]] - ncoff;
/*       nps = 16*(ip/16); */
         nps = (ip >> 4) << 4;
/* loop over particles in this direction in groups of 16 */
         for (j = 0; j < nps; j+=16) {
/* insert incoming particles into holes */
/*          ih += 1;                        */
/*          if (ih <= nh) {                 */
/*             j1 = ihole[2*(ih+noff)] - 1; */
/*          }                               */
/* place overflow at end of array */
/*          else {                          */
/*             j1 = npp;                    */
/*             npp += 1;                    */
/*          }                               */
            v_mm = _mm512_add_epi32(_mm512_set1_epi32(ih),v_it0);
            msk1 = _mm512_cmp_epi32_mask(v_mm,v_0,_MM_CMPINT_LT);
            v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_is);
            v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_mm);
            v_mm = _mm512_mask_add_epi32(v_mm,_mm512_knot(msk1),v_mm,
                                         v_npp);
            v_it = _mm512_mask_i32gather_epi32(v_mm,msk1,v_mm,
                                               (int *)ihole,4);
            v_it = _mm512_mask_sub_epi32(v_it,msk1,v_it,v_1);
            ih += 16;
            nn = ih - nh;
            if (nn > 0) {
               nn = nn < 16 ? nn : 16;
               npp += nn;
            }
            msk1 = _mm512_cmp_epi32_mask(v_it,v_m1,_MM_CMPINT_LT);
            ll = _mm512_mask2int(_mm512_knot(msk1));
            v_it = _mm512_add_epi32(v_it,v_m3);
            for (i = 0; i < idimp; i++) {
/*             if (j1 < nppmx)                     */
/*                ppart[j1+nppmx*i+npoff]          */
/*                = ppbuff[j+ncoff+npbmx*i+nboff]; */
               mm = j + ncoff + npbmx*i + nboff;
               v_x = _mm512_loadunpacklo_ps(v_x,&ppbuff[mm]);
               v_x = _mm512_loadunpackhi_ps(v_x,&ppbuff[mm+16]);
               if (ll==0) {
                  _mm512_i32scatter_ps((float *)ppart,v_it,v_x,4);
               }
               else {
                  _mm512_mask_i32scatter_ps((float *)ppart,msk1,v_it,
                                            v_x,4);
               }
               v_it = _mm512_add_epi32(v_it,v_m1);
            }
            if (ll != 0) {
               ist = 1;
            }
         }
/* loop over remaining particles in this direction */
         for (j = nps; j < ip; j++) {
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
/* holes with locations great than npp-ip do not need to be filled */
      if (ih < nh) {
         ip = nh - ih;
         ii = nh;
         nn = ihole[2*(ii+noff)] - 1;
         v_it0 = _mm512_set1_epi32(nn);
         ih += 1;
         j2 = ihole[2*(ih+noff)] - 1;
         v_m2 = _mm512_sub_epi32(v_m2,v_1);
/* move particles from end into remaining holes */
/* holes are processed in increasing order */
/*       nps = 16*(ip/16); */
         nps = (ip >> 4) << 4;
/* loop over particles in groups of 16 */
         for (j = 0; j < nps; j+=16) {
/*          j2 = ihole[2*(ih+noff)] - 1; */
            v_mm = _mm512_add_epi32(_mm512_set1_epi32(ih),v_ioff);
            v_mm = _mm512_add_epi32(v_mm,v_m2);
            v_mm = _mm512_add_epi32(v_mm,v_mm);
            v_is = _mm512_i32gather_epi32(v_mm,(int *)ihole,4);
            v_is = _mm512_sub_epi32(v_is,v_1);
/*          j1 = npp - j - 1;                      */
/*          if (j1==nn) {                          */
/*             ii -= 1;                            */
/*             nn = ihole[2*(ii+(ntmax+1)*l)] - 1; */
/*          }                                      */
            kk = 0;
            for (ll = 0; ll < 16; ll++) {
               j1 = npp - j - ll - 1;
               if (j1==nn) {
                  ii -= 1;
                  nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
               }
               else {
                  ls[kk] = j1;
                  kk += 1;
               }
            }
            v_it = _mm512_load_epi32(ls);
            v_it0 = _mm512_set1_epi32(kk);
            msk1 = _mm512_cmp_epi32_mask(v_ioff,v_it0,_MM_CMPINT_LT);
            v_is = _mm512_add_epi32(v_is,v_m3);
            v_it = _mm512_add_epi32(v_it,v_m3);
            for (i = 0; i < idimp; i++) {
/*             ppart[j2+nppmx*i+npoff]    */
/*             = ppart[j1+nppmx*i+npoff]; */
               if (kk==16) {
                  v_x = _mm512_i32gather_ps(v_it,(float *)ppart,4);
                  _mm512_i32scatter_ps((float *)ppart,v_is,v_x,4);
               }
               else {
                  v_x = _mm512_mask_i32gather_ps(v_zero,msk1,v_it,
                                                 (float *)ppart,4);
                  _mm512_mask_i32scatter_ps((float *)ppart,msk1,v_is,
                                            v_x,4);
               }
               v_is = _mm512_add_epi32(v_is,v_m1);
               v_it = _mm512_add_epi32(v_it,v_m1);
            }
            ih += kk;
/* holes with locations great than npp-ip do not need to be filled */
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
               j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[l] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncpporderf3lt(float ppart[], float ppbuff[], int kpic[],
                     int ncl[], int ihole[], int idimp, int nppmx,
                     int mx1, int my1, int mz1, int npbmx, int ntmax,
                     int *irc) {
/* this subroutine sorts particles by x,y,z grid in tiles of mx, my, mz
   linear interpolation, with periodic boundary conditions
   tiles are assumed to be arranged in 3D linear memory
   the algorithm has 2 steps.  first, a prefix scan of ncl is performed
   and departing particles are buffered in ppbuff in direction order.
   then we copy the incoming particles from other tiles into ppart.
   it assumes that the number, location, and destination of particles 
   leaving a tile have been previously stored in ncl and ihole by the
   ckncgppushf3lt subroutine.
   input: all except ppbuff, irc
   output: ppart, ppbuff, kpic, ncl, irc
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppbuff[m][i][n] = i co-ordinate of particle n in tile m
   kpic[m] = number of particles in tile m
   ncl[m][i] = number of particles going to destination i, tile m
   ihole[m][:][0] = location of hole in array left by departing particle
   ihole[m][:][1] = direction destination of particle leaving hole
   all for tile m
   ihole[m][0][0] = ih, number of holes left (error, if negative)
   idimp = size of phase space = 6
   nppmx = maximum number of particles in tile
   mx1 = (system length in x direction - 1)/mx + 1
   my1 = (system length in y direction - 1)/my + 1
   mz1 = (system length in z direction - 1)/mz + 1
   npbmx = size of buffer array ppbuff
   ntmax = size of hole array for particles leaving tiles
   irc = maximum overflow, returned only if error occurs, when irc > 0
   requires KNC, ppart, ppbuff need to be 64 byte aligned
   nppmx, npbmx need to be a multiple of 16
local data                                                            */
   int mxy1, mxyz1, noff, npp, npoff, nps, nboff, ncoff;
   int i, j, k, l, ii, kx, ky, kz, ih, nh, ist, nn, mm, ll;
   int ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr;
   int ks[26];
   __m512i v_it, v_0, v_1;
   __m512i v_m1, v_m2, v_m3, v_npp, v_mm, v_is, v_it0, v_ioff;
   __m512 v_x, v_zero;
   __mmask16 msk1;
   __attribute__((aligned(64))) unsigned int ls[32], lm[32];
   mxy1 = mx1*my1;
   mxyz1 = mxy1*mz1;
   v_0 = _mm512_set1_epi32(0);
   v_1 = _mm512_set1_epi32(1);
   v_zero = _mm512_setzero_ps();
/* buffer particles that are leaving tile: update ppbuff, ncl */
/* loop over tiles */
   msk1 = _mm512_int2mask(1023);
   v_m1 = _mm512_set_epi32(11,11,11,11,11,10,9,8,3,3,3,3,3,2,1,0);
   v_m2 = _mm512_set_epi32(7,7,7,7,7,7,7,7,7,6,5,4,3,2,1,0);
#pragma omp parallel for \
private(i,j,l,npoff,nboff,noff,nps,mm,ii,ll,j1,ist,nh,ip,v_it,v_is, \
v_it0,v_ioff,ls,lm)
   for (l = 0; l < mxyz1; l++) {
      npoff = idimp*nppmx*l;
      nboff = idimp*npbmx*l;
      noff = (ntmax+1)*l;
/* find address offset for ordered ppbuff array */
/*    isum = 0;                  */
/*    for (j = 0; j < 26; j++) { */
/*       ist = ncl[j+26*l];      */
/*       ncl[j+26*l] = isum;     */
/*       isum += ist;            */
/*    }                          */
/* perform exclusive prefix scan */
/* load 26 data elements into 32 length vector with zero padding */
      mm = 26*l;
      v_it = _mm512_loadunpacklo_epi32(v_0,&ncl[mm]);
      v_it = _mm512_loadunpackhi_epi32(v_it,&ncl[mm+16]);
      _mm512_store_epi32(ls,v_it);
      v_is = _mm512_mask_loadunpacklo_epi32(v_0,msk1,&ncl[mm+16]);
      v_is = _mm512_mask_loadunpackhi_epi32(v_is,msk1,&ncl[mm+32]);
      _mm512_store_epi32(&ls[16],v_is);
      v_ioff = _mm512_setzero_epi32();
/* vector loop over elements in blocks of 16 */
      for (j = 0; j < 32; j+=16) {
/* load data */
         v_it0 = _mm512_load_epi32(&ls[j]);
/* first pass */
         v_is = _mm512_shuffle_epi32(v_it0,177);
         v_it = _mm512_mask_add_epi32(v_it0,_mm512_int2mask(43690),
                v_it0,v_is);
/* second pass */
         v_is = _mm512_shuffle_epi32(v_it,80);
         v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(52428),v_it,
                v_is);
/* third pass */
         v_is = _mm512_permutevar_epi32(v_m1,v_it);
         v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(61680),v_it,
                v_is);
/* fourth pass */
         v_is = _mm512_permutevar_epi32(v_m2,v_it);
         v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(65280),v_it,
                v_is);
/* add offset */
         v_it = _mm512_add_epi32(v_it,v_ioff);
/* next offset */
         if (j==0) {
            v_ioff = _mm512_shuffle_epi32(v_it,255);
            v_ioff = _mm512_permute4f128_epi32(v_ioff,255);
         }
/* subtract for exclusive scan */
         v_it = _mm512_sub_epi32(v_it,v_it0);
/* write data */
         _mm512_store_epi32(&ls[j],v_it);
      }
      nh = ihole[2*noff];
      nps = 16*(nh/16);
/*    nps = (nh >> 4) << 4; */
      ip = 0;
/* loop over particles leaving tile in groups of 16 */
      for (j = 0; j < nps; j+=16) {
/*       j1 = ihole[2*(j+1+(ntmax+1)*l)] - 1; */
/*       ist = ihole[1+2*(j+1+(ntmax+1)*l)];  */
         mm = 2*(j+1+noff);
         v_it = _mm512_loadunpacklo_epi32(v_0,&ihole[mm]);
         v_it = _mm512_loadunpackhi_epi32(v_it,&ihole[mm+16]);
         _mm512_store_epi32(lm,v_it);
         mm += 16;
         v_is = _mm512_loadunpacklo_epi32(v_0,&ihole[mm]);
         v_is = _mm512_loadunpackhi_epi32(v_is,&ihole[mm+16]);
         _mm512_store_epi32(&lm[16],v_is);
/* buffer particles that are leaving tile, in direction order */
         for (ll = 0; ll < 16; ll++) {
            j1 = lm[2*ll] - 1;
            ist = lm[1+2*ll];
            ii = ls[ist-1];
            if (ii < npbmx) {
               for (i = 0; i < idimp; i++) {
                  ppbuff[ii+npbmx*i+nboff]
                  = ppart[j1+nppmx*i+npoff];
               }
            }
            else {
               ip = 1;
            }
            ls[ist-1] = ii + 1;
         }
      }
/* loop over remaining particles leaving tile */
      for (j = nps; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+noff)] - 1;
         ist = ihole[1+2*(j+1+noff)];
         ii = ls[ist-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*i+nboff]
               = ppart[j1+nppmx*i+npoff];
            }
         }
         else {
            ip = 1;
         }
         ls[ist-1] = ii + 1;
      }
/* store 26 data elements into ncl */
      mm = 26*l;
      v_it = _mm512_load_epi32(ls);
      v_is = _mm512_load_epi32(&ls[16]);
     _mm512_packstorelo_epi32(&ncl[mm],v_it);
     _mm512_packstorehi_epi32(&ncl[mm+16],v_it);
     _mm512_mask_packstorelo_epi32(&ncl[mm+16],msk1,v_is);
     _mm512_mask_packstorehi_epi32(&ncl[mm+32],msk1,v_is);
/* set error */
      if (ip > 0)
         *irc = ncl[25+26*l];
   }
/* ppbuff overflow */
   if (*irc > 0)
      return;

/* copy incoming particles from buffer into ppart: update ppart, kpic */
/* loop over tiles */
   v_ioff = _mm512_set_epi32(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
   v_m1 = _mm512_set1_epi32(nppmx);
#pragma omp parallel for \
private(i,j,k,l,ii,kk,npp,nps,npoff,noff,nboff,kx,ky,kz,kl,kr,kxl,kxr, \
lk,ll,lr,ih,nh,nn,mm,ncoff,ist,j1,j2,ip,v_m2,v_m3,v_it,v_is,v_it0,v_mm, \
v_npp,v_x,msk1,ks,ls)
   for (l = 0; l < mxyz1; l++) {
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      noff = (ntmax+1)*l;
      v_m2 = _mm512_set1_epi32(noff+1);
      v_m3 = _mm512_set1_epi32(npoff);
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
      nh = ihole[2*noff];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      v_it0 = _mm512_set1_epi32(nh);
      v_is = _mm512_add_epi32(v_m2,v_it0);
      v_it0 = _mm512_sub_epi32(v_ioff,v_it0);
      v_npp = _mm512_set1_epi32(npp);
      for (ii = 0; ii < 26; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+26*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+26*ks[ii]] - ncoff;
/*       nps = 16*(ip/16); */
         nps = (ip >> 4) << 4;
/* loop over particles in this direction in groups of 16 */
         for (j = 0; j < nps; j+=16) {
/* insert incoming particles into holes */
/*          ih += 1;                        */
/*          if (ih <= nh) {                 */
/*             j1 = ihole[2*(ih+noff)] - 1; */
/*          }                               */
/* place overflow at end of array */
/*          else {                          */
/*             j1 = npp;                    */
/*             npp += 1;                    */
/*          }                               */
            v_mm = _mm512_add_epi32(_mm512_set1_epi32(ih),v_it0);
            msk1 = _mm512_cmp_epi32_mask(v_mm,v_0,_MM_CMPINT_LT);
            v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_is);
            v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_mm);
            v_mm = _mm512_mask_add_epi32(v_mm,_mm512_knot(msk1),v_mm,
                                         v_npp);
            v_it = _mm512_mask_i32gather_epi32(v_mm,msk1,v_mm,
                                               (int *)ihole,4);
            v_it = _mm512_mask_sub_epi32(v_it,msk1,v_it,v_1);
            ih += 16;
            nn = ih - nh;
            if (nn > 0) {
               nn = nn < 16 ? nn : 16;
               npp += nn;
            }
            msk1 = _mm512_cmp_epi32_mask(v_it,v_m1,_MM_CMPINT_LT);
            ll = _mm512_mask2int(_mm512_knot(msk1));
            v_it = _mm512_add_epi32(v_it,v_m3);
            for (i = 0; i < idimp; i++) {
/*             if (j1 < nppmx)                     */
/*                ppart[j1+nppmx*i+npoff]          */
/*                = ppbuff[j+ncoff+npbmx*i+nboff]; */
               mm = j + ncoff + npbmx*i + nboff;
               v_x = _mm512_loadunpacklo_ps(v_x,&ppbuff[mm]);
               v_x = _mm512_loadunpackhi_ps(v_x,&ppbuff[mm+16]);
               if (ll==0) {
                  _mm512_i32scatter_ps((float *)ppart,v_it,v_x,4);
               }
               else {
                  _mm512_mask_i32scatter_ps((float *)ppart,msk1,v_it,
                                            v_x,4);
               }
               v_it = _mm512_add_epi32(v_it,v_m1);
            }
            if (ll != 0) {
               ist = 1;
            }
         }
/* loop over remaining particles in this direction */
         for (j = nps; j < ip; j++) {
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
/* holes with locations great than npp-ip do not need to be filled */
      if (ih < nh) {
         ip = nh - ih;
         ii = nh;
         nn = ihole[2*(ii+noff)] - 1;
         v_it0 = _mm512_set1_epi32(nn);
         ih += 1;
         j2 = ihole[2*(ih+noff)] - 1;
         v_m2 = _mm512_sub_epi32(v_m2,v_1);
/* move particles from end into remaining holes */
/* holes are processed in increasing order */
/*       nps = 16*(ip/16); */
         nps = (ip >> 4) << 4;
/* loop over particles in groups of 16 */
         for (j = 0; j < nps; j+=16) {
/*          j2 = ihole[2*(ih+noff)] - 1; */
            v_mm = _mm512_add_epi32(_mm512_set1_epi32(ih),v_ioff);
            v_mm = _mm512_add_epi32(v_mm,v_m2);
            v_mm = _mm512_add_epi32(v_mm,v_mm);
            v_is = _mm512_i32gather_epi32(v_mm,(int *)ihole,4);
            v_is = _mm512_sub_epi32(v_is,v_1);
/*          j1 = npp - j - 1;                      */
/*          if (j1==nn) {                          */
/*             ii -= 1;                            */
/*             nn = ihole[2*(ii+(ntmax+1)*l)] - 1; */
/*          }                                      */
            kk = 0;
            for (ll = 0; ll < 16; ll++) {
               j1 = npp - j - ll - 1;
               if (j1==nn) {
                  ii -= 1;
                  nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
               }
               else {
                  ls[kk] = j1;
                  kk += 1;
               }
            }
            v_it = _mm512_load_epi32(ls);
            v_it0 = _mm512_set1_epi32(kk);
            msk1 = _mm512_cmp_epi32_mask(v_ioff,v_it0,_MM_CMPINT_LT);
            v_is = _mm512_add_epi32(v_is,v_m3);
            v_it = _mm512_add_epi32(v_it,v_m3);
            for (i = 0; i < idimp; i++) {
/*             ppart[j2+nppmx*i+npoff]    */
/*             = ppart[j1+nppmx*i+npoff]; */
               if (kk==16) {
                  v_x = _mm512_i32gather_ps(v_it,(float *)ppart,4);
                  _mm512_i32scatter_ps((float *)ppart,v_is,v_x,4);
               }
               else {
                  v_x = _mm512_mask_i32gather_ps(v_zero,msk1,v_it,
                                                 (float *)ppart,4);
                  _mm512_mask_i32scatter_ps((float *)ppart,msk1,v_is,
                                            v_x,4);
               }
               v_is = _mm512_add_epi32(v_is,v_m1);
               v_it = _mm512_add_epi32(v_it,v_m1);
            }
            ih += kk;
/* holes with locations great than npp-ip do not need to be filled */
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
               j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[l] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncpp2order3lt(float ppart[], float ppbuff[], int kpic[],
                     int ncl[], int ihole[], int idimp, int nppmx, 
                     int nx, int ny, int nz, int mx, int my, int mz,
                     int mx1, int my1, int mz1, int npbmx, int ntmax,
                     int *irc) {
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
   ppart[m][0][n] = position x of particle n in tile m
   ppart[m][1][n] = position y of particle n in tile m
   ppart[m][2][n] = position z of particle n in tile m
   ppbuff[m][i][n] = i co-ordinate of particle n in tile m
   kpic[m] = number of particles in tile m
   ncl[m][i] = number of particles going to destination i, tile m
   ihole[m][:][0] = location of hole in array left by departing particle
   ihole[m][:][1] = direction destination of particle leaving hole
   all for tile m
   ihole[m][0][0] = ih, number of holes left (error, if negative)
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
   requires KNC, ppart, ppbuff need to be 64 byte aligned
   nppmx, npbmx need to be a multiple of 16
local data                                                            */
   int mxy1, mxyz1, noff, moff, loff, npoff, npp, nps, nboff, ncoff;
   int i, j, k, l, ii, kx, ky, kz, ih, nh, ist, nn, mm, ll;
   int ip, j1, j2, kxl, kxr, kk, kl, kr, lk, lr;
   float anx, any, anz, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dx, dy, dz;
   int ks[26];
   __m512i v_ist, v_it, v_0, v_1, v_3, v_9;
   __m512i v_m1, v_m2, v_m3, v_m4, v_npp, v_mm, v_is, v_it0, v_ioff;
   __m512 v_anx, v_any, v_anz;
   __m512 v_dx, v_dy, v_dz, v_x;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 v_zero;
   __mmask16 msk1, msk2;
   __attribute__((aligned(64))) unsigned int ls[16], lm[32];
   mxy1 = mx1*my1;
   mxyz1 = mxy1*mz1;
   anx = (float) nx;
   any = (float) ny;
   anz = (float) nz;
   v_0 = _mm512_set1_epi32(0);
   v_1 = _mm512_set1_epi32(1);
   v_3 = _mm512_set1_epi32(3);
   v_9 = _mm512_set1_epi32(9);
   v_anx = _mm512_set1_ps(anx);
   v_any = _mm512_set1_ps(any);
   v_anz = _mm512_set1_ps(anz);
/* find and count particles leaving tiles and determine destination */
/* update ppart, ihole, ncl */
   v_zero = _mm512_setzero_ps();
/* loop over tiles */
#pragma omp parallel for \
private(i,j,k,l,ii,noff,moff,loff,npp,npoff,nps,nn,mm,ll,ih,nh,ist,dx, \
dy,dz,edgelx,edgely,edgelz,edgerx,edgery,edgerz,v_it,v_ist,v_edgelx, \
v_edgely,v_edgelz,v_edgerx,v_edgery,v_edgerz,v_dx,v_dy,v_dz,v_x,msk1, \
msk2,ls,lm)
   for (l = 0; l < mxyz1; l++) {
      loff = l/mxy1;
      k = l - mxy1*loff;
      loff = mz*loff;
      noff = k/mx1;
      moff = my*noff;
      noff = mx*(k - mx1*noff);
      npp = kpic[l];
      npoff = idimp*nppmx*l;
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
      noff = (ntmax+1)*l;
      v_edgelx = _mm512_set1_ps(edgelx);
      v_edgely = _mm512_set1_ps(edgely);
      v_edgelz = _mm512_set1_ps(edgelz);
      v_edgerx = _mm512_set1_ps(edgerx);
      v_edgery = _mm512_set1_ps(edgery);
      v_edgerz = _mm512_set1_ps(edgerz);
/* clear counters */
/*    for (j = 0; j < 26; j++) { */
/*       ncl[j+26*l] = 0;        */
/*    }                          */
      memset((void*)&ncl[26*l],0,26*sizeof(int));
      nps = 16*(npp/16);
/* loop over particles in tile in blocks of 16 */
      for (j = 0; j < nps; j+=16) {
/*       dx = ppart[j+npoff];         */
/*       dy = ppart[j+nppmx+npoff];   */
/*       dz = ppart[j+2*nppmx+npoff]; */
         v_dx = _mm512_load_ps(&ppart[j+npoff]);
         v_dy = _mm512_load_ps(&ppart[j+nppmx+npoff]);
         v_dz = _mm512_load_ps(&ppart[j+2*nppmx+npoff]);
/* find particles going out of bounds */
/*       ist = 0; */
         v_ist = _mm512_setzero_epi32();
/* count how many particles are going in each direction in ncl   */
/* save their address and destination in ihole                   */
/* use periodic boundary conditions and check for roundoff error */
/* ist = direction particle is going                             */
/*       if (dx >= edgerx) {              */
/*          if (dx >= anx)                */
/*             ppart[j+npoff] = dx - anx; */
/*          ist = 2;                      */
/*       }                                */
         msk1 = _mm512_cmp_ps_mask(v_dx,v_edgerx,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dx;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_1,v_1);
               v_ist = _mm512_mask_add_epi32(v_ist,msk1,v_ist,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dx,v_anx,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dx,v_anx);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+npoff],v_x);
            }
/*          if (dx < edgelx) {         */
/*             if (dx < 0.0) {         */
/*                dx += anx;           */
/*                if (dx < anx)        */
/*                   ist = 1;          */
/*                else                 */
/*                   dx = 0.0;         */
/*                ppart[j+npoff] = dx; */
/*             }                       */
/*             else {                  */
/*                ist = 1;             */
/*             }                       */
/*          }                          */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_1);
               msk2 = _mm512_cmp_ps_mask(v_dx,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dx,v_anx);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anx,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_ist = _mm512_add_epi32(v_ist,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+npoff],v_x);
            }
         }
/*       if (dy >= edgery) {                    */
/*          if (dy >= any)                      */
/*             ppart[j+nppmx+npoff] = dy - any; */
/*          ist += 6;                           */
/*       }                                      */
         msk1 = _mm512_cmp_ps_mask(v_dy,v_edgery,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dy;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_3,v_3);
               v_ist = _mm512_mask_add_epi32(v_ist,msk1,v_ist,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dy,v_any,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dy,v_any);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+nppmx+npoff],v_x);
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
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_3);
               msk2 = _mm512_cmp_ps_mask(v_dy,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dy,v_any);
               msk1 = _mm512_cmp_ps_mask(v_x,v_any,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_ist = _mm512_add_epi32(v_ist,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+nppmx+npoff],v_x);
            }
         }
/*       if (dz >= edgerz) {                      */
/*          if (dz >= anz)                        */
/*             ppart[j+2*nppmx+npoff] = dz - anz; */
/*          ist += 18;                            */
/*       }                                        */
         msk1 = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         msk2 = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         ii = _mm512_mask2int(_mm512_kor(msk1,msk2));
/* execute if either test result is true for any particle */
         if (ii != 0) {
            ii = _mm512_mask2int(msk1);
            v_x = v_dz;
/* write output if test result is true for any particle */
            if (ii != 0) {
               v_it = _mm512_add_epi32(v_9,v_9);
               v_ist = _mm512_mask_add_epi32(v_ist,msk1,v_ist,v_it);
               msk1 = _mm512_cmp_ps_mask(v_dz,v_anz,_MM_CMPINT_GE);
               v_x = _mm512_mask_sub_ps(v_x,msk1,v_dz,v_anz);
               ii = _mm512_mask2int(msk1);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_x);
            }
/*          if (dz < edgelz) {                 */
/*             if (dz < 0.0) {                 */
/*                dz += anz;                   */
/*                if (dz < anz)                */
/*                   ist += 9;                 */
/*                else                         */
/*                   dz = 0.0;                 */
/*                ppart[j+2*nppmx+npoff] = dz; */
/*             }                               */
/*             else {                          */
/*                ist += 9;                    */
/*             }                               */
/*          }                                  */
/* write output if test result is true for any particle */
            ii = _mm512_mask2int(msk2);
            if (ii != 0) {
               v_it = _mm512_mask_mov_epi32(v_0,msk2,v_9);
               msk2 = _mm512_cmp_ps_mask(v_dz,v_zero,_MM_CMPINT_LT);
               v_x = _mm512_mask_add_ps(v_x,msk2,v_dz,v_anz);
               msk1 = _mm512_cmp_ps_mask(v_x,v_anz,_MM_CMPINT_GE);
               msk1 = _mm512_kand(msk1,msk2);
               v_x = _mm512_mask_mov_ps(v_x,msk1,v_zero);
               v_it = _mm512_mask_mov_epi32(v_it,msk1,v_0);
               v_ist = _mm512_add_epi32(v_ist,v_it);
               ii = _mm512_mask2int(msk2);
               if (ii != 0)
                  _mm512_store_ps(&ppart[j+2*nppmx+npoff],v_x);
            }
         }
/* increment counters */
/*       if (ist > 0) {                               */
/*          ncl[ist+26*l-1] += 1;                     */
/*          ih += 1;                                  */
/*          if (ih <= ntmax) {                        */
/*             ihole[2*(ih+(ntmax+1)*l)] = j + i + 1; */
/*             ihole[1+2*(ih+(ntmax+1)*l)] = ist;     */
/*          }                                         */
/*          else {                                    */
/*             nh = 1;                                */
/*          }                                         */
/*       }                                            */
         _mm512_store_epi32(ls,v_ist);
/* remove zero ist values and left shift data */
         ll = 0;
         memset((void*)lm,0,32*sizeof(int));
         for (i = 0; i < 16; i++) {
            ist = ls[i];
            if (ist > 0) {
               lm[2*ll] = j + i + 1;
               lm[1+2*ll] = ist;
               ncl[ist+26*l-1] += 1;
               ll += 1;
            }
         }
         if (ll > 0) {
            if ((ih+ll) > ntmax) {
               nh = 1;
            }
            else {
               v_it = _mm512_load_epi32(lm);
               mm = 2*(ih+1+noff);
               _mm512_packstorelo_epi32(&ihole[mm],v_it);
               _mm512_packstorehi_epi32(&ihole[mm+16],v_it);
               if (ll > 8) {
                  v_it = _mm512_load_epi32(&lm[16]);
                  mm += 16;
                  _mm512_packstorelo_epi32(&ihole[mm],v_it);
                  _mm512_packstorehi_epi32(&ihole[mm+16],v_it);
               }
            }
            ih += ll;
         }
      }
/* loop over remaining particles in tile */
      for (j = nps; j < npp; j++) {
         dx = ppart[j+npoff];
         dy = ppart[j+nppmx+npoff];
         dz = ppart[j+2*nppmx+npoff];
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
         if (dz >= edgerz) {
            if (dz >= anz)
               ppart[j+2*nppmx+npoff] = dz - anz;
            ist += 18;
         }
         else if (dz < edgelz) {
            if (dz < 0.0) {
               dz += anz;
               if (dz < anz)
                  ist += 9;
               else
                  dz = 0.0;
               ppart[j+2*nppmx+npoff] = dz;
            }
            else {
               ist += 9;
            }
         }
         if (ist > 0) {
            ncl[ist+26*l-1] += 1;
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
   msk1 = _mm512_int2mask(1023);
   v_m1 = _mm512_set1_epi32(nppmx);
   v_ioff = _mm512_set_epi32(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
#pragma omp parallel for \
private(i,j,l,npoff,nboff,noff,nps,mm,ii,j1,ist,nh,ip,v_it,v_is,v_it0, \
v_mm,v_m2,v_m3,v_m4,lm)
   for (l = 0; l < mxyz1; l++) {
      npoff = idimp*nppmx*l;
      nboff = idimp*npbmx*l;
      noff = (ntmax+1)*l;
      v_m2 = _mm512_set_epi32(11,11,11,11,11,10,9,8,3,3,3,3,3,2,1,0);
      v_m3 = _mm512_set_epi32(7,7,7,7,7,7,7,7,7,6,5,4,3,2,1,0);
/* find address offset for ordered ppbuff array */
/*    isum = 0;                  */
/*    for (j = 0; j < 26; j++) { */
/*       ist = ncl[j+26*l];      */
/*       ncl[j+26*l] = isum;     */
/*       isum += ist;            */
/*    }                          */
/* perform exclusive prefix scan */
/* load 26 data elements into 32 length vector with zero padding */
      mm = 26*l;
      v_it = _mm512_loadunpacklo_epi32(v_0,&ncl[mm]);
      v_it = _mm512_loadunpackhi_epi32(v_it,&ncl[mm+16]);
      _mm512_store_epi32(lm,v_it);
      v_is = _mm512_mask_loadunpacklo_epi32(v_0,msk1,&ncl[mm+16]);
      v_is = _mm512_mask_loadunpackhi_epi32(v_is,msk1,&ncl[mm+32]);
      _mm512_store_epi32(&lm[16],v_is);
      v_mm = _mm512_setzero_epi32();
/* vector loop over elements in blocks of 16 */
      for (j = 0; j < 32; j+=16) {
/* load data */
         v_it0 = _mm512_load_epi32(&lm[j]);
/* first pass */
         v_is = _mm512_shuffle_epi32(v_it0,177);
         v_it = _mm512_mask_add_epi32(v_it0,_mm512_int2mask(43690),
                v_it0,v_is);
/* second pass */
         v_is = _mm512_shuffle_epi32(v_it,80);
         v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(52428),v_it,
                v_is);
/* third pass */
         v_is = _mm512_permutevar_epi32(v_m2,v_it);
         v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(61680),v_it,
                v_is);
/* fourth pass */
         v_is = _mm512_permutevar_epi32(v_m3,v_it);
         v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(65280),v_it,
                v_is);
/* add offset */
         v_it = _mm512_add_epi32(v_it,v_mm);
/* next offset */
         if (j==0) {
            v_mm = _mm512_shuffle_epi32(v_it,255);
            v_mm = _mm512_permute4f128_epi32(v_mm,255);
         }
/* subtract for exclusive scan */
         v_it = _mm512_sub_epi32(v_it,v_it0);
/* write data */
         _mm512_store_epi32(&lm[j],v_it);
      }
/* store 26 data elements into ncl */
      v_it = _mm512_load_epi32(lm);
      v_is = _mm512_load_epi32(&lm[16]);
     _mm512_packstorelo_epi32(&ncl[mm],v_it);
     _mm512_packstorehi_epi32(&ncl[mm+16],v_it);
     _mm512_mask_packstorelo_epi32(&ncl[mm+16],msk1,v_is);
     _mm512_mask_packstorehi_epi32(&ncl[mm+32],msk1,v_is);
      nh = ihole[2*noff];
      nps = 16*(nh/16);
/*    nps = (nh >> 4) << 4; */
      ip = 0;

      v_m2 = _mm512_set1_epi32(noff+1);
      v_m3 = _mm512_set1_epi32(npoff);
      v_m4 = _mm512_set1_epi32(nboff);
      v_it0 = _mm512_set1_epi32(npbmx);

/* loop over particles leaving tile in groups of 16 */
      for (j = 0; j < nps; j+=16) {
/* buffer particles that are leaving tile, in direction order */
/*       j1 = ihole[2*(j+1+(ntmax+1)*l)] - 1; */
/*       ist = ihole[1+2*(j+1+(ntmax+1)*l)];  */
         v_mm = _mm512_add_epi32(_mm512_set1_epi32(j),v_ioff);
         v_mm = _mm512_add_epi32(v_mm,v_m2);
         v_mm = _mm512_add_epi32(v_mm,v_mm);
         v_it = _mm512_i32gather_epi32(v_mm,(int *)ihole,4);
         v_it = _mm512_sub_epi32(v_it,v_1);
         v_mm = _mm512_add_epi32(v_mm,v_1);
         v_is = _mm512_i32gather_epi32(v_mm,(int *)ihole,4);
         _mm512_store_epi32(lm,v_is);
         for (ll = 0; ll < 16; ll++) {
            ist = lm[ll];
            ii = ncl[ist+26*l-1];
            if (ii < npbmx) {
               lm[ll] = ii;
            }
            else {
               ip = 1;
            }
            ncl[ist+26*l-1] = ii + 1;
         }
         v_is = _mm512_load_epi32(lm);
         v_it = _mm512_add_epi32(v_it,v_m3);
         v_is = _mm512_add_epi32(v_is,v_m4);
         if (ip==0) {
            for (i = 0; i < idimp; i++) {
/*             ppbuff[ii+npbmx*i+nboff]   */
/*             = ppart[j1+nppmx*i+npoff]; */
               v_x = _mm512_i32gather_ps(v_it,(float *)ppart,4);
               _mm512_i32scatter_ps((float *)ppbuff,v_is,v_x,4);
               v_it = _mm512_add_epi32(v_it,v_m1);
               v_is = _mm512_add_epi32(v_is,v_it0);
            }
         }

/*
         mm = 2*(j+1+noff);
         v_it = _mm512_loadunpacklo_epi32(v_0,&ihole[mm]);
         v_it = _mm512_loadunpackhi_epi32(v_it,&ihole[mm+16]);
         _mm512_store_epi32(lm,v_it);
         mm += 16;
         v_is = _mm512_loadunpacklo_epi32(v_0,&ihole[mm]);
         v_is = _mm512_loadunpackhi_epi32(v_is,&ihole[mm+16]);
         _mm512_store_epi32(&lm[16],v_is);
         for (ll = 0; ll < 16; ll++) {
            j1 = lm[2*ll] - 1;
            ist = lm[1+2*ll];
            ii = ncl[ist+26*l-1];
            if (ii < npbmx) {
               for (i = 0; i < idimp; i++) {
                  ppbuff[ii+npbmx*i+nboff]
                  = ppart[j1+nppmx*i+npoff];
               }
            }
            else {
               ip = 1;
            }
            ncl[ist+26*l-1] = ii + 1;
         }
*/
      }
/* loop over remaining particles leaving tile */
      for (j = nps; j < nh; j++) {
/* buffer particles that are leaving tile, in direction order */
         j1 = ihole[2*(j+1+noff)] - 1;
         ist = ihole[1+2*(j+1+noff)];
         ii = ncl[ist+26*l-1];
         if (ii < npbmx) {
            for (i = 0; i < idimp; i++) {
               ppbuff[ii+npbmx*i+nboff]
               = ppart[j1+nppmx*i+npoff];
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
   v_ioff = _mm512_set_epi32(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
   v_m1 = _mm512_set1_epi32(nppmx);
#pragma omp parallel for \
private(i,j,k,l,ii,kk,npp,nps,npoff,noff,nboff,kx,ky,kz,kl,kr,kxl,kxr, \
lk,ll,lr,ih,nh,nn,mm,ncoff,ist,j1,j2,ip,v_m2,v_m3,v_it,v_is,v_it0,v_mm, \
v_npp,v_x,msk1,ks,ls)
   for (l = 0; l < mxyz1; l++) {
      npp = kpic[l];
      npoff = idimp*nppmx*l;
      noff = (ntmax+1)*l;
      v_m2 = _mm512_set1_epi32(noff+1);
      v_m3 = _mm512_set1_epi32(npoff);
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
      nh = ihole[2*noff];
      ncoff = 0;
      ih = 0;
      ist = 0;
      j1 = 0;
      v_it0 = _mm512_set1_epi32(nh);
      v_is = _mm512_add_epi32(v_m2,v_it0);
      v_it0 = _mm512_sub_epi32(v_ioff,v_it0);
      v_npp = _mm512_set1_epi32(npp);
      for (ii = 0; ii < 26; ii++) {
         nboff = idimp*npbmx*ks[ii];
         if (ii > 0)
            ncoff = ncl[ii-1+26*ks[ii]];
/* ip = number of particles coming from direction ii */
         ip = ncl[ii+26*ks[ii]] - ncoff;
/*       nps = 16*(ip/16); */
         nps = (ip >> 4) << 4;
/* loop over particles in this direction in groups of 16 */
         for (j = 0; j < nps; j+=16) {
/* insert incoming particles into holes */
/*          ih += 1;                        */
/*          if (ih <= nh) {                 */
/*             j1 = ihole[2*(ih+noff)] - 1; */
/*          }                               */
/* place overflow at end of array */
/*          else {                          */
/*             j1 = npp;                    */
/*             npp += 1;                    */
/*          }                               */
            v_mm = _mm512_add_epi32(_mm512_set1_epi32(ih),v_it0);
            msk1 = _mm512_cmp_epi32_mask(v_mm,v_0,_MM_CMPINT_LT);
            v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_is);
            v_mm = _mm512_mask_add_epi32(v_mm,msk1,v_mm,v_mm);
            v_mm = _mm512_mask_add_epi32(v_mm,_mm512_knot(msk1),v_mm,
                                         v_npp);
            v_it = _mm512_mask_i32gather_epi32(v_mm,msk1,v_mm,
                                               (int *)ihole,4);
            v_it = _mm512_mask_sub_epi32(v_it,msk1,v_it,v_1);
            ih += 16;
            nn = ih - nh;
            if (nn > 0) {
               nn = nn < 16 ? nn : 16;
               npp += nn;
            }
            msk1 = _mm512_cmp_epi32_mask(v_it,v_m1,_MM_CMPINT_LT);
            ll = _mm512_mask2int(_mm512_knot(msk1));
            v_it = _mm512_add_epi32(v_it,v_m3);
            for (i = 0; i < idimp; i++) {
/*             if (j1 < nppmx)                     */
/*                ppart[j1+nppmx*i+npoff]          */
/*                = ppbuff[j+ncoff+npbmx*i+nboff]; */
               mm = j + ncoff + npbmx*i + nboff;
               v_x = _mm512_loadunpacklo_ps(v_x,&ppbuff[mm]);
               v_x = _mm512_loadunpackhi_ps(v_x,&ppbuff[mm+16]);
               if (ll==0) {
                  _mm512_i32scatter_ps((float *)ppart,v_it,v_x,4);
               }
               else {
                  _mm512_mask_i32scatter_ps((float *)ppart,msk1,v_it,
                                            v_x,4);
               }
               v_it = _mm512_add_epi32(v_it,v_m1);
            }
            if (ll != 0) {
               ist = 1;
            }
         }
/* loop over remaining particles in this direction */
         for (j = nps; j < ip; j++) {
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
/* holes with locations great than npp-ip do not need to be filled */
      if (ih < nh) {
         ip = nh - ih;
         ii = nh;
         nn = ihole[2*(ii+noff)] - 1;
         v_it0 = _mm512_set1_epi32(nn);
         ih += 1;
         j2 = ihole[2*(ih+noff)] - 1;
         v_m2 = _mm512_sub_epi32(v_m2,v_1);
/* move particles from end into remaining holes */
/* holes are processed in increasing order */
/*       nps = 16*(ip/16); */
         nps = (ip >> 4) << 4;
/* loop over particles in groups of 16 */
         for (j = 0; j < nps; j+=16) {
/*          j2 = ihole[2*(ih+noff)] - 1; */
            v_mm = _mm512_add_epi32(_mm512_set1_epi32(ih),v_ioff);
            v_mm = _mm512_add_epi32(v_mm,v_m2);
            v_mm = _mm512_add_epi32(v_mm,v_mm);
            v_is = _mm512_i32gather_epi32(v_mm,(int *)ihole,4);
            v_is = _mm512_sub_epi32(v_is,v_1);
/*          j1 = npp - j - 1;                      */
/*          if (j1==nn) {                          */
/*             ii -= 1;                            */
/*             nn = ihole[2*(ii+(ntmax+1)*l)] - 1; */
/*          }                                      */
            kk = 0;
            for (ll = 0; ll < 16; ll++) {
               j1 = npp - j - ll - 1;
               if (j1==nn) {
                  ii -= 1;
                  nn = ihole[2*(ii+(ntmax+1)*l)] - 1;
               }
               else {
                  ls[kk] = j1;
                  kk += 1;
               }
            }
            v_it = _mm512_load_epi32(ls);
            v_it0 = _mm512_set1_epi32(kk);
            msk1 = _mm512_cmp_epi32_mask(v_ioff,v_it0,_MM_CMPINT_LT);
            v_is = _mm512_add_epi32(v_is,v_m3);
            v_it = _mm512_add_epi32(v_it,v_m3);
            for (i = 0; i < idimp; i++) {
/*             ppart[j2+nppmx*i+npoff]    */
/*             = ppart[j1+nppmx*i+npoff]; */
               if (kk==16) {
                  v_x = _mm512_i32gather_ps(v_it,(float *)ppart,4);
                  _mm512_i32scatter_ps((float *)ppart,v_is,v_x,4);
               }
               else {
                  v_x = _mm512_mask_i32gather_ps(v_zero,msk1,v_it,
                                                 (float *)ppart,4);
                  _mm512_mask_i32scatter_ps((float *)ppart,msk1,v_is,
                                            v_x,4);
               }
               v_is = _mm512_add_epi32(v_is,v_m1);
               v_it = _mm512_add_epi32(v_it,v_m1);
            }
            ih += kk;
/* holes with locations great than npp-ip do not need to be filled */
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
               j2 = ihole[2*(ih+(ntmax+1)*l)] - 1;
            }
         }
         npp -= ip;
      }
      kpic[l] = npp;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cknccguard3l(float fxyz[], int nx, int ny, int nz, int nxe,
                  int nye, int nze) {
/* replicate extended periodic vector field fxyz
   linear interpolation
   nx/ny/nz = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   nze = third dimension of field arrays, must be >= nz+1
   requires KNC, fxyz needs to be 64 byte aligned
   nxe needs to be a multiple of 4
local data                                                 */
   int j, k, l, nxs, nxyen, ll;
   nxs = 4*(nx/4);
   nxyen = 4*nxe*nye;
/* copy edges of extended field */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,ll)
      for (l = 0; l < nz; l++) {
         ll = nxyen*l;
         for (k = 0; k < ny; k++) {
            fxyz[4*nx+4*nxe*k+ll] = fxyz[4*nxe*k+ll];
            fxyz[1+4*nx+4*nxe*k+ll] = fxyz[1+4*nxe*k+ll];
            fxyz[2+4*nx+4*nxe*k+ll] = fxyz[2+4*nxe*k+ll];
         }
/* vector loop over elements in blocks of 4 */
         for (j = 0; j < nxs; j+=4) {
            _mm512_mask_store_ps(&fxyz[4*j+4*nxe*ny+ll],
               _mm512_int2mask(30583),_mm512_load_ps(&fxyz[4*j+ll]));
         }
/* loop over remaining elements */
         for (j = nxs; j < nx; j++) {
            fxyz[4*j+4*nxe*ny+ll] = fxyz[4*j+ll];
            fxyz[1+4*j+4*nxe*ny+ll] = fxyz[1+4*j+ll];
            fxyz[2+4*j+4*nxe*ny+ll] = fxyz[2+4*j+ll];
         }
         fxyz[4*nx+4*nxe*ny+ll] = fxyz[ll];
         fxyz[1+4*nx+4*nxe*ny+ll] = fxyz[1+ll];
         fxyz[2+4*nx+4*nxe*ny+ll] = fxyz[2+ll];
      }
#pragma omp for \
private(j,k)
      for (k = 0; k < ny; k++) {
/* vector loop over elements in blocks of 4 */
         for (j = 0; j < nxs; j+=4) {
            _mm512_mask_store_ps(&fxyz[4*j+4*nxe*k+nxyen*nz],
            _mm512_int2mask(30583),_mm512_load_ps(&fxyz[4*j+4*nxe*k]));
         }
/* loop over remaining elements */
         for (j = nxs; j < nx; j++) {
            fxyz[4*j+4*nxe*k+nxyen*nz] = fxyz[4*j+4*nxe*k];
            fxyz[1+4*j+4*nxe*k+nxyen*nz] = fxyz[1+4*j+4*nxe*k];
            fxyz[2+4*j+4*nxe*k+nxyen*nz] = fxyz[2+4*j+4*nxe*k];
         }
         fxyz[4*nx+4*nxe*k+nxyen*nz] = fxyz[4*nxe*k];
         fxyz[1+4*nx+4*nxe*k+nxyen*nz] = fxyz[1+4*nxe*k];
         fxyz[2+4*nx+4*nxe*k+nxyen*nz] = fxyz[2+4*nxe*k];
      }
   }
/* vector loop over elements in blocks of 4 */
   for (j = 0; j < nxs; j+=4) {
      _mm512_mask_store_ps(&fxyz[4*j+4*nxe*ny+nxyen*nz],
      _mm512_int2mask(30583),_mm512_load_ps(&fxyz[4*j]));
   }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
      fxyz[4*j+4*nxe*ny+nxyen*nz] = fxyz[4*j];
      fxyz[1+4*j+4*nxe*ny+nxyen*nz] = fxyz[1+4*j];
      fxyz[2+4*j+4*nxe*ny+nxyen*nz] = fxyz[2+4*j];
   }
   fxyz[4*nx+4*nxe*ny+nxyen*nz] = fxyz[0];
   fxyz[1+4*nx+4*nxe*ny+nxyen*nz] = fxyz[1];
   fxyz[2+4*nx+4*nxe*ny+nxyen*nz] = fxyz[2];
   return;
}

/*--------------------------------------------------------------------*/
void ckncacguard3l(float cu[], int nx, int ny, int nz, int nxe, int nye,
                   int nze) {
/* accumulate extended periodic field cu
   linear interpolation
   nx/ny/nz = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   nze = third dimension of field arrays, must be >= nz+1
   requires KNC, fxyz needs to be 64 byte aligned
   nxe needs to be a multiple of 4
local data                                                 */
   int j, k, l, nxs, nxyen, ll;
   __m512 v_cu, v_zero;
   nxs = 4*(nx/4);
   nxyen = 4*nxe*nye;
   v_zero = _mm512_set1_ps(0.0f);
/* accumulate edges of extended field */
#pragma omp parallel
   {
#pragma omp for \
private(j,k,l,ll,v_cu)
      for (l = 0; l < nz; l++) {
         ll = nxyen*l;
         for (k = 0; k < ny; k++) {
            cu[4*nxe*k+ll] += cu[4*nx+4*nxe*k+ll];
            cu[1+4*nxe*k+ll] += cu[1+4*nx+4*nxe*k+ll];
            cu[2+4*nxe*k+ll] += cu[2+4*nx+4*nxe*k+ll];
            cu[4*nx+4*nxe*k+ll] = 0.0;
            cu[1+4*nx+4*nxe*k+ll] = 0.0;
            cu[2+4*nx+4*nxe*k+ll] = 0.0;
         }
/* vector loop over elements in blocks of 4 */
         for (j = 0; j < nxs; j+=4) {
            v_cu = _mm512_load_ps(&cu[4*j+4*nxe*ny+ll]);
            v_cu = _mm512_add_ps(_mm512_load_ps(&cu[4*j+ll]),v_cu);
            _mm512_store_ps(&cu[4*j+ll],v_cu);
            _mm512_store_ps(&cu[4*j+4*nxe*ny+ll],v_zero);
         }
/* loop over remaining elements */
         for (j = nxs; j < nx; j++) {
            cu[4*j+ll] += cu[4*j+4*nxe*ny+ll];
            cu[1+4*j+ll] += cu[1+4*j+4*nxe*ny+ll];
            cu[2+4*j+ll] += cu[2+4*j+4*nxe*ny+ll];
            cu[4*j+4*nxe*ny+ll] = 0.0;
            cu[1+4*j+4*nxe*ny+ll] = 0.0;
            cu[2+4*j+4*nxe*ny+ll] = 0.0;
         }
         cu[ll] += cu[4*nx+4*nxe*ny+ll];
         cu[1+ll] += cu[1+4*nx+4*nxe*ny+ll];
         cu[2+ll] += cu[2+4*nx+4*nxe*ny+ll];
         cu[4*nx+4*nxe*ny+ll] = 0.0;
         cu[1+4*nx+4*nxe*ny+ll] = 0.0;
         cu[2+4*nx+4*nxe*ny+ll] = 0.0;
      }
#pragma omp for \
private(j,k,v_cu)
      for (k = 0; k < ny; k++) {
/* vector loop over elements in blocks of 4 */
         for (j = 0; j < nxs; j+=4) {
            v_cu = _mm512_load_ps(&cu[4*j+4*nxe*k+nxyen*nz]);
            v_cu = _mm512_add_ps(_mm512_load_ps(&cu[4*j+4*nxe*k]),v_cu);
            _mm512_store_ps(&cu[4*j+4*nxe*k],v_cu);
            _mm512_store_ps(&cu[4*j+4*nxe*k+nxyen*nz],v_zero);
         }
/* loop over remaining elements */
         for (j = nxs; j < nx; j++) {
            cu[4*j+4*nxe*k] += cu[4*j+4*nxe*k+nxyen*nz];
            cu[1+4*j+4*nxe*k] += cu[1+4*j+4*nxe*k+nxyen*nz];
            cu[2+4*j+4*nxe*k] += cu[2+4*j+4*nxe*k+nxyen*nz];
            cu[4*j+4*nxe*k+nxyen*nz] = 0.0;
            cu[1+4*j+4*nxe*k+nxyen*nz] = 0.0;
            cu[2+4*j+4*nxe*k+nxyen*nz] = 0.0;
         }
         cu[4*nxe*k] += cu[4*nx+4*nxe*k+nxyen*nz];
         cu[1+4*nxe*k] += cu[1+4*nx+4*nxe*k+nxyen*nz];
         cu[2+4*nxe*k] += cu[2+4*nx+4*nxe*k+nxyen*nz];
         cu[4*nx+4*nxe*k+nxyen*nz] = 0.0;
         cu[1+4*nx+4*nxe*k+nxyen*nz] = 0.0;
         cu[2+4*nx+4*nxe*k+nxyen*nz] = 0.0;
      }
   }
/* vector loop over elements in blocks of 4 */
   for (j = 0; j < nxs; j+=4) {
      v_cu = _mm512_load_ps(&cu[4*j+4*nxe*ny+nxyen*nz]);
      v_cu = _mm512_add_ps(_mm512_load_ps(&cu[4*j]),v_cu);
      _mm512_store_ps(&cu[4*j],v_cu);
      _mm512_store_ps(&cu[4*j+4*nxe*ny+nxyen*nz],v_zero);
   }
/* loop over remaining elements */
   for (j = nxs; j < nx; j++) {
      cu[4*j] += cu[4*j+4*nxe*ny+nxyen*nz];
      cu[1+4*j] += cu[1+4*j+4*nxe*ny+nxyen*nz];
      cu[2+4*j] += cu[2+4*j+4*nxe*ny+nxyen*nz];
      cu[4*j+4*nxe*ny+nxyen*nz] = 0.0;
      cu[1+4*j+4*nxe*ny+nxyen*nz] = 0.0;
      cu[2+4*j+4*nxe*ny+nxyen*nz] = 0.0;
   }
   cu[0] += cu[4*nx+4*nxe*ny+nxyen*nz];
   cu[1] += cu[1+4*nx+4*nxe*ny+nxyen*nz];
   cu[2] += cu[2+4*nx+4*nxe*ny+nxyen*nz];
   cu[4*nx+4*nxe*ny+nxyen*nz] = 0.0;
   cu[1+4*nx+4*nxe*ny+nxyen*nz] = 0.0;
   cu[2+4*nx+4*nxe*ny+nxyen*nz] = 0.0;
   return;
}

/*--------------------------------------------------------------------*/
void ckncaguard3l(float q[], int nx, int ny, int nz, int nxe, int nye,
                  int nze) {
/* accumulate extended periodic scalar field q
   linear interpolation
   nx/ny/nz = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   nze = third dimension of field arrays, must be >= nz+1
   requires KNC, q needs to be 64 byte aligned
   nxe needs to be a multiple of 16
local data                                                 */
   int j, k, l, nxs, nxye, ll;
   __m512 v_q;
   nxs = 16*(nx/16);
   nxye = nxe*nye;
/* accumulate edges of extended field */
#pragma omp parallel
   {
#pragma omp for \
private(j,k,l,ll,v_q)
      for (l = 0; l < nz; l++) {
         ll = nxye*l;
         for (k = 0; k < ny; k++) {
            q[nxe*k+ll] += q[nx+nxe*k+ll];
            q[nx+nxe*k+ll] = 0.0;
         }
/* vector loop over elements in blocks of 16 */
         for (j = 0; j < nxs; j+=16) {
            v_q = _mm512_load_ps(&q[j+nxe*ny+ll]);
            v_q = _mm512_add_ps(_mm512_load_ps(&q[j+ll]),v_q);
            _mm512_store_ps(&q[j+ll],v_q);
            _mm512_store_ps(&q[j+nxe*ny+ll],_mm512_setzero_ps());
         }
/* loop over remaining elements */
         for (j = nxs; j < nx; j++) {
            q[j+ll] += q[j+nxe*ny+ll];
            q[j+nxe*ny+ll] = 0.0;
         }
         q[ll] += q[nx+nxe*ny+ll];
         q[nx+nxe*ny+ll] = 0.0;
      }
#pragma omp for \
private(j,k,v_q)
      for (k = 0; k < ny; k++) {
/* vector loop over elements in blocks of 16 */
         for (j = 0; j < nxs; j+=16) {
            v_q = _mm512_load_ps(&q[j+nxe*k+nxye*nz]);
            v_q = _mm512_add_ps(_mm512_load_ps(&q[j+nxe*k]),v_q);
            _mm512_store_ps(&q[j+nxe*k],v_q);
            _mm512_store_ps(&q[j+nxe*k+nxye*nz],_mm512_setzero_ps());
         }
/* loop over remaining elements */
         for (j = nxs; j < nx; j++) {
            q[j+nxe*k] += q[j+nxe*k+nxye*nz];
            q[j+nxe*k+nxye*nz] = 0.0;
         }
         q[nxe*k] += q[nx+nxe*k+nxye*nz];
         q[nx+nxe*k+nxye*nz] = 0.0;
      }
   }
/* vector loop over elements in blocks of 16 */
      for (j = 0; j < nxs; j+=16) {
         v_q = _mm512_load_ps(&q[j+nxe*ny+nxye*nz]);
         v_q = _mm512_add_ps(_mm512_load_ps(&q[j]),v_q);
         _mm512_store_ps(&q[j],v_q);
         _mm512_store_ps(&q[j+nxe*ny+nxye*nz],_mm512_setzero_ps());
   }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
      q[j] += q[j+nxe*ny+nxye*nz];
      q[j+nxe*ny+nxye*nz] = 0.0;
   }
   q[0] += q[nx+nxe*ny+nxye*nz];
   q[nx+nxe*ny+nxye*nz] = 0.0;
   return;
}

/*--------------------------------------------------------------------*/
void ckncmpois33(float complex q[], float complex fxyz[], int isign,
                 float complex ffc[], float ax, float ay, float az,
                 float affp, float *we, int nx, int ny, int nz,
                 int nxvh, int nyv, int nzv, int nxhd, int nyhd,
                 int nzhd) {
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
   requires KNC, q, fxy, ffc need to be 64 byte aligned
   nxhd, nxvh need to be a multiple of 8
   fxyz needs to have 4 components
local data                                                 */
   int nxh, nyh, nzh, nxhs, itn, j, k, l, k1, l1, kk, kj, ll, lj;
   int nxyhd, nxvyh; 
   float dnx, dny, dnz, dkx, dky, dkz, at1, at2, at3, at4, at5, at6;
   float complex zero, zt1, zt2;
   double wp, sum1, sum2;
   __m512i v_j, v_it, v_perm;
   __m512 v_dnx, v_dny, v_dnz, v_dky, v_dkz, v_at1, v_at2, v_at3, v_at4;
   __m512 v_zero, v_zt1, v_zt2, v_zt3, v_zt4;
   __m512 a, b, c, d, e, f, g, h;
   __m512d v_wp, v_d;
   __attribute__((aligned(64))) double dd[8];
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxhs = 8*(nxh/8);
   itn = 1 > nxhs ? 1 : nxhs;
   nxyhd = nxhd*nyhd;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   zero = 0.0 + 0.0*_Complex_I;
   v_j = _mm512_set_epi32(7,7,6,6,5,5,4,4,3,3,2,2,1,1,0,0);
   v_dnx = _mm512_set1_ps(dnx);
   v_dny = _mm512_set1_ps(dny);
   v_dnz = _mm512_set1_ps(dnz);
   v_zero = _mm512_setzero_ps();
   v_perm = _mm512_set_epi32(15,14,11,10,7,6,3,2,13,12,9,8,5,4,1,0);
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
private(j,k,l,k1,l1,ll,lj,kk,kj,dky,dkz,at1,at2,at3,at4,zt1,zt2,wp, \
v_it,v_dky,v_dkz,v_at1,v_at2,v_at3,v_at4,v_zt1,v_zt2,v_zt3,v_zt4,a,b, \
c,d,e,f,g,h,v_d,v_wp,dd) \
reduction(+:sum1)
      for (l = 1; l < nzh; l++) {
         dkz = dnz*(float) l;
         v_dkz = _mm512_cvtfxpnt_round_adjustepi32_ps(
               _mm512_set1_epi32(l),_MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkz = _mm512_mul_ps(v_dnz,v_dkz);
         ll = nxyhd*l;
         lj = nxvyh*l;
         l1 = nxvyh*nz - lj;
         wp = 0.0;
         v_wp = _mm512_setzero_pd();
         for (k = 1; k < nyh; k++) {
            dky = dny*(float) k;
            v_it = _mm512_set1_epi32(k);
            v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_dky = _mm512_mul_ps(v_dny,v_dky);
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
/* vector loop over elements in blocks of 8 */
            for (j = 0; j < nxhs; j+=8) {
/*             at1 = crealf(ffc[j+kk+ll])*cimagf(ffc[j+kk+ll]); */
               v_at1 = _mm512_load_ps((float *)&ffc[j+kk+ll]);
               v_at2 = (__m512)_mm512_shuffle_epi32((__m512i)v_at1,177);
               v_at1 = _mm512_mul_ps(v_at1,v_at2);
/*             at2 = at1*dnx*(float) j; */
               v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
               v_at2 = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                       _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
               v_at2 = _mm512_mul_ps(v_at1,_mm512_mul_ps(v_dnx,v_at2));
/*             at3 = dky*at1; */
               v_at3 = _mm512_mul_ps(v_dky,v_at1);
/*             at4 = dkz*at1; */
               v_at4 = _mm512_mul_ps(v_dkz,v_at1);
/*             zt1 = cimagf(q[j+kj+lj]) - crealf(q[j+kj+lj])*_Complex_I; */
               v_zt1 = _mm512_load_ps((float *)&q[j+kj+lj]);
               v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),
                       v_zero,v_zt1);
               v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/*             zt2 = cimagf(q[j+k1+lj]) - crealf(q[j+k1+lj])*_Complex_I; */
               v_zt2 = _mm512_load_ps((float *)&q[j+k1+lj]);
               v_zt2 = _mm512_mask_sub_ps(v_zt2,_mm512_int2mask(21845),
                       v_zero,v_zt2);
               v_zt2 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt2,177);
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),
                                             v_zero);
                  v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(3),
                                             v_zero);
               }
/*             fxyz[4*(j+kj+lj)] = at2*zt1;   */
/*             fxyz[1+4*(j+kj+lj)] = at3*zt1; */
/*             fxyz[2+4*(j+kj+lj)] = at4*zt1; */
               a = _mm512_mul_ps(v_at2,v_zt1);
               b = _mm512_mul_ps(v_at3,v_zt1);
               c = _mm512_mul_ps(v_at4,v_zt1);
/* perform 4x16 transpose for fxyz field components */
               e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),
                    c,78);
               f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),
                   a,78);
               g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                    v_zero,78);
               h = _mm512_mask_permute4f128_ps(v_zero,
                   _mm512_int2mask(255),b,78);
               a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),
                                               g,177);
               b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),
                                               e,177);
               c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),
                                               h,177);
               d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),
                                               f,177);
               a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
               b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
               c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
               d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
               _mm512_store_ps((float *)&fxyz[4*(j+kj+lj)],a);
               _mm512_store_ps((float *)&fxyz[8+4*(j+kj+lj)],b);
               _mm512_store_ps((float *)&fxyz[16+4*(j+kj+lj)],c);
               _mm512_store_ps((float *)&fxyz[24+4*(j+kj+lj)],d);
/*             fxyz[4*(j+k1+lj)] = at2*zt2;    */
/*             fxyz[1+4*(j+k1+lj)] = -at3*zt2; */
/*             fxyz[2+4*(j+k1+lj)] = at4*zt2;  */
               a = _mm512_mul_ps(v_at2,v_zt2);
               b = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at3,v_zt2));
               c = _mm512_mul_ps(v_at4,v_zt2);
/* perform 4x16 transpose for fxyz field components */
               e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),
                    c,78);
               f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),
                   a,78);
               g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                    v_zero,78);
               h = _mm512_mask_permute4f128_ps(v_zero,
                   _mm512_int2mask(255),b,78);
               a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),
                                               g,177);
               b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),
                                               e,177);
               c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),
                                               h,177);
               d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),
                                               f,177);
               a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
               b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
               c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
               d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
               _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],a);
               _mm512_store_ps((float *)&fxyz[8+4*(j+k1+lj)],b);
               _mm512_store_ps((float *)&fxyz[16+4*(j+k1+lj)],c);
               _mm512_store_ps((float *)&fxyz[24+4*(j+k1+lj)],d);
/*             wp += at1*(q[j+kj+lj]*conjf(q[j+kj+lj]) */
/*                + q[j+k1+lj]*conjf(q[j+k1+lj]));     */
               v_zt3 = _mm512_mul_ps(v_zt1,v_zt1);
               v_zt3 = _mm512_add_ps(v_zt3,_mm512_mul_ps(v_zt2,v_zt2));
               v_zt3 = _mm512_mul_ps(v_at1,v_zt3);
/*             zt1 = cimagf(q[j+kj+l1]) - crealf(q[j+kj+l1])*_Complex_I; */
               v_zt1 = _mm512_load_ps((float *)&q[j+kj+l1]);
               v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),
                       v_zero,v_zt1);
               v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/*             zt2 = cimagf(q[j+k1+l1]) - crealf(q[j+k1+l1])*_Complex_I; */
               v_zt2 = _mm512_load_ps((float *)&q[j+k1+l1]);
               v_zt2 = _mm512_mask_sub_ps(v_zt2,_mm512_int2mask(21845),
                       v_zero,v_zt2);
               v_zt2 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt2,177);
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),
                                             v_zero);
                  v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(3),
                                             v_zero);
               }
/*             fxyz[4*(j+kj+l1)] = at2*zt1;    */
/*             fxyz[1+4*(j+kj+l1)] = at3*zt1;  */
/*             fxyz[2+4*(j+kj+l1)] = -at4*zt1; */
               a = _mm512_mul_ps(v_at2,v_zt1);
               b = _mm512_mul_ps(v_at3,v_zt1);
               c = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at4,v_zt1));
/* perform 4x16 transpose for fxyz field components */
               e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),
                                               c,78);
               f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),
                                               a,78);
               g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                                               v_zero,78);
               h = _mm512_mask_permute4f128_ps(v_zero,
                   _mm512_int2mask(255),b,78);
               a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),
                                               g,177);
               b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),
                                               e,177);
               c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),
                                               h,177);
               d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),
                                               f,177);
               a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
               b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
               c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
               d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
               _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],a);
               _mm512_store_ps((float *)&fxyz[8+4*(j+kj+l1)],b);
               _mm512_store_ps((float *)&fxyz[16+4*(j+kj+l1)],c);
               _mm512_store_ps((float *)&fxyz[24+4*(j+kj+l1)],d);
/*             fxyz[4*(j+k1+l1)] = at2*zt2;    */
/*             fxyz[1+4*(j+k1+l1)] = -at3*zt2; */
/*             fxyz[2+4*(j+k1+l1)] = -at4*zt2; */
               a = _mm512_mul_ps(v_at2,v_zt2);
               b = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at3,v_zt2));
               c = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at4,v_zt2));
/* perform 4x16 transpose for fxyz field components */
               e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),
                    c,78);
               f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),
                   a,78);
               g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                    v_zero,78);
               h = _mm512_mask_permute4f128_ps(v_zero,
                   _mm512_int2mask(255),b,78);
               a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),
                                               g,177);
               b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),
                                               e,177);
               c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),
                                               h,177);
               d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),
                                               f,177);
               a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
               b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
               c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
               d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
               _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],a);
               _mm512_store_ps((float *)&fxyz[8+4*(j+k1+l1)],b);
               _mm512_store_ps((float *)&fxyz[16+4*(j+k1+l1)],c);
               _mm512_store_ps((float *)&fxyz[24+4*(j+k1+l1)],d);
/*             wp += at1*(q[j+kj+l1]*conjf(q[j+kj+l1]) */
/*                + q[j+k1+l1]*conjf(q[j+k1+l1]));     */
               v_zt4 = _mm512_mul_ps(v_zt1,v_zt1);
               v_zt4 = _mm512_add_ps(v_zt4,_mm512_mul_ps(v_zt2,v_zt2));
               v_zt3 = _mm512_add_ps(v_zt3,_mm512_mul_ps(v_at1,v_zt4));
/* convert to double precision before accumulating */
               v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt3));
               v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt3,78));
               v_wp = _mm512_add_pd(v_wp,v_d);
            }
/* loop over remaining elements */
            for (j = itn; j < nxh; j++) {
               at1 = crealf(ffc[j+kk+ll])*cimagf(ffc[j+kk+ll]);
               at2 = at1*dnx*(float) j;
               at3 = dky*at1;
               at4 = dkz*at1;
               zt1 = cimagf(q[j+kj+lj]) - crealf(q[j+kj+lj])*_Complex_I;
               zt2 = cimagf(q[j+k1+lj]) - crealf(q[j+k1+lj])*_Complex_I;
               fxyz[4*(j+kj+lj)] = at2*zt1;
               fxyz[1+4*(j+kj+lj)] = at3*zt1;
               fxyz[2+4*(j+kj+lj)] = at4*zt1;
               fxyz[4*(j+k1+lj)] = at2*zt2;
               fxyz[1+4*(j+k1+lj)] = -at3*zt2;
               fxyz[2+4*(j+k1+lj)] = at4*zt2;
               zt1 = cimagf(q[j+kj+l1]) - crealf(q[j+kj+l1])*_Complex_I;
               zt2 = cimagf(q[j+k1+l1]) - crealf(q[j+k1+l1])*_Complex_I;
               fxyz[4*(j+kj+l1)] = at2*zt1;
               fxyz[1+4*(j+kj+l1)] = at3*zt1;
               fxyz[2+4*(j+kj+l1)] = -at4*zt1;
               fxyz[4*(j+k1+l1)] = at2*zt2;
               fxyz[1+4*(j+k1+l1)] = -at3*zt2;
               fxyz[2+4*(j+k1+l1)] = -at4*zt2;
               at1 = at1*(q[j+kj+lj]*conjf(q[j+kj+lj])                 
                   + q[j+k1+lj]*conjf(q[j+k1+lj])
                   + q[j+kj+l1]*conjf(q[j+kj+l1])      
                   + q[j+k1+l1]*conjf(q[j+k1+l1]));
               wp += (double) at1;
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
            fxyz[4*(kj+lj)] = zero;
            fxyz[1+4*(kj+lj)] = at3*zt1;
            fxyz[2+4*(kj+lj)] = at4*zt1;
            fxyz[4*(k1+lj)] = zero;
            fxyz[1+4*(k1+lj)] = zero;
            fxyz[2+4*(k1+lj)] = zero;
            fxyz[4*(kj+l1)] = zero;
            fxyz[1+4*(kj+l1)] = at3*zt2;
            fxyz[2+4*(kj+l1)] = -at4*zt2;
            fxyz[4*(k1+l1)] = zero;
            fxyz[1+4*(k1+l1)] = zero;
            fxyz[2+4*(k1+l1)] = zero;
            at1 = at1*(q[kj+lj]*conjf(q[kj+lj])
                + q[kj+l1]*conjf(q[kj+l1]));
            wp += (double) at1;
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
/* vector loop over elements in blocks of 8 */
         for (j = 0; j < nxhs; j+=8) {
/*          at1 = crealf(ffc[j+ll])*cimagf(ffc[j+ll]); */
            v_at1 = _mm512_load_ps((float *)&ffc[j+ll]);
            v_at2 = (__m512)_mm512_shuffle_epi32((__m512i)v_at1,177);
            v_at1 = _mm512_mul_ps(v_at1,v_at2);
/*          at2 = at1*dnx*(float) j; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_at2 = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_at2 = _mm512_mul_ps(v_at1,_mm512_mul_ps(v_dnx,v_at2));
/*          at4 = dkz*at1; */
            v_at4 = _mm512_mul_ps(v_dkz,v_at1);
/*          zt1 = cimagf(q[j+lj]) - crealf(q[j+lj])*_Complex_I; */
            v_zt1 = _mm512_load_ps((float *)&q[j+lj]);
            v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),
                    v_zero,v_zt1);
            v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/*          zt2 = cimagf(q[j+l1]) - crealf(q[j+l1])*_Complex_I; */
            v_zt2 = _mm512_load_ps((float *)&q[j+l1]);
            v_zt2 = _mm512_mask_sub_ps(v_zt2,_mm512_int2mask(21845),
                    v_zero,v_zt2);
            v_zt2 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt2,177);
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),
                                          v_zero);
               v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(3),
                                          v_zero);
            }
/*          fxyz[4*(j+lj)] = at2*zt1;   */
/*          fxyz[1+4*(j+lj)] = zero;    */
/*          fxyz[2+4*(j+lj)] = at4*zt1; */
            a = _mm512_mul_ps(v_at2,v_zt1);
            b = v_zero;
            c = _mm512_mul_ps(v_at4,v_zt1);
/* perform 4x16 transpose for fxyz field components */
            e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),c,
                                            78);
            f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),a,
                                            78);
            g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                                            v_zero,78);
            h = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(255),
                                            b,78);
            a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),g,
                                            177);
            b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),e,
                                            177);
            c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),h,
                                            177);
            d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),f,
                                            177);
            a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
            b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
            c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
            d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
            _mm512_store_ps((float *)&fxyz[4*(j+lj)],a);
            _mm512_store_ps((float *)&fxyz[8+4*(j+lj)],b);
            _mm512_store_ps((float *)&fxyz[16+4*(j+lj)],c);
            _mm512_store_ps((float *)&fxyz[24+4*(j+lj)],d);
/*          fxyz[4*(j+k1+lj)] = zero;   */
/*          fxyz[1+4*(j+k1+lj)] = zero; */
/*          fxyz[2+4*(j+k1+lj)] = zero; */
            _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],v_zero);
            _mm512_store_ps((float *)&fxyz[8+4*(j+k1+lj)],v_zero);
            _mm512_store_ps((float *)&fxyz[16+4*(j+k1+lj)],v_zero);
            _mm512_store_ps((float *)&fxyz[24+4*(j+k1+lj)],v_zero);
/*          fxyz[4*(j+l1)] = at2*zt2;    */
/*          fxyz[1+4*(j+l1)] = zero;     */
/*          fxyz[2+4*(j+l1)] = -at4*zt2; */
            a = _mm512_mul_ps(v_at2,v_zt2);
            b = v_zero;
            c = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at4,v_zt2));
/* perform 4x16 transpose for fxyz field components */
            e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),c,
                                            78);
            f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),a,
                                            78);
            g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                                            v_zero,78);
            h = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(255),
                                            b,78);
            a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),g,
                                            177);
            b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),e,
                                            177);
            c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),h,
                                            177);
            d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),f,
                                            177);
            a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
            b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
            c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
            d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
            _mm512_store_ps((float *)&fxyz[4*(j+l1)],a);
            _mm512_store_ps((float *)&fxyz[8+4*(j+l1)],b);
            _mm512_store_ps((float *)&fxyz[16+4*(j+l1)],c);
            _mm512_store_ps((float *)&fxyz[24+4*(j+l1)],d);
/*          fxyz[4*(j+k1+l1)] = zero;   */
/*          fxyz[1+4*(j+k1+l1)] = zero; */
/*          fxyz[2+4*(j+k1+l1)] = zero; */
            _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zero);
            _mm512_store_ps((float *)&fxyz[8+4*(j+k1+l1)],v_zero);
            _mm512_store_ps((float *)&fxyz[16+4*(j+k1+l1)],v_zero);
            _mm512_store_ps((float *)&fxyz[24+4*(j+k1+l1)],v_zero);
/*          wp += at1*(q[j+lj]*conjf(q[j+lj]) */
/*             + q[j+l1]*conjf(q[j+l1]));     */
            v_zt3 = _mm512_mul_ps(v_zt1,v_zt1);
            v_zt3 = _mm512_add_ps(v_zt3,_mm512_mul_ps(v_zt2,v_zt2));
            v_zt3 = _mm512_mul_ps(v_at1,v_zt3);
/* convert to double precision before accumulating */
            v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt3));
            v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt3,78));
            v_wp = _mm512_add_pd(v_wp,v_d);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            at1 = crealf(ffc[j+ll])*cimagf(ffc[j+ll]);
            at2 = at1*dnx*(float) j;  
            at4 = dkz*at1;
            zt1 = cimagf(q[j+lj]) - crealf(q[j+lj])*_Complex_I;
            zt2 = cimagf(q[j+l1]) - crealf(q[j+l1])*_Complex_I;
            fxyz[4*(j+lj)] = at2*zt1;
            fxyz[1+4*(j+lj)] = zero;
            fxyz[2+4*(j+lj)] = at4*zt1;
            fxyz[4*(j+k1+lj)] = zero;
            fxyz[1+4*(j+k1+lj)] = zero;
            fxyz[2+4*(j+k1+lj)] = zero;
            fxyz[4*(j+l1)] = at2*zt2;
            fxyz[1+4*(j+l1)] = zero;
            fxyz[2+4*(j+l1)] = -at4*zt2;
            fxyz[4*(j+k1+l1)] = zero;
            fxyz[1+4*(j+k1+l1)] = zero;
            fxyz[2+4*(j+k1+l1)] = zero;
            at1 = at1*(q[j+lj]*conjf(q[j+lj])                           
                + q[j+l1]*conjf(q[j+l1]));
            wp += (double) at1;
         }
/* mode numbers kx = 0, nx/2 */
         at1 = crealf(ffc[ll])*cimagf(ffc[ll]);
         at4 = dkz*at1;
         zt1 = cimagf(q[lj]) - crealf(q[lj])*_Complex_I;
         fxyz[4*lj] = zero;
         fxyz[1+4*lj] = zero;
         fxyz[2+4*lj] = at4*zt1;
         fxyz[4*(k1+lj)] = zero;
         fxyz[1+4*(k1+lj)] = zero;
         fxyz[2+4*(k1+lj)] = zero;
         fxyz[4*l1] = zero;
         fxyz[1+4*l1] = zero;
         fxyz[2+4*l1] = zero;
         fxyz[4*(k1+l1)] = zero;
         fxyz[1+4*(k1+l1)] = zero;
         fxyz[2+4*(k1+l1)] = zero;
         at1 = at1*(q[lj]*conjf(q[lj]));
         wp += (double) at1;
/*       sum1 += wp; */
         _mm512_store_pd(&dd[0],v_wp);
         for (j = 1; j < 8; j++) {
            dd[0] += dd[j];
         }
         sum1 += (wp + dd[0]);
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   sum2 = 0.0;
#pragma omp parallel for \
private(j,k,k1,kk,kj,dky,at1,at2,at3,zt1,zt2,wp,v_it,v_dky,v_at1, \
v_at2,v_at3,v_zt1,v_zt2,v_zt3,a,b,c,d,e,f,g,h,v_d,v_wp,dd) \
reduction(+:sum2)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      v_it = _mm512_set1_epi32(k);
      v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dky = _mm512_mul_ps(v_dny,v_dky);
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      wp = 0.0;
      v_wp = _mm512_set1_pd(0.0);
/* vector loop over elements in blocks of 8 */
      for (j = 0; j < nxhs; j+=8) {
/*       at1 = crealf(ffc[j+kk])*cimagf(ffc[j+kk]); */
         v_at1 = _mm512_load_ps((float *)&ffc[j+kk]);
         v_at2 = (__m512)_mm512_shuffle_epi32((__m512i)v_at1,177);
         v_at1 = _mm512_mul_ps(v_at1,v_at2);
/*       at2 = at1*dnx*(float) j; */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_at2 = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_at2 = _mm512_mul_ps(v_at1,_mm512_mul_ps(v_dnx,v_at2));
/*       at3 = dky*at1; */
         v_at3 = _mm512_mul_ps(v_dky,v_at1);
/*       zt1 = cimagf(q[j+kj]) - crealf(q[j+kj])*_Complex_I; */
         v_zt1 = _mm512_load_ps((float *)&q[j+kj]);
         v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),
                 v_zero,v_zt1);
         v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/*       zt2 = cimagf(q[j+k1]) - crealf(q[j+k1])*_Complex_I; */
         v_zt2 = _mm512_load_ps((float *)&q[j+k1]);
         v_zt2 = _mm512_mask_sub_ps(v_zt2,_mm512_int2mask(21845),
                 v_zero,v_zt2);
         v_zt2 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt2,177);
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),
                                       v_zero);
            v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(3),
                                       v_zero);
         }
/*       fxyz[4*(j+kj)] = at2*zt1;   */
/*       fxyz[1+4*(j+kj)] = at3*zt1; */
/*       fxyz[2+4*(j+kj)] = zero;    */
         a = _mm512_mul_ps(v_at2,v_zt1);
         b = _mm512_mul_ps(v_at3,v_zt1);
         c = v_zero;
 /* perform 4x16 transpose for fxyz field components */
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),c,78);
         f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),a,78);
         g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                                         v_zero,78);
         h = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(255),b,
                                         78);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),g,
                                         177);
         b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),e,177);
         c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),h,
                                         177);
         d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),f,177);
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         _mm512_store_ps((float *)&fxyz[4*(j+kj)],a);
         _mm512_store_ps((float *)&fxyz[8+4*(j+kj)],b);
         _mm512_store_ps((float *)&fxyz[16+4*(j+kj)],c);
         _mm512_store_ps((float *)&fxyz[24+4*(j+kj)],d);
/*       fxyz[4*(j+k1)] = at2*zt2;    */
/*       fxyz[1+4*(j+k1)] = -at3*zt2; */
/*       fxyz[2+4*(j+k1)] = zero;     */
         a = _mm512_mul_ps(v_at2,v_zt2);
         b = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at3,v_zt2));
         c = v_zero;;
/* perform 4x16 transpose for fxyz field components */
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),c,78);
         f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),a,78);
         g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                                         v_zero,78);
         h = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(255),b,
                                         78);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),g,
                                         177);
         b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),e,177);
         c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),h,
                                         177);
         d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),f,177);
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         _mm512_store_ps((float *)&fxyz[4*(j+k1)],a);
         _mm512_store_ps((float *)&fxyz[8+4*(j+k1)],b);
         _mm512_store_ps((float *)&fxyz[16+4*(j+k1)],c);
         _mm512_store_ps((float *)&fxyz[24+4*(j+k1)],d);
/*       fxyz[4*(j+kj+l1)] = zero;    */
/*       fxyz[1+4*(j+kj+l1)] = zero;  */
/*       fxyz[2+4*(j+kj+l1)] = zero;  */
         _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[8+4*(j+kj+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[16+4*(j+kj+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[24+4*(j+kj+l1)],v_zero);
/*       fxyz[4*(j+k1+l1)] = zero;    */
/*       fxyz[1+4*(j+k1+l1)] = zero;  */
/*       fxyz[2+4*(j+k1+l1)] = zero;  */
         _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[8+4*(j+k1+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[16+4*(j+k1+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[24+4*(j+k1+l1)],v_zero);
/*       at1 = at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1])); */
         v_zt3 = _mm512_mul_ps(v_zt1,v_zt1);
         v_zt3 = _mm512_add_ps(v_zt3,_mm512_mul_ps(v_zt2,v_zt2));
         v_zt3 = _mm512_mul_ps(v_at1,v_zt3);
/* convert to double precision before accumulating */
         v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt3));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt3,78));
         v_wp = _mm512_add_pd(v_wp,v_d);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         at1 = crealf(ffc[j+kk])*cimagf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         zt1 = cimagf(q[j+kj]) - crealf(q[j+kj])*_Complex_I;
         zt2 = cimagf(q[j+k1]) - crealf(q[j+k1])*_Complex_I;
         fxyz[4*(j+kj)] = at2*zt1;
         fxyz[1+4*(j+kj)] = at3*zt1;
         fxyz[2+4*(j+kj)] = zero;
         fxyz[4*(j+k1)] = at2*zt2;
         fxyz[1+4*(j+k1)] = -at3*zt2;
         fxyz[2+4*(j+k1)] = zero;
         fxyz[4*(j+kj+l1)] = zero;
         fxyz[1+4*(j+kj+l1)] = zero;
         fxyz[2+4*(j+kj+l1)] = zero;
         fxyz[4*(j+k1+l1)] = zero;
         fxyz[1+4*(j+k1+l1)] = zero;
         fxyz[2+4*(j+k1+l1)] = zero;
         at1 = at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1]));
         wp += (double) at1;
      }
/*    sum2 += wp; */
      _mm512_store_pd(&dd[0],v_wp);
      for (j = 1; j < 8; j++) {
          dd[0] += dd[j];
      }
      sum2 += (wp + dd[0]);
   }
/* mode numbers kx = 0, nx/2 */
   wp = 0.0;
   v_wp = _mm512_setzero_pd();
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      at1 = crealf(ffc[kk])*cimagf(ffc[kk]);
      at3 = at1*dny*(float) k;
      zt1 = cimagf(q[kj]) - crealf(q[kj])*_Complex_I;
      fxyz[4*kj] = zero;
      fxyz[1+4*kj] = at3*zt1;
      fxyz[2+4*kj] = zero;
      fxyz[4*k1] = zero;
      fxyz[1+4*k1] = zero;
      fxyz[2+4*k1] = zero;
      fxyz[4*(kj+l1)] = zero;
      fxyz[1+4*(kj+l1)] = zero;
      fxyz[2+4*(kj+l1)] = zero;
      fxyz[4*(k1+l1)] = zero;
      fxyz[1+4*(k1+l1)] = zero;
      fxyz[2+4*(k1+l1)] = zero;
      at1 = at1*(q[kj]*conjf(q[kj]));
      wp += (double) at1;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
/* vector loop over elements in blocks of 8 */
   for (j = 0; j < nxhs; j+=8) {
/*    at1 = crealf(ffc[j])*cimagf(ffc[j]); */
      v_at1 = _mm512_load_ps((float *)&ffc[j]);
      v_at2 = (__m512)_mm512_shuffle_epi32((__m512i)v_at1,177);
      v_at1 = _mm512_mul_ps(v_at1,v_at2);
/*    at2 = at1*dnx*(float) j; */
      v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
      v_at2 = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_at2 = _mm512_mul_ps(v_at1,_mm512_mul_ps(v_dnx,v_at2));
/*    zt1 = cimagf(q[j]) - crealf(q[j])*_Complex_I; */
      v_zt1 = _mm512_load_ps((float *)&q[j]);
      v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),v_zero,
              v_zt1);
      v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/* zero out kx = 0 mode */
      if (j==0) {
         v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),v_zero);
      }
/*    fxyz[4*j] = at2*zt1; */
/*    fxyz[1+4*j] = zero;  */
/*    fxyz[2+4*j] = zero;  */
      a = _mm512_mul_ps(v_at2,v_zt1);
      b = v_zero;
      c = v_zero;
/* perform 4x16 transpose for fxyz field components */
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),c,78);
      f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),a,78);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),v_zero,
          177);
      b = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(3855),e,
          177);
      c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),v_zero,
          177);
      d = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(3855),f,
          177);
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      _mm512_store_ps((float *)&fxyz[4*j],a);
      _mm512_store_ps((float *)&fxyz[8+4*j],b);
      _mm512_store_ps((float *)&fxyz[16+4*j],c);
      _mm512_store_ps((float *)&fxyz[24+4*j],d);
/*    fxyz[4*(j+k1)] = zero;   */
/*    fxyz[1+4*(j+k1)] = zero; */
/*    fxyz[2+4*(j+k1)] = zero; */
      _mm512_store_ps((float *)&fxyz[4*(j+k1)],v_zero);
      _mm512_store_ps((float *)&fxyz[8+4*(j+k1)],v_zero);
      _mm512_store_ps((float *)&fxyz[16+4*(j+k1)],v_zero);
      _mm512_store_ps((float *)&fxyz[24+4*(j+k1)],v_zero);
/*    fxyz[4*(j+l1)] = zero;    */
/*    fxyz[1+4*(j+l1)] = zero;  */
/*    fxyz[2+4*(j+l1)] = zero; */
      _mm512_store_ps((float *)&fxyz[4*(j+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[8+4*(j+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[16+4*(j+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[24+4*(j+l1)],v_zero);
/*    fxyz[4*(j+k1+l1)] = zero;   */
/*    fxyz[1+4*(j+k1+l1)] = zero; */
/*    fxyz[2+4*(j+k1+l1)] = zero; */
      _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[8+4*(j+k1+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[16+4*(j+k1+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[24+4*(j+k1+l1)],v_zero);
/*    wp += at1*(q[j]*conjf(q[j])); */
      v_zt3 = _mm512_mul_ps(v_at1,_mm512_mul_ps(v_zt1,v_zt1));
/* convert to double precision before accumulating */
      v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt3));
      v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt3,78));
      v_wp = _mm512_add_pd(v_wp,v_d);
   }
/* loop over remaining elements */
   for (j = itn; j < nxh; j++) {
      at1 = crealf(ffc[j])*cimagf(ffc[j]);
      at2 = at1*dnx*(float) j;
      zt1 = cimagf(q[j]) - crealf(q[j])*_Complex_I;
      fxyz[4*j] = at2*zt1;
      fxyz[1+4*j] = zero;
      fxyz[2+4*j] = zero;
      fxyz[4*(j+k1)] = zero;
      fxyz[1+4*(j+k1)] = zero;
      fxyz[2+4*(j+k1)] = zero;
      fxyz[4*(j+l1)] = zero;
      fxyz[1+4*(j+l1)] = zero;
      fxyz[2+4*(j+l1)] = zero;
      fxyz[4*(j+k1+l1)] = zero;
      fxyz[1+4*(j+k1+l1)] = zero;
      fxyz[2+4*(j+k1+l1)] = zero;
      at1 = at1*(q[j]*conjf(q[j]));
      wp += (double) at1;
   }
   fxyz[0] = zero;
   fxyz[1] = zero;
   fxyz[2] = zero;
   fxyz[4*k1] = zero;
   fxyz[1+4*k1] = zero;
   fxyz[2+4*k1] = zero;
   fxyz[4*l1] = zero;
   fxyz[1+4*l1] = zero;
   fxyz[2+4*l1] = zero;
   fxyz[4*(k1+l1)] = zero;
   fxyz[1+4*(k1+l1)] = zero;
   fxyz[2+4*(k1+l1)] = zero;
/* sum2 += wp; */
   _mm512_store_pd(&dd[0],v_wp);
   for (j = 1; j < 8; j++) {
      dd[0] += dd[j];
   }
   sum2 += (wp + dd[0]);
/* *we = wp*((float) nx)*((float) ny)*((float) nz); */
   *we = (sum1 + sum2)*((float) nx)*((float) ny)*((float) nz);
   return;
}

/*--------------------------------------------------------------------*/
void ckncmcuperp3(float complex cu[], int nx, int ny, int nz, int nxvh,
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
   cuy([kz][ky][kx] = cuy[kz][ky][kx]
                    - ky*(kx*cux[kz][ky][kx]+ky*cuy[kz][ky][kx]
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
   requires KNC, cu need to be 64 byte aligned
   nxhd needs to be a multiple of 8
   nxvh needs to be a multiple of 2
   cu needs to have 4 components
local data                                                 */
   int nxh, nyh, nzh, nxhs, itn, j, k, l, k1, l1, kj, lj, nxvyh;
   float dnx, dny, dnz, dkx, dky, dkz, dky2, dkz2, dkyz2, at1;
   float complex zero, zt1;
   __m512i v_j, v_it;
   __m512 v_dnx, v_dny, v_dnz, v_dkx, v_dky, v_dkz, v_dkz2, v_dkyz2;
   __m512 v_dk, v_at1, v_zt1, v_zt2, v_zero, v_one, v_at, v_as;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxhs = 2*(nxh/2);
   itn = 1 > nxhs ? 1 : nxhs;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   zero = 0.0 + 0.0*_Complex_I;
   v_j = _mm512_set_epi32(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0);
   v_dnx = _mm512_set1_ps(dnx);
   v_dny = _mm512_set1_ps(dny);
   v_dnz = _mm512_set1_ps(dnz);
   v_zero = _mm512_setzero_ps();
   v_one = _mm512_set1_ps(1.0f);
/* calculate transverse part of current */
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,k1,l1,lj,kj,dkx,dky,dkz,dkz2,dkyz2,at1,zt1,v_it,v_dk, \
v_dkx,v_dky,v_dkz,v_dkz2,v_dkyz2,v_at1,v_zt1,v_zt2,v_at,v_as)
      for (l = 1; l < nzh; l++) {
         dkz = dnz*(float) l;
         v_dkz = _mm512_cvtfxpnt_round_adjustepi32_ps(_mm512_set1_epi32(l),
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkz = _mm512_mul_ps(v_dnz,v_dkz);
         lj = nxvyh*l;
         l1 = nxvyh*nz - lj;
         dkz2 = dkz*dkz;
         v_dkz2 = _mm512_set1_ps(dkz2);
/* add kz to gradient operator */
         v_dk  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(12336),v_dkz);
         for (k = 1; k < nyh; k++) {
            dky = dny*(float) k;
            v_it = _mm512_set1_epi32(k);
            v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_dky = _mm512_mul_ps(v_dny,v_dky);
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            dkyz2 = dky*dky + dkz2;
            v_dkyz2 = _mm512_fmadd_ps(v_dky,v_dky,v_dkz2);
/* add ky to gradient operator */
            v_dk  = _mm512_mask_mov_ps(v_dk,_mm512_int2mask(3084),v_dky);
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nxhs; j+=2) {
/*             dkx = dnx*(float) j; */
               v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
               v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                       _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
               v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/*             at1 = 1.0/(dkx*dkx + dkyz2); */
               v_at1 = _mm512_fmadd_ps(v_dkx,v_dkx,v_dkyz2);
               v_at1 = _mm512_div_ps(v_one,v_at1);
/* add kx to gradient operator */
               v_dk  = _mm512_mask_mov_ps(v_dk,_mm512_int2mask(771),v_dkx);
/*             zt1 = at1*(dkx*cu[4*(j+kj+lj)] + dky*cu[1+4*(j+kj+lj)] */
/*                      + dkz*cu[2+4*(j+kj+lj)]);                     */
               v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj+lj)]);
               v_zt1 = _mm512_mul_ps(v_dk,v_zt2);
               v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
               v_zt1 = _mm512_add_ps(v_at,v_zt1);
               v_at = _mm512_permute4f128_ps(v_zt1,177);
               v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                          v_zero);
               }
/*             cu[4*(j+kj+lj)] -= dkx*zt1;   */
/*             cu[1+4*(j+kj+lj)] -= dky*zt1; */
/*             cu[2+4*(j+kj+lj)] -= dkz*zt1; */
               v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_dk,v_zt1));
               _mm512_store_ps((float *)&cu[4*(j+kj+lj)],v_zt2);
/*             zt1 = at1*(dkx*cu[4*(j+k1+lj)] - dky*cu[1+4*(j+k1+lj)] */
/*                      + dkz*cu[2+4*(j+k1+lj)]);                     */
               v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1+lj)]);
               v_as = _mm512_mask_sub_ps(v_dk,_mm512_int2mask(3084),v_zero,
                      v_dk);
               v_zt1 = _mm512_mul_ps(v_as,v_zt2);
               v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
               v_zt1 = _mm512_add_ps(v_at,v_zt1);
               v_at = _mm512_permute4f128_ps(v_zt1,177);
               v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                          v_zero);
               }
/*             cu[4*(j+k1+lj)] -= dkx*zt1;   */
/*             cu[1+4*(j+k1+lj)] += dky*zt1; */
/*             cu[2+4*(j+k1+lj)] -= dkz*zt1; */
               v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_as,v_zt1));
               _mm512_store_ps((float *)&cu[4*(j+k1+lj)],v_zt2);
/*             zt1 = at1*(dkx*cu[4*(j+kj+l1)] + dky*cu[1+4*(j+kj+l1)] */
/*                      - dkz*cu[2+4*(j+kj+l1)]);                     */
               v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj+l1)]);
               v_as = _mm512_mask_sub_ps(v_dk,_mm512_int2mask(12336),
                      v_zero,v_dk);
               v_zt1 = _mm512_mul_ps(v_as,v_zt2);
               v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
               v_zt1 = _mm512_add_ps(v_at,v_zt1);
               v_at = _mm512_permute4f128_ps(v_zt1,177);
               v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                          v_zero);
               }
/*             cu[4*(j+kj+l1)] -= dkx*zt1;   */
/*             cu[1+4*(j+kj+l1)] -= dky*zt1; */
/*             cu[2+4*(j+kj+l1)] += dkz*zt1; */
               v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_as,v_zt1));
               _mm512_store_ps((float *)&cu[4*(j+kj+l1)],v_zt2);
/*             zt1 = at1*(dkx*cu[4*(j+k1+l1)] - dky*cu[1+4*(j+k1+l1)] */
/*                      - dkz*cu[2+4*(j+k1+l1)]);                     */
               v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1+l1)]);
               v_as = _mm512_mask_sub_ps(v_dk,_mm512_int2mask(15420),
                      v_zero,v_dk);
               v_zt1 = _mm512_mul_ps(v_as,v_zt2);
               v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
               v_zt1 = _mm512_add_ps(v_at,v_zt1);
               v_at = _mm512_permute4f128_ps(v_zt1,177);
               v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                          v_zero);
               }
/*             cu[4*(j+k1+l1)] -= dkx*zt1;   */
/*             cu[1+4*(j+k1+l1)] += dky*zt1; */
/*             cu[2+4*(j+k1+l1)] += dkz*zt1; */
               v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_as,v_zt1));
               _mm512_store_ps((float *)&cu[4*(j+k1+l1)],v_zt2);
            }
/* loop over remaining elements */
            for (j = itn; j < nxh; j++) {
               dkx = dnx*(float) j;
               at1 = 1.0/(dkx*dkx + dkyz2);
               zt1 = at1*(dkx*cu[4*(j+kj+lj)] + dky*cu[1+4*(j+kj+lj)]
                        + dkz*cu[2+4*(j+kj+lj)]);
               cu[4*(j+kj+lj)] -= dkx*zt1;
               cu[1+4*(j+kj+lj)] -= dky*zt1;
               cu[2+4*(j+kj+lj)] -= dkz*zt1;
               zt1 = at1*(dkx*cu[4*(j+k1+lj)] - dky*cu[1+4*(j+k1+lj)]
                        + dkz*cu[2+4*(j+k1+lj)]);
               cu[4*(j+k1+lj)] -= dkx*zt1;
               cu[1+4*(j+k1+lj)] += dky*zt1;
               cu[2+4*(j+k1+lj)] -= dkz*zt1;
               zt1 = at1*(dkx*cu[4*(j+kj+l1)] + dky*cu[1+4*(j+kj+l1)]
                        - dkz*cu[2+4*(j+kj+l1)]);
               cu[4*(j+kj+l1)] -= dkx*zt1;
               cu[1+4*(j+kj+l1)] -= dky*zt1;
               cu[2+4*(j+kj+l1)] += dkz*zt1;
               zt1 = at1*(dkx*cu[4*(j+k1+l1)] - dky*cu[1+4*(j+k1+l1)]
                        - dkz*cu[2+4*(j+k1+l1)]);
               cu[4*(j+k1+l1)] -= dkx*zt1;
               cu[1+4*(j+k1+l1)] += dky*zt1;
               cu[2+4*(j+k1+l1)] += dkz*zt1;
            }
         }
/* mode numbers kx = 0, nx/2 */
         for (k = 1; k < nyh; k++) {
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            dky = dny*(float) k;
            at1 = 1.0/(dky*dky + dkz2);
            zt1 = at1*(dky*cu[1+4*(kj+lj)] + dkz*cu[2+4*(kj+lj)]);
            cu[1+4*(kj+lj)] -= dky*zt1;
            cu[2+4*(kj+lj)] -= dkz*zt1;
            cu[4*(k1+lj)] = zero;
            cu[1+4*(k1+lj)] = zero;
            cu[2+4*(k1+lj)] = zero;
            zt1 = at1*(dky*cu[1+4*(kj+l1)] - dkz*cu[2+4*(kj+l1)]);
            cu[1+4*(kj+l1)] -= dky*zt1;
            cu[2+4*(kj+l1)] += dkz*zt1;
            cu[4*(k1+l1)] = zero;
            cu[1+4*(k1+l1)] = zero;
            cu[2+4*(k1+l1)] = zero;
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
/* add ky to gradient operator */
         v_dk  = _mm512_mask_mov_ps(v_dk,_mm512_int2mask(3084),v_zero);
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          dkx = dnx*(float) j; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/*          at1 = 1.0/(dkx*dkx + dkz2); */
            v_at1 = _mm512_fmadd_ps(v_dkx,v_dkx,v_dkz2);
            v_at1 = _mm512_div_ps(v_one,v_at1);
/* add kx to gradient operator */
            v_dk  = _mm512_mask_mov_ps(v_dk,_mm512_int2mask(771),v_dkx);
/*          zt1 = at1*(dkx*cu[4*(j+lj)] + dkz*cu[2+4*(j+lj)]); */
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+lj)]);
            v_zt1 = _mm512_mul_ps(v_dk,v_zt2);
            v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
            v_zt1 = _mm512_add_ps(v_at,v_zt1);
            v_at = _mm512_permute4f128_ps(v_zt1,177);
            v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                       v_zero);
            }
/*          cu[4*(j+lj)] -= dkx*zt1;   */
/*          cu[2+4*(j+lj)] -= dkz*zt1; */
            v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_dk,v_zt1));
            _mm512_store_ps((float *)&cu[4*(j+lj)],v_zt2);
/*          cu[4*(j+k1+lj)] = zero;    */
/*          cu[1+4*(j+k1+lj)] = zero; */
/*          cu[2+4*(j+k1+lj)] = zero; */
            _mm512_store_ps((float *)&cu[4*(j+k1+lj)],v_zero);
/*          zt1 = at1*(dkx*cu[4*(j+l1)] - dkz*cu[2+4*(j+l1)]); */
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+l1)]);
            v_as = _mm512_mask_sub_ps(v_dk,_mm512_int2mask(12336),
                   v_zero,v_dk);
            v_zt1 = _mm512_mul_ps(v_as,v_zt2);
            v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
            v_zt1 = _mm512_add_ps(v_at,v_zt1);
            v_at = _mm512_permute4f128_ps(v_zt1,177);
            v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);
            }
/*          cu[4*(j+l1)] -= dkx*zt1;   */
/*          cu[2+4*(j+l1)] += dkz*zt1; */
            v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_as,v_zt1));
            _mm512_store_ps((float *)&cu[4*(j+l1)],v_zt2);
/*          cu[4*(j+k1+l1)] = zero;   */
/*          cu[1+4*(j+k1+l1)] = zero; */
/*          cu[2+4*(j+k1+l1)] = zero; */
            _mm512_store_ps((float *)&cu[4*(j+k1+l1)],v_zero);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            dkx = dnx*(float) j;
            at1 = 1.0/(dkx*dkx + dkz2);
            zt1 = at1*(dkx*cu[4*(j+lj)] + dkz*cu[2+4*(j+lj)]);
            cu[4*(j+lj)] -= dkx*zt1;
            cu[2+4*(j+lj)] -= dkz*zt1;
            cu[4*(j+k1+lj)] = zero;
            cu[1+4*(j+k1+lj)] = zero;
            cu[2+4*(j+k1+lj)] = zero;
            zt1 = at1*(dkx*cu[4*(j+l1)] - dkz*cu[2+4*(j+l1)]);
            cu[4*(j+l1)] -= dkx*zt1;
            cu[2+4*(j+l1)] += dkz*zt1;
            cu[4*(j+k1+l1)] = zero;
            cu[1+4*(j+k1+l1)] = zero;
            cu[2+4*(j+k1+l1)] = zero;
         }
/* mode numbers kx = 0, nx/2 */
         cu[2+4*lj] = zero;
         cu[4*(k1+lj)] = zero;
         cu[1+4*(k1+lj)] = zero;
         cu[2+4*(k1+lj)] = zero;
         cu[4*l1] = zero;
         cu[1+4*l1] = zero;
         cu[2+4*l1] = zero;
         cu[4*(k1+l1)] = zero;
         cu[1+4*(k1+l1)] = zero;
         cu[2+4*(k1+l1)] = zero;
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
#pragma omp parallel for \
private(j,k,k1,kj,dky,dky2,dkx,at1,zt1,v_it,v_dk,v_dkx,v_dky,v_dkyz2, \
v_at1,v_zt1,v_zt2,v_at,v_as)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      v_it = _mm512_set1_epi32(k);
      v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dky = _mm512_mul_ps(v_dny,v_dky);
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      dky2 = dky*dky;
      v_dkyz2 = _mm512_mul_ps(v_dky,v_dky);
/* add ky to gradient operator */
      v_dk  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(3084),v_dky);
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       dkx = dnx*(float) j; */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/*       at1 = 1.0/(dkx*dkx + dky2); */
         v_at1 = _mm512_fmadd_ps(v_dkx,v_dkx,v_dkyz2);
         v_at1 = _mm512_div_ps(v_one,v_at1);
/* add kx to gradient operator */
         v_dk  = _mm512_mask_mov_ps(v_dk,_mm512_int2mask(771),v_dkx);
/*       zt1 = at1*(dkx*cu[4*(j+kj)] + dky*cu[1+4*(j+kj)]); */
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj)]);
         v_zt1 = _mm512_mul_ps(v_dk,v_zt2);
         v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
         v_zt1 = _mm512_add_ps(v_at,v_zt1);
         v_at = _mm512_permute4f128_ps(v_zt1,177);
         v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);
         }
/*       cu[4*(j+kj)] -= dkx*zt1;   */
/*       cu[1+4*(j+kj)] -= dky*zt1; */
         v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_dk,v_zt1));
         _mm512_store_ps((float *)&cu[4*(j+kj)],v_zt2);
/*       zt1 = at1*(dkx*cu[4*(j+k1)]- dky*cu[1+4*(j+k1)]); */
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1)]);
         v_as = _mm512_mask_sub_ps(v_dk,_mm512_int2mask(3084),v_zero,
                v_dk);
         v_zt1 = _mm512_mul_ps(v_as,v_zt2);
         v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
         v_zt1 = _mm512_add_ps(v_at,v_zt1);
         v_at = _mm512_permute4f128_ps(v_zt1,177);
         v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);;
         }
/*       cu[4*(j+k1)] -= dkx*zt1;   */
/*       cu[1+4*(j+k1)] += dky*zt1; */
         v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_as,v_zt1));
         _mm512_store_ps((float *)&cu[4*(j+k1)],v_zt2);
/*       cu[4*(j+kj+l1)] = zero;   */
/*       cu[1+4*(j+kj+l1)] = zero; */
/*       cu[2+4*(j+kj+l1)] = zero; */
         _mm512_store_ps((float *)&cu[4*(j+kj+l1)],v_zero);
/*       cu[4*(j+k1+l1)] = zero;   */
/*       cu[1+4*(j+k1+l1)] = zero; */
/*       cu[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&cu[4*(j+k1+l1)],v_zero);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         dkx = dnx*(float) j;
         at1 = 1.0/(dkx*dkx + dky2);
         zt1 = at1*(dkx*cu[4*(j+kj)] + dky*cu[1+4*(j+kj)]);
         cu[4*(j+kj)] -= dkx*zt1;
         cu[1+4*(j+kj)] -= dky*zt1;
         zt1 = at1*(dkx*cu[4*(j+k1)]- dky*cu[1+4*(j+k1)]);
         cu[4*(j+k1)] -= dkx*zt1;
         cu[1+4*(j+k1)] += dky*zt1;
         cu[4*(j+kj+l1)] = zero;
         cu[1+4*(j+kj+l1)] = zero;
         cu[2+4*(j+kj+l1)] = zero;
         cu[4*(j+k1+l1)] = zero;
         cu[1+4*(j+k1+l1)] = zero;
         cu[2+4*(j+k1+l1)] = zero;
      }

   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      cu[1+4*kj] = zero;
      cu[4*k1] = zero;
      cu[1+4*k1] = zero;
      cu[2+4*k1] = zero;
      cu[4*(kj+l1)] = zero;
      cu[1+4*(kj+l1)] = zero;
      cu[2+4*(kj+l1)] = zero;
      cu[4*(k1+l1)] = zero;
      cu[1+4*(k1+l1)] = zero;
      cu[2+4*(k1+l1)] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
   for (j = 0; j < nxhs; j+=2) {
      v_zt2 = _mm512_load_ps((float *)&cu[4*j]);
/* zero out kx = 0 mode */
      if (j==0) {
         v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(255),v_zero);
      }
/*    cu[4*j] = zero;      */
      v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(771),v_zero);
      _mm512_store_ps((float *)&cu[4*j],v_zt2);
/*    cu[4*(j+k1)] = zero; */
/*    cu[1+4*(j+k1)] = zero; */
/*    cu[2+4*(j+k1)] = zero; */
      _mm512_store_ps((float *)&cu[4*(j+k1)],v_zero);
/*    cu[4*(j+l1)] = zero;   */
/*    cu[1+4*(j+l1)] = zero; */
/*    cu[2+4*(j+l1)] = zero; */
      _mm512_store_ps((float *)&cu[4*(j+l1)],v_zero);
/*    cu[4*(j+k1+l1)] = zero;   */
/*    cu[1+4*(j+k1+l1)] = zero; */
/*    cu[2+4*(j+k1+l1)] = zero; */
      _mm512_store_ps((float *)&cu[4*(j+k1+l1)],v_zero);
   }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
      cu[4*j] = zero;
      cu[4*(j+k1)] = zero;
      cu[1+4*(j+k1)] = zero;
      cu[2+4*(j+k1)] = zero;
      cu[4*(j+l1)] = zero;
      cu[1+4*(j+l1)] = zero;
      cu[2+4*(j+l1)] = zero;
      cu[4*(j+k1+l1)] = zero;
      cu[1+4*(j+k1+l1)] = zero;
      cu[2+4*(j+k1+l1)] = zero;
   }
   cu[0] = zero;
   cu[1] = zero;
   cu[2] = zero;
   cu[4*k1] = zero;
   cu[1+4*k1] = zero;
   cu[2+4*k1] = zero;
   cu[4*l1] = zero;
   cu[1+4*l1] = zero;
   cu[2+4*l1] = zero;
   cu[4*(k1+l1)] = zero;
   cu[1+4*(k1+l1)] = zero;
   cu[2+4*(k1+l1)] = zero;
   return;
}

/*--------------------------------------------------------------------*/
void ckncmibpois33(float complex cu[], float complex bxyz[],
                   float complex ffc[], float ci, float *wm, int nx,
                   int ny, int nz, int nxvh, int nyv, int nzv, int nxhd,
                   int nyhd, int nzhd) {
/* this subroutine solves 3d poisson's equation in fourier space for
   magnetic field with periodic boundary conditions.
   input: cu,ffc,ci,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
   output: bxyz, wm
   approximate flop count is:
   193*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   the magnetic field is calculated using the equations:
   bx[kz][ky][kx] = ci*ci*sqrt(-1)*g[kz][ky][kx]*
                  (ky*cuz[kz][ky][kx]-kz*cuy[kz][ky][kx]),
   by[kz][ky][kx] = ci*ci*sqrt(-1)*g[kz][ky][kx]*
                  (kz*cux[kz][ky][kx]-kx*cuz[kz][ky][kx]),
   bz[kz][ky][kx] = ci*ci*sqrt(-1)*g[kz][ky][kx]*
                  (kx*cuy[kz][ky][kx]-ky*cux[kz][ky][kx]),
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
   j,k,l = fourier mode numbers,
   g[kz][ky][kx] = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
   s[kz][ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
   bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
   bx(ky=pi) = by(ky=pi) = bx(ky=pi) = 0,
   bx(kz=pi) = by(kz=pi) = bz(kz=pi) = 0,
   bx(kx=0,ky=0,kz=0) = by(kx=0,ky=0,kz=0) = bz(kx=0,ky=0,kz=0) = 0.
   cu[l][k][j][i] = complex current density for fourier mode (j,k,l)
   bxyz[l][k][j][i] = i component of complex magnetic field
   all for fourier mode (j,k,l)
   aimag(ffc(j,k,l)) = finite-size particle shape factor s
   for fourier mode (j,k,l)
   real(ffc(j,k,l)) = potential green's function g
   for fourier mode (j,k,l)
   ci = reciprocal of velocity of light
   magnetic field energy is also calculated, using
   wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
      |cu[kz][ky][kx]*s[kz][ky][kx]|**2)
   this expression is valid only if the current is divergence-free
   nx/ny/nz = system length in x/y/z direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nzv = fourth dimension of field arrays, must be >= nz
   nxhd = dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
   nzhd = third dimension of form factor array, must be >= nzh
   requires KNC, cu, bxyz, ffc need to be 64 byte aligned
   nxhd needs to be a multiple of 8
   nxvh needs to be a multiple of 2
   cu, bxyz need to have 4 components
local data                                                 */
   int nxh, nyh, nzh, nxhs, itn, j, k, l, k1, l1, kk, kj, ll, lj;
   int nxyhd, nxvyh;
   float dnx, dny, dnz, dky, dkz, ci2, at1, at2, at3, at4;
   float complex zero, zt1, zt2, zt3;
   double wp, sum1, sum2;
   __m512i v_j, v_it, v_n, v_m;
   __m512 v_dnx, v_dny, v_dnz, v_dkx, v_dky, v_dkz, v_ci2;
   __m512 v_dk1, v_dk2, v_at1, v_at2, v_at3, v_at4, v_zero;
   __m512 v_zt1, v_zt2, v_zt3, v_zt4;
   __m512d v_wp, v_d;
   __attribute__((aligned(64))) double dd[8];
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxhs = 2*(nxh/2);
   itn = 1 > nxhs ? 1 : nxhs;
   nxyhd = nxhd*nyhd;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   zero = 0.0 + 0.0*_Complex_I;
   ci2 = ci*ci;
   v_j = _mm512_set_epi32(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0);
   v_n =  _mm512_set_epi32(15,14,11,10,9,8,13,12,7,6,3,2,1,0,5,4);
   v_m =  _mm512_set_epi32(15,14,9,8,13,12,11,10,7,6,1,0,5,4,3,2);
   v_dnx = _mm512_set1_ps(dnx);
   v_dny = _mm512_set1_ps(dny);
   v_dnz = _mm512_set1_ps(dnz);
   v_zero = _mm512_setzero_ps();
   v_ci2 = _mm512_set1_ps(ci2);
/* calculate magnetic field and sum field energy */
   sum1 = 0.0;
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,k1,l1,ll,lj,kk,kj,dky,dkz,at1,at2,at3,at4,zt1,zt2,zt3, \
wp,v_it,v_dkx,v_dky,v_dkz,v_dk1,v_dk2,v_at1,v_at2,v_at3,v_at4,v_zt1, \
v_zt2,v_zt3,v_zt4,v_d,v_wp,dd) \
reduction(+:sum1)
      for (l = 1; l < nzh; l++) {
         dkz = dnz*(float) l;
         v_dkz = _mm512_cvtfxpnt_round_adjustepi32_ps(_mm512_set1_epi32(l),
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkz = _mm512_mul_ps(v_dnz,v_dkz);
         ll = nxyhd*l;
         lj = nxvyh*l;
         l1 = nxvyh*nz - lj;
         wp = 0.0;
         v_wp = _mm512_set1_pd(0.0);
/* add kz to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(771),v_dkz);
         v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(3084),v_dkz);
         for (k = 1; k < nyh; k++) {
            dky = dny*(float) k;
            v_it = _mm512_set1_epi32(k);
            v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_dky = _mm512_mul_ps(v_dny,v_dky);
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
/* add ky to curl operators */
            v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(12336),
                     v_dky);
            v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(771),
                     v_dky);
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nxhs; j+=2) {
/*             at1 = ci2*crealf(ffc[j+kk+ll]); */
               v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                       _mm512_int2mask(15),(float *)&ffc[j+kk+ll]);
               v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                       _mm512_int2mask(15),(float *)&ffc[j+kk+ll+8]);
               v_at1 = _mm512_permute4f128_ps(v_at1,0);
               v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                       _mm512_int2mask(13260),(__m512i)v_at1,78);
               v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                       _mm512_int2mask(43690),(__m512i)v_at4,177);
               v_at1 = _mm512_mul_ps(v_ci2,v_at1);
/*             at2 = at1*dnx*(float) j; */
/*             at3 = dky*at1;           */
/*             at4 = dkz*at1;           */
               v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
               v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                       _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
               v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
               v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                        v_dkx);
               v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                        v_dkx);
/* normalize curl operators */
               v_at2 = _mm512_mul_ps(v_at1,v_dk1);
               v_at3 = _mm512_mul_ps(v_at1,v_dk2);
/*             at1 = at1*cimagf(ffc[j+kk+ll]); */
               v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                       _mm512_int2mask(21845),(__m512i)v_at4,177);
               v_at1 = _mm512_mul_ps(v_at1,v_at4);
/*             zt1 = -cimagf(cu[2+4*(j+kj+lj)])              */
/*                  + crealf(cu[2+4*(j+kj+lj)])*_Complex_I;/ */
/*             zt2 = -cimagf(cu[1+4*(j+kj+lj)])              */
/*                  + crealf(cu[1+4*(j+kj+lj)])*_Complex_I;  */
/*             zt3 = -cimagf(cu[4*(j+kj+lj)])                */
/*                  + crealf(cu[4*(j+kj+lj)])*_Complex_I;    */
               v_zt3 = _mm512_load_ps((float *)&cu[4*(j+kj+lj)]);
               v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),
                       v_zero,v_zt3);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*             bxyz[4*(j+kj+lj)] = at3*zt1 - at4*zt2;   */
/*             bxyz[1+4*(j+kj+lj)] = at4*zt3 - at2*zt1; */
/*             bxyz[2+4*(j+kj+lj)] = at2*zt2 - at3*zt3; */
               v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
               v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                          v_zero);
                  v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                          v_zero);
               }
               _mm512_store_ps((float *)&bxyz[4*(j+kj+lj)],v_zt1);
/*             wp += at1*(cu[4*(j+kj+lj)]*conjf(cu[4*(j+kj+lj)]) */
/*                + cu[1+4*(j+kj+lj)]*conjf(cu[1+4*(j+kj+lj)])   */
/*                + cu[2+4*(j+kj+lj)]*conjf(cu[2+4*(j+kj+lj)])); */
               v_zt4 = _mm512_mul_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                       _mm512_int2mask(16191),v_zt3,v_zt3));
/*             zt1 = -cimagf(cu[2+4*(j+k1+lj)])             */
/*                  + crealf(cu[2+4*(j+k1+lj)])*_Complex_I; */
/*             zt2 = -cimagf(cu[1+4*(j+k1+lj)])             */
/*                  + crealf(cu[1+4*(j+k1+lj)])*_Complex_I; */
/*             zt3 = -cimagf(cu[4*(j+k1+lj)])               */
/*                  + crealf(cu[4*(j+k1+lj)])*_Complex_I;   */
               v_zt3 = _mm512_load_ps((float *)&cu[4*(j+k1+lj)]);
               v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),
                       v_zero,v_zt3);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
               v_zt1 = _mm512_mask_sub_ps(v_at2,_mm512_int2mask(12336),
                       v_zero,v_at2);
               v_zt2 = _mm512_mask_sub_ps(v_at3,_mm512_int2mask(771),
                       v_zero,v_at3);
/*             bxyz[4*(j+k1+lj)] = -at3*zt1 - at4*zt2;  */
/*             bxyz[1+4*(j+k1+lj)] = at4*zt3 - at2*zt1; */
/*             bxyz[2+4*(j+k1+lj)] = at2*zt2 + at3*zt3; */
               v_zt1 = _mm512_mul_ps(v_zt1,v_zt3);
               v_zt2 = _mm512_mul_ps(v_zt2,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                          v_zero);
                  v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                          v_zero);
               }
               _mm512_store_ps((float *)&bxyz[4*(j+k1+lj)],v_zt1);
/*             wp += at1*(cu[4*(j+k1+lj)]*conjf(cu[4*(j+k1+lj)]) */
/*                + cu[1+4*(j+k1+lj)]*conjf(cu[1+4*(j+k1+lj)])   */
/*                + cu[2+4*(j+k1+lj)]*conjf(cu[2+4*(j+k1+lj)])); */
               v_zt4 = _mm512_fmadd_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                       _mm512_int2mask(16191),v_zt3,v_zt3),v_zt4);
/*             zt1 = -cimagf(cu[2+4*(j+kj+l1)])             */
/*                  + crealf(cu[2+4*(j+kj+l1)])*_Complex_I; */
/*             zt2 = -cimagf(cu[1+4*(j+kj+l1)])             */
/*                  + crealf(cu[1+4*(j+kj+l1)])*_Complex_I; */
/*             zt3 = -cimagf(cu[4*(j+kj+l1)])               */
/*                  + crealf(cu[4*(j+kj+l1)])*_Complex_I;   */
               v_zt3 = _mm512_load_ps((float *)&cu[4*(j+kj+l1)]);
               v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),
                       v_zero,v_zt3);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
               v_zt1 = _mm512_mask_sub_ps(v_at2,_mm512_int2mask(771),
                       v_zero,v_at2);
               v_zt2 = _mm512_mask_sub_ps(v_at3,_mm512_int2mask(3084),
                       v_zero,v_at3);
/*             bxyz[4*(j+kj+l1)] = at3*zt1 + at4*zt2;    */
/*             bxyz[1+4*(j+kj+l1)] = -at4*zt3 - at2*zt1; */
/*             bxyz[2+4*(j+kj+l1)] = at2*zt2 - at3*zt3;  */
               v_zt1 = _mm512_mul_ps(v_zt1,v_zt3);
               v_zt2 = _mm512_mul_ps(v_zt2,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                          v_zero);
                  v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                          v_zero);
               }
               _mm512_store_ps((float *)&bxyz[4*(j+kj+l1)],v_zt1);
/*             wp += at1*(cu[4*(j+kj+l1)]*conjf(cu[4*(j+kj+l1)]) */
/*                + cu[1+4*(j+kj+l1)]*conjf(cu[1+4*(j+kj+l1)])   */
/*                + cu[2+4*(j+kj+l1)]*conjf(cu[2+4*(j+kj+l1)])); */
               v_zt4 = _mm512_fmadd_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                       _mm512_int2mask(16191),v_zt3,v_zt3),v_zt4);
/*             zt1 = -cimagf(cu[2+4*(j+k1+l1)])             */
/*                  + crealf(cu[2+4*(j+k1+l1)])*_Complex_I; */
/*             zt2 = -cimagf(cu[1+4*(j+k1+l1)])             */
/*                  + crealf(cu[1+4*(j+k1+l1)])*_Complex_I; */
/*             zt3 = -cimagf(cu[4*(j+k1+l1)])               */
/*                  + crealf(cu[4*(j+k1+l1)])*_Complex_I;   */
               v_zt3 = _mm512_load_ps((float *)&cu[4*(j+k1+l1)]);
               v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),
                       v_zero,v_zt3);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
               v_zt1 = _mm512_mask_sub_ps(v_at2,_mm512_int2mask(13107),
                       v_zero,v_at2);
               v_zt2 = _mm512_mask_sub_ps(v_at3,_mm512_int2mask(3855),
                       v_zero,v_at3);
/*             bxyz[4*(j+k1+l1)] = -at3*zt1 + at4*zt2;   */
/*             bxyz[1+4*(j+k1+l1)] = -at4*zt3 - at2*zt1; */
/*             bxyz[2+4*(j+k1+l1)] = at2*zt2 + at3*zt3;  */
               v_zt1 = _mm512_mul_ps(v_zt1,v_zt3);
               v_zt2 = _mm512_mul_ps(v_zt2,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                          v_zero);
                  v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                          v_zero);
               }
               _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zt1);
/*             wp += at1*(cu[4*(j+k1+l1)]*conjf(cu[4*(j+k1+l1)]) */
/*                + cu[1+4*(j+k1+l1)]*conjf(cu[1+4*(j+k1+l1)])   */
/*                + cu[2+4*(j+k1+l1)]*conjf(cu[2+4*(j+k1+l1)])); */
               v_zt4 = _mm512_fmadd_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                       _mm512_int2mask(16191),v_zt3,v_zt3),v_zt4);
/* convert to double precision before accumulating */
               v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt4));
               v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt4,78));
               v_wp = _mm512_add_pd(v_wp,v_d);
            }
/* loop over remaining elements */
            for (j = itn; j < nxh; j++) {
               at1 = ci2*crealf(ffc[j+kk+ll]);
               at2 = at1*dnx*(float) j;
               at3 = dky*at1;
               at4 = dkz*at1;
               at1 = at1*cimagf(ffc[j+kk+ll]);
               zt1 = -cimagf(cu[2+4*(j+kj+lj)])
                    + crealf(cu[2+4*(j+kj+lj)])*_Complex_I;
               zt2 = -cimagf(cu[1+4*(j+kj+lj)])
                    + crealf(cu[1+4*(j+kj+lj)])*_Complex_I;
               zt3 = -cimagf(cu[4*(j+kj+lj)])
                    + crealf(cu[4*(j+kj+lj)])*_Complex_I;
               bxyz[4*(j+kj+lj)] = at3*zt1 - at4*zt2;
               bxyz[1+4*(j+kj+lj)] = at4*zt3 - at2*zt1;
               bxyz[2+4*(j+kj+lj)] = at2*zt2 - at3*zt3;
               zt1 = -cimagf(cu[2+4*(j+k1+lj)])
                    + crealf(cu[2+4*(j+k1+lj)])*_Complex_I;
               zt2 = -cimagf(cu[1+4*(j+k1+lj)])
                    + crealf(cu[1+4*(j+k1+lj)])*_Complex_I;
               zt3 = -cimagf(cu[4*(j+k1+lj)])
                    + crealf(cu[4*(j+k1+lj)])*_Complex_I;
               bxyz[4*(j+k1+lj)] = -at3*zt1 - at4*zt2;
               bxyz[1+4*(j+k1+lj)] = at4*zt3 - at2*zt1;
               bxyz[2+4*(j+k1+lj)] = at2*zt2 + at3*zt3;
               zt1 = -cimagf(cu[2+4*(j+kj+l1)])
                    + crealf(cu[2+4*(j+kj+l1)])*_Complex_I;
               zt2 = -cimagf(cu[1+4*(j+kj+l1)])
                    + crealf(cu[1+4*(j+kj+l1)])*_Complex_I;
               zt3 = -cimagf(cu[4*(j+kj+l1)])
                    + crealf(cu[4*(j+kj+l1)])*_Complex_I;
               bxyz[4*(j+kj+l1)] = at3*zt1 + at4*zt2;
               bxyz[1+4*(j+kj+l1)] = -at4*zt3 - at2*zt1;
               bxyz[2+4*(j+kj+l1)] = at2*zt2 - at3*zt3;
               zt1 = -cimagf(cu[2+4*(j+k1+l1)])
                    + crealf(cu[2+4*(j+k1+l1)])*_Complex_I;
               zt2 = -cimagf(cu[1+4*(j+k1+l1)])
                    + crealf(cu[1+4*(j+k1+l1)])*_Complex_I;
               zt3 = -cimagf(cu[4*(j+k1+l1)])
                    + crealf(cu[4*(j+k1+l1)])*_Complex_I;
               bxyz[4*(j+k1+l1)] = -at3*zt1 + at4*zt2;
               bxyz[1+4*(j+k1+l1)] = -at4*zt3 - at2*zt1;
               bxyz[2+4*(j+k1+l1)] = at2*zt2 + at3*zt3;
               at1 = at1*(cu[4*(j+kj+lj)]*conjf(cu[4*(j+kj+lj)])
                   + cu[1+4*(j+kj+lj)]*conjf(cu[1+4*(j+kj+lj)])
                   + cu[2+4*(j+kj+lj)]*conjf(cu[2+4*(j+kj+lj)])
                   + cu[4*(j+k1+lj)]*conjf(cu[4*(j+k1+lj)])
                   + cu[1+4*(j+k1+lj)]*conjf(cu[1+4*(j+k1+lj)])
                   + cu[2+4*(j+k1+lj)]*conjf(cu[2+4*(j+k1+lj)])
                   + cu[4*(j+kj+l1)]*conjf(cu[4*(j+kj+l1)])
                   + cu[1+4*(j+kj+l1)]*conjf(cu[1+4*(j+kj+l1)])
                   + cu[2+4*(j+kj+l1)]*conjf(cu[2+4*(j+kj+l1)])
                   + cu[4*(j+k1+l1)]*conjf(cu[4*(j+k1+l1)])
                   + cu[1+4*(j+k1+l1)]*conjf(cu[1+4*(j+k1+l1)])
                   + cu[2+4*(j+k1+l1)]*conjf(cu[2+4*(j+k1+l1)]));
               wp += (double) at1;
           }
         }
/* mode numbers kx = 0, nx/2 */
         for (k = 1; k < nyh; k++) {
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            at1 = ci2*crealf(ffc[kk+ll]);
            at3 = at1*dny*(float) k;
            at4 = dkz*at1;
            at1 = at1*cimagf(ffc[kk+ll]);
            zt1 = -cimagf(cu[2+4*(kj+lj)])
                 + crealf(cu[2+4*(kj+lj)])*_Complex_I;
            zt2 = -cimagf(cu[1+4*(kj+lj)])
                 + crealf(cu[1+4*(kj+lj)])*_Complex_I;
            zt3 = -cimagf(cu[4*(kj+lj)])
                 + crealf(cu[4*(kj+lj)])*_Complex_I;
            bxyz[4*(kj+lj)] = at3*zt1 - at4*zt2;
            bxyz[1+4*(kj+lj)] = at4*zt3;
            bxyz[2+4*(kj+lj)] = -at3*zt3;
            bxyz[4*(k1+lj)] = zero;
            bxyz[1+4*(k1+lj)] = zero;
            bxyz[2+4*(k1+lj)] = zero;
            zt1 = -cimagf(cu[2+4*(kj+l1)])
                 + crealf(cu[2+4*(kj+l1)])*_Complex_I;
            zt2 = -cimagf(cu[1+4*(kj+l1)])
                 + crealf(cu[1+4*(kj+l1)])*_Complex_I;
            zt3 = -cimagf(cu[4*(kj+l1)])
                 + crealf(cu[4*(kj+l1)])*_Complex_I;
            bxyz[4*(kj+l1)] = at3*zt1 + at4*zt2;
            bxyz[1+4*(kj+l1)] = -at4*zt3;
            bxyz[2+4*(kj+l1)] = -at3*zt3;
            bxyz[4*(k1+l1)] = zero;
            bxyz[1+4*(k1+l1)]  = zero;
            bxyz[2+4*(k1+l1)] = zero;
            at1 = at1*(cu[4*(kj+lj)]*conjf(cu[4*(kj+lj)])
                + cu[1+4*(kj+lj)]*conjf(cu[1+4*(kj+lj)])
                + cu[2+4*(kj+lj)]*conjf(cu[2+4*(kj+lj)])
                + cu[4*(kj+l1)]*conjf(cu[4*(kj+l1)])
                + cu[1+4*(kj+l1)]*conjf(cu[1+4*(kj+l1)])
                + cu[2+4*(kj+l1)]*conjf(cu[2+4*(kj+l1)]));
            wp += (double) at1;
           }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
/* add ky to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(12336),v_zero);
         v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(771),v_zero);
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          at1 = ci2*crealf(ffc[j+ll]); */
            v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,_mm512_int2mask(15),
                    (float *)&ffc[j+ll]);
            v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,_mm512_int2mask(15),
                    (float *)&ffc[j+ll+8]);
            v_at1 = _mm512_permute4f128_ps(v_at1,0);
            v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                    _mm512_int2mask(13260),(__m512i)v_at1,78);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                    _mm512_int2mask(43690),(__m512i)v_at4,177);
            v_at1 = _mm512_mul_ps(v_ci2,v_at1);
/*          at2 = at1*dnx*(float) j;     */
/*          at4 = dkz*at1;               */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
            v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                     v_dkx);
            v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                     v_dkx);
/* normalize curl operators */
            v_at2 = _mm512_mul_ps(v_at1,v_dk1);
            v_at3 = _mm512_mul_ps(v_at1,v_dk2);
/*          at1 = at1*cimagf(ffc[j+ll]); */
            v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                    _mm512_int2mask(21845),(__m512i)v_at4,177);
            v_at1 = _mm512_mul_ps(v_at1,v_at4);
/*          zt1 = -cimagf(cu[2+4*(j+lj)])             */
/*               + crealf(cu[2+4*(j+lj)])*_Complex_I; */
/*          zt2 = -cimagf(cu[1+4*(j+lj)])             */
/*               + crealf(cu[1+4*(j+lj)])*_Complex_I; */
/*          zt3 = -cimagf(cu[4*(j+lj)])               */
/*               + crealf(cu[4*(j+lj)])*_Complex_I;   */
            v_zt3 = _mm512_load_ps((float *)&cu[4*(j+lj)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),v_zero,
                    v_zt3);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          bxyz[4*(j+lj)] = -at4*zt2;            */
/*          bxyz[1+4*(j+lj)] = at4*zt3 - at2*zt1; */
/*          bxyz[2+4*(j+lj)] = at2*zt2;           */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                       v_zero);
               v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                       v_zero);
            }
            _mm512_store_ps((float *)&bxyz[4*(j+lj)],v_zt1);
/*          wp += at1*(cu[4*(j+lj)]*conjf(cu[4*(j+lj)]) */
/*             + cu[1+4*(j+lj)]*conjf(cu[1+4*(j+lj)])   */
/*             + cu[2+4*(j+lj)]*conjf(cu[2+4*(j+lj)])   */
            v_zt4 = _mm512_mul_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                    _mm512_int2mask(16191),v_zt3,v_zt3));
/*          bxyz[4*(j+k1+lj)] = zero;   */
/*          bxyz[1+4*(j+k1+lj)] = zero; */
/*          bxyz[2+4*(j+k1+lj)] = zero; */
            _mm512_store_ps((float *)&bxyz[4*(j+k1+lj)],v_zero);
/*          zt1 = -cimagf(cu[2+4*(j+l1)])             */
/*               + crealf(cu[2+4*(j+l1)])*_Complex_I; */
/*          zt2 = -cimagf(cu[1+4*(j+l1)])             */
/*               + crealf(cu[1+4*(j+l1)])*_Complex_I; */
/*          zt3 = -cimagf(cu[4*(j+l1)])               */
/*               + crealf(cu[4*(j+l1)])*_Complex_I;   */
            v_zt3 = _mm512_load_ps((float *)&cu[4*(j+l1)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),v_zero,
                    v_zt3);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
            v_zt1 = _mm512_mask_sub_ps(v_at2,_mm512_int2mask(771),v_zero,
                    v_at2);
            v_zt2 = _mm512_mask_sub_ps(v_at3,_mm512_int2mask(3084),v_zero,
                    v_at3);
/*          bxyz[4*(j+l1)] = at4*zt2;              */
/*          bxyz[1+4*(j+l1)] = -at4*zt3 - at2*zt1; */
/*          bxyz[2+4*(j+l1)] = at2*zt2;            */
            v_zt1 = _mm512_mul_ps(v_zt1,v_zt3);
            v_zt2 = _mm512_mul_ps(v_zt2,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                       v_zero);
               v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                       v_zero);
            }
            _mm512_store_ps((float *)&bxyz[4*(j+l1)],v_zt1);
/*          wp += at1*(cu[4*(j+l1)]*conjf(cu[4*(j+l1)]) */ 
/*             + cu[1+4*(j+l1)]*conjf(cu[1+4*(j+l1)])   */
/*             + cu[2+4*(j+l1)]*conjf(cu[2+4*(j+l1)])); */
            v_zt4 = _mm512_fmadd_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                    _mm512_int2mask(16191),v_zt3,v_zt3),v_zt4);
/* convert to double precision before accumulating */
            v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt4));
            v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt4,78));
            v_wp = _mm512_add_pd(v_wp,v_d);
/*          bxyz[4*(j+k1+l1)] = zero;   */
/*          bxyz[1+4*(j+k1+l1)] = zero; */
/*          bxyz[2+4*(j+k1+l1)] = zero; */
            _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            at1 = ci2*crealf(ffc[j+ll]);
            at2 = at1*dnx*(float) j;  
            at4 = dkz*at1;
            at1 = at1*cimagf(ffc[j+ll]);
            zt1 = -cimagf(cu[2+4*(j+lj)])
                 + crealf(cu[2+4*(j+lj)])*_Complex_I;
            zt2 = -cimagf(cu[1+4*(j+lj)])
                 + crealf(cu[1+4*(j+lj)])*_Complex_I;
            zt3 = -cimagf(cu[4*(j+lj)])
                 + crealf(cu[4*(j+lj)])*_Complex_I;
            bxyz[4*(j+lj)] = -at4*zt2;
            bxyz[1+4*(j+lj)] = at4*zt3 - at2*zt1;
            bxyz[2+4*(j+lj)] = at2*zt2;
            bxyz[4*(j+k1+lj)] = zero;
            bxyz[1+4*(j+k1+lj)] = zero;
            bxyz[2+4*(j+k1+lj)] = zero;
            zt1 = -cimagf(cu[2+4*(j+l1)])
                 + crealf(cu[2+4*(j+l1)])*_Complex_I;
            zt2 = -cimagf(cu[1+4*(j+l1)])
                 + crealf(cu[1+4*(j+l1)])*_Complex_I;
            zt3 = -cimagf(cu[4*(j+l1)])
                 + crealf(cu[4*(j+l1)])*_Complex_I;
            bxyz[4*(j+l1)] = at4*zt2;
            bxyz[1+4*(j+l1)] = -at4*zt3 - at2*zt1;
            bxyz[2+4*(j+l1)] = at2*zt2;
            bxyz[4*(j+k1+l1)] = zero;
            bxyz[1+4*(j+k1+l1)] = zero;
            bxyz[2+4*(j+k1+l1)] = zero;
            at1 = at1*(cu[4*(j+lj)]*conjf(cu[4*(j+lj)])
                + cu[1+4*(j+lj)]*conjf(cu[1+4*(j+lj)])
                + cu[2+4*(j+lj)]*conjf(cu[2+4*(j+lj)])
                + cu[4*(j+l1)]*conjf(cu[4*(j+l1)])
                + cu[1+4*(j+l1)]*conjf(cu[1+4*(j+l1)])
                + cu[2+4*(j+l1)]*conjf(cu[2+4*(j+l1)]));
            wp += (double) at1;
         }
/* mode numbers kx = 0, nx/2 */
         at1 = ci2*crealf(ffc[ll]);
         at4 = dkz*at1;
         at1 = at1*cimagf(ffc[ll]);
         zt2 = -cimagf(cu[1+4*(lj)]) + crealf(cu[1+4*(lj)])*_Complex_I;
         zt3 = -cimagf(cu[4*(lj)]) + crealf(cu[4*(lj)])*_Complex_I;
         bxyz[4*lj] = -at4*zt2;
         bxyz[1+4*lj] = at4*zt3;
         bxyz[2+4*lj] = zero;
         bxyz[4*(k1+lj)] = zero;
         bxyz[1+4*(k1+lj)] = zero;
         bxyz[2+4*(k1+lj)] = zero;
         bxyz[4*l1] = zero;
         bxyz[1+4*l1] = zero;
         bxyz[2+4*l1] = zero;
         bxyz[4*(k1+l1)] = zero;
         bxyz[1+4*(k1+l1)] = zero;
         bxyz[2+4*(k1+l1)] = zero;
         at1 = at1*(cu[4*lj]*conjf(cu[4*lj])
             + cu[1+4*lj]*conjf(cu[1+4*lj])
             + cu[2+4*lj]*conjf(cu[2+4*lj]));
         wp += (double) at1;
/*       sum1 += wp; */
         _mm512_store_pd(&dd[0],v_wp);
         for (j = 1; j < 8; j++) {
            dd[0] += dd[j];
         }
         sum1 += (wp + dd[0]);
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   sum2 = 0.0;
#pragma omp parallel for \
private(j,k,k1,kk,kj,dky,at1,at2,at3,zt1,zt2,zt3,wp,v_it,v_dkx,v_dky, \
v_dk1,v_dk2,v_at1,v_at2,v_at3,v_at4,v_zt1,v_zt2,v_zt3,v_zt4,v_d,v_wp, \
dd) \
reduction(+:sum2)
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      v_it = _mm512_set1_epi32(k);
      v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dky = _mm512_mul_ps(v_dny,v_dky);
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      wp = 0.0;
      v_wp = _mm512_set1_pd(0.0);
/* add ky to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(12336),v_dky);
      v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(771),v_dky);
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       at1 = ci2*crealf(ffc[j+kk]); */
         v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,_mm512_int2mask(15),
                 (float *)&ffc[j+kk]);
         v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,_mm512_int2mask(15),
                 (float *)&ffc[j+kk+8]);
         v_at1 = _mm512_permute4f128_ps(v_at1,0);
         v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                 _mm512_int2mask(13260),(__m512i)v_at1,78);
         v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                 _mm512_int2mask(43690),(__m512i)v_at4,177);
         v_at1 = _mm512_mul_ps(v_ci2,v_at1);
/*       at2 = at1*dnx*(float) j; */
/*       at3 = dky*at1;           */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                  v_dkx);
         v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                  v_dkx);
/* normalize curl operators */
         v_at2 = _mm512_mul_ps(v_at1,v_dk1);
         v_at3 = _mm512_mul_ps(v_at1,v_dk2);
/*       at1 = at1*cimagf(ffc[j+kk]); */
         v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                    _mm512_int2mask(21845),(__m512i)v_at4,177);
         v_at1 = _mm512_mul_ps(v_at1,v_at4);
/*       zt1 = -cimagf(cu[2+4*(j+kj)])             */
/*            + crealf(cu[2+4*(j+kj)])*_Complex_I; */
/*       zt2 = -cimagf(cu[1+4*(j+kj)])             */
/*            + crealf(cu[1+4*(j+kj)])*_Complex_I; */
/*       zt3 = -cimagf(cu[4*(j+kj)])               */
/*            + crealf(cu[4*(j+kj)])*_Complex_I;   */
         v_zt3 = _mm512_load_ps((float *)&cu[4*(j+kj)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),v_zero,
                 v_zt3);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       bxyz[4*(j+kj)] = at3*zt1;             */
/*       bxyz[1+4*(j+kj)] = -at2*zt1;          */
/*       bxyz[2+4*(j+kj)] = at2*zt2 - at3*zt3; */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);
            v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                    v_zero);
         }
         _mm512_store_ps((float *)&bxyz[4*(j+kj)],v_zt1);
/*       wp += at1*(cu[4*(j+kj)]*conjf(cu[4*(j+kj)]) */
/*          + cu[1+4*(j+kj)]*conjf(cu[1+4*(j+kj)])   */
/*          + cu[2+4*(j+kj)]*conjf(cu[2+4*(j+kj)])); */
         v_zt4 = _mm512_mul_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                 _mm512_int2mask(16191),v_zt3,v_zt3));
/*       zt1 = -cimagf(cu[2+4*(j+k1)])             */
/*            + crealf(cu[2+4*(j+k1)])*_Complex_I; */
/*       zt2 = -cimagf(cu[1+4*(j+k1)])             */
/*            + crealf(cu[1+4*(j+k1)])*_Complex_I; */
/*       zt3 = -cimagf(cu[4*(j+k1)])               */
/*            + crealf(cu[4*(j+k1)])*_Complex_I;   */
         v_zt3 = _mm512_load_ps((float *)&cu[4*(j+k1)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),v_zero,
                 v_zt3);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
         v_zt1 = _mm512_mask_sub_ps(v_at2,_mm512_int2mask(12336),v_zero,
                 v_at2);
         v_zt2 = _mm512_mask_sub_ps(v_at3,_mm512_int2mask(771),v_zero,
                 v_at3);
/*       bxyz[4*(j+k1)] = -at3*zt1;            */
/*       bxyz[1+4*(j+k1)] = -at2*zt1;          */
/*       bxyz[2+4*(j+k1)] = at2*zt2 + at3*zt3; */
         v_zt1 = _mm512_mul_ps(v_zt1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_zt2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);
            v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                    v_zero);
         }
         _mm512_store_ps((float *)&bxyz[4*(j+k1)],v_zt1);
/*       wp += at1*(cu[4*(j+k1)]*conjf(cu[4*(j+k1)]) */
/*          + cu[1+4*(j+k1)]*conjf(cu[1+4*(j+k1)])   */
/*          + cu[2+4*(j+k1)]*conjf(cu[2+4*(j+k1)])); */
         v_zt4 = _mm512_fmadd_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                 _mm512_int2mask(16191),v_zt3,v_zt3),v_zt4);
/* convert to double precision before accumulating */
         v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt4));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt4,78));
         v_wp = _mm512_add_pd(v_wp,v_d);
/*       bxyz[4*(j+kj+l1)] = zero;   */
/*       bxyz[1+4*(j+kj+l1)] = zero; */
/*       bxyz[2+4*(j+kj+l1)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+kj+l1)],v_zero);
/*       bxyz[4*(j+k1+l1)] = zero;   */
/*       bxyz[1+4*(j+k1+l1)] = zero; */
/*       bxyz[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         at1 = ci2*crealf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         at1 = at1*cimagf(ffc[j+kk]);
         zt1 = -cimagf(cu[2+4*(j+kj)])
              + crealf(cu[2+4*(j+kj)])*_Complex_I;
         zt2 = -cimagf(cu[1+4*(j+kj)])
              + crealf(cu[1+4*(j+kj)])*_Complex_I;
         zt3 = -cimagf(cu[4*(j+kj)])
              + crealf(cu[4*(j+kj)])*_Complex_I;
         bxyz[4*(j+kj)] = at3*zt1;
         bxyz[1+4*(j+kj)] = -at2*zt1;
         bxyz[2+4*(j+kj)] = at2*zt2 - at3*zt3;
         zt1 = -cimagf(cu[2+4*(j+k1)])
              + crealf(cu[2+4*(j+k1)])*_Complex_I;
         zt2 = -cimagf(cu[1+4*(j+k1)])
              + crealf(cu[1+4*(j+k1)])*_Complex_I;
         zt3 = -cimagf(cu[4*(j+k1)])
              + crealf(cu[4*(j+k1)])*_Complex_I;
         bxyz[4*(j+k1)] = -at3*zt1;
         bxyz[1+4*(j+k1)] = -at2*zt1;
         bxyz[2+4*(j+k1)] = at2*zt2 + at3*zt3;
         bxyz[4*(j+kj+l1)] = zero;
         bxyz[1+4*(j+kj+l1)] = zero;
         bxyz[2+4*(j+kj+l1)] = zero;
         bxyz[4*(j+k1+l1)] = zero;
         bxyz[1+4*(j+k1+l1)] = zero;
         bxyz[2+4*(j+k1+l1)] = zero;
         at1 = at1*(cu[4*(j+kj)]*conjf(cu[4*(j+kj)])
             + cu[1+4*(j+kj)]*conjf(cu[1+4*(j+kj)])
             + cu[2+4*(j+kj)]*conjf(cu[2+4*(j+kj)])
             + cu[4*(j+k1)]*conjf(cu[4*(j+k1)])
             + cu[1+4*(j+k1)]*conjf(cu[1+4*(j+k1)])
             + cu[2+4*(j+k1)]*conjf(cu[2+4*(j+k1)]));
         wp += (double) at1;
      }
/*    sum2 += wp; */
      _mm512_store_pd(&dd[0],v_wp);
      for (j = 1; j < 8; j++) {
          dd[0] += dd[j];
      }
      sum2 += (wp + dd[0]);
   }
/* mode numbers kx = 0, nx/2 */
   wp = 0.0;
   v_wp = _mm512_setzero_pd();
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      at1 = ci2*crealf(ffc[kk]);
      at3 = at1*dny*(float) k;
      at1 = at1*cimagf(ffc[kk]);
      zt1 = -cimagf(cu[2+4*(kj)]) + crealf(cu[2+4*(kj)])*_Complex_I;
      zt3 = -cimagf(cu[4*(kj)]) + crealf(cu[4*(kj)])*_Complex_I;
      bxyz[4*kj] = at3*zt1;
      bxyz[1+4*kj] = zero;
      bxyz[2+4*kj] = -at3*zt3;
      bxyz[4*k1] = zero;
      bxyz[1+4*k1] = zero;
      bxyz[2+4*k1] = zero;
      bxyz[4*(kj+l1)] = zero;
      bxyz[1+4*(kj+l1)] = zero;
      bxyz[2+4*(kj+l1)] = zero;
      bxyz[4*(k1+l1)] = zero;
      bxyz[1+4*(k1+l1)] = zero;
      bxyz[2+4*(k1+l1)] = zero;
      at1 = at1*(cu[4*kj]*conjf(cu[4*kj])
          + cu[1+4*kj]*conjf(cu[1+4*kj])
          + cu[2+4*kj]*conjf(cu[2+4*kj]));
      wp += (double) at1;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
   for (j = 0; j < nxhs; j+=2) {
/*    at1 = ci2*crealf(ffc[j]); */
      v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,_mm512_int2mask(15),
              (float *)&ffc[j]);
      v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,_mm512_int2mask(15),
              (float *)&ffc[j+8]);
      v_at1 = _mm512_permute4f128_ps(v_at1,0);
      v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
              _mm512_int2mask(13260),(__m512i)v_at1,78);
      v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
               _mm512_int2mask(43690),(__m512i)v_at4,177);
      v_at1 = _mm512_mul_ps(v_ci2,v_at1);
/*    at2 = at1*dnx*(float) j; */
      v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
      v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(3084),v_dkx);
      v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(12336),v_dkx);
/* normalize curl operators */
      v_at2 = _mm512_mul_ps(v_at1,v_dk1);
      v_at3 = _mm512_mul_ps(v_at1,v_dk2);
/*    at1 = at1*cimagf(ffc[j]); */
      v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
              _mm512_int2mask(21845),(__m512i)v_at4,177);
      v_at1 = _mm512_mul_ps(v_at1,v_at4);
/*    zt1 = -cimagf(cu[2+4*j]) + crealf(cu[2+4*j])*_Complex_I; */
/*    zt2 = -cimagf(cu[1+4*j]) + crealf(cu[1+4*j])*_Complex_I; */
      v_zt3 = _mm512_load_ps((float *)&cu[4*j]);
      v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),v_zero,
              v_zt3);
      v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*    bxyz[4*j] = zero;       */
/*    bxyz[1+4*j] = -at2*zt1; */
/*    bxyz[2+4*j] = at2*zt2;  */
      v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
      v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
      v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
      v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
      v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
      if (j==0) {
         v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),v_zero);
         v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),v_zero);
      }
      _mm512_store_ps((float *)&bxyz[4*j],v_zt1);
/*    wp += at1*(cu[4*j]*conjf(cu[4*j]) */
/*       + cu[1+4*j]*conjf(cu[1+4*j])   */
/*       + cu[2+4*j]*conjf(cu[2+4*j])); */
      v_zt4 = _mm512_mul_ps(v_at1,_mm512_mask_mul_ps(v_zero,
              _mm512_int2mask(16191),v_zt3,v_zt3));
/* convert to double precision before accumulating */
      v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt4));
      v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt4,78));
      v_wp = _mm512_add_pd(v_wp,v_d);
/*    bxyz[4*(j+k1)] = zero;   */
/*    bxyz[1+4*(j+k1)] = zero; */
/*    bxyz[2+4*(j+k1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+k1)],v_zero);
/*    bxyz[4*(j+l1)] = zero;   */
/*    bxyz[1+4*(j+l1)] = zero; */
/*    bxyz[2+4*(j+l1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+l1)],v_zero);
/*    bxyz[4*(j+k1+l1)] = zero;    */
/*    bxyz[1+4*(j+k1+l1)] = zero; */
/*    bxyz[2+4*(j+k1+l1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
   }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
      at1 = ci2*crealf(ffc[j]);
      at2 = at1*dnx*(float) j;
      at1 = at1*cimagf(ffc[j]);
      zt1 = -cimagf(cu[2+4*j]) + crealf(cu[2+4*j])*_Complex_I;
      zt2 = -cimagf(cu[1+4*j]) + crealf(cu[1+4*j])*_Complex_I;
      bxyz[4*j] = zero;
      bxyz[1+4*j] = -at2*zt1;
      bxyz[2+4*j] = at2*zt2;
      bxyz[4*(j+k1)] = zero;
      bxyz[1+4*(j+k1)] = zero;
      bxyz[2+4*(j+k1)] = zero;
      bxyz[4*(j+l1)] = zero;
      bxyz[1+4*(j+l1)] = zero;
      bxyz[2+4*(j+l1)] = zero;
      bxyz[4*(j+k1+l1)] = zero;
      bxyz[1+4*(j+k1+l1)] = zero;
      bxyz[2+4*(j+k1+l1)] = zero;
      at1 = at1*(cu[4*j]*conjf(cu[4*j])
          + cu[1+4*j]*conjf(cu[1+4*j])
          + cu[2+4*j]*conjf(cu[2+4*j]));
      wp += (double) at1;
   }
   bxyz[0] = zero;
   bxyz[1] = zero;
   bxyz[2] = zero;
   bxyz[4*k1] = zero;
   bxyz[1+4*k1] = zero;
   bxyz[2+4*k1] = zero;
   bxyz[4*l1] = zero;
   bxyz[1+4*l1] = zero;
   bxyz[2+4*l1] = zero;
   bxyz[4*(k1+l1)] = zero;
   bxyz[1+4*(k1+l1)] = zero;
   bxyz[2+4*(k1+l1)] = zero;
/* sum2 += wp; */
   _mm512_store_pd(&dd[0],v_wp);
   for (j = 1; j < 8; j++) {
      dd[0] += dd[j];
   }
   sum2 += (wp + dd[0]);
/* *wm = wp*((float) nx)*((float) ny)*((float) nz); */
   *wm = (sum1 + sum2)*((float) nx)*((float) ny)*((float) nz);
   return;
}

/*--------------------------------------------------------------------*/
void ckncmmaxwel3(float complex exyz[], float complex bxyz[],
                  float complex cu[], float complex ffc[], float ci,
                  float dt, float *wf, float *wm, int nx, int ny,
                  int nz, int nxvh, int nyv, int nzv, int nxhd,
                  int nyhd, int nzhd) {
/* this subroutine solves 3d maxwell's equation in fourier space for
   transverse electric and magnetic fields with periodic boundary
   conditions.
   input: all, output: wf, wm, exyz, bxyz
   approximate flop count is:
   680*nxc*nyc*nzc + 149*(nxc*nyc + nxc*nzc + nyc*nzc)
   plus nxc*nyc*nzc divides
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   the magnetic field is first updated half a step using the equations:
   bx[kz][ky][kx] = bx[kz][ky][kx] - .5*dt*sqrt(-1)*
                    (ky*ez[kz][ky][kx]-kz*ey[kz][ky][kx])
   by[kz][ky][kx] = by[kz][ky][kx] - .5*dt*sqrt(-1)*
                   (kz*ex[kz][ky][kx]-kx*ez[kz][ky][kx])
   bz[kz][ky][kx] = bz[kz][ky][kx] - .5*dt*sqrt(-1)*
                   (kx*ey[kz][ky][kx]-ky*ex[kz][ky][kx])
   the electric field is then updated a whole step using the equations:
   ex[kz][ky][kx] = ex[kz][ky][kx] + c2*dt*sqrt(-1)
                    *(ky*bz[kz][ky][kx]-kz*by[kz][ky][kx])
                    - affp*dt*cux[kz][ky][kx]*s[kz][ky][kx]
   ey[kz][ky][kx] = ey[kz][ky][kx] + c2*dt*sqrt(-1)*
                    *(kz*bx[kz][ky][kx]-kx*bz[kz][ky][kx])
                    - affp*dt*cuy[kz][ky][kx]*s[kz][ky][kx]
   ez[kz][ky][kx] = ez[kz][ky][kx] + c2*dt*sqrt(-1)
                    *(kx*by[kz][ky][kx]-ky*bx[kz][ky][kx])
                    - affp*dt*cuz[kz][ky][kx]*s[kz][ky][kx]
   the magnetic field is finally updated the remaining half step with
   the new electric field and the previous magnetic field equations.
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, c2 = 1./(ci*ci)
   and s[kz][ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)
   j,k,l = fourier mode numbers, except for
   ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
   ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
   ex(kz=pi) = ey(kz=pi) = ez(kz=pi) = 0,
   ex(kx=0,ky=0,kz=0) = ey(kx=0,ky=0,kz=0) = ez(kx=0,ky=0,kz=0) = 0.
   and similarly for bx, by, bz.
   cu[l][k][j][i] = complex current density
   exyz[l][k][j][i] = complex transverse electric field
   bxyz[l][k][j][i] = complex magnetic field
   for component i, all for fourier mode (j1,k,l)
   real(ffc[0][0][0]) = affp = normalization constant = nx*ny*nz/np,
   where np=number of particles
   aimag(ffc[l][k][j]) = finite-size particle shape factor s,
   s[kz][ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
   for fourier mode (j,k,l)
   ci = reciprocal of velocity of light
   dt = time interval between successive calculations
   transverse electric field energy is also calculated, using
   wf = nx*ny*nz**sum((1/affp)*|exyz[kz][ky][kx]|**2)
   magnetic field energy is also calculated, using
   wm = nx*ny*nz**sum((c2/affp)*|bxyz[kz][ky][kx]|**2)
   nx/ny/nz = system length in x/y/z direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nzv = fourth dimension of field arrays, must be >= nz
   nxhd = second dimension of form factor array, must be >= nxh
   nyhd = third dimension of form factor array, must be >= nyh
   nzhd = fourth dimension of form factor array, must be >= nzh
   requires KNC, cu, exyz, bxyz, ffc need to be 64 byte aligned
   nxhd needs to be a multiple of 8
   nxvh needs to be a multiple of 2
   cu, exyz, bxyz needs to have 4 components
local data                                                 */
   int nxh, nyh, nzh, nxhs, itn, j, k, l, k1, l1, kk, kj, ll, lj;
   int nxyhd, nxvyh;
   float dnx, dny, dnz, dth, c2, cdt, affp, anorm, dkx, dky, dkz;
   float adt, afdt;
   float at1;
   float complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9;
   double wp, ws, sum1, sum2, sum3, sum4;
   __m512i v_j, v_it, v_n, v_m;
   __m512 v_dnx, v_dny, v_dnz, v_dkx, v_dky, v_dkz;
   __m512 v_zero, v_cdt, v_adt, v_afdt, v_dth, v_anorm;
   __m512 v_dk1, v_dk2, v_at2, v_at3;
   __m512 v_zt1, v_zt2, v_zt3, v_zt4, v_zt5, v_zt6, v_zt7;
   __m512d v_wp, v_ws, v_d;
   __attribute__((aligned(64))) double dd[8];
   if (ci <= 0.0)
      return;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxhs = 2*(nxh/2);
   itn = 1 > nxhs ? 1 : nxhs;
   nxyhd = nxhd*nyhd;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   dth = 0.5*dt;
   c2 = 1.0/(ci*ci);
   cdt = c2*dt;
   affp = creal(ffc[0]);
   adt = affp*dt;
   zero = 0.0 + 0.0*_Complex_I;
   anorm = 1.0/affp;
   v_j = _mm512_set_epi32(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0);
   v_n =  _mm512_set_epi32(15,14,11,10,9,8,13,12,7,6,3,2,1,0,5,4);
   v_m =  _mm512_set_epi32(15,14,9,8,13,12,11,10,7,6,1,0,5,4,3,2);
   v_dnx = _mm512_set1_ps(dnx);
   v_dny = _mm512_set1_ps(dny);
   v_dnz = _mm512_set1_ps(dnz);
   v_zero = _mm512_setzero_ps();
   v_cdt = _mm512_set1_ps(cdt);
   v_adt = _mm512_set1_ps(adt);
   v_dth = _mm512_set1_ps(dth);
   v_anorm = _mm512_set1_ps(anorm);
/* update electromagnetic field and sum field energies */
   sum1 = 0.0;
   sum2 = 0.0;
/* calculate the electromagnetic fields */
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
#pragma omp parallel
   {
#pragma omp for nowait \
private(j,k,l,k1,l1,ll,lj,kk,kj,dkz,dky,dkx,afdt,at1,zt1,zt2,zt3,zt4, \
zt5,zt6,zt7,zt8,zt9,ws,wp,v_it,v_dkx,v_dky,v_dkz,v_dk1,v_dk2,v_afdt, \
v_at2,v_at3,v_zt1,v_zt2,v_zt3,v_zt4,v_zt5,v_zt6,v_zt7,v_d,v_ws,v_wp,dd) \
reduction(+:sum1,sum2)
      for (l = 1; l < nzh; l++) {
         dkz = dnz*(float) l;
         v_dkz = _mm512_cvtfxpnt_round_adjustepi32_ps(_mm512_set1_epi32(l),
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkz = _mm512_mul_ps(v_dnz,v_dkz);
         ll = nxyhd*l;
         lj = nxvyh*l;
         l1 = nxvyh*nz - lj;
         ws = 0.0;
         wp = 0.0;
         v_ws = _mm512_set1_pd(0.0);
         v_wp = _mm512_set1_pd(0.0);
/* add kz to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(771),v_dkz);
         v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(3084),v_dkz);
         for (k = 1; k < nyh; k++) {
            dky = dny*(float) k;
            v_it = _mm512_set1_epi32(k);
            v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_dky = _mm512_mul_ps(v_dny,v_dky);
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
/* add ky to curl operators */
            v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(12336),
                     v_dky);
            v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(771),
                     v_dky);
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nxhs; j+=2) {
/*             dkx = dnx*(float) j; */
               v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
               v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                       _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
               v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
               v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                        v_dkx);
               v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                        v_dkx);
/*             afdt = adt*cimagf(ffc[j+kk+ll]); */
               v_afdt = _mm512_mask_loadunpacklo_ps(v_zero,
                        _mm512_int2mask(15),(float *)&ffc[j+kk+ll]);
               v_afdt = _mm512_mask_loadunpackhi_ps(v_afdt,
                        _mm512_int2mask(15),(float *)&ffc[j+kk+ll+8]);
               v_afdt = _mm512_permute4f128_ps(v_afdt,0);
               v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                       _mm512_int2mask(13260),(__m512i)v_afdt,78);
               v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                        _mm512_int2mask(21845),(__m512i)v_afdt,177);
               v_afdt = _mm512_mul_ps(v_adt,v_afdt);
/* update magnetic field half time step, ky > 0, kz > 0 */
/*             zt1 = -cimagf(exyz[2+4*(j+kj+lj)])             */
/*                  + crealf(exyz[2+4*(j+kj+lj)])*_Complex_I; */
/*             zt2 = -cimagf(exyz[1+4*(j+kj+lj)])             */
/*                  + crealf(exyz[1+4*(j+kj+lj)])*_Complex_I; */
/*             zt3 = -cimagf(exyz[4*(j+kj+lj)])               */
/*                  + crealf(exyz[4*(j+kj+lj)])*_Complex_I;   */
               v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+kj+lj)]);
               v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                       v_zero,v_zt4);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*             zt4 = bxyz[4*(j+kj+lj)] - dth*(dky*zt1 - dkz*zt2);   */
/*             zt5 = bxyz[1+4*(j+kj+lj)] - dth*(dkz*zt3 - dkx*zt1); */
/*             zt6 = bxyz[2+4*(j+kj+lj)] - dth*(dkx*zt2 - dky*zt3); */
               v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
               v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
               v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+kj+lj)]);
               v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*             zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*             zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*             zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
               v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),
                       v_zero,v_zt5);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*             zt7 = exyz[4*(j+kj+lj)] + cdt*(dky*zt1 - dkz*zt2)   */
/*                 - afdt*cu[4*(j+kj+lj)];                         */
/*             zt8 = exyz[1+4*(j+kj+lj)] + cdt*(dkz*zt3 - dkx*zt1) */
/*                 - afdt*cu[1+4*(j+kj+lj)];                       */
/*             zt9 = exyz[2+4*(j+kj+lj)] + cdt*(dkx*zt2 - dky*zt3) */
/*                  - afdt*cu[2+4*(j+kj+lj)];                       */
               v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
               v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                        v_zt4);
               v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj+lj)]);
               v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                       v_afdt,v_zt2);
               v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*             zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*             zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*             zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
               v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                       v_zero,v_zt4);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*             exyz[4*(j+kj+lj)] = zt7;   */
/*             exyz[1+4*(j+kj+lj)] = zt8; */
/*             exyz[2+4*(j+kj+lj)] = zt9; */
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                          v_zero);
                  _mm512_mask_store_ps((float *)&exyz[4*(j+kj+lj)],
                  _mm512_int2mask(65280),v_zt4);
               }
               else {
                  _mm512_store_ps((float *)&exyz[4*(j+kj+lj)],v_zt4);
               }
/*             ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) */
/*                + zt9*conjf(zt9));                        */
               v_zt6 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4));
/*             zt4 -= dth*(dky*zt1 - dkz*zt2); */
/*             zt5 -= dth*(dkz*zt3 - dkx*zt1); */
/*             zt6 -= dth*(dkx*zt2 - dky*zt3); */
               v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
               v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
               v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*             bxyz[4*(j+kj+lj)] = zt4;   */
/*             bxyz[1+4*(j+kj+lj)] = zt5; */
/*             bxyz[2+4*(j+kj+lj)] = zt6; */
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                          v_zero);
                  _mm512_mask_store_ps((float *)&bxyz[4*(j+kj+lj)],
                  _mm512_int2mask(65280),v_zt5);
               }
               else {
                  _mm512_store_ps((float *)&bxyz[4*(j+kj+lj)],v_zt5);
               }
/*             wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) */
/*                + zt6*conjf(zt6));                        */
               v_zt7 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5));
/* update magnetic field half time step, ky < 0, kz > 0 */
/*             zt1 = -cimagf(exyz[2+4*(j+k1+lj)])             */
/*                  + crealf(exyz[2+4*(j+k1+lj)])*_Complex_I; */
/*             zt2 = -cimagf(exyz[1+4*(j+k1+lj)])             */
/*                  + crealf(exyz[1+4*(j+k1+lj)])*_Complex_I; */
/*             zt3 = -cimagf(exyz[4*(j+k1+lj)])               */
/*                  + crealf(exyz[4*(j+k1+lj)])*_Complex_I;   */
               v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+k1+lj)]);
               v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                       v_zero,v_zt4);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
               v_at2 = _mm512_mask_sub_ps(v_dk1,_mm512_int2mask(12336),
                       v_zero,v_dk1);
               v_at3 = _mm512_mask_sub_ps(v_dk2,_mm512_int2mask(771),
                       v_zero,v_dk2);
/*             zt4 = bxyz[4*(j+k1+lj)] + dth*(dky*zt1 + dkz*zt2);   */
/*             zt5 = bxyz[1+4*(j+k1+lj)] - dth*(dkz*zt3 - dkx*zt1); */
/*             zt6 = bxyz[2+4*(j+k1+lj)] - dth*(dkx*zt2 + dky*zt3); */
               v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
               v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
               v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+k1+lj)]);
               v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*             zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*             zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*             zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
               v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),
                       v_zero,v_zt5);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*             zt7 = exyz[4*(j+k1+lj)] - cdt*(dky*zt1 + dkz*zt2)   */
/*                 - afdt*cu[4*(j+k1+lj)];                         */
/*             zt8 = exyz[1+4*(j+k1+lj)] + cdt*(dkz*zt3 - dkx*zt1) */
/*                 - afdt*cu[1+4*(j+k1+lj)];                       */
/*             zt9 = exyz[2+4*(j+k1+lj)] + cdt*(dkx*zt2 + dky*zt3) */
/*                 - afdt*cu[2+4*(j+k1+lj)];                       */
               v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
               v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                        v_zt4);
               v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1+lj)]);
               v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                       v_afdt,v_zt2);
               v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*             zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*             zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*             zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
               v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                       v_zero,v_zt4);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*             exyz[4*(j+k1+lj)] = zt7;   */
/*             exyz[1+4*(j+k1+lj)] = zt8; */
/*             exyz[2+4*(j+k1+lj)] = zt9; */
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                          v_zero);
                  _mm512_mask_store_ps((float *)&exyz[4*(j+k1+lj)],
                  _mm512_int2mask(65280),v_zt4);
              }
               else {
                  _mm512_store_ps((float *)&exyz[4*(j+k1+lj)],v_zt4);
               }
/*             ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) */
/*                + zt9*conjf(zt9));                        */
               v_zt6 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4),
                       v_zt6);
/*             zt4 += dth*(dky*zt1 + dkz*zt2); */
/*             zt5 -= dth*(dkz*zt3 - dkx*zt1); */
/*             zt6 -= dth*(dkx*zt2 + dky*zt3); */
               v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
               v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
               v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*             bxyz[4*(j+k1+lj)] = zt4;   */
/*             bxyz[1+4*(j+k1+lj)] = zt5; */
/*             bxyz[2+4*(j+k1+lj)] = zt6; */
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                          v_zero);
                  _mm512_mask_store_ps((float *)&bxyz[4*(j+k1+lj)],
                  _mm512_int2mask(65280),v_zt5);
               }
               else {
                  _mm512_store_ps((float *)&bxyz[4*(j+k1+lj)],v_zt5);
               }
/*             wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) */
/*                        + zt6*conjf(zt6));                */
               v_zt7 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5),
                       v_zt7);
/* update magnetic field half time step, ky > 0, kz < 0 */
/*             zt1 = -cimagf(exyz[2+4*(j+kj+l1)])             */
/*                  + crealf(exyz[2+4*(j+kj+l1)])*_Complex_I; */
/*             zt2 = -cimagf(exyz[1+4*(j+kj+l1)])             */
/*                  + crealf(exyz[1+4*(j+kj+l1)])*_Complex_I; */
/*             zt3 = -cimagf(exyz[4*(j+kj+l1)])               */
/*                  + crealf(exyz[4*(j+kj+l1)])*_Complex_I;   */
               v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+kj+l1)]);
               v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                       v_zero,v_zt4);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
               v_at2 = _mm512_mask_sub_ps(v_dk1,_mm512_int2mask(771),
                       v_zero,v_dk1);
               v_at3 = _mm512_mask_sub_ps(v_dk2,_mm512_int2mask(3084),
                       v_zero,v_dk2);
/*             zt4 = bxyz[4*(j+kj+l1)] - dth*(dky*zt1 + dkz*zt2);   */
/*             zt5 = bxyz[1+4*(j+kj+l1)] + dth*(dkz*zt3 + dkx*zt1); */
/*             zt6 = bxyz[2+4*(j+kj+l1)] - dth*(dkx*zt2 - dky*zt3); */
               v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
               v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
               v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+kj+l1)]);
               v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*             zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*             zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*             zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
               v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),
                       v_zero,v_zt5);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*             zt7 = exyz[4*(j+kj+l1)] + cdt*(dky*zt1 + dkz*zt2)   */
/*                 - afdt*cu[4*(j+kj+l1)];                         */
/*             zt8 = exyz[1+4*(j+kj+l1)] - cdt*(dkz*zt3 + dkx*zt1) */
/*                 - afdt*cu[1+4*(j+kj+l1)];                       */
/*             zt9 = exyz[2+4*(j+kj+l1)] + cdt*(dkx*zt2 - dky*zt3) */
/*                 - afdt*cu[2+4*(j+kj+l1)];                       */
               v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
               v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                        v_zt4);
               v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj+l1)]);
               v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                       v_afdt,v_zt2);
               v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*             zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*             zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*             zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
               v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                       v_zero,v_zt4);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*             exyz[4*(j+kj+l1)] = zt7;   */
/*             exyz[1+4*(j+kj+l1)] = zt8; */
/*             exyz[2+4*(j+kj+l1)] = zt9; */
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                          v_zero);
                  _mm512_mask_store_ps((float *)&exyz[4*(j+kj+l1)],
                  _mm512_int2mask(65280),v_zt4);
               }
               else {
                  _mm512_store_ps((float *)&exyz[4*(j+kj+l1)],v_zt4);
               }
/*             ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) */
/*                        + zt9*conjf(zt9));                */
               v_zt6 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4),
                       v_zt6);
/*             zt4 -= dth*(dky*zt1 + dkz*zt2); */
/*             zt5 += dth*(dkz*zt3 + dkx*zt1); */
/*             zt6 -= dth*(dkx*zt2 - dky*zt3); */
               v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
               v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
               v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*             bxyz[4*(j+kj+l1)] = zt4;   */
/*             bxyz[1+4*(j+kj+l1)] = zt5; */
/*             bxyz[2+4*(j+kj+l1)] = zt6; */
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                          v_zero);
                  _mm512_mask_store_ps((float *)&bxyz[4*(j+kj+l1)],
                  _mm512_int2mask(65280),v_zt5);
               }
               else {
                  _mm512_store_ps((float *)&bxyz[4*(j+kj+l1)],v_zt5);
               }
/*             wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) */
/*                        + zt6*conjf(zt6));                */
               v_zt7 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5),
                       v_zt7);
/* update magnetic field half time step, ky < 0, kz < 0 */
/*             zt1 = -cimagf(exyz[2+4*(j+k1+l1)])             */
/*                  + crealf(exyz[2+4*(j+k1+l1)])*_Complex_I; */
/*             zt2 = -cimagf(exyz[1+4*(j+k1+l1)])             */
/*                  + crealf(exyz[1+4*(j+k1+l1)])*_Complex_I; */
/*             zt3 = -cimagf(exyz[4*(j+k1+l1)])              */
/*                  + crealf(exyz[4*(j+k1+l1)])*_Complex_I;  */
               v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
               v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                       v_zero,v_zt4);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
               v_at2 = _mm512_mask_sub_ps(v_dk1,_mm512_int2mask(13107),
                       v_zero,v_dk1);
               v_at3 = _mm512_mask_sub_ps(v_dk2,_mm512_int2mask(3855),
                       v_zero,v_dk2);
/*             zt4 = bxyz[4*(j+k1+l1)] + dth*(dky*zt1 - dkz*zt2);   */
/*             zt5 = bxyz[1+4*(j+k1+l1)] + dth*(dkz*zt3 + dkx*zt1); */
/*             zt6 = bxyz[2+4*(j+k1+l1)] - dth*(dkx*zt2 + dky*zt3); */
               v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
               v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
               v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+k1+l1)]);
               v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*             zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*             zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*             zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
               v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),
                       v_zero,v_zt5);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*             zt7 = exyz[4*(j+k1+l1)] - cdt*(dky*zt1 - dkz*zt2)   */
/*                 - afdt*cu[4*(j+k1+l1)];                         */
/*             zt8 = exyz[1+4*(j+k1+l1)] - cdt*(dkz*zt3 + dkx*zt1) */
/*                 - afdt*cu[1+4*(j+k1+l1)];                       */
/*             zt9 = exyz[2+4*(j+k1+l1)] + cdt*(dkx*zt2 + dky*zt3) */
/*                 - afdt*cu[2+4*(j+k1+l1)];                       */
               v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
               v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                        v_zt4);
               v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1+l1)]);
               v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                       v_afdt,v_zt2);
               v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*             zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*             zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*             zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
               v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                       v_zero,v_zt4);
               v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*             exyz[4*(j+k1+l1)] = zt7;   */
/*             exyz[1+4*(j+k1+l1)] = zt8; */
/*             exyz[2+4*(j+k1+l1)] = zt9; */
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                          v_zero);
                  _mm512_mask_store_ps((float *)&exyz[4*(j+k1+l1)],
                  _mm512_int2mask(65280),v_zt4);
               }
               else {
                  _mm512_store_ps((float *)&exyz[4*(j+k1+l1)],v_zt4);
               }
/*             ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) */
/*                        + zt9*conjf(zt9));                */
               v_zt6 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4),
                       v_zt6);
/*             zt4 += dth*(dky*zt1 - dkz*zt2); */
/*             zt5 += dth*(dkz*zt3 + dkx*zt1); */
/*             zt6 -= dth*(dkx*zt2 + dky*zt3); */
               v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
               v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
               v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
               v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
               v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
               v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*             bxyz[4*(j+k1+l1)] = zt4;   */
/*             bxyz[1+4*(j+k1+l1)] = zt5; */
/*             bxyz[2+4*(j+k1+l1)] = zt6; */
/* zero out kx = 0 mode */
               if (j==0) {
                  v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                          v_zero);
                  _mm512_mask_store_ps((float *)&bxyz[4*(j+k1+l1)],
                  _mm512_int2mask(65280),v_zt5);
               }
               else {
                  _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zt5);
               }
/*             wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) */
/*                        + zt6*conjf(zt6));                */
               v_zt7 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5),
                       v_zt7);
/* convert to double precision before accumulating */
               v_ws = _mm512_add_pd(v_ws,_mm512_cvtpslo_pd(v_zt6));
               v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt6,78));
               v_ws = _mm512_add_pd(v_ws,v_d);
               v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt7));
               v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt7,78));
               v_wp = _mm512_add_pd(v_wp,v_d);
            }
/* loop over remaining elements */
            for (j = itn; j < nxh; j++) {
               dkx = dnx*(float) j;
               afdt = adt*cimagf(ffc[j+kk+ll]);
/* update magnetic field half time step, ky > 0, kz > 0 */
               zt1 = -cimagf(exyz[2+4*(j+kj+lj)])
                    + crealf(exyz[2+4*(j+kj+lj)])*_Complex_I;
               zt2 = -cimagf(exyz[1+4*(j+kj+lj)])
                    + crealf(exyz[1+4*(j+kj+lj)])*_Complex_I;
               zt3 = -cimagf(exyz[4*(j+kj+lj)])
                    + crealf(exyz[4*(j+kj+lj)])*_Complex_I;
               zt4 = bxyz[4*(j+kj+lj)] - dth*(dky*zt1 - dkz*zt2);
               zt5 = bxyz[1+4*(j+kj+lj)] - dth*(dkz*zt3 - dkx*zt1);
               zt6 = bxyz[2+4*(j+kj+lj)] - dth*(dkx*zt2 - dky*zt3);
/* update electric field whole time step */
               zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
               zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
               zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
               zt7 = exyz[4*(j+kj+lj)] + cdt*(dky*zt1 - dkz*zt2)
                   - afdt*cu[4*(j+kj+lj)];
               zt8 = exyz[1+4*(j+kj+lj)] + cdt*(dkz*zt3 - dkx*zt1)
                   - afdt*cu[1+4*(j+kj+lj)];
               zt9 = exyz[2+4*(j+kj+lj)] + cdt*(dkx*zt2 - dky*zt3)
                   - afdt*cu[2+4*(j+kj+lj)];
/* update magnetic field half time step and store electric field */
               zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
               zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
               zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
               exyz[4*(j+kj+lj)] = zt7;
               exyz[1+4*(j+kj+lj)] = zt8;
               exyz[2+4*(j+kj+lj)] = zt9;
               at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8)
                          + zt9*conjf(zt9));
               ws += (double) at1;
               zt4 -= dth*(dky*zt1 - dkz*zt2);
               zt5 -= dth*(dkz*zt3 - dkx*zt1);
               zt6 -= dth*(dkx*zt2 - dky*zt3);
               bxyz[4*(j+kj+lj)] = zt4;
               bxyz[1+4*(j+kj+lj)] = zt5;
               bxyz[2+4*(j+kj+lj)] = zt6;
               at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
                          + zt6*conjf(zt6));
               wp += (double) at1;
/* update magnetic field half time step, ky < 0, kz > 0 */
               zt1 = -cimagf(exyz[2+4*(j+k1+lj)])
                    + crealf(exyz[2+4*(j+k1+lj)])*_Complex_I;
               zt2 = -cimagf(exyz[1+4*(j+k1+lj)])
                    + crealf(exyz[1+4*(j+k1+lj)])*_Complex_I;
               zt3 = -cimagf(exyz[4*(j+k1+lj)])
                    + crealf(exyz[4*(j+k1+lj)])*_Complex_I;
               zt4 = bxyz[4*(j+k1+lj)] + dth*(dky*zt1 + dkz*zt2);
               zt5 = bxyz[1+4*(j+k1+lj)] - dth*(dkz*zt3 - dkx*zt1);
               zt6 = bxyz[2+4*(j+k1+lj)] - dth*(dkx*zt2 + dky*zt3);
/* update electric field whole time step */
               zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
               zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
               zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
               zt7 = exyz[4*(j+k1+lj)] - cdt*(dky*zt1 + dkz*zt2)
                   - afdt*cu[4*(j+k1+lj)];
               zt8 = exyz[1+4*(j+k1+lj)] + cdt*(dkz*zt3 - dkx*zt1)
                   - afdt*cu[1+4*(j+k1+lj)];
               zt9 = exyz[2+4*(j+k1+lj)] + cdt*(dkx*zt2 + dky*zt3)
                   - afdt*cu[2+4*(j+k1+lj)];
/* update magnetic field half time step and store electric field */
               zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
               zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
               zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
               exyz[4*(j+k1+lj)] = zt7;
               exyz[1+4*(j+k1+lj)] = zt8;
               exyz[2+4*(j+k1+lj)] = zt9;
               at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8)
                          + zt9*conjf(zt9));
               ws += (double) at1;
               zt4 += dth*(dky*zt1 + dkz*zt2);
               zt5 -= dth*(dkz*zt3 - dkx*zt1);
               zt6 -= dth*(dkx*zt2 + dky*zt3);
               bxyz[4*(j+k1+lj)] = zt4;
               bxyz[1+4*(j+k1+lj)] = zt5;
               bxyz[2+4*(j+k1+lj)] = zt6;
               at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
                         + zt6*conjf(zt6));
               wp += (double) at1;
/* update magnetic field half time step, ky > 0, kz < 0 */
               zt1 = -cimagf(exyz[2+4*(j+kj+l1)])
                    + crealf(exyz[2+4*(j+kj+l1)])*_Complex_I;
               zt2 = -cimagf(exyz[1+4*(j+kj+l1)])
                    + crealf(exyz[1+4*(j+kj+l1)])*_Complex_I;
               zt3 = -cimagf(exyz[4*(j+kj+l1)])
                    + crealf(exyz[4*(j+kj+l1)])*_Complex_I;
               zt4 = bxyz[4*(j+kj+l1)] - dth*(dky*zt1 + dkz*zt2);
               zt5 = bxyz[1+4*(j+kj+l1)] + dth*(dkz*zt3 + dkx*zt1);
               zt6 = bxyz[2+4*(j+kj+l1)] - dth*(dkx*zt2 - dky*zt3);
/* update electric field whole time step */
               zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
               zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
               zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
               zt7 = exyz[4*(j+kj+l1)] + cdt*(dky*zt1 + dkz*zt2)
                   - afdt*cu[4*(j+kj+l1)];
               zt8 = exyz[1+4*(j+kj+l1)] - cdt*(dkz*zt3 + dkx*zt1)
                   - afdt*cu[1+4*(j+kj+l1)];
               zt9 = exyz[2+4*(j+kj+l1)] + cdt*(dkx*zt2 - dky*zt3)
                   - afdt*cu[2+4*(j+kj+l1)];
/* update magnetic field half time step and store electric field */
               zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
               zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
               zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
               exyz[4*(j+kj+l1)] = zt7;
               exyz[1+4*(j+kj+l1)] = zt8;
               exyz[2+4*(j+kj+l1)] = zt9;
               at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8)
                          + zt9*conjf(zt9));
               ws += (double) at1;
               zt4 -= dth*(dky*zt1 + dkz*zt2);
               zt5 += dth*(dkz*zt3 + dkx*zt1);
               zt6 -= dth*(dkx*zt2 - dky*zt3);
               bxyz[4*(j+kj+l1)] = zt4;
               bxyz[1+4*(j+kj+l1)] = zt5;
               bxyz[2+4*(j+kj+l1)] = zt6;
               at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
                          + zt6*conjf(zt6));
               wp += (double) at1;
/* update magnetic field half time step, ky < 0, kz < 0 */
               zt1 = -cimagf(exyz[2+4*(j+k1+l1)])
                    + crealf(exyz[2+4*(j+k1+l1)])*_Complex_I;
               zt2 = -cimagf(exyz[1+4*(j+k1+l1)])
                    + crealf(exyz[1+4*(j+k1+l1)])*_Complex_I;
               zt3 = -cimagf(exyz[4*(j+k1+l1)])
                    + crealf(exyz[4*(j+k1+l1)])*_Complex_I;
               zt4 = bxyz[4*(j+k1+l1)] + dth*(dky*zt1 - dkz*zt2);
               zt5 = bxyz[1+4*(j+k1+l1)] + dth*(dkz*zt3 + dkx*zt1);
               zt6 = bxyz[2+4*(j+k1+l1)] - dth*(dkx*zt2 + dky*zt3);
/* update electric field whole time step */
               zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
               zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
               zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
               zt7 = exyz[4*(j+k1+l1)] - cdt*(dky*zt1 - dkz*zt2)
                   - afdt*cu[4*(j+k1+l1)];
               zt8 = exyz[1+4*(j+k1+l1)] - cdt*(dkz*zt3 + dkx*zt1)
                   - afdt*cu[1+4*(j+k1+l1)];
               zt9 = exyz[2+4*(j+k1+l1)] + cdt*(dkx*zt2 + dky*zt3)
                   - afdt*cu[2+4*(j+k1+l1)];
/* update magnetic field half time step and store electric field */
               zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
               zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
               zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
               exyz[4*(j+k1+l1)] = zt7;
               exyz[1+4*(j+k1+l1)] = zt8;
               exyz[2+4*(j+k1+l1)] = zt9;
               at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8)
                          + zt9*conjf(zt9));
               ws += (double) at1;
               zt4 += dth*(dky*zt1 - dkz*zt2);
               zt5 += dth*(dkz*zt3 + dkx*zt1);
               zt6 -= dth*(dkx*zt2 + dky*zt3);
               bxyz[4*(j+k1+l1)] = zt4;
               bxyz[1+4*(j+k1+l1)] = zt5;
               bxyz[2+4*(j+k1+l1)] = zt6;
               at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
                          + zt6*conjf(zt6));
               wp += (double) at1;
            }
         }
/* mode numbers kx = 0, nx/2 */
         for (k = 1; k < nyh; k++) {
            dky = dny*(float) k;
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
            afdt = adt*cimagf(ffc[kk+ll]);
/* update magnetic field half time step, kz > 0 */
            zt1 = -cimagf(exyz[2+4*(kj+lj)])
                 + crealf(exyz[2+4*(kj+lj)])*_Complex_I;
            zt2 = -cimagf(exyz[1+4*(kj+lj)])
                 + crealf(exyz[1+4*(kj+lj)])*_Complex_I;
            zt3 = -cimagf(exyz[4*(kj+lj)])
                 + crealf(exyz[4*(kj+lj)])*_Complex_I;
            zt4 = bxyz[4*(kj+lj)] - dth*(dky*zt1 - dkz*zt2);
            zt5 = bxyz[1+4*(kj+lj)] - dth*(dkz*zt3);
            zt6 = bxyz[2+4*(kj+lj)] + dth*(dky*zt3);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exyz[4*(kj+lj)] + cdt*(dky*zt1 - dkz*zt2)
                - afdt*cu[4*(kj+lj)];
            zt8 = exyz[1+4*(kj+lj)] + cdt*(dkz*zt3) - afdt*cu[1+4*(kj+lj)];
            zt9 = exyz[2+4*(kj+lj)] - cdt*(dky*zt3) - afdt*cu[2+4*(kj+lj)];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exyz[4*(kj+lj)] = zt7;
            exyz[1+4*(kj+lj)] = zt8;
            exyz[2+4*(kj+lj)] = zt9;
            at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
            ws += (double) at1;
            zt4 -= dth*(dky*zt1 - dkz*zt2);
            zt5 -= dth*(dkz*zt3);
            zt6 += dth*(dky*zt3);
            bxyz[4*(kj+lj)] = zt4;
            bxyz[1+4*(kj+lj)] = zt5;
            bxyz[2+4*(kj+lj)] = zt6;
            at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
            wp += (double) at1;
            bxyz[4*(k1+lj)] = zero;
            bxyz[1+4*(k1+lj)] = zero;
            bxyz[2+4*(k1+lj)] = zero;
            exyz[4*(k1+lj)] = zero;
            exyz[1+4*(k1+lj)] = zero;
            exyz[2+4*(k1+lj)] = zero;
/* update magnetic field half time step, kz < 0 */
            zt1 = -cimagf(exyz[2+4*(kj+l1)])
                 + crealf(exyz[2+4*(kj+l1)])*_Complex_I;
            zt2 = -cimagf(exyz[1+4*(kj+l1)])
                 + crealf(exyz[1+4*(kj+l1)])*_Complex_I;
            zt3 = -cimagf(exyz[4*(kj+l1)])
                 + crealf(exyz[4*(kj+l1)])*_Complex_I;
            zt4 = bxyz[4*(kj+l1)] - dth*(dky*zt1 + dkz*zt2);
            zt5 = bxyz[1+4*(kj+l1)] + dth*(dkz*zt3);
            zt6 = bxyz[2+4*(kj+l1)] + dth*(dky*zt3);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exyz[4*(kj+l1)] + cdt*(dky*zt1 + dkz*zt2)
                - afdt*cu[4*(kj+l1)];
            zt8 = exyz[1+4*(kj+l1)] - cdt*(dkz*zt3) - afdt*cu[1+4*(kj+l1)];
            zt9 = exyz[2+4*(kj+l1)] - cdt*(dky*zt3) - afdt*cu[2+4*(kj+l1)];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exyz[4*(kj+l1)] = zt7;
            exyz[1+4*(kj+l1)] = zt8;
            exyz[2+4*(kj+l1)] = zt9;
            at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
            ws += (double) at1;
            zt4 -= dth*(dky*zt1 + dkz*zt2);
            zt5 += dth*(dkz*zt3);
            zt6 += dth*(dky*zt3);
            bxyz[4*(kj+l1)] = zt4;
            bxyz[1+4*(kj+l1)] = zt5;
            bxyz[2+4*(kj+l1)] = zt6;
            at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
            wp += (double) at1;
            bxyz[4*(k1+l1)] = zero;
            bxyz[1+4*(k1+l1)] = zero;
            bxyz[2+4*(k1+l1)] = zero;
            exyz[4*(k1+l1)] = zero;
            exyz[1+4*(k1+l1)] = zero;
            exyz[2+4*(k1+l1)] = zero;
         }
/* mode numbers ky = 0, ny/2 */
         k1 = nxvh*nyh;
/* add ky to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(12336),v_zero);
         v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(771),v_zero);
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          dkx = dnx*(float) j; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
            v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                     v_dkx);
            v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                     v_dkx);
/*          afdt = adt*cimagf(ffc[j+ll]); */
            v_afdt = _mm512_mask_loadunpacklo_ps(v_zero,
                     _mm512_int2mask(15),(float *)&ffc[j+ll]);
            v_afdt = _mm512_mask_loadunpackhi_ps(v_afdt,
                     _mm512_int2mask(15),(float *)&ffc[j+ll+8]);
            v_afdt = _mm512_permute4f128_ps(v_afdt,0);
            v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                     _mm512_int2mask(13260),(__m512i)v_afdt,78);
            v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                      _mm512_int2mask(21845),(__m512i)v_afdt,177);
            v_afdt = _mm512_mul_ps(v_adt,v_afdt);
/* update magnetic field half time step, kz > 0 */
/*          zt1 = -cimagf(exyz[2+4*(j+lj)])             */
/*               + crealf(exyz[2+4*(j+lj)])*_Complex_I; */
/*          zt2 = -cimagf(exyz[1+4*(j+lj)])             */
/*               + crealf(exyz[1+4*(j+lj)])*_Complex_I; */
/*          zt3 = -cimagf(exyz[4*(j+lj)])               */
/*               + crealf(exyz[4*(j+lj)])*_Complex_I;   */
            v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+lj)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                    v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          zt4 = bxyz[4*(j+lj)] + dth*(dkz*zt2);             */
/*          zt5 = bxyz[1+4*(j+lj)] - dth*(dkz*zt3 - dkx*zt1); */
/*          zt6 = bxyz[2+4*(j+lj)] - dth*(dkx*zt2);           */
            v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
            v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+lj)]);
            v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*          zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*          zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*          zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),v_zero,
                    v_zt5);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          zt7 = exyz[4*(j+lj)] - cdt*(dkz*zt2) - afdt*cu[4*(j+lj)];     */
/*          zt8 = exyz[1+4*(j+lj)] + cdt*(dkz*zt3 - dkx*zt1)              */
/*              - afdt*cu[1+4*(j+lj)];                                    */
/*          zt9 = exyz[2+4*(j+lj)] + cdt*(dkx*zt2) - afdt*cu[2+4*(j+lj)]; */
            v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
            v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                    v_zt4);
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+lj)]);
            v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                    v_afdt,v_zt2);
            v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*          zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*          zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*          zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                    v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          exyz[4*(j+lj)] = zt7;   */
/*          exyz[1+4*(j+lj)] = zt8; */
/*          exyz[2+4*(j+lj)] = zt9; */
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                        v_zero);
               _mm512_mask_store_ps((float *)&exyz[4*(j+lj)],
               _mm512_int2mask(65280),v_zt4);
            }
            else {
               _mm512_store_ps((float *)&exyz[4*(j+lj)],v_zt4);
            }
/*          ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
            v_zt6 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4));
/*          zt4 += dth*(dkz*zt2);           */
/*          zt5 -= dth*(dkz*zt3 - dkx*zt1); */
/*          zt6 -= dth*(dkx*zt2);           */
            v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
            v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*          bxyz[4*(j+lj)] = zt4;   */
/*          bxyz[1+4*(j+lj)] = zt5; */
/*          bxyz[2+4*(j+lj)] = zt6; */
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                       v_zero);
               _mm512_mask_store_ps((float *)&bxyz[4*(j+lj)],
               _mm512_int2mask(65280),v_zt5);
            }
            else {
               _mm512_store_ps((float *)&bxyz[4*(j+lj)],v_zt5);
            }
/*          wp +=  anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
            v_zt7 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5));
/*          bxyz[4*(j+k1+lj)] = zero;   */
/*          bxyz[1+4*(j+k1+lj)] = zero; */
/*          bxyz[2+4*(j+k1+lj)] = zero; */
            _mm512_store_ps((float *)&bxyz[4*(j+k1+lj)],v_zero);
/*          exyz[4*(j+k1+lj)] = zero;   */
/*          exyz[1+4*(j+k1+lj)] = zero; */
/*          exyz[2+4*(j+k1+lj)] = zero; */
            _mm512_store_ps((float *)&exyz[4*(j+k1+lj)],v_zero);
/* update magnetic field half time step, kz > 0 */
/*          zt1 = -cimagf(exyz[2+4*(j+l1)])             */
/*               + crealf(exyz[2+4*(j+l1)])*_Complex_I; */
/*          zt2 = -cimagf(exyz[1+4*(j+l1)])             */
/*               + crealf(exyz[1+4*(j+l1)])*_Complex_I; */
/*          zt3 = -cimagf(exyz[4*(j+l1)])               */
/*               + crealf(exyz[4*(j+l1)])*_Complex_I;   */
            v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+l1)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                     v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
            v_at2 = _mm512_mask_sub_ps(v_dk1,_mm512_int2mask(771),v_zero,
                     v_dk1);
            v_at3 = _mm512_mask_sub_ps(v_dk2,_mm512_int2mask(3084),v_zero,
                    v_dk2);
/*          zt4 = bxyz[4*(j+l1)] - dth*(dkz*zt2);             */
/*          zt5 = bxyz[1+4*(j+l1)] + dth*(dkz*zt3 + dkx*zt1); */
/*          zt6 = bxyz[2+4*(j+l1)] - dth*(dkx*zt2);           */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+l1)]);
            v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*          zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*          zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*          zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),v_zero,
                    v_zt5);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          zt7 = exyz[4*(j+l1)] + cdt*(dkz*zt2) - afdt*cu[4*(j+l1)];     */
/*          zt8 = exyz[1+4*(j+l1)] - cdt*(dkz*zt3 + dkx*zt1)              */
/*              - afdt*cu[1+4*(j+l1)];                                    */
/*          zt9 = exyz[2+4*(j+l1)] + cdt*(dkx*zt2) - afdt*cu[2+4*(j+l1)]; */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                    v_zt4);
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+l1)]);
            v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                     v_afdt,v_zt2);
            v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*          zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*          zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*          zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                    v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          exyz[4*(j+l1)] = zt7;   */
/*          exyz[1+4*(j+l1)] = zt8; */
/*          exyz[2+4*(j+l1)] = zt9; */
/* zero out kx = 0 mode */
            if (j==0) {
                v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                        v_zero);
               _mm512_mask_store_ps((float *)&exyz[4*(j+l1)],
               _mm512_int2mask(65280),v_zt4);
            }
            else {
               _mm512_store_ps((float *)&exyz[4*(j+l1)],v_zt4);
            }
/*          ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
            v_zt6 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4),
                    v_zt6);
/*          zt4 -= dth*(dkz*zt2);           */
/*          zt5 += dth*(dkz*zt3 + dkx*zt1); */
/*          zt6 -= dth*(dkx*zt2);           */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*          bxyz[4*(j+l1)] = zt4;   */
/*          bxyz[1+4*(j+l1)] = zt5; */
/*          bxyz[2+4*(j+l1)] = zt6; */
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                       v_zero);
               _mm512_mask_store_ps((float *)&bxyz[4*(j+l1)],
               _mm512_int2mask(65280),v_zt5);
            }
            else {
               _mm512_store_ps((float *)&bxyz[4*(j+l1)],v_zt5);
            }
/*          wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
            v_zt7 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5),
                     v_zt7);
/* convert to double precision before accumulating */
            v_ws = _mm512_add_pd(v_ws,_mm512_cvtpslo_pd(v_zt6));
            v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt6,78));
            v_ws = _mm512_add_pd(v_ws,v_d);
            v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt7));
            v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt7,78));
            v_wp = _mm512_add_pd(v_wp,v_d);
/*          bxyz[4*(j+k1+l1)] = zero;   */
/*          bxyz[1+4*(j+k1+l1)] = zero; */
/*          bxyz[2+4*(j+k1+l1)] = zero; */
            _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
/*          exyz[4*(j+k1+l1)] = zero;   */
/*          exyz[1+4*(j+k1+l1)] = zero; */
/*          exyz[2+4*(j+k1+l1)] = zero; */
            _mm512_store_ps((float *)&exyz[4*(j+k1+l1)],v_zero);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            dkx = dnx*(float) j;  
            afdt = adt*cimagf(ffc[j+ll]);
/* update magnetic field half time step, kz > 0 */
            zt1 = -cimagf(exyz[2+4*(j+lj)])
                 + crealf(exyz[2+4*(j+lj)])*_Complex_I;
            zt2 = -cimagf(exyz[1+4*(j+lj)])
                 + crealf(exyz[1+4*(j+lj)])*_Complex_I;
            zt3 = -cimagf(exyz[4*(j+lj)])
                 + crealf(exyz[4*(j+lj)])*_Complex_I;
            zt4 = bxyz[4*(j+lj)] + dth*(dkz*zt2);
            zt5 = bxyz[1+4*(j+lj)] - dth*(dkz*zt3 - dkx*zt1);
            zt6 = bxyz[2+4*(j+lj)] - dth*(dkx*zt2);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exyz[4*(j+lj)] - cdt*(dkz*zt2) - afdt*cu[4*(j+lj)];
            zt8 = exyz[1+4*(j+lj)] + cdt*(dkz*zt3 - dkx*zt1)
                - afdt*cu[1+4*(j+lj)];
            zt9 = exyz[2+4*(j+lj)] + cdt*(dkx*zt2) - afdt*cu[2+4*(j+lj)];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exyz[4*(j+lj)] = zt7;
            exyz[1+4*(j+lj)] = zt8;
            exyz[2+4*(j+lj)] = zt9;
            at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
            ws += (double) at1;
            zt4 += dth*(dkz*zt2);
            zt5 -= dth*(dkz*zt3 - dkx*zt1);
            zt6 -= dth*(dkx*zt2);
            bxyz[4*(j+lj)] = zt4;
            bxyz[1+4*(j+lj)] = zt5;
            bxyz[2+4*(j+lj)] = zt6;
            at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
            wp += (double) at1;
            bxyz[4*(j+k1+lj)] = zero;
            bxyz[1+4*(j+k1+lj)] = zero;
            bxyz[2+4*(j+k1+lj)] = zero;
            exyz[4*(j+k1+lj)] = zero;
            exyz[1+4*(j+k1+lj)] = zero;
            exyz[2+4*(j+k1+lj)] = zero;
/* update magnetic field half time step, kz > 0 */
            zt1 = -cimagf(exyz[2+4*(j+l1)])
                 + crealf(exyz[2+4*(j+l1)])*_Complex_I;
            zt2 = -cimagf(exyz[1+4*(j+l1)])
                 + crealf(exyz[1+4*(j+l1)])*_Complex_I;
            zt3 = -cimagf(exyz[4*(j+l1)])
                 + crealf(exyz[4*(j+l1)])*_Complex_I;
            zt4 = bxyz[4*(j+l1)] - dth*(dkz*zt2);
            zt5 = bxyz[1+4*(j+l1)] + dth*(dkz*zt3 + dkx*zt1);
            zt6 = bxyz[2+4*(j+l1)] - dth*(dkx*zt2);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exyz[4*(j+l1)] + cdt*(dkz*zt2) - afdt*cu[4*(j+l1)];
            zt8 = exyz[1+4*(j+l1)] - cdt*(dkz*zt3 + dkx*zt1)
                - afdt*cu[1+4*(j+l1)];
            zt9 = exyz[2+4*(j+l1)] + cdt*(dkx*zt2) - afdt*cu[2+4*(j+l1)];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exyz[4*(j+l1)] = zt7;
            exyz[1+4*(j+l1)] = zt8;
            exyz[2+4*(j+l1)] = zt9;
            at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
            ws += (double) at1;
            zt4 -= dth*(dkz*zt2);
            zt5 += dth*(dkz*zt3 + dkx*zt1);
            zt6 -= dth*(dkx*zt2);
            bxyz[4*(j+l1)] = zt4;
            bxyz[1+4*(j+l1)] = zt5;
            bxyz[2+4*(j+l1)] = zt6;
            at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
            wp += (double) at1;
            bxyz[4*(j+k1+l1)] = zero;
            bxyz[1+4*(j+k1+l1)] = zero;
            bxyz[2+4*(j+k1+l1)] = zero;
            exyz[4*(j+k1+l1)] = zero;
            exyz[1+4*(j+k1+l1)] = zero;
            exyz[2+4*(j+k1+l1)] = zero;
         }
/* mode numbers kx = 0, nx/2 */
         afdt = adt*cimagf(ffc[ll]);
/* update magnetic field half time step */
         zt2 = -cimagf(exyz[1+4*(lj)]) + crealf(exyz[1+4*(lj)])*_Complex_I;
         zt3 = -cimagf(exyz[4*(lj)]) + crealf(exyz[4*(lj)])*_Complex_I;
         zt4 = bxyz[4*lj] + dth*(dkz*zt2);
         zt5 = bxyz[1+4*lj] - dth*(dkz*zt3);
/* update electric field whole time step */
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exyz[4*lj] - cdt*(dkz*zt2) - afdt*cu[4*lj];
         zt8 = exyz[1+4*lj] + cdt*(dkz*zt3) - afdt*cu[1+4*lj];
/* update magnetic field half time step and store electric field */
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exyz[4*lj] = zt7;
         exyz[1+4*lj] = zt8;
         exyz[2+4*lj] = zero;
         at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8));
         ws += (double) at1;
         zt4 += dth*(dkz*zt2);
         zt5 -= dth*(dkz*zt3);
         bxyz[4*lj] = zt4;
         bxyz[1+4*lj] = zt5;
         bxyz[2+4*lj] = zero;
         at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5));
         wp += (double) at1;
         bxyz[4*(k1+lj)] = zero;
         bxyz[1+4*(k1+lj)] = zero;
         bxyz[2+4*(k1+lj)] = zero;
         exyz[4*(k1+lj)] = zero;
         exyz[1+4*(k1+lj)] = zero;
         exyz[2+4*(k1+lj)] = zero;
         bxyz[4*l1] = zero;
         bxyz[1+4*l1] = zero;
         bxyz[2+4*l1] = zero;
         exyz[4*l1] = zero;
         exyz[1+4*l1] = zero;
         exyz[2+4*l1] = zero;
         bxyz[4*(k1+l1)] = zero;
         bxyz[1+4*(k1+l1)] = zero;
         bxyz[2+4*(k1+l1)] = zero;
         exyz[4*(k1+l1)] = zero;
         exyz[1+4*(k1+l1)] = zero;
/*       sum1 += ws; */
         _mm512_store_pd(&dd[0],v_ws);
         for (j = 1; j < 8; j++) {
            dd[0] += dd[j];
         }
         sum1 += (ws + dd[0]);
/*       sum2 += wp; */
         _mm512_store_pd(&dd[0],v_wp);
         for (j = 1; j < 8; j++) {
            dd[0] += dd[j];
         }
         sum2 += (wp + dd[0]);
      }
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   sum3 = 0.0;
   sum4 = 0.0;
#pragma omp parallel for \
private(j,k,k1,kk,kj,dky,dkx,afdt,at1,zt1,zt2,zt3,zt4,zt5,zt6,zt7,zt8, \
zt9,ws,wp,v_it,v_dkx,v_dky,v_dk1,v_dk2,v_afdt,v_at2,v_at3,v_zt1,v_zt2, \
v_zt3,v_zt4,v_zt5,v_zt6,v_zt7,v_d,v_ws,v_wp,dd) \
reduction(+:sum3,sum4)
   for (k = 1; k < nyh; k++) {
/*    dky = dny*(float) k; */
      v_it = _mm512_set1_epi32(k);
      v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dky = _mm512_mul_ps(v_dny,v_dky);
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      ws = 0.0;
      wp = 0.0;
      v_ws = _mm512_set1_pd(0.0);
      v_wp = _mm512_set1_pd(0.0);
/* add ky to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(12336),v_dky);
      v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(771),v_dky);
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       dkx = dnx*(float) j; */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                  v_dkx);
         v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                  v_dkx);
/*       afdt = adt*cimagf(ffc[j+kk]); */
         v_afdt = _mm512_mask_loadunpacklo_ps(v_zero,
                  _mm512_int2mask(15),(float *)&ffc[j+kk]);
         v_afdt = _mm512_mask_loadunpackhi_ps(v_afdt,
                  _mm512_int2mask(15),(float *)&ffc[j+kk+8]);
         v_afdt = _mm512_permute4f128_ps(v_afdt,0);
         v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                  _mm512_int2mask(13260),(__m512i)v_afdt,78);
         v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                  _mm512_int2mask(21845),(__m512i)v_afdt,177);
         v_afdt = _mm512_mul_ps(v_adt,v_afdt);
/* update magnetic field half time step, ky > 0 */
/*       zt1 = -cimagf(exyz[2+4*(j+kj)])             */
/*            + crealf(exyz[2+4*(j+kj)])*_Complex_I; */
/*       zt2 = -cimagf(exyz[1+4*(j+kj)])             */
/*            + crealf(exyz[1+4*(j+kj)])*_Complex_I; */
/*       zt3 = -cimagf(exyz[4*(j+kj)])               */
/*            + crealf(exyz[4*(j+kj)])*_Complex_I;   */
         v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+kj)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                 v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       zt4 = bxyz[4*(j+kj)] - dth*(dky*zt1);             */
/*       zt5 = bxyz[1+4*(j+kj)] + dth*(dkx*zt1);           */
/*       zt6 = bxyz[2+4*(j+kj)] - dth*(dkx*zt2 - dky*zt3); */
         v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+kj)]);
         v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*       zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*       zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*       zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),v_zero,
                 v_zt5);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       zt7 = exyz[4*(j+kj)] + cdt*(dky*zt1) - afdt*cu[4*(j+kj)];     */
/*       zt8 = exyz[1+4*(j+kj)] - cdt*(dkx*zt1) - afdt*cu[1+4*(j+kj)]; */
/*       zt9 = exyz[2+4*(j+kj)] + cdt*(dkx*zt2 - dky*zt3)              */
/*           - afdt*cu[2+4*(j+kj)];                                    */
         v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                 v_zt4);
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj)]);
         v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                 v_afdt,v_zt2);
         v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*       zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*       zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*       zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                 v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       exyz[4*(j+kj)] = zt7;   */
/*       exyz[1+4*(j+kj)] = zt8; */
/*       exyz[2+4*(j+kj)] = zt9; */
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                    v_zero);
            _mm512_mask_store_ps((float *)&exyz[4*(j+kj)],
            _mm512_int2mask(65280),v_zt4);
         }
         else {
            _mm512_store_ps((float *)&exyz[4*(j+kj)],v_zt4);
         }
/*       ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
         v_zt6 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4));
/*       zt4 -= dth*(dky*zt1);           */
/*       zt5 += dth*(dkx*zt1);           */
/*       zt6 -= dth*(dkx*zt2 - dky*zt3); */
         v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*       bxyz[4*(j+kj)] = zt4;   */
/*       bxyz[1+4*(j+kj)] = zt5; */
/*       bxyz[2+4*(j+kj)] = zt6; */
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                    v_zero);
            _mm512_mask_store_ps((float *)&bxyz[4*(j+kj)],
            _mm512_int2mask(65280),v_zt5);
         }
         else {
            _mm512_store_ps((float *)&bxyz[4*(j+kj)],v_zt5);
         }
/*       wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
         v_zt7 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5));
/* update magnetic field half time step, ky < 0 */
/*       zt1 = -cimagf(exyz[2+4*(j+k1)])             */
/*            + crealf(exyz[2+4*(j+k1)])*_Complex_I; */
/*       zt2 = -cimagf(exyz[1+4*(j+k1)])             */
/*            + crealf(exyz[1+4*(j+k1)])*_Complex_I; */
/*       zt3 = -cimagf(exyz[4*(j+k1)])               */
/*            + crealf(exyz[4*(j+k1)])*_Complex_I;   */
         v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+k1)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                 v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
         v_at2 = _mm512_mask_sub_ps(v_dk1,_mm512_int2mask(12336),v_zero,
                  v_dk1);
         v_at3 = _mm512_mask_sub_ps(v_dk2,_mm512_int2mask(771),v_zero,
                 v_dk2);
/*       zt4 = bxyz[4*(j+k1)] + dth*(dky*zt1);             */
/*       zt5 = bxyz[1+4*(j+k1)] + dth*(dkx*zt1);           */
/*       zt6 = bxyz[2+4*(j+k1)] - dth*(dkx*zt2 + dky*zt3); */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+k1)]);
         v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*       zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*       zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*       zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),v_zero,
                 v_zt5);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       zt7 = exyz[4*(j+k1)] - cdt*(dky*zt1) - afdt*cu[4*(j+k1)];     */
/*       zt8 = exyz[1+4*(j+k1)] - cdt*(dkx*zt1) - afdt*cu[1+4*(j+k1)]; */
/*       zt9 = exyz[2+4*(j+k1)] + cdt*(dkx*zt2 + dky*zt3)              */
/*           - afdt*cu[2+4*(j+k1)];                                    */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),v_zt4);
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1)]);
         v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                 v_afdt,v_zt2);
         v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*       zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*       zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*       zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                 v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       exyz[4*(j+k1)] = zt7;   */
/*       exyz[1+4*(j+k1)] = zt8; */
/*       exyz[2+4*(j+k1)] = zt9; */
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                    v_zero);
            _mm512_mask_store_ps((float *)&exyz[4*(j+k1)],
            _mm512_int2mask(65280),v_zt4);
         }
         else {
            _mm512_store_ps((float *)&exyz[4*(j+k1)],v_zt4);
         }
/*       ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
         v_zt6 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4),
                 v_zt6);
/*       zt4 += dth*(dky*zt1);           */
/*       zt5 += dth*(dkx*zt1);           */
/*       zt6 -= dth*(dkx*zt2 + dky*zt3); */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*       bxyz[4*(j+k1)] = zt4;   */
/*       bxyz[1+4*(j+k1)] = zt5; */
/*       bxyz[2+4*(j+k1)] = zt6; */
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                    v_zero);
            _mm512_mask_store_ps((float *)&bxyz[4*(j+k1)],
            _mm512_int2mask(65280),v_zt5);
         }
         else {
            _mm512_store_ps((float *)&bxyz[4*(j+k1)],v_zt5);
         }
/*       wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
         v_zt7 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5),
                 v_zt7);
/* convert to double precision before accumulating */
         v_ws = _mm512_add_pd(v_ws,_mm512_cvtpslo_pd(v_zt6));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt6,78));
         v_ws = _mm512_add_pd(v_ws,v_d);
         v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt7));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt7,78));
         v_wp = _mm512_add_pd(v_wp,v_d);
/*       bxyz[4*(j+kj+l1)] = zero; */
/*       bxyz[1+4*(j+kj+l1)] = zero; */
/*       bxyz[2+4*(j+kj+l1)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+kj+l1)],v_zero);
/*       exyz[4*(j+kj+l1)] = zero; */
/*       exyz[1+4*(j+kj+l1)] = zero; */
/*       exyz[2+4*(j+kj+l1)] = zero; */
         _mm512_store_ps((float *)&exyz[4*(j+kj+l1)],v_zero);
/*       bxyz[4*(j+k1+l1)] = zero; */
/*       bxyz[1+4*(j+k1+l1)] = zero; */
/*       bxyz[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
/*       exyz[4*(j+k1+l1)] = zero; */
/*       exyz[1+4*(j+k1+l1)] = zero; */
/*       exyz[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&exyz[4*(j+k1+l1)],v_zero);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         dkx = dnx*(float) j;
         afdt = adt*cimagf(ffc[j+kk]);
/* update magnetic field half time step, ky > 0 */
         zt1 = -cimagf(exyz[2+4*(j+kj)])
              + crealf(exyz[2+4*(j+kj)])*_Complex_I;
         zt2 = -cimagf(exyz[1+4*(j+kj)])
              + crealf(exyz[1+4*(j+kj)])*_Complex_I;
         zt3 = -cimagf(exyz[4*(j+kj)])
              + crealf(exyz[4*(j+kj)])*_Complex_I;
         zt4 = bxyz[4*(j+kj)] - dth*(dky*zt1);
         zt5 = bxyz[1+4*(j+kj)] + dth*(dkx*zt1);
         zt6 = bxyz[2+4*(j+kj)] - dth*(dkx*zt2 - dky*zt3);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exyz[4*(j+kj)] + cdt*(dky*zt1) - afdt*cu[4*(j+kj)];
         zt8 = exyz[1+4*(j+kj)] - cdt*(dkx*zt1) - afdt*cu[1+4*(j+kj)];
         zt9 = exyz[2+4*(j+kj)] + cdt*(dkx*zt2 - dky*zt3) 
             - afdt*cu[2+4*(j+kj)];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exyz[4*(j+kj)] = zt7;
         exyz[1+4*(j+kj)] = zt8;
         exyz[2+4*(j+kj)] = zt9;
         at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         ws += (double) at1;
         zt4 -= dth*(dky*zt1);
         zt5 += dth*(dkx*zt1);
         zt6 -= dth*(dkx*zt2 - dky*zt3);
         bxyz[4*(j+kj)] = zt4;
         bxyz[1+4*(j+kj)] = zt5;
         bxyz[2+4*(j+kj)] = zt6;
         at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
         wp += (double) at1;
/* update magnetic field half time step, ky < 0 */
         zt1 = -cimagf(exyz[2+4*(j+k1)])
              + crealf(exyz[2+4*(j+k1)])*_Complex_I;
         zt2 = -cimagf(exyz[1+4*(j+k1)])
              + crealf(exyz[1+4*(j+k1)])*_Complex_I;
         zt3 = -cimagf(exyz[4*(j+k1)])
              + crealf(exyz[4*(j+k1)])*_Complex_I;
         zt4 = bxyz[4*(j+k1)] + dth*(dky*zt1);
         zt5 = bxyz[1+4*(j+k1)] + dth*(dkx*zt1);
         zt6 = bxyz[2+4*(j+k1)] - dth*(dkx*zt2 + dky*zt3);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exyz[4*(j+k1)] - cdt*(dky*zt1) - afdt*cu[4*(j+k1)];
         zt8 = exyz[1+4*(j+k1)] - cdt*(dkx*zt1) - afdt*cu[1+4*(j+k1)];
         zt9 = exyz[2+4*(j+k1)] + cdt*(dkx*zt2 + dky*zt3)
             - afdt*cu[2+4*(j+k1)];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exyz[4*(j+k1)] = zt7;
         exyz[1+4*(j+k1)] = zt8;
         exyz[2+4*(j+k1)] = zt9;
         at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         ws += (double) at1;
         zt4 += dth*(dky*zt1);
         zt5 += dth*(dkx*zt1);
         zt6 -= dth*(dkx*zt2 + dky*zt3);
         bxyz[4*(j+k1)] = zt4;
         bxyz[1+4*(j+k1)] = zt5;
         bxyz[2+4*(j+k1)] = zt6;
         at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
         wp += (double) at1;
         bxyz[4*(j+kj+l1)] = zero;
         bxyz[1+4*(j+kj+l1)] = zero;
         bxyz[2+4*(j+kj+l1)] = zero;
         exyz[4*(j+kj+l1)] = zero;
         exyz[1+4*(j+kj+l1)] = zero;
         exyz[2+4*(j+kj+l1)] = zero;
         bxyz[4*(j+k1+l1)] = zero;
         bxyz[1+4*(j+k1+l1)] = zero;
         bxyz[2+4*(j+k1+l1)] = zero;
         exyz[4*(j+k1+l1)] = zero;
         exyz[1+4*(j+k1+l1)] = zero;
         exyz[2+4*(j+k1+l1)] = zero;
      }
/*    sum3 += ws; */
      _mm512_store_pd(&dd[0],v_ws);
      for (j = 1; j < 8; j++) {
         dd[0] += dd[j];
      }
      sum3 += (ws + dd[0]);
/*    sum4 += wp; */
      _mm512_store_pd(&dd[0],v_wp);
      for (j = 1; j < 8; j++) {
          dd[0] += dd[j];
      }
      sum4 += (wp + dd[0]);
   }
/* mode numbers kx = 0, nx/2 */
   ws = 0.0;
   wp = 0.0;
   v_ws = _mm512_set1_pd(0.0);
   v_wp = _mm512_setzero_pd();
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      afdt = adt*cimagf(ffc[kk]);
/* update magnetic field half time step */
      zt1 = -cimagf(exyz[2+4*(kj)]) + crealf(exyz[2+4*(kj)])*_Complex_I;
      zt3 = -cimagf(exyz[4*(kj)]) + crealf(exyz[4*(kj)])*_Complex_I;
      zt4 = bxyz[4*kj] - dth*(dky*zt1);
      zt6 = bxyz[2+4*kj] + dth*(dky*zt3);
/* update electric field whole time step */
      zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
      zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
      zt7 = exyz[4*kj] + cdt*(dky*zt1) - afdt*cu[4*kj];
      zt9 = exyz[2+4*kj] - cdt*(dky*zt3) - afdt*cu[2+4*kj];
/* update magnetic field half time step and store electric field */
      zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
      zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
      exyz[4*kj] = zt7;
      exyz[1+4*kj] = zero;
      exyz[2+4*kj] = zt9;
      at1 = anorm*(zt7*conjf(zt7) + zt9*conjf(zt9));
      ws += (double) at1;
      zt4 -= dth*(dky*zt1);
      zt6 += dth*(dky*zt3);
      bxyz[4*kj] = zt4;
      bxyz[1+4*kj] = zero;
      bxyz[2+4*kj] = zt6;
      at1 = anorm*(zt4*conjf(zt4) + zt6*conjf(zt6));
      wp += (double) at1;
      bxyz[4*k1] = zero;
      bxyz[1+4*k1] = zero;
      bxyz[2+4*k1] = zero;
      exyz[4*k1] = zero;
      exyz[1+4*k1] = zero;
      exyz[2+4*k1] = zero;
      bxyz[4*(kj+l1)] = zero;
      bxyz[1+4*(kj+l1)] = zero;
      bxyz[2+4*(kj+l1)]= zero;
      exyz[4*(kj+l1)] = zero;
      exyz[1+4*(kj+l1)] = zero;
      exyz[2+4*(kj+l1)] = zero;
      bxyz[4*(k1+l1)] = zero;
      bxyz[1+4*(k1+l1)] = zero;
      bxyz[2+4*(k1+l1)] = zero;
      exyz[4*(k1+l1)] = zero;
      exyz[1+4*(k1+l1)] = zero;
      exyz[2+4*(k1+l1)] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
   for (j = 0; j < nxhs; j+=2) {
/*    dkx = dnx*(float) j; */
      v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
      v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(3084),v_dkx);
      v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(12336),v_dkx);
/*    afdt = adt*cimagf(ffc[j]); */
      v_afdt = _mm512_mask_loadunpacklo_ps(v_zero,
                _mm512_int2mask(15),(float *)&ffc[j]);
      v_afdt = _mm512_mask_loadunpackhi_ps(v_afdt,
               _mm512_int2mask(15),(float *)&ffc[j+8]);
      v_afdt = _mm512_permute4f128_ps(v_afdt,0);
      v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
               _mm512_int2mask(13260),(__m512i)v_afdt,78);
      v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
               _mm512_int2mask(21845),(__m512i)v_afdt,177);
      v_afdt = _mm512_mul_ps(v_adt,v_afdt);
/* update magnetic field half time step */
/*    zt1 = -cimagf(exyz[2+4*j]) + crealf(exyz[2+4*j])*_Complex_I; */
/*    zt2 = -cimagf(exyz[1+4*j]) + crealf(exyz[1+4*j])*_Complex_I; */
      v_zt4 = _mm512_load_ps((float *)&exyz[4*j]);
      v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
               v_zt4);
      v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*    zt5 = bxyz[1+4*j] + dth*(dkx*zt1); */
/*    zt6 = bxyz[2+4*j] - dth*(dkx*zt2); */
      v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
      v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
      v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
      v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
      v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
      v_zt2 = _mm512_load_ps((float *)&bxyz[4*j]);
      v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*    zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*    zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
      v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),v_zero,
              v_zt5);
      v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*    zt8 = exyz[1+4*j] - cdt*(dkx*zt1) - afdt*cu[1+4*j]; */
/*    zt9 = exyz[2+4*j] + cdt*(dkx*zt2) - afdt*cu[2+4*j]; */
      v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
      v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
      v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
      v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
      v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),v_zt4);
      v_zt2 = _mm512_load_ps((float *)&cu[4*j]);
      v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),v_afdt,
              v_zt2);
      v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*    zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*    zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
      v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
              v_zt4);
      v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*    exyz[4*j] = zero;  */
/*    exyz[1+4*j] = zt8; */
/*    exyz[2+4*j] = zt9; */
/* zero out kx = 0 mode */
      if (j==0) {
         v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),v_zero);
         _mm512_mask_store_ps((float *)&exyz[4*j],
         _mm512_int2mask(65280),v_zt4);
      }
      else {
         _mm512_store_ps((float *)&exyz[4*j],v_zt4);
      }
/*    ws += anorm*(zt8*conjf(zt8) + zt9*conjf(zt9)); */
      v_zt6 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4));
/*    zt5 += dth*(dkx*zt1); */
/*    zt6 -= dth*(dkx*zt2); */
      v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
      v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
      v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
      v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
      v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
      v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*    bxyz[4*j] = zero;  */
/*    bxyz[1+4*j] = zt5; */
/*    bxyz[2+4*j] = zt6; */
/* zero out kx = 0 mode */
      if (j==0) {
         v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),v_zero);
         _mm512_mask_store_ps((float *)&bxyz[4*j],
         _mm512_int2mask(65280),v_zt5);
      }
      else {
         _mm512_store_ps((float *)&bxyz[4*j],v_zt5);
      }
/*    wp += anorm*(zt5*conjf(zt5) + zt6*conjf(zt6)); */
      v_zt7 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5));
/* convert to double precision before accumulating */
      v_ws = _mm512_add_pd(v_ws,_mm512_cvtpslo_pd(v_zt6));
      v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt6,78));
      v_ws = _mm512_add_pd(v_ws,v_d);
      v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt7));
      v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt7,78));
      v_wp = _mm512_add_pd(v_wp,v_d);
/*    bxyz[4*(j+k1)] = zero;   */
/*    bxyz[1+4*(j+k1)] = zero; */
/*    bxyz[2+4*(j+k1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+k1)],v_zero);
/*    exyz[4*(j+k1)] = zero;   */
/*    exyz[1+4*(j+k1)] = zero; */
/*    exyz[2+4*(j+k1)] = zero; */
      _mm512_store_ps((float *)&exyz[4*(j+k1)],v_zero);
/*    bxyz[4*(j+l1)] = zero;   */
/*    bxyz[1+4*(j+l1)] = zero; */
/*    bxyz[2+4*(j+l1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+l1)],v_zero);
/*    exyz[4*(j+l1)] = zero;   */
/*    exyz[1+4*(j+l1)] = zero; */
/*    exyz[2+4*(j+l1)] = zero; */
      _mm512_store_ps((float *)&exyz[4*(j+l1)],v_zero);
/*    bxyz[4*(j+k1+l1)] = zero;   */
/*    bxyz[1+4*(j+k1+l1)] = zero; */
/*    bxyz[2+4*(j+k1+l1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
/*    exyz[4*(j+k1+l1)] = zero;   */
/*    exyz[1+4*(j+k1+l1)] = zero; */
/*    exyz[2+4*(j+k1+l1)] = zero; */
      _mm512_store_ps((float *)&exyz[4*(j+k1+l1)],v_zero);
   }
/* loop over remaining elements */
   for (j = itn; j < nxh; j++) {
      dkx = dnx*(float) j;
      afdt = adt*cimagf(ffc[j]);
/* update magnetic field half time step */
      zt1 = -cimagf(exyz[2+4*j]) + crealf(exyz[2+4*j])*_Complex_I;
      zt2 = -cimagf(exyz[1+4*j]) + crealf(exyz[1+4*j])*_Complex_I;
      zt5 = bxyz[1+4*j] + dth*(dkx*zt1);
      zt6 = bxyz[2+4*j] - dth*(dkx*zt2);
/* update electric field whole time step */
      zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
      zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
      zt8 = exyz[1+4*j] - cdt*(dkx*zt1) - afdt*cu[1+4*j];
      zt9 = exyz[2+4*j] + cdt*(dkx*zt2) - afdt*cu[2+4*j];
/* update magnetic field half time step and store electric field */
      zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
      zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
      exyz[4*j] = zero;
      exyz[1+4*j] = zt8;
      exyz[2+4*j] = zt9;
      at1 = anorm*(zt8*conjf(zt8) + zt9*conjf(zt9));
      ws += (double) at1;
      zt5 += dth*(dkx*zt1);
      zt6 -= dth*(dkx*zt2);
      bxyz[4*j] = zero;
      bxyz[1+4*j] = zt5;
      bxyz[2+4*j] = zt6;
      at1 = anorm*(zt5*conjf(zt5) + zt6*conjf(zt6));
      wp += (double) at1;
      bxyz[4*(j+k1)] = zero;
      bxyz[1+4*(j+k1)] = zero;
      bxyz[2+4*(j+k1)] = zero;
      exyz[4*(j+k1)] = zero;
      exyz[1+4*(j+k1)] = zero;
      exyz[2+4*(j+k1)] = zero;
      bxyz[4*(j+l1)] = zero;
      bxyz[1+4*(j+l1)] = zero;
      bxyz[2+4*(j+l1)] = zero;
      exyz[4*(j+l1)] = zero;
      exyz[1+4*(j+l1)] = zero;
      exyz[2+4*(j+l1)] = zero;
      bxyz[4*(j+k1+l1)] = zero;
      bxyz[1+4*(j+k1+l1)] = zero;
      bxyz[2+4*(j+k1+l1)] = zero;
      exyz[4*(j+k1+l1)] = zero;
      exyz[1+4*(j+k1+l1)] = zero;
      exyz[2+4*(j+k1+l1)] = zero;
   }
   bxyz[0] = zero;
   bxyz[1] = zero;
   bxyz[2] = zero;
   exyz[0] = zero;
   exyz[1] = zero;
   exyz[2]= zero;
   bxyz[4*k1] = zero;
   bxyz[1+4*k1] = zero;
   bxyz[2+4*k1] = zero;
   exyz[4*k1] = zero;
   exyz[1+4*k1] = zero;
   exyz[2+4*k1] = zero;
   bxyz[4*l1] = zero;
   bxyz[1+4*l1] = zero;
   bxyz[2+4*l1] = zero;
   exyz[4*l1] = zero;
   exyz[1+4*l1] = zero;
   exyz[2+4*l1] = zero;
   bxyz[4*(k1+l1)] = zero;
   bxyz[1+4*(k1+l1)] = zero;
   bxyz[2+4*(k1+l1)] = zero;
   exyz[4*(k1+l1)] = zero;
   exyz[1+4*(k1+l1)] = zero;
   exyz[2+4*(k1+l1)] = zero;
/* *wf = ws*((float) nx)*((float) ny)*((float) nz);    */
/* sum3 += ws; */
   _mm512_store_pd(&dd[0],v_ws);
   for (j = 1; j < 8; j++) {
      dd[0] += dd[j];
   }
   sum3 += (ws + dd[0]);
   *wf = (sum1 + sum3)*((float) nx)*((float) ny)*((float) nz);
/* *wm = c2*wp*((float) nx)*((float) ny)*((float) nz); */

/* sum4 += wp; */
   _mm512_store_pd(&dd[0],v_wp);
   for (j = 1; j < 8; j++) {
      dd[0] += dd[j];
   }
   sum4 += (wp + dd[0]);
   *wm = c2*(sum2 + sum4)*((float) nx)*((float) ny)*((float) nz);
   return;
}

/*--------------------------------------------------------------------*/
void ckncmemfield3(float complex fxyz[], float complex exyz[],
                   float complex ffc[], int isign, int nx, int ny,
                   int nz, int nxvh, int nyv, int nzv, int nxhd,
                   int nyhd, int nzhd) {
/* this subroutine either adds complex vector fields if isign > 0
   or copies complex vector fields if isign < 0
   includes additional smoothing
   requires KNC, fxyz, exyz, ffc need to be 64 byte aligned
   nxhd needs to be a multiple of 8
   nxvh needs to be a multiple of 2
   fxyz, exyz needs to have 4 components
local data                                                 */
   int j, k, l, nxh, nyh, nzh, nxhs, itn, k1, l1, kk, kj, ll, lj;
   int nxyhd, nxvyh; 
   float at1;
   __m512 v_at1, v_zero, v_zt1, v_zt2;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxhs = 2*(nxh/2);
   itn = 1 > nxhs ? 1 : nxhs;
   nxyhd = nxhd*nyhd;
   nxvyh = nxvh*nyv;
   v_zero = _mm512_setzero_ps();
/* add the fields */
   if (isign > 0) {
#pragma omp parallel
      {
#pragma omp for nowait \
private(j,k,l,k1,l1,kk,kj,ll,lj,at1,v_at1,v_zt1,v_zt2)
         for (l = 1; l < nzh; l++) {
            ll = nxyhd*l;
            lj = nxvyh*l;
            l1 = nxvyh*nz - lj;
            for (k = 1; k < nyh; k++) {
               kk = nxhd*k;
               kj = nxvh*k;
               k1 = nxvh*ny - kj;
/* vector loop over elements in blocks of 2 */
               for (j = 0; j < nxhs; j+=2) {
/*                at1 = cimagf(ffc[j+kk+ll]); */
                  v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                          _mm512_int2mask(15),(float *)&ffc[j+kk+ll]);
                  v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                          _mm512_int2mask(15),(float *)&ffc[j+kk+ll+8]);
                  v_at1 = _mm512_permute4f128_ps(v_at1,0);
                  v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                         _mm512_int2mask(13260),(__m512i)v_at1,78);
                  v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                          _mm512_int2mask(21845),(__m512i)v_at1,177);
/*                fxyz[4*(j+kj+lj)] += exyz[4*(j+kj+lj)]*at1;     */
/*                fxyz[1+4*(j+kj+lj)] += exyz[1+4*(j+kj+lj)]*at1; */
/*                fxyz[2+4*(j+kj+lj)] += exyz[2+4*(j+kj+lj)]*at1; */
                  v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+lj)]);
                  v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+kj+lj)]);
                  v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
                  _mm512_store_ps((float *)&fxyz[4*(j+kj+lj)],v_zt2);
/*                fxyz[4*(j+k1+lj)] += exyz[4*(j+k1+lj)]*at1;     */
/*                fxyz[1+4*(j+k1+lj)] += exyz[1+4*(j+k1+lj)]*at1; */
/*                fxyz[2+4*(j+k1+lj)] += exyz[2+4*(j+k1+lj)]*at1; */
                  v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+lj)]);
                  v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+lj)]);
                  v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
                  _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],v_zt2);
/*                fxyz[4*(j+kj+l1)] += exyz[4*(j+kj+l1)]*at1;     */
/*                fxyz[1+4*(j+kj+l1)] += exyz[1+4*(j+kj+l1)]*at1; */
/*                fxyz[2+4*(j+kj+l1)] += exyz[2+4*(j+kj+l1)]*at1; */
                  v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+l1)]);
                  v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+kj+l1)]);
                  v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
                  _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],v_zt2);
/*                fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;     */
/*                fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1; */
/*                fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1; */
                  v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
                  v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+l1)]);
                  v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
                  _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
               }
/* loop over remaining elements */
               for (j = itn; j < nxh; j++) {
                  at1 = cimagf(ffc[j+kk+ll]);
                  fxyz[4*(j+kj+lj)] += exyz[4*(j+kj+lj)]*at1;
                  fxyz[1+4*(j+kj+lj)] += exyz[1+4*(j+kj+lj)]*at1;
                  fxyz[2+4*(j+kj+lj)] += exyz[2+4*(j+kj+lj)]*at1;
                  fxyz[4*(j+k1+lj)] += exyz[4*(j+k1+lj)]*at1;
                  fxyz[1+4*(j+k1+lj)] += exyz[1+4*(j+k1+lj)]*at1;
                  fxyz[2+4*(j+k1+lj)] += exyz[2+4*(j+k1+lj)]*at1;
                  fxyz[4*(j+kj+l1)] += exyz[4*(j+kj+l1)]*at1;
                  fxyz[1+4*(j+kj+l1)] += exyz[1+4*(j+kj+l1)]*at1;
                  fxyz[2+4*(j+kj+l1)] += exyz[2+4*(j+kj+l1)]*at1;
                  fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;
                  fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1;
                  fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1;
               }
            }
            k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nxhs; j+=2) {
/*             at1 = cimagf(ffc[j+ll]); */
               v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                       _mm512_int2mask(15),(float *)&ffc[j+ll]);
               v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                       _mm512_int2mask(15),(float *)&ffc[j+kk+8]);
               v_at1 = _mm512_permute4f128_ps(v_at1,0);
               v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                       _mm512_int2mask(13260),(__m512i)v_at1,78);
               v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                       _mm512_int2mask(21845),(__m512i)v_at1,177);
/*             fxyz[4*(j+lj)] += exyz[4*(j+lj)]*at1;     */
/*             fxyz[1+4*(j+lj)] += exyz[1+4*(j+lj)]*at1; */
/*             fxyz[2+4*(j+lj)] += exyz[2+4*(j+lj)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+lj)]);
               v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+lj)]);
               v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
               _mm512_store_ps((float *)&fxyz[4*(j+lj)],v_zt2);
/*             fxyz[4*(j+k1+lj)] += exyz[4*(j+k1+lj)]*at1;     */
/*             fxyz[1+4*(j+k1+lj)] += exyz[1+4*(j+k1+lj)]*at1; */
/*             fxyz[2+4*(j+k1+lj)] += exyz[2+4*(j+k1+lj)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+lj)]);
               v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+lj)]);
               v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
               _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],v_zt2);
/*             fxyz[4*(j+l1)] += exyz[4*(j+l1)]*at1;     */
/*             fxyz[1+4*(j+l1)] += exyz[1+4*(j+l1)]*at1; */
/*             fxyz[2+4*(j+l1)] += exyz[2+4*(j+l1)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+l1)]);
               v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+l1)]);
               v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
               _mm512_store_ps((float *)&fxyz[4*(j+l1)],v_zt2);
/*             fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;     */
/*             fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1; */
/*             fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
               v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+l1)]);
               v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
               _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
            }
/* loop over remaining elements */
            for (j = itn; j < nxh; j++) {
               at1 = cimagf(ffc[j+ll]);
               fxyz[4*(j+lj)] += exyz[4*(j+lj)]*at1;
               fxyz[1+4*(j+lj)] += exyz[1+4*(j+lj)]*at1;
               fxyz[2+4*(j+lj)] += exyz[2+4*(j+lj)]*at1;
               fxyz[4*(j+k1+lj)] += exyz[4*(j+k1+lj)]*at1;
               fxyz[1+4*(j+k1+lj)] += exyz[1+4*(j+k1+lj)]*at1;
               fxyz[2+4*(j+k1+lj)] += exyz[2+4*(j+k1+lj)]*at1;
               fxyz[4*(j+l1)] += exyz[4*(j+l1)]*at1;
               fxyz[1+4*(j+l1)] += exyz[1+4*(j+l1)]*at1;
               fxyz[2+4*(j+l1)] += exyz[2+4*(j+l1)]*at1;
               fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;
               fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1;
               fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1;
            }
         }
      }
      l1 = nxvyh*nzh;
#pragma omp parallel for private(j,k,k1,kk,kj,at1,v_at1,v_zt1,v_zt2)
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          at1 = cimagf(ffc[j+kk]); */
            v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                    _mm512_int2mask(15),(float *)&ffc[j+kk]);
            v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                    _mm512_int2mask(15),(float *)&ffc[j+kk+8]);
            v_at1 = _mm512_permute4f128_ps(v_at1,0);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                    _mm512_int2mask(13260),(__m512i)v_at1,78);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                     _mm512_int2mask(21845),(__m512i)v_at1,177);
/*          fxyz[4*(j+kj)] += exyz[4*(j+kj)]*at1;     */
/*          fxyz[1+4*(j+kj)] += exyz[1+4*(j+kj)]*at1; */
/*          fxyz[2+4*(j+kj)] += exyz[2+4*(j+kj)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+kj)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+kj)],v_zt2);
/*          fxyz[4*(j+k1)] += exyz[4*(j+k1)]*at1;     */
/*          fxyz[1+4*(j+k1)] += exyz[1+4*(j+k1)]*at1; */
/*          fxyz[2+4*(j+k1)] += exyz[2+4*(j+k1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+k1)],v_zt2);
/*          fxyz[4*(j+kj+l1)] += exyz[4*(j+kj+l1)]*at1;     */
/*          fxyz[1+4*(j+kj+l1)] += exyz[1+4*(j+kj+l1)]*at1; */
/*          fxyz[2+4*(j+kj+l1)] += exyz[2+4*(j+kj+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+l1)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+kj+l1)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],v_zt2);
/*          fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;     */
/*          fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1; */
/*          fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+l1)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            at1 = cimagf(ffc[j+kk]);
            fxyz[4*(j+kj)] += exyz[4*(j+kj)]*at1;
            fxyz[1+4*(j+kj)] += exyz[1+4*(j+kj)]*at1;
            fxyz[2+4*(j+kj)] += exyz[2+4*(j+kj)]*at1;
            fxyz[4*(j+k1)] += exyz[4*(j+k1)]*at1;
            fxyz[1+4*(j+k1)] += exyz[1+4*(j+k1)]*at1;
            fxyz[2+4*(j+k1)] += exyz[2+4*(j+k1)]*at1;
            fxyz[4*(j+kj+l1)] += exyz[4*(j+kj+l1)]*at1;
            fxyz[1+4*(j+kj+l1)] += exyz[1+4*(j+kj+l1)]*at1;
            fxyz[2+4*(j+kj+l1)] += exyz[2+4*(j+kj+l1)]*at1;
            fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;
            fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1;
            fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1;
         }
      }
      k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       at1 = cimagf(ffc[j]); */
         v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,_mm512_int2mask(15),
                 (float *)&ffc[j]);
         v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                 _mm512_int2mask(15),(float *)&ffc[j+8]);
         v_at1 = _mm512_permute4f128_ps(v_at1,0);
         v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                 _mm512_int2mask(13260),(__m512i)v_at1,78);
         v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                  _mm512_int2mask(21845),(__m512i)v_at1,177);
/*       fxyz[4*j] += exyz[4*j]*at1;     */
/*       fxyz[1+4*j] += exyz[1+4*j]*at1; */
/*       fxyz[2+4*j] += exyz[2+4*j]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*j]);
         v_zt2 = _mm512_load_ps((float *)&fxyz[4*j]);
         v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
         _mm512_store_ps((float *)&fxyz[4*j],v_zt2);
/*       fxyz[4*(j+k1)] += exyz[4*(j+k1)]*at1;     */
/*       fxyz[1+4*(j+k1)] += exyz[1+4*(j+k1)]*at1; */
/*       fxyz[2+4*(j+k1)] += exyz[2+4*(j+k1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1)]);
         v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1)]);
         v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
         _mm512_store_ps((float *)&fxyz[4*(j+k1)],v_zt2);
/*       fxyz[4*(j+l1)] += exyz[4*(j+l1)]*at1;     */
/*       fxyz[1+4*(j+l1)] += exyz[1+4*(j+l1)]*at1; */
/*       fxyz[2+4*(j+l1)] += exyz[2+4*(j+l1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+l1)]);
         v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+l1)]);
         v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
         _mm512_store_ps((float *)&fxyz[4*(j+l1)],v_zt2);
/*       fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;     */
/*       fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1; */
/*       fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
         v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+l1)]);
         v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
         _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         at1 = cimagf(ffc[j]);
         fxyz[4*j] += exyz[4*j]*at1;
         fxyz[1+4*j] += exyz[1+4*j]*at1;
         fxyz[2+4*j] += exyz[2+4*j]*at1;
         fxyz[4*(j+k1)] += exyz[4*(j+k1)]*at1;
         fxyz[1+4*(j+k1)] += exyz[1+4*(j+k1)]*at1;
         fxyz[2+4*(j+k1)] += exyz[2+4*(j+k1)]*at1;
         fxyz[4*(j+l1)] += exyz[4*(j+l1)]*at1;
         fxyz[1+4*(j+l1)] += exyz[1+4*(j+l1)]*at1;
         fxyz[2+4*(j+l1)] += exyz[2+4*(j+l1)]*at1;
         fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;
         fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1;
         fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1;
      }
   }
/* copy the fields */
   else if (isign < 0) {
#pragma omp parallel
      {
#pragma omp for nowait \
private(j,k,l,k1,l1,kk,kj,ll,lj,at1,v_at1,v_zt1,v_zt2)
         for (l = 1; l < nzh; l++) {
            ll = nxyhd*l;
            lj = nxvyh*l;
            l1 = nxvyh*nz - lj;
            for (k = 1; k < nyh; k++) {
               kk = nxhd*k;
               kj = nxvh*k;
               k1 = nxvh*ny - kj;
/* vector loop over elements in blocks of 2 */
               for (j = 0; j < nxhs; j+=2) {
/*                at1 = cimagf(ffc[j+kk+ll]); */
                  v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                          _mm512_int2mask(15),(float *)&ffc[j+kk+ll]);
                  v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                          _mm512_int2mask(15),(float *)&ffc[j+kk+ll+8]);
                  v_at1 = _mm512_permute4f128_ps(v_at1,0);
                  v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                         _mm512_int2mask(13260),(__m512i)v_at1,78);
                  v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                          _mm512_int2mask(21845),(__m512i)v_at1,177);
/*                fxyz[4*(j+kj+lj)] = exyz[4*(j+kj+lj)]*at1;     */
/*                fxyz[1+4*(j+kj+lj)] = exyz[1+4*(j+kj+lj)]*at1; */
/*                fxyz[2+4*(j+kj+lj)] = exyz[2+4*(j+kj+lj)]*at1; */
                  v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+lj)]);
                  v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
                  _mm512_store_ps((float *)&fxyz[4*(j+kj+lj)],v_zt2);
/*                fxyz[4*(j+k1+lj)] = exyz[4*(j+k1+lj)]*at1;     */
/*                fxyz[1+4*(j+k1+lj)] = exyz[1+4*(j+k1+lj)]*at1; */
/*                fxyz[2+4*(j+k1+lj)] = exyz[2+4*(j+k1+lj)]*at1; */
                  v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+lj)]);
                  v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
                  _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],v_zt2);
/*                fxyz[4*(j+kj+l1)] = exyz[4*(j+kj+l1)]*at1;     */
/*                fxyz[1+4*(j+kj+l1)] = exyz[1+4*(j+kj+l1)]*at1; */
/*                fxyz[2+4*(j+kj+l1)] = exyz[2+4*(j+kj+l1)]*at1; */
                  v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+l1)]);
                  v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
                  _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],v_zt2);
/*                fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;     */
/*                fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1; */
/*                fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1; */
                  v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
                  v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
                  _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
               }
/* loop over remaining elements */
               for (j = itn; j < nxh; j++) {
                  at1 = cimagf(ffc[j+kk+ll]);
                  fxyz[4*(j+kj+lj)] = exyz[4*(j+kj+lj)]*at1;
                  fxyz[1+4*(j+kj+lj)] = exyz[1+4*(j+kj+lj)]*at1;
                  fxyz[2+4*(j+kj+lj)] = exyz[2+4*(j+kj+lj)]*at1;
                  fxyz[4*(j+k1+lj)] = exyz[4*(j+k1+lj)]*at1;
                  fxyz[1+4*(j+k1+lj)] = exyz[1+4*(j+k1+lj)]*at1;
                  fxyz[2+4*(j+k1+lj)] = exyz[2+4*(j+k1+lj)]*at1;
                  fxyz[4*(j+kj+l1)] = exyz[4*(j+kj+l1)]*at1;
                  fxyz[1+4*(j+kj+l1)] = exyz[1+4*(j+kj+l1)]*at1;
                  fxyz[2+4*(j+kj+l1)] = exyz[2+4*(j+kj+l1)]*at1;
                  fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;
                  fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1;
                  fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1;
               }
            }
            k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nxhs; j+=2) {
/*             at1 = cimagf(ffc[j+ll]); */
               v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                       _mm512_int2mask(15),(float *)&ffc[j+ll]);
               v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                       _mm512_int2mask(15),(float *)&ffc[j+kk+8]);
               v_at1 = _mm512_permute4f128_ps(v_at1,0);
               v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                       _mm512_int2mask(13260),(__m512i)v_at1,78);
               v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                       _mm512_int2mask(21845),(__m512i)v_at1,177);
/*             fxyz[4*(j+lj)] = exyz[4*(j+lj)]*at1;     */
/*             fxyz[1+4*(j+lj)] = exyz[1+4*(j+lj)]*at1; */
/*             fxyz[2+4*(j+lj)] = exyz[2+4*(j+lj)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+lj)]);
               v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
               _mm512_store_ps((float *)&fxyz[4*(j+lj)],v_zt2);
/*             fxyz[4*(j+k1+lj)] = exyz[4*(j+k1+lj)]*at1;     */
/*             fxyz[1+4*(j+k1+lj)] = exyz[1+4*(j+k1+lj)]*at1; */
/*             fxyz[2+4*(j+k1+lj)] = exyz[2+4*(j+k1+lj)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+lj)]);
               v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
               _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],v_zt2);
/*             fxyz[4*(j+l1)] = exyz[4*(j+l1)]*at1;     */
/*             fxyz[1+4*(j+l1)] = exyz[1+4*(j+l1)]*at1; */
/*             fxyz[2+4*(j+l1)] = exyz[2+4*(j+l1)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+l1)]);
               v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
               _mm512_store_ps((float *)&fxyz[4*(j+l1)],v_zt2);
/*             fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;     */
/*             fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1; */
/*             fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
               v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
               _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
            }
/* loop over remaining elements */
            for (j = itn; j < nxh; j++) {
               at1 = cimagf(ffc[j+ll]);
               fxyz[4*(j+lj)] = exyz[4*(j+lj)]*at1;
               fxyz[1+4*(j+lj)] = exyz[1+4*(j+lj)]*at1;
               fxyz[2+4*(j+lj)] = exyz[2+4*(j+lj)]*at1;
               fxyz[4*(j+k1+lj)] = exyz[4*(j+k1+lj)]*at1;
               fxyz[1+4*(j+k1+lj)] = exyz[1+4*(j+k1+lj)]*at1;
               fxyz[2+4*(j+k1+lj)] = exyz[2+4*(j+k1+lj)]*at1;
               fxyz[4*(j+l1)] = exyz[4*(j+l1)]*at1;
               fxyz[1+4*(j+l1)] = exyz[1+4*(j+l1)]*at1;
               fxyz[2+4*(j+l1)] = exyz[2+4*(j+l1)]*at1;
               fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;
               fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1;
               fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1;
            }
         }
      }
      l1 = nxvyh*nzh;
#pragma omp parallel for private(j,k,k1,kk,kj,at1,v_at1,v_zt1,v_zt2)
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          at1 = cimagf(ffc[j+kk]); */
            v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                    _mm512_int2mask(15),(float *)&ffc[j+kk]);
            v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                    _mm512_int2mask(15),(float *)&ffc[j+kk+8]);
            v_at1 = _mm512_permute4f128_ps(v_at1,0);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                    _mm512_int2mask(13260),(__m512i)v_at1,78);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                     _mm512_int2mask(21845),(__m512i)v_at1,177);
/*          fxyz[4*(j+kj)] = exyz[4*(j+kj)]*at1;     */
/*          fxyz[1+4*(j+kj)] = exyz[1+4*(j+kj)]*at1; */
/*          fxyz[2+4*(j+kj)] = exyz[2+4*(j+kj)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+kj)],v_zt2);
/*          fxyz[4*(j+k1)] = exyz[4*(j+k1)]*at1;     */
/*          fxyz[1+4*(j+k1)] = exyz[1+4*(j+k1)]*at1; */
/*          fxyz[2+4*(j+k1)] = exyz[2+4*(j+k1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+k1)],v_zt2);
/*          fxyz[4*(j+kj+l1)] = exyz[4*(j+kj+l1)]*at1;     */
/*          fxyz[1+4*(j+kj+l1)] = exyz[1+4*(j+kj+l1)]*at1; */
/*          fxyz[2+4*(j+kj+l1)] = exyz[2+4*(j+kj+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+l1)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],v_zt2);
/*          fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;     */
/*          fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1; */
/*          fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            at1 = cimagf(ffc[j+kk]);
            fxyz[4*(j+kj)] = exyz[4*(j+kj)]*at1;
            fxyz[1+4*(j+kj)] = exyz[1+4*(j+kj)]*at1;
            fxyz[2+4*(j+kj)] = exyz[2+4*(j+kj)]*at1;
            fxyz[4*(j+k1)] = exyz[4*(j+k1)]*at1;
            fxyz[1+4*(j+k1)] = exyz[1+4*(j+k1)]*at1;
            fxyz[2+4*(j+k1)] = exyz[2+4*(j+k1)]*at1;
            fxyz[4*(j+kj+l1)] = exyz[4*(j+kj+l1)]*at1;
            fxyz[1+4*(j+kj+l1)] = exyz[1+4*(j+kj+l1)]*at1;
            fxyz[2+4*(j+kj+l1)] = exyz[2+4*(j+kj+l1)]*at1;
            fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;
            fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1;
            fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1;
         }
      }
      k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       at1 = cimagf(ffc[j]); */
         v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,_mm512_int2mask(15),
                 (float *)&ffc[j]);
         v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                 _mm512_int2mask(15),(float *)&ffc[j+8]);
         v_at1 = _mm512_permute4f128_ps(v_at1,0);
         v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                 _mm512_int2mask(13260),(__m512i)v_at1,78);
         v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                  _mm512_int2mask(21845),(__m512i)v_at1,177);
/*       fxyz[4*j] = exyz[4*j]*at1;     */
/*       fxyz[1+4*j] = exyz[1+4*j]*at1; */
/*       fxyz[2+4*j] = exyz[2+4*j]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*j]);
         v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
         _mm512_store_ps((float *)&fxyz[4*j],v_zt2);
/*       fxyz[4*(j+k1)] = exyz[4*(j+k1)]*at1;     */
/*       fxyz[1+4*(j+k1)] = exyz[1+4*(j+k1)]*at1; */
/*       fxyz[2+4*(j+k1)] = exyz[2+4*(j+k1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1)]);
         v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
         _mm512_store_ps((float *)&fxyz[4*(j+k1)],v_zt2);
/*       fxyz[4*(j+l1)] = exyz[4*(j+l1)]*at1;     */
/*       fxyz[1+4*(j+l1)] = exyz[1+4*(j+l1)]*at1; */
/*       fxyz[2+4*(j+l1)] = exyz[2+4*(j+l1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+l1)]);
         v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
         _mm512_store_ps((float *)&fxyz[4*(j+l1)],v_zt2);
/*       fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;     */
/*       fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1; */
/*       fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
         v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
         _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         at1 = cimagf(ffc[j]);
         fxyz[4*j] = exyz[4*j]*at1;
         fxyz[1+4*j] = exyz[1+4*j]*at1;
         fxyz[2+4*j] = exyz[2+4*j]*at1;
         fxyz[4*(j+k1)] = exyz[4*(j+k1)]*at1;
         fxyz[1+4*(j+k1)] = exyz[1+4*(j+k1)]*at1;
         fxyz[2+4*(j+k1)] = exyz[2+4*(j+k1)]*at1;
         fxyz[4*(j+l1)] = exyz[4*(j+l1)]*at1;
         fxyz[1+4*(j+l1)] = exyz[1+4*(j+l1)]*at1;
         fxyz[2+4*(j+l1)] = exyz[2+4*(j+l1)]*at1;
         fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;
         fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1;
         fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncfft3rmxy(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nzi, int nzp, int nxhd, int nyd, int nzd,
                  int nxhyzd, int nxyzhd) {
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
   requires KNC, f needs to be 64 byte aligned
   nxhd need to be a multiple of 8
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, nxhh, ny, nyh;
   int nz, nxyz, nxhyz, nzt, nrx, nry, nrxb, nryb, nxhyd;
   int i, j, k, l, n, nn, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
   int nss, nxhs, nxhhs, itn;
   float ani;
   float complex t1, t2, t3;
   __m512i v_j, v_kmr, v_m, v_n, v_it;
   __m512 v_zero, v_t1, v_t2, v_t3, v_t4, v_t5, v_ani;
   v_j = _mm512_set_epi32(7,7,6,6,5,5,4,4,3,3,2,2,1,1,0,0);
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
   nxhs = 8*(nxh/8);
   nxhhs = 8*(nxhh/8);
   itn = 1 > nxhhs ? 1 : nxhhs;
   v_m = _mm512_set_epi32(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0);
   v_n = _mm512_set_epi32(1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14);
   v_zero = _mm512_setzero_ps();
   v_t1 = _mm512_setzero_ps();
   v_t2 = _mm512_setzero_ps();
   v_t3 = _mm512_setzero_ps();
   v_t4 = _mm512_setzero_ps();
   if (isign > 0)
      goto L180;
/* inverse fourier transform */
   nrxb = nxhyz/nxh;
   nrx = nxyz/nxh;
   nryb = nxhyz/ny;
   nry = nxyz/ny;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,nss,km,kmr,k1,k2,j1,j2,nn,joff,ani,t1,t2,t3, \
v_it,v_kmr,v_t1,v_ani,v_t2,v_t3,v_t4,v_t5)
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
         nss = 8*(ns/8);
         v_kmr = _mm512_set1_epi32(2*kmr);
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (i = 0; i < ny; i++) {
               joff = nxhd*i + nn;
/* vector loop over elements in blocks of 8 */
               for (j = 0; j < nss; j+=8) {
/*                t1 = sct[kmr*j]; */
                  v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
                  v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
                  v_t1 = _mm512_i32gather_ps(v_it,(float *)sct,4);
/*                t2 = t1*f[j+k2+joff]; */
                  v_t2 = _mm512_load_ps((float *)&f[j+k2+joff]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[j+k2+joff] = f[j+k1+joff] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[j+k1+joff]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[j+k2+joff],v_t4);
/*                f[j+k1+joff] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[j+k1+joff],v_t4);
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
      kmr = nxyz/nx;
      ani = 0.5/(((float) nx)*((float) ny)*((float) nz));
      v_ani = _mm512_set1_ps(ani);
      v_kmr = _mm512_set1_epi32(2*kmr);
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
/* vector loop over elements in blocks of 8 */
         for (j = 0; j < nxhhs; j+=8) {
/*          t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
            v_t3 = _mm512_i32gather_ps(v_it,(float *)sct,4);
            v_t3 = _mm512_mask_sub_ps(v_t3,_mm512_int2mask(21845),
                   v_zero,v_t3);
            v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,177);
/*          t2 = conjf(f[nxh-j+joff]); */
            v_t2 = _mm512_loadunpacklo_ps(v_t2,
                   (float *)&f[nxh-j+joff-7]);
            v_t2 = _mm512_loadunpackhi_ps(v_t2,
                   (float *)&f[nxh-j+joff+1]);
/* reverse data */
            v_t2 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t2);
            v_t2 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(43690),
                   v_zero,v_t2);
/*          t1 = f[j+joff] + t2; */
            v_t4 = _mm512_load_ps((float *)&f[j+joff]);
            v_t1 = _mm512_add_ps(v_t4,v_t2);
/*          t2 = (f[j+joff] - t2)*t3; */
            v_t2 = _mm512_sub_ps(v_t4,v_t2);
            v_t5 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,160);
            v_t5 = _mm512_mul_ps(v_t2,v_t5);
            v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
            v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,245);
            v_t4 = _mm512_mul_ps(v_t2,v_t4);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                   v_zero,v_t4);
            v_t2 = _mm512_add_ps(v_t5,v_t4);
/*          f[j+joff] = ani*(t1 + t2); */
            v_t3 = _mm512_mul_ps(v_ani,_mm512_add_ps(v_t1,v_t2));
/*          f[nxh-j+joff] = ani*conjf(t1 - t2); */
            v_t4 = _mm512_sub_ps(v_t1,v_t2);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(43690),
                   v_zero,v_t4);
            v_t4 = _mm512_mul_ps(v_ani,v_t4);
/* reverse data */
            v_t4 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t4);
            if (j==0) {
               _mm512_mask_store_ps((float *)&f[j+joff],
                  _mm512_int2mask(65532),v_t3);
               _mm512_mask_packstorelo_ps((float *)&f[nxh-j+joff-7],
                  _mm512_int2mask(16383),v_t4);
               _mm512_mask_packstorehi_ps((float *)&f[nxh-j+joff+1],
                  _mm512_int2mask(16383),v_t4);
            }
            else {
               _mm512_store_ps((float *)&f[j+joff],v_t3);
               _mm512_packstorelo_ps((float *)&f[nxh-j+joff-7],v_t4);
               _mm512_packstorehi_ps((float *)&f[nxh-j+joff+1],v_t4);
            }
         }
/* loop over remaining elements */
         for (j = itn; j < nxhh; j++) {
            t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
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
/* vector loop over elements in blocks of 8 */
            for (i = 0; i < nxhs; i+=8) {
/*             t1 = f[i+k1]; */
               v_t1 = _mm512_load_ps((float *)&f[i+k1]);
/*             f[i+k1] = f[i+joff]; */
               v_t2 = _mm512_load_ps((float *)&f[i+joff]);
               _mm512_store_ps((float *)&f[i+k1],v_t2);
/*             f[i+joff] = t1; */
               _mm512_store_ps((float *)&f[i+joff],v_t1);
            }
/* loop over remaining elements */
            for (i = nxhs; i < nxh; i++) {
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
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                      crealf(t1));
/* vector loop over elements in blocks of 8 */
               for (i = 0; i < nxhs; i+=8) {
/*                t2 = t1*f[i+j2]; */
                  v_t2 = _mm512_load_ps((float *)&f[i+j2]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[i+j2] = f[i+j1] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[i+j1]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+j2],v_t4);
/*                f[i+j1] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+j1],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
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
private(i,j,k,l,n,ns,ns2,nss,km,kmr,k1,k2,j1,j2,nn,joff,t1,t2,t3,v_it, \
v_kmr,v_t1,v_t2,v_t3,v_t4,v_t5)
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
/* vector loop over elements in blocks of 8 */
            for (i = 0; i < nxhs; i+=8) {
/*             t1 = f[i+k1]; */
               v_t1 = _mm512_load_ps((float *)&f[i+k1]);
/*             f[i+k1] = f[i+joff]; */
               v_t2 = _mm512_load_ps((float *)&f[i+joff]);
               _mm512_store_ps((float *)&f[i+k1],v_t2);
/*             f[i+joff] = t1; */
               _mm512_store_ps((float *)&f[i+joff],v_t1);
            }
/* loop over remaining elements */
            for (i = nxhs; i < nxh; i++) {
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
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                      crealf(t1));
/* vector loop over elements in blocks of 8 */
               for (i = 0; i < nxhs; i+=8) {
/*                t2 = t1*f[i+j2]; */
                  v_t2 = _mm512_load_ps((float *)&f[i+j2]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[i+j2] = f[i+j1] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[i+j1]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+j2],v_t4);
/*                f[i+j1] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+j1],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
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
      v_kmr = _mm512_set1_epi32(2*kmr);
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
/* vector loop over elements in blocks of 8 */
         for (j = 0; j < nxhhs; j+=8) {
/*          t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
            v_t3 = _mm512_i32gather_ps(v_it,(float *)sct,4);
            v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,177);
/*          t2 = conjf(f[nxh-j+joff]); */
            v_t2 = _mm512_loadunpacklo_ps(v_t2,
                   (float *)&f[nxh-j+joff-7]);
            v_t2 = _mm512_loadunpackhi_ps(v_t2,
                   (float *)&f[nxh-j+joff+1]);
/* reverse data */
            v_t2 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t2);
            v_t2 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(43690),
                   v_zero,v_t2);
/*          t1 = f[j+joff] + t2; */
            v_t4 = _mm512_load_ps((float *)&f[j+joff]);
            v_t1 = _mm512_add_ps(v_t4,v_t2);
/*          t2 = (f[j+joff] - t2)*t3; */
            v_t2 = _mm512_sub_ps(v_t4,v_t2);
            v_t5 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,160);
            v_t5 = _mm512_mul_ps(v_t2,v_t5);
            v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
            v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,245);
            v_t4 = _mm512_mul_ps(v_t2,v_t4);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                   v_zero,v_t4);
            v_t2 = _mm512_add_ps(v_t5,v_t4);
/*          f[j+joff] = t1 + t2; */
            v_t3 = _mm512_add_ps(v_t1,v_t2);
/*          f[nxh-j+joff] = conjf(t1 - t2); */
            v_t4 = _mm512_sub_ps(v_t1,v_t2);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(43690),
                   v_zero,v_t4);
/* reverse data */
            v_t4 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t4);
            if (j==0) {
               _mm512_mask_store_ps((float *)&f[j+joff],
                  _mm512_int2mask(65532),v_t3);
               _mm512_mask_packstorelo_ps((float *)&f[nxh-j+joff-7],
                  _mm512_int2mask(16383),v_t4);
               _mm512_mask_packstorehi_ps((float *)&f[nxh-j+joff+1],
                  _mm512_int2mask(16383),v_t4);
            }
            else {
               _mm512_store_ps((float *)&f[j+joff],v_t3);
               _mm512_packstorelo_ps((float *)&f[nxh-j+joff-7],v_t4);
               _mm512_packstorehi_ps((float *)&f[nxh-j+joff+1],v_t4);
            }
         }
/* loop over remaining elements */
         for (j = itn; j < nxhh; j++) {
            t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
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
         nss = 8*(ns/8);
         v_kmr = _mm512_set1_epi32(2*kmr);
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (i = 0; i < ny; i++) {
               joff = nxhd*i + nn;
/* vector loop over elements in blocks of 8 */
               for (j = 0; j < nss; j+=8) {
/*                t1 = conjf(sct[kmr*j]); */
                  v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
                  v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
                  v_t1 = _mm512_i32gather_ps(v_it,(float *)sct,4);
                  v_t1 = _mm512_mask_sub_ps(v_t1,_mm512_int2mask(43690),
                         v_zero,v_t1);
/*                t2 = t1*f[j+k2+joff]; */
                  v_t2 = _mm512_load_ps((float *)&f[j+k2+joff]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[j+k2+joff] = f[j+k1+joff] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[j+k1+joff]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[j+k2+joff],v_t4);
/*                f[j+k1+joff] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[j+k1+joff],v_t4);
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
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncfft3rmz(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int nyi, int nyp, int nxhd, int nyd, int nzd,
                 int nxhyzd, int nxyzhd) {
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
   requires KNC, f needs to be 64 byte aligned
   nxhd need to be a multiple of 8
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, ny, nyh;
   int nz, nzh, nxyz, nxhyz, nyt, nrz, nrzb, nxhyd, ioff;
   int i, j, k, l, n, ll, j1, j2, k1, k2, l1, ns, ns2, km, kmr, i0, i1;
   int nss, nxhs;
   float complex t1, t2;
   __m512 v_zero, v_t1, v_t2, v_t3, v_t4;
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
   nxhs = 8*(nxh/8);
   v_zero = _mm512_setzero_ps();
   v_t1 = _mm512_setzero_ps();
   v_t2 = _mm512_setzero_ps();
   v_t3 = _mm512_setzero_ps();
   v_t4 = _mm512_setzero_ps();
   if (isign > 0)
      goto L90;
/* inverse fourier transform */
   nrzb = nxhyz/nz;
   nrz = nxyz/nz;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,ll,l1,i0,i1,ioff,t1,t2, \
v_t1,v_t2,v_t3,v_t4)
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
/* vector loop over elements in blocks of 8 */
            for (i = 0; i < nxhs; i+=8) {
/*             t1 = f[i+i1]; */
               v_t1 = _mm512_load_ps((float *)&f[i+i1]);
/*             f[i+i1] = f[i+i0]; */
               v_t2 = _mm512_load_ps((float *)&f[i+i0]);
               _mm512_store_ps((float *)&f[i+i1],v_t2);
/*             f[i+i0] = t1; */
               _mm512_store_ps((float *)&f[i+i0],v_t1);
            }
/* loop over remaining elements */
            for (i = nxhs; i < nxh; i++) {
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
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                   crealf(t1));
               i0 = ioff + j1;
               i1 = ioff + j2;
/* vector loop over elements in blocks of 8 */
               for (i = 0; i < nxhs; i+=8) {
/*                t2 = t1*f[i+i1]; */
                  v_t2 = _mm512_load_ps((float *)&f[i+i1]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[i+i1] = f[i+i0] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[i+i0]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+i1],v_t4);
/*                f[i+i0] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+i0],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
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
   if (nyi==1) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         t1 = f[l1];
         f[l1] = 0.5*(cimagf(f[ll] + t1)
                    + crealf(f[ll] - t1)*_Complex_I);
         f[ll] = 0.5*(crealf(f[ll] + t1)
                    + cimagf(f[ll] - t1)*_Complex_I);
      }
   }
   if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
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
   if (nyi==1) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         t1 = cimagf(f[l1]) + crealf(f[l1])*_Complex_I;
         f[l1] = conjf(f[ll] - t1);
         f[ll] += t1;
      }
   }
   if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
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
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,ll,l1,i0,i1,ioff,t1,t2, \
v_t1,v_t2,v_t3,v_t4)
   for (n = nyi-1; n < nyt; n++) {
      ioff = nxhd*n;
      for (l = 0; l < nz; l++) {
         ll = nxhyd*l;
         l1 = (mixup[l] - 1)/nrzb;
         if (l < l1) {
            l1 = nxhyd*l1;
            i0 = ioff + ll;
            i1 = ioff + l1;
/* vector loop over elements in blocks of 8 */
            for (i = 0; i < nxhs; i+=8) {
/*             t1 = f[i+i1]; */
               v_t1 = _mm512_load_ps((float *)&f[i+i1]);
/*             f[i+i1] = f[i+i0]; */
               v_t2 = _mm512_load_ps((float *)&f[i+i0]);
               _mm512_store_ps((float *)&f[i+i1],v_t2);
/*             f[i+i0] = t1; */
               _mm512_store_ps((float *)&f[i+i0],v_t1);
            }
/* loop over remaining elements */
            for (i = nxhs; i < nxh; i++) {
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
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                      crealf(t1));
               i0 = ioff + j1;
               i1 = ioff + j2;
/* vector loop over elements in blocks of 8 */
               for (i = 0; i < nxhs; i+=8) {
/*                t2 = t1*f[i+i1]; */
                  v_t2 = _mm512_load_ps((float *)&f[i+i1]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[i+i1] = f[i+i0] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[i+i0]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+i1],v_t4);
/*                f[i+i0] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+i0],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
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
void ckncfft3rm3xy(float complex f[], int isign, int mixup[],
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
   f[i][m][n][0:2] = (1/nx*ny*nz)*sum(f[i][k][j][0:2]*
         exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
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
   requires KNC, f needs to be 64 byte aligned
   nxhd need to be a multiple of 2
   f needs to have 4 components
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, nxhh, ny, nyh;
   int nz, nxyz, nxhyz, nzt, nrx, nry, nrxb, nryb, nxhd4, nxhyd;
   int i, j, k, l, n, nn, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
   int nss, nxhs, nxhhs, itn;
   float at1, at2, ani;
   float complex t1, t2, t3, t4;
   __m512i v_j, v_kmr, v_m, v_n, v_l, v_it;
   __m512 v_zero, v_t1, v_t2, v_t3, v_t4, v_t5, v_ani, v_half;
   v_j = _mm512_set_epi32(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0);
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
   nxhd4 = 4*nxhd;
   nxhyd = nxhd4*nyd;
   nxhs = 2*(nxh/2);
   nxhhs = 2*(nxhh/2);
   itn = 1 > nxhhs ? 1 : nxhhs;
   v_m = _mm512_set_epi32(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0);
   v_n = _mm512_set_epi32(7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8);
   v_zero = _mm512_setzero_ps();
   v_t1 = _mm512_setzero_ps();
   v_t2 = _mm512_setzero_ps();
   v_t3 = _mm512_setzero_ps();
   v_t4 = _mm512_setzero_ps();
   v_half = _mm512_set1_ps(0.5f);
   if (isign > 0)
      goto L230;
/* inverse fourier transform */
   nrxb = nxhyz/nxh;
   nrx = nxyz/nxh;
   nryb = nxhyz/ny;
   nry = nxyz/ny;
   v_l = _mm512_set_epi32(15,11,14,10,13,9,12,8,7,3,6,2,5,1,4,0);
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,nss,km,kmr,k1,k2,jj,j1,j2,nn,joff,at1,at2, \
ani,t1,t2,t3,t4,v_it,v_kmr,v_t1,v_ani,v_t2,v_t3,v_t4,v_t5)
   for (n = nzi-1; n < nzt; n++) {
      nn = nxhyd*n;
/* swap complex components */
      for (i = 0; i < ny; i++) {
         joff = nxhd4*i + nn;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          at1 = cimagf(f[2+4*j+joff]);                          */
/*          at2 = crealf(f[2+4*j+joff]);                          */
/*          f[2+4*j+joff] = crealf(f[1+4*j+joff])                 */
/*                          + crealf(f[3+4*j+joff])*_Complex_I;   */
/*          f[1+4*j+joff] = cimagf(f[4*j+joff]) + at1*_Complex_I; */
/*          f[4*j+joff] = crealf(f[4*j+joff]) + at2*_Complex_I;   */
            v_t1 = _mm512_load_ps((float *)&f[4*j+joff]);
            v_t1 = (__m512)_mm512_permutevar_epi32(v_l,(__m512i)v_t1);
            _mm512_store_ps((float *)&f[4*j+joff],v_t1);
          }
/* loop over remaining elements */
         for (j = nxhs; j < nxh; j++) {
            at1 = cimagf(f[2+4*j+joff]);
            at2 = crealf(f[2+4*j+joff]);
            f[2+4*j+joff] = crealf(f[1+4*j+joff])
                            + crealf(f[3+4*j+joff])*_Complex_I;
            f[1+4*j+joff] = cimagf(f[4*j+joff]) + at1*_Complex_I;
            f[4*j+joff] = crealf(f[4*j+joff]) + at2*_Complex_I;
          }
      }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            for (i = 0; i < ny; i++) {
               joff = nxhd4*i + nn;
/*             t1 = f[4*j1+joff];   */
/*             t2 = f[1+4*j1+joff]; */
/*             t3 = f[2+4*j1+joff]; */
               v_t1 = _mm512_mask_loadunpacklo_ps(v_t1,
                      _mm512_int2mask(255),(float *)&f[4*j1+joff]);
               v_t1 = _mm512_mask_loadunpackhi_ps(v_t1,
                      _mm512_int2mask(255),(float *)&f[4*j1+joff+8]);
/*             f[4*j1+joff] = f[4*j+joff];     */
/*             f[1+4*j1+joff] = f[1+4*j+joff]; */
/*             f[2+4*j1+joff] = f[2+4*j+joff]; */
               v_t2 = _mm512_mask_loadunpacklo_ps(v_t2,
                      _mm512_int2mask(255),(float *)&f[4*j+joff]);
               v_t2 = _mm512_mask_loadunpackhi_ps(v_t2,
                      _mm512_int2mask(255),(float *)&f[4*j+joff+8]);
               _mm512_mask_packstorelo_ps((float *)&f[4*j1+joff],
                  _mm512_int2mask(255),v_t2);
               _mm512_mask_packstorehi_ps((float *)&f[4*j1+joff+8],
                  _mm512_int2mask(255),v_t2);
/*             f[4*j+joff] = t1;   */
/*             f[1+4*j+joff] = t2; */
/*             f[2+4*j+joff] = t3; */
               _mm512_mask_packstorelo_ps((float *)&f[4*j+joff],
                  _mm512_int2mask(255),v_t1);
               _mm512_mask_packstorehi_ps((float *)&f[4*j+joff+8],
                  _mm512_int2mask(255),v_t1);
            }
         }
      }
/* first transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         nss = 2*(ns/2);
         v_kmr = _mm512_set1_epi32(2*kmr);
         for (k = 0; k < km; k++) {
            k1 = 4*ns2*k;
            k2 = k1 + 4*ns;
            for (i = 0; i < ny; i++) {
               joff = nxhd4*i + nn;
/* vector loop over elements in blocks of 2 */
               for (j = 0; j < nss; j+=2) {
/*                t1 = sct[kmr*j]; */
                  v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
                  v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
                  v_t1 = _mm512_i32gather_ps(v_it,(float *)sct,4);
/*                t2 = t1*f[4*j+k2+joff];   */
/*                t3 = t1*f[1+4*j+k2+joff]; */
/*                t4 = t1*f[2+4*j+k2+joff]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*j+k2+joff]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*j+k2+joff] = f[4*j+k1+joff] - t2;     */
/*                f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3; */
/*                f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*j+k1+joff]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*j+k2+joff],v_t4);
/*                f[4*j+k1+joff] += t2;   */
/*                f[1+4*j+k1+joff] += t3; */
/*                f[2+4*j+k1+joff] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*j+k1+joff],v_t4);
               }
/* loop over remaining elements */
               for (j = nss; j < ns; j++) {
                  t1 = sct[kmr*j];
                  t2 = t1*f[4*j+k2+joff];
                  t3 = t1*f[1+4*j+k2+joff];
                  t4 = t1*f[2+4*j+k2+joff];
                  f[4*j+k2+joff] = f[4*j+k1+joff] - t2;
                  f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3;
                  f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4;
                  f[4*j+k1+joff] += t2;
                  f[1+4*j+k1+joff] += t3;
                  f[2+4*j+k1+joff] += t4;
               }
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxyz/nx;
      ani = 0.5/(((float) nx)*((float) ny)*((float) nz));
      v_ani = _mm512_set1_ps(ani);
      v_kmr = _mm512_set1_epi32(2*kmr);
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhhs; j+=2) {
/*          t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
            v_t3 = _mm512_i32gather_ps(v_it,(float *)sct,4);
            v_t3 = _mm512_mask_sub_ps(v_t3,_mm512_int2mask(21845),
                   v_zero,v_t3);
            v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,177);
/*          for (jj = 0; jj < 3; jj++) {      */
/*          t2 = conjf(f[jj+4*(nxh-j)+joff]); */
            v_t2 = _mm512_loadunpacklo_ps(v_t2,
                   (float *)&f[4*(nxh-j-1)+joff]);
            v_t2 = _mm512_loadunpackhi_ps(v_t2,
                   (float *)&f[4*(nxh-j-1)+joff+8]);
/* reverse data */
            v_t2 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t2);
            v_t2 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(43690),
                   v_zero,v_t2);
/*          t1 = f[jj+4*j+joff] + t2; */
            v_t4 = _mm512_load_ps((float *)&f[4*j+joff]);
            v_t1 = _mm512_add_ps(v_t4,v_t2);
/*          t2 = (f[jj+4*j+joff] - t2)*t3; */
            v_t2 = _mm512_sub_ps(v_t4,v_t2);
            v_t5 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,160);
            v_t5 = _mm512_mul_ps(v_t2,v_t5);
            v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
            v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,245);
            v_t4 = _mm512_mul_ps(v_t2,v_t4);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                   v_zero,v_t4);
            v_t2 = _mm512_add_ps(v_t5,v_t4);

/*          f[jj+4*j+joff] = ani*(t1 + t2); */
            v_t3 = _mm512_mul_ps(v_ani,_mm512_add_ps(v_t1,v_t2));
/*          f[jj+4*(nxh-j)+joff] = ani*conjf(t1 - t2); */
/*          }                                          */
            v_t4 = _mm512_sub_ps(v_t1,v_t2);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(43690),
                   v_zero,v_t4);
            v_t4 = _mm512_mul_ps(v_ani,v_t4);
/* reverse data */
            v_t4 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t4);
            if (j==0) {
               _mm512_mask_store_ps((float *)&f[4*j+joff],
                  _mm512_int2mask(65280),v_t3);
               _mm512_mask_packstorelo_ps((float *)&f[4*(nxh-j-1)+joff],
                  _mm512_int2mask(255),v_t4);
               _mm512_mask_packstorehi_ps((float *)&f[4*(nxh-j-1)+joff+8],
                  _mm512_int2mask(255),v_t4);
            }
            else {
               _mm512_store_ps((float *)&f[4*j+joff],v_t3);
               _mm512_packstorelo_ps((float *)&f[4*(nxh-j-1)+joff],v_t4);
               _mm512_packstorehi_ps((float *)&f[4*(nxh-j-1)+joff+8],v_t4);
            }
         }
/* loop over remaining elements */
         for (j = itn; j < nxhh; j++) {
            t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
            for (jj = 0; jj < 3; jj++) {
               t2 = conjf(f[jj+4*(nxh-j)+joff]);
               t1 = f[jj+4*j+joff] + t2;
               t2 = (f[jj+4*j+joff] - t2)*t3;
               f[jj+4*j+joff] = ani*(t1 + t2);
               f[jj+4*(nxh-j)+joff] = ani*conjf(t1 - t2);
            }
         }
      }
/*    ani = 2.0*ani; */
      v_ani = _mm512_add_ps(v_ani,v_ani);
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
/*       for (jj = 0; jj < 3; jj++) {                      */
/*       f[jj+4*nxhh+joff] = ani*conjf(f[jj+4*nxhh+joff]); */
         v_t1 = _mm512_mask_load_ps(v_t1,_mm512_int2mask(63),
                (float *)&f[4*nxhh+joff]);
         v_t1 = _mm512_mask_sub_ps(v_t1,_mm512_int2mask(42),v_zero,
                v_t1);
         v_t1 = _mm512_mul_ps(v_ani,v_t1);
         _mm512_mask_store_ps((float *)&f[4*nxhh+joff],
               _mm512_int2mask(63),v_t1);
/*       f[jj+joff] = ani*((crealf(f[jj+joff])            */
/*                     + cimagf(f[jj+joff]))              */
/*                     + (crealf(f[jj+joff])              */
/*                     - cimagf(f[jj+joff]))*_Complex_I); */
/*       }                                                */
         v_t2 = _mm512_mask_load_ps(v_t2,_mm512_int2mask(63),
                (float *)&f[joff]);
         v_t1 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
         v_t3 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(42),v_t1,v_t2);
         v_t3 = _mm512_mask_add_ps(v_t3,_mm512_int2mask(21),v_t1,v_t2);
         v_t3 = _mm512_mul_ps(v_ani,v_t3);
         _mm512_mask_store_ps((float *)&f[joff],_mm512_int2mask(63),
         v_t3);
      }
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nxhd4*k1 + nn;
/* vector loop over elements in blocks of 2 */
            for (i = 0; i < nxhs; i+=2) {
/*             t1 = f[4*i+k1];   */
/*             t2 = f[1+4*i+k1]; */
/*             t3 = f[2+4*i+k1]; */
               v_t1 = _mm512_load_ps((float *)&f[4*i+k1]);
/*             f[4*i+k1] = f[4*i+joff];     */
/*             f[1+4*i+k1] = f[1+4*i+joff]; */
/*             f[2+4*i+k1] = f[2+4*i+joff]; */
               v_t2 = _mm512_load_ps((float *)&f[4*i+joff]);
               _mm512_store_ps((float *)&f[4*i+k1],v_t2);
/*             f[4*i+joff] = t1;   */
/*             f[1+4*i+joff] = t2; */
/*             f[2+4*i+joff] = t3; */
               _mm512_store_ps((float *)&f[4*i+joff],v_t1);
            }
/* loop over remaining elements */
            for (i = nxhs; i < nxh; i++) {
               t1 = f[4*i+k1];
               t2 = f[1+4*i+k1];
               t3 = f[2+4*i+k1];
               f[4*i+k1] = f[4*i+joff];
               f[1+4*i+k1] = f[1+4*i+joff];
               f[2+4*i+k1] = f[2+4*i+joff];
               f[4*i+joff] = t1;
               f[1+4*i+joff] = t2;
               f[2+4*i+joff] = t3;
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
               j1 = nxhd4*(j + k1) + nn;
               j2 = nxhd4*(j + k2) + nn;
               t1 = sct[kmr*j];
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                      crealf(t1));
/* vector loop over elements in blocks of 2 */
               for (i = 0; i < nxhs; i+=2) {
/*                t2 = t1*f[4*i+j2];   */
/*                t3 = t1*f[1+4*i+j2]; */
/*                t4 = t1*f[2+4*i+j2]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*i+j2]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*i+j2] = f[4*i+j1] - t2;     */
/*                f[1+4*i+j2] = f[1+4*i+j1] - t3; */
/*                f[2+4*i+j2] = f[2+4*i+j1] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*i+j1]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+j2],v_t4);
/*                f[4*i+j1] += t2;   */
/*                f[1+4*i+j1] += t3; */
/*                f[2+4*i+j1] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+j1],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[4*i+j2];
                  t3 = t1*f[1+4*i+j2];
                  t4 = t1*f[2+4*i+j2];
                  f[4*i+j2] = f[4*i+j1] - t2;
                  f[1+4*i+j2] = f[1+4*i+j1] - t3;
                  f[2+4*i+j2] = f[2+4*i+j1] - t4;
                  f[4*i+j1] += t2;
                  f[1+4*i+j1] += t3;
                  f[2+4*i+j1] += t4;
               }
            }
         }
         ns = ns2;
      }
/* unscramble modes kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         joff = nxhd4*k;
         k1 = nxhd4*ny - joff + nn;
         joff += nn;
/*       for (jj = 0; jj < 3; jj++) { */
/*       t1 = f[jj+k1];               */
         v_t1 = _mm512_mask_load_ps(v_t1,_mm512_int2mask(63),
                (float *)&f[k1]);
/*       f[jj+k1] = 0.5*(cimagf(f[jj+joff] + t1)            */
/*                   + crealf(f[jj+joff] - t1)*_Complex_I); */
         v_t2 = _mm512_mask_load_ps(v_t2,_mm512_int2mask(63),
                (float *)&f[joff]);
         v_t3 = _mm512_mask_add_ps(v_t3,_mm512_int2mask(42),v_t2,v_t1);
         v_t3 = _mm512_mask_sub_ps(v_t3,_mm512_int2mask(21),v_t2,v_t1);
         v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,177);
         v_t3 = _mm512_mul_ps(v_half,v_t3);
         _mm512_mask_store_ps((float *)&f[k1],_mm512_int2mask(63),v_t3);
/*       f[jj+joff] = 0.5*(crealf(f[jj+joff] + t1)            */
/*                     + cimagf(f[jj+joff] - t1)*_Complex_I); */
/*       }                                                    */
         v_t2 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(42),v_t2,v_t1);
         v_t2 = _mm512_mask_add_ps(v_t2,_mm512_int2mask(21),v_t2,v_t1);
         v_t2 = _mm512_mul_ps(v_half,v_t2);
         _mm512_mask_store_ps((float *)&f[joff],_mm512_int2mask(63),v_t2);
      }
   }
   return;
/* forward fourier transform */
L230: nryb = nxhyz/ny;
   nry = nxyz/ny;
   nrxb = nxhyz/nxh;
   nrx = nxyz/nxh;
   v_l = _mm512_set_epi32(15,13,11,9,14,12,10,8,7,5,3,1,6,4,2,0);
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,nss,km,kmr,k1,k2,jj,j1,j2,nn,joff,at1,at2, \
t1,t2,t3,t4,v_it,v_kmr,v_t1,v_t2,v_t3,v_t4,v_t5)
   for (n = nzi-1; n < nzt; n++) {
      nn = nxhyd*n;
/* scramble modes kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         joff = nxhd4*k;
         k1 = nxhd4*ny - joff + nn;
         joff += nn;
/*       for (jj = 0; jj < 3; jj++) {                            */
/*       t1 = cimagf(f[jj+k1]) + crealf(f[jj+k1])*_Complex_I; */
         v_t1 = _mm512_mask_load_ps(v_t1,_mm512_int2mask(63),
                (float *)&f[k1]);
         v_t1 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,177);
/*       f[jj+k1] = conjf(f[jj+joff] - t1); */
         v_t2 = _mm512_mask_load_ps(v_t2,_mm512_int2mask(63),
                (float *)&f[joff]);
         v_t3 = _mm512_mask_sub_ps(v_t3,_mm512_int2mask(63),v_t2,v_t1);
         v_t3 = _mm512_mask_sub_ps(v_t3,_mm512_int2mask(42),
                v_zero,v_t3);
         _mm512_mask_store_ps((float *)&f[k1],_mm512_int2mask(63),v_t3);
/*       f[jj+joff] += t1; */
/*       }                 */
         v_t2 = _mm512_mask_add_ps(v_t2,_mm512_int2mask(63),v_t2,v_t1);
         _mm512_mask_store_ps((float *)&f[joff],_mm512_int2mask(63),
         v_t2);
      }
/* bit-reverse array elements in y */
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
         k1 = (mixup[k] - 1)/nryb;
         if (k < k1) {
            k1 = nxhd4*k1 + nn;
/* vector loop over elements in blocks of 2 */
            for (i = 0; i < nxhs; i+=2) {
/*             t1 = f[4*i+k1];   */
/*             t2 = f[1+4*i+k1]; */
/*             t3 = f[2+4*i+k1]; */
               v_t1 = _mm512_load_ps((float *)&f[4*i+k1]);
/*             f[4*i+k1] = f[4*i+joff];     */
/*             f[1+4*i+k1] = f[1+4*i+joff]; */
/*             f[2+4*i+k1] = f[2+4*i+joff]; */
               v_t2 = _mm512_load_ps((float *)&f[4*i+joff]);
               _mm512_store_ps((float *)&f[4*i+k1],v_t2);
/*             f[4*i+joff] = t1;   */
/*             f[1+4*i+joff] = t2; */
/*             f[2+4*i+joff] = t3; */
               _mm512_store_ps((float *)&f[4*i+joff],v_t1);
            }
/* loop over remaining elements */
            for (i = nxhs; i < nxh; i++) {
               t1 = f[4*i+k1];
               t2 = f[1+4*i+k1];
               t3 = f[2+4*i+k1];
               f[4*i+k1] = f[4*i+joff];
               f[1+4*i+k1] = f[1+4*i+joff];
               f[2+4*i+k1] = f[2+4*i+joff];
               f[4*i+joff] = t1;
               f[1+4*i+joff] = t2;
               f[2+4*i+joff] = t3;
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
               j1 = nxhd4*(j + k1) + nn;
               j2 = nxhd4*(j + k2) + nn;
               t1 = conjf(sct[kmr*j]);
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                      crealf(t1));
/* vector loop over elements in blocks of 2 */
               for (i = 0; i < nxhs; i+=2) {
/*                t2 = t1*f[4*i+j2];   */
/*                t3 = t1*f[1+4*i+j2]; */
/*                t4 = t1*f[2+4*i+j2]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*i+j2]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*i+j2] = f[4*i+j1] - t2;     */
/*                f[1+4*i+j2] = f[1+4*i+j1] - t3; */
/*                f[2+4*i+j2] = f[2+4*i+j1] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*i+j1]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+j2],v_t4);
/*                f[4*i+j1] += t2;   */
/*                f[1+4*i+j1] += t3; */
/*                f[2+4*i+j1] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+j1],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[4*i+j2];
                  t3 = t1*f[1+4*i+j2];
                  t4 = t1*f[2+4*i+j2];
                  f[4*i+j2] = f[4*i+j1] - t2;
                  f[1+4*i+j2] = f[1+4*i+j1] - t3;
                  f[2+4*i+j2] = f[2+4*i+j1] - t4;
                  f[4*i+j1] += t2;
                  f[1+4*i+j1] += t3;
                  f[2+4*i+j1] += t4;
               }
            }
         }
         ns = ns2;
      }
/* scramble coefficients */
      kmr = nxyz/nx;
      v_kmr = _mm512_set1_epi32(2*kmr);
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhhs; j+=2) {
/*          t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
            v_t3 = _mm512_i32gather_ps(v_it,(float *)sct,4);
            v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,177);
/*          for (jj = 0; jj < 3; jj++) {      */
/*          t2 = conjf(f[jj+4*(nxh-j)+joff]); */
            v_t2 = _mm512_loadunpacklo_ps(v_t2,
                   (float *)&f[4*(nxh-j-1)+joff]);
            v_t2 = _mm512_loadunpackhi_ps(v_t2,
                   (float *)&f[4*(nxh-j-1)+joff+8]);
/* reverse data */
            v_t2 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t2);
            v_t2 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(43690),
                   v_zero,v_t2);
/*          t1 = f[jj+4*j+joff] + t2; */
            v_t4 = _mm512_load_ps((float *)&f[4*j+joff]);
            v_t1 = _mm512_add_ps(v_t4,v_t2);
/*          t2 = (f[jj+4*j+joff] - t2)*t3; */
            v_t2 = _mm512_sub_ps(v_t4,v_t2);
            v_t5 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,160);
            v_t5 = _mm512_mul_ps(v_t2,v_t5);
            v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
            v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,245);
            v_t4 = _mm512_mul_ps(v_t2,v_t4);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                   v_zero,v_t4);
            v_t2 = _mm512_add_ps(v_t5,v_t4);
/*          f[jj+4*j+joff] = t1 + t2; */
            v_t3 = _mm512_add_ps(v_t1,v_t2);
/*          f[jj+4*(nxh-j)+joff] = conjf(t1 - t2); */
/*          }                                      */
            v_t4 = _mm512_sub_ps(v_t1,v_t2);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(43690),
                   v_zero,v_t4);
/* reverse data */
            v_t4 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t4);
            if (j==0) {
               _mm512_mask_store_ps((float *)&f[4*j+joff],
                  _mm512_int2mask(65280),v_t3);
               _mm512_mask_packstorelo_ps((float *)&f[4*(nxh-j-1)+joff],
                  _mm512_int2mask(255),v_t4);
               _mm512_mask_packstorehi_ps((float *)&f[4*(nxh-j-1)+joff+8],
                  _mm512_int2mask(255),v_t4);
            }
            else {
               _mm512_store_ps((float *)&f[4*j+joff],v_t3);
               _mm512_packstorelo_ps((float *)&f[4*(nxh-j-1)+joff],v_t4);
               _mm512_packstorehi_ps((float *)&f[4*(nxh-j-1)+joff+8],v_t4);
            }
         }
/* loop over remaining elements */
         for (j = itn; j < nxhh; j++) {
            t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
            for (jj = 0; jj < 3; jj++) {
               t2 = conjf(f[jj+4*(nxh-j)+joff]);
               t1 = f[jj+4*j+joff] + t2;
               t2 = (f[jj+4*j+joff] - t2)*t3;
               f[jj+4*j+joff] = t1 + t2;
               f[jj+4*(nxh-j)+joff] = conjf(t1 - t2);
            }
         }
      }
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
/*       for (jj = 0; jj < 3; jj++) { */
/*       f[jj+4*nxhh+joff] = 2.0*conjf(f[jj+4*nxhh+joff]); */
         v_t1 = _mm512_mask_load_ps(v_t1,_mm512_int2mask(63),
                (float *)&f[4*nxhh+joff]);
         v_t1 = _mm512_mask_sub_ps(v_t1,_mm512_int2mask(42),v_zero,
                v_t1);
         v_t1 = _mm512_add_ps(v_t1,v_t1);
         _mm512_mask_store_ps((float *)&f[4*nxhh+joff],
               _mm512_int2mask(63),v_t1);
/*       f[jj+joff] = (crealf(f[jj+joff]) + cimagf(f[jj+joff])) */
/*                  + (crealf(f[jj+joff])                       */
/*                  - cimagf(f[jj+joff]))*_Complex_I;           */
/*       }                                                      */
         v_t2 = _mm512_mask_load_ps(v_t2,_mm512_int2mask(63),
                (float *)&f[joff]);
         v_t1 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
         v_t3 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(42),v_t1,v_t2);
         v_t3 = _mm512_mask_add_ps(v_t3,_mm512_int2mask(21),v_t1,v_t2);
         _mm512_mask_store_ps((float *)&f[joff],_mm512_int2mask(63),
         v_t3);
      }
/* bit-reverse array elements in x */
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrxb;
         if (j < j1) {
            for (i = 0; i < ny; i++) {
               joff = nxhd4*i + nn;
/*             t1 = f[4*j1+joff];   */
/*             t2 = f[1+4*j1+joff]; */
/*             t3 = f[2+4*j1+joff]; */
               v_t1 = _mm512_mask_loadunpacklo_ps(v_t1,
                      _mm512_int2mask(255),(float *)&f[4*j1+joff]);
               v_t1 = _mm512_mask_loadunpackhi_ps(v_t1,
                      _mm512_int2mask(255),(float *)&f[4*j1+joff+8]);
/*             f[4*j1+joff] = f[4*j+joff];     */
/*             f[1+4*j1+joff] = f[1+4*j+joff]; */
/*             f[2+4*j1+joff] = f[2+4*j+joff]; */
               v_t2 = _mm512_mask_loadunpacklo_ps(v_t2,
                      _mm512_int2mask(255),(float *)&f[4*j+joff]);
               v_t2 = _mm512_mask_loadunpackhi_ps(v_t2,
                      _mm512_int2mask(255),(float *)&f[4*j+joff+8]);
               _mm512_mask_packstorelo_ps((float *)&f[4*j1+joff],
                  _mm512_int2mask(255),v_t2);
               _mm512_mask_packstorehi_ps((float *)&f[4*j1+joff+8],
                  _mm512_int2mask(255),v_t2);
/*             f[4*j+joff] = t1;   */
/*             f[1+4*j+joff] = t2; */
/*             f[2+4*j+joff] = t3; */
               _mm512_mask_packstorelo_ps((float *)&f[4*j+joff],
                  _mm512_int2mask(255),v_t1);
               _mm512_mask_packstorehi_ps((float *)&f[4*j+joff+8],
                  _mm512_int2mask(255),v_t1);
            }
         }
      }
/* finally transform in x */
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         nss = 2*(ns/2);
         v_kmr = _mm512_set1_epi32(2*kmr);
         for (k = 0; k < km; k++) {
            k1 = 4*ns2*k;
            k2 = k1 + 4*ns;
            for (i = 0; i < ny; i++) {
               joff = nxhd4*i + nn;
/* vector loop over elements in blocks of 2 */
               for (j = 0; j < nss; j+=2) {
/*                t1 = conjf(sct[kmr*j]); */
                  v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
                  v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
                  v_t1 = _mm512_i32gather_ps(v_it,(float *)sct,4);
                  v_t1 = _mm512_mask_sub_ps(v_t1,_mm512_int2mask(43690),
                         v_zero,v_t1);
/*                t2 = t1*f[4*j+k2+joff];   */
/*                t3 = t1*f[1+4*j+k2+joff]; */
/*                t4 = t1*f[2+4*j+k2+joff]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*j+k2+joff]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*j+k2+joff] = f[4*j+k1+joff] - t2;     */
/*                f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3; */
/*                f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*j+k1+joff]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*j+k2+joff],v_t4);
/*                f[4*j+k1+joff] += t2;   */
/*                f[1+4*j+k1+joff] += t3; */
/*                f[2+4*j+k1+joff] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*j+k1+joff],v_t4);
               }
/* loop over remaining elements */
               for (j = nss; j < ns; j++) {
                  t1 = conjf(sct[kmr*j]);
                  t2 = t1*f[4*j+k2+joff];
                  t3 = t1*f[1+4*j+k2+joff];
                  t4 = t1*f[2+4*j+k2+joff];
                  f[4*j+k2+joff] = f[4*j+k1+joff] - t2;
                  f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3;
                  f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4;
                  f[4*j+k1+joff] += t2;
                  f[1+4*j+k1+joff] += t3;
                  f[2+4*j+k1+joff] += t4;
               }
            }
         }
         ns = ns2;
      }
/* swap complex components */
      for (i = 0; i < ny; i++) {
         joff = nxhd4*i + nn;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          f[3+4*j+joff] = cimagf(f[2+4*j+joff]) */
/*                          + cimagf(f[3+4*j+joff])*_Complex_I; */
/*          at1 = crealf(f[2+4*j+joff]); */
/*          f[2+4*j+joff] = cimagf(f[4*j+joff]) */
/*                          + cimagf(f[1+4*j+joff])*_Complex_I; */
/*          at2 = crealf(f[1+4*j+joff]); */
/*          f[1+4*j+joff] = at1 + 0.0*_Complex_I; */
/*          f[4*j+joff] = crealf(f[4*j+joff]) + at2*_Complex_I; */
            v_t1 = _mm512_load_ps((float *)&f[4*j+joff]);
            v_t1 = (__m512)_mm512_permutevar_epi32(v_l,(__m512i)v_t1);
            _mm512_store_ps((float *)&f[4*j+joff],v_t1);
         }
/* loop over remaining elements */
         for (j = nxhs; j < nxh; j++) {
            f[3+4*j+joff] = cimagf(f[2+4*j+joff])
                            + cimagf(f[3+4*j+joff])*_Complex_I;
            at1 = crealf(f[2+4*j+joff]);
            f[2+4*j+joff] = cimagf(f[4*j+joff])
                            + cimagf(f[1+4*j+joff])*_Complex_I;
            at2 = crealf(f[1+4*j+joff]);
            f[1+4*j+joff] = at1 + 0.0*_Complex_I;
            f[4*j+joff] = crealf(f[4*j+joff]) + at2*_Complex_I;
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncfft3rm3z(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nyi, int nyp, int nxhd, int nyd, int nzd,
                  int nxhyzd, int nxyzhd) {
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
   requires KNC, f needs to be 64 byte aligned
   nxhd need to be a multiple of 2
   f needs to have 4 components
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, ny, nyh;
   int nz, nzh, nxyz, nxhyz, nyt, nrz, nrzb, nxhd4, nxhyd, ioff;
   int i, j, k, l, n, ll, jj, j1, j2, k1, k2, l1, ns, ns2, km, kmr;
   int i0, i1;
   int nxhs;
   float complex t1, t2, t3, t4;
   __m512 v_zero, v_t1, v_t2, v_t3, v_t4;
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
   nxhd4 = 4*nxhd;
   nxhyd = nxhd4*nyd;
   nxhs = 2*(nxh/2);
   v_zero = _mm512_setzero_ps();
   v_t1 = _mm512_setzero_ps();
   v_t2 = _mm512_setzero_ps();
   v_t3 = _mm512_setzero_ps();
   v_t4 = _mm512_setzero_ps();
   if (isign > 0)
      goto L110;
/* inverse fourier transform */
   nrzb = nxhyz/nz;
   nrz = nxyz/nz;
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,ll,l1,i0,i1,ioff,t1,t2,t3, \
t4,v_t1,v_t2,v_t3,v_t4)
   for (n = nyi-1; n < nyt; n++) {
      ioff = nxhd4*n;
/* bit-reverse array elements in z */
      for (l = 0; l < nz; l++) {
         ll = nxhyd*l;
         l1 = (mixup[l] - 1)/nrzb;
         if (l < l1) {
            l1 = nxhyd*l1;
            i0 = ioff + ll;
            i1 = ioff + l1;
/* vector loop over elements in blocks of 2 */
            for (i = 0; i < nxhs; i+=2) {
/*             t1 = f[4*i+i1];   */
/*             t2 = f[1+4*i+i1]; */
/*             t3 = f[2+4*i+i1]; */
               v_t1 = _mm512_load_ps((float *)&f[4*i+i1]);
/*             f[4*i+i1] = f[4*i+i0];     */
/*             f[1+4*i+i1] = f[1+4*i+i0]; */
/*             f[2+4*i+i1] = f[2+4*i+i0]; */
               v_t2 = _mm512_load_ps((float *)&f[4*i+i0]);
               _mm512_store_ps((float *)&f[4*i+i1],v_t2);
/*             f[4*i+i0] = t1;   */
/*             f[1+4*i+i0] = t2; */
/*             f[2+4*i+i0] = t3; */
               _mm512_store_ps((float *)&f[4*i+i0],v_t1);
            }
/* loop over remaining elements */
            for (i = nxhs; i < nxh; i++) {
               t1 = f[4*i+i1];
               t2 = f[1+4*i+i1];
               t3 = f[2+4*i+i1];
               f[4*i+i1] = f[4*i+i0];
               f[1+4*i+i1] = f[1+4*i+i0];
               f[2+4*i+i1] = f[2+4*i+i0];
               f[4*i+i0] = t1;
               f[1+4*i+i0] = t2;
               f[2+4*i+i0] = t3;
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
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                      crealf(t1));
               i0 = ioff + j1;
               i1 = ioff + j2;
/* vector loop over elements in blocks of 2 */
               for (i = 0; i < nxhs; i+=2) {
/*                t2 = t1*f[4*i+i1];   */
/*                t3 = t1*f[1+4*i+i1]; */
/*                t4 = t1*f[2+4*i+i1]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*i+i1]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*i+i1] = f[4*i+i0] - t2;     */
/*                f[1+4*i+i1] = f[1+4*i+i0] - t3; */
/*                f[2+4*i+i1] = f[2+4*i+i0] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*i+i0]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+i1],v_t4);
/*                f[4*i+i0] += t2;   */
/*                f[1+4*i+i0] += t3; */
/*                f[2+4*i+i0] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+i0],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[4*i+i1];
                  t3 = t1*f[1+4*i+i1];
                  t4 = t1*f[2+4*i+i1];
                  f[4*i+i1] = f[4*i+i0] - t2;
                  f[1+4*i+i1] = f[1+4*i+i0] - t3;
                  f[2+4*i+i1] = f[2+4*i+i0] - t4;
                  f[4*i+i0] += t2;
                  f[1+4*i+i0] += t3;
                  f[2+4*i+i0] += t4;
               }
            }
         }
         ns = ns2;
      }
   }
/* unscramble modes kx = 0, nx/2 */
   if (nyi==1) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         for (jj = 0; jj < 3; jj++) {
            t1 = f[jj+l1];
            f[jj+l1] = 0.5*(cimagf(f[jj+ll] + t1)
                          + crealf(f[jj+ll] - t1)*_Complex_I);
            f[jj+ll] = 0.5*(crealf(f[jj+ll] + t1)
                          + cimagf(f[jj+ll] - t1)*_Complex_I);
         }
      }
   }
   if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         i1 = nxhd4*nyh;
         i0 = i1 + ll;
         i1 += l1;
         for (jj = 0; jj < 3; jj++) {
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
   if (nyi==1) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         for (jj = 0; jj < 3; jj++) {
            t1 = cimagf(f[jj+l1]) + crealf(f[jj+l1])*_Complex_I;
            f[jj+l1] = conjf(f[jj+ll] - t1);
            f[jj+ll] += t1;
         }
      }
   }
   if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
      for (n = 1; n < nzh; n++) {
         ll = nxhyd*n;
         l1 = nxhyd*nz - ll;
         i1 = nxhd4*nyh;
         i0 = i1 + ll;
         i1 += l1;
         for (jj = 0; jj < 3; jj++) {
            t1 = cimagf(f[jj+i1]) + crealf(f[jj+i1])*_Complex_I;
            f[jj+i1] = conjf(f[jj+i0] - t1);
            f[jj+i0] += t1;
         }
      }
   }
#pragma omp parallel for \
private(i,j,k,l,n,ns,ns2,km,kmr,k1,k2,j1,j2,ll,l1,i0,i1,ioff,t1,t2,t3, \
t4,v_t1,v_t2,v_t3,v_t4)
   for (n = nyi-1; n < nyt; n++) {
      ioff = nxhd4*n;
/* bit-reverse array elements in z */
      for (l = 0; l < nz; l++) {
         ll = nxhyd*l;
         l1 = (mixup[l] - 1)/nrzb;
         if (l < l1) {
            l1 = nxhyd*l1;
            i0 = ioff + ll;
            i1 = ioff + l1;
/* vector loop over elements in blocks of 2 */
            for (i = 0; i < nxhs; i+=2) {
/*             t1 = f[4*i+i1];   */
/*             t2 = f[1+4*i+i1]; */
/*             t3 = f[2+4*i+i1]; */
               v_t1 = _mm512_load_ps((float *)&f[4*i+i1]);
/*             f[4*i+i1] = f[4*i+i0];     */
/*             f[1+4*i+i1] = f[1+4*i+i0]; */
/*             f[2+4*i+i1] = f[2+4*i+i0]; */
               v_t2 = _mm512_load_ps((float *)&f[4*i+i0]);
               _mm512_store_ps((float *)&f[4*i+i1],v_t2);
/*             f[4*i+i0] = t1;   */
/*             f[1+4*i+i0] = t2; */
/*             f[2+4*i+i0] = t3; */
               _mm512_store_ps((float *)&f[4*i+i0],v_t1);
            }
/* loop over remaining elements */
            for (i = nxhs; i < nxh; i++) {
               t1 = f[4*i+i1];
               t2 = f[1+4*i+i1];
               t3 = f[2+4*i+i1];
               f[4*i+i1] = f[4*i+i0];
               f[1+4*i+i1] = f[1+4*i+i0];
               f[2+4*i+i1] = f[2+4*i+i0];
               f[4*i+i0] = t1;
               f[1+4*i+i0] = t2;
               f[2+4*i+i0] = t3;
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
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                      crealf(t1));
               i0 = ioff + j1;
               i1 = ioff + j2;
/* vector loop over elements in blocks of 2 */
               for (i = 0; i < nxhs; i+=2) {
/*                t2 = t1*f[4*i+i1];   */
/*                t3 = t1*f[1+4*i+i1]; */
/*                t4 = t1*f[2+4*i+i1]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*i+i1]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*i+i1] = f[4*i+i0] - t2;     */
/*                f[1+4*i+i1] = f[1+4*i+i0] - t3; */
/*                f[2+4*i+i1] = f[2+4*i+i0] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*i+i0]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+i1],v_t4);
/*                f[4*i+i0] += t2;   */
/*                f[1+4*i+i0] += t3; */
/*                f[2+4*i+i0] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+i0],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[4*i+i1];
                  t3 = t1*f[1+4*i+i1];
                  t4 = t1*f[2+4*i+i1];
                  f[4*i+i1] = f[4*i+i0] - t2;
                  f[1+4*i+i1] = f[1+4*i+i0] - t3;
                  f[2+4*i+i1] = f[2+4*i+i0] - t4;
                  f[4*i+i0] += t2;
                  f[1+4*i+i0] += t3;
                  f[2+4*i+i0] += t4;
               }
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rmx(float complex f[], int isign, int mixup[],
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
      ckncfft3rmxy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                   nxhyzd,nxyzhd);
/* perform z fft */
      ckncfft3rmz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                  nxhyzd,nxyzhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform z fft */
      ckncfft3rmz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                  nxhyzd,nxyzhd);
/* perform xy fft */
      ckncfft3rmxy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                   nxhyzd,nxyzhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rm3(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
/* wrapper function for 3 2d real to complex ffts, with packed data */
/* local data */
   int ny, nz;
   static int nyi = 1, nzi = 1;
/* calculate range of indices */
   ny = 1L<<indy;
   nz = 1L<<indz;
/* inverse fourier transform */
   if (isign < 0) {
/* perform xy fft */
      ckncfft3rm3xy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
                    nzd,nxhyzd,nxyzhd);
/* perform z fft */
      ckncfft3rm3z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                   nxhyzd,nxyzhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform z fft */
      ckncfft3rm3z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                  nxhyzd,nxyzhd);
/* perform xy fft */
      ckncfft3rm3xy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
                    nzd,nxhyzd,nxyzhd);
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void ckncgbppush3lt_(float *ppart, float *fxyz, float *bxyz ,int *kpic,
                     float *qbm, float *dt, float *dtc, float *ek,
                     int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                     int *mx, int *my, int *mz, int *nxv, int *nyv,
                     int *nzv, int *mx1, int *my1, int *mxyz1,
                     int *ipbc) {
   ckncgbppush3lt(ppart,fxyz,bxyz,kpic,*qbm,*dt,*dtc,ek,*idimp,*nppmx,
                  *nx,*ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,
                  *mxyz1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgbppushf3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                      int *ncl, int *ihole, float *qbm, float *dt,
                      float *dtc, float *ek, int *idimp, int *nppmx,
                      int *nx, int *ny, int *nz, int *mx, int *my,
                      int *mz, int *nxv, int *nyv, int *nzv, int *mx1,
                      int *my1, int *mxyz1, int *ntmax, int *irc) {
   ckncgbppushf3lt(ppart,fxyz,bxyz,kpic,ncl,ihole,*qbm,*dt,*dtc,ek,
                   *idimp,*nppmx,*nx,*ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,
                   *mx1,*my1,*mxyz1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgrbppush3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                      float *qbm, float *dt, float *dtc, float *ci,
                      float *ek, int *idimp, int *nppmx, int *nx,
                      int *ny, int *nz, int *mx, int *my, int *mz,
                      int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                      int *mxyz1, int *ipbc) {
   ckncgrbppush3lt(ppart,fxyz,bxyz,kpic,*qbm,*dt,*dtc,*ci,ek,*idimp,
                   *nppmx,*nx,*ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,
                   *my1,*mxyz1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgrbppushf3lt_(float *ppart, float *fxyz, float *bxyz,
                       int *kpic, int *ncl, int *ihole, float *qbm,
                       float *dt, float *dtc, float *ci, float *ek,
                       int *idimp, int *nppmx, int *nx, int *ny,
                       int *nz, int *mx, int *my, int *mz, int *nxv, 
                       int *nyv, int *nzv, int *mx1, int *my1,
                       int *mxyz1, int *ntmax, int *irc) {
   ckncgrbppushf3lt(ppart,fxyz,bxyz,kpic,ncl,ihole,*qbm,*dt,*dtc,*ci,ek,
                    *idimp,*nppmx,*nx,*ny,*nz,*mx,*my,*mz,*nxv,*nyv,
                    *nzv,*mx1,*my1,*mxyz1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgppost3lt_(float *ppart, float *q, int *kpic, float *qm,
                    int *nppmx, int *idimp, int *mx, int *my, int *mz,
                    int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                    int *mxyz1) {
   ckncgppost3lt(ppart,q,kpic,*qm,*nppmx,*idimp,*mx,*my,*mz,*nxv,*nyv,
                 *nzv,*mx1,*my1,*mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cknc2gppost3lt_(float *ppart, float *q, int *kpic, float *qm,
                     int *nppmx, int *idimp, int *mx, int *my, int *mz,
                     int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                     int *mxyz1) {
   cknc2gppost3lt(ppart,q,kpic,*qm,*nppmx,*idimp,*mx,*my,*mz,*nxv,*nyv,
                  *nzv,*mx1,*my1,*mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgjppost3lt_(float *ppart, float *cu, int *kpic, float *qm,
                     float *dt, int *nppmx, int *idimp, int *nx,
                     int *ny, int *nz, int *mx, int *my, int *mz,
                     int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                     int *mxyz1, int *ipbc) {
   ckncgjppost3lt(ppart,cu,kpic,*qm,*dt,*nppmx,*idimp,*nx,*ny,*nz,*mx,
                  *my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgjppostf3lt_(float *ppart, float *cu, int *kpic, int *ncl,
                      int *ihole, float *qm, float *dt, int *nppmx,
                      int *idimp, int *nx, int *ny, int *nz, int *mx,
                      int *my, int *mz, int *nxv, int *nyv, int *nzv,
                      int *mx1, int *my1, int *mxyz1, int *ntmax,
                      int *irc) {
   ckncgjppostf3lt(ppart,cu,kpic,ncl,ihole,*qm,*dt,*nppmx,*idimp,*nx,
                   *ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1,
                   *ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgrjppost3lt_(float *ppart, float *cu, int *kpic, float *qm,
                      float *dt, float *ci, int *nppmx, int *idimp,
                      int *nx, int *ny, int *nz, int *mx, int *my,
                      int *mz, int *nxv, int *nyv, int *nzv, int *mx1, 
                      int *my1, int *mxyz1, int *ipbc) {
   ckncgrjppost3lt(ppart,cu,kpic,*qm,*dt,*ci,*nppmx,*idimp,*nx,*ny,*nz,
                   *mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgrjppostf3lt_(float *ppart, float *cu, int *kpic, int *ncl,
                       int *ihole, float *qm, float *dt, float *ci,
                       int *nppmx, int *idimp, int *nx, int *ny,
                       int *nz, int *mx, int *my, int *mz, int *nxv,
                       int *nyv, int *nzv, int *mx1, int *my1,
                       int *mxyz1, int *ntmax, int *irc) {
   ckncgrjppostf3lt(ppart,cu,kpic,ncl,ihole,*qm,*dt,*ci,*nppmx,*idimp,
                    *nx,*ny,*nz,*mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,
                    *mxyz1,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cknc2gjppost3lt_(float *ppart, float *cu, int *kpic, float *qm,
                      float *dt, int *nppmx, int *idimp, int *nx,
                      int *ny, int *nz, int *mx, int *my, int *mz,
                      int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                      int *mxyz1, int *ipbc) {
   cknc2gjppost3lt(ppart,cu,kpic,*qm,*dt,*nppmx,*idimp,*nx,*ny,*nz,*mx,
                   *my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cknc2grjppost3lt_(float *ppart, float *cu, int *kpic, float *qm,
                       float *dt, float *ci, int *nppmx, int *idimp,
                       int *nx, int *ny, int *nz, int *mx, int *my,
                       int *mz, int *nxv, int *nyv, int *nzv, int *mx1,
                       int *my1, int *mxyz1, int *ipbc) {
   cknc2grjppost3lt(ppart,cu,kpic,*qm,*dt,*ci,*nppmx,*idimp,*nx,*ny,*nz,
                    *mx,*my,*mz,*nxv,*nyv,*nzv,*mx1,*my1,*mxyz1,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncpporder3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                     int *ihole, int *idimp, int *nppmx, int *nx,
                     int *ny, int *nz, int *mx, int *my, int *mz,
                     int *mx1, int *my1, int *mz1, int *npbmx,
                     int *ntmax, int *irc) {
   ckncpporder3lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*nx,*ny,*nz,
                  *mx,*my,*mz,*mx1,*my1,*mz1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncpporderf3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                      int *ihole, int *idimp, int *nppmx, int *mx1,
                      int *my1, int *mz1, int *npbmx, int *ntmax,
                      int *irc) {
   ckncpporderf3lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*mx1,*my1,
                   *mz1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncpp2order3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                      int *ihole, int *idimp, int *nppmx, int *nx,
                      int *ny, int *nz, int *mx, int *my, int *mz,
                      int *mx1, int *my1, int *mz1, int *npbmx,
                      int *ntmax, int *irc) {
   ckncpp2order3lt(ppart,ppbuff,kpic,ncl,ihole,*idimp,*nppmx,*nx,*ny,
                   *nz,*mx,*my,*mz,*mx1,*my1,*mz1,*npbmx,*ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cknccguard3l_(float *fxyz, int *nx, int *ny, int *nz, int *nxe,
                   int *nye, int *nze) {
   cknccguard3l(fxyz,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void ckncacguard3l_(float *cu, int *nx, int *ny, int *nz, int *nxe,
                    int *nye, int *nze) {
   ckncacguard3l(cu,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void ckncaguard3l_(float *q, int *nx, int *ny, int *nz, int *nxe,
                   int *nye, int *nze) {
   ckncaguard3l(q,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void ckncmpois33_(float complex *q, float complex *fxyz, int *isign,
                  float complex *ffc, float *ax, float *ay, float *az,
                  float *affp, float *we, int *nx, int *ny, int *nz,
                  int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
                  int *nzhd) {
   ckncmpois33(q,fxyz,*isign,ffc,*ax,*ay,*az,*affp,we,*nx,*ny,*nz,*nxvh,
               *nyv,*nzv,*nxhd,*nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ckncmcuperp3_(float complex *cu, int *nx, int *ny, int *nz,
                   int *nxvh, int *nyv, int *nzv) {
   ckncmcuperp3(cu,*nx,*ny,*nz,*nxvh,*nyv,*nzv);
   return;
}

/*--------------------------------------------------------------------*/
void ckncmibpois33_(float complex *cu, float complex *bxyz,
                    float complex *ffc, float *ci, float *wm, int *nx,
                    int *ny, int *nz, int *nxvh, int *nyv, int *nzv,
                    int *nxhd, int *nyhd, int *nzhd) {
   ckncmibpois33(cu,bxyz,ffc,*ci,wm,*nx,*ny,*nz,*nxvh,*nyv,*nzv,*nxhd,
                 *nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ckncmmaxwel3_(float complex *exyz, float complex *bxyz,
                   float complex *cu, float complex *ffc, float *ci,
                   float *dt, float *wf, float *wm, int *nx, int *ny,
                   int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd,
                   int *nyhd, int *nzhd) {
   ckncmmaxwel3(exyz,bxyz,cu,ffc,*ci,*dt,wf,wm,*nx,*ny,*nz,*nxvh,*nyv,
                *nzv,*nxhd,*nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ckncmemfield3_(float complex *fxyz, float complex *exyz,
                    float complex *ffc, int *isign, int *nx, int *ny,
                    int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd,
                    int *nyhd, int *nzhd) {
   ckncmemfield3(fxyz,exyz,ffc,*isign,*nx,*ny,*nz,*nxvh,*nyv,*nzv,*nxhd,
                 *nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rmx_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *indz,
                   int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                   int *nxyzhd) {
   ckncwfft3rmx(f,*isign,mixup,sct,*indx,*indy,*indz,*nxhd,*nyd,*nzd,
                *nxhyzd,*nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rm3_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *indz,
                   int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                   int *nxyzhd) {
   ckncwfft3rm3(f,*isign,mixup,sct,*indx,*indy,*indz,*nxhd,*nyd,*nzd,
                *nxhyzd,*nxyzhd);
   return;
}
