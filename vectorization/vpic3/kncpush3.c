/* KNC C Library for Skeleton 3D Electrostatic Vector PIC Code */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <immintrin.h>
#include "kncpush3.h"

/*--------------------------------------------------------------------*/
void ckncxiscan2(int *isdata, int nths) {
/* performs local prefix reduction of integer data shared by threads */
/* using binary tree method, exclusive scan. */
/* requires KNC, isdata needs to be 64 byte aligned */
/* local data */
   int j, ns, isum, ist;
   __m512i v_m1, v_m2, v_it0, v_it, v_is, v_ioff;
   __attribute__((aligned(64))) int ll[16];
   ns = 16*(nths/16);
   v_m1 = _mm512_set_epi32(11,11,11,11,11,10,9,8,3,3,3,3,3,2,1,0);
   v_m2 = _mm512_set_epi32(7,7,7,7,7,7,7,7,7,6,5,4,3,2,1,0);
   isum = 0;
   v_ioff = _mm512_setzero_epi32();
/* vector loop over elements in blocks of 16 */
   for (j = 0; j < ns; j+=16) {
/* load data */
      v_it0 = _mm512_load_epi32(&isdata[j]);
/* first pass */
      v_is = _mm512_shuffle_epi32(v_it0,177);
      v_it = _mm512_mask_add_epi32(v_it0,_mm512_int2mask(43690),v_it0,
             v_is);
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
      v_ioff = _mm512_shuffle_epi32(v_it,255);
      v_ioff = _mm512_permute4f128_epi32(v_ioff,255);
/* subtract for exclusive scan */
      v_it = _mm512_sub_epi32(v_it,v_it0);
/* write data */
      _mm512_store_epi32(&isdata[j],v_it);
   }
   if (ns > 0) {
      _mm512_store_epi32(ll,v_ioff);
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
void ckncgpush3lt(float part[], float fxyz[], float qbm, float dt,
                  float *ek, int idimp, int nop, int npe, int nx,
                  int ny, int nz, int nxv, int nyv, int nzv, int ipbc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space
   scalar version using guard cells
   94 flops/particle, 30 loads, 6 stores
   input: all, output: part, ek
   equations used are:
   vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
   vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
   vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
   where q/m is charge/mass, and
   x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
   z(t+dt) = z(t) + vz(t+dt/2)*dt
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
                  + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
                  + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
   fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
                  + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
                  + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = position z of particle n
   part[3][n] = velocity vx of particle n
   part[4][n] = velocity vy of particle n
   part[5][n] = velocity vz of particle n
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   qbm = particle charge/mass ratio
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
   (vz(t+dt/2)+vz(t-dt/2))**2)
   idimp = size of phase space = 6
   npe = first dimension of particle array
   nx/ny/nz = system length in x/y/z direction
   nxv = second dimension of field array, must be >= nx+1
   nyv = third dimension of field array, must be >= ny+1
   nzv = fourth dimension of field array, must be >= nz+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
   requires KNC, part needs to be 64 byte aligned
   npe needs to be a multiple of 16
   fxyz needs to have 4 components, although one is not used
local data                                                            */
   int j, nps, nn, mm, ll, nxyv;
   float qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx1, x, y, z, dx, dy, dz;
   float vx, vy, vz;
   double sum1, d0;
   __m512i v_nxv4, v_nxyv4;
   __m512i v_nn, v_mm, v_ll, v_it, v_perm;
   __m512 v_qtm, v_dt, v_one, v_zero;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_at, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, p, q, r, s;
   __m512d v_sum1, v_d;
   __mmask16 msk;
   __attribute__((aligned(64))) unsigned int kk[16];
   nxyv = nxv*nyv;
   nps = 16*(nop/16);
   qtm = qbm*dt;
   sum1 = 0.0;
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
   v_nxv4 = _mm512_set1_epi32(4*nxv);
   v_nxyv4 = _mm512_set1_epi32(4*nxyv);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qtm = _mm512_set1_ps(qtm);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_dt = _mm512_set1_ps(dt);
   v_edgelx = _mm512_set1_ps(edgelx);
   v_edgely = _mm512_set1_ps(edgely);
   v_edgelz = _mm512_set1_ps(edgelz);
   v_edgerx = _mm512_set1_ps(edgerx);
   v_edgery = _mm512_set1_ps(edgery);
   v_edgerz = _mm512_set1_ps(edgerz);
   v_sum1 = _mm512_set1_pd(0.0);
/* vector loop over particles in blocks of 16 */
   for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*    x = part[j];       */
/*    y = part[j+npe];   */
/*    z = part[j+2*npe]; */
      v_x = _mm512_load_ps(&part[j]);
      v_y = _mm512_load_ps(&part[j+npe]);
      v_z = _mm512_load_ps(&part[j+2*npe]);
/*    nn = x; */
/*    mm = y; */
/*    ll = z; */
      v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    dxp = x - (float) nn; */
/*    dyp = y - (float) mm; */
/*    dzp = z - (float) ll; */
      v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dxp = _mm512_sub_ps(v_x,v_dxp);
      v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dyp = _mm512_sub_ps(v_y,v_dyp);
      v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*    nn = 4*(nn + nxv*mm + nxyv*ll); */
      v_it = _mm512_mullo_epi32(v_nxyv4,v_ll);
      v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_nxv4,v_mm));
      v_nn = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*    amx = 1.0f - dxp; */
/*    amy = 1.0f - dyp; */
/*    amz = 1.0f - dzp; */
      v_amx = _mm512_sub_ps(v_one,v_dxp);
      v_amy = _mm512_sub_ps(v_one,v_dyp);
      v_amz = _mm512_sub_ps(v_one,v_dzp);
/*    dx1 = dxp*dyp; */
/*    dyp = amx*dyp; */
/*    amx = amx*amy; */
/*    amy = dxp*amy; */
      v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
      v_dyp = _mm512_mul_ps(v_amx,v_dyp);
      v_amx = _mm512_mul_ps(v_amx,v_amy);
      v_amy = _mm512_mul_ps(v_dxp,v_amy);
      _mm512_store_epi32(kk,v_nn);
/* load fxyz[nn:nn+3] and fxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
      mm = kk[0];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[nn:nn+3] field components */
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
/* perform 16x3 transpose for fxyz[nn+4:nn+7] field components */
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
/* find first part of acceleration */
/*    dx = amx*fxyz[nn] + amy*fxyz[nn+4];                     */
      v_dx = _mm512_mul_ps(v_amx,a);
      v_dx = _mm512_fmadd_ps(v_amy,p,v_dx);
/*    dy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];                 */
      v_dy = _mm512_mul_ps(v_amx,b);
      v_dy = _mm512_fmadd_ps(v_amy,q,v_dy);
/*    dz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];                 */
      v_dz = _mm512_mul_ps(v_amx,c);
      v_dz = _mm512_fmadd_ps(v_amy,r,v_dz);
/*    mm = nn + 4*nxv;                                     */
/* load fxyz[mm:mm+3] and fxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
      mm = kk[0] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[mm:mm+3] field components */
/* where mm = nn + 4*nxv;                                    */
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
/* perform 16x3 transpose for fxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*nxv;                                      */
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
/* find second part of acceleration */
/*    dx = amz*(dx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);          */
      v_dx = _mm512_fmadd_ps(v_dyp,a,v_dx);
      v_dx = _mm512_fmadd_ps(v_dx1,p,v_dx);
      v_dx = _mm512_mul_ps(v_amz,v_dx);
/*    dy = amz*(dy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]);      */
      v_dy = _mm512_fmadd_ps(v_dyp,b,v_dy);
      v_dy = _mm512_fmadd_ps(v_dx1,q,v_dy);
      v_dy = _mm512_mul_ps(v_amz,v_dy);
/*    dz = amz*(dz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]);      */
      v_dz = _mm512_fmadd_ps(v_dyp,c,v_dz);
      v_dz = _mm512_fmadd_ps(v_dx1,r,v_dz);
      v_dz = _mm512_mul_ps(v_amz,v_dz);
/*    nn += 4*nxyv;                                           */
      v_nn = _mm512_add_epi32(v_nn,v_nxyv4);
      _mm512_store_epi32(kk,v_nn);
/* load fxyz[nn:nn+3] and fxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
      mm = kk[0];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[nn:nn+3] field components */
/* where nn = nn + 4*nxyv;                                   */
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
/* perform 16x3 transpose for fxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*nxyv;                                     */
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
/* find third part of acceleration */
/*    vx = amx*fxyz[nn] + amy*fxyz[nn+4];                     */
      v_vx = _mm512_mul_ps(v_amx,a);
      v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*    vy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];                 */
      v_vy = _mm512_mul_ps(v_amx,b);
      v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*    vz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];                 */
      v_vz = _mm512_mul_ps(v_amx,c);
      v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*    mm = nn + 4*nxv;                                        */
/* load fxyz[mm:mm+3] and fxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
      mm = kk[0] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[mm:mm+3] field components */
/* where mm = nn + 4*nxyv + 4*nxv;                           */
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
/* perform 16x3 transpose for fxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*nxyv + 4*nxv;                           */
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
/* find fourth part of acceleration */
/*    dx = dx + dzp*(vx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);     */
      v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
      v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
      v_dx = _mm512_fmadd_ps(v_dzp,v_vx,v_dx);
/*    dy = dy + dzp*(vy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]); */
      v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
      v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
      v_dy = _mm512_fmadd_ps(v_dzp,v_vy,v_dy);
/*    dz = dz + dzp*(vz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]); */
      v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
      v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
      v_dz = _mm512_fmadd_ps(v_dzp,v_vz,v_dz);
/* new velocity */
/*    dxp = part[j+3*npe]; */
/*    dyp = part[j+4*npe]; */
/*    dzp = part[j+5*npe]; */
      v_dxp = _mm512_load_ps(&part[j+3*npe]);
      v_dyp = _mm512_load_ps(&part[j+4*npe]);
      v_dzp = _mm512_load_ps(&part[j+5*npe]);
/*    vx = dxp + qtm*dx; */
/*    vy = dyp + qtm*dy; */
/*    vz = dzp + qtm*dz; */
      v_vx = _mm512_fmadd_ps(v_qtm,v_dx,v_dxp);
      v_vy = _mm512_fmadd_ps(v_qtm,v_dy,v_dyp);
      v_vz = _mm512_fmadd_ps(v_qtm,v_dz,v_dzp);
/* average kinetic energy */
/*    dxp += vx; */
/*    dyp += vy; */
/*    dzp += vz; */
      v_dxp = _mm512_add_ps(v_dxp,v_vx);
      v_dyp = _mm512_add_ps(v_dyp,v_vy);
      v_dzp = _mm512_add_ps(v_dzp,v_vz);
/*    sum1 += dxp*dxp + dyp*dyp + dzp*dzp; */
      v_at = _mm512_mul_ps(v_dxp,v_dxp);
      v_at = _mm512_add_ps(v_at,_mm512_mul_ps(v_dyp,v_dyp));
      v_at = _mm512_add_ps(v_at,_mm512_mul_ps(v_dzp,v_dzp));
/* convert to double precision before accumulating */
      v_sum1 = _mm512_add_pd(v_sum1,_mm512_cvtpslo_pd(v_at));
      v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_at,78));
      v_sum1 = _mm512_add_pd(v_sum1,v_d);
/* new position */
/*    dx = x + vx*dt; */
/*    dy = y + vy*dt; */
/*    dz = z + vz*dt; */
      v_dx = _mm512_fmadd_ps(v_vx,v_dt,v_x);
      v_dy = _mm512_fmadd_ps(v_vy,v_dt,v_y);
      v_dz = _mm512_fmadd_ps(v_vz,v_dt,v_z);
/* periodic boundary conditions */
      if (ipbc==1) {
/*       if (dx < edgelx) dx += edgerx; */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         v_dx = _mm512_mask_add_ps(v_dx,msk,v_dx,v_edgerx);
/*       if (dx >= edgerx) dx -= edgerx; */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgerx,_MM_CMPINT_GE);
         v_dx = _mm512_mask_sub_ps(v_dx,msk,v_dx,v_edgerx);
/*       if (dy < edgely) dy += edgery; */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         v_dy = _mm512_mask_add_ps(v_dy,msk,v_dy,v_edgery);
/*       if (dy >= edgery) dy -= edgery; */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgery,_MM_CMPINT_GE);
         v_dy = _mm512_mask_sub_ps(v_dy,msk,v_dy,v_edgery);
/*       if (dz < edgelz) dz += edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         v_dz = _mm512_mask_add_ps(v_dz,msk,v_dz,v_edgerz);
/*       if (dz >= edgerz) dz -= edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         v_dz = _mm512_mask_sub_ps(v_dz,msk,v_dz,v_edgerz);
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          vx = -vx;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
               _MM_CMPINT_GE));
         v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
         v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          vy = -vy;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
               _MM_CMPINT_GE));
         v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
         v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/*       if ((dz < edgelz) || (dz >= edgerz)) { */
/*          dz = z;                             */
/*          vz = -vz;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dz,v_edgerz,
               _MM_CMPINT_GE));
         v_dz = _mm512_mask_blend_ps(msk,v_dz,v_z);
         v_vz = _mm512_mask_sub_ps(v_vz,msk,v_zero,v_vz);
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          vx = -vx;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
               _MM_CMPINT_GE));
         v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
         v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          vy = -vy;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
               _MM_CMPINT_GE));
         v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
         v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/*       if (dz < edgelz) dz += edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         v_dz = _mm512_mask_add_ps(v_dz,msk,v_dz,v_edgerz);
/*       if (dz >= edgerz) dz -= edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         v_dz = _mm512_mask_sub_ps(v_dz,msk,v_dz,v_edgerz);
      }
/* set new position */
/*    part[j] = dx;       */
/*    part[j+npe] = dy;   */
/*    part[j+2*npe] = dz; */
      _mm512_store_ps(&part[j],v_dx);
      _mm512_store_ps(&part[j+npe],v_dy);
      _mm512_store_ps(&part[j+2*npe],v_dz);
/* set new velocity */
/*    part[j+3*npe] = vx; */
/*    part[j+4*npe] = vy; */
/*    part[j+5*npe] = vz; */
      _mm512_store_ps(&part[j+3*npe],v_vx);
      _mm512_store_ps(&part[j+4*npe],v_vy);
      _mm512_store_ps(&part[j+5*npe],v_vz);
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      z = part[j+2*npe];
      nn = x;
      mm = y;
      ll = z;
      dxp = x - (float) nn;
      dyp = y - (float) mm;
      dzp = z - (float) ll;
      nn = 4*(nn + nxv*mm + nxyv*ll);
      amx = 1.0f - dxp;
      amy = 1.0f - dyp;
      dx1 = dxp*dyp;
      dyp = amx*dyp;
      amx = amx*amy;
      amz = 1.0f - dzp;
      amy = dxp*amy;
/* find acceleration */
      dx = amx*fxyz[nn] + amy*fxyz[nn+4];
      dy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];
      dz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];
      mm = nn + 4*nxv;
      dx = amz*(dx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);
      dy = amz*(dy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]);
      dz = amz*(dz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]);
      nn += 4*nxyv;
      vx = amx*fxyz[nn] + amy*fxyz[nn+4];
      vy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];
      vz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];
      mm = nn + 4*nxv;
      dx = dx + dzp*(vx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);
      dy = dy + dzp*(vy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]);
      dz = dz + dzp*(vz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]);
/* new velocity */
      dxp = part[j+3*npe];
      dyp = part[j+4*npe];
      dzp = part[j+5*npe];
      vx = dxp + qtm*dx;
      vy = dyp + qtm*dy;
      vz = dzp + qtm*dz;
/* average kinetic energy */
      dxp += vx;
      dyp += vy;
      dzp += vz;
      sum1 += dxp*dxp + dyp*dyp + dzp*dzp;
/* new position */
      dx = x + vx*dt;
      dy = y + vy*dt;
      dz = z + vz*dt;
/* periodic boundary conditions */
      if (ipbc==1) {
         if (dx < edgelx) dx += edgerx;
         if (dx >= edgerx) dx -= edgerx;
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
         if (dz < edgelz) dz += edgerz;
         if (dz >= edgerz) dz -= edgerz;
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
         if (dz < edgelz) dz += edgerz;
         if (dz >= edgerz) dz -= edgerz;
      }
/* set new position */
      part[j] = dx;
      part[j+npe] = dy;
      part[j+2*npe] = dz;
/* set new velocity */
      part[j+3*npe] = vx;
      part[j+4*npe] = vy;
      part[j+5*npe] = vz;
   }
/* normalize kinetic energy */
/* *ek += 0.125f*sum1; */
   d0 = _mm512_reduce_add_pd(v_sum1);
   *ek += 0.125f*(sum1 + d0);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgpost3lt(float part[], float q[], float qm, int nop, int npe,
                  int idimp, int nxv, int nyv, int nzv) {
/* for 3d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   scalar version using guard cells
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
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = position z of particle n
   q[l][k][j] = charge density at grid point j,k,l
   qm = charge on particle, in units of e
   npe = first dimension of particle array
   idimp = size of phase space = 6
   nxv = first dimension of charge array, must be >= nx+1
   nyv = second dimension of charge array, must be >= ny+1
   nzv = third dimension of charge array, must be >= nz+1
   requires KNC, part needs to be 64 byte aligned
   npe needs to be a multiple of 16
local data                                                            */
   int i, j, nps, nn, mm, ll, nxyv;
   float x, y, z, w, dx1, dxp, dyp, dzp, amx, amy, amz;
   __m512i v_nxv, v_nxyv;
   __m512i v_nn, v_mm, v_ll, v_it;
   __m512 v_qm, v_one;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_as, v_at;
   __m512 a, b, c, d, e, f, g, h, qp, qr;
   __mmask16 msk, msks;
   __attribute__((aligned(64))) unsigned int kk[16];
   nxyv = nxv*nyv;
   nps = 16*(nop/16);
   v_nxv = _mm512_set1_epi32(nxv);
   v_nxyv = _mm512_set1_epi32(nxyv);
   v_qm = _mm512_set1_ps(qm);
   v_one = _mm512_set1_ps(1.0f);
/* vector loop over particles in blocks of 16 */
   for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*    x = part[j];       */
/*    y = part[j+npe];   */
/*    z = part[j+2*npe]; */
      v_x = _mm512_load_ps(&part[j]);
      v_y = _mm512_load_ps(&part[j+npe]);
      v_z = _mm512_load_ps(&part[j+2*npe]);
/*    nn = x; */
/*    mm = y; */
/*    ll = z; */
      v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    dxp = qm*(x - (float) nn); */
/*    dyp = y - (float) mm; */
/*    dzp = z - (float) ll; */
      v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
      v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dyp = _mm512_sub_ps(v_y,v_dyp);
      v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*    nn = nn + nxv*mm + nxyv*ll; */
      v_it = _mm512_mullo_epi32(v_nxyv,v_ll);
      v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_nxv,v_mm));
      v_nn = _mm512_add_epi32(v_nn,v_it);
/*    amx = qm - dxp;   */
/*    amy = 1.0f - dyp; */
/*    amz = 1.0f - dzp; */
      v_amx = _mm512_sub_ps(v_qm,v_dxp);
      v_amy = _mm512_sub_ps(v_one,v_dyp);
      v_amz = _mm512_sub_ps(v_one,v_dzp);
/*    dx1 = dxp*dyp; */
/*    dyp = amx*dyp; */
/*    amx = amx*amy; */
/*    amy = dxp*amy; */
      v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
      v_dyp = _mm512_mul_ps(v_amx,v_dyp);
      v_amx = _mm512_mul_ps(v_amx,v_amy);
      v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*    a = amx*amz; */
/*    b = amy*amz; */
/*    d = dyp*amz; */
/*    d = dx1*amz; */
      a = _mm512_mul_ps(v_amx,v_amz);
      b = _mm512_mul_ps(v_amy,v_amz);
      c = _mm512_mul_ps(v_dyp,v_amz);
      d = _mm512_mul_ps(v_dx1,v_amz);
/*    e = amx*dzp; */
/*    f = amy*dzp; */
/*    g = dyp*dzp; */
/*    h = dx1*dzp; */
      e = _mm512_mul_ps(v_amx,v_dzp);
      f = _mm512_mul_ps(v_amy,v_dzp);
      g = _mm512_mul_ps(v_dyp,v_dzp);
      h = _mm512_mul_ps(v_dx1,v_dzp);
      _mm512_store_epi32(kk,v_nn);
/* deposit charge */
/*    x = q[nn] + amx*amz;       */
/*    y = q[nn+1] + amy*amz;     */
/*    z = q[nn+nxv] + dyp*amz;   */
/*    w = q[nn+1+nxv] + dx1*amz; */
/*    q[nn] = x;                 */
/*    q[nn+1] = y;               */
/*    q[nn+nxv] = z;             */
/*    q[nn+1+nxv] = w;           */
/*    mm = nn + nxyv;            */
/*    x = q[mm] + amx*dzp;       */
/*    y = q[mm+1] + amy*dzp;     */
/*    z = q[mm+nxv] + dyp*dzp;   */
/*    w = q[mm+1+nxv] + dx1*dzp; */
/*    q[mm] = x;                 */
/*    q[mm+1] = y;               */
/*    q[mm+nxv] = z;             */
/*    q[mm+1+nxv] = w;           */
/* deposit charge for two particles at a time */
      for (i = 0; i < 8; i++) {
/* first particle */
         mm = kk[2*i];
         msk = _mm512_int2mask(3<<(2*i));
         msks = _mm512_int2mask(2<<(2*i));
         qp = _mm512_mask_loadunpacklo_ps(qp,msk,&q[mm]);
         qp = _mm512_mask_loadunpackhi_ps(qp,msk,&q[mm+16]);
         v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)a,msks,
                (__m512i)b,177);
         qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
         _mm512_mask_packstorelo_ps(&q[mm],msk,qp);
         _mm512_mask_packstorehi_ps(&q[mm+16],msk,qp);
         ll = mm + nxv;
         qr = _mm512_mask_loadunpacklo_ps(qr,msk,&q[ll]);
         qr = _mm512_mask_loadunpackhi_ps(qr,msk,&q[ll+16]);
         v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)c,msks,
                (__m512i)d,177);
         qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
         _mm512_mask_packstorelo_ps(&q[ll],msk,qr);
         _mm512_mask_packstorehi_ps(&q[ll+16],msk,qr);
         mm = mm + nxyv;
         qp = _mm512_mask_loadunpacklo_ps(qp,msk,&q[mm]);
         qp = _mm512_mask_loadunpackhi_ps(qp,msk,&q[mm+16]);
         v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)e,msks,
                (__m512i)f,177);
         qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
         _mm512_mask_packstorelo_ps(&q[mm],msk,qp);
         _mm512_mask_packstorehi_ps(&q[mm+16],msk,qp);
         ll = mm + nxv;
         qr = _mm512_mask_loadunpacklo_ps(qr,msk,&q[ll]);
         qr = _mm512_mask_loadunpackhi_ps(qr,msk,&q[ll+16]);
         v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)g,msks,
                (__m512i)h,177);
         qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
         _mm512_mask_packstorelo_ps(&q[ll],msk,qr);
         _mm512_mask_packstorehi_ps(&q[ll+16],msk,qr);
/* second particle */
         mm = kk[2*i+1];
         msks = _mm512_int2mask(1<<(2*i));
         qp = _mm512_mask_loadunpacklo_ps(qp,msk,&q[mm]);
         qp = _mm512_mask_loadunpackhi_ps(qp,msk,&q[mm+16]);
         v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)b,msks,
                (__m512i)a,177);
         qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
         _mm512_mask_packstorelo_ps(&q[mm],msk,qp);
         _mm512_mask_packstorehi_ps(&q[mm+16],msk,qp);
         ll = mm + nxv;
         qr = _mm512_mask_loadunpacklo_ps(qr,msk,&q[ll]);
         qr = _mm512_mask_loadunpackhi_ps(qr,msk,&q[ll+16]);
         v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)d,msks,
                (__m512i)c,177);
         qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
         _mm512_mask_packstorelo_ps(&q[ll],msk,qr);
         _mm512_mask_packstorehi_ps(&q[ll+16],msk,qr);
         mm = mm + nxyv;
         qp = _mm512_mask_loadunpacklo_ps(qp,msk,&q[mm]);
         qp = _mm512_mask_loadunpackhi_ps(qp,msk,&q[mm+16]);
         v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)f,msks,
                (__m512i)e,177);
         qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
         _mm512_mask_packstorelo_ps(&q[mm],msk,qp);
         _mm512_mask_packstorehi_ps(&q[mm+16],msk,qp);
         ll = mm + nxv;
         qr = _mm512_mask_loadunpacklo_ps(qr,msk,&q[ll]);
         qr = _mm512_mask_loadunpackhi_ps(qr,msk,&q[ll+16]);
         v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)h,msks,
                (__m512i)g,177);
         qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
         _mm512_mask_packstorelo_ps(&q[ll],msk,qr);
         _mm512_mask_packstorehi_ps(&q[ll+16],msk,qr);
      }
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      z = part[j+2*npe];
      nn = x;
      mm = y;
      ll = z;
      dxp = qm*(x - (float) nn);
      dyp = y - (float) mm;
      dzp = z - (float) ll;
      nn = nn + nxv*mm + nxyv*ll;
      amx = qm - dxp;
      amy = 1.0f - dyp;
      amz = 1.0f - dzp;
      dx1 = dxp*dyp;
      dyp = amx*dyp;
      amx = amx*amy;
      amy = dxp*amy;
/* deposit charge */
      x = q[nn] + amx*amz;
      y = q[nn+1] + amy*amz;
      z = q[nn+nxv] + dyp*amz;
      w = q[nn+1+nxv] + dx1*amz;
      q[nn] = x;
      q[nn+1] = y;
      q[nn+nxv] = z;
      q[nn+1+nxv] = w;
      mm = nn + nxyv;
      x = q[mm] + amx*dzp;
      y = q[mm+1] + amy*dzp;
      z = q[mm+nxv] + dyp*dzp;
      w = q[mm+1+nxv] + dx1*dzp;
      q[mm] = x;
      q[mm+1] = y;
      q[mm+nxv] = z;
      q[mm+1+nxv] = w;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cknc2gpost3lt(float part[], float q[], float qm, int nop, int npe,
                   int idimp, int nxv, int nyv, int nzv) {
/* for 3d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   scalar version using guard cells
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
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = position z of particle n
   q[l][k][j] = charge density at grid point j,k,l
   qm = charge on particle, in units of e
   npe = first dimension of particle array
   idimp = size of phase space = 6
   nxv = first dimension of charge array, must be >= nx+1
   nyv = second dimension of charge array, must be >= ny+1
   nzv = third dimension of charge array, must be >= nz+1
   requires KNC, part needs to be 64 byte aligned
   npe needs to be a multiple of 16
local data                                                            */
   int i, j, nps, nn, mm, ll, nxyv;
   float x, y, z, w, dx1, dxp, dyp, dzp, amx, amy, amz;
   __m512i v_nxv, v_nxyv;
   __m512i v_nn, v_mm, v_ll, v_it;
   __m512 v_qm, v_one;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1;
   __attribute__((aligned(64))) unsigned int kk[16];
   typedef union vfloat {float v[16]; __m512 v16;} vf;
   vf vv[8];
   nxyv = nxv*nyv;
   nps = 16*(nop/16);
   v_nxv = _mm512_set1_epi32(nxv);
   v_nxyv = _mm512_set1_epi32(nxyv);
   v_qm = _mm512_set1_ps(qm);
   v_one = _mm512_set1_ps(1.0f);
/* vector loop over particles in blocks of 16 */
   for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*    x = part[j];       */
/*    y = part[j+npe];   */
/*    z = part[j+2*npe]; */
      v_x = _mm512_load_ps(&part[j]);
      v_y = _mm512_load_ps(&part[j+npe]);
      v_z = _mm512_load_ps(&part[j+2*npe]);
/*    nn = x; */
/*    mm = y; */
/*    ll = z; */
      v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    dxp = qm*(x - (float) nn); */
/*    dyp = y - (float) mm; */
/*    dzp = z - (float) ll; */
      v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
      v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dyp = _mm512_sub_ps(v_y,v_dyp);
      v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*    nn = nn + nxv*mm + nxyv*ll; */
      v_it = _mm512_mullo_epi32(v_nxyv,v_ll);
      v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_nxv,v_mm));
      v_nn = _mm512_add_epi32(v_nn,v_it);
/*    amx = qm - dxp;   */
/*    amy = 1.0f - dyp; */
/*    amz = 1.0f - dzp; */
      v_amx = _mm512_sub_ps(v_qm,v_dxp);
      v_amy = _mm512_sub_ps(v_one,v_dyp);
      v_amz = _mm512_sub_ps(v_one,v_dzp);
/*    dx1 = dxp*dyp; */
/*    dyp = amx*dyp; */
/*    amx = amx*amy; */
/*    amy = dxp*amy; */
      v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
      v_dyp = _mm512_mul_ps(v_amx,v_dyp);
      v_amx = _mm512_mul_ps(v_amx,v_amy);
      v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*    x = amx*amz; */
/*    y = amy*amz; */
/*    z = dyp*amz; */
/*    w = dx1*amz; */
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
/*    x = q[nn] + amx*amz;       */
/*    y = q[nn+1] + amy*amz;     */
/*    z = q[nn+nxv] + dyp*amz;   */
/*    w = q[nn+1+nxv] + dx1*amz; */
/*    q[nn] = x;                 */
/*    q[nn+1] = y;               */
/*    q[nn+nxv] = z;             */
/*    q[nn+1+nxv] = w;           */
/*    mm = nn + nxyv;            */
/*    x = q[mm] + amx*dzp;       */
/*    y = q[mm+1] + amy*dzp;     */
/*    z = q[mm+nxv] + dyp*dzp;   */
/*    w = q[mm+1+nxv] + dx1*dzp; */
/*    q[mm] = x;                 */
/*    q[mm+1] = y;               */
/*    q[mm+nxv] = z;             */
/*    q[mm+1+nxv] = w;           */
      for (i = 0; i < 16; i++) {
         nn = kk[i];
         x = q[nn] + vv[0].v[i];
         y = q[nn+1] + vv[1].v[i];
         z = q[nn+nxv] + vv[2].v[i];
         w = q[nn+1+nxv] + vv[3].v[i];
         q[nn] = x;
         q[nn+1] = y;
         q[nn+nxv] = z;
         q[nn+1+nxv] = w;
         mm = nn + nxyv;
         x = q[mm] + vv[4].v[i];
         y = q[mm+1] + vv[5].v[i];
         z = q[mm+nxv] + vv[6].v[i];
         w = q[mm+1+nxv] + vv[7].v[i];
         q[mm] = x;
         q[mm+1] = y;
         q[mm+nxv] = z;
         q[mm+1+nxv] = w;
      }
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      z = part[j+2*npe];
      nn = x;
      mm = y;
      ll = z;
      dxp = qm*(x - (float) nn);
      dyp = y - (float) mm;
      dzp = z - (float) ll;
      nn = nn + nxv*mm + nxyv*ll;
      amx = qm - dxp;
      amy = 1.0f - dyp;
      amz = 1.0f - dzp;
      dx1 = dxp*dyp;
      dyp = amx*dyp;
      amx = amx*amy;
      amy = dxp*amy;
/* deposit charge */
      x = q[nn] + amx*amz;
      y = q[nn+1] + amy*amz;
      z = q[nn+nxv] + dyp*amz;
      w = q[nn+1+nxv] + dx1*amz;
      q[nn] = x;
      q[nn+1] = y;
      q[nn+nxv] = z;
      q[nn+1+nxv] = w;
      mm = nn + nxyv;
      x = q[mm] + amx*dzp;
      y = q[mm+1] + amy*dzp;
      z = q[mm+nxv] + dyp*dzp;
      w = q[mm+1+nxv] + dx1*dzp;
      q[mm] = x;
      q[mm+1] = y;
      q[mm+nxv] = z;
      q[mm+1+nxv] = w;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncdsortp3yzlt(float parta[], float partb[], int npic[],
                     int idimp, int nop, int npe, int ny1, int nyz1) {
/* this subroutine sorts particles by y,z grid
   linear interpolation
   part = particle array
   part[1][n] = position y of particle n
   part[2][n] = position z of particle n
   npic = address offset for reordering particles
   idimp = size of phase space = 6
   npe = first dimension of particle array
   ny1 = system length in y direction + 1
   nyz1 = ny1*nz1, where nz1 = system length in z direction + 1
   requires KNC, part needs to be 64 byte aligned
   npe needs to be a multiple of 16
local data                                                            */
   int i, j, k, m, l, nps, ip;
/* int isum, ist; */
   __m512i v_ny1, v_m, v_l, v_it;
   __m512 v_y, v_z, v_at;
   __attribute__((aligned(64))) unsigned int pp[16];
   typedef union vint {int ll[16]; __m512i v_l;} vi;
   vi ii;
   nps = 16*(nop/16);
   v_ny1 = _mm512_set1_epi32(ny1);
/* clear counter array */
/* for (k = 0; k < nyz1; k++) { */
/*    npic[k] = 0;              */
/* }                            */
   memset((void *)npic,0,nyz1*sizeof(int));
/* find how many particles in each grid */
/* vector loop over particles in blocks of 4 */
   for (j = 0; j < nps; j+=16) {
/*    m = parta[j+npe];   */
/*    l = parta[j+2*npe]; */
      v_y = _mm512_load_ps(&parta[j+npe]);
      v_z = _mm512_load_ps(&parta[j+2*npe]);
      v_m = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
            _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_l = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
            _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    l = m + ny1*l; */
      ii.v_l =_mm512_fmadd_epi32(v_ny1,v_l,v_m);
      for (k = 0; k < 16; k++) {
         l = ii.ll[k];
         npic[l] += 1;
      }
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
      m = parta[j+npe];
      l = parta[j+2*npe];
      l = m + ny1*l;
      npic[l] += 1;
   }
/* find address offset */
/* isum = 0;                    */
/* for (k = 0; k < nyz1; k++) { */
/*    ist = npic[k];            */
/*    npic[k] = isum;           */
/*    isum += ist;              */
/* }                            */
   ckncxiscan2(npic,nyz1);
/* find addresses of particles at each grid and reorder particles */
/* vector loop over particles in blocks of 4 */
   for (j = 0; j < nps; j+=16) {
/*    m = parta[j+npe];   */
/*    l = parta[j+2*npe]; */
      v_y = _mm512_load_ps(&parta[j+npe]);
      v_z = _mm512_load_ps(&parta[j+2*npe]);
      v_m = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
            _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_l = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
            _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    l = m + ny1*l; */
      ii.v_l = _mm512_fmadd_epi32(v_ny1,v_l,v_m);
      for (k = 0; k < 16; k++) {
         l = ii.ll[k];
         ip = npic[l];
         npic[l] = ip + 1;
         pp[k] = ip;
      }
      v_it = _mm512_load_epi32(pp);
      for (i = 0; i < idimp; i++) {
/*          partb[ip+npe*i] = parta[j+npe*i]; */
            v_at = _mm512_load_ps(&parta[j+npe*i]);
            _mm512_i32scatter_ps(&partb[npe*i],v_it,v_at,4);
      }
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
      m = parta[j+npe];
      l = parta[j+2*npe];
      l = m + ny1*l;
      ip = npic[l];
      for (i = 0; i < idimp; i++) {
         partb[ip+npe*i] = parta[j+npe*i];
      }
      npic[l] = ip + 1;
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
#define N 4
   int j, k, l, nxs, nxyen, ll;
   nxs = 4*(nx/4);
   nxyen = N*nxe*nye;
/* copy edges of extended field */
   for (l = 0; l < nz; l++) {
      ll = nxyen*l;
      for (k = 0; k < ny; k++) {
         fxyz[N*nx+N*nxe*k+ll] = fxyz[N*nxe*k+ll];
         fxyz[1+N*nx+N*nxe*k+ll] = fxyz[1+N*nxe*k+ll];
         fxyz[2+N*nx+N*nxe*k+ll] = fxyz[2+N*nxe*k+ll];
      }
/* vector loop over elements in blocks of 4 */
      for (j = 0; j < nxs; j+=4) {
         _mm512_mask_store_ps(&fxyz[N*j+N*nxe*ny+ll],
            _mm512_int2mask(30583),_mm512_load_ps(&fxyz[N*j+ll]));
      }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
         fxyz[N*j+N*nxe*ny+ll] = fxyz[N*j+ll];
         fxyz[1+N*j+N*nxe*ny+ll] = fxyz[1+N*j+ll];
         fxyz[2+N*j+N*nxe*ny+ll] = fxyz[2+N*j+ll];
      }
      fxyz[N*nx+N*nxe*ny+ll] = fxyz[ll];
      fxyz[1+N*nx+N*nxe*ny+ll] = fxyz[1+ll];
      fxyz[2+N*nx+N*nxe*ny+ll] = fxyz[2+ll];
   }
   for (k = 0; k < ny; k++) {
/* vector loop over elements in blocks of 4 */
      for (j = 0; j < nxs; j+=4) {
         _mm512_mask_store_ps(&fxyz[N*j+N*nxe*k+nxyen*nz],
            _mm512_int2mask(30583),_mm512_load_ps(&fxyz[N*j+N*nxe*k]));
      }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
         fxyz[N*j+N*nxe*k+nxyen*nz] = fxyz[N*j+N*nxe*k];
         fxyz[1+N*j+N*nxe*k+nxyen*nz] = fxyz[1+N*j+N*nxe*k];
         fxyz[2+N*j+N*nxe*k+nxyen*nz] = fxyz[2+N*j+N*nxe*k];
      }
      fxyz[N*nx+N*nxe*k+nxyen*nz] = fxyz[N*nxe*k];
      fxyz[1+N*nx+N*nxe*k+nxyen*nz] = fxyz[1+N*nxe*k];
      fxyz[2+N*nx+N*nxe*k+nxyen*nz] = fxyz[2+N*nxe*k];
   }
/* vector loop over elements in blocks of 4 */
      for (j = 0; j < nxs; j+=4) {
         _mm512_mask_store_ps(&fxyz[N*j+N*nxe*ny+nxyen*nz],
            _mm512_int2mask(30583),_mm512_load_ps(&fxyz[N*j]));
   }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
      fxyz[N*j+N*nxe*ny+nxyen*nz] = fxyz[N*j];
      fxyz[1+N*j+N*nxe*ny+nxyen*nz] = fxyz[1+N*j];
      fxyz[2+N*j+N*nxe*ny+nxyen*nz] = fxyz[2+N*j];
   }
   fxyz[N*nx+N*nxe*ny+nxyen*nz] = fxyz[0];
   fxyz[1+N*nx+N*nxe*ny+nxyen*nz] = fxyz[1];
   fxyz[2+N*nx+N*nxe*ny+nxyen*nz] = fxyz[2];
   return;
#undef N
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
void ckncpois33(float complex q[], float complex fxyz[], int isign,
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
   requires KNC, q, fxy, ffc need to be 64 byte aligned
   nxhd, nxvh need to be a multiple of 8
   fxyz needs to have 4 components
local data                                                 */
   int nxh, nyh, nzh, nxhs, itn, j, k, l, k1, l1, kk, kj, ll, lj;
   int nxyhd, nxvyh; 
   float dnx, dny, dnz, dkx, dky, dkz, at1, at2, at3, at4, at5, at6;
   float complex zero, zt1, zt2;
   double wp, d0;
   __m512i v_j, v_it, v_perm;
   __m512 v_dnx, v_dny, v_dnz, v_dky, v_dkz, v_at1, v_at2, v_at3, v_at4;
   __m512 v_zero, v_zt1, v_zt2, v_zt3, v_zt4;
   __m512 a, b, c, d, e, f, g, h;
   __m512d v_wp, v_d;
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
L40: wp = 0.0;
   v_wp = _mm512_setzero_pd();
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
   for (l = 1; l < nzh; l++) {
      dkz = dnz*(float) l;
      v_dkz = _mm512_cvtfxpnt_round_adjustepi32_ps(_mm512_set1_epi32(l),
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dkz = _mm512_mul_ps(v_dnz,v_dkz);
      ll = nxyhd*l;
      lj = nxvyh*l;
      l1 = nxvyh*nz - lj;
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
/*          at1 = crealf(ffc[j+kk+ll])*cimagf(ffc[j+kk+ll]); */
            v_at1 = _mm512_load_ps((float *)&ffc[j+kk+ll]);
            v_at2 = (__m512)_mm512_shuffle_epi32((__m512i)v_at1,177);
            v_at1 = _mm512_mul_ps(v_at1,v_at2);
/*          at2 = at1*dnx*(float) j; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_at2 = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_at2 = _mm512_mul_ps(v_at1,_mm512_mul_ps(v_dnx,v_at2));
/*          at3 = dky*at1; */
            v_at3 = _mm512_mul_ps(v_dky,v_at1);
/*          at4 = dkz*at1; */
            v_at4 = _mm512_mul_ps(v_dkz,v_at1);
/*          zt1 = cimagf(q[j+kj+lj]) - crealf(q[j+kj+lj])*_Complex_I; */
            v_zt1 = _mm512_load_ps((float *)&q[j+kj+lj]);
            v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),
                    v_zero,v_zt1);
            v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/*          zt2 = cimagf(q[j+k1+lj]) - crealf(q[j+k1+lj])*_Complex_I; */
            v_zt2 = _mm512_load_ps((float *)&q[j+k1+lj]);
            v_zt2 = _mm512_mask_sub_ps(v_zt2,_mm512_int2mask(21845),
                    v_zero,v_zt2);
            v_zt2 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt2,177);
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),v_zero);
               v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(3),v_zero);
            }
/*          fxyz[4*(j+kj+lj)] = at2*zt1;   */
/*          fxyz[1+4*(j+kj+lj)] = at3*zt1; */
/*          fxyz[2+4*(j+kj+lj)] = at4*zt1; */
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
            _mm512_store_ps((float *)&fxyz[4*(j+kj+lj)],a);
            _mm512_store_ps((float *)&fxyz[8+4*(j+kj+lj)],b);
            _mm512_store_ps((float *)&fxyz[16+4*(j+kj+lj)],c);
            _mm512_store_ps((float *)&fxyz[24+4*(j+kj+lj)],d);
/*          fxyz[4*(j+k1+lj)] = at2*zt2;    */
/*          fxyz[1+4*(j+k1+lj)] = -at3*zt2; */
/*          fxyz[2+4*(j+k1+lj)] = at4*zt2;  */
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
            _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],a);
            _mm512_store_ps((float *)&fxyz[8+4*(j+k1+lj)],b);
            _mm512_store_ps((float *)&fxyz[16+4*(j+k1+lj)],c);
            _mm512_store_ps((float *)&fxyz[24+4*(j+k1+lj)],d);
/*          wp += at1*(q[j+kj+lj]*conjf(q[j+kj+lj]) */
/*             + q[j+k1+lj]*conjf(q[j+k1+lj]));     */
            v_zt3 = _mm512_mul_ps(v_zt1,v_zt1);
            v_zt3 = _mm512_add_ps(v_zt3,_mm512_mul_ps(v_zt2,v_zt2));
            v_zt3 = _mm512_mul_ps(v_at1,v_zt3);
/*          zt1 = cimagf(q[j+kj+l1]) - crealf(q[j+kj+l1])*_Complex_I; */
            v_zt1 = _mm512_load_ps((float *)&q[j+kj+l1]);
            v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),
                    v_zero,v_zt1);
            v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/*          zt2 = cimagf(q[j+k1+l1]) - crealf(q[j+k1+l1])*_Complex_I; */
            v_zt2 = _mm512_load_ps((float *)&q[j+k1+l1]);
            v_zt2 = _mm512_mask_sub_ps(v_zt2,_mm512_int2mask(21845),
                    v_zero,v_zt2);
            v_zt2 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt2,177);
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),v_zero);
               v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(3),v_zero);
            }
/*          fxyz[4*(j+kj+l1)] = at2*zt1;    */
/*          fxyz[1+4*(j+kj+l1)] = at3*zt1;  */
/*          fxyz[2+4*(j+kj+l1)] = -at4*zt1; */
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
            _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],a);
            _mm512_store_ps((float *)&fxyz[8+4*(j+kj+l1)],b);
            _mm512_store_ps((float *)&fxyz[16+4*(j+kj+l1)],c);
            _mm512_store_ps((float *)&fxyz[24+4*(j+kj+l1)],d);
/*          fxyz[4*(j+k1+l1)] = at2*zt2;    */
/*          fxyz[1+4*(j+k1+l1)] = -at3*zt2; */
/*          fxyz[2+4*(j+k1+l1)] = -at4*zt2; */
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
            _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],a);
            _mm512_store_ps((float *)&fxyz[8+4*(j+k1+l1)],b);
            _mm512_store_ps((float *)&fxyz[16+4*(j+k1+l1)],c);
            _mm512_store_ps((float *)&fxyz[24+4*(j+k1+l1)],d);
/*          wp += at1*(q[j+kj+l1]*conjf(q[j+kj+l1]) */
/*             + q[j+k1+l1]*conjf(q[j+k1+l1]));     */
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
/*       at1 = crealf(ffc[j+ll])*cimagf(ffc[j+ll]); */
         v_at1 = _mm512_load_ps((float *)&ffc[j+ll]);
         v_at2 = (__m512)_mm512_shuffle_epi32((__m512i)v_at1,177);
         v_at1 = _mm512_mul_ps(v_at1,v_at2);
/*       at2 = at1*dnx*(float) j; */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_at2 = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_at2 = _mm512_mul_ps(v_at1,_mm512_mul_ps(v_dnx,v_at2));
/*       at4 = dkz*at1; */
         v_at4 = _mm512_mul_ps(v_dkz,v_at1);
/*       zt1 = cimagf(q[j+lj]) - crealf(q[j+lj])*_Complex_I; */
         v_zt1 = _mm512_load_ps((float *)&q[j+lj]);
         v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),
                 v_zero,v_zt1);
         v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/*       zt2 = cimagf(q[j+l1]) - crealf(q[j+l1])*_Complex_I; */
         v_zt2 = _mm512_load_ps((float *)&q[j+l1]);
         v_zt2 = _mm512_mask_sub_ps(v_zt2,_mm512_int2mask(21845),
                 v_zero,v_zt2);
         v_zt2 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt2,177);
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),v_zero);
            v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(3),v_zero);
         }
/*       fxyz[4*(j+lj)] = at2*zt1;   */
/*       fxyz[1+4*(j+lj)] = zero;    */
/*       fxyz[2+4*(j+lj)] = at4*zt1; */
         a = _mm512_mul_ps(v_at2,v_zt1);
         b = v_zero;
         c = _mm512_mul_ps(v_at4,v_zt1);
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
         _mm512_store_ps((float *)&fxyz[4*(j+lj)],a);
         _mm512_store_ps((float *)&fxyz[8+4*(j+lj)],b);
         _mm512_store_ps((float *)&fxyz[16+4*(j+lj)],c);
         _mm512_store_ps((float *)&fxyz[24+4*(j+lj)],d);
/*       fxyz[4*(j+k1+lj)] = zero;   */
/*       fxyz[1+4*(j+k1+lj)] = zero; */
/*       fxyz[2+4*(j+k1+lj)] = zero; */
         _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],v_zero);
         _mm512_store_ps((float *)&fxyz[8+4*(j+k1+lj)],v_zero);
         _mm512_store_ps((float *)&fxyz[16+4*(j+k1+lj)],v_zero);
         _mm512_store_ps((float *)&fxyz[24+4*(j+k1+lj)],v_zero);
/*       fxyz[4*(j+l1)] = at2*zt2;    */
/*       fxyz[1+4*(j+l1)] = zero;     */
/*       fxyz[2+4*(j+l1)] = -at4*zt2; */
         a = _mm512_mul_ps(v_at2,v_zt2);
         b = v_zero;
         c = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at4,v_zt2));
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
         _mm512_store_ps((float *)&fxyz[4*(j+l1)],a);
         _mm512_store_ps((float *)&fxyz[8+4*(j+l1)],b);
         _mm512_store_ps((float *)&fxyz[16+4*(j+l1)],c);
         _mm512_store_ps((float *)&fxyz[24+4*(j+l1)],d);
/*       fxyz[4*(j+k1+l1)] = zero;   */
/*       fxyz[1+4*(j+k1+l1)] = zero; */
/*       fxyz[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[8+4*(j+k1+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[16+4*(j+k1+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[24+4*(j+k1+l1)],v_zero);
/*       wp += at1*(q[j+lj]*conjf(q[j+lj]) */
/*          + q[j+l1]*conjf(q[j+l1]));     */
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
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      for (j = 1; j < nxh; j++) {
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
   }
/* mode numbers kx = 0, nx/2 */
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
/*     fxyz[2+4*(j+l1)] = zero; */
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
/* *we = wp*((float) nx)*((float) ny)*((float) nz); */
   d0 = _mm512_reduce_add_pd(v_wp);
   *we = (wp + d0)*((float) nx)*((float) ny)*((float) nz);
   return;
}

/*--------------------------------------------------------------------*/
void ckncfft3rvxy(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nzi, int nzp, int nxhd, int nyd, int nzd,
                  int nxhyzd, int nxyzhd) {
/* this subroutine performs the x-y part of a three dimensional real to
   complex fast fourier transform and its inverse, for a subset of z,
   using complex arithmetic
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, an inverse fourier transform is performed
   f[i][m][n] = (1/nx*ny*nz)*sum(f[i][k][j]*exp(-sqrt(-1)*2pi*n*j/nx)*
         exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, a forward fourier transform is performed
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
   int nz, nxyz, nxhyz, nzt, nrx, nry, nxhyd;
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
   for (n = nzi-1; n < nzt; n++) {
      nn = nxhyd*n;
/* bit-reverse array elements in x */
      nrx = nxhyz/nxh;
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrx;
         if (j >= j1)
            continue;
         for (i = 0; i < ny; i++) {
            joff = nxhd*i + nn;
            t1 = f[j1+joff];
            f[j1+joff] = f[j+joff];
            f[j+joff] = t1;
         }
      }
/* first transform in x */
      nrx = nxyz/nxh;
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
      nry = nxhyz/ny;
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
         k1 = (mixup[k] - 1)/nry;
         if (k >= k1)
            continue;
         k1 = nxhd*k1 + nn;
/* vector loop over elements in blocks of 8 */
         for (i = 0; i < nxhs; i+=8) {
/*          t1 = f[i+k1]; */
            v_t1 = _mm512_load_ps((float *)&f[i+k1]);
/*          f[i+k1] = f[i+joff]; */
            v_t2 = _mm512_load_ps((float *)&f[i+joff]);
            _mm512_store_ps((float *)&f[i+k1],v_t2);
/*          f[i+joff] = t1; */
            _mm512_store_ps((float *)&f[i+joff],v_t1);
         }
/* loop over remaining elements */
         for (i = nxhs; i < nxh; i++) {
            t1 = f[i+k1];
            f[i+k1] = f[i+joff];
            f[i+joff] = t1;
         }
      }
/* then transform in y */
      nry = nxyz/ny;
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
L180: for (n = nzi-1; n < nzt; n++) {
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
      nry = nxhyz/ny;
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
         k1 = (mixup[k] - 1)/nry;
         if (k >= k1)
            continue;
         k1 = nxhd*k1 + nn;
/* vector loop over elements in blocks of 8 */
         for (i = 0; i < nxhs; i+=8) {
/*          t1 = f[i+k1]; */
            v_t1 = _mm512_load_ps((float *)&f[i+k1]);
/*          f[i+k1] = f[i+joff]; */
            v_t2 = _mm512_load_ps((float *)&f[i+joff]);
            _mm512_store_ps((float *)&f[i+k1],v_t2);
/*          f[i+joff] = t1; */
            _mm512_store_ps((float *)&f[i+joff],v_t1);
         }
/* loop over remaining elements */
         for (i = nxhs; i < nxh; i++) {
            t1 = f[i+k1];
            f[i+k1] = f[i+joff];
            f[i+joff] = t1;
         }
      }
/* then transform in y */
      nry = nxyz/ny;
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
      nrx = nxhyz/nxh;
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrx;
         if (j >= j1)
            continue;
         for (i = 0; i < ny; i++) {
            joff = nxhd*i + nn;
            t1 = f[j1+joff];
            f[j1+joff] = f[j+joff];
            f[j+joff] = t1;
         }
      }
/* finally transform in x */
      nrx = nxyz/nxh;
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
void ckncfft3rxz(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int nyi, int nyp, int nxhd, int nyd, int nzd,
                 int nxhyzd, int nxyzhd) {
/* this subroutine performs the z part of a three dimensional real to
   complex fast fourier transform and its inverse, for a subset of y,
   using complex arithmetic
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, an inverse fourier transform is performed
   f[l][k][j] = sum(f[i][k][j]*exp(-sqrt(-1)*2pi*l*i/nz))
   if isign = 1, a forward fourier transform is performed
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
   int nz, nzh, nxyz, nxhyz, nyt, nrz, nxhyd;
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
      goto L100;
/* inverse fourier transform */
/* bit-reverse array elements in z */
   nrz = nxhyz/nz;
   for (l = 0; l < nz; l++) {
      ll = nxhyd*l;
      l1 = (mixup[l] - 1)/nrz;
      if (l >= l1)
         continue;
      l1 = nxhyd*l1;
      for (n = nyi-1; n < nyt; n++) {
         i1 = nxhd*n;
         i0 = i1 + ll;
         i1 += l1;
/* vector loop over elements in blocks of 8 */
         for (i = 0; i < nxhs; i+=8) {
/*          t1 = f[i+i1]; */
            v_t1 = _mm512_load_ps((float *)&f[i+i1]);
/*          f[i+i1] = f[i+i0]; */
            v_t2 = _mm512_load_ps((float *)&f[i+i0]);
            _mm512_store_ps((float *)&f[i+i1],v_t2);
/*          f[i+i0] = t1; */
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
   nrz = nxyz/nz;
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
            for (n = nyi-1; n < nyt; n++) {
               i1 = nxhd*n;
               i0 = i1 + j1;
               i1 += j2;
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
      }
      ns = ns2;
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
/* scramble modes kx = 0, nx/2 */
L100: for (n = 1; n < nzh; n++) {
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
   nrz = nxhyz/nz;
   for (l = 0; l < nz; l++) {
      ll = nxhyd*l;
      l1 = (mixup[l] - 1)/nrz;
      if (l >= l1)
         continue;
      l1 = nxhyd*l1;
      for (n = nyi-1; n < nyt; n++) {
         i1 = nxhd*n;
         i0 = i1 + ll;
         i1 += l1;
/* vector loop over elements in blocks of 8 */
         for (i = 0; i < nxhs; i+=8) {
/*          t1 = f[i+i1]; */
            v_t1 = _mm512_load_ps((float *)&f[i+i1]);
/*          f[i+i1] = f[i+i0]; */
            v_t2 = _mm512_load_ps((float *)&f[i+i0]);
            _mm512_store_ps((float *)&f[i+i1],v_t2);
/*          f[i+i0] = t1; */
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
   nrz = nxyz/nz;
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
            for (n = nyi-1; n < nyt; n++) {
               i1 = nxhd*n;
               i0 = i1 + j1;
               i1 += j2;
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
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncfft3rv3xy(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int nzi, int nzp, int nxhd, int nyd, int nzd,
                   int nxhyzd, int nxyzhd) {
/* this subroutine performs the x-y part of 3 three dimensional complex
   to real fast fourier transforms and their inverses, for a subset of z,
   using complex arithmetic
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, three inverse fourier transforms is performed
   f[i][m][n][0:2] = (1/nx*ny*nz)*sum(f[i][k][j][0:2]*exp(-sqrt(-1)*2pi*n*j/nx)
         *exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, three forward fourier transforms are performed
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
   int nz, nxyz, nxhyz, nzt, nrx, nry, nxhd4, nxhyd;
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
   v_l = _mm512_set_epi32(15,11,14,10,13,9,12,8,7,3,6,2,5,1,4,0);
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
      nrx = nxhyz/nxh;
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrx;
         if (j >= j1)
            continue;
         for (i = 0; i < ny; i++) {
            joff = nxhd4*i + nn;
/*          t1 = f[4*j1+joff];   */
/*          t2 = f[1+4*j1+joff]; */
/*          t3 = f[2+4*j1+joff]; */
            v_t1 = _mm512_mask_loadunpacklo_ps(v_t1,
                   _mm512_int2mask(255),(float *)&f[4*j1+joff]);
            v_t1 = _mm512_mask_loadunpackhi_ps(v_t1,
                   _mm512_int2mask(255),(float *)&f[4*j1+joff+8]);
/*          f[4*j1+joff] = f[4*j+joff];     */
/*          f[1+4*j1+joff] = f[1+4*j+joff]; */
/*          f[2+4*j1+joff] = f[2+4*j+joff]; */
            v_t2 = _mm512_mask_loadunpacklo_ps(v_t2,
                   _mm512_int2mask(255),(float *)&f[4*j+joff]);
            v_t2 = _mm512_mask_loadunpackhi_ps(v_t2,
                   _mm512_int2mask(255),(float *)&f[4*j+joff+8]);
            _mm512_mask_packstorelo_ps((float *)&f[4*j1+joff],
               _mm512_int2mask(255),v_t2);
            _mm512_mask_packstorehi_ps((float *)&f[4*j1+joff+8],
               _mm512_int2mask(255),v_t2);
/*          f[4*j+joff] = t1;   */
/*          f[1+4*j+joff] = t2; */
/*          f[2+4*j+joff] = t3; */
            _mm512_mask_packstorelo_ps((float *)&f[4*j+joff],
               _mm512_int2mask(255),v_t1);
            _mm512_mask_packstorehi_ps((float *)&f[4*j+joff+8],
               _mm512_int2mask(255),v_t1);
         }
      }
/* first transform in x */
      nrx = nxyz/nxh;
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
      nry = nxhyz/ny;
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
         k1 = (mixup[k] - 1)/nry;
         if (k >= k1)
            continue;
         k1 = nxhd4*k1 + nn;
/* vector loop over elements in blocks of 2 */
         for (i = 0; i < nxhs; i+=2) {
/*          t1 = f[4*i+k1];   */
/*          t2 = f[1+4*i+k1]; */
/*          t3 = f[2+4*i+k1]; */
            v_t1 = _mm512_load_ps((float *)&f[4*i+k1]);
/*          f[4*i+k1] = f[4*i+joff];     */
/*          f[1+4*i+k1] = f[1+4*i+joff]; */
/*          f[2+4*i+k1] = f[2+4*i+joff]; */
            v_t2 = _mm512_load_ps((float *)&f[4*i+joff]);
            _mm512_store_ps((float *)&f[4*i+k1],v_t2);
/*          f[4*i+joff] = t1;   */
/*          f[1+4*i+joff] = t2; */
/*          f[2+4*i+joff] = t3; */
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
/* then transform in y */
      nry = nxyz/ny;
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
L230: v_l = _mm512_set_epi32(15,13,11,9,14,12,10,8,7,5,3,1,6,4,2,0);
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
      nry = nxhyz/ny;
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
         k1 = (mixup[k] - 1)/nry;
         if (k >= k1)
            continue;
         k1 = nxhd4*k1 + nn;
/* vector loop over elements in blocks of 2 */
         for (i = 0; i < nxhs; i+=2) {
/*          t1 = f[4*i+k1];   */
/*          t2 = f[1+4*i+k1]; */
/*          t3 = f[2+4*i+k1]; */
            v_t1 = _mm512_load_ps((float *)&f[4*i+k1]);
/*          f[4*i+k1] = f[4*i+joff];     */
/*          f[1+4*i+k1] = f[1+4*i+joff]; */
/*          f[2+4*i+k1] = f[2+4*i+joff]; */
            v_t2 = _mm512_load_ps((float *)&f[4*i+joff]);
            _mm512_store_ps((float *)&f[4*i+k1],v_t2);
/*          f[4*i+joff] = t1;   */
/*          f[1+4*i+joff] = t2; */
/*          f[2+4*i+joff] = t3; */
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
/* then transform in y */
      nry = nxyz/ny;
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
      nrx = nxhyz/nxh;
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrx;
         if (j >= j1)
            continue;
         for (i = 0; i < ny; i++) {
            joff = nxhd4*i + nn;
/*          t1 = f[4*j1+joff];   */
/*          t2 = f[1+4*j1+joff]; */
/*          t3 = f[2+4*j1+joff]; */
            v_t1 = _mm512_mask_loadunpacklo_ps(v_t1,
                   _mm512_int2mask(255),(float *)&f[4*j1+joff]);
            v_t1 = _mm512_mask_loadunpackhi_ps(v_t1,
                   _mm512_int2mask(255),(float *)&f[4*j1+joff+8]);
/*          f[4*j1+joff] = f[4*j+joff];     */
/*          f[1+4*j1+joff] = f[1+4*j+joff]; */
/*          f[2+4*j1+joff] = f[2+4*j+joff]; */
            v_t2 = _mm512_mask_loadunpacklo_ps(v_t2,
                   _mm512_int2mask(255),(float *)&f[4*j+joff]);
            v_t2 = _mm512_mask_loadunpackhi_ps(v_t2,
                   _mm512_int2mask(255),(float *)&f[4*j+joff+8]);
            _mm512_mask_packstorelo_ps((float *)&f[4*j1+joff],
               _mm512_int2mask(255),v_t2);
            _mm512_mask_packstorehi_ps((float *)&f[4*j1+joff+8],
               _mm512_int2mask(255),v_t2);
/*          f[4*j+joff] = t1;   */
/*          f[1+4*j+joff] = t2; */
/*          f[2+4*j+joff] = t3; */
            _mm512_mask_packstorelo_ps((float *)&f[4*j+joff],
               _mm512_int2mask(255),v_t1);
            _mm512_mask_packstorehi_ps((float *)&f[4*j+joff+8],
               _mm512_int2mask(255),v_t1);
         }
      }
/* finally transform in x */
      nrx = nxyz/nxh;
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
void ckncfft3rv3z(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nyi, int nyp, int nxhd, int nyd, int nzd,
                  int nxhyzd, int nxyzhd) {
/* this subroutine performs the z part of 3 three dimensional complex to
   real fast fourier transforms and their inverses, for a subset of y,
   using complex arithmetic
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, three inverse fourier transforms is performed
   f[l][k][j][0:2] = sum(f[i][k][j][0:2]*exp(-sqrt(-1)*2pi*l*i/nz))
   if isign = 1, three forward fourier transforms are performed
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
   int nz, nzh, nxyz, nxhyz, nyt, nrz, nxhd4, nxhyd;
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
      goto L120;
/* inverse fourier transform */
/* bit-reverse array elements in z */
   nrz = nxhyz/nz;
   for (l = 0; l < nz; l++) {
      ll = nxhyd*l;
      l1 = (mixup[l] - 1)/nrz;
      if (l >= l1)
         continue;
      l1 = nxhyd*l1;
      for (n = nyi-1; n < nyt; n++) {
         i1 = nxhd4*n;
         i0 = i1 + ll;
         i1 += l1;
/* vector loop over elements in blocks of 2 */
         for (i = 0; i < nxhs; i+=2) {
/*          t1 = f[4*i+i1];   */
/*          t2 = f[1+4*i+i1]; */
/*          t3 = f[2+4*i+i1]; */
            v_t1 = _mm512_load_ps((float *)&f[4*i+i1]);
/*          f[4*i+i1] = f[4*i+i0];     */
/*          f[1+4*i+i1] = f[1+4*i+i0]; */
/*          f[2+4*i+i1] = f[2+4*i+i0]; */
            v_t2 = _mm512_load_ps((float *)&f[4*i+i0]);
            _mm512_store_ps((float *)&f[4*i+i1],v_t2);
/*          f[4*i+i0] = t1;   */
/*          f[1+4*i+i0] = t2; */
/*          f[2+4*i+i0] = t3; */
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
   nrz = nxyz/nz;
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
            for (n = nyi-1; n < nyt; n++) {
               i1 = nxhd4*n;
               i0 = i1 + j1;
               i1 += j2;
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
      }
      ns = ns2;
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
            i1 = nxhd4*nyh;
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
/* scramble modes kx = 0, nx/2 */
L120: for (n = 1; n < nzh; n++) {
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
            i1 = nxhd4*nyh;
            i0 = i1 + ll;
            i1 += l1;
            t1 = cimagf(f[jj+i1]) + crealf(f[jj+i1])*_Complex_I;
            f[jj+i1] = conjf(f[jj+i0] - t1);
            f[jj+i0] += t1;
         }
      }
   }
/* bit-reverse array elements in z */
   nrz = nxhyz/nz;
   for (l = 0; l < nz; l++) {
      ll = nxhyd*l;
      l1 = (mixup[l] - 1)/nrz;
      if (l >= l1)
         continue;
      l1 = nxhyd*l1;
      for (n = nyi-1; n < nyt; n++) {
         i1 = nxhd4*n;
         i0 = i1 + ll;
         i1 += l1;
/* vector loop over elements in blocks of 2 */
         for (i = 0; i < nxhs; i+=2) {
/*          t1 = f[4*i+i1];   */
/*          t2 = f[1+4*i+i1]; */
/*          t3 = f[2+4*i+i1]; */
            v_t1 = _mm512_load_ps((float *)&f[4*i+i1]);
/*          f[4*i+i1] = f[4*i+i0];     */
/*          f[1+4*i+i1] = f[1+4*i+i0]; */
/*          f[2+4*i+i1] = f[2+4*i+i0]; */
            v_t2 = _mm512_load_ps((float *)&f[4*i+i0]);
            _mm512_store_ps((float *)&f[4*i+i1],v_t2);
/*          f[4*i+i0] = t1;   */
/*          f[1+4*i+i0] = t2; */
/*          f[2+4*i+i0] = t3; */
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
nrz = nxyz/nz;
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
            for (n = nyi-1; n < nyt; n++) {
               i1 = nxhd4*n;
               i0 = i1 + j1;
               i1 += j2;
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
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rvx(float complex f[], int isign, int mixup[],
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
      ckncfft3rvxy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                   nxhyzd,nxyzhd);
/* perform z fft */
      ckncfft3rxz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                  nxhyzd,nxyzhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform z fft */
      ckncfft3rxz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                  nxhyzd,nxyzhd);
/* perform xy fft */
      ckncfft3rvxy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                   nxhyzd,nxyzhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rv3(float complex f[], int isign, int mixup[],
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
      ckncfft3rv3xy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
                    nzd,nxhyzd,nxyzhd);
/* perform z fft */
      ckncfft3rv3z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                   nxhyzd,nxyzhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform z fft */
      ckncfft3rv3z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                  nxhyzd,nxyzhd);
/* perform xy fft */
      ckncfft3rv3xy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
                    nzd,nxhyzd,nxyzhd);
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void ckncxiscan2_(int *isdata, int *nths) {
   ckncxiscan2(isdata,*nths);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgpush3lt_(float *part, float *fxyz, float *qbm, float *dt,
                   float *ek, int *idimp, int *nop, int *npe, int *nx,
                   int *ny, int *nz, int *nxv, int *nyv, int *nzv,
                   int *ipbc) {
   ckncgpush3lt(part,fxyz,*qbm,*dt,ek,*idimp,*nop,*npe,*nx,*ny,*nz,*nxv,
                *nyv,*nzv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgpost3lt_(float *part, float *q, float *qm, int *nop, int *npe,
                   int *idimp, int *nxv, int *nyv, int *nzv) {
   ckncgpost3lt(part,q,*qm,*nop,*npe,*idimp,*nxv,*nyv,*nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cknc2gpost3lt_(float *part, float *q, float *qm, int *nop,
                    int *npe, int *idimp, int *nxv, int *nyv,
                    int *nzv) {
   cknc2gpost3lt(part,q,*qm,*nop,*npe,*idimp,*nxv,*nyv,*nzv);
   return;
}

/*--------------------------------------------------------------------*/
void ckncdsortp3yzlt_(float *parta, float *partb, int *npic, int *idimp,
                      int *nop, int *npe, int *ny1, int *nyz1) {
   ckncdsortp3yzlt(parta,partb,npic,*idimp,*nop,*npe,*ny1,*nyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cknccguard3l_(float *fxyz, int *nx, int *ny, int *nz, int *nxe,
                   int *nye, int *nze) {
   cknccguard3l(fxyz,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void ckncaguard3l_(float *q, int *nx, int *ny, int *nz, int *nxe,
                   int *nye, int *nze) {
   ckncaguard3l(q,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void ckncpois33_(float complex *q, float complex *fxyz, int *isign,
                 float complex *ffc, float *ax, float *ay, float *az,
                 float *affp, float *we, int *nx, int *ny, int *nz,
                 int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
                 int *nzhd) {
   ckncpois33(q,fxyz,*isign,ffc,*ax,*ay,*az,*affp,we,*nx,*ny,*nz,*nxvh,
              *nyv,*nzv,*nxhd,*nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rvx_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *indz,
                   int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                   int *nxyzhd) {
   ckncwfft3rvx(f,*isign,mixup,sct,*indx,*indy,*indz,*nxhd,*nyd,*nzd,
                *nxhyzd,*nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rv3_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *indz,
                   int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                   int *nxyzhd) {
   ckncwfft3rv3(f,*isign,mixup,sct,*indx,*indy,*indz,*nxhd,*nyd,*nzd,
                *nxhyzd,*nxyzhd);
   return;
}
